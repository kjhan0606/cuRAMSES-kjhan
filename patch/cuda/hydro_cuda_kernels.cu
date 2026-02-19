// ==========================================================================
// RAMSES Hydro Godunov Unsplit Solver -- CUDA Implementation
// Adapted from RAMSES-DARWIN-v2 for cuRAMSES (hydro-only, no RT/cooling)
//
// Pipeline:
//   1. ctoprim   : conservative -> primitive + sound speed
//   2. uslope    : TVD slope limiter (3D)
//   3. trace3d   : MUSCL-Hancock half-step predictor
//   4. flux      : Riemann solves (X, Y, Z) -> Godunov fluxes
//   5. difmag    : optional artificial viscosity (cmpdivu + consup)
//
// Memory layout follows Fortran column-major order so that the
// host arrays can be memcpy'd without transposition.
// ==========================================================================

#include "cuda_stream_pool.h"
#include <cstdio>
#include <cfloat>
#include <cstring>

// ---- Riemann solver IDs ----
#define RIEMANN_LLF  0
#define RIEMANN_HLL  1
#define RIEMANN_HLLC 2

// ====================================================================
// Constant-memory hydro parameters (set once at init)
// ====================================================================
struct HydroParams {
    double gamma;
    double smallr;
    double smallc;
    double smallp;
    double smalle;
    double entho;       // 1/(gamma-1)
    int    slope_type;
    double slope_theta;
    int    nvar;
    int    ndim;
    int    riemann_solver;
    int    pressure_fix;
    double difmag;
};

__constant__ HydroParams d_hp;

// ====================================================================
// Index macros -- Fortran column-major with first dim = ngrid
//
// Fortran: arr(1:ngrid, iu1:iu2, ju1:ju2, ku1:ku2, 1:nvar [, 1:ndim])
// C index: g + ngrid * ( (i-IU1) + NI * ( (j-JU1) + NJ * ( (k-KU1) + NK * (n-1 [+ nvar*(d-1)]) ) ) )
// ====================================================================

// uloc / q / c  -- shape (ngrid, NI, NJ, NK, NVAR)
#define IDX5(g, i, j, k, n, ng)                                        \
    ((g) + (long)(ng) * ((long)((i) - IU1)                             \
        + (long)STENCIL_NI * ((long)((j) - JU1)                       \
        + (long)STENCIL_NJ * ((long)((k) - KU1)                       \
        + (long)STENCIL_NK * (long)((n) - 1)))))

// c (sound speed) -- shape (ngrid, NI, NJ, NK)
#define IDX4(g, i, j, k, ng)                                           \
    ((g) + (long)(ng) * ((long)((i) - IU1)                             \
        + (long)STENCIL_NI * ((long)((j) - JU1)                       \
        + (long)STENCIL_NJ * (long)((k) - KU1))))

// gloc (gravity) -- shape (ngrid, NI, NJ, NK, NDIM)
#define IDX5G(g, i, j, k, d, ng) IDX5(g, i, j, k, d, ng)

// dq / qm / qp -- shape (ngrid, NI, NJ, NK, NVAR, NDIM)
#define IDX6(g, i, j, k, n, d, ng)                                     \
    ((g) + (long)(ng) * ((long)((i) - IU1)                             \
        + (long)STENCIL_NI * ((long)((j) - JU1)                       \
        + (long)STENCIL_NJ * ((long)((k) - KU1)                       \
        + (long)STENCIL_NK * ((long)((n) - 1)                         \
        + (long)NVAR * (long)((d) - 1))))))

// flux -- shape (ngrid, FLUX_NI, FLUX_NJ, FLUX_NK, NVAR, NDIM)
#define IDX_FLUX(g, i, j, k, n, d, ng)                                 \
    ((g) + (long)(ng) * ((long)((i) - IF1)                             \
        + (long)FLUX_NI * ((long)((j) - JF1)                          \
        + (long)FLUX_NJ * ((long)((k) - KF1)                          \
        + (long)FLUX_NK * ((long)((n) - 1)                            \
        + (long)NVAR * (long)((d) - 1))))))

// tmp -- shape (ngrid, FLUX_NI, FLUX_NJ, FLUX_NK, 2, NDIM)
#define IDX_TMP(g, i, j, k, n, d, ng)                                  \
    ((g) + (long)(ng) * ((long)((i) - IF1)                             \
        + (long)FLUX_NI * ((long)((j) - JF1)                          \
        + (long)FLUX_NJ * ((long)((k) - KF1)                          \
        + (long)FLUX_NK * ((long)((n) - 1)                            \
        + 2L * (long)((d) - 1))))))

// ====================================================================
// Helper: map linear thread index to (i,j,k) within the 6x6x6 stencil
// ====================================================================
__device__ __forceinline__
void thread_to_ijk(int tid, int &i, int &j, int &k) {
    i = IU1 + (tid % STENCIL_NI);
    j = JU1 + ((tid / STENCIL_NI) % STENCIL_NJ);
    k = KU1 + (tid / (STENCIL_NI * STENCIL_NJ));
}

// ====================================================================
// Kernel 1: ctoprim  (conservative -> primitive)
// One block per grid, 216 threads per block (one per cell in stencil)
// ====================================================================
__global__ void hydro_ctoprim_kernel(
    const double* __restrict__ uloc,
    const double* __restrict__ gloc,
    double* __restrict__ q,
    double* __restrict__ c,
    double dt,
    int ngrid)
{
    const int g = blockIdx.x;
    if (g >= ngrid) return;
    const int tid = threadIdx.x;
    if (tid >= NCELLS_PER_GRID) return;

    int i, j, k;
    thread_to_ijk(tid, i, j, k);

    const double gamma_   = d_hp.gamma;
    const double smallr_  = d_hp.smallr;
    const double smalle_  = d_hp.smalle;
    const double dtxhalf  = dt * 0.5;
    const int    nvar_    = d_hp.nvar;
    const int    ndim_    = d_hp.ndim;

    // Density
    double rho = uloc[IDX5(g, i, j, k, 1, ngrid)];
    rho = fmax(rho, smallr_);
    q[IDX5(g, i, j, k, 1, ngrid)] = rho;

    double oneoverrho = 1.0 / rho;

    // Velocities  (u,v,w -> vars 2,3,4 for NDIM=3)
    double vel[3];
    for (int d = 0; d < ndim_; d++) {
        vel[d] = uloc[IDX5(g, i, j, k, d + 2, ngrid)] * oneoverrho;
        q[IDX5(g, i, j, k, d + 2, ngrid)] = vel[d];
    }

    // Kinetic energy per unit mass
    double eken = 0.0;
    for (int d = 0; d < ndim_; d++) {
        eken += 0.5 * vel[d] * vel[d];
    }

    // Non-thermal energy (NENER not used: erad=0)
    double erad = 0.0;

    // Thermal specific internal energy
    double eint = uloc[IDX5(g, i, j, k, ndim_ + 2, ngrid)] * oneoverrho - eken - erad;
    eint = fmax(eint, smalle_);

    // Thermal pressure
    double pressure = (gamma_ - 1.0) * rho * eint;
    q[IDX5(g, i, j, k, ndim_ + 2, ngrid)] = pressure;

    // Sound speed  c = sqrt(gamma * P / rho)
    double cs2 = gamma_ * pressure * oneoverrho;
    c[IDX4(g, i, j, k, ngrid)] = sqrt(cs2);

    // Gravity predictor  v += g * dt/2
    for (int d = 0; d < ndim_; d++) {
        q[IDX5(g, i, j, k, d + 2, ngrid)] += gloc[IDX5G(g, i, j, k, d + 1, ngrid)] * dtxhalf;
    }

    // Passive scalars  (indices ndim+3 .. nvar, 1-based)
    for (int n = ndim_ + 3; n <= nvar_; n++) {
        q[IDX5(g, i, j, k, n, ngrid)] = uloc[IDX5(g, i, j, k, n, ngrid)] * oneoverrho;
    }
}

// ====================================================================
// Kernel 2: uslope  (TVD slope computation, 3D)
// Interior cells only: i in [0,3], j in [0,3], k in [0,3]
// ====================================================================
__global__ void hydro_uslope_kernel(
    const double* __restrict__ q,
    double* __restrict__ dq,
    int ngrid)
{
    const int g = blockIdx.x;
    if (g >= ngrid) return;
    const int tid = threadIdx.x;
    if (tid >= NCELLS_PER_GRID) return;

    int i, j, k;
    thread_to_ijk(tid, i, j, k);

    // Interior range for slopes: [iu1+1, iu2-1] = [0, 3]
    const int ilo = 0, ihi = 3;
    const int jlo = 0, jhi = 3;
    const int klo = 0, khi = 3;

    if (i < ilo || i > ihi || j < jlo || j > jhi || k < klo || k > khi) {
        // Zero slopes for boundary cells
        for (int n = 1; n <= NVAR; n++) {
            dq[IDX6(g, i, j, k, n, 1, ngrid)] = 0.0;
            dq[IDX6(g, i, j, k, n, 2, ngrid)] = 0.0;
            dq[IDX6(g, i, j, k, n, 3, ngrid)] = 0.0;
        }
        return;
    }

    const int    stype  = d_hp.slope_type;
    const double stheta = d_hp.slope_theta;
    const int    nvar_  = d_hp.nvar;

    if (stype == 0) {
        for (int n = 1; n <= nvar_; n++) {
            dq[IDX6(g, i, j, k, n, 1, ngrid)] = 0.0;
            dq[IDX6(g, i, j, k, n, 2, ngrid)] = 0.0;
            dq[IDX6(g, i, j, k, n, 3, ngrid)] = 0.0;
        }
        return;
    }

    for (int n = 1; n <= nvar_; n++) {
        double qc = q[IDX5(g, i, j, k, n, ngrid)];

        if (stype == 1) {
            // -- minmod limiter --
            // X direction
            {
                double dlft = qc - q[IDX5(g, i-1, j, k, n, ngrid)];
                double drgt = q[IDX5(g, i+1, j, k, n, ngrid)] - qc;
                double val;
                if (dlft * drgt <= 0.0) val = 0.0;
                else if (dlft > 0.0) val = fmin(dlft, drgt);
                else val = fmax(dlft, drgt);
                dq[IDX6(g, i, j, k, n, 1, ngrid)] = val;
            }
            // Y direction
            {
                double dlft = qc - q[IDX5(g, i, j-1, k, n, ngrid)];
                double drgt = q[IDX5(g, i, j+1, k, n, ngrid)] - qc;
                double val;
                if (dlft * drgt <= 0.0) val = 0.0;
                else if (dlft > 0.0) val = fmin(dlft, drgt);
                else val = fmax(dlft, drgt);
                dq[IDX6(g, i, j, k, n, 2, ngrid)] = val;
            }
            // Z direction
            {
                double dlft = qc - q[IDX5(g, i, j, k-1, n, ngrid)];
                double drgt = q[IDX5(g, i, j, k+1, n, ngrid)] - qc;
                double val;
                if (dlft * drgt <= 0.0) val = 0.0;
                else if (dlft > 0.0) val = fmin(dlft, drgt);
                else val = fmax(dlft, drgt);
                dq[IDX6(g, i, j, k, n, 3, ngrid)] = val;
            }

        } else if (stype == 2) {
            // -- moncen limiter --
            // X direction
            {
                double dlft = 2.0 * (qc - q[IDX5(g, i-1, j, k, n, ngrid)]);
                double drgt = 2.0 * (q[IDX5(g, i+1, j, k, n, ngrid)] - qc);
                double dcen = 0.5 * (dlft + drgt) / 2.0;
                double dsgn = (dcen >= 0.0) ? 1.0 : -1.0;
                double slop = fmin(fabs(dlft), fabs(drgt));
                double dlim = (dlft * drgt > 0.0) ? slop : 0.0;
                dq[IDX6(g, i, j, k, n, 1, ngrid)] = dsgn * fmin(dlim, fabs(dcen));
            }
            // Y direction
            {
                double dlft = 2.0 * (qc - q[IDX5(g, i, j-1, k, n, ngrid)]);
                double drgt = 2.0 * (q[IDX5(g, i, j+1, k, n, ngrid)] - qc);
                double dcen = 0.5 * (dlft + drgt) / 2.0;
                double dsgn = (dcen >= 0.0) ? 1.0 : -1.0;
                double slop = fmin(fabs(dlft), fabs(drgt));
                double dlim = (dlft * drgt > 0.0) ? slop : 0.0;
                dq[IDX6(g, i, j, k, n, 2, ngrid)] = dsgn * fmin(dlim, fabs(dcen));
            }
            // Z direction
            {
                double dlft = 2.0 * (qc - q[IDX5(g, i, j, k-1, n, ngrid)]);
                double drgt = 2.0 * (q[IDX5(g, i, j, k+1, n, ngrid)] - qc);
                double dcen = 0.5 * (dlft + drgt) / 2.0;
                double dsgn = (dcen >= 0.0) ? 1.0 : -1.0;
                double slop = fmin(fabs(dlft), fabs(drgt));
                double dlim = (dlft * drgt > 0.0) ? slop : 0.0;
                dq[IDX6(g, i, j, k, n, 3, ngrid)] = dsgn * fmin(dlim, fabs(dcen));
            }

        } else if (stype == 3) {
            // -- positivity-preserving 3D unsplit slope --
            double vmin =  1e99;
            double vmax = -1e99;
            for (int dk = -1; dk <= 1; dk++) {
                for (int dj = -1; dj <= 1; dj++) {
                    for (int di = -1; di <= 1; di++) {
                        double diff = q[IDX5(g, i+di, j+dj, k+dk, n, ngrid)] - qc;
                        vmin = fmin(vmin, diff);
                        vmax = fmax(vmax, diff);
                    }
                }
            }
            double dfx = 0.5 * (q[IDX5(g, i+1, j, k, n, ngrid)] - q[IDX5(g, i-1, j, k, n, ngrid)]);
            double dfy = 0.5 * (q[IDX5(g, i, j+1, k, n, ngrid)] - q[IDX5(g, i, j-1, k, n, ngrid)]);
            double dfz = 0.5 * (q[IDX5(g, i, j, k+1, n, ngrid)] - q[IDX5(g, i, j, k-1, n, ngrid)]);
            double dff = 0.5 * (fabs(dfx) + fabs(dfy) + fabs(dfz));
            double slop;
            if (dff > 0.0) slop = fmin(1.0, fmin(fabs(vmin), fabs(vmax)) / dff);
            else slop = 1.0;
            dq[IDX6(g, i, j, k, n, 1, ngrid)] = slop * dfx;
            dq[IDX6(g, i, j, k, n, 2, ngrid)] = slop * dfy;
            dq[IDX6(g, i, j, k, n, 3, ngrid)] = slop * dfz;

        } else if (stype == 7) {
            // -- van Leer --
            // X
            {
                double dlft = qc - q[IDX5(g, i-1, j, k, n, ngrid)];
                double drgt = q[IDX5(g, i+1, j, k, n, ngrid)] - qc;
                dq[IDX6(g, i, j, k, n, 1, ngrid)] =
                    (dlft * drgt <= 0.0) ? 0.0 : (2.0 * dlft * drgt / (dlft + drgt));
            }
            // Y
            {
                double dlft = qc - q[IDX5(g, i, j-1, k, n, ngrid)];
                double drgt = q[IDX5(g, i, j+1, k, n, ngrid)] - qc;
                dq[IDX6(g, i, j, k, n, 2, ngrid)] =
                    (dlft * drgt <= 0.0) ? 0.0 : (2.0 * dlft * drgt / (dlft + drgt));
            }
            // Z
            {
                double dlft = qc - q[IDX5(g, i, j, k-1, n, ngrid)];
                double drgt = q[IDX5(g, i, j, k+1, n, ngrid)] - qc;
                dq[IDX6(g, i, j, k, n, 3, ngrid)] =
                    (dlft * drgt <= 0.0) ? 0.0 : (2.0 * dlft * drgt / (dlft + drgt));
            }

        } else if (stype == 8) {
            // -- generalized moncen (slope_theta parameterisation) --
            // X
            {
                double dlft = qc - q[IDX5(g, i-1, j, k, n, ngrid)];
                double drgt = q[IDX5(g, i+1, j, k, n, ngrid)] - qc;
                double dcen = 0.5 * (dlft + drgt);
                double dsgn = (dcen >= 0.0) ? 1.0 : -1.0;
                double slop = fmin(stheta * fabs(dlft), stheta * fabs(drgt));
                double dlim = (dlft * drgt > 0.0) ? slop : 0.0;
                dq[IDX6(g, i, j, k, n, 1, ngrid)] = dsgn * fmin(dlim, fabs(dcen));
            }
            // Y
            {
                double dlft = qc - q[IDX5(g, i, j-1, k, n, ngrid)];
                double drgt = q[IDX5(g, i, j+1, k, n, ngrid)] - qc;
                double dcen = 0.5 * (dlft + drgt);
                double dsgn = (dcen >= 0.0) ? 1.0 : -1.0;
                double slop = fmin(stheta * fabs(dlft), stheta * fabs(drgt));
                double dlim = (dlft * drgt > 0.0) ? slop : 0.0;
                dq[IDX6(g, i, j, k, n, 2, ngrid)] = dsgn * fmin(dlim, fabs(dcen));
            }
            // Z
            {
                double dlft = qc - q[IDX5(g, i, j, k-1, n, ngrid)];
                double drgt = q[IDX5(g, i, j, k+1, n, ngrid)] - qc;
                double dcen = 0.5 * (dlft + drgt);
                double dsgn = (dcen >= 0.0) ? 1.0 : -1.0;
                double slop = fmin(stheta * fabs(dlft), stheta * fabs(drgt));
                double dlim = (dlft * drgt > 0.0) ? slop : 0.0;
                dq[IDX6(g, i, j, k, n, 3, ngrid)] = dsgn * fmin(dlim, fabs(dcen));
            }

        } else {
            // Unsupported slope type: set to zero
            dq[IDX6(g, i, j, k, n, 1, ngrid)] = 0.0;
            dq[IDX6(g, i, j, k, n, 2, ngrid)] = 0.0;
            dq[IDX6(g, i, j, k, n, 3, ngrid)] = 0.0;
        }
    }
}

// ====================================================================
// Kernel 3: trace3d  (MUSCL-Hancock half-step predictor, 3D)
// Interior range: [0,3] in each direction
//
// Variable ordering (1-based, NDIM=3):
//   ir=1, iu=2, iv=3, iw=4, ip=5
// ====================================================================
__global__ void hydro_trace3d_kernel(
    const double* __restrict__ q,
    const double* __restrict__ dq,
    double* __restrict__ qm,
    double* __restrict__ qp,
    double dtdx, double dtdy, double dtdz,
    int ngrid)
{
    const int g = blockIdx.x;
    if (g >= ngrid) return;
    const int tid = threadIdx.x;
    if (tid >= NCELLS_PER_GRID) return;

    int i, j, k;
    thread_to_ijk(tid, i, j, k);

    const int ilo = 0, ihi = 3;
    const int jlo = 0, jhi = 3;
    const int klo = 0, khi = 3;

    if (i < ilo || i > ihi || j < jlo || j > jhi || k < klo || k > khi) return;

    const double smallr_ = d_hp.smallr;
    const double gamma_  = d_hp.gamma;
    const int    nvar_   = d_hp.nvar;
    const int    ndim_   = d_hp.ndim;

    // Primitive variable indices (1-based)
    const int ir = 1;
    const int iu = 2;
    const int iv = 3;
    const int iw = 4;
    const int ip = 5;  // ndim+2 for ndim=3

    // Cell-centered primitive values
    double r = q[IDX5(g, i, j, k, ir, ngrid)];
    double u = q[IDX5(g, i, j, k, iu, ngrid)];
    double v = q[IDX5(g, i, j, k, iv, ngrid)];
    double w = q[IDX5(g, i, j, k, iw, ngrid)];
    double p = q[IDX5(g, i, j, k, ip, ngrid)];

    // TVD slopes in X
    double drx = dq[IDX6(g, i, j, k, ir, 1, ngrid)];
    double dpx = dq[IDX6(g, i, j, k, ip, 1, ngrid)];
    double dux = dq[IDX6(g, i, j, k, iu, 1, ngrid)];
    double dvx = dq[IDX6(g, i, j, k, iv, 1, ngrid)];
    double dwx = dq[IDX6(g, i, j, k, iw, 1, ngrid)];

    // TVD slopes in Y
    double dry = dq[IDX6(g, i, j, k, ir, 2, ngrid)];
    double dpy = dq[IDX6(g, i, j, k, ip, 2, ngrid)];
    double duy = dq[IDX6(g, i, j, k, iu, 2, ngrid)];
    double dvy = dq[IDX6(g, i, j, k, iv, 2, ngrid)];
    double dwy = dq[IDX6(g, i, j, k, iw, 2, ngrid)];

    // TVD slopes in Z
    double drz = dq[IDX6(g, i, j, k, ir, 3, ngrid)];
    double dpz = dq[IDX6(g, i, j, k, ip, 3, ngrid)];
    double duz = dq[IDX6(g, i, j, k, iu, 3, ngrid)];
    double dvz = dq[IDX6(g, i, j, k, iv, 3, ngrid)];
    double dwz = dq[IDX6(g, i, j, k, iw, 3, ngrid)];

    // Source terms (including ALL transverse derivatives)
    double divv = dux + dvy + dwz;
    double sr0 = -u * drx - v * dry - w * drz - divv * r;
    double sp0 = -u * dpx - v * dpy - w * dpz - divv * gamma_ * p;
    double su0 = -u * dux - v * duy - w * duz - dpx / r;
    double sv0 = -u * dvx - v * dvy - w * dvz - dpy / r;
    double sw0 = -u * dwx - v * dwy - w * dwz - dpz / r;

    // ---- X direction (dim=1) ----
    // Right state at left interface (qp)
    double qp_r = r - 0.5 * drx + sr0 * dtdx * 0.5;
    double qp_p = p - 0.5 * dpx + sp0 * dtdx * 0.5;
    double qp_u = u - 0.5 * dux + su0 * dtdx * 0.5;
    double qp_v = v - 0.5 * dvx + sv0 * dtdx * 0.5;
    double qp_w = w - 0.5 * dwx + sw0 * dtdx * 0.5;
    if (qp_r < smallr_) qp_r = r;

    qp[IDX6(g, i, j, k, ir, 1, ngrid)] = qp_r;
    qp[IDX6(g, i, j, k, ip, 1, ngrid)] = qp_p;
    qp[IDX6(g, i, j, k, iu, 1, ngrid)] = qp_u;
    qp[IDX6(g, i, j, k, iv, 1, ngrid)] = qp_v;
    qp[IDX6(g, i, j, k, iw, 1, ngrid)] = qp_w;

    // Left state at right interface (qm)
    double qm_r = r + 0.5 * drx + sr0 * dtdx * 0.5;
    double qm_p = p + 0.5 * dpx + sp0 * dtdx * 0.5;
    double qm_u = u + 0.5 * dux + su0 * dtdx * 0.5;
    double qm_v = v + 0.5 * dvx + sv0 * dtdx * 0.5;
    double qm_w = w + 0.5 * dwx + sw0 * dtdx * 0.5;
    if (qm_r < smallr_) qm_r = r;

    qm[IDX6(g, i, j, k, ir, 1, ngrid)] = qm_r;
    qm[IDX6(g, i, j, k, ip, 1, ngrid)] = qm_p;
    qm[IDX6(g, i, j, k, iu, 1, ngrid)] = qm_u;
    qm[IDX6(g, i, j, k, iv, 1, ngrid)] = qm_v;
    qm[IDX6(g, i, j, k, iw, 1, ngrid)] = qm_w;

    // ---- Y direction (dim=2) ----
    qp_r = r - 0.5 * dry + sr0 * dtdy * 0.5;
    qp_p = p - 0.5 * dpy + sp0 * dtdy * 0.5;
    qp_u = u - 0.5 * duy + su0 * dtdy * 0.5;
    qp_v = v - 0.5 * dvy + sv0 * dtdy * 0.5;
    qp_w = w - 0.5 * dwy + sw0 * dtdy * 0.5;
    if (qp_r < smallr_) qp_r = r;

    qp[IDX6(g, i, j, k, ir, 2, ngrid)] = qp_r;
    qp[IDX6(g, i, j, k, ip, 2, ngrid)] = qp_p;
    qp[IDX6(g, i, j, k, iu, 2, ngrid)] = qp_u;
    qp[IDX6(g, i, j, k, iv, 2, ngrid)] = qp_v;
    qp[IDX6(g, i, j, k, iw, 2, ngrid)] = qp_w;

    qm_r = r + 0.5 * dry + sr0 * dtdy * 0.5;
    qm_p = p + 0.5 * dpy + sp0 * dtdy * 0.5;
    qm_u = u + 0.5 * duy + su0 * dtdy * 0.5;
    qm_v = v + 0.5 * dvy + sv0 * dtdy * 0.5;
    qm_w = w + 0.5 * dwy + sw0 * dtdy * 0.5;
    if (qm_r < smallr_) qm_r = r;

    qm[IDX6(g, i, j, k, ir, 2, ngrid)] = qm_r;
    qm[IDX6(g, i, j, k, ip, 2, ngrid)] = qm_p;
    qm[IDX6(g, i, j, k, iu, 2, ngrid)] = qm_u;
    qm[IDX6(g, i, j, k, iv, 2, ngrid)] = qm_v;
    qm[IDX6(g, i, j, k, iw, 2, ngrid)] = qm_w;

    // ---- Z direction (dim=3) ----
    qp_r = r - 0.5 * drz + sr0 * dtdz * 0.5;
    qp_p = p - 0.5 * dpz + sp0 * dtdz * 0.5;
    qp_u = u - 0.5 * duz + su0 * dtdz * 0.5;
    qp_v = v - 0.5 * dvz + sv0 * dtdz * 0.5;
    qp_w = w - 0.5 * dwz + sw0 * dtdz * 0.5;
    if (qp_r < smallr_) qp_r = r;

    qp[IDX6(g, i, j, k, ir, 3, ngrid)] = qp_r;
    qp[IDX6(g, i, j, k, ip, 3, ngrid)] = qp_p;
    qp[IDX6(g, i, j, k, iu, 3, ngrid)] = qp_u;
    qp[IDX6(g, i, j, k, iv, 3, ngrid)] = qp_v;
    qp[IDX6(g, i, j, k, iw, 3, ngrid)] = qp_w;

    qm_r = r + 0.5 * drz + sr0 * dtdz * 0.5;
    qm_p = p + 0.5 * dpz + sp0 * dtdz * 0.5;
    qm_u = u + 0.5 * duz + su0 * dtdz * 0.5;
    qm_v = v + 0.5 * dvz + sv0 * dtdz * 0.5;
    qm_w = w + 0.5 * dwz + sw0 * dtdz * 0.5;
    if (qm_r < smallr_) qm_r = r;

    qm[IDX6(g, i, j, k, ir, 3, ngrid)] = qm_r;
    qm[IDX6(g, i, j, k, ip, 3, ngrid)] = qm_p;
    qm[IDX6(g, i, j, k, iu, 3, ngrid)] = qm_u;
    qm[IDX6(g, i, j, k, iv, 3, ngrid)] = qm_v;
    qm[IDX6(g, i, j, k, iw, 3, ngrid)] = qm_w;

    // ---- Passive scalars (ndim+3 .. nvar) ----
    for (int n = ndim_ + 3; n <= nvar_; n++) {
        double a   = q[IDX5(g, i, j, k, n, ngrid)];
        double uu  = q[IDX5(g, i, j, k, iu, ngrid)];
        double vv  = q[IDX5(g, i, j, k, iv, ngrid)];
        double ww  = q[IDX5(g, i, j, k, iw, ngrid)];
        double dax = dq[IDX6(g, i, j, k, n, 1, ngrid)];
        double day = dq[IDX6(g, i, j, k, n, 2, ngrid)];
        double daz = dq[IDX6(g, i, j, k, n, 3, ngrid)];
        double sa0 = -uu * dax - vv * day - ww * daz;

        qp[IDX6(g, i, j, k, n, 1, ngrid)] = a - 0.5 * dax + sa0 * dtdx * 0.5;
        qm[IDX6(g, i, j, k, n, 1, ngrid)] = a + 0.5 * dax + sa0 * dtdx * 0.5;
        qp[IDX6(g, i, j, k, n, 2, ngrid)] = a - 0.5 * day + sa0 * dtdy * 0.5;
        qm[IDX6(g, i, j, k, n, 2, ngrid)] = a + 0.5 * day + sa0 * dtdy * 0.5;
        qp[IDX6(g, i, j, k, n, 3, ngrid)] = a - 0.5 * daz + sa0 * dtdz * 0.5;
        qm[IDX6(g, i, j, k, n, 3, ngrid)] = a + 0.5 * daz + sa0 * dtdz * 0.5;
    }
}

// ====================================================================
// Device Riemann Solvers
// ====================================================================

// --- LLF (Local Lax-Friedrichs) ---
__device__ void riemann_llf_device(
    const double* __restrict__ ql,
    const double* __restrict__ qr,
    double* __restrict__ fgdnv,
    int nvar_)
{
    const double gamma_  = d_hp.gamma;
    const double smallr_ = d_hp.smallr;
    const double smallc_ = d_hp.smallc;
    const double entho   = d_hp.entho;
    const double smallp  = smallc_ * smallc_ / gamma_;
    const int    ndim_   = d_hp.ndim;

    double rl = fmax(ql[1], smallr_);
    double ul = ql[2];
    double pl = fmax(ql[3], rl * smallp);
    double cl = sqrt(gamma_ * pl / rl);

    double rr = fmax(qr[1], smallr_);
    double ur = qr[2];
    double pr = fmax(qr[3], rr * smallp);
    double cr = sqrt(gamma_ * pr / rr);

    double cmax = fmax(fabs(ul) + cl, fabs(ur) + cr);

    double ucons_l[NVAR + 1];
    double ucons_r[NVAR + 1];
    double fl[NVAR + 1];
    double fr[NVAR + 1];

    ucons_l[1] = rl;
    ucons_r[1] = rr;
    ucons_l[2] = rl * ul;
    ucons_r[2] = rr * ur;
    ucons_l[3] = pl * entho + 0.5 * rl * ul * ul;
    ucons_r[3] = pr * entho + 0.5 * rr * ur * ur;
    if (ndim_ > 1) {
        ucons_l[3] += 0.5 * rl * ql[4] * ql[4];
        ucons_r[3] += 0.5 * rr * qr[4] * qr[4];
    }
    if (ndim_ > 2) {
        ucons_l[3] += 0.5 * rl * ql[5] * ql[5];
        ucons_r[3] += 0.5 * rr * qr[5] * qr[5];
    }
    for (int n = 4; n <= ndim_ + 2; n++) {
        ucons_l[n] = rl * ql[n];
        ucons_r[n] = rr * qr[n];
    }
    for (int n = ndim_ + 3; n <= nvar_; n++) {
        ucons_l[n] = rl * ql[n];
        ucons_r[n] = rr * qr[n];
    }
    ucons_l[nvar_ + 1] = pl * entho;
    ucons_r[nvar_ + 1] = pr * entho;

    fl[1] = ul * ucons_l[1];
    fr[1] = ur * ucons_r[1];
    fl[2] = ul * ucons_l[2] + pl;
    fr[2] = ur * ucons_r[2] + pr;
    fl[3] = ul * (ucons_l[3] + pl);
    fr[3] = ur * (ucons_r[3] + pr);
    for (int n = 4; n <= nvar_ + 1; n++) {
        fl[n] = ul * ucons_l[n];
        fr[n] = ur * ucons_r[n];
    }

    for (int n = 1; n <= nvar_ + 1; n++) {
        fgdnv[n] = 0.5 * (fl[n] + fr[n] - cmax * (ucons_r[n] - ucons_l[n]));
    }
}

// --- HLL ---
__device__ void riemann_hll_device(
    const double* __restrict__ ql,
    const double* __restrict__ qr,
    double* __restrict__ fgdnv,
    int nvar_)
{
    const double gamma_  = d_hp.gamma;
    const double smallr_ = d_hp.smallr;
    const double smallc_ = d_hp.smallc;
    const double entho   = d_hp.entho;
    const double smallp  = smallc_ * smallc_ / gamma_;
    const int    ndim_   = d_hp.ndim;

    double rl = fmax(ql[1], smallr_);
    double ul = ql[2];
    double pl = fmax(ql[3], rl * smallp);
    double cl = sqrt(gamma_ * pl / rl);

    double rr = fmax(qr[1], smallr_);
    double ur = qr[2];
    double pr = fmax(qr[3], rr * smallp);
    double cr = sqrt(gamma_ * pr / rr);

    double SL = fmin(fmin(ul, ur) - fmax(cl, cr), 0.0);
    double SR = fmax(fmax(ul, ur) + fmax(cl, cr), 0.0);

    double ucons_l[NVAR + 1], ucons_r[NVAR + 1];
    double fl[NVAR + 1], fr[NVAR + 1];

    ucons_l[1] = rl;               ucons_r[1] = rr;
    ucons_l[2] = rl * ul;          ucons_r[2] = rr * ur;
    ucons_l[3] = pl * entho + 0.5 * rl * ul * ul;
    ucons_r[3] = pr * entho + 0.5 * rr * ur * ur;
    if (ndim_ > 1) {
        ucons_l[3] += 0.5 * rl * ql[4] * ql[4];
        ucons_r[3] += 0.5 * rr * qr[4] * qr[4];
    }
    if (ndim_ > 2) {
        ucons_l[3] += 0.5 * rl * ql[5] * ql[5];
        ucons_r[3] += 0.5 * rr * qr[5] * qr[5];
    }
    for (int n = 4; n <= ndim_ + 2; n++) {
        ucons_l[n] = rl * ql[n];
        ucons_r[n] = rr * qr[n];
    }
    for (int n = ndim_ + 3; n <= nvar_; n++) {
        ucons_l[n] = rl * ql[n];
        ucons_r[n] = rr * qr[n];
    }
    ucons_l[nvar_ + 1] = pl * entho;
    ucons_r[nvar_ + 1] = pr * entho;

    fl[1] = ucons_l[2];            fr[1] = ucons_r[2];
    fl[2] = pl + ucons_l[2] * ul;  fr[2] = pr + ucons_r[2] * ur;
    fl[3] = ul * (ucons_l[3] + pl); fr[3] = ur * (ucons_r[3] + pr);
    for (int n = 4; n <= nvar_ + 1; n++) {
        fl[n] = ul * ucons_l[n];
        fr[n] = ur * ucons_r[n];
    }

    double dSR = 1.0 / (SR - SL);
    for (int n = 1; n <= nvar_ + 1; n++) {
        fgdnv[n] = (SR * fl[n] - SL * fr[n] + SR * SL * (ucons_r[n] - ucons_l[n])) * dSR;
    }
}

// --- HLLC (Toro) ---
__device__ void riemann_hllc_device(
    const double* __restrict__ ql,
    const double* __restrict__ qr,
    double* __restrict__ fgdnv,
    int nvar_)
{
    const double gamma_  = d_hp.gamma;
    const double smallr_ = d_hp.smallr;
    const double smallc_ = d_hp.smallc;
    const double entho   = d_hp.entho;
    const double smallp  = smallc_ * smallc_ / gamma_;
    const int    ndim_   = d_hp.ndim;

    double rl = fmax(ql[1], smallr_);
    double Pl = fmax(ql[3], rl * smallp);
    double ul = ql[2];

    double el     = Pl * entho;
    double ecinl  = 0.5 * rl * ul * ul;
    if (ndim_ > 1) ecinl += 0.5 * rl * ql[4] * ql[4];
    if (ndim_ > 2) ecinl += 0.5 * rl * ql[5] * ql[5];
    double etotl  = el + ecinl;
    double Ptotl  = Pl;

    double rr = fmax(qr[1], smallr_);
    double Pr = fmax(qr[3], rr * smallp);
    double ur = qr[2];

    double er     = Pr * entho;
    double ecinr  = 0.5 * rr * ur * ur;
    if (ndim_ > 1) ecinr += 0.5 * rr * qr[4] * qr[4];
    if (ndim_ > 2) ecinr += 0.5 * rr * qr[5] * qr[5];
    double etotr  = er + ecinr;
    double Ptotr  = Pr;

    double cfastl = sqrt(fmax(gamma_ * Pl / rl, smallc_ * smallc_));
    double cfastr = sqrt(fmax(gamma_ * Pr / rr, smallc_ * smallc_));

    double SL = fmin(ul, ur) - fmax(cfastl, cfastr);
    double SR = fmax(ul, ur) + fmax(cfastl, cfastr);

    double rcl = rl * (ul - SL);
    double rcr = rr * (SR - ur);

    double ustar    = (rcr * ur + rcl * ul + (Ptotl - Ptotr)) / (rcr + rcl);
    double Ptotstar = (rcr * Ptotl + rcl * Ptotr + rcl * rcr * (ul - ur)) / (rcr + rcl);

    double rstarl     = rl * (SL - ul) / (SL - ustar);
    double etotstarl  = ((SL - ul) * etotl - Ptotl * ul + Ptotstar * ustar) / (SL - ustar);
    double estarl     = el * (SL - ul) / (SL - ustar);

    double rstarr     = rr * (SR - ur) / (SR - ustar);
    double etotstarr  = ((SR - ur) * etotr - Ptotr * ur + Ptotstar * ustar) / (SR - ustar);
    double estarr     = er * (SR - ur) / (SR - ustar);

    double ro, uo, Ptoto, etoto, eo;

    if (SL > 0.0) {
        ro = rl; uo = ul; Ptoto = Ptotl; etoto = etotl; eo = el;
    } else if (ustar > 0.0) {
        ro = rstarl; uo = ustar; Ptoto = Ptotstar; etoto = etotstarl; eo = estarl;
    } else if (SR > 0.0) {
        ro = rstarr; uo = ustar; Ptoto = Ptotstar; etoto = etotstarr; eo = estarr;
    } else {
        ro = rr; uo = ur; Ptoto = Ptotr; etoto = etotr; eo = er;
    }

    fgdnv[1] = ro * uo;
    fgdnv[2] = ro * uo * uo + Ptoto;
    fgdnv[3] = (etoto + Ptoto) * uo;

    for (int n = 4; n <= ndim_ + 2; n++) {
        fgdnv[n] = (ustar > 0.0) ? ro * uo * ql[n] : ro * uo * qr[n];
    }
    for (int n = ndim_ + 3; n <= nvar_; n++) {
        fgdnv[n] = (ustar > 0.0) ? ro * uo * ql[n] : ro * uo * qr[n];
    }
    fgdnv[nvar_ + 1] = uo * eo;
}

// ====================================================================
// Kernel 4: hydro_flux_kernel
// Combines cmpflxm for X, Y, Z directions.
// ====================================================================

// Flux active ranges for non-face dimensions
#define FLUX_ILO 1
#define FLUX_IHI 2
#define FLUX_JLO 1
#define FLUX_JHI 2
#define FLUX_KLO 1
#define FLUX_KHI 2

__global__ void hydro_flux_kernel(
    const double* __restrict__ qm,
    const double* __restrict__ qp,
    double* __restrict__ flux,
    double* __restrict__ tmp,
    double dtdx, double dtdy, double dtdz,
    int ngrid)
{
    const int g = blockIdx.x;
    if (g >= ngrid) return;
    const int tid = threadIdx.x;
    if (tid >= NCELLS_PER_GRID) return;

    int i, j, k;
    thread_to_ijk(tid, i, j, k);

    const int    nvar_    = d_hp.nvar;
    const int    ndim_    = d_hp.ndim;
    const int    solver   = d_hp.riemann_solver;

    double qleft[NVAR + 2];
    double qright[NVAR + 2];
    double fgdnv[NVAR + 2];

    // ---- X-direction flux ----
    if (i >= IF1 && i <= IF2 && j >= FLUX_JLO && j <= FLUX_JHI && k >= FLUX_KLO && k <= FLUX_KHI) {
        const int xdim = 1;
        const int ln = 2, lt1 = 3, lt2 = 4;

        qleft[1] = qm[IDX6(g, i-1, j, k, 1,   xdim, ngrid)];
        qleft[2] = qm[IDX6(g, i-1, j, k, ln,  xdim, ngrid)];
        qleft[3] = qm[IDX6(g, i-1, j, k, ndim_ + 2, xdim, ngrid)];
        qleft[4] = qm[IDX6(g, i-1, j, k, lt1, xdim, ngrid)];
        qleft[5] = qm[IDX6(g, i-1, j, k, lt2, xdim, ngrid)];

        qright[1] = qp[IDX6(g, i, j, k, 1,   xdim, ngrid)];
        qright[2] = qp[IDX6(g, i, j, k, ln,  xdim, ngrid)];
        qright[3] = qp[IDX6(g, i, j, k, ndim_ + 2, xdim, ngrid)];
        qright[4] = qp[IDX6(g, i, j, k, lt1, xdim, ngrid)];
        qright[5] = qp[IDX6(g, i, j, k, lt2, xdim, ngrid)];

        for (int n = ndim_ + 3; n <= nvar_; n++) {
            qleft[n]  = qm[IDX6(g, i-1, j, k, n, xdim, ngrid)];
            qright[n] = qp[IDX6(g, i, j, k, n, xdim, ngrid)];
        }

        if (solver == RIEMANN_HLLC) riemann_hllc_device(qleft, qright, fgdnv, nvar_);
        else if (solver == RIEMANN_HLL) riemann_hll_device(qleft, qright, fgdnv, nvar_);
        else riemann_llf_device(qleft, qright, fgdnv, nvar_);

        flux[IDX_FLUX(g, i, j, k, 1, 1, ngrid)]           = fgdnv[1] * dtdx;
        flux[IDX_FLUX(g, i, j, k, ln, 1, ngrid)]           = fgdnv[2] * dtdx;
        flux[IDX_FLUX(g, i, j, k, lt1, 1, ngrid)]          = fgdnv[4] * dtdx;
        flux[IDX_FLUX(g, i, j, k, lt2, 1, ngrid)]          = fgdnv[5] * dtdx;
        flux[IDX_FLUX(g, i, j, k, ndim_ + 2, 1, ngrid)]   = fgdnv[3] * dtdx;
        for (int n = ndim_ + 3; n <= nvar_; n++) {
            flux[IDX_FLUX(g, i, j, k, n, 1, ngrid)] = fgdnv[n] * dtdx;
        }
        tmp[IDX_TMP(g, i, j, k, 1, 1, ngrid)] = 0.5 * (qleft[2] + qright[2]) * dtdx;
        tmp[IDX_TMP(g, i, j, k, 2, 1, ngrid)] = fgdnv[nvar_ + 1] * dtdx;
    }

    // ---- Y-direction flux ----
    if (i >= FLUX_ILO && i <= FLUX_IHI && j >= JF1 && j <= JF2 && k >= FLUX_KLO && k <= FLUX_KHI) {
        const int xdim = 2;
        const int ln = 3, lt1 = 2, lt2 = 4;

        qleft[1] = qm[IDX6(g, i, j-1, k, 1,   xdim, ngrid)];
        qleft[2] = qm[IDX6(g, i, j-1, k, ln,  xdim, ngrid)];
        qleft[3] = qm[IDX6(g, i, j-1, k, ndim_ + 2, xdim, ngrid)];
        qleft[4] = qm[IDX6(g, i, j-1, k, lt1, xdim, ngrid)];
        qleft[5] = qm[IDX6(g, i, j-1, k, lt2, xdim, ngrid)];

        qright[1] = qp[IDX6(g, i, j, k, 1,   xdim, ngrid)];
        qright[2] = qp[IDX6(g, i, j, k, ln,  xdim, ngrid)];
        qright[3] = qp[IDX6(g, i, j, k, ndim_ + 2, xdim, ngrid)];
        qright[4] = qp[IDX6(g, i, j, k, lt1, xdim, ngrid)];
        qright[5] = qp[IDX6(g, i, j, k, lt2, xdim, ngrid)];

        for (int n = ndim_ + 3; n <= nvar_; n++) {
            qleft[n]  = qm[IDX6(g, i, j-1, k, n, xdim, ngrid)];
            qright[n] = qp[IDX6(g, i, j, k, n, xdim, ngrid)];
        }

        if (solver == RIEMANN_HLLC) riemann_hllc_device(qleft, qright, fgdnv, nvar_);
        else if (solver == RIEMANN_HLL) riemann_hll_device(qleft, qright, fgdnv, nvar_);
        else riemann_llf_device(qleft, qright, fgdnv, nvar_);

        flux[IDX_FLUX(g, i, j, k, 1, 2, ngrid)]           = fgdnv[1] * dtdy;
        flux[IDX_FLUX(g, i, j, k, ln, 2, ngrid)]           = fgdnv[2] * dtdy;
        flux[IDX_FLUX(g, i, j, k, lt1, 2, ngrid)]          = fgdnv[4] * dtdy;
        flux[IDX_FLUX(g, i, j, k, lt2, 2, ngrid)]          = fgdnv[5] * dtdy;
        flux[IDX_FLUX(g, i, j, k, ndim_ + 2, 2, ngrid)]   = fgdnv[3] * dtdy;
        for (int n = ndim_ + 3; n <= nvar_; n++) {
            flux[IDX_FLUX(g, i, j, k, n, 2, ngrid)] = fgdnv[n] * dtdy;
        }
        tmp[IDX_TMP(g, i, j, k, 1, 2, ngrid)] = 0.5 * (qleft[2] + qright[2]) * dtdy;
        tmp[IDX_TMP(g, i, j, k, 2, 2, ngrid)] = fgdnv[nvar_ + 1] * dtdy;
    }

    // ---- Z-direction flux ----
    if (i >= FLUX_ILO && i <= FLUX_IHI && j >= FLUX_JLO && j <= FLUX_JHI && k >= KF1 && k <= KF2) {
        const int xdim = 3;
        const int ln = 4, lt1 = 2, lt2 = 3;

        qleft[1] = qm[IDX6(g, i, j, k-1, 1,   xdim, ngrid)];
        qleft[2] = qm[IDX6(g, i, j, k-1, ln,  xdim, ngrid)];
        qleft[3] = qm[IDX6(g, i, j, k-1, ndim_ + 2, xdim, ngrid)];
        qleft[4] = qm[IDX6(g, i, j, k-1, lt1, xdim, ngrid)];
        qleft[5] = qm[IDX6(g, i, j, k-1, lt2, xdim, ngrid)];

        qright[1] = qp[IDX6(g, i, j, k, 1,   xdim, ngrid)];
        qright[2] = qp[IDX6(g, i, j, k, ln,  xdim, ngrid)];
        qright[3] = qp[IDX6(g, i, j, k, ndim_ + 2, xdim, ngrid)];
        qright[4] = qp[IDX6(g, i, j, k, lt1, xdim, ngrid)];
        qright[5] = qp[IDX6(g, i, j, k, lt2, xdim, ngrid)];

        for (int n = ndim_ + 3; n <= nvar_; n++) {
            qleft[n]  = qm[IDX6(g, i, j, k-1, n, xdim, ngrid)];
            qright[n] = qp[IDX6(g, i, j, k, n, xdim, ngrid)];
        }

        if (solver == RIEMANN_HLLC) riemann_hllc_device(qleft, qright, fgdnv, nvar_);
        else if (solver == RIEMANN_HLL) riemann_hll_device(qleft, qright, fgdnv, nvar_);
        else riemann_llf_device(qleft, qright, fgdnv, nvar_);

        flux[IDX_FLUX(g, i, j, k, 1, 3, ngrid)]           = fgdnv[1] * dtdz;
        flux[IDX_FLUX(g, i, j, k, ln, 3, ngrid)]           = fgdnv[2] * dtdz;
        flux[IDX_FLUX(g, i, j, k, lt1, 3, ngrid)]          = fgdnv[4] * dtdz;
        flux[IDX_FLUX(g, i, j, k, lt2, 3, ngrid)]          = fgdnv[5] * dtdz;
        flux[IDX_FLUX(g, i, j, k, ndim_ + 2, 3, ngrid)]   = fgdnv[3] * dtdz;
        for (int n = ndim_ + 3; n <= nvar_; n++) {
            flux[IDX_FLUX(g, i, j, k, n, 3, ngrid)] = fgdnv[n] * dtdz;
        }
        tmp[IDX_TMP(g, i, j, k, 1, 3, ngrid)] = 0.5 * (qleft[2] + qright[2]) * dtdz;
        tmp[IDX_TMP(g, i, j, k, 2, 3, ngrid)] = fgdnv[nvar_ + 1] * dtdz;
    }
}

// ====================================================================
// Kernel 5: cmpdivu + consup  (artificial viscosity)
// Only adds diffusive fluxes when difmag > 0.
// ====================================================================
__global__ void hydro_difmag_kernel(
    const double* __restrict__ q,
    const double* __restrict__ uloc,
    double* __restrict__ flux,
    double dx, double dy, double dz, double dt,
    int ngrid)
{
    const int g = blockIdx.x;
    if (g >= ngrid) return;
    const int tid = threadIdx.x;
    if (tid >= NCELLS_PER_GRID) return;

    int i, j, k;
    thread_to_ijk(tid, i, j, k);

    const double difmag_ = d_hp.difmag;
    const int    nvar_    = d_hp.nvar;

    double factorx = 0.25 / dx;
    double factory = 0.25 / dy;
    double factorz = 0.25 / dz;

    // Shared memory for divergence: 3x3x3 per grid block
    __shared__ double s_div[FLUX_NI * FLUX_NJ * FLUX_NK];

    // Compute divergence for flux-region vertices
    if (i >= IF1 && i <= IF2 && j >= JF1 && j <= JF2 && k >= KF1 && k <= KF2) {
        int div_idx = (i - IF1) + FLUX_NI * ((j - JF1) + FLUX_NJ * (k - KF1));

        double ux = factorx * (
            (q[IDX5(g, i, j, k, 2, ngrid)]     - q[IDX5(g, i-1, j, k, 2, ngrid)])
          + (q[IDX5(g, i, j-1, k, 2, ngrid)]   - q[IDX5(g, i-1, j-1, k, 2, ngrid)])
          + (q[IDX5(g, i, j, k-1, 2, ngrid)]   - q[IDX5(g, i-1, j, k-1, 2, ngrid)])
          + (q[IDX5(g, i, j-1, k-1, 2, ngrid)] - q[IDX5(g, i-1, j-1, k-1, 2, ngrid)])
        );

        double vy = factory * (
            (q[IDX5(g, i, j, k, 3, ngrid)]     - q[IDX5(g, i, j-1, k, 3, ngrid)])
          + (q[IDX5(g, i-1, j, k, 3, ngrid)]   - q[IDX5(g, i-1, j-1, k, 3, ngrid)])
          + (q[IDX5(g, i, j, k-1, 3, ngrid)]   - q[IDX5(g, i, j-1, k-1, 3, ngrid)])
          + (q[IDX5(g, i-1, j, k-1, 3, ngrid)] - q[IDX5(g, i-1, j-1, k-1, 3, ngrid)])
        );

        double wz = factorz * (
            (q[IDX5(g, i, j, k, 4, ngrid)]     - q[IDX5(g, i, j, k-1, 4, ngrid)])
          + (q[IDX5(g, i, j-1, k, 4, ngrid)]   - q[IDX5(g, i, j-1, k-1, 4, ngrid)])
          + (q[IDX5(g, i-1, j, k, 4, ngrid)]   - q[IDX5(g, i-1, j, k-1, 4, ngrid)])
          + (q[IDX5(g, i-1, j-1, k, 4, ngrid)] - q[IDX5(g, i-1, j-1, k-1, 4, ngrid)])
        );

        s_div[div_idx] = ux + vy + wz;
    }
    __syncthreads();

    #define SDIV(ii, jj, kk) s_div[((ii)-IF1) + FLUX_NI * (((jj)-JF1) + FLUX_NJ * ((kk)-KF1))]

    double factor = 0.25;  // half^(ndim-1) for ndim=3

    // X-direction faces
    if (i >= IF1 && i <= IF2 && j >= JF1 && j <= (JF2-1) && k >= KF1 && k <= (KF2-1)) {
        double div1 = factor * (SDIV(i, j, k) + SDIV(i, j+1, k) + SDIV(i, j, k+1) + SDIV(i, j+1, k+1));
        div1 = difmag_ * fmin(0.0, div1);
        for (int n = 1; n <= nvar_; n++) {
            flux[IDX_FLUX(g, i, j, k, n, 1, ngrid)] +=
                dt * div1 * (uloc[IDX5(g, i, j, k, n, ngrid)] - uloc[IDX5(g, i-1, j, k, n, ngrid)]);
        }
    }

    // Y-direction faces
    if (i >= 1 && i <= 2 && j >= JF1 && j <= JF2 && k >= KF1 && k <= (KF2-1)) {
        double div1 = factor * (SDIV(i, j, k) + SDIV(i+1, j, k) + SDIV(i, j, k+1) + SDIV(i+1, j, k+1));
        div1 = difmag_ * fmin(0.0, div1);
        for (int n = 1; n <= nvar_; n++) {
            flux[IDX_FLUX(g, i, j, k, n, 2, ngrid)] +=
                dt * div1 * (uloc[IDX5(g, i, j, k, n, ngrid)] - uloc[IDX5(g, i, j-1, k, n, ngrid)]);
        }
    }

    // Z-direction faces
    if (i >= 1 && i <= 2 && j >= 1 && j <= 2 && k >= KF1 && k <= KF2) {
        double div1 = factor * (SDIV(i, j, k) + SDIV(i+1, j, k) + SDIV(i, j+1, k) + SDIV(i+1, j+1, k));
        div1 = difmag_ * fmin(0.0, div1);
        for (int n = 1; n <= nvar_; n++) {
            flux[IDX_FLUX(g, i, j, k, n, 3, ngrid)] +=
                dt * div1 * (uloc[IDX5(g, i, j, k, n, ngrid)] - uloc[IDX5(g, i, j, k-1, n, ngrid)]);
        }
    }

    #undef SDIV
}

// ====================================================================
// Kernel 6: GPU Gather -- fills d_uloc/d_gloc from mesh + stencil indices
// One block per grid, 216 threads per block.
// stencil_idx > 0 => direct cell access into mesh_uold
// stencil_idx < 0 => interpolated value from interp_vals (AoS layout)
// stencil_grav > 0 => cell index into mesh_f; == 0 => zero gravity
// ====================================================================
__global__ void hydro_gather_kernel(
    const int* __restrict__ stencil_idx,
    const int* __restrict__ stencil_grav,
    const double* __restrict__ interp_vals,
    const double* __restrict__ mesh_uold,
    const double* __restrict__ mesh_f,
    double* __restrict__ uloc,
    double* __restrict__ gloc,
    long long ncell,
    int ngrid,
    int stride)
{
    const int g = blockIdx.x;
    if (g >= ngrid) return;
    const int tid = threadIdx.x;
    if (tid >= NCELLS_PER_GRID) return;

    int i, j, k;
    thread_to_ijk(tid, i, j, k);

    const int nvar_ = d_hp.nvar;
    const int ndim_ = d_hp.ndim;

    // Hydro data gather
    int sidx = stencil_idx[IDX4(g, i, j, k, stride)];
    if (sidx > 0) {
        // Direct cell: sidx is 1-based Fortran index
        long long cell = (long long)(sidx - 1);
        for (int n = 1; n <= nvar_; n++) {
            uloc[IDX5(g, i, j, k, n, stride)] =
                mesh_uold[cell + ncell * (long long)(n - 1)];
        }
    } else if (sidx < 0) {
        // Interpolated: -sidx is 1-based slot, AoS layout: base = (slot-1)*NVAR
        int base = (-sidx - 1) * NVAR;
        for (int n = 1; n <= nvar_; n++) {
            uloc[IDX5(g, i, j, k, n, stride)] = interp_vals[base + n - 1];
        }
    } else {
        for (int n = 1; n <= nvar_; n++) {
            uloc[IDX5(g, i, j, k, n, stride)] = 0.0;
        }
    }

    // Gravity gather
    int gidx = stencil_grav[IDX4(g, i, j, k, stride)];
    if (gidx > 0) {
        long long cell = (long long)(gidx - 1);
        for (int d = 1; d <= ndim_; d++) {
            gloc[IDX5G(g, i, j, k, d, stride)] =
                mesh_f[cell + ncell * (long long)(d - 1)];
        }
    } else {
        for (int d = 1; d <= ndim_; d++) {
            gloc[IDX5G(g, i, j, k, d, stride)] = 0.0;
        }
    }
}

// ====================================================================
// Host-callable functions (C API for Fortran ISO_C_BINDING)
// ====================================================================
// ====================================================================
// Kernel 6: scatter_reduce
// Combines flux reset + Level L scatter + Level L-1 scatter into compact
// output arrays, eliminating 98 MB D2H transfer of raw flux/tmp.
//
// Thread assignment (NVAR=11):
//   tid 0..87   : Level L (8 children × 11 vars)
//   tid 88..153 : Level L-1 (6 faces × 11 vars)
//   tid 154..161: pressure_fix Level L divu+enew (8 children)
//   tid 162..167: pressure_fix Level L-1 divu+enew (6 faces)
//   tid 0..215  : shared memory loading (ok array)
// Block size: 256 threads (8 warps)
// ====================================================================
__global__ void hydro_scatter_reduce_kernel(
    const double* __restrict__ flux,
    const double* __restrict__ tmp,
    const int* __restrict__ ok_arr,
    double* __restrict__ add_unew,
    double* __restrict__ add_lm1,
    double* __restrict__ add_divu_l,
    double* __restrict__ add_enew_l,
    double* __restrict__ add_divu_lm1,
    double* __restrict__ add_enew_lm1,
    double oneontwotondim,
    int ngrid, int stride)
{
    const int g = blockIdx.x;
    if (g >= ngrid) return;
    const int tid = threadIdx.x;

    const int nvar_ = d_hp.nvar;
    const int pfix  = d_hp.pressure_fix;

    // Load ok array into shared memory (216 ints = 864 bytes)
    __shared__ int s_ok[STENCIL_NI * STENCIL_NJ * STENCIL_NK];
    if (tid < NCELLS_PER_GRID) {
        s_ok[tid] = ok_arr[g + (long)stride * (long)tid];
    }
    __syncthreads();

    // ok(i,j,k) from shared memory (stencil indices -1..4)
    #define LOK(i,j,k) s_ok[((i)-IU1) + STENCIL_NI*(((j)-JU1) + STENCIL_NJ*((k)-KU1))]

    // Inline flux read with flux-reset check:
    // flux at position (i3,j3,k3,ivar,idim) is zero if
    //   ok(i3-i0, j3-j0, k3-k0) || ok(i3, j3, k3)
    // where (i0,j0,k0) is the direction offset for idim.

    // === Level L: threads 0..8*nvar-1 ===
    if (tid < 8 * nvar_) {
        int child = tid / nvar_;      // 0..7 (ind_son-1)
        int ivar  = tid % nvar_ + 1;  // 1..nvar

        int i2 = child & 1;
        int j2 = (child >> 1) & 1;
        int k2 = (child >> 2) & 1;
        int i3 = 1 + i2, j3 = 1 + j2, k3 = 1 + k2;

        double acc = 0.0;

        // idim=1 (i0=1)
        {
            double fl = (LOK(i3-1,j3,k3) || LOK(i3,j3,k3)) ? 0.0
                        : flux[IDX_FLUX(g, i3, j3, k3, ivar, 1, stride)];
            double fr = (LOK(i3,j3,k3) || LOK(i3+1,j3,k3)) ? 0.0
                        : flux[IDX_FLUX(g, i3+1, j3, k3, ivar, 1, stride)];
            acc += fl - fr;
        }
        // idim=2 (j0=1)
        {
            double fl = (LOK(i3,j3-1,k3) || LOK(i3,j3,k3)) ? 0.0
                        : flux[IDX_FLUX(g, i3, j3, k3, ivar, 2, stride)];
            double fr = (LOK(i3,j3,k3) || LOK(i3,j3+1,k3)) ? 0.0
                        : flux[IDX_FLUX(g, i3, j3+1, k3, ivar, 2, stride)];
            acc += fl - fr;
        }
        // idim=3 (k0=1)
        {
            double fl = (LOK(i3,j3,k3-1) || LOK(i3,j3,k3)) ? 0.0
                        : flux[IDX_FLUX(g, i3, j3, k3, ivar, 3, stride)];
            double fr = (LOK(i3,j3,k3) || LOK(i3,j3,k3+1)) ? 0.0
                        : flux[IDX_FLUX(g, i3, j3, k3+1, ivar, 3, stride)];
            acc += fl - fr;
        }

        add_unew[g + (long)stride * (child + 8 * (ivar - 1))] = acc;
    }

    // === Level L-1: threads 8*nvar..14*nvar-1 ===
    {
        int lm1_tid = tid - 8 * nvar_;
        if (lm1_tid >= 0 && lm1_tid < 6 * nvar_) {
            int face = lm1_tid / nvar_;     // 0..5
            int ivar = lm1_tid % nvar_ + 1; // 1..nvar

            double acc = 0.0;

            // face 0: -x  flux(g, 1, j, k, ivar, 1),  j∈{1,2}, k∈{1,2}
            // face 1: +x  flux(g, 3, j, k, ivar, 1)
            // face 2: -y  flux(g, i, 1, k, ivar, 2),  i∈{1,2}, k∈{1,2}
            // face 3: +y  flux(g, i, 3, k, ivar, 2)
            // face 4: -z  flux(g, i, j, 1, ivar, 3),  i∈{1,2}, j∈{1,2}
            // face 5: +z  flux(g, i, j, 3, ivar, 3)
            switch (face) {
            case 0:
                for (int kk=1; kk<=2; kk++)
                for (int jj=1; jj<=2; jj++) {
                    double f = (LOK(0,jj,kk) || LOK(1,jj,kk)) ? 0.0
                               : flux[IDX_FLUX(g, 1, jj, kk, ivar, 1, stride)];
                    acc -= f * oneontwotondim;
                } break;
            case 1:
                for (int kk=1; kk<=2; kk++)
                for (int jj=1; jj<=2; jj++) {
                    double f = (LOK(2,jj,kk) || LOK(3,jj,kk)) ? 0.0
                               : flux[IDX_FLUX(g, 3, jj, kk, ivar, 1, stride)];
                    acc += f * oneontwotondim;
                } break;
            case 2:
                for (int kk=1; kk<=2; kk++)
                for (int ii=1; ii<=2; ii++) {
                    double f = (LOK(ii,0,kk) || LOK(ii,1,kk)) ? 0.0
                               : flux[IDX_FLUX(g, ii, 1, kk, ivar, 2, stride)];
                    acc -= f * oneontwotondim;
                } break;
            case 3:
                for (int kk=1; kk<=2; kk++)
                for (int ii=1; ii<=2; ii++) {
                    double f = (LOK(ii,2,kk) || LOK(ii,3,kk)) ? 0.0
                               : flux[IDX_FLUX(g, ii, 3, kk, ivar, 2, stride)];
                    acc += f * oneontwotondim;
                } break;
            case 4:
                for (int jj=1; jj<=2; jj++)
                for (int ii=1; ii<=2; ii++) {
                    double f = (LOK(ii,jj,0) || LOK(ii,jj,1)) ? 0.0
                               : flux[IDX_FLUX(g, ii, jj, 1, ivar, 3, stride)];
                    acc -= f * oneontwotondim;
                } break;
            case 5:
                for (int jj=1; jj<=2; jj++)
                for (int ii=1; ii<=2; ii++) {
                    double f = (LOK(ii,jj,2) || LOK(ii,jj,3)) ? 0.0
                               : flux[IDX_FLUX(g, ii, jj, 3, ivar, 3, stride)];
                    acc += f * oneontwotondim;
                } break;
            }

            add_lm1[g + (long)stride * (face + 6 * (ivar - 1))] = acc;
        }
    }

    // === pressure_fix: divu and enew ===
    if (pfix) {
        // Level L: threads 14*nvar..14*nvar+7
        int pfix_l = tid - 14 * nvar_;
        if (pfix_l >= 0 && pfix_l < 8) {
            int child = pfix_l;
            int i2 = child & 1;
            int j2 = (child >> 1) & 1;
            int k2 = (child >> 2) & 1;
            int i3 = 1+i2, j3 = 1+j2, k3 = 1+k2;

            double acc_d = 0.0, acc_e = 0.0;
            // idim=1
            {
                int rl = LOK(i3-1,j3,k3) || LOK(i3,j3,k3);
                int rr = LOK(i3,j3,k3) || LOK(i3+1,j3,k3);
                acc_d += (rl?0.0:tmp[IDX_TMP(g,i3,j3,k3,1,1,stride)])
                       - (rr?0.0:tmp[IDX_TMP(g,i3+1,j3,k3,1,1,stride)]);
                acc_e += (rl?0.0:tmp[IDX_TMP(g,i3,j3,k3,2,1,stride)])
                       - (rr?0.0:tmp[IDX_TMP(g,i3+1,j3,k3,2,1,stride)]);
            }
            // idim=2
            {
                int rl = LOK(i3,j3-1,k3) || LOK(i3,j3,k3);
                int rr = LOK(i3,j3,k3) || LOK(i3,j3+1,k3);
                acc_d += (rl?0.0:tmp[IDX_TMP(g,i3,j3,k3,1,2,stride)])
                       - (rr?0.0:tmp[IDX_TMP(g,i3,j3+1,k3,1,2,stride)]);
                acc_e += (rl?0.0:tmp[IDX_TMP(g,i3,j3,k3,2,2,stride)])
                       - (rr?0.0:tmp[IDX_TMP(g,i3,j3+1,k3,2,2,stride)]);
            }
            // idim=3
            {
                int rl = LOK(i3,j3,k3-1) || LOK(i3,j3,k3);
                int rr = LOK(i3,j3,k3) || LOK(i3,j3,k3+1);
                acc_d += (rl?0.0:tmp[IDX_TMP(g,i3,j3,k3,1,3,stride)])
                       - (rr?0.0:tmp[IDX_TMP(g,i3,j3,k3+1,1,3,stride)]);
                acc_e += (rl?0.0:tmp[IDX_TMP(g,i3,j3,k3,2,3,stride)])
                       - (rr?0.0:tmp[IDX_TMP(g,i3,j3,k3+1,2,3,stride)]);
            }

            add_divu_l[g + (long)stride * child] = acc_d;
            add_enew_l[g + (long)stride * child] = acc_e;
        }

        // Level L-1: threads 14*nvar+8..14*nvar+13
        int pfix_lm1 = tid - 14 * nvar_ - 8;
        if (pfix_lm1 >= 0 && pfix_lm1 < 6) {
            int face = pfix_lm1;
            double acc_d = 0.0, acc_e = 0.0;

            switch (face) {
            case 0:
                for (int kk=1; kk<=2; kk++)
                for (int jj=1; jj<=2; jj++) {
                    if (!(LOK(0,jj,kk)||LOK(1,jj,kk))) {
                        acc_d -= tmp[IDX_TMP(g,1,jj,kk,1,1,stride)] * oneontwotondim;
                        acc_e -= tmp[IDX_TMP(g,1,jj,kk,2,1,stride)] * oneontwotondim;
                    }
                } break;
            case 1:
                for (int kk=1; kk<=2; kk++)
                for (int jj=1; jj<=2; jj++) {
                    if (!(LOK(2,jj,kk)||LOK(3,jj,kk))) {
                        acc_d += tmp[IDX_TMP(g,3,jj,kk,1,1,stride)] * oneontwotondim;
                        acc_e += tmp[IDX_TMP(g,3,jj,kk,2,1,stride)] * oneontwotondim;
                    }
                } break;
            case 2:
                for (int kk=1; kk<=2; kk++)
                for (int ii=1; ii<=2; ii++) {
                    if (!(LOK(ii,0,kk)||LOK(ii,1,kk))) {
                        acc_d -= tmp[IDX_TMP(g,ii,1,kk,1,2,stride)] * oneontwotondim;
                        acc_e -= tmp[IDX_TMP(g,ii,1,kk,2,2,stride)] * oneontwotondim;
                    }
                } break;
            case 3:
                for (int kk=1; kk<=2; kk++)
                for (int ii=1; ii<=2; ii++) {
                    if (!(LOK(ii,2,kk)||LOK(ii,3,kk))) {
                        acc_d += tmp[IDX_TMP(g,ii,3,kk,1,2,stride)] * oneontwotondim;
                        acc_e += tmp[IDX_TMP(g,ii,3,kk,2,2,stride)] * oneontwotondim;
                    }
                } break;
            case 4:
                for (int jj=1; jj<=2; jj++)
                for (int ii=1; ii<=2; ii++) {
                    if (!(LOK(ii,jj,0)||LOK(ii,jj,1))) {
                        acc_d -= tmp[IDX_TMP(g,ii,jj,1,1,3,stride)] * oneontwotondim;
                        acc_e -= tmp[IDX_TMP(g,ii,jj,1,2,3,stride)] * oneontwotondim;
                    }
                } break;
            case 5:
                for (int jj=1; jj<=2; jj++)
                for (int ii=1; ii<=2; ii++) {
                    if (!(LOK(ii,jj,2)||LOK(ii,jj,3))) {
                        acc_d += tmp[IDX_TMP(g,ii,jj,3,1,3,stride)] * oneontwotondim;
                        acc_e += tmp[IDX_TMP(g,ii,jj,3,2,3,stride)] * oneontwotondim;
                    }
                } break;
            }

            add_divu_lm1[g + (long)stride * face] = acc_d;
            add_enew_lm1[g + (long)stride * face] = acc_e;
        }
    }

    #undef LOK
}

// ====================================================================
// Kernel 7: GPU-gather ok flags from mesh_son + stencil_idx
// Computes ok_int for scatter_reduce_kernel from GPU-resident son array.
// ok = 1 if stencil_idx > 0 (fine cell) AND mesh_son[cell] > 0 (has children)
// ok = 0 otherwise (interpolated or leaf cell)
// Layout: ok_arr[g + stride * tid] matches scatter_reduce_kernel's s_ok loading
// ====================================================================
__global__ void hydro_gather_ok_kernel(
    const int* __restrict__ stencil_idx,
    const int* __restrict__ mesh_son,
    int* __restrict__ ok_arr,
    int ngrid, int stride)
{
    const int g = blockIdx.x;
    if (g >= ngrid) return;
    const int tid = threadIdx.x;
    if (tid >= NCELLS_PER_GRID) return;

    int sidx = stencil_idx[g + (long)stride * (long)tid];
    int ok = 0;
    if (sidx > 0) {
        ok = (mesh_son[sidx - 1] > 0) ? 1 : 0;
    }
    ok_arr[g + (long)stride * (long)tid] = ok;
}

extern "C" {

void hydro_cuda_init(
    double gamma_, double smallr_, double smallc_,
    int slope_type_, double slope_theta_,
    int nvar_, int ndim_, int riemann_solver_,
    int pressure_fix_, double difmag_)
{
    HydroParams hp;
    hp.gamma          = gamma_;
    hp.smallr         = smallr_;
    hp.smallc         = smallc_;
    hp.smallp         = smallc_ * smallc_ / gamma_;
    hp.smalle         = smallc_ * smallc_ / gamma_ / (gamma_ - 1.0);
    hp.entho          = 1.0 / (gamma_ - 1.0);
    hp.slope_type     = slope_type_;
    hp.slope_theta    = slope_theta_;
    hp.nvar           = nvar_;
    hp.ndim           = ndim_;
    hp.riemann_solver = riemann_solver_;
    hp.pressure_fix   = pressure_fix_;
    hp.difmag         = difmag_;

    cudaMemcpyToSymbol(d_hp, &hp, sizeof(HydroParams));
}

void hydro_cuda_unsplit_async(
    const double* h_uloc,
    const double* h_gloc,
    double* h_flux,
    double* h_tmp,
    const int* h_ok,
    double dx, double dy, double dz, double dt,
    int ngrid,
    int stride,
    int stream_slot)
{

    // Ensure device buffers are large enough (use stride, not ngrid)
    pool_ensure_hydro_buffers(stream_slot, stride);
    pool_ensure_hydro_inter_buffers(stream_slot, stride);

    StreamSlot* s = &get_pool()[stream_slot];
    cudaStream_t strm = s->stream;

    // Sizes: always use NVECTOR stride to match Fortran array layout
    size_t uloc_bytes = (size_t)stride * STENCIL_NI * STENCIL_NJ * STENCIL_NK * NVAR * sizeof(double);
    size_t gloc_bytes = (size_t)stride * STENCIL_NI * STENCIL_NJ * STENCIL_NK * NDIM * sizeof(double);
    size_t flux_bytes = (size_t)stride * FLUX_NI * FLUX_NJ * FLUX_NK * NVAR * NDIM * sizeof(double);
    size_t tmp_bytes  = (size_t)stride * FLUX_NI * FLUX_NJ * FLUX_NK * 2 * NDIM * sizeof(double);

    // Profile events: [0]H2D[1]ctoprim[2]uslope[3]trace3d[4]flux[5]difmag[6]D2H[7]
    cudaEvent_t* ev = s->ev_profile;

    // H -> D: upload input arrays (full NVECTOR-strided arrays)
    cudaEventRecord(ev[0], strm);
    cudaMemcpyAsync(s->d_uloc, h_uloc, uloc_bytes, cudaMemcpyHostToDevice, strm);
    cudaMemcpyAsync(s->d_gloc, h_gloc, gloc_bytes, cudaMemcpyHostToDevice, strm);
    cudaEventRecord(ev[1], strm);

    // Kernel launch: one block per ACTIVE grid, 216 threads per block.
    // Grids g=0..ngrid-1 are active; g=ngrid..stride-1 are unused padding.
    // The stride parameter tells kernels the first-dimension stride.
    int nblocks  = ngrid;
    int nthreads = NCELLS_PER_GRID;  // 216

    double dtdx = dt / dx;
    double dtdy = dt / dy;
    double dtdz = dt / dz;

    // Kernel 1: Conservative to primitive
    hydro_ctoprim_kernel<<<nblocks, nthreads, 0, strm>>>(
        s->d_uloc, s->d_gloc, s->d_q, s->d_c, dt, stride);
    cudaEventRecord(ev[2], strm);

    // Kernel 2: TVD slopes
    hydro_uslope_kernel<<<nblocks, nthreads, 0, strm>>>(
        s->d_q, s->d_dq, stride);
    cudaEventRecord(ev[3], strm);

    // Kernel 3: MUSCL-Hancock trace
    hydro_trace3d_kernel<<<nblocks, nthreads, 0, strm>>>(
        s->d_q, s->d_dq, s->d_qm, s->d_qp, dtdx, dtdy, dtdz, stride);
    cudaEventRecord(ev[4], strm);

    // Zero flux and tmp output arrays
    cudaMemsetAsync(s->d_flux, 0, flux_bytes, strm);
    cudaMemsetAsync(s->d_tmp,  0, tmp_bytes,  strm);

    // Kernel 4: Riemann solves -> fluxes
    hydro_flux_kernel<<<nblocks, nthreads, 0, strm>>>(
        s->d_qm, s->d_qp, s->d_flux, s->d_tmp, dtdx, dtdy, dtdz, stride);
    cudaEventRecord(ev[5], strm);

    // Kernel 5: Artificial viscosity (difmag kernel)
    hydro_difmag_kernel<<<nblocks, nthreads,
        FLUX_NI * FLUX_NJ * FLUX_NK * sizeof(double), strm>>>(
        s->d_q, s->d_uloc, s->d_flux, dx, dy, dz, dt, stride);
    cudaEventRecord(ev[6], strm);

    // D -> H: download results (full NVECTOR-strided arrays)
    cudaMemcpyAsync(h_flux, s->d_flux, flux_bytes, cudaMemcpyDeviceToHost, strm);
    cudaMemcpyAsync(h_tmp,  s->d_tmp,  tmp_bytes,  cudaMemcpyDeviceToHost, strm);
    cudaEventRecord(ev[7], strm);
}

void hydro_cuda_unsplit_sync(
    double* h_flux,
    double* h_tmp,
    int ngrid,
    int stream_slot)
{
    StreamSlot* s = &get_pool()[stream_slot];
    cudaStreamSynchronize(s->stream);
    cudaError_t err = cudaGetLastError();
    if (err != cudaSuccess) {
        fprintf(stderr, "CUDA hydro error (slot %d, ngrid %d): %s\n",
                stream_slot, ngrid, cudaGetErrorString(err));
    }
}

// ====================================================================
// GPU-gather variant: upload stencil indices, gather on GPU, compute
// Eliminates 97 MB/chunk uloc/gloc PCIe transfer
// ====================================================================
void hydro_cuda_gather_unsplit_async(
    const int* h_stencil_idx,
    const int* h_stencil_grav,
    const double* h_interp_vals,
    double* h_flux,
    double* h_tmp,
    double dx, double dy, double dz, double dt,
    int ngrid,
    int stride,
    int n_interp,
    int stream_slot)
{
    pool_ensure_hydro_buffers(stream_slot, stride);
    pool_ensure_hydro_inter_buffers(stream_slot, stride);
    pool_ensure_stencil_buffers(stream_slot, stride, n_interp);

    StreamSlot* s = &get_pool()[stream_slot];
    cudaStream_t strm = s->stream;
    long long ncell = cuda_get_mesh_ncell();

    // Profile events: [0]H2D[1]gather+ctoprim[2]uslope[3]trace3d[4]flux[5]difmag[6]D2H[7]
    cudaEvent_t* ev = s->ev_profile;

    // Upload stencil indices (small: ~7 MB for stride=4096)
    size_t idx_bytes = (size_t)stride * STENCIL_NI * STENCIL_NJ * STENCIL_NK * sizeof(int);
    cudaEventRecord(ev[0], strm);
    cudaMemcpyAsync(s->d_stencil_idx,  h_stencil_idx,  idx_bytes, cudaMemcpyHostToDevice, strm);
    cudaMemcpyAsync(s->d_stencil_grav, h_stencil_grav, idx_bytes, cudaMemcpyHostToDevice, strm);

    // Upload interpolation values (AoS: n_interp * NVAR doubles)
    if (n_interp > 0 && h_interp_vals) {
        size_t interp_bytes = (size_t)n_interp * NVAR * sizeof(double);
        cudaMemcpyAsync(s->d_interp_vals, h_interp_vals, interp_bytes, cudaMemcpyHostToDevice, strm);
    }
    cudaEventRecord(ev[1], strm);

    int nblocks  = ngrid;
    int nthreads = NCELLS_PER_GRID;

    // Kernel 0: GPU gather (fills d_uloc and d_gloc from mesh)
    hydro_gather_kernel<<<nblocks, nthreads, 0, strm>>>(
        s->d_stencil_idx, s->d_stencil_grav, s->d_interp_vals,
        cuda_get_mesh_uold(), cuda_get_mesh_f(),
        s->d_uloc, s->d_gloc,
        ncell, ngrid, stride);

    double dtdx = dt / dx;
    double dtdy = dt / dy;
    double dtdz = dt / dz;

    // Kernel 1: Conservative to primitive
    hydro_ctoprim_kernel<<<nblocks, nthreads, 0, strm>>>(
        s->d_uloc, s->d_gloc, s->d_q, s->d_c, dt, stride);
    cudaEventRecord(ev[2], strm);

    // Kernel 2: TVD slopes
    hydro_uslope_kernel<<<nblocks, nthreads, 0, strm>>>(
        s->d_q, s->d_dq, stride);
    cudaEventRecord(ev[3], strm);

    // Kernel 3: MUSCL-Hancock trace
    hydro_trace3d_kernel<<<nblocks, nthreads, 0, strm>>>(
        s->d_q, s->d_dq, s->d_qm, s->d_qp, dtdx, dtdy, dtdz, stride);
    cudaEventRecord(ev[4], strm);

    // Zero flux and tmp output arrays
    size_t flux_bytes = (size_t)stride * FLUX_NI * FLUX_NJ * FLUX_NK * NVAR * NDIM * sizeof(double);
    size_t tmp_bytes  = (size_t)stride * FLUX_NI * FLUX_NJ * FLUX_NK * 2 * NDIM * sizeof(double);
    cudaMemsetAsync(s->d_flux, 0, flux_bytes, strm);
    cudaMemsetAsync(s->d_tmp,  0, tmp_bytes,  strm);

    // Kernel 4: Riemann solves -> fluxes
    hydro_flux_kernel<<<nblocks, nthreads, 0, strm>>>(
        s->d_qm, s->d_qp, s->d_flux, s->d_tmp, dtdx, dtdy, dtdz, stride);
    cudaEventRecord(ev[5], strm);

    // Kernel 5: Artificial viscosity
    hydro_difmag_kernel<<<nblocks, nthreads,
        FLUX_NI * FLUX_NJ * FLUX_NK * sizeof(double), strm>>>(
        s->d_q, s->d_uloc, s->d_flux, dx, dy, dz, dt, stride);
    cudaEventRecord(ev[6], strm);

    // D -> H: download results
    cudaMemcpyAsync(h_flux, s->d_flux, flux_bytes, cudaMemcpyDeviceToHost, strm);
    cudaMemcpyAsync(h_tmp,  s->d_tmp,  tmp_bytes,  cudaMemcpyDeviceToHost, strm);
    cudaEventRecord(ev[7], strm);
}

void hydro_cuda_gather_unsplit_sync(
    double* h_flux,
    double* h_tmp,
    int ngrid,
    int stream_slot)
{
    StreamSlot* s = &get_pool()[stream_slot];
    cudaStreamSynchronize(s->stream);
    cudaError_t err = cudaGetLastError();
    if (err != cudaSuccess) {
        fprintf(stderr, "CUDA gather-hydro error (slot %d, ngrid %d): %s\n",
                stream_slot, ngrid, cudaGetErrorString(err));
    }
}

// ====================================================================
// Scatter-reduce: 5 standard kernels + scatter_reduce → compact D2H
// Eliminates 98 MB D2H transfer, replaces with ~5 MB.
// ====================================================================
void hydro_cuda_unsplit_reduce_async(
    const double* h_uloc, const double* h_gloc, const int* h_ok,
    double* h_add_unew, double* h_add_lm1,
    double* h_add_divu_l, double* h_add_enew_l,
    double* h_add_divu_lm1, double* h_add_enew_lm1,
    double dx, double dy, double dz, double dt,
    int ngrid, int stride, int stream_slot)
{
    pool_ensure_hydro_buffers(stream_slot, stride);
    pool_ensure_hydro_inter_buffers(stream_slot, stride);
    pool_ensure_reduce_buffers(stream_slot, stride);

    StreamSlot* s = &get_pool()[stream_slot];
    cudaStream_t strm = s->stream;

    size_t uloc_bytes = (size_t)stride * STENCIL_NI * STENCIL_NJ * STENCIL_NK * NVAR * sizeof(double);
    size_t gloc_bytes = (size_t)stride * STENCIL_NI * STENCIL_NJ * STENCIL_NK * NDIM * sizeof(double);
    size_t ok_bytes   = (size_t)stride * STENCIL_NI * STENCIL_NJ * STENCIL_NK * sizeof(int);
    size_t flux_bytes = (size_t)stride * FLUX_NI * FLUX_NJ * FLUX_NK * NVAR * NDIM * sizeof(double);
    size_t tmp_bytes  = (size_t)stride * FLUX_NI * FLUX_NJ * FLUX_NK * 2 * NDIM * sizeof(double);

    cudaEvent_t* ev = s->ev_profile;

    // H2D: direct DMA from host (caller must ensure pinned/registered memory)
    cudaEventRecord(ev[0], strm);
    cudaMemcpyAsync(s->d_uloc,   h_uloc, uloc_bytes, cudaMemcpyHostToDevice, strm);
    cudaMemcpyAsync(s->d_gloc,   h_gloc, gloc_bytes, cudaMemcpyHostToDevice, strm);
    cudaMemcpyAsync(s->d_ok_int, h_ok,   ok_bytes,   cudaMemcpyHostToDevice, strm);
    cudaEventRecord(ev[1], strm);

    int nblocks  = ngrid;
    int nthreads = NCELLS_PER_GRID;  // 216
    double dtdx = dt/dx, dtdy = dt/dy, dtdz = dt/dz;

    // Kernel 1: ctoprim
    hydro_ctoprim_kernel<<<nblocks, nthreads, 0, strm>>>(
        s->d_uloc, s->d_gloc, s->d_q, s->d_c, dt, stride);
    cudaEventRecord(ev[2], strm);

    // Kernel 2: uslope
    hydro_uslope_kernel<<<nblocks, nthreads, 0, strm>>>(
        s->d_q, s->d_dq, stride);
    cudaEventRecord(ev[3], strm);

    // Kernel 3: trace3d
    hydro_trace3d_kernel<<<nblocks, nthreads, 0, strm>>>(
        s->d_q, s->d_dq, s->d_qm, s->d_qp, dtdx, dtdy, dtdz, stride);
    cudaEventRecord(ev[4], strm);

    // Zero flux/tmp
    cudaMemsetAsync(s->d_flux, 0, flux_bytes, strm);
    cudaMemsetAsync(s->d_tmp,  0, tmp_bytes,  strm);

    // Kernel 4: Riemann fluxes
    hydro_flux_kernel<<<nblocks, nthreads, 0, strm>>>(
        s->d_qm, s->d_qp, s->d_flux, s->d_tmp, dtdx, dtdy, dtdz, stride);
    cudaEventRecord(ev[5], strm);

    // Kernel 5: difmag
    hydro_difmag_kernel<<<nblocks, nthreads,
        FLUX_NI * FLUX_NJ * FLUX_NK * sizeof(double), strm>>>(
        s->d_q, s->d_uloc, s->d_flux, dx, dy, dz, dt, stride);
    cudaEventRecord(ev[6], strm);

    // Kernel 6: scatter_reduce (flux+tmp+ok → add_unew+add_lm1)
    double oneontwotondim = 0.125;  // 1/8 for 3D
    hydro_scatter_reduce_kernel<<<nblocks, 256, 0, strm>>>(
        s->d_flux, s->d_tmp, s->d_ok_int,
        s->d_add_unew, s->d_add_lm1,
        s->d_add_divu_l, s->d_add_enew_l,
        s->d_add_divu_lm1, s->d_add_enew_lm1,
        oneontwotondim, ngrid, stride);

    // D2H: compact output only (~5 MB vs 98 MB)
    size_t unew_bytes     = (size_t)stride * 8 * NVAR * sizeof(double);
    size_t lm1_bytes      = (size_t)stride * 6 * NVAR * sizeof(double);
    size_t divu_l_bytes   = (size_t)stride * 8 * sizeof(double);
    size_t enew_l_bytes   = (size_t)stride * 8 * sizeof(double);
    size_t divu_lm1_bytes = (size_t)stride * 6 * sizeof(double);
    size_t enew_lm1_bytes = (size_t)stride * 6 * sizeof(double);

    cudaMemcpyAsync(h_add_unew, s->d_add_unew, unew_bytes, cudaMemcpyDeviceToHost, strm);
    cudaMemcpyAsync(h_add_lm1,  s->d_add_lm1,  lm1_bytes,  cudaMemcpyDeviceToHost, strm);

    if (h_add_divu_l) {
        cudaMemcpyAsync(h_add_divu_l,   s->d_add_divu_l,   divu_l_bytes,   cudaMemcpyDeviceToHost, strm);
        cudaMemcpyAsync(h_add_enew_l,   s->d_add_enew_l,   enew_l_bytes,   cudaMemcpyDeviceToHost, strm);
        cudaMemcpyAsync(h_add_divu_lm1, s->d_add_divu_lm1, divu_lm1_bytes, cudaMemcpyDeviceToHost, strm);
        cudaMemcpyAsync(h_add_enew_lm1, s->d_add_enew_lm1, enew_lm1_bytes, cudaMemcpyDeviceToHost, strm);
    }

    cudaEventRecord(ev[7], strm);
}

void hydro_cuda_unsplit_reduce_sync(int ngrid, int stream_slot)
{
    StreamSlot* s = &get_pool()[stream_slot];
    cudaStreamSynchronize(s->stream);
    cudaError_t err = cudaGetLastError();
    if (err != cudaSuccess) {
        fprintf(stderr, "CUDA hydro reduce error (slot %d, ngrid %d): %s\n",
                stream_slot, ngrid, cudaGetErrorString(err));
    }
}

// ====================================================================
// GPU-gather + scatter-reduce pipeline:
// H2D stencil idx -> gather kernel -> compute kernels -> scatter_reduce
// -> D2H compact output.
// Eliminates 98 MB H2D of uloc/gloc, replaces with ~7 MB stencil indices.
// Requires cuda_mesh_upload() to have been called first.
// ====================================================================
void hydro_cuda_gather_reduce_async(
    const int* h_stencil_idx, const int* h_stencil_grav,
    const double* h_interp_vals,
    double* h_add_unew, double* h_add_lm1,
    double* h_add_divu_l, double* h_add_enew_l,
    double* h_add_divu_lm1, double* h_add_enew_lm1,
    double dx, double dy, double dz, double dt,
    int ngrid, int stride, int n_interp, int stream_slot)
{
    pool_ensure_hydro_buffers(stream_slot, stride);
    pool_ensure_hydro_inter_buffers(stream_slot, stride);
    pool_ensure_stencil_buffers(stream_slot, stride, n_interp);
    pool_ensure_reduce_buffers(stream_slot, stride);

    StreamSlot* s = &get_pool()[stream_slot];
    cudaStream_t strm = s->stream;
    long long ncell = cuda_get_mesh_ncell();

    cudaEvent_t* ev = s->ev_profile;

    // H2D: stencil indices (small: ~7 MB for stride=4096)
    size_t idx_bytes = (size_t)stride * STENCIL_NI * STENCIL_NJ * STENCIL_NK * sizeof(int);
    cudaEventRecord(ev[0], strm);
    cudaMemcpyAsync(s->d_stencil_idx,  h_stencil_idx,  idx_bytes, cudaMemcpyHostToDevice, strm);
    cudaMemcpyAsync(s->d_stencil_grav, h_stencil_grav, idx_bytes, cudaMemcpyHostToDevice, strm);

    // Upload interpolation values (AoS: n_interp * NVAR doubles)
    if (n_interp > 0 && h_interp_vals) {
        size_t interp_bytes = (size_t)n_interp * NVAR * sizeof(double);
        cudaMemcpyAsync(s->d_interp_vals, h_interp_vals, interp_bytes, cudaMemcpyHostToDevice, strm);
    }
    cudaEventRecord(ev[1], strm);

    int nblocks  = ngrid;
    int nthreads = NCELLS_PER_GRID;  // 216

    // Kernel 0a: GPU gather (fills d_uloc and d_gloc from mesh)
    hydro_gather_kernel<<<nblocks, nthreads, 0, strm>>>(
        s->d_stencil_idx, s->d_stencil_grav, s->d_interp_vals,
        cuda_get_mesh_uold(), cuda_get_mesh_f(),
        s->d_uloc, s->d_gloc,
        ncell, ngrid, stride);

    // Kernel 0b: Compute ok flags from mesh_son via stencil_idx
    hydro_gather_ok_kernel<<<nblocks, nthreads, 0, strm>>>(
        s->d_stencil_idx, cuda_get_mesh_son(),
        s->d_ok_int,
        ngrid, stride);

    double dtdx = dt / dx;
    double dtdy = dt / dy;
    double dtdz = dt / dz;

    // Kernel 1: ctoprim
    hydro_ctoprim_kernel<<<nblocks, nthreads, 0, strm>>>(
        s->d_uloc, s->d_gloc, s->d_q, s->d_c, dt, stride);
    cudaEventRecord(ev[2], strm);

    // Kernel 2: uslope
    hydro_uslope_kernel<<<nblocks, nthreads, 0, strm>>>(
        s->d_q, s->d_dq, stride);
    cudaEventRecord(ev[3], strm);

    // Kernel 3: trace3d
    hydro_trace3d_kernel<<<nblocks, nthreads, 0, strm>>>(
        s->d_q, s->d_dq, s->d_qm, s->d_qp, dtdx, dtdy, dtdz, stride);
    cudaEventRecord(ev[4], strm);

    // Zero flux and tmp output arrays
    size_t flux_bytes = (size_t)stride * FLUX_NI * FLUX_NJ * FLUX_NK * NVAR * NDIM * sizeof(double);
    size_t tmp_bytes  = (size_t)stride * FLUX_NI * FLUX_NJ * FLUX_NK * 2 * NDIM * sizeof(double);
    cudaMemsetAsync(s->d_flux, 0, flux_bytes, strm);
    cudaMemsetAsync(s->d_tmp,  0, tmp_bytes,  strm);

    // Kernel 4: Riemann fluxes
    hydro_flux_kernel<<<nblocks, nthreads, 0, strm>>>(
        s->d_qm, s->d_qp, s->d_flux, s->d_tmp, dtdx, dtdy, dtdz, stride);
    cudaEventRecord(ev[5], strm);

    // Kernel 5: difmag
    hydro_difmag_kernel<<<nblocks, nthreads,
        FLUX_NI * FLUX_NJ * FLUX_NK * sizeof(double), strm>>>(
        s->d_q, s->d_uloc, s->d_flux, dx, dy, dz, dt, stride);
    cudaEventRecord(ev[6], strm);

    // Kernel 6: scatter_reduce (flux+tmp+ok -> add_unew+add_lm1)
    double oneontwotondim = 0.125;  // 1/8 for 3D
    hydro_scatter_reduce_kernel<<<nblocks, 256, 0, strm>>>(
        s->d_flux, s->d_tmp, s->d_ok_int,
        s->d_add_unew, s->d_add_lm1,
        s->d_add_divu_l, s->d_add_enew_l,
        s->d_add_divu_lm1, s->d_add_enew_lm1,
        oneontwotondim, ngrid, stride);

    // D2H: compact output only (~5 MB vs 98 MB)
    size_t unew_bytes = (size_t)stride * 8 * NVAR * sizeof(double);
    size_t lm1_bytes  = (size_t)stride * 6 * NVAR * sizeof(double);

    cudaMemcpyAsync(h_add_unew, s->d_add_unew, unew_bytes, cudaMemcpyDeviceToHost, strm);
    cudaMemcpyAsync(h_add_lm1,  s->d_add_lm1,  lm1_bytes,  cudaMemcpyDeviceToHost, strm);

    if (h_add_divu_l) {
        size_t divu_l_bytes   = (size_t)stride * 8 * sizeof(double);
        size_t enew_l_bytes   = (size_t)stride * 8 * sizeof(double);
        size_t divu_lm1_bytes = (size_t)stride * 6 * sizeof(double);
        size_t enew_lm1_bytes = (size_t)stride * 6 * sizeof(double);
        cudaMemcpyAsync(h_add_divu_l,   s->d_add_divu_l,   divu_l_bytes,   cudaMemcpyDeviceToHost, strm);
        cudaMemcpyAsync(h_add_enew_l,   s->d_add_enew_l,   enew_l_bytes,   cudaMemcpyDeviceToHost, strm);
        cudaMemcpyAsync(h_add_divu_lm1, s->d_add_divu_lm1, divu_lm1_bytes, cudaMemcpyDeviceToHost, strm);
        cudaMemcpyAsync(h_add_enew_lm1, s->d_add_enew_lm1, enew_lm1_bytes, cudaMemcpyDeviceToHost, strm);
    }

    cudaEventRecord(ev[7], strm);
}

void hydro_cuda_gather_reduce_sync(int ngrid, int stream_slot)
{
    StreamSlot* s = &get_pool()[stream_slot];
    cudaStreamSynchronize(s->stream);
    cudaError_t err = cudaGetLastError();
    if (err != cudaSuccess) {
        fprintf(stderr, "CUDA gather-reduce error (slot %d, ngrid %d): %s\n",
                stream_slot, ngrid, cudaGetErrorString(err));
    }
}

} // extern "C"

// ====================================================================
// Auxiliary CUDA Kernels: gradient_phi and upload_fine
// For dynamic OMP/CUDA hybrid dispatch
// ====================================================================

// gradient_phi: 5-point FDA force from pre-gathered phi stencil
// phi_buf layout: phi_buf[idx*4 + 0..3] = phi1,phi2,phi3,phi4
// f_buf layout: f_buf[idx] = a*(phi1-phi2) - b*(phi3-phi4)
// total = ngrid * 8 * 3 (cells * dims)
__global__ void gradient_phi_kernel(
    const double* __restrict__ phi_buf,
    double* __restrict__ f_buf,
    double a, double b,
    int total)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= total) return;
    int base = idx * 4;
    f_buf[idx] = a * (phi_buf[base] - phi_buf[base+1])
               - b * (phi_buf[base+2] - phi_buf[base+3]);
}

// upload_fine: restriction (average 8 children per parent)
// child_buf layout (Fortran col-major): child_buf[v + nvar*(c + 8*p)]
// parent_buf layout: parent_buf[v + (nvar+1)*p]
//   v=0..nvar-1: averaged conservative variables
//   v=nvar: averaged internal energy (if do_eint)
__global__ void upload_fine_kernel(
    const double* __restrict__ child_buf,
    double* __restrict__ parent_buf,
    int nsplit, int nvar, int ndim, double smallr, int do_eint)
{
    int p = blockIdx.x * blockDim.x + threadIdx.x;
    if (p >= nsplit) return;

    int out_stride = nvar + 1;

    // Average all conservative variables
    for (int v = 0; v < nvar; v++) {
        double sum = 0;
        for (int c = 0; c < 8; c++) {
            sum += child_buf[v + nvar * (c + 8 * p)];
        }
        parent_buf[v + out_stride * p] = sum * 0.125;
    }

    // Compute averaged internal energy if needed
    if (do_eint) {
        double avg_eint = 0;
        for (int c = 0; c < 8; c++) {
            int base = nvar * (c + 8 * p);
            double rho = fmax(child_buf[base], smallr);
            double ekin = 0;
            for (int d = 0; d < ndim; d++) {
                double mom = child_buf[base + 1 + d];
                ekin += 0.5 * mom * mom / rho;
            }
            double etot = child_buf[base + ndim + 1];
            // erad = 0 for NENER=0
            avg_eint += (etot - ekin);
        }
        parent_buf[nvar + out_stride * p] = avg_eint * 0.125;
    }
}

// Per-stream auxiliary buffers for gradient_phi and upload_fine
struct AuxBuffers {
    double *d_phi_buf, *d_f_buf;
    int phi_cap;
    double *d_child_buf, *d_parent_buf;
    int child_cap, child_nvar;
};

static AuxBuffers g_aux[N_STREAMS];
static bool g_aux_initialized = false;

static void init_aux_buffers() {
    if (g_aux_initialized) return;
    memset(g_aux, 0, sizeof(g_aux));
    g_aux_initialized = true;
}

extern "C" {

void gradient_phi_cuda_async(
    const double* h_phi_buf, double* h_f_buf,
    double a, double b, int ngrid, int stream_slot)
{
    if (!is_pool_initialized()) return;
    init_aux_buffers();
    AuxBuffers& aux = g_aux[stream_slot];
    cudaStream_t strm = cuda_get_stream_internal(stream_slot);

    if (ngrid > aux.phi_cap) {
        int cap = ngrid * 2;
        if (aux.d_phi_buf) cudaFree(aux.d_phi_buf);
        if (aux.d_f_buf)   cudaFree(aux.d_f_buf);
        cudaMalloc(&aux.d_phi_buf, (size_t)cap * 24 * 4 * sizeof(double));
        cudaMalloc(&aux.d_f_buf,   (size_t)cap * 24 * sizeof(double));
        aux.phi_cap = cap;
    }

    int total = ngrid * 24;
    size_t phi_bytes = (size_t)total * 4 * sizeof(double);
    size_t f_bytes   = (size_t)total * sizeof(double);

    cudaMemcpyAsync(aux.d_phi_buf, h_phi_buf, phi_bytes, cudaMemcpyHostToDevice, strm);
    int nth = 256;
    gradient_phi_kernel<<<(total+nth-1)/nth, nth, 0, strm>>>(
        aux.d_phi_buf, aux.d_f_buf, a, b, total);
    cudaMemcpyAsync(h_f_buf, aux.d_f_buf, f_bytes, cudaMemcpyDeviceToHost, strm);
}

void gradient_phi_cuda_sync(int stream_slot) {
    if (!is_pool_initialized()) return;
    cudaStreamSynchronize(cuda_get_stream_internal(stream_slot));
    cudaError_t err = cudaGetLastError();
    if (err != cudaSuccess)
        fprintf(stderr, "CUDA gradient_phi error (slot %d): %s\n",
                stream_slot, cudaGetErrorString(err));
}

void upload_fine_cuda_async(
    const double* h_child_buf, double* h_parent_buf,
    int nsplit, int nvar, int ndim,
    double smallr, int do_eint, int stream_slot)
{
    if (!is_pool_initialized()) return;
    init_aux_buffers();
    AuxBuffers& aux = g_aux[stream_slot];
    cudaStream_t strm = cuda_get_stream_internal(stream_slot);

    if (nsplit > aux.child_cap || nvar != aux.child_nvar) {
        int cap = (nsplit > 512) ? nsplit * 2 : 1024;
        if (aux.d_child_buf)  cudaFree(aux.d_child_buf);
        if (aux.d_parent_buf) cudaFree(aux.d_parent_buf);
        cudaMalloc(&aux.d_child_buf,  (size_t)cap * 8 * nvar * sizeof(double));
        cudaMalloc(&aux.d_parent_buf, (size_t)cap * (nvar + 1) * sizeof(double));
        aux.child_cap = cap;
        aux.child_nvar = nvar;
    }

    size_t child_bytes  = (size_t)nsplit * 8 * nvar * sizeof(double);
    size_t parent_bytes = (size_t)nsplit * (nvar + 1) * sizeof(double);

    cudaMemcpyAsync(aux.d_child_buf, h_child_buf, child_bytes, cudaMemcpyHostToDevice, strm);
    int nth = 256;
    upload_fine_kernel<<<(nsplit+nth-1)/nth, nth, 0, strm>>>(
        aux.d_child_buf, aux.d_parent_buf, nsplit, nvar, ndim, smallr, do_eint);
    cudaMemcpyAsync(h_parent_buf, aux.d_parent_buf, parent_bytes, cudaMemcpyDeviceToHost, strm);
}

void upload_fine_cuda_sync(int stream_slot) {
    if (!is_pool_initialized()) return;
    cudaStreamSynchronize(cuda_get_stream_internal(stream_slot));
    cudaError_t err = cudaGetLastError();
    if (err != cudaSuccess)
        fprintf(stderr, "CUDA upload_fine error (slot %d): %s\n",
                stream_slot, cudaGetErrorString(err));
}

// ====================================================================
// synchro_hydro: gravity momentum + energy update
// buf layout (Fortran column-major, SoA):
//   buf[k + 0*stride] = rho, buf[k + d*stride] = momentum_d (d=1..ndim),
//   buf[k + (ndim+1)*stride] = E_total,
//   buf[k + (ndim+2+d-1)*stride] = f_grav_d (d=0..ndim-1)
// In-place update: momentum += rho*f*dt, energy adjusted for kinetic
// ====================================================================
__global__ void synchro_hydro_kernel(
    double* __restrict__ buf,
    double smallr, double dteff,
    int ncell, int ndim, int stride)
{
    int k = blockIdx.x * blockDim.x + threadIdx.x;
    if (k >= ncell) return;

    double rho = buf[k];
    if (rho < smallr) rho = smallr;
    double E = buf[k + (ndim + 1) * stride];

    // Remove kinetic energy
    for (int d = 0; d < ndim; d++) {
        double mom = buf[k + (d + 1) * stride];
        E -= 0.5 * mom * mom / rho;
    }

    // Update momentum
    for (int d = 0; d < ndim; d++) {
        double fg = buf[k + (ndim + 2 + d) * stride];
        buf[k + (d + 1) * stride] += rho * fg * dteff;
    }

    // Add back kinetic energy with updated momentum
    for (int d = 0; d < ndim; d++) {
        double mom = buf[k + (d + 1) * stride];
        E += 0.5 * mom * mom / rho;
    }

    buf[k + (ndim + 1) * stride] = E;
}

// Per-stream generic device buffers for synchro/cmpdt
static double* g_generic_buf[N_STREAMS] = {};
static size_t  g_generic_cap[N_STREAMS] = {};

static void ensure_generic_buf(int slot, size_t needed) {
    if (needed > g_generic_cap[slot]) {
        if (g_generic_buf[slot]) cudaFree(g_generic_buf[slot]);
        cudaError_t err = cudaMalloc(&g_generic_buf[slot], needed * 2);
        if (err != cudaSuccess) {
            fprintf(stderr, "ensure_generic_buf: cudaMalloc(%zu) failed for slot %d: %s\n",
                    needed * 2, slot, cudaGetErrorString(err));
            g_generic_buf[slot] = nullptr;
            g_generic_cap[slot] = 0;
            return;
        }
        g_generic_cap[slot] = needed * 2;
    }
}

void synchro_cuda_async(
    double* h_buf,
    double smallr, double dteff,
    int ncell, int ndim, int stride, int stream_slot)
{
    if (!is_pool_initialized()) return;
    cudaStream_t strm = cuda_get_stream_internal(stream_slot);

    int nprops = ndim + 2 + ndim;
    size_t total_bytes = (size_t)stride * nprops * sizeof(double);

    ensure_generic_buf(stream_slot, total_bytes);
    if (!g_generic_buf[stream_slot]) return;

    cudaError_t err;
    err = cudaMemcpyAsync(g_generic_buf[stream_slot], h_buf, total_bytes, cudaMemcpyHostToDevice, strm);
    if (err != cudaSuccess) {
        fprintf(stderr, "synchro H2D memcpy error: %s (slot=%d bytes=%zu h_buf=%p d_buf=%p)\n",
                cudaGetErrorString(err), stream_slot, total_bytes, (void*)h_buf, (void*)g_generic_buf[stream_slot]);
        return;
    }
    int nth = 256;
    synchro_hydro_kernel<<<(ncell+nth-1)/nth, nth, 0, strm>>>(
        g_generic_buf[stream_slot], smallr, dteff, ncell, ndim, stride);
    err = cudaMemcpyAsync(h_buf, g_generic_buf[stream_slot], total_bytes, cudaMemcpyDeviceToHost, strm);
    if (err != cudaSuccess) {
        fprintf(stderr, "synchro D2H memcpy error: %s (slot=%d bytes=%zu)\n",
                cudaGetErrorString(err), stream_slot, total_bytes);
    }
}

void synchro_cuda_sync(int stream_slot) {
    if (!is_pool_initialized()) return;
    cudaStreamSynchronize(cuda_get_stream_internal(stream_slot));
    cudaError_t err = cudaGetLastError();
    if (err != cudaSuccess)
        fprintf(stderr, "CUDA synchro error (slot %d): %s\n",
                stream_slot, cudaGetErrorString(err));
}

// ====================================================================
// cmpdt: CFL timestep computation per cell
// Input buf (SoA): uold(1:nvar) + f(1:ndim) per cell
// Output dt_buf: per-cell dt value (min-reduced on CPU)
// ====================================================================
__global__ void cmpdt_cuda_kernel(
    const double* __restrict__ buf,
    double* __restrict__ dt_buf,
    double dx, double smallr, double smallc, double gamma_val, double courant_factor,
    int ncell, int nvar, int ndim, int stride)
{
    int k = blockIdx.x * blockDim.x + threadIdx.x;
    if (k >= ncell) return;

    double rho = buf[k];
    if (rho < smallr) rho = smallr;

    // Internal energy = E_total
    double eint = buf[k + (ndim + 1) * stride];
    // Subtract kinetic energy
    for (int d = 0; d < ndim; d++) {
        double v = buf[k + (d + 1) * stride] / rho;
        eint -= 0.5 * rho * v * v;
    }

    // Pressure
    double smallp = smallc * smallc / gamma_val;
    double p = (gamma_val - 1.0) * eint;
    if (p < rho * smallp) p = rho * smallp;

    // Sound speed
    double cs = sqrt(gamma_val * p / rho);

    // Wave speed = ndim*cs + sum(|v|)
    double wavespeed = (double)ndim * cs;
    for (int d = 0; d < ndim; d++) {
        double v = buf[k + (d + 1) * stride] / rho;
        wavespeed += fabs(v);
    }

    // Gravity strength ratio
    double grav_ratio = 0.0;
    for (int d = 0; d < ndim; d++) {
        grav_ratio += fabs(buf[k + (nvar + d) * stride]);
    }
    grav_ratio = grav_ratio * dx / (wavespeed * wavespeed);
    if (grav_ratio < 0.0001) grav_ratio = 0.0001;

    // CFL timestep
    double dtcell = dx / wavespeed * (sqrt(1.0 + 2.0 * courant_factor * grav_ratio) - 1.0) / grav_ratio;
    dt_buf[k] = dtcell;
}

void cmpdt_cuda_async(
    double* h_buf, double* h_dt_buf,
    double dx, double smallr, double smallc, double gamma_val, double courant_factor,
    int ncell, int nvar, int ndim, int stride, int stream_slot)
{
    if (!is_pool_initialized()) return;
    cudaStream_t strm = cuda_get_stream_internal(stream_slot);

    int nprops = nvar + ndim;
    size_t buf_bytes = (size_t)stride * nprops * sizeof(double);
    size_t dt_bytes  = (size_t)ncell * sizeof(double);

    ensure_generic_buf(stream_slot, buf_bytes + dt_bytes);
    if (!g_generic_buf[stream_slot]) return;
    double* d_buf = g_generic_buf[stream_slot];
    double* d_dt  = d_buf + (size_t)stride * nprops;

    cudaMemcpyAsync(d_buf, h_buf, buf_bytes, cudaMemcpyHostToDevice, strm);
    int nth = 256;
    cmpdt_cuda_kernel<<<(ncell+nth-1)/nth, nth, 0, strm>>>(
        d_buf, d_dt, dx, smallr, smallc, gamma_val, courant_factor,
        ncell, nvar, ndim, stride);
    cudaMemcpyAsync(h_dt_buf, d_dt, dt_bytes, cudaMemcpyDeviceToHost, strm);
}

void cmpdt_cuda_sync(int stream_slot) {
    if (!is_pool_initialized()) return;
    cudaStreamSynchronize(cuda_get_stream_internal(stream_slot));
    cudaError_t err = cudaGetLastError();
    if (err != cudaSuccess)
        fprintf(stderr, "CUDA cmpdt error (slot %d): %s\n",
                stream_slot, cudaGetErrorString(err));
}

} // extern "C" -- auxiliary kernels
