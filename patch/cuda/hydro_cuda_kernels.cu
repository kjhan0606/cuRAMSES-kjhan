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

    // H -> D: upload input arrays (full NVECTOR-strided arrays)
    cudaMemcpyAsync(s->d_uloc, h_uloc, uloc_bytes, cudaMemcpyHostToDevice, strm);
    cudaMemcpyAsync(s->d_gloc, h_gloc, gloc_bytes, cudaMemcpyHostToDevice, strm);

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

    // Kernel 2: TVD slopes
    hydro_uslope_kernel<<<nblocks, nthreads, 0, strm>>>(
        s->d_q, s->d_dq, stride);

    // Kernel 3: MUSCL-Hancock trace
    hydro_trace3d_kernel<<<nblocks, nthreads, 0, strm>>>(
        s->d_q, s->d_dq, s->d_qm, s->d_qp, dtdx, dtdy, dtdz, stride);

    // Zero flux and tmp output arrays
    cudaMemsetAsync(s->d_flux, 0, flux_bytes, strm);
    cudaMemsetAsync(s->d_tmp,  0, tmp_bytes,  strm);

    // Kernel 4: Riemann solves -> fluxes
    hydro_flux_kernel<<<nblocks, nthreads, 0, strm>>>(
        s->d_qm, s->d_qp, s->d_flux, s->d_tmp, dtdx, dtdy, dtdz, stride);

    // Kernel 5: Artificial viscosity (difmag kernel)
    hydro_difmag_kernel<<<nblocks, nthreads,
        FLUX_NI * FLUX_NJ * FLUX_NK * sizeof(double), strm>>>(
        s->d_q, s->d_uloc, s->d_flux, dx, dy, dz, dt, stride);

    // D -> H: download results (full NVECTOR-strided arrays)
    cudaMemcpyAsync(h_flux, s->d_flux, flux_bytes, cudaMemcpyDeviceToHost, strm);
    cudaMemcpyAsync(h_tmp,  s->d_tmp,  tmp_bytes,  cudaMemcpyDeviceToHost, strm);
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

    // Upload stencil indices (small: ~7 MB for stride=4096)
    size_t idx_bytes = (size_t)stride * STENCIL_NI * STENCIL_NJ * STENCIL_NK * sizeof(int);
    cudaMemcpyAsync(s->d_stencil_idx,  h_stencil_idx,  idx_bytes, cudaMemcpyHostToDevice, strm);
    cudaMemcpyAsync(s->d_stencil_grav, h_stencil_grav, idx_bytes, cudaMemcpyHostToDevice, strm);

    // Upload interpolation values (AoS: n_interp * NVAR doubles)
    if (n_interp > 0 && h_interp_vals) {
        size_t interp_bytes = (size_t)n_interp * NVAR * sizeof(double);
        cudaMemcpyAsync(s->d_interp_vals, h_interp_vals, interp_bytes, cudaMemcpyHostToDevice, strm);
    }

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

    // Kernel 2: TVD slopes
    hydro_uslope_kernel<<<nblocks, nthreads, 0, strm>>>(
        s->d_q, s->d_dq, stride);

    // Kernel 3: MUSCL-Hancock trace
    hydro_trace3d_kernel<<<nblocks, nthreads, 0, strm>>>(
        s->d_q, s->d_dq, s->d_qm, s->d_qp, dtdx, dtdy, dtdz, stride);

    // Zero flux and tmp output arrays
    size_t flux_bytes = (size_t)stride * FLUX_NI * FLUX_NJ * FLUX_NK * NVAR * NDIM * sizeof(double);
    size_t tmp_bytes  = (size_t)stride * FLUX_NI * FLUX_NJ * FLUX_NK * 2 * NDIM * sizeof(double);
    cudaMemsetAsync(s->d_flux, 0, flux_bytes, strm);
    cudaMemsetAsync(s->d_tmp,  0, tmp_bytes,  strm);

    // Kernel 4: Riemann solves -> fluxes
    hydro_flux_kernel<<<nblocks, nthreads, 0, strm>>>(
        s->d_qm, s->d_qp, s->d_flux, s->d_tmp, dtdx, dtdy, dtdz, stride);

    // Kernel 5: Artificial viscosity
    hydro_difmag_kernel<<<nblocks, nthreads,
        FLUX_NI * FLUX_NJ * FLUX_NK * sizeof(double), strm>>>(
        s->d_q, s->d_uloc, s->d_flux, dx, dy, dz, dt, stride);

    // D -> H: download results
    cudaMemcpyAsync(h_flux, s->d_flux, flux_bytes, cudaMemcpyDeviceToHost, strm);
    cudaMemcpyAsync(h_tmp,  s->d_tmp,  tmp_bytes,  cudaMemcpyDeviceToHost, strm);
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

} // extern "C"
