// ==========================================================================
// GPU-accelerated AGN feedback kernels for cuRAMSES
// average_AGN + AGN_blast with flat cell arrays and spatial binning
// ==========================================================================

#include "cuda_stream_pool.h"
#include <cstdio>
#include <cmath>

// ==========================================================================
// Static GPU arrays (persistent across calls, resized as needed)
// ==========================================================================

// Flat cell arrays
static double *d_cell_x = nullptr, *d_cell_y = nullptr, *d_cell_z = nullptr;
static double *d_cell_vol = nullptr, *d_cell_d = nullptr, *d_cell_dx = nullptr;
static int    *d_cell_idx = nullptr;
static double *d_cell_uold = nullptr;  // [ncells * nvar] for AGN_blast

// AGN data arrays
static double *d_agn_buf = nullptr;     // packed AGN doubles
static int    *d_agn_ibuf = nullptr;    // packed AGN ints

// Spatial hash
static int    *d_bin_head = nullptr;
static int    *d_agn_next = nullptr;

// Per-AGN output arrays
static double *d_vol_gas = nullptr, *d_mass_gas = nullptr, *d_psy_norm = nullptr;
static double *d_vol_blast = nullptr, *d_mass_blast = nullptr;
static int    *d_ind_blast = nullptr;
static double *d_Esave_out = nullptr;

// Allocation tracking
static int g_cell_cap = 0;
static int g_agn_cap = 0;
static int g_bin_cap = 0;
static int g_uold_cap = 0;

// Dedicated stream
static cudaStream_t g_sink_stream = nullptr;

// ==========================================================================
// average_AGN GPU kernel
// Thread per leaf cell: 27-bin search, atomicAdd to per-AGN accumulators
// ==========================================================================
__global__ void kernel_average_AGN(
    const double* __restrict__ cell_x,
    const double* __restrict__ cell_y,
    const double* __restrict__ cell_z,
    const double* __restrict__ cell_vol,
    const double* __restrict__ cell_d,
    const double* __restrict__ cell_dx,
    const int*    __restrict__ cell_idx,
    int ncells,
    const double* __restrict__ agn_x,       // SoA: [x|y|z] x nAGN
    const double* __restrict__ agn_j,       // SoA: [jx|jy|jz] x nAGN
    const double* __restrict__ agn_dMBH,
    const double* __restrict__ agn_dMEd,
    const double* __restrict__ agn_Esave,
    const int*    __restrict__ agn_ok_blast,
    int nAGN,
    const int*    __restrict__ bin_head,     // flat [nbin^3], Fortran col-major
    const int*    __restrict__ agn_next,     // [nAGN], 1-based linked list
    int nbin, double inv_bin_size, double rmax, double rmax2, double X_floor,
    double* vol_gas, double* mass_gas, double* psy_norm,
    double* vol_blast, double* mass_blast, int* ind_blast)
{
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if(tid >= ncells) return;

    double x = cell_x[tid];
    double y = cell_y[tid];
    double z = cell_z[tid];
    double vloc = cell_vol[tid];
    double dval = cell_d[tid];
    double dxloc = cell_dx[tid];
    int    idx  = cell_idx[tid];

    int ibx = max(1, min(nbin, (int)(x * inv_bin_size) + 1));
    int iby = max(1, min(nbin, (int)(y * inv_bin_size) + 1));
    int ibz = max(1, min(nbin, (int)(z * inv_bin_size) + 1));

    for(int jbz = max(1, ibz-1); jbz <= min(nbin, ibz+1); jbz++)
    for(int jby = max(1, iby-1); jby <= min(nbin, iby+1); jby++)
    for(int jbx = max(1, ibx-1); jbx <= min(nbin, ibx+1); jbx++) {
        int bidx = (jbx-1) + (jby-1)*nbin + (jbz-1)*nbin*nbin;
        int iA = bin_head[bidx];  // 1-based
        while(iA > 0) {
            int ai = iA - 1;  // 0-based
            double dxx = x - agn_x[ai];
            double dyy = y - agn_x[nAGN + ai];
            double dzz = z - agn_x[2*nAGN + ai];
            double dr2 = dxx*dxx + dyy*dyy + dzz*dzz;
            double Xr  = agn_dMBH[ai] / agn_dMEd[ai];

            if(agn_Esave[ai] > 0.0) {
                // Case 0: saved energy from previous step
                if(dr2 <= rmax2) {
                    atomicAdd(&vol_gas[ai], vloc * dval);
                }
                double drc = fmax(fabs(dxx), fmax(fabs(dyy), fabs(dzz)));
                if(drc <= dxloc * 0.5) {
                    ind_blast[ai] = idx;
                    vol_blast[ai] = vloc;
                }
            } else {
                if(Xr < X_floor) {
                    // Jet case
                    if(agn_ok_blast[ai]) {
                        double jx = agn_j[ai];
                        double jy = agn_j[nAGN + ai];
                        double jz = agn_j[2*nAGN + ai];
                        double jtot = sqrt(jx*jx + jy*jy + jz*jz);
                        if(jtot > 0.0) {
                            double nj_x = jx/jtot, nj_y = jy/jtot, nj_z = jz/jtot;
                            double dzjet = dxx*nj_x + dyy*nj_y + dzz*nj_z;
                            double drj2 = dr2 - dzjet*dzjet;
                            double drjet = (drj2 > 0.0) ? sqrt(drj2) : 0.0;
                            if(drjet <= rmax && fabs(dzjet) <= rmax) {
                                atomicAdd(&vol_gas[ai], vloc);
                                double psy = exp(-drjet*drjet / (2.0*rmax*rmax));
                                atomicAdd(&psy_norm[ai], psy);
                            }
                            double drc = fmax(fabs(dxx), fmax(fabs(dyy), fabs(dzz)));
                            if(drc <= dxloc * 0.5) {
                                ind_blast[ai] = idx;
                                vol_blast[ai] = vloc;
                                mass_blast[ai] = vloc * dval;
                            }
                        } else {
                            // No spin → spherical
                            if(dr2 <= rmax2) {
                                atomicAdd(&vol_gas[ai], vloc);
                            }
                        }
                    }
                } else {
                    // Thermal case
                    if(dr2 <= rmax2) {
                        atomicAdd(&vol_gas[ai], vloc * dval);
                        atomicAdd(&mass_gas[ai], vloc * dval);
                    }
                    double drc = fmax(fabs(dxx), fmax(fabs(dyy), fabs(dzz)));
                    if(drc <= dxloc * 0.5) {
                        ind_blast[ai] = idx;
                        vol_blast[ai] = vloc;
                        mass_blast[ai] = vloc * dval;
                    }
                }
            }
            iA = agn_next[ai];  // next 1-based value
        }
    }
}

// ==========================================================================
// AGN_blast GPU kernel
// Thread per leaf cell: 27-bin search, modify cell uold, atomicAdd Esave
// cell_uold is column-major: [ncells, nvar], access: cell_uold[tid + ivar*ncells]
// ==========================================================================
__global__ void kernel_AGN_blast(
    const double* __restrict__ cell_x,
    const double* __restrict__ cell_y,
    const double* __restrict__ cell_z,
    const double* __restrict__ cell_vol,
    const double* __restrict__ cell_dx,
    double* cell_uold,      // [ncells * nvar], modified in place
    int ncells, int nvar,
    // AGN data (SoA, nAGN-sized)
    const double* __restrict__ agn_x,
    const double* __restrict__ agn_j,
    const double* __restrict__ agn_v,
    const double* __restrict__ agn_p_gas,
    const double* __restrict__ agn_uBlast,
    const double* __restrict__ agn_mAGN,
    const double* __restrict__ agn_ZAGN,
    const double* __restrict__ agn_psy_norm,
    const double* __restrict__ agn_vol_gas,
    const double* __restrict__ agn_dMBH,
    const double* __restrict__ agn_dMEd,
    const int*    __restrict__ agn_ok_blast,
    const int*    __restrict__ agn_ok_save,
    int nAGN,
    const int*    __restrict__ bin_head,
    const int*    __restrict__ agn_next,
    int nbin, double inv_bin_size, double rmax, double rmax2, double X_floor,
    double scale_T2, double T2maxAGNz, int imetal_0, double gamma_val,
    double* Esave_out)
{
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if(tid >= ncells) return;

    double x = cell_x[tid];
    double y = cell_y[tid];
    double z = cell_z[tid];
    double vloc = cell_vol[tid];

    int ibx = max(1, min(nbin, (int)(x * inv_bin_size) + 1));
    int iby = max(1, min(nbin, (int)(y * inv_bin_size) + 1));
    int ibz = max(1, min(nbin, (int)(z * inv_bin_size) + 1));

    double gm1 = gamma_val - 1.0;

    for(int jbz = max(1, ibz-1); jbz <= min(nbin, ibz+1); jbz++)
    for(int jby = max(1, iby-1); jby <= min(nbin, iby+1); jby++)
    for(int jbx = max(1, ibx-1); jbx <= min(nbin, ibx+1); jbx++) {
        int bidx = (jbx-1) + (jby-1)*nbin + (jbz-1)*nbin*nbin;
        int iA = bin_head[bidx];
        while(iA > 0) {
            int ai = iA - 1;

            double dxx = x - agn_x[ai];
            double dyy = y - agn_x[nAGN + ai];
            double dzz = z - agn_x[2*nAGN + ai];
            double dr2 = dxx*dxx + dyy*dyy + dzz*dzz;

            if(agn_ok_save[ai]) {
                // Case 0: saved energy → thermal deposit
                if(dr2 <= rmax2) {
                    double d = cell_uold[tid + 0*ncells];
                    double ekk = 0.0;
                    for(int dm=0; dm<3; dm++) {
                        double mom = cell_uold[tid + (dm+1)*ncells];
                        ekk += 0.5*mom*mom/d;
                    }
                    double etot = cell_uold[tid + 4*ncells];
                    double eint = etot - ekk;
                    double T2_1 = gm1*eint/d*scale_T2;
                    double pg = agn_p_gas[ai];
                    if(T2_1 < T2maxAGNz) {
                        eint += pg * d;
                        double T2_2 = gm1*eint/d*scale_T2;
                        if(T2_2 <= T2maxAGNz) {
                            cell_uold[tid + 4*ncells] += pg * d;
                        } else {
                            cell_uold[tid + 4*ncells] += T2maxAGNz/scale_T2/gm1*d;
                            atomicAdd(&Esave_out[ai], (T2_2-T2maxAGNz)/scale_T2/gm1*d*vloc);
                        }
                    } else {
                        atomicAdd(&Esave_out[ai], pg*d*vloc);
                    }
                }

            } else if(agn_ok_blast[ai]) {
                double Xr = agn_dMBH[ai] / agn_dMEd[ai];

                if(Xr < X_floor) {
                    // Jet case
                    double jx = agn_j[ai], jy = agn_j[nAGN+ai], jz = agn_j[2*nAGN+ai];
                    double jtot = sqrt(jx*jx + jy*jy + jz*jz);
                    if(jtot > 0.0) {
                        double nj_x=jx/jtot, nj_y=jy/jtot, nj_z=jz/jtot;
                        double dzjet = dxx*nj_x + dyy*nj_y + dzz*nj_z;
                        double drj2 = dr2 - dzjet*dzjet;
                        double drjet = (drj2 > 0.0) ? sqrt(drj2) : 0.0;

                        if(drjet <= rmax && fabs(dzjet) <= rmax) {
                            double pnorm = agn_psy_norm[ai];
                            if(pnorm <= 0.0) { iA = agn_next[ai]; continue; }
                            double psy = exp(-drjet*drjet/(2.0*rmax*rmax)) / pnorm;

                            double u_jet=0, v_jet=0, w_jet=0;
                            if(dzjet < 0.0) {
                                u_jet = -nj_x*agn_uBlast[ai];
                                v_jet = -nj_y*agn_uBlast[ai];
                                w_jet = -nj_z*agn_uBlast[ai];
                            } else if(dzjet > 0.0) {
                                u_jet =  nj_x*agn_uBlast[ai];
                                v_jet =  nj_y*agn_uBlast[ai];
                                w_jet =  nj_z*agn_uBlast[ai];
                            }

                            double d_old = cell_uold[tid + 0*ncells];
                            double ekk = 0.0;
                            for(int dm=0; dm<3; dm++) {
                                double mom = cell_uold[tid + (dm+1)*ncells];
                                ekk += 0.5*mom*mom/d_old;
                            }
                            double etot = cell_uold[tid + 4*ncells];
                            double eint = etot - ekk;
                            double T2_1 = gm1*eint/d_old*scale_T2;

                            double vg = agn_vol_gas[ai];
                            double d_gas = (vg > 0.0) ? agn_mAGN[ai]/vg : 0.0;

                            // Update density
                            cell_uold[tid + 0*ncells] += d_gas * psy;
                            double d_new = cell_uold[tid + 0*ncells];

                            // Update metallicity
                            if(imetal_0 >= 0) {
                                cell_uold[tid + imetal_0*ncells] += agn_ZAGN[ai]*d_gas*psy;
                            }

                            // Velocity = jet + AGN bulk
                            double u = u_jet + agn_v[ai];
                            double v = v_jet + agn_v[nAGN+ai];
                            double w = w_jet + agn_v[2*nAGN+ai];

                            // Update momentum
                            cell_uold[tid + 1*ncells] += d_gas*u*psy;
                            cell_uold[tid + 2*ncells] += d_gas*v*psy;
                            cell_uold[tid + 3*ncells] += d_gas*w*psy;

                            // New kinetic energy
                            ekk = 0.0;
                            for(int dm=0; dm<3; dm++) {
                                double mom = cell_uold[tid + (dm+1)*ncells];
                                ekk += 0.5*mom*mom/d_new;
                            }

                            double pg = agn_p_gas[ai];
                            if(T2_1 >= T2maxAGNz) {
                                double etot_new = etot + 0.5*d_gas*(u*u+v*v+w*w)*psy + pg;
                                cell_uold[tid + 4*ncells] = ekk + eint;
                                double T2_2 = gm1*(etot_new - cell_uold[tid+4*ncells])/d_new*scale_T2;
                                atomicAdd(&Esave_out[ai], T2_2/scale_T2/gm1*d_new*vloc);
                            } else {
                                double etot_new = etot + 0.5*d_gas*(u*u+v*v+w*w)*psy + pg;
                                eint = etot_new - ekk;
                                double T2_2 = gm1*eint/d_new*scale_T2;
                                if(T2_2 <= T2maxAGNz) {
                                    cell_uold[tid + 4*ncells] = ekk + T2_2/scale_T2/gm1*d_new;
                                } else {
                                    cell_uold[tid + 4*ncells] = ekk + T2maxAGNz/scale_T2/gm1*d_new;
                                    atomicAdd(&Esave_out[ai], (T2_2-T2maxAGNz)/scale_T2/gm1*d_new*vloc);
                                }
                            }
                        }
                    } else {
                        // No spin: just add p_gas to energy
                        if(dr2 <= rmax2) {
                            cell_uold[tid + 4*ncells] += agn_p_gas[ai];
                        }
                    }
                } else {
                    // Thermal case
                    if(dr2 <= rmax2) {
                        double d = cell_uold[tid + 0*ncells];
                        double ekk = 0.0;
                        for(int dm=0; dm<3; dm++) {
                            double mom = cell_uold[tid + (dm+1)*ncells];
                            ekk += 0.5*mom*mom/d;
                        }
                        double etot = cell_uold[tid + 4*ncells];
                        double eint = etot - ekk;
                        double T2_1 = gm1*eint/d*scale_T2;
                        double pg = agn_p_gas[ai];
                        if(T2_1 < T2maxAGNz) {
                            eint += pg * d;
                            double T2_2 = gm1*eint/d*scale_T2;
                            if(T2_2 <= T2maxAGNz) {
                                cell_uold[tid + 4*ncells] += pg * d;
                            } else {
                                cell_uold[tid + 4*ncells] += T2maxAGNz/scale_T2/gm1*d;
                                atomicAdd(&Esave_out[ai], (T2_2-T2maxAGNz)/scale_T2/gm1*d*vloc);
                            }
                        } else {
                            atomicAdd(&Esave_out[ai], pg*d*vloc);
                        }
                    }
                }
            }

            iA = agn_next[ai];
        }
    }
}

// ==========================================================================
// Helper: ensure device arrays have enough capacity
// ==========================================================================
static void ensure_cell_arrays(int ncells) {
    if(ncells <= g_cell_cap) return;
    int cap = ncells * 2;  // 2x growth
    if(d_cell_x)   cudaFree(d_cell_x);
    if(d_cell_y)   cudaFree(d_cell_y);
    if(d_cell_z)   cudaFree(d_cell_z);
    if(d_cell_vol) cudaFree(d_cell_vol);
    if(d_cell_d)   cudaFree(d_cell_d);
    if(d_cell_dx)  cudaFree(d_cell_dx);
    if(d_cell_idx) cudaFree(d_cell_idx);
    bool ok = true;
    ok = ok && (cudaMalloc(&d_cell_x,   cap * sizeof(double)) == cudaSuccess);
    ok = ok && (cudaMalloc(&d_cell_y,   cap * sizeof(double)) == cudaSuccess);
    ok = ok && (cudaMalloc(&d_cell_z,   cap * sizeof(double)) == cudaSuccess);
    ok = ok && (cudaMalloc(&d_cell_vol, cap * sizeof(double)) == cudaSuccess);
    ok = ok && (cudaMalloc(&d_cell_d,   cap * sizeof(double)) == cudaSuccess);
    ok = ok && (cudaMalloc(&d_cell_dx,  cap * sizeof(double)) == cudaSuccess);
    ok = ok && (cudaMalloc(&d_cell_idx, cap * sizeof(int))    == cudaSuccess);
    if(!ok) {
        fprintf(stderr, "CUDA SINK: cell array alloc FAILED\n");
        g_cell_cap = 0;
        return;
    }
    g_cell_cap = cap;
}

static void ensure_agn_arrays(int nAGN) {
    if(nAGN <= g_agn_cap) return;
    int cap = max(nAGN * 2, 4096);
    if(d_agn_buf)   cudaFree(d_agn_buf);
    if(d_agn_ibuf)  cudaFree(d_agn_ibuf);
    if(d_agn_next)  cudaFree(d_agn_next);
    if(d_vol_gas)   cudaFree(d_vol_gas);
    if(d_mass_gas)  cudaFree(d_mass_gas);
    if(d_psy_norm)  cudaFree(d_psy_norm);
    if(d_vol_blast) cudaFree(d_vol_blast);
    if(d_mass_blast)cudaFree(d_mass_blast);
    if(d_ind_blast) cudaFree(d_ind_blast);
    if(d_Esave_out) cudaFree(d_Esave_out);
    // AGN doubles: 3*nAGN (x) + 3*nAGN (j) + 3*nAGN (v) + 10*nAGN (scalars) = 19*nAGN
    bool ok = true;
    ok = ok && (cudaMalloc(&d_agn_buf,    cap * 19 * sizeof(double)) == cudaSuccess);
    ok = ok && (cudaMalloc(&d_agn_ibuf,   cap * 3  * sizeof(int))   == cudaSuccess);
    ok = ok && (cudaMalloc(&d_agn_next,   cap * sizeof(int))        == cudaSuccess);
    ok = ok && (cudaMalloc(&d_vol_gas,    cap * sizeof(double))     == cudaSuccess);
    ok = ok && (cudaMalloc(&d_mass_gas,   cap * sizeof(double))     == cudaSuccess);
    ok = ok && (cudaMalloc(&d_psy_norm,   cap * sizeof(double))     == cudaSuccess);
    ok = ok && (cudaMalloc(&d_vol_blast,  cap * sizeof(double))     == cudaSuccess);
    ok = ok && (cudaMalloc(&d_mass_blast, cap * sizeof(double))     == cudaSuccess);
    ok = ok && (cudaMalloc(&d_ind_blast,  cap * sizeof(int))        == cudaSuccess);
    ok = ok && (cudaMalloc(&d_Esave_out,  cap * sizeof(double))     == cudaSuccess);
    if(!ok) {
        fprintf(stderr, "CUDA SINK: AGN array alloc FAILED\n");
        g_agn_cap = 0;
        return;
    }
    g_agn_cap = cap;
}

static void ensure_bin_array(int nbin3) {
    if(nbin3 <= g_bin_cap) return;
    if(d_bin_head) cudaFree(d_bin_head);
    cudaError_t e = cudaMalloc(&d_bin_head, (size_t)nbin3 * sizeof(int));
    if(e != cudaSuccess) {
        fprintf(stderr, "CUDA SINK: bin_head alloc FAILED\n");
        g_bin_cap = 0;
        return;
    }
    g_bin_cap = nbin3;
}

static void ensure_uold_array(int ncells, int nvar) {
    int need = ncells * nvar;
    if(need <= g_uold_cap) return;
    int cap = need * 2;
    if(d_cell_uold) cudaFree(d_cell_uold);
    cudaError_t e = cudaMalloc(&d_cell_uold, (size_t)cap * sizeof(double));
    if(e != cudaSuccess) {
        fprintf(stderr, "CUDA SINK: uold array alloc FAILED (%.1f MB)\n",
                cap * 8.0 / 1e6);
        g_uold_cap = 0;
        return;
    }
    g_uold_cap = cap;
}

// ==========================================================================
// C API: average_AGN on GPU
// ==========================================================================
extern "C" {

int cuda_sink_average_agn(
    const double* cell_x, const double* cell_y, const double* cell_z,
    const double* cell_vol, const double* cell_d, const double* cell_dx,
    const int* cell_idx, int ncells,
    const double* agn_x, const double* agn_j,
    const double* agn_dMBH, const double* agn_dMEd,
    const double* agn_Esave, const int* agn_ok_blast,
    int nAGN,
    const int* bin_head, const int* agn_next,
    int nbin, double inv_bin_size, double rmax, double rmax2, double X_floor,
    double* h_vol_gas, double* h_mass_gas, double* h_psy_norm,
    double* h_vol_blast, double* h_mass_blast, int* h_ind_blast)
{
    if(!is_pool_initialized()) return 0;
    if(ncells <= 0 || nAGN <= 0) return 0;

    // Check available GPU memory
    size_t free_mem = 0, total_mem = 0;
    cudaMemGetInfo(&free_mem, &total_mem);
    size_t need = (size_t)ncells * 7 * 8 + (size_t)nAGN * 25 * 8
                + (size_t)nbin * nbin * nbin * 4;
    if(need > free_mem * 8 / 10) {
        fprintf(stderr, "CUDA SINK avg: SKIP — need %.1f MB, free %.1f MB\n",
                need/1e6, free_mem/1e6);
        return 0;
    }

    // Create stream on first call
    if(!g_sink_stream)
        cudaStreamCreateWithFlags(&g_sink_stream, cudaStreamNonBlocking);

    // Ensure capacity
    int nbin3 = nbin * nbin * nbin;
    ensure_cell_arrays(ncells);
    ensure_agn_arrays(nAGN);
    ensure_bin_array(nbin3);
    if(g_cell_cap < ncells || g_agn_cap < nAGN || g_bin_cap < nbin3)
        return 0;

    // Upload cell data
    cudaMemcpyAsync(d_cell_x,   cell_x,   ncells*sizeof(double), cudaMemcpyHostToDevice, g_sink_stream);
    cudaMemcpyAsync(d_cell_y,   cell_y,   ncells*sizeof(double), cudaMemcpyHostToDevice, g_sink_stream);
    cudaMemcpyAsync(d_cell_z,   cell_z,   ncells*sizeof(double), cudaMemcpyHostToDevice, g_sink_stream);
    cudaMemcpyAsync(d_cell_vol, cell_vol,  ncells*sizeof(double), cudaMemcpyHostToDevice, g_sink_stream);
    cudaMemcpyAsync(d_cell_d,   cell_d,   ncells*sizeof(double), cudaMemcpyHostToDevice, g_sink_stream);
    cudaMemcpyAsync(d_cell_dx,  cell_dx,  ncells*sizeof(double), cudaMemcpyHostToDevice, g_sink_stream);
    cudaMemcpyAsync(d_cell_idx, cell_idx, ncells*sizeof(int),    cudaMemcpyHostToDevice, g_sink_stream);

    // Upload AGN data: x[3*nAGN], j[3*nAGN], dMBH[nAGN], dMEd[nAGN], Esave[nAGN]
    // Pack into d_agn_buf: offsets 0..3n-1=x, 3n..6n-1=j, 6n..7n-1=dMBH, 7n..8n-1=dMEd, 8n..9n-1=Esave
    double* dbuf = d_agn_buf;
    cudaMemcpyAsync(dbuf,             agn_x,     3*nAGN*sizeof(double), cudaMemcpyHostToDevice, g_sink_stream);
    cudaMemcpyAsync(dbuf + 3*nAGN,    agn_j,     3*nAGN*sizeof(double), cudaMemcpyHostToDevice, g_sink_stream);
    cudaMemcpyAsync(dbuf + 6*nAGN,    agn_dMBH,  nAGN*sizeof(double),   cudaMemcpyHostToDevice, g_sink_stream);
    cudaMemcpyAsync(dbuf + 7*nAGN,    agn_dMEd,  nAGN*sizeof(double),   cudaMemcpyHostToDevice, g_sink_stream);
    cudaMemcpyAsync(dbuf + 8*nAGN,    agn_Esave, nAGN*sizeof(double),   cudaMemcpyHostToDevice, g_sink_stream);
    cudaMemcpyAsync(d_agn_ibuf,       agn_ok_blast, nAGN*sizeof(int),   cudaMemcpyHostToDevice, g_sink_stream);

    // Upload spatial hash
    cudaMemcpyAsync(d_bin_head, bin_head, nbin3*sizeof(int),    cudaMemcpyHostToDevice, g_sink_stream);
    cudaMemcpyAsync(d_agn_next, agn_next, nAGN*sizeof(int),    cudaMemcpyHostToDevice, g_sink_stream);

    // Zero output arrays
    cudaMemsetAsync(d_vol_gas,   0, nAGN*sizeof(double), g_sink_stream);
    cudaMemsetAsync(d_mass_gas,  0, nAGN*sizeof(double), g_sink_stream);
    cudaMemsetAsync(d_psy_norm,  0, nAGN*sizeof(double), g_sink_stream);
    cudaMemsetAsync(d_vol_blast, 0, nAGN*sizeof(double), g_sink_stream);
    cudaMemsetAsync(d_mass_blast,0, nAGN*sizeof(double), g_sink_stream);
    // Init ind_blast to -1
    cudaMemsetAsync(d_ind_blast, 0xFF, nAGN*sizeof(int), g_sink_stream);  // -1 in two's complement

    // Launch kernel
    int block = 256;
    int grid  = (ncells + block - 1) / block;
    kernel_average_AGN<<<grid, block, 0, g_sink_stream>>>(
        d_cell_x, d_cell_y, d_cell_z, d_cell_vol, d_cell_d, d_cell_dx, d_cell_idx,
        ncells,
        dbuf, dbuf + 3*nAGN,  // agn_x, agn_j
        dbuf + 6*nAGN, dbuf + 7*nAGN, dbuf + 8*nAGN,  // dMBH, dMEd, Esave
        d_agn_ibuf,  // ok_blast
        nAGN,
        d_bin_head, d_agn_next,
        nbin, inv_bin_size, rmax, rmax2, X_floor,
        d_vol_gas, d_mass_gas, d_psy_norm,
        d_vol_blast, d_mass_blast, d_ind_blast);

    // Download results
    cudaMemcpyAsync(h_vol_gas,   d_vol_gas,   nAGN*sizeof(double), cudaMemcpyDeviceToHost, g_sink_stream);
    cudaMemcpyAsync(h_mass_gas,  d_mass_gas,  nAGN*sizeof(double), cudaMemcpyDeviceToHost, g_sink_stream);
    cudaMemcpyAsync(h_psy_norm,  d_psy_norm,  nAGN*sizeof(double), cudaMemcpyDeviceToHost, g_sink_stream);
    cudaMemcpyAsync(h_vol_blast, d_vol_blast, nAGN*sizeof(double), cudaMemcpyDeviceToHost, g_sink_stream);
    cudaMemcpyAsync(h_mass_blast,d_mass_blast,nAGN*sizeof(double), cudaMemcpyDeviceToHost, g_sink_stream);
    cudaMemcpyAsync(h_ind_blast, d_ind_blast, nAGN*sizeof(int),    cudaMemcpyDeviceToHost, g_sink_stream);

    cudaStreamSynchronize(g_sink_stream);

    cudaError_t err = cudaGetLastError();
    if(err != cudaSuccess) {
        fprintf(stderr, "CUDA SINK avg error: %s\n", cudaGetErrorString(err));
        return 0;
    }
    return 1;
}

// ==========================================================================
// C API: AGN_blast on GPU
// ==========================================================================
int cuda_sink_agn_blast(
    const double* cell_x, const double* cell_y, const double* cell_z,
    const double* cell_vol, const double* cell_dx,
    double* cell_uold, int ncells, int nvar,
    const double* agn_x, const double* agn_j, const double* agn_v,
    const double* agn_p_gas, const double* agn_uBlast,
    const double* agn_mAGN, const double* agn_ZAGN,
    const double* agn_psy_norm, const double* agn_vol_gas,
    const double* agn_dMBH, const double* agn_dMEd,
    const int* agn_ok_blast, const int* agn_ok_save,
    int nAGN,
    const int* bin_head, const int* agn_next,
    int nbin, double inv_bin_size, double rmax, double rmax2, double X_floor,
    double scale_T2, double T2maxAGNz, int imetal, double gamma_val,
    double* h_Esave_out)
{
    if(!is_pool_initialized()) return 0;
    if(ncells <= 0 || nAGN <= 0) return 0;

    // Check available GPU memory
    size_t free_mem = 0, total_mem = 0;
    cudaMemGetInfo(&free_mem, &total_mem);
    long long need_uold = (long long)ncells * nvar * 8;
    size_t need = need_uold + (size_t)ncells * 5 * 8 + (size_t)nAGN * 25 * 8
                + (size_t)nbin * nbin * nbin * 4;
    if(need > free_mem * 8 / 10) {
        fprintf(stderr, "CUDA SINK blast: SKIP — need %.1f MB, free %.1f MB\n",
                need/1e6, free_mem/1e6);
        return 0;
    }

    if(!g_sink_stream)
        cudaStreamCreateWithFlags(&g_sink_stream, cudaStreamNonBlocking);

    int nbin3 = nbin * nbin * nbin;
    ensure_cell_arrays(ncells);
    ensure_agn_arrays(nAGN);
    ensure_bin_array(nbin3);
    ensure_uold_array(ncells, nvar);
    if(g_cell_cap < ncells || g_agn_cap < nAGN || g_bin_cap < nbin3 ||
       g_uold_cap < ncells * nvar)
        return 0;

    // Upload cell position data
    cudaMemcpyAsync(d_cell_x,   cell_x,   ncells*sizeof(double), cudaMemcpyHostToDevice, g_sink_stream);
    cudaMemcpyAsync(d_cell_y,   cell_y,   ncells*sizeof(double), cudaMemcpyHostToDevice, g_sink_stream);
    cudaMemcpyAsync(d_cell_z,   cell_z,   ncells*sizeof(double), cudaMemcpyHostToDevice, g_sink_stream);
    cudaMemcpyAsync(d_cell_vol, cell_vol,  ncells*sizeof(double), cudaMemcpyHostToDevice, g_sink_stream);
    cudaMemcpyAsync(d_cell_dx,  cell_dx,  ncells*sizeof(double), cudaMemcpyHostToDevice, g_sink_stream);

    // Upload cell uold
    cudaMemcpyAsync(d_cell_uold, cell_uold, (long long)ncells*nvar*sizeof(double),
                    cudaMemcpyHostToDevice, g_sink_stream);

    // Upload AGN data: x[3n], j[3n], v[3n], p_gas[n], uBlast[n], mAGN[n], ZAGN[n],
    //                  psy_norm[n], vol_gas[n], dMBH[n], dMEd[n] = 19n doubles
    double* dbuf = d_agn_buf;
    int n = nAGN;
    cudaMemcpyAsync(dbuf,         agn_x,        3*n*sizeof(double), cudaMemcpyHostToDevice, g_sink_stream);
    cudaMemcpyAsync(dbuf + 3*n,   agn_j,        3*n*sizeof(double), cudaMemcpyHostToDevice, g_sink_stream);
    cudaMemcpyAsync(dbuf + 6*n,   agn_v,        3*n*sizeof(double), cudaMemcpyHostToDevice, g_sink_stream);
    cudaMemcpyAsync(dbuf + 9*n,   agn_p_gas,    n*sizeof(double),   cudaMemcpyHostToDevice, g_sink_stream);
    cudaMemcpyAsync(dbuf + 10*n,  agn_uBlast,   n*sizeof(double),   cudaMemcpyHostToDevice, g_sink_stream);
    cudaMemcpyAsync(dbuf + 11*n,  agn_mAGN,     n*sizeof(double),   cudaMemcpyHostToDevice, g_sink_stream);
    cudaMemcpyAsync(dbuf + 12*n,  agn_ZAGN,     n*sizeof(double),   cudaMemcpyHostToDevice, g_sink_stream);
    cudaMemcpyAsync(dbuf + 13*n,  agn_psy_norm, n*sizeof(double),   cudaMemcpyHostToDevice, g_sink_stream);
    cudaMemcpyAsync(dbuf + 14*n,  agn_vol_gas,  n*sizeof(double),   cudaMemcpyHostToDevice, g_sink_stream);
    cudaMemcpyAsync(dbuf + 15*n,  agn_dMBH,     n*sizeof(double),   cudaMemcpyHostToDevice, g_sink_stream);
    cudaMemcpyAsync(dbuf + 16*n,  agn_dMEd,     n*sizeof(double),   cudaMemcpyHostToDevice, g_sink_stream);

    // Upload AGN int arrays: ok_blast[n], ok_save[n]
    cudaMemcpyAsync(d_agn_ibuf,       agn_ok_blast, n*sizeof(int), cudaMemcpyHostToDevice, g_sink_stream);
    cudaMemcpyAsync(d_agn_ibuf + n,   agn_ok_save,  n*sizeof(int), cudaMemcpyHostToDevice, g_sink_stream);

    // Upload spatial hash
    cudaMemcpyAsync(d_bin_head, bin_head, nbin3*sizeof(int), cudaMemcpyHostToDevice, g_sink_stream);
    cudaMemcpyAsync(d_agn_next, agn_next, n*sizeof(int),    cudaMemcpyHostToDevice, g_sink_stream);

    // Zero Esave output
    cudaMemsetAsync(d_Esave_out, 0, n*sizeof(double), g_sink_stream);

    // Launch kernel
    int block = 256;
    int grid  = (ncells + block - 1) / block;
    // imetal is 1-based Fortran, convert to 0-based for kernel
    int imetal_0 = imetal - 1;  // -1 if imetal=0 (no metallicity)

    kernel_AGN_blast<<<grid, block, 0, g_sink_stream>>>(
        d_cell_x, d_cell_y, d_cell_z, d_cell_vol, d_cell_dx,
        d_cell_uold, ncells, nvar,
        dbuf,               // agn_x
        dbuf + 3*n,         // agn_j
        dbuf + 6*n,         // agn_v
        dbuf + 9*n,         // agn_p_gas
        dbuf + 10*n,        // agn_uBlast
        dbuf + 11*n,        // agn_mAGN
        dbuf + 12*n,        // agn_ZAGN
        dbuf + 13*n,        // agn_psy_norm
        dbuf + 14*n,        // agn_vol_gas
        dbuf + 15*n,        // agn_dMBH
        dbuf + 16*n,        // agn_dMEd
        d_agn_ibuf,         // agn_ok_blast
        d_agn_ibuf + n,     // agn_ok_save
        nAGN,
        d_bin_head, d_agn_next,
        nbin, inv_bin_size, rmax, rmax2, X_floor,
        scale_T2, T2maxAGNz, imetal_0, gamma_val,
        d_Esave_out);

    // Download modified uold
    cudaMemcpyAsync(cell_uold, d_cell_uold, (long long)ncells*nvar*sizeof(double),
                    cudaMemcpyDeviceToHost, g_sink_stream);

    // Download Esave
    cudaMemcpyAsync(h_Esave_out, d_Esave_out, n*sizeof(double),
                    cudaMemcpyDeviceToHost, g_sink_stream);

    cudaStreamSynchronize(g_sink_stream);

    cudaError_t err = cudaGetLastError();
    if(err != cudaSuccess) {
        fprintf(stderr, "CUDA SINK blast error: %s\n", cudaGetErrorString(err));
        return 0;
    }
    return 1;
}

// ==========================================================================
// C API: cleanup
// ==========================================================================
void cuda_sink_finalize(void) {
    if(d_cell_x)     { cudaFree(d_cell_x);     d_cell_x = nullptr; }
    if(d_cell_y)     { cudaFree(d_cell_y);     d_cell_y = nullptr; }
    if(d_cell_z)     { cudaFree(d_cell_z);     d_cell_z = nullptr; }
    if(d_cell_vol)   { cudaFree(d_cell_vol);   d_cell_vol = nullptr; }
    if(d_cell_d)     { cudaFree(d_cell_d);     d_cell_d = nullptr; }
    if(d_cell_dx)    { cudaFree(d_cell_dx);    d_cell_dx = nullptr; }
    if(d_cell_idx)   { cudaFree(d_cell_idx);   d_cell_idx = nullptr; }
    if(d_cell_uold)  { cudaFree(d_cell_uold);  d_cell_uold = nullptr; }
    if(d_agn_buf)    { cudaFree(d_agn_buf);    d_agn_buf = nullptr; }
    if(d_agn_ibuf)   { cudaFree(d_agn_ibuf);   d_agn_ibuf = nullptr; }
    if(d_agn_next)   { cudaFree(d_agn_next);   d_agn_next = nullptr; }
    if(d_bin_head)   { cudaFree(d_bin_head);   d_bin_head = nullptr; }
    if(d_vol_gas)    { cudaFree(d_vol_gas);    d_vol_gas = nullptr; }
    if(d_mass_gas)   { cudaFree(d_mass_gas);   d_mass_gas = nullptr; }
    if(d_psy_norm)   { cudaFree(d_psy_norm);   d_psy_norm = nullptr; }
    if(d_vol_blast)  { cudaFree(d_vol_blast);  d_vol_blast = nullptr; }
    if(d_mass_blast) { cudaFree(d_mass_blast); d_mass_blast = nullptr; }
    if(d_ind_blast)  { cudaFree(d_ind_blast);  d_ind_blast = nullptr; }
    if(d_Esave_out)  { cudaFree(d_Esave_out);  d_Esave_out = nullptr; }
    g_cell_cap = 0; g_agn_cap = 0; g_bin_cap = 0; g_uold_cap = 0;
    if(g_sink_stream) {
        cudaStreamSynchronize(g_sink_stream);
        cudaStreamDestroy(g_sink_stream);
        g_sink_stream = nullptr;
    }
}

} // extern "C"
