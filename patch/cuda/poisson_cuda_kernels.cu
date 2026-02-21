// ==========================================================================
// GPU-accelerated Poisson Multigrid kernels for cuRAMSES
// Gauss-Seidel red-black smoother + residual computation with L2 norm
// ==========================================================================

#include "cuda_stream_pool.h"
#include <cufft.h>
#include <cufftXt.h>
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <cmath>

#ifdef USE_CUFFTMP
#include <mpi.h>
// Declare cuFFTMp and NVSHMEM symbols directly to avoid header version mismatch
// (cuFFTMp 11.4 headers expect CUFFT_VERSION 11400, system cuFFT has 12000)
#define CUFFT_COMM_MPI_VAL 0x00
extern "C" cufftResult cufftMpAttachComm(cufftHandle plan, int comm_type, void *comm_handle);

// NVSHMEM removed — causes PMI hang with Intel MPI
#endif

// ==========================================================================
// Lookup tables in constant memory
// For 3D: iii(idim,inbor,ind), jjj(idim,inbor,ind)
// Indices: idim=0..2, inbor=0..1, ind=0..7 (C 0-based)
// iii maps to neighbor grid shift slot (0-6), jjj maps to child cell (1-8)
// ==========================================================================

// iii[idim][inbor][ind] — neighbor grid slot for direction (idim,inbor) at cell ind
// Fortran: iii(1,1,1:8)=(/1,0,1,0,1,0,1,0/) etc.
// C: iii[0][0][0..7] = {1,0,1,0,1,0,1,0}  (0-based idim,inbor; 0-based ind)
__constant__ int d_iii[3][2][8] = {
    {{1,0,1,0,1,0,1,0}, {0,2,0,2,0,2,0,2}},  // idim=0 (x)
    {{3,3,0,0,3,3,0,0}, {0,0,4,4,0,0,4,4}},  // idim=1 (y)
    {{5,5,5,5,0,0,0,0}, {0,0,0,0,6,6,6,6}}   // idim=2 (z)
};

// jjj[idim][inbor][ind] — neighbor child cell index (1-based Fortran)
__constant__ int d_jjj[3][2][8] = {
    {{2,1,4,3,6,5,8,7}, {2,1,4,3,6,5,8,7}},  // idim=0
    {{3,4,1,2,7,8,5,6}, {3,4,1,2,7,8,5,6}},  // idim=1
    {{5,6,7,8,1,2,3,4}, {5,6,7,8,1,2,3,4}}   // idim=2
};

// Red/black cell indices for 3D (1-based Fortran)
__constant__ int d_ired[4]   = {1, 4, 6, 7};
__constant__ int d_iblack[4] = {2, 3, 5, 8};

// ==========================================================================
// GPU-resident Poisson MG arrays
// ==========================================================================

static double* d_mg_phi    = nullptr;  // potential
static double* d_mg_f1     = nullptr;  // residual output (f(:,1))
static double* d_mg_f2     = nullptr;  // BC-modified RHS (f(:,2))
static double* d_mg_f3     = nullptr;  // mask (f(:,3))
static int*    d_mg_flag2  = nullptr;  // scan flag
static int*    d_mg_nbor   = nullptr;  // nbor_grid_fine(0:6, 1:ngrid)
static int*    d_mg_igrid  = nullptr;  // active(ilevel)%igrid(1:ngrid)
static double* d_mg_partial_norm2 = nullptr; // partial norm2 for reduction
static long long g_mg_ncell = 0;
static int     g_mg_alloc_ngrid = 0;  // allocated capacity for grid arrays
static bool    g_mg_ready  = false;

// Pinned host buffer for partial norm2 reduction
static double* h_mg_partial_norm2 = nullptr;
static int     h_mg_partial_cap   = 0;

// Dedicated stream for MG operations
static cudaStream_t g_mg_stream = nullptr;

// Forward declarations for restrict/interp arrays (used by cuda_mg_finalize)
// Constant memory for interpolation coefficients
__constant__ double d_bbb_interp[8];
__constant__ int d_ccc_interp[8][8];  // [ind_average][ind_f], 1-based values

// Device arrays for restrict/interp mappings
static int*    d_restrict_target    = nullptr;
static int*    d_interp_nbor_flat   = nullptr;
static double* d_coarse_rhs_flat    = nullptr;
static double* d_coarse_phi_flat    = nullptr;

// Pinned host buffers for coarse data transfer
static double* h_coarse_rhs_pinned  = nullptr;
static double* h_coarse_phi_pinned  = nullptr;

static int     g_ri_ngrid  = 0;
static int     g_ri_ncells = 0;

// ==========================================================================
// GPU-resident halo exchange arrays for MG phi
// ==========================================================================
static int*    d_halo_emit_cells = nullptr;  // emission cell indices (1-based)
static int*    d_halo_recv_cells = nullptr;  // reception cell indices (1-based)
static double* d_halo_emit_buf   = nullptr;  // GPU staging for gather
static double* d_halo_recv_buf   = nullptr;  // GPU staging for scatter
static double* h_halo_emit_buf   = nullptr;  // pinned host staging for gather
static double* h_halo_recv_buf   = nullptr;  // pinned host staging for scatter
static int     g_halo_n_emit = 0;
static int     g_halo_n_recv = 0;

// ==========================================================================
// cuFFT Poisson solver state (forward declarations for cuda_mg_finalize)
// ==========================================================================
static cufftHandle g_fft_plan_d2z = 0;
static cufftHandle g_fft_plan_z2d = 0;
static cufftDoubleReal*    d_fft_real    = nullptr;
static cufftDoubleComplex* d_fft_complex = nullptr;
static double*             d_fft_green   = nullptr;
static int*                d_fft_map     = nullptr;
static int g_fft_Nx = 0, g_fft_Ny = 0, g_fft_Nz = 0;

// ==========================================================================
// cuFFTMp distributed Poisson solver state
// ==========================================================================
#ifdef USE_CUFFTMP
static cufftHandle    g_fftmp_plan = 0;
static cudaLibXtDesc *g_fftmp_desc = nullptr;
static int g_fftmp_Nx = 0, g_fftmp_Ny = 0, g_fftmp_Nz = 0;
static int g_fftmp_initialized = 0;
static double g_fftmp_dx2 = 0.0;
#endif

// ==========================================================================
// Gauss-Seidel kernel: one thread per active grid, processes 4 cells
// ==========================================================================
__global__ void gauss_seidel_mg_fine_kernel(
    double* __restrict__ phi,
    const double* __restrict__ f2,
    const double* __restrict__ f3,
    const int* __restrict__ flag2,
    const int* __restrict__ nbor_grid,  // (7, ngrid) — nbor_grid[igrid_mg*7 + j]
    const int* __restrict__ igrid_arr,
    int ngrid, int ngridmax, int ncoarse,
    double dx2, int color, int safe_mode_flag)
{
    // NOTE: All cell/grid indices (icell, igrid_amr, etc.) are 1-based Fortran indices.
    // C arrays are 0-based, so we subtract 1 when accessing: arr[idx - 1].

    int igrid_mg = blockIdx.x * blockDim.x + threadIdx.x;
    if (igrid_mg >= ngrid) return;

    int igrid_amr = igrid_arr[igrid_mg]; // 1-based Fortran grid index

    // Load pre-computed neighbor grids into registers (0=self, 1..6=neighbors)
    int nbors[7];
    for (int j = 0; j < 7; j++)
        nbors[j] = nbor_grid[igrid_mg * 7 + j]; // 1-based or 0 (no neighbor)

    // Process 4 cells (red or black)
    for (int ind0 = 0; ind0 < 4; ind0++) {
        int ind = (color == 0) ? d_ired[ind0] : d_iblack[ind0]; // 1-based
        int iskip = ncoarse + (ind - 1) * ngridmax;
        int icell = iskip + igrid_amr; // 1-based Fortran cell index

        double nb_sum = 0.0;

        // Check scan flag
        int scan_val = flag2[icell - 1] / ngridmax;
        if (scan_val == 0) {
            // Fast path: interior cell
            for (int idim = 0; idim < 3; idim++) {
                for (int inbor = 0; inbor < 2; inbor++) {
                    int igshift = d_iii[idim][inbor][ind-1];
                    int jj = d_jjj[idim][inbor][ind-1]; // 1-based
                    int igrid_nbor = nbors[igshift]; // 1-based
                    int icell_nbor = igrid_nbor + ncoarse + (jj - 1) * ngridmax; // 1-based
                    nb_sum += phi[icell_nbor - 1];
                }
            }
            phi[icell - 1] = (nb_sum - dx2 * f2[icell - 1]) / 6.0;
        } else {
            // Slow path: boundary cell
            double mask_c = f3[icell - 1];
            if (mask_c <= 0.0) continue;
            if (safe_mode_flag && mask_c < 1.0) continue;

            double weight = 0.0;
            for (int idim = 0; idim < 3; idim++) {
                for (int inbor = 0; inbor < 2; inbor++) {
                    int igshift = d_iii[idim][inbor][ind-1];
                    int jj = d_jjj[idim][inbor][ind-1];
                    int igrid_nbor = nbors[igshift];

                    if (igrid_nbor == 0) {
                        // No neighbor cell
                        weight -= 1.0 / mask_c;
                    } else {
                        int icell_nbor = igrid_nbor + ncoarse + (jj - 1) * ngridmax;
                        double mask_nbor = f3[icell_nbor - 1];
                        if (mask_nbor <= 0.0) {
                            weight += mask_nbor / mask_c;
                        } else {
                            nb_sum += phi[icell_nbor - 1];
                        }
                    }
                }
            }
            phi[icell - 1] = (nb_sum - dx2 * f2[icell - 1]) / (6.0 - weight);
        }
    }
}

// ==========================================================================
// Residual kernel: one thread per active grid, processes 8 cells
// Optional L2 norm reduction via shared memory
// ==========================================================================
__global__ void cmp_residual_mg_fine_kernel(
    const double* __restrict__ phi,
    double* __restrict__ f1,
    const double* __restrict__ f2,
    const double* __restrict__ f3,
    const int* __restrict__ flag2,
    const int* __restrict__ nbor_grid,
    const int* __restrict__ igrid_arr,
    int ngrid, int ngridmax, int ncoarse,
    double oneoverdx2, double dtwondim,
    double* __restrict__ partial_norm2,
    int compute_norm)
{
    extern __shared__ double sdata[];

    int igrid_mg = blockIdx.x * blockDim.x + threadIdx.x;
    double my_norm2 = 0.0;

    // NOTE: All cell/grid indices (icell, igrid_amr, etc.) are 1-based Fortran indices.
    // C arrays are 0-based, so we subtract 1 when accessing: arr[idx - 1].

    if (igrid_mg < ngrid) {
        int igrid_amr = igrid_arr[igrid_mg]; // 1-based Fortran grid index

        // Load neighbor grids
        int nbors[7];
        for (int j = 0; j < 7; j++)
            nbors[j] = nbor_grid[igrid_mg * 7 + j]; // 1-based or 0 (no neighbor)

        // Process all 8 cells
        for (int ind = 1; ind <= 8; ind++) {
            int iskip = ncoarse + (ind - 1) * ngridmax;
            int icell = iskip + igrid_amr; // 1-based Fortran cell index

            double phi_c = phi[icell - 1];
            double nb_sum = 0.0;
            double mask_c = f3[icell - 1];

            if (flag2[icell - 1] / ngridmax == 0) {
                // Fast path: interior cell
                for (int idim = 0; idim < 3; idim++) {
                    for (int inbor = 0; inbor < 2; inbor++) {
                        int igshift = d_iii[idim][inbor][ind-1];
                        int jj = d_jjj[idim][inbor][ind-1]; // 1-based
                        int igrid_nbor = nbors[igshift]; // 1-based
                        int icell_nbor = igrid_nbor + ncoarse + (jj - 1) * ngridmax; // 1-based
                        nb_sum += phi[icell_nbor - 1];
                    }
                }
            } else {
                // Slow path: boundary cell
                if (mask_c <= 0.0) {
                    f1[icell - 1] = 0.0;
                    continue;
                }
                for (int idim = 0; idim < 3; idim++) {
                    for (int inbor = 0; inbor < 2; inbor++) {
                        int igshift = d_iii[idim][inbor][ind-1];
                        int jj = d_jjj[idim][inbor][ind-1];
                        int igrid_nbor = nbors[igshift];

                        if (igrid_nbor == 0) {
                            nb_sum -= phi_c / mask_c;
                        } else {
                            int icell_nbor = igrid_nbor + ncoarse + (jj - 1) * ngridmax; // 1-based
                            double mask_nbor = f3[icell_nbor - 1];
                            if (mask_nbor <= 0.0) {
                                nb_sum += phi_c * (mask_nbor / mask_c);
                            } else {
                                nb_sum += phi[icell_nbor - 1];
                            }
                        }
                    }
                }
            }

            double res = -oneoverdx2 * (nb_sum - dtwondim * phi_c) + f2[icell - 1];
            f1[icell - 1] = res;

            if (compute_norm && mask_c > 0.0) {
                my_norm2 += res * res;
            }
        }
    }

    // Block-level reduction for norm2
    if (compute_norm) {
        sdata[threadIdx.x] = my_norm2;
        __syncthreads();
        for (int s = blockDim.x / 2; s > 0; s >>= 1) {
            if (threadIdx.x < s)
                sdata[threadIdx.x] += sdata[threadIdx.x + s];
            __syncthreads();
        }
        if (threadIdx.x == 0)
            partial_norm2[blockIdx.x] = sdata[0];
    }
}

// ==========================================================================
// Halo gather kernel: extract emission cells from phi into compact buffer
// ==========================================================================
__global__ void mg_halo_gather_kernel(
    const double* __restrict__ phi,
    double* __restrict__ emit_buf,
    const int* __restrict__ emit_cells,
    int n)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < n) {
        emit_buf[i] = phi[emit_cells[i] - 1]; // 1-based Fortran → 0-based C
    }
}

// ==========================================================================
// Halo scatter kernel: write received values into phi at reception cells
// ==========================================================================
__global__ void mg_halo_scatter_kernel(
    double* __restrict__ phi,
    const double* __restrict__ recv_buf,
    const int* __restrict__ recv_cells,
    int n)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < n) {
        phi[recv_cells[i] - 1] = recv_buf[i]; // 1-based Fortran → 0-based C
    }
}

// ==========================================================================
// C API for Fortran binding
// ==========================================================================

extern "C" {

void cuda_mg_upload(const double* phi, const double* f,
                    const int* flag2, long long ncell,
                    const int* nbor_grid, const int* igrid, int ngrid)
{
    if (!is_pool_initialized()) return;

    // Create dedicated MG stream on first call
    if (!g_mg_stream) {
        cudaStreamCreateWithFlags(&g_mg_stream, cudaStreamNonBlocking);
    }

    // Reallocate cell arrays if needed
    if (ncell != g_mg_ncell) {
        if (d_mg_phi)   { cudaFree(d_mg_phi);   d_mg_phi   = nullptr; }
        if (d_mg_f1)    { cudaFree(d_mg_f1);    d_mg_f1    = nullptr; }
        if (d_mg_f2)    { cudaFree(d_mg_f2);    d_mg_f2    = nullptr; }
        if (d_mg_f3)    { cudaFree(d_mg_f3);    d_mg_f3    = nullptr; }
        if (d_mg_flag2) { cudaFree(d_mg_flag2); d_mg_flag2 = nullptr; }

        // Check available GPU memory
        size_t free_mem = 0, total_mem = 0;
        cudaMemGetInfo(&free_mem, &total_mem);
        size_t need = (size_t)ncell * 4 * sizeof(double) + (size_t)ncell * sizeof(int);
        if (need > free_mem * 9 / 10) {  // Leave 10% margin
            fprintf(stderr, "CUDA MG: SKIP — need %.1f GB but only %.1f GB free\n",
                    (double)need / (1024.0*1024.0*1024.0),
                    (double)free_mem / (1024.0*1024.0*1024.0));
            g_mg_ncell = 0;
            g_mg_ready = false;
            return;
        }

        cudaError_t e1 = cudaMalloc(&d_mg_phi,   (size_t)ncell * sizeof(double));
        cudaError_t e2 = cudaMalloc(&d_mg_f1,    (size_t)ncell * sizeof(double));
        cudaError_t e3 = cudaMalloc(&d_mg_f2,    (size_t)ncell * sizeof(double));
        cudaError_t e4 = cudaMalloc(&d_mg_f3,    (size_t)ncell * sizeof(double));
        cudaError_t e5 = cudaMalloc(&d_mg_flag2, (size_t)ncell * sizeof(int));

        if (e1 || e2 || e3 || e4 || e5) {
            fprintf(stderr, "CUDA MG: cell array allocation FAILED (%.1f GB)\n",
                    (double)need / (1024.0*1024.0*1024.0));
            if (d_mg_phi)   { cudaFree(d_mg_phi);   d_mg_phi   = nullptr; }
            if (d_mg_f1)    { cudaFree(d_mg_f1);    d_mg_f1    = nullptr; }
            if (d_mg_f2)    { cudaFree(d_mg_f2);    d_mg_f2    = nullptr; }
            if (d_mg_f3)    { cudaFree(d_mg_f3);    d_mg_f3    = nullptr; }
            if (d_mg_flag2) { cudaFree(d_mg_flag2); d_mg_flag2 = nullptr; }
            g_mg_ncell = 0;
            g_mg_ready = false;
            return;
        }
        g_mg_ncell = ncell;
    }

    // Reallocate grid arrays if needed
    if (ngrid > g_mg_alloc_ngrid) {
        if (d_mg_nbor)  { cudaFree(d_mg_nbor);  d_mg_nbor  = nullptr; }
        if (d_mg_igrid) { cudaFree(d_mg_igrid); d_mg_igrid = nullptr; }

        int cap = ngrid * 2;  // 2x over-allocation
        cudaError_t e1 = cudaMalloc(&d_mg_nbor,  (size_t)cap * 7 * sizeof(int));
        cudaError_t e2 = cudaMalloc(&d_mg_igrid, (size_t)cap * sizeof(int));
        if (e1 || e2) {
            fprintf(stderr, "CUDA MG: grid array allocation FAILED\n");
            if (d_mg_nbor)  { cudaFree(d_mg_nbor);  d_mg_nbor  = nullptr; }
            if (d_mg_igrid) { cudaFree(d_mg_igrid); d_mg_igrid = nullptr; }
            g_mg_alloc_ngrid = 0;
            g_mg_ready = false;
            return;
        }
        g_mg_alloc_ngrid = cap;
    }
    // Upload cell data (phi, f(:,1..3), flag2)
    cudaMemcpy(d_mg_phi,   phi,   (size_t)ncell * sizeof(double), cudaMemcpyHostToDevice);
    // f is Fortran f(1:ncell, 1:3) — column-major
    // f(:,1) = f[0..ncell-1], f(:,2) = f[ncell..2*ncell-1], f(:,3) = f[2*ncell..3*ncell-1]
    cudaMemcpy(d_mg_f1,    f,             (size_t)ncell * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_mg_f2,    f + ncell,     (size_t)ncell * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_mg_f3,    f + 2 * ncell, (size_t)ncell * sizeof(double), cudaMemcpyHostToDevice);
    // flag2 is Fortran flag2(0:ncell) — starts at index 0, unlike phi/f which start at 1.
    // Skip element 0 so d_mg_flag2[icell-1] = Fortran flag2(icell), matching phi/f indexing.
    cudaMemcpy(d_mg_flag2, flag2 + 1, (size_t)ncell * sizeof(int), cudaMemcpyHostToDevice);

    // (debug prints removed — use cuda_mg_check_flag2 from Fortran instead)

    // Upload grid data (nbor_grid_fine, active igrid)
    // nbor_grid_fine is Fortran (0:twondim, 1:ngrid) = (0:6, 1:ngrid) -> 7 × ngrid
    // Fortran column-major: nbor_grid_fine(j, igrid) at offset j + igrid_mg*7
    // We need row-major: nbor_grid[igrid_mg * 7 + j]
    // Since Fortran stores (0:6, 1:ngrid) column-major, element (j, igrid) is at j + (igrid-1)*7
    // This matches our C layout if we copy directly: nbor_grid[igrid_mg*7 + j] = nbor_grid_fine(j, igrid_mg+1)
    // BUT Fortran column-major stores as: (0,1), (1,1), ..., (6,1), (0,2), (1,2), ..., (6,2), ...
    // which maps to C as: nbor_grid[(igrid_mg)*7 + j] = Fortran nbor_grid_fine(j, igrid_mg+1)
    // So direct memcpy works!
    cudaMemcpy(d_mg_nbor,  nbor_grid, (size_t)ngrid * 7 * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_mg_igrid, igrid,     (size_t)ngrid * sizeof(int), cudaMemcpyHostToDevice);

    g_mg_ready = true;
    cudaDeviceSynchronize();

    cudaError_t err = cudaGetLastError();
    if (err != cudaSuccess) {
        fprintf(stderr, "CUDA MG upload error: %s\n", cudaGetErrorString(err));
        g_mg_ready = false;
    }
}

void cuda_mg_download_phi(double* phi, long long ncell)
{
    if (!g_mg_ready || ncell != g_mg_ncell) return;
    cudaMemcpy(phi, d_mg_phi, (size_t)ncell * sizeof(double), cudaMemcpyDeviceToHost);
}

void cuda_mg_upload_phi(const double* phi, long long ncell)
{
    if (!g_mg_ready || ncell != g_mg_ncell) return;
    cudaMemcpy(d_mg_phi, phi, (size_t)ncell * sizeof(double), cudaMemcpyHostToDevice);
}


void cuda_mg_download_f1(double* f1, long long ncell)
{
    if (!g_mg_ready || ncell != g_mg_ncell) return;
    cudaMemcpy(f1, d_mg_f1, (size_t)ncell * sizeof(double), cudaMemcpyDeviceToHost);
}

void cuda_mg_gauss_seidel(int ngrid, int ngridmax, int ncoarse,
                          double dx2, int color, int safe_mode)
{
    if (!g_mg_ready || ngrid <= 0) return;

    int block = 256;
    int grid = (ngrid + block - 1) / block;

    gauss_seidel_mg_fine_kernel<<<grid, block, 0, g_mg_stream>>>(
        d_mg_phi, d_mg_f2, d_mg_f3, d_mg_flag2,
        d_mg_nbor, d_mg_igrid,
        ngrid, ngridmax, ncoarse,
        dx2, color, safe_mode);

    // Sync after kernel (phi is needed by next color sweep)
    cudaStreamSynchronize(g_mg_stream);
}

void cuda_mg_residual(int ngrid, int ngridmax, int ncoarse,
                      double oneoverdx2, double dtwondim, double dx2_norm,
                      double* norm2, int compute_norm)
{
    if (!g_mg_ready || ngrid <= 0) return;

    int block = 256;
    int grid = (ngrid + block - 1) / block;

    // Ensure partial norm2 buffers
    if (compute_norm) {
        if (grid > h_mg_partial_cap) {
            if (d_mg_partial_norm2) cudaFree(d_mg_partial_norm2);
            if (h_mg_partial_norm2) cudaFreeHost(h_mg_partial_norm2);
            int cap = grid * 2;
            cudaMalloc(&d_mg_partial_norm2, (size_t)cap * sizeof(double));
            cudaMallocHost(&h_mg_partial_norm2, (size_t)cap * sizeof(double));
            h_mg_partial_cap = cap;
        }
    }

    size_t smem = compute_norm ? block * sizeof(double) : 0;

    cmp_residual_mg_fine_kernel<<<grid, block, smem, g_mg_stream>>>(
        d_mg_phi, d_mg_f1, d_mg_f2, d_mg_f3, d_mg_flag2,
        d_mg_nbor, d_mg_igrid,
        ngrid, ngridmax, ncoarse,
        oneoverdx2, dtwondim,
        d_mg_partial_norm2, compute_norm);

    cudaStreamSynchronize(g_mg_stream);

    if (compute_norm && norm2) {
        cudaMemcpy(h_mg_partial_norm2, d_mg_partial_norm2,
                   (size_t)grid * sizeof(double), cudaMemcpyDeviceToHost);
        double total = 0.0;
        for (int i = 0; i < grid; i++)
            total += h_mg_partial_norm2[i];
        *norm2 = dx2_norm * total;
    }
}

void cuda_mg_free(void)
{
    // Don't free cell/grid arrays — keep for reuse across levels/steps
    // Only reset state
    g_mg_ready = false;
    // Note: actual GPU memory freed in cuda_mg_finalize or pool_finalize
}

void cuda_mg_finalize(void)
{
    // Free halo exchange arrays
    if (d_halo_emit_cells) { cudaFree(d_halo_emit_cells); d_halo_emit_cells = nullptr; }
    if (d_halo_recv_cells) { cudaFree(d_halo_recv_cells); d_halo_recv_cells = nullptr; }
    if (d_halo_emit_buf)   { cudaFree(d_halo_emit_buf);   d_halo_emit_buf   = nullptr; }
    if (d_halo_recv_buf)   { cudaFree(d_halo_recv_buf);   d_halo_recv_buf   = nullptr; }
    if (h_halo_emit_buf)   { cudaFreeHost(h_halo_emit_buf); h_halo_emit_buf = nullptr; }
    if (h_halo_recv_buf)   { cudaFreeHost(h_halo_recv_buf); h_halo_recv_buf = nullptr; }
    g_halo_n_emit = 0;
    g_halo_n_recv = 0;

    // Free restrict/interp arrays
    if (d_restrict_target)   { cudaFree(d_restrict_target);   d_restrict_target   = nullptr; }
    if (d_interp_nbor_flat)  { cudaFree(d_interp_nbor_flat);  d_interp_nbor_flat  = nullptr; }
    if (d_coarse_rhs_flat)   { cudaFree(d_coarse_rhs_flat);   d_coarse_rhs_flat   = nullptr; }
    if (d_coarse_phi_flat)   { cudaFree(d_coarse_phi_flat);   d_coarse_phi_flat   = nullptr; }
    if (h_coarse_rhs_pinned) { cudaFreeHost(h_coarse_rhs_pinned); h_coarse_rhs_pinned = nullptr; }
    if (h_coarse_phi_pinned) { cudaFreeHost(h_coarse_phi_pinned); h_coarse_phi_pinned = nullptr; }
    g_ri_ngrid  = 0;
    g_ri_ncells = 0;

    // Free cuFFT arrays
    if (g_fft_plan_d2z) { cufftDestroy(g_fft_plan_d2z); g_fft_plan_d2z = 0; }
    if (g_fft_plan_z2d) { cufftDestroy(g_fft_plan_z2d); g_fft_plan_z2d = 0; }
    if (d_fft_real)    { cudaFree(d_fft_real);    d_fft_real    = nullptr; }
    if (d_fft_complex) { cudaFree(d_fft_complex); d_fft_complex = nullptr; }
    if (d_fft_green)   { cudaFree(d_fft_green);   d_fft_green   = nullptr; }
    if (d_fft_map)     { cudaFree(d_fft_map);     d_fft_map     = nullptr; }
    g_fft_Nx = 0; g_fft_Ny = 0; g_fft_Nz = 0;

#ifdef USE_CUFFTMP
    // Free cuFFTMp arrays
    if (g_fftmp_desc) { cufftXtFree(g_fftmp_desc); g_fftmp_desc = nullptr; }
    if (g_fftmp_plan) { cufftDestroy(g_fftmp_plan); g_fftmp_plan = 0; }
    g_fftmp_Nx = 0; g_fftmp_Ny = 0; g_fftmp_Nz = 0;
    g_fftmp_initialized = 0;
#endif

    if (d_mg_phi)   { cudaFree(d_mg_phi);   d_mg_phi   = nullptr; }
    if (d_mg_f1)    { cudaFree(d_mg_f1);    d_mg_f1    = nullptr; }
    if (d_mg_f2)    { cudaFree(d_mg_f2);    d_mg_f2    = nullptr; }
    if (d_mg_f3)    { cudaFree(d_mg_f3);    d_mg_f3    = nullptr; }
    if (d_mg_flag2) { cudaFree(d_mg_flag2); d_mg_flag2 = nullptr; }
    if (d_mg_nbor)  { cudaFree(d_mg_nbor);  d_mg_nbor  = nullptr; }
    if (d_mg_igrid) { cudaFree(d_mg_igrid); d_mg_igrid = nullptr; }
    if (d_mg_partial_norm2) { cudaFree(d_mg_partial_norm2); d_mg_partial_norm2 = nullptr; }
    if (h_mg_partial_norm2) { cudaFreeHost(h_mg_partial_norm2); h_mg_partial_norm2 = nullptr; }
    if (g_mg_stream) {
        cudaStreamSynchronize(g_mg_stream);
        cudaStreamDestroy(g_mg_stream);
        g_mg_stream = nullptr;
    }
    g_mg_ncell = 0;
    g_mg_alloc_ngrid = 0;
    g_mg_ready = false;
    h_mg_partial_cap = 0;
}

int cuda_mg_is_ready(void)
{
    return g_mg_ready ? 1 : 0;
}

// ==========================================================================
// Halo exchange C API
// ==========================================================================

void cuda_mg_halo_setup(const int* emit_cells, int n_emit,
                        const int* recv_cells, int n_recv)
{
    // Free previous allocations
    if (d_halo_emit_cells) { cudaFree(d_halo_emit_cells); d_halo_emit_cells = nullptr; }
    if (d_halo_recv_cells) { cudaFree(d_halo_recv_cells); d_halo_recv_cells = nullptr; }
    if (d_halo_emit_buf)   { cudaFree(d_halo_emit_buf);   d_halo_emit_buf   = nullptr; }
    if (d_halo_recv_buf)   { cudaFree(d_halo_recv_buf);   d_halo_recv_buf   = nullptr; }
    if (h_halo_emit_buf)   { cudaFreeHost(h_halo_emit_buf); h_halo_emit_buf = nullptr; }
    if (h_halo_recv_buf)   { cudaFreeHost(h_halo_recv_buf); h_halo_recv_buf = nullptr; }

    g_halo_n_emit = n_emit;
    g_halo_n_recv = n_recv;

    if (n_emit > 0) {
        cudaMalloc(&d_halo_emit_cells, (size_t)n_emit * sizeof(int));
        cudaMalloc(&d_halo_emit_buf,   (size_t)n_emit * sizeof(double));
        cudaMallocHost(&h_halo_emit_buf, (size_t)n_emit * sizeof(double));
        cudaMemcpy(d_halo_emit_cells, emit_cells,
                   (size_t)n_emit * sizeof(int), cudaMemcpyHostToDevice);
    }
    if (n_recv > 0) {
        cudaMalloc(&d_halo_recv_cells, (size_t)n_recv * sizeof(int));
        cudaMalloc(&d_halo_recv_buf,   (size_t)n_recv * sizeof(double));
        cudaMallocHost(&h_halo_recv_buf, (size_t)n_recv * sizeof(double));
        cudaMemcpy(d_halo_recv_cells, recv_cells,
                   (size_t)n_recv * sizeof(int), cudaMemcpyHostToDevice);
    }

    printf("MG halo setup: n_emit=%d n_recv=%d (%.1f MB each)\n",
           n_emit, n_recv,
           (double)(n_emit > n_recv ? n_emit : n_recv) * sizeof(double) / (1024.0*1024.0));
    fflush(stdout);
}

void cuda_mg_halo_gather(double* h_emit_buf, int n)
{
    if (!g_mg_ready || n <= 0 || !d_halo_emit_cells) return;

    int block = 256;
    int grid = (n + block - 1) / block;

    // All operations on g_mg_stream — strict ordering guaranteed
    mg_halo_gather_kernel<<<grid, block, 0, g_mg_stream>>>(
        d_mg_phi, d_halo_emit_buf, d_halo_emit_cells, n);
    cudaMemcpyAsync(h_halo_emit_buf, d_halo_emit_buf,
                    (size_t)n * sizeof(double), cudaMemcpyDeviceToHost, g_mg_stream);
    cudaStreamSynchronize(g_mg_stream);

    // Copy pinned → user buffer
    memcpy(h_emit_buf, h_halo_emit_buf, (size_t)n * sizeof(double));
}

void cuda_mg_halo_scatter(const double* h_recv_buf, int n)
{
    if (!g_mg_ready || n <= 0 || !d_halo_recv_cells) return;

    // Copy user buffer → pinned
    memcpy(h_halo_recv_buf, h_recv_buf, (size_t)n * sizeof(double));

    int block = 256;
    int grid_blocks = (n + block - 1) / block;

    // All operations on g_mg_stream — strict ordering guaranteed
    cudaMemcpyAsync(d_halo_recv_buf, h_halo_recv_buf,
                    (size_t)n * sizeof(double), cudaMemcpyHostToDevice, g_mg_stream);
    mg_halo_scatter_kernel<<<grid_blocks, block, 0, g_mg_stream>>>(
        d_mg_phi, d_halo_recv_buf, d_halo_recv_cells, n);
    cudaStreamSynchronize(g_mg_stream);
}

void cuda_mg_halo_free(void)
{
    if (d_halo_emit_cells) { cudaFree(d_halo_emit_cells); d_halo_emit_cells = nullptr; }
    if (d_halo_recv_cells) { cudaFree(d_halo_recv_cells); d_halo_recv_cells = nullptr; }
    if (d_halo_emit_buf)   { cudaFree(d_halo_emit_buf);   d_halo_emit_buf   = nullptr; }
    if (d_halo_recv_buf)   { cudaFree(d_halo_recv_buf);   d_halo_recv_buf   = nullptr; }
    if (h_halo_emit_buf)   { cudaFreeHost(h_halo_emit_buf); h_halo_emit_buf = nullptr; }
    if (h_halo_recv_buf)   { cudaFreeHost(h_halo_recv_buf); h_halo_recv_buf = nullptr; }
    g_halo_n_emit = 0;
    g_halo_n_recv = 0;
}

} // extern "C"

// ==========================================================================
// GPU Restriction + Interpolation for MG V-cycle
// Eliminates full phi/f1 D2H/H2D transfers by doing restrict/interp on GPU
// Only small coarse-level data (~20 MB) transferred via PCIe
// ==========================================================================

// Static arrays d_restrict_target, d_interp_nbor_flat, d_coarse_rhs_flat,
// d_coarse_phi_flat, h_coarse_rhs_pinned, h_coarse_phi_pinned,
// g_ri_ngrid, g_ri_ncells declared at top of file (before cuda_mg_finalize).

// ==========================================================================
// Restriction kernel: accumulate fine residual into coarse RHS
// 1 thread per fine grid. Each fine grid has a unique father → no atomics.
// ==========================================================================
__global__ void restrict_residual_fine_kernel(
    const double* __restrict__ f1,
    const double* __restrict__ f3,
    const int* __restrict__ igrid_arr,
    const int* __restrict__ restrict_target,
    double* __restrict__ coarse_rhs_flat,
    int ngrid, int ngridmax, int ncoarse)
{
    int igrid_f_mg = blockIdx.x * blockDim.x + threadIdx.x;
    if (igrid_f_mg >= ngrid) return;

    int target = restrict_target[igrid_f_mg];
    if (target < 0) return;  // coarse cell masked or not in MG

    int igrid_amr = igrid_arr[igrid_f_mg]; // 1-based Fortran grid index

    double sum = 0.0;
    for (int ind_f = 0; ind_f < 8; ind_f++) {
        int icell_f = ncoarse + ind_f * ngridmax + igrid_amr; // 1-based
        if (f3[icell_f - 1] > 0.0) {
            sum += f1[icell_f - 1];
        }
    }

    coarse_rhs_flat[target] = sum / 8.0;
}

// ==========================================================================
// Interpolation kernel: correct fine phi from coarse solution
// 1 thread per fine grid. Loads 27 coarse phi values, applies bbb/ccc weights.
// ==========================================================================
__global__ void interpolate_correct_fine_kernel(
    double* __restrict__ phi,
    const double* __restrict__ f3,
    const int* __restrict__ igrid_arr,
    const int* __restrict__ interp_nbor_flat,
    const double* __restrict__ coarse_phi_flat,
    int ngrid, int ngridmax, int ncoarse)
{
    int igrid_f_mg = blockIdx.x * blockDim.x + threadIdx.x;
    if (igrid_f_mg >= ngrid) return;

    int igrid_amr = igrid_arr[igrid_f_mg]; // 1-based

    // Load 27 coarse phi values from flat array
    double coarse_vals[27];
    for (int j = 0; j < 27; j++) {
        int idx = interp_nbor_flat[igrid_f_mg * 27 + j];
        coarse_vals[j] = (idx >= 0) ? coarse_phi_flat[idx] : 0.0;
    }

    // Process 8 children
    for (int ind_f = 0; ind_f < 8; ind_f++) {
        int icell_f = ncoarse + ind_f * ngridmax + igrid_amr; // 1-based
        if (f3[icell_f - 1] <= 0.0) continue;  // fine cell masked

        double corr = 0.0;
        for (int ind_avg = 0; ind_avg < 8; ind_avg++) {
            int j = d_ccc_interp[ind_avg][ind_f] - 1;  // 0-based index into 27
            corr += d_bbb_interp[ind_avg] * coarse_vals[j];
        }
        phi[icell_f - 1] += corr;
    }
}

// ==========================================================================
// C API for restrict/interp
// ==========================================================================

extern "C" {

void cuda_mg_ri_setup(const int* restrict_target, const int* interp_nbor_flat,
                      int ngrid_fine, int total_coarse_cells,
                      const double* bbb, const int* ccc)
{
    if (!g_mg_stream || !g_mg_ready) return;

    // Free previous allocations
    if (d_restrict_target)   { cudaFree(d_restrict_target);   d_restrict_target   = nullptr; }
    if (d_interp_nbor_flat)  { cudaFree(d_interp_nbor_flat);  d_interp_nbor_flat  = nullptr; }
    if (d_coarse_rhs_flat)   { cudaFree(d_coarse_rhs_flat);   d_coarse_rhs_flat   = nullptr; }
    if (d_coarse_phi_flat)   { cudaFree(d_coarse_phi_flat);   d_coarse_phi_flat   = nullptr; }
    if (h_coarse_rhs_pinned) { cudaFreeHost(h_coarse_rhs_pinned); h_coarse_rhs_pinned = nullptr; }
    if (h_coarse_phi_pinned) { cudaFreeHost(h_coarse_phi_pinned); h_coarse_phi_pinned = nullptr; }

    g_ri_ngrid  = ngrid_fine;
    g_ri_ncells = total_coarse_cells;

    if (ngrid_fine <= 0 || total_coarse_cells <= 0) {
        g_ri_ngrid = 0;
        return;
    }

    // Check GPU memory
    size_t free_mem = 0, total_mem = 0;
    cudaMemGetInfo(&free_mem, &total_mem);
    size_t need = (size_t)ngrid_fine * sizeof(int)                   // restrict_target
                + (size_t)ngrid_fine * 27 * sizeof(int)              // interp_nbor_flat
                + (size_t)total_coarse_cells * 2 * sizeof(double);   // rhs + phi

    if (need > free_mem * 9 / 10) {
        fprintf(stderr, "CUDA MG RI: need %.1f MB but only %.1f MB free — SKIP\n",
                (double)need / (1024.0*1024.0), (double)free_mem / (1024.0*1024.0));
        g_ri_ngrid = 0;
        return;
    }

    // Allocate device arrays
    cudaError_t e1 = cudaMalloc(&d_restrict_target,  (size_t)ngrid_fine * sizeof(int));
    cudaError_t e2 = cudaMalloc(&d_interp_nbor_flat, (size_t)ngrid_fine * 27 * sizeof(int));
    cudaError_t e3 = cudaMalloc(&d_coarse_rhs_flat,  (size_t)total_coarse_cells * sizeof(double));
    cudaError_t e4 = cudaMalloc(&d_coarse_phi_flat,  (size_t)total_coarse_cells * sizeof(double));

    if (e1 || e2 || e3 || e4) {
        fprintf(stderr, "CUDA MG RI: device allocation FAILED (%.1f MB)\n",
                (double)need / (1024.0*1024.0));
        if (d_restrict_target)  { cudaFree(d_restrict_target);  d_restrict_target  = nullptr; }
        if (d_interp_nbor_flat) { cudaFree(d_interp_nbor_flat); d_interp_nbor_flat = nullptr; }
        if (d_coarse_rhs_flat)  { cudaFree(d_coarse_rhs_flat);  d_coarse_rhs_flat  = nullptr; }
        if (d_coarse_phi_flat)  { cudaFree(d_coarse_phi_flat);  d_coarse_phi_flat  = nullptr; }
        g_ri_ngrid = 0;
        return;
    }

    // Allocate pinned host buffers
    cudaMallocHost(&h_coarse_rhs_pinned, (size_t)total_coarse_cells * sizeof(double));
    cudaMallocHost(&h_coarse_phi_pinned, (size_t)total_coarse_cells * sizeof(double));

    // Upload mappings
    cudaMemcpy(d_restrict_target, restrict_target,
               (size_t)ngrid_fine * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_interp_nbor_flat, interp_nbor_flat,
               (size_t)ngrid_fine * 27 * sizeof(int), cudaMemcpyHostToDevice);

    // Upload interpolation coefficients to constant memory
    cudaMemcpyToSymbol(d_bbb_interp, bbb, 8 * sizeof(double));
    cudaMemcpyToSymbol(d_ccc_interp, ccc, 64 * sizeof(int));

    cudaDeviceSynchronize();

    cudaError_t err = cudaGetLastError();
    if (err != cudaSuccess) {
        fprintf(stderr, "CUDA MG RI setup error: %s\n", cudaGetErrorString(err));
        g_ri_ngrid = 0;
        return;
    }

    printf("MG RI setup: ngrid_fine=%d total_coarse=%d (%.1f MB)\n",
           ngrid_fine, total_coarse_cells, (double)need / (1024.0*1024.0));
    fflush(stdout);
}

void cuda_mg_restrict_execute(int ngrid, int ngridmax, int ncoarse)
{
    if (g_ri_ngrid <= 0 || !g_mg_ready || !d_coarse_rhs_flat) return;

    // Clear coarse RHS
    cudaMemsetAsync(d_coarse_rhs_flat, 0,
                    (size_t)g_ri_ncells * sizeof(double), g_mg_stream);

    int block = 256;
    int grid = (ngrid + block - 1) / block;

    restrict_residual_fine_kernel<<<grid, block, 0, g_mg_stream>>>(
        d_mg_f1, d_mg_f3, d_mg_igrid, d_restrict_target,
        d_coarse_rhs_flat,
        ngrid, ngridmax, ncoarse);

    cudaStreamSynchronize(g_mg_stream);
}

void cuda_mg_restrict_download(double* h_coarse_rhs, int total_cells)
{
    if (g_ri_ngrid <= 0 || !d_coarse_rhs_flat || total_cells <= 0) return;

    cudaMemcpyAsync(h_coarse_rhs_pinned, d_coarse_rhs_flat,
                    (size_t)total_cells * sizeof(double),
                    cudaMemcpyDeviceToHost, g_mg_stream);
    cudaStreamSynchronize(g_mg_stream);

    memcpy(h_coarse_rhs, h_coarse_rhs_pinned, (size_t)total_cells * sizeof(double));
}

void cuda_mg_interp_upload(const double* h_coarse_phi, int total_cells)
{
    if (g_ri_ngrid <= 0 || !d_coarse_phi_flat || total_cells <= 0) return;

    memcpy(h_coarse_phi_pinned, h_coarse_phi, (size_t)total_cells * sizeof(double));

    cudaMemcpyAsync(d_coarse_phi_flat, h_coarse_phi_pinned,
                    (size_t)total_cells * sizeof(double),
                    cudaMemcpyHostToDevice, g_mg_stream);
    cudaStreamSynchronize(g_mg_stream);
}

void cuda_mg_interp_execute(int ngrid, int ngridmax, int ncoarse)
{
    if (g_ri_ngrid <= 0 || !g_mg_ready || !d_coarse_phi_flat) return;

    int block = 256;
    int grid = (ngrid + block - 1) / block;

    interpolate_correct_fine_kernel<<<grid, block, 0, g_mg_stream>>>(
        d_mg_phi, d_mg_f3, d_mg_igrid, d_interp_nbor_flat,
        d_coarse_phi_flat,
        ngrid, ngridmax, ncoarse);

    cudaStreamSynchronize(g_mg_stream);
}

int cuda_mg_ri_is_ready(void)
{
    return (g_ri_ngrid > 0) ? 1 : 0;
}

void cuda_mg_ri_free(void)
{
    if (d_restrict_target)   { cudaFree(d_restrict_target);   d_restrict_target   = nullptr; }
    if (d_interp_nbor_flat)  { cudaFree(d_interp_nbor_flat);  d_interp_nbor_flat  = nullptr; }
    if (d_coarse_rhs_flat)   { cudaFree(d_coarse_rhs_flat);   d_coarse_rhs_flat   = nullptr; }
    if (d_coarse_phi_flat)   { cudaFree(d_coarse_phi_flat);   d_coarse_phi_flat   = nullptr; }
    if (h_coarse_rhs_pinned) { cudaFreeHost(h_coarse_rhs_pinned); h_coarse_rhs_pinned = nullptr; }
    if (h_coarse_phi_pinned) { cudaFreeHost(h_coarse_phi_pinned); h_coarse_phi_pinned = nullptr; }
    g_ri_ngrid  = 0;
    g_ri_ncells = 0;
}

} // extern "C" — restrict/interp

// ==========================================================================
// cuFFT Poisson solver for uniform base level
// Direct spectral solve: FFT → Green's function multiply → IFFT
// Replaces iterative MG V-cycle for fully uniform levels with periodic BC.
// ==========================================================================

// (Static cuFFT state declared at top of file with other statics)

// ==========================================================================
// Kernel: gather RHS from RAMSES cells via fft_map into 3D real array
// d_fft_real[idx_3d] = d_rhs[fft_map[idx_3d] - 1]  (fft_map is 1-based)
// ==========================================================================
__global__ void fft_gather_kernel(
    cufftDoubleReal* __restrict__ d_real,
    const double* __restrict__ d_rhs,
    const int* __restrict__ fft_map,
    int N)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= N) return;
    int icell = fft_map[idx];  // 1-based Fortran index
    d_real[idx] = (icell > 0) ? d_rhs[icell - 1] : 0.0;
}

// ==========================================================================
// Kernel: multiply complex spectrum by real Green's function
// d_complex[i] *= d_green[i]  (element-wise complex × real)
// ==========================================================================
__global__ void fft_green_kernel(
    cufftDoubleComplex* __restrict__ d_cplx,
    const double* __restrict__ d_green,
    int N)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= N) return;
    double g = d_green[idx];
    d_cplx[idx].x *= g;
    d_cplx[idx].y *= g;
}

// ==========================================================================
// Kernel: scatter solved phi from 3D real array back to RAMSES cells
// d_phi[fft_map[idx_3d] - 1] = d_fft_real[idx_3d]
// ==========================================================================
__global__ void fft_scatter_kernel(
    double* __restrict__ d_phi,
    const cufftDoubleReal* __restrict__ d_real,
    const int* __restrict__ fft_map,
    int N)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= N) return;
    int icell = fft_map[idx];  // 1-based Fortran index
    if (icell > 0) d_phi[icell - 1] = d_real[idx];
}

// ==========================================================================
// Host helper: compute scaled Green's function on CPU, upload to GPU
// G_scaled[kx,ky,kz] = dx2 / (N^3 * L_k)
// where L_k = 2cos(2pi*kx/Nx) + 2cos(2pi*ky/Ny) + 2cos(2pi*kz/Nz) - 6
// k=0: G_scaled = 0 (zero-mean)
// ==========================================================================
static void compute_green_function(int Nx, int Ny, int Nz, double dx2)
{
    int Nz_complex = Nz / 2 + 1;
    long long N_complex = (long long)Nx * Ny * Nz_complex;
    double N3 = (double)Nx * (double)Ny * (double)Nz;
    double twopi_x = 2.0 * M_PI / (double)Nx;
    double twopi_y = 2.0 * M_PI / (double)Ny;
    double twopi_z = 2.0 * M_PI / (double)Nz;

    double* h_green = (double*)malloc((size_t)N_complex * sizeof(double));
    if (!h_green) {
        fprintf(stderr, "cuFFT: failed to allocate host Green's function\n");
        return;
    }

    for (int kx = 0; kx < Nx; kx++) {
        double cx = 2.0 * cos(twopi_x * kx);
        for (int ky = 0; ky < Ny; ky++) {
            double cxy = cx + 2.0 * cos(twopi_y * ky);
            for (int kz = 0; kz < Nz_complex; kz++) {
                long long idx = (long long)kx * Ny * Nz_complex
                              + (long long)ky * Nz_complex + kz;
                double Lk = cxy + 2.0 * cos(twopi_z * kz) - 6.0;
                if (kx == 0 && ky == 0 && kz == 0) {
                    h_green[idx] = 0.0;  // zero mode
                } else {
                    h_green[idx] = dx2 / (N3 * Lk);
                }
            }
        }
    }

    cudaMemcpy(d_fft_green, h_green, (size_t)N_complex * sizeof(double),
               cudaMemcpyHostToDevice);
    free(h_green);
}

// ==========================================================================
// C API: cuda_fft_poisson_setup
// Upload fft_map, create cuFFT plans, compute Green's function
// Called once per level (or when grid changes)
// ==========================================================================
extern "C" {

void cuda_fft_poisson_setup(const int* fft_map, int Nx, int Ny, int Nz, double dx2)
{
    if (!g_mg_stream) {
        cudaStreamCreateWithFlags(&g_mg_stream, cudaStreamNonBlocking);
    }

    long long N_real = (long long)Nx * Ny * Nz;
    int Nz_complex = Nz / 2 + 1;
    long long N_complex = (long long)Nx * Ny * Nz_complex;

    // Check if plans need recreation (grid size changed)
    if (Nx != g_fft_Nx || Ny != g_fft_Ny || Nz != g_fft_Nz) {
        // Destroy old plans
        if (g_fft_plan_d2z) { cufftDestroy(g_fft_plan_d2z); g_fft_plan_d2z = 0; }
        if (g_fft_plan_z2d) { cufftDestroy(g_fft_plan_z2d); g_fft_plan_z2d = 0; }
        if (d_fft_real)    { cudaFree(d_fft_real);    d_fft_real    = nullptr; }
        if (d_fft_complex) { cudaFree(d_fft_complex); d_fft_complex = nullptr; }
        if (d_fft_green)   { cudaFree(d_fft_green);   d_fft_green   = nullptr; }
        if (d_fft_map)     { cudaFree(d_fft_map);     d_fft_map     = nullptr; }

        // Check GPU memory
        size_t free_mem = 0, total_mem = 0;
        cudaMemGetInfo(&free_mem, &total_mem);
        size_t need = (size_t)N_real * sizeof(cufftDoubleReal)
                    + (size_t)N_complex * sizeof(cufftDoubleComplex)
                    + (size_t)N_complex * sizeof(double)
                    + (size_t)N_real * sizeof(int);
        if (need > free_mem * 8 / 10) {
            fprintf(stderr, "cuFFT: need %.1f MB but only %.1f MB free — SKIP\n",
                    (double)need / (1024.0*1024.0), (double)free_mem / (1024.0*1024.0));
            g_fft_Nx = 0; g_fft_Ny = 0; g_fft_Nz = 0;
            return;
        }

        // Allocate arrays
        cudaError_t e1 = cudaMalloc(&d_fft_real,    (size_t)N_real * sizeof(cufftDoubleReal));
        cudaError_t e2 = cudaMalloc(&d_fft_complex, (size_t)N_complex * sizeof(cufftDoubleComplex));
        cudaError_t e3 = cudaMalloc(&d_fft_green,   (size_t)N_complex * sizeof(double));
        cudaError_t e4 = cudaMalloc(&d_fft_map,     (size_t)N_real * sizeof(int));
        if (e1 || e2 || e3 || e4) {
            fprintf(stderr, "cuFFT: allocation FAILED\n");
            if (d_fft_real)    { cudaFree(d_fft_real);    d_fft_real    = nullptr; }
            if (d_fft_complex) { cudaFree(d_fft_complex); d_fft_complex = nullptr; }
            if (d_fft_green)   { cudaFree(d_fft_green);   d_fft_green   = nullptr; }
            if (d_fft_map)     { cudaFree(d_fft_map);     d_fft_map     = nullptr; }
            g_fft_Nx = 0; g_fft_Ny = 0; g_fft_Nz = 0;
            return;
        }

        // Create cuFFT plans
        int dims[3] = {Nx, Ny, Nz};  // row-major: x varies slowest
        cufftResult r1 = cufftPlanMany(&g_fft_plan_d2z, 3, dims,
                                        NULL, 1, 0, NULL, 1, 0,
                                        CUFFT_D2Z, 1);
        cufftResult r2 = cufftPlanMany(&g_fft_plan_z2d, 3, dims,
                                        NULL, 1, 0, NULL, 1, 0,
                                        CUFFT_Z2D, 1);
        if (r1 != CUFFT_SUCCESS || r2 != CUFFT_SUCCESS) {
            fprintf(stderr, "cuFFT: plan creation FAILED (r1=%d r2=%d)\n", r1, r2);
            if (g_fft_plan_d2z) { cufftDestroy(g_fft_plan_d2z); g_fft_plan_d2z = 0; }
            if (g_fft_plan_z2d) { cufftDestroy(g_fft_plan_z2d); g_fft_plan_z2d = 0; }
            g_fft_Nx = 0; g_fft_Ny = 0; g_fft_Nz = 0;
            return;
        }

        // Set plans to use MG stream
        cufftSetStream(g_fft_plan_d2z, g_mg_stream);
        cufftSetStream(g_fft_plan_z2d, g_mg_stream);

        g_fft_Nx = Nx;
        g_fft_Ny = Ny;
        g_fft_Nz = Nz;

        // Compute and upload Green's function
        compute_green_function(Nx, Ny, Nz, dx2);

        printf("cuFFT setup: %d x %d x %d, dx2=%.6e, GPU mem=%.1f MB\n",
               Nx, Ny, Nz, dx2, (double)need / (1024.0*1024.0));
        fflush(stdout);
    }

    // Always upload fft_map (may change after load balancing)
    cudaMemcpy(d_fft_map, fft_map, (size_t)N_real * sizeof(int),
               cudaMemcpyHostToDevice);
}

// ==========================================================================
// C API: cuda_fft_poisson_solve
// Assumes d_mg_f2 has RHS already uploaded (via cuda_mg_upload),
// and d_fft_map is set via cuda_fft_poisson_setup.
// After solve, d_mg_phi contains the solution.
// ==========================================================================
void cuda_fft_poisson_solve(const double* h_rhs_3d, int N_real_in)
{
    if (g_fft_Nx <= 0 || !d_fft_real || !g_mg_stream) return;

    long long N_real = (long long)g_fft_Nx * g_fft_Ny * g_fft_Nz;
    int Nz_complex = g_fft_Nz / 2 + 1;
    long long N_complex = (long long)g_fft_Nx * g_fft_Ny * Nz_complex;

    // Upload global RHS 3D array to d_fft_real
    cudaMemcpyAsync(d_fft_real, h_rhs_3d, (size_t)N_real * sizeof(double),
                    cudaMemcpyHostToDevice, g_mg_stream);

    // Forward FFT: R2C
    cufftResult r = cufftExecD2Z(g_fft_plan_d2z, d_fft_real, d_fft_complex);
    if (r != CUFFT_SUCCESS) {
        fprintf(stderr, "cuFFT forward FAILED: %d\n", r);
        return;
    }

    // Apply Green's function
    int block = 256;
    int grid = ((int)N_complex + block - 1) / block;
    fft_green_kernel<<<grid, block, 0, g_mg_stream>>>(
        d_fft_complex, d_fft_green, (int)N_complex);

    // Inverse FFT: C2R
    r = cufftExecZ2D(g_fft_plan_z2d, d_fft_complex, d_fft_real);
    if (r != CUFFT_SUCCESS) {
        fprintf(stderr, "cuFFT inverse FAILED: %d\n", r);
        return;
    }

    // Scatter from 3D array to RAMSES phi via fft_map
    grid = ((int)N_real + block - 1) / block;
    fft_scatter_kernel<<<grid, block, 0, g_mg_stream>>>(
        d_mg_phi, d_fft_real, d_fft_map, (int)N_real);

    cudaStreamSynchronize(g_mg_stream);
}

// ==========================================================================
// C API: cuda_fft_poisson_solve_host
// Like cuda_fft_poisson_solve but returns result in host h_phi_3d array.
// Does NOT use d_mg_phi — fully independent of MG GPU arrays.
// ==========================================================================
void cuda_fft_poisson_solve_host(const double* h_rhs_3d, double* h_phi_3d, int N_real_in)
{
    if (g_fft_Nx <= 0 || !d_fft_real || !g_mg_stream) return;

    long long N_real = (long long)g_fft_Nx * g_fft_Ny * g_fft_Nz;
    int Nz_complex = g_fft_Nz / 2 + 1;
    long long N_complex = (long long)g_fft_Nx * g_fft_Ny * Nz_complex;

    // Upload global RHS 3D array to d_fft_real
    cudaMemcpyAsync(d_fft_real, h_rhs_3d, (size_t)N_real * sizeof(double),
                    cudaMemcpyHostToDevice, g_mg_stream);

    // Forward FFT: R2C
    cufftResult r = cufftExecD2Z(g_fft_plan_d2z, d_fft_real, d_fft_complex);
    if (r != CUFFT_SUCCESS) {
        fprintf(stderr, "cuFFT forward FAILED: %d\n", r);
        return;
    }

    // Apply Green's function
    int block = 256;
    int grid = ((int)N_complex + block - 1) / block;
    fft_green_kernel<<<grid, block, 0, g_mg_stream>>>(
        d_fft_complex, d_fft_green, (int)N_complex);

    // Inverse FFT: C2R
    r = cufftExecZ2D(g_fft_plan_z2d, d_fft_complex, d_fft_real);
    if (r != CUFFT_SUCCESS) {
        fprintf(stderr, "cuFFT inverse FAILED: %d\n", r);
        return;
    }

    // Download result to host (no scatter to d_mg_phi)
    cudaMemcpyAsync(h_phi_3d, d_fft_real, (size_t)N_real * sizeof(double),
                    cudaMemcpyDeviceToHost, g_mg_stream);
    cudaStreamSynchronize(g_mg_stream);
}

// ==========================================================================
// C API: cuda_fft_poisson_free
// Free cuFFT plans and arrays. Called at finalize.
// ==========================================================================
void cuda_fft_poisson_free(void)
{
    if (g_fft_plan_d2z) { cufftDestroy(g_fft_plan_d2z); g_fft_plan_d2z = 0; }
    if (g_fft_plan_z2d) { cufftDestroy(g_fft_plan_z2d); g_fft_plan_z2d = 0; }
    if (d_fft_real)    { cudaFree(d_fft_real);    d_fft_real    = nullptr; }
    if (d_fft_complex) { cudaFree(d_fft_complex); d_fft_complex = nullptr; }
    if (d_fft_green)   { cudaFree(d_fft_green);   d_fft_green   = nullptr; }
    if (d_fft_map)     { cudaFree(d_fft_map);     d_fft_map     = nullptr; }
    g_fft_Nx = 0; g_fft_Ny = 0; g_fft_Nz = 0;
}

} // extern "C" — cuFFT Poisson

// ==========================================================================
// cuFFTMp distributed Poisson solver
// Uses NVIDIA cuFFTMp for multi-GPU distributed 3D FFT.
// Eliminates MPI_ALLREDUCE bottleneck by distributing data across GPUs.
// ==========================================================================
#ifdef USE_CUFFTMP

// ==========================================================================
// Kernel: Apply Green's function on distributed complex data (Y-slab)
// After forward D2Z, cuFFTMp shuffles data into Y-slab decomposition.
// Each rank processes its local Y-range of the complex array.
// Layout: [Nx][local_Ny][Nz/2+1] in row-major
// ==========================================================================
__global__ void fftmp_green_kernel(
    cufftDoubleComplex* __restrict__ d_cplx,
    int Nx, int Ny, int Nz,
    int local_Ny, int y_start,
    double dx2)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    int Nz_complex = Nz / 2 + 1;
    long long total = (long long)Nx * local_Ny * Nz_complex;
    if (idx >= total) return;

    // Decompose linear index into (kx, local_ky, kz)
    int kz = idx % Nz_complex;
    int rem = idx / Nz_complex;
    int local_ky = rem % local_Ny;
    int kx = rem / local_Ny;

    int ky = local_ky + y_start;

    // Compute Green's function: G = dx2 / (N^3 * Lk)
    // Lk = 2*cos(2*pi*kx/Nx) + 2*cos(2*pi*ky/Ny) + 2*cos(2*pi*kz/Nz) - 6
    double N3 = (double)Nx * (double)Ny * (double)Nz;
    double twopi_x = 2.0 * M_PI / (double)Nx;
    double twopi_y = 2.0 * M_PI / (double)Ny;
    double twopi_z = 2.0 * M_PI / (double)Nz;

    double Lk = 2.0 * cos(twopi_x * kx)
              + 2.0 * cos(twopi_y * ky)
              + 2.0 * cos(twopi_z * kz) - 6.0;

    double g;
    if (kx == 0 && ky == 0 && kz == 0) {
        g = 0.0;  // zero mode
    } else {
        g = dx2 / (N3 * Lk);
    }

    d_cplx[idx].x *= g;
    d_cplx[idx].y *= g;
}

// ==========================================================================
// C API: cuda_fftmp_poisson_setup
// Create cuFFTMp plan with MPI communicator, allocate distributed memory.
// ==========================================================================
extern "C" {

void cuda_fftmp_poisson_setup(int f_comm, int Nx, int Ny, int Nz, double dx2)
{
    if (g_fftmp_initialized && Nx == g_fftmp_Nx && Ny == g_fftmp_Ny && Nz == g_fftmp_Nz) {
        return;  // already set up for this grid size
    }

    // Cleanup previous state
    if (g_fftmp_initialized) {
        if (g_fftmp_desc) { cufftXtFree(g_fftmp_desc); g_fftmp_desc = nullptr; }
        if (g_fftmp_plan) { cufftDestroy(g_fftmp_plan); g_fftmp_plan = 0; }
        g_fftmp_initialized = 0;
    }

    // Convert Fortran MPI communicator to C
    MPI_Comm c_comm = MPI_Comm_f2c(f_comm);

    int myrank, nranks;
    MPI_Comm_rank(c_comm, &myrank);
    MPI_Comm_size(c_comm, &nranks);

    // Create cuFFT plan (NVSHMEM init removed — causes PMI hang)
    cufftResult r = cufftCreate(&g_fftmp_plan);
    if (r != CUFFT_SUCCESS) {
        if (myrank == 0)
            fprintf(stderr, "cuFFTMp: cufftCreate FAILED: %d\n", r);
        g_fftmp_initialized = -1;  // mark as permanently failed
        return;
    }

    // Attach MPI communicator (from libcufftMp)
    r = cufftMpAttachComm(g_fftmp_plan, CUFFT_COMM_MPI_VAL, &c_comm);
    if (r != CUFFT_SUCCESS) {
        if (myrank == 0)
            fprintf(stderr, "cuFFTMp: cufftMpAttachComm FAILED: %d "
                    "(cuFFT version mismatch? system=%d, cuFFTMp expects 11400)\n",
                    r, CUFFT_VERSION);
        cufftDestroy(g_fftmp_plan); g_fftmp_plan = 0;
        g_fftmp_initialized = -1;  // mark as permanently failed
        return;
    }

    // Make 3D plan (D2Z = real-to-complex)
    size_t workSize = 0;
    r = cufftMakePlan3d(g_fftmp_plan, Nx, Ny, Nz, CUFFT_D2Z, &workSize);
    if (r != CUFFT_SUCCESS) {
        if (myrank == 0)
            fprintf(stderr, "cuFFTMp: cufftMakePlan3d FAILED: %d\n", r);
        cufftDestroy(g_fftmp_plan); g_fftmp_plan = 0;
        g_fftmp_initialized = -1;
        return;
    }

    // Allocate distributed memory (cuFFTMp manages the slab decomposition)
    r = cufftXtMalloc(g_fftmp_plan, &g_fftmp_desc, CUFFT_XT_FORMAT_INPLACE);
    if (r != CUFFT_SUCCESS) {
        if (myrank == 0)
            fprintf(stderr, "cuFFTMp: cufftXtMalloc FAILED: %d\n", r);
        cufftDestroy(g_fftmp_plan); g_fftmp_plan = 0;
        g_fftmp_initialized = -1;
        return;
    }

    g_fftmp_Nx = Nx;
    g_fftmp_Ny = Ny;
    g_fftmp_Nz = Nz;
    g_fftmp_dx2 = dx2;
    g_fftmp_initialized = 1;

    if (myrank == 0) {
        printf("cuFFTMp setup: %d x %d x %d, dx2=%.6e, workSize=%.1f MB, nranks=%d\n",
               Nx, Ny, Nz, dx2, (double)workSize / (1024.0*1024.0), nranks);
        fflush(stdout);
    }
}

// ==========================================================================
// C API: cuda_fftmp_poisson_solve
// Upload local X-slab RHS, forward FFT, Green, inverse FFT, download phi.
// h_rhs_slab: host array [local_Nx][Ny][Nz] (real, row-major)
// h_phi_slab: host array [local_Nx][Ny][Nz] (real, row-major, output)
// local_Nx, Ny, Nz: slab dimensions
// ==========================================================================
void cuda_fftmp_poisson_solve(const double* h_rhs_slab, double* h_phi_slab,
                              int local_Nx, int Ny, int Nz,
                              int y_start_shuffled, int local_Ny_shuffled)
{
    if (!g_fftmp_initialized || !g_fftmp_desc) return;

    // Upload host RHS to distributed descriptor (INPLACE = X-slab layout)
    cufftResult r = cufftXtMemcpy(g_fftmp_plan, g_fftmp_desc,
                                  (void*)h_rhs_slab, CUFFT_COPY_HOST_TO_DEVICE);
    if (r != CUFFT_SUCCESS) {
        fprintf(stderr, "cuFFTMp: H2D memcpy FAILED: %d\n", r);
        return;
    }

    // Forward FFT: D2Z (real-to-complex)
    // After exec, data is in INPLACE_SHUFFLED format (Y-slab decomposition)
    r = cufftXtExecDescriptorD2Z(g_fftmp_plan, g_fftmp_desc, g_fftmp_desc);
    if (r != CUFFT_SUCCESS) {
        fprintf(stderr, "cuFFTMp: forward D2Z FAILED: %d\n", r);
        return;
    }

    // Apply Green's function on GPU (data is now on local GPU in Y-slab)
    {
        int Nz_complex = Nz / 2 + 1;
        long long total = (long long)g_fftmp_Nx * local_Ny_shuffled * Nz_complex;
        int block = 256;
        int grid = ((int)total + block - 1) / block;

        // Get pointer to local GPU data from descriptor
        cudaXtDesc *xtdesc = g_fftmp_desc->descriptor;
        cufftDoubleComplex *d_local = (cufftDoubleComplex*)xtdesc->data[0];

        fftmp_green_kernel<<<grid, block>>>(
            d_local, g_fftmp_Nx, Ny, Nz,
            local_Ny_shuffled, y_start_shuffled, g_fftmp_dx2);
        cudaDeviceSynchronize();
    }

    // Inverse FFT: Z2D (complex-to-real)
    // cuFFTMp inverse with D2Z plan: use cufftXtExecDescriptorZ2D
    r = cufftXtExecDescriptorZ2D(g_fftmp_plan, g_fftmp_desc, g_fftmp_desc);
    if (r != CUFFT_SUCCESS) {
        fprintf(stderr, "cuFFTMp: inverse Z2D FAILED: %d\n", r);
        return;
    }

    // Download result (INPLACE = X-slab) back to host
    r = cufftXtMemcpy(g_fftmp_plan, (void*)h_phi_slab,
                      g_fftmp_desc, CUFFT_COPY_DEVICE_TO_HOST);
    if (r != CUFFT_SUCCESS) {
        fprintf(stderr, "cuFFTMp: D2H memcpy FAILED: %d\n", r);
        return;
    }
}

// ==========================================================================
// C API: cuda_fftmp_poisson_free
// Free cuFFTMp plan and distributed memory.
// ==========================================================================
void cuda_fftmp_poisson_free(void)
{
    if (g_fftmp_desc) { cufftXtFree(g_fftmp_desc); g_fftmp_desc = nullptr; }
    if (g_fftmp_plan) { cufftDestroy(g_fftmp_plan); g_fftmp_plan = 0; }
    g_fftmp_Nx = 0; g_fftmp_Ny = 0; g_fftmp_Nz = 0;
    g_fftmp_initialized = 0;
}

// ==========================================================================
// C API: cuda_fftmp_is_ready
// Returns 1 if cuFFTMp is initialized and ready, 0 otherwise, -1 if failed
// ==========================================================================
int cuda_fftmp_is_ready(void)
{
    return g_fftmp_initialized;
}

// ==========================================================================
// C API: cuda_fftmp_get_local_sizes
// Query the local slab sizes from the plan descriptor.
// For INPLACE (input, X-slab): local_Nx = Nx/nranks (approx)
// For INPLACE_SHUFFLED (output, Y-slab): local_Ny = Ny/nranks (approx)
// Returns via pointers: local_Nx, x_start, local_Ny_shuffled, y_start_shuffled
// ==========================================================================
void cuda_fftmp_get_local_sizes(int* local_Nx_out, int* x_start_out,
                                int* local_Ny_out, int* y_start_out)
{
    if (!g_fftmp_initialized || !g_fftmp_desc) {
        *local_Nx_out = 0; *x_start_out = 0;
        *local_Ny_out = 0; *y_start_out = 0;
        return;
    }

    // For D2Z plan, INPLACE allocation is max of:
    //   real:    local_Nx * Ny * Nz * sizeof(double)
    //   complex: Nx * local_Ny * (Nz/2+1) * sizeof(cufftDoubleComplex)
    // We determine local_Nx and local_Ny from the grid dimensions and rank.

    int myrank, nranks;
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    MPI_Comm_size(MPI_COMM_WORLD, &nranks);

    // cuFFTMp slabs: first dimension (X) for input, second dimension (Y) for output
    // Default decomposition divides evenly with remainder to last ranks
    int Nx = g_fftmp_Nx, Ny = g_fftmp_Ny;

    int base_Nx = Nx / nranks;
    int rem_Nx  = Nx % nranks;
    // Ranks 0..rem_Nx-1 get (base_Nx+1), rest get base_Nx
    if (myrank < rem_Nx) {
        *local_Nx_out = base_Nx + 1;
        *x_start_out  = myrank * (base_Nx + 1);
    } else {
        *local_Nx_out = base_Nx;
        *x_start_out  = rem_Nx * (base_Nx + 1) + (myrank - rem_Nx) * base_Nx;
    }

    int base_Ny = Ny / nranks;
    int rem_Ny  = Ny % nranks;
    if (myrank < rem_Ny) {
        *local_Ny_out = base_Ny + 1;
        *y_start_out  = myrank * (base_Ny + 1);
    } else {
        *local_Ny_out = base_Ny;
        *y_start_out  = rem_Ny * (base_Ny + 1) + (myrank - rem_Ny) * base_Ny;
    }
}

} // extern "C" — cuFFTMp
#endif
