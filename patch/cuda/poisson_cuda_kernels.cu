// ==========================================================================
// GPU-accelerated Poisson Multigrid kernels for cuRAMSES
// Gauss-Seidel red-black smoother + residual computation with L2 norm
// ==========================================================================

#include "cuda_stream_pool.h"
#include <cstdio>
#include <cstring>
#include <cstdlib>

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
    // Free halo exchange arrays (inline to avoid forward declaration)
    if (d_halo_emit_cells) { cudaFree(d_halo_emit_cells); d_halo_emit_cells = nullptr; }
    if (d_halo_recv_cells) { cudaFree(d_halo_recv_cells); d_halo_recv_cells = nullptr; }
    if (d_halo_emit_buf)   { cudaFree(d_halo_emit_buf);   d_halo_emit_buf   = nullptr; }
    if (d_halo_recv_buf)   { cudaFree(d_halo_recv_buf);   d_halo_recv_buf   = nullptr; }
    if (h_halo_emit_buf)   { cudaFreeHost(h_halo_emit_buf); h_halo_emit_buf = nullptr; }
    if (h_halo_recv_buf)   { cudaFreeHost(h_halo_recv_buf); h_halo_recv_buf = nullptr; }
    g_halo_n_emit = 0;
    g_halo_n_recv = 0;

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
