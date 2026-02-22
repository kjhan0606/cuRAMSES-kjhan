// ==========================================================================
// Unified CUDA stream pool for cuRAMSES
// Shared header for all CUDA kernel files (hydro, etc.)
// Thread-safe acquire/release for OpenMP concurrent access
// ==========================================================================
#ifndef CUDA_STREAM_POOL_H
#define CUDA_STREAM_POOL_H

#include <cuda_runtime.h>

// Compile-time constants (passed via -D flags from Makefile)
#ifndef NDIM
#define NDIM 3
#endif
#ifndef NVAR
#define NVAR 11
#endif
#ifndef NVECTOR
#define NVECTOR 32
#endif

#ifndef N_STREAMS
#define N_STREAMS 1
#endif

// Stencil ranges (matches hydro_parameters.f90: iu1:iu2 = -1:4)
#define IU1 (-1)
#define IU2 (4)
#define JU1 (-1)
#define JU2 (4)
#define KU1 (-1)
#define KU2 (4)
#define STENCIL_NI (IU2 - IU1 + 1)   // 6
#define STENCIL_NJ (JU2 - JU1 + 1)   // 6
#define STENCIL_NK (KU2 - KU1 + 1)   // 6

// Flux ranges (if1:if2 = 1:3)
#define IF1 1
#define IF2 3
#define JF1 1
#define JF2 3
#define KF1 1
#define KF2 3
#define FLUX_NI (IF2 - IF1 + 1)   // 3
#define FLUX_NJ (JF2 - JF1 + 1)   // 3
#define FLUX_NK (KF2 - KF1 + 1)   // 3

#define NCELLS_PER_GRID (STENCIL_NI * STENCIL_NJ * STENCIL_NK)  // 216
#define NFLUX_PER_GRID  (FLUX_NI * FLUX_NJ * FLUX_NK)           // 27

// Profiling: 8 events mark boundaries of 7 phases
// [0]H2D[1]ctoprim[2]uslope[3]trace3d[4]flux[5]difmag[6]D2H[7]
#define N_PROFILE_EVENTS 8

// Per-stream buffer slot for hydro computation
typedef struct {
    cudaStream_t stream;
    volatile int busy;

    // Hydro I/O buffers (host ↔ device transfer)
    double *d_uloc, *d_gloc, *d_flux, *d_tmp;
    int *d_ok;
    int hydro_cap;       // allocated grid capacity

    // Hydro intermediate arrays (device only)
    double *d_q, *d_c, *d_dq, *d_qm, *d_qp;
    int hydro_inter_cap; // allocated grid capacity

    // GPU-gather stencil index buffers (device only)
    int *d_stencil_idx;     // (stride, NI, NJ, NK)
    int *d_stencil_grav;    // (stride, NI, NJ, NK)
    double *d_interp_vals;  // (n_interp, NVAR)
    int stencil_cap;        // allocated stride capacity
    int interp_cap;         // allocated interp capacity

    // Scatter-reduce output buffers (device only)
    int *d_ok_int;          // (stride, NI, NJ, NK) - ok flags for flux reset
    double *d_add_unew;     // (stride, 8, NVAR) - Level L scatter result
    double *d_add_lm1;      // (stride, 6, NVAR) - Level L-1 scatter result
    double *d_add_divu_l;   // (stride, 8) - pressure_fix divu Level L
    double *d_add_enew_l;   // (stride, 8) - pressure_fix enew Level L
    double *d_add_divu_lm1; // (stride, 6) - pressure_fix divu Level L-1
    double *d_add_enew_lm1; // (stride, 6) - pressure_fix enew Level L-1
    int reduce_cap;         // allocated grid capacity for reduce buffers

    // Pinned host staging buffers for fast H2D DMA
    double *h_uloc_pin;     // pinned mirror of uloc for H2D
    double *h_gloc_pin;     // pinned mirror of gloc for H2D
    int    *h_ok_pin;       // pinned mirror of ok for H2D
    int pin_cap;            // allocated pinned capacity

    // Profiling events
    cudaEvent_t ev_profile[N_PROFILE_EVENTS];
    int ev_initialized;
} StreamSlot;

// C API (callable from Fortran via ISO_C_BINDING)
#ifdef __cplusplus
extern "C" {
#endif
    void cuda_pool_init(int local_rank);
    int  cuda_acquire_stream(void);
    void cuda_release_stream(int slot);
    void cuda_stream_sync(int slot);
    int  cuda_stream_query(int slot);
    int  cuda_get_n_streams(void);
    void cuda_pool_finalize(void);
    int  cuda_pool_is_initialized(void);

    // Per-kernel profiling
    void hydro_cuda_profile_accumulate(int stream_slot, int ngrid);
    void hydro_cuda_profile_report(void);
    void hydro_cuda_profile_reset(void);

    // Host memory pinning
    int  cuda_host_register(void* ptr, long long nbytes);
    void cuda_host_unregister(void* ptr);
    void cuda_h2d_bandwidth_test(void);

    // GPU-resident mesh arrays for GPU-gather
    void cuda_mesh_upload(const double* uold, const double* f_grav,
                          const int* son, long long ncell, int nvar, int ndim);
    void cuda_mesh_free(void);
    int  cuda_mesh_is_ready(void);  // returns 1 if mesh allocated successfully
    cudaEvent_t cuda_get_upload_event(void);  // event signaling async upload done

    // Scatter-reduce: 5 kernels + scatter_reduce kernel, D2H compact output
    void hydro_cuda_unsplit_reduce_async(
        const double* h_uloc, const double* h_gloc, const int* h_ok,
        double* h_add_unew, double* h_add_lm1,
        double* h_add_divu_l, double* h_add_enew_l,
        double* h_add_divu_lm1, double* h_add_enew_lm1,
        double dx, double dy, double dz, double dt,
        int ngrid, int stride, int stream_slot);
    void hydro_cuda_unsplit_reduce_sync(int ngrid, int stream_slot);

    // GPU-gather + scatter-reduce: stencil idx -> GPU gather -> compute -> scatter_reduce -> compact D2H
    void hydro_cuda_gather_reduce_async(
        const int* h_stencil_idx, const int* h_stencil_grav,
        const double* h_interp_vals,
        double* h_add_unew, double* h_add_lm1,
        double* h_add_divu_l, double* h_add_enew_l,
        double* h_add_divu_lm1, double* h_add_enew_lm1,
        double dx, double dy, double dz, double dt,
        int ngrid, int stride, int n_interp, int stream_slot);
    void hydro_cuda_gather_reduce_sync(int ngrid, int stream_slot);

    // GPU-accelerated Poisson MG
    void cuda_mg_upload(const double* phi, const double* f,
                        const int* flag2, long long ncell,
                        const int* nbor_grid, const int* igrid, int ngrid);
    void cuda_mg_download_phi(double* phi, long long ncell);
    void cuda_mg_upload_phi(const double* phi, long long ncell);
    void cuda_mg_download_f1(double* f1, long long ncell);
    void cuda_mg_gauss_seidel(int ngrid, int ngridmax, int ncoarse,
                              double dx2, int color, int safe_mode);
    void cuda_mg_residual(int ngrid, int ngridmax, int ncoarse,
                          double oneoverdx2, double dtwondim, double dx2_norm,
                          double* norm2, int compute_norm);
    void cuda_mg_free(void);
    void cuda_mg_finalize(void);
    int  cuda_mg_is_ready(void);
#ifdef __cplusplus
}
#endif

// Internal helpers for .cu files
#ifdef __cplusplus
StreamSlot* get_pool();
bool is_pool_initialized();
bool pool_ensure_hydro_buffers(int slot, int ngrid);
bool pool_ensure_hydro_inter_buffers(int slot, int ngrid);
bool pool_ensure_stencil_buffers(int slot, int ngrid, int n_interp);
bool pool_ensure_reduce_buffers(int slot, int ngrid);
void pool_ensure_pinned_buffers(int slot, int ngrid);
cudaStream_t cuda_get_stream_internal(int slot);
// GPU-resident mesh array accessors
double* cuda_get_mesh_uold();
double* cuda_get_mesh_f();
int*    cuda_get_mesh_son();
long long cuda_get_mesh_ncell();
#endif

#endif // CUDA_STREAM_POOL_H
