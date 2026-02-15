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

#define N_STREAMS 8

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
#ifdef __cplusplus
}
#endif

// Internal helpers for .cu files
#ifdef __cplusplus
StreamSlot* get_pool();
bool is_pool_initialized();
void pool_ensure_hydro_buffers(int slot, int ngrid);
void pool_ensure_hydro_inter_buffers(int slot, int ngrid);
cudaStream_t cuda_get_stream_internal(int slot);
#endif

#endif // CUDA_STREAM_POOL_H
