// ==========================================================================
// Unified CUDA stream pool for cuRAMSES
// Thread-safe acquire/release for OpenMP concurrent access
// Hydro-only version (no RT/cooling)
// ==========================================================================

#include "cuda_stream_pool.h"
#include <cstdio>
#include <cstring>

// Global pool
static StreamSlot g_pool[N_STREAMS];
static bool pool_initialized = false;
static int g_device_id = 0;

// ============================================================================
// Buffer management helpers (2x over-allocation)
// ============================================================================

static void ensure_hydro_buffers(int slot, int ngrid) {
    StreamSlot& s = g_pool[slot];
    if (ngrid <= s.hydro_cap) return;
    int cap = ngrid * 2;

    size_t uloc_sz = (size_t)cap * STENCIL_NI * STENCIL_NJ * STENCIL_NK * NVAR * sizeof(double);
    size_t gloc_sz = (size_t)cap * STENCIL_NI * STENCIL_NJ * STENCIL_NK * NDIM * sizeof(double);
    size_t flux_sz = (size_t)cap * FLUX_NI * FLUX_NJ * FLUX_NK * NVAR * NDIM * sizeof(double);
    size_t tmp_sz  = (size_t)cap * FLUX_NI * FLUX_NJ * FLUX_NK * 2 * NDIM * sizeof(double);
    size_t ok_sz   = (size_t)cap * STENCIL_NI * STENCIL_NJ * STENCIL_NK * sizeof(int);

    if (s.d_uloc) cudaFree(s.d_uloc);
    if (s.d_gloc) cudaFree(s.d_gloc);
    if (s.d_flux) cudaFree(s.d_flux);
    if (s.d_tmp)  cudaFree(s.d_tmp);
    if (s.d_ok)   cudaFree(s.d_ok);

    cudaMalloc(&s.d_uloc, uloc_sz);
    cudaMalloc(&s.d_gloc, gloc_sz);
    cudaMalloc(&s.d_flux, flux_sz);
    cudaMalloc(&s.d_tmp,  tmp_sz);
    cudaMalloc(&s.d_ok,   ok_sz);
    s.hydro_cap = cap;
}

static void ensure_hydro_inter_buffers(int slot, int ngrid) {
    StreamSlot& s = g_pool[slot];
    if (ngrid <= s.hydro_inter_cap) return;
    int cap = ngrid * 2;

    size_t q_sz  = (size_t)cap * STENCIL_NI * STENCIL_NJ * STENCIL_NK * NVAR * sizeof(double);
    size_t c_sz  = (size_t)cap * STENCIL_NI * STENCIL_NJ * STENCIL_NK * sizeof(double);
    size_t dq_sz = (size_t)cap * STENCIL_NI * STENCIL_NJ * STENCIL_NK * NVAR * NDIM * sizeof(double);

    if (s.d_q)  cudaFree(s.d_q);
    if (s.d_c)  cudaFree(s.d_c);
    if (s.d_dq) cudaFree(s.d_dq);
    if (s.d_qm) cudaFree(s.d_qm);
    if (s.d_qp) cudaFree(s.d_qp);

    cudaMalloc(&s.d_q,  q_sz);
    cudaMalloc(&s.d_c,  c_sz);
    cudaMalloc(&s.d_dq, dq_sz);
    cudaMalloc(&s.d_qm, dq_sz);
    cudaMalloc(&s.d_qp, dq_sz);
    s.hydro_inter_cap = cap;
}

// ============================================================================
// Public C API
// ============================================================================

extern "C" {

void cuda_pool_init(int local_rank) {
    int device_count = 0;
    cudaGetDeviceCount(&device_count);
    if (device_count <= 0) {
        printf("CUDA pool: No GPU devices found. Running CPU-only.\n");
        pool_initialized = false;
        return;
    }
    g_device_id = local_rank % device_count;
    cudaSetDevice(g_device_id);

    cudaDeviceProp prop;
    cudaGetDeviceProperties(&prop, g_device_id);
    printf("CUDA pool: MPI local rank %d -> GPU %d (%s, %.1f GB, SM %d.%d)\n",
           local_rank, g_device_id, prop.name,
           prop.totalGlobalMem / 1073741824.0,
           prop.major, prop.minor);

    for (int s = 0; s < N_STREAMS; s++) {
        cudaStreamCreateWithFlags(&g_pool[s].stream, cudaStreamNonBlocking);
        g_pool[s].busy = 0;
        g_pool[s].d_uloc = nullptr; g_pool[s].d_gloc = nullptr;
        g_pool[s].d_flux = nullptr; g_pool[s].d_tmp  = nullptr;
        g_pool[s].d_ok   = nullptr; g_pool[s].hydro_cap = 0;
        g_pool[s].d_q  = nullptr; g_pool[s].d_c  = nullptr;
        g_pool[s].d_dq = nullptr; g_pool[s].d_qm = nullptr;
        g_pool[s].d_qp = nullptr; g_pool[s].hydro_inter_cap = 0;
    }

    pool_initialized = true;
    cudaDeviceSynchronize();
    cudaError_t err = cudaGetLastError();
    if (err != cudaSuccess) {
        fprintf(stderr, "CUDA pool init error: %s\n", cudaGetErrorString(err));
        pool_initialized = false;
    }
}

int cuda_acquire_stream(void) {
    if (!pool_initialized) return -1;
    for (int s = 0; s < N_STREAMS; s++) {
        if (!g_pool[s].busy) {
            if (__sync_lock_test_and_set(&g_pool[s].busy, 1) == 0) {
                return s;
            }
        }
        if (cudaStreamQuery(g_pool[s].stream) == cudaSuccess) {
            if (__sync_lock_test_and_set(&g_pool[s].busy, 1) == 0) {
                return s;
            }
        }
    }
    return -1;
}

void cuda_release_stream(int slot) {
    if (slot < 0 || slot >= N_STREAMS) return;
    __sync_lock_release(&g_pool[slot].busy);
}

void cuda_stream_sync(int slot) {
    if (slot < 0 || slot >= N_STREAMS) return;
    cudaStreamSynchronize(g_pool[slot].stream);
    __sync_lock_release(&g_pool[slot].busy);
}

int cuda_stream_query(int slot) {
    if (slot < 0 || slot >= N_STREAMS) return 0;
    return (cudaStreamQuery(g_pool[slot].stream) == cudaSuccess) ? 1 : 0;
}

int cuda_get_n_streams(void) {
    if (!pool_initialized) return 0;
    return N_STREAMS;
}

void cuda_pool_finalize(void) {
    if (!pool_initialized) return;
    for (int s = 0; s < N_STREAMS; s++) {
        cudaStreamSynchronize(g_pool[s].stream);
        cudaStreamDestroy(g_pool[s].stream);
        if (g_pool[s].d_uloc) cudaFree(g_pool[s].d_uloc);
        if (g_pool[s].d_gloc) cudaFree(g_pool[s].d_gloc);
        if (g_pool[s].d_flux) cudaFree(g_pool[s].d_flux);
        if (g_pool[s].d_tmp)  cudaFree(g_pool[s].d_tmp);
        if (g_pool[s].d_ok)   cudaFree(g_pool[s].d_ok);
        if (g_pool[s].d_q)    cudaFree(g_pool[s].d_q);
        if (g_pool[s].d_c)    cudaFree(g_pool[s].d_c);
        if (g_pool[s].d_dq)   cudaFree(g_pool[s].d_dq);
        if (g_pool[s].d_qm)   cudaFree(g_pool[s].d_qm);
        if (g_pool[s].d_qp)   cudaFree(g_pool[s].d_qp);
    }
    pool_initialized = false;
    printf("CUDA pool: finalized.\n");
}

int cuda_pool_is_initialized(void) {
    return pool_initialized ? 1 : 0;
}

} // extern "C"

// Internal helpers for other .cu files
StreamSlot* get_pool() { return g_pool; }
bool is_pool_initialized() { return pool_initialized; }
void pool_ensure_hydro_buffers(int slot, int ngrid) { ensure_hydro_buffers(slot, ngrid); }
void pool_ensure_hydro_inter_buffers(int slot, int ngrid) { ensure_hydro_inter_buffers(slot, ngrid); }
cudaStream_t cuda_get_stream_internal(int slot) {
    if (slot < 0 || slot >= N_STREAMS) return 0;
    return g_pool[slot].stream;
}
