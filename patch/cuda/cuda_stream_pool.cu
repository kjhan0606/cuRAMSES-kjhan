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

// Per-kernel profiling accumulators (ms -> converted to seconds for report)
// Phases: 0=H2D, 1=ctoprim, 2=uslope, 3=trace3d, 4=flux, 5=difmag, 6=D2H
static double g_profile_ms[7] = {0,0,0,0,0,0,0};
static long long g_profile_grids = 0;
static int g_profile_flushes = 0;
static const char* g_profile_names[7] = {
    "H2D transfer", "ctoprim", "uslope", "trace3d", "flux+memset", "difmag", "D2H transfer"
};

// Persistent GPU-resident mesh arrays (uploaded once per step)
static double*    d_mesh_uold = nullptr;
static double*    d_mesh_f    = nullptr;
static int*       d_mesh_son  = nullptr;
static long long  g_mesh_ncell = 0;

// Async mesh upload: dedicated stream + event for overlap with CPU work
static cudaStream_t g_upload_stream = nullptr;
static cudaEvent_t  g_upload_done_event = nullptr;
static bool         g_mesh_host_pinned = false;

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

static void ensure_stencil_buffers(int slot, int ngrid, int n_interp) {
    StreamSlot& s = g_pool[slot];
    if (ngrid > s.stencil_cap) {
        int cap = ngrid * 2;
        if (s.d_stencil_idx)  cudaFree(s.d_stencil_idx);
        if (s.d_stencil_grav) cudaFree(s.d_stencil_grav);
        size_t idx_sz = (size_t)cap * STENCIL_NI * STENCIL_NJ * STENCIL_NK * sizeof(int);
        cudaMalloc(&s.d_stencil_idx,  idx_sz);
        cudaMalloc(&s.d_stencil_grav, idx_sz);
        s.stencil_cap = cap;
    }
    if (n_interp > s.interp_cap) {
        int cap = (n_interp > 0) ? n_interp * 2 : 1024;
        if (s.d_interp_vals) cudaFree(s.d_interp_vals);
        cudaMalloc(&s.d_interp_vals, (size_t)cap * NVAR * sizeof(double));
        s.interp_cap = cap;
    }
}

static void ensure_reduce_buffers(int slot, int ngrid) {
    StreamSlot& s = g_pool[slot];
    if (ngrid <= s.reduce_cap) return;
    int cap = ngrid * 2;

    if (s.d_ok_int)       cudaFree(s.d_ok_int);
    if (s.d_add_unew)     cudaFree(s.d_add_unew);
    if (s.d_add_lm1)      cudaFree(s.d_add_lm1);
    if (s.d_add_divu_l)   cudaFree(s.d_add_divu_l);
    if (s.d_add_enew_l)   cudaFree(s.d_add_enew_l);
    if (s.d_add_divu_lm1) cudaFree(s.d_add_divu_lm1);
    if (s.d_add_enew_lm1) cudaFree(s.d_add_enew_lm1);

    size_t ok_sz = (size_t)cap * STENCIL_NI * STENCIL_NJ * STENCIL_NK * sizeof(int);
    cudaMalloc(&s.d_ok_int, ok_sz);
    cudaMalloc(&s.d_add_unew,     (size_t)cap * 8 * NVAR * sizeof(double));
    cudaMalloc(&s.d_add_lm1,      (size_t)cap * 6 * NVAR * sizeof(double));
    cudaMalloc(&s.d_add_divu_l,   (size_t)cap * 8 * sizeof(double));
    cudaMalloc(&s.d_add_enew_l,   (size_t)cap * 8 * sizeof(double));
    cudaMalloc(&s.d_add_divu_lm1, (size_t)cap * 6 * sizeof(double));
    cudaMalloc(&s.d_add_enew_lm1, (size_t)cap * 6 * sizeof(double));
    s.reduce_cap = cap;
}

static void ensure_pinned_buffers(int slot, int ngrid) {
    StreamSlot& s = g_pool[slot];
    if (ngrid <= s.pin_cap) return;
    int cap = ngrid * 2;

    if (s.h_uloc_pin) cudaFreeHost(s.h_uloc_pin);
    if (s.h_gloc_pin) cudaFreeHost(s.h_gloc_pin);
    if (s.h_ok_pin)   cudaFreeHost(s.h_ok_pin);

    size_t uloc_sz = (size_t)cap * STENCIL_NI * STENCIL_NJ * STENCIL_NK * NVAR * sizeof(double);
    size_t gloc_sz = (size_t)cap * STENCIL_NI * STENCIL_NJ * STENCIL_NK * NDIM * sizeof(double);
    size_t ok_sz   = (size_t)cap * STENCIL_NI * STENCIL_NJ * STENCIL_NK * sizeof(int);

    cudaMallocHost(&s.h_uloc_pin, uloc_sz);
    cudaMallocHost(&s.h_gloc_pin, gloc_sz);
    cudaMallocHost(&s.h_ok_pin,   ok_sz);
    s.pin_cap = cap;
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
        g_pool[s].d_stencil_idx = nullptr; g_pool[s].d_stencil_grav = nullptr;
        g_pool[s].d_interp_vals = nullptr;
        g_pool[s].stencil_cap = 0; g_pool[s].interp_cap = 0;
        g_pool[s].d_ok_int = nullptr;
        g_pool[s].d_add_unew = nullptr; g_pool[s].d_add_lm1 = nullptr;
        g_pool[s].d_add_divu_l = nullptr; g_pool[s].d_add_enew_l = nullptr;
        g_pool[s].d_add_divu_lm1 = nullptr; g_pool[s].d_add_enew_lm1 = nullptr;
        g_pool[s].reduce_cap = 0;
        g_pool[s].h_uloc_pin = nullptr; g_pool[s].h_gloc_pin = nullptr;
        g_pool[s].h_ok_pin = nullptr; g_pool[s].pin_cap = 0;
        // Initialize profiling events
        for (int e = 0; e < N_PROFILE_EVENTS; e++) {
            cudaEventCreate(&g_pool[s].ev_profile[e]);
        }
        g_pool[s].ev_initialized = 1;
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
    // Ensure this OMP thread uses the correct GPU device (per-thread context)
    cudaSetDevice(g_device_id);
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
        if (g_pool[s].d_stencil_idx)  cudaFree(g_pool[s].d_stencil_idx);
        if (g_pool[s].d_stencil_grav) cudaFree(g_pool[s].d_stencil_grav);
        if (g_pool[s].d_interp_vals)  cudaFree(g_pool[s].d_interp_vals);
        if (g_pool[s].d_ok_int)       cudaFree(g_pool[s].d_ok_int);
        if (g_pool[s].d_add_unew)     cudaFree(g_pool[s].d_add_unew);
        if (g_pool[s].d_add_lm1)      cudaFree(g_pool[s].d_add_lm1);
        if (g_pool[s].d_add_divu_l)   cudaFree(g_pool[s].d_add_divu_l);
        if (g_pool[s].d_add_enew_l)   cudaFree(g_pool[s].d_add_enew_l);
        if (g_pool[s].d_add_divu_lm1) cudaFree(g_pool[s].d_add_divu_lm1);
        if (g_pool[s].d_add_enew_lm1) cudaFree(g_pool[s].d_add_enew_lm1);
        if (g_pool[s].h_uloc_pin) cudaFreeHost(g_pool[s].h_uloc_pin);
        if (g_pool[s].h_gloc_pin) cudaFreeHost(g_pool[s].h_gloc_pin);
        if (g_pool[s].h_ok_pin)   cudaFreeHost(g_pool[s].h_ok_pin);
        if (g_pool[s].ev_initialized) {
            for (int e = 0; e < N_PROFILE_EVENTS; e++) {
                cudaEventDestroy(g_pool[s].ev_profile[e]);
            }
        }
    }
    // Free persistent mesh arrays
    if (d_mesh_uold) { cudaFree(d_mesh_uold); d_mesh_uold = nullptr; }
    if (d_mesh_f)    { cudaFree(d_mesh_f);    d_mesh_f    = nullptr; }
    if (d_mesh_son)  { cudaFree(d_mesh_son);  d_mesh_son  = nullptr; }
    g_mesh_ncell = 0;
    pool_initialized = false;
    printf("CUDA pool: finalized.\n");
}

int cuda_pool_is_initialized(void) {
    return pool_initialized ? 1 : 0;
}

void cuda_mesh_upload(const double* uold, const double* f_grav,
                      const int* son, long long ncell, int nvar, int ndim) {
    if (!pool_initialized) return;

    // Create dedicated upload stream + event on first call
    if (!g_upload_stream) {
        cudaStreamCreateWithFlags(&g_upload_stream, cudaStreamNonBlocking);
        cudaEventCreateWithFlags(&g_upload_done_event, cudaEventDisableTiming);
    }

    // Reallocate if size changed
    if (ncell != g_mesh_ncell) {
        // Unpin previous host arrays if pinned
        if (g_mesh_host_pinned) {
            cudaHostUnregister((void*)uold);  // same pointer, just unpin
            if (f_grav) cudaHostUnregister((void*)f_grav);
            cudaHostUnregister((void*)son);
            g_mesh_host_pinned = false;
        }
        if (d_mesh_uold) { cudaFree(d_mesh_uold); d_mesh_uold = nullptr; }
        if (d_mesh_f)    { cudaFree(d_mesh_f);    d_mesh_f    = nullptr; }
        if (d_mesh_son)  { cudaFree(d_mesh_son);  d_mesh_son  = nullptr; }

        double gb = ((double)ncell * (nvar + ndim) * sizeof(double) +
                     (double)ncell * sizeof(int)) / (1024.0*1024.0*1024.0);

        // Check available GPU memory before allocating
        size_t free_mem = 0, total_mem = 0;
        cudaMemGetInfo(&free_mem, &total_mem);
        size_t need = (size_t)ncell * nvar * sizeof(double)
                    + (size_t)ncell * ndim * sizeof(double)
                    + (size_t)ncell * sizeof(int);
        if (need > free_mem) {
            fprintf(stderr, "CUDA mesh: SKIP — need %.1f GB but only %.1f GB free (total %.1f GB)\n",
                    (double)need / (1024.0*1024.0*1024.0),
                    (double)free_mem / (1024.0*1024.0*1024.0),
                    (double)total_mem / (1024.0*1024.0*1024.0));
            g_mesh_ncell = 0;
            return;
        }

        cudaError_t e1 = cudaMalloc(&d_mesh_uold, (size_t)ncell * nvar * sizeof(double));
        cudaError_t e2 = cudaMalloc(&d_mesh_f,    (size_t)ncell * ndim * sizeof(double));
        cudaError_t e3 = cudaMalloc(&d_mesh_son,  (size_t)ncell * sizeof(int));
        if (e1 != cudaSuccess || e2 != cudaSuccess || e3 != cudaSuccess) {
            fprintf(stderr, "CUDA mesh: allocation FAILED (%.1f GB). Falling back to CPU gather.\n", gb);
            if (d_mesh_uold) { cudaFree(d_mesh_uold); d_mesh_uold = nullptr; }
            if (d_mesh_f)    { cudaFree(d_mesh_f);    d_mesh_f    = nullptr; }
            if (d_mesh_son)  { cudaFree(d_mesh_son);  d_mesh_son  = nullptr; }
            g_mesh_ncell = 0;
            return;
        }
        g_mesh_ncell = ncell;
        printf("CUDA mesh: allocated %.1f GB (ncell=%lld, nvar=%d, free=%.1f/%.1f GB)\n",
               gb, ncell, nvar,
               (double)(free_mem - need) / (1024.0*1024.0*1024.0),
               (double)total_mem / (1024.0*1024.0*1024.0));
    }

    // Pin host arrays for fast async DMA (one-time cost per allocation)
    if (!g_mesh_host_pinned && g_mesh_ncell > 0) {
        cudaError_t ep = cudaHostRegister((void*)uold,
            (size_t)ncell * nvar * sizeof(double), cudaHostRegisterDefault);
        if (ep == cudaSuccess) {
            if (f_grav)
                cudaHostRegister((void*)f_grav,
                    (size_t)ncell * ndim * sizeof(double), cudaHostRegisterDefault);
            cudaHostRegister((void*)son,
                (size_t)ncell * sizeof(int), cudaHostRegisterDefault);
            g_mesh_host_pinned = true;
            double gb_pin = ((double)ncell * (nvar + ndim) * sizeof(double) +
                             (double)ncell * sizeof(int)) / (1024.0*1024.0*1024.0);
            printf("CUDA mesh: pinned %.1f GB host memory for async DMA\n", gb_pin);
        }
    }

    // Async upload on dedicated stream (non-blocking to CPU)
    if (g_mesh_host_pinned) {
        cudaMemcpyAsync(d_mesh_uold, uold,
            (size_t)ncell * nvar * sizeof(double), cudaMemcpyHostToDevice, g_upload_stream);
        if (f_grav) {
            cudaMemcpyAsync(d_mesh_f, f_grav,
                (size_t)ncell * ndim * sizeof(double), cudaMemcpyHostToDevice, g_upload_stream);
        } else {
            cudaMemsetAsync(d_mesh_f, 0,
                (size_t)ncell * ndim * sizeof(double), g_upload_stream);
        }
        cudaMemcpyAsync(d_mesh_son, son,
            (size_t)ncell * sizeof(int), cudaMemcpyHostToDevice, g_upload_stream);
        cudaEventRecord(g_upload_done_event, g_upload_stream);
    } else {
        // Sync fallback if pinning failed
        cudaMemcpy(d_mesh_uold, uold,
            (size_t)ncell * nvar * sizeof(double), cudaMemcpyHostToDevice);
        if (f_grav) {
            cudaMemcpy(d_mesh_f, f_grav,
                (size_t)ncell * ndim * sizeof(double), cudaMemcpyHostToDevice);
        } else {
            cudaMemset(d_mesh_f, 0, (size_t)ncell * ndim * sizeof(double));
        }
        cudaMemcpy(d_mesh_son, son,
            (size_t)ncell * sizeof(int), cudaMemcpyHostToDevice);
    }
}

void cuda_mesh_free(void) {
    if (g_upload_stream) {
        cudaStreamSynchronize(g_upload_stream);
        cudaStreamDestroy(g_upload_stream);  g_upload_stream = nullptr;
    }
    if (g_upload_done_event) {
        cudaEventDestroy(g_upload_done_event); g_upload_done_event = nullptr;
    }
    // Note: host arrays are still in use by Fortran, don't unregister here
    // (they will be freed by Fortran deallocate which handles unpin)
    g_mesh_host_pinned = false;
    if (d_mesh_uold) { cudaFree(d_mesh_uold); d_mesh_uold = nullptr; }
    if (d_mesh_f)    { cudaFree(d_mesh_f);    d_mesh_f    = nullptr; }
    if (d_mesh_son)  { cudaFree(d_mesh_son);  d_mesh_son  = nullptr; }
    g_mesh_ncell = 0;
}

void hydro_cuda_profile_accumulate(int stream_slot, int ngrid) {
    if (stream_slot < 0 || stream_slot >= N_STREAMS) return;
    StreamSlot& s = g_pool[stream_slot];
    if (!s.ev_initialized) return;
    // Events must already be recorded and stream synchronized
    for (int p = 0; p < 7; p++) {
        float ms = 0;
        cudaEventElapsedTime(&ms, s.ev_profile[p], s.ev_profile[p+1]);
        g_profile_ms[p] += (double)ms;
    }
    g_profile_grids += ngrid;
    g_profile_flushes++;
}

void hydro_cuda_profile_report(void) {
    double total = 0;
    for (int p = 0; p < 7; p++) total += g_profile_ms[p];
    if (total <= 0 || g_profile_flushes == 0) return;

    printf(" === GPU Hydro Kernel Profiling ===\n");
    printf("   Flushes: %d, Total grids: %lld\n", g_profile_flushes, g_profile_grids);
    printf("   %-14s %10s %6s\n", "Phase", "Time(s)", "%");
    for (int p = 0; p < 7; p++) {
        printf("   %-14s %10.3f %5.1f%%\n",
               g_profile_names[p], g_profile_ms[p] / 1000.0,
               g_profile_ms[p] / total * 100.0);
    }
    printf("   %-14s %10.3f\n", "Total", total / 1000.0);
}

void hydro_cuda_profile_reset(void) {
    for (int p = 0; p < 7; p++) g_profile_ms[p] = 0;
    g_profile_grids = 0;
    g_profile_flushes = 0;
}

int cuda_host_register(void* ptr, long long nbytes) {
    if (!pool_initialized || !ptr || nbytes <= 0) return -1;
    cudaError_t err = cudaHostRegister(ptr, (size_t)nbytes, cudaHostRegisterDefault);
    if (err != cudaSuccess) {
        fprintf(stderr, "cudaHostRegister(%p, %lld bytes) failed: %s\n",
                ptr, nbytes, cudaGetErrorString(err));
        return -1;
    }
    return 0;
}

void cuda_host_unregister(void* ptr) {
    if (!pool_initialized || !ptr) return;
    cudaHostUnregister(ptr);
}

void cuda_h2d_bandwidth_test(void) {
    if (!pool_initialized) return;
    const size_t test_sz = 100 * 1024 * 1024; // 100 MB
    double *h_pin = nullptr, *d_buf = nullptr;
    double *h_page = (double*)malloc(test_sz);
    cudaMallocHost(&h_pin, test_sz);
    cudaMalloc(&d_buf, test_sz);
    if (!h_pin || !d_buf || !h_page) {
        printf("CUDA BW test: alloc failed\n");
        if (h_pin) cudaFreeHost(h_pin); if (d_buf) cudaFree(d_buf);
        free(h_page);
        return;
    }
    memset(h_pin, 0, test_sz);
    memset(h_page, 0, test_sz);
    cudaDeviceSynchronize();

    // Test 1: pinned H2D
    cudaEvent_t e0, e1;
    cudaEventCreate(&e0); cudaEventCreate(&e1);
    cudaEventRecord(e0);
    for (int i = 0; i < 10; i++)
        cudaMemcpy(d_buf, h_pin, test_sz, cudaMemcpyHostToDevice);
    cudaEventRecord(e1);
    cudaEventSynchronize(e1);
    float ms1 = 0;
    cudaEventElapsedTime(&ms1, e0, e1);
    double bw_pin = 10.0 * test_sz / (ms1 / 1000.0) / 1e9;

    // Test 2: pageable H2D
    cudaEventRecord(e0);
    for (int i = 0; i < 10; i++)
        cudaMemcpy(d_buf, h_page, test_sz, cudaMemcpyHostToDevice);
    cudaEventRecord(e1);
    cudaEventSynchronize(e1);
    float ms2 = 0;
    cudaEventElapsedTime(&ms2, e0, e1);
    double bw_page = 10.0 * test_sz / (ms2 / 1000.0) / 1e9;

    // Test 3: registered H2D
    cudaHostRegister(h_page, test_sz, cudaHostRegisterDefault);
    cudaEventRecord(e0);
    for (int i = 0; i < 10; i++)
        cudaMemcpy(d_buf, h_page, test_sz, cudaMemcpyHostToDevice);
    cudaEventRecord(e1);
    cudaEventSynchronize(e1);
    float ms3 = 0;
    cudaEventElapsedTime(&ms3, e0, e1);
    double bw_reg = 10.0 * test_sz / (ms3 / 1000.0) / 1e9;
    cudaHostUnregister(h_page);

    printf(" CUDA H2D bandwidth test (100 MB x 10):\n");
    printf("   Pinned:     %.1f ms -> %.1f GB/s\n", ms1, bw_pin);
    printf("   Pageable:   %.1f ms -> %.1f GB/s\n", ms2, bw_page);
    printf("   Registered: %.1f ms -> %.1f GB/s\n", ms3, bw_reg);

    cudaEventDestroy(e0); cudaEventDestroy(e1);
    cudaFreeHost(h_pin); cudaFree(d_buf); free(h_page);
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
void pool_ensure_stencil_buffers(int slot, int ngrid, int n_interp) {
    ensure_stencil_buffers(slot, ngrid, n_interp);
}
void pool_ensure_reduce_buffers(int slot, int ngrid) {
    ensure_reduce_buffers(slot, ngrid);
}
void pool_ensure_pinned_buffers(int slot, int ngrid) {
    ensure_pinned_buffers(slot, ngrid);
}
double*   cuda_get_mesh_uold()  { return d_mesh_uold; }
double*   cuda_get_mesh_f()     { return d_mesh_f; }
int*      cuda_get_mesh_son()   { return d_mesh_son; }
long long cuda_get_mesh_ncell() { return g_mesh_ncell; }
int       cuda_mesh_is_ready()  { return (g_mesh_ncell > 0 && d_mesh_uold && d_mesh_son) ? 1 : 0; }
cudaEvent_t cuda_get_upload_event() { return g_upload_done_event; }
