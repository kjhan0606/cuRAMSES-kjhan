#!/usr/bin/env python3
"""Generate OMP scaling figure for cuRAMSES paper."""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# OMP scaling data: ncpu=4 fixed, varying OMP_NUM_THREADS
omp = np.array([1, 2, 4, 6, 8, 10, 12, 16, 20, 24, 30])
total_cores = 4 * omp

# Elapsed wall-clock time (seconds)
elapsed = np.array([2147.59, 1234.68, 736.68, 583.42, 492.53,
                    441.60, 407.30, 369.75, 374.30, 363.09, 365.59])

# Per-rank average timers (seconds)
data = {
    'MG':        np.array([1085.44, 623.45, 344.85, 258.76, 202.49,
                           169.65, 148.89, 121.24, 124.44, 114.31, 103.19]),
    'Godunov':   np.array([610.29, 306.07, 155.53, 109.07, 83.98,
                           70.14, 60.85, 53.61, 56.78, 59.02, 65.38]),
    'Sinks':     np.array([312.36, 209.03, 149.92, 130.33, 120.29,
                           114.24, 110.26, 105.91, 105.78, 102.99, 105.16]),
    'Poisson':   np.array([297.85, 193.33, 139.66, 122.71, 112.80,
                           108.85, 105.01, 100.72, 99.71, 97.50, 97.99]),
    'Particles': np.array([135.30, 90.40, 65.61, 58.06, 53.87,
                           51.54, 49.63, 47.87, 47.71, 46.38, 47.35]),
    'Feedback':  np.array([53.31, 30.24, 18.52, 15.00, 12.72,
                           11.58, 10.82, 9.73, 10.00, 9.53, 9.45]),
    'Flag':      np.array([41.34, 25.09, 15.17, 12.08, 10.49,
                           9.28, 8.60, 7.72, 7.84, 7.39, 6.99]),
}

# Compute speedup relative to OMP=1
ideal = omp.astype(float)

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4.5))

# --- Left panel: Component times ---
markers = ['o', 's', '^', 'D', 'v', 'P', 'X']
colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2']
labels = list(data.keys())

for i, label in enumerate(labels):
    ax1.plot(omp, data[label], marker=markers[i], color=colors[i],
             label=label, linewidth=1.5, markersize=5)

# Elapsed time
ax1.plot(omp, elapsed, 'k-', marker='*', markersize=8, linewidth=2,
         label='Elapsed', zorder=10)

# Ideal scaling reference (from OMP=1 elapsed)
ideal_time = elapsed[0] / omp
ax1.plot(omp, ideal_time, 'k--', linewidth=1, alpha=0.4, label='Ideal')

# Vertical line at physical core limit (64 cores / 4 MPI = 16 threads)
ax1.axvline(x=16, color='gray', linestyle=':', linewidth=1, alpha=0.5)
ax1.text(17, 2000, '64 cores', fontsize=7, color='gray', alpha=0.7)

ax1.set_xscale('log', base=2)
ax1.set_yscale('log')
ax1.set_xlabel(r'$N_{\rm omp}$ (threads per rank)')
ax1.set_ylabel('Time (s)')
ax1.set_xticks(omp)
ax1.set_xticklabels([str(n) for n in omp], fontsize=7)
ax1.set_ylim(3, 4000)
ax1.legend(fontsize=6.5, ncol=2, loc='upper right')
ax1.set_title(r'(a) Per-component timers ($N_{\rm cpu}=4$)')
ax1.grid(True, alpha=0.3, which='both')

# --- Right panel: Speedup ---
for i, label in enumerate(labels):
    speedup = data[label][0] / data[label]
    ax2.plot(omp, speedup, marker=markers[i], color=colors[i],
             label=label, linewidth=1.5, markersize=5)

# Elapsed speedup
speedup_elapsed = elapsed[0] / elapsed
ax2.plot(omp, speedup_elapsed, 'k-', marker='*', markersize=8, linewidth=2,
         label='Elapsed', zorder=10)

# Ideal scaling line
ax2.plot(omp, ideal, 'k--', linewidth=1, alpha=0.4, label='Ideal')

# Vertical line at physical core limit
ax2.axvline(x=16, color='gray', linestyle=':', linewidth=1, alpha=0.5)

ax2.set_xscale('log', base=2)
ax2.set_yscale('log', base=2)
ax2.set_xlabel(r'$N_{\rm omp}$ (threads per rank)')
ax2.set_ylabel('Speedup relative to 1 thread')
ax2.set_xticks(omp)
ax2.set_xticklabels([str(n) for n in omp], fontsize=7)
ax2.set_yticks([1, 2, 4, 8, 16, 32])
ax2.set_yticklabels(['1', '2', '4', '8', '16', '32'])
ax2.set_ylim(0.8, 40)
ax2.legend(fontsize=6.5, ncol=2, loc='upper left')
ax2.set_title('(b) Speedup')
ax2.grid(True, alpha=0.3, which='both')

plt.tight_layout()
plt.savefig('/gpfs/kjhan/Hydro/CUBE_HR5/code_cube/misc/omp_scaling.pdf',
            bbox_inches='tight', dpi=150)
plt.savefig('/gpfs/kjhan/Hydro/CUBE_HR5/code_cube/misc/omp_scaling.png',
            bbox_inches='tight', dpi=150)
print('Saved omp_scaling.pdf and omp_scaling.png')
