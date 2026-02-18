#!/usr/bin/env python3
"""Generate strong scaling figure for cuRAMSES paper."""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# Data from Table 5 (tab:scaling)
ncpu = np.array([4, 8, 16, 24, 32, 48, 64])

# Per-rank average timers (seconds)
data = {
    'Particles':  np.array([7.16, 3.61, 1.98, 1.56, 1.15, 0.98, 0.76]),
    'Poisson':    np.array([6.11, 3.29, 2.05, 1.53, 1.38, 1.22, 1.14]),
    'MG':         np.array([18.77, 9.91, 7.01, 5.99, 5.43, 5.66, 5.69]),
    'Godunov':    np.array([3.50, 2.32, 3.38, 6.22, 8.38, 7.93, 8.64]),
    'Flag':       np.array([4.70, 2.46, 1.43, 0.97, 0.84, 0.65, 0.59]),
    'Hydro-GZ':   np.array([0.58, 0.64, 0.55, 0.52, 0.70, 0.74, 0.95]),
}
elapsed = np.array([43.4, 25.3, 20.1, 21.5, 23.5, 23.4, 25.6])
load_bal = np.array([12.2, 11.9, 16.4, 23.3, 33.1, 55.4, 83.1])

# Compute speedup relative to 4 ranks
ncpu_ref = 4
ideal = ncpu / ncpu_ref

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4.5))

# --- Left panel: Component times ---
markers = ['o', 's', '^', 'D', 'v', 'P']
colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b']
labels = ['Particles', 'Poisson', 'MG', 'Godunov', 'Flag', 'Hydro-GZ']

for i, label in enumerate(labels):
    ax1.plot(ncpu, data[label], marker=markers[i], color=colors[i],
             label=label, linewidth=1.5, markersize=5)

# Elapsed time
ax1.plot(ncpu, elapsed, 'k-', marker='*', markersize=8, linewidth=2,
         label='Elapsed', zorder=10)

# Ideal scaling reference (from 4-rank elapsed)
ideal_time = elapsed[0] * (ncpu_ref / ncpu)
ax1.plot(ncpu, ideal_time, 'k--', linewidth=1, alpha=0.4, label='Ideal')

ax1.set_xscale('log', base=2)
ax1.set_yscale('log')
ax1.set_xlabel(r'$N_{\rm cpu}$')
ax1.set_ylabel('Time (s)')
ax1.set_xticks(ncpu)
ax1.set_xticklabels([str(n) for n in ncpu])
ax1.set_ylim(0.3, 60)
ax1.legend(fontsize=7, ncol=2, loc='upper right')
ax1.set_title('(a) Per-component timers')
ax1.grid(True, alpha=0.3, which='both')

# --- Right panel: Speedup ---
for i, label in enumerate(labels):
    speedup = data[label][0] / data[label]
    ax2.plot(ncpu, speedup, marker=markers[i], color=colors[i],
             label=label, linewidth=1.5, markersize=5)

# Elapsed speedup
speedup_elapsed = elapsed[0] / elapsed
ax2.plot(ncpu, speedup_elapsed, 'k-', marker='*', markersize=8, linewidth=2,
         label='Elapsed', zorder=10)

# Ideal scaling line
ax2.plot(ncpu, ideal, 'k--', linewidth=1, alpha=0.4, label='Ideal')

ax2.set_xscale('log', base=2)
ax2.set_yscale('log', base=2)
ax2.set_xlabel(r'$N_{\rm cpu}$')
ax2.set_ylabel('Speedup relative to 4 ranks')
ax2.set_xticks(ncpu)
ax2.set_xticklabels([str(n) for n in ncpu])
ax2.set_yticks([1, 2, 4, 8, 16])
ax2.set_yticklabels(['1', '2', '4', '8', '16'])
ax2.legend(fontsize=7, ncol=2, loc='upper left')
ax2.set_title('(b) Speedup')
ax2.grid(True, alpha=0.3, which='both')

plt.tight_layout()
plt.savefig('/gpfs/kjhan/Hydro/CUBE_HR5/code_cube/misc/strong_scaling.pdf',
            bbox_inches='tight', dpi=150)
plt.savefig('/gpfs/kjhan/Hydro/CUBE_HR5/code_cube/misc/strong_scaling.png',
            bbox_inches='tight', dpi=150)
print('Saved strong_scaling.pdf and strong_scaling.png')
