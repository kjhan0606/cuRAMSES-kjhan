#!/usr/bin/env python3
"""Generate strong scaling figure for cuRAMSES paper."""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# Data from Table 5 (tab:scaling) — 54M grids, 136M particles, full physics
ncpu = np.array([1, 2, 4, 8, 12, 16, 24, 32, 48, 64])

# Elapsed wall-clock time (seconds)
elapsed = np.array([8181.5, 4253.1, 2150.9, 1138.6, 793.3, 613.4, 425.3, 357.2, 271.0, 241.7])

# Per-rank average timers (seconds)
data = {
    'MG':         np.array([4034.1, 2152.5, 1084.6, 558.4, 387.3, 296.2, 205.6, 169.0, 121.9, 105.5]),
    'Godunov':    np.array([2536.4, 1218.0,  612.1, 293.9, 198.6, 144.1,  98.6,  73.9,  52.9,  43.4]),
    'Sinks':      np.array([1345.5,  656.1,  312.4, 171.9, 131.3,  94.9,  69.2,  53.7,  47.5,  41.5]),
    'Poisson':    np.array([1130.1,  584.7,  295.0, 159.1, 107.4,  85.3,  57.5,  49.4,  35.4,  31.8]),
    'Particles':  np.array([ 538.3,  272.0,  134.9,  73.9,  51.6,  38.6,  26.4,  21.5,  16.8,  13.9]),
    'Feedback':   np.array([ 210.7,  107.3,   53.3,  27.8,  20.6,  15.4,  11.3,   9.4,   7.7,   7.5]),
    'Flag':       np.array([ 161.3,   83.5,   41.5,  21.2,  14.0,  10.9,   7.8,   6.7,   5.1,   4.8]),
}
load_bal = np.array([26.4, 83.7, 50.9, 33.0, 26.0, 21.2, 18.9, 16.5, 16.1, 16.0])

# Compute speedup relative to 1 rank
ncpu_ref = 1
ideal = ncpu / ncpu_ref

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4.5))

# --- Left panel: Component times ---
markers = ['o', 's', '^', 'D', 'v', 'P', 'X']
colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2']
labels = list(data.keys())

for i, label in enumerate(labels):
    ax1.plot(ncpu, data[label], marker=markers[i], color=colors[i],
             label=label, linewidth=1.5, markersize=5)

# Elapsed time
ax1.plot(ncpu, elapsed, 'k-', marker='*', markersize=8, linewidth=2,
         label='Elapsed', zorder=10)

# Ideal scaling reference (from 1-rank elapsed)
ideal_time = elapsed[0] / ncpu
ax1.plot(ncpu, ideal_time, 'k--', linewidth=1, alpha=0.4, label='Ideal')

ax1.set_xscale('log', base=2)
ax1.set_yscale('log')
ax1.set_xlabel(r'$N_{\rm cpu}$')
ax1.set_ylabel('Time (s)')
ax1.set_xticks(ncpu)
ax1.set_xticklabels([str(n) for n in ncpu], fontsize=7)
ax1.set_ylim(3, 15000)
ax1.legend(fontsize=6.5, ncol=2, loc='upper right')
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
ax2.set_ylabel('Speedup relative to 1 rank')
ax2.set_xticks(ncpu)
ax2.set_xticklabels([str(n) for n in ncpu], fontsize=7)
ax2.set_yticks([1, 2, 4, 8, 16, 32, 64])
ax2.set_yticklabels(['1', '2', '4', '8', '16', '32', '64'])
ax2.set_ylim(0.8, 100)
ax2.legend(fontsize=6.5, ncol=2, loc='upper left')
ax2.set_title('(b) Speedup')
ax2.grid(True, alpha=0.3, which='both')

plt.tight_layout()
plt.savefig('/gpfs/kjhan/Hydro/CUBE_HR5/code_cube/misc/strong_scaling.pdf',
            bbox_inches='tight', dpi=150)
plt.savefig('/gpfs/kjhan/Hydro/CUBE_HR5/code_cube/misc/strong_scaling.png',
            bbox_inches='tight', dpi=150)
print('Saved strong_scaling.pdf and strong_scaling.png')
