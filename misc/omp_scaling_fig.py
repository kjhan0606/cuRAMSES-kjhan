#!/usr/bin/env python3
"""Generate hybrid MPI/OpenMP scaling figure for cuRAMSES paper (Cosmo1024).

Two test series on Grammar cluster (8 nodes, 64 cores/node, 512 cores total):
  - Hybrid scaling: ncpu * OMP = 512, OMP = 1,2,4,8,16,32,64
  - Pure OMP scaling: ncpu = 64 fixed, OMP = 1,2,4 (+ OMP=8 from hybrid)

The figure uses the hybrid series (constant total resources) as it matches
the paper text (Section 8.5.2).
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# === Hybrid scaling data: ncpu * OMP = 512 (8 nodes) ===
omp_hybrid = np.array([1,    2,    4,    8,   16,   32,   64])
ncpu_hybrid = 512 // omp_hybrid  # [512, 256, 128, 64, 32, 16, 8]

# Per-step wall-clock (average µs/pt from logs, excluding first step)
uspt_hybrid = np.array([70.59, 24.85, 11.61, 7.16, 4.92, 4.17, 4.60])

# Component average times (seconds) from final TIMER block
timer_hybrid = {
    'MG AMR':    np.array([564.26, 613.37, 649.75, 798.67, 1012.48, 1512.65, 3434.85]),
    'MG base':   np.array([27.89,  38.53,  52.90,  90.57,  157.13,  326.81,  676.45]),
    'Godunov':   np.array([154.67, 156.31, 161.80, 169.72, 175.49,  268.86,  728.85]),
    'Sinks':     np.array([626.84, 373.59, 278.30, 354.09, 538.92,  956.31, 2250.60]),
    'Poisson':   np.array([74.93,  94.90, 130.82, 206.09, 347.52,  657.63, 1316.91]),
    'Particles': np.array([60.02,  62.56,  73.34, 113.62, 186.64,  364.14,  795.47]),
    'Feedback':  np.array([268.43, 97.16,  51.75,  44.09,  48.16,   78.50,  125.04]),
    'Flag':      np.array([7.99,    8.97,  10.37,  14.26,  21.15,   34.74,   70.79]),
}
total_hybrid = np.array([2752.10, 2119.10, 2067.79, 1905.52, 3121.54, 4989.99, 10583.00])

# Elapsed per coarse step (seconds) = total / nsteps (5 steps each)
elapsed_hybrid = total_hybrid / 5.0

# === Figure ===
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4.5))

# Color scheme
comp_order = ['MG AMR', 'MG base', 'Godunov', 'Sinks', 'Poisson', 'Particles', 'Feedback', 'Flag']
colors = ['#d62728', '#ff9896', '#1f77b4', '#2ca02c', '#ff7f0e', '#9467bd', '#8c564b', '#e377c2']
markers = ['o', 'o', 's', '^', 'D', 'v', 'P', 'X']
fills   = ['full','none','full','full','full','full','full','full']

# --- Left panel: stacked bar of component times ---
bar_width = 0.7
x_pos = np.arange(len(omp_hybrid))
bottom = np.zeros(len(omp_hybrid))

for i, comp in enumerate(comp_order):
    vals = timer_hybrid[comp] / 5.0  # per-step
    ax1.bar(x_pos, vals, bar_width, bottom=bottom, color=colors[i],
            label=comp, edgecolor='white', linewidth=0.3)
    bottom += vals

ax1.set_xticks(x_pos)
ax1.set_xticklabels([f'{n}\n({c})' for n, c in zip(omp_hybrid, ncpu_hybrid)], fontsize=7.5)
ax1.set_xlabel(r'$N_{\rm thread}$ ($N_{\rm rank}$)')
ax1.set_ylabel('Wall-clock time per step (s)')
_bbox = dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.85, edgecolor='none')
ax1.text(0.03, 0.97, r'(a) Per-step time ($N_{\rm rank} \times N_{\rm thread} = 512$)',
         transform=ax1.transAxes, fontsize=8, fontweight='bold',
         va='top', ha='left', bbox=_bbox, zorder=20)
ax1.legend(fontsize=6.5, ncol=2, loc='upper left',
           bbox_to_anchor=(0.0, 0.87))
ax1.set_ylim(0, bottom.max() * 1.15)
ax1.grid(True, alpha=0.2, axis='y')

# Highlight production config (OMP=8)
prod_idx = np.where(omp_hybrid == 8)[0][0]
ax1.axvline(x=prod_idx, color='red', linestyle='--', linewidth=1.2, alpha=0.6)
ax1.annotate('Production', xy=(prod_idx, bottom[prod_idx]*1.02),
             fontsize=7, color='red', ha='center', va='bottom')

# --- Right panel: speedup relative to OMP=1 (pure MPI) ---
# Use µs/pt for speedup (lower is better → invert)
speedup = uspt_hybrid[0] / uspt_hybrid
ideal = omp_hybrid.astype(float)

ax2.plot(omp_hybrid, speedup, 'ko-', linewidth=2, markersize=7, label='Measured', zorder=10)
ax2.plot(omp_hybrid, ideal, 'k--', linewidth=1, alpha=0.4, label='Ideal')

# Fill between measured and ideal
ax2.fill_between(omp_hybrid, speedup, ideal, alpha=0.1, color='blue')

# Annotate efficiency
for i, (o, s) in enumerate(zip(omp_hybrid, speedup)):
    eff = s / o * 100
    if o > 1:
        ax2.annotate(f'{eff:.0f}%', xy=(o, s), xytext=(5, -12),
                     textcoords='offset points', fontsize=7, color='#333')

# Physical core limit: 64 cores/node, 8 ranks/node → 8 threads = physical
ax2.axvline(x=8, color='gray', linestyle=':', linewidth=1, alpha=0.6)
ax2.text(9, 1.5, 'Physical\ncores', fontsize=6.5, color='gray', alpha=0.8)

# Production config
ax2.axvline(x=8, color='red', linestyle='--', linewidth=1.2, alpha=0.3)

ax2.set_xscale('log', base=2)
ax2.set_yscale('log', base=2)
ax2.set_xlabel(r'$N_{\rm thread}$ (threads per rank)')
ax2.set_ylabel(r'Speedup ($\mu$s/pt relative to $N_{\rm thread}=1$)')
ax2.set_xticks(omp_hybrid)
ax2.set_xticklabels([str(n) for n in omp_hybrid])
ax2.set_yticks([1, 2, 4, 8, 16, 32, 64])
ax2.set_yticklabels(['1', '2', '4', '8', '16', '32', '64'])
ax2.set_ylim(0.8, 80)
ax2.set_xlim(0.7, 90)
ax2.legend(fontsize=8, loc='upper left',
           bbox_to_anchor=(0.0, 0.87))
ax2.text(0.03, 0.97, r'(b) Speedup vs $N_{\rm thread}=1$',
         transform=ax2.transAxes, fontsize=8, fontweight='bold',
         va='top', ha='left', bbox=_bbox, zorder=20)
ax2.grid(True, alpha=0.2, which='both')

# Upper x-axis: ncpu
ax2t = ax2.twiny()
ax2t.set_xscale('log', base=2)
ax2t.set_xlim(ax2.get_xlim())
ax2t.set_xticks(omp_hybrid)
ax2t.set_xticklabels([str(n) for n in ncpu_hybrid], fontsize=7)
ax2t.set_xlabel(r'$N_{\rm rank}$ (MPI ranks)', fontsize=8)

plt.tight_layout()
outbase = '/gpfs/kjhan/Hydro/CUBE_HR5/code_cube/misc/fig_cosmo1024_omp_scaling'
plt.savefig(outbase + '.pdf', bbox_inches='tight', dpi=150)
plt.savefig(outbase + '.png', bbox_inches='tight', dpi=150)
print(f'Saved {outbase}.pdf and .png')

# Print summary table
print('\n=== Hybrid Scaling Summary (ncpu × OMP = 512, 8 nodes) ===')
print(f'{"OMP":>4} {"ncpu":>5} {"µs/pt":>7} {"Speedup":>8} {"Eff%":>6}')
for o, n, u in zip(omp_hybrid, ncpu_hybrid, uspt_hybrid):
    sp = uspt_hybrid[0] / u
    eff = sp / o * 100
    print(f'{o:4d} {n:5d} {u:7.2f} {sp:8.2f}x {eff:5.1f}%')

# =====================================================================
# Cosmo512 OMP scaling: ncpu=4 fixed, OMP = 1..30 (single node)
# =====================================================================
omp_512 = np.array([1, 2, 4, 6, 8, 10, 12, 16, 20, 24, 30])
elapsed_512 = np.array([2147.59, 1234.68, 736.68, 583.42, 492.53,
                        441.60, 407.30, 369.75, 374.30, 363.09, 365.59])
data_512 = {
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

fig2, (bx1, bx2) = plt.subplots(1, 2, figsize=(10, 4.5))

markers_512 = ['o', 's', '^', 'D', 'v', 'P', 'X']
colors_512 = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2']

for i, label in enumerate(data_512):
    bx1.plot(omp_512, data_512[label], marker=markers_512[i], color=colors_512[i],
             label=label, linewidth=1.5, markersize=5)
bx1.plot(omp_512, elapsed_512, 'k-', marker='*', markersize=8, linewidth=2,
         label='Elapsed', zorder=10)
ideal_time_512 = elapsed_512[0] / omp_512
bx1.plot(omp_512, ideal_time_512, 'k--', linewidth=1, alpha=0.4, label='Ideal')
bx1.axvline(x=16, color='gray', linestyle=':', linewidth=1, alpha=0.5)
bx1.text(17, 2000, '64 cores', fontsize=7, color='gray', alpha=0.7)
bx1.set_xscale('log', base=2)
bx1.set_yscale('log')
bx1.set_xlabel(r'$N_{\rm thread}$ (threads per rank)')
bx1.set_ylabel('Time (s)')
bx1.set_xticks(omp_512)
bx1.set_xticklabels([str(n) for n in omp_512], fontsize=7)
bx1.set_ylim(3, 4000)
bx1.legend(fontsize=6.5, ncol=2, loc='upper right')
bx1.text(0.03, 0.97, r'(a) Per-component timers ($N_{\rm rank}=4$)',
         transform=bx1.transAxes, fontsize=8, fontweight='bold',
         va='top', ha='left', bbox=_bbox, zorder=20)
bx1.grid(True, alpha=0.3, which='both')

for i, label in enumerate(data_512):
    speedup_512 = data_512[label][0] / data_512[label]
    bx2.plot(omp_512, speedup_512, marker=markers_512[i], color=colors_512[i],
             label=label, linewidth=1.5, markersize=5)
speedup_elapsed_512 = elapsed_512[0] / elapsed_512
bx2.plot(omp_512, speedup_elapsed_512, 'k-', marker='*', markersize=8, linewidth=2,
         label='Elapsed', zorder=10)
bx2.plot(omp_512, omp_512.astype(float), 'k--', linewidth=1, alpha=0.4, label='Ideal')
bx2.axvline(x=16, color='gray', linestyle=':', linewidth=1, alpha=0.5)
bx2.set_xscale('log', base=2)
bx2.set_yscale('log', base=2)
bx2.set_xlabel(r'$N_{\rm thread}$ (threads per rank)')
bx2.set_ylabel('Speedup relative to 1 thread')
bx2.set_xticks(omp_512)
bx2.set_xticklabels([str(n) for n in omp_512], fontsize=7)
bx2.set_yticks([1, 2, 4, 8, 16, 32])
bx2.set_yticklabels(['1', '2', '4', '8', '16', '32'])
bx2.set_ylim(0.8, 40)
bx2.legend(fontsize=6.5, ncol=2, loc='upper left',
           bbox_to_anchor=(0.0, 0.87))
bx2.text(0.03, 0.97, '(b) Speedup',
         transform=bx2.transAxes, fontsize=8, fontweight='bold',
         va='top', ha='left', bbox=_bbox, zorder=20)
bx2.grid(True, alpha=0.3, which='both')

plt.tight_layout()
outbase2 = '/gpfs/kjhan/Hydro/CUBE_HR5/code_cube/misc/omp_scaling'
plt.savefig(outbase2 + '.pdf', bbox_inches='tight', dpi=150)
plt.savefig(outbase2 + '.png', bbox_inches='tight', dpi=150)
print(f'Saved {outbase2}.pdf and .png')
