#!/usr/bin/env python3
"""Generate combined strong scaling figure (Cosmo512 + Cosmo1024) for cuRAMSES paper.

Figure layout: 2x2 panels
  (a) Cosmo512 per-component timers   (b) Cosmo512 speedup
  (c) Cosmo1024 per-component timers  (d) Cosmo1024 speedup
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# =========================================================================
# Cosmo512: single-node, ncpu = 1..64, OMP = 1 (pure MPI)
# =========================================================================
ncpu_512 = np.array([1, 2, 4, 8, 12, 16, 24, 32, 48, 64])
omp_512 = 1
elapsed_512 = np.array([8181.5, 4253.1, 2150.9, 1138.6, 793.3, 613.4, 425.3, 357.2, 271.0, 241.7])
data_512 = {
    'MG AMR':     np.array([4034.1, 2152.5, 1084.6, 558.4, 387.3, 296.2, 205.6, 169.0, 121.9, 105.5]),
    'Godunov':    np.array([2536.4, 1218.0,  612.1, 293.9, 198.6, 144.1,  98.6,  73.9,  52.9,  43.4]),
    'Sinks':      np.array([1345.5,  656.1,  312.4, 171.9, 131.3,  94.9,  69.2,  53.7,  47.5,  41.5]),
    'Poisson':    np.array([1130.1,  584.7,  295.0, 159.1, 107.4,  85.3,  57.5,  49.4,  35.4,  31.8]),
    'Particles':  np.array([ 538.3,  272.0,  134.9,  73.9,  51.6,  38.6,  26.4,  21.5,  16.8,  13.9]),
    'Feedback':   np.array([ 210.7,  107.3,   53.3,  27.8,  20.6,  15.4,  11.3,   9.4,   7.7,   7.5]),
    'Flag':       np.array([ 161.3,   83.5,   41.5,  21.2,  14.0,  10.9,   7.8,   6.7,   5.1,   4.8]),
}

# =========================================================================
# Cosmo1024: multi-node Grammar, ncpu = 8..256, OMP = 8
# =========================================================================
ncpu_1024 = np.array([8, 16, 32, 64, 128, 256])
omp_1024 = 8
nodes_1024 = ncpu_1024 // 8
cores_1024 = ncpu_1024 * omp_1024
# Per-step elapsed (average of 4 measured coarse steps, excluding warmup)
elapsed_1024 = np.array([2858.0, 1432.0, 721.0, 388.0, 182.0, 108.0])
# TOTAL timer from finalize_timer (accumulated over all 5 steps + startup)
total_1024 = np.array([18048.5, 9236.9, 4809.0, 2218.0, 1446.0, 885.1])
# Raw component data (accumulated totals from finalize_timer, average across ranks)
data_1024_raw = {
    'MG AMR':     np.array([8218.07, 4079.30, 1972.73,  965.11,  378.06,  185.26]),
    'Godunov':    np.array([1471.26,  675.78,  336.89,  170.27,   81.50,   39.85]),
    'Sinks':      np.array([2609.54, 1349.74,  737.36,  457.42,  205.43,  143.56]),
    'Poisson':    np.array([1516.98,  780.12,  391.40,  210.31,  108.33,   56.52]),
    'MG base':    np.array([ 635.44,  345.35,  166.82,   91.24,   46.51,   25.34]),
    'Particles':  np.array([ 880.24,  452.57,  228.56,  126.74,   61.36,   38.58]),
    'Feedback':   np.array([ 240.34,  122.62,   69.95,   45.49,   31.14,   29.53]),
    'Flag':       np.array([ 133.50,   70.03,   34.91,   17.29,    7.34,    3.83]),
}
# Convert to per-step by dividing uniformly by nsteps=5 (all coarse steps).
# A uniform divisor preserves the raw speedup ratios across rank counts;
# a per-rank scaling factor (elapsed/total) would distort them because the
# warmup-step overhead (load_balance) varies with rank count.
nsteps_1024 = 5
data_1024 = {comp: vals / nsteps_1024 for comp, vals in data_1024_raw.items()}

# =========================================================================
# Shared style
# =========================================================================
comp_order = ['MG AMR', 'Godunov', 'Sinks', 'Poisson', 'Particles', 'Feedback', 'Flag']
colors =     ['#d62728', '#1f77b4', '#2ca02c', '#ff7f0e', '#9467bd', '#8c564b', '#e377c2']
markers =    ['o',       's',       '^',       'D',       'v',       'P',       'X']

fig, axes = plt.subplots(2, 2, figsize=(10, 8.5))

def plot_panel(ax_time, ax_spd, ncpu, elapsed, data, omp, title_prefix, xlabel):
    """Plot (left) component timer lines + elapsed, (right) speedup."""
    # --- Component timers ---
    for i, comp in enumerate(comp_order):
        if comp not in data:
            continue
        ax_time.plot(ncpu, data[comp], marker=markers[i], color=colors[i],
                     label=comp, linewidth=1.5, markersize=5)

    # MG base (only Cosmo1024)
    if 'MG base' in data:
        ax_time.plot(ncpu, data['MG base'], marker='o', color='#ff9896',
                     fillstyle='none', label='MG base', linewidth=1.2, markersize=5)

    # Elapsed
    ax_time.plot(ncpu, elapsed, 'k-', marker='*', markersize=8, linewidth=2,
                 label='Elapsed', zorder=10)
    # Ideal
    ideal_time = elapsed[0] / (ncpu / ncpu[0])
    ax_time.plot(ncpu, ideal_time, 'k--', linewidth=1, alpha=0.4, label='Ideal')

    ax_time.set_xscale('log', base=2)
    ax_time.set_yscale('log')
    ax_time.set_xlabel(xlabel)
    ax_time.set_ylabel('Time (s)')
    ax_time.set_xticks(ncpu)
    ax_time.set_xticklabels([str(n) for n in ncpu], fontsize=7)
    ax_time.legend(fontsize=6, ncol=2, loc='upper right')
    ax_time.grid(True, alpha=0.3, which='both')

    # --- Speedup ---
    for i, comp in enumerate(comp_order):
        if comp not in data:
            continue
        spd = data[comp][0] / data[comp]
        ax_spd.plot(ncpu, spd, marker=markers[i], color=colors[i],
                    label=comp, linewidth=1.5, markersize=5)
    if 'MG base' in data:
        spd_mgb = data['MG base'][0] / data['MG base']
        ax_spd.plot(ncpu, spd_mgb, marker='o', color='#ff9896',
                    fillstyle='none', label='MG base', linewidth=1.2, markersize=5)

    spd_elapsed = elapsed[0] / elapsed
    ax_spd.plot(ncpu, spd_elapsed, 'k-', marker='*', markersize=8, linewidth=2,
                label='Elapsed', zorder=10)
    ideal_spd = ncpu / ncpu[0]
    ax_spd.plot(ncpu, ideal_spd, 'k--', linewidth=1, alpha=0.4, label='Ideal')

    ax_spd.set_xscale('log', base=2)
    ax_spd.set_yscale('log', base=2)
    ax_spd.set_xlabel(xlabel)
    ax_spd.set_ylabel('Speedup')
    ax_spd.set_xticks(ncpu)
    ax_spd.set_xticklabels([str(n) for n in ncpu], fontsize=7)

    max_ratio = int(ncpu[-1] / ncpu[0])
    yticks = [2**i for i in range(0, int(np.log2(max_ratio))+2)]
    ax_spd.set_yticks(yticks)
    ax_spd.set_yticklabels([str(y) for y in yticks])
    ax_spd.set_ylim(0.8, max_ratio * 2)
    ax_spd.legend(fontsize=6, ncol=2, loc='upper left',
                   bbox_to_anchor=(0.0, 0.87))
    ax_spd.grid(True, alpha=0.3, which='both')

# --- Cosmo512 (top row) ---
plot_panel(axes[0,0], axes[0,1], ncpu_512, elapsed_512, data_512,
           omp_512, '(a)', r'$N_\mathrm{rank}$')

_bbox = dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.85, edgecolor='none')
axes[0,0].text(0.03, 0.97, r'(a) Cosmo512, $N_\mathrm{thread}=1$: timers',
               transform=axes[0,0].transAxes, fontsize=8, fontweight='bold',
               va='top', ha='left', bbox=_bbox, zorder=20)
axes[0,1].text(0.03, 0.97, r'(b) Cosmo512, $N_\mathrm{thread}=1$: speedup',
               transform=axes[0,1].transAxes, fontsize=8, fontweight='bold',
               va='top', ha='left', bbox=_bbox, zorder=20)

# --- Cosmo1024 (bottom row) ---
plot_panel(axes[1,0], axes[1,1], ncpu_1024, elapsed_1024, data_1024,
           omp_1024, '(c)', r'$N_\mathrm{rank}$')

axes[1,0].text(0.03, 0.97, r'(c) Cosmo1024, $N_\mathrm{thread}=8$: timers',
               transform=axes[1,0].transAxes, fontsize=8, fontweight='bold',
               va='top', ha='left', bbox=_bbox, zorder=20)
axes[1,1].text(0.03, 0.97, r'(d) Cosmo1024, $N_\mathrm{thread}=8$: speedup',
               transform=axes[1,1].transAxes, fontsize=8, fontweight='bold',
               va='top', ha='left', bbox=_bbox, zorder=20)

plt.tight_layout()
outbase = '/gpfs/kjhan/Hydro/CUBE_HR5/code_cube/misc/fig_strong_scaling'
plt.savefig(outbase + '.pdf', bbox_inches='tight', dpi=150)
plt.savefig(outbase + '.png', bbox_inches='tight', dpi=150)
print(f'Saved {outbase}.pdf and .png')
