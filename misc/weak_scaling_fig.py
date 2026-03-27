#!/usr/bin/env python3
"""Generate weak scaling figure for cuRAMSES paper.

Weak scaling: 256^3 (2 ranks, 1 node), 512^3 (16 ranks, 2 nodes),
              1024^3 (128 ranks, 16 nodes).
Same cosmology, same AMR depth (levelmax = levelmin+5), OMP=8.
~8.4M cells/rank (constant work per rank).
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import re, os

# ============================================================
# Parse per-step wall-clock times and load_balance times from logs
# ============================================================
weakdir = '/gpfs/kjhan/Hydro/CUBE_HR5/code_cube/test_cosmo1024/scaling/weak'

def parse_step_times(logfile):
    """Extract per-coarse-step wall-clock times and mus/pt."""
    times = []
    muspt = []
    with open(logfile) as f:
        for line in f:
            m = re.search(r'Time elapsed since last coarse step:\s+([\d.]+)\s+s\s+([\d.]+)\s+mus/pt', line)
            if m:
                times.append(float(m.group(1)))
                muspt.append(float(m.group(2)))
    return np.array(times), np.array(muspt)

def parse_remap_times(logfile):
    """Extract per-event load_balance times."""
    remap = []
    with open(logfile) as f:
        for line in f:
            m = re.search(r'load_balance total:\s+([\d.]+)\s+s', line)
            if m:
                remap.append(float(m.group(1)))
    return np.array(remap)

def subtract_remap(times, remap_times):
    """Subtract remap time from steps that contain a remap event.

    Remap steps are detected as those with time > 1.5 * median(times).
    The remap_times are matched to remap steps in order.
    """
    adjusted = times.copy()
    median_t = np.median(times[1:])  # exclude first step (init)
    remap_mask = times > 1.5 * median_t
    remap_indices = np.where(remap_mask)[0]
    n = min(len(remap_indices), len(remap_times))
    for i in range(n):
        adjusted[remap_indices[i]] -= remap_times[i]
    return adjusted

t256, mu256 = parse_step_times(os.path.join(weakdir, 'log_256.log'))
t512, mu512 = parse_step_times(os.path.join(weakdir, 'log_512.log'))
t1024, mu1024 = parse_step_times(os.path.join(weakdir, 'log_1024.log'))

r256 = parse_remap_times(os.path.join(weakdir, 'log_256.log'))
r512 = parse_remap_times(os.path.join(weakdir, 'log_512.log'))
r1024 = parse_remap_times(os.path.join(weakdir, 'log_1024.log'))

t256_adj = subtract_remap(t256, r256)
t512_adj = subtract_remap(t512, r512)
t1024_adj = subtract_remap(t1024, r1024)

# Configuration
configs = {
    r'$256^3$': dict(ncpu=2, nodes=1, cores=16, times=t256, times_adj=t256_adj, muspt=mu256),
    r'$512^3$': dict(ncpu=16, nodes=2, cores=128, times=t512, times_adj=t512_adj, muspt=mu512),
    r'$1024^3$': dict(ncpu=128, nodes=16, cores=1024, times=t1024, times_adj=t1024_adj, muspt=mu1024),
}

# Skip first step (includes initialization overhead)
skip = 1

# ============================================================
# Figure: 2-panel
# Left: wall-clock per step vs step number (with remap subtracted)
# Right: bar chart of average mus/pt
# ============================================================
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4.5))
_bbox = dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.85, edgecolor='none')

colors = ['#1f77b4', '#ff7f0e', '#2ca02c']
markers = ['o', 's', '^']

# --- Left panel: wall-clock per step (remap subtracted) ---
for i, (label, cfg) in enumerate(configs.items()):
    steps = np.arange(1, len(cfg['times']) + 1)
    # Show remap-subtracted times as solid line
    ax1.plot(steps[skip:], cfg['times_adj'][skip:],
             marker=markers[i], color=colors[i], label=label,
             linewidth=1.5, markersize=4, markevery=2)
    # Show original times as faint line for reference
    ax1.plot(steps[skip:], cfg['times'][skip:],
             color=colors[i], linewidth=0.6, alpha=0.3, linestyle='--')

ax1.set_xlabel('Coarse step')
ax1.set_ylabel('Wall-clock time per step (s)')
ax1.legend(frameon=True, fontsize=9, loc='upper left',
           bbox_to_anchor=(0.0, 0.87))
ax1.text(0.03, 0.97, '(a) Per-step time (remap subtracted)',
         transform=ax1.transAxes, fontsize=8, fontweight='bold',
         va='top', ha='left', bbox=_bbox, zorder=20)
ax1.grid(True, alpha=0.3)

# --- Right panel: average mus/pt (efficiency metric) ---
labels_list = []
avg_muspt = []
avg_muspt_std = []
ncores_list = []

for label, cfg in configs.items():
    labels_list.append(label)
    mu = cfg['muspt'][skip:]
    avg_muspt.append(np.mean(mu))
    avg_muspt_std.append(np.std(mu))
    ncores_list.append(cfg['cores'])

x = np.arange(len(labels_list))
bars = ax2.bar(x, avg_muspt, yerr=avg_muspt_std, width=0.5,
               color=colors[:len(labels_list)], edgecolor='black',
               capsize=4, alpha=0.85)

# Add core count annotations
for xi, nc, mu in zip(x, ncores_list, avg_muspt):
    ax2.text(xi, mu + 0.3, f'{nc} cores', ha='center', va='bottom', fontsize=8)

# Ideal line (constant mus/pt = perfect weak scaling)
ax2.axhline(y=avg_muspt[0], color='k', linestyle='--', alpha=0.4,
            label=f'Ideal ({avg_muspt[0]:.2f} $\\mu$s/pt)')

ax2.set_xticks(x)
ax2.set_xticklabels(labels_list)
ax2.set_ylabel(r'Average $\mu$s / grid-point')
ax2.text(0.03, 0.97, '(b) Weak scaling efficiency',
         transform=ax2.transAxes, fontsize=8, fontweight='bold',
         va='top', ha='left', bbox=_bbox, zorder=20)
ax2.legend(frameon=True, fontsize=9, loc='upper left',
           bbox_to_anchor=(0.0, 0.87))
ax2.grid(True, alpha=0.3, axis='y')

# Compute and print efficiency
print("Weak scaling summary:")
print(f"  256^3:  {avg_muspt[0]:.2f} +/- {avg_muspt_std[0]:.2f} mus/pt  ({ncores_list[0]} cores)")
print(f"  512^3:  {avg_muspt[1]:.2f} +/- {avg_muspt_std[1]:.2f} mus/pt  ({ncores_list[1]} cores)")
print(f"  1024^3: {avg_muspt[2]:.2f} +/- {avg_muspt_std[2]:.2f} mus/pt  ({ncores_list[2]} cores)")
eff_512 = avg_muspt[0] / avg_muspt[1] * 100
eff_1024 = avg_muspt[0] / avg_muspt[2] * 100
print(f"  Efficiency 512^3:  {eff_512:.1f}%")
print(f"  Efficiency 1024^3: {eff_1024:.1f}%")

plt.tight_layout()
outfile = '/gpfs/kjhan/Hydro/CUBE_HR5/code_cube/misc/fig_weak_scaling.png'
plt.savefig(outfile, dpi=200, bbox_inches='tight')
print(f"Saved: {outfile}")

# Also save PDF for paper
outpdf = outfile.replace('.png', '.pdf')
plt.savefig(outpdf, bbox_inches='tight')
print(f"Saved: {outpdf}")
