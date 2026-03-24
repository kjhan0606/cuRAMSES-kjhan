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
# Parse per-step wall-clock times from logs
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

t256, mu256 = parse_step_times(os.path.join(weakdir, 'log_256.log'))
t512, mu512 = parse_step_times(os.path.join(weakdir, 'log_512.log'))
t1024, mu1024 = parse_step_times(os.path.join(weakdir, 'log_1024.log'))

# Configuration
configs = {
    r'$256^3$': dict(ncpu=2, nodes=1, cores=16, times=t256, muspt=mu256),
    r'$512^3$': dict(ncpu=16, nodes=2, cores=128, times=t512, muspt=mu512),
    r'$1024^3$': dict(ncpu=128, nodes=16, cores=1024, times=t1024, muspt=mu1024),
}

# Skip first step (includes initialization overhead)
skip = 1

# ============================================================
# Figure: 2-panel
# Left: wall-clock per step vs step number
# Right: bar chart of average mus/pt
# ============================================================
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4.5))

colors = ['#1f77b4', '#ff7f0e', '#2ca02c']
markers = ['o', 's', '^']

# --- Left panel: wall-clock per step ---
for i, (label, cfg) in enumerate(configs.items()):
    steps = np.arange(1, len(cfg['times']) + 1)
    ax1.plot(steps[skip:], cfg['times'][skip:],
             marker=markers[i], color=colors[i], label=label,
             linewidth=1.5, markersize=4, markevery=2)

ax1.set_xlabel('Coarse step')
ax1.set_ylabel('Wall-clock time per step (s)')
ax1.legend(frameon=True, fontsize=9)
ax1.set_title('Per-step wall-clock time')
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
ax2.set_title('Weak scaling efficiency')
ax2.legend(frameon=True, fontsize=9)
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
