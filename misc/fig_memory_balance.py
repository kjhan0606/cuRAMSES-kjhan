import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

ncpu = np.array([2, 4, 8, 12, 16, 24, 32, 48, 64])
m_min = np.array([87.0, 50.1, 28.0, 20.9, 16.7, 12.5, 10.7, 8.2, 8.2])
m_max = np.array([87.4, 50.5, 28.6, 21.4, 17.1, 13.1, 11.1, 8.6, 8.6])
m_ratio = np.array([1.005, 1.008, 1.021, 1.024, 1.024, 1.048, 1.037, 1.049, 1.049])
l10_ratio = np.array([1.13, 1.15, 1.22, 1.33, 1.35, 1.40, 1.36, 1.56, 1.61])

c_blue = '#3B82F6'
c_red = '#EF4444'
c_gray = '#9CA3AF'

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7.0, 3.0))

# === Panel (a): Memory per rank ===
ax1.fill_between(ncpu, m_min, m_max, alpha=0.3, color=c_blue,
                 label=r'$M_{\rm min}$' + u'\u2013' + r'$M_{\rm max}$')
ax1.plot(ncpu, m_min, 'o-', color=c_blue, markersize=4, linewidth=1.5,
         label=r'$M_{\rm min}$')
ax1.plot(ncpu, m_max, 's-', color='#1D4ED8', markersize=4, linewidth=1.5,
         label=r'$M_{\rm max}$')

ncpu_ideal = np.array([1, 2, 4, 8, 16, 32, 64])
m_ideal = 147.2 / ncpu_ideal
ax1.plot(ncpu_ideal, m_ideal, '--', color=c_gray, linewidth=1,
         label=r'Ideal $\propto 1/N_{\rm cpu}$')

ax1.axhline(y=4.0, color=c_red, linewidth=0.8, linestyle=':', alpha=0.6)
ax1.text(50, 4.5, 'Fixed overhead\n~4 GB', fontsize=5.5, color=c_red, ha='center')

ax1.set_xscale('log', base=2)
ax1.set_yscale('log')
ax1.set_xlabel(r'$N_{\rm cpu}$', fontsize=9)
ax1.set_ylabel('Memory per rank (GB)', fontsize=9)
ax1.set_title('(a) Memory per rank', fontsize=9, fontweight='bold')
ax1.set_xticks(ncpu)
ax1.set_xticklabels([str(n) for n in ncpu], fontsize=7)
ax1.tick_params(axis='y', labelsize=7)
ax1.legend(fontsize=5.5, loc='upper right', frameon=True, fancybox=True)
ax1.grid(True, alpha=0.2)

# === Panel (b): Imbalance ratios ===
ln1 = ax2.plot(ncpu, m_ratio, 'o-', color=c_blue, markersize=5, linewidth=1.5,
               label=r'Memory $M_{\rm max}/M_{\rm min}$', zorder=5)
ax2.axhline(y=1.05, color=c_blue, linewidth=0.8, linestyle='--', alpha=0.5)
ax2.text(3, 1.053, '5% imbalance', fontsize=5.5, color=c_blue, alpha=0.7)
ax2.set_ylabel('Memory imbalance ratio', fontsize=8, color=c_blue)
ax2.tick_params(axis='y', labelcolor=c_blue, labelsize=7)
ax2.set_ylim(0.99, 1.07)

ax2_r = ax2.twinx()
ln2 = ax2_r.plot(ncpu, l10_ratio, 's-', color=c_red, markersize=5, linewidth=1.5,
                  label='Level-10 grids max/min', zorder=4)
ax2_r.set_ylabel('Grid count imbalance ratio', fontsize=8, color=c_red)
ax2_r.tick_params(axis='y', labelcolor=c_red, labelsize=7)
ax2_r.set_ylim(1.0, 1.75)

ax2.set_xscale('log', base=2)
ax2.set_xlabel(r'$N_{\rm cpu}$', fontsize=9)
ax2.set_title('(b) Load balance quality', fontsize=9, fontweight='bold')
ax2.set_xticks(ncpu)
ax2.set_xticklabels([str(n) for n in ncpu], fontsize=7)
ax2.grid(True, alpha=0.2)

lns = ln1 + ln2
labs = [l.get_label() for l in lns]
ax2.legend(lns, labs, fontsize=5.5, loc='upper left', frameon=True, fancybox=True)

plt.tight_layout()
plt.savefig('/gpfs/kjhan/Hydro/CUBE_HR5/code_cube/misc/fig_memory_balance.pdf',
            dpi=300, bbox_inches='tight', pad_inches=0.05)
plt.close()
print("fig_memory_balance.pdf generated successfully")
