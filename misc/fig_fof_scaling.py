import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

# Full data from benchmark (12 OMP threads)
N    = np.array([50, 100, 200, 500, 1000, 2000, 5000, 10000, 50000, 100000])
seq  = np.array([0.007, 0.024, 0.089, 0.543, 2.188, 8.681, 54.284, 218.233, 5666.133, 22584.649])
tree = np.array([0.019, 0.053, 0.118, 0.361, 0.805, 1.902, 5.356, 11.802, 75.802, 163.410])
omp  = np.array([0.134, 0.153, 0.135, 0.210, 0.369, 0.779, 2.069, 3.239, 25.651, 37.917])

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(3.5, 6.0))

# --- Panel (a): Wall-clock time (full range) ---
ax1.loglog(N, seq,  'o-', color='#1D4ED8', ms=5, lw=1.5, label='Sequential brute-force')
ax1.loglog(N, tree, 's-', color='#EA580C', ms=5, lw=1.5, label='Serial oct-tree')
ax1.loglog(N, omp,  '^-', color='#16A34A', ms=5, lw=1.5, label='OpenMP oct-tree (12 thr)')

# O(N^2) reference
N_ref = np.array([30, 200000])
t_ref = 2.19 * (N_ref / 1000.0)**2
ax1.loglog(N_ref, t_ref, '--', color='#1D4ED8', lw=0.8, alpha=0.5, label=r'$O(N^2)$')

# O(N log N) reference
t_nlogn = 0.82 * (N_ref / 1000.0) * np.log2(N_ref) / np.log2(1000)
ax1.loglog(N_ref, t_nlogn, '--', color='#EA580C', lw=0.8, alpha=0.5, label=r'$O(N\log N)$')

# Shade the overhead region (N < 500)
ax1.axvspan(30, 500, alpha=0.08, color='red')
ax1.text(120, 8000, 'tree\noverhead', fontsize=7, color='red', alpha=0.7, ha='center')

ax1.set_xlabel(r'$N_{\rm sink}$', fontsize=10)
ax1.set_ylabel('Wall-clock time (ms)', fontsize=10)
ax1.set_xlim(30, 200000)
ax1.set_ylim(0.003, 100000)
ax1.legend(fontsize=5.5, loc='upper left', framealpha=0.9)
ax1.tick_params(labelsize=8)
ax1.grid(True, alpha=0.3, which='both')
ax1.text(0.02, 0.95, '(a)', transform=ax1.transAxes, fontsize=10, fontweight='bold', va='top')

# --- Panel (b): Speedup ratio ---
speedup_tree = seq / tree
speedup_omp  = seq / omp

ax2.loglog(N, speedup_tree, 's-', color='#EA580C', ms=5, lw=1.5, label='Serial tree / brute-force')
ax2.loglog(N, speedup_omp,  '^-', color='#16A34A', ms=5, lw=1.5, label='OMP tree / brute-force')
ax2.axhline(y=1.0, color='gray', ls='--', lw=0.8, alpha=0.7)

# Shade overhead region
ax2.axvspan(30, 500, alpha=0.08, color='red')

# Annotate crossover points
ax2.annotate(f'OMP crossover\n$N_{{\\rm sink}}\\approx 300$',
             xy=(300, 1.0), xytext=(80, 100),
             fontsize=6.5, ha='center',
             arrowprops=dict(arrowstyle='->', color='#16A34A', lw=1.0),
             color='#16A34A')

ax2.set_xlabel(r'$N_{\rm sink}$', fontsize=10)
ax2.set_ylabel('Speedup over brute-force', fontsize=10)
ax2.set_xlim(30, 200000)
ax2.set_ylim(0.03, 1000)
ax2.legend(fontsize=6.5, loc='upper left', framealpha=0.9)
ax2.tick_params(labelsize=8)
ax2.grid(True, alpha=0.3, which='both')
ax2.text(0.02, 0.95, '(b)', transform=ax2.transAxes, fontsize=10, fontweight='bold', va='top')

plt.tight_layout()
plt.savefig('/gpfs/kjhan/Hydro/CUBE_HR5/code_cube/misc/fig_fof_scaling.pdf',
            dpi=300, bbox_inches='tight')
plt.savefig('/gpfs/kjhan/Hydro/CUBE_HR5/code_cube/misc/fig_fof_scaling.png',
            dpi=150, bbox_inches='tight')
plt.close()
print("fig_fof_scaling.pdf/png generated successfully")
