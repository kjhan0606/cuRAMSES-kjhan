import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

# Data from Table 4
N = np.array([500, 1000, 5000, 10000, 50000, 100000])
seq  = np.array([0.54, 2.19, 54.96, 220.26, 5762.31, 22584.65])
tree = np.array([0.37, 0.82, 5.42, 11.92, 75.10, 163.41])
omp  = np.array([0.21, 0.37, 2.07, 3.24, 25.65, 37.92])

fig, ax = plt.subplots(1, 1, figsize=(3.5, 3.0))

ax.loglog(N, seq,  'o-', color='#1D4ED8', ms=5, lw=1.5, label='Sequential brute-force')
ax.loglog(N, tree, 's-', color='#EA580C', ms=5, lw=1.5, label=r'Serial oct-tree')
ax.loglog(N, omp,  '^-', color='#16A34A', ms=5, lw=1.5, label=r'OpenMP oct-tree (12 threads)')

# O(N^2) reference dashed line, anchored to sequential data
N_ref = np.array([300, 200000])
# Anchor at N=1000: seq=2.19ms, scale as (N/1000)^2
t_ref = 2.19 * (N_ref / 1000.0)**2
ax.loglog(N_ref, t_ref, '--', color='#1D4ED8', lw=0.8, alpha=0.5, label=r'$O(N^2)$')

# O(N log N) reference dashed line near tree data
t_nlogn = 0.82 * (N_ref / 1000.0) * np.log2(N_ref) / np.log2(1000)
ax.loglog(N_ref, t_nlogn, '--', color='#EA580C', lw=0.8, alpha=0.5, label=r'$O(N\log N)$')

ax.set_xlabel(r'$N_{\rm sink}$', fontsize=10)
ax.set_ylabel('Wall-clock time (ms)', fontsize=10)
ax.set_xlim(300, 200000)
ax.set_ylim(0.1, 100000)
ax.legend(fontsize=6, loc='upper left', framealpha=0.9)
ax.tick_params(labelsize=8)
ax.grid(True, alpha=0.3, which='both')

plt.tight_layout()
plt.savefig('/gpfs/kjhan/Hydro/CUBE_HR5/code_cube/misc/fig_fof_scaling.pdf',
            dpi=300, bbox_inches='tight')
plt.close()
print("fig_fof_scaling.pdf generated successfully")
