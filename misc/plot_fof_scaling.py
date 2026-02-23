#!/usr/bin/env python3
"""Generate FoF scaling figure for paper."""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

# Data from benchmark (3 MPI x 4 OMP threads for OMP)
nsink = np.array([50, 100, 200, 500, 1000, 2000, 5000, 10000, 50000, 100000])
seq   = np.array([0.007, 0.024, 0.089, 0.543, 2.188, 8.681, 54.284, 218.233, 5666.133, 22584.649])
tree  = np.array([0.019, 0.053, 0.118, 0.361, 0.805, 1.902, 5.356, 11.802, 75.802, 163.410])
omp   = np.array([0.134, 0.153, 0.135, 0.210, 0.369, 0.779, 2.069, 3.239, 25.651, 37.917])

fig, ax = plt.subplots(1, 1, figsize=(4.5, 3.5))

ax.loglog(nsink, seq,  'o-', color='#2166ac', markersize=4, linewidth=1.5, label='Sequential $O(N^2)$')
ax.loglog(nsink, tree, 's-', color='#d95f02', markersize=4, linewidth=1.5, label='Serial oct-tree')
ax.loglog(nsink, omp,  '^-', color='#1b9e77', markersize=4, linewidth=1.5, label=r'MPI+OMP tree (3$\times$4)')

# Reference slopes
n_ref = np.array([500, 100000])
ax.loglog(n_ref, 0.002*n_ref**2/500**2, '--', color='#2166ac', alpha=0.3, linewidth=1)
ax.loglog(n_ref, 0.3*n_ref*np.log2(n_ref)/(500*np.log2(500)), '--', color='#d95f02', alpha=0.3, linewidth=1)

ax.set_xlabel(r'$N_{\rm sink}$', fontsize=11)
ax.set_ylabel('Wall-clock time (ms)', fontsize=11)
ax.legend(fontsize=8, loc='upper left')
ax.set_xlim(30, 2e5)
ax.set_ylim(0.005, 5e4)
ax.grid(True, which='both', alpha=0.2)
ax.tick_params(labelsize=9)

plt.tight_layout()
plt.savefig('fig_fof_scaling.pdf', bbox_inches='tight')
plt.savefig('fig_fof_scaling.png', dpi=150, bbox_inches='tight')
print('Saved fig_fof_scaling.pdf and fig_fof_scaling.png')
