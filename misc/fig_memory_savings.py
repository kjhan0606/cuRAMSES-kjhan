import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np

fig, ax = plt.subplots(figsize=(7.0, 2.5))

categories = ['nbor\nelimination', 'hilbert_key\nelimination', 'bisec_ind_cell\non-demand', 
              'cell_level\non-demand', 'Defrag\nscratch']
values = [240, 640, 160, 160, 40]
colors = ['#3B82F6', '#8B5CF6', '#10B981', '#F59E0B', '#EC4899']
hatches = ['', '', '//', '//', '']

# Horizontal stacked bar
left = 0
bars = []
for i, (v, c, h) in enumerate(zip(values, colors, hatches)):
    bar = ax.barh(0, v, left=left, color=c, edgecolor='white', linewidth=1.5, 
                  height=0.5, alpha=0.85, hatch=h)
    # Label inside bar
    cx = left + v/2
    if v >= 100:
        ax.text(cx, 0, f'{v}\nMB', ha='center', va='center', fontsize=7, 
                fontweight='bold', color='white')
    else:
        ax.text(cx, 0, f'{v}', ha='center', va='center', fontsize=6, 
                fontweight='bold', color='white')
    left += v

# Hash table overhead (negative)
overhead = 50
ax.barh(0, -overhead, left=left, color='#EF4444', edgecolor='white', linewidth=1.5,
        height=0.5, alpha=0.7)
ax.text(left - overhead/2, 0, f'-{overhead}\nMB', ha='center', va='center', fontsize=6, 
        fontweight='bold', color='white')

# Net savings line
net = sum(values) - overhead
ax.axvline(x=net, color='#DC2626', linewidth=1.5, linestyle='--', alpha=0.8)
ax.text(net, 0.35, f'Net savings\n> {net} MB', ha='center', va='bottom', fontsize=7, 
        color='#DC2626', fontweight='bold')

# Legend
legend_labels = ['nbor array (always)', 'hilbert_key (k-section)', 
                 'bisec_ind_cell (between LB)', 'cell_level (between LB)',
                 'Defrag scratch', 'Hash table overhead']
legend_colors = colors + ['#EF4444']
legend_patches = [patches.Patch(facecolor=c, alpha=0.85, label=l) 
                  for c, l in zip(legend_colors, legend_labels)]
ax.legend(handles=legend_patches, loc='upper center', bbox_to_anchor=(0.5, -0.25),
          ncol=3, fontsize=6, frameon=False)

ax.set_xlabel('Memory (MB)', fontsize=9)
ax.set_xlim(-70, 1300)
ax.set_yticks([])
ax.set_title(r'Memory savings per MPI rank ($N_{\rm gridmax} = 5\,\mathrm{M}$)', 
             fontsize=10, fontweight='bold', pad=10)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['left'].set_visible(False)

plt.tight_layout()
plt.savefig('/gpfs/kjhan/Hydro/CUBE_HR5/code_cube/misc/fig_memory_savings.pdf',
            dpi=300, bbox_inches='tight', pad_inches=0.05)
plt.close()
print('fig_memory_savings.pdf generated successfully')
