import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np

fig, axes = plt.subplots(2, 1, figsize=(7.0, 3.0), gridspec_kw={'hspace': 0.6})

# Colors
c_red = '#EF4444'
c_blue = '#3B82F6'
c_darkblue = '#1D4ED8'
c_orange = '#F59E0B'
c_green = '#10B981'
c_gray = '#9CA3AF'
c_purple = '#8B5CF6'
c_dark = '#1E293B'

def draw_pipeline(ax, steps, title, n_exch, label_color, show_savings=False):
    ax.set_xlim(-0.5, len(steps)*1.1 + 0.5)
    ax.set_ylim(-0.5, 1.5)
    ax.axis('off')
    
    ax.text(-0.3, 1.2, title, fontsize=9, fontweight='bold', color=c_dark, va='center')
    
    x = 0
    for i, (name, stype, color) in enumerate(steps):
        if stype == 'compute':
            w = max(0.6, len(name)*0.08)
            rect = patches.FancyBboxPatch((x, 0.3), w, 0.7, boxstyle='round,pad=0.05',
                                           facecolor=color, edgecolor='white', linewidth=1, alpha=0.85)
            ax.add_patch(rect)
            ax.text(x + w/2, 0.65, name, fontsize=5.5, ha='center', va='center', 
                    color='white', fontweight='bold')
            x += w + 0.15
        elif stype == 'exchange':
            diamond = patches.RegularPolygon((x + 0.25, 0.65), 4, radius=0.25,
                                              facecolor=c_orange, edgecolor='#D97706', 
                                              linewidth=1, alpha=0.9)
            ax.add_patch(diamond)
            ax.text(x + 0.25, 0.65, 'E', fontsize=6, ha='center', va='center', 
                    color='white', fontweight='bold')
            ax.text(x + 0.25, 0.05, name, fontsize=4.5, ha='center', va='center', 
                    color='#D97706', rotation=0)
            x += 0.65
        elif stype == 'removed':
            diamond = patches.RegularPolygon((x + 0.25, 0.65), 4, radius=0.25,
                                              facecolor='#FEF3C7', edgecolor=c_gray, 
                                              linewidth=1, alpha=0.5, linestyle='--')
            ax.add_patch(diamond)
            ax.plot([x+0.05, x+0.45], [0.45, 0.85], color=c_red, linewidth=1.5, alpha=0.7)
            ax.plot([x+0.05, x+0.45], [0.85, 0.45], color=c_red, linewidth=1.5, alpha=0.7)
            x += 0.65
        
        if i < len(steps) - 1:
            ax.annotate('', xy=(x+0.05, 0.65), xytext=(x-0.1, 0.65),
                        arrowprops=dict(arrowstyle='->', color=c_gray, lw=0.8))
            x += 0.1
    
    ax.text(x + 0.3, 0.65, f'{n_exch} exchanges', fontsize=8, va='center',
            color=label_color, fontweight='bold',
            bbox=dict(boxstyle='round,pad=0.2', facecolor='white', edgecolor=label_color, linewidth=1.2))
    if show_savings:
        ax.text(x + 0.3, 0.2, '(−44%)', fontsize=7, va='center', color=c_green, fontweight='bold')

orig_simple = [
    ('Red', 'compute', c_red),
    ('φ', 'exchange', None),
    ('Black', 'compute', c_darkblue),
    ('φ', 'exchange', None),
    ('Resid.', 'compute', c_purple),
    ('r', 'exchange', None),
    ('Norm', 'compute', '#6B7280'),
    ('r', 'exchange', None),
    ('Restr.', 'compute', c_green),
    ('φ', 'exchange', None),
    ('Prolog.', 'compute', '#EC4899'),
    ('φ', 'exchange', None),
    (' ... ', 'compute', '#6B7280'),
]

opt_simple = [
    ('Red', 'compute', c_red),
    ('', 'removed', None),
    ('Black', 'compute', c_darkblue),
    ('φ', 'exchange', None),
    ('Resid.+Norm', 'compute', c_purple),
    ('', 'removed', None),
    ('r', 'exchange', None),
    ('Restr.', 'compute', c_green),
    ('φ', 'exchange', None),
    ('Prolog.', 'compute', '#EC4899'),
    ('φ', 'exchange', None),
    ('φ', 'exchange', None),
]

draw_pipeline(axes[0], orig_simple, '(a) Original', 9, c_orange)
draw_pipeline(axes[1], opt_simple, '(b) Optimized', 5, c_green, show_savings=True)

fig.suptitle('Multigrid V-cycle exchange optimization', fontsize=10, fontweight='bold', y=0.98)

plt.savefig('/gpfs/kjhan/Hydro/CUBE_HR5/code_cube/misc/fig_mg_redblack.pdf',
            dpi=300, bbox_inches='tight', pad_inches=0.05)
plt.close()
print('fig_mg_redblack.pdf generated successfully')
