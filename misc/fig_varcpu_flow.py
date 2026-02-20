import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np

fig, ax = plt.subplots(figsize=(7.0, 3.5))
ax.set_xlim(0, 10)
ax.set_ylim(0, 4.5)
ax.axis('off')

# Colors
c_blue = '#3B82F6'
c_green = '#10B981'
c_purple = '#8B5CF6'
c_dark = '#1E293B'
c_light_blue = '#DBEAFE'
c_light_green = '#D1FAE5'
c_light_purple = '#EDE9FE'
c_orange = '#F59E0B'
c_gray = '#9CA3AF'

def draw_stage(ax, x, y, w, h, title, items, color, light_color, stage_num):
    # Main box
    rect = patches.FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.1",
                                   facecolor=light_color, edgecolor=color, linewidth=2)
    ax.add_patch(rect)

    # Stage number circle
    circle = patches.Circle((x + 0.35, y + h - 0.3), 0.22, facecolor=color,
                              edgecolor='white', linewidth=1.5, zorder=5)
    ax.add_patch(circle)
    ax.text(x + 0.35, y + h - 0.3, str(stage_num), fontsize=9, ha='center', va='center',
            color='white', fontweight='bold', zorder=6)

    # Title (placed right of the stage number circle)
    ax.text(x + 0.65, y + h - 0.3, title, fontsize=7, ha='left', va='center',
            fontweight='bold', color=color)

    # Items
    for i, item in enumerate(items):
        ax.text(x + 0.3, y + h - 0.8 - i*0.4, '\u2022', fontsize=8, va='center', color=color)
        ax.text(x + 0.5, y + h - 0.8 - i*0.4, item, fontsize=6, va='center', color=c_dark)

# Title
ax.text(5.0, 4.2, 'Binary variable-$N_{\\rm cpu}$ restart: 3-stage distributed I/O',
        fontsize=10, ha='center', fontweight='bold', color=c_dark)

# Stage dimensions
sw, sh = 2.7, 3.2
sy = 0.5
gap = 0.5

# Stage 1
s1x = 0.3
draw_stage(ax, s1x, sy, sw, sh,
           'Detection & Assignment',
           ['Probe header: read $N_{\\rm cpu}^{\\rm file}$',
            'Round-robin file mapping',
            '$(f-1)\\;\\mathrm{mod}\\;N_{\\rm cpu} = r-1$',
            'Each rank: \u2264 2 files'],
           c_blue, c_light_blue, 1)

# Small file icons for Stage 1
for i in range(4):
    fx = s1x + 0.3 + i*0.35
    fy = sy + 0.25
    frect = patches.Rectangle((fx, fy), 0.2, 0.25, facecolor='white',
                                edgecolor=c_blue, linewidth=0.5)
    ax.add_patch(frect)
    ax.text(fx+0.1, fy+0.12, 'f', fontsize=4, ha='center', va='center', color=c_blue)

# Arrow 1->2
ax.annotate('', xy=(s1x + sw + gap, sy + sh/2), xytext=(s1x + sw + 0.1, sy + sh/2),
            arrowprops=dict(arrowstyle='->', color=c_dark, lw=1.5))

# Stage 2
s2x = s1x + sw + gap
draw_stage(ax, s2x, sy, sw, sh,
           'AMR Reconstruction',
           ['Read grid positions $\\mathbf{x}_g$',
            'Ownership via cpumap',
            'MPI_ALLTOALLV exchange',
            'Build local AMR tree',
            'Store exchange metadata'],
           c_green, c_light_green, 2)

# Arrow 2->3
ax.annotate('', xy=(s2x + sw + gap, sy + sh/2), xytext=(s2x + sw + 0.1, sy + sh/2),
            arrowprops=dict(arrowstyle='->', color=c_dark, lw=1.5))

# Reuse arrow (curved, from stage 2 bottom to stage 3)
ax.annotate('reuse\nmetadata', xy=(s2x + sw + gap + 0.3, sy + 0.5),
            xytext=(s2x + sw/2, sy + 0.1),
            fontsize=5.5, color=c_orange, fontweight='bold', ha='center',
            arrowprops=dict(arrowstyle='->', color=c_orange, lw=1.2,
                          connectionstyle='arc3,rad=-0.3'))

# Stage 3
s3x = s2x + sw + gap
draw_stage(ax, s3x, sy, sw, sh,
           'Hydro & Gravity Scatter',
           ['Read hydro/gravity files',
            'Reuse send ordering',
            'Reuse ALLTOALLV counts',
            'Scatter to local cells',
            'Primitive \u2192 conservative'],
           c_purple, c_light_purple, 3)

plt.savefig('/gpfs/kjhan/Hydro/CUBE_HR5/code_cube/misc/fig_varcpu_flow.pdf',
            dpi=300, bbox_inches='tight', pad_inches=0.05)
plt.close()
print("fig_varcpu_flow.pdf generated successfully")
