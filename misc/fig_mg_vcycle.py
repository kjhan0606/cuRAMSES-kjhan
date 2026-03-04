import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np

fig, ax = plt.subplots(1, 1, figsize=(3.5, 3.0))

# V-cycle levels (4 levels: fine at top, coarse at bottom)
n_levels = 4
level_labels = [f'Level {n_levels - i}' for i in range(n_levels)]
y_positions = np.linspace(2.8, 0.4, n_levels)

# V-shape x positions
x_left = np.linspace(0.5, 2.0, n_levels)   # descending (restriction)
x_right = np.linspace(4.5, 3.0, n_levels)  # ascending (prolongation)

c_restrict = '#3B82F6'   # blue for restriction
c_prolong  = '#10B981'   # green for prolongation
c_smooth   = '#EF4444'   # red for smoothing boxes
c_solve    = '#8B5CF6'   # purple for direct solve
c_dark     = '#1E293B'

box_w, box_h = 0.55, 0.22

# Draw restriction arrows (left leg, top to bottom)
for i in range(n_levels - 1):
    ax.annotate('', xy=(x_left[i+1] + box_w/2, y_positions[i+1] + box_h),
                xytext=(x_left[i] + box_w/2, y_positions[i]),
                arrowprops=dict(arrowstyle='->', color=c_restrict, lw=1.5))

# Draw prolongation arrows (right leg, bottom to top)
for i in range(n_levels - 1):
    ax.annotate('', xy=(x_right[i] + box_w/2, y_positions[i]),
                xytext=(x_right[i+1] + box_w/2, y_positions[i+1] + box_h),
                arrowprops=dict(arrowstyle='->', color=c_prolong, lw=1.5))

# Draw smoothing boxes on both legs
for i in range(n_levels):
    # Left leg (pre-smoothing)
    if i < n_levels - 1:
        rect = patches.FancyBboxPatch((x_left[i], y_positions[i]),
                                       box_w, box_h, boxstyle='round,pad=0.03',
                                       facecolor='#FEE2E2', edgecolor=c_smooth, linewidth=1.2)
        ax.add_patch(rect)
        ax.text(x_left[i] + box_w/2, y_positions[i] + box_h/2, 'Smooth',
                ha='center', va='center', fontsize=6, color=c_smooth, fontweight='bold')
    else:
        # Direct solve at coarsest
        rect = patches.FancyBboxPatch((x_left[i], y_positions[i]),
                                       box_w, box_h, boxstyle='round,pad=0.03',
                                       facecolor='#EDE9FE', edgecolor=c_solve, linewidth=1.2)
        ax.add_patch(rect)
        ax.text(x_left[i] + box_w/2, y_positions[i] + box_h/2, 'Solve',
                ha='center', va='center', fontsize=6, color=c_solve, fontweight='bold')

    # Right leg (post-smoothing), skip coarsest
    if i < n_levels - 1:
        rect = patches.FancyBboxPatch((x_right[i], y_positions[i]),
                                       box_w, box_h, boxstyle='round,pad=0.03',
                                       facecolor='#FEE2E2', edgecolor=c_smooth, linewidth=1.2)
        ax.add_patch(rect)
        ax.text(x_right[i] + box_w/2, y_positions[i] + box_h/2, 'Smooth',
                ha='center', va='center', fontsize=6, color=c_smooth, fontweight='bold')

# Level labels on the right
for i in range(n_levels):
    ax.text(5.3, y_positions[i] + box_h/2, level_labels[i],
            ha='left', va='center', fontsize=7, color=c_dark)

# Leg labels
ax.text(1.0, 3.1, 'Restrict', ha='center', fontsize=7, color=c_restrict, fontweight='bold')
ax.text(4.0, 3.1, 'Prolongate', ha='center', fontsize=7, color=c_prolong, fontweight='bold')

# Exchange annotation
ax.annotate('5 exchanges\n(was 9)',
            xy=(x_left[0] + box_w + 0.05, y_positions[0] + box_h/2),
            xytext=(x_left[0] + box_w + 0.5, y_positions[0] + box_h/2 + 0.35),
            fontsize=5.5, color='#6B7280', ha='center',
            arrowprops=dict(arrowstyle='->', color='#9CA3AF', lw=0.8))

ax.set_xlim(-0.1, 6.2)
ax.set_ylim(-0.1, 3.4)
ax.set_aspect('equal')
ax.axis('off')

plt.tight_layout()
plt.savefig('/gpfs/kjhan/Hydro/CUBE_HR5/code_cube/misc/fig_mg_vcycle.pdf',
            dpi=300, bbox_inches='tight', pad_inches=0.05)
plt.close()
print("fig_mg_vcycle.pdf generated successfully")
