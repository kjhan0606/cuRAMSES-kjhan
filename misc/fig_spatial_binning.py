import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np

np.random.seed(42)
fig, axes = plt.subplots(1, 2, figsize=(7.0, 3.5))

n_sn = 40
sn_x = np.random.uniform(0.05, 0.95, n_sn)
sn_y = np.random.uniform(0.05, 0.95, n_sn)

# Target cell position
cx, cy = 0.45, 0.55

c_blue = '#3B82F6'
c_red = '#EF4444'
c_green = '#10B981'
c_orange = '#F59E0B'
c_gray = '#D1D5DB'
c_dark = '#1E293B'

# === Panel (a): Brute force ===
ax = axes[0]
ax.set_xlim(0, 1)
ax.set_ylim(0, 1)
ax.set_aspect('equal')
ax.set_title('(a) Brute force', fontsize=9, fontweight='bold', pad=8)

# Draw lines from cell to ALL SNe
for i in range(n_sn):
    ax.plot([cx, sn_x[i]], [cy, sn_y[i]], color=c_red, linewidth=0.4, alpha=0.4)

# Plot SN events
ax.scatter(sn_x, sn_y, marker='*', s=50, c=c_orange, edgecolors='#B45309', 
           linewidth=0.5, zorder=5)

# Highlight target cell
ax.scatter([cx], [cy], marker='s', s=80, c=c_blue, edgecolors='white', 
           linewidth=1.5, zorder=6)

# Complexity label
ax.text(0.5, -0.08, r'$\mathcal{O}(N_{\rm cells} \times N_{\rm SN})$', fontsize=9,
        ha='center', color=c_red, fontweight='bold', transform=ax.transAxes)

ax.set_xticks([])
ax.set_yticks([])
for spine in ax.spines.values():
    spine.set_edgecolor(c_gray)

# === Panel (b): Spatial binning ===
ax = axes[1]
ax.set_xlim(0, 1)
ax.set_ylim(0, 1)
ax.set_aspect('equal')
ax.set_title('(b) Spatial hash binning', fontsize=9, fontweight='bold', pad=8)

# Draw bin grid
n_bins = 8
for i in range(n_bins + 1):
    ax.axhline(y=i/n_bins, color=c_gray, linewidth=0.5, alpha=0.5)
    ax.axvline(x=i/n_bins, color=c_gray, linewidth=0.5, alpha=0.5)

# Find which bin the target cell is in
target_bx = int(cx * n_bins)
target_by = int(cy * n_bins)

# Shade 3x3 neighboring bins
for dx in range(-1, 2):
    for dy in range(-1, 2):
        bx = target_bx + dx
        by = target_by + dy
        if 0 <= bx < n_bins and 0 <= by < n_bins:
            color = '#DBEAFE' if (dx == 0 and dy == 0) else '#EFF6FF'
            rect = patches.Rectangle((bx/n_bins, by/n_bins), 1/n_bins, 1/n_bins,
                                      facecolor=color, edgecolor=c_blue, 
                                      linewidth=0.8, alpha=0.7)
            ax.add_patch(rect)

# Determine which SNe are in the 3x3 neighborhood
nearby_mask = np.zeros(n_sn, dtype=bool)
for i in range(n_sn):
    sbx = int(sn_x[i] * n_bins)
    sby = int(sn_y[i] * n_bins)
    if abs(sbx - target_bx) <= 1 and abs(sby - target_by) <= 1:
        nearby_mask[i] = True

# Draw lines only to nearby SNe
for i in range(n_sn):
    if nearby_mask[i]:
        ax.plot([cx, sn_x[i]], [cy, sn_y[i]], color=c_green, linewidth=1.0, alpha=0.7)

# Plot distant SNe (grayed out)
ax.scatter(sn_x[~nearby_mask], sn_y[~nearby_mask], marker='*', s=30, 
           c=c_gray, edgecolors='#9CA3AF', linewidth=0.3, zorder=4, alpha=0.5)

# Plot nearby SNe (highlighted)
ax.scatter(sn_x[nearby_mask], sn_y[nearby_mask], marker='*', s=60, 
           c=c_orange, edgecolors='#B45309', linewidth=0.5, zorder=5)

# Highlight target cell
ax.scatter([cx], [cy], marker='s', s=80, c=c_blue, edgecolors='white', 
           linewidth=1.5, zorder=6)

# Label the 3x3 region
ax.text(0.5, -0.08, r'$\mathcal{O}(N_{\rm cells} \times 27\,\bar{n}_{\rm SN/bin})$', 
        fontsize=9, ha='center', color=c_green, fontweight='bold', transform=ax.transAxes)

# Add "3x3 bins" annotation
mid_x = (target_bx + 0.5) / n_bins
mid_y = (target_by + 1.5) / n_bins + 0.02
ax.annotate('9 bins\n(27 in 3D)', xy=(mid_x, mid_y), fontsize=6, ha='center',
            color=c_blue, fontweight='bold')

ax.set_xticks([])
ax.set_yticks([])
for spine in ax.spines.values():
    spine.set_edgecolor(c_gray)

plt.tight_layout()
plt.savefig('/gpfs/kjhan/Hydro/CUBE_HR5/code_cube/misc/fig_spatial_binning.pdf',
            dpi=300, bbox_inches='tight', pad_inches=0.05)
plt.close()
print("fig_spatial_binning.pdf generated successfully")
