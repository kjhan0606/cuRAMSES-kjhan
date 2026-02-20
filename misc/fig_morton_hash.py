import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np

fig, ax = plt.subplots(1, 1, figsize=(7.0, 3.0))
ax.set_xlim(0, 10)
ax.set_ylim(0, 4)
ax.axis('off')

# Colors
c_blue = '#3B82F6'
c_light = '#DBEAFE'
c_red = '#EF4444'
c_gray = '#9CA3AF'
c_dark = '#1E3A5A'
c_green = '#10B981'

# === Left: Original nbor array (eliminated) ===
nbor_x, nbor_y = 0.3, 0.5
nbor_w, nbor_h = 1.2, 2.8
rect = patches.FancyBboxPatch((nbor_x, nbor_y), nbor_w, nbor_h, 
                                boxstyle="round,pad=0.05", 
                                facecolor=c_light, edgecolor=c_gray, linewidth=1.5)
ax.add_patch(rect)
ax.text(nbor_x + nbor_w/2, nbor_y + nbor_h + 0.15, 'Original', fontsize=8, 
        ha='center', fontweight='bold', color=c_dark)
ax.text(nbor_x + nbor_w/2, nbor_y + nbor_h/2 + 0.3, 'nbor', fontsize=9, 
        ha='center', fontfamily='monospace', color=c_dark)
ax.text(nbor_x + nbor_w/2, nbor_y + nbor_h/2 - 0.1, r'$(N_{\rm gridmax}$', fontsize=7, 
        ha='center', color=c_gray)
ax.text(nbor_x + nbor_w/2, nbor_y + nbor_h/2 - 0.45, r'$\times\;6)$', fontsize=7, 
        ha='center', color=c_gray)
ax.text(nbor_x + nbor_w/2, nbor_y + 0.25, '240 MB', fontsize=7, 
        ha='center', color=c_red, fontweight='bold')

# X mark
ax.plot([nbor_x+0.15, nbor_x+nbor_w-0.15], [nbor_y+0.7, nbor_y+nbor_h-0.3], 
        color=c_red, linewidth=3, alpha=0.6)
ax.plot([nbor_x+0.15, nbor_x+nbor_w-0.15], [nbor_y+nbor_h-0.3, nbor_y+0.7], 
        color=c_red, linewidth=3, alpha=0.6)

# === Right: Morton hash table flow ===
# Step 1: Grid cell with coordinates
gx, gy = 2.5, 2.2
rect1 = patches.FancyBboxPatch((gx, gy), 1.0, 0.7, boxstyle="round,pad=0.05",
                                 facecolor='#E0F2FE', edgecolor=c_blue, linewidth=1.2)
ax.add_patch(rect1)
# Draw mini 2x2 grid inside
for i in range(3):
    ax.plot([gx+0.25+i*0.17, gx+0.25+i*0.17], [gy+0.15, gy+0.55], color=c_blue, linewidth=0.5)
    ax.plot([gx+0.25, gx+0.25+0.34], [gy+0.15+i*0.2, gy+0.15+i*0.2], color=c_blue, linewidth=0.5)
# Highlight one cell
highlight = patches.Rectangle((gx+0.42, gy+0.35), 0.17, 0.2, facecolor=c_blue, alpha=0.4)
ax.add_patch(highlight)
ax.text(gx+0.5, gy-0.15, r'$(i_x, i_y, i_z)$', fontsize=7, ha='center', color=c_dark)

# Arrow 1
ax.annotate('', xy=(4.0, gy+0.35), xytext=(gx+1.05, gy+0.35),
            arrowprops=dict(arrowstyle='->', color=c_dark, lw=1.2))

# Step 2: Bit interleave
bx, by = 4.0, gy-0.0
rect2 = patches.FancyBboxPatch((bx, by), 1.6, 0.7, boxstyle="round,pad=0.05",
                                 facecolor='#FEF3C7', edgecolor='#F59E0B', linewidth=1.2)
ax.add_patch(rect2)
ax.text(bx+0.8, by+0.42, 'Bit interleave', fontsize=7, ha='center', color='#92400E', fontweight='bold')
ax.text(bx+0.8, by+0.15, r'$x_0 y_0 z_0 x_1 y_1 z_1 \ldots$', fontsize=6, ha='center', color='#92400E')

# Arrow 2
ax.annotate('', xy=(6.1, by+0.35), xytext=(bx+1.65, by+0.35),
            arrowprops=dict(arrowstyle='->', color=c_dark, lw=1.2))

# Step 3: Morton key
mx, my = 6.1, by-0.0
rect3 = patches.FancyBboxPatch((mx, my), 1.3, 0.7, boxstyle="round,pad=0.05",
                                 facecolor='#ECFDF5', edgecolor=c_green, linewidth=1.2)
ax.add_patch(rect3)
ax.text(mx+0.65, my+0.42, 'Morton key', fontsize=7, ha='center', color='#065F46', fontweight='bold')
ax.text(mx+0.65, my+0.15, r'$M$ (64-bit)', fontsize=7, ha='center', color='#065F46')

# Arrow 3 (down)
ax.annotate('', xy=(mx+0.65, my-0.25), xytext=(mx+0.65, my-0.05),
            arrowprops=dict(arrowstyle='->', color=c_dark, lw=1.2))

# Step 4: Hash function
hx, hy = 5.8, 0.8
rect4 = patches.FancyBboxPatch((hx, hy), 1.8, 0.65, boxstyle="round,pad=0.05",
                                 facecolor='#EDE9FE', edgecolor='#7C3AED', linewidth=1.2)
ax.add_patch(rect4)
ax.text(hx+0.9, hy+0.4, r'$h(M)$', fontsize=8, ha='center', color='#5B21B6', fontweight='bold')
ax.text(hx+0.9, hy+0.12, 'multiplicative hash', fontsize=6, ha='center', color='#7C3AED')

# Arrow 4 (left to hash table)
ax.annotate('', xy=(hx-0.05, hy+0.32), xytext=(hx-0.5, hy+0.32),
            arrowprops=dict(arrowstyle='<-', color=c_dark, lw=1.2))

# Step 5: Hash table visualization
tx, ty = 2.3, 0.4
table_w, table_h = 3.0, 1.0
rect5 = patches.FancyBboxPatch((tx, ty), table_w, table_h, boxstyle="round,pad=0.05",
                                 facecolor='white', edgecolor=c_blue, linewidth=1.2)
ax.add_patch(rect5)
ax.text(tx+table_w/2, ty+table_h+0.1, 'Per-level hash table', fontsize=7, 
        ha='center', color=c_dark, fontweight='bold')

# Draw table slots
n_slots = 12
slot_w = 0.2
slot_h = 0.35
start_x = tx + 0.15
start_y = ty + 0.35
filled = [0, 2, 3, 5, 7, 9, 10]  # filled slots
for i in range(n_slots):
    color = c_blue if i in filled else '#F3F4F6'
    edge = c_dark if i in filled else c_gray
    slot = patches.Rectangle((start_x + i*(slot_w+0.03), start_y), slot_w, slot_h,
                               facecolor=color, edgecolor=edge, linewidth=0.5, alpha=0.7)
    ax.add_patch(slot)

ax.text(tx+table_w/2, ty+0.15, r'capacity = $2^k$, load factor $< 0.7$', fontsize=5.5, 
        ha='center', color=c_gray)

# Arrow from hash table to output
ax.annotate('', xy=(tx+table_w+0.4, ty+0.65), xytext=(tx+table_w+0.05, ty+0.65),
            arrowprops=dict(arrowstyle='->', color=c_dark, lw=1.2))
# Output: igrid
ox, oy = tx+table_w+0.4, ty+0.35
rect6 = patches.FancyBboxPatch((ox, oy), 0.9, 0.6, boxstyle="round,pad=0.05",
                                 facecolor='#FEE2E2', edgecolor='#DC2626', linewidth=1.2)
ax.add_patch(rect6)
ax.text(ox+0.45, oy+0.35, 'igrid', fontsize=8, ha='center', fontfamily='monospace', 
        color='#991B1B', fontweight='bold')

# Memory annotation
ax.text(tx+table_w/2, ty-0.15, r'$\approx 17\,N_{\rm grids}$ bytes (typically $<$ 50 MB)', 
        fontsize=6.5, ha='center', color=c_green, fontweight='bold')

# Title
ax.text(5.0, 3.85, 'Morton key hash table replaces nbor array', fontsize=9, 
        ha='center', fontweight='bold', color=c_dark)

plt.tight_layout()
plt.savefig('/gpfs/kjhan/Hydro/CUBE_HR5/code_cube/misc/fig_morton_hash.pdf', 
            dpi=300, bbox_inches='tight', pad_inches=0.05)
plt.close()
print("fig_morton_hash.pdf generated successfully")
