import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np

fig, ax = plt.subplots(1, 1, figsize=(3.5, 7.2))

# ---- Colors ----
C_RED   = '#DC2626'
C_BLUE  = '#1D4ED8'
C_RESTR = '#3B82F6'
C_PROL  = '#10B981'
C_SOLVE = '#7C3AED'
C_DARK  = '#1E293B'
C_EXCH  = '#F59E0B'
C_GRAY  = '#9CA3AF'
C_ELIM  = '#EF4444'

# ==================================================================
#  UPPER SECTION: V-Cycle Schematic
# ==================================================================
n_levels = 4
y_pos = [5.5, 4.5, 3.5, 2.5]
box_w, box_h = 0.62, 0.38

x_left  = [0.05, 0.43, 0.81, 1.19]
x_right = [2.58, 2.20, 1.82, 1.19]


def draw_rb(x, y):
    """Draw a smoothing box as a 2x2 checkerboard (R/B pattern)."""
    pad = 0.02
    # Outer border
    border = mpatches.FancyBboxPatch((x, y), box_w, box_h,
        boxstyle='round,pad=0.02', fc='white', ec=C_DARK, lw=0.8)
    ax.add_patch(border)
    # 2x2 checkerboard cells
    cw = (box_w - 2 * pad) / 2
    ch = (box_h - 2 * pad) / 2
    for row in range(2):
        for col in range(2):
            is_red = (row + col) % 2 == 0
            fc = '#FEE2E2' if is_red else '#DBEAFE'
            ec_c = C_RED if is_red else C_BLUE
            label = 'R' if is_red else 'B'
            cx = x + pad + col * cw
            cy = y + pad + row * ch
            cell = mpatches.Rectangle((cx, cy), cw, ch,
                fc=fc, ec=ec_c, lw=0.6)
            ax.add_patch(cell)
            ax.text(cx + cw / 2, cy + ch / 2, label, ha='center',
                    va='center', fontsize=5, color=ec_c, fontweight='bold')


def draw_exch(x, y, active=True, sz=0.045):
    """Draw a small exchange diamond (filled or crossed-out)."""
    if active:
        d = mpatches.RegularPolygon((x, y), 4, radius=sz,
            fc=C_EXCH, ec='#D97706', lw=0.5)
        ax.add_patch(d)
    else:
        d = mpatches.RegularPolygon((x, y), 4, radius=sz,
            fc='#FEF3C7', ec='#D1D5DB', lw=0.5, alpha=0.4)
        ax.add_patch(d)
        ax.plot([x - 0.03, x + 0.03], [y - 0.03, y + 0.03],
                color=C_ELIM, lw=0.8)
        ax.plot([x - 0.03, x + 0.03], [y + 0.03, y - 0.03],
                color=C_ELIM, lw=0.8)


# ---- Boxes ----
for i in range(n_levels):
    if i < n_levels - 1:
        draw_rb(x_left[i], y_pos[i])
        draw_rb(x_right[i], y_pos[i])
    else:
        r = mpatches.FancyBboxPatch((x_left[i], y_pos[i]), box_w, box_h,
            boxstyle='round,pad=0.04', fc='#EDE9FE', ec=C_SOLVE, lw=1.5)
        ax.add_patch(r)
        ax.text(x_left[i] + box_w / 2, y_pos[i] + box_h / 2, 'Solve',
                ha='center', va='center', fontsize=8, color=C_SOLVE,
                fontweight='bold')

# ---- Restrict arrows (left leg) ----
for i in range(n_levels - 1):
    ax.annotate('',
        xy=(x_left[i + 1] + box_w / 2, y_pos[i + 1] + box_h),
        xytext=(x_left[i] + box_w / 2, y_pos[i]),
        arrowprops=dict(arrowstyle='->', color=C_RESTR, lw=1.5))

# ---- Prolongate arrows (right leg) ----
for i in range(n_levels - 1):
    src_x = x_left[i + 1] if i == n_levels - 2 else x_right[i + 1]
    ax.annotate('',
        xy=(x_right[i] + box_w / 2, y_pos[i]),
        xytext=(src_x + box_w / 2, y_pos[i + 1] + box_h),
        arrowprops=dict(arrowstyle='->', color=C_PROL, lw=1.5))

# ---- Prolongation exchange diamonds (on arrows) ----
for i in range(n_levels - 1):
    src_x = x_left[i + 1] if i == n_levels - 2 else x_right[i + 1]
    mx = 0.5 * (src_x + box_w / 2 + x_right[i] + box_w / 2)
    my = 0.5 * (y_pos[i + 1] + box_h + y_pos[i])
    draw_exch(mx + 0.08, my, active=True, sz=0.04)

# ---- Level labels ----
for i in range(n_levels):
    ax.text(3.5, y_pos[i] + box_h / 2, f'Level {n_levels - i}',
            ha='left', va='center', fontsize=7, color=C_DARK)

# ---- Leg labels ----
ax.text(x_left[0] + box_w / 2, y_pos[0] + box_h + 0.15, 'Restrict',
        ha='center', fontsize=7, color=C_RESTR, fontweight='bold')
ax.text(x_right[0] + box_w / 2, y_pos[0] + box_h + 0.15, 'Prolongate',
        ha='center', fontsize=7, color=C_PROL, fontweight='bold')

# ---- Compact legend ----
ly = 2.05
ax.text(0.0, ly, 'R', fontsize=8, color=C_RED, fontweight='bold', va='center')
ax.text(0.15, ly, '= Red GS', fontsize=5.5, color=C_GRAY, va='center')
ax.text(1.05, ly, 'B', fontsize=8, color=C_BLUE, fontweight='bold', va='center')
ax.text(1.20, ly, '= Black GS', fontsize=5.5, color=C_GRAY, va='center')
d = mpatches.RegularPolygon((2.25, ly), 4, radius=0.05,
    fc=C_EXCH, ec='#D97706', lw=0.5)
ax.add_patch(d)
ax.text(2.40, ly, '= Exchange', fontsize=5.5, color=C_GRAY, va='center')

# ==================================================================
#  LOWER SECTION: Per-level exchange detail
# ==================================================================
ax.text(0.0, 1.65, 'Exchanges per smoothing step:', fontsize=6.5,
        color=C_DARK, fontweight='bold')

# Helper: draw small operation box
def draw_op(x, y, w, h, label, fc, ec, fs=5.5):
    r = mpatches.FancyBboxPatch((x, y), w, h,
        boxstyle='round,pad=0.015', fc=fc, ec=ec, lw=0.7)
    ax.add_patch(r)
    ax.text(x + w / 2, y + h / 2, label, ha='center', va='center',
            fontsize=fs, color=ec, fontweight='bold')

# Dimensions for small boxes
bw = 0.30
bh = 0.22
dsp = 0.14          # space for diamond between boxes
x0 = 0.62           # left margin for operation boxes

# ---- Original (4 exchanges per smooth) ----
oy = 1.15
ax.text(0.0, oy + bh / 2, 'Original:', fontsize=5.5, color=C_GRAY,
        va='center')

x = x0
draw_op(x, oy, bw, bh, 'R', '#FEE2E2', C_RED)
x += bw
draw_exch(x + dsp / 2, oy + bh / 2, active=True, sz=0.04)
x += dsp
draw_op(x, oy, bw, bh, 'B', '#DBEAFE', C_BLUE)
x += bw
draw_exch(x + dsp / 2, oy + bh / 2, active=True, sz=0.04)
x += dsp
draw_op(x, oy, bw + 0.02, bh, 'Res', '#F3F4F6', C_DARK)
x += bw + 0.02
draw_exch(x + dsp / 2, oy + bh / 2, active=True, sz=0.04)
x += dsp
draw_op(x, oy, bw + 0.06, bh, u'\u2016r\u2016\u2082', '#F3F4F6', C_DARK, fs=5)
x += bw + 0.06
draw_exch(x + dsp / 2, oy + bh / 2, active=True, sz=0.04)
x_end_orig = x + dsp

# Count annotation
ax.text(x_end_orig + 0.15, oy + bh / 2, u'\u00d72 +1 = 9',
        fontsize=5.5, color='#D97706', fontweight='bold', va='center')

# ---- Optimized (2 exchanges per smooth) ----
oy2 = 0.60
ax.text(0.0, oy2 + bh / 2, 'Optimized:', fontsize=5.5, color=C_GRAY,
        va='center')

x = x0
draw_op(x, oy2, bw, bh, 'R', '#FEE2E2', C_RED)
x += bw
draw_exch(x + dsp / 2, oy2 + bh / 2, active=False, sz=0.04)   # eliminated
x += dsp
draw_op(x, oy2, bw, bh, 'B', '#DBEAFE', C_BLUE)
x += bw
draw_exch(x + dsp / 2, oy2 + bh / 2, active=True, sz=0.04)    # kept
x += dsp
# Fused Res+norm
draw_op(x, oy2, bw + 0.28, bh, u'Res+\u2016r\u2016\u2082', '#F3F4F6',
        C_DARK, fs=5)
x += bw + 0.28
draw_exch(x + dsp / 2, oy2 + bh / 2, active=True, sz=0.04)    # kept (fused)
x_end_opt = x + dsp

# Count annotation
ax.text(x_end_opt + 0.15, oy2 + bh / 2, u'\u00d72 +1 = 5',
        fontsize=5.5, color=C_PROL, fontweight='bold', va='center')
ax.text(x_end_opt + 0.15, oy2 + bh / 2 - 0.15, u'(\u221244%)',
        fontsize=5, color=C_PROL, va='center')

# ---- Annotations: "merged" and "fused" ----
# "merged" between R and B in optimized (bracket pointing to eliminated diamond)
xm = x0 + bw + dsp / 2
ax.text(xm, oy2 - 0.06, 'merged', fontsize=4.2, color=C_ELIM,
        ha='center', va='top', style='italic')

# "fused" under combined Res+norm box
xf = x0 + bw + dsp + bw + dsp + (bw + 0.28) / 2
ax.text(xf, oy2 - 0.06, 'fused', fontsize=4.2, color=C_PROL,
        ha='center', va='top', style='italic')

# ---- Note about ×2 + 1 ----
ax.text(x0, 0.28, u'\u00d72 (pre + post-smooth) + 1 (prolongation)',
        fontsize=4.5, color=C_GRAY, style='italic')

# ---- Axis ----
ax.set_xlim(-0.15, 4.3)
ax.set_ylim(0.10, 6.3)
ax.axis('off')

plt.tight_layout()
plt.savefig('/gpfs/kjhan/Hydro/CUBE_HR5/code_cube/misc/fig_mg_vcycle.pdf',
            dpi=300, bbox_inches='tight', pad_inches=0.05)
plt.close()
print("fig_mg_vcycle.pdf generated successfully")
