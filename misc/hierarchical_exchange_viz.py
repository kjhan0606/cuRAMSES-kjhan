#!/usr/bin/env python3
"""
Hierarchical exchange communication pattern for N_cpu=12 (=3x2x2).
Shows the progression of communication through tree levels and stages.
Node groupings shown as solid bordered rectangles with inter-group gaps.
"""

import colorsys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches


def assign_colors_12():
    """HSV-based colors for 12 ranks, matching ksection_progressive_viz.py."""
    base_hues = [0.0, 0.30, 0.60]
    nx, ny, nz = 3, 2, 2
    colors = []
    for ix in range(nx):
        hue = base_hues[ix]
        for iy in range(ny):
            sat = 0.55 + 0.40 * iy / max(ny - 1, 1)
            for iz in range(nz):
                val = 0.65 + 0.30 * iz / max(nz - 1, 1)
                r, g, b = colorsys.hsv_to_rgb(hue, sat, val)
                colors.append((r, g, b))
    return colors


# ── Position mapping with inter-group gaps ────────────────────────────
# Level-1 groups: {0-3}, {4-7}, {8-11} → add gap between groups
GAP = 0.5  # gap between level-1 groups
def rank_x(rank_idx):
    """Map 0-indexed rank to x-position with inter-group gaps."""
    group = rank_idx // 4
    within = rank_idx % 4
    return group * (4 + GAP) + within + 0.5


def draw_rank_boxes(ax, colors, box_w=0.7, box_h=0.5, y0=0.0):
    """Draw 12 colored rank boxes at y=y0 with group gaps."""
    for i in range(12):
        cx = rank_x(i)
        rect = mpatches.FancyBboxPatch(
            (cx - box_w / 2, y0 - box_h / 2), box_w, box_h,
            boxstyle="round,pad=0.06",
            facecolor=colors[i], edgecolor='#333333', linewidth=1.0,
            zorder=10)
        ax.add_patch(rect)
        r, g, b = colors[i]
        lum = 0.299 * r + 0.587 * g + 0.114 * b
        tc = 'white' if lum < 0.55 else 'black'
        ax.text(cx, y0, str(i + 1), ha='center', va='center',
                fontsize=9, fontweight='bold', color=tc, zorder=11)


def draw_node_groups(ax, groups, y0=0.0, box_h=0.5, pad=0.15,
                     edgecolor='#666666', lw=1.2, label_y=None):
    """Draw solid transparent rectangles around rank groups (tree nodes)."""
    for start, end, label in groups:
        x_left = rank_x(start) - 0.35 - pad
        x_right = rank_x(end) + 0.35 + pad
        y_bot = y0 - box_h / 2 - pad
        h = box_h + 2 * pad
        rect = mpatches.FancyBboxPatch(
            (x_left, y_bot), x_right - x_left, h,
            boxstyle="round,pad=0.05",
            facecolor='none', edgecolor=edgecolor,
            linewidth=lw, linestyle='-', zorder=3)
        ax.add_patch(rect)
        if label and label_y is not None:
            ax.text((rank_x(start) + rank_x(end)) / 2, label_y, label,
                    ha='center', va='top', fontsize=7, color=edgecolor)


def draw_arc_simple(ax, r1, r2, y_base, arc_color='#555555', lw=1.8):
    """Draw a clean arc with small arrow tips between rank indices."""
    x1c = rank_x(r1)
    x2c = rank_x(r2)
    cx = (x1c + x2c) / 2
    width = abs(x2c - x1c)
    height = np.sqrt(width) * 0.45 + 0.3

    theta = np.linspace(0.08, np.pi - 0.08, 80)
    xs = cx + (width / 2) * np.cos(theta)
    ys = y_base + height * np.sin(theta)
    ax.plot(xs, ys, color=arc_color, linewidth=lw, solid_capstyle='round',
            zorder=5)

    tri_size = 0.13
    for end_idx, xend in [(0, x1c), (-1, x2c)]:
        dx = xs[end_idx] - xs[end_idx + (1 if end_idx == 0 else -1)]
        dy = ys[end_idx] - ys[end_idx + (1 if end_idx == 0 else -1)]
        norm = np.sqrt(dx**2 + dy**2)
        if norm > 0:
            dx, dy = dx / norm * tri_size, dy / norm * tri_size
        px, py = -dy * 0.6, dx * 0.6
        tip_x, tip_y = xs[end_idx] + dx * 0.5, ys[end_idx] + dy * 0.5
        base_x, base_y = xs[end_idx] - dx * 0.8, ys[end_idx] - dy * 0.8
        tri = plt.Polygon(
            [(tip_x, tip_y),
             (base_x + px, base_y + py),
             (base_x - px, base_y - py)],
            closed=True, facecolor=arc_color, edgecolor='none', zorder=6)
        ax.add_patch(tri)


# ── Communication stages ──────────────────────────────────────────────

node_groups_lv1 = [
    (0, 3, 'child 1'), (4, 7, 'child 2'), (8, 11, 'child 3'),
]
node_groups_lv2 = [
    (0, 1, ''), (2, 3, ''), (4, 5, ''), (6, 7, ''), (8, 9, ''), (10, 11, ''),
]

stages = [
    {
        'label': 'Level 1, step 1\n'
                 '($k_1{=}3$, siblings 1$\\leftrightarrow$2)',
        'pairs': [
            {'pair': (0, 4), 'color': '#AA4444'},
            {'pair': (1, 5), 'color': '#AA4444'},
            {'pair': (2, 6), 'color': '#AA4444'},
            {'pair': (3, 7), 'color': '#AA4444'},
        ],
        'groups': node_groups_lv1,
    },
    {
        'label': 'Level 1, step 2\n'
                 '($k_1{=}3$, siblings $\\{$1,2$\\}$$\\leftrightarrow$3)',
        'pairs': [
            {'pair': (0, 8),  'color': '#BB3333'},
            {'pair': (1, 9),  'color': '#BB3333'},
            {'pair': (2, 10), 'color': '#BB3333'},
            {'pair': (3, 11), 'color': '#BB3333'},
            {'pair': (4, 8),  'color': '#DD8822'},
            {'pair': (5, 9),  'color': '#DD8822'},
            {'pair': (6, 10), 'color': '#DD8822'},
            {'pair': (7, 11), 'color': '#DD8822'},
        ],
        'legend': [
            ('#BB3333', 'child 1$\\leftrightarrow$3'),
            ('#DD8822', 'child 2$\\leftrightarrow$3'),
        ],
        'groups': node_groups_lv1,
    },
    {
        'label': 'Level 2\n($k_2{=}2$)',
        'pairs': [
            {'pair': (0, 2),  'color': '#44AA44'},
            {'pair': (1, 3),  'color': '#44AA44'},
            {'pair': (4, 6),  'color': '#44AA44'},
            {'pair': (5, 7),  'color': '#44AA44'},
            {'pair': (8, 10), 'color': '#44AA44'},
            {'pair': (9, 11), 'color': '#44AA44'},
        ],
        'groups': node_groups_lv2,
    },
    {
        'label': 'Level 3\n($k_3{=}2$)',
        'pairs': [
            {'pair': (0, 1),   'color': '#4444AA'},
            {'pair': (2, 3),   'color': '#4444AA'},
            {'pair': (4, 5),   'color': '#4444AA'},
            {'pair': (6, 7),   'color': '#4444AA'},
            {'pair': (8, 9),   'color': '#4444AA'},
            {'pair': (10, 11), 'color': '#4444AA'},
        ],
        'groups': [],
    },
]


# ── Main figure ───────────────────────────────────────────────────────

colors = assign_colors_12()

box_y = 0.0
arc_y = 0.30
x_max = rank_x(11) + 0.8  # rightmost rank + margin

ylim_ranges = []
for stage in stages:
    max_width = max(abs(rank_x(p['pair'][1]) - rank_x(p['pair'][0]))
                    for p in stage['pairs'])
    max_arc_h = np.sqrt(max_width) * 0.45 + 0.3
    ylo, yhi = -0.65, arc_y + max_arc_h + 0.3
    ylim_ranges.append((ylo, yhi, max_arc_h))

height_ratios = [yhi - ylo for ylo, yhi, _ in ylim_ranges]

fig, axes = plt.subplots(4, 1, figsize=(10, 11.5),
                         gridspec_kw={'hspace': 0.15,
                                      'height_ratios': height_ratios})

fig.suptitle(
    'Hierarchical exchange communication pattern for '
    '$N_{\\rm cpu}=12\\;(=3\\times 2\\times 2)$',
    fontsize=13, y=0.97, fontweight='bold')

for ax_idx, (ax, stage) in enumerate(zip(axes, stages)):
    ylo, yhi, max_arc_h = ylim_ranges[ax_idx]

    ax.set_xlim(-0.5, x_max + 0.5)
    ax.set_ylim(ylo, yhi)
    ax.set_aspect('equal')
    ax.axis('off')

    # Draw node grouping boxes (behind everything)
    if stage.get('groups'):
        draw_node_groups(ax, stage['groups'], y0=box_y, box_h=0.5,
                         pad=0.15, edgecolor='#666666', lw=1.2,
                         label_y=-0.52)

    # Draw rank boxes
    draw_rank_boxes(ax, colors, y0=box_y)

    # Draw arcs
    for p in stage['pairs']:
        a, b = p['pair']
        draw_arc_simple(ax, a, b, y_base=arc_y,
                        arc_color=p['color'], lw=1.6)

    # Stage label on the left
    ax.text(-0.4, (arc_y + max_arc_h) * 0.45, stage['label'],
            ha='right', va='center', fontsize=9,
            bbox=dict(boxstyle='round,pad=0.3', facecolor='#f0f0f0',
                      edgecolor='#cccccc', alpha=0.9))

    # Legend for multi-color stages
    if 'legend' in stage:
        legend_y = arc_y + max_arc_h + 0.05
        legend_x = x_max - 1.5
        for li, (lc, ltxt) in enumerate(stage['legend']):
            ax.plot([legend_x, legend_x + 0.5], [legend_y - li * 0.35] * 2,
                    color=lc, linewidth=2.5, solid_capstyle='round')
            ax.text(legend_x + 0.7, legend_y - li * 0.35, ltxt,
                    ha='left', va='center', fontsize=8)

    # Separator line below (except last)
    if ax_idx < 3:
        ax.axhline(y=-0.62, color='#dddddd', linewidth=0.5,
                   xmin=0.02, xmax=0.98)

fig.text(0.5, 0.02,
         '$N_{\\rm msg} = (k_1{-}1) + (k_2{-}1) + (k_3{-}1) = 2+1+1 = 4$'
         '  messages per rank per exchange',
         ha='center', va='center', fontsize=11,
         bbox=dict(boxstyle='round,pad=0.4', facecolor='#ffffdd',
                   edgecolor='#cccc88'))

plt.savefig('/gpfs/kjhan/Hydro/CUBE_HR5/code_cube/misc/hierarchical_exchange.pdf',
            bbox_inches='tight', dpi=300)
plt.savefig('/gpfs/kjhan/Hydro/CUBE_HR5/code_cube/misc/hierarchical_exchange.png',
            bbox_inches='tight', dpi=200)
print("Saved: hierarchical_exchange.pdf / .png")
