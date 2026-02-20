#!/usr/bin/env python3
"""
Hierarchical exchange communication pattern for N_cpu=12 (=3x2x2).
Shows the progression of communication through tree levels and stages.
All panels use identical rank positions for uniform layout.
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


# ── Unified position mapping ────────────────────────────────────────
GAP_LV1 = 0.6   # gap between level-1 groups (groups of 4)
GAP_LV2 = 0.35  # gap between level-2 sub-groups (pairs) within a group of 4


def rank_x(rank_idx):
    """Unified position for all panels: gaps between groups of 4 AND pairs."""
    group4 = rank_idx // 4
    pair_in_group = (rank_idx % 4) // 2
    within_pair = rank_idx % 2
    group4_width = 4 + GAP_LV2
    return (group4 * (group4_width + GAP_LV1)
            + pair_in_group * (2 + GAP_LV2)
            + within_pair + 0.5)


def draw_rank_boxes(ax, colors, box_w=0.7, box_h=0.5, y0=0.0):
    """Draw 12 colored rank boxes at y=y0."""
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
                     edgecolor='#666666', lw=1.2):
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


# ── Node groupings ────────────────────────────────────────────────────

# Level 1 groups (3 children of root, 4 ranks each)
groups_lv1 = [
    (0, 3, 'child 1'), (4, 7, 'child 2'), (8, 11, 'child 3'),
]
# Level 2 groups (6 children, 2 ranks each)
groups_lv2 = [
    (0, 1, ''), (2, 3, ''), (4, 5, ''), (6, 7, ''), (8, 9, ''), (10, 11, ''),
]

# ── Communication stages ──────────────────────────────────────────────

stages = [
    {
        'label': 'Level 1, step 1  ($k_1{=}3$, siblings 1$\\leftrightarrow$2)',
        'pairs': [
            {'pair': (0, 4), 'color': '#AA4444'},
            {'pair': (1, 5), 'color': '#AA4444'},
            {'pair': (2, 6), 'color': '#AA4444'},
            {'pair': (3, 7), 'color': '#AA4444'},
        ],
        'groups': groups_lv1,
    },
    {
        'label': 'Level 1, step 2  ($k_1{=}3$, siblings $\\{$1,2$\\}$$\\leftrightarrow$3)',
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
        'groups': groups_lv1,
    },
    {
        'label': 'Level 2  ($k_2{=}2$)',
        'pairs': [
            {'pair': (0, 2),  'color': '#44AA44'},
            {'pair': (1, 3),  'color': '#44AA44'},
            {'pair': (4, 6),  'color': '#44AA44'},
            {'pair': (5, 7),  'color': '#44AA44'},
            {'pair': (8, 10), 'color': '#44AA44'},
            {'pair': (9, 11), 'color': '#44AA44'},
        ],
        'groups': groups_lv2,
    },
    {
        'label': 'Level 3  ($k_3{=}2$)',
        'pairs': [
            {'pair': (0, 1),   'color': '#4444AA'},
            {'pair': (2, 3),   'color': '#4444AA'},
            {'pair': (4, 5),   'color': '#4444AA'},
            {'pair': (6, 7),   'color': '#4444AA'},
            {'pair': (8, 9),   'color': '#4444AA'},
            {'pair': (10, 11), 'color': '#4444AA'},
        ],
        'groups': groups_lv2,
    },
]


# ── Main figure ───────────────────────────────────────────────────────

colors = assign_colors_12()

box_y = 0.0
arc_y = 0.30

# Compute max arc height across ALL panels for uniform sizing
all_max_arc_h = 0
for stage in stages:
    for p in stage['pairs']:
        w = abs(rank_x(p['pair'][1]) - rank_x(p['pair'][0]))
        h = np.sqrt(w) * 0.45 + 0.3
        all_max_arc_h = max(all_max_arc_h, h)

# Uniform y-limits for all panels
ylo = -1.45   # room for label + legend below
yhi = arc_y + all_max_arc_h + 0.15
panel_height = yhi - ylo

# Uniform x-limits
x_min = rank_x(0) - 1.0
x_max = rank_x(11) + 1.0

fig, axes = plt.subplots(4, 1, figsize=(10, 10.0),
                         gridspec_kw={'hspace': 0.08,
                                      'height_ratios': [1, 1, 1, 1]})

fig.suptitle(
    'Hierarchical exchange communication pattern for '
    '$N_{\\rm cpu}=12\\;(=3\\times 2\\times 2)$',
    fontsize=13, y=0.97, fontweight='bold')

for ax_idx, (ax, stage) in enumerate(zip(axes, stages)):
    ax.set_xlim(x_min, x_max)
    ax.set_ylim(ylo, yhi)
    ax.set_aspect('equal')
    ax.axis('off')

    # Panel border (rectangular frame)
    # Use axes coordinates for the border
    border = mpatches.FancyBboxPatch(
        (x_min + 0.05, ylo + 0.05),
        (x_max - x_min) - 0.10,
        (yhi - ylo) - 0.10,
        boxstyle="round,pad=0.05",
        facecolor='none', edgecolor='#aaaaaa',
        linewidth=1.0, linestyle='-', zorder=1,
        transform=ax.transData)
    ax.add_patch(border)

    # Panel label (a), (b), (c), (d) at top-left
    panel_labels = ['(a)', '(b)', '(c)', '(d)']
    ax.text(x_min + 0.3, yhi - 0.15, panel_labels[ax_idx],
            ha='left', va='top', fontsize=11, fontweight='bold',
            zorder=12)

    # Draw node grouping boxes
    if stage.get('groups'):
        draw_node_groups(ax, stage['groups'], y0=box_y, box_h=0.5,
                         pad=0.15, edgecolor='#666666', lw=1.2)

    # Draw rank boxes
    draw_rank_boxes(ax, colors, y0=box_y)

    # Draw arcs
    for p in stage['pairs']:
        a, b = p['pair']
        draw_arc_simple(ax, a, b, y_base=arc_y,
                        arc_color=p['color'], lw=1.6)

    # Stage label at center-bottom (well below node grouping boxes at y=-0.40)
    label_y = -0.80
    ax.text((x_min + x_max) / 2, label_y, stage['label'],
            ha='center', va='center', fontsize=9,
            bbox=dict(boxstyle='round,pad=0.3', facecolor='#f0f0f0',
                      edgecolor='#cccccc', alpha=0.9),
            zorder=12)

    # Legend for multi-color stages (center-bottom, below label)
    if 'legend' in stage:
        legend_cx = (x_min + x_max) / 2
        legend_y = label_y - 0.35
        total_w = len(stage['legend']) * 3.5
        start_x = legend_cx - total_w / 2
        for li, (lc, ltxt) in enumerate(stage['legend']):
            lx = start_x + li * 3.5
            ax.plot([lx, lx + 0.5], [legend_y, legend_y],
                    color=lc, linewidth=2.5, solid_capstyle='round',
                    zorder=12)
            ax.text(lx + 0.7, legend_y, ltxt,
                    ha='left', va='center', fontsize=8, zorder=12)

fig.text(0.5, 0.015,
         '$N_{\\rm msg} = (k_1{-}1) + (k_2{-}1) + (k_3{-}1) = 2+1+1 = 4$'
         '  messages per rank per exchange',
         ha='center', va='bottom', fontsize=11,
         bbox=dict(boxstyle='round,pad=0.4', facecolor='#ffffdd',
                   edgecolor='#cccc88'))

plt.subplots_adjust(bottom=0.06)

plt.savefig('/gpfs/kjhan/Hydro/CUBE_HR5/code_cube/misc/hierarchical_exchange.pdf',
            bbox_inches='tight', dpi=300)
plt.savefig('/gpfs/kjhan/Hydro/CUBE_HR5/code_cube/misc/hierarchical_exchange.png',
            bbox_inches='tight', dpi=200)
print("Saved: hierarchical_exchange.pdf / .png")
