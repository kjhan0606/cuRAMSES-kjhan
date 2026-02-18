#!/usr/bin/env python3
"""
Progressive K-Section visualization: 3 -> 3x2 -> 3x2x2
Three isometric panels showing step-by-step domain decomposition.
Single-pass I/O over the TSC density file for all stages.
"""

import struct, math, colorsys
import numpy as np
from scipy.ndimage import zoom
from PIL import Image, ImageDraw, ImageFont


# ─── Data I/O ───────────────────────────────────────────────────────────────

def make_uniform_boxes(nx, ny, nz):
    """Create uniform subbox coordinates in [0,1]^3."""
    boxes = []
    for ix in range(nx):
        for iy in range(ny):
            for iz in range(nz):
                boxes.append((
                    ix / nx, (ix + 1) / nx,
                    iy / ny, (iy + 1) / ny,
                    iz / nz, (iz + 1) / nz,
                ))
    return boxes


def make_progressive_ksection_boxes():
    """
    Build progressive k-section boxes with non-uniform wall positions.
    Volume variation ~10-30% between domains at each split level.

    Stage 1: k=3 along x  (walls shifted from uniform 1/3, 2/3)
    Stage 2: k=2 along y  (per x-strip, wall shifted from 0.5)
    Stage 3: k=2 along z  (per xy-cell, wall shifted from 0.5)
    """
    # Stage 1: split x into 3 non-uniform strips
    # Uniform: [0, 0.333, 0.667, 1.0] → each 33.3%
    # Non-uniform: 27%, 34%, 39%
    x_walls = [0.0, 0.27, 0.61, 1.0]
    stage1 = []
    for i in range(3):
        stage1.append((x_walls[i], x_walls[i+1], 0.0, 1.0, 0.0, 1.0))

    # Stage 2: split each x-strip along y into 2
    # Different y-split for each x-strip (shifted ±5-10% from 0.5)
    y_splits = [0.44, 0.56, 0.41]
    stage2 = []
    for i, box in enumerate(stage1):
        xlo, xhi, ylo, yhi, zlo, zhi = box
        ymid = ylo + y_splits[i] * (yhi - ylo)
        stage2.append((xlo, xhi, ylo, ymid, zlo, zhi))
        stage2.append((xlo, xhi, ymid, yhi, zlo, zhi))

    # Stage 3: split each xy-cell along z into 2
    # Different z-split for each of the 6 xy-cells
    z_splits = [0.47, 0.54, 0.43, 0.57, 0.46, 0.53]
    stage3 = []
    for i, box in enumerate(stage2):
        xlo, xhi, ylo, yhi, zlo, zhi = box
        zmid = zlo + z_splits[i] * (zhi - zlo)
        stage3.append((xlo, xhi, ylo, yhi, zlo, zmid))
        stage3.append((xlo, xhi, ylo, yhi, zmid, zhi))

    return [stage1, stage2, stage3]


def precompute_all_projections(fname, all_boxes, target_ng=1024):
    """
    Single-pass projection computation for all box configurations combined.
    Reads TSC density file plane-by-plane, computes projections for all boxes.
    """
    with open(fname, 'rb') as f:
        nx, ny, nz = struct.unpack('iii', f.read(12))
    stride = nx // target_ng
    ng = target_ng
    print(f"  TSC file: {nx}x{ny}x{nz}, stride={stride}, effective {ng}^3")

    plane_bytes = nx * ny * 4
    header_bytes = 12

    box_ranges = []
    for box in all_boxes:
        xmin, xmax, ymin, ymax, zmin, zmax = box
        i0, i1 = max(0, int(xmin * ng)), min(ng, int(xmax * ng))
        j0, j1 = max(0, int(ymin * ng)), min(ng, int(ymax * ng))
        k0, k1 = max(0, int(zmin * ng)), min(ng, int(zmax * ng))
        box_ranges.append((i0, i1, j0, j1, k0, k1))

    projs = {}
    for bi, (i0, i1, j0, j1, k0, k1) in enumerate(box_ranges):
        projs[(bi, 0)] = np.zeros((j1 - j0, k1 - k0), dtype=np.float64)
        projs[(bi, 1)] = np.zeros((i1 - i0, k1 - k0), dtype=np.float64)
        projs[(bi, 2)] = np.zeros((i1 - i0, j1 - j0), dtype=np.float64)

    with open(fname, 'rb') as f:
        for kds in range(ng):
            k_file = kds * stride
            f.seek(header_bytes + k_file * plane_bytes)
            raw = f.read(plane_bytes)
            plane_full = np.frombuffer(raw, dtype=np.float32).reshape(ny, nx)
            plane_ds = plane_full[::stride, ::stride].T[:ng, :ng]

            for bi, (i0, i1, j0, j1, k0, k1) in enumerate(box_ranges):
                if kds < k0 or kds >= k1:
                    continue
                sub = plane_ds[i0:i1, j0:j1]
                klocal = kds - k0
                projs[(bi, 0)][:, klocal] += sub.sum(axis=0)
                projs[(bi, 1)][:, klocal] += sub.sum(axis=1)
                projs[(bi, 2)] += sub

            if kds % 200 == 0:
                print(f"    Plane {kds}/{ng}")

    print(f"    All {ng} planes done")
    return projs, box_ranges, ng


# ─── Coloring ──────────────────────────────────────────────────────────────

def assign_progressive_colors(nx, ny, nz):
    """
    HSV-based colors for progressive k-section.
    Hue varies with x-index (3 families: red, green, blue).
    Saturation varies with y-index, value with z-index.
    """
    colors = []
    base_hues = [0.0, 0.30, 0.60]

    for ix in range(nx):
        hue = base_hues[ix] if nx <= 3 else ix / nx
        for iy in range(ny):
            sat = 0.55 + 0.40 * iy / max(ny - 1, 1) if ny > 1 else 0.75
            for iz in range(nz):
                val = 0.65 + 0.30 * iz / max(nz - 1, 1) if nz > 1 else 0.85
                r, g, b = colorsys.hsv_to_rgb(hue, sat, val)
                colors.append((r, g, b))
    return colors


# ─── Rendering ─────────────────────────────────────────────────────────────

def make_subbox_rgba(density_2d, base_rgb, vmin, vmax, face_res=512):
    """Convert projected density to RGBA image with position-based coloring."""
    h_in, w_in = density_2d.shape
    if h_in == 0 or w_in == 0:
        return Image.new('RGBA', (face_res, face_res), (0, 0, 0, 255))

    resized = zoom(density_2d.astype(np.float64),
                   (face_res / h_in, face_res / w_in), order=1)
    resized = resized[:face_res, :face_res]

    log_d = np.log10(1.0 + resized)
    t = np.clip((log_d - vmin) / max(vmax - vmin, 1e-30), 0.0, 1.0)
    t = t ** 0.6

    brightness = 0.03 + 0.97 * t
    r = base_rgb[0] * brightness
    g = base_rgb[1] * brightness
    b = base_rgb[2] * brightness

    rgb = np.stack([r, g, b], axis=-1)
    rgba = np.zeros((*rgb.shape[:2], 4), dtype=np.uint8)
    rgba[:, :, :3] = np.clip(rgb * 255, 0, 255).astype(np.uint8)
    rgba[:, :, 3] = 255
    return Image.fromarray(rgba, 'RGBA')


def iso_project(x, y, z, scale, cx, cy):
    """Isometric projection: 3D -> 2D screen coordinates."""
    cos30 = math.cos(math.radians(30))
    sin30 = math.sin(math.radians(30))
    sx = (x - y) * cos30 * scale + cx
    sy = -(x + y) * sin30 * scale - z * scale + cy
    return sx, sy


def find_perspective_coeffs(src_pts, dst_pts):
    """Find perspective transform coefficients for PIL."""
    matrix = []
    for (xs, ys), (xd, yd) in zip(src_pts, dst_pts):
        matrix.append([xs, ys, 1, 0, 0, 0, -xd * xs, -xd * ys])
        matrix.append([0, 0, 0, xs, ys, 1, -yd * xs, -yd * ys])
    A = np.array(matrix, dtype=np.float64)
    b = np.array([c for pt in dst_pts for c in pt], dtype=np.float64)
    coeffs = np.linalg.lstsq(A, b, rcond=None)[0]
    return tuple(coeffs)


def warp_face_to_canvas(face_img, quad_corners, canvas_size):
    """Warp face image onto canvas at given quad corners."""
    w, h = face_img.size
    src_pts = [(0, 0), (w, 0), (w, h), (0, h)]
    dst_pts = list(quad_corners)

    xs = [p[0] for p in dst_pts]
    ys = [p[1] for p in dst_pts]
    x0, x1 = int(math.floor(min(xs))) - 2, int(math.ceil(max(xs))) + 2
    y0, y1 = int(math.floor(min(ys))) - 2, int(math.ceil(max(ys))) + 2
    x0, y0 = max(0, x0), max(0, y0)
    x1, y1 = min(canvas_size[0], x1), min(canvas_size[1], y1)

    cw, ch = x1 - x0, y1 - y0
    if cw <= 0 or ch <= 0:
        return Image.new('RGBA', canvas_size, (0, 0, 0, 0))

    shifted_dst = [(px - x0, py - y0) for px, py in dst_pts]
    coeffs_shifted = find_perspective_coeffs(shifted_dst, src_pts)

    warped = face_img.transform((cw, ch), Image.PERSPECTIVE, coeffs_shifted,
                                Image.BILINEAR)

    mask = Image.new('L', (cw, ch), 0)
    mask_draw = ImageDraw.Draw(mask)
    mask_draw.polygon([(int(round(px)), int(round(py))) for px, py in shifted_dst],
                      fill=255)
    warped.putalpha(mask)

    result = Image.new('RGBA', canvas_size, (0, 0, 0, 0))
    result.paste(warped, (x0, y0))
    return result


def render_panel(projs, boxes, colors, vmin, vmax,
                 canvas_w=2000, canvas_h=2000, gap=0.015, face_res=512):
    """Render one isometric panel for a single decomposition stage."""
    canvas = Image.new('RGBA', (canvas_w, canvas_h), (255, 255, 255, 255))

    margin = 0.03
    cos30 = math.cos(math.radians(30))
    sin30 = math.sin(math.radians(30))
    total_x = 2 * cos30
    total_y = 1 + 2 * sin30
    scale = min(canvas_w * (1 - 2 * margin) / total_x,
                canvas_h * (1 - 2 * margin) / total_y)

    scx, scy = iso_project(0.5, 0.5, 0.5, scale, 0, 0)
    cx = canvas_w / 2 - scx
    cy = canvas_h / 2 - scy

    ncpu = len(boxes)

    def box_depth(bi):
        b = boxes[bi]
        return (b[0]+b[1])/2 + (b[2]+b[3])/2 - (b[4]+b[5])/2

    order = sorted(range(ncpu), key=lambda i: box_depth(i), reverse=True)

    edge_color = (255, 255, 255, 255)
    edge_width = 2

    for bi in order:
        box = boxes[bi]
        bx0, bx1, by0, by1, bz0, bz1 = box

        mx, my, mz = (bx0+bx1)/2, (by0+by1)/2, (bz0+bz1)/2
        hx = (bx1 - bx0) / 2 - gap
        hy = (by1 - by0) / 2 - gap
        hz = (bz1 - bz0) / 2 - gap
        x0, x1 = mx - hx, mx + hx
        y0, y1 = my - hy, my + hy
        z0, z1 = mz - hz, mz + hz

        base_rgb = colors[bi]

        # Left face: -X
        proj = projs[(bi, 0)]
        fi = make_subbox_rgba(proj.T, base_rgb, vmin, vmax, face_res)
        fi = fi.transpose(Image.FLIP_TOP_BOTTOM)
        quad_left = [
            iso_project(x0, y0, z1, scale, cx, cy),
            iso_project(x0, y1, z1, scale, cx, cy),
            iso_project(x0, y1, z0, scale, cx, cy),
            iso_project(x0, y0, z0, scale, cx, cy),
        ]
        warped = warp_face_to_canvas(fi, quad_left, (canvas_w, canvas_h))
        canvas = Image.alpha_composite(canvas, warped)

        # Right face: -Y
        proj = projs[(bi, 1)]
        fi = make_subbox_rgba(proj.T, base_rgb, vmin, vmax, face_res)
        fi = fi.transpose(Image.FLIP_TOP_BOTTOM)
        quad_right = [
            iso_project(x0, y0, z1, scale, cx, cy),
            iso_project(x1, y0, z1, scale, cx, cy),
            iso_project(x1, y0, z0, scale, cx, cy),
            iso_project(x0, y0, z0, scale, cx, cy),
        ]
        warped = warp_face_to_canvas(fi, quad_right, (canvas_w, canvas_h))
        canvas = Image.alpha_composite(canvas, warped)

        # Top face: +Z
        proj = projs[(bi, 2)]
        fi = make_subbox_rgba(proj.T, base_rgb, vmin, vmax, face_res)
        quad_top = [
            iso_project(x0, y0, z1, scale, cx, cy),
            iso_project(x1, y0, z1, scale, cx, cy),
            iso_project(x1, y1, z1, scale, cx, cy),
            iso_project(x0, y1, z1, scale, cx, cy),
        ]
        warped = warp_face_to_canvas(fi, quad_top, (canvas_w, canvas_h))
        canvas = Image.alpha_composite(canvas, warped)

        # Wireframe edges
        draw = ImageDraw.Draw(canvas)
        for i in range(4):
            draw.line([quad_top[i], quad_top[(i+1)%4]], fill=edge_color, width=edge_width)
        draw.line([quad_left[1], quad_left[2]], fill=edge_color, width=edge_width)
        draw.line([quad_left[2], quad_left[3]], fill=edge_color, width=edge_width)
        draw.line([quad_right[1], quad_right[2]], fill=edge_color, width=edge_width)
        draw.line([quad_right[2], quad_right[3]], fill=edge_color, width=edge_width)
        draw.line([quad_left[3], quad_left[0]], fill=edge_color, width=edge_width)

    return canvas


def draw_arrow(draw, x0, y0, x1, y1, color=(120, 120, 120, 255), width=4, head_len=20):
    """Draw an arrow from (x0,y0) to (x1,y1) with a triangular head."""
    draw.line([(x0, y0), (x1, y1)], fill=color, width=width)
    angle = math.atan2(y1 - y0, x1 - x0)
    lx = x1 - head_len * math.cos(angle - 0.4)
    ly = y1 - head_len * math.sin(angle - 0.4)
    rx = x1 - head_len * math.cos(angle + 0.4)
    ry = y1 - head_len * math.sin(angle + 0.4)
    draw.polygon([(x1, y1), (int(lx), int(ly)), (int(rx), int(ry))], fill=color)


def draw_ksection_tree(draw, stage_idx, domain_colors,
                        x0, y0, width, height, font_small=None):
    """
    Draw a k-section tree diagram for a given decomposition stage.
    stage 0: single root node (undivided)
    stage 1: root → 3 leaves
    stage 2: root → 3 → 6 leaves
    stage 3: root → 3 → 6 → 12 leaves
    """
    branching = [3, 2, 2]

    # Build tree level counts
    levels = [1]
    for s in range(stage_idx):
        levels.append(levels[-1] * branching[s])

    n_levels = len(levels)
    max_leaves = levels[-1]

    node_r = max(8, min(15, int(width / (max(max_leaves, 1) * 5))))

    # Compute node positions
    positions = []
    for l in range(n_levels):
        n = levels[l]
        if n_levels == 1:
            y = y0 + height // 2
        else:
            y = y0 + node_r + int(l / (n_levels - 1) * (height - 2 * node_r))
        pts = []
        for i in range(n):
            x = x0 + int((i + 0.5) / n * width)
            pts.append((x, y))
        positions.append(pts)

    # Convert domain colors to 0-255
    leaf_colors_255 = [tuple(int(v * 255) for v in c[:3]) for c in domain_colors]

    # Compute node colors: leaves = domain colors, internal = avg of children
    node_colors = [None] * n_levels
    node_colors[-1] = leaf_colors_255[:levels[-1]]

    for l in range(n_levels - 2, -1, -1):
        k = branching[l] if l < len(branching) else 1
        colors = []
        for i in range(levels[l]):
            child_start = i * k
            child_end = min(child_start + k, levels[l + 1])
            avg = [0, 0, 0]
            count = 0
            for ci in range(child_start, child_end):
                if ci < len(node_colors[l + 1]):
                    for ch in range(3):
                        avg[ch] += node_colors[l + 1][ci][ch]
                    count += 1
            avg = tuple(int(a / count) for a in avg) if count > 0 else (160, 160, 160)
            colors.append(avg)
        node_colors[l] = colors

    # Root = gray
    node_colors[0] = [(180, 180, 180)]

    # Draw edges
    for l in range(1, n_levels):
        k = branching[l - 1]
        for i in range(levels[l]):
            parent_i = i // k
            px, py = positions[l - 1][parent_i]
            cx, cy = positions[l][i]
            draw.line([(px, py), (cx, cy)], fill=(160, 160, 160, 255), width=2)

    # Draw nodes
    for l in range(n_levels):
        for i in range(levels[l]):
            x, y = positions[l][i]
            r, g, b = node_colors[l][i][:3]
            draw.ellipse([x - node_r, y - node_r, x + node_r, y + node_r],
                         fill=(r, g, b, 255),
                         outline=(60, 60, 60, 255), width=1)

    # Branching factor labels on the left side
    if font_small and n_levels > 1:
        for l in range(1, n_levels):
            k = branching[l - 1]
            label = f"\u00d7{k}"
            mid_y = (positions[l - 1][0][1] + positions[l][0][1]) // 2 - 8
            draw.text((x0 - 40, mid_y), label, fill=(120, 120, 120, 255),
                      font=font_small)


# ─── Main ────────────────────────────────────────────────────────────────────

def main():
    tsc_file = '/gpfs/kjhan/Hydro/CUBE_HR5/code_cube/misc/TSC_Density.02981.dat'
    print(f"Using TSC density: {tsc_file}")

    target_ng = 1024

    # 4 stages: undivided + 3 progressive k-section stages
    stage_labels = [
        (1, 1, 1, "(a) Full domain"),
        (3, 1, 1, "(b) k = 3"),
        (3, 2, 1, "(c) k = 3 \u00d7 2 = 6"),
        (3, 2, 2, "(d) k = 3 \u00d7 2 \u00d7 2 = 12"),
    ]

    # Build progressive k-section boxes
    ksec_stages = make_progressive_ksection_boxes()

    # Stage 0: single undivided box
    stage0_boxes = [(0.0, 1.0, 0.0, 1.0, 0.0, 1.0)]
    all_stages = [stage0_boxes] + ksec_stages

    # Collect all boxes from all stages
    all_boxes = []
    stage_info = []
    for si, (nx, ny, nz, name) in enumerate(stage_labels):
        offset = len(all_boxes)
        boxes = all_stages[si]
        stage_info.append((offset, len(boxes), nx, ny, nz, name))
        all_boxes.extend(boxes)
    total_boxes = len(all_boxes)
    print(f"Total boxes: {total_boxes} (1 + 3 + 6 + 12 = 22)")

    # Single-pass projection computation
    print("Computing projections...")
    projs, box_ranges, ng = precompute_all_projections(tsc_file, all_boxes, target_ng)

    # Global normalization
    print("Computing global normalization...")
    all_vals = []
    for bi in range(total_boxes):
        for axis in range(3):
            proj = projs[(bi, axis)]
            if proj.size > 0:
                all_vals.append(np.log10(1.0 + proj).ravel())
    all_vals = np.concatenate(all_vals)
    vmin = np.percentile(all_vals, 2)
    vmax = np.percentile(all_vals, 99)
    print(f"  vmin={vmin:.4f}, vmax={vmax:.4f}")
    del all_vals

    # Render each stage
    panel_w, panel_h = 1500, 1500
    cube_gap = 0.015
    face_res = 512
    tree_h = 400

    panels = []
    stage_colors = []
    for si, (offset, count, nx, ny, nz, name) in enumerate(stage_info):
        print(f"Rendering stage {si}: {name} ({count} boxes)...")
        stage_boxes = all_boxes[offset:offset+count]
        stage_projs = {}
        for bi_local in range(count):
            for axis in range(3):
                stage_projs[(bi_local, axis)] = projs[(offset + bi_local, axis)]

        if si == 0:
            # Undivided: neutral blue-gray tint
            colors = [(0.75, 0.80, 0.90)]
        else:
            colors = assign_progressive_colors(nx, ny, nz)

        panel = render_panel(stage_projs, stage_boxes, colors, vmin, vmax,
                             panel_w, panel_h, cube_gap, face_res)
        panels.append((panel, name))
        stage_colors.append(colors)

    # Combine into 2x2 layout with trees below each cube
    gap = 30          # narrow gap between cells
    label_h = 60
    cell_h = panel_h + label_h + tree_h
    margin_bottom = 20
    total_w = panel_w * 2 + gap
    total_h = cell_h * 2 + gap + margin_bottom

    final = Image.new('RGBA', (total_w, total_h), (255, 255, 255, 255))

    try:
        font_label = ImageFont.truetype(
            '/usr/share/fonts/urw-base35/NimbusSans-Regular.otf', 42)
        font_small = ImageFont.truetype(
            '/usr/share/fonts/urw-base35/NimbusSans-Regular.otf', 28)
    except OSError:
        font_label = ImageFont.load_default()
        font_small = font_label

    draw = ImageDraw.Draw(final)

    # 2x2 grid positions
    positions = [
        (0, 0),                        # (a) top-left
        (panel_w + gap, 0),            # (b) top-right
        (0, cell_h + gap),             # (c) bottom-left
        (panel_w + gap, cell_h + gap), # (d) bottom-right
    ]

    for si, ((panel, name), (x_off, y_off)) in enumerate(zip(panels, positions)):
        # Paste cube
        final.paste(panel, (x_off, y_off))

        # Label below cube
        bbox = draw.textbbox((0, 0), name, font=font_label)
        lw = bbox[2] - bbox[0]
        lx = x_off + (panel_w - lw) // 2
        ly = y_off + panel_h + 5
        draw.text((lx, ly), name, fill=(0, 0, 0, 255), font=font_label)

        # Tree below label
        tree_x0 = x_off + panel_w // 6 + 40
        tree_y0 = y_off + panel_h + label_h + 5
        tree_w = panel_w * 2 // 3
        draw_ksection_tree(draw, si, stage_colors[si],
                           tree_x0, tree_y0, tree_w, tree_h - 20,
                           font_small=font_small)

    # Short bold arrows showing progression: (a)→(b), (b)→(c), (c)→(d)
    arrow_color = (80, 80, 80, 255)
    aw = 7       # bold width
    ahl = 18     # arrowhead length
    aL = 50      # arrow shaft length

    # (a) → (b): short arrow in gap between top panels
    mid_x = panel_w + gap // 2
    mid_y = panel_h // 2
    draw_arrow(draw, mid_x - aL // 2, mid_y, mid_x + aL // 2, mid_y,
               color=arrow_color, width=aw, head_len=ahl)

    # (b) → (c): short diagonal arrow (top-right → bottom-left)
    mid_x = total_w // 2
    mid_y = cell_h + gap // 2
    dx, dy = aL // 2, aL // 3
    draw_arrow(draw, mid_x + dx, mid_y - dy, mid_x - dx, mid_y + dy,
               color=arrow_color, width=aw, head_len=ahl)

    # (c) → (d): short arrow in gap between bottom panels
    mid_x = panel_w + gap // 2
    mid_y = cell_h + gap + panel_h // 2
    draw_arrow(draw, mid_x - aL // 2, mid_y, mid_x + aL // 2, mid_y,
               color=arrow_color, width=aw, head_len=ahl)

    # Save
    outdir = '/gpfs/kjhan/Hydro/CUBE_HR5/code_cube/misc'

    out_png = f'{outdir}/ksection_progressive.png'
    final_rgb = final.convert('RGB')
    final_rgb.save(out_png, dpi=(300, 300))
    print(f"Saved: {out_png}")

    out_pdf = f'{outdir}/ksection_progressive.pdf'
    final_rgb.save(out_pdf, dpi=(300, 300))
    print(f"Saved: {out_pdf}")


if __name__ == '__main__':
    main()
