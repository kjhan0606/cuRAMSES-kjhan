#!/usr/bin/env python3
"""
Visualize 3x3x2 uniform domain decomposition with TSC density on surfaces.
PIL-based isometric projection for publication-quality rendering.
On-the-fly projection: reads z-planes sequentially, accumulates projections
without storing the full 3D array. Supports 2048^3 effective resolution.
"""

import struct
import math
import numpy as np
from scipy.ndimage import zoom
from PIL import Image, ImageDraw, ImageFont


# ─── Data I/O ───────────────────────────────────────────────────────────────

def make_uniform_boxes(nx, ny, nz):
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


def precompute_projections(fname, boxes, target_ng=2048):
    """
    Read TSC density file plane-by-plane with stride, compute all box projections.
    Returns dict: projections[(box_idx, proj_axis)] = 2D array
    proj_axis: 0=sum over x (yz-face), 1=sum over y (xz-face), 2=sum over z (xy-face)
    """
    with open(fname, 'rb') as f:
        nx, ny, nz = struct.unpack('iii', f.read(12))
    stride = nx // target_ng
    ng = target_ng
    print(f"TSC file: {nx}x{ny}x{nz}, stride={stride}, effective {ng}^3")

    plane_bytes = nx * ny * 4
    header_bytes = 12

    # Compute box pixel ranges on the downsampled grid
    box_ranges = []
    for box in boxes:
        xmin, xmax, ymin, ymax, zmin, zmax = box
        i0, i1 = max(0, int(xmin * ng)), min(ng, int(xmax * ng))
        j0, j1 = max(0, int(ymin * ng)), min(ng, int(ymax * ng))
        k0, k1 = max(0, int(zmin * ng)), min(ng, int(zmax * ng))
        box_ranges.append((i0, i1, j0, j1, k0, k1))

    # Allocate projection arrays
    projs = {}
    for bi, (i0, i1, j0, j1, k0, k1) in enumerate(box_ranges):
        projs[(bi, 0)] = np.zeros((j1 - j0, k1 - k0), dtype=np.float64)  # x-proj: (ny, nz)
        projs[(bi, 1)] = np.zeros((i1 - i0, k1 - k0), dtype=np.float64)  # y-proj: (nx, nz)
        projs[(bi, 2)] = np.zeros((i1 - i0, j1 - j0), dtype=np.float64)  # z-proj: (nx, ny)

    # Single pass through z-planes
    with open(fname, 'rb') as f:
        for kds in range(ng):
            k_file = kds * stride
            f.seek(header_bytes + k_file * plane_bytes)
            raw = f.read(plane_bytes)
            # Fortran column-major: den(x,y) at this z-plane
            plane_full = np.frombuffer(raw, dtype=np.float32).reshape(ny, nx)
            # Subsample and transpose: plane_ds[ix, iy]
            plane_ds = plane_full[::stride, ::stride].T[:ng, :ng]

            for bi, (i0, i1, j0, j1, k0, k1) in enumerate(box_ranges):
                if kds < k0 or kds >= k1:
                    continue
                sub = plane_ds[i0:i1, j0:j1]  # (nx_sub, ny_sub)
                klocal = kds - k0
                projs[(bi, 0)][:, klocal] += sub.sum(axis=0)  # sum over x -> (ny_sub,)
                projs[(bi, 1)][:, klocal] += sub.sum(axis=1)  # sum over y -> (nx_sub,)
                projs[(bi, 2)] += sub                          # accumulate z-proj

            if kds % 200 == 0:
                print(f"  Processed plane {kds}/{ng}")

    print(f"  All {ng} planes processed")
    return projs, box_ranges, ng


# ─── RGB Positional Coloring ─────────────────────────────────────────────────

def make_subbox_rgba(density_2d, base_rgb, vmin, vmax, face_res=512):
    """
    Convert projected column density to RGBA image.
    base_rgb: (R, G, B) tuple in [0,1] — subbox's position-based color.
    Low density -> near black, high density -> full base_rgb color.
    """
    h_in, w_in = density_2d.shape
    if h_in == 0 or w_in == 0:
        return Image.new('RGBA', (face_res, face_res), (0, 0, 0, 255))

    resized = zoom(density_2d.astype(np.float64),
                   (face_res / h_in, face_res / w_in), order=1)
    resized = resized[:face_res, :face_res]

    # Log scaling + normalize to [0, 1]
    log_d = np.log10(1.0 + resized)
    t = np.clip((log_d - vmin) / max(vmax - vmin, 1e-30), 0.0, 1.0)

    # Gamma stretch for better contrast
    t = t ** 0.6

    # RGB = base_color * brightness(t)
    brightness = 0.03 + 0.97 * t  # 0.03 at void, 1.0 at peak
    r = base_rgb[0] * brightness
    g = base_rgb[1] * brightness
    b = base_rgb[2] * brightness

    rgb = np.stack([r, g, b], axis=-1)
    rgba = np.zeros((*rgb.shape[:2], 4), dtype=np.uint8)
    rgba[:, :, :3] = np.clip(rgb * 255, 0, 255).astype(np.uint8)
    rgba[:, :, 3] = 255

    return Image.fromarray(rgba, 'RGBA')


def assign_colors(nx_dom, ny_dom, nz_dom):
    """
    Assign RGB base colors to subboxes: ix→R, iy→G, iz→B.
    Adjacent subboxes differ in only one channel → smooth spatial gradient.
    Returns list of (R, G, B) tuples in [0,1].
    """
    colors = []
    for ix in range(nx_dom):
        for iy in range(ny_dom):
            for iz in range(nz_dom):
                r = 0.45 + 0.55 * ix / max(nx_dom - 1, 1)
                g = 0.45 + 0.55 * iy / max(ny_dom - 1, 1)
                b = 0.45 + 0.55 * iz / max(nz_dom - 1, 1)
                colors.append((r, g, b))
    return colors


# ─── Isometric Projection ───────────────────────────────────────────────────

def iso_project(x, y, z, scale, cx, cy):
    """
    Isometric projection: 3D (x, y, z) -> 2D screen (sx, sy).
    +x → screen right-down, +y → screen left-down, +z → screen up.
    Camera at (-inf, -inf, +inf); visible faces: -X, -Y, +Z.
    """
    cos30 = math.cos(math.radians(30))
    sin30 = math.sin(math.radians(30))
    sx = (x - y) * cos30 * scale + cx
    sy = -(x + y) * sin30 * scale - z * scale + cy
    return sx, sy


def find_perspective_coeffs(src_pts, dst_pts):
    """
    Find perspective transform coefficients mapping src_pts -> dst_pts.
    Each is list of 4 (x, y) tuples.
    Returns 8 coefficients for PIL Image.transform(PERSPECTIVE).
    """
    matrix = []
    for (xs, ys), (xd, yd) in zip(src_pts, dst_pts):
        matrix.append([xs, ys, 1, 0, 0, 0, -xd * xs, -xd * ys])
        matrix.append([0, 0, 0, xs, ys, 1, -yd * xs, -yd * ys])
    A = np.array(matrix, dtype=np.float64)
    b = np.array([c for pt in dst_pts for c in pt], dtype=np.float64)
    coeffs = np.linalg.lstsq(A, b, rcond=None)[0]
    return tuple(coeffs)


def warp_face_to_canvas(face_img, quad_corners, canvas_size):
    """
    Warp a face image onto the canvas at the given quad corners.
    quad_corners: 4 (sx, sy) tuples for TL, TR, BR, BL of the face.
    Returns RGBA image of canvas_size with the warped face, masked to quad.
    """
    w, h = face_img.size
    src_pts = [(0, 0), (w, 0), (w, h), (0, h)]
    dst_pts = list(quad_corners)

    # Compute bounding box of destination quad
    xs = [p[0] for p in dst_pts]
    ys = [p[1] for p in dst_pts]
    x0, x1 = int(math.floor(min(xs))) - 2, int(math.ceil(max(xs))) + 2
    y0, y1 = int(math.floor(min(ys))) - 2, int(math.ceil(max(ys))) + 2

    # Clamp to canvas
    x0 = max(0, x0)
    y0 = max(0, y0)
    x1 = min(canvas_size[0], x1)
    y1 = min(canvas_size[1], y1)

    cw, ch = x1 - x0, y1 - y0
    if cw <= 0 or ch <= 0:
        return Image.new('RGBA', canvas_size, (0, 0, 0, 0))

    # Shift dst_pts so bounding box starts at (0,0)
    shifted_dst = [(px - x0, py - y0) for px, py in dst_pts]
    coeffs_shifted = find_perspective_coeffs(shifted_dst, src_pts)

    warped = face_img.transform((cw, ch), Image.PERSPECTIVE, coeffs_shifted,
                                Image.BILINEAR)

    # Create polygon mask for the quad (clip outside quad region)
    mask = Image.new('L', (cw, ch), 0)
    mask_draw = ImageDraw.Draw(mask)
    mask_draw.polygon([(int(round(px)), int(round(py))) for px, py in shifted_dst],
                      fill=255)
    warped.putalpha(mask)

    # Place on full canvas
    result = Image.new('RGBA', canvas_size, (0, 0, 0, 0))
    result.paste(warped, (x0, y0))
    return result


def draw_line_on_canvas(draw, p0, p1, color=(255, 255, 255, 255), width=2):
    """Draw a line on the canvas."""
    draw.line([p0, p1], fill=color, width=width)


def render_isometric(projs, box_ranges, ng, boxes, colors, vmin, vmax,
                     canvas_w=4000, canvas_h=4000, gap=0.015, face_res=512):
    """
    Render all subboxes with isometric projection using PIL.
    Painter's algorithm: farthest subbox first, closest subbox last.

    Isometric null-space analysis:
      M = [[cos30, -cos30, 0], [-sin30, -sin30, -1]]
      null(M) = (1, 1, -1)  →  camera at (-inf, -inf, +inf)
      Visible faces: -X (x=xmin), -Y (y=ymin), +Z (z=zmax)
      Depth along view axis: d = x + y - z  (larger = farther from camera)
    """
    canvas = Image.new('RGBA', (canvas_w, canvas_h), (255, 255, 255, 255))

    # Isometric scale and centering
    margin = 0.12
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
        """Depth along viewing direction (1,1,-1). Larger = farther from camera."""
        b = boxes[bi]
        cx_b = (b[0] + b[1]) / 2
        cy_b = (b[2] + b[3]) / 2
        cz_b = (b[4] + b[5]) / 2
        return cx_b + cy_b - cz_b

    # Painter's algorithm: draw farthest first (largest depth first = descending)
    order = sorted(range(ncpu), key=lambda i: box_depth(i), reverse=True)
    print(f"  Render order: depth {box_depth(order[0]):.3f}(far) → {box_depth(order[-1]):.3f}(near)")

    edge_color = (255, 255, 255, 255)
    edge_width = 2

    for bi in order:
        box = boxes[bi]
        bx0, bx1, by0, by1, bz0, bz1 = box

        # Apply gap
        mx, my, mz = (bx0+bx1)/2, (by0+by1)/2, (bz0+bz1)/2
        hx = (bx1 - bx0) / 2 - gap
        hy = (by1 - by0) / 2 - gap
        hz = (bz1 - bz0) / 2 - gap
        x0, x1 = mx - hx, mx + hx
        y0, y1 = my - hy, my + hy
        z0, z1 = mz - hz, mz + hz

        base_rgb = colors[bi]

        # --- Left face: -X (x=x0), shows yz-plane, projection axis 0 ---
        proj = projs[(bi, 0)]  # (ny_sub, nz_sub)
        fi = make_subbox_rgba(proj.T, base_rgb, vmin, vmax, face_res)
        # proj.T: rows=z(min→max), cols=y(min→max)
        # Need row0=z_max for TL → FLIP_TOP_BOTTOM
        fi = fi.transpose(Image.FLIP_TOP_BOTTOM)
        quad_left = [
            iso_project(x0, y0, z1, scale, cx, cy),  # TL (y_min, z_max)
            iso_project(x0, y1, z1, scale, cx, cy),  # TR (y_max, z_max)
            iso_project(x0, y1, z0, scale, cx, cy),  # BR (y_max, z_min)
            iso_project(x0, y0, z0, scale, cx, cy),  # BL (y_min, z_min)
        ]
        warped = warp_face_to_canvas(fi, quad_left, (canvas_w, canvas_h))
        canvas = Image.alpha_composite(canvas, warped)

        # --- Right face: -Y (y=y0), shows xz-plane, projection axis 1 ---
        proj = projs[(bi, 1)]  # (nx_sub, nz_sub)
        fi = make_subbox_rgba(proj.T, base_rgb, vmin, vmax, face_res)
        # proj.T: rows=z(min→max), cols=x(min→max)
        # Need row0=z_max → FLIP_TOP_BOTTOM
        fi = fi.transpose(Image.FLIP_TOP_BOTTOM)
        quad_right = [
            iso_project(x0, y0, z1, scale, cx, cy),  # TL (x_min, z_max)
            iso_project(x1, y0, z1, scale, cx, cy),  # TR (x_max, z_max)
            iso_project(x1, y0, z0, scale, cx, cy),  # BR (x_max, z_min)
            iso_project(x0, y0, z0, scale, cx, cy),  # BL (x_min, z_min)
        ]
        warped = warp_face_to_canvas(fi, quad_right, (canvas_w, canvas_h))
        canvas = Image.alpha_composite(canvas, warped)

        # --- Top face: +Z (z=z1), shows xy-plane, projection axis 2 ---
        proj = projs[(bi, 2)]  # (nx_sub, ny_sub)
        fi = make_subbox_rgba(proj.T, base_rgb, vmin, vmax, face_res)
        # proj.T: rows=y(min→max), cols=x(min→max)
        # TL=(x_min,y_min) matches row0,col0 → no flip needed
        quad_top = [
            iso_project(x0, y0, z1, scale, cx, cy),  # TL (x_min, y_min)
            iso_project(x1, y0, z1, scale, cx, cy),  # TR (x_max, y_min)
            iso_project(x1, y1, z1, scale, cx, cy),  # BR (x_max, y_max)
            iso_project(x0, y1, z1, scale, cx, cy),  # BL (x_min, y_max)
        ]
        warped = warp_face_to_canvas(fi, quad_top, (canvas_w, canvas_h))
        canvas = Image.alpha_composite(canvas, warped)

        # --- Wireframe edges (visible silhouette) ---
        draw = ImageDraw.Draw(canvas)
        # Top face diamond (4 edges)
        for i in range(4):
            draw.line([quad_top[i], quad_top[(i + 1) % 4]], fill=edge_color, width=edge_width)
        # Left face: far-left vertical + bottom
        draw.line([quad_left[1], quad_left[2]], fill=edge_color, width=edge_width)
        draw.line([quad_left[2], quad_left[3]], fill=edge_color, width=edge_width)
        # Right face: far-right vertical + bottom
        draw.line([quad_right[1], quad_right[2]], fill=edge_color, width=edge_width)
        draw.line([quad_right[2], quad_right[3]], fill=edge_color, width=edge_width)
        # Front vertical (shared edge between left and right)
        draw.line([quad_left[3], quad_left[0]], fill=edge_color, width=edge_width)

    return canvas


# ─── Main ────────────────────────────────────────────────────────────────────

def main():
    tsc_file = '/gpfs/kjhan/Hydro/CUBE_HR5/code_cube/misc/TSC_Density.02981.dat'
    print(f"Using TSC density: {tsc_file}")

    nx_dom, ny_dom, nz_dom = 3, 3, 2
    ncpu = nx_dom * ny_dom * nz_dom
    target_ng = 2048

    boxes = make_uniform_boxes(nx_dom, ny_dom, nz_dom)
    print(f"{ncpu} uniform subboxes: {nx_dom}x{ny_dom}x{nz_dom}")

    # Pre-compute all projections in single file pass
    print(f"Computing projections at {target_ng}^3 resolution...")
    projs, box_ranges, ng = precompute_projections(tsc_file, boxes, target_ng=target_ng)

    # Compute global normalization (log10(1+density))
    print("Computing global normalization...")
    all_vals = []
    for bi in range(ncpu):
        for axis in range(3):
            proj = projs[(bi, axis)]
            if proj.size > 0:
                log_p = np.log10(1.0 + proj)
                all_vals.append(log_p.ravel())
    all_vals = np.concatenate(all_vals)
    vmin = np.percentile(all_vals, 2)
    vmax = np.percentile(all_vals, 99)
    print(f"  Global range: vmin={vmin:.4f}, vmax={vmax:.4f}")
    del all_vals

    # Assign RGB colors: ix→R, iy→G, iz→B
    colors = assign_colors(nx_dom, ny_dom, nz_dom)
    print(f"  Colors: {['({:.2f},{:.2f},{:.2f})'.format(*c) for c in colors]}")

    # Render
    gap = 0.012
    face_res = 768
    canvas_w, canvas_h = 4000, 4000

    print("Rendering isometric view...")
    canvas = render_isometric(projs, box_ranges, ng, boxes, colors, vmin, vmax,
                              canvas_w=canvas_w, canvas_h=canvas_h,
                              gap=gap, face_res=face_res)

    # Add title
    try:
        font = ImageFont.truetype('/usr/share/fonts/urw-base35/NimbusRoman-Regular.otf', 72)
    except OSError:
        font = ImageFont.load_default()

    draw = ImageDraw.Draw(canvas)
    title_lines = [
        f"Domain Decomposition {nx_dom}x{ny_dom}x{nz_dom} = {ncpu} subboxes",
        f"4096^3 TSC density  |  log10(1+Sigma) projected column density",
    ]
    y_text = 60
    for line in title_lines:
        bbox = draw.textbbox((0, 0), line, font=font)
        tw = bbox[2] - bbox[0]
        x_text = (canvas_w - tw) // 2
        draw.text((x_text, y_text), line, fill=(0, 0, 0, 255), font=font)
        y_text += bbox[3] - bbox[1] + 15

    # Save
    outdir = '/gpfs/kjhan/Hydro/CUBE_HR5/code_cube/misc'

    out_png = f'{outdir}/domain_decomposition_improved.png'
    canvas_rgb = canvas.convert('RGB')
    canvas_rgb.save(out_png, dpi=(300, 300))
    print(f"Saved: {out_png}")

    out_pdf = f'{outdir}/domain_decomposition_improved.pdf'
    canvas_rgb.save(out_pdf, dpi=(300, 300))
    print(f"Saved: {out_pdf}")


if __name__ == '__main__':
    main()
