#!/usr/bin/env python3
"""
Visualize 3x3x2 uniform domain decomposition with TSC density on surfaces.
On-the-fly projection: reads z-planes sequentially, accumulates projections
without storing the full 3D array. Supports 2048^3 effective resolution.
"""

import struct
import math
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize
from scipy.ndimage import zoom


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
                projs[(bi, 0)][:, klocal] += sub.sum(axis=0)  # sum over x → (ny_sub,)
                projs[(bi, 1)][:, klocal] += sub.sum(axis=1)  # sum over y → (nx_sub,)
                projs[(bi, 2)] += sub                          # accumulate z-proj

            if kds % 200 == 0:
                print(f"  Processed plane {kds}/{ng}")

    print(f"  All {ng} planes processed")
    return projs, box_ranges, ng


def sample_face_from_proj(projs, box_ranges, ng, box_idx, box, face, npts=96):
    """
    Get face meshgrid and projected density from pre-computed projections.
    face: 0=x_min, 1=x_max, 2=y_min, 3=y_max, 4=z_min, 5=z_max
    """
    xmin, xmax, ymin, ymax, zmin, zmax = box

    def resize(arr, ny_t, nx_t):
        if arr.size == 0:
            return np.zeros((ny_t, nx_t))
        s = arr.astype(float)
        out = zoom(s, (ny_t / max(s.shape[0], 1), nx_t / max(s.shape[1], 1)), order=1)
        return out[:ny_t, :nx_t]

    if face in (0, 1):  # x-face: show y-z projection (sum over x)
        proj = projs[(box_idx, 0)]  # (ny_sub, nz_sub)
        u = np.linspace(ymin, ymax, npts)
        v = np.linspace(zmin, zmax, npts)
        U, V = np.meshgrid(u, v)
        X = np.full_like(U, xmin if face == 0 else xmax)
        d2d = resize(proj, npts, npts)
        return X, U, V, d2d

    elif face in (2, 3):  # y-face: show x-z projection (sum over y)
        proj = projs[(box_idx, 1)]  # (nx_sub, nz_sub)
        u = np.linspace(xmin, xmax, npts)
        v = np.linspace(zmin, zmax, npts)
        U, V = np.meshgrid(u, v)
        Y = np.full_like(U, ymin if face == 2 else ymax)
        d2d = resize(proj, npts, npts)
        return U, Y, V, d2d

    else:  # z-face: show x-y projection (sum over z)
        proj = projs[(box_idx, 2)]  # (nx_sub, ny_sub)
        u = np.linspace(xmin, xmax, npts)
        v = np.linspace(ymin, ymax, npts)
        U, V = np.meshgrid(u, v)
        Z = np.full_like(U, zmin if face == 4 else zmax)
        d2d = resize(proj, npts, npts)
        return U, V, Z, d2d


def density_to_facecolors(d2d, cmap):
    """Map projected column density to RGBA using log10 scaling."""
    ny, nx = d2d.shape
    d_face = (d2d[:-1,:-1] + d2d[1:,:-1] + d2d[:-1,1:] + d2d[1:,1:]) / 4.0
    log_face = np.log10(1.0 + d_face)
    lmin, lmax = log_face.min(), log_face.max()
    if lmax > lmin:
        norm = Normalize(vmin=lmin, vmax=lmax)
    else:
        norm = Normalize(vmin=0, vmax=1)
    sm = ScalarMappable(norm=norm, cmap=cmap)
    rgba = sm.to_rgba(log_face)
    return rgba


def get_visible_faces(elev, azim):
    """Return 3 visible face indices for the given view angle."""
    cam_x = math.cos(math.radians(elev)) * math.cos(math.radians(azim))
    cam_y = math.cos(math.radians(elev)) * math.sin(math.radians(azim))
    cam_z = math.sin(math.radians(elev))
    faces = []
    faces.append(1 if cam_x > 0 else 0)
    faces.append(3 if cam_y > 0 else 2)
    faces.append(5 if cam_z > 0 else 4)
    return faces, (cam_x, cam_y, cam_z)


def draw_visible_edges(ax, box, visible_faces, lw=1.0, zorder=0):
    """Draw only edges adjacent to at least one visible face."""
    xmin, xmax, ymin, ymax, zmin, zmax = box
    c = [(xmin,ymin,zmin), (xmax,ymin,zmin), (xmin,ymax,zmin), (xmax,ymax,zmin),
         (xmin,ymin,zmax), (xmax,ymin,zmax), (xmin,ymax,zmax), (xmax,ymax,zmax)]
    edges = [
        (0,1, 2,4), (2,3, 3,4), (4,5, 2,5), (6,7, 3,5),
        (0,2, 0,4), (1,3, 1,4), (4,6, 0,5), (5,7, 1,5),
        (0,4, 0,2), (1,5, 1,2), (2,6, 0,3), (3,7, 1,3),
    ]
    vf = set(visible_faces)
    for ci, cj, f1, f2 in edges:
        if f1 in vf or f2 in vf:
            p0, p1 = c[ci], c[cj]
            ax.plot([p0[0],p1[0]], [p0[1],p1[1]], [p0[2],p1[2]],
                    'k-', lw=lw, zorder=zorder)


def render_scene(ax, projs, box_ranges, ng, boxes, cmaps, elev, azim, gap, npts, title_str=None):
    """Render all boxes with proper occlusion for given view angle."""
    ncpu = len(boxes)
    visible_faces, cam = get_visible_faces(elev, azim)
    print(f"  View elev={elev}, azim={azim}: visible faces={visible_faces}")

    def box_depth(box):
        cx = (box[0]+box[1])/2
        cy = (box[2]+box[3])/2
        cz = (box[4]+box[5])/2
        return cx*cam[0] + cy*cam[1] + cz*cam[2]

    order = sorted(range(ncpu), key=lambda i: box_depth(boxes[i]))

    for rank, i in enumerate(order):
        box = boxes[i]
        cx = (box[0]+box[1])/2; cy = (box[2]+box[3])/2; cz = (box[4]+box[5])/2
        sx = (box[1]-box[0])/2 - gap
        sy = (box[3]-box[2])/2 - gap
        sz = (box[5]-box[4])/2 - gap
        sbox = (cx-sx, cx+sx, cy-sy, cy+sy, cz-sz, cz+sz)

        cmap = cmaps[i % len(cmaps)]
        zorder_face = rank * 2
        zorder_edge = rank * 2 + 1

        for face in visible_faces:
            X, Y, Z, d2d = sample_face_from_proj(projs, box_ranges, ng, i, sbox, face, npts)
            fc = density_to_facecolors(d2d, cmap)
            ax.plot_surface(X, Y, Z, facecolors=fc, shade=False,
                            edgecolor='none', rstride=1, cstride=1,
                            zorder=zorder_face)

        draw_visible_edges(ax, sbox, visible_faces, lw=0.8, zorder=zorder_edge)

    ax.view_init(elev=elev, azim=azim)
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.set_zlim(0, 1)
    ax.set_axis_off()
    if title_str:
        ax.set_title(title_str, fontsize=14, pad=15)


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

    # 18 distinct sequential colormaps
    cmap_names = [
        'Reds', 'Blues', 'Greens', 'Purples', 'Oranges', 'YlOrBr',
        'RdPu', 'BuGn', 'PuBu', 'OrRd', 'GnBu', 'YlGn',
        'BuPu', 'YlOrRd', 'PuRd', 'Greys', 'bone_r', 'pink_r',
    ]
    cmaps = [plt.colormaps[name] for name in cmap_names]

    gap = 0.015
    npts = 96

    title = (f'Domain Decomposition {nx_dom}x{ny_dom}x{nz_dom} = {ncpu} subboxes\n'
             f'4096$^3$ TSC density  |  '
             r'log$_{10}$(1+$\Sigma$) projected column density')

    outdir = '/gpfs/kjhan/Hydro/CUBE_HR5/code_cube/doc/manual'

    print("Rendering view 1...")
    fig1 = plt.figure(figsize=(14, 12))
    ax1 = fig1.add_subplot(111, projection='3d', computed_zorder=False)
    render_scene(ax1, projs, box_ranges, ng, boxes, cmaps,
                 elev=30, azim=45, gap=gap, npts=npts, title_str=title)
    plt.tight_layout()
    outfile1 = f'{outdir}/ksection_domains_18ranks.png'
    fig1.savefig(outfile1, dpi=200, bbox_inches='tight')
    print(f"Saved: {outfile1}")
    plt.close(fig1)

    print("Rendering view 2...")
    fig2 = plt.figure(figsize=(14, 12))
    ax2 = fig2.add_subplot(111, projection='3d', computed_zorder=False)
    render_scene(ax2, projs, box_ranges, ng, boxes, cmaps,
                 elev=20, azim=55, gap=gap, npts=npts, title_str=title)
    plt.tight_layout()
    outfile2 = f'{outdir}/ksection_domains_18ranks_wide.png'
    fig2.savefig(outfile2, dpi=200, bbox_inches='tight')
    print(f"Saved: {outfile2}")
    plt.close(fig2)


if __name__ == '__main__':
    main()
