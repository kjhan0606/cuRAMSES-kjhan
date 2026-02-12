#!/usr/bin/env python3
"""
Visualize 3x2x2 uniform domain decomposition with density slice images on surfaces.
Each subbox uses a different sequential colormap.
High density = dark/saturated, low density = white/light.
Proper occlusion: only camera-facing faces drawn, back-to-front painter's algorithm.
"""

import struct
import glob
import os
import math
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize


def read_fortran_record(f):
    raw = f.read(4)
    if len(raw) < 4:
        return None
    nbytes = struct.unpack('<I', raw)[0]
    data = f.read(nbytes)
    f.read(4)
    return data


def read_particle_positions(output_dir):
    part_files = sorted(glob.glob(os.path.join(output_dir, 'part_*.out*')))
    if not part_files:
        raise FileNotFoundError(f"No particle files in {output_dir}")
    all_x, all_y, all_z = [], [], []
    for pf in part_files:
        with open(pf, 'rb') as f:
            read_fortran_record(f)
            read_fortran_record(f)
            npart_data = read_fortran_record(f)
            npart = struct.unpack('<I', npart_data)[0]
            for _ in range(5):
                read_fortran_record(f)
            if npart == 0:
                continue
            pos = []
            for idim in range(3):
                rec = read_fortran_record(f)
                pos.append(np.frombuffer(rec, dtype='<f8'))
            all_x.append(pos[0])
            all_y.append(pos[1])
            all_z.append(pos[2])
    x = np.concatenate(all_x)
    y = np.concatenate(all_y)
    z = np.concatenate(all_z)
    print(f"Read {len(x):,} particles from {len(part_files)} files")
    return x, y, z


def count_in_cell(x, y, z, ngrid=512):
    ix = np.clip((x * ngrid).astype(int), 0, ngrid - 1)
    iy = np.clip((y * ngrid).astype(int), 0, ngrid - 1)
    iz = np.clip((z * ngrid).astype(int), 0, ngrid - 1)
    density = np.zeros((ngrid, ngrid, ngrid), dtype=np.float64)
    np.add.at(density, (ix, iy, iz), 1)
    return density


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


def sample_face(density, box, face, npts=48):
    """
    Project density through the subbox onto a face (column density).
    face: 0=x_min, 1=x_max, 2=y_min, 3=y_max, 4=z_min, 5=z_max
    Returns meshgrids and 2D projected density.
    """
    from scipy.ndimage import zoom
    ng = density.shape[0]
    xmin, xmax, ymin, ymax, zmin, zmax = box

    i0, i1 = max(0, int(xmin * ng)), min(ng, int(xmax * ng))
    j0, j1 = max(0, int(ymin * ng)), min(ng, int(ymax * ng))
    k0, k1 = max(0, int(zmin * ng)), min(ng, int(zmax * ng))
    sub = density[i0:i1, j0:j1, k0:k1]  # (nx, ny, nz)

    def resize(arr, ny_t, nx_t):
        if arr.size == 0:
            return np.zeros((ny_t, nx_t))
        s = arr.astype(float)
        out = zoom(s, (ny_t / max(s.shape[0], 1), nx_t / max(s.shape[1], 1)), order=1)
        return out[:ny_t, :nx_t]

    if face in (0, 1):  # project along x -> (y, z)
        proj = sub.sum(axis=0)  # (ny, nz)
        u = np.linspace(ymin, ymax, npts)
        v = np.linspace(zmin, zmax, npts)
        U, V = np.meshgrid(u, v)
        X = np.full_like(U, xmin if face == 0 else xmax)
        d2d = resize(proj, npts, npts)
        return X, U, V, d2d

    elif face in (2, 3):  # project along y -> (x, z)
        proj = sub.sum(axis=1)  # (nx, nz)
        u = np.linspace(xmin, xmax, npts)
        v = np.linspace(zmin, zmax, npts)
        U, V = np.meshgrid(u, v)
        Y = np.full_like(U, ymin if face == 2 else ymax)
        d2d = resize(proj, npts, npts)
        return U, Y, V, d2d

    else:  # project along z -> (x, y)
        proj = sub.sum(axis=2)  # (nx, ny)
        u = np.linspace(xmin, xmax, npts)
        v = np.linspace(ymin, ymax, npts)
        U, V = np.meshgrid(u, v)
        Z = np.full_like(U, zmin if face == 4 else zmax)
        d2d = resize(proj, npts, npts)
        return U, V, Z, d2d


def density_to_facecolors(d2d, cmap):
    """
    Map projected column density to RGBA using log10 scaling.
    Returns (ny-1, nx-1, 4) array for plot_surface facecolors.
    """
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
    faces.append(1 if cam_x > 0 else 0)  # x_max or x_min
    faces.append(3 if cam_y > 0 else 2)  # y_max or y_min
    faces.append(5 if cam_z > 0 else 4)  # z_max or z_min
    return faces, (cam_x, cam_y, cam_z)


def draw_visible_edges(ax, box, visible_faces, lw=1.0, zorder=0):
    """Draw only edges adjacent to at least one visible face (9 of 12)."""
    xmin, xmax, ymin, ymax, zmin, zmax = box
    # Corner indices: bit0=x(0=min,1=max), bit1=y, bit2=z
    c = [(xmin,ymin,zmin), (xmax,ymin,zmin), (xmin,ymax,zmin), (xmax,ymax,zmin),
         (xmin,ymin,zmax), (xmax,ymin,zmax), (xmin,ymax,zmax), (xmax,ymax,zmax)]
    # (corner_i, corner_j, adjacent_face1, adjacent_face2)
    edges = [
        (0,1, 2,4), (2,3, 3,4), (4,5, 2,5), (6,7, 3,5),  # along x
        (0,2, 0,4), (1,3, 1,4), (4,6, 0,5), (5,7, 1,5),  # along y
        (0,4, 0,2), (1,5, 1,2), (2,6, 0,3), (3,7, 1,3),  # along z
    ]
    vf = set(visible_faces)
    for ci, cj, f1, f2 in edges:
        if f1 in vf or f2 in vf:
            p0, p1 = c[ci], c[cj]
            ax.plot([p0[0],p1[0]], [p0[1],p1[1]], [p0[2],p1[2]],
                    'k-', lw=lw, zorder=zorder)


def render_scene(ax, density, boxes, cmaps, elev, azim, gap, npts, npart):
    """Render all boxes with proper occlusion for given view angle."""
    ncpu = len(boxes)
    visible_faces, cam = get_visible_faces(elev, azim)
    print(f"  View elev={elev}, azim={azim}: visible faces={visible_faces}, cam={cam}")

    # Sort back-to-front using camera direction dot product
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

        # Draw only camera-facing faces (opaque)
        for face in visible_faces:
            X, Y, Z, d2d = sample_face(density, sbox, face, npts)
            fc = density_to_facecolors(d2d, cmap)
            ax.plot_surface(X, Y, Z, facecolors=fc, shade=False,
                            edgecolor='none', rstride=1, cstride=1,
                            zorder=zorder_face)

        # Draw visible edges (on top of same box's faces, under next box's faces)
        draw_visible_edges(ax, sbox, visible_faces, lw=0.8, zorder=zorder_edge)

    ax.view_init(elev=elev, azim=azim)
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.set_zlim(0, 1)
    ax.set_axis_off()
    ax.set_title(
        f'Domain Decomposition 3x2x2 = {ncpu} subboxes\n'
        f'{npart:,} particles  |  512$^3$ CIC  |  '
        r'log$_{10}$(1+$\Sigma$) projected column density',
        fontsize=14, pad=15
    )


def main():
    output_dir = '/gpfs/jaehyun/Darwin/Darwin2/01_test/output_00066'
    print(f"Using: {output_dir}")

    nx, ny, nz = 3, 2, 2
    ncpu = nx * ny * nz
    ngrid_cic = 512

    # Read particles
    x, y, z = read_particle_positions(output_dir)
    npart = len(x)

    # Build density field
    print(f"Building {ngrid_cic}^3 CIC density field...")
    density = count_in_cell(x, y, z, ngrid_cic)
    print(f"Density: min={density.min():.0f}, max={density.max():.0f}, mean={density.mean():.2f}")

    # Uniform boxes
    boxes = make_uniform_boxes(nx, ny, nz)
    print(f"{ncpu} uniform subboxes: {nx}x{ny}x{nz}")

    # 12 distinct sequential colormaps (white=low, dark=high)
    cmap_names = [
        'Reds', 'Blues', 'Greens', 'Purples', 'Oranges', 'YlOrBr',
        'RdPu', 'BuGn', 'PuBu', 'OrRd', 'GnBu', 'YlGn',
    ]
    cmaps = [plt.colormaps[name] for name in cmap_names]

    gap = 0.018
    npts = 48

    # View 1: elev=30, azim=45
    print("Rendering view 1...")
    fig1 = plt.figure(figsize=(14, 11))
    ax1 = fig1.add_subplot(111, projection='3d', computed_zorder=False)
    render_scene(ax1, density, boxes, cmaps, elev=30, azim=45, gap=gap, npts=npts, npart=npart)
    plt.tight_layout()
    outfile1 = '/gpfs/kjhan/Hydro/CUBE_HR5/code_cube/misc/ksection_domains_18ranks.png'
    fig1.savefig(outfile1, dpi=200, bbox_inches='tight')
    print(f"Saved: {outfile1}")
    plt.close(fig1)

    # View 2: elev=20, azim=55
    print("Rendering view 2...")
    fig2 = plt.figure(figsize=(14, 11))
    ax2 = fig2.add_subplot(111, projection='3d', computed_zorder=False)
    render_scene(ax2, density, boxes, cmaps, elev=20, azim=55, gap=gap, npts=npts, npart=npart)
    plt.tight_layout()
    outfile2 = '/gpfs/kjhan/Hydro/CUBE_HR5/code_cube/misc/ksection_domains_18ranks_wide.png'
    fig2.savefig(outfile2, dpi=200, bbox_inches='tight')
    print(f"Saved: {outfile2}")
    plt.close(fig2)


if __name__ == '__main__':
    main()
