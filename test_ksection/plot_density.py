#!/usr/bin/env python3
"""
Plot 2D gas density projection from RAMSES output.
Custom binary reader for cuRAMSES (ksection ordering).
Cell-size-aware AMR projection: each cell fills its actual pixel footprint.
Author: Juhan Kim
"""
import sys
import os
import numpy as np
import struct
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def read_record(f):
    """Read one Fortran unformatted record, return raw bytes."""
    raw = f.read(4)
    if len(raw) < 4:
        return None
    n = struct.unpack('i', raw)[0]
    data = f.read(n)
    f.read(4)
    return data


def ri(f):
    """Read a single integer from one record."""
    d = read_record(f)
    if d is None:
        raise EOFError("Unexpected end of file")
    return struct.unpack('i', d[:4])[0]


def rd(f):
    """Read double(s) from one record."""
    d = read_record(f)
    n = len(d) // 8
    vals = struct.unpack(f'{n}d', d)
    return vals[0] if n == 1 else vals


def skip(f, n=1):
    """Skip n records."""
    for _ in range(n):
        read_record(f)


def read_info(output_dir):
    """Parse info_NNNNN.txt."""
    outnum = os.path.basename(output_dir).split('_')[1]
    info = {}
    with open(os.path.join(output_dir, f'info_{outnum}.txt')) as fh:
        for line in fh:
            if '=' in line:
                k, v = line.split('=', 1)
                k, v = k.strip(), v.strip()
                try:
                    info[k] = float(v) if ('.' in v or 'E' in v.upper()) else int(v)
                except ValueError:
                    info[k] = v
    return info


def project_cells_amr(proj_map, res, x_min, y_min, x_ext, y_ext,
                       cell_x, cell_y, cell_rho, dx):
    """Project cells onto map with cell-size-aware area fill.

    Each cell fills all pixels within its spatial footprint [cx-dx/2, cx+dx/2].
    Column density contribution per pixel = rho * dz (where dz = dx).
    """
    dpix_x = x_ext / res
    dpix_y = y_ext / res

    # Filter cells that overlap the view
    in_view = ((cell_x + dx / 2 > x_min) & (cell_x - dx / 2 < x_min + x_ext) &
               (cell_y + dx / 2 > y_min) & (cell_y - dx / 2 < y_min + y_ext))
    if not np.any(in_view):
        return 0

    cx = cell_x[in_view]
    cy = cell_y[in_view]
    rho = cell_rho[in_view]
    ncells = len(cx)
    contribution = rho * dx  # column density dz = dx

    npix_x = max(1, int(round(dx / dpix_x)))
    npix_y = max(1, int(round(dx / dpix_y)))

    if npix_x <= 1 and npix_y <= 1:
        # Sub-pixel cells: point projection (fast vectorized)
        ix = ((cx - x_min) / x_ext * res).astype(int)
        iy = ((cy - y_min) / y_ext * res).astype(int)
        valid = (ix >= 0) & (ix < res) & (iy >= 0) & (iy < res)
        np.add.at(proj_map, (iy[valid], ix[valid]), contribution[valid])
    else:
        # Super-pixel cells: fill the cell's pixel footprint
        ix_lo = np.floor((cx - dx / 2 - x_min) / x_ext * res).astype(int)
        iy_lo = np.floor((cy - dx / 2 - y_min) / y_ext * res).astype(int)
        for di in range(npix_x):
            for dj in range(npix_y):
                ix = ix_lo + di
                iy = iy_lo + dj
                valid = (ix >= 0) & (ix < res) & (iy >= 0) & (iy < res)
                np.add.at(proj_map, (iy[valid], ix[valid]), contribution[valid])

    return ncells


def plot_density(output_dir, full_res=512, zoom_res=1024,
                 zoom_center=None, zoom_width=None):
    outnum = os.path.basename(output_dir).split('_')[1]
    info = read_info(output_dir)
    ncpu = info['ncpu']
    ndim = info['ndim']
    levelmin = info['levelmin']
    nlevelmax = info.get('levelmax', 20)
    aexp = info['aexp']
    z = 1.0 / aexp - 1.0
    boxlen = info['boxlen']
    unit_l = info['unit_l']
    unit_d = info['unit_d']
    boxlen_mpc = boxlen * unit_l / 3.0857e24

    twotondim = 2 ** ndim
    twondim = 2 * ndim

    print(f"Output: {output_dir}")
    print(f"  ncpu={ncpu}, ndim={ndim}, levelmin={levelmin}, levelmax={nlevelmax}")
    print(f"  aexp={aexp:.4f}, z={z:.2f}, box={boxlen_mpc:.2f} Mpc/h")

    # Projection grids
    proj_full = np.zeros((full_res, full_res), dtype=np.float64)
    proj_zoom = np.zeros((zoom_res, zoom_res), dtype=np.float64)

    # Auto-detect zoom region from first pass if not specified
    if zoom_center is None:
        zoom_center = (boxlen / 2, boxlen / 2)
    if zoom_width is None:
        zoom_width = boxlen * 0.3

    zx_min = zoom_center[0] - zoom_width / 2
    zy_min = zoom_center[1] - zoom_width / 2
    zx_ext = zoom_width
    zy_ext = zoom_width

    print(f"  Full box: {full_res}x{full_res}")
    print(f"  Zoom: {zoom_res}x{zoom_res}, center=({zoom_center[0]:.3f},{zoom_center[1]:.3f}), "
          f"width={zoom_width:.3f} [{zoom_width*boxlen_mpc:.3f} Mpc/h]")

    total_leaf = 0

    for icpu in range(1, ncpu + 1):
        amr_f = os.path.join(output_dir, f'amr_{outnum}.out{icpu:05d}')
        hyd_f = os.path.join(output_dir, f'hydro_{outnum}.out{icpu:05d}')

        # ===================== READ HYDRO =====================
        hydro_data = {}
        with open(hyd_f, 'rb') as fh:
            ncpu_h = ri(fh)
            nvar_h = ri(fh)
            ndim_h = ri(fh)
            nlevelmax_h = ri(fh)
            nboundary_h = ri(fh)
            gamma_h = rd(fh)

            for ilevel in range(nlevelmax_h):
                for ibound in range(ncpu_h + nboundary_h):
                    _ilev = ri(fh)
                    ncache = ri(fh)

                    if ncache == 0:
                        continue

                    is_active = (ibound + 1 == icpu)

                    for ind in range(twotondim):
                        for ivar in range(nvar_h):
                            if is_active and ivar == 0:
                                rec = read_record(fh)
                                arr = np.frombuffer(rec, dtype=np.float64).copy()
                                hydro_data[(ilevel, ind)] = arr
                            else:
                                skip(fh, 1)

        # ===================== READ AMR =====================
        with open(amr_f, 'rb') as fa:
            ncpu_a = ri(fa)
            ndim_a = ri(fa)
            rec = read_record(fa)
            nx, ny, nz = struct.unpack('iii', rec)
            nlevelmax_a = ri(fa)
            ngridmax_a = ri(fa)
            nboundary_a = ri(fa)
            ngrid_current = ri(fa)
            boxlen_a = rd(fa)

            skip(fa, 3)  # noutput/iout/ifout, tout, aout
            skip(fa, 3)  # t, dtold, dtnew
            skip(fa, 2)  # nstep/nstep_coarse, const/mass_tot_0/rho_tot
            skip(fa, 2)  # omega_m..boxlen_ini, aexp..epot_tot_old
            skip(fa, 1)  # mass_sph

            skip(fa, 1)  # headl
            skip(fa, 1)  # taill
            numbl_rec = read_record(fa)
            numbl_all = np.frombuffer(numbl_rec, dtype=np.int32).reshape(
                ncpu_a, nlevelmax_a, order='F')
            skip(fa, 1)  # numbtot

            if nboundary_a > 0:
                skip(fa, 3)

            skip(fa, 1)  # headf,tailf,numbf,used_mem,used_mem_tot

            ordering_rec = read_record(fa)
            ordering = ordering_rec.decode('ascii', errors='ignore').strip()

            if 'ksection' in ordering.lower():
                skip(fa, 10)
            elif 'bisection' in ordering.lower():
                skip(fa, 5)
            else:
                skip(fa, 1)

            skip(fa, 3)  # coarse: son, flag1, cpu_map

            nrec_per_ibound = 3 + ndim_a + 1 + 2 * ndim_a + 3 * twotondim

            for ilevel in range(nlevelmax_a):
                for ibound in range(ncpu_a + nboundary_a):
                    if ibound < ncpu_a:
                        ncache = int(numbl_all[ibound, ilevel])
                    else:
                        ncache = 0

                    if ncache == 0:
                        continue

                    is_active = (ibound + 1 == icpu)

                    if not is_active:
                        skip(fa, nrec_per_ibound)
                        continue

                    skip(fa, 3)  # ind_grid, next, prev

                    xg = []
                    for idim in range(ndim_a):
                        rec = read_record(fa)
                        arr = np.frombuffer(rec, dtype=np.float64).copy()
                        xg.append(arr)

                    skip(fa, 1)            # father
                    skip(fa, 2 * ndim_a)   # nbor

                    son = []
                    for ind in range(twotondim):
                        rec = read_record(fa)
                        arr = np.frombuffer(rec, dtype=np.int32).copy()
                        son.append(arr)

                    skip(fa, twotondim)  # cpu_map
                    skip(fa, twotondim)  # ref_map

                    dx = boxlen * 0.5 ** (ilevel + 1)

                    for ind in range(twotondim):
                        iz = ind // 4
                        iy = (ind - 4 * iz) // 2
                        ix = ind - 2 * iy - 4 * iz

                        ox = (float(ix) - 0.5) * dx
                        oy = (float(iy) - 0.5) * dx

                        is_leaf = (son[ind] == 0)
                        if not np.any(is_leaf):
                            continue

                        key = (ilevel, ind)
                        if key not in hydro_data:
                            continue

                        rho = hydro_data[key][is_leaf]
                        cell_x = xg[0][is_leaf] + ox
                        cell_y = xg[1][is_leaf] + oy

                        # Project to full-box map
                        project_cells_amr(proj_full, full_res,
                                          0, 0, boxlen, boxlen,
                                          cell_x, cell_y, rho, dx)

                        # Project to zoom map
                        project_cells_amr(proj_zoom, zoom_res,
                                          zx_min, zy_min, zx_ext, zy_ext,
                                          cell_x, cell_y, rho, dx)

                        total_leaf += int(np.sum(is_leaf))

        if icpu % 3 == 0 or icpu == ncpu:
            print(f"  CPU {icpu}/{ncpu}, leaf cells: {total_leaf:,}")

    print(f"  Total leaf cells: {total_leaf:,}")

    # ===================== PLOTTING =====================
    fig, axes = plt.subplots(1, 2, figsize=(16, 7))

    # --- Full box ---
    ax = axes[0]
    nonzero = proj_full[proj_full > 0]
    if len(nonzero) > 0:
        vmin_f = np.percentile(nonzero, 1)
        vmax_f = np.percentile(nonzero, 99.9)
        proj_full[proj_full <= 0] = vmin_f
        im = ax.imshow(np.log10(proj_full), origin='lower',
                       extent=[0, boxlen_mpc, 0, boxlen_mpc],
                       cmap='inferno',
                       vmin=np.log10(vmin_f), vmax=np.log10(vmax_f))
        plt.colorbar(im, ax=ax, label=r'log$_{10}$ $\Sigma$ [code units]')
    ax.set_xlabel('x [Mpc/h comoving]')
    ax.set_ylabel('y [Mpc/h comoving]')
    ax.set_title(f'Gas Column Density (full box)\nz = {z:.1f}, a = {aexp:.4f}')

    # Mark zoom region on full box
    zx0 = zx_min * boxlen_mpc
    zy0 = zy_min * boxlen_mpc
    zw = zoom_width * boxlen_mpc
    rect = plt.Rectangle((zx0, zy0), zw, zw,
                          linewidth=1.5, edgecolor='cyan', facecolor='none',
                          linestyle='--')
    ax.add_patch(rect)

    # --- Zoom region (cell-size-aware) ---
    ax = axes[1]
    nonzero_z = proj_zoom[proj_zoom > 0]
    if len(nonzero_z) > 0:
        vmin_z = np.percentile(nonzero_z, 1)
        vmax_z = np.percentile(nonzero_z, 99.9)
        proj_zoom[proj_zoom <= 0] = vmin_z
        zx0_mpc = zx_min * boxlen_mpc
        zy0_mpc = zy_min * boxlen_mpc
        zx1_mpc = (zx_min + zx_ext) * boxlen_mpc
        zy1_mpc = (zy_min + zy_ext) * boxlen_mpc
        im2 = ax.imshow(np.log10(proj_zoom), origin='lower',
                        extent=[zx0_mpc, zx1_mpc, zy0_mpc, zy1_mpc],
                        cmap='inferno',
                        vmin=np.log10(vmin_z), vmax=np.log10(vmax_z))
        plt.colorbar(im2, ax=ax, label=r'log$_{10}$ $\Sigma$ [code units]')

        # Annotate cell sizes
        for lev in range(levelmin, nlevelmax + 1):
            dx_lev = boxlen * 0.5 ** lev
            dx_pix = dx_lev / zoom_width * zoom_res
            if dx_pix >= 1.5:
                dx_kpc = dx_lev * boxlen_mpc * 1000  # kpc/h
                ax.text(0.02, 0.97 - 0.04 * (lev - levelmin), f'L{lev}: {dx_kpc:.1f} kpc/h ({dx_pix:.0f} pix)',
                        transform=ax.transAxes, fontsize=7,
                        color='white', va='top',
                        bbox=dict(boxstyle='round,pad=0.15', facecolor='black', alpha=0.5))

    ax.set_xlabel('x [Mpc/h comoving]')
    ax.set_ylabel('y [Mpc/h comoving]')
    ax.set_title(f'Gas Column Density (zoom, AMR cell-size)\nz = {z:.1f}, a = {aexp:.4f}')

    plt.tight_layout()
    outfile = os.path.join(os.path.dirname(output_dir),
                           f'density_{os.path.basename(output_dir)}.png')
    plt.savefig(outfile, dpi=150, bbox_inches='tight')
    print(f"\nSaved: {outfile}")
    plt.close()


if __name__ == '__main__':
    if len(sys.argv) < 2:
        d = '/gpfs/kjhan/Hydro/CUBE_HR5/code_cube/test_ksection/run_zoomin_lv11/output_00002'
    else:
        d = sys.argv[1]
    plot_density(d)
