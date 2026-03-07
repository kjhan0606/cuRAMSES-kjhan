#!/usr/bin/env python3
"""
Generate a multi-redshift Grackle cooling table.

Calls grackle_cooling_grid for each redshift snapshot and packs
all results into a single binary file for RAMSES runtime interpolation.

Multi-z binary format:
  int32:           n_z
  int32:           L, M, N
  float64[n_z]:    z_values (sorted ascending)
  float64[7]:      log_dmin, log_dmax, log_zmin, log_zmax, log_tmin, log_tmax, (unused=0)
  float64[N]:      log10(temperature) grid
  For iz = 0..n_z-1:
    float64[L*M*N]:  net_rate  [erg/cm^3/s]
    float64[L*M*N]:  mmw       [dimensionless]
"""
import argparse
import struct
import subprocess
import tempfile
import os
import sys
import numpy as np


def read_single_z_binary(path):
    """Read a single-z Grackle binary table."""
    with open(path, 'rb') as f:
        L, M, N = struct.unpack('iii', f.read(12))
        hdr = struct.unpack('7d', f.read(56))  # z, log_d{min,max}, log_z{min,max}, log_t{min,max}
        log_T = np.frombuffer(f.read(N * 8), dtype=np.float64).copy()
        net_rate = np.frombuffer(f.read(L * M * N * 8), dtype=np.float64).copy()
        mmw = np.frombuffer(f.read(L * M * N * 8), dtype=np.float64).copy()
    return L, M, N, hdr, log_T, net_rate, mmw


def main():
    parser = argparse.ArgumentParser(description='Generate multi-z Grackle cooling table')
    parser.add_argument('--output', default='grackle_multi_z.bin',
                        help='Output binary file (default: grackle_multi_z.bin)')
    parser.add_argument('--generator',
                        default='/home/kjhan/BACKUP/Eunha.A1/Prerequisites/Cooling_Grackle/grackle_cooling_grid',
                        help='Path to grackle_cooling_grid executable')
    parser.add_argument('--cloudy',
                        default='/home/kjhan/grackle/grackle_data_files/input/CloudyData_UVB=HM2012.h5',
                        help='Path to Cloudy HM2012 data file')
    parser.add_argument('--zmin', type=float, default=0.0, help='Minimum redshift (for linspace mode)')
    parser.add_argument('--zmax', type=float, default=10.0, help='Maximum redshift (for linspace mode)')
    parser.add_argument('--nz', type=int, default=21, help='Number of redshift snapshots (for linspace mode)')
    parser.add_argument('--z_list', type=str, default=None,
                        help='Comma-separated z values (overrides zmin/zmax/nz). '
                             'Default: 0.0,0.5,1.0,1.5,2.0,3.0,4.0,5.0,7.0,10.0')
    args = parser.parse_args()

    # Default: low-z dense, high-z sparse
    if args.z_list is not None:
        z_values = np.array(sorted(float(z) for z in args.z_list.split(',')))
    else:
        z_values = np.array([0.0, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0, 7.0, 10.0])
    nz_total = len(z_values)
    print(f'Generating {nz_total} redshift snapshots: z = {z_values}')

    all_rates = []
    all_mmw = []
    L_ref = M_ref = N_ref = None
    hdr_ref = log_T_ref = None

    for iz, z in enumerate(z_values):
        print(f'  [{iz+1}/{nz_total}] z = {z:.2f} ...', end=' ', flush=True)
        with tempfile.NamedTemporaryFile(suffix='.bin', delete=False) as tmp:
            tmppath = tmp.name

        try:
            result = subprocess.run(
                [args.generator, tmppath, f'{z:.6f}', args.cloudy],
                capture_output=True, text=True, timeout=300
            )
            if result.returncode != 0:
                print(f'FAILED (rc={result.returncode})')
                print(result.stderr)
                sys.exit(1)

            L, M, N, hdr, log_T, net_rate, mmw = read_single_z_binary(tmppath)

            if L_ref is None:
                L_ref, M_ref, N_ref = L, M, N
                hdr_ref = hdr
                log_T_ref = log_T
            else:
                assert (L, M, N) == (L_ref, M_ref, N_ref), \
                    f'Grid mismatch at z={z}: ({L},{M},{N}) vs ({L_ref},{M_ref},{N_ref})'

            all_rates.append(net_rate)
            all_mmw.append(mmw)
            print(f'OK  (L={L}, M={M}, N={N})')
        finally:
            os.unlink(tmppath)

    # Write multi-z binary
    nz = len(z_values)
    total_bytes = 4 + 12 + nz*8 + 7*8 + N_ref*8 + nz*2*L_ref*M_ref*N_ref*8
    print(f'\nWriting {args.output} ({total_bytes/1e6:.1f} MB) ...')

    with open(args.output, 'wb') as f:
        f.write(struct.pack('i', nz))
        f.write(struct.pack('iii', L_ref, M_ref, N_ref))
        f.write(z_values.astype(np.float64).tobytes())
        f.write(struct.pack('7d', *hdr_ref))
        f.write(log_T_ref.astype(np.float64).tobytes())
        for iz in range(nz):
            f.write(all_rates[iz].astype(np.float64).tobytes())
            f.write(all_mmw[iz].astype(np.float64).tobytes())

    print(f'Done. {nz} z-snapshots, grid {L_ref}x{M_ref}x{N_ref}')


if __name__ == '__main__':
    main()
