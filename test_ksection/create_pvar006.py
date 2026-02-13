#!/usr/bin/env python3
"""
Create ic_pvar_00006 (zoom geometry scalar) for RAMSES zoom-in simulation.
- Level 7 (base): read ic_refmap to get zoom mask (1.0 for zoom, 0.0 for background)
- Levels 8-11 (zoom sub-levels): all cells = 1.0 (entire grid is zoom region)

GRAFIC2 binary format:
  Record 1 (header): [rec_marker] n1,n2,n3(int), dx(float), x1o,x2o,x3o(float),
                      astart(float), omegam,omegav,h0(float) [rec_marker]
  Records 2..n3+1: [rec_marker] n1*n2 float32 values [rec_marker]
"""
import struct
import numpy as np
import os
import sys

IC_DIR = "IC_zoomin_lv11"
LEVELS = [7, 8, 9, 10, 11]
LEVEL_DIRS = {
    7:  "level_007",
    8:  "level_008",
    9:  "level_009",
    10: "level_010",
    11: "level_011",
}

def read_grafic2_header(fname):
    """Read GRAFIC2 header, return (n1,n2,n3, header_bytes)"""
    with open(fname, 'rb') as f:
        rec_len = struct.unpack('i', f.read(4))[0]
        assert rec_len == 44, f"Expected 44-byte header, got {rec_len}"
        n1, n2, n3 = struct.unpack('3i', f.read(12))
        rest = f.read(32)  # dx, x1o,x2o,x3o, astart, omegam,omegav,h0
        rec_end = struct.unpack('i', f.read(4))[0]
        assert rec_end == 44
        # Return full header bytes for rewriting
        f.seek(0)
        header_bytes = f.read(4 + 44 + 4)  # rec_start + data + rec_end
    return n1, n2, n3, header_bytes

def read_grafic2_data(fname, n1, n2, n3):
    """Read full 3D data from GRAFIC2 file"""
    data = np.zeros((n3, n2, n1), dtype=np.float32)
    with open(fname, 'rb') as f:
        f.read(4 + 44 + 4)  # skip header
        for i3 in range(n3):
            rec_len = struct.unpack('i', f.read(4))[0]
            assert rec_len == n1 * n2 * 4
            plane = np.frombuffer(f.read(n1 * n2 * 4), dtype=np.float32).reshape(n2, n1)
            data[i3] = plane
            rec_end = struct.unpack('i', f.read(4))[0]
    return data

def write_grafic2(fname, header_bytes, data):
    """Write GRAFIC2 file with given header and 3D float32 data"""
    n3, n2, n1 = data.shape
    with open(fname, 'wb') as f:
        f.write(header_bytes)
        for i3 in range(n3):
            plane = data[i3].astype(np.float32).tobytes()
            rec_len = struct.pack('i', n1 * n2 * 4)
            f.write(rec_len)
            f.write(plane)
            f.write(rec_len)

for level in LEVELS:
    level_dir = os.path.join(IC_DIR, LEVEL_DIRS[level])
    deltab_file = os.path.join(level_dir, "ic_deltab")
    refmap_file = os.path.join(level_dir, "ic_refmap")
    pvar_file = os.path.join(level_dir, "ic_pvar_00006")

    # Read header from ic_deltab
    n1, n2, n3, header_bytes = read_grafic2_header(deltab_file)
    print(f"Level {level}: {n1}x{n2}x{n3}")

    if level == 7:
        # Base level: read ic_refmap to get zoom mask
        refmap_data = read_grafic2_data(refmap_file, n1, n2, n3)
        # Convert to zoom scalar: nonzero -> 1.0, zero -> 0.0
        pvar_data = np.where(refmap_data != 0, 1.0, 0.0).astype(np.float32)
        nonzero = np.count_nonzero(pvar_data)
        total = n1 * n2 * n3
        print(f"  Zoom cells: {nonzero}/{total} ({100*nonzero/total:.1f}%)")
    else:
        # Zoom sub-levels: all cells are in zoom region
        pvar_data = np.ones((n3, n2, n1), dtype=np.float32)
        print(f"  All {n1*n2*n3} cells set to 1.0 (zoom region)")

    write_grafic2(pvar_file, header_bytes, pvar_data)
    print(f"  Written: {pvar_file}")

print("\nDone! ic_pvar_00006 created for all levels.")
