#!/usr/bin/env python3
"""
Townsend (2009) exact cooling validation test.

Two tests:
  Test 1: Piecewise power-law cooling function (Townsend 2009, Table 1)
          — reproduces the isochoric cooling test from the paper
  Test 2: Grackle table-based cooling
          — validates our Y-function implementation against tiny-dt reference

Method:
  - Y-function exact integration (Townsend 2009, Eq 6-10)
  - Reference: RK4 numerical integration with very small dt
  - Compares T(t) for various dt sizes to verify dt-independence
"""

import struct
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


# ── Physical constants (CGS) ──
k_B = 1.3806504e-16      # Boltzmann [erg/K]
m_H = 1.6726219e-24      # proton mass [g]
GAMMA_AD = 5.0/3.0
X_H = 0.76
Y_He = 0.24
Z_SOLAR = 0.02
Myr = 3.15576e13         # 1 Myr in seconds


# ═══════════════════════════════════════════════════════════════
# TEST 1: Townsend (2009) piecewise power-law
# ═══════════════════════════════════════════════════════════════

# Table 1 from Townsend (2009, ApJS 181, 391)
# Piecewise power-law: Lambda_N = C_j * (T / T_j)^alpha_j  [erg cm^3 s^-1]
PL_SEGMENTS = [
    # T_j [K],     C [erg cm^3/s],  alpha
    (1.0e2,   3.720e-26,  7.000),   # segment 0
    (2.0e3,   3.720e-26,  7.000),   # segment 1
    (8.0e3,   1.000e-24,  0.000),   # segment 2 (flat)
    (1.0e5,   4.684e-23, -1.000),   # segment 3
    (4.0e7,   1.170e-27,  0.500),   # segment 4
]
PL_T_BOUNDS = [1e2, 2e3, 8e3, 1e5, 4e7, 1e9]


def pl_cooling_rate(T):
    """Piecewise power-law Lambda_N(T) [erg cm^3 s^-1]."""
    T = np.atleast_1d(np.asarray(T, dtype=np.float64))
    result = np.zeros_like(T)
    for iseg, (T_j, C_j, alpha_j) in enumerate(PL_SEGMENTS):
        T_lo = PL_T_BOUNDS[iseg]
        T_hi = PL_T_BOUNDS[iseg + 1]
        if iseg == len(PL_SEGMENTS) - 1:
            mask = (T >= T_lo)
        else:
            mask = (T >= T_lo) & (T < T_hi)
        result[mask] = C_j * (T[mask] / T_j) ** alpha_j
    return result


def numerical_integrate_pl(T_init, nH, dt_total, nsteps):
    """
    RK4 numerical integration of isochoric cooling.
    Townsend convention: e = nH * kT / (gamma-1)
    de/dt = -nH^2 * Lambda_N(T)
    => dT/dt = -(gamma-1) * nH * Lambda_N(T) / kB
    """
    dt = dt_total / nsteps
    T = float(T_init)

    t_arr = [0.0]
    T_arr = [T]
    T_floor = PL_T_BOUNDS[0]

    for step in range(nsteps):
        def dTdt(Tv):
            Tv = max(Tv, T_floor)
            LN = float(pl_cooling_rate(np.array([Tv]))[0])
            return -(GAMMA_AD - 1.0) * nH * LN / k_B

        k1 = dTdt(T) * dt
        k2 = dTdt(T + 0.5 * k1) * dt
        k3 = dTdt(T + 0.5 * k2) * dt
        k4 = dTdt(T + k3) * dt
        T = T + (k1 + 2*k2 + 2*k3 + k4) / 6.0
        T = max(T, T_floor)

        t_arr.append((step + 1) * dt)
        T_arr.append(T)

    return np.array(t_arr), np.array(T_arr)


def build_Y_from_rates(T_grid, L_spec):
    """Build Y-function from T grid and specific cooling rates (Townsend Eq 6-9)."""
    N = len(T_grid)
    u_grid = k_B * T_grid / (GAMMA_AD - 1.0)  # e/nH [erg]

    k_ref = N // 2
    u_ref = u_grid[k_ref]
    L_ref = L_spec[k_ref]

    # Find non-zero reference
    if abs(L_ref) < 1e-30:
        for k in range(N):
            if abs(L_spec[k]) > 1e-30:
                k_ref = k
                u_ref = u_grid[k_ref]
                L_ref = L_spec[k_ref]
                break
    if abs(L_ref) < 1e-30:
        return np.zeros(N), k_ref, u_ref, 1e-30

    norm = L_ref / u_ref
    Y = np.zeros(N)

    def calc_dY(u1, u2, L1, L2):
        du = u2 - u1
        dL = L2 - L1
        if abs(L1) < 1e-30 or abs(L2) < 1e-30 or L1 * L2 < 0:
            return 0.0
        if abs(dL) < 1e-8 * abs(L1 + L2) * 0.5 + 1e-30:
            L_avg = 0.5 * (L1 + L2)
            return norm * du / L_avg if abs(L_avg) > 1e-30 else 0.0
        A = du / dL
        B = (u2 * L1 - u1 * L2) / dL
        arg1 = B + u1
        arg2 = B + u2
        if arg1 * arg2 <= 0 or abs(arg1) < 1e-30 or abs(arg2) < 1e-30:
            L_avg = 0.5 * (L1 + L2)
            return norm * du / L_avg if abs(L_avg) > 1e-30 else 0.0
        return norm * A * np.log(arg2 / arg1)

    # Upward
    for k in range(k_ref + 1, N):
        dY = calc_dY(u_grid[k-1], u_grid[k], L_spec[k-1], L_spec[k])
        Y[k] = Y[k-1] + dY

    # Downward
    for k in range(k_ref - 1, -1, -1):
        dY = calc_dY(u_grid[k+1], u_grid[k], L_spec[k+1], L_spec[k])
        Y[k] = Y[k+1] + dY

    return Y, k_ref, u_ref, L_ref


def exact_integrate_generic(T_grid, Y_grid, L_spec, u_ref, L_ref, T_init, dt):
    """
    Townsend exact integration using precomputed Y-function.

    Returns T_final after isochoric cooling for time dt.
    """
    log_T = np.log10(T_grid)
    N = len(T_grid)
    u_grid = k_B * T_grid / (GAMMA_AD - 1.0)

    # Get initial specific cooling rate
    log_T_init = np.log10(T_init)
    Lambda_init = np.interp(log_T_init, log_T, L_spec)

    if abs(Lambda_init) < 1e-30:
        return T_init

    is_cooling = (Lambda_init > 0)

    # Y at T_init
    Y_init = np.interp(log_T_init, log_T, Y_grid)

    # Delta_Y = (|L_ref| / u_ref) * dt
    Delta_Y = (abs(L_ref) / u_ref) * dt

    if is_cooling:
        Y_target = Y_init - Delta_Y
    else:
        Y_target = Y_init + Delta_Y

    # Find T where Y = Y_target
    if is_cooling:
        # Search downward: Y generally decreases with decreasing T
        idx_init = np.searchsorted(log_T, log_T_init)
        idx_init = max(1, min(idx_init, N - 1))

        for k in range(idx_init, 0, -1):
            Y_lo = Y_grid[k - 1]
            Y_hi = Y_grid[k]
            if (Y_target - Y_lo) * (Y_target - Y_hi) <= 0:
                if abs(Y_hi - Y_lo) > 1e-30:
                    frac = (Y_target - Y_lo) / (Y_hi - Y_lo)
                    frac = max(0.0, min(1.0, frac))
                else:
                    frac = 0.5
                T_new = 10.0 ** (log_T[k-1] + frac * (log_T[k] - log_T[k-1]))
                return max(T_grid[0], min(T_new, T_grid[-1]))

        # Check equilibrium crossings
        for k in range(idx_init, 0, -1):
            if L_spec[k] > 0 and L_spec[k-1] <= 0:
                if abs(L_spec[k] - L_spec[k-1]) > 1e-30:
                    frac = L_spec[k] / (L_spec[k] - L_spec[k-1])
                    frac = max(0.0, min(1.0, frac))
                    return 10.0 ** (log_T[k] + frac * (log_T[k-1] - log_T[k]))
                return T_grid[k]
        return T_grid[0]  # floor
    else:
        # Heating: search upward
        idx_init = np.searchsorted(log_T, log_T_init)
        idx_init = max(0, min(idx_init, N - 2))

        for k in range(idx_init, N - 1):
            Y_lo = Y_grid[k]
            Y_hi = Y_grid[k + 1]
            if (Y_target - Y_lo) * (Y_target - Y_hi) <= 0:
                if abs(Y_hi - Y_lo) > 1e-30:
                    frac = (Y_target - Y_lo) / (Y_hi - Y_lo)
                    frac = max(0.0, min(1.0, frac))
                else:
                    frac = 0.5
                T_new = 10.0 ** (log_T[k] + frac * (log_T[k+1] - log_T[k]))
                return max(T_grid[0], min(T_new, T_grid[-1]))
        return T_grid[-1]


def test1_piecewise_powerlaw():
    """Test 1: Townsend piecewise power-law isochoric cooling."""
    print("=" * 70)
    print("TEST 1: Piecewise power-law (Townsend 2009)")
    print("=" * 70)

    nH = 1.0  # cm^-3
    T_init = 1e8  # K

    # Build fine T grid and cooling rates
    N_T = 4000
    T_grid = np.logspace(np.log10(PL_T_BOUNDS[0]), np.log10(PL_T_BOUNDS[-1]), N_T)
    LN_grid = pl_cooling_rate(T_grid)

    # Specific rate in Townsend convention: L_spec = nH * Lambda_N
    # Because dT/dt = -(gamma-1) * nH * Lambda_N / kB
    # and u = kT / (gamma-1), so du/dt = -nH * Lambda_N ≡ -L_spec
    L_spec = nH * LN_grid

    # t_cool = kT / ((gamma-1) * nH * Lambda_N) = u / L_spec
    u_init = k_B * T_init / (GAMMA_AD - 1.0)
    L_init = np.interp(np.log10(T_init), np.log10(T_grid), L_spec)
    t_cool = u_init / L_init

    t_total = 10.0 * t_cool

    print(f"  nH = {nH:.1f} cm^-3,  T_init = {T_init:.1e} K")
    print(f"  t_cool = {t_cool:.3e} s = {t_cool/Myr:.3f} Myr")
    print(f"  t_total = 10 * t_cool = {t_total:.3e} s = {t_total/Myr:.3f} Myr")

    # Build Y-function
    Y_grid, k_ref, u_ref, L_ref = build_Y_from_rates(T_grid, L_spec)

    # Reference: RK4 with 500,000 steps
    n_ref = 500000
    print(f"\n  Computing RK4 reference ({n_ref} steps)...", end=" ", flush=True)
    t_ref, T_ref = numerical_integrate_pl(T_init, nH, t_total, n_ref)
    print(f"T_final = {T_ref[-1]:.6e} K")

    # Exact integration with various dt
    dt_factors = [0.01, 0.1, 0.5, 1.0, 2.0, 5.0, 10.0]
    print(f"\n  {'dt/t_cool':>12s}  {'N_steps':>8s}  {'T_exact':>14s}  {'T_RK4_ref':>14s}  {'rel_error':>12s}")
    print("  " + "-" * 70)

    exact_trajectories = {}
    for dtf in dt_factors:
        dt_step = dtf * t_cool
        n_steps = max(1, int(np.ceil(t_total / dt_step)))
        T = T_init
        t_traj = [0.0]
        T_traj = [T_init]
        for step in range(n_steps):
            dt_actual = min(dt_step, t_total - step * dt_step)
            if dt_actual <= 0:
                break
            T = exact_integrate_generic(T_grid, Y_grid, L_spec, u_ref, L_ref, T, dt_actual)
            t_traj.append((step + 1) * dt_step)
            T_traj.append(T)
        exact_trajectories[dtf] = (np.array(t_traj), np.array(T_traj))

        T_ref_final = T_ref[-1]
        if T_ref_final > PL_T_BOUNDS[0] * 1.01:
            rel_err = abs(T - T_ref_final) / abs(T_ref_final)
        else:
            rel_err = abs(T - T_ref_final) / max(abs(T_ref_final), 1.0)
        print(f"  {dtf:12.2f}  {n_steps:8d}  {T:14.6e}  {T_ref_final:14.6e}  {rel_err:12.4e}")

    return {
        't_ref': t_ref, 'T_ref': T_ref, 't_cool': t_cool,
        'exact': exact_trajectories, 'T_grid': T_grid, 'LN_grid': LN_grid
    }


# ═══════════════════════════════════════════════════════════════
# TEST 2: Grackle table-based cooling
# ═══════════════════════════════════════════════════════════════

def read_grackle_table(path):
    """Read single-z or multi-z Grackle binary table (first z-snapshot)."""
    with open(path, 'rb') as f:
        first_int = struct.unpack('i', f.read(4))[0]

    is_multi_z = False
    if first_int <= 200 and first_int != 60 and first_int != 40:
        # Could be multi-z: check next ints
        with open(path, 'rb') as f:
            nz_cand = struct.unpack('i', f.read(4))[0]
            L_cand, M_cand, N_cand = struct.unpack('iii', f.read(12))
            if L_cand >= 10 and M_cand >= 10 and N_cand >= 10 and nz_cand != L_cand:
                is_multi_z = True

    if is_multi_z:
        with open(path, 'rb') as f:
            nz = struct.unpack('i', f.read(4))[0]
            L, M, N = struct.unpack('iii', f.read(12))
            z_values = np.frombuffer(f.read(nz * 8), dtype=np.float64).copy()
            # Fixed header: log_dmin, log_dmax, log_zmin, log_zmax, log_tmin, log_tmax, unused
            log_dmin, log_dmax, log_zmin, log_zmax, log_tmin, log_tmax, _ = \
                struct.unpack('7d', f.read(56))
            log_T = np.frombuffer(f.read(N * 8), dtype=np.float64).copy()
            # Read first z snapshot
            raw_rate = np.frombuffer(f.read(L * M * N * 8), dtype=np.float64).copy().reshape(L, M, N)
            raw_mmw = np.frombuffer(f.read(L * M * N * 8), dtype=np.float64).copy().reshape(L, M, N)
            z_val = z_values[0]
    else:
        with open(path, 'rb') as f:
            L, M, N = struct.unpack('iii', f.read(12))
            z_val, log_dmin, log_dmax, log_zmin, log_zmax, log_tmin, log_tmax = \
                struct.unpack('7d', f.read(56))
            log_T = np.frombuffer(f.read(N * 8), dtype=np.float64).copy()
            raw_rate = np.frombuffer(f.read(L * M * N * 8), dtype=np.float64).copy().reshape(L, M, N)
            raw_mmw = np.frombuffer(f.read(L * M * N * 8), dtype=np.float64).copy().reshape(L, M, N)

    # Build grids
    density = np.logspace(log_dmin, log_dmax, L)
    metallicity = np.logspace(log_zmin, log_zmax, M)
    temperature = 10.0 ** log_T

    # Convert net_rate: erg/cm^3/s -> erg/g/s
    # Negate: Grackle negative=cooling -> positive=cooling
    net_rate = np.zeros_like(raw_rate)
    for i in range(L):
        nH_val = density[i]
        for j in range(M):
            Z_metal = metallicity[j] * Z_SOLAR
            Z_metal = max(0.0, min(1.0, Z_metal))
            f_scale = (1.0 - Z_metal) / (X_H + Y_He)
            mf_X = X_H * f_scale
            rho = nH_val * m_H / mf_X
            net_rate[i, j, :] = -raw_rate[i, j, :] / rho

    return {
        'L': L, 'M': M, 'N': N, 'z': z_val,
        'density': density, 'metallicity': metallicity, 'temperature': temperature,
        'log_T': log_T, 'net_rate': net_rate, 'mmw': raw_mmw,
        'is_multi_z': is_multi_z
    }


def interp_table_3d(table, nH, Z_met, T, field='net_rate'):
    """Trilinear interpolation in log-space."""
    data = table[field]  # shape (L, M, N)
    den, met, temp = table['density'], table['metallicity'], table['temperature']

    def find_idx(arr, val):
        n = len(arr)
        if val <= arr[0]: return 0, 0.0
        if val >= arr[-1]: return n - 2, 1.0
        idx = np.searchsorted(arr, val) - 1
        idx = max(0, min(idx, n - 2))
        log_val = np.log10(max(val, 1e-30))
        log_lo = np.log10(max(arr[idx], 1e-30))
        log_hi = np.log10(max(arr[idx+1], 1e-30))
        if abs(log_hi - log_lo) > 1e-30:
            frac = (log_val - log_lo) / (log_hi - log_lo)
        else:
            frac = 0.0
        return idx, max(0.0, min(1.0, frac))

    i, fn = find_idx(den, nH)
    j, fz = find_idx(met, Z_met)
    k, ft = find_idx(temp, T)

    c000 = data[i,   j,   k  ];  c001 = data[i,   j,   k+1]
    c010 = data[i,   j+1, k  ];  c011 = data[i,   j+1, k+1]
    c100 = data[i+1, j,   k  ];  c101 = data[i+1, j,   k+1]
    c110 = data[i+1, j+1, k  ];  c111 = data[i+1, j+1, k+1]

    c00 = c000 + ft * (c001 - c000);  c01 = c010 + ft * (c011 - c010)
    c10 = c100 + ft * (c101 - c100);  c11 = c110 + ft * (c111 - c110)
    c0 = c00 + fz * (c01 - c00);  c1 = c10 + fz * (c11 - c10)

    return c0 + fn * (c1 - c0)


def build_Y_grackle(table, i_d, j_z):
    """Build Y-function for one (nH, Z) grid cell using Grackle table."""
    N = table['N']
    T_arr = table['temperature']
    L_arr = table['net_rate'][i_d, j_z, :]  # erg/g/s
    mu_arr = table['mmw'][i_d, j_z, :]

    u_arr = k_B * T_arr / ((GAMMA_AD - 1.0) * mu_arr * m_H)

    k_ref = N // 2
    u_ref = u_arr[k_ref]
    L_ref = L_arr[k_ref]
    if abs(L_ref) < 1e-30:
        for k in range(N):
            if abs(L_arr[k]) > 1e-30:
                k_ref = k; u_ref = u_arr[k]; L_ref = L_arr[k]; break
    if abs(L_ref) < 1e-30:
        return np.zeros(N), k_ref, u_ref, 1e-30

    norm = L_ref / u_ref
    Y = np.zeros(N)

    def calc_dY(u1, u2, L1, L2):
        du = u2 - u1; dL = L2 - L1
        if abs(L1) < 1e-30 or abs(L2) < 1e-30 or L1 * L2 < 0: return 0.0
        if abs(dL) < 1e-8 * abs(L1 + L2) * 0.5 + 1e-30:
            L_avg = 0.5 * (L1 + L2)
            return norm * du / L_avg if abs(L_avg) > 1e-30 else 0.0
        A = du / dL; B = (u2 * L1 - u1 * L2) / dL
        arg1 = B + u1; arg2 = B + u2
        if arg1 * arg2 <= 0 or abs(arg1) < 1e-30 or abs(arg2) < 1e-30:
            L_avg = 0.5 * (L1 + L2)
            return norm * du / L_avg if abs(L_avg) > 1e-30 else 0.0
        return norm * A * np.log(arg2 / arg1)

    for k in range(k_ref + 1, N):
        Y[k] = Y[k-1] + calc_dY(u_arr[k-1], u_arr[k], L_arr[k-1], L_arr[k])
    for k in range(k_ref - 1, -1, -1):
        Y[k] = Y[k+1] + calc_dY(u_arr[k+1], u_arr[k], L_arr[k+1], L_arr[k])

    return Y, k_ref, u_ref, L_ref


def exact_integrate_grackle(table, nH, Z_met, T_init, dt):
    """Townsend exact integration using Grackle table with bilinear Y interpolation."""
    T_arr = table['temperature']
    log_T = np.log10(T_arr)
    N = table['N']
    den, met = table['density'], table['metallicity']

    def find_idx(arr, val):
        n = len(arr)
        if val <= arr[0]: return 0, 0.0
        if val >= arr[-1]: return n - 2, 1.0
        idx = np.searchsorted(arr, val) - 1
        idx = max(0, min(idx, n - 2))
        log_val = np.log10(max(val, 1e-30))
        log_lo = np.log10(max(arr[idx], 1e-30))
        log_hi = np.log10(max(arr[idx+1], 1e-30))
        if abs(log_hi - log_lo) > 1e-30:
            frac = (log_val - log_lo) / (log_hi - log_lo)
        else:
            frac = 0.0
        return idx, max(0.0, min(1.0, frac))

    i_d, f_n = find_idx(den, nH)
    j_z, f_z = find_idx(met, Z_met)

    # Build Y for 4 corners, bilinear interpolate
    corners = []
    i_list = [i_d, min(i_d + 1, table['L'] - 1)]
    j_list = [j_z, min(j_z + 1, table['M'] - 1)]
    for ii in i_list:
        for jj in j_list:
            Y, k_ref, u_ref, L_ref = build_Y_grackle(table, ii, jj)
            corners.append((Y, k_ref, u_ref, L_ref))

    Y00, Y01, Y10, Y11 = [c[0] for c in corners]
    Y_interp = np.zeros(N)
    for k in range(N):
        y0 = Y00[k] + f_z * (Y01[k] - Y00[k])
        y1 = Y10[k] + f_z * (Y11[k] - Y10[k])
        Y_interp[k] = y0 + f_n * (y1 - y0)

    k_ref, u_ref, L_ref = corners[0][1], corners[0][2], corners[0][3]

    Lambda_init = interp_table_3d(table, nH, Z_met, T_init, 'net_rate')
    if abs(Lambda_init) < 1e-30:
        return T_init

    is_cooling = (Lambda_init > 0)

    # Y at T_init
    log_T_init = np.log10(max(T_init, T_arr[0]))
    Y_init = np.interp(log_T_init, log_T, Y_interp)

    Delta_Y = (abs(L_ref) / u_ref) * dt
    Y_target = Y_init - Delta_Y if is_cooling else Y_init + Delta_Y

    # Search for T where Y = Y_target
    idx_init = np.searchsorted(log_T, log_T_init)
    idx_init = max(1, min(idx_init, N - 1))

    if is_cooling:
        for k in range(idx_init, 0, -1):
            y0, y1 = Y_interp[k-1], Y_interp[k]
            if (Y_target - y0) * (Y_target - y1) <= 0:
                if abs(y1 - y0) > 1e-30:
                    frac = (Y_target - y0) / (y1 - y0)
                    frac = max(0.0, min(1.0, frac))
                else:
                    frac = 0.5
                T_new = 10.0 ** (log_T[k-1] + frac * (log_T[k] - log_T[k-1]))
                return max(T_arr[0], min(T_new, T_arr[-1]))
            # Equilibrium crossing check
            if table['net_rate'][i_d, j_z, k] > 0 and table['net_rate'][i_d, j_z, k-1] <= 0:
                Lk = table['net_rate'][i_d, j_z, k]
                Lk1 = table['net_rate'][i_d, j_z, k-1]
                if abs(Lk - Lk1) > 1e-30:
                    fr = Lk / (Lk - Lk1)
                    fr = max(0.0, min(1.0, fr))
                    return 10.0 ** (log_T[k] + fr * (log_T[k-1] - log_T[k]))
                return T_arr[k]
        return T_arr[0]
    else:
        for k in range(idx_init, N - 1):
            y0, y1 = Y_interp[k], Y_interp[k+1]
            if (Y_target - y0) * (Y_target - y1) <= 0:
                if abs(y1 - y0) > 1e-30:
                    frac = (Y_target - y0) / (y1 - y0)
                    frac = max(0.0, min(1.0, frac))
                else:
                    frac = 0.5
                T_new = 10.0 ** (log_T[k] + frac * (log_T[k+1] - log_T[k]))
                return max(T_arr[0], min(T_new, T_arr[-1]))
        return T_arr[-1]


def numerical_integrate_grackle(table, nH, Z_met, T_init, dt_total, nsteps):
    """RK4 numerical integration using Grackle table."""
    dt = dt_total / nsteps
    T = float(T_init)
    T_min = table['temperature'][0]
    T_max = table['temperature'][-1]

    t_arr = [0.0]
    T_arr_out = [T]

    for step in range(nsteps):
        def dTdt(Tv):
            Tv = max(T_min, min(Tv, T_max))
            L_spec = interp_table_3d(table, nH, Z_met, Tv, 'net_rate')
            mu_v = interp_table_3d(table, nH, Z_met, Tv, 'mmw')
            return -(GAMMA_AD - 1.0) * mu_v * m_H / k_B * L_spec

        k1 = dTdt(T) * dt
        k2 = dTdt(T + 0.5 * k1) * dt
        k3 = dTdt(T + 0.5 * k2) * dt
        k4 = dTdt(T + k3) * dt
        T = T + (k1 + 2*k2 + 2*k3 + k4) / 6.0
        T = max(T_min, min(T, T_max))

        t_arr.append((step + 1) * dt)
        T_arr_out.append(T)

    return np.array(t_arr), np.array(T_arr_out)


def test2_grackle(table_path):
    """Test 2: Grackle table-based Y-function validation."""
    print("\n" + "=" * 70)
    print("TEST 2: Grackle table cooling (Y-function validation)")
    print("=" * 70)

    import os
    if not os.path.exists(table_path):
        print(f"  Table not found: {table_path}")
        return None

    table = read_grackle_table(table_path)
    fmt = "multi-z" if table.get('is_multi_z') else "single-z"
    print(f"  Loaded {fmt} table, z={table['z']:.1f}")
    print(f"  Grid: L={table['L']} M={table['M']} N={table['N']}")
    print(f"  Density:     {table['density'][0]:.2e} - {table['density'][-1]:.2e} cm^-3")
    print(f"  Z_solar:     {table['metallicity'][0]:.2e} - {table['metallicity'][-1]:.2e}")
    print(f"  Temperature: {table['temperature'][0]:.2e} - {table['temperature'][-1]:.2e} K")

    # Test conditions (within table range)
    test_cases = [
        (1.0,   1.0,  1e7,  "nH=1, Z=Zsun, T=10^7 K"),
        (0.01,  0.1,  1e6,  "nH=0.01, Z=0.1Zsun, T=10^6 K"),
        (100.0, 1.0,  1e5,  "nH=100, Z=Zsun, T=10^5 K"),
        (1.0,   0.01, 1e8,  "nH=1, Z=0.01Zsun, T=10^8 K"),
    ]

    all_results = []

    for nH, Z_met, T_init, label in test_cases:
        print(f"\n  --- {label} ---")

        # Check table range
        if nH < table['density'][0] or nH > table['density'][-1]:
            print(f"    nH={nH} outside table range, skipping")
            continue

        L_init = interp_table_3d(table, nH, Z_met, T_init, 'net_rate')
        mu_init = interp_table_3d(table, nH, Z_met, T_init, 'mmw')

        if abs(L_init) < 1e-30:
            print(f"    Negligible cooling rate (L_spec={L_init:.3e}), skipping")
            continue

        u_init = k_B * T_init / ((GAMMA_AD - 1.0) * mu_init * m_H)
        t_cool = abs(u_init / L_init)
        t_total = 10.0 * t_cool

        print(f"    L_spec={L_init:.3e} erg/g/s, mu={mu_init:.4f}, t_cool={t_cool:.3e} s ({t_cool/Myr:.2f} Myr)")

        # Reference: RK4
        n_ref = 100000
        print(f"    Computing RK4 reference ({n_ref} steps)...", end=" ", flush=True)
        t_ref, T_ref = numerical_integrate_grackle(table, nH, Z_met, T_init, t_total, n_ref)
        print(f"T_final = {T_ref[-1]:.6e} K")

        # Exact integration with various dt
        dt_factors = [0.01, 0.1, 0.5, 1.0, 2.0, 5.0, 10.0]
        print(f"    {'dt/t_cool':>12s}  {'N_steps':>8s}  {'T_exact':>14s}  {'T_RK4':>14s}  {'rel_err':>12s}")
        print("    " + "-" * 70)

        case_results = []
        for dtf in dt_factors:
            dt_step = dtf * t_cool
            n_steps = max(1, int(np.ceil(t_total / dt_step)))
            T = T_init
            for step in range(n_steps):
                dt_actual = min(dt_step, t_total - step * dt_step)
                if dt_actual <= 0: break
                T = exact_integrate_grackle(table, nH, Z_met, T, dt_actual)

            T_ref_final = T_ref[-1]
            if abs(T_ref_final) > 1.0:
                rel_err = abs(T - T_ref_final) / abs(T_ref_final)
            else:
                rel_err = abs(T - T_ref_final)
            print(f"    {dtf:12.2f}  {n_steps:8d}  {T:14.6e}  {T_ref_final:14.6e}  {rel_err:12.4e}")
            case_results.append((dtf, T, rel_err))

        all_results.append({
            'label': label, 'nH': nH, 'Z_met': Z_met, 'T_init': T_init,
            't_cool': t_cool, 't_ref': t_ref, 'T_ref': T_ref,
            'case_results': case_results
        })

    return table, all_results


# ═══════════════════════════════════════════════════════════════
# PLOTTING
# ═══════════════════════════════════════════════════════════════

def make_plots(test1_data, test2_data, outdir):
    """Generate validation plots."""
    import os
    os.makedirs(outdir, exist_ok=True)

    # ── Test 1: Power-law cooling ──
    if test1_data is not None:
        fig, axes = plt.subplots(1, 2, figsize=(14, 5))

        t_ref = test1_data['t_ref']
        T_ref = test1_data['T_ref']
        tc = test1_data['t_cool']

        # Left: T(t) cooling curves
        ax = axes[0]
        ax.plot(t_ref / tc, T_ref, 'k-', lw=2.5, label='RK4 ref (500k steps)', zorder=10)

        colors = plt.cm.tab10(np.linspace(0, 0.7, len(test1_data['exact'])))
        markers = ['o', 's', '^', 'D', 'v', '<', '>']
        for idx, (dtf, (t_ex, T_ex)) in enumerate(test1_data['exact'].items()):
            ax.plot(t_ex / tc, T_ex, marker=markers[idx % len(markers)],
                    color=colors[idx], ms=5, lw=1, alpha=0.8,
                    label=f'Exact dt={dtf}tc ({len(t_ex)-1} steps)')

        ax.set_yscale('log')
        ax.set_xlabel(r'$t / t_{\rm cool}$', fontsize=12)
        ax.set_ylabel('T [K]', fontsize=12)
        ax.set_title('Piecewise Power-Law Cooling (Townsend 2009)', fontsize=12)
        ax.legend(fontsize=7, loc='lower left')
        ax.set_ylim(bottom=80)
        ax.grid(True, alpha=0.3)

        # Right: Cooling function
        ax2 = axes[1]
        T_plot = test1_data['T_grid']
        LN_plot = test1_data['LN_grid']
        ax2.loglog(T_plot, LN_plot, 'b-', lw=2)
        ax2.set_xlabel('T [K]', fontsize=12)
        ax2.set_ylabel(r'$\Lambda_N$ [erg cm$^3$ s$^{-1}$]', fontsize=12)
        ax2.set_title('Cooling Function', fontsize=12)
        ax2.grid(True, alpha=0.3)
        # Mark segment boundaries
        for Tb in PL_T_BOUNDS[1:-1]:
            ax2.axvline(Tb, color='gray', ls=':', alpha=0.4)

        plt.tight_layout()
        fig.savefig(os.path.join(outdir, 'test1_powerlaw.png'), dpi=150, bbox_inches='tight')
        print(f"\n  Saved: {outdir}/test1_powerlaw.png")
        plt.close(fig)

    # ── Test 2: Grackle table ──
    if test2_data is not None:
        table, all_results = test2_data
        if not all_results:
            print("  No valid test cases for Test 2, skipping plots.")
            return

        n_cases = len(all_results)

        # Plot 1: T(t) for each test case
        fig, axes = plt.subplots(1, n_cases, figsize=(5*n_cases, 5), squeeze=False)
        for idx, res in enumerate(all_results):
            ax = axes[0][idx]
            tc = res['t_cool']
            ax.plot(res['t_ref'] / tc, res['T_ref'], 'k-', lw=2, label='RK4 ref', zorder=10)

            for dtf, T_f, err in res['case_results']:
                ax.axhline(T_f, ls=':', alpha=0.3, color='gray')
                ax.plot(10.0, T_f, 'o', ms=6,
                        label=f'dt={dtf}tc err={err:.1e}')

            ax.set_yscale('log')
            ax.set_xlabel(r'$t / t_{\rm cool}$')
            ax.set_ylabel('T [K]')
            ax.set_title(res['label'], fontsize=10)
            ax.legend(fontsize=6, loc='upper right')
            ax.grid(True, alpha=0.3)

        plt.tight_layout()
        fig.savefig(os.path.join(outdir, 'test2_T_evolution.png'), dpi=150, bbox_inches='tight')
        print(f"  Saved: {outdir}/test2_T_evolution.png")
        plt.close(fig)

        # Plot 2: Error convergence
        fig, ax = plt.subplots(figsize=(8, 5))
        for res in all_results:
            dtfs = [cr[0] for cr in res['case_results']]
            errors = [max(cr[2], 1e-16) for cr in res['case_results']]
            ax.semilogy(dtfs, errors, 'o-', ms=8, lw=2, label=res['label'])

        ax.set_xlabel(r'$\Delta t / t_{\rm cool}$', fontsize=12)
        ax.set_ylabel('Relative Error vs RK4', fontsize=12)
        ax.set_title('Exact Integration Error Convergence', fontsize=12)
        ax.axhline(0.01, color='gray', ls='--', alpha=0.5, label='1% threshold')
        ax.legend(fontsize=8)
        ax.grid(True, alpha=0.3)
        fig.savefig(os.path.join(outdir, 'test2_error_convergence.png'), dpi=150, bbox_inches='tight')
        print(f"  Saved: {outdir}/test2_error_convergence.png")
        plt.close(fig)

        # Plot 3: Cooling curves
        fig, ax = plt.subplots(figsize=(8, 5))
        T_plot = table['temperature']
        for res in all_results:
            nH, Z_met = res['nH'], res['Z_met']
            rates = np.array([interp_table_3d(table, nH, Z_met, T, 'net_rate') for T in T_plot])
            ax.semilogy(np.log10(T_plot), np.abs(rates), lw=2,
                        label=f'nH={nH}, Z={Z_met}Zsun')
            # Mark heating (negative rate) regions
            heating = rates < 0
            if np.any(heating):
                ax.semilogy(np.log10(T_plot[heating]), np.abs(rates[heating]),
                            'x', ms=3, alpha=0.5)

        ax.set_xlabel(r'$\log_{10}(T)$ [K]', fontsize=12)
        ax.set_ylabel(r'$|\Lambda_{\rm spec}|$ [erg/g/s]', fontsize=12)
        ax.set_title(f'Grackle Specific Cooling Rate (z={table["z"]:.1f})', fontsize=12)
        ax.legend(fontsize=8)
        ax.grid(True, alpha=0.3)
        fig.savefig(os.path.join(outdir, 'test2_cooling_curves.png'), dpi=150, bbox_inches='tight')
        print(f"  Saved: {outdir}/test2_cooling_curves.png")
        plt.close(fig)


# ═══════════════════════════════════════════════════════════════
# MAIN
# ═══════════════════════════════════════════════════════════════

if __name__ == '__main__':
    import sys, os

    code_dir = os.path.dirname(os.path.abspath(__file__))
    base_dir = os.path.join(code_dir, '..')

    # Use single-z table for reliability (multi-z also works now)
    grackle_single = os.path.join(base_dir, 'bin', 'grackle_z5.1.bin')
    grackle_multi = os.path.join(base_dir, 'bin', 'grackle_multi_z.bin')

    table_path = grackle_multi if os.path.exists(grackle_multi) else grackle_single
    if len(sys.argv) > 1:
        table_path = sys.argv[1]

    # Run tests
    test1_data = test1_piecewise_powerlaw()
    test2_data = test2_grackle(table_path) if table_path and os.path.exists(table_path) else None

    # Generate plots
    outdir = os.path.join(base_dir, 'misc')
    make_plots(test1_data, test2_data, outdir)

    print("\nAll tests complete.")
