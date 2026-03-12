#!/usr/bin/env python3
"""
Cooling Solver Validation: Townsend (2009) exact vs RAMSES N-R subcycling.

Reads the multi-z Grackle binary cooling table and compares three solvers:
  1. RK4 reference  — adaptive-step 4th-order Runge-Kutta (gold standard)
  2. Townsend exact  — Phi-function exact integration (dt-independent)
  3. RAMSES N-R      — semi-implicit subcycling with varmax parameter

Generates four publication-quality figures:
  fig_cooling_dtindep.pdf      — T_final vs dt (dt-independence)
  fig_cooling_error.pdf        — relative error vs varmax
  fig_cooling_nsub.pdf         — number of subcycles vs varmax
  fig_cooling_trajectory.pdf   — T(t) cooling trajectories

Usage:
  python3 test_exact_cooling.py --table ../bin/grackle_multi_z.bin --redshift 5.0
"""
import argparse
import struct
import sys
import os

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import LogLocator

# ============================================================
# Physical constants (CGS)
# ============================================================
kB   = 1.3806504e-16   # erg/K
mH   = 1.6726219e-24   # g
X_H  = 0.76
Y_He = 0.24
Z_SOLAR = 0.02
GAMMA = 5.0 / 3.0
yr_s = 3.15576e7       # seconds per year
Myr  = 3.15576e13      # seconds per Myr


# ============================================================
# GrackleTable: multi-z reader with z interpolation
# ============================================================
class GrackleTable:
    """Read multi-z or single-z Grackle binary cooling table."""

    def __init__(self, filename):
        self.filename = filename
        self._read(filename)

    def _read(self, filename):
        with open(filename, 'rb') as f:
            first_int = struct.unpack('i', f.read(4))[0]

        # Detect multi-z vs single-z
        is_multi_z = False
        if first_int <= 200 and first_int not in (60, 40):
            with open(filename, 'rb') as f:
                nz_cand = struct.unpack('i', f.read(4))[0]
                L_cand, M_cand, N_cand = struct.unpack('iii', f.read(12))
                if L_cand >= 10 and M_cand >= 10 and N_cand >= 10 and nz_cand != L_cand:
                    is_multi_z = True

        if is_multi_z:
            self._read_multi_z(filename)
        else:
            self._read_single_z(filename)

        self.is_multi_z = is_multi_z
        self.temperature = 10.0 ** self.log_T
        self.density = np.logspace(self.log_dmin, self.log_dmax, self.L)
        self.metallicity = np.logspace(self.log_zmin, self.log_zmax, self.M)

    def _read_multi_z(self, filename):
        with open(filename, 'rb') as f:
            self.nz = struct.unpack('i', f.read(4))[0]
            self.L, self.M, self.N = struct.unpack('iii', f.read(12))
            self.z_values = np.frombuffer(f.read(self.nz * 8), dtype=np.float64).copy()
            self.log_dmin, self.log_dmax, self.log_zmin, self.log_zmax, \
                self.log_tmin, self.log_tmax, _ = struct.unpack('7d', f.read(56))
            self.log_T = np.frombuffer(f.read(self.N * 8), dtype=np.float64).copy()

            sz = self.L * self.M * self.N
            self.all_net_rate = []
            self.all_mmw = []
            for iz in range(self.nz):
                raw_rate = np.frombuffer(f.read(sz * 8), dtype=np.float64).copy().reshape(self.L, self.M, self.N)
                raw_mmw = np.frombuffer(f.read(sz * 8), dtype=np.float64).copy().reshape(self.L, self.M, self.N)
                self.all_net_rate.append(raw_rate)
                self.all_mmw.append(raw_mmw)

    def _read_single_z(self, filename):
        with open(filename, 'rb') as f:
            self.L, self.M, self.N = struct.unpack('iii', f.read(12))
            z_val, self.log_dmin, self.log_dmax, self.log_zmin, self.log_zmax, \
                self.log_tmin, self.log_tmax = struct.unpack('7d', f.read(56))
            self.log_T = np.frombuffer(f.read(self.N * 8), dtype=np.float64).copy()
            sz = self.L * self.M * self.N
            raw_rate = np.frombuffer(f.read(sz * 8), dtype=np.float64).copy().reshape(self.L, self.M, self.N)
            raw_mmw = np.frombuffer(f.read(sz * 8), dtype=np.float64).copy().reshape(self.L, self.M, self.N)
            self.nz = 1
            self.z_values = np.array([z_val])
            self.all_net_rate = [raw_rate]
            self.all_mmw = [raw_mmw]

    def interp_z(self, z_target):
        """Interpolate tables to a target redshift. Returns (net_rate, mmw) arrays.

        Convention: net_rate is converted to erg/g/s and negated so positive = cooling.
        """
        if self.nz == 1:
            iz0, wz = 0, 0.0
        else:
            # Find bracketing z indices
            if z_target <= self.z_values[0]:
                iz0, wz = 0, 0.0
            elif z_target >= self.z_values[-1]:
                iz0, wz = self.nz - 2, 1.0
            else:
                iz0 = np.searchsorted(self.z_values, z_target) - 1
                iz0 = max(0, min(iz0, self.nz - 2))
                wz = (z_target - self.z_values[iz0]) / \
                     (self.z_values[iz0 + 1] - self.z_values[iz0])

        raw_rate = (1 - wz) * self.all_net_rate[iz0]
        raw_mmw = (1 - wz) * self.all_mmw[iz0]
        if self.nz > 1 and wz > 0:
            raw_rate = raw_rate + wz * self.all_net_rate[iz0 + 1]
            raw_mmw = raw_mmw + wz * self.all_mmw[iz0 + 1]

        # Convert to erg/g/s: divide by rho, negate (positive = cooling)
        net_rate = np.zeros_like(raw_rate)
        density = np.logspace(self.log_dmin, self.log_dmax, self.L)
        metallicity = np.logspace(self.log_zmin, self.log_zmax, self.M)
        for i in range(self.L):
            for j in range(self.M):
                Z_metal = metallicity[j] * Z_SOLAR
                Z_metal = max(0.0, min(1.0, Z_metal))
                f_scale = (1.0 - Z_metal) / (X_H + Y_He)
                mf_X = X_H * f_scale
                rho = density[i] * mH / mf_X
                net_rate[i, j, :] = -raw_rate[i, j, :] / rho

        return net_rate, raw_mmw

    def get_table_at_z(self, z_target):
        """Return a dict with interpolated table at z_target."""
        net_rate, mmw = self.interp_z(z_target)
        return {
            'L': self.L, 'M': self.M, 'N': self.N,
            'density': self.density, 'metallicity': self.metallicity,
            'temperature': self.temperature, 'log_T': self.log_T,
            'net_rate': net_rate, 'mmw': mmw, 'z': z_target,
        }


# ============================================================
# Interpolation helper
# ============================================================
def interp_table_3d(tab, nH_val, Z_val, T_val, field='net_rate'):
    """Trilinear interpolation of field at (nH, Z, T). Scalar T only."""
    data = tab[field]
    den, met, temp = tab['density'], tab['metallicity'], tab['temperature']

    def find_idx(arr, val):
        n = len(arr)
        if val <= arr[0]: return 0, 0.0
        if val >= arr[-1]: return n - 2, 1.0
        idx = np.searchsorted(arr, val) - 1
        idx = max(0, min(idx, n - 2))
        log_val = np.log10(max(val, 1e-30))
        log_lo = np.log10(max(arr[idx], 1e-30))
        log_hi = np.log10(max(arr[idx + 1], 1e-30))
        if abs(log_hi - log_lo) > 1e-30:
            frac = (log_val - log_lo) / (log_hi - log_lo)
        else:
            frac = 0.0
        return idx, max(0.0, min(1.0, frac))

    i, fn = find_idx(den, nH_val)
    j, fz = find_idx(met, Z_val)
    k, ft = find_idx(temp, T_val)

    c000 = data[i,   j,   k  ];  c001 = data[i,   j,   k+1]
    c010 = data[i,   j+1, k  ];  c011 = data[i,   j+1, k+1]
    c100 = data[i+1, j,   k  ];  c101 = data[i+1, j,   k+1]
    c110 = data[i+1, j+1, k  ];  c111 = data[i+1, j+1, k+1]

    c00 = c000 + ft * (c001 - c000);  c01 = c010 + ft * (c011 - c010)
    c10 = c100 + ft * (c101 - c100);  c11 = c110 + ft * (c111 - c110)
    c0 = c00 + fz * (c01 - c00);  c1 = c10 + fz * (c11 - c10)
    return c0 + fn * (c1 - c0)


# ============================================================
# CoolingProfile: 1D Lambda(T) and mu(T) at fixed (nH, Z)
# ============================================================
class CoolingProfile:
    """Pre-computed 1D profiles at fixed (nH, Z)."""

    def __init__(self, tab, nH, Z, NT=4096):
        self.nH = nH
        self.Z = Z
        self.NT = NT

        log_T_min = tab['log_T'][0]
        log_T_max = tab['log_T'][-1]
        self.log_T_grid = np.linspace(log_T_min, log_T_max, NT)
        self.T_grid = 10.0 ** self.log_T_grid
        self.dlogT = (log_T_max - log_T_min) / (NT - 1)

        # Vectorized build using trilinear interp
        self.Lambda = np.array([
            interp_table_3d(tab, nH, Z, T, 'net_rate') for T in self.T_grid
        ])
        self.mu = np.array([
            interp_table_3d(tab, nH, Z, T, 'mmw') for T in self.T_grid
        ])

    def get_Lambda(self, T):
        log_T = np.log10(max(T, 10.0))
        return float(np.interp(log_T, self.log_T_grid, self.Lambda))

    def get_mu(self, T):
        log_T = np.log10(max(T, 10.0))
        return float(np.interp(log_T, self.log_T_grid, self.mu))

    def get_Lambda_and_deriv(self, T):
        log_T = np.log10(max(T, 10.0))
        L = float(np.interp(log_T, self.log_T_grid, self.Lambda))
        Lp = float(np.interp(log_T + self.dlogT, self.log_T_grid, self.Lambda))
        Lm = float(np.interp(log_T - self.dlogT, self.log_T_grid, self.Lambda))
        Tp = 10.0 ** (log_T + self.dlogT)
        Tm = 10.0 ** (log_T - self.dlogT)
        dLdT = (Lp - Lm) / (Tp - Tm)
        return L, dLdT

    def find_equilibria(self):
        L = self.Lambda
        equil = []
        for j in range(len(L) - 1):
            if L[j] * L[j + 1] < 0:
                frac = abs(L[j]) / (abs(L[j]) + abs(L[j + 1]))
                log_T_eq = self.log_T_grid[j] + frac * (self.log_T_grid[j + 1] - self.log_T_grid[j])
                equil.append(10.0 ** log_T_eq)
        return equil

    def cooling_time(self, T):
        """Cooling timescale t_cool = |T / (dT/dt)| [seconds].

        dT/dt = -(gamma-1)*mu*mH/kB * Lambda_spec  (Lambda_spec in erg/g/s).
        No nH factor because Lambda_spec already incorporates density.
        """
        L = self.get_Lambda(T)
        mu = self.get_mu(T)
        if abs(L) < 1e-50:
            return np.inf
        C = mu * mH * (GAMMA - 1.0) / kB
        return abs(T / (C * L))


# ============================================================
# Townsend (2009) Exact Integration Solver
# ============================================================
class TownsendSolver:
    """Phi-function exact integration for isochoric cooling/heating."""

    def __init__(self, prof):
        self.prof = prof
        self.T_eq_list = prof.find_equilibria()
        self._build_Phi_table()

    def _build_Phi_table(self):
        NT = self.prof.NT
        T = self.prof.T_grid
        L = self.prof.Lambda
        self.Phi = np.zeros(NT)
        for j in range(1, NT):
            Lavg = 0.5 * (L[j - 1] + L[j])
            dT = T[j] - T[j - 1]
            if abs(Lavg) > 1e-50:
                self.Phi[j] = self.Phi[j - 1] + dT / Lavg
            else:
                self.Phi[j] = self.Phi[j - 1]

    def solve(self, T_init, dt):
        prof = self.prof
        L0 = prof.get_Lambda(T_init)
        mu0 = prof.get_mu(T_init)
        if abs(L0) < 1e-50:
            return T_init

        # dT/dt = -C * Lambda_spec; Lambda_spec in erg/g/s (no nH needed)
        C = mu0 * mH * (GAMMA - 1.0) / kB
        t_cool = abs(T_init / (C * L0))

        # Short timestep: 2nd-order Taylor
        if dt < 1e-4 * t_cool:
            _, dLdT = prof.get_Lambda_and_deriv(T_init)
            dT = -C * L0 * dt + 0.5 * C**2 * L0 * dLdT * dt**2
            return max(T_init + dT, 10.0)

        # Phi-table integration
        Phi_init = float(np.interp(np.log10(max(T_init, 10.0)),
                                   prof.log_T_grid, self.Phi))
        Phi_target = Phi_init - C * dt

        direction = 'down' if L0 > 0 else 'up'
        T_eq = self._find_nearest_eq(T_init, direction)
        return self._inverse_Phi(Phi_target, T_init, direction, T_eq)

    def _find_nearest_eq(self, T_init, direction):
        candidates = [Teq for Teq in self.T_eq_list
                      if (direction == 'down' and Teq < T_init) or
                         (direction == 'up' and Teq > T_init)]
        if not candidates:
            return None
        return max(candidates) if direction == 'down' else min(candidates)

    def _inverse_Phi(self, Phi_target, T_init, direction, T_eq):
        NT = self.prof.NT
        Phi = self.Phi
        T = self.prof.T_grid
        logT = self.prof.log_T_grid

        idx = np.searchsorted(logT, np.log10(max(T_init, 10.0)))
        idx = np.clip(idx, 1, NT - 2)

        if direction == 'down':
            for j in range(idx, 0, -1):
                if T_eq is not None and T[j] <= T_eq:
                    return T_eq
                if (Phi[j] - Phi_target) * (Phi[j - 1] - Phi_target) <= 0:
                    dPhi = Phi[j] - Phi[j - 1]
                    if abs(dPhi) > 1e-50:
                        frac = np.clip((Phi_target - Phi[j - 1]) / dPhi, 0, 1)
                        return 10.0 ** (logT[j - 1] + frac * (logT[j] - logT[j - 1]))
                    return T[j]
            return T_eq if T_eq is not None else T[0]
        else:
            for j in range(idx, NT - 1):
                if T_eq is not None and T[j] >= T_eq:
                    return T_eq
                if (Phi[j] - Phi_target) * (Phi[j + 1] - Phi_target) <= 0:
                    dPhi = Phi[j + 1] - Phi[j]
                    if abs(dPhi) > 1e-50:
                        frac = np.clip((Phi_target - Phi[j]) / dPhi, 0, 1)
                        return 10.0 ** (logT[j] + frac * (logT[j + 1] - logT[j]))
                    return T[j]
            return T_eq if T_eq is not None else T[-1]


# ============================================================
# RAMSES Semi-Implicit Subcycling Solver
# ============================================================
def ramses_nr_solve(prof, T_init, dt, varmax=4.0, iter_max=5000):
    """
    RAMSES N-R subcycling solver with configurable varmax.

    Fortran code (cooling_module.f90:630):
        wcool = MAX(abs(lambda)/tau * varmax, wmax, -lambda_prime * varmax)
        tau_new = tau * (1 + L'/w - L/(tau*w)) / (1 + L'/w)
        time += 1/wcool

    Returns (T_final, n_subcycles).
    """
    nH = prof.nH
    precoeff = 2.0 * X_H / (3.0 * kB)
    time_max = dt * precoeff * nH
    if time_max <= 0:
        return T_init, 0

    wmax = 1.0 / time_max
    mu_init = prof.get_mu(T_init)
    tau = T_init / mu_init
    # Convert Lambda_spec [erg/g/s] to solver's lambda:
    #   Solver: d(tau)/d(time) = -lam, time = t*2*X*nH/(3*kB), tau = T/mu
    #   Physics: dT/dt = -(gamma-1)*mu*mH/kB * Lambda_spec
    #   => lam = 3*(gamma-1)*mH / (2*X_H*nH) * Lambda_spec
    conv = 3.0 * (GAMMA - 1.0) * mH / (2.0 * X_H * nH)

    time_acc = 0.0
    tau_old = tau
    time_old = 0.0
    n_sub = 0

    for it in range(iter_max):
        T_curr = tau * mu_init
        Lambda_curr, dLdT = prof.get_Lambda_and_deriv(T_curr)

        lam = conv * Lambda_curr
        lam_prime = conv * dLdT * mu_init

        wcool = max(abs(lam) / max(tau, 1e-30) * varmax, wmax)
        if lam_prime < 0:
            wcool = max(wcool, -lam_prime * varmax)

        tau_old = tau
        time_old = time_acc
        denom = 1.0 + lam_prime / wcool
        if abs(denom) < 1e-30:
            denom = 1e-30

        tau = tau * (1.0 + lam_prime / wcool - lam / (tau * wcool)) / denom
        time_acc += 1.0 / wcool
        n_sub += 1
        if tau < 1.0:
            tau = 1.0
        if time_acc >= time_max:
            break

    dt_sub = time_acc - time_old
    if dt_sub > 1e-50:
        tau_final = tau * (time_max - time_old) / dt_sub + \
                    tau_old * (time_acc - time_max) / dt_sub
    else:
        tau_final = tau
    return max(tau_final * mu_init, 1.0), n_sub


# ============================================================
# Reference solver: adaptive RK4
# ============================================================
def rk4_solve(prof, T_init, dt, tol=1e-8):
    """Adaptive RK4 reference solver. Returns T_final."""
    mu = prof.get_mu(T_init)
    C = mu * mH * (GAMMA - 1.0) / kB  # No nH; Lambda_spec in erg/g/s
    T_floor = 10.0
    T_eq_list = prof.find_equilibria()
    L0 = prof.get_Lambda(T_init)

    if abs(L0) < 1e-50:
        return T_init

    if L0 > 0:
        candidates = [Teq for Teq in T_eq_list if Teq < T_init]
        T_eq = max(candidates) if candidates else None
    else:
        candidates = [Teq for Teq in T_eq_list if Teq > T_init]
        T_eq = min(candidates) if candidates else None

    T = T_init
    t = 0.0
    t_cool = abs(T_init / (C * L0))
    h = min(dt * 1e-3, 0.01 * t_cool)
    h = max(h, dt * 1e-14)

    for _ in range(500000):
        if t >= dt:
            break
        if t + h > dt:
            h = dt - t
        if h <= 0:
            break

        k1 = -C * prof.get_Lambda(max(T, T_floor))
        k2 = -C * prof.get_Lambda(max(T + 0.5 * h * k1, T_floor))
        k3 = -C * prof.get_Lambda(max(T + 0.5 * h * k2, T_floor))
        k4 = -C * prof.get_Lambda(max(T + h * k3, T_floor))
        T_full = T + h * (k1 + 2 * k2 + 2 * k3 + k4) / 6.0

        hh = 0.5 * h
        k1a = k1
        k2a = -C * prof.get_Lambda(max(T + 0.5 * hh * k1a, T_floor))
        k3a = -C * prof.get_Lambda(max(T + 0.5 * hh * k2a, T_floor))
        k4a = -C * prof.get_Lambda(max(T + hh * k3a, T_floor))
        T_half = T + hh * (k1a + 2 * k2a + 2 * k3a + k4a) / 6.0

        k1b = -C * prof.get_Lambda(max(T_half, T_floor))
        k2b = -C * prof.get_Lambda(max(T_half + 0.5 * hh * k1b, T_floor))
        k3b = -C * prof.get_Lambda(max(T_half + 0.5 * hh * k2b, T_floor))
        k4b = -C * prof.get_Lambda(max(T_half + hh * k3b, T_floor))
        T_dbl = T_half + hh * (k1b + 2 * k2b + 2 * k3b + k4b) / 6.0

        err = abs(T_full - T_dbl) / 15.0
        T_scale = max(abs(T), T_floor)

        if err < tol * T_scale or h < dt * 1e-14:
            T_new = max(T_dbl + (T_dbl - T_full) / 15.0, T_floor)
            if T_eq is not None:
                if (L0 > 0 and T_new <= T_eq) or (L0 < 0 and T_new >= T_eq):
                    T = T_eq
                    break
            T = T_new
            t += h
            if err > 1e-50:
                h_new = h * min((tol * T_scale / err) ** 0.2, 5.0)
            else:
                h_new = h * 5.0
            h = min(h_new, dt - t)
            h = max(h, dt * 1e-14)

            L_curr = prof.get_Lambda(T)
            if abs(L_curr) < 1e-50:
                break
            if (L0 > 0 and L_curr < 0) or (L0 < 0 and L_curr > 0):
                if T_eq is not None:
                    T = T_eq
                break
            if T_eq is not None and abs(T - T_eq) / T_eq < 1e-3:
                T = T_eq
                break
        else:
            h = h * max((tol * T_scale / err) ** 0.25, 0.1)
            h = max(h, dt * 1e-14)

    return max(T, T_floor)


# ============================================================
# RK4 trajectory solver (returns arrays for trajectory plot)
# ============================================================
def rk4_trajectory(prof, T_init, t_total, n_output=500, tol=1e-8):
    """Adaptive RK4 returning (t_arr, T_arr) for plotting."""
    mu = prof.get_mu(T_init)
    C = mu * mH * (GAMMA - 1.0) / kB  # No nH; Lambda_spec in erg/g/s
    T_floor = 10.0
    T_eq_list = prof.find_equilibria()
    L0 = prof.get_Lambda(T_init)

    if abs(L0) < 1e-50:
        return np.array([0, t_total]), np.array([T_init, T_init])

    if L0 > 0:
        candidates = [Teq for Teq in T_eq_list if Teq < T_init]
        T_eq = max(candidates) if candidates else None
    else:
        candidates = [Teq for Teq in T_eq_list if Teq > T_init]
        T_eq = min(candidates) if candidates else None

    dt_out = t_total / n_output
    t_arr = [0.0]
    T_arr = [T_init]

    T = T_init
    t = 0.0
    t_cool = abs(T_init / (C * L0))
    h = min(t_total * 1e-4, 0.01 * t_cool)
    h = max(h, t_total * 1e-14)
    next_out = dt_out

    for _ in range(2000000):
        if t >= t_total:
            break
        if t + h > t_total:
            h = t_total - t
        if h <= 0:
            break

        k1 = -C * prof.get_Lambda(max(T, T_floor))
        k2 = -C * prof.get_Lambda(max(T + 0.5 * h * k1, T_floor))
        k3 = -C * prof.get_Lambda(max(T + 0.5 * h * k2, T_floor))
        k4 = -C * prof.get_Lambda(max(T + h * k3, T_floor))
        T_full = T + h * (k1 + 2 * k2 + 2 * k3 + k4) / 6.0

        hh = 0.5 * h
        k2a = -C * prof.get_Lambda(max(T + 0.5 * hh * k1, T_floor))
        k3a = -C * prof.get_Lambda(max(T + 0.5 * hh * k2a, T_floor))
        k4a = -C * prof.get_Lambda(max(T + hh * k3a, T_floor))
        T_half = T + hh * (k1 + 2 * k2a + 2 * k3a + k4a) / 6.0

        k1b = -C * prof.get_Lambda(max(T_half, T_floor))
        k2b = -C * prof.get_Lambda(max(T_half + 0.5 * hh * k1b, T_floor))
        k3b = -C * prof.get_Lambda(max(T_half + 0.5 * hh * k2b, T_floor))
        k4b = -C * prof.get_Lambda(max(T_half + hh * k3b, T_floor))
        T_dbl = T_half + hh * (k1b + 2 * k2b + 2 * k3b + k4b) / 6.0

        err = abs(T_full - T_dbl) / 15.0
        T_scale = max(abs(T), T_floor)

        if err < tol * T_scale or h < t_total * 1e-14:
            T_new = max(T_dbl + (T_dbl - T_full) / 15.0, T_floor)
            if T_eq is not None:
                if (L0 > 0 and T_new <= T_eq) or (L0 < 0 and T_new >= T_eq):
                    T = T_eq
                    t += h
                    t_arr.append(t)
                    T_arr.append(T)
                    # Fill remaining
                    while t < t_total:
                        t = min(t + dt_out, t_total)
                        t_arr.append(t)
                        T_arr.append(T)
                    break
            T = T_new
            t += h

            if t >= next_out or t >= t_total:
                t_arr.append(t)
                T_arr.append(T)
                next_out = t + dt_out

            if err > 1e-50:
                h_new = h * min((tol * T_scale / err) ** 0.2, 5.0)
            else:
                h_new = h * 5.0
            h = min(h_new, t_total - t)
            h = max(h, t_total * 1e-14)

            L_curr = prof.get_Lambda(T)
            if abs(L_curr) < 1e-50:
                break
            if (L0 > 0 and L_curr < 0) or (L0 < 0 and L_curr > 0):
                if T_eq is not None:
                    T = T_eq
                t_arr.append(t)
                T_arr.append(T)
                break
            if T_eq is not None and abs(T - T_eq) / T_eq < 1e-3:
                T = T_eq
                t_arr.append(t)
                T_arr.append(T)
                break
        else:
            h = h * max((tol * T_scale / err) ** 0.25, 0.1)
            h = max(h, t_total * 1e-14)

    return np.array(t_arr), np.array(T_arr)


# ============================================================
# N-R trajectory solver (for trajectory plot)
# ============================================================
def nr_trajectory(prof, T_init, t_total, n_output=500, varmax=4.0):
    """RAMSES N-R solver returning (t_arr, T_arr) trajectory."""
    dt_out = t_total / n_output
    t_arr = [0.0]
    T_arr = [T_init]
    T = T_init
    t = 0.0
    next_out = dt_out

    while t < t_total:
        dt_step = min(dt_out, t_total - t)
        T, _ = ramses_nr_solve(prof, T, dt_step, varmax=varmax)
        t += dt_step
        t_arr.append(t)
        T_arr.append(T)

    return np.array(t_arr), np.array(T_arr)


# ============================================================
# Townsend trajectory solver
# ============================================================
def townsend_trajectory(solver, T_init, t_total, n_output=500):
    """Townsend exact solver returning (t_arr, T_arr) trajectory."""
    dt_out = t_total / n_output
    t_arr = [0.0]
    T_arr = [T_init]
    T = T_init
    t = 0.0

    while t < t_total:
        dt_step = min(dt_out, t_total - t)
        T = solver.solve(T, dt_step)
        t += dt_step
        t_arr.append(t)
        T_arr.append(T)

    return np.array(t_arr), np.array(T_arr)


# ============================================================
# Test conditions
# ============================================================
TEST_CONDITIONS = [
    # (name, nH, T0, Z, description)
    ('IGM',       1e-4,  1e6,   1e-3, r'IGM: $n_{\rm H}=10^{-4}$, $T=10^6\,$K'),
    ('Filament',  1e-2,  1e5,   1e-2, r'Filament: $n_{\rm H}=10^{-2}$, $T=10^5\,$K'),
    ('ISM',       1e0,   3e4,   1e0,  r'ISM: $n_{\rm H}=1$, $T=3\times10^4\,$K'),
    ('Dense_hot', 1e2,   1e7,   1e0,  r'Dense: $n_{\rm H}=10^2$, $T=10^7\,$K'),
    ('Strong',    1e0,   1e7,   1e0,  r'Strong: $n_{\rm H}=1$, $T=10^7\,$K'),
    ('Dense_ISM', 1e1,   3e4,   1e0,  r'Dense ISM: $n_{\rm H}=10$, $T=3\times10^4\,$K'),
]


# ============================================================
# Figure A: T_final vs dt (dt-independence)
# ============================================================
def fig_dtindep(profiles, townsend_solvers, outdir):
    """2x3 panel: T_final vs dt for each condition."""
    print('\n[Figure A] T_final vs dt (dt-independence)...')

    n_dt = 12
    fig, axes = plt.subplots(2, 3, figsize=(7.0, 4.5))
    axes_flat = axes.flatten()

    for ip, (name, nH, T0, Z, label) in enumerate(TEST_CONDITIONS):
        ax = axes_flat[ip]
        prof = profiles[(nH, Z)]
        ts = townsend_solvers[(nH, Z)]
        t_cool = prof.cooling_time(T0)

        dt_ratios = np.logspace(-2, 3, n_dt)
        dt_vals = dt_ratios * t_cool

        T_rk4 = np.array([rk4_solve(prof, T0, dt) for dt in dt_vals])
        T_tw  = np.array([ts.solve(T0, dt) for dt in dt_vals])
        T_nr4  = np.array([ramses_nr_solve(prof, T0, dt, varmax=4.0)[0] for dt in dt_vals])
        T_nr20 = np.array([ramses_nr_solve(prof, T0, dt, varmax=20.0)[0] for dt in dt_vals])
        T_nr100 = np.array([ramses_nr_solve(prof, T0, dt, varmax=100.0)[0] for dt in dt_vals])

        ax.semilogx(dt_ratios, T_rk4, 'k-', lw=1.5, label='RK4 ref', zorder=5)
        ax.semilogx(dt_ratios, T_tw, 'b-o', lw=1.2, ms=3, label='Townsend', zorder=4)
        ax.semilogx(dt_ratios, T_nr4, 'r--s', lw=1.0, ms=3, label=r'N-R $\Delta_{\rm var}$=4')
        ax.semilogx(dt_ratios, T_nr20, 'g-.^', lw=1.0, ms=3, label=r'N-R $\Delta_{\rm var}$=20')
        ax.semilogx(dt_ratios, T_nr100, 'm:D', lw=1.0, ms=3, label=r'N-R $\Delta_{\rm var}$=100')

        equil = prof.find_equilibria()
        for Teq in equil:
            yvals = np.concatenate([T_rk4, T_tw, T_nr4])
            if min(yvals) * 0.5 < Teq < max(yvals) * 2:
                ax.axhline(Teq, color='orange', ls='-.', lw=0.8, alpha=0.7)

        ax.set_xlabel(r'$\Delta t / t_{\rm cool}$', fontsize=8)
        ax.set_ylabel(r'$T_{\rm final}$ [K]', fontsize=8)
        ax.set_title(label, fontsize=7)
        ax.tick_params(labelsize=7)
        if ip == 0:
            ax.legend(fontsize=5, loc='best')
        ax.grid(True, alpha=0.3, which='both')

        print(f'  {name}: t_cool = {t_cool/yr_s:.2e} yr, T_eq = {equil}')

    # All 6 panels used (no empty legend panel)

    plt.tight_layout()
    outfile = os.path.join(outdir, 'fig_cooling_dtindep.pdf')
    fig.savefig(outfile, dpi=200, bbox_inches='tight')
    outfile_png = outfile.replace('.pdf', '.png')
    fig.savefig(outfile_png, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f'  Saved: {outfile}')


# ============================================================
# Figure B: Relative error vs varmax
# ============================================================
def fig_error_vs_varmax(profiles, townsend_solvers, outdir):
    """2-panel: error vs varmax (full range + practical zoom)."""
    print('\n[Figure B] Relative error vs varmax...')

    varmax_vals = np.logspace(np.log10(0.5), np.log10(200), 30)
    dt_factor = 10.0  # dt = 10 * t_cool (strongly nonlinear regime)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7.0, 3.2))

    colors = plt.cm.tab10(np.linspace(0, 0.5, len(TEST_CONDITIONS)))

    for ip, (name, nH, T0, Z, label) in enumerate(TEST_CONDITIONS):
        prof = profiles[(nH, Z)]
        ts = townsend_solvers[(nH, Z)]
        t_cool = prof.cooling_time(T0)
        dt = dt_factor * t_cool

        T_ref = rk4_solve(prof, T0, dt)
        T_tw = ts.solve(T0, dt)

        errors = []
        for vm in varmax_vals:
            T_nr, _ = ramses_nr_solve(prof, T0, dt, varmax=vm)
            if abs(T_ref) > 1.0:
                err = abs(T_nr - T_ref) / abs(T_ref)
            else:
                err = abs(T_nr - T_ref)
            errors.append(max(err, 1e-15))

        errors = np.array(errors)
        short_label = name

        # Townsend error
        tw_err = abs(T_tw - T_ref) / abs(T_ref) if abs(T_ref) > 1 else abs(T_tw - T_ref)

        for ax in (ax1, ax2):
            ax.loglog(varmax_vals, errors, '-o', color=colors[ip], lw=1.2,
                      ms=3, label=short_label)
            # Townsend as horizontal line
            if tw_err > 1e-15:
                ax.axhline(tw_err, color=colors[ip], ls=':', lw=0.6, alpha=0.5)

        print(f'  {name}: T_ref={T_ref:.4e}, Townsend_err={tw_err:.2e}')

    for ax in (ax1, ax2):
        ax.axhline(0.01, color='gray', ls='--', lw=1, alpha=0.7, label='1%')
        ax.axhline(0.001, color='gray', ls=':', lw=0.8, alpha=0.5, label='0.1%')
        ax.axvline(4.0, color='red', ls='--', lw=0.8, alpha=0.7, label=r'$\Delta_{\rm var}$=4 (default)')
        ax.set_xlabel(r'$\Delta_{\rm var}$', fontsize=9)
        ax.set_ylabel('Relative error vs RK4', fontsize=9)
        ax.tick_params(labelsize=7)
        ax.grid(True, alpha=0.3, which='both')

    ax1.set_title(r'Full range ($\Delta t = 10\,t_{\rm cool}$)', fontsize=9)
    ax1.legend(fontsize=5.5, loc='upper right', ncol=2)
    ax2.set_title('Practical range', fontsize=9)
    ax2.set_xlim(1, 50)
    ax2.set_ylim(1e-5, 1)

    plt.tight_layout()
    outfile = os.path.join(outdir, 'fig_cooling_error.pdf')
    fig.savefig(outfile, dpi=200, bbox_inches='tight')
    fig.savefig(outfile.replace('.pdf', '.png'), dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f'  Saved: {outfile}')


# ============================================================
# Figure C: N_subcycles vs varmax
# ============================================================
def fig_nsub_vs_varmax(profiles, outdir):
    """1 panel: subcycle count vs varmax."""
    print('\n[Figure C] N_subcycles vs varmax...')

    varmax_vals = np.logspace(np.log10(0.5), np.log10(200), 30)
    dt_factor = 10.0

    fig, ax = plt.subplots(figsize=(3.5, 3.2))
    colors = plt.cm.tab10(np.linspace(0, 0.5, len(TEST_CONDITIONS)))

    for ip, (name, nH, T0, Z, label) in enumerate(TEST_CONDITIONS):
        prof = profiles[(nH, Z)]
        t_cool = prof.cooling_time(T0)
        dt = dt_factor * t_cool

        nsubs = []
        for vm in varmax_vals:
            _, nsub = ramses_nr_solve(prof, T0, dt, varmax=vm)
            nsubs.append(nsub)

        nsubs = np.array(nsubs)
        ax.loglog(varmax_vals, nsubs, '-o', color=colors[ip], lw=1.2,
                  ms=3, label=name)

        # Theoretical estimate: N ~ varmax * dt/t_cool
        nsub_theory = varmax_vals * dt_factor
        if ip == 0:
            ax.loglog(varmax_vals, nsub_theory, 'k:', lw=0.8, alpha=0.5,
                      label=r'$\sim \Delta_{\rm var} \cdot \Delta t/t_{\rm cool}$')

        print(f'  {name}: nsub range = [{nsubs.min()}, {nsubs.max()}]')

    ax.axhline(500, color='red', ls='--', lw=1, alpha=0.7, label='iter_max=500')
    ax.axvline(4.0, color='red', ls=':', lw=0.8, alpha=0.5)
    ax.set_xlabel(r'$\Delta_{\rm var}$', fontsize=9)
    ax.set_ylabel(r'$N_{\rm sub}$', fontsize=9)
    ax.set_title(r'Subcycle count ($\Delta t = 10\,t_{\rm cool}$)', fontsize=9)
    ax.legend(fontsize=5.5, loc='upper left')
    ax.tick_params(labelsize=7)
    ax.grid(True, alpha=0.3, which='both')

    plt.tight_layout()
    outfile = os.path.join(outdir, 'fig_cooling_nsub.pdf')
    fig.savefig(outfile, dpi=200, bbox_inches='tight')
    fig.savefig(outfile.replace('.pdf', '.png'), dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f'  Saved: {outfile}')


# ============================================================
# Figure D: T(t) cooling trajectories
# ============================================================
def fig_trajectory(profiles, townsend_solvers, outdir):
    """2x2 panels: T=1.5e4 K cooling trajectories at various n_H.

    Key: N-R uses COARSE steps (dt = 1 t_cool per hydro step) to reveal
    overshoot behavior. Townsend uses same coarse steps to show dt-independence.
    RK4 reference uses fine adaptive steps.
    """
    print('\n[Figure D] T(t) cooling trajectories (T=1.5e4 K, multiple nH)...')

    T0 = 1.5e4
    Z = 1.0  # solar metallicity
    traj_cases = [
        (1e-2, r'$n_{\rm H}=10^{-2}\,{\rm cm}^{-3}$'),
        (1e-1, r'$n_{\rm H}=10^{-1}\,{\rm cm}^{-3}$'),
        (1e0,  r'$n_{\rm H}=1\,{\rm cm}^{-3}$'),
        (1e1,  r'$n_{\rm H}=10\,{\rm cm}^{-3}$'),
    ]

    fig, axes = plt.subplots(2, 2, figsize=(7.0, 5.5))
    axes_flat = axes.flatten()

    for panel_idx, (nH, label) in enumerate(traj_cases):
        ax = axes_flat[panel_idx]
        key = (nH, Z)

        # Get or build profile
        if key in profiles:
            prof = profiles[key]
            ts = townsend_solvers[key]
        else:
            # Need to build from the table used in main
            prof = profiles.get(key)
            ts = townsend_solvers.get(key)
            if prof is None:
                print(f'  WARNING: profile ({nH}, {Z}) not found, skipping')
                continue

        t_cool = prof.cooling_time(T0)
        t_total = 5.0 * t_cool

        # RK4 reference trajectory (fine steps)
        t_rk4, T_rk4 = rk4_trajectory(prof, T0, t_total, n_output=500)

        # Townsend and N-R with coarse steps (markers only, no lines)
        n_coarse = 12
        t_tw, T_tw = townsend_trajectory(ts, T0, t_total, n_output=n_coarse)

        # N-R with coarse steps to reveal overshoot
        t_nr4, T_nr4 = nr_trajectory(prof, T0, t_total, n_output=n_coarse, varmax=4.0)
        t_nr20, T_nr20 = nr_trajectory(prof, T0, t_total, n_output=n_coarse, varmax=20.0)

        ax.semilogy(t_rk4 / t_cool, T_rk4, 'k-', lw=2, label='RK4 ref', zorder=5)
        ax.semilogy(t_tw / t_cool, T_tw, 'bo', ms=4, markerfacecolor='none', lw=0,
                    label=r'Townsend', zorder=4)
        ax.semilogy(t_nr4 / t_cool, T_nr4, 'rs', ms=4, markerfacecolor='none', lw=0,
                    label=r'N-R $\Delta_{\rm var}$=4', zorder=3)
        ax.semilogy(t_nr20 / t_cool, T_nr20, 'g^', ms=4, markerfacecolor='none', lw=0,
                    label=r'N-R $\Delta_{\rm var}$=20', zorder=3)

        # Equilibrium lines
        equil = prof.find_equilibria()
        for Teq in equil:
            yvals = np.concatenate([T_rk4, T_tw])
            if min(yvals) * 0.1 < Teq < max(yvals) * 10:
                ax.axhline(Teq, color='orange', ls='-.', lw=1, alpha=0.7,
                           label=f'$T_{{\\rm eq}}$ = {Teq:.0f} K')

        ax.set_xlabel(r'$t / t_{\rm cool}$', fontsize=9)
        ax.set_ylabel('T [K]', fontsize=9)
        ax.set_title(label, fontsize=8)
        if panel_idx == 0:
            ax.legend(fontsize=5.5, loc='best')
        ax.tick_params(labelsize=7)
        ax.grid(True, alpha=0.3, which='both')
        ax.set_xlim(0, 5)

        print(f'  nH={nH:.0e}: t_cool = {t_cool/yr_s:.2e} yr, equil = {equil}')

    plt.tight_layout()
    outfile = os.path.join(outdir, 'fig_cooling_trajectory.pdf')
    fig.savefig(outfile, dpi=200, bbox_inches='tight')
    fig.savefig(outfile.replace('.pdf', '.png'), dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f'  Saved: {outfile}')


# ============================================================
# Summary table
# ============================================================
def print_summary_table(profiles, townsend_solvers):
    """Print summary: t_cool, T_eq, varmax needed for 1% accuracy."""
    print('\n' + '=' * 90)
    print('SUMMARY TABLE')
    print('=' * 90)
    print(f'  {"Condition":<12s}  {"nH":>8s}  {"T0":>10s}  {"Z":>8s}  '
          f'{"t_cool [yr]":>12s}  {"T_eq [K]":>12s}  '
          f'{"varmax(1%)":>12s}  {"Tw_err":>10s}')
    print(f'  {"-" * 88}')

    dt_factor = 10.0
    varmax_scan = np.logspace(np.log10(0.5), np.log10(200), 60)

    for name, nH, T0, Z, _ in TEST_CONDITIONS:
        prof = profiles[(nH, Z)]
        ts = townsend_solvers[(nH, Z)]
        t_cool = prof.cooling_time(T0)
        dt = dt_factor * t_cool
        equil = prof.find_equilibria()
        T_eq_str = f'{equil[-1]:.2e}' if equil else 'none'

        T_ref = rk4_solve(prof, T0, dt)
        T_tw = ts.solve(T0, dt)
        tw_err = abs(T_tw - T_ref) / abs(T_ref) if abs(T_ref) > 1 else abs(T_tw - T_ref)

        # Find minimum varmax for 1% accuracy
        varmax_1pct = '>200'
        for vm in varmax_scan:
            T_nr, _ = ramses_nr_solve(prof, T0, dt, varmax=vm)
            err = abs(T_nr - T_ref) / abs(T_ref) if abs(T_ref) > 1 else abs(T_nr - T_ref)
            if err < 0.01:
                varmax_1pct = f'{vm:.1f}'
                break

        print(f'  {name:<12s}  {nH:8.1e}  {T0:10.1e}  {Z:8.1e}  '
              f'{t_cool/yr_s:12.2e}  {T_eq_str:>12s}  '
              f'{varmax_1pct:>12s}  {tw_err:10.2e}')


# ============================================================
# Main
# ============================================================
def main():
    parser = argparse.ArgumentParser(description='Cooling solver validation')
    parser.add_argument('--table', '-t', default=None,
                        help='Path to Grackle binary table')
    parser.add_argument('--redshift', '-z', type=float, default=5.0,
                        help='Target redshift (default: 5.0)')
    args = parser.parse_args()

    # Find table
    basedir = os.path.dirname(os.path.abspath(__file__))
    if args.table:
        table_file = args.table
    else:
        table_file = os.path.join(basedir, '..', 'bin', 'grackle_multi_z.bin')
        if not os.path.exists(table_file):
            table_file = os.path.join(basedir, '..', 'bin', 'grackle_z5.1.bin')

    if not os.path.exists(table_file):
        print(f'ERROR: Table not found: {table_file}')
        sys.exit(1)

    print(f'Reading Grackle table: {table_file}')
    gtab = GrackleTable(table_file)
    print(f'  Multi-z: {gtab.is_multi_z}, nz={gtab.nz}, '
          f'z_values={gtab.z_values}')
    print(f'  Grid: L={gtab.L}, M={gtab.M}, N={gtab.N}')
    print(f'  log(nH): [{gtab.log_dmin:.1f}, {gtab.log_dmax:.1f}]')
    print(f'  log(Z):  [{gtab.log_zmin:.1f}, {gtab.log_zmax:.1f}]')
    print(f'  log(T):  [{gtab.log_T[0]:.2f}, {gtab.log_T[-1]:.2f}]')

    z_target = args.redshift
    print(f'\nInterpolating to z = {z_target:.1f}...')
    tab = gtab.get_table_at_z(z_target)

    # Build profiles and solvers
    print('\nBuilding cooling profiles and solvers...')
    profiles = {}
    townsend_solvers = {}
    for name, nH, T0, Z, label in TEST_CONDITIONS:
        key = (nH, Z)
        if key not in profiles:
            profiles[key] = CoolingProfile(tab, nH, Z)
            townsend_solvers[key] = TownsendSolver(profiles[key])
            equil = profiles[key].find_equilibria()
            t_cool = profiles[key].cooling_time(T0)
            L0 = profiles[key].get_Lambda(T0)
            print(f'  ({nH:.0e}, Z={Z:.0e}): T_eq={equil}, '
                  f'Lambda(T0)={L0:.3e}, t_cool={t_cool/yr_s:.2e} yr')

    # Also build profiles for trajectory figure (T=1.5e4, various nH, Z=1)
    traj_nH_list = [1e-2, 1e-1, 1e0, 1e1]
    for nH_traj in traj_nH_list:
        key = (nH_traj, 1.0)
        if key not in profiles:
            profiles[key] = CoolingProfile(tab, nH_traj, 1.0)
            townsend_solvers[key] = TownsendSolver(profiles[key])
            equil = profiles[key].find_equilibria()
            t_cool = profiles[key].cooling_time(1.5e4)
            L0 = profiles[key].get_Lambda(1.5e4)
            print(f'  ({nH_traj:.0e}, Z=1e+00): T_eq={equil}, '
                  f'Lambda(T0=1.5e4)={L0:.3e}, t_cool={t_cool/yr_s:.2e} yr')

    outdir = basedir

    # Generate figures
    fig_dtindep(profiles, townsend_solvers, outdir)
    fig_error_vs_varmax(profiles, townsend_solvers, outdir)
    fig_nsub_vs_varmax(profiles, outdir)
    fig_trajectory(profiles, townsend_solvers, outdir)

    # Summary table
    print_summary_table(profiles, townsend_solvers)

    print('\nAll figures generated successfully.')


if __name__ == '__main__':
    main()
