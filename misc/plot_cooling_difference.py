#!/usr/bin/env python3
"""
Compare RAMSES cooling module vs Grackle cooling table on a (nH, T) diagram.
Implements the RAMSES cooling physics in Python for direct comparison.
"""
import struct
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import SymLogNorm
import os

# ============================================================
# Physical constants (same as cooling_module.f90)
# ============================================================
kB = 1.38062e-16      # erg/K
mH = 1.66e-24         # g
eV = 1.6022e-12       # erg
hplanck = 6.6262e-27  # erg s
clight = 2.99793e10   # cm/s
twopi = 6.2831853
X_H = 0.76
Y_He = 0.24
dumfac_ion = 2.0      # dumfac_ion_theuns
dumfac_rec = 0.75     # dumfac_rec_theuns
alpha_spec = 1.0      # J(nu) ~ nu^-alpha
zreioniz = 10.0       # from namelist z_reion
J0min_ref = 2.77168510365299962e-25
aexp_ref = 0.0001

# Courty polynomial coefficients (0:7, 1:6)
# columns: taux_rad(HI,HEI,HEII), heat_rad(HI,HEI,HEII)
coefcourty = np.array([
    [-13.5857,  1.24475,    0.187739, -0.430409, 0.152544,  -0.0246448,  0.00192622, -5.89772e-05],
    [-14.0242,  1.99211,   -0.490766, -0.122646, 0.0776501, -0.0146310,  0.00123335, -3.96066e-05],
    [-15.6627,  0.128240,   1.65633,  -1.23799,  0.372157,  -0.0561687,  0.00422696, -0.000126344],
    [-24.8422,  1.50750,   -0.0699428,-0.308682, 0.122196,  -0.0205179,  0.00163695, -5.08050e-05],
    [-25.0252,  1.79577,   -0.159054, -0.300924, 0.125343,  -0.0214598,  0.00173377, -5.43576e-05],
    [-26.4168,  0.0479454,  1.70948,  -1.26395,  0.378922,  -0.0570957,  0.00428897, -0.000127909],
])  # shape (6, 8)
coef_fit = np.array([20., 20., 20., 20., 20., 20.])
beta_fit = np.array([6, 6, 8, 6, 6, 8])

# Cloudy 07 metal cooling table
temperature_cc07 = np.array([
    3.9684,4.0187,4.0690,4.1194,4.1697,4.2200,4.2703,
    4.3206,4.3709,4.4212,4.4716,4.5219,4.5722,4.6225,
    4.6728,4.7231,4.7734,4.8238,4.8741,4.9244,4.9747,
    5.0250,5.0753,5.1256,5.1760,5.2263,5.2766,5.3269,
    5.3772,5.4275,5.4778,5.5282,5.5785,5.6288,5.6791,
    5.7294,5.7797,5.8300,5.8804,5.9307,5.9810,6.0313,
    6.0816,6.1319,6.1822,6.2326,6.2829,6.3332,6.3835,
    6.4338,6.4841,6.5345,6.5848,6.6351,6.6854,6.7357,
    6.7860,6.8363,6.8867,6.9370,6.9873,7.0376,7.0879,
    7.1382,7.1885,7.2388,7.2892,7.3395,7.3898,7.4401,
    7.4904,7.5407,7.5911,7.6414,7.6917,7.7420,7.7923,
    7.8426,7.8929,7.9433,7.9936,8.0439,8.0942,8.1445,
    8.1948,8.2451,8.2955,8.3458,8.3961,8.4464,8.4967])
excess_cooling_cc07 = np.array([
    -24.9082, -24.9082, -24.5503, -24.0898, -23.5328, -23.0696, -22.7758,
    -22.6175, -22.5266, -22.4379, -22.3371, -22.2289, -22.1181, -22.0078,
    -21.8992, -21.7937, -21.6921, -21.5961, -21.5089, -21.4343, -21.3765,
    -21.3431, -21.3274, -21.3205, -21.3142, -21.3040, -21.2900, -21.2773,
    -21.2791, -21.3181, -21.4006, -21.5045, -21.6059, -21.6676, -21.6877,
    -21.6934, -21.7089, -21.7307, -21.7511, -21.7618, -21.7572, -21.7532,
    -21.7668, -21.7860, -21.8129, -21.8497, -21.9035, -21.9697, -22.0497,
    -22.1327, -22.2220, -22.3057, -22.3850, -22.4467, -22.4939, -22.5205,
    -22.5358, -22.5391, -22.5408, -22.5408, -22.5475, -22.5589, -22.5813,
    -22.6122, -22.6576, -22.7137, -22.7838, -22.8583, -22.9348, -23.0006,
    -23.0547, -23.0886, -23.1101, -23.1139, -23.1147, -23.1048, -23.1017,
    -23.0928, -23.0969, -23.0968, -23.1105, -23.1191, -23.1388, -23.1517,
    -23.1717, -23.1837, -23.1986, -23.2058, -23.2134, -23.2139, -23.2107])

# Courty phi table for UV suppression of metal cooling
z_courty = np.array([
    0.00000,0.04912,0.10060,0.15470,0.21140,0.27090,0.33330,0.39880,
    0.46750,0.53960,0.61520,0.69450,0.77780,0.86510,0.95670,1.05300,
    1.15400,1.25900,1.37000,1.48700,1.60900,1.73700,1.87100,2.01300,
    2.16000,2.31600,2.47900,2.64900,2.82900,3.01700,3.21400,3.42100,
    3.63800,3.86600,4.10500,4.35600,4.61900,4.89500,5.18400,5.48800,
    5.80700,6.14100,6.49200,6.85900,7.24600,7.65000,8.07500,8.52100,
    8.98900,9.50000])
phi_courty = np.array([
    0.0499886,0.0582622,0.0678333,0.0788739,0.0915889,0.1061913,0.1229119,
    0.1419961,0.1637082,0.1883230,0.2161014,0.2473183,0.2822266,0.3210551,
    0.3639784,0.4111301,0.4623273,0.5172858,0.5752659,0.6351540,0.6950232,
    0.7529284,0.8063160,0.8520859,0.8920522,0.9305764,0.9682031,1.0058810,
    1.0444020,1.0848160,1.1282190,1.1745120,1.2226670,1.2723200,1.3231350,
    1.3743020,1.4247480,1.4730590,1.5174060,1.5552610,1.5833640,1.5976390,
    1.5925270,1.5613110,1.4949610,1.3813710,1.2041510,0.9403100,0.5555344,
    0.0000000])


# ============================================================
# RAMSES Courty UV rates
# ============================================================
def courty_rates(z):
    """Compute Courty UV photo-ionization and photo-heating rates at redshift z."""
    aexp = 1.0 / (1.0 + z)
    J0min = J0min_ref / (aexp / aexp_ref)**2

    t_rad = np.zeros(3)  # HI, HeI, HeII
    h_rad = np.zeros(3)

    if z < zreioniz:
        zz = max(z, 1e-15)
        for ispec in range(6):
            hh = sum(coefcourty[ispec, i] * zz**i for i in range(8))
            hhreion = coef_fit[ispec] * (zz / zreioniz)**beta_fit[ispec]
            val = 10.0**(hh - hhreion)
            if ispec < 3:
                t_rad[ispec] = val
            else:
                h_rad[ispec - 3] = val

    # Apply Theuns floor
    for ispec in range(3):
        # taux_rad_theuns floor
        tt_floor = _taux_rad_theuns(ispec, J0min)
        t_rad[ispec] = max(t_rad[ispec], tt_floor)
        ht_floor = _heat_rad_theuns(ispec, J0min)
        h_rad[ispec] = max(h_rad[ispec], ht_floor)

    return t_rad, h_rad


def _taux_rad_theuns(ispec, J0):
    """Theuns photoionization rate [s^-1]."""
    Jloc = max(J0, 1e-29)
    # These are integrals of J(nu)*sigma(nu)/h/nu dnu
    if ispec == 0:    # HI
        return dumfac_ion * 1.26e-11 * Jloc / 1e-22 / (2.0 + alpha_spec) / (1.0 + alpha_spec)
    elif ispec == 1:  # HeI
        return dumfac_ion * 1.48e-11 * Jloc / 1e-22 * 0.553**alpha_spec * \
               (1.66 / (alpha_spec + 1.05) - 0.66 / (alpha_spec + 2.05))
    else:             # HeII
        return dumfac_ion * 3.16e-12 * Jloc / 1e-22 * 0.249**alpha_spec / \
               (2.0 + alpha_spec) / (1.0 + alpha_spec)


def _heat_rad_theuns(ispec, J0):
    """Theuns photo-heating rate [erg s^-1]."""
    Jloc = max(J0, 1e-29)
    if ispec == 0:    # HI
        return 2.91e-1 * Jloc / (2.0 + alpha_spec) / (3.0 + alpha_spec)
    elif ispec == 1:  # HeI
        return 5.84e-1 * Jloc * 0.553**alpha_spec * \
               (1.66 / (alpha_spec + 1.05) - 2.32 / (alpha_spec + 2.05) + 0.66 / (alpha_spec + 3.05))
    else:             # HeII
        return 2.92e-1 * Jloc * 0.249**alpha_spec / (2.0 + alpha_spec) / (3.0 + alpha_spec)


# ============================================================
# H/He atomic rates
# ============================================================
def cool_bre(ispec, T):
    """Bremsstrahlung cooling [erg cm^3/s]."""
    fac = 1.42e-27 * np.sqrt(T) * (1.1 + 0.34 * np.exp(-(5.5 - np.log10(T))**2 / 3.0))
    if ispec == 2:  # HeII
        return 4.0 * fac
    return fac

def cool_exc(ispec, T):
    """Collisional excitation cooling [erg cm^3/s]."""
    T5 = 1e-5 * T
    sq = 1.0 + np.sqrt(T5)
    if ispec == 0:    return 7.50e-19 / sq * np.exp(-118348.0 / T)
    elif ispec == 1:  return 9.10e-27 / sq / T**0.1687 * np.exp(-13179.0 / T)
    else:             return 5.54e-17 / sq / T**0.397 * np.exp(-473638.0 / T)

def cool_rec(ispec, T):
    """Recombination cooling [erg cm^3/s]."""
    T3, T6 = 1e-3 * T, 1e-6 * T
    if ispec == 0:    return 8.70e-27 * np.sqrt(T) / T3**0.2 / (1.0 + T6**0.7)
    elif ispec == 1:  return 1.55e-26 * T**0.3647
    else:             return 3.48e-26 * np.sqrt(T) / T3**0.2 / (1.0 + T6**0.7)

def cool_die(T):
    """Dielectric recombination cooling [erg cm^3/s]."""
    return 1.24e-13 * T**(-1.5) * np.exp(-470000.0 / T) * (1.0 + 0.3 * np.exp(-94000.0 / T))

def cool_ion(ispec, T):
    """Collisional ionization cooling [erg cm^3/s]."""
    T5 = 1e-5 * T
    sq = 1.0 + np.sqrt(T5)
    if ispec == 0:    return dumfac_ion * 1.27e-21 * np.sqrt(T) / sq * np.exp(-157809.1 / T)
    elif ispec == 1:  return dumfac_ion * 9.38e-22 * np.sqrt(T) / sq * np.exp(-285335.4 / T)
    else:             return dumfac_ion * 4.95e-22 * np.sqrt(T) / sq * np.exp(-631515.0 / T)

def taux_rec(ispec, T):
    """Recombination rate [cm^3/s]."""
    T3, T6 = 1e-3 * T, 1e-6 * T
    if ispec == 0:    return dumfac_rec * 8.40e-11 / np.sqrt(T) / T3**0.2 / (1.0 + T6**0.7)
    elif ispec == 1:
        tdie = 1.9e-3 * T**(-1.5) * np.exp(-470000.0 / T) * (1.0 + 0.3 * np.exp(-94000.0 / T))
        return 1.50e-10 / T**0.6353 + tdie
    else:             return 3.36e-10 / np.sqrt(T) / T3**0.2 / (1.0 + T6**0.7)

def taux_ion(ispec, T):
    """Collisional ionization rate [cm^3/s]."""
    T5 = 1e-5 * T
    sq = 1.0 + np.sqrt(T5)
    if ispec == 0:    return dumfac_ion * 5.85e-11 * np.sqrt(T) / sq * np.exp(-157809.1 / T)
    elif ispec == 1:  return dumfac_ion * 2.38e-11 * np.sqrt(T) / sq * np.exp(-285335.4 / T)
    else:             return dumfac_ion * 5.68e-12 * np.sqrt(T) / sq * np.exp(-631515.0 / T)


# ============================================================
# Chemical equilibrium
# ============================================================
def cmp_chem_eq(T, nH, t_rad):
    """Compute ionization equilibrium for H+He at given T, nH, UV rates."""
    xx = 1.0 - Y_He
    yy = Y_He / (1.0 - Y_He) / 4.0

    t_rec_HI = taux_rec(0, T)
    t_rec_HEI = taux_rec(1, T)
    t_rec_HEII = taux_rec(2, T)
    t_ion_HI = taux_ion(0, T)
    t_ion_HEI = taux_ion(1, T)
    t_ion_HEII = taux_ion(2, T)

    n_E = nH  # initial guess
    for _ in range(200):
        n_E_safe = max(n_E, 1e-15 * nH)
        t2_HI = t_ion_HI + t_rad[0] / n_E_safe
        t2_HEI = t_ion_HEI + t_rad[1] / n_E_safe
        t2_HEII = t_ion_HEII + t_rad[2] / n_E_safe

        n_HI = t_rec_HI / (t2_HI + t_rec_HI) * nH
        n_HII = t2_HI / (t2_HI + t_rec_HI) * nH

        x1 = t_rec_HEII * t_rec_HEI + t2_HEI * t_rec_HEII + t2_HEII * t2_HEI
        n_HEIII = yy * t2_HEII * t2_HEI / x1 * nH
        n_HEII = yy * t2_HEI * t_rec_HEII / x1 * nH
        n_HEI = yy * t_rec_HEII * t_rec_HEI / x1 * nH

        n_E_new = n_HII + n_HEII + 2.0 * n_HEIII
        err = abs(n_E - n_E_new) / nH
        n_E = 0.5 * n_E + 0.5 * n_E_new
        if err < 1e-8:
            break

    n_TOT = n_E + n_HI + n_HII + n_HEI + n_HEII + n_HEIII
    mu = nH / xx / n_TOT
    return n_E, n_HI, n_HII, n_HEI, n_HEII, n_HEIII, mu


# ============================================================
# RAMSES net cooling rate computation
# ============================================================
def ramses_net_rate(nH, T, z, Zsolar):
    """Compute RAMSES net cooling rate [erg/cm^3/s] at given (nH, T, z, Z)."""
    aexp = 1.0 / (1.0 + z)
    t_rad, h_rad = courty_rates(z)

    # Chemical equilibrium
    n_E, n_HI, n_HII, n_HEI, n_HEII, n_HEIII, mu = cmp_chem_eq(T, nH, t_rad)

    # H/He cooling per volume [erg/cm^3/s]
    cool_tot = (
        cool_bre(0, T) * n_E * n_HII +
        cool_bre(1, T) * n_E * n_HEII +
        cool_bre(2, T) * n_E * n_HEIII +
        cool_ion(0, T) * n_E * n_HI +
        cool_ion(1, T) * n_E * n_HEI +
        cool_ion(2, T) * n_E * n_HEII +
        cool_rec(0, T) * n_E * n_HII +
        cool_rec(1, T) * n_E * n_HEII +
        cool_rec(2, T) * n_E * n_HEIII +
        cool_die(T) * n_E * n_HEII +
        cool_exc(0, T) * n_E * n_HI +
        cool_exc(1, T) * n_E * n_HEI +
        cool_exc(2, T) * n_E * n_HEII
    )

    # UV photo-heating per volume [erg/cm^3/s]
    heat_tot = (
        h_rad[0] * n_HI +
        h_rad[1] * n_HEI +
        h_rad[2] * n_HEII
    )

    # Compton cooling/heating per volume [erg/cm^3/s]
    cool_compton = 5.406e-36 * T / aexp**4 * n_E
    heat_compton = 5.406e-36 * 2.726 / aexp**5 * n_E

    # Metal cooling [erg/cm^3/s]
    metal = _metal_cooling(T, nH, mu, z, aexp) * Zsolar * nH * nH

    # Net: positive = cooling, negative = heating
    net = cool_tot - heat_tot + (cool_compton - heat_compton) + metal
    return net, mu


def _metal_cooling(T, nH, mu, z, aexp):
    """Metal cooling rate per nH^2 per Z_solar [erg cm^3/s], with UV suppression."""
    lTT = np.log10(T)

    # UV suppression factor (Courty model)
    if z <= 0.0 or z >= z_courty[-1]:
        ux = 0.0
    else:
        iZ = min(max(int(z / z_courty[-1] * 49.0), 1), 49)
        deltaZ = z_courty[iZ] - z_courty[iZ - 1]
        ux = 1e-4 * (phi_courty[iZ] * (z - z_courty[iZ - 1]) / deltaZ +
                      phi_courty[iZ - 1] * (z_courty[iZ] - z) / deltaZ) / nH

    c1, c2, TT0, TTC, alpha1 = 0.4, 10.0, 1e5, 1e6, 0.15
    g_courty = c1 * (T / TT0)**alpha1 + c2 * np.exp(-TTC / T)
    f_courty = 1.0 / (1.0 + ux / max(g_courty, 1e-30))

    if lTT >= temperature_cc07[-1]:
        return 1e-100
    elif lTT >= 1.0:
        lcool1 = -100.0
        if lTT >= temperature_cc07[0]:
            iT = min(max(int((lTT - temperature_cc07[0]) /
                             (temperature_cc07[-1] - temperature_cc07[0]) * 90.0), 0), 89)
            deltaT = temperature_cc07[iT + 1] - temperature_cc07[iT]
            lcool1 = (excess_cooling_cc07[iT + 1] * (lTT - temperature_cc07[iT]) / deltaT +
                      excess_cooling_cc07[iT] * (temperature_cc07[iT + 1] - lTT) / deltaT)

        # Fine structure IR cooling
        lcool2 = -31.522879 + 2.0 * lTT - 20.0 / T - T * 4.342944e-5
        metal_tot = (10.0**lcool1 + 10.0**lcool2) * f_courty
        return metal_tot
    else:
        return 1e-100


# ============================================================
# Read Grackle binary table
# ============================================================
def read_grackle_table(filename):
    with open(filename, 'rb') as f:
        L, M, N = struct.unpack('iii', f.read(12))
        z_val, log_dmin, log_dmax, log_zmin, log_zmax, log_tmin, log_tmax = \
            struct.unpack('7d', f.read(56))
        log_T = np.array(struct.unpack(f'{N}d', f.read(N * 8)))
        total = L * M * N
        net_rate = np.array(struct.unpack(f'{total}d', f.read(total * 8)))
        mmw = np.array(struct.unpack(f'{total}d', f.read(total * 8)))
    net_rate = net_rate.reshape(L, M, N)
    mmw = mmw.reshape(L, M, N)
    density = np.logspace(log_dmin, log_dmax, L)
    metallicity = np.logspace(log_zmin, log_zmax, M)
    temperature = 10.0**log_T
    return {
        'L': L, 'M': M, 'N': N, 'redshift': z_val,
        'density': density, 'metallicity': metallicity, 'temperature': temperature,
        'net_rate': net_rate, 'mmw': mmw, 'log_T': log_T,
    }


# ============================================================
# Main: compute RAMSES rates and compare with Grackle
# ============================================================
if __name__ == '__main__':
    basedir = '/gpfs/kjhan/Hydro/CUBE_HR5/code_cube'
    outdir = os.path.join(basedir, 'test_ksection/run_compare_cooling')
    table_file = os.path.join(outdir, 'grackle_z5.1.bin')

    if not os.path.exists(table_file):
        print(f'Table not found: {table_file}')
        exit(1)

    tab = read_grackle_table(table_file)
    z_val = tab['redshift']
    print(f'Grackle table: L={tab["L"]}, M={tab["M"]}, N={tab["N"]}, z={z_val:.1f}')

    # Choose metallicity: z_ave = 1e-3 Z_solar (from namelist)
    Zsolar = 1e-3
    Z_idx = np.argmin(np.abs(np.log10(tab['metallicity']) - np.log10(Zsolar)))
    Z_actual = tab['metallicity'][Z_idx]
    print(f'Using Z = {Z_actual:.2e} Z_solar (index {Z_idx})')

    # Define the (nH, T) grid — use a subset for computation speed
    nH_grid = tab['density']    # 60 points
    T_grid = tab['temperature']  # 256 points
    Nn, Nt = len(nH_grid), len(T_grid)

    # Grackle net rate [erg/cm^3/s] on (nH, T) — already in the table
    # Negate: Grackle convention (negative=cooling) -> RAMSES (positive=cooling)
    grackle_rate = -tab['net_rate'][:, Z_idx, :]  # shape (Nn, Nt)
    grackle_mu = tab['mmw'][:, Z_idx, :]          # shape (Nn, Nt)

    # Compute RAMSES rate on the same grid
    print(f'Computing RAMSES cooling rates on {Nn}x{Nt} grid...')
    ramses_rate = np.zeros((Nn, Nt))
    ramses_mu = np.zeros((Nn, Nt))

    for i in range(Nn):
        for j in range(Nt):
            nH_val = nH_grid[i]
            T_val = T_grid[j]
            try:
                rate, mu = ramses_net_rate(nH_val, T_val, z_val, Z_actual)
                ramses_rate[i, j] = rate
                ramses_mu[i, j] = mu
            except Exception:
                ramses_rate[i, j] = 0.0
                ramses_mu[i, j] = 0.6

    print('Done computing RAMSES rates.')

    # ============================================================
    # Plot 1: 4-panel comparison on (nH, T) diagram
    # ============================================================
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    log_T = np.log10(T_grid)
    log_nH = np.log10(nH_grid)

    # Panel 1: RAMSES net cooling rate
    ax = axes[0, 0]
    r = ramses_rate.copy()
    vmax = np.max(np.abs(r[r != 0])) if np.any(r != 0) else 1e-20
    im = ax.pcolormesh(log_nH, log_T, r.T, cmap='RdBu_r',
                       norm=SymLogNorm(linthresh=1e-30, vmin=-vmax, vmax=vmax),
                       shading='nearest')
    plt.colorbar(im, ax=ax, label=r'$\Lambda_{\rm net}$ [erg/cm$^3$/s]')
    ax.set_title(f'RAMSES (Courty UV, CC07 metal)', fontsize=11)
    ax.set_xlabel(r'log$_{10}$(n$_H$ [cm$^{-3}$])')
    ax.set_ylabel(r'log$_{10}$(T [K])')
    ax.set_xlim(-8, 4)
    ax.set_ylim(1, 9)

    # Panel 2: Grackle net cooling rate
    ax = axes[0, 1]
    g = grackle_rate.copy()
    im = ax.pcolormesh(log_nH, log_T, g.T, cmap='RdBu_r',
                       norm=SymLogNorm(linthresh=1e-30, vmin=-vmax, vmax=vmax),
                       shading='nearest')
    plt.colorbar(im, ax=ax, label=r'$\Lambda_{\rm net}$ [erg/cm$^3$/s]')
    ax.set_title(f'Grackle (HM2012 UV, non-LTE)', fontsize=11)
    ax.set_xlabel(r'log$_{10}$(n$_H$ [cm$^{-3}$])')
    ax.set_ylabel(r'log$_{10}$(T [K])')
    ax.set_xlim(-8, 4)
    ax.set_ylim(1, 9)

    # Panel 3: log10 ratio |RAMSES/Grackle| — clipped
    ax = axes[1, 0]
    # Avoid division by zero
    eps = 1e-40
    with np.errstate(divide='ignore', invalid='ignore'):
        ratio = np.log10(np.abs(ramses_rate + eps) / np.abs(grackle_rate + eps))
    ratio = np.clip(ratio, -3, 3)
    im = ax.pcolormesh(log_nH, log_T, ratio.T, cmap='coolwarm',
                       vmin=-3, vmax=3, shading='nearest')
    cb = plt.colorbar(im, ax=ax, label=r'log$_{10}$(|$\Lambda_{\rm RAMSES}$| / |$\Lambda_{\rm Grackle}$|)')
    ax.set_title('Cooling Rate Ratio (red = RAMSES stronger)', fontsize=11)
    ax.set_xlabel(r'log$_{10}$(n$_H$ [cm$^{-3}$])')
    ax.set_ylabel(r'log$_{10}$(T [K])')
    ax.set_xlim(-8, 4)
    ax.set_ylim(1, 9)
    # Equilibrium temperature contour (net=0 lines)
    try:
        ax.contour(log_nH, log_T, grackle_rate.T, levels=[0], colors='lime',
                   linewidths=1.5, linestyles='--')
        ax.contour(log_nH, log_T, ramses_rate.T, levels=[0], colors='yellow',
                   linewidths=1.5, linestyles='-')
    except Exception:
        pass

    # Panel 4: mu difference
    ax = axes[1, 1]
    mu_diff = ramses_mu - grackle_mu
    vabs = max(np.max(np.abs(mu_diff)), 0.01)
    im = ax.pcolormesh(log_nH, log_T, mu_diff.T, cmap='coolwarm',
                       vmin=-vabs, vmax=vabs, shading='nearest')
    plt.colorbar(im, ax=ax, label=r'$\mu_{\rm RAMSES} - \mu_{\rm Grackle}$')
    ax.set_title(r'Mean Molecular Weight Difference $\Delta\mu$', fontsize=11)
    ax.set_xlabel(r'log$_{10}$(n$_H$ [cm$^{-3}$])')
    ax.set_ylabel(r'log$_{10}$(T [K])')
    ax.set_xlim(-8, 4)
    ax.set_ylim(1, 9)

    fig.suptitle(f'RAMSES vs Grackle Cooling Comparison  (z = {z_val:.1f}, '
                 f'Z = {Z_actual:.1e} Z$_\\odot$)',
                 fontsize=14, fontweight='bold')
    plt.tight_layout()
    outfile = os.path.join(basedir, 'misc', 'cooling_rate_difference.png')
    plt.savefig(outfile, dpi=150, bbox_inches='tight')
    plt.close()
    print(f'Saved {outfile}')

    # ============================================================
    # Plot 2: Also show at solar metallicity Z=1
    # ============================================================
    Z_idx_sol = np.argmin(np.abs(np.log10(tab['metallicity'])))
    Z_sol = tab['metallicity'][Z_idx_sol]
    print(f'\nAlso computing for Z = {Z_sol:.2f} Z_solar (index {Z_idx_sol})...')

    # Negate: Grackle convention (negative=cooling) -> RAMSES (positive=cooling)
    grackle_rate_sol = -tab['net_rate'][:, Z_idx_sol, :]
    grackle_mu_sol = tab['mmw'][:, Z_idx_sol, :]

    ramses_rate_sol = np.zeros((Nn, Nt))
    ramses_mu_sol = np.zeros((Nn, Nt))
    for i in range(Nn):
        for j in range(Nt):
            try:
                rate, mu = ramses_net_rate(nH_grid[i], T_grid[j], z_val, Z_sol)
                ramses_rate_sol[i, j] = rate
                ramses_mu_sol[i, j] = mu
            except Exception:
                ramses_rate_sol[i, j] = 0.0
                ramses_mu_sol[i, j] = 0.6

    fig, axes = plt.subplots(2, 2, figsize=(16, 12))

    vmax_sol = np.max(np.abs(ramses_rate_sol[ramses_rate_sol != 0]))
    ax = axes[0, 0]
    im = ax.pcolormesh(log_nH, log_T, ramses_rate_sol.T, cmap='RdBu_r',
                       norm=SymLogNorm(linthresh=1e-30, vmin=-vmax_sol, vmax=vmax_sol),
                       shading='nearest')
    plt.colorbar(im, ax=ax, label=r'$\Lambda_{\rm net}$ [erg/cm$^3$/s]')
    ax.set_title('RAMSES (Courty UV, CC07 metal)', fontsize=11)
    ax.set_xlabel(r'log$_{10}$(n$_H$ [cm$^{-3}$])')
    ax.set_ylabel(r'log$_{10}$(T [K])')
    ax.set_xlim(-8, 4)
    ax.set_ylim(1, 9)

    ax = axes[0, 1]
    im = ax.pcolormesh(log_nH, log_T, grackle_rate_sol.T, cmap='RdBu_r',
                       norm=SymLogNorm(linthresh=1e-30, vmin=-vmax_sol, vmax=vmax_sol),
                       shading='nearest')
    plt.colorbar(im, ax=ax, label=r'$\Lambda_{\rm net}$ [erg/cm$^3$/s]')
    ax.set_title('Grackle (HM2012 UV, non-LTE)', fontsize=11)
    ax.set_xlabel(r'log$_{10}$(n$_H$ [cm$^{-3}$])')
    ax.set_ylabel(r'log$_{10}$(T [K])')
    ax.set_xlim(-8, 4)
    ax.set_ylim(1, 9)

    ax = axes[1, 0]
    with np.errstate(divide='ignore', invalid='ignore'):
        ratio_sol = np.log10(np.abs(ramses_rate_sol + eps) / np.abs(grackle_rate_sol + eps))
    ratio_sol = np.clip(ratio_sol, -3, 3)
    im = ax.pcolormesh(log_nH, log_T, ratio_sol.T, cmap='coolwarm',
                       vmin=-3, vmax=3, shading='nearest')
    plt.colorbar(im, ax=ax, label=r'log$_{10}$(|$\Lambda_{\rm RAMSES}$| / |$\Lambda_{\rm Grackle}$|)')
    ax.set_title('Cooling Rate Ratio (red = RAMSES stronger)', fontsize=11)
    ax.set_xlabel(r'log$_{10}$(n$_H$ [cm$^{-3}$])')
    ax.set_ylabel(r'log$_{10}$(T [K])')
    ax.set_xlim(-8, 4)
    ax.set_ylim(1, 9)
    try:
        ax.contour(log_nH, log_T, grackle_rate_sol.T, levels=[0], colors='lime',
                   linewidths=1.5, linestyles='--')
        ax.contour(log_nH, log_T, ramses_rate_sol.T, levels=[0], colors='yellow',
                   linewidths=1.5, linestyles='-')
    except Exception:
        pass

    ax = axes[1, 1]
    mu_diff_sol = ramses_mu_sol - grackle_mu_sol
    vabs_sol = max(np.max(np.abs(mu_diff_sol)), 0.01)
    im = ax.pcolormesh(log_nH, log_T, mu_diff_sol.T, cmap='coolwarm',
                       vmin=-vabs_sol, vmax=vabs_sol, shading='nearest')
    plt.colorbar(im, ax=ax, label=r'$\mu_{\rm RAMSES} - \mu_{\rm Grackle}$')
    ax.set_title(r'Mean Molecular Weight Difference $\Delta\mu$', fontsize=11)
    ax.set_xlabel(r'log$_{10}$(n$_H$ [cm$^{-3}$])')
    ax.set_ylabel(r'log$_{10}$(T [K])')
    ax.set_xlim(-8, 4)
    ax.set_ylim(1, 9)

    fig.suptitle(f'RAMSES vs Grackle Cooling Comparison  (z = {z_val:.1f}, '
                 f'Z = {Z_sol:.2f} Z$_\\odot$)',
                 fontsize=14, fontweight='bold')
    plt.tight_layout()
    outfile2 = os.path.join(basedir, 'misc', 'cooling_rate_difference_Zsolar.png')
    plt.savefig(outfile2, dpi=150, bbox_inches='tight')
    plt.close()
    print(f'Saved {outfile2}')
