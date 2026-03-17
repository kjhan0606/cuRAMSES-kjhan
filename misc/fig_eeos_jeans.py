#!/usr/bin/env python3
"""
Figure: eEOS Polytropic Floor and the Jeans Length
- Panel (a): Temperature vs density (cooling equilibrium + polytropic floor + eEOS floor)
- Panel (b): Jeans length vs density (with/without floor, cell sizes for each AMR level)
- Panel (c): dt/dt_0 ratio showing CFL collapse when eint < 0
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch

plt.rcParams.update({
    'font.size': 10,
    'font.family': 'serif',
    'axes.labelsize': 11,
    'legend.fontsize': 8,
    'xtick.labelsize': 9,
    'ytick.labelsize': 9,
    'figure.dpi': 150,
    'text.usetex': False,
})

# Physical constants (CGS)
k_B = 1.3807e-16       # Boltzmann constant
m_H = 1.6726e-24       # Hydrogen mass
G = 6.674e-8            # Gravitational constant
gamma = 5.0/3.0
mu = 0.6                # mean molecular weight (ionised)
mu_neutral = 1.22       # neutral gas

# Simulation parameters
n_star = 0.1            # cm^-3, star formation threshold
T2_star = 750.0         # K, polytropic normalization temperature (T/mu)
g_star = 5.0/3.0        # polytropic index
boxlen = 142.0          # Mpc/h comoving
h0 = 0.6777
aexp = 0.21             # scale factor at crash epoch (z~3.7)

# Physical box size at aexp
Lbox_phys = aexp * boxlen / h0 * 3.086e24  # cm (physical)
levelmax = 16
dx = np.array([Lbox_phys / 2**ell for ell in range(9, levelmax+1)])  # levels 9-16
dx_pc = dx / 3.086e18   # convert to pc

# Density range
nH = np.logspace(-6, 5, 500)  # cm^-3
rho = nH * m_H / 0.76   # total mass density

# ---- Temperature models ----
# Approximate cooling equilibrium T(n) for z~3.7 (HM12 UV background)
# Low density: photoionisation equilibrium ~1e4 K
# High density: approaches molecular cooling ~100 K without floor
T_eq = np.where(nH < 1e-2,
                1.0e4 * (nH / 1e-2)**(-0.05),   # warm photoionised
                1.0e4 * (nH / 1e-2)**(-0.7))     # cooling toward cold phase
T_eq = np.maximum(T_eq, 50.0)  # minimum ~50 K from CMB at z~3.7

# Polytropic floor
T_poly = T2_star * mu * (nH / n_star)**(g_star - 1)
T_poly_active = np.where(nH >= n_star, T_poly, 0)

# eEOS floor (same polytropic form but applied at all densities via energy floor)
# In code: e_int >= eeos_poly_coeff * rho^alpha
# This corresponds to T >= T2_star * mu * (n/n_star)^(g_star-1) everywhere
T_eeos = T2_star * mu * (nH / n_star)**(g_star - 1)

# Effective temperature with standard polytropic floor (only above n_star)
T_eff_standard = np.where(nH >= n_star, np.maximum(T_eq, T_poly), T_eq)

# Effective temperature with eEOS floor (everywhere)
T_eff_eeos = np.maximum(T_eq, T_eeos)

# ---- Jeans length ----
def jeans_length(T, rho_val):
    """Jeans length in cm"""
    cs2 = gamma * k_B * T / (mu * m_H)
    cs2 = np.maximum(cs2, 1e-10)  # avoid sqrt of negative
    return np.sqrt(np.pi * cs2 / (G * rho_val))

lJ_eq = jeans_length(T_eq, rho)
lJ_standard = jeans_length(T_eff_standard, rho)
lJ_eeos = jeans_length(T_eff_eeos, rho)

# Convert to pc
lJ_eq_pc = lJ_eq / 3.086e18
lJ_standard_pc = lJ_standard / 3.086e18
lJ_eeos_pc = lJ_eeos / 3.086e18

# ---- Figure ----
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5.5))

# === Panel (a): Temperature vs density ===
ax1.loglog(nH, T_eq, 'b-', lw=1.5, label=r'Cooling equilibrium $T_{\rm eq}(n)$', alpha=0.7)
ax1.loglog(nH[nH >= n_star], T_poly[nH >= n_star], 'r--', lw=2,
           label=r'Polytropic floor ($n \geq n_\star$)')
ax1.loglog(nH, T_eeos, 'r-', lw=1.0, alpha=0.4,
           label=r'eEOS floor (all $n$)')
ax1.loglog(nH, T_eff_standard, 'k-', lw=2, label=r'$T_{\rm eff}$ (standard)', zorder=5)
ax1.loglog(nH, T_eff_eeos, color='darkred', ls='-', lw=2,
           label=r'$T_{\rm eff}$ (eEOS floor)', zorder=4)

# Mark n_star
ax1.axvline(n_star, color='gray', ls=':', lw=1, alpha=0.7)
ax1.text(n_star * 1.3, 3e6, r'$n_\star$', fontsize=10, color='gray')

# Shading: region where eEOS floor is active below n_star
mask_active = (nH < n_star) & (T_eeos > T_eq)
if np.any(mask_active):
    ax1.fill_between(nH[mask_active], T_eq[mask_active], T_eeos[mask_active],
                     color='red', alpha=0.08, label='eEOS protection zone')

# Region where crash occurs (low T, moderate n)
crash_n = np.logspace(-3, -1, 50)
crash_T = 50 * np.ones_like(crash_n)
ax1.fill_between(crash_n, 10, 200, color='blue', alpha=0.08)
ax1.text(3e-3, 30, r'$e_{\rm int} < 0$ crash zone', fontsize=8,
         color='blue', style='italic')

ax1.set_xlabel(r'$n_{\rm H}$ [cm$^{-3}$]')
ax1.set_ylabel(r'Temperature [K]')
ax1.set_xlim(1e-6, 1e5)
ax1.set_ylim(10, 1e8)
ax1.legend(loc='upper left', framealpha=0.9, fontsize=8)
ax1.set_title(r'(a) Temperature--Density Relation', fontsize=11)

# === Panel (b): Jeans length vs density with AMR levels ===
ax2.loglog(nH, lJ_eq_pc, 'b-', lw=1.5, alpha=0.7,
           label=r'$\lambda_J$ (cooling eq. only)')
ax2.loglog(nH, lJ_standard_pc, 'k-', lw=2,
           label=r'$\lambda_J$ (standard polytrope)')
ax2.loglog(nH, lJ_eeos_pc, color='darkred', ls='-', lw=2,
           label=r'$\lambda_J$ (eEOS floor)')

# AMR level cell sizes (4 dx for N_J = 4)
N_J = 4
colors_lv = plt.cm.viridis(np.linspace(0.2, 0.9, len(dx)))
for i, (d, ell) in enumerate(zip(dx_pc, range(9, levelmax+1))):
    if ell in [10, 12, 14, 16]:
        ax2.axhline(N_J * d, color=colors_lv[i], ls='--', lw=0.8, alpha=0.6)
        ax2.text(2e-6, N_J * d * 1.15, f'$4\\Delta x_{{\\ell={ell}}}$',
                fontsize=7, color=colors_lv[i])

# Mark n_star
ax2.axvline(n_star, color='gray', ls=':', lw=1, alpha=0.7)
ax2.text(n_star * 1.3, 5e5, r'$n_\star$', fontsize=10, color='gray')

# Annotation: refinement cascade region
ax2.annotate('Refinement\ncascade\nzone',
             xy=(1e-2, 30), fontsize=8, color='blue',
             ha='center', style='italic',
             bbox=dict(boxstyle='round,pad=0.3', facecolor='lightblue', alpha=0.3))

# Annotation: eEOS-protected region
ax2.annotate('eEOS floor\nprevents\ncascade',
             xy=(1e-4, 5e2), fontsize=8, color='darkred',
             ha='center', style='italic',
             bbox=dict(boxstyle='round,pad=0.3', facecolor='lightyellow', alpha=0.5))

ax2.set_xlabel(r'$n_{\rm H}$ [cm$^{-3}$]')
ax2.set_ylabel(r'Jeans length [pc]')
ax2.set_xlim(1e-6, 1e5)
ax2.set_ylim(1, 1e6)
ax2.legend(loc='upper right', framealpha=0.9, fontsize=8)
ax2.set_title(r'(b) Jeans Length and AMR Resolution', fontsize=11)

plt.tight_layout()
plt.savefig('fig_eeos_jeans.pdf', bbox_inches='tight', dpi=300)
plt.savefig('fig_eeos_jeans.png', bbox_inches='tight', dpi=150)
print("Saved fig_eeos_jeans.pdf and fig_eeos_jeans.png")


# ========== Figure 2: dt collapse ==========
fig2, ax3 = plt.subplots(1, 1, figsize=(8, 4.5))

# Read actual dt data from the crash run
# Data extracted from the analysis above
steps_crash = np.arange(980, 1020)
dt_before = np.full(12, 1.76e-4)
dt_collapse = [7.004e-17, 3.730e-10, 5.427e-10, 8.713e-10, 1.003e-9, 1.098e-9,
               1.232e-9, 1.371e-9]
dt_crash = np.concatenate([dt_before, dt_collapse])

# Normalize to dt_0
dt0 = 1.76e-4
dt_ratio = dt_crash / dt0

fine_steps = np.arange(11560, 11560 + len(dt_crash))

ax3.semilogy(fine_steps, dt_ratio, 'r.-', lw=1.5, markersize=4,
             label='Without eEOS floor (crash)')
ax3.axhline(1.0, color='gray', ls='--', lw=0.8)

# Diagnostic run: stable dt
dt_diag_steps = np.arange(11560, 11610)
dt_diag_ratio = np.ones(50) * (1.75e-4 / dt0)  # ~stable
ax3.semilogy(dt_diag_steps, dt_diag_ratio, 'g-', lw=2, alpha=0.8,
             label='With eEOS floor (stable)')

# Mark the crash point
ax3.axvline(11572, color='red', ls=':', lw=1, alpha=0.7)
ax3.text(11573, 1e-6, r'$e_{\rm int} < 0$' + '\n' + r'$\downarrow$',
         fontsize=9, color='red', ha='left')

# Annotation
ax3.annotate(r'$\Delta t / \Delta t_0 \approx 10^{-12}$',
             xy=(11572, 7e-17/dt0), xytext=(11580, 1e-10),
             fontsize=9, color='red',
             arrowprops=dict(arrowstyle='->', color='red', lw=1.5))

ax3.set_xlabel('Fine step number')
ax3.set_ylabel(r'$\Delta t / \Delta t_0$')
ax3.set_xlim(11555, 11600)
ax3.set_ylim(1e-14, 10)
ax3.legend(loc='lower right', fontsize=9)
ax3.set_title(r'CFL Timestep Collapse from Negative Internal Energy', fontsize=11)

plt.tight_layout()
plt.savefig('fig_dt_collapse_eeos.pdf', bbox_inches='tight', dpi=300)
plt.savefig('fig_dt_collapse_eeos.png', bbox_inches='tight', dpi=150)
print("Saved fig_dt_collapse_eeos.pdf and fig_dt_collapse_eeos.png")
