#!/usr/bin/env python3
"""
Generate dark energy transfer function ratio table R_DE(k,a) for the
cuRAMSES linear response DE perturbation implementation.

Uses CAMB Weyl potential method to extract the NEWTONIAN GAUGE effective
R_DE that correctly enters the modified Poisson equation:
    nabla^2 phi = 1.5*Omega_m*a*delta_m * [1 + (Omega_DE(a)*a^3/Omega_m)*R_DE]

Method:
    1. Run CAMB with target cs2_de (clustering DE)
    2. Run CAMB with cs2=1 (smooth DE, reference — no DE clustering)
    3. R_DE_eff = [(Weyl_target/delta_m_target)/(Weyl_smooth/delta_m_smooth) - 1]
                  * Omega_m / (Omega_L * f_de(a) * a^3)

The old method (extracting delta_DE from delta_tot_de) gave synchronous gauge
values, which are ~26x too large and cause over-correction in the N-body code.

Output: ASCII table readable by dark_energy_commons.f90

Usage:
    python generate_de_table.py [options] [-o output_file]

Examples:
    python generate_de_table.py --w0=-0.8 --wa=0 --cs2=0.0 -o de_table_wcdm.dat
    python generate_de_table.py --w0=-0.9 --wa=0.1 --cs2=0.0 -o de_table_thaw.dat
"""

import sys
import argparse
import numpy as np

try:
    import camb
    from camb.dark_energy import DarkEnergyFluid
except ImportError:
    print("ERROR: CAMB not installed. Install with: pip install camb")
    sys.exit(1)


def generate_de_table(
    output_file="de_table.dat",
    # DE parameters
    w0=-0.9,
    wa=0.1,
    cs2_de=0.0,
    # Cosmological parameters (Planck 2018)
    omega_b=0.0486,
    omega_m=0.3089,
    h0=0.6774,
    n_s=0.9667,
    A_s=2.142e-9,
    # Table grid
    nk=200,
    na=50,
    k_min=1.0e-4,      # h/Mpc
    k_max=10.0,         # h/Mpc
    a_min=0.01,         # scale factor (z=99)
    a_max=1.0,          # scale factor (z=0)
):
    """Generate R_DE(k,a) table using Weyl potential method (Newtonian gauge)."""

    print(f"CAMB version: {camb.__version__}")
    print(f"DE parameters: w0={w0}, wa={wa}, cs2_de={cs2_de}")
    print(f"Cosmology: omega_m={omega_m}, omega_b={omega_b}, h0={h0}")
    print(f"Table: {nk} k-points x {na} a-points")
    print(f"Method: Weyl potential (Newtonian gauge effective R_DE)")

    # k and a grids (log-spaced)
    k_arr = np.logspace(np.log10(k_min), np.log10(k_max), nk)  # h/Mpc
    a_arr = np.logspace(np.log10(a_min), np.log10(a_max), na)
    z_arr = 1.0 / a_arr - 1.0  # redshifts

    omega_cdm = omega_m - omega_b
    omega_l = 1.0 - omega_m  # flat universe
    H0_full = h0 * 100  # km/s/Mpc
    print(f"  omega_cdm={omega_cdm:.6f}, omega_l={omega_l:.6f}")

    # Check w range — DarkEnergyFluid requires w >= -1
    w_min = min(w0, w0 + wa)
    if w_min < -1.0:
        print(f"  ERROR: w_min = {w_min:.3f} < -1. DarkEnergyFluid cannot handle phantom.")
        print(f"  For phantom models, use de_linear_test_setup.py with PPF.")
        sys.exit(1)

    q_arr = k_arr * h0  # Mpc^-1 (CAMB uses Mpc^-1, not h/Mpc)
    z_sorted = sorted(z_arr.tolist(), reverse=True)

    # ================================================================
    # Run CAMB: target model (cs2=cs2_de)
    # ================================================================
    def run_camb(cs2):
        pars = camb.CAMBparams()
        pars.set_cosmology(
            H0=H0_full, ombh2=omega_b * h0**2, omch2=omega_cdm * h0**2,
            mnu=0.0,
        )
        pars.InitPower.set_params(ns=n_s, As=A_s)
        pars.DarkEnergy = DarkEnergyFluid(w=w0, wa=wa, cs2=cs2)
        pars.set_matter_power(redshifts=z_sorted, kmax=k_max * 1.5)
        pars.WantTransfer = True
        return camb.get_results(pars)

    print(f"Running CAMB (target cs2={cs2_de})...")
    results_target = run_camb(cs2_de)

    print(f"Running CAMB (reference cs2=1.0)...")
    results_smooth = run_camb(1.0)

    # ================================================================
    # Compute R_DE from Weyl potential ratio
    # ================================================================
    def f_de(a):
        """DE density ratio: rho_DE(a)/rho_DE(1)"""
        return a**(-3*(1+w0+wa)) * np.exp(-3*wa*(1-a))

    R_DE = np.zeros((nk, na))

    for ia, z_val in enumerate(z_arr):
        a_val = a_arr[ia]

        # Get Weyl potential and matter density from target (clustering DE)
        ev_t = results_target.get_redshift_evolution(
            q_arr, [z_val], ['Weyl', 'delta_cdm', 'delta_baryon'])
        weyl_t    = ev_t[:, 0, 0]
        delta_cdm_t = ev_t[:, 0, 1]
        delta_b_t   = ev_t[:, 0, 2]
        delta_m_t = (omega_cdm * delta_cdm_t + omega_b * delta_b_t) / omega_m

        # Get Weyl potential and matter density from reference (smooth DE)
        ev_s = results_smooth.get_redshift_evolution(
            q_arr, [z_val], ['Weyl', 'delta_cdm', 'delta_baryon'])
        weyl_s    = ev_s[:, 0, 0]
        delta_cdm_s = ev_s[:, 0, 1]
        delta_b_s   = ev_s[:, 0, 2]
        delta_m_s = (omega_cdm * delta_cdm_s + omega_b * delta_b_s) / omega_m

        # Effective de_factor from Weyl potential ratio:
        # de_factor_eff = (Weyl_target/delta_m_target) / (Weyl_smooth/delta_m_smooth)
        # This measures how much the gravitational potential per unit matter
        # perturbation changes due to DE clustering.
        omega_de_a = omega_l * f_de(a_val) * a_val**3

        mask = (np.abs(delta_m_t) > 1e-30) & (np.abs(delta_m_s) > 1e-30) & \
               (np.abs(weyl_s) > 1e-30)
        if np.sum(mask) > 0 and omega_de_a > 1e-30:
            pot_per_dm_t = weyl_t[mask] / delta_m_t[mask]
            pot_per_dm_s = weyl_s[mask] / delta_m_s[mask]
            de_factor_eff = pot_per_dm_t / pot_per_dm_s
            R_DE[mask, ia] = (de_factor_eff - 1.0) * omega_m / omega_de_a

    # ================================================================
    # Write output table (same format as neutrino_table.dat)
    # ================================================================
    print(f"Writing table to {output_file}")
    with open(output_file, "w") as f:
        f.write(f"# Dark energy effective R_DE(k,a) for Newtonian gauge Poisson equation\n")
        f.write(f"# Generated by generate_de_table.py (Weyl method) using CAMB {camb.__version__}\n")
        f.write(f"# DE: w0={w0} wa={wa} cs2_de={cs2_de}\n")
        f.write(f"# Cosmology: omega_m={omega_m} omega_b={omega_b} h0={h0}\n")
        f.write(f"# Method: Newtonian gauge from Weyl potential (NOT sync gauge delta_tot_de)\n")
        f.write(f"# nk  na\n")
        f.write(f"{nk}  {na}\n")
        f.write(f"# k values (h/Mpc), log-spaced\n")
        for i, k in enumerate(k_arr):
            f.write(f"{k:.6e}")
            if (i + 1) % 10 == 0 or i == nk - 1:
                f.write("\n")
            else:
                f.write("  ")
        f.write(f"# a values (scale factor), log-spaced\n")
        for i, a in enumerate(a_arr):
            f.write(f"{a:.6e}")
            if (i + 1) % 10 == 0 or i == na - 1:
                f.write("\n")
            else:
                f.write("  ")
        f.write(f"# R_DE(k,a) matrix: nk rows x na columns\n")
        for ik in range(nk):
            for ia in range(na):
                f.write(f"{R_DE[ik, ia]:.6e}")
                if (ia + 1) % 10 == 0 or ia == na - 1:
                    f.write("\n")
                else:
                    f.write("  ")

    print(f"Done. Table shape: {nk} x {na}")
    print(f"k range: [{k_arr[0]:.2e}, {k_arr[-1]:.2e}] h/Mpc")
    print(f"a range: [{a_arr[0]:.4f}, {a_arr[-1]:.4f}]")
    print(f"R_DE range: [{R_DE.min():.6f}, {R_DE.max():.6f}]")

    # Diagnostic: print de_factor at key epochs
    print("\nde_factor at reference epochs:")
    for z_check in [5.0, 2.0, 1.0, 0.5, 0.0]:
        a_check = 1.0 / (1.0 + z_check)
        ia_check = np.argmin(np.abs(a_arr - a_check))
        omega_de_a_check = omega_l * f_de(a_arr[ia_check]) * a_arr[ia_check]**3
        R_mid = np.median(R_DE[:, ia_check])
        de_f = 1.0 + omega_de_a_check / omega_m * R_mid
        print(f"  z={z_check:.1f}: R_DE={R_mid:.6f}, "
              f"omega_DE_a/omega_m={omega_de_a_check/omega_m:.4f}, "
              f"de_factor={de_f:.6f}")

    # Verify: for cs2_de=0, R_DE should be k-independent at each a
    if cs2_de == 0.0:
        for ia in range(na):
            vals = R_DE[:, ia]
            mask_v = np.abs(vals) > 1e-10
            if np.sum(mask_v) > 2:
                std_k = np.std(vals[mask_v])
                mean_k = np.mean(np.abs(vals[mask_v]))
                if mean_k > 1e-10 and std_k / mean_k > 0.05:
                    print(f"  WARNING: cs2=0 but R_DE varies with k at a={a_arr[ia]:.4f}: "
                          f"std/mean={std_k/mean_k:.4f}")
        print("  cs2=0 check: R_DE should be approximately k-independent")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate DE transfer function table")
    parser.add_argument("--w0", type=float, default=-0.9, help="DE EOS w0 (default: -0.9)")
    parser.add_argument("--wa", type=float, default=0.1, help="DE EOS wa (default: 0.1)")
    parser.add_argument("--cs2", type=float, default=0.0, help="DE sound speed squared (default: 0.0)")
    parser.add_argument("--omega_m", type=float, default=0.3089, help="Omega_m (default: 0.3089)")
    parser.add_argument("--omega_b", type=float, default=0.0486, help="Omega_b (default: 0.0486)")
    parser.add_argument("--h0", type=float, default=0.6774, help="h (default: 0.6774)")
    parser.add_argument("-o", "--output", type=str, default="de_table.dat", help="Output file")
    parser.add_argument("--nk", type=int, default=200, help="Number of k points")
    parser.add_argument("--na", type=int, default=50, help="Number of a points")
    args = parser.parse_args()

    generate_de_table(
        output_file=args.output,
        w0=args.w0, wa=args.wa, cs2_de=args.cs2,
        omega_m=args.omega_m, omega_b=args.omega_b, h0=args.h0,
        nk=args.nk, na=args.na,
    )
