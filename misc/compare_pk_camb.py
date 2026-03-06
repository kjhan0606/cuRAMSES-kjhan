#!/usr/bin/env python3
"""
Compare simulation P(k) with CAMB linear predictions.

Reads pk_NNNNN.dat from RAMSES output directories, computes CAMB P(k) at
the EXACT output redshift, then produces comparison plots.

Usage:
    python3 compare_pk_camb.py --basedir test_ksection/run_de_pk_test/

Figures produced:
    1. pk_ratio_by_redshift.png  — P_sim/P_CAMB per redshift (4 models per panel)
    2. pk_ratio_by_model.png     — P_sim/P_CAMB per model (6 redshifts per panel)
    3. growth_factor.png         — D^2(z) ratio (sim vs CAMB)
"""

import os
import sys
import argparse
import numpy as np

try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
except ImportError:
    print("ERROR: matplotlib not installed. Install with: pip install matplotlib")
    sys.exit(1)

try:
    import camb
    HAS_CAMB = True
except ImportError:
    HAS_CAMB = False
    print("WARNING: camb not installed. Will use pre-computed P(k) files only.")


# ===================================================================
# Configuration
# ===================================================================
MODELS = ['lcdm', 'wcdm', 'cpl_thaw', 'cpl_phantom']
MODEL_PARAMS = {
    'lcdm':        {'w0': -1.0, 'wa': 0.0},
    'wcdm':        {'w0': -0.8, 'wa': 0.0},
    'cpl_thaw':    {'w0': -0.9, 'wa': 0.2},
    'cpl_phantom': {'w0': -1.1, 'wa': -0.2},
}
MODEL_LABELS = {
    'lcdm': r'$\Lambda$CDM',
    'wcdm': r'$w$CDM ($w_0=-0.8$)',
    'cpl_thaw': r'CPL thaw ($w_0=-0.9, w_a=0.2$)',
    'cpl_phantom': r'CPL phantom ($w_0=-1.1, w_a=-0.2$)',
}
MODEL_COLORS = {
    'lcdm': 'black',
    'wcdm': 'C0',
    'cpl_thaw': 'C1',
    'cpl_phantom': 'C2',
}

# Cosmological parameters (must match de_linear_test_setup.py)
OMEGA_M = 0.3111
OMEGA_B = 0.0490
H0 = 67.66
N_S = 0.9665
A_S = 2.1e-9

Z_OUTPUTS = [5.0, 2.0, 1.0, 0.5, 0.2, 0.0]

# Expected output directories: output_00002 through output_00007
# (output_00001 is the initial snapshot at z_start)
OUTPUT_NUMS = ['00002', '00003', '00004', '00005', '00006', '00007']


def camb_pk_at_z(model_name, z_val):
    """Compute CAMB linear matter P(k) at exact redshift z_val."""
    if not HAS_CAMB:
        return None, None
    # Clamp z to non-negative (sim may overshoot a=1)
    z_use = max(z_val, 0.0)
    params = MODEL_PARAMS[model_name]
    pars = camb.CAMBparams()
    h = H0 / 100.0
    pars.set_cosmology(
        H0=H0,
        ombh2=OMEGA_B * h**2,
        omch2=(OMEGA_M - OMEGA_B) * h**2,
        mnu=0.0,
        omk=0,
    )
    w0, wa = params['w0'], params['wa']
    # Check if w ever crosses or stays below -1 (DarkEnergyFluid can't handle)
    # w(a) = w0 + wa*(1-a); check at a=0 and a=1
    w_min = min(w0, w0 + wa)
    if model_name == 'lcdm' or w_min >= -1.0:
        if model_name == 'lcdm':
            pars.set_dark_energy(w=-1.0, wa=0.0)
        else:
            pars.DarkEnergy = camb.dark_energy.DarkEnergyFluid(
                w=w0, wa=wa, cs2=0.0)
    else:
        # w < -1 at some epoch: use PPF (less accurate but stable)
        pars.set_dark_energy(w=w0, wa=wa, dark_energy_model='ppf')
    pars.InitPower.set_params(As=A_S, ns=N_S)
    pars.set_matter_power(redshifts=[z_use], kmax=50.0, nonlinear=False)
    pars.NonLinear = camb.model.NonLinear_none
    pars.WantTransfer = True
    results = camb.get_results(pars)
    kh, _, pk = results.get_linear_matter_power_spectrum()
    return kh, pk[0]  # pk[0] = P(k) at the single redshift


def read_pk_sim(filepath):
    """Read simulation P(k) file."""
    k, pk, pk_shot, nmodes = [], [], [], []
    aexp = None
    with open(filepath) as f:
        for line in f:
            if line.startswith('#'):
                if 'a_exp' in line:
                    aexp = float(line.split('=')[-1].strip())
                continue
            vals = line.split()
            if len(vals) >= 4:
                k.append(float(vals[0]))
                pk.append(float(vals[1]))
                pk_shot.append(float(vals[2]))
                nmodes.append(float(vals[3]))
    return np.array(k), np.array(pk), np.array(pk_shot), np.array(nmodes), aexp


def read_pk_camb(filepath):
    """Read CAMB linear P(k) file."""
    k, pk = [], []
    with open(filepath) as f:
        for line in f:
            if line.startswith('#'):
                continue
            vals = line.split()
            if len(vals) >= 2:
                k.append(float(vals[0]))
                pk.append(float(vals[1]))
    return np.array(k), np.array(pk)


def interpolate_pk(k_target, k_source, pk_source):
    """Log-log interpolation of P(k)."""
    mask = (pk_source > 0) & (k_source > 0)
    if np.sum(mask) < 2:
        return np.ones_like(k_target) * np.nan
    log_interp = np.interp(
        np.log10(k_target),
        np.log10(k_source[mask]),
        np.log10(pk_source[mask]),
        left=np.nan, right=np.nan,
    )
    return 10.0**log_interp


def find_sim_pk_files(basedir, model_name):
    """Find all pk_*.dat files for a model in output directories."""
    pk_files = {}
    for inum, nchar in enumerate(OUTPUT_NUMS):
        candidates = [
            os.path.join(basedir, model_name, f'output_{nchar}', f'pk_{nchar}.dat'),
            os.path.join(basedir, f'output_{nchar}', f'pk_{nchar}.dat'),
        ]
        for cand in candidates:
            if os.path.exists(cand):
                pk_files[inum] = cand
                break
    return pk_files


def plot_ratio_by_redshift(all_data, basedir):
    """Figure 1: P_sim/P_CAMB per redshift, 4 models in each panel."""
    fig, axes = plt.subplots(2, 3, figsize=(16, 10))
    axes = axes.flatten()

    for iz, z in enumerate(Z_OUTPUTS):
        ax = axes[iz]
        ax.set_title(f'z = {z:.1f}', fontsize=14)
        ax.axhline(1.0, color='gray', ls='--', lw=0.8)
        ax.axhspan(0.95, 1.05, color='gray', alpha=0.1)

        for model in MODELS:
            if model not in all_data or iz not in all_data[model]:
                continue
            data = all_data[model][iz]
            k = data['k_sim']
            ratio = data['ratio']
            mask = np.isfinite(ratio)
            if np.sum(mask) == 0:
                continue
            ax.semilogx(k[mask], ratio[mask], '-', color=MODEL_COLORS[model],
                        label=MODEL_LABELS[model], lw=1.5)

        ax.set_xlim(0.01, 5.0)
        ax.set_ylim(0.5, 2.0)
        ax.set_xlabel(r'$k$ [h/Mpc]')
        ax.set_ylabel(r'$P_{\rm sim}(k) / P_{\rm CAMB}(k)$')
        if iz == 0:
            ax.legend(fontsize=8, loc='upper left')

    plt.tight_layout()
    outfile = os.path.join(basedir, 'pk_ratio_by_redshift.png')
    fig.savefig(outfile, dpi=150)
    print(f"  Saved: {outfile}")
    plt.close(fig)


def plot_ratio_by_model(all_data, basedir):
    """Figure 2: P_sim/P_CAMB per model, 6 redshifts in each panel."""
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    axes = axes.flatten()

    z_colors = plt.cm.viridis(np.linspace(0.1, 0.9, len(Z_OUTPUTS)))

    for im, model in enumerate(MODELS):
        ax = axes[im]
        ax.set_title(MODEL_LABELS[model], fontsize=13)
        ax.axhline(1.0, color='gray', ls='--', lw=0.8)
        ax.axhspan(0.95, 1.05, color='gray', alpha=0.1)

        if model not in all_data:
            continue

        for iz, z in enumerate(Z_OUTPUTS):
            if iz not in all_data[model]:
                continue
            data = all_data[model][iz]
            k = data['k_sim']
            ratio = data['ratio']
            mask = np.isfinite(ratio)
            if np.sum(mask) == 0:
                continue
            aexp = data.get('aexp', None)
            z_actual = 1.0/aexp - 1 if aexp and aexp > 0 else z
            ax.semilogx(k[mask], ratio[mask], '-', color=z_colors[iz],
                        label=f'z={z_actual:.2f}', lw=1.5)

        ax.set_xlim(0.01, 5.0)
        ax.set_ylim(0.5, 2.0)
        ax.set_xlabel(r'$k$ [h/Mpc]')
        ax.set_ylabel(r'$P_{\rm sim}(k) / P_{\rm CAMB}(k)$')
        ax.legend(fontsize=8, loc='upper left')

    plt.tight_layout()
    outfile = os.path.join(basedir, 'pk_ratio_by_model.png')
    fig.savefig(outfile, dpi=150)
    print(f"  Saved: {outfile}")
    plt.close(fig)


def plot_growth_factor(all_data, basedir):
    """Figure 3: Growth factor D^2(z) = P(k,z)/P(k,z_ref) ratio."""
    z_ref_idx = 0  # z=5 as reference

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    axes = axes.flatten()

    for im, model in enumerate(MODELS):
        ax = axes[im]
        ax.set_title(MODEL_LABELS[model], fontsize=13)
        ax.axhline(1.0, color='gray', ls='--', lw=0.8)

        if model not in all_data or z_ref_idx not in all_data[model]:
            continue

        k_ref = all_data[model][z_ref_idx]['k_sim']
        pk_ref_sim = all_data[model][z_ref_idx]['pk_sim']
        pk_ref_camb = all_data[model][z_ref_idx]['pk_camb_interp']

        z_colors = plt.cm.plasma(np.linspace(0.1, 0.9, len(Z_OUTPUTS)-1))
        for iz in range(1, len(Z_OUTPUTS)):
            if iz not in all_data[model]:
                continue
            z = Z_OUTPUTS[iz]
            data = all_data[model][iz]

            pk_ref_sim_interp = interpolate_pk(data['k_sim'], k_ref, pk_ref_sim)
            pk_ref_camb_interp = interpolate_pk(data['k_sim'], k_ref, pk_ref_camb)

            d2_sim = data['pk_sim'] / pk_ref_sim_interp
            d2_camb = data['pk_camb_interp'] / pk_ref_camb_interp
            ratio = d2_sim / d2_camb

            k = data['k_sim']
            mask = np.isfinite(ratio) & (ratio > 0)
            if np.sum(mask) == 0:
                continue
            aexp = data.get('aexp', None)
            z_actual = 1.0/aexp - 1 if aexp and aexp > 0 else z
            ax.semilogx(k[mask], ratio[mask], '-', color=z_colors[iz-1],
                        label=f'z={z_actual:.2f}', lw=1.5)

        ax.set_xlim(0.01, 5.0)
        ax.set_ylim(0.7, 2.0)
        ax.set_xlabel(r'$k$ [h/Mpc]')
        ax.set_ylabel(r'$D^2_{\rm sim}(z) / D^2_{\rm CAMB}(z)$')
        ax.legend(fontsize=8, loc='upper left')

    plt.tight_layout()
    outfile = os.path.join(basedir, 'growth_factor.png')
    fig.savefig(outfile, dpi=150)
    print(f"  Saved: {outfile}")
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description="Compare sim P(k) with CAMB")
    parser.add_argument('--basedir', type=str, default='.',
                        help='Base directory containing model outputs')
    args = parser.parse_args()
    basedir = args.basedir

    print(f"Base directory: {os.path.abspath(basedir)}")
    print(f"Models: {MODELS}")
    print(f"Redshifts (nominal): {Z_OUTPUTS}")
    print(f"CAMB on-the-fly: {HAS_CAMB}")

    # ===================================================================
    # Load all data
    # ===================================================================
    all_data = {}

    for model in MODELS:
        print(f"\n--- {model} ---")

        sim_pk_files = find_sim_pk_files(basedir, model)
        if not sim_pk_files:
            print(f"  WARNING: No simulation P(k) files found for {model}")
            continue

        all_data[model] = {}
        for iz, z_nominal in enumerate(Z_OUTPUTS):
            if iz not in sim_pk_files:
                print(f"  z={z_nominal:.1f}: no sim P(k)")
                continue

            # Read sim P(k)
            k_sim, pk_sim, pk_shot, nmodes, aexp_sim = read_pk_sim(sim_pk_files[iz])
            z_actual = 1.0 / aexp_sim - 1.0 if aexp_sim and aexp_sim > 0 else z_nominal
            print(f"  z_nom={z_nominal:.1f} z_act={z_actual:.3f} (a={aexp_sim:.4f}): {len(k_sim)} k-bins")

            # Get CAMB P(k) at EXACT output redshift
            pk_camb_interp = None
            if HAS_CAMB:
                try:
                    k_camb, pk_camb = camb_pk_at_z(model, z_actual)
                    pk_camb_interp = interpolate_pk(k_sim, k_camb, pk_camb)
                    print(f"    CAMB P(k) computed at z={z_actual:.3f}")
                except Exception as e:
                    print(f"    CAMB computation failed: {e}")

            if pk_camb_interp is None:
                # Fallback: use pre-computed files at nominal z
                camb_file = os.path.join(basedir, f'camb_pk_{model}_z{z_nominal:.1f}.dat')
                if not os.path.exists(camb_file):
                    print(f"    WARNING: {camb_file} not found")
                    continue
                k_camb, pk_camb = read_pk_camb(camb_file)
                pk_camb_interp = interpolate_pk(k_sim, k_camb, pk_camb)
                print(f"    Using pre-computed CAMB P(k) at z={z_nominal:.1f} (redshift mismatch!)")

            # Compute ratio (use shot-noise corrected P(k))
            ratio = np.where(
                (pk_camb_interp > 0) & np.isfinite(pk_camb_interp),
                pk_shot / pk_camb_interp,
                np.nan,
            )

            all_data[model][iz] = {
                'k_sim': k_sim,
                'pk_sim': pk_sim,
                'pk_shot': pk_shot,
                'pk_camb_interp': pk_camb_interp,
                'ratio': ratio,
                'aexp': aexp_sim,
                'z_actual': z_actual,
            }

            # Print linear-scale ratio
            lin_mask = (k_sim > 0.01) & (k_sim < 0.1) & np.isfinite(ratio)
            if np.sum(lin_mask) > 0:
                mean_ratio = np.mean(ratio[lin_mask])
                print(f"    Linear ratio (0.01<k<0.1): {mean_ratio:.4f}")

    if not all_data:
        print("\nNo data loaded. Check file paths.")
        return

    # ===================================================================
    # Generate plots
    # ===================================================================
    print("\n=== Generating plots ===")
    plot_ratio_by_redshift(all_data, basedir)
    plot_ratio_by_model(all_data, basedir)
    plot_growth_factor(all_data, basedir)

    # ===================================================================
    # Summary statistics
    # ===================================================================
    print("\n=== Summary ===")
    print(f"{'Model':<15} {'z_actual':<10} {'<P_sim/P_CAMB>':<15} {'(k<0.1 h/Mpc)'}")
    print("-" * 55)
    for model in MODELS:
        if model not in all_data:
            continue
        for iz, z in enumerate(Z_OUTPUTS):
            if iz not in all_data[model]:
                continue
            data = all_data[model][iz]
            z_act = data.get('z_actual', z)
            lin_mask = (data['k_sim'] > 0.01) & (data['k_sim'] < 0.1) & \
                       np.isfinite(data['ratio'])
            if np.sum(lin_mask) > 0:
                mean_r = np.mean(data['ratio'][lin_mask])
                std_r = np.std(data['ratio'][lin_mask])
                print(f"{model:<15} {z_act:<10.3f} {mean_r:.4f} +/- {std_r:.4f}")

    print("\nDone.")


if __name__ == "__main__":
    main()
