#!/usr/bin/env python3
"""
DE Perturbation Linear P(k) Validation Test Setup

Generates CAMB transfer functions, linear P(k), DE tables, MUSIC IC configs,
RAMSES namelists, and SLURM scripts for 4 cosmological models:
  1. LCDM       (w0=-1.0, wa=0.0, de_perturb=F)
  2. wCDM       (w0=-0.8, wa=0.0, de_perturb=T, cs2=0)
  3. CPL-thaw   (w0=-0.9, wa=0.2, de_perturb=T, cs2=0)
  4. CPL-phantom (w0=-1.1, wa=-0.2, de_perturb=T, cs2=0)

Usage:
    cd test_ksection/run_de_pk_test/
    python3 ../../misc/de_linear_test_setup.py
"""

import os
import sys
import numpy as np

try:
    import camb
    from camb.dark_energy import DarkEnergyFluid, DarkEnergyPPF
except ImportError:
    print("ERROR: CAMB not installed. Install with: pip install camb")
    sys.exit(1)

# ===================================================================
# Common cosmological parameters (Planck 2018)
# ===================================================================
OMEGA_M = 0.3111
OMEGA_B = 0.049
H0 = 67.66
N_S = 0.9665
A_S = 2.1e-9

# Box and grid
BOXLEN = 500.0  # Mpc/h
LEVELMIN = 8    # 256^3
Z_START = 17
SEED = 12345

# Output redshifts
Z_OUTPUTS = [5.0, 2.0, 1.0, 0.5, 0.2, 0.0]
A_OUTPUTS = [1.0/(1.0+z) for z in Z_OUTPUTS]

# Paths
MUSIC_BIN = "/gpfs/kjhan/Hydro/MUSIC/MUSIC2/build/MUSIC"
RAMSES_BIN = "../../bin/ramses_final3d"
DE_TABLE_SCRIPT = "../../misc/generate_de_table.py"

# ===================================================================
# Model definitions
# ===================================================================
MODELS = {
    'lcdm': {
        'w0': -1.0, 'wa': 0.0, 'cs2_de': 0.0,
        'de_perturb': False,
        'desc': 'LCDM baseline',
    },
    'wcdm': {
        'w0': -0.8, 'wa': 0.0, 'cs2_de': 0.0,
        'de_perturb': True,
        'desc': 'Constant w=-0.8, clustering DE',
    },
    'cpl_thaw': {
        'w0': -0.9, 'wa': 0.2, 'cs2_de': 0.0,
        'de_perturb': True,
        'desc': 'CPL thawing (w0=-0.9, wa=0.2)',
    },
    'cpl_phantom': {
        'w0': -1.1, 'wa': -0.2, 'cs2_de': 0.0,
        'de_perturb': True,
        'desc': 'CPL phantom crossing (w0=-1.1, wa=-0.2)',
    },
}


def setup_camb_params(w0, wa, cs2_de, redshifts, kmax=50.0):
    """Create CAMB parameters for given DE model."""
    h = H0 / 100.0
    omega_cdm = OMEGA_M - OMEGA_B

    pars = camb.CAMBparams()
    pars.set_cosmology(
        H0=H0,
        ombh2=OMEGA_B * h**2,
        omch2=omega_cdm * h**2,
        mnu=0.0,
    )
    pars.InitPower.set_params(ns=N_S, As=A_S)
    # Use PPF if w ever goes below -1 (DarkEnergyFluid rejects w<-1)
    w_today = w0
    w_early = w0 + wa
    needs_ppf = (w_today < -1) or (w_early < -1)
    if needs_ppf:
        pars.DarkEnergy = DarkEnergyPPF(w=w0, wa=wa)
    else:
        pars.DarkEnergy = DarkEnergyFluid(w=w0, wa=wa, cs2=cs2_de)
    pars.set_matter_power(
        redshifts=sorted(redshifts, reverse=True),
        kmax=kmax,
    )
    pars.WantTransfer = True
    return pars


def generate_camb_transfer(model_name, model, outdir='.'):
    """Generate CAMB transfer function file for MUSIC IC generation."""
    print(f"\n--- {model_name}: CAMB transfer at z={Z_START} ---")

    pars = setup_camb_params(model['w0'], model['wa'], model['cs2_de'],
                              [Z_START, 0.0], kmax=50.0)
    results = camb.get_results(pars)

    # Get transfer data (z=17 is the first redshift in sorted list)
    trans = results.get_matter_transfer_data()
    # transfer_data shape: (nvar, nk, nz) with nvar=13
    # CAMB indices: 0=k/h, 1=cdm, 2=b, 3=photon, 4=nu_massless, 5=nu_massive,
    #               6=total, 7=no_nu, 8=total_de, 9=Weyl, 10=v_cdm, 11=v_b, 12=v_b-v_cdm
    # MUSIC2 expects all 13 columns (reads up to v_b column 12)
    kh = trans.transfer_data[0, :, 0]

    transfer_file = os.path.join(outdir, f'transfer_{model_name}_z{Z_START}.dat')
    nvar = trans.transfer_data.shape[0]
    with open(transfer_file, 'w') as f:
        nk = len(kh)
        iz = 0  # z=17 is first (CAMB descending sort)
        for ik in range(nk):
            vals = [trans.transfer_data[iv, ik, iz] for iv in range(nvar)]
            f.write(' '.join(f'{v:15.8e}' for v in vals) + '\n')

    print(f"  Transfer file: {transfer_file} ({nk} k-points)")

    # Also get sigma_8 for MUSIC config
    sigma8 = results.get_sigma8_0()
    print(f"  sigma_8(z=0) = {sigma8:.6f}")

    return transfer_file, sigma8


def generate_camb_pk(model_name, model, outdir='.'):
    """Generate CAMB linear P(k) at 6 output redshifts for comparison."""
    print(f"\n--- {model_name}: CAMB linear P(k) at z={Z_OUTPUTS} ---")

    pars = setup_camb_params(model['w0'], model['wa'], model['cs2_de'],
                              Z_OUTPUTS, kmax=10.0)
    pars.NonLinear = camb.model.NonLinear_none  # linear only
    results = camb.get_results(pars)

    for z in Z_OUTPUTS:
        # Get linear matter power spectrum
        kh, z_out, pk = results.get_linear_matter_power_spectrum()
        # Find the index for this redshift
        iz = np.argmin(np.abs(np.array(z_out) - z))

        pk_file = os.path.join(outdir, f'camb_pk_{model_name}_z{z:.1f}.dat')
        with open(pk_file, 'w') as f:
            f.write(f"# CAMB linear P(k) at z={z_out[iz]:.4f} (a={1/(1+z_out[iz]):.6f})\n")
            f.write(f"# w0={model['w0']} wa={model['wa']} cs2_de={model['cs2_de']}\n")
            f.write(f"# k [h/Mpc]    P(k) [(Mpc/h)^3]\n")
            for ik in range(len(kh)):
                f.write(f"{kh[ik]:15.8e} {pk[iz, ik]:15.8e}\n")
        print(f"  z={z:.1f}: {pk_file} ({len(kh)} k-points)")


def generate_de_table(model_name, model, outdir='.'):
    """Generate DE transfer function table using Weyl potential method.

    For phantom models (w < -1), DE perturbation is disabled since
    DarkEnergyFluid/PPF doesn't produce significant sub-Hubble clustering.
    """
    if not model['de_perturb']:
        return None

    w0, wa = model['w0'], model['wa']
    w_min = min(w0, w0 + wa)
    if w_min < -1.0:
        print(f"\n--- {model_name}: Phantom model (w_min={w_min:.2f}), "
              f"DE perturbation disabled ---")
        # Disable de_perturb for phantom — PPF doesn't cluster on sub-Hubble scales
        model['de_perturb'] = False
        return None

    de_file = os.path.join(outdir, f'de_table_{model_name}.dat')
    print(f"\n--- {model_name}: DE table → {de_file} ---")

    cs2_de = model['cs2_de']
    h = H0 / 100.0

    sys.path.insert(0, os.path.dirname(os.path.abspath(DE_TABLE_SCRIPT)))
    from generate_de_table import generate_de_table as gen_de
    gen_de(
        output_file=de_file,
        w0=w0, wa=wa, cs2_de=cs2_de,
        omega_m=OMEGA_M, omega_b=OMEGA_B, h0=h,
        nk=200, na=50,
    )

    return de_file

def write_music_config(model_name, model, sigma8, transfer_file, outdir='.'):
    """Write MUSIC IC configuration file."""
    music_file = os.path.join(outdir, f'music_{model_name}.conf')
    print(f"\n--- {model_name}: MUSIC config → {music_file} ---")

    h = H0 / 100.0
    omega_l = 1.0 - OMEGA_M  # flat universe

    config = f"""[setup]
boxlength       = {BOXLEN}
zstart          = {Z_START}
levelmin        = {LEVELMIN}
levelmax        = {LEVELMIN}
baryons         = no
use_2LPT        = yes
use_LLA         = no
padding         = 8

[cosmology]
Omega_m         = {OMEGA_M}
Omega_L         = {omega_l}
Omega_b         = {OMEGA_B}
H0              = {H0}
sigma_8         = {sigma8:.6f}
n_s             = {N_S}
w_0             = {model['w0']}
w_a             = {model['wa']}
transfer        = camb_file
transfer_file   = {os.path.basename(transfer_file)}

[random]
seed[{LEVELMIN}]         = {SEED}

[output]
format          = grafic2
filename        = ics_{model_name}

[poisson]
fft_fine        = yes
accuracy        = 1e-5
pre_smooth      = 3
post_smooth     = 3
smoother        = gs
laplace_order   = 6
grad_order      = 6
"""
    with open(music_file, 'w') as f:
        f.write(config)


def write_ramses_namelist(model_name, model, de_table_file, outdir='.'):
    """Write RAMSES namelist file."""
    nml_file = os.path.join(outdir, f'cosmo_{model_name}.nml')
    print(f"\n--- {model_name}: RAMSES namelist → {nml_file} ---")

    de_perturb_str = '.true.' if model['de_perturb'] else '.false.'
    aout_str = ','.join(f'{a:.6f}' for a in A_OUTPUTS)

    nml = f"""&RUN_PARAMS
cosmo=.true.
pic=.true.
poisson=.true.
hydro=.false.
sink=.false.
nrestart=0
nsubcycle=1,1,1
ordering='ksection'
memory_balance=.true.
use_fftw=.true.
de_perturb={de_perturb_str}
dump_pk=.true.
/

&OUTPUT_PARAMS
noutput=6
aout={aout_str}
foutput=100000
/

&COSMO_PARAMS
omega_m={OMEGA_M}
omega_l={1.0 - OMEGA_M}
omega_b=0.0  ! gravity-only: particles represent total matter
h0={H0}
/
"""

    # Always include CPL_PARAMS for non-LCDM models (background expansion)
    # de_perturb only controls Poisson correction, not background cosmology
    is_non_lcdm = (model['w0'] != -1.0) or (model['wa'] != 0.0)
    if is_non_lcdm:
        nml += f"""
&CPL_PARAMS
w0={model['w0']}
wa={model['wa']}
"""
        if model['de_perturb'] and de_table_file:
            nml += f"""cs2_de={model['cs2_de']}
de_table='{os.path.basename(de_table_file)}'
"""
        nml += "/\n"

    nml += f"""
&INIT_PARAMS
filetype='grafic'
initfile(1)='ics_{model_name}/level_008'
/

&AMR_PARAMS
levelmin={LEVELMIN}
levelmax={LEVELMIN}
nexpand=1
ngridtot=30000000
nparttot=30000000
/

&POISSON_PARAMS
epsilon=1d-4
/
"""

    with open(nml_file, 'w') as f:
        f.write(nml)


def write_slurm_script(model_name, model, de_table_file, outdir='.'):
    """Write SLURM job submission script.

    Each model runs in its own subdirectory to avoid output conflicts.
    """
    # Create model subdirectory
    model_dir = os.path.join(outdir, model_name)
    os.makedirs(model_dir, exist_ok=True)

    script_file = os.path.join(outdir, f'run_{model_name}.sh')
    print(f"\n--- {model_name}: SLURM script → {script_file} ---")

    # Symlinks needed in model subdir (relative paths)
    symlinks = [f'../cosmo_{model_name}.nml', f'../ics_{model_name}']
    if de_table_file:
        symlinks.append(f'../{os.path.basename(de_table_file)}')

    symlink_cmds = '\n'.join(
        f'ln -sf {s} .' for s in symlinks
    )

    script = f"""#!/bin/bash
#SBATCH -J de_pk_{model_name}
#SBATCH -p a10
#SBATCH -N 1
#SBATCH --ntasks-per-node=12
#SBATCH --cpus-per-task=4
#SBATCH -t 04:00:00
#SBATCH -o de_pk_{model_name}_%j.log
#SBATCH -e de_pk_{model_name}_%j.err

export OMP_NUM_THREADS=4
export OMP_STACKSIZE=256M
export LD_LIBRARY_PATH=/home/kjhan/local/hdf5/lib:$LD_LIBRARY_PATH

echo "=== DE P(k) test: {model_name} ==="
echo "Start: $(date)"
echo "Nodes: $SLURM_NNODES, Tasks: $SLURM_NTASKS, OMP: $OMP_NUM_THREADS"

# Run in model subdirectory so outputs don't collide
cd {model_name}
{symlink_cmds}

mpirun -np 12 ../../../bin/ramses_final3d cosmo_{model_name}.nml

echo "End: $(date)"
"""

    with open(script_file, 'w') as f:
        f.write(script)
    os.chmod(script_file, 0o755)


def write_run_all_script(outdir='.'):
    """Write master script to generate ICs and run all models."""
    script_file = os.path.join(outdir, 'run_all.sh')
    print(f"\n--- Master script → {script_file} ---")

    models_str = ' '.join(MODELS.keys())
    script = f"""#!/bin/bash
# DE P(k) Validation Test - Master Run Script
# Usage: ./run_all.sh [ic|run|compare|all]

set -e

MUSIC="{MUSIC_BIN}"
MODELS="{models_str}"

case "${{1:-all}}" in
  ic)
    echo "=== Generating ICs with MUSIC ==="
    for model in $MODELS; do
      echo "--- MUSIC: $model ---"
      $MUSIC music_${{model}}.conf
    done
    ;;

  run)
    echo "=== Submitting RAMSES jobs ==="
    for model in $MODELS; do
      echo "--- Submitting: $model ---"
      sbatch run_${{model}}.sh
    done
    ;;

  compare)
    echo "=== Comparing P(k) with CAMB ==="
    python3 ../../misc/compare_pk_camb.py --basedir .
    ;;

  all)
    echo "=== Step 1: Generate ICs ==="
    $0 ic
    echo ""
    echo "=== Step 2: Submit jobs ==="
    $0 run
    echo ""
    echo "After jobs finish, run: $0 compare"
    ;;

  *)
    echo "Usage: $0 [ic|run|compare|all]"
    exit 1
    ;;
esac
"""

    with open(script_file, 'w') as f:
        f.write(script)
    os.chmod(script_file, 0o755)


def main():
    outdir = '.'
    print(f"CAMB version: {camb.__version__}")
    print(f"Output directory: {os.path.abspath(outdir)}")
    print(f"Models: {list(MODELS.keys())}")

    for model_name, model in MODELS.items():
        print(f"\n{'='*60}")
        print(f"Model: {model_name} — {model['desc']}")
        print(f"  w0={model['w0']}, wa={model['wa']}, cs2_de={model['cs2_de']}")
        print(f"  de_perturb={model['de_perturb']}")
        print(f"{'='*60}")

        # 1a. CAMB transfer function for MUSIC
        transfer_file, sigma8 = generate_camb_transfer(model_name, model, outdir)

        # 1b. CAMB linear P(k) at output redshifts
        generate_camb_pk(model_name, model, outdir)

        # 1c. DE table (only for non-LCDM)
        de_table_file = generate_de_table(model_name, model, outdir)

        # 1d. MUSIC config
        write_music_config(model_name, model, sigma8, transfer_file, outdir)

        # 1e. RAMSES namelist
        write_ramses_namelist(model_name, model, de_table_file, outdir)

        # 1f. SLURM script
        write_slurm_script(model_name, model, de_table_file, outdir)

    # Master run script
    write_run_all_script(outdir)

    print(f"\n{'='*60}")
    print("Setup complete! Next steps:")
    print(f"  1. cd test_ksection/run_de_pk_test/")
    print(f"  2. ./run_all.sh ic    # Generate ICs with MUSIC")
    print(f"  3. ./run_all.sh run   # Submit RAMSES jobs")
    print(f"  4. ./run_all.sh compare  # Compare P(k) after jobs finish")
    print(f"{'='*60}")


if __name__ == "__main__":
    main()
