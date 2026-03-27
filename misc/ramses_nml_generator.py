#!/usr/bin/env python3
"""
RAMSES Namelist Generator
=========================
Interactive and CLI tool for generating RAMSES namelist (.nml) files.
Supports presets, parameter editing, validation, import/export.

Usage:
  python3 ramses_nml_generator.py                                    # Interactive
  python3 ramses_nml_generator.py --preset hr5_production -o out.nml # Non-interactive
  python3 ramses_nml_generator.py --import old.nml -o new.nml        # Import+edit
  python3 ramses_nml_generator.py --list-presets                     # Show presets

Requires: Python 3.6+ (no external dependencies)
"""
from __future__ import print_function
import argparse
import copy
import os
import re
import sys
from collections import OrderedDict

# ============================================================================
# Parameter Definition
# ============================================================================
class ParamDef:
    """Definition of a single RAMSES namelist parameter."""
    __slots__ = ('name', 'ftype', 'default', 'group', 'section', 'desc',
                 'visible_when', 'choices', 'validate')
    def __init__(self, name, ftype, default, group, section, desc,
                 visible_when='', choices=None, validate=''):
        self.name = name
        self.ftype = ftype          # 'int','real','bool','str','int_arr','real_arr','str_arr'
        self.default = default
        self.group = group          # Fortran namelist group
        self.section = section      # UI section for interactive display
        self.desc = desc
        self.visible_when = visible_when
        self.choices = choices
        self.validate = validate

# --- Section names ---
S_SIMTYPE  = '1. Simulation Type'
S_COSMO    = '2. Cosmology'
S_AMR      = '3. AMR Grid'
S_IC       = '4. Initial Conditions'
S_OUTPUT   = '5. Output'
S_TIME     = '6. Time Stepping'
S_DOMAIN   = '7. Domain Decomposition'
S_HYDRO    = '8. Hydro Solver'
S_REFINE   = '9. Refinement'
S_FPR      = '10. FPR / Holdback'
S_COOL     = '11. Cooling & Star Formation'
S_FEED     = '12. Feedback (SN/AGN)'
S_GPU      = '13. GPU & FFTW'
S_MODGRAV  = '14. Modified Gravity'
S_DE       = '15. Dark Energy'
S_NU       = '16. Neutrino'
S_SIDM     = '17. SIDM'
S_SGS      = '18. SGS Turbulence'
S_POISSON  = '19. Poisson Solver'
S_LIGHTCONE = '20. Lightcone'

# --- Full parameter database ---
PARAMS = [
    # ====== RUN_PARAMS ======
    ParamDef('cosmo',   'bool', False, 'RUN_PARAMS', S_SIMTYPE, 'Cosmological simulation'),
    ParamDef('pic',     'bool', False, 'RUN_PARAMS', S_SIMTYPE, 'Particle-in-cell (N-body)'),
    ParamDef('poisson', 'bool', False, 'RUN_PARAMS', S_SIMTYPE, 'Self-gravity (Poisson solver)'),
    ParamDef('hydro',   'bool', False, 'RUN_PARAMS', S_SIMTYPE, 'Hydrodynamics'),
    ParamDef('sink',    'bool', False, 'RUN_PARAMS', S_SIMTYPE, 'Sink particles (BH/stars)'),
    ParamDef('rt',      'bool', False, 'RUN_PARAMS', S_SIMTYPE, 'Radiative transfer'),
    ParamDef('lightcone','bool',False, 'RUN_PARAMS', S_SIMTYPE, 'Lightcone output'),
    ParamDef('clumpfind','bool',False, 'RUN_PARAMS', S_SIMTYPE, 'Clump finder'),
    ParamDef('verbose', 'bool', False, 'RUN_PARAMS', S_SIMTYPE, 'Verbose output'),
    ParamDef('debug',   'bool', False, 'RUN_PARAMS', S_SIMTYPE, 'Debug mode'),
    ParamDef('static',  'bool', False, 'RUN_PARAMS', S_SIMTYPE, 'Static grid (no refinement)'),

    # Time stepping
    ParamDef('nrestart',   'int',  0,      'RUN_PARAMS', S_TIME, 'Restart from output number (0=fresh)'),
    ParamDef('ncontrol',   'int',  1,      'RUN_PARAMS', S_TIME, 'Steps between control output'),
    ParamDef('nstepmax',   'int',  1000000,'RUN_PARAMS', S_TIME, 'Maximum number of coarse steps'),
    ParamDef('nsubcycle',  'int_arr','1,1,2','RUN_PARAMS',S_TIME, 'Subcycling per level (csv)'),

    # Domain decomposition
    ParamDef('nremap',     'int',  0,      'RUN_PARAMS', S_DOMAIN, 'Steps between load rebalancing'),
    ParamDef('ordering',   'str',  'hilbert','RUN_PARAMS',S_DOMAIN,'Domain ordering',
             choices=['hilbert','morton','ksection']),
    ParamDef('memory_balance','bool',False,'RUN_PARAMS', S_DOMAIN, 'Memory-based load balancing'),
    ParamDef('mem_weight_grid','int',270,  'RUN_PARAMS', S_DOMAIN, 'Memory weight per grid cell (bytes)',
             visible_when='memory_balance==True'),
    ParamDef('mem_weight_part','int',12,   'RUN_PARAMS', S_DOMAIN, 'Memory weight per particle (bytes)',
             visible_when='memory_balance==True'),
    ParamDef('mem_weight_sink','int',500,  'RUN_PARAMS', S_DOMAIN, 'Memory weight per sink (bytes)',
             visible_when='memory_balance==True'),
    ParamDef('exchange_method','str','auto','RUN_PARAMS', S_DOMAIN,'MPI exchange method',
             choices=['auto','alltoall','ksection','p2p']),
    ParamDef('overload',   'real', 1.0,    'RUN_PARAMS', S_DOMAIN, 'Load imbalance tolerance factor'),
    ParamDef('cost_weighting','bool',False,'RUN_PARAMS', S_DOMAIN, 'Cost-based weighting'),
    ParamDef('bisec_tol',  'real', 0.05,   'RUN_PARAMS', S_DOMAIN, 'Bisection tolerance'),

    # GPU & FFTW
    ParamDef('gpu_hydro',  'bool', False,  'RUN_PARAMS', S_GPU, 'GPU-accelerated hydro'),
    ParamDef('gpu_poisson','bool', False,  'RUN_PARAMS', S_GPU, 'GPU-accelerated Poisson MG'),
    ParamDef('gpu_fft',    'bool', False,  'RUN_PARAMS', S_GPU, 'GPU FFT for direct Poisson'),
    ParamDef('gpu_sink',   'bool', False,  'RUN_PARAMS', S_GPU, 'GPU sink particle'),
    ParamDef('gpu_auto_tune','bool',True,  'RUN_PARAMS', S_GPU, 'Auto-tune CPU vs GPU (disable for benchmarks)'),
    ParamDef('use_fftw',   'bool', False,  'RUN_PARAMS', S_GPU, 'FFTW3 CPU direct Poisson solver'),

    # Misc run_params
    ParamDef('jobcontrolfile','str','',    'RUN_PARAMS', S_TIME, 'Job control file path'),
    ParamDef('dump_pk',    'bool', False,  'RUN_PARAMS', S_OUTPUT, 'Dump power spectrum'),
    ParamDef('sinkprops',  'bool', False,  'RUN_PARAMS', S_SIMTYPE, 'Output sink properties'),

    # Physics toggles in RUN_PARAMS
    ParamDef('use_neutrino','bool',False,  'RUN_PARAMS', S_NU,  'Massive neutrino linear response'),
    ParamDef('sidm',       'bool', False,  'RUN_PARAMS', S_SIDM,'Self-interacting dark matter'),
    ParamDef('de_perturb', 'bool', False,  'RUN_PARAMS', S_DE,  'Dark energy perturbations'),
    ParamDef('use_mond',   'bool', False,  'RUN_PARAMS', S_MODGRAV,'MOND (QUMOND/AQUAL)'),
    ParamDef('use_fR',     'bool', False,  'RUN_PARAMS', S_MODGRAV,'f(R) Hu-Sawicki gravity'),
    ParamDef('use_nDGP',   'bool', False,  'RUN_PARAMS', S_MODGRAV,'nDGP braneworld gravity'),
    ParamDef('use_symmetron','bool',False, 'RUN_PARAMS', S_MODGRAV,'Symmetron scalar field'),
    ParamDef('use_dilaton','bool', False,  'RUN_PARAMS', S_MODGRAV,'Dilaton scalar field'),
    ParamDef('use_galileon','bool',False,  'RUN_PARAMS', S_MODGRAV,'Galileon scalar field'),
    ParamDef('use_coupled_de','bool',False,'RUN_PARAMS', S_DE,  'Coupled dark energy'),
    ParamDef('use_ede',    'bool', False,  'RUN_PARAMS', S_DE,  'Early dark energy'),
    ParamDef('use_sgs',    'bool', False,  'RUN_PARAMS', S_SGS, 'SGS turbulence model'),

    # ====== COSMO_PARAMS ======
    ParamDef('omega_b',  'real', 0.045,  'COSMO_PARAMS', S_COSMO, 'Baryon density parameter'),
    ParamDef('omega_m',  'real', 0.3,    'COSMO_PARAMS', S_COSMO, 'Matter density parameter'),
    ParamDef('omega_l',  'real', 0.7,    'COSMO_PARAMS', S_COSMO, 'Dark energy density parameter'),
    ParamDef('h0',       'real', 0.7,    'COSMO_PARAMS', S_COSMO, 'Hubble parameter (H0/100)'),

    # ====== OUTPUT_PARAMS ======
    ParamDef('noutput',   'int',  1,     'OUTPUT_PARAMS', S_OUTPUT, 'Number of output snapshots'),
    ParamDef('foutput',   'int',  10000, 'OUTPUT_PARAMS', S_OUTPUT, 'Steps between forced outputs'),
    ParamDef('fbackup',   'int',  10000, 'OUTPUT_PARAMS', S_OUTPUT, 'Steps between backup outputs'),
    ParamDef('aout',      'real_arr','',  'OUTPUT_PARAMS', S_OUTPUT, 'Output scale factors (csv)',
             visible_when='cosmo==True'),
    ParamDef('tout',      'real_arr','',  'OUTPUT_PARAMS', S_OUTPUT, 'Output times (csv)',
             visible_when='cosmo==False'),
    ParamDef('outformat', 'str',  '',     'OUTPUT_PARAMS', S_OUTPUT, 'Output format',
             choices=['original','hdf5']),
    ParamDef('output_mode','str', '',     'OUTPUT_PARAMS', S_OUTPUT, 'Output mode'),
    ParamDef('walltime_hrs','real',0.0,   'OUTPUT_PARAMS', S_OUTPUT, 'Walltime limit (hours, 0=off)'),
    ParamDef('minutes_dump','real',0.0,   'OUTPUT_PARAMS', S_OUTPUT, 'Minutes before walltime to dump'),

    # ====== INIT_PARAMS ======
    ParamDef('filetype',  'str',  'grafic','INIT_PARAMS', S_IC, 'IC file type',
             choices=['grafic','ascii','gadget','one','two']),
    ParamDef('initfile(1)','str', '',      'INIT_PARAMS', S_IC, 'IC directory path'),
    ParamDef('multiple',  'bool', False,   'INIT_PARAMS', S_IC, 'Multiple IC files'),
    ParamDef('aexp_ini',  'real', 0.0,     'INIT_PARAMS', S_IC, 'Initial expansion factor',
             visible_when='cosmo==True'),

    # ====== AMR_PARAMS ======
    ParamDef('levelmin',  'int',  7,       'AMR_PARAMS', S_AMR, 'Minimum refinement level'),
    ParamDef('levelmax',  'int',  7,       'AMR_PARAMS', S_AMR, 'Maximum refinement level'),
    ParamDef('nexpand',   'int',  1,       'AMR_PARAMS', S_AMR, 'Buffer cells for refinement'),
    ParamDef('ngridtot',  'int',  0,       'AMR_PARAMS', S_AMR, 'Total grid allocation (0=auto)'),
    ParamDef('nparttot',  'int',  0,       'AMR_PARAMS', S_AMR, 'Total particle allocation (0=auto)'),
    ParamDef('ngridmax',  'int',  0,       'AMR_PARAMS', S_AMR, 'Max grids per CPU (0=auto)'),
    ParamDef('npartmax',  'int',  0,       'AMR_PARAMS', S_AMR, 'Max particles per CPU (0=auto)'),
    ParamDef('boxlen',    'real', 1.0,     'AMR_PARAMS', S_AMR, 'Box length (code units)'),
    ParamDef('nsinkmax',  'int',  2000,    'AMR_PARAMS', S_AMR, 'Max number of sink particles',
             visible_when='sink==True'),

    # ====== LIGHTCONE_PARAMS ======
    ParamDef('zmax_cone', 'real', 0.0,     'LIGHTCONE_PARAMS', S_LIGHTCONE, 'Max redshift for lightcone',
             visible_when='lightcone==True'),

    # ====== REFINE_PARAMS ======
    ParamDef('m_refine',  'real_arr','',   'REFINE_PARAMS', S_REFINE, 'Mass refinement threshold per level'),
    ParamDef('ivar_refine','int', -1,      'REFINE_PARAMS', S_REFINE, 'Variable for refinement (0=mass)'),
    ParamDef('interpol_var','int', 0,      'REFINE_PARAMS', S_REFINE, 'Interpolation variable (0=cons,1=prim)'),
    ParamDef('interpol_type','int',1,      'REFINE_PARAMS', S_REFINE, 'Interpolation type (0=minmod,1=linear)'),
    ParamDef('jeans_refine','real_arr','',  'REFINE_PARAMS', S_REFINE, 'Jeans refinement criteria per level'),
    ParamDef('mass_sph',  'real', 0.0,     'REFINE_PARAMS', S_REFINE, 'SPH kernel mass for refinement'),
    ParamDef('err_grad_d','real', -1.0,    'REFINE_PARAMS', S_REFINE, 'Density gradient error threshold'),
    ParamDef('err_grad_p','real', -1.0,    'REFINE_PARAMS', S_REFINE, 'Pressure gradient error threshold'),
    ParamDef('err_grad_u','real', -1.0,    'REFINE_PARAMS', S_REFINE, 'Velocity gradient error threshold'),
    ParamDef('sink_refine','bool', False,  'REFINE_PARAMS', S_REFINE, 'Refine around sink particles',
             visible_when='sink==True'),
    ParamDef('var_cut_refine','real',-1.0, 'REFINE_PARAMS', S_REFINE, 'Variable cut for refinement'),

    # FPR / Holdback
    ParamDef('q_refine_holdback','bool',False,'REFINE_PARAMS',S_FPR, 'Graduated levelmax (holdback)'),
    ParamDef('dr_proper', 'real', 0.0,     'REFINE_PARAMS', S_FPR, 'Proper resolution target [kpc]'),
    ParamDef('jeans_ncells','int', 0,       'PHYSICS_PARAMS',S_FPR, 'Jeans floor: N cells per Jeans length',
             visible_when='dr_proper>0'),
    ParamDef('d_jeans_thre','real',0.0,    'REFINE_PARAMS', S_FPR, 'Jeans density threshold'),

    # ====== HYDRO_PARAMS ======
    ParamDef('gamma',     'real', 1.4,     'HYDRO_PARAMS', S_HYDRO, 'Adiabatic index',
             visible_when='hydro==True'),
    ParamDef('courant_factor','real',0.8,  'HYDRO_PARAMS', S_HYDRO, 'Courant number (CFL)',
             visible_when='hydro==True', validate='0<x<=1'),
    ParamDef('scheme',    'str',  'muscl', 'HYDRO_PARAMS', S_HYDRO, 'Hydro scheme',
             visible_when='hydro==True', choices=['muscl','plmde','induction']),
    ParamDef('slope_type','int',  1,       'HYDRO_PARAMS', S_HYDRO, 'Slope limiter (1=minmod,2=MC,3=superbee)',
             visible_when='hydro==True'),
    ParamDef('riemann',   'str',  'hllc',  'HYDRO_PARAMS', S_HYDRO, 'Riemann solver',
             visible_when='hydro==True', choices=['exact','hll','hllc','hlld','roe','llf','acoustic']),
    ParamDef('smallr',    'real', 1e-10,   'HYDRO_PARAMS', S_HYDRO, 'Minimum density',
             visible_when='hydro==True'),
    ParamDef('smallc',    'real', 1e-10,   'HYDRO_PARAMS', S_HYDRO, 'Minimum sound speed',
             visible_when='hydro==True'),
    ParamDef('difmag',    'real', 0.0,     'HYDRO_PARAMS', S_HYDRO, 'Magnetic diffusivity',
             visible_when='hydro==True'),
    ParamDef('pressure_fix','bool',False,  'HYDRO_PARAMS', S_HYDRO, 'Pressure fix for low-p regions',
             visible_when='hydro==True'),
    ParamDef('beta_fix',  'real', 0.0,     'HYDRO_PARAMS', S_HYDRO, 'Beta for pressure fix',
             visible_when='hydro==True'),
    ParamDef('niter_riemann','int',10,     'HYDRO_PARAMS', S_HYDRO, 'Iterations for exact Riemann',
             visible_when='hydro==True'),

    # ====== POISSON_PARAMS ======
    ParamDef('epsilon',   'real', 1e-4,    'POISSON_PARAMS', S_POISSON, 'Poisson solver tolerance',
             visible_when='poisson==True'),
    ParamDef('gravity_type','int',0,       'POISSON_PARAMS', S_POISSON, 'Gravity type (0=self,1-3=ext)',
             visible_when='poisson==True'),
    ParamDef('cg_levelmin','int', 0,       'POISSON_PARAMS', S_POISSON, 'CG solver min level',
             visible_when='poisson==True'),

    # ====== PHYSICS_PARAMS ======
    # Cooling & SF
    ParamDef('cooling',   'bool', False,   'PHYSICS_PARAMS', S_COOL, 'Radiative cooling'),
    ParamDef('haardt_madau','bool',False,   'PHYSICS_PARAMS', S_COOL, 'Haardt-Madau UV background'),
    ParamDef('metal',     'bool', False,   'PHYSICS_PARAMS', S_COOL, 'Metal cooling'),
    ParamDef('isothermal','bool', False,   'PHYSICS_PARAMS', S_COOL, 'Isothermal EOS'),
    ParamDef('cooling_method','str','original','PHYSICS_PARAMS',S_COOL,'Cooling method',
             choices=['original','exact']),
    ParamDef('grackle_table','str','',     'PHYSICS_PARAMS', S_COOL, 'Grackle/Eunha cooling table file',
             visible_when="cooling_method=='exact'"),
    ParamDef('z_reion',   'real', 8.5,     'PHYSICS_PARAMS', S_COOL, 'Reionization redshift',
             visible_when='cooling==True'),
    ParamDef('z_ave',     'real', 0.0,     'PHYSICS_PARAMS', S_COOL, 'Average metallicity'),
    ParamDef('J21',       'real', 0.0,     'PHYSICS_PARAMS', S_COOL, 'LW background intensity'),
    ParamDef('a_spec',    'real', 1.0,     'PHYSICS_PARAMS', S_COOL, 'UV spectral index'),
    ParamDef('self_shielding','bool',False,'PHYSICS_PARAMS', S_COOL, 'Self-shielding from UV'),
    ParamDef('T2max',     'real', 0.0,     'PHYSICS_PARAMS', S_COOL, 'Max temperature for cooling'),

    # Star formation
    ParamDef('eps_star',  'real', 0.0,     'PHYSICS_PARAMS', S_COOL, 'Star formation efficiency'),
    ParamDef('n_star',    'real', 0.1,     'PHYSICS_PARAMS', S_COOL, 'SF density threshold [H/cc]'),
    ParamDef('m_star',    'real', 0.0,     'PHYSICS_PARAMS', S_COOL, 'Star particle mass (code units, 0=auto)'),
    ParamDef('t_star',    'real', 0.0,     'PHYSICS_PARAMS', S_COOL, 'SF timescale [Gyr]'),
    ParamDef('T2_star',   'real', 0.0,     'PHYSICS_PARAMS', S_COOL, 'Polytropic floor temperature [K]'),
    ParamDef('g_star',    'real', 1.0,     'PHYSICS_PARAMS', S_COOL, 'Polytropic index for SF EOS'),
    ParamDef('del_star',  'real', 200.0,   'PHYSICS_PARAMS', S_COOL, 'SF overdensity threshold'),
    ParamDef('T2thres_SF','real', 0.0,     'PHYSICS_PARAMS', S_COOL, 'Temperature threshold for SF'),
    ParamDef('sf_virial', 'int',  0,       'PHYSICS_PARAMS', S_COOL, 'Virial SF criterion (0=off)'),
    ParamDef('sf_model',  'int',  0,       'PHYSICS_PARAMS', S_COOL, 'SF model selection'),
    ParamDef('star_maker','int',  0,       'PHYSICS_PARAMS', S_COOL, 'Star maker algorithm'),
    ParamDef('star_imf',  'str',  '',      'PHYSICS_PARAMS', S_COOL, 'IMF type'),

    # Feedback (SN)
    ParamDef('f_w',       'real', 0.0,     'PHYSICS_PARAMS', S_FEED, 'SN mass loading factor'),
    ParamDef('f_ek',      'real', 0.0,     'PHYSICS_PARAMS', S_FEED, 'SN kinetic energy fraction'),
    ParamDef('delayed_cooling','bool',False,'PHYSICS_PARAMS',S_FEED, 'Delayed cooling after SN'),
    ParamDef('rbubble',   'real', 0.0,     'PHYSICS_PARAMS', S_FEED, 'SN bubble radius [cells]'),
    ParamDef('ndebris',   'int',  0,       'PHYSICS_PARAMS', S_FEED, 'Number of SN debris particles'),
    ParamDef('mass_gmc',  'real', 0.0,     'PHYSICS_PARAMS', S_FEED, 'GMC mass [Msun]'),
    ParamDef('kappa_IR',  'real', 0.0,     'PHYSICS_PARAMS', S_FEED, 'IR opacity [cm^2/g]'),
    ParamDef('eps_sn1',   'real', 0.0,     'PHYSICS_PARAMS', S_FEED, 'SNIa energy parameter'),
    ParamDef('eps_sn2',   'real', 0.0,     'PHYSICS_PARAMS', S_FEED, 'SNII energy parameter'),
    ParamDef('tol',       'real', 0.0,     'PHYSICS_PARAMS', S_FEED, 'Tolerance for feedback'),
    ParamDef('E_SNII',    'real', 0.0,     'PHYSICS_PARAMS', S_FEED, 'SNII energy [erg]'),
    ParamDef('M_SNII',    'real', 0.0,     'PHYSICS_PARAMS', S_FEED, 'SNII progenitor mass [Msun]'),
    ParamDef('SF_kick_kms','real',0.0,     'PHYSICS_PARAMS', S_FEED, 'SF natal kick [km/s]'),
    ParamDef('SN_dT2_min','real', 0.0,     'PHYSICS_PARAMS', S_FEED, 'SN minimum temperature boost'),
    ParamDef('ir_feedback','bool',False,   'PHYSICS_PARAMS', S_FEED, 'IR radiation feedback'),

    # AGN/Sink feedback
    ParamDef('n_sink',    'real', 0.0,     'PHYSICS_PARAMS', S_FEED, 'Sink formation density [H/cc]'),
    ParamDef('sink_AGN',  'bool', False,   'PHYSICS_PARAMS', S_FEED, 'AGN feedback from sinks',
             visible_when='sink==True'),
    ParamDef('bondi',     'bool', False,   'PHYSICS_PARAMS', S_FEED, 'Bondi accretion',
             visible_when='sink==True'),
    ParamDef('drag',      'bool', False,   'PHYSICS_PARAMS', S_FEED, 'Dynamical friction on BH',
             visible_when='sink==True'),
    ParamDef('spin_bh',   'bool', False,   'PHYSICS_PARAMS', S_FEED, 'BH spin evolution',
             visible_when='sink==True'),
    ParamDef('mad_jet',   'bool', False,   'PHYSICS_PARAMS', S_FEED, 'MAD jet model',
             visible_when='sink==True'),
    ParamDef('rAGN',      'real', 0.0,     'PHYSICS_PARAMS', S_FEED, 'AGN feedback radius [cells]',
             visible_when='sink_AGN==True'),
    ParamDef('eAGN_K',    'real', 0.0,     'PHYSICS_PARAMS', S_FEED, 'AGN kinetic efficiency',
             visible_when='sink_AGN==True'),
    ParamDef('eAGN_T',    'real', 0.0,     'PHYSICS_PARAMS', S_FEED, 'AGN thermal efficiency',
             visible_when='sink_AGN==True'),
    ParamDef('TAGN',      'real', 0.0,     'PHYSICS_PARAMS', S_FEED, 'AGN temperature threshold [K]',
             visible_when='sink_AGN==True'),
    ParamDef('X_floor',   'real', 0.0,     'PHYSICS_PARAMS', S_FEED, 'Accretion rate floor',
             visible_when='sink_AGN==True'),
    ParamDef('r_gal',     'real', 0.0,     'PHYSICS_PARAMS', S_FEED, 'Galaxy radius [kpc]',
             visible_when='sink_AGN==True'),
    ParamDef('boost_acc', 'real', 0.0,     'PHYSICS_PARAMS', S_FEED, 'Bondi accretion boost factor',
             visible_when='sink_AGN==True'),
    ParamDef('boost_drag','real', 0.0,     'PHYSICS_PARAMS', S_FEED, 'Drag boost factor',
             visible_when='sink_AGN==True'),
    ParamDef('Mseed',     'real', 0.0,     'PHYSICS_PARAMS', S_FEED, 'BH seed mass [Msun]',
             visible_when='sink==True'),
    ParamDef('T2maxAGN',  'real', 0.0,     'PHYSICS_PARAMS', S_FEED, 'AGN max temperature [K]',
             visible_when='sink_AGN==True'),
    ParamDef('vrel_merge','bool', False,   'PHYSICS_PARAMS', S_FEED, 'Relative velocity for BH merger',
             visible_when='sink==True'),
    ParamDef('rmerge',    'real', 0.0,     'PHYSICS_PARAMS', S_FEED, 'BH merger radius [cells]',
             visible_when='sink==True'),
    ParamDef('random_jet','bool', False,   'PHYSICS_PARAMS', S_FEED, 'Randomized jet direction',
             visible_when='sink_AGN==True'),

    # Units
    ParamDef('units_density','real',1.0,   'PHYSICS_PARAMS', S_COOL, 'Density unit [g/cm^3]'),
    ParamDef('units_time',   'real',1.0,   'PHYSICS_PARAMS', S_COOL, 'Time unit [s]'),
    ParamDef('units_length', 'real',1.0,   'PHYSICS_PARAMS', S_COOL, 'Length unit [cm]'),

    # Misc physics
    ParamDef('omega_b',   'real', 0.0,     'PHYSICS_PARAMS', S_COSMO,'Baryon fraction (physics_params)',
             visible_when='cosmo==True'),
    ParamDef('yieldtablefilename','str','', 'PHYSICS_PARAMS', S_FEED, 'Yield table filename'),

    # ====== CPL_PARAMS (Dark Energy) ======
    ParamDef('w0',   'real', -1.0, 'CPL_PARAMS', S_DE, 'DE equation of state w0',
             visible_when='cosmo==True'),
    ParamDef('wa',   'real',  0.0, 'CPL_PARAMS', S_DE, 'DE equation of state wa',
             visible_when='cosmo==True'),
    ParamDef('cs2_de','real', 1.0, 'CPL_PARAMS', S_DE, 'DE sound speed squared',
             visible_when='de_perturb==True'),
    ParamDef('de_table','str','',  'CPL_PARAMS', S_DE, 'DE table file',
             visible_when='de_perturb==True'),

    # ====== NEUTRINO_PARAMS ======
    ParamDef('omega_nu',     'real', 0.0,  'NEUTRINO_PARAMS', S_NU, 'Neutrino density parameter',
             visible_when='use_neutrino==True'),
    ParamDef('neutrino_table','str', '',    'NEUTRINO_PARAMS', S_NU, 'Neutrino transfer function table',
             visible_when='use_neutrino==True'),

    # ====== FR_PARAMS ======
    ParamDef('fR0',     'real', 1e-6,  'FR_PARAMS', S_MODGRAV, 'f(R) amplitude |fR0|',
             visible_when='use_fR==True'),
    ParamDef('fR_n',    'int',  1,     'FR_PARAMS', S_MODGRAV, 'f(R) power-law index',
             visible_when='use_fR==True'),
    ParamDef('n_iter_fR','int', 10,    'FR_PARAMS', S_MODGRAV, 'f(R) solver iterations',
             visible_when='use_fR==True'),
    ParamDef('fR_eps',  'real', 1e-6,  'FR_PARAMS', S_MODGRAV, 'f(R) solver tolerance',
             visible_when='use_fR==True'),

    # ====== NDGP_PARAMS ======
    ParamDef('omega_rc',   'real', 0.25, 'NDGP_PARAMS', S_MODGRAV, 'nDGP crossover parameter',
             visible_when='use_nDGP==True'),
    ParamDef('nDGP_branch','int',  1,    'NDGP_PARAMS', S_MODGRAV, 'nDGP branch (+1=normal,-1=self-acc)',
             visible_when='use_nDGP==True'),
    ParamDef('n_iter_nDGP','int',  10,   'NDGP_PARAMS', S_MODGRAV, 'nDGP solver iterations',
             visible_when='use_nDGP==True'),
    ParamDef('nDGP_eps',   'real', 1e-6, 'NDGP_PARAMS', S_MODGRAV, 'nDGP solver tolerance',
             visible_when='use_nDGP==True'),

    # ====== SYMMETRON_PARAMS ======
    ParamDef('a_ssb',       'real', 0.5,  'SYMMETRON_PARAMS', S_MODGRAV, 'Symmetron SSB scale factor',
             visible_when='use_symmetron==True'),
    ParamDef('beta_symmetron','real',1.0, 'SYMMETRON_PARAMS', S_MODGRAV, 'Symmetron coupling',
             visible_when='use_symmetron==True'),
    ParamDef('L_symmetron',  'real',1.0,  'SYMMETRON_PARAMS', S_MODGRAV, 'Symmetron Compton wavelength',
             visible_when='use_symmetron==True'),
    ParamDef('n_iter_symmetron','int',10, 'SYMMETRON_PARAMS', S_MODGRAV, 'Symmetron solver iterations',
             visible_when='use_symmetron==True'),
    ParamDef('symmetron_eps','real',1e-6, 'SYMMETRON_PARAMS', S_MODGRAV, 'Symmetron solver tolerance',
             visible_when='use_symmetron==True'),

    # ====== DILATON_PARAMS ======
    ParamDef('beta_dilaton','real',1.0,   'DILATON_PARAMS', S_MODGRAV, 'Dilaton coupling',
             visible_when='use_dilaton==True'),
    ParamDef('L_dilaton',   'real',1.0,   'DILATON_PARAMS', S_MODGRAV, 'Dilaton Compton wavelength',
             visible_when='use_dilaton==True'),
    ParamDef('a0_dilaton',  'real',0.5,   'DILATON_PARAMS', S_MODGRAV, 'Dilaton a0 parameter',
             visible_when='use_dilaton==True'),
    ParamDef('n_iter_dilaton','int',10,   'DILATON_PARAMS', S_MODGRAV, 'Dilaton solver iterations',
             visible_when='use_dilaton==True'),
    ParamDef('dilaton_eps', 'real',1e-6,  'DILATON_PARAMS', S_MODGRAV, 'Dilaton solver tolerance',
             visible_when='use_dilaton==True'),

    # ====== GALILEON_PARAMS ======
    ParamDef('c2_galileon','real',0.0,    'GALILEON_PARAMS', S_MODGRAV, 'Galileon c2 coefficient',
             visible_when='use_galileon==True'),
    ParamDef('c3_galileon','real',0.0,    'GALILEON_PARAMS', S_MODGRAV, 'Galileon c3 coefficient',
             visible_when='use_galileon==True'),
    ParamDef('n_iter_galileon','int',10,  'GALILEON_PARAMS', S_MODGRAV, 'Galileon solver iterations',
             visible_when='use_galileon==True'),
    ParamDef('galileon_eps','real',1e-6,  'GALILEON_PARAMS', S_MODGRAV, 'Galileon solver tolerance',
             visible_when='use_galileon==True'),

    # ====== COUPLED_DE_PARAMS ======
    ParamDef('beta_cde', 'real', 0.0,     'COUPLED_DE_PARAMS', S_DE, 'Coupled DE coupling constant',
             visible_when='use_coupled_de==True'),

    # ====== EDE_PARAMS ======
    ParamDef('omega_ede','real', 0.0,     'EDE_PARAMS', S_DE, 'Early dark energy density',
             visible_when='use_ede==True'),
    ParamDef('z_ede',    'real', 0.0,     'EDE_PARAMS', S_DE, 'EDE transition redshift',
             visible_when='use_ede==True'),
    ParamDef('w_ede',    'real', 0.0,     'EDE_PARAMS', S_DE, 'EDE equation of state',
             visible_when='use_ede==True'),

    # ====== MOND_PARAMS ======
    ParamDef('a0_mond',    'real', 1.2e-8,'MOND_PARAMS', S_MODGRAV, 'MOND acceleration scale [cm/s^2]',
             visible_when='use_mond==True'),
    ParamDef('mond_mu_type','int', 1,     'MOND_PARAMS', S_MODGRAV, 'MOND mu function type',
             visible_when='use_mond==True'),
    ParamDef('mond_type',  'int',  0,     'MOND_PARAMS', S_MODGRAV, 'MOND type (0=algebraic,1=QUMOND,2=AQUAL)',
             visible_when='use_mond==True'),
    ParamDef('n_iter_mond','int',  10,    'MOND_PARAMS', S_MODGRAV, 'MOND solver iterations',
             visible_when='use_mond==True'),
    ParamDef('mond_eps',   'real', 1e-6,  'MOND_PARAMS', S_MODGRAV, 'MOND solver tolerance',
             visible_when='use_mond==True'),
    ParamDef('g_ext_mond', 'real_arr','0,0,0','MOND_PARAMS', S_MODGRAV, 'External field [cm/s^2] (3 components)',
             visible_when='use_mond==True'),

    # ====== SGS_PARAMS ======
    ParamDef('sgs_C_prod','real', 0.04,   'SGS_PARAMS', S_SGS, 'SGS production coefficient',
             visible_when='use_sgs==True'),
    ParamDef('sgs_C_diss','real', 0.5,    'SGS_PARAMS', S_SGS, 'SGS dissipation coefficient',
             visible_when='use_sgs==True'),
    ParamDef('sgs_C_smag','real', 0.18,   'SGS_PARAMS', S_SGS, 'SGS Smagorinsky coefficient',
             visible_when='use_sgs==True'),
    ParamDef('sgs_floor', 'real', 1e-30,  'SGS_PARAMS', S_SGS, 'SGS energy floor',
             visible_when='use_sgs==True'),
    ParamDef('sgs_cap',   'real', 1e30,   'SGS_PARAMS', S_SGS, 'SGS energy cap',
             visible_when='use_sgs==True'),
    ParamDef('sgs_e_init','real', 0.0,    'SGS_PARAMS', S_SGS, 'SGS initial energy',
             visible_when='use_sgs==True'),
    ParamDef('sgs_hydro', 'bool', True,   'SGS_PARAMS', S_SGS, 'SGS affects hydro',
             visible_when='use_sgs==True'),

    # ====== SIDM_PARAMS ======
    ParamDef('sidm_cross_section','real',1.0,'SIDM_PARAMS',S_SIDM,'SIDM cross section [cm^2/g]',
             visible_when='sidm==True'),
    ParamDef('sidm_npart_min','int',3,    'SIDM_PARAMS', S_SIDM, 'Min particles per cell for scattering',
             visible_when='sidm==True'),
    ParamDef('sidm_type',  'int',  0,     'SIDM_PARAMS', S_SIDM, 'SIDM type (0=const,1=vel-dep,2=yukawa)',
             visible_when='sidm==True'),
    ParamDef('sidm_v0',    'real', 0.0,   'SIDM_PARAMS', S_SIDM, 'SIDM velocity scale [km/s]',
             visible_when='sidm==True'),
    ParamDef('sidm_power', 'real', 0.0,   'SIDM_PARAMS', S_SIDM, 'SIDM velocity power-law index',
             visible_when='sidm==True'),
    ParamDef('sidm_courant','real',0.1,   'SIDM_PARAMS', S_SIDM, 'SIDM timestep Courant factor',
             visible_when='sidm==True'),
    ParamDef('sidm_angular','bool',False, 'SIDM_PARAMS', S_SIDM, 'Anisotropic scattering',
             visible_when='sidm==True'),
    ParamDef('sidm_epsilon','real',0.0,   'SIDM_PARAMS', S_SIDM, 'Anisotropy parameter',
             visible_when='sidm==True'),
    ParamDef('sidm_inelastic','bool',False,'SIDM_PARAMS',S_SIDM, 'Inelastic scattering (iSIDM)',
             visible_when='sidm==True'),
    ParamDef('sidm_delta', 'real', 0.0,   'SIDM_PARAMS', S_SIDM, 'Inelastic mass splitting [keV]',
             visible_when='sidm==True'),
    ParamDef('sidm_frac_excited','real',0.0,'SIDM_PARAMS',S_SIDM,'Initial excited state fraction',
             visible_when='sidm==True'),
]

# Build lookup dict (case-insensitive keys)
PARAM_BY_NAME = OrderedDict()
for _p in PARAMS:
    PARAM_BY_NAME[_p.name.lower()] = _p

# ============================================================================
# Presets
# ============================================================================
PRESETS = OrderedDict([
    ('cosmo_uniform', {
        'desc': 'Cosmological uniform box (no zoom-in)',
        'values': {
            'cosmo': True, 'pic': True, 'poisson': True, 'hydro': True,
            'levelmin': 9, 'levelmax': 9, 'nexpand': 1,
            'filetype': 'grafic', 'ordering': 'hilbert',
            'gamma': 1.6666667, 'courant_factor': 0.8,
            'scheme': 'muscl', 'slope_type': 1,
        }
    }),
    ('cosmo_zoomin', {
        'desc': 'Cosmological zoom-in with AMR + SF + cooling',
        'values': {
            'cosmo': True, 'pic': True, 'poisson': True, 'hydro': True,
            'sink': True,
            'levelmin': 9, 'levelmax': 16, 'nexpand': 1,
            'filetype': 'grafic', 'ordering': 'ksection',
            'memory_balance': True, 'use_fftw': True,
            'nsubcycle': '1,1,2',
            'gamma': 1.6666667, 'courant_factor': 0.8,
            'scheme': 'muscl', 'slope_type': 1,
            'm_refine': '9*8.',
            'q_refine_holdback': True,
            'cooling': True, 'haardt_madau': True, 'metal': True,
            'eps_star': 0.02, 'del_star': 55.0, 'n_star': 0.1,
            'T2_star': 750, 'g_star': 1.66667, 'jeans_ncells': 4,
        }
    }),
    ('restart', {
        'desc': 'Restart from existing output (nrestart>0)',
        'values': {
            'nrestart': 1,
        }
    }),
    ('gravity_only', {
        'desc': 'N-body only (no hydro)',
        'values': {
            'cosmo': True, 'pic': True, 'poisson': True, 'hydro': False,
            'levelmin': 9, 'levelmax': 9,
            'filetype': 'grafic', 'ordering': 'hilbert',
        }
    }),
    ('hr5_production', {
        'desc': 'HR5 production run (full physics)',
        'values': {
            'cosmo': True, 'pic': True, 'poisson': True, 'hydro': True,
            'sink': True,
            'levelmin': 9, 'levelmax': 16, 'nexpand': 1,
            'filetype': 'grafic', 'ordering': 'ksection',
            'memory_balance': True, 'use_fftw': True,
            'nsubcycle': '1,1,2', 'ncontrol': 1, 'nremap': 5,
            'gamma': 1.6666667, 'courant_factor': 0.8,
            'scheme': 'muscl', 'slope_type': 1,
            'm_refine': '9*8.', 'ivar_refine': 0,
            'interpol_var': 1, 'interpol_type': 0,
            'q_refine_holdback': True, 'dr_proper': 2.27,
            'jeans_ncells': 4,
            'cooling': True, 'haardt_madau': True, 'metal': True,
            'cooling_method': 'exact',
            'grackle_table': 'grackle_multi_z.bin',
            'eps_star': 0.02, 'del_star': 55.0, 'n_star': 0.1,
            'm_star': 0.2, 'T2_star': 750, 'g_star': 1.66667,
            'z_ave': 1.0e-3, 'z_reion': 10.0,
            'f_w': 3.0, 'f_ek': 0.3, 'delayed_cooling': True,
            'eps_sn1': 2.0, 'eps_sn2': 2.0, 'tol': 1e-3,
            'Mseed': 1e4,
            'sink_AGN': True, 'bondi': True, 'drag': True,
            'rAGN': 4.0, 'X_floor': 1e-2,
            'eAGN_K': 1.0, 'eAGN_T': 0.15, 'TAGN': 0.0,
            'r_gal': 50.0, 'T2maxAGN': 1e8,
            'boost_acc': 2.0, 'boost_drag': 2.0,
            'vrel_merge': True, 'rmerge': 4.0,
            'spin_bh': True, 'mad_jet': True,
            'omega_b': 0.04,
            'yieldtablefilename': 'yield_table.asc',
        }
    }),
])

# ============================================================================
# Validation
# ============================================================================
class ValidationMsg:
    """A validation message (error or warning)."""
    def __init__(self, level, msg):
        self.level = level  # 'ERROR' or 'WARNING'
        self.msg = msg
    def __str__(self):
        return '[{}] {}'.format(self.level, self.msg)

def validate_params(values):
    """Run all validation rules. Returns list of ValidationMsg."""
    values = _normalize_values(values)
    msgs = []

    # --- Mutual exclusion: modified gravity (max 1) ---
    mg_flags = ['use_fr', 'use_ndgp', 'use_mond',
                'use_symmetron', 'use_dilaton', 'use_galileon']
    active_mg = [f for f in mg_flags if values.get(f, False)]
    if len(active_mg) > 1:
        msgs.append(ValidationMsg('ERROR',
            'Modified gravity mutual exclusion: only 1 allowed, got: {}'.format(
                ', '.join(active_mg))))

    # --- Dependencies ---
    if values.get('hydro') and not values.get('gamma'):
        msgs.append(ValidationMsg('WARNING', 'hydro=True but gamma not set'))

    if values.get('cooling') and not values.get('z_reion'):
        msgs.append(ValidationMsg('WARNING', 'cooling=True but z_reion not set'))

    if values.get('sink') and not values.get('mseed'):
        msgs.append(ValidationMsg('WARNING', 'sink=True but Mseed not set'))

    dr = values.get('dr_proper', 0)
    if isinstance(dr, (int, float)) and dr > 0:
        jn = values.get('jeans_ncells', 0)
        if not jn or (isinstance(jn, (int, float)) and jn <= 0):
            msgs.append(ValidationMsg('WARNING',
                'FPR active (dr_proper={}) but jeans_ncells not set (eEOS floor needed)'.format(dr)))

    if values.get('use_sgs') and not values.get('hydro'):
        msgs.append(ValidationMsg('ERROR', 'use_sgs=True requires hydro=True (also build with NVAR=12)'))

    if values.get('sidm') and not values.get('pic'):
        msgs.append(ValidationMsg('ERROR', 'sidm=True requires pic=True'))

    if values.get('use_neutrino') and not values.get('use_fftw'):
        msgs.append(ValidationMsg('WARNING', 'use_neutrino=True typically requires use_fftw=True'))

    # --- Range checks ---
    cf = values.get('courant_factor')
    if cf is not None and isinstance(cf, (int, float)):
        if cf <= 0 or cf > 1:
            msgs.append(ValidationMsg('ERROR',
                'courant_factor={} out of range (0,1]'.format(cf)))

    lmin = values.get('levelmin', 1)
    lmax = values.get('levelmax', 1)
    if isinstance(lmin, int) and isinstance(lmax, int):
        if lmin < 1 or lmin > 30:
            msgs.append(ValidationMsg('ERROR', 'levelmin={} out of range [1,30]'.format(lmin)))
        if lmax < 1 or lmax > 30:
            msgs.append(ValidationMsg('ERROR', 'levelmax={} out of range [1,30]'.format(lmax)))
        if lmax < lmin:
            msgs.append(ValidationMsg('ERROR',
                'levelmax={} < levelmin={}'.format(lmax, lmin)))

    return msgs

# ============================================================================
# Fortran Namelist Formatter
# ============================================================================

# Canonical group ordering
GROUP_ORDER = [
    'RUN_PARAMS', 'COSMO_PARAMS', 'OUTPUT_PARAMS', 'INIT_PARAMS',
    'AMR_PARAMS', 'LIGHTCONE_PARAMS', 'REFINE_PARAMS',
    'HYDRO_PARAMS', 'POISSON_PARAMS', 'PHYSICS_PARAMS',
    # Optional groups (only emitted when relevant)
    'CPL_PARAMS', 'NEUTRINO_PARAMS',
    'FR_PARAMS', 'NDGP_PARAMS', 'SYMMETRON_PARAMS',
    'DILATON_PARAMS', 'GALILEON_PARAMS',
    'COUPLED_DE_PARAMS', 'EDE_PARAMS',
    'MOND_PARAMS', 'SGS_PARAMS', 'SIDM_PARAMS',
]

# Groups that are always emitted (even if empty)
ALWAYS_EMIT = {
    'RUN_PARAMS', 'OUTPUT_PARAMS', 'INIT_PARAMS', 'AMR_PARAMS',
    'LIGHTCONE_PARAMS', 'REFINE_PARAMS', 'HYDRO_PARAMS',
    'POISSON_PARAMS', 'PHYSICS_PARAMS',
}

# Groups that require a toggle to be True
GROUP_REQUIRES = {
    'COSMO_PARAMS': 'cosmo',
    'CPL_PARAMS': None,  # emitted if w0 or wa differ from defaults
    'NEUTRINO_PARAMS': 'use_neutrino',
    'FR_PARAMS': 'use_fR',
    'NDGP_PARAMS': 'use_nDGP',
    'SYMMETRON_PARAMS': 'use_symmetron',
    'DILATON_PARAMS': 'use_dilaton',
    'GALILEON_PARAMS': 'use_galileon',
    'COUPLED_DE_PARAMS': 'use_coupled_de',
    'EDE_PARAMS': 'use_ede',
    'MOND_PARAMS': 'use_mond',
    'SGS_PARAMS': 'use_sgs',
    'SIDM_PARAMS': 'sidm',
}


def _fmt_fortran_value(val, ftype):
    """Convert Python value to Fortran namelist string."""
    if ftype == 'bool':
        if isinstance(val, str):
            val = val.strip().lower() in ('true', '.true.', '1', 'yes', 't')
        return '.true.' if val else '.false.'
    elif ftype == 'str':
        s = str(val).strip()
        if not s:
            return None  # skip empty strings
        # Don't double-quote if already quoted
        if s.startswith("'") and s.endswith("'"):
            return s
        return "'{}'".format(s)
    elif ftype in ('int_arr', 'real_arr', 'str_arr'):
        s = str(val).strip()
        if not s:
            return None  # skip empty arrays
        return s
    elif ftype == 'int':
        if isinstance(val, float):
            val = int(val)
        return str(int(val))
    elif ftype == 'real':
        return _fmt_real(val)
    return str(val)


def _fmt_real(val):
    """Format a real number for Fortran, preserving Fortran-style notation."""
    if isinstance(val, str):
        val = val.strip()
        # Already in Fortran notation? Pass through
        if 'd' in val.lower() or 'e' in val.lower():
            return val
        try:
            val = float(val)
        except ValueError:
            return val
    if isinstance(val, int):
        val = float(val)
    if val == 0.0:
        return '0.0'
    # For very large/small values, use Fortran d notation
    av = abs(val)
    if av >= 1e6 or av < 1e-3:
        s = '{:.6e}'.format(val)
        # Clean up trailing zeros in mantissa
        m, e = s.split('e')
        m = m.rstrip('0').rstrip('.')
        if '.' not in m:
            m += '.0'
        exp = int(e)
        return '{}d{}'.format(m, exp)
    else:
        # Plain float
        s = '{:.10g}'.format(val)
        if '.' not in s:
            s += '.'
        return s


def _normalize_values(values):
    """Normalize all keys to lowercase for consistent lookup."""
    return OrderedDict((k.lower(), v) for k, v in values.items())


def format_namelist(values, emit_defaults=False):
    """Generate complete Fortran namelist string from values dict."""
    values = _normalize_values(values)
    # Collect parameters by group
    group_params = OrderedDict()
    for g in GROUP_ORDER:
        group_params[g] = []

    for p in PARAMS:
        if p.group not in group_params:
            group_params[p.group] = []
        val = values.get(p.name.lower())
        if val is None:
            if emit_defaults and p.default is not None:
                val = p.default
            else:
                continue
        # Skip default values unless emit_defaults
        if not emit_defaults:
            if val == p.default:
                continue
        # Format the value
        fval = _fmt_fortran_value(val, p.ftype)
        if fval is None:
            continue
        group_params[p.group].append((p.name, fval))

    # Build output
    lines = []
    first = True
    for g in GROUP_ORDER:
        params = group_params.get(g, [])

        # Should we emit this group?
        if g in ALWAYS_EMIT:
            emit = True
        elif g in GROUP_REQUIRES:
            req = GROUP_REQUIRES[g]
            if req is None:
                # Emit only if there are params to write
                emit = len(params) > 0
            else:
                emit = bool(values.get(req.lower(), False))
        else:
            emit = len(params) > 0

        if not emit:
            continue

        if not first:
            lines.append('')
        first = False

        lines.append('&{}'.format(g))
        for name, fval in params:
            lines.append('{}={}'.format(name, fval))
        lines.append('/')

    return '\n'.join(lines) + '\n'

# ============================================================================
# Namelist Parser (Import)
# ============================================================================

def parse_namelist(text):
    """Parse a Fortran namelist file into a dict of {param_name: raw_string_value}.
    Returns (values_dict, group_order_list)."""
    values = OrderedDict()
    groups_seen = []

    # Remove comments (but preserve strings)
    cleaned_lines = []
    for line in text.split('\n'):
        # Simple comment removal: strip everything after ! not inside quotes
        in_str = False
        result = []
        for ch in line:
            if ch == "'" and not in_str:
                in_str = True
            elif ch == "'" and in_str:
                in_str = False
            elif ch == '!' and not in_str:
                break
            result.append(ch)
        cleaned_lines.append(''.join(result))

    content = '\n'.join(cleaned_lines)

    # Find all namelist groups
    pattern = r'&(\w+)\s*(.*?)\s*/'
    for m in re.finditer(pattern, content, re.DOTALL):
        group_name = m.group(1).upper()
        body = m.group(2)
        groups_seen.append(group_name)

        # Parse assignments: name=value (possibly multi-line, comma-separated values)
        # We need to handle things like: m_refine=9*8.,\nivar_refine=0
        # Strategy: split on pattern `name=` but be careful with strings
        assignments = _parse_assignments(body)
        for name, raw_val in assignments:
            values[name.strip().lower()] = raw_val.strip()

    return values, groups_seen


def _parse_assignments(body):
    """Parse 'name=value' pairs from a namelist group body."""
    assignments = []
    # Tokenize: find all `identifier = ` patterns
    # We look for pattern: word = value, where value extends until next word= or end
    tokens = []
    i = 0
    chars = body
    n = len(chars)

    while i < n:
        # Skip whitespace
        while i < n and chars[i] in ' \t\n\r':
            i += 1
        if i >= n:
            break

        # Try to match identifier=
        m = re.match(r'([a-zA-Z_]\w*(?:\(\d+\))?)\s*=', chars[i:])
        if m:
            name = m.group(1)
            i += m.end()
            # Now collect value until next `identifier=` or end
            val_start = i
            while i < n:
                # Check if we're at a new assignment
                m2 = re.match(r'\s+([a-zA-Z_]\w*(?:\(\d+\))?)\s*=', chars[i:])
                if m2:
                    # Check it's not inside a string
                    # Simple heuristic: count quotes before this position
                    preceding = chars[val_start:i]
                    if preceding.count("'") % 2 == 0:
                        break
                i += 1
            val = chars[val_start:i].strip().rstrip(',').strip()
            assignments.append((name, val))
        else:
            i += 1

    return assignments


def import_to_values(raw_parsed):
    """Convert raw parsed string values to typed Python values using PARAM_BY_NAME."""
    values = OrderedDict()
    for name, raw in raw_parsed.items():
        pdef = PARAM_BY_NAME.get(name)
        if pdef is None:
            # Unknown parameter - keep as string
            values[name] = raw
            continue
        values[name] = _parse_value(raw, pdef.ftype)
    return values


def _parse_value(raw, ftype):
    """Parse a raw string value into a Python type."""
    raw = raw.strip()
    if ftype == 'bool':
        return raw.lower() in ('.true.', 'true', '1', 't', 'yes')
    elif ftype == 'int':
        try:
            return int(float(raw))
        except (ValueError, TypeError):
            return raw
    elif ftype == 'real':
        # Handle Fortran d notation
        raw_clean = raw.lower().replace('d', 'e')
        try:
            return float(raw_clean)
        except (ValueError, TypeError):
            return raw
    elif ftype == 'str':
        return raw.strip("'\"")
    elif ftype in ('int_arr', 'real_arr', 'str_arr'):
        return raw
    return raw

# ============================================================================
# Interactive TUI
# ============================================================================

def get_sections():
    """Get ordered list of unique sections."""
    seen = set()
    sections = []
    for p in PARAMS:
        if p.section not in seen:
            seen.add(p.section)
            sections.append(p.section)
    return sections


def _eval_condition(cond, values):
    """Evaluate a simple visibility condition against current values."""
    if not cond:
        return True
    cond = cond.strip()

    # Handle simple comparisons: "param==value", "param>0", "param=='str'"
    for op in ['==', '!=', '>=', '<=', '>', '<']:
        if op in cond:
            parts = cond.split(op, 1)
            if len(parts) == 2:
                param_name = parts[0].strip()
                expected = parts[1].strip()
                actual = values.get(param_name.lower())
                if actual is None:
                    # Check defaults
                    pdef = PARAM_BY_NAME.get(param_name.lower())
                    if pdef:
                        actual = pdef.default
                    else:
                        return False

                # Parse expected value
                if expected.lower() in ('true', '.true.'):
                    expected_val = True
                elif expected.lower() in ('false', '.false.'):
                    expected_val = False
                elif expected.startswith("'") and expected.endswith("'"):
                    expected_val = expected.strip("'")
                else:
                    try:
                        expected_val = float(expected)
                        if expected_val == int(expected_val):
                            expected_val = int(expected_val)
                    except ValueError:
                        expected_val = expected

                if op == '==':
                    return actual == expected_val
                elif op == '!=':
                    return actual != expected_val
                elif op == '>':
                    try:
                        return float(actual) > float(expected_val)
                    except (ValueError, TypeError):
                        return False
                elif op == '<':
                    try:
                        return float(actual) < float(expected_val)
                    except (ValueError, TypeError):
                        return False
                elif op == '>=':
                    try:
                        return float(actual) >= float(expected_val)
                    except (ValueError, TypeError):
                        return False
                elif op == '<=':
                    try:
                        return float(actual) <= float(expected_val)
                    except (ValueError, TypeError):
                        return False
    return True


def interactive_edit(values):
    """Interactive section-by-section parameter editing."""
    sections = get_sections()

    print('\n--- Editing parameters (section by section) ---')
    print("  Type 'param=value' to change, Enter to keep, 'n' for next section, 'q' to finish\n")

    for section in sections:
        # Collect visible params for this section
        visible = []
        for p in PARAMS:
            if p.section != section:
                continue
            if not _eval_condition(p.visible_when, values):
                continue
            visible.append(p)

        if not visible:
            continue

        print('--- {} ---'.format(section))
        for p in visible:
            val = values.get(p.name.lower())
            if val is None:
                val = p.default
            # Display
            display_val = _fmt_display(val, p.ftype)
            desc = '  # {}'.format(p.desc) if p.desc else ''
            if p.choices:
                desc += '  [{}]'.format('|'.join(str(c) for c in p.choices))
            print('  {:<25s} = {}{}'.format(p.name, display_val, desc))

        while True:
            try:
                inp = input('> ').strip()
            except (EOFError, KeyboardInterrupt):
                print()
                return values

            if inp == '' or inp.lower() == 'n':
                break
            if inp.lower() == 'q':
                return values

            # Parse param=value
            if '=' not in inp:
                print("  (Use 'param=value' format, 'n' for next, 'q' to finish)")
                continue

            eq_pos = inp.index('=')
            pname = inp[:eq_pos].strip().lower()
            pval = inp[eq_pos+1:].strip()

            # Find param definition
            pdef = PARAM_BY_NAME.get(pname)
            if pdef is None:
                # Allow arbitrary params
                values[pname] = pval
                print('  {} = {} (custom parameter)'.format(pname, pval))
                continue

            # Validate choices
            if pdef.choices:
                clean_val = pval.strip("'\"")
                if clean_val not in [str(c) for c in pdef.choices]:
                    print('  Invalid choice. Options: {}'.format(pdef.choices))
                    continue

            # Parse and set
            typed_val = _parse_value(pval, pdef.ftype)
            values[pname] = typed_val
            print('  {} = {}'.format(pname, _fmt_display(typed_val, pdef.ftype)))

        print()

    return values


def _fmt_display(val, ftype):
    """Format a value for display."""
    if ftype == 'bool':
        if isinstance(val, str):
            val = val.strip().lower() in ('true', '.true.', '1', 'yes', 't')
        return '.true.' if val else '.false.'
    elif ftype == 'str':
        return "'{}'".format(val) if val else "''"
    elif ftype in ('int_arr', 'real_arr', 'str_arr'):
        return str(val) if val else '(empty)'
    elif ftype == 'real':
        return _fmt_real(val) if val is not None else '0.0'
    elif ftype == 'int':
        return str(int(val)) if val is not None else '0'
    return str(val)


def show_review(values):
    """Show final namelist review before writing."""
    print('\n' + '='*60)
    print('FINAL NAMELIST REVIEW')
    print('='*60)
    output = format_namelist(values)
    print(output)

    # Run validation
    msgs = validate_params(values)
    if msgs:
        print('--- Validation Messages ---')
        for m in msgs:
            print('  {}'.format(m))
        print()

    return output

# ============================================================================
# CLI
# ============================================================================

def list_presets():
    """Print available presets."""
    print('\nAvailable presets:')
    for i, (name, info) in enumerate(PRESETS.items(), 1):
        print('  [{}] {:<20s} - {}'.format(i, name, info['desc']))
    print()


def main():
    parser = argparse.ArgumentParser(
        description='RAMSES Namelist Generator',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""Examples:
  %(prog)s                                         # Interactive mode
  %(prog)s --preset hr5_production -o cosmo.nml     # Non-interactive
  %(prog)s --import old.nml -o new.nml              # Import & re-export
  %(prog)s --import old.nml                         # Import & interactive edit
  %(prog)s --list-presets                            # Show presets
""")
    parser.add_argument('--preset', '-p', help='Use preset (name or number)')
    parser.add_argument('--output', '-o', help='Output namelist file')
    parser.add_argument('--import-nml', '--import', dest='import_nml',
                        help='Import existing namelist file')
    parser.add_argument('--list-presets', action='store_true',
                        help='List available presets')
    parser.add_argument('--non-interactive', '-n', action='store_true',
                        help='Skip interactive editing')
    parser.add_argument('--emit-defaults', action='store_true',
                        help='Emit all parameters including defaults')

    args = parser.parse_args()

    if args.list_presets:
        list_presets()
        return 0

    values = OrderedDict()

    # --- Import existing namelist ---
    if args.import_nml:
        path = args.import_nml
        if not os.path.isfile(path):
            print('Error: file not found: {}'.format(path), file=sys.stderr)
            return 1
        with open(path, 'r') as f:
            text = f.read()
        raw_parsed, _ = parse_namelist(text)
        values = import_to_values(raw_parsed)
        print('Imported {} parameters from {}'.format(len(values), path))

    # --- Apply preset ---
    if args.preset:
        preset_name = args.preset
        # Allow numeric selection
        if preset_name.isdigit():
            idx = int(preset_name) - 1
            names = list(PRESETS.keys())
            if 0 <= idx < len(names):
                preset_name = names[idx]
            else:
                print('Error: preset index {} out of range'.format(args.preset),
                      file=sys.stderr)
                return 1
        if preset_name not in PRESETS:
            print('Error: unknown preset "{}". Use --list-presets'.format(preset_name),
                  file=sys.stderr)
            return 1
        preset = PRESETS[preset_name]
        print('Using preset: {} - {}'.format(preset_name, preset['desc']))
        for k, v in preset['values'].items():
            values[k.lower()] = v

    # --- Interactive mode ---
    is_interactive = not args.non_interactive and sys.stdin.isatty()

    if is_interactive and not args.preset and not args.import_nml:
        # Show preset selection
        print('\n=== RAMSES Namelist Generator ===')
        print('Select preset:')
        for i, (name, info) in enumerate(PRESETS.items(), 1):
            print('  [{}] {:<20s} - {}'.format(i, name, info['desc']))
        print('  [{}] {:<20s} - {}'.format(len(PRESETS)+1, 'blank', 'Start from scratch'))

        try:
            choice = input('> ').strip()
        except (EOFError, KeyboardInterrupt):
            print()
            return 0

        if choice.isdigit():
            idx = int(choice) - 1
            names = list(PRESETS.keys())
            if 0 <= idx < len(names):
                preset_name = names[idx]
                preset = PRESETS[preset_name]
                print('Using preset: {} - {}'.format(preset_name, preset['desc']))
                for k, v in preset['values'].items():
                    values[k] = v
            elif idx == len(names):
                print('Starting from scratch...')
            else:
                print('Invalid choice, starting from scratch...')
        elif choice in PRESETS:
            preset = PRESETS[choice]
            print('Using preset: {} - {}'.format(choice, preset['desc']))
            for k, v in preset['values'].items():
                values[k] = v

    if is_interactive and not args.non_interactive:
        values = interactive_edit(values)

    # --- Generate output ---
    output = show_review(values)

    # --- Write file ---
    if args.output:
        outpath = args.output
    elif is_interactive:
        default_name = 'cosmo.nml'
        try:
            inp = input('Write to file [{}]: '.format(default_name)).strip()
        except (EOFError, KeyboardInterrupt):
            print()
            return 0
        outpath = inp if inp else default_name
    else:
        # Non-interactive, no output specified -> stdout
        return 0

    with open(outpath, 'w') as f:
        f.write(output)
    print('Written: {}'.format(outpath))
    return 0


if __name__ == '__main__':
    sys.exit(main() or 0)
