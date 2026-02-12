module amr_parameters

  ! Define real types
  integer,parameter::sp=kind(1.0E0)
#ifndef NPRE
  integer,parameter::dp=kind(1.0E0) ! default
#else
#if NPRE==4
  integer,parameter::dp=kind(1.0E0) ! real*4
#else
  integer,parameter::dp=kind(1.0D0) ! real*8
#endif
#endif
#ifdef QUADHILBERT
  integer,parameter::qdp=kind(1.0_16) ! real*16
#else
  integer,parameter::qdp=kind(1.0_8) ! real*8
#endif
  integer,parameter::MAXOUT=1000
  integer,parameter::MAXLEVEL=100
  
  ! Define integer types (for particle IDs mostly)
  integer,parameter::i4b=4
#ifndef LONGINT
  integer,parameter::i8b=4  ! default long int are short int
#else
  integer,parameter::i8b=8  ! long int are long int
#endif

  ! Number of dimensions
#ifndef NDIM
  integer,parameter::ndim=1
#else
  integer,parameter::ndim=NDIM
#endif
  integer,parameter::twotondim=2**ndim
  integer,parameter::threetondim=3**ndim
  integer,parameter::twondim=2*ndim

  ! Vectorization parameter
#ifndef NVECTOR
  integer,parameter::nvector=500  ! Size of vector sweeps
#else
  integer,parameter::nvector=NVECTOR
#endif

  integer, parameter :: nstride = 65536

  ! Run control
  logical::verbose =.false.   ! Write everything
  logical::hydro   =.false.   ! Hydro activated
  logical::pic     =.false.   ! Particle In Cell activated
  logical::poisson =.false.   ! Poisson solver activated
  logical::cosmo   =.false.   ! Cosmology activated
  logical::star    =.false.   ! Star formation activated
  logical::sink    =.false.   ! Sink particles activated
  logical::sinkprops=.false.  ! Write sink properties at each coarse step
  logical::rt      =.false.   ! Radiative transfer activated
  logical::debug   =.false.   ! Debug mode activated
  logical::static  =.false.   ! Static mode activated
  logical::tracer  =.false.   ! Tracer particles activated
  logical::lightcone=.false.  ! Enable lightcone generation
  logical::clumpfind=.false.  ! Enable clump finder
  logical::aton=.false.       ! Enable ATON coarse grid radiation transfer

  ! Mesh parameters
  integer::geom=1             ! 1: cartesian, 2: cylindrical, 3: spherical
  integer::nx=1,ny=1,nz=1     ! Number of coarse cells in each dimension
  integer::levelmin=1         ! Full refinement up to levelmin
  integer::nlevelmax=1        ! Maximum number of level
  integer::ngridmax=0         ! Maximum number of grids
  integer,dimension(1:MAXLEVEL)::nexpand=1 ! Number of mesh expansion
  integer::nexpand_bound=1    ! Number of mesh expansion for virtual boundaries
  real(dp)::boxlen=1.0D0      ! Box length along x direction
  character(len=128)::ordering='hilbert'
  logical::cost_weighting=.true. ! Activate load balancing according to cpu time
  ! Recursive bisection tree parameters
  integer::nbilevelmax=1      ! Max steps of bisection partitioning
  integer::nbinodes=3         ! Max number of internal nodes
  integer::nbileafnodes=2     ! Max number of leaf (terminal) nodes
  real(dp)::bisec_tol=0.05d0  ! Tolerance for bisection load balancing

  ! Step parameters
  integer::nrestart=0         ! New run or backup file number
!jhshin1
  integer::nrestart_quad=0    ! Restart with double precision Hilbert keys
  logical::restart_remap=.false. ! Force load balance on restart
!jhshin2
  integer::nstepmax=1000000   ! Maximum number of time steps
  integer::ncontrol=1         ! Write control variables
  integer::fbackup=1000000    ! Backup data to disk
  integer::nremap=0           ! Load balancing frequency (0: never)

  ! Output parameters
  integer::iout=1             ! Increment for output times
  integer::ifout=1            ! Increment for output files
  integer::iback=1            ! Increment for backup files
  integer::noutput=1          ! Total number of outputs
  integer::foutput=1000000    ! Frequency of outputs
  integer::output_mode=0      ! Output mode (for hires runs)
  logical::gadget_output=.false. ! Output in gadget format
  logical::output_now=.false. ! write output next step

  ! Lightcone parameters
  real(dp)::thetay_cone=12.5
  real(dp)::thetaz_cone=12.5
  real(dp)::zmax_cone=2.0

  ! Cosmology and physical parameters
  real(dp)::boxlen_ini        ! Box size in h-1 Mpc
  real(dp)::omega_b=0.0D0     ! Omega Baryon
  real(dp)::omega_m=1.0D0     ! Omega Matter
  real(dp)::omega_l=0.0D0     ! Omega Lambda
  real(dp)::omega_k=0.0D0     ! Omega Curvature
  real(dp)::h0     =1.0D0     ! Hubble constant in km/s/Mpc
  real(dp)::aexp   =1.0D0     ! Current expansion factor
  real(dp)::hexp   =0.0D0     ! Current Hubble parameter
  real(dp)::texp   =0.0D0     ! Current proper time
!jhshin1
  real(dp)::n_sink =-1D0      ! Sink particle density threshold in H/cc
!jhshin2
  real(dp)::rho_sink = -1.D0  ! Sink particle density threshold in g/cc
  real(dp)::d_sink = -1.D0    ! Sink particle density threshold in user units
  real(dp)::m_star =-1.0      ! Star particle mass in units of mass_sph
  real(dp)::n_star =0.1D0     ! Star formation density threshold in H/cc
  real(dp)::t_star =0.0D0     ! Star formation time scale in Gyr
  real(dp)::eps_star=0.0D0    ! Star formation efficiency (0.02 at n_star=0.1 gives t_star=8 Gyr)
  real(dp)::T2_star=0.0D0     ! Typical ISM polytropic temperature
  real(dp)::g_star =1.6D0     ! Typical ISM polytropic index
  real(dp)::jeans_ncells=-1   ! Jeans polytropic EOS
  real(dp)::del_star=2.D2     ! Minimum overdensity to define ISM
  real(dp)::f_ek   =1.0D0     ! Supernovae kinetic energy fraction (only between 0 and 1)
  real(dp)::rbubble=0.0D0     ! Supernovae superbubble radius in pc
  real(dp)::f_w    =0.0D0     ! Supernovae mass loading factor
  integer ::ndebris=1         ! Supernovae debris particle number
  real(dp)::mass_gmc=-1.0     ! Stochastic exploding GMC mass
  real(dp)::z_ave  =0.0D0     ! Average metal abundance
  real(dp)::B_ave  =0.0D0     ! Average magnetic field
  real(dp)::z_reion=8.5D0     ! Reionization redshift
  real(dp)::T2_start          ! Starting gas temperature
!jhshin1
  real(dp)::T2max= 1d50       ! Temperature ceiling for cooling_fine
!jhshin2
  real(dp)::t_diss =20.0D0    ! Dissipation timescale for feedback
  real(dp)::t_sne =10.0D0     ! Supernova blast time
  real(dp)::J21    =0.0D0     ! UV flux at threshold in 10^21 units
  real(dp)::a_spec =1.0D0     ! Slope of the UV spectrum
  real(dp)::beta_fix=0.0D0    ! Pressure fix parameter
  real(dp)::kappa_IR=0d0      ! IR dust opacity
  real(dp)::ind_rsink=4.0d0   ! Number of cells defining the radius of the sphere where AGN feedback is active
  real(dp)::ir_eff=0.75       ! efficiency of the IR feedback (only when ir_feedback=.true.)
!jhshin1
  real(dp)::sf_trelax=0.0D0  ! Relaxation time for star formation (cosmo=.false.only)
  integer::sf_model=3         ! Virial star formation model
!jhshin2
  real(dp)::f_ekAGN=1.0D0     ! AGN kinetic energy fraction (only between 0 and 1)
  real(dp)::rAGN   =0.0D0     ! AGN superbubble radius in kpc
  real(dp)::eAGN_K =1d0       ! AGN energy efficiency in Radio mode 
  real(dp)::eAGN_T =0.15d0    ! AGN energy efficiency in Quasar mode
  real(dp)::TAGN   =0.0d0     ! AGN temperature floor for energy release in Quasar mode
  real(dp)::T2maxAGN=1d10     ! Maximum temperature allowed after AGN energy input 
  real(dp)::jetfrac=0.0d0     ! Fraction of accreted mass before releasing jet AGN energy
  real(dp)::Mseed  =1.0d5     ! Mass of the initial sink particle in solar mass
  real(dp)::boost_acc=2.0d0   ! Boost power factor for the accretion rate
  real(dp)::boost_drag=2.0d0  ! Boost power factor for drag force
  real(dp)::r_gal  =1.0d2     ! SMBH sphere radius of influence in kpc
  real(dp)::sigmav_max=10d15  ! Maximum relative velocity in the Bondi accretion rate in kpc
  real(dp)::mloadAGN=1d2      ! Mass loading factor of the jet
  real(dp)::f_bondi=1d0       ! Fraction of the Bondi rate
  real(dp)::maxspin=0.998d0   ! Maximum spin amplitude
  real(dp)::X_floor=1.0d-2    ! radio/quasar floor
  real(dp)::rmerge=1.0d0      ! Number of dx_min to allow for BH coalescence
  real(dp)::star_ratio_floor=0.25d0
  real(dp)::d_jeans_thre=0.d0 ! Gas density threshold to trigger Jeans-based refinement criterion

  logical ::self_shielding=.false.
  logical ::pressure_fix=.false.
  logical ::nordlund_fix=.true.
  logical ::cooling=.false.
  logical ::neq_chem=.false.  ! Non-equilbrium chemistry activated
  logical ::isothermal=.false.
  logical ::metal=.false.
  logical ::haardt_madau=.false.
  logical ::delayed_cooling=.false.
  logical ::smbh=.false.
  logical ::agn=.false.
  logical::convert_birth_times=.false. ! Convert stellar birthtimes: conformal-> proper
!jhshin2
  logical ::ir_feedback=.false. ! Activate ir feedback from accreting sinks
!jhshin1
  logical ::sf_virial=.false.   ! Activate SF Virial criterion
  logical ::sf_birth_properties=.true. ! Output birth properties of stars
!jhshin2
  logical ::bondi=.true.      ! Activate Bondi accretion onto sink particle 
  logical ::sink_AGN=.false.  ! Activate AGN feedback on sink particles
  logical ::drag=.false.      ! Activate drag force
  logical ::random_jet=.false.
  logical ::spin_bh=.true.
  logical ::vrel_merge=.true.
  logical ::bhspinmerge=.true.
  logical ::selfgrav=.true.
  logical ::mad_jet=.false.

!chemo flags and variables
  real(dp)::rcell  =2.0D0     ! Supernovae superbubble radius in cells
  real(dp)::tol    =5.D-6     ! Tolerance factor used in the chemo version
  real(dp)::eps_sn2=1.D-1     ! SNII efficiency           
  real(dp)::eps_sn1=1.D-1     ! SNIa efficiency     
  integer ::nelt   =0         ! Number of chemical elements
  character(len=200) :: yieldtablefilename


  ! Output times
  real(dp),dimension(1:MAXOUT)::aout=1.1       ! Output expansion factors
  real(dp),dimension(1:MAXOUT)::tout=0.0       ! Output times

  ! Movie
  integer::imovout=0             ! Increment for output times
  integer::imov=1                ! Initialize
!jhshin1
  real(kind=8)::tstartmov=0.,astartmov=0.
!jhshin2
  real(kind=8)::tendmov=0.,aendmov=0.
  real(kind=8),allocatable,dimension(:)::amovout,tmovout
  logical::movie=.false.
!jhshin1
  logical::zoom_only=.false.
!jhshin2
  integer::nw_frame=512 ! prev: nx_frame, width of frame in pixels
  integer::nh_frame=512 ! prev: ny_frame, height of frame in pixels
  integer::levelmax_frame=0
  integer::ivar_frame=1
  real(kind=8),dimension(1:20)::xcentre_frame=0d0
  real(kind=8),dimension(1:20)::ycentre_frame=0d0
  real(kind=8),dimension(1:20)::zcentre_frame=0d0
  real(kind=8),dimension(1:10)::deltax_frame=0d0
  real(kind=8),dimension(1:10)::deltay_frame=0d0
  real(kind=8),dimension(1:10)::deltaz_frame=0d0
!jhshin1
  real(kind=8),dimension(1:5)::dtheta_camera=0d0
  real(kind=8),dimension(1:5)::dphi_camera=0d0
  real(kind=8),dimension(1:5)::theta_camera=0d0
  real(kind=8),dimension(1:5)::phi_camera=0d0
  real(kind=8),dimension(1:5)::tstart_theta_camera=0d0
  real(kind=8),dimension(1:5)::tstart_phi_camera=0d0
  real(kind=8),dimension(1:5)::tend_theta_camera=0d0
  real(kind=8),dimension(1:5)::tend_phi_camera=0d0
  real(kind=8),dimension(1:5)::focal_camera=0d0
  real(kind=8),dimension(1:5)::smooth_frame=1d0
  logical,dimension(1:5)::perspective_camera=.false.
!jhshin2
  character(LEN=5)::proj_axis='z' ! x->x, y->y, projection along z
!jhshin1
  character(LEN=6),dimension(1:5)::shader_frame='square'
!jhshin2
#ifdef SOLVERmhd
  integer,dimension(0:NVAR+6)::movie_vars=0
  character(len=5),dimension(0:NVAR+6)::movie_vars_txt=''
#else
  integer,dimension(0:NVAR+2)::movie_vars=0
  character(len=5),dimension(0:NVAR+2)::movie_vars_txt=''
#endif

  ! Refinement parameters for each level
  real(dp),dimension(1:MAXLEVEL)::m_refine =-1.0 ! Lagrangian threshold
  real(dp),dimension(1:MAXLEVEL)::r_refine =-1.0 ! Radius of refinement region
  real(dp),dimension(1:MAXLEVEL)::x_refine = 0.0 ! Center of refinement region
  real(dp),dimension(1:MAXLEVEL)::y_refine = 0.0 ! Center of refinement region
  real(dp),dimension(1:MAXLEVEL)::z_refine = 0.0 ! Center of refinement region
  real(dp),dimension(1:MAXLEVEL)::exp_refine = 2.0 ! Exponent for distance
  real(dp),dimension(1:MAXLEVEL)::a_refine = 1.0 ! Ellipticity (Y/X)
  real(dp),dimension(1:MAXLEVEL)::b_refine = 1.0 ! Ellipticity (Z/X)
  real(dp)::var_cut_refine=-1.0 ! Threshold for variable-based refinement
  real(dp)::mass_cut_refine=-1.0 ! Mass threshold for particle-based refinement
  integer::ivar_refine=-1 ! Variable index for refinement
  logical::sink_refine=.false. ! Fully refine on sink particles

  ! Initial condition files for each level
  logical::multiple=.false.
  character(LEN=80),dimension(1:MAXLEVEL)::initfile=' '
  character(LEN=20)::filetype='ascii'

  ! Initial condition regions parameters
  integer,parameter::MAXREGION=100
  integer                           ::nregion=0
  character(LEN=10),dimension(1:MAXREGION)::region_type='square'
  real(dp),dimension(1:MAXREGION)   ::x_center=0.
  real(dp),dimension(1:MAXREGION)   ::y_center=0.
  real(dp),dimension(1:MAXREGION)   ::z_center=0.
  real(dp),dimension(1:MAXREGION)   ::length_x=1.E10
  real(dp),dimension(1:MAXREGION)   ::length_y=1.E10
  real(dp),dimension(1:MAXREGION)   ::length_z=1.E10
  real(dp),dimension(1:MAXREGION)   ::exp_region=2.0

  ! Boundary conditions parameters
  integer,parameter::MAXBOUND=100
  logical                           ::simple_boundary=.false.
  integer                           ::nboundary=0
  integer                           ::icoarse_min=0
  integer                           ::icoarse_max=0
  integer                           ::jcoarse_min=0
  integer                           ::jcoarse_max=0
  integer                           ::kcoarse_min=0
  integer                           ::kcoarse_max=0
  integer ,dimension(1:MAXBOUND)    ::boundary_type=0
  integer ,dimension(1:MAXBOUND)    ::ibound_min=0
  integer ,dimension(1:MAXBOUND)    ::ibound_max=0
  integer ,dimension(1:MAXBOUND)    ::jbound_min=0
  integer ,dimension(1:MAXBOUND)    ::jbound_max=0
  integer ,dimension(1:MAXBOUND)    ::kbound_min=0
  integer ,dimension(1:MAXBOUND)    ::kbound_max=0
  logical                           ::no_inflow=.false.

  !Number of processes sharing one token
  !Only one process can write at a time in an I/O group
  integer::IOGROUPSIZE=0           ! Main snapshot
  integer::IOGROUPSIZECONE=0       ! Lightcone
  integer::IOGROUPSIZEREP=0        ! Subfolder size
  logical::withoutmkdir=.false.    !If true mkdir should be done before the run
  logical::print_when_io=.false.   !If true print when IO
  logical::synchro_when_io=.false. !If true synchronize when IO

  ! Stochastic feedback:
  real(dp)::SN_dT2_min = 0d0        ! Temperature increase for stochastic SNe
  ! Velocity kick maximum for new stellar particles (in random direction)
  real(dp)::SF_kick_kms = -1d0
  ! Keep track of densities at which stellar particles are born and go SN:
  logical::write_stellar_densities=.false.
  ! Kimm feedback stuff:
  ! Efficiency of stellar feedback energy used to heat up/blow out the gas:
!  real(dp)::eff_sfbk=1.0D0
  real(dp)::E_SNII=1d51        ! different from ESN used in feedback.f90 
  integer ::loading_type = 0 ! 0: uniform, 1: turbulence-based
  ! Realistic time delay for individual star particle. t_delay is oldest age to consider if set:
  logical ::sn2_real_delay=.false.
  logical ::use_initial_mass=.false. ! read/write initial mass of particles
  real(dp) :: A_SN=3d5
  real(dp) :: expN_SN=-2d0/17d0
  real(dp) :: A_SN_geen = 5d5 
  real(dp) :: expN_SN_boost = -0.15    ! parameter needed to make sure 
  real(dp) :: expE_SN_boost = 16./17.  ! that pSN_geen = A_SN_Geen during an adiabatic phase
  real(dp) :: porosity = 1.0d0

  ! Star formation stuff:
  character(len=10)::star_maker='density' ! density,hopkins, cen, padoan
  real(dp)::T2thres_SF=1d10  ! Temperature threshold
  real(dp)::fstar_min=1d0     ! Mstar,min = nH*dx_min^3*fstar_min
  real(dp)::M_SNII=10d0       ! Mean progenitor mass of the TypeII SNe
  real(dp)::sf_lamjt=1d0      ! Number of cells below which the Jeans length trigger SF
  character(len=8 )::star_imf=''          ! salpeter,kroupa,(chabrier)
  ! Density above which SF will occur regardless of the kinematic condition (dark cloud):
  real(dp)::n_dc=1d10
  ! Density above which Hopkins and Cen SF routines are evaluated:
  real(dp)::n_gmc=1D-1
  ! Star particle mass in units of the number of SN:
  integer ::nsn2mass=-1

end module amr_parameters
