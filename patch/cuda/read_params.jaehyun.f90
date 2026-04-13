subroutine read_params
  use amr_commons
  use pm_parameters
  use poisson_parameters
  use hydro_parameters
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  !--------------------------------------------------
  ! Local variables
  !--------------------------------------------------
  integer::i,narg,iargc,ierr,levelmax
  character(LEN=80)::infile
  character(LEN=80)::cmdarg
  integer(kind=8)::ngridtot=0
  integer(kind=8)::nparttot=0
  real(kind=8)::delta_tout=0,tend=0
  real(kind=8)::delta_aout=0,aend=0
  logical::nml_ok
  integer,parameter::tag=1134
  integer::dummy_io,info2
  !--------------------------------------------------
  ! Namelist definitions
  !--------------------------------------------------
!jhshin1
  namelist/run_params/clumpfind,cosmo,pic,sink,sinkprops,lightcone,poisson,hydro,rt,verbose,debug &
       & ,nrestart,ncontrol,nstepmax,nsubcycle,nremap,remap_thresh,ordering &
       & ,bisec_tol,static,geom,overload,cost_weighting,aton,varcpu_chunk_nfile &
       & ,memory_balance,mem_weight_grid,mem_weight_part,mem_weight_sink &
       & ,time_balance_alpha &
       & ,jobcontrolfile &
       & ,gpu_hydro,gpu_poisson,gpu_fft,gpu_sink,gpu_auto_tune,n_cuda_streams &
       & ,use_fftw &
       & ,dump_pk &
       & ,exchange_method &
       & ,use_neutrino &
       & ,sidm &
       & ,de_perturb &
       & ,use_mond &
       & ,use_fR &
       & ,use_nDGP &
       & ,use_symmetron &
       & ,use_dilaton &
       & ,use_galileon &
       & ,use_coupled_de &
       & ,use_ede &
       & ,use_sgs &
       & ,use_adm
  ! Non-standard model namelists (read only when enabled)
  namelist/cpl_params/w0,wa,cs2_de,de_table
  namelist/neutrino_params/omega_nu,neutrino_table
  namelist/fR_params/fR0,fR_n,n_iter_fR,fR_eps
  namelist/nDGP_params/omega_rc,nDGP_branch,n_iter_nDGP,nDGP_eps
  namelist/symmetron_params/a_ssb,beta_symmetron,L_symmetron, &
       & n_iter_symmetron,symmetron_eps
  namelist/dilaton_params/beta_dilaton,L_dilaton,a0_dilaton, &
       & n_iter_dilaton,dilaton_eps
  namelist/galileon_params/c2_galileon,c3_galileon, &
       & n_iter_galileon,galileon_eps
  namelist/coupled_de_params/beta_cde
  namelist/ede_params/omega_ede,z_ede,w_ede
  namelist/sgs_params/sgs_C_prod,sgs_C_diss,sgs_C_smag,sgs_floor,sgs_cap,sgs_e_init,sgs_hydro
  namelist/sidm_params/sidm_cross_section,sidm_npart_min, &
       & sidm_type,sidm_v0,sidm_power, &
       & sidm_courant, &
       & sidm_angular,sidm_epsilon, &
       & sidm_inelastic,sidm_delta,sidm_frac_excited, &
       & sidm_nstates,sidm_energy,sidm_frac_init, &
       & sidm_a_type,sidm_a_transition,sidm_sigma_ratio,sidm_a_width, &
       & sidm_fdiss, &
       & sidm_baryon,sidm_baryon_sigma,sidm_baryon_power
namelist/adm_params/adm_alpha,adm_mp,adm_me_ratio,adm_xi, &
       & adm_cross_section
  namelist/fdm_params/m_axion,fdm_courant,fdm_nrefine_dB,fdm_hybrid
  namelist/mond_params/a0_mond,mond_mu_type,mond_type, &
       & n_iter_mond,mond_eps,g_ext_mond
  namelist/cosmo_params/omega_b,omega_m,omega_l,h0
  namelist/output_params/noutput,foutput,fbackup,aout,tout,output_mode &
       & ,tend,delta_tout,aend,delta_aout,gadget_output,walltime_hrs,minutes_dump &
       & ,informat,outformat
  namelist/amr_params/levelmin,levelmax,ngridmax,ngridtot &
       & ,npartmax,nparttot,nexpand,boxlen,nsinkmax
  namelist/poisson_params/epsilon,gravity_type,gravity_params &
       & ,cg_levelmin,cic_levelmax
  namelist/lightcone_params/zmax_cone  &
       & ,elongated_axis_cone1,observer_cone1,minboxr_cone1,maxboxr_cone1 &
       & ,elongated_axis_cone2,observer_cone2,minboxr_cone2,maxboxr_cone2
  namelist/spherical_region_params/spherical_region &
       & ,scenter1,scenter2,scenter3,scenter4,scenter5 &
       & ,sradius1,sradius2,sradius3,sradius4,sradius5
!yonghwi (changed all from previous namelist/movie_params/)
  namelist/movie_params/levelmax_frame,nw_frame,nh_frame,ivar_frame &
       & ,xcentre_frame,ycentre_frame,zcentre_frame &
       & ,deltax_frame,deltay_frame,deltaz_frame,movie,zoom_only &
       & ,imovout,imov,tstartmov,astartmov,tendmov,aendmov &
       & ,proj_axis,movie_vars,movie_vars_txt &
       & ,theta_camera,phi_camera,dtheta_camera,dphi_camera,focal_camera &
       & ,perspective_camera,smooth_frame,shader_frame &
       & ,tstart_theta_camera,tstart_phi_camera &
       & ,tend_theta_camera,tend_phi_camera,dist_camera,ddist_camera
!yonghwi

!jhshin2
  ! MPI initialization
#ifndef WITHOUTMPI
  call MPI_INIT(ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,ncpu,ierr)
  myid=myid+1 ! Careful with this...
#endif
#ifdef WITHOUTMPI
  ncpu=1
  myid=1
#endif
  !--------------------------------------------------
  ! Advertise cuRAMSES
  !--------------------------------------------------
  if(myid==1)then
  write(*,*)'  _/_/_/   _/    _/   _/_/_/         _/_/     _/    _/    _/_/_/    _/_/_/_/    _/_/_/ '
  write(*,*)' _/        _/    _/    _/    _/     _/  _/    _/_/_/_/   _/    _/   _/         _/    _/'
  write(*,*)'_/         _/    _/    _/    _/    _/    _/   _/ _/ _/   _/         _/         _/      '
  write(*,*)'_/         _/    _/    _/_/_/     _/_/_/_/    _/    _/     _/_/     _/_/_/       _/_/  '
  write(*,*)'_/         _/    _/    _/    _/   _/    _/    _/    _/         _/   _/               _/'
  write(*,*)' _/        _/    _/    _/    _/   _/    _/    _/    _/   _/    _/   _/         _/    _/'
  write(*,*)'  _/_/_/    _/_/_/     _/    _/   _/    _/    _/    _/    _/_/_/    _/_/_/_/    _/_/_/ '
  write(*,*)'                              Version 3.0                                             '
  write(*,*)'             written by Romain Teyssier (University of Zurich)                        '
  write(*,*)'             GPU acceleration by Juhan Kim (KIAS)                                       '
  write(*,*)'                     (c) CEA 1999-2007, UZH 2008-2014                                 '
  write(*,*)'        GPU & optimization by Juhan Kim (KIAS) 2026                                   '
  write(*,*)' '
  write(*,'(" Working with nproc = ",I4," for ndim = ",I1)')ncpu,ndim
  ! Check nvar is not too small
#ifdef SOLVERhydro
  write(*,'(" Using solver = hydro with nvar = ",I2)')nvar
  if(nvar<ndim+2)then
     write(*,*)'You should have: nvar>=ndim+2'
     write(*,'(" Please recompile with -DNVAR=",I2)')ndim+2
     call clean_stop
  endif
#endif
#ifdef SOLVERmhd
  write(*,'(" Using solver = mhd with nvar = ",I2)')nvar
  if(nvar<8)then
     write(*,*)'You should have: nvar>=8'
     write(*,'(" Please recompile with -DNVAR=8")')
     call clean_stop
  endif
#endif
  
  !Write I/O group size information
  if(IOGROUPSIZE>0.or.IOGROUPSIZECONE>0.or.IOGROUPSIZEREP>0)write(*,*)' '
  if(IOGROUPSIZE>0) write(*,*)'IOGROUPSIZE=',IOGROUPSIZE
  if(IOGROUPSIZECONE>0) write(*,*)'IOGROUPSIZECONE=',IOGROUPSIZECONE
  if(IOGROUPSIZEREP>0) write(*,*)'IOGROUPSIZEREP=',IOGROUPSIZEREP
  if(IOGROUPSIZE>0.or.IOGROUPSIZECONE>0.or.IOGROUPSIZEREP>0)write(*,*)' '

  ! Write information about git version
  call write_gitinfo

  ! Read namelist filename from command line argument
  narg = iargc()
  IF(narg .LT. 1)THEN
     write(*,*)'You should type: ramses3d input.nml [nrestart]'
     write(*,*)'File input.nml should contain a parameter namelist'
     write(*,*)'nrestart is optional'
     call clean_stop
  END IF
  CALL getarg(1,infile)
  endif
#ifndef WITHOUTMPI
  call MPI_BCAST(infile,80,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
#endif

  !-------------------------------------------------
  ! Read the namelist
  !-------------------------------------------------

  ! Wait for the token                                                                                                                                                                                
#ifndef WITHOUTMPI
     if(IOGROUPSIZE>0) then
        if (mod(myid-1,IOGROUPSIZE)/=0) then
           call MPI_RECV(dummy_io,1,MPI_INTEGER,myid-1-1,tag,&
                & MPI_COMM_WORLD,MPI_STATUS_IGNORE,info2)
        end if
     endif
#endif


  namelist_file=TRIM(infile)
  INQUIRE(file=infile,exist=nml_ok)
  if(.not. nml_ok)then
     if(myid==1)then
        write(*,*)'File '//TRIM(infile)//' does not exist'
     endif
     call clean_stop
  end if

  open(1,file=infile)
  rewind(1)
  read(1,NML=run_params)
  rewind(1)
  read(1,NML=output_params)
  rewind(1)
  read(1,NML=amr_params)
  rewind(1)
  read(1,NML=lightcone_params,END=84)
84 continue
  rewind(1)
  read(1,NML=spherical_region_params,END=83)
83 continue
  rewind(1)
  read(1,NML=movie_params,END=82)
82 continue
  rewind(1)
  read(1,NML=poisson_params,END=81)
81 continue
  rewind(1)
  read(1,NML=cosmo_params,END=80)
80 continue
  rewind(1)
  read(1,NML=cpl_params,END=79)
79 continue
  rewind(1)
  read(1,NML=neutrino_params,END=78)
78 continue
  rewind(1)
  read(1,NML=sidm_params,END=77)
77 continue
  ! aDM parameters
  if(use_adm) then
     rewind(1)
     read(1,NML=adm_params,END=76)
76   continue
  end if
  ! FDM parameters
  if(use_fdm) then
     rewind(1)
     read(1,NML=fdm_params,END=176)
176  continue
  end if
  rewind(1)
  read(1,NML=mond_params,END=75)
75 continue
  ! f(R) parameters
  if(use_fR) then
     rewind(1)
     read(1,NML=fR_params,END=74)
74   continue
  end if
  ! nDGP parameters
  if(use_nDGP) then
     rewind(1)
     read(1,NML=nDGP_params,END=73)
73   continue
  end if
  ! Symmetron parameters
  if(use_symmetron) then
     rewind(1)
     read(1,NML=symmetron_params,END=72)
72   continue
  end if
  ! Dilaton parameters
  if(use_dilaton) then
     rewind(1)
     read(1,NML=dilaton_params,END=71)
71   continue
  end if
  ! Galileon parameters
  if(use_galileon) then
     rewind(1)
     read(1,NML=galileon_params,END=70)
70   continue
  end if
  ! Coupled DE parameters
  if(use_coupled_de) then
     rewind(1)
     read(1,NML=coupled_de_params,END=69)
69   continue
  end if
  ! Early DE parameters
  if(use_ede) then
     rewind(1)
     read(1,NML=ede_params,END=68)
68   continue
  end if
  ! SGS turbulence parameters
  if(use_sgs) then
     rewind(1)
     read(1,NML=sgs_params,END=67)
67   continue
  end if

  !-------------------------------------------------
  ! Read optional nrestart command-line argument
  !-------------------------------------------------
  if (myid==1 .and. narg == 2) then
     CALL getarg(2,cmdarg)
     read(cmdarg,*) nrestart
  endif

#ifndef WITHOUTMPI
  call MPI_BCAST(nrestart,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
#endif

  !-------------------------------------------------
  ! GPU acceleration: disable if not compiled with USE_CUDA
  !-------------------------------------------------
#ifndef HYDRO_CUDA
  if(gpu_hydro .or. gpu_poisson .or. gpu_fft .or. gpu_sink) then
     if(myid==1) write(*,*) 'WARNING: gpu_* options ignored (not compiled with USE_CUDA)'
     gpu_hydro = .false.
     gpu_poisson = .false.
     gpu_fft = .false.
     gpu_sink = .false.
  end if
#else
  if(myid==1 .and. (gpu_hydro .or. gpu_poisson .or. gpu_fft .or. gpu_sink)) then
     write(*,'(A,L1,A,L1,A,L1,A,L1,A,I0)') &
          ' GPU acceleration: hydro=',gpu_hydro, &
          ' poisson=',gpu_poisson,' fft=',gpu_fft,' sink=',gpu_sink, &
          ' streams=',n_cuda_streams
  end if
#endif

  !-------------------------------------------------
  ! FFTW3 CPU Poisson solver: disable if not compiled with USE_FFTW
  !-------------------------------------------------
#ifndef USE_FFTW
  if(use_fftw) then
     if(myid==1) write(*,*) 'WARNING: use_fftw ignored (not compiled with USE_FFTW)'
     use_fftw = .false.
  end if
#else
  if(myid==1 .and. use_fftw) then
     write(*,'(A)') ' FFTW3 CPU direct Poisson solver enabled (use_fftw=T)'
  end if
#endif

  !-------------------------------------------------
  ! Auto-compute mem_weight_grid from nvar if sentinel (0)
  ! Per cell (×twotondim per grid):
  !   Hydro:   2*nvar*8 (uold+unew)
  !   Topo:    5*4 (son,flag1,flag2,cpu_map,cpu_map2) + 8 (hilbert_key)
  !   PFix:    2*8 (enew,divu)  [pressure_fix, default for cosmo]
  !   Poisson: 7*8 (rho,rho_star,phi,phi_old,f*3)
  ! Per grid:
  !   AMR:  3*8+3*4+6*4 (xg,father/next/prev,nbor) = 48
  !   Part: 3*4 (headp,tailp,numbp) = 12
  !-------------------------------------------------
  if(memory_balance .and. mem_weight_grid <= 0) then
     mem_weight_grid = twotondim * (2*nvar*8 + 28 + 16 + 56) + 48 + 12
     if(myid==1) write(*,'(A,I6,A,I3,A)') &
          ' Memory balance: mem_weight_grid=',mem_weight_grid,' (nvar=',nvar,')'
  end if

  !-------------------------------------------------
  ! Exchange method auto-tune
  !-------------------------------------------------
  if(ordering=='ksection' .and. myid==1) then
     write(*,'(A,A)') ' Exchange method: ', trim(exchange_method)
  end if

  !-------------------------------------------------
  ! DE perturbation (CPL)
  ! Three modes:
  !   1) de_table provided → table-based linear response (any cs2_de)
  !   2) no de_table, cs2_de>0 → kappa2/alpha quasi-static method
  !   3) no de_table, cs2_de<=0 → unsupported, disable
  !-------------------------------------------------
  if(de_perturb) then
     if(.not. cosmo) then
        if(myid==1) write(*,*) 'WARNING: de_perturb=T but not cosmo run, disabling'
        de_perturb = .false.
     else if(len_trim(de_table) > 0) then
        ! Table-based linear response (works for any cs2_de)
        if(myid==1) then
           write(*,'(A,A)') ' DE perturbation (table): ', trim(de_table)
           write(*,'(A,ES10.3,A,F6.3,A,F6.3)') &
                '   cs2_de=', cs2_de, ' w0=', w0, ' wa=', wa
        end if
     else if(cs2_de <= 0.0d0) then
        if(myid==1) write(*,*) 'WARNING: de_perturb=T, no de_table, cs2_de<=0 -> disabling'
        de_perturb = .false.
     else
        ! Fallback: kappa2/alpha quasi-static method
        if(myid==1) then
           write(*,'(A,ES10.3,A,F6.3,A,F6.3)') &
                ' DE perturbation (kappa2/alpha): cs2_de=', cs2_de, &
                ' w0=', w0, ' wa=', wa
        end if
     end if
  end if

  !-------------------------------------------------
  ! Neutrino linear response
  !-------------------------------------------------
  if(use_neutrino) then
     if(omega_nu <= 0.0d0) then
        if(myid==1) write(*,*) 'WARNING: use_neutrino=T but omega_nu<=0, disabling'
        use_neutrino = .false.
     else if(len_trim(neutrino_table) == 0) then
        if(myid==1) write(*,*) 'WARNING: use_neutrino=T but neutrino_table not set, disabling'
        use_neutrino = .false.
     else
        if(myid==1) then
           write(*,'(A,F7.4,A,F7.4)') &
                ' Neutrino linear response: omega_nu=', omega_nu, &
                ' omega_cb=', omega_m - omega_nu
           write(*,'(A,A)') '   table: ', trim(neutrino_table)
        end if
     end if
  end if

  !-------------------------------------------------
  ! SIDM (Self-Interacting Dark Matter) scattering
  !-------------------------------------------------
  if(sidm) then
     if(.not. pic) then
        if(myid==1) write(*,*) 'ERROR: sidm=T requires pic=T'
        call clean_stop
     end if
     if(sidm_cross_section <= 0.0d0) then
        if(myid==1) write(*,*) 'ERROR: sidm=T but sidm_cross_section<=0'
        call clean_stop
     end if
     if(sidm_inelastic .and. sidm_delta <= 0.0d0 &
          .and. sidm_energy(1) <= 0.0d0) then
        if(myid==1) write(*,*) 'ERROR: sidm_inelastic=T but no energy splitting set'
        call clean_stop
     end if
     ! Auto-populate multi-state arrays from 2-state shortcut
     if(sidm_inelastic .and. sidm_nstates==2 &
          .and. sidm_energy(1)==0.0d0 .and. sidm_delta>0.0d0) then
        sidm_energy(0) = 0.0d0
        sidm_energy(1) = sidm_delta
     end if
     if(sidm_inelastic .and. sidm_nstates==2 &
          .and. sidm_frac_init(0)==0.0d0 .and. sidm_frac_init(1)==0.0d0 &
          .and. sidm_frac_excited>0.0d0) then
        sidm_frac_init(0) = 1.0d0 - sidm_frac_excited
        sidm_frac_init(1) = sidm_frac_excited
     end if
     if(myid==1) then
        write(*,'(A,ES10.3,A)') ' SIDM enabled: sigma/m=', sidm_cross_section, ' cm^2/g'
        write(*,'(A,A)')        '   cross-section type: ', trim(sidm_type)
        if(trim(sidm_type) /= 'constant') then
           write(*,'(A,F8.1,A)') '   v0=', sidm_v0, ' km/s'
           if(trim(sidm_type) == 'power_law') &
                write(*,'(A,F6.2)') '   power=', sidm_power
        end if
        write(*,'(A,I4)')        '   npart_min=', sidm_npart_min
        write(*,'(A,F5.2)')      '   courant=', sidm_courant
        write(*,'(A,A)')         '   angular: ', trim(sidm_angular)
        if(trim(sidm_angular) == 'rutherford') &
             write(*,'(A,ES10.3)') '   epsilon=', sidm_epsilon
        if(sidm_inelastic) then
           write(*,'(A,I3)')       '   iSIDM: nstates=', sidm_nstates
           write(*,'(A,10ES10.3)') '   energies [keV]=', &
                sidm_energy(0:sidm_nstates-1)
           write(*,'(A,10F7.3)')   '   frac_init=', &
                sidm_frac_init(0:sidm_nstates-1)
        end if
        if(trim(sidm_a_type) /= 'none') then
           write(*,'(A,A)')         '   Phase transition type: ', trim(sidm_a_type)
           write(*,'(A,F8.4)')      '   a_transition=', sidm_a_transition
           write(*,'(A,F8.2)')      '   sigma_ratio=', sidm_sigma_ratio
           if(trim(sidm_a_type) == 'sigmoid') &
                write(*,'(A,F8.4)') '   a_width=', sidm_a_width
        end if
        if(sidm_fdiss > 0.0d0) then
           write(*,'(A,F6.3)')     '   dSIDM: fdiss=', sidm_fdiss
        end if
        if(sidm_baryon) then
           write(*,'(A,ES10.3,A,F5.1)') &
                '   IDM: sigma_DM-b/m=', sidm_baryon_sigma, &
                ' cm^2/g  power=', sidm_baryon_power
        end if
     end if
     if(sidm_fdiss < 0.0d0 .or. sidm_fdiss >= 1.0d0) then
        if(myid==1) write(*,*) 'ERROR: sidm_fdiss must be in [0,1)'
        call clean_stop
     end if
     if(sidm_baryon .and. .not.hydro) then
        if(myid==1) write(*,*) 'ERROR: sidm_baryon requires hydro=T'
        call clean_stop
     end if
     if(sidm_baryon .and. sidm_baryon_sigma <= 0.0d0) then
        if(myid==1) write(*,*) 'ERROR: sidm_baryon_sigma must be > 0'
        call clean_stop
     end if
  end if

  !-------------------------------------------------
  ! Atomic Dark Matter (aDM)
  !-------------------------------------------------
  if(use_adm) then
     if(.not. sidm) then
        if(myid==1) write(*,*) 'ERROR: use_adm=T requires sidm=T'
        call clean_stop
     end if
     if(adm_alpha <= 0.0d0 .or. adm_alpha >= 1.0d0) then
        if(myid==1) write(*,*) 'ERROR: adm_alpha must be in (0,1)'
        call clean_stop
     end if
     if(adm_mp <= 0.0d0) then
        if(myid==1) write(*,*) 'ERROR: adm_mp must be > 0'
        call clean_stop
     end if
     if(adm_me_ratio <= 0.0d0 .or. adm_me_ratio >= 1.0d0) then
        if(myid==1) write(*,*) 'ERROR: adm_me_ratio must be in (0,1)'
        call clean_stop
     end if
     if(adm_xi <= 0.0d0) then
        if(myid==1) write(*,*) 'ERROR: adm_xi must be > 0'
        call clean_stop
     end if
     if(myid==1) then
        write(*,'(A)') ' Atomic Dark Matter (aDM) enabled:'
        write(*,'(A,ES10.3)') '   alpha_D  =', adm_alpha
        write(*,'(A,F8.3,A)') '   m_p''     =', adm_mp, ' GeV'
        write(*,'(A,ES10.3)') '   m_e''/m_p''=', adm_me_ratio
        write(*,'(A,F6.3)')   '   xi       =', adm_xi
        write(*,'(A,ES10.3,A)') '   sigma/m  =', adm_cross_section, ' cm^2/g'
     end if
  end if

  !-------------------------------------------------
  ! Fuzzy Dark Matter (FDM)
  !-------------------------------------------------
  if(use_fdm) then
     if(.not. poisson) then
        if(myid==1) write(*,*) 'ERROR: use_fdm=T requires poisson=T'
        call clean_stop
     end if
     if(m_axion <= 0.0d0) then
        if(myid==1) write(*,*) 'ERROR: m_axion must be > 0'
        call clean_stop
     end if
     if(sidm) then
        if(myid==1) write(*,*) 'ERROR: FDM and SIDM are mutually exclusive'
        call clean_stop
     end if
     if(myid==1) then
        write(*,'(A)')         ' Fuzzy Dark Matter (FDM) enabled:'
        write(*,'(A,ES10.3,A)') '   m_axion  =', m_axion, ' eV'
        write(*,'(A,F5.2)')    '   courant  =', fdm_courant
        write(*,'(A,I3)')      '   nrefine_dB=', fdm_nrefine_dB
        write(*,'(A,L1)')      '   hybrid   =', fdm_hybrid
     end if
  end if

  !-------------------------------------------------
  ! MOND (Modified Newtonian Dynamics)
  !-------------------------------------------------
  if(use_mond) then
     if(.not. poisson) then
        if(myid==1) write(*,*) 'ERROR: use_mond=T requires poisson=T'
        call clean_stop
     end if
     if(a0_mond <= 0.0d0) then
        if(myid==1) write(*,*) 'ERROR: use_mond=T but a0_mond<=0'
        call clean_stop
     end if
     if(mond_mu_type < 1 .or. mond_mu_type > 2) then
        if(myid==1) write(*,*) 'ERROR: mond_mu_type must be 1 (simple) or 2 (standard)'
        call clean_stop
     end if
     if(mond_type < 0 .or. mond_type > 2) then
        if(myid==1) write(*,*) 'ERROR: mond_type must be 0, 1 (QUMOND), or 2 (AQUAL)'
        call clean_stop
     end if
     if(mond_type == 2) then
        if(n_iter_mond < 1) then
           if(myid==1) write(*,*) 'ERROR: n_iter_mond must be >= 1'
           call clean_stop
        end if
        if(mond_eps <= 0d0) then
           if(myid==1) write(*,*) 'ERROR: mond_eps must be > 0'
           call clean_stop
        end if
     end if
     ! Warn if DM particles likely present
     if(pic) then
        if(myid==1) then
           write(*,'(A)') ' WARNING: use_mond=T with pic=T (DM particles likely present)'
           write(*,'(A)') '   MOND replaces DM — running both gives excess gravity.'
           write(*,'(A)') '   Use DM-free ICs or set pic=.false. for pure MOND.'
        end if
     end if
     if(myid==1) then
        write(*,'(A)') ' MOND (QUMOND) enabled:'
        write(*,'(A,ES12.4,A)') '   a0 = ', a0_mond, ' cm/s^2'
        if(mond_mu_type==1) then
           write(*,'(A)') '   mu-function: simple [mu=x/(1+x)]'
        else
           write(*,'(A)') '   mu-function: standard [mu=x/sqrt(1+x^2)]'
        end if
        if(mond_type==0) then
           write(*,'(A)') '   mode: algebraic QUMOND (Phase 0)'
        else if(mond_type==1) then
           write(*,'(A)') '   mode: full QUMOND with phantom density (Phase 1)'
        else
           write(*,'(A,I3,A,ES10.3)') &
                '   mode: AQUAL iterative (Phase 2), max_iter=', n_iter_mond, &
                ' eps=', mond_eps
        end if
        if(g_ext_mond(1)**2+g_ext_mond(2)**2+g_ext_mond(3)**2 > 0d0) then
           write(*,'(A,3ES12.4,A)') '   g_ext = (', g_ext_mond, ') cm/s^2'
        end if
     end if
  end if

  !-------------------------------------------------
  ! f(R) Hu-Sawicki gravity
  !-------------------------------------------------
  if(use_fR) then
     ! Mutual exclusion checks
     if(use_nDGP) then
        if(myid==1) write(*,*) 'ERROR: Cannot use both f(R) and nDGP simultaneously'
        call clean_stop
     end if
     if(use_mond) then
        if(myid==1) write(*,*) 'ERROR: Cannot use both f(R) and MOND simultaneously'
        call clean_stop
     end if
     if(.not. cosmo) then
        if(myid==1) write(*,*) 'ERROR: f(R) gravity requires cosmo=.true.'
        call clean_stop
     end if
     if(.not. poisson) then
        if(myid==1) write(*,*) 'ERROR: use_fR=T requires poisson=T'
        call clean_stop
     end if
     if(fR0 >= 0d0) then
        if(myid==1) write(*,*) 'ERROR: fR0 must be negative'
        call clean_stop
     end if
     if(fR_n < 1) then
        if(myid==1) write(*,*) 'ERROR: fR_n must be >= 1'
        call clean_stop
     end if
     if(myid==1) then
        write(*,'(A)') ' f(R) Hu-Sawicki gravity enabled'
        write(*,'(A,ES10.3,A,I2)') '   fR0=', fR0, '  n=', fR_n
        write(*,'(A,I3,A,ES10.3)') '   max_iter=', n_iter_fR, '  eps=', fR_eps
     end if
  end if

  !-------------------------------------------------
  ! nDGP gravity
  !-------------------------------------------------
  if(use_nDGP) then
     if(use_mond) then
        if(myid==1) write(*,*) 'ERROR: Cannot use both nDGP and MOND simultaneously'
        call clean_stop
     end if
     if(.not. cosmo) then
        if(myid==1) write(*,*) 'ERROR: nDGP gravity requires cosmo=.true.'
        call clean_stop
     end if
     if(.not. poisson) then
        if(myid==1) write(*,*) 'ERROR: use_nDGP=T requires poisson=T'
        call clean_stop
     end if
     if(omega_rc <= 0d0) then
        if(myid==1) write(*,*) 'ERROR: omega_rc must be > 0'
        call clean_stop
     end if
     if(abs(nDGP_branch) /= 1) then
        if(myid==1) write(*,*) 'ERROR: nDGP_branch must be 1 or -1'
        call clean_stop
     end if
     if(myid==1) then
        write(*,'(A)') ' nDGP gravity enabled'
        write(*,'(A,ES10.3,A,I2)') '   omega_rc=', omega_rc, '  branch=', nDGP_branch
        write(*,'(A,I3,A,ES10.3)') '   max_iter=', n_iter_nDGP, '  eps=', nDGP_eps
     end if
  end if

  !-------------------------------------------------
  ! Scalar field mutual exclusion check
  ! fR, nDGP, MOND, symmetron, dilaton, galileon: at most 1
  !-------------------------------------------------
  i = 0
  if(use_fR)        i = i + 1
  if(use_nDGP)      i = i + 1
  if(use_mond)      i = i + 1
  if(use_symmetron)  i = i + 1
  if(use_dilaton)    i = i + 1
  if(use_galileon)   i = i + 1
  if(i > 1) then
     if(myid==1) write(*,*) 'ERROR: Only one scalar field gravity model allowed at a time'
     if(myid==1) write(*,*) '  (fR, nDGP, MOND, symmetron, dilaton, galileon)'
     call clean_stop
  end if

  !-------------------------------------------------
  ! Symmetron gravity
  !-------------------------------------------------
  if(use_symmetron) then
     if(.not. cosmo) then
        if(myid==1) write(*,*) 'ERROR: Symmetron requires cosmo=.true.'
        call clean_stop
     end if
     if(.not. poisson) then
        if(myid==1) write(*,*) 'ERROR: use_symmetron=T requires poisson=T'
        call clean_stop
     end if
     if(a_ssb <= 0d0 .or. a_ssb >= 1d0) then
        if(myid==1) write(*,*) 'ERROR: a_ssb must be in (0,1)'
        call clean_stop
     end if
     if(L_symmetron <= 0d0) then
        if(myid==1) write(*,*) 'ERROR: L_symmetron must be > 0'
        call clean_stop
     end if
     if(myid==1) then
        write(*,'(A)') ' Symmetron gravity enabled'
        write(*,'(A,F6.3,A,F6.3,A,F8.3,A)') &
             '   a_ssb=', a_ssb, '  beta=', beta_symmetron, '  L=', L_symmetron, ' Mpc/h'
        write(*,'(A,I3,A,ES10.3)') '   max_iter=', n_iter_symmetron, '  eps=', symmetron_eps
     end if
  end if

  !-------------------------------------------------
  ! Dilaton gravity
  !-------------------------------------------------
  if(use_dilaton) then
     if(.not. cosmo) then
        if(myid==1) write(*,*) 'ERROR: Dilaton requires cosmo=.true.'
        call clean_stop
     end if
     if(.not. poisson) then
        if(myid==1) write(*,*) 'ERROR: use_dilaton=T requires poisson=T'
        call clean_stop
     end if
     if(L_dilaton <= 0d0) then
        if(myid==1) write(*,*) 'ERROR: L_dilaton must be > 0'
        call clean_stop
     end if
     if(myid==1) then
        write(*,'(A)') ' Dilaton gravity enabled'
        write(*,'(A,F6.3,A,F6.3,A,F8.3,A)') &
             '   a0=', a0_dilaton, '  beta=', beta_dilaton, '  L=', L_dilaton, ' Mpc/h'
        write(*,'(A,I3,A,ES10.3)') '   max_iter=', n_iter_dilaton, '  eps=', dilaton_eps
     end if
  end if

  !-------------------------------------------------
  ! Galileon (cubic) gravity
  !-------------------------------------------------
  if(use_galileon) then
     if(.not. cosmo) then
        if(myid==1) write(*,*) 'ERROR: Galileon requires cosmo=.true.'
        call clean_stop
     end if
     if(.not. poisson) then
        if(myid==1) write(*,*) 'ERROR: use_galileon=T requires poisson=T'
        call clean_stop
     end if
     if(abs(c3_galileon) < 1d-30) then
        if(myid==1) write(*,*) 'ERROR: c3_galileon must be non-zero'
        call clean_stop
     end if
     if(myid==1) then
        write(*,'(A)') ' Cubic Galileon gravity enabled'
        write(*,'(A,ES10.3,A,ES10.3)') '   c2=', c2_galileon, '  c3=', c3_galileon
        write(*,'(A,I3,A,ES10.3)') '   max_iter=', n_iter_galileon, '  eps=', galileon_eps
     end if
  end if

  !-------------------------------------------------
  ! Coupled Dark Energy
  !-------------------------------------------------
  if(use_coupled_de) then
     if(.not. poisson) then
        if(myid==1) write(*,*) 'ERROR: use_coupled_de=T requires poisson=T'
        call clean_stop
     end if
     if(myid==1) then
        write(*,'(A,F8.4)') ' Coupled Dark Energy enabled: beta_cde=', beta_cde
        write(*,'(A,F10.6)') '   G_eff/G = ', 1d0 + 2d0*beta_cde**2
     end if
  end if

  !-------------------------------------------------
  ! Early Dark Energy (Doran-Robbers)
  !-------------------------------------------------
  if(use_ede) then
     if(.not. cosmo) then
        if(myid==1) write(*,*) 'ERROR: EDE requires cosmo=.true.'
        call clean_stop
     end if
     if(omega_ede < 0d0 .or. omega_ede >= 1d0) then
        if(myid==1) write(*,*) 'ERROR: omega_ede must be in [0,1)'
        call clean_stop
     end if
     if(z_ede <= 0d0) then
        if(myid==1) write(*,*) 'ERROR: z_ede must be > 0'
        call clean_stop
     end if
     if(myid==1) then
        write(*,'(A)') ' Early Dark Energy (Doran-Robbers) enabled'
        write(*,'(A,F8.4,A,F10.1,A,F6.3)') &
             '   omega_ede=', omega_ede, '  z_ede=', z_ede, '  w_ede=', w_ede
     end if
  end if

  !-------------------------------------------------
  ! SGS (Sub-Grid Scale) Turbulence model
  !-------------------------------------------------
  if(use_sgs) then
     if(.not. hydro) then
        if(myid==1) write(*,*) 'ERROR: use_sgs=T requires hydro=T'
        call clean_stop
     end if
     if(sgs_C_diss <= 0d0) then
        if(myid==1) write(*,*) 'ERROR: sgs_C_diss must be > 0'
        call clean_stop
     end if
     if(myid==1) then
        write(*,'(A)') ' SGS turbulence model enabled'
        write(*,'(A,F6.3,A,F6.3,A,F6.3)') &
             '   C_prod=', sgs_C_prod, '  C_diss=', sgs_C_diss, '  C_smag=', sgs_C_smag
        write(*,'(A,ES10.3,A,F6.2,A,ES10.3,A,L1)') '   floor=', sgs_floor, &
             '  cap=', sgs_cap, '  e_init=', sgs_e_init, '  hydro=', sgs_hydro
        write(*,'(A,I3,A,I3)') '   isgs=', isgs, '  nvar=', nvar
     end if
  end if

  !-------------------------------------------------
  ! Compute time step for outputs
  !-------------------------------------------------
  if(tend>0)then
     if(delta_tout==0)delta_tout=tend
     noutput=MIN(int(tend/delta_tout),MAXOUT)
     do i=1,noutput
        tout(i)=dble(i)*delta_tout
     end do
  else if(aend>0)then
     if(delta_aout==0)delta_aout=aend
     noutput=MIN(int(aend/delta_aout),MAXOUT)
     do i=1,noutput
        aout(i)=dble(i)*delta_aout
     end do
  endif
  noutput=MIN(noutput,MAXOUT)
  if(imovout>0) then
     allocate(tmovout(1:imovout))
     allocate(amovout(1:imovout))
     tmovout=1d100
     amovout=1d100
     if(tendmov>0)then
        do i=1,imovout
           tmovout(i)=tendmov*dble(i)/dble(imovout)
        enddo
     endif
     if(aendmov>0)then
        do i=1,imovout
           amovout(i)=aendmov*dble(i)/dble(imovout)
        enddo
     endif
     if(tendmov==0.and.aendmov==0)movie=.false.
  endif
  !--------------------------------------------------
  ! Check for errors in the namelist so far
  !--------------------------------------------------
  levelmin=MAX(levelmin,1)
  nlevelmax=levelmax
  nml_ok=.true.
  if(levelmin<1)then
     if(myid==1)write(*,*)'Error in the namelist:'
     if(myid==1)write(*,*)'levelmin should not be lower than 1 !!!'
     nml_ok=.false.
  end if
  if(nlevelmax<levelmin)then
     if(myid==1)write(*,*)'Error in the namelist:'
     if(myid==1)write(*,*)'levelmax should not be lower than levelmin'
     nml_ok=.false.
  end if
  if(ngridmax==0)then
     if(ngridtot==0)then
        if(myid==1)write(*,*)'Error in the namelist:'
        if(myid==1)write(*,*)'Allocate some space for refinements !!!'
        nml_ok=.false.
     else
        ngridmax=ngridtot/int(ncpu,kind=8)
     endif
  end if
  if(npartmax==0)then
     npartmax=nparttot/int(ncpu,kind=8)
  endif
  if(myid>1)verbose=.false.
  if(sink.and.(.not.pic))then
     pic=.true.
  endif
  if(clumpfind.and.(.not.pic))then
     pic=.true.
  endif
  !if(pic.and.(.not.poisson))then
  !   poisson=.true.
  !endif

  call read_hydro_params(nml_ok)
#ifdef RT
  call rt_read_hydro_params(nml_ok)
#endif
  if (clumpfind)call read_clumpfind_params
  if (movie)call set_movie_vars


  close(1)

  !-------------------------------------------------
  ! FPR (Fixed Proper Resolution, Gnedin 2016)
  ! Must be after read_hydro_params (reads dr_proper)
  !-------------------------------------------------
  if(dr_proper > 0.0d0) then
     if(.not. cosmo) then
        if(myid==1) write(*,*) 'WARNING: dr_proper>0 but not cosmo run, disabling FPR'
        dr_proper = 0.0d0
     else
        if(myid==1) then
           write(*,'(A)') ' FPR (Fixed Proper Resolution, Gnedin 2016) enabled'
           write(*,'(A,F8.3,A)') '   dr_proper = ', dr_proper, ' kpc'
           if(q_refine_holdback) then
              write(*,'(A)') '   mode: FPR + binary holdback'
           else
              write(*,'(A)') '   mode: FPR only (no holdback)'
           end if
        end if
     end if
  end if

  ! Send the token
#ifndef WITHOUTMPI
  if(IOGROUPSIZE>0) then
     if(mod(myid,IOGROUPSIZE)/=0 .and.(myid.lt.ncpu))then
        dummy_io=1
        call MPI_SEND(dummy_io,1,MPI_INTEGER,myid-1+1,tag, &
             & MPI_COMM_WORLD,info2)
     end if
  endif
#endif
  


  !-----------------
  ! Max size checks
  !-----------------
  if(nlevelmax>MAXLEVEL)then
     write(*,*) 'Error: nlevelmax>MAXLEVEL'
     call clean_stop
  end if
  if(nregion>MAXREGION)then
     write(*,*) 'Error: nregion>MAXREGION'
     call clean_stop
  end if
  
  !-----------------------------------
  ! Rearrange level dependent arrays
  !-----------------------------------
  do i=nlevelmax,levelmin,-1
     nexpand   (i)=nexpand   (i-levelmin+1)
     nsubcycle (i)=nsubcycle (i-levelmin+1)
     r_refine  (i)=r_refine  (i-levelmin+1)
     a_refine  (i)=a_refine  (i-levelmin+1)
     b_refine  (i)=b_refine  (i-levelmin+1)
     x_refine  (i)=x_refine  (i-levelmin+1)
     y_refine  (i)=y_refine  (i-levelmin+1)
     z_refine  (i)=z_refine  (i-levelmin+1)
     m_refine  (i)=m_refine  (i-levelmin+1)
     m_basic_refine(i) = m_refine  (i) !(ONS)
     exp_refine(i)=exp_refine(i-levelmin+1)
     initfile  (i)=initfile  (i-levelmin+1)
  end do
  do i=1,levelmin-1
     nexpand   (i)= 1
     nsubcycle (i)= 1
     r_refine  (i)=-1.0
     a_refine  (i)= 1.0
     b_refine  (i)= 1.0
     x_refine  (i)= 0.0
     y_refine  (i)= 0.0
     z_refine  (i)= 0.0
     m_refine  (i)=-1.0
     exp_refine(i)= 2.0
     initfile  (i)= ' '
  end do
  ! Initialize m_refine_eff from m_refine (FPR adjusts at runtime)
  m_refine_eff = m_refine

  if(.not. nml_ok)then
     if(myid==1)write(*,*)'Too many errors in the namelist'
     if(myid==1)write(*,*)'Aborting...'
     call clean_stop
  end if

#ifndef WITHOUTMPI
  call MPI_BARRIER(MPI_COMM_WORLD,ierr)
#endif

end subroutine read_params

