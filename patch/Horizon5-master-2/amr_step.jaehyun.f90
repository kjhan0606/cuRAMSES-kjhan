recursive subroutine amr_step(ilevel,icount)
  use amr_commons
  use pm_commons
  use hydro_commons
  use poisson_commons
#ifdef RT
  use rt_hydro_commons
  use SED_module
  use UV_module
  use coolrates_module, only: update_coolrates_tables
  use rt_cooling_module, only: update_UVrates
#endif
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer::ilevel,icount
  !-------------------------------------------------------------------!
  ! This routine is the adaptive-mesh/adaptive-time-step main driver. !
  ! Each routine is called using a specific order, don't change it,   !
  ! unless you check all consequences first                           !
  !-------------------------------------------------------------------!
  integer::i,idim,ivar,mpi_err
  logical::ok_defrag,output_now_all
!!$  integer::i,idim,ivar
!!$  logical::ok_defrag
  logical,save::first_step=.true.
  integer:: info

  real(kind=4):: real_mem, real_mem_tot

  ! Particle sub-timers
  integer(kind=8) :: pt_t1, pt_t2, pt_rate
  real(dp), save :: pt_mktree=0, pt_killtree=0, pt_synchro=0, pt_move=0, pt_merge=0

  ! Sink sub-timers
  integer(kind=8) :: sk_t1, sk_t2
  real(dp), save :: sk_agn_fb=0, sk_create_sink=0, sk_grow=0, sk_bondi_hoyle=0

#ifdef HYDRO_CUDA
  ! GPU auto-tuning framework
  ! Phase 0: first call → force CPU, record time
  ! Phase 1: second call → force GPU, record time
  ! Phase 2+: use faster path, keep booking times
  !
  ! AGN feedback (sink particles)
  integer, save :: sk_auto_phase = 0
  real(dp), save :: sk_cpu_ref = 0d0, sk_gpu_ref = 0d0
  logical, save :: sk_use_gpu = .false.
  logical, save :: sk_auto_init = .false.
  real(dp) :: sk_dt_agn
  !
  ! Hydro Godunov
  integer, save :: hy_auto_phase = 0
  real(dp), save :: hy_cpu_ref = 0d0, hy_gpu_ref = 0d0
  logical, save :: hy_use_gpu = .false.
  logical, save :: hy_auto_init = .false.
  integer(kind=8) :: hy_t1, hy_t2
  real(dp) :: hy_dt
  !
  ! Poisson MG (gpu_poisson only; gpu_fft excluded from auto-tuning)
  integer, save :: mg_auto_phase = 0
  real(dp), save :: mg_cpu_ref = 0d0, mg_gpu_ref = 0d0
  logical, save :: mg_use_gpu = .false.
  logical, save :: mg_auto_init = .false.
  logical, save :: mg_orig_poisson = .false.
  integer(kind=8) :: mg_t1, mg_t2
  real(dp) :: mg_dt
#endif

  if(numbtot(1,ilevel)==0)return

  if(verbose)write(*,999)icount,ilevel, levelmin

#ifdef HYDRO_CUDA
  !---------------------------------------------------
  ! GPU auto-tuning: set flags at start of coarse step
  !---------------------------------------------------
  if(ilevel==levelmin) then
     ! Hydro auto-tuning: init + flag setting
     if(.not. hy_auto_init .and. gpu_hydro) then
        hy_auto_init = .true.
        hy_auto_phase = 0
     endif
     if(hy_auto_init) then
        if(hy_auto_phase == 0) then
           gpu_hydro = .false.
        else if(hy_auto_phase == 1) then
           gpu_hydro = .true.
        else
           gpu_hydro = hy_use_gpu
        endif
     endif
     ! Poisson auto-tuning: init is done near multigrid_fine call
     ! (mg_auto_init set there), flags set at levelmin entry
  endif
#endif

  !-------------------------------------------
  ! Make new refinements and update boundaries
  !-------------------------------------------
                               call timer('refine','start')
  if(levelmin.lt.nlevelmax .and..not. static)then
     if(ilevel==levelmin.or.icount>1)then
        do i=ilevel,nlevelmax
           if(i>levelmin)then

              !--------------------------
              ! Build communicators
              !--------------------------
              call build_comm(i)

              !--------------------------
              ! Update boundaries
              !--------------------------
              call make_virtual_fine_int(cpu_map(1),i)
              if(hydro)then
#ifdef SOLVERmhd
                 call make_virtual_fine_dp_bulk(uold,nvar+3,i)
#else
                 call make_virtual_fine_dp_bulk(uold,nvar,i)
#endif
                 if(simple_boundary)call make_boundary_hydro(i)
              end if
#ifdef RT
              if(rt)then
                 call make_virtual_fine_dp_bulk(rtuold,nrtvar,i)
                 if(simple_boundary)call rt_make_boundary_hydro(i)
              end if
#endif
              if(poisson)then
                 call make_virtual_fine_dp(phi(1),i)
                 call make_virtual_fine_dp_bulk(f,ndim,i)
                 if(simple_boundary)call make_boundary_force(i)
              end if
           end if

           !--------------------------
           ! Refine grids
           !--------------------------
           call refine_fine(i)
        end do
     end if
  end if

  !--------------------------
  ! Load balance
  !--------------------------
                               call timer('loadbalance','start')
  ok_defrag=.false.
  if(levelmin.lt.nlevelmax)then
     if(ilevel==levelmin)then
        if(varcpu_restart_done)then
           ! Force load balance after variable-ncpu restart
           if(myid==1) write(*,*) 'Forcing load_balance after variable-ncpu restart'
           call load_balance
           call defrag
           ok_defrag=.true.
           varcpu_restart_done=.false.
           first_step=.false.
        else if(nremap>0)then
           ! Skip first load balance because it has been performed before file dump
           if(nrestart>0.and.first_step)then
              first_step=.false.
           else
              if(MOD(nstep_coarse,nremap)==0)then
                 call load_balance
                 call defrag
                 ok_defrag=.true.
              endif
           end if
        end if
     endif
  end if

  !-----------------
  ! Particle leakage
  !-----------------
                               call timer('particles','start')
  call system_clock(pt_t1, pt_rate)
  if(pic)call make_tree_fine(ilevel)
  call system_clock(pt_t2)
  pt_mktree = pt_mktree + dble(pt_t2-pt_t1)/dble(pt_rate)
  
  !------------------------
  ! Output results to files
  !------------------------
   if(ilevel==levelmin) then
   endif


  if(ilevel==levelmin)then
     ! check if any of the processes received a signal for output
     call MPI_BARRIER(MPI_COMM_WORLD,mpi_err)
     call MPI_ALLREDUCE(output_now,output_now_all,1,MPI_LOGICAL,MPI_LOR,MPI_COMM_WORLD,mpi_err)




     if(mod(nstep_coarse,foutput)==0.or.aexp>=aout(iout).or.t>=tout(iout).or.output_now_all.EQV..true.)then
                               call timer('io','start')
        if(.not.ok_defrag)then
           call defrag
        endif

        call dump_all

        ! Run the clumpfinder, (produce output, don't keep arrays alive on output)
        if(clumpfind .and. ndim==3) call clump_finder(.true.,.false.)

        ! Dump lightcone
!jhshin1
!        if(lightcone) call output_cone()
!jhshin2
        if (output_now_all.EQV..true.) then
          output_now=.false.
        endif

     endif

  endif

  !----------------------------
  ! Output frame to movie dump (without synced levels)
  !----------------------------
  if(movie) then
     if(imov.le.imovout)then ! ifort returns error for next statement if looking
                             ! beyond what is allocated as an array amovout/tmovout
        if(aexp>=amovout(imov).or.t>=tmovout(imov))then
                               call timer('io','start')
           call output_frame()
        endif
     endif
  end if

  !-----------------------------------------------------------
  ! Put here all stuffs that are done only at coarse time step
  !-----------------------------------------------------------
  if(ilevel==levelmin)then

     !jaehyun

        if(aexp_old2 .le. 0.0001) aexp_old2=aexp_old
	 
	 if(spherical_region) then
        do i=1,5
	   if(hydro)   call output_sphere_hydro(i)
	   call output_sphere_part(i)
	   if(sink)   call output_sphere_sink(i)
        enddo
     endif

	 if(lightcone) then
        do i=1,2 
           if(hydro) call output_cone_hydro(i)
           call output_cone_part(i)
           if(sink) call output_cone_sink(i)
        enddo
	 endif
	 aexp_old2=aexp

     !jaehyun

                               call timer('feedback','start')
        !----------------------------------------------------
        ! Kinetic feedback
        !----------------------------------------------------
     if(hydro.and.star.and.f_w>0.)call kinetic_feedback
     
     call timer('sinks','start')
#ifdef HYDRO_CUDA
     ! --- GPU auto-tuning: set gpu_sink for this step ---
     if(.not. sk_auto_init .and. gpu_sink) then
        sk_auto_init = .true.
        sk_auto_phase = 0
     endif
     if(sk_auto_init) then
        if(sk_auto_phase == 0) then
           gpu_sink = .false.  ! Phase 0: force CPU path
        else if(sk_auto_phase == 1) then
           gpu_sink = .true.   ! Phase 1: force GPU path
        else
           gpu_sink = sk_use_gpu  ! Phase 2+: use decided path
        endif
     endif
#endif
     call system_clock(sk_t1)
     if(sink .and. sink_AGN)call AGN_feedback
     call system_clock(sk_t2)
     sk_agn_fb = sk_agn_fb + dble(sk_t2-sk_t1)/dble(pt_rate)
#ifdef HYDRO_CUDA
     ! --- GPU auto-tuning: record time and update decision ---
     if(sk_auto_init .and. sink .and. sink_AGN) then
        sk_dt_agn = dble(sk_t2-sk_t1)/dble(pt_rate)
        if(sk_auto_phase == 0) then
           ! Phase 0 done: recorded CPU time
           sk_cpu_ref = sk_dt_agn
           sk_auto_phase = 1
           if(myid==1) write(*,'(A,F8.3,A)') &
                ' [GPU auto-tune] AGN_feedback CPU: ', sk_cpu_ref, ' s'
        else if(sk_auto_phase == 1) then
           ! Phase 1 done: recorded GPU time, decide
           sk_gpu_ref = sk_dt_agn
           sk_use_gpu = (sk_gpu_ref < sk_cpu_ref)
           sk_auto_phase = 2
           gpu_sink = sk_use_gpu
           if(myid==1) then
              write(*,'(A,F8.3,A,F8.3,A)') &
                   ' [GPU auto-tune] AGN_feedback GPU: ', sk_gpu_ref, &
                   ' s  (CPU was ', sk_cpu_ref, ' s)'
              if(sk_use_gpu) then
                 write(*,'(A)') ' [GPU auto-tune] Decision: GPU (faster)'
              else
                 write(*,'(A)') ' [GPU auto-tune] Decision: CPU (faster)'
              endif
           endif
        else
           ! Phase 2+: keep booking, check for switch
           if(sk_use_gpu) then
              sk_gpu_ref = sk_dt_agn
           else
              sk_cpu_ref = sk_dt_agn
              ! Switch to GPU if CPU became slower
              if(sk_cpu_ref > sk_gpu_ref .and. sk_gpu_ref > 0d0) then
                 sk_use_gpu = .true.
                 gpu_sink = .true.
                 if(myid==1) write(*,'(A,F8.3,A,F8.3,A)') &
                      ' [GPU auto-tune] Switching to GPU: CPU=', &
                      sk_cpu_ref, ' s > GPU=', sk_gpu_ref, ' s'
              endif
           endif
        endif
     endif
#endif
     !-----------------------------------------------------
     ! Create sink particles and associated cloud particles
     !-----------------------------------------------------
     call system_clock(sk_t1)
     if(sink)call create_sink
     call system_clock(sk_t2)
     sk_create_sink = sk_create_sink + dble(sk_t2-sk_t1)/dble(pt_rate)

  endif



  !--------------------
  ! Poisson source term
  !--------------------
  if(poisson)then
                               call timer('poisson','start')
     !save old potential for time-extrapolation at level boundaries
     call save_phi_old(ilevel)
     call rho_fine(ilevel,icount)
  endif

  !-------------------------------------------
  ! Sort particles between ilevel and ilevel+1
  !-------------------------------------------
  if(pic)then
     ! Remove particles to finer levels
                               call timer('particles','start')
     call system_clock(pt_t1)
!    call MPI_BARRIER(MPI_COMM_WORLD,mpi_err)
     call kill_tree_fine(ilevel)
     ! Update boundary conditions for remaining particles
!    call MPI_BARRIER(MPI_COMM_WORLD,mpi_err)
     call virtual_tree_fine(ilevel)
     call system_clock(pt_t2)
     pt_killtree = pt_killtree + dble(pt_t2-pt_t1)/dble(pt_rate)
  end if
!############################################################
!       call getmem(real_mem)
!       call MPI_ALLREDUCE(real_mem,real_mem_tot,1,MPI_REAL,MPI_MAX,MPI_COMM_WORLD,info)
!       if(myid==1) then
!           call writemem(real_mem_tot)
!       endif
!############################################################

  !---------------
  ! Gravity update
  !---------------
  if(poisson)then
     call timer('poisson','start')
 
     ! Remove gravity source term with half time step and old force
!    call MPI_BARRIER(MPI_COMM_WORLD,mpi_err)
     if(hydro)then
        call synchro_hydro_fine(ilevel,-0.5*dtnew(ilevel))
     endif
!############################################################
!       call getmem(real_mem)
!       call MPI_ALLREDUCE(real_mem,real_mem_tot,1,MPI_REAL,MPI_MAX,MPI_COMM_WORLD,info)
!       if(myid==1) then
!           call writemem(real_mem_tot)
!       endif
!############################################################

     ! Compute gravitational potential
!    call MPI_BARRIER(MPI_COMM_WORLD,mpi_err)
!jhshin1
#ifdef HYDRO_CUDA
     ! --- GPU auto-tuning for Poisson MG (NOT fft) ---
     ! gpu_fft is excluded from auto-tuning because cuFFT direct solve
     ! and MG V-cycle produce different potential scales, so switching
     ! between them mid-run causes energy conservation failure.
     if(.not. mg_auto_init .and. gpu_poisson) then
        mg_auto_init = .true.
        mg_auto_phase = 0
        mg_orig_poisson = gpu_poisson
     endif
     if(mg_auto_init .and. ilevel==levelmin) then
        if(mg_auto_phase == 0) then
           gpu_poisson = .false.
        else if(mg_auto_phase == 1) then
           gpu_poisson = mg_orig_poisson
        else
           if(.not. mg_use_gpu) then
              gpu_poisson = .false.
           endif
        endif
     endif
#endif
#ifdef HYDRO_CUDA
     call system_clock(mg_t1)
#endif
     if(ilevel>levelmin)then
        if(ilevel .ge. cg_levelmin) then
           call timer('poisson - cg', 'start')
           call phi_fine_cg(ilevel,icount)
        else
           call timer('poisson - mg AMR', 'start')
           call multigrid_fine(ilevel,icount)
        end if
     else
#ifdef USE_FFTW
        if(use_fftw) then
           call timer('poisson-fftw3 base','start')
        else
#endif
#ifdef HYDRO_CUDA
        if(gpu_fft .and. cuda_pool_is_initialized_c()/=0) then
           call timer('poisson-cuFFT base','start')
        else
#endif
           call timer('poisson - mg base', 'start')
#ifdef HYDRO_CUDA
        end if
#endif
#ifdef USE_FFTW
        end if
#endif
        call multigrid_fine(levelmin,icount)
     end if
#ifdef HYDRO_CUDA
     call system_clock(mg_t2)
     if(mg_auto_init) then
        mg_dt = dble(mg_t2-mg_t1)/dble(pt_rate)
        if(mg_auto_phase == 0) then
           mg_cpu_ref = mg_cpu_ref + mg_dt
        else if(mg_auto_phase == 1) then
           mg_gpu_ref = mg_gpu_ref + mg_dt
        endif
     endif
#endif
!############################################################
!       call getmem(real_mem)
!       call MPI_ALLREDUCE(real_mem,real_mem_tot,1,MPI_REAL,MPI_MAX,MPI_COMM_WORLD,info)
!       if(myid==1) then
!           call writemem(real_mem_tot)
!       endif
!############################################################
      call timer('poisson', 'start')
     !when there is no old potential...
!jhshin2      
     
     if (nstep==0)call save_phi_old(ilevel)
!############################################################
!       call getmem(real_mem)
!       call MPI_ALLREDUCE(real_mem,real_mem_tot,1,MPI_REAL,MPI_MAX,MPI_COMM_WORLD,info)
!       if(myid==1) then
!           call writemem(real_mem_tot)
!       endif
!############################################################

     ! Compute gravitational acceleration
!    call MPI_BARRIER(MPI_COMM_WORLD,mpi_err)
     call force_fine(ilevel,icount)
!############################################################
!       call getmem(real_mem)
!       call MPI_ALLREDUCE(real_mem,real_mem_tot,1,MPI_REAL,MPI_MAX,MPI_COMM_WORLD,info)
!       if(myid==1) then
!           call writemem(real_mem_tot)
!       endif
!############################################################

     ! Synchronize remaining particles for gravity
     if(pic)then
                               call timer('particles','start')
        call system_clock(pt_t1)
!    call MPI_BARRIER(MPI_COMM_WORLD,mpi_err)
        call synchro_fine(ilevel)
        call system_clock(pt_t2)
        pt_synchro = pt_synchro + dble(pt_t2-pt_t1)/dble(pt_rate)
     end if
!############################################################
!       call getmem(real_mem)
!       call MPI_ALLREDUCE(real_mem,real_mem_tot,1,MPI_REAL,MPI_MAX,MPI_COMM_WORLD,info)
!       if(myid==1) then
!           call writemem(real_mem_tot)
!       endif
!############################################################

     if(hydro)then
                               call timer('poisson','start')

!    call MPI_BARRIER(MPI_COMM_WORLD,mpi_err)
        ! Add gravity source term with half time step and new force
        call synchro_hydro_fine(ilevel,+0.5*dtnew(ilevel))
!############################################################
!       call getmem(real_mem)
!       call MPI_ALLREDUCE(real_mem,real_mem_tot,1,MPI_REAL,MPI_MAX,MPI_COMM_WORLD,info)
!       if(myid==1) then
!           call writemem(real_mem_tot)
!       endif
!############################################################


        ! Density threshold and/or Bondi accretion onto sink particle
                               call timer('sinks','start')
        call system_clock(sk_t1)
        if(sink)then
           if(bondi)then
              call grow_bondi(ilevel)
           else
              call grow_jeans(ilevel)
           endif
        endif
        call system_clock(sk_t2)
        sk_grow = sk_grow + dble(sk_t2-sk_t1)/dble(pt_rate)
!############################################################
!       call getmem(real_mem)
!       call MPI_ALLREDUCE(real_mem,real_mem_tot,1,MPI_REAL,MPI_MAX,MPI_COMM_WORLD,info)
!       if(myid==1) then
!           call writemem(real_mem_tot)
!       endif
!############################################################

        ! Update boundaries
                               call timer('hydro - ghostzones','start')
#ifdef SOLVERmhd
        call make_virtual_fine_dp_bulk(uold,nvar+3,ilevel)
#else
        call make_virtual_fine_dp_bulk(uold,nvar,ilevel)
#endif
        if(simple_boundary)call make_boundary_hydro(ilevel)
        
     end if
  end if

#ifdef RT
  ! Turn on RT in case of rt_stars and first stars just created:
  ! Update photon packages according to star particles
                               call timer('radiative transfer','start')
  if(rt .and. rt_star) call update_star_RT_feedback(ilevel)
#endif

  !----------------------
  ! Compute new time step
  !----------------------
                               call timer('courant','start')
  call newdt_fine(ilevel)
  if(ilevel>levelmin)then
     dtnew(ilevel)=MIN(dtnew(ilevel-1)/real(nsubcycle(ilevel-1)),dtnew(ilevel))
  end if

  ! Set unew equal to uold
                               call timer('hydro - set unew','start')
  if(hydro)call set_unew(ilevel)

#ifdef RT
  ! Set rtunew equal to rtuold
                               call timer('radiative transfer','start')
  if(rt)call rt_set_unew(ilevel)
#endif

  !---------------------------
  ! Recursive call to amr_step
  !---------------------------
  if(ilevel<nlevelmax)then
     if(numbtot(1,ilevel+1)>0)then
        if(nsubcycle(ilevel)==2)then
           call amr_step(ilevel+1,1)
           call amr_step(ilevel+1,2)
        else
           call amr_step(ilevel+1,1)
        endif
     else 
        ! Otherwise, update time and finer level time-step
        dtold(ilevel+1)=dtnew(ilevel)/dble(nsubcycle(ilevel))
        dtnew(ilevel+1)=dtnew(ilevel)/dble(nsubcycle(ilevel))
        call update_time(ilevel)
     end if
  else
     call update_time(ilevel)
  end if

  ! Thermal feedback from stars (also call if no feedback, for bookkeeping)
  if(hydro.and.star) then
                               call timer('feedback','start')
     call thermal_feedback(ilevel)
  endif
  
  !-----------
  ! Hydro step
  !-----------
  if(hydro)then

     ! Hyperbolic solver
                               call timer('hydro - godunov','start')
#ifdef HYDRO_CUDA
     ! --- GPU auto-tuning for hydro: flags set at levelmin entry ---
     ! (actual flag setting is done at start of amr_step(levelmin))
#endif
#ifdef HYDRO_CUDA
     call system_clock(hy_t1)
#endif
     call godunov_fine(ilevel)
#ifdef HYDRO_CUDA
     call system_clock(hy_t2)
     if(hy_auto_init) then
        hy_dt = dble(hy_t2-hy_t1)/dble(pt_rate)
        if(hy_auto_phase == 0) then
           hy_cpu_ref = hy_cpu_ref + hy_dt
        else if(hy_auto_phase == 1) then
           hy_gpu_ref = hy_gpu_ref + hy_dt
        endif
     endif
#endif

     ! Reverse update boundaries
                               call timer('hydro - rev ghostzones','start')
#ifdef SOLVERmhd
     call make_virtual_reverse_dp_bulk(unew,nvar+3,ilevel)
#else
     call make_virtual_reverse_dp_bulk(unew,nvar,ilevel)
#endif
     if(pressure_fix)then
        call make_virtual_reverse_dp(enew(1),ilevel)
        call make_virtual_reverse_dp(divu(1),ilevel)
     endif

     ! Set uold equal to unew
                               call timer('hydro - set uold','start')
     call set_uold(ilevel)

     ! Add gravity source term with half time step and old force
     ! in order to complete the time step 
!    call MPI_BARRIER(MPI_COMM_WORLD,mpi_err)
                               call timer('poisson','start')
     if(poisson)call synchro_hydro_fine(ilevel,+0.5*dtnew(ilevel))

     ! Restriction operator
                               call timer('hydro upload fine','start')
     call upload_fine(ilevel)

  endif
 

  !---------------------
  ! Do RT/Chemistry step
  !---------------------
#ifdef RT
  if(rt .and. rt_advect) then  
                               call timer('radiative transfer','start')
     call rt_step(ilevel)
  else
     ! Still need a chemistry call if RT is defined but not
     ! actually doing radiative transfer (i.e. rt==false):
                               call timer('cooling','start')
     if(neq_chem.or.cooling.or.T2_star>0.0)call cooling_fine(ilevel)
  endif
  ! Regular updates and book-keeping:
  if(ilevel==levelmin) then
                               call timer('radiative transfer','start')
     if(cosmo) call update_rt_c
     if(cosmo .and. haardt_madau) call update_UVrates(aexp)
     if(cosmo .and. rt_isDiffuseUVsrc) call update_UVsrc
                               call timer('cooling','start')
     if(cosmo) call update_coolrates_tables(dble(aexp))
                               call timer('radiative transfer','start')
     if(ilevel==levelmin) call output_rt_stats
  endif
#else
                               call timer('cooling','start')
  if(neq_chem.or.cooling.or.T2_star>0.0)call cooling_fine(ilevel)
#endif
  
  !---------------
  ! Move particles
  !---------------
  if(pic)then
                               call timer('particles','start')
     call system_clock(pt_t1)
     call move_fine(ilevel) ! Only remaining particles
     call system_clock(pt_t2)
     pt_move = pt_move + dble(pt_t2-pt_t1)/dble(pt_rate)
  end if

  !----------------------------------
  ! Star formation in leaf cells only
  !----------------------------------
                               call timer('feedback','start')
  if(hydro.and.star)call star_formation(ilevel)

  ! Compute Bondi-Hoyle accretion parameters
                               call timer('sinks','start')
  call system_clock(sk_t1)
  if(sink.and.bondi)call bondi_hoyle(ilevel)
  call system_clock(sk_t2)
  sk_bondi_hoyle = sk_bondi_hoyle + dble(sk_t2-sk_t1)/dble(pt_rate)

  !---------------------------------------
  ! Update physical and virtual boundaries
  !---------------------------------------
  if(hydro)then
                               call timer('hydro - ghostzones','start')
#ifdef SOLVERmhd
     call make_virtual_fine_dp_bulk(uold,nvar+3,ilevel)
#else
     call make_virtual_fine_dp_bulk(uold,nvar,ilevel)
#endif
     if(simple_boundary)call make_boundary_hydro(ilevel)
  endif

#ifdef SOLVERmhd
  ! Magnetic diffusion step
 if(hydro)then
     if(eta_mag>0d0.and.ilevel==levelmin)then
                               call timer('hydro - diffusion','start')
        call diffusion
     endif
  end if
#endif

  !-----------------------
  ! Compute refinement map
  !-----------------------
                               call timer('flag','start')
  if(.not.static) call flag_fine(ilevel,icount)


  !----------------------------
  ! Merge finer level particles
  !----------------------------
                               call timer('particles','start')
  call system_clock(pt_t1)
  if(pic)call merge_tree_fine(ilevel)
  call system_clock(pt_t2)
  pt_merge = pt_merge + dble(pt_t2-pt_t1)/dble(pt_rate)

  !---------------
  ! Radiation step
  !---------------
#ifdef ATON
  if(aton.and.ilevel==levelmin)then
                               call timer('aton','start')
     call rad_step(dtnew(ilevel))
  endif
#endif

  !-------------------------------
  ! Update coarser level time-step
  !-------------------------------
  if(ilevel>levelmin)then
     if(nsubcycle(ilevel-1)==1)dtnew(ilevel-1)=dtnew(ilevel)
     if(icount==2)dtnew(ilevel-1)=dtold(ilevel)+dtnew(ilevel)
  end if

#ifdef HYDRO_CUDA
  ! --- GPU auto-tuning: end-of-coarse-step decision for hydro & poisson ---
  if(ilevel==levelmin) then
     ! Hydro auto-tuning decision
     if(hy_auto_init) then
        if(hy_auto_phase == 0) then
           hy_auto_phase = 1
           if(myid==1) write(*,'(A,F8.3,A)') &
                ' [GPU auto-tune] Hydro CPU: ', hy_cpu_ref, ' s'
        else if(hy_auto_phase == 1) then
           hy_use_gpu = (hy_gpu_ref < hy_cpu_ref)
           hy_auto_phase = 2
           gpu_hydro = hy_use_gpu
           if(myid==1) then
              write(*,'(A,F8.3,A,F8.3,A)') &
                   ' [GPU auto-tune] Hydro GPU: ', hy_gpu_ref, &
                   ' s  (CPU was ', hy_cpu_ref, ' s)'
              if(hy_use_gpu) then
                 write(*,'(A)') ' [GPU auto-tune] Hydro decision: GPU (faster)'
              else
                 write(*,'(A)') ' [GPU auto-tune] Hydro decision: CPU (faster)'
              endif
           endif
        else
           ! Phase 2+: monitor accumulated time per coarse step
           if(hy_use_gpu) then
              hy_gpu_ref = 0d0  ! will re-accumulate next step
           else
              hy_cpu_ref = 0d0
           endif
        endif
     endif
     ! Poisson MG auto-tuning decision (gpu_fft excluded)
     if(mg_auto_init) then
        if(mg_auto_phase == 0) then
           mg_auto_phase = 1
           if(myid==1) write(*,'(A,F8.3,A)') &
                ' [GPU auto-tune] Poisson MG CPU: ', mg_cpu_ref, ' s'
        else if(mg_auto_phase == 1) then
           mg_use_gpu = (mg_gpu_ref < mg_cpu_ref)
           mg_auto_phase = 2
           if(.not. mg_use_gpu) then
              gpu_poisson = .false.
           endif
           if(myid==1) then
              write(*,'(A,F8.3,A,F8.3,A)') &
                   ' [GPU auto-tune] Poisson MG GPU: ', mg_gpu_ref, &
                   ' s  (CPU was ', mg_cpu_ref, ' s)'
              if(mg_use_gpu) then
                 write(*,'(A)') ' [GPU auto-tune] Poisson MG decision: GPU (faster)'
              else
                 write(*,'(A)') ' [GPU auto-tune] Poisson MG decision: CPU (faster)'
              endif
           endif
        else
           ! Phase 2+: monitor
           if(mg_use_gpu) then
              mg_gpu_ref = 0d0
           else
              mg_cpu_ref = 0d0
           endif
        endif
     endif
  endif
#endif

  ! Print particle sub-timers at last step
  if(ilevel==levelmin .and. nstep>=nstepmax .and. myid==1) then
     write(*,'(A)') ' === Particle sub-timers ==='
     write(*,'(A,F8.3,A)') '   make_tree  : ', pt_mktree, ' s'
     write(*,'(A,F8.3,A)') '   kill+virt  : ', pt_killtree, ' s'
     write(*,'(A,F8.3,A)') '   synchro    : ', pt_synchro, ' s'
     write(*,'(A,F8.3,A)') '   move       : ', pt_move, ' s'
     write(*,'(A,F8.3,A)') '   merge      : ', pt_merge, ' s'
     write(*,'(A,F8.3,A)') '   TOTAL      : ', &
          pt_mktree+pt_killtree+pt_synchro+pt_move+pt_merge, ' s'
     write(*,'(A)') ' === Sink sub-timers ==='
     write(*,'(A,F8.3,A)') '   AGN_feedback : ', sk_agn_fb, ' s'
     write(*,'(A,F8.3,A)') '   create_sink  : ', sk_create_sink, ' s'
     write(*,'(A,F8.3,A)') '   grow_bondi   : ', sk_grow, ' s'
     write(*,'(A,F8.3,A)') '   bondi_hoyle  : ', sk_bondi_hoyle, ' s'
     write(*,'(A,F8.3,A)') '   TOTAL        : ', &
          sk_agn_fb+sk_create_sink+sk_grow+sk_bondi_hoyle, ' s'
  end if

999 format(' Entering amr_step',i1,' for level',i2, '  for a levelmin ',i3)

end subroutine amr_step

!##########################################################################
!##########################################################################
!##########################################################################
!##########################################################################

#ifdef RT
subroutine rt_step(ilevel)
  use amr_parameters, only: dp
  use amr_commons,    only: levelmin, t, dtnew, myid
  use rt_parameters, only: rt_isDiffuseUVsrc
  use rt_cooling_module, only: update_UVrates
  use rt_hydro_commons
  use UV_module
  use SED_module,     only: star_RT_feedback
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer, intent(in) :: ilevel

!--------------------------------------------------------------------------
!  Radiative transfer and chemistry step. Either do one step on ilevel,
!  with radiation field updates in coarser level neighbours, or, if
!  rt_nsubsteps>1, do many substeps in ilevel only, using Dirichlet
!  boundary conditions for the level boundaries. 
!--------------------------------------------------------------------------

  real(dp) :: dt_hydro, t_left, dt_rt, t_save
  integer  :: i_substep, ivar

  dt_hydro = dtnew(ilevel)                   ! Store hydro timestep length
  t_left = dt_hydro
  ! We shift the time backwards one hydro-dt, to get evolution of stellar
  ! ages within the hydro timestep, in the case of rt subcycling:
  t_save=t ; t=t-t_left
  
  i_substep = 0
  do while (t_left > 0)                      !                RT sub-cycle
     i_substep = i_substep + 1
     call get_rt_courant_coarse(dt_rt)
     ! Temporarily change timestep length to rt step:
     dtnew(ilevel) = MIN(t_left, dt_rt/2.0**(ilevel-levelmin))
     t = t + dtnew(ilevel) ! Shift the time forwards one dt_rt

     ! If (myid==1) write(*,900) dt_hydro, dtnew(ilevel), i_substep, ilevel    
     if (i_substep > 1) call rt_set_unew(ilevel)

     if(rt_star) call star_RT_feedback(ilevel,dtnew(ilevel))

     ! Hyperbolic solver
     if(rt_advect) call rt_godunov_fine(ilevel,dtnew(ilevel))

     call add_rt_sources(ilevel,dtnew(ilevel))

     ! Reverse update boundaries
     do ivar=1,nrtvar
        call make_virtual_reverse_dp(rtunew(1,ivar),ilevel)
     end do

     ! Set rtuold equal to rtunew
     call rt_set_uold(ilevel)

     if(neq_chem.or.cooling.or.T2_star>0.0)call cooling_fine(ilevel)
     
     do ivar=1,nrtvar
        call make_virtual_fine_dp(rtuold(1,ivar),ilevel)
     end do
     if(simple_boundary)call rt_make_boundary_hydro(ilevel)

     t_left = t_left - dtnew(ilevel)
  end do                                   !          End RT subcycle loop
  dtnew(ilevel) = dt_hydro                 ! Restore hydro timestep length
  t = t_save       ! Restore original time (otherwise tiny roundoff error)
  
  ! Restriction operator to update coarser level split cells
  call rt_upload_fine(ilevel)

  if (myid==1 .and. rt_nsubcycle .gt. 1) write(*,901) ilevel, i_substep

900 format (' dt_hydro=', 1pe12.3, ' dt_rt=', 1pe12.3, ' i_sub=', I5, ' level=', I5)
901 format (' Performed level', I3, ' RT-step with ', I5, ' subcycles')
  
end subroutine rt_step
#endif
