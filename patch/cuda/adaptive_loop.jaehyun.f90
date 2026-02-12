subroutine adaptive_loop
  use amr_commons
  use hydro_commons
  use pm_commons
  use poisson_commons
  use cooling_module
#ifdef RT
  use rt_hydro_commons
#endif
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer(kind=8)::n_step
  integer::ilevel,idim,ivar,info,tot_pt
!jhshin1
  real(kind=8)::tt1,tt2,muspt,muspt_this_step,wallsec,dumpsec
!jhshin2
  real(kind=4)::real_mem,real_mem_tot,real_mem_min,real_mem_max
!jhshin1
  real(kind=8),save::tstart=0.0
!jhshin2

#ifndef WITHOUTMPI
  tt1=MPI_WTIME()
!jhshin1
  ! for calculating total run time
  if (tstart.eq.0.0) then
     tstart = MPI_WTIME()
  end if
!jhshin2
#endif

  call init_amr                      ! Initialize AMR variables
  call init_time                     ! Initialize time variables
  if(hydro)call init_hydro           ! Initialize hydro variables
#ifdef RT
  if(rt.or.neq_chem) &
       & call rt_init_hydro          ! Initialize radiation variables
#endif
  if(poisson)call init_poisson       ! Initialize poisson variables
#ifdef ATON
  if(aton)call init_radiation        ! Initialize radiation variables
#endif
  if(nrestart==0)call init_refine    ! Build initial AMR grid

#ifdef grackle
  if(cosmo)then
     ! Compute cooling table at current aexp
  endif
#else  
  if(cooling.and..not.neq_chem) &
       call set_table(dble(aexp))    ! Initialize cooling look up table
#endif
  if(pic)call init_part              ! Initialize particle variables
  if(pic)call init_tree              ! Initialize particle tree
  if(nrestart==0)call init_refine_2  ! Build initial AMR grid again

  ! Rebuild Morton hash tables after full AMR grid is available
  if(nrestart==0) call morton_hash_build_and_verify()

#ifndef WITHOUTMPI
  muspt=0.
  tot_pt=-1
  tt2=MPI_WTIME()
  if(myid==1)write(*,*)'Time elapsed since startup:',tt2-tt1
#endif


  if(myid==1)then
     write(*,*)'Initial mesh structure'
     do ilevel=1,nlevelmax
        if(numbtot(1,ilevel)>0)write(*,999)ilevel,numbtot(1:4,ilevel)
     end do
  end if

  nstep_coarse_old=nstep_coarse

  if(myid==1)write(*,*)'Starting time integration' 

  do ! Main time loop
                               call timer('coarse levels','start')

#ifndef WITHOUTMPI
     tt1=MPI_WTIME()
#endif

     if(verbose)write(*,*)'Entering amr_step_coarse'

     epot_tot=0.0D0  ! Reset total potential energy
     ekin_tot=0.0D0  ! Reset total kinetic energy
     mass_tot=0.0D0  ! Reset total mass
     eint_tot=0.0D0  ! Reset total internal energy
#ifdef SOLVERmhd
     emag_tot=0.0D0  ! Reset total magnetic energy
#endif

     ! Make new refinements
     if(levelmin.lt.nlevelmax.and.(.not.static.or.(nstep_coarse_old.eq.nstep_coarse.and.restart_remap)))then
        call refine_coarse
        do ilevel=1,levelmin
           call build_comm(ilevel)
           call make_virtual_fine_int(cpu_map(1),ilevel)
           if(hydro)then
#ifdef SOLVERmhd
              call make_virtual_fine_dp_bulk(uold,nvar+3,ilevel)
#else
              call make_virtual_fine_dp_bulk(uold,nvar,ilevel)
#endif
              if(simple_boundary)call make_boundary_hydro(ilevel)
           endif
#ifdef RT
           if(rt)then
              call make_virtual_fine_dp_bulk(rtuold,nrtvar,ilevel)
              if(simple_boundary)call rt_make_boundary_hydro(ilevel)
           endif
#endif
           if(poisson)then
              call make_virtual_fine_dp(phi(1),ilevel)
              call make_virtual_fine_dp_bulk(f,ndim,ilevel)
           end if
           if(ilevel<levelmin)call refine_fine(ilevel)
        end do
     endif

     ! Call base level
     call amr_step(levelmin,1)

     ! Rebuild and verify Morton hash tables after all levels updated
     call morton_hash_rebuild
     call morton_hash_verify('step')

                               call timer('coarse levels','start')

     if(levelmin.lt.nlevelmax.and.(.not.static.or.(nstep_coarse_old.eq.nstep_coarse.and.restart_remap)))then

        do ilevel=levelmin-1,1,-1
           ! Hydro book-keeping
           if(hydro)then
              call upload_fine(ilevel)
#ifdef SOLVERmhd
              call make_virtual_fine_dp_bulk(uold,nvar+3,ilevel)
#else
              call make_virtual_fine_dp_bulk(uold,nvar,ilevel)
#endif
              if(simple_boundary)call make_boundary_hydro(ilevel)
           end if
#ifdef RT
           ! Radiation book-keeping
           if(rt)then
              call rt_upload_fine(ilevel)
              call make_virtual_fine_dp_bulk(rtuold,nrtvar,ilevel)
              if(simple_boundary)call rt_make_boundary_hydro(ilevel)
           end if
#endif
           ! Gravity book-keeping
           if(poisson)then
              call make_virtual_fine_dp(phi(1),ilevel)
              call make_virtual_fine_dp_bulk(f,ndim,ilevel)
           end if
        end do
        
        ! Build refinement map
        do ilevel=levelmin-1,1,-1
           call flag_fine(ilevel,2)
        end do
        call flag_coarse
!		if(lightcone) then
!			call output_cone_hydro(1)
!			call output_cone_part(1)
!			call output_cone_sink(1)
!			call output_cone_hydro(2)
!			call output_cone_part(2)
!			call output_cone_sink(2)
!		endif

 !       aexp_old2=aexp_old
!		if(myid .eq.1) write(*,*)aexp_old
!		if(myid .eq.1) write(*,*)aexp
     endif

     ! New coarse time-step
     nstep_coarse=nstep_coarse+1

!jhshin1 -walltime setting
#ifndef WITHOUTMPI
     tt2=MPI_WTIME()
     if(mod(nstep_coarse,ncontrol)==0)then
        call getmem(real_mem)
        call MPI_ALLREDUCE(real_mem,real_mem_max,1,MPI_REAL,MPI_MAX,MPI_COMM_WORLD,info)
        call MPI_ALLREDUCE(real_mem,real_mem_min,1,MPI_REAL,MPI_MIN,MPI_COMM_WORLD,info)
        real_mem_tot=real_mem_max
        if(myid==1)then
           if (tot_pt==0) muspt=0. ! dont count first timestep
           n_step = int(numbtot(1,levelmin),kind=8)*twotondim
           do ilevel=levelmin+1,nlevelmax
             n_step = n_step + int(numbtot(1,ilevel),kind=8)*product(nsubcycle(levelmin:ilevel-1))*(twotondim-1)
           enddo
           muspt_this_step = (tt2-tt1)*1e6/n_step*ncpu
           muspt = muspt + muspt_this_step
           tot_pt = tot_pt + 1
           write(*,'(a,f8.2,a,f12.2,a,f12.2,a)')' Time elapsed since last coarse step:',tt2-tt1 &
          ,' s',muspt_this_step,' mus/pt'  &
          ,muspt / max(tot_pt,1), ' mus/pt (av)'
           call writemem(real_mem_tot)
           call writemem_minmax(real_mem_min,real_mem_max)
           write(*,*)'Total running time:', NINT((tt2-tstart)*100.0)*0.01,'s'
        endif
        if(walltime_hrs.gt.0d0) then
           wallsec = walltime_hrs*3600.     ! Convert from hours to seconds
           dumpsec = minutes_dump*60.       ! Convert minutes before end to seconds
           if(wallsec-dumpsec.lt.tt2-tstart) then
              output_now=.true.
              if(myid==1) write(*,*) 'Dumping snapshot before walltime runs out'
              ! Now set walltime to a negative number so we don't keep printing
              ! outputs
              walltime_hrs = -1d0
           endif
        endif
     endif
#endif
!jhshin2

  end do

999 format(' Level ',I2,' has ',I10,' grids (',3(I8,','),')')

end subroutine adaptive_loop
