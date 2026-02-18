subroutine init_amr
  use amr_commons
  use hydro_commons
  use pm_commons  
  use poisson_commons
  use bisection
  use ksection
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'  
#endif
  integer::i,idim,ncell,iskip,ind,ncache,ilevel,ibound,nboundary2
  integer::ncpu2,ndim2,nx2,ny2,nz2,ngridmax2,nlevelmax2
  integer::noutput2,iout2,ifout2,ilun,info
  integer::ix,iy,iz,ix_max,iy_max,iz_max,nxny,nx_loc
  real(dp)::mass_sph2 
  integer,dimension(:),allocatable::ind_grid,iig,pos,grid
  real(dp),dimension(1:MAXOUT)::aout2=1.1d0 
  real(dp),dimension(1:MAXOUT)::tout2=0.0d0 
  real(dp),dimension(:),allocatable::xxg
  integer ,dimension(1:nvector)::c
  real(dp),dimension(1:nvector,1:ndim)::x
  real(qdp),dimension(1:nvector)::order_min,order_max
  logical::ok
  real(dp)::dx_loc,scale
  character(LEN=128)::ordering2
  character(LEN=80)::fileloc
  character(LEN=5)::nchar,ncharcpu
  real(dp),allocatable,dimension(:)::bxmin,bxmax
  integer,parameter::tag=1100
  integer::dummy_io,info2
  real(kind=8),allocatable,dimension(:)::bound_key_restart
  
  if(verbose.and.myid==1)write(*,*)'Entering init_amr'

  ! Constants
  ncoarse=nx*ny*nz
  ncell=ncoarse+twotondim*ngridmax
  nxny=nx*ny
  ix_max=0; iy_max=0; iz_max=0
  if(ndim>0)ix_max=1 
  if(ndim>1)iy_max=1
  if(ndim>2)iz_max=1
  nx_loc=icoarse_max-icoarse_min+1
  scale=boxlen/dble(nx_loc)

  ! Initial time step for each level
  dtold=0.0D0
  dtnew=0.0D0

  ! Allocate AMR cell-based arrays
  allocate(flag1(0:ncell)) ! Note: starting from 0
  allocate(flag2(0:ncell)) ! Note: starting from 0
  allocate(son  (1:ncell)) ! Son index
  flag1=0; flag2=0; son=0

  ! Allocate MPI cell-based arrays
  allocate(cpu_map    (1:ncell)) ! Cpu map
  allocate(cpu_map2   (1:ncell)) ! New cpu map for load balance
  cpu_map=0; cpu_map2=0
  if(ordering=='ksection') then
     allocate(hilbert_key(1:1))     ! Defrag uses local scratch
  else
     allocate(hilbert_key(1:ncell)) ! Ordering key
  end if
  hilbert_key=0.0d0

  ! Bisection ordering: compute array boundaries and
  ! allocate arrays if needed
  nbilevelmax=ceiling(log(dble(ncpu))/log(2.0))
  nbinodes=2**(nbilevelmax+1)-1
  nbileafnodes=2**nbilevelmax
  bisec_nres=2**(nlevelmax+1)
  bisec_res=scale/dble(bisec_nres)

  if(ordering=='bisection') then
    use_cpubox_decomp=.true.
    ! allocate bisection tree structure
    allocate(bisec_wall(1:nbinodes))
    allocate(bisec_next(1:nbinodes,1:2))
    allocate(bisec_indx(1:nbinodes))
    bisec_wall=0.0d0; bisec_next=0; bisec_indx=0; bisec_root=0
    ! allocate some other bisection stuff
    allocate(bisec_cpubox_min (1:ncpu,1:ndim))
    allocate(bisec_cpubox_max (1:ncpu,1:ndim))
    allocate(bisec_cpubox_min2(1:ncpu,1:ndim))
    allocate(bisec_cpubox_max2(1:ncpu,1:ndim))
    allocate(bisec_cpu_load(1:ncpu))
    bisec_cpubox_min=0;  bisec_cpubox_max=0;
    bisec_cpubox_min2=0; bisec_cpubox_max2=0;
    bisec_cpu_load=0;
    ! allocate histograms
    allocate(bisec_hist(1:nbileafnodes,1:bisec_nres))
    allocate(bisec_hist_bounds(1:(nbileafnodes+1)))
    allocate(new_hist_bounds  (1:(nbileafnodes+1)))
    allocate(bisec_ind_cell(1:ncell))    ! big array
    allocate(cell_level    (1:ncell))    ! big array
    bisec_hist=0
    bisec_hist_bounds=0; new_hist_bounds=0
    bisec_ind_cell=0; cell_level=0

  else if(ordering=='ksection') then

    use_cpubox_decomp=.true.
    ! Compute prime factorization and split sequence
    call compute_ksection_factorization(ncpu)
    ! Compute split directions (longest axis per level)
    call compute_ksection_directions(scale)
    ! Compute tree node count
    ksec_nbinodes=compute_ksec_nbinodes()
    if(myid==1) then
       write(*,*) 'K-section decomposition:'
       write(*,*) '   nksec_levels  =', nksec_levels
       write(*,*) '   ksec_kmax     =', ksec_kmax
       write(*,*) '   ksec_nbinodes =', ksec_nbinodes
       write(*,*) '   ksec_factor   =', ksec_factor(1:nksec_levels)
       write(*,*) '   ksec_dir      =', ksec_dir(1:nksec_levels)
    end if
    ! Allocate k-section tree structure
    allocate(ksec_wall(1:ksec_nbinodes,1:max(ksec_kmax-1,1)))
    allocate(ksec_next(1:ksec_nbinodes,1:ksec_kmax))
    allocate(ksec_indx(1:ksec_nbinodes))
    ksec_wall=0.0d0; ksec_next=0; ksec_indx=0; ksec_root=1
    ! Allocate k-section tree navigation arrays
    allocate(ksec_cpumin(1:ksec_nbinodes))
    allocate(ksec_cpumax(1:ksec_nbinodes))
    ksec_cpumin=0; ksec_cpumax=0
    ! Allocate cpubox arrays (reused for authorize_fine compatibility)
    allocate(bisec_cpubox_min (1:ncpu,1:ndim))
    allocate(bisec_cpubox_max (1:ncpu,1:ndim))
    allocate(bisec_cpubox_min2(1:ncpu,1:ndim))
    allocate(bisec_cpubox_max2(1:ncpu,1:ndim))
    allocate(bisec_cpu_load(1:ncpu))
    bisec_cpubox_min=0;  bisec_cpubox_max=0;
    bisec_cpubox_min2=0; bisec_cpubox_max2=0;
    bisec_cpu_load=0;
    ! Allocate histograms (reuse bisection histogram arrays)
    allocate(bisec_hist(1:nbileafnodes,1:bisec_nres))
    allocate(bisec_hist_bounds(1:(nbileafnodes+1)))
    allocate(new_hist_bounds  (1:(nbileafnodes+1)))
    ! bisec_ind_cell and cell_level: on-demand in init_bisection_histogram
    bisec_hist=0
    bisec_hist_bounds=0; new_hist_bounds=0
  end if

 bisection_or_ordering:if(.not.use_cpubox_decomp) then ! use usual ordering machinery

    ! Cpu boundaries in chosen ordering
    ndomain=ncpu*overload
    allocate(bound_key (0:ndomain))
    allocate(bound_key2(0:ndomain))

    ! Compute minimum and maximum ordering key
    dx_loc=scale
    x(1,1)=0.5*scale
#if NDIM>1
    x(1,2)=0.5*scale
#endif
#if NDIM>2
    x(1,3)=0.5*scale
#endif
    call cmp_minmaxorder(x,order_min,order_max,dx_loc,1)
    order_all_min=order_min(1)
    order_all_max=order_max(1)
    do iz=kcoarse_min,kcoarse_max
       do iy=jcoarse_min,jcoarse_max
          do ix=icoarse_min,icoarse_max
             ind=1+ix+iy*nx+iz*nxny
             x(1,1)=(dble(ix)+0.5d0-dble(icoarse_min))*scale
#if NDIM>1
             x(1,2)=(dble(iy)+0.5d0-dble(jcoarse_min))*scale
#endif
#if NDIM>2
             x(1,3)=(dble(iz)+0.5d0-dble(kcoarse_min))*scale
#endif
             call cmp_minmaxorder(x,order_min,order_max,dx_loc,1)
             order_all_min=min(order_all_min,order_min(1))
             order_all_max=max(order_all_max,order_max(1))
          end do
       end do
    end do

    ! Set initial cpu boundaries
    do i=0,ndomain-1
#ifdef QUADHILBERT
       bound_key(i)=order_all_min+real(i,16)/real(ndomain,16)* &
            & (order_all_max-order_all_min)
#else
       bound_key(i)=order_all_min+real(i,8)/real(ndomain,8)* &
            & (order_all_max-order_all_min)
#endif
    end do
    bound_key(ndomain)=order_all_max

      else if(ordering=='bisection') then ! Init bisection balancing

       call build_bisection(update=.false.)

      else ! Init ksection balancing

       call build_ksection(update=.false.)

  end if bisection_or_ordering

  ! Compute coarse cpu map
  do iz=kcoarse_min,kcoarse_max
  do iy=jcoarse_min,jcoarse_max
  do ix=icoarse_min,icoarse_max
     ind=1+ix+iy*nx+iz*nxny
     x(1,1)=(dble(ix)+0.5d0-dble(icoarse_min))*scale
#if NDIM>1
     x(1,2)=(dble(iy)+0.5d0-dble(jcoarse_min))*scale
#endif
#if NDIM>2
     x(1,3)=(dble(iz)+0.5d0-dble(kcoarse_min))*scale
#endif
     call cmp_cpumap(x,c,1)
     cpu_map(ind)=c(1)
  end do
  end do
  end do

  ! Allocate linked list for each level
  allocate(headl(1:ncpu,1:nlevelmax))
  allocate(taill(1:ncpu,1:nlevelmax))
  allocate(numbl(1:ncpu,1:nlevelmax))
  allocate(numbtot(1:10,1:nlevelmax))
  headl=0    ! Head grid in the level
  taill=0    ! Tail grid in the level
  numbl=0    ! Number of grids in the level
  numbtot=0  ! Total number of grids in the level

  ! Allocate communicators
  allocate(active(1:nlevelmax))
  allocate(emission(1:ncpu,1:nlevelmax))
  allocate(reception(1:ncpu,1:nlevelmax))
  do ilevel=1,nlevelmax
     active(ilevel)%ngrid=0
     do i=1,ncpu
        emission (i,ilevel)%ngrid=0
        emission (i,ilevel)%npart=0
        reception(i,ilevel)%ngrid=0
        reception(i,ilevel)%npart=0
     end do
  end do
  ! Allocate lookup array for multigrid fine
  if(poisson)then
     allocate(lookup_mg(1:ngridmax))
     lookup_mg=0
  endif

  ! Allocate physical boundary for each level
  allocate(headb   (1:MAXBOUND,1:nlevelmax))
  allocate(tailb   (1:MAXBOUND,1:nlevelmax))
  allocate(numbb   (1:MAXBOUND,1:nlevelmax))
  allocate(boundary(1:MAXBOUND,1:nlevelmax))
  do i=1,MAXBOUND
     do ilevel=1,nlevelmax
        headb   (i,ilevel)=0       ! Head grid in boundary
        tailb   (i,ilevel)=0       ! Tail grid in boundary
        numbb   (i,ilevel)=0       ! Number of grids in boundary
        boundary(i,ilevel)%ngrid=0 ! Communicators
     end do
  end do

  ! Allocate grid center coordinates
  allocate(xg(1:ngridmax,1:ndim))
  xg=0.0D0

  ! Allocate tree arrays
  allocate(father(1:ngridmax))
  allocate(nbor  (1:1,1:1))  ! Minimal — nbor computed from Morton keys
  allocate(next  (1:ngridmax))
  allocate(prev  (1:ngridmax))
  father=0; nbor=0; next=0; prev=0

  ! Allocate pointer to particles linked lists
  if(pic)then
     allocate(headp(1:ngridmax))
     allocate(tailp(1:ngridmax))
     allocate(numbp(1:ngridmax))
     headp=0; tailp=0; numbp=0
  endif

  ! Initialize AMR grid linked list
  do i=1,ngridmax-1
     next(i)=i+1
  end do
  do i=2,ngridmax
     prev(i)=i-1
  end do
  headf=1          ! Pointer to first grid in free memory
  tailf=ngridmax   ! Pointer to last grid in free memory
  prev(headf)=0; next(tailf)=0
  numbf=ngridmax   ! Number of grids in free memory
  used_mem=ngridmax-numbf

  !----------------------------
  ! Read amr file for a restart
  !----------------------------
  if(nrestart>0)then
#ifdef HDF5
     if(informat == 'hdf5') then
        call restore_amr_hdf5()
        if(myid==1) write(*,*) 'HDF5 AMR backup files read completed'
        goto 998
     end if
#endif
     ! Early varcpu detection: rank 1 probes file 00001 for ncpu
#ifndef WITHOUTMPI
     call title(nrestart,nchar)
     if(myid==1) then
        if(IOGROUPSIZEREP>0)then
           call title(1,ncharcpu)
           fileloc='output_'//TRIM(nchar)//'/group_'//TRIM(ncharcpu)//'/amr_'//TRIM(nchar)//'.out'
        else
           fileloc='output_'//TRIM(nchar)//'/amr_'//TRIM(nchar)//'.out'
        endif
        call title(1,ncharcpu)
        fileloc=TRIM(fileloc)//TRIM(ncharcpu)
        open(unit=99,file=fileloc,form='unformatted')
        read(99) ncpu2
        read(99) ndim2
        read(99) nx2,ny2,nz2
        read(99) nlevelmax2
        close(99)
     end if
     call MPI_BCAST(ncpu2,1,MPI_INTEGER,0,MPI_COMM_WORLD,info2)
     call MPI_BCAST(nlevelmax2,1,MPI_INTEGER,0,MPI_COMM_WORLD,info2)
     if(ncpu2.ne.ncpu)then
        if(ordering /= 'ksection') then
           if(myid==1) write(*,*)'Variable-ncpu restart only supported with ksection ordering'
           call clean_stop
        end if
        call restore_amr_binary_varcpu(ncpu2,nlevelmax2)
        goto 998
     end if
#endif
     ! Wait for the token
#ifndef WITHOUTMPI
     if(IOGROUPSIZE>0) then
        if (mod(myid-1,IOGROUPSIZE)/=0) then
           call MPI_RECV(dummy_io,1,MPI_INTEGER,myid-1-1,tag,&
                & MPI_COMM_WORLD,MPI_STATUS_IGNORE,info2)
        end if
     endif
#endif

#ifdef QUADHILBERT
    if(nrestart_quad.eq.nrestart) allocate(bound_key_restart(0:ndomain))
#endif

     ilun=myid+10
     call title(nrestart,nchar)

     if(IOGROUPSIZEREP>0)then
        call title(((myid-1)/IOGROUPSIZEREP)+1,ncharcpu)
        fileloc='output_'//TRIM(nchar)//'/group_'//TRIM(ncharcpu)//'/amr_'//TRIM(nchar)//'.out'
     else
        fileloc='output_'//TRIM(nchar)//'/amr_'//TRIM(nchar)//'.out'
     endif

     call title(myid,nchar)
     fileloc=TRIM(fileloc)//TRIM(nchar)
     inquire(file=fileloc, exist=ok)
     if(.not. ok)then
        write(*,*)'Restart failed:'
        write(*,*)'File '//TRIM(fileloc)//' not found'
        call clean_stop
     end if
     if(debug)write(*,*)'amr.tmp opened for processor ',myid
     open(unit=ilun,file=fileloc,form='unformatted')
     ! Read grid variables
     read(ilun)ncpu2
     read(ilun)ndim2
     read(ilun)nx2,ny2,nz2
     read(ilun)nlevelmax2
     read(ilun)ngridmax2
     read(ilun)nboundary2
     read(ilun)ngrid_current
     read(ilun)boxlen

     if(ncpu2.ne.ncpu)then
        ! Should not reach here (caught by early detection above)
        close(ilun)
        call restore_amr_binary_varcpu(ncpu2,nlevelmax2)
        goto 998
     end if
     ! Read time variables
     read(ilun)noutput2,iout2,ifout2
     if(noutput2>MAXOUT)then
       write(*,*) 'Error: noutput>MAXOUT'
       call clean_stop
     end if
     read(ilun)tout2(1:noutput2)
     read(ilun)aout2(1:noutput2)
     ! Check compatibility with current parameters
     if((ndim2.ne.ndim).or.(nx2.ne.nx).or.(ny2.ne.ny).or.(nz2.ne.nz).or.&
          & (nboundary2.ne.nboundary).or.(nlevelmax2>nlevelmax).or.&
          & (ngrid_current>ngridmax).or.(noutput2>noutput) )then
        write(*,*)'File amr.tmp is not compatible with namelist'
        write(*,*)'         ndim   nx   ny   nz nlevelmax noutput   ngridmax nboundary'
        write(*,'("amr.tmp  =",4(I4,1x),5x,I4,4x,I4,3x,I8)')&
             & ndim2,nx2,ny2,nz2,nlevelmax2,noutput2,ngrid_current,nboundary2
        write(*,'("namelist =",4(I4,1x),5x,I4,4x,I4,3x,I8)')&
             & ndim ,nx ,ny ,nz ,nlevelmax ,noutput, ngridmax     ,nboundary
        !jhshin1
        !if(myid==1)write(*,*)'Restart failed'
        !call clean_stop 
        !jhshin2
     end if
     !jhshin1
     ! Old output times
     !tout(1:noutput2)=tout2(1:noutput2)
     !aout(1:noutput2)=aout2(1:noutput2)
     !jhshin2
     iout=iout2
     ifout=ifout2
     if(ifout.gt.nrestart+1) ifout=nrestart+1
     read(ilun)t
     read(ilun)dtold(1:nlevelmax2)
     read(ilun)dtnew(1:nlevelmax2)
     read(ilun)nstep,nstep_coarse
     nstep_coarse_old=nstep_coarse
     read(ilun)const,mass_tot_0,rho_tot
     read(ilun)omega_m,omega_l,omega_k,omega_b,h0,aexp_ini,boxlen_ini
     read(ilun)aexp,hexp,aexp_old,epot_tot_int,epot_tot_old
     if(cosmo)then
        read(ilun)mass_sph
     else
        read(ilun)mass_sph2
     endif
     if(myid==1)write(*,*)'Restarting at t=',t,' nstep_coarse=',nstep_coarse

     ! Compute movie frame number if applicable
     if(imovout>0) then
        do i=2,imovout
           if(aendmov>0)then
              if(aexp>amovout(i-1).and.aexp<amovout(i)) then
                 imov=i
              endif
           else
              if(t>tmovout(i-1).and.t<tmovout(i)) then
                 imov=i
              endif
           endif
        enddo
        if(aendmov>0)then
           if(myid==1)write(*,*) "Frame number, aexp ",imov, amovout(imov)
        else
           if(myid==1)write(*,*) "Frame number, t ",imov, tmovout(imov)
        endif
     endif

     ! Read levels variables
     read(ilun)headl(1:ncpu,1:nlevelmax2)
     read(ilun)taill(1:ncpu,1:nlevelmax2)
     read(ilun)numbl(1:ncpu,1:nlevelmax2)
     read(ilun)numbtot(1:10,1:nlevelmax2)
     ! Read boundary linked list
     if(simple_boundary)then
        read(ilun)headb(1:nboundary,1:nlevelmax2)
        read(ilun)tailb(1:nboundary,1:nlevelmax2)
        read(ilun)numbb(1:nboundary,1:nlevelmax2)
     end if
     ! Read free memory
     read(ilun)headf,tailf,numbf,used_mem,used_mem_tot
     headf=ngrid_current+1
     tailf=ngridmax
     numbf=ngridmax-ngrid_current
     prev(headf)=0
     next(tailf)=0
     ! Read cpu boundaries
     read(ilun)ordering2
     if(ordering2.ne.ordering)then
        if(myid==1)write(*,*)'Ordering is uncompatible'
        call clean_stop
     endif
     if(ordering=='bisection') then
        read(ilun)bisec_wall(1:nbinodes)
        read(ilun)bisec_next(1:nbinodes,1:2)
        read(ilun)bisec_indx(1:nbinodes)
        read(ilun)bisec_cpubox_min(1:ncpu,1:ndim)
        read(ilun)bisec_cpubox_max(1:ncpu,1:ndim)
     else if(ordering=='ksection') then
        read(ilun)nksec_levels
        read(ilun)ksec_kmax
        read(ilun)ksec_nbinodes
        read(ilun)ksec_factor(1:nksec_levels)
        read(ilun)ksec_dir(1:nksec_levels)
        read(ilun)ksec_wall(1:ksec_nbinodes,1:max(ksec_kmax-1,1))
        read(ilun)ksec_next(1:ksec_nbinodes,1:ksec_kmax)
        read(ilun)ksec_indx(1:ksec_nbinodes)
        read(ilun)bisec_cpubox_min(1:ncpu,1:ndim)
        read(ilun)bisec_cpubox_max(1:ncpu,1:ndim)
        ! Rebuild tree navigation arrays from restored tree
        call rebuild_ksec_cpuranges()
        call compute_ksec_cpu_path()
     else
#ifdef QUADHILBERT
        if(nrestart_quad.eq.nrestart) then
           read(ilun)bound_key_restart(0:ndomain)
           bound_key(0:ndomain)=bound_key_restart(0:ndomain)
        else
           read(ilun)bound_key(0:ndomain)
        endif
#else
        read(ilun)bound_key(0:ndomain)
#endif
     endif
     ! Read coarse level
     read(ilun)son(1:ncoarse)
     read(ilun)flag1(1:ncoarse)
     read(ilun)cpu_map(1:ncoarse)
     ! Read fine levels
     do ilevel=1,nlevelmax2
        do ibound=1,nboundary+ncpu
           if(ibound<=ncpu)then
              ncache=numbl(ibound,ilevel)
           else
              ncache=numbb(ibound-ncpu,ilevel)
           end if
           if(ncache>0)then
              allocate(ind_grid(1:ncache))
              allocate(xxg(1:ncache))
              allocate(iig(1:ncache))
              allocate(pos(1:ncache))
              allocate(grid(1:ncache))
              ! Read grid index
              read(ilun)ind_grid
              ! Read next index
              read(ilun)iig
              do i=1,ncache
                 next(ind_grid(i))=iig(i)
              end do
              ! Read prev index
              read(ilun)iig
              do i=1,ncache
                 prev(ind_grid(i))=iig(i)
              end do              
              ! Read grid center
              do idim=1,ndim
                 read(ilun)xxg
                 do i=1,ncache
                    xg(ind_grid(i),idim)=xxg(i)
                 end do
              end do
              ! Read father index
              read(ilun)iig
              if(ngridmax.ne.ngridmax2.and.ilevel>1)then
                 do i=1,ncache
                    pos(i)=(iig(i)-ncoarse-1)/ngridmax2
                 end do
                 do i=1,ncache
                    grid(i)=iig(i)-ncoarse-pos(i)*ngridmax2
                 end do
                 do i=1,ncache
                    iig(i)=ncoarse+pos(i)*ngridmax+grid(i)
                 end do
              end if
              do i=1,ncache
                 father(ind_grid(i))=iig(i)
              end do
              ! Read nbor index (read but ignore — computed from Morton keys)
              do ind=1,twondim
                 read(ilun)iig
              end do
              ! Read son index
              do ind=1,twotondim
                 iskip=ncoarse+(ind-1)*ngridmax
                 read(ilun)iig
                 do i=1,ncache
                    son(ind_grid(i)+iskip)=iig(i)
                 end do
              end do
              ! Read cpu map
              do ind=1,twotondim
                 iskip=ncoarse+(ind-1)*ngridmax
                 read(ilun)iig
                 do i=1,ncache
                    cpu_map(ind_grid(i)+iskip)=iig(i)
                 end do
              end do
              ! Read refinement map
              do ind=1,twotondim
                 iskip=ncoarse+(ind-1)*ngridmax
                 read(ilun)iig
                 do i=1,ncache
                    flag1(ind_grid(i)+iskip)=iig(i)
                 end do
              end do
              deallocate(xxg,iig,pos,grid,ind_grid)
           end if
        end do
     end do
     close(ilun)
#ifdef QUADHILBERT
     if(nrestart_quad.eq.nrestart) deallocate(bound_key_restart)
#endif

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



#ifndef WITHOUTMPI
     if(debug)write(*,*)'amr.tmp read for processor ',myid
     call MPI_BARRIER(MPI_COMM_WORLD,info)
#endif
     if(verbose)write(*,*)'AMR backup files read completed'

     ! Build communicators
     do ilevel=1,nlevelmax
        call build_comm(ilevel)
     end do

  end if

998 continue  ! HDF5 restore jumps here

  ! Test ksection exchange subroutines
  if(ordering=='ksection') call test_ksection_exchange()

  ! Build Morton hash tables
  call morton_hash_build_and_verify()

end subroutine init_amr


!###########################################################################
! Binary variable-ncpu AMR restore (distributed I/O version)
! Called when ncpu_file != ncpu and ordering=='ksection'
! Each rank reads only its assigned files, then MPI_ALLTOALLV exchanges
! grid data to the correct owners level by level.
! Exchange metadata is stored for hydro/poisson reuse.
!###########################################################################
subroutine restore_amr_binary_varcpu(ncpu2_in, nlevelmax2_in)
  use amr_commons
  use hydro_commons
  use pm_commons
  use poisson_commons
  use ksection
  use morton_keys
  use morton_hash
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer, intent(in) :: ncpu2_in, nlevelmax2_in

  integer :: icpu_file, ilevel, ibound, i, j, idim, ind, info
  integer :: ncache, ilun, nrec_skip, ifile
  integer :: nlevelmax_file, ngridmax_file, nboundary_file
  integer, allocatable :: numbl_file(:,:)  ! (ncpu_file, nlevelmax_file) per-file numbl
  integer, allocatable :: numbb_file(:,:)  ! (nboundary, nlevelmax_file)
  integer, allocatable :: nactive_local(:,:)

  ! Per-level local grid buffers (only grids from assigned files)
  type local_level_t
     real(dp), allocatable :: xg(:,:)        ! (ngrids, ndim)
     integer, allocatable :: son_flag(:)     ! (ngrids * twotondim)
     integer :: ngrids = 0
  end type
  type(local_level_t) :: my_lvl(1:MAXLEVEL)
  integer :: lvl_offset(1:MAXLEVEL)

  ! File reading buffers
  integer, allocatable :: iig(:)
  real(dp), allocatable :: xxg(:)

  ! MPI exchange variables
  integer :: nlocal, nsend, nrecv, npack, d
  integer, allocatable :: dest(:), cursor(:)
  real(dp), allocatable :: sendbuf(:), recvbuf(:)
  integer, allocatable :: scount(:), rcount(:), sdispl(:), rdispl(:)
  integer, allocatable :: scount_p(:), rcount_p(:), sdispl_p(:), rdispl_p(:)
  integer :: ipacked

  ! Grid creation variables
  integer :: igrid_new, igrid_prev_cpu, igrid_father, ind_cell
  integer :: ix, iy, iz, ix_p, iy_p, iz_p, nxny, nx_loc, iskip
  integer(mkb) :: mkey
  real(dp) :: twotol, dx, scale, xx_cell(1,3), xx_father(1,3), xg_recv(3)
  integer :: c_tmp(1)
  real(dp) :: mass_sph2

  character(LEN=128) :: ordering_file
  character(LEN=80) :: fileloc
  character(LEN=5) :: nchar, ncharcpu

  ! ============================================================
  ! Step 1: Setup
  ! ============================================================
  ncpu_file = ncpu2_in
  nlevelmax_file = nlevelmax2_in
  varcpu_restart = .true.
  varcpu_restart_done = .true.

  if(myid==1) then
     write(*,*) '============================================================'
     write(*,*) '=====  VARIABLE NCPU BINARY RESTART (distributed I/O)  ====='
     write(*,*) '=====  File ncpu = ', ncpu_file
     write(*,*) '=====  New  ncpu = ', ncpu
     write(*,*) '============================================================'
  end if

  allocate(varcpu_ngrid_file(1:nlevelmax))
  varcpu_ngrid_file = 0

  nxny = nx * ny
  nx_loc = icoarse_max - icoarse_min + 1
  scale = boxlen / dble(nx_loc)

  ! Assign files to ranks (round-robin)
  varcpu_nfiles_local = 0
  do icpu_file = 1, ncpu_file
     if(mod(icpu_file-1, ncpu) + 1 == myid) &
          varcpu_nfiles_local = varcpu_nfiles_local + 1
  end do
  allocate(varcpu_my_files(max(varcpu_nfiles_local,1)))
  j = 0
  do icpu_file = 1, ncpu_file
     if(mod(icpu_file-1, ncpu) + 1 == myid) then
        j = j + 1
        varcpu_my_files(j) = icpu_file
     end if
  end do

  ! ============================================================
  ! Step 2: Read header from file 00001 (all ranks read — small file)
  ! ============================================================
  ilun = 99
  call title(nrestart, nchar)
  if(IOGROUPSIZEREP>0)then
     call title(1, ncharcpu)
     fileloc='output_'//TRIM(nchar)//'/group_'//TRIM(ncharcpu)//'/amr_'//TRIM(nchar)//'.out'
  else
     fileloc='output_'//TRIM(nchar)//'/amr_'//TRIM(nchar)//'.out'
  endif
  call title(1, ncharcpu)
  fileloc=TRIM(fileloc)//TRIM(ncharcpu)

  open(unit=ilun, file=fileloc, form='unformatted')
  read(ilun)  ! ncpu
  read(ilun)  ! ndim
  read(ilun)  ! nx,ny,nz
  read(ilun)  ! nlevelmax
  read(ilun) ngridmax_file
  read(ilun) nboundary_file
  read(ilun) ngrid_current
  read(ilun)  ! boxlen

  ! Time variables
  block
    integer :: nout_tmp, iout_tmp, ifout_tmp
    real(dp), allocatable :: dp_tmp(:)
    read(ilun) nout_tmp, iout_tmp, ifout_tmp
    iout = iout_tmp; ifout = ifout_tmp
    if(ifout.gt.nrestart+1) ifout=nrestart+1
    allocate(dp_tmp(max(nout_tmp,1)))
    read(ilun) dp_tmp(1:nout_tmp)  ! tout
    read(ilun) dp_tmp(1:nout_tmp)  ! aout
    deallocate(dp_tmp)
  end block
  read(ilun) t
  read(ilun) dtold(1:nlevelmax_file)
  read(ilun) dtnew(1:nlevelmax_file)
  read(ilun) nstep, nstep_coarse
  nstep_coarse_old = nstep_coarse
  read(ilun) const, mass_tot_0, rho_tot
  read(ilun) omega_m, omega_l, omega_k, omega_b, h0, aexp_ini, boxlen_ini
  read(ilun) aexp, hexp, aexp_old, epot_tot_int, epot_tot_old
  if(cosmo)then
     read(ilun) mass_sph
  else
     read(ilun) mass_sph2
  endif
  if(myid==1) write(*,*) 'Restarting at t=', t, ' nstep_coarse=', nstep_coarse

  ! Level linked lists (ncpu_file × nlevelmax_file) — just for file 00001 header
  allocate(numbl_file(1:ncpu_file, 1:nlevelmax_file))
  read(ilun)  ! headl
  read(ilun)  ! taill
  read(ilun) numbl_file
  read(ilun)  ! numbtot

  ! Boundary linked lists
  if(simple_boundary)then
     allocate(numbb_file(1:nboundary_file, 1:nlevelmax_file))
     read(ilun)  ! headb
     read(ilun)  ! tailb
     read(ilun) numbb_file
  else
     allocate(numbb_file(1:1, 1:1))
     numbb_file = 0
  end if

  ! Free memory — skip
  read(ilun)

  ! Ordering — skip tree records (will rebuild)
  read(ilun) ordering_file
  if(trim(ordering_file) == 'ksection') then
     do i = 1, 10; read(ilun); end do
  else if(trim(ordering_file) == 'bisection') then
     do i = 1, 5; read(ilun); end do
  else
     read(ilun)
  end if

  ! Coarse grid
  read(ilun) son(1:ncoarse)
  read(ilun) flag1(1:ncoarse)
  read(ilun) cpu_map(1:ncoarse)
  close(ilun)

  son(1:ncoarse) = 0  ! Reset (set when level-1 grids are created)

  ! ============================================================
  ! Step 3: Rebuild domain decomposition for new ncpu
  ! ============================================================
  call build_ksection(update=.false.)
  call rebuild_ksec_cpuranges()
  call compute_ksec_cpu_path()
  bisec_cpubox_min2 = bisec_cpubox_min
  bisec_cpubox_max2 = bisec_cpubox_max

  ! Recompute coarse cpu_map
  do iz = kcoarse_min, kcoarse_max
  do iy = jcoarse_min, jcoarse_max
  do ix = icoarse_min, icoarse_max
     ind = 1 + ix + iy * nx + iz * nxny
     xx_cell(1,1) = (dble(ix) + 0.5d0 - dble(icoarse_min)) * scale
     xx_cell(1,2) = (dble(iy) + 0.5d0 - dble(jcoarse_min)) * scale
     xx_cell(1,3) = (dble(iz) + 0.5d0 - dble(kcoarse_min)) * scale
     call cmp_ksection_cpumap(xx_cell, c_tmp, 1)
     cpu_map(ind) = c_tmp(1)
  end do
  end do
  end do
  cpu_map2(1:ncoarse) = cpu_map(1:ncoarse)
  if(myid==1) write(*,*) 'Binary varcpu: ksection tree rebuilt, coarse cpu_map recomputed'

  ! ============================================================
  ! Step 4: Initialize Morton hash tables
  ! ============================================================
  if(.not. allocated(mort_table)) allocate(mort_table(1:nlevelmax))
  if(.not. allocated(grid_level)) then
     allocate(grid_level(1:ngridmax))
     grid_level = 0
  end if
  do ilevel = 1, nlevelmax
     call morton_hash_init(mort_table(ilevel), 16)
  end do

  ! ============================================================
  ! Step 5: Read numbl from assigned files, ALLREDUCE for totals
  ! ============================================================
  allocate(nactive_local(1:ncpu_file, 1:nlevelmax_file))
  allocate(varcpu_nactive(1:ncpu_file, 1:nlevelmax_file))
  nactive_local = 0

  do ifile = 1, varcpu_nfiles_local
     icpu_file = varcpu_my_files(ifile)
     if(IOGROUPSIZEREP>0)then
        call title(((icpu_file-1)/IOGROUPSIZEREP)+1, ncharcpu)
        fileloc='output_'//TRIM(nchar)//'/group_'//TRIM(ncharcpu)//'/amr_'//TRIM(nchar)//'.out'
     else
        fileloc='output_'//TRIM(nchar)//'/amr_'//TRIM(nchar)//'.out'
     endif
     call title(icpu_file, ncharcpu)
     fileloc=TRIM(fileloc)//TRIM(ncharcpu)
     open(unit=ilun, file=fileloc, form='unformatted')
     do i = 1, 8; read(ilun); end do
     block
       integer :: nout_tmp, dum1, dum2
       read(ilun) nout_tmp, dum1, dum2
       read(ilun); read(ilun)
     end block
     do i = 1, 8; read(ilun); end do
     read(ilun); read(ilun)
     read(ilun) numbl_file
     close(ilun)
     nactive_local(icpu_file, 1:nlevelmax_file) = numbl_file(icpu_file, 1:nlevelmax_file)
  end do

  call MPI_ALLREDUCE(nactive_local, varcpu_nactive, ncpu_file*nlevelmax_file, &
       MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, info)
  deallocate(nactive_local)

  ! Compute per-level totals
  do ilevel = 1, nlevelmax_file
     varcpu_ngrid_file(ilevel) = sum(varcpu_nactive(:, ilevel))
  end do

  if(myid==1) then
     do ilevel = 1, nlevelmax_file
        if(varcpu_ngrid_file(ilevel) > 0) then
           write(*,'(A,I3,A,I10)') &
                ' Binary varcpu level ', ilevel, ' total grids: ', varcpu_ngrid_file(ilevel)
        end if
     end do
  end if

  ! ============================================================
  ! Step 6: Read active grids from assigned files only
  ! ============================================================
  ! Compute local grid counts per level
  do ilevel = 1, nlevelmax_file
     nlocal = 0
     do ifile = 1, varcpu_nfiles_local
        nlocal = nlocal + varcpu_nactive(varcpu_my_files(ifile), ilevel)
     end do
     my_lvl(ilevel)%ngrids = nlocal
     if(nlocal > 0) then
        allocate(my_lvl(ilevel)%xg(nlocal, ndim))
        allocate(my_lvl(ilevel)%son_flag(nlocal * twotondim))
     end if
  end do

  lvl_offset = 0
  do ifile = 1, varcpu_nfiles_local
     icpu_file = varcpu_my_files(ifile)
     if(IOGROUPSIZEREP>0)then
        call title(((icpu_file-1)/IOGROUPSIZEREP)+1, ncharcpu)
        fileloc='output_'//TRIM(nchar)//'/group_'//TRIM(ncharcpu)//'/amr_'//TRIM(nchar)//'.out'
     else
        fileloc='output_'//TRIM(nchar)//'/amr_'//TRIM(nchar)//'.out'
     endif
     call title(icpu_file, ncharcpu)
     fileloc=TRIM(fileloc)//TRIM(ncharcpu)
     open(unit=ilun, file=fileloc, form='unformatted')

     ! Skip 8 header records
     do i = 1, 8; read(ilun); end do
     ! Skip noutput/tout/aout
     block
       integer :: nout_tmp, dum1, dum2
       read(ilun) nout_tmp, dum1, dum2
       read(ilun); read(ilun)
     end block
     ! Skip 8 time records
     do i = 1, 8; read(ilun); end do
     ! Read THIS file's numbl
     read(ilun); read(ilun)
     read(ilun) numbl_file
     read(ilun)  ! numbtot
     if(simple_boundary) then
        do i = 1, 3; read(ilun); end do
     end if
     read(ilun)  ! free memory
     block
       character(LEN=128) :: ord_tmp
       read(ilun) ord_tmp
       if(trim(ord_tmp) == 'ksection') then
          do i = 1, 10; read(ilun); end do
       else if(trim(ord_tmp) == 'bisection') then
          do i = 1, 5; read(ilun); end do
       else
          read(ilun)
       end if
     end block
     do i = 1, 3; read(ilun); end do  ! coarse

     ! Read fine levels — extract active blocks (ibound == icpu_file)
     do ilevel = 1, nlevelmax_file
        do ibound = 1, nboundary_file + ncpu_file
           if(ibound <= ncpu_file) then
              ncache = numbl_file(ibound, ilevel)
           else
              ncache = numbb_file(ibound - ncpu_file, ilevel)
           end if

           if(ncache > 0) then
              if(ibound == icpu_file) then
                 allocate(xxg(1:ncache), iig(1:ncache))
                 ind = lvl_offset(ilevel)

                 read(ilun); read(ilun); read(ilun)  ! ind_grid, next, prev
                 do idim = 1, ndim
                    read(ilun) xxg
                    my_lvl(ilevel)%xg(ind+1:ind+ncache, idim) = xxg(1:ncache)
                 end do
                 read(ilun)  ! father
                 do i = 1, twondim; read(ilun); end do  ! nbor
                 do iskip = 1, twotondim
                    read(ilun) iig
                    do i = 1, ncache
                       if(iig(i) > 0) then
                          my_lvl(ilevel)%son_flag((ind + i - 1)*twotondim + iskip) = 1
                       else
                          my_lvl(ilevel)%son_flag((ind + i - 1)*twotondim + iskip) = 0
                       end if
                    end do
                 end do
                 do i = 1, 2*twotondim; read(ilun); end do  ! cpu_map + flag1

                 lvl_offset(ilevel) = ind + ncache
                 deallocate(xxg, iig)
              else
                 nrec_skip = 4 + ndim + twondim + 3*twotondim
                 do i = 1, nrec_skip; read(ilun); end do
              end if
           end if
        end do
     end do
     close(ilun)
  end do

  ! ============================================================
  ! Step 7: Level-by-level MPI exchange + grid creation
  ! ============================================================
  balance = .true.
  shrink = .false.
  npack = 4  ! xg(1:3) + packed_son_flag

  allocate(varcpu_exc(1:nlevelmax))

  do ilevel = 1, nlevelmax_file
     if(varcpu_ngrid_file(ilevel) == 0) cycle

     nlocal = my_lvl(ilevel)%ngrids

     ! Resize Morton hash table for this level
     call morton_hash_init(mort_table(ilevel), &
          max(4 * (varcpu_ngrid_file(ilevel) / max(ncpu, 1) + 1), 16))

     dx = 0.5d0**ilevel

     ! 7a. Determine owner for each local grid via cmp_ksection_cpumap
     allocate(dest(max(nlocal,1)))
     do i = 1, nlocal
        twotol = 2.0d0**(ilevel-1)
        ix = int(my_lvl(ilevel)%xg(i, 1) * twotol)
        iy = int(my_lvl(ilevel)%xg(i, 2) * twotol)
        iz = int(my_lvl(ilevel)%xg(i, 3) * twotol)
        xx_father(1,1) = ((dble(ix) + 0.5d0) / twotol - dble(icoarse_min)) * scale
        xx_father(1,2) = ((dble(iy) + 0.5d0) / twotol - dble(jcoarse_min)) * scale
        xx_father(1,3) = ((dble(iz) + 0.5d0) / twotol - dble(kcoarse_min)) * scale
        call cmp_ksection_cpumap(xx_father, c_tmp, 1)
        dest(i) = c_tmp(1)
     end do

     ! 7b. Count sends per destination
     allocate(scount(0:ncpu-1), rcount(0:ncpu-1))
     scount = 0
     do i = 1, nlocal
        scount(dest(i) - 1) = scount(dest(i) - 1) + 1
     end do

     ! 7c. Exchange counts
     call MPI_ALLTOALL(scount, 1, MPI_INTEGER, rcount, 1, MPI_INTEGER, &
          MPI_COMM_WORLD, info)
     nsend = sum(scount)
     nrecv = sum(rcount)

     ! 7d. Compute displacements
     allocate(sdispl(0:ncpu-1), rdispl(0:ncpu-1))
     sdispl(0) = 0; rdispl(0) = 0
     do i = 1, ncpu-1
        sdispl(i) = sdispl(i-1) + scount(i-1)
        rdispl(i) = rdispl(i-1) + rcount(i-1)
     end do

     ! 7e. Pack send buffer (npack dp per grid) and record send_order
     allocate(sendbuf(npack * max(nsend,1)))
     allocate(cursor(0:ncpu-1))
     cursor = 0

     ! Store exchange metadata
     allocate(varcpu_exc(ilevel)%scount(0:ncpu-1))
     allocate(varcpu_exc(ilevel)%rcount(0:ncpu-1))
     allocate(varcpu_exc(ilevel)%sdispl(0:ncpu-1))
     allocate(varcpu_exc(ilevel)%rdispl(0:ncpu-1))
     allocate(varcpu_exc(ilevel)%send_order(max(nsend,1)))
     varcpu_exc(ilevel)%scount = scount
     varcpu_exc(ilevel)%rcount = rcount
     varcpu_exc(ilevel)%sdispl = sdispl
     varcpu_exc(ilevel)%rdispl = rdispl
     varcpu_exc(ilevel)%nsend = nsend
     varcpu_exc(ilevel)%nrecv = nrecv

     do i = 1, nlocal
        d = dest(i) - 1
        j = sdispl(d) + cursor(d)  ! 0-based index in send buffer
        sendbuf(npack*j + 1) = my_lvl(ilevel)%xg(i, 1)
        sendbuf(npack*j + 2) = my_lvl(ilevel)%xg(i, 2)
        sendbuf(npack*j + 3) = my_lvl(ilevel)%xg(i, 3)
        ! Pack twotondim son_flag bits into one integer
        ipacked = 0
        do iskip = 1, twotondim
           if(my_lvl(ilevel)%son_flag((i-1)*twotondim + iskip) > 0) &
                ipacked = ipacked + 2**(iskip-1)
        end do
        sendbuf(npack*j + 4) = dble(ipacked)
        varcpu_exc(ilevel)%send_order(j + 1) = i  ! 1-based
        cursor(d) = cursor(d) + 1
     end do
     deallocate(cursor, dest)

     ! 7f. MPI_ALLTOALLV exchange
     allocate(recvbuf(npack * max(nrecv,1)))
     allocate(scount_p(0:ncpu-1), rcount_p(0:ncpu-1))
     allocate(sdispl_p(0:ncpu-1), rdispl_p(0:ncpu-1))
     scount_p = scount * npack
     rcount_p = rcount * npack
     sdispl_p = sdispl * npack
     rdispl_p = rdispl * npack
     call MPI_ALLTOALLV(sendbuf, scount_p, sdispl_p, MPI_DOUBLE_PRECISION, &
          recvbuf, rcount_p, rdispl_p, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, info)
     deallocate(sendbuf, scount_p, rcount_p, sdispl_p, rdispl_p)
     deallocate(scount, rcount, sdispl, rdispl)

     ! 7g. Create grids from received data (all owned by myid)
     allocate(varcpu_exc(ilevel)%recv_igrid(max(nrecv,1)))
     varcpu_exc(ilevel)%recv_igrid = 0
     igrid_prev_cpu = 0

     do i = 1, nrecv
        xg_recv(1) = recvbuf(npack*(i-1) + 1)
        xg_recv(2) = recvbuf(npack*(i-1) + 2)
        xg_recv(3) = recvbuf(npack*(i-1) + 3)
        ipacked = int(recvbuf(npack*(i-1) + 4))

        ! Find father cell
        if(ilevel == 1) then
           ix = int(xg_recv(1))
           iy = int(xg_recv(2))
           iz = int(xg_recv(3))
           igrid_father = 1 + ix + iy * nx + iz * nxny
        else
           twotol = 2.0d0**(ilevel-1)
           ix = int(xg_recv(1) * twotol)
           iy = int(xg_recv(2) * twotol)
           iz = int(xg_recv(3) * twotol)
           ix_p = ix / 2
           iy_p = iy / 2
           iz_p = iz / 2
           mkey = morton_encode(ix_p, iy_p, iz_p)
           igrid_father = morton_hash_lookup(mort_table(ilevel-1), mkey)
           if(igrid_father == 0) then
              varcpu_exc(ilevel)%recv_igrid(i) = 0
              cycle  ! Father not in hash (virtual); handled by refine_fine
           end if
           ind_cell = 1 + mod(ix, 2) + 2 * mod(iy, 2) + 4 * mod(iz, 2)
           igrid_father = ncoarse + (ind_cell - 1) * ngridmax + igrid_father
        end if

        ! Allocate grid from free list
        igrid_new = headf
        if(igrid_new == 0) then
           write(*,*) 'ERROR: out of free grids, myid=', myid, &
                ' level=', ilevel, ' recv=', i
           call clean_stop
        end if
        headf = next(igrid_new)
        if(headf > 0) prev(headf) = 0
        numbf = numbf - 1

        ! Set xg
        xg(igrid_new, 1:ndim) = xg_recv(1:ndim)

        ! Unpack son_flag → flag1
        do iskip = 1, twotondim
           ind_cell = ncoarse + (iskip - 1) * ngridmax + igrid_new
           son(ind_cell) = 0
           if(btest(ipacked, iskip-1)) then
              flag1(ind_cell) = 1
           else
              flag1(ind_cell) = 0
           end if
        end do

        ! Set cpu_map for each child cell
        do iskip = 1, twotondim
           ind_cell = ncoarse + (iskip - 1) * ngridmax + igrid_new
           iz = (iskip - 1) / 4
           iy = (iskip - 1 - 4 * iz) / 2
           ix = (iskip - 1 - 2 * iy - 4 * iz)
           xx_cell(1,1) = (xg_recv(1) + (dble(ix) - 0.5d0) * dx &
                - dble(icoarse_min)) * scale
           xx_cell(1,2) = (xg_recv(2) + (dble(iy) - 0.5d0) * dx &
                - dble(jcoarse_min)) * scale
           xx_cell(1,3) = (xg_recv(3) + (dble(iz) - 0.5d0) * dx &
                - dble(kcoarse_min)) * scale
           call cmp_ksection_cpumap(xx_cell, c_tmp, 1)
           cpu_map(ind_cell) = c_tmp(1)
           cpu_map2(ind_cell) = c_tmp(1)
        end do

        ! Set father + son in parent
        father(igrid_new) = igrid_father
        son(igrid_father) = igrid_new

        ! Insert into Morton hash
        mkey = grid_to_morton(igrid_new, ilevel)
        call morton_hash_insert(mort_table(ilevel), mkey, igrid_new)
        grid_level(igrid_new) = ilevel

        ! Append to linked list for myid
        if(igrid_prev_cpu == 0) then
           headl(myid, ilevel) = igrid_new
        else
           next(igrid_prev_cpu) = igrid_new
        end if
        prev(igrid_new) = igrid_prev_cpu
        next(igrid_new) = 0
        taill(myid, ilevel) = igrid_new
        numbl(myid, ilevel) = numbl(myid, ilevel) + 1
        igrid_prev_cpu = igrid_new

        varcpu_exc(ilevel)%recv_igrid(i) = igrid_new
     end do
     deallocate(recvbuf)

     ! Free level data
     if(allocated(my_lvl(ilevel)%xg)) deallocate(my_lvl(ilevel)%xg)
     if(allocated(my_lvl(ilevel)%son_flag)) deallocate(my_lvl(ilevel)%son_flag)

     ! Create virtual grids via RAMSES refine mechanism
     if(ilevel == 1) then
        call flag_coarse
        call refine_coarse
     else
        call refine_fine(ilevel - 1)
     end if

     call build_comm(ilevel)

     ! Exchange flag1/cpu_map to virtual grids (needed for next level's refine)
     call make_virtual_fine_int(flag1(1), ilevel)
     call make_virtual_fine_int(cpu_map(1), ilevel)
     call make_virtual_fine_int(cpu_map2(1), ilevel)

     if(myid==1) write(*,'(A,I3,A,I10)') &
          ' Binary varcpu level ', ilevel, ' active local: ', numbl(myid, ilevel)
  end do

  balance = .false.
  shrink = .false.

  deallocate(numbl_file, numbb_file)

  if(myid==1) write(*,*) 'Binary variable-ncpu AMR restore done.'

end subroutine restore_amr_binary_varcpu

