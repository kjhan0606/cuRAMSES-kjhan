subroutine init_poisson
  use pm_commons
  use amr_commons
  use poisson_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer::ncell,ncache,iskip,igrid,i,ilevel,ind,ivar
  integer::nvar2,ilevel2,numbl2,ilun,ibound,istart,info
  integer::ncpu2,ndim2,nlevelmax2,nboundary2
  integer ,dimension(:),allocatable::ind_grid
  real(dp),dimension(:),allocatable::xx
  character(LEN=80)::fileloc
  character(LEN=5)::nchar,ncharcpu
  integer,parameter::tag=1114
  integer::dummy_io,info2

  if(verbose)write(*,*)'Entering init_poisson'

  !------------------------------------------------------
  ! Allocate cell centered variables arrays
  !------------------------------------------------------
  ncell=ncoarse+twotondim*ngridmax
  allocate(rho (1:ncell))
  allocate(rho_star (1:ncell))
  allocate(phi (1:ncell))
  allocate(phi_old (1:ncell))
  allocate(f   (1:ncell,1:3))
  rho=0.0D0; rho_star=0d0; phi=0.0D0; f=0.0D0
  if(cic_levelmax>0)then
     allocate(rho_top(1:ncell))
     rho_top=0d0
  endif

  !------------------------------------------------------
  ! Allocate multigrid variables
  !------------------------------------------------------
  ! Allocate communicators for coarser multigrid levels
  allocate(active_mg    (1:ncpu,1:nlevelmax-1))
  allocate(emission_mg  (1:ncpu,1:nlevelmax-1))
  do ilevel=1,nlevelmax-1
     do i=1,ncpu
        active_mg   (i,ilevel)%ngrid=0
        active_mg   (i,ilevel)%npart=0
        emission_mg (i,ilevel)%ngrid=0
        emission_mg (i,ilevel)%npart=0
     end do
  end do
  allocate(safe_mode(1:nlevelmax))
  safe_mode = .false.

  !--------------------------------
  ! For a restart, read poisson file
  !--------------------------------
  if(nrestart>0)then
#ifdef HDF5
     if(informat == 'hdf5') then
        call restore_poisson_hdf5()
        if(verbose)write(*,*)'HDF5 POISSON backup files read completed'
        return
     end if
#endif
     if(varcpu_restart) then
        call restore_poisson_binary_varcpu()
#ifndef WITHOUTMPI
        call MPI_BARRIER(MPI_COMM_WORLD,info)
#endif
        if(verbose)write(*,*)'Binary varcpu POISSON backup files read completed'
        return
     end if
     ilun=ncpu+myid+10
     call title(nrestart,nchar)
     if(IOGROUPSIZEREP>0)then
        call title(((myid-1)/IOGROUPSIZEREP)+1,ncharcpu)
        fileloc='output_'//TRIM(nchar)//'/group_'//TRIM(ncharcpu)//'/grav_'//TRIM(nchar)//'.out'
     else
        fileloc='output_'//TRIM(nchar)//'/grav_'//TRIM(nchar)//'.out'
     endif
     call title(myid,nchar)
     fileloc=TRIM(fileloc)//TRIM(nchar)

     ! Wait for the token
#ifndef WITHOUTMPI
     if(IOGROUPSIZE>0) then
        if (mod(myid-1,IOGROUPSIZE)/=0) then
           call MPI_RECV(dummy_io,1,MPI_INTEGER,myid-1-1,tag,&
                & MPI_COMM_WORLD,MPI_STATUS_IGNORE,info2)
        end if
     endif
#endif

     open(unit=ilun,file=fileloc,form='unformatted')
     read(ilun)ncpu2
     read(ilun)ndim2
     read(ilun)nlevelmax2
     read(ilun)nboundary2
     if(ndim2.ne.ndim+1)then
        write(*,*)'File poisson.tmp is not compatible'
        write(*,*)'Found   =',ndim2
        write(*,*)'Expected=',ndim+1
        call clean_stop
     end if
     do ilevel=1,nlevelmax2
        do ibound=1,nboundary+ncpu
           if(ibound<=ncpu)then
              ncache=numbl(ibound,ilevel)
              istart=headl(ibound,ilevel)
           else
              ncache=numbb(ibound-ncpu,ilevel)
              istart=headb(ibound-ncpu,ilevel)
           end if
           read(ilun)ilevel2
           read(ilun)numbl2
           if(numbl2.ne.ncache)then
              write(*,*)'File poisson.tmp is not compatible'
              write(*,*)'Found   =',numbl2,' for level ',ilevel2
              write(*,*)'Expected=',ncache,' for level ',ilevel
           end if
           if(ncache>0)then
              allocate(ind_grid(1:ncache))
              allocate(xx(1:ncache))
              ! Loop over level grids
              igrid=istart
              do i=1,ncache
                 ind_grid(i)=igrid
                 igrid=next(igrid)
              end do
              ! Loop over cells
              do ind=1,twotondim
                 iskip=ncoarse+(ind-1)*ngridmax
                 ! Read potential
                 read(ilun)xx
                 do i=1,ncache
                    phi(ind_grid(i)+iskip)=xx(i)
                 end do
                 ! Read force
                 do ivar=1,ndim
                    read(ilun)xx
                    do i=1,ncache
                       f(ind_grid(i)+iskip,ivar)=xx(i)
                    end do
                 end do
              end do
              deallocate(ind_grid,xx)
           end if
        end do
     end do
     close(ilun)

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
     if(debug)write(*,*)'poisson.tmp read for processor ',myid
     call MPI_BARRIER(MPI_COMM_WORLD,info)
#endif
     if(verbose)write(*,*)'POISSON backup files read completed'
  end if

end subroutine init_poisson

subroutine restore_poisson_binary_varcpu
  !--------------------------------------------------------------
  ! Distributed I/O version: each rank reads its assigned grav files,
  ! then MPI_ALLTOALLV exchanges data using stored AMR exchange metadata.
  ! Frees all exchange metadata after completion.
  !--------------------------------------------------------------
  use amr_commons
  use poisson_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer :: icpu_file, ilevel, ibound, i, j, ind, ivar, info, ifile
  integer :: ncache, ilun, ndim2, nlevelmax2, nboundary2, ncpu2
  integer :: igrid, iskip, icell, nval_per_grid, base
  integer :: nlocal, nsend, nrecv, ngrav_per_oct
  real(dp), allocatable :: xx(:)
  character(LEN=80) :: fileloc
  character(LEN=5) :: nchar, ncharcpu

  ! Per-level local gravity buffers: phi + force from assigned files
  type grav_level_t
     real(dp), allocatable :: gdata(:,:)  ! (ngrids, twotondim * ngrav_per_oct)
  end type
  type(grav_level_t) :: my_gvl(1:MAXLEVEL)
  integer :: lvl_offset(1:MAXLEVEL)

  ! Exchange buffers
  real(dp), allocatable :: sendbuf(:), recvbuf(:)
  integer, allocatable :: scount_v(:), rcount_v(:), sdispl_v(:), rdispl_v(:)

  if(myid==1) write(*,*) 'Binary varcpu poisson restore (distributed): ncpu_file=', ncpu_file

  ilun = 99
  call title(nrestart, nchar)
  ngrav_per_oct = 1 + ndim  ! phi + force(1:ndim)
  nval_per_grid = twotondim * ngrav_per_oct

  ! Allocate per-level buffers
  do ilevel = 1, nlevelmax
     nlocal = 0
     do ifile = 1, varcpu_nfiles_local
        nlocal = nlocal + varcpu_nactive(varcpu_my_files(ifile), ilevel)
     end do
     if(nlocal > 0) allocate(my_gvl(ilevel)%gdata(nlocal, nval_per_grid))
  end do

  ! Read assigned grav files
  lvl_offset = 0
  do ifile = 1, varcpu_nfiles_local
     icpu_file = varcpu_my_files(ifile)
     if(IOGROUPSIZEREP>0) then
        call title(((icpu_file-1)/IOGROUPSIZEREP)+1, ncharcpu)
        fileloc='output_'//TRIM(nchar)//'/group_'//TRIM(ncharcpu)//'/grav_'//TRIM(nchar)//'.out'
     else
        fileloc='output_'//TRIM(nchar)//'/grav_'//TRIM(nchar)//'.out'
     end if
     call title(icpu_file, ncharcpu)
     fileloc=TRIM(fileloc)//TRIM(ncharcpu)

     open(unit=ilun, file=fileloc, form='unformatted')
     read(ilun) ncpu2
     read(ilun) ndim2
     read(ilun) nlevelmax2
     read(ilun) nboundary2

     do ilevel = 1, nlevelmax2
        do ibound = 1, nboundary2 + ncpu2
           read(ilun)  ! ilevel2
           read(ilun) ncache

           if(ncache > 0) then
              if(ibound == icpu_file) then
                 ! Active block: read phi + force into buffer
                 allocate(xx(1:ncache))
                 ind = lvl_offset(ilevel)
                 do iskip = 1, twotondim
                    ! phi (1 record)
                    read(ilun) xx
                    do i = 1, ncache
                       my_gvl(ilevel)%gdata(ind+i, (iskip-1)*ngrav_per_oct+1) = xx(i)
                    end do
                    ! force (ndim records)
                    do ivar = 1, ndim
                       read(ilun) xx
                       do i = 1, ncache
                          my_gvl(ilevel)%gdata(ind+i, (iskip-1)*ngrav_per_oct+1+ivar) = xx(i)
                       end do
                    end do
                 end do
                 lvl_offset(ilevel) = ind + ncache
                 deallocate(xx)
              else
                 ! Virtual: skip twotondim * (1+ndim) records
                 do i = 1, twotondim * ngrav_per_oct
                    read(ilun)
                 end do
              end if
           end if
        end do
     end do
     close(ilun)
  end do

  ! For each level: pack, exchange, scatter
  do ilevel = 1, nlevelmax
     if(.not. allocated(varcpu_exc(ilevel)%scount)) cycle

     nsend = varcpu_exc(ilevel)%nsend
     nrecv = varcpu_exc(ilevel)%nrecv

     ! Pack send buffer using send_order from AMR exchange
     allocate(sendbuf(nval_per_grid * max(nsend, 1)))
     do j = 1, nsend
        i = varcpu_exc(ilevel)%send_order(j)
        sendbuf((j-1)*nval_per_grid+1 : j*nval_per_grid) = &
             my_gvl(ilevel)%gdata(i, 1:nval_per_grid)
     end do

     ! MPI_ALLTOALLV exchange
     allocate(recvbuf(nval_per_grid * max(nrecv, 1)))
     allocate(scount_v(0:ncpu-1), rcount_v(0:ncpu-1))
     allocate(sdispl_v(0:ncpu-1), rdispl_v(0:ncpu-1))
     scount_v = varcpu_exc(ilevel)%scount * nval_per_grid
     rcount_v = varcpu_exc(ilevel)%rcount * nval_per_grid
     sdispl_v = varcpu_exc(ilevel)%sdispl * nval_per_grid
     rdispl_v = varcpu_exc(ilevel)%rdispl * nval_per_grid
     call MPI_ALLTOALLV(sendbuf, scount_v, sdispl_v, MPI_DOUBLE_PRECISION, &
          recvbuf, rcount_v, rdispl_v, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, info)
     deallocate(sendbuf, scount_v, rcount_v, sdispl_v, rdispl_v)

     ! Scatter to local grids
     do i = 1, nrecv
        igrid = varcpu_exc(ilevel)%recv_igrid(i)
        if(igrid == 0) cycle

        do iskip = 1, twotondim
           icell = igrid + ncoarse + (iskip-1)*ngridmax
           base = (i-1)*nval_per_grid + (iskip-1)*ngrav_per_oct
           phi(icell) = recvbuf(base + 1)
           do ivar = 1, ndim
              f(icell, ivar) = recvbuf(base + 1 + ivar)
           end do
        end do
     end do
     deallocate(recvbuf)

     ! Free level buffer
     if(allocated(my_gvl(ilevel)%gdata)) deallocate(my_gvl(ilevel)%gdata)
  end do

  ! Free all exchange metadata — no longer needed
  do ilevel = 1, nlevelmax
     if(allocated(varcpu_exc(ilevel)%scount)) deallocate(varcpu_exc(ilevel)%scount)
     if(allocated(varcpu_exc(ilevel)%rcount)) deallocate(varcpu_exc(ilevel)%rcount)
     if(allocated(varcpu_exc(ilevel)%sdispl)) deallocate(varcpu_exc(ilevel)%sdispl)
     if(allocated(varcpu_exc(ilevel)%rdispl)) deallocate(varcpu_exc(ilevel)%rdispl)
     if(allocated(varcpu_exc(ilevel)%send_order)) deallocate(varcpu_exc(ilevel)%send_order)
     if(allocated(varcpu_exc(ilevel)%recv_igrid)) deallocate(varcpu_exc(ilevel)%recv_igrid)
  end do
  if(allocated(varcpu_exc)) deallocate(varcpu_exc)
  if(allocated(varcpu_nactive)) deallocate(varcpu_nactive)
  if(allocated(varcpu_my_files)) deallocate(varcpu_my_files)
  if(allocated(varcpu_ngrid_file)) deallocate(varcpu_ngrid_file)

  if(myid==1) write(*,*) 'Binary varcpu poisson restore done.'

end subroutine restore_poisson_binary_varcpu

