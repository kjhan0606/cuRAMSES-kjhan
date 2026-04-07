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
  ! rho/phi/f: recomputed every timestep by rho_fine/Poisson solver.
  ! Free-list cells get zero from mmap(MAP_ANONYMOUS) lazy page allocation.
  ! Skip full-array zeroing to avoid paging in ~9 GB at startup.
  if(use_fR .or. use_nDGP .or. use_symmetron .or. use_dilaton .or. use_galileon) then
     allocate(scalar_gr    (1:ncell))
     allocate(scalar_gr_old(1:ncell))
     ! scalar_gr: recomputed by MG solver. Mmap provides zero pages.
  end if
  if(use_fdm) then
     allocate(psi_re(1:ncell))
     allocate(psi_im(1:ncell))
     ! psi: mmap provides zero pages; initialized by fdm_init_psi or restore
  end if
  if(cic_levelmax>0)then
     allocate(rho_top(1:ncell))
     ! rho_top: set by rho_fine. Mmap provides zero pages.
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
  ! Chunked distributed I/O version: reads grav files in chunks,
  ! exchanges via ksection_exchange_dp (O(log_k ncpu) memory),
  ! uses Morton hash lookup for position→igrid mapping.
  ! Frees all varcpu metadata after completion.
  !--------------------------------------------------------------
  use amr_commons
  use poisson_commons
  use ksection
  use morton_keys
  use morton_hash
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer :: icpu_file, ilevel, ibound, i, j, k, ind, ivar, info, ifile
  integer :: ncache, ilun, ndim2, nlevelmax2, nboundary2, ncpu2
  integer :: igrid, iskip, icell, nval_per_grid, ngrav_per_oct, nprops, base
  integer :: nlocal, nrecv
  real(dp), allocatable :: xx(:)
  character(LEN=80) :: fileloc
  character(LEN=5) :: nchar, ncharcpu

  ! Chunking variables
  integer :: nchunk, ichunk, chunk_file_lo, chunk_file_hi
  integer :: chunk_ngrids(1:MAXLEVEL), chunk_lvl_offset(1:MAXLEVEL)
  integer :: n_from_file, xg_base

  ! Per-level chunk gravity buffers
  type grav_level_t
     real(dp), allocatable :: gdata(:,:)  ! (ngrids, nval_per_grid)
  end type
  type(grav_level_t) :: chunk_gvl(1:MAXLEVEL)

  ! Ksection exchange buffers
  real(dp), allocatable :: sendbuf_2d(:,:), recvbuf_2d(:,:)
  integer, allocatable :: dest(:)

  ! Morton hash lookup
  integer(8) :: ixm, iym, izm
  type(mkey_t) :: mkey
  real(dp) :: xg_recv(3), xx_father(1:1, 1:ndim), scale, twotol
  integer :: c_tmp(1:1), nx_loc, nxny

  if(myid==1) write(*,*) 'Binary varcpu poisson restore (chunked ksection): ncpu_file=', ncpu_file

  ilun = 99
  call title(nrestart, nchar)
  nx_loc = icoarse_max - icoarse_min + 1
  scale = boxlen / dble(nx_loc)
  nxny = nx * ny
  ngrav_per_oct = 1 + ndim  ! phi + force(1:ndim)
  nval_per_grid = twotondim * ngrav_per_oct
  nprops = ndim + nval_per_grid  ! xg(1:3) prepended to data

  ! Compute chunk boundaries
  if(varcpu_chunk_nfile <= 0) then
     nchunk = 1
  else
     nchunk = (ncpu_file + varcpu_chunk_nfile - 1) / varcpu_chunk_nfile
  end if

  ! Main chunked loop
  do ichunk = 1, nchunk
     if(varcpu_chunk_nfile <= 0) then
        chunk_file_lo = 1
        chunk_file_hi = ncpu_file
     else
        chunk_file_lo = (ichunk - 1) * varcpu_chunk_nfile + 1
        chunk_file_hi = min(ichunk * varcpu_chunk_nfile, ncpu_file)
     end if

     ! Compute per-level grid counts for this chunk (my files only)
     chunk_ngrids = 0
     do j = 1, varcpu_nfiles_local
        if(varcpu_my_files(j) < chunk_file_lo .or. &
           varcpu_my_files(j) > chunk_file_hi) cycle
        do ilevel = 1, nlevelmax
           chunk_ngrids(ilevel) = chunk_ngrids(ilevel) + &
                varcpu_nactive(varcpu_my_files(j), ilevel)
        end do
     end do

     ! Allocate chunk-level gravity buffers
     do ilevel = 1, nlevelmax
        if(chunk_ngrids(ilevel) > 0) &
             allocate(chunk_gvl(ilevel)%gdata(chunk_ngrids(ilevel), nval_per_grid))
     end do

     ! Read assigned grav files in this chunk
     chunk_lvl_offset = 0
     do ifile = 1, varcpu_nfiles_local
        icpu_file = varcpu_my_files(ifile)
        if(icpu_file < chunk_file_lo .or. icpu_file > chunk_file_hi) cycle

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
                    allocate(xx(1:ncache))
                    ind = chunk_lvl_offset(ilevel)
                    do iskip = 1, twotondim
                       ! phi (1 record)
                       read(ilun) xx
                       do i = 1, ncache
                          chunk_gvl(ilevel)%gdata(ind+i, (iskip-1)*ngrav_per_oct+1) = xx(i)
                       end do
                       ! force (ndim records)
                       do ivar = 1, ndim
                          read(ilun) xx
                          do i = 1, ncache
                             chunk_gvl(ilevel)%gdata(ind+i, (iskip-1)*ngrav_per_oct+1+ivar) = xx(i)
                          end do
                       end do
                    end do
                    chunk_lvl_offset(ilevel) = ind + ncache
                    deallocate(xx)
                 else
                    do i = 1, twotondim * ngrav_per_oct
                       read(ilun)
                    end do
                 end if
              end if
           end do
        end do
        close(ilun)
     end do

     ! Exchange and scatter level by level
     do ilevel = 1, nlevelmax
        if(varcpu_ngrid_file(ilevel) == 0) cycle

        nlocal = chunk_ngrids(ilevel)
        twotol = 2.0d0**(ilevel-1)

        ! Pack sendbuf: xg(1:ndim) + gdata(1:nval_per_grid)
        allocate(sendbuf_2d(1:nprops, 1:max(nlocal,1)))
        allocate(dest(max(nlocal,1)))
        k = 0
        do j = 1, varcpu_nfiles_local
           if(varcpu_my_files(j) < chunk_file_lo .or. &
              varcpu_my_files(j) > chunk_file_hi) cycle
           n_from_file = varcpu_nactive(varcpu_my_files(j), ilevel)
           xg_base = varcpu_file_start(j-1, ilevel)
           do i = 1, n_from_file
              k = k + 1
              sendbuf_2d(1, k) = varcpu_lvl(ilevel)%xg(xg_base + i, 1)
              sendbuf_2d(2, k) = varcpu_lvl(ilevel)%xg(xg_base + i, 2)
              sendbuf_2d(3, k) = varcpu_lvl(ilevel)%xg(xg_base + i, 3)
              sendbuf_2d(ndim+1:nprops, k) = chunk_gvl(ilevel)%gdata(k, 1:nval_per_grid)
              ! Determine owner CPU
              ixm = int(sendbuf_2d(1, k) * twotol, 8)
              iym = int(sendbuf_2d(2, k) * twotol, 8)
              izm = int(sendbuf_2d(3, k) * twotol, 8)
              xx_father(1,1) = ((dble(ixm) + 0.5d0) / twotol - dble(icoarse_min)) * scale
              xx_father(1,2) = ((dble(iym) + 0.5d0) / twotol - dble(jcoarse_min)) * scale
              xx_father(1,3) = ((dble(izm) + 0.5d0) / twotol - dble(kcoarse_min)) * scale
              if(ordering == 'ksection') then
                 call cmp_ksection_cpumap(xx_father, c_tmp, 1)
              else
                 call cmp_cpumap(xx_father, c_tmp, 1)
              end if
              dest(k) = c_tmp(1)
           end do
        end do

        ! Ksection hierarchical exchange
        call ksection_exchange_dp(sendbuf_2d, nlocal, dest, nprops, recvbuf_2d, nrecv)
        deallocate(sendbuf_2d, dest)

        ! Scatter to local grids
        do i = 1, nrecv
           xg_recv(1:ndim) = recvbuf_2d(1:ndim, i)

           ! Morton hash lookup → igrid
           ixm = int(xg_recv(1) * twotol, 8)
           iym = int(xg_recv(2) * twotol, 8)
           izm = int(xg_recv(3) * twotol, 8)
           mkey = morton_encode(ixm, iym, izm)
           igrid = morton_hash_lookup(mort_table(ilevel), mkey)
           if(igrid == 0) cycle

           do iskip = 1, twotondim
              icell = igrid + ncoarse + (iskip-1)*ngridmax
              base = ndim + (iskip-1)*ngrav_per_oct
              phi(icell) = recvbuf_2d(base + 1, i)
              do ivar = 1, ndim
                 f(icell, ivar) = recvbuf_2d(base + 1 + ivar, i)
              end do
           end do
        end do
        deallocate(recvbuf_2d)

        ! Free chunk level buffer
        if(allocated(chunk_gvl(ilevel)%gdata)) deallocate(chunk_gvl(ilevel)%gdata)
     end do
  end do  ! ichunk

  ! Free all varcpu metadata — no longer needed
  do ilevel = 1, nlevelmax
     if(allocated(varcpu_lvl(ilevel)%xg)) deallocate(varcpu_lvl(ilevel)%xg)
  end do
  if(allocated(varcpu_lvl)) deallocate(varcpu_lvl)
  if(allocated(varcpu_file_start)) deallocate(varcpu_file_start)
  if(allocated(varcpu_nactive)) deallocate(varcpu_nactive)
  if(allocated(varcpu_my_files)) deallocate(varcpu_my_files)
  if(allocated(varcpu_ngrid_file)) deallocate(varcpu_ngrid_file)

  if(myid==1) write(*,*) 'Binary varcpu poisson restore done.'

end subroutine restore_poisson_binary_varcpu

