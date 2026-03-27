subroutine init_hydro
  use amr_commons
  use hydro_commons
#ifdef RT      
  use rt_parameters,only: convert_birth_times
#endif
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer::ncell,ncache,iskip,igrid,i,ilevel,ind,ivar,irad
  integer::nvar2,ilevel2,numbl2,ilun,ibound,istart,info
  integer::ncpu2,ndim2,nlevelmax2,nboundary2
  integer ,dimension(:),allocatable::ind_grid
  real(dp),dimension(:),allocatable::xx
  real(dp)::gamma2
  character(LEN=80)::fileloc
  character(LEN=5)::nchar,ncharcpu
  integer,parameter::tag=1108
  integer::dummy_io,info2

  if(verbose)write(*,*)'Entering init_hydro'
  
  !------------------------------------------------------
  ! Allocate conservative, cell-centered variables arrays
  !------------------------------------------------------
  ncell=ncoarse+twotondim*ngridmax
  allocate(uold(1:ncell,1:nvar))
  allocate(unew(1:ncell,1:nvar))
  ! uold/unew: Active cells initialized by restart reader or init_flow_fine.
  ! Free-list cells get zero from mmap(MAP_ANONYMOUS) lazy page allocation.
  ! Skip full-array zeroing to avoid paging in 18 GB at startup.
  if(pressure_fix)then
     allocate(divu(1:ncell))
     allocate(enew(1:ncell))
     divu=0.0d0; enew=0.0d0
  end if

  !--------------------------------
  ! For a restart, read hydro file
  !--------------------------------
  if(nrestart>0)then
#ifdef HDF5
     if(informat == 'hdf5') then
        call restore_hydro_hdf5()
        if(verbose)write(*,*)'HDF5 HYDRO backup files read completed'
        call sgs_init_restart
        return
     end if
#endif
     if(varcpu_restart) then
        call restore_hydro_binary_varcpu()
#ifndef WITHOUTMPI
        call MPI_BARRIER(MPI_COMM_WORLD,info)
#endif
        if(verbose)write(*,*)'Binary varcpu HYDRO backup files read completed'
        call sgs_init_restart
        return
     end if
     ilun=ncpu+myid+10
     call title(nrestart,nchar)

     if(IOGROUPSIZEREP>0)then
        call title(((myid-1)/IOGROUPSIZEREP)+1,ncharcpu)
        fileloc='output_'//TRIM(nchar)//'/group_'//TRIM(ncharcpu)//'/hydro_'//TRIM(nchar)//'.out'
     else
        fileloc='output_'//TRIM(nchar)//'/hydro_'//TRIM(nchar)//'.out'
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
     read(ilun)nvar2
     read(ilun)ndim2
     read(ilun)nlevelmax2
     read(ilun)nboundary2
     read(ilun)gamma2
     if(.not.(neq_chem.or.rt) .and. nvar2.ne.nvar)then
        write(*,*)'File hydro.tmp is not compatible'
        write(*,*)'Found   =',nvar2
        write(*,*)'Expected=',nvar
        call clean_stop
     end if
#ifdef RT
     if((neq_chem.or.rt).and.nvar2.lt.nvar)then ! OK to add ionization fraction vars
        ! Convert birth times for RT postprocessing:
        if(rt.and.static) convert_birth_times=.true.
        if(myid==1) write(*,*)'File hydro.tmp is not compatible'
        if(myid==1) write(*,*)'Found nvar2  =',nvar2
        if(myid==1) write(*,*)'Expected=',nvar
        if(myid==1) write(*,*)'..so only reading first ',nvar2, &
                  'variables and setting the rest to zero'
     end if
     if((neq_chem.or.rt).and.nvar2.gt.nvar)then ! Not OK to drop variables 
        if(myid==1) write(*,*)'File hydro.tmp is not compatible'
        if(myid==1) write(*,*)'Found   =',nvar2
        if(myid==1) write(*,*)'Expected=',nvar
        call clean_stop
     end if
#endif
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
              write(*,*)'File hydro.tmp is not compatible'
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

                 ! Read density and velocities --> density and momenta
                 do ivar=1,ndim+1
                    read(ilun)xx
                    if(ivar==1)then
                       do i=1,ncache
                          uold(ind_grid(i)+iskip,1)=xx(i)
                       end do
                    else if(ivar>=2.and.ivar<=ndim+1)then
                       do i=1,ncache
                          uold(ind_grid(i)+iskip,ivar)=xx(i)*max(uold(ind_grid(i)+iskip,1),smallr)
                       end do
                    endif
                 end do

#if NENER>0
                 ! Read non-thermal pressures --> non-thermal energies
                 do ivar=ndim+3,ndim+2+nener
                    read(ilun)xx
                    do i=1,ncache
                       uold(ind_grid(i)+iskip,ivar)=xx(i)/(gamma_rad(ivar-ndim-2)-1d0)
                    end do
                 end do
#endif
                 ! Read thermal pressure --> total fluid energy
                 read(ilun)xx
                 do i=1,ncache
                    xx(i)=xx(i)/(gamma-1d0)
                    if (uold(ind_grid(i)+iskip,1)>0.)then
                    xx(i)=xx(i)+0.5d0*uold(ind_grid(i)+iskip,2)**2/max(uold(ind_grid(i)+iskip,1),smallr)
#if NDIM>1
                    xx(i)=xx(i)+0.5d0*uold(ind_grid(i)+iskip,3)**2/max(uold(ind_grid(i)+iskip,1),smallr)
#endif
#if NDIM>2
                    xx(i)=xx(i)+0.5d0*uold(ind_grid(i)+iskip,4)**2/max(uold(ind_grid(i)+iskip,1),smallr)
#endif
#if NENER>0
                    do irad=1,nener
                       xx(i)=xx(i)+uold(ind_grid(i)+iskip,ndim+2+irad)
                    end do
#endif
                 else
                    xx(i)=0.
                 end if
                    uold(ind_grid(i)+iskip,ndim+2)=xx(i)
                 end do
#if NVAR>NDIM+2+NENER
                 ! Read passive scalars
                 do ivar=ndim+3+nener,min(nvar,nvar2)
                    read(ilun)xx
                    do i=1,ncache
                       uold(ind_grid(i)+iskip,ivar)=xx(i)*max(uold(ind_grid(i)+iskip,1),smallr)
                    end do
                 end do
#endif
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
     if(debug)write(*,*)'hydro.tmp read for processor ',myid
     call MPI_BARRIER(MPI_COMM_WORLD,info)
#endif
     if(verbose)write(*,*)'HYDRO backup files read completed'
  end if

end subroutine init_hydro

subroutine restore_hydro_binary_varcpu
  !--------------------------------------------------------------
  ! Chunked distributed I/O version: reads hydro files in chunks,
  ! exchanges via ksection_exchange_dp (O(log_k ncpu) memory),
  ! and uses Morton hash lookup for position→igrid mapping.
  !--------------------------------------------------------------
  use amr_commons
  use hydro_commons
  use ksection
  use morton_keys
  use morton_hash
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer :: icpu_file, ilevel, ibound, i, j, k, ind, ivar, irad, info, ifile
  integer :: ncache, ilun, nvar2, nvar_send, ndim2, nlevelmax2, nboundary2, ncpu2
  integer :: igrid, iskip, icell, nval_per_grid, nprops, base
  integer :: nlocal, nrecv
  real(dp) :: gamma2, eval, rho_val, twotol
  real(dp), allocatable :: xx(:)
  character(LEN=80) :: fileloc
  character(LEN=5) :: nchar, ncharcpu

  ! Chunking variables
  integer :: nchunk, ichunk, chunk_file_lo, chunk_file_hi
  integer :: chunk_ngrids(1:MAXLEVEL), chunk_lvl_offset(1:MAXLEVEL)
  integer :: n_from_file, xg_base

  ! Per-level chunk hydro buffers
  type hydro_level_t
     real(dp), allocatable :: udata(:,:)  ! (ngrids, nval_per_grid)
  end type
  type(hydro_level_t) :: chunk_hvl(1:MAXLEVEL)

  ! Ksection exchange buffers
  real(dp), allocatable :: sendbuf_2d(:,:), recvbuf_2d(:,:)
  integer, allocatable :: dest(:)

  ! Morton hash lookup
  integer(8) :: ixm, iym, izm
  type(mkey_t) :: mkey
  real(dp) :: xg_recv(3), xx_father(1:1, 1:ndim), scale
  integer :: c_tmp(1:1), nx_loc, nxny

  if(myid==1) write(*,*) 'Binary varcpu hydro restore (chunked ksection): ncpu_file=', ncpu_file

  ilun = 99
  call title(nrestart, nchar)
  nx_loc = icoarse_max - icoarse_min + 1
  scale = boxlen / dble(nx_loc)
  nxny = nx * ny

  ! Get nvar2 from file 00001
  if(myid == 1) then
     if(IOGROUPSIZEREP>0) then
        call title(1, ncharcpu)
        fileloc='output_'//TRIM(nchar)//'/group_'//TRIM(ncharcpu)//'/hydro_'//TRIM(nchar)//'.out'
     else
        fileloc='output_'//TRIM(nchar)//'/hydro_'//TRIM(nchar)//'.out'
     end if
     call title(1, ncharcpu)
     fileloc=TRIM(fileloc)//TRIM(ncharcpu)
     open(unit=ilun, file=fileloc, form='unformatted')
     read(ilun)  ! ncpu
     read(ilun) nvar2
     close(ilun)
  end if
  call MPI_BCAST(nvar2, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, info)

  nvar_send = min(nvar, nvar2)
  nval_per_grid = twotondim * nvar_send
  nprops = ndim + nval_per_grid  ! xg(1:3) prepended to data

  ! Compute chunk boundaries on global file indices
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

     ! Allocate chunk-level hydro buffers
     do ilevel = 1, nlevelmax
        if(chunk_ngrids(ilevel) > 0) &
             allocate(chunk_hvl(ilevel)%udata(chunk_ngrids(ilevel), nval_per_grid))
     end do

     ! Read assigned hydro files in this chunk
     chunk_lvl_offset = 0
     do ifile = 1, varcpu_nfiles_local
        icpu_file = varcpu_my_files(ifile)
        if(icpu_file < chunk_file_lo .or. icpu_file > chunk_file_hi) cycle

        if(IOGROUPSIZEREP>0) then
           call title(((icpu_file-1)/IOGROUPSIZEREP)+1, ncharcpu)
           fileloc='output_'//TRIM(nchar)//'/group_'//TRIM(ncharcpu)//'/hydro_'//TRIM(nchar)//'.out'
        else
           fileloc='output_'//TRIM(nchar)//'/hydro_'//TRIM(nchar)//'.out'
        end if
        call title(icpu_file, ncharcpu)
        fileloc=TRIM(fileloc)//TRIM(ncharcpu)

        open(unit=ilun, file=fileloc, form='unformatted')
        read(ilun) ncpu2
        read(ilun) nvar2
        read(ilun) ndim2
        read(ilun) nlevelmax2
        read(ilun) nboundary2
        read(ilun) gamma2

        do ilevel = 1, nlevelmax2
           do ibound = 1, nboundary2 + ncpu2
              read(ilun)  ! ilevel2
              read(ilun) ncache

              if(ncache > 0) then
                 if(ibound == icpu_file) then
                    allocate(xx(1:ncache))
                    ind = chunk_lvl_offset(ilevel)
                    do iskip = 1, twotondim
                       do ivar = 1, nvar_send
                          read(ilun) xx
                          do i = 1, ncache
                             chunk_hvl(ilevel)%udata(ind+i, (iskip-1)*nvar_send+ivar) = xx(i)
                          end do
                       end do
                       do ivar = nvar_send+1, nvar2
                          read(ilun)  ! skip extra variables
                       end do
                    end do
                    chunk_lvl_offset(ilevel) = ind + ncache
                    deallocate(xx)
                 else
                    do i = 1, twotondim * nvar2
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

        ! Pack sendbuf: xg(1:ndim) + udata(1:nval_per_grid)
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
              ! Grid position from varcpu_lvl (saved during AMR restore)
              sendbuf_2d(1, k) = varcpu_lvl(ilevel)%xg(xg_base + i, 1)
              sendbuf_2d(2, k) = varcpu_lvl(ilevel)%xg(xg_base + i, 2)
              sendbuf_2d(3, k) = varcpu_lvl(ilevel)%xg(xg_base + i, 3)
              ! Hydro data
              sendbuf_2d(ndim+1:nprops, k) = chunk_hvl(ilevel)%udata(k, 1:nval_per_grid)
              ! Determine owner CPU (same father-level lookup as AMR)
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

        ! Scatter to local grids with primitive → conservative conversion
        do i = 1, nrecv
           xg_recv(1:ndim) = recvbuf_2d(1:ndim, i)

           ! Morton hash lookup → igrid
           ixm = int(xg_recv(1) * twotol, 8)
           iym = int(xg_recv(2) * twotol, 8)
           izm = int(xg_recv(3) * twotol, 8)
           mkey = morton_encode(ixm, iym, izm)
           igrid = morton_hash_lookup(mort_table(ilevel), mkey)
           if(igrid == 0) cycle  ! virtual grid, skip

           do iskip = 1, twotondim
              icell = igrid + ncoarse + (iskip-1)*ngridmax
              base = ndim + (iskip-1)*nvar_send

              ! Density (file record 1)
              rho_val = recvbuf_2d(base + 1, i)
              uold(icell, 1) = rho_val

              ! Velocities → momenta (file records 2..ndim+1)
              do ivar = 2, ndim+1
                 uold(icell, ivar) = recvbuf_2d(base + ivar, i) * max(rho_val, smallr)
              end do

#if NENER>0
              ! Non-thermal pressures → energies
              do irad = 1, nener
                 uold(icell, ndim+2+irad) = recvbuf_2d(base + ndim+1+irad, i) / &
                      (gamma_rad(irad) - 1d0)
              end do
#endif

              ! Thermal pressure → total energy
              eval = recvbuf_2d(base + ndim+2+nener, i) / (gamma - 1d0)
              if(rho_val > 0d0) then
                 do ivar = 2, ndim+1
                    eval = eval + 0.5d0*uold(icell,ivar)**2 / max(rho_val, smallr)
                 end do
#if NENER>0
                 do irad = 1, nener
                    eval = eval + uold(icell, ndim+2+irad)
                 end do
#endif
              else
                 eval = 0d0
              end if
              uold(icell, ndim+2) = eval

#if NVAR>NDIM+2+NENER
              ! Passive scalars
              do ivar = ndim+3+nener, nvar_send
                 uold(icell, ivar) = recvbuf_2d(base + ivar, i) * max(rho_val, smallr)
              end do
#endif
           end do
        end do
        deallocate(recvbuf_2d)

        ! Free chunk level buffer
        if(allocated(chunk_hvl(ilevel)%udata)) deallocate(chunk_hvl(ilevel)%udata)
     end do
  end do  ! ichunk

  if(myid==1) write(*,*) 'Binary varcpu hydro restore done.'

end subroutine restore_hydro_binary_varcpu

