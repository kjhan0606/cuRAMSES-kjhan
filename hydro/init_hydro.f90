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
  uold=0.0d0; unew=0.0d0
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
        return
     end if
#endif
     if(varcpu_restart) then
        call restore_hydro_binary_varcpu()
#ifndef WITHOUTMPI
        call MPI_BARRIER(MPI_COMM_WORLD,info)
#endif
        if(verbose)write(*,*)'Binary varcpu HYDRO backup files read completed'
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
  ! Distributed I/O version: each rank reads its assigned hydro files,
  ! then MPI_ALLTOALLV exchanges data to correct owners using stored
  ! AMR exchange metadata. Primitive → conservative conversion on receive.
  !--------------------------------------------------------------
  use amr_commons
  use hydro_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer :: icpu_file, ilevel, ibound, i, j, ind, ivar, irad, info, ifile
  integer :: ncache, ilun, nvar2, nvar_send, ndim2, nlevelmax2, nboundary2, ncpu2
  integer :: igrid, iskip, icell, nval_per_grid, base
  integer :: nlocal, nsend, nrecv
  real(dp) :: gamma2, eval, rho_val
  real(dp), allocatable :: xx(:)
  character(LEN=80) :: fileloc
  character(LEN=5) :: nchar, ncharcpu

  ! Per-level local hydro buffers: raw primitive data from assigned files
  type hydro_level_t
     real(dp), allocatable :: udata(:,:)  ! (ngrids, twotondim * nvar_send)
  end type
  type(hydro_level_t) :: my_hvl(1:MAXLEVEL)
  integer :: lvl_offset(1:MAXLEVEL)
  integer :: nlvl_grids(1:MAXLEVEL)

  ! Exchange buffers
  real(dp), allocatable :: sendbuf(:), recvbuf(:)
  integer, allocatable :: scount_v(:), rcount_v(:), sdispl_v(:), rdispl_v(:)

  if(myid==1) write(*,*) 'Binary varcpu hydro restore (distributed): ncpu_file=', ncpu_file

  ilun = 99
  call title(nrestart, nchar)

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

  ! Allocate per-level buffers
  nlvl_grids = 0
  do ilevel = 1, nlevelmax
     nlocal = 0
     do ifile = 1, varcpu_nfiles_local
        nlocal = nlocal + varcpu_nactive(varcpu_my_files(ifile), ilevel)
     end do
     nlvl_grids(ilevel) = nlocal
     if(nlocal > 0) allocate(my_hvl(ilevel)%udata(nlocal, nval_per_grid))
  end do

  ! Read assigned hydro files
  lvl_offset = 0
  do ifile = 1, varcpu_nfiles_local
     icpu_file = varcpu_my_files(ifile)
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
                 ! Active block: read raw primitive data into buffer
                 allocate(xx(1:ncache))
                 ind = lvl_offset(ilevel)
                 do iskip = 1, twotondim
                    do ivar = 1, nvar_send
                       read(ilun) xx
                       do i = 1, ncache
                          my_hvl(ilevel)%udata(ind+i, (iskip-1)*nvar_send+ivar) = xx(i)
                       end do
                    end do
                    do ivar = nvar_send+1, nvar2
                       read(ilun)  ! skip extra variables
                    end do
                 end do
                 lvl_offset(ilevel) = ind + ncache
                 deallocate(xx)
              else
                 ! Virtual: skip
                 do i = 1, twotondim * nvar2
                    read(ilun)
                 end do
              end if
           end if
        end do
     end do
     close(ilun)
  end do

  ! For each level: pack, exchange, scatter with primitive→conservative conversion
  do ilevel = 1, nlevelmax
     if(.not. allocated(varcpu_exc(ilevel)%scount)) cycle

     nsend = varcpu_exc(ilevel)%nsend
     nrecv = varcpu_exc(ilevel)%nrecv

     ! Pack send buffer using send_order from AMR exchange
     allocate(sendbuf(nval_per_grid * max(nsend, 1)))
     do j = 1, nsend
        i = varcpu_exc(ilevel)%send_order(j)
        sendbuf((j-1)*nval_per_grid+1 : j*nval_per_grid) = &
             my_hvl(ilevel)%udata(i, 1:nval_per_grid)
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

     ! Scatter to local grids with primitive → conservative conversion
     do i = 1, nrecv
        igrid = varcpu_exc(ilevel)%recv_igrid(i)
        if(igrid == 0) cycle

        do iskip = 1, twotondim
           icell = igrid + ncoarse + (iskip-1)*ngridmax
           base = (i-1)*nval_per_grid + (iskip-1)*nvar_send

           ! Density (file record 1)
           rho_val = recvbuf(base + 1)
           uold(icell, 1) = rho_val

           ! Velocities → momenta (file records 2..ndim+1)
           do ivar = 2, ndim+1
              uold(icell, ivar) = recvbuf(base + ivar) * max(rho_val, smallr)
           end do

#if NENER>0
           ! Non-thermal pressures → energies (file records ndim+2..ndim+1+nener)
           do irad = 1, nener
              uold(icell, ndim+2+irad) = recvbuf(base + ndim+1+irad) / &
                   (gamma_rad(irad) - 1d0)
           end do
#endif

           ! Thermal pressure → total energy (file record ndim+2+nener)
           eval = recvbuf(base + ndim+2+nener) / (gamma - 1d0)
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
           ! Passive scalars (file records ndim+3+nener..nvar_send)
           do ivar = ndim+3+nener, nvar_send
              uold(icell, ivar) = recvbuf(base + ivar) * max(rho_val, smallr)
           end do
#endif
        end do
     end do
     deallocate(recvbuf)

     ! Free level buffer
     if(allocated(my_hvl(ilevel)%udata)) deallocate(my_hvl(ilevel)%udata)
  end do

  if(myid==1) write(*,*) 'Binary varcpu hydro restore done.'

end subroutine restore_hydro_binary_varcpu

