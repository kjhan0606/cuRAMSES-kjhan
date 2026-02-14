!###########################################################################
! backup_hdf5.f90 — HDF5 parallel output for cuRAMSES
!
! dump_all_hdf5(filedir) creates a single HDF5 file per output:
!   output_NNNNN/data_NNNNN.h5
!###########################################################################
#ifdef HDF5
subroutine dump_all_hdf5(filedir, nchar)
  use amr_commons
  use pm_commons
  use hydro_commons
  use ramses_hdf5_io
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  character(len=*), intent(in) :: filedir
  character(len=5), intent(in) :: nchar
  character(len=200) :: h5filename
  integer :: info

  h5filename = trim(filedir)//'data_'//trim(nchar)//'.h5'

  if(myid==1) write(*,*) 'HDF5 output: ', trim(h5filename)

  ! Create file with parallel access
  call hdf5_create_parallel(h5filename, MPI_COMM_WORLD)

  ! Write sections
  call backup_header_hdf5()
  call backup_amr_hdf5()
  if(hydro) call backup_hydro_hdf5()
  if(poisson) call backup_poisson_hdf5()
  if(pic) call backup_part_hdf5()
  if(sink) call backup_sink_hdf5()

  call hdf5_close_file()

#ifndef WITHOUTMPI
  call MPI_BARRIER(MPI_COMM_WORLD, info)
#endif
  if(myid==1) write(*,*) 'HDF5 output done.'

end subroutine dump_all_hdf5

!###########################################################################
! Header: simulation parameters as attributes
!###########################################################################
subroutine backup_header_hdf5()
  use amr_commons
  use hydro_commons
  use hydro_parameters, only: gamma
  use pm_commons
  use ramses_hdf5_io
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer(HID_T) :: grp_id

  call hdf5_create_group('/header', grp_id)

  ! Grid parameters
  call hdf5_write_attr_int(grp_id, 'ncpu', ncpu)
  call hdf5_write_attr_int(grp_id, 'ndim', ndim)
  call hdf5_write_attr_int(grp_id, 'nx', nx)
  call hdf5_write_attr_int(grp_id, 'ny', ny)
  call hdf5_write_attr_int(grp_id, 'nz', nz)
  call hdf5_write_attr_int(grp_id, 'nlevelmax', nlevelmax)
  call hdf5_write_attr_int(grp_id, 'ngridmax', ngridmax)
  call hdf5_write_attr_int(grp_id, 'nboundary', nboundary)
  call hdf5_write_attr_int(grp_id, 'ncoarse', ncoarse)
  call hdf5_write_attr_int(grp_id, 'ngrid_current', ngrid_current)
  call hdf5_write_attr_dp(grp_id, 'boxlen', boxlen)

  ! Hydro
  call hdf5_write_attr_int(grp_id, 'nvar', nvar)
  call hdf5_write_attr_dp(grp_id, 'gamma', gamma)

  ! Time
  call hdf5_write_attr_dp(grp_id, 'time', t)
  call hdf5_write_attr_int(grp_id, 'nstep', nstep)
  call hdf5_write_attr_int(grp_id, 'nstep_coarse', nstep_coarse)
  call hdf5_write_attr_int(grp_id, 'noutput', noutput)
  call hdf5_write_attr_int(grp_id, 'iout', iout)
  call hdf5_write_attr_int(grp_id, 'ifout', ifout)
  call hdf5_write_attr_1d_dp(grp_id, 'tout', tout, noutput)
  call hdf5_write_attr_1d_dp(grp_id, 'aout', aout, noutput)
  call hdf5_write_attr_1d_dp(grp_id, 'dtold', dtold, nlevelmax)
  call hdf5_write_attr_1d_dp(grp_id, 'dtnew', dtnew, nlevelmax)

  ! Cosmology
  call hdf5_write_attr_dp(grp_id, 'const', const)
  call hdf5_write_attr_dp(grp_id, 'mass_tot_0', mass_tot_0)
  call hdf5_write_attr_dp(grp_id, 'rho_tot', rho_tot)
  call hdf5_write_attr_dp(grp_id, 'omega_m', omega_m)
  call hdf5_write_attr_dp(grp_id, 'omega_l', omega_l)
  call hdf5_write_attr_dp(grp_id, 'omega_k', omega_k)
  call hdf5_write_attr_dp(grp_id, 'omega_b', omega_b)
  call hdf5_write_attr_dp(grp_id, 'h0', h0)
  call hdf5_write_attr_dp(grp_id, 'aexp_ini', aexp_ini)
  call hdf5_write_attr_dp(grp_id, 'boxlen_ini', boxlen_ini)
  call hdf5_write_attr_dp(grp_id, 'aexp', aexp)
  call hdf5_write_attr_dp(grp_id, 'hexp', hexp)
  call hdf5_write_attr_dp(grp_id, 'aexp_old', aexp_old)
  call hdf5_write_attr_dp(grp_id, 'epot_tot_int', epot_tot_int)
  call hdf5_write_attr_dp(grp_id, 'epot_tot_old', epot_tot_old)
  call hdf5_write_attr_dp(grp_id, 'mass_sph', mass_sph)

  ! Ordering
  call hdf5_write_attr_string(grp_id, 'ordering', ordering)

  ! Flags
  call hdf5_write_attr_int(grp_id, 'levelmin', levelmin)

  call hdf5_close_group(grp_id)

end subroutine backup_header_hdf5

!###########################################################################
! AMR tree: coarse grid + per-level grid data
! Each level stores: xg, son_flag, cpu_map, ngrid_per_cpu
!###########################################################################
subroutine backup_amr_hdf5()
  use amr_commons
  use ksection
  use ramses_hdf5_io
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer :: ilevel, i, igrid, ind, iskip, idim, info
  integer :: ngrid_loc, ncache, nlevelmax_file
  integer, allocatable :: ngrid_all(:)
  integer(i8b) :: ngrid_total, offset_grid
  integer(i8b) :: ncells_total, offset_cells
  integer, allocatable :: son_flag_buf(:), cpu_map_buf(:)
  real(dp), allocatable :: xg_buf(:)
  integer(HID_T) :: grp_id, lvl_grp_id
  character(len=40) :: grp_name
  character(len=10) :: lvl_str

  allocate(ngrid_all(ncpu))

  ! Write coarse-level data (all ranks write their portion)
  call hdf5_create_group('/coarse', grp_id)
  call hdf5_write_dataset_1d_int(grp_id, 'son', son(1:ncoarse), ncoarse, &
       int(0, i8b), int(ncoarse, i8b))
  call hdf5_write_dataset_1d_int(grp_id, 'cpu_map', cpu_map(1:ncoarse), ncoarse, &
       int(0, i8b), int(ncoarse, i8b))
  call hdf5_write_dataset_1d_int(grp_id, 'flag1', flag1(1:ncoarse), ncoarse, &
       int(0, i8b), int(ncoarse, i8b))
  call hdf5_close_group(grp_id)

  ! Write domain decomposition info
  call hdf5_create_group('/domain', grp_id)
  call hdf5_write_attr_string(grp_id, 'ordering', ordering)
  if(ordering=='bisection') then
     ! Store bisection tree info (serial from rank 0)
     call hdf5_write_dataset_serial_dp(grp_id, 'bisec_wall', bisec_wall, nbinodes, myid)
  else if(ordering=='ksection') then
     call hdf5_write_attr_int(grp_id, 'nksec_levels', nksec_levels)
     call hdf5_write_attr_int(grp_id, 'ksec_kmax', ksec_kmax)
     call hdf5_write_attr_int(grp_id, 'ksec_nbinodes', ksec_nbinodes)
     call hdf5_write_dataset_serial_int(grp_id, 'ksec_factor', ksec_factor, nksec_levels, myid)
     call hdf5_write_dataset_serial_int(grp_id, 'ksec_dir', ksec_dir, nksec_levels, myid)
     call hdf5_write_dataset_serial_int(grp_id, 'ksec_indx', ksec_indx, ksec_nbinodes, myid)
     ! Flatten 2D arrays column-major: ksec_wall(nbinodes, ncols), ksec_next(nbinodes, ncols)
     block
        integer :: nw_cols, nn_cols, nw_total, nn_total
        real(dp), allocatable :: wall_flat(:)
        integer, allocatable :: next_flat(:)
        nw_cols = max(ksec_kmax - 1, 1)
        nn_cols = ksec_kmax
        nw_total = ksec_nbinodes * nw_cols
        nn_total = ksec_nbinodes * nn_cols
        allocate(wall_flat(nw_total))
        allocate(next_flat(nn_total))
        wall_flat = reshape(ksec_wall(1:ksec_nbinodes, 1:nw_cols), (/nw_total/))
        next_flat = reshape(ksec_next(1:ksec_nbinodes, 1:nn_cols), (/nn_total/))
        call hdf5_write_dataset_serial_dp(grp_id, 'ksec_wall', wall_flat, nw_total, myid)
        call hdf5_write_dataset_serial_int(grp_id, 'ksec_next', next_flat, nn_total, myid)
        deallocate(wall_flat, next_flat)
     end block
     ! CPU bounding boxes
     block
        real(dp), allocatable :: cpubox_flat(:)
        integer :: ncb
        ncb = ncpu * ndim
        allocate(cpubox_flat(ncb))
        cpubox_flat = reshape(bisec_cpubox_min(1:ncpu, 1:ndim), (/ncb/))
        call hdf5_write_dataset_serial_dp(grp_id, 'bisec_cpubox_min', cpubox_flat, ncb, myid)
        cpubox_flat = reshape(bisec_cpubox_max(1:ncpu, 1:ndim), (/ncb/))
        call hdf5_write_dataset_serial_dp(grp_id, 'bisec_cpubox_max', cpubox_flat, ncb, myid)
        deallocate(cpubox_flat)
     end block
  else
     ! Hilbert ordering: write bound_key (all ranks have same copy)
     ! bound_key is real(qdp), write as dp (sufficient for restart)
     block
        real(dp), allocatable :: bkey_dp(:)
        integer :: nd
        nd = ndomain + 1
        allocate(bkey_dp(nd))
        do i = 0, ndomain
           bkey_dp(i+1) = real(bound_key(i), dp)
        end do
        call hdf5_write_dataset_serial_dp(grp_id, 'bound_key', bkey_dp, nd, myid)
        deallocate(bkey_dp)
     end block
  end if
  call hdf5_close_group(grp_id)

  ! Compute nlevelmax_file (max level with grids > 0)
  nlevelmax_file = 0
  do ilevel = 1, nlevelmax
     block
        integer :: ngtmp
        call MPI_ALLREDUCE(numbl(myid, ilevel), ngtmp, 1, MPI_INTEGER, &
             MPI_SUM, MPI_COMM_WORLD, info)
        if(ngtmp > 0) nlevelmax_file = ilevel
     end block
  end do

  ! Per-level AMR data
  call hdf5_create_group('/amr', grp_id)
  call hdf5_write_attr_int(grp_id, 'nlevelmax_file', nlevelmax_file)
  call hdf5_close_group(grp_id)

  do ilevel = 1, nlevelmax
     ! Count local grids at this level (only active, not boundary)
     ngrid_loc = numbl(myid, ilevel)

     ! Allgather to compute global total and offset
     call MPI_ALLGATHER(ngrid_loc, 1, MPI_INTEGER, ngrid_all, 1, MPI_INTEGER, &
          MPI_COMM_WORLD, info)
     ngrid_total = 0
     offset_grid = 0
     do i = 1, ncpu
        if(i < myid) offset_grid = offset_grid + ngrid_all(i)
        ngrid_total = ngrid_total + ngrid_all(i)
     end do

     if(ngrid_total == 0) cycle

     ! Create level group
     write(lvl_str, '(I0)') ilevel
     grp_name = '/amr/level_'//trim(lvl_str)
     call hdf5_create_group(grp_name, lvl_grp_id)

     ! Store ngrid_per_cpu as attribute
     call hdf5_write_attr_int(lvl_grp_id, 'ngrid_total', int(ngrid_total))

     ! Allocate and fill buffers
     ! xg: grid positions [ngrid_loc * ndim], stored as ndim separate datasets
     if(ngrid_loc > 0) then
        allocate(xg_buf(ngrid_loc))
        allocate(son_flag_buf(ngrid_loc * twotondim))
        allocate(cpu_map_buf(ngrid_loc * twotondim))
     else
        allocate(xg_buf(1))
        allocate(son_flag_buf(1))
        allocate(cpu_map_buf(1))
     end if

     ! Write grid positions (ndim datasets)
     do idim = 1, ndim
        igrid = headl(myid, ilevel)
        do i = 1, ngrid_loc
           xg_buf(i) = xg(igrid, idim)
           igrid = next(igrid)
        end do
        write(lvl_str, '(I0)') idim
        call hdf5_write_dataset_1d_dp(lvl_grp_id, 'xg_'//trim(lvl_str), &
             xg_buf, ngrid_loc, offset_grid, ngrid_total)
     end do

     ! Write son_flag (0=leaf, >0=refined) and cpu_map per cell
     ! Layout: twotondim contiguous blocks, each of ngrid_total
     ncells_total = ngrid_total * twotondim
     offset_cells = offset_grid * twotondim

     igrid = headl(myid, ilevel)
     do i = 1, ngrid_loc
        do ind = 1, twotondim
           iskip = ncoarse + (ind - 1) * ngridmax
           ! son_flag: 0 if no son, 1 if has son
           if(son(igrid + iskip) > 0) then
              son_flag_buf((i-1)*twotondim + ind) = 1
           else
              son_flag_buf((i-1)*twotondim + ind) = 0
           end if
           cpu_map_buf((i-1)*twotondim + ind) = cpu_map(igrid + iskip)
        end do
        igrid = next(igrid)
     end do

     call hdf5_write_dataset_1d_int(lvl_grp_id, 'son_flag', &
          son_flag_buf, ngrid_loc * twotondim, offset_cells, ncells_total)
     call hdf5_write_dataset_1d_int(lvl_grp_id, 'cpu_map', &
          cpu_map_buf, ngrid_loc * twotondim, offset_cells, ncells_total)

     ! Write ngrid_per_cpu as dataset (so we know how to split on read)
     call hdf5_write_dataset_serial_int(lvl_grp_id, 'ngrid_per_cpu', ngrid_all, ncpu, myid)

     deallocate(xg_buf, son_flag_buf, cpu_map_buf)
     call hdf5_close_group(lvl_grp_id)
  end do

  deallocate(ngrid_all)

end subroutine backup_amr_hdf5

!###########################################################################
! Hydro: per-level uold (raw conservative variables)
!###########################################################################
subroutine backup_hydro_hdf5()
  use amr_commons
  use hydro_commons
  use ramses_hdf5_io
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer :: ilevel, i, igrid, ind, iskip, ivar, info
  integer :: ngrid_loc
  integer, allocatable :: ngrid_all(:)
  integer(i8b) :: ncells_total, offset_cells, ngrid_total
  real(dp), allocatable :: ubuf(:)
  integer(HID_T) :: grp_id, lvl_grp_id
  character(len=40) :: grp_name
  character(len=10) :: lvl_str, var_str

  allocate(ngrid_all(ncpu))

  call hdf5_create_group('/hydro', grp_id)
  call hdf5_close_group(grp_id)

  do ilevel = 1, nlevelmax
     ngrid_loc = numbl(myid, ilevel)
     call MPI_ALLGATHER(ngrid_loc, 1, MPI_INTEGER, ngrid_all, 1, MPI_INTEGER, &
          MPI_COMM_WORLD, info)
     ngrid_total = 0
     do i = 1, ncpu
        ngrid_total = ngrid_total + ngrid_all(i)
     end do
     if(ngrid_total == 0) cycle

     ncells_total = ngrid_total * twotondim
     offset_cells = 0
     do i = 1, myid - 1
        offset_cells = offset_cells + int(ngrid_all(i), i8b) * twotondim
     end do

     write(lvl_str, '(I0)') ilevel
     grp_name = '/hydro/level_'//trim(lvl_str)
     call hdf5_create_group(grp_name, lvl_grp_id)

     ! Write uold: one dataset per variable, each of size ncells_total
     ! Data layout: grid-major ordering (all cells of grid i, then grid i+1)
     if(ngrid_loc > 0) then
        allocate(ubuf(ngrid_loc * twotondim))
     else
        allocate(ubuf(1))
     end if

     do ivar = 1, nvar
        igrid = headl(myid, ilevel)
        do i = 1, ngrid_loc
           do ind = 1, twotondim
              iskip = ncoarse + (ind - 1) * ngridmax
              ubuf((i-1)*twotondim + ind) = uold(igrid + iskip, ivar)
           end do
           igrid = next(igrid)
        end do
        write(var_str, '(I0)') ivar
        call hdf5_write_dataset_1d_dp(lvl_grp_id, 'uold_'//trim(var_str), &
             ubuf, ngrid_loc * twotondim, offset_cells, ncells_total)
     end do

     deallocate(ubuf)
     call hdf5_close_group(lvl_grp_id)
  end do

  deallocate(ngrid_all)

end subroutine backup_hydro_hdf5

!###########################################################################
! Poisson: per-level phi and f
!###########################################################################
subroutine backup_poisson_hdf5()
  use amr_commons
  use poisson_commons
  use ramses_hdf5_io
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer :: ilevel, i, igrid, ind, iskip, idim, info
  integer :: ngrid_loc
  integer, allocatable :: ngrid_all(:)
  integer(i8b) :: ncells_total, offset_cells, ngrid_total
  real(dp), allocatable :: pbuf(:)
  integer(HID_T) :: grp_id, lvl_grp_id
  character(len=40) :: grp_name
  character(len=10) :: lvl_str, dim_str

  allocate(ngrid_all(ncpu))

  call hdf5_create_group('/gravity', grp_id)
  call hdf5_close_group(grp_id)

  do ilevel = 1, nlevelmax
     ngrid_loc = numbl(myid, ilevel)
     call MPI_ALLGATHER(ngrid_loc, 1, MPI_INTEGER, ngrid_all, 1, MPI_INTEGER, &
          MPI_COMM_WORLD, info)
     ngrid_total = 0
     do i = 1, ncpu
        ngrid_total = ngrid_total + ngrid_all(i)
     end do
     if(ngrid_total == 0) cycle

     ncells_total = ngrid_total * twotondim
     offset_cells = 0
     do i = 1, myid - 1
        offset_cells = offset_cells + int(ngrid_all(i), i8b) * twotondim
     end do

     write(lvl_str, '(I0)') ilevel
     grp_name = '/gravity/level_'//trim(lvl_str)
     call hdf5_create_group(grp_name, lvl_grp_id)

     if(ngrid_loc > 0) then
        allocate(pbuf(ngrid_loc * twotondim))
     else
        allocate(pbuf(1))
     end if

     ! Write phi
     igrid = headl(myid, ilevel)
     do i = 1, ngrid_loc
        do ind = 1, twotondim
           iskip = ncoarse + (ind - 1) * ngridmax
           pbuf((i-1)*twotondim + ind) = phi(igrid + iskip)
        end do
        igrid = next(igrid)
     end do
     call hdf5_write_dataset_1d_dp(lvl_grp_id, 'phi', &
          pbuf, ngrid_loc * twotondim, offset_cells, ncells_total)

     ! Write force components
     do idim = 1, ndim
        igrid = headl(myid, ilevel)
        do i = 1, ngrid_loc
           do ind = 1, twotondim
              iskip = ncoarse + (ind - 1) * ngridmax
              pbuf((i-1)*twotondim + ind) = f(igrid + iskip, idim)
           end do
           igrid = next(igrid)
        end do
        write(dim_str, '(I0)') idim
        call hdf5_write_dataset_1d_dp(lvl_grp_id, 'f_'//trim(dim_str), &
             pbuf, ngrid_loc * twotondim, offset_cells, ncells_total)
     end do

     deallocate(pbuf)
     call hdf5_close_group(lvl_grp_id)
  end do

  deallocate(ngrid_all)

end subroutine backup_poisson_hdf5

!###########################################################################
! Particles: parallel hyperslab write
!###########################################################################
subroutine backup_part_hdf5()
  use amr_commons
  use pm_commons
  use ramses_hdf5_io
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer :: i, idim, ipart, info
  integer :: npart_loc
  integer(i8b) :: npart_total, offset_part, tmp_long
  integer, allocatable :: npart_all(:)
  real(dp), allocatable :: dbuf(:)
  integer(i8b), allocatable :: ibuf8(:)
  integer, allocatable :: ibuf(:)
  integer(HID_T) :: grp_id

  allocate(npart_all(ncpu))

  ! Gather particle counts
  npart_loc = npart
#ifndef WITHOUTMPI
  call MPI_ALLGATHER(npart_loc, 1, MPI_INTEGER, npart_all, 1, MPI_INTEGER, &
       MPI_COMM_WORLD, info)
#else
  npart_all(1) = npart_loc
#endif

  npart_total = 0
  offset_part = 0
  do i = 1, ncpu
     if(i < myid) offset_part = offset_part + npart_all(i)
     npart_total = npart_total + npart_all(i)
  end do

  call hdf5_create_group('/particles', grp_id)
  call hdf5_write_attr_int8(grp_id, 'npart_total', npart_total)
  call hdf5_write_attr_int8(grp_id, 'nstar_tot', int(nstar_tot, i8b))
  call hdf5_write_attr_dp(grp_id, 'mstar_tot', mstar_tot)
  call hdf5_write_attr_dp(grp_id, 'mstar_lost', mstar_lost)
  ! Store npart_per_cpu for fast same-rank restart
  call hdf5_write_dataset_serial_int(grp_id, 'npart_per_cpu', npart_all, ncpu, myid)

  ! Pack active particles into contiguous buffers
  allocate(dbuf(npart_loc))
  allocate(ibuf8(npart_loc))
  allocate(ibuf(npart_loc))

  ! Position
  do idim = 1, ndim
     ipart = 0
     do i = 1, npartmax
        if(levelp(i) > 0) then
           ipart = ipart + 1
           dbuf(ipart) = xp(i, idim)
        end if
     end do
     block
        character(len=10) :: dstr
        write(dstr, '(I0)') idim
        call hdf5_write_dataset_1d_dp(grp_id, 'x_'//trim(dstr), &
             dbuf, npart_loc, offset_part, npart_total)
     end block
  end do

  ! Velocity
  do idim = 1, ndim
     ipart = 0
     do i = 1, npartmax
        if(levelp(i) > 0) then
           ipart = ipart + 1
           dbuf(ipart) = vp(i, idim)
        end if
     end do
     block
        character(len=10) :: dstr
        write(dstr, '(I0)') idim
        call hdf5_write_dataset_1d_dp(grp_id, 'v_'//trim(dstr), &
             dbuf, npart_loc, offset_part, npart_total)
     end block
  end do

  ! Mass
  ipart = 0
  do i = 1, npartmax
     if(levelp(i) > 0) then
        ipart = ipart + 1
        dbuf(ipart) = mp(i)
     end if
  end do
  call hdf5_write_dataset_1d_dp(grp_id, 'mass', &
       dbuf, npart_loc, offset_part, npart_total)

  ! Identity
  ipart = 0
  do i = 1, npartmax
     if(levelp(i) > 0) then
        ipart = ipart + 1
        ibuf8(ipart) = idp(i)
     end if
  end do
  call hdf5_write_dataset_1d_int8(grp_id, 'identity', &
       ibuf8, npart_loc, offset_part, npart_total)

  ! Level
  ipart = 0
  do i = 1, npartmax
     if(levelp(i) > 0) then
        ipart = ipart + 1
        ibuf(ipart) = levelp(i)
     end if
  end do
  call hdf5_write_dataset_1d_int(grp_id, 'levelp', &
       ibuf, npart_loc, offset_part, npart_total)

  ! Birth epoch (tp) and metallicity (zp)
  if(star .or. sink) then
     ipart = 0
     do i = 1, npartmax
        if(levelp(i) > 0) then
           ipart = ipart + 1
           dbuf(ipart) = tp(i)
        end if
     end do
     call hdf5_write_dataset_1d_dp(grp_id, 'birth_epoch', &
          dbuf, npart_loc, offset_part, npart_total)

     if(metal) then
        ipart = 0
        do i = 1, npartmax
           if(levelp(i) > 0) then
              ipart = ipart + 1
              dbuf(ipart) = zp(i)
           end if
        end do
        call hdf5_write_dataset_1d_dp(grp_id, 'metallicity', &
             dbuf, npart_loc, offset_part, npart_total)
     end if
  end if

  deallocate(dbuf, ibuf8, ibuf)
  call hdf5_close_group(grp_id)
  deallocate(npart_all)

end subroutine backup_part_hdf5

!###########################################################################
! Sinks: rank 0 only write
!###########################################################################
subroutine backup_sink_hdf5()
  use amr_commons
  use pm_commons
  use ramses_hdf5_io
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer :: idim, ilevel
  integer(HID_T) :: grp_id
  real(dp), allocatable :: dbuf(:)
  integer, allocatable :: ibuf(:)
  character(len=10) :: dstr

  call hdf5_create_group('/sinks', grp_id)
  call hdf5_write_attr_int(grp_id, 'nsink', nsink)
  call hdf5_write_attr_int(grp_id, 'nindsink', nindsink)
  call hdf5_write_attr_int(grp_id, 'levelmin', levelmin)
  call hdf5_write_attr_int(grp_id, 'nlevelmax', nlevelmax)

  if(nsink > 0) then
     allocate(ibuf(nsink))
     ibuf(1:nsink) = idsink(1:nsink)
     call hdf5_write_dataset_serial_int(grp_id, 'idsink', ibuf, nsink, myid)
     deallocate(ibuf)

     call hdf5_write_dataset_serial_dp(grp_id, 'msink', msink, nsink, myid)

     allocate(dbuf(nsink))
     ! Position
     do idim = 1, ndim
        dbuf(1:nsink) = xsink(1:nsink, idim)
        write(dstr, '(I0)') idim
        call hdf5_write_dataset_serial_dp(grp_id, 'xsink_'//trim(dstr), dbuf, nsink, myid)
     end do
     ! Velocity
     do idim = 1, ndim
        dbuf(1:nsink) = vsink(1:nsink, idim)
        write(dstr, '(I0)') idim
        call hdf5_write_dataset_serial_dp(grp_id, 'vsink_'//trim(dstr), dbuf, nsink, myid)
     end do
     ! Time
     call hdf5_write_dataset_serial_dp(grp_id, 'tsink', tsink, nsink, myid)
     ! Accretion
     dbuf(1:nsink) = dMsmbh(1:nsink)
     call hdf5_write_dataset_serial_dp(grp_id, 'dMsmbh', dbuf, nsink, myid)
     dbuf(1:nsink) = dMBH_coarse(1:nsink)
     call hdf5_write_dataset_serial_dp(grp_id, 'dMBH_coarse', dbuf, nsink, myid)
     dbuf(1:nsink) = dMEd_coarse(1:nsink)
     call hdf5_write_dataset_serial_dp(grp_id, 'dMEd_coarse', dbuf, nsink, myid)
     ! Esave
     dbuf(1:nsink) = Esave(1:nsink)
     call hdf5_write_dataset_serial_dp(grp_id, 'Esave', dbuf, nsink, myid)
     ! Gas angular momentum (jsink)
     do idim = 1, ndim
        dbuf(1:nsink) = jsink(1:nsink, idim)
        write(dstr, '(I0)') idim
        call hdf5_write_dataset_serial_dp(grp_id, 'jsink_'//trim(dstr), dbuf, nsink, myid)
     end do
     ! BH spin axis
     do idim = 1, ndim
        dbuf(1:nsink) = bhspin(1:nsink, idim)
        write(dstr, '(I0)') idim
        call hdf5_write_dataset_serial_dp(grp_id, 'bhspin_'//trim(dstr), dbuf, nsink, myid)
     end do
     ! BH spin magnitude
     dbuf(1:nsink) = spinmag(1:nsink)
     call hdf5_write_dataset_serial_dp(grp_id, 'spinmag', dbuf, nsink, myid)
     ! BH efficiency
     dbuf(1:nsink) = eps_sink(1:nsink)
     call hdf5_write_dataset_serial_dp(grp_id, 'eps_sink', dbuf, nsink, myid)
     ! Sink statistics
     do idim = 1, ndim*2+1
        do ilevel = levelmin, nlevelmax
           dbuf(1:nsink) = sink_stat(1:nsink, ilevel, idim)
           write(dstr, '(I0,"_",I0)') idim, ilevel
           call hdf5_write_dataset_serial_dp(grp_id, 'sink_stat_'//trim(dstr), dbuf, nsink, myid)
        end do
     end do
     deallocate(dbuf)
  end if

  call hdf5_close_group(grp_id)

end subroutine backup_sink_hdf5

#endif
