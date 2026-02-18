!###########################################################################
! restore_hdf5.f90 — HDF5 parallel input for cuRAMSES
!
! Reads from output_NNNNN/data_NNNNN.h5 and restores:
!   - AMR tree (grid positions, son_flag, cpu_map)
!   - Hydro (uold conservative variables)
!   - Gravity (phi, f)
!   - Particles (xp, vp, mp, idp, levelp, tp, zp)
!   - Sinks
!
! Supports variable-ncpu restart (ncpu_file != ncpu) with ksection ordering:
!   - Builds uniform ksection tree for new ncpu
!   - Recomputes cpu_map via cmp_ksection_cpumap
!   - ALL ranks read ALL data, scatter to locally owned grids
!   - First coarse step forces load_balance for optimal rebalancing
!###########################################################################
#ifdef HDF5

!###########################################################################
! AMR restore from HDF5
!###########################################################################
subroutine restore_amr_hdf5()
  use amr_commons
  use hydro_commons
  use ksection
  use morton_keys
  use morton_hash
  use ramses_hdf5_io
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer :: ilevel, i, igrid, ind, iskip, idim, info
  integer :: nlevelmax_file
  integer(i8b) :: ngrid_total
  integer, allocatable :: ngrid_per_cpu(:)
  integer, allocatable :: son_flag_buf(:), cpu_map_buf(:)
  real(dp), allocatable :: xg_all(:,:)
  integer(HID_T) :: grp_id, lvl_grp_id
  character(len=200) :: h5filename
  character(len=40) :: grp_name
  character(len=10) :: lvl_str
  character(len=5) :: nchar
  character(len=128) :: ordering_file

  ! Grid creation variables
  integer :: igrid_new, igrid_prev_cpu, igrid_father, ind_cell
  integer :: ix, iy, iz, ix_p, iy_p, iz_p, nxny
  integer :: ngrid_all_int, icpu, grid_offset
  integer(mkb) :: mkey
  real(dp) :: twotol

  ! Variable-ncpu specific variables
  integer :: father_cell
  real(dp) :: scale, dx, xx_cell(1,3), xc_off(3)
  integer :: c_tmp(1)
  integer :: nx_loc

  call title(nrestart, nchar)
  h5filename = 'output_'//trim(nchar)//'/data_'//trim(nchar)//'.h5'

  if(myid==1) write(*,*) 'HDF5 restore AMR: ', trim(h5filename)

  call hdf5_open_parallel(h5filename, MPI_COMM_WORLD)

  nxny = nx * ny

  !=====================================================
  ! Step 1: Read header (time, cosmology, nstep, etc.)
  !=====================================================
  call hdf5_open_group('/header', grp_id)
  call hdf5_read_attr_int(grp_id, 'ncpu', ncpu_file)
  call hdf5_read_attr_dp(grp_id, 'time', t)
  call hdf5_read_attr_int(grp_id, 'nstep', nstep)
  call hdf5_read_attr_int(grp_id, 'nstep_coarse', nstep_coarse)
  ! Output scheduling (noutput/tout/aout) is set by the namelist.
  ! Only restore the output counters (iout/ifout) from the file.
  ! Read noutput_file to skip the stored arrays.
  block
    integer :: noutput_file
    real(dp), allocatable :: dummy_dp(:)
    call hdf5_read_attr_int(grp_id, 'noutput', noutput_file)
    call hdf5_read_attr_int(grp_id, 'iout', iout)
    call hdf5_read_attr_int(grp_id, 'ifout', ifout)
    allocate(dummy_dp(noutput_file))
    call hdf5_read_attr_1d_dp(grp_id, 'tout', dummy_dp, noutput_file)
    call hdf5_read_attr_1d_dp(grp_id, 'aout', dummy_dp, noutput_file)
    deallocate(dummy_dp)
  end block
  call hdf5_read_attr_1d_dp(grp_id, 'dtold', dtold, nlevelmax)
  call hdf5_read_attr_1d_dp(grp_id, 'dtnew', dtnew, nlevelmax)
  call hdf5_read_attr_dp(grp_id, 'const', const)
  call hdf5_read_attr_dp(grp_id, 'mass_tot_0', mass_tot_0)
  call hdf5_read_attr_dp(grp_id, 'rho_tot', rho_tot)
  call hdf5_read_attr_dp(grp_id, 'omega_m', omega_m)
  call hdf5_read_attr_dp(grp_id, 'omega_l', omega_l)
  call hdf5_read_attr_dp(grp_id, 'omega_k', omega_k)
  call hdf5_read_attr_dp(grp_id, 'omega_b', omega_b)
  call hdf5_read_attr_dp(grp_id, 'h0', h0)
  call hdf5_read_attr_dp(grp_id, 'aexp_ini', aexp_ini)
  call hdf5_read_attr_dp(grp_id, 'boxlen_ini', boxlen_ini)
  call hdf5_read_attr_dp(grp_id, 'aexp', aexp)
  call hdf5_read_attr_dp(grp_id, 'hexp', hexp)
  call hdf5_read_attr_dp(grp_id, 'aexp_old', aexp_old)
  call hdf5_read_attr_dp(grp_id, 'epot_tot_int', epot_tot_int)
  call hdf5_read_attr_dp(grp_id, 'epot_tot_old', epot_tot_old)
  call hdf5_read_attr_dp(grp_id, 'mass_sph', mass_sph)
  call hdf5_read_attr_string(grp_id, 'ordering', ordering_file)
  call hdf5_close_group(grp_id)

  nstep_coarse_old = nstep_coarse

  if(myid==1) write(*,*) 'Restarting at t=', t, ' nstep_coarse=', nstep_coarse

  !=====================================================
  ! Step 2: Check ncpu_file vs ncpu
  !=====================================================
  if(ncpu_file /= ncpu) then
     if(ordering /= 'ksection') then
        if(myid==1) then
           write(*,*) 'ERROR: Variable-ncpu restart only supported with ksection ordering.'
           write(*,*) '       File ncpu=', ncpu_file, ', current ncpu=', ncpu
        end if
        call clean_stop
     end if
     varcpu_restart = .true.
     varcpu_restart_done = .true.   ! flag in amr_commons: force load_balance on first step
     if(myid==1) then
        write(*,*) '============================================================'
        write(*,*) '============================================================'
        write(*,*) '=====  VARIABLE NCPU RESTART                           ====='
        write(*,*) '=====  File ncpu = ', ncpu_file
        write(*,*) '=====  New  ncpu = ', ncpu
        write(*,*) '=====  Rebuilding domain decomposition for new ncpu    ====='
        write(*,*) '=====  Load balance will be forced on first coarse step ====='
        write(*,*) '============================================================'
        write(*,*) '============================================================'
     end if
     ! Allocate mapping arrays
     allocate(varcpu_grid_file_idx(1:ngridmax))
     varcpu_grid_file_idx = 0
     allocate(varcpu_ngrid_file(1:nlevelmax))
     varcpu_ngrid_file = 0
  else
     varcpu_restart = .false.
  end if

  !=====================================================
  ! Step 2b: Read coarse grid (all ranks read full coarse)
  !=====================================================
  call hdf5_open_group('/coarse', grp_id)
  call hdf5_read_dataset_all_int(grp_id, 'son', son(1:ncoarse), ncoarse)
  call hdf5_read_dataset_all_int(grp_id, 'cpu_map', cpu_map(1:ncoarse), ncoarse)
  ! Try to read flag1 (new format); fall back to computing from son
  block
     integer :: h5err
     logical :: dset_exists
     call h5lexists_f(grp_id, 'flag1', dset_exists, h5err)
     if(dset_exists .and. h5err == 0) then
        call hdf5_read_dataset_all_int(grp_id, 'flag1', flag1(1:ncoarse), ncoarse)
     else
        do i = 1, ncoarse
           if(son(i) > 0) then
              flag1(i) = 1
           else
              flag1(i) = 0
           end if
        end do
     end if
  end block
  call hdf5_close_group(grp_id)

  ! Reset coarse son — will be set when level 1 grids are created
  son(1:ncoarse) = 0

  !=====================================================
  ! Step 3: Domain decomposition
  !=====================================================
  if(varcpu_restart) then
     !--- Variable ncpu: rebuild ksection tree for new ncpu ---
     ! Arrays already allocated in init_amr for ncpu
     ! Build uniform tree (volume-balanced, no load data needed)
     nx_loc = icoarse_max - icoarse_min + 1
     scale = boxlen / dble(nx_loc)
     call build_ksection(update=.false.)
     ! Rebuild tree navigation
     call rebuild_ksec_cpuranges()
     call compute_ksec_cpu_path()
     bisec_cpubox_min2 = bisec_cpubox_min
     bisec_cpubox_max2 = bisec_cpubox_max
     if(myid==1) write(*,*) 'HDF5: ksection tree rebuilt for ncpu=', ncpu

     ! Recompute coarse cpu_map via ksection tree
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
     if(myid==1) write(*,*) 'HDF5: coarse cpu_map recomputed'
  else
     !--- Same ncpu: read ksection tree from file ---
     call hdf5_open_group('/domain', grp_id)
     if(ordering=='ksection') then
        call hdf5_read_attr_int(grp_id, 'ksec_kmax', ksec_kmax)
        call hdf5_read_attr_int(grp_id, 'ksec_nbinodes', ksec_nbinodes)
        call hdf5_read_dataset_all_int(grp_id, 'ksec_factor', ksec_factor, nksec_levels)
        call hdf5_read_dataset_all_int(grp_id, 'ksec_dir', ksec_dir, nksec_levels)
        call hdf5_read_dataset_all_int(grp_id, 'ksec_indx', ksec_indx, ksec_nbinodes)
        ! Read flattened 2D arrays
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
           call hdf5_read_dataset_all_dp(grp_id, 'ksec_wall', wall_flat, nw_total)
           call hdf5_read_dataset_all_int(grp_id, 'ksec_next', next_flat, nn_total)
           ksec_wall(1:ksec_nbinodes, 1:nw_cols) = &
                reshape(wall_flat, (/ksec_nbinodes, nw_cols/))
           ksec_next(1:ksec_nbinodes, 1:nn_cols) = &
                reshape(next_flat, (/ksec_nbinodes, nn_cols/))
           deallocate(wall_flat, next_flat)
        end block
        ! Read CPU bounding boxes
        block
           real(dp), allocatable :: cpubox_flat(:)
           integer :: ncb
           ncb = ncpu * ndim
           allocate(cpubox_flat(ncb))
           call hdf5_read_dataset_all_dp(grp_id, 'bisec_cpubox_min', cpubox_flat, ncb)
           bisec_cpubox_min(1:ncpu, 1:ndim) = reshape(cpubox_flat, (/ncpu, ndim/))
           call hdf5_read_dataset_all_dp(grp_id, 'bisec_cpubox_max', cpubox_flat, ncb)
           bisec_cpubox_max(1:ncpu, 1:ndim) = reshape(cpubox_flat, (/ncpu, ndim/))
           deallocate(cpubox_flat)
        end block
        bisec_cpubox_min2 = bisec_cpubox_min
        bisec_cpubox_max2 = bisec_cpubox_max
        ! Rebuild tree navigation arrays
        call rebuild_ksec_cpuranges()
        call compute_ksec_cpu_path()
        if(myid==1) write(*,*) 'HDF5: ksection tree restored'
     end if
     call hdf5_close_group(grp_id)
  end if

  ! Set cpu_map2 = cpu_map for coarse
  cpu_map2(1:ncoarse) = cpu_map(1:ncoarse)

  !=====================================================
  ! Step 4: Read nlevelmax_file from /amr group
  !=====================================================
  call hdf5_open_group('/amr', grp_id)
  call hdf5_read_attr_int(grp_id, 'nlevelmax_file', nlevelmax_file)
  call hdf5_close_group(grp_id)
  if(myid==1) write(*,*) 'HDF5: nlevelmax_file =', nlevelmax_file

  !=====================================================
  ! Step 5: Initialize Morton hash tables
  !=====================================================
  if(.not. allocated(mort_table)) allocate(mort_table(1:nlevelmax))
  if(.not. allocated(grid_level)) then
     allocate(grid_level(1:ngridmax))
     grid_level = 0
  end if
  ! Initialize empty hash tables for all levels
  do ilevel = 1, nlevelmax
     call morton_hash_init(mort_table(ilevel), 16)
  end do

  !=====================================================
  ! Step 6: Per-level grid creation
  !=====================================================
  ! Suppress HDF5 error messages for missing level groups
  call hdf5_suppress_errors()

  ! Compute scale for cell coordinate computation
  nx_loc = icoarse_max - icoarse_min + 1
  scale = boxlen / dble(nx_loc)

  if(varcpu_restart) then
     !===================================================
     ! VARIABLE NCPU PATH: distributed active-only grid creation
     ! Each rank creates only its active grids (cpu_map(father)==myid)
     ! Virtual grids are created via refine_coarse/refine_fine
     ! This reduces per-rank memory by ~ncpu× vs reading ALL grids
     !===================================================
     balance = .true.
     shrink = .false.

     do ilevel = 1, nlevelmax_file
        ! Check if level group exists
        write(lvl_str, '(I0)') ilevel
        grp_name = '/amr/level_'//trim(lvl_str)

        block
           integer(HID_T) :: test_grp_id
           integer :: h5err
           logical :: grp_exists
           call h5lexists_f(hdf5_file_id, trim(grp_name), grp_exists, h5err)
           if(.not. grp_exists .or. h5err /= 0) cycle
           call h5gopen_f(hdf5_file_id, trim(grp_name), test_grp_id, h5err)
           if(h5err /= 0) cycle
           lvl_grp_id = test_grp_id
        end block

        ! Read ngrid_per_cpu from file (ncpu_file entries)
        allocate(ngrid_per_cpu(ncpu_file))
        call hdf5_read_dataset_all_int(lvl_grp_id, 'ngrid_per_cpu', ngrid_per_cpu, ncpu_file)

        ! Compute total grids at this level
        ngrid_total = 0
        do i = 1, ncpu_file
           ngrid_total = ngrid_total + ngrid_per_cpu(i)
        end do
        ngrid_all_int = int(ngrid_total)
        varcpu_ngrid_file(ilevel) = ngrid_all_int

        if(myid==1) write(*,'(A,I3,A,I10)') &
             ' HDF5 level ', ilevel, ' total grids: ', ngrid_total

        deallocate(ngrid_per_cpu)

        if(ngrid_all_int == 0) then
           call hdf5_close_group(lvl_grp_id)
           cycle
        end if

        ! ALL ranks read ALL grids' data at this level (temporary buffer)
        allocate(xg_all(ngrid_all_int, ndim))
        allocate(son_flag_buf(ngrid_all_int * twotondim))

        ! Read ALL grid positions
        do idim = 1, ndim
           write(lvl_str, '(I0)') idim
           call hdf5_read_dataset_all_dp(lvl_grp_id, 'xg_'//trim(lvl_str), &
                xg_all(:, idim), ngrid_all_int)
        end do

        ! Read ALL son flags (for flag1)
        call hdf5_read_dataset_all_int(lvl_grp_id, 'son_flag', &
             son_flag_buf, ngrid_all_int * twotondim)

        ! Initialize hash table for this level (sized for local grids)
        call morton_hash_init(mort_table(ilevel), &
             max(4 * (ngrid_all_int / max(ncpu, 1) + 1), 16))

        ! Subcell offsets at this level
        dx = 0.5d0**ilevel

        !--- Phase 1: Create ACTIVE grids only (cpu_map(father)==myid) ---
        igrid_prev_cpu = 0
        do i = 1, ngrid_all_int
           ! Compute father cell index
           if(ilevel == 1) then
              twotol = 1.0d0
              ix = int(xg_all(i, 1) * twotol)
              iy = int(xg_all(i, 2) * twotol)
              iz = int(xg_all(i, 3) * twotol)
              father_cell = 1 + ix + iy * nx + iz * nxny
           else
              twotol = 2.0d0**(ilevel-1)
              ix = int(xg_all(i, 1) * twotol)
              iy = int(xg_all(i, 2) * twotol)
              iz = int(xg_all(i, 3) * twotol)
              ix_p = ix / 2
              iy_p = iy / 2
              iz_p = iz / 2
              mkey = morton_encode(ix_p, iy_p, iz_p)
              igrid_father = morton_hash_lookup(mort_table(ilevel-1), mkey)
              if(igrid_father == 0) cycle  ! Father not on this rank → not active
              ind_cell = 1 + mod(ix, 2) + 2 * mod(iy, 2) + 4 * mod(iz, 2)
              father_cell = ncoarse + (ind_cell - 1) * ngridmax + igrid_father
           end if

           ! Check if this grid is active (owned by myid)
           if(cpu_map(father_cell) /= myid) cycle

           ! Allocate grid from free list
           igrid_new = headf
           if(igrid_new == 0) then
              write(*,*) 'ERROR: out of free grids, myid=', myid, &
                   ' level=', ilevel, ' grid=', i
              call clean_stop
           end if
           headf = next(igrid_new)
           if(headf > 0) prev(headf) = 0
           numbf = numbf - 1

           ! Set xg
           xg(igrid_new, 1:ndim) = xg_all(i, 1:ndim)

           ! Set son_flag → flag1, and init son to 0
           do iskip = 1, twotondim
              ind_cell = ncoarse + (iskip - 1) * ngridmax + igrid_new
              son(ind_cell) = 0
              flag1(ind_cell) = son_flag_buf((i-1)*twotondim + iskip)
           end do

           ! Compute cpu_map for each cell via cmp_ksection_cpumap
           do iskip = 1, twotondim
              ind_cell = ncoarse + (iskip - 1) * ngridmax + igrid_new
              iz = (iskip - 1) / 4
              iy = (iskip - 1 - 4 * iz) / 2
              ix = (iskip - 1 - 2 * iy - 4 * iz)
              xx_cell(1,1) = (xg_all(i, 1) + (dble(ix) - 0.5d0) * dx &
                   - dble(icoarse_min)) * scale
              xx_cell(1,2) = (xg_all(i, 2) + (dble(iy) - 0.5d0) * dx &
                   - dble(jcoarse_min)) * scale
              xx_cell(1,3) = (xg_all(i, 3) + (dble(iz) - 0.5d0) * dx &
                   - dble(kcoarse_min)) * scale
              call cmp_ksection_cpumap(xx_cell, c_tmp, 1)
              cpu_map(ind_cell) = c_tmp(1)
              cpu_map2(ind_cell) = c_tmp(1)
           end do

           ! Set father pointer
           father(igrid_new) = father_cell

           ! Set son in parent cell
           son(father_cell) = igrid_new

           ! Insert into Morton hash
           mkey = grid_to_morton(igrid_new, ilevel)
           call morton_hash_insert(mort_table(ilevel), mkey, igrid_new)
           grid_level(igrid_new) = ilevel

           ! Append to linked list for myid at this level
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

           ! Save file→igrid mapping for hydro/poisson restore
           varcpu_grid_file_idx(igrid_new) = i
        end do

        deallocate(xg_all, son_flag_buf)
        call hdf5_close_group(lvl_grp_id)

        !--- Phase 2: Create virtual grids via RAMSES refine mechanism ---
        ! Active grids already set son(father)>0, so refine skips those cells
        ! Only cells with flag1==1 AND son==0 get new (virtual) grids
        if(ilevel == 1) then
           call flag_coarse
           call refine_coarse
        else
           call refine_fine(ilevel - 1)
        end if

        ! Build communicators for this level
        call build_comm(ilevel)

        ! Exchange flag1/cpu_map/cpu_map2 from active to virtual grids
        ! (needed for refine_fine at next level)
        call make_virtual_fine_int(flag1(1), ilevel)
        call make_virtual_fine_int(cpu_map(1), ilevel)
        call make_virtual_fine_int(cpu_map2(1), ilevel)

        if(myid==1) write(*,'(A,I3,A,I10,A,I10)') &
             ' HDF5 level ', ilevel, ' active: ', numbl(myid, ilevel), &
             ' total: ', numbtot(1, ilevel)
     end do

     balance = .false.
     shrink = .false.

  else
     !===================================================
     ! SAME NCPU PATH: original code (grids ordered by CPU)
     !===================================================
     allocate(ngrid_per_cpu(ncpu))

     do ilevel = 1, nlevelmax_file
        ! Check if level group exists
        write(lvl_str, '(I0)') ilevel
        grp_name = '/amr/level_'//trim(lvl_str)

        block
           integer(HID_T) :: test_grp_id
           integer :: h5err
           logical :: grp_exists
           call h5lexists_f(hdf5_file_id, trim(grp_name), grp_exists, h5err)
           if(.not. grp_exists .or. h5err /= 0) cycle
           call h5gopen_f(hdf5_file_id, trim(grp_name), test_grp_id, h5err)
           if(h5err /= 0) cycle
           lvl_grp_id = test_grp_id
        end block

        ! Read ngrid_per_cpu (all ranks read full array)
        call hdf5_read_dataset_all_int(lvl_grp_id, 'ngrid_per_cpu', ngrid_per_cpu, ncpu)

        ! Compute total grids at this level
        ngrid_total = 0
        do i = 1, ncpu
           ngrid_total = ngrid_total + ngrid_per_cpu(i)
        end do
        ngrid_all_int = int(ngrid_total)

        if(myid==1) write(*,'(A,I3,A,I10)') &
             ' HDF5 level ', ilevel, ' total grids: ', ngrid_total

        if(ngrid_all_int == 0) then
           call hdf5_close_group(lvl_grp_id)
           cycle
        end if

        ! ALL ranks read ALL grids' data at this level
        allocate(xg_all(ngrid_all_int, ndim))
        allocate(son_flag_buf(ngrid_all_int * twotondim))
        allocate(cpu_map_buf(ngrid_all_int * twotondim))

        ! Read ALL grid positions (every rank reads full datasets)
        do idim = 1, ndim
           write(lvl_str, '(I0)') idim
           call hdf5_read_dataset_all_dp(lvl_grp_id, 'xg_'//trim(lvl_str), &
                xg_all(:, idim), ngrid_all_int)
        end do

        ! Read ALL cell data
        call hdf5_read_dataset_all_int(lvl_grp_id, 'son_flag', &
             son_flag_buf, ngrid_all_int * twotondim)
        call hdf5_read_dataset_all_int(lvl_grp_id, 'cpu_map', &
             cpu_map_buf, ngrid_all_int * twotondim)

        ! Initialize hash table for this level
        call morton_hash_init(mort_table(ilevel), max(2 * ngrid_all_int, 16))

        !--- Create grid slots for ALL CPUs' grids ---
        ! Grids are ordered by CPU: CPU 1 first, then CPU 2, etc.
        grid_offset = 0
        do icpu = 1, ncpu
           igrid_prev_cpu = 0
           do i = 1, ngrid_per_cpu(icpu)
              ! Allocate grid from free list
              igrid_new = headf
              if(igrid_new == 0) then
                 write(*,*) 'ERROR: out of free grids, myid=', myid, &
                      ' level=', ilevel, ' cpu=', icpu
                 call clean_stop
              end if
              headf = next(igrid_new)
              if(headf > 0) prev(headf) = 0
              numbf = numbf - 1

              ! Set linked list for this CPU at this level
              if(igrid_prev_cpu == 0) then
                 headl(icpu, ilevel) = igrid_new
              else
                 next(igrid_prev_cpu) = igrid_new
              end if
              prev(igrid_new) = igrid_prev_cpu
              next(igrid_new) = 0
              taill(icpu, ilevel) = igrid_new
              numbl(icpu, ilevel) = numbl(icpu, ilevel) + 1

              ! Set xg, cpu_map, cpu_map2, son, flag1
              ind = grid_offset + i
              xg(igrid_new, 1:ndim) = xg_all(ind, 1:ndim)
              do iskip = 1, twotondim
                 ind_cell = ncoarse + (iskip - 1) * ngridmax + igrid_new
                 cpu_map(ind_cell) = cpu_map_buf((ind-1)*twotondim + iskip)
                 cpu_map2(ind_cell) = cpu_map_buf((ind-1)*twotondim + iskip)
                 son(ind_cell) = 0
                 flag1(ind_cell) = son_flag_buf((ind-1)*twotondim + iskip)
              end do

              ! Set father pointer
              if(ilevel == 1) then
                 twotol = 1.0d0
                 ix = int(xg_all(ind, 1) * twotol)
                 iy = int(xg_all(ind, 2) * twotol)
                 iz = int(xg_all(ind, 3) * twotol)
                 father(igrid_new) = 1 + ix + iy * nx + iz * nxny
              else
                 twotol = 2.0d0**(ilevel-1)
                 ix = int(xg_all(ind, 1) * twotol)
                 iy = int(xg_all(ind, 2) * twotol)
                 iz = int(xg_all(ind, 3) * twotol)
                 ix_p = ix / 2
                 iy_p = iy / 2
                 iz_p = iz / 2
                 mkey = morton_encode(ix_p, iy_p, iz_p)
                 igrid_father = morton_hash_lookup(mort_table(ilevel-1), mkey)
                 if(igrid_father == 0) then
                    write(*,*) 'ERROR: father not found! myid=', myid, &
                         ' level=', ilevel, ' cpu=', icpu, ' grid=', i, &
                         ' ix_p=', ix_p, ' iy_p=', iy_p, ' iz_p=', iz_p
                    call clean_stop
                 end if
                 ind_cell = 1 + mod(ix, 2) + 2 * mod(iy, 2) + 4 * mod(iz, 2)
                 father(igrid_new) = ncoarse + (ind_cell - 1) * ngridmax + igrid_father
              end if

              ! Set son in parent cell
              son(father(igrid_new)) = igrid_new

              ! Insert into Morton hash
              mkey = grid_to_morton(igrid_new, ilevel)
              call morton_hash_insert(mort_table(ilevel), mkey, igrid_new)
              grid_level(igrid_new) = ilevel

              igrid_prev_cpu = igrid_new
           end do
           grid_offset = grid_offset + ngrid_per_cpu(icpu)
        end do

        deallocate(xg_all, son_flag_buf, cpu_map_buf)
        call hdf5_close_group(lvl_grp_id)
     end do

     deallocate(ngrid_per_cpu)
  end if

  ! Restore HDF5 error printing
  call hdf5_restore_errors()

  !=====================================================
  ! Step 7: Set ngrid_current, update used_mem
  !=====================================================
  ngrid_current = ngridmax - numbf
  used_mem = ngrid_current
  if(myid==1) write(*,*) 'HDF5: ngrid_current =', ngrid_current

  !=====================================================
  ! Step 8: Compute numbtot via MPI_ALLREDUCE
  !=====================================================
  do ilevel = 1, nlevelmax
     numbtot(1:4, ilevel) = 0
     ! numbtot(1,ilevel) = total active grids across all CPUs
     call MPI_ALLREDUCE(numbl(myid, ilevel), numbtot(1, ilevel), 1, &
          MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, info)
  end do

  !=====================================================
  ! Step 9: Build communicators for each level
  !=====================================================
  do ilevel = 1, nlevelmax
     call build_comm(ilevel)
  end do
  if(myid==1) write(*,*) 'HDF5: build_comm completed for all levels'

  !=====================================================
  ! Step 10: Close HDF5 file
  !=====================================================
  call hdf5_close_file()

  if(myid==1) write(*,*) 'HDF5 AMR restore done.'

end subroutine restore_amr_hdf5

!###########################################################################
! Hydro restore from HDF5
!###########################################################################
subroutine restore_hydro_hdf5()
  use amr_commons
  use hydro_commons
  use ramses_hdf5_io
  ! varcpu variables now in amr_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer :: ilevel, i, igrid, ind, iskip, ivar, info
  integer :: ngrid_loc, nvar_file, fidx
  integer, allocatable :: ngrid_all(:)
  integer(i8b) :: ncells_total, offset_cells, ngrid_total
  real(dp), allocatable :: ubuf(:), ubuf_all(:)
  integer(HID_T) :: grp_id, lvl_grp_id, hdr_grp_id
  character(len=200) :: h5filename
  character(len=40) :: grp_name
  character(len=10) :: lvl_str, var_str
  character(len=5) :: nchar

  call title(nrestart, nchar)
  h5filename = 'output_'//trim(nchar)//'/data_'//trim(nchar)//'.h5'

  if(myid==1) write(*,*) 'HDF5 restore hydro: ', trim(h5filename)

  call hdf5_open_parallel(h5filename, MPI_COMM_WORLD)

  ! Read header to get ncpu_file and nvar_file
  call hdf5_open_group('/header', hdr_grp_id)
  call hdf5_read_attr_int(hdr_grp_id, 'ncpu', ncpu_file)
  call hdf5_read_attr_int(hdr_grp_id, 'nvar', nvar_file)
  call hdf5_close_group(hdr_grp_id)

  if(nvar_file /= nvar) then
     if(myid==1) write(*,*) 'WARNING: HDF5 nvar mismatch, file=', nvar_file, ' expected=', nvar
  end if

  ! Suppress HDF5 error messages for missing level groups
  call hdf5_suppress_errors()

  if(varcpu_restart) then
     !===================================================
     ! VARIABLE NCPU PATH: ALL ranks read ALL data, scatter
     !===================================================
     do ilevel = 1, nlevelmax
        ngrid_loc = numbl(myid, ilevel)
        if(varcpu_ngrid_file(ilevel) == 0) cycle

        write(lvl_str, '(I0)') ilevel
        grp_name = '/hydro/level_'//trim(lvl_str)

        block
           integer :: h5err
           logical :: grp_exists
           call h5lexists_f(hdf5_file_id, trim(grp_name), grp_exists, h5err)
           if(.not. grp_exists .or. h5err /= 0) cycle
           call h5gopen_f(hdf5_file_id, trim(grp_name), lvl_grp_id, h5err)
           if(h5err /= 0) cycle
        end block

        ! ALL ranks read ALL hydro data for this level
        allocate(ubuf_all(varcpu_ngrid_file(ilevel) * twotondim))

        do ivar = 1, min(nvar, nvar_file)
           write(var_str, '(I0)') ivar
           call hdf5_read_dataset_all_dp(lvl_grp_id, 'uold_'//trim(var_str), &
                ubuf_all, varcpu_ngrid_file(ilevel) * twotondim)

           ! Scatter to active grids owned by this rank
           igrid = headl(myid, ilevel)
           do while(igrid > 0)
              fidx = varcpu_grid_file_idx(igrid)
              do ind = 1, twotondim
                 iskip = ncoarse + (ind - 1) * ngridmax
                 uold(igrid + iskip, ivar) = ubuf_all((fidx-1)*twotondim + ind)
              end do
              igrid = next(igrid)
           end do
        end do

        deallocate(ubuf_all)
        call hdf5_close_group(lvl_grp_id)
     end do
  else
     !===================================================
     ! SAME NCPU PATH: parallel I/O with per-rank offsets
     !===================================================
     allocate(ngrid_all(ncpu))

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

        block
           integer :: h5err
           logical :: grp_exists
           call h5lexists_f(hdf5_file_id, trim(grp_name), grp_exists, h5err)
           if(.not. grp_exists .or. h5err /= 0) cycle
           call h5gopen_f(hdf5_file_id, trim(grp_name), lvl_grp_id, h5err)
           if(h5err /= 0) cycle
        end block

        if(ngrid_loc > 0) then
           allocate(ubuf(ngrid_loc * twotondim))
        else
           allocate(ubuf(1))
        end if

        ! Read uold: raw conservative variables
        do ivar = 1, min(nvar, nvar_file)
           write(var_str, '(I0)') ivar
           call hdf5_read_dataset_1d_dp(lvl_grp_id, 'uold_'//trim(var_str), &
                ubuf, ngrid_loc * twotondim, offset_cells)

           ! Scatter to uold array
           igrid = headl(myid, ilevel)
           do i = 1, ngrid_loc
              do ind = 1, twotondim
                 iskip = ncoarse + (ind - 1) * ngridmax
                 uold(igrid + iskip, ivar) = ubuf((i-1)*twotondim + ind)
              end do
              igrid = next(igrid)
           end do
        end do

        deallocate(ubuf)
        call hdf5_close_group(lvl_grp_id)
     end do

     deallocate(ngrid_all)
  end if

  call hdf5_restore_errors()
  call hdf5_close_file()

  ! Virtual exchange to populate ghost cells
  do ilevel = 1, nlevelmax
     if(numbtot(1, ilevel) == 0) cycle
     do ivar = 1, nvar
        call make_virtual_fine_dp(uold(1,ivar), ilevel)
     end do
  end do

  if(myid==1) write(*,*) 'HDF5 hydro restore done.'

end subroutine restore_hydro_hdf5

!###########################################################################
! Poisson restore from HDF5
!###########################################################################
subroutine restore_poisson_hdf5()
  use amr_commons
  use poisson_commons
  use ramses_hdf5_io
  ! varcpu variables now in amr_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer :: ilevel, i, igrid, ind, iskip, idim, info, fidx
  integer :: ngrid_loc
  integer, allocatable :: ngrid_all(:)
  integer(i8b) :: ncells_total, offset_cells, ngrid_total
  real(dp), allocatable :: pbuf(:), pbuf_all(:)
  integer(HID_T) :: lvl_grp_id
  character(len=200) :: h5filename
  character(len=40) :: grp_name
  character(len=10) :: lvl_str, dim_str
  character(len=5) :: nchar

  call title(nrestart, nchar)
  h5filename = 'output_'//trim(nchar)//'/data_'//trim(nchar)//'.h5'

  if(myid==1) write(*,*) 'HDF5 restore poisson: ', trim(h5filename)

  call hdf5_open_parallel(h5filename, MPI_COMM_WORLD)

  ! Suppress HDF5 error messages for missing level groups
  call hdf5_suppress_errors()

  if(varcpu_restart) then
     !===================================================
     ! VARIABLE NCPU PATH: ALL ranks read ALL data, scatter
     !===================================================
     do ilevel = 1, nlevelmax
        ngrid_loc = numbl(myid, ilevel)
        if(varcpu_ngrid_file(ilevel) == 0) cycle

        write(lvl_str, '(I0)') ilevel
        grp_name = '/gravity/level_'//trim(lvl_str)

        block
           integer :: h5err
           logical :: grp_exists
           call h5lexists_f(hdf5_file_id, trim(grp_name), grp_exists, h5err)
           if(.not. grp_exists .or. h5err /= 0) cycle
           call h5gopen_f(hdf5_file_id, trim(grp_name), lvl_grp_id, h5err)
           if(h5err /= 0) cycle
        end block

        ! ALL ranks read ALL data for this level
        allocate(pbuf_all(varcpu_ngrid_file(ilevel) * twotondim))

        ! Read phi
        call hdf5_read_dataset_all_dp(lvl_grp_id, 'phi', &
             pbuf_all, varcpu_ngrid_file(ilevel) * twotondim)
        igrid = headl(myid, ilevel)
        do while(igrid > 0)
           fidx = varcpu_grid_file_idx(igrid)
           do ind = 1, twotondim
              iskip = ncoarse + (ind - 1) * ngridmax
              phi(igrid + iskip) = pbuf_all((fidx-1)*twotondim + ind)
           end do
           igrid = next(igrid)
        end do

        ! Read force
        do idim = 1, ndim
           write(dim_str, '(I0)') idim
           call hdf5_read_dataset_all_dp(lvl_grp_id, 'f_'//trim(dim_str), &
                pbuf_all, varcpu_ngrid_file(ilevel) * twotondim)
           igrid = headl(myid, ilevel)
           do while(igrid > 0)
              fidx = varcpu_grid_file_idx(igrid)
              do ind = 1, twotondim
                 iskip = ncoarse + (ind - 1) * ngridmax
                 f(igrid + iskip, idim) = pbuf_all((fidx-1)*twotondim + ind)
              end do
              igrid = next(igrid)
           end do
        end do

        deallocate(pbuf_all)
        call hdf5_close_group(lvl_grp_id)
     end do

     ! Free mapping arrays — no longer needed
     if(allocated(varcpu_grid_file_idx)) deallocate(varcpu_grid_file_idx)
     if(allocated(varcpu_ngrid_file)) deallocate(varcpu_ngrid_file)
  else
     !===================================================
     ! SAME NCPU PATH: parallel I/O with per-rank offsets
     !===================================================
     allocate(ngrid_all(ncpu))

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

        block
           integer :: h5err
           logical :: grp_exists
           call h5lexists_f(hdf5_file_id, trim(grp_name), grp_exists, h5err)
           if(.not. grp_exists .or. h5err /= 0) cycle
           call h5gopen_f(hdf5_file_id, trim(grp_name), lvl_grp_id, h5err)
           if(h5err /= 0) cycle
        end block

        if(ngrid_loc > 0) then
           allocate(pbuf(ngrid_loc * twotondim))
        else
           allocate(pbuf(1))
        end if

        ! Read phi
        call hdf5_read_dataset_1d_dp(lvl_grp_id, 'phi', &
             pbuf, ngrid_loc * twotondim, offset_cells)
        igrid = headl(myid, ilevel)
        do i = 1, ngrid_loc
           do ind = 1, twotondim
              iskip = ncoarse + (ind - 1) * ngridmax
              phi(igrid + iskip) = pbuf((i-1)*twotondim + ind)
           end do
           igrid = next(igrid)
        end do

        ! Read force
        do idim = 1, ndim
           write(dim_str, '(I0)') idim
           call hdf5_read_dataset_1d_dp(lvl_grp_id, 'f_'//trim(dim_str), &
                pbuf, ngrid_loc * twotondim, offset_cells)
           igrid = headl(myid, ilevel)
           do i = 1, ngrid_loc
              do ind = 1, twotondim
                 iskip = ncoarse + (ind - 1) * ngridmax
                 f(igrid + iskip, idim) = pbuf((i-1)*twotondim + ind)
              end do
              igrid = next(igrid)
           end do
        end do

        deallocate(pbuf)
        call hdf5_close_group(lvl_grp_id)
     end do

     deallocate(ngrid_all)
  end if

  call hdf5_restore_errors()
  call hdf5_close_file()

  ! Virtual exchange to populate ghost cells
  do ilevel = 1, nlevelmax
     if(numbtot(1, ilevel) == 0) cycle
     call make_virtual_fine_dp(phi(1), ilevel)
     do idim = 1, ndim
        call make_virtual_fine_dp(f(1,idim), ilevel)
     end do
  end do

  if(myid==1) write(*,*) 'HDF5 poisson restore done.'

end subroutine restore_poisson_hdf5

!###########################################################################
! Particle restore from HDF5
!###########################################################################
subroutine restore_part_hdf5()
  use amr_commons
  use pm_commons
  use ramses_hdf5_io
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer :: i, idim, info
  integer :: npart_loc
  integer(i8b) :: npart_total_file, offset_part, tmp_long
  integer, allocatable :: npart_per_cpu(:)
  real(dp), allocatable :: dbuf(:)
  integer(i8b), allocatable :: ibuf8(:)
  integer, allocatable :: ibuf(:)
  integer(HID_T) :: grp_id, hdr_grp_id
  character(len=200) :: h5filename
  character(len=10) :: dstr
  character(len=5) :: nchar

  call title(nrestart, nchar)
  h5filename = 'output_'//trim(nchar)//'/data_'//trim(nchar)//'.h5'

  if(myid==1) write(*,*) 'HDF5 restore particles: ', trim(h5filename)

  call hdf5_open_parallel(h5filename, MPI_COMM_WORLD)

  ! Read header
  call hdf5_open_group('/header', hdr_grp_id)
  call hdf5_read_attr_int(hdr_grp_id, 'ncpu', ncpu_file)
  call hdf5_close_group(hdr_grp_id)

  ! Read particle group
  call hdf5_open_group('/particles', grp_id)
  call hdf5_read_attr_int8(grp_id, 'npart_total', npart_total_file)
  block
     integer(i8b) :: nstar_tot_tmp
     call hdf5_read_attr_int8(grp_id, 'nstar_tot', nstar_tot_tmp)
     nstar_tot = nstar_tot_tmp
  end block
  call hdf5_read_attr_dp(grp_id, 'mstar_tot', mstar_tot)
  call hdf5_read_attr_dp(grp_id, 'mstar_lost', mstar_lost)

  ! Read npart_per_cpu
  allocate(npart_per_cpu(ncpu_file))
  call hdf5_read_dataset_all_int(grp_id, 'npart_per_cpu', npart_per_cpu, ncpu_file)

  ! Determine local particle count and offset
  if(ncpu_file == ncpu) then
     ! Same rank: direct mapping
     npart_loc = npart_per_cpu(myid)
     offset_part = 0
     do i = 1, myid - 1
        offset_part = offset_part + npart_per_cpu(i)
     end do
  else
     ! Different rank: distribute evenly
     block
        integer(i8b) :: npp, remainder
        npp = npart_total_file / ncpu
        remainder = npart_total_file - npp * ncpu
        if(myid <= int(remainder)) then
           npart_loc = int(npp + 1)
           offset_part = int(npp + 1, i8b) * (myid - 1)
        else
           npart_loc = int(npp)
           offset_part = (npp + 1) * remainder + npp * (myid - 1 - remainder)
        end if
     end block
  end if

  if(npart_loc > npartmax) then
     write(*,*) 'ERROR: npart_loc > npartmax', npart_loc, npartmax
     call clean_stop
  end if

  ! Allocate buffers
  allocate(dbuf(npart_loc))
  allocate(ibuf8(npart_loc))
  allocate(ibuf(npart_loc))

  ! Read positions
  do idim = 1, ndim
     write(dstr, '(I0)') idim
     call hdf5_read_dataset_1d_dp(grp_id, 'x_'//trim(dstr), dbuf, npart_loc, offset_part)
     xp(1:npart_loc, idim) = dbuf(1:npart_loc)
  end do

  ! Read velocities
  do idim = 1, ndim
     write(dstr, '(I0)') idim
     call hdf5_read_dataset_1d_dp(grp_id, 'v_'//trim(dstr), dbuf, npart_loc, offset_part)
     vp(1:npart_loc, idim) = dbuf(1:npart_loc)
  end do

  ! Read mass
  call hdf5_read_dataset_1d_dp(grp_id, 'mass', dbuf, npart_loc, offset_part)
  mp(1:npart_loc) = dbuf(1:npart_loc)

  ! Read identity
  call hdf5_read_dataset_1d_int8(grp_id, 'identity', ibuf8, npart_loc, offset_part)
  idp(1:npart_loc) = ibuf8(1:npart_loc)

  ! Read level
  call hdf5_read_dataset_1d_int(grp_id, 'levelp', ibuf, npart_loc, offset_part)
  levelp(1:npart_loc) = ibuf(1:npart_loc)

  ! Read birth epoch and metallicity if star/sink
  if(star .or. sink) then
     call hdf5_read_dataset_1d_dp(grp_id, 'birth_epoch', dbuf, npart_loc, offset_part)
     tp(1:npart_loc) = dbuf(1:npart_loc)
     if(metal) then
        call hdf5_read_dataset_1d_dp(grp_id, 'metallicity', dbuf, npart_loc, offset_part)
        zp(1:npart_loc) = dbuf(1:npart_loc)
     end if
  end if

  npart = npart_loc

  deallocate(dbuf, ibuf8, ibuf, npart_per_cpu)
  call hdf5_close_group(grp_id)

  ! Read sinks
  if(sink) then
     call hdf5_open_group('/sinks', grp_id)
     call hdf5_read_attr_int(grp_id, 'nsink', nsink)
     call hdf5_read_attr_int(grp_id, 'nindsink', nindsink)

     if(nsink > 0) then
        block
           real(dp), allocatable :: sbuf(:)
           integer, allocatable :: sibuf(:)
           integer :: ilevel
           character(len=20) :: stat_name

           allocate(sibuf(nsink))
           call hdf5_read_dataset_all_int(grp_id, 'idsink', sibuf, nsink)
           idsink(1:nsink) = sibuf(1:nsink)
           nindsink = maxval(idsink(1:nsink))
           deallocate(sibuf)

           allocate(sbuf(nsink))
           call hdf5_read_dataset_all_dp(grp_id, 'msink', sbuf, nsink)
           msink(1:nsink) = sbuf(1:nsink)

           ! Position
           do idim = 1, ndim
              write(dstr, '(I0)') idim
              call hdf5_read_dataset_all_dp(grp_id, 'xsink_'//trim(dstr), sbuf, nsink)
              xsink(1:nsink, idim) = sbuf(1:nsink)
           end do
           ! Velocity
           do idim = 1, ndim
              write(dstr, '(I0)') idim
              call hdf5_read_dataset_all_dp(grp_id, 'vsink_'//trim(dstr), sbuf, nsink)
              vsink(1:nsink, idim) = sbuf(1:nsink)
           end do
           ! Time
           call hdf5_read_dataset_all_dp(grp_id, 'tsink', sbuf, nsink)
           tsink(1:nsink) = sbuf(1:nsink)
           ! Accretion
           call hdf5_read_dataset_all_dp(grp_id, 'dMsmbh', sbuf, nsink)
           dMsmbh(1:nsink) = sbuf(1:nsink)
           call hdf5_read_dataset_all_dp(grp_id, 'dMBH_coarse', sbuf, nsink)
           dMBH_coarse(1:nsink) = sbuf(1:nsink)
           call hdf5_read_dataset_all_dp(grp_id, 'dMEd_coarse', sbuf, nsink)
           dMEd_coarse(1:nsink) = sbuf(1:nsink)
           ! Esave
           call hdf5_read_dataset_all_dp(grp_id, 'Esave', sbuf, nsink)
           Esave(1:nsink) = sbuf(1:nsink)
           ! Gas angular momentum (jsink)
           do idim = 1, ndim
              write(dstr, '(I0)') idim
              call hdf5_read_dataset_all_dp(grp_id, 'jsink_'//trim(dstr), sbuf, nsink)
              jsink(1:nsink, idim) = sbuf(1:nsink)
           end do
           ! BH spin axis
           do idim = 1, ndim
              write(dstr, '(I0)') idim
              call hdf5_read_dataset_all_dp(grp_id, 'bhspin_'//trim(dstr), sbuf, nsink)
              bhspin(1:nsink, idim) = sbuf(1:nsink)
           end do
           ! BH spin magnitude
           call hdf5_read_dataset_all_dp(grp_id, 'spinmag', sbuf, nsink)
           spinmag(1:nsink) = sbuf(1:nsink)
           ! BH efficiency
           call hdf5_read_dataset_all_dp(grp_id, 'eps_sink', sbuf, nsink)
           eps_sink(1:nsink) = sbuf(1:nsink)
           ! Sink statistics
           do idim = 1, ndim*2+1
              do ilevel = levelmin, nlevelmax
                 write(stat_name, '(I0,"_",I0)') idim, ilevel
                 call hdf5_read_dataset_all_dp(grp_id, 'sink_stat_'//trim(stat_name), sbuf, nsink)
                 sink_stat(1:nsink, ilevel, idim) = sbuf(1:nsink)
              end do
           end do
           deallocate(sbuf)
        end block
     end if

     call hdf5_close_group(grp_id)
  end if

  call hdf5_close_file()

  if(myid==1) write(*,*) 'HDF5 particle restore done.'

end subroutine restore_part_hdf5

#endif
