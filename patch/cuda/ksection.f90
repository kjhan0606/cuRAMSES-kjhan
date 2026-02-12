! Recursive k-section domain decomposition module
! Generalizes bisection (k=2) to arbitrary k per level.
! For N_rank = p^m * q^n * ... (p>q>...), splits by p (m times),
! then q (n times), etc., always along the longest axis.

module ksection
   use amr_parameters
   use amr_commons
   use bisection, only: round_to_bisec_res, &
        init_bisection_histogram, build_bisection_histogram

   implicit none

contains

   !================================================================
   ! Prime factorization: compute split sequence for ncpu
   ! Output stored in module variables: ksec_factor, nksec_levels, ksec_kmax
   !================================================================
   subroutine compute_ksection_factorization(n)
      integer, intent(in) :: n

      integer :: i, j, temp, nfactors
      integer :: primes(64), mults(64)
      integer :: itmp

      ! Factor n into primes
      temp = n
      nfactors = 0
      i = 2
      do while (i * i <= temp)
         if (mod(temp, i) == 0) then
            nfactors = nfactors + 1
            primes(nfactors) = i
            mults(nfactors) = 0
            do while (mod(temp, i) == 0)
               mults(nfactors) = mults(nfactors) + 1
               temp = temp / i
            end do
         end if
         i = i + 1
      end do
      if (temp > 1) then
         nfactors = nfactors + 1
         primes(nfactors) = temp
         mults(nfactors) = 1
      end if

      ! Sort primes descending (selection sort, nfactors is small)
      do i = 1, nfactors - 1
         do j = i + 1, nfactors
            if (primes(j) > primes(i)) then
               itmp = primes(i); primes(i) = primes(j); primes(j) = itmp
               itmp = mults(i);  mults(i)  = mults(j);  mults(j)  = itmp
            end if
         end do
      end do

      ! Count total levels
      nksec_levels = 0
      do i = 1, nfactors
         nksec_levels = nksec_levels + mults(i)
      end do

      ! Max factor
      ksec_kmax = primes(1)

      ! Build factor sequence: [p, p, ...(m times), q, q, ...(n times), ...]
      if(allocated(ksec_factor)) deallocate(ksec_factor)
      allocate(ksec_factor(1:nksec_levels))
      j = 0
      do i = 1, nfactors
         do temp = 1, mults(i)
            j = j + 1
            ksec_factor(j) = primes(i)
         end do
      end do

   end subroutine compute_ksection_factorization


   !================================================================
   ! Compute split direction per level using longest-axis selection
   !================================================================
   subroutine compute_ksection_directions(scale_val)
      real(dp), intent(in) :: scale_val

      real(dp) :: dims(1:3)
      integer :: i, dir, maxdir

      if(allocated(ksec_dir)) deallocate(ksec_dir)
      allocate(ksec_dir(1:nksec_levels))

      dims = 0.0d0
      do i = 1, ndim
         dims(i) = scale_val
      end do

      do i = 1, nksec_levels
         ! Find longest dimension
         maxdir = 1
         do dir = 2, ndim
            if (dims(dir) > dims(maxdir)) maxdir = dir
         end do
         ksec_dir(i) = maxdir
         dims(maxdir) = dims(maxdir) / dble(ksec_factor(i))
      end do

   end subroutine compute_ksection_directions


   !================================================================
   ! Compute total number of tree nodes for the k-ary tree
   !================================================================
   function compute_ksec_nbinodes() result(ntot)
      integer :: ntot
      integer :: lvl, nc

      ntot = 1   ! root
      nc = 1
      do lvl = 1, nksec_levels
         nc = nc * ksec_factor(lvl)
         ntot = ntot + nc
      end do
   end function compute_ksec_nbinodes


   !================================================================
   ! Map points to CPU ids by walking the k-ary tree
   !================================================================
   subroutine cmp_ksection_cpumap(x, c, nn)
      real(dp), intent(in), dimension(:,:) :: x
      integer, intent(out), dimension(:)   :: c
      integer, intent(in) :: nn

      integer :: p, lvl, cur, j, k, dir, child_idx

      do p = 1, nn
         cur = ksec_root
         lvl = 0
         do while (ksec_next(cur, 1) > 0)  ! not a leaf
            lvl = lvl + 1
            k   = ksec_factor(lvl)
            dir = ksec_dir(lvl)
            ! Find which partition point falls in
            child_idx = k  ! default: last partition
            do j = 1, k - 1
               if (x(p, dir) <= ksec_wall(cur, j)) then
                  child_idx = j
                  exit
               end if
            end do
            cur = ksec_next(cur, child_idx)
         end do
         c(p) = ksec_indx(cur)
      end do

   end subroutine cmp_ksection_cpumap


   !================================================================
   ! K-way histogram splitsort
   ! Generalizes splitsort_bisection_histogram from 2-way to k-way
   !================================================================
   subroutine splitsort_ksection_histogram(lev, dir, all_walls, kfac)
      integer, intent(in) :: lev, dir, kfac
      real(dp), intent(in), dimension(:,:) :: all_walls  ! (nc, 1:kfac-1)

      integer :: i, nc, j, lmost, rmost, tmp, nxny
      integer :: part, nparts_done
      real(dp) :: xc_lmost, subcell_c, dx, scale

      integer :: ix, iy, iz
      integer :: icell, igrid, isubcell, nx_loc
      integer, dimension(1:3) :: iarray, icoarse_array

      if(verbose) print *,'entering splitsort_ksection_histogram'

      iarray = 0
      icoarse_array = (/ icoarse_min, jcoarse_min, kcoarse_min /)
      nx_loc = icoarse_max - icoarse_min + 1
      scale = boxlen / dble(nx_loc)
      nxny = nx * ny

      ! Compute nc = nodes at level lev
      nc = 1
      do i = 1, lev
         nc = nc * ksec_factor(i)
      end do

      new_hist_bounds = 0

      do i = 1, nc
         ! For each partition boundary (k-1 walls), split the domain
         ! We do k-way split by iterating: for each wall, move cells
         ! past the wall to the right end, then advance
         nparts_done = 0
         lmost = bisec_hist_bounds(i)

         do part = 1, kfac - 1
            ! Split at wall part: cells < wall go left, >= wall go right
            rmost = bisec_hist_bounds(i + 1) - 1

            ! But we only sort within the remaining unsorted range
            ! lmost is the start of the current unsorted region
            ! We want to separate cells < walls(i,part) to the left
            rmost = bisec_hist_bounds(i + 1) - 1
            ! Actually, count from lmost
            ! Partition: all cells in [lmost, rmost] such that
            ! coord < walls(i,part) go to [lmost, split_point-1]
            ! coord >= walls(i,part) go to [split_point, rmost]

            ! Find the split point using in-place partitioning
            rmost = bisec_hist_bounds(i + 1) - 1
            do while (rmost - lmost > 0)
               ! Compute coordinate of leftmost cell
               icell = bisec_ind_cell(lmost)
               dx = 0.5d0**cell_level(icell)
               if (icell <= ncoarse) then
                  iz = (icell - 1) / nxny
                  iy = ((icell - 1) - iz * nxny) / nx
                  ix = ((icell - 1) - iy * nx - iz * nxny)
                  iarray = (/ ix, iy, iz /)
                  xc_lmost = scale * (dble(iarray(dir)) + 0.5d0 - dble(icoarse_array(dir)))
               else
                  isubcell = ((icell - ncoarse) / ngridmax) + 1
                  igrid = icell - ncoarse - ngridmax * (isubcell - 1)
                  iz = (isubcell - 1) / 4
                  iy = (isubcell - 1 - 4 * iz) / 2
                  ix = (isubcell - 1 - 2 * iy - 4 * iz)
                  iarray = (/ ix, iy, iz /)
                  subcell_c = (dble(iarray(dir)) - 0.5d0) * dx - dble(icoarse_array(dir))
                  xc_lmost = scale * (xg(igrid, dir) + subcell_c)
               end if

               if (xc_lmost < all_walls(i, part)) then
                  lmost = lmost + 1
               else
                  tmp = bisec_ind_cell(lmost)
                  bisec_ind_cell(lmost) = bisec_ind_cell(rmost)
                  bisec_ind_cell(rmost) = tmp
                  rmost = rmost - 1
               end if
            end do

            ! Record split point as end of partition 'part' / start of partition 'part+1'
            new_hist_bounds(kfac * (i - 1) + part + 1) = lmost
            ! lmost is already at the right position for the next partition
         end do

         ! First bound for domain i (start of partition 1)
         new_hist_bounds(kfac * (i - 1) + 1) = bisec_hist_bounds(i)
         ! Last bound for domain i (end of last partition)
         new_hist_bounds(kfac * i + 1) = bisec_hist_bounds(i + 1)
      end do

      bisec_hist_bounds = new_hist_bounds

   end subroutine splitsort_ksection_histogram


   !================================================================
   ! MAIN K-SECTION TREE CREATION/UPDATING ROUTINE
   !================================================================
   subroutine build_ksection(update)
#ifndef WITHOUTMPI
      include 'mpif.h'
#endif
      logical, intent(in) :: update

      ! Tree-wide variables
      integer,  allocatable, dimension(:) :: tmp_imin, tmp_imax
      integer(i8b),  allocatable, dimension(:) :: tmp_load
      real(dp), allocatable, dimension(:,:) :: tmp_bxmin, tmp_bxmax

      ! Level-wide variables for dichotomy (flattened: nc*(k-1) walls)
      integer :: nwalls_level, nwalls_max
      logical,  allocatable, dimension(:) :: skip_wall
      real(dp), allocatable, dimension(:) :: l_limit, u_limit
      real(dp), allocatable, dimension(:) :: last_wall_arr, best_wall_arr, best_score_arr
      real(dp), allocatable, dimension(:) :: target_cum_frac
      integer(i8b),  allocatable, dimension(:) :: myload, cum_load_global, totload_arr
      logical,  allocatable, dimension(:) :: skip_node

      ! Wall storage per level for splitsort
      real(dp), allocatable, dimension(:,:) :: walls_2d

      logical :: all_skip, start_bisec
      integer(i8b) :: mytmp, tottmp
      real(dp) :: scale, mean, var, stdev
      real(dp) :: target_load_val, actual_load_val, score

      integer :: nc, dir, i, j, lvl, ierr, iter, k, fw
      integer :: lncpu, cpuid, child
      integer :: cur_levelstart, cur_cell
      integer :: base_ncpu, rem_ncpu, cum_ncpu
      integer, allocatable, dimension(:) :: ncpu_part

      scale = boxlen / dble(icoarse_max - icoarse_min + 1)

      if(verbose) print *,'entering build_ksection with update = ', update

      ! Allocate tree-wide temporaries
      allocate(tmp_imin(1:ksec_nbinodes))
      allocate(tmp_imax(1:ksec_nbinodes))
      allocate(tmp_load(1:ksec_nbinodes))
      allocate(tmp_bxmin(1:ksec_nbinodes, 1:ndim))
      allocate(tmp_bxmax(1:ksec_nbinodes, 1:ndim))

      ! Compute max walls at any level for allocation
      nwalls_max = 0
      nc = 1
      do lvl = 1, nksec_levels
         nwalls_max = max(nwalls_max, nc * (ksec_factor(lvl) - 1))
         nc = nc * ksec_factor(lvl)
      end do

      allocate(skip_wall(1:nwalls_max))
      allocate(l_limit(1:nwalls_max))
      allocate(u_limit(1:nwalls_max))
      allocate(last_wall_arr(1:nwalls_max))
      allocate(best_wall_arr(1:nwalls_max))
      allocate(best_score_arr(1:nwalls_max))
      allocate(target_cum_frac(1:nwalls_max))
      allocate(myload(1:nwalls_max))
      allocate(cum_load_global(1:nwalls_max))
      allocate(totload_arr(1:nwalls_max))
      allocate(skip_node(1:ncpu))

      ! TREE INIT
      ksec_root = 1
      ksec_next = 0

      tmp_imin = 0; tmp_imax = 0
      tmp_imin(1) = 1; tmp_imax(1) = ncpu
      tmp_bxmin(1,:) = 0.0d0; tmp_bxmax(1,:) = scale

      cur_levelstart = 1

      if (update) then
         call init_bisection_histogram
         dir = ksec_dir(1)
         call build_bisection_histogram(0, dir, 1)
         mytmp = bisec_hist(1, bisec_nres)
#ifndef WITHOUTMPI
         call MPI_ALLREDUCE(mytmp, tottmp, 1, MPI_INTEGER8, MPI_SUM, MPI_COMM_WORLD, ierr)
#else
         tottmp = mytmp
#endif
         tmp_load(1) = tottmp
      end if

      ! Loop through levels
      nc = 1
      level_loop: do lvl = 0, nksec_levels - 1

         k   = ksec_factor(lvl + 1)
         dir = ksec_dir(lvl + 1)
         nwalls_level = nc * (k - 1)

         ! Rebuild histograms
         if (update .and. lvl > 0) call build_bisection_histogram(lvl, dir, nc)

         ! Allocate per-node partition info
         allocate(ncpu_part(1:k))
         allocate(walls_2d(1:nc, 1:max(k-1,1)))
         walls_2d = 0.0d0

         ! WALL-FINDING
         start_bisec = .true.
         skip_wall(1:nwalls_level) = .false.
         skip_node(1:nc) = .false.

         ! Check for skip conditions and compute ncpu_part
         do i = 1, nc
            cur_cell = cur_levelstart + (i - 1)
            if (tmp_imax(cur_cell) == 0) then
               skip_node(i) = .true.
               do j = 1, k - 1
                  skip_wall((i-1)*(k-1) + j) = .true.
               end do
               cycle
            end if
            lncpu = tmp_imax(cur_cell) - tmp_imin(cur_cell) + 1
            if (lncpu == 1) then
               skip_node(i) = .true.
               do j = 1, k - 1
                  skip_wall((i-1)*(k-1) + j) = .true.
               end do
               cycle
            end if
         end do

         ! BUILD FROM SCRATCH (not update)
         build_from_scratch: if (.not. update) then
            do i = 1, nc
               cur_cell = cur_levelstart + (i - 1)
               if (skip_node(i)) cycle
               lncpu = tmp_imax(cur_cell) - tmp_imin(cur_cell) + 1
               base_ncpu = lncpu / k
               rem_ncpu  = mod(lncpu, k)
               cum_ncpu = 0
               do j = 1, k - 1
                  if (j <= rem_ncpu) then
                     cum_ncpu = cum_ncpu + base_ncpu + 1
                  else
                     cum_ncpu = cum_ncpu + base_ncpu
                  end if
                  ksec_wall(cur_cell, j) = round_to_bisec_res( &
                       tmp_bxmin(cur_cell, dir) + &
                       (tmp_bxmax(cur_cell, dir) - tmp_bxmin(cur_cell, dir)) * &
                       dble(cum_ncpu) / dble(lncpu) )
                  ! Safety check
                  if (ksec_wall(cur_cell, j) <= tmp_bxmin(cur_cell, dir) .or. &
                      ksec_wall(cur_cell, j) >= tmp_bxmax(cur_cell, dir)) then
                     if(myid==1) print *, "Problem in ksection tree creation: insufficient resolution"
#ifndef WITHOUTMPI
                     call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
#endif
                     stop
                  end if
                  walls_2d(i, j) = ksec_wall(cur_cell, j)
               end do
            end do
         end if build_from_scratch

         ! UPDATE: iterative dichotomy for load balancing
         update_dichotomy: if (update) then

            ! Compute target cumulative fractions
            do i = 1, nc
               cur_cell = cur_levelstart + (i - 1)
               if (skip_node(i)) cycle
               lncpu = tmp_imax(cur_cell) - tmp_imin(cur_cell) + 1
               base_ncpu = lncpu / k
               rem_ncpu  = mod(lncpu, k)
               cum_ncpu = 0
               do j = 1, k - 1
                  fw = (i - 1) * (k - 1) + j
                  if (j <= rem_ncpu) then
                     cum_ncpu = cum_ncpu + base_ncpu + 1
                  else
                     cum_ncpu = cum_ncpu + base_ncpu
                  end if
                  target_cum_frac(fw) = dble(cum_ncpu) / dble(lncpu)
               end do
            end do

            iter = 0
            dichotomy_loop: do
               iter = iter + 1

               ! Check all done
               all_skip = .true.
               do fw = 1, nwalls_level
                  all_skip = all_skip .and. skip_wall(fw)
               end do
               if (all_skip) exit

               if (start_bisec) then
                  ! Initial: extract cumulative loads at wall positions
                  do i = 1, nc
                     cur_cell = cur_levelstart + (i - 1)
                     if (skip_node(i)) cycle
                     do j = 1, k - 1
                        fw = (i - 1) * (k - 1) + j
                        ! Check wall position is within bounds
                        if (ksec_wall(cur_cell, j) <= tmp_bxmin(cur_cell, dir) .or. &
                            ksec_wall(cur_cell, j) >= tmp_bxmax(cur_cell, dir)) then
                           ksec_wall(cur_cell, j) = round_to_bisec_res( &
                                tmp_bxmin(cur_cell, dir) + &
                                target_cum_frac(fw) * &
                                (tmp_bxmax(cur_cell, dir) - tmp_bxmin(cur_cell, dir)))
                        end if
                        myload(fw) = bisec_hist(i, floor(ksec_wall(cur_cell, j) / bisec_res) + 1)
                     end do
                  end do

#ifndef WITHOUTMPI
                  call MPI_ALLREDUCE(myload, cum_load_global, nwalls_level, &
                       MPI_INTEGER8, MPI_SUM, MPI_COMM_WORLD, ierr)
#else
                  cum_load_global(1:nwalls_level) = myload(1:nwalls_level)
#endif

                  ! Initialize search bounds and best values
                  do i = 1, nc
                     cur_cell = cur_levelstart + (i - 1)
                     if (skip_node(i)) cycle
                     do j = 1, k - 1
                        fw = (i - 1) * (k - 1) + j
                        best_wall_arr(fw) = ksec_wall(cur_cell, j)
#if NPRE==4
                        best_score_arr(fw) = huge(1.0e0)
#else
                        best_score_arr(fw) = huge(1.0d0)
#endif
                        l_limit(fw) = tmp_bxmin(cur_cell, dir)
                        u_limit(fw) = tmp_bxmax(cur_cell, dir)
                     end do
                  end do
                  start_bisec = .false.

               else
                  ! Differential loads
                  do i = 1, nc
                     cur_cell = cur_levelstart + (i - 1)
                     if (skip_node(i)) cycle
                     do j = 1, k - 1
                        fw = (i - 1) * (k - 1) + j
                        if (skip_wall(fw)) cycle
                        myload(fw) = abs( &
                             bisec_hist(i, floor(max(ksec_wall(cur_cell,j), last_wall_arr(fw)) / bisec_res) + 1) &
                           - bisec_hist(i, floor(min(ksec_wall(cur_cell,j), last_wall_arr(fw)) / bisec_res) + 1))
                     end do
                  end do

#ifndef WITHOUTMPI
                  call MPI_ALLREDUCE(myload, totload_arr, nwalls_level, &
                       MPI_INTEGER8, MPI_SUM, MPI_COMM_WORLD, ierr)
#else
                  totload_arr(1:nwalls_level) = myload(1:nwalls_level)
#endif

                  ! Update cumulative loads
                  do i = 1, nc
                     cur_cell = cur_levelstart + (i - 1)
                     if (skip_node(i)) cycle
                     do j = 1, k - 1
                        fw = (i - 1) * (k - 1) + j
                        if (skip_wall(fw)) cycle
                        if (ksec_wall(cur_cell, j) > last_wall_arr(fw)) then
                           cum_load_global(fw) = cum_load_global(fw) + totload_arr(fw)
                        else
                           cum_load_global(fw) = cum_load_global(fw) - totload_arr(fw)
                        end if
                     end do
                  end do
               end if

               ! Score computation and convergence
               do i = 1, nc
                  cur_cell = cur_levelstart + (i - 1)
                  if (skip_node(i)) cycle
                  do j = 1, k - 1
                     fw = (i - 1) * (k - 1) + j
                     if (skip_wall(fw)) cycle

                     ! Target and actual cumulative loads
                     target_load_val = target_cum_frac(fw) * dble(tmp_load(cur_cell))
                     actual_load_val = dble(cum_load_global(fw))

                     ! Score: relative imbalance
                     if (target_load_val > 0.0d0) then
                        score = abs(actual_load_val - target_load_val) / target_load_val
                     else
                        score = 0.0d0
                     end if

                     ! Tolerance met?
                     if (score < bisec_tol) then
                        skip_wall(fw) = .true.
                        cycle
                     end if

                     ! Store best wall
                     if (score < best_score_arr(fw)) then
                        best_score_arr(fw) = score
                        best_wall_arr(fw) = ksec_wall(cur_cell, j)
                     end if

                     ! Adjust wall position
                     if (actual_load_val > target_load_val) then
                        u_limit(fw) = ksec_wall(cur_cell, j)
                     else
                        l_limit(fw) = ksec_wall(cur_cell, j)
                     end if

                     last_wall_arr(fw) = ksec_wall(cur_cell, j)
                     ksec_wall(cur_cell, j) = round_to_bisec_res(0.5d0 * (u_limit(fw) + l_limit(fw)))

                     ! Check resolution limit
                     if (abs(ksec_wall(cur_cell, j) - u_limit(fw)) < 0.5d0 * bisec_res .or. &
                         abs(ksec_wall(cur_cell, j) - l_limit(fw)) < 0.5d0 * bisec_res) then
                        ksec_wall(cur_cell, j) = best_wall_arr(fw)
                        skip_wall(fw) = .true.
                     end if
                  end do
               end do

            end do dichotomy_loop

            ! Copy final walls for splitsort
            do i = 1, nc
               cur_cell = cur_levelstart + (i - 1)
               if (skip_node(i)) cycle
               do j = 1, k - 1
                  walls_2d(i, j) = ksec_wall(cur_cell, j)
               end do
            end do

         end if update_dichotomy

         ! CHILDREN CREATION AND LEAF PROCESSING
         do i = 1, nc
            cur_cell = cur_levelstart + (i - 1)
            if (tmp_imax(cur_cell) == 0) cycle

            ! Is current cell a leaf?
            if (tmp_imin(cur_cell) == tmp_imax(cur_cell)) then
               cpuid = tmp_imin(cur_cell)
               ksec_indx(cur_cell) = cpuid
               if (update) then
                  bisec_cpubox_min2(cpuid, :) = tmp_bxmin(cur_cell, :)
                  bisec_cpubox_max2(cpuid, :) = tmp_bxmax(cur_cell, :)
               else
                  bisec_cpubox_min(cpuid, :) = tmp_bxmin(cur_cell, :)
                  bisec_cpubox_max(cpuid, :) = tmp_bxmax(cur_cell, :)
               end if
               bisec_cpu_load(cpuid) = tmp_load(cur_cell)
               ksec_next(cur_cell, :) = 0
               cycle
            end if

            ! Create k children
            lncpu = tmp_imax(cur_cell) - tmp_imin(cur_cell) + 1
            base_ncpu = lncpu / k
            rem_ncpu  = mod(lncpu, k)
            cum_ncpu = 0

            do j = 1, k
               child = cur_levelstart + nc + (i - 1) * k + (j - 1)
               ksec_next(cur_cell, j) = child

               ! Bounding box: copy parent, then adjust split direction
               tmp_bxmin(child, :) = tmp_bxmin(cur_cell, :)
               tmp_bxmax(child, :) = tmp_bxmax(cur_cell, :)

               if (j > 1) then
                  tmp_bxmin(child, dir) = ksec_wall(cur_cell, j - 1)
               end if
               if (j < k) then
                  tmp_bxmax(child, dir) = ksec_wall(cur_cell, j)
               end if

               ! CPU index ranges
               if (j <= rem_ncpu) then
                  ncpu_part(j) = base_ncpu + 1
               else
                  ncpu_part(j) = base_ncpu
               end if

               if (j == 1) then
                  tmp_imin(child) = tmp_imin(cur_cell)
               else
                  tmp_imin(child) = tmp_imin(cur_cell) + cum_ncpu
               end if
               cum_ncpu = cum_ncpu + ncpu_part(j)
               tmp_imax(child) = tmp_imin(cur_cell) + cum_ncpu - 1

               ! Load estimate from histogram (if update)
               if (update) then
                  if (j == 1) then
                     tmp_load(child) = cum_load_global((i-1)*(k-1) + 1)
                  else if (j < k) then
                     tmp_load(child) = cum_load_global((i-1)*(k-1) + j) - &
                                       cum_load_global((i-1)*(k-1) + j - 1)
                  else
                     tmp_load(child) = tmp_load(cur_cell) - &
                                       cum_load_global((i-1)*(k-1) + k - 1)
                  end if
               end if
            end do
         end do

         ! Splitsort for next histogram computation
         if (update) call splitsort_ksection_histogram(lvl, dir, walls_2d, k)

         deallocate(ncpu_part)
         deallocate(walls_2d)

         ! Advance level start cursor
         cur_levelstart = cur_levelstart + nc

         ! Next level: nc grows by factor k
         nc = nc * k

      end do level_loop

      ! Process final level leaves
      do i = 1, nc
         cur_cell = cur_levelstart + (i - 1)
         if (tmp_imax(cur_cell) == 0) cycle
         if (tmp_imin(cur_cell) == tmp_imax(cur_cell)) then
            cpuid = tmp_imin(cur_cell)
            ksec_indx(cur_cell) = cpuid
            if (update) then
               bisec_cpubox_min2(cpuid, :) = tmp_bxmin(cur_cell, :)
               bisec_cpubox_max2(cpuid, :) = tmp_bxmax(cur_cell, :)
            else
               bisec_cpubox_min(cpuid, :) = tmp_bxmin(cur_cell, :)
               bisec_cpubox_max(cpuid, :) = tmp_bxmax(cur_cell, :)
            end if
            bisec_cpu_load(cpuid) = tmp_load(cur_cell)
            ksec_next(cur_cell, :) = 0
         end if
      end do

      ! Statistics
      mean  = sum(dble(bisec_cpu_load)) / ncpu
      var   = sum(dble(bisec_cpu_load) * dble(bisec_cpu_load)) / ncpu
      stdev = sqrt(max(var - mean * mean, 0.0d0))

      if (verbose .and. update) then
         print *, "K-section load balancing report ..."
         print *, "   Average CPU load     :", mean
         print *, "   Standard deviation   :", stdev
         print *, "   Balancing accuracy   :", stdev / mean
         print *, "   Requested tolerance  :", bisec_tol
      end if

      ! Store CPU ranges for hierarchical exchange
      ksec_cpumin(1:ksec_nbinodes) = tmp_imin(1:ksec_nbinodes)
      ksec_cpumax(1:ksec_nbinodes) = tmp_imax(1:ksec_nbinodes)
      call compute_ksec_cpu_path()

      ! Cleanup
      deallocate(tmp_imin, tmp_imax, tmp_load, tmp_bxmin, tmp_bxmax)
      deallocate(skip_wall, l_limit, u_limit, last_wall_arr)
      deallocate(best_wall_arr, best_score_arr, target_cum_frac)
      deallocate(myload, cum_load_global, totload_arr, skip_node)

      if(verbose) print *,'done with build_ksection'

   end subroutine build_ksection

   !================================================================
   ! Rebuild ksec_cpumin/cpumax from the tree structure (for restart)
   ! Walks the tree recursively and sets min/max CPU from ksec_indx
   !================================================================
   subroutine rebuild_ksec_cpuranges()
      implicit none

      ksec_cpumin = 0
      ksec_cpumax = 0
      call rebuild_cpuranges_node(ksec_root)

   contains

      recursive subroutine rebuild_cpuranges_node(nd)
         integer, intent(in) :: nd
         integer :: j

         if(ksec_next(nd, 1) == 0) then
            ! Leaf node
            ksec_cpumin(nd) = ksec_indx(nd)
            ksec_cpumax(nd) = ksec_indx(nd)
         else
            ! Internal node: recurse into children
            ksec_cpumin(nd) = ncpu + 1
            ksec_cpumax(nd) = 0
            do j = 1, ksec_kmax
               if(ksec_next(nd, j) == 0) exit
               call rebuild_cpuranges_node(ksec_next(nd, j))
               ksec_cpumin(nd) = min(ksec_cpumin(nd), ksec_cpumin(ksec_next(nd, j)))
               ksec_cpumax(nd) = max(ksec_cpumax(nd), ksec_cpumax(ksec_next(nd, j)))
            end do
         end if
      end subroutine rebuild_cpuranges_node

   end subroutine rebuild_ksec_cpuranges


   !================================================================
   ! Compute ksec_cpu_path: for each CPU, record its child index
   ! at each level of the k-section tree
   !================================================================
   subroutine compute_ksec_cpu_path()
      implicit none
      integer :: icpu, lvl, node, j, k

      if(allocated(ksec_cpu_path)) deallocate(ksec_cpu_path)
      allocate(ksec_cpu_path(1:ncpu, 1:nksec_levels))
      ksec_cpu_path = 0

      do icpu = 1, ncpu
         node = ksec_root
         do lvl = 1, nksec_levels
            k = ksec_factor(lvl)
            ! Find which child contains this CPU
            do j = 1, k
               if(ksec_next(node, j) == 0) exit
               if(icpu >= ksec_cpumin(ksec_next(node, j)) .and. &
                  icpu <= ksec_cpumax(ksec_next(node, j))) then
                  ksec_cpu_path(icpu, lvl) = j
                  node = ksec_next(node, j)
                  exit
               end if
            end do
         end do
      end do

   end subroutine compute_ksec_cpu_path


   !================================================================
   ! Hierarchical k-section exchange for double precision data
   ! Level-by-level correspondent exchange using k-section tree
   !================================================================
   subroutine ksection_exchange_dp(sendbuf, nitem, dest_cpu, nprops, recvbuf, nrecv)
#ifndef WITHOUTMPI
      include 'mpif.h'
#endif
      integer, intent(in) :: nitem, nprops
      real(dp), intent(in) :: sendbuf(1:nprops, 1:nitem)
      integer, intent(in) :: dest_cpu(1:nitem)
      real(dp), allocatable, intent(out) :: recvbuf(:,:)
      integer, intent(out) :: nrecv

      ! Working buffers (dynamic, managed with move_alloc)
      real(dp), allocatable :: wbuf(:,:)
      integer, allocatable :: wdest(:)
      real(dp), allocatable :: wbuf_new(:,:)
      integer, allocatable :: wdest_new(:)
      integer :: wnitem

      ! Pre-allocated per-level arrays (persistent across calls)
      integer, allocatable, save :: child_count(:), child_offset(:)
      integer, allocatable, save :: peer_list(:), peer_child(:)
      integer, allocatable, save :: send_count(:), recv_count(:)
      integer, allocatable, save :: peer_recv_offset(:)
#ifndef WITHOUTMPI
      integer, allocatable, save :: req_send_cnt(:), req_recv_cnt(:)
      integer, allocatable, save :: req_send_data(:), req_recv_data(:)
      integer, allocatable, save :: req_send_dest(:), req_recv_dest(:)
      integer, allocatable, save :: mpi_stat(:,:)
#endif
      integer, save :: kmax_alloc = 0

      ! Pre-allocated peer receive buffers (grow-only across calls)
      real(dp), allocatable, save :: peer_recv_data(:,:)
      integer, allocatable, save :: peer_recv_dest(:)
      integer, save :: prv_cap = 0, prv_np = 0

      integer :: lvl, k, node, my_child, j, c, p
      integer :: my_pos, child_size, peer, npeers
      integer :: isend_offset, irecv_total, irecv_offset
      integer :: i, idx, ierr

      ! Ensure per-level arrays are allocated to ksec_kmax
      if(ksec_kmax > kmax_alloc) then
         kmax_alloc = ksec_kmax
         if(allocated(child_count)) deallocate(child_count, child_offset)
         allocate(child_count(kmax_alloc), child_offset(kmax_alloc))
         if(allocated(peer_list)) deallocate(peer_list, peer_child, &
              send_count, recv_count, peer_recv_offset)
         allocate(peer_list(kmax_alloc-1), peer_child(kmax_alloc-1))
         allocate(send_count(kmax_alloc-1), recv_count(kmax_alloc-1))
         allocate(peer_recv_offset(kmax_alloc-1))
#ifndef WITHOUTMPI
         if(allocated(req_send_cnt)) deallocate(req_send_cnt, req_recv_cnt)
         if(allocated(req_send_data)) deallocate(req_send_data, req_recv_data, &
              req_send_dest, req_recv_dest)
         if(allocated(mpi_stat)) deallocate(mpi_stat)
         allocate(req_send_cnt(kmax_alloc-1), req_recv_cnt(kmax_alloc-1))
         allocate(req_send_data(kmax_alloc-1), req_recv_data(kmax_alloc-1))
         allocate(req_send_dest(kmax_alloc-1), req_recv_dest(kmax_alloc-1))
         allocate(mpi_stat(MPI_STATUS_SIZE, 4*(kmax_alloc-1)))
#endif
      end if

      ! Copy input to working buffers
      wnitem = nitem
      allocate(wbuf(1:nprops, 1:max(wnitem,1)))
      allocate(wdest(1:max(wnitem,1)))
      if(wnitem > 0) then
         wbuf(1:nprops, 1:wnitem) = sendbuf(1:nprops, 1:wnitem)
         wdest(1:wnitem) = dest_cpu(1:wnitem)
      end if

      ! Walk down the tree level by level
      node = ksec_root
      do lvl = 1, nksec_levels
         k = ksec_factor(lvl)
         my_child = ksec_cpu_path(myid, lvl)

         ! --- 1. Classify items by child index ---
         child_count(1:k) = 0

         ! Count items per child
         do i = 1, wnitem
            c = ksec_cpu_path(wdest(i), lvl)
            child_count(c) = child_count(c) + 1
         end do

         ! Compute offsets (exclusive prefix sum)
         child_offset(1) = 0
         do j = 2, k
            child_offset(j) = child_offset(j-1) + child_count(j-1)
         end do

         ! Counting sort into temporary arrays
         allocate(wbuf_new(1:nprops, 1:max(wnitem,1)))
         allocate(wdest_new(1:max(wnitem,1)))
         ! Reset child_count for use as running index
         child_count(1:k) = 0
         do i = 1, wnitem
            c = ksec_cpu_path(wdest(i), lvl)
            child_count(c) = child_count(c) + 1
            idx = child_offset(c) + child_count(c)
            wbuf_new(1:nprops, idx) = wbuf(1:nprops, i)
            wdest_new(idx) = wdest(i)
         end do

         ! Swap working buffers (zero-copy)
         call move_alloc(wbuf_new, wbuf)
         call move_alloc(wdest_new, wdest)

         ! --- 2. Determine correspondents ---
         npeers = k - 1

         ! My position within my child group (0-based)
         my_pos = myid - ksec_cpumin(ksec_next(node, my_child))

         p = 0
         do j = 1, k
            if(j == my_child) cycle
            p = p + 1
            child_size = ksec_cpumax(ksec_next(node, j)) - &
                         ksec_cpumin(ksec_next(node, j)) + 1
            peer_list(p) = ksec_cpumin(ksec_next(node, j)) + &
                           min(my_pos, child_size - 1)
            peer_child(p) = j
            send_count(p) = child_count(j)
         end do

#ifndef WITHOUTMPI
         ! --- 3. Exchange counts ---
         do p = 1, npeers
            call MPI_ISEND(send_count(p), 1, MPI_INTEGER, &
                 peer_list(p)-1, 100+lvl, MPI_COMM_WORLD, req_send_cnt(p), ierr)
            call MPI_IRECV(recv_count(p), 1, MPI_INTEGER, &
                 peer_list(p)-1, 100+lvl, MPI_COMM_WORLD, req_recv_cnt(p), ierr)
         end do
         if(npeers > 0) then
            call MPI_WAITALL(npeers, req_send_cnt, mpi_stat, ierr)
            call MPI_WAITALL(npeers, req_recv_cnt, mpi_stat, ierr)
         end if

         ! --- 4. Exchange data ---
         irecv_total = sum(recv_count(1:npeers))
         peer_recv_offset(1) = 0
         do p = 2, npeers
            peer_recv_offset(p) = peer_recv_offset(p-1) + recv_count(p-1)
         end do

         ! Ensure peer receive buffers are large enough (grow-only on capacity)
         ! First dimension must match nprops exactly (MPI stride requirement)
         if(irecv_total > 0) then
            if(irecv_total > prv_cap .or. nprops /= prv_np) then
               prv_cap = max(irecv_total, prv_cap)
               prv_np = nprops
               if(allocated(peer_recv_data)) deallocate(peer_recv_data, peer_recv_dest)
               allocate(peer_recv_data(1:nprops, 1:prv_cap))
               allocate(peer_recv_dest(1:prv_cap))
            end if
         end if

         do p = 1, npeers
            ! Send data for child peer_child(p)
            j = peer_child(p)

            if(send_count(p) > 0) then
               call MPI_ISEND(wbuf(1, child_offset(j)+1), send_count(p)*nprops, &
                    MPI_DOUBLE_PRECISION, peer_list(p)-1, 200+lvl, &
                    MPI_COMM_WORLD, req_send_data(p), ierr)
               call MPI_ISEND(wdest(child_offset(j)+1), send_count(p), &
                    MPI_INTEGER, peer_list(p)-1, 300+lvl, &
                    MPI_COMM_WORLD, req_send_dest(p), ierr)
            else
               req_send_data(p) = MPI_REQUEST_NULL
               req_send_dest(p) = MPI_REQUEST_NULL
            end if

            ! Recv data
            if(recv_count(p) > 0) then
               call MPI_IRECV(peer_recv_data(1, peer_recv_offset(p)+1), &
                    recv_count(p)*nprops, MPI_DOUBLE_PRECISION, &
                    peer_list(p)-1, 200+lvl, MPI_COMM_WORLD, req_recv_data(p), ierr)
               call MPI_IRECV(peer_recv_dest(peer_recv_offset(p)+1), &
                    recv_count(p), MPI_INTEGER, peer_list(p)-1, 300+lvl, &
                    MPI_COMM_WORLD, req_recv_dest(p), ierr)
            else
               req_recv_data(p) = MPI_REQUEST_NULL
               req_recv_dest(p) = MPI_REQUEST_NULL
            end if
         end do

         if(npeers > 0) then
            call MPI_WAITALL(npeers, req_send_data, mpi_stat, ierr)
            call MPI_WAITALL(npeers, req_send_dest, mpi_stat, ierr)
            call MPI_WAITALL(npeers, req_recv_data, mpi_stat, ierr)
            call MPI_WAITALL(npeers, req_recv_dest, mpi_stat, ierr)
         end if
#else
         irecv_total = 0
#endif

         ! --- 5. Merge: my_child items + received items ---
         ! My child items are at offset child_offset(my_child)+1 .. +child_count(my_child)
         idx = child_count(my_child) + irecv_total
         allocate(wbuf_new(1:nprops, 1:max(idx,1)))
         allocate(wdest_new(1:max(idx,1)))

         ! Copy my_child items
         do i = 1, child_count(my_child)
            wbuf_new(1:nprops, i) = wbuf(1:nprops, child_offset(my_child)+i)
            wdest_new(i) = wdest(child_offset(my_child)+i)
         end do

         ! Copy received items
#ifndef WITHOUTMPI
         if(irecv_total > 0) then
            do i = 1, irecv_total
               wbuf_new(1:nprops, child_count(my_child)+i) = peer_recv_data(1:nprops, i)
               wdest_new(child_count(my_child)+i) = peer_recv_dest(i)
            end do
         end if
#endif

         ! Update working set (zero-copy)
         wnitem = idx
         call move_alloc(wbuf_new, wbuf)
         call move_alloc(wdest_new, wdest)

         ! Advance to my child node for next level
         node = ksec_next(node, my_child)
      end do

      ! Output (zero-copy)
      nrecv = wnitem
      deallocate(wdest)
      if(nrecv > 0) then
         call move_alloc(wbuf, recvbuf)
      else
         deallocate(wbuf)
         allocate(recvbuf(1:nprops, 1:1))
      end if

   end subroutine ksection_exchange_dp


   !================================================================
   ! Hierarchical k-section exchange for overlapping regions
   ! Items are routed based on spatial bounding box overlap
   ! with k-section tree children (an item can go to multiple CPUs)
   !================================================================
   subroutine ksection_exchange_dp_overlap(sendbuf, nitem, item_xmin, item_xmax, &
        nprops, recvbuf, nrecv, periodic)
#ifndef WITHOUTMPI
      include 'mpif.h'
#endif
      integer, intent(in) :: nitem, nprops
      real(dp), intent(in) :: sendbuf(1:nprops, 1:nitem)
      real(dp), intent(in) :: item_xmin(1:ndim, 1:nitem)
      real(dp), intent(in) :: item_xmax(1:ndim, 1:nitem)
      real(dp), allocatable, intent(out) :: recvbuf(:,:)
      integer, intent(out) :: nrecv
      logical, intent(in), optional :: periodic

      ! Packed working buffer: data(1:nprops) + xmin(1:ndim) + xmax(1:ndim)
      integer :: npack
      real(dp), allocatable :: wpacked(:,:), wpacked_new(:,:)
      integer :: wnitem

      ! Per-child / per-peer arrays (dynamic per level)
      integer, allocatable :: my_child_idx(:)
      integer :: my_count
      real(dp), allocatable :: send_buf(:,:)
      integer :: total_send, total_recv

      ! Pre-allocated per-level arrays (persistent across calls)
      integer, allocatable, save :: peer_list(:), peer_child_idx(:)
      integer, allocatable, save :: send_count(:), recv_count(:)
      integer, allocatable, save :: send_offset(:), recv_offset(:)
#ifndef WITHOUTMPI
      integer, allocatable, save :: req_sc(:), req_rc(:)
      integer, allocatable, save :: req_sd(:), req_rd(:)
      integer, allocatable, save :: mpi_stat(:,:)
#endif
      integer, save :: kmax_alloc = 0

      ! Pre-allocated receive buffer (grow-only across calls)
      real(dp), allocatable, save :: recv_buf(:,:)
      integer, save :: rv_cap = 0, rv_np = 0

      integer :: lvl, k, dir, node, my_child, j, p
      integer :: my_pos, child_size, npeers
      integer :: i, idx, ierr
      real(dp) :: cur_bxmin(1:3), cur_bxmax(1:3)
      real(dp) :: clo, chi, xlo, xhi
      integer :: nx_loc
      real(dp) :: scale
      logical :: do_periodic
      integer :: nextra, nwrap_dims, nsubsets, mask, bit, d
      integer :: wrap_d(1:3)
      real(dp) :: wrap_shift(1:3)

      ! Ensure per-level arrays are allocated to ksec_kmax
      if(ksec_kmax > kmax_alloc) then
         kmax_alloc = ksec_kmax
         if(allocated(peer_list)) deallocate(peer_list, peer_child_idx, &
              send_count, recv_count, send_offset, recv_offset)
         allocate(peer_list(kmax_alloc-1), peer_child_idx(kmax_alloc-1))
         allocate(send_count(kmax_alloc-1), recv_count(kmax_alloc-1))
         allocate(send_offset(kmax_alloc-1), recv_offset(kmax_alloc-1))
#ifndef WITHOUTMPI
         if(allocated(req_sc)) deallocate(req_sc, req_rc, req_sd, req_rd)
         if(allocated(mpi_stat)) deallocate(mpi_stat)
         allocate(req_sc(kmax_alloc-1), req_rc(kmax_alloc-1))
         allocate(req_sd(kmax_alloc-1), req_rd(kmax_alloc-1))
         allocate(mpi_stat(MPI_STATUS_SIZE, 2*(kmax_alloc-1)))
#endif
      end if

      npack = nprops + 2 * ndim
      nx_loc = icoarse_max - icoarse_min + 1
      scale = boxlen / dble(nx_loc)

      ! Pack input into working buffer
      wnitem = nitem
      allocate(wpacked(1:npack, 1:max(wnitem,1)))
      if(wnitem > 0) then
         wpacked(1:nprops, 1:wnitem) = sendbuf(1:nprops, 1:wnitem)
         do i = 1, wnitem
            wpacked(nprops+1:nprops+ndim, i) = item_xmin(1:ndim, i)
            wpacked(nprops+ndim+1:nprops+2*ndim, i) = item_xmax(1:ndim, i)
         end do
      end if

      ! Periodic BC: create shifted copies for items that wrap around domain
      do_periodic = .false.
      if(present(periodic)) do_periodic = periodic

      if(do_periodic .and. wnitem > 0) then
         ! First pass: count extra items needed
         nextra = 0
         do i = 1, wnitem
            nsubsets = 1
            do d = 1, ndim
               xlo = wpacked(nprops + d, i)
               xhi = wpacked(nprops + ndim + d, i)
               if(xlo < 0.0d0 .or. xhi > scale) nsubsets = nsubsets * 2
            end do
            nextra = nextra + nsubsets - 1
         end do

         if(nextra > 0) then
            allocate(wpacked_new(1:npack, 1:wnitem + nextra))
            wpacked_new(1:npack, 1:wnitem) = wpacked(1:npack, 1:wnitem)
            idx = wnitem

            do i = 1, wnitem
               nwrap_dims = 0
               do d = 1, ndim
                  xlo = wpacked(nprops + d, i)
                  xhi = wpacked(nprops + ndim + d, i)
                  if(xlo < 0.0d0) then
                     nwrap_dims = nwrap_dims + 1
                     wrap_d(nwrap_dims) = d
                     wrap_shift(nwrap_dims) = scale
                  else if(xhi > scale) then
                     nwrap_dims = nwrap_dims + 1
                     wrap_d(nwrap_dims) = d
                     wrap_shift(nwrap_dims) = -scale
                  end if
               end do

               if(nwrap_dims == 0) cycle

               ! Generate all non-empty subsets of wrapping dimensions
               nsubsets = 2**nwrap_dims - 1
               do mask = 1, nsubsets
                  idx = idx + 1
                  wpacked_new(1:npack, idx) = wpacked(1:npack, i)
                  do bit = 1, nwrap_dims
                     if(iand(mask, 2**(bit-1)) /= 0) then
                        wpacked_new(nprops + wrap_d(bit), idx) = &
                             wpacked_new(nprops + wrap_d(bit), idx) + wrap_shift(bit)
                        wpacked_new(nprops + ndim + wrap_d(bit), idx) = &
                             wpacked_new(nprops + ndim + wrap_d(bit), idx) + wrap_shift(bit)
                     end if
                  end do
               end do
            end do

            wnitem = wnitem + nextra
            call move_alloc(wpacked_new, wpacked)
         end if
      end if

      ! Initialize node bounding box
      cur_bxmin(1:ndim) = 0.0d0
      cur_bxmax(1:ndim) = scale

      ! Walk down the tree level by level
      node = ksec_root
      do lvl = 1, nksec_levels
         k = ksec_factor(lvl)
         dir = ksec_dir(lvl)
         my_child = ksec_cpu_path(myid, lvl)

         ! --- Determine correspondents ---
         npeers = k - 1

         my_pos = myid - ksec_cpumin(ksec_next(node, my_child))
         p = 0
         do j = 1, k
            if(j == my_child) cycle
            p = p + 1
            child_size = ksec_cpumax(ksec_next(node,j)) - &
                         ksec_cpumin(ksec_next(node,j)) + 1
            peer_list(p) = ksec_cpumin(ksec_next(node,j)) + &
                           min(my_pos, child_size - 1)
            peer_child_idx(p) = j
         end do

         ! --- First pass: count overlaps ---
         my_count = 0
         send_count(1:npeers) = 0
         do i = 1, wnitem
            xlo = wpacked(nprops + dir, i)
            xhi = wpacked(nprops + ndim + dir, i)
            do j = 1, k
               ! Child j boundary along split direction
               if(j == 1) then; clo = cur_bxmin(dir)
               else;            clo = ksec_wall(node, j-1); end if
               if(j == k) then; chi = cur_bxmax(dir)
               else;            chi = ksec_wall(node, j);   end if
               ! Check overlap (strict inequalities)
               if(xhi > clo .and. xlo < chi) then
                  if(j == my_child) then
                     my_count = my_count + 1
                  else
                     do p = 1, npeers
                        if(peer_child_idx(p) == j) then
                           send_count(p) = send_count(p) + 1
                           exit
                        end if
                     end do
                  end if
               end if
            end do
         end do

         ! --- Allocate buffers ---
         total_send = sum(send_count(1:npeers))
         send_offset(1) = 0
         do p = 2, npeers
            send_offset(p) = send_offset(p-1) + send_count(p-1)
         end do
         if(total_send > 0) allocate(send_buf(1:npack, 1:total_send))
         allocate(my_child_idx(1:max(my_count,1)))

         ! --- Second pass: fill buffers ---
         my_count = 0
         send_count(1:npeers) = 0
         do i = 1, wnitem
            xlo = wpacked(nprops + dir, i)
            xhi = wpacked(nprops + ndim + dir, i)
            do j = 1, k
               if(j == 1) then; clo = cur_bxmin(dir)
               else;            clo = ksec_wall(node, j-1); end if
               if(j == k) then; chi = cur_bxmax(dir)
               else;            chi = ksec_wall(node, j);   end if
               if(xhi > clo .and. xlo < chi) then
                  if(j == my_child) then
                     my_count = my_count + 1
                     my_child_idx(my_count) = i
                  else
                     do p = 1, npeers
                        if(peer_child_idx(p) == j) then
                           send_count(p) = send_count(p) + 1
                           idx = send_offset(p) + send_count(p)
                           send_buf(1:npack, idx) = wpacked(1:npack, i)
                           exit
                        end if
                     end do
                  end if
               end if
            end do
         end do

#ifndef WITHOUTMPI
         ! --- Exchange counts ---
         do p = 1, npeers
            call MPI_ISEND(send_count(p), 1, MPI_INTEGER, &
                 peer_list(p)-1, 400+lvl, MPI_COMM_WORLD, req_sc(p), ierr)
            call MPI_IRECV(recv_count(p), 1, MPI_INTEGER, &
                 peer_list(p)-1, 400+lvl, MPI_COMM_WORLD, req_rc(p), ierr)
         end do
         if(npeers > 0) then
            call MPI_WAITALL(npeers, req_sc, mpi_stat, ierr)
            call MPI_WAITALL(npeers, req_rc, mpi_stat, ierr)
         end if

         ! --- Exchange data ---
         total_recv = sum(recv_count(1:npeers))
         recv_offset(1) = 0
         do p = 2, npeers
            recv_offset(p) = recv_offset(p-1) + recv_count(p-1)
         end do

         ! Ensure receive buffer is large enough (grow-only on capacity)
         ! First dimension must match npack exactly (MPI stride requirement)
         if(total_recv > 0) then
            if(total_recv > rv_cap .or. npack /= rv_np) then
               rv_cap = max(total_recv, rv_cap)
               rv_np = npack
               if(allocated(recv_buf)) deallocate(recv_buf)
               allocate(recv_buf(1:npack, 1:rv_cap))
            end if
         end if

         do p = 1, npeers
            if(send_count(p) > 0) then
               call MPI_ISEND(send_buf(1, send_offset(p)+1), &
                    send_count(p)*npack, MPI_DOUBLE_PRECISION, &
                    peer_list(p)-1, 500+lvl, MPI_COMM_WORLD, req_sd(p), ierr)
            else
               req_sd(p) = MPI_REQUEST_NULL
            end if
            if(recv_count(p) > 0) then
               call MPI_IRECV(recv_buf(1, recv_offset(p)+1), &
                    recv_count(p)*npack, MPI_DOUBLE_PRECISION, &
                    peer_list(p)-1, 500+lvl, MPI_COMM_WORLD, req_rd(p), ierr)
            else
               req_rd(p) = MPI_REQUEST_NULL
            end if
         end do
         if(npeers > 0) then
            call MPI_WAITALL(npeers, req_sd, mpi_stat, ierr)
            call MPI_WAITALL(npeers, req_rd, mpi_stat, ierr)
         end if
#else
         total_recv = 0
#endif

         ! --- Merge: my_child items + received items ---
         idx = my_count + total_recv
         allocate(wpacked_new(1:npack, 1:max(idx,1)))
         do i = 1, my_count
            wpacked_new(1:npack, i) = wpacked(1:npack, my_child_idx(i))
         end do
#ifndef WITHOUTMPI
         if(total_recv > 0) then
            do i = 1, total_recv
               wpacked_new(1:npack, my_count+i) = recv_buf(1:npack, i)
            end do
         end if
#endif

         ! Update working set (zero-copy)
         wnitem = idx
         call move_alloc(wpacked_new, wpacked)

         ! Update node bounding box for my_child
         if(my_child > 1) cur_bxmin(dir) = ksec_wall(node, my_child - 1)
         if(my_child < k) cur_bxmax(dir) = ksec_wall(node, my_child)

         ! Cleanup (only dynamic per-level arrays)
         deallocate(my_child_idx)
         if(allocated(send_buf)) deallocate(send_buf)

         ! Advance to my child node
         node = ksec_next(node, my_child)
      end do

      ! Output: extract data portion only
      nrecv = wnitem
      if(nrecv > 0) then
         allocate(recvbuf(1:nprops, 1:nrecv))
         recvbuf(1:nprops, 1:nrecv) = wpacked(1:nprops, 1:nrecv)
      else
         allocate(recvbuf(1:nprops, 1:1))
      end if
      deallocate(wpacked)

   end subroutine ksection_exchange_dp_overlap


   !================================================================
   ! Test subroutine for ksection exchange (exclusive + overlap)
   !================================================================
   subroutine test_ksection_exchange()
#ifndef WITHOUTMPI
      include 'mpif.h'
#endif
      integer, parameter :: np = 3  ! sender, target, index
      integer :: nitem, nrecv, ierr
      integer :: i, d, pass_local, pass_global
      real(dp), allocatable :: sbuf(:,:), rbuf(:,:)
      integer, allocatable :: dcpu(:)
      real(dp), allocatable :: xmn(:,:), xmx(:,:)
      real(dp) :: ctr
      integer :: nx_loc
      real(dp) :: scale
      integer :: nexpected
      logical :: overlap_all, dim_ok
      real(dp) :: ictr, ixmin_d, ixmax_d, mylo, myhi

      nx_loc = icoarse_max - icoarse_min + 1
      scale = boxlen / dble(nx_loc)

      ! ===== Test 1: Exclusive exchange =====
      ! Each rank sends 1 item to every rank (including self)
      nitem = ncpu
      allocate(sbuf(1:np, 1:nitem), dcpu(1:nitem))
      do i = 1, ncpu
         dcpu(i) = i
         sbuf(1, i) = dble(myid)  ! sender
         sbuf(2, i) = dble(i)     ! target
         sbuf(3, i) = dble(i)     ! index
      end do

      call ksection_exchange_dp(sbuf, nitem, dcpu, np, rbuf, nrecv)

      pass_local = 1
      if(nrecv /= ncpu) then
         pass_local = 0
         write(*,*) 'EXCL FAIL: rank', myid, 'nrecv=', nrecv, 'expected', ncpu
      end if
      do i = 1, nrecv
         if(nint(rbuf(2, i)) /= myid) then
            pass_local = 0
            write(*,*) 'EXCL FAIL: rank', myid, 'got dest=', nint(rbuf(2,i))
         end if
      end do

#ifndef WITHOUTMPI
      call MPI_ALLREDUCE(pass_local, pass_global, 1, MPI_INTEGER, &
           MPI_MIN, MPI_COMM_WORLD, ierr)
#else
      pass_global = pass_local
#endif
      if(myid==1) then
         if(pass_global==1) then
            write(*,*) 'ksection_exchange (exclusive) test: PASSED'
         else
            write(*,*) 'ksection_exchange (exclusive) test: FAILED'
         end if
      end if
      deallocate(sbuf, dcpu, rbuf)

      ! ===== Test 2: Overlap exchange (point boxes at domain centers) =====
      ! Each rank sends 1 item per CPU, with a tiny box at that CPU's domain center
      ! Result: each rank should receive ncpu items targeted at it
      nitem = ncpu
      allocate(sbuf(1:np, 1:nitem))
      allocate(xmn(1:ndim, 1:nitem), xmx(1:ndim, 1:nitem))
      do i = 1, ncpu
         sbuf(1, i) = dble(myid)
         sbuf(2, i) = dble(i)
         sbuf(3, i) = dble(i)
         do d = 1, ndim
            ctr = 0.5d0 * (bisec_cpubox_min(i, d) + bisec_cpubox_max(i, d))
            xmn(d, i) = ctr - 1.0d-10
            xmx(d, i) = ctr + 1.0d-10
         end do
      end do

      call ksection_exchange_dp_overlap(sbuf, nitem, xmn, xmx, np, rbuf, nrecv)

      pass_local = 1
      if(nrecv /= ncpu) then
         pass_local = 0
         write(*,*) 'OVLP FAIL: rank', myid, 'nrecv=', nrecv, 'expected', ncpu
      end if
      do i = 1, nrecv
         if(nint(rbuf(2, i)) /= myid) then
            pass_local = 0
            write(*,*) 'OVLP FAIL: rank', myid, 'got target=', nint(rbuf(2,i))
         end if
      end do

#ifndef WITHOUTMPI
      call MPI_ALLREDUCE(pass_local, pass_global, 1, MPI_INTEGER, &
           MPI_MIN, MPI_COMM_WORLD, ierr)
#else
      pass_global = pass_local
#endif
      if(myid==1) then
         if(pass_global==1) then
            write(*,*) 'ksection_exchange (overlap point) test: PASSED'
         else
            write(*,*) 'ksection_exchange (overlap point) test: FAILED'
         end if
      end if
      deallocate(sbuf, xmn, xmx, rbuf)

      ! ===== Test 3: Overlap exchange (full-domain box) =====
      ! Each rank sends 1 item covering entire domain
      ! Result: every rank should receive ncpu items (one from each sender)
      nitem = 1
      allocate(sbuf(1:np, 1:1))
      allocate(xmn(1:ndim, 1:1), xmx(1:ndim, 1:1))
      sbuf(1, 1) = dble(myid)
      sbuf(2, 1) = 0.0d0
      sbuf(3, 1) = 1.0d0
      do d = 1, ndim
         xmn(d, 1) = 0.0d0
         xmx(d, 1) = scale
      end do

      call ksection_exchange_dp_overlap(sbuf, nitem, xmn, xmx, np, rbuf, nrecv)

      pass_local = 1
      if(nrecv /= ncpu) then
         pass_local = 0
         write(*,*) 'FULL FAIL: rank', myid, 'nrecv=', nrecv, 'expected', ncpu
      end if

#ifndef WITHOUTMPI
      call MPI_ALLREDUCE(pass_local, pass_global, 1, MPI_INTEGER, &
           MPI_MIN, MPI_COMM_WORLD, ierr)
#else
      pass_global = pass_local
#endif
      if(myid==1) then
         if(pass_global==1) then
            write(*,*) 'ksection_exchange (full overlap) test: PASSED'
         else
            write(*,*) 'ksection_exchange (full overlap) test: FAILED'
         end if
      end if
      deallocate(sbuf, xmn, xmx, rbuf)

      ! ===== Test 4: Periodic overlap exchange =====
      ! Each rank sends 1 item at domain center with radius 20% of scale
      nitem = 1
      allocate(sbuf(1:np, 1:1))
      allocate(xmn(1:ndim, 1:1), xmx(1:ndim, 1:1))
      sbuf(1, 1) = dble(myid)
      sbuf(2, 1) = 0.0d0
      sbuf(3, 1) = 1.0d0
      do d = 1, ndim
         ctr = 0.5d0 * (bisec_cpubox_min(myid, d) + bisec_cpubox_max(myid, d))
         xmn(d, 1) = ctr - 0.2d0 * scale
         xmx(d, 1) = ctr + 0.2d0 * scale
      end do

      call ksection_exchange_dp_overlap(sbuf, nitem, xmn, xmx, np, rbuf, nrecv, periodic=.true.)

      ! Verify: compute expected receives by brute-force periodic overlap check
      pass_local = 1
      nexpected = 0
      do i = 1, ncpu
         overlap_all = .true.
         do d = 1, ndim
            ictr = 0.5d0 * (bisec_cpubox_min(i, d) + bisec_cpubox_max(i, d))
            ixmin_d = ictr - 0.2d0 * scale
            ixmax_d = ictr + 0.2d0 * scale
            mylo = bisec_cpubox_min(myid, d)
            myhi = bisec_cpubox_max(myid, d)
            ! Check overlap with any periodic image
            dim_ok = .false.
            if(ixmax_d > mylo .and. ixmin_d < myhi) dim_ok = .true.
            if(ixmax_d + scale > mylo .and. ixmin_d + scale < myhi) dim_ok = .true.
            if(ixmax_d - scale > mylo .and. ixmin_d - scale < myhi) dim_ok = .true.
            if(.not. dim_ok) overlap_all = .false.
         end do
         if(overlap_all) nexpected = nexpected + 1
      end do

      if(nrecv /= nexpected) then
         pass_local = 0
         write(*,*) 'PERI FAIL: rank', myid, 'nrecv=', nrecv, 'expected=', nexpected
      end if

#ifndef WITHOUTMPI
      call MPI_ALLREDUCE(pass_local, pass_global, 1, MPI_INTEGER, &
           MPI_MIN, MPI_COMM_WORLD, ierr)
#else
      pass_global = pass_local
#endif
      if(myid==1) then
         if(pass_global==1) then
            write(*,*) 'ksection_exchange (periodic overlap) test: PASSED'
         else
            write(*,*) 'ksection_exchange (periodic overlap) test: FAILED'
         end if
      end if
      deallocate(sbuf, xmn, xmx, rbuf)

   end subroutine test_ksection_exchange


end module ksection
