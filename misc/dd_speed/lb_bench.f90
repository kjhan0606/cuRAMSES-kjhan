!=====================================================================
! lb_bench.f90 - Load Balance Quality Benchmark
!
! Compares load balance quality of Hilbert vs K-Section DD using a
! non-uniform weight field with Gaussian clusters.
!
! Metrics: weight imbalance, ghost cell ratio, neighbor partner count
!
! Plane-by-plane ownership computation: O(N^2) memory, not O(N^3).
!
! Usage: mpirun -np P ./lb_bench [ngrid] [nclusters]
!        Default: ngrid=256, nclusters=20
!=====================================================================

!---------------------------------------------------------------------
! Hilbert curve module (from ghost_bench.f90)
!---------------------------------------------------------------------
module hilbert_mod
  implicit none
contains

  subroutine hilbert3d(x, y, z, order, bit_length, npoint)
    implicit none
    integer, intent(in) :: bit_length, npoint
    integer, intent(in) :: x(npoint), y(npoint), z(npoint)
    real(8), intent(out) :: order(npoint)
    logical :: ib(0:3*bit_length-1)
    logical :: xb(0:bit_length-1), yb(0:bit_length-1), zb(0:bit_length-1)
    integer :: sd(0:7,0:1,0:11)
    integer :: i, ip, cs, ns, b0, b1, b2, si, hi

    sd = RESHAPE((/ &
      1,2,3,2,4,5,3,5, 0,1,3,2,7,6,4,5, &
      2,6,0,7,8,8,0,7, 0,7,1,6,3,4,2,5, &
      0,9,10,9,1,1,11,11, 0,3,7,4,1,2,6,5, &
      6,0,6,11,9,0,9,8, 2,3,1,0,5,4,6,7, &
      11,11,0,7,5,9,0,7, 4,3,5,2,7,0,6,1, &
      4,4,8,8,0,6,10,6, 6,5,1,2,7,4,0,3, &
      5,7,5,3,1,1,11,11, 4,7,3,0,5,6,2,1, &
      6,1,6,10,9,4,9,10, 6,7,5,4,1,0,2,3, &
      10,3,1,1,10,3,5,9, 2,5,3,4,1,6,0,7, &
      4,4,8,8,2,7,2,3, 2,1,5,6,3,0,4,7, &
      7,2,11,2,7,5,8,5, 4,5,7,6,3,2,0,1, &
      10,3,2,6,10,3,4,4, 6,1,7,0,5,2,4,3/), (/8,2,12/))

    do ip = 1, npoint
      do i = 0, bit_length-1
        xb(i) = btest(x(ip),i); yb(i) = btest(y(ip),i); zb(i) = btest(z(ip),i)
      end do
      do i = 0, bit_length-1
        ib(3*i+2) = xb(i); ib(3*i+1) = yb(i); ib(3*i) = zb(i)
      end do
      cs = 0
      do i = bit_length-1, 0, -1
        b2=0; if(ib(3*i+2)) b2=1
        b1=0; if(ib(3*i+1)) b1=1
        b0=0; if(ib(3*i  )) b0=1
        si = b2*4+b1*2+b0
        ns = sd(si,0,cs); hi = sd(si,1,cs)
        ib(3*i+2) = btest(hi,2); ib(3*i+1) = btest(hi,1); ib(3*i) = btest(hi,0)
        cs = ns
      end do
      order(ip) = 0.0d0
      do i = 0, 3*bit_length-1
        b0=0; if(ib(i)) b0=1
        order(ip) = order(ip) + real(b0,kind=8)*real(2,kind=8)**i
      end do
    end do
  end subroutine hilbert3d

end module hilbert_mod

!---------------------------------------------------------------------
! K-Section module: tree DD (no exchange needed for this benchmark)
!---------------------------------------------------------------------
module ksection_mod
  use mpi
  implicit none

  integer, parameter :: KS_MAXLEV = 32
  integer :: nklev, kfac(KS_MAXLEV), kkmax

  integer :: knn
  integer, allocatable :: kcmin(:), kcmax(:)
  real(8), allocatable :: kbmin(:,:), kbmax(:,:)
  integer, allocatable :: knc(:), kco(:), kdir(:)

  integer :: kpnode(KS_MAXLEV+1), kpchild(KS_MAXLEV)

contains

  subroutine ks_factorize(n)
    implicit none
    integer, intent(in) :: n
    integer :: m, d, nf, i, j, tmp

    kfac = 0; nf = 0; m = n; d = 2
    do while(m > 1)
      do while(mod(m,d) == 0)
        nf = nf + 1; kfac(nf) = d; m = m / d
      end do
      d = d + 1
      if(d*d > m .and. m > 1) then
        nf = nf + 1; kfac(nf) = m; m = 1
      end if
    end do
    do i = 1, nf-1
      do j = i+1, nf
        if(kfac(j) > kfac(i)) then
          tmp = kfac(i); kfac(i) = kfac(j); kfac(j) = tmp
        end if
      end do
    end do
    nklev = nf
    kkmax = 2
    if(nf > 0) kkmax = kfac(1)
  end subroutine

  subroutine ks_build(ncpu)
    implicit none
    integer, intent(in) :: ncpu
    integer :: mxn, lvl, k, nd, cn, j, dir
    integer :: ls, le, ns, cppc, cs, ce

    call ks_factorize(ncpu)

    mxn = 2*ncpu + 16
    if(allocated(kcmin)) deallocate(kcmin,kcmax,kbmin,kbmax,knc,kco,kdir)
    allocate(kcmin(mxn), kcmax(mxn), kbmin(3,mxn), kbmax(3,mxn))
    allocate(knc(mxn), kco(mxn), kdir(mxn))

    knn = 1
    kcmin(1) = 0; kcmax(1) = ncpu - 1
    kbmin(:,1) = 0.0d0; kbmax(:,1) = 1.0d0
    knc(1) = 0; kco(1) = 0; kdir(1) = 0

    ls = 1; le = 1
    do lvl = 1, nklev
      k = kfac(lvl); ns = knn + 1
      do nd = ls, le
        dir = 1
        if(kbmax(2,nd)-kbmin(2,nd) > kbmax(dir,nd)-kbmin(dir,nd)) dir = 2
        if(kbmax(3,nd)-kbmin(3,nd) > kbmax(dir,nd)-kbmin(dir,nd)) dir = 3
        kdir(nd) = dir; kco(nd) = knn; knc(nd) = k
        cppc = (kcmax(nd) - kcmin(nd) + 1) / k; cs = kcmin(nd)
        do j = 1, k
          knn = knn + 1; cn = knn
          kbmin(:,cn) = kbmin(:,nd); kbmax(:,cn) = kbmax(:,nd)
          kbmin(dir,cn) = kbmin(dir,nd) + dble(j-1)/dble(k)*(kbmax(dir,nd)-kbmin(dir,nd))
          kbmax(dir,cn) = kbmin(dir,nd) + dble(j)  /dble(k)*(kbmax(dir,nd)-kbmin(dir,nd))
          ce = cs + cppc - 1; kcmin(cn) = cs; kcmax(cn) = ce; cs = ce + 1
          knc(cn) = 0; kco(cn) = 0; kdir(cn) = 0
        end do
      end do
      ls = ns; le = knn
    end do
  end subroutine

  subroutine ks_path(myid)
    implicit none
    integer, intent(in) :: myid
    integer :: lvl, nd, j, cn

    nd = 1; kpnode(1) = 1
    do lvl = 1, nklev
      do j = 1, knc(nd)
        cn = kco(nd) + j
        if(myid >= kcmin(cn) .and. myid <= kcmax(cn)) then
          kpchild(lvl) = j; kpnode(lvl+1) = cn; nd = cn; exit
        end if
      end do
    end do
  end subroutine

  subroutine ks_dd(pts, np, dest)
    implicit none
    integer, intent(in) :: np
    real(8), intent(in) :: pts(3, np)
    integer, intent(out) :: dest(np)
    integer :: i, nd, lvl, j, k, dir
    real(8) :: c, lo, hi

    do i = 1, np
      nd = 1
      do lvl = 1, nklev
        k = knc(nd); if(k == 0) exit
        dir = kdir(nd)
        c = pts(dir, i); lo = kbmin(dir, nd); hi = kbmax(dir, nd)
        j = int((c - lo) / (hi - lo) * dble(k)) + 1
        j = max(1, min(k, j))
        nd = kco(nd) + j
      end do
      dest(i) = kcmin(nd)
    end do
  end subroutine

end module ksection_mod

!=====================================================================
! Main benchmark program
!=====================================================================
program lb_bench
  use mpi
  use hilbert_mod
  use ksection_mod
  implicit none

  integer :: myid, ncpu, ierr, comm
  integer :: ngrid, nclusters
  character(len=32) :: arg
  integer(8) :: ncells8

  ! Cluster parameters
  real(8), allocatable :: cx(:), cy(:), cz(:), amp(:), sigma(:)

  ! Hilbert parameters
  integer :: bit_length, hgrid
  real(8) :: max_key

  ! Work arrays for plane computation (N^2 size)
  integer :: nplane
  integer, allocatable :: work_ix(:), work_iy(:), work_iz(:), work_dd(:)
  real(8), allocatable :: work_hk(:), work_pts(:,:)

  ! DD statistics (local)
  real(8) :: my_weight
  integer :: my_ghost, my_owned
  integer, allocatable :: partner_flag(:)
  integer :: my_partners

  ! DD statistics (global)
  real(8) :: wt_min, wt_max, wt_sum, wt_mean, wt_imbal
  integer :: gh_min, gh_max, gh_sum
  real(8) :: gh_mean, gh_ratio
  integer :: pt_min, pt_max, pt_sum
  real(8) :: pt_mean

  real(8) :: t0, t1, t_dd
  real(8) :: total_weight

  integer :: i, ic, idd
  character(len=10) :: dd_name(2)

  call MPI_INIT(ierr)
  comm = MPI_COMM_WORLD
  call MPI_COMM_RANK(comm, myid, ierr)
  call MPI_COMM_SIZE(comm, ncpu, ierr)

  ! Parse command line
  ngrid = 256; nclusters = 20
  if(command_argument_count() >= 1) then
    call get_command_argument(1, arg); read(arg,*) ngrid
  end if
  if(command_argument_count() >= 2) then
    call get_command_argument(2, arg); read(arg,*) nclusters
  end if

  ncells8 = int(ngrid, 8)**3
  nplane = ngrid * ngrid

  ! Generate deterministic cluster parameters (same on all ranks)
  allocate(cx(nclusters), cy(nclusters), cz(nclusters))
  allocate(amp(nclusters), sigma(nclusters))
  do ic = 1, nclusters
    cx(ic)    = mod(ic*1327 + 853, 10000) / 10000.0d0
    cy(ic)    = mod(ic*2659 + 1741, 10000) / 10000.0d0
    cz(ic)    = mod(ic*3967 + 2111, 10000) / 10000.0d0
    amp(ic)   = 10.0d0 + mod(ic*4273, 90)
    sigma(ic) = 0.02d0 + mod(ic*5381, 80) * 0.001d0
  end do

  ! Compute total weight (distributed across ranks, then reduced)
  total_weight = compute_total_weight(ngrid, nclusters, cx, cy, cz, amp, sigma, &
                                       myid, ncpu, comm)

  ! Pre-compute Hilbert parameters
  bit_length = 1
  do while(2**bit_length < ngrid)
    bit_length = bit_length + 1
  end do
  hgrid = 2**bit_length
  max_key = 2.0d0**(3*bit_length)

  ! Build K-Section tree + path
  call ks_build(ncpu)
  call ks_path(myid)

  ! Print header
  if(myid == 0) then
    write(*,'(A)') '=== Load Balance Benchmark ==='
    write(*,'(A,I0,A,I0,A,I0)') 'ncpu=', ncpu, ', N=', ngrid, ', nclusters=', nclusters
    write(*,*)

    write(*,'(A)') 'Cluster parameters:'
    write(*,'(A6, A10, A10, A10, A8, A10)') &
         'ic', 'cx', 'cy', 'cz', 'amp', 'sigma'
    do ic = 1, nclusters
      write(*,'(I6, F10.4, F10.4, F10.4, F8.1, F10.4)') &
           ic, cx(ic), cy(ic), cz(ic), amp(ic), sigma(ic)
    end do
    write(*,*)

    write(*,'(A,ES14.6)') 'Total weight = ', total_weight
    write(*,'(A,ES14.6)') 'Mean weight per CPU = ', total_weight / dble(ncpu)
    write(*,*)

    write(*,'(A)', advance='no') 'K-Section factors: '
    do i = 1, nklev
      if(i > 1) write(*,'(A)', advance='no') ' x '
      write(*,'(I0)', advance='no') kfac(i)
    end do
    write(*,'(A,I0,A)') '  (', nklev, ' levels)'
    write(*,*)
  end if

  ! Allocate work arrays
  allocate(work_ix(nplane), work_iy(nplane), work_iz(nplane))
  allocate(work_hk(nplane), work_pts(3, nplane), work_dd(nplane))
  allocate(partner_flag(0:ncpu-1))

  dd_name(1) = 'Hilbert  '; dd_name(2) = 'K-Section'

  ! ===== Benchmark each DD method =====
  do idd = 1, 2
    call MPI_BARRIER(comm, ierr)
    t0 = MPI_WTIME()

    call compute_dd_stats(idd, ngrid, nclusters, cx, cy, cz, amp, sigma, &
                          myid, ncpu, my_weight, my_ghost, my_owned, my_partners)

    call MPI_BARRIER(comm, ierr)
    t1 = MPI_WTIME()
    t_dd = t1 - t0

    ! Global reductions for weight
    wt_min = my_weight; wt_max = my_weight; wt_sum = my_weight
    call MPI_ALLREDUCE(MPI_IN_PLACE, wt_min, 1, MPI_DOUBLE_PRECISION, MPI_MIN, comm, ierr)
    call MPI_ALLREDUCE(MPI_IN_PLACE, wt_max, 1, MPI_DOUBLE_PRECISION, MPI_MAX, comm, ierr)
    call MPI_ALLREDUCE(MPI_IN_PLACE, wt_sum, 1, MPI_DOUBLE_PRECISION, MPI_SUM, comm, ierr)
    wt_mean = wt_sum / dble(ncpu)
    wt_imbal = wt_max / wt_mean

    ! Global reductions for ghost count
    gh_min = my_ghost; gh_max = my_ghost; gh_sum = my_ghost
    call MPI_ALLREDUCE(MPI_IN_PLACE, gh_min, 1, MPI_INTEGER, MPI_MIN, comm, ierr)
    call MPI_ALLREDUCE(MPI_IN_PLACE, gh_max, 1, MPI_INTEGER, MPI_MAX, comm, ierr)
    call MPI_ALLREDUCE(MPI_IN_PLACE, gh_sum, 1, MPI_INTEGER, MPI_SUM, comm, ierr)
    gh_mean = dble(gh_sum) / dble(ncpu)

    ! Ghost ratio: use global mean ghost / mean owned
    gh_ratio = gh_mean / (dble(ncells8) / dble(ncpu))

    ! Global reductions for partner count
    pt_min = my_partners; pt_max = my_partners; pt_sum = my_partners
    call MPI_ALLREDUCE(MPI_IN_PLACE, pt_min, 1, MPI_INTEGER, MPI_MIN, comm, ierr)
    call MPI_ALLREDUCE(MPI_IN_PLACE, pt_max, 1, MPI_INTEGER, MPI_MAX, comm, ierr)
    call MPI_ALLREDUCE(MPI_IN_PLACE, pt_sum, 1, MPI_INTEGER, MPI_SUM, comm, ierr)
    pt_mean = dble(pt_sum) / dble(ncpu)

    if(myid == 0) then
      write(*,'(A,A,A,F8.1,A)') '--- ', trim(dd_name(idd)), &
           ' DD (time: ', t_dd, 's) ---'
      write(*,'(A,ES12.4,A,ES12.4,A,ES12.4,A,F6.3)') &
           '  Weight:   min=', wt_min, '  max=', wt_max, &
           '  mean=', wt_mean, '  imbalance=', wt_imbal
      write(*,'(A,I0,A,I0,A,F12.1,A,F8.5)') &
           '  Ghost:    min=', gh_min, '  max=', gh_max, &
           '  mean=', gh_mean, '  ratio=', gh_ratio
      write(*,'(A,I0,A,I0,A,F6.1)') &
           '  Partners: min=', pt_min, '  max=', pt_max, &
           '  mean=', pt_mean
      write(*,*)
    end if
  end do

  ! Cleanup
  deallocate(cx, cy, cz, amp, sigma)
  deallocate(work_ix, work_iy, work_iz, work_hk, work_pts, work_dd)
  deallocate(partner_flag)

  call MPI_FINALIZE(ierr)

contains

  !-------------------------------------------------------------------
  ! Compute cell weight on-the-fly from cluster parameters
  !-------------------------------------------------------------------
  pure real(8) function cell_weight(i, j, k, N, cx_in, cy_in, cz_in, &
                                     amp_in, sigma_in, nc)
    integer, intent(in) :: i, j, k, N, nc
    real(8), intent(in) :: cx_in(nc), cy_in(nc), cz_in(nc)
    real(8), intent(in) :: amp_in(nc), sigma_in(nc)
    real(8) :: x, y, z, dx, dy, dz, r2, w
    integer :: ic2

    x = (dble(i) - 0.5d0) / dble(N)
    y = (dble(j) - 0.5d0) / dble(N)
    z = (dble(k) - 0.5d0) / dble(N)
    w = 1.0d0
    do ic2 = 1, nc
      dx = min(abs(x - cx_in(ic2)), 1.0d0 - abs(x - cx_in(ic2)))
      dy = min(abs(y - cy_in(ic2)), 1.0d0 - abs(y - cy_in(ic2)))
      dz = min(abs(z - cz_in(ic2)), 1.0d0 - abs(z - cz_in(ic2)))
      r2 = dx*dx + dy*dy + dz*dz
      w = w + amp_in(ic2) * exp(-r2 / (2.0d0 * sigma_in(ic2)**2))
    end do
    cell_weight = w
  end function

  !-------------------------------------------------------------------
  ! Compute total weight (distributed across ranks, then reduced)
  !-------------------------------------------------------------------
  real(8) function compute_total_weight(N, nc, cx_in, cy_in, cz_in, &
                                         amp_in, sigma_in, myid_in, ncpu_in, comm_in)
    integer, intent(in) :: N, nc, myid_in, ncpu_in, comm_in
    real(8), intent(in) :: cx_in(nc), cy_in(nc), cz_in(nc)
    real(8), intent(in) :: amp_in(nc), sigma_in(nc)
    real(8) :: local_sum
    integer :: i, j, k, ierr2

    local_sum = 0.0d0
    do k = 1, N
      if(mod(k-1, ncpu_in) /= myid_in) cycle
      do j = 1, N
        do i = 1, N
          local_sum = local_sum + cell_weight(i, j, k, N, cx_in, cy_in, cz_in, &
                                               amp_in, sigma_in, nc)
        end do
      end do
    end do
    call MPI_ALLREDUCE(local_sum, compute_total_weight, 1, &
         MPI_DOUBLE_PRECISION, MPI_SUM, comm_in, ierr2)
  end function

  !-------------------------------------------------------------------
  ! Compute ownership of all N^2 cells in z-plane kplane (1..N)
  ! Uses host-associated: bit_length, hgrid, max_key, ncpu, work_*
  !-------------------------------------------------------------------
  subroutine compute_plane_owner(idd_type, kplane, N, owners)
    integer, intent(in) :: idd_type, kplane, N
    integer, intent(out) :: owners(N, N)
    integer :: ii, jj, pp, np, iz_val

    np = N * N

    if(idd_type == 1) then
      ! Hilbert DD
      iz_val = min(int((dble(kplane)-0.5d0)/dble(N)*dble(hgrid)), hgrid-1)
      pp = 0
      do jj = 1, N; do ii = 1, N
        pp = pp + 1
        work_ix(pp) = min(int((dble(ii)-0.5d0)/dble(N)*dble(hgrid)), hgrid-1)
        work_iy(pp) = min(int((dble(jj)-0.5d0)/dble(N)*dble(hgrid)), hgrid-1)
        work_iz(pp) = iz_val
      end do; end do
      call hilbert3d(work_ix, work_iy, work_iz, work_hk, bit_length, np)
      pp = 0
      do jj = 1, N; do ii = 1, N
        pp = pp + 1
        owners(ii, jj) = min(ncpu-1, int(work_hk(pp)/max_key*dble(ncpu)))
      end do; end do
    else
      ! K-Section DD
      pp = 0
      do jj = 1, N; do ii = 1, N
        pp = pp + 1
        work_pts(1, pp) = (dble(ii)-0.5d0)/dble(N)
        work_pts(2, pp) = (dble(jj)-0.5d0)/dble(N)
        work_pts(3, pp) = (dble(kplane)-0.5d0)/dble(N)
      end do; end do
      call ks_dd(work_pts, np, work_dd)
      pp = 0
      do jj = 1, N; do ii = 1, N
        pp = pp + 1
        owners(ii, jj) = work_dd(pp)
      end do; end do
    end if
  end subroutine

  !-------------------------------------------------------------------
  ! Get owner CPU of the d-th face neighbor of cell (i,j) in current plane
  !-------------------------------------------------------------------
  function nbr_owner(i, j, d, planes, N, ip, ic, in) result(q)
    integer, intent(in) :: i, j, d, N, ip, ic, in
    integer, intent(in) :: planes(N, N, 0:2)
    integer :: q, ni, nj

    select case(d)
    case(1)  ! i+1
      ni = modulo(i, N) + 1
      q = planes(ni, j, ic)
    case(2)  ! i-1
      ni = modulo(i-2, N) + 1
      q = planes(ni, j, ic)
    case(3)  ! j+1
      nj = modulo(j, N) + 1
      q = planes(i, nj, ic)
    case(4)  ! j-1
      nj = modulo(j-2, N) + 1
      q = planes(i, nj, ic)
    case(5)  ! k+1
      q = planes(i, j, in)
    case(6)  ! k-1
      q = planes(i, j, ip)
    case default
      q = -1
    end select
  end function

  !-------------------------------------------------------------------
  ! Compute DD statistics using 3-plane sliding window
  ! For each cell owned by myid: accumulate weight, count ghost faces
  !-------------------------------------------------------------------
  subroutine compute_dd_stats(idd_type, N, nc, cx_in, cy_in, cz_in, &
                               amp_in, sigma_in, myid_in, ncpu_in, &
                               out_weight, out_ghost, out_owned, out_partners)
    integer, intent(in) :: idd_type, N, nc, myid_in, ncpu_in
    real(8), intent(in) :: cx_in(nc), cy_in(nc), cz_in(nc)
    real(8), intent(in) :: amp_in(nc), sigma_in(nc)
    real(8), intent(out) :: out_weight
    integer, intent(out) :: out_ghost, out_owned, out_partners

    integer, allocatable :: planes(:,:,:)
    integer :: idx_prev, idx_curr, idx_next, itmp
    integer :: i, j, k, kn, d, q
    logical :: is_ghost

    allocate(planes(N, N, 0:2))

    out_weight = 0.0d0
    out_ghost = 0
    out_owned = 0
    partner_flag(:) = 0

    ! Initialize 3-plane window: prev=plane(N), curr=plane(1)
    idx_prev = 0; idx_curr = 1; idx_next = 2
    call compute_plane_owner(idd_type, N, N, planes(:,:,idx_prev))
    call compute_plane_owner(idd_type, 1, N, planes(:,:,idx_curr))

    do k = 1, N
      kn = modulo(k, N) + 1
      call compute_plane_owner(idd_type, kn, N, planes(:,:,idx_next))

      do j = 1, N
        do i = 1, N
          if(planes(i, j, idx_curr) == myid_in) then
            ! My cell: accumulate weight
            out_owned = out_owned + 1
            out_weight = out_weight + cell_weight(i, j, k, N, cx_in, cy_in, cz_in, &
                                                   amp_in, sigma_in, nc)

            ! Check if this cell is a ghost cell (has a face neighbor on another rank)
            is_ghost = .false.
            do d = 1, 6
              q = nbr_owner(i, j, d, planes, N, idx_prev, idx_curr, idx_next)
              if(q /= myid_in) then
                is_ghost = .true.
                partner_flag(q) = 1
              end if
            end do
            if(is_ghost) out_ghost = out_ghost + 1
          end if
        end do
      end do

      ! Rotate plane indices
      itmp = idx_prev; idx_prev = idx_curr; idx_curr = idx_next; idx_next = itmp
    end do

    ! Count unique neighbor ranks
    out_partners = 0
    do i = 0, ncpu_in-1
      if(partner_flag(i) == 1) out_partners = out_partners + 1
    end do

    deallocate(planes)
  end subroutine

end program lb_bench
