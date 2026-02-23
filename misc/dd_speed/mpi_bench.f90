program mpi_bench
  use mpi
  implicit none

  ! ---- parameters ----
  integer, parameter :: dp = kind(1.0d0)
  integer, parameter :: nsizes = 8
  integer            :: sizes(nsizes)
  data sizes / 1, 10, 100, 1000, 10000, 100000, 1000000, 10000000 /

  integer :: myid, ncpu, ierr, nrepeat
  character(len=32) :: arg

  ! ---- init ----
  call MPI_INIT(ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, ncpu, ierr)

  ! ---- parse command line ----
  nrepeat = 100
  if (command_argument_count() >= 1) then
     call get_command_argument(1, arg)
     read(arg, *) nrepeat
  end if

  if (myid == 0) then
     write(*,'(A)')        '=== MPI Communication Benchmark ==='
     write(*,'(A,I0,A,I0)') 'ncpu=', ncpu, ', nrepeat=', nrepeat
     write(*,*)
  end if

  call test_allreduce()
  call test_binary_tree_allreduce()
  call test_bcast()
  call test_neighbor_p2p()
  call test_alltoall()

  call MPI_FINALIZE(ierr)

contains

  !=====================================================================
  ! Test 1: MPI_ALLREDUCE
  !=====================================================================
  subroutine test_allreduce()
    implicit none
    integer :: isz, irep, i, idx
    real(dp), allocatable :: sbuf(:), rbuf(:)
    real(dp) :: t0, t1
    real(dp) :: times(nsizes), stds(nsizes)
    real(dp), allocatable :: tlog(:)

    if (myid == 0) write(*,'(A)') '--- Test 1: MPI_ALLREDUCE ---'

    allocate(tlog(nrepeat))

    do idx = 1, nsizes
       isz = sizes(idx)
       allocate(sbuf(isz), rbuf(isz))

       ! fill
       do i = 1, isz
          sbuf(i) = dble(myid + i)
       end do

       ! warmup
       call MPI_ALLREDUCE(sbuf, rbuf, isz, MPI_DOUBLE_PRECISION, &
            MPI_SUM, MPI_COMM_WORLD, ierr)

       ! timed
       do irep = 1, nrepeat
          do i = 1, isz
             sbuf(i) = dble(myid + i)
          end do
          call MPI_BARRIER(MPI_COMM_WORLD, ierr)
          t0 = MPI_WTIME()
          call MPI_ALLREDUCE(sbuf, rbuf, isz, MPI_DOUBLE_PRECISION, &
               MPI_SUM, MPI_COMM_WORLD, ierr)
          t1 = MPI_WTIME()
          tlog(irep) = t1 - t0
       end do

       call compute_stats(tlog, nrepeat, times(idx), stds(idx))

       deallocate(sbuf, rbuf)
    end do

    if (myid == 0) then
       call print_table(sizes, nsizes, times, stds)
       write(*,*)
    end if

    deallocate(tlog)
  end subroutine test_allreduce

  !=====================================================================
  ! Test 2: Binary Tree ALLREDUCE
  !=====================================================================
  subroutine test_binary_tree_allreduce()
    implicit none
    integer :: isz, irep, i, idx
    real(dp), allocatable :: sbuf(:), rbuf(:)
    real(dp) :: t0, t1
    real(dp) :: times(nsizes), stds(nsizes)
    real(dp), allocatable :: tlog(:)

    if (myid == 0) write(*,'(A)') '--- Test 2: Binary Tree ALLREDUCE ---'

    allocate(tlog(nrepeat))

    do idx = 1, nsizes
       isz = sizes(idx)
       allocate(sbuf(isz), rbuf(isz))

       ! warmup
       do i = 1, isz
          sbuf(i) = dble(myid + i)
       end do
       call binary_tree_allreduce(sbuf, rbuf, isz)

       ! timed
       do irep = 1, nrepeat
          do i = 1, isz
             sbuf(i) = dble(myid + i)
          end do
          call MPI_BARRIER(MPI_COMM_WORLD, ierr)
          t0 = MPI_WTIME()
          call binary_tree_allreduce(sbuf, rbuf, isz)
          t1 = MPI_WTIME()
          tlog(irep) = t1 - t0
       end do

       call compute_stats(tlog, nrepeat, times(idx), stds(idx))

       deallocate(sbuf, rbuf)
    end do

    if (myid == 0) then
       call print_table(sizes, nsizes, times, stds)
       write(*,*)
    end if

    deallocate(tlog)
  end subroutine test_binary_tree_allreduce

  !=====================================================================
  ! Test 3: MPI_BCAST
  !=====================================================================
  subroutine test_bcast()
    implicit none
    integer :: isz, irep, i, idx
    real(dp), allocatable :: buf(:)
    real(dp) :: t0, t1
    real(dp) :: times(nsizes), stds(nsizes)
    real(dp), allocatable :: tlog(:)

    if (myid == 0) write(*,'(A)') '--- Test 3: MPI_BCAST ---'

    allocate(tlog(nrepeat))

    do idx = 1, nsizes
       isz = sizes(idx)
       allocate(buf(isz))

       ! warmup
       if (myid == 0) then
          do i = 1, isz
             buf(i) = dble(i)
          end do
       end if
       call MPI_BCAST(buf, isz, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

       ! timed
       do irep = 1, nrepeat
          if (myid == 0) then
             do i = 1, isz
                buf(i) = dble(i + irep)
             end do
          end if
          call MPI_BARRIER(MPI_COMM_WORLD, ierr)
          t0 = MPI_WTIME()
          call MPI_BCAST(buf, isz, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
          t1 = MPI_WTIME()
          tlog(irep) = t1 - t0
       end do

       call compute_stats(tlog, nrepeat, times(idx), stds(idx))

       deallocate(buf)
    end do

    if (myid == 0) then
       call print_table(sizes, nsizes, times, stds)
       write(*,*)
    end if

    deallocate(tlog)
  end subroutine test_bcast

  !=====================================================================
  ! Test 4: Neighbor P2P Exchange (1D ring via MPI_SENDRECV)
  !=====================================================================
  subroutine test_neighbor_p2p()
    implicit none
    integer :: isz, irep, i, idx
    integer :: left, right
    real(dp), allocatable :: sbuf_r(:), rbuf_r(:)
    real(dp), allocatable :: sbuf_l(:), rbuf_l(:)
    integer :: stat(MPI_STATUS_SIZE)
    real(dp) :: t0, t1
    real(dp) :: times(nsizes), stds(nsizes)
    real(dp), allocatable :: tlog(:)

    if (myid == 0) write(*,'(A)') '--- Test 4: Neighbor P2P Exchange (1D ring) ---'

    allocate(tlog(nrepeat))

    ! 1D ring neighbors
    left  = mod(myid - 1 + ncpu, ncpu)
    right = mod(myid + 1, ncpu)

    do idx = 1, nsizes
       isz = sizes(idx)
       allocate(sbuf_r(isz), rbuf_r(isz), sbuf_l(isz), rbuf_l(isz))

       do i = 1, isz
          sbuf_r(i) = dble(myid * 1000 + i)
          sbuf_l(i) = dble(myid * 2000 + i)
       end do

       ! warmup
       call MPI_SENDRECV(sbuf_r, isz, MPI_DOUBLE_PRECISION, right, 100, &
            rbuf_r, isz, MPI_DOUBLE_PRECISION, left, 100, &
            MPI_COMM_WORLD, stat, ierr)
       call MPI_SENDRECV(sbuf_l, isz, MPI_DOUBLE_PRECISION, left, 200, &
            rbuf_l, isz, MPI_DOUBLE_PRECISION, right, 200, &
            MPI_COMM_WORLD, stat, ierr)

       ! timed
       do irep = 1, nrepeat
          call MPI_BARRIER(MPI_COMM_WORLD, ierr)
          t0 = MPI_WTIME()
          ! send right, recv from left
          call MPI_SENDRECV(sbuf_r, isz, MPI_DOUBLE_PRECISION, right, 100, &
               rbuf_r, isz, MPI_DOUBLE_PRECISION, left, 100, &
               MPI_COMM_WORLD, stat, ierr)
          ! send left, recv from right
          call MPI_SENDRECV(sbuf_l, isz, MPI_DOUBLE_PRECISION, left, 200, &
               rbuf_l, isz, MPI_DOUBLE_PRECISION, right, 200, &
               MPI_COMM_WORLD, stat, ierr)
          t1 = MPI_WTIME()
          tlog(irep) = t1 - t0
       end do

       call compute_stats(tlog, nrepeat, times(idx), stds(idx))

       deallocate(sbuf_r, rbuf_r, sbuf_l, rbuf_l)
    end do

    if (myid == 0) then
       call print_table(sizes, nsizes, times, stds)
       write(*,*)
    end if

    deallocate(tlog)
  end subroutine test_neighbor_p2p

  !=====================================================================
  ! Test 5: MPI_ALLTOALL (uniform)
  !=====================================================================
  subroutine test_alltoall()
    implicit none
    integer :: isz, irep, i, idx
    real(dp), allocatable :: sbuf(:), rbuf(:)
    real(dp) :: t0, t1
    real(dp) :: times(nsizes), stds(nsizes)
    real(dp), allocatable :: tlog(:)
    integer(8) :: mem_per_rank
    logical :: skip(nsizes)
    integer(8) :: one_gb

    if (myid == 0) write(*,'(A)') '--- Test 5: MPI_ALLTOALL ---'

    allocate(tlog(nrepeat))
    one_gb = 1073741824_8  ! 1 GB in bytes

    do idx = 1, nsizes
       isz = sizes(idx)
       ! Memory check: each rank needs isz*ncpu doubles for send + recv
       mem_per_rank = int(isz, 8) * int(ncpu, 8) * 8_8 * 2_8
       if (mem_per_rank > one_gb) then
          skip(idx) = .true.
          times(idx) = -1.0d0
          stds(idx)  = -1.0d0
          cycle
       end if
       skip(idx) = .false.

       allocate(sbuf(int(isz, 8) * int(ncpu, 8)))
       allocate(rbuf(int(isz, 8) * int(ncpu, 8)))

       do i = 1, isz * ncpu
          sbuf(i) = dble(myid * 1000 + i)
       end do

       ! warmup
       call MPI_ALLTOALL(sbuf, isz, MPI_DOUBLE_PRECISION, &
            rbuf, isz, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr)

       ! timed
       do irep = 1, nrepeat
          call MPI_BARRIER(MPI_COMM_WORLD, ierr)
          t0 = MPI_WTIME()
          call MPI_ALLTOALL(sbuf, isz, MPI_DOUBLE_PRECISION, &
               rbuf, isz, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr)
          t1 = MPI_WTIME()
          tlog(irep) = t1 - t0
       end do

       call compute_stats(tlog, nrepeat, times(idx), stds(idx))

       deallocate(sbuf, rbuf)
    end do

    if (myid == 0) then
       call print_table_alltoall(sizes, nsizes, times, stds, skip)
       write(*,*)
    end if

    deallocate(tlog)
  end subroutine test_alltoall

  !=====================================================================
  ! Binary tree allreduce: reduce to root + broadcast
  !=====================================================================
  subroutine binary_tree_allreduce(sbuf, rbuf, n)
    implicit none
    integer, intent(in) :: n
    real(dp), intent(inout) :: sbuf(n)
    real(dp), intent(out)   :: rbuf(n)

    integer :: stride, partner, i, ierr2
    integer :: stat(MPI_STATUS_SIZE)
    real(dp), allocatable :: tmp(:)

    allocate(tmp(n))

    ! --- Reduce phase (bottom-up) ---
    stride = 1
    do while (stride < ncpu)
       if (mod(myid, 2 * stride) == 0) then
          partner = myid + stride
          if (partner < ncpu) then
             call MPI_RECV(tmp, n, MPI_DOUBLE_PRECISION, partner, 301, &
                  MPI_COMM_WORLD, stat, ierr2)
             do i = 1, n
                sbuf(i) = sbuf(i) + tmp(i)
             end do
          end if
       else if (mod(myid, 2 * stride) == stride) then
          partner = myid - stride
          call MPI_SEND(sbuf, n, MPI_DOUBLE_PRECISION, partner, 301, &
               MPI_COMM_WORLD, ierr2)
       end if
       stride = stride * 2
    end do

    ! Root has full sum in sbuf
    if (myid == 0) then
       do i = 1, n
          rbuf(i) = sbuf(i)
       end do
    end if

    ! --- Broadcast phase (top-down) ---
    stride = 1
    do while (stride * 2 < ncpu)
       stride = stride * 2
    end do

    do while (stride >= 1)
       if (mod(myid, 2 * stride) == 0) then
          partner = myid + stride
          if (partner < ncpu) then
             call MPI_SEND(rbuf, n, MPI_DOUBLE_PRECISION, partner, 302, &
                  MPI_COMM_WORLD, ierr2)
          end if
       else if (mod(myid, 2 * stride) == stride) then
          partner = myid - stride
          call MPI_RECV(rbuf, n, MPI_DOUBLE_PRECISION, partner, 302, &
               MPI_COMM_WORLD, stat, ierr2)
       end if
       stride = stride / 2
    end do

    deallocate(tmp)
  end subroutine binary_tree_allreduce

  !=====================================================================
  ! Compute statistics: max across ranks per repeat, then mean +/- std
  !=====================================================================
  subroutine compute_stats(tlog, nrep, tmean, tstd)
    implicit none
    integer, intent(in)  :: nrep
    real(dp), intent(in) :: tlog(nrep)
    real(dp), intent(out) :: tmean, tstd

    real(dp) :: tmax(nrep), sumv, sumv2
    integer :: irep, ierr2

    call MPI_ALLREDUCE(tlog, tmax, nrep, MPI_DOUBLE_PRECISION, &
         MPI_MAX, MPI_COMM_WORLD, ierr2)

    sumv = 0.0d0
    do irep = 1, nrep
       sumv = sumv + tmax(irep)
    end do
    tmean = sumv / dble(nrep)

    sumv2 = 0.0d0
    do irep = 1, nrep
       sumv2 = sumv2 + (tmax(irep) - tmean)**2
    end do
    if (nrep > 1) then
       tstd = sqrt(sumv2 / dble(nrep - 1))
    else
       tstd = 0.0d0
    end if
  end subroutine compute_stats

  !=====================================================================
  ! Print table for Tests 1-4
  !=====================================================================
  subroutine print_table(sz, ns, times, stds)
    implicit none
    integer, intent(in) :: ns, sz(ns)
    real(dp), intent(in) :: times(ns), stds(ns)

    integer :: idx
    real(dp) :: time_us, bw_gbs, msg_bytes

    write(*,'(A12,A16,A12,A16)') 'Elements', 'Time(us)', 'Std(us)', 'BW(GB/s)'
    do idx = 1, ns
       time_us = times(idx) * 1.0d6
       msg_bytes = dble(sz(idx)) * 8.0d0
       if (times(idx) > 0.0d0) then
          bw_gbs = msg_bytes / times(idx) / 1.0d9
       else
          bw_gbs = 0.0d0
       end if
       write(*,'(I12, F16.1, F12.1, F16.4)') &
            sz(idx), time_us, stds(idx)*1.0d6, bw_gbs
    end do
  end subroutine print_table

  !=====================================================================
  ! Print table for Test 5 (Alltoall) with skip handling
  !=====================================================================
  subroutine print_table_alltoall(sz, ns, times, stds, skip)
    implicit none
    integer, intent(in) :: ns, sz(ns)
    real(dp), intent(in) :: times(ns), stds(ns)
    logical, intent(in) :: skip(ns)

    integer :: idx
    real(dp) :: time_us, bw_gbs, msg_bytes

    write(*,'(A12,A16,A12,A16)') 'Elements', 'Time(us)', 'Std(us)', 'BW(GB/s)'
    do idx = 1, ns
       if (skip(idx)) then
          write(*,'(I12, A16)') sz(idx), '       [skipped]'
          cycle
       end if
       time_us = times(idx) * 1.0d6
       msg_bytes = dble(sz(idx)) * dble(ncpu) * 8.0d0
       if (times(idx) > 0.0d0) then
          bw_gbs = msg_bytes / times(idx) / 1.0d9
       else
          bw_gbs = 0.0d0
       end if
       write(*,'(I12, F16.1, F12.1, F16.4)') &
            sz(idx), time_us, stds(idx)*1.0d6, bw_gbs
    end do
  end subroutine print_table_alltoall

end program mpi_bench
