! Test program for FoF: Sequential vs Tree vs OMP+MPI Tree
! Build: mpiifx -O2 -qopenmp -fpp -D_KJHAN_TEST -DNDIM=3 -o test_fof test_fof.f90 ../patch/Horizon5-master-2/kjhan.f90
! Run:   OMP_NUM_THREADS=N mpirun -np P ./test_fof

program test_fof
  use AA
  use omp_lib
  implicit none
  include 'mpif.h'

  integer, parameter :: NTESTS = 12
  integer :: test_nsink(NTESTS)
  real(dp) :: test_foflink(NTESTS)
  integer :: itest, nsink, i, j, seed_size, ierr, myid_mpi, nprocs
  integer, allocatable :: seed(:)
  real(dp) :: foflink, dx_min2, factG
  integer, allocatable :: psink_seq(:), gsink_seq(:)
  integer, allocatable :: psink_tree(:), gsink_tree(:)
  integer, allocatable :: psink_omp(:), gsink_omp(:)
  integer, allocatable :: psink_auto(:), gsink_auto(:)
  integer :: igrp_seq, igrp_tree, igrp_omp, igrp_auto
  integer :: nfail, nthr
  real(dp) :: t1, t2, t3, t4, t5, t6, t7, t8
  integer, external :: Do_FoF_Sequential, Do_FoF_Tree, Omp_Do_FoF_Tree, Auto_Do_FoF
  logical :: match_tree, match_omp, match_auto

  call MPI_INIT(ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD, myid_mpi, ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr)

  test_nsink   = (/ 1, 2, 50, 100, 200, 500, 1000, 2000, 5000, 10000, 50000, 100000 /)
  test_foflink = (/ 0.50d0, 0.30d0, 0.10d0, 0.08d0, 0.06d0, 0.05d0, 0.04d0, 0.03d0, 0.02d0, 0.015d0, 0.007d0, 0.005d0 /)

  nthr = 1
  !$ nthr = omp_get_max_threads()

  if(myid_mpi == 0) then
    write(*,'(A,I0,A,I0,A)') '=== FoF Benchmark (', nprocs, ' MPI x ', nthr, ' threads) ==='
    write(*,*)
    write(*,'(A)') '  nsink  foflink  seq_grp tree_grp  omp_grp auto_grp   seq(ms)  tree(ms)   omp(ms)  auto(ms) tree_ok  omp_ok auto_ok'
    write(*,'(A)') '  -----  -------  ------- --------  ------- --------  --------  --------  --------  -------- -------  ------ -------'
  endif

  nfail = 0

  do itest = 1, NTESTS
    nsink = test_nsink(itest)
    foflink = test_foflink(itest)

    if(allocated(xsink)) deallocate(xsink)
    allocate(xsink(1:nsink, 1:3))

    ! Generate reproducible random particles with clusters
    call random_seed(size=seed_size)
    allocate(seed(seed_size))
    seed = 42 + itest * 7
    call random_seed(put=seed)
    deallocate(seed)

    ! First half: random positions
    do i = 1, nsink/2
      call random_number(xsink(i,1))
      call random_number(xsink(i,2))
      call random_number(xsink(i,3))
    enddo

    ! Second half: clusters (groups of 4 particles close together)
    i = nsink/2 + 1
    do while(i <= nsink)
      call random_number(xsink(i,1))
      call random_number(xsink(i,2))
      call random_number(xsink(i,3))
      do j = 1, min(3, nsink - i)
        xsink(i+j,1) = xsink(i,1) + dble(j)*0.3d0*foflink
        xsink(i+j,2) = xsink(i,2) + dble(j)*0.2d0*foflink
        xsink(i+j,3) = xsink(i,3) + dble(j)*0.1d0*foflink
        xsink(i+j,1) = modulo(xsink(i+j,1), 1.0d0)
        xsink(i+j,2) = modulo(xsink(i+j,2), 1.0d0)
        xsink(i+j,3) = modulo(xsink(i+j,3), 1.0d0)
      enddo
      i = i + 4
    enddo

    dx_min2 = foflink * foflink / (rmerge * rmerge)
    factG = 1.0d0

    allocate(psink_seq(1:nsink), gsink_seq(1:nsink))
    allocate(psink_tree(1:nsink), gsink_tree(1:nsink))
    allocate(psink_omp(1:nsink), gsink_omp(1:nsink))
    allocate(psink_auto(1:nsink), gsink_auto(1:nsink))

    ! --- Method 1: Sequential brute-force FoF ---
    do i = 1, nsink
      psink_seq(i) = i
      gsink_seq(i) = 0
    enddo
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    t1 = MPI_WTIME()
    igrp_seq = Do_FoF_Sequential(nsink, psink_seq, gsink_seq, dx_min2, factG)
    t2 = MPI_WTIME()

    ! --- Method 2: Serial tree-based Do_FoF_Tree ---
    do i = 1, nsink
      psink_tree(i) = i
      gsink_tree(i) = 0
    enddo
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    t3 = MPI_WTIME()
    igrp_tree = Do_FoF_Tree(nsink, psink_tree, gsink_tree, dx_min2, factG)
    t4 = MPI_WTIME()

    ! --- Method 3: OMP+MPI Omp_Do_FoF_Tree ---
    do i = 1, nsink
      psink_omp(i) = i
      gsink_omp(i) = 0
    enddo
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    t5 = MPI_WTIME()
    igrp_omp = Omp_Do_FoF_Tree(nsink, psink_omp, gsink_omp, dx_min2, factG)
    t6 = MPI_WTIME()

    ! --- Method 4: Auto_Do_FoF ---
    do i = 1, nsink
      psink_auto(i) = i
      gsink_auto(i) = 0
    enddo
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    t7 = MPI_WTIME()
    igrp_auto = Auto_Do_FoF(nsink, psink_auto, gsink_auto, dx_min2, factG)
    t8 = MPI_WTIME()

    ! --- Compare partitions (O(n) via group mapping) ---
    match_tree = .true.
    match_omp = .true.
    match_auto = .true.
    if(myid_mpi == 0) then
      call check_partition_match(nsink, igrp_seq, gsink_seq, igrp_tree, gsink_tree, match_tree)
      call check_partition_match(nsink, igrp_seq, gsink_seq, igrp_omp, gsink_omp, match_omp)
      call check_partition_match(nsink, igrp_seq, gsink_seq, igrp_auto, gsink_auto, match_auto)
    endif

    if(myid_mpi == 0) then
      write(*,'(I7,F9.4,I9,I9,I9,I9,F10.3,F10.3,F10.3,F10.3,A,A,A)') &
        nsink, foflink, igrp_seq, igrp_tree, igrp_omp, igrp_auto, &
        (t2-t1)*1000.0d0, (t4-t3)*1000.0d0, (t6-t5)*1000.0d0, (t8-t7)*1000.0d0, &
        merge('     OK','   FAIL', match_tree), &
        merge('     OK','   FAIL', match_omp), &
        merge('     OK','   FAIL', match_auto)
      if(.not. match_tree) nfail = nfail + 1
      if(.not. match_omp) nfail = nfail + 1
      if(.not. match_auto) nfail = nfail + 1
    endif

    deallocate(psink_seq, gsink_seq)
    deallocate(psink_tree, gsink_tree)
    deallocate(psink_omp, gsink_omp)
    deallocate(psink_auto, gsink_auto)
  enddo

  if(myid_mpi == 0) then
    write(*,*)
    if(nfail == 0) then
      write(*,'(A)') 'All tests PASSED!'
    else
      write(*,'(A,I0,A)') 'FAILED: ', nfail, ' tests'
    endif
  endif

  call MPI_FINALIZE(ierr)

contains

  subroutine check_partition_match(n, ngrp_a, gsink_a, ngrp_b, gsink_b, match)
    ! O(n) partition equivalence check via group mapping
    implicit none
    integer, intent(in) :: n, ngrp_a, ngrp_b
    integer, intent(in) :: gsink_a(n), gsink_b(n)
    logical, intent(out) :: match
    integer, allocatable :: map_a2b(:), map_b2a(:)
    integer :: i, ga, gb

    match = .true.
    if(ngrp_a == 0 .and. ngrp_b == 0) return

    allocate(map_a2b(1:max(ngrp_a,1)), map_b2a(1:max(ngrp_b,1)))
    map_a2b = 0
    map_b2a = 0

    do i = 1, n
      ga = gsink_a(i)
      gb = gsink_b(i)
      if(ga < 1 .or. gb < 1) then
        match = .false.; exit
      endif
      ! Check a->b mapping consistency
      if(map_a2b(ga) == 0) then
        map_a2b(ga) = gb
      else if(map_a2b(ga) /= gb) then
        match = .false.; exit
      endif
      ! Check b->a mapping consistency
      if(map_b2a(gb) == 0) then
        map_b2a(gb) = ga
      else if(map_b2a(gb) /= ga) then
        match = .false.; exit
      endif
    enddo

    deallocate(map_a2b, map_b2a)
  end subroutine

end program test_fof
