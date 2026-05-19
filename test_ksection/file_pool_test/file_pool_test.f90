! ----------------------------------------------------------------------
! file_pool_test — unit test for Phase D file_pool LRU cache
!
! Creates 5 trivial HDF5 files, exercises the pool with a scripted access
! sequence, and verifies (hits, misses, evictions, n_open, handle reuse)
! against expected values. Single-rank MPI program (suffices for the LRU
! logic; collective semantics inside HDF5 still go through the parallel
! plist setup path used by file_pool_get).
!
! Run:  mpirun -n 1 ./file_pool_test
! ----------------------------------------------------------------------
program file_pool_test
  use hdf5
  use ramses_hdf5_io
  use amr_parameters, only: i8b
  implicit none
  include 'mpif.h'

  integer, parameter :: NFILE = 5
  character(len=64)  :: fname(NFILE)
  integer            :: ierr, i, ifail
  integer(HID_T)     :: fid, fid_a, fid_b
  integer            :: n_open, cap
  integer(i8b)       :: n_hits, n_misses, n_evictions

  ! Expected counters at each checkpoint
  integer            :: e_open
  integer(i8b)       :: e_hits, e_misses, e_evictions

  ifail = 0
  call MPI_Init(ierr)

  ! ---- step 0: create 5 trivial HDF5 files (each just an empty root) ----
  call h5open_f(ierr)
  do i = 1, NFILE
     write(fname(i),'("fp_test_",I0,".h5")') i
     call create_empty_file(trim(fname(i)))
  end do

  ! ---- step 1: init pool with cap=3 ----
  call file_pool_init(3)
  call file_pool_stats(n_open, n_hits, n_misses, n_evictions, cap)
  call check_eq_i('cap after init',   cap,    3, ifail)
  call check_eq_i('n_open after init',n_open, 0, ifail)

  ! ---- step 2: get file1, file2, file3 (3 cold misses, no evictions) ----
  call file_pool_get(trim(fname(1)), MPI_COMM_WORLD, fid)
  call file_pool_get(trim(fname(2)), MPI_COMM_WORLD, fid)
  call file_pool_get(trim(fname(3)), MPI_COMM_WORLD, fid)
  call file_pool_stats(n_open, n_hits, n_misses, n_evictions, cap)
  e_open=3; e_hits=0_i8b; e_misses=3_i8b; e_evictions=0_i8b
  call check_eq_i ('A: n_open',      n_open,      e_open,      ifail)
  call check_eq_i8('A: n_hits',      n_hits,      e_hits,      ifail)
  call check_eq_i8('A: n_misses',    n_misses,    e_misses,    ifail)
  call check_eq_i8('A: n_evictions', n_evictions, e_evictions, ifail)

  ! ---- step 3: re-get file1 → cache hit, file1 bumps to most-recent ----
  call file_pool_get(trim(fname(1)), MPI_COMM_WORLD, fid_a)
  call file_pool_get(trim(fname(1)), MPI_COMM_WORLD, fid_b)
  call check_eq_hid('B: handle reuse fid_a==fid_b', fid_a, fid_b, ifail)
  call file_pool_stats(n_open, n_hits, n_misses, n_evictions, cap)
  e_open=3; e_hits=2_i8b; e_misses=3_i8b; e_evictions=0_i8b
  call check_eq_i ('B: n_open',      n_open,      e_open,      ifail)
  call check_eq_i8('B: n_hits',      n_hits,      e_hits,      ifail)
  call check_eq_i8('B: n_misses',    n_misses,    e_misses,    ifail)
  call check_eq_i8('B: n_evictions', n_evictions, e_evictions, ifail)

  ! ---- step 4: get file4 → miss, must evict file2 (LRU at this point) ----
  call file_pool_get(trim(fname(4)), MPI_COMM_WORLD, fid)
  call file_pool_stats(n_open, n_hits, n_misses, n_evictions, cap)
  e_open=3; e_hits=2_i8b; e_misses=4_i8b; e_evictions=1_i8b
  call check_eq_i ('C: n_open',      n_open,      e_open,      ifail)
  call check_eq_i8('C: n_hits',      n_hits,      e_hits,      ifail)
  call check_eq_i8('C: n_misses',    n_misses,    e_misses,    ifail)
  call check_eq_i8('C: n_evictions', n_evictions, e_evictions, ifail)
  ! Touching file3 must still hit (file3 was evicted? no — file2 was LRU).
  ! Pool now holds {file1, file3, file4}. Re-getting file3 must be a hit.
  call file_pool_get(trim(fname(3)), MPI_COMM_WORLD, fid)
  call file_pool_stats(n_open, n_hits, n_misses, n_evictions, cap)
  e_open=3; e_hits=3_i8b; e_misses=4_i8b; e_evictions=1_i8b
  call check_eq_i8('D: n_hits after re-get file3', n_hits, e_hits, ifail)

  ! ---- step 5: re-get file2 → MISS (was evicted in step 4), evicts file1
  !   Tick trace at this point: file1.lru=5 (last bumped in B),
  !   file4.lru=6 (added in step 4), file3.lru=7 (re-bumped in step 4).
  !   So LRU = file1 (smallest tick) → file1 is the victim.
  call file_pool_get(trim(fname(2)), MPI_COMM_WORLD, fid)
  call file_pool_stats(n_open, n_hits, n_misses, n_evictions, cap)
  e_open=3; e_hits=3_i8b; e_misses=5_i8b; e_evictions=2_i8b
  call check_eq_i ('E: n_open',      n_open,      e_open,      ifail)
  call check_eq_i8('E: n_hits',      n_hits,      e_hits,      ifail)
  call check_eq_i8('E: n_misses',    n_misses,    e_misses,    ifail)
  call check_eq_i8('E: n_evictions', n_evictions, e_evictions, ifail)
  ! Confirm file1 is gone: re-get must be a miss (evicts file4 now LRU)
  call file_pool_get(trim(fname(1)), MPI_COMM_WORLD, fid)
  call file_pool_stats(n_open, n_hits, n_misses, n_evictions, cap)
  e_open=3; e_hits=3_i8b; e_misses=6_i8b; e_evictions=3_i8b
  call check_eq_i8('F: n_misses after re-get file1', n_misses, e_misses, ifail)
  call check_eq_i8('F: n_evictions',                 n_evictions, e_evictions, ifail)

  ! ---- step 6: capacity stays at 3 across heavy churn ----
  do i = 1, 50
     call file_pool_get(trim(fname(1 + mod(i,NFILE))), MPI_COMM_WORLD, fid)
  end do
  call file_pool_stats(n_open, n_hits, n_misses, n_evictions, cap)
  call check_eq_i ('G: n_open <= cap', n_open, 3, ifail)
  call check_eq_i ('G: cap unchanged', cap,    3, ifail)

  ! ---- step 7: clean shutdown ----
  call file_pool_close_all()
  call file_pool_stats(n_open, n_hits, n_misses, n_evictions, cap)
  call check_eq_i('H: n_open after close_all', n_open, 0, ifail)
  call check_eq_i('H: cap after close_all',    cap,    0, ifail)

  ! ---- step 8: re-init must work (idempotent) ----
  call file_pool_init(2)
  call file_pool_get(trim(fname(1)), MPI_COMM_WORLD, fid)
  call file_pool_get(trim(fname(2)), MPI_COMM_WORLD, fid)
  call file_pool_get(trim(fname(3)), MPI_COMM_WORLD, fid)  ! evict file1
  call file_pool_stats(n_open, n_hits, n_misses, n_evictions, cap)
  call check_eq_i ('I: re-init cap=2',           cap,    2, ifail)
  call check_eq_i ('I: re-init n_open=2',        n_open, 2, ifail)
  call check_eq_i8('I: re-init n_evictions=1',   n_evictions, 1_i8b, ifail)
  call file_pool_close_all()

  ! ---- step 9: cleanup files ----
  do i = 1, NFILE
     call delete_file(trim(fname(i)))
  end do
  call h5close_f(ierr)

  if(ifail == 0) then
     write(*,'(A)') ''
     write(*,'(A)') '========================================'
     write(*,'(A)') ' file_pool_test: PASS (all assertions ok)'
     write(*,'(A)') '========================================'
  else
     write(*,'(A,I0,A)') 'file_pool_test: FAIL (', ifail, ' assertion(s) failed)'
  end if

  call MPI_Finalize(ierr)
  if(ifail /= 0) call exit(1)

contains

  subroutine create_empty_file(filename)
    character(len=*), intent(in) :: filename
    integer(HID_T) :: fid_loc, plist_loc
    integer        :: ierr_loc
    call h5pcreate_f(H5P_FILE_ACCESS_F, plist_loc, ierr_loc)
    call h5pset_fapl_mpio_f(plist_loc, MPI_COMM_WORLD, MPI_INFO_NULL, ierr_loc)
    call h5fcreate_f(trim(filename), H5F_ACC_TRUNC_F, fid_loc, ierr_loc, &
         access_prp=plist_loc)
    call h5pclose_f(plist_loc, ierr_loc)
    call h5fclose_f(fid_loc, ierr_loc)
  end subroutine create_empty_file

  subroutine delete_file(filename)
    character(len=*), intent(in) :: filename
    integer :: u, ios
    open(newunit=u, file=trim(filename), status='old', iostat=ios)
    if(ios == 0) close(u, status='delete')
  end subroutine delete_file

  subroutine check_eq_i(label, got, expected, ifailcnt)
    character(len=*), intent(in) :: label
    integer,          intent(in) :: got, expected
    integer,          intent(inout) :: ifailcnt
    if(got == expected) then
       write(*,'("  [PASS] ",A,"  =",I0)') label, got
    else
       write(*,'("  [FAIL] ",A,"  got=",I0," expected=",I0)') label, got, expected
       ifailcnt = ifailcnt + 1
    end if
  end subroutine check_eq_i

  subroutine check_eq_i8(label, got, expected, ifailcnt)
    character(len=*), intent(in) :: label
    integer(i8b),     intent(in) :: got, expected
    integer,          intent(inout) :: ifailcnt
    if(got == expected) then
       write(*,'("  [PASS] ",A,"  =",I0)') label, got
    else
       write(*,'("  [FAIL] ",A,"  got=",I0," expected=",I0)') label, got, expected
       ifailcnt = ifailcnt + 1
    end if
  end subroutine check_eq_i8

  subroutine check_eq_hid(label, a, b, ifailcnt)
    character(len=*), intent(in) :: label
    integer(HID_T),   intent(in) :: a, b
    integer,          intent(inout) :: ifailcnt
    if(a == b) then
       write(*,'("  [PASS] ",A,"  =",I0)') label, a
    else
       write(*,'("  [FAIL] ",A,"  a=",I0," b=",I0)') label, a, b
       ifailcnt = ifailcnt + 1
    end if
  end subroutine check_eq_hid

end program file_pool_test
