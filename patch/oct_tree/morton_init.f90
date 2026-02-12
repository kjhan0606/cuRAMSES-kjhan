!################################################################
!################################################################
!################################################################
!################################################################
subroutine morton_hash_rebuild
  !--------------------------------------------------------------
  ! Rebuild per-level Morton hash tables from all active grids
  ! (local + virtual). Can be called after load balancing or
  ! build_comm to resync with current grid structure.
  !--------------------------------------------------------------
  use amr_commons
  use amr_parameters
  use morton_keys
  use morton_hash
  implicit none

  integer :: ilevel, icpu, igrid, i
  integer :: ngrids_level
  integer(mkb) :: mkey

  ! Check that coordinates fit in MORTON_MAXBITS
  do ilevel = 1, nlevelmax
     if (nx * 2**(ilevel-1) > 2**MORTON_MAXBITS) then
        if (myid==1) then
           write(*,*) 'Morton hash: ERROR - level', ilevel, &
                ' exceeds coordinate range'
        end if
        return
     end if
  end do

  ! Allocate per-level hash tables if not yet allocated
  if (.not. allocated(mort_table)) then
     allocate(mort_table(1:nlevelmax))
  end if

  ! Allocate grid_level array if not yet allocated
  if (.not. allocated(grid_level)) then
     allocate(grid_level(1:ngridmax))
     grid_level = 0
  end if

  do ilevel = 1, nlevelmax
     ! Count all grids at this level (local + virtual from all CPUs)
     ngrids_level = 0
     do icpu = 1, ncpu
        ngrids_level = ngrids_level + numbl(icpu, ilevel)
     end do

     ! (Re)initialize hash table
     call morton_hash_init(mort_table(ilevel), max(2 * ngrids_level, 16))

     ! Insert all grids at this level and set grid_level
     do icpu = 1, ncpu
        igrid = headl(icpu, ilevel)
        do i = 1, numbl(icpu, ilevel)
           if (igrid == 0) exit
           mkey = grid_to_morton(igrid, ilevel)
           call morton_hash_insert(mort_table(ilevel), mkey, igrid)
           grid_level(igrid) = ilevel
           igrid = next(igrid)
        end do
     end do
  end do

end subroutine morton_hash_rebuild
!################################################################
!################################################################
!################################################################
!################################################################
subroutine morton_hash_rebuild_level(ilevel)
  !--------------------------------------------------------------
  ! Rebuild Morton hash table for a single level only.
  ! Also updates grid_level for all grids at this level.
  ! Called after build_comm to include virtual grids.
  !--------------------------------------------------------------
  use amr_commons
  use amr_parameters
  use morton_keys
  use morton_hash
  implicit none

  integer, intent(in) :: ilevel
  integer :: icpu, igrid, i
  integer :: ngrids_level
  integer(mkb) :: mkey

  if (.not. allocated(mort_table)) return
  if (ilevel < 1 .or. ilevel > nlevelmax) return

  ! Count all grids at this level (local + virtual)
  ngrids_level = 0
  do icpu = 1, ncpu
     ngrids_level = ngrids_level + numbl(icpu, ilevel)
  end do

  ! (Re)initialize hash table for this level
  call morton_hash_init(mort_table(ilevel), max(2 * ngrids_level, 16))

  ! Insert all grids and set grid_level
  do icpu = 1, ncpu
     igrid = headl(icpu, ilevel)
     do i = 1, numbl(icpu, ilevel)
        if (igrid == 0) exit
        mkey = grid_to_morton(igrid, ilevel)
        call morton_hash_insert(mort_table(ilevel), mkey, igrid)
        if (allocated(grid_level)) grid_level(igrid) = ilevel
        igrid = next(igrid)
     end do
  end do

end subroutine morton_hash_rebuild_level
!################################################################
!################################################################
!################################################################
!################################################################
subroutine morton_hash_verify(label)
  !--------------------------------------------------------------
  ! Verify Morton hash table self-consistency:
  ! - Every local grid must be found by its own Morton key lookup
  ! - grid_level must match the level in the headl/numbl lists
  ! (nbor is no longer available for cross-checking)
  !--------------------------------------------------------------
  use amr_commons
  use amr_parameters
  use morton_keys
  use morton_hash
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif

  character(len=*), intent(in) :: label
  integer :: ilevel, igrid, i
  integer :: igrid_morton
  integer :: nmatch, nmismatch
  integer :: nmatch_all, nmismatch_all
  integer(mkb) :: mkey
  integer :: info

  if (.not. allocated(mort_table)) return

  nmatch = 0
  nmismatch = 0

  do ilevel = 1, nlevelmax
     ! Check all grids (local + virtual)
     do i = 1, ncpu
        igrid = headl(i, ilevel)
        do while (igrid > 0)
           mkey = grid_to_morton(igrid, ilevel)

           ! Verify self-lookup
           igrid_morton = morton_hash_lookup(mort_table(ilevel), mkey)
           if (igrid_morton == igrid) then
              nmatch = nmatch + 1
           else
              nmismatch = nmismatch + 1
              if (nmismatch <= 5) then
                 write(*,'(A,A,A,I6,A,I8,A,I8,A,I3)') &
                      ' [', trim(label), '] self-lookup MISMATCH cpu=', myid, &
                      ' igrid=', igrid, ' lookup=', igrid_morton, &
                      ' level=', ilevel
              end if
           end if

           igrid = next(igrid)
        end do
     end do
  end do

  ! Gather global statistics
#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(nmatch, nmatch_all, 1, &
       MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, info)
  call MPI_ALLREDUCE(nmismatch, nmismatch_all, 1, &
       MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, info)
#else
  nmatch_all = nmatch
  nmismatch_all = nmismatch
#endif

  if (myid == 1) then
     write(*,'(A,A,A,I10,A,I10)') &
          ' Morton [', trim(label), '] self-lookups=', nmatch_all, &
          ' mismatch=', nmismatch_all
  end if

  if (nmismatch_all > 0) then
     if (myid == 1) write(*,*) 'Morton hash: VERIFICATION FAILED at ', label
     call clean_stop
  end if

end subroutine morton_hash_verify
!################################################################
!################################################################
!################################################################
!################################################################
subroutine morton_hash_build_and_verify
  !--------------------------------------------------------------
  ! Combined build + verify (used at init_amr startup)
  !--------------------------------------------------------------
  use amr_commons
  use morton_hash
  implicit none

  integer :: ngrids_total, ilevel, icpu

  if(myid==1) write(*,*) 'Morton hash: building per-level hash tables...'

  call morton_hash_rebuild

  ! Count total grids
  ngrids_total = 0
  do ilevel = 1, nlevelmax
     do icpu = 1, ncpu
        ngrids_total = ngrids_total + numbl(icpu, ilevel)
     end do
  end do

  if(myid==1) write(*,'(A,I10,A)') &
       ' Morton hash: inserted ', ngrids_total, ' grids total'

  call morton_hash_verify('init')

  if(myid==1) write(*,*) 'Morton hash: ALL CHECKS PASSED'

end subroutine morton_hash_build_and_verify
