!################################################################
!################################################################
!################################################################
!################################################################
module morton_hash
  !--------------------------------------------------------------
  ! Open-addressing hash table for Morton key → igrid mapping.
  ! Per-level hash tables: mort_table(1:nlevelmax)
  !
  ! Each table maps morton_key (integer(mkb)) → igrid (integer).
  ! Uses linear probing with power-of-2 capacity.
  ! Sentinel values: MHASH_EMPTY = -1, MHASH_DELETED = -2
  !--------------------------------------------------------------
  use morton_keys, only: mkb
  implicit none

  integer(mkb), parameter :: MHASH_EMPTY   = -1_mkb
  integer(mkb), parameter :: MHASH_DELETED = -2_mkb

  type :: morton_hash_table
     integer(mkb), allocatable :: keys(:)    ! Morton keys
     integer, allocatable      :: igrids(:)  ! Grid indices
     integer :: capacity  ! Always power of 2
     integer :: count     ! Number of entries
  end type morton_hash_table

  ! Global per-level hash tables
  type(morton_hash_table), allocatable :: mort_table(:)

  ! Per-grid level array: grid_level(igrid) = AMR level of grid igrid
  ! Maintained by make_grid_coarse/fine, kill_grid, and morton_hash_rebuild
  integer, allocatable :: grid_level(:)

contains

  !--------------------------------------------------------------
  ! Initialize hash table with given capacity (rounded up to power of 2)
  !--------------------------------------------------------------
  subroutine morton_hash_init(table, cap)
    type(morton_hash_table), intent(inout) :: table
    integer, intent(in) :: cap
    integer :: actual_cap

    ! Round up to next power of 2, minimum 16
    actual_cap = 16
    do while (actual_cap < cap)
       actual_cap = actual_cap * 2
    end do

    if (allocated(table%keys)) deallocate(table%keys)
    if (allocated(table%igrids)) deallocate(table%igrids)
    allocate(table%keys(actual_cap))
    allocate(table%igrids(actual_cap))
    table%keys = MHASH_EMPTY
    table%igrids = 0
    table%capacity = actual_cap
    table%count = 0
  end subroutine morton_hash_init

  !--------------------------------------------------------------
  ! Destroy hash table
  !--------------------------------------------------------------
  subroutine morton_hash_destroy(table)
    type(morton_hash_table), intent(inout) :: table

    if (allocated(table%keys)) deallocate(table%keys)
    if (allocated(table%igrids)) deallocate(table%igrids)
    table%capacity = 0
    table%count = 0
  end subroutine morton_hash_destroy

  !--------------------------------------------------------------
  ! Hash function: multiplicative hashing for 64-bit keys
  ! Returns 1-based slot index in [1, capacity]
  !--------------------------------------------------------------
  pure function morton_hash_func(key, capacity) result(slot)
    integer(mkb), intent(in) :: key
    integer, intent(in) :: capacity
    integer :: slot
    integer(mkb) :: h

    ! Multiplicative hash (Knuth's golden ratio method)
    h = key * 2654435761_mkb
    h = ieor(h, ishft(h, -16))
    h = h * (-7046029254386353131_mkb)  ! 0x9E3779B97F4A7C15
    h = ieor(h, ishft(h, -13))

    ! Map to [1, capacity]  (capacity is power of 2)
    slot = int(iand(h, int(capacity - 1, mkb))) + 1
  end function morton_hash_func

  !--------------------------------------------------------------
  ! Insert key → igrid into hash table
  ! If key already exists, update igrid
  !--------------------------------------------------------------
  subroutine morton_hash_insert(table, key, igrid)
    type(morton_hash_table), intent(inout) :: table
    integer(mkb), intent(in) :: key
    integer, intent(in) :: igrid
    integer :: slot, i

    ! Rehash if load factor > 0.7
    if (table%count * 10 > table%capacity * 7) then
       call morton_hash_rehash(table)
    end if

    slot = morton_hash_func(key, table%capacity)

    do i = 0, table%capacity - 1
       slot = int(iand(int(slot - 1 + i, mkb), int(table%capacity - 1, mkb))) + 1

       if (table%keys(slot) == MHASH_EMPTY .or. &
           table%keys(slot) == MHASH_DELETED) then
          ! Empty or deleted slot: insert here
          table%keys(slot) = key
          table%igrids(slot) = igrid
          table%count = table%count + 1
          return
       else if (table%keys(slot) == key) then
          ! Key already exists: update
          table%igrids(slot) = igrid
          return
       end if
    end do

    ! Table full (should not happen due to rehashing)
    write(*,*) 'morton_hash_insert: FATAL - hash table full!'
    stop
  end subroutine morton_hash_insert

  !--------------------------------------------------------------
  ! Delete key from hash table
  !--------------------------------------------------------------
  subroutine morton_hash_delete(table, key)
    type(morton_hash_table), intent(inout) :: table
    integer(mkb), intent(in) :: key
    integer :: slot, i

    slot = morton_hash_func(key, table%capacity)

    do i = 0, table%capacity - 1
       slot = int(iand(int(slot - 1 + i, mkb), int(table%capacity - 1, mkb))) + 1

       if (table%keys(slot) == MHASH_EMPTY) then
          ! Key not found
          return
       else if (table%keys(slot) == key) then
          ! Found: mark as deleted
          table%keys(slot) = MHASH_DELETED
          table%igrids(slot) = 0
          table%count = table%count - 1
          return
       end if
    end do
  end subroutine morton_hash_delete

  !--------------------------------------------------------------
  ! Lookup key in hash table
  ! Returns igrid if found, 0 if not found
  !--------------------------------------------------------------
  pure function morton_hash_lookup(table, key) result(igrid)
    type(morton_hash_table), intent(in) :: table
    integer(mkb), intent(in) :: key
    integer :: igrid
    integer :: slot, i

    igrid = 0
    if (table%capacity == 0) return

    slot = morton_hash_func(key, table%capacity)

    do i = 0, table%capacity - 1
       slot = int(iand(int(slot - 1 + i, mkb), int(table%capacity - 1, mkb))) + 1

       if (table%keys(slot) == MHASH_EMPTY) then
          ! Not found
          return
       else if (table%keys(slot) == key) then
          ! Found
          igrid = table%igrids(slot)
          return
       end if
       ! MHASH_DELETED: continue probing
    end do
  end function morton_hash_lookup

  !--------------------------------------------------------------
  ! Rehash: double the capacity and re-insert all entries
  !--------------------------------------------------------------
  subroutine morton_hash_rehash(table)
    type(morton_hash_table), intent(inout) :: table
    integer(mkb), allocatable :: old_keys(:)
    integer, allocatable :: old_igrids(:)
    integer :: old_cap, i

    old_cap = table%capacity
    allocate(old_keys(old_cap))
    allocate(old_igrids(old_cap))
    old_keys = table%keys
    old_igrids = table%igrids

    ! Double capacity
    deallocate(table%keys, table%igrids)
    table%capacity = old_cap * 2
    allocate(table%keys(table%capacity))
    allocate(table%igrids(table%capacity))
    table%keys = MHASH_EMPTY
    table%igrids = 0
    table%count = 0

    ! Re-insert all entries
    do i = 1, old_cap
       if (old_keys(i) /= MHASH_EMPTY .and. &
           old_keys(i) /= MHASH_DELETED) then
          call morton_hash_insert(table, old_keys(i), old_igrids(i))
       end if
    end do

    deallocate(old_keys, old_igrids)
  end subroutine morton_hash_rehash

end module morton_hash
