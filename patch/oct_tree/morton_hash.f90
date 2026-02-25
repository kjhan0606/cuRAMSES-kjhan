!################################################################
!################################################################
!################################################################
!################################################################
module morton_hash
  !--------------------------------------------------------------
  ! Open-addressing hash table for Morton key → igrid mapping.
  ! Per-level hash tables: mort_table(1:nlevelmax)
  !
  ! Each table maps type(mkey_t) → igrid (integer).
  ! Uses linear probing with power-of-2 capacity.
  ! Sentinel values: MHASH_EMPTY = MK_MINUS1, MHASH_DELETED = MK_MINUS2
  !--------------------------------------------------------------
  use morton_keys, only: mkey_t, MK_ZERO, MK_MINUS1, MK_MINUS2, &
       MORTON_MAXBITS, grid_to_morton, morton_neighbor, &
       morton_decode, morton_encode, &
       operator(==), operator(/=)
  implicit none

  type(mkey_t), parameter :: MHASH_EMPTY   = MK_MINUS1
  type(mkey_t), parameter :: MHASH_DELETED = MK_MINUS2

  type :: morton_hash_table
     type(mkey_t), allocatable :: keys(:)    ! Morton keys
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
  ! Hash function: multiplicative hashing with XOR mixing
  ! Returns 1-based slot index in [1, capacity]
  !--------------------------------------------------------------
  pure function morton_hash_func(key, capacity) result(slot)
    type(mkey_t), intent(in) :: key
    integer, intent(in) :: capacity
    integer :: slot
    integer(8) :: h

#ifdef MORTON128
    ! XOR fold 128-bit key to 64-bit
    h = ieor(key%lo, key%hi)
#else
    h = key%lo
#endif

    ! Multiplicative hash (Knuth's golden ratio method)
    h = h * 2654435761_8
    h = ieor(h, ishft(h, -16))
    h = h * (-7046029254386353131_8)  ! 0x9E3779B97F4A7C15
    h = ieor(h, ishft(h, -13))

    ! Map to [1, capacity]  (capacity is power of 2)
    slot = int(iand(h, int(capacity - 1, 8))) + 1
  end function morton_hash_func

  !--------------------------------------------------------------
  ! Insert key → igrid into hash table
  ! If key already exists, update igrid
  !--------------------------------------------------------------
  subroutine morton_hash_insert(table, key, igrid)
    type(morton_hash_table), intent(inout) :: table
    type(mkey_t), intent(in) :: key
    integer, intent(in) :: igrid
    integer :: slot, i

    ! Rehash if load factor > 0.7
    if (table%count * 10 > table%capacity * 7) then
       call morton_hash_rehash(table)
    end if

    slot = morton_hash_func(key, table%capacity)

    do i = 0, table%capacity - 1
       slot = iand(slot - 1 + i, table%capacity - 1) + 1

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
    type(mkey_t), intent(in) :: key
    integer :: slot, i

    slot = morton_hash_func(key, table%capacity)

    do i = 0, table%capacity - 1
       slot = iand(slot - 1 + i, table%capacity - 1) + 1

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
    type(mkey_t), intent(in) :: key
    integer :: igrid
    integer :: slot, i

    igrid = 0
    if (table%capacity == 0) return

    slot = morton_hash_func(key, table%capacity)

    do i = 0, table%capacity - 1
       slot = iand(slot - 1 + i, table%capacity - 1) + 1

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
    type(mkey_t), allocatable :: old_keys(:)
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

  !--------------------------------------------------------------
  ! morton_nbor_grid: replaces son(nbor(igrid, j))
  ! Returns same-level neighbor grid index, or 0 if not found
  !--------------------------------------------------------------
  function morton_nbor_grid(igrid, ilevel, j) result(igridn)
    use amr_commons, only: nx, ny, nz
    integer, intent(in) :: igrid, ilevel, j
    integer :: igridn
    integer(8) :: nmax_x, nmax_y, nmax_z
    type(mkey_t) :: mkey, nkey

    nmax_x = int(nx, 8) * 2_8**(ilevel-1)
    nmax_y = int(ny, 8) * 2_8**(ilevel-1)
    nmax_z = int(nz, 8) * 2_8**(ilevel-1)

    mkey = grid_to_morton(igrid, ilevel)
    nkey = morton_neighbor(mkey, j, nmax_x, nmax_y, nmax_z)
    if (nkey == MK_MINUS1) then
       igridn = 0
    else
       igridn = morton_hash_lookup(mort_table(ilevel), nkey)
    end if
  end function morton_nbor_grid

  !--------------------------------------------------------------
  ! morton_nbor_cell: replaces nbor(igrid, j)
  ! Returns father cell index of the neighbor in direction j
  ! Level 1: coarse cell index (1..nx*ny*nz)
  ! Level >= 2: fine cell = ncoarse + (ind_oct-1)*ngridmax + igrid_parent
  !--------------------------------------------------------------
  function morton_nbor_cell(igrid, ilevel, j) result(icell)
    use amr_commons, only: nx, ny, nz, ncoarse, ngridmax
    integer, intent(in) :: igrid, ilevel, j
    integer :: icell
    integer(8) :: nmax_x, nmax_y, nmax_z
    integer(8) :: nix, niy, niz, pix, piy, piz
    integer :: ind_oct, igrid_parent
    type(mkey_t) :: mkey, nkey, pkey

    nmax_x = int(nx, 8) * 2_8**(ilevel-1)
    nmax_y = int(ny, 8) * 2_8**(ilevel-1)
    nmax_z = int(nz, 8) * 2_8**(ilevel-1)

    mkey = grid_to_morton(igrid, ilevel)
    nkey = morton_neighbor(mkey, j, nmax_x, nmax_y, nmax_z)

    if (nkey == MK_MINUS1) then
       icell = 0
       return
    end if

    call morton_decode(nkey, nix, niy, niz)

    if (ilevel == 1) then
       ! Coarse cell (safe downcast: level 1 coords always small)
       icell = 1 + int(nix) + int(niy)*nx + int(niz)*nx*ny
    else
       ! Fine cell: find parent grid at level ilevel-1
       pix = nix / 2_8
       piy = niy / 2_8
       piz = niz / 2_8
       ind_oct = 1 + int(mod(nix, 2_8)) + 2*int(mod(niy, 2_8)) + 4*int(mod(niz, 2_8))
       pkey = morton_encode(pix, piy, piz)
       igrid_parent = morton_hash_lookup(mort_table(ilevel-1), pkey)
       if (igrid_parent > 0) then
          icell = ncoarse + (ind_oct-1)*ngridmax + igrid_parent
       else
          icell = 0
       end if
    end if
  end function morton_nbor_cell

end module morton_hash
