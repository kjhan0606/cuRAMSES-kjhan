!################################################################
!################################################################
!################################################################
!################################################################
module morton_keys
  !--------------------------------------------------------------
  ! Morton (Z-order) key encoding for 3D AMR octree navigation.
  ! Provides O(1) bit-operation neighbor finding to replace
  ! the nbor() array indirection.
  !
  ! Morton keys use type(mkey_t) — a dual-64-bit struct:
  !   64-bit mode:  lo only (21 bits/coord, level 22)
  !   128-bit mode: lo + hi (42 bits/coord, level 43)
  !
  ! Coordinates (ix, iy, iz, nmax) are always integer(8).
  !
  ! Compile with -DMORTON128 for 128-bit keys.
  ! Both modes work with ifx and gfortran.
  !--------------------------------------------------------------
  use amr_parameters, only: dp, ndim
  implicit none

  ! Morton key width parameters
#ifdef MORTON128
  integer, parameter :: MORTON_MAXBITS = 42 ! 42 bits/coord, level 43 (nx=1)
#else
  integer, parameter :: MORTON_MAXBITS = 21 ! 21 bits/coord, level 22 (nx=1)
#endif

  ! Dual-64-bit Morton key type
  type :: mkey_t
     integer(8) :: lo        ! bits 0-63
#ifdef MORTON128
     integer(8) :: hi        ! bits 64-125
#endif
  end type mkey_t

  ! Constants
#ifdef MORTON128
  type(mkey_t), parameter :: MK_ZERO    = mkey_t(0_8, 0_8)
  type(mkey_t), parameter :: MK_MINUS1  = mkey_t(-1_8, -1_8)
  type(mkey_t), parameter :: MK_MINUS2  = mkey_t(-2_8, -1_8)
#else
  type(mkey_t), parameter :: MK_ZERO    = mkey_t(0_8)
  type(mkey_t), parameter :: MK_MINUS1  = mkey_t(-1_8)
  type(mkey_t), parameter :: MK_MINUS2  = mkey_t(-2_8)
#endif

  ! Operator overloads
  interface operator(==)
     module procedure mkey_eq
  end interface
  interface operator(/=)
     module procedure mkey_ne
  end interface

contains

  !--------------------------------------------------------------
  ! Operator: equality
  !--------------------------------------------------------------
  pure function mkey_eq(a, b) result(r)
    type(mkey_t), intent(in) :: a, b
    logical :: r
#ifdef MORTON128
    r = (a%lo == b%lo) .and. (a%hi == b%hi)
#else
    r = (a%lo == b%lo)
#endif
  end function mkey_eq

  !--------------------------------------------------------------
  ! Operator: inequality
  !--------------------------------------------------------------
  pure function mkey_ne(a, b) result(r)
    type(mkey_t), intent(in) :: a, b
    logical :: r
    r = .not. mkey_eq(a, b)
  end function mkey_ne

  !--------------------------------------------------------------
  ! Set bit at position pos in Morton key
  !--------------------------------------------------------------
  pure function mk_setbit(key, pos) result(r)
    type(mkey_t), intent(in) :: key
    integer, intent(in) :: pos
    type(mkey_t) :: r
    r = key
#ifdef MORTON128
    if (pos < 64) then
       r%lo = ibset(r%lo, pos)
    else
       r%hi = ibset(r%hi, pos - 64)
    end if
#else
    r%lo = ibset(r%lo, pos)
#endif
  end function mk_setbit

  !--------------------------------------------------------------
  ! Test bit at position pos in Morton key
  !--------------------------------------------------------------
  pure function mk_btest(key, pos) result(r)
    type(mkey_t), intent(in) :: key
    integer, intent(in) :: pos
    logical :: r
#ifdef MORTON128
    if (pos < 64) then
       r = btest(key%lo, pos)
    else
       r = btest(key%hi, pos - 64)
    end if
#else
    r = btest(key%lo, pos)
#endif
  end function mk_btest

  !--------------------------------------------------------------
  ! Shift Morton key by shift bits (positive=left, negative=right)
  !--------------------------------------------------------------
  pure function mk_ishft(key, shift) result(r)
    type(mkey_t), intent(in) :: key
    integer, intent(in) :: shift
    type(mkey_t) :: r
#ifdef MORTON128
    integer :: n
    if (shift == 0) then
       r = key
    else if (shift > 0) then
       ! Left shift
       n = shift
       if (n >= 64) then
          r%lo = 0_8
          r%hi = ishft(key%lo, n - 64)
       else
          r%lo = ishft(key%lo, n)
          r%hi = ior(ishft(key%hi, n), ishft(key%lo, n - 64))
       end if
    else
       ! Right shift
       n = -shift
       if (n >= 64) then
          r%lo = ishft(key%hi, -(n - 64))
          r%hi = 0_8
       else
          r%lo = ior(ishft(key%lo, -n), ishft(key%hi, 64 - n))
          r%hi = ishft(key%hi, -n)
       end if
    end if
#else
    r%lo = ishft(key%lo, shift)
#endif
  end function mk_ishft

  !--------------------------------------------------------------
  ! Bitwise OR of two Morton keys
  !--------------------------------------------------------------
  pure function mk_ior(a, b) result(r)
    type(mkey_t), intent(in) :: a, b
    type(mkey_t) :: r
    r%lo = ior(a%lo, b%lo)
#ifdef MORTON128
    r%hi = ior(a%hi, b%hi)
#endif
  end function mk_ior

  !--------------------------------------------------------------
  ! Bitwise OR of Morton key with small integer (affects lo only)
  !--------------------------------------------------------------
  pure function mk_ior_int(a, n) result(r)
    type(mkey_t), intent(in) :: a
    integer, intent(in) :: n
    type(mkey_t) :: r
    r%lo = ior(a%lo, int(n, 8))
#ifdef MORTON128
    r%hi = a%hi
#endif
  end function mk_ior_int

  !--------------------------------------------------------------
  ! Interleave bits of (ix, iy, iz) into a Morton key
  !--------------------------------------------------------------
  pure function morton_encode(ix, iy, iz) result(key)
    integer(8), intent(in) :: ix, iy, iz
    type(mkey_t) :: key
    integer :: i

    key = MK_ZERO
    do i = 0, MORTON_MAXBITS - 1
       if (btest(ix, i)) key = mk_setbit(key, 3*i)
       if (btest(iy, i)) key = mk_setbit(key, 3*i+1)
       if (btest(iz, i)) key = mk_setbit(key, 3*i+2)
    end do
  end function morton_encode

  !--------------------------------------------------------------
  ! Decode Morton key back to (ix, iy, iz)
  !--------------------------------------------------------------
  pure subroutine morton_decode(key, ix, iy, iz)
    type(mkey_t), intent(in) :: key
    integer(8), intent(out) :: ix, iy, iz
    integer :: i

    ix = 0_8; iy = 0_8; iz = 0_8
    do i = 0, MORTON_MAXBITS - 1
       if (mk_btest(key, 3*i))   ix = ibset(ix, i)
       if (mk_btest(key, 3*i+1)) iy = ibset(iy, i)
       if (mk_btest(key, 3*i+2)) iz = ibset(iz, i)
    end do
  end subroutine morton_decode

  !--------------------------------------------------------------
  ! Compute Morton key of neighbor grid in given direction
  ! RAMSES convention: 1=-x, 2=+x, 3=-y, 4=+y, 5=-z, 6=+z
  ! Returns MK_MINUS1 if neighbor is out of bounds
  !--------------------------------------------------------------
  pure function morton_neighbor(key, dir, nmax_x, nmax_y, nmax_z) result(nkey)
    type(mkey_t), intent(in) :: key
    integer, intent(in) :: dir
    integer(8), intent(in) :: nmax_x, nmax_y, nmax_z
    type(mkey_t) :: nkey
    integer(8) :: ix, iy, iz

    call morton_decode(key, ix, iy, iz)

    select case(dir)
    case(1); ix = ix - 1  ! -x
    case(2); ix = ix + 1  ! +x
    case(3); iy = iy - 1  ! -y
    case(4); iy = iy + 1  ! +y
    case(5); iz = iz - 1  ! -z
    case(6); iz = iz + 1  ! +z
    end select

    ! Periodic boundary: wrap around
    if (ix < 0) then
       ix = ix + nmax_x
    else if (ix >= nmax_x) then
       ix = ix - nmax_x
    end if
    if (iy < 0) then
       iy = iy + nmax_y
    else if (iy >= nmax_y) then
       iy = iy - nmax_y
    end if
    if (iz < 0) then
       iz = iz + nmax_z
    else if (iz >= nmax_z) then
       iz = iz - nmax_z
    end if

    ! Check bounds (for non-periodic: caller should not get here)
    if (ix < 0 .or. ix >= nmax_x .or. &
        iy < 0 .or. iy >= nmax_y .or. &
        iz < 0 .or. iz >= nmax_z) then
       nkey = MK_MINUS1
       return
    end if

    nkey = morton_encode(ix, iy, iz)
  end function morton_neighbor

  !--------------------------------------------------------------
  ! Parent Morton key (remove lowest 3 bits = one level up)
  !--------------------------------------------------------------
  pure function morton_parent(key) result(pkey)
    type(mkey_t), intent(in) :: key
    type(mkey_t) :: pkey
    pkey = mk_ishft(key, -3)
  end function morton_parent

  !--------------------------------------------------------------
  ! Child Morton key (add 3 bits for subcell)
  ! ichild: 0-7  (ix + 2*iy + 4*iz within parent oct)
  !--------------------------------------------------------------
  pure function morton_child(key, ichild) result(ckey)
    type(mkey_t), intent(in) :: key
    integer, intent(in) :: ichild
    type(mkey_t) :: ckey
    ckey = mk_ior_int(mk_ishft(key, 3), ichild)
  end function morton_child

  !--------------------------------------------------------------
  ! Compute Morton key from grid's xg coordinates and level
  !--------------------------------------------------------------
  function grid_to_morton(igrid, ilevel) result(key)
    use amr_commons, only: xg
    integer, intent(in) :: igrid, ilevel
    type(mkey_t) :: key
    integer(8) :: ix, iy, iz
    real(dp) :: twotol

    twotol = 2.0d0**(ilevel-1)
    ix = int(xg(igrid, 1) * twotol, 8)
    iy = int(xg(igrid, 2) * twotol, 8)
    iz = int(xg(igrid, 3) * twotol, 8)

    key = morton_encode(ix, iy, iz)
  end function grid_to_morton

end module morton_keys
