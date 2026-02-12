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
  ! Morton key layout (64-bit signed integer, using bits 0-62):
  !   bit 3*i+0 = x coordinate bit i
  !   bit 3*i+1 = y coordinate bit i
  !   bit 3*i+2 = z coordinate bit i
  ! for i = 0..20  (21 bits per coordinate, 63 bits total)
  !
  ! Maximum coordinate value: 2^21 - 1 = 2,097,151
  ! For nx=1: supports up to level 22
  ! For nx=4: supports up to level 20
  !--------------------------------------------------------------
  use amr_parameters, only: dp, ndim
  implicit none

  ! Always use 8-byte integers for Morton keys
  integer, parameter :: mkb = 8

  ! Maximum bits per coordinate (21 bits → 63 total, fits in i8b)
  integer, parameter :: MORTON_MAXBITS = 21

contains

  !--------------------------------------------------------------
  ! Interleave bits of (ix, iy, iz) into a 63-bit Morton key
  !--------------------------------------------------------------
  pure function morton_encode(ix, iy, iz) result(key)
    integer, intent(in) :: ix, iy, iz
    integer(mkb) :: key
    integer(mkb) :: x, y, z
    integer :: i

    key = 0_mkb
    x = int(ix, mkb)
    y = int(iy, mkb)
    z = int(iz, mkb)

    do i = 0, MORTON_MAXBITS - 1
       ! x bit i → position 3*i
       key = ior(key, ishft(iand(ishft(x, -i), 1_mkb), 3*i))
       ! y bit i → position 3*i+1
       key = ior(key, ishft(iand(ishft(y, -i), 1_mkb), 3*i+1))
       ! z bit i → position 3*i+2
       key = ior(key, ishft(iand(ishft(z, -i), 1_mkb), 3*i+2))
    end do
  end function morton_encode

  !--------------------------------------------------------------
  ! Decode Morton key back to (ix, iy, iz)
  !--------------------------------------------------------------
  pure subroutine morton_decode(key, ix, iy, iz)
    integer(mkb), intent(in) :: key
    integer, intent(out) :: ix, iy, iz
    integer(mkb) :: x, y, z
    integer :: i

    x = 0_mkb; y = 0_mkb; z = 0_mkb
    do i = 0, MORTON_MAXBITS - 1
       x = ior(x, ishft(iand(ishft(key, -(3*i  )), 1_mkb), i))
       y = ior(y, ishft(iand(ishft(key, -(3*i+1)), 1_mkb), i))
       z = ior(z, ishft(iand(ishft(key, -(3*i+2)), 1_mkb), i))
    end do
    ix = int(x)
    iy = int(y)
    iz = int(z)
  end subroutine morton_decode

  !--------------------------------------------------------------
  ! Compute Morton key of neighbor grid in given direction
  ! RAMSES convention (matching nbor array):
  !   dir: 1=-x, 2=+x, 3=-y, 4=+y, 5=-z, 6=+z
  ! nmax_x/y/z: number of grids per dimension at this level
  ! Returns -1 if neighbor is out of bounds
  !--------------------------------------------------------------
  pure function morton_neighbor(key, dir, nmax_x, nmax_y, nmax_z) result(nkey)
    integer(mkb), intent(in) :: key
    integer, intent(in) :: dir, nmax_x, nmax_y, nmax_z
    integer(mkb) :: nkey
    integer :: ix, iy, iz

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
       nkey = -1_mkb
       return
    end if

    nkey = morton_encode(ix, iy, iz)
  end function morton_neighbor

  !--------------------------------------------------------------
  ! Parent Morton key (remove lowest 3 bits = one level up)
  !--------------------------------------------------------------
  pure function morton_parent(key) result(pkey)
    integer(mkb), intent(in) :: key
    integer(mkb) :: pkey
    pkey = ishft(key, -3)
  end function morton_parent

  !--------------------------------------------------------------
  ! Child Morton key (add 3 bits for subcell)
  ! ichild: 0-7  (ix + 2*iy + 4*iz within parent oct)
  !--------------------------------------------------------------
  pure function morton_child(key, ichild) result(ckey)
    integer(mkb), intent(in) :: key
    integer, intent(in) :: ichild
    integer(mkb) :: ckey
    ckey = ior(ishft(key, 3), int(ichild, mkb))
  end function morton_child

  !--------------------------------------------------------------
  ! Compute Morton key from grid's xg coordinates and level
  ! At level l, grid centers are at (m + 0.5) * dx
  ! where dx = 0.5^(l-1) and m is the integer coordinate.
  ! So m = floor(xg * 2^(l-1))
  !--------------------------------------------------------------
  function grid_to_morton(igrid, ilevel) result(key)
    use amr_commons, only: xg
    integer, intent(in) :: igrid, ilevel
    integer(mkb) :: key
    integer :: ix, iy, iz
    real(dp) :: twotol

    twotol = 2.0d0**(ilevel-1)
    ix = int(xg(igrid, 1) * twotol)
    iy = int(xg(igrid, 2) * twotol)
    iz = int(xg(igrid, 3) * twotol)

    key = morton_encode(ix, iy, iz)
  end function grid_to_morton

end module morton_keys
