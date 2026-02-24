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
  ! Compile with -DMORTON128 for 128-bit keys (42 bits/coord,
  ! level 43, requires gfortran), default is 64-bit keys
  ! (21 bits/coord, level 22, works with ifx and gfortran).
  !--------------------------------------------------------------
  use amr_parameters, only: dp, ndim
  implicit none

#ifdef MORTON128
  integer, parameter :: mkb = 16            ! 128-bit Morton key
  integer, parameter :: MORTON_MAXBITS = 42 ! 42 bits/coord, level 43 (nx=1)
#else
  integer, parameter :: mkb = 8             ! 64-bit Morton key
  integer, parameter :: MORTON_MAXBITS = 21 ! 21 bits/coord, level 22 (nx=1)
#endif

contains

  !--------------------------------------------------------------
  ! Interleave bits of (ix, iy, iz) into a 126-bit Morton key
  !--------------------------------------------------------------
  pure function morton_encode(ix, iy, iz) result(key)
    integer(mkb), intent(in) :: ix, iy, iz
    integer(mkb) :: key
    integer :: i

    key = 0_mkb

    do i = 0, MORTON_MAXBITS - 1
       ! x bit i → position 3*i
       key = ior(key, ishft(iand(ishft(ix, -i), 1_mkb), 3*i))
       ! y bit i → position 3*i+1
       key = ior(key, ishft(iand(ishft(iy, -i), 1_mkb), 3*i+1))
       ! z bit i → position 3*i+2
       key = ior(key, ishft(iand(ishft(iz, -i), 1_mkb), 3*i+2))
    end do
  end function morton_encode

  !--------------------------------------------------------------
  ! Decode Morton key back to (ix, iy, iz)
  !--------------------------------------------------------------
  pure subroutine morton_decode(key, ix, iy, iz)
    integer(mkb), intent(in) :: key
    integer(mkb), intent(out) :: ix, iy, iz
    integer :: i

    ix = 0_mkb; iy = 0_mkb; iz = 0_mkb
    do i = 0, MORTON_MAXBITS - 1
       ix = ior(ix, ishft(iand(ishft(key, -(3*i  )), 1_mkb), i))
       iy = ior(iy, ishft(iand(ishft(key, -(3*i+1)), 1_mkb), i))
       iz = ior(iz, ishft(iand(ishft(key, -(3*i+2)), 1_mkb), i))
    end do
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
    integer, intent(in) :: dir
    integer(mkb), intent(in) :: nmax_x, nmax_y, nmax_z
    integer(mkb) :: nkey
    integer(mkb) :: ix, iy, iz

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
    integer(mkb) :: ix, iy, iz
    real(dp) :: twotol

    twotol = 2.0d0**(ilevel-1)
    ix = int(xg(igrid, 1) * twotol, mkb)
    iy = int(xg(igrid, 2) * twotol, mkb)
    iz = int(xg(igrid, 3) * twotol, mkb)

    key = morton_encode(ix, iy, iz)
  end function grid_to_morton

end module morton_keys
