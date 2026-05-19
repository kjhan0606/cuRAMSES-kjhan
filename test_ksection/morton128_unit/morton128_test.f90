program morton128_test
  use morton_keys
  use morton_hash
  implicit none

  integer :: nfail, ntest
  integer(8) :: ix, iy, iz, jx, jy, jz, nmax
  type(mkey_t) :: k, k2, kp, kc, kn
  type(morton_hash_table) :: tab
  integer :: igrid, ilevel, ic
  integer(8) :: c

  nfail = 0
  ntest = 0

  write(*,*) '======================================================'
  write(*,*) ' MORTON128 dual-int8 unit test'
  write(*,*) ' MORTON_MAXBITS = ', MORTON_MAXBITS
  write(*,*) '======================================================'

  ! ---- Test 1: encode/decode round trip at various levels ----
  call check_roundtrip(0_8, 0_8, 0_8, 'origin')
  call check_roundtrip(1_8, 1_8, 1_8, 'coord 1 (lvl 2)')
  call check_roundtrip(1023_8, 1023_8, 1023_8, '2^10-1 (lvl 11)')
  call check_roundtrip(2097151_8, 2097151_8, 2097151_8, '2^21-1 (lvl 22, max lo)')
  call check_roundtrip(2097152_8, 0_8, 0_8, '2^21 in x ONLY (first hi-bit)')
  call check_roundtrip(0_8, 2097152_8, 0_8, '2^21 in y ONLY')
  call check_roundtrip(0_8, 0_8, 2097152_8, '2^21 in z ONLY')
  call check_roundtrip(4194303_8, 4194303_8, 4194303_8, '2^22-1 (lvl 23 deep)')
  call check_roundtrip(4194304_8, 4194304_8, 4194304_8, '2^22 (lvl 24)')
  call check_roundtrip(2147483647_8, 2147483647_8, 2147483647_8, '2^31-1 (lvl 32)')
  if (MORTON_MAXBITS >= 42) then
    call check_roundtrip(4398046511103_8, 4398046511103_8, 4398046511103_8, '2^42-1 (lvl 43 deepest)')
  end if

  ! ---- Test 2: bit position 66 (3*22) lands in hi field ----
  k = morton_encode(4194304_8, 0_8, 0_8)  ! ix=2^22 → bit 22 of ix → key bit pos 66
  ntest = ntest + 1
  if (k%lo /= 0_8) then
    write(*,*) 'FAIL bit-66 hi check: lo should be 0 but got ', k%lo
    nfail = nfail + 1
#ifdef MORTON128
  else if (k%hi == 0_8) then
    write(*,*) 'FAIL bit-66 hi check: hi should be nonzero but got 0'
    nfail = nfail + 1
  else
    write(*,*) 'PASS  bit-66 lands in hi: lo=', k%lo, ' hi=', k%hi
#else
  else
    write(*,*) 'PASS  bit-66 in 64-bit (lo): lo=', k%lo
#endif
  end if

  ! ---- Test 3: parent/child round trip at deep level ----
  call check_parent_child(2097151_8, 1234567_8, 999999_8, 'lvl 22 boundary')
  call check_parent_child(4194303_8, 4194300_8, 4194301_8, 'lvl 23')
  if (MORTON_MAXBITS >= 42) then
    call check_parent_child(4398046511100_8, 4398046511101_8, 4398046511102_8, 'lvl 43 deep')
  end if

  ! ---- Test 4: child(parent(key)) recovers key ----
  do ilevel = 10, MORTON_MAXBITS, 8
     ix = ishft(1_8, ilevel-1) - 7_8
     iy = ishft(1_8, ilevel-1) - 3_8
     iz = ishft(1_8, ilevel-1) - 5_8
     if (ilevel == MORTON_MAXBITS) then
        ix = ishft(1_8, ilevel-1) - 7_8
     end if
     if (ix < 0_8 .or. iy < 0_8 .or. iz < 0_8) cycle
     k = morton_encode(ix, iy, iz)
     kp = morton_parent(k)
     ! ix,iy,iz lowest 3 bits → child index
     ic = 0
     if (btest(ix, 0)) ic = ic + 1
     if (btest(iy, 0)) ic = ic + 2
     if (btest(iz, 0)) ic = ic + 4
     kc = morton_child(kp, ic)
     ntest = ntest + 1
     if (kc == k) then
        write(*,'(A,I0,A,3I20)') ' PASS  parent/child @ lvl ', ilevel, '  coord=', ix, iy, iz
     else
        write(*,'(A,I0)') ' FAIL  parent/child @ lvl ', ilevel
        nfail = nfail + 1
     end if
  end do

  ! ---- Test 5: neighbor walk at deep coords ----
  call check_neighbor(4194304_8, 4194304_8, 4194304_8, 8388608_8, 'deep +/- neighbors')
  if (MORTON_MAXBITS >= 42) then
    call check_neighbor(1099511627776_8, 1099511627776_8, 1099511627776_8, 2199023255552_8, 'lvl 41 neighbors')
  end if

  ! ---- Test 6: hash insert/lookup at deep keys ----
  call morton_hash_init(tab, 32)
  ntest = ntest + 1
  do ic = 0, 99
    ix = 4194304_8 + int(ic, 8) * 17_8
    iy = 4194304_8 + int(ic, 8) * 31_8
    iz = 4194304_8 + int(ic, 8) * 7_8
    k = morton_encode(ix, iy, iz)
    call morton_hash_insert(tab, k, ic + 1000)
  end do
  ! Verify all lookups
  do ic = 0, 99
    ix = 4194304_8 + int(ic, 8) * 17_8
    iy = 4194304_8 + int(ic, 8) * 31_8
    iz = 4194304_8 + int(ic, 8) * 7_8
    k = morton_encode(ix, iy, iz)
    igrid = morton_hash_lookup(tab, k)
    if (igrid /= ic + 1000) then
       write(*,*) 'FAIL hash lookup @ deep coord ic=', ic, ' got=', igrid
       nfail = nfail + 1
       exit
    end if
  end do
  if (nfail == 0 .or. tab%count == 100) then
    write(*,*) 'PASS  hash insert/lookup 100 deep keys (cap=', tab%capacity, ')'
  end if
  call morton_hash_destroy(tab)

  ! ---- Test 7: bitshift across the 64-bit boundary ----
  k = morton_encode(2097152_8, 2097152_8, 2097152_8)  ! all coords = 2^21 (bit 21)
  ! Each axis's bit 21 sits at key bit positions 3*21=63, 3*21+1=64, 3*21+2=65
  ! So lo should have bit 63 set, hi should have bits 0 and 1 set
  ntest = ntest + 1
#ifdef MORTON128
  if (k%lo == ishft(1_8, 63) .and. k%hi == 3_8) then
    write(*,*) 'PASS  bit 63/64/65 boundary: lo=2^63 hi=3'
  else
    write(*,'(A,I20,A,I20)') ' FAIL  bit boundary: lo=', k%lo, ' hi=', k%hi
    nfail = nfail + 1
  end if
#else
  if (k%lo /= 0_8) then
    write(*,*) 'PASS  64-bit mode boundary: lo=', k%lo
  else
    write(*,*) 'FAIL  64-bit mode: lo=0 (lost bits)'
    nfail = nfail + 1
  end if
#endif

  ! ---- Test 8: ishft round trip across boundary ----
  k = morton_encode(2097151_8, 2097151_8, 2097151_8)  ! max lo coord (bits 0..20)
  k2 = mk_ishft(k, 3)   ! shift left by 3 — should partly cross into hi
  ix = -1_8; iy = -1_8; iz = -1_8
  call morton_decode(k, ix, iy, iz)
  call morton_decode(mk_ishft(k2, -3), jx, jy, jz)
  ntest = ntest + 1
  if (ix == jx .and. iy == jy .and. iz == jz) then
    write(*,*) 'PASS  ishft(+3) then ishft(-3) preserves key'
  else
    write(*,*) 'FAIL  ishft round trip: orig=', ix, iy, iz, ' got=', jx, jy, jz
    nfail = nfail + 1
  end if

  ! ---- Summary ----
  write(*,*) '======================================================'
  write(*,'(A,I0,A,I0,A)') '  Total: ', ntest, ' tests, ', nfail, ' failures'
  write(*,*) '======================================================'

  if (nfail > 0) then
     stop 1
  end if

contains

  subroutine check_roundtrip(ix, iy, iz, label)
    integer(8), intent(in) :: ix, iy, iz
    character(len=*), intent(in) :: label
    integer(8) :: jx, jy, jz
    type(mkey_t) :: k

    ntest = ntest + 1
    k = morton_encode(ix, iy, iz)
    call morton_decode(k, jx, jy, jz)
    if (ix == jx .and. iy == jy .and. iz == jz) then
#ifdef MORTON128
       write(*,'(A,A,A,3I16,A,Z16,Z16)') ' PASS  ', trim(label), '  in=', ix,iy,iz, ' key=', k%lo, k%hi
#else
       write(*,'(A,A,A,3I16,A,Z16)') ' PASS  ', trim(label), '  in=', ix,iy,iz, ' key=', k%lo
#endif
    else
       write(*,'(A,A,A,3I16,A,3I16)') ' FAIL  ', trim(label), ' in=', ix,iy,iz, ' decoded=', jx,jy,jz
       nfail = nfail + 1
    end if
  end subroutine check_roundtrip

  subroutine check_parent_child(ix, iy, iz, label)
    integer(8), intent(in) :: ix, iy, iz
    character(len=*), intent(in) :: label
    type(mkey_t) :: k, kp, kc
    integer :: ic
    integer(8) :: pix, piy, piz, cix, ciy, ciz

    ntest = ntest + 1
    k = morton_encode(ix, iy, iz)
    kp = morton_parent(k)
    call morton_decode(kp, pix, piy, piz)
    ! Parent coord should be (ix>>1, iy>>1, iz>>1)
    if (pix /= ishft(ix,-1) .or. piy /= ishft(iy,-1) .or. piz /= ishft(iz,-1)) then
       write(*,*) 'FAIL parent ', trim(label), '  got=', pix,piy,piz, ' expected=', ishft(ix,-1),ishft(iy,-1),ishft(iz,-1)
       nfail = nfail + 1
       return
    end if
    ! Reconstruct child index
    ic = 0
    if (btest(ix,0)) ic = ic + 1
    if (btest(iy,0)) ic = ic + 2
    if (btest(iz,0)) ic = ic + 4
    kc = morton_child(kp, ic)
    call morton_decode(kc, cix, ciy, ciz)
    if (cix == ix .and. ciy == iy .and. ciz == iz) then
       write(*,'(A,A)') ' PASS  parent/child ', trim(label)
    else
       write(*,*) 'FAIL parent/child ', trim(label), ' got=', cix,ciy,ciz
       nfail = nfail + 1
    end if
  end subroutine check_parent_child

  subroutine check_neighbor(ix, iy, iz, nmax, label)
    integer(8), intent(in) :: ix, iy, iz, nmax
    character(len=*), intent(in) :: label
    type(mkey_t) :: k, kn
    integer(8) :: nx, ny, nz
    integer :: dir
    integer(8), dimension(6,3) :: exp
    logical :: all_ok

    all_ok = .true.
    k = morton_encode(ix, iy, iz)
    exp(1,:) = (/ix-1, iy, iz/)
    exp(2,:) = (/ix+1, iy, iz/)
    exp(3,:) = (/ix, iy-1, iz/)
    exp(4,:) = (/ix, iy+1, iz/)
    exp(5,:) = (/ix, iy, iz-1/)
    exp(6,:) = (/ix, iy, iz+1/)
    do dir = 1, 6
      kn = morton_neighbor(k, dir, nmax, nmax, nmax)
      call morton_decode(kn, nx, ny, nz)
      if (nx /= exp(dir,1) .or. ny /= exp(dir,2) .or. nz /= exp(dir,3)) then
         write(*,*) 'FAIL neighbor dir=', dir, ' got=', nx,ny,nz, ' exp=', exp(dir,:)
         all_ok = .false.
      end if
    end do
    ntest = ntest + 1
    if (all_ok) then
       write(*,'(A,A,A,3I16)') ' PASS  neighbor walk ', trim(label), '  coord=', ix,iy,iz
    else
       nfail = nfail + 1
    end if
  end subroutine check_neighbor

end program morton128_test
