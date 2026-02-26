module neutrino_commons
  use amr_parameters, only: dp
  implicit none

  ! Transfer function table: R_nu(k,a) = T_nu/T_cb from CAMB
  integer :: nu_nk = 0, nu_na = 0
  real(dp), allocatable :: nu_log_k(:)    ! log10(k) in h/Mpc, size nu_nk
  real(dp), allocatable :: nu_log_a(:)    ! log10(a), size nu_na
  real(dp), allocatable :: nu_ratio(:,:)  ! R_nu(ik, ia), size (nu_nk, nu_na)
  logical :: nu_table_loaded = .false.

contains

  subroutine read_neutrino_table(filename)
    !-----------------------------------------------------------
    ! Read CAMB transfer function ratio table R_nu(k,a)
    ! Format: see misc/generate_neutrino_table.py
    !-----------------------------------------------------------
    implicit none
#ifndef WITHOUTMPI
    include 'mpif.h'
#endif
    character(len=*), intent(in) :: filename
    integer :: iu, ik, ia, ios
    character(len=512) :: line
    real(dp), allocatable :: k_arr(:), a_arr(:)
    integer :: info

    if(nu_table_loaded) return

    ! Only rank 0 reads, then broadcast
    iu = 79  ! free unit number

    ! Read on all ranks (file is accessible via shared filesystem)
    open(iu, file=trim(filename), status='old', action='read', iostat=ios)
    if(ios /= 0) then
       write(*,*) 'ERROR: Cannot open neutrino table: ', trim(filename)
       stop
    end if

    ! Skip comment lines starting with #
    do
       read(iu, '(A)', iostat=ios) line
       if(ios /= 0) then
          write(*,*) 'ERROR: Unexpected end of neutrino table'
          stop
       end if
       line = adjustl(line)
       if(line(1:1) /= '#') exit
    end do

    ! First non-comment line: nk na
    read(line, *) nu_nk, nu_na

    ! Allocate arrays
    allocate(k_arr(nu_nk))
    allocate(a_arr(nu_na))
    allocate(nu_log_k(nu_nk))
    allocate(nu_log_a(nu_na))
    allocate(nu_ratio(nu_nk, nu_na))

    ! Skip comment line for k values
    do
       read(iu, '(A)', iostat=ios) line
       if(ios /= 0) exit
       line = adjustl(line)
       if(line(1:1) /= '#') then
          backspace(iu)
          exit
       end if
    end do

    ! Read k values
    read(iu, *, iostat=ios) (k_arr(ik), ik=1,nu_nk)
    if(ios /= 0) then
       write(*,*) 'ERROR: Failed reading k values from neutrino table'
       stop
    end if

    ! Skip comment line for a values
    do
       read(iu, '(A)', iostat=ios) line
       if(ios /= 0) exit
       line = adjustl(line)
       if(line(1:1) /= '#') then
          backspace(iu)
          exit
       end if
    end do

    ! Read a values
    read(iu, *, iostat=ios) (a_arr(ia), ia=1,nu_na)
    if(ios /= 0) then
       write(*,*) 'ERROR: Failed reading a values from neutrino table'
       stop
    end if

    ! Skip comment line for R_nu matrix
    do
       read(iu, '(A)', iostat=ios) line
       if(ios /= 0) exit
       line = adjustl(line)
       if(line(1:1) /= '#') then
          backspace(iu)
          exit
       end if
    end do

    ! Read R_nu matrix: nk rows x na columns
    do ik = 1, nu_nk
       read(iu, *, iostat=ios) (nu_ratio(ik, ia), ia=1,nu_na)
       if(ios /= 0) then
          write(*,*) 'ERROR: Failed reading R_nu row', ik
          stop
       end if
    end do

    close(iu)

    ! Convert to log10
    do ik = 1, nu_nk
       nu_log_k(ik) = log10(k_arr(ik))
    end do
    do ia = 1, nu_na
       nu_log_a(ia) = log10(a_arr(ia))
    end do

    deallocate(k_arr, a_arr)
    nu_table_loaded = .true.

  end subroutine read_neutrino_table


  function get_nu_ratio(k_hmpc, aexp_val) result(R_nu)
    !-----------------------------------------------------------
    ! Interpolate R_nu at given k (h/Mpc) and a (scale factor)
    ! Bilinear interpolation in (log10 k, log10 a)
    ! Clamp to table boundaries
    !-----------------------------------------------------------
    implicit none
    real(dp), intent(in) :: k_hmpc, aexp_val
    real(dp) :: R_nu
    real(dp) :: logk, loga, tk, ta
    integer :: ik, ia

    if(.not. nu_table_loaded .or. nu_nk < 2 .or. nu_na < 2) then
       R_nu = 0.0d0
       return
    end if

    logk = log10(max(k_hmpc, 1.0d-30))
    loga = log10(max(aexp_val, 1.0d-30))

    ! Clamp to table range
    if(logk <= nu_log_k(1)) then
       ik = 1; tk = 0.0d0
    else if(logk >= nu_log_k(nu_nk)) then
       ik = nu_nk - 1; tk = 1.0d0
    else
       ! Binary search for k index
       call bisect_search(nu_log_k, nu_nk, logk, ik)
       tk = (logk - nu_log_k(ik)) / (nu_log_k(ik+1) - nu_log_k(ik))
    end if

    if(loga <= nu_log_a(1)) then
       ia = 1; ta = 0.0d0
    else if(loga >= nu_log_a(nu_na)) then
       ia = nu_na - 1; ta = 1.0d0
    else
       call bisect_search(nu_log_a, nu_na, loga, ia)
       ta = (loga - nu_log_a(ia)) / (nu_log_a(ia+1) - nu_log_a(ia))
    end if

    ! Bilinear interpolation
    R_nu = (1.0d0-tk)*(1.0d0-ta) * nu_ratio(ik,   ia)   &
         + tk        *(1.0d0-ta) * nu_ratio(ik+1, ia)   &
         + (1.0d0-tk)*ta         * nu_ratio(ik,   ia+1) &
         + tk        *ta         * nu_ratio(ik+1, ia+1)

  end function get_nu_ratio


  subroutine bisect_search(arr, n, val, idx)
    !-----------------------------------------------------------
    ! Binary search: find idx such that arr(idx) <= val < arr(idx+1)
    ! arr must be sorted in ascending order, 1 <= idx <= n-1
    !-----------------------------------------------------------
    implicit none
    integer, intent(in) :: n
    real(dp), intent(in) :: arr(n), val
    integer, intent(out) :: idx
    integer :: lo, hi, mid

    lo = 1
    hi = n
    do while(hi - lo > 1)
       mid = (lo + hi) / 2
       if(arr(mid) <= val) then
          lo = mid
       else
          hi = mid
       end if
    end do
    idx = lo

  end subroutine bisect_search

end module neutrino_commons
