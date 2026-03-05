module dark_energy_commons
  !--------------------------------------------------------------
  ! Dark energy perturbation support for CPL w(a) = w0 + wa*(1-a)
  !
  ! Provides k-space correction factor for FFT Poisson solver:
  !   de_factor(k) = (k_tilde2 + kappa2) / (k_tilde2 + kappa2 + alpha)
  !
  ! where k_tilde2 = (2*pi)^2 * (kx^2 + ky^2 + kz^2)
  !       kappa2 = a^2 * E^2(a) * (boxlen_ini/C_H100)^2 / cs2_de
  !       alpha  = 1.5 * (boxlen_ini/C_H100)^2 * omega_l * f_de(a) * (1+w(a)) / cs2_de
  !
  ! Physics: combined Poisson + quasi-static Helmholtz system
  !   Poisson:   nabla^2 Phi = fourpi_m * delta_m + fourpi_de * delta_de
  !   Helmholtz: (nabla^2 - kappa^2) delta_de = [(1+w)/cs2] * nabla^2 Phi
  !
  ! For LCDM (w0=-1, wa=0): (1+w)=0, alpha=0, de_factor=1 (bit-wise)
  !--------------------------------------------------------------
  use amr_parameters, only: dp
  implicit none

  ! c / (H0/100) in Mpc  (= 299792.458 / 100)
  real(dp), parameter :: C_H100 = 2997.92458d0

  ! CAMB transfer function table: R_DE(k,a) = delta_DE/delta_m
  integer :: de_nk = 0, de_na = 0
  real(dp), allocatable :: de_log_k(:)    ! log10(k) in h/Mpc, size de_nk
  real(dp), allocatable :: de_log_a(:)    ! log10(a), size de_na
  real(dp), allocatable :: de_ratio(:,:)  ! R_DE(ik, ia), size (de_nk, de_na)
  logical :: de_table_loaded = .false.

contains

  !--------------------------------------------------------------
  ! f_de: DE density ratio rho_de(a)/rho_de(1)
  ! f_de = a^{-3(1+w0+wa)} * exp(-3*wa*(1-a))
  ! Matches f_de() in init_time.f90
  !--------------------------------------------------------------
  function f_de_val(a) result(fde)
    use amr_parameters, only: w0, wa
    real(dp), intent(in) :: a
    real(dp) :: fde
    if (wa == 0.0d0 .and. w0 == -1.0d0) then
       fde = 1.0d0
    else if (wa == 0.0d0) then
       fde = a**(-3.0d0*(1.0d0 + w0))
    else
       fde = a**(-3.0d0*(1.0d0 + w0 + wa)) * exp(-3.0d0*wa*(1.0d0 - a))
    end if
  end function f_de_val

  !--------------------------------------------------------------
  ! compute_de_kspace_params: precompute kappa2 and alpha
  ! for DE k-space correction in FFT Poisson solver
  !
  ! kappa2 = a^2 * E^2(a) * (boxlen_ini/C_H100)^2 / cs2_de
  !   where E^2(a) = omega_m/a^3 + omega_l*f_de(a) + omega_k/a^2
  !
  ! alpha = 1.5 * (boxlen_ini/C_H100)^2 * omega_l * f_de(a) * (1+w(a)) / cs2_de
  !--------------------------------------------------------------
  subroutine compute_de_kspace_params(aexp_val, kappa2_out, alpha_out)
    use amr_parameters, only: omega_m, omega_l, omega_k, w0, wa, &
         cs2_de, boxlen_ini
    real(dp), intent(in)  :: aexp_val
    real(dp), intent(out) :: kappa2_out, alpha_out
    real(dp) :: fde, E2_a, w_a, boxratio_sq

    fde = f_de_val(aexp_val)

    ! Friedmann equation: E^2(a) = H^2(a)/H0^2
    E2_a = omega_m / aexp_val**3 + omega_l * fde + omega_k / aexp_val**2

    ! Current DE equation of state
    w_a = w0 + wa * (1.0d0 - aexp_val)

    ! (boxlen_ini / (c/H0))^2  [dimensionless]
    boxratio_sq = (boxlen_ini / C_H100)**2

    ! Helmholtz mass term: kappa^2 = (aH/c)^2 / cs2_de  [in box^-2 units]
    kappa2_out = aexp_val**2 * E2_a * boxratio_sq / cs2_de

    ! DE coupling term: alpha = fourpi_de_phys * (1+w) / cs2_de  [in box^-2 units]
    ! fourpi_de_phys = (3/2) * H0^2 * omega_l * f_de(a)
    alpha_out = 1.5d0 * boxratio_sq * omega_l * fde * (1.0d0 + w_a) / cs2_de

  end subroutine compute_de_kspace_params


  subroutine read_de_table(filename)
    !-----------------------------------------------------------
    ! Read CAMB DE transfer function ratio table R_DE(k,a)
    ! Format: same as neutrino_table.dat (see generate_de_table.py)
    !-----------------------------------------------------------
    implicit none
#ifndef WITHOUTMPI
    include 'mpif.h'
#endif
    character(len=*), intent(in) :: filename
    integer :: iu, ik, ia, ios
    character(len=512) :: line
    real(dp), allocatable :: k_arr(:), a_arr(:)

    if(de_table_loaded) return

    iu = 78  ! free unit number

    open(iu, file=trim(filename), status='old', action='read', iostat=ios)
    if(ios /= 0) then
       write(*,*) 'ERROR: Cannot open DE table: ', trim(filename)
       stop
    end if

    ! Skip comment lines starting with #
    do
       read(iu, '(A)', iostat=ios) line
       if(ios /= 0) then
          write(*,*) 'ERROR: Unexpected end of DE table'
          stop
       end if
       line = adjustl(line)
       if(line(1:1) /= '#') exit
    end do

    ! First non-comment line: nk na
    read(line, *) de_nk, de_na

    ! Allocate arrays
    allocate(k_arr(de_nk))
    allocate(a_arr(de_na))
    allocate(de_log_k(de_nk))
    allocate(de_log_a(de_na))
    allocate(de_ratio(de_nk, de_na))

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
    read(iu, *, iostat=ios) (k_arr(ik), ik=1,de_nk)
    if(ios /= 0) then
       write(*,*) 'ERROR: Failed reading k values from DE table'
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
    read(iu, *, iostat=ios) (a_arr(ia), ia=1,de_na)
    if(ios /= 0) then
       write(*,*) 'ERROR: Failed reading a values from DE table'
       stop
    end if

    ! Skip comment line for R_DE matrix
    do
       read(iu, '(A)', iostat=ios) line
       if(ios /= 0) exit
       line = adjustl(line)
       if(line(1:1) /= '#') then
          backspace(iu)
          exit
       end if
    end do

    ! Read R_DE matrix: nk rows x na columns
    do ik = 1, de_nk
       read(iu, *, iostat=ios) (de_ratio(ik, ia), ia=1,de_na)
       if(ios /= 0) then
          write(*,*) 'ERROR: Failed reading R_DE row', ik
          stop
       end if
    end do

    close(iu)

    ! Convert to log10
    do ik = 1, de_nk
       de_log_k(ik) = log10(k_arr(ik))
    end do
    do ia = 1, de_na
       de_log_a(ia) = log10(a_arr(ia))
    end do

    deallocate(k_arr, a_arr)
    de_table_loaded = .true.

  end subroutine read_de_table


  function get_de_ratio(k_hmpc, aexp_val) result(R_DE)
    !-----------------------------------------------------------
    ! Interpolate R_DE at given k (h/Mpc) and a (scale factor)
    ! Bilinear interpolation in (log10 k, log10 a)
    ! Clamp to table boundaries
    !-----------------------------------------------------------
    implicit none
    real(dp), intent(in) :: k_hmpc, aexp_val
    real(dp) :: R_DE
    real(dp) :: logk, loga, tk, ta
    integer :: ik, ia

    if(.not. de_table_loaded .or. de_nk < 2 .or. de_na < 2) then
       R_DE = 0.0d0
       return
    end if

    logk = log10(max(k_hmpc, 1.0d-30))
    loga = log10(max(aexp_val, 1.0d-30))

    ! Clamp to table range
    if(logk <= de_log_k(1)) then
       ik = 1; tk = 0.0d0
    else if(logk >= de_log_k(de_nk)) then
       ik = de_nk - 1; tk = 1.0d0
    else
       call bisect_search_de(de_log_k, de_nk, logk, ik)
       tk = (logk - de_log_k(ik)) / (de_log_k(ik+1) - de_log_k(ik))
    end if

    if(loga <= de_log_a(1)) then
       ia = 1; ta = 0.0d0
    else if(loga >= de_log_a(de_na)) then
       ia = de_na - 1; ta = 1.0d0
    else
       call bisect_search_de(de_log_a, de_na, loga, ia)
       ta = (loga - de_log_a(ia)) / (de_log_a(ia+1) - de_log_a(ia))
    end if

    ! Bilinear interpolation
    R_DE = (1.0d0-tk)*(1.0d0-ta) * de_ratio(ik,   ia)   &
         + tk        *(1.0d0-ta) * de_ratio(ik+1, ia)   &
         + (1.0d0-tk)*ta         * de_ratio(ik,   ia+1) &
         + tk        *ta         * de_ratio(ik+1, ia+1)

  end function get_de_ratio


  subroutine bisect_search_de(arr, n, val, idx)
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

  end subroutine bisect_search_de

end module dark_energy_commons
