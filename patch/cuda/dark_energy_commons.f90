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

end module dark_energy_commons
