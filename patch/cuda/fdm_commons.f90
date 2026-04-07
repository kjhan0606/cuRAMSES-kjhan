!################################################################
! Fuzzy Dark Matter (FDM) commons module
! Stores parameters and shared state for the Schrödinger-Poisson solver
!################################################################
module fdm_commons
  use amr_commons, only: dp
  implicit none

  ! Computed at runtime from m_axion and code units
  real(dp)::hbar_code = 0.0d0   ! Effective hbar in code units

  ! Physical constants (CGS)
  real(dp),parameter::hbar_cgs = 1.0545718d-27   ! erg s
  real(dp),parameter::eV_to_erg = 1.602176634d-12 ! erg/eV
  real(dp),parameter::c_light = 2.99792458d10     ! cm/s

end module fdm_commons
