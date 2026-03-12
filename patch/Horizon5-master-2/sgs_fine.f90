!###########################################################
!###########################################################
! SGS (Sub-Grid Scale) Turbulence Source Terms
!
! Source terms for the SGS turbulent kinetic energy equation:
!   d(rho*e_sgs)/dt = -P_turb*div(v) + Production - Dissipation
!
! where:
!   P_turb = (2/3)*rho*e_sgs        (turbulent pressure)
!   Production ~ C_prod * rho * e_sgs * sqrt(factG*rho)  (gravitational stirring)
!   Dissipation = C_diss * rho * (2*e_sgs/3)^(3/2) / dx  (cascade to thermal)
!
! The advection of rho*e_sgs and P_turb in momentum/energy fluxes
! are handled by the Riemann solver (see godunov_utils).
!
! Dissipation uses an analytical implicit solve (unconditionally stable):
!   sigma_new = sigma_old / (1 + C_diss*sigma_old*dt/(2*dx))
!
! References:
!   Schmidt & Federrath (2011), Scannapieco & Bruggen (2008)
!###########################################################
subroutine sgs_fine(ilevel)
  use amr_commons
  use hydro_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer::ilevel

  integer::igrid,ngrid,ncache

  if(.not. use_sgs) return
  if(isgs <= 0) return
  if(numbtot(1,ilevel)==0) return

  ncache=active(ilevel)%ngrid
!$omp parallel do private(igrid,ngrid) schedule(dynamic)
  do igrid=1,ncache,nvector
     ngrid=MIN(nvector,ncache-igrid+1)
     call sgs_source(ilevel,igrid,ngrid)
  end do

end subroutine sgs_fine
!###########################################################
!###########################################################
subroutine sgs_source(ilevel,istart,ngrid)
  use amr_commons
  use hydro_commons
  use cooling_module, ONLY: twopi
  implicit none

  integer::ilevel,istart,ngrid

  ! Local variables
  integer::i,ind,iskip,nleaf
  integer,dimension(1:nvector)::ind_grid,ind_cell,ind_leaf
  real(dp)::dx,dx_loc,scale,dt_loc
  real(dp)::pi,factG_loc
  real(dp)::rho,esgs,sigma,sigma_new,esgs_new
  real(dp)::diss_rate,prod_rate,pdv_rate
  real(dp)::esgs_dissipated,pturb
  integer::nx_loc

  ! Mesh spacing in that level
  dx=0.5D0**ilevel
  nx_loc=(icoarse_max-icoarse_min+1)
  scale=boxlen/dble(nx_loc)
  dx_loc=dx*scale

  ! Time step at this level
  dt_loc=dtnew(ilevel)

  ! Gravitational constant factor (same as jeans_length_refine)
  pi = twopi / 2d0
  factG_loc=1d0
  if(cosmo)factG_loc=3d0/8d0/pi*omega_m*aexp

  ! Gather grid indices
  do i=1,ngrid
     ind_grid(i)=active(ilevel)%igrid(istart+i-1)
  end do

  ! Loop over cells
  do ind=1,twotondim
     iskip=ncoarse+(ind-1)*ngridmax
     do i=1,ngrid
        ind_cell(i)=iskip+ind_grid(i)
     end do

     ! Gather leaf cells
     nleaf=0
     do i=1,ngrid
        if(son(ind_cell(i))==0)then
           nleaf=nleaf+1
           ind_leaf(nleaf)=ind_cell(i)
        end if
     end do
     if(nleaf==0)cycle

     ! Apply source terms to each leaf cell
     do i=1,nleaf
        rho  = max(uold(ind_leaf(i),1), smallr)
        esgs = uold(ind_leaf(i),isgs) / rho  ! specific SGS energy
        esgs = max(esgs, sgs_floor)

        ! Turbulent velocity dispersion: sigma^2 = 2*e_sgs/3
        sigma = sqrt(2d0*esgs/3d0)

        ! =============================================
        ! 1. PdV coupling: -P_turb * div(v)
        !    Uses divu array from godunov_fine if available
        !    Convention: divu(cell) = -div(v)*dt
        ! =============================================
        pdv_rate = 0d0
        if(pressure_fix) then
           ! divu = -div(v)*dt, so -P_turb*div(v)*dt = P_turb*divu
           pturb = (2d0/3d0)*rho*esgs
           pdv_rate = pturb * divu(ind_leaf(i)) / dt_loc  ! rate per unit time
        end if

        ! =============================================
        ! 2. Gravitational production
        !    P_grav = C_prod * rho * e_sgs * sqrt(factG*rho)
        !    Timescale ~ free-fall time t_ff = 1/sqrt(4*pi*G*rho)
        !    In code units: sqrt(factG*rho) is the free-fall rate
        ! =============================================
        prod_rate = sgs_C_prod * rho * esgs * sqrt(factG_loc * rho)

        ! =============================================
        ! 3. Dissipation: D = C_diss * rho * sigma^3 / dx
        !    Uses analytical implicit solution (unconditionally stable):
        !    sigma_new = sigma / (1 + C_diss*sigma*dt/(2*dx))
        ! =============================================
        ! First apply production + PdV explicitly
        esgs = esgs + (prod_rate/rho + pdv_rate/rho) * dt_loc
        esgs = max(esgs, sgs_floor)
        sigma = sqrt(2d0*esgs/3d0)

        ! Then apply dissipation analytically (implicit, stable)
        sigma_new = sigma / (1d0 + sgs_C_diss*sigma*dt_loc/(2d0*dx_loc))
        esgs_new = 1.5d0 * sigma_new**2

        ! Safety floor
        esgs_new = max(esgs_new, sgs_floor)

        ! Dissipated energy (to be added to thermal energy)
        esgs_dissipated = esgs - esgs_new
        esgs_dissipated = max(esgs_dissipated, 0d0)

        ! =============================================
        ! Update conserved variables
        ! =============================================
        ! Update SGS energy (rho*e_sgs)
        uold(ind_leaf(i),isgs) = rho * esgs_new

        ! Add dissipated SGS energy to thermal energy
        ! (total internal energy = uold(ndim+2) includes kinetic)
        ! The thermal part: E_total - E_kin - E_rad
        ! Adding dissipation: E_total += rho * esgs_dissipated
        uold(ind_leaf(i),ndim+2) = uold(ind_leaf(i),ndim+2) &
             & + rho * esgs_dissipated

     end do  ! nleaf

  end do  ! twotondim

end subroutine sgs_source
