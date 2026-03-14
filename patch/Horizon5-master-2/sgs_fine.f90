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
!   sigma_new = sigma_old / (1 + C_diss*sigma_old*dt/(3*dx))
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
  integer::i,j,ind,iskip
  integer,dimension(1:nvector)::ind_grid,ind_cell
  integer,dimension(1:nvector,0:twondim)::igridn
  integer,dimension(1:nvector,1:twondim)::indn
  real(dp)::dx,dx_loc,scale,dt_loc
  real(dp)::pi,factG_loc
  real(dp)::rho,esgs,sigma,sigma_new,esgs_new
  real(dp)::eint_cell,esgs_cap,esgs_dissipated_cap
  real(dp)::prod_rate,pdv_rate,shear_prod
  real(dp)::esgs_dissipated,pturb
  real(dp)::inv2dx
  real(dp)::gradPx,gradPy,gradPz  ! grad(P_turb) for source term
  real(dp)::dvx_dx,dvx_dy,dvx_dz
  real(dp)::dvy_dx,dvy_dy,dvy_dz
  real(dp)::dvz_dx,dvz_dy,dvz_dz
  real(dp)::S11,S22,S33,S12,S13,S23
  real(dp)::divv,Sstar2
  logical::ok_nbor
  integer::nx_loc

  ! Mesh spacing in that level
  dx=0.5D0**ilevel
  nx_loc=(icoarse_max-icoarse_min+1)
  scale=boxlen/dble(nx_loc)
  dx_loc=dx*scale
  inv2dx=0.5d0/dx_loc

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

  ! Get neighboring grids (once per batch)
  call getnborgrids(ind_grid,igridn,ngrid)

  ! Loop over cells
  do ind=1,twotondim
     iskip=ncoarse+(ind-1)*ngridmax
     do i=1,ngrid
        ind_cell(i)=iskip+ind_grid(i)
     end do

     ! Get face-neighbor cells for this oct position
     call getnborcells(igridn,ind,indn,ngrid)

     ! Apply source terms to each cell (skip non-leaf)
     do i=1,ngrid
        if(son(ind_cell(i))/=0) cycle

        rho  = max(uold(ind_cell(i),1), smallr)
        esgs = uold(ind_cell(i),isgs) / rho  ! specific SGS energy
        esgs = max(esgs, sgs_floor)

        ! Turbulent velocity dispersion: sigma^2 = 2*e_sgs/3
        sigma = sqrt(2d0*esgs/3d0)

        ! Check all 6 face neighbors exist (needed for PdV and shear)
        ok_nbor = .true.
        do j=1,twondim
           if(indn(i,j)==0) ok_nbor = .false.
        end do

        ! =============================================
        ! Compute velocity gradients (shared by PdV and shear)
        ! =============================================
        divv = 0d0
        Sstar2 = 0d0
        if(ok_nbor) then
              ! Velocity gradients from face neighbors
              ! j=1:-X, j=2:+X, j=3:-Y, j=4:+Y, j=5:-Z, j=6:+Z
              dvx_dx = (uold(indn(i,2),2)/uold(indn(i,2),1) &
                   & - uold(indn(i,1),2)/uold(indn(i,1),1)) * inv2dx
              dvx_dy = (uold(indn(i,4),2)/uold(indn(i,4),1) &
                   & - uold(indn(i,3),2)/uold(indn(i,3),1)) * inv2dx
              dvx_dz = (uold(indn(i,6),2)/uold(indn(i,6),1) &
                   & - uold(indn(i,5),2)/uold(indn(i,5),1)) * inv2dx
              dvy_dx = (uold(indn(i,2),3)/uold(indn(i,2),1) &
                   & - uold(indn(i,1),3)/uold(indn(i,1),1)) * inv2dx
              dvy_dy = (uold(indn(i,4),3)/uold(indn(i,4),1) &
                   & - uold(indn(i,3),3)/uold(indn(i,3),1)) * inv2dx
              dvy_dz = (uold(indn(i,6),3)/uold(indn(i,6),1) &
                   & - uold(indn(i,5),3)/uold(indn(i,5),1)) * inv2dx
              dvz_dx = (uold(indn(i,2),4)/uold(indn(i,2),1) &
                   & - uold(indn(i,1),4)/uold(indn(i,1),1)) * inv2dx
              dvz_dy = (uold(indn(i,4),4)/uold(indn(i,4),1) &
                   & - uold(indn(i,3),4)/uold(indn(i,3),1)) * inv2dx
              dvz_dz = (uold(indn(i,6),4)/uold(indn(i,6),1) &
                   & - uold(indn(i,5),4)/uold(indn(i,5),1)) * inv2dx

              ! Divergence of velocity
              divv = dvx_dx + dvy_dy + dvz_dz

              ! Strain rate tensor S_ij = (dv_i/dx_j + dv_j/dx_i)/2
              S11 = dvx_dx
              S22 = dvy_dy
              S33 = dvz_dz
              S12 = 0.5d0*(dvx_dy + dvy_dx)
              S13 = 0.5d0*(dvx_dz + dvz_dx)
              S23 = 0.5d0*(dvy_dz + dvz_dy)

              ! Traceless strain magnitude: |S*|^2 = S_ij*S_ij - (1/3)*(div v)^2
              Sstar2 = S11**2 + S22**2 + S33**2 &
                   & + 2d0*(S12**2 + S13**2 + S23**2) &
                   & - (1d0/3d0)*divv**2
        end if

        ! =============================================
        ! 1. PdV coupling: -P_turb * div(v)
        !    Only when sgs_hydro=T (P_turb is in Riemann flux).
        !    When sgs_hydro=F, P_turb is not in the flux, so PdV
        !    would create energy from nothing.
        ! =============================================
        pdv_rate = 0d0
        if(sgs_hydro .and. ok_nbor) then
           pturb = (2d0/3d0)*rho*esgs
           pdv_rate = -pturb * divv  ! -P_turb * div(v)
        end if

        ! =============================================
        ! 2. Gravitational production
        !    P_grav = C_prod * rho * e_sgs * sqrt(factG*rho)
        ! =============================================
        prod_rate = sgs_C_prod * rho * esgs * sqrt(factG_loc * rho)

        ! =============================================
        ! 2b. Smagorinsky shear production
        !     P_shear = 2*C_smag*rho*sigma*dx*|S*|^2
        ! =============================================
        shear_prod = 0d0
        if(sgs_C_smag > 0d0 .and. ok_nbor) then
              shear_prod = 2d0 * sgs_C_smag * rho * sigma * dx_loc * Sstar2
        end if

        prod_rate = prod_rate + shear_prod

        ! =============================================
        ! 3. Dissipation: D = C_diss * rho * sigma^3 / dx
        !    Implicit: sigma_new = sigma / (1 + C_diss*sigma*dt/(3*dx))
        ! =============================================
        ! First apply production + PdV explicitly
        esgs = esgs + (prod_rate/rho + pdv_rate/rho) * dt_loc
        esgs = max(esgs, sgs_floor)
        sigma = sqrt(2d0*esgs/3d0)

        ! Then apply dissipation analytically (implicit, stable)
        sigma_new = sigma / (1d0 + sgs_C_diss*sigma*dt_loc/(3d0*dx_loc))
        esgs_new = 1.5d0 * sigma_new**2

        ! =============================================
        ! Cap: P_turb <= sgs_cap * P_thermal
        !   (2/3)*rho*esgs <= sgs_cap * (gamma-1)*rho*e_int
        !   => esgs <= sgs_cap * (3/2)*(gamma-1)*e_int
        ! =============================================
        if(sgs_cap > 0d0) then
           eint_cell = uold(ind_cell(i),ndim+2)/rho &
                & - 0.5d0*(uold(ind_cell(i),2)**2+uold(ind_cell(i),3)**2 &
                &         +uold(ind_cell(i),4)**2)/rho**2
           eint_cell = max(eint_cell, smallc**2)
           esgs_cap = sgs_cap * 1.5d0 * (gamma-1d0) * eint_cell
           if(esgs_new > esgs_cap) then
              esgs_new = esgs_cap
              ! NOTE: Do NOT transfer excess to thermal.
              ! For gamma=5/3, gamma*(gamma-1)=10/9 > 2/3, so transferring
              ! SGS to thermal INCREASES c_eff, making the cap counterproductive.
              ! The excess energy is simply removed (non-conservative but stable).
           end if
        end if

        ! Safety floor
        esgs_new = max(esgs_new, sgs_floor)

        ! Dissipated energy (to be added to thermal energy)
        esgs_dissipated = esgs - esgs_new
        esgs_dissipated = max(esgs_dissipated, 0d0)

        ! =============================================
        ! Update conserved variables
        ! =============================================
        ! Update SGS energy (rho*e_sgs)
        uold(ind_cell(i),isgs) = rho * esgs_new

        ! Add dissipated SGS energy to thermal energy
        ! (only when sgs_hydro=T; when F, energy simply removed)
        if(sgs_hydro) then
           uold(ind_cell(i),ndim+2) = uold(ind_cell(i),ndim+2) &
                & + rho * esgs_dissipated
        end if

        ! =============================================
        ! 4. Momentum source: -grad(P_turb) now handled in Riemann solver
        !    P_turb = (2/3)*rho*e_sgs added to Ptot in all Riemann solvers
        !    and CFL via c_eff^2 = gamma*P/rho + (2/3)*e_sgs
        ! =============================================

     end do  ! ngrid

  end do  ! twotondim

end subroutine sgs_source
!###########################################################
!###########################################################
! Initialize SGS energy at restart from non-SGS output
! Sets e_sgs = sgs_e_init * e_thermal for cells where isgs=0
!###########################################################
subroutine sgs_init_restart
  use amr_commons
  use hydro_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif

  integer::ilevel,igrid,ind,iskip,icell,ncache,info
  real(dp)::rho,ek,eint_loc,esgs_init_val
  integer(kind=8)::ninit_loc,ninit_glob

  if(.not. use_sgs) return
  if(isgs <= 0) return
  if(sgs_e_init <= 0d0) return

  ninit_loc = 0

  do ilevel = levelmin, nlevelmax
     ncache = active(ilevel)%ngrid
     do igrid = 1, ncache
        do ind = 1, twotondim
           iskip = ncoarse + (ind-1)*ngridmax
           icell = iskip + active(ilevel)%igrid(igrid)
           if(son(icell) /= 0) cycle

           ! Only initialize cells where e_sgs is effectively zero
           rho = max(uold(icell,1), smallr)
           if(uold(icell,isgs)/rho > 10d0*sgs_floor) cycle

           ! Compute specific internal energy
           ek = 0.5d0*(uold(icell,2)**2 + uold(icell,3)**2 + uold(icell,4)**2)/rho
           eint_loc = max(uold(icell,ndim+2) - ek, smallr*smallc**2) / rho

           ! Set e_sgs = sgs_e_init * e_thermal
           esgs_init_val = sgs_e_init * eint_loc
           uold(icell,isgs) = rho * esgs_init_val
           ninit_loc = ninit_loc + 1
        end do
     end do
  end do

#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(ninit_loc, ninit_glob, 1, MPI_INTEGER8, MPI_SUM, MPI_COMM_WORLD, info)
#else
  ninit_glob = ninit_loc
#endif
  if(myid==1 .and. ninit_glob > 0) then
     write(*,'(A,I12,A,ES10.3)') ' SGS init: seeded ', ninit_glob, &
          ' cells with e_sgs/e_th=', sgs_e_init
  end if

end subroutine sgs_init_restart
