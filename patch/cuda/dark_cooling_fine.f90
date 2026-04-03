!################################################################
!################################################################
! Dark Cooling Fine: Apply dark-sector cooling to DM particles
!
! For Atomic Dark Matter (aDM): each DM particle carries a dark
! internal energy edp. Per-particle cooling via backward Euler.
!
! Pattern follows cooling_fine.kjhan.f90 (grid traversal)
!################################################################
subroutine dark_cooling_fine(ilevel)
  use amr_commons
  use pm_commons
  use dark_cooling_mod, only: dark_cool_implicit
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer,intent(in)::ilevel

  integer::ncache,igrid,ngrid

  if(.not. use_adm) return
  if(numbtot(1,ilevel)==0) return
  if(verbose) write(*,111) ilevel

  ncache = active(ilevel)%ngrid
!$omp parallel do private(igrid,ngrid) schedule(dynamic)
  do igrid = 1, ncache, nvector
     ngrid = MIN(nvector, ncache - igrid + 1)
     call sub_dark_cooling_fine(ilevel, igrid, ngrid)
  end do

111 format('   Entering dark_cooling_fine for level ',I2)
end subroutine dark_cooling_fine
!################################################################
!################################################################
subroutine sub_dark_cooling_fine(ilevel, igrid_start, ngrid)
  use amr_commons
  use pm_commons
  use dark_cooling_mod, only: dark_cool_implicit
  implicit none

  integer,intent(in)::ilevel, igrid_start, ngrid

  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  real(dp)::dx,dx_loc,scale,vol_phys,dt_phys
  real(dp)::rho_D,n_D,mp_D_g,mp_phys
  integer::nx_loc,i,ipart,jpart
  integer,dimension(1:nvector)::ind_grid
  real(dp)::edp_new

  ! Physical constants
  real(dp),parameter::GeV_to_g = 1.78266192d-24

  ! Unit conversions
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

  ! Cell size at this level
  dx = 0.5D0**ilevel
  nx_loc = (icoarse_max-icoarse_min+1)
  scale = boxlen/dble(nx_loc)
  dx_loc = dx*scale

  ! Physical cell volume [cm^3]
  vol_phys = (dx_loc*scale_l/aexp)**3

  ! Physical timestep [s]
  dt_phys = dtnew(ilevel)*scale_t

  ! Dark proton mass [g]
  mp_D_g = adm_mp * GeV_to_g

  ! Fill grid index array
  do i = 1, ngrid
     ind_grid(i) = active(ilevel)%igrid(igrid_start + i - 1)
  end do

  ! Traverse particles per grid, apply per-particle dark cooling
  do i = 1, ngrid
     ipart = headp(ind_grid(i))
     do jpart = 1, numbp(ind_grid(i))
        ! Only DM particles (idp > 0, tp <= 0: ground or excited state)
        if(idp(ipart) > 0 .and. tp(ipart) <= 0.0d0) then
           ! Particle physical mass [g]
           ! In RAMSES: mp [code] = rho * vol in code units
           ! mp_phys = mp * scale_d * scale_l^3
           mp_phys = mp(ipart) * scale_d * scale_l**3

           ! Estimate local dark density: particle mass / cell volume [g/cm^3]
           rho_D = mp_phys / vol_phys
           n_D = rho_D / mp_D_g

           if(n_D > 0.0d0) then
              ! Apply dark cooling (edp in code energy/mass units)
              edp_new = dark_cool_implicit( &
                   edp(ipart) * scale_v**2, &  ! to physical [erg/g]
                   rho_D, n_D, dt_phys, aexp)
              edp(ipart) = edp_new / scale_v**2  ! back to code units
           end if
        end if

        ipart = nextp(ipart)
     end do
  end do

end subroutine sub_dark_cooling_fine
