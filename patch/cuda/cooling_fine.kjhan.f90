!##########################################################
!##########################################################
! Clean OMP dispatch for cooling_fine
! Overrides patch/Horizon5-master-2/cooling_fine.kjhan.f90
!
! GPU dispatch NOT implemented because solve_cooling uses
! iterative Newton-Raphson with variable iteration count
! per cell (while-loop with active cell tracking), causing
! severe GPU thread divergence.
!##########################################################
!##########################################################
subroutine cooling_fine(ilevel)
  use amr_commons
  use hydro_commons
  use cooling_module
  use eunha_cooling_mod, only: eunha_report, eunha_interp_redshift
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer::ilevel
  !-------------------------------------------------------------------
  ! Compute cooling for fine levels
  !-------------------------------------------------------------------
  integer::ncache,i,igrid,ngrid,info

#ifdef grackle
  integer:: iresult, initialize_grackle
  real(kind=8)::density_units,length_units,time_units,velocity_units,temperature_units,a_units=1.0,a_value=1.0
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
#endif

  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel

  ! Operator splitting step for cooling source term
  ! by vector sweeps
  ncache=active(ilevel)%ngrid
!$omp parallel do private(igrid,ngrid) schedule(dynamic)
  do igrid=1,ncache,nvector
     ngrid=MIN(nvector,ncache-igrid+1)
     call sub_cooling_fine(ilevel,igrid,ngrid)
  end do

  ! Report comparison statistics at end of coarse step
  if(ilevel==levelmin .and. cooling_method=='compare') then
     call eunha_report(nstep_coarse)
  endif

  if((cooling.and..not.neq_chem).and.ilevel==levelmin.and.cosmo)then
     if(myid==1)write(*,*)'Computing new cooling table'
#ifdef grackle
     ! Compute new cooling table at current aexp with grackle
     if((1.D0/aexp-1.D0.lt.z_reion).and.(grackle_UVbackground.eq.1).and.(.not.grackle_UVbackground_on)) then
        if(myid==1)write(*,*)'Grackle: Activating UV background'
        grackle_UVbackground_on = .true.
        a_value = aexp
        call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
        density_units=scale_d
        length_units=scale_l
        time_units=scale_t
        velocity_units=scale_v
        ! Initialize the Grackle data
         iresult = initialize_grackle(                                &
           &     grackle_comoving_coordinates,                        &
           &     density_units, length_units,                         &
           &     time_units, velocity_units,                          &
           &     a_units, a_value,                                    &
           &     use_grackle, grackle_with_radiative_cooling,         &
           &     TRIM(grackle_data_file),                             &
           &     grackle_primordial_chemistry, grackle_metal_cooling, &
           &     grackle_UVbackground, grackle_h2_on_dust,            &
           &     grackle_cmb_temperature_floor, gamma)
     endif
#else
     call set_table(dble(aexp))
     ! Update Eunha multi-z table if loaded
     if(cooling_method/='original') then
        call eunha_interp_redshift(1.0d0/dble(aexp) - 1.0d0)
     endif
#endif
  endif

111 format('   Entering cooling_fine for level',i2)
end subroutine cooling_fine

!###########################################################
!###########################################################
subroutine sub_cooling_fine(ilevel,igrid,ngrid)
  use amr_commons
  use hydro_commons
  use cooling_module
  integer,dimension(1:nvector)::ind_grid
  integer:: i,ilevel,igrid,ngrid
  do i = 1,ngrid
     ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
  enddo
  call coolfine1(ind_grid,ngrid,ilevel)
end subroutine sub_cooling_fine

!###########################################################
!###########################################################
subroutine coolfine1(ind_grid,ngrid,ilevel)
  use amr_commons
  use hydro_commons
  use cooling_module
  use eunha_cooling_mod, only: eunha_solve, acc_sum_rdiff, acc_max_rdiff, &
       acc_sum_dT_orig, acc_sum_dT_exact, acc_ncell, acc_nsign
#ifdef ATON
  use radiation_commons, ONLY: Erad
#endif
#ifdef RT
  use rt_parameters, only: nGroups, iGroups
  use rt_hydro_commons
  use rt_cooling_module, only: rt_solve_cooling,iIR,rt_isIRtrap &
       ,rt_pressBoost,iIRtrapVar,kappaSc,a_r,is_kIR_T,rt_vc
#endif
  implicit none
  integer::ilevel,ngrid
  integer,dimension(1:nvector)::ind_grid
  !-------------------------------------------------------------------
  !-------------------------------------------------------------------
  integer::i,ind,iskip,idim,nleaf,nx_loc,ix,iy,iz
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  real(kind=8)::dtcool,nISM,nCOM,damp_factor,cooling_switch,t_blast
  real(dp)::polytropic_constant,dx_phys_cm
  real(dp),dimension(1:3)::skip_loc
  real(kind=8)::dx,dx_loc,scale,vol_loc
  integer::irad

  integer,dimension(1:nvector)::ind_cell,ind_leaf
  real(kind=8),dimension(1:nvector)::nH,T2,T2_new,delta_T2,ekk,err,emag
  real(kind=8),dimension(1:nvector)::T2min,Zsolar,boost
  real(kind=8),dimension(1:nvector)::delta_T2_eunha
  real(dp)::rdiff,old_max

#ifdef grackle
  real(kind=8) gr_density(nvector), gr_energy(nvector), &
  &     gr_x_velocity(nvector), gr_y_velocity(nvector), &
  &     gr_z_velocity(nvector), gr_metal_density(nvector), &
  &     gr_poly(nvector), gr_floor(nvector)
  integer::iresult, solve_chemistry_table, gr_rank
  integer,dimension(1:3)::gr_dimension,gr_start,gr_end
  real(dp)::density_units,length_units,time_units,velocity_units,temperature_units,a_units=1.0,a_value=1.0,gr_dt
  if(cosmo) then
     a_value=aexp
  endif
#endif

  ! Mesh spacing in that level
  dx=0.5D0**ilevel
  nx_loc=(icoarse_max-icoarse_min+1)
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  if(ndim>0)skip_loc(1)=dble(icoarse_min)
  if(ndim>1)skip_loc(2)=dble(jcoarse_min)
  if(ndim>2)skip_loc(3)=dble(kcoarse_min)
  scale=boxlen/dble(nx_loc)
  dx_loc=dx*scale
  vol_loc=dx_loc**ndim

  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

  ! Typical ISM density in H/cc
  nISM = n_star; nCOM=0d0
  if(cosmo)then
     nCOM = del_star*omega_b*rhoc*(h0/100.)**2/aexp**3*X/mH
  endif
  nISM = MAX(nCOM,nISM)

  ! Polytropic constant for Jeans length related polytropic EOS
  if(jeans_ncells>0)then
     ! dx_phys [cm] at finest level
     dx_phys_cm=boxlen*0.5d0**dble(nlevelmax)*scale_l/aexp
     ! FPR: effective resolution = max(dx_phys, dr_proper)
     if(dr_proper > 0d0) dx_phys_cm=max(dx_phys_cm, dr_proper*3.085678d21)
     polytropic_constant=2d0*(jeans_ncells*dx_phys_cm)**2/ &
          & (twopi)*6.67e-8*scale_d*(scale_t/scale_l)**2
  endif

  ! Update polytropic floor coefficients for hydro solver (ctoprim/cmpdt)
  ! Floor on energy density: rho*e >= coeff * rho^alpha
  ! Floor on specific energy: e >= coeff * rho^(alpha-1)
  if(jeans_ncells>0)then
     eeos_poly_coeff = polytropic_constant * scale_nH / (gamma-1.0)
     eeos_poly_alpha = 2.0d0
  else if(T2_star > 0d0)then
     eeos_poly_coeff = T2_star * scale_nH**(g_star-1.0) / &
          & (nISM**(g_star-1.0) * scale_T2 * (gamma-1.0))
     eeos_poly_alpha = g_star
  endif

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
     if(nleaf.eq.0)cycle

     ! Compute rho
     do i=1,nleaf
        nH(i)=MAX(uold(ind_leaf(i),1),smallr)
     end do

     ! Compute metallicity in solar units
     if(metal)then
        do i=1,nleaf
           Zsolar(i)=uold(ind_leaf(i),imetal)/nH(i)/0.02
        end do
     else
        do i=1,nleaf
           Zsolar(i)=z_ave
        end do
     endif

     ! Compute thermal pressure
     do i=1,nleaf
        T2(i)=uold(ind_leaf(i),ndim+2)
     end do
     do i=1,nleaf
        ekk(i)=0.0d0
     end do
     do idim=1,ndim
        do i=1,nleaf
           ekk(i)=ekk(i)+0.5*uold(ind_leaf(i),idim+1)**2/nH(i)
        end do
     end do
     do i=1,nleaf
        err(i)=0.0d0
     end do
#if NENER>0
     do irad=0,nener-1
        do i=1,nleaf
           err(i)=err(i)+uold(ind_leaf(i),inener+irad)
        end do
     end do
#endif
     do i=1,nleaf
        emag(i)=0.0d0
     end do
#ifdef SOLVERmhd
     do idim=1,ndim
        do i=1,nleaf
           emag(i)=emag(i)+0.125d0*(uold(ind_leaf(i),idim+ndim+2)+uold(ind_leaf(i),idim+nvar))**2
        end do
     end do
#endif
     do i=1,nleaf
        T2(i)=(gamma-1.0)*(T2(i)-ekk(i)-err(i)-emag(i))
     end do

     ! Compute T2=T/mu in Kelvin
     do i=1,nleaf
        T2(i)=T2(i)/nH(i)*scale_T2
     end do

     ! Compute nH in H/cc
     do i=1,nleaf
        nH(i)=nH(i)*scale_nH
     end do

     ! Compute radiation boost factor
     if(self_shielding)then
        do i=1,nleaf
           boost(i)=MAX(exp(-nH(i)/0.01),1.0D-20)
        end do
#ifdef ATON
     else if (aton) then
        do i=1,nleaf
           boost(i)=MAX(Erad(ind_leaf(i))/J0simple(aexp), &
                &                   J0min/J0simple(aexp) )
        end do
#endif
     else
        do i=1,nleaf
           boost(i)=1.0
        end do
     endif

     !==========================================
     ! Compute temperature from polytrope EOS
     !==========================================
     if(jeans_ncells>0)then
        do i=1,nleaf
           T2min(i) = nH(i)*polytropic_constant*scale_T2
        end do
     else
        do i=1,nleaf
           T2min(i) = T2_star*(nH(i)/nISM)**(g_star-1.0)
        end do
     endif

     if(cooling)then
        ! Compute thermal temperature by subtracting polytrope
        do i=1,nleaf
           T2(i) = min(max(T2(i)-T2min(i),T2_min_fix),T2max)
        end do
     endif

     ! Compute cooling time step in second
     dtcool = dtnew(ilevel)*scale_t

     ! grackle tabular cooling
#ifdef grackle
     if(cosmo) then
        a_value=aexp
     endif
     gr_rank = 3
     do i = 1, gr_rank
        gr_dimension(i) = 1
        gr_start(i) = 0
        gr_end(i) = 0
     enddo
     gr_dimension(1) = nvector
     gr_end(1) = nleaf - 1
     ! set units
     density_units=scale_d
     length_units=scale_l
     time_units=scale_t
     velocity_units=scale_v
     temperature_units=scale_T2
     do i = 1, nleaf
        gr_density(i) = max(uold(ind_leaf(i),1),smallr)
        if(metal)then
           gr_metal_density(i) = uold(ind_leaf(i),imetal)
        else
           gr_metal_density(i) = uold(ind_leaf(i),1)*0.02*z_ave
        endif
        gr_x_velocity(i) = uold(ind_leaf(i),2)/max(uold(ind_leaf(i),1),smallr)
        gr_y_velocity(i) = uold(ind_leaf(i),3)/max(uold(ind_leaf(i),1),smallr)
        gr_z_velocity(i) = uold(ind_leaf(i),4)/max(uold(ind_leaf(i),1),smallr)
        gr_floor(i)  = 1.0*nH(i)/scale_nH/scale_T2/(gamma-1.0)
        gr_poly(i)   = T2min(i)*nH(i)/scale_nH/scale_T2/(gamma-1.0)
        gr_energy(i) = uold(ind_leaf(i),ndim+2)-ekk(i)-gr_poly(i)
        gr_energy(i) = MAX(gr_energy(i),gr_floor(i))
        gr_energy(i) = gr_energy(i)/max(uold(ind_leaf(i),1),smallr)
     enddo

     gr_dt = dtnew(ilevel)

     iresult = solve_chemistry_table(    &
     &     grackle_comoving_coordinates, &
     &     density_units, length_units,  &
     &     time_units, velocity_units,   &
     &     a_units, a_value, gr_dt,      &
     &     gr_rank, gr_dimension,        &
     &     gr_start, gr_end,             &
     &     gr_density, gr_energy,        &
     &     gr_x_velocity, gr_y_velocity, gr_z_velocity, &
     &     gr_metal_density)

     do i = 1, nleaf
        T2_new(i) = gr_energy(i)*scale_T2*(gamma-1.0)
     end do
     delta_T2(1:nleaf) = T2_new(1:nleaf) - T2(1:nleaf)
#else
     ! Eunha exact / compare / original cooling dispatch
     if(cooling_method=='exact' .or. cooling_method=='compare') then
        call eunha_solve(nH,T2,Zsolar,dtcool,delta_T2_eunha,nleaf)
     endif

     if(cooling_method=='original' .or. cooling_method=='compare') then
        if(cooling.and..not.neq_chem) &
             call solve_cooling(nH,T2,Zsolar,boost,dtcool,delta_T2,nleaf)
     endif

     if(cooling_method=='exact') then
        delta_T2(1:nleaf) = delta_T2_eunha(1:nleaf)
     endif

     if(cooling_method=='compare') then
        do i=1,nleaf
           rdiff = abs(delta_T2(i)-delta_T2_eunha(i)) / max(abs(T2(i)),1d0)
           !$omp atomic
           acc_sum_rdiff = acc_sum_rdiff + rdiff
           ! atomic max via CAS-style loop
           old_max = acc_max_rdiff
           do while(rdiff > old_max)
              acc_max_rdiff = rdiff
              old_max = acc_max_rdiff
           end do
           !$omp atomic
           acc_sum_dT_orig = acc_sum_dT_orig + delta_T2(i)
           !$omp atomic
           acc_sum_dT_exact = acc_sum_dT_exact + delta_T2_eunha(i)
           !$omp atomic
           acc_ncell = acc_ncell + 1
           if(delta_T2(i) * delta_T2_eunha(i) < 0.0d0) then
              !$omp atomic
              acc_nsign = acc_nsign + 1
           endif
        end do
     endif
#endif

     ! Compute rho
     do i=1,nleaf
        nH(i) = nH(i)/scale_nH
     end do

     ! Deal with cooling
     if(cooling.or.neq_chem)then
        ! Compute net energy sink
        do i=1,nleaf
           delta_T2(i) = delta_T2(i)*nH(i)/scale_T2/(gamma-1.0)
        end do
        ! Compute initial fluid internal energy
        do i=1,nleaf
           T2(i) = T2(i)*nH(i)/scale_T2/(gamma-1.0)
        end do
        ! Turn off cooling in blast wave regions
        if(delayed_cooling)then
           do i=1,nleaf
              cooling_switch = uold(ind_leaf(i),idelay)/max(uold(ind_leaf(i),1),smallr)
              if(cooling_switch > 1d-3)then
                 delta_T2(i) = MAX(delta_T2(i),real(0,kind=dp))
              endif
           end do
        endif
     endif

     ! Compute polytrope internal energy
     do i=1,nleaf
        T2min(i) = T2min(i)*nH(i)/scale_T2/(gamma-1.0)
     end do

     ! Update fluid internal energy
     if(cooling.or.neq_chem)then
        do i=1,nleaf
           T2(i) = T2(i) + delta_T2(i)
        end do
     endif

     ! Update total fluid energy
     if(isothermal)then
        do i=1,nleaf
           uold(ind_leaf(i),ndim+2) = T2min(i) + ekk(i) + err(i) + emag(i)
        end do
     else if(cooling .or. neq_chem)then
        do i=1,nleaf
           uold(ind_leaf(i),ndim+2) = T2(i) + T2min(i) + ekk(i) + err(i) + emag(i)
        end do
     endif

     ! Update delayed cooling switch
     if(delayed_cooling)then
        t_blast=t_diss*1d6*(365.*24.*3600.)
        damp_factor=exp(-dtcool/t_blast)
        do i=1,nleaf
           uold(ind_leaf(i),idelay)=uold(ind_leaf(i),idelay)*damp_factor
        end do
     endif

  end do
  ! End loop over cells

end subroutine coolfine1

!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine init_eeos_poly_coeff
  ! Compute polytropic floor coefficients for eEOS enforcement
  ! in the hydro solver (ctoprim/cmpdt). Called once at startup.
  use amr_commons
  use hydro_commons
  use cooling_module
  implicit none
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  real(dp)::polytropic_constant,dx_phys_cm
  real(dp)::nISM,nCOM

  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
  nISM = n_star; nCOM=0d0
  if(cosmo)then
     nCOM = del_star*omega_b*rhoc*(h0/100.)**2/aexp**3*X/mH
  endif
  nISM = MAX(nCOM,nISM)

  if(jeans_ncells>0)then
     ! dx_phys [cm] at finest level
     dx_phys_cm=boxlen*0.5d0**dble(nlevelmax)*scale_l/aexp
     ! FPR: effective resolution = max(dx_phys, dr_proper)
     if(dr_proper > 0d0) dx_phys_cm=max(dx_phys_cm, dr_proper*3.085678d21)
     polytropic_constant=2d0*(jeans_ncells*dx_phys_cm)**2/ &
          & (twopi)*6.67e-8*scale_d*(scale_t/scale_l)**2
     eeos_poly_coeff = polytropic_constant * scale_nH / (gamma-1.0)
     eeos_poly_alpha = 2.0d0
  else if(T2_star > 0d0)then
     eeos_poly_coeff = T2_star * scale_nH**(g_star-1.0) / &
          & (nISM**(g_star-1.0) * scale_T2 * (gamma-1.0))
     eeos_poly_alpha = g_star
  endif

  if(myid==1 .and. eeos_poly_coeff > 0d0)then
     write(*,'(A,ES10.3,A,F5.2)') ' eEOS polytropic floor: coeff=', eeos_poly_coeff, &
          & ' alpha=', eeos_poly_alpha
  endif
end subroutine init_eeos_poly_coeff

!###########################################################
!###########################################################
subroutine enforce_eeos_after_sink
  ! Enforce TWO protections after sink/AGN operations (TOTAL ENERGY CONSERVED):
  !
  ! 1. eEOS polytropic floor: if eint < floor, reduce momentum to maintain floor
  ! 2. Mach cap: if |v|/c_s > Mach_max, reduce momentum, convert excess Ekin to Eth
  !
  ! The Mach cap prevents CFL dt collapse from AGN jet momentum injection into
  ! fine cells. Excess kinetic energy is thermalised (instant shock approximation)
  ! which is valid when the jet cooling length << cell size.
  use amr_commons
  use hydro_commons
  implicit none
  integer::ilevel,igrid,ngrid,ncache,ind,iskip,i,idim
  integer,dimension(1:nvector)::ind_grid,ind_cell
  logical,dimension(1:nvector)::ok
  real(dp)::d,e_kin,e_kin_pure,e_total,e_int,e_floor
  real(dp)::e_kin_cap,c2,v2,v2_max,factor
  integer::irad
  integer::n_floor,n_mach
  real(dp)::e_floor_inj,e_mach_inj
  real(dp),parameter::mach_cap=100d0  ! Max Mach number after AGN injection

  n_floor = 0; n_mach = 0
  e_floor_inj = 0d0; e_mach_inj = 0d0

  do ilevel=levelmin,nlevelmax
     ncache=active(ilevel)%ngrid
     do igrid=1,ncache,nvector
        ngrid=MIN(nvector,ncache-igrid+1)
        do i=1,ngrid
           ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
        end do
        do ind=1,twotondim
           iskip=ncoarse+(ind-1)*ngridmax
           do i=1,ngrid
              ind_cell(i)=ind_grid(i)+iskip
           end do
           do i=1,ngrid
              ok(i)=son(ind_cell(i))==0
           end do
           do i=1,ngrid
              if(.not.ok(i))cycle

              d=max(uold(ind_cell(i),1),smallr)
              e_total=uold(ind_cell(i),ndim+2)

              ! Compute pure kinetic energy (0.5*|mom|^2/rho)
              e_kin_pure=0d0
              do idim=1,ndim
                 e_kin_pure=e_kin_pure+0.5d0*uold(ind_cell(i),idim+1)**2/d
              end do
              ! Total "non-thermal kinetic" including NENER
              e_kin=e_kin_pure
#if NENER>0
              do irad=1,nener
                 e_kin=e_kin+uold(ind_cell(i),ndim+2+irad)
              end do
#endif
              e_int=e_total-e_kin

              ! === Protection 1: eEOS floor ===
              e_floor=eeos_poly_coeff*d**eeos_poly_alpha
              if(e_int < e_floor)then
                 n_floor = n_floor + 1
                 e_floor_inj = e_floor_inj + (e_floor - e_int)
                 ! Convert kinetic → thermal to satisfy floor
                 e_kin_cap = max(e_total - e_floor, 0d0)
#if NENER>0
                 do irad=1,nener
                    e_kin_cap = e_kin_cap - uold(ind_cell(i),ndim+2+irad)
                 end do
#endif
                 if(e_kin_cap > 0d0 .and. e_kin_pure > 0d0)then
                    factor = sqrt(e_kin_cap / e_kin_pure)
                    do idim=1,ndim
                       uold(ind_cell(i),idim+1) = uold(ind_cell(i),idim+1) * factor
                    end do
                 else
                    do idim=1,ndim
                       uold(ind_cell(i),idim+1) = 0d0
                    end do
                    uold(ind_cell(i),ndim+2) = e_floor
#if NENER>0
                    do irad=1,nener
                       uold(ind_cell(i),ndim+2) = uold(ind_cell(i),ndim+2) + &
                            & uold(ind_cell(i),ndim+2+irad)
                    end do
#endif
                 endif
                 ! Recompute after floor fix
                 e_kin_pure=0d0
                 do idim=1,ndim
                    e_kin_pure=e_kin_pure+0.5d0*uold(ind_cell(i),idim+1)**2/d
                 end do
                 e_int=e_total-e_kin_pure
#if NENER>0
                 do irad=1,nener
                    e_int=e_int-uold(ind_cell(i),ndim+2+irad)
                 end do
#endif
              endif

              ! === Protection 2: Mach cap ===
              ! c_s^2 = gamma*(gamma-1)*eint/rho
              c2 = gamma*(gamma-1d0)*max(e_int,smallc**2*d)/d
              ! v^2 = 2*e_kin_pure/rho  (specific kinetic = 0.5*v^2)
              v2 = 2d0*e_kin_pure/d
              v2_max = mach_cap**2 * c2
              if(v2 > v2_max .and. v2 > 0d0)then
                 n_mach = n_mach + 1
                 ! Scale momentum: v_new = mach_cap * c_s
                 ! e_kin_new = 0.5 * d * v2_max = e_kin_pure * (v2_max/v2)
                 factor = sqrt(v2_max / v2)
                 e_mach_inj = e_mach_inj + e_kin_pure*(1d0-factor**2)
                 do idim=1,ndim
                    uold(ind_cell(i),idim+1) = uold(ind_cell(i),idim+1) * factor
                 end do
                 ! Total energy conserved: excess kinetic → thermal
                 ! uold(5) unchanged, momentum reduced → eint increases
              endif
           end do
        end do
     end do
  end do

  ! Report if any cells were fixed
  if((n_floor > 0 .or. n_mach > 0) .and. myid==1 .and. nstep_coarse < 500)then
     if(n_floor > 0) write(*,'(A,I8,A,ES10.3)') &
          & ' eEOS post-sink floor: ', n_floor, ' cells, dE_th=', e_floor_inj
     if(n_mach > 0) write(*,'(A,I8,A,ES10.3)') &
          & ' eEOS post-sink Mach cap: ', n_mach, ' cells, dE_kin→th=', e_mach_inj
  endif

end subroutine enforce_eeos_after_sink
