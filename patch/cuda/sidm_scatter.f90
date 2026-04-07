!################################################################
!################################################################
!################################################################
! SIDM (Self-Interacting Dark Matter) Monte Carlo Scattering
!
! Features:
! - Velocity-dependent cross-section (constant/yukawa/power_law)
! - Isotropic or anisotropic (Rutherford-like) scattering angle
! - Inelastic scattering with mass splitting (iSIDM)
! - Timestep constraint via P_max tracking
!################################################################
subroutine sidm_scatter(ilevel)
  use pm_commons
  use pm_parameters, only: iseed
  use amr_commons
  use random
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer,intent(in)::ilevel

  integer::icpu,igrid,subnump,info,ith
  integer::mythread,nthreads
  integer,dimension(:),allocatable::nparticles,ptrhead
  integer::sidm_n_scatter,sidm_n_pairs,sidm_n_up,sidm_n_down
  real(dp)::sidm_dp(3),sidm_dEk,sidm_dp_max,sidm_dEk_max
  real(dp)::sidm_Pmax_level
  real(dp)::sidm_dEdiss  ! Cumulative dissipated energy (dSIDM)
  real(dp)::R_dummy
  integer,dimension(1:ncpu,1:IRandNumSize)::allseed
  ! Per-thread seed array (deterministic, avoids localseed race)
  integer,dimension(:,:),allocatable::sidm_seeds
  logical,save::sidm_excited_init = .false.
  common /sidm_omp/ mythread
  common /sidm_diag/ sidm_n_scatter,sidm_n_pairs,sidm_n_up,sidm_n_down
  common /sidm_cons/ sidm_dp,sidm_dEk,sidm_dp_max,sidm_dEk_max
  common /sidm_pmax_c/ sidm_Pmax_level
  common /sidm_diss/ sidm_dEdiss
!$omp threadprivate(/sidm_omp/)

  if(.not.sidm) return
  if(numbtot(1,ilevel)==0) return
  if(verbose) write(*,111) ilevel

  ! Initialize excited DM fraction (once, for iSIDM)
  if(sidm_inelastic .and. .not.sidm_excited_init) then
     call sidm_init_excited()
     sidm_excited_init = .true.
  end if

!$omp parallel
  mythread = omp_get_thread_num()
  if(mythread==0) nthreads = omp_get_num_threads()
!$omp end parallel
  allocate(ptrhead(0:nthreads-1), nparticles(0:nthreads-1))

  ! Initialize localseed if first call
  if(localseed(1)==-1) then
     call rans(ncpu,iseed,allseed)
     localseed = allseed(myid,1:IRandNumSize)
  end if

  ! Initialize per-thread seeds from localseed (deterministic)
  allocate(sidm_seeds(IRandNumSize, 0:nthreads-1))
  sidm_seeds(:,0) = localseed
  do ith=1,nthreads-1
     sidm_seeds(:,ith) = sidm_seeds(:,ith-1)
     ! Advance seed 100 steps to decorrelate streams
     do igrid=1,100
        call ranf(sidm_seeds(:,ith), R_dummy)
     end do
  end do

  ! Reset level-wide counters
  sidm_n_scatter  = 0
  sidm_n_pairs    = 0
  sidm_n_up       = 0
  sidm_n_down     = 0
  sidm_dp(:)      = 0.0d0
  sidm_dEk        = 0.0d0
  sidm_dp_max     = 0.0d0
  sidm_dEk_max    = 0.0d0
  sidm_Pmax_level = 0.0d0
  sidm_dEdiss     = 0.0d0

#if NDIM==3
  do icpu=1,ncpu
     if(numbl(icpu,ilevel)<=0) cycle
     call pthreadLinkedList(headl(icpu,ilevel),numbl(icpu,ilevel), &
          nthreads,nparticles,ptrhead,next)
!$omp parallel private(subnump,igrid)
     subnump = nparticles(mythread)
     igrid   = ptrhead(mythread)
     call sub_sidm_scatter(ilevel,icpu,igrid,subnump,sidm_seeds)
!$omp end parallel
  end do
#endif

  ! Update localseed deterministically: advance last thread's seed
  localseed = sidm_seeds(:,nthreads-1)
  do igrid=1,100
     call ranf(localseed, R_dummy)
  end do

  deallocate(ptrhead, nparticles, sidm_seeds)

  ! MPI reduce all diagnostics
  call MPI_ALLREDUCE(MPI_IN_PLACE,sidm_n_scatter,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(MPI_IN_PLACE,sidm_n_pairs,  1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(MPI_IN_PLACE,sidm_n_up,     1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(MPI_IN_PLACE,sidm_n_down,   1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(MPI_IN_PLACE,sidm_dp,   3,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(MPI_IN_PLACE,sidm_dEk,  1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(MPI_IN_PLACE,sidm_dp_max, 1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(MPI_IN_PLACE,sidm_dEk_max,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(MPI_IN_PLACE,sidm_Pmax_level,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(MPI_IN_PLACE,sidm_dEdiss,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)

  ! Store P_max for timestep constraint in newdt_fine
  sidm_Pmax(ilevel) = sidm_Pmax_level

  ! Report statistics
  if(myid==1 .and. sidm_n_pairs>0) then
     write(*,'(A,I2,A,I8,A,I10,A,ES9.2)') &
          ' SIDM level ',ilevel,': scattered=',sidm_n_scatter, &
          ' pairs=',sidm_n_pairs,' Pmax=',sidm_Pmax_level
     if(sidm_n_scatter>0) then
        write(*,'(A,3ES12.4)') '   dp(x,y,z)=', sidm_dp(1:3)
        write(*,'(A,ES12.4,A,ES12.4)') '   dEk_total=', sidm_dEk, &
             '  dp_max=', sidm_dp_max
        write(*,'(A,ES12.4)') '   dEk_max  =', sidm_dEk_max
     end if
     if(sidm_inelastic .and. (sidm_n_up>0 .or. sidm_n_down>0)) then
        write(*,'(A,I6,A,I6,A,I6)') '   iSIDM: up=', &
             sidm_n_up,' down=',sidm_n_down, &
             ' net=',sidm_n_up-sidm_n_down
     end if
     if(sidm_fdiss > 0.0d0 .and. sidm_dEdiss /= 0.0d0) then
        write(*,'(A,ES12.4)') '   dEdiss_total=', sidm_dEdiss
     end if
  end if

111 format('   Entering sidm_scatter for level ',I2)
end subroutine sidm_scatter
!################################################################
!################################################################
subroutine sub_sidm_scatter(ilevel,icpu,kgrid,subnump,thread_seeds)
  use pm_commons
  use amr_commons
  use random
  implicit none
  integer,intent(in)::ilevel,icpu,kgrid,subnump
  integer,dimension(IRandNumSize,0:*)::thread_seeds

  ! External function
  integer,external::cell_index_from_part

  ! Shared counters (common with sidm_scatter)
  integer::sidm_n_scatter,sidm_n_pairs,sidm_n_up,sidm_n_down
  real(dp)::sidm_dp(3),sidm_dEk,sidm_dp_max,sidm_dEk_max
  real(dp)::sidm_Pmax_level
  real(dp)::sidm_dEdiss
  common /sidm_diag/ sidm_n_scatter,sidm_n_pairs,sidm_n_up,sidm_n_down
  common /sidm_cons/ sidm_dp,sidm_dEk,sidm_dp_max,sidm_dEk_max
  common /sidm_pmax_c/ sidm_Pmax_level
  common /sidm_diss/ sidm_dEdiss

  ! Thread ID from sidm_scatter wrapper
  integer::mythread_loc
  common /sidm_omp/ mythread_loc
!$omp threadprivate(/sidm_omp/)

  ! Local variables
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  real(dp)::dx,dx_loc,scale,vol_phys,dt_phys,sigma_over_m
  real(dp)::mp_phys,P_scatter,twopi
  real(dp)::cos_theta,sin_theta,phi_rand,v_rel_mag,v_rel_kms
  real(dp),dimension(1:3)::v1,v2,v_cm,v_rel_vec,nhat,v_rel_new
  real(dp),dimension(1:3)::p_before,p_after,v1_new,v2_new
  real(dp)::Ek_before,Ek_after,dp_mag,dEk_evt
  real(dp)::m1,m2,mtot,R1,R2
  integer::igrid,jgrid,ind,iskip,icell,nx_loc
  integer::ipart,jpart,npart1,ndm_cell,ip,jp,npairs
  integer,dimension(1:nvector)::ind_dm
  integer::itemp,ipair
  integer,dimension(IRandNumSize)::seed_loc

  ! Scattering counters (thread-local)
  integer::n_scatter_loc,n_pairs_loc,n_up_loc,n_down_loc
  ! Conservation diagnostics (thread-local)
  real(dp)::dp_loc(3),dEk_loc,dp_max_loc,dEk_max_loc
  real(dp)::P_max_loc
  real(dp)::dEdiss_loc,KE_cm_before_diss

  ! Anisotropic scattering variables
  real(dp)::eps2,a_ruth,b_ruth,e1_mag
  real(dp),dimension(1:3)::v_hat,e1,e2

  ! Inelastic scattering variables
  real(dp)::mu,KE_cm,KE_cm_new,delta_KE,v_rel_new_mag
  logical::do_transition
  integer::istate_1,istate_2,new_state_1,new_state_2
  integer::n_channels,ichan,j1,j2
  real(dp)::R_inel,Ei_1,Ei_2,dE
  ! Multi-state channel arrays (max N*(N+1)/2 = 55 for N=10)
  integer::chan_j1(55),chan_j2(55)
  real(dp)::chan_dE(55)

  if(subnump<=0) return

  ! Unit conversions
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

  twopi = 2.0d0*ACOS(-1.0d0)

  ! Cell size at this level
  dx = 0.5D0**ilevel
  nx_loc = (icoarse_max-icoarse_min+1)
  scale = boxlen/dble(nx_loc)
  dx_loc = dx*scale

  ! Physical cell volume [cm^3]
  vol_phys = (dx_loc*scale_l/aexp)**3

  ! Physical timestep [s]
  dt_phys = dtnew(ilevel)*scale_t

  ! (Multi-state energies are accessed directly as sidm_energy(:) [keV]
  !  and converted to erg inline: E [erg] = E [keV] * 1.602d-9)

  ! Precompute Rutherford parameters
  if(trim(sidm_angular) == 'rutherford') then
     eps2 = 2.0d0*sidm_epsilon
     a_ruth = 1.0d0/eps2
     b_ruth = 1.0d0/(2.0d0+eps2)
  end if

  ! Use per-thread seed from thread_seeds array
  seed_loc = thread_seeds(:,mythread_loc)

  n_scatter_loc = 0
  n_pairs_loc   = 0
  n_up_loc      = 0
  n_down_loc    = 0
  dp_loc(:)     = 0.0d0
  dEk_loc       = 0.0d0
  dp_max_loc    = 0.0d0
  dEk_max_loc   = 0.0d0
  dEdiss_loc    = 0.0d0
  P_max_loc     = 0.0d0

  ! Loop over grids assigned to this thread
  igrid = kgrid
  do jgrid=1,subnump

     npart1 = numbp(igrid)
     if(npart1<sidm_npart_min) then
        igrid = next(igrid)
        cycle
     end if

     ! Loop over 8 subcells (twotondim=8 in 3D)
     do ind=1,twotondim
        iskip = ncoarse+(ind-1)*ngridmax
        icell = iskip+igrid

        ! Only process leaf cells
        if(son(icell)/=0) cycle

        ! Collect DM particles in this subcell
        ! DM: idp>0 and tp<=0 (ground: tp=0, excited: tp=-1)
        ndm_cell = 0
        ipart = headp(igrid)
        do jp=1,npart1
           ! Check if this particle belongs to subcell 'ind'
           if(cell_index_from_part(ipart,igrid,ilevel)==ind) then
              if(idp(ipart)>0 .and. tp(ipart)<=0.0d0) then
                 ndm_cell = ndm_cell+1
                 if(ndm_cell<=nvector) ind_dm(ndm_cell) = ipart
              end if
           end if
           ipart = nextp(ipart)
        end do

        ! Need at least sidm_npart_min DM particles
        if(ndm_cell<sidm_npart_min) cycle

        ! Clamp to nvector
        if(ndm_cell>nvector) ndm_cell = nvector

        ! Fisher-Yates shuffle
        do ip=ndm_cell,2,-1
           call ranf(seed_loc, R1)
           jp = 1 + int(R1*dble(ip))
           if(jp>ip) jp = ip  ! safety clamp
           if(jp/=ip) then
              itemp = ind_dm(ip)
              ind_dm(ip) = ind_dm(jp)
              ind_dm(jp) = itemp
           end if
        end do

        ! Sequential pairing: (1,2), (3,4), ...
        npairs = ndm_cell/2
        n_pairs_loc = n_pairs_loc + npairs

        do ip=1,npairs
           ! Pair index into ind_dm array
           ipair = 2*ip-1
           ! Particle masses (physical) [g]
           m1 = mp(ind_dm(ipair))   * scale_d * scale_l**3
           m2 = mp(ind_dm(ipair+1)) * scale_d * scale_l**3
           mtot = m1 + m2

           ! Particle velocities (code units -> physical [cm/s])
           v1(1:3) = vp(ind_dm(ipair),   1:3) * scale_v / aexp
           v2(1:3) = vp(ind_dm(ipair+1), 1:3) * scale_v / aexp

           ! Relative velocity
           v_rel_vec(1:3) = v1(1:3) - v2(1:3)
           v_rel_mag = sqrt(v_rel_vec(1)**2 + v_rel_vec(2)**2 &
                          + v_rel_vec(3)**2)

           if(v_rel_mag==0.0d0) cycle

           ! --- Velocity-dependent cross-section ---
           v_rel_kms = v_rel_mag * 1.0d-5  ! cm/s -> km/s
           select case(trim(sidm_type))
           case('yukawa')
              ! sigma(v) = sigma_0 / (1 + (v/v0)^2)^2
              sigma_over_m = sidm_cross_section &
                   / (1.0d0 + (v_rel_kms/sidm_v0)**2)**2
           case('power_law')
              ! sigma(v) = sigma_0 * (v/v0)^n, capped
              if(v_rel_kms > 1.0d0) then
                 sigma_over_m = sidm_cross_section &
                      * (v_rel_kms/sidm_v0)**sidm_power
              else
                 sigma_over_m = sidm_cross_section &
                      * (1.0d0/sidm_v0)**sidm_power
              end if
              sigma_over_m = min(sigma_over_m, &
                   1.0d4*sidm_cross_section)
           case default  ! 'constant'
              sigma_over_m = sidm_cross_section
           end select

           ! --- Dark phase transition: sigma(v,a) = sigma(v) * g(a) ---
           select case(trim(sidm_a_type))
           case('step')
              if(aexp > sidm_a_transition) then
                 sigma_over_m = sigma_over_m * sidm_sigma_ratio
              end if
           case('sigmoid')
              sigma_over_m = sigma_over_m * (1.0d0 + &
                   (sidm_sigma_ratio - 1.0d0) * 0.5d0 * &
                   (1.0d0 + tanh((aexp - sidm_a_transition)/sidm_a_width)))
           case default  ! 'none'
              ! no modification
           end select

           ! Scattering probability
           ! P = (sigma/m) * m_p * v_rel * dt / V_cell * (N_dm - 1)
           mp_phys = 0.5d0*(m1+m2)  ! representative particle mass
           P_scatter = sigma_over_m * mp_phys * v_rel_mag * dt_phys &
                     / vol_phys * dble(ndm_cell-1)

           ! Track P_max for timestep constraint
           if(P_scatter > P_max_loc) P_max_loc = P_scatter

           ! Warning if P > 1 (timestep too large for SIDM)
           if(P_scatter > 1.0d0 .and. myid==1) then
              write(*,'(A,ES10.3,A,I2,A,I10)') &
                   ' WARNING: SIDM P=', P_scatter, &
                   ' > 1 at level ', ilevel, ' cell ', icell
           end if

           ! Monte Carlo: scatter if random < P
           call ranf(seed_loc, R1)
           if(R1 >= P_scatter) cycle

           ! === SCATTER EVENT ===
           n_scatter_loc = n_scatter_loc + 1

           ! Before: momentum and kinetic energy (physical)
           p_before(1:3) = m1*v1(1:3) + m2*v2(1:3)
           Ek_before = 0.5d0*(m1*(v1(1)**2+v1(2)**2+v1(3)**2) &
                             +m2*(v2(1)**2+v2(2)**2+v2(3)**2))

           ! Center-of-mass velocity (physical)
           v_cm(1:3) = (m1*v1(1:3) + m2*v2(1:3)) / mtot

           ! --- Inelastic energy change (multi-state iSIDM) ---
           delta_KE = 0.0d0
           do_transition = .false.
           new_state_1 = 0
           new_state_2 = 0
           if(sidm_inelastic) then
              ! Current states: tp=0 -> state 0, tp=-1 -> state 1, tp=-2 -> state 2, etc.
              istate_1 = nint(-tp(ind_dm(ipair)))
              istate_2 = nint(-tp(ind_dm(ipair+1)))
              istate_1 = max(0, min(istate_1, sidm_nstates-1))
              istate_2 = max(0, min(istate_2, sidm_nstates-1))
              Ei_1 = sidm_energy(istate_1) * 1.602d-9  ! keV -> erg
              Ei_2 = sidm_energy(istate_2) * 1.602d-9
              mu = m1*m2/mtot
              KE_cm = 0.5d0 * mu * v_rel_mag**2

              ! Enumerate all energetically accessible final channels
              n_channels = 0
              do j1=0,sidm_nstates-1
                 do j2=j1,sidm_nstates-1
                    ! Energy change: dE = (E_j1+E_j2) - (E_i1+E_i2)
                    ! dE>0 = endothermic (costs KE), dE<0 = exothermic (releases KE)
                    dE = (sidm_energy(j1)+sidm_energy(j2) &
                         -sidm_energy(istate_1)-sidm_energy(istate_2)) * 1.602d-9
                    if(j1==istate_1 .and. j2==istate_2) cycle  ! skip elastic
                    if(j1==istate_2 .and. j2==istate_1) cycle  ! skip trivial swap
                    if(dE > 0.0d0 .and. KE_cm < dE) cycle     ! not enough energy
                    n_channels = n_channels + 1
                    chan_j1(n_channels) = j1
                    chan_j2(n_channels) = j2
                    chan_dE(n_channels) = dE
                 end do
              end do

              if(n_channels > 0) then
                 call ranf(seed_loc, R_inel)
                 ichan = min(int(R_inel*dble(n_channels))+1, n_channels)
                 new_state_1 = chan_j1(ichan)
                 new_state_2 = chan_j2(ichan)
                 delta_KE = -chan_dE(ichan)  ! negative dE = KE gained
                 do_transition = .true.
                 ! Count up/down transitions
                 if(new_state_1 > istate_1) n_up_loc = n_up_loc + 1
                 if(new_state_1 < istate_1) n_down_loc = n_down_loc + 1
                 if(new_state_2 > istate_2) n_up_loc = n_up_loc + 1
                 if(new_state_2 < istate_2) n_down_loc = n_down_loc + 1
              end if
           end if

           ! New relative velocity magnitude
           mu = m1*m2/mtot
           if(delta_KE == 0.0d0) then
              v_rel_new_mag = v_rel_mag  ! elastic
           else
              KE_cm_new = 0.5d0*mu*v_rel_mag**2 + delta_KE
              if(KE_cm_new < 0.0d0) cycle  ! safety
              v_rel_new_mag = sqrt(2.0d0*KE_cm_new/mu)
           end if

           ! Dissipative SIDM: remove fraction f_diss of CM kinetic energy
           if(sidm_fdiss > 0.0d0) then
              KE_cm_before_diss = 0.5d0*mu*v_rel_new_mag**2
              v_rel_new_mag = v_rel_new_mag * sqrt(1.0d0 - sidm_fdiss)
              dEdiss_loc = dEdiss_loc + sidm_fdiss*KE_cm_before_diss
           end if

           ! --- Scattering angle ---
           select case(trim(sidm_angular))
           case('rutherford')
              ! Rutherford-like: dsigma/dOmega ~ 1/(1-cos_theta+2*eps)^2
              ! CDF inversion: cos_theta = 1+2*eps - 1/(b + R*(a-b))
              call ranf(seed_loc, R1)
              call ranf(seed_loc, R2)
              cos_theta = 1.0d0 + eps2 &
                   - 1.0d0/(b_ruth + R1*(a_ruth - b_ruth))
              cos_theta = max(-1.0d0, min(1.0d0, cos_theta))
              sin_theta = sqrt(max(0.0d0, 1.0d0 - cos_theta**2))
              phi_rand = twopi*R2

              ! Rotate scattering angle from v_rel frame to lab frame
              v_hat(1:3) = v_rel_vec(1:3) / v_rel_mag
              ! Orthonormal basis: e1 perp to v_hat
              if(abs(v_hat(3)) < 0.9d0) then
                 e1(1) =  v_hat(2)
                 e1(2) = -v_hat(1)
                 e1(3) =  0.0d0
              else
                 e1(1) =  0.0d0
                 e1(2) =  v_hat(3)
                 e1(3) = -v_hat(2)
              end if
              e1_mag = sqrt(e1(1)**2 + e1(2)**2 + e1(3)**2)
              e1(1:3) = e1(1:3) / e1_mag
              ! e2 = v_hat x e1
              e2(1) = v_hat(2)*e1(3) - v_hat(3)*e1(2)
              e2(2) = v_hat(3)*e1(1) - v_hat(1)*e1(3)
              e2(3) = v_hat(1)*e1(2) - v_hat(2)*e1(1)
              ! Scattered direction in lab frame
              nhat(1:3) = sin_theta*cos(phi_rand)*e1(1:3) &
                        + sin_theta*sin(phi_rand)*e2(1:3) &
                        + cos_theta*v_hat(1:3)

           case default  ! 'isotropic'
              ! Random unit vector (isotropic)
              call ranf(seed_loc, R1)
              call ranf(seed_loc, R2)
              cos_theta = 2.0d0*R1 - 1.0d0
              sin_theta = sqrt(1.0d0 - cos_theta**2)
              phi_rand  = twopi*R2
              nhat(1) = sin_theta*cos(phi_rand)
              nhat(2) = sin_theta*sin(phi_rand)
              nhat(3) = cos_theta
           end select

           ! New relative velocity vector
           v_rel_new(1:3) = v_rel_new_mag * nhat(1:3)

           ! New velocities in COM frame
           v1_new(1:3) = v_cm(1:3) + (m2/mtot)*v_rel_new(1:3)
           v2_new(1:3) = v_cm(1:3) - (m1/mtot)*v_rel_new(1:3)

           ! Apply state transitions (multi-state iSIDM)
           if(do_transition) then
              tp(ind_dm(ipair))   = -dble(new_state_1)
              tp(ind_dm(ipair+1)) = -dble(new_state_2)
           end if

           ! Conservation diagnostics (momentum: exact, energy: exact
           ! for elastic; intentional change for inelastic)
           p_after(1:3) = m1*v1_new(1:3) + m2*v2_new(1:3)
           Ek_after = 0.5d0*(m1*(v1_new(1)**2+v1_new(2)**2 &
                                +v1_new(3)**2) &
                            +m2*(v2_new(1)**2+v2_new(2)**2 &
                                +v2_new(3)**2))
           dp_loc(1:3) = dp_loc(1:3) + (p_after(1:3) - p_before(1:3))
           dEk_evt = Ek_after - Ek_before
           dEk_loc = dEk_loc + dEk_evt
           dp_mag = sqrt((p_after(1)-p_before(1))**2 &
                        +(p_after(2)-p_before(2))**2 &
                        +(p_after(3)-p_before(3))**2)
           if(dp_mag > dp_max_loc) dp_max_loc = dp_mag
           if(abs(dEk_evt) > dEk_max_loc) dEk_max_loc = abs(dEk_evt)

           ! Write new velocities (physical -> code units)
           vp(ind_dm(ipair),   1:3) = v1_new(1:3) * aexp / scale_v
           vp(ind_dm(ipair+1), 1:3) = v2_new(1:3) * aexp / scale_v
        end do  ! pairs

     end do  ! subcells (ind)

     igrid = next(igrid)
  end do  ! grids

  ! Save seed back to per-thread array (no race)
  thread_seeds(:,mythread_loc) = seed_loc

  ! Accumulate thread-local counters to shared counters
  !$omp critical
  sidm_n_scatter = sidm_n_scatter + n_scatter_loc
  sidm_n_pairs   = sidm_n_pairs + n_pairs_loc
  sidm_n_up      = sidm_n_up + n_up_loc
  sidm_n_down    = sidm_n_down + n_down_loc
  sidm_dp(1:3)   = sidm_dp(1:3) + dp_loc(1:3)
  sidm_dEk       = sidm_dEk + dEk_loc
  if(dp_max_loc > sidm_dp_max) sidm_dp_max = dp_max_loc
  if(dEk_max_loc > sidm_dEk_max) sidm_dEk_max = dEk_max_loc
  if(P_max_loc > sidm_Pmax_level) sidm_Pmax_level = P_max_loc
  sidm_dEdiss = sidm_dEdiss + dEdiss_loc
  !$omp end critical

end subroutine sub_sidm_scatter
!################################################################
!################################################################
! Initialize excited DM fraction for iSIDM
!################################################################
subroutine sidm_init_excited()
  use pm_commons
  use pm_parameters, only: iseed
  use amr_commons
  use random
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer::ilevel,icpu,igrid,jgrid,ipart,jpart,npart1,is
  integer::ndm_total,ndm_existing,info
  integer,dimension(0:9)::ndm_state
  integer,dimension(1:ncpu,1:IRandNumSize)::allseed
  integer,dimension(1:IRandNumSize)::seed_init
  real(dp)::R1,cum_frac
  real(dp),dimension(0:9)::cum_prob  ! Cumulative probability for N states

  ! Check if any non-ground population is requested
  if(sum(sidm_frac_init(1:sidm_nstates-1)) <= 0.0d0) return

  ! Build cumulative probability
  cum_prob(0) = sidm_frac_init(0)
  do is=1,sidm_nstates-1
     cum_prob(is) = cum_prob(is-1) + sidm_frac_init(is)
  end do
  ! Normalize (in case fractions don't sum to 1)
  if(cum_prob(sidm_nstates-1) > 0.0d0) then
     cum_prob = cum_prob / cum_prob(sidm_nstates-1)
  end if

  ! Check if already initialized (restart case: some DM have tp<0)
  ndm_existing = 0
  do ilevel=levelmin,nlevelmax
     igrid = headl(myid,ilevel)
     do jgrid=1,numbl(myid,ilevel)
        ipart = headp(igrid)
        do jpart=1,numbp(igrid)
           if(idp(ipart)>0 .and. tp(ipart)<-0.5d0) then
              ndm_existing = ndm_existing + 1
           end if
           ipart = nextp(ipart)
        end do
        igrid = next(igrid)
     end do
  end do
  call MPI_ALLREDUCE(MPI_IN_PLACE,ndm_existing,1,MPI_INTEGER, &
       MPI_SUM,MPI_COMM_WORLD,info)
  if(ndm_existing > 0) then
     if(myid==1) write(*,'(A,I10,A)') &
          ' iSIDM: found ', ndm_existing, &
          ' pre-existing excited DM (skipping init)'
     return
  end if

  ! Initialize seed (offset to avoid correlation with main RNG)
  call rans(ncpu,iseed+12345,allseed)
  seed_init = allseed(myid,1:IRandNumSize)

  ndm_total = 0
  ndm_state = 0

  ! Loop over own grids and randomly assign states
  do ilevel=levelmin,nlevelmax
     igrid = headl(myid,ilevel)
     do jgrid=1,numbl(myid,ilevel)
        ipart = headp(igrid)
        do jpart=1,numbp(igrid)
           ! DM: idp>0, tp==0 (ground state, not yet assigned)
           if(idp(ipart)>0 .and. tp(ipart)==0.0d0) then
              ndm_total = ndm_total + 1
              call ranf(seed_init, R1)
              ! Assign state by cumulative probability
              do is=0,sidm_nstates-1
                 if(R1 < cum_prob(is)) then
                    tp(ipart) = -dble(is)  ! state 0 -> tp=0, state 1 -> tp=-1, etc.
                    ndm_state(is) = ndm_state(is) + 1
                    exit
                 end if
              end do
           end if
           ipart = nextp(ipart)
        end do
        igrid = next(igrid)
     end do
  end do

  ! Report
  call MPI_ALLREDUCE(MPI_IN_PLACE,ndm_total,1,MPI_INTEGER, &
       MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(MPI_IN_PLACE,ndm_state,sidm_nstates,MPI_INTEGER, &
       MPI_SUM,MPI_COMM_WORLD,info)
  if(myid==1) then
     write(*,'(A,I10,A,I3,A)') &
          ' iSIDM: initialized ', ndm_total, ' DM into ', sidm_nstates, ' states:'
     do is=0,sidm_nstates-1
        write(*,'(A,I2,A,I10,A,F7.4)') &
             '   state ', is, ': N=', ndm_state(is), &
             '  frac=', dble(ndm_state(is))/dble(max(1,ndm_total))
     end do
  end if

end subroutine sidm_init_excited
!################################################################
!################################################################
! Helper: determine subcell index (1..8) for a particle in a grid
!################################################################
integer function cell_index_from_part(ipart,igrid,ilevel)
  use pm_commons, only: xp
  use amr_commons, only: xg,dp
  implicit none
  integer,intent(in)::ipart,igrid,ilevel
  integer::ix,iy,iz

  ! Subcell indices (0 or 1 in each dimension)
  ix = 0; iy = 0; iz = 0
  if(xp(ipart,1) > xg(igrid,1)) ix = 1
  if(xp(ipart,2) > xg(igrid,2)) iy = 1
#if NDIM==3
  if(xp(ipart,3) > xg(igrid,3)) iz = 1
#endif

  ! RAMSES subcell ordering: ind = 1 + ix + 2*iy + 4*iz
  cell_index_from_part = 1 + ix + 2*iy + 4*iz

end function cell_index_from_part
!################################################################
!################################################################
! Report excited DM fraction (called from adaptive_loop)
!################################################################
subroutine sidm_report_excited_fraction()
  use pm_commons
  use amr_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer::ilevel,igrid,jgrid,ipart,jpart,info,is,istate
  integer::ndm_total
  integer,dimension(0:9)::ndm_state

  if(.not.sidm_inelastic) return

  ndm_total = 0
  ndm_state = 0
  do ilevel=levelmin,nlevelmax
     igrid = headl(myid,ilevel)
     do jgrid=1,numbl(myid,ilevel)
        ipart = headp(igrid)
        do jpart=1,numbp(igrid)
           if(idp(ipart)>0 .and. tp(ipart)<=0.0d0) then
              ndm_total = ndm_total + 1
              istate = nint(-tp(ipart))
              istate = max(0, min(istate, 9))
              ndm_state(istate) = ndm_state(istate) + 1
           end if
           ipart = nextp(ipart)
        end do
        igrid = next(igrid)
     end do
  end do
  call MPI_ALLREDUCE(MPI_IN_PLACE,ndm_total,1,MPI_INTEGER, &
       MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(MPI_IN_PLACE,ndm_state,sidm_nstates,MPI_INTEGER, &
       MPI_SUM,MPI_COMM_WORLD,info)

  if(myid==1 .and. ndm_total>0) then
     if(sidm_nstates==2) then
        ! Compact format for 2-state (backward compatible output)
        write(*,'(A,I10,A,I10,A,F8.5)') &
             ' iSIDM excited: ', ndm_state(1), ' / ', ndm_total, &
             '  f_exc=', dble(ndm_state(1))/dble(ndm_total)
     else
        write(*,'(A,I10,A,I3,A)') &
             ' iSIDM populations: N_DM=', ndm_total, '  (', sidm_nstates, ' states)'
        do is=0,sidm_nstates-1
           write(*,'(A,I2,A,I10,A,F7.4)') &
                '   state ', is, ': ', ndm_state(is), &
                '  f=', dble(ndm_state(is))/dble(ndm_total)
        end do
     end if
  end if

end subroutine sidm_report_excited_fraction
!################################################################
!################################################################
! DM-baryon drag force: continuous momentum exchange
!
! Physics: F/m_DM = -(sigma_db/m_DM) * rho_b * |v_rel| * v_rel
! Velocity-dependent: sigma(v) = sigma_0 * (v/100 km/s)^n
! Backreaction on gas conserves total momentum.
!################################################################
subroutine sidm_baryon_drag(ilevel)
  use pm_commons
  use amr_commons
  use hydro_commons, only: uold,nvar,gamma
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer,intent(in)::ilevel

  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  real(dp)::dx,dx_loc,scale,vol_phys,vol_code,dt_phys
  real(dp)::rho_b_phys,sigma_eff
  real(dp),dimension(1:3)::v_gas_phys,v_dm_phys,v_rel,dv_dm
  real(dp)::v_rel_mag,v_rel_kms,drag_coeff
  real(dp)::m_dm_phys,dv_max,dv_mag
  integer::igrid,jgrid,ind,iskip,icell,nx_loc
  integer::ipart,jpart,npart1,info
  ! Diagnostics
  integer::n_drag_loc,n_drag_tot
  real(dp)::dp_drag_loc(3),dp_drag_tot(3),dv_max_loc,dv_max_tot

  if(.not.sidm_baryon) return
  if(.not.hydro) return
  if(numbtot(1,ilevel)==0) return

  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

  dx = 0.5D0**ilevel
  nx_loc = (icoarse_max-icoarse_min+1)
  scale = boxlen/dble(nx_loc)
  dx_loc = dx*scale
  vol_code = dx_loc**3
  vol_phys = (dx_loc*scale_l/aexp)**3
  dt_phys  = dtnew(ilevel)*scale_t

  n_drag_loc  = 0
  dp_drag_loc = 0.0d0
  dv_max_loc  = 0.0d0

  ! Loop over local grids
  igrid = headl(myid,ilevel)
  do jgrid=1,numbl(myid,ilevel)

     npart1 = numbp(igrid)
     if(npart1==0) then
        igrid = next(igrid)
        cycle
     end if

     ! Loop over 8 subcells
     do ind=1,twotondim
        iskip = ncoarse+(ind-1)*ngridmax
        icell = iskip+igrid

        ! Only leaf cells with gas
        if(son(icell)/=0) cycle

        ! Gas density and velocity (code units)
        rho_b_phys = uold(icell,1) * scale_d / aexp**3
        if(uold(icell,1) <= 0.0d0) cycle

        v_gas_phys(1) = uold(icell,2)/uold(icell,1) * scale_v / aexp
        v_gas_phys(2) = uold(icell,3)/uold(icell,1) * scale_v / aexp
        v_gas_phys(3) = uold(icell,4)/uold(icell,1) * scale_v / aexp

        ! Loop over particles in this grid
        ipart = headp(igrid)
        do jpart=1,npart1
           ! DM only (idp>0, tp<=0 for non-stars)
           if(idp(ipart)>0 .and. tp(ipart)<=0.0d0) then
              ! Check subcell match
              if(cell_index_from_part_inline(ipart,igrid)==ind) then
                 ! DM velocity (physical)
                 v_dm_phys(1) = vp(ipart,1) * scale_v / aexp
                 v_dm_phys(2) = vp(ipart,2) * scale_v / aexp
                 v_dm_phys(3) = vp(ipart,3) * scale_v / aexp

                 ! Relative velocity
                 v_rel(1:3) = v_dm_phys(1:3) - v_gas_phys(1:3)
                 v_rel_mag = sqrt(v_rel(1)**2+v_rel(2)**2+v_rel(3)**2)
                 if(v_rel_mag == 0.0d0) then
                    ipart = nextp(ipart)
                    cycle
                 end if

                 ! Velocity-dependent cross-section
                 v_rel_kms = v_rel_mag * 1.0d-5
                 if(sidm_baryon_power == 0.0d0) then
                    sigma_eff = sidm_baryon_sigma
                 else
                    sigma_eff = sidm_baryon_sigma &
                         * (v_rel_kms/100.0d0)**sidm_baryon_power
                    sigma_eff = min(sigma_eff, 1.0d4*sidm_baryon_sigma)
                 end if

                 ! Drag: dv = -(sigma/m) * rho_b * |v_rel| * v_rel * dt
                 drag_coeff = sigma_eff * rho_b_phys * v_rel_mag * dt_phys
                 dv_dm(1:3) = -drag_coeff * v_rel(1:3)

                 ! Limit kick to avoid overshoot (cap at 50% of v_rel)
                 dv_mag = sqrt(dv_dm(1)**2+dv_dm(2)**2+dv_dm(3)**2)
                 dv_max = 0.5d0 * v_rel_mag
                 if(dv_mag > dv_max) then
                    dv_dm(1:3) = dv_dm(1:3) * (dv_max/dv_mag)
                    dv_mag = dv_max
                 end if

                 ! Apply to DM particle (physical -> code)
                 vp(ipart,1) = vp(ipart,1) + dv_dm(1)*aexp/scale_v
                 vp(ipart,2) = vp(ipart,2) + dv_dm(2)*aexp/scale_v
                 vp(ipart,3) = vp(ipart,3) + dv_dm(3)*aexp/scale_v

                 ! Backreaction on gas: dp_gas = -m_DM * dv_DM
                 ! Update momentum density: d(rho*v) = m_DM * (-dv_DM) / V_cell
                 m_dm_phys = mp(ipart) * scale_d * scale_l**3
                 uold(icell,2) = uold(icell,2) &
                      - m_dm_phys*dv_dm(1) / (vol_phys*scale_d/aexp**3) &
                      * (aexp/scale_v)
                 uold(icell,3) = uold(icell,3) &
                      - m_dm_phys*dv_dm(2) / (vol_phys*scale_d/aexp**3) &
                      * (aexp/scale_v)
                 uold(icell,4) = uold(icell,4) &
                      - m_dm_phys*dv_dm(3) / (vol_phys*scale_d/aexp**3) &
                      * (aexp/scale_v)

                 ! Diagnostics
                 n_drag_loc = n_drag_loc + 1
                 dp_drag_loc(1:3) = dp_drag_loc(1:3) + m_dm_phys*dv_dm(1:3)
                 if(dv_mag > dv_max_loc) dv_max_loc = dv_mag
              end if
           end if
           ipart = nextp(ipart)
        end do  ! particles

     end do  ! subcells

     igrid = next(igrid)
  end do  ! grids

  ! MPI reduce diagnostics
  call MPI_ALLREDUCE(n_drag_loc, n_drag_tot, 1, MPI_INTEGER, &
       MPI_SUM, MPI_COMM_WORLD, info)
  call MPI_ALLREDUCE(dp_drag_loc, dp_drag_tot, 3, MPI_DOUBLE_PRECISION, &
       MPI_SUM, MPI_COMM_WORLD, info)
  call MPI_ALLREDUCE(dv_max_loc, dv_max_tot, 1, MPI_DOUBLE_PRECISION, &
       MPI_MAX, MPI_COMM_WORLD, info)

  if(myid==1 .and. n_drag_tot>0) then
     write(*,'(A,I2,A,I10,A,ES9.2,A)') &
          ' IDM level ',ilevel,': n_drag=',n_drag_tot, &
          ' dv_max=',dv_max_tot*1d-5,' km/s'
  end if

contains
  ! Inline subcell index (avoids external function in tight loop)
  integer function cell_index_from_part_inline(ip,ig)
    integer,intent(in)::ip,ig
    integer::ixx,iyy,izz
    ixx=0; iyy=0; izz=0
    if(xp(ip,1) > xg(ig,1)) ixx=1
    if(xp(ip,2) > xg(ig,2)) iyy=1
    if(xp(ip,3) > xg(ig,3)) izz=1
    cell_index_from_part_inline = 1 + ixx + 2*iyy + 4*izz
  end function cell_index_from_part_inline

end subroutine sidm_baryon_drag
