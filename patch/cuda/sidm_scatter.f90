!################################################################
!################################################################
!################################################################
! SIDM (Self-Interacting Dark Matter) Monte Carlo Scattering
!
! Cell-based random pairing with isotropic elastic scattering.
! Constant cross-section (sigma/m) in cm^2/g.
!################################################################
subroutine sidm_scatter(ilevel)
  use pm_commons
  use amr_commons
  implicit none
  integer,intent(in)::ilevel

  integer::icpu,igrid,subnump
  integer::mythread,nthreads
  integer,dimension(:),allocatable::nparticles,ptrhead
  common /sidm_omp/ mythread
!$omp threadprivate(/sidm_omp/)

  if(.not.sidm) return
  if(numbtot(1,ilevel)==0) return
  if(verbose) write(*,111) ilevel

!$omp parallel
  mythread = omp_get_thread_num()
  if(mythread==0) nthreads = omp_get_num_threads()
!$omp end parallel
  allocate(ptrhead(0:nthreads-1), nparticles(0:nthreads-1))

#if NDIM==3
  do icpu=1,ncpu
     if(numbl(icpu,ilevel)<=0) cycle
     call pthreadLinkedList(headl(icpu,ilevel),numbl(icpu,ilevel), &
          nthreads,nparticles,ptrhead,next)
!$omp parallel private(subnump,igrid)
     subnump = nparticles(mythread)
     igrid   = ptrhead(mythread)
     call sub_sidm_scatter(ilevel,icpu,igrid,subnump)
!$omp end parallel
  end do
#endif

  deallocate(ptrhead, nparticles)

111 format('   Entering sidm_scatter for level ',I2)
end subroutine sidm_scatter
!################################################################
!################################################################
subroutine sub_sidm_scatter(ilevel,icpu,kgrid,subnump)
  use pm_commons
  use pm_parameters, only: iseed
  use amr_commons
  use random
  implicit none
  integer,intent(in)::ilevel,icpu,kgrid,subnump

  ! External function
  integer,external::cell_index_from_part

  ! Local variables
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  real(dp)::dx,dx_loc,scale,vol_phys,dt_phys,sigma_over_m
  real(dp)::mp_phys,P_scatter,twopi
  real(dp)::cos_theta,sin_theta,phi_rand,v_rel_mag
  real(dp),dimension(1:3)::v1,v2,v_cm,v_rel_vec,nhat,v_rel_new
  real(dp)::m1,m2,mtot,R1,R2
  integer::igrid,jgrid,ind,iskip,icell,nx_loc
  integer::ipart,jpart,npart1,ndm_cell,ip,jp,npairs
  integer,dimension(1:nvector)::ind_dm
  integer::itemp,ipair
  integer,dimension(IRandNumSize)::seed_loc
  integer,dimension(1:ncpu,1:IRandNumSize)::allseed

  ! Scattering counters (thread-local)
  integer::n_scatter_loc,n_pairs_loc

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

  ! sigma/m in CGS [cm^2/g]
  sigma_over_m = sidm_cross_section

  ! Initialize thread-local random seed from global localseed
  if(localseed(1)==-1) then
     call rans(ncpu,iseed,allseed)
     localseed = allseed(myid,1:IRandNumSize)
  end if
  seed_loc = localseed

  n_scatter_loc = 0
  n_pairs_loc = 0

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
        ndm_cell = 0
        ipart = headp(igrid)
        do jp=1,npart1
           ! Check if this particle belongs to subcell 'ind'
           if(cell_index_from_part(ipart,igrid,ilevel)==ind) then
              ! DM particle: idp>0 and tp==0 (not star, not sink)
              if(idp(ipart)>0 .and. tp(ipart)==0.0d0) then
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
           v_rel_mag = sqrt(v_rel_vec(1)**2 + v_rel_vec(2)**2 + v_rel_vec(3)**2)

           if(v_rel_mag==0.0d0) cycle

           ! Scattering probability
           ! P = (sigma/m) * m_p * v_rel * dt / V_cell * (N_dm - 1)
           mp_phys = 0.5d0*(m1+m2)  ! representative particle mass
           P_scatter = sigma_over_m * mp_phys * v_rel_mag * dt_phys &
                     / vol_phys * dble(ndm_cell-1)

           ! Warning if P > 1 (timestep too large for SIDM)
           if(P_scatter > 1.0d0 .and. myid==1) then
              write(*,'(A,ES10.3,A,I2,A,I10)') &
                   ' WARNING: SIDM P=', P_scatter, ' > 1 at level ', ilevel, &
                   ' cell ', icell
           end if

           ! Monte Carlo: scatter if random < P
           call ranf(seed_loc, R1)
           if(R1 >= P_scatter) cycle

           ! --- Isotropic elastic scattering ---
           n_scatter_loc = n_scatter_loc + 1

           ! Center-of-mass velocity (physical)
           v_cm(1:3) = (m1*v1(1:3) + m2*v2(1:3)) / mtot

           ! Random unit vector (isotropic)
           call ranf(seed_loc, R1)
           call ranf(seed_loc, R2)
           cos_theta = 2.0d0*R1 - 1.0d0
           sin_theta = sqrt(1.0d0 - cos_theta**2)
           phi_rand  = twopi*R2
           nhat(1) = sin_theta*cos(phi_rand)
           nhat(2) = sin_theta*sin(phi_rand)
           nhat(3) = cos_theta

           ! New relative velocity: same magnitude, random direction
           v_rel_new(1:3) = v_rel_mag * nhat(1:3)

           ! New velocities in COM frame (physical -> code units)
           ! v1_new = v_cm + (m2/mtot)*v_rel_new
           ! v2_new = v_cm - (m1/mtot)*v_rel_new
           vp(ind_dm(ipair),   1:3) = (v_cm(1:3) + (m2/mtot)*v_rel_new(1:3)) &
                                      * aexp / scale_v
           vp(ind_dm(ipair+1), 1:3) = (v_cm(1:3) - (m1/mtot)*v_rel_new(1:3)) &
                                      * aexp / scale_v
        end do  ! pairs

     end do  ! subcells (ind)

     igrid = next(igrid)
  end do  ! grids

  ! Update global seed
  !$omp critical
  localseed = seed_loc
  !$omp end critical

end subroutine sub_sidm_scatter
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
