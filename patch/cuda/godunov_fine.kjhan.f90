#ifdef HYDRO_CUDA
!###########################################################
! Module for hybrid CPU/GPU hydro dispatch data structures
!###########################################################
module hydro_hybrid_commons
  use amr_parameters, only: dp
  implicit none

  integer, parameter :: HYBRID_SUPER_SIZE = 4096
  integer, parameter :: SBUF_INIT_CAP = 65536

  type scatter_buf_t
     integer :: count = 0
     integer :: capacity = 0
     integer, allocatable :: icell(:)
     real(dp), allocatable :: dunew(:,:)
     real(dp), allocatable :: ddivu(:)
     real(dp), allocatable :: denew(:)
  end type scatter_buf_t

  type gpu_state_t
     real(dp), allocatable :: super_uloc(:,:,:,:,:)
     real(dp), allocatable :: super_gloc(:,:,:,:,:)
     real(dp), allocatable :: super_flux(:,:,:,:,:,:)
     real(dp), allocatable :: super_tmp(:,:,:,:,:,:)
     logical,  allocatable :: super_ok(:,:,:,:)
     integer,  allocatable :: super_ind_grid(:)
     integer :: off = 0
  end type gpu_state_t

  type(scatter_buf_t), allocatable, save :: scatter_bufs(:)
  type(gpu_state_t), allocatable, save :: gpu_states(:)
  integer, save :: hybrid_max_threads = 0

  ! Timing accumulators
  real(dp), save :: acc_gather = 0, acc_gpu = 0, acc_scatter_l = 0, acc_merge = 0
  integer, save :: n_gpu_flushes = 0, n_cpu_batches = 0

contains

  subroutine scatter_buf_init(buf, cap, nv)
    type(scatter_buf_t), intent(inout) :: buf
    integer, intent(in) :: cap, nv
    buf%count = 0
    buf%capacity = cap
    allocate(buf%icell(cap))
    allocate(buf%dunew(cap, nv))
    allocate(buf%ddivu(cap))
    allocate(buf%denew(cap))
  end subroutine scatter_buf_init

  subroutine scatter_buf_grow(buf)
    type(scatter_buf_t), intent(inout) :: buf
    integer :: new_cap, nv
    integer, allocatable :: tmp_icell(:)
    real(dp), allocatable :: tmp_dunew(:,:), tmp_ddivu(:), tmp_denew(:)

    new_cap = buf%capacity * 2
    nv = size(buf%dunew, 2)

    allocate(tmp_icell(new_cap))
    allocate(tmp_dunew(new_cap, nv))
    allocate(tmp_ddivu(new_cap))
    allocate(tmp_denew(new_cap))

    tmp_icell(1:buf%capacity) = buf%icell(1:buf%capacity)
    tmp_dunew(1:buf%capacity,:) = buf%dunew(1:buf%capacity,:)
    tmp_ddivu(1:buf%capacity) = buf%ddivu(1:buf%capacity)
    tmp_denew(1:buf%capacity) = buf%denew(1:buf%capacity)

    call move_alloc(tmp_icell, buf%icell)
    call move_alloc(tmp_dunew, buf%dunew)
    call move_alloc(tmp_ddivu, buf%ddivu)
    call move_alloc(tmp_denew, buf%denew)

    buf%capacity = new_cap
  end subroutine scatter_buf_grow

end module hydro_hybrid_commons
#endif
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine godunov_fine(ilevel)
  use amr_commons
  use hydro_commons
#ifdef HYDRO_CUDA
  use cuda_commons
  use hydro_cuda_interface
#endif
  implicit none
  integer::ilevel
  !--------------------------------------------------------------------------
  ! This routine is a wrapper to the second order Godunov solver.
  ! Small grids (2x2x2) are gathered from level ilevel and sent to the
  ! hydro solver. On entry, hydro variables are gathered from array uold.
  ! On exit, unew has been updated.
  !
  ! HYDRO_CUDA: OMP threads dynamically dispatch batches to GPU streams
  ! or CPU fallback depending on stream availability.
  !--------------------------------------------------------------------------
  integer::i,ivar,igrid,ncache,ngrid
#ifdef HYDRO_CUDA
  logical, save :: cuda_init_done = .false.
#endif

  if(numbtot(1,ilevel)==0)return
  if(static)return
  if(verbose)write(*,111)ilevel

#ifdef HYDRO_CUDA
  ! One-time CUDA initialization (after hydro parameters are read)
  if (.not. cuda_init_done) then
     !$omp single
     call cuda_pool_init_f()
     if (cuda_available) then
        call hydro_cuda_init_f()
        if(myid==1) write(*,'(A,I2,A)') &
             ' CUDA hydro: initialized with ', cuda_n_streams, ' GPU streams'
     else
        if(myid==1) write(*,'(A)') &
             ' CUDA hydro: no GPU available, using CPU-only mode'
     end if
     cuda_init_done = .true.
     !$omp end single
  end if
#endif

  ! Loop over active grids by vector sweeps
  ncache=active(ilevel)%ngrid

#ifdef HYDRO_CUDA
  if (cuda_available .and. hydro_cuda_initialized) then
     call godunov_fine_hybrid(ilevel, ncache)
  else
#endif
#ifndef OMP_TEST
!$omp parallel do private(igrid,ngrid) shared(ncache, ilevel)
#endif
  do igrid=1,ncache,nvector
     ngrid=MIN(nvector,ncache-igrid+1)
     call godfine1(ilevel, igrid,ngrid)
  end do
#ifdef HYDRO_CUDA
  end if
#endif

111 format('  +Entering godunov_fine for level ',i2)
end subroutine godunov_fine
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine set_unew(ilevel)
  use amr_commons
  use hydro_commons
  implicit none
  integer::ilevel
  !--------------------------------------------------------------------------
  ! This routine sets array unew to its initial value uold before calling
  ! the hydro scheme. unew is set to zero in virtual boundaries.
  !--------------------------------------------------------------------------
  integer::i,ivar,irad,ind,icpu,iskip
  real(dp)::d,u,v,w,e

  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel

  ! Set unew to uold for myid cells
!$omp parallel do private(ind,iskip,ivar,i,d,u,v,w,e,irad) schedule(dynamic)
  do ind=1,twotondim
     iskip=ncoarse+(ind-1)*ngridmax
     do ivar=1,nvar
        do i=1,active(ilevel)%ngrid
           unew(active(ilevel)%igrid(i)+iskip,ivar) = uold(active(ilevel)%igrid(i)+iskip,ivar)
        end do
     end do
     if(pressure_fix)then
        do i=1,active(ilevel)%ngrid
           divu(active(ilevel)%igrid(i)+iskip) = 0.0
        end do
        do i=1,active(ilevel)%ngrid
           d=max(uold(active(ilevel)%igrid(i)+iskip,1),smallr)
           u=0.0; v=0.0; w=0.0
           if(ndim>0)u=uold(active(ilevel)%igrid(i)+iskip,2)/d
           if(ndim>1)v=uold(active(ilevel)%igrid(i)+iskip,3)/d
           if(ndim>2)w=uold(active(ilevel)%igrid(i)+iskip,4)/d
           e=uold(active(ilevel)%igrid(i)+iskip,ndim+2)-0.5*d*(u**2+v**2+w**2)
#if NENER>0
           do irad=1,nener
              e=e-uold(active(ilevel)%igrid(i)+iskip,ndim+2+irad)
           end do
#endif
           enew(active(ilevel)%igrid(i)+iskip)=e
        end do
     end if
  end do

  ! Set unew to 0 for virtual boundary cells
  do icpu=1,ncpu
!$omp parallel do private(ind,iskip,ivar,i) schedule(dynamic)
  do ind=1,twotondim
     iskip=ncoarse+(ind-1)*ngridmax
     do ivar=1,nvar
        do i=1,reception(icpu,ilevel)%ngrid
           unew(reception(icpu,ilevel)%igrid(i)+iskip,ivar)=0.0
        end do
     end do
     if(pressure_fix)then
        do i=1,reception(icpu,ilevel)%ngrid
           divu(reception(icpu,ilevel)%igrid(i)+iskip) = 0.0
           enew(reception(icpu,ilevel)%igrid(i)+iskip) = 0.0
        end do
     end if
  end do
  end do

111 format('  +Entering set_unew for level ',i2)

end subroutine set_unew
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine set_uold(ilevel)
  use amr_commons
  use hydro_commons
  use poisson_commons
  implicit none
  integer::ilevel
  !---------------------------------------------------------
  ! This routine sets array uold to its new value unew
  ! after the hydro step.
  !---------------------------------------------------------
  integer::i,ivar,irad,ind,iskip,nx_loc,ind_cell
  real(dp)::scale,d,u,v,w
  real(dp)::e_kin,e_cons,e_prim,e_trunc,div,dx,fact,d_old

  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel

  nx_loc=icoarse_max-icoarse_min+1
  scale=boxlen/dble(nx_loc)
  dx=0.5d0**ilevel*scale

  ! Add gravity source terms to unew
  if(poisson)then
     call add_gravity_source_terms(ilevel)
  end if

  ! Add non conservative pdV terms to unew
  ! for thermal and/or non-thermal energies
  if(pressure_fix.OR.nener>0)then
     call add_pdv_source_terms(ilevel)
  endif

  ! Set uold to unew for myid cells
!$omp parallel do private(ind,iskip,ivar,i,ind_cell,d,u,v,w,e_kin,e_cons,e_prim,div,e_trunc)
  do ind=1,twotondim
     iskip=ncoarse+(ind-1)*ngridmax
     do ivar=1,nvar
        do i=1,active(ilevel)%ngrid
           uold(active(ilevel)%igrid(i)+iskip,ivar) = unew(active(ilevel)%igrid(i)+iskip,ivar)
        end do
     end do
     if(pressure_fix)then
        ! Correct total energy if internal energy is too small
        do i=1,active(ilevel)%ngrid
           ind_cell=active(ilevel)%igrid(i)+iskip
           d=max(uold(ind_cell,1),smallr)
           u=0.0; v=0.0; w=0.0
           if(ndim>0)u=uold(ind_cell,2)/d
           if(ndim>1)v=uold(ind_cell,3)/d
           if(ndim>2)w=uold(ind_cell,4)/d
           e_kin=0.5*d*(u**2+v**2+w**2)
#if NENER>0
           do irad=1,nener
              e_kin=e_kin+uold(ind_cell,ndim+2+irad)
           end do
#endif
           e_cons=uold(ind_cell,ndim+2)-e_kin
           e_prim=enew(ind_cell)
           ! Note: here divu=-div.u*dt
           div=abs(divu(ind_cell))*dx/dtnew(ilevel)
           e_trunc=beta_fix*d*max(div,3.0*hexp*dx)**2
           if(e_cons<e_trunc)then
              uold(ind_cell,ndim+2)=e_prim+e_kin
           end if
        end do
     end if
  end do

111 format('  +Entering set_uold for level ',i2)

end subroutine set_uold
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine add_gravity_source_terms(ilevel)
  use amr_commons
  use hydro_commons
  use poisson_commons
  implicit none
  integer::ilevel
  !--------------------------------------------------------------------------
  ! This routine adds to unew the gravity source terms
  ! with only half a time step. Only the momentum and the
  ! total energy are modified in array unew.
  !--------------------------------------------------------------------------
  integer::i,ivar,ind,iskip,nx_loc,ind_cell
  real(dp)::d,u,v,w,e_kin,e_prim,d_old,fact

  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel

  ! Add gravity source term at time t with half time step
! do ind=1,twotondim
!    iskip=ncoarse+(ind-1)*ngridmax
!    do i=1,active(ilevel)%ngrid
!$omp parallel do private(ind,iskip,i,ind_cell,d,u,v,w,e_kin,d_old,e_prim,fact)
  do i=1,active(ilevel)%ngrid
     do ind=1,twotondim
        iskip=ncoarse+(ind-1)*ngridmax
        ind_cell=active(ilevel)%igrid(i)+iskip
        d=max(unew(ind_cell,1),smallr)
        u=0.0; v=0.0; w=0.0
        if(ndim>0)u=unew(ind_cell,2)/d
        if(ndim>1)v=unew(ind_cell,3)/d
        if(ndim>2)w=unew(ind_cell,4)/d
        e_kin=0.5*d*(u**2+v**2+w**2)
        e_prim=unew(ind_cell,ndim+2)-e_kin
        d_old=max(uold(ind_cell,1),smallr)
        fact=d_old/d*0.5*dtnew(ilevel)
        if(ndim>0)then
           u=u+f(ind_cell,1)*fact
           unew(ind_cell,2)=d*u
        endif
        if(ndim>1)then
           v=v+f(ind_cell,2)*fact
           unew(ind_cell,3)=d*v
        end if
        if(ndim>2)then
           w=w+f(ind_cell,3)*fact
           unew(ind_cell,4)=d*w
        endif
        e_kin=0.5*d*(u**2+v**2+w**2)
        unew(ind_cell,ndim+2)=e_prim+e_kin
     end do
  end do

111 format('  +Entering add_gravity_source_terms for level ',i2)

end subroutine add_gravity_source_terms
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine add_pdv_source_terms(ilevel)
  use amr_commons
  use hydro_commons
  implicit none
  integer:: ncache, igrid,ngrid,ilevel

  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel

  ncache = active(ilevel)%ngrid
!$omp parallel do private(igrid,ngrid) schedule(dynamic)
  do igrid=1, ncache, nvector
     ngrid = MIN(nvector,ncache-igrid+1)
     call sub_add_pdv_source_terms(ilevel,igrid,ngrid)
  enddo
111 format('  +Entering add_pdv_source_terms for level ',i2)
end subroutine add_pdv_source_terms

subroutine sub_add_pdv_source_terms(ilevel,igrid,ngrid)
  use amr_commons
  use hydro_commons
  use morton_hash
  implicit none
  integer::ilevel
  !---------------------------------------------------------
  ! This routine adds the pdV source term to the internal
  ! energy equation and to the non-thermal energy equations.
  !---------------------------------------------------------
  integer::i,ivar,irad,ind,iskip,nx_loc,ind_cell1
  integer::ncache,igrid,ngrid,idim,id1,ig1,ih1,id2,ig2,ih2
  integer,dimension(1:3,1:2,1:8)::iii,jjj
  real(dp)::scale,dx,dx_loc,d,u,v,w,eold

  integer ,dimension(1:nvector)::ind_grid,ind_cell
  integer ,dimension(1:nvector,0:twondim)::igridn
  integer ,dimension(1:nvector,1:ndim)::ind_left,ind_right
  real(dp),dimension(1:nvector,1:ndim,1:ndim)::velg,veld
  real(dp),dimension(1:nvector,1:ndim)::dx_g,dx_d
  real(dp),dimension(1:nvector)::divu_loc


  nx_loc=icoarse_max-icoarse_min+1
  scale=boxlen/dble(nx_loc)
  dx=0.5d0**ilevel
  dx_loc=dx*scale

  iii(1,1,1:8)=(/1,0,1,0,1,0,1,0/); jjj(1,1,1:8)=(/2,1,4,3,6,5,8,7/)
  iii(1,2,1:8)=(/0,2,0,2,0,2,0,2/); jjj(1,2,1:8)=(/2,1,4,3,6,5,8,7/)
  iii(2,1,1:8)=(/3,3,0,0,3,3,0,0/); jjj(2,1,1:8)=(/3,4,1,2,7,8,5,6/)
  iii(2,2,1:8)=(/0,0,4,4,0,0,4,4/); jjj(2,2,1:8)=(/3,4,1,2,7,8,5,6/)
  iii(3,1,1:8)=(/5,5,5,5,0,0,0,0/); jjj(3,1,1:8)=(/5,6,7,8,1,2,3,4/)
  iii(3,2,1:8)=(/0,0,0,0,6,6,6,6/); jjj(3,2,1:8)=(/5,6,7,8,1,2,3,4/)

  ! Loop over myid grids by vector sweeps
! do igrid=1,ncache,nvector

     ! Gather nvector grids
!    ngrid=MIN(nvector,ncache-igrid+1)
     do i=1,ngrid
        ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
        ! Gather neighboring grids
        igridn(i,0)=ind_grid(i)
     end do
     do idim=1,ndim
        do i=1,ngrid
           igridn(i,2*idim-1)=morton_nbor_grid(ind_grid(i),ilevel,2*idim-1)
           igridn(i,2*idim  )=morton_nbor_grid(ind_grid(i),ilevel,2*idim  )
           ind_left (i,idim)=morton_nbor_cell(ind_grid(i),ilevel,2*idim-1)
           ind_right(i,idim)=morton_nbor_cell(ind_grid(i),ilevel,2*idim  )
        end do
     end do

     ! Loop over cells
     do ind=1,twotondim

        ! Compute central cell index
        iskip=ncoarse+(ind-1)*ngridmax
        do i=1,ngrid
           ind_cell(i)=iskip+ind_grid(i)
        end do

        ! Gather all neighboring velocities
        do idim=1,ndim
           id1=jjj(idim,1,ind); ig1=iii(idim,1,ind)
           ih1=ncoarse+(id1-1)*ngridmax
           do i=1,ngrid
              if(igridn(i,ig1)>0)then
                 velg(i,idim,1:ndim) = uold(igridn(i,ig1)+ih1,2:ndim+1)/max(uold(igridn(i,ig1)+ih1,1),smallr)
                 dx_g(i,idim) = dx_loc
              else
                 velg(i,idim,1:ndim) = uold(ind_left(i,idim),2:ndim+1)/max(uold(ind_left(i,idim),1),smallr)
                 dx_g(i,idim) = dx_loc*1.5_dp
              end if
           enddo
           id2=jjj(idim,2,ind); ig2=iii(idim,2,ind)
           ih2=ncoarse+(id2-1)*ngridmax
           do i=1,ngrid
              if(igridn(i,ig2)>0)then
                 veld(i,idim,1:ndim)= uold(igridn(i,ig2)+ih2,2:ndim+1)/max(uold(igridn(i,ig2)+ih2,1),smallr)
                 dx_d(i,idim)=dx_loc
              else
                 veld(i,idim,1:ndim)= uold(ind_right(i,idim),2:ndim+1)/max(uold(ind_right(i,idim),1),smallr)
                 dx_d(i,idim)=dx_loc*1.5_dp
              end if
           enddo
        end do
        ! End loop over dimensions

        ! Compute divu = Trace G
        divu_loc(1:ngrid)=0.0d0
        do i=1,ngrid
           do idim=1,ndim
              divu_loc(i) = divu_loc(i) + (veld(i,idim,idim)-velg(i,idim,idim)) &
                   &                    / (dx_g(i,idim)     +dx_d(i,idim))
           enddo
        end do

        ! Update thermal internal energy
        if(pressure_fix)then
           do i=1,ngrid
              ! Compute old thermal energy
              d=max(uold(ind_cell(i),1),smallr)
              u=0.0; v=0.0; w=0.0
              if(ndim>0)u=uold(ind_cell(i),2)/d
              if(ndim>1)v=uold(ind_cell(i),3)/d
              if(ndim>2)w=uold(ind_cell(i),4)/d
              eold=uold(ind_cell(i),ndim+2)-0.5*d*(u**2+v**2+w**2)
#if NENER>0
              do irad=1,nener
                 eold=eold-uold(ind_cell(i),ndim+2+irad)
              end do
#endif
              ! Add -pdV term
              enew(ind_cell(i))=enew(ind_cell(i)) &
                   & -(gamma-1.0d0)*eold*divu_loc(i)*dtnew(ilevel)
           end do
        end if

#if NENER>0
        do irad=1,nener
           do i=1,ngrid
              ! Add -pdV term
              unew(ind_cell(i),ndim+2+irad)=unew(ind_cell(i),ndim+2+irad) &
                & -(gamma_rad(irad)-1.0d0)*uold(ind_cell(i),ndim+2+irad)*divu_loc(i)*dtnew(ilevel)
           end do
        end do
#endif

     enddo
     ! End loop over cells
! end do
  ! End loop over grids
end subroutine sub_add_pdv_source_terms





!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine godfine1(ilevel, jgrid,mgrid)
  use amr_commons
  use hydro_commons
  use poisson_commons
  use morton_hash
#ifdef HYDRO_CUDA
  use iso_c_binding, only: c_int, c_double
  use cuda_commons
  use hydro_cuda_interface
#endif
  implicit none
  integer, intent(in) ::ilevel,jgrid,mgrid
  integer,dimension(1:nvector)::ind_grid
  integer:: ngrid
  !-------------------------------------------------------------------
  ! This routine gathers first hydro variables from neighboring grids
  ! to set initial conditions in a 6x6x6 grid. It interpolate from
  ! coarser level missing grid variables. It then calls the
  ! Godunov solver that computes fluxes. These fluxes are zeroed at
  ! coarse-fine boundaries, since contribution from finer levels has
  ! already been taken into account. Conservative variables are updated
  ! and stored in array unew(:), both at the current level and at the
  ! coarser level if necessary.
  !-------------------------------------------------------------------
  integer ,dimension(1:nvector,1:threetondim     )::nbors_father_cells
  integer ,dimension(1:nvector,1:twotondim       )::nbors_father_grids
  integer ,dimension(1:nvector,0:twondim         )::ibuffer_father
  real(dp),dimension(1:nvector,0:twondim  ,1:nvar)::u1
  real(dp),dimension(1:nvector,1:twotondim,1:nvar)::u2
! real(dp),dimension(1:nvector,0:twondim  ,1:ndim)::g1=0.0d0
! real(dp),dimension(1:nvector,1:twotondim,1:ndim)::g2=0.0d0

  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar)::uloc
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:ndim)::gloc
  real(dp),dimension(1:nvector,if1:if2,jf1:jf2,kf1:kf2,1:nvar,1:ndim)::flux
  real(dp),dimension(1:nvector,if1:if2,jf1:jf2,kf1:kf2,1:2,1:ndim)::tmp
  logical ,dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2)::ok

  integer,dimension(1:nvector)::igrid_nbor,ind_cell,ind_cell2,ind_buffer, ind_buffer2,ind_exist,ind_nexist

  integer::i,j,ivar,idim,ind_son,ind_father,iskip,nbuffer,ibuffer
  integer::i0,j0,k0,i1,j1,k1,i2,j2,k2,i3,j3,k3,nx_loc,nb_noneigh,nb_noneigh2,nexist
  integer::igridn_tmp
  integer::i1min,i1max,j1min,j1max,k1min,k1max
  integer::i2min,i2max,j2min,j2max,k2min,k2max
  integer::i3min,i3max,j3min,j3max,k3min,k3max
  real(dp)::dx,scale,oneontwotondim

#ifdef HYDRO_CUDA
  integer :: stream_slot
  logical :: use_gpu
#endif


  gloc = 0.0d0
  ngrid = mgrid

  do i=1,ngrid
     ind_grid(i)=active(ilevel)%igrid(jgrid+i-1)
  end do

  oneontwotondim = 1.d0/dble(twotondim)

  ! Mesh spacing in that level
  nx_loc=icoarse_max-icoarse_min+1
  scale=boxlen/dble(nx_loc)
  dx=0.5D0**ilevel*scale

  ! Integer constants
  i1min=0; i1max=0; i2min=0; i2max=0; i3min=1; i3max=1
  j1min=0; j1max=0; j2min=0; j2max=0; j3min=1; j3max=1
  k1min=0; k1max=0; k2min=0; k2max=0; k3min=1; k3max=1
  if(ndim>0)then
     i1max=2; i2max=1; i3max=2
  end if
  if(ndim>1)then
     j1max=2; j2max=1; j3max=2
  end if
  if(ndim>2)then
     k1max=2; k2max=1; k3max=2
  end if

  !------------------------------------------
  ! Gather 3^ndim neighboring father cells
  !------------------------------------------
  do i=1,ngrid
     ind_cell(i)=father(ind_grid(i))
  end do
  call get3cubefather(ind_cell,nbors_father_cells,nbors_father_grids,ngrid,ilevel)

  !---------------------------
  ! Gather 6x6x6 cells stencil
  !---------------------------
  ! Loop over 3x3x3 neighboring father cells
  do k1=k1min,k1max
  do j1=j1min,j1max
  do i1=i1min,i1max

     ! Check if neighboring grid exists
     nbuffer=0
     nexist=0
     ind_father=1+i1+3*j1+9*k1
     do i=1,ngrid
        igrid_nbor(i)=son(nbors_father_cells(i,ind_father))
        if(igrid_nbor(i)>0) then
           nexist=nexist+1
           ind_exist(nexist)=i
        else
          nbuffer=nbuffer+1
          ind_nexist(nbuffer)=i
          ind_buffer(nbuffer)=nbors_father_cells(i,ind_father)
        end if
     end do

     ! If not, interpolate hydro variables from parent cells
     if(nbuffer>0)then
        call getnborfather(ind_buffer,ibuffer_father,nbuffer,ilevel)
        do j=0,twondim
           do ivar=1,nvar
              do i=1,nbuffer
                 u1(i,j,ivar)=uold(ibuffer_father(i,j),ivar)
              end do
           end do
        end do
        call interpol_hydro(u1,u2,nbuffer)
     endif

     ! Loop over 2x2x2 cells
     do k2=k2min,k2max
     do j2=j2min,j2max
     do i2=i2min,i2max

        ind_son=1+i2+2*j2+4*k2
        iskip=ncoarse+(ind_son-1)*ngridmax
        do i=1,nexist
           ind_cell(i)=iskip+igrid_nbor(ind_exist(i))
        end do

        i3=1; j3=1; k3=1
        if(ndim>0)i3=1+2*(i1-1)+i2
        if(ndim>1)j3=1+2*(j1-1)+j2
        if(ndim>2)k3=1+2*(k1-1)+k2

        ! Gather hydro variables
        do ivar=1,nvar
           do i=1,nexist
              uloc(ind_exist(i),i3,j3,k3,ivar)=uold(ind_cell(i),ivar)
           end do
           do i=1,nbuffer
              uloc(ind_nexist(i),i3,j3,k3,ivar)=u2(i,ind_son,ivar)
           end do
        end do

        ! Gather gravitational acceleration
        if(poisson)then
           do idim=1,ndim
              do i=1,nexist
                 gloc(ind_exist(i),i3,j3,k3,idim)=f(ind_cell(i),idim)
              end do
              ! Use straight injection for buffer cells
              do i=1,nbuffer
                 gloc(ind_nexist(i),i3,j3,k3,idim)=f(ibuffer_father(i,0),idim)
              end do
           end do
        end if

        ! Gather refinement flag
        do i=1,nexist
           ok(ind_exist(i),i3,j3,k3)=son(ind_cell(i))>0
        end do
        do i=1,nbuffer
           ok(ind_nexist(i),i3,j3,k3)=.false.
        end do

     end do
     end do
     end do
     ! End loop over cells

  end do
  end do
  end do
  ! End loop over neighboring grids

  !-----------------------------------------------
  ! Compute flux using second-order Godunov method
  !-----------------------------------------------
#ifdef HYDRO_CUDA
  use_gpu = .false.
  if (cuda_available .and. hydro_cuda_initialized) then
     stream_slot = int(cuda_acquire_stream_c())
     if (stream_slot >= 0) then
        use_gpu = .true.
        call hydro_cuda_unsplit_async_f( &
             uloc, gloc, flux, tmp, &
             real(dx, c_double), real(dx, c_double), real(dx, c_double), &
             real(dtnew(ilevel), c_double), &
             int(ngrid, c_int), int(nvector, c_int), int(stream_slot, c_int))
        call hydro_cuda_unsplit_sync_f( &
             flux, tmp, int(ngrid, c_int), int(stream_slot, c_int))
        call cuda_release_stream_c(int(stream_slot, c_int))
     end if
  end if
  if (.not. use_gpu) then
     call unsplit(uloc,gloc,flux,tmp,dx,dx,dx,dtnew(ilevel),ngrid)
  end if
#else
  call unsplit(uloc,gloc,flux,tmp,dx,dx,dx,dtnew(ilevel),ngrid)
#endif

  !------------------------------------------------
  ! Reset flux along direction at refined interface
  !------------------------------------------------
  do idim=1,ndim
     i0=0; j0=0; k0=0
     if(idim==1)i0=1
     if(idim==2)j0=1
     if(idim==3)k0=1
     do k3=k3min,k3max+k0
     do j3=j3min,j3max+j0
     do i3=i3min,i3max+i0
        do ivar=1,nvar
           do i=1,ngrid
              if(ok(i,i3-i0,j3-j0,k3-k0) .or. ok(i,i3,j3,k3))then
                 flux(i,i3,j3,k3,ivar,idim)=0.0d0
              end if
           end do
        end do
        if(pressure_fix)then
        do ivar=1,2
           do i=1,ngrid
              if(ok(i,i3-i0,j3-j0,k3-k0) .or. ok(i,i3,j3,k3))then
                 tmp (i,i3,j3,k3,ivar,idim)=0.0d0
              end if
           end do
        end do
        end if
     end do
     end do
     end do
  end do
  !--------------------------------------
  ! Conservative update at level ilevel
  !--------------------------------------
#ifndef OMP_TEST
!$omp flush
!$omp critical (godunov)
#endif
  do idim=1,ndim
     i0=0; j0=0; k0=0
     if(idim==1)i0=1
     if(idim==2)j0=1
     if(idim==3)k0=1
     do k2=k2min,k2max
     do j2=j2min,j2max
     do i2=i2min,i2max
        ind_son=1+i2+2*j2+4*k2
        iskip=ncoarse+(ind_son-1)*ngridmax
        do i=1,ngrid
           ind_cell(i)=iskip+ind_grid(i)
        end do
        i3=1+i2
        j3=1+j2
        k3=1+k2
        ! Update conservative variables new state vector
        do ivar=1,nvar
           do i=1,ngrid
              unew(ind_cell(i),ivar)=unew(ind_cell(i),ivar)+ &
                   & (flux(i,i3   ,j3   ,k3   ,ivar,idim) &
                   & -flux(i,i3+i0,j3+j0,k3+k0,ivar,idim))
           end do
        end do
        if(pressure_fix)then
        ! Update velocity divergence
        do i=1,ngrid
           divu(ind_cell(i))=divu(ind_cell(i))+ &
                & (tmp(i,i3   ,j3   ,k3   ,1,idim) &
                & -tmp(i,i3+i0,j3+j0,k3+k0,1,idim))
        end do
        ! Update internal energy
        do i=1,ngrid
           enew(ind_cell(i))=enew(ind_cell(i))+ &
                & (tmp(i,i3   ,j3   ,k3   ,2,idim) &
                & -tmp(i,i3+i0,j3+j0,k3+k0,2,idim))
        end do
        end if
     end do
     end do
     end do
  end do


  !--------------------------------------
  ! Conservative update at level ilevel-1
  !--------------------------------------
  ! Loop over dimensions
  do idim=1,ndim
     i0=0; j0=0; k0=0
     if(idim==1)i0=1
     if(idim==2)j0=1
     if(idim==3)k0=1

     !----------------------
     ! Left flux at boundary
     !----------------------
     ! Check if grids sits near left boundary
     ! and gather neighbor father cells index
     nb_noneigh=0
     do i=1,ngrid
        igridn_tmp = morton_nbor_grid(ind_grid(i),ilevel,2*idim-1)
        if (igridn_tmp==0) then
           nb_noneigh = nb_noneigh + 1
           ind_buffer(nb_noneigh) = morton_nbor_cell(ind_grid(i),ilevel,2*idim-1)
           ind_cell(nb_noneigh) = i
        end if
     end do

     !-----------------------
     ! Right flux at boundary
     !-----------------------
     ! Check if grids sits near right boundary
     ! and gather neighbor father cells index
     nb_noneigh2=0
     do i=1,ngrid
        igridn_tmp = morton_nbor_grid(ind_grid(i),ilevel,2*idim)
        if (igridn_tmp==0) then
           nb_noneigh2 = nb_noneigh2 + 1
           ind_buffer2(nb_noneigh2) = morton_nbor_cell(ind_grid(i),ilevel,2*idim)
           ind_cell2(nb_noneigh2) = i
        end if
     end do
     ! Conservative update of new state variables
     do ivar=1,nvar
        ! Loop over boundary cells
        do k3=k3min,k3max-k0
        do j3=j3min,j3max-j0
        do i3=i3min,i3max-i0
           do i=1,nb_noneigh
              unew(ind_buffer(i),ivar)=unew(ind_buffer(i),ivar) &
                   & -flux(ind_cell(i),i3,j3,k3,ivar,idim)*oneontwotondim
           end do
        end do
        end do
        end do
        ! Loop over boundary cells
        do k3=k3min+k0,k3max
        do j3=j3min+j0,j3max
        do i3=i3min+i0,i3max
           do i=1,nb_noneigh2
              unew(ind_buffer2(i),ivar)=unew(ind_buffer2(i),ivar) &
                   & +flux(ind_cell2(i),i3+i0,j3+j0,k3+k0,ivar,idim)*oneontwotondim
           end do
        end do
        end do
        end do
     end do
     if(pressure_fix)then
        ! Update velocity divergence
        do k3=k3min,k3max-k0
        do j3=j3min,j3max-j0
        do i3=i3min,i3max-i0
           do i=1,nb_noneigh
              divu(ind_buffer(i))=divu(ind_buffer(i)) &
                   & -tmp(ind_cell(i),i3,j3,k3,1,idim)*oneontwotondim
           end do
        end do
        end do
        end do
        ! Update velocity divergence
        do k3=k3min+k0,k3max
        do j3=j3min+j0,j3max
        do i3=i3min+i0,i3max
           do i=1,nb_noneigh2
              divu(ind_buffer2(i))=divu(ind_buffer2(i)) &
                   & +tmp(ind_cell2(i),i3+i0,j3+j0,k3+k0,1,idim)*oneontwotondim
           end do
        end do
        end do
        end do
        ! Update internal energy
        do k3=k3min,k3max-k0
        do j3=j3min,j3max-j0
        do i3=i3min,i3max-i0
           do i=1,nb_noneigh
              enew(ind_buffer(i))=enew(ind_buffer(i)) &
                   & -tmp(ind_cell(i),i3,j3,k3,2,idim)*oneontwotondim
           end do
        end do
        end do
        end do
        ! Update internal energy
        do k3=k3min+k0,k3max
        do j3=j3min+j0,j3max
        do i3=i3min+i0,i3max
           do i=1,nb_noneigh2
              enew(ind_buffer2(i))=enew(ind_buffer2(i)) &
                   & +tmp(ind_cell2(i),i3+i0,j3+j0,k3+k0,2,idim)*oneontwotondim
           end do
        end do
        end do
        end do
     end if

  end do
#ifndef OMP_TEST
!$omp flush
!$omp end critical (godunov)
#endif
  ! End loop over dimensions

end subroutine godfine1
#ifdef HYDRO_CUDA
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine godunov_fine_hybrid(ilevel, ncache)
  use omp_lib
  use amr_commons
  use hydro_commons
  use hydro_parameters
  use hydro_hybrid_commons
  use iso_c_binding, only: c_int
  use cuda_commons
  use hydro_cuda_interface
  implicit none
  integer, intent(in) :: ilevel, ncache
  !-------------------------------------------------------------------
  ! Dynamic hybrid CPU/GPU dispatcher for Godunov solver.
  ! OMP threads acquire GPU streams dynamically; fallback to CPU.
  ! Level L scatter: direct to unew (conflict-free).
  ! Level L-1 scatter: per-thread buffer, serial merge after barrier.
  !-------------------------------------------------------------------
  integer :: tid, nthreads, stream_slot, igrid, ngrid
  integer :: it, i, ic, ivar
  integer(kind=8) :: t_start, t_now, clock_rate

  call system_clock(count_rate=clock_rate)

  ! First-call initialization (serial)
  if (hybrid_max_threads == 0) then
     nthreads = 1
     !$ nthreads = omp_get_max_threads()
     hybrid_max_threads = nthreads
     allocate(scatter_bufs(0:nthreads-1))
     do it = 0, nthreads-1
        call scatter_buf_init(scatter_bufs(it), SBUF_INIT_CAP, nvar)
     end do
     allocate(gpu_states(0:max(cuda_n_streams,1)-1))
  end if

  ! Reset scatter buffer counts (serial)
  do it = 0, hybrid_max_threads-1
     scatter_bufs(it)%count = 0
  end do

  !$omp parallel private(tid, stream_slot, igrid, ngrid)
  tid = 0
  !$ tid = omp_get_thread_num()
  stream_slot = int(cuda_acquire_stream_c())

  ! Ensure GPU state allocated (lazy, per-stream-slot)
  if (stream_slot >= 0) then
     if (.not. allocated(gpu_states(stream_slot)%super_uloc)) then
        allocate(gpu_states(stream_slot)%super_uloc(1:HYBRID_SUPER_SIZE,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar))
        allocate(gpu_states(stream_slot)%super_gloc(1:HYBRID_SUPER_SIZE,iu1:iu2,ju1:ju2,ku1:ku2,1:ndim))
        allocate(gpu_states(stream_slot)%super_flux(1:HYBRID_SUPER_SIZE,if1:if2,jf1:jf2,kf1:kf2,1:nvar,1:ndim))
        allocate(gpu_states(stream_slot)%super_tmp (1:HYBRID_SUPER_SIZE,if1:if2,jf1:jf2,kf1:kf2,1:2,1:ndim))
        allocate(gpu_states(stream_slot)%super_ok  (1:HYBRID_SUPER_SIZE,iu1:iu2,ju1:ju2,ku1:ku2))
        allocate(gpu_states(stream_slot)%super_ind_grid(1:HYBRID_SUPER_SIZE))
        gpu_states(stream_slot)%super_gloc = 0.0d0
     end if
     gpu_states(stream_slot)%off = 0
  end if

  !$omp do schedule(dynamic)
  do igrid = 1, ncache, nvector
     ngrid = MIN(nvector, ncache - igrid + 1)
     if (stream_slot >= 0) then
        call hybrid_gpu_gather_batch(gpu_states(stream_slot), ilevel, igrid, ngrid, &
             scatter_bufs(tid), stream_slot)
     else
        call hybrid_cpu_process_batch(ilevel, igrid, ngrid, scatter_bufs(tid))
     end if
  end do
  !$omp end do nowait

  ! GPU threads: flush remaining super buffer
  if (stream_slot >= 0) then
     if (gpu_states(stream_slot)%off > 0) then
        call hybrid_gpu_flush_scatter(gpu_states(stream_slot), stream_slot, &
             ilevel, scatter_bufs(tid))
     end if
     call cuda_release_stream_c(int(stream_slot, c_int))
  end if
  !$omp end parallel

  ! Serial merge of L-1 scatter buffers
  call system_clock(t_start)
  do it = 0, hybrid_max_threads-1
     do i = 1, scatter_bufs(it)%count
        ic = scatter_bufs(it)%icell(i)
        do ivar = 1, nvar
           unew(ic, ivar) = unew(ic, ivar) + scatter_bufs(it)%dunew(i, ivar)
        end do
        if (pressure_fix) then
           divu(ic) = divu(ic) + scatter_bufs(it)%ddivu(i)
           enew(ic) = enew(ic) + scatter_bufs(it)%denew(i)
        end if
     end do
     scatter_bufs(it)%count = 0
  end do
  call system_clock(t_now)
  acc_merge = acc_merge + dble(t_now - t_start) / dble(clock_rate)

  ! Print timing at last step
  if (nstep >= nstepmax .and. myid == 1) then
     write(*,'(A)') ' === Hybrid CPU/GPU timing ==='
     write(*,'(A,F8.3,A)') '   gather      : ', acc_gather, ' s'
     write(*,'(A,F8.3,A)') '   GPU compute : ', acc_gpu, ' s'
     write(*,'(A,F8.3,A)') '   L scatter   : ', acc_scatter_l, ' s'
     write(*,'(A,F8.3,A)') '   L-1 merge   : ', acc_merge, ' s'
     write(*,'(A,I6)') '   GPU flushes : ', n_gpu_flushes
     write(*,'(A,I6)') '   CPU batches : ', n_cpu_batches
  end if

end subroutine godunov_fine_hybrid
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine hybrid_cpu_process_batch(ilevel, igrid_start, ngrid, sbuf)
  use amr_commons
  use hydro_commons
  use hydro_parameters
  use poisson_commons
  use morton_hash
  use hydro_hybrid_commons
  implicit none
  integer, intent(in) :: ilevel, igrid_start, ngrid
  type(scatter_buf_t), intent(inout) :: sbuf
  !-------------------------------------------------------------------
  ! CPU path: gather + unsplit + scatter.
  ! Level L scatter: direct to unew (conflict-free, no critical).
  ! Level L-1 scatter: pre-summed entries appended to per-thread buffer.
  !-------------------------------------------------------------------
  integer ,dimension(1:nvector)::ind_grid
  integer ,dimension(1:nvector,1:threetondim     )::nbors_father_cells
  integer ,dimension(1:nvector,1:twotondim       )::nbors_father_grids
  integer ,dimension(1:nvector,0:twondim         )::ibuffer_father
  real(dp),dimension(1:nvector,0:twondim  ,1:nvar)::u1
  real(dp),dimension(1:nvector,1:twotondim,1:nvar)::u2
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar)::uloc
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:ndim)::gloc
  real(dp),dimension(1:nvector,if1:if2,jf1:jf2,kf1:kf2,1:nvar,1:ndim)::flux
  real(dp),dimension(1:nvector,if1:if2,jf1:jf2,kf1:kf2,1:2,1:ndim)::tmp
  logical ,dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2)::ok
  integer ,dimension(1:nvector)::igrid_nbor,ind_cell,ind_buffer,ind_exist,ind_nexist

  integer :: i,j,ivar,idim,ind_son,ind_father,iskip,nbuffer,nexist
  integer :: i0,j0,k0,i1,j1,k1,i2,j2,k2,i3,j3,k3,nx_loc
  integer :: igridn_tmp, icell_nbor, idx
  integer :: i1min,i1max,j1min,j1max,k1min,k1max
  integer :: i2min,i2max,j2min,j2max,k2min,k2max
  integer :: i3min,i3max,j3min,j3max,k3min,k3max
  real(dp) :: dx,scale,oneontwotondim
  real(dp) :: acc_unew(1:nvar), acc_ddivu, acc_denew

  gloc = 0.0d0

  do i=1,ngrid
     ind_grid(i) = active(ilevel)%igrid(igrid_start+i-1)
  end do

  oneontwotondim = 1.d0/dble(twotondim)
  nx_loc = icoarse_max-icoarse_min+1
  scale = boxlen/dble(nx_loc)
  dx = 0.5D0**ilevel*scale

  ! Integer constants
  i1min=0; i1max=0; i2min=0; i2max=0; i3min=1; i3max=1
  j1min=0; j1max=0; j2min=0; j2max=0; j3min=1; j3max=1
  k1min=0; k1max=0; k2min=0; k2max=0; k3min=1; k3max=1
  if(ndim>0)then
     i1max=2; i2max=1; i3max=2
  end if
  if(ndim>1)then
     j1max=2; j2max=1; j3max=2
  end if
  if(ndim>2)then
     k1max=2; k2max=1; k3max=2
  end if

  !------------------------------------------
  ! Gather 3^ndim neighboring father cells
  !------------------------------------------
  do i=1,ngrid
     ind_cell(i) = father(ind_grid(i))
  end do
  call get3cubefather(ind_cell,nbors_father_cells,nbors_father_grids,ngrid,ilevel)

  !---------------------------
  ! Gather 6x6x6 cells stencil
  !---------------------------
  do k1=k1min,k1max
  do j1=j1min,j1max
  do i1=i1min,i1max

     nbuffer=0
     nexist=0
     ind_father=1+i1+3*j1+9*k1
     do i=1,ngrid
        igrid_nbor(i) = son(nbors_father_cells(i,ind_father))
        if(igrid_nbor(i)>0) then
           nexist=nexist+1
           ind_exist(nexist)=i
        else
           nbuffer=nbuffer+1
           ind_nexist(nbuffer)=i
           ind_buffer(nbuffer)=nbors_father_cells(i,ind_father)
        end if
     end do

     if(nbuffer>0)then
        call getnborfather(ind_buffer,ibuffer_father,nbuffer,ilevel)
        do j=0,twondim
           do ivar=1,nvar
              do i=1,nbuffer
                 u1(i,j,ivar)=uold(ibuffer_father(i,j),ivar)
              end do
           end do
        end do
        call interpol_hydro(u1,u2,nbuffer)
     endif

     do k2=k2min,k2max
     do j2=j2min,j2max
     do i2=i2min,i2max

        ind_son=1+i2+2*j2+4*k2
        iskip=ncoarse+(ind_son-1)*ngridmax
        do i=1,nexist
           ind_cell(i)=iskip+igrid_nbor(ind_exist(i))
        end do

        i3=1; j3=1; k3=1
        if(ndim>0)i3=1+2*(i1-1)+i2
        if(ndim>1)j3=1+2*(j1-1)+j2
        if(ndim>2)k3=1+2*(k1-1)+k2

        do ivar=1,nvar
           do i=1,nexist
              uloc(ind_exist(i),i3,j3,k3,ivar)=uold(ind_cell(i),ivar)
           end do
           do i=1,nbuffer
              uloc(ind_nexist(i),i3,j3,k3,ivar)=u2(i,ind_son,ivar)
           end do
        end do

        if(poisson)then
           do idim=1,ndim
              do i=1,nexist
                 gloc(ind_exist(i),i3,j3,k3,idim)=f(ind_cell(i),idim)
              end do
              do i=1,nbuffer
                 gloc(ind_nexist(i),i3,j3,k3,idim)=f(ibuffer_father(i,0),idim)
              end do
           end do
        end if

        do i=1,nexist
           ok(ind_exist(i),i3,j3,k3)=son(ind_cell(i))>0
        end do
        do i=1,nbuffer
           ok(ind_nexist(i),i3,j3,k3)=.false.
        end do

     end do
     end do
     end do

  end do
  end do
  end do

  !-----------------------------------------------
  ! Compute flux using second-order Godunov method
  !-----------------------------------------------
  call unsplit(uloc,gloc,flux,tmp,dx,dx,dx,dtnew(ilevel),ngrid)

  !------------------------------------------------
  ! Reset flux along direction at refined interface
  !------------------------------------------------
  do idim=1,ndim
     i0=0; j0=0; k0=0
     if(idim==1)i0=1
     if(idim==2)j0=1
     if(idim==3)k0=1
     do k3=k3min,k3max+k0
     do j3=j3min,j3max+j0
     do i3=i3min,i3max+i0
        do ivar=1,nvar
           do i=1,ngrid
              if(ok(i,i3-i0,j3-j0,k3-k0) .or. ok(i,i3,j3,k3))then
                 flux(i,i3,j3,k3,ivar,idim)=0.0d0
              end if
           end do
        end do
        if(pressure_fix)then
        do ivar=1,2
           do i=1,ngrid
              if(ok(i,i3-i0,j3-j0,k3-k0) .or. ok(i,i3,j3,k3))then
                 tmp (i,i3,j3,k3,ivar,idim)=0.0d0
              end if
           end do
        end do
        end if
     end do
     end do
     end do
  end do

  !--------------------------------------
  ! Level L scatter (direct, conflict-free)
  !--------------------------------------
  do idim=1,ndim
     i0=0; j0=0; k0=0
     if(idim==1)i0=1
     if(idim==2)j0=1
     if(idim==3)k0=1
     do k2=k2min,k2max
     do j2=j2min,j2max
     do i2=i2min,i2max
        ind_son=1+i2+2*j2+4*k2
        iskip=ncoarse+(ind_son-1)*ngridmax
        do i=1,ngrid
           ind_cell(i)=iskip+ind_grid(i)
        end do
        i3=1+i2
        j3=1+j2
        k3=1+k2
        do ivar=1,nvar
           do i=1,ngrid
              unew(ind_cell(i),ivar)=unew(ind_cell(i),ivar)+ &
                   & (flux(i,i3   ,j3   ,k3   ,ivar,idim) &
                   & -flux(i,i3+i0,j3+j0,k3+k0,ivar,idim))
           end do
        end do
        if(pressure_fix)then
        do i=1,ngrid
           divu(ind_cell(i))=divu(ind_cell(i))+ &
                & (tmp(i,i3   ,j3   ,k3   ,1,idim) &
                & -tmp(i,i3+i0,j3+j0,k3+k0,1,idim))
        end do
        do i=1,ngrid
           enew(ind_cell(i))=enew(ind_cell(i))+ &
                & (tmp(i,i3   ,j3   ,k3   ,2,idim) &
                & -tmp(i,i3+i0,j3+j0,k3+k0,2,idim))
        end do
        end if
     end do
     end do
     end do
  end do

  !--------------------------------------
  ! Level L-1 scatter (to per-thread buffer)
  !--------------------------------------
  do i=1,ngrid
     do idim=1,ndim
        i0=0; j0=0; k0=0
        if(idim==1)i0=1
        if(idim==2)j0=1
        if(idim==3)k0=1

        ! Left boundary
        igridn_tmp = morton_nbor_grid(ind_grid(i),ilevel,2*idim-1)
        if (igridn_tmp==0) then
           icell_nbor = morton_nbor_cell(ind_grid(i),ilevel,2*idim-1)
           acc_unew(1:nvar) = 0.0d0
           acc_ddivu = 0.0d0; acc_denew = 0.0d0
           do k3=k3min,k3max-k0
           do j3=j3min,j3max-j0
           do i3=i3min,i3max-i0
              do ivar=1,nvar
                 acc_unew(ivar) = acc_unew(ivar) &
                      & - flux(i,i3,j3,k3,ivar,idim)*oneontwotondim
              end do
              if(pressure_fix)then
                 acc_ddivu = acc_ddivu - tmp(i,i3,j3,k3,1,idim)*oneontwotondim
                 acc_denew = acc_denew - tmp(i,i3,j3,k3,2,idim)*oneontwotondim
              end if
           end do
           end do
           end do
           sbuf%count = sbuf%count + 1
           idx = sbuf%count
           if (idx > sbuf%capacity) call scatter_buf_grow(sbuf)
           sbuf%icell(idx) = icell_nbor
           sbuf%dunew(idx, 1:nvar) = acc_unew(1:nvar)
           if(pressure_fix) then
              sbuf%ddivu(idx) = acc_ddivu
              sbuf%denew(idx) = acc_denew
           end if
        end if

        ! Right boundary
        igridn_tmp = morton_nbor_grid(ind_grid(i),ilevel,2*idim)
        if (igridn_tmp==0) then
           icell_nbor = morton_nbor_cell(ind_grid(i),ilevel,2*idim)
           acc_unew(1:nvar) = 0.0d0
           acc_ddivu = 0.0d0; acc_denew = 0.0d0
           do k3=k3min+k0,k3max
           do j3=j3min+j0,j3max
           do i3=i3min+i0,i3max
              do ivar=1,nvar
                 acc_unew(ivar) = acc_unew(ivar) &
                      & + flux(i,i3+i0,j3+j0,k3+k0,ivar,idim)*oneontwotondim
              end do
              if(pressure_fix)then
                 acc_ddivu = acc_ddivu + tmp(i,i3+i0,j3+j0,k3+k0,1,idim)*oneontwotondim
                 acc_denew = acc_denew + tmp(i,i3+i0,j3+j0,k3+k0,2,idim)*oneontwotondim
              end if
           end do
           end do
           end do
           sbuf%count = sbuf%count + 1
           idx = sbuf%count
           if (idx > sbuf%capacity) call scatter_buf_grow(sbuf)
           sbuf%icell(idx) = icell_nbor
           sbuf%dunew(idx, 1:nvar) = acc_unew(1:nvar)
           if(pressure_fix) then
              sbuf%ddivu(idx) = acc_ddivu
              sbuf%denew(idx) = acc_denew
           end if
        end if

     end do
  end do

  !$omp atomic
  n_cpu_batches = n_cpu_batches + 1

end subroutine hybrid_cpu_process_batch
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine hybrid_gpu_gather_batch(gstate, ilevel, igrid_start, ngrid, sbuf, stream_slot)
  use amr_commons
  use hydro_commons
  use hydro_parameters
  use poisson_commons
  use morton_hash
  use hydro_hybrid_commons
  use iso_c_binding, only: c_int
  implicit none
  type(gpu_state_t), intent(inout) :: gstate
  integer, intent(in) :: ilevel, igrid_start, ngrid, stream_slot
  type(scatter_buf_t), intent(inout) :: sbuf
  !-------------------------------------------------------------------
  ! GPU path: gather stencil into super buffer, flush when full.
  !-------------------------------------------------------------------
  integer ,dimension(1:nvector)::ind_grid_loc,igrid_nbor,ind_cell
  integer ,dimension(1:nvector)::ind_exist,ind_nexist,ind_buffer
  integer ,dimension(1:nvector,1:threetondim)::nbors_father_cells
  integer ,dimension(1:nvector,1:twotondim)::nbors_father_grids
  integer ,dimension(1:nvector,0:twondim)::ibuffer_father
  real(dp),dimension(1:nvector,0:twondim,1:nvar)::u1
  real(dp),dimension(1:nvector,1:twotondim,1:nvar)::u2

  integer :: off, i, j, ivar, idim, ind_son, iskip, ind_father, nbuffer, nexist
  integer :: i1,j1,k1,i2,j2,k2,i3,j3,k3
  integer :: i1min,i1max,j1min,j1max,k1min,k1max
  integer :: i2min,i2max,j2min,j2max,k2min,k2max
  integer :: i3min,i3max,j3min,j3max,k3min,k3max
  integer(kind=8) :: t_start, t_now, clock_rate

  call system_clock(count_rate=clock_rate)
  call system_clock(t_start)

  off = gstate%off

  ! Integer constants
  i1min=0; i1max=0; i2min=0; i2max=0; i3min=1; i3max=1
  j1min=0; j1max=0; j2min=0; j2max=0; j3min=1; j3max=1
  k1min=0; k1max=0; k2min=0; k2max=0; k3min=1; k3max=1
  if(ndim>0)then
     i1max=2; i2max=1; i3max=2
  end if
  if(ndim>1)then
     j1max=2; j2max=1; j3max=2
  end if
  if(ndim>2)then
     k1max=2; k2max=1; k3max=2
  end if

  do i=1,ngrid
     ind_grid_loc(i) = active(ilevel)%igrid(igrid_start+i-1)
  end do

  ! Gather father cells
  do i=1,ngrid
     ind_cell(i) = father(ind_grid_loc(i))
  end do
  call get3cubefather(ind_cell,nbors_father_cells,nbors_father_grids,ngrid,ilevel)

  ! Gather 6x6x6 stencil into super buffer at offset 'off'
  do k1=k1min,k1max
  do j1=j1min,j1max
  do i1=i1min,i1max
     nbuffer=0
     nexist=0
     ind_father=1+i1+3*j1+9*k1
     do i=1,ngrid
        igrid_nbor(i)=son(nbors_father_cells(i,ind_father))
        if(igrid_nbor(i)>0)then
           nexist=nexist+1
           ind_exist(nexist)=i
        else
           nbuffer=nbuffer+1
           ind_nexist(nbuffer)=i
           ind_buffer(nbuffer)=nbors_father_cells(i,ind_father)
        end if
     end do

     if(nbuffer>0)then
        call getnborfather(ind_buffer,ibuffer_father,nbuffer,ilevel)
        do j=0,twondim
           do ivar=1,nvar
              do i=1,nbuffer
                 u1(i,j,ivar)=uold(ibuffer_father(i,j),ivar)
              end do
           end do
        end do
        call interpol_hydro(u1,u2,nbuffer)
     endif

     do k2=k2min,k2max
     do j2=j2min,j2max
     do i2=i2min,i2max
        ind_son=1+i2+2*j2+4*k2
        iskip=ncoarse+(ind_son-1)*ngridmax
        do i=1,nexist
           ind_cell(i)=iskip+igrid_nbor(ind_exist(i))
        end do

        i3=1; j3=1; k3=1
        if(ndim>0)i3=1+2*(i1-1)+i2
        if(ndim>1)j3=1+2*(j1-1)+j2
        if(ndim>2)k3=1+2*(k1-1)+k2

        do ivar=1,nvar
           do i=1,nexist
              gstate%super_uloc(off+ind_exist(i),i3,j3,k3,ivar)=uold(ind_cell(i),ivar)
           end do
           do i=1,nbuffer
              gstate%super_uloc(off+ind_nexist(i),i3,j3,k3,ivar)=u2(i,ind_son,ivar)
           end do
        end do

        if(poisson)then
           do idim=1,ndim
              do i=1,nexist
                 gstate%super_gloc(off+ind_exist(i),i3,j3,k3,idim)=f(ind_cell(i),idim)
              end do
              do i=1,nbuffer
                 gstate%super_gloc(off+ind_nexist(i),i3,j3,k3,idim)=f(ibuffer_father(i,0),idim)
              end do
           end do
        end if

        do i=1,nexist
           gstate%super_ok(off+ind_exist(i),i3,j3,k3)=son(ind_cell(i))>0
        end do
        do i=1,nbuffer
           gstate%super_ok(off+ind_nexist(i),i3,j3,k3)=.false.
        end do

     end do
     end do
     end do
  end do
  end do
  end do

  ! Record grid indices
  gstate%super_ind_grid(off+1:off+ngrid) = ind_grid_loc(1:ngrid)
  gstate%off = off + ngrid

  call system_clock(t_now)
  !$omp atomic
  acc_gather = acc_gather + dble(t_now - t_start) / dble(clock_rate)

  ! Flush if buffer full
  if (gstate%off + nvector > HYBRID_SUPER_SIZE) then
     call hybrid_gpu_flush_scatter(gstate, stream_slot, ilevel, sbuf)
  end if

end subroutine hybrid_gpu_gather_batch
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine hybrid_gpu_flush_scatter(gstate, stream_slot, ilevel, sbuf)
  use amr_commons
  use hydro_commons
  use hydro_parameters
  use morton_hash
  use hydro_hybrid_commons
  use iso_c_binding, only: c_int, c_double
  use cuda_commons
  use hydro_cuda_interface
  implicit none
  type(gpu_state_t), intent(inout) :: gstate
  integer, intent(in) :: stream_slot, ilevel
  type(scatter_buf_t), intent(inout) :: sbuf
  !-------------------------------------------------------------------
  ! GPU compute + scatter for accumulated super buffer.
  ! Level L scatter: direct to unew (conflict-free).
  ! Level L-1 scatter: pre-summed entries to per-thread buffer.
  !-------------------------------------------------------------------
  integer :: super_total, g, ig, ic, igridn_tmp, icell_nbor, idx
  integer :: i0, j0, k0, i2, j2, k2, i3, j3, k3, ind_son, iskip, ivar, idim
  integer :: i2min,i2max,j2min,j2max,k2min,k2max
  integer :: i3min,i3max,j3min,j3max,k3min,k3max
  integer :: nx_loc
  real(dp) :: dx, scale, oneontwotondim, dt_val
  real(dp) :: acc_unew(1:nvar), acc_ddivu, acc_denew
  integer(kind=8) :: t_start, t_now, clock_rate

  call system_clock(count_rate=clock_rate)

  super_total = gstate%off
  if (super_total == 0) return

  nx_loc = icoarse_max-icoarse_min+1
  scale = boxlen/dble(nx_loc)
  dx = 0.5D0**ilevel*scale
  oneontwotondim = 1.d0/dble(twotondim)

  i2min=0; i2max=0; i3min=1; i3max=1
  j2min=0; j2max=0; j3min=1; j3max=1
  k2min=0; k2max=0; k3min=1; k3max=1
  if(ndim>0)then
     i2max=1; i3max=2
  end if
  if(ndim>1)then
     j2max=1; j3max=2
  end if
  if(ndim>2)then
     k2max=1; k3max=2
  end if

  ! GPU compute
  call system_clock(t_start)
  call hydro_cuda_unsplit_async_f( &
       gstate%super_uloc, gstate%super_gloc, gstate%super_flux, gstate%super_tmp, &
       real(dx, c_double), real(dx, c_double), real(dx, c_double), &
       real(dtnew(ilevel), c_double), &
       int(super_total, c_int), int(HYBRID_SUPER_SIZE, c_int), int(stream_slot, c_int))
  call hydro_cuda_unsplit_sync_f( &
       gstate%super_flux, gstate%super_tmp, int(super_total, c_int), int(stream_slot, c_int))
  call system_clock(t_now)
  dt_val = dble(t_now - t_start) / dble(clock_rate)
  !$omp atomic
  acc_gpu = acc_gpu + dt_val

  ! Scatter phase
  call system_clock(t_start)

  ! Flux reset at refined interfaces
  do g=1,super_total
     do idim=1,ndim
        i0=0; j0=0; k0=0
        if(idim==1)i0=1
        if(idim==2)j0=1
        if(idim==3)k0=1
        do k3=k3min,k3max+k0
        do j3=j3min,j3max+j0
        do i3=i3min,i3max+i0
           if(gstate%super_ok(g,i3-i0,j3-j0,k3-k0) .or. gstate%super_ok(g,i3,j3,k3))then
              do ivar=1,nvar
                 gstate%super_flux(g,i3,j3,k3,ivar,idim)=0.0d0
              end do
              if(pressure_fix)then
                 gstate%super_tmp(g,i3,j3,k3,1,idim)=0.0d0
                 gstate%super_tmp(g,i3,j3,k3,2,idim)=0.0d0
              end if
           end if
        end do
        end do
        end do
     end do
  end do

  ! Level L scatter (direct, conflict-free)
  do g=1,super_total
     ig = gstate%super_ind_grid(g)
     do idim=1,ndim
        i0=0; j0=0; k0=0
        if(idim==1)i0=1
        if(idim==2)j0=1
        if(idim==3)k0=1
        do k2=k2min,k2max
        do j2=j2min,j2max
        do i2=i2min,i2max
           ind_son=1+i2+2*j2+4*k2
           iskip=ncoarse+(ind_son-1)*ngridmax
           ic = iskip + ig
           i3=1+i2
           j3=1+j2
           k3=1+k2
           do ivar=1,nvar
              unew(ic,ivar)=unew(ic,ivar)+ &
                   & (gstate%super_flux(g,i3   ,j3   ,k3   ,ivar,idim) &
                   & -gstate%super_flux(g,i3+i0,j3+j0,k3+k0,ivar,idim))
           end do
           if(pressure_fix)then
              divu(ic)=divu(ic)+ &
                   & (gstate%super_tmp(g,i3   ,j3   ,k3   ,1,idim) &
                   & -gstate%super_tmp(g,i3+i0,j3+j0,k3+k0,1,idim))
              enew(ic)=enew(ic)+ &
                   & (gstate%super_tmp(g,i3   ,j3   ,k3   ,2,idim) &
                   & -gstate%super_tmp(g,i3+i0,j3+j0,k3+k0,2,idim))
           end if
        end do
        end do
        end do
     end do
  end do

  ! Level L-1 scatter (to per-thread buffer)
  do g=1,super_total
     ig = gstate%super_ind_grid(g)
     do idim=1,ndim
        i0=0; j0=0; k0=0
        if(idim==1)i0=1
        if(idim==2)j0=1
        if(idim==3)k0=1

        ! Left boundary
        igridn_tmp = morton_nbor_grid(ig,ilevel,2*idim-1)
        if (igridn_tmp==0) then
           icell_nbor = morton_nbor_cell(ig,ilevel,2*idim-1)
           acc_unew(1:nvar) = 0.0d0
           acc_ddivu = 0.0d0; acc_denew = 0.0d0
           do k3=k3min,k3max-k0
           do j3=j3min,j3max-j0
           do i3=i3min,i3max-i0
              do ivar=1,nvar
                 acc_unew(ivar) = acc_unew(ivar) &
                      & - gstate%super_flux(g,i3,j3,k3,ivar,idim)*oneontwotondim
              end do
              if(pressure_fix)then
                 acc_ddivu = acc_ddivu - gstate%super_tmp(g,i3,j3,k3,1,idim)*oneontwotondim
                 acc_denew = acc_denew - gstate%super_tmp(g,i3,j3,k3,2,idim)*oneontwotondim
              end if
           end do
           end do
           end do
           sbuf%count = sbuf%count + 1
           idx = sbuf%count
           if (idx > sbuf%capacity) call scatter_buf_grow(sbuf)
           sbuf%icell(idx) = icell_nbor
           sbuf%dunew(idx, 1:nvar) = acc_unew(1:nvar)
           if(pressure_fix) then
              sbuf%ddivu(idx) = acc_ddivu
              sbuf%denew(idx) = acc_denew
           end if
        end if

        ! Right boundary
        igridn_tmp = morton_nbor_grid(ig,ilevel,2*idim)
        if (igridn_tmp==0) then
           icell_nbor = morton_nbor_cell(ig,ilevel,2*idim)
           acc_unew(1:nvar) = 0.0d0
           acc_ddivu = 0.0d0; acc_denew = 0.0d0
           do k3=k3min+k0,k3max
           do j3=j3min+j0,j3max
           do i3=i3min+i0,i3max
              do ivar=1,nvar
                 acc_unew(ivar) = acc_unew(ivar) &
                      & + gstate%super_flux(g,i3+i0,j3+j0,k3+k0,ivar,idim)*oneontwotondim
              end do
              if(pressure_fix)then
                 acc_ddivu = acc_ddivu + gstate%super_tmp(g,i3+i0,j3+j0,k3+k0,1,idim)*oneontwotondim
                 acc_denew = acc_denew + gstate%super_tmp(g,i3+i0,j3+j0,k3+k0,2,idim)*oneontwotondim
              end if
           end do
           end do
           end do
           sbuf%count = sbuf%count + 1
           idx = sbuf%count
           if (idx > sbuf%capacity) call scatter_buf_grow(sbuf)
           sbuf%icell(idx) = icell_nbor
           sbuf%dunew(idx, 1:nvar) = acc_unew(1:nvar)
           if(pressure_fix) then
              sbuf%ddivu(idx) = acc_ddivu
              sbuf%denew(idx) = acc_denew
           end if
        end if

     end do
  end do

  call system_clock(t_now)
  dt_val = dble(t_now - t_start) / dble(clock_rate)
  !$omp atomic
  acc_scatter_l = acc_scatter_l + dt_val
  !$omp atomic
  n_gpu_flushes = n_gpu_flushes + 1

  gstate%off = 0

end subroutine hybrid_gpu_flush_scatter
#endif
