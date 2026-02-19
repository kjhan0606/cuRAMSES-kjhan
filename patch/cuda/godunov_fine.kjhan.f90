!###########################################################
! Module for scatter buffer (always compiled, no #ifdef)
!###########################################################
module hydro_scatter_commons
  use amr_parameters, only: dp
  implicit none
  integer, parameter :: SBUF_INIT_CAP = 65536

  type scatter_buf_t
     integer :: count = 0
     integer :: capacity = 0
     integer, allocatable :: icell(:)
     real(dp), allocatable :: dunew(:,:)
     real(dp), allocatable :: ddivu(:)
     real(dp), allocatable :: denew(:)
  end type scatter_buf_t

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

end module hydro_scatter_commons
#ifdef HYDRO_CUDA
!###########################################################
! Module for hybrid CPU/GPU hydro dispatch data structures
!###########################################################
module hydro_hybrid_commons
  use hydro_scatter_commons
  use amr_parameters, only: dp
  implicit none

  integer, parameter :: HYBRID_SUPER_SIZE = 4096
  integer, parameter :: INTERP_INIT_CAP = 131072   ! 128K interpolated cells

  type gpu_state_t
     ! GPU-gather stencil index buffers (H2D per flush, lightweight)
     integer,  allocatable :: super_stencil_idx(:,:,:,:)  ! (SUPER_SIZE, iu1:iu2, ju1:ju2, ku1:ku2)
     integer,  allocatable :: super_stencil_grav(:,:,:,:) ! (SUPER_SIZE, iu1:iu2, ju1:ju2, ku1:ku2)
     real(dp), allocatable :: super_interp_vals(:,:)      ! (nvar, interp_cap) AoS for GPU
     integer :: super_n_interp = 0
     integer :: interp_cap = 0
     ! D2H compact output buffers (scatter-reduce results)
     real(dp), allocatable :: super_add_unew(:,:,:)       ! (SUPER_SIZE, 8, nvar)
     real(dp), allocatable :: super_add_lm1(:,:,:)        ! (SUPER_SIZE, 6, nvar)
     real(dp), allocatable :: super_add_divu_l(:,:)       ! (SUPER_SIZE, 8) [pressure_fix]
     real(dp), allocatable :: super_add_enew_l(:,:)       ! (SUPER_SIZE, 8)
     real(dp), allocatable :: super_add_divu_lm1(:,:)     ! (SUPER_SIZE, 6)
     real(dp), allocatable :: super_add_enew_lm1(:,:)     ! (SUPER_SIZE, 6)
     ! Grid info
     integer,  allocatable :: super_ind_grid(:)           ! (SUPER_SIZE)
     integer :: off = 0
  end type gpu_state_t

  type(scatter_buf_t), allocatable, save :: scatter_bufs(:)
  type(gpu_state_t), allocatable, save :: gpu_states(:)
  integer, save :: hybrid_max_threads = 0

  ! Timing accumulators
  real(dp), save :: acc_gather = 0, acc_gpu = 0, acc_scatter_l = 0, acc_merge = 0
  real(dp), save :: acc_mesh_upload = 0, acc_omp_total = 0
  integer, save :: n_gpu_flushes = 0, n_cpu_batches = 0, n_hybrid_calls = 0

contains

  subroutine pin_4d_int(arr)
    use iso_c_binding, only: c_loc, c_long_long
    use hydro_cuda_interface, only: cuda_host_register_c
    integer, intent(in), target :: arr(:,:,:,:)
    integer :: irc
    irc = cuda_host_register_c(c_loc(arr), int(size(arr)*4, c_long_long))
  end subroutine

  subroutine pin_3d_dp(arr)
    use iso_c_binding, only: c_loc, c_long_long
    use hydro_cuda_interface, only: cuda_host_register_c
    real(dp), intent(in), target :: arr(:,:,:)
    integer :: irc
    irc = cuda_host_register_c(c_loc(arr), int(size(arr)*8, c_long_long))
  end subroutine

  subroutine grow_interp_vals(gstate, nvar_in)
    type(gpu_state_t), intent(inout) :: gstate
    integer, intent(in) :: nvar_in
    real(dp), allocatable :: tmp(:,:)
    integer :: new_cap
    new_cap = gstate%interp_cap * 2
    allocate(tmp(nvar_in, new_cap))
    tmp(:, 1:gstate%interp_cap) = gstate%super_interp_vals(:, 1:gstate%interp_cap)
    call move_alloc(tmp, gstate%super_interp_vals)
    gstate%interp_cap = new_cap
  end subroutine grow_interp_vals

end module hydro_hybrid_commons
#endif
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine godunov_fine(ilevel)
  use amr_commons
  use hydro_commons
  use hydro_scatter_commons
#ifdef HYDRO_CUDA
  use cuda_commons
  use hydro_cuda_interface
#endif
!$ use omp_lib
  implicit none
  integer::ilevel
  !--------------------------------------------------------------------------
  ! This routine is a wrapper to the second order Godunov solver.
  ! Small grids (2x2x2) are gathered from level ilevel and sent to the
  ! hydro solver. On entry, hydro variables are gathered from array uold.
  ! On exit, unew has been updated.
  !--------------------------------------------------------------------------
  integer::i,ivar,igrid,ncache,ngrid,ic
  integer::nthreads,tid,it
  type(scatter_buf_t), allocatable :: cpu_scatter_bufs(:)
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

  ! Allocate per-thread scatter buffers for Level L-1
  nthreads = 1
  !$ nthreads = omp_get_max_threads()
  allocate(cpu_scatter_bufs(0:nthreads-1))
  do it = 0, nthreads-1
     call scatter_buf_init(cpu_scatter_bufs(it), SBUF_INIT_CAP, nvar)
  end do

  ! Dynamic scheduling with per-thread scatter buffers
  !$omp parallel do private(igrid,ngrid,tid) schedule(dynamic)
  do igrid=1,ncache,nvector
     tid = 0
     !$ tid = omp_get_thread_num()
     ngrid=MIN(nvector,ncache-igrid+1)
     call godfine1(ilevel, igrid, ngrid, cpu_scatter_bufs(tid))
  end do

  ! Serial merge: L-1 scatter buffers -> unew/divu/enew
  do it = 0, nthreads-1
     do i = 1, cpu_scatter_bufs(it)%count
        ic = cpu_scatter_bufs(it)%icell(i)
        do ivar = 1, nvar
           unew(ic, ivar) = unew(ic, ivar) + cpu_scatter_bufs(it)%dunew(i, ivar)
        end do
        if (pressure_fix) then
           divu(ic) = divu(ic) + cpu_scatter_bufs(it)%ddivu(i)
           enew(ic) = enew(ic) + cpu_scatter_bufs(it)%denew(i)
        end if
     end do
  end do

  ! Cleanup scatter buffers
  do it = 0, nthreads-1
     deallocate(cpu_scatter_bufs(it)%icell, cpu_scatter_bufs(it)%dunew, &
                cpu_scatter_bufs(it)%ddivu, cpu_scatter_bufs(it)%denew)
  end do
  deallocate(cpu_scatter_bufs)

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
subroutine godfine1(ilevel, jgrid, mgrid, sbuf)
  use amr_commons
  use hydro_commons
  use hydro_scatter_commons
  use poisson_commons
  use morton_hash
#ifdef HYDRO_CUDA
  use iso_c_binding, only: c_int, c_double
  use cuda_commons
  use hydro_cuda_interface
#endif
  implicit none
  integer, intent(in) :: ilevel, jgrid, mgrid
  type(scatter_buf_t), intent(inout), optional :: sbuf
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
  integer::icell_nbor,idx
  real(dp)::acc_unew(1:nvar),acc_ddivu,acc_denew
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
  ! (conflict-free: each thread writes to its own grid cells)
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
  if (present(sbuf)) then
     ! Scatter buffer path: accumulate per-grid L-1 contributions
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
  else
     ! Fallback: original critical section path
     !$omp critical (godunov)
     do idim=1,ndim
        i0=0; j0=0; k0=0
        if(idim==1)i0=1
        if(idim==2)j0=1
        if(idim==3)k0=1
        nb_noneigh=0
        do i=1,ngrid
           igridn_tmp = morton_nbor_grid(ind_grid(i),ilevel,2*idim-1)
           if (igridn_tmp==0) then
              nb_noneigh = nb_noneigh + 1
              ind_buffer(nb_noneigh) = morton_nbor_cell(ind_grid(i),ilevel,2*idim-1)
              ind_cell(nb_noneigh) = i
           end if
        end do
        nb_noneigh2=0
        do i=1,ngrid
           igridn_tmp = morton_nbor_grid(ind_grid(i),ilevel,2*idim)
           if (igridn_tmp==0) then
              nb_noneigh2 = nb_noneigh2 + 1
              ind_buffer2(nb_noneigh2) = morton_nbor_cell(ind_grid(i),ilevel,2*idim)
              ind_cell2(nb_noneigh2) = i
           end if
        end do
        do ivar=1,nvar
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
     !$omp end critical (godunov)
  end if

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
  use poisson_commons
  use hydro_hybrid_commons
  use iso_c_binding, only: c_int, c_long_long, c_loc
  use cuda_commons
  use hydro_cuda_interface
  implicit none
  integer, intent(in) :: ilevel, ncache
  !-------------------------------------------------------------------
  ! Dynamic hybrid CPU/GPU dispatcher for Godunov solver.
  ! GPU-gather path: upload mesh once, send stencil indices to GPU,
  ! GPU does gather + compute + scatter_reduce, D2H compact output.
  !-------------------------------------------------------------------
  integer :: tid, nthreads, stream_slot, igrid, ngrid
  integer :: it, i, ic, ivar, irc
  integer(kind=8) :: t_start, t_now, clock_rate
  integer(c_long_long) :: ncell_total
  logical :: mesh_on_gpu

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

  n_hybrid_calls = n_hybrid_calls + 1

  ! Reset scatter buffer counts (serial)
  do it = 0, hybrid_max_threads-1
     scatter_bufs(it)%count = 0
  end do

  ! Upload mesh arrays to GPU (serial, once per godunov_fine call)
  call system_clock(t_start)
  ncell_total = int(ncoarse, c_long_long) + int(twotondim, c_long_long) * int(ngridmax, c_long_long)
  call cuda_mesh_upload_f(uold, f, son, ncell_total, int(nvar, c_int), int(ndim, c_int), poisson)
  mesh_on_gpu = (cuda_mesh_is_ready_c() == 1)
  call system_clock(t_now)
  acc_mesh_upload = acc_mesh_upload + dble(t_now - t_start) / dble(clock_rate)

  call system_clock(t_start)  ! OMP total timer
  !$omp parallel private(tid, stream_slot, igrid, ngrid)
  tid = 0
  !$ tid = omp_get_thread_num()
  if (mesh_on_gpu) then
     stream_slot = int(cuda_acquire_stream_c())
  else
     stream_slot = -1  ! GPU mesh allocation failed, fall back to CPU
  end if

  ! Ensure GPU state allocated (lazy, per-stream-slot)
  if (stream_slot >= 0) then
     if (.not. allocated(gpu_states(stream_slot)%super_stencil_idx)) then
        allocate(gpu_states(stream_slot)%super_stencil_idx(1:HYBRID_SUPER_SIZE,iu1:iu2,ju1:ju2,ku1:ku2))
        allocate(gpu_states(stream_slot)%super_stencil_grav(1:HYBRID_SUPER_SIZE,iu1:iu2,ju1:ju2,ku1:ku2))
        allocate(gpu_states(stream_slot)%super_interp_vals(1:nvar,1:INTERP_INIT_CAP))
        gpu_states(stream_slot)%interp_cap = INTERP_INIT_CAP
        allocate(gpu_states(stream_slot)%super_add_unew(1:HYBRID_SUPER_SIZE,1:8,1:nvar))
        allocate(gpu_states(stream_slot)%super_add_lm1(1:HYBRID_SUPER_SIZE,1:6,1:nvar))
        allocate(gpu_states(stream_slot)%super_add_divu_l(1:HYBRID_SUPER_SIZE,1:8))
        allocate(gpu_states(stream_slot)%super_add_enew_l(1:HYBRID_SUPER_SIZE,1:8))
        allocate(gpu_states(stream_slot)%super_add_divu_lm1(1:HYBRID_SUPER_SIZE,1:6))
        allocate(gpu_states(stream_slot)%super_add_enew_lm1(1:HYBRID_SUPER_SIZE,1:6))
        allocate(gpu_states(stream_slot)%super_ind_grid(1:HYBRID_SUPER_SIZE))
        ! Pin stencil index buffers for fast H2D DMA
        call pin_4d_int(gpu_states(stream_slot)%super_stencil_idx)
        call pin_4d_int(gpu_states(stream_slot)%super_stencil_grav)
        ! Pin D2H buffers
        call pin_3d_dp(gpu_states(stream_slot)%super_add_unew)
        call pin_3d_dp(gpu_states(stream_slot)%super_add_lm1)
        ! Run bandwidth test once on first allocation
        !$omp critical (bwtest)
        if (n_hybrid_calls == 1) call cuda_h2d_bandwidth_test_c()
        !$omp end critical (bwtest)
     end if
     gpu_states(stream_slot)%off = 0
     gpu_states(stream_slot)%super_n_interp = 0
  end if

  !$omp do schedule(dynamic)
  do igrid = 1, ncache, nvector
     ngrid = MIN(nvector, ncache - igrid + 1)
     if (stream_slot >= 0) then
        call hybrid_gpu_stencil_batch(gpu_states(stream_slot), ilevel, igrid, ngrid, &
             scatter_bufs(tid), stream_slot)
     else
        call hybrid_cpu_process_batch(ilevel, igrid, ngrid, scatter_bufs(tid))
     end if
  end do
  !$omp end do nowait

  ! GPU threads: flush remaining super buffer
  if (stream_slot >= 0) then
     if (gpu_states(stream_slot)%off > 0) then
        call hybrid_gpu_gather_flush_scatter(gpu_states(stream_slot), stream_slot, &
             ilevel, scatter_bufs(tid))
     end if
     call cuda_release_stream_c(int(stream_slot, c_int))
  end if
  !$omp end parallel
  call system_clock(t_now)
  acc_omp_total = acc_omp_total + dble(t_now - t_start) / dble(clock_rate)

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
     write(*,'(A,F8.3,A)') '   mesh upload : ', acc_mesh_upload, ' s'
     write(*,'(A,F8.3,A)') '   OMP total   : ', acc_omp_total, ' s'
     write(*,'(A,F8.3,A)') '   gather      : ', acc_gather, ' s'
     write(*,'(A,F8.3,A)') '   GPU compute : ', acc_gpu, ' s'
     write(*,'(A,F8.3,A)') '   L scatter   : ', acc_scatter_l, ' s'
     write(*,'(A,F8.3,A)') '   L-1 merge   : ', acc_merge, ' s'
     write(*,'(A,I6)') '   hybrid calls: ', n_hybrid_calls
     write(*,'(A,I6)') '   GPU flushes : ', n_gpu_flushes
     write(*,'(A,I6)') '   CPU batches : ', n_cpu_batches
     call hydro_cuda_profile_report_c()
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
subroutine hybrid_gpu_stencil_batch(gstate, ilevel, igrid_start, ngrid, sbuf, stream_slot)
  use amr_commons
  use hydro_commons
  use hydro_parameters
  use poisson_commons
  use hydro_hybrid_commons
  implicit none
  type(gpu_state_t), intent(inout) :: gstate
  integer, intent(in) :: ilevel, igrid_start, ngrid, stream_slot
  type(scatter_buf_t), intent(inout) :: sbuf
  !-------------------------------------------------------------------
  ! GPU-gather path: store stencil INDICES into super buffer instead
  ! of copying data. GPU will gather from mesh using these indices.
  ! For interpolated (coarse) cells, compute interpol_hydro on CPU and
  ! store result in interp_vals buffer, referenced by negative index.
  !-------------------------------------------------------------------
  integer ,dimension(1:nvector)::ind_grid_loc,igrid_nbor,ind_cell
  integer ,dimension(1:nvector)::ind_exist,ind_nexist,ind_buffer
  integer ,dimension(1:nvector,1:threetondim)::nbors_father_cells
  integer ,dimension(1:nvector,1:twotondim)::nbors_father_grids
  integer ,dimension(1:nvector,0:twondim)::ibuffer_father
  real(dp),dimension(1:nvector,0:twondim,1:nvar)::u1
  real(dp),dimension(1:nvector,1:twotondim,1:nvar)::u2

  integer :: off, i, j, ivar, ind_son, iskip, ind_father, nbuffer, nexist
  integer :: i1,j1,k1,i2,j2,k2,i3,j3,k3, cell_idx, n_interp_slot
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

  ! Initialize stencil_grav to zero (default: no gravity)
  gstate%super_stencil_grav(off+1:off+ngrid,iu1:iu2,ju1:ju2,ku1:ku2) = 0

  ! Gather father cells
  do i=1,ngrid
     ind_cell(i) = father(ind_grid_loc(i))
  end do
  call get3cubefather(ind_cell,nbors_father_cells,nbors_father_grids,ngrid,ilevel)

  ! Build stencil index arrays for 6x6x6 cells
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

     ! Interpolate hydro variables for coarse cells (same as before)
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

        i3=1; j3=1; k3=1
        if(ndim>0)i3=1+2*(i1-1)+i2
        if(ndim>1)j3=1+2*(j1-1)+j2
        if(ndim>2)k3=1+2*(k1-1)+k2

        ! Existing fine grids: store cell index (positive)
        do i=1,nexist
           cell_idx = iskip + igrid_nbor(ind_exist(i))
           gstate%super_stencil_idx(off+ind_exist(i),i3,j3,k3) = cell_idx
           if(poisson) gstate%super_stencil_grav(off+ind_exist(i),i3,j3,k3) = cell_idx
        end do

        ! Interpolated coarse cells: store interp values, negative index
        do i=1,nbuffer
           gstate%super_n_interp = gstate%super_n_interp + 1
           n_interp_slot = gstate%super_n_interp
           if (n_interp_slot > gstate%interp_cap) call grow_interp_vals(gstate, nvar)
           gstate%super_interp_vals(1:nvar, n_interp_slot) = u2(i, ind_son, 1:nvar)
           gstate%super_stencil_idx(off+ind_nexist(i),i3,j3,k3) = -n_interp_slot
           ! Gravity: straight injection from father cell
           if(poisson) gstate%super_stencil_grav(off+ind_nexist(i),i3,j3,k3) = ibuffer_father(i,0)
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
     call hybrid_gpu_gather_flush_scatter(gstate, stream_slot, ilevel, sbuf)
  end if

end subroutine hybrid_gpu_stencil_batch
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine hybrid_gpu_gather_flush_scatter(gstate, stream_slot, ilevel, sbuf)
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
  ! GPU-gather + scatter-reduce pipeline for accumulated super buffer.
  ! H2D: stencil indices (~7 MB) + interp_vals (tiny)
  ! GPU: gather from mesh → compute → scatter_reduce
  ! D2H: compact add_unew/add_lm1 (~5 MB)
  ! CPU applies results to unew/sbuf.
  !-------------------------------------------------------------------
  integer :: super_total, g, ig, ic, igridn_tmp, icell_nbor, idx
  integer :: ind_son, iskip, ivar, face
  integer :: nx_loc
  real(dp) :: dx, scale, dt_val
  integer(kind=8) :: t_start, t_now, clock_rate

  call system_clock(count_rate=clock_rate)

  super_total = gstate%off
  if (super_total == 0) return

  nx_loc = icoarse_max-icoarse_min+1
  scale = boxlen/dble(nx_loc)
  dx = 0.5D0**ilevel*scale

  ! GPU compute: H2D stencil idx → gather → 5 kernels + scatter_reduce → D2H compact
  call system_clock(t_start)
  call hydro_cuda_gather_reduce_async_f( &
       gstate%super_stencil_idx, gstate%super_stencil_grav, &
       gstate%super_interp_vals, &
       gstate%super_add_unew, gstate%super_add_lm1, &
       gstate%super_add_divu_l, gstate%super_add_enew_l, &
       gstate%super_add_divu_lm1, gstate%super_add_enew_lm1, &
       real(dx, c_double), real(dx, c_double), real(dx, c_double), &
       real(dtnew(ilevel), c_double), &
       int(super_total, c_int), int(HYBRID_SUPER_SIZE, c_int), &
       int(gstate%super_n_interp, c_int), int(stream_slot, c_int), pressure_fix)
  call hydro_cuda_gather_reduce_sync_f( &
       int(super_total, c_int), int(stream_slot, c_int))
  call hydro_cuda_profile_accumulate_c(int(stream_slot, c_int), int(super_total, c_int))
  call system_clock(t_now)
  dt_val = dble(t_now - t_start) / dble(clock_rate)
  !$omp atomic
  acc_gpu = acc_gpu + dt_val

  ! Apply compact results to unew
  call system_clock(t_start)

  ! Level L scatter (direct, conflict-free): apply add_unew to unew
  do g = 1, super_total
     ig = gstate%super_ind_grid(g)
     do ind_son = 1, 8
        iskip = ncoarse + (ind_son-1)*ngridmax
        ic = iskip + ig
        do ivar = 1, nvar
           unew(ic, ivar) = unew(ic, ivar) + gstate%super_add_unew(g, ind_son, ivar)
        end do
        if (pressure_fix) then
           divu(ic) = divu(ic) + gstate%super_add_divu_l(g, ind_son)
           enew(ic) = enew(ic) + gstate%super_add_enew_l(g, ind_son)
        end if
     end do
  end do

  ! Level L-1 scatter: apply add_lm1 entries to per-thread buffer
  ! face 1=-x, 2=+x, 3=-y, 4=+y, 5=-z, 6=+z (RAMSES convention)
  do g = 1, super_total
     ig = gstate%super_ind_grid(g)
     do face = 1, 2*ndim
        igridn_tmp = morton_nbor_grid(ig, ilevel, face)
        if (igridn_tmp == 0) then
           icell_nbor = morton_nbor_cell(ig, ilevel, face)
           sbuf%count = sbuf%count + 1
           idx = sbuf%count
           if (idx > sbuf%capacity) call scatter_buf_grow(sbuf)
           sbuf%icell(idx) = icell_nbor
           do ivar = 1, nvar
              sbuf%dunew(idx, ivar) = gstate%super_add_lm1(g, face, ivar)
           end do
           if (pressure_fix) then
              sbuf%ddivu(idx) = gstate%super_add_divu_lm1(g, face)
              sbuf%denew(idx) = gstate%super_add_enew_lm1(g, face)
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
  gstate%super_n_interp = 0

end subroutine hybrid_gpu_gather_flush_scatter
#endif
