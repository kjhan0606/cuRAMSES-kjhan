!###########################################################
!###########################################################
! Dynamic OMP/CUDA hybrid dispatch for upload_fine
! Overrides patch/Horizon5-master-2/interpol_hydro.kjhan.f90
!###########################################################
!###########################################################

#ifdef HYDRO_CUDA
!=========================================================
! Module: GPU state for upload_fine superbatch
!=========================================================
module upload_hybrid_commons
  use amr_parameters, only: dp, nvector, ndim, twotondim, ngridmax
  implicit none

  integer, parameter :: UPLOAD_SUPER_SIZE = 8192

  type upload_gpu_state_t
     ! child_buf(nvar, 8, cap) — uold at 8 children per parent
     real(dp), allocatable :: child_buf(:,:,:)
     ! parent_buf(nvar+1, cap) — averaged parent values + avg internal energy
     real(dp), allocatable :: parent_buf(:,:)
     ! parent_idx(cap) — global cell index for scatter
     integer, allocatable :: parent_idx(:)
     integer :: nsplit = 0   ! number of split cells accumulated
     integer :: cap = 0      ! current capacity
     integer :: nv = 0       ! nvar used for allocation
  end type

  type(upload_gpu_state_t), allocatable, save, target :: upload_gstates(:)
  integer, save :: upload_hybrid_inited = 0

contains

  subroutine upload_gstate_ensure(gs, needed, nvar_in)
    type(upload_gpu_state_t), intent(inout) :: gs
    integer, intent(in) :: needed, nvar_in
    if (needed <= gs%cap .and. nvar_in == gs%nv) return
    gs%cap = max(needed, UPLOAD_SUPER_SIZE)
    gs%nv = nvar_in
    if (allocated(gs%child_buf)) deallocate(gs%child_buf)
    if (allocated(gs%parent_buf)) deallocate(gs%parent_buf)
    if (allocated(gs%parent_idx)) deallocate(gs%parent_idx)
    allocate(gs%child_buf(nvar_in, 8, gs%cap))
    allocate(gs%parent_buf(nvar_in + 1, gs%cap))
    allocate(gs%parent_idx(gs%cap))
  end subroutine

end module upload_hybrid_commons
#endif

!###########################################################
!###########################################################
subroutine upload_fine(ilevel)
  use amr_commons
  use hydro_commons
  implicit none
  integer::ilevel
  !----------------------------------------------------------------------
  ! This routine performs a restriction operation (averaging down)
  ! for the hydro variables.
  !----------------------------------------------------------------------
  integer::i,ncache,igrid,ngrid

  if(ilevel==nlevelmax)return
  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel

  ncache=active(ilevel)%ngrid

#ifdef HYDRO_CUDA
  call upload_fine_hybrid(ilevel, ncache)
#else
!$omp parallel do private(igrid,ngrid) schedule(dynamic)
  do igrid=1,ncache,nvector
     ngrid=MIN(nvector,ncache-igrid+1)
     call sub_upload_fine(ilevel, igrid, ngrid)
  end do
#endif

111 format('   Entering upload_fine for level',i2)
end subroutine upload_fine

#ifdef HYDRO_CUDA
!###########################################################
!###########################################################
! Hybrid dispatch for upload_fine
!###########################################################
!###########################################################
subroutine upload_fine_hybrid(ilevel, ncache)
  use amr_commons
  use hydro_commons
  use upload_hybrid_commons
  use cuda_commons
  use hydro_cuda_interface
  use iso_c_binding
  implicit none
  integer, intent(in) :: ilevel, ncache

  integer :: igrid, ngrid, stream_slot
  type(upload_gpu_state_t), pointer :: gs

  ! First-call: allocate GPU states
  if (upload_hybrid_inited == 0) then
     allocate(upload_gstates(0:7))
     upload_hybrid_inited = 1
  end if

  !$omp parallel private(igrid, ngrid, stream_slot, gs)
  stream_slot = cuda_acquire_stream_c()
  if (stream_slot >= 0) then
     gs => upload_gstates(stream_slot)
     call upload_gstate_ensure(gs, UPLOAD_SUPER_SIZE, nvar)
     gs%nsplit = 0
  end if

  !$omp do schedule(dynamic)
  do igrid = 1, ncache, nvector
     ngrid = MIN(nvector, ncache - igrid + 1)
     if (stream_slot >= 0) then
        call upload_gpu_gather_batch(gs, ilevel, igrid, ngrid, stream_slot)
     else
        call sub_upload_fine(ilevel, igrid, ngrid)
     end if
  end do
  !$omp end do nowait

  if (stream_slot >= 0) then
     if (gs%nsplit > 0) call upload_gpu_flush_scatter(gs, stream_slot)
     call cuda_release_stream_c(stream_slot)
  end if
  !$omp end parallel

end subroutine upload_fine_hybrid

!###########################################################
! GPU gather: collect child cell data for superbatch
!###########################################################
subroutine upload_gpu_gather_batch(gs, ilevel, igrid_start, ngrid, stream_slot)
  use amr_commons
  use hydro_commons
  use upload_hybrid_commons
  implicit none

  type(upload_gpu_state_t), intent(inout) :: gs
  integer, intent(in) :: ilevel, igrid_start, ngrid, stream_slot

  integer :: i, ind, ind_son, iskip, iskip_son, ivar
  integer, dimension(1:nvector) :: ind_grid, ind_cell
  integer :: igrid_son, ns

  ! Load grid indices
  do i = 1, ngrid
     ind_grid(i) = active(ilevel)%igrid(igrid_start + i - 1)
  end do

  ! Loop over cells in each grid
  do ind = 1, twotondim
     iskip = ncoarse + (ind - 1) * ngridmax
     do i = 1, ngrid
        ind_cell(i) = iskip + ind_grid(i)
     end do

     ! Process split cells (cells with children)
     do i = 1, ngrid
        if (son(ind_cell(i)) > 0) then
           ns = gs%nsplit + 1
           gs%parent_idx(ns) = ind_cell(i)
           igrid_son = son(ind_cell(i))

           ! Gather 8 children's uold values
           do ind_son = 1, twotondim
              iskip_son = ncoarse + (ind_son - 1) * ngridmax
              do ivar = 1, nvar
                 gs%child_buf(ivar, ind_son, ns) = &
                      uold(iskip_son + igrid_son, ivar)
              end do
           end do

           gs%nsplit = ns

           ! Flush if buffer is near full
           if (gs%nsplit + nvector * twotondim > gs%cap) then
              call upload_gpu_flush_scatter(gs, stream_slot)
           end if
        end if
     end do
  end do

end subroutine upload_gpu_gather_batch

!###########################################################
! GPU flush: launch kernel, scatter results
!###########################################################
subroutine upload_gpu_flush_scatter(gs, stream_slot)
  use amr_commons
  use hydro_commons
  use hydro_parameters, only: smallr
  use upload_hybrid_commons
  use hydro_cuda_interface
  use iso_c_binding
  implicit none

  type(upload_gpu_state_t), intent(inout) :: gs
  integer, intent(in) :: stream_slot

  integer :: p, ivar, idim, do_eint
  real(dp) :: ekin, erad, avg_eint
#if NENER>0
  integer :: irad
#endif

  if (gs%nsplit == 0) return

  ! Determine if internal energy correction is needed
  do_eint = 0
  if (interpol_var == 1 .or. interpol_var == 2) do_eint = 1

  ! Launch GPU kernel
  call upload_fine_cuda_async_f(gs%child_buf, gs%parent_buf, &
       int(gs%nsplit, c_int), int(nvar, c_int), int(ndim, c_int), &
       real(smallr, c_double), int(do_eint, c_int), int(stream_slot, c_int))
  call upload_fine_cuda_sync_f(int(stream_slot, c_int))

  ! Apply internal energy correction on CPU if needed
  if (do_eint == 1) then
     do p = 1, gs%nsplit
        avg_eint = gs%parent_buf(nvar + 1, p)
        ! Compute new kinetic energy from averaged momenta
        ekin = 0.0d0
        do idim = 1, ndim
           ekin = ekin + 0.5d0 * gs%parent_buf(1+idim, p)**2 &
                / max(gs%parent_buf(1, p), smallr)
        end do
        ! Compute new radiative energy
        erad = 0.0d0
#if NENER>0
        do irad = 1, nener
           erad = erad + gs%parent_buf(ndim+2+irad, p)
        end do
#endif
        ! Correct total energy
        gs%parent_buf(ndim + 2, p) = avg_eint + ekin + erad
     end do
  end if

  ! Scatter averaged values to global uold array
  do p = 1, gs%nsplit
     do ivar = 1, nvar
        uold(gs%parent_idx(p), ivar) = gs%parent_buf(ivar, p)
     end do
  end do

  gs%nsplit = 0

end subroutine upload_gpu_flush_scatter
#endif

!##########################################################################
!##########################################################################
! sub_upload_fine (CPU fallback, unchanged)
!##########################################################################
!##########################################################################
subroutine sub_upload_fine(ilevel,igrid,ngrid)
  use amr_commons
  use hydro_commons
  implicit none
  integer::ilevel
  integer::i,ncache,igrid,ngrid,ind,iskip,nsplit,icell
  integer,dimension(1:nvector)::ind_grid,ind_cell,ind_split
  logical,dimension(1:nvector)::ok
     do i=1,ngrid
        ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
     end do

     ! Loop over cells
     do ind=1,twotondim
        iskip=ncoarse+(ind-1)*ngridmax
        do i=1,ngrid
           ind_cell(i)=iskip+ind_grid(i)
        end do

        ! Gather split cells
        do i=1,ngrid
           ok(i)=son(ind_cell(i))>0
        end do

        ! Count split cells
        nsplit=0
        do i=1,ngrid
           if(ok(i))nsplit=nsplit+1
        end do

        ! Upload for selected cells
        if(nsplit>0)then
           icell=0
           do i=1,ngrid
              if(ok(i))then
                 icell=icell+1
                 ind_split(icell)=ind_cell(i)
              end if
           end do
           call upl(ind_split,nsplit)
        end if

     end do

end subroutine sub_upload_fine
!##########################################################################
!##########################################################################
subroutine upl(ind_cell,ncell)
  use amr_commons
  use hydro_commons
  implicit none
  integer::ncell
  integer,dimension(1:nvector)::ind_cell
  !---------------------------------------------------------------------
  ! This routine performs a restriction operation (averaging down)
  !---------------------------------------------------------------------
  integer ::ivar,irad,i,idim,ind_son,iskip_son
  integer ,dimension(1:nvector)::igrid_son,ind_cell_son
  real(dp),dimension(1:nvector)::getx,ekin,erad

  ! Get child oct index
  do i=1,ncell
     igrid_son(i)=son(ind_cell(i))
  end do

  ! Average conservative variables
  do ivar=1,nvar

     getx(1:ncell)=0.0d0
     do ind_son=1,twotondim
        iskip_son=ncoarse+(ind_son-1)*ngridmax
        do i=1,ncell
           ind_cell_son(i)=iskip_son+igrid_son(i)
        end do
        do i=1,ncell
           getx(i)=getx(i)+uold(ind_cell_son(i),ivar)
        end do
     end do

     do i=1,ncell
        uold(ind_cell(i),ivar)=getx(i)/dble(twotondim)
     end do

  end do

  ! Average internal energy instead of total energy
  if(interpol_var==1 .or. interpol_var==2)then

     getx(1:ncell)=0.0d0
     do ind_son=1,twotondim
        iskip_son=ncoarse+(ind_son-1)*ngridmax
        do i=1,ncell
           ind_cell_son(i)=iskip_son+igrid_son(i)
        end do
        ekin(1:ncell)=0.0d0
        do idim=1,ndim
           do i=1,ncell
              ekin(i)=ekin(i)+0.5d0*uold(ind_cell_son(i),1+idim)**2 &
                   &               /max(uold(ind_cell_son(i),1),smallr)
           end do
        end do
        erad(1:ncell)=0.0d0
#if NENER>0
        do irad=1,nener
           do i=1,ncell
              erad(i)=erad(i)+uold(ind_cell_son(i),ndim+2+irad)
           end do
        end do
#endif
        do i=1,ncell
           getx(i)=getx(i)+uold(ind_cell_son(i),ndim+2)-ekin(i)-erad(i)
        end do
     end do

     ekin(1:ncell)=0.0d0
     do idim=1,ndim
        do i=1,ncell
           ekin(i)=ekin(i)+0.5d0*uold(ind_cell(i),1+idim)**2 &
                &               /max(uold(ind_cell(i),1),smallr)
        end do
     end do
     erad(1:ncell)=0.0d0
#if NENER>0
     do irad=1,nener
        do i=1,ncell
           erad(i)=erad(i)+uold(ind_cell(i),ndim+2+irad)
        end do
     end do
#endif

     do i=1,ncell
        uold(ind_cell(i),ndim+2)=getx(i)/dble(twotondim)+ekin(i)+erad(i)
     end do

  end if

end subroutine upl
!###########################################################
!###########################################################
subroutine interpol_hydro(u1,u2,nn)
  use amr_commons
  use hydro_commons
  use poisson_commons
  implicit none
  integer::nn
  real(dp),dimension(1:nvector,0:twondim  ,1:nvar)::u1
  real(dp),dimension(1:nvector,1:twotondim,1:nvar)::u2
  !----------------------------------------------------------
  ! This routine performs a prolongation (interpolation)
  ! operation for newly refined cells or buffer cells.
  !----------------------------------------------------------
  integer::i,j,ivar,irad,idim,ind,ix,iy,iz,ind2
  real(dp)::oneover_twotondim
  real(dp),dimension(1:twotondim,1:3)::xc
  real(dp),dimension(1:nvector,0:twondim)::a
  real(dp),dimension(1:nvector,1:ndim)::w
  real(dp),dimension(1:nvector)::ekin,mom
  real(dp),dimension(1:nvector)::erad

  oneover_twotondim=1.D0/dble(twotondim)

  do ind=1,twotondim
     iz=(ind-1)/4
     iy=(ind-1-4*iz)/2
     ix=(ind-1-2*iy-4*iz)
     if(ndim>0)xc(ind,1)=(dble(ix)-0.5D0)
     if(ndim>1)xc(ind,2)=(dble(iy)-0.5D0)
     if(ndim>2)xc(ind,3)=(dble(iz)-0.5D0)
  end do

  if(interpol_var==1 .or. interpol_var==2)then
     do j=0,twondim
        ekin(1:nn)=0.0d0
        do idim=1,ndim
           do i=1,nn
              ekin(i)=ekin(i)+0.5d0*u1(i,j,idim+1)**2/max(u1(i,j,1),smallr)
           end do
        end do
        erad(1:nn)=0.0d0
#if NENER>0
        do irad=1,nener
           do i=1,nn
              erad(i)=erad(i)+u1(i,j,ndim+2+irad)
           end do
        end do
#endif
        do i=1,nn
           u1(i,j,ndim+2)=u1(i,j,ndim+2)-ekin(i)-erad(i)
        end do
        if(interpol_var==2)then
           do idim=1,ndim
              do i=1,nn
                 u1(i,j,idim+1)=u1(i,j,idim+1)/max(u1(i,j,1),smallr)
              end do
           end do
        end if
     end do
  end if

  do ivar=1,nvar
     do j=0,twondim
        do i=1,nn
           a(i,j)=u1(i,j,ivar)
        end do
     end do
     w(1:nn,1:ndim)=0.0D0
     if(interpol_type==1)call compute_limiter_minmod(a,w,nn)
     if(interpol_type==2)call compute_limiter_central(a,w,nn)
     if(interpol_type==3)call compute_central(a,w,nn)
     if(interpol_type==4)then
        if (interpol_var .ne. 2)then
           write(*,*)'interpol_type=4 is designed for interpol_var=2'
           call clean_stop
        end if
        if (ivar>1 .and. (ivar <= 1+ndim))then
           call compute_central(a,w,nn)
        else
           call compute_limiter_central(a,w,nn)
        end if
     end if
     do ind=1,twotondim
        u2(1:nn,ind,ivar)=a(1:nn,0)
        do idim=1,ndim
           do i=1,nn
              u2(i,ind,ivar)=u2(i,ind,ivar)+w(i,idim)*xc(ind,idim)
           end do
        end do
     end do
  end do

  if(interpol_var==1 .or. interpol_var==2)then
     if(interpol_var==2)then
        do ind=1,twotondim
           do idim=1,ndim
              do i=1,nn
                 u2(i,ind,idim+1)=u2(i,ind,idim+1)*u2(i,ind,1)
              end do
           end do
        end do
        do idim=1,ndim
           mom(1:nn)=0.
           do ind=1,twotondim
              do i=1,nn
                 mom(i)=mom(i)+u2(i,ind,idim+1)*oneover_twotondim
              end do
           end do
           do i=1,nn
              mom(i)=mom(i)-u1(i,0,idim+1)*u1(i,0,1)
              u2(i,1:twotondim,idim+1)=u2(i,1:twotondim,idim+1)-mom(i)
           end do
        end do
     end if
     do ind=1,twotondim
        ekin(1:nn)=0.0d0
        do idim=1,ndim
           do i=1,nn
              ekin(i)=ekin(i)+0.5d0*u2(i,ind,idim+1)**2/max(u2(i,ind,1),smallr)
           end do
        end do
        erad(1:nn)=0.0d0
#if NENER>0
        do irad=1,nener
           do i=1,nn
              erad(i)=erad(i)+u2(i,ind,ndim+2+irad)
           end do
        end do
#endif
        do i=1,nn
           u2(i,ind,ndim+2)=u2(i,ind,ndim+2)+ekin(i)+erad(i)
        end do
     end do
  end if

end subroutine interpol_hydro
!###########################################################
!###########################################################
subroutine compute_limiter_minmod(a,w,nn)
  use amr_commons
  use hydro_commons
  implicit none
  integer::nn
  real(dp),dimension(1:nvector,0:twondim)::a
  real(dp),dimension(1:nvector,1:ndim)::w
  integer::i,idim
  real(dp)::diff_left,diff_right,minmod

  do idim=1,ndim
     do i=1,nn
        diff_left=0.5*(a(i,2*idim)-a(i,0))
        diff_right=0.5*(a(i,0)-a(i,2*idim-1))
        if(diff_left*diff_right<=0.0)then
           minmod=0.0
        else
           minmod=MIN(ABS(diff_left),ABS(diff_right)) &
                &   *diff_left/ABS(diff_left)
        end if
        w(i,idim)=minmod
     end do
  end do

end subroutine compute_limiter_minmod
!###########################################################
!###########################################################
subroutine compute_limiter_central(a,w,nn)
  use amr_commons
  use hydro_commons
  implicit none
  integer::nn
  real(dp),dimension(1:nvector,0:twondim)::a
  real(dp),dimension(1:nvector,1:ndim)::w
  integer::i,j,idim,ind,ix,iy,iz
  real(dp),dimension(1:twotondim,1:3)::xc
  real(dp)::xxc
  real(dp),dimension(1:nvector,1:twotondim)::ac
  real(dp),dimension(1:nvector)::corner,kernel,diff_corner,diff_kernel
  real(dp),dimension(1:nvector)::max_limiter,min_limiter,limiter

  do ind=1,twotondim
     iz=(ind-1)/4
     iy=(ind-1-4*iz)/2
     ix=(ind-1-2*iy-4*iz)
     if(ndim>0)xc(ind,1)=(dble(ix)-0.5D0)
     if(ndim>1)xc(ind,2)=(dble(iy)-0.5D0)
     if(ndim>2)xc(ind,3)=(dble(iz)-0.5D0)
  end do

  do idim=1,ndim
     do i=1,nn
        w(i,idim)=0.25D0*(a(i,2*idim)-a(i,2*idim-1))
     end do
  end do
  do ind=1,twotondim
     do i=1,nn
        ac(i,ind)=a(i,0)
     end do
  end do
  do idim=1,ndim
     do ind=1,twotondim
        xxc = xc(ind,idim)
        do i=1,nn
           corner(i)=ac(i,ind)+2.D0*w(i,idim)*xxc
        end do
        do i=1,nn
           ac(i,ind)=corner(i)
        end do
     end do
  end do
  do i=1,nn
     corner(i)=ac(i,1)
  end do
  do j=2,twotondim
     do i=1,nn
        corner(i)=MAX(corner(i),ac(i,j))
     end do
  end do
  do i=1,nn
     kernel(i)=a(i,1)
  end do
  do j=2,twondim
     do i=1,nn
        kernel(i)=MAX(kernel(i),a(i,j))
     end do
  end do
  do i=1,nn
     diff_kernel(i)=a(i,0)-kernel(i)
     diff_corner(i)=a(i,0)-corner(i)
  end do
  max_limiter=0.0D0
  do i=1,nn
     if(diff_kernel(i)*diff_corner(i) > 0.0D0)then
        max_limiter(i)=MIN(1.0_dp,diff_kernel(i)/diff_corner(i))
     end if
  end do
  do i=1,nn
     corner(i)=ac(i,1)
  end do
  do j=2,twotondim
     do i=1,nn
        corner(i)=MIN(corner(i),ac(i,j))
     end do
  end do
  do i=1,nn
     kernel(i)=a(i,1)
  end do
  do j=2,twondim
     do i=1,nn
        kernel(i)=MIN(kernel(i),a(i,j))
     end do
  end do
  do i=1,nn
     diff_kernel(i)=a(i,0)-kernel(i)
     diff_corner(i)=a(i,0)-corner(i)
  end do
  min_limiter=0.0D0
  do i=1,nn
     if(diff_kernel(i)*diff_corner(i) > 0.0D0)then
        min_limiter(i)=MIN(1.0_dp,diff_kernel(i)/diff_corner(i))
     end if
  end do
  do i=1,nn
     limiter(i)=MIN(min_limiter(i),max_limiter(i))
  end do
  do idim=1,ndim
     do i=1,nn
        w(i,idim)=w(i,idim)*limiter(i)
     end do
  end do

end subroutine compute_limiter_central
!###########################################################
!###########################################################
subroutine compute_central(a,w,nn)
  use amr_commons
  use hydro_commons
  implicit none
  integer::nn
  real(dp),dimension(1:nvector,0:twondim)::a
  real(dp),dimension(1:nvector,1:ndim)::w
  integer::i,idim

  do idim=1,ndim
     do i=1,nn
        w(i,idim)=0.25D0*(a(i,2*idim)-a(i,2*idim-1))
     end do
  end do

end subroutine compute_central
