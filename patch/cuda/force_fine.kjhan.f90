!#########################################################
!#########################################################
! Dynamic OMP/CUDA hybrid dispatch for force_fine
! Overrides patch/Horizon5-master-2/force_fine.kjhan.f90
!#########################################################
!#########################################################

#ifdef HYDRO_CUDA
!=========================================================
! Module: GPU state for gradient_phi superbatch
!=========================================================
module force_hybrid_commons
  use amr_parameters, only: dp, nvector, ndim, twotondim, twondim
  implicit none

  integer, parameter :: FORCE_SUPER_SIZE = 4096

  type force_gpu_state_t
     ! phi stencil buffer: phi_buf(4, cap*24)
     ! Layout: phi_buf(1:4, slot) where slot=(g-1)*24+(ind-1)*3+idim
     real(dp), allocatable :: phi_buf(:,:)
     ! Force output buffer: f_buf(cap*24)
     real(dp), allocatable :: f_buf(:)
     ! Cell index buffer for scatter: cell_buf(cap, 8)
     integer, allocatable :: cell_buf(:,:)
     integer :: off = 0     ! number of grids accumulated
     integer :: cap = 0     ! current capacity
  end type

  type(force_gpu_state_t), allocatable, save, target :: force_gstates(:)
  integer, save :: force_hybrid_inited = 0

contains

  subroutine force_gstate_ensure(gs, needed)
    type(force_gpu_state_t), intent(inout) :: gs
    integer, intent(in) :: needed
    if (needed <= gs%cap) return
    gs%cap = max(needed, FORCE_SUPER_SIZE)
    if (allocated(gs%phi_buf)) deallocate(gs%phi_buf)
    if (allocated(gs%f_buf)) deallocate(gs%f_buf)
    if (allocated(gs%cell_buf)) deallocate(gs%cell_buf)
    allocate(gs%phi_buf(4, gs%cap * 24))
    allocate(gs%f_buf(gs%cap * 24))
    allocate(gs%cell_buf(gs%cap, 8))
  end subroutine

end module force_hybrid_commons
#endif

!#########################################################
!#########################################################
subroutine force_fine(ilevel,icount)
  use amr_commons
  use pm_commons
  use poisson_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer::ilevel,icount
  !----------------------------------------------------------
  ! This routine computes the gravitational acceleration,
  ! the maximum density rho_max, and the potential energy
  !----------------------------------------------------------
  integer::igrid,ngrid,ncache,i,ind,iskip,ix,iy,iz
  integer::info,ibound,nx_loc,idim
  real(dp)::dx,dx_loc,scale,fact,fourpi
  real(kind=8)::rho_loc,rho_all,epot_loc,epot_all
  real(dp),dimension(1:twotondim,1:3)::xc
  real(dp),dimension(1:3)::skip_loc

  ! Work arrays (thread-private via OMP private clause)
  integer ,dimension(1:nvector)::ind_grid_w,ind_cell_w
  real(dp),dimension(1:nvector,1:ndim)::xx_w,ff_w

  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel

  ! Mesh size at level ilevel in coarse cell units
  dx=0.5D0**ilevel

  ! Rescaling factors
  nx_loc=(icoarse_max-icoarse_min+1)
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  if(ndim>0)skip_loc(1)=dble(icoarse_min)
  if(ndim>1)skip_loc(2)=dble(jcoarse_min)
  if(ndim>2)skip_loc(3)=dble(kcoarse_min)
  scale=boxlen/dble(nx_loc)
  dx_loc=dx*scale

  ! Set position of cell centers relative to grid center
  do ind=1,twotondim
     iz=(ind-1)/4
     iy=(ind-1-4*iz)/2
     ix=(ind-1-2*iy-4*iz)
     if(ndim>0)xc(ind,1)=(dble(ix)-0.5D0)*dx
     if(ndim>1)xc(ind,2)=(dble(iy)-0.5D0)*dx
     if(ndim>2)xc(ind,3)=(dble(iz)-0.5D0)*dx
  end do

  !-------------------------------------
  ! Compute analytical gravity force
  !-------------------------------------
  if(gravity_type>0)then

     ncache=active(ilevel)%ngrid
!$omp parallel do private(igrid,ngrid,i,ind,iskip,idim,ind_grid_w,ind_cell_w,xx_w,ff_w) schedule(dynamic)
     do igrid=1,ncache,nvector
        ngrid=MIN(nvector,ncache-igrid+1)
        do i=1,ngrid
           ind_grid_w(i)=active(ilevel)%igrid(igrid+i-1)
        end do
        do ind=1,twotondim
           iskip=ncoarse+(ind-1)*ngridmax
           do i=1,ngrid
              ind_cell_w(i)=iskip+ind_grid_w(i)
           end do
           do idim=1,ndim
              do i=1,ngrid
                 xx_w(i,idim)=xg(ind_grid_w(i),idim)+xc(ind,idim)
              end do
           end do
           do idim=1,ndim
              do i=1,ngrid
                 xx_w(i,idim)=(xx_w(i,idim)-skip_loc(idim))*scale
              end do
           end do
           call gravana(xx_w,ff_w,dx_loc,ngrid)
           do idim=1,ndim
              do i=1,ngrid
                 f(ind_cell_w(i),idim)=ff_w(i,idim)
              end do
           end do
        end do
     end do

     do idim=1,ndim
        call make_virtual_fine_dp(f(1,idim),ilevel)
     end do
     if(simple_boundary)call make_boundary_force(ilevel)

  !------------------------------
  ! Compute gradient of potential
  !------------------------------
  else
     call make_boundary_phi(ilevel)

     ncache=active(ilevel)%ngrid
#ifdef HYDRO_CUDA
     call force_gradient_hybrid(ilevel, icount, ncache)
#else
!$omp parallel do private(igrid,ngrid) schedule(dynamic)
     do igrid=1,ncache,nvector
        ngrid=MIN(nvector,ncache-igrid+1)
        call gradient_phi(ilevel,igrid,ngrid,icount)
     end do
#endif

     do idim=1,ndim
        call make_virtual_fine_dp(f(1,idim),ilevel)
     end do
     if(simple_boundary)call make_boundary_force(ilevel)

  endif

  !----------------------------------------------
  ! Compute gravity potential and maximum density
  !----------------------------------------------
  rho_loc =0.0; rho_all =0.0
  epot_loc=0.0; epot_all=0.0
  fourpi=4.0D0*ACOS(-1.0D0)
  if(cosmo)fourpi=1.5D0*omega_m*aexp
  fact=-dx_loc**ndim/fourpi/2.0D0

  ncache=active(ilevel)%ngrid
!$omp parallel do private(igrid,ngrid,i,ind,iskip,idim,ind_grid_w,ind_cell_w) &
!$omp& reduction(+:epot_loc) reduction(max:rho_loc)
  do igrid=1,ncache,nvector
     ngrid=MIN(nvector,ncache-igrid+1)
     do i=1,ngrid
        ind_grid_w(i)=active(ilevel)%igrid(igrid+i-1)
     end do
     do ind=1,twotondim
        iskip=ncoarse+(ind-1)*ngridmax
        do i=1,ngrid
           ind_cell_w(i)=iskip+ind_grid_w(i)
        end do
        do idim=1,ndim
           do i=1,ngrid
              if(son(ind_cell_w(i))==0)then
                 epot_loc=epot_loc+fact*f(ind_cell_w(i),idim)**2
              end if
           end do
        end do
        do i=1,ngrid
           rho_loc=MAX(rho_loc,dble(abs(rho(ind_cell_w(i)))))
        end do
     end do
  end do

#ifndef WITHOUTMPI
     call MPI_ALLREDUCE(epot_loc,epot_all,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
     call MPI_ALLREDUCE(rho_loc ,rho_all ,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,info)
     epot_loc=epot_all
     rho_loc =rho_all
#endif
     epot_tot=epot_tot+epot_loc
     rho_max(ilevel)=rho_loc

111 format('   Entering force_fine for level ',I2)
end subroutine force_fine

#ifdef HYDRO_CUDA
!#########################################################
!#########################################################
! Hybrid dispatch for gradient_phi
!#########################################################
!#########################################################
subroutine force_gradient_hybrid(ilevel, icount, ncache)
  use amr_commons
  use poisson_commons
  use force_hybrid_commons
  use cuda_commons
  use hydro_cuda_interface
  use iso_c_binding
  implicit none
  integer, intent(in) :: ilevel, icount, ncache

  integer :: igrid, ngrid, stream_slot
  type(force_gpu_state_t), pointer :: gs

  ! First-call: allocate GPU states
  if (force_hybrid_inited == 0) then
     allocate(force_gstates(0:7))
     force_hybrid_inited = 1
  end if

  !$omp parallel private(igrid, ngrid, stream_slot, gs)
  stream_slot = cuda_acquire_stream_c()
  if (stream_slot >= 0) then
     gs => force_gstates(stream_slot)
     call force_gstate_ensure(gs, FORCE_SUPER_SIZE)
     gs%off = 0
  end if

  !$omp do schedule(dynamic)
  do igrid = 1, ncache, nvector
     ngrid = MIN(nvector, ncache - igrid + 1)
     if (stream_slot >= 0) then
        call force_gpu_gather_batch(gs, ilevel, icount, igrid, ngrid, stream_slot)
     else
        call gradient_phi(ilevel, igrid, ngrid, icount)
     end if
  end do
  !$omp end do nowait

  if (stream_slot >= 0) then
     if (gs%off > 0) call force_gpu_flush_scatter(gs, stream_slot, ilevel)
     call cuda_release_stream_c(stream_slot)
  end if
  !$omp end parallel

end subroutine force_gradient_hybrid

!#########################################################
! GPU gather: fill phi stencil buffer for superbatch
!#########################################################
subroutine force_gpu_gather_batch(gs, ilevel, icount, igrid_start, ngrid, stream_slot)
  use amr_commons
  use poisson_commons
  use morton_hash
  use force_hybrid_commons
  implicit none

  type(force_gpu_state_t), intent(inout) :: gs
  integer, intent(in) :: ilevel, icount, igrid_start, ngrid, stream_slot

  integer :: i, ind, idim, nx_loc, iskip
  integer :: id1, id2, id3, id4
  integer :: ig1, ig2, ig3, ig4
  integer :: ih1, ih2, ih3, ih4
  real(dp) :: dx, scale, dx_loc
  integer, dimension(1:3,1:4,1:8) :: ggg, hhh

  integer, dimension(1:nvector) :: ind_grid
  integer, dimension(1:nvector, 0:twondim) :: igridn
  integer, dimension(1:nvector, 1:ndim) :: ind_left, ind_right
  real(dp), dimension(1:nvector, 1:twotondim, 1:ndim) :: phi_left, phi_right

  integer :: g, slot, off

  ! Load grid indices
  do i = 1, ngrid
     ind_grid(i) = active(ilevel)%igrid(igrid_start + i - 1)
  end do

  ! Mesh size
  dx = 0.5D0**ilevel
  nx_loc = icoarse_max - icoarse_min + 1
  scale = boxlen / dble(nx_loc)
  dx_loc = dx * scale

  ! Stencil lookup tables
  ggg(1,1,1:8)=(/1,0,1,0,1,0,1,0/); hhh(1,1,1:8)=(/2,1,4,3,6,5,8,7/)
  ggg(1,2,1:8)=(/0,2,0,2,0,2,0,2/); hhh(1,2,1:8)=(/2,1,4,3,6,5,8,7/)
  ggg(1,3,1:8)=(/1,1,1,1,1,1,1,1/); hhh(1,3,1:8)=(/1,2,3,4,5,6,7,8/)
  ggg(1,4,1:8)=(/2,2,2,2,2,2,2,2/); hhh(1,4,1:8)=(/1,2,3,4,5,6,7,8/)
  ggg(2,1,1:8)=(/3,3,0,0,3,3,0,0/); hhh(2,1,1:8)=(/3,4,1,2,7,8,5,6/)
  ggg(2,2,1:8)=(/0,0,4,4,0,0,4,4/); hhh(2,2,1:8)=(/3,4,1,2,7,8,5,6/)
  ggg(2,3,1:8)=(/3,3,3,3,3,3,3,3/); hhh(2,3,1:8)=(/1,2,3,4,5,6,7,8/)
  ggg(2,4,1:8)=(/4,4,4,4,4,4,4,4/); hhh(2,4,1:8)=(/1,2,3,4,5,6,7,8/)
  ggg(3,1,1:8)=(/5,5,5,5,0,0,0,0/); hhh(3,1,1:8)=(/5,6,7,8,1,2,3,4/)
  ggg(3,2,1:8)=(/0,0,0,0,6,6,6,6/); hhh(3,2,1:8)=(/5,6,7,8,1,2,3,4/)
  ggg(3,3,1:8)=(/5,5,5,5,5,5,5,5/); hhh(3,3,1:8)=(/1,2,3,4,5,6,7,8/)
  ggg(3,4,1:8)=(/6,6,6,6,6,6,6,6/); hhh(3,4,1:8)=(/1,2,3,4,5,6,7,8/)

  ! Gather neighboring grids
  do i = 1, ngrid
     igridn(i, 0) = ind_grid(i)
  end do
  do idim = 1, ndim
     do i = 1, ngrid
        igridn(i, 2*idim-1) = morton_nbor_grid(ind_grid(i), ilevel, 2*idim-1)
        igridn(i, 2*idim  ) = morton_nbor_grid(ind_grid(i), ilevel, 2*idim  )
        ind_left (i, idim) = morton_nbor_cell(ind_grid(i), ilevel, 2*idim-1)
        ind_right(i, idim) = morton_nbor_cell(ind_grid(i), ilevel, 2*idim  )
     end do
  end do

  ! Interpolate phi from upper level
  if (ilevel > levelmin) then
     do idim = 1, ndim
        call interpol_phi(ind_left(1,idim), phi_left(1,1,idim), ngrid, ilevel, icount)
        call interpol_phi(ind_right(1,idim), phi_right(1,1,idim), ngrid, ilevel, icount)
     end do
  end if

  ! Fill phi_buf and cell_buf
  off = gs%off
  do ind = 1, twotondim
     iskip = ncoarse + (ind - 1) * ngridmax
     do i = 1, ngrid
        gs%cell_buf(off + i, ind) = iskip + ind_grid(i)
     end do

     do idim = 1, ndim
        id1 = hhh(idim,1,ind); ig1 = ggg(idim,1,ind); ih1 = ncoarse+(id1-1)*ngridmax
        id2 = hhh(idim,2,ind); ig2 = ggg(idim,2,ind); ih2 = ncoarse+(id2-1)*ngridmax
        id3 = hhh(idim,3,ind); ig3 = ggg(idim,3,ind); ih3 = ncoarse+(id3-1)*ngridmax
        id4 = hhh(idim,4,ind); ig4 = ggg(idim,4,ind); ih4 = ncoarse+(id4-1)*ngridmax

        do i = 1, ngrid
           g = off + i
           slot = (g - 1) * 24 + (ind - 1) * 3 + idim

           ! phi1
           if (igridn(i, ig1) > 0) then
              gs%phi_buf(1, slot) = phi(igridn(i, ig1) + ih1)
           else
              gs%phi_buf(1, slot) = phi_left(i, id1, idim)
           end if
           ! phi2
           if (igridn(i, ig2) > 0) then
              gs%phi_buf(2, slot) = phi(igridn(i, ig2) + ih2)
           else
              gs%phi_buf(2, slot) = phi_right(i, id2, idim)
           end if
           ! phi3
           if (igridn(i, ig3) > 0) then
              gs%phi_buf(3, slot) = phi(igridn(i, ig3) + ih3)
           else
              gs%phi_buf(3, slot) = phi_left(i, id3, idim)
           end if
           ! phi4
           if (igridn(i, ig4) > 0) then
              gs%phi_buf(4, slot) = phi(igridn(i, ig4) + ih4)
           else
              gs%phi_buf(4, slot) = phi_right(i, id4, idim)
           end if
        end do
     end do
  end do

  gs%off = off + ngrid

  ! Flush if buffer is near full
  if (gs%off + nvector > gs%cap) then
     call force_gpu_flush_scatter(gs, stream_slot, ilevel)
  end if

end subroutine force_gpu_gather_batch

!#########################################################
! GPU flush: launch kernel, scatter results
!#########################################################
subroutine force_gpu_flush_scatter(gs, stream_slot, ilevel)
  use amr_commons
  use poisson_commons
  use force_hybrid_commons
  use hydro_cuda_interface
  use iso_c_binding
  implicit none

  type(force_gpu_state_t), intent(inout) :: gs
  integer, intent(in) :: stream_slot, ilevel

  real(dp) :: dx, a, b
  integer :: g, ind, idim, slot

  if (gs%off == 0) return

  dx = 0.5D0**ilevel
  a = 0.50D0 * 4.0D0 / 3.0D0 / dx
  b = 0.25D0 * 1.0D0 / 3.0D0 / dx

  ! Launch GPU kernel
  call gradient_phi_cuda_async_f(gs%phi_buf, gs%f_buf, a, b, &
       int(gs%off, c_int), int(stream_slot, c_int))
  call gradient_phi_cuda_sync_f(int(stream_slot, c_int))

  ! Scatter force values to global f() array
  do g = 1, gs%off
     do ind = 1, twotondim
        do idim = 1, ndim
           slot = (g - 1) * 24 + (ind - 1) * 3 + idim
           f(gs%cell_buf(g, ind), idim) = gs%f_buf(slot)
        end do
     end do
  end do

  gs%off = 0

end subroutine force_gpu_flush_scatter
#endif

!#########################################################
!#########################################################
! gradient_phi (unchanged from base version)
!#########################################################
!#########################################################
subroutine gradient_phi(ilevel,igrid,ngrid,icount)
  use amr_commons
  use pm_commons
  use hydro_commons
  use poisson_commons
  use morton_hash
  implicit none
  integer::ngrid,ilevel,icount
  integer,dimension(1:nvector)::ind_grid
  !-------------------------------------------------
  ! This routine compute the 3-force for all cells
  ! in grids ind_grid(:) at level ilevel, using a
  ! 5 nodes kernel (5 points FDA).
  !-------------------------------------------------
  integer::i,idim,ind,iskip,nx_loc,igrid
  integer::id1,id2,id3,id4
  integer::ig1,ig2,ig3,ig4
  integer::ih1,ih2,ih3,ih4
  real(dp)::dx,a,b,scale,dx_loc
  integer,dimension(1:3,1:4,1:8)::ggg,hhh
  integer ,dimension(1:nvector)::ind_cell
  integer ,dimension(1:nvector,1:ndim)::ind_left,ind_right
  integer ,dimension(1:nvector,0:twondim)::igridn
  real(dp),dimension(1:nvector)::phi1,phi2,phi3,phi4
  real(dp),dimension(1:nvector,1:twotondim,1:ndim)::phi_left,phi_right

  do i=1,ngrid
     ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
  end do

  ! Mesh size at level ilevel
  dx=0.5D0**ilevel

  ! Rescaling factor
  nx_loc=icoarse_max-icoarse_min+1
  scale=boxlen/dble(nx_loc)
  dx_loc=dx*scale

  a=0.50D0*4.0D0/3.0D0/dx
  b=0.25D0*1.0D0/3.0D0/dx
  ggg(1,1,1:8)=(/1,0,1,0,1,0,1,0/); hhh(1,1,1:8)=(/2,1,4,3,6,5,8,7/)
  ggg(1,2,1:8)=(/0,2,0,2,0,2,0,2/); hhh(1,2,1:8)=(/2,1,4,3,6,5,8,7/)
  ggg(1,3,1:8)=(/1,1,1,1,1,1,1,1/); hhh(1,3,1:8)=(/1,2,3,4,5,6,7,8/)
  ggg(1,4,1:8)=(/2,2,2,2,2,2,2,2/); hhh(1,4,1:8)=(/1,2,3,4,5,6,7,8/)
  ggg(2,1,1:8)=(/3,3,0,0,3,3,0,0/); hhh(2,1,1:8)=(/3,4,1,2,7,8,5,6/)
  ggg(2,2,1:8)=(/0,0,4,4,0,0,4,4/); hhh(2,2,1:8)=(/3,4,1,2,7,8,5,6/)
  ggg(2,3,1:8)=(/3,3,3,3,3,3,3,3/); hhh(2,3,1:8)=(/1,2,3,4,5,6,7,8/)
  ggg(2,4,1:8)=(/4,4,4,4,4,4,4,4/); hhh(2,4,1:8)=(/1,2,3,4,5,6,7,8/)
  ggg(3,1,1:8)=(/5,5,5,5,0,0,0,0/); hhh(3,1,1:8)=(/5,6,7,8,1,2,3,4/)
  ggg(3,2,1:8)=(/0,0,0,0,6,6,6,6/); hhh(3,2,1:8)=(/5,6,7,8,1,2,3,4/)
  ggg(3,3,1:8)=(/5,5,5,5,5,5,5,5/); hhh(3,3,1:8)=(/1,2,3,4,5,6,7,8/)
  ggg(3,4,1:8)=(/6,6,6,6,6,6,6,6/); hhh(3,4,1:8)=(/1,2,3,4,5,6,7,8/)

  ! Gather neighboring grids
  do i=1,ngrid
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

  ! Interpolate potential from upper level
  if (ilevel>levelmin)then
     do idim=1,ndim
        call interpol_phi(ind_left (1,idim),phi_left (1,1,idim),ngrid,ilevel,icount)
        call interpol_phi(ind_right(1,idim),phi_right(1,1,idim),ngrid,ilevel,icount)
     end do
  end if
  ! Loop over cells
  do ind=1,twotondim
     iskip=ncoarse+(ind-1)*ngridmax
     do i=1,ngrid
        ind_cell(i)=iskip+ind_grid(i)
     end do

     ! Loop over dimensions
     do idim=1,ndim

        ! Loop over nodes
        id1=hhh(idim,1,ind); ig1=ggg(idim,1,ind); ih1=ncoarse+(id1-1)*ngridmax
        id2=hhh(idim,2,ind); ig2=ggg(idim,2,ind); ih2=ncoarse+(id2-1)*ngridmax
        id3=hhh(idim,3,ind); ig3=ggg(idim,3,ind); ih3=ncoarse+(id3-1)*ngridmax
        id4=hhh(idim,4,ind); ig4=ggg(idim,4,ind); ih4=ncoarse+(id4-1)*ngridmax

        ! Gather potential
        do i=1,ngrid
           if(igridn(i,ig1)>0)then
              phi1(i)=phi(igridn(i,ig1)+ih1)
           else
              phi1(i)=phi_left(i,id1,idim)
           end if
        end do
        do i=1,ngrid
           if(igridn(i,ig2)>0)then
              phi2(i)=phi(igridn(i,ig2)+ih2)
           else
              phi2(i)=phi_right(i,id2,idim)
           end if
        end do
        do i=1,ngrid
           if(igridn(i,ig3)>0)then
              phi3(i)=phi(igridn(i,ig3)+ih3)
           else
              phi3(i)=phi_left(i,id3,idim)
           end if
        end do
        do i=1,ngrid
           if(igridn(i,ig4)>0)then
              phi4(i)=phi(igridn(i,ig4)+ih4)
           else
              phi4(i)=phi_right(i,id4,idim)
           end if
        end do
        do i=1,ngrid
           f(ind_cell(i),idim)=a*(phi1(i)-phi2(i)) &
                &             -b*(phi3(i)-phi4(i))
        end do
     end do
  end do

end subroutine gradient_phi
