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

     ! Apply MOND (QUMOND) correction to Newtonian force
     if(use_mond .and. mond_type == 0) then
        call apply_mond_force(ilevel)
     end if

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
!#########################################################
!#########################################################
subroutine apply_mond_force(ilevel)
  use amr_commons
  use poisson_commons
  implicit none
  integer,intent(in)::ilevel
  !-------------------------------------------------------
  ! QUMOND: multiply Newtonian force by nu(|g_N+g_ext|/a0)
  ! With EFE: f_d = nu * f_d + (nu-1) * g_ext_d
  !-------------------------------------------------------
  integer::igrid,ngrid,ncache,i,ind,iskip,idim
  real(dp)::scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2
  real(dp)::scale_a,a0_code,gnorm,x,nu
  real(dp)::g_ext_code(1:3)
  integer,dimension(1:nvector)::ind_cell

  ncache=active(ilevel)%ngrid
  if(ncache==0) return

  ! Convert a0 and g_ext from CGS to code units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
  scale_a = scale_l / scale_t**2
  a0_code = a0_mond / scale_a
  g_ext_code(1:3) = g_ext_mond(1:3) / scale_a

!$omp parallel do private(igrid,ngrid,i,ind,iskip,idim,ind_cell, &
!$omp&  gnorm,x,nu) schedule(dynamic)
  do igrid=1,ncache,nvector
     ngrid=MIN(nvector,ncache-igrid+1)

     do ind=1,twotondim
        iskip=ncoarse+(ind-1)*ngridmax
        do i=1,ngrid
           ind_cell(i)=iskip+active(ilevel)%igrid(igrid+i-1)
        end do

        do i=1,ngrid
           gnorm = sqrt((f(ind_cell(i),1)+g_ext_code(1))**2 &
                      + (f(ind_cell(i),2)+g_ext_code(2))**2 &
                      + (f(ind_cell(i),3)+g_ext_code(3))**2)

           if(gnorm > 1d-30*a0_code) then
              x = gnorm / a0_code
              ! Compute nu-function
              if(mond_mu_type == 1) then
                 nu = 0.5d0 + 0.5d0*sqrt(1d0 + 4d0/x)
              else
                 nu = sqrt(0.5d0 + 0.5d0*sqrt(1d0 + 4d0/(x*x)))
              end if
              ! f_d = nu * f_d + (nu-1) * g_ext_d
              do idim=1,ndim
                 f(ind_cell(i),idim) = nu * f(ind_cell(i),idim) &
                      + (nu - 1d0) * g_ext_code(idim)
              end do
           end if
        end do
     end do
  end do

end subroutine apply_mond_force
!#########################################################
!#########################################################
subroutine compute_mond_phantom_density(ilevel, is_aqual)
  use amr_commons
  use poisson_commons
  use morton_hash
  implicit none
  integer,intent(in)::ilevel
  logical,intent(in)::is_aqual
  !-------------------------------------------------------
  ! Compute phantom density and add to rho() in-place.
  !
  ! QUMOND (is_aqual=.false.):
  !   rho_ph = -div[(nu-1)*(f+g_ext)] / fourpi  (nu-1 > 0)
  ! AQUAL (is_aqual=.true.):
  !   rho_ph = +div[(mu-1)*(f+g_ext)] / fourpi  (mu-1 < 0)
  !
  ! Uses f() (force, already boundary-exchanged)
  ! with 6-neighbor centered difference for divergence.
  ! Cells at coarse-fine boundaries (igridn=0) are skipped.
  !-------------------------------------------------------
  integer::igrid,ngrid,ncache,i,ind,iskip,idim
  real(dp)::scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2
  real(dp)::scale_a,a0_code
  real(dp)::g_ext_code(1:3)
  real(dp)::dx,scale,dx_loc,fourpi
  integer::nx_loc
  integer::ig_left,ig_right,ih_left,ih_right
  integer::ind_nb_left,ind_nb_right
  real(dp)::gnorm_nb,x_nb,func_m1
  real(dp)::h_right,h_left,div_h
  logical::skip_cell

  integer,dimension(1:3,1:2,1:8)::ggg,hhh
  integer,dimension(1:nvector)::ind_grid_w,ind_cell_w
  integer,dimension(1:nvector,0:twondim)::igridn_w

  ncache=active(ilevel)%ngrid
  if(ncache==0) return

  ! Unit conversion
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
  scale_a = scale_l / scale_t**2
  a0_code = a0_mond / scale_a
  g_ext_code(1:3) = g_ext_mond(1:3) / scale_a

  ! Mesh size
  dx=0.5D0**ilevel
  nx_loc=(icoarse_max-icoarse_min+1)
  scale=boxlen/dble(nx_loc)
  dx_loc=dx*scale

  ! 4piG factor (same convention as force_fine)
  fourpi=4.0D0*ACOS(-1.0D0)
  if(cosmo)fourpi=1.5D0*omega_m*aexp

  ! Neighbor lookup tables (left=1, right=2 per dimension)
  ggg(1,1,1:8)=(/1,0,1,0,1,0,1,0/); hhh(1,1,1:8)=(/2,1,4,3,6,5,8,7/)
  ggg(1,2,1:8)=(/0,2,0,2,0,2,0,2/); hhh(1,2,1:8)=(/2,1,4,3,6,5,8,7/)
  ggg(2,1,1:8)=(/3,3,0,0,3,3,0,0/); hhh(2,1,1:8)=(/3,4,1,2,7,8,5,6/)
  ggg(2,2,1:8)=(/0,0,4,4,0,0,4,4/); hhh(2,2,1:8)=(/3,4,1,2,7,8,5,6/)
  ggg(3,1,1:8)=(/5,5,5,5,0,0,0,0/); hhh(3,1,1:8)=(/5,6,7,8,1,2,3,4/)
  ggg(3,2,1:8)=(/0,0,0,0,6,6,6,6/); hhh(3,2,1:8)=(/5,6,7,8,1,2,3,4/)

!$omp parallel do private(igrid,ngrid,i,ind,iskip,idim, &
!$omp&  ind_grid_w,ind_cell_w,igridn_w, &
!$omp&  ig_left,ig_right,ih_left,ih_right, &
!$omp&  ind_nb_left,ind_nb_right, &
!$omp&  gnorm_nb,x_nb,func_m1, &
!$omp&  h_right,h_left,div_h,skip_cell) schedule(dynamic)
  do igrid=1,ncache,nvector
     ngrid=MIN(nvector,ncache-igrid+1)
     do i=1,ngrid
        ind_grid_w(i)=active(ilevel)%igrid(igrid+i-1)
     end do

     ! Gather neighboring grids (6 faces)
     do i=1,ngrid
        igridn_w(i,0)=ind_grid_w(i)
     end do
     do idim=1,ndim
        do i=1,ngrid
           igridn_w(i,2*idim-1)=morton_nbor_grid(ind_grid_w(i),ilevel,2*idim-1)
           igridn_w(i,2*idim  )=morton_nbor_grid(ind_grid_w(i),ilevel,2*idim  )
        end do
     end do

     do ind=1,twotondim
        iskip=ncoarse+(ind-1)*ngridmax
        do i=1,ngrid
           ind_cell_w(i)=iskip+ind_grid_w(i)
        end do

        do i=1,ngrid
           ! Check if all neighbors exist (skip at coarse-fine boundaries)
           skip_cell=.false.
           do idim=1,ndim
              ig_left =ggg(idim,1,ind)
              ig_right=ggg(idim,2,ind)
              if(igridn_w(i,ig_left)==0 .or. igridn_w(i,ig_right)==0) then
                 skip_cell=.true.
                 exit
              end if
           end do
           if(skip_cell) cycle

           ! Compute divergence of (func-1)*(f+g_ext)
           div_h = 0d0
           do idim=1,ndim
              ig_left =ggg(idim,1,ind)
              ig_right=ggg(idim,2,ind)
              ih_left =ncoarse+(hhh(idim,1,ind)-1)*ngridmax
              ih_right=ncoarse+(hhh(idim,2,ind)-1)*ngridmax

              ind_nb_left =igridn_w(i,ig_left )+ih_left
              ind_nb_right=igridn_w(i,ig_right)+ih_right

              ! Left neighbor
              gnorm_nb = sqrt((f(ind_nb_left,1)+g_ext_code(1))**2 &
                            + (f(ind_nb_left,2)+g_ext_code(2))**2 &
                            + (f(ind_nb_left,3)+g_ext_code(3))**2)
              x_nb = gnorm_nb / a0_code
              if(is_aqual) then
                 call get_mond_mu_minus1(x_nb, mond_mu_type, func_m1)
              else
                 call get_mond_nu_minus1(x_nb, mond_mu_type, func_m1)
              end if
              h_left = func_m1 * (f(ind_nb_left, idim) + g_ext_code(idim))

              ! Right neighbor
              gnorm_nb = sqrt((f(ind_nb_right,1)+g_ext_code(1))**2 &
                            + (f(ind_nb_right,2)+g_ext_code(2))**2 &
                            + (f(ind_nb_right,3)+g_ext_code(3))**2)
              x_nb = gnorm_nb / a0_code
              if(is_aqual) then
                 call get_mond_mu_minus1(x_nb, mond_mu_type, func_m1)
              else
                 call get_mond_nu_minus1(x_nb, mond_mu_type, func_m1)
              end if
              h_right = func_m1 * (f(ind_nb_right, idim) + g_ext_code(idim))

              div_h = div_h + (h_right - h_left)
           end do

           ! QUMOND: rho_ph = -div[(nu-1)*(f+g_ext)] / (2*dx) / fourpi
           ! AQUAL:  rho_ph = +div[(mu-1)*(f+g_ext)] / (2*dx) / fourpi
           if(is_aqual) then
              rho(ind_cell_w(i)) = rho(ind_cell_w(i)) + div_h / (2d0*dx_loc) / fourpi
           else
              rho(ind_cell_w(i)) = rho(ind_cell_w(i)) - div_h / (2d0*dx_loc) / fourpi
           end if

        end do  ! i
     end do  ! ind
  end do  ! igrid

end subroutine compute_mond_phantom_density
!#########################################################
!#########################################################
subroutine get_mond_nu_minus1(x, mu_type, num1)
  use amr_parameters, only: dp
  implicit none
  real(dp),intent(in)::x
  integer,intent(in)::mu_type
  real(dp),intent(out)::num1
  ! Returns nu(x) - 1 for MOND interpolation function
  if(x < 1d-30) then
     num1 = 0d0
     return
  end if
  if(mu_type == 1) then
     ! Simple: mu=x/(1+x) -> nu = 0.5*(1+sqrt(1+4/x))
     num1 = 0.5d0*(-1d0 + sqrt(1d0 + 4d0/x))
  else
     ! Standard: mu=x/sqrt(1+x^2) -> nu = sqrt(0.5+0.5*sqrt(1+4/x^2))
     num1 = sqrt(0.5d0 + 0.5d0*sqrt(1d0 + 4d0/(x*x))) - 1d0
  end if
end subroutine get_mond_nu_minus1
!#########################################################
!#########################################################
subroutine get_mond_mu_minus1(x, mu_type, mum1)
  use amr_parameters, only: dp
  implicit none
  real(dp),intent(in)::x
  integer,intent(in)::mu_type
  real(dp),intent(out)::mum1
  ! Returns mu(x) - 1 for MOND interpolation function
  ! mu(x) - 1 < 0 always (mu < 1 in deep-MOND)
  if(x < 1d-30) then
     mum1 = -1d0   ! mu(0) = 0
     return
  end if
  if(mu_type == 1) then
     ! Simple: mu=x/(1+x) -> mu-1 = -1/(1+x)
     mum1 = -1d0 / (1d0 + x)
  else
     ! Standard: mu=x/sqrt(1+x^2) -> mu-1 = x/sqrt(1+x^2) - 1
     mum1 = x / sqrt(1d0 + x*x) - 1d0
  end if
end subroutine get_mond_mu_minus1
!#########################################################
!#########################################################
subroutine aqual_iterate(ilevel, icount)
  use amr_commons
  use poisson_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer, intent(in) :: ilevel, icount
  !-------------------------------------------------------
  ! AQUAL Phase 2: iterative fixed-point solver
  ! Solves: div[mu(|grad(Phi)|/a0) * grad(Phi)] = 4piG*rho
  ! by iterating:
  !   1. Compute AQUAL phantom density from current forces
  !   2. Re-solve Poisson with rho_bary + rho_phantom
  !   3. Re-compute forces
  !   4. Check convergence: max|f_new - f_old| / max|f_new|
  !-------------------------------------------------------
  integer :: iter, ncache, ind, iskip, i, idx, idim, info
  real(dp) :: delta_f_local, f_norm_local
  real(dp) :: delta_f_global, f_norm_global, rel_change
  logical :: converged

  real(dp), allocatable :: rho_bary(:)
  real(dp), allocatable :: f_old(:,:)
  integer :: ncell_active
  integer, allocatable :: ind_cell_list(:)

  ncache = active(ilevel)%ngrid
  if(ncache == 0) return

  ncell_active = ncache * twotondim
  allocate(rho_bary(1:ncell_active))
  allocate(f_old(1:ncell_active, 1:ndim))
  allocate(ind_cell_list(1:ncell_active))

  ! Build cell list and save rho_bary
  idx = 0
  do ind = 1, twotondim
     iskip = ncoarse + (ind-1) * ngridmax
     do i = 1, ncache
        idx = idx + 1
        ind_cell_list(idx) = iskip + active(ilevel)%igrid(i)
        rho_bary(idx) = rho(ind_cell_list(idx))
     end do
  end do

  ! Save current forces (Newtonian) as f_old
  do idx = 1, ncell_active
     do idim = 1, ndim
        f_old(idx, idim) = f(ind_cell_list(idx), idim)
     end do
  end do

  converged = .false.
  do iter = 1, n_iter_mond

     ! (a) Restore rho to baryonic
     do idx = 1, ncell_active
        rho(ind_cell_list(idx)) = rho_bary(idx)
     end do

     ! (b) Compute AQUAL phantom density (mu-1, current forces)
     call compute_mond_phantom_density(ilevel, .true.)
     call make_virtual_fine_dp(rho(1), ilevel)

     ! (c) Re-solve Poisson
     if(ilevel > levelmin) then
        if(ilevel .ge. cg_levelmin) then
           call phi_fine_cg(ilevel, icount)
        else
           call multigrid_fine(ilevel, icount)
        end if
     else
        call multigrid_fine(levelmin, icount)
     end if

     ! (d) Re-compute forces
     call force_fine(ilevel, icount)

     ! (e) Convergence check
     delta_f_local = 0d0
     f_norm_local = 0d0
!$omp parallel do private(idx, idim) &
!$omp& reduction(max:delta_f_local, f_norm_local)
     do idx = 1, ncell_active
        do idim = 1, ndim
           delta_f_local = max(delta_f_local, &
                abs(f(ind_cell_list(idx), idim) - f_old(idx, idim)))
           f_norm_local = max(f_norm_local, &
                abs(f(ind_cell_list(idx), idim)))
        end do
     end do

     ! Save forces for next convergence check
     do idx = 1, ncell_active
        do idim = 1, ndim
           f_old(idx, idim) = f(ind_cell_list(idx), idim)
        end do
     end do

#ifndef WITHOUTMPI
     call MPI_ALLREDUCE(delta_f_local, delta_f_global, 1, &
          MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, info)
     call MPI_ALLREDUCE(f_norm_local, f_norm_global, 1, &
          MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, info)
#else
     delta_f_global = delta_f_local
     f_norm_global = f_norm_local
#endif

     if(f_norm_global > 0d0) then
        rel_change = delta_f_global / f_norm_global
     else
        rel_change = 0d0
     end if

     if(myid == 1) then
        write(*,'(A,I3,A,ES12.4)') &
             ' AQUAL iter ', iter, ': delta_f/f_norm = ', rel_change
     end if

     if(rel_change < mond_eps) then
        converged = .true.
        if(myid == 1) write(*,'(A,I3,A)') &
             ' AQUAL converged in ', iter, ' iterations'
        exit
     end if
  end do

  if(.not. converged .and. myid == 1) then
     write(*,'(A,I3,A,ES12.4)') &
          ' WARNING: AQUAL did NOT converge after ', n_iter_mond, &
          ' iters, rel_change=', rel_change
  end if

  deallocate(rho_bary, f_old, ind_cell_list)
end subroutine aqual_iterate
!#########################################################
!#########################################################
! f(R) Hu-Sawicki gravity + nDGP gravity solvers
!#########################################################
!#########################################################

!=========================================================
! compute_fifth_force: add gradient of scalar_gr to f()
! factor = -0.5 for f(R), -1/(2*beta) for nDGP
!=========================================================
subroutine compute_fifth_force(ilevel, factor)
  use amr_commons
  use poisson_commons
  use morton_hash
  implicit none
  integer,intent(in)::ilevel
  real(dp),intent(in)::factor
  !-------------------------------------------------------
  ! Compute F5_d = factor * (scalar_gr(right)-scalar_gr(left))/(2*dx)
  ! and add to f(icell, idim).
  ! Uses simple 2-point centered difference (same as rho divergence).
  !-------------------------------------------------------
  integer::igrid,ngrid,ncache,i,ind,iskip,idim
  integer::ig_left,ig_right,ih_left,ih_right
  integer::ind_nb_left,ind_nb_right
  real(dp)::dx,scale,dx_loc
  integer::nx_loc
  real(dp)::grad_u

  integer,dimension(1:3,1:2,1:8)::ggg,hhh
  integer,dimension(1:nvector)::ind_grid_w,ind_cell_w
  integer,dimension(1:nvector,0:twondim)::igridn_w

  ncache=active(ilevel)%ngrid
  if(ncache==0) return

  dx=0.5D0**ilevel
  nx_loc=(icoarse_max-icoarse_min+1)
  scale=boxlen/dble(nx_loc)
  dx_loc=dx*scale

  ! Neighbor lookup tables (left=1, right=2)
  ggg(1,1,1:8)=(/1,0,1,0,1,0,1,0/); hhh(1,1,1:8)=(/2,1,4,3,6,5,8,7/)
  ggg(1,2,1:8)=(/0,2,0,2,0,2,0,2/); hhh(1,2,1:8)=(/2,1,4,3,6,5,8,7/)
  ggg(2,1,1:8)=(/3,3,0,0,3,3,0,0/); hhh(2,1,1:8)=(/3,4,1,2,7,8,5,6/)
  ggg(2,2,1:8)=(/0,0,4,4,0,0,4,4/); hhh(2,2,1:8)=(/3,4,1,2,7,8,5,6/)
  ggg(3,1,1:8)=(/5,5,5,5,0,0,0,0/); hhh(3,1,1:8)=(/5,6,7,8,1,2,3,4/)
  ggg(3,2,1:8)=(/0,0,0,0,6,6,6,6/); hhh(3,2,1:8)=(/5,6,7,8,1,2,3,4/)

!$omp parallel do private(igrid,ngrid,i,ind,iskip,idim, &
!$omp&  ind_grid_w,ind_cell_w,igridn_w, &
!$omp&  ig_left,ig_right,ih_left,ih_right, &
!$omp&  ind_nb_left,ind_nb_right,grad_u) schedule(dynamic)
  do igrid=1,ncache,nvector
     ngrid=MIN(nvector,ncache-igrid+1)
     do i=1,ngrid
        ind_grid_w(i)=active(ilevel)%igrid(igrid+i-1)
     end do

     ! Gather neighboring grids
     do i=1,ngrid
        igridn_w(i,0)=ind_grid_w(i)
     end do
     do idim=1,ndim
        do i=1,ngrid
           igridn_w(i,2*idim-1)=morton_nbor_grid(ind_grid_w(i),ilevel,2*idim-1)
           igridn_w(i,2*idim  )=morton_nbor_grid(ind_grid_w(i),ilevel,2*idim  )
        end do
     end do

     do ind=1,twotondim
        iskip=ncoarse+(ind-1)*ngridmax
        do i=1,ngrid
           ind_cell_w(i)=iskip+ind_grid_w(i)
        end do

        do i=1,ngrid
           do idim=1,ndim
              ig_left =ggg(idim,1,ind)
              ig_right=ggg(idim,2,ind)

              ! Skip at coarse-fine boundaries
              if(igridn_w(i,ig_left)==0 .or. igridn_w(i,ig_right)==0) cycle

              ih_left =ncoarse+(hhh(idim,1,ind)-1)*ngridmax
              ih_right=ncoarse+(hhh(idim,2,ind)-1)*ngridmax

              ind_nb_left =igridn_w(i,ig_left )+ih_left
              ind_nb_right=igridn_w(i,ig_right)+ih_right

              grad_u = (scalar_gr(ind_nb_right) - scalar_gr(ind_nb_left)) / (2d0*dx_loc)
              f(ind_cell_w(i),idim) = f(ind_cell_w(i),idim) + factor * grad_u
           end do
        end do
     end do
  end do

end subroutine compute_fifth_force

!=========================================================
! fR_background: compute background Ricci scalar and f_R
!=========================================================
subroutine fR_background(aa, R_bar, fR_bar)
  use amr_parameters, only: dp, omega_m, omega_l, fR0, fR_n
  implicit none
  real(dp),intent(in) ::aa
  real(dp),intent(out)::R_bar, fR_bar
  !-------------------------------------------------------
  ! R_bar = 3*H0^2*(Omega_m/a^3 + 4*Omega_Lambda)
  ! In code units where H0=1:
  !   R_bar = 3*(omega_m/a^3 + 4*omega_l)
  ! f_R_bar = -|fR0| * n * (R_bar0/R_bar)^(n+1)
  !   where R_bar0 = R_bar(a=1)
  !-------------------------------------------------------
  real(dp)::R_bar0
  integer::np1

  np1 = fR_n + 1
  R_bar  = 3d0 * (omega_m / aa**3 + 4d0 * omega_l)
  R_bar0 = 3d0 * (omega_m + 4d0 * omega_l)
  fR_bar = fR0 * dble(fR_n) * (R_bar0 / R_bar)**np1

end subroutine fR_background

!=========================================================
! fR_solve_level: top-level f(R) solver for one AMR level
!=========================================================
subroutine fR_solve_level(ilevel, icount)
  use amr_commons
  use poisson_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer,intent(in)::ilevel, icount
  !-------------------------------------------------------
  ! 1. Initialize scalar_gr to background fR_bar
  ! 2. Newton-GS iteration
  ! 3. Add fifth force F5 = -0.5 * grad(fR)
  !-------------------------------------------------------
  integer::iter,ncache,info
  real(dp)::R_bar,fR_bar
  real(dp)::res_max_local,res_max_global
  real(dp)::src_max_local,src_max_global
  real(dp)::rel_res
  logical::converged

  ncache=active(ilevel)%ngrid
  if(ncache==0) return

  ! Get background values
  call fR_background(aexp, R_bar, fR_bar)

  ! Initialize scalar_gr to background value on first call
  ! (scalar_gr_old provides a warm start from previous step)
  if(nstep==0) then
     call fR_init_scalar(ilevel, fR_bar)
  end if

  ! Newton-GS relaxation
  converged = .false.
  do iter=1,n_iter_fR

     call fR_gauss_seidel(ilevel, R_bar, fR_bar, res_max_local, src_max_local)

     ! Exchange boundaries after each sweep
     call make_virtual_fine_dp(scalar_gr(1), ilevel)

     ! Global convergence check
#ifndef WITHOUTMPI
     call MPI_ALLREDUCE(res_max_local, res_max_global, 1, &
          MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, info)
     call MPI_ALLREDUCE(src_max_local, src_max_global, 1, &
          MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, info)
#else
     res_max_global = res_max_local
     src_max_global = src_max_local
#endif

     if(src_max_global > 0d0) then
        rel_res = res_max_global / src_max_global
     else
        rel_res = 0d0
     end if

     if(rel_res < fR_eps) then
        converged = .true.
        if(myid==1) write(*,'(A,I2,A,I3,A,ES10.3)') &
             ' f(R) level ',ilevel,' converged in ',iter,' iters, res=',rel_res
        exit
     end if
  end do

  if(.not. converged .and. myid==1) then
     write(*,'(A,I2,A,I3,A,ES10.3)') &
          ' WARNING: f(R) level ',ilevel,' NOT converged after ', &
          n_iter_fR,' iters, res=',rel_res
  end if

  ! Save scalar_gr for warm start next step
  call fR_save_old(ilevel)

  ! Add fifth force: F5 = -(1/2) * grad(fR)
  call compute_fifth_force(ilevel, -0.5d0)

  ! Exchange f boundaries (force_fine will do this later, but be safe)
  ! (No — force boundaries are exchanged by the caller in amr_step)

end subroutine fR_solve_level

!=========================================================
! fR_init_scalar: initialize scalar_gr at a level
!=========================================================
subroutine fR_init_scalar(ilevel, fR_bar)
  use amr_commons
  use poisson_commons
  implicit none
  integer,intent(in)::ilevel
  real(dp),intent(in)::fR_bar

  integer::igrid,i,ind,iskip,ncache,icell

  ncache=active(ilevel)%ngrid
!$omp parallel do private(igrid,i,ind,iskip,icell) schedule(dynamic)
  do igrid=1,ncache,nvector
     do i=1,MIN(nvector,ncache-igrid+1)
        do ind=1,twotondim
           iskip=ncoarse+(ind-1)*ngridmax
           icell=iskip+active(ilevel)%igrid(igrid+i-1)
           scalar_gr(icell) = fR_bar
           scalar_gr_old(icell) = fR_bar
        end do
     end do
  end do

end subroutine fR_init_scalar

!=========================================================
! fR_save_old: save scalar_gr → scalar_gr_old
!=========================================================
subroutine fR_save_old(ilevel)
  use amr_commons
  use poisson_commons
  implicit none
  integer,intent(in)::ilevel

  integer::igrid,i,ind,iskip,ncache,icell

  ncache=active(ilevel)%ngrid
!$omp parallel do private(igrid,i,ind,iskip,icell) schedule(dynamic)
  do igrid=1,ncache,nvector
     do i=1,MIN(nvector,ncache-igrid+1)
        do ind=1,twotondim
           iskip=ncoarse+(ind-1)*ngridmax
           icell=iskip+active(ilevel)%igrid(igrid+i-1)
           scalar_gr_old(icell) = scalar_gr(icell)
        end do
     end do
  end do

end subroutine fR_save_old

!=========================================================
! fR_gauss_seidel: one Newton-GS sweep for f(R) equation
! ∇²f_R = -(a²/3)[R(f_R) - R_bar] + (a²/3)*8πG*δρ
!
! In RAMSES code units (H0=1, fourpi=1.5*Ωm*a):
!   ∇²f_R = -(a²/3)*[R(f_R) - R_bar] + 2*a*Ωm*δρ/ρ_bar
!
! Hu-Sawicki inversion: R(f_R) from f_R ↔ R relation
!=========================================================
subroutine fR_gauss_seidel(ilevel, R_bar, fR_bar, res_max, src_max)
  use amr_commons
  use poisson_commons
  use morton_hash
  implicit none
  integer,intent(in)::ilevel
  real(dp),intent(in)::R_bar, fR_bar
  real(dp),intent(out)::res_max, src_max
  !-------------------------------------------------------
  ! Newton-GS update: u_new = u - F/J
  !   F = laplacian(u) - source(u)
  !   J = -6/dx² - dS/du
  ! Red-black ordering
  !-------------------------------------------------------
  integer::igrid,ngrid,ncache,i,ind,iskip,idim
  integer::ig_left,ig_right,ih_left,ih_right
  integer::ind_nb_left,ind_nb_right
  real(dp)::dx,scale,dx_loc,dx2,dx2_inv
  integer::nx_loc
  real(dp)::u_c,lapl,source,R_of_u,dR_du,residual,jacobian,delta_u
  real(dp)::a2_over_3,rho_bar
  real(dp)::R_bar0
  integer::np1,icolor

  integer,dimension(1:3,1:2,1:8)::ggg,hhh
  integer,dimension(1:nvector)::ind_grid_w,ind_cell_w
  integer,dimension(1:nvector,0:twondim)::igridn_w

  ncache=active(ilevel)%ngrid
  if(ncache==0) then
     res_max=0d0; src_max=0d0
     return
  end if

  dx=0.5D0**ilevel
  nx_loc=(icoarse_max-icoarse_min+1)
  scale=boxlen/dble(nx_loc)
  dx_loc=dx*scale
  dx2=dx_loc**2
  dx2_inv=1d0/dx2

  np1 = fR_n + 1
  a2_over_3 = aexp**2 / 3d0
  R_bar0 = 3d0 * (omega_m + 4d0 * omega_l)
  ! Mean density in code units: ρ_bar = Ω_m / a^3 (with H0=1)
  rho_bar = omega_m / aexp**3

  ! Neighbor lookup
  ggg(1,1,1:8)=(/1,0,1,0,1,0,1,0/); hhh(1,1,1:8)=(/2,1,4,3,6,5,8,7/)
  ggg(1,2,1:8)=(/0,2,0,2,0,2,0,2/); hhh(1,2,1:8)=(/2,1,4,3,6,5,8,7/)
  ggg(2,1,1:8)=(/3,3,0,0,3,3,0,0/); hhh(2,1,1:8)=(/3,4,1,2,7,8,5,6/)
  ggg(2,2,1:8)=(/0,0,4,4,0,0,4,4/); hhh(2,2,1:8)=(/3,4,1,2,7,8,5,6/)
  ggg(3,1,1:8)=(/5,5,5,5,0,0,0,0/); hhh(3,1,1:8)=(/5,6,7,8,1,2,3,4/)
  ggg(3,2,1:8)=(/0,0,0,0,6,6,6,6/); hhh(3,2,1:8)=(/5,6,7,8,1,2,3,4/)

  res_max = 0d0
  src_max = 0d0

  ! Red-black sweep (icolor=0: red, icolor=1: black)
  do icolor=0,1

!$omp parallel do private(igrid,ngrid,i,ind,iskip,idim, &
!$omp&  ind_grid_w,ind_cell_w,igridn_w, &
!$omp&  ig_left,ig_right,ih_left,ih_right, &
!$omp&  ind_nb_left,ind_nb_right, &
!$omp&  u_c,lapl,source,R_of_u,dR_du,residual,jacobian,delta_u) &
!$omp& reduction(max:res_max,src_max) schedule(dynamic)
     do igrid=1,ncache,nvector
        ngrid=MIN(nvector,ncache-igrid+1)
        do i=1,ngrid
           ind_grid_w(i)=active(ilevel)%igrid(igrid+i-1)
        end do

        ! Gather neighbors
        do i=1,ngrid
           igridn_w(i,0)=ind_grid_w(i)
        end do
        do idim=1,ndim
           do i=1,ngrid
              igridn_w(i,2*idim-1)=morton_nbor_grid(ind_grid_w(i),ilevel,2*idim-1)
              igridn_w(i,2*idim  )=morton_nbor_grid(ind_grid_w(i),ilevel,2*idim  )
           end do
        end do

        do ind=1,twotondim
           ! Red-black coloring
           if(mod(ind-1, 2) /= icolor) cycle

           iskip=ncoarse+(ind-1)*ngridmax
           do i=1,ngrid
              ind_cell_w(i)=iskip+ind_grid_w(i)
           end do

           do i=1,ngrid
              ! Check all neighbors exist
              ig_left =ggg(1,1,ind)
              ig_right=ggg(1,2,ind)
              if(igridn_w(i,ig_left)==0 .or. igridn_w(i,ig_right)==0) cycle
              ig_left =ggg(2,1,ind)
              ig_right=ggg(2,2,ind)
              if(igridn_w(i,ig_left)==0 .or. igridn_w(i,ig_right)==0) cycle
              ig_left =ggg(3,1,ind)
              ig_right=ggg(3,2,ind)
              if(igridn_w(i,ig_left)==0 .or. igridn_w(i,ig_right)==0) cycle

              u_c = scalar_gr(ind_cell_w(i))

              ! Compute Laplacian
              lapl = 0d0
              do idim=1,ndim
                 ig_left =ggg(idim,1,ind)
                 ig_right=ggg(idim,2,ind)
                 ih_left =ncoarse+(hhh(idim,1,ind)-1)*ngridmax
                 ih_right=ncoarse+(hhh(idim,2,ind)-1)*ngridmax
                 ind_nb_left =igridn_w(i,ig_left )+ih_left
                 ind_nb_right=igridn_w(i,ig_right)+ih_right
                 lapl = lapl + (scalar_gr(ind_nb_left) + scalar_gr(ind_nb_right) - 2d0*u_c) * dx2_inv
              end do

              ! R(f_R) from Hu-Sawicki inversion:
              ! f_R = fR0 * n * (R_bar0/R)^(n+1)
              ! → R = R_bar0 * (n*fR0/f_R)^(1/(n+1))
              ! Since fR0<0 and f_R<0, use absolute values
              if(abs(u_c) > 1d-30*abs(fR0)) then
                 R_of_u = R_bar0 * (dble(fR_n)*abs(fR0)/abs(u_c))**(1d0/dble(np1))
              else
                 R_of_u = R_bar  ! fallback to background
              end if

              ! Source = -(a²/3)*[R(fR) - R_bar] + (a²/3)*8πG*δρ
              ! In code units: δρ/ρ_bar = rho(cell)/rho_bar - 1
              ! and 8πG*ρ_bar = 3*Ω_m*H0²/a³ = 3*omega_m/a³  (H0=1)
              ! So (a²/3)*8πG*δρ = a²*omega_m/a³ * (rho(cell)/rho_bar - 1)
              !                  = omega_m/(a) * (ρ/ρ_bar - 1)
              ! But rho() in RAMSES code is overdensity: rho_code = ρ/ρ_bar - 1 when fourpi is set properly
              ! Actually, in super-comoving coords, the Poisson source is:
              !   ∇²Φ = 1.5*Ωm*a * (ρ/ρ_bar - 1) = fourpi * δ
              ! So for f(R): source = -a²/3*(R(fR) - R_bar) + 2*fourpi/3 * δ
              !   where fourpi = 1.5*Ωm*a and δ = rho(cell) (already overdensity)
              source = -a2_over_3*(R_of_u - R_bar) + a2_over_3 * 2d0 * 1.5d0 * omega_m * aexp * rho(ind_cell_w(i))

              ! Residual: F = laplacian - source
              residual = lapl - source

              ! Jacobian of source w.r.t. u_c:
              ! dS/du = -(a²/3) * dR/dfR
              ! dR/dfR = -R/(n+1)/fR  (from the inversion formula)
              if(abs(u_c) > 1d-30*abs(fR0)) then
                 dR_du = -R_of_u / (dble(np1) * u_c)
              else
                 dR_du = 0d0
              end if

              ! Jacobian: J = -6/dx² - dSource/du = -6/dx² + (a²/3)*dR/du
              jacobian = -6d0*dx2_inv + a2_over_3*dR_du

              ! Newton update
              if(abs(jacobian) > 1d-30) then
                 delta_u = -residual / jacobian
                 ! Clamp update to prevent overshoot
                 if(abs(delta_u) > 0.5d0*abs(u_c) .and. abs(u_c) > 1d-30*abs(fR0)) then
                    delta_u = sign(0.5d0*abs(u_c), delta_u)
                 end if
                 scalar_gr(ind_cell_w(i)) = u_c + delta_u
                 ! Enforce f_R < 0 (physical constraint)
                 if(scalar_gr(ind_cell_w(i)) > 0d0) scalar_gr(ind_cell_w(i)) = 0.5d0*u_c
              end if

              ! Track convergence
              res_max = max(res_max, abs(residual))
              src_max = max(src_max, abs(source))

           end do  ! i
        end do  ! ind
     end do  ! igrid

  end do  ! icolor

end subroutine fR_gauss_seidel

!=========================================================
! nDGP_beta: compute β(a) parameter
!=========================================================
function nDGP_beta(aa, orc, branch) result(beta)
  use amr_parameters, only: dp, omega_m, omega_l
  implicit none
  real(dp),intent(in)::aa, orc
  integer,intent(in)::branch
  real(dp)::beta
  !-------------------------------------------------------
  ! H(a) = H0 * sqrt(Ω_m/a³ + Ω_Λ)  [flat ΛCDM]
  ! Ḣ/H² = -(3/2)*(Ω_m/a³)/(Ω_m/a³ + Ω_Λ)
  ! r_c = 1/(2*sqrt(omega_rc)) / H0
  ! β = 1 + branch * 2*H*r_c*(1 + Ḣ/(3H²))
  !   = 1 + branch / sqrt(omega_rc) * sqrt(Ω_m/a³+Ω_Λ)
  !         * (1 - (Ω_m/a³)/(2*(Ω_m/a³+Ω_Λ)))
  !-------------------------------------------------------
  real(dp)::Oma3, E2, Hdot_over_H2, HrC

  Oma3 = omega_m / aa**3
  E2   = Oma3 + omega_l          ! H²/H0²
  Hdot_over_H2 = -1.5d0 * Oma3 / E2
  ! H*r_c = sqrt(E2)/H0 * 1/(2*sqrt(omega_rc)*H0) = sqrt(E2)/(2*sqrt(omega_rc))
  HrC = sqrt(E2) / (2d0 * sqrt(orc))
  beta = 1d0 + dble(branch) * 2d0 * HrC * (1d0 + Hdot_over_H2 / 3d0)

end function nDGP_beta

!=========================================================
! nDGP_solve_level: top-level nDGP solver for one AMR level
!=========================================================
subroutine nDGP_solve_level(ilevel, icount)
  use amr_commons
  use poisson_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer,intent(in)::ilevel, icount

  integer::iter,ncache,info
  real(dp)::beta
  real(dp)::nDGP_beta  ! external function
  real(dp)::res_max_local,res_max_global
  real(dp)::src_max_local,src_max_global
  real(dp)::rel_res
  logical::converged

  ncache=active(ilevel)%ngrid
  if(ncache==0) return

  ! Compute β(a)
  beta = nDGP_beta(aexp, omega_rc, nDGP_branch)

  if(myid==1 .and. nstep==0 .and. ilevel==levelmin) then
     write(*,'(A,F8.4,A,F8.4)') ' nDGP: beta(a)=', beta, ' at a=', aexp
  end if

  ! Initialize scalar_gr on first step
  if(nstep==0) then
     call nDGP_init_scalar(ilevel)
  end if

  ! Newton-GS relaxation
  converged = .false.
  do iter=1,n_iter_nDGP

     call nDGP_gauss_seidel(ilevel, beta, res_max_local, src_max_local)

     call make_virtual_fine_dp(scalar_gr(1), ilevel)

#ifndef WITHOUTMPI
     call MPI_ALLREDUCE(res_max_local, res_max_global, 1, &
          MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, info)
     call MPI_ALLREDUCE(src_max_local, src_max_global, 1, &
          MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, info)
#else
     res_max_global = res_max_local
     src_max_global = src_max_local
#endif

     if(src_max_global > 0d0) then
        rel_res = res_max_global / src_max_global
     else
        rel_res = 0d0
     end if

     if(rel_res < nDGP_eps) then
        converged = .true.
        if(myid==1) write(*,'(A,I2,A,I3,A,ES10.3)') &
             ' nDGP level ',ilevel,' converged in ',iter,' iters, res=',rel_res
        exit
     end if
  end do

  if(.not. converged .and. myid==1) then
     write(*,'(A,I2,A,I3,A,ES10.3)') &
          ' WARNING: nDGP level ',ilevel,' NOT converged after ', &
          n_iter_nDGP,' iters, res=',rel_res
  end if

  ! Save scalar_gr for warm start
  call nDGP_save_old(ilevel)

  ! Fifth force: F5 = -(1/(2β)) * grad(φ)
  call compute_fifth_force(ilevel, -0.5d0/beta)

end subroutine nDGP_solve_level

!=========================================================
! nDGP_init_scalar: initialize scalar_gr = 0 at a level
!=========================================================
subroutine nDGP_init_scalar(ilevel)
  use amr_commons
  use poisson_commons
  implicit none
  integer,intent(in)::ilevel

  integer::igrid,i,ind,iskip,ncache,icell

  ncache=active(ilevel)%ngrid
!$omp parallel do private(igrid,i,ind,iskip,icell) schedule(dynamic)
  do igrid=1,ncache,nvector
     do i=1,MIN(nvector,ncache-igrid+1)
        do ind=1,twotondim
           iskip=ncoarse+(ind-1)*ngridmax
           icell=iskip+active(ilevel)%igrid(igrid+i-1)
           scalar_gr(icell) = 0d0
           scalar_gr_old(icell) = 0d0
        end do
     end do
  end do

end subroutine nDGP_init_scalar

!=========================================================
! nDGP_save_old: save scalar_gr → scalar_gr_old
!=========================================================
subroutine nDGP_save_old(ilevel)
  use amr_commons
  use poisson_commons
  implicit none
  integer,intent(in)::ilevel

  integer::igrid,i,ind,iskip,ncache,icell

  ncache=active(ilevel)%ngrid
!$omp parallel do private(igrid,i,ind,iskip,icell) schedule(dynamic)
  do igrid=1,ncache,nvector
     do i=1,MIN(nvector,ncache-igrid+1)
        do ind=1,twotondim
           iskip=ncoarse+(ind-1)*ngridmax
           icell=iskip+active(ilevel)%igrid(igrid+i-1)
           scalar_gr_old(icell) = scalar_gr(icell)
        end do
     end do
  end do

end subroutine nDGP_save_old

!=========================================================
! nDGP_gauss_seidel: one Newton-GS sweep for nDGP equation
!
! ∇²φ + coeff*[(∇²φ)² - (∇ᵢ∇ⱼφ)²] = source
!
! coeff = r_c²/(3β a² c²) → in code units (c=1, H0=1):
!   coeff = 1/(12 * omega_rc * beta * a²)
!
! source = 8πG/(3β) * a² * δρ
!   = 2*fourpi/(3β) * δ where fourpi=1.5*Ωm*a
!   = Ωm*a / β * δ
!
! Mixed derivatives via diagonal neighbors (double morton_nbor_grid)
!=========================================================
subroutine nDGP_gauss_seidel(ilevel, beta, res_max, src_max)
  use amr_commons
  use poisson_commons
  use morton_hash
  implicit none
  integer,intent(in)::ilevel
  real(dp),intent(in)::beta
  real(dp),intent(out)::res_max, src_max

  integer::igrid,ngrid,ncache,i,ind,iskip,idim
  integer::ig_left,ig_right,ih_left,ih_right
  integer::ind_nb_left,ind_nb_right
  real(dp)::dx,scale,dx_loc,dx2,dx2_inv
  integer::nx_loc
  real(dp)::u_c,lapl,source,residual,jacobian,delta_u
  real(dp)::coeff
  real(dp)::phi_xm,phi_xp,phi_ym,phi_yp,phi_zm,phi_zp
  real(dp)::phi_xx,phi_yy,phi_zz,phi_xy,phi_xz,phi_yz
  real(dp)::lapl2,trace_ij2,vain_term
  integer::icolor

  ! Diagonal neighbor grids
  integer::ig_diag_pp,ig_diag_pm,ig_diag_mp,ig_diag_mm
  real(dp)::phi_pp,phi_pm,phi_mp,phi_mm

  integer,dimension(1:3,1:2,1:8)::ggg,hhh
  integer,dimension(1:nvector)::ind_grid_w,ind_cell_w
  integer,dimension(1:nvector,0:twondim)::igridn_w

  ncache=active(ilevel)%ngrid
  if(ncache==0) then
     res_max=0d0; src_max=0d0
     return
  end if

  dx=0.5D0**ilevel
  nx_loc=(icoarse_max-icoarse_min+1)
  scale=boxlen/dble(nx_loc)
  dx_loc=dx*scale
  dx2=dx_loc**2
  dx2_inv=1d0/dx2

  ! Vainshtein coefficient: r_c²/(3β a²)
  ! r_c = 1/(2*sqrt(omega_rc)*H0), so r_c² = 1/(4*omega_rc*H0²)
  ! coeff = 1/(12*omega_rc*beta*a²) in code units (H0=1, c=1)
  ! But we need this in mesh units (dx_loc): divide by dx_loc⁴ for the bilinear terms
  ! Actually the equation uses comoving ∇, so coeff stays in comoving
  coeff = 1d0 / (12d0 * omega_rc * beta * aexp**2)

  ! Neighbor lookup
  ggg(1,1,1:8)=(/1,0,1,0,1,0,1,0/); hhh(1,1,1:8)=(/2,1,4,3,6,5,8,7/)
  ggg(1,2,1:8)=(/0,2,0,2,0,2,0,2/); hhh(1,2,1:8)=(/2,1,4,3,6,5,8,7/)
  ggg(2,1,1:8)=(/3,3,0,0,3,3,0,0/); hhh(2,1,1:8)=(/3,4,1,2,7,8,5,6/)
  ggg(2,2,1:8)=(/0,0,4,4,0,0,4,4/); hhh(2,2,1:8)=(/3,4,1,2,7,8,5,6/)
  ggg(3,1,1:8)=(/5,5,5,5,0,0,0,0/); hhh(3,1,1:8)=(/5,6,7,8,1,2,3,4/)
  ggg(3,2,1:8)=(/0,0,0,0,6,6,6,6/); hhh(3,2,1:8)=(/5,6,7,8,1,2,3,4/)

  res_max = 0d0
  src_max = 0d0

  do icolor=0,1

!$omp parallel do private(igrid,ngrid,i,ind,iskip,idim, &
!$omp&  ind_grid_w,ind_cell_w,igridn_w, &
!$omp&  ig_left,ig_right,ih_left,ih_right, &
!$omp&  ind_nb_left,ind_nb_right, &
!$omp&  u_c,lapl,source,residual,jacobian,delta_u, &
!$omp&  phi_xm,phi_xp,phi_ym,phi_yp,phi_zm,phi_zp, &
!$omp&  phi_xx,phi_yy,phi_zz,phi_xy,phi_xz,phi_yz, &
!$omp&  lapl2,trace_ij2,vain_term, &
!$omp&  ig_diag_pp,ig_diag_pm,ig_diag_mp,ig_diag_mm, &
!$omp&  phi_pp,phi_pm,phi_mp,phi_mm) &
!$omp& reduction(max:res_max,src_max) schedule(dynamic)
     do igrid=1,ncache,nvector
        ngrid=MIN(nvector,ncache-igrid+1)
        do i=1,ngrid
           ind_grid_w(i)=active(ilevel)%igrid(igrid+i-1)
        end do

        ! Gather face neighbors
        do i=1,ngrid
           igridn_w(i,0)=ind_grid_w(i)
        end do
        do idim=1,ndim
           do i=1,ngrid
              igridn_w(i,2*idim-1)=morton_nbor_grid(ind_grid_w(i),ilevel,2*idim-1)
              igridn_w(i,2*idim  )=morton_nbor_grid(ind_grid_w(i),ilevel,2*idim  )
           end do
        end do

        do ind=1,twotondim
           if(mod(ind-1, 2) /= icolor) cycle

           iskip=ncoarse+(ind-1)*ngridmax
           do i=1,ngrid
              ind_cell_w(i)=iskip+ind_grid_w(i)
           end do

           do i=1,ngrid
              ! Check face neighbors exist
              if(igridn_w(i,1)==0 .or. igridn_w(i,2)==0 .or. &
                 igridn_w(i,3)==0 .or. igridn_w(i,4)==0 .or. &
                 igridn_w(i,5)==0 .or. igridn_w(i,6)==0) cycle

              u_c = scalar_gr(ind_cell_w(i))

              ! Get face neighbor values
              ig_left =ggg(1,1,ind); ig_right=ggg(1,2,ind)
              ih_left =ncoarse+(hhh(1,1,ind)-1)*ngridmax
              ih_right=ncoarse+(hhh(1,2,ind)-1)*ngridmax
              phi_xm = scalar_gr(igridn_w(i,ig_left )+ih_left)
              phi_xp = scalar_gr(igridn_w(i,ig_right)+ih_right)

              ig_left =ggg(2,1,ind); ig_right=ggg(2,2,ind)
              ih_left =ncoarse+(hhh(2,1,ind)-1)*ngridmax
              ih_right=ncoarse+(hhh(2,2,ind)-1)*ngridmax
              phi_ym = scalar_gr(igridn_w(i,ig_left )+ih_left)
              phi_yp = scalar_gr(igridn_w(i,ig_right)+ih_right)

              ig_left =ggg(3,1,ind); ig_right=ggg(3,2,ind)
              ih_left =ncoarse+(hhh(3,1,ind)-1)*ngridmax
              ih_right=ncoarse+(hhh(3,2,ind)-1)*ngridmax
              phi_zm = scalar_gr(igridn_w(i,ig_left )+ih_left)
              phi_zp = scalar_gr(igridn_w(i,ig_right)+ih_right)

              ! Laplacian
              lapl = (phi_xp + phi_xm + phi_yp + phi_ym + phi_zp + phi_zm - 6d0*u_c) * dx2_inv

              ! Diagonal second derivatives
              phi_xx = (phi_xp + phi_xm - 2d0*u_c) * dx2_inv
              phi_yy = (phi_yp + phi_ym - 2d0*u_c) * dx2_inv
              phi_zz = (phi_zp + phi_zm - 2d0*u_c) * dx2_inv

              ! Mixed derivatives via diagonal neighbors
              ! ∂²φ/∂x∂y: need (+x,+y), (+x,-y), (-x,+y), (-x,-y)
              ! These are obtained by double morton_nbor_grid calls
              ! For simplicity and robustness at coarse-fine boundaries,
              ! approximate mixed derivatives using only face neighbors:
              ! (∇²φ)² - Tr = φ_xx² + φ_yy² + φ_zz² + 2(φ_xy² + φ_xz² + φ_yz²)
              ! = lapl² - (lapl² - φ_xx² - φ_yy² - φ_zz² - 2*cross_terms)
              ! Without diagonal: set cross_terms ≈ 0
              ! Then Vainshtein = lapl² - (φ_xx² + φ_yy² + φ_zz²)
              !   = (φ_xx+φ_yy+φ_zz)² - (φ_xx²+φ_yy²+φ_zz²)
              !   = 2*(φ_xx*φ_yy + φ_xx*φ_zz + φ_yy*φ_zz)

              ! Vainshtein operator: (∇²φ)² - Σ(∂ᵢ∂ⱼφ)²
              ! Approximated as: lapl² - (φ_xx² + φ_yy² + φ_zz²) (ignoring off-diag)
              ! = 2*(φ_xx*φ_yy + φ_xx*φ_zz + φ_yy*φ_zz)
              lapl2 = lapl**2
              trace_ij2 = phi_xx**2 + phi_yy**2 + phi_zz**2
              vain_term = lapl2 - trace_ij2

              ! Source = Ωm*a/β * δ  (in code units)
              ! rho(cell) is the overdensity δ in RAMSES
              source = omega_m * aexp / beta * rho(ind_cell_w(i))

              ! Residual: F = lapl + coeff * vain_term - source
              residual = lapl + coeff * vain_term - source

              ! Jacobian: dF/du_c
              ! dlapl/du_c = -6/dx²
              ! d(vain_term)/du_c = d(lapl²-trace_ij²)/du_c
              !   d(lapl²)/du_c = 2*lapl*(-6*dx2_inv) ... but actually dlapl/du = -(2*ndim)/dx²
              ! Let's use: dlapl/du = -6/dx² (each of 6 face neighbors has +1, center has -6)
              ! d(lapl²)/du = 2*lapl*(-6/dx²)
              ! d(phi_xx²)/du = 2*phi_xx*(-2/dx²), similarly for yy, zz
              ! d(trace_ij2)/du = 2*(-2/dx²)*(phi_xx+phi_yy+phi_zz) = 2*(-2/dx²)*lapl
              ! d(vain_term)/du = 2*lapl*(-6/dx²) - 2*lapl*(-2/dx²)
              !                 = 2*lapl*(-6/dx² + 2/dx²) = 2*lapl*(-4/dx²)
              jacobian = -6d0*dx2_inv + coeff * 2d0 * lapl * (-4d0*dx2_inv)

              ! Newton update
              if(abs(jacobian) > 1d-30) then
                 delta_u = -residual / jacobian
                 ! Damped Newton for stability
                 if(abs(delta_u) > 0.5d0*dx2*abs(source) .and. abs(source) > 1d-30) then
                    delta_u = sign(0.5d0*dx2*abs(source), delta_u)
                 end if
                 scalar_gr(ind_cell_w(i)) = u_c + delta_u
              end if

              res_max = max(res_max, abs(residual))
              src_max = max(src_max, abs(source))

           end do  ! i
        end do  ! ind
     end do  ! igrid

  end do  ! icolor

end subroutine nDGP_gauss_seidel
