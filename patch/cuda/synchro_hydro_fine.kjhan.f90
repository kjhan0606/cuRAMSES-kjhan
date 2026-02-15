!################################################################
!################################################################
! Dynamic OMP/CUDA hybrid dispatch for synchro_hydro_fine
! Overrides patch/Horizon5-master-2/synchro_hydro_fine.kjhan.f90
!################################################################
!################################################################
subroutine synchro_hydro_fine(ilevel,dteff)
  use amr_commons
  use hydro_commons
  implicit none
  integer::ilevel
  real(dp)::dteff
  !-------------------------------------------------------------------
  ! Update velocity from gravitational acceleration
  !-------------------------------------------------------------------
  integer::ncache,ngrid,i,igrid,iskip,ind

  if(.not. poisson)return
  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel

  ncache=active(ilevel)%ngrid

#ifdef HYDRO_CUDA
  call synchro_hydro_hybrid(ilevel, ncache, dteff)
#else
!$omp parallel do private(igrid,ngrid,i,ind,iskip)
  do igrid=1,ncache,nvector
     ngrid=MIN(nvector,ncache-igrid+1)
     call sub_synchro_hydro_fine(ilevel, igrid, ngrid, dteff)
  end do
#endif

111 format('   Entering synchro_hydro_fine for level',i2)
end subroutine synchro_hydro_fine

!################################################################
!################################################################
subroutine sub_synchro_hydro_fine(ilevel, igrid, ngrid, dteff)
  use amr_commons
  use hydro_commons
  implicit none
  integer, intent(in) :: ilevel, igrid, ngrid
  real(dp), intent(in) :: dteff
  integer :: i, ind, iskip
  integer, dimension(1:nvector) :: ind_grid, ind_cell

  do i=1,ngrid
     ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
  end do
  do ind=1,twotondim
     iskip=ncoarse+(ind-1)*ngridmax
     do i=1,ngrid
        ind_cell(i)=ind_grid(i)+iskip
     end do
     call synchydrofine1(ind_cell,ngrid,dteff)
  end do
end subroutine sub_synchro_hydro_fine

#ifdef HYDRO_CUDA
!################################################################
! Hybrid dispatch: OMP threads with dynamic GPU/CPU fallback
!################################################################
subroutine synchro_hydro_hybrid(ilevel, ncache, dteff)
  use amr_commons
  use hydro_commons
  use poisson_commons
  use hydro_parameters, only: smallr
  use cuda_commons
  use hydro_cuda_interface
  use iso_c_binding
  implicit none
  integer, intent(in) :: ilevel, ncache
  real(dp), intent(in) :: dteff

  integer, parameter :: SYNC_SUPER_SIZE = 16384
  integer :: igrid, ngrid, stream_slot
  integer :: i, ind, iskip, nprops
  integer, dimension(1:nvector) :: ind_grid, ind_cell

  ! Superbatch buffers (per-thread via OMP private)
  real(dp), allocatable :: sbuf(:,:)
  integer, allocatable :: sidx(:)
  integer :: scount, scap

  nprops = ndim + 2 + ndim  ! rho + mom(ndim) + E + f(ndim)

!$omp parallel private(igrid, ngrid, stream_slot, i, ind, iskip, &
!$omp&   ind_grid, ind_cell, sbuf, sidx, scount, scap)
  stream_slot = cuda_acquire_stream_c()

  if (stream_slot >= 0) then
     scap = SYNC_SUPER_SIZE
     allocate(sbuf(scap, nprops), sidx(scap))
     scount = 0
  end if

!$omp do schedule(dynamic)
  do igrid = 1, ncache, nvector
     ngrid = MIN(nvector, ncache - igrid + 1)

     if (stream_slot < 0) then
        ! CPU path
        call sub_synchro_hydro_fine(ilevel, igrid, ngrid, dteff)
     else
        ! GPU path: gather cells into superbatch
        do i = 1, ngrid
           ind_grid(i) = active(ilevel)%igrid(igrid+i-1)
        end do
        do ind = 1, twotondim
           iskip = ncoarse + (ind - 1) * ngridmax
           do i = 1, ngrid
              scount = scount + 1
              ind_cell(i) = ind_grid(i) + iskip
              sidx(scount) = ind_cell(i)
              ! Gather uold(1:ndim+2)
              sbuf(scount, 1) = uold(ind_cell(i), 1)
              sbuf(scount, 2) = uold(ind_cell(i), 2)
              sbuf(scount, 3) = uold(ind_cell(i), 3)
              sbuf(scount, 4) = uold(ind_cell(i), 4)
              sbuf(scount, 5) = uold(ind_cell(i), 5)
              ! Gather f(1:ndim)
              sbuf(scount, 6) = f(ind_cell(i), 1)
              sbuf(scount, 7) = f(ind_cell(i), 2)
              sbuf(scount, 8) = f(ind_cell(i), 3)
           end do
        end do

        ! Flush when near full
        if (scount + nvector * twotondim > scap) then
           call synchro_gpu_flush(sbuf, scap, nprops, sidx, scount, &
                stream_slot, dteff)
        end if
     end if
  end do
!$omp end do nowait

  ! Flush remaining
  if (stream_slot >= 0) then
     if (scount > 0) call synchro_gpu_flush(sbuf, scap, nprops, sidx, &
          scount, stream_slot, dteff)
     deallocate(sbuf, sidx)
     call cuda_release_stream_c(stream_slot)
  end if
!$omp end parallel

end subroutine synchro_hydro_hybrid

!################################################################
! GPU flush: explicit-shape args (no assumed-shape, no descriptor)
!################################################################
subroutine synchro_gpu_flush(sbuf, sbuf_ld, nprops, sidx, scount, &
     stream_slot, dteff)
  use amr_commons, only: dp, ndim
  use hydro_commons, only: uold
  use hydro_parameters, only: smallr
  use hydro_cuda_interface
  use iso_c_binding
  implicit none
  integer, intent(in) :: sbuf_ld, nprops
  real(dp), intent(inout), target :: sbuf(sbuf_ld, nprops)
  integer, intent(in) :: sidx(*)
  integer, intent(inout) :: scount
  integer, intent(in) :: stream_slot
  real(dp), intent(in) :: dteff
  integer :: k, idim

  call synchro_cuda_async_f(sbuf, &
       real(smallr, c_double), real(dteff, c_double), &
       int(scount, c_int), int(ndim, c_int), &
       int(sbuf_ld, c_int), int(stream_slot, c_int))
  call synchro_cuda_sync_f(int(stream_slot, c_int))

  ! Scatter updated values back
  do k = 1, scount
     do idim = 1, ndim
        uold(sidx(k), idim + 1) = sbuf(k, idim + 1)
     end do
     uold(sidx(k), ndim + 2) = sbuf(k, ndim + 2)
  end do

  scount = 0
end subroutine synchro_gpu_flush
#endif

!################################################################
!################################################################
subroutine synchydrofine1(ind_cell,ncell,dteff)
  use amr_commons
  use hydro_commons
  use poisson_commons
  implicit none
  integer::ncell
  real(dp)::dteff
  integer,dimension(1:nvector)::ind_cell
  !-------------------------------------------------------------------
  ! Gravity update for hydro variables
  !-------------------------------------------------------------------
  integer::i,idim,neul=ndim+2,nndim=ndim
  real(dp),dimension(1:nvector)::pp

  ! Compute internal + magnetic + radiative energy
  do i=1,ncell
     pp(i)=uold(ind_cell(i),neul)
  end do
  do idim=1,nndim
     do i=1,ncell
        pp(i)=pp(i)-0.5*uold(ind_cell(i),idim+1)**2/max(uold(ind_cell(i),1),smallr)
     end do
  end do
  do i=1,ncell
     uold(ind_cell(i),neul)=pp(i)
  end do

  ! Update momentum
  do idim=1,ndim
     do i=1,ncell
        pp(i)=uold(ind_cell(i),idim+1)+ &
             & max(uold(ind_cell(i),1),smallr)*f(ind_cell(i),idim)*dteff
     end do
     do i=1,ncell
        uold(ind_cell(i),idim+1)=pp(i)
     end do
  end do

  ! Update total energy
  do i=1,ncell
     pp(i)=uold(ind_cell(i),neul)
  end do
  do idim=1,nndim
     do i=1,ncell
        pp(i)=pp(i)+0.5*uold(ind_cell(i),idim+1)**2/max(uold(ind_cell(i),1),smallr)
     end do
  end do
  do i=1,ncell
     uold(ind_cell(i),neul)=pp(i)
  end do

end subroutine synchydrofine1
