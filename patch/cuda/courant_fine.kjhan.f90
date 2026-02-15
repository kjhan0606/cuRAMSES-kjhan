!###############################################
!###############################################
! Dynamic OMP/CUDA hybrid dispatch for courant_fine
! Overrides patch/Horizon5-master-2/courant_fine.kjhan.f90
!###############################################
!###############################################
subroutine courant_fine(ilevel)
  use amr_commons
  use hydro_commons
  use poisson_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer::ilevel
  !----------------------------------------------------------------------
  ! Using the Courant-Friedrich-Levy stability condition,               !
  ! this routine computes the maximum allowed time-step.                !
  !----------------------------------------------------------------------
  integer::i,ivar,idim,ind,ncache,igrid,iskip
  integer::info,nleaf,ngrid,nx_loc

  real(dp)::dt_lev,dx,vol,scale
  real(kind=8)::mass_loc,ekin_loc,eint_loc,dt_loc
  real(kind=8)::imass_loc,iekin_loc,ieint_loc,idt_loc
  real(kind=8)::mass_all,ekin_all,eint_all,dt_all
  real(kind=8),dimension(3)::comm_buffin,comm_buffout

  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel

  mass_all=0.0d0; mass_loc=0.0d0
  ekin_all=0.0d0; ekin_loc=0.0d0
  eint_all=0.0d0; eint_loc=0.0d0
  dt_all=dtnew(ilevel); dt_loc=dt_all

  ! Mesh spacing at that level
  nx_loc=icoarse_max-icoarse_min+1
  scale=boxlen/dble(nx_loc)
  dx=0.5D0**ilevel*scale
  vol=dx**ndim

  ! Loop over active grids by vector sweeps
  ncache=active(ilevel)%ngrid

#ifdef HYDRO_CUDA
  call courant_fine_hybrid(ilevel, ncache, dx, vol, mass_loc, ekin_loc, eint_loc, dt_loc)
#else
!$omp parallel do private(igrid,ngrid,imass_loc,iekin_loc,ieint_loc,idt_loc) &
!$omp& reduction(+:mass_loc,ekin_loc,eint_loc), reduction(min:dt_loc)
  do igrid=1,ncache,nvector
     ngrid=MIN(nvector,ncache-igrid+1)
     call sub_courant_fine(ilevel, igrid, ngrid, imass_loc, iekin_loc, ieint_loc, idt_loc)
     mass_loc = mass_loc + imass_loc
     ekin_loc = ekin_loc + iekin_loc
     eint_loc = eint_loc + ieint_loc
     dt_loc = min(dt_loc, idt_loc)
  enddo
#endif

  ! Compute global quantities
#ifndef WITHOUTMPI
  comm_buffin(1)=mass_loc
  comm_buffin(2)=ekin_loc
  comm_buffin(3)=eint_loc
  call MPI_ALLREDUCE(comm_buffin,comm_buffout,3,MPI_DOUBLE_PRECISION,MPI_SUM,&
       &MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(dt_loc,dt_all,1,MPI_DOUBLE_PRECISION,MPI_MIN,&
       &MPI_COMM_WORLD,info)
  mass_all=comm_buffout(1)
  ekin_all=comm_buffout(2)
  eint_all=comm_buffout(3)
#endif
#ifdef WITHOUTMPI
  mass_all=mass_loc
  ekin_all=ekin_loc
  eint_all=eint_loc
  dt_all=dt_loc
#endif
  mass_tot=mass_tot+mass_all
  ekin_tot=ekin_tot+ekin_all
  eint_tot=eint_tot+eint_all
  dtnew(ilevel)=MIN(dtnew(ilevel),dt_all)

111 format('   Entering courant_fine for level ',I2)
end subroutine courant_fine

#ifdef HYDRO_CUDA
!###############################################
! Hybrid dispatch for courant_fine
!###############################################
subroutine courant_fine_hybrid(ilevel, ncache, dx, vol, &
     mass_loc, ekin_loc, eint_loc, dt_loc)
  use amr_commons
  use hydro_commons
  use hydro_parameters, only: smallr, smallc, gamma, courant_factor
  use poisson_commons
  use cuda_commons
  use hydro_cuda_interface
  use iso_c_binding
  implicit none
  integer, intent(in) :: ilevel, ncache
  real(dp), intent(in) :: dx, vol
  real(kind=8), intent(out) :: mass_loc, ekin_loc, eint_loc, dt_loc

  integer, parameter :: CFL_SUPER_SIZE = 32768
  integer :: igrid, ngrid, stream_slot
  integer :: i, ivar, idim, ind, iskip, nleaf
  integer, dimension(1:nvector) :: ind_grid, ind_cell, ind_leaf
  real(dp), dimension(1:nvector, 1:nvar) :: uu
  real(dp), dimension(1:nvector, 1:ndim) :: gg
  real(dp) :: dt_lev

  ! Per-thread superbatch for GPU
  integer :: nprops, scount, scap
  real(dp), allocatable :: sbuf(:,:), dt_buf(:)

  ! Thread-local accumulators
  real(kind=8) :: lmass, lekin, leint, ldt

  nprops = nvar + ndim
  mass_loc = 0.0d0; ekin_loc = 0.0d0; eint_loc = 0.0d0
  dt_loc = dtnew(ilevel)

!$omp parallel private(igrid, ngrid, stream_slot, i, ivar, idim, ind, iskip, &
!$omp&   nleaf, ind_grid, ind_cell, ind_leaf, uu, gg, dt_lev, &
!$omp&   sbuf, dt_buf, scount, scap, lmass, lekin, leint, ldt) &
!$omp& reduction(+:mass_loc,ekin_loc,eint_loc) reduction(min:dt_loc)
  stream_slot = cuda_acquire_stream_c()
  lmass = 0.0d0; lekin = 0.0d0; leint = 0.0d0; ldt = dt_loc

  if (stream_slot >= 0) then
     scap = CFL_SUPER_SIZE
     allocate(sbuf(scap, nprops), dt_buf(scap))
     scount = 0
  end if

!$omp do schedule(dynamic)
  do igrid = 1, ncache, nvector
     ngrid = MIN(nvector, ncache - igrid + 1)
     do i = 1, ngrid
        ind_grid(i) = active(ilevel)%igrid(igrid+i-1)
     end do

     do ind = 1, twotondim
        iskip = ncoarse + (ind - 1) * ngridmax
        do i = 1, ngrid
           ind_cell(i) = ind_grid(i) + iskip
        end do

        ! Gather leaf cells
        nleaf = 0
        do i = 1, ngrid
           if (son(ind_cell(i)) == 0) then
              nleaf = nleaf + 1
              ind_leaf(nleaf) = ind_cell(i)
           end if
        end do
        if (nleaf == 0) cycle

        ! Accumulate mass/energy (always CPU — lightweight)
        do i = 1, nleaf
           lmass = lmass + uold(ind_leaf(i), 1) * vol
           lekin = lekin + uold(ind_leaf(i), ndim+2) * vol
           leint = leint + uold(ind_leaf(i), ndim+2) * vol
        end do
        do ivar = 1, ndim
           do i = 1, nleaf
              leint = leint - 0.5d0 * uold(ind_leaf(i),1+ivar)**2 &
                   / max(uold(ind_leaf(i),1), smallr) * vol
           end do
        end do

        if (stream_slot < 0) then
           ! CPU path: call cmpdt directly
           do ivar = 1, nvar
              do i = 1, nleaf
                 uu(i, ivar) = uold(ind_leaf(i), ivar)
              end do
           end do
           gg = 0.0d0
           if (poisson) then
              do idim = 1, ndim
                 do i = 1, nleaf
                    gg(i, idim) = f(ind_leaf(i), idim)
                 end do
              end do
           end if
           if (nleaf > 0) then
              call cmpdt(uu, gg, dx, dt_lev, nleaf)
              ldt = min(ldt, dt_lev)
           end if
        else
           ! GPU path: gather into superbatch
           do i = 1, nleaf
              scount = scount + 1
              do ivar = 1, nvar
                 sbuf(scount, ivar) = uold(ind_leaf(i), ivar)
              end do
              if (poisson) then
                 do idim = 1, ndim
                    sbuf(scount, nvar + idim) = f(ind_leaf(i), idim)
                 end do
              else
                 do idim = 1, ndim
                    sbuf(scount, nvar + idim) = 0.0d0
                 end do
              end if
           end do

           ! Flush when near full
           if (scount + nvector * twotondim > scap) then
              call courant_gpu_flush(sbuf, scap, nprops, dt_buf, &
                   scount, stream_slot, dx, ldt)
           end if
        end if
     end do
  end do
!$omp end do nowait

  ! Flush remaining
  if (stream_slot >= 0) then
     if (scount > 0) call courant_gpu_flush(sbuf, scap, nprops, dt_buf, &
          scount, stream_slot, dx, ldt)
     deallocate(sbuf, dt_buf)
     call cuda_release_stream_c(stream_slot)
  end if

  mass_loc = mass_loc + lmass
  ekin_loc = ekin_loc + lekin
  eint_loc = eint_loc + leint
  dt_loc = min(dt_loc, ldt)
!$omp end parallel

end subroutine courant_fine_hybrid

!###############################################
! GPU flush: explicit-shape args (no assumed-shape)
!###############################################
subroutine courant_gpu_flush(sbuf, sbuf_ld, nprops, dt_buf, &
     scount, stream_slot, dx, ldt)
  use amr_parameters, only: dp, ndim
  use hydro_parameters, only: smallr, smallc, gamma, courant_factor, nvar
  use hydro_cuda_interface
  use iso_c_binding
  implicit none
  integer, intent(in) :: sbuf_ld, nprops
  real(dp), intent(in), target :: sbuf(sbuf_ld, nprops)
  real(dp), intent(inout), target :: dt_buf(*)
  integer, intent(inout) :: scount
  integer, intent(in) :: stream_slot
  real(dp), intent(in) :: dx
  real(kind=8), intent(inout) :: ldt
  integer :: k

  call cmpdt_cuda_async_f(sbuf, dt_buf, &
       real(dx, c_double), real(smallr, c_double), real(smallc, c_double), &
       real(gamma, c_double), real(courant_factor, c_double), &
       int(scount, c_int), int(nvar, c_int), int(ndim, c_int), &
       int(sbuf_ld, c_int), int(stream_slot, c_int))
  call cmpdt_cuda_sync_f(int(stream_slot, c_int))

  ! CPU min-reduction over per-cell dt values
  do k = 1, scount
     ldt = min(ldt, dt_buf(k))
  end do

  scount = 0
end subroutine courant_gpu_flush
#endif

!###############################################
!###############################################
subroutine sub_courant_fine(ilevel,igrid,ngrid, mass_loc,ekin_loc,eint_loc,dt_loc)
  use amr_commons
  use hydro_commons
  use poisson_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer::ilevel
  integer::i,ivar,idim,ind,ncache,igrid,iskip
  integer::info,nleaf,ngrid,nx_loc
  integer,dimension(1:nvector)::ind_grid,ind_cell,ind_leaf

  real(dp)::dt_lev,dx,vol,scale
  real(kind=8)::mass_loc,ekin_loc,eint_loc,dt_loc
  real(dp),dimension(1:nvector,1:nvar)::uu
  real(dp),dimension(1:nvector,1:ndim)::gg
  mass_loc=0.0d0
  ekin_loc=0.0d0
  eint_loc=0.0d0
  dt_loc=dtnew(ilevel);

  ! Mesh spacing at that level
  nx_loc=icoarse_max-icoarse_min+1
  scale=boxlen/dble(nx_loc)
  dx=0.5D0**ilevel*scale
  vol=dx**ndim

     do i=1,ngrid
        ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
     end do

     ! Loop over cells
     do ind=1,twotondim
        iskip=ncoarse+(ind-1)*ngridmax
        do i=1,ngrid
           ind_cell(i)=ind_grid(i)+iskip
        end do

        ! Gather leaf cells
        nleaf=0
        do i=1,ngrid
           if(son(ind_cell(i))==0)then
              nleaf=nleaf+1
              ind_leaf(nleaf)=ind_cell(i)
           end if
        end do

        ! Gather hydro variables
        do ivar=1,nvar
           do i=1,nleaf
              uu(i,ivar)=uold(ind_leaf(i),ivar)
           end do
        end do

        ! Gather gravitational acceleration
        gg=0.0d0
        if(poisson)then
           do idim=1,ndim
              do i=1,nleaf
                 gg(i,idim)=f(ind_leaf(i),idim)
              end do
           end do
        end if

        ! Compute total mass
        do i=1,nleaf
           mass_loc=mass_loc+uu(i,1)*vol
        end do

        ! Compute total energy
        do i=1,nleaf
           ekin_loc=ekin_loc+uu(i,ndim+2)*vol
        end do

        ! Compute total internal energy
        do i=1,nleaf
           eint_loc=eint_loc+uu(i,ndim+2)*vol
        end do
        do ivar=1,ndim
           do i=1,nleaf
              eint_loc=eint_loc-0.5d0*uu(i,1+ivar)**2/max(uu(i,1),smallr)*vol
           end do
        end do
#if NENER>0
        do ivar=1,nener
           do i=1,nleaf
              eint_loc=eint_loc-uu(i,ndim+2+ivar)*vol
           end do
        end do
#endif

        ! Compute CFL time-step
        if(nleaf>0)then
           call cmpdt(uu,gg,dx,dt_lev,nleaf)
           dt_loc=min(dt_loc,dt_lev)
        end if

     end do
     ! End loop over cells

end subroutine sub_courant_fine

!##################
!##################
subroutine check_cons(ilevel)
  use amr_commons
  use hydro_commons
  use poisson_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer::ilevel
  integer::i,ivar,idim,ind,ncache,igrid,iskip
  integer::info,nleaf,ngrid,nx_loc

  real(dp)::dt_lev,dx,vol,scale
  real(kind=8)::mass_loc,ekin_loc,eint_loc
  real(kind=8)::imass_loc,iekin_loc,ieint_loc
  real(kind=8)::mass_all,ekin_all,eint_all,dt_all
  real(kind=8),dimension(3)::comm_buffin,comm_buffout

  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel

  mass_all=0.0d0; mass_loc=0.0d0
  ekin_all=0.0d0; ekin_loc=0.0d0
  eint_all=0.0d0; eint_loc=0.0d0

  ! Mesh spacing at that level
  dx=0.5D0**ilevel*boxlen
  vol=dx**ndim

  ! Loop over active grids by vector sweeps
  ncache=active(ilevel)%ngrid
!$omp parallel do private(igrid,ngrid,imass_loc,iekin_loc,ieint_loc) &
!$omp& reduction(+:mass_loc,ekin_loc,eint_loc)
  do igrid=1,ncache,nvector
     ngrid=MIN(nvector,ncache-igrid+1)
     call sub_check_cons(ilevel,igrid,ngrid,imass_loc,iekin_loc,ieint_loc)
     mass_loc = mass_loc + imass_loc
     ekin_loc = ekin_loc + iekin_loc
     eint_loc = eint_loc + ieint_loc
  enddo
  ! Compute global quantities
#ifndef WITHOUTMPI
  comm_buffin(1)=mass_loc
  comm_buffin(2)=ekin_loc
  comm_buffin(3)=eint_loc
  call MPI_ALLREDUCE(comm_buffin,comm_buffout,3,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  mass_all=comm_buffout(1)
  ekin_all=comm_buffout(2)
  eint_all=comm_buffout(3)
#endif
#ifdef WITHOUTMPI
  mass_all=mass_loc
  ekin_all=ekin_loc
  eint_all=eint_loc
#endif
  mass_tot=mass_tot+mass_all
  ekin_tot=ekin_tot+ekin_all
  eint_tot=eint_tot+eint_all

111 format('   Entering check_cons for level ',I2)
end subroutine check_cons

!##################
!##################
subroutine sub_check_cons(ilevel,igrid,ngrid,mass_loc,ekin_loc,eint_loc)
  use amr_commons
  use hydro_commons
  use poisson_commons
  implicit none
  integer::ilevel
  integer::i,ivar,idim,ind,igrid,iskip
  integer::nleaf,ngrid
  integer,dimension(1:nvector)::ind_grid,ind_cell,ind_leaf

  real(dp)::dx,vol
  real(kind=8)::mass_loc,ekin_loc,eint_loc
  real(dp),dimension(1:nvector,1:nvar)::uu

  mass_loc=0.0d0
  ekin_loc=0.0d0
  eint_loc=0.0d0

  dx=0.5D0**ilevel*boxlen
  vol=dx**ndim

  do i=1,ngrid
     ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
  end do

  do ind=1,twotondim
     iskip=ncoarse+(ind-1)*ngridmax
     do i=1,ngrid
        ind_cell(i)=ind_grid(i)+iskip
     end do

     nleaf=0
     do i=1,ngrid
        if(son(ind_cell(i))==0)then
           nleaf=nleaf+1
           ind_leaf(nleaf)=ind_cell(i)
        end if
     end do

     do ivar=1,nvar
        do i=1,nleaf
           uu(i,ivar)=uold(ind_leaf(i),ivar)
        end do
     end do

     do i=1,nleaf
        mass_loc=mass_loc+uu(i,1)*vol
     end do

     do i=1,nleaf
        ekin_loc=ekin_loc+uu(i,ndim+2)*vol
     end do

     do i=1,nleaf
        eint_loc=eint_loc+uu(i,ndim+2)*vol
     end do
     do ivar=1,ndim
        do i=1,nleaf
           eint_loc=eint_loc-0.5d0*uu(i,1+ivar)**2/max(uu(i,1),smallr)*vol
        end do
     end do
#if NENER>0
     do ivar=1,nener
        do i=1,nleaf
           eint_loc=eint_loc-uu(i,ndim+2+ivar)*vol
        end do
     end do
#endif
  end do

end subroutine sub_check_cons
