! ------------------------------------------------------------------------
! Multigrid Poisson solver for refined AMR levels
! ------------------------------------------------------------------------
! This file contains all generic fine multigrid routines, such as
!   * multigrid iterations @ fine and coarse MG levels
!   * communicator building
!   * MPI routines
!   * helper functions
!
! Used variables:
!                       finest(AMR)level     coarse(MG)levels
!     -----------------------------------------------------------------
!     potential            phi            active_mg(myid,ilevel)%u(:,1)
!     physical RHS         rho            active_mg(myid,ilevel)%u(:,2)
!     residual             f(:,1)         active_mg(myid,ilevel)%u(:,3)
!     BC-modified RHS      f(:,2)                  N/A
!     mask                 f(:,3)         active_mg(myid,ilevel)%u(:,4)
!
! ------------------------------------------------------------------------


! ------------------------------------------------------------------------
! Main multigrid routine, called by amr_step
! ------------------------------------------------------------------------

subroutine multigrid_fine(ilevel,icount)
   use amr_commons
   use poisson_commons
   use poisson_parameters
#ifdef HYDRO_CUDA
   use poisson_cuda_interface
   use iso_c_binding
   use cuda_commons, only: cuda_pool_is_initialized_c
#endif

   implicit none
#ifndef WITHOUTMPI
   include "mpif.h"
#endif

   integer, intent(in) :: ilevel,icount

   interface
      subroutine cmp_residual_mg_fine(ilevel, norm2)
         use amr_commons, only: dp
         integer, intent(in) :: ilevel
         real(dp), intent(out), optional :: norm2
      end subroutine
   end interface

   integer, parameter  :: MAXITER  = 10
   real(dp), parameter :: SAFE_FACTOR = 0.5

   integer  :: ifine, i, iter, info, icpu
   real(kind=8) :: res_norm2, i_res_norm2, i_res_norm2_tot, res_norm2_tot
   real(kind=8) :: debug_norm2, debug_norm2_tot
   real(kind=8) :: err, last_err

   logical :: allmasked, allmasked_tot

   ! GPU Poisson MG variables
   logical :: use_mg_gpu, use_ri_gpu
#ifdef HYDRO_CUDA
   integer(c_long_long) :: ncell_tot_c
   integer :: ncell_tot, safe_int, ri_flag
   real(kind=8) :: dx_mg, dx2_mg, oneoverdx2_mg, dx2_norm_mg
   real(kind=8) :: gpu_norm2, dummy_norm2
   ! cuFFT variables
   logical :: is_uniform_fft
   integer(i8b) :: expected_grids
#endif

   if(gravity_type>0)return
   if(numbtot(1,ilevel)==0)return

   if(verbose) print '(A,I2)','Entering fine multigrid at level ',ilevel

   ! ---------------------------------------------------------------------
   ! Prepare first guess, mask and BCs at finest level
   ! ---------------------------------------------------------------------

   if(ilevel>levelmin)then
      call make_initial_phi(ilevel,icount)         ! Interpolate phi down
   else
      call make_multipole_phi(ilevel)       ! Fill with simple initial guess
   endif
   call make_virtual_fine_dp(phi(1),ilevel) ! Update boundaries
   call make_boundary_phi(ilevel)           ! Update physical boundaries

   call make_fine_mask  (ilevel)            ! Fill the fine mask
   call make_virtual_fine_dp(f(:,3),ilevel) ! Communicate mask
   call make_boundary_mask(ilevel)          ! Set mask to -1 in phys bounds

   ! Pre-compute neighbor grids BEFORE bc_rhs so it can use the array
   call precompute_nbor_grid_fine(ilevel)

   call make_fine_bc_rhs(ilevel,icount)            ! Fill BC-modified RHS

   ! ---------------------------------------------------------------------
   ! Build communicators up
   ! ---------------------------------------------------------------------

   ! @ finer level
   call build_parent_comms_mg(active(ilevel),ilevel)
   ! @ coarser levels
   do ifine=(ilevel-1),2,-1
      call build_parent_comms_mg(active_mg(myid,ifine),ifine)
   end do

   ! ---------------------------------------------------------------------
   ! Restrict mask up, then set scan flag
   ! ---------------------------------------------------------------------
   ! @ finer level

   if(ilevel>1) then
      ! Restrict and communicate mask
      call restrict_mask_fine_reverse(ilevel)
      call make_reverse_mg_dp(4,ilevel-1)
      call make_virtual_mg_dp(4,ilevel-1)

      ! Convert volume fraction to mask value
      do icpu=1,ncpu
         if(active_mg(icpu,ilevel-1)%ngrid==0) cycle
         active_mg(icpu,ilevel-1)%u(:,4)=2d0*active_mg(icpu,ilevel-1)%u(:,4)-1d0
      end do

      ! Check active mask state
      if(active_mg(myid,ilevel-1)%ngrid>0) then
         allmasked=(maxval(active_mg(myid,ilevel-1)%u(:,4))<=0d0)
      else
         allmasked=.true.
      end if

      ! Allreduce on mask state
#ifndef WITHOUTMPI
      call MPI_ALLREDUCE(allmasked, allmasked_tot, 1, MPI_LOGICAL, &
           & MPI_LAND, MPI_COMM_WORLD, info)
      allmasked=allmasked_tot
#endif
   else
      allmasked=.true.
   endif

   ! @ coarser levels
   ! Restrict mask and compute levelmin_mg in the process
   if (.not. allmasked) then
      levelmin_mg=1
      do ifine=(ilevel-1),2,-1

         ! Restrict and communicate mask
         call restrict_mask_coarse_reverse(ifine)
         call make_reverse_mg_dp(4,ifine-1)
         call make_virtual_mg_dp(4,ifine-1)

         ! Convert volume fraction to mask value
         do icpu=1,ncpu
            if(active_mg(icpu,ifine-1)%ngrid==0) cycle
            active_mg(icpu,ifine-1)%u(:,4)=2d0*active_mg(icpu,ifine-1)%u(:,4)-1d0
         end do

         ! Check active mask state
         if(active_mg(myid,ifine-1)%ngrid>0) then
            allmasked=(maxval(active_mg(myid,ifine-1)%u(:,4))<=0d0)
         else
            allmasked=.true.
         end if

         ! Allreduce on mask state
#ifndef WITHOUTMPI
         call MPI_ALLREDUCE(allmasked,allmasked_tot,1,MPI_LOGICAL, &
                 & MPI_LAND,MPI_COMM_WORLD,info)
         allmasked=allmasked_tot
#endif

         if(allmasked) then ! Coarser level is fully masked: stop here
            levelmin_mg=ifine
            exit
         end if
      end do
   else
      levelmin_mg=ilevel
   end if
   if(nboundary>0)levelmin_mg=max(levelmin_mg,2)

   ! nbor_grid_fine already precomputed before make_fine_bc_rhs

   ! Update flag with scan flag (uses nbor_grid_fine if available)
   call set_scan_flag_fine(ilevel)
   do ifine=levelmin_mg,ilevel-1
      call set_scan_flag_coarse(ifine)
   end do

   ! Precompute neighbor grids for coarse MG levels
   if(ilevel>1 .and. levelmin_mg < ilevel) then
      call precompute_nbor_grid_coarse(levelmin_mg, ilevel-1)
   end if

   ! ---------------------------------------------------------------------
   ! Initiate solve at fine level
   ! ---------------------------------------------------------------------

   ! GPU setup: upload MG arrays to GPU if CUDA available
   use_mg_gpu = .false.
   use_ri_gpu = .false.
#ifdef HYDRO_CUDA
   if(cuda_pool_is_initialized_c() /= 0) then
      ncell_tot = ncoarse + twotondim*ngridmax
      ncell_tot_c = int(ncell_tot, c_long_long)
      dx_mg  = 0.5d0**ilevel
      dx2_mg = dx_mg*dx_mg
      oneoverdx2_mg = 1.0d0/dx2_mg
      dx2_norm_mg = dx_mg**(ndim)
      call cuda_mg_upload_c( &
           phi, f, flag2, ncell_tot_c, &
           nbor_grid_fine, active(ilevel)%igrid, &
           int(active(ilevel)%ngrid, c_int))
      use_mg_gpu = (cuda_mg_is_ready_c() /= 0)
      if(use_mg_gpu) call build_mg_halo_indices(ilevel)
      ! Setup GPU restrict/interp if MG GPU is ready and coarse levels exist
      if(use_mg_gpu .and. ilevel > 1) then
         call precompute_mg_gpu_restrict_interp(ilevel)
         use_ri_gpu = (cuda_mg_ri_is_ready_c() /= 0)
      end if
      ! Synchronize use_ri_gpu across all ranks (MPI collective paths must match)
#ifndef WITHOUTMPI
      ri_flag = 0; if(use_ri_gpu) ri_flag = 1
      call MPI_ALLREDUCE(MPI_IN_PLACE, ri_flag, 1, &
           MPI_INTEGER, MPI_MIN, MPI_COMM_WORLD, info)
      use_ri_gpu = (ri_flag == 1)
#endif
      if(myid==1) write(*,'(A,I3,A,L1,A,L1,A,I12,A,I10)') &
           ' MG GPU: level=',ilevel,' ready=',use_mg_gpu, &
           ' ri=',use_ri_gpu, &
           ' ncell=',ncell_tot,' ngrid=',active(ilevel)%ngrid
   end if

   ! -----------------------------------------------------------------
   ! cuFFT direct solve for fully uniform level (periodic BC)
   ! -----------------------------------------------------------------
   is_uniform_fft = .false.
   if(use_mg_gpu) then
      expected_grids = int(nx,i8b)*int(ny,i8b)*int(nz,i8b) &
                     * 8_i8b**(ilevel-1)
      if(numbtot(1,ilevel) == expected_grids) is_uniform_fft = .true.
   end if

   if(is_uniform_fft) then
      call fft_poisson_solve_uniform(ilevel, icount)
      ! Download phi from GPU
      call cuda_mg_download_phi_c(phi, ncell_tot_c)
      call make_virtual_fine_dp(phi(1), ilevel)
      ! Cleanup GPU state
      call cuda_mg_halo_free_c()
      call cuda_mg_free_c()
      if(allocated(mg_halo_emit_cells)) deallocate(mg_halo_emit_cells)
      if(allocated(mg_halo_recv_cells)) deallocate(mg_halo_recv_cells)
      if(allocated(mg_halo_emit_buf))   deallocate(mg_halo_emit_buf)
      if(allocated(mg_halo_recv_buf))   deallocate(mg_halo_recv_buf)
      mg_halo_n_emit = 0
      mg_halo_n_recv = 0
      if(use_ri_gpu) call cuda_mg_ri_free_c()
      if(allocated(mg_ri_flat_offset)) deallocate(mg_ri_flat_offset)
      if(allocated(mg_ri_coarse_rhs))  deallocate(mg_ri_coarse_rhs)
      if(allocated(mg_ri_coarse_phi))  deallocate(mg_ri_coarse_phi)
      mg_ri_total_coarse = 0
      ! Free precomputed neighbors
      if(allocated(nbor_grid_fine)) deallocate(nbor_grid_fine)
      if(ilevel>1 .and. levelmin_mg < ilevel) then
         call cleanup_nbor_grid_coarse(levelmin_mg, ilevel-1)
      end if
      ! Cleanup MG levels
      do ifine=1,ilevel-1
         call cleanup_mg_level(ifine)
      end do
      if(myid==1) write(*,'(A,I5,A)') &
           '   ==> Level=',ilevel,' cuFFT direct solve DONE'
      return
   end if
#endif

   iter = 0
   err = 1.0d0
   main_iteration_loop: do
      iter=iter+1

      ! Pre-smoothing
      do i=1,ngs_fine
#ifdef HYDRO_CUDA
         if(use_mg_gpu) then
            safe_int = 0
            if(safe_mode(ilevel)) safe_int = 1
            call cuda_mg_gauss_seidel_c(int(active(ilevel)%ngrid,c_int), &
                 int(ngridmax,c_int), int(ncoarse,c_int), dx2_mg, 0, safe_int)
            call cuda_mg_gauss_seidel_c(int(active(ilevel)%ngrid,c_int), &
                 int(ngridmax,c_int), int(ncoarse,c_int), dx2_mg, 1, safe_int)
            call make_virtual_fine_dp_gpu(ilevel)
         else
#endif
            call gauss_seidel_mg_fine(ilevel,.true. )  ! Red step
            call gauss_seidel_mg_fine(ilevel,.false.)  ! Black step
            call make_virtual_fine_dp(phi(1),ilevel)   ! Communicate phi
#ifdef HYDRO_CUDA
         end if
#endif
      end do

      ! Compute residual and restrict into upper level RHS
#ifdef HYDRO_CUDA
      if(use_mg_gpu) then
         gpu_norm2 = 0.0d0
         if(iter==1) then
            call cuda_mg_residual_c(int(active(ilevel)%ngrid,c_int), &
                 int(ngridmax,c_int), int(ncoarse,c_int), &
                 oneoverdx2_mg, dble(twondim), dx2_norm_mg, &
                 gpu_norm2, 1)
            i_res_norm2 = gpu_norm2
#ifndef WITHOUTMPI
            call MPI_ALLREDUCE(i_res_norm2,i_res_norm2_tot,1, &
                    & MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
            i_res_norm2=i_res_norm2_tot
#endif
         else
            call cuda_mg_residual_c(int(active(ilevel)%ngrid,c_int), &
                 int(ngridmax,c_int), int(ncoarse,c_int), &
                 oneoverdx2_mg, dble(twondim), dx2_norm_mg, &
                 dummy_norm2, 0)
         end if
         if(.not. use_ri_gpu) call cuda_mg_download_f1_c(f, ncell_tot_c)
      else
#endif
         if(iter==1) then
            call cmp_residual_mg_fine(ilevel, i_res_norm2)
#ifndef WITHOUTMPI
            call MPI_ALLREDUCE(i_res_norm2,i_res_norm2_tot,1, &
                    & MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
            i_res_norm2=i_res_norm2_tot
#endif
         else
            call cmp_residual_mg_fine(ilevel)
         end if
#ifdef HYDRO_CUDA
      end if
#endif

      ! Restrict residual into upper level RHS
#ifdef HYDRO_CUDA
      if(use_ri_gpu) then
         ! GPU restriction: d_mg_f1 → d_coarse_rhs_flat → host → active_mg
         call cuda_mg_restrict_execute_c(int(active(ilevel)%ngrid,c_int), &
              int(ngridmax,c_int), int(ncoarse,c_int))
         call cuda_mg_restrict_download_c(mg_ri_coarse_rhs, &
              int(mg_ri_total_coarse,c_int))
         call scatter_coarse_rhs_from_flat(ilevel)
      else
#endif
         ! First clear the rhs in coarser reception comms
         do icpu=1,ncpu
            if(active_mg(icpu,ilevel-1)%ngrid==0) cycle
            active_mg(icpu,ilevel-1)%u(:,2)=0.0d0
         end do
         ! Restrict
         call restrict_residual_fine_reverse(ilevel)
#ifdef HYDRO_CUDA
      end if
#endif
      call make_reverse_mg_dp(2,ilevel-1) ! communicate rhs

      if(ilevel>1) then
         ! Reset correction at upper level before solve
         do icpu=1,ncpu
            if(active_mg(icpu,ilevel-1)%ngrid==0) cycle
            active_mg(icpu,ilevel-1)%u(:,1)=0.0d0
         end do

         ! Multigrid-solve the upper level
         call recursive_multigrid_coarse(ilevel-1, safe_mode(ilevel))

         ! Interpolate coarse solution and correct fine solution
#ifdef HYDRO_CUDA
         if(use_ri_gpu) then
            ! GPU interpolation: active_mg → flat → GPU → correct d_mg_phi
            call gather_coarse_phi_to_flat(ilevel)
            call cuda_mg_interp_upload_c(mg_ri_coarse_phi, &
                 int(mg_ri_total_coarse,c_int))
            call cuda_mg_interp_execute_c(int(active(ilevel)%ngrid,c_int), &
                 int(ngridmax,c_int), int(ncoarse,c_int))
            call make_virtual_fine_dp_gpu(ilevel)
         else
#endif
#ifdef HYDRO_CUDA
            if(use_mg_gpu) then
               call cuda_mg_download_phi_c(phi, ncell_tot_c)
            end if
#endif
            call interpolate_and_correct_fine(ilevel)
            call make_virtual_fine_dp(phi(1),ilevel)
#ifdef HYDRO_CUDA
            if(use_mg_gpu) then
               call cuda_mg_upload_phi_c(phi, ncell_tot_c)
            end if
         end if
#endif
      end if

      ! Post-smoothing
      do i=1,ngs_fine
#ifdef HYDRO_CUDA
         if(use_mg_gpu) then
            safe_int = 0
            if(safe_mode(ilevel)) safe_int = 1
            call cuda_mg_gauss_seidel_c(int(active(ilevel)%ngrid,c_int), &
                 int(ngridmax,c_int), int(ncoarse,c_int), dx2_mg, 0, safe_int)
            call cuda_mg_gauss_seidel_c(int(active(ilevel)%ngrid,c_int), &
                 int(ngridmax,c_int), int(ncoarse,c_int), dx2_mg, 1, safe_int)
            call make_virtual_fine_dp_gpu(ilevel)
         else
#endif
            call gauss_seidel_mg_fine(ilevel,.true. )  ! Red step
            call gauss_seidel_mg_fine(ilevel,.false.)  ! Black step
            call make_virtual_fine_dp(phi(1),ilevel)   ! Communicate phi
#ifdef HYDRO_CUDA
         end if
#endif
      end do

      ! Update fine residual (fused with norm computation)
#ifdef HYDRO_CUDA
      if(use_mg_gpu) then
         gpu_norm2 = 0.0d0
         call cuda_mg_residual_c(int(active(ilevel)%ngrid,c_int), &
              int(ngridmax,c_int), int(ncoarse,c_int), &
              oneoverdx2_mg, dble(twondim), dx2_norm_mg, &
              gpu_norm2, 1)
         res_norm2 = gpu_norm2
      else
#endif
         call cmp_residual_mg_fine(ilevel, res_norm2)
#ifdef HYDRO_CUDA
      end if
#endif
#ifndef WITHOUTMPI
      call MPI_ALLREDUCE(res_norm2,res_norm2_tot,1, &
              & MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
      res_norm2=res_norm2_tot
#endif

      last_err = err
      err = sqrt(res_norm2/(i_res_norm2+1d-20*rho_tot**2))

      ! Verbosity
      if(verbose .or. use_mg_gpu) then
         if(myid==1) print '(A,I5,A,1pE10.3,A,I3,A,L1)', &
              '   ==> Step=',iter,' Error=',err, &
              ' level=',ilevel,' gpu=',use_mg_gpu
      end if

      ! Converged?
      if(err<epsilon .or. iter>=MAXITER) exit

      ! Not converged, check error and possibly enable safe mode for the level
      if(err > last_err*SAFE_FACTOR .and. (.not. safe_mode(ilevel))) then
         if(verbose)print *,'CAUTION: Switching to safe MG mode for level ',ilevel
         safe_mode(ilevel) = .true.
      end if

   end do main_iteration_loop

   ! Download final phi from GPU to CPU before cleanup
#ifdef HYDRO_CUDA
   if(use_mg_gpu) then
      call cuda_mg_download_phi_c(phi, ncell_tot_c)
   end if
#endif

   ! Cleanup GPU MG state
#ifdef HYDRO_CUDA
   if(use_ri_gpu) then
      call cuda_mg_ri_free_c()
   end if
   if(allocated(mg_ri_flat_offset))  deallocate(mg_ri_flat_offset)
   if(allocated(mg_ri_coarse_rhs))   deallocate(mg_ri_coarse_rhs)
   if(allocated(mg_ri_coarse_phi))   deallocate(mg_ri_coarse_phi)
   mg_ri_total_coarse = 0
   if(use_mg_gpu) then
      call cuda_mg_halo_free_c()
      call cuda_mg_free_c()
   end if
   if(allocated(mg_halo_emit_cells)) deallocate(mg_halo_emit_cells)
   if(allocated(mg_halo_recv_cells)) deallocate(mg_halo_recv_cells)
   if(allocated(mg_halo_emit_buf))   deallocate(mg_halo_emit_buf)
   if(allocated(mg_halo_recv_buf))   deallocate(mg_halo_recv_buf)
   mg_halo_n_emit = 0
   mg_halo_n_recv = 0
#endif

   ! Free pre-computed neighbor grids
   if(allocated(nbor_grid_fine)) deallocate(nbor_grid_fine)

   ! Free pre-computed coarse neighbor cache
   if(ilevel>1 .and. levelmin_mg < ilevel) then
      call cleanup_nbor_grid_coarse(levelmin_mg, ilevel-1)
   end if

   if(myid==1) print '(A,I5,A,I5,A,1pE10.3)','   ==> Level=',ilevel, ' Step=', &
            iter,' Error=',err
   if(myid==1 .and. iter==MAXITER) print *,'WARN: Fine multigrid &
      &Poisson failed to converge...'

   ! ---------------------------------------------------------------------
   ! Cleanup MG levels after solve complete
   ! ---------------------------------------------------------------------
   do ifine=1,ilevel-1
      call cleanup_mg_level(ifine)
   end do

end subroutine multigrid_fine


! ########################################################################
! ########################################################################
! ########################################################################
! ########################################################################

! ------------------------------------------------------------------------
! Build flat halo cell index arrays for GPU MG phi exchange
! Enumerates emission and reception cell indices in the same order
! as make_virtual_fine_dp packing.
! ------------------------------------------------------------------------
#ifdef HYDRO_CUDA
subroutine build_mg_halo_indices(ilevel)
   use amr_commons
   use poisson_commons
   use poisson_cuda_interface
   use iso_c_binding
   implicit none

   integer, intent(in) :: ilevel

   integer :: icpu, i, j, idx, iskip

   ! Count total emission and reception cells
   mg_halo_n_emit = 0
   mg_halo_n_recv = 0
   do icpu = 1, ncpu
      mg_halo_n_emit = mg_halo_n_emit + emission(icpu,ilevel)%ngrid * twotondim
      mg_halo_n_recv = mg_halo_n_recv + reception(icpu,ilevel)%ngrid * twotondim
   end do

   ! Allocate flat arrays
   if(allocated(mg_halo_emit_cells)) deallocate(mg_halo_emit_cells)
   if(allocated(mg_halo_recv_cells)) deallocate(mg_halo_recv_cells)
   if(allocated(mg_halo_emit_buf))   deallocate(mg_halo_emit_buf)
   if(allocated(mg_halo_recv_buf))   deallocate(mg_halo_recv_buf)

   allocate(mg_halo_emit_cells(1:max(mg_halo_n_emit,1)))
   allocate(mg_halo_recv_cells(1:max(mg_halo_n_recv,1)))
   allocate(mg_halo_emit_buf(1:max(mg_halo_n_emit,1)))
   allocate(mg_halo_recv_buf(1:max(mg_halo_n_recv,1)))

   ! Build emission cell indices
   ! Order: for each CPU, for each child cell j, for each grid i
   ! This matches the packing order in make_virtual_fine_dp
   idx = 0
   do icpu = 1, ncpu
      if(emission(icpu,ilevel)%ngrid > 0) then
         do j = 1, twotondim
            iskip = ncoarse + (j-1)*ngridmax
            do i = 1, emission(icpu,ilevel)%ngrid
               idx = idx + 1
               mg_halo_emit_cells(idx) = emission(icpu,ilevel)%igrid(i) + iskip
            end do
         end do
      end if
   end do

   ! Build reception cell indices (same order)
   idx = 0
   do icpu = 1, ncpu
      if(reception(icpu,ilevel)%ngrid > 0) then
         do j = 1, twotondim
            iskip = ncoarse + (j-1)*ngridmax
            do i = 1, reception(icpu,ilevel)%ngrid
               idx = idx + 1
               mg_halo_recv_cells(idx) = reception(icpu,ilevel)%igrid(i) + iskip
            end do
         end do
      end if
   end do

   if(myid==1) write(*,'(A,I3,A,I8,A,I8)') &
        ' MG halo indices: level=',ilevel, &
        ' n_emit=',mg_halo_n_emit,' n_recv=',mg_halo_n_recv

   ! Upload cell indices to GPU
   call cuda_mg_halo_setup_c( &
        mg_halo_emit_cells, int(mg_halo_n_emit, c_int), &
        mg_halo_recv_cells, int(mg_halo_n_recv, c_int))

end subroutine build_mg_halo_indices
#endif


! ------------------------------------------------------------------------
! GPU phi exchange for MG smoothing
! Full D2H → MPI exchange → full H2D
! ------------------------------------------------------------------------
#ifdef HYDRO_CUDA
subroutine make_virtual_fine_dp_gpu(ilevel)
   use amr_commons
   use poisson_commons
   use poisson_cuda_interface
   use iso_c_binding
   implicit none

   integer, intent(in) :: ilevel
   integer :: i

   if(mg_halo_n_emit == 0 .and. mg_halo_n_recv == 0) return

   ! Step 1: GPU gather — emission cell values from GPU phi → flat host buffer
   if(mg_halo_n_emit > 0) then
      call cuda_mg_halo_gather_c(mg_halo_emit_buf, int(mg_halo_n_emit, c_int))
      ! Step 2: Scatter flat buffer → host phi at emission positions
      do i = 1, mg_halo_n_emit
         phi(mg_halo_emit_cells(i)) = mg_halo_emit_buf(i)
      end do
   end if

   ! Step 3: MPI exchange (packs phi at emission cells, unpacks to reception cells)
   call make_virtual_fine_dp(phi(1), ilevel)

   ! Step 4: Gather host phi at reception positions → flat buffer
   if(mg_halo_n_recv > 0) then
      do i = 1, mg_halo_n_recv
         mg_halo_recv_buf(i) = phi(mg_halo_recv_cells(i))
      end do
      ! Step 5: GPU scatter — flat host buffer → GPU phi at reception positions
      call cuda_mg_halo_scatter_c(mg_halo_recv_buf, int(mg_halo_n_recv, c_int))
   end if

end subroutine make_virtual_fine_dp_gpu
#endif


! ########################################################################
! ########################################################################
! ########################################################################
! ########################################################################

! ------------------------------------------------------------------------
! Precompute GPU restrict/interp mapping arrays for V-cycle
! Builds restrict_target and interp_nbor_flat, then uploads to GPU.
! Also allocates flat host buffers for coarse data transfer.
! ------------------------------------------------------------------------
#ifdef HYDRO_CUDA
subroutine precompute_mg_gpu_restrict_interp(ilevel)
   use amr_commons
   use poisson_commons
   use poisson_cuda_interface
   use iso_c_binding
   implicit none

   integer, intent(in) :: ilevel

   integer :: icoarselevel, ngrid_fine
   integer :: igrid_f_mg, igrid_f_amr, icell_c_amr, ind_c_cell
   integer :: igrid_c_amr, igrid_c_mg, cpu_amr, icell_c_mg
   integer :: i, j, istart, nbatch, icpu, ngrid_c

   integer, allocatable :: h_restrict_target(:)
   integer, allocatable :: h_interp_nbor_flat(:,:)

   real(dp) :: a, b, c, d
   real(dp) :: bbb(8)
   integer :: ccc(8,8), ccc_gpu(8,8)

   ! Work arrays for get3cubefather
   integer :: ind_cell_father_loc(1:nvector)
   integer :: nbors_father_cells_loc(1:nvector, 1:threetondim)
   integer :: nbors_father_grids_loc(1:nvector, 1:twotondim)

   icoarselevel = ilevel - 1
   ngrid_fine = active(ilevel)%ngrid

   if(ngrid_fine == 0 .or. icoarselevel < 1) return

   ! Compute flat_offset for all CPUs
   if(allocated(mg_ri_flat_offset)) deallocate(mg_ri_flat_offset)
   allocate(mg_ri_flat_offset(1:ncpu))
   mg_ri_flat_offset(1) = 0
   do icpu = 2, ncpu
      mg_ri_flat_offset(icpu) = mg_ri_flat_offset(icpu-1) + &
           active_mg(icpu-1,icoarselevel)%ngrid * twotondim
   end do
   mg_ri_total_coarse = mg_ri_flat_offset(ncpu) + &
        active_mg(ncpu,icoarselevel)%ngrid * twotondim

   if(mg_ri_total_coarse == 0) return

   ! Allocate host flat buffers for coarse data
   if(allocated(mg_ri_coarse_rhs)) deallocate(mg_ri_coarse_rhs)
   if(allocated(mg_ri_coarse_phi)) deallocate(mg_ri_coarse_phi)
   allocate(mg_ri_coarse_rhs(1:mg_ri_total_coarse))
   allocate(mg_ri_coarse_phi(1:mg_ri_total_coarse))

   ! Set interpolation coefficients (from interpolate_and_correct_fine)
   a = 1.0D0/4.0D0**ndim
   b = 3.0D0*a
   c = 9.0D0*a
   d = 27.0D0*a
   bbb(:) = (/a, b, b, c, b, c, c, d/)

   ccc(:,1)=(/1 ,2 ,4 ,5 ,10,11,13,14/)
   ccc(:,2)=(/3 ,2 ,6 ,5 ,12,11,15,14/)
   ccc(:,3)=(/7 ,8 ,4 ,5 ,16,17,13,14/)
   ccc(:,4)=(/9 ,8 ,6 ,5 ,18,17,15,14/)
   ccc(:,5)=(/19,20,22,23,10,11,13,14/)
   ccc(:,6)=(/21,20,24,23,12,11,15,14/)
   ccc(:,7)=(/25,26,22,23,16,17,13,14/)
   ccc(:,8)=(/27,26,24,23,18,17,15,14/)

   ! Transpose ccc for GPU C row-major layout:
   ! Fortran ccc_gpu(f,a) stored at offset (a-1)*8+(f-1) matches C d_ccc[a-1][f-1]
   do i = 1, 8
      do j = 1, 8
         ccc_gpu(j, i) = ccc(i, j)
      end do
   end do

   ! === Compute restrict_target ===
   allocate(h_restrict_target(1:ngrid_fine))

   do igrid_f_mg = 1, ngrid_fine
      igrid_f_amr = active(ilevel)%igrid(igrid_f_mg)
      icell_c_amr = father(igrid_f_amr)
      ind_c_cell  = (icell_c_amr - ncoarse - 1) / ngridmax + 1
      igrid_c_amr = icell_c_amr - ncoarse - (ind_c_cell - 1) * ngridmax
      cpu_amr     = cpu_map(father(igrid_c_amr))
      igrid_c_mg  = lookup_mg(igrid_c_amr)

      if(igrid_c_mg <= 0) then
         h_restrict_target(igrid_f_mg) = -1
      else
         ngrid_c = active_mg(cpu_amr,icoarselevel)%ngrid
         icell_c_mg = (ind_c_cell - 1) * ngrid_c + igrid_c_mg
         ! Check coarse mask
         if(active_mg(cpu_amr,icoarselevel)%u(icell_c_mg,4) <= 0d0) then
            h_restrict_target(igrid_f_mg) = -1
         else
            h_restrict_target(igrid_f_mg) = mg_ri_flat_offset(cpu_amr) + &
                 icell_c_mg - 1  ! 0-based for C
         end if
      end if
   end do

   ! === Compute interp_nbor_flat ===
   allocate(h_interp_nbor_flat(1:27, 1:ngrid_fine))

   do istart = 1, ngrid_fine, nvector
      nbatch = min(nvector, ngrid_fine - istart + 1)

      ! Get father cells
      do i = 1, nbatch
         ind_cell_father_loc(i) = father(active(ilevel)%igrid(istart+i-1))
      end do

      ! Get 3x3x3 neighbor cells
      call get3cubefather(ind_cell_father_loc, nbors_father_cells_loc, &
           nbors_father_grids_loc, nbatch, ilevel)

      ! Convert to flat indices
      do j = 1, threetondim  ! 27
         do i = 1, nbatch
            igrid_f_mg  = istart + i - 1
            icell_c_amr = nbors_father_cells_loc(i, j)
            ind_c_cell  = (icell_c_amr - ncoarse - 1) / ngridmax + 1
            igrid_c_amr = icell_c_amr - ncoarse - (ind_c_cell - 1) * ngridmax
            cpu_amr     = cpu_map(father(igrid_c_amr))
            igrid_c_mg  = lookup_mg(igrid_c_amr)

            if(igrid_c_mg <= 0) then
               h_interp_nbor_flat(j, igrid_f_mg) = -1
            else
               ngrid_c = active_mg(cpu_amr,icoarselevel)%ngrid
               icell_c_mg = (ind_c_cell - 1) * ngrid_c + igrid_c_mg
               h_interp_nbor_flat(j, igrid_f_mg) = &
                    mg_ri_flat_offset(cpu_amr) + icell_c_mg - 1  ! 0-based
            end if
         end do
      end do
   end do

   ! Upload to GPU
   call cuda_mg_ri_setup_c(h_restrict_target, h_interp_nbor_flat, &
        int(ngrid_fine, c_int), int(mg_ri_total_coarse, c_int), &
        bbb, ccc_gpu)

   deallocate(h_restrict_target, h_interp_nbor_flat)

   if(myid==1) write(*,'(A,I3,A,I10,A,I10)') &
        ' MG RI precompute: level=',ilevel, &
        ' ngrid_fine=',ngrid_fine,' total_coarse=',mg_ri_total_coarse

end subroutine precompute_mg_gpu_restrict_interp
#endif


! ------------------------------------------------------------------------
! Scatter coarse RHS from flat array into active_mg structure
! Called after GPU restrict download
! ------------------------------------------------------------------------
#ifdef HYDRO_CUDA
subroutine scatter_coarse_rhs_from_flat(ilevel)
   use amr_commons
   use poisson_commons
   implicit none

   integer, intent(in) :: ilevel
   integer :: icoarselevel, icpu, ngrid_c, ncells_c, offset

   icoarselevel = ilevel - 1

   ! Scatter flat data into active_mg%u(:,2)
   ! The flat buffer replaces (not accumulates) since it was zeroed on GPU
   do icpu = 1, ncpu
      ngrid_c = active_mg(icpu,icoarselevel)%ngrid
      if(ngrid_c == 0) cycle
      ncells_c = ngrid_c * twotondim
      offset = mg_ri_flat_offset(icpu)
      active_mg(icpu,icoarselevel)%u(1:ncells_c,2) = &
           mg_ri_coarse_rhs(offset+1:offset+ncells_c)
   end do

end subroutine scatter_coarse_rhs_from_flat
#endif


! ------------------------------------------------------------------------
! Gather coarse phi from active_mg structure into flat array
! Called before GPU interp upload
! ------------------------------------------------------------------------
#ifdef HYDRO_CUDA
subroutine gather_coarse_phi_to_flat(ilevel)
   use amr_commons
   use poisson_commons
   implicit none

   integer, intent(in) :: ilevel
   integer :: icoarselevel, icpu, ngrid_c, ncells_c, offset

   icoarselevel = ilevel - 1

   ! Gather active_mg%u(:,1) into flat array
   do icpu = 1, ncpu
      ngrid_c = active_mg(icpu,icoarselevel)%ngrid
      if(ngrid_c == 0) cycle
      ncells_c = ngrid_c * twotondim
      offset = mg_ri_flat_offset(icpu)
      mg_ri_coarse_phi(offset+1:offset+ncells_c) = &
           active_mg(icpu,icoarselevel)%u(1:ncells_c,1)
   end do

end subroutine gather_coarse_phi_to_flat
#endif


! ########################################################################
! ########################################################################
! ########################################################################
! ########################################################################

! ------------------------------------------------------------------------
! Recursive multigrid routine for coarse MG levels
! ------------------------------------------------------------------------
recursive subroutine recursive_multigrid_coarse(ifinelevel, safe)
   use amr_commons
   use poisson_commons
   implicit none
#ifndef WITHOUTMPI
   include "mpif.h"
#endif

   integer, intent(in) :: ifinelevel
   logical, intent(in) :: safe

   real(dp) :: debug_norm2, debug_norm2_tot
   integer :: i, icpu, info, icycle, ncycle

   if(ifinelevel<=levelmin_mg) then
      ! Solve 'directly'
      do i=1,2*ngs_coarse
         call gauss_seidel_mg_coarse(ifinelevel,safe,.true. )  ! Red step
         call gauss_seidel_mg_coarse(ifinelevel,safe,.false.)  ! Black step
         call make_virtual_mg_dp(1,ifinelevel)  ! Communicate solution
      end do
      return
   end if

   if(safe) then
      ncycle=ncycles_coarse_safe
   else
      ncycle=1
   endif

   do icycle=1,ncycle

      ! Pre-smoothing
      do i=1,ngs_coarse
         call gauss_seidel_mg_coarse(ifinelevel,safe,.true. )  ! Red step
         call gauss_seidel_mg_coarse(ifinelevel,safe,.false.)  ! Black step
         call make_virtual_mg_dp(1,ifinelevel)  ! Communicate solution
      end do

      ! Compute residual and restrict into upper level RHS
      call cmp_residual_mg_coarse(ifinelevel)

      ! First clear the rhs in coarser reception comms
      do icpu=1,ncpu
         if(active_mg(icpu,ifinelevel-1)%ngrid==0) cycle
         active_mg(icpu,ifinelevel-1)%u(:,2)=0.0d0
      end do

      ! Restrict and do reverse-comm
      call restrict_residual_coarse_reverse(ifinelevel)
      call make_reverse_mg_dp(2,ifinelevel-1) ! communicate rhs

      ! Reset correction from upper level before solve
      do icpu=1,ncpu
         if(active_mg(icpu,ifinelevel-1)%ngrid==0) cycle
         active_mg(icpu,ifinelevel-1)%u(:,1)=0.0d0
      end do

      ! Multigrid-solve the upper level
      call recursive_multigrid_coarse(ifinelevel-1, safe)

      ! Interpolate coarse solution and correct back into fine solution
      call interpolate_and_correct_coarse(ifinelevel)
      call make_virtual_mg_dp(1,ifinelevel)  ! Communicate solution

      ! Post-smoothing
      do i=1,ngs_coarse
         call gauss_seidel_mg_coarse(ifinelevel,safe,.true. )  ! Red step
         call gauss_seidel_mg_coarse(ifinelevel,safe,.false.)  ! Black step
         call make_virtual_mg_dp(1,ifinelevel)  ! Communicate solution
      end do

   end do

end subroutine recursive_multigrid_coarse

! ########################################################################
! ########################################################################
! ########################################################################
! ########################################################################

! ------------------------------------------------------------------------
! Multigrid communicator building
! ------------------------------------------------------------------------
subroutine build_parent_comms_mg(active_f_comm, ifinelevel)
   use amr_commons
   use poisson_commons
   use ksection
   implicit none

#ifndef WITHOUTMPI
   include "mpif.h"
   integer, dimension (MPI_STATUS_SIZE, ncpu) :: statuses
   integer :: ntotal_ksec, nrecv_ksec, idx_ksec
   real(dp), allocatable :: sbuf_ksec(:,:), rbuf_ksec(:,:)
   integer, allocatable :: dcpu_ksec(:)
#endif

   integer, intent(in) :: ifinelevel
   type(communicator), intent(in) :: active_f_comm

   integer :: icoarselevel
   integer :: ngrids, cur_grid, cur_cpu, cur_cell, newgrids
   integer :: i, nbatch, ind, icpu, istart, info

   integer :: nact_tot, nreq_tot, nreq_tot2
   integer, dimension(1:ncpu) :: nreq, nreq2

   ! Per-thread work arrays for OpenMP (Stages 1 & 4)
   integer, dimension(:,:), target, allocatable :: P_icf_bp
   integer, dimension(:,:,:), target, allocatable :: P_nfg_bp, P_nfc_bp
   integer, dimension(:), pointer :: ind_cell_father
   integer, dimension(:,:), pointer :: nbors_father_grids
   integer, dimension(:,:), pointer :: nbors_father_cells
   common /omp_build_parent_comms/ ind_cell_father, nbors_father_cells, nbors_father_grids
!$omp threadprivate(/omp_build_parent_comms/)
   integer :: mythread, nthreads
   common /openmpthreads/ mythread, nthreads
!$omp threadprivate(/openmpthreads/)

   type(communicator), dimension(1:ncpu) :: comm_send, comm_receive
   type(communicator), dimension(1:ncpu) :: comm_send2, comm_receive2

   integer, dimension(1:ncpu) :: indx, recvbuf, recvbuf2
   integer, dimension(1:ncpu) :: reqsend, reqrecv
   integer :: countrecv, countsend
   integer :: tag = 777


   icoarselevel=ifinelevel-1

   nact_tot=0
   nreq_tot=0; nreq=0
   indx=0; recvbuf=0

   ! Setup per-thread work arrays for OpenMP
!$omp parallel
   mythread = omp_get_thread_num()
   nthreads = omp_get_num_threads()
!$omp end parallel
   allocate(P_icf_bp(1:nvector, 0:nthreads-1))
   allocate(P_nfg_bp(1:nvector, 1:twotondim, 0:nthreads-1))
   allocate(P_nfc_bp(1:nvector, 1:threetondim, 0:nthreads-1))
!$omp parallel
   ind_cell_father => P_icf_bp(:, mythread)
   nbors_father_grids => P_nfg_bp(:, :, mythread)
   nbors_father_cells => P_nfc_bp(:, :, mythread)
!$omp end parallel

   ! ---------------------------------------------------------------------
   ! STAGE 1 : Coarse grid MG activation for local grids (1st pass)
   ! ---------------------------------------------------------------------

   ! Loop over the AMR active communicator first
   ngrids = active_f_comm%ngrid
!$omp parallel do private(istart,nbatch,i,ind,cur_grid,cur_cpu) schedule(dynamic,128)
   do istart=1,ngrids,nvector
      nbatch=min(nvector,ngrids-istart+1)
      ! Gather grid indices and retrieve parent cells
      do i=1,nbatch
         ind_cell_father(i)=father( active_f_comm%igrid(istart+i-1) )
      end do

      ! Compute neighbouring father cells and grids
      call get3cubefather(ind_cell_father,nbors_father_cells, &
         & nbors_father_grids,nbatch,ifinelevel)

      ! Now process the twotondim father grids
      do ind=1,twotondim
         do i=1,nbatch
            cur_grid = nbors_father_grids(i,ind)
            if(lookup_mg(cur_grid)>0) cycle ! Grid already active (pre-check)

            cur_cpu=cpu_map(father(cur_grid))
            if(cur_cpu==0) cycle

            !$omp critical(stage1_update)
            if(lookup_mg(cur_grid)<=0) then  ! Definitive check under lock
               if(cur_cpu==myid) then
                  ! Stack grid for local activation
                  nact_tot=nact_tot+1
                  flag2(nact_tot)=cur_grid
                  lookup_mg(cur_grid)=nact_tot
               else
                  ! Stack grid for remote activation
                  nreq_tot=nreq_tot+1
                  nreq(cur_cpu)=nreq(cur_cpu)+1
                  flag2(ngridmax+nreq_tot)=cur_grid
                  lookup_mg(cur_grid)=abs(lookup_mg(cur_grid))
               end if
            end if
            !$omp end critical(stage1_update)
         end do
      end do
   end do


   ! ---------------------------------------------------------------------
   ! STAGE 2 : Coarse grid MG activation request
   ! ---------------------------------------------------------------------

#ifndef WITHOUTMPI
   ! Share number of requests and replies
   if(ordering=='ksection') then
      ntotal_ksec = 0
      do icpu = 1, ncpu
         if(nreq(icpu) > 0) ntotal_ksec = ntotal_ksec + 1
      end do
      allocate(sbuf_ksec(1:2, 1:max(ntotal_ksec,1)))
      allocate(dcpu_ksec(1:max(ntotal_ksec,1)))
      idx_ksec = 0
      do icpu = 1, ncpu
         if(nreq(icpu) > 0) then
            idx_ksec = idx_ksec + 1
            dcpu_ksec(idx_ksec) = icpu
            sbuf_ksec(1, idx_ksec) = dble(myid)
            sbuf_ksec(2, idx_ksec) = dble(nreq(icpu))
         end if
      end do
      call ksection_exchange_dp(sbuf_ksec, ntotal_ksec, dcpu_ksec, 2, rbuf_ksec, nrecv_ksec)
      recvbuf = 0
      do idx_ksec = 1, nrecv_ksec
         recvbuf(nint(rbuf_ksec(1, idx_ksec))) = nint(rbuf_ksec(2, idx_ksec))
      end do
      deallocate(sbuf_ksec, dcpu_ksec, rbuf_ksec)
   else
      call MPI_ALLTOALL(nreq, 1, MPI_INTEGER, recvbuf, 1, MPI_INTEGER, &
         & MPI_COMM_WORLD, info)
   end if

   ! Allocate inbound comms
   do icpu=1,ncpu
      comm_receive(icpu)%ngrid=recvbuf(icpu)
      if(recvbuf(icpu)>0) allocate(comm_receive(icpu)%igrid(1:recvbuf(icpu)))
   end do

   ! Receive to-be-activated grids
   countrecv=0; reqrecv=0
   do icpu=1,ncpu
      if(icpu==myid) cycle
      ngrids=comm_receive(icpu)%ngrid
      if(ngrids>0) then
         countrecv=countrecv+1
         call MPI_IRECV(comm_receive(icpu)%igrid, ngrids, MPI_INTEGER, &
            & icpu-1, tag, MPI_COMM_WORLD, reqrecv(countrecv), info)
      end if
   end do

   ! Allocate and then fill outbound (activation request) communicators
   do icpu=1,ncpu
      comm_send(icpu)%ngrid=nreq(icpu)
      if(nreq(icpu)>0) allocate(comm_send(icpu)%igrid(1:nreq(icpu)))
   end do
   nreq=0
   do i=1,nreq_tot
      cur_grid=flag2(ngridmax+i) ! Local AMR index
      cur_cpu =cpu_map(father(cur_grid))
      nreq(cur_cpu)=nreq(cur_cpu)+1
      comm_send(cur_cpu)%igrid(nreq(cur_cpu))=lookup_mg(cur_grid) ! Remote
   end do

   ! Send to-be-activated grids
   countsend=0; reqsend=0
   do icpu=1,ncpu
      if(icpu==myid) cycle
      ngrids=comm_send(icpu)%ngrid
      if(ngrids>0) then
         countsend=countsend+1
         call MPI_ISEND(comm_send(icpu)%igrid, ngrids, MPI_INTEGER, &
            & icpu-1, tag, MPI_COMM_WORLD, reqsend(countsend), info)
      end if
   end do

   ! Wait for completion of receives
   call MPI_WAITALL(countrecv, reqrecv, statuses, info)

   ! Activate requested grids
   do icpu=1,ncpu
      if(icpu==myid) cycle
      ngrids=comm_receive(icpu)%ngrid
      if(ngrids>0) then
         do i=1,ngrids
            cur_grid=comm_receive(icpu)%igrid(i) ! Local AMR index
            if(lookup_mg(cur_grid)>0) cycle      ! Already active: cycle
            ! Activate grid
            nact_tot=nact_tot+1
            flag2(nact_tot)=cur_grid
            lookup_mg(cur_grid)=nact_tot
         end do
      end if
   end do

   ! Wait for completion of sends
   call MPI_WAITALL(countsend, reqsend, statuses, info)
#endif

   ! ---------------------------------------------------------------------
   ! STAGE 3 : Coarse grid MG active comm gathering
   ! ---------------------------------------------------------------------

   active_mg(myid,icoarselevel)%ngrid=nact_tot
   if(nact_tot>0) then
      allocate( active_mg(myid,icoarselevel)%igrid(1:nact_tot) )
      allocate( active_mg(myid,icoarselevel)%u(1:nact_tot*twotondim,1:4) )
      allocate( active_mg(myid,icoarselevel)%f(1:nact_tot*twotondim,1:1) )
      active_mg(myid,icoarselevel)%igrid=0
      active_mg(myid,icoarselevel)%u=0.0d0
      active_mg(myid,icoarselevel)%f=0
   end if
   do i=1,nact_tot
      active_mg(myid,icoarselevel)%igrid(i)=flag2(i)
   end do

   ! ---------------------------------------------------------------------
   ! STAGE 4 : Screen active grid neighbors for new reception grids
   ! ---------------------------------------------------------------------
   ngrids = active_mg(myid,icoarselevel)%ngrid
   nreq2 = 0
   nreq_tot2 = 0
!$omp parallel do private(istart,nbatch,i,ind,cur_cell,cur_cpu,cur_grid) schedule(dynamic,128)
   do istart=1,ngrids,nvector
      nbatch=min(nvector,ngrids-istart+1)
      ! Gather grid indices and retrieve parent cells
      do i=1,nbatch
         ind_cell_father(i)=father( active_mg(myid,icoarselevel)%igrid(istart+i-1) )
      end do

      ! Compute neighbouring father cells
      call get3cubefather(ind_cell_father,nbors_father_cells,nbors_father_grids,nbatch,icoarselevel)

      ! Now process the father grids
      do ind=1,threetondim
         do i=1,nbatch
            cur_cell = nbors_father_cells(i,ind)
            cur_cpu  = cpu_map(cur_cell)
            if(cur_cpu==0) cycle
            cur_grid = son(cur_cell)
            if(cur_cpu/=myid) then
               ! Neighbor cell is not managed by current CPU
               if (cur_grid==0) cycle              ! No grid there
               if (lookup_mg(cur_grid)>0) cycle    ! Already selected (pre-check)
               !$omp critical(stage4_update)
               if(lookup_mg(cur_grid)<=0) then  ! Definitive check under lock
                  nreq_tot2=nreq_tot2+1
                  nreq2(cur_cpu)=nreq2(cur_cpu)+1
                  flag2(ngridmax+nreq_tot+nreq_tot2)=cur_grid
                  lookup_mg(cur_grid)=abs(lookup_mg(cur_grid))
               end if
               !$omp end critical(stage4_update)
            end if
         end do
      end do
   end do

   ! ---------------------------------------------------------------------
   ! STAGE 5 : Share new reception grid requests, build emission comms
   ! ---------------------------------------------------------------------

#ifndef WITHOUTMPI
   ! Share number of requests and replies
   recvbuf2=0
   if(ordering=='ksection') then
      ntotal_ksec = 0
      do icpu = 1, ncpu
         if(nreq2(icpu) > 0) ntotal_ksec = ntotal_ksec + 1
      end do
      allocate(sbuf_ksec(1:2, 1:max(ntotal_ksec,1)))
      allocate(dcpu_ksec(1:max(ntotal_ksec,1)))
      idx_ksec = 0
      do icpu = 1, ncpu
         if(nreq2(icpu) > 0) then
            idx_ksec = idx_ksec + 1
            dcpu_ksec(idx_ksec) = icpu
            sbuf_ksec(1, idx_ksec) = dble(myid)
            sbuf_ksec(2, idx_ksec) = dble(nreq2(icpu))
         end if
      end do
      call ksection_exchange_dp(sbuf_ksec, ntotal_ksec, dcpu_ksec, 2, rbuf_ksec, nrecv_ksec)
      do idx_ksec = 1, nrecv_ksec
         recvbuf2(nint(rbuf_ksec(1, idx_ksec))) = nint(rbuf_ksec(2, idx_ksec))
      end do
      deallocate(sbuf_ksec, dcpu_ksec, rbuf_ksec)
   else
      call MPI_ALLTOALL(nreq2, 1, MPI_INTEGER, recvbuf2, 1, MPI_INTEGER, &
         & MPI_COMM_WORLD, info)
   end if

   ! Allocate inbound comms
   do icpu=1,ncpu
      comm_receive2(icpu)%ngrid=recvbuf2(icpu)
      if(recvbuf2(icpu)>0) allocate(comm_receive2(icpu)%igrid(1:recvbuf2(icpu)))
   end do

   ! Receive potential reception grids
   countrecv=0; reqrecv=0
   do icpu=1,ncpu
      if(icpu==myid) cycle
      ngrids=comm_receive2(icpu)%ngrid
      if(ngrids>0) then
         countrecv=countrecv+1
         call MPI_IRECV(comm_receive2(icpu)%igrid, ngrids, MPI_INTEGER, &
            & icpu-1, tag, MPI_COMM_WORLD, reqrecv(countrecv), info)
      end if
   end do

   ! Allocate and then fill outbound (reception request) communicators
   do icpu=1,ncpu
      comm_send2(icpu)%ngrid=nreq2(icpu)
      if(nreq2(icpu)>0) allocate(comm_send2(icpu)%igrid(1:nreq2(icpu)))
   end do
   nreq2=0
   do i=1,nreq_tot2
      cur_grid=flag2(ngridmax+nreq_tot+i) ! Local AMR index
      cur_cpu =cpu_map(father(cur_grid))
      nreq2(cur_cpu)=nreq2(cur_cpu)+1
      comm_send2(cur_cpu)%igrid(nreq2(cur_cpu))=lookup_mg(cur_grid) ! Remote AMR index
      ! Restore negative lookup_mg
      lookup_mg(cur_grid)=-abs(lookup_mg(cur_grid))
   end do

   ! Send reception request grids
   countsend=0; reqsend=0
   do icpu=1,ncpu
      if(icpu==myid) cycle
      ngrids=comm_send2(icpu)%ngrid
      if(ngrids>0) then
         countsend=countsend+1
         call MPI_ISEND(comm_send2(icpu)%igrid,ngrids,MPI_INTEGER,icpu-1, &
            & tag, MPI_COMM_WORLD, reqsend(countsend), info)
      end if
   end do

   ! Wait for completion of receives
   call MPI_WAITALL(countrecv, reqrecv, statuses, info)

   ! Compute local MG indices of inbound grids, alloc & fill emission comms
   do icpu=1,ncpu
      if(icpu==myid) cycle
      newgrids=0
      do i=1,recvbuf2(icpu)
         ! MAP AMR -> MG INDICES IN PLACE
         comm_receive2(icpu)%igrid(i)=lookup_mg(comm_receive2(icpu)%igrid(i))
         if(comm_receive2(icpu)%igrid(i)>0) newgrids=newgrids+1
      end do
      ! Allocate emission communicators
      ngrids=recvbuf(icpu)+newgrids
      emission_mg(icpu,icoarselevel)%ngrid=ngrids
      if(ngrids>0) then
         allocate(emission_mg(icpu,icoarselevel)%igrid(1:ngrids))
         allocate(emission_mg(icpu,icoarselevel)%u(1:ngrids*twotondim,1:4) )
         allocate(emission_mg(icpu,icoarselevel)%f(1:ngrids*twotondim,1:1))
         emission_mg(icpu,icoarselevel)%igrid=0
         emission_mg(icpu,icoarselevel)%u=0.0d0
         emission_mg(icpu,icoarselevel)%f=0
      end if
      ! First part: activation request emission grids
      do i=1,recvbuf(icpu)
         emission_mg(icpu,icoarselevel)%igrid(i)=lookup_mg(comm_receive(icpu)%igrid(i))
      end do
      ! Second part: new emission grids
      cur_grid=recvbuf(icpu)
      do i=1,recvbuf2(icpu)
         if(comm_receive2(icpu)%igrid(i)>0) then
            cur_grid=cur_grid+1
            emission_mg(icpu,icoarselevel)%igrid(cur_grid)=comm_receive2(icpu)%igrid(i)
         end if
      end do
   end do

   ! Wait for completion of sends
   call MPI_WAITALL(countsend, reqsend, statuses, info)


   ! ---------------------------------------------------------------------
   ! STAGE 6 : Reply with local MG grid status and build reception comms
   ! ---------------------------------------------------------------------
   ! Receive MG mappings from other CPUs back into comm_send2
   countrecv=0; reqrecv=0
   do icpu=1,ncpu
      if(icpu==myid) cycle
      ngrids=nreq2(icpu)
      if(ngrids>0) then
         countrecv=countrecv+1
         call MPI_IRECV(comm_send2(icpu)%igrid,ngrids,MPI_INTEGER,icpu-1, &
            & tag, MPI_COMM_WORLD, reqrecv(countrecv), info)
      end if
   end do

   ! Send local MG mappings to other CPUs from comm_receive
   countsend=0; reqsend=0
   do icpu=1,ncpu
      if(icpu==myid) cycle
      ngrids=recvbuf2(icpu)
      if(ngrids>0) then
         countsend=countsend+1
         call MPI_ISEND(comm_receive2(icpu)%igrid,ngrids,MPI_INTEGER, &
            & icpu-1, tag, MPI_COMM_WORLD, reqsend(countsend), info)
      end if
   end do

   ! Wait for full completion of receives
   call MPI_WAITALL(countrecv, reqrecv, statuses, info)

   ! Count remotely active MG grids, and allocate and fill reception comms
   do icpu=1,ncpu
      if(icpu==myid) cycle
      ! Count requested grids which are MG-active remotely
      newgrids=0
      do i=1,nreq2(icpu)
         if(comm_send2(icpu)%igrid(i)>0) newgrids=newgrids+1
      end do
      ! Allocate and fill reception communicators on the fly
      ngrids=nreq(icpu)+newgrids
      active_mg(icpu,icoarselevel)%ngrid=ngrids
      if(ngrids>0) then
         allocate(active_mg(icpu,icoarselevel)%igrid(1:ngrids))
         allocate(active_mg(icpu,icoarselevel)%u(1:ngrids*twotondim,1:4))
         allocate(active_mg(icpu,icoarselevel)%f(1:ngrids*twotondim,1:1))
         active_mg(icpu,icoarselevel)%igrid=0
         active_mg(icpu,icoarselevel)%u=0.0d0
         active_mg(icpu,icoarselevel)%f=0
      end if
   end do
   
   nreq=0
   do i=1,nreq_tot
      cur_grid=flag2(ngridmax+i)
      cur_cpu =cpu_map(father(cur_grid))
      nreq(cur_cpu)=nreq(cur_cpu)+1
      ! Add to reception comm
      active_mg(cur_cpu,icoarselevel)%igrid(nreq(cur_cpu))=cur_grid
      ! Backup lookup_mg into flag2
      flag2(cur_grid)=lookup_mg(cur_grid)
      ! Update lookup_mg
      lookup_mg(cur_grid)=nreq(cur_cpu)
   end do

   nreq2=0; indx=nreq
   do i=1,nreq_tot2
      cur_grid=flag2(ngridmax+nreq_tot+i)
      cur_cpu =cpu_map(father(cur_grid))
      nreq2(cur_cpu)=nreq2(cur_cpu)+1
      if(comm_send2(cur_cpu)%igrid(nreq2(cur_cpu))>0) then
         indx(cur_cpu)=indx(cur_cpu)+1
         ! Add to reception comm
         active_mg(cur_cpu,icoarselevel)%igrid(indx(cur_cpu))=cur_grid
         ! Backup lookup_mg
         flag2(cur_grid)=-lookup_mg(cur_grid)
         ! Update lookup_mg
         lookup_mg(cur_grid)=indx(cur_cpu)
      end if
   end do

   ! Wait for full completion of sends
   call MPI_WAITALL(countsend, reqsend, statuses, info)


   ! Cleanup
   do icpu=1,ncpu
      if(comm_send (icpu)%ngrid>0) deallocate(comm_send (icpu)%igrid)
      if(comm_send2(icpu)%ngrid>0) deallocate(comm_send2(icpu)%igrid)
      if(comm_receive (icpu)%ngrid>0) deallocate(comm_receive (icpu)%igrid)
      if(comm_receive2(icpu)%ngrid>0) deallocate(comm_receive2(icpu)%igrid)
   end do
#endif

   deallocate(P_icf_bp, P_nfg_bp, P_nfc_bp)

end subroutine build_parent_comms_mg


! ########################################################################
! ########################################################################
! ########################################################################
! ########################################################################

! ------------------------------------------------------------------------
! Multigrid level cleanup
! ------------------------------------------------------------------------
subroutine cleanup_mg_level(ilevel)
   use amr_commons
   use pm_commons
   use poisson_commons
   implicit none

   integer, intent(in) :: ilevel

   integer :: igrid, icpu, cur_grid, cur_cpu

   ! ---------------------------------------------------------------------
   ! Cleanup lookup table
   ! ---------------------------------------------------------------------
   do icpu=1,ncpu
      do igrid=1,active_mg(icpu,ilevel)%ngrid
         cur_grid=active_mg(icpu,ilevel)%igrid(igrid)
         cur_cpu=cpu_map(father(cur_grid))
         if(cur_cpu==myid) then
            lookup_mg(cur_grid)=0
         else
            lookup_mg(cur_grid)=-mod(flag2(cur_grid),ngridmax)
         end if
      end do
   end do

   ! ---------------------------------------------------------------------
   ! Deallocate communicators
   ! ---------------------------------------------------------------------
   do icpu=1,ncpu
      if(active_mg(icpu,ilevel)%ngrid>0)then
         deallocate(active_mg(icpu,ilevel)%igrid)
         deallocate(active_mg(icpu,ilevel)%u)
         deallocate(active_mg(icpu,ilevel)%f)
      endif
      active_mg(icpu,ilevel)%ngrid=0
      if(emission_mg(icpu,ilevel)%ngrid>0)then
         deallocate(emission_mg(icpu,ilevel)%igrid)
         deallocate(emission_mg(icpu,ilevel)%u)
         deallocate(emission_mg(icpu,ilevel)%f)
      endif
      emission_mg(icpu,ilevel)%ngrid=0
   end do

end subroutine cleanup_mg_level

! ########################################################################
! ########################################################################
! ########################################################################
! ########################################################################

! ------------------------------------------------------------------------
! Initialize mask at fine level into f(:,3)
! ------------------------------------------------------------------------
subroutine make_fine_mask(ilevel)

   use amr_commons
   use pm_commons
   use poisson_commons
   implicit none
   integer, intent(in) :: ilevel

   integer  :: ngrid
   integer  :: ind, igrid_mg, icpu, ibound
   integer  :: igrid_amr, icell_amr, iskip_amr

   ngrid=active(ilevel)%ngrid
   do ind=1,twotondim
      iskip_amr = ncoarse+(ind-1)*ngridmax
      do igrid_mg=1,ngrid
         igrid_amr = active(ilevel)%igrid(igrid_mg)
         icell_amr = iskip_amr + igrid_amr
         ! Init mask to 1.0 on active cells :
         f(icell_amr,3) = 1.0d0
      end do
   end do

   do icpu=1,ncpu
      ngrid=reception(icpu,ilevel)%ngrid
      do ind=1,twotondim
         iskip_amr = ncoarse+(ind-1)*ngridmax
         do igrid_mg=1,ngrid
            igrid_amr = reception(icpu,ilevel)%igrid(igrid_mg)
            icell_amr = iskip_amr + igrid_amr
            ! Init mask to 1.0 on virtual cells :
            f(icell_amr,3) = 1.0d0
         end do
      end do
   end do

   do ibound=1,nboundary
      ngrid=boundary(ibound,ilevel)%ngrid
      do ind=1,twotondim
         iskip_amr=ncoarse+(ind-1)*ngridmax
         do igrid_mg=1,ngrid
            igrid_amr = boundary(ibound,ilevel)%igrid(igrid_mg)
            icell_amr = iskip_amr + igrid_amr
            ! Init mask to -1.0 on boundary cells :
            f(icell_amr,3) = -1.0d0
         end do
      end do
   end do

end subroutine make_fine_mask

! ########################################################################
! ########################################################################
! ########################################################################
! ########################################################################

! ------------------------------------------------------------------------
! Preprocess the fine (AMR) level RHS to account for boundary conditions
!
!  _____#_____
! |     #     |      Cell I is INSIDE active domain (mask > 0)
! |  I  #  O  |      Cell O is OUTSIDE (mask <= 0 or nonexistent cell)
! |_____#_____|      # is the boundary
!       #
!
! phi(I) and phi(O) must BOTH be set at call time, if applicable
! phi(#) is computed from phi(I), phi(O) and the mask values
! If AMR cell O does not exist, phi(O) is computed by interpolation
!
! Sets BC-modified RHS    into f(:,2)
!
! ------------------------------------------------------------------------
subroutine make_fine_bc_rhs(ilevel,icount)

   use amr_commons
   use pm_commons
   use poisson_commons
   use morton_hash
   implicit none
   integer, intent(in) :: ilevel,icount

   integer, dimension(1:3,1:2,1:8) :: iii, jjj

   real(dp) :: dx, oneoverdx2

   integer  :: ngrid
   integer  :: ind, igrid_mg, idim, inbor
   integer  :: igrid_amr, icell_amr, iskip_amr
   integer  :: igshift, igrid_nbor_amr, icell_nbor_amr
   integer  :: ifathercell_nbor_amr

   ! Thread-private variables for OpenMP
   real(dp) :: phi_b, nb_mask, nb_phi, w
   real(dp), dimension(1:nvector,1:twotondim) :: phi_int
   integer,  dimension(1:nvector) :: ind_cell

   integer  :: nx_loc
   real(dp) :: scale, fourpi

   ! Set constants
   nx_loc = icoarse_max-icoarse_min+1
   scale  = boxlen/dble(nx_loc)
   fourpi = 4.D0*ACOS(-1.0D0)*scale
   if(cosmo) fourpi = 1.5D0*omega_m*aexp*scale

   dx  = 0.5d0**ilevel
   oneoverdx2 = 1.0d0/(dx*dx)

   iii(1,1,1:8)=(/1,0,1,0,1,0,1,0/); jjj(1,1,1:8)=(/2,1,4,3,6,5,8,7/)
   iii(1,2,1:8)=(/0,2,0,2,0,2,0,2/); jjj(1,2,1:8)=(/2,1,4,3,6,5,8,7/)
   iii(2,1,1:8)=(/3,3,0,0,3,3,0,0/); jjj(2,1,1:8)=(/3,4,1,2,7,8,5,6/)
   iii(2,2,1:8)=(/0,0,4,4,0,0,4,4/); jjj(2,2,1:8)=(/3,4,1,2,7,8,5,6/)
   iii(3,1,1:8)=(/5,5,5,5,0,0,0,0/); jjj(3,1,1:8)=(/5,6,7,8,1,2,3,4/)
   iii(3,2,1:8)=(/0,0,0,0,6,6,6,6/); jjj(3,2,1:8)=(/5,6,7,8,1,2,3,4/)

   ngrid=active(ilevel)%ngrid

   ! Loop over cells
   do ind=1,twotondim
      iskip_amr = ncoarse+(ind-1)*ngridmax

      ! Loop over active grids — OpenMP parallelized
      ! Each igrid_mg writes to unique icell_amr, so no race condition
!$omp parallel do private(igrid_mg,igrid_amr,icell_amr,idim,inbor, &
!$omp    igshift,igrid_nbor_amr,icell_nbor_amr,ifathercell_nbor_amr, &
!$omp    nb_mask,nb_phi,w,phi_b,phi_int,ind_cell) &
!$omp schedule(dynamic,1024)
      do igrid_mg=1,ngrid
         igrid_amr = active(ilevel)%igrid(igrid_mg)
         icell_amr = iskip_amr + igrid_amr

         ! Init BC-modified RHS to rho - rho_tot :
         f(icell_amr,2) = fourpi*(rho(icell_amr) - rho_tot)

         if(f(icell_amr,3)<=0.0) cycle ! Do not process masked cells

         ! Separate directions
         do idim=1,ndim
            ! Loop over the 2 neighbors
            do inbor=1,2
               ! Get neighbor grid shift
               igshift = iii(idim,inbor,ind)

               ! Get neighbor grid using precomputed array
               if(igshift==0) then
                  igrid_nbor_amr = igrid_amr
               else
                  igrid_nbor_amr = nbor_grid_fine(igshift, igrid_mg)
               end if

               if(igrid_nbor_amr==0) then
                  ! No neighbor (rare boundary case): interp. phi
                  nb_mask = -1.0d0
                  ! Only call morton_nbor_cell for this rare case
                  ifathercell_nbor_amr = morton_nbor_cell(igrid_amr,ilevel,igshift)
                  ind_cell(1)=ifathercell_nbor_amr
                  call interpol_phi(ind_cell,phi_int,1,ilevel,icount)
                  nb_phi = phi_int(1,jjj(idim,inbor,ind))
               else
                  ! Fetch neighbor cell id
                  icell_nbor_amr = igrid_nbor_amr + (ncoarse + (jjj(idim,inbor,ind)-1)*ngridmax)
                  ! Check neighbor cell mask
                  nb_mask = f(icell_nbor_amr,3)
                  if(nb_mask>0) cycle ! Neighbor cell is active too: cycle
                  nb_phi  = phi(icell_nbor_amr)
               end if
               ! phi(#) interpolated with mask:
               w = nb_mask/(nb_mask-f(icell_amr,3)) ! Linear parameter
               phi_b = ((1.0d0-w)*nb_phi + w*phi(icell_amr))

               ! Increment correction for current cell
               f(icell_amr,2) = f(icell_amr,2) - 2.0d0*oneoverdx2*phi_b
            end do
         end do
      end do
!$omp end parallel do
   end do

end subroutine make_fine_bc_rhs


! ########################################################################
! ########################################################################
! ########################################################################
! ########################################################################

! ------------------------------------------------------------------------
! MPI routines for MG communication for CPU boundaries,
! Those are the MG versions of the make_virtual_* AMR routines
! ------------------------------------------------------------------------

subroutine make_virtual_mg_dp(ivar,ilevel)
  use amr_commons
  use poisson_commons
  use ksection

  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
  integer,dimension(MPI_STATUS_SIZE,ncpu)::statuses
#endif
  integer::ilevel,ivar,icell
  integer::icpu,i,j,ncache,iskip,step
  integer::countsend,countrecv
  integer::info,tag=101
  integer,dimension(ncpu)::reqsend,reqrecv

#ifndef WITHOUTMPI
  if(ordering=='ksection') then
     call make_virtual_mg_dp_ksec(ivar,ilevel)
     return
  end if

  ! Receive all messages
  countrecv=0
  do icpu=1,ncpu
     if(icpu==myid)cycle
     ncache=active_mg(icpu,ilevel)%ngrid
     if(ncache>0) then
       countrecv=countrecv+1
       call MPI_IRECV(active_mg(icpu,ilevel)%u(1,ivar),ncache*twotondim, &
            & MPI_DOUBLE_PRECISION,icpu-1,tag,MPI_COMM_WORLD,reqrecv(countrecv),info)
     end if
  end do

  ! Gather emission array
  do icpu=1,ncpu
     if (emission_mg(icpu,ilevel)%ngrid>0) then
        do j=1,twotondim
           step=(j-1)*emission_mg(icpu,ilevel)%ngrid
           iskip=(j-1)*active_mg(myid,ilevel)%ngrid
           do i=1,emission_mg(icpu,ilevel)%ngrid
              icell=emission_mg(icpu,ilevel)%igrid(i)+iskip
              emission_mg(icpu,ilevel)%u(i+step,1)=active_mg(myid,ilevel)%u(icell,ivar)
           end do
        end do
     end if
  end do

  ! Send all messages
  countsend=0
  do icpu=1,ncpu
     ncache=emission_mg(icpu,ilevel)%ngrid
     if(ncache>0) then
       countsend=countsend+1
       call MPI_ISEND(emission_mg(icpu,ilevel)%u,ncache*twotondim, &
            & MPI_DOUBLE_PRECISION,icpu-1,tag,MPI_COMM_WORLD,reqsend(countsend),info)
     end if
  end do

  ! Wait for full completion of receives
  call MPI_WAITALL(countrecv,reqrecv,statuses,info)

  ! Wait for full completion of sends
  call MPI_WAITALL(countsend,reqsend,statuses,info)

#endif

111 format('   Entering make_virtual_mg for level ',I2)

end subroutine make_virtual_mg_dp

! ########################################################################
! ########################################################################

subroutine make_virtual_mg_int(ilevel)
  use amr_commons
  use poisson_commons
  use ksection
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
  integer,dimension(MPI_STATUS_SIZE,ncpu)::statuses
#endif
  integer::ilevel
  integer::icpu,i,j,ncache,iskip,step,icell
  integer::countsend,countrecv
  integer::info,tag=101
  integer,dimension(ncpu)::reqsend,reqrecv

#ifndef WITHOUTMPI
  if(ordering=='ksection') then
     call make_virtual_mg_int_ksec(ilevel)
     return
  end if

  ! Receive all messages
  countrecv=0
  do icpu=1,ncpu
     if(icpu==myid)cycle
     ncache=active_mg(icpu,ilevel)%ngrid
     if(ncache>0) then
       countrecv=countrecv+1
       call MPI_IRECV(active_mg(icpu,ilevel)%f(1,1),ncache*twotondim, &
            & MPI_INTEGER,icpu-1,tag,MPI_COMM_WORLD,reqrecv(countrecv),info)
     end if
  end do

  ! Gather emission array
  do icpu=1,ncpu
     if (emission_mg(icpu,ilevel)%ngrid>0) then
        do j=1,twotondim
           step=(j-1)*emission_mg(icpu,ilevel)%ngrid
           iskip=(j-1)*active_mg(myid,ilevel)%ngrid
           do i=1,emission_mg(icpu,ilevel)%ngrid
              icell=emission_mg(icpu,ilevel)%igrid(i)+iskip
              emission_mg(icpu,ilevel)%f(i+step,1)=active_mg(myid,ilevel)%f(icell,1)
           end do
        end do
     end if
  end do

  ! Send all messages
  countsend=0
  do icpu=1,ncpu
     ncache=emission_mg(icpu,ilevel)%ngrid
     if(ncache>0) then
       countsend=countsend+1
       call MPI_ISEND(emission_mg(icpu,ilevel)%f,ncache*twotondim, &
            & MPI_INTEGER,icpu-1,tag,MPI_COMM_WORLD,reqsend(countsend),info)
     end if
  end do

  ! Wait for full completion of receives
  call MPI_WAITALL(countrecv,reqrecv,statuses,info)

  ! Wait for full completion of sends
  call MPI_WAITALL(countsend,reqsend,statuses,info)

#endif

111 format('   Entering make_virtual_mg for level ',I2)

end subroutine make_virtual_mg_int

! ########################################################################
! ########################################################################

subroutine make_reverse_mg_dp(ivar,ilevel)
  use amr_commons
  use poisson_commons
  use ksection
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
  integer,dimension(MPI_STATUS_SIZE,ncpu)::statuses
#endif
  integer::ilevel,ivar,icell
  integer::icpu,i,j,ncache,iskip,step
  integer::countsend,countrecv
  integer::info,tag=101
  integer,dimension(ncpu)::reqsend,reqrecv

#ifndef WITHOUTMPI
  if(ordering=='ksection') then
     call make_reverse_mg_dp_ksec(ivar,ilevel)
     return
  end if

  ! Receive all messages
  countrecv=0
  do icpu=1,ncpu
     ncache=emission_mg(icpu,ilevel)%ngrid
     if(ncache>0) then
       countrecv=countrecv+1
       call MPI_IRECV(emission_mg(icpu,ilevel)%u,ncache*twotondim, &
            & MPI_DOUBLE_PRECISION,icpu-1,tag,MPI_COMM_WORLD,reqrecv(countrecv),info)
     end if
  end do

  ! Send all messages
  countsend=0
  do icpu=1,ncpu
     if(icpu==myid)cycle
     ncache=active_mg(icpu,ilevel)%ngrid
     if(ncache>0) then
       countsend=countsend+1
       call MPI_ISEND(active_mg(icpu,ilevel)%u(1,ivar),ncache*twotondim, &
            & MPI_DOUBLE_PRECISION,icpu-1,tag,MPI_COMM_WORLD,reqsend(countsend),info)
     end if
  end do

  ! Wait for full completion of receives
  call MPI_WAITALL(countrecv,reqrecv,statuses,info)

  ! Gather emission array
  do icpu=1,ncpu
     if (emission_mg(icpu,ilevel)%ngrid>0) then
        do j=1,twotondim
           step=(j-1)*emission_mg(icpu,ilevel)%ngrid
           iskip=(j-1)*active_mg(myid,ilevel)%ngrid
           do i=1,emission_mg(icpu,ilevel)%ngrid
              icell=emission_mg(icpu,ilevel)%igrid(i)+iskip
              active_mg(myid,ilevel)%u(icell,ivar)=active_mg(myid,ilevel)%u(icell,ivar)+ &
                   & emission_mg(icpu,ilevel)%u(i+step,1)
           end do
        end do
     end if
  end do

  ! Wait for full completion of sends
  call MPI_WAITALL(countsend,reqsend,statuses,info)

#endif

111 format('   Entering make_reverse_mg for level ',I2)

end subroutine make_reverse_mg_dp

! ########################################################################
! ########################################################################

subroutine make_reverse_mg_int(ilevel)
  use amr_commons
  use poisson_commons
  use ksection
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
  integer,dimension(MPI_STATUS_SIZE,ncpu)::statuses
#endif
  integer::ilevel,icell
  integer::icpu,i,j,ncache,iskip,step
  integer::countsend,countrecv
  integer::info,tag=101
  integer,dimension(ncpu)::reqsend,reqrecv

#ifndef WITHOUTMPI
  if(ordering=='ksection') then
     call make_reverse_mg_int_ksec(ilevel)
     return
  end if

  ! Receive all messages
  countrecv=0
  do icpu=1,ncpu
     ncache=emission_mg(icpu,ilevel)%ngrid
     if(ncache>0) then
       countrecv=countrecv+1
       call MPI_IRECV(emission_mg(icpu,ilevel)%f,ncache*twotondim, &
            & MPI_INTEGER,icpu-1,tag,MPI_COMM_WORLD,reqrecv(countrecv),info)
     end if
  end do

  ! Send all messages
  countsend=0
  do icpu=1,ncpu
     if(icpu==myid)cycle
     ncache=active_mg(icpu,ilevel)%ngrid
     if(ncache>0) then
       countsend=countsend+1
       call MPI_ISEND(active_mg(icpu,ilevel)%f,ncache*twotondim, &
            & MPI_INTEGER,icpu-1,tag,MPI_COMM_WORLD,reqsend(countsend),info)
     end if
  end do

  ! Wait for full completion of receives
  call MPI_WAITALL(countrecv,reqrecv,statuses,info)

  ! Gather emission array
  do icpu=1,ncpu
     if (emission_mg(icpu,ilevel)%ngrid>0) then
        do j=1,twotondim
           step=(j-1)*emission_mg(icpu,ilevel)%ngrid
           iskip=(j-1)*active_mg(myid,ilevel)%ngrid
           do i=1,emission_mg(icpu,ilevel)%ngrid
              icell=emission_mg(icpu,ilevel)%igrid(i)+iskip
              active_mg(myid,ilevel)%f(icell,1)=active_mg(myid,ilevel)%f(icell,1)+&
                 & emission_mg(icpu,ilevel)%f(i+step,1)
           end do
        end do
     end if
  end do

  ! Wait for full completion of sends
  call MPI_WAITALL(countsend,reqsend,statuses,info)

#endif

111 format('   Entering make_reverse_mg for level ',I2)

end subroutine make_reverse_mg_int

! ########################################################################
! ########################################################################
! Ksection-based MG communication routines
! ########################################################################
! ########################################################################

subroutine make_virtual_mg_dp_ksec(ivar,ilevel)
  use amr_commons
  use poisson_commons
  use ksection
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer,intent(in)::ivar,ilevel
  ! -------------------------------------------------------------------
  ! Ksection-based forward MG ghost zone exchange for double precision.
  ! Packs emission_mg data with metadata, exchanges via ksection tree,
  ! then scatters to active_mg(sender) grids on the receiver.
  ! -------------------------------------------------------------------
  integer::icpu,i,j,idx,ntotal,nrecv,nprops_ksec,sender,ridx,step
  real(dp),allocatable::sendbuf(:,:),recvbuf(:,:)
  integer,allocatable::dest_cpu(:)

#ifndef WITHOUTMPI
  nprops_ksec = twotondim + 2

  ! Count total emission items
  ntotal = 0
  do icpu = 1, ncpu
     ntotal = ntotal + emission_mg(icpu,ilevel)%ngrid
  end do

  ! Pack sendbuf + dest_cpu
  allocate(sendbuf(1:nprops_ksec, 1:max(ntotal,1)))
  allocate(dest_cpu(1:max(ntotal,1)))
  idx = 0
  do icpu = 1, ncpu
     do i = 1, emission_mg(icpu,ilevel)%ngrid
        idx = idx + 1
        dest_cpu(idx) = icpu
        do j = 1, twotondim
           step = (j-1)*active_mg(myid,ilevel)%ngrid
           sendbuf(j, idx) = active_mg(myid,ilevel)%u( &
                & emission_mg(icpu,ilevel)%igrid(i) + step, ivar)
        end do
        sendbuf(twotondim+1, idx) = dble(myid)
        sendbuf(twotondim+2, idx) = dble(i)
     end do
  end do

  ! Exchange via ksection tree
  call ksection_exchange_dp(sendbuf, ntotal, dest_cpu, nprops_ksec, &
       & recvbuf, nrecv)

  ! Scatter received data to active_mg(sender) on this CPU
  do i = 1, nrecv
     sender = nint(recvbuf(twotondim+1, i))
     ridx   = nint(recvbuf(twotondim+2, i))
     do j = 1, twotondim
        step = (j-1)*active_mg(sender,ilevel)%ngrid
        active_mg(sender,ilevel)%u(ridx + step, ivar) = recvbuf(j, i)
     end do
  end do

  deallocate(sendbuf, dest_cpu, recvbuf)
#endif

end subroutine make_virtual_mg_dp_ksec

! ########################################################################
! ########################################################################

subroutine make_virtual_mg_int_ksec(ilevel)
  use amr_commons
  use poisson_commons
  use ksection
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer,intent(in)::ilevel
  ! -------------------------------------------------------------------
  ! Ksection-based forward MG ghost zone exchange for integer arrays.
  ! Converts int to dp, exchanges via ksection tree, converts back.
  ! -------------------------------------------------------------------
  integer::icpu,i,j,idx,ntotal,nrecv,nprops_ksec,sender,ridx,step
  real(dp),allocatable::sendbuf(:,:),recvbuf(:,:)
  integer,allocatable::dest_cpu(:)

#ifndef WITHOUTMPI
  nprops_ksec = twotondim + 2

  ! Count total emission items
  ntotal = 0
  do icpu = 1, ncpu
     ntotal = ntotal + emission_mg(icpu,ilevel)%ngrid
  end do

  ! Pack sendbuf (int->dp) + dest_cpu
  allocate(sendbuf(1:nprops_ksec, 1:max(ntotal,1)))
  allocate(dest_cpu(1:max(ntotal,1)))
  idx = 0
  do icpu = 1, ncpu
     do i = 1, emission_mg(icpu,ilevel)%ngrid
        idx = idx + 1
        dest_cpu(idx) = icpu
        do j = 1, twotondim
           step = (j-1)*active_mg(myid,ilevel)%ngrid
           sendbuf(j, idx) = dble(active_mg(myid,ilevel)%f( &
                & emission_mg(icpu,ilevel)%igrid(i) + step, 1))
        end do
        sendbuf(twotondim+1, idx) = dble(myid)
        sendbuf(twotondim+2, idx) = dble(i)
     end do
  end do

  ! Exchange via ksection tree
  call ksection_exchange_dp(sendbuf, ntotal, dest_cpu, nprops_ksec, &
       & recvbuf, nrecv)

  ! Scatter received data (dp->int) to active_mg(sender)%f
  do i = 1, nrecv
     sender = nint(recvbuf(twotondim+1, i))
     ridx   = nint(recvbuf(twotondim+2, i))
     do j = 1, twotondim
        step = (j-1)*active_mg(sender,ilevel)%ngrid
        active_mg(sender,ilevel)%f(ridx + step, 1) = nint(recvbuf(j, i))
     end do
  end do

  deallocate(sendbuf, dest_cpu, recvbuf)
#endif

end subroutine make_virtual_mg_int_ksec

! ########################################################################
! ########################################################################

subroutine make_reverse_mg_dp_ksec(ivar,ilevel)
  use amr_commons
  use poisson_commons
  use ksection
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer,intent(in)::ivar,ilevel
  ! -------------------------------------------------------------------
  ! Ksection-based reverse MG ghost zone exchange for double precision.
  ! Packs active_mg(icpu) data, exchanges via ksection tree,
  ! then accumulates (+=) into local active_mg(myid) using emission_mg.
  ! -------------------------------------------------------------------
  integer::icpu,i,j,idx,ntotal,nrecv,nprops_ksec,sender,ridx,step,icell
  real(dp),allocatable::sendbuf(:,:),recvbuf(:,:)
  integer,allocatable::dest_cpu(:)

#ifndef WITHOUTMPI
  nprops_ksec = twotondim + 2

  ! Count total items from remote active_mg
  ntotal = 0
  do icpu = 1, ncpu
     if(icpu == myid) cycle
     ntotal = ntotal + active_mg(icpu,ilevel)%ngrid
  end do

  ! Pack sendbuf from active_mg(icpu) + dest_cpu
  allocate(sendbuf(1:nprops_ksec, 1:max(ntotal,1)))
  allocate(dest_cpu(1:max(ntotal,1)))
  idx = 0
  do icpu = 1, ncpu
     if(icpu == myid) cycle
     do i = 1, active_mg(icpu,ilevel)%ngrid
        idx = idx + 1
        dest_cpu(idx) = icpu
        do j = 1, twotondim
           step = (j-1)*active_mg(icpu,ilevel)%ngrid
           sendbuf(j, idx) = active_mg(icpu,ilevel)%u(i + step, ivar)
        end do
        sendbuf(twotondim+1, idx) = dble(myid)
        sendbuf(twotondim+2, idx) = dble(i)
     end do
  end do

  ! Exchange via ksection tree
  call ksection_exchange_dp(sendbuf, ntotal, dest_cpu, nprops_ksec, &
       & recvbuf, nrecv)

  ! Accumulate received data into local active_mg(myid) using emission_mg
  do i = 1, nrecv
     sender = nint(recvbuf(twotondim+1, i))
     ridx   = nint(recvbuf(twotondim+2, i))
     do j = 1, twotondim
        step  = (j-1)*active_mg(myid,ilevel)%ngrid
        icell = emission_mg(sender,ilevel)%igrid(ridx) + step
        active_mg(myid,ilevel)%u(icell, ivar) = &
             & active_mg(myid,ilevel)%u(icell, ivar) + recvbuf(j, i)
     end do
  end do

  deallocate(sendbuf, dest_cpu, recvbuf)
#endif

end subroutine make_reverse_mg_dp_ksec

! ########################################################################
! ########################################################################

subroutine make_reverse_mg_int_ksec(ilevel)
  use amr_commons
  use poisson_commons
  use ksection
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer,intent(in)::ilevel
  ! -------------------------------------------------------------------
  ! Ksection-based reverse MG ghost zone exchange for integer arrays.
  ! Packs active_mg(icpu)%f, exchanges via ksection tree,
  ! then accumulates (+=) into local active_mg(myid)%f using emission_mg.
  ! -------------------------------------------------------------------
  integer::icpu,i,j,idx,ntotal,nrecv,nprops_ksec,sender,ridx,step,icell
  real(dp),allocatable::sendbuf(:,:),recvbuf(:,:)
  integer,allocatable::dest_cpu(:)

#ifndef WITHOUTMPI
  nprops_ksec = twotondim + 2

  ! Count total items from remote active_mg
  ntotal = 0
  do icpu = 1, ncpu
     if(icpu == myid) cycle
     ntotal = ntotal + active_mg(icpu,ilevel)%ngrid
  end do

  ! Pack sendbuf from active_mg(icpu)%f (int->dp) + dest_cpu
  allocate(sendbuf(1:nprops_ksec, 1:max(ntotal,1)))
  allocate(dest_cpu(1:max(ntotal,1)))
  idx = 0
  do icpu = 1, ncpu
     if(icpu == myid) cycle
     do i = 1, active_mg(icpu,ilevel)%ngrid
        idx = idx + 1
        dest_cpu(idx) = icpu
        do j = 1, twotondim
           step = (j-1)*active_mg(icpu,ilevel)%ngrid
           sendbuf(j, idx) = dble(active_mg(icpu,ilevel)%f(i + step, 1))
        end do
        sendbuf(twotondim+1, idx) = dble(myid)
        sendbuf(twotondim+2, idx) = dble(i)
     end do
  end do

  ! Exchange via ksection tree
  call ksection_exchange_dp(sendbuf, ntotal, dest_cpu, nprops_ksec, &
       & recvbuf, nrecv)

  ! Accumulate received data (dp->int) into active_mg(myid)%f
  do i = 1, nrecv
     sender = nint(recvbuf(twotondim+1, i))
     ridx   = nint(recvbuf(twotondim+2, i))
     do j = 1, twotondim
        step  = (j-1)*active_mg(myid,ilevel)%ngrid
        icell = emission_mg(sender,ilevel)%igrid(ridx) + step
        active_mg(myid,ilevel)%f(icell, 1) = &
             & active_mg(myid,ilevel)%f(icell, 1) + nint(recvbuf(j, i))
     end do
  end do

  deallocate(sendbuf, dest_cpu, recvbuf)
#endif

end subroutine make_reverse_mg_int_ksec

! ########################################################################
! ########################################################################
! ########################################################################
! ########################################################################

subroutine dump_mg_levels(ilevel,idout)
   use amr_commons
   use poisson_commons
   implicit none
#ifndef WITHOUTMPI
   include 'mpif.h'  
#endif
   integer, intent(in) :: idout, ilevel

   character(len=24)  :: cfile
   character(len=5)   :: ccpu='00000'
   character(len=5)   :: cout='00000'

   integer :: i, ngrids, igrid, icpu, idim
   
   integer,parameter::tag=1119
   integer::dummy_io,info2

   write(ccpu,'(I5.5)') myid
   write(cout,'(I5.5)') idout
   cfile='multigrid_'//cout//'.out'//ccpu
   
   ! Wait for the token
#ifndef WITHOUTMPI
   if(IOGROUPSIZE>0) then
      if (mod(myid-1,IOGROUPSIZE)/=0) then
         call MPI_RECV(dummy_io,1,MPI_INTEGER,myid-1-1,tag,&
              & MPI_COMM_WORLD,MPI_STATUS_IGNORE,info2)
      end if
   endif
#endif
   
   open(unit=10,file=cfile,status='unknown',form='formatted')

   write(10,'(I1)') ndim
   write(10,'(I1)') myid
   write(10,'(I1)') ncpu
   write(10,'(I2)') ilevel

   ! Loop over levels
   do i=1,ilevel-1
      ! Active grids
      ngrids=active_mg(myid,i)%ngrid
      write(10,*) ngrids
      do igrid=1,ngrids
         do idim=1,ndim
            write(10,*) xg(active_mg(myid,i)%igrid(igrid),idim)
         end do
      end do

      ! Reception grids
      do icpu=1,ncpu
         if(icpu==myid)cycle
         ngrids=active_mg(icpu,i)%ngrid
         write(10,*) ngrids
         do igrid=1,ngrids
            do idim=1,ndim
               write(10,*) xg(active_mg(icpu,i)%igrid(igrid),idim)
            end do
         end do
      end do

   end do

   close(10)

        ! Send the token
#ifndef WITHOUTMPI
   if(IOGROUPSIZE>0) then
      if(mod(myid,IOGROUPSIZE)/=0 .and.(myid.lt.ncpu))then
         dummy_io=1
         call MPI_SEND(dummy_io,1,MPI_INTEGER,myid-1+1,tag, &
              & MPI_COMM_WORLD,info2)
      end if
   endif
#endif

end subroutine dump_mg_levels


#ifdef HYDRO_CUDA
! ########################################################################
! cuFFT direct Poisson solver for fully uniform levels (periodic BC)
! Replaces MG V-cycle iteration with single FFT solve: O(N log N) vs O(N * niter)
! ########################################################################

subroutine fft_poisson_solve_uniform(ilevel, icount)
   use amr_commons
   use poisson_commons
   use poisson_parameters
   use poisson_cuda_interface
   use iso_c_binding

   implicit none
#ifndef WITHOUTMPI
   include "mpif.h"
#endif

   integer, intent(in) :: ilevel, icount

   ! Grid dimensions
   integer :: fft_Nx, fft_Ny, fft_Nz
   integer(i8b) :: N_total
   real(dp) :: dx_fft, dx2_fft

   ! Persistent arrays (saved across calls to avoid reallocation)
   integer,  allocatable, save :: fft_map(:)
   real(dp), allocatable, save :: rhs_3d(:), rhs_local(:)
   integer(i8b), save :: saved_N_total = 0

   ! Loop variables
   integer :: igrid, ngrid_loc, ind, iskip
   integer :: icell_amr, igrid_amr
   integer :: ix, iy, iz, idx_3d
   integer :: Kx, Ky, Kz
   integer :: info

   ! Grid dimensions for this level (cell grid, not oct grid)
   fft_Nx = nx * 2**ilevel
   fft_Ny = ny * 2**ilevel
   fft_Nz = nz * 2**ilevel
   N_total = int(fft_Nx, i8b) * int(fft_Ny, i8b) * int(fft_Nz, i8b)

   dx_fft  = 0.5d0**ilevel
   dx2_fft = dx_fft * dx_fft

   if(myid==1) write(*,'(A,I3,A,I5,A,I5,A,I5,A,I15)') &
        ' cuFFT Poisson: level=', ilevel, &
        ' grid=', fft_Nx, 'x', fft_Ny, 'x', fft_Nz, &
        ' N=', N_total

   ! ------------------------------------------------------------------
   ! Persistent allocation (only allocate on first call or size change)
   ! ------------------------------------------------------------------
   if(N_total /= saved_N_total) then
      if(allocated(fft_map))   deallocate(fft_map)
      if(allocated(rhs_local)) deallocate(rhs_local)
      if(allocated(rhs_3d))    deallocate(rhs_3d)
      allocate(fft_map(0:N_total-1))
      allocate(rhs_local(0:N_total-1))
      allocate(rhs_3d(0:N_total-1))
      saved_N_total = N_total
   end if

   ! ------------------------------------------------------------------
   ! Step 1: Build fft_map and gather local RHS
   ! fft_map(idx_3d) = icell_amr (1-based), idx_3d is 0-based C index
   ! ------------------------------------------------------------------
   fft_map   = 0
   rhs_local = 0.0d0

   ngrid_loc = active(ilevel)%ngrid
   do igrid = 1, ngrid_loc
      igrid_amr = active(ilevel)%igrid(igrid)

      ! Grid center coordinates -> integer cell coords
      Kx = nint(xg(igrid_amr, 1) * dble(fft_Nx))
      Ky = nint(xg(igrid_amr, 2) * dble(fft_Ny))
      Kz = nint(xg(igrid_amr, 3) * dble(fft_Nz))

      do ind = 1, twotondim
         ! Child cell offset within oct (0-based)
         ix = Kx - 1 + mod(ind-1, 2)
         iy = Ky - 1 + mod((ind-1)/2, 2)
         iz = Kz - 1 + (ind-1)/4

         ! Periodic wrapping
         ix = modulo(ix, fft_Nx)
         iy = modulo(iy, fft_Ny)
         iz = modulo(iz, fft_Nz)

         ! C row-major index: idx = ix * Ny * Nz + iy * Nz + iz
         idx_3d = ix * fft_Ny * fft_Nz + iy * fft_Nz + iz

         ! RAMSES cell index
         iskip = ncoarse + (ind - 1) * ngridmax
         icell_amr = iskip + igrid_amr

         fft_map(idx_3d) = icell_amr
         rhs_local(idx_3d) = f(icell_amr, 2)
      end do
   end do

   ! ------------------------------------------------------------------
   ! Step 2: MPI_Allreduce to get global RHS (all ranks contribute local parts)
   ! ------------------------------------------------------------------
#ifndef WITHOUTMPI
   call MPI_ALLREDUCE(rhs_local, rhs_3d, int(N_total), &
        MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, info)
#else
   rhs_3d = rhs_local
#endif

   ! Note: fft_map stays LOCAL (each rank's own cell indices)
   ! The scatter kernel uses local fft_map to write to local d_mg_phi

   ! ------------------------------------------------------------------
   ! Step 3: cuFFT setup (plans + Green's function, only on first call or grid change)
   ! ------------------------------------------------------------------
   call cuda_fft_poisson_setup_c(fft_map, &
        int(fft_Nx, c_int), int(fft_Ny, c_int), int(fft_Nz, c_int), dx2_fft)

   ! ------------------------------------------------------------------
   ! Step 4: cuFFT solve (upload RHS, FFT, Green, IFFT, scatter to d_mg_phi)
   ! ------------------------------------------------------------------
   call cuda_fft_poisson_solve_c(rhs_3d, int(N_total, c_int))

end subroutine fft_poisson_solve_uniform
#endif

! ########################################################################
! ########################################################################
! ########################################################################
! ########################################################################
