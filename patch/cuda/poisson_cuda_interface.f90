!###########################################################
! Fortran interface to CUDA Poisson MG kernels for cuRAMSES
! Provides ISO_C_BINDING wrappers for poisson_cuda_kernels.cu
!###########################################################
module poisson_cuda_interface
  use iso_c_binding
  implicit none

  ! C interface to Poisson CUDA kernels
  ! Arrays passed by reference (no VALUE) so Fortran passes pointer to C
  interface
     subroutine cuda_mg_upload_c(phi, f, flag2, ncell, nbor, igrid, ngrid) &
          bind(C, name='cuda_mg_upload')
       import :: c_double, c_int, c_long_long
       real(c_double) :: phi(*)
       real(c_double) :: f(*)
       integer(c_int) :: flag2(*)
       integer(c_long_long), value :: ncell
       integer(c_int) :: nbor(*)
       integer(c_int) :: igrid(*)
       integer(c_int), value :: ngrid
     end subroutine

     subroutine cuda_mg_download_phi_c(phi, ncell) &
          bind(C, name='cuda_mg_download_phi')
       import :: c_double, c_long_long
       real(c_double) :: phi(*)
       integer(c_long_long), value :: ncell
     end subroutine

     subroutine cuda_mg_upload_phi_c(phi, ncell) &
          bind(C, name='cuda_mg_upload_phi')
       import :: c_double, c_long_long
       real(c_double) :: phi(*)
       integer(c_long_long), value :: ncell
     end subroutine

     subroutine cuda_mg_download_f1_c(f1, ncell) &
          bind(C, name='cuda_mg_download_f1')
       import :: c_double, c_long_long
       real(c_double) :: f1(*)
       integer(c_long_long), value :: ncell
     end subroutine

     subroutine cuda_mg_gauss_seidel_c(ngrid, ngridmax_c, ncoarse_c, &
          dx2, color, safe_mode) &
          bind(C, name='cuda_mg_gauss_seidel')
       import :: c_int, c_double
       integer(c_int), value :: ngrid, ngridmax_c, ncoarse_c
       real(c_double), value :: dx2
       integer(c_int), value :: color, safe_mode
     end subroutine

     subroutine cuda_mg_residual_c(ngrid, ngridmax_c, ncoarse_c, &
          oneoverdx2, dtwondim, dx2_norm, norm2, compute_norm) &
          bind(C, name='cuda_mg_residual')
       import :: c_int, c_double
       integer(c_int), value :: ngrid, ngridmax_c, ncoarse_c
       real(c_double), value :: oneoverdx2, dtwondim, dx2_norm
       real(c_double) :: norm2
       integer(c_int), value :: compute_norm
     end subroutine

     subroutine cuda_mg_free_c() bind(C, name='cuda_mg_free')
     end subroutine

     subroutine cuda_mg_finalize_c() bind(C, name='cuda_mg_finalize')
     end subroutine

     integer(c_int) function cuda_mg_is_ready_c() &
          bind(C, name='cuda_mg_is_ready')
       import :: c_int
     end function

     subroutine cuda_mg_halo_setup_c(emit_cells, n_emit, &
          recv_cells, n_recv) &
          bind(C, name='cuda_mg_halo_setup')
       import :: c_int
       integer(c_int) :: emit_cells(*)
       integer(c_int), value :: n_emit
       integer(c_int) :: recv_cells(*)
       integer(c_int), value :: n_recv
     end subroutine

     subroutine cuda_mg_halo_gather_c(h_emit_buf, n) &
          bind(C, name='cuda_mg_halo_gather')
       import :: c_double, c_int
       real(c_double) :: h_emit_buf(*)
       integer(c_int), value :: n
     end subroutine

     subroutine cuda_mg_halo_scatter_c(h_recv_buf, n) &
          bind(C, name='cuda_mg_halo_scatter')
       import :: c_double, c_int
       real(c_double) :: h_recv_buf(*)
       integer(c_int), value :: n
     end subroutine

     subroutine cuda_mg_halo_free_c() &
          bind(C, name='cuda_mg_halo_free')
     end subroutine

     ! Restrict/Interp GPU setup: upload mapping arrays + bbb/ccc
     subroutine cuda_mg_ri_setup_c(restrict_target, interp_nbor_flat, &
          ngrid_fine, total_coarse_cells, bbb, ccc) &
          bind(C, name='cuda_mg_ri_setup')
       import :: c_int, c_double
       integer(c_int) :: restrict_target(*)
       integer(c_int) :: interp_nbor_flat(*)
       integer(c_int), value :: ngrid_fine, total_coarse_cells
       real(c_double) :: bbb(*)
       integer(c_int) :: ccc(*)
     end subroutine

     ! Restrict: memset coarse RHS + run kernel
     subroutine cuda_mg_restrict_execute_c(ngrid, ngridmax_c, ncoarse_c) &
          bind(C, name='cuda_mg_restrict_execute')
       import :: c_int
       integer(c_int), value :: ngrid, ngridmax_c, ncoarse_c
     end subroutine

     ! Restrict: download coarse RHS to host
     subroutine cuda_mg_restrict_download_c(h_coarse_rhs, total_cells) &
          bind(C, name='cuda_mg_restrict_download')
       import :: c_double, c_int
       real(c_double) :: h_coarse_rhs(*)
       integer(c_int), value :: total_cells
     end subroutine

     ! Interp: upload coarse phi to GPU
     subroutine cuda_mg_interp_upload_c(h_coarse_phi, total_cells) &
          bind(C, name='cuda_mg_interp_upload')
       import :: c_double, c_int
       real(c_double) :: h_coarse_phi(*)
       integer(c_int), value :: total_cells
     end subroutine

     ! Interp: run kernel
     subroutine cuda_mg_interp_execute_c(ngrid, ngridmax_c, ncoarse_c) &
          bind(C, name='cuda_mg_interp_execute')
       import :: c_int
       integer(c_int), value :: ngrid, ngridmax_c, ncoarse_c
     end subroutine

     ! Check if restrict/interp GPU is ready
     integer(c_int) function cuda_mg_ri_is_ready_c() &
          bind(C, name='cuda_mg_ri_is_ready')
       import :: c_int
     end function

     ! Free restrict/interp GPU arrays
     subroutine cuda_mg_ri_free_c() &
          bind(C, name='cuda_mg_ri_free')
     end subroutine

     ! cuFFT Poisson: setup plans, Green's function, fft_map
     subroutine cuda_fft_poisson_setup_c(fft_map, Nx, Ny, Nz, dx2) &
          bind(C, name='cuda_fft_poisson_setup')
       import :: c_int, c_double
       integer(c_int) :: fft_map(*)
       integer(c_int), value :: Nx, Ny, Nz
       real(c_double), value :: dx2
     end subroutine

     ! cuFFT Poisson: solve (h_rhs_3d → d_mg_phi)
     subroutine cuda_fft_poisson_solve_c(h_rhs_3d, N_real) &
          bind(C, name='cuda_fft_poisson_solve')
       import :: c_double, c_int
       real(c_double) :: h_rhs_3d(*)
       integer(c_int), value :: N_real
     end subroutine

     ! cuFFT Poisson: free plans and arrays
     subroutine cuda_fft_poisson_free_c() &
          bind(C, name='cuda_fft_poisson_free')
     end subroutine
  end interface

end module poisson_cuda_interface
