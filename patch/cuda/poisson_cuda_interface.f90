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
  end interface

end module poisson_cuda_interface
