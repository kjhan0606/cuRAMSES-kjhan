!###########################################################
! Fortran interface to CUDA hydro kernels for cuRAMSES
! Provides ISO_C_BINDING wrappers for hydro_cuda_kernels.cu
!###########################################################
module hydro_cuda_interface
  use iso_c_binding
  implicit none

  logical :: hydro_cuda_initialized = .false.

  ! Riemann solver IDs (must match hydro_cuda_kernels.cu)
  integer, parameter :: RIEMANN_LLF  = 0
  integer, parameter :: RIEMANN_HLL  = 1
  integer, parameter :: RIEMANN_HLLC = 2

  ! C interface to hydro CUDA kernels
  interface
     subroutine hydro_cuda_init_c( &
          gamma, smallr, smallc, &
          slope_type, slope_theta, &
          nvar, ndim, riemann_solver, &
          pressure_fix, difmag) &
          bind(C, name='hydro_cuda_init')
       import :: c_double, c_int
       real(c_double), value :: gamma, smallr, smallc
       integer(c_int), value :: slope_type
       real(c_double), value :: slope_theta
       integer(c_int), value :: nvar, ndim, riemann_solver
       integer(c_int), value :: pressure_fix
       real(c_double), value :: difmag
     end subroutine

     subroutine hydro_cuda_unsplit_async_c( &
          h_uloc, h_gloc, h_flux, h_tmp, h_ok, &
          dx, dy, dz, dt, ngrid, stream_slot) &
          bind(C, name='hydro_cuda_unsplit_async')
       import :: c_double, c_int, c_ptr
       type(c_ptr), value :: h_uloc, h_gloc, h_flux, h_tmp, h_ok
       real(c_double), value :: dx, dy, dz, dt
       integer(c_int), value :: ngrid, stream_slot
     end subroutine

     subroutine hydro_cuda_unsplit_sync_c( &
          h_flux, h_tmp, ngrid, stream_slot) &
          bind(C, name='hydro_cuda_unsplit_sync')
       import :: c_double, c_int, c_ptr
       type(c_ptr), value :: h_flux, h_tmp
       integer(c_int), value :: ngrid, stream_slot
     end subroutine
  end interface

contains

  !-----------------------------------------------------------
  ! Initialize CUDA hydro solver with current parameters
  ! Called once after hydro_parameters are read
  !-----------------------------------------------------------
  subroutine hydro_cuda_init_f()
    use hydro_parameters, only: gamma, smallr, smallc, &
         slope_type, slope_theta, difmag, riemann, nvar, ndim
    use amr_parameters, only: pressure_fix
    implicit none
    integer :: riemann_id, pfix_int

    ! Map Fortran riemann string to C integer ID
    select case(trim(riemann))
    case('llf')
       riemann_id = RIEMANN_LLF
    case('hll')
       riemann_id = RIEMANN_HLL
    case('hllc')
       riemann_id = RIEMANN_HLLC
    case default
       ! Unsupported solver on GPU — fall back to LLF
       riemann_id = RIEMANN_LLF
    end select

    ! Convert logical to integer for C interface
    pfix_int = 0
    if (pressure_fix) pfix_int = 1

    call hydro_cuda_init_c( &
         real(gamma, c_double), &
         real(smallr, c_double), &
         real(smallc, c_double), &
         int(slope_type, c_int), &
         real(slope_theta, c_double), &
         int(nvar, c_int), &
         int(ndim, c_int), &
         int(riemann_id, c_int), &
         int(pfix_int, c_int), &
         real(difmag, c_double))

    hydro_cuda_initialized = .true.
  end subroutine hydro_cuda_init_f

  !-----------------------------------------------------------
  ! Launch async GPU hydro solve for a batch of grids
  ! Arrays must be contiguous (Fortran default)
  ! Note: ok is not used by CUDA kernels, c_null_ptr is passed
  !-----------------------------------------------------------
  subroutine hydro_cuda_unsplit_async_f( &
       uloc, gravin, flux, tmp, &
       dx, dy, dz, dt, ngrid, stream_slot)
    implicit none
    real(c_double), intent(in), target :: uloc(*)
    real(c_double), intent(in), target :: gravin(*)
    real(c_double), intent(inout), target :: flux(*)
    real(c_double), intent(inout), target :: tmp(*)
    real(c_double), intent(in) :: dx, dy, dz, dt
    integer(c_int), intent(in) :: ngrid, stream_slot

    call hydro_cuda_unsplit_async_c( &
         c_loc(uloc(1)), c_loc(gravin(1)), &
         c_loc(flux(1)), c_loc(tmp(1)), c_null_ptr, &
         dx, dy, dz, dt, ngrid, stream_slot)
  end subroutine hydro_cuda_unsplit_async_f

  !-----------------------------------------------------------
  ! Wait for GPU hydro solve to complete
  ! After return, flux and tmp contain valid results
  !-----------------------------------------------------------
  subroutine hydro_cuda_unsplit_sync_f( &
       flux, tmp, ngrid, stream_slot)
    implicit none
    real(c_double), intent(inout), target :: flux(*)
    real(c_double), intent(inout), target :: tmp(*)
    integer(c_int), intent(in) :: ngrid, stream_slot

    call hydro_cuda_unsplit_sync_c( &
         c_loc(flux(1)), c_loc(tmp(1)), ngrid, stream_slot)
  end subroutine hydro_cuda_unsplit_sync_f

end module hydro_cuda_interface
