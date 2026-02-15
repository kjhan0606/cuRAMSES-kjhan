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
          dx, dy, dz, dt, ngrid, stride, stream_slot) &
          bind(C, name='hydro_cuda_unsplit_async')
       import :: c_double, c_int, c_ptr
       type(c_ptr), value :: h_uloc, h_gloc, h_flux, h_tmp, h_ok
       real(c_double), value :: dx, dy, dz, dt
       integer(c_int), value :: ngrid, stride, stream_slot
     end subroutine

     subroutine hydro_cuda_unsplit_sync_c( &
          h_flux, h_tmp, ngrid, stream_slot) &
          bind(C, name='hydro_cuda_unsplit_sync')
       import :: c_double, c_int, c_ptr
       type(c_ptr), value :: h_flux, h_tmp
       integer(c_int), value :: ngrid, stream_slot
     end subroutine

     ! GPU-gather variant: stencil indices -> GPU gather -> compute
     subroutine cuda_mesh_upload_c(uold, f_grav, son, ncell, nvar, ndim) &
          bind(C, name='cuda_mesh_upload')
       import :: c_ptr, c_long_long, c_int
       type(c_ptr), value :: uold, f_grav, son
       integer(c_long_long), value :: ncell
       integer(c_int), value :: nvar, ndim
     end subroutine

     subroutine cuda_mesh_free_c() bind(C, name='cuda_mesh_free')
     end subroutine

     subroutine hydro_cuda_gather_unsplit_async_c( &
          h_stencil_idx, h_stencil_grav, h_interp_vals, &
          h_flux, h_tmp, &
          dx, dy, dz, dt, ngrid, stride, n_interp, stream_slot) &
          bind(C, name='hydro_cuda_gather_unsplit_async')
       import :: c_double, c_int, c_ptr
       type(c_ptr), value :: h_stencil_idx, h_stencil_grav, h_interp_vals
       type(c_ptr), value :: h_flux, h_tmp
       real(c_double), value :: dx, dy, dz, dt
       integer(c_int), value :: ngrid, stride, n_interp, stream_slot
     end subroutine

     subroutine hydro_cuda_gather_unsplit_sync_c( &
          h_flux, h_tmp, ngrid, stream_slot) &
          bind(C, name='hydro_cuda_gather_unsplit_sync')
       import :: c_double, c_int, c_ptr
       type(c_ptr), value :: h_flux, h_tmp
       integer(c_int), value :: ngrid, stream_slot
     end subroutine
  end interface

contains

  !-----------------------------------------------------------
  ! Initialize CUDA hydro solver with current parameters
  !-----------------------------------------------------------
  subroutine hydro_cuda_init_f()
    use hydro_parameters, only: gamma, smallr, smallc, &
         slope_type, slope_theta, difmag, riemann, nvar, ndim
    use amr_parameters, only: pressure_fix
    implicit none
    integer :: riemann_id, pfix_int

    select case(trim(riemann))
    case('llf')
       riemann_id = RIEMANN_LLF
    case('hll')
       riemann_id = RIEMANN_HLL
    case('hllc')
       riemann_id = RIEMANN_HLLC
    case default
       riemann_id = RIEMANN_LLF
    end select

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
  ! Launch async GPU hydro solve (CPU-gathered data path)
  !-----------------------------------------------------------
  subroutine hydro_cuda_unsplit_async_f( &
       uloc, gravin, flux, tmp, &
       dx, dy, dz, dt, ngrid, stride, stream_slot)
    implicit none
    real(c_double), intent(in), target :: uloc(*)
    real(c_double), intent(in), target :: gravin(*)
    real(c_double), intent(inout), target :: flux(*)
    real(c_double), intent(inout), target :: tmp(*)
    real(c_double), intent(in) :: dx, dy, dz, dt
    integer(c_int), intent(in) :: ngrid, stride, stream_slot

    call hydro_cuda_unsplit_async_c( &
         c_loc(uloc(1)), c_loc(gravin(1)), &
         c_loc(flux(1)), c_loc(tmp(1)), c_null_ptr, &
         dx, dy, dz, dt, ngrid, stride, stream_slot)
  end subroutine hydro_cuda_unsplit_async_f

  !-----------------------------------------------------------
  ! Wait for GPU hydro solve to complete
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

  !-----------------------------------------------------------
  ! Upload mesh arrays to GPU (once per step)
  !-----------------------------------------------------------
  subroutine cuda_mesh_upload_f(uold_arr, f_arr, son_arr, ncell, nvar_in, ndim_in, has_grav)
    implicit none
    real(c_double), intent(in), target :: uold_arr(*)
    real(c_double), intent(in), target :: f_arr(*)
    integer(c_int), intent(in), target :: son_arr(*)
    integer(c_long_long), intent(in) :: ncell
    integer(c_int), intent(in) :: nvar_in, ndim_in
    logical, intent(in) :: has_grav
    type(c_ptr) :: f_ptr

    if (has_grav) then
       f_ptr = c_loc(f_arr(1))
    else
       f_ptr = c_null_ptr
    end if

    call cuda_mesh_upload_c(c_loc(uold_arr(1)), f_ptr, c_loc(son_arr(1)), &
         ncell, nvar_in, ndim_in)
  end subroutine cuda_mesh_upload_f

  !-----------------------------------------------------------
  ! GPU-gather async: upload stencil indices, GPU gather+compute
  !-----------------------------------------------------------
  subroutine hydro_cuda_gather_unsplit_async_f( &
       stencil_idx, stencil_grav, interp_vals, flux, tmp, &
       dx, dy, dz, dt, ngrid, stride, n_interp, stream_slot)
    implicit none
    integer(c_int), intent(in), target :: stencil_idx(*), stencil_grav(*)
    real(c_double), intent(in), target :: interp_vals(*)
    real(c_double), intent(inout), target :: flux(*), tmp(*)
    real(c_double), intent(in) :: dx, dy, dz, dt
    integer(c_int), intent(in) :: ngrid, stride, n_interp, stream_slot
    type(c_ptr) :: interp_ptr

    if (n_interp > 0) then
       interp_ptr = c_loc(interp_vals(1))
    else
       interp_ptr = c_null_ptr
    end if

    call hydro_cuda_gather_unsplit_async_c( &
         c_loc(stencil_idx(1)), c_loc(stencil_grav(1)), interp_ptr, &
         c_loc(flux(1)), c_loc(tmp(1)), &
         dx, dy, dz, dt, ngrid, stride, n_interp, stream_slot)
  end subroutine hydro_cuda_gather_unsplit_async_f

  !-----------------------------------------------------------
  ! GPU-gather sync: wait for completion
  !-----------------------------------------------------------
  subroutine hydro_cuda_gather_unsplit_sync_f( &
       flux, tmp, ngrid, stream_slot)
    implicit none
    real(c_double), intent(inout), target :: flux(*)
    real(c_double), intent(inout), target :: tmp(*)
    integer(c_int), intent(in) :: ngrid, stream_slot

    call hydro_cuda_gather_unsplit_sync_c( &
         c_loc(flux(1)), c_loc(tmp(1)), ngrid, stream_slot)
  end subroutine hydro_cuda_gather_unsplit_sync_f

end module hydro_cuda_interface
