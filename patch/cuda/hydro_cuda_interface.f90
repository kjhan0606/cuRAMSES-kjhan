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

     ! Scatter-reduce: 5 kernels + scatter_reduce, compact D2H
     subroutine hydro_cuda_unsplit_reduce_async_c( &
          h_uloc, h_gloc, h_ok, &
          h_add_unew, h_add_lm1, &
          h_add_divu_l, h_add_enew_l, &
          h_add_divu_lm1, h_add_enew_lm1, &
          dx, dy, dz, dt, ngrid, stride, stream_slot) &
          bind(C, name='hydro_cuda_unsplit_reduce_async')
       import :: c_double, c_int, c_ptr
       type(c_ptr), value :: h_uloc, h_gloc, h_ok
       type(c_ptr), value :: h_add_unew, h_add_lm1
       type(c_ptr), value :: h_add_divu_l, h_add_enew_l
       type(c_ptr), value :: h_add_divu_lm1, h_add_enew_lm1
       real(c_double), value :: dx, dy, dz, dt
       integer(c_int), value :: ngrid, stride, stream_slot
     end subroutine

     subroutine hydro_cuda_unsplit_reduce_sync_c(ngrid, stream_slot) &
          bind(C, name='hydro_cuda_unsplit_reduce_sync')
       import :: c_int
       integer(c_int), value :: ngrid, stream_slot
     end subroutine

     ! gradient_phi CUDA kernel
     subroutine gradient_phi_cuda_async_c( &
          h_phi_buf, h_f_buf, a, b, ngrid, stream_slot) &
          bind(C, name='gradient_phi_cuda_async')
       import :: c_double, c_int, c_ptr
       type(c_ptr), value :: h_phi_buf, h_f_buf
       real(c_double), value :: a, b
       integer(c_int), value :: ngrid, stream_slot
     end subroutine

     subroutine gradient_phi_cuda_sync_c(stream_slot) &
          bind(C, name='gradient_phi_cuda_sync')
       import :: c_int
       integer(c_int), value :: stream_slot
     end subroutine

     ! upload_fine CUDA kernel
     subroutine upload_fine_cuda_async_c( &
          h_child_buf, h_parent_buf, &
          nsplit, nvar, ndim, smallr, do_eint, stream_slot) &
          bind(C, name='upload_fine_cuda_async')
       import :: c_double, c_int, c_ptr
       type(c_ptr), value :: h_child_buf, h_parent_buf
       integer(c_int), value :: nsplit, nvar, ndim, do_eint, stream_slot
       real(c_double), value :: smallr
     end subroutine

     subroutine upload_fine_cuda_sync_c(stream_slot) &
          bind(C, name='upload_fine_cuda_sync')
       import :: c_int
       integer(c_int), value :: stream_slot
     end subroutine

     ! synchro_hydro CUDA kernel
     subroutine synchro_cuda_async_c( &
          h_buf, smallr, dteff, ncell, ndim, stride, stream_slot) &
          bind(C, name='synchro_cuda_async')
       import :: c_double, c_int, c_ptr
       type(c_ptr), value :: h_buf
       real(c_double), value :: smallr, dteff
       integer(c_int), value :: ncell, ndim, stride, stream_slot
     end subroutine

     subroutine synchro_cuda_sync_c(stream_slot) &
          bind(C, name='synchro_cuda_sync')
       import :: c_int
       integer(c_int), value :: stream_slot
     end subroutine

     ! cmpdt CUDA kernel
     subroutine cmpdt_cuda_async_c( &
          h_buf, h_dt_buf, &
          dx, smallr, smallc, gamma_val, courant_factor, &
          ncell, nvar, ndim, stride, stream_slot) &
          bind(C, name='cmpdt_cuda_async')
       import :: c_double, c_int, c_ptr
       type(c_ptr), value :: h_buf, h_dt_buf
       real(c_double), value :: dx, smallr, smallc, gamma_val, courant_factor
       integer(c_int), value :: ncell, nvar, ndim, stride, stream_slot
     end subroutine

     subroutine cmpdt_cuda_sync_c(stream_slot) &
          bind(C, name='cmpdt_cuda_sync')
       import :: c_int
       integer(c_int), value :: stream_slot
     end subroutine

     ! Host memory pinning
     function cuda_host_register_c(ptr, nbytes) result(rc) &
          bind(C, name='cuda_host_register')
       import :: c_ptr, c_long_long, c_int
       type(c_ptr), value :: ptr
       integer(c_long_long), value :: nbytes
       integer(c_int) :: rc
     end function

     subroutine cuda_host_unregister_c(ptr) &
          bind(C, name='cuda_host_unregister')
       import :: c_ptr
       type(c_ptr), value :: ptr
     end subroutine

     subroutine cuda_h2d_bandwidth_test_c() &
          bind(C, name='cuda_h2d_bandwidth_test')
     end subroutine

     ! Per-kernel profiling
     subroutine hydro_cuda_profile_accumulate_c(stream_slot, ngrid) &
          bind(C, name='hydro_cuda_profile_accumulate')
       import :: c_int
       integer(c_int), value :: stream_slot, ngrid
     end subroutine

     subroutine hydro_cuda_profile_report_c() &
          bind(C, name='hydro_cuda_profile_report')
     end subroutine

     subroutine hydro_cuda_profile_reset_c() &
          bind(C, name='hydro_cuda_profile_reset')
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
  ! Scatter-reduce: 5 kernels + scatter_reduce, compact D2H
  ! H2D: uloc+gloc+ok (~101 MB), D2H: add_unew+add_lm1 (~5 MB)
  !-----------------------------------------------------------
  subroutine hydro_cuda_unsplit_reduce_async_f( &
       uloc, gravin, ok_int, &
       add_unew, add_lm1, &
       add_divu_l, add_enew_l, add_divu_lm1, add_enew_lm1, &
       dx, dy, dz, dt, ngrid, stride, stream_slot, has_pfix)
    implicit none
    real(c_double), intent(in), target :: uloc(*), gravin(*)
    integer(c_int), intent(in), target :: ok_int(*)
    real(c_double), intent(inout), target :: add_unew(*), add_lm1(*)
    real(c_double), intent(inout), target :: add_divu_l(*), add_enew_l(*)
    real(c_double), intent(inout), target :: add_divu_lm1(*), add_enew_lm1(*)
    real(c_double), intent(in) :: dx, dy, dz, dt
    integer(c_int), intent(in) :: ngrid, stride, stream_slot
    logical, intent(in) :: has_pfix
    type(c_ptr) :: pfix_ptr_dl, pfix_ptr_el, pfix_ptr_dlm1, pfix_ptr_elm1

    if (has_pfix) then
       pfix_ptr_dl   = c_loc(add_divu_l(1))
       pfix_ptr_el   = c_loc(add_enew_l(1))
       pfix_ptr_dlm1 = c_loc(add_divu_lm1(1))
       pfix_ptr_elm1 = c_loc(add_enew_lm1(1))
    else
       pfix_ptr_dl   = c_null_ptr
       pfix_ptr_el   = c_null_ptr
       pfix_ptr_dlm1 = c_null_ptr
       pfix_ptr_elm1 = c_null_ptr
    end if

    call hydro_cuda_unsplit_reduce_async_c( &
         c_loc(uloc(1)), c_loc(gravin(1)), c_loc(ok_int(1)), &
         c_loc(add_unew(1)), c_loc(add_lm1(1)), &
         pfix_ptr_dl, pfix_ptr_el, pfix_ptr_dlm1, pfix_ptr_elm1, &
         dx, dy, dz, dt, ngrid, stride, stream_slot)
  end subroutine hydro_cuda_unsplit_reduce_async_f

  subroutine hydro_cuda_unsplit_reduce_sync_f(ngrid, stream_slot)
    implicit none
    integer(c_int), intent(in) :: ngrid, stream_slot
    call hydro_cuda_unsplit_reduce_sync_c(ngrid, stream_slot)
  end subroutine hydro_cuda_unsplit_reduce_sync_f

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

  !-----------------------------------------------------------
  ! gradient_phi CUDA: async launch
  !-----------------------------------------------------------
  subroutine gradient_phi_cuda_async_f(phi_buf, f_buf, a, b, ngrid, stream_slot)
    implicit none
    real(c_double), intent(in), target :: phi_buf(*)
    real(c_double), intent(inout), target :: f_buf(*)
    real(c_double), intent(in) :: a, b
    integer(c_int), intent(in) :: ngrid, stream_slot
    call gradient_phi_cuda_async_c(c_loc(phi_buf(1)), c_loc(f_buf(1)), &
         a, b, ngrid, stream_slot)
  end subroutine gradient_phi_cuda_async_f

  !-----------------------------------------------------------
  ! gradient_phi CUDA: sync
  !-----------------------------------------------------------
  subroutine gradient_phi_cuda_sync_f(stream_slot)
    implicit none
    integer(c_int), intent(in) :: stream_slot
    call gradient_phi_cuda_sync_c(stream_slot)
  end subroutine gradient_phi_cuda_sync_f

  !-----------------------------------------------------------
  ! upload_fine CUDA: async launch
  !-----------------------------------------------------------
  subroutine upload_fine_cuda_async_f(child_buf, parent_buf, &
       nsplit, nvar_in, ndim_in, smallr_in, do_eint, stream_slot)
    implicit none
    real(c_double), intent(in), target :: child_buf(*)
    real(c_double), intent(inout), target :: parent_buf(*)
    integer(c_int), intent(in) :: nsplit, nvar_in, ndim_in, do_eint, stream_slot
    real(c_double), intent(in) :: smallr_in
    call upload_fine_cuda_async_c(c_loc(child_buf(1)), c_loc(parent_buf(1)), &
         nsplit, nvar_in, ndim_in, smallr_in, do_eint, stream_slot)
  end subroutine upload_fine_cuda_async_f

  !-----------------------------------------------------------
  ! upload_fine CUDA: sync
  !-----------------------------------------------------------
  subroutine upload_fine_cuda_sync_f(stream_slot)
    implicit none
    integer(c_int), intent(in) :: stream_slot
    call upload_fine_cuda_sync_c(stream_slot)
  end subroutine upload_fine_cuda_sync_f

  !-----------------------------------------------------------
  ! synchro_hydro CUDA: async launch
  !-----------------------------------------------------------
  subroutine synchro_cuda_async_f(buf, smallr_in, dteff, ncell, ndim_in, stride, stream_slot)
    implicit none
    real(c_double), intent(inout), target :: buf(*)
    real(c_double), intent(in) :: smallr_in, dteff
    integer(c_int), intent(in) :: ncell, ndim_in, stride, stream_slot
    call synchro_cuda_async_c(c_loc(buf(1)), smallr_in, dteff, ncell, ndim_in, stride, stream_slot)
  end subroutine synchro_cuda_async_f

  !-----------------------------------------------------------
  ! synchro_hydro CUDA: sync
  !-----------------------------------------------------------
  subroutine synchro_cuda_sync_f(stream_slot)
    implicit none
    integer(c_int), intent(in) :: stream_slot
    call synchro_cuda_sync_c(stream_slot)
  end subroutine synchro_cuda_sync_f

  !-----------------------------------------------------------
  ! cmpdt CUDA: async launch
  !-----------------------------------------------------------
  subroutine cmpdt_cuda_async_f(buf, dt_buf, dx, smallr_in, smallc_in, &
       gamma_in, courant_in, ncell, nvar_in, ndim_in, stride, stream_slot)
    implicit none
    real(c_double), intent(in), target :: buf(*)
    real(c_double), intent(out), target :: dt_buf(*)
    real(c_double), intent(in) :: dx, smallr_in, smallc_in, gamma_in, courant_in
    integer(c_int), intent(in) :: ncell, nvar_in, ndim_in, stride, stream_slot
    call cmpdt_cuda_async_c(c_loc(buf(1)), c_loc(dt_buf(1)), &
         dx, smallr_in, smallc_in, gamma_in, courant_in, &
         ncell, nvar_in, ndim_in, stride, stream_slot)
  end subroutine cmpdt_cuda_async_f

  !-----------------------------------------------------------
  ! cmpdt CUDA: sync
  !-----------------------------------------------------------
  subroutine cmpdt_cuda_sync_f(stream_slot)
    implicit none
    integer(c_int), intent(in) :: stream_slot
    call cmpdt_cuda_sync_c(stream_slot)
  end subroutine cmpdt_cuda_sync_f

end module hydro_cuda_interface
