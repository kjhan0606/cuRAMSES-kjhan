!###########################################################
! Fortran interface to CUDA AGN feedback kernels
! Provides ISO_C_BINDING wrappers for sink_cuda_kernels.cu
!###########################################################
module sink_cuda_interface
  use iso_c_binding
  implicit none

  interface
     ! average_AGN on GPU: returns 1 on success, 0 on failure
     integer(c_int) function cuda_sink_average_agn_c( &
          cell_x, cell_y, cell_z, cell_vol, cell_d, cell_dx, cell_idx, ncells, &
          agn_x, agn_j, agn_dMBH, agn_dMEd, agn_Esave, agn_ok_blast, nAGN, &
          bin_head, agn_next, &
          nbin, inv_bin_size, rmax, rmax2, X_floor, &
          vol_gas, mass_gas, psy_norm, vol_blast, mass_blast, ind_blast) &
          bind(C, name='cuda_sink_average_agn')
       import :: c_double, c_int
       real(c_double) :: cell_x(*), cell_y(*), cell_z(*)
       real(c_double) :: cell_vol(*), cell_d(*), cell_dx(*)
       integer(c_int) :: cell_idx(*)
       integer(c_int), value :: ncells
       real(c_double) :: agn_x(*), agn_j(*)
       real(c_double) :: agn_dMBH(*), agn_dMEd(*), agn_Esave(*)
       integer(c_int) :: agn_ok_blast(*)
       integer(c_int), value :: nAGN
       integer(c_int) :: bin_head(*), agn_next(*)
       integer(c_int), value :: nbin
       real(c_double), value :: inv_bin_size, rmax, rmax2, X_floor
       real(c_double) :: vol_gas(*), mass_gas(*), psy_norm(*)
       real(c_double) :: vol_blast(*), mass_blast(*)
       integer(c_int) :: ind_blast(*)
     end function

     ! AGN_blast on GPU: returns 1 on success, 0 on failure
     integer(c_int) function cuda_sink_agn_blast_c( &
          cell_x, cell_y, cell_z, cell_vol, cell_dx, &
          cell_uold, ncells, nvar, &
          agn_x, agn_j, agn_v, &
          agn_p_gas, agn_uBlast, agn_mAGN, agn_ZAGN, &
          agn_psy_norm, agn_vol_gas, &
          agn_dMBH, agn_dMEd, &
          agn_ok_blast, agn_ok_save, nAGN, &
          bin_head, agn_next, &
          nbin, inv_bin_size, rmax, rmax2, X_floor, &
          scale_T2, T2maxAGNz, imetal, gamma_val, &
          Esave_out) &
          bind(C, name='cuda_sink_agn_blast')
       import :: c_double, c_int
       real(c_double) :: cell_x(*), cell_y(*), cell_z(*)
       real(c_double) :: cell_vol(*), cell_dx(*)
       real(c_double) :: cell_uold(*)
       integer(c_int), value :: ncells, nvar
       real(c_double) :: agn_x(*), agn_j(*), agn_v(*)
       real(c_double) :: agn_p_gas(*), agn_uBlast(*), agn_mAGN(*), agn_ZAGN(*)
       real(c_double) :: agn_psy_norm(*), agn_vol_gas(*)
       real(c_double) :: agn_dMBH(*), agn_dMEd(*)
       integer(c_int) :: agn_ok_blast(*), agn_ok_save(*)
       integer(c_int), value :: nAGN
       integer(c_int) :: bin_head(*), agn_next(*)
       integer(c_int), value :: nbin
       real(c_double), value :: inv_bin_size, rmax, rmax2, X_floor
       real(c_double), value :: scale_T2, T2maxAGNz
       integer(c_int), value :: imetal
       real(c_double), value :: gamma_val
       real(c_double) :: Esave_out(*)
     end function

     subroutine cuda_sink_finalize_c() bind(C, name='cuda_sink_finalize')
     end subroutine
  end interface

end module sink_cuda_interface
