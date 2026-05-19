! Minimal amr_parameters stub for file_pool_test (decouples ramses_hdf5_io
! from the rest of the RAMSES build).
module amr_parameters
  implicit none
  integer, parameter :: dp  = kind(1.0d0)
  integer, parameter :: i8b = selected_int_kind(15)
end module amr_parameters
