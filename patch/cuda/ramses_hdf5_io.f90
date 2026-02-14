!###########################################################################
! ramses_hdf5_io.f90 — HDF5 wrapper module for cuRAMSES
!
! Provides helper routines for parallel HDF5 file create/open/close,
! attribute and dataset I/O.
! All routines are guarded by #ifdef HDF5.
!###########################################################################
#ifdef HDF5
module ramses_hdf5_io
  use amr_parameters, only: dp, i8b
  use hdf5
  implicit none

  integer(HID_T) :: hdf5_file_id    ! current file handle
  integer(HID_T) :: hdf5_plist_id   ! file access property list (MPI-IO)

contains

  !=========================================================================
  ! File-level operations
  !=========================================================================
  subroutine hdf5_create_parallel(filename, comm)
    implicit none
    include 'mpif.h'
    character(len=*), intent(in) :: filename
    integer, intent(in) :: comm
    integer :: ierr

    call h5open_f(ierr)
    ! Create file access property list for parallel I/O
    call h5pcreate_f(H5P_FILE_ACCESS_F, hdf5_plist_id, ierr)
    call h5pset_fapl_mpio_f(hdf5_plist_id, comm, MPI_INFO_NULL, ierr)
    ! Create file
    call h5fcreate_f(trim(filename), H5F_ACC_TRUNC_F, hdf5_file_id, ierr, &
         access_prp=hdf5_plist_id)
    call h5pclose_f(hdf5_plist_id, ierr)
  end subroutine hdf5_create_parallel

  subroutine hdf5_open_parallel(filename, comm)
    implicit none
    include 'mpif.h'
    character(len=*), intent(in) :: filename
    integer, intent(in) :: comm
    integer :: ierr

    call h5open_f(ierr)
    call h5pcreate_f(H5P_FILE_ACCESS_F, hdf5_plist_id, ierr)
    call h5pset_fapl_mpio_f(hdf5_plist_id, comm, MPI_INFO_NULL, ierr)
    call h5fopen_f(trim(filename), H5F_ACC_RDONLY_F, hdf5_file_id, ierr, &
         access_prp=hdf5_plist_id)
    call h5pclose_f(hdf5_plist_id, ierr)
  end subroutine hdf5_open_parallel

  subroutine hdf5_close_file()
    implicit none
    integer :: ierr
    call h5fclose_f(hdf5_file_id, ierr)
    call h5close_f(ierr)
  end subroutine hdf5_close_file

  !=========================================================================
  ! Group operations
  !=========================================================================
  subroutine hdf5_create_group(grp_name, grp_id)
    implicit none
    character(len=*), intent(in) :: grp_name
    integer(HID_T), intent(out) :: grp_id
    integer :: ierr
    call h5gcreate_f(hdf5_file_id, trim(grp_name), grp_id, ierr)
  end subroutine hdf5_create_group

  subroutine hdf5_open_group(grp_name, grp_id)
    implicit none
    character(len=*), intent(in) :: grp_name
    integer(HID_T), intent(out) :: grp_id
    integer :: ierr
    call h5gopen_f(hdf5_file_id, trim(grp_name), grp_id, ierr)
  end subroutine hdf5_open_group

  subroutine hdf5_close_group(grp_id)
    implicit none
    integer(HID_T), intent(in) :: grp_id
    integer :: ierr
    call h5gclose_f(grp_id, ierr)
  end subroutine hdf5_close_group

  !=========================================================================
  ! Attribute write helpers
  !=========================================================================
  subroutine hdf5_write_attr_int(loc_id, name, val)
    implicit none
    integer(HID_T), intent(in) :: loc_id
    character(len=*), intent(in) :: name
    integer, intent(in) :: val
    integer(HID_T) :: space_id, attr_id
    integer(HSIZE_T), dimension(1) :: dims = (/1/)
    integer :: ierr
    call h5screate_simple_f(1, dims, space_id, ierr)
    call h5acreate_f(loc_id, trim(name), H5T_NATIVE_INTEGER, space_id, attr_id, ierr)
    call h5awrite_f(attr_id, H5T_NATIVE_INTEGER, val, dims, ierr)
    call h5aclose_f(attr_id, ierr)
    call h5sclose_f(space_id, ierr)
  end subroutine hdf5_write_attr_int

  subroutine hdf5_write_attr_int8(loc_id, name, val)
    implicit none
    integer(HID_T), intent(in) :: loc_id
    character(len=*), intent(in) :: name
    integer(i8b), intent(in) :: val
    integer(HID_T) :: space_id, attr_id, type_id
    integer(HSIZE_T), dimension(1) :: dims = (/1/)
    integer :: ierr
    call h5screate_simple_f(1, dims, space_id, ierr)
    call h5tcopy_f(H5T_NATIVE_INTEGER, type_id, ierr)
    call h5tset_size_f(type_id, int(8, SIZE_T), ierr)
    call h5acreate_f(loc_id, trim(name), type_id, space_id, attr_id, ierr)
    call h5awrite_f(attr_id, type_id, val, dims, ierr)
    call h5aclose_f(attr_id, ierr)
    call h5tclose_f(type_id, ierr)
    call h5sclose_f(space_id, ierr)
  end subroutine hdf5_write_attr_int8

  subroutine hdf5_write_attr_dp(loc_id, name, val)
    implicit none
    integer(HID_T), intent(in) :: loc_id
    character(len=*), intent(in) :: name
    real(dp), intent(in) :: val
    integer(HID_T) :: space_id, attr_id
    integer(HSIZE_T), dimension(1) :: dims = (/1/)
    integer :: ierr
    call h5screate_simple_f(1, dims, space_id, ierr)
    call h5acreate_f(loc_id, trim(name), H5T_NATIVE_DOUBLE, space_id, attr_id, ierr)
    call h5awrite_f(attr_id, H5T_NATIVE_DOUBLE, val, dims, ierr)
    call h5aclose_f(attr_id, ierr)
    call h5sclose_f(space_id, ierr)
  end subroutine hdf5_write_attr_dp

  subroutine hdf5_write_attr_string(loc_id, name, val)
    implicit none
    integer(HID_T), intent(in) :: loc_id
    character(len=*), intent(in) :: name
    character(len=*), intent(in) :: val
    integer(HID_T) :: space_id, attr_id, type_id
    integer(HSIZE_T), dimension(1) :: dims = (/1/)
    integer(SIZE_T) :: slen
    integer :: ierr
    slen = len_trim(val)
    call h5screate_simple_f(1, dims, space_id, ierr)
    call h5tcopy_f(H5T_NATIVE_CHARACTER, type_id, ierr)
    call h5tset_size_f(type_id, slen, ierr)
    call h5acreate_f(loc_id, trim(name), type_id, space_id, attr_id, ierr)
    call h5awrite_f(attr_id, type_id, trim(val), dims, ierr)
    call h5aclose_f(attr_id, ierr)
    call h5tclose_f(type_id, ierr)
    call h5sclose_f(space_id, ierr)
  end subroutine hdf5_write_attr_string

  !=========================================================================
  ! Attribute read helpers
  !=========================================================================
  subroutine hdf5_read_attr_int(loc_id, name, val)
    implicit none
    integer(HID_T), intent(in) :: loc_id
    character(len=*), intent(in) :: name
    integer, intent(out) :: val
    integer(HID_T) :: attr_id
    integer(HSIZE_T), dimension(1) :: dims = (/1/)
    integer :: ierr
    call h5aopen_f(loc_id, trim(name), attr_id, ierr)
    call h5aread_f(attr_id, H5T_NATIVE_INTEGER, val, dims, ierr)
    call h5aclose_f(attr_id, ierr)
  end subroutine hdf5_read_attr_int

  subroutine hdf5_read_attr_int8(loc_id, name, val)
    implicit none
    integer(HID_T), intent(in) :: loc_id
    character(len=*), intent(in) :: name
    integer(i8b), intent(out) :: val
    integer(HID_T) :: attr_id, type_id
    integer(HSIZE_T), dimension(1) :: dims = (/1/)
    integer :: ierr
    call h5aopen_f(loc_id, trim(name), attr_id, ierr)
    call h5tcopy_f(H5T_NATIVE_INTEGER, type_id, ierr)
    call h5tset_size_f(type_id, int(8, SIZE_T), ierr)
    call h5aread_f(attr_id, type_id, val, dims, ierr)
    call h5aclose_f(attr_id, ierr)
    call h5tclose_f(type_id, ierr)
  end subroutine hdf5_read_attr_int8

  subroutine hdf5_read_attr_dp(loc_id, name, val)
    implicit none
    integer(HID_T), intent(in) :: loc_id
    character(len=*), intent(in) :: name
    real(dp), intent(out) :: val
    integer(HID_T) :: attr_id
    integer(HSIZE_T), dimension(1) :: dims = (/1/)
    integer :: ierr
    call h5aopen_f(loc_id, trim(name), attr_id, ierr)
    call h5aread_f(attr_id, H5T_NATIVE_DOUBLE, val, dims, ierr)
    call h5aclose_f(attr_id, ierr)
  end subroutine hdf5_read_attr_dp

  subroutine hdf5_read_attr_string(loc_id, name, val)
    implicit none
    integer(HID_T), intent(in) :: loc_id
    character(len=*), intent(in) :: name
    character(len=*), intent(out) :: val
    integer(HID_T) :: attr_id, type_id
    integer(HSIZE_T), dimension(1) :: dims = (/1/)
    integer(SIZE_T) :: slen
    integer :: ierr
    slen = len(val)
    call h5aopen_f(loc_id, trim(name), attr_id, ierr)
    call h5tcopy_f(H5T_NATIVE_CHARACTER, type_id, ierr)
    call h5tset_size_f(type_id, slen, ierr)
    val = ''
    call h5aread_f(attr_id, type_id, val, dims, ierr)
    call h5aclose_f(attr_id, ierr)
    call h5tclose_f(type_id, ierr)
  end subroutine hdf5_read_attr_string

  !=========================================================================
  ! 1D array attribute write/read
  !=========================================================================
  subroutine hdf5_write_attr_1d_dp(loc_id, name, arr, n)
    implicit none
    integer(HID_T), intent(in) :: loc_id
    character(len=*), intent(in) :: name
    integer, intent(in) :: n
    real(dp), intent(in) :: arr(n)
    integer(HID_T) :: space_id, attr_id
    integer(HSIZE_T), dimension(1) :: dims
    integer :: ierr
    dims(1) = n
    call h5screate_simple_f(1, dims, space_id, ierr)
    call h5acreate_f(loc_id, trim(name), H5T_NATIVE_DOUBLE, space_id, attr_id, ierr)
    call h5awrite_f(attr_id, H5T_NATIVE_DOUBLE, arr, dims, ierr)
    call h5aclose_f(attr_id, ierr)
    call h5sclose_f(space_id, ierr)
  end subroutine hdf5_write_attr_1d_dp

  subroutine hdf5_read_attr_1d_dp(loc_id, name, arr, n)
    implicit none
    integer(HID_T), intent(in) :: loc_id
    character(len=*), intent(in) :: name
    integer, intent(in) :: n
    real(dp), intent(out) :: arr(n)
    integer(HID_T) :: attr_id
    integer(HSIZE_T), dimension(1) :: dims
    integer :: ierr
    dims(1) = n
    call h5aopen_f(loc_id, trim(name), attr_id, ierr)
    call h5aread_f(attr_id, H5T_NATIVE_DOUBLE, arr, dims, ierr)
    call h5aclose_f(attr_id, ierr)
  end subroutine hdf5_read_attr_1d_dp

  !=========================================================================
  ! Dataset write: all ranks collectively write with hyperslab selection
  !=========================================================================
  subroutine hdf5_write_dataset_1d_dp(grp_id, name, data, nlocal, offset_global, ntotal)
    ! Write a 1D double-precision dataset with hyperslab parallel write
    implicit none
    include 'mpif.h'
    integer(HID_T), intent(in) :: grp_id
    character(len=*), intent(in) :: name
    integer, intent(in) :: nlocal
    integer(i8b), intent(in) :: offset_global, ntotal
    real(dp), intent(in) :: data(nlocal)
    integer(HID_T) :: dspace_id, dset_id, memspace_id, plist_id
    integer(HSIZE_T), dimension(1) :: gdims, ldims, offset, count
    integer :: ierr

    gdims(1) = ntotal
    call h5screate_simple_f(1, gdims, dspace_id, ierr)
    call h5dcreate_f(grp_id, trim(name), H5T_NATIVE_DOUBLE, dspace_id, dset_id, ierr)

    ! Select hyperslab in file space
    offset(1) = offset_global
    count(1) = nlocal
    call h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F, offset, count, ierr)

    ! Create memory space
    ldims(1) = nlocal
    call h5screate_simple_f(1, ldims, memspace_id, ierr)

    ! Collective write
    call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, ierr)
    call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, ierr)
    call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, data, ldims, ierr, &
         mem_space_id=memspace_id, file_space_id=dspace_id, xfer_prp=plist_id)

    call h5pclose_f(plist_id, ierr)
    call h5sclose_f(memspace_id, ierr)
    call h5dclose_f(dset_id, ierr)
    call h5sclose_f(dspace_id, ierr)
  end subroutine hdf5_write_dataset_1d_dp

  subroutine hdf5_write_dataset_1d_int(grp_id, name, data, nlocal, offset_global, ntotal)
    implicit none
    include 'mpif.h'
    integer(HID_T), intent(in) :: grp_id
    character(len=*), intent(in) :: name
    integer, intent(in) :: nlocal
    integer(i8b), intent(in) :: offset_global, ntotal
    integer, intent(in) :: data(nlocal)
    integer(HID_T) :: dspace_id, dset_id, memspace_id, plist_id
    integer(HSIZE_T), dimension(1) :: gdims, ldims, offset, count
    integer :: ierr

    gdims(1) = ntotal
    call h5screate_simple_f(1, gdims, dspace_id, ierr)
    call h5dcreate_f(grp_id, trim(name), H5T_NATIVE_INTEGER, dspace_id, dset_id, ierr)

    offset(1) = offset_global
    count(1) = nlocal
    call h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F, offset, count, ierr)

    ldims(1) = nlocal
    call h5screate_simple_f(1, ldims, memspace_id, ierr)

    call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, ierr)
    call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, ierr)
    call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, data, ldims, ierr, &
         mem_space_id=memspace_id, file_space_id=dspace_id, xfer_prp=plist_id)

    call h5pclose_f(plist_id, ierr)
    call h5sclose_f(memspace_id, ierr)
    call h5dclose_f(dset_id, ierr)
    call h5sclose_f(dspace_id, ierr)
  end subroutine hdf5_write_dataset_1d_int

  subroutine hdf5_write_dataset_1d_int8(grp_id, name, data, nlocal, offset_global, ntotal)
    implicit none
    include 'mpif.h'
    integer(HID_T), intent(in) :: grp_id
    character(len=*), intent(in) :: name
    integer, intent(in) :: nlocal
    integer(i8b), intent(in) :: offset_global, ntotal
    integer(i8b), intent(in) :: data(nlocal)
    integer(HID_T) :: dspace_id, dset_id, memspace_id, plist_id, type_id
    integer(HSIZE_T), dimension(1) :: gdims, ldims, offset, count
    integer :: ierr

    gdims(1) = ntotal
    call h5screate_simple_f(1, gdims, dspace_id, ierr)
    call h5tcopy_f(H5T_NATIVE_INTEGER, type_id, ierr)
    call h5tset_size_f(type_id, int(8, SIZE_T), ierr)
    call h5dcreate_f(grp_id, trim(name), type_id, dspace_id, dset_id, ierr)

    offset(1) = offset_global
    count(1) = nlocal
    call h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F, offset, count, ierr)

    ldims(1) = nlocal
    call h5screate_simple_f(1, ldims, memspace_id, ierr)

    call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, ierr)
    call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, ierr)
    call h5dwrite_f(dset_id, type_id, data, ldims, ierr, &
         mem_space_id=memspace_id, file_space_id=dspace_id, xfer_prp=plist_id)

    call h5pclose_f(plist_id, ierr)
    call h5tclose_f(type_id, ierr)
    call h5sclose_f(memspace_id, ierr)
    call h5dclose_f(dset_id, ierr)
    call h5sclose_f(dspace_id, ierr)
  end subroutine hdf5_write_dataset_1d_int8

  !=========================================================================
  ! Dataset read: all ranks collectively read with hyperslab selection
  !=========================================================================
  subroutine hdf5_read_dataset_1d_dp(grp_id, name, data, nlocal, offset_global)
    implicit none
    include 'mpif.h'
    integer(HID_T), intent(in) :: grp_id
    character(len=*), intent(in) :: name
    integer, intent(in) :: nlocal
    integer(i8b), intent(in) :: offset_global
    real(dp), intent(out) :: data(nlocal)
    integer(HID_T) :: dset_id, dspace_id, memspace_id, plist_id
    integer(HSIZE_T), dimension(1) :: ldims, offset, count
    integer :: ierr

    call h5dopen_f(grp_id, trim(name), dset_id, ierr)
    call h5dget_space_f(dset_id, dspace_id, ierr)

    offset(1) = offset_global
    count(1) = nlocal
    call h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F, offset, count, ierr)

    ldims(1) = nlocal
    call h5screate_simple_f(1, ldims, memspace_id, ierr)

    call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, ierr)
    call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, ierr)
    call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, data, ldims, ierr, &
         mem_space_id=memspace_id, file_space_id=dspace_id, xfer_prp=plist_id)

    call h5pclose_f(plist_id, ierr)
    call h5sclose_f(memspace_id, ierr)
    call h5sclose_f(dspace_id, ierr)
    call h5dclose_f(dset_id, ierr)
  end subroutine hdf5_read_dataset_1d_dp

  subroutine hdf5_read_dataset_1d_int(grp_id, name, data, nlocal, offset_global)
    implicit none
    include 'mpif.h'
    integer(HID_T), intent(in) :: grp_id
    character(len=*), intent(in) :: name
    integer, intent(in) :: nlocal
    integer(i8b), intent(in) :: offset_global
    integer, intent(out) :: data(nlocal)
    integer(HID_T) :: dset_id, dspace_id, memspace_id, plist_id
    integer(HSIZE_T), dimension(1) :: ldims, offset, count
    integer :: ierr

    call h5dopen_f(grp_id, trim(name), dset_id, ierr)
    call h5dget_space_f(dset_id, dspace_id, ierr)

    offset(1) = offset_global
    count(1) = nlocal
    call h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F, offset, count, ierr)

    ldims(1) = nlocal
    call h5screate_simple_f(1, ldims, memspace_id, ierr)

    call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, ierr)
    call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, ierr)
    call h5dread_f(dset_id, H5T_NATIVE_INTEGER, data, ldims, ierr, &
         mem_space_id=memspace_id, file_space_id=dspace_id, xfer_prp=plist_id)

    call h5pclose_f(plist_id, ierr)
    call h5sclose_f(memspace_id, ierr)
    call h5sclose_f(dspace_id, ierr)
    call h5dclose_f(dset_id, ierr)
  end subroutine hdf5_read_dataset_1d_int

  subroutine hdf5_read_dataset_1d_int8(grp_id, name, data, nlocal, offset_global)
    implicit none
    include 'mpif.h'
    integer(HID_T), intent(in) :: grp_id
    character(len=*), intent(in) :: name
    integer, intent(in) :: nlocal
    integer(i8b), intent(in) :: offset_global
    integer(i8b), intent(out) :: data(nlocal)
    integer(HID_T) :: dset_id, dspace_id, memspace_id, plist_id, type_id
    integer(HSIZE_T), dimension(1) :: ldims, offset, count
    integer :: ierr

    call h5dopen_f(grp_id, trim(name), dset_id, ierr)
    call h5dget_space_f(dset_id, dspace_id, ierr)

    offset(1) = offset_global
    count(1) = nlocal
    call h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F, offset, count, ierr)

    ldims(1) = nlocal
    call h5screate_simple_f(1, ldims, memspace_id, ierr)

    call h5tcopy_f(H5T_NATIVE_INTEGER, type_id, ierr)
    call h5tset_size_f(type_id, int(8, SIZE_T), ierr)

    call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, ierr)
    call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, ierr)
    call h5dread_f(dset_id, type_id, data, ldims, ierr, &
         mem_space_id=memspace_id, file_space_id=dspace_id, xfer_prp=plist_id)

    call h5pclose_f(plist_id, ierr)
    call h5tclose_f(type_id, ierr)
    call h5sclose_f(memspace_id, ierr)
    call h5sclose_f(dspace_id, ierr)
    call h5dclose_f(dset_id, ierr)
  end subroutine hdf5_read_dataset_1d_int8

  !=========================================================================
  ! Rank-0-only dataset write/read (for sinks, small metadata)
  !=========================================================================
  subroutine hdf5_write_dataset_serial_dp(grp_id, name, data, n, myid)
    implicit none
    integer(HID_T), intent(in) :: grp_id
    character(len=*), intent(in) :: name
    integer, intent(in) :: n, myid
    real(dp), intent(in) :: data(n)
    integer(HID_T) :: dspace_id, dset_id, memspace_id, plist_id
    integer(HSIZE_T), dimension(1) :: gdims, ldims, offset, count
    integer :: ierr

    gdims(1) = n
    call h5screate_simple_f(1, gdims, dspace_id, ierr)
    call h5dcreate_f(grp_id, trim(name), H5T_NATIVE_DOUBLE, dspace_id, dset_id, ierr)

    ! Only rank 1 (myid==1) writes; others select none
    if(myid == 1) then
       offset(1) = 0; count(1) = n
    else
       offset(1) = 0; count(1) = 0
    end if
    call h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F, offset, count, ierr)

    ldims(1) = count(1)
    call h5screate_simple_f(1, ldims, memspace_id, ierr)
    if(myid /= 1) call h5sselect_none_f(memspace_id, ierr)

    call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, ierr)
    call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, ierr)
    call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, data, ldims, ierr, &
         mem_space_id=memspace_id, file_space_id=dspace_id, xfer_prp=plist_id)

    call h5pclose_f(plist_id, ierr)
    call h5sclose_f(memspace_id, ierr)
    call h5dclose_f(dset_id, ierr)
    call h5sclose_f(dspace_id, ierr)
  end subroutine hdf5_write_dataset_serial_dp

  subroutine hdf5_write_dataset_serial_int(grp_id, name, data, n, myid)
    implicit none
    integer(HID_T), intent(in) :: grp_id
    character(len=*), intent(in) :: name
    integer, intent(in) :: n, myid
    integer, intent(in) :: data(n)
    integer(HID_T) :: dspace_id, dset_id, memspace_id, plist_id
    integer(HSIZE_T), dimension(1) :: gdims, ldims, offset, count
    integer :: ierr

    gdims(1) = n
    call h5screate_simple_f(1, gdims, dspace_id, ierr)
    call h5dcreate_f(grp_id, trim(name), H5T_NATIVE_INTEGER, dspace_id, dset_id, ierr)

    if(myid == 1) then
       offset(1) = 0; count(1) = n
    else
       offset(1) = 0; count(1) = 0
    end if
    call h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F, offset, count, ierr)

    ldims(1) = count(1)
    call h5screate_simple_f(1, ldims, memspace_id, ierr)
    if(myid /= 1) call h5sselect_none_f(memspace_id, ierr)

    call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, ierr)
    call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, ierr)
    call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, data, ldims, ierr, &
         mem_space_id=memspace_id, file_space_id=dspace_id, xfer_prp=plist_id)

    call h5pclose_f(plist_id, ierr)
    call h5sclose_f(memspace_id, ierr)
    call h5dclose_f(dset_id, ierr)
    call h5sclose_f(dspace_id, ierr)
  end subroutine hdf5_write_dataset_serial_int

  subroutine hdf5_read_dataset_all_dp(grp_id, name, data, n)
    ! All ranks read the entire dataset (for small arrays like sinks)
    implicit none
    include 'mpif.h'
    integer(HID_T), intent(in) :: grp_id
    character(len=*), intent(in) :: name
    integer, intent(in) :: n
    real(dp), intent(out) :: data(n)
    integer(HID_T) :: dset_id, plist_id
    integer(HSIZE_T), dimension(1) :: dims
    integer :: ierr

    dims(1) = n
    call h5dopen_f(grp_id, trim(name), dset_id, ierr)
    call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, ierr)
    call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, ierr)
    call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, data, dims, ierr, xfer_prp=plist_id)
    call h5pclose_f(plist_id, ierr)
    call h5dclose_f(dset_id, ierr)
  end subroutine hdf5_read_dataset_all_dp

  subroutine hdf5_read_dataset_all_int(grp_id, name, data, n)
    implicit none
    include 'mpif.h'
    integer(HID_T), intent(in) :: grp_id
    character(len=*), intent(in) :: name
    integer, intent(in) :: n
    integer, intent(out) :: data(n)
    integer(HID_T) :: dset_id, plist_id
    integer(HSIZE_T), dimension(1) :: dims
    integer :: ierr

    dims(1) = n
    call h5dopen_f(grp_id, trim(name), dset_id, ierr)
    call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, ierr)
    call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, ierr)
    call h5dread_f(dset_id, H5T_NATIVE_INTEGER, data, dims, ierr, xfer_prp=plist_id)
    call h5pclose_f(plist_id, ierr)
    call h5dclose_f(dset_id, ierr)
  end subroutine hdf5_read_dataset_all_int

  subroutine hdf5_suppress_errors()
    implicit none
    integer :: ierr
    call h5eset_auto_f(0, ierr)
  end subroutine hdf5_suppress_errors

  subroutine hdf5_restore_errors()
    implicit none
    integer :: ierr
    call h5eset_auto_f(1, ierr)
  end subroutine hdf5_restore_errors

end module ramses_hdf5_io
#endif
