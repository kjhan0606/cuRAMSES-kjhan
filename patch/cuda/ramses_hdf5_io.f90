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

  ! ----- Phase D: file pool (LRU cache of open HDF5 file handles) -----
  ! Future-proofing for multi-file outputs (one file per writer-group).
  ! Single-shared-file callers (current production path) bypass the pool
  ! entirely; nothing changes for them. Callers that want pooling use
  !   call file_pool_init(max_open)
  !   call file_pool_get(filename, comm, file_id)
  !   ... (handle is owned by the pool; do not close it directly)
  !   call file_pool_close_all()
  integer, parameter :: FP_NAME_LEN    = 512
  integer, parameter :: FP_MAX_DEFAULT = 8
  integer, parameter :: FP_MAX_LIMIT   = 1024
  integer                    :: fp_cap   = 0
  integer                    :: fp_count = 0
  integer(i8b)               :: fp_tick  = 0
  character(len=FP_NAME_LEN), allocatable :: fp_name(:)
  integer(HID_T),             allocatable :: fp_file_id(:)
  integer(i8b),               allocatable :: fp_lru(:)
  ! Diagnostic counters
  integer(i8b) :: fp_hits = 0, fp_misses = 0, fp_evictions = 0

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
    ! Replicate dataset to ALL ranks via rank-0 streamed read + MPI_Bcast.
    ! Was: every rank collectively re-read the full dataset → O(N²) GPFS traffic.
    ! Now: rank 0 reads n/nproc-sized chunks sequentially, broadcasting each chunk
    !      so the GPFS volume sees ~1× n bytes total instead of nproc × n bytes.
    implicit none
    include 'mpif.h'
    integer(HID_T), intent(in) :: grp_id
    character(len=*), intent(in) :: name
    integer, intent(in) :: n
    real(dp), intent(out) :: data(n)
    integer(HID_T) :: dset_id, dspace_id, memspace_id, plist_id
    integer(HSIZE_T), dimension(1) :: ldims, fooffset, count
    integer :: ierr, myrank, nproc, ichunk, nchunk, chunk_size, this_chunk, istart

    call MPI_Comm_rank(MPI_COMM_WORLD, myrank, ierr)
    call MPI_Comm_size(MPI_COMM_WORLD, nproc, ierr)

    if (n <= 0) return

    if (n < nproc) then
       chunk_size = n
    else
       chunk_size = (n + nproc - 1) / nproc
    end if
    nchunk = (n + chunk_size - 1) / chunk_size

    call h5dopen_f(grp_id, trim(name), dset_id, ierr)
    if (myrank == 0) call h5dget_space_f(dset_id, dspace_id, ierr)
    call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, ierr)
    call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_INDEPENDENT_F, ierr)

    do ichunk = 1, nchunk
       istart = (ichunk-1) * chunk_size + 1
       this_chunk = min(chunk_size, n - istart + 1)
       if (this_chunk <= 0) exit
       if (myrank == 0) then
          fooffset(1) = int(istart - 1, HSIZE_T)
          count(1) = int(this_chunk, HSIZE_T)
          call h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F, fooffset, count, ierr)
          ldims(1) = int(this_chunk, HSIZE_T)
          call h5screate_simple_f(1, ldims, memspace_id, ierr)
          call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, data(istart:istart+this_chunk-1), &
               ldims, ierr, mem_space_id=memspace_id, file_space_id=dspace_id, &
               xfer_prp=plist_id)
          call h5sclose_f(memspace_id, ierr)
       end if
       call MPI_Bcast(data(istart), this_chunk, MPI_DOUBLE_PRECISION, 0, &
            MPI_COMM_WORLD, ierr)
    end do

    call h5pclose_f(plist_id, ierr)
    if (myrank == 0) call h5sclose_f(dspace_id, ierr)
    call h5dclose_f(dset_id, ierr)
  end subroutine hdf5_read_dataset_all_dp

  subroutine hdf5_read_dataset_all_int(grp_id, name, data, n)
    ! Integer counterpart of hdf5_read_dataset_all_dp (rank-0 chunk read + Bcast).
    implicit none
    include 'mpif.h'
    integer(HID_T), intent(in) :: grp_id
    character(len=*), intent(in) :: name
    integer, intent(in) :: n
    integer, intent(out) :: data(n)
    integer(HID_T) :: dset_id, dspace_id, memspace_id, plist_id
    integer(HSIZE_T), dimension(1) :: ldims, fooffset, count
    integer :: ierr, myrank, nproc, ichunk, nchunk, chunk_size, this_chunk, istart

    call MPI_Comm_rank(MPI_COMM_WORLD, myrank, ierr)
    call MPI_Comm_size(MPI_COMM_WORLD, nproc, ierr)

    if (n <= 0) return

    if (n < nproc) then
       chunk_size = n
    else
       chunk_size = (n + nproc - 1) / nproc
    end if
    nchunk = (n + chunk_size - 1) / chunk_size

    call h5dopen_f(grp_id, trim(name), dset_id, ierr)
    if (myrank == 0) call h5dget_space_f(dset_id, dspace_id, ierr)
    call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, ierr)
    call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_INDEPENDENT_F, ierr)

    do ichunk = 1, nchunk
       istart = (ichunk-1) * chunk_size + 1
       this_chunk = min(chunk_size, n - istart + 1)
       if (this_chunk <= 0) exit
       if (myrank == 0) then
          fooffset(1) = int(istart - 1, HSIZE_T)
          count(1) = int(this_chunk, HSIZE_T)
          call h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F, fooffset, count, ierr)
          ldims(1) = int(this_chunk, HSIZE_T)
          call h5screate_simple_f(1, ldims, memspace_id, ierr)
          call h5dread_f(dset_id, H5T_NATIVE_INTEGER, data(istart:istart+this_chunk-1), &
               ldims, ierr, mem_space_id=memspace_id, file_space_id=dspace_id, &
               xfer_prp=plist_id)
          call h5sclose_f(memspace_id, ierr)
       end if
       call MPI_Bcast(data(istart), this_chunk, MPI_INTEGER, 0, &
            MPI_COMM_WORLD, ierr)
    end do

    call h5pclose_f(plist_id, ierr)
    if (myrank == 0) call h5sclose_f(dspace_id, ierr)
    call h5dclose_f(dset_id, ierr)
  end subroutine hdf5_read_dataset_all_int

  subroutine hdf5_read_dataset_chunk_dp(grp_id, name, data, n, offset_global)
    ! Replicate a chunk [offset_global+1 .. offset_global+n] of dataset to all
    ! ranks via rank-0 read + MPI_Bcast. Used by streaming varcpu restore to
    ! cap per-rank transient buffers (vs. full-N allocation in _all variants).
    implicit none
    include 'mpif.h'
    integer(HID_T), intent(in) :: grp_id
    character(len=*), intent(in) :: name
    integer, intent(in) :: n
    integer(i8b), intent(in) :: offset_global
    real(dp), intent(out) :: data(n)
    integer(HID_T) :: dset_id, dspace_id, memspace_id, plist_id
    integer(HSIZE_T), dimension(1) :: ldims, fooffset, count
    integer :: ierr, myrank

    call MPI_Comm_rank(MPI_COMM_WORLD, myrank, ierr)
    if (n <= 0) return

    call h5dopen_f(grp_id, trim(name), dset_id, ierr)
    if (myrank == 0) then
       call h5dget_space_f(dset_id, dspace_id, ierr)
       fooffset(1) = int(offset_global, HSIZE_T)
       count(1) = int(n, HSIZE_T)
       call h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F, fooffset, count, ierr)
       ldims(1) = int(n, HSIZE_T)
       call h5screate_simple_f(1, ldims, memspace_id, ierr)
       call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, ierr)
       call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_INDEPENDENT_F, ierr)
       call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, data, ldims, ierr, &
            mem_space_id=memspace_id, file_space_id=dspace_id, xfer_prp=plist_id)
       call h5pclose_f(plist_id, ierr)
       call h5sclose_f(memspace_id, ierr)
       call h5sclose_f(dspace_id, ierr)
    end if
    call h5dclose_f(dset_id, ierr)
    call MPI_Bcast(data, n, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  end subroutine hdf5_read_dataset_chunk_dp

  subroutine hdf5_read_dataset_chunk_int(grp_id, name, data, n, offset_global)
    ! Integer counterpart of hdf5_read_dataset_chunk_dp.
    implicit none
    include 'mpif.h'
    integer(HID_T), intent(in) :: grp_id
    character(len=*), intent(in) :: name
    integer, intent(in) :: n
    integer(i8b), intent(in) :: offset_global
    integer, intent(out) :: data(n)
    integer(HID_T) :: dset_id, dspace_id, memspace_id, plist_id
    integer(HSIZE_T), dimension(1) :: ldims, fooffset, count
    integer :: ierr, myrank

    call MPI_Comm_rank(MPI_COMM_WORLD, myrank, ierr)
    if (n <= 0) return

    call h5dopen_f(grp_id, trim(name), dset_id, ierr)
    if (myrank == 0) then
       call h5dget_space_f(dset_id, dspace_id, ierr)
       fooffset(1) = int(offset_global, HSIZE_T)
       count(1) = int(n, HSIZE_T)
       call h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F, fooffset, count, ierr)
       ldims(1) = int(n, HSIZE_T)
       call h5screate_simple_f(1, ldims, memspace_id, ierr)
       call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, ierr)
       call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_INDEPENDENT_F, ierr)
       call h5dread_f(dset_id, H5T_NATIVE_INTEGER, data, ldims, ierr, &
            mem_space_id=memspace_id, file_space_id=dspace_id, xfer_prp=plist_id)
       call h5pclose_f(plist_id, ierr)
       call h5sclose_f(memspace_id, ierr)
       call h5sclose_f(dspace_id, ierr)
    end if
    call h5dclose_f(dset_id, ierr)
    call MPI_Bcast(data, n, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  end subroutine hdf5_read_dataset_chunk_int

  subroutine hdf5_read_dataset_collective_dp(grp_id, name, data, n, offset_global)
    ! Tier 2 parallel restore primitive: each rank reads its own disjoint
    ! hyperslab [offset_global+1 .. offset_global+n] in a single collective
    ! H5Dread. No post-read Bcast: each rank gets only its slice.
    ! Ranks with n==0 must still participate (H5Sselect_none) — collective.
    implicit none
    include 'mpif.h'
    integer(HID_T), intent(in) :: grp_id
    character(len=*), intent(in) :: name
    integer, intent(in) :: n
    integer(i8b), intent(in) :: offset_global
    real(dp), intent(out) :: data(*)
    integer(HID_T) :: dset_id, dspace_id, memspace_id, plist_id
    integer(HSIZE_T), dimension(1) :: ldims, fooffset, count
    integer :: ierr

    call h5dopen_f(grp_id, trim(name), dset_id, ierr)
    call h5dget_space_f(dset_id, dspace_id, ierr)
    if (n > 0) then
       fooffset(1) = int(offset_global, HSIZE_T)
       count(1) = int(n, HSIZE_T)
       call h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F, fooffset, count, ierr)
       ldims(1) = int(n, HSIZE_T)
       call h5screate_simple_f(1, ldims, memspace_id, ierr)
    else
       call h5sselect_none_f(dspace_id, ierr)
       ldims(1) = 0_HSIZE_T
       call h5screate_simple_f(1, ldims, memspace_id, ierr)
       call h5sselect_none_f(memspace_id, ierr)
    end if
    call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, ierr)
    call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, ierr)
    call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, data, ldims, ierr, &
         mem_space_id=memspace_id, file_space_id=dspace_id, xfer_prp=plist_id)
    call h5pclose_f(plist_id, ierr)
    call h5sclose_f(memspace_id, ierr)
    call h5sclose_f(dspace_id, ierr)
    call h5dclose_f(dset_id, ierr)
  end subroutine hdf5_read_dataset_collective_dp

  subroutine hdf5_read_dataset_collective_int(grp_id, name, data, n, offset_global)
    implicit none
    include 'mpif.h'
    integer(HID_T), intent(in) :: grp_id
    character(len=*), intent(in) :: name
    integer, intent(in) :: n
    integer(i8b), intent(in) :: offset_global
    integer, intent(out) :: data(*)
    integer(HID_T) :: dset_id, dspace_id, memspace_id, plist_id
    integer(HSIZE_T), dimension(1) :: ldims, fooffset, count
    integer :: ierr

    call h5dopen_f(grp_id, trim(name), dset_id, ierr)
    call h5dget_space_f(dset_id, dspace_id, ierr)
    if (n > 0) then
       fooffset(1) = int(offset_global, HSIZE_T)
       count(1) = int(n, HSIZE_T)
       call h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F, fooffset, count, ierr)
       ldims(1) = int(n, HSIZE_T)
       call h5screate_simple_f(1, ldims, memspace_id, ierr)
    else
       call h5sselect_none_f(dspace_id, ierr)
       ldims(1) = 0_HSIZE_T
       call h5screate_simple_f(1, ldims, memspace_id, ierr)
       call h5sselect_none_f(memspace_id, ierr)
    end if
    call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, ierr)
    call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, ierr)
    call h5dread_f(dset_id, H5T_NATIVE_INTEGER, data, ldims, ierr, &
         mem_space_id=memspace_id, file_space_id=dspace_id, xfer_prp=plist_id)
    call h5pclose_f(plist_id, ierr)
    call h5sclose_f(memspace_id, ierr)
    call h5sclose_f(dspace_id, ierr)
    call h5dclose_f(dset_id, ierr)
  end subroutine hdf5_read_dataset_collective_int

  !=========================================================================
  ! Phase D: file pool (LRU cache of open HDF5 file handles)
  !
  ! Rationale: at exascale, restore may need to read from O(10^3-10^5)
  ! per-writer files. Keeping all of them open exhausts file descriptors;
  ! open/close per access wastes O(ms) HDF5 setup per dataset. An LRU pool
  ! of size K caps open descriptors at K while amortising open overhead
  ! over many dataset accesses to the same file.
  !
  ! Implementation: parallel array of (name, file_id, last_use_tick).
  ! get() does linear search (K small, typically 8-64) → returns cached
  ! file_id on hit; on miss opens with parallel MPI-IO, evicting the
  ! lowest-tick slot if full. Tick is monotonically incremented on every
  ! get() to define LRU order.
  !=========================================================================
  subroutine file_pool_init(max_open)
    implicit none
    integer, intent(in), optional :: max_open
    integer :: cap, ierr
    if(present(max_open)) then
       cap = max(1, min(max_open, FP_MAX_LIMIT))
    else
       cap = FP_MAX_DEFAULT
    end if
    ! Idempotent: re-init flushes whatever was there
    if(allocated(fp_file_id)) call file_pool_close_all()
    fp_cap = cap
    allocate(fp_name(cap), fp_file_id(cap), fp_lru(cap))
    fp_name = ''
    fp_file_id = -1_HID_T
    fp_lru = 0_i8b
    fp_count = 0
    fp_tick = 0_i8b
    fp_hits = 0_i8b
    fp_misses = 0_i8b
    fp_evictions = 0_i8b
    call h5open_f(ierr)
  end subroutine file_pool_init

  subroutine file_pool_get(filename, comm, file_id)
    implicit none
    include 'mpif.h'
    character(len=*), intent(in)  :: filename
    integer,          intent(in)  :: comm
    integer(HID_T),   intent(out) :: file_id
    integer :: i, victim, ierr
    integer(HID_T) :: plist_id

    ! Auto-init if caller forgot (uses default capacity)
    if(.not. allocated(fp_file_id)) call file_pool_init()

    fp_tick = fp_tick + 1_i8b

    ! Cache lookup
    do i = 1, fp_count
       if(trim(fp_name(i)) == trim(filename)) then
          file_id = fp_file_id(i)
          fp_lru(i) = fp_tick
          fp_hits = fp_hits + 1_i8b
          return
       end if
    end do

    fp_misses = fp_misses + 1_i8b

    ! Evict LRU slot if pool is full
    if(fp_count >= fp_cap) then
       victim = 1
       do i = 2, fp_count
          if(fp_lru(i) < fp_lru(victim)) victim = i
       end do
       call h5fclose_f(fp_file_id(victim), ierr)
       fp_evictions = fp_evictions + 1_i8b
       ! Compact: move the last live slot into the victim's hole
       if(victim /= fp_count) then
          fp_name(victim)    = fp_name(fp_count)
          fp_file_id(victim) = fp_file_id(fp_count)
          fp_lru(victim)     = fp_lru(fp_count)
       end if
       fp_name(fp_count)    = ''
       fp_file_id(fp_count) = -1_HID_T
       fp_lru(fp_count)     = 0_i8b
       fp_count = fp_count - 1
    end if

    ! Open new
    call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, ierr)
    call h5pset_fapl_mpio_f(plist_id, comm, MPI_INFO_NULL, ierr)
    call h5fopen_f(trim(filename), H5F_ACC_RDONLY_F, file_id, ierr, &
         access_prp=plist_id)
    call h5pclose_f(plist_id, ierr)

    fp_count = fp_count + 1
    fp_name(fp_count)    = trim(filename)
    fp_file_id(fp_count) = file_id
    fp_lru(fp_count)     = fp_tick
  end subroutine file_pool_get

  subroutine file_pool_close_all()
    implicit none
    integer :: i, ierr
    if(.not. allocated(fp_file_id)) return
    do i = 1, fp_count
       call h5fclose_f(fp_file_id(i), ierr)
    end do
    deallocate(fp_name, fp_file_id, fp_lru)
    fp_cap = 0
    fp_count = 0
    fp_tick = 0_i8b
  end subroutine file_pool_close_all

  subroutine file_pool_stats(n_open, n_hits, n_misses, n_evictions, cap)
    implicit none
    integer,      intent(out) :: n_open, cap
    integer(i8b), intent(out) :: n_hits, n_misses, n_evictions
    n_open      = fp_count
    cap         = fp_cap
    n_hits      = fp_hits
    n_misses    = fp_misses
    n_evictions = fp_evictions
  end subroutine file_pool_stats

  ! Convenience wrapper: open via pool AND publish handle into module's
  ! hdf5_file_id global so callers using the single-file API keep working.
  ! Pair with hdf5_release_pooled() (does NOT close — pool owns the handle).
  subroutine hdf5_open_pooled(filename, comm)
    implicit none
    character(len=*), intent(in) :: filename
    integer,          intent(in) :: comm
    integer(HID_T) :: fid
    call file_pool_get(filename, comm, fid)
    hdf5_file_id = fid
  end subroutine hdf5_open_pooled

  subroutine hdf5_release_pooled()
    implicit none
    hdf5_file_id = -1_HID_T
  end subroutine hdf5_release_pooled

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
