!###########################################################
! CUDA commons module for cuRAMSES
! Provides Fortran interface to the unified CUDA stream pool
!###########################################################
module cuda_commons
  use iso_c_binding
  implicit none

  logical :: cuda_available = .false.
  integer :: cuda_n_streams = 0

  ! C interface to stream pool
  interface
     subroutine cuda_pool_init_c(local_rank) bind(C, name='cuda_pool_init')
       import :: c_int
       integer(c_int), value :: local_rank
     end subroutine

     function cuda_acquire_stream_c() result(slot) bind(C, name='cuda_acquire_stream')
       import :: c_int
       integer(c_int) :: slot
     end function

     subroutine cuda_release_stream_c(slot) bind(C, name='cuda_release_stream')
       import :: c_int
       integer(c_int), value :: slot
     end subroutine

     subroutine cuda_stream_sync_c(slot) bind(C, name='cuda_stream_sync')
       import :: c_int
       integer(c_int), value :: slot
     end subroutine

     function cuda_stream_query_c(slot) result(ready) bind(C, name='cuda_stream_query')
       import :: c_int
       integer(c_int), value :: slot
       integer(c_int) :: ready
     end function

     function cuda_get_n_streams_c() result(n) bind(C, name='cuda_get_n_streams')
       import :: c_int
       integer(c_int) :: n
     end function

     subroutine cuda_pool_finalize_c() bind(C, name='cuda_pool_finalize')
     end subroutine

     function cuda_pool_is_initialized_c() result(flag) bind(C, name='cuda_pool_is_initialized')
       import :: c_int
       integer(c_int) :: flag
     end function
  end interface

contains

  !-----------------------------------------------------------
  ! Initialize CUDA pool using MPI node-local rank
  !-----------------------------------------------------------
  subroutine cuda_pool_init_f()
    use amr_commons, only: myid
    implicit none
    include 'mpif.h'
    integer :: node_comm, local_rank, ierr

    ! Get node-local rank for GPU assignment
    call MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, &
         myid, MPI_INFO_NULL, node_comm, ierr)
    call MPI_Comm_rank(node_comm, local_rank, ierr)
    call MPI_Comm_free(node_comm, ierr)

    call cuda_pool_init_c(int(local_rank, c_int))

    cuda_n_streams = int(cuda_get_n_streams_c())
    cuda_available = (cuda_n_streams > 0)
  end subroutine cuda_pool_init_f

  !-----------------------------------------------------------
  ! Finalize CUDA pool
  !-----------------------------------------------------------
  subroutine cuda_pool_finalize_f()
    implicit none
    call cuda_pool_finalize_c()
    cuda_available = .false.
    cuda_n_streams = 0
  end subroutine cuda_pool_finalize_f

end module cuda_commons
