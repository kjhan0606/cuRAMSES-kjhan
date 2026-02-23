!----------------------------------------------------------------------
! poisson_bench.f90 - Distributed Poisson Solver Benchmark
!
! Compares three iterative solvers for the 3D Poisson equation:
!   1) Jacobi
!   2) Red-Black Gauss-Seidel
!   3) Conjugate Gradient (CG)
!
! Usage: mpirun -np P ./poisson_bench [ngrid] [max_iter]
!        Defaults: ngrid=256, max_iter=100
!
! Slab decomposition in z. Requires N >= ncpu and N mod ncpu == 0.
!----------------------------------------------------------------------
program poisson_bench
  use mpi
  implicit none

  integer :: N, max_iter
  integer :: myid, ncpu, ierr, nloc
  integer :: left, right
  real(8) :: h, h2
  character(len=32) :: arg

  real(8), allocatable :: phi(:,:,:)       ! (N, N, 0:nloc+1)
  real(8), allocatable :: phi_new(:,:,:)   ! (N, N, nloc) - Jacobi scratch
  real(8), allocatable :: rhs(:,:,:)       ! (N, N, nloc)
  real(8), allocatable :: phi_exact(:,:,:) ! (N, N, nloc)

  ! CG work arrays
  real(8), allocatable :: r_cg(:,:,:)      ! (N, N, nloc)
  real(8), allocatable :: p_cg(:,:,:)      ! (N, N, 0:nloc+1)
  real(8), allocatable :: Ap_cg(:,:,:)     ! (N, N, nloc)

  real(8) :: t_total, t_comp, t_comm
  real(8) :: t0, t1
  real(8) :: residual, error_linf
  integer :: iters_done
  logical :: converged

  integer :: i, j, k, k_global
  real(8), parameter :: PI = 3.14159265358979323846d0
  real(8), parameter :: TWO_PI = 2.0d0 * PI
  real(8) :: x, y, z
  real(8) :: f_norm2, f_norm2_local

  call MPI_INIT(ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, ncpu, ierr)

  N = 256
  max_iter = 100

  if (command_argument_count() >= 1) then
    call get_command_argument(1, arg)
    read(arg, *) N
  end if
  if (command_argument_count() >= 2) then
    call get_command_argument(2, arg)
    read(arg, *) max_iter
  end if

  if (N < ncpu) then
    if (myid == 0) write(*,*) 'ERROR: N must be >= ncpu'
    call MPI_FINALIZE(ierr); stop
  end if
  if (mod(N, ncpu) /= 0) then
    if (myid == 0) write(*,*) 'ERROR: N must be divisible by ncpu'
    call MPI_FINALIZE(ierr); stop
  end if

  nloc = N / ncpu
  h  = 1.0d0 / dble(N)
  h2 = h * h

  left  = mod(myid - 1 + ncpu, ncpu)
  right = mod(myid + 1, ncpu)

  allocate(phi(N, N, 0:nloc+1))
  allocate(phi_new(N, N, nloc))
  allocate(rhs(N, N, nloc))
  allocate(phi_exact(N, N, nloc))
  allocate(r_cg(N, N, nloc))
  allocate(p_cg(N, N, 0:nloc+1))
  allocate(Ap_cg(N, N, nloc))

  ! Initialize RHS and exact solution
  do k = 1, nloc
    k_global = myid * nloc + k
    z = (dble(k_global) - 0.5d0) * h
    do j = 1, N
      y = (dble(j) - 0.5d0) * h
      do i = 1, N
        x = (dble(i) - 0.5d0) * h
        rhs(i, j, k) = -12.0d0 * PI**2 * sin(TWO_PI*x) * sin(TWO_PI*y) * sin(TWO_PI*z)
        phi_exact(i, j, k) = sin(TWO_PI*x) * sin(TWO_PI*y) * sin(TWO_PI*z)
      end do
    end do
  end do

  ! Compute ||f||_2 for relative residual in CG
  f_norm2_local = 0.0d0
  do k = 1, nloc
    do j = 1, N
      do i = 1, N
        f_norm2_local = f_norm2_local + rhs(i,j,k)**2
      end do
    end do
  end do
  call MPI_ALLREDUCE(f_norm2_local, f_norm2, 1, MPI_DOUBLE_PRECISION, &
                     MPI_SUM, MPI_COMM_WORLD, ierr)
  f_norm2 = sqrt(f_norm2)

  if (myid == 0) then
    write(*,'(A)')        '=== Poisson Solver Benchmark ==='
    write(*,'(A,I6,A,I6,A,I6,A,ES10.3,A,I6)') &
      'ncpu=', ncpu, ', N=', N, ', nloc=', nloc, ', h=', h, ', max_iter=', max_iter
    write(*,*)
  end if

  call MPI_BARRIER(MPI_COMM_WORLD, ierr)

  ! ===== 1) Jacobi =====
  phi = 0.0d0
  t_comp = 0.0d0; t_comm = 0.0d0

  t0 = MPI_WTIME()
  call solve_jacobi(phi, phi_new, rhs, N, nloc, h2, left, right, &
                    MPI_COMM_WORLD, max_iter, t_comp, t_comm, residual, iters_done)
  t1 = MPI_WTIME()
  t_total = t1 - t0

  call compute_error(phi, phi_exact, N, nloc, MPI_COMM_WORLD, error_linf)
  if (myid == 0) call print_result('Jacobi', iters_done, .false., residual, error_linf, &
                                    t_total, t_comp, t_comm)
  call MPI_BARRIER(MPI_COMM_WORLD, ierr)

  ! ===== 2) Red-Black GS =====
  phi = 0.0d0
  t_comp = 0.0d0; t_comm = 0.0d0

  t0 = MPI_WTIME()
  call solve_rbgs(phi, rhs, N, nloc, h2, left, right, &
                  MPI_COMM_WORLD, myid, max_iter, t_comp, t_comm, residual, iters_done)
  t1 = MPI_WTIME()
  t_total = t1 - t0

  call compute_error(phi, phi_exact, N, nloc, MPI_COMM_WORLD, error_linf)
  if (myid == 0) call print_result('Red-Black GS', iters_done, .false., residual, error_linf, &
                                    t_total, t_comp, t_comm)
  call MPI_BARRIER(MPI_COMM_WORLD, ierr)

  ! ===== 3) CG =====
  phi = 0.0d0
  t_comp = 0.0d0; t_comm = 0.0d0

  t0 = MPI_WTIME()
  call solve_cg(phi, r_cg, p_cg, Ap_cg, rhs, N, nloc, h2, left, right, &
                MPI_COMM_WORLD, max_iter, f_norm2, t_comp, t_comm, &
                residual, iters_done, converged)
  t1 = MPI_WTIME()
  t_total = t1 - t0

  call compute_error(phi, phi_exact, N, nloc, MPI_COMM_WORLD, error_linf)
  if (myid == 0) call print_result('Conjugate Gradient', iters_done, converged, residual, &
                                    error_linf, t_total, t_comp, t_comm)

  deallocate(phi, phi_new, rhs, phi_exact)
  deallocate(r_cg, p_cg, Ap_cg)

  call MPI_FINALIZE(ierr)

contains

  subroutine exchange_ghosts(phi, N, nloc, left, right, comm)
    implicit none
    integer, intent(in) :: N, nloc, left, right, comm
    real(8), intent(inout) :: phi(N, N, 0:nloc+1)
    integer :: ierr2, count
    count = N * N
    call MPI_SENDRECV(phi(1,1,1),      count, MPI_DOUBLE_PRECISION, left,  100, &
                      phi(1,1,nloc+1), count, MPI_DOUBLE_PRECISION, right, 100, &
                      comm, MPI_STATUS_IGNORE, ierr2)
    call MPI_SENDRECV(phi(1,1,nloc), count, MPI_DOUBLE_PRECISION, right, 200, &
                      phi(1,1,0),    count, MPI_DOUBLE_PRECISION, left,  200, &
                      comm, MPI_STATUS_IGNORE, ierr2)
  end subroutine exchange_ghosts

  subroutine apply_laplacian(phi, Aphi, N, nloc, h2)
    implicit none
    integer, intent(in) :: N, nloc
    real(8), intent(in) :: phi(N, N, 0:nloc+1)
    real(8), intent(out) :: Aphi(N, N, nloc)
    real(8), intent(in) :: h2
    integer :: i, j, k, ip, im, jp, jm
    do k = 1, nloc
      do j = 1, N
        jp = mod(j, N) + 1; jm = mod(j - 2 + N, N) + 1
        do i = 1, N
          ip = mod(i, N) + 1; im = mod(i - 2 + N, N) + 1
          Aphi(i,j,k) = ( phi(ip,j,k) + phi(im,j,k) &
                        + phi(i,jp,k) + phi(i,jm,k) &
                        + phi(i,j,k+1) + phi(i,j,k-1) &
                        - 6.0d0 * phi(i,j,k) ) / h2
        end do
      end do
    end do
  end subroutine apply_laplacian

  function global_dot(a, b, N, nloc, comm) result(dot)
    implicit none
    integer, intent(in) :: N, nloc, comm
    real(8), intent(in) :: a(N, N, nloc), b(N, N, nloc)
    real(8) :: dot
    real(8) :: dot_local
    integer :: i, j, k, ierr2
    dot_local = 0.0d0
    do k = 1, nloc
      do j = 1, N
        do i = 1, N
          dot_local = dot_local + a(i,j,k) * b(i,j,k)
        end do
      end do
    end do
    call MPI_ALLREDUCE(dot_local, dot, 1, MPI_DOUBLE_PRECISION, MPI_SUM, comm, ierr2)
  end function global_dot

  function compute_residual_norm(phi, rhs, N, nloc, h2, left, right, comm) result(rnorm)
    implicit none
    integer, intent(in) :: N, nloc, left, right, comm
    real(8), intent(inout) :: phi(N, N, 0:nloc+1)
    real(8), intent(in) :: rhs(N, N, nloc)
    real(8), intent(in) :: h2
    real(8) :: rnorm
    real(8) :: lap, rnorm_local
    integer :: i, j, k, ip, im, jp, jm, ierr2

    rnorm_local = 0.0d0
    do k = 1, nloc
      do j = 1, N
        jp = mod(j, N) + 1; jm = mod(j - 2 + N, N) + 1
        do i = 1, N
          ip = mod(i, N) + 1; im = mod(i - 2 + N, N) + 1
          lap = ( phi(ip,j,k) + phi(im,j,k) &
                + phi(i,jp,k) + phi(i,jm,k) &
                + phi(i,j,k+1) + phi(i,j,k-1) &
                - 6.0d0 * phi(i,j,k) ) / h2
          rnorm_local = rnorm_local + (rhs(i,j,k) - lap)**2
        end do
      end do
    end do
    call MPI_ALLREDUCE(rnorm_local, rnorm, 1, MPI_DOUBLE_PRECISION, MPI_SUM, comm, ierr2)
    rnorm = sqrt(rnorm)
  end function compute_residual_norm

  subroutine compute_error(phi, phi_exact, N, nloc, comm, error_linf)
    implicit none
    integer, intent(in) :: N, nloc, comm
    real(8), intent(in) :: phi(N, N, 0:nloc+1)
    real(8), intent(in) :: phi_exact(N, N, nloc)
    real(8), intent(out) :: error_linf
    real(8) :: err_local
    integer :: i, j, k, ierr2
    err_local = 0.0d0
    do k = 1, nloc
      do j = 1, N
        do i = 1, N
          err_local = max(err_local, abs(phi(i,j,k) - phi_exact(i,j,k)))
        end do
      end do
    end do
    call MPI_ALLREDUCE(err_local, error_linf, 1, MPI_DOUBLE_PRECISION, MPI_MAX, comm, ierr2)
  end subroutine compute_error

  subroutine print_result(name, iters, conv, residual, error, t_total, t_comp, t_comm)
    implicit none
    character(len=*), intent(in) :: name
    integer, intent(in) :: iters
    logical, intent(in) :: conv
    real(8), intent(in) :: residual, error, t_total, t_comp, t_comm
    real(8) :: pct_comp, pct_comm, per_iter
    pct_comp = 0.0d0; pct_comm = 0.0d0; per_iter = 0.0d0
    if (t_total > 0.0d0) then
      pct_comp = 100.0d0 * t_comp / t_total
      pct_comm = 100.0d0 * t_comm / t_total
    end if
    if (iters > 0) per_iter = t_total / dble(iters)

    write(*,'(A,A,A)') '--- ', trim(name), ' ---'
    if (conv) then
      write(*,'(A,I6,A,ES10.3,A,ES10.3)') &
        '  Iterations: ', iters, ' (converged),  Residual: ', residual, ',  Error: ', error
    else
      write(*,'(A,I6,A,ES10.3,A,ES10.3)') &
        '  Iterations: ', iters, ',  Residual: ', residual, ',  Error: ', error
    end if
    write(*,'(A,F8.3,A,F8.3,A,F5.1,A,F8.3,A,F5.1,A)') &
      '  Time: total=', t_total, 's  comp=', t_comp, 's(', pct_comp, &
      '%)  comm=', t_comm, 's(', pct_comm, '%)'
    write(*,'(A,F8.5,A)') '  Per iteration: ', per_iter, 's'
    write(*,*)
  end subroutine print_result

  subroutine solve_jacobi(phi, phi_new, rhs, N, nloc, h2, left, right, &
                          comm, max_iter, t_comp, t_comm, residual, iters_done)
    implicit none
    integer, intent(in) :: N, nloc, left, right, comm, max_iter
    real(8), intent(inout) :: phi(N, N, 0:nloc+1)
    real(8), intent(inout) :: phi_new(N, N, nloc)
    real(8), intent(in) :: rhs(N, N, nloc)
    real(8), intent(in) :: h2
    real(8), intent(inout) :: t_comp, t_comm
    real(8), intent(out) :: residual
    integer, intent(out) :: iters_done

    integer :: iter, i, j, k, ip, im, jp, jm
    real(8) :: tc0, tc1

    residual = 0.0d0
    do iter = 1, max_iter
      tc0 = MPI_WTIME()
      call exchange_ghosts(phi, N, nloc, left, right, comm)
      tc1 = MPI_WTIME()
      t_comm = t_comm + (tc1 - tc0)

      tc0 = MPI_WTIME()
      do k = 1, nloc
        do j = 1, N
          jp = mod(j, N) + 1; jm = mod(j - 2 + N, N) + 1
          do i = 1, N
            ip = mod(i, N) + 1; im = mod(i - 2 + N, N) + 1
            phi_new(i,j,k) = ( phi(ip,j,k) + phi(im,j,k) &
                             + phi(i,jp,k) + phi(i,jm,k) &
                             + phi(i,j,k+1) + phi(i,j,k-1) &
                             + h2 * rhs(i,j,k) ) / 6.0d0
          end do
        end do
      end do
      phi(:,:,1:nloc) = phi_new(:,:,1:nloc)
      tc1 = MPI_WTIME()
      t_comp = t_comp + (tc1 - tc0)

      if (mod(iter, 10) == 0 .or. iter == max_iter) then
        tc0 = MPI_WTIME()
        call exchange_ghosts(phi, N, nloc, left, right, comm)
        residual = compute_residual_norm(phi, rhs, N, nloc, h2, left, right, comm)
        tc1 = MPI_WTIME()
        t_comm = t_comm + (tc1 - tc0)
      end if
      iters_done = iter
    end do
  end subroutine solve_jacobi

  subroutine solve_rbgs(phi, rhs, N, nloc, h2, left, right, &
                        comm, myid, max_iter, t_comp, t_comm, residual, iters_done)
    implicit none
    integer, intent(in) :: N, nloc, left, right, comm, myid, max_iter
    real(8), intent(inout) :: phi(N, N, 0:nloc+1)
    real(8), intent(in) :: rhs(N, N, nloc)
    real(8), intent(in) :: h2
    real(8), intent(inout) :: t_comp, t_comm
    real(8), intent(out) :: residual
    integer, intent(out) :: iters_done

    integer :: iter, i, j, k, ip, im, jp, jm, kg, color
    real(8) :: tc0, tc1

    residual = 0.0d0
    do iter = 1, max_iter
      ! Red sweep
      tc0 = MPI_WTIME()
      call exchange_ghosts(phi, N, nloc, left, right, comm)
      tc1 = MPI_WTIME()
      t_comm = t_comm + (tc1 - tc0)

      tc0 = MPI_WTIME()
      do k = 1, nloc
        kg = myid * nloc + k
        do j = 1, N
          jp = mod(j, N) + 1; jm = mod(j - 2 + N, N) + 1
          do i = 1, N
            color = mod(i + j + kg, 2)
            if (color == 0) then
              ip = mod(i, N) + 1; im = mod(i - 2 + N, N) + 1
              phi(i,j,k) = ( phi(ip,j,k) + phi(im,j,k) &
                           + phi(i,jp,k) + phi(i,jm,k) &
                           + phi(i,j,k+1) + phi(i,j,k-1) &
                           + h2 * rhs(i,j,k) ) / 6.0d0
            end if
          end do
        end do
      end do
      tc1 = MPI_WTIME()
      t_comp = t_comp + (tc1 - tc0)

      ! Black sweep
      tc0 = MPI_WTIME()
      call exchange_ghosts(phi, N, nloc, left, right, comm)
      tc1 = MPI_WTIME()
      t_comm = t_comm + (tc1 - tc0)

      tc0 = MPI_WTIME()
      do k = 1, nloc
        kg = myid * nloc + k
        do j = 1, N
          jp = mod(j, N) + 1; jm = mod(j - 2 + N, N) + 1
          do i = 1, N
            color = mod(i + j + kg, 2)
            if (color == 1) then
              ip = mod(i, N) + 1; im = mod(i - 2 + N, N) + 1
              phi(i,j,k) = ( phi(ip,j,k) + phi(im,j,k) &
                           + phi(i,jp,k) + phi(i,jm,k) &
                           + phi(i,j,k+1) + phi(i,j,k-1) &
                           + h2 * rhs(i,j,k) ) / 6.0d0
            end if
          end do
        end do
      end do
      tc1 = MPI_WTIME()
      t_comp = t_comp + (tc1 - tc0)

      if (mod(iter, 10) == 0 .or. iter == max_iter) then
        tc0 = MPI_WTIME()
        call exchange_ghosts(phi, N, nloc, left, right, comm)
        residual = compute_residual_norm(phi, rhs, N, nloc, h2, left, right, comm)
        tc1 = MPI_WTIME()
        t_comm = t_comm + (tc1 - tc0)
      end if
      iters_done = iter
    end do
  end subroutine solve_rbgs

  subroutine solve_cg(phi, r, p, Ap, rhs, N, nloc, h2, left, right, &
                      comm, max_iter, f_norm2, t_comp, t_comm, &
                      residual, iters_done, converged)
    implicit none
    integer, intent(in) :: N, nloc, left, right, comm, max_iter
    real(8), intent(inout) :: phi(N, N, 0:nloc+1)
    real(8), intent(inout) :: r(N, N, nloc)
    real(8), intent(inout) :: p(N, N, 0:nloc+1)
    real(8), intent(inout) :: Ap(N, N, nloc)
    real(8), intent(in) :: rhs(N, N, nloc)
    real(8), intent(in) :: h2, f_norm2
    real(8), intent(inout) :: t_comp, t_comm
    real(8), intent(out) :: residual
    integer, intent(out) :: iters_done
    logical, intent(out) :: converged

    integer :: iter, i, j, k
    real(8) :: rr, rr_new, pAp, alpha, beta
    real(8) :: tc0, tc1
    real(8), parameter :: TOL = 1.0d-10

    converged = .false.

    tc0 = MPI_WTIME()
    call exchange_ghosts(phi, N, nloc, left, right, comm)
    tc1 = MPI_WTIME()
    t_comm = t_comm + (tc1 - tc0)

    tc0 = MPI_WTIME()
    call apply_laplacian(phi, Ap, N, nloc, h2)
    do k = 1, nloc
      do j = 1, N
        do i = 1, N
          r(i,j,k) = rhs(i,j,k) - Ap(i,j,k)
        end do
      end do
    end do
    p(:,:,1:nloc) = r(:,:,1:nloc)
    tc1 = MPI_WTIME()
    t_comp = t_comp + (tc1 - tc0)

    tc0 = MPI_WTIME()
    rr = global_dot(r, r, N, nloc, comm)
    tc1 = MPI_WTIME()
    t_comm = t_comm + (tc1 - tc0)

    do iter = 1, max_iter
      tc0 = MPI_WTIME()
      call exchange_ghosts(p, N, nloc, left, right, comm)
      tc1 = MPI_WTIME()
      t_comm = t_comm + (tc1 - tc0)

      tc0 = MPI_WTIME()
      call apply_laplacian(p, Ap, N, nloc, h2)
      tc1 = MPI_WTIME()
      t_comp = t_comp + (tc1 - tc0)

      tc0 = MPI_WTIME()
      pAp = global_dot(p(:,:,1:nloc), Ap, N, nloc, comm)
      tc1 = MPI_WTIME()
      t_comm = t_comm + (tc1 - tc0)

      alpha = rr / pAp

      tc0 = MPI_WTIME()
      do k = 1, nloc
        do j = 1, N
          do i = 1, N
            phi(i,j,k) = phi(i,j,k) + alpha * p(i,j,k)
            r(i,j,k)   = r(i,j,k)   - alpha * Ap(i,j,k)
          end do
        end do
      end do
      tc1 = MPI_WTIME()
      t_comp = t_comp + (tc1 - tc0)

      tc0 = MPI_WTIME()
      rr_new = global_dot(r, r, N, nloc, comm)
      tc1 = MPI_WTIME()
      t_comm = t_comm + (tc1 - tc0)

      residual = sqrt(rr_new)
      iters_done = iter

      if (f_norm2 > 0.0d0) then
        if (sqrt(rr_new) / f_norm2 < TOL) then
          converged = .true.; exit
        end if
      else
        if (sqrt(rr_new) < TOL) then
          converged = .true.; exit
        end if
      end if

      beta = rr_new / rr

      tc0 = MPI_WTIME()
      do k = 1, nloc
        do j = 1, N
          do i = 1, N
            p(i,j,k) = r(i,j,k) + beta * p(i,j,k)
          end do
        end do
      end do
      tc1 = MPI_WTIME()
      t_comp = t_comp + (tc1 - tc0)

      rr = rr_new
    end do
  end subroutine solve_cg

end program poisson_bench
