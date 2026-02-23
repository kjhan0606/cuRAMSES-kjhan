!=====================================================================
! ghost_bench.f90 - Ghost Zone Exchange Benchmark (v2)
!
! Compares ghost zone communication patterns:
!   Hilbert vs K-Section domain decomposition
!   P2P ISEND/IRECV vs MPI_ALLTOALLV vs K-Section hierarchical exchange
!
! 6 combinations: 2 DD x 3 Exchange methods
!
! v2: Plane-by-plane ownership computation for large grid support.
!     Memory: O(N^2) per rank instead of O(N^3).
!
! Usage: mpirun -np P ./ghost_bench [ngrid] [nvar] [nrepeat]
!=====================================================================

!---------------------------------------------------------------------
! Hilbert curve module (from dd_bench.f90)
!---------------------------------------------------------------------
module hilbert_mod
  implicit none
contains

  subroutine hilbert3d(x, y, z, order, bit_length, npoint)
    implicit none
    integer, intent(in) :: bit_length, npoint
    integer, intent(in) :: x(npoint), y(npoint), z(npoint)
    real(8), intent(out) :: order(npoint)
    logical :: ib(0:3*bit_length-1)
    logical :: xb(0:bit_length-1), yb(0:bit_length-1), zb(0:bit_length-1)
    integer :: sd(0:7,0:1,0:11)
    integer :: i, ip, cs, ns, b0, b1, b2, si, hi

    sd = RESHAPE((/ &
      1,2,3,2,4,5,3,5, 0,1,3,2,7,6,4,5, &
      2,6,0,7,8,8,0,7, 0,7,1,6,3,4,2,5, &
      0,9,10,9,1,1,11,11, 0,3,7,4,1,2,6,5, &
      6,0,6,11,9,0,9,8, 2,3,1,0,5,4,6,7, &
      11,11,0,7,5,9,0,7, 4,3,5,2,7,0,6,1, &
      4,4,8,8,0,6,10,6, 6,5,1,2,7,4,0,3, &
      5,7,5,3,1,1,11,11, 4,7,3,0,5,6,2,1, &
      6,1,6,10,9,4,9,10, 6,7,5,4,1,0,2,3, &
      10,3,1,1,10,3,5,9, 2,5,3,4,1,6,0,7, &
      4,4,8,8,2,7,2,3, 2,1,5,6,3,0,4,7, &
      7,2,11,2,7,5,8,5, 4,5,7,6,3,2,0,1, &
      10,3,2,6,10,3,4,4, 6,1,7,0,5,2,4,3/), (/8,2,12/))

    do ip = 1, npoint
      do i = 0, bit_length-1
        xb(i) = btest(x(ip),i); yb(i) = btest(y(ip),i); zb(i) = btest(z(ip),i)
      end do
      do i = 0, bit_length-1
        ib(3*i+2) = xb(i); ib(3*i+1) = yb(i); ib(3*i) = zb(i)
      end do
      cs = 0
      do i = bit_length-1, 0, -1
        b2=0; if(ib(3*i+2)) b2=1
        b1=0; if(ib(3*i+1)) b1=1
        b0=0; if(ib(3*i  )) b0=1
        si = b2*4+b1*2+b0
        ns = sd(si,0,cs); hi = sd(si,1,cs)
        ib(3*i+2) = btest(hi,2); ib(3*i+1) = btest(hi,1); ib(3*i) = btest(hi,0)
        cs = ns
      end do
      order(ip) = 0.0d0
      do i = 0, 3*bit_length-1
        b0=0; if(ib(i)) b0=1
        order(ip) = order(ip) + real(b0,kind=8)*real(2,kind=8)**i
      end do
    end do
  end subroutine hilbert3d

end module hilbert_mod

!---------------------------------------------------------------------
! K-Section module: tree DD + generalized hierarchical exchange
!---------------------------------------------------------------------
module ksection_mod
  use mpi
  implicit none

  integer, parameter :: KS_MAXLEV = 32
  integer :: nklev, kfac(KS_MAXLEV), kkmax

  integer :: knn
  integer, allocatable :: kcmin(:), kcmax(:)
  real(8), allocatable :: kbmin(:,:), kbmax(:,:)
  integer, allocatable :: knc(:), kco(:), kdir(:)

  integer :: kpnode(KS_MAXLEV+1), kpchild(KS_MAXLEV)

contains

  subroutine ks_factorize(n)
    implicit none
    integer, intent(in) :: n
    integer :: m, d, nf, i, j, tmp

    kfac = 0; nf = 0; m = n; d = 2
    do while(m > 1)
      do while(mod(m,d) == 0)
        nf = nf + 1; kfac(nf) = d; m = m / d
      end do
      d = d + 1
      if(d*d > m .and. m > 1) then
        nf = nf + 1; kfac(nf) = m; m = 1
      end if
    end do
    do i = 1, nf-1
      do j = i+1, nf
        if(kfac(j) > kfac(i)) then
          tmp = kfac(i); kfac(i) = kfac(j); kfac(j) = tmp
        end if
      end do
    end do
    nklev = nf
    kkmax = 2
    if(nf > 0) kkmax = kfac(1)
  end subroutine

  subroutine ks_build(ncpu)
    implicit none
    integer, intent(in) :: ncpu
    integer :: mxn, lvl, k, nd, cn, j, dir
    integer :: ls, le, ns, cppc, cs, ce

    call ks_factorize(ncpu)

    mxn = 2*ncpu + 16
    if(allocated(kcmin)) deallocate(kcmin,kcmax,kbmin,kbmax,knc,kco,kdir)
    allocate(kcmin(mxn), kcmax(mxn), kbmin(3,mxn), kbmax(3,mxn))
    allocate(knc(mxn), kco(mxn), kdir(mxn))

    knn = 1
    kcmin(1) = 0; kcmax(1) = ncpu - 1
    kbmin(:,1) = 0.0d0; kbmax(:,1) = 1.0d0
    knc(1) = 0; kco(1) = 0; kdir(1) = 0

    ls = 1; le = 1
    do lvl = 1, nklev
      k = kfac(lvl); ns = knn + 1
      do nd = ls, le
        dir = 1
        if(kbmax(2,nd)-kbmin(2,nd) > kbmax(dir,nd)-kbmin(dir,nd)) dir = 2
        if(kbmax(3,nd)-kbmin(3,nd) > kbmax(dir,nd)-kbmin(dir,nd)) dir = 3
        kdir(nd) = dir; kco(nd) = knn; knc(nd) = k
        cppc = (kcmax(nd) - kcmin(nd) + 1) / k; cs = kcmin(nd)
        do j = 1, k
          knn = knn + 1; cn = knn
          kbmin(:,cn) = kbmin(:,nd); kbmax(:,cn) = kbmax(:,nd)
          kbmin(dir,cn) = kbmin(dir,nd) + dble(j-1)/dble(k)*(kbmax(dir,nd)-kbmin(dir,nd))
          kbmax(dir,cn) = kbmin(dir,nd) + dble(j)  /dble(k)*(kbmax(dir,nd)-kbmin(dir,nd))
          ce = cs + cppc - 1; kcmin(cn) = cs; kcmax(cn) = ce; cs = ce + 1
          knc(cn) = 0; kco(cn) = 0; kdir(cn) = 0
        end do
      end do
      ls = ns; le = knn
    end do
  end subroutine

  subroutine ks_path(myid)
    implicit none
    integer, intent(in) :: myid
    integer :: lvl, nd, j, cn

    nd = 1; kpnode(1) = 1
    do lvl = 1, nklev
      do j = 1, knc(nd)
        cn = kco(nd) + j
        if(myid >= kcmin(cn) .and. myid <= kcmax(cn)) then
          kpchild(lvl) = j; kpnode(lvl+1) = cn; nd = cn; exit
        end if
      end do
    end do
  end subroutine

  subroutine ks_dd(pts, np, dest)
    implicit none
    integer, intent(in) :: np
    real(8), intent(in) :: pts(3, np)
    integer, intent(out) :: dest(np)
    integer :: i, nd, lvl, j, k, dir
    real(8) :: c, lo, hi

    do i = 1, np
      nd = 1
      do lvl = 1, nklev
        k = knc(nd); if(k == 0) exit
        dir = kdir(nd)
        c = pts(dir, i); lo = kbmin(dir, nd); hi = kbmax(dir, nd)
        j = int((c - lo) / (hi - lo) * dble(k)) + 1
        j = max(1, min(k, j))
        nd = kco(nd) + j
      end do
      dest(i) = kcmin(nd)
    end do
  end subroutine

  subroutine ks_exchange_gen(sbuf, sdest, ns, rbuf, nr, ndata, myid, ncpu, comm)
    implicit none
    integer, intent(in) :: ns, ndata, myid, ncpu, comm
    real(8), intent(in) :: sbuf(ndata, ns)
    integer, intent(in) :: sdest(ns)
    real(8), allocatable, intent(out) :: rbuf(:,:)
    integer, intent(out) :: nr

    real(8), allocatable :: wb(:,:), w2(:,:), tb(:,:)
    integer, allocatable :: wd(:), d2(:), td(:)
    integer :: nw, nw2
    integer :: lvl, k, mc, nd, cn, mp, cppc
    integer :: j, p, np2, roff, nrq, ierr, myc, trv
    integer, allocatable :: cc(:), of(:)
    integer, allocatable :: pl(:), psc(:), prc(:)
    integer, allocatable :: rq(:), st(:,:)

    nw = ns
    allocate(wb(ndata, max(1,nw)), wd(max(1,nw)))
    if(nw > 0) then
      wb(:,1:nw) = sbuf(:,1:nw)
      wd(1:nw) = sdest(1:nw)
    end if

    allocate(cc(kkmax), of(kkmax+1))
    allocate(pl(kkmax), psc(kkmax), prc(kkmax))
    allocate(rq(6*kkmax), st(MPI_STATUS_SIZE, 6*kkmax))

    nd = 1
    do lvl = 1, nklev
      k = knc(nd)
      mc = kpchild(lvl)
      cn = kco(nd) + mc
      mp = myid - kcmin(cn)
      cppc = (kcmax(nd) - kcmin(nd) + 1) / k

      cc(1:k) = 0
      do j = 1, nw
        p = (wd(j) - kcmin(nd)) / cppc + 1
        p = max(1, min(k, p))
        cc(p) = cc(p) + 1
      end do

      of(1) = 0
      do j = 1, k
        of(j+1) = of(j) + cc(j)
      end do

      allocate(w2(ndata, max(1,nw)), d2(max(1,nw)))
      cc(1:k) = 0
      do j = 1, nw
        p = (wd(j) - kcmin(nd)) / cppc + 1
        p = max(1, min(k, p))
        cc(p) = cc(p) + 1
        w2(:, of(p)+cc(p)) = wb(:, j)
        d2(of(p)+cc(p)) = wd(j)
      end do

      np2 = 0
      do j = 1, k
        if(j == mc) cycle
        np2 = np2 + 1
        pl(np2) = kcmin(kco(nd)+j) + min(mp, kcmax(kco(nd)+j)-kcmin(kco(nd)+j))
        psc(np2) = cc(j)
      end do

      nrq = 0
      do p = 1, np2
        nrq = nrq + 1
        call MPI_ISEND(psc(p), 1, MPI_INTEGER, pl(p), 100+lvl, comm, rq(nrq), ierr)
        nrq = nrq + 1
        call MPI_IRECV(prc(p), 1, MPI_INTEGER, pl(p), 100+lvl, comm, rq(nrq), ierr)
      end do
      if(nrq > 0) call MPI_WAITALL(nrq, rq, st, ierr)

      trv = 0
      do p = 1, np2
        trv = trv + prc(p)
      end do
      myc = of(mc+1) - of(mc)

      nw2 = myc + trv
      allocate(tb(ndata, max(1,nw2)), td(max(1,nw2)))
      if(myc > 0) then
        tb(:, 1:myc) = w2(:, of(mc)+1:of(mc)+myc)
        td(1:myc) = d2(of(mc)+1:of(mc)+myc)
      end if

      nrq = 0
      p = 0
      do j = 1, k
        if(j == mc) cycle
        p = p + 1
        if(psc(p) > 0) then
          nrq = nrq + 1
          call MPI_ISEND(w2(1,of(j)+1), ndata*psc(p), MPI_DOUBLE_PRECISION, &
               pl(p), 200+lvl, comm, rq(nrq), ierr)
          nrq = nrq + 1
          call MPI_ISEND(d2(of(j)+1), psc(p), MPI_INTEGER, &
               pl(p), 300+lvl, comm, rq(nrq), ierr)
        end if
      end do

      roff = myc
      do p = 1, np2
        if(prc(p) > 0) then
          nrq = nrq + 1
          call MPI_IRECV(tb(1,roff+1), ndata*prc(p), MPI_DOUBLE_PRECISION, &
               pl(p), 200+lvl, comm, rq(nrq), ierr)
          nrq = nrq + 1
          call MPI_IRECV(td(roff+1), prc(p), MPI_INTEGER, &
               pl(p), 300+lvl, comm, rq(nrq), ierr)
        end if
        roff = roff + prc(p)
      end do
      if(nrq > 0) call MPI_WAITALL(nrq, rq, st, ierr)

      deallocate(wb, wd, w2, d2)
      nw = nw2
      allocate(wb(ndata, max(1,nw)), wd(max(1,nw)))
      if(nw > 0) then
        wb(:, 1:nw) = tb(:, 1:nw)
        wd(1:nw) = td(1:nw)
      end if
      deallocate(tb, td)

      nd = cn
    end do

    nr = nw
    allocate(rbuf(ndata, max(1,nr)))
    if(nr > 0) rbuf(:, 1:nr) = wb(:, 1:nr)
    deallocate(wb, wd, cc, of, pl, psc, prc, rq, st)
  end subroutine

end module ksection_mod

!---------------------------------------------------------------------
! Ghost zone communication module
!---------------------------------------------------------------------
module ghost_mod
  use mpi
  use ksection_mod, only: ks_exchange_gen
  implicit none

  type ghost_comm
    integer :: ncpu, myid, ngrid
    integer, allocatable :: send_cnt(:), recv_cnt(:)   ! (0:ncpu-1)
    integer, allocatable :: send_off(:), recv_off(:)   ! (0:ncpu-1)
    integer, allocatable :: send_cells(:,:)             ! (3, nsend_tot)
    integer, allocatable :: recv_cells(:,:)             ! (3, nrecv_tot)
    integer :: nsend_tot, nrecv_tot, n_partner
  end type

contains

  pure real(8) function cell_value(i, j, k, ivar, N)
    integer, intent(in) :: i, j, k, ivar, N
    cell_value = dble(i) + dble(j)*dble(N) + dble(k)*dble(N)**2 + dble(ivar)*dble(N)**3
  end function

  subroutine ghost_p2p_exchange(gc, nvar, comm, nerr)
    type(ghost_comm), intent(in) :: gc
    integer, intent(in) :: nvar, comm
    integer, intent(out) :: nerr

    real(8), allocatable :: sendbuf(:,:), recvbuf(:,:)
    integer, allocatable :: rq(:), st(:,:)
    integer :: nrq, p, ii, iv, ierr, N

    N = gc%ngrid; nerr = 0
    allocate(sendbuf(nvar, max(1,gc%nsend_tot)))
    allocate(recvbuf(nvar, max(1,gc%nrecv_tot)))
    allocate(rq(2*gc%ncpu), st(MPI_STATUS_SIZE, 2*gc%ncpu))

    do ii = 1, gc%nsend_tot
      do iv = 1, nvar
        sendbuf(iv,ii) = cell_value(gc%send_cells(1,ii), &
             gc%send_cells(2,ii), gc%send_cells(3,ii), iv, N)
      end do
    end do

    nrq = 0
    do p = 0, gc%ncpu-1
      if(gc%recv_cnt(p) == 0) cycle
      nrq = nrq + 1
      call MPI_IRECV(recvbuf(1,gc%recv_off(p)+1), nvar*gc%recv_cnt(p), &
           MPI_DOUBLE_PRECISION, p, 500, comm, rq(nrq), ierr)
    end do

    do p = 0, gc%ncpu-1
      if(gc%send_cnt(p) == 0) cycle
      nrq = nrq + 1
      call MPI_ISEND(sendbuf(1,gc%send_off(p)+1), nvar*gc%send_cnt(p), &
           MPI_DOUBLE_PRECISION, p, 500, comm, rq(nrq), ierr)
    end do

    if(nrq > 0) call MPI_WAITALL(nrq, rq, st, ierr)

    do ii = 1, gc%nrecv_tot
      do iv = 1, nvar
        if(recvbuf(iv,ii) /= cell_value(gc%recv_cells(1,ii), &
             gc%recv_cells(2,ii), gc%recv_cells(3,ii), iv, N)) nerr = nerr + 1
      end do
    end do

    deallocate(sendbuf, recvbuf, rq, st)
  end subroutine

  subroutine ghost_a2av_exchange(gc, nvar, comm, nerr)
    type(ghost_comm), intent(in) :: gc
    integer, intent(in) :: nvar, comm
    integer, intent(out) :: nerr

    real(8), allocatable :: sendbuf(:,:), recvbuf(:,:)
    integer, allocatable :: scnt(:), rcnt(:), sdsp(:), rdsp(:)
    integer :: p, ii, iv, ierr, N

    N = gc%ngrid; nerr = 0
    allocate(sendbuf(nvar, max(1,gc%nsend_tot)))
    allocate(recvbuf(nvar, max(1,gc%nrecv_tot)))

    do ii = 1, gc%nsend_tot
      do iv = 1, nvar
        sendbuf(iv,ii) = cell_value(gc%send_cells(1,ii), &
             gc%send_cells(2,ii), gc%send_cells(3,ii), iv, N)
      end do
    end do

    allocate(scnt(0:gc%ncpu-1), rcnt(0:gc%ncpu-1))
    allocate(sdsp(0:gc%ncpu-1), rdsp(0:gc%ncpu-1))
    do p = 0, gc%ncpu-1
      scnt(p) = gc%send_cnt(p) * nvar
      rcnt(p) = gc%recv_cnt(p) * nvar
      sdsp(p) = gc%send_off(p) * nvar
      rdsp(p) = gc%recv_off(p) * nvar
    end do

    call MPI_ALLTOALLV(sendbuf, scnt, sdsp, MPI_DOUBLE_PRECISION, &
                       recvbuf, rcnt, rdsp, MPI_DOUBLE_PRECISION, comm, ierr)

    do ii = 1, gc%nrecv_tot
      do iv = 1, nvar
        if(recvbuf(iv,ii) /= cell_value(gc%recv_cells(1,ii), &
             gc%recv_cells(2,ii), gc%recv_cells(3,ii), iv, N)) nerr = nerr + 1
      end do
    end do

    deallocate(sendbuf, recvbuf, scnt, rcnt, sdsp, rdsp)
  end subroutine

  ! K-Section hierarchical ghost exchange (integer(8) safe for large N)
  subroutine ghost_ks_exchange(gc, nvar, myid, ncpu, comm, nerr)
    type(ghost_comm), intent(in) :: gc
    integer, intent(in) :: nvar, myid, ncpu, comm
    integer, intent(out) :: nerr

    real(8), allocatable :: sendbuf(:,:), recvbuf(:,:)
    integer, allocatable :: sdest(:)
    integer :: N, ii, iv, p, q, nr
    integer(8) :: lid8, N8

    N = gc%ngrid; nerr = 0; N8 = int(N, 8)

    allocate(sendbuf(nvar+1, max(1, gc%nsend_tot)))
    allocate(sdest(max(1, gc%nsend_tot)))

    do p = 0, gc%ncpu-1
      do q = 1, gc%send_cnt(p)
        ii = gc%send_off(p) + q
        do iv = 1, nvar
          sendbuf(iv, ii) = cell_value(gc%send_cells(1,ii), &
               gc%send_cells(2,ii), gc%send_cells(3,ii), iv, N)
        end do
        ! Encode linear cell index as real(8) — avoids integer overflow for large N
        sendbuf(nvar+1, ii) = dble(gc%send_cells(1,ii)-1) + &
             dble(gc%send_cells(2,ii)-1)*dble(N) + &
             dble(gc%send_cells(3,ii)-1)*dble(N)*dble(N)
        sdest(ii) = p
      end do
    end do

    call ks_exchange_gen(sendbuf, sdest, gc%nsend_tot, recvbuf, nr, &
         nvar+1, myid, ncpu, comm)

    do ii = 1, nr
      lid8 = nint(recvbuf(nvar+1, ii), kind=8)
      do iv = 1, nvar
        if(recvbuf(iv, ii) /= cell_value( &
             int(mod(lid8, N8)) + 1, &
             int(mod(lid8/N8, N8)) + 1, &
             int(lid8/(N8*N8)) + 1, iv, N)) nerr = nerr + 1
      end do
    end do

    deallocate(sendbuf, sdest)
    if(allocated(recvbuf)) deallocate(recvbuf)
  end subroutine

  subroutine ghost_free(gc)
    type(ghost_comm), intent(inout) :: gc
    if(allocated(gc%send_cnt))   deallocate(gc%send_cnt)
    if(allocated(gc%recv_cnt))   deallocate(gc%recv_cnt)
    if(allocated(gc%send_off))   deallocate(gc%send_off)
    if(allocated(gc%recv_off))   deallocate(gc%recv_off)
    if(allocated(gc%send_cells)) deallocate(gc%send_cells)
    if(allocated(gc%recv_cells)) deallocate(gc%recv_cells)
  end subroutine

end module ghost_mod

!=====================================================================
! Main benchmark program
!=====================================================================
program ghost_bench
  use mpi
  use hilbert_mod
  use ksection_mod
  use ghost_mod
  implicit none

  integer :: myid, ncpu, ierr, comm
  integer :: ngrid, nvar, nrepeat
  character(len=32) :: arg
  integer(8) :: ncells8

  type(ghost_comm) :: gc(2)
  integer :: i, idd, iex, ir, nerr, nerr_total

  real(8) :: t0, t1, t_build(2)
  real(8), allocatable :: dt_exch(:)
  real(8) :: m_exch, s_exch

  ! DD statistics
  integer :: sum_send, sum_recv, sum_partner
  real(8) :: avg_owned, avg_send, avg_recv, avg_partner, sv_ratio

  ! Hilbert parameters
  integer :: bit_length, hgrid, nplane
  real(8) :: max_key

  ! Work arrays for plane computation (host-associated by contained subs)
  integer, allocatable :: work_ix(:), work_iy(:), work_iz(:), work_dd(:)
  real(8), allocatable :: work_hk(:), work_pts(:,:)

  character(len=10) :: dd_name(2), ex_name(3)

  call MPI_INIT(ierr)
  comm = MPI_COMM_WORLD
  call MPI_COMM_RANK(comm, myid, ierr)
  call MPI_COMM_SIZE(comm, ncpu, ierr)

  ! Parse command line arguments
  ngrid = 128; nvar = 5; nrepeat = 5
  if(command_argument_count() >= 1) then
    call get_command_argument(1, arg); read(arg,*) ngrid
  end if
  if(command_argument_count() >= 2) then
    call get_command_argument(2, arg); read(arg,*) nvar
  end if
  if(command_argument_count() >= 3) then
    call get_command_argument(3, arg); read(arg,*) nrepeat
  end if

  ncells8 = int(ngrid, 8)**3
  nplane = ngrid * ngrid

  if(myid == 0) then
    write(*,'(A)') '=== Ghost Zone Exchange Benchmark ==='
    write(*,'(A,I0,A,I0,A,I0,A,I0,A,I0)') &
         'ncpu=', ncpu, ', N=', ngrid, ' (N^3=', ncells8, &
         '), nvar=', nvar, ', nrepeat=', nrepeat
  end if

  ! Pre-compute Hilbert parameters
  bit_length = 1
  do while(2**bit_length < ngrid)
    bit_length = bit_length + 1
  end do
  hgrid = 2**bit_length
  max_key = 2.0d0**(3*bit_length)

  ! Build K-Section tree + path
  call ks_build(ncpu)
  call ks_path(myid)

  if(myid == 0) then
    write(*,'(A)', advance='no') 'K-Section factors: '
    do i = 1, nklev
      if(i > 1) write(*,'(A)', advance='no') ' x '
      write(*,'(I0)', advance='no') kfac(i)
    end do
    write(*,'(A,I0,A)') '  (', nklev, ' levels)'
    write(*,*)
  end if

  ! Allocate work arrays for plane computation
  allocate(work_ix(nplane), work_iy(nplane), work_iz(nplane))
  allocate(work_hk(nplane), work_pts(3, nplane), work_dd(nplane))
  allocate(dt_exch(nrepeat))

  dd_name(1) = 'Hilbert  '; dd_name(2) = 'K-Section'
  ex_name(1) = 'P2P      '; ex_name(2) = 'ALLTOALLV'; ex_name(3) = 'K-Section'

  ! ===== Build ghost_comm, time it, print DD stats =====
  do idd = 1, 2
    if(myid == 0) then
      write(*,'(A,A,A)', advance='no') 'Building ', trim(dd_name(idd)), ' ghost lists...'
      call flush(6)
    end if

    call MPI_BARRIER(comm, ierr); t0 = MPI_WTIME()
    call build_ghost_comm(idd, ngrid, myid, ncpu, gc(idd))
    call MPI_BARRIER(comm, ierr); t1 = MPI_WTIME()
    t_build(idd) = t1 - t0

    if(myid == 0) write(*,'(A,F10.2,A)') ' done (', t_build(idd), 's)'

    sum_send = gc(idd)%nsend_tot; sum_recv = gc(idd)%nrecv_tot
    sum_partner = gc(idd)%n_partner
    call MPI_ALLREDUCE(MPI_IN_PLACE, sum_send,    1, MPI_INTEGER, MPI_SUM, comm, ierr)
    call MPI_ALLREDUCE(MPI_IN_PLACE, sum_recv,    1, MPI_INTEGER, MPI_SUM, comm, ierr)
    call MPI_ALLREDUCE(MPI_IN_PLACE, sum_partner, 1, MPI_INTEGER, MPI_SUM, comm, ierr)

    if(myid == 0) then
      avg_owned   = dble(ncells8) / dble(ncpu)
      avg_send    = dble(sum_send) / dble(ncpu)
      avg_recv    = dble(sum_recv) / dble(ncpu)
      avg_partner = dble(sum_partner) / dble(ncpu)
      sv_ratio    = avg_send / avg_owned

      write(*,'(A,A,A)') '--- ', trim(dd_name(idd)), ' DD ---'
      write(*,'(A,F14.1,A,F14.1,A,F14.1)') &
           'Owned cells: ', avg_owned, ' (avg), Ghost send: ', avg_send, &
           ', Ghost recv: ', avg_recv
      write(*,'(A,F5.1,A,F7.4)') &
           'Neighbor CPUs: ', avg_partner, ' (avg), Surface/Volume: ', sv_ratio
      write(*,*)
    end if
  end do

  ! ===== Table header =====
  if(myid == 0) then
    write(*,'(A21, A14, A18, A20)') &
         'Combination', 'Build(s)', 'Exch(s)', 'Total(s)'
    write(*,'(A)') repeat('-', 73)
  end if

  ! ===== Benchmark 6 combinations =====
  do idd = 1, 2
    do iex = 1, 3
      ! Warmup + verify
      call run_exchange(gc(idd), iex, nerr)
      nerr_total = nerr
      call MPI_ALLREDUCE(MPI_IN_PLACE, nerr_total, 1, MPI_INTEGER, MPI_SUM, comm, ierr)
      if(nerr_total > 0 .and. myid == 0) then
        write(*,'(A,A,A,A,A,I0)') '  *** VERIFY FAIL: ', trim(dd_name(idd)), &
             ' + ', trim(ex_name(iex)), ', errors=', nerr_total
      end if

      ! Timed repeats
      do ir = 1, nrepeat
        call MPI_BARRIER(comm, ierr); t0 = MPI_WTIME()
        call run_exchange(gc(idd), iex, nerr)
        call MPI_BARRIER(comm, ierr); t1 = MPI_WTIME()
        dt_exch(ir) = t1 - t0
      end do
      call compute_stats(dt_exch, nrepeat, m_exch, s_exch, comm)

      if(myid == 0) then
        write(*,'(A9,A3,A9, F14.4, F10.4,A2,F6.4, F12.4,A2,F6.4)') &
             dd_name(idd), ' + ', ex_name(iex), &
             t_build(idd), &
             m_exch, '+-', s_exch, &
             t_build(idd)+m_exch, '+-', s_exch
      end if
    end do
  end do

  if(myid == 0) write(*,'(A)') repeat('-', 73)

  ! Cleanup
  call ghost_free(gc(1)); call ghost_free(gc(2))
  deallocate(work_ix, work_iy, work_iz, work_hk, work_pts, work_dd)
  deallocate(dt_exch)
  call MPI_FINALIZE(ierr)

contains

  !-------------------------------------------------------------------
  ! Compute ownership of all N^2 cells in z-plane kplane (1..N)
  ! Uses host-associated: bit_length, hgrid, max_key, ncpu, work_*
  !-------------------------------------------------------------------
  subroutine compute_plane_owner(idd_type, kplane, N, owners)
    integer, intent(in) :: idd_type, kplane, N
    integer, intent(out) :: owners(N, N)
    integer :: ii, jj, pp, np, iz_val
    real(8) :: zcoord

    np = N * N

    if(idd_type == 1) then
      ! Hilbert DD
      iz_val = min(int((dble(kplane)-0.5d0)/dble(N)*dble(hgrid)), hgrid-1)
      zcoord = (dble(kplane)-0.5d0)/dble(N)
      pp = 0
      do jj = 1, N; do ii = 1, N
        pp = pp + 1
        work_ix(pp) = min(int((dble(ii)-0.5d0)/dble(N)*dble(hgrid)), hgrid-1)
        work_iy(pp) = min(int((dble(jj)-0.5d0)/dble(N)*dble(hgrid)), hgrid-1)
        work_iz(pp) = iz_val
      end do; end do
      call hilbert3d(work_ix, work_iy, work_iz, work_hk, bit_length, np)
      pp = 0
      do jj = 1, N; do ii = 1, N
        pp = pp + 1
        owners(ii, jj) = min(ncpu-1, int(work_hk(pp)/max_key*dble(ncpu)))
      end do; end do
    else
      ! K-Section DD
      zcoord = (dble(kplane)-0.5d0)/dble(N)
      pp = 0
      do jj = 1, N; do ii = 1, N
        pp = pp + 1
        work_pts(1, pp) = (dble(ii)-0.5d0)/dble(N)
        work_pts(2, pp) = (dble(jj)-0.5d0)/dble(N)
        work_pts(3, pp) = zcoord
      end do; end do
      call ks_dd(work_pts, np, work_dd)
      pp = 0
      do jj = 1, N; do ii = 1, N
        pp = pp + 1
        owners(ii, jj) = work_dd(pp)
      end do; end do
    end if
  end subroutine

  !-------------------------------------------------------------------
  ! Build ghost zone communication lists using 3-plane sliding window
  ! Two passes: count (pass 1) then fill (pass 2)
  ! Memory: O(N^2) per rank instead of O(N^3)
  !-------------------------------------------------------------------
  subroutine build_ghost_comm(idd_type, N, myid_in, ncpu_in, gc_out)
    integer, intent(in) :: idd_type, N, myid_in, ncpu_in
    type(ghost_comm), intent(out) :: gc_out

    integer, allocatable :: planes(:,:,:)
    integer :: idx_prev, idx_curr, idx_next, itmp
    integer :: i, j, k, kn, d, p, q, m, n_nbr, ni, nj
    integer :: nbr_cpu(6)
    logical :: is_new
    integer, allocatable :: stmp(:), rtmp(:)
    integer :: ipass

    allocate(planes(N, N, 0:2))

    gc_out%ncpu = ncpu_in; gc_out%myid = myid_in; gc_out%ngrid = N
    allocate(gc_out%send_cnt(0:ncpu_in-1), gc_out%recv_cnt(0:ncpu_in-1))
    allocate(gc_out%send_off(0:ncpu_in-1), gc_out%recv_off(0:ncpu_in-1))
    gc_out%send_cnt = 0; gc_out%recv_cnt = 0

    do ipass = 1, 2
      ! Initialize 3-plane window: prev=plane(N), curr=plane(1)
      idx_prev = 0; idx_curr = 1; idx_next = 2
      call compute_plane_owner(idd_type, N, N, planes(:,:,idx_prev))
      call compute_plane_owner(idd_type, 1, N, planes(:,:,idx_curr))

      do k = 1, N
        kn = modulo(k, N) + 1
        call compute_plane_owner(idd_type, kn, N, planes(:,:,idx_next))

        do j = 1, N
          do i = 1, N
            p = planes(i, j, idx_curr)
            if(p == myid_in) then
              ! My cell: find unique neighbor CPUs for ghost send
              n_nbr = 0
              do d = 1, 6
                q = nbr_owner(i, j, d, planes, N, idx_prev, idx_curr, idx_next)
                if(q == myid_in) cycle
                is_new = .true.
                do m = 1, n_nbr
                  if(nbr_cpu(m) == q) then; is_new = .false.; exit; end if
                end do
                if(is_new) then
                  n_nbr = n_nbr + 1; nbr_cpu(n_nbr) = q
                  if(ipass == 1) then
                    gc_out%send_cnt(q) = gc_out%send_cnt(q) + 1
                  else
                    stmp(q) = stmp(q) + 1
                    gc_out%send_cells(:, gc_out%send_off(q)+stmp(q)) = (/i,j,k/)
                  end if
                end if
              end do
            else
              ! Not my cell: check if any face neighbor is mine (ghost recv)
              do d = 1, 6
                q = nbr_owner(i, j, d, planes, N, idx_prev, idx_curr, idx_next)
                if(q == myid_in) then
                  if(ipass == 1) then
                    gc_out%recv_cnt(p) = gc_out%recv_cnt(p) + 1
                  else
                    rtmp(p) = rtmp(p) + 1
                    gc_out%recv_cells(:, gc_out%recv_off(p)+rtmp(p)) = (/i,j,k/)
                  end if
                  exit  ! count once per ghost cell
                end if
              end do
            end if
          end do
        end do

        ! Rotate plane indices (no data copy)
        itmp = idx_prev; idx_prev = idx_curr; idx_curr = idx_next; idx_next = itmp
      end do

      if(ipass == 1) then
        ! Compute offsets
        gc_out%send_off(0) = 0; gc_out%recv_off(0) = 0
        do p = 1, ncpu_in-1
          gc_out%send_off(p) = gc_out%send_off(p-1) + gc_out%send_cnt(p-1)
          gc_out%recv_off(p) = gc_out%recv_off(p-1) + gc_out%recv_cnt(p-1)
        end do
        gc_out%nsend_tot = gc_out%send_off(ncpu_in-1) + gc_out%send_cnt(ncpu_in-1)
        gc_out%nrecv_tot = gc_out%recv_off(ncpu_in-1) + gc_out%recv_cnt(ncpu_in-1)

        gc_out%n_partner = 0
        do p = 0, ncpu_in-1
          if(gc_out%send_cnt(p)>0 .or. gc_out%recv_cnt(p)>0) gc_out%n_partner = gc_out%n_partner+1
        end do

        allocate(gc_out%send_cells(3, max(1,gc_out%nsend_tot)))
        allocate(gc_out%recv_cells(3, max(1,gc_out%nrecv_tot)))
        allocate(stmp(0:ncpu_in-1), rtmp(0:ncpu_in-1))
        stmp = 0; rtmp = 0
      end if
    end do

    deallocate(planes, stmp, rtmp)
  end subroutine

  !-------------------------------------------------------------------
  ! Get owner CPU of the d-th face neighbor of cell (i,j) in current plane
  !-------------------------------------------------------------------
  function nbr_owner(i, j, d, planes, N, ip, ic, in) result(q)
    integer, intent(in) :: i, j, d, N, ip, ic, in
    integer, intent(in) :: planes(N, N, 0:2)
    integer :: q, ni, nj

    select case(d)
    case(1)  ! i+1
      ni = modulo(i, N) + 1
      q = planes(ni, j, ic)
    case(2)  ! i-1
      ni = modulo(i-2, N) + 1
      q = planes(ni, j, ic)
    case(3)  ! j+1
      nj = modulo(j, N) + 1
      q = planes(i, nj, ic)
    case(4)  ! j-1
      nj = modulo(j-2, N) + 1
      q = planes(i, nj, ic)
    case(5)  ! k+1
      q = planes(i, j, in)
    case(6)  ! k-1
      q = planes(i, j, ip)
    case default
      q = -1
    end select
  end function

  ! Dispatch to the correct exchange routine
  subroutine run_exchange(gc_in, iex_type, nerr_out)
    type(ghost_comm), intent(in) :: gc_in
    integer, intent(in) :: iex_type
    integer, intent(out) :: nerr_out

    if(iex_type == 1) then
      call ghost_p2p_exchange(gc_in, nvar, comm, nerr_out)
    else if(iex_type == 2) then
      call ghost_a2av_exchange(gc_in, nvar, comm, nerr_out)
    else
      call ghost_ks_exchange(gc_in, nvar, myid, ncpu, comm, nerr_out)
    end if
  end subroutine

  ! Compute stats: max across ranks, then mean+-std
  subroutine compute_stats(t, n, mean_out, std_out, comm2)
    real(8), intent(in) :: t(:)
    integer, intent(in) :: n, comm2
    real(8), intent(out) :: mean_out, std_out
    real(8) :: tmax(n), mn, s2
    integer :: ii, ierr2

    tmax(1:n) = t(1:n)
    call MPI_ALLREDUCE(MPI_IN_PLACE, tmax, n, MPI_DOUBLE_PRECISION, MPI_MAX, comm2, ierr2)
    mn = sum(tmax(1:n)) / dble(n)
    s2 = 0.0d0
    if(n > 1) s2 = sqrt(sum((tmax(1:n) - mn)**2) / dble(n-1))
    mean_out = mn; std_out = s2
  end subroutine

end program ghost_bench
