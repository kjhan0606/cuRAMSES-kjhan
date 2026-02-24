!@################################################
! (C)  Copyright 2023-29 Dr. Juhan Kim's Software 

!#define _KJHAN_TEST 

#ifdef _KJHAN_TEST
module AA
   integer, parameter:: dp=kind(1.0D0)
   logical:: vrel_merge=.false.
   real(dp),allocatable, dimension(:,:):: xsink
   real(dp),allocatable, dimension(:):: msink
   real(dp),allocatable, dimension(:,:):: vsink
   real(dp)::rmerge=1.0d0
end module AA
#endif
module Tree_commons
#ifdef _KJHAN_TEST
  use AA
#else
  use amr_commons
#endif
   real(dp):: MINCELLWIDTH = 1.e-6
   integer:: MIN_CELL_PARTICLE_NUM = 5
   TYPE Tree_Node
     integer:: sibling
    logical:: included
     real(dp):: x0,y0,z0,R
     integer:: daughter, nump
   end TYPE Tree_Node
   TYPE Tree_Particle
     integer:: sibling
     integer:: gid
     logical:: included
     real(dp):: x,y,z
   end TYPE Tree_Particle

end module Tree_commons
logical function Division(thisnode, pnum)
  use Tree_commons
  implicit none
  type(Tree_Node):: thisnode
  integer:: pnum
  if(thisnode%R .gt. 0.5*MINCELLWIDTH .and. pnum .ge. MIN_CELL_PARTICLE_NUM) then
     Division = .true.
  else
     Division = .false.
  endif
end function Division

subroutine EraseFromTree(node, particle,np, optr, ptr, nptr)
  use Tree_commons
  implicit none
  integer:: np
  TYPE(Tree_Node), dimension(1:np):: node
  TYPE(Tree_Particle),  dimension(1:np):: particle
   integer:: optr, ptr, nptr
  if(optr .gt. 0) then
     particle(optr)%sibling = nptr
  else
     if(node(-optr)%daughter .eq. ptr) then
        node(-optr)%daughter = nptr
     else
        node(-optr)%sibling = nptr
     endif
  endif

end subroutine EraseFromTree

logical function pfof_open(node, particle,np, inow, ptr, foflink, boxsize)
  use Tree_commons
  implicit none
  integer:: np
  TYPE(Tree_Node), dimension(1:np):: node
  TYPE(Tree_Particle),  dimension(1:np):: particle
  real(dp), dimension(1:3):: boxsize
  integer:: inow, ptr
  real(dp):: foflink,tmpx,tmpy,tmpz,dist,diff

! tmpx = abs(particle(inow)%x - node(-ptr)%x0)
! tmpy = abs(particle(inow)%y - node(-ptr)%y0)
! tmpz = abs(particle(inow)%z - node(-ptr)%z0)
! if(tmpx .ge. boxsize(1)*0.5) tmpx = boxsize(1) - tmpx
! if(tmpy .ge. boxsize(2)*0.5) tmpy = boxsize(2) - tmpy
! if(tmpz .ge. boxsize(3)*0.5) tmpz = boxsize(3) - tmpz
  tmpx = (particle(inow)%x - node(-ptr)%x0)
  tmpy = (particle(inow)%y - node(-ptr)%y0)
  tmpz = (particle(inow)%z - node(-ptr)%z0)
  if(tmpx .ge. boxsize(1)*0.5) then
    tmpx = tmpx - boxsize(1)
  endif
  if(tmpx .le. -boxsize(1)*0.5) then
    tmpx = tmpx + boxsize(1)
  endif
  if(tmpy .ge. boxsize(2)*0.5)then
    tmpy = tmpy - boxsize(2)
  endif
  if(tmpy .le. -boxsize(2)*0.5)then
    tmpy = tmpy + boxsize(2)
  endif
  if(tmpz .ge. boxsize(3)*0.5)then
    tmpz = tmpz - boxsize(3)
  endif
  if(tmpz .le. -boxsize(3)*0.5)then
    tmpz = tmpz + boxsize(3)
  endif

  dist = dsqrt(tmpx*tmpx + tmpy*tmpy + tmpz*tmpz)
  diff = dist - node(-ptr)%R
  pfof_open = diff <= foflink
end function pfof_open


integer function Omp_FoF_link(gndx, node, particle, np, foflink , linked, gindx, boxsize, factG,included,gid) result(ncount)
#ifdef _KJHAN_TEST
  use AA
#else
  use pm_commons
  use amr_commons
#endif
  use Tree_commons
  implicit none
  logical, dimension(1:np):: included
  integer, dimension(1:np):: gid
  integer:: i,j,k,ptr,optr, nptr,now,gndx,gindx,np
  TYPE(Tree_Node), dimension(1:np), intent(inout):: node
  TYPE(Tree_Particle),  dimension(1:np), intent(inout):: particle
  integer, dimension(1:np), intent(out):: linked
  real(dp), dimension(3):: boxsize
  real(dp):: foflink2,tmpx,tmpy,tmpz, foflink,rr
  real(dp):: egrav,ekin,uxcom,uycom,uzcom,v2rel1,v2rel2,factG
  logical, external:: pfof_open

  foflink2 = foflink*foflink
  ncount = 0
  now = 0

  do while(.true.)
     optr =  -1
     ptr = -1;
     do while(ptr .ne. 0)
        if(ptr.lt. 0) then
           if(node(-ptr)%sibling .eq. node(-ptr)%daughter) then
!             call EraseFromTree(node, particle,np, optr,ptr, node(-ptr)%sibling)
              ptr = node(-ptr)%sibling
           else
              if(pfof_open(node, particle,np,gndx, ptr, foflink, boxsize).eqv. .true.) then
                 optr = ptr
                 ptr = node(-ptr)%daughter
              else
                 optr = ptr;
                 ptr = node(-ptr)%sibling
              endif
           endif
        else if(ptr .gt.0) then
!          if(particle(ptr)%included .eqv. .true.) then
           if(included(ptr) .eqv. .true.) then
              nptr = particle(ptr)%sibling
!             call EraseFromTree(node, particle,np, optr, ptr, nptr)
              ptr = nptr
           else
!             tmpx = abs(particle(gndx)%x - particle(ptr)%x)
!             tmpy = abs(particle(gndx)%y - particle(ptr)%y)
!             tmpz = abs(particle(gndx)%z - particle(ptr)%z)
!             if(tmpx.ge.boxsize(1)*0.5) tmpx = boxsize(1) - tmpx
!             if(tmpy.ge.boxsize(2)*0.5) tmpy = boxsize(2) - tmpy
!             if(tmpz.ge.boxsize(3)*0.5) tmpz = boxsize(3) - tmpz
              tmpx = (particle(gndx)%x - particle(ptr)%x)
              tmpy = (particle(gndx)%y - particle(ptr)%y)
              tmpz = (particle(gndx)%z - particle(ptr)%z)
              if(tmpx.ge.boxsize(1)*0.5) then
                tmpx = tmpx - boxsize(1)
              endif
              if(tmpx.le.-boxsize(1)*0.5) then
                tmpx = tmpx + boxsize(1)
              endif
              rr = tmpx*tmpx
#if NDIM>1
              if(tmpy.ge.boxsize(2)*0.5) then
                tmpy = tmpy - boxsize(2)
              endif
              if(tmpy.le.-boxsize(2)*0.5) then
                tmpy = tmpy + boxsize(2)
              endif
              rr = rr+ tmpy*tmpy
#endif
#if NDIM>2
              if(tmpz.ge.boxsize(3)*0.5) then
                tmpz = tmpz - boxsize(3)
              endif
              if(tmpz.le.-boxsize(3)*0.5) then
                tmpz = tmpz + boxsize(3)
              endif
              rr = rr+ tmpz*tmpz
#endif
              if(rr .le. foflink2) then
#ifdef _KJHAN_TEST
                    ncount = ncount + 1
                    linked(ncount) = ptr
                    gid(ptr) = gindx
                    included(ptr) = .true.
                    nptr = particle(ptr)%sibling
!                   call EraseFromTree(node, particle,np, optr, ptr, nptr)
#else
                 if(vrel_merge) then
                    egrav=msink(ptr)*msink(gndx)/(rr+tiny(0.d0))*factG
                    uxcom=(msink(ptr)*vsink(ptr,1)+msink(gndx)*vsink(gndx,1))/(msink(ptr)+msink(gndx))
                    uycom=(msink(ptr)*vsink(ptr,2)+msink(gndx)*vsink(gndx,2))/(msink(ptr)+msink(gndx))
                    uzcom=(msink(ptr)*vsink(ptr,3)+msink(gndx)*vsink(gndx,3))/(msink(ptr)+msink(gndx))
                    v2rel1=(vsink(ptr,1)-uxcom)**2+(vsink(ptr,2)-uycom)**2+(vsink(ptr,3)-uzcom)**2
                    v2rel2=(vsink(gndx,1)-uxcom)**2+(vsink(gndx,2)-uycom)**2+(vsink(gndx,3)-uzcom)**2
                    ekin=0.5d0*(msink(ptr)*v2rel1+msink(gndx)*v2rel2)
                    if(ekin .lt. egrav) then
                       ncount = ncount + 1
                       linked(ncount) = ptr
                       gid(ptr) = gindx
                       included(ptr) = .true.
                       nptr = particle(ptr)%sibling
!                      call EraseFromTree(node, particle,np, optr, ptr, nptr)
                    else
                       optr = ptr
                    endif
                 else
                    ncount = ncount + 1
                    linked(ncount) = ptr
                    gid(ptr) = gindx
                    included(ptr) = .true.
                    nptr = particle(ptr)%sibling
!                   call EraseFromTree(node, particle,np, optr, ptr, nptr)
                 endif
#endif
              else
                 optr = ptr
              endif
              ptr = particle(ptr)%sibling
           endif
        endif
     enddo
     now = now + 1
     if(now .gt. ncount) exit
     gndx = linked(now)
  enddo
end function Omp_FoF_link


integer function Omp_Do_FoF_Tree(numsink, psink, gsink, dx_min2, factG) result(gindx)
#ifdef _OPENMP
  use OMP_LIB
#endif
#ifdef _KJHAN_TEST
  use AA
#else
  use pm_commons
  use amr_commons
#endif
  use Tree_commons
  implicit none
  include 'mpif.h'
  integer:: i,j,k, indx,numsink,kindx
  TYPE(Tree_Node),  dimension(:), allocatable:: node
  TYPE(Tree_Particle),  dimension(:), allocatable:: particle
  INTEGER mpistatus(MPI_STATUS_SIZE)
  integer,dimension(1:numsink)::gsink
  integer, dimension(:), allocatable:: linked
  logical, dimension(:), allocatable:: included
  integer, dimension(:), allocatable:: gid
  integer,dimension(1:numsink)::psink
  integer, target, allocatable:: Pgindx(:),Picount(:)
  integer, target, allocatable:: Pptmp(:,:),Pgtmp(:,:)
  integer, pointer ::ptmp(:),gtmp(:)
  logical:: crossflag
  integer:: gnum,ii,icount,ticount,ipos,jcount,gjndx
  real(dp), dimension(3):: boxsize
  real(dp):: foflink
  real(dp):: dx_min2
  real(dp):: factG, scale
  real(dp),dimension(1:3)::xbound,skip_loc
  integer, external:: Omp_FoF_link
  integer:: nx_loc
  integer:: nid, myrank, nthreads, mythread, ierror
  integer:: ssink,fsink, issink, ifsink, mynumsink,njobs
#ifndef _OPENMP
  integer, external:: OMP_get_num_threads,OMP_get_thread_num
#endif

#ifdef _KJHAN_TEST
  xbound(1:3) = (/1,1,1/)
  scale = 1
#else
  xbound(1:3)=(/dble(nx),dble(ny),dble(nz)/)
  nx_loc=(icoarse_max-icoarse_min+1)
  scale=boxlen/dble(nx_loc)
#endif



  allocate(node(1:2*numsink))
  allocate(particle(1:numsink))

  foflink = dsqrt(rmerge**2*dx_min2)
  boxsize(1) = scale*xbound(1)
  boxsize(2) = scale*xbound(2)
  boxsize(3) = scale*xbound(3)
  do i = 1, numsink
     indx = psink(i)
     particle(i)%x = xsink(indx,1)
     particle(i)%y = xsink(indx,2)
     particle(i)%z = xsink(indx,3)
!    particle(i)%included = .false.
  enddo
  call Build_Tree(node, particle, numsink)

! inserted for the MPI+OpenMP implementatio for MPI+OpenMP
  call MPI_Comm_size(MPI_COMM_WORLD, nid, ierror)
  call MPI_Comm_rank(MPI_COMM_WORLD, myrank, ierror)


  mynumsink = (numsink+nid-1)/nid
  ssink = mynumsink*myrank  + 1
  fsink = ssink + mynumsink - 1
  fsink = MIN(fsink, numsink)
  mynumsink = fsink - ssink + 1

! inserted for the MPI+OpenMP implementatio for MPI+OpenMP
!$OMP PARALLEL 
!$OMP MASTER 
   nthreads = OMP_get_num_threads()
!$OMP END MASTER 
!$OMP END PARALLEL

   allocate(Pptmp(1:numsink, 0:nthreads-1))
   allocate(Pgtmp(1:numsink, 0:nthreads-1))
   allocate(Pgindx(0:nthreads-1))
   allocate(Picount(0:nthreads-1))

!$OMP PARALLEL PRIVATE(mythread,gindx,icount, i, ii,linked,ptmp,gtmp,included,gid,issink, ifsink,njobs,jcount, gjndx,crossflag,gnum)
   mythread = OMP_get_thread_num()
   allocate(included(1:numsink))
   allocate(gid(1:numsink))
   allocate(linked(1:numsink))
   gtmp => Pgtmp(:, mythread)
   ptmp => Pptmp(:, mythread)

  do i = 1, numsink
     included(i) = .false.
  enddo

  gindx = 0
  icount = 0
  njobs = (mynumsink+nthreads-1)/nthreads

  issink = njobs*mythread + 1
  ifsink = issink + njobs - 1
  ifsink = MIN(ifsink, mynumsink)
  njobs = ifsink - issink + 1
! print *, myid, mythread, njobs, ssink,ssink-1+issink, ssink-1+ifsink
  issink = issink + ssink -1
  ifsink = ifsink + ssink - 1
  do i = issink, ifsink
     if(included(i) .eqv. .false.) then
        ii = i
        gindx = gindx + 1
        gnum = Omp_FoF_link(ii, node, particle, numsink, foflink, linked, gindx, boxsize, factG,included,gid)
        crossflag = .false.
        do j = 1, gnum
           if(linked(j) .lt. issink) then
              crossflag = .true.
           endif
        enddo
        if(crossflag .eqv. .false.) then
           do j = 1, gnum
              ptmp(icount+j) = linked(j)
              gtmp(icount+j) = gindx
!             psink(icount+j) = linked(j)
!             gsink(linked(j)) = gindx
           enddo
           icount = icount + gnum
        else
           do j = 1, gnum
              gid(linked(j)) = -1
              included(linked(j)) = .false.
           enddo
           gindx = gindx - 1
        endif
     endif
  enddo
  deallocate(included,gid,linked)
  Pgindx(mythread) = gindx
  Picount(mythread) = icount
!#ifdef _KJHAN_TEST
!  print *, myid, mythread,gindx, icount,issink, ifsink
!#endif
!$OMP barrier
  jcount = 0
  gjndx = 0
  do i = 0, mythread-1
     jcount = jcount + Picount(i)
     gjndx = gjndx + Pgindx(i)
  enddo
  do i = 1, icount
     gtmp(i) = gtmp(i) + gjndx
  enddo
!$OMP barrier
!$OMP MASTER
  ipos = icount
  do i = 1, nthreads-1
     do j = 1, Picount(i)
        ipos = ipos + 1
        gtmp(ipos) = Pgtmp(j, i)
        ptmp(ipos) = Pptmp(j, i)
     enddo
  enddo
  ticount = ipos
!$OMP END MASTER

!$OMP END PARALLEL
  if(ticount > 0) then
     gindx = Pgtmp(ticount,0)
  else
     gindx = 0
  endif
  icount = ticount
  ptmp => Pptmp(:,0)
  gtmp => Pgtmp(:,0)


  kindx = 0
  if(myrank .ne.0) then
     call MPI_Recv(kindx, 1, MPI_INTEGER, myrank-1, 0, MPI_COMM_WORLD, mpistatus, ierror)
  endif

  do i = 1, icount
     gtmp(i) = gtmp(i) + kindx
  enddo
  if(icount > 0) then
     kindx = gtmp(icount)
  else
     kindx = kindx + gindx
  endif
  if(myrank.ne.nid-1) then
     call MPI_Send(kindx, 1, MPI_INTEGER, myrank+1, 0,MPI_COMM_WORLD, mpistatus, ierror)
  endif

  if(myrank.eq.0) then
     do i = 1, nid-1
        call MPI_Recv(jcount, 1, MPI_INTEGER, i, 0, MPI_COMM_WORLD, mpistatus, ierror)
        call MPI_Recv(ptmp(icount+1), jcount, MPI_INTEGER, i, 1, MPI_COMM_WORLD, mpistatus, ierror)
        call MPI_Recv(gtmp(icount+1), jcount, MPI_INTEGER, i, 2, MPI_COMM_WORLD, mpistatus, ierror)
        icount = icount + jcount
     enddo
     do i = 1, icount
       psink(i) = ptmp(i)
       gsink(psink(i)) = gtmp(i)
     enddo
  else
     call MPI_Send(icount, 1, MPI_INTEGER, 0, 0, MPI_COMM_WORLD, mpistatus, ierror)
     call MPI_Send(ptmp, icount, MPI_INTEGER, 0, 1, MPI_COMM_WORLD, mpistatus, ierror)
     call MPI_Send(gtmp, icount, MPI_INTEGER, 0, 2, MPI_COMM_WORLD, mpistatus, ierror)
  endif
  gindx = kindx
  call MPI_Bcast(gindx, 1,MPI_INTEGER, nid-1, MPI_COMM_WORLD,ierror)
  call MPI_Bcast(psink, numsink,MPI_INTEGER, 0, MPI_COMM_WORLD,ierror)
  call MPI_Bcast(gsink, numsink,MPI_INTEGER, 0, MPI_COMM_WORLD,ierror)

  deallocate(Pgtmp,Pptmp,Pgindx,Picount)
  deallocate(node, particle)

end function Omp_Do_FoF_Tree

integer function FoF_link(gndx, node, particle, np, foflink , linked, gindx, boxsize, factG) result(ncount)
#ifdef _KJHAN_TEST
  use AA
#else
  use pm_commons
  use amr_commons
#endif
  use Tree_commons
  implicit none
  integer:: i,j,k,ptr,optr, nptr,now,gndx,gindx,np
  TYPE(Tree_Node), dimension(1:np), intent(inout):: node
  TYPE(Tree_Particle),  dimension(1:np), intent(inout):: particle
  integer, dimension(1:np), intent(out):: linked
  real(dp), dimension(3):: boxsize
  real(dp):: foflink2,tmpx,tmpy,tmpz, foflink,rr
  real(dp):: egrav,ekin,uxcom,uycom,uzcom,v2rel1,v2rel2,factG
  logical:: pfof_open

  foflink2 = foflink*foflink
  ncount = 0
  now = 0

  do while(.true.)
     optr =  -1
     ptr = -1;
     do while(ptr .ne. 0)
        if(ptr.lt. 0) then
           if(node(-ptr)%sibling .eq. node(-ptr)%daughter) then
              call EraseFromTree(node, particle,np, optr,ptr, node(-ptr)%sibling)
              ptr = node(-ptr)%sibling
           else
              if(pfof_open(node, particle,np,gndx, ptr, foflink, boxsize).eqv. .true.) then
                 optr = ptr
                 ptr = node(-ptr)%daughter
              else
                 optr = ptr;
                 ptr = node(-ptr)%sibling
              endif
           endif
        else if(ptr .gt.0) then
           if(particle(ptr)%included .eqv. .true.) then
              nptr = particle(ptr)%sibling
              call EraseFromTree(node, particle,np, optr, ptr, nptr)
              ptr = nptr
           else
!             tmpx = abs(particle(gndx)%x - particle(ptr)%x)
!             tmpy = abs(particle(gndx)%y - particle(ptr)%y)
!             tmpz = abs(particle(gndx)%z - particle(ptr)%z)
!             if(tmpx.ge.boxsize(1)*0.5) tmpx = boxsize(1) - tmpx
!             if(tmpy.ge.boxsize(2)*0.5) tmpy = boxsize(2) - tmpy
!             if(tmpz.ge.boxsize(3)*0.5) tmpz = boxsize(3) - tmpz
              tmpx = (particle(gndx)%x - particle(ptr)%x)
              tmpy = (particle(gndx)%y - particle(ptr)%y)
              tmpz = (particle(gndx)%z - particle(ptr)%z)
              if(tmpx.ge.boxsize(1)*0.5) then
                tmpx = tmpx - boxsize(1)
              endif
              if(tmpx.le.-boxsize(1)*0.5) then
                tmpx = tmpx + boxsize(1)
              endif
              rr = tmpx*tmpx
#if NDIM>1
              if(tmpy.ge.boxsize(2)*0.5) then
                tmpy = tmpy - boxsize(2)
              endif
              if(tmpy.le.-boxsize(2)*0.5) then
                tmpy = tmpy + boxsize(2)
              endif
              rr = rr+ tmpy*tmpy
#endif
#if NDIM>1
              if(tmpz.ge.boxsize(3)*0.5) then
                tmpz = tmpz - boxsize(3)
              endif
              if(tmpz.le.-boxsize(3)*0.5) then
                tmpz = tmpz + boxsize(3)
              endif
              rr = rr+ tmpz*tmpz
#endif
              if(rr .le. foflink2) then
                 if(vrel_merge) then
                    egrav=msink(ptr)*msink(gndx)/(rr+tiny(0.d0))*factG
                    uxcom=(msink(ptr)*vsink(ptr,1)+msink(gndx)*vsink(gndx,1))/(msink(ptr)+msink(gndx))
                    uycom=(msink(ptr)*vsink(ptr,2)+msink(gndx)*vsink(gndx,2))/(msink(ptr)+msink(gndx))
                    uzcom=(msink(ptr)*vsink(ptr,3)+msink(gndx)*vsink(gndx,3))/(msink(ptr)+msink(gndx))
                    v2rel1=(vsink(ptr,1)-uxcom)**2+(vsink(ptr,2)-uycom)**2+(vsink(ptr,3)-uzcom)**2
                    v2rel2=(vsink(gndx,1)-uxcom)**2+(vsink(gndx,2)-uycom)**2+(vsink(gndx,3)-uzcom)**2
                    ekin=0.5d0*(msink(ptr)*v2rel1+msink(gndx)*v2rel2)
                    if(ekin .lt. egrav) then
                       ncount = ncount + 1
                       linked(ncount) = ptr
                       particle(ptr)%gid = gindx
                       particle(ptr)%included = .true.
                       nptr = particle(ptr)%sibling
                       call EraseFromTree(node, particle,np, optr, ptr, nptr)
                    else
                       optr = ptr
                    endif
                 else
                    ncount = ncount + 1
                    linked(ncount) = ptr
                    particle(ptr)%gid = gindx
                    particle(ptr)%included = .true.
                    nptr = particle(ptr)%sibling
                    call EraseFromTree(node, particle,np, optr, ptr, nptr)
                 endif
              else
                 optr = ptr
              endif
              ptr = particle(ptr)%sibling
           endif
        endif
     enddo
     now = now + 1
     if(now .gt. ncount) exit
     gndx = linked(now)
  enddo
end function FoF_link


function Auto_Do_FoF(numsink, psink, gsink, dx_min2, factG) result(gindx)
! Auto-tuning FoF: picks fastest among Sequential, Tree, OMP+MPI Tree
! Phase 0: run all 3, pick winner
! Phase 1+: use winner
#ifdef _KJHAN_TEST
  use AA
#else
  use pm_commons
  use amr_commons
#endif
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer, intent(in):: numsink
  integer, dimension(1:numsink), intent(inout):: psink, gsink
  real(dp), intent(in):: dx_min2, factG
  integer:: gindx

  integer, external:: Do_FoF_Sequential, Do_FoF_Tree, Omp_Do_FoF_Tree

  ! Auto-tune state (saved between calls)
  integer, save:: fof_phase = 0    ! 0=probe, 1=decided
  integer, save:: fof_winner = 0   ! 1=seq, 2=tree, 3=omp_tree
  integer, save:: ncalls = 0

  integer:: i, ierr, myid_loc
  integer, allocatable:: psink_tmp(:), gsink_tmp(:)
  real(dp):: t1, t2, t_seq, t_tree, t_omp
  real(dp):: t_min

  myid_loc = 1
#ifndef WITHOUTMPI
  call MPI_COMM_RANK(MPI_COMM_WORLD, myid_loc, ierr)
  myid_loc = myid_loc + 1
#endif

  ncalls = ncalls + 1

  if(fof_phase == 0) then
    ! Phase 0: probe all 3 methods
    allocate(psink_tmp(1:numsink), gsink_tmp(1:numsink))

    ! --- Method 1: Sequential ---
    do i = 1, numsink
      psink_tmp(i) = psink(i)
      gsink_tmp(i) = gsink(i)
    enddo
#ifndef WITHOUTMPI
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    t1 = MPI_WTIME()
#endif
    gindx = Do_FoF_Sequential(numsink, psink_tmp, gsink_tmp, dx_min2, factG)
#ifndef WITHOUTMPI
    t2 = MPI_WTIME()
    t_seq = t2 - t1
#else
    t_seq = 1.0d10
#endif

    ! --- Method 2: Tree ---
    do i = 1, numsink
      psink_tmp(i) = psink(i)
      gsink_tmp(i) = gsink(i)
    enddo
#ifndef WITHOUTMPI
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    t1 = MPI_WTIME()
#endif
    gindx = Do_FoF_Tree(numsink, psink_tmp, gsink_tmp, dx_min2, factG)
#ifndef WITHOUTMPI
    t2 = MPI_WTIME()
    t_tree = t2 - t1
#else
    t_tree = 1.0d10
#endif

    ! --- Method 3: OMP+MPI Tree ---
    do i = 1, numsink
      psink_tmp(i) = psink(i)
      gsink_tmp(i) = gsink(i)
    enddo
#ifndef WITHOUTMPI
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    t1 = MPI_WTIME()
#endif
    gindx = Omp_Do_FoF_Tree(numsink, psink_tmp, gsink_tmp, dx_min2, factG)
#ifndef WITHOUTMPI
    t2 = MPI_WTIME()
    t_omp = t2 - t1
#else
    t_omp = 1.0d10
#endif

    deallocate(psink_tmp, gsink_tmp)

    ! Pick winner (all ranks agree since timings are global)
    t_min = min(t_seq, t_tree, t_omp)
    if(t_min == t_seq) then
      fof_winner = 1
    else if(t_min == t_tree) then
      fof_winner = 2
    else
      fof_winner = 3
    endif

#ifndef WITHOUTMPI
    ! Ensure all ranks agree
    call MPI_BCAST(fof_winner, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
#endif

    fof_phase = 1
    if(myid_loc == 1) then
      write(*,'(A,F8.3,A,F8.3,A,F8.3,A)') &
        ' [FoF auto-tune] seq=', t_seq*1000, 'ms tree=', t_tree*1000, &
        'ms omp=', t_omp*1000, 'ms'
      if(fof_winner == 1) write(*,'(A)') ' [FoF auto-tune] → Sequential'
      if(fof_winner == 2) write(*,'(A)') ' [FoF auto-tune] → Tree (serial)'
      if(fof_winner == 3) write(*,'(A)') ' [FoF auto-tune] → OMP+MPI Tree'
    endif
  endif

  ! Phase 1: run winner
  if(fof_winner == 1) then
    gindx = Do_FoF_Sequential(numsink, psink, gsink, dx_min2, factG)
  else if(fof_winner == 2) then
    gindx = Do_FoF_Tree(numsink, psink, gsink, dx_min2, factG)
  else
    gindx = Omp_Do_FoF_Tree(numsink, psink, gsink, dx_min2, factG)
  endif

end function Auto_Do_FoF


function Do_FoF_Sequential(numsink, psink, gsink, dx_min2, factG) result(gindx)
! Brute-force O(n^2) FoF — same as sink_particle NO_OMP path
#ifdef _KJHAN_TEST
  use AA
#else
  use pm_commons
  use amr_commons
#endif
  implicit none
  integer, intent(in):: numsink
  integer, dimension(1:numsink), intent(inout):: psink
  integer, dimension(1:numsink), intent(inout):: gsink
  real(dp), intent(in):: dx_min2, factG
  integer:: gindx

  integer:: icomp, ifirst, ilast, gndx, indx
  real(dp):: xx, yy, zz, rr, foflink, foflink2, scale
  real(dp), dimension(1:3):: boxsize, xbound
  integer:: nx_loc

#ifdef _KJHAN_TEST
  xbound(1:3) = (/1,1,1/)
  scale = 1
#else
  xbound(1:3)=(/dble(nx),dble(ny),dble(nz)/)
  nx_loc=(icoarse_max-icoarse_min+1)
  scale=boxlen/dble(nx_loc)
#endif

  foflink = dsqrt(rmerge**2 * dx_min2)
  foflink2 = foflink * foflink
  boxsize(1) = scale*xbound(1)
  boxsize(2) = scale*xbound(2)
  boxsize(3) = scale*xbound(3)

  gindx = 0
  icomp = 1
  ifirst = 2

  do while(icomp <= numsink)
    gndx = psink(icomp)
    if(gsink(gndx) == 0) then
      gindx = gindx + 1
      gsink(gndx) = gindx
    endif
    ilast = numsink
    do while((ilast - ifirst + 1) > 0)
      indx = psink(ifirst)
      xx = xsink(indx,1) - xsink(gndx,1)
      if(xx >= boxsize(1)*0.5d0) xx = xx - boxsize(1)
      if(xx <= -boxsize(1)*0.5d0) xx = xx + boxsize(1)
      rr = xx*xx
#if NDIM>1
      yy = xsink(indx,2) - xsink(gndx,2)
      if(yy >= boxsize(2)*0.5d0) yy = yy - boxsize(2)
      if(yy <= -boxsize(2)*0.5d0) yy = yy + boxsize(2)
      rr = rr + yy*yy
#endif
#if NDIM>2
      zz = xsink(indx,3) - xsink(gndx,3)
      if(zz >= boxsize(3)*0.5d0) zz = zz - boxsize(3)
      if(zz <= -boxsize(3)*0.5d0) zz = zz + boxsize(3)
      rr = rr + zz*zz
#endif
      if(rr <= foflink2) then
#ifndef _KJHAN_TEST
        if(vrel_merge) then
          call check_vrel_merge(indx, gndx, rr, factG, ifirst, ilast, gindx, psink, gsink, numsink)
        else
#endif
          ifirst = ifirst + 1
          gsink(indx) = gindx
#ifndef _KJHAN_TEST
        endif
#endif
      else
        psink(ifirst) = psink(ilast)
        psink(ilast) = indx
        ilast = ilast - 1
      endif
    enddo
    icomp = icomp + 1
  enddo

contains
  subroutine check_vrel_merge(indx_in, gndx_in, rr_in, factG_in, &
       ifirst_io, ilast_io, igrp_in, psink_io, gsink_io, ns)
#ifndef _KJHAN_TEST
    use pm_commons
#endif
    implicit none
    integer, intent(in):: indx_in, gndx_in, igrp_in, ns
    real(dp), intent(in):: rr_in, factG_in
    integer, intent(inout):: ifirst_io, ilast_io
    integer, intent(inout):: psink_io(ns), gsink_io(ns)
    real(dp):: egrav, ekin, uxcom, uycom, uzcom, v2rel1, v2rel2
#ifndef _KJHAN_TEST
    egrav = msink(indx_in)*msink(gndx_in)/(rr_in+tiny(0.d0))*factG_in
    uxcom = (msink(indx_in)*vsink(indx_in,1)+msink(gndx_in)*vsink(gndx_in,1)) &
         / (msink(indx_in)+msink(gndx_in))
    uycom = (msink(indx_in)*vsink(indx_in,2)+msink(gndx_in)*vsink(gndx_in,2)) &
         / (msink(indx_in)+msink(gndx_in))
    uzcom = (msink(indx_in)*vsink(indx_in,3)+msink(gndx_in)*vsink(gndx_in,3)) &
         / (msink(indx_in)+msink(gndx_in))
    v2rel1 = (vsink(indx_in,1)-uxcom)**2+(vsink(indx_in,2)-uycom)**2+(vsink(indx_in,3)-uzcom)**2
    v2rel2 = (vsink(gndx_in,1)-uxcom)**2+(vsink(gndx_in,2)-uycom)**2+(vsink(gndx_in,3)-uzcom)**2
    ekin = 0.5d0*(msink(indx_in)*v2rel1+msink(gndx_in)*v2rel2)
    if(ekin < egrav) then
      ifirst_io = ifirst_io + 1
      gsink_io(indx_in) = igrp_in
    else
      psink_io(ifirst_io) = psink_io(ilast_io)
      psink_io(ilast_io) = indx_in
      ilast_io = ilast_io - 1
    endif
#endif
  end subroutine
end function Do_FoF_Sequential

function Do_FoF_Tree(numsink, psink, gsink, dx_min2, factG) result(gindx)
#ifdef _KJHAN_TEST
  use AA
#else
  use pm_commons
  use amr_commons
#endif
  use Tree_commons
  implicit none
  integer:: i,j,k, indx,numsink, gindx
! TYPE(Tree_Node),  dimension(1:2*numsink):: node
! TYPE(Tree_Particle),  dimension(1:numsink):: particle
! integer, dimension(1:numsink):: linked

  TYPE(Tree_Node),  dimension(:), allocatable:: node
  TYPE(Tree_Particle),  dimension(:),allocatable:: particle
  integer, dimension(:),allocatable:: linked



  integer, dimension(1:numsink)::psink,gsink
  integer:: gnum,ii,icount
  real(dp), dimension(1:3):: boxsize
  real(dp):: foflink
  real(dp):: dx_min2
  real(dp):: factG, scale
  real(dp),dimension(1:3)::xbound,skip_loc
  integer, external:: FoF_link
  integer:: nx_loc


  allocate(node(1:2*numsink))
  allocate(particle(1:numsink))
  allocate(linked(1:numsink))

#ifdef _KJHAN_TEST
    xbound(1:3) = (/1,1,1/)
    scale = 1
#else
    xbound(1:3)=(/dble(nx),dble(ny),dble(nz)/)
    nx_loc=(icoarse_max-icoarse_min+1)
    scale=boxlen/dble(nx_loc)
#endif


  foflink = sqrt(rmerge**2*dx_min2)
  MINCELLWIDTH = 5*foflink
  boxsize(1) = scale*xbound(1)
  boxsize(2) = scale*xbound(2)
  boxsize(3) = scale*xbound(3)
  do i = 1, numsink
     particle(i)%x = xsink(i,1)
     particle(i)%y = xsink(i,2)
     particle(i)%z = xsink(i,3)
     particle(i)%included = .false.
  enddo
  call Build_Tree(node, particle, numsink)

  gindx = 0
  icount = 0
  do i = 1, numsink
     if(particle(i)%included .eqv. .false.) then
        ii = i
        gindx = gindx + 1
        gnum = FoF_link(ii, node, particle, numsink, foflink, linked, gindx, boxsize, factG)
!       if(gnum .eq.0 .or. particle(i)%included .ne. .true.) then
!          print *, 'Missing finding', gnum, particle(i)%included
!          stop
!       endif
        do j = 1, gnum
!          if(linked(j) .le.0 .or. linked(j).gt.numsink) then
!             print *,'wrong detection of linked', j, linked(j)
!             stop
!          endif
           psink(icount+j) = linked(j)
           gsink(linked(j)) = gindx
        enddo
        icount = icount + gnum
     endif
  enddo
! if(icount .ne. numsink) then
!    print *, 'Missing particles are detected ', icount, numsink
! endif
  deallocate(node,particle,linked)
end function Do_FoF_Tree

subroutine Build_Tree(node, particle, np)
  use Tree_commons
  implicit none
  integer:: np,i,j,k,iSpareNode
  TYPE(Tree_Node), dimension(1:np):: node
  TYPE(Tree_Particle), dimension(1:np):: particle
  interface
     recursive function Make_Tree(node,particle,np,iworknode,iSpareNode0) result(iNextSpareNode)
       use Tree_commons
       integer:: np
       TYPE(Tree_Node), target, dimension(1:np):: node
       TYPE(Tree_Particle), dimension(1:np):: particle
       integer:: iNextSpareNode
       integer, intent(in):: iSpareNode0, iworknode
     end function Make_Tree
  end interface
  integer:: iworknode, iSpareNode0


  do i = 1, np
     particle(i)%sibling = i + 1
     node(i)%sibling = 0
     node(i)%daughter = 0
  enddo
  particle(np)%sibling = 0
  node(1)%daughter = 1
  node(1)%nump = np

  iworknode = -1
  iSpareNode0 = -2

  iSpareNode = Make_Tree(node, particle,np, iworknode, iSpareNode0)

end subroutine Build_Tree

recursive function Make_Tree(node,particle,np,iworknode, iSpareNode0) result(iNextSpareNode)
  use Tree_commons
  implicit none
  integer:: np
  TYPE(Tree_Node),  target, dimension(1:np):: node
  TYPE(Tree_Node), pointer :: worknode
  TYPE(Tree_Particle),  dimension(1:np):: particle
  TYPE(Tree_Node),   dimension(0:7):: tmpnode
  integer:: iSpareNode,iNextSpareNode
  integer, intent(in)::iSpareNode0, iworknode
  integer:: i,j,k,ix,iy,iz,inextjobnode
  integer:: ntmp,mx,my,mz,mm
  integer:: iFirstDaughter, iNowDaughter,iLeftSibling
  real(dp):: tmpx,tmpy,tmpz,d2,maxd2
  logical, external:: Division


  iSpareNode = iSpareNode0

  worknode => node(-iworknode)

  worknode%x0 = 0
  worknode%y0 = 0
  worknode%z0 = 0
  ntmp = worknode%daughter
  worknode%nump = 0
  do while(ntmp .gt.0)
     worknode%x0 = worknode%x0 + particle(ntmp)%x
     worknode%y0 = worknode%y0 + particle(ntmp)%y
     worknode%z0 = worknode%z0 + particle(ntmp)%z
     worknode%nump = worknode%nump + 1
     ntmp = particle(ntmp)%sibling
  enddo
  worknode%x0 = worknode%x0/worknode%nump
  worknode%y0 = worknode%y0/worknode%nump
  worknode%z0 = worknode%z0/worknode%nump
  maxd2 = -1E20
  ntmp = worknode%daughter
  do while(ntmp .gt. 0)
     tmpx = particle(ntmp)%x - worknode%x0
     tmpy = particle(ntmp)%y - worknode%y0
     tmpz = particle(ntmp)%z - worknode%z0
     d2 = tmpx*tmpx + tmpy*tmpy + tmpz*tmpz
     maxd2 = MAX(maxd2, d2)
     ntmp = particle(ntmp)%sibling
  enddo
  worknode%R = dsqrt(maxd2)
  do i = 0, 7
     tmpnode(i)%sibling = 0
     tmpnode(i)%daughter = 0
     tmpnode(i)%nump = 0
  enddo
  ntmp = worknode%daughter
  do i = 1, worknode%nump
     if(particle(ntmp)%x .ge. worknode%x0) then
        mx = 1
     else
        mx = 0
     endif
     if(particle(ntmp)%y .ge. worknode%y0) then
        my = 1
     else
        my = 0
     endif
     if(particle(ntmp)%z .ge. worknode%z0) then
        mz = 1
     else
        mz = 0
     endif
     mm = mx+2*(my+2*mz)
     tmpnode(mm)%nump = tmpnode(mm)%nump + 1
     j = tmpnode(mm)%daughter
     k = particle(ntmp)%sibling
     tmpnode(mm)%daughter = ntmp
     particle(ntmp)%sibling = j
     ntmp = k
  enddo
  do i = 0, 7
     if(tmpnode(i)%nump .gt.0) goto 20
  enddo
20 continue
  iFirstDaughter = iSpareNode
  iNowDaughter = iSpareNode
  if(Division(worknode, tmpnode(i)%nump)) then
     worknode%daughter = iFirstDaughter
  else
     worknode%daughter = tmpnode(i)%daughter
  endif
  iLeftSibling = 0
  do i = 0, 7
     if(Division(worknode, tmpnode(i)%nump)) then
        node(-iNowDaughter)%daughter = tmpnode(i)%daughter
        node(-iNowDaughter)%nump = tmpnode(i)%nump
        if(iLeftSibling .gt. 0) then
           particle(iLeftSibling)%sibling = iNowDaughter
        else if (iLeftSibling .lt. 0) then
           node(-iLeftSibling)%sibling = iNowDaughter
        endif
        iLeftSibling = iNowDaughter
        iNowDaughter = iNowDaughter -1

     else if (tmpnode(i)%nump .gt.0) then
        ntmp = tmpnode(i)%daughter
        if(iLeftSibling .gt. 0) then
           particle(iLeftSibling)%sibling = ntmp
        else if (iLeftSibling .lt. 0) then
           node(-iLeftSibling)%sibling = ntmp
        endif
        do while(ntmp .ne.0)
           iLeftSibling = ntmp
           ntmp = particle(ntmp)%sibling
        enddo
     endif
  enddo
  if(iLeftSibling .gt.0) then
     particle(iLeftSibling)%sibling = worknode%sibling
  else if (iLeftSibling .lt.0) then
     node(-iLeftSibling)%sibling = worknode%sibling
  endif

  iSpareNode = iNowDaughter


  inextjobnode = iFirstDaughter
  do i = 0, 7
     if(Division(worknode, tmpnode(i)%nump)) then
        j = inextjobnode
        k = iSpareNode
        iSpareNode = Make_Tree(node, particle,np, j, k)
        inextjobnode = inextjobnode -1
     endif
  enddo
  iNextSpareNode = iSpareNode

end function Make_Tree


!@################################################
subroutine pthreadLinkedList(ipart,nworks,nthreads,nparticles,ptrhead,nextlink)
   integer, intent(in):: ipart, nworks
   integer:: eqlwrks, nthreads,icount,ithread,i,jpart
   integer, dimension(1:*), intent(in):: nextlink
   integer, dimension(0:nthreads-1), intent(out):: nparticles, ptrhead
   integer, dimension(0:nthreads-1):: beforenwrks
   integer:: runningwrks,remainwrks,allocworks,k
   eqlwrks = (nworks+nthreads-1)/ nthreads
   runningwrks = 0
   remainwrks = nworks
   jpart = ipart

!   k = 1000
!   do while( k.ge.0)
!      k = 1000
!   enddo

   do i = 0, nthreads-1
       if(remainwrks .ge. eqlwrks) then
          allocworks = eqlwrks
       else
          allocworks = remainwrks
       endif
       nparticles(i) = allocworks
       runningwrks = runningwrks + allocworks
       remainwrks = remainwrks - allocworks
   enddo
   beforenwrks(0) = 0
   do i = 1, nthreads-1
      beforenwrks(i) = nparticles(i-1) + beforenwrks(i-1)
   enddo
!  ptrhead(0) = jpart
!  ithread = 1;
   ptrhead = 0

   ithread = 0

   do i=1,nworks
       if(i .eq. beforenwrks(ithread)+1) then
          ptrhead(ithread) = jpart
          ithread = ithread + 1
          if(ithread .ge. nthreads) then
             return
          endif
       endif
       jpart=nextlink(jpart)
   enddo
   return
end subroutine pthreadLinkedList

