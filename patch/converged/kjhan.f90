!@################################################
! (C)  Copyright 2018-24 Dr. Juhan Kim's Software 


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
              if(pfof_open(node, particle,np,gndx, ptr, foflink, boxsize).eq. .true.) then
                 optr = ptr
                 ptr = node(-ptr)%daughter
              else
                 optr = ptr;
                 ptr = node(-ptr)%sibling
              endif
           endif
        else if(ptr .gt.0) then
           if(particle(ptr)%included .eq. .true.) then
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


integer function Do_FoF_Tree(numsink, psink, gsink, dx_min2, factG) result(gindx)
#ifdef _KJHAN_TEST
  use AA
#else
  use pm_commons
  use amr_commons
#endif
  use Tree_commons
  implicit none
  integer:: i,j,k, indx,numsink
  TYPE(Tree_Node),  dimension(1:numsink):: node
  TYPE(Tree_Particle),  dimension(1:numsink):: particle
  integer, dimension(1:numsink):: linked
  integer,dimension(1:numsink)::psink,gsink
  integer:: gnum,ii,icount
  real(dp), dimension(3):: boxsize
  real(dp):: foflink
  real(dp):: dx_min2
  real(dp):: factG, scale
  real(dp),dimension(1:3)::xbound,skip_loc
  integer, external:: FoF_link
  integer:: nx_loc


  xbound(1:3)=(/dble(nx),dble(ny),dble(nz)/)
  nx_loc=(icoarse_max-icoarse_min+1)
  scale=boxlen/dble(nx_loc)



  foflink = dsqrt(rmerge**2*dx_min2)
  boxsize(1) = scale*xbound(1)
  boxsize(2) = scale*xbound(2)
  boxsize(3) = scale*xbound(3)
  do i = 1, numsink
     indx = psink(i)
     particle(i)%x = xsink(indx,1)
     particle(i)%y = xsink(indx,2)
     particle(i)%z = xsink(indx,3)
     particle(i)%included = .false.
  enddo
  call Build_Tree(node, particle, numsink)

  gindx = 0
  icount = 0
  do i = 1, numsink
     if(particle(i)%included .eq. .false.) then
        ii = i
        gindx = gindx + 1
        gnum = FoF_link(ii, node, particle, numsink, foflink, linked, gindx, boxsize, factG)
        do j = 1, gnum
           psink(icount+j) = linked(j)
           gsink(linked(j)) = gindx
        enddo
        icount = icount + gnum
     endif
  enddo


end function Do_FoF_Tree

subroutine Build_Tree(node, particle, np)
  use Tree_commons
  implicit none
  integer:: np,i,j,k,iSpareNode
  TYPE(Tree_Node), dimension(1:np):: node
  TYPE(Tree_Particle), dimension(1:np):: particle
  integer, external:: Make_Tree
  integer:: iworknode, iSpareNode0


  do i = 1, np
     particle(i)%sibling = i + 1
     node(i)%sibling = 0
     node(i)%daughter = 0
  enddo
  particle(np).sibling = 0
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
           particle(iLeftSibling).sibling = iNowDaughter
        else if (iLeftSibling .lt. 0) then
           node(-iLeftSibling).sibling = iNowDaughter
        endif
        iLeftSibling = iNowDaughter
        iNowDaughter = iNowDaughter -1 

     else if (tmpnode(i)%nump .gt.0) then
        ntmp = tmpnode(i)%daughter
        if(iLeftSibling .gt. 0) then
           particle(iLeftSibling).sibling = ntmp
        else if (iLeftSibling .lt. 0) then
           node(-iLeftSibling).sibling = ntmp
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
     if(Division(worknode, tmpnode(i).nump)) then
        j = inextjobnode
        k = iSpareNode
        iSpareNode = Make_Tree(node, particle,np, j, k)
        inextjobnode = inextjobnode -1
     endif
  enddo
  iNextSpareNode = iSpareNode
  
end function Make_Tree


!@################################################
!@################################################
!@################################################
!@################################################
!@################################################
!@################################################
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
