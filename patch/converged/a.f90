subroutine merge_sink(ilevel)
  use pm_commons
  use amr_commons
  use cooling_module, ONLY: XH=>X, rhoc, mH, twopi
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer::ilevel,ilun
  !------------------------------------------------------------------------
  ! This routine merges sink usink the FOF algorithm.
  ! It keeps only the group centre of mass and remove other sinks.
  !------------------------------------------------------------------------
  integer::j,isink,ii,jj,kk,ind,idim,new_sink
  real(dp)::dx_loc,scale,dx_min,xx,yy,zz,rr,rmax2,rmax
  integer::igrid,jgrid,ipart,jpart,next_part,info
  integer::i,ig,ip,npart1,npart2,icpu,nx_loc
  integer::igrp,icomp,gndx,ifirst,ilast,indx
  integer,dimension(1:nvector),save::ind_grid,ind_part,ind_grid_part
  integer,dimension(:),allocatable::psink,gsink
  real(dp),dimension(1:3)::xbound,skip_loc
  real(dp)::dx,vol_min

  integer,dimension(:),allocatable::rank_old,idsink_old
  real(dp),dimension(:),allocatable::tsink_old
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  real(dp)::d0,mstar,nISM,nCOM,vr,cs,factG
  real(dp)::egrav,ekin,uxcom,uycom,uzcom,v2rel1,v2rel2,dx_min2
  real(dp)::xc,yc,zc,vxc,vyc,vzc,xx1,yy1,zz1,xx2,yy2,zz2,vvx1,vvy1,vvz1,vvx2,vvy2,vvz2,Lx1,Ly1,Lz1,Lx2,Ly2,Lz2
  real(dp)::q,M1,M2,a1,a2,a1a2,a1L,a2L,Lx,Ly,Lz,Lmod,a1mod,a2mod,mu,af,Lmodana
  real(dp)::ax1,ay1,az1,ax2,ay2,az2,x1,y1,z1,x2,y2,z2,vx1,vy1,vz1,vx2,vy2,vz2
  real(dp)::s4=-0.129d0,s5=-0.384d0,t0=-2.686d0,t2=-3.454d0,t3=+2.353d0,pi

  pi=twopi/2.0d0
  factG=1d0
  if(cosmo)factG=3d0/8d0/pi*omega_m*aexp

  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

  if(numbtot(1,ilevel)==0)return
  if(nsink==0)return
  if(verbose)write(*,111)ilevel

  ! Mesh spacing in that level
  dx_loc=0.5D0**ilevel
  xbound(1:3)=(/dble(nx),dble(ny),dble(nz)/)
  nx_loc=(icoarse_max-icoarse_min+1)
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  if(ndim>0)skip_loc(1)=dble(icoarse_min)
  if(ndim>1)skip_loc(2)=dble(jcoarse_min)
  if(ndim>2)skip_loc(3)=dble(kcoarse_min)
  scale=boxlen/dble(nx_loc)
  dx_min=scale*0.5D0**nlevelmax/aexp
  dx_min2=dx_min*dx_min
  rmax=4.0D0*dx_min
  rmax2=rmax*rmax

  allocate(psink(1:nsink),gsink(1:nsink))

  allocate(rank_old(1:nsink),idsink_old(1:nsink),tsink_old(1:nsink))

  !-------------------------------
  ! Merge sinks using FOF
  !-------------------------------
  do isink=1,nsink
     psink(isink)=isink
     gsink(isink)=0
  end do

  igrp=0
  icomp=1
  ifirst=2
  do while(icomp.le.nsink)
     gndx=psink(icomp)
     if(gsink(gndx)==0)then
    igrp=igrp+1
    gsink(gndx)=igrp
     endif
     ilast=nsink
     do while((ilast-ifirst+1)>0)
    indx=psink(ifirst)
    xx=xsink(indx,1)-xsink(gndx,1)
    if(xx>scale*xbound(1)/2.0)then
       xx=xx-scale*xbound(1)
    endif
    if(xx<-scale*xbound(1)/2.0)then
       xx=xx+scale*xbound(1)
    endif
    rr=xx**2
#if NDIM>1
    yy=xsink(indx,2)-xsink(gndx,2)
    if(yy>scale*xbound(2)/2.0)then
       yy=yy-scale*xbound(2)
    endif
    if(yy<-scale*xbound(2)/2.0)then
       yy=yy+scale*xbound(2)
    endif
    rr=yy**2+rr
#endif
#if NDIM>2
    zz=xsink(indx,3)-xsink(gndx,3)
    if(zz>scale*xbound(3)/2.0)then
       zz=zz-scale*xbound(3)
    endif
    if(zz<-scale*xbound(3)/2.0)then
       zz=zz+scale*xbound(3)
    endif
    rr=zz**2+rr
#endif
        if(rr.le.rmerge**2.*dx_min2)then
           if(vrel_merge)then
              egrav=msink(indx)*msink(gndx)/(rr+tiny(0.d0))*factG
              uxcom=(msink(indx)*vsink(indx,1)+msink(gndx)*vsink(gndx,1))/(msink(indx)+msink(gndx))
              uycom=(msink(indx)*vsink(indx,2)+msink(gndx)*vsink(gndx,2))/(msink(indx)+msink(gndx))
              uzcom=(msink(indx)*vsink(indx,3)+msink(gndx)*vsink(gndx,3))/(msink(indx)+msink(gndx))
              v2rel1=(vsink(indx,1)-uxcom)**2+(vsink(indx,2)-uycom)**2+(vsink(indx,3)-uzcom)**2
              v2rel2=(vsink(gndx,1)-uxcom)**2+(vsink(gndx,2)-uycom)**2+(vsink(gndx,3)-uzcom)**2
              ekin=0.5d0*(msink(indx)*v2rel1+msink(gndx)*v2rel2)
              if(ekin.lt.egrav)then
                 ifirst=ifirst+1
                 gsink(indx)=igrp
              else
                 psink(ifirst)=psink(ilast)
                 psink(ilast)=indx
                 ilast=ilast-1
              endif
           else
              ifirst=ifirst+1
              gsink(indx)=igrp
           endif
        else
           psink(ifirst)=psink(ilast)
           psink(ilast)=indx
           ilast=ilast-1
        endif
     end do
     icomp=icomp+1
  end do
  new_sink=igrp
  if(myid==1)then
     write(*,*)'Found ',new_sink,' groups'
     !do isink=1,nsink
     !     write(*,'(3(I4,1x),3(1PE10.3))')isink,psink(isink),gsink(isink),xsink(isink,1:ndim)
     !end do
  endif

  !----------------------------------------------------
  ! Compute group centre of mass and average velocity
  !----------------------------------------------------
  xsink_new=0d0; vsink_new=0d0; msink_new=0d0; dMsmbh_new=0d0; Esave_new=0d0; idsink_new=0
  oksink_all=0d0; oksink_new=0d0; tsink_new=0d0
  rank_old=0d0; idsink_old=0d0; tsink_old=0d0
  bhspin_new=0d0; spinmag_new=0d0
  do isink=1,nsink
     igrp=gsink(isink)

     !----------------------------------------------------
     ! This is done to keep track of the most massive sink
     ! after a merger with a companion
     !----------------------------------------------------
     if ( rank_old(igrp) .eq. 0)then
    rank_old(igrp)=isink
    idsink_old(igrp)=idsink(isink)
    tsink_old(igrp) =tsink(isink)
     endif
     if ( msink(isink) .gt. msink(rank_old(igrp)) )then
    rank_old(igrp)=isink
    idsink_new(igrp)=idsink(isink)
    idsink_old(igrp)=idsink(isink)
    tsink_new(igrp) =tsink(isink)
    tsink_old(igrp) =tsink(isink)
     else
    idsink_new(igrp)=idsink_old(igrp)
    tsink_new(igrp) =tsink_old(igrp)
     endif
     !----------------------------------------------------
     !----------------------------------------------------

     !idsink_new(igrp)=idsink(isink)
     if(oksink_new(igrp)==0d0)then
    oksink_all(isink)=igrp
    oksink_new(igrp)=isink
     endif
     !YDspin----------------------------------------------------
     if(msink_new(igrp).ge.msink(isink))then
        if(msink(isink).eq.0d0)then
           write(*,*)'Problem in merge_sink for spins, msink=0'
           stop
        endif
        M1=msink_new(igrp)
        M2=msink(isink)
        a1=ABS(spinmag_new(igrp))
        a2=ABS(spinmag(isink))
        ax1=bhspin_new(igrp,1)
        ay1=bhspin_new(igrp,2)
        az1=bhspin_new(igrp,3)
        ax2=bhspin(isink,1)
        ay2=bhspin(isink,2)
        az2=bhspin(isink,3)
        xx1=xsink_new(igrp,1)
        yy1=xsink_new(igrp,2)
        zz1=xsink_new(igrp,3)
        xx2=xsink(isink,1)
        yy2=xsink(isink,2)
        zz2=xsink(isink,3)
        vvx1=vsink_new(igrp,1)
        vvy1=vsink_new(igrp,2)
        vvz1=vsink_new(igrp,3)
        vvx2=vsink(isink,1)
        vvy2=vsink(isink,2)
        vvz2=vsink(isink,3)
     else
        if(msink_new(igrp).gt.0d0)then
           M1=msink(isink)
           M2=msink_new(igrp)
           a1=ABS(spinmag(isink))
           a2=ABS(spinmag_new(igrp))
           ax1=bhspin(isink,1)
           ay1=bhspin(isink,2)
           az1=bhspin(isink,3)
           ax2=bhspin_new(igrp,1)
           ay2=bhspin_new(igrp,2)
           az2=bhspin_new(igrp,3)
           xx1=xsink(isink,1)
           yy1=xsink(isink,2)
           zz1=xsink(isink,3)
           xx2=xsink_new(igrp,1)
           yy2=xsink_new(igrp,2)
           zz2=xsink_new(igrp,3)
           vvx1=vsink(isink,1)
           vvy1=vsink(isink,2)
           vvz1=vsink(isink,3)
           vvx2=vsink_new(igrp,1)
           vvy2=vsink_new(igrp,2)
           vvz2=vsink_new(igrp,3)
        else
           M2=0d0
        endif
     endif
     if(M2.ne.0d0.and.bhspinmerge)then
        q=M2/M1
        mu=q/(1d0+q)**2
        a1mod=SQRT(ax1**2+ay1**2+az1**2)
        a2mod=SQRT(ax2**2+ay2**2+az2**2)
        if(a1mod>0.and.a2mod>0)then
           a1a2=(ax1*ax2+ay1*ay2+az1*az2)/(a1mod*a2mod)
        else
           a1a2=0d0
        endif        
        ! C.O.M.
        ! Check for periodicity
        xx=xx2-xx1
        if(xx>scale*xbound(1)/2.0)then
           xx2=xx2-scale*xbound(1)
        endif
        if(xx<-scale*xbound(1)/2.0)then
           xx1=xx1+scale*xbound(1)
        endif        
        yy=yy2-yy1
        if(yy>scale*xbound(2)/2.0)then
           yy2=yy2-scale*xbound(2)
        endif
        if(yy<-scale*xbound(2)/2.0)then
           yy1=yy1+scale*xbound(2)
        endif
        zz=zz2-zz1
        if(zz>scale*xbound(3)/2.0)then
           zz2=zz2-scale*xbound(3)
        endif
        if(zz<-scale*xbound(3)/2.0)then
           zz1=zz1+scale*xbound(3)
        endif
        xc =(M1*xx1+M2*xx2)/(M1+M2)
        yc =(M1*yy1+M2*yy2)/(M1+M2)
        zc =(M1*zz1+M2*zz2)/(M1+M2)
        vxc=(M1*vvx1+M2*vvx2)/(M1+M2)
        vyc=(M1*vvy1+M2*vvy2)/(M1+M2)
        vzc=(M1*vvz1+M2*vvz2)/(M1+M2)
        x1 =xx1-xc
        y1 =yy1-yc
        z1 =zz1-zc
        x2 =xx2-xc
        y2 =yy2-yc
        z2 =zz2-zc        
        vx1=vvx1-vxc
        vy1=vvy1-vyc
        vz1=vvz1-vzc
        vx2=vvx2-vxc
        vy2=vvy2-vyc
        vz2=vvz2-vzc
        Lx1=M1*(y1*vz1-z1*vy1)
        Ly1=M1*(z1*vx1-x1*vz1)
        Lz1=M1*(x1*vy1-y1*vx1)
        Lx2=M2*(y2*vz2-z2*vy2)
        Ly2=M2*(z2*vx2-x2*vz2)
        Lz2=M2*(x2*vy2-y2*vx2)
        Lx=Lx1+Lx2
        Ly=Ly1+Ly2
        Lz=Lz1+Lz2
        Lmod =SQRT(Lx**2+Ly**2+Lz**2)
        if(a1mod>0.and.a2mod>0.and.Lmod>0.)then
           a1L=(ax1*Lx+ay1*Ly+az1*Lz)/(a1mod*Lmod)
           a2L=(ax2*Lx+ay2*Ly+az2*Lz)/(a2mod*Lmod)
        else
           a1L=0d0
           a2L=0d0
        endif
        Lmodana = s4 / ( 1d0 + q**2 )**2 * ( a1**2 + a2**2*q**4 + 2d0*a1*a2*q**2*a1a2 ) &
             & + ( s5*mu + t0 + 2d0 ) / ( 1d0 + q**2 ) * ( a1*a1L + a2*q**2*a2L ) &
             & + 2d0*SQRT(3d0) + t2*mu + t3*mu**2     
        af = 1d0 / ( 1d0 + q )**2 * SQRT( a1**2 + a2**2*q**4 + 2d0*a1*a2*q**2*a1a2 &
             & + 2d0 * ( a1*a1L + a2*q**2*a2L ) * Lmodana*q + Lmodana**2*q**2 )     
        bhspin_new(igrp,1)=1d0/(1d0+q)**2 * ( a1*ax1+a2*ax2*q**2+(Lx/Lmod)*Lmodana*q )
        bhspin_new(igrp,2)=1d0/(1d0+q)**2 * ( a1*ay1+a2*ay2*q**2+(Ly/Lmod)*Lmodana*q )
        bhspin_new(igrp,3)=1d0/(1d0+q)**2 * ( a1*az1+a2*az2*q**2+(Lz/Lmod)*Lmodana*q )
        spinmag_new(igrp)=SQRT(bhspin_new(igrp,1)**2 + bhspin_new(igrp,2)**2 + bhspin_new(igrp,3)**2 )
        if(spinmag_new(igrp).gt.+maxspin) spinmag_new(igrp)=+maxspin
        if(spinmag_new(igrp).lt.-maxspin) spinmag_new(igrp)=-maxspin     
        ax1=bhspin_new(igrp,1)
        ay1=bhspin_new(igrp,2)
        az1=bhspin_new(igrp,3)
        a1mod=SQRT(ax1**2+ay1**2+az1**2)
        bhspin_new(igrp,1)=ax1/a1mod
        bhspin_new(igrp,2)=ay1/a1mod
        bhspin_new(igrp,3)=az1/a1mod
     else
        if ( msink(isink) .gt. msink_new(igrp) )then        
           ! Case it's not a merger
           bhspin_new(igrp,1:ndim)=bhspin(isink,1:ndim)
           spinmag_new(igrp)=spinmag(isink)
        endif
     endif
     !YDspin----------------------------------------------------

     msink_new (igrp)=msink_new (igrp)+msink (isink)
     dMsmbh_new(igrp)=dMsmbh_new(igrp)+dMsmbh(isink)
     Esave_new (igrp)=Esave_new (igrp)+Esave (isink)
     xx=xsink(isink,1)-xsink(int(oksink_new(igrp)),1)
     if(xx>scale*xbound(1)/2.0)then
    xx=xx-scale*xbound(1)
     endif
     if(xx<-scale*xbound(1)/2.0)then
    xx=xx+scale*xbound(1)
     endif
     xsink_new(igrp,1)=xsink_new(igrp,1)+msink(isink)*xx
     vsink_new(igrp,1)=vsink_new(igrp,1)+msink(isink)*vsink(isink,1)
#if NDIM>1
     yy=xsink(isink,2)-xsink(int(oksink_new(igrp)),2)
     if(yy>scale*xbound(2)/2.0)then
    yy=yy-scale*xbound(2)
     endif
     if(yy<-scale*xbound(2)/2.0)then
    yy=yy+scale*xbound(2)
     endif
     xsink_new(igrp,2)=xsink_new(igrp,2)+msink(isink)*yy
     vsink_new(igrp,2)=vsink_new(igrp,2)+msink(isink)*vsink(isink,2)
#endif
#if NDIM>2
     zz=xsink(isink,3)-xsink(int(oksink_new(igrp)),3)
     if(zz>scale*xbound(3)/2.0)then
    zz=zz-scale*xbound(3)
     endif
     if(zz<-scale*xbound(3)/2.0)then
    zz=zz+scale*xbound(3)
     endif
     xsink_new(igrp,3)=xsink_new(igrp,3)+msink(isink)*zz
     vsink_new(igrp,3)=vsink_new(igrp,3)+msink(isink)*vsink(isink,3)
#endif
  end do
  do isink=1,new_sink
     xsink_new(isink,1)=xsink_new(isink,1)/msink_new(isink)+xsink(int(oksink_new(isink)),1)
     vsink_new(isink,1)=vsink_new(isink,1)/msink_new(isink)
#if NDIM>1
     xsink_new(isink,2)=xsink_new(isink,2)/msink_new(isink)+xsink(int(oksink_new(isink)),2)
     vsink_new(isink,2)=vsink_new(isink,2)/msink_new(isink)
#endif
#if NDIM>2
     xsink_new(isink,3)=xsink_new(isink,3)/msink_new(isink)+xsink(int(oksink_new(isink)),3)
     vsink_new(isink,3)=vsink_new(isink,3)/msink_new(isink)
#endif
  end do
  nsink=new_sink
  msink (1:nsink)=msink_new (1:nsink)
  dMsmbh(1:nsink)=dMsmbh_new(1:nsink)
  Esave (1:nsink)=Esave_new (1:nsink)
  idsink(1:nsink)=idsink_new(1:nsink)
  tsink (1:nsink)=tsink_new (1:nsink)
  xsink(1:nsink,1:ndim)=xsink_new(1:nsink,1:ndim)
  vsink(1:nsink,1:ndim)=vsink_new(1:nsink,1:ndim)
  bhspin(1:nsink,1:ndim)=bhspin_new(1:nsink,1:ndim)
  spinmag(1:nsink)=spinmag_new(1:nsink)
  ! Periodic boundary conditions
  do isink=1,nsink
     xx=xsink(isink,1)
     if(xx<-scale*skip_loc(1))then
    xx=xx+scale*(xbound(1)-skip_loc(1))
     endif
     if(xx>scale*(xbound(1)-skip_loc(1)))then
    xx=xx-scale*(xbound(1)-skip_loc(1))
     endif
     xsink(isink,1)=xx
#if NDIM>1
     yy=xsink(isink,2)
     if(yy<-scale*skip_loc(2))then
    yy=yy+scale*(xbound(2)-skip_loc(2))
     endif
     if(yy>scale*(xbound(2)-skip_loc(2)))then
    yy=yy-scale*(xbound(2)-skip_loc(2))
     endif
     xsink(isink,2)=yy
#endif
#if NDIM>2
     zz=xsink(isink,3)
     if(zz<-scale*skip_loc(3))then
    zz=zz+scale*(xbound(3)-skip_loc(3))
     endif
     if(zz>scale*(xbound(3)-skip_loc(3)))then
    zz=zz-scale*(xbound(3)-skip_loc(3))
     endif
     xsink(isink,3)=zz
#endif
  enddo

  deallocate(psink,gsink)

  deallocate(rank_old,idsink_old,tsink_old)

  !-----------------------------------------------------
  ! Remove sink particles that are part of a FOF group.
  !-----------------------------------------------------
  ! Loop over cpus
  do icpu=1,ncpu
     igrid=headl(icpu,ilevel)
     ig=0
     ip=0
     ! Loop over grids
     do jgrid=1,numbl(icpu,ilevel)
    npart1=numbp(igrid)  ! Number of particles in the grid
    npart2=0

    ! Count sink particles
    if(npart1>0)then
       ipart=headp(igrid)
       ! Loop over particles
       do jpart=1,npart1
          ! Save next particle   <--- Very important !!!
          next_part=nextp(ipart)
          if(idp(ipart).lt.0 .and. tp(ipart).eq.0.d0)then
         npart2=npart2+1
          endif
          ipart=next_part  ! Go to next particle
       end do
    endif

    ! Gather sink particles
    if(npart2>0)then
       ig=ig+1
       ind_grid(ig)=igrid
       ipart=headp(igrid)
       ! Loop over particles
       do jpart=1,npart1
          ! Save next particle   <--- Very important !!!
          next_part=nextp(ipart)
          ! Select only sink particles
          if(idp(ipart).lt.0 .and. tp(ipart).eq.0.d0)then
         if(ig==0)then
            ig=1
            ind_grid(ig)=igrid
         end if
         ip=ip+1
         ind_part(ip)=ipart
         ind_grid_part(ip)=ig
          endif
          if(ip==nvector)then
         call kill_sink(ind_grid,ind_part,ind_grid_part,ig,ip,ilevel)
         ip=0
         ig=0
          end if
          ipart=next_part  ! Go to next particle
       end do
       ! End loop over particles
    end if

    igrid=next(igrid)   ! Go to next grid
     end do

     ! End loop over grids
     if(ip>0)call kill_sink(ind_grid,ind_part,ind_grid_part,ig,ip,ilevel)
  end do
  ! End loop over cpus

111 format('   Entering merge_sink for level ',I2)

end subroutine merge_sink
