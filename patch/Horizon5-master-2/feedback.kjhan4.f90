!################################################################
!################################################################
!################################################################
!################################################################
subroutine thermal_feedback(ilevel)
  use pm_commons
  use amr_commons
  implicit none
  integer::ilevel
  !------------------------------------------------------------------------
  ! This routine computes the thermal energy, the kinetic energy and 
  ! the metal mass dumped in the gas by stars (SNII, SNIa, winds).
  ! This routine is called every fine time step.
  !------------------------------------------------------------------------
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v, &
       t1,xx,xx0,yy,ESN1,ESN2

  real(dp)::t0,scale,dx_min,vsn,rdebris,ethermal
  real::wxa,wxb,wya,wyb,wxa0,wxb0,rmin,rmax,dx1,yieldval, &
       yieldval1,val1,val2,valsn1,zmin,zmax,dy1,valsn2


  integer::igrid,jgrid,ipart,jpart,next_part,iastar
  integer::i,ig,ip,npart1,npart2,icpu,nx_loc
  integer::info,iSN,ivar,ielt,n11,n22,ihx, &
       ihy,ii,ihx0,indt1,indt2,som
  logical ::ok_free
! real(dp),dimension(1:3)::skip_loc
! integer,dimension(1:nvector)::ind_grid,ind_part,ind_grid_part
! real(dp),dimension(1:nvector)::mej,zej,nbsn1,nbsn2
! real(dp),dimension(1:nvector,1:nelt)::ceej

  integer mythread, nthreads,nwork,icount,jcount, npart3, subnump
  common /thermal_fb/ mythread
!$omp threadprivate(/thermal_fb/)
  integer, dimension(:), allocatable:: nparticles, ptrhead

! common /thermal_feedback_units/ scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2

  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel

! ESN1=eps_sn1*(1d51/2.d33)/scale_v**2
! ESN2=eps_sn2*(1d51/2.d33)/scale_v**2
! 

  n11=yieldtab%na
  n22=yieldtab%nz



!$omp parallel
  mythread = omp_get_thread_num()
  if(mythread.eq.0) nthreads = omp_get_num_threads()
!$omp end parallel
  allocate(ptrhead(0:nthreads-1), nparticles(0:nthreads-1))


#if NDIM==3
  do icpu=1,ncpu
  if(numbl(icpu,ilevel) .le.0) goto 12
  call pthreadLinkedList(headl(icpu,ilevel),numbl(icpu,ilevel),nthreads,nparticles,ptrhead,next)
!$omp parallel private(subnump,igrid) 
  subnump = nparticles(mythread) 
  igrid = ptrhead(mythread) 
  call sub_thermal_feedback(ilevel,icpu, igrid,subnump,n11,n22)
!$omp end parallel
12 end do
#endif
  deallocate(ptrhead, nparticles)

111 format('   Entering thermal_feedback for level ',I2)
end subroutine thermal_feedback


subroutine sub_thermal_feedback(ilevel,icpu, kgrid,subnump,n11,n22)
  use pm_commons
  use amr_commons
  implicit none
  integer,intent(in)::ilevel,icpu,kgrid,subnump,n11,n22
  !------------------------------------------------------------------------
  ! This routine computes the thermal energy, the kinetic energy and 
  ! the metal mass dumped in the gas by stars (SNII, SNIa, winds).
  ! This routine is called every fine time step.
  !------------------------------------------------------------------------
  integer::info,iSN,ivar,ielt,ihx, &
       ihy,ii,ihx0,indt1,indt2,som
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v, &
       t1,xx,xx0,yy,ESN1,ESN2
  real(dp)::t0,scale,dx_min,vsn,rdebris,ethermal
  real::wxa,wxb,wya,wyb,wxa0,wxb0,rmin,rmax,dx1,yieldval, &
       yieldval1,val1,val2,valsn1,zmin,zmax,dy1,valsn2
  real,dimension(1:nelt)::val3
  real,dimension(1:n11)::astarproper
  real,dimension(1:n11,1:n22)::Xeject
  integer::igrid,jgrid,ipart,jpart,next_part,iastar
  integer::i,ig,ip,npart1,npart2,nx_loc
  logical ::ok_free
  real(dp),dimension(1:3)::skip_loc
  integer,dimension(1:nvector)::ind_grid,ind_part,ind_grid_part
  real(dp),dimension(1:nvector)::mej,zej,nbsn1,nbsn2
  real(dp),dimension(1:nvector,1:nelt)::ceej
! common /thermal_feedback_units/ scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2


  igrid = kgrid

  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

  ESN1=eps_sn1*(1d51/2.d33)/scale_v**2
  ESN2=eps_sn2*(1d51/2.d33)/scale_v**2
  
!for test 3  
  do iastar=1,n11
     astarproper(iastar) = yieldtab%astar(iastar)*aexp**2*(1d9*365.*24.*3600./scale_t) !bring it from supercomoving to proper time in Gyr by multiplying by aexp**2
  end do

  rmin=log10(astarproper(1))
  rmax=log10(astarproper(n11))

  dx1=real(n11-1)/(rmax-rmin)

  zmin=log10(yieldtab%zstar(2))
  zmax=log10(yieldtab%zstar(n22))
  dy1=real(n22-2)/(zmax-zmin)

  ! Gather star particles eligible for a SN event 

#if NDIM==3
  ! Loop over cpus
  ig=0
  ip=0
  ! Loop over grids
  do jgrid=1,subnump
     npart1=numbp(igrid)  ! Number of particles in the grid
     npart2=0
        
     ! Count star particles
     if(npart1>0)then
        ipart=headp(igrid)
        ! Loop over particles
        do jpart=1,npart1
           ! Save next particle   <--- Very important !!!
           next_part=nextp(ipart)

           xx=texp-tpp(ipart)
           yy=zp(ipart)
           ihx=1+int(dx1*(log10(xx)-rmin))
           ihy=2+int(dy1*(log10(yy)-zmin))
           if (yy<yieldtab%zstar(2)) ihy=1
              
           if( tp(ipart)/=0. .and. xx>0. .and. ihx>0 )then
              wxa=(xx-astarproper(ihx))  /(astarproper(ihx+1)-astarproper(ihx))
              wxb = 1.-wxa

              wya=(yy-yieldtab%zstar(ihy))  /(yieldtab%zstar(ihy+1)-yieldtab%zstar(ihy))
              wyb=1.-wya

              xx0=indtab(ipart)
              if (xx0>0.d0) then
                 ihx0=1+int(dx1*(log10(xx0)-rmin))
                 wxa0=(xx0-astarproper(ihx0))    /(astarproper(ihx0+1)-astarproper(ihx0))
                 wxb0=1.-wxa0
              else
                 ihx0=1
                 wxa0=0.
                 wxb0=0.
              endif

              ok_free=.true.

              !Total ejected mass                                                       
              Xeject(:n11,:n22)=yieldtab%Meject(:n11,:n22)
              yieldval=Xeject(ihx,ihy)*wxb*wyb   + Xeject(ihx,ihy+1)*wxb*wya   + &
                   Xeject(ihx+1,ihy)*wxa*wyb + Xeject(ihx+1,ihy+1)*wxa*wya
              yieldval1=Xeject(ihx0,ihy)*wxb0*wyb   + Xeject(ihx0,ihy+1)*wxb0*wya   + &
                   Xeject(ihx0+1,ihy)*wxa0*wyb + Xeject(ihx0+1,ihy+1)*wxa0*wya
              if((yieldval-yieldval1)/(yieldval1+1.d-99)<tol)ok_free=.false.

              !Amount of ejected metals                                                 
              Xeject(:n11,:n22)=yieldtab%Zeject(:n11,:n22)
              yieldval=Xeject(ihx,ihy)*wxb*wyb   + Xeject(ihx,ihy+1)*wxb*wya   + &
                   Xeject(ihx+1,ihy)*wxa*wyb + Xeject(ihx+1,ihy+1)*wxa*wya
              yieldval1=Xeject(ihx0,ihy)*wxb0*wyb   + Xeject(ihx0,ihy+1)*wxb0*wya   + &
                   Xeject(ihx0+1,ihy)*wxa0*wyb + Xeject(ihx0+1,ihy+1)*wxa0*wya
              if((yieldval-yieldval1)/(yieldval1+1.d-99)<tol)ok_free=.false.

              !Nb of SNII                                                               
              valsn2=0.
              Xeject(:n11,:n22)=yieldtab%NSN2(:n11,:n22)
              yieldval=Xeject(ihx,ihy)*wxb*wyb   + Xeject(ihx,ihy+1)*wxb*wya   + &
                   Xeject(ihx+1,ihy)*wxa*wyb + Xeject(ihx+1,ihy+1)*wxa*wya
              yieldval1=Xeject(ihx0,ihy)*wxb0*wyb   + Xeject(ihx0,ihy+1)*wxb0*wya   + &
                   Xeject(ihx0+1,ihy)*wxa0*wyb + Xeject(ihx0+1,ihy+1)*wxa0*wya
              if((yieldval-yieldval1)/(yieldval1+1.d-99)>=tol)valsn2=yieldval-yieldval1
              !local feedback for AGB/SNIa                                              
              !local feedback for SNII only if f_w=0.                                   
              !otherwise, for f_w>0., spherical-feedback scheme (-->routine kinetic_feedback)                                             
              if( f_w>0. .and. valsn2>0. )ok_free=.false.

              if (ok_free) npart2=npart2+1


           endif
           ipart=next_part  ! Go to next particle
        end do
     endif
        
        ! Gather star particles
        if(npart2>0)then        
           ig=ig+1
           ind_grid(ig)=igrid
           ipart=headp(igrid)
           ! Loop over particles
           do jpart=1,npart1
              ! Save next particle   <--- Very important !!!
              next_part=nextp(ipart)
              ! Select only star particles
              xx=texp-tpp(ipart)
              yy=zp(ipart)
              ihx=1+int(dx1*(log10(xx)-rmin))
              ihy=2+int(dy1*(log10(yy)-zmin))
              if (yy<yieldtab%zstar(2)) ihy=1

              if( tp(ipart)/=0. .and. xx>0. .and. ihx>0 )then

                 wxa=(xx-astarproper(ihx))  /(astarproper(ihx+1)-astarproper(ihx))
                 wxb=1.-wxa

                 wya=(yy-yieldtab%zstar(ihy))  /(yieldtab%zstar(ihy+1)-yieldtab%zstar(ihy))
                 wyb=1.-wya

                 xx0=indtab(ipart)
                 if (xx0>0.d0) then
                    ihx0=1+int(dx1*(log10(xx0)-rmin))
                    wxa0=(xx0-astarproper(ihx0))    /(astarproper(ihx0+1)-astarproper(ihx0))
                    wxb0=1.-wxa0
                 else
                    ihx0=1
                    wxa0=0.
                    wxb0=0.
                 endif

                 ok_free=.true.

                 !Total ejected mass                                                                 
                 Xeject(:n11,:n22)=yieldtab%Meject(:n11,:n22)
                 yieldval=Xeject(ihx,ihy)*wxb*wyb   + Xeject(ihx,ihy+1)*wxb*wya   + &
                      Xeject(ihx+1,ihy)*wxa*wyb + Xeject(ihx+1,ihy+1)*wxa*wya
                 yieldval1=Xeject(ihx0,ihy)*wxb0*wyb   + Xeject(ihx0,ihy+1)*wxb0*wya   + &
                      Xeject(ihx0+1,ihy)*wxa0*wyb + Xeject(ihx0+1,ihy+1)*wxa0*wya
                 if((yieldval-yieldval1)/(yieldval1+1.d-99)<tol)ok_free=.false.
                 if((yieldval-yieldval1)/(yieldval1+1.d-99)>=tol)val1=yieldval-yieldval1

                 !Amount of ejected metals                                                           
                 Xeject(:n11,:n22)=yieldtab%Zeject(:n11,:n22)
                 yieldval=Xeject(ihx,ihy)*wxb*wyb   + Xeject(ihx,ihy+1)*wxb*wya   + &
                      Xeject(ihx+1,ihy)*wxa*wyb + Xeject(ihx+1,ihy+1)*wxa*wya
                 yieldval1=Xeject(ihx0,ihy)*wxb0*wyb   + Xeject(ihx0,ihy+1)*wxb0*wya   + &
                      Xeject(ihx0+1,ihy)*wxa0*wyb + Xeject(ihx0+1,ihy+1)*wxa0*wya
                 if((yieldval-yieldval1)/(yieldval1+1.d-99)<tol)ok_free=.false.
                 if((yieldval-yieldval1)/(yieldval1+1.d-99)>=tol)val2=yieldval-yieldval1

                 !Nb of SNII                                                                         
                 valsn2=0.
                 Xeject(:n11,:n22)=yieldtab%NSN2(:n11,:n22)
                 yieldval=Xeject(ihx,ihy)*wxb*wyb   + Xeject(ihx,ihy+1)*wxb*wya   + &
                      Xeject(ihx+1,ihy)*wxa*wyb + Xeject(ihx+1,ihy+1)*wxa*wya
                 yieldval1=Xeject(ihx0,ihy)*wxb0*wyb   + Xeject(ihx0,ihy+1)*wxb0*wya   + &
                      Xeject(ihx0+1,ihy)*wxa0*wyb + Xeject(ihx0+1,ihy+1)*wxa0*wya
                 if((yieldval-yieldval1)/(yieldval1+1.d-99)>=tol)valsn2=yieldval-yieldval1
                 if( f_w>0. .and. valsn2>0. )ok_free=.false.

                 if (ok_free) then

                    if(ig==0)then
                       ig=1
                       ind_grid(ig)=igrid
                    end if

                    !Nb of SNIa                                                                      
                    valsn1=0.
                    Xeject(:n11,:n22)=yieldtab%NSN1(:n11,:n22)
                    if (Xeject(ihx,ihy)>0.d0) then
                       yieldval=Xeject(ihx,ihy)*wxb*wyb   + Xeject(ihx,ihy+1)*wxb*wya   + &
                            Xeject(ihx+1,ihy)*wxa*wyb + Xeject(ihx+1,ihy+1)*wxa*wya
                       yieldval1=Xeject(ihx0,ihy)*wxb0*wyb   + Xeject(ihx0,ihy+1)*wxb0*wya   + &
                            Xeject(ihx0+1,ihy)*wxa0*wyb + Xeject(ihx0+1,ihy+1)*wxa0*wya
                       if((yieldval-yieldval1)/(yieldval1+1.d-99)>=tol)valsn1=yieldval-yieldval1
                    endif


                    !Amount of ejected chemical elements                                            
                    val3=0.
                    do ielt=1,nelt
                       Xeject(:n11,:n22)=yieldtab%Eeject(:n11,:n22,ielt)
                       yieldval=Xeject(ihx,ihy)*wxb*wyb   + Xeject(ihx,ihy+1)*wxb*wya   + &
                            Xeject(ihx+1,ihy)*wxa*wyb + Xeject(ihx+1,ihy+1)*wxa*wya
                       yieldval1=Xeject(ihx0,ihy)*wxb0*wyb   + Xeject(ihx0,ihy+1)*wxb0*wya   + &
                            Xeject(ihx0+1,ihy)*wxa0*wyb + Xeject(ihx0+1,ihy+1)*wxa0*wya
                       if((yieldval-yieldval1)/(yieldval1+1.d-99)>=tol)val3(ielt)=yieldval-yieldval1
                    enddo

                    ip=ip+1
                    ind_part(ip)=ipart
                    ind_grid_part(ip)=ig

                    zej(ip)=val2/val1
                    do ielt=1,nelt
                       ceej(ip,ielt)=val3(ielt)/val1
                    enddo
                    mej(ip)=val1*mp0(ipart)
                    nbsn1(ip)=valsn1*mp0(ipart)
                    nbsn2(ip)=valsn2*mp0(ipart)
                    mp(ipart)=mp(ipart)-mej(ip)
                    indtab(ipart)=xx

                 endif
              endif

              if(ip==nvector)then
                 call feedbk(ind_grid,ind_part,ind_grid_part,mej,nbsn1,nbsn2,zej,ceej,ESN1,ESN2,ig,ip,ilevel)
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
     if(ip>0)call feedbk(ind_grid,ind_part,ind_grid_part,mej,nbsn1,nbsn2,zej,ceej,ESN1,ESN2,ig,ip,ilevel)


#endif
!111 format('   Entering thermal_feedback for level ',I2)

end subroutine sub_thermal_feedback
!################################################################
!################################################################
!################################################################
!################################################################
subroutine feedbk(ind_grid,ind_part,ind_grid_part,mej,nbsn1,nbsn2,zej,ceej,ESN1,ESN2,ng,np,ilevel)
  use amr_commons
  use pm_commons
  use hydro_commons
  use random
  implicit none
  integer::ng,np,ilevel
  integer,dimension(1:nvector)::ind_grid
  integer,dimension(1:nvector)::ind_grid_part,ind_part
  real(dp),dimension(1:nvector)::mej,nbsn1,nbsn2,zej
  real(dp),dimension(1:nvector,1:nelt)::ceej
  !-----------------------------------------------------------------------
  ! This routine is called by subroutine feedback. Each stellar particle
  ! dumps mass, momentum and energy in the nearest grid cell using array
  ! unew.
  !-----------------------------------------------------------------------
  integer::i,j,idim,nx_loc,ielt
  real(kind=8)::RandNum
  real(dp)::SN_BOOST,dx_min
  real(dp)::ESN1,ESN2,vol_loc
  real(dp)::ERAD,RAD_BOOST,tauIR,eta_sig
  real(dp)::sigma_d,delta_x,tau_factor,rad_factor
  real(dp)::dx,dx_loc,scale,current_time
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  logical::error
  ! Grid based arrays
  real(dp),dimension(1:nvector,1:ndim)::x0
  integer ,dimension(1:nvector)::ind_cell
  integer ,dimension(1:nvector,1:threetondim)::nbors_father_cells
  integer ,dimension(1:nvector,1:twotondim)::nbors_father_grids
  ! Particle based arrays
  integer,dimension(1:nvector)::igrid_son,ind_son
  integer,dimension(1:nvector)::list1
  logical,dimension(1:nvector)::ok
  real(dp),dimension(1:nvector)::mloss,ekinetic
  real(dp),dimension(1:nvector,1:ndim):: x
  integer ,dimension(1:nvector,1:ndim)::id,igd,icd
  integer ,dimension(1:nvector)::igrid,icell,indp,kg
  real(dp),dimension(1:3)::skip_loc


  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

  ! Mesh spacing in that level
  dx=0.5D0**ilevel
  nx_loc=(icoarse_max-icoarse_min+1)
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  if(ndim>0)skip_loc(1)=dble(icoarse_min)
  if(ndim>1)skip_loc(2)=dble(jcoarse_min)
  if(ndim>2)skip_loc(3)=dble(kcoarse_min)
  scale=boxlen/dble(nx_loc)
  dx_loc=dx*scale
  vol_loc=dx_loc**ndim

#if NDIM==3
  ! Lower left corner of 3x3x3 grid-cube
  do idim=1,ndim
     do i=1,ng
        x0(i,idim)=xg(ind_grid(i),idim)-3.0D0*dx
     end do
  end do

  ! Gather 27 neighboring father cells (should be present anytime !)
  do i=1,ng
     ind_cell(i)=father(ind_grid(i))
  end do
  call get3cubefather(ind_cell,nbors_father_cells,nbors_father_grids,ng,ilevel)

  ! Rescale position at level ilevel
  do idim=1,ndim
     do j=1,np
        x(j,idim)=xp(ind_part(j),idim)/scale+skip_loc(idim)
     end do
  end do
  do idim=1,ndim
     do j=1,np
        x(j,idim)=x(j,idim)-x0(ind_grid_part(j),idim)
     end do
  end do
  do idim=1,ndim
     do j=1,np
        x(j,idim)=x(j,idim)/dx
     end do
  end do

  ! Removed since this is done right after anyway (in move_particles)
  !! Check for illegal moves
  !error=.false.
  !do idim=1,ndim
  !   do j=1,np
  !      if(x(j,idim)<=0.5D0.or.x(j,idim)>=5.5D0)error=.true.
  !   end do
  !end do
  !if(error)then
  !   write(*,*)'problem in sn2'
  !   write(*,*)ilevel,ng,np
  !end if

  ! NGP at level ilevel
  do idim=1,ndim
     do j=1,np
        id(j,idim)=x(j,idim)
     end do
  end do

   ! Compute parent grids
  do idim=1,ndim
     do j=1,np
        igd(j,idim)=id(j,idim)/2
     end do
  end do
  do j=1,np
     kg(j)=1+igd(j,1)+3*igd(j,2)+9*igd(j,3)
  end do
  do j=1,np
     igrid(j)=son(nbors_father_cells(ind_grid_part(j),kg(j)))
  end do

  ! Check if particles are entirely in level ilevel
  ok(1:np)=.true.
  do j=1,np
     ok(j)=ok(j).and.igrid(j)>0
  end do

  ! Compute parent cell position
  do idim=1,ndim
     do j=1,np
        if(ok(j))then
           icd(j,idim)=id(j,idim)-2*igd(j,idim)
        end if
     end do
  end do
  do j=1,np
     if(ok(j))then
        icell(j)=1+icd(j,1)+2*icd(j,2)+4*icd(j,3)
     end if
  end do

  ! Compute parent cell adresses
  do j=1,np
     if(ok(j))then
        indp(j)=ncoarse+(icell(j)-1)*ngridmax+igrid(j)
!     else THIS SECTION DOES EXIST IN THE LAST VERSION BUT WAS REMOVED IN THE CHEMO PATCH
!        indp(j) = nbors_father_cells(ind_grid_part(j),kg(j))
!        vol_loc(j)=vol_loc(j)*2**ndim ! ilevel-1 cell volume
     end if
  end do

  ! Reset ejected mass, metallicity, thermal energy
  do j=1,np
     mloss(j)=mej(j)/vol_loc
  end do

  ! Update hydro variables due to feedback

  ! For IR radiation trapping,
  ! we use a fixed length to estimate the column density of gas
!  delta_x=200.*3d18
!  if(metal)then
!     tau_factor=kappa_IR*delta_x*scale_d/0.02
!  else
!     tau_factor=kappa_IR*delta_x*scale_d*z_ave
!  endif
!  rad_factor=ERAD/ESN
!  do j=1,np

     ! Infrared photon trapping boost
!     if(metal)then
!        tauIR=tau_factor*max(uold(indp(j),imetal),smallr)
!     else
!        tauIR=tau_factor*max(uold(indp(j),1),smallr)
!     endif
!     if(uold(indp(j),1)*scale_nH > 10.)then
!        RAD_BOOST=rad_factor*(1d0-exp(-tauIR))
!     else
!        RAD_BOOST=0.0
!     endif


  ! Update hydro variables due to feedback                                                                                                
  do j=1,np
     if(ok(j))then
        ! Specific kinetic energy of the star                                                                                             
        ekinetic(j)=0.5*(vp(ind_part(j),1)**2 &
             &          +vp(ind_part(j),2)**2 &
             &          +vp(ind_part(j),3)**2)
        ! Update hydro variable in NGP cell                                                                                               
        unew(indp(j),1)=unew(indp(j),1)+mloss(j)
        unew(indp(j),2)=unew(indp(j),2)+mloss(j)*vp(ind_part(j),1)
        unew(indp(j),3)=unew(indp(j),3)+mloss(j)*vp(ind_part(j),2)
        unew(indp(j),4)=unew(indp(j),4)+mloss(j)*vp(ind_part(j),3)
        unew(indp(j),5)=unew(indp(j),5)+mloss(j)*ekinetic(j)+nbsn1(j)*ESN1/vol_loc+nbsn2(j)*ESN2/vol_loc
     endif
  end do
  ! Add metals
  if(metal)then
     do j=1,np
        if(ok(j))then
           unew(indp(j),imetal)=unew(indp(j),imetal)+mloss(j)*zej(j)
        endif
     end do
     do j=1,np
        if(ok(j))then
           do ielt=1,nelt
              unew(indp(j),ichem+ielt-1)=unew(indp(j),ichem+ielt-1)+mloss(j)*ceej(j,ielt)
           enddo
        endif
     enddo
  endif

  ! Add delayed cooling switch variable                                                                                                   
  if(delayed_cooling)then
     do j=1,np
        if(ok(j))then
           unew(indp(j),idelay)=unew(indp(j),idelay)+mloss(j)
        endif
     end do
  endif
#endif
  
end subroutine feedbk
!################################################################
!################################################################
!################################################################
!################################################################
subroutine sub2_kinetic_feedback(icpu,kgrid,subnump,n11,n22,iSN,nSN_tot,xSN_tot,vSN_tot,mSN_tot,ZSN_tot,NbSN_tot,ceSN_tot)
  use amr_commons
  use pm_commons
  use hydro_commons
  use cooling_module, ONLY: XH=>X, rhoc, mH 
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
  integer::nSN_tot_all,subnump,kgrid
  integer,dimension(1:ncpu)::nSN_icpu_all
  real(dp),dimension(:),allocatable::mSN_all,ZSN_all,NbSN_all
  real(dp),dimension(:,:),allocatable::xSN_all,vSN_all,ceSN_all
#endif
  !----------------------------------------------------------------------
  ! This subroutine compute the kinetic feedback due to SNII and
  ! imolement this using exploding GMC particles. 
  ! This routine is called only at coarse time step.
  ! Yohan Dubois
  !----------------------------------------------------------------------
  ! local constants
  integer::icpu,igrid,jgrid,npart1,npart2,ipart,jpart,next_part
  integer::nSN,nSN_loc,nSN_tot,info,iSN,ilevel,ivar,ielt,iSN2,iSN1,n11,n22,ihx, &
       ihy,ii,ihx0,i,nSN1,nSN2,nSN3,iSN3,indt1,indt2,nSN_tot2,npart22
  integer,dimension(1:ncpu)::nSN_icpu
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v,t0, &
       dxx,dyy,t1,xx,xx0,yy,ESN,t2,mass_load
  real::wxa,wxb,wya,wyb,wxa0,wxb0,rmin,rmax,dx1,yieldval,yieldval1, &
       val1,val2,valsn,zmin,zmax,dy1
  real,dimension(:),allocatable::val3
  real(dp)::current_time
  real(dp)::scale,dx_min,vol_min,nISM,nCOM,d0,mstar
  integer::nx_loc
!  integer,dimension(:),allocatable::ind_part,ind_grid
  logical ::ok_free
  integer ,dimension(:),allocatable::indSN,iSN_myid
  real(dp),dimension(:),allocatable::mSN,ZSN,vol_gas,ekBlast,NbSN,mloadSN,ZloadSN
  real(dp),dimension(:,:),allocatable::xSN,vSN,dq,ceSN,vloadSN,celoadSN
  real,dimension(:,:),allocatable::Xeject
  integer ,dimension(:),allocatable::indSN_tot,itemp
  real(dp),dimension(1:nSN_tot)::mSN_tot,ZSN_tot,NbSN_tot
  real(dp),dimension(1:nSN_tot,1:3)::xSN_tot,vSN_tot
  real(dp),dimension(1:nSN_tot,1:nelt)::ceSN_tot
  integer::isort,iastar
  real,dimension(:),allocatable::astarproper

  if(.not. hydro)return
  if(ndim.ne.3)return

 ! if(verbose .and. myid .eq. 1)write(*,*)'Entering make_sn'
  
  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

  igrid = kgrid

  ESN=eps_sn2*(1d51/2.d33)/scale_v**2

  n11=yieldtab%na
  n22=yieldtab%nz
  allocate(Xeject(1:n11,1:n22),val3(1:nelt),astarproper(1:n11))


!for test 3  
  do iastar=1,n11
     astarproper(iastar) = yieldtab%astar(iastar)*aexp**2*(1d9*365.*24.*3600./scale_t) !bring it from supercomoving to proper time in Gyr by multiplying by aexp**2
  end do

  rmin=log10(astarproper(1))
  rmax=log10(astarproper(n11))

  dx1=real(n11-1)/(rmax-rmin)

  zmin=log10(yieldtab%zstar(2))
  zmax=log10(yieldtab%zstar(n22))
  dy1=real(n22-2)/(zmax-zmin)

  !------------------------------------------------------                          
  ! Gather star particles in a SNII phase                                          
  !------------------------------------------------------
     do jgrid=1,subnump
        npart1=numbp(igrid)  ! Number of particles in the grid
        ! Count old enough star particles that have not exploded
        if(npart1>0)then
           ipart=headp(igrid)
           ! Loop over particles
           do jpart=1,npart1
              ! Save next particle   <--- Very important !!!
              next_part=nextp(ipart)

              xx=texp-tpp(ipart)
              yy=zp(ipart)
              ihx=1+int(dx1*(log10(xx)-rmin))
              ihy=2+int(dy1*(log10(yy)-zmin))
              if (yy<yieldtab%zstar(2)) ihy=1

              if( tp(ipart)/=0. .and. xx>0. .and. ihx>0 )then
                 wxa=(xx-astarproper(ihx))  /(astarproper(ihx+1)-astarproper(ihx))
                 wxb=1.-wxa

                 wya=(yy-yieldtab%zstar(ihy))  /(yieldtab%zstar(ihy+1)-yieldtab%zstar(ihy))
                 wyb = 1.-wya

                 xx0=indtab(ipart)
                 if (xx0>0.d0) then
                    ihx0=1+int(dx1*(log10(xx0)-rmin))
                    wxa0=(xx0-astarproper(ihx0))    /(astarproper(ihx0+1)-astarproper(ihx0))
                    wxb0 = 1.-wxa0
                 else
                    ihx0=1
                    wxa0=0.
                    wxb0=0.
                 endif

                 ok_free=.true.

                 !Total ejected mass                                                                                                                         
                 Xeject(:n11,:n22)=yieldtab%Meject(:n11,:n22)
                 yieldval=Xeject(ihx,ihy)*wxb*wyb   + Xeject(ihx,ihy+1)*wxb*wya   + &
                      Xeject(ihx+1,ihy)*wxa*wyb + Xeject(ihx+1,ihy+1)*wxa*wya
                 yieldval1=Xeject(ihx0,ihy)*wxb0*wyb   + Xeject(ihx0,ihy+1)*wxb0*wya   + &
                      Xeject(ihx0+1,ihy)*wxa0*wyb + Xeject(ihx0+1,ihy+1)*wxa0*wya
                 if((yieldval-yieldval1)/(yieldval1+1.d-99)<tol)ok_free=.false.
                 if((yieldval-yieldval1)/(yieldval1+1.d-99)>=tol)val1=yieldval-yieldval1

                 !Amount of ejected metals                                                                                                                   
                 Xeject(:n11,:n22)=yieldtab%Zeject(:n11,:n22)
                 yieldval=Xeject(ihx,ihy)*wxb*wyb   + Xeject(ihx,ihy+1)*wxb*wya   + &
                      Xeject(ihx+1,ihy)*wxa*wyb + Xeject(ihx+1,ihy+1)*wxa*wya
                 yieldval1=Xeject(ihx0,ihy)*wxb0*wyb   + Xeject(ihx0,ihy+1)*wxb0*wya   + &
                      Xeject(ihx0+1,ihy)*wxa0*wyb + Xeject(ihx0+1,ihy+1)*wxa0*wya
                 if((yieldval-yieldval1)/(yieldval1+1.d-99)<tol)ok_free=.false.
                 if((yieldval-yieldval1)/(yieldval1+1.d-99)>=tol)val2=yieldval-yieldval1

                 !Nb of SNII                                                                                                                                 
                 Xeject(:n11,:n22)=yieldtab%NSN2(:n11,:n22)
                 yieldval=Xeject(ihx,ihy)*wxb*wyb   + Xeject(ihx,ihy+1)*wxb*wya   + &
                      Xeject(ihx+1,ihy)*wxa*wyb + Xeject(ihx+1,ihy+1)*wxa*wya
                 yieldval1=Xeject(ihx0,ihy)*wxb0*wyb   + Xeject(ihx0,ihy+1)*wxb0*wya   + &
                      Xeject(ihx0+1,ihy)*wxa0*wyb + Xeject(ihx0+1,ihy+1)*wxa0*wya
                 if((yieldval-yieldval1)/(yieldval1+1.d-99)<tol)ok_free=.false.
                 if((yieldval-yieldval1)/(yieldval1+1.d-99)>=tol)valsn=yieldval-yieldval1

                 if (ok_free) then

                    !Amount of ejected chemical elements                                                                                                     
                    val3=0.
                    do ielt=1,nelt
                       Xeject(:n11,:n22)=yieldtab%Eeject(:n11,:n22,ielt)
                       yieldval=Xeject(ihx,ihy)*wxb*wyb   + Xeject(ihx,ihy+1)*wxb*wya   + &
                            Xeject(ihx+1,ihy)*wxa*wyb + Xeject(ihx+1,ihy+1)*wxa*wya
                       yieldval1=Xeject(ihx0,ihy)*wxb0*wyb   + Xeject(ihx0,ihy+1)*wxb0*wya   + &
                            Xeject(ihx0+1,ihy)*wxa0*wyb + Xeject(ihx0+1,ihy+1)*wxa0*wya
                              if((yieldval-yieldval1)/(yieldval1+1.d-99)>=tol) val3(ielt)=yieldval-yieldval1
                    enddo

                    iSN=iSN+1
                    xSN_tot(iSN,1)=xp(ipart,1)
                    xSN_tot(iSN,2)=xp(ipart,2)
                    xSN_tot(iSN,3)=xp(ipart,3)
                    vSN_tot(iSN,1)=vp(ipart,1)
                    vSN_tot(iSN,2)=vp(ipart,2)
                    vSN_tot(iSN,3)=vp(ipart,3)

                    ZSN_tot(iSN)=val2/val1
                    do ielt=1,nelt
                       ceSN_tot(iSN,ielt)=val3(ielt)/val1
                    enddo
                    mSN_tot(iSN)=val1*mp0(ipart)
                    NbSN_tot(iSN)=valsn*mp0(ipart)
                    mp(ipart)=mp(ipart)-mSN_tot(iSN)
                    indtab(ipart)=xx

                 endif        !ok_free 
              endif
              ipart=next_part  ! Go to next particle
           end do
        endif
        igrid=next(igrid)   ! Go to next grid
     end do
  deallocate(Xeject, val3, astarproper)
end subroutine sub2_kinetic_feedback
!################################################################
!################################################################
!################################################################
!################################################################
subroutine sub1_kinetic_feedback(icpu,kgrid,subnump,n11,n22,nSN_loc)
  use amr_commons
  use pm_commons
  use hydro_commons
  use cooling_module, ONLY: XH=>X, rhoc, mH 
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
  integer::nSN_tot_all,subnump,kgrid
  integer,dimension(1:ncpu)::nSN_icpu_all
  real(dp),dimension(:),allocatable::mSN_all,ZSN_all,NbSN_all
  real(dp),dimension(:,:),allocatable::xSN_all,vSN_all,ceSN_all
#endif
  !----------------------------------------------------------------------
  ! This subroutine compute the kinetic feedback due to SNII and
  ! imolement this using exploding GMC particles. 
  ! This routine is called only at coarse time step.
  ! Yohan Dubois
  !----------------------------------------------------------------------
  ! local constants
  integer::icpu,igrid,jgrid,npart1,npart2,ipart,jpart,next_part
  integer::nSN,nSN_loc,nSN_tot,info,iSN,ilevel,ivar,ielt,iSN2,iSN1,n11,n22,ihx, &
       ihy,ii,ihx0,i,nSN1,nSN2,nSN3,iSN3,indt1,indt2,nSN_tot2,npart22
  integer,dimension(1:ncpu)::nSN_icpu
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v,t0, &
       dxx,dyy,t1,xx,xx0,yy,ESN,t2,mass_load
  real::wxa,wxb,wya,wyb,wxa0,wxb0,rmin,rmax,dx1,yieldval,yieldval1, &
       val1,val2,valsn,zmin,zmax,dy1
  real,dimension(:),allocatable::val3
  real(dp)::current_time
  real(dp)::scale,dx_min,vol_min,nISM,nCOM,d0,mstar
  integer::nx_loc
!  integer,dimension(:),allocatable::ind_part,ind_grid
  logical ::ok_free
  integer ,dimension(:),allocatable::indSN,iSN_myid
  real(dp),dimension(:),allocatable::mSN,ZSN,vol_gas,ekBlast,NbSN,mloadSN,ZloadSN
  real(dp),dimension(:,:),allocatable::xSN,vSN,dq,ceSN,vloadSN,celoadSN
  real,dimension(:,:),allocatable::Xeject
  integer ,dimension(:),allocatable::indSN_tot,itemp
  real(dp),dimension(:),allocatable::mSN_tot,ZSN_tot,NbSN_tot
  real(dp),dimension(:,:),allocatable::xSN_tot,vSN_tot,ceSN_tot
  integer::isort,iastar
  real,dimension(:),allocatable::astarproper

  if(.not. hydro)return
  if(ndim.ne.3)return

  !if(verbose .and. icpu .eq. 1)write(*,*)'Entering make_sn'
  
  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

  igrid = kgrid

  ESN=eps_sn2*(1d51/2.d33)/scale_v**2

  n11=yieldtab%na
  n22=yieldtab%nz
  allocate(Xeject(1:n11,1:n22),val3(1:nelt),astarproper(1:n11))


!for test 3  
  do iastar=1,n11
     astarproper(iastar) = yieldtab%astar(iastar)*aexp**2*(1d9*365.*24.*3600./scale_t) !bring it from supercomoving to proper time in Gyr by multiplying by aexp**2
  end do

  rmin=log10(astarproper(1))
  rmax=log10(astarproper(n11))

  dx1=real(n11-1)/(rmax-rmin)

  zmin=log10(yieldtab%zstar(2))
  zmax=log10(yieldtab%zstar(n22))
  dy1=real(n22-2)/(zmax-zmin)

  !------------------------------------------------------                          
  ! Gather star particles in a SNII phase                                          
  !------------------------------------------------------
  nSN_loc=0
  ! Loop over levels
     do jgrid=1,subnump
        npart1=numbp(igrid)  ! Number of particles in the grid
        npart2=0
        ! Count old enough GMC particles
        if(npart1>0)then
           ipart=headp(igrid)
           ! Loop over particles
           do jpart=1,npart1
              ! Save next particle   <--- Very important !!!
              next_part=nextp(ipart)
              
              xx=texp-tpp(ipart)
              yy=zp(ipart)
              ihx=1+int(dx1*(log10(xx)-rmin))
              ihy=2+int(dy1*(log10(yy)-zmin))
              if (yy<yieldtab%zstar(2)) ihy=1

              if( tp(ipart)/=0. .and. xx>0. .and. ihx>0 )then

                 wxa=(xx-astarproper(ihx))  /(astarproper(ihx+1)-astarproper(ihx))
                 wxb=1.-wxa

                 wya=(yy-yieldtab%zstar(ihy))  /(yieldtab%zstar(ihy+1)-yieldtab%zstar(ihy))
                 wyb=1.-wya
                 
                 xx0=indtab(ipart)
                 if (xx0>0.d0) then
                    ihx0=1+int(dx1*(log10(xx0)-rmin))
                    wxa0=(xx0-astarproper(ihx0))    /(astarproper(ihx0+1)-astarproper(ihx0))
                    wxb0=1.-wxa0
                 else
                    ihx0=1
                    wxa0=0.
                    wxb0=0.
                 endif

                 ok_free=.true.

                 !Total ejected mass                                               
                 Xeject(:n11,:n22)=yieldtab%Meject(:n11,:n22)
                 yieldval=Xeject(ihx,ihy)*wxb*wyb   + Xeject(ihx,ihy+1)*wxb*wya + &
                      Xeject(ihx+1,ihy)*wxa*wyb + Xeject(ihx+1,ihy+1)*wxa*wya
                 yieldval1=Xeject(ihx0,ihy)*wxb0*wyb   + Xeject(ihx0,ihy+1)*wxb0*wya + &
                      Xeject(ihx0+1,ihy)*wxa0*wyb + Xeject(ihx0+1,ihy+1)*wxa0*wya
                 if((yieldval-yieldval1)/(yieldval1+1.d-99)<tol)ok_free=.false.
                 if((yieldval-yieldval1)/(yieldval1+1.d-99)>=tol)val1=yieldval-yieldval1

                 !Amount of ejected metals                                         
                 Xeject(:n11,:n22)=yieldtab%Zeject(:n11,:n22)
                 yieldval=Xeject(ihx,ihy)*wxb*wyb   + Xeject(ihx,ihy+1)*wxb*wya + &
                      Xeject(ihx+1,ihy)*wxa*wyb + Xeject(ihx+1,ihy+1)*wxa*wya
                 yieldval1=Xeject(ihx0,ihy)*wxb0*wyb   + Xeject(ihx0,ihy+1)*wxb0*wya + &
                      Xeject(ihx0+1,ihy)*wxa0*wyb + Xeject(ihx0+1,ihy+1)*wxa0*wya
                 if((yieldval-yieldval1)/(yieldval1+1.d-99)<tol)ok_free=.false.

                 !Nb of SNII                                                       
                 Xeject(:n11,:n22)=yieldtab%NSN2(:n11,:n22)
                 yieldval=Xeject(ihx,ihy)*wxb*wyb   + Xeject(ihx,ihy+1)*wxb*wya + &
                      Xeject(ihx+1,ihy)*wxa*wyb + Xeject(ihx+1,ihy+1)*wxa*wya
                 yieldval1=Xeject(ihx0,ihy)*wxb0*wyb   + Xeject(ihx0,ihy+1)*wxb0*wya   + &
                      Xeject(ihx0+1,ihy)*wxa0*wyb + Xeject(ihx0+1,ihy+1)*wxa0*wya
                 if((yieldval-yieldval1)/(yieldval1+1.d-99)<tol)ok_free=.false.

                 if (ok_free) npart2=npart2+1
              endif
              ipart=next_part  ! Go to next particle
            end do
        endif
        nSN_loc=nSN_loc+npart2   ! Add SNe to the total
        igrid=next(igrid)   ! Go to next grid
     end do
! end do
  ! End loop over levels
  deallocate(Xeject, val3, astarproper)
end subroutine sub1_kinetic_feedback

!################################################################
!################################################################
!################################################################
!################################################################
subroutine kinetic_feedback
  use amr_commons
  use pm_commons
  use hydro_commons
  use cooling_module, ONLY: XH=>X, rhoc, mH 
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
  integer::nSN_tot_all
  integer,dimension(1:ncpu)::nSN_icpu_all
  real(dp),dimension(:),allocatable::mSN_all,ZSN_all,NbSN_all
  real(dp),dimension(:,:),allocatable::xSN_all,vSN_all,ceSN_all
  real(dp),dimension(:),allocatable::SN_allreduce1, SN_allreduce2
#endif
  !----------------------------------------------------------------------
  ! This subroutine compute the kinetic feedback due to SNII and
  ! imolement this using exploding GMC particles. 
  ! This routine is called only at coarse time step.
  ! Yohan Dubois
  !----------------------------------------------------------------------
  ! local constants
  integer::icpu,igrid,jgrid,npart1,npart2,ipart,jpart,next_part
  integer::nSN,nSN_loc,nSN_tot,info,iSN,ilevel,ivar,ielt,iSN2,iSN1,n11,n22,ihx, &
       ihy,ii,ihx0,i,nSN1,nSN2,nSN3,iSN3,indt1,indt2,nSN_tot2,npart22
  integer:: jSN,kSN
  integer,dimension(1:ncpu)::nSN_icpu
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v,t0, &
       dxx,dyy,t1,xx,xx0,yy,ESN,t2,mass_load
  real::wxa,wxb,wya,wyb,wxa0,wxb0,rmin,rmax,dx1,yieldval,yieldval1, &
       val1,val2,valsn,zmin,zmax,dy1
  real,dimension(:),allocatable::val3
  real(dp)::current_time
  real(dp)::scale,dx_min,vol_min,nISM,nCOM,d0,mstar
  integer::nx_loc
!  integer,dimension(:),allocatable::ind_part,ind_grid
  logical ::ok_free
  integer ,dimension(:),allocatable::indSN,iSN_myid
  real(dp),dimension(:),allocatable::mSN,ZSN,vol_gas,ekBlast,NbSN,mloadSN,ZloadSN
  real(dp),dimension(:,:),allocatable::xSN,vSN,dq,ceSN,vloadSN,celoadSN
  real,dimension(:,:),allocatable::Xeject
  integer ,dimension(:),allocatable::indSN_tot,itemp
  real(dp),dimension(:),allocatable::mSN_tot,ZSN_tot,NbSN_tot
  real(dp),dimension(:,:),allocatable::xSN_tot,vSN_tot,ceSN_tot
  integer::isort,iastar
  real,dimension(:),allocatable::astarproper
  integer mythread, nthreads, icount,jcount, subnump
  common /kinetic_fb/ mythread,jSN
!$omp threadprivate(/kinetic_fb/)
  integer, dimension(:), allocatable:: nparticles, ptrhead,mynewSN
  integer:: start_mynewSN 

  if(.not. hydro)return
  if(ndim.ne.3)return

!$omp parallel
  mythread = omp_get_thread_num()
  if(mythread .eq.0) nthreads = omp_get_num_threads()
!$omp end parallel
  allocate(ptrhead(0:nthreads-1), nparticles(0:nthreads-1), mynewSN(0:nthreads-1))
  mynewSN = 0

  if(verbose .and. myid .eq. 1) write(*,*)'Entering make_sn'
  
  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

  ESN=eps_sn2*(1d51/2.d33)/scale_v**2

  n11=yieldtab%na
  n22=yieldtab%nz
! allocate(Xeject(1:n11,1:n22),val3(1:nelt),astarproper(1:n11))
  allocate(astarproper(1:n11))


!for test 3  
  do iastar=1,n11
     astarproper(iastar) = yieldtab%astar(iastar)*aexp**2*(1d9*365.*24.*3600./scale_t) !bring it from supercomoving to proper time in Gyr by multiplying by aexp**2
  end do

  rmin=log10(astarproper(1))
  rmax=log10(astarproper(n11))

  dx1=real(n11-1)/(rmax-rmin)

  zmin=log10(yieldtab%zstar(2))
  zmax=log10(yieldtab%zstar(n22))
  dy1=real(n22-2)/(zmax-zmin)

  !------------------------------------------------------                          
  ! Gather star particles in a SNII phase                                          
  !------------------------------------------------------
  nSN_loc=0
  ! Loop over levels
  do icpu=1,ncpu
        if(numbl(icpu,levelmin) .le.0) goto 13
        call pthreadLinkedList(headl(icpu,levelmin),numbl(icpu,levelmin),nthreads,nparticles,ptrhead,next)
!$omp parallel private(subnump,igrid) reduction(+:nSN_loc)
        subnump = nparticles(mythread)
        igrid = ptrhead(mythread)
        call sub1_kinetic_feedback(icpu,igrid,subnump,n11,n22,nSN_loc)
        mynewSN(mythread) = mynewSN(mythread) + nSN_loc
!$omp end parallel
13      continue
  end do
  ! End loop over levels
  nSN_icpu=0
  nSN_icpu(myid)=nSN_loc
#ifndef WITHOUTMPI
  ! Give an array of number of SN on each cpu available to all cpus
  call MPI_ALLREDUCE(nSN_icpu,nSN_icpu_all,ncpu,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
  nSN_icpu=nSN_icpu_all
#endif

  nSN_tot=sum(nSN_icpu(1:ncpu))

  if (nSN_tot .eq. 0) return

  if(myid==1)then
     write(*,*)'-----------------------------------------------'
     write(*,*)'Number of particles with SNII ',nSN_tot
     write(*,*)'-----------------------------------------------'
  endif

  ! Allocate arrays for the position and the mass of the SN
  allocate(xSN_tot(1:nSN_tot,1:3),vSN_tot(1:nSN_tot,1:3))
  allocate(mSN_tot(1:nSN_tot),ZSN_tot(1:nSN_tot),ceSN_tot(1:nSN_tot,1:nelt),NbSN_tot(1:nSN_tot),itemp(1:nSN_tot))
  xSN_tot=0.;vSN_tot=0.;mSN_tot=0.;ZSN_tot=0.;ceSN_tot=0.;NbSN_tot=0.

  !------------------------------------------------------
  ! Store position and mass of the GMC into the SN array
  !------------------------------------------------------
  if(myid==1)then
     iSN=0
  else
     iSN=sum(nSN_icpu(1:myid-1))
  endif

!$omp parallel
  if(mythread==0) then
     jSN = iSN
  else
     jSN = iSN + sum(mynewSN(0:mythread-1))
  endif
!$omp end parallel

  ! Loop over levels
  do icpu=1,ncpu
     if(numbl(icpu,levelmin) .le.0) goto 14
     call pthreadLinkedList(headl(icpu,levelmin),numbl(icpu,levelmin),nthreads,nparticles,ptrhead,next)
!$omp parallel private(subnump,igrid) 
     subnump = nparticles(mythread)
     igrid = ptrhead(mythread)
     call sub2_kinetic_feedback(icpu,igrid,subnump,n11,n22,jSN,nSN_tot,xSN_tot,vSN_tot,mSN_tot,ZSN_tot,NbSN_tot,ceSN_tot)
!$omp end parallel
14   continue
  end do 
  ! End loop over levels

  deallocate(nparticles,ptrhead)
  deallocate(mynewSN)

#ifndef WITHOUTMPI
  if(nSN_tot*(3+3+2+1+nelt) .gt. 1000000000 .or. nid .lt. 100) then
    allocate(xSN_all(1:nSN_tot,1:3),vSN_all(1:nSN_tot,1:3),mSN_all(1:nSN_tot),ZSN_all(1:nSN_tot))
    allocate(ceSN_all(1:nSN_tot,1:nelt),NbSN_all(1:nSN_tot))
    call MPI_ALLREDUCE(xSN_tot,xSN_all,nSN_tot*3,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
    call MPI_ALLREDUCE(vSN_tot,vSN_all,nSN_tot*3,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
    call MPI_ALLREDUCE(mSN_tot,mSN_all,nSN_tot  ,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
    call MPI_ALLREDUCE(ZSN_tot,ZSN_all,nSN_tot  ,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
    call MPI_ALLREDUCE(NbSN_tot,NbSN_all,nSN_tot,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
    call MPI_ALLREDUCE(ceSN_tot,ceSN_all,nSN_tot*nelt,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
    xSN_tot=xSN_all
    vSN_tot=vSN_all
    mSN_tot=mSN_all
    ZSN_tot=ZSN_all
    NbSN_tot=NbSN_all
    ceSN_tot=ceSN_all
    deallocate(xSN_all,vSN_all,mSN_all,ZSN_all,NbSN_all,ceSN_all)
  else
    allocate(SN_allreduce1(1:nSN_tot*(3+3+2+1+nelt))
    allocate(SN_allreduce2(1:nSN_tot*(3+3+2+1+nelt))
    SN_allreduce1(1:nSN_tot*3) = xSN_tot(1:nSN_tot*3);
    SN_allreduce1(nSN_tot*3+1:nSN_tot*6) = vSN_tot(1:nSN_tot*3);
    SN_allreduce1(nSN_tot*6+1:nSN_tot*(6+1)) = mSN_tot(1:nSN_tot);
    SN_allreduce1(nSN_tot*(6+1)+1:nSN_tot*(6+2)) = ZSN_tot(1:nSN_tot);
    SN_allreduce1(nSN_tot*(6+2)+1:nSN_tot*(6+3)) = NbSN_tot(1:nSN_tot);
    SN_allreduce1(nSN_tot*(6+3)+1:nSN_tot*(6+3+nelt)) = ceSN_tot(1:nSN_tot*nelt);
    call MPI_ALLREDUCE(SN_allreduce1, SN_allreduce2, nSN_tot*(3+3+2+1+nelt), &
                     MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
    xSN_tot(1:nSN_tot*3)     = SN_allreduce2(            1:nSN_tot*3       )
    vSN_tot(1:nSN_tot*3)     = SN_allreduce2(nSN_tot*3 + 1:nSN_tot*6       )
    mSN_tot(1:nSN_tot*1)     = SN_allreduce2(nSN_tot*6 + 1:nSN_tot*7       )
    ZSN_tot(1:nSN_tot*1)     = SN_allreduce2(nSN_tot*7 + 1:nSN_tot*8       )
    NbSN_tot(1:nSN_tot*1)    = SN_allreduce2(nSN_tot*8 + 1:nSN_tot*9       )
    ceSN_tot(1:nSN_tot*nelt) = SN_allreduce2(nSN_tot*9 + 1:nSN_tot*(9+nelt))
    deallocate(SN_allreduce1, SN_allreduce2)
  endif
#endif

  call getSNonmyid(itemp,nSN,xSN_tot,nSN_tot)

  ! Allocate the arrays for the position and the mass of the SN                                                                                              
  allocate(xSN(1:nSN,1:3),vSN(1:nSN,1:3),mSN(1:nSN),ZSN(1:nSN),NbSN(1:nSN),ceSN(1:nSN,1:nelt),iSN_myid(1:nSN))
  xSN=0d0; vSN=0d0; mSN=0d0; ZSN=0d0; NbSN=0d0; ceSN=0d0; iSN_myid=0

!$omp parallel do private(iSN,isort,ielt)
  do iSN=1,nSN
     isort=itemp(iSN)
     iSN_myid(iSN)=isort
     xSN(iSN,1)=xSN_tot(isort,1)
     xSN(iSN,2)=xSN_tot(isort,2)
     xSN(iSN,3)=xSN_tot(isort,3)
     vSN(iSN,1)=vSN_tot(isort,1)
     vSN(iSN,2)=vSN_tot(isort,2)
     vSN(iSN,3)=vSN_tot(isort,3)
     mSN(iSN)  =mSN_tot(isort)
     ZSN(iSN)  =ZSN_tot(isort)
     NbSN(iSN) =NbSN_tot(isort)
     do ielt=1,nelt
        ceSN(iSN,ielt) =ceSN_tot(isort,ielt)
     enddo
  enddo
  deallocate(xSN_tot,vSN_tot,mSN_tot,ZSN_tot,NbSN_tot,ceSN_tot,itemp,astarproper)

  allocate(vol_gas(1:nSN),dq(1:nSN,1:3),ekBlast(1:nSN))
  allocate(indSN(1:nSN))
  allocate(mloadSN(1:nSN),ZloadSN(1:nSN),celoadSN(1:nSN,1:nelt),vloadSN(1:nSN,1:3))

  ! Compute the grid discretization effects
  call average_SN(xSN,vSN,NbSN,vol_gas,dq,ekBlast,indSN,nSN,nSN_tot,iSN_myid,mSN,mloadSN,ZSN,ZloadSN,ceSN,celoadSN,vloadSN)
  ! Modify hydro quantities to account for a Sedov blast wave
  call Sedov_blast(xSN,mSN,NbSN,indSN,vol_gas,dq,ekBlast,nSN,mloadSN,ZloadSN,celoadSN,vloadSN)
  deallocate(xSN,vSN,mSN,ZSN,indSN,vol_gas,dq,ekBlast,ceSN,NbSN,iSN_myid)
  deallocate(mloadSN,ZloadSN,celoadSN,vloadSN)

  ! Update hydro quantities for split cells
  do ilevel=nlevelmax,levelmin,-1
     call upload_fine(ilevel)
     do ivar=1,nvar
        call make_virtual_fine_dp(uold(1,ivar),ilevel)
     enddo
  enddo
end subroutine kinetic_feedback
!################################################################
!################################################################
!################################################################
!################################################################
subroutine average_SN(xSN,vSN,NbSN,vol_gas,dq,ekBlast,ind_blast,nSN,nSN_tot,iSN_myid,mSN,mloadSN,ZSN,ZloadSN,ceSN,celoadSN,vloadSN)
  use pm_commons
  use amr_commons
  use hydro_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  !------------------------------------------------------------------------
  ! This routine average the hydro quantities inside the SN bubble
  !------------------------------------------------------------------------
  integer::ilevel,ncache,nSN,j,iSN,ind,ix,iy,iz,ngrid,iskip,nSN_tot,ielt
  integer::i,nx_loc,igrid,info,jSN,kSN,iSNsize
  integer,dimension(1:nvector)::ind_grid,ind_cell
  real(dp)::x,y,z,dr_SN,d,u,v,w,ek,u2,v2,w2,dr_cell
  real(dp)::scale,dx,dxx,dyy,dzz,dx_min,dx_loc,vol_loc,rmax2,rmax
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  real(dp)::eint,ekk,ekk1,ekk2,mload,Zload,ceload
  real(dp),dimension(1:3)::skip_loc
  real(dp),dimension(1:twotondim,1:3)::xc
  integer ,dimension(1:nSN)::ind_blast
  real(dp),dimension(1:nSN)::mSN,m_gas,vol_gas,ekBlast,ZSN,NbSN,mloadSN,ZloadSN
  real(dp),dimension(1:nSN,1:3)::xSN,vSN,dq,u2Blast,vloadSN
  real(dp),dimension(1:nSN,1:nelt)::celoadSN,ceSN
#ifndef WITHOUTMPI
! real(dp),dimension(1:nSN_tot)::vol_gas_all,ekBlast_all
! real(dp),dimension(1:nSN_tot,1:3)::dq_all,u2Blast_all,vloadSN_all
! real(dp),dimension(1:nSN_tot)::vol_gas_mpi,ekBlast_mpi
! real(dp),dimension(1:nSN_tot,1:3)::dq_mpi,u2Blast_mpi,vloadSN_mpi
! real(dp),dimension(1:nSN_tot)::mloadSN_mpi,mloadSN_all,ZloadSN_mpi,ZloadSN_all
! real(dp),dimension(1:nSN_tot,1:nelt)::celoadSN_mpi,celoadSN_all
  real(dp),dimension(:), allocatable::vol_gas_all,ekBlast_all
  real(dp),dimension(:,:), allocatable::dq_all,u2Blast_all,vloadSN_all
  real(dp),dimension(:), allocatable::vol_gas_mpi,ekBlast_mpi
  real(dp),dimension(:,:), allocatable::dq_mpi,u2Blast_mpi,vloadSN_mpi
  real(dp),dimension(:), allocatable::mloadSN_mpi,mloadSN_all,ZloadSN_mpi,ZloadSN_all
  real(dp),dimension(:,:), allocatable::celoadSN_mpi,celoadSN_all
  real(dp),dimension(:,:), allocatable::commall1, commall2
#endif
  logical ,dimension(1:nvector)::ok
  integer ,dimension(1:nSN)::iSN_myid
  integer::ind_SN
  integer mythread, nthreads
  common /average_SN_common/ mythread
!$omp threadprivate(/average_SN_common/)

  if(verbose .and. myid .eq. 1)write(*,*)'Entering average_SN'
!$omp parallel
  mythread = omp_get_thread_num()
  if(mythread.eq.0) nthreads = omp_get_num_threads()
!$omp end parallel

  allocate(vol_gas_all(1:nSN_tot),ekBlast_all(1:nSN_tot))
  allocate(dq_all(1:nSN_tot,1:3),u2Blast_all(1:nSN_tot,1:3),vloadSN_all(1:nSN_tot,1:3))
  allocate(vol_gas_mpi(1:nSN_tot),ekBlast_mpi(1:nSN_tot))
  allocate(dq_mpi(1:nSN_tot,1:3),u2Blast_mpi(1:nSN_tot,1:3),vloadSN_mpi(1:nSN_tot,1:3))
  allocate(mloadSN_mpi(1:nSN_tot),mloadSN_all(1:nSN_tot))
  allocate(ZloadSN_mpi(1:nSN_tot),ZloadSN_all(1:nSN_tot))
  allocate(celoadSN_mpi(1:nSN_tot,1:nelt),celoadSN_all(1:nSN_tot,1:nelt))


  ! Mesh spacing in that level
  nx_loc=(icoarse_max-icoarse_min+1)
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  skip_loc(1)=dble(icoarse_min)
  skip_loc(2)=dble(jcoarse_min)
  skip_loc(3)=dble(kcoarse_min)
  scale=boxlen/dble(nx_loc)
  dx_min=scale*0.5D0**nlevelmax

  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

  ! Maximum radius of the ejecta
  rmax=MAX(rcell*dx_min*scale_l/aexp,rbubble*3.08d18)
  rmax=rmax/scale_l
  rmax2=rmax*rmax

  ! Initialize the averaged variables
  vol_gas=0.0;dq=0.0;u2Blast=0.0;ekBlast=0.0;ind_blast=-1
  mloadSN=0.0;ZloadSN=0.0;celoadSN=0.0;vloadSN=0.0


  ! Loop over levels
  do ilevel=levelmin,nlevelmax
     ! Computing local volume (important for averaging hydro quantities) 
     dx=0.5D0**ilevel 
     dx_loc=dx*scale
     vol_loc=dx_loc**ndim
     ! Cells center position relative to grid center position
     do ind=1,twotondim  
        iz=(ind-1)/4
        iy=(ind-1-4*iz)/2
        ix=(ind-1-2*iy-4*iz)
        xc(ind,1)=(dble(ix)-0.5D0)*dx
        xc(ind,2)=(dble(iy)-0.5D0)*dx
        xc(ind,3)=(dble(iz)-0.5D0)*dx
     end do

     ! Loop over grids
     ncache=active(ilevel)%ngrid
!$omp parallel private(iSN,igrid,ngrid,i,ind_grid,ind,iskip,ind_cell,ok,x,y,z,&
!$omp       dxx,dyy,dzz,dr_SN,dr_cell, &
!$omp       d,u,v,w,ekk,eint,mload,Zload,ceload,ielt,jSN,kSN,iSNsize)
     iSNsize = (nSN + nthreads-1)/nthreads
	 jSN = iSNsize*mythread + 1
	 kSN = iSNsize*(mythread+1) 
	 kSN = min(kSN,nSN)
     do igrid=1,ncache,nvector
        ngrid=MIN(nvector,ncache-igrid+1)
        do i=1,ngrid
           ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
        end do
        !    Loop over cells
        do ind=1,twotondim  
           iskip=ncoarse+(ind-1)*ngridmax
           do i=1,ngrid
              ind_cell(i)=iskip+ind_grid(i)
           end do
           ! Flag leaf cells
           do i=1,ngrid
              ok(i)=son(ind_cell(i))==0
           end do

           do i=1,ngrid
              if(ok(i))then
                 ! Get gas cell position
                 x=(xg(ind_grid(i),1)+xc(ind,1)-skip_loc(1))*scale
                 y=(xg(ind_grid(i),2)+xc(ind,2)-skip_loc(2))*scale
                 z=(xg(ind_grid(i),3)+xc(ind,3)-skip_loc(3))*scale
			     do iSN=jSN,kSN
                 	! Check if the cell lies within the SN radius
                    dxx=x-xSN(iSN,1)
                    dyy=y-xSN(iSN,2)
                    dzz=z-xSN(iSN,3)
                    dr_SN=dxx**2+dyy**2+dzz**2
                    dr_cell=MAX(ABS(dxx),ABS(dyy),ABS(dzz))
                    if(dr_SN.lt.rmax2)then
                       vol_gas(iSN)=vol_gas(iSN)+vol_loc
                       ! Take account for grid effects on the conservation of the
                       ! normalized linear momentum
                       u=dxx/rmax
                       v=dyy/rmax
                       w=dzz/rmax
                       ! Add the local normalized linear momentum to the total linear
                       ! momentum of the blast wave (should be zero with no grid effect)
                       dq(iSN,1)=dq(iSN,1)+u*vol_loc
                       dq(iSN,2)=dq(iSN,2)+v*vol_loc
                       dq(iSN,3)=dq(iSN,3)+w*vol_loc
                       u2Blast(iSN,1)=u2Blast(iSN,1)+u*u*vol_loc
                       u2Blast(iSN,2)=u2Blast(iSN,2)+v*v*vol_loc
                       u2Blast(iSN,3)=u2Blast(iSN,3)+w*w*vol_loc
                    endif
                    if(dr_cell.le.dx_loc/2.0)then !same as above?
                       ind_blast(iSN)=ind_cell(i)
                       ekBlast  (iSN)=vol_loc
                       d=uold(ind_blast(iSN),1)                                                                     
                       u=uold(ind_blast(iSN),2)/d                                                                   
                       v=uold(ind_blast(iSN),3)/d                                                                   
                       w=uold(ind_blast(iSN),4)/d                                                                   
                       ekk=0.5d0*d*(u*u+v*v+w*w)                                                                    
                       eint=uold(ind_blast(iSN),5)-ekk                                                              
                       ! Mass loading factor of the Sedov explosion                                                 
                       ! Ensure that no more that 25% of the gas content is removed                                 
                       mload=min(f_w*mSN(iSN),0.25d0*d*vol_loc)                                                     
                       mloadSN(iSN)=mSN(iSN)+mload                                                                  
                       ! Update gas mass and metal content in the cell                                              
                       if(metal)then                                                                                
                          Zload=uold(ind_blast(iSN),imetal)/d                                                       
                          ZloadSN(iSN)=( mload*Zload + ZSN(iSN)*mSN(iSN) ) / mloadSN(iSN)                           
                          uold(ind_blast(iSN),imetal)=uold(ind_blast(iSN),imetal)-Zload*mload/vol_loc
                       endif                                                                                        
                       do ielt=1,nelt                                                                               
                          ceload=uold(ind_blast(iSN),ichem+ielt-1)/d                                                 
                          celoadSN(iSN,ielt)=( mload*ceload + ceSN(iSN,ielt)*mSN(iSN) ) / mloadSN(iSN)              
                          uold(ind_blast(iSN),ichem+ielt-1)=uold(ind_blast(iSN),ichem+ielt-1)-ceload*mload/vol_loc
                       enddo                                                                                        
                       d=uold(ind_blast(iSN),1)-mload/vol_loc                                                       
                                                                                                                    
                       uold(ind_blast(iSN),1)=d                                                                     
                       uold(ind_blast(iSN),2)=d*u                                                                   
                       uold(ind_blast(iSN),3)=d*v                                                                   
                       uold(ind_blast(iSN),4)=d*w                                                                   
                       uold(ind_blast(iSN),5)=eint+0.5d0*d*(u*u+v*v+w*w)                                            

                       vloadSN(iSN,1)=(mSN(iSN)*vSN(iSN,1)+mload*u)/mloadSN(iSN)
                       vloadSN(iSN,2)=(mSN(iSN)*vSN(iSN,2)+mload*v)/mloadSN(iSN)
                       vloadSN(iSN,3)=(mSN(iSN)*vSN(iSN,3)+mload*w)/mloadSN(iSN)
                    endif
       			end do ! End loop over nSN
              endif
           end do
           
        end do
           ! End loop over cells
     end do
     ! End loop over grids
!$omp end parallel
  end do
  ! End loop over levels

#ifndef WITHOUTMPI
  vol_gas_mpi=0d0; dq_mpi=0d0; u2Blast_mpi=0d0
  mloadSN_mpi=0d0; ZloadSN_mpi=0d0; celoadSN_mpi=0d0; vloadSN_mpi=0d0
  ! Put the nSN size arrays into nSN_tot size arrays to synchronize processors                                
  if(nSN .gt. 1000000000 .or. nid .lt. 100) then
!$omp parallel do private(iSN,ind_SN,ielt)
     do iSN=1,nSN
        ind_SN=iSN_myid(iSN)
        vol_gas_mpi(ind_SN)=vol_gas(iSN)
        mloadSN_mpi(ind_SN)=mloadSN(iSN)
        ZloadSN_mpi(ind_SN)=ZloadSN(iSN)
        do ielt=1,nelt
           celoadSN_mpi(ind_SN,ielt)=celoadSN(iSN,ielt)
        enddo
        vloadSN_mpi(ind_SN,1)=vloadSN(iSN,1)
        vloadSN_mpi(ind_SN,2)=vloadSN(iSN,2)
        vloadSN_mpi(ind_SN,3)=vloadSN(iSN,3)
        dq_mpi     (ind_SN,1)=dq     (iSN,1)
        dq_mpi     (ind_SN,2)=dq     (iSN,2)
        dq_mpi     (ind_SN,3)=dq     (iSN,3)
        u2Blast_mpi(ind_SN,1)=u2Blast(iSN,1)
        u2Blast_mpi(ind_SN,2)=u2Blast(iSN,2)
        u2Blast_mpi(ind_SN,3)=u2Blast(iSN,3)
     enddo
     call MPI_ALLREDUCE(vol_gas_mpi,vol_gas_all,nSN_tot  ,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
     call MPI_ALLREDUCE(mloadSN_mpi,mloadSN_all,nSN_tot  ,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
     call MPI_ALLREDUCE(ZloadSN_mpi,ZloadSN_all,nSN_tot  ,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
     call MPI_ALLREDUCE(vloadSN_mpi,vloadSN_all,nSN_tot*3,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
     call MPI_ALLREDUCE(dq_mpi     ,dq_all     ,nSN_tot*3,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
     call MPI_ALLREDUCE(u2Blast_mpi,u2Blast_all,nSN_tot*3,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
     call MPI_ALLREDUCE(celoadSN_mpi,celoadSN_all,nSN_tot*nelt,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
     vol_gas_mpi=vol_gas_all
     mloadSN_mpi=mloadSN_all
     ZloadSN_mpi=ZloadSN_all
     vloadSN_mpi=vloadSN_all
     dq_mpi     =dq_all
     u2Blast_mpi=u2Blast_all
     celoadSN_mpi=celoadSN_all
  ! Put the nSN_tot size arrays into nSN size arrays                                                          
!$omp parallel do private(iSN,ind_SN,ielt)
     do iSN=1,nSN
        ind_SN=iSN_myid(iSN)
        vol_gas(iSN)=vol_gas_mpi(ind_SN)
        mloadSN(iSN)=mloadSN_mpi(ind_SN)
        ZloadSN(iSN)=ZloadSN_mpi(ind_SN)
        vloadSN(iSN,1)=vloadSN_mpi(ind_SN,1)
        vloadSN(iSN,2)=vloadSN_mpi(ind_SN,2)
        vloadSN(iSN,3)=vloadSN_mpi(ind_SN,3)
        dq     (iSN,1)=dq_mpi     (ind_SN,1)
        dq     (iSN,2)=dq_mpi     (ind_SN,2)
        dq     (iSN,3)=dq_mpi     (ind_SN,3)
        u2Blast(iSN,1)=u2Blast_mpi(ind_SN,1)
        u2Blast(iSN,2)=u2Blast_mpi(ind_SN,2)
        u2Blast(iSN,3)=u2Blast_mpi(ind_SN,3)
        do ielt=1,nelt
           celoadSN(iSN,ielt)=celoadSN_mpi(ind_SN,ielt)
        enddo
     enddo
  else
     allocate(commall1(1:nSN*(9+nelt+3)))
     allocate(commall2(1:nSN*(9+nelt+3)))
!$omp parallel do private(iSN,ind_SN,ielt)
     do iSN=1,nSN
        ind_SN=iSN_myid(iSN)
        commall1(ind_SN)=vol_gas(iSN)
        commall1(nSN+ind_SN)=mloadSN(iSN)
        commall1(nSN*2+ind_SN)=ZloadSN(iSN)
        commall1(nSN*3+ind_SN)=vloadSN(iSN,1)
        commall1(nSN*4+ind_SN)=vloadSN(iSN,2)
        commall1(nSN*5+ind_SN)=vloadSN(iSN,3)
        commall1(nSN*6+ind_SN)=dq     (iSN,1)
        commall1(nSN*7+ind_SN)=dq     (iSN,2)
        commall1(nSN*8+ind_SN)=dq     (iSN,3)
        commall1(nSN*9+ind_SN)=u2Blast(iSN,1)
        commall1(nSN*10+ind_SN)=u2Blast(iSN,2)
        commall1(nSN*11+ind_SN)=u2Blast(iSN,3)
        do ielt=1,nelt
           commall1(nSN*(11+ielt)+ind_SN)=celoadSN(iSN,ielt)
        enddo
     enddo
     call MPI_ALLREDUCE(commall1,commall2,nSN_tot*(12+nelt),MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
!$omp parallel do private(iSN,ind_SN,ielt)
     do iSN=1,nSN
        ind_SN=iSN_myid(iSN)
        vol_gas(iSN)=commall2(ind_SN)
        mloadSN(iSN)=commall2(nSN+ind_SN)
        ZloadSN(iSN)=commall2(nSN*2+ind_SN)
        vloadSN(iSN,1)=commall2(nSN*3+ind_SN)
        vloadSN(iSN,2)=commall2(nSN*4+ind_SN)
        vloadSN(iSN,3)=commall2(nSN*5+ind_SN)
        dq     (iSN,1)=commall2(nSN*6+ind_SN)
        dq     (iSN,2)=commall2(nSN*7+ind_SN)
        dq     (iSN,3)=commall2(nSN*8+ind_SN)
        u2Blast(iSN,1)=commall2(nSN*9+ind_SN)
        u2Blast(iSN,2)=commall2(nSN*10+ind_SN)
        u2Blast(iSN,3)=commall2(nSN*11+ind_SN)
        do ielt=1,nelt
           celoadSN(iSN,ielt)=commall2(nSN*(11+ielt)+ind_SN)
        enddo
     enddo

     deallocate(commall1,commall2)
     
  endif
#endif
!$omp parallel do private(iSN,u2,v2,w2)
  do iSN=1,nSN
     if(vol_gas(iSN)>0d0)then
        dq(iSN,1)=dq(iSN,1)/vol_gas(iSN)
        dq(iSN,2)=dq(iSN,2)/vol_gas(iSN)
        dq(iSN,3)=dq(iSN,3)/vol_gas(iSN)
        u2Blast(iSN,1)=u2Blast(iSN,1)/vol_gas(iSN)
        u2Blast(iSN,2)=u2Blast(iSN,2)/vol_gas(iSN)
        u2Blast(iSN,3)=u2Blast(iSN,3)/vol_gas(iSN)
        u2=u2Blast(iSN,1)-dq(iSN,1)**2
        v2=u2Blast(iSN,2)-dq(iSN,2)**2
        w2=u2Blast(iSN,3)-dq(iSN,3)**2
        ekBlast(iSN)=max(0.5d0*(u2+v2+w2),0.0d0)
     endif
  end do


  deallocate(vol_gas_all,ekBlast_all,dq_all,u2Blast_all,vloadSN_all,vol_gas_mpi,ekBlast_mpi)
  deallocate(dq_mpi,u2Blast_mpi,vloadSN_mpi,mloadSN_mpi,mloadSN_all,ZloadSN_mpi,ZloadSN_all)
  deallocate(celoadSN_mpi,celoadSN_all)

  if(verbose .and. myid .eq. 1)write(*,*)'Exiting average_SN'

end subroutine average_SN
!################################################################
!################################################################
!################################################################
!################################################################
subroutine Sedov_blast(xSN,mSN,NbSN,indSN,vol_gas,dq,ekBlast,nSN,mloadSN,ZloadSN,celoadSN,vloadSN)
  use pm_commons
  use amr_commons
  use hydro_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  !------------------------------------------------------------------------
  ! This routine merges SN using the FOF algorithm.
  !------------------------------------------------------------------------
  integer::ilevel,j,iSN,nSN,ind,ix,iy,iz,ngrid,iskip,ielt
  integer::i,nx_loc,igrid,info,ncache
  integer,dimension(1:nvector)::ind_grid,ind_cell
  real(dp)::x,y,z,dx,dxx,dyy,dzz,dr_SN,d,u,v,w,ek,u_r,ESN,d_gas
  real(dp)::scale,dx_min,dx_loc,vol_loc,rmax2,rmax
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  real(dp),dimension(1:3)::skip_loc
  real(dp),dimension(1:twotondim,1:3)::xc
  real(dp),dimension(1:nSN)::mSN,ZloadSN,p_gas,vol_gas,uSedov,ekBlast,NbSN,mloadSN
  real(dp),dimension(1:nSN,1:3)::xSN,vloadSN,dq
  real(dp),dimension(1:nSN,1:nelt)::celoadSN
  integer ,dimension(1:nSN)::indSN
  logical ,dimension(1:nvector)::ok

  if(verbose .and. myid .eq. 1)write(*,*)'Entering Sedov_blast'

  ! Mesh spacing in that level
  nx_loc=(icoarse_max-icoarse_min+1)
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  skip_loc(1)=dble(icoarse_min)
  skip_loc(2)=dble(jcoarse_min)
  skip_loc(3)=dble(kcoarse_min)
  scale=boxlen/dble(nx_loc)
  dx_min=scale*0.5D0**nlevelmax

  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

  ! Maximum radius of the ejecta
  rmax=MAX(rcell*dx_min*scale_l/aexp,rbubble*3.08d18)
  rmax=rmax/scale_l
  rmax2=rmax*rmax
  
  ! Supernova specific energy from cgs to code units
  ESN=(eps_sn2/(10d0*2d33))/scale_v**2

!$omp parallel do private(iSN)
  do iSN=1,nSN
     if(vol_gas(iSN)>0d0)then
        if(ekBlast(iSN)==0d0)then
           p_gas(iSN)=NbSN(iSN)*ESN/vol_gas(iSN)
           uSedov(iSN)=0d0
        else
           p_gas(iSN)=(1d0-f_ek)*NbSN(iSN)*ESN/vol_gas(iSN)
           uSedov(iSN)=sqrt(f_ek*NbSN(iSN)*ESN/mloadSN(iSN)/ekBlast(iSN))
        endif
     else
        p_gas(iSN)=NbSN(iSN)*ESN/ekBlast(iSN)
     endif
  end do

  ! Loop over levels
  do ilevel=levelmin,nlevelmax
     ! Computing local volume (important for averaging hydro quantities) 
     dx=0.5D0**ilevel 
     dx_loc=dx*scale
     vol_loc=dx_loc**ndim
     ! Cells center position relative to grid center position
     do ind=1,twotondim  
        iz=(ind-1)/4
        iy=(ind-1-4*iz)/2
        ix=(ind-1-2*iy-4*iz)
        xc(ind,1)=(dble(ix)-0.5D0)*dx
        xc(ind,2)=(dble(iy)-0.5D0)*dx
        xc(ind,3)=(dble(iz)-0.5D0)*dx
     end do

     ! Loop over grids
     ncache=active(ilevel)%ngrid
!$omp parallel do private(igrid,ngrid,i,ind_grid,ind,iskip,ind_cell,&
!$omp             ok,x,y,z,iSN,dxx,dyy,dzz,dr_SN,d_gas,ielt,u,v,w)
     do igrid=1,ncache,nvector
        ngrid=MIN(nvector,ncache-igrid+1)
        do i=1,ngrid
           ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
        end do

        ! Loop over cells
        do ind=1,twotondim  
           iskip=ncoarse+(ind-1)*ngridmax
           do i=1,ngrid
              ind_cell(i)=iskip+ind_grid(i)
           end do

           ! Flag leaf cells
           do i=1,ngrid
              ok(i)=son(ind_cell(i))==0
           end do

           do i=1,ngrid
              if(ok(i))then
                 ! Get gas cell position
                 x=(xg(ind_grid(i),1)+xc(ind,1)-skip_loc(1))*scale
                 y=(xg(ind_grid(i),2)+xc(ind,2)-skip_loc(2))*scale
                 z=(xg(ind_grid(i),3)+xc(ind,3)-skip_loc(3))*scale
                 do iSN=1,nSN
                    ! Check if the cell lies within the SN radius
                    dxx=x-xSN(iSN,1)
                    dyy=y-xSN(iSN,2)
                    dzz=z-xSN(iSN,3)
                    dr_SN=dxx**2+dyy**2+dzz**2
                    if(dr_SN.lt.rmax2)then
                       d_gas=mloadSN(iSN)/vol_gas(iSN)
                       ! Compute the mass density in the cell
                       uold(ind_cell(i),1)=uold(ind_cell(i),1)+d_gas
                       ! Compute the metal density in the cell
                       if(metal)uold(ind_cell(i),imetal)=uold(ind_cell(i),imetal)+d_gas*ZloadSN(iSN)
                       do ielt=1,nelt
                          if (d_gas*celoadSN(iSN,ielt)>0d0) uold(ind_cell(i),ichem+ielt-1)=uold(ind_cell(i),ichem+ielt-1)+d_gas*celoadSN(iSN,ielt)
                       enddo
                      ! Velocity at a given dr_SN linearly interpolated between zero and uSedov
                       u=uSedov(iSN)*(dxx/rmax-dq(iSN,1))+vloadSN(iSN,1)
                       v=uSedov(iSN)*(dyy/rmax-dq(iSN,2))+vloadSN(iSN,2)
                       w=uSedov(iSN)*(dzz/rmax-dq(iSN,3))+vloadSN(iSN,3)
                       ! Add each momentum component of the blast wave to the gas
                       uold(ind_cell(i),2)=uold(ind_cell(i),2)+d_gas*u
                       uold(ind_cell(i),3)=uold(ind_cell(i),3)+d_gas*v
                       uold(ind_cell(i),4)=uold(ind_cell(i),4)+d_gas*w
                       ! Finally update the total energy of the gas
                       uold(ind_cell(i),5)=uold(ind_cell(i),5)+0.5*d_gas*(u*u+v*v+w*w)+p_gas(iSN)
                    endif
                 end do
              endif
           end do
           
        end do
        ! End loop over cells
     end do
     ! End loop over grids
  end do
  ! End loop over levels

  do iSN=1,nSN
     if(vol_gas(iSN)==0d0)then
        d_gas=mloadSN(iSN)/ekBlast(iSN)
        u=vloadSN(iSN,1)
        v=vloadSN(iSN,2)
        w=vloadSN(iSN,3)
        if(indSN(iSN)>0)then
           uold(indSN(iSN),1)=uold(indSN(iSN),1)+d_gas
           uold(indSN(iSN),2)=uold(indSN(iSN),2)+d_gas*u
           uold(indSN(iSN),3)=uold(indSN(iSN),3)+d_gas*v
           uold(indSN(iSN),4)=uold(indSN(iSN),4)+d_gas*w
           uold(indSN(iSN),5)=uold(indSN(iSN),5)+d_gas*0.5*(u*u+v*v+w*w)+p_gas(iSN)
           if(metal)uold(indSN(iSN),imetal)=uold(indSN(iSN),imetal)+d_gas*ZloadSN(iSN)
           do ielt=1,nelt
              if (d_gas*celoadSN(iSN,ielt)>0d0) uold(indSN(iSN),ichem+ielt-1)=uold(indSN(iSN),ichem+ielt-1)+d_gas*celoadSN(iSN,ielt)
           enddo

        endif
     endif
  end do

  if(verbose)write(*,*)'Exiting Sedov_blast'

end subroutine Sedov_blast
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine getSNonmyid(iSN_myid,nSN_myid,xSN,nSN)
  use amr_commons
  use pm_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer,dimension(1:nSN)::iSN_myid
  integer::nSN_myid,ii,info,nSN
  integer,dimension(1:8)::idom,jdom,kdom,cpu_min,cpu_max
  integer::lmin,iSN,nx_loc,ilevel,lmax,bit_length,maxdom,icpu
  integer::imin,jmin,kmin,imax,jmax,kmax,ndom,impi,i,j,k,ncpu_read
  integer,dimension(1:ncpu)::cpu_list
  logical,dimension(1:ncpu)::cpu_read
  real(dp)::scale,dx,dx_min,drSN
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  real(dp),dimension(1:3)::skip_loc
  real(dp)::xxmin,yymin,zzmin,xxmax,yymax,zzmax,dmax
  real(qdp),dimension(1:8)::bounding_min,bounding_max
  real(qdp)::dkey,order_min,oneqdp=1.0
  real(dp),dimension(1:nSN,1:3)::xSN
  integer, dimension(:), allocatable:: iflag

  allocate(iflag(1:nSN));
  iflag = 0

  ! Mesh spacing in that level
  nx_loc=(icoarse_max-icoarse_min+1)
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  skip_loc(1)=dble(icoarse_min)
  skip_loc(2)=dble(jcoarse_min)
  skip_loc(3)=dble(kcoarse_min)
  scale=boxlen/dble(nx_loc)
  dx_min=scale*0.5D0**nlevelmax

  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

  ! Maximum radius of the ejecta
  drSN=2d0*MAX(rcell*dx_min*scale_l/aexp,rbubble*3.08d18)
  drSN=drSN/scale_l

  !-----------------------
  ! Map parameters
  !-----------------------
  lmax=nlevelmax
  iSN_myid=0
  ii=0
!$omp parallel do private(iSN,cpu_read,xxmin,yymin,zzmin,xxmax,yymax,zzmax,dmax,ilevel,dx, &
!$omp        lmin,bit_length,maxdom,imin,imax,jmin,jmax,kmin,kmax, &
!$omp        dkey,ndom,idom,jdom,kdom,i,order_min,bounding_min,bounding_max,&
!$omp        cpu_min,cpu_max,impi,ncpu_read, cpu_list, j,k,icpu,ii) &
!$omp        schedule(guided)
!!$omp        schedule(dynamic,200)
  do iSN=1,nSN

     cpu_read=.false.
     ! Compute boundaries for the SN cube of influence
     xxmin=(xSN(iSN,1)-drSN)/scale ; xxmax=(xSN(iSN,1)+drSN)/scale
     yymin=(xSN(iSN,2)-drSN)/scale ; yymax=(xSN(iSN,2)+drSN)/scale
     zzmin=(xSN(iSN,3)-drSN)/scale ; zzmax=(xSN(iSN,3)+drSN)/scale

     if(TRIM(ordering).eq.'hilbert')then
        dmax=max(xxmax-xxmin,yymax-yymin,zzmax-zzmin)
        do ilevel=1,lmax
           dx=0.5d0**ilevel
           if(dx.lt.dmax)exit
        end do
        lmin=ilevel
        bit_length=lmin-1
        maxdom=2**bit_length
        imin=0; imax=0; jmin=0; jmax=0; kmin=0; kmax=0
        if(bit_length>0)then
           imin=int(xxmin*dble(maxdom))
           imax=imin+1
           jmin=int(yymin*dble(maxdom))
           jmax=jmin+1
           kmin=int(zzmin*dble(maxdom))
           kmax=kmin+1
        endif

        dkey=(real(2**(nlevelmax+1),kind=qdp)/real(maxdom,kind=qdp))**ndim
        ndom=1
        if(bit_length>0)ndom=8
        idom(1)=imin; idom(2)=imax
        idom(3)=imin; idom(4)=imax
        idom(5)=imin; idom(6)=imax
        idom(7)=imin; idom(8)=imax
        jdom(1)=jmin; jdom(2)=jmin
        jdom(3)=jmax; jdom(4)=jmax
        jdom(5)=jmin; jdom(6)=jmin
        jdom(7)=jmax; jdom(8)=jmax
        kdom(1)=kmin; kdom(2)=kmin
        kdom(3)=kmin; kdom(4)=kmin
        kdom(5)=kmax; kdom(6)=kmax
        kdom(7)=kmax; kdom(8)=kmax
        do i=1,ndom
           if(bit_length>0)then
              call hilbert3d(idom(i),jdom(i),kdom(i),order_min,bit_length,1)
           else
              order_min=0.0d0
           endif
           bounding_min(i)=(order_min)*dkey
           bounding_max(i)=(order_min+oneqdp)*dkey
        end do

        cpu_min=0; cpu_max=0
        do impi=1,ncpu
           do i=1,ndom
              if (   bound_key(impi-1).le.bounding_min(i).and.&
                   & bound_key(impi  ).gt.bounding_min(i))then
                 cpu_min(i)=impi
              endif
              if (   bound_key(impi-1).lt.bounding_max(i).and.&
                   & bound_key(impi  ).ge.bounding_max(i))then
                 cpu_max(i)=impi
              endif
           end do
        end do
        ncpu_read=0
        do i=1,ndom
           do j=cpu_min(i),cpu_max(i)
              if(.not. cpu_read(j))then
                 ncpu_read=ncpu_read+1
                 cpu_list(ncpu_read)=j
                 cpu_read(j)=.true.
              endif
           enddo
        enddo
     else
        ncpu_read=ncpu
        do j=1,ncpu
           cpu_list(j)=j
        end do
     end  if

     ! Create the index array for SN in processor myid                                  
     do k=1,ncpu_read
        icpu=cpu_list(k)
        if(icpu==myid)then
!          ii=ii+1
!          iSN_myid(ii)=iSN
           iflag(iSN) = 1
        endif
     enddo

  enddo
  ii = 0
  do iSN = 1, nSN
     if(iflag(iSN) .eq.1) then
	    ii = ii + 1
		iSN_myid(ii) = iSN
	 endif
  enddo

  ! Number of SN in processor myid                                                      
  nSN_myid=ii
  deallocate(iflag)

end subroutine getSNonmyid
!################################################################                       
!################################################################                       
!################################################################                       
!################################################################
