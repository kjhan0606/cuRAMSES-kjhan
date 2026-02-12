!=======================================================================
!=======================================================================
!=======================================================================
!=======================================================================
subroutine output_cone_hydro(obs)
  use amr_commons
  use hydro_commons
  use pm_commons
  use poisson_commons

  implicit none

#ifndef WITHOUTMPI
#include "mpif.h"
#endif
  
  integer::info,dummy_io,info2
  integer,parameter::tag=1118

  character(len=5) :: istep_str
  character(len=100) :: conedir, conecmd, conefile, infofile
  
  character(LEN=80)::fileloc
  character(LEN=5)::nchar
  real(kind=8) :: z1,z2,om0in,omLin,hubin,Lbox
  real(kind=8) :: Lobserver(3)
  integer::igrid,idim,icpu,ilevel
  integer::i


  !jaehyun1
  integer,dimension(:),allocatable::ind_grid
  integer(kind=8),dimension(:),allocatable::ii8,ii8out,ii8_out
  real(kind=8),dimension(:,:),allocatable::gpos,gvel,getc,tmp_gpos,tmp_gvel,tmp_getc
  real(kind=8),dimension(:,:),allocatable::gposout,gvelout,getcout
  real(kind=8),dimension(:,:),allocatable::gpos_out,gvel_out,getc_out
  real(kind=8),dimension(:),allocatable::gzout
  integer::nleaf1,ncache,ibound,iskip,ngout,j,istart,iglun,ind,nhvar,ivar,irad,k,end_tag,zoomed_level,print_mark,mncell,l,tngout,elongated_axis_cone,obs
  real(kind=8) :: cpi(3),dx,coord_distance,Omega0,OmegaL,OmegaR,coverH0,dist1,dist2,lboxz(3),minboxr(3),maxboxr(3),observer_cone(3),minboxr_cone(3),maxboxr_cone(3)
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  !jaehyun2

  logical::opened
  opened=.false.
  end_tag=-1
  minboxr=1
  maxboxr=0
  tngout=0

  if(nstep_coarse.lt.2) return

  z2=1./aexp_old-1.
  z1=1./aexp-1.
  if(z1<0.) z1=0. 
  if(z2.gt.zmax_cone) return
  if(abs(z2-z1)<1d-6) return
  
  nhvar=nvar-ndim+2
  om0in=omega_m
  omLin=omega_l
  hubin=h0/100.
  Lbox=boxlen_ini/hubin
 ! observer=(/Lbox/2.0,Lbox/2.0,Lbox/2.0/)
!  observer=(/Lbox/2.0,Lbox*0.0,Lbox/2.0/)
  if(obs==1) then
	elongated_axis_cone=elongated_axis_cone1
    do idim=1,3
      Lobserver(idim)=observer_cone1(idim)*Lbox
	  minboxr_cone(idim)=minboxr_cone1(idim)
	  maxboxr_cone(idim)=maxboxr_cone1(idim)
    enddo
  else
	elongated_axis_cone=elongated_axis_cone2
    do idim=1,3
      Lobserver(idim)=observer_cone2(idim)*Lbox
	  minboxr_cone(idim)=minboxr_cone2(idim)
	  maxboxr_cone(idim)=maxboxr_cone2(idim)
    enddo
  endif



  if(myid==1)write(*,*)'Computing and dumping lightcone (hydro of leaf cells)'
  call init_cosmo_cone(om0in,omLin,hubin,Omega0,OmegaL,OmegaR,coverH0)
  dist1=coord_distance(z1,Omega0,OmegaL,OmegaR,coverH0)
  dist2=coord_distance(z2,Omega0,OmegaL,OmegaR,coverH0)
  if(myid==1)  write(*,*)'lightcone redshifts',z1,z2
  if(myid==1)  write(*,*)'Distance (Mpc)',dist1,dist2



  iglun=3*ncpu+myid+10
  
  ! Determine the filename, dir, etc

  call title(nstep_coarse, istep_str)
  if(obs==1) then
    conedir = "cone_" // trim(istep_str) // "/observer1/"
  else
    conedir = "cone_" // trim(istep_str) // "/observer2/"
  endif
  conecmd = "mkdir -p " // trim(conedir)
  if(.not.withoutmkdir) then
     if (myid==1) call system(conecmd)
  endif
 
#ifndef WITHOUTMPI
  call MPI_BARRIER(MPI_COMM_WORLD, info)
#endif

  call title(myid,nchar)

  conefile = trim(conedir)//'cone_hydro_'//trim(istep_str)//'.out'
  fileloc=TRIM(conefile)//TRIM(nchar)



  ! Wait for the token
#ifndef WITHOUTMPI
  if(IOGROUPSIZECONE>0) then
     if (mod(myid-1,IOGROUPSIZECONE)/=0) then
        call MPI_RECV(dummy_io,1,MPI_INTEGER,myid-1-1,tag,&
             & MPI_COMM_WORLD,MPI_STATUS_IGNORE,info2)
     end if
  endif
#endif

   !jaehyun1

   mncell=65536
   allocate(gpos(1:ndim,1:mncell))
   allocate(gvel(1:ndim,1:mncell))
   allocate(getc(1:nhvar,1:mncell))
   allocate(ii8(1:mncell))



   do idim=1,ndim
	  maxboxr(idim)=maxboxr_cone(idim)*Lbox
	  minboxr(idim)=minboxr_cone(idim)*Lbox
      lboxz(idim)=maxboxr(idim)-minboxr(idim)
	  Lobserver(idim)=Lobserver(idim)-minboxr(idim)
   enddo


  do ibound=1,nboundary+ncpu
     do ilevel=levelmin_cone,nlevelmax
        dx=0.5**ilevel
	    nleaf1=0
        if(ibound<=ncpu) then
	      ncache=numbl(ibound,ilevel)
	      istart=headl(ibound,ilevel)
		else
	      ncache=numbb(ibound-ncpu,ilevel)
	      istart=headb(ibound-ncpu,ilevel)
		endif
		if(ncache>0) then
		  allocate(ind_grid(1:ncache))
		  igrid=istart
		  do i=1,ncache
		    ind_grid(i)=igrid
			igrid=next(igrid)
		  enddo
			
		  do ind=1,twotondim
		    iskip=ncoarse+(ind-1)*ngridmax
			do i=1,ncache
			  if(son(ind_grid(i)+iskip)==0) then
				nleaf1=nleaf1+1
			  endif
			enddo
		  enddo
		
		  if(nleaf1>0) then
		    j=0
			print_mark=0
		    do ind=1,twotondim
		      iskip=ncoarse+(ind-1)*ngridmax
              if(ind==1) cpi=(/0,0,0/)
              if(ind==2) cpi=(/1,0,0/)
              if(ind==3) cpi=(/0,1,0/)
              if(ind==4) cpi=(/1,1,0/)
              if(ind==5) cpi=(/0,0,1/)
              if(ind==6) cpi=(/1,0,1/)
              if(ind==7) cpi=(/0,1,1/)
              if(ind==8) cpi=(/1,1,1/)
		  	  do i=1,ncache
			    if(son(ind_grid(i)+iskip)==0) then
			      j=j+1
				  do idim=1,ndim
				    gpos(idim,j)=(xg(ind_grid(i),idim)+(cpi(idim)-0.5)*dx)*Lbox-minboxr(idim)
				    gvel(idim,j)=uold(ind_grid(i)+iskip,idim+1)/max(uold(ind_grid(i)+iskip,1),smallr)
				  enddo
				  getc(1,j)=dx
				  getc(2,j)=uold(ind_grid(i)+iskip,1)
                  k=2
#if NENER>0
				  do ivar=ndim+3,ndim+2+nener
				    getc(3+ivar-ndim-3,j)=(gamma_rad(ivar-ndim-2)-1d0)*uold(ind_grid(i)+iskip,ivar)
				  enddo
				  k=k+nener
#endif
				  getc(k+1,j)=uold(ind_grid(i)+iskip,ndim+2)
				  getc(k+1,j)=getc(k+1,j)-0.5d0*uold(ind_grid(i)+iskip,2)**2/max(uold(ind_grid(i)+iskip,1),smallr)
#if NDIM>1
                  getc(k+1,j)=getc(k+1,j)-0.5d0*uold(ind_grid(i)+iskip,3)**2/max(uold(ind_grid(i)+iskip,1),smallr)
#endif
#if NDIM>2
                  getc(k+1,j)=getc(k+1,j)-0.5d0*uold(ind_grid(i)+iskip,4)**2/max(uold(ind_grid(i)+iskip,1),smallr)
#endif
#if NENER>0
                  do irad=1,nener
                    getc(k+1,j)=getc(k+1,j)-uold(ind_grid(i)+iskip,ndim+2+irad)
                  end do
#endif
                  getc(k+1,j)=(gamma-1d0)*getc(k+1,j)
#if NVAR>NDIM+2+NENER
                  do ivar=ndim+3+nener,nvar
                    getc(k+2+ivar-ndim-3-nener,j)=uold(ind_grid(i)+iskip,ivar)/max(uold(ind_grid(i)+iskip,1),smallr)
			      enddo
#endif
				  getc(nhvar,j)=phi(ind_grid(i)+iskip)

				  if(j==mncell) then
                    print_mark=1
                  endif
                  
                  if(ind==twotondim .and. i==ncache .and. j<mncell) then
                    print_mark=1
                    allocate(tmp_gpos(1:ndim,1:j))
                    allocate(tmp_gvel(1:ndim,1:j))
                    allocate(tmp_getc(1:nhvar,1:j))
					do idim=1,ndim
                       tmp_gpos(idim,1:j)=gpos(idim,1:j)
                       tmp_gvel(idim,1:j)=gvel(idim,1:j)
					enddo
                    do k=1,nhvar
					   tmp_getc(k,1:j)=getc(k,1:j)
					enddo
                    deallocate(gpos)
                    deallocate(gvel)
                    deallocate(getc)
                    deallocate(ii8)
                    allocate(gpos(1:ndim,1:j))
                    allocate(gvel(1:ndim,1:j))
                    allocate(getc(1:nhvar,1:j))
                    allocate(ii8(1:j))
					do idim=1,ndim
                       gpos(idim,1:j)=tmp_gpos(idim,1:j)
                       gvel(idim,1:j)=tmp_gvel(idim,1:j)
					enddo
                    do k=1,nhvar
					   getc(k,1:j)=tmp_getc(k,1:j)
					enddo

                  endif
 
                  if (print_mark==1) then
                    !===========================================================================
                    ! Count number of leaf cells within a redshift range
                     call perform_my_selection(.true.,z1,z2, &
                     &                           om0in,omLin,hubin,lboxz, &
                     &                           Lobserver,elongated_axis_cone, &
        		     &	    					 j, nhvar, &
                     &                           ii8,gpos,gvel,getc,j, &
                     &                           ii8out,gposout,gvelout,getcout,gzout,ngout,.false.)

            
		             if (ngout > 0) then
       		            allocate(gposout(1:ndim,1:ngout))
		                allocate(gvelout(1:ndim,1:ngout))
		                allocate(getcout(1:nhvar,1:ngout))
		                allocate(ii8out(1:ngout))
		                allocate(gpos_out(1:ngout,1:ndim))
		                allocate(gvel_out(1:ngout,1:ndim))
		                allocate(getc_out(1:ngout,1:nhvar))
		                allocate(gzout(1:ngout))
		                allocate(ii8_out(1:ngout))
                        tngout=tngout+ngout

!            call extend_arrays_if_needed()
            ! Perform actual selection
                        call perform_my_selection(.false.,z1,z2, &
                        &                           om0in,omLin,hubin,lboxz, &
                        &                           Lobserver,elongated_axis_cone, &
		   	            &						    j, nhvar, &
                        &                           ii8,gpos,gvel,getc,j, &
                        &                           ii8out,gposout,gvelout,getcout,gzout,ngout,.false.)
          !===========================================================================

			            do idim=1,ndim
			               do k=1,ngout
				              gpos_out(k,idim)=gposout(idim,k)/Lbox
   				              gvel_out(k,idim)=gvelout(idim,k)
				              if(idim==1) then
				                  do l=1,nhvar
			                         getc_out(k,l)=getcout(l,k)
  					              enddo
	   			               endif
			                enddo
	  		             enddo

                         if(.not.opened) then
                            open(iglun,file=TRIM(fileloc),form='unformatted')
                            rewind(iglun)  
                            write(iglun)ncpu
                            write(iglun)ndim
				            write(iglun)nhvar
				            write(iglun)Lbox
                            opened=.true.
                         endif
	   	                 write(iglun)ngout
                         do idim=1,ndim
                   	        write(iglun)gpos_out(1:ngout,idim)
                         end do
                         do idim=1,ndim
                     	    write(iglun)gvel_out(1:ngout,idim)
                         end do
                         write(iglun)gzout(1:ngout)
                         do k=1,nhvar
                     	    write(iglun)getc_out(1:ngout,k)
                         end do
                       endif
            		   if(allocated(gposout)) deallocate(gposout)
		               if(allocated(gvelout)) deallocate(gvelout)
		               if(allocated(getcout)) deallocate(getcout)
		               if(allocated(ii8out)) deallocate(ii8out)
		               if(allocated(gpos_out)) deallocate(gpos_out)
		               if(allocated(gvel_out)) deallocate(gvel_out)
		               if(allocated(getc_out)) deallocate(getc_out)
		               if(allocated(ii8_out)) deallocate(ii8_out)
		               if(allocated(gzout)) deallocate(gzout)

                       print_mark=0
				       j=0

		   		     endif		  
				  endif
			   enddo
		    enddo

		    if(allocated(gpos)) deallocate(gpos)
		    if(allocated(gvel)) deallocate(gvel)
		    if(allocated(getc)) deallocate(getc)
		    if(allocated(ii8)) deallocate(ii8)
		    if(allocated(tmp_gpos)) deallocate(tmp_gpos)
		    if(allocated(tmp_gvel)) deallocate(tmp_gvel)
		    if(allocated(tmp_getc)) deallocate(tmp_getc)
					
  		    allocate(gpos(1:ndim,1:mncell))
		    allocate(gvel(1:ndim,1:mncell))
		    allocate(getc(1:nhvar,1:mncell))
		    allocate(ii8(1:mncell))

		  endif
		  deallocate(ind_grid)
		endif
     enddo
  enddo
  if(opened)write(iglun)end_tag
  if(opened)close(iglun)

  if (verbose)write(*,*)'cone hydro output=',myid,tngout
  if (tngout>0) then
     open(iglun,file=TRIM(fileloc)//".txt",form='formatted')
	 rewind(iglun)
	 write(iglun,*) ncpu
	 write(iglun,*) tngout
!     write(iglun,*) aexp_old
!     write(iglun,*) 1./(1.+z1)
     close(iglun)
   endif

   if (myid == 1) then
     call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
     infofile = trim(conedir)//'info_'//trim(istep_str)
     open(iglun,file=TRIM(infofile)//".txt",form='formatted')
	 rewind(iglun)
	 write(iglun,'("a           =",E23.15)') aexp_old
	 write(iglun,'("a+da        =",E23.15)') 1./(1.+z1)
	 write(iglun,'("unit_l      =",E23.15)') scale_l 
	 write(iglun,'("unit_t      =",E23.15)') scale_t
     write(iglun,'("unit_d      =",E23.15)') scale_d
     write(iglun,'("unit_v      =",E23.15)') scale_v
     write(iglun,'("unit_nH     =",E23.15)') scale_nH
     write(iglun,'("unit_T2     =",E23.15)') scale_T2
     close(iglun)
   endif



#ifndef WITHOUTMPI
   if(IOGROUPSIZECONE>0) then
      if(mod(myid,IOGROUPSIZECONE)/=0 .and.(myid.lt.ncpu))then
         dummy_io=1
         call MPI_SEND(dummy_io,1,MPI_INTEGER,myid-1+1,tag, &
              & MPI_COMM_WORLD,info2)
      end if
   endif
#endif
 
 
   if((opened.and.(tngout==0)).or.((.not.opened).and.(tngout>0))) then
      write(*,*)'Error in output_cone_cell'
      write(*,*)'tngout=',tngout,'opened=',opened
      stop
   endif


!jaehyun2


end subroutine output_cone_hydro

!jhshin1
subroutine perform_my_selection(justcount,z1,z2, &
     &                          om0in,omLin,hubin,lboxz, &
     &                          observer,elongated_axis_cone,&
	 &						    nvector, nvar,&
     &                          pid,pos,vel,etc,npart, &
     &                          pidout,posout,velout,etcout,zout,npartout,verbose)
!jhshin2
  !===========================================================================
  ! All the quantities below are real*8 except
  !      juscount : logical
  !      npart,npartout : integer*4
  !
  ! juscount : .true. to just count the particles to be selected. Needed
  !            to allocate appropriately the arrays posout,velout and zout
  !            The parameter npartout becomes then an output given the number
  !            of particles selected.
  !            .false. to perform the actual selection: npartout is an input
  !            then  posout, velout, zout are appropriate outputs.
  !
  ! z1,z2    : the lightcone part of interest is in between z1 and z2, 
  !            with z1 < z2. If we consider a redshift z(t) where all the 
  !            particles are synchrone, and if coarse timestep is
  !            a fixed dt, it is most appropriate to choose z1 and z2 such
  !            that z1=z(t+dt/2) and z2=z(t-dt/2) to have best accuracy.
  !
  ! om0in    : the value of the cosmological parameter omega0 (typically 0.3)
  !
  ! omLin    : the value of the cosmological constant lambda (typically 0.7)
  !
  ! hubin    : the value of H0/100, where H0 is the present Hubble constant
  !            in km/s/Mpc (typically 0.7)
  !
  ! Lbox     : the comoving size of the simulation box in Mpc (NOT in Mpc/h)
  !
  ! observer(3) : the observer position in the box in Mpc, assuming that
  !            coordinates are in [0,Lbox[
  !
  ! thetay   : half the opening angle in degrees of the lightcone along y direction
  !            (it should be obviously smaller than 90 degrees to avoid catastrophic 
  !            behavior). The lightcone is assume to be aligned with x axis (after
  !            appropriates rotations given by angles theta and phi)
  !
  ! thetaz   : half the opening angle in degrees of the lightcone along z direction
  !            Given thetay and thetaz, the area of the survey is thus 4.thetay.thetaz
  !
  ! theta, phi : 2 angles in degrees defining a rotation to avoid alignement of
  !            the lightcone with the major axes of the simulation box.
  !            Example : theta=21, phi=17.
  !
  ! pos(3,npart) : comoving positions of the input particles in Mpc, assumed to be
  !            in [0,Lbox[.
  !
  ! vel(3,npart) : velocities of the input particles (in any unit, it does not
  !            matter)
  !
  ! npart    : number of input particles to be treated
  !
  ! posout(3,npartout) : output comoving positions of selected particles in Mpc.
  !
  ! velout(3,npartout) : output velocities of selected particles
  !
  ! zout(npartout) : output redshift of selected particles
  !
  ! npartout : number of selected particles. To be computed appropriately,
  !            this routine must be called with juscount=.true., which will give
  !            npartout as an output. Then this routine must be called with 
  !            juscount=.false. with the correct value of npartout, after having
  !            allocated correctly arrays posout,velout,zout.
  !===========================================================================
  !use amr_parameters, ONLY: nvector
  implicit none
  logical :: justcount,verbose
  integer :: npart,npartout
  real(kind=8) :: z1,z2,om0in,omLin,hubin
  real(kind=8) :: Omega0,OmegaL,OmegaR,coverH0
  real(kind=8) :: observer(3)
  !jaehyun1
  integer(kind=8) :: pid(1:nvector)
  integer(kind=8) :: pidout(1:npartout)
  integer :: nvector,nvar,elongated_axis_cone,irange_min(3),irange_max(3)
  real(kind=8) :: lboxz(3),coord(3)
  !jaehyun2
  !jhshin1
  real(kind=8) :: pos(1:3,1:nvector),vel(1:3,1:nvector),etc(1:nvar,1:nvector)
  real(kind=8) :: posout(3,npartout),velout(3,npartout),zout(npartout)
  real(kind=8) :: etcout(1:nvar,npartout)
  !jhshin2
  real(kind=8) :: coord_distance
  real(kind=8) :: dist1,dist2
  real(kind=8) :: dist,cdist,dxtest1,dxtest2,facnorm
  real(kind=8) :: small=1d-6
  integer,dimension(:,:,:),allocatable::ulbox
  
  integer :: i,j,k,l,np,npartcount
  
  if (verbose) write(*,*) 'Entering perform_my_selection'
  
  
  ! Initialize cosmological parameters
  call init_cosmo_cone(om0in,omLin,hubin,Omega0,OmegaL,OmegaR,coverH0)
  
  
  
  ! Compute comoving distance of the photon planes from the observer
  ! dist1,dist2=integral of c.dt/a between zero and z1,z2
  dist1=coord_distance(z1,Omega0,OmegaL,OmegaR,coverH0)
  dist2=coord_distance(z2,Omega0,OmegaL,OmegaR,coverH0)
  
  
  facnorm=1.0d0/(dist2-dist1)
 
  call compute_replica_box_index(dist2,lboxz,observer,elongated_axis_cone,&
       &                       irange_min,irange_max)    

  allocate(ulbox(irange_min(1):irange_max(1),irange_min(2):irange_max(2),irange_min(3):irange_max(3)))
  ulbox=0
  call compute_replica_box(dist1,dist2,lboxz,observer,irange_min,irange_max,&
       &                       ulbox)    
  npartcount=0    
  ! loop on all the replica of potential interest
  do k=irange_min(3),irange_max(3),1
     do j=irange_min(2),irange_max(2),1
        do i=irange_min(1),irange_max(1),1
           if(ulbox(i,j,k)==1) then
             do np=1,npart
				if((pos(1,np)>=0.) .and. (pos(1,np)<=lboxz(1)) .and. (pos(2,np)>=0.) .and. (pos(2,np)<=lboxz(2)) .and. (pos(3,np)>=0.) .and. (pos(3,np)<=lboxz(3))) then
                   coord(1)=pos(1,np)+lboxz(1)*dble(i)-observer(1)
                   coord(2)=pos(2,np)+lboxz(2)*dble(j)-observer(2)
                   coord(3)=pos(3,np)+lboxz(3)*dble(k)-observer(3)
              
!                if (xcoord > small) then ! To avoid divergences near the origin
                   dist=sqrt(coord(1)**2+coord(2)**2+coord(3)**2)


				   cdist=abs(coord(abs(elongated_axis_cone)))
				   if(cdist > small) then
				      if (cdist > dist1 .and. cdist <= dist2) then
                        npartcount=npartcount+1

                        if (.not. justcount) then
							do l=1,3
                              posout(l,npartcount)=coord(l)
    	   			  	      velout(l,npartcount)=vel(l,np)
                            enddo
                        !jhshin1
						    do l=1,nvar
                              etcout(l,npartcount)=etc(l,np)
						    enddo
							pidout(npartcount)=pid(np)
                        !jhshin2

                        ! Compute the redshift of the particle using linear
                        ! interpolation
                            dxtest1=dist-dist1
                            dxtest2=dist2-dist
                            zout(npartcount)=(dxtest1*z2+dxtest2*z1)*facnorm
                        endif
                     endif
                  endif
               endif
             enddo
           endif
        enddo
     enddo
  enddo
  npartout=npartcount
  if (verbose) write(*,*) 'End of perform_my_selection'
  if (allocated(ulbox)) deallocate(ulbox)
end subroutine perform_my_selection









!jaehyun1
subroutine compute_replica_box_index(dist2,lboxz,observer, elongated_axis_cone,&
     &                           irange_min,irange_max)
  !===========================================================================
  ! 2*theta1 and 2*theta2 are the opening angles of the lightcone in degrees.
  ! The observer position is expressed in comoving Mpc, as well as the simulation
  ! box size Lbox. Furthermore, the positions of particles inside the simulation
  ! are supposed to be in [0,Lbox[.
  ! z1 and z2 are the redshifts of the successive photon planes, z1 < z2
  !===========================================================================
  implicit none
  real(kind=8) :: observer(3),lboxz(3),dist2,dx
  integer ::myint,elongated_axis_cone,irange_min(3),irange_max(3),i

  if(elongated_axis_cone .lt. 0) then
    do i=1,3
      irange_max(i)=0
      if(i==abs(elongated_axis_cone)) then
        dx=observer(i)-dist2
  	    irange_min(i)=myint(dx/lboxz(i))
      else
        irange_min(i)=0
      endif
    enddo

  else
    do i=1,3
      irange_min(i)=0
      if(i==elongated_axis_cone) then
        dx=observer(i)+dist2
	    irange_max(i)=myint(dx/lboxz(i))
      else
        irange_max(i)=0
      endif
    enddo
  endif
!  dx=observer(2)+dist
!  nrepyp=myint(dx/lboxz(2))
!  dx=observer(3)-dist2
!  nrepzm=myint(dx/lboxz(3)
end subroutine compute_replica_box_index
  


subroutine compute_replica_box(dist1,dist2,lboxz,observer,irange_min,irange_max,&
     &                           ulbox)
  !===========================================================================
  ! 2*theta1 and 2*theta2 are the opening angles of the lightcone in degrees.
  ! The observer position is expressed in comoving Mpc, as well as the simulation
  ! box size Lbox. Furthermore, the positions of particles inside the simulation
  ! are supposed to be in [0,Lbox[.
  ! z1 and z2 are the redshifts of the successive photon planes, z1 < z2
  !===========================================================================
  implicit none
  real(kind=8) :: observer(3),lboxz(3),dist1,dist2,i00,i01,i1,j00,j01,j1,k00,k01,k1,dx,dy,dz,ndist,fdist
  integer :: irange_min(3),irange_max(3),i,j,k
  integer :: ulbox(irange_min(1):irange_max(1),irange_min(2):irange_max(2),irange_min(3):irange_max(3))


  do i=irange_min(1),irange_max(1),1
    do j=irange_min(2),irange_max(2),1
      do k=irange_min(3),irange_max(3),1
        if(i==0 .and. j==0 .and. k==0) then
          i00=real(int(observer(1)/0.5/lboxz(1)))
          j00=real(int(observer(2)/0.5/lboxz(2)))
          k00=real(int(observer(3)/0.5/lboxz(3)))
          i1=1.
		  j1=1.
		  k1=1.

		  if(i00 .ge. 1.0) then
            i1=0.
		  endif
		  if(j00 .ge. 1.0) then
            j1=0.
		  endif
		  if(k00 .ge. 1.0) then
            k1=0.
		  endif
          
		  fdist=sqrt((lboxz(1)*i1-observer(1))**2.+(lboxz(2)*j1-observer(2))**2.+(lboxz(3)*k1-observer(3))**2.)
          if(fdist > dist1) then
            ulbox(i,j,k)=1
		  endif
		else
		  i00=0.
		  i01=dble(abs(i))
		  i1=0.
          if(i>0) then
          	i00=dble(abs(i)-1)
            i1=1.
			dx=lboxz(1)-observer(1)
		  elseif(i<0) then
          	i00=dble(abs(i)-1)
            i1=1.
			dx=observer(1)
		  endif

		  j00=0.
		  j01=dble(abs(j))
		  j1=0.
          if(j>0) then
          	j00=dble(abs(j)-1)
            j1=1.
			dy=lboxz(2)-observer(2)
		  elseif(j<0) then
          	j00=dble(abs(j)-1)
            j1=1.
			dy=observer(2)
		  endif

		  k00=0.
		  k01=dble(abs(k))
		  k1=0.
          if(k>0) then
          	k00=dble(abs(k)-1)
            k1=1.
			dz=lboxz(3)-observer(3)
		  elseif(k<0) then
          	k00=dble(abs(k)-1)
            k1=1.
			dz=observer(3)
		  endif

          ndist=sqrt((dx*i1+lboxz(1)*i00)**2.+(dy*j1+lboxz(2)*j00)**2.+(dz*k1+lboxz(3)*k00)**2.)
          fdist=sqrt((dx*i1+lboxz(1)*i01)**2.+(dy*j1+lboxz(2)*j01)**2.+(dz*k1+lboxz(3)*k01)**2.)
          if(ndist<dist2 .and. fdist>dist1) then
            ulbox(i,j,k)=1
          endif
		endif


      enddo
	enddo
  enddo
  
end subroutine compute_replica_box

!jaehyun2



 
!===================
!cone cosmo routines
!===================


!===========================================================================
subroutine init_cosmo_cone(om0in,omLin,hubin,Omega0,OmegaL,OmegaR,coverH0)
  !===========================================================================
  ! om0in : the value of omega0
  ! omLin : the value of Lambda
  !         We MUST have omega0+Lambda=1.0d0
  ! hubin : the value of H0/100 where H0 is the present Hubble constant 
  !         in km/s/Mpc
  !===========================================================================
  implicit none
  real(kind=8) :: om0in,omLin,hubin
  real(kind=8) :: Omega0,OmegaL,OmegaR,coverH0
  real(kind=8) :: verysmall=1d-6
  omega0=om0in
  omegaL=omLin
  omegaR=1.0d0-omega0-omegaL
  if (abs(omegaR) > verysmall) then
     write(*,*) 'ERROR in propagate_photons, init_cosmo.'
     write(*,*) 'This routine works only for flat universes, omega0+Lambda=1.'
     STOP
  endif
  coverH0=299792.5d0/(100.0d0*hubin)
end subroutine init_cosmo_cone


!===========================================================================
function coord_distance(zz,Omega0,OmegaL,OmegaR,coverH0)
  !===========================================================================
  implicit none
  real(kind=8) :: z,res,coord_distance,zz
  real(kind=8) :: Omega0,OmegaL,OmegaR,coverH0
  z=abs(zz)
  call qromb(0.d0,z,res,omega0,omegaL,OmegaR)
  coord_distance=coverH0*res
  if (zz.lt.0) coord_distance=-coord_distance
end function coord_distance

!===========================================================================
function funcE(z,Omega0,OmegaL,OmegaR)
  !===========================================================================
  implicit none
  real(kind=8) :: funcE,z,HsurH0
  real(kind=8) :: omega0,omegaL,OmegaR
  
  funcE=1.d0/HsurH0(z,Omega0,OmegaL,OmegaR)
end function funcE

!===========================================================================
function HsurH0(z,omega0,omegaL,OmegaR)
  !===========================================================================
  implicit none
  real(kind=8) :: z,omega0,omegaL,OmegaR,HsurH0
  HsurH0=sqrt(Omega0*(1.d0+z)**3+OmegaR*(1.d0+z)**2+OmegaL)
end function HsurH0


!===========================================================================
SUBROUTINE qromb(a,b,ss,omega0,omegaL,OmegaR)
  !===========================================================================
  implicit none
  INTEGER :: JMAX,JMAXP,K,KM
  REAL(kind=8) :: a,b,ss,EPS,omega0,omegaL,OmegaR
  PARAMETER (EPS=1.d-6, JMAX=20, JMAXP=JMAX+1, K=5, KM=K-1)
  !  USES polint,trapzd
  INTEGER :: j
  REAL(kind=8) :: dss,h(JMAXP),s(JMAXP)
  h(1)=1.
  do j=1,JMAX
     call trapzd(a,b,s(j),j,omega0,omegaL,OmegaR)
     if (j.ge.K) then
        call polint(h(j-KM),s(j-KM),K,0.d0,ss,dss)
        if (abs(dss).le.EPS*abs(ss)) return
     endif
     s(j+1)=s(j)
         h(j+1)=0.25*h(j)
      enddo

  print *, 'too many steps in qromb'
END SUBROUTINE qromb

!===========================================================================
SUBROUTINE polint(xa,ya,n,x,y,dy)
  !===========================================================================
  implicit none
  INTEGER :: n,NMAX
  REAL(kind=8) :: dy,x,y,xa(n),ya(n)
  PARAMETER (NMAX=10)
  INTEGER :: i,m,ns
  REAL(kind=8) ::den,dif,dift,ho,hp,w,c(NMAX),d(NMAX)
  ns=1
  dif=abs(x-xa(1))
  do i=1,n
     dift=abs(x-xa(i))
     if (dift.lt.dif) then
        ns=i
        dif=dift
     endif
     c(i)=ya(i)
     d(i)=ya(i)
  enddo
  y=ya(ns)
  ns=ns-1
  do m=1,n-1
     do i=1,n-m
        ho=xa(i)-x
        hp=xa(i+m)-x
        w=c(i+1)-d(i)
        den=ho-hp
        if(den.eq.0.) print *, 'failure in polint'
        den=w/den
        d(i)=hp*den
        c(i)=ho*den
     enddo
     if (2*ns.lt.n-m)then
        dy=c(ns+1)
     else
        dy=d(ns)
        ns=ns-1
     endif
     y=y+dy
  enddo
  return
END SUBROUTINE polint

!===========================================================================
SUBROUTINE trapzd(a,b,s,n,omega0,omegaL,OmegaR)
  !===========================================================================
  implicit none
  INTEGER :: n
  REAL(kind=8) :: a,b,s,funcE,omega0,omegaL,OmegaR
  INTEGER :: it,j
  REAL(kind=8) :: del,sum,tnm,x
  if (n.eq.1) then
     s=0.5*(b-a)*(funcE(a,omega0,omegaL,OmegaR)+funcE(b,omega0,omegaL,OmegaR))
  else
     it=2**(n-2)
     tnm=it
     del=(b-a)/tnm
     x=a+0.5*del
     sum=0.
     do j=1,it
        sum=sum+funcE(x,omega0,omegaL,OmegaR)
        x=x+del
     enddo
     s=0.5*(s+(b-a)*sum/tnm)
  endif
  return
END SUBROUTINE trapzd
!=======================================================================
function myint(x)
  !=======================================================================
  ! The REAL int function
  !=======================================================================
  real(kind=8) :: x
  integer :: myint
  
  if (x >= 0.0d0) then
     myint=int(x)
  else
     myint=int(x)-1
  endif
end function myint
