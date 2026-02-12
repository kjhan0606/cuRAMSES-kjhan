!=======================================================================
!=======================================================================
!=======================================================================
!=======================================================================
subroutine output_sphere_hydro()
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
  character(len=100) :: spheredir, spherecmd, spherefile, infofile
  
  character(LEN=80)::fileloc
  character(LEN=5)::nchar
  integer::igrid,idim,icpu,ilevel
  integer::i,j,k,l,m


  !jaehyun1
  integer,dimension(:),allocatable::ind_grid, stag, ttag
  integer(kind=8),dimension(:),allocatable::ii8,ii8out,ii8_out
  real(kind=8),dimension(:,:),allocatable::gpos,gvel,getc,tmp_gpos,tmp_gvel,tmp_getc
  real(kind=8),dimension(:,:),allocatable::gposout,gvelout,getcout
  real(kind=8),dimension(:,:),allocatable::gpos_out,gvel_out,getc_out
  integer::nleaf1,ncache,ibound,iskip,ngout,istart,iglun,ind,nhvar,ivar,irad,end_tag,print_mark,mncell,tngout,nprint,ngout0,tnprint
  real(kind=8) :: cpi(3),dx,ccenter(3),cradius
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  !jaehyun2

  logical::opened
  opened=.false.
  end_tag=-1
  tngout=0

  
  nhvar=nvar-ndim+2



  if(myid==1) write(*,*)'Dumping data in a spherical region (leaf cells)'



  iglun=3*ncpu+myid+10
  
  ! Determine the filename, dir, etc

  call title(nstep_coarse, istep_str)
  spheredir = "sphere/sphere_" // trim(istep_str)
  spherecmd = "mkdir -p " // trim(spheredir)
  if(.not.withoutmkdir) then
     if (myid==1) call system(spherecmd)
  endif
 
#ifndef WITHOUTMPI
  call MPI_BARRIER(MPI_COMM_WORLD, info)
#endif

  call title(myid,nchar)

  spherefile = trim(spheredir)//'sphere_hydro_'//trim(istep_str)//'.out'
  fileloc=TRIM(spherefile)//TRIM(nchar)


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
   allocate(gpos_out(1:mncell,1:ndim))
   allocate(gvel_out(1:mncell,1:ndim))
   allocate(getc_out(1:mncell,1:nhvar))
   allocate(ii8_out(1:mncell))
   allocate(stag(1:mncell))
   allocate(ttag(1:mncell))
!   allocate(gposout(1:ndim,1:mncell))
!   allocate(gvelout(1:ndim,1:mncell))
!   allocate(getcout(1:nhvar,1:mncell))
!   allocate(ii8out(1:mncell))


 
  stag=0
  ttag=0


  j=0
  nprint=0
  do ibound=1,nboundary+ncpu
     if(ibound == myid) then
        do ilevel=levelmin_sphere,nlevelmax
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
				             gpos(idim,j)=xg(ind_grid(i),idim)+(cpi(idim)-0.5)*dx
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
                  
                          if(ind==twotondim .and. i==ncache) then
                             print_mark=1
                          endif
 
                          if(print_mark==1) then
                          !===========================================================================
                          ! Count number of leaf cells within a spherical region
                             call selecting_cells_in_sphere(sradius, scenter, mncell, j, gpos, ngout, stag, .false.)
  !      		             &                      mncell, nhvar, &
  !                           &                      ii8,gpos,gvel,getc, j&
  !                           &                      ii8out,gposout,gvelout,getcout,ngout,.false.)


		                     if(ngout > 0) then
                                tngout=tngout+ngout

            			  ! Perform actual selection
 !                               call selecting_cells_in_sphere(.false.,sradius, scenter,&
 !		   	                 &                         mncell, nhvar, &
 !                            &                         ii8,gpos,gvel,getc, j&
 !                            &                         ii8out,gposout,gvelout,getcout,ngout,.false.)
                         !===========================================================================
				                l=0
                                do k=1,mncell
								   if(stag[k] .eq. 1) then
								      l=l+1
								      ttag[l]=k									
								   endif
								enddo

  						        tnprint=nprint+ngout
						        ngout0=ngout
						        if(tnprint .ge. mncell) then
                                   ngout0=ngout-tnprint+mncell
						           print_mark=3
						        endif

						        do k=1,ngout0
						           nprint=nprint+1
  								   l=ttag[k]
						           do idim=1,ndim
						              gpos_out(nprint,idim)=gpos(l,idim)
						              gvel_out(nprint,idim)=gvel(l,idim)	  
                                   enddo
						           do m=1,nhvar
						             getc_out(nprint,m)=getc(l,m)
                                   enddo
						        enddo
							 endif
							 j=0
					      endif
						

					      if(print_mark==3) then

                             if(.not.opened) then
                                open(iglun,file=TRIM(fileloc),form='unformatted')
                                rewind(iglun)  
                                write(iglun)ncpu
                                write(iglun)ndim
				                write(iglun)nhvar
                                opened=.true.
                             endif
	   	                     write(iglun)nprint
                             do idim=1,ndim
                   	            write(iglun)gpos_out(1:nprint,idim)
                             enddo
                             do idim=1,ndim
                     	        write(iglun)gvel_out(1:nprint,idim)
                             enddo
                             do k=1,nhvar
                                write(iglun)getc_out(1:nprint,k)
                             enddo
						   
                             nprint=0
                             if(ngout0 .lt. ngout) then
						  
                                do idim=1,ndim
                                   do k=ngout0+1,ngout
                                      nprint=nprint+1
									  l=ttag[k]
       					              gpos_out(nprint,idim)=gpos(l,idim)
							          gvel_out(nrpint,idim)=gvel(l,idim)
					                  if(idim==1) then
					                     do m=1,nhvar
			                                getc_out(nprint,m)=getc(l,m)
  					                     enddo
	   			                      endif
							       enddo
                                enddo
                             endif
                          endif
                          print_mark=0
			           endif
		            enddo
		         enddo
		      endif
		      deallocate(ind_grid)
           endif 
	    enddo
        if(nprint .gt. 0) then
           if(.not.opened) then
               open(iglun,file=TRIM(fileloc),form='unformatted')
               rewind(iglun)  
               write(iglun)ncpu
               write(iglun)ndim
               write(iglun)nhvar
               opened=.true.
            endif
            write(iglun)nprint
            do idim=1,ndim
                write(iglun)gpos_out(1:nprint,idim)
            enddo
            do idim=1,ndim
               write(iglun)gvel_out(1:nprint,idim)
            enddo
            do k=1,nhvar
               write(iglun)getc_out(1:nprint,k)
            enddo	 
            nprint=0
        endif
     endif
  enddo
  if(opened)write(iglun)end_tag
  if(opened)close(iglun)


  deallocate(gpos)
  deallocate(gvel)
  deallocate(getc)
  deallocate(ii8)
  deallocate(gpos_out)
  deallocate(gvel_out)
  deallocate(getc_out)
  deallocate(ii8_out)
  deallocate(stag)
!  deallocate(gposout)
!  deallocate(gvelout)
!  deallocate(getcout)
!  deallocate(ii8out)



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
     infofile = trim(spheredir)//'info_'//trim(istep_str)
     open(iglun,file=TRIM(infofile)//".txt",form='formatted')
     rewind(iglun)
     write(iglun,'("a           =",E23.15)') aexp_old
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
      write(*,*)'Error in output_sphere_cell'
      write(*,*)'tngout=',tngout,'opened=',opened
      stop
   endif


!jaehyun2


end subroutine output_sphere_hydro

!jhshin1
subroutine selecting_cells_in_sphere(r, center, nvector, npart, pos, npartout, stag, verbose)
!     &                          center,&
!	 &						    nvector, nvar,&
!     &                          pid,pos,vel,etc,npart, npart, &
!     &                          pidout,posout,velout,etcout,npartout,verbose)
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
  !
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
  !
  ! npartout : number of selected particles. To be computed appropriately,
  !            this routine must be called with juscount=.true., which will give
  !            npartout as an output. Then this routine must be called with 
  !            juscount=.false. with the correct value of npartout, after having
  !            allocated correctly arrays posout,velout,zout.
  !===========================================================================
  !use amr_parameters, ONLY: nvector
  implicit none
 ! logical :: justcount,verbose
  logical :: verbose
  integer :: npart,npartout,stag(1:nvector)
  real(kind=8) :: center(3),r
  !jaehyun1
!  integer(kind=8) :: pid(1:nvector)
 ! integer(kind=8) :: pidout(1:nvector)
  integer :: nvector!,nvar
  real(kind=8) :: dx,dy,dz
  !jaehyun2
  !jhshin1
  real(kind=8) :: pos(1:3,1:nvector)!,vel(1:3,1:nvector),etc(1:nvar,1:nvector)
!  real(kind=8) :: posout(3,nvector),velout(3,nvector)
 ! real(kind=8) :: etcout(1:nvar,nvector)
  !jhshin2
  real(kind=8) :: dist
  
  integer :: i,np,npartcount
  
  if (verbose) write(*,*) 'Entering perform_my_selection'
  
  
  stag=0
  
   
  npartcount=0    
  do np=1,npart
     dx=pos(1,np)-center(1)
     dy=pos(2,np)-center(2)
     dz=pos(3,np)-center(3)
              
     dist=sqrt(dx**2+dy**2+dz**2)
   

	 if(dist .le. r) then
        npartcount=npartcount+1
        stag(np)=1
     endif
  enddo
  !      if (.not. justcount) then
!	       do i=1,3
 !             posout(i,npartcount)=pos(i,np)
  !     	      velout(i,npartcount)=vel(i,np)
   !        enddo
	!	   do i=1,nvar
    !          etcout(i,npartcount)=etc(i,np)
!		   enddo
!		   pidout(npartcount)=pid(np)
 !       endif
 !    endif
 ! enddo
  npartout=npartcount
  if (verbose) write(*,*) 'End of perform_my_selection'
end subroutine selecting_cells_in_sphere
