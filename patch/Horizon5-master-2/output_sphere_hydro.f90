subroutine output_sphere_hydro(isphere)
  use amr_commons
  use hydro_commons
  use poisson_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'  
#endif


  integer::i,ivar,ncache,ind,ilevel,igrid,iskip,ilun,istart,ibound,irad,nprint,j,nleaf,cpi(8,3),idim,k,isphere
  integer,allocatable,dimension(:)::ind_grid,ttag
  integer(kind=2),allocatable,dimension(:)::stag
  real(dp),allocatable,dimension(:)::xdp
  real(kind=4),allocatable,dimension(:)::fxdp
  character(LEN=2)::sisphere
  character(LEN=5)::nchar
  character(LEN=5)::istep_str
  character(LEN=150)::fileloc
  character(len=150) :: spheredir, spherecmd, spherefile, infofile
  integer,parameter::tag=1121
  integer::info,dummy_io,info2,endtag
  real(dp)::dx,sradius,scenter(3)
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  logical::opened
  opened=.false.
  endtag=0
  ilun=ncpu+myid+10
     

  ilun=3*ncpu+myid+10
  
  if(isphere .eq. 1) then
    sradius=sradius1
    scenter=scenter1
    sisphere="01"
    if(myid==1)write(*,*)'Dumping a spherical region 1'
  else if(isphere .eq. 2) then
    sradius=sradius2
    scenter=scenter2
    sisphere="02"
    if(myid==1)write(*,*)'Dumping a spherical region 2'
  else if(isphere .eq. 3) then
    sradius=sradius3
    scenter=scenter3
    sisphere="03"
    if(myid==1)write(*,*)'Dumping a spherical region 3'
  else if(isphere .eq. 4) then
    sradius=sradius4
    scenter=scenter4
    sisphere="04"
    if(myid==1)write(*,*)'Dumping a spherical region 4'
  else if(isphere .eq. 5) then
    sradius=sradius5
    scenter=scenter5
    sisphere="05"
    if(myid==1)write(*,*)'Dumping a spherical region 5'
  endif



  ! Determine the filename, dir, etc

  call title(nstep_coarse, istep_str)

  spheredir = "sphere/sphere_" // trim(sisphere) // "/nstep_" // trim(istep_str) 

  spherecmd = "mkdir -p " // trim(spheredir)
  if(.not.withoutmkdir) then
!     if (myid==1) call system(spherecmd)
     if (myid==1) call execute_command_line(spherecmd,wait=.true.)
  endif

#ifndef WITHOUTMPI
  call MPI_BARRIER(MPI_COMM_WORLD, info)
#endif
  spherefile = trim(spheredir)//'/sphere_hydro_'//trim(istep_str)//'.out'
  call title(myid,nchar)
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
 
  cpi(1,1:3)=(/0,0,0/)
  cpi(2,1:3)=(/1,0,0/)
  cpi(3,1:3)=(/0,1,0/)
  cpi(4,1:3)=(/1,1,0/)
  cpi(5,1:3)=(/0,0,1/)
  cpi(6,1:3)=(/1,0,1/)
  cpi(7,1:3)=(/0,1,1/)
  cpi(8,1:3)=(/1,1,1/)



  do ilevel=1,nlevelmax
     do ibound=1,nboundary+ncpu
	    if(ibound==myid) then
		   dx=0.5d0**ilevel
		   nleaf=0
           if(ibound<=ncpu)then 
              ncache=numbl(ibound,ilevel)
              istart=headl(ibound,ilevel)
           else
              ncache=numbb(ibound-ncpu,ilevel)
              istart=headb(ibound-ncpu,ilevel)
           end if
           if(ncache>0)then
              allocate(ind_grid(1:ncache))
              ! Loop over level grids
              igrid=istart
              do i=1,ncache
                 ind_grid(i)=igrid
                 igrid=next(igrid)
              end do

              do ind=1,twotondim
                 iskip=ncoarse+(ind-1)*ngridmax
                 do i=1,ncache
                    if(son(ind_grid(i)+iskip)==0)then
				       nleaf=nleaf+1
					endif
				 enddo
			  enddo
              if(nleaf .gt. 0) then
                 allocate(xdp(1:nleaf))
                 allocate(stag(1:nleaf))
				 xdp=0.
				 stag=0
				 nprint=0
                 j=0
                 do ind=1, twotondim
	                iskip=ncoarse+(ind-1)*ngridmax
                    do i=1,ncache
                       if(son(ind_grid(i)+iskip)==0)then
			              j=j+1
					      do idim=1,ndim
					         xdp(j)=xdp(j)+(xg(ind_grid(i),idim)+(cpi(ind,idim)-0.5)*dx-scenter(idim))**2
					      enddo
						  xdp(j)=dsqrt(xdp(j))
						  if(xdp(j) .le. sradius) then
						     stag(j)=1
						     nprint=nprint+1
						  endif
					   endif
				    enddo
				 enddo
                 
				 deallocate(xdp)
				 if(nprint .gt. 0) then
				    allocate(ttag(1:nprint))
				    allocate(xdp(1:nprint))
				    allocate(fxdp(1:nprint))
					if(.not.opened) then
                       open(unit=ilun,file=fileloc,form='unformatted')
					   rewind(ilun)
                       write(ilun)ndim
                       write(ilun)nvar+ndim+1
                       write(ilun)nlevelmax
                       write(ilun)nstep_coarse
					   opened=.true.
                    endif
					write(ilun)ilevel
					write(ilun)nprint
					do idim=1,ndim
					   j=0
					   k=0
					   do ind=1, twotondim
				          iskip=ncoarse+(ind-1)*ngridmax
					      do i=1,ncache
					         if(son(ind_grid(i)+iskip)==0) then
							    k=k+1
							    if(stag(k) .eq. 1) then
								   j=j+1
							       if(idim .eq. 1) ttag(j)=ind_grid(i)+iskip
								   xdp(j)=xg(ind_grid(i),idim)+(cpi(ind,idim)-0.5)*dx
							    endif
						     endif
					      enddo
					   enddo
					   write(ilun)xdp
                    enddo
                    deallocate(stag)
                    deallocate(ind_grid)
					do ivar=1,ndim+1
					   if(ivar==1) then
                          do i=1,nprint
					         j=ttag(i)
					         fxdp(i)=uold(j,1)
					      enddo
					   else if(ivar>=2 .and. ivar<=ndim+1) then
					      do i=1,nprint
						     j=ttag(i)
						     fxdp(i)=uold(j,ivar)/max(uold(j,1),smallr)
					      enddo
					   endif
					   write(ilun)fxdp
                    enddo
#if NENER>0
              ! Write non-thermal pressures
                    do ivar=ndim+3,ndim+2+nener
                       do i=1,nprint
					      j=ttag(i)
                          fxdp(i)=(gamma_rad(ivar-ndim-2)-1d0)*uold(j,ivar)
                       enddo
                       write(ilun)fxdp
                    enddo
#endif
              ! Write thermal pressure
                    do i=1,nprint
					   j=ttag(i)
                       fxdp(i)=uold(j,ndim+2)
                       fxdp(i)=fxdp(i)-0.5d0*uold(j,2)**2/max(uold(j,1),smallr)
#if NDIM>1
                       fxdp(i)=fxdp(i)-0.5d0*uold(j,3)**2/max(uold(j,1),smallr)
#endif
#if NDIM>2
                       fxdp(i)=fxdp(i)-0.5d0*uold(j,4)**2/max(uold(j,1),smallr)
#endif
#if NENER>0
                       do irad=1,nener
                          fxdp(i)=fxdp(i)-uold(j,ndim+2+irad)
                       enddo
#endif
                       fxdp(i)=(gamma-1d0)*fxdp(i)
                    enddo
                    write(ilun)fxdp
#if NVAR>NDIM+2+NENER
              ! Write passive scalars
                    do ivar=ndim+3+nener,nvar
                       do i=1,nprint
						  j=ttag(i)
                          fxdp(i)=uold(j,ivar)/max(uold(j,1),smallr)
                       end do
                       write(ilun)fxdp
                    enddo
#endif
                    do i=1,nprint
					   j=ttag(i)
					   fxdp(i)=phi(j)
					enddo
					write(ilun)fxdp

				    deallocate(xdp)
				    deallocate(fxdp)
				    deallocate(ttag)
                 endif
              endif
			  if(allocated(ind_grid)) deallocate(ind_grid)
			  if(allocated(stag)) deallocate(stag)
           endif
        endif
     end do
  end do
  if(opened)write(ilun)endtag
  if(opened)write(ilun)endtag
  if(opened)close(ilun)


   if (myid == 1) then
      call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
      infofile = trim(spheredir)//'/info_'//trim(istep_str)
      open(ilun,file=TRIM(infofile)//".txt",form='formatted')
      rewind(ilun)
      write(ilun,'("a           =",E23.15)') aexp
      write(ilun,'("unit_l      =",E23.15)') scale_l
      write(ilun,'("unit_t      =",E23.15)') scale_t
      write(ilun,'("unit_d      =",E23.15)') scale_d
      write(ilun,'("unit_v      =",E23.15)') scale_v
      write(ilun,'("unit_nH     =",E23.15)') scale_nH
      write(ilun,'("unit_T2     =",E23.15)') scale_T2
      close(ilun)
   endif





  ! Send the token
#ifndef WITHOUTMPI
  if(IOGROUPSIZECONE>0) then
     if(mod(myid,IOGROUPSIZECONE)/=0 .and.(myid.lt.ncpu))then
        dummy_io=1
        call MPI_SEND(dummy_io,1,MPI_INTEGER,myid-1+1,tag, &
             & MPI_COMM_WORLD,info2)
     end if
  endif
#endif
  
#ifndef WITHOUTMPI
  call MPI_BARRIER(MPI_COMM_WORLD, info)
#endif


  
end subroutine output_sphere_hydro
