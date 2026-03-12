!=======================================================================
!=======================================================================
!=======================================================================
!=======================================================================
subroutine output_sphere_part(isphere)
  use amr_commons
  use pm_commons
  implicit none

#ifndef WITHOUTMPI
#include "mpif.h"
#endif
  
  integer::info,dummy_io,info2,nprint,ipart,ilun,idim,i,j,isphere
  integer,parameter::tag=1118

  character(len=2) :: sisphere
  character(len=5) :: istep_str
  character(len=150) :: spheredir, spherecmd, spherefile, infofile
  character(LEN=150)::fileloc
  character(LEN=5)::nchar

  integer(kind=2),allocatable,dimension(:)::stag
  integer(kind=4),allocatable,dimension(:)::ll,ttag
  integer(kind=8),dimension(:),allocatable::ii8
  real(dp),allocatable,dimension(:)::xdp
  real(kind=4),allocatable,dimension(:)::fxdp
  real(dp)::sradius,scenter(3)
  real(dp)::scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2
	


  ilun=3*ncpu+myid+10
  
  if(isphere .eq. 1) then
    sradius=sradius1
    scenter=scenter1
    sisphere="01"
    if(myid==1 .and. .not. hydro)write(*,*)'Dumping a spherical region 1'
  else if(isphere .eq. 2) then
    sradius=sradius2
    scenter=scenter2
    sisphere="02"
    if(myid==1 .and. .not. hydro)write(*,*)'Dumping a spherical region 2'
  else if(isphere .eq. 3) then
    sradius=sradius3
    scenter=scenter3
    sisphere="03"
    if(myid==1 .and. .not. hydro)write(*,*)'Dumping a spherical region 3'
  else if(isphere .eq. 4) then
    sradius=sradius4
    scenter=scenter4
    sisphere="04"
    if(myid==1 .and. .not. hydro)write(*,*)'Dumping a spherical region 4'
  else if(isphere .eq. 5) then
    sradius=sradius5
    scenter=scenter5
    sisphere="05"
    if(myid==1 .and. .not. hydro)write(*,*)'Dumping a spherical region 5'
  endif


  ! Determine the filename, dir, etc

  call title(nstep_coarse, istep_str)

  spheredir = "sphere/sphere_" // trim(sisphere) // "/nstep_" // trim(istep_str)

  spherecmd = "mkdir -p " // trim(spheredir)
  if(.not.withoutmkdir) then
    ! if (myid==1) call system(spherecmd)
     if (myid==1) call execute_command_line(spherecmd,wait=.true.)
  endif
  
#ifndef WITHOUTMPI
  call MPI_BARRIER(MPI_COMM_WORLD, info)
#endif

  spherefile = trim(spheredir)//'/sphere_part_'//trim(istep_str)//'.out'
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
  

   allocate(stag(1:npartmax))
   allocate(xdp(1:npart))

   stag=0
   nprint=0
   xdp=0.

      ipart=0
   do i=1,npartmax
	  if(levelp(i)>0) then
		 ipart=ipart+1
		 do idim=1,ndim
	        xdp(ipart)=xdp(ipart)+(xp(i,idim)-scenter(idim))**2
         enddo
	   	 xdp(ipart)=dsqrt(xdp(ipart))
		 if(xdp(ipart) .le. sradius) then
            stag(i)=1
			nprint=nprint+1
		 endif
      endif
   enddo
   

   deallocate(xdp)
   if(nprint .gt. 0) then

      allocate(ttag(1:nprint))
      j=0
      do i=1,npartmax
         if(stag(i) .eq. 1) then
	        j=j+1
            ttag(j)=i
	     endif
      enddo

      deallocate(stag)

      allocate(xdp(1:nprint))
      allocate(fxdp(1:nprint))

      open(ilun,file=TRIM(fileloc),form='unformatted')
      rewind(ilun)  
      write(ilun)ndim
      write(ilun)nstep_coarse
      write(ilun)nprint
   
      do idim=1,ndim
	     do i=1,nprint
	        j=ttag(i)
	        xdp(i)=xp(j,idim)
	     enddo
         write(ilun)xdp   
      enddo


      deallocate(xdp)
      do idim=1,ndim
	     do i=1,nprint
	        j=ttag(i)
	        fxdp(i)=vp(j,idim)
	     enddo
         write(ilun)fxdp   
      enddo

      do i=1,nprint
         j=ttag(i)
	     fxdp(i)=mp(j)
	  enddo
      write(ilun)fxdp   

      deallocate(fxdp)
	  allocate(ii8(1:nprint))

      do i=1,nprint
         j=ttag(i)
	     ii8(i)=idp(j)
	  enddo
      write(ilun)ii8
	  deallocate(ii8)

      allocate(ll(1:nprint))

      do i=1,nprint
         j=ttag(i)
	     ll(i)=levelp(j)
	  enddo
      write(ilun)ll
	  deallocate(ll)

#ifdef OUTPUT_PARTICLE_POTENTIAL
      allocate(fxdp(1:nprint))
      do i=1, nprint
         j=ttag(i)
         fxdp(i)=ptcl_phi(j)
      end do
      write(ilun)fxdp
	  deallocate(fxdp)
#endif

      if(star .or. sink) then
         allocate(fxdp(1:nprint))
         do i=1, nprint
            j=ttag(i)
            fxdp(i)=tp(j)
         end do
         write(ilun)fxdp
         
		 if(metal) then
            do i=1, nprint
               j=ttag(i)
               fxdp(i)=zp(j)
            end do
            write(ilun)fxdp
	     endif

         do i=1, nprint
            j=ttag(i)
            fxdp(i)=tpp(j)
         end do
         write(ilun)fxdp

         do i=1, nprint
            j=ttag(i)
            fxdp(i)=mp0(j)
         end do
         write(ilun)fxdp

         do i=1, nprint
            j=ttag(i)
            fxdp(i)=indtab(j)
         end do
         write(ilun)fxdp
		 deallocate(fxdp)
		 deallocate(ttag)
     endif

     close(ilun)

   endif

   if(allocated(stag)) deallocate(stag)


   if (myid == 1 .and. .not. hydro) then
      call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
      infofile = trim(spheredir)//'/info_'//trim(istep_str)
      open(ilun,file=TRIM(infofile)//".txt",form='formatted')
      rewind(ilun)
      write(ilun,'("a           =",E23.15)') aexp
      write(ilun,'("unit_l      =",E23.15)') scale_l
      write(ilun,'("unit_t      =",E23.15)') scale_t
      write(ilun,'("unit_d      =",E23.15)') scale_d
      write(ilun,'("unit_v      =",E23.15)') scale_v
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



end subroutine output_sphere_part
