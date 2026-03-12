!=======================================================================
!=======================================================================
!=======================================================================
!=======================================================================
subroutine output_sphere_sink(isphere)
  use amr_commons
  use pm_commons
  implicit none

#ifndef WITHOUTMPI
#include "mpif.h"
#endif
  
  integer::info,dummy_io,info2,nprint,ilun,idim,i,j,ilevel,isphere
  integer,parameter::tag=1118

  character(len=2) :: sisphere
  character(len=5) :: istep_str
  character(len=150) :: spheredir, spherecmd, spherefile
  character(LEN=150)::fileloc
  character(LEN=5)::nchar

  integer(kind=2),allocatable,dimension(:)::stag
  integer(kind=4),allocatable,dimension(:)::ll,ttag
  real(dp),allocatable,dimension(:)::xdp

  real(dp)::sradius,scenter(3)


	
  if(.not. sink) return

  ilun=3*ncpu+myid+10
 
  if(isphere .eq. 1) then
    sradius=sradius1
    scenter=scenter1
    sisphere="01"
  else if(isphere .eq. 2) then
    sradius=sradius2
    scenter=scenter2
    sisphere="02"
  else if(isphere .eq. 3) then
    sradius=sradius3
    scenter=scenter3
    sisphere="03"
  else if(isphere .eq. 4) then
    sradius=sradius4
    scenter=scenter4
    sisphere="04"
  else if(isphere .eq. 5) then
    sradius=sradius5
    scenter=scenter5
    sisphere="05"
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

  fileloc = trim(spheredir)//'/sphere_sink_'//trim(istep_str)//'.out'
!  call title(myid,nchar)
!  fileloc=TRIM(spherefile)//TRIM(nchar)

 

  ! Wait for the token
!#ifndef WITHOUTMPI
!  if(IOGROUPSIZEOUT>0) then
!     if (mod(myid-1,IOGROUPSIZEOUT)/=0) then
!        call MPI_RECV(dummy_io,1,MPI_INTEGER,myid-1-1,tag,&
!             & MPI_COMM_WORLD,MPI_STATUS_IGNORE,info2)
!     end if
!  endif
!#endif
  
   if(myid == 1) then 
      allocate(stag(1:nsink))
      allocate(xdp(1:nsink))

      stag=0
      nprint=0
      xdp=0.

      do i=1,nsink
	     do idim=1,ndim
	        xdp(i)=xdp(i)+(xsink(i,idim)-scenter(idim))**2
         enddo
   	     xdp(i)=dsqrt(xdp(i))
	     if(xdp(i) .le. sradius) then
            stag(i)=1
	        nprint=nprint+1
	     endif
      enddo

      deallocate(xdp)
      if(nprint .gt. 0) then

         allocate(ttag(1:nprint))
         j=0
         do i=1,nsink
            if(stag(i) .eq. 1) then
	           j=j+1
               ttag(j)=i
	        endif
         enddo

         deallocate(stag)


         open(ilun,file=TRIM(fileloc),form='unformatted')
         rewind(ilun)  
		 write(ilun)ndim
		 write(ilun)nstep_coarse
         write(ilun)nprint
   
         allocate(ll(1:nprint))

         do i=1,nprint
            j=ttag(i)
	        ll(i)=idsink(j)
	     enddo
         write(ilun)ll
	     deallocate(ll)
	 
         allocate(xdp(1:nprint))
         do i=1,nprint
            j=ttag(i)
	        xdp(i)=msink(j)
	     enddo
         write(ilun)xdp   


         do idim=1,ndim
	        do i=1,nprint
	           j=ttag(i)
	           xdp(i)=xsink(j,idim)
	        enddo
            write(ilun)xdp   
         enddo


         do idim=1,ndim
	        do i=1,nprint
	           j=ttag(i)
	           xdp(i)=vsink(j,idim)
	        enddo
            write(ilun)xdp   
         enddo


         do i=1, nprint
            j=ttag(i)
            xdp(i)=tsink(j)
         end do
         write(ilun)xdp

         do i=1, nprint
            j=ttag(i)
            xdp(i)=dMsmbh(j)
         end do
         write(ilun)xdp
        
         do i=1, nprint
            j=ttag(i)
            xdp(i)=dMBH_coarse(j)
         end do
         write(ilun)xdp

         do i=1, nprint
            j=ttag(i)
            xdp(i)=dMEd_coarse(j)
         end do
         write(ilun)xdp

         do idim=1,ndim
	        do i=1,nprint
	           j=ttag(i)
	           xdp(i)=jsink(j,idim)
	        enddo
            write(ilun)xdp   
         enddo

         do idim=1,ndim
	        do i=1,nprint
	           j=ttag(i)
	           xdp(i)=bhspin(j,idim)
	        enddo
            write(ilun)xdp   
         enddo

         do i=1, nprint
            j=ttag(i)
            xdp(i)=spinmag(j)
         end do
         write(ilun)xdp
         do i=1, nprint
            j=ttag(i)
            xdp(i)=eps_sink(j)
         end do
         write(ilun)xdp

         do idim=1,ndim*2+1
		    do ilevel=levelmin,nlevelmax
	           do i=1,nprint
	              j=ttag(i)
	              xdp(i)=sink_stat(j,ilevel,idim)
			   enddo
               write(ilun)xdp   
	        enddo
         enddo

		 deallocate(xdp)
         close(ilun)
	  endif
      if(allocated(stag)) deallocate(stag)
   endif

#ifndef WITHOUTMPI
  call MPI_BARRIER(MPI_COMM_WORLD, info)
#endif



     ! Send the token
!#ifndef WITHOUTMPI
!  if(IOGROUPSIZEOUT>0) then
!     if(mod(myid,IOGROUPSIZEOUT)/=0 .and.(myid.lt.ncpu))then
!        dummy_io=1
!        call MPI_SEND(dummy_io,1,MPI_INTEGER,myid-1+1,tag, &
!             & MPI_COMM_WORLD,info2)
!     end if
!  endif
!#endif
  



end subroutine output_sphere_sink
