!=======================================================================
!=======================================================================
!=======================================================================
!=======================================================================
subroutine output_cone_part(obs)
  use amr_commons
  use pm_commons
  implicit none

#ifndef WITHOUTMPI
#include "mpif.h"
#endif
  
  integer::info,dummy_io,info2
  integer,parameter::tag=1118

  character(len=5) :: istep_str
  character(len=100) :: conedir, conecmd, conefile
  
  integer::ilun,nx_loc,ipout,npout,npart_out
  character(LEN=80)::fileloc
  character(LEN=5)::nchar
  real(kind=8),dimension(1:3,1:nvector),save::pos,vel
!jhshin1
  real(kind=8),dimension(1:3,1:nvector),save::etc
  real(kind=8),dimension(:,:),allocatable::etcout
  real(sp),dimension(:,:),allocatable::etc_out
!jhshin2    
  real(kind=8),dimension(:,:),allocatable::posout,velout
  real(kind=8),dimension(:),allocatable::zout
  real(kind=8),dimension(:,:),allocatable::tmparr
  real(sp),dimension(:,:),allocatable::xp_out,vp_out
  real(sp),dimension(:),allocatable::zp_out
  real(kind=8) :: z1,z2,om0in,omLin,hubin,Lbox
  real(kind=8) :: Lobserver(3)
  integer::igrid,jgrid,ipart,jpart,idim,icpu,ilevel
  integer::i,ip,npart1
  integer::nalloc1,nalloc2

  integer,dimension(1:nvector),save::ind_part
!jaehyun1
  integer(kind=8),dimension(1:nvector),save::ii8
  integer(kind=8),dimension(:),allocatable::ii8out,ii8_out
  integer::elongated_axis_cone,obs
  real(kind=8) :: lboxz(3),minboxr(3),maxboxr(3),observer_cone(3),minboxr_cone(3),maxboxr_cone(3)
!jaehyun2
  
  logical::opened
  opened=.false.


  if(nstep_coarse.lt.2) return

  z2=1./aexp_old-1.
  z1=1./aexp-1.

  if(z1<0.) z1=0.
  if(z2.gt.zmax_cone) return
  if(abs(z2-z1)<1d-6) return
  
  
  om0in=omega_m
  omLin=omega_l
  hubin=h0/100.
  Lbox=boxlen_ini/hubin
  !observer=(/Lbox/2.0,Lbox*0.,Lbox/2.0/)

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
	


  ilun=3*ncpu+myid+10
  
  ! Determine the filename, dir, etc
  if(myid==1)write(*,*)'Computing and dumping lightcone (star particles)'

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

  conefile = trim(conedir)//'cone_part_'//trim(istep_str)//'.out'
  call title(myid,nchar)
  fileloc=TRIM(conefile)//TRIM(nchar)

 
  npart_out=0
  ipout=0
  npout=0

  ! Pre-allocate arrays for particle selection -----
  nalloc1=nvector
  allocate(posout(1:3, 1:nalloc1))
  allocate(velout(1:3, 1:nalloc1))
  !jhshin1
  allocate(etcout(1:3, 1:nalloc1))
  !jhshin2
  allocate(zout(1:nalloc1))
  !jaehyun1
  allocate(ii8out(1:nalloc1))
  !jaehyun2

  nalloc2=nvector+nstride
  allocate(xp_out(1:nalloc2,1:3))
  allocate(vp_out(1:nalloc2,1:3))
  !jhshin1
  allocate(etc_out(1:nalloc2,1:3))
  !jhshin2
  allocate(zp_out(1:nalloc2))
  allocate(tmparr(1:3, 1:nalloc2))
  !jaehyun1
  allocate(ii8_out(1:nalloc2))
  !jaehyun2
  ! ------------------------------------------------

  ! Wait for the token
#ifndef WITHOUTMPI
  if(IOGROUPSIZECONE>0) then
     if (mod(myid-1,IOGROUPSIZECONE)/=0) then
        call MPI_RECV(dummy_io,1,MPI_INTEGER,myid-1-1,tag,&
             & MPI_COMM_WORLD,MPI_STATUS_IGNORE,info2)
     end if
  endif
#endif
  

  

  do idim=1,ndim
    maxboxr(idim)=maxboxr_cone(idim)*Lbox
	minboxr(idim)=minboxr_cone(idim)*Lbox
	lboxz(idim)=maxboxr(idim)-minboxr(idim)
	Lobserver(idim)=Lobserver(idim)-minboxr(idim)
  enddo





  ilevel=levelmin_cone
  ! Loop over cpus
  do icpu=1,ncpu
     ! Loop over grids
     igrid=headl(icpu,ilevel)
     ip=0
     do jgrid=1,numbl(icpu,ilevel)
        npart1=numbp(igrid)  ! Number of particles in the grid
		

        if(npart1>0)then        
           ipart=headp(igrid)
           
           ! Loop over particles
           do jpart=1,npart1
              ip=ip+1
              ind_part(ip)=ipart
              if(ip==nvector)then
                 ! Lower left corner of 3x3x3 grid-cube
                 do idim=1,ndim
                    do i=1,ip
                       pos(idim,i)=xp(ind_part(i),idim)*Lbox-minboxr(idim)
                       vel(idim,i)=vp(ind_part(i),idim)
                    end do
                 end do

                 !jhshin1
				 do i=1,ip
  			  	    ii8(i)=idp(ind_part(i))
                    etc(1,i)=mp(ind_part(i))
                    etc(2,i)=tp(ind_part(i))
                    etc(3,i)=zp(ind_part(i))
				 end do
                !jhshin2

                 !===========================================================================
                 ! Count selection particles
                 !jhshin1
                 call perform_my_selection(.true.,z1,z2, &
                      &                           om0in,omLin,hubin,lboxz, &
                      &                           Lobserver,elongated_axis_cone, &
			 		  & 	        	     	 nvector, 3, &
                      &                           ii8,pos,vel,etc,ip, &
                      &                           ii8out,posout,velout,etcout,zout,npout,.false.)


                 call extend_arrays_if_needed()
                 ! Perform actual selection
                 call perform_my_selection(.false.,z1,z2, &
                      &                           om0in,omLin,hubin,lboxz, &
                      &                           Lobserver,elongated_axis_cone, &
			 		  & 	        	     	 nvector, 3, &
                      &                           ii8,pos,vel,etc,ip, &
                      &                           ii8out,posout,velout,etcout,zout,npout,.false.)
                 !jhshin2
                 !===========================================================================
                 if(npout>0)then
                    do idim=1,ndim
                       do i=1,npout
                          xp_out(ipout+i,idim)=posout(idim,i)/Lbox
                          vp_out(ipout+i,idim)=velout(idim,i)
                       end do
                    end do
                    do i=1,npout
                       zp_out(ipout+i)=zout(i)
                       !jhshin1
					   ii8_out(ipout+i)=ii8out(i)
                       etc_out(ipout+i,1)=etcout(1,i)
                       etc_out(ipout+i,2)=etcout(2,i)
                       etc_out(ipout+i,3)=etcout(3,i)
                       !jhshin2

                    end do
                    ipout=ipout+npout
                    npart_out=npart_out+npout
                 endif

                 ip=0
              end if
              if(ipout>=nstride)then
                 if(.not.opened) then
                    open(ilun,file=TRIM(fileloc),form='unformatted')
                    rewind(ilun)  
                    write(ilun)ncpu
            		write(ilun)nstride
		            write(ilun)npart
                    opened=.true.
                 endif
				 !jaehyun1
				 write(ilun)ii8_out(1:nstride)
				 !jaehyun2
                 do idim=1,ndim
                    write(ilun)xp_out(1:nstride,idim)
                    write(ilun)vp_out(1:nstride,idim)
                 end do
                 !jhshin1
                 do idim=1,3
                    write(ilun)etc_out(1:nstride,idim)
                 end do
                 !jhshin2
                 write(ilun)zp_out(1:nstride)
                 do idim=1,ndim
                    do i=1,ipout-nstride
                       xp_out(i,idim)=xp_out(i+nstride,idim)
                       vp_out(i,idim)=vp_out(i+nstride,idim)
                    end do
                 end do
                 !jhshin1
                 do idim=1,3
                   do i=1,ipout-nstride
                       etc_out(i,idim)=etc_out(i+nstride,idim)
                   end do
                 end do
                 !jhshin2
                 do i=1,ipout-nstride
                    zp_out(i)=zp_out(i+nstride)
					!jaehyun1
			        ii8_out(i)=ii8_out(i+nstride)
					!jaehyun2
                 end do
                 ipout=ipout-nstride
              endif
              ipart=nextp(ipart)  ! Go to next particle
           end do
           ! End loop over particles           
        end if
        igrid=next(igrid)   ! Go to next grid
     end do
     ! End loop over grids

     if(ip>0)then
        ! Lower left corner of 3x3x3 grid-cube
        do idim=1,ndim
           do i=1,ip
              pos(idim,i)=xp(ind_part(i),idim)*Lbox-minboxr(idim)
              vel(idim,i)=vp(ind_part(i),idim)
           end do
        end do
        !jhshin1
		do i=1,ip
           ii8(i)=idp(ind_part(i))
           etc(1,i)=mp(ind_part(i))
           etc(2,i)=tp(ind_part(i))
           etc(3,i)=zp(ind_part(i))
		end do
        !jhshin2

        !===========================================================================
        !jhshin1
        ! Count selection particles
        call perform_my_selection(.true.,z1,z2, &
             &                           om0in,omLin,hubin,lboxz, &
             &                           Lobserver,elongated_axis_cone, &
			 &	 	        	     	 nvector, 3, &
             &                           ii8,pos,vel,etc,ip, &
             &                           ii8out,posout,velout,etcout,zout,npout,.false.)

        call extend_arrays_if_needed()
            
        ! Perform actual selection
        call perform_my_selection(.false.,z1,z2, &
             &                           om0in,omLin,hubin,lboxz, &
             &                           Lobserver,elongated_axis_cone, &
			 & 		        	     	 nvector, 3, &
             &                           ii8,pos,vel,etc,ip, &
             &                           ii8out,posout,velout,etcout,zout,npout,.false.)

        !jhshin2
        !===========================================================================
        if(npout>0)then
           do idim=1,ndim
              do i=1,npout
                 xp_out(ipout+i,idim)=posout(idim,i)/Lbox
                 vp_out(ipout+i,idim)=velout(idim,i)
              end do
           end do
           do i=1,npout
              zp_out(ipout+i)=zout(i)
			  !jaehyun1
			  ii8_out(ipout+i)=ii8out(i)
			  !jaehyun2
              !jhshin1
              etc_out(ipout+i,1)=etcout(1,i)
              etc_out(ipout+i,2)=etcout(2,i)
              etc_out(ipout+i,3)=etcout(3,i)
              !jhshin2
           end do
           ipout=ipout+npout
           npart_out=npart_out+npout
        endif
     endif
     if(ipout>=nstride)then
        if(.not.opened) then
           open(ilun,file=TRIM(fileloc),form='unformatted')
           rewind(ilun)  
           write(ilun)ncpu
           write(ilun)nstride
           write(ilun)npart
           opened=.true.
        endif
	 !jaehyun1
		write(ilun)ii8_out(1:nstride)
	 !jaehyun2
        do idim=1,ndim
           write(ilun)xp_out(1:nstride,idim)
           write(ilun)vp_out(1:nstride,idim)
        end do
        !jhshin1
        do idim=1,3
           write(ilun)etc_out(1:nstride,idim)
        end do
        !jhshin2
        write(ilun)zp_out(1:nstride)
        do idim=1,ndim
           do i=1,ipout-nstride
              xp_out(i,idim)=xp_out(i+nstride,idim)
              vp_out(i,idim)=vp_out(i+nstride,idim)
           end do
        end do
        !jhshin1
        do idim=1,ndim
           do i=1,ipout-nstride
              etc_out(i,idim)=etc_out(i+nstride,idim)
            end do
        end do
        !jhshin2
        do i=1,ipout-nstride
           zp_out(i)=zp_out(i+nstride)
		   !jaehyun1
		   ii8_out(i)=ii8_out(i+nstride)
		   !jaehyun2
        end do
        ipout=ipout-nstride
     endif
  end do
  ! End loop over cpus

  if(ipout>0)then
     if(.not.opened) then
        open(ilun,file=TRIM(fileloc),form='unformatted')
        rewind(ilun)  
        write(ilun)ncpu
        write(ilun)nstride
        write(ilun)npart
        opened=.true.
     endif
	 !jaehyun1
  	 write(ilun)ii8_out(1:nstride)
	 !jaehyun2
     do idim=1,ndim
        write(ilun)xp_out(1:ipout,idim)
        write(ilun)vp_out(1:ipout,idim)
     end do
     !jhshin1
     do idim=1,ndim
        write(ilun)etc_out(1:ipout,idim)
     end do
     !jhshin2
     write(ilun)zp_out(1:ipout)
  endif

  if(opened)close(ilun)
  
  if (verbose)write(*,*)'cone output=',myid,npart_out

  if(npart_out>0) then
     open(ilun,file=TRIM(fileloc)//".txt",form='formatted')
     rewind(ilun)
     write(ilun,*) ncpu
     write(ilun,*) nstride
     write(ilun,*) npart_out
!     write(ilun,*) aexp_old
!     write(ilun,*) 1./(1.+z1)
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
  

   if((opened.and.(npart_out==0)).or.((.not.opened).and.(npart_out>0))) then
     write(*,*)'Error in output_cone'
     write(*,*)'npart_out=',npart_out,'opened=',opened
     stop
  endif


contains

    ! Extends (deallocates and reallocates) the arrays
    ! posout, velout, zout, xp_out, vp_out and zp_out
    ! after npout has been updated, so they can hold enough particles 
    !
    ! Reallocation is done in chunks of size alloc_chunk_size, to avoid
    ! reallocating too frequently.

    subroutine extend_arrays_if_needed()

        ! Allocation chunk size
        integer, parameter :: alloc_chunk_size = 100
        integer :: new_nalloc1, new_nalloc2
        integer :: nchunks1, nchunks2

        if (nalloc1 >= npout .and. nalloc2 >= npout+nstride) return


        ! Compute new array sizes
        nchunks1 = npout / alloc_chunk_size
        if (mod(npout, alloc_chunk_size) > 0) nchunks1=nchunks1+1

        nchunks2 = (npout+nstride) / alloc_chunk_size
        if (mod(npout+nstride, alloc_chunk_size) > 0) nchunks2=nchunks2+1

        new_nalloc1 = nchunks1 * alloc_chunk_size
        new_nalloc2 = nchunks2 * alloc_chunk_size

        ! Resize temp array
        deallocate(tmparr)
        allocate(tmparr(1:3,1:max(new_nalloc1,new_nalloc2)))


        ! Resize xp_out, vp_out, zp_out
        do idim=1,ndim
            tmparr(idim,1:nalloc2)=xp_out(1:nalloc2,idim)
        end do
        deallocate(xp_out); allocate(xp_out(1:new_nalloc2,1:3))
        do idim=1,ndim
            xp_out(1:nalloc2,idim)=tmparr(idim,1:nalloc2)
        end do

        do idim=1,ndim
            tmparr(idim,1:nalloc2)=vp_out(1:nalloc2,idim) 
        end do
        deallocate(vp_out); allocate(vp_out(1:new_nalloc2,1:3))
        do idim=1,ndim
            vp_out(1:nalloc2,idim)=tmparr(idim,1:nalloc2)
        end do

        !jhshin1
        do idim=1,ndim
            tmparr(idim,1:nalloc2)=etc_out(1:nalloc2,idim) 
        end do
        deallocate(etc_out); allocate(etc_out(1:new_nalloc2,1:3))
        do idim=1,ndim
            etc_out(1:nalloc2,idim)=tmparr(idim,1:nalloc2)
        end do
        !jhshin2

		!jaehyun1
		tmparr(1,1:nalloc2)=ii8_out(1:nalloc2)
        deallocate(ii8_out); allocate(ii8_out(1:new_nalloc2))
		ii8_out(1:nalloc2)=tmparr(1,1:nalloc2)
		!jaehyun2

        tmparr(1,1:nalloc2)=zp_out(1:nalloc2) 
        deallocate(zp_out); allocate(zp_out(1:new_nalloc2))
        zp_out(1:nalloc2)=tmparr(1,1:nalloc2)

        nalloc2 = new_nalloc2


        ! Resize posout, velout, zout
        do idim=1,ndim
            tmparr(idim,1:nalloc1)=posout(idim,1:nalloc1) 
        deallocate(posout); allocate(posout(1:3,1:new_nalloc1))
        end do
        do idim=1,ndim
            posout(idim,1:nalloc1)=tmparr(idim,1:nalloc1)
        end do

        do idim=1,ndim
            tmparr(idim,1:nalloc1)=velout(idim,1:nalloc1) 
        end do
        deallocate(velout); allocate(velout(1:3,1:new_nalloc1))
        do idim=1,ndim
            velout(idim,1:nalloc1)=tmparr(idim,1:nalloc1)
        end do

        !jhshin1
        do idim=1,ndim
            tmparr(idim,1:nalloc1)=etcout(idim,1:nalloc1) 
        end do
        deallocate(etcout); allocate(etcout(1:3,1:new_nalloc1))
        do idim=1,ndim
            etcout(idim,1:nalloc1)=tmparr(idim,1:nalloc1)
        end do
        !jhshin2

		!jaehyun1
        tmparr(1,1:nalloc1)=ii8out(1:nalloc1) 
        deallocate(ii8out); allocate(ii8out(1:new_nalloc1))
        ii8out(1:nalloc1)=tmparr(1,1:nalloc1)
		!jaehyun2

        tmparr(1,1:nalloc1)=zout(1:nalloc1) 
        deallocate(zout); allocate(zout(1:new_nalloc1))
        zout(1:nalloc1)=tmparr(1,1:nalloc1)

        nalloc1 = new_nalloc1

    end subroutine extend_arrays_if_needed
end subroutine output_cone_part
