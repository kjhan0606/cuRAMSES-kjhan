!=======================================================================
!=======================================================================
!=======================================================================
!=======================================================================
subroutine output_cone_sink(obs)
  use amr_commons
  use pm_commons
  implicit none

#ifndef WITHOUTMPI
#include "mpif.h"
#endif
  
  integer::info,dummy_io,info2
  integer,parameter::tag=1118

  character(len=5) :: istep_str
  character(len=150) :: conedir, conecmd, conefile
  
  integer::ilun,npout,ns,i,j,k,idim,mnsink,nsink_out,end_tag,ilevel
  real(kind=8),dimension(:,:),allocatable::pos,vel,etc
  real(kind=8),dimension(:,:),allocatable::posout,velout,etcout
  real(kind=8),dimension(:),allocatable::zout
  real(kind=8),dimension(:,:),allocatable::pos_out,vel_out,etc_out
  real(kind=8) :: z1,z2,om0in,omLin,hubin,Lbox
  real(kind=8) :: Lobserver(3),nfine,daexp
  integer(kind=8),dimension(:),allocatable::ii8
  integer(kind=8),dimension(:),allocatable::ii8out
  integer::elongated_axis_cone,obs
  real(kind=8) :: lboxz(3),minboxr(3),maxboxr(3),observer_cone(3),minboxr_cone(3),maxboxr_cone(3)
  logical::opened
  opened=.false.


  if(nstep_coarse.lt.2) return

  nfine=1.
  do ilevel=levelmin,nlevelmax
     if(numbtot(1,ilevel) .gt. 0) then
        nfine=nfine*nsubcycle(ilevel-1)
     endif
  enddo


  daexp=max((aexp-aexp_old2)/nfine,aexp-aexp_old_fine)

  z2=1./(aexp_old2-daexp)-1.
  z1=1./(aexp+daexp)-1.


  if(z1<0.) z1=0.
  if(z2.gt.zmax_cone) return
  if(abs(z2-z1)<1d-6) return
  

  om0in=omega_m
  omLin=omega_l
  hubin=h0/100.
  Lbox=boxlen_ini/hubin
!  Lobserver=(/Lbox/2.0,Lbox*0.,Lbox/2.0/)
!  observer=(/Lbox*0.5,Lbox*0.5,Lbox*0.5/)


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
  if(myid==1 .and. obs==1 )write(*,*)'Computing and dumping lightcone (sink particles)'

  call title(nstep_coarse, istep_str)

  if(obs==1) then
    conedir = "light_cone/cone_" // trim(istep_str) // "/observer1/"
  else
    conedir = "light_cone/cone_" // trim(istep_str) // "/observer2/"
  endif

  conecmd = "mkdir -p " // trim(conedir)
  if(.not.withoutmkdir) then
    ! if (myid==1) call system(conecmd)
     if (myid==1) call execute_command_line(conecmd,wait=.true.)
  endif
  
#ifndef WITHOUTMPI
  call MPI_BARRIER(MPI_COMM_WORLD, info)
#endif

  conefile = trim(conedir)//'cone_sink_'//trim(istep_str)//'.out'
!  call title(myid,nchar)
!  fileloc=TRIM(conefile)//TRIM(nchar)

  end_tag=-1
  npout=0
  j=0
  mnsink=65536
  nsink_out=0

  allocate(pos(1:3, 1:mnsink))
  allocate(vel(1:3, 1:mnsink))
  allocate(etc(1:13, 1:mnsink))
  allocate(ii8(1:mnsink))


  ! Wait for the token
!#ifndef WITHOUTMPI
!  if(IOGROUPSIZECONE>0) then
!     if (mod(myid-1,IOGROUPSIZECONE)/=0) then
!        call MPI_RECV(dummy_io,1,MPI_INTEGER,myid-1-1,tag,&
!             & MPI_COMM_WORLD,MPI_STATUS_IGNORE,info2)
!     end if
!  endif
!#endif


      
   do idim=1,ndim
     maxboxr(idim)=maxboxr_cone(idim)*Lbox
     minboxr(idim)=minboxr_cone(idim)*Lbox
     lboxz(idim)=maxboxr(idim)-minboxr(idim)
     Lobserver(idim)=Lobserver(idim)-minboxr(idim)
   enddo


  if (myid == 1) then 

    if (nsink>0) then
       if(obs==1) write(*,*)'Number of sink particles = ',nsink
    ! Loop over particles
       do ns=1,nsink
          j=j+1
          do idim=1,ndim
             pos(idim,j)=xsink(ns,idim)*Lbox-minboxr(idim)
             vel(idim,j)=vsink(ns,idim)
          enddo
          ii8(j)=idsink(ns)
          etc(1,j)=msink(ns)
          etc(2,j)=tsink(ns)
          etc(3,j)=dMsmbh(ns)
          etc(4,j)=dMBH_coarse(ns)
          etc(5,j)=dMEd_coarse(ns)
		  do idim=1,ndim
             etc(6+idim-1,j)=jsink(ns,idim)
			 etc(9+idim-1,j)=bhspin(ns,idim)
		  enddo
          etc(12,j)=spinmag(ns)
          etc(13,j)=eps_sink(ns)

          if (j==mnsink .or. ns==nsink) then
             !===========================================================================
             ! Count selection particles
             call perform_my_selection(.true.,z1,z2, &
                  &                           om0in,omLin,hubin,lboxz, &
                  &                           Lobserver,elongated_axis_cone, &
                  &                           j, 13, &
                  &                           ii8,pos,vel,etc,j, &
                  &                           ii8out,posout,velout,etcout,zout,npout,.false.)


             if (npout > 0) then
                allocate(ii8out(1:npout))
                allocate(posout(1:ndim,1:npout))
                allocate(velout(1:ndim,1:npout))
                allocate(etcout(1:13,1:npout))
                allocate(zout(1:npout))
                allocate(pos_out(1:npout,1:ndim))
                allocate(vel_out(1:npout,1:ndim))
                allocate(etc_out(1:npout,1:13))
                

!          call extend_arrays_if_needed_2()
          ! Perform actual selection
                call perform_my_selection(.false.,z1,z2, &
                     &                           om0in,omLin,hubin,lboxz, &
                     &                           Lobserver,elongated_axis_cone, &
                     &                           j, 13, &
                     &                           ii8,pos,vel,etc,j, &
                     &                           ii8out,posout,velout,etcout,zout,npout,.false.)
             !===========================================================================
                nsink_out=nsink_out+npout
                do idim=1,ndim
                   do i=1,npout
                      pos_out(i,idim)=posout(idim,i)/Lbox
                      vel_out(i,idim)=velout(idim,i)
                      if (idim==1) then
						 do k=1,13
                            etc_out(i,k)=etcout(k,i)
						 enddo
                      endif
                   end do
                end do

                if (.not.opened) then
                   open(ilun,file=TRIM(conefile),form='unformatted')
                   rewind(ilun)  
                   opened=.true.
                endif
                write(ilun)npout
                write(ilun)ii8out(1:npout)
                do idim=1,ndim
                   write(ilun)pos_out(1:npout,idim)
                   write(ilun)vel_out(1:npout,idim)
                end do
                do idim=1,13
                   write(ilun)etc_out(1:npout,idim)
                end do
                write(ilun)zout(1:npout)

                if(allocated(ii8out)) deallocate(ii8out)
                if(allocated(posout)) deallocate(posout)
                if(allocated(velout)) deallocate(velout)
                if(allocated(etcout)) deallocate(etcout)
                if(allocated(zout)) deallocate(zout)
                if(allocated(pos_out)) deallocate(pos_out)
                if(allocated(vel_out)) deallocate(vel_out)
                if(allocated(etc_out)) deallocate(etc_out)
              endif
              j=0
           endif
        enddo
     endif
     if(opened) write(ilun)end_tag
     if(opened) close(ilun)


    if(nsink_out>0) then
       open(ilun,file=TRIM(conefile)//".txt",form='formatted')
       rewind(ilun)
       write(ilun,*) nsink_out
       close(ilun)
     endif
  endif







     ! Send the token
!#ifndef WITHOUTMPI
!  if(IOGROUPSIZECONE>0) then
!     if(mod(myid,IOGROUPSIZECONE)/=0 .and.(myid.lt.ncpu))then
!        dummy_io=1
!        call MPI_SEND(dummy_io,1,MPI_INTEGER,myid-1+1,tag, &
!             & MPI_COMM_WORLD,info2)
!     end if
!  endif
!#endif
  

   if((opened.and.(nsink_out==0)).or.((.not.opened).and.(nsink_out>0))) then
     write(*,*)'Error in output_cone'
     write(*,*)'npart_out=',nsink_out,'opened=',opened
     stop
  endif

#ifndef WITHOUTMPI
  call MPI_BARRIER(MPI_COMM_WORLD, info)
#endif



end subroutine output_cone_sink
