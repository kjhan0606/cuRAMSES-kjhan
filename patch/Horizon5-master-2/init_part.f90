subroutine init_part
  use amr_commons
  use pm_commons
  use hydro_parameters, ONLY: ichem
  use clfind_commons
  use ksection

#ifdef RT
  use rt_parameters,only: convert_birth_times
#endif
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  !------------------------------------------------------------
  ! Allocate particle-based arrays.
  ! Read particles positions and velocities from grafic files
  !------------------------------------------------------------
  integer::npart2,ndim2,ncpu2
  integer::ipart,jpart,ipart_old,ilevel,idim
  integer::i,igrid,ncache,ngrid,iskip,isink
  integer::ind,ix,iy,iz,ilun,info,icpu,nx_loc
  integer::i1,i2,i3,i1_min,i1_max,i2_min,i2_max,i3_min,i3_max
  integer(kind=8)::byte_pos,hdr_bytes,plane_bytes
  integer::buf_count,indglob,npart_new,nDMloc
  integer(i8b)::tmp_long
  real(dp)::dx,xx1,xx2,xx3,vv1,vv2,vv3,mm1,ll1,ll2,ll3
  real(dp)::scale,dx_loc,rr,rmax,dx_min
  integer::ncode,bit_length,temp
  real(kind=8)::bscale
  real(dp),dimension(1:twotondim,1:3)::xc
  integer ,dimension(1:nvector)::ind_grid,ind_cell,cc,ii
  integer(i8b),dimension(1:ncpu)::npart_cpu,npart_all
  real(dp),allocatable,dimension(:)::xdp
  integer,allocatable,dimension(:)::isp
  integer(i8b),allocatable,dimension(:)::isp8
  logical,allocatable,dimension(:)::nb
  real(kind=4),allocatable,dimension(:,:)::init_plane,init_plane_x
  real(kind=4),allocatable,dimension(:,:,:)::init_array,init_array_x
  real(kind=8),dimension(1:nvector,1:3)::xx,vv,xs
  real(dp),dimension(1:nvector,1:3)::xx_dp
  integer,dimension(1:nvector)::ixx,iyy,izz
  real(qdp),dimension(1:nvector)::order
  real(kind=8),dimension(1:nvector)::mm
  real(kind=8)::dispmax=0.0,dispall
  real(dp),dimension(1:3)::skip_loc
  real(dp),dimension(1:3)::centerofmass

  integer::ibuf,tag=101,tagf=102,tagu=102
  integer::countsend,countrecv
#ifndef WITHOUTMPI
  integer,dimension(MPI_STATUS_SIZE,2*ncpu)::statuses
  integer,dimension(2*ncpu)::reqsend,reqrecv
  integer,dimension(ncpu)::sendbuf,recvbuf
  integer::ntotal_ksec,nrecv_ksec,idx_ksec
  real(dp),allocatable::sbuf_ksec(:,:),rbuf_ksec(:,:)
  integer,allocatable::dcpu_ksec(:)
#endif

  logical::error,keep_part,eof,jumped,ic_sink=.false.,read_pos=.false.,ok
  character(LEN=80)::filename,filename_x
  character(LEN=80)::fileloc
  character(LEN=20)::filetype_loc
  character(LEN=5)::nchar,ncharcpu
  integer,parameter::tagg=1109,tagg2=1110,tagg3=1111
  integer::dummy_io,info2

  integer::nsteps,nmetal,j,k,ielt
  real::xc1!,xc2,xc3,xc4,xc5,xc6,xc7,xc8,xc9,xc10,xc11,xc12,xc13,xc14,xc15,xc16
  real,allocatable,dimension(:)::yield_table_readline
  real(dp)::astarconv
  real(dp) :: scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2
  character(LEN=1)::a1
  character(LEN=2)::nelt_str
  character(LEN=11)::format_string
  logical::file_exist

  real(kind=4):: real_mem,real_mem_tot


  if(verbose)write(*,*)'Entering init_part'

  if(allocated(xp))then
     if(verbose)write(*,*)'Initial conditions already set'
     return
  end if

  ! Allocate particle variables
  allocate(xp    (npartmax,ndim))
  allocate(vp    (npartmax,ndim))
  allocate(mp    (npartmax))
  allocate(nextp (npartmax))
  allocate(prevp (npartmax))
  allocate(levelp(npartmax))
  allocate(idp   (npartmax))
#ifdef OUTPUT_PARTICLE_POTENTIAL
  allocate(ptcl_phi(npartmax))
#endif
  xp=0.0; vp=0.0; mp=0.0; levelp=0; idp=0
  if(star.or.sink)then
     allocate(tp(npartmax))
     tp=0.0
     if(metal)then
        allocate(zp(npartmax))
        zp=0.0
     end if


     !Read the yield table
     open(33,file=trim(yieldtablefilename),status='old', form='formatted')
     read(33,'(a8,I)')a1,nmetal
     read(33,'(a8,I)')a1,nsteps
     read(33,'(a11,I)')a1,nelt
     yieldtab%na=nsteps
     yieldtab%nz=nmetal
     allocate(elem_list(nelt))
     write(nelt_str,'(I2.2)') nelt
     if(nelt.gt.0)then
        format_string = '(a15,'//trim(nelt_str)//'a3)'
     else
        format_string = '(a15)'
     endif
     read(33,format_string)a1,elem_list

     allocate(tpp(npartmax))
     tpp=0.0
     allocate(mp0(npartmax))
     mp0=0.d0
     allocate(indtab(npartmax))
     indtab=0.d0
!     allocate(cep(npartmax,nelt))
!     cep=0.d0
     call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
     astarconv = 1d9*365.*24.*3600.*scale_t
     allocate(yieldtab%astar(nsteps),yieldtab%zstar(nmetal))
     allocate(yieldtab%Eeject(nsteps,nmetal,nelt))
     allocate(yieldtab%Zeject(nsteps,nmetal),yieldtab%Meject(nsteps,nmetal))
     allocate(yieldtab%NSN1(nsteps,nmetal),yieldtab%NSN2(nsteps,nmetal))
     allocate(yield_table_readline(nelt+5))
     write(nelt_str,'(I2.2)') nelt+5
     format_string = '('//trim(nelt_str)//'e19.11)'
     do i=1,nmetal
        read(33,'(a1,e16.5)')a1,xc1
        yieldtab%zstar(i)=xc1
        do j=1,nsteps
           !time, NSNII, NSNIa, Mej, MZ, and 8 elements (the names and order are stored in elem_list)
           read(33,format_string) yield_table_readline
           yieldtab%astar(j)=yield_table_readline(1)
           yieldtab%NSN2(j,i)=yield_table_readline(2)
           yieldtab%NSN1(j,i)=yield_table_readline(3)
           yieldtab%Meject(j,i)=yield_table_readline(4)
           yieldtab%Zeject(j,i)=yield_table_readline(5)
           do ielt=1,nelt
              yieldtab%Eeject(j,i,ielt)=yield_table_readline(5+ielt)
           enddo
        enddo
     enddo
     close(33)

     ilun=2*ncpu+myid+10

     if(nrestart>0) then
        call title(nrestart,nchar) !I think if this is nrestart = 0 or 1 we need to reset the output1 ones 
     else
        call title(nrestart+1,nchar) !I think if this is nrestart = 0 or 1 we need to reset the output1 ones 
     endif

     if(IOGROUPSIZEREP>0)then
        call title(((myid-1)/IOGROUPSIZEREP)+1,ncharcpu)
        fileloc = 'output_'//TRIM(nchar)//'/group_'//TRIM(ncharcpu)//'/stars_'//TRIM(nchar)//'.out'
     else
        fileloc = 'output_'//TRIM(nchar)//'/stars_'//TRIM(nchar)//'.out'
     endif

     call title(myid,nchar)
     fileloc = TRIM(fileloc)//TRIM(nchar)
     inquire(file=fileloc,exist=file_exist)
     if(file_exist) then
     !resetting the star files on a restart
        open(unit=ilun, file=fileloc, form='formatted')
        close(ilun,STATUS='DELETE')
     endif

  end if

  !--------------------
  ! Read part.tmp file
  !--------------------
  if(nrestart>0)then

     ilun=2*ncpu+myid+10
     call title(nrestart,nchar)

     if(IOGROUPSIZEREP>0)then
        call title(((myid-1)/IOGROUPSIZEREP)+1,ncharcpu)
        fileloc='output_'//TRIM(nchar)//'/group_'//TRIM(ncharcpu)//'/part_'//TRIM(nchar)//'.out'
     else
        fileloc='output_'//TRIM(nchar)//'/part_'//TRIM(nchar)//'.out'
     endif

     call title(myid,nchar)
     fileloc=TRIM(fileloc)//TRIM(nchar)
     ! Wait for the token                                                                                                                                                                    
#ifndef WITHOUTMPI
     if(IOGROUPSIZE>0) then
        if (mod(myid-1,IOGROUPSIZE)/=0) then
           call MPI_RECV(dummy_io,1,MPI_INTEGER,myid-1-1,tagg,&
                & MPI_COMM_WORLD,MPI_STATUS_IGNORE,info2)
        end if
     endif
#endif

     open(unit=ilun,file=fileloc,form='unformatted')
     rewind(ilun)
     read(ilun)ncpu2
     read(ilun)ndim2
     read(ilun)npart2
     read(ilun)localseed
     read(ilun)nstar_tot
     read(ilun)mstar_tot
     read(ilun)mstar_lost
     read(ilun)nsink     
     if(ncpu2.ne.ncpu.or.ndim2.ne.ndim.or.npart2.gt.npartmax)then
        write(*,*)'File part.tmp not compatible'
        write(*,*)'Found   =',ncpu2,ndim2,npart2
        write(*,*)'Expected=',ncpu,ndim,npartmax
        call clean_stop
     end if
     ! Read position
     allocate(xdp(1:npart2))
     do idim=1,ndim
        read(ilun)xdp
        xp(1:npart2,idim)=xdp
     end do
     ! Read velocity
     do idim=1,ndim
        read(ilun)xdp
        vp(1:npart2,idim)=xdp
     end do
     ! Read mass
     read(ilun)xdp
     mp(1:npart2)=xdp
     deallocate(xdp)
     ! Read identity
     allocate(isp8(1:npart2))
     read(ilun)isp8
     idp(1:npart2)=isp8
     deallocate(isp8)
     ! Read level
     allocate(isp(1:npart2))
     read(ilun)isp
     levelp(1:npart2)=isp
     deallocate(isp)
#ifdef OUTPUT_PARTICLE_POTENTIAL
      read(ilun)
#endif
     
     
     if(star.or.sink)then
        ! Read birth epoch
        allocate(xdp(1:npart2))
        read(ilun)xdp
        tp(1:npart2)=xdp
        if(convert_birth_times) then
           do i = 1, npart2 ! Convert birth time to proper for RT postpr.
              call getProperTime(tp(i),tp(i))
           enddo
        endif
        if(metal)then
           ! Read metallicity
           read(ilun)xdp
           zp(1:npart2)=xdp
        end if

        read(ilun)xdp
        tpp(1:npart2)=xdp
        read(ilun)xdp
        mp0(1:npart2)=xdp
        read(ilun)xdp
        indtab(1:npart2)=xdp
        deallocate(xdp)
     end if
     close(ilun)

     !determine NDM
     nDMloc=0
     if(star .or. sink) then
        do i = 1, npart2
           if(idp(i)>0.and.tp(i)==0.) then 
              nDMloc=nDMloc+1
           endif
        enddo
      else 
        nDMloc=npart2
      endif
#ifndef LONGINT
     call MPI_ALLREDUCE(nDMloc,nDM,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
#else
     tmp_long=nDMloc
     call MPI_ALLREDUCE(tmp_long,nDM,1,MPI_INTEGER8,MPI_SUM,MPI_COMM_WORLD,info)
#endif


     ! Send the token      
#ifndef WITHOUTMPI
     if(IOGROUPSIZE>0) then
        if(mod(myid,IOGROUPSIZE)/=0 .and.(myid.lt.ncpu))then
           dummy_io=1
           call MPI_SEND(dummy_io,1,MPI_INTEGER,myid-1+1,tagg, &
                & MPI_COMM_WORLD,info2)
        end if
     endif
#endif


     if(debug)write(*,*)'part.tmp read for processor ',myid
     npart=npart2

  else     

     filetype_loc=filetype
     if(.not. cosmo)filetype_loc='ascii'

     select case (filetype_loc)

     case ('grafic')

        !----------------------------------------------------
        ! Reading initial conditions GRAFIC2 multigrid arrays  
        !----------------------------------------------------
        ipart=0
        ! Loop over initial condition levels
        do ilevel=levelmin,nlevelmax
           
           if(initfile(ilevel)==' ')cycle
           
           ! Mesh size at level ilevel in coarse cell units
           dx=0.5D0**ilevel
           
           ! Set position of cell centers relative to grid center
           do ind=1,twotondim
              iz=(ind-1)/4
              iy=(ind-1-4*iz)/2
              ix=(ind-1-2*iy-4*iz)
              if(ndim>0)xc(ind,1)=(dble(ix)-0.5D0)*dx
              if(ndim>1)xc(ind,2)=(dble(iy)-0.5D0)*dx
              if(ndim>2)xc(ind,3)=(dble(iz)-0.5D0)*dx
           end do
           
           !--------------------------------------------------------------
           ! First step: compute level boundaries and particle positions
           !--------------------------------------------------------------
           i1_min=n1(ilevel)+1; i1_max=0
           i2_min=n2(ilevel)+1; i2_max=0
           i3_min=n3(ilevel)+1; i3_max=0
           ipart_old=ipart
           
           ! Loop over grids by vector sweeps
           ncache=active(ilevel)%ngrid
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
                 do i=1,ngrid
                    xx1=xg(ind_grid(i),1)+xc(ind,1)
                    xx1=(xx1*(dxini(ilevel)/dx)-xoff1(ilevel))/dxini(ilevel)
                    xx2=xg(ind_grid(i),2)+xc(ind,2)
                    xx2=(xx2*(dxini(ilevel)/dx)-xoff2(ilevel))/dxini(ilevel)
                    xx3=xg(ind_grid(i),3)+xc(ind,3)
                    xx3=(xx3*(dxini(ilevel)/dx)-xoff3(ilevel))/dxini(ilevel)
                    i1_min=MIN(i1_min,int(xx1)+1)
                    i1_max=MAX(i1_max,int(xx1)+1)
                    i2_min=MIN(i2_min,int(xx2)+1)
                    i2_max=MAX(i2_max,int(xx2)+1)
                    i3_min=MIN(i3_min,int(xx3)+1)
                    i3_max=MAX(i3_max,int(xx3)+1)
                    keep_part=son(ind_cell(i))==0
                    if(keep_part)then
                       ipart=ipart+1
                       if(ipart>npartmax)then
                          write(*,*)'Maximum number of particles incorrect'
                          write(*,*)'npartmax should be greater than',ipart
                          call clean_stop
                       endif
                       if(ndim>0)xp(ipart,1)=xg(ind_grid(i),1)+xc(ind,1)
                       if(ndim>1)xp(ipart,2)=xg(ind_grid(i),2)+xc(ind,2)
                       if(ndim>2)xp(ipart,3)=xg(ind_grid(i),3)+xc(ind,3)
                       mp(ipart)=0.5d0**(3*ilevel)*(1.0d0-omega_b/omega_m)
                    end if
                 end do
              end do
              ! End loop over cells
           end do
           ! End loop over grids
           
           ! Check that all grids are within initial condition region
           error=.false.
           if(active(ilevel)%ngrid>0)then
              if(i1_min<1.or.i1_max>n1(ilevel))error=.true.
              if(i2_min<1.or.i2_max>n2(ilevel))error=.true.
              if(i3_min<1.or.i3_max>n3(ilevel))error=.true.
           end if
           if(error) then
              write(*,*)'Some grid are outside initial conditions sub-volume'
              write(*,*)'for ilevel=',ilevel
              write(*,*)i1_min,i1_max
              write(*,*)i2_min,i2_max
              write(*,*)i3_min,i3_max
              write(*,*)n1(ilevel),n2(ilevel),n3(ilevel)
!jhshin1
!              call clean_stop
!jhshin2
           end if
           if(debug)then
              write(*,*)myid,i1_min,i1_max,i2_min,i2_max,i3_min,i3_max
           endif
           
           !---------------------------------------------------------------------
           ! Second step: read initial condition file and set particle velocities
           !---------------------------------------------------------------------
           ! Allocate initial conditions array
           if(active(ilevel)%ngrid>0)then
              allocate(init_array(i1_min:i1_max,i2_min:i2_max,i3_min:i3_max))
              allocate(init_array_x(i1_min:i1_max,i2_min:i2_max,i3_min:i3_max))
              init_array=0d0
              init_array_x=0d0
           end if
           allocate(init_plane(1:n1(ilevel),1:n2(ilevel)))
           allocate(init_plane_x(1:n1(ilevel),1:n2(ilevel)))
           
           ! Loop over input variables
           do idim=1,ndim
              
              ! Read dark matter initial displacement field
              if(multiple)then
                 call title(myid,nchar)
                 if(idim==1)filename=TRIM(initfile(ilevel))//'/dir_velcx/ic_velcx.'//TRIM(nchar)
                 if(idim==2)filename=TRIM(initfile(ilevel))//'/dir_velcy/ic_velcy.'//TRIM(nchar)
                 if(idim==3)filename=TRIM(initfile(ilevel))//'/dir_velcz/ic_velcz.'//TRIM(nchar)
              else
                 if(idim==1)filename=TRIM(initfile(ilevel))//'/ic_velcx'
                 if(idim==2)filename=TRIM(initfile(ilevel))//'/ic_velcy'
                 if(idim==3)filename=TRIM(initfile(ilevel))//'/ic_velcz'

                 if(idim==1)filename_x=TRIM(initfile(ilevel))//'/ic_poscx'
                 if(idim==2)filename_x=TRIM(initfile(ilevel))//'/ic_poscy'
                 if(idim==3)filename_x=TRIM(initfile(ilevel))//'/ic_poscz'

                 INQUIRE(file=filename_x,exist=ok)
                 if(.not.ok)then
                    read_pos = .false.
                 else
                    read_pos = .true.
                    if(myid==1)write(*,*)'Reading file '//TRIM(filename_x)
                 end if

              endif

              if(myid==1)write(*,*)'Reading file '//TRIM(filename)
                               
              if(multiple)then
                 ilun=myid+10
                 ! Wait for the token                                                                                                                                                        
#ifndef WITHOUTMPI
                 if(IOGROUPSIZE>0) then
                    if (mod(myid-1,IOGROUPSIZE)/=0) then
                       call MPI_RECV(dummy_io,1,MPI_INTEGER,myid-1-1,tagg2,&
                            & MPI_COMM_WORLD,MPI_STATUS_IGNORE,info2)
                    end if
                 endif
#endif
                 open(ilun,file=filename,form='unformatted')
                 rewind ilun
                 read(ilun) ! skip first line
                 do i3=1,n3(ilevel)
                    read(ilun)((init_plane(i1,i2),i1=1,n1(ilevel)),i2=1,n2(ilevel))
                    if(active(ilevel)%ngrid>0)then
                       if(i3.ge.i3_min.and.i3.le.i3_max)then
                          init_array(i1_min:i1_max,i2_min:i2_max,i3) = &
                               & init_plane(i1_min:i1_max,i2_min:i2_max)
                       end if
                    endif
                 end do
                 close(ilun)
                 ! Send the token                                                                                                                                                            
#ifndef WITHOUTMPI
                 if(IOGROUPSIZE>0) then
                    if(mod(myid,IOGROUPSIZE)/=0 .and.(myid.lt.ncpu))then
                       dummy_io=1
                       call MPI_SEND(dummy_io,1,MPI_INTEGER,myid-1+1,tagg2, &
                            & MPI_COMM_WORLD,info2)
                    end if
                 endif
#endif

              else
                 ! Parallel IC reading with stream access: jump directly to needed planes
                 if(active(ilevel)%ngrid>0)then
                    open(10,file=filename,access='stream',form='unformatted',status='old')
                    hdr_bytes = 52_8   ! 4 + 44 + 4 (GRAFIC2 header record)
                    plane_bytes = int(n1(ilevel),8) * int(n2(ilevel),8) * 4_8 + 8_8
                    do i3=i3_min,i3_max
                       byte_pos = hdr_bytes + int(i3-1,8)*plane_bytes + 5_8
                       read(10, pos=byte_pos)((init_plane(i1,i2),i1=1,n1(ilevel)),i2=1,n2(ilevel))
                       init_array(i1_min:i1_max,i2_min:i2_max,i3) = &
                            & init_plane(i1_min:i1_max,i2_min:i2_max)
                    end do
                    close(10)
                 endif

                 if(read_pos) then
                    if(active(ilevel)%ngrid>0)then
                       open(10,file=filename_x,access='stream',form='unformatted',status='old')
                       hdr_bytes = 52_8
                       plane_bytes = int(n1(ilevel),8) * int(n2(ilevel),8) * 4_8 + 8_8
                       do i3=i3_min,i3_max
                          byte_pos = hdr_bytes + int(i3-1,8)*plane_bytes + 5_8
                          read(10, pos=byte_pos)((init_plane_x(i1,i2),i1=1,n1(ilevel)),i2=1,n2(ilevel))
                          init_array_x(i1_min:i1_max,i2_min:i2_max,i3) = &
                               & init_plane_x(i1_min:i1_max,i2_min:i2_max)
                       end do
                       close(10)
                    endif
                 end if

              endif
              
              if(active(ilevel)%ngrid>0)then
                 ! Rescale initial displacement field to code units
                 init_array=dfact(ilevel)*dx/dxini(ilevel)*init_array/vfact(ilevel)
                 if(read_pos)then
                    init_array_x = init_array_x/boxlen_ini
                 endif
                 ! Loop over grids by vector sweeps
                 ipart=ipart_old
                 ncache=active(ilevel)%ngrid
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
                       do i=1,ngrid
                          xx1=xg(ind_grid(i),1)+xc(ind,1)
                          xx1=(xx1*(dxini(ilevel)/dx)-xoff1(ilevel))/dxini(ilevel)
                          xx2=xg(ind_grid(i),2)+xc(ind,2)
                          xx2=(xx2*(dxini(ilevel)/dx)-xoff2(ilevel))/dxini(ilevel)
                          xx3=xg(ind_grid(i),3)+xc(ind,3)
                          xx3=(xx3*(dxini(ilevel)/dx)-xoff3(ilevel))/dxini(ilevel)
                          i1=int(xx1)+1
                          i1=int(xx1)+1
                          i2=int(xx2)+1
                          i2=int(xx2)+1
                          i3=int(xx3)+1
                          i3=int(xx3)+1
                          keep_part=son(ind_cell(i))==0
                          if(keep_part)then
                             ipart=ipart+1
                             vp(ipart,idim)=init_array(i1,i2,i3)
                             if(.not. read_pos)then
                                dispmax=max(dispmax,abs(init_array(i1,i2,i3)/dx))
                             else
                                xp(ipart,idim)=xg(ind_grid(i),idim)+xc(ind,idim)+init_array_x(i1,i2,i3)
                                dispmax=max(dispmax,abs(init_array_x(i1,i2,i3)/dx))
                             endif
                          end if
                       end do
                    end do
                    ! End loop over cells
                 end do
                 ! End loop over grids
              endif

           end do
           ! End loop over input variables
           
           ! Deallocate initial conditions array
           if(active(ilevel)%ngrid>0)then
              deallocate(init_array,init_array_x)
           end if
           deallocate(init_plane,init_plane_x)
           
           if(debug)write(*,*)'npart=',ipart,'/',npartmax,' forPE=',myid,dispmax
      
           
        end do
        ! End loop over levels
        
        ! Initial particle number
        npart=ipart

        !has to be added because we don't track clouds and can't retrieve DM numbers from the total
#ifndef LONGINT
        call MPI_ALLREDUCE(npart,nDM,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
#else
        tmp_long=npart
        call MPI_ALLREDUCE(tmp_long,nDM,1,MPI_INTEGER8,MPI_SUM,MPI_COMM_WORLD,info)
#endif


        ! Move particle according to Zeldovich approximation
        if(.not. read_pos)then
           xp(1:npart,1:ndim)=xp(1:npart,1:ndim)+vp(1:npart,1:ndim)
        endif

        ! Scale displacement to velocity
        vp(1:npart,1:ndim)=vfact(1)*vp(1:npart,1:ndim)
        
        ! Periodic box
        do ipart=1,npart
#if NDIM>0
           if(xp(ipart,1)<  0.0d0  )xp(ipart,1)=xp(ipart,1)+dble(nx)
           if(xp(ipart,1)>=dble(nx))xp(ipart,1)=xp(ipart,1)-dble(nx)
#endif
#if NDIM>1
           if(xp(ipart,2)<  0.0d0  )xp(ipart,2)=xp(ipart,2)+dble(ny)
           if(xp(ipart,2)>=dble(ny))xp(ipart,2)=xp(ipart,2)-dble(ny)
#endif
#if NDIM>2
           if(xp(ipart,3)<  0.0d0  )xp(ipart,3)=xp(ipart,3)+dble(nz)
           if(xp(ipart,3)>=dble(nz))xp(ipart,3)=xp(ipart,3)-dble(nz)
#endif
        end do
        
#ifndef WITHOUTMPI        

        ! Compute particle Hilbert ordering
        sendbuf=0
        do ipart=1,npart
           xx(1,1:3)=xp(ipart,1:3)
           xx_dp(1,1:3)=xx(1,1:3)
           call cmp_cpumap(xx_dp,cc,1)
           if(cc(1).ne.myid)sendbuf(cc(1))=sendbuf(cc(1))+1
        end do
           
        ! Allocate communication buffer in emission
        do icpu=1,ncpu
           ncache=sendbuf(icpu)
           if(ncache>0)then
              allocate(emission(icpu,1)%up(1:ncache,1:twondim+1))
           end if
        end do

        ! Fill communicators
        jpart=0
        sendbuf=0
        do ipart=1,npart
           xx(1,1:3)=xp(ipart,1:3)
           xx_dp(1,1:3)=xx(1,1:3)
           call cmp_cpumap(xx_dp,cc,1)
           if(cc(1).ne.myid)then
              icpu=cc(1)
              sendbuf(icpu)=sendbuf(icpu)+1
              ibuf=sendbuf(icpu)
              emission(icpu,1)%up(ibuf,1)=xp(ipart,1)
              emission(icpu,1)%up(ibuf,2)=xp(ipart,2)
              emission(icpu,1)%up(ibuf,3)=xp(ipart,3)
              emission(icpu,1)%up(ibuf,4)=vp(ipart,1)
              emission(icpu,1)%up(ibuf,5)=vp(ipart,2)
              emission(icpu,1)%up(ibuf,6)=vp(ipart,3)
              emission(icpu,1)%up(ibuf,7)=mp(ipart)
           else
              jpart=jpart+1
              xp(jpart,1:3)=xp(ipart,1:3)
              vp(jpart,1:3)=vp(ipart,1:3)
              mp(jpart)    =mp(ipart)
           endif
        end do

        ! Communicate virtual particle number to parent cpu
        if(ordering=='ksection') then
           ntotal_ksec = 0
           do icpu = 1, ncpu
              if(sendbuf(icpu) > 0) ntotal_ksec = ntotal_ksec + 1
           end do
           allocate(sbuf_ksec(1:2, 1:max(ntotal_ksec,1)))
           allocate(dcpu_ksec(1:max(ntotal_ksec,1)))
           idx_ksec = 0
           do icpu = 1, ncpu
              if(sendbuf(icpu) > 0) then
                 idx_ksec = idx_ksec + 1
                 dcpu_ksec(idx_ksec) = icpu
                 sbuf_ksec(1, idx_ksec) = dble(myid)
                 sbuf_ksec(2, idx_ksec) = dble(sendbuf(icpu))
              end if
           end do
           call ksection_exchange_dp(sbuf_ksec, ntotal_ksec, dcpu_ksec, 2, rbuf_ksec, nrecv_ksec)
           recvbuf = 0
           do idx_ksec = 1, nrecv_ksec
              recvbuf(nint(rbuf_ksec(1, idx_ksec))) = nint(rbuf_ksec(2, idx_ksec))
           end do
           deallocate(sbuf_ksec, dcpu_ksec, rbuf_ksec)
        else
           call MPI_ALLTOALL(sendbuf,1,MPI_INTEGER,recvbuf,1,MPI_INTEGER,MPI_COMM_WORLD,info)
        end if

        ! Compute total number of newly created particles
        npart_new=0
        do icpu=1,ncpu
           npart_new=npart_new+recvbuf(icpu)
        end do

        if(jpart+npart_new.gt.npartmax)then
           write(*,*)'No more free memory for particles'
           write(*,*)'Increase npartmax'
           write(*,*)myid
           write(*,*)jpart,npart_new
           write(*,*)bound_key
           call MPI_ABORT(MPI_COMM_WORLD,1,info)
        end if

        ! Allocate communication buffer in reception
        do icpu=1,ncpu
           ncache=recvbuf(icpu)
           if(ncache>0)then
              allocate(reception(icpu,1)%up(1:ncache,1:twondim+1))
           end if
        end do

        ! Receive particles
        countrecv=0
        do icpu=1,ncpu
           ncache=recvbuf(icpu)
           if(ncache>0)then
              buf_count=ncache*(twondim+1)
              countrecv=countrecv+1
              call MPI_IRECV(reception(icpu,1)%up,buf_count, &
                   & MPI_DOUBLE_PRECISION,icpu-1,&
                   & tagu,MPI_COMM_WORLD,reqrecv(countrecv),info)
           end if
        end do

        ! Send particles
        countsend=0
        do icpu=1,ncpu
           ncache=sendbuf(icpu)
           if(ncache>0)then
              buf_count=ncache*(twondim+1)
              countsend=countsend+1
              call MPI_ISEND(emission(icpu,1)%up,buf_count, &
                   & MPI_DOUBLE_PRECISION,icpu-1,&
                   & tagu,MPI_COMM_WORLD,reqsend(countsend),info)
           end if
        end do

        ! Wait for full completion of receives
        call MPI_WAITALL(countrecv,reqrecv,statuses,info)
        
        ! Wait for full completion of sends
        call MPI_WAITALL(countsend,reqsend,statuses,info)

        ! Create new particles
        do icpu=1,ncpu
           do ibuf=1,recvbuf(icpu)
              jpart=jpart+1
              xp(jpart,1)=reception(icpu,1)%up(ibuf,1)
              xp(jpart,2)=reception(icpu,1)%up(ibuf,2)
              xp(jpart,3)=reception(icpu,1)%up(ibuf,3)
              vp(jpart,1)=reception(icpu,1)%up(ibuf,4)
              vp(jpart,2)=reception(icpu,1)%up(ibuf,5)
              vp(jpart,3)=reception(icpu,1)%up(ibuf,6)
              mp(jpart)  =reception(icpu,1)%up(ibuf,7)
           end do
        end do
        
        ! Erase old particles
        do ipart=jpart+1,npart
           xp(ipart,1)=0d0
           xp(ipart,2)=0d0
           xp(ipart,3)=0d0
           vp(ipart,1)=0d0
           vp(ipart,2)=0d0
           vp(ipart,3)=0d0
           mp(ipart)  =0d0
        end do
        npart=jpart

        ! Deallocate communicators
        do icpu=1,ncpu
           if(sendbuf(icpu)>0)deallocate(emission(icpu,1)%up)
           if(recvbuf(icpu)>0)deallocate(reception(icpu,1)%up)
        end do

        write(*,*)'npart=',ipart,'/',npartmax,' for PE=',myid
#endif
        

        ! Compute particle initial level
        do ipart=1,npart
           levelp(ipart)=levelmin
        end do

        ! Compute particle initial age and metallicity
        if(star.or.sink)then
           do ipart=1,npart
              tp(ipart)=0d0
              if(metal)then
                 zp(ipart)=0d0
              end if
              tpp(ipart)=0d0
              mp0(ipart)=0d0
              indtab(ipart)=0d0
!              cep(ipart,:nelt)=0d0
           end do
        end if

        ! Compute particle initial identity
        npart_cpu=0; npart_all=0
        npart_cpu(myid)=npart
#ifndef WITHOUTMPI
#ifndef LONGINT
        call MPI_ALLREDUCE(npart_cpu,npart_all,ncpu,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
#else
        call MPI_ALLREDUCE(npart_cpu,npart_all,ncpu,MPI_INTEGER8,MPI_SUM,MPI_COMM_WORLD,info)
#endif
        npart_cpu(1)=npart_all(1)
#endif
        do icpu=2,ncpu
           npart_cpu(icpu)=npart_cpu(icpu-1)+npart_all(icpu)
        end do
        if(myid==1)then
           do ipart=1,npart
              idp(ipart)=ipart
           end do
        else
           do ipart=1,npart
              idp(ipart)=npart_cpu(myid-1)+ipart
           end do
        end if

     case ('ascii')

        ! Local particle count
        ipart=0

        if(TRIM(initfile(levelmin)).NE.' ')then

        filename=TRIM(initfile(levelmin))//'/ic_part'
        if(myid==1)then
           open(10,file=filename,form='formatted')
           indglob=0
        end if
        eof=.false.

        do while (.not.eof)
           xx=0.0
           if(myid==1)then
              jpart=0
              do i=1,nvector
                 read(10,*,end=100)xx1,xx2,xx3,vv1,vv2,vv3,mm1
                 jpart=jpart+1
                 indglob=indglob+1
                 xx(i,1)=xx1+boxlen/2.0
                 xx(i,2)=xx2+boxlen/2.0
                 xx(i,3)=xx3+boxlen/2.0
                 vv(i,1)=vv1
                 vv(i,2)=vv2
                 vv(i,3)=vv3
                 mm(i  )=mm1
                 ii(i  )=indglob
              end do
100           continue
              if(jpart<nvector)eof=.true.
           endif
           buf_count=nvector*3
#ifndef WITHOUTMPI
           call MPI_BCAST(xx,buf_count,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,info)
           call MPI_BCAST(vv,buf_count,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,info)
           call MPI_BCAST(mm,nvector  ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,info)
           call MPI_BCAST(ii,nvector  ,MPI_INTEGER         ,0,MPI_COMM_WORLD,info)
           call MPI_BCAST(eof,1       ,MPI_LOGICAL         ,0,MPI_COMM_WORLD,info)
           call MPI_BCAST(jpart,1     ,MPI_INTEGER         ,0,MPI_COMM_WORLD,info)
           call cmp_cpumap(xx,cc,jpart)
#endif


           do i=1,jpart
#ifndef WITHOUTMPI
              if(cc(i)==myid)then
#endif
                 ipart=ipart+1
                 if(ipart>npartmax)then
                    write(*,*)'Maximum number of particles incorrect'
                    write(*,*)'npartmax should be greater than',ipart
                    call clean_stop
                 endif
                 xp(ipart,1:3)=xx(i,1:3)
                 vp(ipart,1:3)=vv(i,1:3)
                 mp(ipart)    =mm(i)
                 levelp(ipart)=levelmin
                 idp(ipart)   =ii(i)
#ifndef WITHOUTMPI
              endif
#endif
           enddo

        end do

        if(myid==1)close(10)

        end if
        npart=ipart
        
        !has to be added because we don't track clouds and can't retrieve DM numbers from the total
#ifndef LONGINT
        call MPI_ALLREDUCE(npart,nDM,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
#else
        tmp_long=npart
        call MPI_ALLREDUCE(tmp_long,nDM,1,MPI_INTEGER8,MPI_SUM,MPI_COMM_WORLD,info)
#endif

        ! Compute total number of particle
        npart_cpu=0; npart_all=0
        npart_cpu(myid)=npart
#ifndef WITHOUTMPI
#ifndef LONGINT
        call MPI_ALLREDUCE(npart_cpu,npart_all,ncpu,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
#else
        call MPI_ALLREDUCE(npart_cpu,npart_all,ncpu,MPI_INTEGER8,MPI_SUM,MPI_COMM_WORLD,info)
#endif
        npart_cpu(1)=npart_all(1)

#endif
        do icpu=2,ncpu
           npart_cpu(icpu)=npart_cpu(icpu-1)+npart_all(icpu)
        end do
        if(debug)write(*,*)'npart=',npart,'/',npart_cpu(ncpu)

     case ('gadget')
        call load_gadget

     case DEFAULT
        write(*,*) 'Unsupported format file ' // filetype
        call clean_stop

     end select
  end if

  call getmem(real_mem)
  call MPI_ALLREDUCE(real_mem,real_mem_tot,1,MPI_REAL,MPI_MAX,MPI_COMM_WORLD,info)
  if(myid==1) then
      call writemem(real_mem_tot)
  endif

  if(sink)call init_sink

end subroutine init_part
#define TIME_START(cs) call SYSTEM_CLOCK(COUNT=cs)
#define TIME_END(ce) call SYSTEM_CLOCK(COUNT=ce)
#define TIME_SPENT(cs,ce,cr) REAL((ce-cs)/cr)
subroutine load_gadget
  use amr_commons
  use pm_commons
  use gadgetreadfilemod
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif

  logical::ok
  TYPE(gadgetheadertype) :: gadgetheader
  integer::numfiles
  integer::ifile
  real(dp),dimension(1:nvector,1:3)::xx_dp
  real, dimension(:, :), allocatable:: pos, vel
  real(dp)::massparticles
  integer(kind=8)::allparticles
  integer(i8b) :: tmp_long
  integer(i8b), dimension(:), allocatable:: ids  
  integer::nparticles, arraysize
  integer::i, icpu, ipart, info, np, start
  integer(i8b),dimension(1:ncpu)::npart_cpu,npart_all
  character(LEN=256)::filename
  integer ,dimension(1:nvector)::cc
  integer :: clock_start, clock_end, clock_rate
  integer :: mpi_cs, mpi_ce
  real:: gadgetvfact
  ! Local particle count
  ipart=0
  call SYSTEM_CLOCK(COUNT_RATE=clock_rate)

  if(TRIM(initfile(levelmin)).NE.' ')then
     filename=TRIM(initfile(levelmin))
     ! read first header to get information
     call gadgetreadheader(filename, 0, gadgetheader, ok)
     if(.not.ok) call clean_stop
     numfiles = gadgetheader%numfiles
     gadgetvfact = SQRT(aexp) / gadgetheader%boxsize * aexp / 100.
#ifndef LONGINT
     allparticles=int(gadgetheader%nparttotal(2),kind=8)
#else
     allparticles=int(gadgetheader%nparttotal(2),kind=8) &
          & +int(gadgetheader%totalhighword(2),kind=8)*4294967296 !2^32
#endif
     massparticles=1d0/dble(allparticles)
     do ifile=0,numfiles-1
        call gadgetreadheader(filename, ifile, gadgetheader, ok)
        nparticles = gadgetheader%npart(2)
        allocate(pos(3,nparticles))
        allocate(vel(3,nparticles))
        allocate(ids(nparticles))
        TIME_START(clock_start)
        call gadgetreadfile(filename,ifile,gadgetheader, pos, vel, ids)
        TIME_END(clock_end)
        if(debug) write(*,*) myid, ':Read ', nparticles, ' from gadget file ', ifile, ' in ', &
        TIME_SPENT(clock_start, clock_end, clock_rate)
        start = 1
        TIME_START(clock_start)
#ifndef WITHOUTMPI
        do i=1,nparticles
           xx_dp(1,1) = pos(1,i)/gadgetheader%boxsize
           xx_dp(1,2) = pos(2,i)/gadgetheader%boxsize
           xx_dp(1,3) = pos(3,i)/gadgetheader%boxsize
           call cmp_cpumap(xx_dp,cc,1)
           if(cc(1)==myid)then
#endif
              ipart=ipart+1
              if (ipart .ge. size(mp)) then
                 write(*,*) "For ", myid, ipart, " exceeds ", size(mp)
                 call clean_stop
              end if
              xp(ipart,1:3)=xx_dp(1,1:3)
              vp(ipart,1)  =vel(1, i) * gadgetvfact
              vp(ipart,2)  =vel(2, i) * gadgetvfact
              vp(ipart,3)  =vel(3, i) * gadgetvfact
              mp(ipart)    = massparticles
              levelp(ipart)=levelmin
              idp(ipart)   =ids(i)
#ifndef WITHOUTMPI
            endif
        enddo
        TIME_END(clock_end)
        if(debug) write(*,*) myid, ':Processed ', nparticles, ' in ',&
             &  TIME_SPENT(clock_start, clock_end, clock_rate), " ipart now ", ipart
#endif
        deallocate(pos,vel,ids)
     end do

  end if
  npart=ipart

  !has to be added because we don't track clouds and can't retrieve DM numbers from the total
#ifndef LONGINT
  call MPI_ALLREDUCE(npart,nDM,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
#else
  tmp_long=npart
  call MPI_ALLREDUCE(tmp_long,nDM,1,MPI_INTEGER8,MPI_SUM,MPI_COMM_WORLD,info)
#endif

  ! Compute total number of particleclock_rate
  npart_cpu=0; npart_all=0
  npart_cpu(myid)=npart
#ifndef WITHOUTMPI
#ifndef LONGINT
  call MPI_ALLREDUCE(npart_cpu,npart_all,ncpu,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
#else
  call MPI_ALLREDUCE(npart_cpu,npart_all,ncpu,MPI_INTEGER8,MPI_SUM,MPI_COMM_WORLD,info)
#endif
  npart_cpu(1)=npart_all(1)
#endif
  do icpu=2,ncpu
     npart_cpu(icpu)=npart_cpu(icpu-1)+npart_all(icpu)
  end do
  write(*,*)'npart=',npart,'/',npartmax

end subroutine load_gadget
