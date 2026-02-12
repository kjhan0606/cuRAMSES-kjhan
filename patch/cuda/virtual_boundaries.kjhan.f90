!################################################################ 
!################################################################
!################################################################
!################################################################
subroutine authorize_coarse
  use amr_commons
  implicit none
  !----------------------------------------------------------------------
  ! This routine authorizes all base cells for refinement.
  ! This duplicates all base grids over all cpu's.
  !----------------------------------------------------------------------
  integer(i8b)::nxny,i,j,k,ind

  if(verbose)write(*,*)'  Entering authorize_coarse'
  ! Constants
  nxny=nx*ny
  ! Initialize flag2(0) to zero
  flag2(0)=0
  ! Duplicate full domain over cpus
!$omp parallel do private(i,j,k,ind)
  do k=0,nz-1
  do j=0,ny-1
  do i=0,nx-1
     ind=1+i+j*nx+k*nxny
     flag2(ind)=1
  end do
  end do
  end do
  
end subroutine authorize_coarse
!################################################################
!################################################################
!################################################################
!################################################################
subroutine sub2_authorize_fine(ibound,ismooth,icpu,ilevel, igrid,ngrid)
  use amr_commons
  use bisection
  implicit none
  integer::ilevel
  ! -------------------------------------------------------------------
  ! This routine computes the authorization map (flag2) for level ilevel.
  ! All myid cells are first marked for authorization.
  ! All virtual cells that intersect the local ordering domain are
  ! also marked for authorization. Finally, the routine performs
  ! a dilatation of the authorization map of one cell width.
  ! Array flag1 for virtual cells is used as temporary work space.
  ! -------------------------------------------------------------------
  integer::ismooth,ibound,ngrid,i,ncache,iskip,igrid,ind,icpu
  integer::ix,iy,iz,idim,nx_loc,isub
  integer,dimension(1:3)::n_nbor
  integer,dimension(1:nvector)::ind_grid,ind_cell
  integer,dimension(1:nvector,0:twondim)::igridn
  real(dp)::dx,dx_loc,scale
  real(dp),dimension(1:3)::skip_loc
  real(dp),dimension(1:twotondim,1:3)::xc
  real(dp),dimension(1:nvector,1:ndim)::xx
  real(qdp),dimension(1:nvector)::order_min,order_max
  logical::test
  real(dp),dimension(1:ndim)::xmin,xmax
  n_nbor(1:3)=(/1,2,3/)
  do i=1,ngrid
     ind_grid(i)=reception(icpu,ilevel)%igrid(igrid+i-1)
  end do
  do ind=1,twotondim
     iskip=ncoarse+(ind-1)*ngridmax
     do i=1,ngrid
        ind_cell(i)=iskip+ind_grid(i)
     end do
     do i=1,ngrid
        flag1(ind_cell(i))=0
     end do
  end do
  call getnborgrids(ind_grid,igridn,ngrid)
  do ind=1,twotondim
     call count_nbors2(igridn,ind,n_nbor(ismooth),ngrid)
  end do
  do ind=1,twotondim
     iskip=ncoarse+(ind-1)*ngridmax
     do i=1,ngrid
        ind_cell(i)=iskip+ind_grid(i)
     end do
     do i=1,ngrid
        if(flag1(ind_cell(i))==1)flag2(ind_cell(i))=1
     end do
  end do
end subroutine sub2_authorize_fine
!################################################################
!################################################################
!################################################################
!################################################################
subroutine sub1_authorize_fine(icpu,ilevel, igrid,ngrid, dx_loc,xc,skip_loc, scale)
  use amr_commons
  use bisection
  implicit none
  integer::ilevel
  ! -------------------------------------------------------------------
  ! This routine computes the authorization map (flag2) for level ilevel.
  ! All myid cells are first marked for authorization.
  ! All virtual cells that intersect the local ordering domain are
  ! also marked for authorization. Finally, the routine performs
  ! a dilatation of the authorization map of one cell width.
  ! Array flag1 for virtual cells is used as temporary work space.
  ! -------------------------------------------------------------------
  integer::ismooth,ibound,ngrid,i,ncache,iskip,igrid,ind,icpu
  integer::ix,iy,iz,idim,nx_loc,isub
  integer,dimension(1:3)::n_nbor
  integer,dimension(1:nvector)::ind_grid,ind_cell
  integer,dimension(1:nvector,0:twondim)::igridn
  real(dp)::dx,dx_loc,scale
  real(dp),dimension(1:3)::skip_loc
  real(dp),dimension(1:twotondim,1:3)::xc
  real(dp),dimension(1:nvector,1:ndim)::xx
  real(qdp),dimension(1:nvector)::order_min,order_max
  logical::test
  real(dp),dimension(1:ndim)::xmin,xmax

  do i=1,ngrid
     ind_grid(i)=reception(icpu,ilevel)%igrid(igrid+i-1)
  end do
  ! Loop over cells
  do ind=1,twotondim
     ! Gather cell indices
     iskip=ncoarse+(ind-1)*ngridmax
     do i=1,ngrid
        ind_cell(i)=iskip+ind_grid(i)
     end do
     ! Gather cell centre positions
     do idim=1,ndim
        do i=1,ngrid
           xx(i,idim)=xg(ind_grid(i),idim)+xc(ind,idim)
        end do
     end do
     ! Rescale position from code units to user units
     do idim=1,ndim
        do i=1,ngrid
           xx(i,idim)=(xx(i,idim)-skip_loc(idim))*scale
        end do
     end do
     ! Reset flag2
     do i=1,ngrid
        flag2(ind_cell(i))=0
     end do

     if (.not.use_cpubox_decomp) then
        ! Compute minimum and maximum ordering key
        call cmp_minmaxorder(xx,order_min,order_max,dx_loc,ngrid)
        ! Determine if cell is authorized
        do isub=1,overload
           do i=1,ngrid
              if(    order_max(i)>bound_key(myid-1+(isub-1)*ncpu).and.&
                   & order_min(i)<bound_key(myid  +(isub-1)*ncpu) )then
                 flag2(ind_cell(i))=1
              endif
           end do
        end do
     else ! recursive bisection method                                                          
         do i=1,ngrid
            ! Test if cell overlaps the cpu                                                     
            test=.true.
            xmin=xx(i,:)-0.5*dx_loc
            xmax=xx(i,:)+0.5*dx_loc
            do idim=1,ndim
               ! This needs to be a >=, not a >, to precisely match the                         
               ! ordering/=case for refinement flagging                                         
               test=test .and. (bisec_cpubox_max(myid,idim).ge.xmin(idim) &
                                    .and. bisec_cpubox_min(myid,idim).le.xmax(idim))
            end do
            if(test) flag2(ind_cell(i))=1
         end do
     endif

     ! For load balancing operations
     if(balance)then
        if(.not.use_cpubox_decomp) then
           do isub=1,overload
              do i=1,ngrid
                 if(    order_max(i)>bound_key2(myid-1+(isub-1)*ncpu).and.&
                      & order_min(i)<bound_key2(myid  +(isub-1)*ncpu) )then
                    flag2(ind_cell(i))=1
                 endif
              end do
           end do
        else
           do i=1,ngrid
              ! Test if cell overlaps the cpu with new cpu map                                  
              test=.true.
              xmin=xx(i,:)-0.5*dx_loc
              xmax=xx(i,:)+0.5*dx_loc
              do idim=1,ndim
                 ! This needs to be a >=, not a >, to precisely match the                       
                 ! ordering/=case for refinement flagging                                       
                 test=test .and. (bisec_cpubox_max2(myid,idim).ge.xmin(idim) &
                      .and. bisec_cpubox_min2(myid,idim).le.xmax(idim))
              end do
              if(test) flag2(ind_cell(i))=1
           end do
        end if
        do i=1,ngrid
           if(cpu_map2(father(ind_grid(i)))==myid)then
              flag2(ind_cell(i))=1
           endif
        end do
     end if
  end do
  ! End loop over cells
end subroutine sub1_authorize_fine
!################################################################
!################################################################
!################################################################
!################################################################
subroutine authorize_fine(ilevel)
  use amr_commons
  use bisection
  implicit none
  integer::ilevel
  ! -------------------------------------------------------------------
  ! This routine computes the authorization map (flag2) for level ilevel.
  ! All myid cells are first marked for authorization.
  ! All virtual cells that intersect the local ordering domain are
  ! also marked for authorization. Finally, the routine performs
  ! a dilatation of the authorization map of one cell width.
  ! Array flag1 for virtual cells is used as temporary work space.
  ! -------------------------------------------------------------------
  integer::ismooth,ibound,ngrid,i,ncache,iskip,igrid,ind,icpu
  integer::ix,iy,iz,idim,nx_loc,isub,indgrid
  integer,dimension(1:3)::n_nbor
  integer,dimension(1:nvector)::ind_grid,ind_cell
  integer,dimension(1:nvector,0:twondim)::igridn
  real(dp)::dx,dx_loc,scale
  real(dp),dimension(1:3)::skip_loc
  real(dp),dimension(1:twotondim,1:3)::xc
  real(dp),dimension(1:nvector,1:ndim)::xx
  real(qdp),dimension(1:nvector)::order_min,order_max
  logical::test
  real(dp),dimension(1:ndim)::xmin,xmax

  if(ilevel==nlevelmax)return
  if(verbose)write(*,111)ilevel

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
  
  ! Scaling factor
  nx_loc=(icoarse_max-icoarse_min+1)
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  if(ndim>0)skip_loc(1)=dble(icoarse_min)
  if(ndim>1)skip_loc(2)=dble(jcoarse_min)
  if(ndim>2)skip_loc(3)=dble(kcoarse_min)
  scale=boxlen/dble(nx_loc)
  dx_loc=dx*scale

  ! Authorize all myid grids (needed for uploads)
  ncache=active(ilevel)%ngrid
  ! Loop over grids by vector sweeps

!$omp parallel do private(igrid,ngrid,i,ind_grid,ind,iskip,ind_cell)
  do igrid=1,ncache,nvector
     ! Gather nvector grids
     ngrid=MIN(nvector,ncache-igrid+1)
     do i=1,ngrid
        ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
     end do
     ! Loop over cells
     do ind=1,twotondim
        ! Gather cell indices
        iskip=ncoarse+(ind-1)*ngridmax
        do i=1,ngrid
           ind_cell(i)=iskip+ind_grid(i)
        end do
        do i=1,ngrid
           flag2(ind_cell(i))=1
        end do
     end do
     ! End loop over cells
  end do
  ! End loop over grids

  ! Authorize virtual cells that contains myid children cells
  do icpu=1,ncpu
     ncache=reception(icpu,ilevel)%ngrid
     ! Loop over grids by vector sweeps
!$omp parallel do private(igrid,ngrid)
     do igrid=1,ncache,nvector
        ! Gather nvector grids
        ngrid=MIN(nvector,ncache-igrid+1)
        call sub1_authorize_fine(icpu,ilevel, igrid,ngrid,dx_loc,xc,skip_loc, scale)
     enddo
  enddo
  ! End loop over cpus

  ! Apply dilatation operator over flag2 cells on virtual cells only
     
  flag2(0)=0
  ! Set flag2 to 0 for physical boundary grids
!$omp parallel do private(ibound,ind,iskip,i)
  do ibound=1,nboundary
  do ind=1,twotondim
     iskip=ncoarse+(ind-1)*ngridmax
     do i=1,boundary(ibound,ilevel)%ngrid
        flag2(boundary(ibound,ilevel)%igrid(i)+iskip)=0
     end do
  end do
  end do

  ! Loop over steps
  do ibound=1,nexpand_bound
  n_nbor(1:3)=(/1,2,3/)
  do ismooth=1,ndim
     ! Initialize flag1 to 0 in virtual cells
!$omp parallel do firstprivate(icpu,ncache,igrid,ngrid,i,ind_grid,ind,iskip,ind_cell)
     do icpu=1,ncpu
        ncache=reception(icpu,ilevel)%ngrid
        do igrid=1,ncache,nvector
           ngrid=MIN(nvector,ncache-igrid+1)
           do i=1,ngrid
              ind_grid(i)=reception(icpu,ilevel)%igrid(igrid+i-1)
           end do
           do ind=1,twotondim
              iskip=ncoarse+(ind-1)*ngridmax
              do i=1,ngrid
                 ind_cell(i)=iskip+ind_grid(i)
              end do
              do i=1,ngrid
                 flag1(ind_cell(i))=0
              end do
           end do
        end do
     end do

     ! Count neighbors and set flag2 accordingly
!$omp parallel do firstprivate(icpu,ncache,igrid,ngrid,i,ind_grid,ind,igridn)
     do icpu=1,ncpu
        ncache=reception(icpu,ilevel)%ngrid
        do igrid=1,ncache,nvector
           ngrid=MIN(nvector,ncache-igrid+1)
           do i=1,ngrid
              ind_grid(i)=reception(icpu,ilevel)%igrid(igrid+i-1)
           end do
           call getnborgrids(ind_grid,igridn,ngrid)
           do ind=1,twotondim
              call count_nbors2(igridn,ind,n_nbor(ismooth),ngrid)
           end do
        end do
     end do

     ! Set flag2=1 for cells with flag1=1
!$omp parallel do firstprivate(icpu,ncache,igrid,ngrid,i,ind_grid,ind,iskip,ind_cell)
     do icpu=1,ncpu
        ncache=reception(icpu,ilevel)%ngrid
        do igrid=1,ncache,nvector
           ngrid=MIN(nvector,ncache-igrid+1)
           do i=1,ngrid
              ind_grid(i)=reception(icpu,ilevel)%igrid(igrid+i-1)
           end do
           do ind=1,twotondim
              iskip=ncoarse+(ind-1)*ngridmax
              do i=1,ngrid
                 ind_cell(i)=iskip+ind_grid(i)
              end do
              do i=1,ngrid
                 if(flag1(ind_cell(i))==1)flag2(ind_cell(i))=1
              end do
           end do
        end do
     end do

  end do
  ! End loop over steps
  end do

  ! Compute authorization map for physical boundaries
! if(verbose) write(*,*) ' Entering init_boundary_fine'
  if(simple_boundary)call init_boundary_fine(ilevel)

  ! Restore boundaries for flag1
! if(verbose) write(*,*) ' Entering make_virtual_fine_int'
  call make_virtual_fine_int(flag1(1),ilevel)
  if(simple_boundary .and. verbose) write(*,*) ' Entering make_boundary_flag'
  if(simple_boundary)call make_boundary_flag(ilevel)

111 format('  +Entering authorize_fine for level ',I2)

end subroutine authorize_fine
!################################################################
!################################################################
!################################################################
!################################################################
subroutine make_virtual_coarse_int(xx)
  use amr_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer,dimension(1:ncoarse+ngridmax*twotondim)::xx
  !-----------------------------------------------------------
  ! This routine communicates virtual boundary conditions 
  ! at the coarse level for integer arrays.
  !-----------------------------------------------------------
  integer::nxny,ncell
  integer::i,j,k
  integer::icell,info
  integer,dimension(:),allocatable::ind_cell,fff,ffg

  ! Constants
  nxny=nx*ny
  ncell=  (icoarse_max-icoarse_min+1) &
       & *(jcoarse_max-jcoarse_min+1) &
       & *(kcoarse_max-kcoarse_min+1)

#ifndef WITHOUTMPI
  ! Allocate local arrays
  allocate(ind_cell(1:ncell),fff(1:ncell),ffg(1:ncell))

  ! Compute cell indices
  icell=0
  do k=kcoarse_min,kcoarse_max
  do j=jcoarse_min,jcoarse_max
  do i=icoarse_min,icoarse_max
     icell=icell+1
     ind_cell(icell)=1+i+j*nx+k*nxny
  end do
  end do
  end do
    
  ! Communications
  fff=0; ffg=0
  do icell=1,ncell
     if(cpu_map(ind_cell(icell))==myid)fff(icell)=xx(ind_cell(icell))
  end do
  call MPI_ALLREDUCE(fff,ffg,ncell,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
  do icell=1,ncell
     xx(ind_cell(icell))=ffg(icell)
  end do
     
  ! Dealocate local arrays
  deallocate(ind_cell,fff,ffg)
#endif

end subroutine make_virtual_coarse_int
!################################################################
!################################################################
!################################################################
!################################################################
subroutine make_virtual_fine_dp(xx,ilevel)
  use amr_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
  integer,dimension(MPI_STATUS_SIZE,ncpu)::statuses
#endif
  integer::ilevel
  real(dp),dimension(1:ncoarse+ngridmax*twotondim)::xx
  ! -------------------------------------------------------------------
  ! This routine communicates virtual boundaries among all cpu's.
  ! at level ilevel for any double precision array in the AMR grid.
  ! -------------------------------------------------------------------
  integer::icpu,i,j,ncache,ind,iskip,step
  integer::countsend,countrecv
  integer::info,buf_count,tag=101
  integer,dimension(ncpu)::reqsend,reqrecv

  if(numbtot(1,ilevel)==0)return

  if(ordering=='ksection') then
     call make_virtual_fine_dp_ksec(xx,ilevel)
     return
  end if

#ifndef WITHOUTMPI
  ! Receive all messages
  countrecv=0
  do icpu=1,ncpu
     ncache=reception(icpu,ilevel)%ngrid
     if(ncache>0) then
       countrecv=countrecv+1
       call MPI_IRECV(reception(icpu,ilevel)%u,ncache*twotondim, &
            & MPI_DOUBLE_PRECISION,icpu-1,tag,MPI_COMM_WORLD,reqrecv(countrecv),info)
     end if
  end do
  
  ! Gather emission array
  do icpu=1,ncpu
    if (emission(icpu,ilevel)%ngrid>0) then
      do j=1,twotondim
        step=(j-1)*emission(icpu,ilevel)%ngrid
        iskip=ncoarse+(j-1)*ngridmax
        do i=1,emission(icpu,ilevel)%ngrid
          emission(icpu,ilevel)%u(i+step,1)=xx(emission(icpu,ilevel)%igrid(i)+iskip)
        end do
      end do
    end if
  end do

  ! Send all messages
  countsend=0
  do icpu=1,ncpu
     ncache=emission(icpu,ilevel)%ngrid
     if(ncache>0) then
       countsend=countsend+1
       call MPI_ISEND(emission(icpu,ilevel)%u,ncache*twotondim, &
            & MPI_DOUBLE_PRECISION,icpu-1,tag,MPI_COMM_WORLD,reqsend(countsend),info)
     end if
  end do

  ! Wait for full completion of receives
  call MPI_WAITALL(countrecv,reqrecv,statuses,info)

  ! Scatter reception array
  do icpu=1,ncpu
    if (reception(icpu,ilevel)%ngrid>0) then
      do j=1,twotondim
        step=(j-1)*reception(icpu,ilevel)%ngrid
        iskip=ncoarse+(j-1)*ngridmax
        do i=1,reception(icpu,ilevel)%ngrid
          xx(reception(icpu,ilevel)%igrid(i)+iskip)=reception(icpu,ilevel)%u(i+step,1)
        end do
      end do 
    end if
  end do

  ! Wait for full completion of sends
  call MPI_WAITALL(countsend,reqsend,statuses,info)
#endif

111 format('   Entering make_virtual_fine for level ',I2)

end subroutine make_virtual_fine_dp
!################################################################
!################################################################
!################################################################
!################################################################
subroutine make_virtual_fine_int(xx,ilevel)
  use amr_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
  integer,dimension(MPI_STATUS_SIZE,ncpu)::statuses
#endif
  integer::ilevel
  integer,dimension(1:ncoarse+ngridmax*twotondim)::xx
  ! -------------------------------------------------------------------
  ! This routine communicates virtual boundaries among all cpu's.
  ! at level ilevel for any integer array in the AMR grid.
  ! -------------------------------------------------------------------
  integer::icpu,i,j,ncache,ind,iskip,step
  integer::countsend,countrecv
  integer::info,buf_count,tag=101
  integer,dimension(ncpu)::reqsend,reqrecv

  if(numbtot(1,ilevel)==0)return

  if(ordering=='ksection') then
     call make_virtual_fine_int_ksec(xx,ilevel)
     return
  end if

#ifndef WITHOUTMPI
  ! Receive all messages
  countrecv=0
  do icpu=1,ncpu
     ncache=reception(icpu,ilevel)%ngrid
     if(ncache>0) then
       countrecv=countrecv+1
       call MPI_IRECV(reception(icpu,ilevel)%f,ncache*twotondim, &
            & MPI_INTEGER,icpu-1,tag,MPI_COMM_WORLD,reqrecv(countrecv),info)
     end if
  end do
  
  ! Gather emission array
  do icpu=1,ncpu
    if (emission(icpu,ilevel)%ngrid>0) then
      do j=1,twotondim
        step=(j-1)*emission(icpu,ilevel)%ngrid
        iskip=ncoarse+(j-1)*ngridmax
        do i=1,emission(icpu,ilevel)%ngrid
          emission(icpu,ilevel)%f(i+step,1)=xx(emission(icpu,ilevel)%igrid(i)+iskip)
        end do
      end do
    end if
  end do

  ! Send all messages
  countsend=0
  do icpu=1,ncpu
     ncache=emission(icpu,ilevel)%ngrid
     if(ncache>0) then
       countsend=countsend+1
       call MPI_ISEND(emission(icpu,ilevel)%f,ncache*twotondim, &
            & MPI_INTEGER,icpu-1,tag,MPI_COMM_WORLD,reqsend(countsend),info)
     end if
  end do

  ! Wait for full completion of receives
  call MPI_WAITALL(countrecv,reqrecv,statuses,info)

  ! Scatter reception array
  do icpu=1,ncpu
    if (reception(icpu,ilevel)%ngrid>0) then
      do j=1,twotondim
        step=(j-1)*reception(icpu,ilevel)%ngrid
        iskip=ncoarse+(j-1)*ngridmax
        do i=1,reception(icpu,ilevel)%ngrid
          xx(reception(icpu,ilevel)%igrid(i)+iskip)=reception(icpu,ilevel)%f(i+step,1)
        end do
      end do 
    end if
  end do

  ! Wait for full completion of sends
  call MPI_WAITALL(countsend,reqsend,statuses,info)
#endif

111 format('   Entering make_virtual_fine for level ',I2)

end subroutine make_virtual_fine_int
!################################################################
!################################################################
!################################################################
!################################################################
subroutine make_virtual_reverse_dp(xx,ilevel)
  use amr_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer::ilevel
  real(dp),dimension(1:ncoarse+ngridmax*twotondim)::xx
  ! -------------------------------------------------------------------
  ! This routine communicates virtual boundaries among all cpu's.
  ! at level ilevel in a reverse way for double precision arrays.
  ! -------------------------------------------------------------------
  integer::icpu,i,j,ncache,ind,iskip,step,icell,ibuf
  integer::countsend,countrecv
  integer::info,buf_count,tag=101
  integer,dimension(ncpu)::reqsend,reqrecv
#ifndef WITHOUTMPI
  integer,dimension(MPI_STATUS_SIZE,ncpu)::statuses
  integer::switchlevel=3
#endif

  if(numbtot(1,ilevel)==0)return

  if(ordering=='ksection') then
     call make_virtual_reverse_dp_ksec(xx,ilevel)
     return
  end if

#ifndef WITHOUTMPI
  if(ilevel.LE.switchlevel)then

 ! Gather emission array
  do icpu=1,ncpu
     if (reception(icpu,ilevel)%ngrid>0) then
        do j=1,twotondim
           step=(j-1)*reception(icpu,ilevel)%ngrid
           iskip=ncoarse+(j-1)*ngridmax
           do i=1,reception(icpu,ilevel)%ngrid
              icell=reception(icpu,ilevel)%igrid(i)+iskip
              ibuf=i+step
              reception(icpu,ilevel)%u(ibuf,1)=xx(icell)
           end do
        end do
     end if
  end do

  ! Receive all messages
  countrecv=0
  do icpu=1,myid-1
     ncache=emission(icpu,ilevel)%ngrid
     if(ncache>0) then
        countrecv=countrecv+1
        ! request to send
        call MPI_SEND(countrecv,0, MPI_INTEGER, icpu-1,101,MPI_COMM_WORLD,info)
        call MPI_RECV(emission(icpu,ilevel)%u,ncache*twotondim, &
             & MPI_DOUBLE_PRECISION,icpu-1,tag,MPI_COMM_WORLD,MPI_STATUS_IGNORE,info)
     end if
  end do

  ! Send all messages
  countsend=0
  do icpu=1,ncpu
     ncache=reception(icpu,ilevel)%ngrid
     if(ncache>0) then
        countsend=countsend+1
        ! wait for request to send
        call MPI_RECV(countrecv,0, MPI_INTEGER, icpu-1,101,MPI_COMM_WORLD, &
             & MPI_STATUS_IGNORE, info)
        call MPI_SEND(reception(icpu,ilevel)%u,ncache*twotondim, &
             & MPI_DOUBLE_PRECISION,icpu-1,tag,MPI_COMM_WORLD,info)
     end if
  end do
  
  ! Receive all messages
  countrecv=0
  do icpu=myid+1,ncpu
     ncache=emission(icpu,ilevel)%ngrid
     if(ncache>0) then
        countrecv=countrecv+1
        ! request to send
        call MPI_SEND(countrecv,0, MPI_INTEGER, icpu-1,101,MPI_COMM_WORLD,info)
        call MPI_RECV(emission(icpu,ilevel)%u,ncache*twotondim, &
             & MPI_DOUBLE_PRECISION,icpu-1,tag,MPI_COMM_WORLD,MPI_STATUS_IGNORE,info)
     end if
  end do
  
  ! Scatter reception array
  do icpu=1,ncpu
     if (emission(icpu,ilevel)%ngrid>0) then
        do j=1,twotondim
           step=(j-1)*emission(icpu,ilevel)%ngrid
           iskip=ncoarse+(j-1)*ngridmax
           do i=1,emission(icpu,ilevel)%ngrid
              xx(emission(icpu,ilevel)%igrid(i)+iskip)= &
                   & xx(emission(icpu,ilevel)%igrid(i)+iskip) + emission(icpu,ilevel)%u(i+step,1)
           end do
        end do
     end if
  end do

  else

  ! Receive all messages
  countrecv=0
  do icpu=1,ncpu
     ncache=emission(icpu,ilevel)%ngrid
     if(ncache>0) then
        countrecv=countrecv+1
        call MPI_IRECV(emission(icpu,ilevel)%u,ncache*twotondim, &
             & MPI_DOUBLE_PRECISION,icpu-1,tag,MPI_COMM_WORLD,reqrecv(countrecv),info)
     end if
  end do

  ! Gather emission array
  do icpu=1,ncpu
     if (reception(icpu,ilevel)%ngrid>0) then
        do j=1,twotondim
           step=(j-1)*reception(icpu,ilevel)%ngrid
           iskip=ncoarse+(j-1)*ngridmax
           do i=1,reception(icpu,ilevel)%ngrid
              reception(icpu,ilevel)%u(i+step,1)=xx(reception(icpu,ilevel)%igrid(i)+iskip)
           end do
        end do
     end if
  end do

  ! Send all messages
  countsend=0
  do icpu=1,ncpu
     ncache=reception(icpu,ilevel)%ngrid
     if(ncache>0) then
       countsend=countsend+1
        call MPI_ISEND(reception(icpu,ilevel)%u,ncache*twotondim, &
             & MPI_DOUBLE_PRECISION,icpu-1,tag,MPI_COMM_WORLD,reqsend(countsend),info)
     end if
  end do

  ! Wait for full completion of receives
  call MPI_WAITALL(countrecv,reqrecv,statuses,info)

  ! Scatter reception array
  do icpu=1,ncpu
    if (emission(icpu,ilevel)%ngrid>0) then
      do j=1,twotondim
        step=(j-1)*emission(icpu,ilevel)%ngrid
        iskip=ncoarse+(j-1)*ngridmax
        do i=1,emission(icpu,ilevel)%ngrid
           xx(emission(icpu,ilevel)%igrid(i)+iskip)= &
                & xx(emission(icpu,ilevel)%igrid(i)+iskip) + emission(icpu,ilevel)%u(i+step,1)
        end do
      end do 
    end if
  end do

  ! Wait for full completion of sends
  call MPI_WAITALL(countsend,reqsend,statuses,info)

  endif

#endif

111 format('   Entering make_virtual_reverse for level ',I2)

end subroutine make_virtual_reverse_dp
!################################################################
!################################################################
!################################################################
!################################################################
subroutine make_virtual_reverse_int(xx,ilevel)
  use amr_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer::ilevel
  integer,dimension(1:ncoarse+ngridmax*twotondim)::xx
  ! -------------------------------------------------------------------
  ! This routine communicates virtual boundaries among all cpu's.
  ! at level ilevel in a reverse way for integer arrays.
  ! -------------------------------------------------------------------
  integer::icpu,i,j,ncache,ind,iskip,step,icell,ibuf
  integer::countsend,countrecv
  integer::info,buf_count,tag=101
  integer,dimension(ncpu)::reqsend,reqrecv
#ifndef WITHOUTMPI
  integer,dimension(MPI_STATUS_SIZE,ncpu)::statuses
  integer::switchlevel=3
#endif

  if(numbtot(1,ilevel)==0)return

  if(ordering=='ksection') then
     call make_virtual_reverse_int_ksec(xx,ilevel)
     return
  end if

#ifndef WITHOUTMPI

  if(ilevel.le.switchlevel) then

  ! Gather emission array
  do icpu=1,ncpu
     if (reception(icpu,ilevel)%ngrid>0) then
        do j=1,twotondim
           step=(j-1)*reception(icpu,ilevel)%ngrid
           iskip=ncoarse+(j-1)*ngridmax
           do i=1,reception(icpu,ilevel)%ngrid
              icell=reception(icpu,ilevel)%igrid(i)+iskip
              ibuf=i+step
              reception(icpu,ilevel)%f(ibuf,1)=xx(icell)
           end do
        end do
     end if
  end do
  
  ! Receive all messages
  countrecv=0
  do icpu=1,myid-1
     ncache=emission(icpu,ilevel)%ngrid
     if(ncache>0) then
        countrecv=countrecv+1
        ! request to send 
        call MPI_SEND(countrecv,0, MPI_INTEGER, icpu-1,101,MPI_COMM_WORLD,info)
        call MPI_RECV(emission(icpu,ilevel)%f,ncache*twotondim, &
             & MPI_INTEGER,icpu-1,tag,MPI_COMM_WORLD,MPI_STATUS_IGNORE,info)
     end if
  end do
  
  ! Send all messages
  countsend=0
  do icpu=1,ncpu
     ncache=reception(icpu,ilevel)%ngrid
     if(ncache>0) then
        countsend=countsend+1
        ! wait for request to send
        call MPI_RECV(countrecv,0, MPI_INTEGER, icpu-1,101,MPI_COMM_WORLD, &
             & MPI_STATUS_IGNORE, info)
        call MPI_SEND(reception(icpu,ilevel)%f,ncache*twotondim, &
             & MPI_INTEGER,icpu-1,tag,MPI_COMM_WORLD,info)
     end if
  end do
  
  ! Receive all messages
  countrecv=0
  do icpu=myid+1,ncpu
     ncache=emission(icpu,ilevel)%ngrid
     if(ncache>0) then
        countrecv=countrecv+1
        ! request to send
        call MPI_SEND(countrecv,0, MPI_INTEGER, icpu-1,101,MPI_COMM_WORLD,info)
        call MPI_RECV(emission(icpu,ilevel)%f,ncache*twotondim, &
             & MPI_INTEGER,icpu-1,tag,MPI_COMM_WORLD,MPI_STATUS_IGNORE,info)
     end if
  end do
  
  ! Scatter reception array
  do icpu=1,ncpu
     if (emission(icpu,ilevel)%ngrid>0) then
        do j=1,twotondim
           step=(j-1)*emission(icpu,ilevel)%ngrid
           iskip=ncoarse+(j-1)*ngridmax
           do i=1,emission(icpu,ilevel)%ngrid
              xx(emission(icpu,ilevel)%igrid(i)+iskip)= &
                   & xx(emission(icpu,ilevel)%igrid(i)+iskip) + emission(icpu,ilevel)%f(i+step,1)
           end do
        end do
     end if
  end do
  
  else

  ! Receive all messages
  countrecv=0
  do icpu=1,ncpu
     ncache=emission(icpu,ilevel)%ngrid
     if(ncache>0) then
        countrecv=countrecv+1
        call MPI_IRECV(emission(icpu,ilevel)%f,ncache*twotondim, &
             & MPI_INTEGER,icpu-1,tag,MPI_COMM_WORLD,reqrecv(countrecv),info)
     end if
  end do

  ! Gather emission array
  do icpu=1,ncpu
     if (reception(icpu,ilevel)%ngrid>0) then
        do j=1,twotondim
           step=(j-1)*reception(icpu,ilevel)%ngrid
           iskip=ncoarse+(j-1)*ngridmax
           do i=1,reception(icpu,ilevel)%ngrid
              reception(icpu,ilevel)%f(i+step,1)=xx(reception(icpu,ilevel)%igrid(i)+iskip)
           end do
        end do
     end if
  end do
  
  ! Send all messages
  countsend=0
  do icpu=1,ncpu
     ncache=reception(icpu,ilevel)%ngrid
     if(ncache>0) then
       countsend=countsend+1
        call MPI_ISEND(reception(icpu,ilevel)%f,ncache*twotondim, &
             & MPI_INTEGER,icpu-1,tag,MPI_COMM_WORLD,reqsend(countsend),info)
     end if
  end do

  ! Wait for full completion of receives
  call MPI_WAITALL(countrecv,reqrecv,statuses,info)

  ! Scatter reception array
  do icpu=1,ncpu
     if (emission(icpu,ilevel)%ngrid>0) then
        do j=1,twotondim
           step=(j-1)*emission(icpu,ilevel)%ngrid
           iskip=ncoarse+(j-1)*ngridmax
           do i=1,emission(icpu,ilevel)%ngrid
              xx(emission(icpu,ilevel)%igrid(i)+iskip)= &
                   & xx(emission(icpu,ilevel)%igrid(i)+iskip) + emission(icpu,ilevel)%f(i+step,1)
           end do
        end do
     end if
  end do
  
  ! Wait for full completion of sends
  call MPI_WAITALL(countsend,reqsend,statuses,info)

  endif

#endif

111 format('   Entering make_virtual_reverse for level ',I2)

end subroutine make_virtual_reverse_int
!################################################################
!################################################################
!################################################################
!################################################################
subroutine build_comm(ilevel)
  use amr_commons
  use poisson_commons, only: lookup_mg
  use ksection
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer::ilevel
  ! -------------------------------------------------------------------
  ! This routine builds the communication structure for level ilevel.
  ! Array flag2 is used as temporary work space.
  ! -------------------------------------------------------------------
  integer::icpu,ibound
  integer::ncache,ind,iskip
  integer::i,j,k,nxny
  integer::igrid,jgrid,ngrid
  integer::info,tag=101
  integer,dimension(ncpu)::reqsend,reqrecv
#ifndef WITHOUTMPI
  integer,dimension(ncpu)::sendbuf,recvbuf
  integer,dimension(MPI_STATUS_SIZE,ncpu)::statuses
  integer::countsend,countrecv
  ! Ksection variables
  integer::ntotal_bc,nrecv_bc,idx_bc,sender_bc
  real(dp),allocatable::sendbuf_bc(:,:),recvbuf_bc(:,:)
  integer,allocatable::destcpu_bc(:),emit_count(:)
#endif
  integer,dimension(1:nvector),save::ind_grid,ind_cell

  if(verbose)write(*,111)ilevel
  nxny=nx*ny

  !----------------------------------------------------------------
  ! Compute grids global adress using flag2 array at level ilevel-1
  !----------------------------------------------------------------
  if(ilevel==1)then
     do k=kcoarse_min,kcoarse_max
     do j=jcoarse_min,jcoarse_max
     do i=icoarse_min,icoarse_max
        ind=1+i+j*nx+k*nxny
        if(cpu_map(ind)==myid)then
           flag2(ind)=son(ind)
        else
           flag2(ind)=0
        end if
     end do
     end do
     end do    
     call make_virtual_coarse_int(flag2(1))
  else
     ! Initialize flag2 to local adress for cpu map = myid cells
     ncache=active(ilevel-1)%ngrid
     do igrid=1,ncache,nvector
        ngrid=MIN(nvector,ncache-igrid+1)
        do i=1,ngrid
           ind_grid(i)=active(ilevel-1)%igrid(igrid+i-1)
        end do
        do ind=1,twotondim
           iskip=ncoarse+(ind-1)*ngridmax
           do i=1,ngrid
              ind_cell(i)=iskip+ind_grid(i)
           end do
           do i=1,ngrid
              if(cpu_map(ind_cell(i))==myid)then
                 flag2(ind_cell(i))=son(ind_cell(i))
              else
                 flag2(ind_cell(i))=0
              end if
           end do
        end do
     end do
     do icpu=1,ncpu
        ncache=reception(icpu,ilevel-1)%ngrid
        do igrid=1,ncache,nvector
           ngrid=MIN(nvector,ncache-igrid+1)
           do i=1,ngrid
              ind_grid(i)=reception(icpu,ilevel-1)%igrid(igrid+i-1)
           end do
           do ind=1,twotondim
              iskip=ncoarse+(ind-1)*ngridmax
              do i=1,ngrid
                 ind_cell(i)=iskip+ind_grid(i)
              end do
              do i=1,ngrid
                 if(cpu_map(ind_cell(i))==myid)then
                    flag2(ind_cell(i))=son(ind_cell(i))
                 else
                    flag2(ind_cell(i))=0
                 end if
              end do
           end do
        end do
     end do
     call make_virtual_reverse_int(flag2(1),ilevel-1)
     call make_virtual_fine_int   (flag2(1),ilevel-1)
  end if

  !--------------------------------------------------------
  ! Compute number and index of active grid at level ilevel
  !--------------------------------------------------------
  ncache=numbl(myid,ilevel)
  ! Reset old communicator
  if(active(ilevel)%ngrid>0)then
     active(ilevel)%ngrid=0
     deallocate(active(ilevel)%igrid)
  end if
  if(ncache>0)then     
     ! Allocate grid index to new communicator
     active(ilevel)%ngrid=ncache
     allocate(active(ilevel)%igrid(1:ncache))
     ! Gather all grids
     igrid=headl(myid,ilevel)
     do jgrid=1,numbl(myid,ilevel)
        active(ilevel)%igrid(jgrid)=igrid
        igrid=next(igrid)
     end do
  end if
  ! Fill up lookup_mg for active
  if(poisson)then
     igrid=headl(myid,ilevel)
     do jgrid=1,numbl(myid,ilevel)
        lookup_mg(igrid)=0
        igrid=next(igrid)
     end do
  end if

  !----------------------------------------------------
  ! Compute number and index of physical boundary grids
  !----------------------------------------------------
  do ibound=1,nboundary
     ncache=numbb(ibound,ilevel)
     ! Reset old communicator
     if(boundary(ibound,ilevel)%ngrid>0)then
        boundary(ibound,ilevel)%ngrid=0
        deallocate(boundary(ibound,ilevel)%igrid)
     end if
     if(ncache>0)then     
        ! Allocate grid index to new communicator
        boundary(ibound,ilevel)%ngrid=ncache
        allocate(boundary(ibound,ilevel)%igrid(1:ncache))
        ! Gather all grids
        igrid=headb(ibound,ilevel)
        do jgrid=1,numbb(ibound,ilevel)
           boundary(ibound,ilevel)%igrid(jgrid)=igrid
           igrid=next(igrid)
        end do
     end if
  end do

  !----------------------------------------------------
  ! Compute number and index of virtual boundary grids
  !----------------------------------------------------
#ifndef WITHOUTMPI
   do icpu=1,ncpu
      ncache=0
      if(icpu.ne.myid)ncache=numbl(icpu,ilevel)
      ! Reset old communicators
      if(emission(icpu,ilevel)%ngrid>0)then
         emission(icpu,ilevel)%ngrid=0
         deallocate(emission(icpu,ilevel)%igrid)
         deallocate(emission(icpu,ilevel)%u)
         deallocate(emission(icpu,ilevel)%f)
      end if
      if(reception(icpu,ilevel)%ngrid>0)then
         reception(icpu,ilevel)%ngrid=0
         deallocate(reception(icpu,ilevel)%igrid)
         deallocate(reception(icpu,ilevel)%u)
         deallocate(reception(icpu,ilevel)%f)
      end if
      if(ncache>0)then
         ! Allocate grid index to new communicator
         reception(icpu,ilevel)%ngrid=ncache
         allocate(reception(icpu,ilevel)%igrid(1:ncache))
         ! Gather all grids
         igrid=headl(icpu,ilevel)
         do jgrid=1,numbl(icpu,ilevel)
            reception(icpu,ilevel)%igrid(jgrid)=igrid
            igrid=next(igrid)
         end do
         ! Allocate temporary communication buffer
         allocate(reception(icpu,ilevel)%f(1:ncache,1:1))
         do i=1,ncache
            reception(icpu,ilevel)%f(i,1) = &
            & flag2(father(reception(icpu,ilevel)%igrid(i)))
         end do
         ! Fill up lookup_mg for reception
         if(poisson)then
            do i=1,ncache
               lookup_mg(reception(icpu,ilevel)%igrid(i))= -reception(icpu,ilevel)%f(i,1)
            end do
         end if
      end if
      sendbuf(icpu)=reception(icpu,ilevel)%ngrid
   end do

  !--------------------------------------------------------
  ! Communicate virtual grid number and index to parent cpu
  !--------------------------------------------------------
  if(ordering=='ksection') then

  ! Ksection version: exchange grid indices via tree routing
  ! Pack: (sender_id, reception_index, global_grid_address) for each reception grid
  ntotal_bc = 0
  do icpu = 1, ncpu
     ntotal_bc = ntotal_bc + reception(icpu,ilevel)%ngrid
  end do

  allocate(sendbuf_bc(1:3, 1:max(ntotal_bc,1)))
  allocate(destcpu_bc(1:max(ntotal_bc,1)))
  idx_bc = 0
  do icpu = 1, ncpu
     do i = 1, reception(icpu,ilevel)%ngrid
        idx_bc = idx_bc + 1
        destcpu_bc(idx_bc) = icpu
        sendbuf_bc(1, idx_bc) = dble(myid)
        sendbuf_bc(2, idx_bc) = dble(i)
        sendbuf_bc(3, idx_bc) = dble(reception(icpu,ilevel)%f(i,1))
     end do
  end do

  call ksection_exchange_dp(sendbuf_bc, ntotal_bc, destcpu_bc, 3, &
       & recvbuf_bc, nrecv_bc)
  deallocate(sendbuf_bc, destcpu_bc)

  ! Count emission grids per sender
  allocate(emit_count(1:ncpu))
  emit_count = 0
  do i = 1, nrecv_bc
     sender_bc = nint(recvbuf_bc(1, i))
     emit_count(sender_bc) = emit_count(sender_bc) + 1
  end do

  ! Allocate emission igrid
  do icpu = 1, ncpu
     emission(icpu, ilevel)%ngrid = emit_count(icpu)
     if(emit_count(icpu) > 0) &
          & allocate(emission(icpu, ilevel)%igrid(1:emit_count(icpu)))
  end do
  deallocate(emit_count)

  ! Fill emission igrid using reception index as position
  do i = 1, nrecv_bc
     sender_bc = nint(recvbuf_bc(1, i))
     j = nint(recvbuf_bc(2, i))
     emission(sender_bc, ilevel)%igrid(j) = nint(recvbuf_bc(3, i))
  end do
  deallocate(recvbuf_bc)

  else

  ! Original MPI_ALLTOALL version
  call MPI_ALLTOALL(sendbuf,1,MPI_INTEGER,recvbuf,1,MPI_INTEGER,MPI_COMM_WORLD,info)

  ! Allocate grid index
  do icpu=1,ncpu
     emission(icpu,ilevel)%ngrid=recvbuf(icpu)
     ncache=emission(icpu,ilevel)%ngrid
     if(ncache>0)allocate(emission(icpu,ilevel)%igrid(1:ncache))
  end do

  ! Receive grid list
  countrecv=0
  do icpu=1,ncpu
     ncache=emission(icpu,ilevel)%ngrid
     if(ncache>0) then
        countrecv=countrecv+1
        call MPI_IRECV(emission(icpu,ilevel)%igrid,ncache, &
             & MPI_INTEGER,icpu-1,tag,MPI_COMM_WORLD,reqrecv(countrecv),info)
     end if
  end do

  ! Send global index
  countsend=0
  do icpu=1,ncpu
     ncache=reception(icpu,ilevel)%ngrid
     if(ncache>0) then
        countsend=countsend+1
        call MPI_ISEND(reception(icpu,ilevel)%f,ncache, &
             & MPI_INTEGER,icpu-1,tag,MPI_COMM_WORLD,reqsend(countsend),info)
     end if
  end do

  ! Wait for full completion of sends
  call MPI_WAITALL(countsend,reqsend,statuses,info)

  ! Wait for full completion of receives
  call MPI_WAITALL(countrecv,reqrecv,statuses,info)

  end if

  ! Deallocate temporary communication buffers
  do icpu=1,ncpu
     ncache=reception(icpu,ilevel)%ngrid
     if(ncache>0)deallocate(reception(icpu,ilevel)%f)
  end do

  ! Allocate temporary communication buffers
  do icpu=1,ncpu
     ncache=emission(icpu,ilevel)%ngrid
     if(ncache>0)then
        allocate(emission(icpu,ilevel)%u(1:ncache*twotondim,1:1))
        allocate(emission(icpu,ilevel)%f(1:ncache*twotondim,1:1))
     endif
     ncache=reception(icpu,ilevel)%ngrid
     if(ncache>0)then
        allocate(reception(icpu,ilevel)%u(1:ncache*twotondim,1:1))
        allocate(reception(icpu,ilevel)%f(1:ncache*twotondim,1:1))
     endif
  end do

#endif

  ! Update Morton hash table for this level (includes new virtual grids)
  call morton_hash_rebuild_level(ilevel)

111 format('   Entering build_comm for level ',I2)

end subroutine build_comm
!################################################################
!################################################################
!################################################################
!################################################################
subroutine make_virtual_fine_dp_ksec(xx,ilevel)
  use amr_commons
  use ksection
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer::ilevel
  real(dp),dimension(1:ncoarse+ngridmax*twotondim)::xx
  ! -------------------------------------------------------------------
  ! Ksection-based forward ghost zone exchange for double precision.
  ! Packs emission data with metadata, exchanges via ksection tree,
  ! then scatters to reception grids using metadata.
  ! -------------------------------------------------------------------
  integer::icpu,i,j,idx,ntotal,nrecv,nprops_ksec,sender,ridx,igrid
  real(dp),allocatable::sendbuf(:,:),recvbuf(:,:)
  integer,allocatable::dest_cpu(:)

#ifndef WITHOUTMPI
  nprops_ksec = twotondim + 2

  ! Count total emission items
  ntotal = 0
  do icpu = 1, ncpu
     ntotal = ntotal + emission(icpu,ilevel)%ngrid
  end do

  ! Pack sendbuf + dest_cpu
  allocate(sendbuf(1:nprops_ksec, 1:max(ntotal,1)))
  allocate(dest_cpu(1:max(ntotal,1)))
  idx = 0
  do icpu = 1, ncpu
     do i = 1, emission(icpu,ilevel)%ngrid
        idx = idx + 1
        dest_cpu(idx) = icpu
        do j = 1, twotondim
           sendbuf(j, idx) = xx(emission(icpu,ilevel)%igrid(i) &
                & + ncoarse + (j-1)*ngridmax)
        end do
        sendbuf(twotondim+1, idx) = dble(myid)
        sendbuf(twotondim+2, idx) = dble(i)
     end do
  end do

  ! Exchange via ksection tree
  call ksection_exchange_dp(sendbuf, ntotal, dest_cpu, nprops_ksec, &
       & recvbuf, nrecv)

  ! Scatter received data to reception grids
  do i = 1, nrecv
     sender = nint(recvbuf(twotondim+1, i))
     ridx   = nint(recvbuf(twotondim+2, i))
     igrid  = reception(sender, ilevel)%igrid(ridx)
     do j = 1, twotondim
        xx(igrid + ncoarse + (j-1)*ngridmax) = recvbuf(j, i)
     end do
  end do

  deallocate(sendbuf, dest_cpu, recvbuf)
#endif

end subroutine make_virtual_fine_dp_ksec
!################################################################
!################################################################
!################################################################
!################################################################
subroutine make_virtual_reverse_dp_ksec(xx,ilevel)
  use amr_commons
  use ksection
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer::ilevel
  real(dp),dimension(1:ncoarse+ngridmax*twotondim)::xx
  ! -------------------------------------------------------------------
  ! Ksection-based reverse ghost zone exchange for double precision.
  ! Packs reception data with metadata, exchanges via ksection tree,
  ! then accumulates (+=) into owner's emission grids.
  ! -------------------------------------------------------------------
  integer::icpu,i,j,idx,ntotal,nrecv,nprops_ksec,sender,eidx,igrid
  real(dp),allocatable::sendbuf(:,:),recvbuf(:,:)
  integer,allocatable::dest_cpu(:)

#ifndef WITHOUTMPI
  nprops_ksec = twotondim + 2

  ! Count total reception items
  ntotal = 0
  do icpu = 1, ncpu
     ntotal = ntotal + reception(icpu,ilevel)%ngrid
  end do

  ! Pack sendbuf from reception grids + dest_cpu
  allocate(sendbuf(1:nprops_ksec, 1:max(ntotal,1)))
  allocate(dest_cpu(1:max(ntotal,1)))
  idx = 0
  do icpu = 1, ncpu
     do i = 1, reception(icpu,ilevel)%ngrid
        idx = idx + 1
        dest_cpu(idx) = icpu
        do j = 1, twotondim
           sendbuf(j, idx) = xx(reception(icpu,ilevel)%igrid(i) &
                & + ncoarse + (j-1)*ngridmax)
        end do
        sendbuf(twotondim+1, idx) = dble(myid)
        sendbuf(twotondim+2, idx) = dble(i)
     end do
  end do

  ! Exchange via ksection tree
  call ksection_exchange_dp(sendbuf, ntotal, dest_cpu, nprops_ksec, &
       & recvbuf, nrecv)

  ! Accumulate received data into emission grids
  do i = 1, nrecv
     sender = nint(recvbuf(twotondim+1, i))
     eidx   = nint(recvbuf(twotondim+2, i))
     igrid  = emission(sender, ilevel)%igrid(eidx)
     do j = 1, twotondim
        xx(igrid + ncoarse + (j-1)*ngridmax) = &
             & xx(igrid + ncoarse + (j-1)*ngridmax) + recvbuf(j, i)
     end do
  end do

  deallocate(sendbuf, dest_cpu, recvbuf)
#endif

end subroutine make_virtual_reverse_dp_ksec
!################################################################
!################################################################
!################################################################
!################################################################
subroutine make_virtual_fine_int_ksec(xx,ilevel)
  use amr_commons
  use ksection
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer::ilevel
  integer,dimension(1:ncoarse+ngridmax*twotondim)::xx
  ! -------------------------------------------------------------------
  ! Ksection-based forward ghost zone exchange for integer arrays.
  ! Converts int to dp, exchanges via ksection tree, converts back.
  ! -------------------------------------------------------------------
  integer::icpu,i,j,idx,ntotal,nrecv,nprops_ksec,sender,ridx,igrid
  real(dp),allocatable::sendbuf(:,:),recvbuf(:,:)
  integer,allocatable::dest_cpu(:)

#ifndef WITHOUTMPI
  nprops_ksec = twotondim + 2

  ! Count total emission items
  ntotal = 0
  do icpu = 1, ncpu
     ntotal = ntotal + emission(icpu,ilevel)%ngrid
  end do

  ! Pack sendbuf (int->dp) + dest_cpu
  allocate(sendbuf(1:nprops_ksec, 1:max(ntotal,1)))
  allocate(dest_cpu(1:max(ntotal,1)))
  idx = 0
  do icpu = 1, ncpu
     do i = 1, emission(icpu,ilevel)%ngrid
        idx = idx + 1
        dest_cpu(idx) = icpu
        do j = 1, twotondim
           sendbuf(j, idx) = dble(xx(emission(icpu,ilevel)%igrid(i) &
                & + ncoarse + (j-1)*ngridmax))
        end do
        sendbuf(twotondim+1, idx) = dble(myid)
        sendbuf(twotondim+2, idx) = dble(i)
     end do
  end do

  ! Exchange via ksection tree
  call ksection_exchange_dp(sendbuf, ntotal, dest_cpu, nprops_ksec, &
       & recvbuf, nrecv)

  ! Scatter received data (dp->int) to reception grids
  do i = 1, nrecv
     sender = nint(recvbuf(twotondim+1, i))
     ridx   = nint(recvbuf(twotondim+2, i))
     igrid  = reception(sender, ilevel)%igrid(ridx)
     do j = 1, twotondim
        xx(igrid + ncoarse + (j-1)*ngridmax) = nint(recvbuf(j, i))
     end do
  end do

  deallocate(sendbuf, dest_cpu, recvbuf)
#endif

end subroutine make_virtual_fine_int_ksec
!################################################################
!################################################################
!################################################################
!################################################################
subroutine make_virtual_reverse_int_ksec(xx,ilevel)
  use amr_commons
  use ksection
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer::ilevel
  integer,dimension(1:ncoarse+ngridmax*twotondim)::xx
  ! -------------------------------------------------------------------
  ! Ksection-based reverse ghost zone exchange for integer arrays.
  ! Converts int to dp, exchanges via ksection tree, then accumulates.
  ! -------------------------------------------------------------------
  integer::icpu,i,j,idx,ntotal,nrecv,nprops_ksec,sender,eidx,igrid
  real(dp),allocatable::sendbuf(:,:),recvbuf(:,:)
  integer,allocatable::dest_cpu(:)

#ifndef WITHOUTMPI
  nprops_ksec = twotondim + 2

  ! Count total reception items
  ntotal = 0
  do icpu = 1, ncpu
     ntotal = ntotal + reception(icpu,ilevel)%ngrid
  end do

  ! Pack sendbuf from reception grids (int->dp) + dest_cpu
  allocate(sendbuf(1:nprops_ksec, 1:max(ntotal,1)))
  allocate(dest_cpu(1:max(ntotal,1)))
  idx = 0
  do icpu = 1, ncpu
     do i = 1, reception(icpu,ilevel)%ngrid
        idx = idx + 1
        dest_cpu(idx) = icpu
        do j = 1, twotondim
           sendbuf(j, idx) = dble(xx(reception(icpu,ilevel)%igrid(i) &
                & + ncoarse + (j-1)*ngridmax))
        end do
        sendbuf(twotondim+1, idx) = dble(myid)
        sendbuf(twotondim+2, idx) = dble(i)
     end do
  end do

  ! Exchange via ksection tree
  call ksection_exchange_dp(sendbuf, ntotal, dest_cpu, nprops_ksec, &
       & recvbuf, nrecv)

  ! Accumulate received data (dp->int) into emission grids
  do i = 1, nrecv
     sender = nint(recvbuf(twotondim+1, i))
     eidx   = nint(recvbuf(twotondim+2, i))
     igrid  = emission(sender, ilevel)%igrid(eidx)
     do j = 1, twotondim
        xx(igrid + ncoarse + (j-1)*ngridmax) = &
             & xx(igrid + ncoarse + (j-1)*ngridmax) + nint(recvbuf(j, i))
     end do
  end do

  deallocate(sendbuf, dest_cpu, recvbuf)
#endif

end subroutine make_virtual_reverse_int_ksec
!################################################################
!################################################################
!################################################################
!################################################################
