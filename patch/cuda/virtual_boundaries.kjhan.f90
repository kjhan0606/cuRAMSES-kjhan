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
!$omp parallel do private(icpu) firstprivate(ncache,igrid,ngrid,i,ind_grid,ind,iskip,ind_cell)
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
!$omp parallel do private(icpu) firstprivate(ncache,igrid,ngrid,i,ind_grid,ind,igridn)
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
!$omp parallel do private(icpu) firstprivate(ncache,igrid,ngrid,i,ind_grid,ind,iskip,ind_cell)
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
  logical::use_ksec
  real(dp)::t1

  if(numbtot(1,ilevel)==0)return

  ! Auto-tune dispatch: P2P vs K-Section (comp 1: fine_dp)
  use_ksec = .false.
  if(ordering=='ksection') then
     if(exchange_method=='ksection') then
        use_ksec = .true.
     else if(exchange_method=='auto') then
        use_ksec = (xchg_phase(1)==1) .or. &
                   (xchg_phase(1)==2 .and. xchg_chosen(1)==1) .or. &
                   (xchg_phase(1)==3 .and. xchg_chosen(1)==0)
     end if
  end if
  if(use_ksec) then
#ifndef WITHOUTMPI
     t1 = MPI_WTIME()
#endif
     call make_virtual_fine_dp_ksec(xx,ilevel)
#ifndef WITHOUTMPI
     if(exchange_method=='auto') call xchg_autotune_update(1, MPI_WTIME()-t1)
#endif
     return
  end if

#ifndef WITHOUTMPI
  if(ordering=='ksection') t1 = MPI_WTIME()
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
  if(ordering=='ksection' .and. exchange_method=='auto') &
     call xchg_autotune_update(1, MPI_WTIME()-t1)
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
  logical::use_ksec
  real(dp)::t1

  if(numbtot(1,ilevel)==0)return

  ! Auto-tune dispatch: P2P vs K-Section (comp 2: fine_int)
  use_ksec = .false.
  if(ordering=='ksection') then
     if(exchange_method=='ksection') then
        use_ksec = .true.
     else if(exchange_method=='auto') then
        use_ksec = (xchg_phase(2)==1) .or. &
                   (xchg_phase(2)==2 .and. xchg_chosen(2)==1) .or. &
                   (xchg_phase(2)==3 .and. xchg_chosen(2)==0)
     end if
  end if
  if(use_ksec) then
#ifndef WITHOUTMPI
     t1 = MPI_WTIME()
#endif
     call make_virtual_fine_int_ksec(xx,ilevel)
#ifndef WITHOUTMPI
     if(exchange_method=='auto') call xchg_autotune_update(2, MPI_WTIME()-t1)
#endif
     return
  end if

#ifndef WITHOUTMPI
  if(ordering=='ksection') t1 = MPI_WTIME()
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
  if(ordering=='ksection' .and. exchange_method=='auto') &
     call xchg_autotune_update(2, MPI_WTIME()-t1)
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
  logical::use_ksec
  real(dp)::t1
#ifndef WITHOUTMPI
  integer,dimension(MPI_STATUS_SIZE,ncpu)::statuses
  integer::switchlevel=3
#endif

  if(numbtot(1,ilevel)==0)return

  ! Auto-tune dispatch: P2P vs K-Section (comp 3: reverse_dp)
  use_ksec = .false.
  if(ordering=='ksection') then
     if(exchange_method=='ksection') then
        use_ksec = .true.
     else if(exchange_method=='auto') then
        use_ksec = (xchg_phase(3)==1) .or. &
                   (xchg_phase(3)==2 .and. xchg_chosen(3)==1) .or. &
                   (xchg_phase(3)==3 .and. xchg_chosen(3)==0)
     end if
  end if
  if(use_ksec) then
#ifndef WITHOUTMPI
     t1 = MPI_WTIME()
#endif
     call make_virtual_reverse_dp_ksec(xx,ilevel)
#ifndef WITHOUTMPI
     if(exchange_method=='auto') call xchg_autotune_update(3, MPI_WTIME()-t1)
#endif
     return
  end if

#ifndef WITHOUTMPI
  if(ordering=='ksection') t1 = MPI_WTIME()
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

  if(ordering=='ksection' .and. exchange_method=='auto') &
     call xchg_autotune_update(3, MPI_WTIME()-t1)
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
  logical::use_ksec
  real(dp)::t1
#ifndef WITHOUTMPI
  integer,dimension(MPI_STATUS_SIZE,ncpu)::statuses
  integer::switchlevel=3
#endif

  if(numbtot(1,ilevel)==0)return

  ! Auto-tune dispatch: P2P vs K-Section (comp 4: reverse_int)
  use_ksec = .false.
  if(ordering=='ksection') then
     if(exchange_method=='ksection') then
        use_ksec = .true.
     else if(exchange_method=='auto') then
        use_ksec = (xchg_phase(4)==1) .or. &
                   (xchg_phase(4)==2 .and. xchg_chosen(4)==1) .or. &
                   (xchg_phase(4)==3 .and. xchg_chosen(4)==0)
     end if
  end if
  if(use_ksec) then
#ifndef WITHOUTMPI
     t1 = MPI_WTIME()
#endif
     call make_virtual_reverse_int_ksec(xx,ilevel)
#ifndef WITHOUTMPI
     if(exchange_method=='auto') call xchg_autotune_update(4, MPI_WTIME()-t1)
#endif
     return
  end if

#ifndef WITHOUTMPI
  if(ordering=='ksection') t1 = MPI_WTIME()

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

  if(ordering=='ksection' .and. exchange_method=='auto') &
     call xchg_autotune_update(4, MPI_WTIME()-t1)
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
subroutine make_virtual_fine_int_pair(xx1,xx2,ilevel)
  use amr_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer::ilevel
  integer,dimension(1:ncoarse+ngridmax*twotondim)::xx1,xx2
  ! -------------------------------------------------------------------
  ! Exchange two integer arrays simultaneously in a single communication
  ! round. Used for cpu_map + cpu_map2 in load_balance expand pass.
  ! -------------------------------------------------------------------
  logical::use_ksec
  real(dp)::t1

  if(numbtot(1,ilevel)==0)return

  ! Auto-tune dispatch: K-Section bulk vs P2P fallback (comp 5: pair_int)
  use_ksec = .false.
  if(ordering=='ksection') then
     if(exchange_method=='ksection') then
        use_ksec = .true.
     else if(exchange_method=='auto') then
        use_ksec = (xchg_phase(5)==1) .or. &
                   (xchg_phase(5)==2 .and. xchg_chosen(5)==1) .or. &
                   (xchg_phase(5)==3 .and. xchg_chosen(5)==0)
     end if
  end if
  if(use_ksec) then
#ifndef WITHOUTMPI
     t1 = MPI_WTIME()
#endif
     call make_virtual_fine_int_pair_ksec(xx1,xx2,ilevel)
#ifndef WITHOUTMPI
     if(exchange_method=='auto') call xchg_autotune_update(5, MPI_WTIME()-t1)
#endif
     return
  end if

  ! Fallback: two separate calls (each uses its own auto-tune)
#ifndef WITHOUTMPI
  if(ordering=='ksection') t1 = MPI_WTIME()
#endif
  call make_virtual_fine_int(xx1,ilevel)
  call make_virtual_fine_int(xx2,ilevel)
#ifndef WITHOUTMPI
  if(ordering=='ksection' .and. exchange_method=='auto') &
     call xchg_autotune_update(5, MPI_WTIME()-t1)
#endif

end subroutine make_virtual_fine_int_pair
!################################################################
!################################################################
!################################################################
!################################################################
subroutine make_virtual_fine_int_pair_ksec(xx1,xx2,ilevel)
  use amr_commons
  use ksection
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer::ilevel
  integer,dimension(1:ncoarse+ngridmax*twotondim)::xx1,xx2
  ! -------------------------------------------------------------------
  ! Ksection-based forward exchange for two integer arrays at once.
  ! Packs 2*twotondim cells + 2 metadata per emission grid.
  ! -------------------------------------------------------------------
  integer::icpu,i,j,idx,ntotal,nrecv,nprops_ksec,sender,ridx,igrid
  real(dp),allocatable::sendbuf(:,:),recvbuf(:,:)
  integer,allocatable::dest_cpu(:)

#ifndef WITHOUTMPI
  nprops_ksec = 2 * twotondim + 2

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
           sendbuf(j, idx) = dble(xx1(emission(icpu,ilevel)%igrid(i) &
                & + ncoarse + (j-1)*ngridmax))
           sendbuf(twotondim + j, idx) = dble(xx2(emission(icpu,ilevel)%igrid(i) &
                & + ncoarse + (j-1)*ngridmax))
        end do
        sendbuf(2*twotondim + 1, idx) = dble(myid)
        sendbuf(2*twotondim + 2, idx) = dble(i)
     end do
  end do

  ! Exchange via ksection tree
  call ksection_exchange_dp(sendbuf, ntotal, dest_cpu, nprops_ksec, &
       & recvbuf, nrecv)

  ! Scatter received data (dp->int) to reception grids
  do i = 1, nrecv
     sender = nint(recvbuf(2*twotondim + 1, i))
     ridx   = nint(recvbuf(2*twotondim + 2, i))
     igrid  = reception(sender, ilevel)%igrid(ridx)
     do j = 1, twotondim
        xx1(igrid + ncoarse + (j-1)*ngridmax) = nint(recvbuf(j, i))
        xx2(igrid + ncoarse + (j-1)*ngridmax) = nint(recvbuf(twotondim + j, i))
     end do
  end do

  deallocate(sendbuf, dest_cpu, recvbuf)
#endif

end subroutine make_virtual_fine_int_pair_ksec
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
subroutine make_virtual_fine_dp_bulk(xx,ncols,ilevel)
  use amr_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer,intent(in)::ncols,ilevel
  real(dp),dimension(1:ncoarse+ngridmax*twotondim,1:ncols)::xx
  ! -------------------------------------------------------------------
  ! Bulk forward ghost zone exchange for multiple columns of a 2D array.
  ! Dispatches to ksec version or falls back to per-column calls.
  ! -------------------------------------------------------------------
  integer::ivar
  logical::use_ksec
  real(dp)::t1

  if(numbtot(1,ilevel)==0)return

  ! Auto-tune dispatch: K-Section bulk vs P2P fallback (comp 6: bulk_dp)
  use_ksec = .false.
  if(ordering=='ksection') then
     if(exchange_method=='ksection') then
        use_ksec = .true.
     else if(exchange_method=='auto') then
        use_ksec = (xchg_phase(6)==1) .or. &
                   (xchg_phase(6)==2 .and. xchg_chosen(6)==1) .or. &
                   (xchg_phase(6)==3 .and. xchg_chosen(6)==0)
     end if
  end if
  if(use_ksec) then
#ifndef WITHOUTMPI
     t1 = MPI_WTIME()
#endif
     call make_virtual_fine_dp_bulk_ksec(xx,ncols,ilevel)
#ifndef WITHOUTMPI
     if(exchange_method=='auto') call xchg_autotune_update(6, MPI_WTIME()-t1)
#endif
     return
  end if

  ! Fallback: call single-variable version for each column
#ifndef WITHOUTMPI
  if(ordering=='ksection') t1 = MPI_WTIME()
#endif
  do ivar=1,ncols
     call make_virtual_fine_dp(xx(1,ivar),ilevel)
  end do
#ifndef WITHOUTMPI
  if(ordering=='ksection' .and. exchange_method=='auto') &
     call xchg_autotune_update(6, MPI_WTIME()-t1)
#endif

end subroutine make_virtual_fine_dp_bulk
!################################################################
!################################################################
!################################################################
!################################################################
subroutine make_virtual_fine_dp_bulk_ksec(xx,ncols,ilevel)
  use amr_commons
  use ksection
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer,intent(in)::ncols,ilevel
  real(dp),dimension(1:ncoarse+ngridmax*twotondim,1:ncols)::xx
  ! -------------------------------------------------------------------
  ! Ksection-based bulk forward ghost zone exchange.
  ! Packs all ncols variables per emission grid into a single buffer,
  ! exchanges once via ksection tree, then scatters to reception grids.
  ! -------------------------------------------------------------------
  integer::icpu,i,j,iv,idx,ntotal,nrecv,nprops_ksec,sender,ridx,igrid
  real(dp),allocatable::sendbuf(:,:),recvbuf(:,:)
  integer,allocatable::dest_cpu(:)

#ifndef WITHOUTMPI
  nprops_ksec = ncols * twotondim + 2

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
        do iv = 1, ncols
           do j = 1, twotondim
              sendbuf((iv-1)*twotondim + j, idx) = &
                   & xx(emission(icpu,ilevel)%igrid(i) &
                   & + ncoarse + (j-1)*ngridmax, iv)
           end do
        end do
        sendbuf(ncols*twotondim + 1, idx) = dble(myid)
        sendbuf(ncols*twotondim + 2, idx) = dble(i)
     end do
  end do

  ! Exchange via ksection tree
  call ksection_exchange_dp(sendbuf, ntotal, dest_cpu, nprops_ksec, &
       & recvbuf, nrecv)

  ! Scatter received data to reception grids
  do i = 1, nrecv
     sender = nint(recvbuf(ncols*twotondim + 1, i))
     ridx   = nint(recvbuf(ncols*twotondim + 2, i))
     igrid  = reception(sender, ilevel)%igrid(ridx)
     do iv = 1, ncols
        do j = 1, twotondim
           xx(igrid + ncoarse + (j-1)*ngridmax, iv) = &
                & recvbuf((iv-1)*twotondim + j, i)
        end do
     end do
  end do

  deallocate(sendbuf, dest_cpu, recvbuf)
#endif

end subroutine make_virtual_fine_dp_bulk_ksec
!################################################################
!################################################################
!################################################################
!################################################################
subroutine make_virtual_reverse_dp_bulk(xx,ncols,ilevel)
  use amr_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer,intent(in)::ncols,ilevel
  real(dp),dimension(1:ncoarse+ngridmax*twotondim,1:ncols)::xx
  ! -------------------------------------------------------------------
  ! Bulk reverse ghost zone exchange for multiple columns of a 2D array.
  ! Dispatches to ksec version or falls back to per-column calls.
  ! -------------------------------------------------------------------
  integer::ivar
  logical::use_ksec
  real(dp)::t1

  if(numbtot(1,ilevel)==0)return

  ! Auto-tune dispatch: K-Section bulk vs P2P fallback (comp 7: bulk_rev_dp)
  use_ksec = .false.
  if(ordering=='ksection') then
     if(exchange_method=='ksection') then
        use_ksec = .true.
     else if(exchange_method=='auto') then
        use_ksec = (xchg_phase(7)==1) .or. &
                   (xchg_phase(7)==2 .and. xchg_chosen(7)==1) .or. &
                   (xchg_phase(7)==3 .and. xchg_chosen(7)==0)
     end if
  end if
  if(use_ksec) then
#ifndef WITHOUTMPI
     t1 = MPI_WTIME()
#endif
     call make_virtual_reverse_dp_bulk_ksec(xx,ncols,ilevel)
#ifndef WITHOUTMPI
     if(exchange_method=='auto') call xchg_autotune_update(7, MPI_WTIME()-t1)
#endif
     return
  end if

  ! Fallback: call single-variable version for each column
#ifndef WITHOUTMPI
  if(ordering=='ksection') t1 = MPI_WTIME()
#endif
  do ivar=1,ncols
     call make_virtual_reverse_dp(xx(1,ivar),ilevel)
  end do
#ifndef WITHOUTMPI
  if(ordering=='ksection' .and. exchange_method=='auto') &
     call xchg_autotune_update(7, MPI_WTIME()-t1)
#endif

end subroutine make_virtual_reverse_dp_bulk
!################################################################
!################################################################
!################################################################
!################################################################
subroutine make_virtual_reverse_dp_bulk_ksec(xx,ncols,ilevel)
  use amr_commons
  use ksection
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer,intent(in)::ncols,ilevel
  real(dp),dimension(1:ncoarse+ngridmax*twotondim,1:ncols)::xx
  ! -------------------------------------------------------------------
  ! Ksection-based bulk reverse ghost zone exchange.
  ! Packs all ncols variables per reception grid into a single buffer,
  ! exchanges once via ksection tree, then accumulates (+=) into
  ! owner's emission grids.
  ! -------------------------------------------------------------------
  integer::icpu,i,j,iv,idx,ntotal,nrecv,nprops_ksec,sender,eidx,igrid
  real(dp),allocatable::sendbuf(:,:),recvbuf(:,:)
  integer,allocatable::dest_cpu(:)

#ifndef WITHOUTMPI
  nprops_ksec = ncols * twotondim + 2

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
        do iv = 1, ncols
           do j = 1, twotondim
              sendbuf((iv-1)*twotondim + j, idx) = &
                   & xx(reception(icpu,ilevel)%igrid(i) &
                   & + ncoarse + (j-1)*ngridmax, iv)
           end do
        end do
        sendbuf(ncols*twotondim + 1, idx) = dble(myid)
        sendbuf(ncols*twotondim + 2, idx) = dble(i)
     end do
  end do

  ! Exchange via ksection tree
  call ksection_exchange_dp(sendbuf, ntotal, dest_cpu, nprops_ksec, &
       & recvbuf, nrecv)

  ! Accumulate received data into emission grids
  do i = 1, nrecv
     sender = nint(recvbuf(ncols*twotondim + 1, i))
     eidx   = nint(recvbuf(ncols*twotondim + 2, i))
     igrid  = emission(sender, ilevel)%igrid(eidx)
     do iv = 1, ncols
        do j = 1, twotondim
           xx(igrid + ncoarse + (j-1)*ngridmax, iv) = &
                & xx(igrid + ncoarse + (j-1)*ngridmax, iv) + &
                & recvbuf((iv-1)*twotondim + j, i)
        end do
     end do
  end do

  deallocate(sendbuf, dest_cpu, recvbuf)
#endif

end subroutine make_virtual_reverse_dp_bulk_ksec
!################################################################
!################################################################
!################################################################
!################################################################
subroutine xchg_autotune_update(icomp, elapsed)
  use amr_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer, intent(in) :: icomp
  real(dp), intent(in) :: elapsed
  ! -------------------------------------------------------------------
  ! Continuous adaptive auto-tune for exchange component icomp.
  !
  ! Phase 0: Test P2P for XCHG_NTRIAL calls
  ! Phase 1: Test K-Section for XCHG_NTRIAL calls, then decide
  ! Phase 2: Run with chosen method, track EMA.
  !          Every XCHG_NRECHECK calls, enter Phase 3 to probe.
  ! Phase 3: Probe alternative for XCHG_NPROBE calls.
  !          If alternative beats EMA by XCHG_SWITCH_MARGIN, switch.
  !          Return to Phase 2.
  !
  ! All phase transitions use MPI_ALLREDUCE(SUM) for consistency.
  ! -------------------------------------------------------------------
  real(dp) :: t_p2p_avg, t_ksec_avg, probe_avg
  real(dp) :: local_times(2), global_times(2)
  real(dp) :: local_val, global_val
  integer  :: ierr, old_chosen
  character(len=12) :: comp_names(NXCHG_COMP)
  character(len=8)  :: method_str

  if(exchange_method /= 'auto') return

  comp_names(1) = 'fine_dp'
  comp_names(2) = 'fine_int'
  comp_names(3) = 'reverse_dp'
  comp_names(4) = 'reverse_int'
  comp_names(5) = 'pair_int'
  comp_names(6) = 'bulk_dp'
  comp_names(7) = 'bulk_rev_dp'
  comp_names(8) = 'reserved'

  xchg_ncall(icomp) = xchg_ncall(icomp) + 1

  ! ===== Phase 0: Accumulate P2P timing =====
  if(xchg_phase(icomp) == 0) then
     xchg_time_p2p(icomp) = xchg_time_p2p(icomp) + elapsed
     if(xchg_ncall(icomp) >= XCHG_NTRIAL) then
        xchg_phase(icomp) = 1
        xchg_ncall(icomp) = 0
     end if

  ! ===== Phase 1: Accumulate K-Section timing, then decide =====
  else if(xchg_phase(icomp) == 1) then
     xchg_time_ksec(icomp) = xchg_time_ksec(icomp) + elapsed
     if(xchg_ncall(icomp) >= XCHG_NTRIAL) then
        t_p2p_avg  = xchg_time_p2p(icomp)  / XCHG_NTRIAL
        t_ksec_avg = xchg_time_ksec(icomp) / XCHG_NTRIAL

#ifndef WITHOUTMPI
        local_times(1) = t_p2p_avg
        local_times(2) = t_ksec_avg
        call MPI_ALLREDUCE(local_times, global_times, 2, &
             MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
        if(global_times(2) < global_times(1)) then
           xchg_chosen(icomp) = 1  ! K-Section faster
           xchg_ema(icomp) = global_times(2) / ncpu
        else
           xchg_chosen(icomp) = 0  ! P2P faster
           xchg_ema(icomp) = global_times(1) / ncpu
        end if
#else
        if(t_ksec_avg < t_p2p_avg) then
           xchg_chosen(icomp) = 1
           xchg_ema(icomp) = t_ksec_avg
        else
           xchg_chosen(icomp) = 0
           xchg_ema(icomp) = t_p2p_avg
        end if
#endif
        xchg_phase(icomp) = 2
        xchg_ncall(icomp) = 0
        xchg_run_count(icomp) = 0
        xchg_recheck_interval(icomp) = XCHG_NRECHECK

        if(myid == 1) then
           method_str = merge('ksection','P2P     ', xchg_chosen(icomp)==1)
           write(*,'(A,A,A,ES10.3,A,ES10.3,A,A)') &
                ' [xchg auto-tune] ', trim(comp_names(icomp)), &
                ': P2P=', t_p2p_avg*1d3, 'ms Ksec=', t_ksec_avg*1d3, &
                'ms -> ', trim(method_str)
        end if
     end if

  ! ===== Phase 2: Running with chosen method, tracking EMA =====
  else if(xchg_phase(icomp) == 2) then
     ! Update EMA locally (no MPI in hot path)
     if(xchg_ema(icomp) > 0.0) then
        xchg_ema(icomp) = XCHG_EMA_ALPHA * elapsed + &
             (1.0d0 - XCHG_EMA_ALPHA) * xchg_ema(icomp)
     else
        xchg_ema(icomp) = elapsed
     end if

     xchg_run_count(icomp) = xchg_run_count(icomp) + 1

     ! Time to re-check? Sync EMA and enter Phase 3
     if(xchg_run_count(icomp) >= xchg_recheck_interval(icomp)) then
#ifndef WITHOUTMPI
        local_val = xchg_ema(icomp)
        call MPI_ALLREDUCE(local_val, global_val, 1, &
             MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
        xchg_ema(icomp) = global_val / ncpu
#endif
        xchg_phase(icomp) = 3
        xchg_ncall(icomp) = 0
        xchg_probe_sum(icomp) = 0.0
     end if

  ! ===== Phase 3: Probing alternative method =====
  else if(xchg_phase(icomp) == 3) then
     xchg_probe_sum(icomp) = xchg_probe_sum(icomp) + elapsed
     if(xchg_ncall(icomp) >= XCHG_NPROBE) then
        probe_avg = xchg_probe_sum(icomp) / XCHG_NPROBE

        ! Sync probe result across ranks
#ifndef WITHOUTMPI
        local_val = probe_avg
        call MPI_ALLREDUCE(local_val, global_val, 1, &
             MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
        probe_avg = global_val / ncpu
#endif

        ! Compare probe vs current EMA: switch if probe is faster by margin
        old_chosen = xchg_chosen(icomp)
        if(probe_avg < xchg_ema(icomp) * (1.0d0 - XCHG_SWITCH_MARGIN)) then
           ! Switch to alternative method
           xchg_chosen(icomp) = 1 - xchg_chosen(icomp)
           xchg_ema(icomp) = probe_avg  ! reset EMA to probe value
           xchg_nswitch(icomp) = xchg_nswitch(icomp) + 1
           ! Backoff: double the re-check interval (max 16x base)
           xchg_recheck_interval(icomp) = min( &
                xchg_recheck_interval(icomp) * 2, XCHG_NRECHECK * 16)

           if(myid == 1) then
              method_str = merge('ksection','P2P     ', xchg_chosen(icomp)==1)
              write(*,'(A,A,A,A,A,ES10.3,A,ES10.3,A)') &
                   ' [xchg re-tune] ', trim(comp_names(icomp)), &
                   ': switch to ', trim(method_str), &
                   ' (probe=', probe_avg*1d3, &
                   'ms ema=', xchg_ema(icomp)*1d3, 'ms)'
           end if
        else
           ! No switch: reset interval to base (probe didn't beat current)
           xchg_recheck_interval(icomp) = XCHG_NRECHECK
        end if

        ! Return to Phase 2
        xchg_phase(icomp) = 2
        xchg_ncall(icomp) = 0
        xchg_run_count(icomp) = 0
     end if
  end if

end subroutine xchg_autotune_update
!################################################################
!################################################################
!################################################################
!################################################################
