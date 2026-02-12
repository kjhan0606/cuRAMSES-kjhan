!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine phi_fine_cg(ilevel,icount)
  use amr_commons
  use pm_commons
  use poisson_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer::ilevel,icount
  !=========================================================
  ! Iterative Poisson solver with Conjugate Gradient method
  ! to solve A x = b
  ! r  : stored in f(i,1)
  ! p  : stored in f(i,2)
  ! A p: stored in f(i,3)
  ! x  : stored in phi(i)
  ! b  : stored in rho(i)
  !=========================================================
  integer::i,idim,info,ind,iter,iskip,itermax,nx_loc
  integer::idx
  real(dp)::error,error_ini
  real(dp)::dx2,fourpi,scale,oneoversix,fact,fact2
  real(dp)::r2_old,alpha_cg,beta_cg
  real(kind=8)::r2,pAp,rhs_norm,r2_all,pAp_all,rhs_norm_all

  if(gravity_type>0)return
  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel


  ! Set constants
  dx2=(0.5D0**ilevel)**2
  nx_loc=icoarse_max-icoarse_min+1
  scale=boxlen/dble(nx_loc)
  fourpi=4.D0*ACOS(-1.0D0)*scale
  if(cosmo)fourpi=1.5D0*omega_m*aexp*scale
  oneoversix=1.0D0/dble(twondim)
  fact=oneoversix*fourpi*dx2
  fact2 = fact*fact

  !===============================
  ! Compute initial phi
  !===============================
   if(ilevel>levelmin)then
      call make_initial_phi(ilevel,icount)              ! Interpolate phi down
   else
      call make_multipole_phi(ilevel)            ! Fill up with simple initial guess
   endif
   call make_virtual_fine_dp(phi(1),ilevel)      ! Update boundaries
   call make_boundary_phi(ilevel)                ! Update physical boundaries

  !===============================
  ! Compute right-hand side norm
  !===============================
  rhs_norm=0.d0
!$omp parallel default(firstprivate) shared(active,rho) reduction(+:rhs_norm)
  do ind=1,twotondim
     iskip=ncoarse+(ind-1)*ngridmax
!$omp do
     do i=1,active(ilevel)%ngrid
        idx=active(ilevel)%igrid(i)+iskip
        rhs_norm=rhs_norm+fact2*(rho(idx)-rho_tot)*(rho(idx)-rho_tot)
     end do
!$omp end do nowait
  end do
!$omp end parallel
  ! Compute global norms
#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(rhs_norm,rhs_norm_all,1,MPI_DOUBLE_PRECISION,MPI_SUM,&
       & MPI_COMM_WORLD,info)
  rhs_norm=rhs_norm_all
#endif
  rhs_norm=DSQRT(rhs_norm/dble(twotondim*numbtot(1,ilevel)))

  !==============================================
  ! Compute r = b - Ax and store it into f(i,1)
  ! Also set p = r and store it into f(i,2)
  !==============================================
  call cmp_residual_cg(ilevel,icount)

  !====================================
  ! Main iteration loop
  !====================================
  iter=0; itermax=10000
  error=1.0D0; error_ini=1.0D0
  do while(error>epsilon*error_ini.and.iter<itermax)

     iter=iter+1

     !====================================
     ! Compute residual norm
     !====================================
     r2=0.0d0
!$omp parallel default(firstprivate) shared(active,f) reduction(+:r2)
     do ind=1,twotondim
        iskip=ncoarse+(ind-1)*ngridmax
!$omp do
        do i=1,active(ilevel)%ngrid
           idx=active(ilevel)%igrid(i)+iskip
           r2=r2+f(idx,1)*f(idx,1)
        end do
!$omp end do nowait
     end do
!$omp end parallel
     ! Compute global norm
#ifndef WITHOUTMPI
     call MPI_ALLREDUCE(r2,r2_all,1,MPI_DOUBLE_PRECISION,MPI_SUM,&
          & MPI_COMM_WORLD,info)
     r2=r2_all
#endif

     !====================================
     ! Compute beta factor
     !====================================
     if(iter==1)then
        beta_cg=0.
     else
        beta_cg=r2/r2_old
     end if
     r2_old=r2

     !====================================
     ! Recurrence on p
     !====================================
!$omp parallel default(firstprivate) shared(active,f)
     do ind=1,twotondim
        iskip=ncoarse+(ind-1)*ngridmax
!$omp do
        do i=1,active(ilevel)%ngrid
           idx=active(ilevel)%igrid(i)+iskip
           f(idx,2)=f(idx,1)+beta_cg*f(idx,2)
        end do
!$omp end do nowait
     end do
!$omp end parallel
     ! Update boundaries
     call make_virtual_fine_dp(f(1,2),ilevel)

     !==============================================
     ! Compute z = Ap and store it into f(i,3)
     !==============================================
     call cmp_Ap_cg(ilevel)

     !====================================
     ! Compute p.Ap scalar product
     !====================================
     pAp=0.0d0
!$omp parallel default(firstprivate) shared(active,f) reduction(+:pAp)
     do ind=1,twotondim
        iskip=ncoarse+(ind-1)*ngridmax
!$omp do
        do i=1,active(ilevel)%ngrid
           idx=active(ilevel)%igrid(i)+iskip
           pAp=pAp+f(idx,2)*f(idx,3)
        end do
!$omp end do nowait
     end do
!$omp end parallel
     ! Compute global sum
#ifndef WITHOUTMPI
     call MPI_ALLREDUCE(pAp,pAp_all,1,MPI_DOUBLE_PRECISION,MPI_SUM,&
          & MPI_COMM_WORLD,info)
     pAp=pAp_all
#endif

     !====================================
     ! Compute alpha factor
     !====================================
     alpha_cg = r2/pAp

     !====================================
     ! Recurrence on x
     !====================================
!$omp parallel default(firstprivate) shared(active,f,phi)
     do ind=1,twotondim
        iskip=ncoarse+(ind-1)*ngridmax
!$omp do
        do i=1,active(ilevel)%ngrid
           idx=active(ilevel)%igrid(i)+iskip
           phi(idx)=phi(idx)+alpha_cg*f(idx,2)
        end do
!$omp end do nowait
     end do
!$omp end parallel
     !====================================
     ! Recurrence on r
     !====================================
!$omp parallel default(firstprivate) shared(active,f)
     do ind=1,twotondim
        iskip=ncoarse+(ind-1)*ngridmax
!$omp do
        do i=1,active(ilevel)%ngrid
           idx=active(ilevel)%igrid(i)+iskip
           f(idx,1)=f(idx,1)-alpha_cg*f(idx,3)
        end do
!$omp end do nowait
     end do
!$omp end parallel
     ! Compute error
     error=DSQRT(r2/dble(twotondim*numbtot(1,ilevel)))
     if(iter==1)error_ini=error
     if(verbose)write(*,112)iter,error/rhs_norm,error/error_ini

  end do
  ! End main iteration loop

  if(myid==1)write(*,115)ilevel,iter,error/rhs_norm,error/error_ini
  if(iter >= itermax)then
     if(myid==1)write(*,*)'Poisson failed to converge...'
  end if

  ! Update boundaries
  call make_virtual_fine_dp(phi(1),ilevel)

111 format('   Entering phi_fine_cg for level ',I2)
112 format('   ==> Step=',i5,' Error=',2(1pe10.3,1x))
115 format('   ==> Level=',i5,' Step=',i5,' Error=',2(1pe10.3,1x))

end subroutine phi_fine_cg
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine cmp_residual_cg(ilevel,icount)
  use amr_commons
  use pm_commons
  use hydro_commons
  use poisson_commons
  implicit none
  integer::ilevel,icount
  !------------------------------------------------------------------
  ! This routine computes the residual for the Conjugate Gradient
  ! Poisson solver. The residual is stored in f(i,1).
  !------------------------------------------------------------------
  integer::i,idim,igrid,ngrid,ncache,ind,iskip,nx_loc
  integer::id1,id2,ig1,ig2,ih1,ih2
  real(dp)::dx2,fourpi,scale,oneoversix,fact
  integer,dimension(1:3,1:2,1:8)::iii,jjj

  !integer, parameter::nvector=32
  !integer ,dimension(1:nvector),save::ind_grid,ind_cell
  !integer ,dimension(1:nvector,0:twondim),save::igridn
  !integer ,dimension(1:nvector,1:ndim),save::ind_left,ind_right
  !real(dp),dimension(1:nvector,1:ndim),save::phig,phid
  !real(dp),dimension(1:nvector,1:twotondim,1:ndim),save::phi_left,phi_right
  !real(dp),dimension(1:nvector),save::residu

  integer ,dimension(1:nvector)::ind_grid,ind_cell
  integer ,dimension(1:nvector,0:twondim)::igridn
  integer ,dimension(1:nvector,1:ndim)::ind_left,ind_right
  real(dp),dimension(1:nvector,1:ndim)::phig,phid
  real(dp),dimension(1:nvector,1:twotondim,1:ndim)::phi_left,phi_right
  real(dp),dimension(1:nvector)::residu

  ! Set constants
  dx2=(0.5D0**ilevel)**2
  nx_loc=icoarse_max-icoarse_min+1
  scale=boxlen/dble(nx_loc)
  fourpi=4.D0*ACOS(-1.0D0)*scale
  if(cosmo)fourpi=1.5D0*omega_m*aexp*scale
  oneoversix=1.0D0/dble(twondim)
  fact=oneoversix*fourpi*dx2

  iii(1,1,1:8)=(/1,0,1,0,1,0,1,0/); jjj(1,1,1:8)=(/2,1,4,3,6,5,8,7/)
  iii(1,2,1:8)=(/0,2,0,2,0,2,0,2/); jjj(1,2,1:8)=(/2,1,4,3,6,5,8,7/)
  iii(2,1,1:8)=(/3,3,0,0,3,3,0,0/); jjj(2,1,1:8)=(/3,4,1,2,7,8,5,6/)
  iii(2,2,1:8)=(/0,0,4,4,0,0,4,4/); jjj(2,2,1:8)=(/3,4,1,2,7,8,5,6/)
  iii(3,1,1:8)=(/5,5,5,5,0,0,0,0/); jjj(3,1,1:8)=(/5,6,7,8,1,2,3,4/)
  iii(3,2,1:8)=(/0,0,0,0,6,6,6,6/); jjj(3,2,1:8)=(/5,6,7,8,1,2,3,4/)

  ! Loop over myid grids by vector sweeps
  ncache=active(ilevel)%ngrid
!$omp parallel do default(firstprivate) shared(active,nbor,son,phi,f)
  do igrid=1,ncache,nvector

     ! Gather nvector grids
     ngrid=MIN(nvector,ncache-igrid+1)
     do i=1,ngrid
        ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
     end do

     ! Gather neighboring grids
     do i=1,ngrid
        igridn(i,0)=ind_grid(i)
     end do
     do idim=1,ndim
        do i=1,ngrid
           ind_left (i,idim)=nbor(ind_grid(i),2*idim-1)
           ind_right(i,idim)=nbor(ind_grid(i),2*idim  )
           igridn(i,2*idim-1)=son(ind_left (i,idim))
           igridn(i,2*idim  )=son(ind_right(i,idim))
        end do
     end do

     ! Interpolate potential from upper level
     do idim=1,ndim
        call interpol_phi_cg(ind_left (1,idim),phi_left (1,1,idim),ngrid,ilevel,icount)
        call interpol_phi_cg(ind_right(1,idim),phi_right(1,1,idim),ngrid,ilevel,icount)
     end do

     ! Loop over cells
     do ind=1,twotondim
        ! Gather neighboring potential
        do idim=1,ndim
           id1=jjj(idim,1,ind); ig1=iii(idim,1,ind)
           ih1=ncoarse+(id1-1)*ngridmax
           do i=1,ngrid
              if(igridn(i,ig1)>0)then
                 phig(i,idim)=phi(igridn(i,ig1)+ih1)
              else
                 phig(i,idim)=phi_left(i,id1,idim)
              end if
           end do
           id2=jjj(idim,2,ind); ig2=iii(idim,2,ind)
           ih2=ncoarse+(id2-1)*ngridmax
           do i=1,ngrid
              if(igridn(i,ig2)>0)then
                 phid(i,idim)=phi(igridn(i,ig2)+ih2)
              else
                 phid(i,idim)=phi_right(i,id2,idim)
              end if
           end do
        end do

        ! Compute central cell index
        iskip=ncoarse+(ind-1)*ngridmax
        do i=1,ngrid
           ind_cell(i)=iskip+ind_grid(i)
        end do

        ! Compute residual using 6 neighbors potential
        do i=1,ngrid
           residu(i)=phi(ind_cell(i))
        end do
        do idim=1,ndim
           do i=1,ngrid
              residu(i)=residu(i)-oneoversix*(phig(i,idim)+phid(i,idim))
           end do
        end do
        do i=1,ngrid
           residu(i)=residu(i)+fact*(rho(ind_cell(i))-rho_tot)
        end do

        ! Store results in f(i,1)
        do i=1,ngrid
           f(ind_cell(i),1)=residu(i)
        end do

        ! Store results in f(i,2)
        do i=1,ngrid
           f(ind_cell(i),2)=residu(i)
        end do

     end do
     ! End loop over cells

  end do
  ! End loop over grids
!!$omp end parallel do

end subroutine cmp_residual_cg
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine cmp_Ap_cg(ilevel)
  use amr_commons
  use pm_commons
  use hydro_commons
  use poisson_commons
  implicit none
  integer::ilevel
  !------------------------------------------------------------------
  ! This routine computes Ap for the Conjugate Gradient
  ! Poisson Solver and store the result into f(i,3).
  !------------------------------------------------------------------
  integer::i,idim,igrid,ngrid,ncache,ind,iskip
  integer::id1,id2,ig1,ig2,ih1,ih2
  real(dp)::oneoversix
  integer,dimension(1:3,1:2,1:8)::iii,jjj

  !integer,parameter::nvector=32
  integer,dimension(1:nvector)::ind_grid,ind_cell
  integer,dimension(1:nvector,0:twondim)::igridn
  real(dp),dimension(1:nvector,1:ndim)::phig,phid
  real(dp),dimension(1:nvector)::residu

  ! Set constants
  oneoversix=1.0D0/dble(twondim)

  iii(1,1,1:8)=(/1,0,1,0,1,0,1,0/); jjj(1,1,1:8)=(/2,1,4,3,6,5,8,7/)
  iii(1,2,1:8)=(/0,2,0,2,0,2,0,2/); jjj(1,2,1:8)=(/2,1,4,3,6,5,8,7/)
  iii(2,1,1:8)=(/3,3,0,0,3,3,0,0/); jjj(2,1,1:8)=(/3,4,1,2,7,8,5,6/)
  iii(2,2,1:8)=(/0,0,4,4,0,0,4,4/); jjj(2,2,1:8)=(/3,4,1,2,7,8,5,6/)
  iii(3,1,1:8)=(/5,5,5,5,0,0,0,0/); jjj(3,1,1:8)=(/5,6,7,8,1,2,3,4/)
  iii(3,2,1:8)=(/0,0,0,0,6,6,6,6/); jjj(3,2,1:8)=(/5,6,7,8,1,2,3,4/)

  ! Loop over myid grids by vector sweeps
  ncache=active(ilevel)%ngrid
!$omp parallel do default(firstprivate) shared(active,nbor,son,f)
  do igrid=1,ncache,nvector

     ! Gather nvector grids
     ngrid=MIN(nvector,ncache-igrid+1)
     do i=1,ngrid
        ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
     end do

     ! Gather neighboring grids
     do i=1,ngrid
        igridn(i,0)=ind_grid(i)
     end do
     do idim=1,ndim
        do i=1,ngrid
           igridn(i,2*idim-1)=son(nbor(ind_grid(i),2*idim-1))
           igridn(i,2*idim  )=son(nbor(ind_grid(i),2*idim  ))
        end do
     end do

     ! Loop over cells
     do ind=1,twotondim

        ! Gather neighboring potential
        do idim=1,ndim
           id1=jjj(idim,1,ind); ig1=iii(idim,1,ind)
           ih1=ncoarse+(id1-1)*ngridmax
           do i=1,ngrid
              if(igridn(i,ig1)>0)then
                 phig(i,idim)=f(igridn(i,ig1)+ih1,2)
              else
                 phig(i,idim)=0.
              end if
           end do
           id2=jjj(idim,2,ind); ig2=iii(idim,2,ind)
           ih2=ncoarse+(id2-1)*ngridmax
           do i=1,ngrid
              if(igridn(i,ig2)>0)then
                 phid(i,idim)=f(igridn(i,ig2)+ih2,2)
              else
                 phid(i,idim)=0.
              end if
           end do
        end do

        ! Compute central cell index
        iskip=ncoarse+(ind-1)*ngridmax
        do i=1,ngrid
           ind_cell(i)=iskip+ind_grid(i)
        end do

        ! Compute Ap using neighbors potential
        do i=1,ngrid
           residu(i)=-f(ind_cell(i),2)
        end do
        do idim=1,ndim
           do i=1,ngrid
              residu(i)=residu(i)+oneoversix*(phig(i,idim)+phid(i,idim))
           end do
        end do
        ! Store results in f(i,3)
        do i=1,ngrid
           f(ind_cell(i),3)=residu(i)
        end do

     end do
     ! End loop over cells

  end do
  ! End loop over grids
!$omp end parallel do
end subroutine cmp_Ap_cg
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine make_initial_phi(ilevel,icount)
  use amr_commons
  use pm_commons
  use poisson_commons
  implicit none
  integer::ilevel,icount
  !
  !
  !
  !integer, parameter::nvector=32
  integer::igrid,ncache,i,ngrid,ind,iskip,idim,ibound
  integer ,dimension(1:nvector),save::ind_grid,ind_cell,ind_cell_father
  real(dp),dimension(1:nvector,1:twotondim),save::phi_int

  ! Loop over myid grids by vector sweeps
  ncache=active(ilevel)%ngrid
!$omp parallel do default(firstprivate) shared(active,phi,f)
  do igrid=1,ncache,nvector
     ! Gather nvector grids
     ngrid=MIN(nvector,ncache-igrid+1)
     do i=1,ngrid
        ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
     end do

     if(ilevel==1)then
        ! Loop over cells
        do ind=1,twotondim
           iskip=ncoarse+(ind-1)*ngridmax
           do i=1,ngrid
              ind_cell(i)=iskip+ind_grid(i)
           end do
           do i=1,ngrid
              phi(ind_cell(i))=0.0d0
           end do
           do idim=1,ndim
              do i=1,ngrid
                 f(ind_cell(i),idim)=0.0
              end do
           end do
        end do
        ! End loop over cells
     else
        ! Compute father cell index
        do i=1,ngrid
           ind_cell_father(i)=father(ind_grid(i))
        end do

        ! Interpolate
        call interpol_phi_cg(ind_cell_father,phi_int,ngrid,ilevel,icount)

        ! Loop over cells
        do ind=1,twotondim
           iskip=ncoarse+(ind-1)*ngridmax
           do i=1,ngrid
              ind_cell(i)=iskip+ind_grid(i)
           end do
           do i=1,ngrid
              phi(ind_cell(i))=phi_int(i,ind)
           end do
           do idim=1,ndim
              do i=1,ngrid
                 f(ind_cell(i),idim)=0.0
              end do
           end do
        end do
        ! End loop over cells
     end if

  end do
  ! End loop over grids

end subroutine make_initial_phi
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine make_multipole_phi(ilevel)
  use amr_commons
  use pm_commons
  use poisson_commons
  implicit none
  integer::ilevel
  !
  !
  !
  integer::ibound,boundary_dir,idim,inbor
  integer::i,ncache,ivar,igrid,ngrid,ind
  integer::iskip,iskip_ref,gdim,nx_loc,ix,iy,iz
  integer,dimension(1:nvector),save::ind_grid,ind_cell

  real(dp)::dx,dx_loc,scale,fourpi,boxlen2,eps,r2
  real(dp),dimension(1:3)::skip_loc
  real(dp),dimension(1:twotondim,1:3)::xc
  real(dp),dimension(1:nvector),save::rr,pp
  real(dp),dimension(1:nvector,1:ndim),save::xx
  real(dp),dimension(1:nvector,1:ndim),save::ff


  ! Mesh size at level ilevel
  dx=0.5D0**ilevel

  ! Rescaling factors
  nx_loc=(icoarse_max-icoarse_min+1)
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  if(ndim>0)skip_loc(1)=dble(icoarse_min)
  if(ndim>1)skip_loc(2)=dble(jcoarse_min)
  if(ndim>2)skip_loc(3)=dble(kcoarse_min)
  scale=boxlen/dble(nx_loc)
  dx_loc=dx*scale
  fourpi=4.D0*ACOS(-1.0D0)
  boxlen2=boxlen**2
  eps=dx_loc

  ! Set position of cell centers relative to grid center
  do ind=1,twotondim
     !iz=(ind-1)/4
     !iy=(ind-1-4*iz)/2
     iz=ISHFT((ind-1),-2)
     iy=ISHFT((ind-1-4*iz),-1)
     ix=(ind-1-2*iy-4*iz)
     if(ndim>0)xc(ind,1)=(dble(ix)-0.5D0)*dx
     if(ndim>1)xc(ind,2)=(dble(iy)-0.5D0)*dx
     if(ndim>2)xc(ind,3)=(dble(iz)-0.5D0)*dx
  end do

  ! Loop over myid grids by vector sweeps
  ncache=active(ilevel)%ngrid
!$omp parallel do default(firstprivate) shared(active,xg,skip_loc,multipole,phi)
  do igrid=1,ncache,nvector
     ! Gather nvector grids
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

        if(simple_boundary)then
           ! Compute cell center in code units
           do idim=1,ndim
               do i=1,ngrid
                  xx(i,idim)=xg(ind_grid(i),idim)+xc(ind,idim)
               end do
           end do

           ! Rescale position from code units to user units
           rr(1:ngrid)=0.0d0
           do idim=1,ndim
               do i=1,ngrid
                  xx(i,idim)=(xx(i,idim)-skip_loc(idim))*scale
                  rr(i)=rr(i)+(xx(i,idim)-multipole(idim+1)/multipole(1))**2
               end do
           end do

           do i=1,ngrid
               rr(i)=max(eps,sqrt(rr(i)))       ! Cutoff
           end do

           if(ngrid>0) call phi_ana(rr,pp,ngrid)

           ! Scatter variables
           do i=1,ngrid
               phi(ind_cell(i))=pp(i)/scale
           end do

        else
           do i=1,ngrid
               phi(ind_cell(i))=0d0
           end do
        endif

        ! End loop over cells
     end do

  end do
  ! End loop over grids

end subroutine make_multipole_phi

subroutine interpol_phi_cg(ind_cell,phi_int,ncell,ilevel,icount)
  use amr_commons
  use poisson_commons, only:phi,phi_old
  implicit none
  integer::ncell,ilevel,icount
!  integer, parameter::nvector=32
  integer ,dimension(1:nvector)::ind_cell
  real(dp),dimension(1:nvector,1:twotondim)::phi_int

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Routine for interpolation at level-boundaries. Interpolation is used for
  ! - boundary conditions for solving poisson equation at fine level
  ! - computing force (gradient_phi) at fine level for cells close to boundary
  ! Interpolation is performed in space (CIC) and - if adaptive timestepping is
  ! on -
  ! time (linear extrapolation of the change in phi during the last coarse step 
  ! onto the first fine step)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  integer ,dimension(1:nvector,1:twotondim)::nbors_father_grids
  integer ,dimension(1:nvector,1:threetondim)::nbors_father_cells
  integer::i,ind,indice,ind_average,ind_father
  real(dp)::dx,tfrac
  real(dp)::aa,bb,cc,dd,coeff,add
  integer,dimension(1:8,1:8)::ccc
  real(dp),dimension(1:8)::bbbb

  ! CIC method constants
  aa = 1.0D0/4.0D0**ndim
  bb = 3.0D0*aa
  cc = 9.0D0*aa
  dd = 27.D0*aa
  bbbb(:)  =(/aa ,bb ,bb ,cc ,bb ,cc ,cc ,dd/)

  ! Sampling positions in the 3x3x3 father cell cube
  ccc(:,1)=(/1 ,2 ,4 ,5 ,10,11,13,14/)
  ccc(:,2)=(/3 ,2 ,6 ,5 ,12,11,15,14/)
  ccc(:,3)=(/7 ,8 ,4 ,5 ,16,17,13,14/)
  ccc(:,4)=(/9 ,8 ,6 ,5 ,18,17,15,14/)
  ccc(:,5)=(/19,20,22,23,10,11,13,14/)
  ccc(:,6)=(/21,20,24,23,12,11,15,14/)
  ccc(:,7)=(/25,26,22,23,16,17,13,14/)
  ccc(:,8)=(/27,26,24,23,18,17,15,14/)

  if (icount .ne. 1 .and. icount .ne. 2)then
     write(*,*)'icount has bad value'
     call clean_stop
  endif

  ! Compute fraction of timesteps for interpolation
  if (dtold(ilevel-1)> 0)then
     tfrac=1.0*dtnew(ilevel)/dtold(ilevel-1)*(icount-1)
  else
     tfrac=0.
  end if

  ! Mesh size at level ilevel
  dx=0.5D0**ilevel
  call get3cubefather_cg(ind_cell,nbors_father_cells,nbors_father_grids,ncell,ilevel)

  ! Third order phi interpolation
  do ind=1,twotondim
     do i=1,ncell
        phi_int(i,ind)=0d0
     end do
     do ind_average=1,twotondim
        ind_father=ccc(ind_average,ind)
        coeff=bbbb(ind_average)
        do i=1,ncell
           indice=nbors_father_cells(i,ind_father)
           if (indice==0) then
              add=coeff*(phi(ind_cell(i))+(phi(ind_cell(i))-phi_old(ind_cell(i)))*tfrac)
              ! add=coeff*(-3d0/8d0*dx**2*boxlen*rho(ind_cell(i))+phi(ind_cell(i)))
           else
              add=coeff*(phi(indice)+(phi(indice)-phi_old(indice))*tfrac)
              ! add=coeff*(-3d0/8d0*dx**2*boxlen*rho(indice)+phi(indice)) 
           endif
           phi_int(i,ind)=phi_int(i,ind)+add
        end do
     end do
  end do

end subroutine interpol_phi_cg

subroutine get3cubefather_cg(ind_cell_father,nbors_father_cells,&
     &                    nbors_father_grids,ncell,ilevel)
! Comments by kjhan: This subroutine seems to be thread safe.
  use amr_commons
  implicit none
  integer::ncell,ilevel
!  integer, parameter::nvector=32
  integer,dimension(1:nvector)::ind_cell_father
  integer,dimension(1:nvector,1:threetondim)::nbors_father_cells
  integer,dimension(1:nvector,1:twotondim)::nbors_father_grids
  !------------------------------------------------------------------
  ! This subroutine determines the 3^ndim neighboring father cells 
  ! of the input father cell. According to the refinement rule, 
  ! they should be present anytime.
  !------------------------------------------------------------------
  integer::i,j,nxny,i1,j1,k1,ind,iok
  integer::i1min,i1max,j1min,j1max,k1min,k1max,ind_father
  integer,dimension(1:nvector)::ix,iy,iz,iix,iiy,iiz
  integer,dimension(1:nvector)::pos,ind_grid_father,ind_grid_ok
  integer,dimension(1:nvector,1:threetondim)::nbors_father_ok
  integer,dimension(1:nvector,1:twotondim)::nbors_grids_ok
  logical::oups

  nxny=nx*ny

  if(ilevel==1)then  ! Easy...

     oups=.false.
     do i=1,ncell
        if(ind_cell_father(i)>ncoarse)oups=.true.
     end do
     if(oups)then
        write(*,*)'get3cubefather_cg'
        write(*,*)'oupsssss !'
        call clean_stop
     endif

     do i=1,ncell
        iz(i)=(ind_cell_father(i)-1)/nxny
     end do
     do i=1,ncell
        iy(i)=(ind_cell_father(i)-1-iz(i)*nxny)/nx
     end do
     do i=1,ncell
        ix(i)=(ind_cell_father(i)-1-iy(i)*nx-iz(i)*nxny)
     end do


     i1min=0; i1max=0
     if(ndim > 0)i1max=2
     j1min=0; j1max=0
     if(ndim > 1)j1max=2
     k1min=0; k1max=0
     if(ndim > 2)k1max=2

     ! Loop over 3^ndim neighboring father cells
     do k1=k1min,k1max
        iiz=iz
        if(ndim > 2)then
           do i=1,ncell
              iiz(i)=iz(i)+k1-1
              if(iiz(i) < 0   )iiz(i)=nz-1
              if(iiz(i) > nz-1)iiz(i)=0
           end do
        end if
        do j1=j1min,j1max
           iiy=iy
           if(ndim > 1)then
              do i=1,ncell
                 iiy(i)=iy(i)+j1-1
                 if(iiy(i) < 0   )iiy(i)=ny-1
                 if(iiy(i) > ny-1)iiy(i)=0
              end do
           end if
           do i1=i1min,i1max
              iix=ix
              if(ndim > 0)then
                 do i=1,ncell
                    iix(i)=ix(i)+i1-1
                    if(iix(i) < 0   )iix(i)=nx-1
                    if(iix(i) > nx-1)iix(i)=0
                 end do
              end if
              ind_father=1+i1+3*j1+9*k1
              do i=1,ncell
                 nbors_father_cells(i,ind_father)=1 &
                      & +iix(i) &
                      & +iiy(i)*nx &
                      & +iiz(i)*nxny
              end do
           end do
        end do
     end do

     i1min=0; i1max=0
     if(ndim > 0)i1max=1
     j1min=0; j1max=0
     if(ndim > 1)j1max=1
     k1min=0; k1max=0
     if(ndim > 2)k1max=1
     ! Loop over 2^ndim neighboring father grids
     do k1=k1min,k1max
        iiz=iz
        if(ndim > 2)then
           do i=1,ncell
              iiz(i)=iz(i)+2*k1-1
              if(iiz(i) < 0   )iiz(i)=nz-1
              if(iiz(i) > nz-1)iiz(i)=0
           end do
        end if
        do j1=j1min,j1max
           iiy=iy
           if(ndim > 1)then
              do i=1,ncell
                 iiy(i)=iy(i)+2*j1-1
                 if(iiy(i) < 0   )iiy(i)=ny-1
                 if(iiy(i) > ny-1)iiy(i)=0
              end do
           end if
           do i1=i1min,i1max
              iix=ix
              if(ndim > 0)then
                 do i=1,ncell
                    iix(i)=ix(i)+2*i1-1
                    if(iix(i) < 0   )iix(i)=nx-1
                    if(iix(i) > nx-1)iix(i)=0
                 end do
              end if
              ind_father=1+i1+2*j1+4*k1
              do i=1,ncell
                 !nbors_father_grids(i,ind_father)=1 &
                 !     & +(iix(i)/2) &
                 !     & +(iiy(i)/2)*(nx/2) &
                 !     & +(iiz(i)/2)*(nxny/4)
                 nbors_father_grids(i,ind_father)=1 &
                      & +(ISHFT(iix(i), -1)) &
                      & +(ISHFT(iiy(i), -1))*(ISHFT(nx,-1)) &
                      & +(ISHFT(iiz(i),-1))*(ISHFT(nxny,-2))
              end do
           end do
        end do
     end do

  else    ! else, more complicated...

     ! Get father cell position in the grid
     do i=1,ncell
        pos(i)=(ind_cell_father(i)-ncoarse-1)/ngridmax+1
     end do
     ! Get father grid
     do i=1,ncell
        ind_grid_father(i)=ind_cell_father(i)-ncoarse-(pos(i)-1)*ngridmax
     end do

     ! Loop over position
     do ind=1,twotondim

        ! Select father cells that sit at position ind
        iok=0
        do i=1,ncell
           if(pos(i)==ind)then
              iok=iok+1
              ind_grid_ok(iok)=ind_grid_father(i)
           end if
        end do

        if(iok>0)&
        & call get3cubepos_cg(ind_grid_ok,ind,nbors_father_ok,nbors_grids_ok,iok)

        ! Store neighboring father cells for selected cells
        do j=1,threetondim
           iok=0
           do i=1,ncell
              if(pos(i)==ind)then
                 iok=iok+1
                 nbors_father_cells(i,j)=nbors_father_ok(iok,j)
              end if
           end do
        end do

        ! Store neighboring father grids for selected cells
        do j=1,twotondim
           iok=0
           do i=1,ncell
              if(pos(i)==ind)then
                 iok=iok+1
                 nbors_father_grids(i,j)=nbors_grids_ok(iok,j)
              end if
           end do
        end do

     end do

  end if

end subroutine get3cubefather_cg

subroutine get3cubepos_cg(ind_grid,ind,nbors_father_cells,nbors_father_grids,ng)
! Comments by kjhan: This subroutine seems to be thread safe.
  use amr_commons
  implicit none
  integer::ng,ind
  !integer, parameter::nvector=32
  integer,dimension(1:nvector)::ind_grid
  integer,dimension(1:nvector,1:threetondim)::nbors_father_cells
  integer,dimension(1:nvector,1:twotondim)::nbors_father_grids
  !--------------------------------------------------------------------
  ! This subroutine determines the 3^ndim neighboring father cells 
  ! of the input cell at position ind in grid ind_grid. According to 
  ! the refinements rules and since the input cell is refined, 
  ! they should be present anytime.
  !--------------------------------------------------------------------
  integer::i,j,iskip
  integer::ii,iimin,iimax
  integer::jj,jjmin,jjmax
  integer::kk,kkmin,kkmax
  integer::icell,igrid,inbor
  integer,dimension(1:8)::iii=(/1,2,1,2,1,2,1,2/)
  integer,dimension(1:8)::jjj=(/3,3,4,4,3,3,4,4/)
  integer,dimension(1:8)::kkk=(/5,5,5,5,6,6,6,6/)
  integer,dimension(1:27,1:8,1:3)::lll,mmm
  integer,dimension(1:nvector)::ind_grid1,ind_grid2,ind_grid3
  integer,dimension(1:nvector,1:twotondim)::nbors_grids

  lll=0; mmm=0
  ! -> ndim=1
  ! @ind =1
  lll(1:3,1,1)=(/2,1,1/)
  mmm(1:3,1,1)=(/2,1,2/)
  ! @ind =2
  lll(1:3,2,1)=(/1,1,2/)
  mmm(1:3,2,1)=(/1,2,1/)

  ! -> ndim=2
  ! @ind =1
  lll(1:9,1,2)=(/4,3,3,2,1,1,2,1,1/)
  mmm(1:9,1,2)=(/4,3,4,2,1,2,4,3,4/)
  ! @ind =2
  lll(1:9,2,2)=(/3,3,4,1,1,2,1,1,2/)
  mmm(1:9,2,2)=(/3,4,3,1,2,1,3,4,3/)
  ! @ind =3
  lll(1:9,3,2)=(/2,1,1,2,1,1,4,3,3/)
  mmm(1:9,3,2)=(/2,1,2,4,3,4,2,1,2/)
  ! @ind =4
  lll(1:9,4,2)=(/1,1,2,1,1,2,3,3,4/)
  mmm(1:9,4,2)=(/1,2,1,3,4,3,1,2,1/)

  ! -> ndim= 3
  ! @ind = 1
  lll(1:27,1,3)=(/8,7,7,6,5,5,6,5,5,4,3,3,2,1,1,2,1,1,4,3,3,2,1,1,2,1,1/)
  mmm(1:27,1,3)=(/8,7,8,6,5,6,8,7,8,4,3,4,2,1,2,4,3,4,8,7,8,6,5,6,8,7,8/)
  ! @ind = 2
  lll(1:27,2,3)=(/7,7,8,5,5,6,5,5,6,3,3,4,1,1,2,1,1,2,3,3,4,1,1,2,1,1,2/)
  mmm(1:27,2,3)=(/7,8,7,5,6,5,7,8,7,3,4,3,1,2,1,3,4,3,7,8,7,5,6,5,7,8,7/)
  ! @ind = 3
  lll(1:27,3,3)=(/6,5,5,6,5,5,8,7,7,2,1,1,2,1,1,4,3,3,2,1,1,2,1,1,4,3,3/)
  mmm(1:27,3,3)=(/6,5,6,8,7,8,6,5,6,2,1,2,4,3,4,2,1,2,6,5,6,8,7,8,6,5,6/)
  ! @ind = 4
  lll(1:27,4,3)=(/5,5,6,5,5,6,7,7,8,1,1,2,1,1,2,3,3,4,1,1,2,1,1,2,3,3,4/)
  mmm(1:27,4,3)=(/5,6,5,7,8,7,5,6,5,1,2,1,3,4,3,1,2,1,5,6,5,7,8,7,5,6,5/)
  ! @ind = 5
  lll(1:27,5,3)=(/4,3,3,2,1,1,2,1,1,4,3,3,2,1,1,2,1,1,8,7,7,6,5,5,6,5,5/)
  mmm(1:27,5,3)=(/4,3,4,2,1,2,4,3,4,8,7,8,6,5,6,8,7,8,4,3,4,2,1,2,4,3,4/)
  ! @ind = 6
  lll(1:27,6,3)=(/3,3,4,1,1,2,1,1,2,3,3,4,1,1,2,1,1,2,7,7,8,5,5,6,5,5,6/)
  mmm(1:27,6,3)=(/3,4,3,1,2,1,3,4,3,7,8,7,5,6,5,7,8,7,3,4,3,1,2,1,3,4,3/)
  ! @ind = 7
  lll(1:27,7,3)=(/2,1,1,2,1,1,4,3,3,2,1,1,2,1,1,4,3,3,6,5,5,6,5,5,8,7,7/)
  mmm(1:27,7,3)=(/2,1,2,4,3,4,2,1,2,6,5,6,8,7,8,6,5,6,2,1,2,4,3,4,2,1,2/)
  ! @ind = 8
  lll(1:27,8,3)=(/1,1,2,1,1,2,3,3,4,1,1,2,1,1,2,3,3,4,5,5,6,5,5,6,7,7,8/)
  mmm(1:27,8,3)=(/1,2,1,3,4,3,1,2,1,5,6,5,7,8,7,5,6,5,1,2,1,3,4,3,1,2,1/)

  iimin=0; iimax=0
  if(ndim>0)iimax=1
  jjmin=0; jjmax=0
  if(ndim>1)jjmax=1
  kkmin=0; kkmax=0
  if(ndim>2)kkmax=1

  do kk=kkmin,kkmax
     do i=1,ng
        ind_grid1(i)=ind_grid(i)
     end do
     if(kk>0)then
        inbor=kkk(ind)
        do i=1,ng
           if(ind_grid(i)>0)then
              ind_grid1(i)=son(nbor(ind_grid(i),inbor))
           endif
        end do
     end if

     do jj=jjmin,jjmax
        do i=1,ng
           ind_grid2(i)=ind_grid1(i)
        end do
        if(jj>0)then
           inbor=jjj(ind)
           do i=1,ng
              if(ind_grid1(i)>0)then
                 ind_grid2(i)=son(nbor(ind_grid1(i),inbor))
              endif
           end do
        end if

        do ii=iimin,iimax
           do i=1,ng
              ind_grid3(i)=ind_grid2(i)
           end do
           if(ii>0)then
              inbor=iii(ind)
              do i=1,ng
                 if(ind_grid2(i)>0)then
                    ind_grid3(i)=son(nbor(ind_grid2(i),inbor))
                 endif
              end do
           end if

           inbor=1+ii+2*jj+4*kk
           do i=1,ng
              nbors_grids(i,inbor)=ind_grid3(i)
           end do

        end do
     end do
  end do

  do j=1,twotondim
     do i=1,ng
        nbors_father_grids(i,j)=nbors_grids(i,j)
     end do
  end do

  do j=1,threetondim
     igrid=lll(j,ind,ndim)
     icell=mmm(j,ind,ndim)
     iskip=ncoarse+(icell-1)*ngridmax
     do i=1,ng
        if(nbors_grids(i,igrid)>0)then
           nbors_father_cells(i,j)=iskip+nbors_grids(i,igrid)
        else
           nbors_father_cells(i,j)=0
        endif
     end do
  end do

end subroutine get3cubepos_cg


