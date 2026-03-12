! ------------------------------------------------------------------------
! Multigrid Poisson solver for refined AMR levels
! ------------------------------------------------------------------------
! This file contains all MG-fine-level related routines
!
! Used variables:
!                       finest(AMR)level     coarse(MG)levels
!     -----------------------------------------------------------------
!     potential            phi            active_mg(myid,ilevel)%u(:,1)
!     physical RHS         rho            active_mg(myid,ilevel)%u(:,2)
!     residual             f(:,1)         active_mg(myid,ilevel)%u(:,3)
!     BC-modified RHS      f(:,2)                  N/A
!     mask                 f(:,3)         active_mg(myid,ilevel)%u(:,4)
!
! ------------------------------------------------------------------------


! ------------------------------------------------------------------------
! Mask restriction (top-down, OBSOLETE, UNUSED)
! ------------------------------------------------------------------------

subroutine restrict_mask_fine(ifinelevel,allmasked)
   use amr_commons
   use poisson_commons
   implicit none
   integer, intent(in)  :: ifinelevel
   logical, intent(out) :: allmasked

   integer :: ind_c_cell,ind_f_cell

   integer :: iskip_f_amr, iskip_c_amr, iskip_c_mg
   integer :: igrid_f_amr, igrid_c_amr, igrid_c_mg
   integer :: icell_f_amr, icell_c_amr, icell_c_mg

   real(dp) :: ngpmask
   real(dp) :: dtwotondim = (twotondim)

   integer :: icoarselevel
   icoarselevel=ifinelevel-1
   allmasked=.true.

   if(ifinelevel==1) return

   ! Loop over coarse cells of the coarse active comm for myid
   do ind_c_cell=1,twotondim
      iskip_c_amr=ncoarse+(ind_c_cell-1)*ngridmax
      iskip_c_mg =(ind_c_cell-1)*active_mg(myid,icoarselevel)%ngrid

      ! Loop over coarse grids
      do igrid_c_mg=1,active_mg(myid,icoarselevel)%ngrid
         igrid_c_amr=active_mg(myid,icoarselevel)%igrid(igrid_c_mg)
         icell_c_amr=iskip_c_amr+igrid_c_amr
         icell_c_mg =iskip_c_mg +igrid_c_mg

         if(son(icell_c_amr)==0) then
            ! Cell is not refined
            ngpmask      = -1.0d0
         else
            igrid_f_amr=son(icell_c_amr)
            ngpmask = 0.0d0
            do ind_f_cell=1,twotondim
               iskip_f_amr=ncoarse+(ind_f_cell-1)*ngridmax
               icell_f_amr=iskip_f_amr+igrid_f_amr
               ngpmask=ngpmask+f(icell_f_amr,3)
            end do
            ngpmask=ngpmask/dtwotondim
         end if
         ! Store cell mask
         active_mg(myid,icoarselevel)%u(icell_c_mg,4)=ngpmask
         allmasked=allmasked .and. (ngpmask<=0.0)
      end do
   end do

end subroutine restrict_mask_fine

! ------------------------------------------------------------------------
! Mask restriction (bottom-up)
! ------------------------------------------------------------------------

subroutine restrict_mask_fine_reverse(ifinelevel)
   use amr_commons
   use poisson_commons
   implicit none
   integer, intent(in) :: ifinelevel

   integer :: ind_c_cell, ind_f_cell, cpu_amr

   integer :: iskip_c_mg
   integer :: igrid_c_amr, igrid_c_mg
   integer :: icell_c_amr, icell_c_mg

   integer :: iskip_f_amr
   integer :: igrid_f_amr, igrid_f_mg
   integer :: icell_f_amr

   real(dp) :: ngpmask
   real(dp) :: dtwotondim = (twotondim)

   integer :: icoarselevel
   icoarselevel=ifinelevel-1

   ! Loop over fine cells of the myid active comm
   do ind_f_cell=1,twotondim
      iskip_f_amr=ncoarse+(ind_f_cell-1)*ngridmax

      ! Loop over fine grids of myid
      do igrid_f_mg=1,active(ifinelevel)%ngrid
         igrid_f_amr=active(ifinelevel)%igrid(igrid_f_mg)
         icell_f_amr=igrid_f_amr+iskip_f_amr

         ! Get coarse grid AMR index and CPU id
         icell_c_amr=father(igrid_f_amr)
         ind_c_cell=(icell_c_amr-ncoarse-1)/ngridmax+1
         igrid_c_amr=icell_c_amr-ncoarse-(ind_c_cell-1)*ngridmax
         cpu_amr=cpu_map(father(igrid_c_amr))

         ! Convert to MG index, get MG coarse cell id
         igrid_c_mg=lookup_mg(igrid_c_amr)
         iskip_c_mg=(ind_c_cell-1)*active_mg(cpu_amr,icoarselevel)%ngrid
         icell_c_mg=iskip_c_mg+igrid_c_mg

         ! Stack cell volume fraction in coarse cell
         ngpmask=(1d0+f(icell_f_amr,3))/2d0/dtwotondim
         active_mg(cpu_amr,icoarselevel)%u(icell_c_mg,4)=&
            active_mg(cpu_amr,icoarselevel)%u(icell_c_mg,4)+ngpmask
      end do
   end do
end subroutine restrict_mask_fine_reverse

! ------------------------------------------------------------------------
! Residual computation
! ------------------------------------------------------------------------

subroutine cmp_residual_mg_fine(ilevel)
   ! Computes the residual the fine (AMR) level, and stores it into f(:,1)
   use amr_commons
   use poisson_commons
   implicit none
   integer, intent(in) :: ilevel

   integer, dimension(1:3,1:2,1:8) :: iii, jjj

   real(dp) :: dx, oneoverdx2, phi_c, nb_sum
   integer  :: ngrid
   integer  :: ind, igrid_mg, idim, inbor
   integer  :: igrid_amr, icell_amr, iskip_amr
   integer  :: igshift, igrid_nbor_amr, icell_nbor_amr

   real(dp) :: dtwondim = (twondim)

   ! Set constants
   dx  = 0.5d0**ilevel
   oneoverdx2 = 1.0d0/(dx*dx)

   iii(1,1,1:8)=(/1,0,1,0,1,0,1,0/); jjj(1,1,1:8)=(/2,1,4,3,6,5,8,7/)
   iii(1,2,1:8)=(/0,2,0,2,0,2,0,2/); jjj(1,2,1:8)=(/2,1,4,3,6,5,8,7/)
   iii(2,1,1:8)=(/3,3,0,0,3,3,0,0/); jjj(2,1,1:8)=(/3,4,1,2,7,8,5,6/)
   iii(2,2,1:8)=(/0,0,4,4,0,0,4,4/); jjj(2,2,1:8)=(/3,4,1,2,7,8,5,6/)
   iii(3,1,1:8)=(/5,5,5,5,0,0,0,0/); jjj(3,1,1:8)=(/5,6,7,8,1,2,3,4/)
   iii(3,2,1:8)=(/0,0,0,0,6,6,6,6/); jjj(3,2,1:8)=(/5,6,7,8,1,2,3,4/)

   ngrid=active(ilevel)%ngrid

   ! Loop over cells
!$omp parallel default(firstprivate) shared(active,flag2,son,nbor,phi,f)
   do ind=1,twotondim
      iskip_amr = ncoarse+(ind-1)*ngridmax

      ! Loop over active grids
!$omp do
      do igrid_mg=1,ngrid
         igrid_amr = active(ilevel)%igrid(igrid_mg)
         icell_amr = iskip_amr + igrid_amr

         phi_c = phi(icell_amr)           ! Value of potential on center cell
         nb_sum=0.0d0                     ! Sum of phi on neighbors

         ! SCAN FLAG TEST
         if(flag2(icell_amr)/ngridmax==0) then
            do inbor=1,2
               do idim=1,ndim
                  ! Get neighbor grid
                  igshift = iii(idim,inbor,ind)
                  if(igshift==0) then
                     igrid_nbor_amr = igrid_amr
                  else
                     igrid_nbor_amr = son(nbor(igrid_amr,igshift))
                  end if
                  icell_nbor_amr = igrid_nbor_amr + &
                      (ncoarse + (jjj(idim,inbor,ind)-1)*ngridmax)
                  nb_sum = nb_sum + phi(icell_nbor_amr)
               end do
            end do
         else ! PERFORM SCAN
            if(f(icell_amr,3)<=0.0) then
               f(icell_amr,1)=0.0d0
               cycle
            end if
            do idim=1,ndim
               do inbor=1,2
                  ! Get neighbor grid
                  igshift = iii(idim,inbor,ind)
                  if(igshift==0) then
                     igrid_nbor_amr = igrid_amr
                  else
                     igrid_nbor_amr = son(nbor(igrid_amr,igshift))
                  end if

                  if(igrid_nbor_amr==0) then
                     ! No neighbor cell !
                     ! Virtual phi value on unrefined neighbor cell :
                     ! -phi_c/mask_c
                     ! (simulates mask=-1.0 for the nonexistent refined cell)
                     nb_sum = nb_sum - phi_c/f(icell_amr,3)
                  else
                     ! Fetch neighbor cell
                     icell_nbor_amr = igrid_nbor_amr + &
                         (ncoarse + (jjj(idim,inbor,ind)-1)*ngridmax)
                     if(f(icell_nbor_amr,3)<=0.0) then
                        ! Neighbor cell is masked :
                        ! compute its virtual phi with the mask
                        nb_sum = nb_sum + &
                            phi_c*(f(icell_nbor_amr,3)/f(icell_amr,3))
                     else
                        ! Neighbor cell is active, use its true potential
                        nb_sum = nb_sum + phi(icell_nbor_amr)
                     end if
                  end if
               end do
            end do
         end if ! END SCAN TEST

         ! Store ***MINUS THE RESIDUAL*** in f(:,1), using BC-modified RHS
         f(icell_amr,1) = -oneoverdx2*( nb_sum - dtwondim*phi_c )+f(icell_amr,2)
      end do
!$omp end do
   end do
!$omp end parallel

end subroutine cmp_residual_mg_fine

! ##################################################################
! ##################################################################

subroutine cmp_residual_norm2_fine(ilevel, norm2)
   use amr_commons
   use poisson_commons
   implicit none

   integer,  intent(in)  :: ilevel
   real(kind=8), intent(out) :: norm2

   real(kind=8) :: dx2
   integer  :: ngrid
   integer  :: ind, igrid_mg
   integer  :: igrid_amr, icell_amr, iskip_amr

   ! Set constants
   dx2  = (0.5d0**ilevel)**ndim
   ngrid=active(ilevel)%ngrid

   norm2 = 0.0d0
   ! Loop over cells
!$omp parallel default(firstprivate) shared(active,f) reduction(+:norm2)
   do ind=1,twotondim
      iskip_amr = ncoarse+(ind-1)*ngridmax
      ! Loop over active grids
!$omp do
      do igrid_mg=1,ngrid
         igrid_amr = active(ilevel)%igrid(igrid_mg)
         icell_amr = iskip_amr + igrid_amr
         if(f(icell_amr,3)<=0.0) then      ! Do not count masked cells
            cycle
         end if
         norm2 = norm2 + f(icell_amr,1)**2
      end do
!$omp end do nowait
   end do
!$omp end parallel
   norm2 = dx2*norm2

end subroutine cmp_residual_norm2_fine

subroutine cmp_ivar_norm2_fine(ilevel, ivar, norm2)
   use amr_commons
   use poisson_commons
   implicit none

   integer,  intent(in)  :: ilevel, ivar
   real(kind=8), intent(out) :: norm2

   real(kind=8) :: dx2
   integer  :: ngrid
   integer  :: ind, igrid_mg
   integer  :: igrid_amr, icell_amr, iskip_amr

   ! Set constants
   dx2  = (0.5d0**ilevel)**ndim
   ngrid=active(ilevel)%ngrid

   norm2 = 0.0d0
   ! Loop over cells
   do ind=1,twotondim
      iskip_amr = ncoarse+(ind-1)*ngridmax
      ! Loop over active grids
      do igrid_mg=1,ngrid
         igrid_amr = active(ilevel)%igrid(igrid_mg)
         icell_amr = iskip_amr + igrid_amr
         if(f(icell_amr,3)<=0.0) then      ! Do not count masked cells
            cycle
         end if
         if(ivar>0)then
            norm2 = norm2 + f(icell_amr,ivar)**2
         else
            norm2 = norm2 + phi(icell_amr)**2
         endif
      end do
   end do
   norm2 = dx2*norm2

end subroutine cmp_ivar_norm2_fine

! ------------------------------------------------------------------------
! Gauss-Seidel smoothing
! ------------------------------------------------------------------------

subroutine gauss_seidel_mg_fine(ilevel,redstep)
   use amr_commons
   use pm_commons
   use poisson_commons
   implicit none
   integer, intent(in) :: ilevel
   logical, intent(in) :: redstep

   integer, dimension(1:3,1:2,1:8) :: iii, jjj
   integer, dimension(1:3,1:4)     :: ired, iblack

   real(dp) :: dx2, nb_sum, weight
   integer  :: ngrid
   integer  :: ind, ind0, igrid_mg, idim, inbor
   integer  :: igrid_amr, icell_amr, iskip_amr
   integer  :: igshift, igrid_nbor_amr, icell_nbor_amr

   real(dp) :: dtwondim = (twondim)

   ! Set constants
   dx2  = (0.5d0**ilevel)**2

   ired  (1,1:4)=(/1,0,0,0/)
   iblack(1,1:4)=(/2,0,0,0/)
   ired  (2,1:4)=(/1,4,0,0/)
   iblack(2,1:4)=(/2,3,0,0/)
   ired  (3,1:4)=(/1,4,6,7/)
   iblack(3,1:4)=(/2,3,5,8/)

   iii(1,1,1:8)=(/1,0,1,0,1,0,1,0/); jjj(1,1,1:8)=(/2,1,4,3,6,5,8,7/)
   iii(1,2,1:8)=(/0,2,0,2,0,2,0,2/); jjj(1,2,1:8)=(/2,1,4,3,6,5,8,7/)
   iii(2,1,1:8)=(/3,3,0,0,3,3,0,0/); jjj(2,1,1:8)=(/3,4,1,2,7,8,5,6/)
   iii(2,2,1:8)=(/0,0,4,4,0,0,4,4/); jjj(2,2,1:8)=(/3,4,1,2,7,8,5,6/)
   iii(3,1,1:8)=(/5,5,5,5,0,0,0,0/); jjj(3,1,1:8)=(/5,6,7,8,1,2,3,4/)
   iii(3,2,1:8)=(/0,0,0,0,6,6,6,6/); jjj(3,2,1:8)=(/5,6,7,8,1,2,3,4/)

   ngrid=active(ilevel)%ngrid
   ! Loop over cells, with red/black ordering
!$omp parallel default(firstprivate) shared(active,flag2,son,nbor,phi,f,safe_mode)
   do ind0=1,twotondim/2      ! Only half of the cells for a red or black sweep
      if(redstep) then
         ind = ired  (ndim,ind0)
      else
         ind = iblack(ndim,ind0)
      end if

      iskip_amr = ncoarse+(ind-1)*ngridmax

      ! Loop over active grids
!$omp do
      do igrid_mg=1,ngrid
         igrid_amr = active(ilevel)%igrid(igrid_mg)
         icell_amr = iskip_amr + igrid_amr

         nb_sum=0.0d0                       ! Sum of phi on neighbors

         ! Read scan flag
         if(flag2(icell_amr)/ngridmax==0) then
            ! Use max-speed "dumb" Gauss-Seidel for "inner" cells
            ! Those cells are active, have all their neighbors active
            ! and all neighbors are in the AMR+MG trees
            do inbor=1,2
               do idim=1,ndim
                  ! Get neighbor grid shift
                  igshift = iii(idim,inbor,ind)
                  ! Get neighbor grid
                  if(igshift==0) then
                     igrid_nbor_amr = igrid_amr
                  else
                     igrid_nbor_amr = son(nbor(igrid_amr,igshift))
                  end if
                  icell_nbor_amr = igrid_nbor_amr + &
                      (ncoarse + (jjj(idim,inbor,ind)-1)*ngridmax)
                  nb_sum = nb_sum + phi(icell_nbor_amr)
               end do
            end do
            ! Update the potential, solving for potential on icell_amr
            phi(icell_amr) = (nb_sum - dx2*f(icell_amr,2)) / dtwondim
         else
            ! Use the finer "solve" Gauss-Seidel near boundaries,
            ! with all necessary checks
            if (f(icell_amr,3)<=0.0) cycle
            if (safe_mode(ilevel) .and. f(icell_amr,3)<1.0) cycle

            weight=0.0d0 ! Central weight for "Solve G-S"
            do inbor=1,2
               do idim=1,ndim
                  ! Get neighbor grid shift
                  igshift = iii(idim,inbor,ind)

                  ! Get neighbor grid
                  if(igshift==0) then
                     igrid_nbor_amr = igrid_amr
                  else
                     igrid_nbor_amr = son(nbor(igrid_amr,igshift))
                  end if

                  if(igrid_nbor_amr==0) then
                     ! No neighbor cell,
                     ! set mask=-1 on nonexistent neighbor cell
                     weight = weight - 1.0d0/f(icell_amr,3)
                  else
                     ! Fetch neighbor cell
                     icell_nbor_amr = igrid_nbor_amr + &
                         (ncoarse + (jjj(idim,inbor,ind)-1)*ngridmax)
                     if(f(icell_nbor_amr,3)<=0.0) then
                        ! Neighbor cell is masked
                        weight = weight + f(icell_nbor_amr,3)/f(icell_amr,3)
                     else
                        ! Neighbor cell is active, increment neighbor sum
                        nb_sum = nb_sum + phi(icell_nbor_amr)
                     end if
                  end if
               end do
            end do
            ! Update the potential, solving for potential on icell_amr
            phi(icell_amr) = (nb_sum - dx2*f(icell_amr,2)) / (dtwondim - weight)
         end if
      end do
!$omp end do nowait
   end do
!$omp end parallel
end subroutine gauss_seidel_mg_fine

! ------------------------------------------------------------------------
! Residual restriction (top-down, OBSOLETE, UNUSED)
! ------------------------------------------------------------------------

subroutine restrict_residual_fine(ifinelevel)
   ! Restrict fine (AMR) residual at level ifinelevel using injection
   ! into coarser residual at level ifinelevel-1
   ! Restricted residual is stored into the RHS at the coarser level
   use amr_commons
   use pm_commons
   use poisson_commons
   implicit none
   integer, intent(in) :: ifinelevel

   real(dp) :: val
   real(dp) :: dtwotondim = (twotondim)

   integer  :: icoarselevel
   integer  :: ngrid_c, ind_c, iskip_c_amr, iskip_c_mg
   integer  :: igrid_c_amr, icell_c_amr, icell_c_mg, igrid_c_mg
   integer  :: ind_f, igrid_f_amr, iskip_f_amr, icell_f_amr

   icoarselevel=ifinelevel-1

   ! Loop over coarse MG cells
   ngrid_c=active_mg(myid,icoarselevel)%ngrid
   do ind_c=1,twotondim
      iskip_c_amr = ncoarse + (ind_c-1)*ngridmax
      iskip_c_mg  = (ind_c-1)*ngrid_c

      do igrid_c_mg=1,ngrid_c
         igrid_c_amr = active_mg(myid,icoarselevel)%igrid(igrid_c_mg)
         icell_c_amr = igrid_c_amr + iskip_c_amr
         icell_c_mg  = igrid_c_mg  + iskip_c_mg

         ! Get AMR child grid
         igrid_f_amr = son(icell_c_amr)
         if(igrid_f_amr==0) then
            ! Nullify residual (coarser RHS)
            active_mg(myid,icoarselevel)%u(icell_c_mg,2) = 0.0d0
            cycle
         end if

         val = 0.0d0
         ! Loop over child (fine MG) cells
         do ind_f=1,twotondim
            iskip_f_amr = ncoarse + (ind_f-1)*ngridmax
            icell_f_amr = igrid_f_amr + iskip_f_amr

            if (f(icell_f_amr,3)<=0.0) cycle
            val = val + f(icell_f_amr,1)
         end do
         ! Store restricted residual into RHS of coarse level
         active_mg(myid,icoarselevel)%u(icell_c_mg,2) = val/dtwotondim
      end do
   end do
end subroutine restrict_residual_fine


! ------------------------------------------------------------------------
! Residual restriction (bottom-up)
! ------------------------------------------------------------------------

subroutine restrict_residual_fine_reverse(ifinelevel)
   use amr_commons
   use poisson_commons
   implicit none
   integer, intent(in) :: ifinelevel

   integer :: ind_c_cell, ind_f_cell, cpu_amr

   integer :: iskip_c_mg
   integer :: igrid_c_amr, igrid_c_mg
   integer :: icell_c_amr, icell_c_mg

   integer :: iskip_f_amr
   integer :: igrid_f_amr, igrid_f_mg
   integer :: icell_f_amr

   real(dp) :: res
   real(dp) :: dtwotondim = (twotondim)

   integer :: icoarselevel
   icoarselevel=ifinelevel-1

   ! Loop over fine cells of the myid active comm
!$omp parallel default(firstprivate) shared(active,active_mg,f,cpu_map,father,lookup_mg)
   do ind_f_cell=1,twotondim
      iskip_f_amr=ncoarse+(ind_f_cell-1)*ngridmax

      ! Loop over fine grids of myid
!$omp do
      do igrid_f_mg=1,active(ifinelevel)%ngrid
         igrid_f_amr=active(ifinelevel)%igrid(igrid_f_mg)
         icell_f_amr=igrid_f_amr+iskip_f_amr
         ! Is fine cell masked?
         if(f(icell_f_amr,3)<=0d0) cycle

         ! Get coarse grid AMR index and CPU id
         icell_c_amr=father(igrid_f_amr)
         ind_c_cell=(icell_c_amr-ncoarse-1)/ngridmax+1
         igrid_c_amr=icell_c_amr-ncoarse-(ind_c_cell-1)*ngridmax
         cpu_amr=cpu_map(father(igrid_c_amr))

         ! Convert to MG index, get MG coarse cell id
         igrid_c_mg=lookup_mg(igrid_c_amr)
         iskip_c_mg=(ind_c_cell-1)*active_mg(cpu_amr,icoarselevel)%ngrid
         icell_c_mg=iskip_c_mg+igrid_c_mg

         ! Is coarse cell masked?
         if(active_mg(cpu_amr,icoarselevel)%u(icell_c_mg,4)<=0d0) cycle

         ! Stack fine cell residual in coarse cell rhs
         res=f(icell_f_amr,1)/dtwotondim
         active_mg(cpu_amr,icoarselevel)%u(icell_c_mg,2)=&
            active_mg(cpu_amr,icoarselevel)%u(icell_c_mg,2)+res
      end do
!$omp end do nowait
   end do
!$omp end parallel
end subroutine restrict_residual_fine_reverse

! ------------------------------------------------------------------------
! Interpolation and correction
! ------------------------------------------------------------------------

subroutine interpolate_and_correct_fine(ifinelevel)
   use amr_commons
   use poisson_commons
   implicit none
   integer, intent(in) :: ifinelevel

   integer  :: i, ind_father, ind_average, ind_f, iskip_f_amr
   integer  :: ngrid_f, istart, nbatch
   integer  :: icell_c_amr, igrid_c_amr, igrid_c_mg, icell_c_mg
   integer  :: icoarselevel, ind_c, cpu_amr

   real(dp) :: a, b, c, d, coeff
   real(dp), dimension(1:8)     :: bbb
   integer,  dimension(1:8,1:8) :: ccc

   integer,  dimension(1:nvector)                :: igrid_f_amr, icell_amr
   integer,  dimension(1:nvector,1:threetondim)  :: nbors_father_cells
   integer,  dimension(1:nvector,1:twotondim)    :: nbors_father_grids
   real(dp), dimension(1:nvector)                :: corr

   ! Local constants
   a = 1.0D0/4.0D0**ndim
   b = 3.0D0*a
   c = 9.0D0*a
   d = 27.D0*a
   icoarselevel=ifinelevel-1

   bbb(:)  =(/a ,b ,b ,c ,b ,c ,c ,d/)

   ccc(:,1)=(/1 ,2 ,4 ,5 ,10,11,13,14/)
   ccc(:,2)=(/3 ,2 ,6 ,5 ,12,11,15,14/)
   ccc(:,3)=(/7 ,8 ,4 ,5 ,16,17,13,14/)
   ccc(:,4)=(/9 ,8 ,6 ,5 ,18,17,15,14/)
   ccc(:,5)=(/19,20,22,23,10,11,13,14/)
   ccc(:,6)=(/21,20,24,23,12,11,15,14/)
   ccc(:,7)=(/25,26,22,23,16,17,13,14/)
   ccc(:,8)=(/27,26,24,23,18,17,15,14/)

   ! Loop over fine grids by vector sweeps
   ngrid_f=active(ifinelevel)%ngrid
!$omp parallel do default(firstprivate) shared(active,active_mg,father,cpu_map,lookup_mg,f,phi)
   do istart=1,ngrid_f,nvector

      ! Gather nvector grids
      nbatch=MIN(nvector,ngrid_f-istart+1)
      do i=1,nbatch
         igrid_f_amr(i)=active(ifinelevel)%igrid(istart+i-1)
      end do

      ! Compute father (coarse) cell index
      do i=1,nbatch
         icell_amr(i)=father(igrid_f_amr(i))
      end do

      ! Gather 3x3x3 neighboring parent cells
      call get3cubefather_fine_mg(icell_amr,nbors_father_cells,nbors_father_grids, &
              nbatch,ifinelevel)

      ! Update solution for fine grid cells
      do ind_f=1,twotondim
         iskip_f_amr = ncoarse+(ind_f-1)*ngridmax

         do i=1,nbatch
            ! Compute fine cell indices
            icell_amr(i) = iskip_f_amr + igrid_f_amr(i)
         end do
         corr=0.0d0

         ! Loop over relevant parent cells
         do ind_average=1,twotondim
            ind_father = ccc(ind_average,ind_f)
            coeff      = bbb(ind_average)
            do i=1,nbatch
               if(f(icell_amr(i),3)<=0.0) then
                  corr(i)=0.0d0        ! Fine cell is masked : no correction
                  cycle
               end if
               icell_c_amr = nbors_father_cells(i,ind_father)
               ind_c       = (icell_c_amr-ncoarse-1)/ngridmax + 1
               igrid_c_amr = icell_c_amr - ncoarse - (ind_c-1)*ngridmax
               cpu_amr     = cpu_map(father(igrid_c_amr))
               igrid_c_mg  = lookup_mg(igrid_c_amr)
               if(igrid_c_mg<=0) cycle

               icell_c_mg=(ind_c-1)*active_mg(cpu_amr,icoarselevel)%ngrid+igrid_c_mg
               corr(i)=corr(i)+coeff*active_mg(cpu_amr,icoarselevel)%u(icell_c_mg,1)
            end do
         end do

         ! Correct potential
         do i=1,nbatch
            phi(icell_amr(i))=phi(icell_amr(i))+corr(i)
         end do

      end do
      ! End loop over cells

   end do
!$omp end parallel do
   ! End loop over grids
end subroutine interpolate_and_correct_fine


! ------------------------------------------------------------------------
! Flag setting
! ------------------------------------------------------------------------

subroutine set_scan_flag_fine(ilevel)
   use amr_commons
   use poisson_commons
   implicit none

   integer, intent(in) :: ilevel

   integer :: ind, ngrid, scan_flag
   integer :: igrid_mg, inbor, idim, igshift
   integer :: igrid_amr, igrid_nbor_amr

   integer :: iskip_amr, icell_amr, icell_nbor_amr

   integer, dimension(1:3,1:2,1:8) :: iii, jjj

   iii(1,1,1:8)=(/1,0,1,0,1,0,1,0/); jjj(1,1,1:8)=(/2,1,4,3,6,5,8,7/)
   iii(1,2,1:8)=(/0,2,0,2,0,2,0,2/); jjj(1,2,1:8)=(/2,1,4,3,6,5,8,7/)
   iii(2,1,1:8)=(/3,3,0,0,3,3,0,0/); jjj(2,1,1:8)=(/3,4,1,2,7,8,5,6/)
   iii(2,2,1:8)=(/0,0,4,4,0,0,4,4/); jjj(2,2,1:8)=(/3,4,1,2,7,8,5,6/)
   iii(3,1,1:8)=(/5,5,5,5,0,0,0,0/); jjj(3,1,1:8)=(/5,6,7,8,1,2,3,4/)
   iii(3,2,1:8)=(/0,0,0,0,6,6,6,6/); jjj(3,2,1:8)=(/5,6,7,8,1,2,3,4/)

   ngrid = active(ilevel)%ngrid

   ! Loop over cells and set fine SCAN flag
   do ind=1,twotondim
      iskip_amr = ncoarse+(ind-1)*ngridmax
      do igrid_mg=1,ngrid
         igrid_amr = active(ilevel)%igrid(igrid_mg)
         icell_amr = iskip_amr + igrid_amr

         if(f(icell_amr,3)==1.0) then
            scan_flag=0       ! Init flag to 'no scan needed'
            scan_flag_loop: do inbor=1,2
               do idim=1,ndim
                  igshift = iii(idim,inbor,ind)
                  if(igshift==0) then
                     igrid_nbor_amr = igrid_amr
                  else
                     igrid_nbor_amr = son(nbor(igrid_amr,igshift))
                  end if

                  if(igrid_nbor_amr==0) then
                     scan_flag=1
                     exit scan_flag_loop
                  else
                     icell_nbor_amr = igrid_nbor_amr + &
                           ncoarse+(jjj(idim,inbor,ind)-1)*ngridmax
                     if(f(icell_nbor_amr,3)<=0.0) then
                        scan_flag=1
                        exit scan_flag_loop
                     end if
                  end if
               end do
            end do scan_flag_loop
         else
            scan_flag=1
         end if
         ! Update flag2 with scan flag,
         ! BEWARE as lookup_mg backups are stored in flag2
         ! Safety init:
         if(flag2(icell_amr)>ngridmax .or. flag2(icell_amr)<0) flag2(icell_amr)=0
         ! Do NOT overwrite flag2 !
         flag2(icell_amr)=flag2(icell_amr)+ngridmax*scan_flag
      end do
   end do
end subroutine set_scan_flag_fine

subroutine get3cubefather_fine_mg(ind_cell_father,nbors_father_cells,&
     &                    nbors_father_grids,ncell,ilevel)
! Comments by kjhan: This subroutine seems to be thread safe.
  use amr_commons
  implicit none
  integer::ncell,ilevel
  integer, parameter::nvector_cg=32
  integer,dimension(1:nvector_cg)::ind_cell_father
  integer,dimension(1:nvector_cg,1:threetondim)::nbors_father_cells
  integer,dimension(1:nvector_cg,1:twotondim)::nbors_father_grids
  !------------------------------------------------------------------
  ! This subroutine determines the 3^ndim neighboring father cells 
  ! of the input father cell. According to the refinement rule, 
  ! they should be present anytime.
  !------------------------------------------------------------------
  integer::i,j,nxny,i1,j1,k1,ind,iok
  integer::i1min,i1max,j1min,j1max,k1min,k1max,ind_father
  integer,dimension(1:nvector_cg)::ix,iy,iz,iix,iiy,iiz
  integer,dimension(1:nvector_cg)::pos,ind_grid_father,ind_grid_ok
  integer,dimension(1:nvector_cg,1:threetondim)::nbors_father_ok
  integer,dimension(1:nvector_cg,1:twotondim)::nbors_grids_ok
  logical::oups

  nxny=nx*ny

  if(ilevel==1)then  ! Easy...

     oups=.false.
     do i=1,ncell
        if(ind_cell_father(i)>ncoarse)oups=.true.
     end do
     if(oups)then
        write(*,*)'get3cubefather_fine_mg'
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
        & call get3cubepos_fine_mg(ind_grid_ok,ind,nbors_father_ok,nbors_grids_ok,iok)

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

end subroutine get3cubefather_fine_mg

subroutine get3cubepos_fine_mg(ind_grid,ind,nbors_father_cells,nbors_father_grids,ng)
! Comments by kjhan: This subroutine seems to be thread safe.
  use amr_commons
  implicit none
  integer::ng,ind
  integer, parameter::nvector_cg=32
  integer,dimension(1:nvector_cg)::ind_grid
  integer,dimension(1:nvector_cg,1:threetondim)::nbors_father_cells
  integer,dimension(1:nvector_cg,1:twotondim)::nbors_father_grids
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
  integer,dimension(1:nvector_cg)::ind_grid1,ind_grid2,ind_grid3
  integer,dimension(1:nvector_cg,1:twotondim)::nbors_grids

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
  ! @ind =3
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

end subroutine get3cubepos_fine_mg



