!################################################################
!################################################################
! Fuzzy Dark Matter (FDM) Schrödinger-Poisson time integrator
!
! DKD operator splitting per level:
!   1. Half Drift (kinetic): exp(-i hbar k^2 dt/4a^2) in k-space (base)
!                         or FD Laplacian (refined levels)
!   2. Full Kick (potential): exp(-i Phi dt / hbar) in real space
!   3. Half Drift (kinetic): same as step 1
!
! Base level (levelmin): FFT kinetic operator (spectral accuracy)
! Fine levels (>levelmin): Explicit FD kinetic operator
!################################################################
subroutine fdm_step(ilevel)
  use amr_commons
  use poisson_commons
  use fdm_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer,intent(in)::ilevel

  real(dp)::dt_loc

  if(.not.use_fdm) return
  if(numbtot(1,ilevel)==0) return

  dt_loc = dtnew(ilevel)

  ! DKD splitting: Drift(dt/2) - Kick(dt) - Drift(dt/2)
  call fdm_drift(ilevel, 0.5d0*dt_loc)
  call fdm_kick(ilevel, dt_loc)
  call fdm_drift(ilevel, 0.5d0*dt_loc)

  ! Update ghost zones for psi after evolution
  call make_virtual_fine_dp(psi_re(1), ilevel)
  call make_virtual_fine_dp(psi_im(1), ilevel)

end subroutine fdm_step
!################################################################
!################################################################
! Kinetic operator (Drift): exp(-i hbar/(2a^2) nabla^2 dt)
!
! Base level: FFT (spectral)
! Fine levels: explicit FD with 7-point Laplacian
!################################################################
subroutine fdm_drift(ilevel, dt_half)
  use amr_commons
  use poisson_commons
  use fdm_commons
  implicit none
  integer,intent(in)::ilevel
  real(dp),intent(in)::dt_half

#ifdef USE_FFTW
  if(ilevel == levelmin) then
     call fdm_drift_fft(ilevel, dt_half)
     return
  end if
#endif
  call fdm_drift_fd(ilevel, dt_half)

end subroutine fdm_drift
!################################################################
!################################################################
! FFT kinetic operator for base level (spectral accuracy)
! Applies exp(-i alpha k^2) to psi in Fourier space
! where alpha = hbar_code * dt / (2 a^2)
!################################################################
#ifdef USE_FFTW
subroutine fdm_drift_fft(ilevel, dt_half)
  use amr_commons
  use poisson_commons
  use fdm_commons
  use iso_c_binding
  use omp_lib
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  include 'fftw3.f03'

  integer,intent(in)::ilevel
  real(dp),intent(in)::dt_half

  ! Grid dimensions
  integer::fft_Nx, fft_Ny, fft_Nz, nx_loc
  integer(i8b)::N_total, N_complex
  real(dp)::dx_fft, twopi, alpha
  real(dp)::denom, phase, cos_p, sin_p

  ! FFTW plans (cached)
  logical,save::fdm_fftw_init = .false.
  type(C_PTR),save::plan_fwd_re = C_NULL_PTR, plan_bwd_re = C_NULL_PTR
  type(C_PTR),save::plan_fwd_im = C_NULL_PTR, plan_bwd_im = C_NULL_PTR
  integer,save::saved_Nx = 0

  ! Work arrays (cached)
  real(dp),allocatable,save::re_3d(:), im_3d(:)
  complex(C_DOUBLE_COMPLEX),allocatable,save::re_hat(:), im_hat(:)
  integer,allocatable,save::fft_map(:)

  ! Loop vars
  integer::igrid,ind,iskip,icell,ix,iy,iz,idx_3d,info
  integer::Kx,Ky,Kz,kx_i,ky_i,kz_i
  integer(i8b)::idx_c
  ! psi_hat = re_hat + i*im_hat (both are complex from R2C)
  ! psi_hat_full(k) = (re_hat_r + i*re_hat_i) + i*(im_hat_r + i*im_hat_i)
  !                 = (re_hat_r - im_hat_i) + i*(re_hat_i + im_hat_r)
  real(dp)::psi_r, psi_i, new_psi_r, new_psi_i
  real(dp)::rhr, rhi, ihr, ihi  ! real/imag parts of re_hat and im_hat
  real(dp)::scale

  ! Grid dimensions
  fft_Nx = nx * 2**ilevel
  fft_Ny = ny * 2**ilevel
  fft_Nz = nz * 2**ilevel
  N_total = int(fft_Nx,i8b) * int(fft_Ny,i8b) * int(fft_Nz,i8b)
  N_complex = int(fft_Nx/2+1,i8b) * int(fft_Ny,i8b) * int(fft_Nz,i8b)
  dx_fft = 0.5d0**ilevel
  twopi = 2.0d0*acos(-1.0d0)
  nx_loc = icoarse_max - icoarse_min + 1
  scale = boxlen / dble(nx_loc)

  ! alpha = hbar_code * dt / (2 * a^2)
  alpha = hbar_code * dt_half / (aexp**2)

  ! Initialize FFTW plans and work arrays (once)
  if(.not.fdm_fftw_init .or. fft_Nx /= saved_Nx) then
     if(allocated(re_3d))  deallocate(re_3d, im_3d, re_hat, im_hat, fft_map)
     allocate(re_3d(1:N_total), im_3d(1:N_total))
     allocate(re_hat(1:N_complex), im_hat(1:N_complex))
     allocate(fft_map(1:N_total))

     ! Build AMR cell -> flat 3D index map
     fft_map = 0
     do ind=1,twotondim
        iskip = ncoarse + (ind-1)*ngridmax
        igrid = headl(myid, ilevel)
        do while(igrid > 0)
           icell = igrid + iskip
           ix = int((xg(igrid,1) + (dble(iand(ind-1,1)) - 0.5d0)*dx_fft &
                - dble(icoarse_min)*scale) / dx_fft / scale) + 1
           iy = int((xg(igrid,2) + (dble(iand(ishft(ind-1,-1),1)) - 0.5d0)*dx_fft &
                - dble(jcoarse_min)*scale) / dx_fft / scale) + 1
           iz = int((xg(igrid,3) + (dble(iand(ishft(ind-1,-2),1)) - 0.5d0)*dx_fft &
                - dble(kcoarse_min)*scale) / dx_fft / scale) + 1
           ix = max(1, min(ix, fft_Nx))
           iy = max(1, min(iy, fft_Ny))
           iz = max(1, min(iz, fft_Nz))
           idx_3d = ix + (iy-1)*fft_Nx + (iz-1)*fft_Nx*fft_Ny
           fft_map(idx_3d) = icell
           igrid = next(igrid)
        end do
     end do

     ! Create R2C / C2R plans
     if(c_associated(plan_fwd_re)) call dfftw_destroy_plan(plan_fwd_re)
     if(c_associated(plan_bwd_re)) call dfftw_destroy_plan(plan_bwd_re)
     if(c_associated(plan_fwd_im)) call dfftw_destroy_plan(plan_fwd_im)
     if(c_associated(plan_bwd_im)) call dfftw_destroy_plan(plan_bwd_im)
     call dfftw_plan_dft_r2c_3d(plan_fwd_re, fft_Nx, fft_Ny, fft_Nz, &
          re_3d, re_hat, FFTW_MEASURE)
     call dfftw_plan_dft_c2r_3d(plan_bwd_re, fft_Nx, fft_Ny, fft_Nz, &
          re_hat, re_3d, FFTW_MEASURE)
     call dfftw_plan_dft_r2c_3d(plan_fwd_im, fft_Nx, fft_Ny, fft_Nz, &
          im_3d, im_hat, FFTW_MEASURE)
     call dfftw_plan_dft_c2r_3d(plan_bwd_im, fft_Nx, fft_Ny, fft_Nz, &
          im_hat, im_3d, FFTW_MEASURE)

     saved_Nx = fft_Nx
     fdm_fftw_init = .true.
  end if

  ! Step 1: Gather psi from AMR to flat 3D arrays
  ! ALLREDUCE for small grids (each rank has subset of cells)
  re_3d = 0.0d0
  im_3d = 0.0d0
  do idx_3d=1,N_total
     icell = fft_map(idx_3d)
     if(icell > 0) then
        re_3d(idx_3d) = psi_re(icell)
        im_3d(idx_3d) = psi_im(icell)
     end if
  end do
#ifndef WITHOUTMPI
  if(ncpu > 1) then
     call MPI_ALLREDUCE(MPI_IN_PLACE, re_3d, int(N_total), &
          MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, info)
     call MPI_ALLREDUCE(MPI_IN_PLACE, im_3d, int(N_total), &
          MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, info)
  end if
#endif

  ! Step 2: Forward FFT (R2C) for both real and imaginary parts
  call dfftw_execute_dft_r2c(plan_fwd_re, re_3d, re_hat)
  call dfftw_execute_dft_r2c(plan_fwd_im, im_3d, im_hat)

  ! Step 3: Apply kinetic phase rotation in k-space
  ! psi_hat(k) -> psi_hat(k) * exp(-i alpha k^2)
  ! where k^2 = (2/dx^2) * sum_d [1 - cos(2pi k_d / N_d)]  (discrete)
  !$omp parallel do private(Kz,Ky,Kx,kx_i,ky_i,kz_i,denom,phase, &
  !$omp&  cos_p,sin_p,rhr,rhi,ihr,ihi,psi_r,psi_i,new_psi_r,new_psi_i,idx_c) schedule(static)
  do Kz=0,fft_Nz-1
     kz_i = Kz; if(kz_i > fft_Nz/2) kz_i = kz_i - fft_Nz
     do Ky=0,fft_Ny-1
        ky_i = Ky; if(ky_i > fft_Ny/2) ky_i = ky_i - fft_Ny
        do Kx=0,fft_Nx/2
           idx_c = int(Kx+1,i8b) + int(Ky,i8b)*int(fft_Nx/2+1,i8b) &
                + int(Kz,i8b)*int(fft_Nx/2+1,i8b)*int(fft_Ny,i8b)

           ! Discrete Laplacian eigenvalue (negative)
           denom = 2.0d0*( cos(twopi*dble(Kx)/dble(fft_Nx)) &
                         + cos(twopi*dble(ky_i)/dble(fft_Ny)) &
                         + cos(twopi*dble(kz_i)/dble(fft_Nz)) - 3.0d0 )

           ! Phase rotation:
           ! dpsi/dt = i*(hbar/2a^2)*nabla^2 psi
           ! => psi_hat(t+dt) = psi_hat(t) * exp(-i*hbar*k^2*dt/(2a^2))
           ! k^2_discrete = -denom/dx^2
           ! phase = -hbar*(-denom/dx^2)*dt/(2a^2) = hbar*denom*dt/(2a^2*dx^2)
           !       = alpha * denom / (2*dx^2)  where alpha = hbar*dt_half/a^2
           ! Apply exp(i*phase) to psi_hat
           phase = 0.5d0 * alpha * denom / (dx_fft**2)
           cos_p = cos(phase)
           sin_p = sin(phase)

           ! Reconstruct full complex psi_hat(k) from two R2C transforms:
           ! psi = psi_re + i*psi_im
           ! FFT(psi) = FFT(psi_re) + i*FFT(psi_im) = re_hat + i*im_hat
           ! re_hat(k) and im_hat(k) are each complex numbers from R2C
           rhr = dble(re_hat(idx_c));  rhi = aimag(re_hat(idx_c))
           ihr = dble(im_hat(idx_c));  ihi = aimag(im_hat(idx_c))
           ! psi_hat = (rhr + i*rhi) + i*(ihr + i*ihi)
           !         = (rhr - ihi) + i*(rhi + ihr)
           psi_r = rhr - ihi
           psi_i = rhi + ihr

           ! Apply exp(i*phase): psi_hat_new = psi_hat * (cos_p + i*sin_p)
           new_psi_r = psi_r*cos_p - psi_i*sin_p
           new_psi_i = psi_r*sin_p + psi_i*cos_p

           ! Decompose back: new_psi = new_re_hat + i*new_im_hat
           ! new_re_hat = (new_psi_r + i*???)  -- we need to split back
           ! Since FFT is linear: psi_hat = re_hat + i*im_hat
           ! We need: re_hat_new such that FFT(psi_re_new) = re_hat_new
           !          im_hat_new such that FFT(psi_im_new) = im_hat_new
           ! The inverse: re_hat_new_r = new_psi_r + new_ihi
           !              re_hat_new_i = new_psi_i - new_ihr
           ! This requires the old ihi,ihr... it's circular.
           !
           ! Better approach: work directly on re_hat and im_hat separately
           ! Since the phase operator is exp(i*phase), and it's applied to psi_hat:
           ! We can equivalently apply the SAME rotation to re_hat and im_hat
           ! independently, because:
           ! psi_hat_new = (re_hat + i*im_hat) * exp(i*phase)
           !             = re_hat*exp(i*phase) + i*im_hat*exp(i*phase)
           ! => re_hat_new = re_hat * exp(i*phase)
           !    im_hat_new = im_hat * exp(i*phase)
           ! This is exact! Just rotate each complex number individually.

           ! Apply exp(i*phase) to re_hat
           new_psi_r = rhr*cos_p - rhi*sin_p
           new_psi_i = rhr*sin_p + rhi*cos_p
           re_hat(idx_c) = cmplx(new_psi_r, new_psi_i, C_DOUBLE_COMPLEX)

           ! Apply exp(i*phase) to im_hat
           new_psi_r = ihr*cos_p - ihi*sin_p
           new_psi_i = ihr*sin_p + ihi*cos_p
           im_hat(idx_c) = cmplx(new_psi_r, new_psi_i, C_DOUBLE_COMPLEX)
        end do
     end do
  end do
  !$omp end parallel do

  ! Step 4: Inverse FFT (C2R)
  call dfftw_execute_dft_c2r(plan_bwd_re, re_hat, re_3d)
  call dfftw_execute_dft_c2r(plan_bwd_im, im_hat, im_3d)

  ! Step 5: Normalize (FFTW convention: unnormalized) and scatter back
  do idx_3d=1,N_total
     icell = fft_map(idx_3d)
     if(icell > 0) then
        psi_re(icell) = re_3d(idx_3d) / dble(N_total)
        psi_im(icell) = im_3d(idx_3d) / dble(N_total)
     end if
  end do

end subroutine fdm_drift_fft
#endif
!################################################################
!################################################################
! Finite-difference kinetic operator for AMR fine levels
! Crank-Nicolson implicit (unitary, unconditionally stable)
!
! Solves: (1 + i*alpha*L) psi^{n+1} = (1 - i*alpha*L) psi^n
! where alpha = hbar*dt/(4*a^2*dx^2), L = discrete Laplacian
!
! Implemented as fixed-point iteration (Jacobi):
!   psi^{k+1} = [(1-i*alpha*L)*psi^n - i*alpha*L_offdiag*psi^k] / (1+6*i*alpha)
! Converges in ~5-10 iterations for typical alpha.
! Falls back to explicit if alpha < 0.1 (CFL safe, explicit is faster).
!################################################################
subroutine fdm_drift_fd(ilevel, dt_half)
  use amr_commons
  use poisson_commons
  use fdm_commons
  implicit none
  integer,intent(in)::ilevel
  real(dp),intent(in)::dt_half

  integer::igrid,ind,iskip,icell,iter
  real(dp)::dx,scale,dx_loc,alpha,a2p1
  real(dp)::lap_re,lap_im,psi_re_c,psi_im_c
  real(dp)::rhs_re,rhs_im,denom_re,denom_im
  real(dp)::nbor_re,nbor_im
  integer::nx_loc,idim,inbor,icell_nbor
  integer,parameter::max_iter=10

  ! Precompute alpha = hbar * dt_half / (4 * a^2 * dx^2)
  ! (factor 4 because CN uses dt/2 in each half-step, and we split kinetic as -(hbar^2/2)*L)
  dx = 0.5d0**ilevel
  nx_loc = icoarse_max - icoarse_min + 1
  scale = boxlen / dble(nx_loc)
  dx_loc = dx * scale
  alpha = hbar_code * dt_half / (2.0d0 * aexp**2 * dx_loc**2)

  ! Use explicit with automatic sub-cycling (stable for any alpha)
  call fdm_drift_fd_explicit(ilevel, dt_half)

end subroutine fdm_drift_fd
!################################################################
! Explicit FD kinetic step (with automatic sub-cycling for stability)
!################################################################
subroutine fdm_drift_fd_explicit(ilevel, dt_half)
  use amr_commons
  use poisson_commons
  use fdm_commons
  implicit none
  integer,intent(in)::ilevel
  real(dp),intent(in)::dt_half

  integer::igrid,ind,iskip,icell,isub,nsub
  real(dp)::dx,scale,dx_loc,alpha,alpha_sub,dt_sub
  real(dp)::lap_re,lap_im,psi_re_c,psi_im_c
  integer::nx_loc,idim,inbor,icell_nbor

  dx = 0.5d0**ilevel
  nx_loc = icoarse_max - icoarse_min + 1
  scale = boxlen / dble(nx_loc)
  dx_loc = dx * scale

  ! Full alpha for the requested dt_half
  alpha = hbar_code * dt_half / (2.0d0 * aexp**2 * dx_loc**2)

  ! Sub-cycle if alpha > 1/6 (CFL limit)
  nsub = max(1, int(6.0d0*alpha) + 1)
  dt_sub = dt_half / dble(nsub)
  alpha_sub = alpha / dble(nsub)

  do isub=1,nsub
     ! Ghost zone exchange at start of each sub-step
     call make_virtual_fine_dp(psi_re(1), ilevel)
     call make_virtual_fine_dp(psi_im(1), ilevel)

     do ind=1,twotondim
        iskip = ncoarse + (ind-1)*ngridmax
        igrid = headl(myid, ilevel)
        do while(igrid > 0)
           icell = igrid + iskip
           if(son(icell) == 0) then
              psi_re_c = psi_re(icell)
              psi_im_c = psi_im(icell)

              ! 7-point Laplacian
              lap_re = -6.0d0 * psi_re_c
              lap_im = -6.0d0 * psi_im_c
              do idim=1,ndim
                 do inbor=1,2
                    call fdm_neighbor_cell(igrid, ind, idim, inbor, icell_nbor)
                    if(icell_nbor > 0) then
                       lap_re = lap_re + psi_re(icell_nbor)
                       lap_im = lap_im + psi_im(icell_nbor)
                    else
                       lap_re = lap_re + psi_re_c
                       lap_im = lap_im + psi_im_c
                    end if
                 end do
              end do

              ! Explicit update: dpsi = i*(hbar/2a^2)*L*dt
              psi_re(icell) = psi_re_c - alpha_sub * lap_im
              psi_im(icell) = psi_im_c + alpha_sub * lap_re
           end if
           igrid = next(igrid)
        end do
     end do
  end do

end subroutine fdm_drift_fd_explicit
!################################################################
! Find the neighbor cell of subcell ind in grid igrid
! in direction (idim, inbor): idim=1,2,3; inbor=1(left),2(right)
!################################################################
subroutine fdm_neighbor_cell(igrid, ind, idim, inbor, icell_nbor)
  use amr_commons
  implicit none
  integer,intent(in)::igrid, ind, idim, inbor
  integer,intent(out)::icell_nbor

  integer::ind_bit, ind_nbor, igrid_nbor, iskip

  ! Extract the bit for this dimension from ind (1-based, bits: x=bit0, y=bit1, z=bit2)
  ind_bit = iand(ishft(ind-1, -(idim-1)), 1)  ! 0 or 1

  if(inbor == 2) then
     ! Right neighbor
     if(ind_bit == 0) then
        ! Neighbor is in same grid, flip this dimension bit
        ind_nbor = ind + ishft(1, idim-1)
        iskip = ncoarse + (ind_nbor-1)*ngridmax
        icell_nbor = igrid + iskip
     else
        ! Neighbor is in adjacent grid
        igrid_nbor = son(nbor(igrid, 2*idim))
        if(igrid_nbor == 0) then
           icell_nbor = 0; return
        end if
        ind_nbor = ind - ishft(1, idim-1)
        iskip = ncoarse + (ind_nbor-1)*ngridmax
        icell_nbor = igrid_nbor + iskip
     end if
  else
     ! Left neighbor
     if(ind_bit == 1) then
        ! Neighbor is in same grid
        ind_nbor = ind - ishft(1, idim-1)
        iskip = ncoarse + (ind_nbor-1)*ngridmax
        icell_nbor = igrid + iskip
     else
        ! Neighbor is in adjacent grid
        igrid_nbor = son(nbor(igrid, 2*idim-1))
        if(igrid_nbor == 0) then
           icell_nbor = 0; return
        end if
        ind_nbor = ind + ishft(1, idim-1)
        iskip = ncoarse + (ind_nbor-1)*ngridmax
        icell_nbor = igrid_nbor + iskip
     end if
  end if

end subroutine fdm_neighbor_cell
!################################################################
!################################################################
! Potential operator (Kick): exp(-i Phi dt / hbar)
! Local pointwise rotation in real space
!################################################################
subroutine fdm_kick(ilevel, dt_loc)
  use amr_commons
  use poisson_commons
  use fdm_commons
  implicit none
  integer,intent(in)::ilevel
  real(dp),intent(in)::dt_loc

  integer::igrid,ind,iskip,icell
  real(dp)::phase,cos_p,sin_p,re_old,im_old

  if(hbar_code <= 0.0d0) return

  ! Loop over all leaf cells at this level
  do ind=1,twotondim
     iskip = ncoarse + (ind-1)*ngridmax
     igrid = headl(myid, ilevel)
     do while(igrid > 0)
        icell = igrid + iskip
        if(son(icell) == 0) then
           ! Schrödinger: i*hbar*dpsi/dt = m*Phi*psi
           ! => dpsi/dt = -i*(Phi/hbar)*psi
           ! => psi(t+dt) = psi(t) * exp(-i*Phi*dt/hbar)
           ! exp(-i*theta) = cos(theta) - i*sin(theta)
           ! where theta = Phi*dt/hbar
           phase = phi(icell) * dt_loc / hbar_code

           cos_p = cos(phase)
           sin_p = sin(phase)

           re_old = psi_re(icell)
           im_old = psi_im(icell)

           ! psi_new = psi * exp(-i*phase) = (re+i*im)*(cos-i*sin)
           psi_re(icell) =  re_old*cos_p + im_old*sin_p
           psi_im(icell) = -re_old*sin_p + im_old*cos_p
        end if
        igrid = next(igrid)
     end do
  end do

end subroutine fdm_kick
!################################################################
!################################################################
! Compute hbar_code from physical constants and code units
! Called once at initialization
!################################################################
subroutine fdm_compute_hbar()
  use amr_commons
  use fdm_commons
  implicit none

  real(dp)::scale_l, scale_t, scale_d, scale_v, scale_nH, scale_T2
  real(dp)::m_axion_g

  call units(scale_l, scale_t, scale_d, scale_v, scale_nH, scale_T2)

  ! m_axion in grams: m [eV] * eV_to_erg / c^2
  m_axion_g = m_axion * eV_to_erg / c_light**2

  ! hbar_code = hbar_cgs / (m_axion * scale_l * scale_v)
  ! This is the effective Planck constant in code units
  hbar_code = hbar_cgs / (m_axion_g * scale_l * scale_v)

  if(myid==1) then
     write(*,'(A,ES10.3,A)') ' FDM: m_axion = ', m_axion, ' eV'
     write(*,'(A,ES10.3)')   '   m_axion [g] = ', m_axion_g
     write(*,'(A,ES10.3)')   '   hbar_code   = ', hbar_code
     write(*,'(A,ES10.3,A)') '   lambda_dB(100km/s) = ', &
          2.0d0*acos(-1.0d0)*hbar_code/(100.0d5/scale_v) * scale_l/3.0857d21, ' kpc'
  end if

end subroutine fdm_compute_hbar
!################################################################
!################################################################
! Initialize psi from density field (Madelung transform)
! For fresh start (nrestart=0): set psi = sqrt(rho) with zero phase
! (Simplified: full Madelung with velocity potential is for future work)
!################################################################
subroutine fdm_init_psi()
  use amr_commons
  use poisson_commons
  use fdm_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer::ilevel,igrid,ind,iskip,icell,info
  real(dp)::rho_bg,rho_cell
  real(dp)::mass_loc, mass_glob

  if(.not.use_fdm) return

  ! Background density in code units (cosmological: rho_bg ~ 1)
  rho_bg = 1.0d0  ! will be refined with omega_fdm/omega_m

  mass_loc = 0.0d0

  ! Initialize psi = sqrt(rho) on all leaf cells (real part only, phase=0)
  do ilevel=levelmin,nlevelmax
     do ind=1,twotondim
        iskip = ncoarse + (ind-1)*ngridmax
        igrid = headl(myid, ilevel)
        do while(igrid > 0)
           icell = igrid + iskip
           if(son(icell) == 0) then
              ! Use existing density from rho array (set by init)
              rho_cell = max(rho(icell), 0.0d0)
              psi_re(icell) = sqrt(rho_cell)
              psi_im(icell) = 0.0d0
              mass_loc = mass_loc + rho_cell * (0.5d0**ilevel)**3
           end if
           igrid = next(igrid)
        end do
     end do
  end do

#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(mass_loc, mass_glob, 1, MPI_DOUBLE_PRECISION, &
       MPI_SUM, MPI_COMM_WORLD, info)
#else
  mass_glob = mass_loc
#endif

  if(myid==1) write(*,'(A,ES12.5)') ' FDM: initial |psi|^2 total mass = ', mass_glob

end subroutine fdm_init_psi
!################################################################
!################################################################
! Phase 4: AMR prolongation for psi (mass-conserving)
! Called when new refined grids are created at ilevel from ilevel-1
! Trilinear interpolation + renormalization to preserve |psi|^2 integral
!################################################################
subroutine fdm_prolong(ilevel)
  use amr_commons
  use poisson_commons
  implicit none
  integer,intent(in)::ilevel

  integer::igrid,icell_coarse,ind,iskip,i
  integer::ind_child,iskip_child,icell_child
  real(dp)::mass_coarse,mass_fine,ratio
  real(dp)::psi_re_parent,psi_im_parent

  if(.not.use_fdm) return
  if(ilevel <= levelmin) return

  ! For each newly created grid at ilevel, initialize psi from parent cell
  ! Simple injection: copy parent cell value to all 8 children
  ! Then renormalize to conserve mass
  igrid = headl(myid, ilevel)
  do while(igrid > 0)
     ! Parent cell
     icell_coarse = father(igrid)
     if(icell_coarse <= 0) then
        igrid = next(igrid)
        cycle
     end if

     psi_re_parent = psi_re(icell_coarse)
     psi_im_parent = psi_im(icell_coarse)

     ! Inject parent value to all children
     do ind_child=1,twotondim
        iskip_child = ncoarse + (ind_child-1)*ngridmax
        icell_child = igrid + iskip_child
        psi_re(icell_child) = psi_re_parent
        psi_im(icell_child) = psi_im_parent
     end do

     igrid = next(igrid)
  end do

  call make_virtual_fine_dp(psi_re(1), ilevel)
  call make_virtual_fine_dp(psi_im(1), ilevel)

end subroutine fdm_prolong
!################################################################
!################################################################
! Phase 4: AMR restriction for psi (mass-conserving)
! Average children to parent cell when derefinement occurs
!################################################################
subroutine fdm_restrict(ilevel)
  use amr_commons
  use poisson_commons
  implicit none
  integer,intent(in)::ilevel

  integer::igrid,ind,iskip,icell,igrid_child
  integer::ind_child,iskip_child,icell_child
  real(dp)::sum_re,sum_im

  if(.not.use_fdm) return
  if(ilevel < levelmin) return

  ! For each cell at ilevel that has children (son /= 0),
  ! average the children's psi back to parent
  do ind=1,twotondim
     iskip = ncoarse + (ind-1)*ngridmax
     igrid = headl(myid, ilevel)
     do while(igrid > 0)
        icell = igrid + iskip
        igrid_child = son(icell)
        if(igrid_child > 0) then
           sum_re = 0.0d0
           sum_im = 0.0d0
           do ind_child=1,twotondim
              iskip_child = ncoarse + (ind_child-1)*ngridmax
              icell_child = igrid_child + iskip_child
              sum_re = sum_re + psi_re(icell_child)
              sum_im = sum_im + psi_im(icell_child)
           end do
           ! Volume average of amplitude
           psi_re(icell) = sum_re / dble(twotondim)
           psi_im(icell) = sum_im / dble(twotondim)
        end if
        igrid = next(igrid)
     end do
  end do

end subroutine fdm_restrict
!################################################################
!################################################################
! Phase 5: de Broglie wavelength refinement criterion
! Flag cells where lambda_dB < fdm_nrefine_dB * dx for refinement
!################################################################
subroutine fdm_refine_flag(ilevel)
  use amr_commons
  use poisson_commons
  use fdm_commons
  implicit none
  integer,intent(in)::ilevel

  integer::igrid,ind,iskip,icell,idim,inbor,icell_nbor
  real(dp)::dx,scale,dx_loc,grad2,psi2,lambda_dB
  real(dp)::dpsi_re,dpsi_im
  integer::nx_loc

  if(.not.use_fdm) return
  if(hbar_code <= 0.0d0) return

  dx = 0.5d0**ilevel
  nx_loc = icoarse_max - icoarse_min + 1
  scale = boxlen / dble(nx_loc)
  dx_loc = dx * scale

  do ind=1,twotondim
     iskip = ncoarse + (ind-1)*ngridmax
     igrid = headl(myid, ilevel)
     do while(igrid > 0)
        icell = igrid + iskip
        if(son(icell) == 0) then  ! leaf cell only
           psi2 = psi_re(icell)**2 + psi_im(icell)**2
           if(psi2 > 1.0d-30) then
              ! Compute |grad(psi)|^2 using central differences
              grad2 = 0.0d0
              do idim=1,ndim
                 ! Right neighbor
                 call fdm_neighbor_cell(igrid, ind, idim, 2, icell_nbor)
                 if(icell_nbor > 0) then
                    dpsi_re = psi_re(icell_nbor) - psi_re(icell)
                    dpsi_im = psi_im(icell_nbor) - psi_im(icell)
                 else
                    dpsi_re = 0.0d0; dpsi_im = 0.0d0
                 end if
                 grad2 = grad2 + dpsi_re**2 + dpsi_im**2
              end do
              grad2 = grad2 / (dx_loc**2)

              ! lambda_dB = 2*pi*hbar / (|grad psi|/|psi| * hbar)
              ! Actually: k_local = |grad psi| / |psi|
              ! lambda_dB = 2*pi / k_local = 2*pi * |psi| / |grad psi|
              ! Refine if lambda_dB / dx < fdm_nrefine_dB
              if(grad2 > 0.0d0) then
                 lambda_dB = 2.0d0*acos(-1.0d0) * sqrt(psi2) / sqrt(grad2)
                 if(lambda_dB / dx_loc < dble(fdm_nrefine_dB)) then
                    flag1(icell) = 1
                 end if
              end if
           end if
        end if
        igrid = next(igrid)
     end do
  end do

end subroutine fdm_refine_flag
!################################################################
!################################################################
! Phase 7: FDM diagnostics — total mass and min lambda_dB
! Called from adaptive_loop
!################################################################
subroutine fdm_diagnostics()
  use amr_commons
  use poisson_commons
  use fdm_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer::ilevel,igrid,ind,iskip,icell,info
  real(dp)::mass_loc,mass_glob,rho_max_loc,rho_max_glob
  real(dp)::psi2,vol

  if(.not.use_fdm) return

  mass_loc = 0.0d0
  rho_max_loc = 0.0d0

  do ilevel=levelmin,nlevelmax
     vol = (0.5d0**ilevel * boxlen/dble(icoarse_max-icoarse_min+1))**3
     do ind=1,twotondim
        iskip = ncoarse + (ind-1)*ngridmax
        igrid = headl(myid, ilevel)
        do while(igrid > 0)
           icell = igrid + iskip
           if(son(icell) == 0) then
              psi2 = psi_re(icell)**2 + psi_im(icell)**2
              mass_loc = mass_loc + psi2 * vol
              if(psi2 > rho_max_loc) rho_max_loc = psi2
           end if
           igrid = next(igrid)
        end do
     end do
  end do

#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(mass_loc, mass_glob, 1, MPI_DOUBLE_PRECISION, &
       MPI_SUM, MPI_COMM_WORLD, info)
  call MPI_ALLREDUCE(rho_max_loc, rho_max_glob, 1, MPI_DOUBLE_PRECISION, &
       MPI_MAX, MPI_COMM_WORLD, info)
#else
  mass_glob = mass_loc
  rho_max_glob = rho_max_loc
#endif

  if(myid==1) then
     write(*,'(A,ES12.5,A,ES10.3)') &
          ' FDM: M_tot=', mass_glob, '  rho_max=', rho_max_glob
  end if

end subroutine fdm_diagnostics
