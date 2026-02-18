module poisson_commons 
  use amr_commons
  use poisson_parameters


  real(dp),allocatable,dimension(:)  ::phi,phi_old       ! Potential
  real(dp),allocatable,dimension(:)  ::rho               ! Density
  real(dp),allocatable,dimension(:)  ::rho_star          ! Star density
  real(dp),allocatable,dimension(:,:)::f                 ! 3-force

  real(dp),allocatable,dimension(:)  ::rho_top   ! Density at last CIC level                                 

  ! Multigrid lookup table for amr -> mg index mapping
  integer, allocatable, dimension(:) :: lookup_mg   ! Lookup table

  ! Communicator arrays for multigrid levels
  type(communicator), allocatable, dimension(:,:) :: active_mg
  type(communicator), allocatable, dimension(:,:) :: emission_mg

  ! Minimum MG level
  integer :: levelmin_mg

  ! Multigrid safety switch
  logical, allocatable, dimension(:) :: safe_mode

  ! Pre-computed neighbor grids for fine-level MG solver
  integer, allocatable, dimension(:,:) :: nbor_grid_fine

  ! Pre-computed neighbor grids/cpus for coarse-level MG solver
  type mg_nbor_cache_t
     integer, allocatable, dimension(:,:) :: grid  ! (0:twondim, 1:ngrid)
     integer, allocatable, dimension(:,:) :: cpu   ! (1:twondim, 1:ngrid)
  end type mg_nbor_cache_t
  type(mg_nbor_cache_t), allocatable, dimension(:) :: nbor_mg_cache

  ! Coarse MG timing accumulators
  real(dp) :: mg_c_gs=0d0, mg_c_res=0d0, mg_c_restrict=0d0
  real(dp) :: mg_c_recurse=0d0, mg_c_interp=0d0, mg_c_total=0d0

  ! Multipole coefficients
  real(dp),dimension(1:ndim+1)::multipole

end module poisson_commons
