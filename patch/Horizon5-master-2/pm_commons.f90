! Patch changes:
! - added mp0 variable for particles
module pm_commons
  use amr_parameters
  use pm_parameters
  use random
  ! Sink particle related arrays
  real(dp),allocatable,dimension(:)::msink,r2sink,v2sink,c2sink,oksink_new,oksink_all,tsink
  real(dp),allocatable,dimension(:)::msink_new,msink_all,r2k,v2sink_new,c2sink_new,tsink_new,tsink_all
  real(dp),allocatable,dimension(:)::v2sink_all,c2sink_all
  real(dp),allocatable,dimension(:)::dMBHoverdt,dMEdoverdt,wdens,wvol,wc2
  real(dp),allocatable,dimension(:)::wdens_new,wvol_new,wc2_new,total_volume
  real(dp),allocatable,dimension(:,:)::wmom,wmom_new
  real(dp),allocatable,dimension(:,:)::vsink,vsink_new,vsink_all
  real(dp),allocatable,dimension(:,:)::xsink,xsink_new,xsink_all
  real(dp),allocatable,dimension(:,:)::weighted_density,weighted_volume,weighted_c2
  real(dp),allocatable,dimension(:,:)::jsink,jsink_new,jsink_all
  real(dp),allocatable,dimension(:)::dMBH_coarse,dMEd_coarse,dMsmbh,dMBH_coarse_new
  real(dp),allocatable,dimension(:)::dMEd_coarse_new,dMsmbh_new,dMBH_coarse_all,dMEd_coarse_all,dMsmbh_all
  real(dp),allocatable,dimension(:)::Esave,Esave_new,Esave_all
  real(dp),allocatable,dimension(:,:,:)::weighted_momentum
  real(dp),allocatable,dimension(:,:,:)::sink_stat,sink_stat_all
  real(dp),allocatable,dimension(:)::c_avgptr,v_avgptr,d_avgptr
  real(dp),allocatable,dimension(:)::spinmag,spinmag_new,spinmag_all
  real(dp),allocatable,dimension(:,:)::bhspin,bhspin_new,bhspin_all
  real(dp),allocatable,dimension(:)::eps_sink
  integer ,allocatable,dimension(:)::idsink,idsink_new,idsink_all
  integer::nindsink=0

  ! Particles related arrays
  real(dp),allocatable,dimension(:,:)::xp       ! Positions
  real(dp),allocatable,dimension(:,:)::vp       ! Velocities
  real(dp),allocatable,dimension(:)  ::mp       ! Masses
#ifdef OUTPUT_PARTICLE_POTENTIAL
  real(dp),allocatable,dimension(:)  ::ptcl_phi ! Potential of particle added by AP for output purposes 
#endif
  real(dp),allocatable,dimension(:)  ::tp       ! Birth epoch
  real(dp),allocatable,dimension(:,:)::weightp  ! weight of cloud parts for sink accretion only
  real(dp),allocatable,dimension(:)  ::zp       ! Birth metallicity
  integer ,allocatable,dimension(:)  ::nextp    ! Next particle in list
  integer ,allocatable,dimension(:)  ::prevp    ! Previous particle in list
  integer ,allocatable,dimension(:)  ::levelp   ! Current level of particle
  integer(i8b),allocatable,dimension(:)::idp    ! Identity of particle

  ! Tree related arrays
  integer ,allocatable,dimension(:)  ::headp    ! Head particle in grid
  integer ,allocatable,dimension(:)  ::tailp    ! Tail particle in grid
  integer ,allocatable,dimension(:)  ::numbp    ! Number of particles in grid
  ! Global particle linked lists
  integer::headp_free,tailp_free,numbp_free=0,numbp_free_tot=0
  ! Local and current seed for random number generator
  integer,dimension(IRandNumSize) :: localseed=-1

  !for chemo components
  real(dp),allocatable,dimension(:)  ::tpp, mp0, indtab
!  real(dp),allocatable,dimension(:,:)::cep
  character(len=3),allocatable,dimension(:) ::elem_list

  type yield_table
     integer::na
     integer::nz
     real,dimension(:)      ,pointer::astar
     real,dimension(:)      ,pointer::zstar
     real,dimension(:,:,:)  ,pointer::Eeject
     real,dimension(:,:)    ,pointer::Zeject
     real,dimension(:,:)    ,pointer::Meject
     real,dimension(:,:)    ,pointer::NSN1  
     real,dimension(:,:)    ,pointer::NSN2 
  end type yield_table
  type(yield_table)::yieldtab


  contains
  function cross(a,b)
    use amr_parameters, only:dp
    real(dp),dimension(1:3)::a,b
    real(dp),dimension(1:3)::cross
    !computes the cross product c= a x b
    cross(1)=a(2)*b(3)-a(3)*b(2)
    cross(2)=a(3)*b(1)-a(1)*b(3)
    cross(3)=a(1)*b(2)-a(2)*b(1)
  end function cross
  
end module pm_commons
