program main
! Test program for EnKF analysis
   use mod_dimensions
   use m_pseudo1D
   use m_set_random_seed2
   use m_enkf
   use m_measurements
   use m_ensmean
   use m_ensvar
   use m_tecsol
   use m_residuals
   implicit none
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Analysis scheme to use
   integer mode_analysis                           ! 10 Stochastic EnKF, Exact inversion with diagonal R
                                                   ! 11 Stochastic EnKF, Eigen value decomposition of (SS'+(N-1)R) 
                                                   ! 12 Stochastic EnKF, Eigen value decomposition of (SS'+(N-1)R) 
                                                   ! 13 Stochastic EnKF, Subspace inversion of        (SS'+(N-1)R)
                                                   ! 21 SQRT EnKF, Subspace inversion of        (SS'+EE')
                                                   ! 22 SQRT EnKF, Subspace inversion of        (SS'+(N-1)R)
                                                   ! 23 SQRT EnKF, Subspace inversion of        (SS'+EE')
! Main tests of consistency between schemes: 
! Rexact=.true., nrobs < nrens, truncation=1.00, covmodel=diagonal, lsymsqrt=true, lrandrot=true, inflation=0, localization=0
!   ==> 10, 11, 12, 21, 22 gives exactly same solution for the mean (13 and 23 uses ensemble R)
!   ==> 13, 23 gives exactly same solution for the mean
!   ==> 10, 11, 12 gives exactly same variance
!   ==> 21, 22 gives exactly same variance

! Rexact=.false., nrobs < nrens, truncation=1.00, covmodel=diagonal, lsymsqrt=true, lrandrot=true, inflation=0, localization=0
!   ==> 11, 12, 13, 21, 22, 23 gives exactly same solution for the mean (10 still uses Rexact implicitly)
!   ==> 11, 12, 13 gives exactly same variance
!   ==> 21, 22, 23 gives exactly same variance


! Model for measurements and errors
   integer, parameter :: nrobs=50                  ! Number of measurements
   logical            :: Rexact=.true.             ! Use exact(true) or lowrank(false) R matrix (allows for testing the impact of lowrank R)
   real               :: obsvar=0.25               ! Measurement variance
   character(len=8  ) :: covmodel='diagonal'       ! diagonal or gaussian
   real               :: rd=40.0                   ! Horizontal correlation of observation errors in Gaussian case

! Model ensemble 
   integer, parameter :: nrens=100                 ! ensemble size
   real               :: const=4.0                 ! mean of analytical solution
   real               :: rh=40.0                   ! Horizontal correlation of model fields
   real               :: dx=1.0                    ! horizontal grid spacing
   real               :: inivar=1.0                ! Variance of initial ensemble

! Options
   logical            :: lsymsqrt=.true.           ! Always use Sakovs symmetrical square root rather than one-sided
   logical            :: lrandrot=.true.           ! Introduce a mean-preserving random rotation in SQRT schemes
   real               :: truncation=1.00           ! SVD truncation in inversion

! Inflation
   integer            :: inflate=0                 ! Inflation(0=off, 1=multiplicative, 2=adaptive according to Evensen 2009)
   real               :: infmult=1.0               ! Factor for multiplicative inflation (also used to scale adaptive inflation if /=1.0)

! Localization
   integer            :: local=0                   ! Localization (0=off, 1=distance_based (robs), 2=adaptive(obstreshold))
   real               :: robs=200                  ! Influence radius for localization
   real               :: obstreshold=0.0           ! Correlation threshold when adaptive is used

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Cases run
   integer, parameter :: nc=10
   character(len=12) :: cc(1:nc)  =(/ 'truth       ','prior       ','prior       ',&
                                     &'10          ','11          ','12          ','13          ',& 
                                     &'21          ','22          ','23          '/)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   real mem(nx,nrens)      ! Work ensemble
   real mem0(nx,nrens)     ! Initial ensemble
   real E0(nrobs,nrens)    ! A common matrix of measurment perts
   real ana(nx)            ! analytical solution
   real fg(nx)             ! first guess
   real ave(nx,nc)         ! ensemble average
   real var(nx,nc)         ! ensemble variance
   integer obspos(nrobs)   ! Position of data point
   real    obs(nrobs)      ! Observations
   integer i,j,m,ic        ! Counters etc.
   real dobs

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   call set_random_seed2   ! New random seed every time
   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! The true solution is a perturbation around the value "const" where the 
! perturbation is a smooth pseudo random field drawn from  N(0,1,rh).
   call pseudo1D(ana,nx,1,rh,dx,nx)
   ana=ana+const
   ave(:,1)=ana(:)
   var(:,1)=0.0
   print '(a)','main: ana ok'

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! First guess solution is a random perturbation from N(0,1,rh) added to the analytical truth
   call pseudo1D(fg,nx,1,rh,dx,nx)
   fg=(fg + ana-const)/sqrt(2.0) +  const !+ana    
   ave(:,2)=fg(:)
   var(:,2)=inivar
   print '(a)','main: fg ok'

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Generating the measurements uniformly throughout the domain
   dobs=real(nx)/real(nrobs)
   do m=1,nrobs
      obspos(m)= nint(real(m-1)*dobs + 0.5*dobs)
   enddo
   call measurements(ana,nx,obs,obspos,nrobs,obsvar,covmodel,rd,dx)
   print '(a)','main: measurements ok'

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Initial ensemble is a random perturbation from N(0,inivar,rh) added to the analytical truth
   call pseudo1D(mem0,nx,nrens,rh,dx,nx)
   do j=1,nrens
      mem0(:,j)=fg(:) + sqrt(inivar)*mem0(:,j)
   enddo
   call ensmean(mem0,ave(1,3),nx,nrens)
   call ensvar(mem0,ave(1,3),var(1,3),nx,nrens)
   print '(a)','main: ensemble ok'

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Computation of analyses
   E0(1,1)=0.0 ! used in EnKF to only resample E0 the first time enkf is called

! Stochastic EnKF 
   do i=0,3
      mode_analysis=10+i
      ic=4+i
      print '(a)','++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
      print '(a,i2,a,i2)','main: calling enkf with mode_analysis=',mode_analysis,' ic=',ic
      mem=mem0
      call enkf(mem,nx,nrens,obs,obsvar,obspos,nrobs,1,1,mode_analysis,&
               &truncation,covmodel,dx,rh,Rexact,rd,lrandrot,lsymsqrt,&
               &inflate,infmult,local,robs,obstreshold,E0)
      call ensmean(mem,ave(1,ic),nx,nrens)
      call ensvar(mem,ave(1,ic),var(1,ic),nx,nrens)
   enddo

! SQRT EnKF
   do i=1,3
      mode_analysis=20+i
      ic=7+i
      print '(a)','++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
      print '(a,i2,a,i2)','main: calling enkf with mode_analysis=',mode_analysis,' ic=',ic
      mem=mem0
      call enkf(mem,nx,nrens,obs,obsvar,obspos,nrobs,1,1,mode_analysis,&
               &truncation,covmodel,dx,rh,Rexact,rd,lrandrot,lsymsqrt,&
               &inflate,infmult,local,robs,obstreshold,E0)
      call ensmean(mem,ave(1,ic),nx,nrens)
      call ensvar(mem,ave(1,ic),var(1,ic),nx,nrens)
   enddo
   print '(a)','++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Dumping the solutions in tecplot format
   call tecsol(ave,var,obs,obspos,nx,nrobs,nc)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Dumping a recidual matrix for the different cases
   call residuals(ave,var,cc,nx,nc)

end program main

