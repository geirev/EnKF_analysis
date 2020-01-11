program main


! Test program for EnKF analysis
   use mod_dimensions
   use m_pseudo1D
   use m_fixsample1D
   use m_set_random_seed2
   use m_enkf
   use m_measurements
   use m_ensmean
   use m_ensvar
   use m_obspert
   use m_tecsol
   use m_dumpensemble
   use m_residuals
   implicit none
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Analysis scheme to use
   integer mode_analysis                           ! 10 Stochastic EnKF, Exact inversion with diagonal R
                                                   ! 11 Stochastic EnKF, Eigen value decomposition of (SS'+(N-1)R) 
                                                   ! 12 Stochastic EnKF, Eigen value decomposition of (SS'+(N-1)Re)  Re=EE'
                                                   ! 13 Stochastic EnKF, Subspace inversion of        (SS'+(N-1)EE)
                                                   ! 21 SQRT EnKF, Subspace inversion of        (SS'+R)
                                                   ! 22 SQRT EnKF, Subspace inversion of        (SS'+(N-1)Re)        Re=EE'
                                                   ! 23 SQRT EnKF, Subspace inversion of        (SS'+EE')
! Main tests of consistency between schemes: 
! nrobs < nrens, truncation=1.00, covmodel=diagonal, lsymsqrt=true, lrandrot=true, inflation=0, localization=0
!   ==> 10, 11, 21 gives exactly same solution for the mean (all uses exact diagonal R)
!   ==> 12, 13, 22, 23 gives exactly same solution for the mean (all uses ensemble R)
!   ==> 10, 11, gives exactly same variance (different solvers with exact R=I)
!   ==> 12, 13 gives exactly same variance (subspace solvers with ensemble R)
!   ==> 23, 23 gives exactly same variance (subspace solvers with ensemble R)


! Model for measurements and errors
   integer, parameter :: nrobs=50                  ! Number of measurements
   real               :: obsvar=0.25               ! Measurement variance
   character(len=8  ) :: covmodel='gaussian'       ! diagonal or gaussian
   real               :: rd=20.0                   ! Horizontal correlation of observation errors in Gaussian case
   integer, parameter :: ne=10                     ! scaling size of E used in the analysis scheme R=EE'

! Model ensemble 
   integer, parameter :: nrens=100                 ! ensemble size
   real               :: const=4.0                 ! mean of analytical solution
   real               :: rh=40.0                   ! Horizontal correlation of model fields
   real               :: dx=1.0                    ! horizontal grid spacing
   real               :: inivar=1.0                ! Variance of initial ensemble

! Options
   logical            :: lsymsqrt=.true.           ! Always use Sakovs symmetrical square root rather than one-sided
   logical            :: lrandrot=.false.           ! Introduce a mean-preserving random rotation in SQRT schemes
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
   real E0(nrobs,nrens*ne) ! A common matrix of measurment perts
   real E(nrobs,nrens*ne)  ! measurement pert
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
   print '(a)','main: sampling random truth'
   call pseudo1D(ana,nx,1,rh,dx,nx)
   ana=ana+const
   ic=1 ! Truth
   ave(:,ic)=ana(:)
   var(:,ic)=0.0
   print '(a)','main: truth ok'
   print *

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! First guess solution is a random perturbation from N(0,1,rh) added to the analytical truth
   print '(a)','main: generating random first guess'
   call pseudo1D(fg,nx,1,rh,dx,nx)
   fg=(fg + ana-const)/sqrt(2.0) +  const 
   ic=2 ! first guess
   ave(:,ic)=fg(:)
   var(:,ic)=inivar
   print '(a)','main: first guess ok'
   print *

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Generating the measurements uniformly throughout the domain
   print '(a)','main: generating measurements by first guess'
   dobs=real(nx)/real(nrobs)
   do m=1,nrobs
      obspos(m)= nint(real(m-1)*dobs + 0.5*dobs)
   enddo
   call measurements(ana,nx,obs,obspos,nrobs,obsvar,covmodel,rd,dx)
   print '(a)','main: measurements ok'
   print *

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Construct observation perturbations E
   print '(a)','main: sampling measurement pert into E0'
   call obspert(E0,nrens*ne,nrobs,.true.,dx,rd,covmodel,obspos)
! Introduce correct variances
   E0(:,:)=sqrt(obsvar)*E0(:,:)
   print '(a)','main: measurement pert ok'
   print *

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Initial ensemble is a random perturbation from N(0,inivar,rh) added to the analytical truth
   print '(a)','main: Generating initial ensemble'
   call pseudo1D(mem0,nx,nrens,rh,dx,nx)
   call fixsample1D(mem0,nx,nrens)
   do j=1,nrens
      mem0(:,j)=fg(:) + sqrt(inivar)*mem0(:,j)
   enddo
   ic=3 ! Prior (should be equal to fg)
   call ensmean(mem0,ave(1,ic),nx,nrens)
   call ensvar(mem0,ave(1,ic),var(1,ic),nx,nrens)
   print '(a)','main: ensemble ok'
   print *

! Stochastic EnKF 
   do i=0,3
      mode_analysis=10+i
      ic=ic+1
      print '(a)','++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
      print '(a,i2,a,i2,tr2,a)','main: calling enkf with mode_analysis=',mode_analysis,' ic=',ic,cc(ic)
      mem=mem0
      E=E0
      call enkf(mem,nx,nrens,obs,obsvar,obspos,nrobs,mode_analysis,&
               &truncation,covmodel,dx,rh,rd,lrandrot,lsymsqrt,&
               &inflate,infmult,local,robs,obstreshold,E,ne)
      call ensmean(mem,ave(1,ic),nx,nrens)
      call ensvar(mem,ave(1,ic),var(1,ic),nx,nrens)
      call dumpensemble(mem,ave,var,nrens,nx,ic,cc,nc)


   enddo

! SQRT EnKF
   do i=1,3
      mode_analysis=20+i
      ic=ic+1
      print '(a)','++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
      print '(a,i2,a,i2,tr2,a)','main: calling enkf with mode_analysis=',mode_analysis,' ic=',ic,cc(ic)
      mem=mem0
      E=E0
      call enkf(mem,nx,nrens,obs,obsvar,obspos,nrobs,mode_analysis,&
               &truncation,covmodel,dx,rh,rd,lrandrot,lsymsqrt,&
               &inflate,infmult,local,robs,obstreshold,E,ne)
      call ensmean(mem,ave(1,ic),nx,nrens)
      call ensvar(mem,ave(1,ic),var(1,ic),nx,nrens)
      call dumpensemble(mem,ave,var,nrens,nx,ic,cc,nc)
   enddo
   print '(a)','++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Dumping the solutions in tecplot format
   call tecsol(ave,var,obs,obspos,nx,nrobs,nc)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Dumping a recidual matrix for the different cases
   call residuals(ave,var,cc,nx,nc)

end program main

