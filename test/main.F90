program main
   use mod_dimensions
   use m_sample1D
   use m_set_random_seed2
   use m_enkflocal
   use m_measurements
   use m_ensmean
   use m_ensvar
   implicit none

   integer mode_analysis                           ! 11 12 13 21 22 23
   real               :: rh=40.0                   ! Horizontal correlation
   real               :: const=4.0                 ! mean of analytical solution
   real               :: dx=1.0                    ! horizontal grid spacing
   integer            :: nrens=1000                 ! ensemble size
   logical            :: samp_fix=.true.
   real               :: inivar=1.0 
   real               :: obsvar=0.25
   integer            :: nrobs=200                   ! Number of measurement per assimilation time
   logical            :: mkobs=.true.              ! Create or read measurements
   logical            :: lrandrot=.true.
   logical            :: lupdate_randrot=.true.
   logical            :: lsakov=.true.             ! Always use the symmetrical square root rather than one-sided
   logical            :: linflate=.false.
   logical            :: ladapinf=.true.
   character(len=100) :: covmodel='gaussian'       ! diagonal or gaussian
   logical            :: Rexact=.false.             ! Use exact(true) or lowrank(false) R matrix
   real               :: truncation=0.99          ! SVD truncation in inversion


   logical            :: llocal=.false.             !Localization on or off
   logical            :: distance_based_loc=.true. ! Distance based
   real               :: robs=200                  ! Influence radius for localization

   logical            :: adaptive_based_loc=.false.! Adaptive localization
   real               :: obstreshold=0.0           ! Correlation threshold when adaptive is used


   real infmult




! other variables
   integer i,j,m,nn,ic,jc          ! Counters etc.


   real, allocatable :: mem(:,:)
   real, allocatable :: mem0(:,:)

   real ana(nx)  ! analytical solution
   real fg(nx)   ! first guess
   real ave(nx,0:10)  ! ensemble average
   real var(nx,0:10)  ! ensemble variance
   real resave(6,6)
   real resvar(6,6)


   integer, allocatable :: obspos(:)     ! Position of data point
   real,    allocatable :: obs(:)        ! Observations

   call set_random_seed2

   allocate (mem0(nx,nrens))   ! Initial ensemble
   allocate (mem(nx,nrens))    ! Work ensemble
   allocate(obs(nrobs))
   allocate(obspos(nrobs))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! The true solution is a perturbation around the value "const" where the 
! perturbation is a smooth pseudo random field drawn from  N(0,1,rh).
   call sample1D(ana,nx,1,1,1,dx,rh,.false.,.true.)
   ana=ana+const
   ave(:,1)=ana(:)
   var(:,1)=0.0
   print *,'main: ana ok'

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! First guess solution is a random perturbation from N(0,1,rh) added to the analytical truth
   call sample1D(fg,nx,1,1,1,dx,rh,.false.,.true.)
   fg=(fg + ana-const)/sqrt(2.0) +  const !+ana    
   ave(:,2)=fg(:)
   var(:,2)=inivar
   print *,'main: fg ok'

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Generating the measurements uniformly throughout the domain
   nn=nint(float(nx)/float(nrobs))
   obspos(1)=nint(real(nn)/2.0)
   do m=2,nrobs
      obspos(m)=min(obspos(m-1)+nn,nx)
   enddo

   open(10,file='obspos.dat')
      do m=1,nrobs
         write(10,'(i4,i4)')m,obspos(m)
      enddo
   close(10)
   print *,'main: obs ok'
   call measurements(ana,nx,obs,obspos,nrobs,obsvar,1,mkobs,0.0)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Initial ensemble is a random perturbation from N(0,inivar,rh) added to the analytical truth
   call sample1D(mem0,nx,nrens,1,1,dx,rh,samp_fix,.true.)
   do j=1,nrens
      mem0(:,j)=fg(:) + sqrt(inivar)*mem0(:,j)
   enddo
   call ensmean(mem0,ave(1,3),nx,nrens)
   call ensvar(mem0,ave(1,3),var(1,3),nx,nrens)
   print *,'main: ensemble ok'

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Computation of analyses

! stochastic EnKF
   do i=1,3
      mode_analysis=10+i
      ic=3+i
      print *,'mode_analysis=',mode_analysis,' ic=',ic

      mem=mem0
      call enkflocal(mem,nx,nrens,obs,obsvar,obspos,nrobs,samp_fix,1,1,mode_analysis,&
               &truncation,covmodel,dx,rh,Rexact,lrandrot,lupdate_randrot,lsakov,&
               &linflate,ladapinf,infmult,&
               &llocal,robs,distance_based_loc, adaptive_based_loc, obstreshold)

      call ensmean(mem,ave(1,ic),nx,nrens)
      call ensvar(mem,ave(1,ic),var(1,ic),nx,nrens)
   enddo

! SQRT EnKF
   do i=1,3
      mode_analysis=20+i
      ic=6+i
      print *,'mode_analysis=',mode_analysis,' ic=',ic

      mem=mem0
      call enkflocal(mem,nx,nrens,obs,obsvar,obspos,nrobs,samp_fix,1,1,mode_analysis,&
               &truncation,covmodel,dx,rh,Rexact,lrandrot,lupdate_randrot,lsakov,&
               &linflate,ladapinf,infmult,&
               &llocal,robs,distance_based_loc, adaptive_based_loc, obstreshold)

      call ensmean(mem,ave(1,ic),nx,nrens)
      call ensvar(mem,ave(1,ic),var(1,ic),nx,nrens)
   enddo

   open(10,file='solutions.dat')
   write(10,*)'TITLE = "Solutions"'
   write(10,*)'VARIABLES = "iens" "Truth" "First guess" "Prior" "11" "12" "13" "21" "22" "23"'
   write(10,'(a,i5,a,i5,a)')' ZONE T="Average"  F=BLOCK, I=',nx,', J=1, K=1'
   write(10,'(20I5)')(i         ,i=1,nx)
   do ic=1,9
      write(10,'(20g13.5)')(ave(i,ic) ,i=1,nx)
   enddo

   write(10,'(a,i5,a,i5,a)')' ZONE T="Std Dev"  F=BLOCK, I=',nx,', J=1, K=1'
   write(10,'(20I5)')(i         ,i=1,nx)
   do ic=1,9
      write(10,'(20g13.5)')(var(i,ic) ,i=1,nx)
   enddo

   write(10,'(a,i5,a,i5,a)')' ZONE T="observations"  F=BLOCK, I=',nrobs,', J=1, K=1'
   write(10,'(20I5)')(obspos(i) ,i=1,nrobs)
   do ic=1,9
      write(10,'(20g13.5)')(obs(i) ,i=1,nrobs)
   enddo

   close(10)

   resave=0.0
   resvar=0.0
   var=sqrt(var)
   do jc=4,9
   do ic=4,9
      resave(ic-3,jc-3)=dot_product(ave(1:nx,ic)-ave(1:nx,jc),ave(1:nx,ic)-ave(1:nx,jc))/real(nx)
      resvar(ic-3,jc-3)=dot_product(var(1:nx,ic)-var(1:nx,jc),var(1:nx,ic)-var(1:nx,jc))/real(nx)
   enddo
   enddo

   write(*,'(a)')'res_ave:'
   write(*,'(a5,6i12)')'     ',(jc,jc=1,6)
   do jc=1,6
      write(*,'(i5,6f12.6)')jc,resave(1:jc,jc)
   enddo

   write(*,'(a)')'res_var:'
   do jc=1,6
      write(*,'(i5,6f12.6)')jc,resvar(1:jc,jc)
   enddo


end program main

