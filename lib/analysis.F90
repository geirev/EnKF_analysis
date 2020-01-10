subroutine analysis(A, R, E, S, D, innov, ndim, nrens, nrobs, verbose, truncation, mode, &
                    lrandrot,lupdate_randrot,lsymsqrt,inflate, infmult, ne)
! Computes the analysed ensemble for A using the EnKF or square root schemes.

   use mod_anafunc
   use m_multa
   use m_ensmean
   use m_ensvar
   implicit none
   
   integer, intent(in) :: ndim             ! dimension of model state
   integer, intent(in) :: nrens            ! number of ensemble members
   integer, intent(in) :: nrobs            ! number of observations
   integer, intent(in) :: ne               ! factor of increase ensemble size in E

   
   real, intent(inout) :: A(ndim,nrens)    ! ensemble matrix
   real, intent(in)    :: R(nrobs,nrobs)   ! matrix holding R (only used if mode=?1 or ?2)
   real, intent(in)    :: D(nrobs,nrens)   ! matrix holding perturbed measurments innovation d+E-HA 
   real, intent(in)    :: E(nrobs,nrens*ne)! matrix holding perturbations (only used if mode=?3)
   real, intent(in)    :: S(nrobs,nrens)   ! matrix holding HA` 
   real, intent(in)    :: innov(nrobs)     ! vector holding d-H*mean(A)

   logical, intent(in) :: verbose          ! Printing some diagnostic output

   real, intent(in)    :: truncation       ! The ratio of variance retained in pseudo inversion (0.99)

   integer, intent(in) :: mode             ! first integer means (EnKF=1, SQRT=2)
                                           ! Second integer is pseudo inversion
                                           !  1=eigen value pseudo inversion of SS'+(N-1)R
                                           !  2=SVD subspace pseudo inversion of SS'+(N-1)R
                                           !  3=SVD subspace pseudo inversion of SS'+EE'

   logical, intent(in) :: lrandrot         ! True if additional random rotation is used
   logical, intent(in) :: lupdate_randrot   ! Normally true; false for all but first grid point
                                           ! updates when using local analysis since all grid
                                           ! points need to use the same rotation.

   logical, intent(in) :: lsymsqrt         ! true if symmetrical sqare root of Sakov is used (should be used)

   integer, intent(in) :: inflate          ! Inflation(0=off, 1=multiplicative, 2=adaptive according to Evensen 2009)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   real X5(nrens,nrens)
   integer i,j,nrmin,iblkmax
   logical lreps
   real, intent(in)    :: infmult
   real inffac
   real ave(ndim)

   real, allocatable :: eig(:)
   real, allocatable :: Z(:,:)
   real, allocatable :: X2(:,:)
   real, allocatable :: X3(:,:)
   real, allocatable :: Reps(:,:)

   lreps = .FALSE.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Pseudo inversion of C=SS' +(N-1)*R
   if (nrobs == 1) then
      nrmin=1
      allocate(Z(1,1))
      allocate(eig(1))
      eig(1)=dot_product(S(1,:),S(1,:))+real(nrens-1)*R(1,1)
      eig(1)=1.0/eig(1)
      Z(1,1)=1.0

   else
      select case (mode)
      case(10)
         call exact_diag_inversion(S,D,X5,nrens,nrobs)

      case(11,21)
         nrmin=nrobs
!        Evaluate R= S*S` + (nrens-1)*R
         call dgemm('n','t',nrobs,nrobs,nrens, &
                       1.0, S, nrobs, &
                            S, nrobs, &
             real(nrens-1), R, nrobs)

!        Compute eigenvalue decomposition of R -> Z*eig*Z` 
         allocate(Z(nrobs,nrobs))
         allocate(eig(nrobs))
         call eigC(R,nrobs,Z,eig)
         call eigsign(eig,nrobs,truncation)

      case(12,22)
         nrmin=min(nrobs,nrens)
         allocate(Z(nrobs,nrmin))
         allocate(eig(nrmin))
         call lowrankCinv(S,R,nrobs,nrens,nrmin,Z,eig,truncation)

      case(13,23)
         nrmin=min(nrobs,nrens)
         allocate(Z(nrobs,nrmin))
         allocate(eig(nrmin))
         call lowrankE(S,E,nrobs,nrens,nrmin,Z,eig,truncation,ne)

      case default
         print *,'error analysis: Unknown mode: ',mode
         stop
      end select
   endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Generation of X5 (or representers in EnKF case with few measurements)
!   print *,'      analysis: Generation of X5:'
   select case (mode)
   case(10)

   case(11,12,13)
      allocate(X3(nrobs,nrens))
      if (nrobs > 1) then
         call genX3(nrens,nrobs,nrmin,eig,Z,D,X3)
      else
         X3=D*eig(1)
      endif

      if (2_8*ndim*nrobs < 1_8*nrens*(nrobs+ndim) .and. inflate/=2) then
!        Code for few observations ( m<nN/(2n-N) )
         if (verbose) print '(a)' ,'   analysis: Representer approach is used'
         lreps=.true.
         allocate (Reps(ndim,nrobs))
         call dgemm('n','t',ndim,nrobs,nrens,1.0,A,ndim,S,nrobs,0.0,Reps,ndim)
      else
         if (verbose) print '(a)' ,'   analysis: X5 approach is used'
         call dgemm('t','n',nrens,nrens,nrobs,1.0,S,nrobs,X3,nrobs,0.0,X5,nrens)
         do i=1,nrens
            X5(i,i)=X5(i,i)+1.0
         enddo
      endif

   case(21,22,23)
! Mean part of X5
      call meanX5(nrens,nrobs,nrmin,S,Z,eig,innov,X5)

! Generating X2
      allocate(X2(nrmin,nrens))
      call genX2(nrens,nrobs,nrmin,S,Z,eig,X2)

! Generating X5 matrix
      call X5sqrt(X2,nrobs,nrens,nrmin,X5,lrandrot,lupdate_randrot,mode,lsymsqrt)

   case default
      print *,'error analysis: Unknown flag for mode: ',mode
      stop
   end select


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Final ensemble update
   if (verbose) print '(a)' ,'   analysis: final update'
   if (lreps) then
      !     A=A+matmul(Reps,X3)
      call dgemm('n','n',ndim,nrens,nrobs,1.0,Reps,ndim,X3,nrobs,1.0,A,ndim)
   else
      iblkmax=min(ndim,200)
      call multa(A, X5, ndim, nrens, iblkmax )
   endif

   if (verbose) print '(a)' ,'   analysis: final update done'

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Inflation
   if (inflate==1) then
      inffac=infmult
   elseif (inflate==2) then
      call inflationfactor(X5,nrens,inffac)  ! Adaptive inflation factor
      inffac=1.0+(inffac-1.0)*infmult        ! Adjustment of addaptive
   endif

   if (inflate > 0) then
      print '(a,f10.4)','   analysis: inflation update with inflation factor= ',inffac
      call ensmean(A,ave,ndim,nrens)
      do j=1,nrens
      do i=1,ndim
         A(i,j)=ave(i) + (A(i,j)-ave(i))*inffac
      enddo
      enddo
   endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   if (allocated(X2))    deallocate(X2)
   if (allocated(X3))    deallocate(X3)
   if (allocated(eig))   deallocate(eig)
   if (allocated(Z))     deallocate(Z)
   if (allocated(Reps))  deallocate(Reps)

end subroutine
