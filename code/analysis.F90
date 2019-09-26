subroutine analysis(A, R, E, S, D, innov, ndim, nrens, nrobs, verbose, truncation,mode,update_randrot)
! Computes the analysed ensemble for A using the EnKF or square root schemes.

   use mod_anafunc
   use m_multa
   implicit none
   integer, intent(in) :: ndim             ! dimension of model state
   integer, intent(in) :: nrens            ! number of ensemble members
   integer, intent(in) :: nrobs            ! number of observations
   
   real, intent(inout) :: A(ndim,nrens)    ! ensemble matrix
   real, intent(in)    :: R(nrobs,nrobs)   ! matrix holding R (only used if mode=?1 or ?2)
   real, intent(in)    :: D(nrobs,nrens)   ! matrix holding perturbed measurments
   real, intent(in)    :: E(nrobs,nrens)   ! matrix holding perturbations (only used if mode=?3)
   real, intent(in)    :: S(nrobs,nrens)   ! matrix holding HA` 
   real, intent(in)    :: innov(nrobs)     ! vector holding d-H*mean(A)

   logical, intent(in) :: verbose          ! Printing some diagnostic output

   real, intent(in)    :: truncation       ! The ratio of variaince retained in pseudo inversion (0.99)

   integer, intent(in) :: mode             ! first integer means (EnKF=1, SQRT=2)
                                           ! Second integer is pseudo inversion
                                           !  1=eigen value pseudo inversion of SS'+(N-1)R
                                           !  2=SVD subspace pseudo inversion of SS'+(N-1)R
                                           !  3=SVD subspace pseudo inversion of SS'+EE'

   logical, intent(in) :: update_randrot   ! Normally true; false for all but first grid point
                                           ! updates when using local analysis since all grid
                                           ! points need to use the same rotation.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   real X5(nrens,nrens)
   integer i,nrmin,iblkmax
   logical lreps


   real, allocatable :: eig(:)
   real, allocatable :: W(:,:)
   real, allocatable :: X2(:,:)
   real, allocatable :: X3(:,:)
   real, allocatable :: Reps(:,:)



   lreps=.FALSE.
   if (verbose) print * ,'analysis: verbose is on'
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Pseudo inversion of C=SS' +(N-1)*R
   print *,'      analysis: Inversion of C:'
   if (nrobs == 1) then
      nrmin=1
      allocate(W(1,1))
      allocate(eig(1))
      eig(1)=dot_product(S(1,:),S(1,:))+real(nrens-1)*R(1,1)
      eig(1)=1.0/eig(1)
      W(1,1)=1.0

   else
      select case (mode)
      case(11,21)
         nrmin=nrobs
!        Evaluate R= S*S` + (nrens-1)*R
         call dgemm('n','t',nrobs,nrobs,nrens, &
                       1.0, S, nrobs, &
                            S, nrobs, &
             real(nrens-1), R, nrobs)

!        Compute eigenvalue decomposition of R -> W*eig*W` 
         allocate(W(nrobs,nrobs))
         allocate(eig(nrobs))
         call eigC(R,nrobs,W,eig)
         call eigsign(eig,nrobs,truncation)

      case(12,22)
         nrmin=min(nrobs,nrens)
         allocate(W(nrobs,nrmin))
         allocate(eig(nrmin))
         call lowrankCinv(S,R,nrobs,nrens,nrmin,W,eig,truncation)

      case(13,23)
         nrmin=min(nrobs,nrens)
         allocate(W(nrobs,nrmin))
         allocate(eig(nrmin))
         call lowrankE(S,E,nrobs,nrens,nrmin,W,eig,truncation)

      case default
         print *,'analysis: Unknown mode: ',mode
         stop
      end select
   endif








!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Generation of X5 (or representers in EnKF case with few measurements)
   print *,'      analysis: Generation of X5:'
   select case (mode)
   case(11,12,13)
      allocate(X3(nrobs,nrens))
      if (nrobs > 1) then
         call genX3(nrens,nrobs,nrmin,eig,W,D,X3)
      else
         X3=D*eig(1)
      endif

      if (2_8*ndim*nrobs < 1_8*nrens*(nrobs+ndim)) then
!        Code for few observations ( m<nN/(2n-N) )
         if (verbose) print * ,'analysis: Representer approach is used'
         lreps=.true.
         allocate (Reps(ndim,nrobs))
!        Reps=matmul(A,transpose(S))
         call dgemm('n','t',ndim,nrobs,nrens,1.0,A,ndim,S,nrobs,0.0,Reps,ndim)
      else
         if (verbose) print * ,'analysis: X5 approach is used'
!        X5=matmul(transpose(S),X3)
         call dgemm('t','n',nrens,nrens,nrobs,1.0,S,nrobs,X3,nrobs,0.0,X5,nrens)
         do i=1,nrens
            X5(i,i)=X5(i,i)+1.0
         enddo
      endif

   case(21,22,23)
! Mean part of X5
      call meanX5(nrens,nrobs,nrmin,S,W,eig,innov,X5)

! Generating X2
      allocate(X2(nrmin,nrens))
      call genX2(nrens,nrobs,nrmin,S,W,eig,X2)

! Generating X5 matrix
      call X5sqrt(X2,nrobs,nrens,nrmin,X5,update_randrot,mode)

   case default
      print *,'analysis: Unknown flag for mode: ',mode
      stop
   end select


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Generation of inflation
!   call inflationTEST(X5,nrens)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Final ensemble update
   print *,'      analysis: Final ensemble update:'
   if (lreps) then
!     A=A+matmul(Reps,X3)
      call dgemm('n','n',ndim,nrens,nrobs,1.0,Reps,ndim,X3,nrobs,1.0,A,ndim)
      call dumpX3(X3,S,nrobs,nrens)
   else
      iblkmax=min(ndim,200)
      call multa(A, X5, ndim, nrens, iblkmax )
      call dumpX5(X5,nrens)
   endif



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   if (allocated(X2))    deallocate(X2)
   if (allocated(X3))    deallocate(X3)
   if (allocated(eig))   deallocate(eig)
   if (allocated(W))     deallocate(W)
   if (allocated(Reps))  deallocate(Reps)
end subroutine
