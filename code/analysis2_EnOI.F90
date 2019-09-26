subroutine analysis2_EnOI(A,psi,D, R, S, ndim, nrens, nrobs, verbose)
! Computes the analysed vector psi for the static ensemble A 
   use m_multa
   implicit none
   integer, intent(in) :: ndim             ! dimension of model state
   integer, intent(in) :: nrens            ! number of ensemble members
   integer, intent(in) :: nrobs            ! number of observations
   
   real, intent(in) :: A(ndim,nrens)   ! ensemble matrix
   real, intent(inout) :: psi(ndim)       ! state vector (forecast -> vector)
   real, intent(in)    :: D(nrobs)  ! matrix  holding observation innovations
   real, intent(in)    :: S(nrobs,nrens)  ! matrix holding HA' 
   real, intent(inout) :: R(nrobs,nrobs)  ! Error covariance matrix for observations
   logical, intent(in) :: verbose


   real, allocatable, dimension(:,:) :: X1,X2,U,Reps,V,I_N
   real, allocatable, dimension(:) :: X4,X5
!   real X3(nrobs,1)
   real, allocatable, dimension(:,:) :: X3
   real, allocatable, dimension(:)   :: sig,work

   real sigsum,sigsum1,oneobs(1,1)
   integer ierr,nrsigma,i,j,lwork,m
   integer iblkmax, iens
   integer, parameter :: target= 3*22+1  ! Case dependent: test one value
   character(len=2) tag2



!   if(verbose) then
!            do i=1,nrobs
!            print *, i,R(i,i)*nrens
!           enddo 
!      endif
   if (nrobs > 1) then
!      R=float(nrens)*R+matmul(S,transpose(S))
      call dgemm('n','t', nrobs, nrobs, nrens, 1.0, S, nrobs, S, nrobs, float(nrens), R, nrobs)
!EnKF call dgemm('n','t', nrobs, nrobs, nrens, 1.0, S, nrobs, S, nrobs, float(nrens), R, nrobs)



      allocate (U(nrobs,nrobs) , stat=ierr )
      allocate (V(nrobs,nrobs) , stat=ierr )
      allocate (sig(nrobs) , stat=ierr )
      lwork=2*max(3*nrobs+nrobs,5*nrobs)
      allocate(work(lwork), stat=ierr)
      sig=0.0
      !R=U*sigma*V
      !note that V=transpose(V) ....... Nasty
!$OMP CRITICAL
      call dgesvd('A', 'A', nrobs, nrobs, R, nrobs, sig, U, nrobs, V, nrobs, work, lwork, ierr)
!$OMP END CRITICAL
      deallocate(work)
      if (ierr /= 0) then
         print *,'ierr from call dgesvd= ',ierr
         stop
      endif

      if(verbose) then
         open(10,file='sigma.dat')
            do i=1,nrobs
               write(10,'(i5,g12.3)')i,sig(i)
            enddo
         close(10)
      endif

      sigsum=sum( sig(1:nrobs) )
      sigsum1=0.0
   ! Significant eigenvalues. digma=sigma^-1
      nrsigma=0
      do i=1,nrobs                 ! singular values are in descending order
         if (sigsum1/sigsum < 0.999) then !using all the singular value until we
                                          !0.999 %
            nrsigma=nrsigma+1
            sigsum1=sigsum1+sig(i)
            sig(i) = 1.0/sig(i)
         else
            sig(i:nrobs)=0.0
            exit
         endif
      enddo

      if (verbose) then
         write(*,'(a,i5,g13.5)') ' dominant sing. values and share ',nrsigma,sigsum1/sigsum
         write(*,'(8g12.3)')1./sig(1:nrsigma), sig(nrsigma+1:nrobs)
      endif

      allocate (X1(nrobs,nrobs))
      do i=1,nrobs
      do j=1,nrobs
         X1(i,j)=sig(i)*U(j,i)
      enddo
      enddo
      deallocate(sig)

      allocate (X2(nrobs,1))
!     X2=X1*D : (nrobs,nrobs)*(nrobs,1)
      call dgemm('n','n',nrobs,1,nrobs,1.0,X1,nrobs,D ,nrobs,0.0,X2,nrobs)
      deallocate(X1) 

!     X3=transpose(transpose(V))*X2  : (nrobs,nnrobs)*(nrobs,1)
       allocate (X3(nrobs,1), stat=ierr)
       call dgemm('t','n',nrobs,1    ,nrobs,1.0,V ,nrobs,X2,nrobs,0.0,X3,nrobs)
!ENKF  call dgemm('t','n',nrobs,nrens,nrobs,1.0,V ,nrobs,X2,nrobs,0.0,X3,nrobs)
      deallocate(V)
      deallocate(X2)

   else
!      print *,'BEWARE ; ONLY 1 OBS FOUND !!!!!!!'
       allocate (X3(nrobs,1), stat=ierr)
      oneobs=matmul(S,transpose(S))+R*float(nrens)
      print *,'oneobs: ',oneobs(1,1)
      X3(:,1)=D/oneobs(1,1)
   endif
!
   if (2_8*ndim*nrobs < 1_8*nrens*(nrobs+ndim)) then
!    Code for few observations ( m<nN/(2n-N) )
!    Representer option 
      if (verbose) print * ,'analysis: Representer approach is used'
      allocate (Reps(ndim,nrobs))

!    Reps=matmul(A,transpose(S))  :  (ndim,nrens)*(nrobs,nrens)^T = (ndim,nrobs)
      call dgemm('n','t',ndim,nrobs,nrens,1.0,A,ndim,S,nrobs,0.0,Reps,ndim)
!
! Remove mean of Representers
    allocate (I_N(nrobs,nrobs))
    I_N=-1./float(nrobs)
    do i=1, nrobs
       I_N(i,i)=I_N(i,i)+1.
    end do
!
!   Reps = matmul(Reps,I_N)
!!!!!    call dgemm('n','n',ndim,nrobs,nrobs,1.0,Reps,ndim,I_N,nrobs,0.0,Reps,ndim)
    deallocate(I_N)
!
!    psi=psi+matmul(Reps,X3)
! EnKF      call dgemm('n','n',ndim,nrens,nrobs,1.0,Reps,ndim,X3,nrobs,1.0,A,ndim)
      allocate (X3(nrobs,1), stat=ierr)
      call dgemv('n',ndim,nrobs,1.0,Reps,ndim,X3,1,1.0,psi,1) ! single step operation
      deallocate(Reps)

      tag2(1:2)='X3'
      open(10,file='X.uf',form='unformatted')
         write(10)tag2,nrens,nrobs,X3,S
      close(10)

   else
      allocate(X4(nrens))
!      X4=matmul(transpose(S),X3)    : (nrens,nrobs)*(nrobs,1)
       call dgemm('t','n',nrens,1    ,nrobs,1.0,S,nrobs,X3,nrobs,0.0,X4,nrens)
!EnKF  call dgemm('t','n',nrens,nrens,nrobs,1.0,S,nrobs,X3,nrobs,0.0,X4,nrens)
     allocate (I_N(nrens,nrens))
      I_N=-1.0/float(nrens)
      do i=1,nrens
         I_N(i,i)=I_N(i,i)+1.0
      enddo

      allocate(X5(nrens))
      !in order to have A' we center X4: A'*X4=A*I_N*X4=A*X5
      X5=matmul(I_N,X4)
      deallocate(I_N)
      deallocate(X4)
      deallocate(X3)
      
!psi=psi+A*X5    :  (ndim,nrens)*(nrens)=ndim
!      call dgemm('n','n',ndim,1,nrens,1.0,A,ndim,X5,nrens,1.0,psi,ndim)
      call dgemv('n',ndim,nrens,1.0,A,ndim,X5,1,1.0,psi,1) ! single step operation

   deallocate(X5)
   if(verbose) print *,'Analysis, TEM(1) =', psi(target)
!
   endif 

end subroutine analysis2_EnOI
