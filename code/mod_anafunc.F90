module mod_anafunc
contains

subroutine lowrankE(S,E,nrobs,nrens,nrmin,W,eig,truncation)
   implicit none
   integer, intent(in)  :: nrobs
   integer, intent(in)  :: nrens
   integer, intent(in)  :: nrmin
   real,    intent(in)  :: S(nrobs,nrens)
   real,    intent(in)  :: E(nrobs,nrens)
   real,    intent(out) :: W(nrobs,nrmin)
   real,    intent(out) :: eig(nrmin)
   real,    intent(in)  :: truncation

   real U0(nrobs,nrmin),sig0(nrmin)
   real X0(nrmin,nrens)
   integer i,j

   real U1(nrmin,nrmin),VT1(1,1)
   real, allocatable :: work(:)
   integer lwork
   integer ierr

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Compute SVD of S=HA`  ->  U0, sig0
   call  svdS(S,nrobs,nrens,nrmin,U0,sig0,truncation)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Compute X0=sig0^{*T} U0^T E 

! X0= U0^T R
   call dgemm('t','n',nrmin,nrens,nrobs, 1.0,U0,nrobs, E,nrobs, 0.0,X0,nrmin)


   do j=1,nrens
   do i=1,nrmin
      X0(i,j)=sig0(i)*X0(i,j)
   enddo
   enddo


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Compute singular value decomposition  of X0(nrmin,nrens)
   lwork=2*max(3*nrens+nrobs,5*nrens)
   allocate(work(lwork))
   eig=0.0

   call dgesvd('S', 'N', nrmin, nrens, X0, nrmin, eig, U1, nrmin, VT1, 1, work, lwork, ierr)
   deallocate(work)
   if (ierr /= 0) then
      print *,'mod_anafunc (lowrankE): ierr from call dgesvd 1= ',ierr; stop
   endif

   do i=1,nrmin
      eig(i)=1.0/(1.0+eig(i)**2)
   enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! W = U0 * sig0^{-1} * U1
   do j=1,nrmin
   do i=1,nrmin
      U1(i,j)=sig0(i)*U1(i,j)
   enddo
   enddo

   call dgemm('n','n',nrobs,nrmin,nrmin, 1.0,U0,nrobs, U1,nrmin, 0.0,W,nrobs)


end subroutine


subroutine eigC(R,nrobs,Z,eig)
! Compute eigenvalue decomposition of R -> Z*eig*Z` 
   integer, intent(in) :: nrobs
   real, intent(in)    :: R(nrobs,nrobs)
   real, intent(out) :: Z(nrobs,nrobs)
   real, intent(out)   :: eig(nrobs)

#ifdef IBM
   real, allocatable :: ap(:)
   integer k
#endif
   real  RR(nrobs,nrobs)

   real fwork(8*nrobs)  
   integer iwork(5*nrobs)
   integer ifail(nrobs)
   real abstol,ddum
   integer idum,neig,ierr
   real, external :: DLAMCH

   idum=1

#ifdef IBM
! Upper packed storage as in ESSL manual
   allocate (ap(nrobs*(nrobs+1)/2) )
   k=0
   do j=1,nrobs
   do i=1,j
       k=k+1
       ap(k)=R(i,j)
   enddo
   enddo
   call dspev(21,ap,eig,Z,nrobs,nrobs,fwork,2*nrobs)
   deallocate(ap)
#else
   abstol=2.0*DLAMCH('S')
   RR=R
   call dsyevx('V', 'A', 'U', nrobs, RR, nrobs, ddum, ddum, idum, idum, abstol, &
            neig, eig, Z, nrobs, fwork, 8*nrobs, iwork, ifail, ierr )
#endif


end subroutine



subroutine eigsign(eig,nrobs,truncation)
! Returns the inverse of the truncated eigenvalue spectrum
implicit none
integer, intent(in)    :: nrobs
real,    intent(inout) :: eig(nrobs)
real,    intent(in)    :: truncation

integer i,nrsigma
real sigsum,sigsum1
logical ex

   inquire(file='eigenvalues.dat',exist=ex)
   if (ex) then
      open(10,file='eigenvalues.dat',position='append')
         write(10,'(a,i5,a)')' ZONE  F=POINT, I=',nrobs,' J=1 K=1'
         do i=1,nrobs
            write(10,'(i3,g13.5)')i,eig(nrobs-i+1)
         enddo
      close(10)
   else
      open(10,file='eigenvalues.dat')
         write(10,*)'TITLE = "Eigenvalues of C"'
         write(10,*)'VARIABLES = "obs" "eigenvalues"'
         write(10,'(a,i5,a)')' ZONE  F=POINT, I=',nrobs,' J=1 K=1'
         do i=1,nrobs
            write(10,'(i3,g13.5)')i,eig(nrobs-i+1)
         enddo
      close(10)
   endif

! Significant eigenvalues
   sigsum=sum( eig(1:nrobs) )
   sigsum1=0.0
   nrsigma=0
   do i=nrobs,1,-1
!      print '(a,i5,g13.5)','Eigen values: ',i,eig(i)
      if (sigsum1/sigsum < truncation) then
         nrsigma=nrsigma+1
         sigsum1=sigsum1+eig(i)
         eig(i) = 1.0/eig(i)
      else
         eig(1:i)=0.0
         exit
      endif
   enddo
   write(*,'(2(a,i5))')      '       analysis: Number of dominant eigenvalues: ',nrsigma,' of ',nrobs
   write(*,'(2(a,g13.4),a)') '       analysis: Share (and truncation)        : ',sigsum1/sigsum,' (',truncation,')'


end subroutine



subroutine genX2(nrens,nrobs,idim,S,W,eig,X2)
! Generate X2= (I+eig)^{-0.5} * W^T * S
   implicit none
   integer, intent(in) :: nrens
   integer, intent(in) :: nrobs
   integer, intent(in) :: idim ! idim=nrobs for A4 and nrmin for A5
   real, intent(in)    :: W(idim,nrens)
   real, intent(in)    :: S(nrobs,nrens)
   real, intent(in)    :: eig(idim)
   real, intent(out)   :: X2(idim,nrens)
   integer i,j

   call dgemm('t','n',idim,nrens,nrobs,1.0,W,nrobs, S,nrobs, 0.0,X2,idim)

   do j=1,nrens
   do i=1,idim
      X2(i,j)=sqrt(eig(i))*X2(i,j)
   enddo
   enddo

end subroutine



subroutine genX3(nrens,nrobs,nrmin,eig,W,D,X3)
   implicit none
   integer, intent(in) :: nrens
   integer, intent(in) :: nrobs
   integer, intent(in) :: nrmin
   real,    intent(in) :: eig(nrmin)
   real,    intent(in) :: W(nrobs,nrmin)
   real,    intent(in) :: D(nrobs,nrens)
   real,    intent(out) :: X3(nrobs,nrmin)

   real X1(nrmin,nrobs)
   real X2(nrmin,nrens)
   integer i,j

   do i=1,nrmin
   do j=1,nrobs
      X1(i,j)=eig(i)*W(j,i)
   enddo
   enddo

!     X2=matmul(X1,D)
      call dgemm('n','n',nrmin,nrens,nrobs,1.0,X1,nrmin,D ,nrobs,0.0,X2,nrmin)

!     X3=matmul(W,X2)
      call dgemm('n','n',nrobs,nrens,nrmin,1.0,W ,nrobs,X2,nrmin,0.0,X3,nrobs)

end subroutine



subroutine meanX5(nrens,nrobs,nrmin,S,W,eig,innov,X5)
   implicit none
   integer, intent(in) :: nrens
   integer, intent(in) :: nrobs
   integer, intent(in) :: nrmin
   real, intent(in)    :: W(nrmin,nrmin)
   real, intent(in)    :: S(nrobs,nrens)
   real, intent(in)    :: eig(nrmin)
   real, intent(in)    :: innov(nrobs)
   real, intent(out)   :: X5(nrens,nrens)

   real y1(nrmin)
   real y2(nrmin)
   real y3(nrobs)
   real y4(nrens) 
   integer i

   if (nrobs==1) then
      y1(1)=W(1,1)*innov(1)
      y2(1)=eig(1)*y1(1)
      y3(1)=W(1,1)*y2(1)
      y4(:)=y3(1)*S(1,:)
   else
      call dgemv('t',nrobs,nrmin,1.0,W,nrobs,innov,1,0.0,y1 ,1)
      y2=eig*y1  
      call dgemv('n',nrobs,nrmin,1.0,W ,nrobs,y2,1,0.0,y3 ,1)
      call dgemv('t',nrobs,nrens,1.0,S ,nrobs,y3,1,0.0,y4 ,1)
   endif

   do i=1,nrens
      X5(:,i)=y4(:)
   enddo

! X5=enN + (I - enN) X5  = enN + X5
   X5=1.0/real(nrens) + X5

end subroutine



subroutine X5sqrt(X2,nrobs,nrens,nrmin,X5,update_randrot,mode)
   use m_randrot
   use m_mean_preserving_rotation
   implicit none
   integer, intent(in) :: nrobs
   integer, intent(in) :: nrens
   integer, intent(inout) :: nrmin ! note that nrmin=nrobs in a4
   real, intent(in)    :: X2(nrmin,nrens)
   real, intent(inout) :: X5(nrens,nrens)
   logical, intent(in) :: update_randrot
   integer, intent(in) :: mode

   real X3(nrens,nrens)
   real X33(nrens,nrens)
   real X4(nrens,nrens)
   real IenN(nrens,nrens)
   real, save, allocatable :: ROT(:,:)


   real U(nrmin,1),sig(nrmin),VT(nrens,nrens)
   real, allocatable, dimension(:)   :: work,isigma
   integer i,j,lwork,ierr

   print *,'      analysis (X5sqrt): update_randrot= ',update_randrot
   if (update_randrot) then
      if (allocated(ROT)) deallocate(ROT)
      allocate(ROT(nrens,nrens))
!!      ROT=0.0
!!      do i=1,nrens
!!         ROT(i,i)=1.0
!!      enddo
!      call randrot(ROT,nrens)  ! old non mean preserving random rotation
      call mean_preserving_rotation(ROT,nrens)
   endif

! SVD of X2
   lwork=2*max(3*nrens+nrens,5*nrens); allocate(work(lwork))
   sig=0.0
   call dgesvd('N', 'A', nrmin, nrens, X2, nrmin, sig, U, nrmin, VT, nrens, work, lwork, ierr)
   deallocate(work)
   if (ierr /= 0) then
      print *,'X5sqrt: ierr from call dgesvd = ',ierr
      stop
   endif


   if (mode == 21) nrmin=min(nrens,nrobs)
   allocate(isigma(nrmin))
   isigma=1.0
   do i=1,nrmin
      if ( sig(i) > 1.0 ) print *,'X5sqrt: WARNING (m_X5sqrt): sig > 1',i,sig(i)
      isigma(i)=sqrt( max(1.0-sig(i)**2,0.0) )
   enddo

   do j=1,nrens
      X3(:,j)=VT(j,:)
   enddo


   do j=1,nrmin
      X3(:,j)=X3(:,j)*isigma(j)
   enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Multiply  X3* V' = (V*sqrt(I-sigma*sigma) * V' to ensure symmetric sqrt and 
! mean preserving rotation.   Sakov paper eq 13
   call dgemm('n','n',nrens,nrens,nrens,1.0,X3,nrens,VT,nrens,0.0,X33,nrens)
! Check mean preservation X33*1_N = a* 1_N (Sakov paper eq 15)
!   do i=1,nrens
!      print *,'sum(X33)= ',i,sum(X33(i,:))
!   enddo
   
!!X33=X3

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!   print '(a)','X5sqrt: sig: '
!   print '(5g11.3)',sig(1:min(nrmin,nrens))

   call dgemm('n','n',nrens,nrens,nrens,1.0,X33,nrens,ROT,nrens,0.0,X4,nrens)

   IenN=-1.0/real(nrens)
   do i=1,nrens
      IenN(i,i)=  IenN(i,i) + 1.0
   enddo

   call dgemm('n','n',nrens,nrens,nrens,1.0,IenN,nrens,X4,nrens,1.0,X5,nrens)

   deallocate(isigma)

end subroutine



subroutine dumpX3(X3,S,nrobs,nrens)
   implicit none
   integer, intent(in) :: nrens
   integer, intent(in) :: nrobs
   real,    intent(in) :: X3(nrens,nrens)
   real,    intent(in) :: S(nrobs,nrens)
   character(len=2) :: tag2

   tag2(1:2)='X3'
   open(10,file='X5.uf',form='unformatted')
      write(10)tag2,nrens,nrobs,X3,S
   close(10)

end subroutine



subroutine dumpX5(X5,nrens)
   implicit none
   integer, intent(in) :: nrens
   real,    intent(in) :: X5(nrens,nrens)
   integer j
   character(len=2) :: tag2

   tag2(1:2)='X5'
   open(10,file='X5.uf',form='unformatted')
      write(10)tag2,nrens,X5
   close(10)

   open(10,file='X5col.dat')
      do j=1,nrens
         write(10,'(i5,f10.4)')j,sum(X5(:,j))
      enddo
   close(10)

   open(10,file='X5row.dat')
      do j=1,nrens
         write(10,'(i5,f10.4)')j,sum(X5(j,:))/real(nrens)
       enddo
   close(10)
end subroutine



subroutine lowrankCinv(S,R,nrobs,nrens,nrmin,W,eig,truncation)
   implicit none
   integer, intent(in)  :: nrobs
   integer, intent(in)  :: nrens
   integer, intent(in)  :: nrmin
   real,    intent(in)  :: S(nrobs,nrens)
   real,    intent(in)  :: R(nrobs,nrobs)
   real,    intent(out) :: W(nrobs,nrmin)
   real,    intent(out) :: eig(nrmin)
   real,    intent(in)  :: truncation

   real U0(nrobs,nrmin),sig0(nrmin)
   real B(nrmin,nrmin),Z(nrmin,nrmin)
   integer i,j

! Compute SVD of S=HA`  ->  U0, sig0
   call  svdS(S,nrobs,nrens,nrmin,U0,sig0,truncation)

! Compute B=sig0^{-1} U0^T R U0 sig0^{-1}
   call lowrankCee(B,nrmin,nrobs,nrens,R,U0,sig0)

! Compute eigenvalue decomposition  of B(nrmin,nrmin)
   call eigC(B,nrmin,Z,eig)

!   print *,'eig:',nrmin
!   print '(6g11.3)',eig

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Compute inverse diagonal of (I+Lamda)
   do i=1,nrmin
      eig(i)=1.0/(1.0+eig(i))
   enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! W = U0 * sig0^{-1} * Z
   do j=1,nrmin
   do i=1,nrmin
      Z(i,j)=sig0(i)*Z(i,j)
   enddo
   enddo

   call dgemm('n','n',nrobs,nrmin,nrmin, 1.0,U0,nrobs, Z,nrmin, 0.0,W,nrobs)

end subroutine




subroutine lowrankCee(B,nrmin,nrobs,nrens,R,U0,sig0)
implicit none
integer, intent(in) :: nrmin
integer, intent(in) :: nrobs
integer, intent(in) :: nrens
real, intent(inout) :: B(nrmin,nrmin)
real, intent(in)    :: R(nrobs,nrobs)
real, intent(in)    :: U0(nrobs,nrmin)
real, intent(in)    :: sig0(nrmin)
real X0(nrmin,nrobs)
integer  i,j

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Compute B=sig0^{-1} U0^T R U0 sig0^{-1}

! X0= U0^T R
   call dgemm('t','n',nrmin,nrobs,nrobs, 1.0,U0,nrobs, R,nrobs, 0.0,X0,nrmin)

! B= X0 U0
   call dgemm('n','n',nrmin,nrmin,nrobs, 1.0,X0,nrmin, U0,nrobs, 0.0,B,nrmin)

   do j=1,nrmin
   do i=1,nrmin
      B(i,j)=sig0(i)*B(i,j)
   enddo
   enddo

   do j=1,nrmin
   do i=1,nrmin
      B(i,j)=sig0(j)*B(i,j)
   enddo
   enddo

   B=real(nrens-1)*B

end subroutine


subroutine svdS(S,nrobs,nrens,nrmin,U0,sig0,truncation)
   integer, intent(in)  :: nrobs
   integer, intent(in)  :: nrens
   integer, intent(in)  :: nrmin
   real,    intent(in)  :: S(nrobs,nrens)
   real,    intent(out) :: sig0(nrmin)
   real,    intent(in)  :: U0(nrobs,nrmin)
   real,    intent(in)  :: truncation

   real S0(nrobs,nrens)
   real VT0(1,1)
   integer ierr
   integer lwork
   real, allocatable, dimension(:)   :: work
   integer nrsigma,i

   real sigsum,sigsum1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Compute SVD of S=HA`  ->  U0, sig0
   lwork=2*max(3*nrens+nrobs,5*nrens)
   allocate(work(lwork))

   S0=S
   sig0=0.0
   call dgesvd('S', 'N', nrobs, nrens, S0, nrobs, sig0, U0, nrobs, VT0, nrens, work, lwork, ierr)
   deallocate(work)
   if (ierr /= 0) then
      print *,'svdS: ierr from call dgesvd 0= ',ierr; stop
   endif

   sigsum=0.0
   do i=1,nrmin
      sigsum=sigsum+sig0(i)**2
   enddo

   sigsum1=0.0
! Significant eigenvalues.
   nrsigma=0
   do i=1,nrmin                       
      if (sigsum1/sigsum < truncation) then
         nrsigma=nrsigma+1
         sigsum1=sigsum1+sig0(i)**2
      else
         sig0(i:nrmin)=0.0
         exit
      endif
   enddo

   write(*,'(a,i5,g13.5)') '      analysis svdS: dominant sing. values and share ',nrsigma,sigsum1/sigsum
!   write(*,'(5g11.3)')sig0

   do i=1,nrsigma
       sig0(i) = 1.0/sig0(i)
   enddo

end subroutine


end module
