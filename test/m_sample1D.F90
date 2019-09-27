module m_sample1D
! This routine samples pseudo random fields with improved coonditioning.
! This is done by first drawing a large sample and then construct the
! final sample using the dominant singular vectors of the large sample.

contains
subroutine sample1D(A2,n,nrens,nre,nrr,dx,rh,samp_fix,periodic)
   use m_pseudo1D
   use m_fixsample1D
   use m_randrot
   use mod_anafunc
   implicit none
   integer, intent(in)     ::  n                 ! Dimension of state vector
   integer, intent(in)     ::  nrens             ! Number of realizations to be simulted
   integer, intent(in)     ::  nre               ! Size of starting ensemble is ns=nre*nrens
   integer, intent(in)     ::  nrr               ! Number of singular vectors to use for sampling nsz=nrr*nrens
   real,    intent(in)     ::  rh
   real,    intent(in)     ::  dx
   logical, intent(in)     ::  samp_fix
   logical, intent(in)     ::  periodic
   real,    intent(inout)    ::  A2(n,nrens)

   integer ns,msx,i,j,nsx,n1,nrnn
   integer lwork,ierr
   real summ
   real, allocatable, dimension(:,:) :: A,Atmp,U,VT,VT1,CA,CE,UE,ZE,ZA
   real, allocatable, dimension(:)   :: sig,work,eigE,eigA

   logical :: debug=.false.
   logical :: diag=.false.
   logical leig

   if (periodic) then
      n1=n
   else
      n1=nint(real(n)*1.2)
   endif

#ifdef SGI
   if (mod(n1,2) == 1 ) n1=n1+1
#endif

#if defined(IBM) || defined(LINUX)
   do i=1,100
      if (2**i >= n1) then
         n1=2**i
         exit
      endif
   enddo
#endif
   if (periodic .and. n /= n1) then
      print '(a,i5,a,i5)','m_sampling1D: Modify n=',n,' to n1=',n1
      stop 'm_sample1D: You have to change model grid size for periodic samples'
   endif


   ns=nre*nrens
   msx=min(ns,n)
   nsx=min(nrr*nrens,n) 
   nrnn=min(nrens,n)


! Start with possibly oversized ensemble of ns members
   allocate(A(n,ns))
   call pseudo1D(A,n,ns,rh,dx,n1)

   if (nre > 1) then
      print '(4(a,I5))','Improved sampling : nre=',nre,' nrr=',nrr,' ns=',ns,' nsx=',nsx  
! Compute SVD of oversized ensemble
      lwork=2*max(3*ns+max(n,ns),5*ns)
      allocate(work(lwork), U(n,msx), sig(msx), VT(1,1) )
      call dgesvd('S', 'N', n, ns, A, n, sig, U, n, VT, msx, work, lwork, ierr)
      if (ierr /= 0) print *, 'ierr',ierr

! make an orthogonal VT1 used as linear combination for final ensemble
      allocate ( VT1(nsx,nsx) )
      call randrot(VT1,nsx)

! Generate first nrens members based nsx dominant singular vectors
      A2=0.0
      do j=1,nrens
         do i=1,nsx
            A2(:,j)=A2(:,j)+U(:,i)*sig(i)*sqrt(real(nsx)/real(ns))*VT1(i,j)
         enddo
      enddo
      deallocate(work, A, VT1, U, sig,VT)
   else
      print '(4(a,I5))','Standard sampling : nre=',nre,' nrr=',nrr,' ns=',ns,' nsx=',nsx 
      A2(:,1:nrens)=A(:,1:nrens)
      deallocate(A)
   endif

   print *,'sampling done'

! subtract mean and correct variance
   if (samp_fix) call fixsample1D(A2,n,nrens)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Diagnostic testing of sampling algorithms
   if (debug) then
! Generate an initial ensemble
      allocate(A(n,nrens),Atmp(n,nrens))
      A=A2

! Ensemble covariance matrix
      allocate(CE(n,n))
      call dgemm('n','t',n,n,nrens,1.0/real(nrens-1),A,n,A,n,0.0,CE,n)

!  Analytical covariance matrix 
      allocate(CA(n,n))
      do j=1,n
      do i=1,n
         CA(i,j)=exp(-(min( abs(i-j), min(i,j)+n-max(j,i) )*dx/rh)**2)
      enddo
      enddo

      print *,'C.dat'
      open(10,file='C.dat')
         write(10,*)'TITLE = "covariance functions"'
         write(10,*)'VARIABLES = "X" "CA" "CE"'
         do j=1,nrnn
         write(10,'(a,i3,a,i5,a)')' ZONE T="C(',j,')", F=POINT, I=',n,' J=1 K=1'
         do i=1,n
            write(10,'(i4,10g13.5)')i,CA(i,j),CE(i,j)
         enddo
         enddo
      close(10)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Singular- and eigenvectors
      leig=.true.

! SVD of ensemble
      lwork=2*max(3*nrnn+max(n,nrnn),5*nrnn)
      allocate(UE(n,nrnn), VT(1,1), sig(nrnn), work(lwork))
      Atmp=A
      call dgesvd('S', 'N', n, nrens, Atmp, n, sig, UE, n, VT, nrens, work, lwork, ierr)
      deallocate(VT, work)


! SVD of ensemble covariance
      allocate( ZE(n,n), eigE(n), VT(1,1), work(10*n) )
      if (leig) then
         call eigC(CE,n,ZE,eigE)
      else
         call dgesvd('S', 'N', n, n, CE, n, eigE, ZE, n, VT, n, work, 10*n, ierr)
         if (ierr /= 0) print *, 'ierr',ierr
      endif
      deallocate(VT,work)

! SVD of analytical covariance
      allocate( ZA(n,n), eigA(n), VT(1,1), work(10*n) )
      if (leig) then
         call eigC(CA,n,ZA,eigA)
      else
         call dgesvd('S', 'N', n, n, CA, n, eigA, ZA, n, VT, n, work, 10*n, ierr)
         if (ierr /= 0) print *, 'ierr',ierr
      endif
      deallocate(VT,work)


      print *,'sigma.dat'
      open(10,file='sigma.dat')
         write(10,*)'VARIABLES = "i" "sigE" "eigE" "eigA"'
         write(10,'(a,i5,a)')' ZONE T="Eigenvalues", F=POINT, I=',n,' J=1 K=1'
         summ=0.0
         do i=1,n
            summ=summ+sig(min(i,nrnn))**2
            if (leig) then
               write(10,'(i4,3e12.4)')i,sig(min(i,nrnn))**2/nrens,eigE(n-i+1),eigA(n-i+1)
            else
               write(10,'(i4,3e12.4)')i,sig(min(i,nrnn))**2/nrens,eigE(i),eigA(i)
            endif
         enddo
      close(10)
      deallocate(sig)

! Write all eigen and singular vectors
      print *,'all.dat'
      open(10,file='all.dat')
         write(10,*)'TITLE = "eigen and singular vectors "'
         write(10,*)'VARIABLES = "X" "UE" "ZE" "ZA"'
         do j=1,nrnn
         write(10,'(a,i3,a,i5,a)')' ZONE T="mode(',j,')", F=POINT, I=',n,' J=1 K=1'
         do i=1,n
            if (leig) then
               write(10,'(i4,10g13.5)')i,UE(i,j),ZE(i,n-j+1),ZA(i,n-j+1)
            else
               write(10,'(i4,10g13.5)')i,UE(i,j),ZE(i,j),ZA(i,j)
            endif
         enddo
         enddo
      close(10)

      deallocate(CA, CE, Atmp, A)

   endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


   if (diag) then
! SVD of new ensemble
      allocate (U(n,nrnn))
      allocate (sig(nrnn))
      allocate (VT(nrnn,nrnn))
      sig=0.0
      allocate( work(10*n) )
      call dgesvd('S', 'S', n, nrnn, A2, n,  sig,  U, n, VT, nrnn, work, 10*n, ierr)
      if (ierr /= 0) print *, 'ierr',ierr

      print *,'sigma2.dat'
      open(10,file='sigma2.dat')
         write(10,*)'TITLE = "Normalized singular values"'
         write(10,*)'VARIABLES = "i" "sigma" "sigma2" "accumulated"'
         summ=0.0
         do i=1,nrnn
            summ=summ+sig(i)**2
            write(10,'(i4,3e12.4)')i,sig(i)/sig(1),sig(i)**2/sig(1)**2,summ/real(n*ns)
         enddo
      close(10)
      deallocate(U, VT, sig, work)
      stop
      
   endif




end subroutine sample1D
end module m_sample1D
