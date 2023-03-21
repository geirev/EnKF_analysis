module m_preexact
contains
subroutine preexact(E,S,D,nrobs,nrens,nre)
!   use mod_anafunc
   implicit none
   integer, intent(in) :: nrobs
   integer, intent(in) :: nrens
   integer, intent(in) :: nre
   integer nrmin

   real,    intent(in)    :: E(nrobs,nre)
   real,    intent(inout) :: D(nrobs,nrens)
   real,    intent(inout) :: S(nrobs,nrens)

   integer m
   real E0(nrobs,nre)
   real U(nrobs,min(nrobs,nre))
   real sigma(min(nrobs,nre))
   real S0(min(nrobs,nre),nrens)
   real D0(min(nrobs,nre),nrens)

   real VT(1,1)
   integer ierr
   integer lwork
   real, allocatable, dimension(:)   :: work

   real sigacc,sigsum
   integer nrsigma

   nrmin=min(nrobs,nre)

   print '(tr7,a,4(tr1,i0))','preexact: transformation nrobs,nrens,nre,nrmin',nrobs,nrens,nre,nrmin

   E0=E/sqrt(real(nre-1))
   sigma=0.0
   lwork=2*max(3*nre+nrobs,5*nre)
   allocate(work(lwork))
   sigma=0.0
   call dgesvd('S', 'N', nrobs, nre, E0, nrobs, sigma, U, nrobs, VT, nre, work, lwork, ierr)
   deallocate(work)
   if (ierr /= 0) then
      print '(tr7,a,i2)','prexact: ierr from call dgesvd 0= ',ierr; stop
   endif

   sigsum=sum(sigma)
   sigacc=0.0
   do m=1,nrmin
      sigacc=sigacc+sigma(m)
      if (sigacc/sigsum < 0.99999) then
         sigma(m)=1.0/sigma(m)
      else
         sigma(m)=0.0
      endif
   enddo

   nrsigma=nrmin
   do m=1,nrmin
      if (sigma(m)==0.0) then
         print '(a,I5)','       Number of sigma values: ',m-1
         nrsigma=m-1
         exit
      endif
   enddo
!   print '(a)','sigma:'
!   print '(10f10.4)',1.0/sigma(1:m-1)

!   S0=matmul(transpose(U),S)
   call dgemm('t','n',nrmin,nrens,nrobs,1.0,U,nrobs,S,nrobs,0.0,S0,nrmin)

!   D0=matmul(transpose(U),D)
   call dgemm('t','n',nrmin,nrens,nrobs,1.0,U,nrobs,D,nrobs,0.0,D0,nrmin)

   do m=1,nrmin
      S(m,:)=sigma(m)*S0(m,:)
      D(m,:)=sigma(m)*D0(m,:)
   enddo

   if (nrmin < nrobs) then
      do m=nrmin+1,nrobs
         S(m,:)=0.0
         D(m,:)=0.0
      enddo
   endif

!   print '(a)','S0:'
!   do m=1,nrsigma
!      print '(i5,10f10.4)',m,S(m,1:10)
!   enddo
!   print '(a)','D0:'
!   do m=1,nrsigma
!      print '(i5,10f10.4)',m,D(m,1:10)
!   enddo

end subroutine
end module
