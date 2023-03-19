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

   nrmin=min(nrobs,nre)

   print *,'preexact transformation nrobs,nrens,nrmin',nrobs,nrens,nrmin

   E0=E
   sigma=0.0
   lwork=2*max(3*nre+nrobs,5*nre)
   allocate(work(lwork))
   sigma=0.0
   call dgesvd('S', 'N', nrobs, nre, E0, nrobs, sigma, U, nrobs, VT, nre, work, lwork, ierr)
   deallocate(work)
   if (ierr /= 0) then
      print *,'svdS: ierr from call dgesvd 0= ',ierr; stop
   endif

   S0=matmul(transpose(U),S)
   D0=matmul(transpose(U),D)
   do m=1,nrmin
      S(m,:)=(1.0/sigma(m))*S0(m,:)
      D(m,:)=(1.0/sigma(m))*D0(m,:)
   enddo
!   if (nrmin < nrobs) then
!      do m=nrmin,nrobs
!         S(m,:)=0.0
!         D(m,:)=0.0
!      enddo
!   endif

end subroutine
end module
