module m_randrot
contains
subroutine randrot(Q,nrens)
   implicit none
   integer, intent(in)  :: nrens
   real,    intent(out) :: Q(nrens,nrens)

   real, dimension(nrens,nrens) ::  A, B
   real sigma(nrens), work(10*nrens)
   real, parameter :: pi=3.14159253589
   integer ierr
   real meanB

   call random_number(B)
   call random_number(A)
   Q = sqrt(-2.*log(A+tiny(A))) * cos(2.*pi*B)

!$OMP CRITICAL
! QR factorization
   call dgeqrf(nrens, nrens, Q, nrens, sigma, work, 10*nrens, ierr )
   if (ierr /= 0) print *, 'randrot: dgeqrf ierr=',ierr

! Construction of Q
   call dorgqr(nrens, nrens, nrens, Q, nrens, sigma, work, 10*nrens, ierr )
   if (ierr /= 0) print *, 'randrot: dorgqr ierr=',ierr
!$OMP END CRITICAL


end subroutine randrot
end module m_randrot

