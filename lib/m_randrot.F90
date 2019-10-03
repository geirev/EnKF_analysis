module m_randrot
contains
subroutine randrot(Q,nrens)
! This routine generates a real orthogonal random matrix.
! The algorithm is the one by
!    Francesco Mezzadri (2007), How to generate random matrices from the classical
!    compact groups, Notices of the AMS, Vol. 54, pp 592-604.
! 1. First a matrix with independent random normal numbers are simulated.
! 2. Then the QR decomposition is computed, and Q will then be a random orthogonal matrix.
! 3. The diagonal elements of R are extracted and we construct the diagonal matrix X(j,j)=R(j,j)/|R(j,j)|
! 4. An updated Q'=Q X is computed, and this is now a random orthogonal matrix with a Haar measure.

   implicit none
   integer, intent(in)  :: nrens
   real,    intent(out) :: Q(nrens,nrens)

   real, dimension(nrens,nrens) :: A, B
   real, dimension(nrens)       :: diagR
   real sigma(nrens), work(10*nrens)
   real, parameter :: pi=3.14159253589
   integer i,ierr

   call random_number(B)
   call random_number(A)
   Q = sqrt(-2.*log(A+tiny(A))) * cos(2.*pi*B)

!$OMP CRITICAL
! QR factorization
   call dgeqrf(nrens, nrens, Q, nrens, sigma, work, 10*nrens, ierr )
   if (ierr /= 0) print *, 'randrot: dgeqrf ierr=',ierr

   do i=1,nrens
      diagR(i)=Q(i,i)/abs(Q(i,i))
      print '(a,2f12.4)','R(i,i), diagR(i) ',Q(i,i), diagR(i)
   enddo

! Construction of Q
   call dorgqr(nrens, nrens, nrens, Q, nrens, sigma, work, 10*nrens, ierr )
   if (ierr /= 0) print *, 'randrot: dorgqr ierr=',ierr
!$OMP END CRITICAL

   do i=1,nrens
      if (diagR(i) < 0) Q(:,i)=-Q(:,i)
   enddo


end subroutine randrot
end module m_randrot

