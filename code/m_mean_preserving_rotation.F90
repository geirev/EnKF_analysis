module m_mean_preserving_rotation
contains
subroutine mean_preserving_rotation(Up,nrens)
! Generates the mean preserving random rotation for the EnKF SQRT algorithm
! using the algorithm from Sakov 2006-07.  I.e, generate rotation Up suceh that
! Up*Up^T=I and Up*1=1 (all rows have sum = 1)  see eq 17.
! From eq 18,    Up=B * Upb * B^T 
! B is a random orthonormal basis with the elements in the first column equals 1/sqrt(nrens)
! Upb = | 1  0 |
!       | 0  U |
! where U is an arbitrary orthonormal matrix of dim nrens-1 x nrens-1  (eq. 19)


use m_randrot

implicit none

integer, intent(in)    :: nrens
real,    intent(out)   :: Up(nrens,nrens)

real B(nrens,nrens),Q(nrens,nrens),R(nrens,nrens),U(nrens-1,nrens-1), Upb(nrens,nrens)

integer j,k




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Generating the B matrix
! Starting with a random matrix with the correct 1st column
   call random_number(B)
   B(:,1)=1.0/sqrt(real(nrens))

! modified_gram_schmidt is used to create the orthonormal basis
!   do k=1,nrens
!      R(k,k)=sqrt(dot_product(B(:,k),B(:,k)))
!      Q(:,k)=B(:,k)/R(k,k)
!      do j=k+1,nrens
!         R(k,j)=dot_product(Q(:,k),B(:,j))
!         B(:,j)=B(:,j)- Q(:,k)*R(k,j)
!      enddo
!   enddo
!   B=Q

! with overwriting of B
   do k=1,nrens
      R(k,k)=sqrt(dot_product(B(:,k),B(:,k)))
      B(:,k)=B(:,k)/R(k,k)
      do j=k+1,nrens
         R(k,j)=dot_product(B(:,k),B(:,j))
         B(:,j)=B(:,j)- B(:,k)*R(k,j)
      enddo
   enddo


! Check on orthonormality of B
!   do k=1,nrens
!      do j=k,min(k+14,nrens)
!         write(*,'(15f10.4)',advance='no')dot_product(B(:,j),B(:,k))
!      enddo
!      write(*,*)' '
!   enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Creating the orthonormal nrens-1 x nrens-1 U matrix
   call randrot(U,nrens-1)
!! Check on orthonormality of U
!   do k=1,nrens-1
!      do j=k,min(k+14,nrens-1)
!         write(*,'(15f10.4)',advance='no')dot_product(U(:,j),U(:,k))
!      enddo
!      write(*,*)' '
!   enddo
!   stop

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Creating the orthonormal nrens x nrens Upb matrix
  Upb(2:nrens,2:nrens)=U(1:nrens-1,1:nrens-1)
  Upb(1,1)=1.0
  Upb(2:nrens,1)=0.0
  Upb(1,2:nrens)=0.0

! Check on orthonormality of Upb
!   do k=1,nrens
!      do j=k,min(k+14,nrens)
!         write(*,'(15f10.4)',advance='no')dot_product(Upb(:,j),Upb(:,k))
!      enddo
!      write(*,*)' '
!   enddo
!   stop

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Creating the random orthonormal mean preserving nrens x nrens Upb matrix: Up=B^T Upb B
   call dgemm('n','n',nrens,nrens,nrens,1.0,B,nrens,Upb,nrens,0.0,Q,nrens)
   call dgemm('n','t',nrens,nrens,nrens,1.0,Q,nrens,B,nrens,0.0,Up,nrens)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Checks
!   do k=1,nrens
!      print *,'Up row sum: ',k,sum(Up(k,1:nrens))
!   enddo

! Check on orthonormality of Up
!   do k=1,nrens
!      do j=k,min(k+14,nrens)
!         write(*,'(15f10.4)',advance='no')dot_product(Up(:,j),Up(:,k))
!      enddo
!      write(*,*)' '
!   enddo
!   stop


end subroutine 
end module
