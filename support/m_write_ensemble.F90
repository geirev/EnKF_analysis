module m_write_ensemble
! Writes the ensemble held in A to ensembleA.uf
contains
subroutine write_ensemble(A,nrens)
   use mod_dimensions
   use mod_states
   implicit none
   integer,      intent(in)  :: nrens
   type(states), intent(in)  :: A(nrens)      ! Ensemble matrix

   type(states4) A4

   integer reclA
   integer j

   character(len=9), parameter :: rident='HYCOM_1.2'

   inquire(iolength=reclA)rident,A4
   open(10,file='ensembleA.uf',form='unformatted',access='direct',recl=reclA)
      do j=1,nrens
         A4=A(j)
         print *,'Writing ensemble member: ',j
         write(10,rec=j)rident,A4
      enddo
   close(10)
end subroutine write_ensemble
end module m_write_ensemble

