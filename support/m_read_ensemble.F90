module m_read_ensemble
! Reads the modelled forecast from ensembleF.uf and stores it in the matrix A
contains
subroutine read_ensemble(A,nrens)
   use mod_dimensions
   use mod_states
   implicit none
   integer,      intent(in)  :: nrens
   type(states), intent(out) :: A(nrens)      ! Ensemble matrix

   type(states4) A4

   character(len=9) rident  
   integer reclA
   integer iostat
   integer j

   logical ex

   character(len=9), parameter :: cident4='HYCOM_1.2'

   inquire(iolength=reclA)rident,A4
   open(10,file='ensembleF.uf',form='unformatted',access='direct',recl=reclA)
      read(10,rec=1,err=90,iostat=iostat)rident    
      if (rident == cident4) then
         do j=1,nrens
            read(10,rec=j)rident,A4
            print *,'Reading ensemble member: ',j
            A(j)=A4
         enddo
      else
         print *,'get_mod_forecast ERROR: wrong ident of restart file: ',rident
         stop
      endif
   close(10)

   return
   90 stop 'read_ensemble: error reading ensembleF.uf '
end subroutine read_ensemble
end module m_read_ensemble
