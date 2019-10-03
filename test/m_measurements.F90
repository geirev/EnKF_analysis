module m_measurements
use m_sample1D
contains
subroutine measurements(ana,nx,obs,obspos,nrobs,obsvar,covmodel,rd,dx)
   use m_pseudo1D
   use m_random
   implicit none
   integer, intent(in) :: nx
   real,    intent(in) :: ana(nx)
   integer, intent(in) :: nrobs
   integer, intent(in) :: obspos(nrobs)
   real,    intent(in) :: obsvar
   real,    intent(in) :: rd
   real,    intent(in) :: dx
   real,    intent(out):: obs(nrobs)
   character(len=8), intent(in) :: covmodel

   real :: obserr(nx)

   integer m

   if (covmodel=='gaussian' .and. (rd > 0.0)) then
      call pseudo1D(obserr,nx,1,rd,dx,nx)
   else
      call random(obserr,nx)
   endif

   do m=1,nrobs
      obs(m)= ana(obspos(m)) + sqrt(obsvar)*obserr(obspos(m))
   enddo

end subroutine measurements
end module m_measurements


