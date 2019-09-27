module m_measurements
contains
subroutine measurements(ana,nx,obs,obspos,nrobs,obsvar,iobs,mkobs,time)
   use m_pseudo1D
   use m_random
   implicit none
   integer, intent(in) :: nx
   real,    intent(in) :: ana(nx)
   integer, intent(in) :: nrobs
   integer, intent(in) :: obspos(nrobs)
   real,    intent(in) :: obsvar
   real,    intent(out):: obs(nrobs)
   integer, intent(in) :: iobs
   logical, intent(in) :: mkobs
   real,    intent(in) :: time

   integer m,i,reclA
   real tt
   logical ex

   inquire(iolength=reclA)tt,i,obs
   inquire(file='obs.uf',exist=ex)

   if (.not. mkobs) then
      if (ex) then
         open(10,file='obs.uf',form='unformatted',access='direct',recl=reclA)
            read(10,rec=iobs,err=100)tt,i,obs
         close(10)
         if (tt /= time .or. i /= nrobs) then
            print *,'Problem reading obs.uf :',iobs,tt,time,i,nrobs
         endif
      else
         print *,'file obs.uf does not exist'
         stop
      endif

    elseif (mkobs) then
         
      call random(obs,nrobs)

      do m=1,nrobs
         obs(m)= ana(obspos(m)) + sqrt(obsvar)*obs(m)
      enddo

      open(10,file='obs.uf',form='unformatted',access='direct',recl=reclA)
         write(10,rec=iobs)time,nrobs,obs
      close(10)
   endif

   return
   100 stop 'Error reading obs.uf (record does not exist?)'
end subroutine measurements
end module m_measurements


