module m_tecsol
implicit none
contains
subroutine tecsol(ave,var,obs,obspos,nx,nrobs,nc)
   integer, intent(in)   :: nx
   integer, intent(in)   :: nrobs
   integer, intent(in)   :: nc
   real,    intent(in)   :: ave(nx,nc)
   real,    intent(in)   :: var(nx,nc)
   integer, intent(in)   :: obspos(nrobs)
   real,    intent(in)   :: obs(nrobs)
   integer ic,i
   
   print '(a)','Dumping solutions to solutions.dat'
   open(10,file='solutions.dat')
   write(10,*)'TITLE = "Solutions"'
   write(10,*)'VARIABLES = "i" "Truth" "First guess" "Prior" "10" "11" "12" "13" "21" "22" "23"'
   write(10,'(a,i5,a,i5,a)')' ZONE T="Average"  F=BLOCK, I=',nx,', J=1, K=1'
   write(10,'(20I5)')(i         ,i=1,nx)
   do ic=1,nc
      write(10,'(20g13.5)')(ave(i,ic) ,i=1,nx)
   enddo

   write(10,'(a,i5,a,i5,a)')' ZONE T="Std Dev"  F=BLOCK, I=',nx,', J=1, K=1'
   write(10,'(20I5)')(i         ,i=1,nx)
   do ic=1,nc
      write(10,'(20g13.5)')(var(i,ic) ,i=1,nx)
   enddo

   write(10,'(a,i5,a,i5,a)')' ZONE T="observations"  F=BLOCK, I=',nrobs,', J=1, K=1'
   write(10,'(20I5)')(obspos(i) ,i=1,nrobs)
   do ic=1,nc
      write(10,'(20g13.5)')(obs(i) ,i=1,nrobs)
   enddo
   close(10)
end subroutine
end module

