module m_dumpensemble
contains
subroutine dumpensemble(mem,ave,var,nrens,nx,ic,cc,nc)
   integer, intent(in) :: nrens
   integer, intent(in) :: nx
   integer, intent(in) :: nc
   real, intent(in)    :: mem(nx,nrens)
   real, intent(in)    :: ave(nx,nc)
   real, intent(in)    :: var(nx,nc)
   integer, intent(in) :: ic
   character(len=12), intent(in) :: cc(nc)
   integer, parameter :: nrensmax=100

   integer :: nn,j,k
   nn=min(nrensmax,nrens)
   
   print '(a,a,i2)','Dumping ensemble:',cc(ic)
   open(10,file='ensemble_'//trim(cc(ic))//'.dat')
      write(10,*)'TITLE = "ensemble"'
      write(10,*)'VARIABLES = "i-index" "ave" "var"'
      write(10,'(30(a,i3,a))')(' "',j,'"',j=1,nn)
      write(10,'(a,i5,a)')' ZONE  F=BLOCK, I=',nx
      write(10,'(30I5)')(k,k=1,nx)
      write(10,900)(ave(k,ic),k=1,nx)
      write(10,900)(var(k,ic),k=1,nx)
      do j=1,nn
         write(10,900)(mem(k,j),k=1,nx)
      enddo
  900 format(10(1x,e12.5))
   close(10)

end subroutine
end module
