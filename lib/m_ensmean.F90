module m_ensmean
contains
subroutine ensmean(A,ave,nx,nrens)
   implicit none
   integer, intent(in) :: nx
   integer, intent(in) :: nrens
   real, intent(in)  :: A(nx,nrens)
   real, intent(out) :: ave(nx)
   integer j

   ave(:)=A(:,1)
   do j=2,nrens
      ave(:)=ave(:)+A(:,j)
   enddo
   ave=(1.0/real(nrens))*ave

end subroutine ensmean
end module m_ensmean
