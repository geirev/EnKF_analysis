module m_ensvar
contains
subroutine ensvar(A,ave,var,nx,nrens)
   implicit none
   integer, intent(in) :: nx
   integer, intent(in) :: nrens
   real, intent(in)    :: A(nx,nrens)
   real, intent(in)    :: ave(nx)
   real, intent(out)   :: var(nx)
   integer j

   var=0.0
   do j=1,nrens
      var(:)=var(:)+(A(:,j)-ave(:))*(A(:,j)-ave(:))
   enddo
   var=(1.0/real(nrens-1))*var

end subroutine ensvar
end module m_ensvar
