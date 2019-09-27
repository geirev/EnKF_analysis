module m_newtonfunc1D
contains
subroutine newtonfunc1D(f,f1,r1,n1,dx,rx)
   implicit none
   integer                 n1
   real,    intent(out) :: f,f1
   real,    intent(in)  :: r1
   real,    intent(in)  :: dx
   real,    intent(in)  :: rx
   real, parameter   :: pi=3.141592653589

   integer l
   real pi2,e,kappa,kappa2

   pi2=2.0*pi

   kappa=pi2/(real(n1)*dx)
   kappa2=kappa**2

   f=0.0; f1=0.0

   do l=-n1/2+1,n1/2
      if (l == 0 ) cycle

      e=exp( -2.0*( kappa2*real(l*l)/r1**2 ) )

      f=f + e * ( cos(kappa *real(l)*rx) - exp(-1.0) )

      f1=f1 + e * (4.0*kappa2 *real(l*l)/r1**3) * ( cos(kappa *float(l)*rx) - exp(-1.0) )
      
   enddo

end subroutine newtonfunc1D
end module m_newtonfunc1D
