module m_newton1D
contains
subroutine newton1D(r1,n1,dx,rx,lconv)
   use m_newtonfunc1D
   implicit none
   integer                    n1
   real,    intent(inout)  :: r1
   real,    intent(in)     :: dx
   real,    intent(in)     :: rx

   logical lconv
   real, parameter :: eps=1.0E-05
   real inc1,err1
   integer i,j
   real f,f1,r1ini,gamma

   lconv=.false.
   r1ini=r1
   gamma=1.25

   do j=1,10
      r1=real(j)*r1ini
      gamma=gamma-0.25/real(j)
      do i=1,100

         call newtonfunc1D(f,f1,r1,n1,dx,rx)

         inc1 =  f/f1

         r1 = r1 - gamma*inc1

         r1=max(r1,0.000001)
         r1=min(r1,1.0)

         err1 = inc1/(abs(r1)+eps)

         if (abs(err1) < eps) then
!            print *,' Newton converged in iteration ',j,i
            lconv=.true.
            exit
         endif

      enddo
      if (lconv) exit
   enddo

   if (.not. lconv) then
      print *,'Newton did not converge'
   endif

end subroutine newton1D
end module m_newton1D




