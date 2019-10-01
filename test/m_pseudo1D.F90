module m_pseudo1D
contains
subroutine pseudo1D(A,nx,nrfields,rx,dx,n1)
   use m_random
   use m_newton1D
#ifdef LINUX
   use mod_fftw3
#endif
   implicit none
   integer, intent(in) :: nx                ! Dimension of random vector
   integer, intent(in) :: nrfields          ! Number of random vectors
   real, intent(out)   :: A(nx,nrfields)    ! The random vectors
   real, intent(in)    :: rx                ! characteristic lengthscale
   real, intent(in)    :: dx                ! delta x  ( >1 )
   integer, intent(in) :: n1              ! Dimension of random vector

   real r1,r12
   logical cnv
   integer l,i
   real c
   real kappa2,kappa
   real pi2,deltak,summ

   real fampl(0:n1/2,2)
   real phi(0:n1/2)


   real tt
   logical, save :: diag=.false.

   real, parameter :: pi=3.141592653589

#ifdef SGI
   integer, parameter :: sign=-1 
   real coeff(15+n1)
   complex arrayC(0:n1/2)
   real y(n1+2)
#endif

#ifdef LINUX
   integer(kind=8) plan
   complex arrayC(n1/2+1)
   real y(n1)
#endif

#ifndef SGI
#ifndef LINUX
   print *,'ranfield is only running on the following machines:'
   print *,'   SGI'
   print *,'   LINUX having FFTW3'
   stop
#endif
#endif

   pi2=2.0*pi
   deltak=pi2/(real(n1)*dx)
   kappa=pi2/(real(n1)*dx)
   kappa2=kappa**2

#ifdef LINUX
   call dfftw_plan_dft_c2r_1d(plan,n1,arrayC,y,FFTW_ESTIMATE)
#endif

#ifdef SGI
   call dzfft1dui(n1,coeff)
#endif

   r1=3.0/rx
   if (diag) print '(a,f13.5,i5,2f12.2)','Call newton1D with ',r1,n1,dx,rx
   call newton1D(r1,n1,dx,rx,cnv)
   if (.not.cnv) then
      stop 'pseudo1D: newton did not converge.  Recompile with diag set to true in pseudo1D.'
   endif
   if (diag) print *, 'Newton gave r1= ',r1



   r12=r1**2
   summ=0.0
   do l=1,n1/2
      summ=summ+2.0*exp(-2.0*(kappa2*real(l*l))/r12)
   enddo
   c=sqrt(1.0/(deltak*summ))

   if (diag) then
      print *,'pseudo1D: summ',summ
      print *,'pseudo1D: rx  ',rx
      print *,'pseudo1D: r1  ',r1
      print *,'pseudo1D: c=  ',c
   endif


   do i=1,nrfields
!     Calculating the random wave phases 
      call random_number(phi)
      phi=pi2*phi

      do l=0,n1/2 
         tt=kappa2*real(l*l)/r12
         fampl(l,1)=exp(-tt)*cos(phi(l))*sqrt(deltak)*c
         fampl(l,2)=exp(-tt)*sin(phi(l))*sqrt(deltak)*c
      enddo
      fampl(0,1)=0.0
      fampl(0,2)=0.0


      arrayC(:)=cmplx(fampl(:,1),fampl(:,2))
#ifdef SGI
      call zdfft1du(sign,n1,arrayC,1,coeff)
      if ( (n1/2+1)*2 == n1+1) then
         y(1:n1+1)=transfer(arrayC(0:n1/2),(/0.0/))
      elseif ( (n1/2+1)*2 == n1+2) then
         y(1:n1+2)=transfer(arrayC(0:n1/2),(/0.0/))
      else
         print *,'pseudo1D: n1 problem: ',n1,(n1/2+1)*2,n1+2
      endif
#endif

#ifdef LINUX
      call dfftw_execute(plan)
#endif
      A(1:nx,i)=y(1:nx)

   enddo

#ifdef LINUX
   call dfftw_destroy_plan(plan)
#endif

end subroutine pseudo1D
end module m_pseudo1D
