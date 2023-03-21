module m_residuals
contains
subroutine residuals(ave,var,cc,nx,nc)
   integer, intent(in)           :: nx
   integer, intent(in)           :: nc
   real, intent(in)              :: ave(nx,nc)
   real, intent(inout)           :: var(nx,nc)
   character(len=12), intent(in) :: cc(nc)
   real norm    
   integer ic,jc
   real resave(nc,nc)
   real resstd(nc)
   real std(nx,nc)
   std=sqrt(var)

   print '(a)','++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
   print '(a)','Dumping resuduals between methods etc'
   resave=0.0
   resstd=0.0
   var=sqrt(var)
   norm=sqrt(dot_product(ave(1:nx,1),ave(1:nx,1))/real(nx))
   do jc=1,nc
   do ic=1,nc
      resave(ic,jc)=sqrt(dot_product(ave(1:nx,ic)-ave(1:nx,jc),ave(1:nx,ic)-ave(1:nx,jc))/real(nx))
   enddo
   resstd(jc)=sqrt(dot_product(std(1:nx,jc),std(1:nx,jc))/real(nx))
   enddo

   write(*,'(a)')'residuals of averages:'
!   write(*,'(tr12,100i12)')(jc,jc=1,nc)
   write(*,'(tr23,100a12)')(cc(jc),jc=1,nc)
   do jc=1,nc
      write(*,'(i5,tr2,a12,100f12.6)')jc,cc(jc),resave(1:jc-1,jc)/norm
   enddo

   write(*,'(a)')'Standard deviations:'
   write(*,'(tr19,100f12.6)')resstd(1:nc)
end subroutine
end module
