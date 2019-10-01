module m_getD
contains
subroutine getD(D,subD,nrobs,nrens,lobs,nobs)
! Returns the subD matrix corresponding to active measurements
   implicit none
   integer, intent(in)  :: nrobs
   integer, intent(in)  :: nrens
   integer, intent(in)  :: nobs
   real,    intent(in)  :: D(nrobs,nrens)
   logical, intent(in)  :: lobs(nrobs)
   real,    intent(out) :: subD(nobs,nrens)

   integer j,m

   j=0
   do m=1,nrobs
      if (lobs(m)) then
         j=j+1
         subD(j,:)=D(m,:)
      endif
   enddo

end subroutine getD
end module m_getD
