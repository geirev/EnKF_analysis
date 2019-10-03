module m_obs_pert
contains
subroutine obs_pert(E,nrens,nrobs,fixsamp,nre,nrr,dx,rh,covmodel,obspos)
   use mod_dimensions
   use m_random
   use m_randrot
   use m_pseudo1D
   use m_fixsample1D
   implicit none
   integer, intent(in)    :: nrobs
   integer, intent(in)    :: nrens
   integer, intent(in)    :: nre
   integer, intent(in)    :: nrr
   real,    intent(out)   :: E(nrobs,nrens)
   logical, intent(in)    :: fixsamp
   real,    intent(in)    :: dx
   real,    intent(in)    :: rh
   character(len=8),  intent(in)    :: covmodel
   integer,    intent(in)    :: obspos(nrobs)

   real, allocatable :: work(:)
   real, allocatable :: EE(:,:)
   real, allocatable :: EEfield(:,:)
   integer iens,m

   integer ns,msx,nsx,i,j,lwork,ierr
   real, allocatable :: U(:,:),VT(:,:),VT1(:,:),sig(:)

   ns=nre*nrens
   msx=min(ns,nrobs)
   nsx=min(nrens,nrobs)

! Start with oversized ensemble
   allocate (EE(nrobs,ns))
   if (trim(covmodel) == 'diagonal') then
      call random(EE,nrobs*ns)   
      if (fixsamp) call fixsample1D(EE,nrobs,ns)

   elseif(trim(covmodel) == 'gaussian') then
     
      allocate(EEfield(nx,ns))
      call pseudo1D(EEfield,nx,nrens,rh,dx,nx)

      do iens=1,ns
      do m=1,nrobs
         EE(m,iens)=EEfield(obspos(m),iens)
      enddo
      enddo
      deallocate(EEfield)
   else
      print *,'Problem with covmodel:+++',trim(covmodel),'+++'
      stop
   endif

   if (nre == 1) then
      do iens=1,nrens
      do m=1,nrobs
         E(m,iens)=EE(m,iens)
      enddo
      enddo

   else
      allocate(VT1(nsx,nsx))
      call randrot(VT1,nsx)

! Compute SVD of oversized ensemble
      lwork=2*max(3*ns+max(nrobs,ns),5*ns)
      allocate(work(lwork))
      allocate( U(nrobs,msx), sig(msx), VT(msx,msx) )
      call dgesvd('S', 'N', nrobs, ns, EE, nrobs, sig, U, nrobs, VT, ns, work, lwork, ierr)
      if (ierr /= 0) print *, 'ierr',ierr
      deallocate(work)

! Generate first nrens or nsx members
      E=0.0
      do j=1,nsx
         do i=1,nsx
            E(:,j)=E(:,j)+U(:,i)*sig(i)/sqrt(float(nre))*VT1(i,j)
         enddo
      enddo
      deallocate(U, VT, sig, VT1)

      if (fixsamp) call fixsample1D(E,nrobs,nrens)
   endif

   deallocate(EE)

end subroutine obs_pert
end module m_obs_pert
