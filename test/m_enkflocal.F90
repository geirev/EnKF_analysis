module m_enkflocal
contains
subroutine enkflocal(mem,nx,nrens,obs,obsvar,obspos,nrobs,fixsamp,nre,nrr,mode_analysis,truncation,covmodel,dx,rh,Rexact,&
               &lrandrot,lupdate_randrot,lsakov,&
               &linflate, ladapinf, infmult,&
               &llocal,robs,distance_based_loc, adaptive_based_loc, obstreshold)

   use m_obs_pert
   use m_getD
   implicit none
   integer, intent(in) :: nx
   integer, intent(in) :: nrens
   integer, intent(in) :: nrobs
   integer, intent(in) :: nre
   integer, intent(in) :: nrr
   integer, intent(in) :: mode_analysis

   real,    intent(inout) :: mem(nx,nrens)
   real,    intent(in) :: obs(nrobs)
   real,    intent(in) :: obsvar
   integer, intent(in) :: obspos(nrobs)
   logical, intent(in) :: fixsamp
   logical, intent(in) :: Rexact

   logical, intent(in) :: lrandrot
   logical, intent(in) :: lupdate_randrot
   logical, intent(in) :: lsakov

   logical, intent(in) :: llocal
   real,    intent(in) :: robs          ! influence radii for the measurements
   logical, intent(in) :: distance_based_loc
   logical, intent(in) :: adaptive_based_loc
   real,    intent(in) :: obstreshold

   logical, intent(in) :: linflate
   logical, intent(in) :: ladapinf
   real, intent(in)    :: infmult

   real,    intent(in) :: truncation
   character(len=100), intent(in) :: covmodel
   real,    intent(in) :: dx
   real,    intent(in) :: rh

   real, allocatable :: R(:,:)
   real, allocatable :: E(:,:)
   real, allocatable :: D(:,:) 
   real, allocatable :: S(:,:)
   real, allocatable :: meanS(:)
   real, allocatable :: innovation(:)



   integer iens,m,i,j

! Local analysis variables
   integer l,icall
   logical                lobs(nrobs)           ! which measurements are active
   integer nobs
   logical local_rot
   real corr(nrobs),stdA,aveA,stdS
   real, allocatable, dimension(:,:) :: subS,subE,subD, subR
   real, allocatable, dimension(:,:)   :: submem
   real, allocatable :: subinnovation(:)


   print *,'EnKFlocal: ',llocal

! End Local analysis variables

   allocate(E(nrobs,nrens))
   allocate(D(nrobs,nrens))
   allocate(S(nrobs,nrens))
   allocate(meanS(nrobs))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Observe ensemble to construct the matrix S=HA
   do iens =1, nrens
      do m =1, nrobs
         S(m,iens)      =  mem(obspos(m),iens)
      enddo
   enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Construct observation perturbations E
   call obs_pert(E,nrens,nrobs,fixsamp,nre,nrr,dx,rh,covmodel,obspos)

! Introduce correct variances
   do iens=1,nrens
      do m=1,nrobs
         E(m,iens)=sqrt(obsvar)*E(m,iens)
      enddo
   enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Construct ensemble of measurements D=d+E
   do iens=1,nrens
      do m=1,nrobs
         D(m,iens)=obs(m)+E(m,iens)
      enddo
   enddo


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Compute innovation D'=D-HA
   D=D-S

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Compute mean(HA)
   meanS=0.0
   do iens=1,nrens
   do m=1,nrobs
      meanS(m)=meanS(m)+S(m,iens)
   enddo
   enddo
   meanS=(1.0/float(nrens))*meanS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Compute HA'=HA-mean(HA)
   do iens=1,nrens
      S(:,iens)=S(:,iens)-meanS(:)
   enddo


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   allocate(R(nrobs,nrobs))
   R=0.0

   if (Rexact) then
      print *,'   enkf: Exact R using covariance model: ',trim(covmodel)
      select case (trim(covmodel))
      case ('diagonal')
         do m=1,nrobs
            R(m,m)=obsvar
         enddo
      case ('gaussian')
         do i=1,nrobs
         do j=1,nrobs
            R(i,j)=obsvar*exp(-real(obspos(i)-obspos(j))**2/20.0**2)
         enddo
         enddo
      case default
         print *,'Covmodel is invalid : ',trim(covmodel)
      end select
   else
      print *,'   enkf: Lowrank R using covariance model: ',trim(covmodel)
      R=matmul(E,transpose(E))/float(nrens)
   endif


   allocate(innovation(nrobs))
   do m =1, nrobs
      innovation(m)      = obs(m)- meanS(m)
   enddo


   if (llocal) then
      print *,'EnKFlocal : ',llocal, obspos(:)
!      open(11,file='meascorr.dat',status='REPLACE')
      icall=0
      do i=1,nx
         nobs=0
         lobs(:)=.false.


! Computing correlation functions etc
         stdA=0.0
         aveA=0.0
         do l=1,nrens
            aveA=aveA+mem(i,l)
            stdA=stdA+mem(i,l)*mem(i,l)
         enddo
         aveA=aveA/real(nrens)
         stdA=stdA/real(nrens)
         stdA=stdA-aveA**2
         stdA=sqrt(stdA)

         do m=1,nrobs
            corr(m)=0.0
            stdS=0.0
            do l=1,nrens
               stdS=stdS+S(m,l)*S(m,l)
               corr(m)=corr(m)+S(m,l)*mem(i,l)
            enddo
            stdS=sqrt(stdS/real(nrens))
            corr(m)=corr(m)/real(nrens)
            corr(m)=corr(m)/(stdS*stdA)
         enddo

         if (distance_based_loc) then
            do m=1,nrobs
               if (real(abs(obspos(m)-i)) < robs) then
                  lobs(m)=.true.
                  nobs=nobs+1
               endif
            enddo
         endif

         if(adaptive_based_loc) then
            do m=1,nrobs
               if (abs(corr(m)) > obstreshold) then
                  lobs(m)=.true.
                  nobs=nobs+1
               endif
            enddo
         endif

!         print '(2(a,i5),a,4l1,a,4f10.4)','i=',i,' nobs=',nobs,' lobs=',lobs(:),' corr=',corr(:)
!        write(11,'(i5,4f10.4)')i,corr(:)


         if (nobs > 0) then
            allocate(subD(nobs,nrens))
            allocate(subE(nobs,nrens))
            allocate(subS(nobs,nrens))
            allocate(subR(nobs,nobs))
            call getD(D,subD,nrobs,nrens,lobs,nobs) ! the innovations to use 
            call getD(E,subE,nrobs,nrens,lobs,nobs) ! the observation errors to use
            call getD(S,subS,nrobs,nrens,lobs,nobs) ! the HA' to use
            subR=matmul(subE,transpose(subE))/float(nrens)
            allocate(subinnovation(nobs))
            allocate(submem(1,nrens))
            l=0
            do m =1, nrobs
               if (lobs(m)) then
                  l=l+1
                  subinnovation(l)      = obs(m)- meanS(m)
               endif
            enddo


            submem(1,:)=mem(i,:)
            icall=icall+1
            if (lupdate_randrot .and. icall==1) then
                local_rot=.true.
            else
                local_rot=.false.
            endif
            call analysis(submem, subR, subE, subS, subD, subinnovation, 1, nrens, nobs, .false., truncation, mode_analysis, &
                         lrandrot, local_rot, 1, './', lsakov, linflate, ladapinf, infmult)
            mem(i,:)=submem(1,:)
            deallocate(subD, subE, subS, subR, subinnovation, submem)
         endif
      enddo
!      close(11)
!      pause
   else
      print *,'   enkf: calling global analysis with mode: ',mode_analysis
      call analysis(mem, R, E, S, D, innovation, nx, nrens, nrobs, .true., truncation, mode_analysis, &
                      lrandrot, lupdate_randrot, 1, './', lsakov, linflate, ladapinf, infmult)

   endif

   deallocate(innovation)
   deallocate(R)


end subroutine
end module
