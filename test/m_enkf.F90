module m_enkf
contains
subroutine enkf(mem,nx,nrens,obs,obsvar,obspos,nrobs,mode_analysis,truncation,covmodel,dx,rh,rd,&
               &lrandrot,lsymsqrt,&
               &inflate, infmult,&
               &local,robs, obstreshold,E0)

   use m_obs_pert
   use m_getD
   implicit none
   integer, intent(in) :: nx
   integer, intent(in) :: nrens
   integer, intent(in) :: nrobs
   integer, intent(in) :: mode_analysis

   real,    intent(inout) :: mem(nx,nrens)
   real,    intent(in) :: obs(nrobs)
   real,    intent(in) :: obsvar
   integer, intent(in) :: obspos(nrobs)

   logical, intent(in) :: lrandrot
   logical, intent(in) :: lsymsqrt

   integer, intent(in) :: local
   real,    intent(in) :: robs          ! influence radii for the measurements
   real,    intent(in) :: obstreshold

   integer, intent(in) :: inflate
   real, intent(in)    :: infmult

   real,    intent(in) :: truncation
   character(len=8), intent(in) :: covmodel
   real,    intent(in) :: dx
   real,    intent(in) :: rh
   real,    intent(in) :: rd

   real,    intent(inout) :: E0(nrobs,nrens)   ! meas pert use to comapare between schemes

   real, allocatable :: R(:,:)
   real, allocatable :: E(:,:)
   real, allocatable :: D(:,:) 
   real, allocatable :: S(:,:)
   real, allocatable :: meanS(:)
   real, allocatable :: innovation(:)
   real, allocatable :: scaling(:)

   integer iens,m,i,j
   logical :: lupdate_randrot=.true.

! Local analysis variables
   integer l,icall
   logical lobs(nrobs)           ! which measurements are active
   integer nobs
   logical local_rot
   real corr(nrobs),stdA,aveA,stdS,dist
   real, allocatable, dimension(:,:) :: subS,subE,subD, subR
   real, allocatable, dimension(:,:)   :: submem
   real, allocatable :: subinnovation(:)
   real :: start, finish


! End Local analysis variables

   allocate(E(nrobs,nrens))
   allocate(D(nrobs,nrens))
   allocate(S(nrobs,nrens))
   allocate(meanS(nrobs))
   allocate(scaling(nrobs))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Observe ensemble to construct the matrix S=HA
   do iens =1, nrens
      do m =1, nrobs
         S(m,iens)      =  mem(obspos(m),iens)
      enddo
   enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Construct observation perturbations E
   if (E0(1,1) == 0.0 ) then
      print '(a)','       enkf: sampling measurement pert into E0'
      call obs_pert(E,nrens,nrobs,.true.,dx,rh,covmodel,obspos)
      E0=E
   else
      print '(a)','       enkf: Reusing previous E stored in E0'
      E=E0  ! use the same E for the different experiments
   endif

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
   if (mode_analysis==11 .or. mode_analysis==21) then
      print '(a,a)','       enkf: Exact R using covariance model: ',trim(covmodel)
      select case (trim(covmodel))
      case ('diagonal')
         do m=1,nrobs
            R(m,m)=obsvar
         enddo
      case ('gaussian')
         do i=1,nrobs
         do j=1,nrobs
            dist=abs(obspos(i)-obspos(j))
            if (dist > real(nx)/2.0) dist = real(nx) - dist
            R(i,j)=obsvar*exp(-dist**2/rd**2)
         enddo
         enddo
      case default
         print '(a,a)','       enkf: covmodel is invalid : ',trim(covmodel)
      end select
   else
      print '(a,a)','       enkf: lowrank R using covariance model: ',trim(covmodel)
      R=matmul(E,transpose(E))/float(nrens-1)
   endif


   allocate(innovation(nrobs))
   do m =1, nrobs
      innovation(m)      = obs(m)- meanS(m)
   enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Scaling of matrices
   print '(a)','       enkf: scale matrices'
   do m=1,nrobs
      scaling(m)=1./sqrt(R(m,m)) 
      S(m,:)=scaling(m)*S(m,:)
      E(m,:)=scaling(m)*E(m,:)
      D(m,:)=scaling(m)*D(m,:)
      innovation(m)=scaling(m)*innovation(m)
   enddo

   do j=1,nrobs
   do i=1,nrobs
      R(i,j)=scaling(i)*R(i,j)*scaling(j)
   enddo
   enddo


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Computing analysis
   call cpu_time(start)
   if (local==0) then
      print '(a,i2)','       enkf: calling global analysis with mode: ',mode_analysis
      call analysis(mem, R, E, S, D, innovation, nx, nrens, nrobs, .true., truncation, mode_analysis, &
                      lrandrot, lupdate_randrot, lsymsqrt, inflate, infmult)

   else 
      print '(a,i2)','       enkf: calling local analysis with mode: ',mode_analysis,local
      icall=0
      do i=1,nx
         nobs=0
         lobs(:)=.false.



! Distace based localization
         if (local == 1) then
            do m=1,nrobs
               if (real(abs(obspos(m)-i)) < robs) then
                  lobs(m)=.true.
                  nobs=nobs+1
               endif
            enddo
         endif

! Adaptive localization
         if(local == 2) then
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

            do m=1,nrobs
               if (abs(corr(m)) > obstreshold) then
                  lobs(m)=.true.
                  nobs=nobs+1
               endif
            enddo
         endif


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
                         lrandrot, local_rot, lsymsqrt, inflate, infmult)
            mem(i,:)=submem(1,:)
            deallocate(subD, subE, subS, subR, subinnovation, submem)
         endif
      enddo
   endif
   call cpu_time(finish)
   print '("   analysis: time = ",f6.3," seconds.")',finish-start
   print '(a)','       enkf: done'

   deallocate(innovation)
   deallocate(R)


end subroutine
end module
