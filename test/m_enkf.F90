module m_enkf
contains
integer function printlocal(lobs,nrobs,nobs)
   integer, intent(in) :: nrobs
   integer, intent(in) :: nobs
   logical, intent(in) :: lobs(nrobs)
   integer m
   write(*,'(a,i4)',advance='no')'localization:',nobs
   do m=1,nrobs
      if (lobs(m)) write(*,'(tr1,i4)',advance='no')m
   enddo
   write(*,*)
   printlocal=1
end function

real function measerrinf(cutoff,corr,inf,mininf)
   real, intent(in) :: cutoff
   real, intent(in) :: corr
   real, intent(in) :: mininf
   integer, intent(in) :: inf
   real, parameter :: pi=3.1415927
   real a,b,c,fx
   c=sqrt(sqrt(cutoff))
   a=-c/(1.0-c)
   b=1.0/(1.0-c)
   if (abs(corr) .le. cutoff) then
      fx=0.0
      measerrinf=0.0
   else
      select case (inf)
      case(0)
         measerrinf=1.0
      case(1)
         c=cutoff
         a=-c/(1.0-c)
         b=1.0/(1.0-c)
         fx=0.5-0.5*cos((a+b*abs(corr))*pi)
         measerrinf=1.0/(fx+0.001)
      case(2)
         c=sqrt(cutoff)
         a=-c/(1.0-c)
         b=1.0/(1.0-c)
         fx=0.5-0.5*cos((a+b*sqrt(abs(corr)))*pi)
         measerrinf=1.0/(fx+0.001)
      case(3)
         c=sqrt(sqrt(cutoff))
         a=-c/(1.0-c)
         b=1.0/(1.0-c)
         fx=0.5-0.5*cos((a+b*sqrt(sqrt(abs(corr))))*pi)
         measerrinf=1.0/(fx+0.001)
      case default
         print *,'measerrinf: not define inf=',inf
      end select
      measerrinf=min(mininf,measerrinf)
   endif
end function

subroutine enkf(mem,nx,nrens,obs,obsvar,obspos,nrobs,mode_analysis,truncation,covmodel,dx,rh,rd,&
               &lrandrot,lsymsqrt,inflate,infmult,local,robs,obstreshold,E,ne)

   use m_getD
   use m_obspert
   use m_preexact
   implicit none
   integer, intent(in) :: nx
   integer, intent(in) :: nrens
   integer, intent(in) :: ne
   integer, intent(in) :: nrobs
   integer, intent(inout) :: mode_analysis

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

   real,    intent(inout) :: E(nrobs,nrens*ne)
   integer iprt

!   character(len=2) cmode
!   real, allocatable :: Rsamp(:,:)
   real, allocatable :: R(:,:)
   real, allocatable :: D(:,:)
   real, allocatable :: S(:,:)
   real, allocatable :: meanS(:)
   real, allocatable :: innovation(:)
   real, allocatable :: scaling(:)
   real, allocatable :: ES(:,:)

   integer iens,m,i,j
   logical :: lupdate_randrot=.true.

! Local analysis variables
   integer l,icall
   logical lobs(nrobs)           ! which measurements are active
   integer nobs
   logical local_rot
   real corr(nrobs),stdA,aveA,stdS
!   real dist
   real, allocatable, dimension(:,:) :: subS,subE,subD, subR
   real, allocatable, dimension(:,:)   :: submem
   real, allocatable :: subinnovation(:)
   real, allocatable :: errinf(:)
   real :: start, finish
! End Local analysis variables

   allocate(D(nrobs,nrens))
   allocate(S(nrobs,nrens))
   allocate(meanS(nrobs))
   allocate(scaling(nrobs))
   allocate(errinf(nrobs))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Observe ensemble to construct the matrix S=HA
   do iens =1, nrens
      do m =1, nrobs
         S(m,iens)      =  mem(obspos(m),iens)
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
!   allocate(Rsamp(nrobs,nrobs))
   R=0.0
   if (mode_analysis==11 .or. mode_analysis==21) then
      print '(a,a)','       enkf: Exact R using covariance model: ',trim(covmodel)
      select case (trim(covmodel))
      case ('diagonal')
         do m=1,nrobs
            R(m,m)=obsvar
         enddo
      case ('gaussian')
!         do i=1,nrobs
!         do j=1,nrobs
!            dist=abs(obspos(i)-obspos(j))
!            if (dist > real(nx)/2.0) dist = real(nx) - dist
!            R(i,j)=obsvar*exp(-dist**2/rd**2)
!         enddo
!         enddo
         allocate(ES(nrobs,10000))
         print '(a)','main: sampling measurement pert into ES(nrobs,2000)'
         call obspert(ES,10000,nrobs,.true.,dx,rd,covmodel,obspos)
         ES(:,:)=sqrt(obsvar)*ES(:,:)
         R=matmul(ES,transpose(ES))/real(10000-1)
         deallocate(ES)
      case default
         print '(a,a)','       enkf: covmodel is invalid : ',trim(covmodel)
      end select

   else
      print '(a,a)',   '       enkf: lowrank R using covariance model: ',trim(covmodel)
      print '(a,i6,a)','       enkf: lowrank R generated using ',ne*nrens,' realizations'
      R=matmul(E,transpose(E))/float(nrens*ne-1)
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


   if (mode_analysis==0) then
      call preexact()
      mode_analysis=10
   endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Computing analysis
   call cpu_time(start)
   if (local==0) then
      print '(a,i2)','       enkf: calling global analysis with mode: ',mode_analysis
      call analysis(mem, R, E, S, D, innovation, nx, nrens, nrobs, .true., truncation, mode_analysis, &
                      lrandrot, lupdate_randrot, lsymsqrt, inflate, infmult, ne)

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
                  errinf(nobs)=measerrinf(obstreshold,corr(m),0,10.0)
            !      print *,'errinf:',m,nobs,corr(m),errinf(nobs)
               endif
            enddo
            iprt=printlocal(lobs,nrobs,nobs)
         endif


         if (nobs > 0) then
            allocate(subD(nobs,nrens))
            allocate(subE(nobs,nrens*ne))
            allocate(subS(nobs,nrens))
            allocate(subR(nobs,nobs))
            call getD(D,subD,nrobs,nrens,lobs,nobs) ! the innovations to use
            call getD(E,subE,nrobs,nrens*ne,lobs,nobs) ! the observation errors to use
            subR=matmul(subE,transpose(subE))/float(nrens*ne)
            do m=1,nobs
               subE(m,:)=errinf(m)*subE(m,:)
            enddo
            call getD(S,subS,nrobs,nrens,lobs,nobs) ! the HA' to use
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
                         lrandrot, local_rot, lsymsqrt, inflate, infmult, ne)
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
