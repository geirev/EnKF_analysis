module mod_states
! Modelstate definition for HYCOM
   use mod_dimensions
   real, parameter ::  onem=9806.
! standard model state
   type states
      real u(nx,ny,nz)                         ! 3-D u-velocity
      real v(nx,ny,nz)                         ! 3-D v-velocity
      real d(nx,ny,nz)                         ! 3-D Layer thickness
      real t(nx,ny,nz)                         ! 3-D Temperature
      real s(nx,ny,nz)                         ! 3-D Temperature 
      real ub(nx,ny)                           ! 2-D barotropic u-velocity
      real vb(nx,ny)                           ! 2-D barotropic v-velocity
      real pb(nx,ny)                           ! 2-D barotropic pressure
   end type states
   integer, parameter ::  global_ndim = 5*nx*ny*nz+3*nx*ny  ! Dimension of states

! single precision model state (used for read and write to files)
   type states4
      real*4 u(nx,ny,nz)
      real*4 v(nx,ny,nz)
      real*4 d(nx,ny,nz)
      real*4 t(nx,ny,nz)
      real*4 s(nx,ny,nz)
      real*4 ub(nx,ny)
      real*4 vb(nx,ny)
      real*4 pb(nx,ny)
   end type states4

! model state at one grid point (used in local analysis)
   type sub_states
      real u(nz)                               ! 1-D u-velocity
      real v(nz)                               ! 1-D v-velocity
      real d(nz)                               ! 1-D Layer thickness
      real t(nz)                               ! 1-D Temperature
      real s(nz)                               ! 1-D Temperature
      real ub                                  ! 0-D barotropic u-velocity
      real vb                                  ! 0-D barotropic v-velocity
      real pb                                  ! 0-D barotropic pressure
   end type sub_states
   integer, parameter ::  local_ndim=5*nz+3    ! Dimension of sub_states


! Overloaded and generic operators
   interface operator(+)
      module procedure add_states
   end interface

   interface operator(-)
      module procedure subtract_states
   end interface

   interface operator(*)
      module procedure states_real_mult,&
                       real_states_mult,&
                       states_states_mult
   end interface

!   interface operator(/)
!      module procedure divide_states
!   end interface

   interface assignment(=)
      module procedure assign_states
      module procedure states4to8
      module procedure states8to4
   end interface


contains
   type (sub_states) function getA(A,i,j,m)
      implicit none
      type(states), intent(in)     :: A
      integer, intent(in) :: i,j,m
      getA%u(:)=A%u(i,j,:)
      getA%v(:)=A%v(i,j,:)
      getA%d(:)=A%d(i,j,:)
      getA%t(:)=A%t(i,j,:)
      getA%s(:)=A%s(i,j,:)
      getA%ub=A%ub(i,j)
      getA%vb=A%vb(i,j)
      getA%pb=A%pb(i,j)
   end function getA

   subroutine putA(subA,A,i,j)
      implicit none
      type(sub_states), intent(in) :: subA
      type(states), intent(inout)  :: A
      integer, intent(in) :: i,j
      A%u(i,j,:)=subA%u(:)
      A%v(i,j,:)=subA%v(:)
      A%d(i,j,:)=subA%d(:)
      A%t(i,j,:)=subA%t(:)
      A%s(i,j,:)=subA%s(:)
      A%ub(i,j)=subA%ub
      A%vb(i,j)=subA%vb
      A%pb(i,j)=subA%pb
   end subroutine putA



   function add_states(A,B)
      type(states) add_states
      type(states), intent(in) :: A
      type(states), intent(in) :: B
       add_states%u = A%u + B%u
       add_states%v = A%v + B%v
       add_states%d = A%d + B%d
       add_states%t = A%t + B%t
       add_states%s = A%s + B%s
       add_states%ub = A%ub + B%ub
       add_states%vb = A%vb + B%vb
       add_states%pb = A%pb + B%pb
   end function add_states

   function subtract_states(A,B)
      type(states) subtract_states
      type(states), intent(in) :: A
      type(states), intent(in) :: B
       subtract_states%u = A%u - B%u
       subtract_states%v = A%v - B%v
       subtract_states%d = A%d - B%d
       subtract_states%t = A%t - B%t
       subtract_states%s = A%s - B%s
       subtract_states%ub = A%ub - B%ub
       subtract_states%vb = A%vb - B%vb
       subtract_states%pb = A%pb - B%pb
   end function subtract_states

   function states_real_mult(A,B)
      type(states) states_real_mult
      type(states), intent(in) :: A
      real, intent(in) :: B
       states_real_mult%u = B*A%u
       states_real_mult%v = B*A%v
       states_real_mult%d = B*A%d
       states_real_mult%t = B*A%t
       states_real_mult%s = B*A%s
       states_real_mult%ub = B*A%ub
       states_real_mult%vb = B*A%vb
       states_real_mult%pb = B*A%pb
   end function states_real_mult

   function real_states_mult(B,A)
      type(states) real_states_mult
      type(states), intent(in) :: A
      real, intent(in) :: B
       real_states_mult%u = B*A%u
       real_states_mult%v = B*A%v
       real_states_mult%d = B*A%d
       real_states_mult%t = B*A%t
       real_states_mult%s = B*A%s
       real_states_mult%ub = B*A%ub
       real_states_mult%vb = B*A%vb
       real_states_mult%pb = B*A%pb
   end function real_states_mult

   function states_states_mult(A,B)
      type(states) states_states_mult
      type(states), intent(in) :: A
      type(states), intent(in) :: B
       states_states_mult%u = A%u * B%u
       states_states_mult%v = A%v * B%v
       states_states_mult%d = A%d * B%d
       states_states_mult%t = A%t * B%t
       states_states_mult%s = A%s * B%s
       states_states_mult%ub = A%ub * B%ub
       states_states_mult%vb = A%vb * B%vb
       states_states_mult%pb = A%pb * B%pb
   end function states_states_mult


   subroutine assign_states(A,r)
      type(states), intent(out) :: A
      real, intent(in) :: r
       A%u = r
       A%v = r
       A%d = r
       A%t = r
       A%s = r
       A%ub = r
       A%vb = r
       A%pb = r
   end subroutine assign_states

   subroutine states4to8(A,B)
      type(states), intent(out) :: A
      type(states4), intent(in)  :: B
      A%u=DBLE(B%u)
      A%v=DBLE(B%v)
      A%d=DBLE(B%d)
      A%t=DBLE(B%t)
      A%s=DBLE(B%s)
      A%ub=DBLE(B%ub)
      A%vb=DBLE(B%vb)
      A%pb=DBLE(B%pb)
   end subroutine states4to8

   subroutine states8to4(A,B)
      type(states), intent(in)  :: B
      type(states4),  intent(out) :: A
      A%u=real(B%u)
      A%v=real(B%v)
      A%d=real(B%d)
      A%t=real(B%t)
      A%s=real(B%s)
      A%ub=real(B%ub)
      A%vb=real(B%vb)
      A%pb=real(B%pb)
   end subroutine states8to4

end module mod_states

