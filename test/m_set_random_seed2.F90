module m_set_random_seed2
contains
subroutine set_random_seed2
! Sets a random seed based on the system and wall clock time
   implicit none 

   integer , dimension(8)::val
   integer cnt
   integer sze
   integer, allocatable, dimension(:):: pt

   call DATE_AND_TIME(values=val)
   call SYSTEM_CLOCK(count=cnt)
   call RANDOM_SEED(size=sze)
   allocate(pt(max(sze,2)))
   pt(1) = val(8)*val(3)
   pt(2) = cnt
   call RANDOM_SEED(put=pt)
   deallocate(pt)
end subroutine set_random_seed2
end module m_set_random_seed2
