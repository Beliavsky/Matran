program RdiagDebug

use MatranUtil_m
use Rdiag_m
use Rmat_m
use RmatPrint_m

implicit none

   integer tn, i
   real(wp) :: s = 4
   real(wp), pointer :: ary(:)=>null()
   type(Rdiag) :: D, E
   type(Rmat) :: A

   print *, "Enter test number"
   read *, tn
   print *, 'Test case ', tn

   select case (tn)

   case(1)
      allocate(ary(4))
      do i=1,4
         ary(i) = i
      end do
      D = ary(1:3)
      print *, D%a, D%order, D%na, D%adjustable, associated(D%temporary)
      E = (/5/)
      print *, E%a, E%order, E%na, E%adjustable, associated(E%temporary)
      E = (/5,6/)
      print *, E%a, E%order, E%na, E%adjustable, associated(E%temporary)
      E = .rd. s
      print *, E%a, E%order, E%na, E%adjustable, associated(E%temporary)
      E = D
      print *, E%a, E%order, E%na, E%adjustable, associated(E%temporary)
      D = .rd.ary(1:3)
      print *, D%a, D%order, D%na, D%adjustable, associated(D%temporary)

      D = s
      print *, D%a, D%order, D%na, D%adjustable, associated(D%temporary)
      D = .rd.s
      print *, D%a, D%order, D%na, D%adjustable, associated(D%temporary)
      
      call Clean(E)
      print *, associated(E%a), E%order, E%na, E%adjustable, associated(E%temporary)
      E = D
      call SetTemp(E)
      call CleanTemp(E)
      print *, associated(E%a), E%order, E%na, E%adjustable, associated(E%temporary)

   case(2)

      allocate(ary(4))
      do i=1,4
         ary(i) = i
      end do
      D = ary(1:3)

      A = (/4,5/)
      A = D
      call Print(A, 9, 1)

      call Print(.rm.D, 9, 1)

  case(3)
      D = (/2,3/)
      call SetTemp(D)
      call farb(D)
      print *, associated(D%a), D%order, D%na, D%adjustable, D%temporary
   end select

contains

   subroutine farb(D)
      type(Rdiag), intent(in) :: D;
      call GuardTemp(D)
      call forb(D)
      print *, associated(D%a), D%order, D%na, D%adjustable, D%temporary
      call CleanTemp(D)
   end subroutine farb

   subroutine forb(D)
      type(Rdiag), intent(in) :: D;
      call GuardTemp(D)
      print *, associated(D%a), D%order, D%na, D%adjustable, D%temporary
      call CleanTemp(D)
   end subroutine forb


end program RdiagDebug

