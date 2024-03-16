program RmatTest

use MatranUtil_m
use Rmat_m
use RmatPrint_m

implicit none

   integer :: tn, i, j
   real(wp) :: s = 5 
   real(wp), pointer :: ary(:,:)
   type(Rmat) :: A, B

   print *, ' '
   print *, "Enter test number"
   read *, tn
   print *, 'Test number ', tn


   select case (tn)

   case(1)
      allocate(ary(2,3))
      do j=1,3
         do i=1,2
            ary(i,j) = i+j
         end do
      end do
      B = ary
      call Print(B, 10, 2)
      call Print(.rm.ary, 10,2)
      A = B
      call Print(A, 10, 2)
      B = ary(1:2,1:2)
      call Print(B, 10, 2)
      A = B
      call Print(A, 10, 2)
      call Print(A%a, A%narow, A%nacol, 10, 2)

      call Print(.rm.s, 10, 2)

   case (2)

      A = (/2,3/)
      call Print(A, 10, 2)
      call Clean(A)
      call Print(.rm.(/4,3/), 9, 1)

      A = (/0,0,3,4/)
      print *, associated(A%a)
      call Print(A, 10, 2)
      A = (/2,3/)
      call Print(A, 10, 2)
      call Clean(A)
      print *, associated(A%temporary)
      call Print(A, 10, 2)

      call Print(.rm.(/4,3/), 9, 1)

   case(3)
      A = (/2,3/)
      B = (/4,5/)
      A = B
      call Print(A, 10, 2)

   case(4)
      call farb(.rm.(/2,3/))

   case default
      print *, 'No such test case.'

   end select      

contains

   subroutine farb(A)
      type(Rmat), intent(in) :: A;
      call GuardTemp(A)
      call forb(A)
      call Print(A, 10, 2)
      call CleanTemp(A)
   end subroutine farb

   subroutine forb(A)
      type(Rmat), intent(in) :: A;
      call GuardTemp(A)
      call Print(A, 10, 2)
      call CleanTemp(A)
   end subroutine forb

end program RmatTest



