program RdiagProductDebug

use MatranUtil_m
use Rdiag_m
use RdiagProduct_m
use Rmat_m
use RmatPrint_m

implicit none

   integer :: tn
   type(Rdiag) :: D, E, F
   type(Rmat) :: A, B
   real(wp) :: s=2

   print *, ' '
   print *, "Enter test number"
   read *, tn
   print *, 'Test number ', tn

   select case (tn)

   case(1)

      D = (/4/); E=(/4,5/)
      D%a = (/1,2,3,4/)
      E%a = (/1,2,3,4/)
      F = D*E
      print *, F%a, F%order, F%adjustable, associated(F%temporary)
      F = s*E
      print *, F%a, F%order, F%adjustable, associated(F%temporary)
      F = E*s
      print *, F%a, F%order, F%adjustable, associated(F%temporary)

   case(2)
      B = (/7,9/)
      D = (/4/)
      D%a = (/1,2,3,4/)
      A = (/4,5/)
      A%a = 1
      B = D*A
      call Print(B, 9, 1)
      A = (/5,4/)
      A%a = 1
      B = A*D
      call Print(B, 9, 1)
     
   case(3)
      D = (/0/)
      E = (/0/)
      F = D*E
      print *, F%a, F%order, F%adjustable, associated(F%temporary)

   case default

      print *, 'No such test.'

   end select


end program RdiagProductDebug


