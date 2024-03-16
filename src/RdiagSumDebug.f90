program RdiagSumDebug

use MatranUtil_m
use Rdiag_m
use RdiagSum_m
use Rmat_m
use RmatPrint_m

implicit none

   integer :: tn
   type(Rdiag) :: D, E, F
   type(Rmat) :: A, B

   print *, ' '
   print *, "Enter test number"
   read *, tn
   print *, 'Test number ', tn


   select case (tn)

   case(1)

      D = (/4, 5/); E=(/4/)
      D%a = (/1,2,3,4/)
      E%a = (/1,2,3,4/)
      F = D + E
      print *, F%a, F%order, F%adjustable, associated(F%temporary)
      F = D - E
      print *, F%a, F%order, F%adjustable, associated(F%temporary)
      F = -D
      print *, F%a, F%order, F%adjustable, associated(F%temporary)

   case(2)
      B = (/7,9/)
      D = (/4/)
      D%a = (/1,2,3,4/)
      A = (/4,4/)
      A%a = 1
      B = D + A
      call print(B, 10, 1)
      B = A + D
      call print(B, 10, 1)
      B = D - A
      call print(B, 10, 1)
      B = A - D
      call print(B, 10, 1)

   case default

      print *, 'No such test.'

   end select

end program RdiagSumDebug

