program RdiagSolveDebug

use MatranUtil_m
use Rdiag_m
use RdiagSolve_m
use Rmat_m
use RmatPrint_m

implicit none

   integer :: tn
   type(Rdiag) :: D, E, F
   type(Rmat) :: A
   real(wp) :: s=4

   print *, ' '
   print *, "Enter test number"
   read *, tn
   print *, 'Test number ', tn


   select case (tn)

   case(1)

      D = (/4/)
      D%a = (/1,2,3,4/)

      E = D/s

      print *, E%a, E%order, E%adjustable, associated(E%temporary)


    case(2)


      D = (/4/)
      D%a = (/1,2,3,4/)

      E = (/4/)
      E%a = (/4,3,2,1/)

      F = D.xiy.E

      print *, F%a, F%order, F%adjustable, associated(F%temporary)

      F = D.xyi.E

      print *, F%a, F%order, F%adjustable, associated(F%temporary)


   case(3)

      D = (/4/)
      D%a = (/1,2,3,4/)

      E = (/2/)
      E%a = (/1,2/)

      A = (/4,2/)
      A%a = 1



      call Print(D.xiy.A, 10, 2)
      call Print(A.xyi.E, 10, 2)
      

   case default

      print *, 'No such test.'

   end select



end program RdiagSolveDebug

