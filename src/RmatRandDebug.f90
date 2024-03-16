program RmatRandDebug

use MatranUtil_m
use RmatRand_m
use RmatPrint_m

implicit none

   type(Rmat) :: A
   
   integer :: tn

   A = (/1000, 1/)

   print *, ' '
   print *, "Enter test number"
   read *, tn
   print *, 'Test number ', tn

   select case (tn)

   case(1) 

      call Clean(A)

      call RandU(A, 5)

      call Print(A, 9, 1)

      call Print(RrandU(5), 9, 1)

      call Print(RrandU(5, 3), 9, 1)

      call Print(RrandU(3, 5), 9, 1)


   case (2)

      call Clean(A)

      call RandN(A, 5)

      call Print(A, 9, 1)

      call Print(RrandN(5), 9, 1)

      call Print(RrandN(5, 3), 9, 1)

      call Print(RrandN(3, 5), 9, 1)

   case default

      print *, 'No such case.'

   end select

end program RmatRandDebug
