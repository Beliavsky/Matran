Program RmatSubmatrixDebug

use MatranUtil_m
use Rmat_m
use RmatPrint_m
use RmatSubmatrix_m
implicit none

   type(Rmat) A, C

   integer :: tn, i


   A = (/3, 4, 4, 5/)
   A%a(1:3,1:4) = reshape((/(I, I=1,12)/), (/3,4/))
   call Print(A, 9, 1)

   print *, ' '
   print *, "Enter test number"
   read *, tn
   print *, 'Test number ', tn


   select case (tn)

   case(1)

      call Print(Sbm(A,2,3,1,3), 9, 1)

   case(2)

      call Print(Col(A,2,3), 9, 1)
      call GetCol(C, A, 2, 3)
      call Print(C, 9, 1)

      call Print(Col(A,2), 9, 1)
      call GetCol(C, A, 2)
      call Print(C, 9, 1)


   case(3)

      call Print(Row(A,2,3), 9, 1)
      call GetRow(C, A, 2, 3)
      call Print(C, 9, 1)

      call Print(Row(A,2), 9, 1)
      call GetRow(C, A, 2)
      call Print(C, 9, 1)

   end select
end program RmatSubmatrixDebug
