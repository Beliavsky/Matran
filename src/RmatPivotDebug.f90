program RmatPivotDebug

use MatranUtil_m
use Rmat_m
use RmatPivot_m
use RmatPrint_m

implicit none

   type(Rmat) :: A
   integer :: pvt(4), i, j

   A = (/4,2/)
   do i=1,4
      A%a(i,1:2) = i
      pvt(i) = i+1
   end do
   pvt(4) = 4

   call PivotRow(A, pvt, 4)
   call Print(A, 9, 1)
   call PivotInvRow(A, pvt, 4)
   call Print(A, 9, 1)

   A = (/2,4/)
   do j=1,4
      A%a(1:2,j) = j
      pvt(j) = j+1
   end do
   pvt(4) = 4

   call PivotCol(A, pvt, 4)
   call Print(A, 9, 1)
   call PivotInvCol(A, pvt, 4)
   call Print(A, 9, 1)

end program RmatPivotDebug
