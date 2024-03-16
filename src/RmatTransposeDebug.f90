program RmatTransposeDebug

use MatranUtil_m
use Rmat_m
use RmatTranspose_m
use RmatPrint_m

implicit none

   integer :: i, j
   type(Rmat) :: A

   A = (/3, 4, 4, 5/)
   do i=1,3
      do j=1,4
         A%a(i,j) = i+j**2
      end do
   end do
   call Print(A, 10, 2)
   call Print(.ctp.A, 10, 2)

   A = (/5, 0, 7, 3/)
   call Print(A, 10, 2)
   call Print(.ctp.A, 10, 2)

   A = (/0, 3, 7, 3/)
   call Print(A, 10, 2)
   call Print(.ctp.A, 10, 2)

end program RmatTransposeDebug
