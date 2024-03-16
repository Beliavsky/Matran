program RmatEyeDebug

use MatranUtil_m
use RmatEye_m
use RmatPrint_m

implicit none

   type(Rmat) :: A

   call Clean(A)

   call Eye(A, 5)

   call Print(A, 9, 1)

   call Print(Reye(5), 9, 1)


   call Print(Reye(5, 3), 9, 1)


   call Print(Reye(3, 5), 9, 1)



end program RmatEyeDebug
