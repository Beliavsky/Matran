program RmatNormDebug

use MatranUtil_m
use Rmat_m
use RmatNorm_m
implicit none

   type(Rmat) A

   A = (/4,5,5,6/)

   print *, normf(A)

   A%a(1:4,1:5) = 1

   print *, normf(A), sqrt(20.0_wp)
   print *, norm1(A)
   print *, norminf(A)

   A = (/0, 5, 2, 6/)
   A%a(1:0, 1:5) = 1

   print *, normf(A), norm1(A), norminf(A)

end program RmatNormDebug

