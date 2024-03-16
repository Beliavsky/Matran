program RmatNorm2Debug

use MatranUtil_m
use RmatNorm2_m
use Rmat_m
use Rdiag_m
use RmatQR_m
use RmatTranspose_m
use RmatProduct_m
use RdiagProduct_m

implicit none

   type(Rmat) :: A
   type(Rdiag) :: D
   type(RmatQR) :: qr1, qr2

   integer :: i, m, mm, n, nn

   print *, 'Enter m and n'

   read *, mm, nn
   m = max(mm, nn)
   n = min(mm, nn)
   print *, 'm = ', m, 'n = ', n

   A = (/m,n/)
   call random_number(A%a)
   call QR(qr1, A)
   A = (/n,n/)
   call random_number(A%a)
   call QR(qr2, A)

   D = (/n/)
   do i=1,n
      D%a(i) = i
   end do

   A = qr1%Q * D * qr2%Q

   print *, Norm2(A)
   print *, Norm2(.ctp.A)

end program RmatNorm2Debug
