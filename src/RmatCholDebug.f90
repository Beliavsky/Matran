program RmatCholDebug

use MatranUtil_m
use Rmat_m
use RmatProduct_m
use RmatSum_m
use RmatNorm_m
use RmatChol_m
use RmatPrint_m
implicit none

   type(Rmat) :: A
   type(RmatChol), target :: chl
   type(Rmat), pointer :: R

   integer i, info, j, n, tn

   print *, ' '
   print *, "Enter test number"
   read *, tn
   print *, 'Test number ', tn


   select case (tn)

   case(1)
      A = (/4,4/)
      do j=1,4
         do i=j,4
            A%a(i,j) = j
            A%a(j,i) = A%a(i,j)
         end do
      end do

      A%tag = 'HP'

      call Chol(chl, A); R => chl%R
      call Print(R, 10, 2)
      print *, Normf(A - .xhx.R)/Normf(A)

   case(2)
      A = (/4,4/)
      do j=1,4
         do i=j,4
            A%a(i,j) = j
            A%a(j,i) = A%a(i,j)
         end do
      end do

      A%a(3,3) = -1;

      A%tag = 'HP'

      call Chol(chl, A, info=info); R => chl%R
      print *, info
      call Print(R, 10, 2)
      print *, Normf(A - .xhx.R)/Normf(A)

   case(3)
      print *, "Enter n"
      read *, n
      print *, 'n = ', n      
      A = (/n,n, n+1,n+1/)
      do i=1,n
         do j=1,n
            A%a(i,j) = 1
         end do
         A%a(i,i) = n+2
      end do

      A%tag = 'HP'

      Call Chol(chl, A)
      R => chl%R
      print *, Normf(A - .xhx.R)/Normf(A)

   case default

      print *, 'No such test case.'

   end select

end program RmatCholDebug
