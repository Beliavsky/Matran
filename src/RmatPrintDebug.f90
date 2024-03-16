program RmatPrintTest

use MatranUtil_m
use Rmat_m
use RmatPrint_m

implicit none

   type(Rmat) :: A
   real(wp), pointer :: Array(:,:)
   integer i, j, m, n, tn

   print *, ' '
   print *, "Enter test number"
   read *, tn
   print *, 'Test number ', tn



   m = 5; n = 20
   allocate(Array(m,n))
   do j=1,n
      do i=1,m
         Array(i,j) = i+j
      end do
   end do

   select case(tn)

   case(1)


      call Print(Array(1:3,1:3), 3, 3, 9, 1)
      call Print(Array, m, n, 12, 4, 3)

   case(2)

      A = Array(1:4,1:5)
      call Print(A, 9, 1)

   case(3)
      A =  Array(1:4,1:5)
      call Print(A, 9, 1, "This is a test")

   case(4)
      A =  Array(1:4,1:5)
      call SetTemp(A)
      call Print(A, 9, 1, "Test temp variable")

   case default

      print *, "No such test case."

   end select


end program RmatPrintTest



