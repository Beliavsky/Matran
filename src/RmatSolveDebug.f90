program RmatSolveTest

use MatranUtil_m
use RmatSolve_m
use RmatSum_m
use RmatProduct_m
use RmatNorm_m
use RmatPrint_m

implicit none

   real(wp) :: s
   integer :: tn
   integer :: i, j, n
   type(Rmat) :: A, B, C, X
   type(RmatLudpp) :: lu
   type(RmatChol) :: chl

   print *, ' '
   print *, "Enter test number"
   read *, tn
   print *, 'Test number ', tn


   select case(tn)

   case(1)

      s = 2.
      A = (/4,4/)
      A%a = 1.
      A = A/s
      call Print(A, 10, 2)

   case(2)

      A = (/4,4/)
      
      do i=1,4
         do j=1,4
            if (i>=j) then
               A%a(i,j) = i
            else
               A%a(i,j) = 0
            end if
         end do
      end do
      A%tag = 'LT'

      X = (/4,2/)
      X%a = 1

      B = A*X
      call SolveXiy(C, A, B)
      print *, normf(B - A*C)
      C = A.xiy.B
      print *, normf(B - A*C)

      B = A.xhy.X
      C = A.xihy.B
      print *, normf(B - (A.xhy.C))

      X = (/2,4/)
      X%a = 1

      B = X*A
      C = B.xyi.A
      print *, normf(B - (C*A))

      B = X.xyh.A
      C = B.xyih.A
      print *, normf(B - (C.xyh.A))

   case(3)

      A = (/4,4/)

      do i=1,4
         do j=1,4
            if (i<=j) then
               A%a(i,j) = i
            else
               A%a(i,j) = 0
            end if
         end do
      end do
      A%tag = 'UT'

      X = (/4,2/)
      X%a = 1

      B = A*X
      call SolveXiy(C, A, B)
      print *, normf(B - A*C)
      C = A.xiy.B
      print *, normf(B - A*C)

      B = A.xhy.X
      call SolveXihy(C, A, B)
      print *, normf(B - (A.xhy.C))
      C = A.xihy.B
      print *, normf(B - (A.xhy.C))

      X = (/2,4/)
      X%a = 1

      B = X*A
      C = B.xyi.A
      print *, normf(B - (C*A))
      
      B = X.xyh.A
      C = B.xyih.A
      print *, normf(B - (C.xyh.A))

   case(4)

      n=4
      A = (/n,n/)
      do i=1,n
         do j=1,n
            A%a(i,j) = exp(.1*(i+j))
         end do
      end do

      X = (/n,2/)
      X%a = 1

      B = A*X
      call SolveXiy(C, A, B, lu)
      print *, normf(B - A*C)/(normf(A)*normf(C))

      call SolveXiy(C, A, B, lu)
      print *, normf(B - A*C)/(normf(A)*normf(C))
      C = A.xiy.B
      print *, normf(B - A*C)/(normf(A)*normf(C))

      B = A.xhy.X
      C = A.xihy.B
      print *, normf(B - (A.xhy.C))/(normf(A)*normf(C))

      X = (/2,4/)
      X%a = 1

      B = X*A
      C = B.xyi.A
      print *, normf(B - (C*A))/(normf(A)*normf(C))

      B = X.xyh.A
      C = B.xyih.A
      print *, normf(B - (C.xyh.A))/(normf(A)*normf(C))

   case(5)

      n=100
      A = (/n,n/)
      do i=1,n
         do j=1,n
            A%a(i,j) = 1
         end do
         A%a(i,i) = n
      end do
      A = .xhx.A

      X = (/n,2/)
      X%a = 1

      B = A*X
      call SolveXiy(C, A, B)
      print *, normf(B - A*C)/(normf(A)*normf(C))

      call SolveXiy(C, A, B, chl=chl)
      print *, normf(B - A*C)/(normf(A)*normf(C))
      C = A.xiy.B
      print *, normf(B - A*C)/(normf(A)*normf(C))

      X = (/2,n/)
      X%a = 1

      B = X*A
      C = B.xyi.A
      print *, normf(B - C*A)/(normf(A)*normf(C))

   case default

      print *, 'No such test case.'

   end select
end program RmatSolveTest
