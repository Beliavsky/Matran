program RmatLudppTest

use MatranUtil_m
use RmatProduct_m
use RmatSum_m
use RmatNorm_m
use RmatPivot_m
use RmatLudpp_m
use RmatPrint_m
implicit none

   type(Rmat) :: A
   type(RmatLudpp) :: lu
   
   integer i, info, j, tn, m, n, minmn

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

      call ludpp(lu, A)
      print *, lu%pvt
      call Print(lu%L, 10, 2)
      call Print(lu%U, 10, 2)
      print *, lu%companion
      print *, Normf(A - lu%L*lu%U)/normf(A)

      call Clean(lu)

      print *, associated(lu%pvt)
      call Print(lu%L, 10, 2)
      call Print(lu%U, 10, 2)
      print *, lu%companion

   case(2)

      A = (/4,4,5,5/)
      do j=1,4
         do i=j,4
            A%a(i,j) = j
            A%a(j,i) = A%a(i,j)
         end do
      end do
      A%a(4,4) = 3

      call ludpp(lu, A, info=info)
      print *, info
      print *, lu%pvt
      call Print(lu%L, 10, 2)
      call Print(lu%U, 10, 2)
      print *, lu%companion
      print *, Normf(A - lu%L*lu%U)/normf(A)


   case(3)
      print *, "Enter m and n"
      read *, m, n
      print *, 'm = ', m, 'n = ', n
      A = (/m,n/)
      minmn = min(m, n)
      do i=1,m
         do j=1,n
            A%a(i,j) = 1
         end do
         if (i<=minmn) &
            A%a(i,minmn-i+1) = A%a(i,minmn-i+1) + 2*max(m, n)
      end do

      call ludpp(lu, A)
      print *, lu%pvt
      call PivotRow(A, lu%pvt, lu%npvt)
      print *, Normf(A - lu%L*lu%U)/Normf(A)

  case default

      print *, 'No such test case.'

  end select



end program RmatLudppTest
