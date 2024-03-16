Program RmatJoinDebug

use MatranUtil_m
use RmatJoin_m
use RmatPrint_m
use RmatNorm_m
use RmatSum_m
use Rmat_m
implicit none

   type(Rmat) A, B, C

   integer :: tn, i, j

   print *, ' '
   print *, "Enter test number"
   read *, tn
   print *, 'Test number ', tn


   select case (tn)

   case(1) ! test JoinWE

      A = (/3, 4, 4, 5/)
      A%tag = 'UT'
      A = (/3, 4, 4, 5/)
      do i=1,3
         do j = 1, 4
            A%a(i,j) = i+j
         end do
      end do
      B = (/3, 2, 4, 5/)
      B%a(1:3, 1:2) = 2
      call JoinWE(C, A, B)
      call Print(C, 10, 2)
      print *, Norm1(C-(A.jwe.B))

   case(2) ! test JoinNS

      A = (/3, 4, 4, 5/)
      A%tag = 'LT'
      A%a(1:3, 1:4) = 1
      B = (/2, 4, 4, 5/)
      B%a(1:2, 1:4) = 2
      call JoinNS(C, A, B)
      call Print(C, 10, 2)
      print *, Norm1(C- (A.jns.B))

   case default

      print *, 'No such test.'

   end select
end program RmatJoinDebug
