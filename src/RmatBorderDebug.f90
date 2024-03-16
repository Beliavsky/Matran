Program RmatBorderDebug

use MatranUtil_m
use RmatBorder_m
use RmatJoin_m
use RmatPrint_m
use RmatNorm_m
use RmatSum_m
use RmatSubmatrix_m
use Rmat_m
implicit none

   type(Rmat) A, B, C

   integer :: tn, i, j

   print *, ' '
   print *, "Enter test number"
   read *, tn
   print *, 'Test number ', tn


   select case (tn)

   case(1) ! test borderE

      A = (/3, 4, 4, 5/)
      A%tag = 'UT'
      A = (/3, 4, 4, 5/)
      A%a(1:3, 1:4) = 1
      B = (/3, 2, 4, 5/)
      B%a(1:3, 1:2) = 2
      C = A.jwe.B
      call BorderE(A, B)
      print *, Norm1(C-A)

   case(2) ! test borderN

      A = (/3, 4, 4, 5/)
      A%tag = 'LT'
      do i=1,3
         do j = 1, 4
            A%a(i,j) = i+j
         end do
      end do
      B = (/2, 4, 4, 5/)
      B%a(1:2, 1:4) = 2
      C = B.jns.A
      call BorderN(A, B)
      print *, Norm1(C-A)

   case(3) ! test BorderW

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
      C = B.jwe.A
      call BorderW(A, B)
      print *, Norm1(C-A)

   case(4) ! test Border S

      A = (/3, 4, 4, 5/)
      A%tag = 'LT'
      A%a(1:3, 1:4) = 1
      B = (/2, 4, 4, 5/)
      B%a(1:2, 1:4) = 2
      C = A.jns.B
      call BorderS(A, B)
      print *, Norm1(C-A)

   case(5) ! test border SE
      A = (/3, 4, 4, 5/)
      A%a = reshape((/(I, I=1,12)/), (/3,4/))
      
      B = Sbm(A,1,2,1,2)
      call BorderSE(B, Sbm(A,3,3,1,2), Sbm(A,1,2,3,4), Sbm(A,3,3,3,4))
      print *, Norm1(A-B)

   case(6) ! test border SW
      A = (/3, 4, 4, 5/)
      A%a = reshape((/(I, I=1,12)/), (/3,4/))
      
      B = Sbm(A,1,2,3,4)
      call BorderSW(B, Sbm(A,3,3,3,4),Sbm(A,1,2,1,2),Sbm(A,3,3,1,2))
      print *, Norm1(A-B)


   case(7) ! test border NE
      A = (/3, 4, 4, 5/)
      A%a = reshape((/(I, I=1,12)/), (/3,4/))
      
      B = Sbm(A,3,3,1,2)
      call BorderNE(B, Sbm(A,1,2,1,2), Sbm(A,3,3,3,4), Sbm(A,1,2,3,4))
      print *, Norm1(A-B)

   case(8) ! test border NW
      A = (/3, 4, 4, 5/)
      A%a = reshape((/(I, I=1,12)/), (/3,4/))
      
      B = Sbm(A,3,3,3,4)
      call BorderNW(B, Sbm(A,1,2,3,4), Sbm(A,3,3,1,2),  Sbm(A,1,2,1,2))
      print *, Norm1(A-B)

   case default
      print *, 'No such test.'

   end select
end program RmatBorderDebug

