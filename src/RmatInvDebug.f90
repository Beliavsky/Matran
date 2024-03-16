program RmatInvDebug

use MatranUtil_m
use Rmat_m
use RmatSum_m
use RmatProduct_m
use RmatInv_m
use RmatNorm_m
use RmatEye_m
use RmatLudpp_m
use RmatChol_m

implicit none

   type(Rmat) :: A, C
   integer :: tn, i, j, n
   type(RmatLudpp):: luda
   type(RmatChol):: chola


   n = 20
   A = (/n , n/)
   call Random_Number(A%a)
   A%a = A%a - 0.5

   print *, ' '
   print *, "Enter test number"
   read *, tn
   print *, 'Test number ', tn


   select case (tn)
   case (1) 
      print *, "A%tag = GE"
      call Inv(C, A)
      print *, norm(A*C-Reye(n))/norm(A)

   case (2)
      print *, "A%tag = HE"
      do i = 1, n
         do j = i+1, n
            A%a(i,  j) = A%a(j, i)
         end do
      end do
      A%tag = "HE"
      call Inv(C, A)
      print *, norm(A*C-Reye(n))/norm(A)

   case (3)
      print *, "A%tag = UT"
      do i = 1, n
         do j = i+1, n
            A%a(j, i) = 0
         end do
      A%a(i,i) = 1
      end do

      A%tag = "UT"
      call Inv(C, A)
      print *, norm(A*C-Reye(n))/norm(A)

   case (4)
      print *, "A%tag = LT"
      do i = 1, n
         do j = i+1, n
            A%a(i, j) = 0
         end do
      A%a(i,i) = 1
      end do

      A%tag = "LT"
      call Inv(C, A)
      print *, norm(A*C-Reye(n))/norm(A)

   case (5)
      print *, "A%tag = HP"
      A = A .xhy. A
      A%tag = "HP"
      call Inv(C, A)
      print *, norm(A*C-Reye(n))/norm(A)

   case (6)
      print *, "A%tag = GE, using Ludpp"

      call Ludpp(luda, A)
      call Inv(C, A, luda)
      print *, norm(A*C-Reye(n))/norm(A)

!      call Clean(luda)
      call Inv(C, A, luda)
      print *, norm(A*C-Reye(n))/norm(A)

   case (7)
      print *, "A%tag = HP, using Chol"
      A = A .xhy. A
      A%tag = "HP"

      call Chol(chola, A)
      call Inv(C, A, chola = chola)
      print *, norm(A*C-Reye(n))/norm(A)

!      call Clean(chola)
      call Inv(C, A, chola = chola)
      print *, norm(A*C-Reye(n))/norm(A)
                    
      C = .inv.(A)
      print *, norm(A*C-Reye(n))/norm(A)

   case default

      print *, 'No such test.'

   end select

end program RmatInvDebug
