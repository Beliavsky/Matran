program RmatSumDebug

use MatranUtil_m
use RmatSum_m
use RmatPrint_m

implicit none

   integer :: tn, i, j
   real(wp), pointer :: arya(:,:), aryb(:,:)
   type(Rmat) :: A, B, C

   print *, ' '
   print *, "Enter test number"
   read *, tn
   print *, 'Test number ', tn

   select case (tn)

   case(1)
      allocate(arya(2,3), aryb(2,3))
      do j=1,3
         do i=1,2
            arya(i,j) = i+j
            aryb(i,j) = i-j
         end do
      end do
      A = arya
      call Print(A, 10, 2)
      B = aryb
      call Print(B, 10, 2)
      call Plus(C, A, B)
      call Print(C, 10, 2)

      
      call Minus(C, A, B)
      call Print(C, 10, 2)

      call Minus(C, A)
      call Print(C, 10, 2)

   case(2)
      allocate(arya(2,3), aryb(2,3))
      do j=1,3
         do i=1,2
            arya(i,j) = i+j
            aryb(i,j) = i-j
         end do
      end do
      A = arya
      B = aryb
      C = A + B
      call Print(C, 10, 2)

      C = C - B
      call Print(C, 10, 2)

      C = -A
      call Print(C, 10, 2)

   case(3)
      allocate(arya(2,3), aryb(2,3))
      do j=1,3
         do i=1,2
            arya(i,j) = i+j
            aryb(i,j) = i-j
         end do
      end do
      A = arya
      B = aryb
      C = A + (B - A)
      call Print(C, 10, 2)
      C = C - B
      call Print(C, 10, 2)


   case default
      print *, 'No such test case.'

   end select



end program RmatSumDebug
