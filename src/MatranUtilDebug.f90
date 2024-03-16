program MatranUtilDebug

use MatranUtil_m

implicit none
   integer :: tn
   integer, pointer :: I1(:)=>null(), I2(:,:)=>null()
   real(wp), pointer :: R1(:)=>null(), R2(:,:)=>null()
   complex(wp), pointer :: C1(:)=>null(), C2(:,:)=>null()

   Print *, 'Enter test number'
   read *, tn
   print *, 'Test number ', tn

   select case(tn)

   case(1)

      call MatranError("This is an error message.")

   case(2)

      call SupportError("Error in xxxxxx: info =", 69)

   case(3)

      call ReshapeAry(I1, 5)
      print *, shape(I1)
      call ReshapeAry(I1, 3)
      print *, shape(I1)
      call ReshapeAry(I1, 7)
      print *, shape(I1)

      call ReshapeAry(I2, 3, 5)
      print *, shape(I2)
      call ReshapeAry(I2, 2, 2)
      print *, shape(I2)
      call ReshapeAry(I2, 5, 2)
      print *, shape(I2)
      call ReshapeAry(I2, 2, 5)
      print *, shape(I2)

   case(4)

      call ReshapeAry(R1, 5)
      print *, shape(R1)
      call ReshapeAry(R1, 3)
      print *, shape(R1)
      call ReshapeAry(R1, 7)
      print *, shape(R1)

      call ReshapeAry(R2, 3, 5)
      print *, shape(R2)
      call ReshapeAry(R2, 2, 2)
      print *, shape(R2)
      call ReshapeAry(R2, 5, 2)
      print *, shape(R2)
      call ReshapeAry(R2, 2, 5)
      print *, shape(R2)

   case(5)

      call ReshapeAry(C1, 5)
      print *, shape(C1)
      call ReshapeAry(C1, 3)
      print *, shape(C1)
      call ReshapeAry(C1, 7)
      print *, shape(C1)

      call ReshapeAry(C2, 3, 5)
      print *, shape(C2)
      call ReshapeAry(C2, 2, 2)
      print *, shape(C2)
      call ReshapeAry(C2, 5, 2)
      print *, shape(C2)
      call ReshapeAry(C2, 2, 5)
      print *, shape(C2)

   case default

      print *, "No such test case."

   end select

end program MatranUtilDebug
