program RmatProductDebug
use Rmat_m
use RmatPrint_m
use MatranUtil_m
use RmatProduct_m

implicit none

   integer :: tn, i, j
   type(Rmat) :: A, B, C
   real(wp) :: s=2

   print *, ' '
   print *, "Enter test number"
   read *, tn
   print *, 'Test number ', tn


   select case (tn)

   case(1)
      A = (/2,3,3,4/)
      A%a(1:2,1:3) = 1
      call Print(A, 10, 2)
      call Times(C, s, A)
      call Print(C, 10, 2)
      C = s*A
      call Print(C, 10, 2)
      C = A*s
      call Print(C, 10, 2)


   case(2)

      A = (/2,3,3,4/)
      B = (/3,2,4,4/)
      C = (/3,3,3,3/)

      do j=1,3
         do i=1,2
            A%a(i,j) = i
            B%a(j,i) = j
         end do
      end do

      call Times(C, A, B)
      call Print(C, 10, 2)

      C = A*B
      call Print(C, 10, 2)

      C = B*A
      call Print(C, 10, 2)

      B = (/1,1, 2,4/)
      B%a(1,1) = 3
      C = B*A
      call Print(C, 10, 2)
      C = A*B
      call Print(C, 10, 2)

      A = (/5, 0, 10, 10/)
      B = (/0, 5, 10, 10/)
      C = A*B
      call Print(C, 10, 2)
      C = B*A
      call Print(C, 10, 2)

   case(3)


      A = (/2,3,3,4/)
      B = (/2,3,4,4/)
      C = (/3,3,3,3/)

      do j=1,3
         do i=1,2
            A%a(i,j) = i
            B%a(i,j) = i+1
         end do
      end do

      call TimesXhy(C, A, B)
      call Print(C, 10, 2)

      call TimesXyh(C, A, B)
      call Print(C, 10, 2)

      C = A.xhy.B
      call Print(C, 10, 2)
      
      C = A.xyh.B
      call Print(C, 10, 2)

      A = (/5, 0, 10, 10/)
      B = (/5, 0, 10, 10/)
      C = A.xhy.B
      call Print(C, 10, 2)
      C = A.xyh.B
      call Print(C, 10, 2)

   case(4)

      A = (/2,3,3,4/)
      C = (/3,3,3,3/)

      do j=1,3
         do i=1,2
            A%a(i,j) = i
         end do
      end do

      call TimesXhx(C, A)
      call Print(C, 10, 2)

      call TimesXxh(C, A)
      call Print(C, 10, 2)

      C = .xhx.A
      call Print(C, 10, 2)

      C = .xxh.A
      call Print(C, 10, 2)

      A = (/5, 0, 10, 10/)
      C = .xhx.A
      call Print(C, 10, 2)
      C = .xxh.A
      call Print(C, 10, 2)


   case default

      print *, 'No such test case.'

   end select
   

end program RmatProductDebug
