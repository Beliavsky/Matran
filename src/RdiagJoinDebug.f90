Program RdiagJoinDebug

use MatranUtil_m
use Rdiag_m
use RdiagJoin_m
use Rmat_m
use RmatPrint_m

implicit none

   type(Rmat) :: A
   type(Rdiag) :: D, E

   integer :: tn

   A = (/3, 4, 5, 6/)
   A%a(1:3, 1:4) = 1

   print *, "Enter test number:"
   read *, tn
   print *, 'Test number ', tn

   select case (tn)

   case(1) ! test JoinWE

      D = (/3, 4/)
      D%a(1:3) = 2
      E = (/3, 3/)
      E%a(1:3) = 3

      call Print(D.jwe.E, 9, 1)
      call Print(A.jwe.D, 9, 1)
      call Print(D.jwe.A, 9, 1)

   case(2) ! test JoinNS

      D = (/4, 5/)
      D%a(1:4) = 2
      E = (/4, 4/)
      E%a(1:4) = 3

      call Print(D.jns.E, 9, 1)
      call Print(A.jns.D, 9, 1)
      call Print(D.jns.A, 9, 1)

   case default

      print *, 'No such case.'

   end select

end program RdiagJoinDebug
