module RmatJoin_m

use MatranUtil_m
use Rmat_m
implicit none

#ifdef OVERVIEW

RmatJoin implements the joining of two matrices to create a new
matrix.  Specifically, it consists of the following functions.
The results are expressed in Matlab notation.

   JoinWE
     subroutine RmJoinWE(C, A, B): C = [A, B]
   .jwe.
     function RmJoin_WE(A, B) result(C): C = [A, B]

   JoinNS
     subroutine RmJoinNS(C, A, B): C = [A; B]
   .jns.
     function RmJoin_NS(A, B) result(C): C = [A; B]

Author: Che Rung Lee, Pete Stewart
Jun 17 2003

#endif

   interface JoinWE
      module procedure RmJoinWE
   end interface

   interface operator(.jwe.)
      module procedure RmJoin_WE
   end interface
   
   interface JoinNS
      module procedure RmJoinNS
   end interface
   
   interface operator(.jns.)
      module procedure RmJoin_NS
   end interface
   
contains
   
   ! C = [A, B]

   subroutine RmJoinWE(C, A, B)
      type(Rmat), intent(out) :: C
      type(Rmat), intent(in) :: A
      type(Rmat), intent(in) :: B
      
      call GuardTemp(A)
      call GuardTemp(B)

      ! Check Dimensions.
      if (A%nrow /= B%nrow) &
         call MatranError('RmJoinWE in RmatJoin: &
                               &Incompatible dimensions.')


      ! Prepare the result.
      call ReshapeAry(C, A%nrow, A%ncol+B%ncol)

      ! Adjust the type.

      if (A%tag == 'UT') then
         C%tag = 'UT'
      else
         C%tag = 'GE'
      end if

      ! Join the matrices.

      C%a(1:A%nrow,1:A%ncol) = A%a(1:A%nrow,1:A%ncol)
      C%a(1:A%nrow,A%ncol+1:A%ncol+B%ncol) = B%a(1:B%nrow,1:B%ncol)

      ! Clean up.
      call CleanTemp(A)
      call CleanTemp(B)

   end subroutine RmJoinWE


   ! Implements A.jwe.B = [A, B]

   function RmJoin_WE(A, B) result(C)
      type(Rmat) :: C
      type(Rmat), intent(in) :: A, B

      call GuardTemp(A)
      call GuardTemp(B)

      C%a => null()
      C%temporary => null()
      call Clean(C)

      call JoinWE(C, A, B)

      call SetTemp(C)

      call CleanTemp(A)
      call CleanTemp(B)
   end function RmJoin_WE

   ! C = [A; B]

   subroutine RmJoinNS(C, A, B)
      type(Rmat), intent(out) :: C
      type(Rmat), intent(in) :: A
      type(Rmat), intent(in) :: B
      
      call GuardTemp(A)
      call GuardTemp(B)

      if (A%ncol /= B%ncol) &
         call MatranError('RmJoinTb in RmatJoin: &
                               &Incompatible dimensions.')


      call ReshapeAry(C, A%nrow+B%nrow, A%ncol)

      if (A%tag == 'LT') then
         C%tag = 'LT'
      else
         C%tag = 'GE'
      end if

      C%a(1:A%nrow,1:A%ncol) = A%a(1:A%nrow,1:A%ncol)
      C%a(A%nrow+1:A%nrow+B%nrow, 1:A%ncol) = B%a(1:B%nrow,1:B%ncol)

      call CleanTemp(A)
      call CleanTemp(B)
   end subroutine RmJoinNS

   ! Implements A.jtb.B = [A; B]

   function RmJoin_NS(A, B) result(C)
      type(Rmat) :: C
      type(Rmat), intent(in) :: A, B

      call GuardTemp(A)
      call GuardTemp(B)

      C%a => null()
      C%temporary => null()
      call Clean(C)

      call JoinNS(C, A, B)

      call SetTemp(C)

      call CleanTemp(A)
      call CleanTemp(B)
   end function RmJoin_NS


end module RmatJoin_m
