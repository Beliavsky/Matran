module RmatSum_m
use MatranUtil_m
use Rmat_m
implicit none

#ifdef OVERVIEW

RmatSum implements the sum of matrices.  As is true of all Matran
operations, the sum comes in two forms: an explicit call to a
subroutine and a function implementing a binary or unary operator.
The actual computation takes place in the former, which is called by
the latter.  They also differ in the sizes of the arrays concerned.
If the result is an mxn matrix, the subroutine will try to fit it into
the current array C, reshaping it only if necessary.  The functional
form always returns the result in an mxn array.

The routines are

   Plus
      subroutine RmPlusRmRm(C, A, B): C = A + B
   +
      function RmPlusRmRm_o(A, B) result(C): C = A + B

   Minus
      subroutine RmMinusRmRm(C, A, B): C = A - B
      subroutine RmMinusRmRm(C, A): C = -A
   -
      function RmMinusRmRm_o(A, B) result(C): C= A - B
      function RmMinusRm_o(A, B) result(C): C = -A
module RmatSum_m

Author: Pete Stewart
Oct  1 2003

#endif

   interface Plus
      module procedure RmPlusRmRm
   end interface

   interface operator(+)
      module procedure RmPlusRmRm_o
   end interface

   interface Minus
      module procedure RmMinusRmRm, RmMinusRm
   end interface

   interface operator(-)
      module procedure RmMinusRmRm_o, RmMinusRm_o
   end interface

   contains

   ! C = A + B

   ! This routine is prototypical and its comments serve
   ! for the procedures that follow.

   subroutine RmPlusRmRm(C, A, B)
      type(Rmat), intent(inout) :: C
      type(Rmat), intent(in) :: A, B


      integer ::  m, n

      ! Take care of temporaries.

      call GuardTemp(A)
      call GuardTemp(B)

      m = A%nrow
      n = A%ncol

      if (m/=B%nrow .or. n/=B%ncol) then
         call MatranError("PlusRmatRmat in RmatSum: Matrices not &
                          &conformable for summation.")
      end if

      ! Make sure C is large enough to hold the sum.

      call ReshapeAry(C, m, n)

      ! Adjust the type of C.

      if (A%tag == B%tag) then
         C%tag = A%tag
      end if

      ! Compute the sum.

      if (m/=0 .and. n/=0)&
         C%a(1:m,1:n) = A%a(1:m,1:n) + B%a(1:m,1:n)

      call CleanTemp(A)
      call CleanTemp(B)

   end subroutine RmPlusRmRm

   ! C = A + B, overloading +

   function RmPlusRmRm_o(A, B) result(C)
      type(Rmat) :: C
      type(Rmat), intent(in) :: A, B

      ! This sequence is necessary because the Sun
      ! Fortran 95 6.2 compiler does not initialize
      ! the results of functions properly.

      C%a => null()
      C%temporary => null()
      call Clean(C)

      ! Compute the sum.

      call Plus(C, A, B)

      ! Declare the result temporary.

      call SetTemp(C)

   end function RmPlusRmRm_o

   ! C = A - B 

   subroutine RmMinusRmRm(C, A, B)
      type(Rmat), intent(inout) :: C
      type(Rmat), intent(in) :: A, B

      integer ::  m, n

      call GuardTemp(A)
      call GuardTemp(B)
      m = A%nrow
      n = A%ncol

      if (m/=B%nrow .or. n/=B%ncol) then
         call MatranError("MinusRmatRmat in RmatSum: Matrices not &
                          &conformable for summation.")
      end if

      call ReshapeAry(C, m, n)

      if (A%tag==B%tag .and. A%tag/='HP') then
         C%tag = A%tag
      end if

      if (m/=0 .and. n/=0)&
         C%a(1:m,1:n) = A%a(1:m,1:n) - B%a(1:m,1:n)

      call CleanTemp(A)
      call CleanTemp(B)

   end subroutine RmMinusRmRm


   ! C = A - B, overloading '-'

   function RmMinusRmRm_o(A, B) result(C)
      type(Rmat) :: C
      type(Rmat), intent(in) :: A, B

      C%a => null()
      C%temporary => null()
      call Clean(C)
      call Minus(C, A, B)
      call SetTemp(C)

   end function RmMinusRmRm_o

   ! C = -A

   subroutine RmMinusRm(C, A)
      type(Rmat), intent(inout) :: C
      type(Rmat), intent(in) :: A

      integer ::  m, n

      call GuardTemp(A)

      m = A%nrow
      n = A%ncol

      call ReshapeAry(C, m, n)

      if (A%tag/='SP') then
         C%tag = A%tag
      end if

      if (m/=0 .and. n/=0)&
         C%a(1:m,1:n) = -A%a(1:m,1:n)

      call CleanTemp(A)

   end subroutine RmMinusRm

   ! C = -A, overloading '-'

   function RmMinusRm_o(A) result(C)
      type(Rmat) :: C
      type(Rmat), intent(in) :: A

      C%a => null()
      C%temporary => null()
      call Clean(C)
      call Minus(C, A)
      call SetTemp(C)
   end function RmMinusRm_o

end module RmatSum_m
