module RdiagSum_m

use MatranUtil_m
use Rdiag_m
use Rmat_m

implicit none

#ifdef OVERVIEW

RdiagSum implements the sum of matrices.  As is true of all Matran
operations, the sum comes in two forms: an explicit call to a
subroutine and a function implementing a binary or unary operator.
The actual computation takes place in the former, which is called by
the latter.  They also differ in the sizes of the arrays concerned.
The subroutine attempts to fit the output into the array of the output
variable.  If it is not possible and the array is adjustable, the
subroutine reallocates the array to the right shape.  The function
always allocates the array.

The routines are (here D, DD, E, F are Rdiags and A, AA, B are Rmats)

   Plus
      subroutine RdPlusRdRd(F, D, E): F = D + E
      subroutine RdPlusRdRm(B, D, A): B = D + A
   +
      function RdPlusRdRd_o(D, E) result(F): F = D + E
      function RdPlusRdRm_o(D, A) result(B): B = D + A
      function RdPlusRmRd_o(A, D) result(B): B = A + D

   Minus
      subroutine RdMinusRdRd(F, D, E): F = D - E
      subroutine RdMinusRdRm(B, D, A): B = D - A
      subroutine RdMinusRmRd(B, AA, DD): B = AA - DD
      subroutine RdMinusRd(E, D): E = -D
   -
      function RdMinusRdRd_o(D, E) result(F): F = D - E
      function RdMinusRdRm_o(D, A) result(B): B = D - A
      function RdMinusRmRd_o(A, D) result(B): B = A - D
      function RdMinusRd_o(D) result(E): E = -D 

Author: Pete Stewart
Jun 17 2003


#endif

   interface Plus
      module procedure RdPlusRdRd, RdPlusRdRm
   end interface Plus

   interface operator (+)
      module procedure RdPlusRdRd_o, RdPlusRdRm_o, RdPlusRmRd_o
   end interface

   interface Minus
      Module procedure RdMinusRdRd, RdMinusRdRm, RdMinusRmRd, &
                       RdMinusRd
   end interface

   interface operator (-)
      module procedure RdMinusRdRd_o, RdMinusRdRm_o, RdMinusRmRd_o, &
                       RdMinusRd_o
   end interface

contains

   ! F = D + E

   subroutine RdPlusRdRd(F, D, E)
      type(Rdiag), intent(out) :: F
      type(Rdiag), intent(in)  :: D
      type(Rdiag), intent(in)  :: E

      integer :: n

      call GuardTemp(D)
      call GuardTemp(E)

      n = D%order
      if (n /= E%order)&
         call MatranError("RdPlusRdRd in RdiagSum: &
                           &Nonconforming orders.")

      call ReshapeAry(F, n)

      F%a(1:n) = D%a(1:n) + E%a(1:n)

      call CleanTemp(D)
      call CleanTemp(E)

   end subroutine RdPlusRdRd

   ! F = D + E

   function RdPlusRdRd_o(D, E) result(F)
      type(Rdiag) :: F
      type(Rdiag), intent(in) :: D
      type(Rdiag), intent(in) :: E

      F%a => null()
      F%temporary => null()
      call Clean(F)

      call Plus(F, D, E)

      call Settemp(F)

   end function RdPlusRdRd_o

   ! B = D + A    

   subroutine RdPlusRdRm(B, D, A)
      type(Rmat), intent(out) :: B
      type(Rdiag), intent(in) :: D
      type(Rmat), intent(in) :: A

      integer n, i

      n = D%order

      if (A%nrow/=n .or. A%ncol/=n)&
         call MatranError("RdPlusRdRm in RdiagSum: &
                           &Nonconforming dimenstions.")
      call GuardTemp(D)
      call GuardTemp(A)

      B = A%a(1:n,1:n)

      do i=1,n
        B%a(i,i) = D%a(i) + B%a(i,i)
      end do

      if (A%tag=='PO') then
         B%tag = 'HE'
      else
         B%tag = A%tag
      end if

      call CleanTemp(D)
      call CleanTemp(A)

   end subroutine RdPlusRdRm

   ! B = D + A

   function RdPlusRdRm_o(D, A) result(B)
      type(Rmat) :: B
      type(Rdiag), intent(in) :: D
      type(Rmat), intent(in) :: A

      B%a => null()
      B%temporary => null()
      call Clean(B)

      call RdPlusRdRm(B, D, A)
      call Settemp(B)

   end function RdPlusRdRm_o

   ! B = A + D

   function RdPlusRmRd_o(A, D) result(B)
      type(Rmat) :: B
      type(Rmat), intent(in) :: A
      type(Rdiag), intent(in) :: D

      B%a => null()
      B%temporary => null()
      call Clean(B)

      call RdPlusRdRm(B, D, A)
      call SetTemp(B)

   end function RdPlusRmRd_o

   ! F = D - E

   subroutine RdMinusRdRd(F, D, E)
      type(Rdiag), intent(out) :: F
      type(Rdiag), intent(in)  :: D
      type(Rdiag), intent(in)  :: E

      integer :: n

      call GuardTemp(D)
      call GuardTemp(E)

      n = D%order
      if (n /= E%order)&
         call MatranError("RdMinusRdRd in RdiagSum: &
                           &Nonconforming orders.")

      call ReshapeAry(F, n)

      F%a(1:n) = D%a(1:n) - E%a(1:n)

      call CleanTemp(D)
      call CleanTemp(E)

   end subroutine RdMinusRdRd
      
   ! F = D - E

   function RdMinusRdRd_o(D, E) result(F)
      type(Rdiag) :: F
      type(Rdiag), intent(in) :: D
      type(Rdiag), intent(in) :: E

      F%a => null()
      F%temporary => null()
      call Clean(F)

      call Minus(F, D, E)
      call SetTemp(F)

   end function RdMinusRdRd_o

   ! B = D - A

   subroutine RdMinusRdRm(B, D, A)
      type(Rmat), intent(out) :: B
      type(Rdiag), intent(in) :: D
      type(Rmat), intent(in) :: A

      integer n, i

      n = D%order

      if (A%nrow/=n .or. A%ncol/=n)&
         call MatranError("RdMinusRdRm in RdiagSum: &
                           &Nonconforming dimenstions.")
      call GuardTemp(D)
      call GuardTemp(A)

      B = -A%a(1:n,1:n)

      do i=1,n
        B%a(i,i) = D%a(i) + B%a(i,i)
      end do

      if (A%tag=='P0') then
         B%tag = 'HE'
      else
         B%tag = A%tag
      end if

      call CleanTemp(D)
      call CleanTemp(A)

   end subroutine RdMinusRdRm

   ! B = D - A

   function RdMinusRdRm_o(D, A) result(B)
      type(Rmat) :: B
      type(Rdiag), intent(in) :: D
      type(Rmat), intent(in) :: A

      B%a => null()
      B%temporary => null()
      call Clean(B)

      call RdMinusRdRm(B, D, A)
      call SetTemp(B)

   end function RdMinusRdRm_o

   ! B = AA - DD

   subroutine RdMinusRmRd(B, AA, DD)
      type(Rmat), intent(out) :: B
      type(Rdiag), intent(in) :: DD
      type(Rmat), intent(in) :: AA

      integer n, i

      n = DD%order

      if (AA%nrow/=n .or. AA%ncol/=n)&
         call MatranError("RdMinusRmRd in RdiagSum: &
                           &Nonconforming dimensions.")
      call GuardTemp(DD)
      call GuardTemp(AA)

      B = AA%a(1:n,1:n)

      do i=1,n
        B%a(i,i) =  B%a(i,i) - DD%a(i)
      end do

      if (AA%tag=='PO') then
         B%tag = 'HE'
      else
         B%tag = AA%tag
      end if

      call CleanTemp(AA)
      call CleanTemp(DD)

   end subroutine RdMinusRmRd

   ! B = A - D

   function RdMinusRmRd_o(A, D) result(B)
      type(Rmat) :: B
      type(Rmat), intent(in) :: A
      type(Rdiag), intent(in) :: D

      B%a => null()
      B%temporary => null()
      call Clean(B)

      call RdMinusRmRd(B, A, D)
      call SetTemp(B)

   end function RdMinusRmRd_o

   ! E = -D

   subroutine RdMinusRd(E, D)
      type(Rdiag), intent(inout) :: E
      type(Rdiag), intent(in) :: D

      integer n

      call GuardTemp(D)
      n = D%order
      call ReshapeAry(E, n)

      E%a(1:n) = -D%a(1:n)

      call CleanTemp(D)

   end subroutine RdMinusRd

   ! E = -D

   function RdMinusRd_o(D) result(E)
      type(Rdiag) :: E
      type(Rdiag), intent(in) :: D

      E%a => null()
      E%temporary => null()
      call Clean(E)

      call Minus(E, D)
      call SetTemp(E)

   end function RdMinusRd_o

end module RdiagSum_m
