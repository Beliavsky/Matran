module RdiagProduct_m

use MatranUtil_m
use Rdiag_m
use Rmat_m

implicit none

#ifdef OVERVIEW

RdiagProduct implements the product of matrices, one of which is
diagonal.  As is true of all Matran operations, the product comes in
two forms: an explicit call to a subroutine implementing the generic
function Times and a function implementing the binary operator *.  The
actual computation takes place in the former, which is called by the
latter.

The subroutines are (here D, E, F are Rdiags, A, B, C are Rmats, and
s is a real scalar)

   Times

      subroutine RdTimesRdRd(F, D, E):  F = D*E
      subroutine RdTimesRdRm(B, D, A):  B = D*A
      subroutine RdTimesRmRd(B, A, D):  B = A*D
      subroutine RdTimesRsRd(E, s, D):  E = s*D

   *
      function RdTimesRdRd(D, E) result(F): F = D*E
      function RdTimesRdRm(D, A) result(B): B = D*A
      function RdTimesRmDa(A, D) result(B): B = A*D
      function RdTimesRsRd(s, D) result(E): E = s*D
      function RdTimesRdRs(D, s) result(E): E = D*s

Author: Pete Stewart
Jun 17 2003

#endif

   interface Times
      module procedure RdTimesRdRd, RdTimesRdRm, RdTimesRmRd, RdTimesRsRd
   end interface Times

   interface operator (*)
      module procedure RdTimesRdRd_o, RdTimesRdRm_o, RdTimesRmRd_o, &
                       RdTimesRsRd_o, RdTimesRdRs_o
   end interface


contains

!  F = D*E

   subroutine RdTimesRdRd(F, D, E)
      type(Rdiag), intent(out) :: F
      type(Rdiag), intent(in) :: D
      type(Rdiag), intent(in) :: E

      integer n

      n = D%order
      if (n /= E%order) then
         call MatranError('RdTimesRdRd in RdiagProduct: &
                           &Dimensions incompatible for multiplication.')
      end if

      call GuardTemp(D)
      call GuardTemp(E)

      call ReshapeAry(F, n)
      F%a(1:n) = D%a(1:n)*E%a(1:n)
      
      call CleanTemp(D)
      call CleanTemp(E)

   end subroutine RdTimesRdRd

   function RdTimesRdRd_o(D, E) result(F)
      type(Rdiag) :: F
      type(Rdiag), intent(in) :: D
      type(Rdiag), intent(in) :: E
      
      call GuardTemp(D)
      call GuardTemp(E)

      F%a => null()
      F%temporary => null()
      call Clean(F)
      call Times(F, D, E)
      call Settemp(F)

      call CleanTemp(D)
      call CleanTemp(E)


   end function RdTimesRdRd_o

!  B = D*A

   subroutine RdTimesRdRm(B, D, A)
      type(Rmat), intent(inout) :: B
      type(Rdiag), intent(in) :: D
      type(Rmat), intent(in) :: A

      integer j, m, n

      m = D%order
      if (m /= A%nrow)&
         call MatranError("RdTimesRdRm_o in RdiagProduct: &
                           &Nonconforming dimensions.")

      call GuardTemp(D)
      call GuardTemp(A)

      n = A%ncol

      call ReshapeAry(B, m, n)
      do j=1,n
         B%a(1:m,j) = D%a(1:m)*A%a(1:m,j)
      end do

      if (A%tag=='GE' .or. A%tag=='LT' .or. A%tag=='UT') then
         B%tag = A%tag
      else
         B%tag = 'GE'
      end if

      call CleanTemp(D)
      call CleanTemp(A)

   end subroutine RdTimesRdRm


   function RdTimesRdRm_o(D, A) result(B)
      type(Rmat) :: B
      type(Rdiag), intent(in) :: D
      type(Rmat), intent(in) :: A

      call GuardTemp(D)
      call GuardTemp(A)

      B%a => null()
      B%temporary => null()
      call Clean(B)


      call Times(B, D, A)

      call SetTemp(B)
      call CleanTemp(D)
      call CleanTemp(A)

   end function RdTimesRdRm_o

!  B = AA*DD

   subroutine RdTimesRmRd(B, AA, DD)
      type(Rmat), intent(inout) :: B
      type(Rmat), intent(in) :: AA
      type(Rdiag), intent(in) :: DD

      integer j, m, n

      n = DD%order
      if (n /= AA%ncol)&
         call MatranError("RdTimesRdRm_o in RdiagProduct: &
                           &Nonconforming dimensions.")
      m = AA%nrow

      call GuardTemp(DD)
      call GuardTemp(AA)

      call ReshapeAry(B, m, n)
      do j=1,n
         B%a(1:m,j) = AA%a(1:m,j)*DD%a(j)
      end do

      if (AA%tag=='GE' .or. AA%tag=='LT' .or. AA%tag=='UT') then
         B%tag = AA%tag
      else
         B%tag = 'GE'
      end if

      call CleanTemp(DD)
      call CleanTemp(AA)

   end subroutine RdTimesRmRd

   function RdTimesRmRd_o(A, D) result(B)
      type(Rmat) :: B
      type(Rdiag), intent(in) :: D
      type(Rmat), intent(in) :: A

      call GuardTemp(D)
      call GuardTemp(A)

      B%a => null()
      B%temporary => null()
      call Clean(B)

      call Times(B, A, D)

      call SetTemp(B)
      call CleanTemp(D)
      call CleanTemp(A)

   end function RdTimesRmRd_o

!  E = s*D

   subroutine RdTimesRsRd(E, s, D)
      type(Rdiag), intent(inout) :: E
      real(wp), intent(in) :: s
      type(Rdiag), intent(in) :: D

      integer :: n

      n = D%order

      call GuardTemp(D)

      call ReshapeAry(E, n)
     
      E%a(1:n) = s*D%a(1:n)

      call CleanTemp(D)

   end subroutine RdTimesRsRd

   function RdTimesRsRd_o(s, D) result(E)
      type(Rdiag) :: E
      real(wp), intent(in) :: s
      type(Rdiag), intent(in) :: D

      call GuardTemp(D)

      E%a => null()
      E%temporary => null()
      call Clean(E)

      call Times(E, s, D)

      call SetTemp(E)
      call CleanTemp(D)

   end function RdTimesRsRd_o

   function RdTimesRdRs_o(D, s) result(E)
      type(Rdiag) :: E
      type(Rdiag), intent(in) :: D
      real(wp), intent(in) :: s

      call GuardTemp(D)

      E%a => null()
      E%temporary => null()
      call Clean(E)

      call Times(E, s, D)

      call SetTemp(E)
      call CleanTemp(D)

   end function RdTimesRdRs_o


end module RdiagProduct_m




