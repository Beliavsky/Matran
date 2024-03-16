module RdiagSolve_m

use MatranUtil_m
use Rdiag_m
use Rmat_m

implicit none


#ifdef OVERVIEW

RdiagSolve implements the products of a matrix with its inverse.  As is
true of all Matran operations, the product comes in two forms: an
explicit call to a subroutine and a function implementing a binary or
unary operator.  The actual computation takes place in the former,
which is called by the latter.  They also differ in the sizes of the
arrays concerned.  The subroutine will always try to fit its results
in the array of the current output object, reallocating only if
necessary.  The functional form (necessarily) returns the result in an
array allocated for the purpose.


The routines are

   Solve
      subroutine RdSolveDiv(F, D, s): F = D/s
      subroutine RdSolveRdRdXiy(F, D, E): F = inv(D)*E
      subroutine RdSolveRdRmXiy(B, D, A): F = inv(D)*A
      subroutine RdSolveRdRdXyi(F, D, E): F = D*inv(E)
      subroutine RdSolveRmRdXyi(B, A, D): B = A*inv(D)

   /
      function RdSolve_div(D, s) result(F): F = D/s
   .xiy.
      function RdSolveRdRd_xiy(D, E) result(F): F = inv(D)*E
      function RdSolveRdRm_xiy(D, A) result(B): F = inv(D)*E
   .xyi
      function RdSolveRdRd_xyi(D, E) result(F): F = D*inv(E)
      function RdSolveRmRd_xyi(A, D) result(B): B = A*inv(D)

Note the absence of routines for computing inv(A)*D and the like.
This is just a scaled inverse of A and can be computed using
the Inverse suite: .inv.A*D.

Author Pete Stewart
Jun 17 2003

#endif

   interface SolveDiv
      module procedure RdSolveDiv
   end interface

   interface operator(/)
      module procedure RdSolve_div
   end interface

   interface SolveXiy
      module procedure RdSolveRdRdXiy, RdSolveRdRmXiy
   end interface

   interface operator(.xiy.)
      module procedure RdSolveRdRd_xiy, RdSolveRdRm_xiy
   end interface

   interface SolveXyi
      module procedure RdSolveRdRdXyi, RdSolveRmRdXyi
   end interface

   interface operator(.xyi.)
      module procedure RdSolveRdRd_xyi, RdSolveRmRd_xyi
   end interface

contains

   ! F = D/s

   subroutine RdSolveDiv(F, D, s)
      type(Rdiag), intent(inout) :: F 
      type(Rdiag), intent(in) :: D
      real(wp), intent(in) :: s

      real(wp) :: t

      call GuardTemp(D)

      if (s == 0)&
         call MatranError("RdSolveDiv in RdiagSolve: &
                           &Cannot divide by zero.")
      call GuardTemp(D)

      call ReshapeAry(F, D%order)

      t = 1/s
      F%a(1:D%order) = t*D%a(1:D%order)

      call CleanTemp(D)

   end subroutine RdSolveDiv

   ! F = D/s

   function RdSolve_div(D, s) result(F)
      type(Rdiag) :: F
      type(Rdiag), intent(in) :: D
      real(wp), intent(in) :: s

      call GuardTemp(D)

      F%a => null()
      F%temporary => null()
      call Clean(F)
      call SolveDiv(F, D, s)
      call SetTemp(F)

      call CleanTemp(D)

   end function RdSolve_div

   ! F = inv(D)*E

   subroutine RdSolveRdRdXiy(F, D, E)
      type(Rdiag), intent(inout) :: F
      type(Rdiag), intent(in) :: D, E

      integer :: i, n

      call GuardTemp(D)
      call GuardTemp(E)

      n = D%order

      if (E%order /= n)&
         call MatranError("RdSolveRdRdXiy in RdiagSolve: &
                           &Incompatible orders.")

      do i=1,n
         if (D%a(i) .eq. 0)&
             call MatranError("RdSolveRdRdXiy in RdiagSolve: &
                               &Singular matrix.")
      end do

      call ReshapeAry(F, n)

         F%a(1:n) = E%a(1:n)/D%a(1:n)

      call CleanTemp(D)
      call CleanTemp(E)

   end subroutine RdSolveRdRdXiy


   ! F = inv(D)*E

   function RdSolveRdRd_xiy(D, E) result(F)

      type(Rdiag) :: F
      type(Rdiag), intent(in) :: D, E

      call GuardTemp(D)
      call GuardTemp(E)

      F%a => null()
      F%temporary => null()
      call Clean(F)

      call SolveXiy(F, D, E)

      call SetTemp(F)

      call CleanTemp(D)
      call CleanTemp(E)

   end function RdSolveRdRd_xiy

   ! B = inv(D)*A


   subroutine RdSolveRdRmXiy(B, D, A)
      type(Rmat), intent(inout) :: B
      type(Rdiag), intent(in) :: D
      type(Rmat), intent(in) :: A

      real(wp) :: t
      integer :: i, n

      call GuardTemp(D)
      call GuardTemp(A)

      n = D%order

      if (A%nrow /= n)&
         call MatranError("RdSolveRdRmXiy in RdiagSolve: &
                           &Incompatible dimensionns.")

      do i=1,n
         if (D%a(i) .eq. 0)&
             call MatranError("RdSolveRdRdXiy in RdiagSolve: &
                               &Singular matrix.")
      end do

      call ReshapeAry(B, n, A%ncol)


      do i=1,n
         t = 1./D%a(i)
         B%a(i,1:A%ncol) = t*A%a(i,1:A%ncol)
      end do

      if (A%tag=='GE' .or. A%tag=='LT' .or. A%tag=='UT') then
         B%tag = A%tag
      else
         B%tag = 'GE'
      end if

      call CleanTemp(D)
      call CleanTemp(A)

   end subroutine RdSolveRdRmXiy

   ! B = inv(D)*A


   function RdSolveRdRm_xiy(D, A) result(B)

      type(Rmat) :: B
      type(Rdiag), intent(in) :: D
      type(Rmat),intent(in) :: A

      call GuardTemp(D)
      call GuardTemp(A)


      B%a => null()
      B%temporary => null()
      call Clean(B)

      call SolveXiy(B, D, A)
      call SetTemp(B)

      call CleanTemp(D)
      call CleanTemp(A)

   end function RdSolveRdRm_xiy

   ! F = D*inv(E)

   subroutine RdSolveRdRdXyi(F, D, E)
      type(Rdiag), intent(inout) :: F
      type(Rdiag), intent(in) :: D, E

      integer :: i, n

      call GuardTemp(D)
      call GuardTemp(E)

      n = E%order

      if (D%order /= n)&
         call MatranError("RdSolveRdRdXyi in RdiagSolve: &
                           &Incompatible orders.")

      do i=1,n
         if (E%a(i) .eq. 0)&
             call MatranError("RdSolveRdRdXiy in RdiagSolve: &
                               &Singular matrix.")
      end do


      call ReshapeAry(F, n)

         F%a(1:n) = D%a(1:n)/E%a(1:n)

      call CleanTemp(D)
      call CleanTemp(E)

   end subroutine RdSolveRdRdXyi


   ! F = D*inv(E)

   function RdSolveRdRd_xyi(D, E) result(F)

      type(Rdiag) :: F
      type(Rdiag), intent(in) :: D, E

      call GuardTemp(D)
      call GuardTemp(E)

      F%a => null()
      F%temporary => null()
      call Clean(F)
      call SolveXyi(F, D, E)
      call SetTemp(F)

      call CleanTemp(D)
      call CleanTemp(E)

   end function RdSolveRdRd_xyi

   ! B = A*inv(D)


   subroutine RdSolveRmRdXyi(B, A, D)
      type(Rmat), intent(inout) :: B
      type(Rmat), intent(in) :: A
      type(Rdiag), intent(in) :: D

      real(wp) :: t
      integer :: i, n

      call GuardTemp(D)
      call GuardTemp(A)

      n = D%order

      if (A%ncol /= n)&
         call MatranError("RdSolveRdRmXyi in RdiagSolve: &
                           &Incompatible dimensionns.")

      do i=1,n
         if (D%a(i) .eq. 0)&
             call MatranError("RdSolveRdRdXyi in RdiagSolve: &
                               &Singular matrix.")
      end do

      call ReshapeAry(B, A%nrow, n)

      do i=1,n
         t = 1./D%a(i)
         B%a(1:A%nrow, i) = t*A%a(1:A%nrow, i)
      end do

      if (A%tag=='GE' .or. A%tag=='LT' .or. A%tag=='UT') then
         B%tag = A%tag
      else
         B%tag = 'GE'
      end if

      call CleanTemp(D)
      call CleanTemp(A)

   end subroutine RdSolveRmRdXyi

   ! B = A*inv(D)


   function RdSolveRmRd_xyi(A, D) result(B)

      type(Rmat) :: B
      type(Rmat),intent(in) :: A
      type(Rdiag), intent(in) :: D

      call GuardTemp(D)
      call GuardTemp(A)

      B%a => null()
      B%temporary => null()
      call Clean(B)

      call SolveXyi(B, A, D)

      call SetTemp(B)

      call CleanTemp(D)
      call CleanTemp(A)

   end function RdSolveRmRd_xyi


end module RdiagSolve_m

