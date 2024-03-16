module RmatNorm2_m

use MatranUtil_m
use Rmat_m
use RmatProduct_m
use RmatSpec_m

implicit none

#ifdef OVERVIEW

The 2-norm norm2 of a vector x is defined by

   norm2(x) = sqrt(sum x(i)^2).

The 2-norm of a matrix A is defined by

   norm2(A) =    max     norm2(Ax)
              norm2(x)=1

This module defines a generic function, Norm2, to compute the 2-norm
of a Rmat.  It uses the characterization that the square of the
2-norm is the largest eigenvalue of A'A or AA'.

If A is mxn, the computation of its 2-norm requires
O(max(m,n)*min(m,n)^2), or O(n^3) when m=n.  The computation of the
norms in the Norm suite requires only O(mn) operations.  Therefore,
the 2-norm should be computed only if a less expensive norm will not
do.

#endif

   interface Norm2
      module procedure RmNorm2
   end interface

contains

   function RmNorm2(A) result(norm2)

      real(wp) :: norm2
      type(Rmat), intent(in) :: A

      type(Rmat) :: B
      type(RmatSpec) :: E
      integer :: i, j
      real(wp) :: mx, scale

      call GuardTemp(A)

      B = A

      ! Scale B so that its largest element is 1 in magnitude.

      mx = 0.0
      do j=1,B%ncol
         do i=1,B%nrow
            mx = max(mx, abs(B%a(i,j)))
         end do
      end do

      if (mx == 0.0D0) then
         norm2 = mx
         call CleanTemp(A)
         return
      end if

      scale = 1.0/mx
      B%a(1:B%nrow, 1:B%ncol) = scale*B%a(1:B%nrow, 1:B%ncol)

      ! Compute the largest eigenvalue of B'*B or B*B', whichever
      ! is the smaller.

      if (B%nrow >= B%ncol) then
         call Spec(E, .xhx.B)
      else
         call Spec(E, .xxh.B)
      end if

      ! Return the 2-norm.

      norm2 = mx*sqrt(E%D%a(1))

      ! Clean up.

      call CleanTemp(A)

   end function RmNorm2

end module RmatNorm2_m
