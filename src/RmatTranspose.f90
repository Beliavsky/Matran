module RmatTranspose_m

use MatranUtil_m
use Rmat_m

#ifdef OVERVIEW

RmatTranspose_m implements the transpose and conjugate transpose of a
Rmat.  They are of course are the same; but if there is any chance the
program using them will be converted to complex arithmetic, the
conjugate transpose should be preferred.

As customary, the operations have subroutine and function forms,
the latter implementing an operator.

   Ctp
      subroutine RmCtp(C, A): C = conjugate transpose of A
   Trp
      subroutine RmCtp(C, A): C = transpose of A
   .ctp.
      function RmCtp_o(A) result(C): C = conjugate transpose of A
   .trp.
      function RmCtp_o(A) result(C): C = transpose of A

Note that for Rmats, Ctp and Trp, produce identical results, but they
will be different for Zmats.  It is recommended that one use the Ctp
forms.

Author: Pete Stewart
May  3 2003

#endif

implicit none

   interface Ctp
      module procedure RmCtp
   end interface

   interface operator(.ctp.)
      module procedure RmCtp_o
   end interface

   interface Trp
      module procedure RmCtp
   end interface

   interface operator(.trp.)
      module procedure RmCtp_o
   end interface


contains

   ! Computes C = A'

   subroutine RmCtp(C, A)
      type(Rmat), intent(inout) :: C
      type(Rmat), intent(in) :: A

      integer :: i, j, m, n

      ! Protect temporary.

      call GuardTemp(A)

      m = A%nrow
      n = A%ncol

      ! Insure that there is enough storage.

      Call ReshapeAry(C, n, m)

      ! Compute the transpose.

      do i=1,m
         do j=1,n
            C%a(j,i) = A%a(i,j)
         end do
      end do

      ! Adjust the tag component.

      select case(A%tag)

      case('LT')
         C%tag = 'UT'

      case('UT')
         C%tag = 'LT'

      case default
         C%tag = A%tag

      end select

      call CleanTemp(A)
   end subroutine RmCtp

   function RmCtp_o(A) result(C)
      type(Rmat) :: C
      type(Rmat), intent(in) :: A

      ! This sequence is necessary because the Sun
      ! Fortran 95 6.2 compiler does not initialize
      ! the results of functions properly.

      C%a => null()
      C%temporary => null()
      call Clean(C)

      call RmCtp(C, A)

      call SetTemp(C)
   end function RmCtp_o

end module RmatTranspose_m
