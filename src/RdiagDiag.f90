module RdiagDiag_m

use MatranUtil_m
use Rmat_m
use Rdiag_m

implicit none

#ifdef OVERVIEW

   Rdiag extracts diagonals from a Rmat.  The routines are

   Diag
      subroutine RdDiag(D, A, k): D = k-th diagonal of A (k<0
                                      subdiagonal, k>0 superdiagonal)

   .diag.
       function RdDiag_o(A) result(D): D = the principal diagonal
                                       of A

Author: Pete Stewart
Jun 17 2003

#endif

   interface Diag
      module procedure RdDiag
   end interface

   interface operator (.diag.)
      module procedure RdDiag_o
   end interface

contains

!  Extract the k-th diagonal of A.

   subroutine RdDiag(D, A, k)
      type(Rdiag), intent(inout) :: D
      type(Rmat), intent(in)     :: A
      integer, optional          :: k

      integer i, kd, order

      call GuardTemp(A)

      if (present(k)) then
         kd = k
      else
         kd = 0
      end if

      if (kd >= 0) then
         if (kd >= A%ncol) then
            call MatranError("RdDiag in Rdiag: Illegal diagonal.")
         end if
         order = min(A%nrow, A%ncol-kd)
         call ReshapeAry(D, order)
         do i=1, order
            D%a(i) = A%a(i,i+kd)
         end do
      else
         kd = -kd
         if (kd >= A%nrow) then
            call MatranError("RdDiag in Rdiag: Illegal diagonal.")
         end if

         order = min(A%ncol, A%nrow-kd)
         call ReshapeAry(D, order)
         do i=1,order
            D%a(i) = A%a(i+kd,i)
         end do
      end if

      call CleanTemp(A)
   end subroutine RdDiag

   ! D = the principal diagonal of A

   function RdDiag_o(A) result(D)
      type(Rmat), intent(in) :: A
      type(Rdiag) :: D

      D%a => null()
      D%temporary => null()
      call Clean(D)
      call Diag(D, A)
      call SetTemp(D)
   end function RdDiag_o

end module RdiagDiag_m





