module RdiagJoin_m

use MatranUtil_m
use Rmat_m
use Rdiag_m

implicit none

#ifdef OVERVIEW

RdiagJoin implements the joining of two matrices to create a new
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

Author: Pete Stewart
Jun 17 2003

#endif


   interface JoinWE
      module procedure RdJoinWE_RdRd
      module procedure RdJoinWE_RmRd
      module procedure RdJoinWE_RdRm
   end interface

   interface operator(.jwe.)
      module procedure RdJoin_WE_RdRd
      module procedure RdJoin_WE_RmRd
      module procedure RdJoin_WE_RdRm
   end interface
   
   interface JoinNS
      module procedure RdJoinNS_RdRd
      module procedure RdJoinNS_RmRd
      module procedure RdJoinNS_RdRm
   end interface
   
   interface operator(.jns.)
      module procedure RdJoin_NS_RdRd
      module procedure RdJoin_NS_RmRd
      module procedure RdJoin_NS_RdRm
   end interface
   
contains

   ! C = [D, E]

   subroutine RdJoinWE_RdRd(C, D, E)
      type(Rmat), intent(out) :: C
      type(Rdiag), intent(in) :: D, E

      integer :: i

      call GuardTemp(D)
      call GuardTemp(E)

      ! Check Dimensions

      if (D%order /= E%order)&
         call MatranError('RmJoinWE_RdRd in RdiagJoin: &
                           &Incompatible dimensions.')

      ! Get storage for result.

      call ReshapeAry(C, D%order, D%order+E%order)

      ! Adjust the type

      C%tag = 'UT'

      ! Join the matrices

      do i=1, D%order
         C%a(i,i) = D%a(i)
         C%a(i,i+D%order) = E%a(i)
      end do

      ! Clean up

      call CleanTemp(D)
      call CleanTemp(E)

   end subroutine RdJoinWE_RdRd

   ! Implements D.jwe.E = [D, E]

   function RdJoin_WE_RdRd(D, E) result(C)
      type(Rmat) :: C
      type(Rdiag), intent(in) :: D, E

      call GuardTemp(D)
      call GuardTemp(E)

      C%a => null()
      C%temporary => null()
      call Clean(C)
      call JoinWE(C, D, E)

      call SetTemp(C)

      call CleanTemp(D)
      call CleanTemp(E)
   end function RdJoin_WE_RdRd

   ! C = [A, D]

   subroutine RdJoinWE_RmRd(C, A, D)
      type(Rmat), intent(out) :: C
      type(Rmat), intent(in) :: A
      type(Rdiag), intent(in) :: D

      integer :: i

      call GuardTemp(A)
      call GuardTemp(D)

      ! Check Dimensions

      if (A%nrow /= D%order) &
         call MatranError('RmJoinWE_RmRd in RdiagJoin: &
                           &Incompatible dimensions.')

      ! Get storage for result.

      call ReshapeAry(C, A%nrow, A%ncol+D%order)

      ! Adjust the type

      if (A%tag == 'UT') then
         C%tag = 'UT'
      else
         C%tag = 'GE'
      end if
      
      ! Join the matrices

      C%a(1:A%nrow,1:A%ncol) = A%a(1:A%nrow,1:A%ncol)
      do i=1, D%order
         C%a(i,i+A%ncol) = D%a(i)
      end do

      ! Clean up

      call CleanTemp(A)
      call CleanTemp(D)

   end subroutine RdJoinWE_RmRd

   ! Implements D.jwe.E = [A, D]

   function RdJoin_WE_RmRd(A, D) result(C)
      type(Rmat) :: C
      type(Rmat), intent(in) :: A
      type(Rdiag), intent(in) :: D

      call GuardTemp(A)
      call GuardTemp(D)

      C%a => null()
      C%temporary => null()
      call Clean(C)
      call JoinWE(C, A, D)

      call SetTemp(C)

      call CleanTemp(A)
      call CleanTemp(D)
   end function RdJoin_WE_RmRd

   ! C = [DD, A]

   subroutine RdJoinWE_RdRm(C, DD, A)
      type(Rmat), intent(out) :: C
      type(Rdiag), intent(in) :: DD
      type(Rmat), intent(in) :: A

      integer :: i

      call GuardTemp(A)
      call GuardTemp(DD)

      ! Check Dimensions

      if (A%nrow /= DD%order) &
         call MatranError('RmJoinWE_RdRm in RdiagJoin: &
                           &Incompatible dimensions.')

      ! Get storage for result.

      call ReshapeAry(C, A%nrow, A%ncol+DD%order)

      ! Adjust the type

      C%tag = 'UT'
      
      ! Join the matrices

      do i=1, DD%order
         C%a(i,i) = DD%a(i)
      end do
      C%a(1:A%nrow, 1+DD%order:DD%order+A%ncol) = A%a(1:A%nrow, 1:A%ncol)

      ! Clean up

      call CleanTemp(A)
      call CleanTemp(DD)

   end subroutine RdJoinWE_RdRm

   ! Implements D.jwe.E = [D, A]

   function RdJoin_WE_RdRm(D, A) result(C)
      type(Rmat) :: C
      type(Rmat), intent(in) :: A
      type(Rdiag), intent(in) :: D

      call GuardTemp(A)
      call GuardTemp(D)

      C%a => null()
      C%temporary => null()
      call Clean(C)
      call JoinWE(C, D, A)

      call SetTemp(C)

      call CleanTemp(A)
      call CleanTemp(D)
   end function RdJoin_WE_RdRm

   ! C = [D; E]

   subroutine RdJoinNS_RdRd(C, D, E)
      type(Rmat), intent(out) :: C
      type(Rdiag), intent(in) :: D, E

      integer :: i

      call GuardTemp(D)
      call GuardTemp(E)

      ! Check Dimensions

      if (D%order /= E%order)&
         call MatranError('RmJoinNS_RdRd in RdiagJoin: &
                           &Incompatible dimensions.')

      ! Get storage for result.

      call ReshapeAry(C, D%order+E%order, D%order)

      ! Adjust the type

      C%tag = 'LT'

      ! Join the matrices

      do i=1, D%order
         C%a(i,i) = D%a(i)
         C%a(i+D%order,i) = E%a(i)
      end do

      ! Clean up

      call CleanTemp(D)
      call CleanTemp(E)

   end subroutine RdJoinNS_RdRd

   ! Implements D.jns.E = [D; E]

   function RdJoin_NS_RdRd(D, E) result(C)
      type(Rmat) :: C
      type(Rdiag), intent(in) :: D, E

      call GuardTemp(D)
      call GuardTemp(E)

      C%a => null()
      C%temporary => null()
      call Clean(C)
      call JoinNS(C, D, E)

      call SetTemp(C)

      call CleanTemp(D)
      call CleanTemp(E)
   end function RdJoin_NS_RdRd

   ! C = [A; D]

   subroutine RdJoinNS_RmRd(C, A, D)
      type(Rmat), intent(out) :: C
      type(Rmat), intent(in) :: A
      type(Rdiag), intent(in) :: D

      integer :: i

      call GuardTemp(A)
      call GuardTemp(D)

      ! Check Dimensions

      if (A%ncol /= D%order) &
         call MatranError('RmJoinNS_RmRd in RdiagJoin: &
                           &Incompatible dimensions.')

      ! Get storage for result.

      call ReshapeAry(C, A%nrow+D%order, A%ncol)

      ! Adjust the type

      if (A%tag == 'LT') then
         C%tag = 'LT'
      else
         C%tag = 'GE'
      end if
      
      ! Join the matrices

      C%a(1:A%nrow,1:A%ncol) = A%a(1:A%nrow,1:A%ncol)
      do i=1, D%order
         C%a(i+A%nrow, i) = D%a(i)
      end do

      ! Clean up

      call CleanTemp(A)
      call CleanTemp(D)

   end subroutine RdJoinNS_RmRd

   ! Implements D.jns.E = [A; D]

   function RdJoin_NS_RmRd(A, D) result(C)
      type(Rmat) :: C
      type(Rmat), intent(in) :: A
      type(Rdiag), intent(in) :: D

      call GuardTemp(A)
      call GuardTemp(D)

      C%a => null()
      C%temporary => null()
      call Clean(C)
      call JoinNS(C, A, D)

      call SetTemp(C)

      call CleanTemp(A)
      call CleanTemp(D)
   end function RdJoin_NS_RmRd

   ! C = [DD; A]

   subroutine RdJoinNS_RdRm(C, DD, A)
      type(Rmat), intent(out) :: C
      type(Rdiag), intent(in) :: DD
      type(Rmat), intent(in) :: A

      integer :: i

      call GuardTemp(A)
      call GuardTemp(DD)

      ! Check Dimensions

      if (A%ncol /= DD%order) &
         call MatranError('RmJoinNS_RdRm in RdiagJoin: &
                           &Incompatible dimensions.')

      ! Get storage for result.

      call ReshapeAry(C, A%nrow+DD%order, A%ncol)

      ! Adjust the type

      C%tag = 'LT'
      
      ! Join the matrices

      do i=1, DD%order
         C%a(i,i) = DD%a(i)
      end do
      C%a(1+DD%order:DD%order+A%nrow, 1:A%ncol) = A%a(1:A%nrow, 1:A%ncol)

      ! Clean up

      call CleanTemp(A)
      call CleanTemp(DD)

   end subroutine RdJoinNS_RdRm

   ! Implements D.jns.E = [D, A]

   function RdJoin_NS_RdRm(D, A) result(C)
      type(Rmat) :: C
      type(Rmat), intent(in) :: A
      type(Rdiag), intent(in) :: D

      call GuardTemp(A)
      call GuardTemp(D)

      C%a => null()
      C%temporary => null()
      call Clean(C)
      call JoinNS(C, D, A)

      call SetTemp(C)

      call CleanTemp(A)
      call CleanTemp(D)
   end function RdJoin_NS_RdRm

end module RdiagJoin_m
