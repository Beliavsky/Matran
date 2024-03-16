module RmatEye_m
use MatranUtil_m
use Rmat_m
implicit none

#ifdef OVERVIEW

RmatEye_m produces matrices with ones on their principal diagonals.
The subroutine form has two forms.

   call Eye(A, m)
      On return A is an mxm zero Rmat with ones on its
      principal diagonal.

   call Eye(A, m, n)
      On return A is an mxn zero Rmat with ones on its
      principal diagonal.

There is also a functional form.

   Reye(m)
      Returns an mxm zero Rmat with ones on its principal diagonal.

   Reye(m, n)
      Returns an mxn zero Rmat with ones on its principal diagonal.

Author: Pete Stewart
Oct  1 2003

#endif

   interface Eye
      module procedure RmEye
   end interface

   interface Reye
      module procedure Reye
   end interface

contains
      
   ! Sets A to a Rmat with ones on the principal diagonal.

   subroutine RmEye(A, m, n)
      type(Rmat), intent(inout) :: A
      integer                   :: m
      integer, optional         :: n

      integer :: mm, nn, i , min_mn
      
      ! compute the dimensions of the result.

      mm = m
      if(present(n)) then
         nn = n
      else
         nn = m
      end if

      if (mm < 0 .or. nn < 0)&
         call MatranError("Eye in RmatEye: Negative dimension.")

      call ReshapeAry(A, mm, nn)
      
      ! Place ones on the principal diagonal.

      min_mn = min(mm, nn)

      A%a = 0.0
      do i = 1, min_mn
         A%a(i,i) = 1.0
      end do

      if(mm .eq. nn) then
         A%tag = 'HP'
      else
         A%tag = 'GE'
      end if

   end subroutine RmEye
 

   ! Returns a Rmat with ones on the principal diagonal.

   function Reye(m, n) result(A)

      type(Rmat) :: A
      integer :: m
      integer, optional :: n
      integer nn

      A%a => null()
      A%temporary => null()
      call Clean(A)
      
      if (present(n)) then
         nn = n
      else
         nn = m
      end if

      call Eye(A, m, nn)
      
      call SetTemp(A)
      
   end function Reye

end module RmatEye_m
