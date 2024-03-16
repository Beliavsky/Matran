program RdiagDiagDebug

use MatranUtil_m
use Rdiag_m
use Rmat_m
use RdiagDiag_m

implicit none

      integer :: i, j, m, n
      type(Rmat) :: A
      type(Rdiag) :: D

      m = 4
      n = 5

      A = (/m,n/)

      do j=1,n
         do i=1,m
            A%a(i,j) = i+j-1
         end do
      end do

      D = .diag.A
      print *, D%a, D%order, D%na, D%adjustable, associated(D%temporary)

      call Clean(D)
      call Diag(D, A)
      print *, D%a, D%order, D%na, D%adjustable, associated(D%temporary)

      call Clean(D)
      call Diag(D, A, 2)
      print *, D%a, D%order, D%na, D%adjustable, associated(D%temporary)

      call Clean(D)
      call Diag(D, A, -1)
      print *, D%a, D%order, D%na, D%adjustable, associated(D%temporary)

      call Clean(D)
      call Diag(D, A, 5)
      print *, D%a, D%order, D%na, D%adjustable, associated(D%temporary)

end program RdiagDiagDebug

