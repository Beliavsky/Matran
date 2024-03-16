module RmatSubmatrix_m

use MatranUtil_m
use Rmat_m
implicit none

#ifdef OVERVIEW

This suite provides functions and subroutines for extracting
submatrices from a Rmat.  The generic functions are summarized
in the following table, in which we use Matlab notation to
describe the submatrix to be extracted.

A(i1:i2, j1:j2)
    function Sbm(A, i1, i2, j1, j2)
    subroutine GetSbm(C, A, i1, i2, j1, j2)

A(:, j1, j2)
    function Col(A, j1, j2)
    subroutine GetCol(C, A, j1, j2)

A(:, j)
    function Col(A, j)
    subroutine GetCol(C, A, j)

A(i1:i2, :)
    function Row(A, i1, i2)
    subroutine GetRow(C, A, i1, i2)

A(i, :)
    function Row(A, i)
    subroutine GetRow(C, A, i)

Author: Che Rung Lee
Jun 17 2003


#endif

   interface GetSbm
      module procedure RmGetSbm
   end interface

   interface Sbm
      module procedure RmSbm
   end interface

   interface GetCol
      module procedure RmGetCol
   end interface 


   interface Col
      module procedure RmCol
   end interface

   interface GetRow
      module procedure RmGetRow
   end interface

   interface Row
      module procedure RmRow
   end interface


contains
	
   ! C = A[m1:m2, n1:n2]

   subroutine RmGetSbm(C, A, m1, m2, n1, n2)      
      type(Rmat), intent(out)::  C
      type(Rmat), intent(in)::  A
      integer::  m1, m2, n1, n2

      integer :: m, n

      call GuardTemp(A)

      ! Check indices.

      m = A%nrow
      n = A%ncol

      if (m1<1 .or. m2>m .or. m1>m2 .or. &
          n1<1 .or. n2>n .or. n1>n2) &
         call MatranError("RmatSbm in RmatSubmatrix:&
                            & Indices out of range!")

      ! Get the submatrix

      m = m2 - m1 + 1
      n = n2 - n1 + 1

      call ReshapeAry(C, m, n)
      C%a(1:m,1:n) = A%a(m1:m2, n1:n2)

      ! Adjust the tag

      if(n1 == m1 .and. n2 == m2) then 
         C%tag = A%tag
      else if ((A%tag=='LT'.or.A%tag=='UT') .and. n1==m1) then
         C%tag = A%tag
      else 
         C%tag = 'GE'
      end if

      call CleanTemp(A)
   end subroutine RmGetSbm

   ! C = A[m1:m2, n1:n2]

   function RmSbm(A, m1, m2, n1, n2) result(C)
      type(Rmat)             ::  C
      type(Rmat), intent(in) ::  A
      integer                ::  m1, m2, n1, n2

      call GuardTemp(A)

      C%a => null()
      C%temporary => null()
      call Clean(C)

      call GetSbm(C, A, m1, m2, n1, n2)

      call SetTemp(C)

      call CleanTemp(A)

   end function RmSbm

   ! C = A{m1:m2, :]

   subroutine RmGetRow(C, A, m1, m2)
      type(Rmat), intent(inout)     ::  C
      type(Rmat), intent(in)        ::  A
      integer, intent(in)           ::  m1
      integer, optional, intent(in) ::  m2

      call GuardTemp(A)

      if (present(m2)) then
         call GetSbm(C, A, m1, m2, 1, A%ncol)
      else
         call GetSbm(C, A, m1, m1, 1, A%ncol)
      end if

      call CleanTemp(A)

   end subroutine RmGetRow


   function RmRow(A, m1, m2) result(C)
      type(Rmat)                    ::  C
      type(Rmat), intent(in)        ::  A
      integer, intent(in)           ::  m1
      integer, optional, intent(in) ::  m2

      call GuardTemp(A)

      C%a => null()
      C%temporary => null()
      call Clean(C)

      if (present(m2)) then
         call GetSbm(C, A, m1, m2, 1, A%ncol)
      else
         call GetSbm(C, A, m1, m1, 1, A%ncol)
      end if

      call SetTemp(C)

      call CleanTemp(A)

   end function RmRow

   ! C = A[:, n1:n2]

   subroutine RmGetCol(C, A, n1, n2)
      type(Rmat), intent(inout)     :: C
      type(Rmat), intent(in)        :: A
      integer, intent(in)           :: n1
      integer, optional, intent(in) :: n2

      call GuardTemp(A)

      if (present(n2)) then
         call GetSbm(C, A, 1, A%nrow, n1, n2)
      else
         call GetSbm(C, A, 1, A%nrow, n1, n1)
      end if

      call CleanTemp(A)

   end subroutine RmGetCol


   function RmCol(A, n1, n2) result(C)
      type(Rmat)                    :: C
      type(Rmat), intent(in)        :: A
      integer, intent(in)           :: n1
      integer, optional, intent(in) :: n2

      call GuardTemp(A)

      C%a => null()
      C%temporary => null()
      call Clean(C)

      if (present(n2)) then
         call GetSbm(C, A, 1, A%nrow, n1, n2)
      else
         call GetSbm(C, A, 1, A%nrow, n1, n1)
      end if

      call SetTemp(C)

      call CleanTemp(A)

   end function RmCol


end module RmatSubmatrix_m
