module RmatPivot_m
use MatranUtil_m
use Rmat_m
implicit none

#ifdef OVERVIEW

The pivot suite applies interchanges to the rows or columns of a Rmat,
thus effecting a permutation of the rows or columns.  It also applies
the inverse permutation.  The permutation is specified by an array pvt
of length npvt.  The effect of pivoting and its inverse on an array x
is given by the following fragments of pseudo-code.

   Pivoting                             Inverse pivoting

   do i=1 to npvt                       do i=npvt,1,-1
      swap x(i) and x(pvt(i))              swap x(i) and x(pvt(i))
   end do                               end do

The generic routines and their instantiations are

   PivotRow
      subroutine RmPivotRow(A, pvt, npvt)

   PivotInvRow
      subroutine RmPivotInvRow(A, pvt, npvt)

   PivotCol
      subroutine RmPivotCol(A, pvt, npvt)

   PivotInvCol
      subroutine RmPivotInvCol(A, pvt, npvt)

Author: Pete Stewart
Jun 17 2003

#endif


   interface PivotRow
      module procedure RmPivotRow
   end interface

   interface PivotInvRow
      module procedure RmPivotInvRow
   end interface

   interface PivotCol
      module procedure RmPivotCol
   end interface
   
   interface PivotInvCol
      module procedure RmPivotInvCol
   end interface
   
contains

   ! Applies pvt to the row of A.

   subroutine RmPivotRow(A, pvt, npvt)
      type(Rmat), intent(inout) :: A
      integer, intent(in) :: pvt(:)
      integer, intent(in) :: npvt

      integer :: i

      if (npvt > A%nrow) then
         call MatranError('RmPivotRow in RmatPivot_m: Too many pivots.')
      end if

      call GuardTemp(A)

      do i=1,npvt
         if (pvt(i) > A%nrow) then
            call MatranError('RmPivotRow in RmatPivot_m:&
                                &  Illegal pivot index.')
         end if
         if (i /= pvt(i)) then
#ifdef dbl
            call dswap(A%ncol,&
                       A%a(i,1), A%narow, A%a(pvt(i),1), A%narow)
#endif
#ifdef sngl
            call sswap(A%ncol,&
                       A%a(i,1), A%narow, A%a(pvt(i),1), A%narow)
#endif
         end if
      end do
      
      call CleanTemp(A)
   end subroutine RmPivotRow

   ! Applies the inverse of pvt to the rows of A.

   subroutine RmPivotInvRow(A, pvt, npvt)
      type(Rmat), intent(inout) :: A
      integer, intent(in) :: pvt(:)
      integer, intent(in) :: npvt

      integer :: i

      call GuardTemp(A)
      if (npvt > A%nrow) then
         call MatranError('RmPivotInvRow in RmatPivot_m: Too many pivots.')
      end if

      do i=npvt,1,-1
         if (pvt(i) > A%nrow) then
            call MatranError('RmPivotInvRow in RmatPivot_m:&
                                &  Illegal pivot index.')
         end if
         if (i /= pvt(i)) then
#ifdef dbl
            call dswap(A%ncol,&
                       A%a(i,1), A%narow, A%a(pvt(i),1), A%narow)
#endif
#ifdef sngl
            call sswap(A%ncol,&
                       A%a(i,1), A%narow, A%a(pvt(i),1), A%narow)
#endif
         end if
      end do

      call CleanTemp(A)

   end subroutine RmPivotInvRow

   ! Applies pvt to the columns of A.

   subroutine RmPivotCol(A, pvt, npvt)
      type(Rmat), intent(inout) :: A
      integer, intent(in) :: pvt(:)
      integer, intent(in) :: npvt

      integer :: j

      if (npvt > A%ncol) then
         call MatranError('RmPivotCol in RmatPivot_m: Too many pivots.')
      end if

      call GuardTemp(A)

      do j=1,npvt
         if (pvt(j) > A%ncol) then
            call MatranError('RmPivotCol in RmatPivot_m:&
                                &  Illegal pivot index.')
         end if
         if (j /= pvt(j)) then
#ifdef dbl
            call dswap(A%nrow,&
                       A%a(1,j), 1, A%a(1, pvt(j)), 1)
#endif
#ifdef sngl
            call sswap(A%nrow,&
                       A%a(1,j), 1, A%a(1, pvt(j)), 1)
#endif
         end if
      end do

      call CleanTemp(A)
   end subroutine RmPivotCol

   ! Applies the inverse of pvt to the columns of A.

   subroutine RmPivotInvCol(A, pvt, npvt)
      type(Rmat), intent(inout) :: A
      integer, intent(in) :: pvt(:)
      integer, intent(in) :: npvt

      integer :: j

      if (npvt > A%ncol) then
         call MatranError('RmPivotInvCol in RmatPivot_m: Too many pivots.')
      end if

      call GuardTemp(A)

      do j=npvt,1,-1
         if (pvt(j) > A%ncol) then
            call MatranError('RmPivotInvCol in RmatPivot_m:&
                                &  Illegal pivot index.')
         end if
         if (j /= pvt(j)) then
#ifdef dbl
            call dswap(A%nrow,&
                       A%a(1,j), 1, A%a(1,pvt(j)), 1)
#endif
#ifdef sngl
            call sswap(A%nrow,&
                       A%a(1,j), 1, A%a(1,pvt(j)), 1)
#endif
         end if
      end do

      call CleanTemp(A)
   end subroutine RmPivotInvCol


end module RmatPivot_m

