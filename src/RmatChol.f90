module RmatChol_m
use MatranUtil_m
use Rmat_m
implicit none

#ifdef OVERVIEW

RmatChol is a module to compute the Cholesky factor of a symmetric
positive definite Rmat.  Specifically, given a symmetric positive
definite matrix A of order n, there is an upper trianguarl matrix
R such that

   A = R'R.

The matrix R is called the Cholesky factor of A.

The decomposition is represented by the type RmatChol defined by

   type RmatChol
      type(Rmat) :: R
      logical :: companion
   end type RmatChol

where

   R            is the Cholesky factor
   companion    is true if the decomposition is associated
                with a Rmat of interest.

The Cholesky decomposition of a Rmat of tag HP is computed by the
generic subroutine Chol, whose instantiation is RmGetChol.

A RmatChol is cleaned by the generic subroutine Clean whose
instatiation is RmCholClean.

Author: Pete Stewart
Oct  1 2003

#endif

   type RmatChol
      type(Rmat) :: R                  ! The Cholesky factor
      logical :: companion = .false.   ! True if the decomposition is
                                       ! associated with a Rmat
                                       ! of interest.
   end type RmatChol

   interface Chol
      module procedure RmChol
   end interface Chol

   interface Clean
      module procedure RmCleanChol
   end interface

contains

   subroutine RmChol(chola, A, info)

      type(RmatChol), intent(inout), target :: chola
         ! The Cholesky decomposition of A

      type(Rmat), intent(in)                :: A
         ! The input matrix

      integer, optional, intent(out) :: info
         ! If present, RmGetChol returns the info parameter from
         ! the matlab subroutine dpotrf in lieu of an error
         ! return when it is not zero.

      integer :: i, inf, j, n
      type(Rmat), pointer :: R

      call GuardTemp(A)

      ! Check tag and dimensions.

      if (A%tag /= 'HP')&
         call MatranError('RmGetChol in RmatChol_m:&
                           & Rmat not of type HP' )

      if (A%nrow /= A%ncol)&
         call MatranError('RmGetChol in RmatChol_m:&
                           & Rmat not square' )

      ! Get space for the Cholesky factor and initialize it to A.

      n = A%nrow
      R => chola%R
      call ReshapeAry(R, n, n)

      R%a(1:n,1:n) = A%a(1:n, 1:n)

      do j=1,n-1
         do i=j+1,n
            R%a(i,j) = 0.0
         end do
      end do

      ! Compute the Cholesky factor

#ifdef dbl
      call dpotrf('U', n, R%a, R%narow, inf)
#endif
#ifdef sngl
      call spotrf('U', n, R%a, R%narow, inf)
#endif

      if (inf /=0) then
         if (present(info)) then
            info = inf
            return
         else
            call SupportError('RmChol in RmatChol_m:&
                           & Error in Lapack routine DSPTRF',inf)
         end if
      end if

      ! Set the tag of R and clean up.

      R%tag = 'UT'

      chola%companion = .true.

      call CleanTemp(A)

   end subroutine RmChol

   ! Clean for type RmatChol

   subroutine RmCleanChol(chl)
      type(RmatChol), intent(inout) :: chl

      call Clean(chl%R)

      chl%companion = .false.

   end subroutine RmCleanChol

end module RmatChol_m

