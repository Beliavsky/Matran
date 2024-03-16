module RmatLudpp_m

use MatranUtil_m
use Rmat_m
implicit none

#ifdef OVERVIEW

RmatLudpp is a module to define and compute the partially pivoted LU
dcomposition of a rectangular Rmat.  Specifically, given an mxn
matrix $A$, there is a permutation matrix P and lower and upper
triangular (actually trapezoidal) matrices L and U such that

   PA = LU.

If m<=n, then L is mxm and U is mxn; otherwise L is mxn and U is nxn.
The matrix L is unit lower triangular (i.e., it has ones on its
diagona) and all its subdiagonal elements are less than or equal to
one in magnitude.  The matrix U may have zeros on its diagonal.

The decomposition is represented by the type RmatLudpp defined by

   type RmatLudpp
      type(Rmat) :: L
      type(Rmat) :: U
      integer, pointer :: pvt(:)
      integer :: npvt
      logical :: companion
   end type RmatLudpp

where

   L         is the lower triangular factor.
   U         is the upper triagular factor.
   pvt       defines the permutation matrix P (see below).
   npvt      is the number of interchanges in the representation of P,
             specifically, min(m,n).
   companion is true if the decomposition is that of an existing
             matrix.

The permutation P is defined in terms of interchange (aka pivot)
operations.  Specically, the matrix PA is obtained by swapping rows
A(:,k) and A(:,pvt(k)) for k = 1,...,npvt.  Also see RmatPivot.

The LU decomposition of a Rmat is computed by the generic
subroutine Ludpp, whose instantiation is RmGetLudpp.

A RmatLudpp is cleaned by the generic subroutine Clean whose
instantiation is RmLudppClean.

Author: Pete Stewart
Aug 21 2003

#endif

! LU decomposition of a Rmat.

type RmatLudpp
   type(Rmat) :: L                      ! The L-factor
   type(Rmat) :: U                      ! The U-factor
   integer, pointer :: pvt(:) => null() ! The pivot arry
   integer :: npvt = 0                  ! The number of pivots.
   logical :: companion = .false.       ! True if is the decomposition
                                        ! of an existing Rmat.
end type RmatLudpp

! Ludpp computes the LU factorization.

interface Ludpp
   module procedure RmLudpp
end interface

interface Clean
   module procedure RmLudppClean
end interface

contains

   ! RmGetLudpp computes the LU decomposition of a Rmat using
   ! the Lapack routine dgetrf.  Overloads Ludpp.

   subroutine RmLudpp(lu, A, info)

      type(RmatLudpp), intent(inout), target :: lu
         ! On return contains the LU
         ! decomposition of A.  If on input lu contains storage
         ! for L and U, RmatgetLudpp will attempt to use it. 

      type(Rmat), intent(in) :: A
         ! The Rmat whose whose LU decomposition is to be computed.

      integer, intent(out), optional :: info
         ! The info parameter from the Lapack routine Dgetrf.
         ! info = 0 is the normal return.
         ! If info < 0, there is an error in the calling sequation.
         ! If info > 0, then the info diagonal element of U is
         !    zero and is the first such element.

      integer :: i, j, m, n, inf
      integer, pointer :: npvt
      type(Rmat), pointer :: L, U

      call GuardTemp(A)

      ! Some convenient abbreviations.

      m = A%nrow
      n = A%ncol


      L => lu%L
      U => lu%U
      npvt => lu%npvt

      ! Set type and npvt

      L%tag = 'LT'
      U%tag = 'UT'
      npvt = min(m,n)

      ! There are two cases.
      ! 1. If m>=n copy A into L and perform the
      !    factorization there.
      !
      ! 2. If m>=n copy A into U and perform the
      !    factorization there.


      if (m >= n) then

         ! Get L and U and copy into L.

         call ReshapeAry(L, m, n)
         call ReshapeAry(U, n, n)
         L%a(1:m,1:n) = A%a(1:m,1:n)

         ! Set up the pivot array

         if (associated(lu%pvt)) then
            if (npvt > size(lu%pvt)) then
               deallocate(lu%pvt)
               allocate(lu%pvt(lu%npvt))
            end if
         else
            allocate(lu%pvt(npvt))
         end if

         ! Factor A.

#ifdef dbl
         call dgetrf(m, n, L%a, L%narow, lu%pvt, inf)
#endif
#ifdef sngl
         call sgetrf(m, n, L%a, L%narow, lu%pvt, inf)
#endif

         ! Return or stop if there was an Lapack exception.

         if (present(info)) then
            info = inf
            if (inf < 0) return
         else
            if (inf < 0)&
                 call SupportError('RmGetLudpp in RmatLudpp: Bad argument&
                 & in Lapack routine DGETRF', inf)
         end if

         ! Unbundle L into L and U.

         do j=1,n
            do i=1,j
               U%a(i,j) = L%a(i,j)
               if (i/=j) then
                  L%a(i,j) = 0
               else
                  L%a(i,j) = 1.0D0
               end if
            end do
         end do

      else ! m > n

         ! Get L and U and copy A into U.

         call ReshapeAry(L, m, m)
         call ReshapeAry(U, m, n)
         U%a(1:m,1:n) = A%a(1:m, 1:n)

         ! Set up the pivot array.

         if (associated(lu%pvt)) then
            if (npvt>size(lu%pvt)) then
               deallocate(lu%pvt)
               allocate(lu%pvt(npvt))
            end if
         else
            allocate(lu%pvt(npvt))
         end if

         ! Factor A.


#ifdef dbl
         call dgetrf(m, n, U%a, U%narow, lu%pvt, inf)
#endif
#ifdef sngl
         call sgetrf(m, n, U%a, U%narow, lu%pvt, inf)
#endif

         ! Return or stop if there was an exception.

         if (present(info)) then
            info = inf
            if (inf < 0) return
         else
            if (inf < 0)&
               call SupportError('RmGetLudpp in RmatLudpp: Bad argument&
                            & in Lapack routine DGETRF', inf)
         end if

         ! Unbundle U into L and U.

         do j=1,m
            L%a(j,j) = 1.0D0
            do i=j+1,m
               L%a(i,j) = U%a(i,j)
               U%a(i,j) = 0
            end do
         end do
      end if

      lu%companion = .true.

      call CleanTemp(A)

   end subroutine RmLudpp

   ! Clean for type RmatLudpp

   subroutine RmLudppClean(lupp)
      type(RmatLudpp), intent(inout) :: lupp


      call Clean(lupp%L)
      call Clean(lupp%U)

      if (associated(lupp%pvt)) &
         deallocate(lupp%pvt)


      lupp%npvt = 0
      lupp%companion = .false.

   end subroutine RmLudppClean



end module RmatLudpp_m
