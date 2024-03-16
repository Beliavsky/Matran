module RmatSpec_m

use MatranUtil_m
use Rmat_m
use Rdiag_m
implicit none

#ifdef OVERVIEW

Let A be a symmetric matrix of order n.  Then there is an
orthogonal matrix V such that

   A = VDV`   (*)

where

   D = diag(d1, ..., dn),   d1 >= ... >= dn.

The decomposition is called the spectral decomposition of A.  The
module RmatSpec_m provides means for computing the spectral
decomposition of a Rmat.  The decomposition is contained in a derived
type RmatSpec, which is defined below.  The spectral decomposition
is compute by the generic subroutine Spec, which is also defined
below.

Author: CheRung Lee
Oct  1 2003

#endif

   type RmatSpec
      type(Rdiag) :: D          ! Contains the eigenvalues
      type(Rmat)  :: V          ! Contains the eigenvectors
      logical     :: companion  ! True if the decomposition is
                                ! associated with a matrix of interest
   end type RmatSpec

   interface Spec
      module procedure RmSpec
   end interface

   interface Clean
      module procedure RmCleanSpec
   end interface

contains

   ! Clean for type RmatSpec.

   subroutine RmCleanSpec(s)
      type(RmatSpec), intent(inout) :: s

      call Clean(s%V)
      call Clean(s%D)

      s%companion = .false.
   end subroutine RmCleanSpec



   subroutine RmSpec(E, A, ev, info, mywork)

      type(RmatSpec), intent(out) :: E
         ! Contains the spectral decomposition of A.

      type(Rmat), intent(in) :: A
         ! The symmatric Rmat whose spectral decomposition 
         ! is to be computed.

      logical, optional, intent(in)  :: ev
         ! If present and true, compute eigenvectors.
         ! Otherwise only eigenvalues computed.

      integer, optional, intent(out) :: info
         ! If present, return dgeev's info parameter in lieu of
         ! an error return.

      real(wp), pointer, optional :: mywork(:)
         ! Optional work array for dsyev.  If mywork is present
         ! it will be used, possibly after a reallocation, to
         ! obtain enough memory.  mywork is not deallocated
         ! by RmSpec.


      integer :: i, j, n, inf

      logical :: cev, usemywork, usetemp

      character(1) :: jobz

      real(wp), pointer :: work(:)=>null()

      real(wp) :: ts, twork(1)
      integer  :: lwork


      call GuardTemp(A)

      ! Check the Rmat A.

      n = A%nrow
      if (n /= A%ncol)&
         call MatranError('RmGetSpec in RmatSpec: Matrix not square.')

      if (A%tag /= 'HP' .and. A%tag /= 'HE') &
         call MatranError('RmGetSpec in RmatSpec: Matrix not symmetric.')
      
      ! Set up storage for the eigenvalues/eigenvectors

      call ReshapeAry(E%D, n)

      ! Decide whether to compute V.

      cev  = .false.
      jobz = 'N'
      if (present(ev)) then
         if(ev) then
            cev = .true.
            jobz = 'V'
         end if
      end if
      E%V = A

      ! Set up the work array.

#ifdef dbl
      call dsyev(jobz, 'L', n, E%V%a, E%V%narow, E%D%a, &
                    twork, -1, inf)
#endif
#ifdef sngl
      call ssyev(jobz, 'L', n, E%V%a, E%V%narow, E%D%a, &
                    twork, -1, inf)
#endif

      if (inf /= 0) go to 50    ! Error.
      lwork = twork(1)

      if (present(mywork)) then
         usemywork = .true.
         call ReshapeAry(mywork, lwork)
         work => mywork
      else
         usemywork = .false.
         call ReshapeAry(work, lwork)
      end if


      ! Compute the decomposition.

#ifdef dbl
      call dsyev(jobz, 'L', n, E%V%a, E%V%narow, E%D%a, &
                 work, lwork, inf)
#endif
#ifdef sngl
      call ssyev(jobz, 'L', n, E%V%a, E%V%narow, E%D%a, &
                 work, lwork, inf)
#endif


      if (inf /= 0) go to 50


      ! Arrange eigenvalues in descending order.

      if (n > 1) then
         do i=1,n/2
            ts = E%D%a(i)
            E%D%a(i) = E%D%a(n-i+1)
            E%D%a(n-i+1) = ts
            if (cev) &
#ifdef dbl
               call dswap(n, E%V%a(1,i), 1, E%V%a(1,n-i+1), 1)
#endif
#ifdef sngl
               call sswap(n, E%V%a(1,i), 1, E%V%a(1,n-i+1), 1)
#endif
         end do
      end if

      E%companion = .true.
      go to 100

50    if (.not. present(info)) then
         call SupportError('RmGetSpec in RmatSpec: &
              &Error in Lapack routine dsyev.', inf)
      else
         info = inf
      end if

      ! Clean up.

100   if (.not.cev) call Clean(E%V)
      if (.not. usemywork) deallocate(work)
      call CleanTemp(A)

   end subroutine RmSpec

end module RmatSpec_m
