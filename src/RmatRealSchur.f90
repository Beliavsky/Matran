module RmatRealSchur_m

use MatranUtil_m
use Rmat_m
implicit none

#ifdef OVERVIEW

Let A be a matrix of order n.  Then there is an orthogonal matrix U
such that

   A = UTU'

where T is block upper triangular with 1x1 or 2x2 blocks on its
diagonal.  The 1x1 blocks are the real eigenvalues of A.  The 2x2
blocks contain the complex eigenvalues of A.  Such a decomposition is
called a real Schur decompostion of A.  The 2x2 blocks can be
standardized to have the form

   r b
   c r

where bc < 0.  The Scalar r is the common real part of the eigenvalues
of this block.  The imaginary parts are sqrt(b)*sqrt(c) and its
negative.

The real Schur form of a Rmat A is computed by the generic subroutine

    Schur(S, A, wantu, info, mywork)

whose arguments are described below.  It returns the real Schur
form in a defined type RmatRealSchur, which is defined below.

The order in which eigenvalues appear on the diagonal of T cannot be
predicted.  Thus it may be necessary to reorder the blocks.  The
subroutine SchurReorder.  moves diagonal a block up or down
the diagonal of $T$ by pairwise exchanges.  Its calling sequence is

    ReorderSchur(S, i1, i2, info).

For further information, see below.

Authors: Che Rung Lee, Pete Stewart
Jun 17 2003


#endif

   type RmatRealSchur
      type(Rmat) :: T                ! The block upper triangular matrix
                                     ! of the decomposition.
      type(Rmat) :: U                ! The orthogonal matrix of the
                                     ! decomposition.
      complex(wp), pointer :: D(:)   ! D containes the eigenvalues of T
                                     ! in the order the appear on the
                                     ! diagonal of T.
      logical           :: companion ! True if the decomposition is
                                     ! associated with a Rmat of
                                     ! interest.
   end type RmatRealSchur

   interface Schur
      module procedure RmRealSchur
   end interface

   interface ReorderSchur
      module procedure RmReorderSchur
   end interface

   interface Clean
      module procedure RmCleanRealSchur
   end interface

contains

   ! Clean for type RmatRealSchur.

   subroutine RmCleanRealSchur(s)
      type(RmatRealSchur), intent(inout) :: s

      if (associated(s%D)) deallocate(s%D)
      call clean(s%U)
      call clean(s%T)

      s%companion = .false.
   end subroutine RmCleanRealSchur



   subroutine RmRealSchur(S, A, wantu, info, mywork)

      type(RmatRealSchur), intent(out) :: S
         ! The real Schur decomposition of A.

      type(Rmat), intent(in)     :: A
         ! The Rmat whose real Schur decomposition is to be computed.

      logical, optional:: wantu
         ! If present and true, compute the orthogonal part
         ! of the decomposition.
      
      integer, optional, intent(out) :: info
         ! If present, RmRealSchur returns the info parameter
         ! from DGEES in liew of an error message.

      real(wp), optional, pointer:: mywork(:)
         ! Optional work array for DGEES.  If it is present,
         ! it will be used, possibly after a reallocation.
         ! It is not deallocated by RmRealSchur.

      real(wp), pointer :: work(:)=>null()
      real(wp), target :: twork(1), wr(A%nrow), wi(A%nrow)
      integer  :: lwork, inf, n, sdim
      logical  :: usemywork, dum
      character(1) :: jobu


      call GuardTemp(A)

      ! Check for square matrix.

      n = A%nrow
      if (n /= A%ncol)&
         call MatranError('RmRealSchur in RmatRealSchur:&
                           & Matrix is not square.')
      
      ! Initialize.

      call ReshapeAry(S%T, n, n)
      call ReshapeAry(S%D, n)
      S%T = A

      jobu = 'N'
      if (present(wantu)) then
         if (wantu) then 
            jobu = 'V'
            call ReshapeAry(S%U, n, n)
         end if
      end if

      ! Allocate memory

#ifdef dbl
      call dgees(jobu,'N',dum, n, S%T%a, S%T%narow, sdim,&
           wr, wi, S%U%a, S%U%narow, twork,-1, dum, inf)
#endif
#ifdef sngl
      call sgees(jobu,'N',dum, n, S%T%a, S%T%narow, sdim,&
           wr, wi, S%U%a, S%U%narow, twork,-1, dum, inf)
#endif

      if (inf /= 0) go to 50  ! Error.

      lwork = twork(1)

      if (present(mywork)) then
         usemywork = .true.
         call ReshapeAry(mywork, lwork)
         work => mywork
      else
         usemywork = .false.
         call ReshapeAry(work, lwork)
      end if


      ! compute Schur decomposition

#ifdef dbl
      call dgees(jobu,'N',dum, n, S%T%a, S%T%narow, sdim, &
           wr, wi, S%U%a, S%U%narow, work, lwork, dum, inf)
#endif
#ifdef sngl
      call sgees(jobu,'N',dum, n, S%T%a, S%T%narow, sdim, &
           wr, wi, S%U%a, S%U%narow, work, lwork, dum, inf)
#endif

      if (inf /= 0) go to 50   ! Error

      ! compute the eigenvalues

      S%D = cmplx(wr,wi,wp)

      S%companion = .true.

      go to 100  ! Clean up.
      
      ! Error Handler.

50    if (.not. present(info)) then
         call SupportError('RmRealSchur in RmatRealSchur: &
              &Error in Lapack routine dgees.', inf)
      else
         info = inf 
      end if

      ! Clean up.

100   call CleanTemp(A)
      if (.not. usemywork) deallocate(work)

   end subroutine RmRealSchur

!  RmReorderSchur reorders the blocks in a real Schur decomposition.

   subroutine RmReorderSchur(S, i1, i2, info)
      type(RmatRealSchur), intent(inout) :: S
!         The Schur decomposition whose blocks are to be
!         reordered.

      integer, intent(inout):: i1, i2
!        i1 is the row containing the block to be moved.
!        i2 is row where it ends up.  Since the row of a
!        complex block is ambiguous, i1 and i2 may change
!        plus or minus one on return.

    integer, optional, intent(out):: info
!        If present return the info parameter from dtrex in
!        lieu of an error return.

      integer :: inf, i, iu
      character(1) :: compQ
      real(wp) :: work(S%T%nrow), wi, wr

      ! Initialize.

      compQ = 'N'
      if(associated(S%U%a)) compQ = 'V'

      ! Reorder the Schur form.

#ifdef dbl
      call dtrexc(compQ, S%T%nrow, S%T%a, S%T%narow, &
                  S%U%a, S%U%narow, i1, i2, work, inf)
#endif
#ifdef sngl
      call strexc(compQ, S%T%nrow, S%T%a, S%T%narow, &
                  S%U%a, S%U%narow, i1, i2, work, inf)
#endif


      ! Error handler.

      if(inf /= 0) then
         if (present(info)) then
            info = inf 
         else
            if(inf == 1) &
                 call SupportError('RmRSchurReorder in RmatRealSchur: &
                 &eigenvalues are too closed.', inf)
            if(inf<0) &
                 call SupportError('RmRSchurReorder in RmatRealSchur: &
                 &Error in Lapack routine dtrexc.', inf)
         end if
      end if

      ! recompute the eigenvalues between i1 and i2

      if (i1 <= i2) then
         i = i1; iu = i2
      else
         i = i2; iu = i1
      end if
      do while (i<=iu)
         wr = S%T%a(i,i)
         if (i /= S%T%nrow) then
            if (S%T%a(i+1,i)/=0) then
               wi = sqrt(abs(S%T%a(i,i+1)))*sqrt(abs(S%T%a(i+1,i)))
               S%D(i) = cmplx(wr, wi)
               s%D(i+1) = cmplx(wr,-wi)
               i = i+2
            else
               S%D(i) = cmplx(wr, 0)
               i = i+1
            end if
         else
            S%D(i) = cmplx(wr, 0)
            i = i+1
         end if
      end do

   end subroutine RmReorderSchur

end module RmatRealSchur_m
