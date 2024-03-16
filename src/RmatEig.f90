module RmatEig_m

use MatranUtil_m
use Rmat_m
implicit none

#ifdef OVERVIEW

Let A be a nondefective matrix of order n.  Then there is a nonsingular
X such that if Y` is the inverse of X,

   Y`AX = D = diag(d1, ... , dn)    (*)

(here Y` denotes the conjugate transpose of Y).  The scalars di are
the eigenvalues of A and the corresponding columns of X and Y are the
right and left eigenvectors of A.  These quantities are in general
complex, even when A is real.  The decomposition (*) is called the
eigendecomposition of A.

This module provides a generic subroutine Eig to compute the
eigendecomposition of a Rmat.  The decomposition is contained in a
derived type RmatEig defined below.

Authors: CheRung Lee and Pete Stewart
Jun 17 2003

#endif

   type RmatEig
      integer :: n = 0                       ! The order of the decomposition
      complex(wp), pointer :: D(:)=>null()   ! The eigenvalues
      complex(wp), pointer :: X(:,:)=>null() ! The right eigenvectors
      complex(wp), pointer :: Y(:,:)=>null() ! The left eigenvectors
      logical              :: companion &    ! True if the decomposition is
                               = .false.     ! associated with a Rmat of
                                             ! interest
   end type RmatEig

   interface Eig
      module procedure RmEig
   end interface

   interface Clean
      module procedure RmCleanEig
   end interface

contains

   ! Clean for type RmatEig

   subroutine RmCleanEig(s)
      type(RmatEig), intent(inout) :: s

      if (associated(s%D)) deallocate(s%D)
      if (associated(s%X)) deallocate(s%X)
      if (associated(s%Y)) deallocate(s%Y)

      s%companion = .false.
   end subroutine RmCleanEig



   subroutine RmEig(E, A, rvecs, lvecs, info, wrv, wlv, mywork)

      type(RmatEig), intent(out) :: E
         ! Contains the eigendecomposition of A.

      type(Rmat), intent(in)     :: A
         ! The Rmat whose eigendecomposition is to be computed.

      logical, optional, intent(in) :: rvecs
         ! If present and true, compute right eigenvectors.
         ! Otherwise do not compute right eigenvectors.

      logical, optional, intent(in) :: lvecs
         ! If present and true, compute left eigenvectors.
         ! Otherwise do not compute left eigenvectors.

      integer, optional, intent(out) :: info
         ! If present, return dgeev's info parameter in lieu of
         ! an error return.

      real(wp), pointer, optional :: wrv(:,:), wlv(:,:), mywork(:)
         ! Optional work arrays for dgeev.  If any of
         ! these is present, they will be used, possibly after
         ! a reallocation.  They are not deallocated by RmEig.

      integer :: lwork, inf, j, ldrv, ldlv, n

      logical :: crv, clv, usewrv,usewlv, usemywork

      character(1) :: joblv, jobrv
      
      real(wp) :: wr(A%nrow), wi(A%nrow)
      real(wp) :: Aa(A%nrow, A%ncol)

      real(wp), pointer :: rv(:,:)=>null(), lv(:,:)=>null(), &
                work(:)=>null()
      real(wp) :: twork(1)


      call GuardTemp(A)

      ! Check for square matrix.

      n = A%nrow
      if (n /= A%ncol)&
         call MatranError('RmEig in RmatEig: Matrix not square.')

      call ReshapeAry(E%D, n)

      ! Check to see if right eigenvectors are wanted, and
      ! if so allocate temporary storage for the real representation.

      crv = .false.
      jobrv = 'N'
      ldrv = 1
      usewrv = .false.
      if (present(rvecs)) then
         if (rvecs) then
            crv = .true.
            jobrv = 'V'
            if (present(wrv)) then
               call ReshapeAry(wrv, n, n)
               rv => wrv
               usewrv = .true.
            else
               usewrv = .false.
               call ReshapeAry(rv, n, n)
         end if
               ldrv = size(rv, 1)
         end if
      end if

      ! Check to see if left eigenvectors are wanted, and
      ! if so allocate temporary storage for the real representation.

      clv = .false.
      joblv = 'N'
      ldlv = 1
      usewlv = .false.
      if (present(lvecs)) then
         if (lvecs) then
            clv = .true.
            joblv = 'V'
            if (present(wlv)) then
               call ReshapeAry(wlv, n, n)
               lv => wlv
               usewlv = .true.
            else
               usewlv = .false.
               call ReshapeAry(lv, n, n)
            end if
            ldlv = size(lv, 1)
         end if
      end if

      ! Copy A%a to Aa
      
      Aa = A%a(1:n, 1:n)

      ! Set up the work array.

#ifdef dbl
      call dgeev(joblv, jobrv, n, Aa, n, wr, wi, &
                    lv, ldlv, rv, ldrv, twork, -1, inf)
#endif
#ifdef sngl
      call sgeev(joblv, jobrv, n, Aa, n, wr, wi, &
                    lv, ldlv, rv, ldrv, twork, -1, inf)
#endif

      if (inf /= 0) go to 50

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
      call dgeev(joblv, jobrv, n, Aa, n, wr, wi, &
                 lv, ldlv, rv, ldrv, work, lwork, inf)
#endif
#ifdef sngl
      call sgeev(joblv, jobrv, n, Aa, n, wr, wi, &
                 lv, ldlv, rv, ldrv, work, lwork, inf)
#endif

      if (inf /= 0) go to 50

      ! Move the results into the RmatEig.

      if (crv) call ReshapeAry(E%X,n,n)
      if (clv) call ReshapeAry(E%Y,n,n)

      j = 1; 
      do while (j <= n)
        
         E%D(j) = cmplx(wr(j), wi(j), wp)

         if (wi(j) == 0) then  ! real eig val

            if (crv) then
               E%X(:,j) = cmplx(rv(:,j), 0.0D0, wp)
            end if
            
            if (clv) then
               E%Y(:,j) = cmplx(lv(:,j), 0.0D0, wp)
            end if

            j = j + 1
         else                ! complex eigenvalue

            if (crv) then
                E%X(:,j)   = cmplx(rv(:,j),  rv(:,j+1), wp)
                E%X(:,j+1) = cmplx(rv(:,j), -rv(:,j+1), wp) 
            end if

            if (clv) then
                E%Y(:,j) = cmplx(lv(:,j), lv(:,j+1), wp) 
                E%Y(:,j+1) = cmplx(lv(:,j), -lv(:,j+1), wp)        
            end if

            E%D(j+1) = cmplx(wr(j), -wi(j), wp)
            j = j + 2
         end if
      end do

      E%companion = .true.
      go to 100

      ! error handling
50    if (.not. present(info)) then
         call SupportError('RmEig in RmatEig: &
              &Error in Lapack routine dgeev.',inf)
      else
         info = inf
      end if


      ! Clean up.
100   call CleanTemp(A)
      if (.not.usemywork) deallocate(work)
      if (crv .and. .not.usewrv) deallocate(rv)
      if (clv .and. .not.usewlv) deallocate(lv)

   end subroutine RmEig

end module RmatEig_m
