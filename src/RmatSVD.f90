module RmatSVD_m

use MatranUtil_m
use Rmat_m
use Rdiag_m
use RmatTranspose_m

implicit none

#ifdef OVERVIEW

Let A be an mxn matrix with m>=n.  Then A can be written in the
form

       U|D|V'
   A =  |0|       (*)

where U and V are mxm and nxn orthogonal matrices and D is a diagonal
matrix of order n with diagonal elements

   d1 >= d2 >= ... >= dn >= 0.

The columns of U and V are called the left and right singular vectors
of A, and d1, ..., dn are called the singular values.  If we partition
U = |U1 U2|, where U1 has m columns, we can write

   A = U1 D V'.   (**)

If $m<n$, the decomposition has the form

   A = U|D 0|V'.  (***),

were D is now of order m.  If we partition V = |V1 V2|, where V1 has m
columns, we can write

   A = U D V1'    (****).

The above factoraizations are generically termed singular value
decompositions (SVD) of A.  The SVDs (*) and (***) are called full
decompositions; (**) and (****) are called economy sized
decompositions.

RmatSVD_m uses the lapack routine DGESVD to compute one of the above
SVDs of a Rmat.  The result is returned in defined type RmatSVD defined
below.  The generic subroutine for computing the SVD is

   SVD(svdc, A, wantu, wantv, full, info, mywork)

whose arguments are described belos.

Authors: Che Rung Lee and Pete Stewart
Jun 17 2003

#endif

   type RmatSVD
      type(Rdiag)       :: D         ! The singular values
      type(Rmat)        :: U         ! The left singular vectors
      type(Rmat)        :: V         ! The right singular vectors
      logical           :: companion ! True if the  decomposition is of
                                     ! a Rmat of interest.
   end type RmatSVD

   interface SVD
      module procedure RmSVD
   end interface

contains

   subroutine RmSVD(svdc, A, wantu, wantv, full, info, mywork)

      type(RmatSVD), intent(out), target :: svdc
      ! The singular value decomposition of A.

      type(Rmat), intent(in) :: A
      ! The Rmat whose SVD is to be computed

      logical, optional, intent(in) :: wantu
      ! if present and true, compute left singular vectors.

      logical, optional, intent(in) :: wantv
      ! if present and true, compute right singular vectors.

      logical, intent(in), optional :: full
      ! If present and true compute the full SVD.  Otherwise
      ! compute the economy size SVD.

      integer, optional, intent(out) :: info
      ! If present, return dgesvd's info parameter in lieu of
      ! an error return.

      real(wp), pointer, optional :: mywork(:)
      ! Optional work arrays for dgesvd.  If it
      ! is present, it will be used, possibly after
      ! a reallocation.  It is not deallocated by RmEig.

      ! Internal variables.

      integer :: info_svd, ldu, ldv, lwork, m, min_mn, n, lda

      character(1) :: jobu, jobv

      real(wp) :: temp(A%nrow, A%ncol), twork(1)

      real(wp), pointer :: work(:)=>null()

      logical :: usemywork, wantub, wantvb, fullb


      call GuardTemp(A)

      ! convenient abreviations.

      m = A%nrow
      n = A%ncol
      min_mn = min(m,n)

      ! Transfer A to a working array.

      temp = A%a
      lda = A%narow
      call ReshapeAry(svdc%D, min_mn)

      ! Determine what to compute and how (see the calling sequence
      ! for DGESVD for the meaning of jobu and jobv).

      jobu = 'N'
      jobv = 'N'
      ldv = 1
      ldu = 1

      wantub = .false.
      wantvb = .false.
      fullb  = .false.
      if (present(wantu)) wantub = wantu
      if (present(full))  fullb = full
      if (present(wantv)) wantvb = wantv

      if (m >= n) then

         ! Set up options for U.

         if (wantub) then
            if (fullb) then
               jobu = 'A'
               call ReshapeAry(svdc%U, m, m)
               ldu = svdc%U%narow
            else
               jobu = 'S'
               call ReshapeAry(svdc%U, m, n)
               ldu = svdc%U%narow
            end if
         end if
         
         ! Set up options for V.

         if (wantvb) then
            jobv = 'A'
            call ReshapeAry(svdc%v, n, n)
            ldv = svdc%v%nrow
         end if

      else  ! m < n

         ! Set up options for V.

         if (wantvb) then
            if (fullb) then
               jobv = 'A'
               call ReshapeAry(svdc%V, n, n)
               ldv = svdc%V%narow
            else
               jobv = 'O' ! write into temp then do transpose
               call ReshapeAry(svdc%V, n, m)
            end if
         end if

         ! Set up options for U.

         if (wantub) then
            jobu = 'A'
            call ReshapeAry(svdc%U, m, m)
            ldu = svdc%U%narow
         end if
      end if


      ! Set up work array.

#ifdef dbl
      call dgesvd(jobu, jobv, m, n, temp, lda, svdc%D%a, &
           svdc%U%a, ldu, svdc%V%a, ldv, twork, -1, info_svd)
#endif
#ifdef sngl
      call sgesvd(jobu, jobv, m, n, temp, lda, svdc%D%a, &
           svdc%U%a, ldu, svdc%V%a, ldv, twork, -1, info_svd)
#endif

      if (info_svd /= 0) go to 50        ! error handling
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
      call dgesvd(jobu, jobv, m, n, temp, lda, svdc%D%a, &
           svdc%U%a, ldu, svdc%V%a, ldv, work, lwork, info_svd)
#endif
#ifdef sngl
      call sgesvd(jobu, jobv, m, n, temp, lda, svdc%D%a, &
           svdc%U%a, ldu, svdc%V%a, ldv, work, lwork, info_svd)
#endif

      if (info_svd /= 0) go to 50        ! Error haldling.


      ! DGESVD returns the transpose of V.  If V has been computed
      ! transpose it.

      if (wantvb) then
         if (m>=n) then
            svdc%V = .ctp.svdc%V
         else 
            if (fullb) then
               svdc%V = .ctp.svdc%V
            else
               svdc%V%a(1:n, 1:m) = transpose(temp(1:m,1:n))
            end if
         end if
      end if

      svdc%companion = .true.
      go to 100

50    if (.not. present(info)) then
         call SupportError('RmSVD in RmatSVD: &
              &Error in Lapack routine dgesvd.',info_svd)
      else
         info = info_svd
      end if

      ! Clean up.

100   if(.not. usemywork) deallocate(work)
      call CleanTemp(A)

   end subroutine RmSVD

end module RmatSVD_m
