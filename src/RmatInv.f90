module RmatInv_m

use MatranUtil_m
use Rmat_m
use RmatLudpp_m
use RmatChol_m

implicit none

#ifdef OVERVIEW

This module contains a generic function Inv to compute the inverse of
a square Rmat.  It also defines a unary operator .inv. that returns
the inverse of a Rmat.

It should be stressed that the inverse of a matrix is seldom needed in
practice, and that operations like multiplying by an inverse, can
usually be done more stably and efficiently by the programs in the
Solve suite.

Authors: Che Rung Lee, Pete Stewart
Oct  1 2003

#endif

   interface Inv
      module procedure RmInv
   end interface
   
   interface operator(.inv.)
      module procedure RmInv_o
   end interface

contains

   subroutine RmInv(C, A, luda, chola, info, mywork)
      type(Rmat), intent(out) :: C
         ! On return contains the inverse of A

      type(Rmat), intent(in)  :: A
         ! The Rmat whose inverse is to be computed.

      type(RmatLudpp), optional, intent(inout) :: luda
         ! If present, luda%companion is true, 
         ! and A%type = GE, RmInv
         ! uses the lu dcomposition contained in luda to
         ! compute the inverse of A.
         ! If present, luda%companion is false,
         ! and A%type = GE, RmInv computes an
         ! LU decomposition and passes it back
         ! to the user for later use.


      type(RmatChol), optional, intent(inout)  :: chola
         ! If present, luda%companion is true, 
         ! and A%tag = HP, RmInv
         ! uses the Cholesky dcomposition contained in chola to
         ! compute the inverse of A.
         ! If present, luda%companion is false,
         ! and A%tag = HP, RmInv computes a
         ! Cholesky decomposition and passes it back
         ! to the user for later use.

      integer, optional :: info
         ! If present, RmInv returns the info paramater
         ! from any LAPACK routines in lieu of an
         ! error return.

      real(wp), pointer, optional :: mywork(:)
         ! Optional work array for Lapack routines.
         ! If it is present, it will be used, possibly
         ! after a reallocation.  It is not deallocated
         ! by RmEig.

      real(wp), pointer :: work(:) => null()
      real(wp) :: twork(1)

      integer :: ipiv(A%nrow)
      integer :: n, lwork, inf, i, j
      logical :: usemywork = .true.
      

      n = A%nrow
      if (n /= A%ncol)&
           call MatranError('RmInv in RmatInv: Matrix not square.')

      call GuardTemp(A)
      call ReshapeAry(C, n, n)

      ! Compute the inverse according to the flag.

      select case(A%tag)

      ! Lower triangular matrix.  Use the LAPACK routine
      ! DTRTRI to compute the inverse.

      case("LT")
         C%a(1:n, 1:n) = A%a(1:n, 1:n)

#ifdef dbl
         call dtrtri("L", "N", n, C%a, C%narow, inf)
#endif
#ifdef sngl
         call strtri("L", "N", n, C%a, C%narow, inf)
#endif

         if (inf/=0) then ! error
            if (present(info)) then
               info = inf
               go to 200
            else
               call SupportError('RmInv in RmatInv: &
                       &Error in LAPACK routine DTRTRI, info = ', inf)
            end if
         end if

      ! Upper triangular matrix.  Use the LAPACK routine
      ! DTRTRI to compute the inverse.

      case("UT")
         C%a(1:n, 1:n) = A%a(1:n, 1:n)

#ifdef dbl
         call dtrtri("U", "N", n, C%a, C%narow, inf)
#endif
#ifdef sngl
         call strtri("U", "N", n, C%a, C%narow, inf)
#endif

         if (inf/=0) then ! error
            if (present(info)) then
               info = inf
               go to 200
            else
               call SupportError('RmInv in RmatInv: &
                       &Error in LAPACK routine DTRTRI, info = ', inf)
            end if
         end if

      ! General matrix.  Get a LU decomposition and use it
      ! to compute the inverse via the LAPACK routine DGETRI.

      case("GE")

         ! Get an LU decomposition.  This decomposition must
         ! end up in C in the form required by DGETRI

         if(present(luda)) then

            ! If necessry compute a RmatLudpp of A.

            if(.not. luda%companion) call ludpp(luda, A)

            ! Copy to C.

            C%a(1:n,1:n) = luda%U%a(1:n,1:n)
            do j=1,n
               C%a(j+1:n,j) = luda%L%a(j+1:n,j)
            end do
            ipiv =luda%pvt

         else

            ! luda not present.  Compute an LU decomposition
            ! in C.

            C%a(1:n, 1:n) = A%a(1:n, 1:n)

#ifdef dbl
            call dgetrf(n, n, C%a, C%narow, ipiv,inf)
#endif
#ifdef sngl
            call sgetrf(n, n, C%a, C%narow, ipiv,inf)
#endif

            if (inf/=0) then ! error
               if (present(info)) then
                  info = inf
                  go to 200
               else
                  call SupportError('RmInv in RmatInv: &
                       &Error in LAPACK routine DGETRF, info = ', inf)
               end if
            end if
         end if

         ! Get working storage.

#ifdef dbl
         call dgetri(n, C%a, C%narow, ipiv, twork, -1, inf)
#endif
#ifdef sngl
         call sgetri(n, C%a, C%narow, ipiv, twork, -1, inf)
#endif

         if (inf/=0) then ! error
            if (present(info)) then
               info = inf
               go to 200
            else
               call SupportError('RmInv in RmatInv: &
                       &Error in LAPACK routine DGETRI, info = ', inf)
            end if
         end if

         lwork = twork(1)
            
         ! allocate working space
         if (present(mywork)) then
            call ReshapeAry(mywork, lwork)
            work => mywork
         else
            usemywork = .false.
            call ReshapeAry(work, lwork)
         end if

         ! Compute the inverse.

#ifdef dbl
         call dgetri(n, C%a, C%narow, ipiv, work, lwork, inf)
#endif
#ifdef sngl
         call sgetri(n, C%a, C%narow, ipiv, work, lwork, inf)
#endif

         if (inf/=0) then ! erro
            if (present(info)) then
               info = inf
               go to 200
            else
               call SupportError('RmInv in RmatInv: &
                    &Error in LAPACK routine DGETRI, info = ', inf)
            end if
         end if

      ! Symmmetric matrix.  Use DSYTRF to get a symmetric factorization
      ! of A.  Then invert A using

      case("HE")

         ! Compute the symmetric factoriization

         C%a(1:n, 1:n) = A%a(1:n, 1:n)
         
         ! Get working storage.

#ifdef dbl
         call dsytrf("L", n, C%a, C%narow, ipiv, twork, -1, inf)            
#endif
#ifdef sngl
         call ssytrf("L", n, C%a, C%narow, ipiv, twork, -1, inf)            
#endif

         if (inf/=0) then ! error handle
            if (present(info)) then
               info = inf
               go to 200
            else
               call SupportError('RmInv in RmatInv: &
                       &Error in LAPACK routine DSYTRF, info = ', inf)
            end if
         end if

         
         lwork = twork(1)
         if (lwork < n) lwork = n

         if (present(mywork)) then
            call ReshapeAry(mywork, lwork)
            work => mywork
         else
            usemywork = .false.
            call ReshapeAry(work, lwork)
         end if

         ! Compute the factorization

#ifdef dbl
         call dsytrf("L", n, C%a, C%narow, ipiv, work, lwork, inf)
#endif
#ifdef sngl
         call ssytrf("L", n, C%a, C%narow, ipiv, work, lwork, inf)
#endif

         if (inf/=0) then ! error handle
            if (present(info)) then
               info = inf
               go to 200
            else
               call SupportError('RmInv in RmatInv: &
                       &Error in LAPACK routine DSYTRF, info = ', inf)
            end if
         end if

         ! Compute the inverse.

#ifdef dbl
         call dsytri("L", n, C%a, C%narow, ipiv, work, inf)
#endif
#ifdef sngl
         call ssytri("L", n, C%a, C%narow, ipiv, work, inf)
#endif

         if (inf/=0) then ! error handle
            if (present(info)) then
               info = inf
               go to 200
            else
               call SupportError('RmInv in RmatInv: &
                       &Error in LAPACK routine DSYTRI, info = ', inf)
            end if
         end if

        ! copy the result
         do i = 1, n
            do j = i+1, n
               C%a(i,j) = C%a(j,i)
            end do
         end do

      ! Symmetric positive definite matrix.  Get a Cholesky
      ! factorization and use it to compute the inverse.

      case("HP")

         ! Get a Cholesky decompostion.  It must end up
         ! in C.

         if(present(chola)) then

            ! If necessary compute a RmatChol of A.

            if(.not. chola%companion) call chol(chola,A)

            ! Copy to C.

            C%a(1:n,1:n) = chola%R%a(1:n,1:n)

         else

            ! chola not present.  Compute a Cholesky decompositon
            ! in C.

            C%a(1:n, 1:n) = A%a(1:n, 1:n)

#ifdef dbl
            call dpotrf("U", n, C%a, C%narow, inf)
#endif
#ifdef sngl
            call spotrf("U", n, C%a, C%narow, inf)
#endif

            if (inf/=0) then ! error
               if (present(info)) then
                  info = inf
                  go to 200
               else
                  call SupportError('RmInv in RmatInv: &
                       &Error in LAPACK routine DPOTRF, info = ', inf)
               end if
            end if
         end if

         ! Compute the Inverse

#ifdef dbl
         call dpotri("U", n, C%a, C%narow, inf)
#endif
#ifdef sngl
         call spotri("U", n, C%a, C%narow, inf)
#endif

         if (inf/=0) then ! Error
            if (present(info)) then
               info = inf
               go to 200
            else
               call SupportError('RmInv in RmatInv: &
                       &Error in LAPACK routine DPOTRF, info = ', inf)
            end if
         end if

         do i = 1, n
            do j = i+1, n
               C%a(j, i) = C%a(i, j)
            end do
         end do

      case default

         call MatranError('RmInv in RmatInv: Illegal tag.')

      end select
      
      ! Clean up.

200   if(.not. usemywork) deallocate(work)

      call CleanTemp(A)

   end subroutine RmInv
   
   ! C = .inv.A

   function RmInv_o(A) result(C)

      type(Rmat) :: C
      type(Rmat), intent(in) :: A

      call GuardTemp(A)
      C%a => null()
      C%temporary => null()
      call Clean(C)

      call Inv(C, A)
      call SetTemp(C)

      call CleanTemp(A)
   end function RmInv_o

end module RmatInv_m
