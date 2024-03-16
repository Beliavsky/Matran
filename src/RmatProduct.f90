module RmatProduct_m

use MatranUtil_m
use Rmat_m
implicit none

#ifdef OVERVIEW

RmatProduct implements the product of matrices.  As is true of all
Matran operations, the product comes in two forms: an explicit call
to a subroutine and a function implementing a binary or unary
operator.  The actual computation takes place in the former, which is
called by the latter.  They also differ in the sizes of the arrays
concerned.  If the result is an mxn matrix, the subroutine will try to
fit it into the current output array, reshaping it only if necessary.
The functional form always returns the result in an mxn array.

The routines are (here s is a scalar)

   Times
      subroutine RmTimesScalarRm(C, s, A): C = s*A
      subroutine RmTimesRmRm(C, A, B) : C = A*B
      subroutine RmTimesXhy(C, A, B): C = A'*B
      subroutine RmTimesXyh(C, A, B): C = A*B'
      subroutine RmTimesXhx(C, A): C = A'*A
      subroutine RmTimesXxh(C, A): C = A*A'
   *
      function RmTimesScalarRm_o(s, A) result(C): C = s*A
      function RmTimesRmScalar_o(A, s) result(C): C =  A*s
      function RmTimesRmRm_o(A, B) result(C): C = A*B
   .xhy.
      function RmTimes_xhy(A, B) result(C): C = A'*B   
   .xyh.
      function RmTimes_xyh(A, B) result(C): C = A*B'
   .xhx.
      function RmTimes_xhx(A) result(C): C = A'*A
   .xxh.
      function RmTimes_xxh(A) result(C): C = A*A'

In the names of the procedures

   x is the first argument
   y is the second argument
   h is the conjugate transpose operator

Thus xyh stands for x times y conjugate transpose.

In RmTimesRmRm and the corresponding operator *, if one of the
operands is a 1x1 matrix, it is treated as a scalar.

Note: When one of A or B is triangular, RmatProduct attempts to
perform the multiplication using the LAPACK routine dtrmm.  However,
this routine will not support the following combinations, which are
therefore performed by the routine dgemm.

   A*B'   A triangular, B general
   A'*B   A general, B triangular


Author: Pete Stewart
Aug 20 2003

#endif

   real(wp), parameter :: ZERO = 0, ONE = 1

   interface Times
      module procedure RmTimesScalarRm, RmTimesRmRm
   end interface Times

   interface operator(*)
      module procedure RmTimesScalarRm_o, RmTimesRmScalar_o, &
                       RmTimesRmRm_o
   end interface

   interface TimesXhy
      module procedure RmTimesXhy
   end interface TimesXhy

   interface operator(.xhy.)
      module procedure RmTimes_xhy
   end interface

   interface TimesXyh
      module procedure RmTimesXyh
   end interface TimesXyh
   
   interface operator(.xyh.)
      module procedure RmTimes_xyh
   end interface

   interface TimesXhx
      module procedure RmTimesXhx
   end interface

   interface operator(.xhx.)
      module procedure RmTimes_xhx
   end interface

   interface TimesXxh
      module procedure RmTimesXxh
   end interface

   interface operator(.xxh.)
      module procedure RmTimes_xxh
   end interface

contains

   ! Scalar times Rmat

   subroutine RmTimesScalarRm(C, s, A)
      type(Rmat), intent(out) :: C
      real(wp), intent(in) :: s
      type(Rmat), intent(in) :: A

      integer m, n

      m = A%nrow
      n = A%ncol

      call GuardTemp(A)
      call ReshapeAry(C, m, n)

      if (m==0 .or. n==0) then
         C%tag = 'GE'

      else

         if (s>=0 .OR. A%tag/='HP') then
            C%tag = A%tag
         else
            C%tag = 'HE'
         end if
         C%a(1:m,1:n) = s*A%a(1:m,1:n)

      end if

      call CleanTemp(A)

   end subroutine RmTimesScalarRm

   ! s*A

   function RmTimesScalarRm_o(s, A) result(C)
      type(Rmat) :: C
      real(wp), intent(in) :: s
      type(Rmat), intent(in) :: A

      C%a => null()
      C%temporary => null()
      call Clean(C)

      call Times(C, s, A)
      call SetTemp(C)

   end function RmTimesScalarRm_o

   ! A*s

   function RmTimesRmScalar_o(A, s) result(C)
      type(Rmat) :: C
      real(wp), intent(in) :: s
      type(Rmat), intent(in) :: A

      C%a => null()
      C%temporary => null()
      call Clean(C)
      call Times(C, s, A)
      call SetTemp(C)

   end function RmTimesRmScalar_o

   ! Rmat times Rmat

   subroutine RmTimesRmRm(C, A, B)
      type(Rmat), intent(inout) :: C
      type(Rmat), intent(in)  :: A, B

      integer :: k, m, n

      call GuardTemp(A)
      call GuardTemp(B)

      ! Treat 1x1 matrices as scalars


      if (A%nrow==1 .and. A%ncol==1) then
         call RmTimesScalarRm(C, A%a(1,1), B)

      else if (B%nrow==1 .and. B%ncol==1) then
         call RmTimesScalarRm(C, B%a(1,1), A)

      else

         ! Check Dimenstions

         m = A%nrow
         k = A%ncol
         n = B%ncol

         if (k /= B%nrow) then
            call MatranError('RmTimesRmRm in RmatProduct: &
                             &Dimensions incompatible for multiplication.')
         end if

         call ReshapeAry(C, m, n)

         if (m==0 .or. k==0 .or. m==0) then

            ! Null matrix.  Set the tag and get out.

            C%tag = 'GE'

         else

            ! Set the tag.

            if (A%tag=='LT' .and. B%tag=='LT' .or. &
                A%tag=='UT' .and. B%tag=='UT' ) then
               C%tag = A%tag
            else
               C%tag = 'GE'
            end if

            ! Perform the multiplication.

            if (A%tag=='UT' .and. m==k) then

               C%a(1:m,1:n) = B%a(1:m,1:n)
#ifdef dbl               
               call dtrmm('L', 'U', 'N', 'N', m, n, ONE, &
                                            A%a, A%narow, &
                                            C%a, B%narow)
#endif
#ifdef sngl               
               call strmm('L', 'U', 'N', 'N', m, n, ONE, &
                                            A%a, A%narow, &
                                            C%a, B%narow)
#endif

            else if (A%tag == 'LT' .and.  m==k) then

               C%a(1:m,1:n) = B%a(1:m,1:n)
#ifdef dbl
               call dtrmm('L', 'L', 'N', 'N', m, n, ONE, &
                                            A%a, A%narow, &
                                            C%a, B%narow)
#endif
#ifdef sngl
               call strmm('L', 'L', 'N', 'N', m, n, ONE, &
                                            A%a, A%narow, &
                                            C%a, B%narow)
#endif

            else if (B%tag == 'UT' .and. n==k) then

               C%a(1:m,1:n) = A%a(1:m,1:n)
#ifdef dbl
               call dtrmm('R', 'U', 'N', 'N', m, n, ONE, &
                                            B%a, B%narow, &
                                            C%a, A%narow)
#endif
#ifdef sngl
               call strmm('R', 'U', 'N', 'N', m, n, ONE, &
                                            B%a, B%narow, &
                                            C%a, A%narow)
#endif

            else if (B%tag == 'LT' .and. n==k) then

               C%a(1:m,1:n) = A%a(1:m,1:n)
#ifdef dbl
               call dtrmm('R', 'L', 'N', 'N', m, n, ONE, &
                                            B%a, B%narow, &
                                            C%a, A%narow)
#endif
#ifdef sngl
               call strmm('R', 'L', 'N', 'N', m, n, ONE, &
                                            B%a, B%narow, &
                                            C%a, A%narow)
#endif
            else

#ifdef dbl
               call dgemm('N', 'N', m, n, k, ONE, A%a, A%narow, &
                                                   B%a, B%narow, &
                                            ZERO,  C%a, C%narow)
#endif
#ifdef sngl
               call sgemm('N', 'N', m, n, k, ONE, A%a, A%narow, &
                                                   B%a, B%narow, &
                                            ZERO,  C%a, C%narow)
#endif
            end if
         end if
      end if

      call CleanTemp(A)
      call CleanTemp(B)

   end subroutine RmTimesRmRm

   ! C = A*B


   function RmTimesRmRm_o(A, B) result(C)
      type(Rmat) :: C
      type(Rmat), intent(in) :: A, B

      C%a => null()
      C%temporary => null()
      call Clean(C)
      call Times(C, A, B)
      call SetTemp(C)

   end function RmTimesRmRm_o

   ! Conjugate transpose of a Rmat times a Rmat

   subroutine RmTimesXhy(C, A, B)
      type(Rmat), intent(inout) :: C
      type(Rmat), intent(in)  :: A, B

      integer :: k, m, n

      call GuardTemp(A)
      call GuardTemp(B)

      m = A%ncol
      k = A%nrow
      n = B%ncol
      if (k /= B%nrow)&
         call MatranError('RmTimesXhy in RmatProduct: &
                          &Dimensions incompatible for multiplication.')

      call ReshapeAry(C, m, n)

      if (m==0 .or. k==0 .or. m==0) then

         ! Null matrix.  Set the tag and get out.

         C%tag = 'GE'
      else

         if (A%tag=='LT' .and. B%tag=='UT' .or. &
             A%tag=='UT' .and. B%tag=='LT' ) then
            C%tag = B%tag
         else
            C%tag = 'GE'
         end if

         if (A%tag=='UT' .and. m==k) then

            C%a(1:m,1:n) = B%a(1:m,1:n)
#ifdef dbl
            call dtrmm('L', 'U', 'T', 'N', m, n, ONE, &
                                         A%a, A%narow, &
                                         C%a, B%narow)
#endif
#ifdef sngl
            call strmm('L', 'U', 'T', 'N', m, n, ONE, &
                                         A%a, A%narow, &
                                         C%a, B%narow)
#endif

         else if (A%tag == 'LT' .and. m==k) then


            C%a(1:m,1:n) = B%a(1:m,1:n)
#ifdef dbl
            call dtrmm('L', 'L', 'T', 'N', m, n, ONE, &
                                         A%a, A%narow, &
                                         C%a, B%narow)
#endif
#ifdef sngl
            call strmm('L', 'L', 'T', 'N', m, n, ONE, &
                                         A%a, A%narow, &
                                         C%a, B%narow)
#endif

         else
#ifdef dbl
            call dgemm('T', 'N', m, n, k, ONE, A%a, A%narow, &
                                                B%a, B%narow, &
                                          ZERO, C%a, C%narow)
#endif
#ifdef sngl
            call sgemm('T', 'N', m, n, k, ONE, A%a, A%narow, &
                                                B%a, B%narow, &
                                          ZERO, C%a, C%narow)
#endif
         end if
      end if

      call CleanTemp(A)
      call CleanTemp(B)

   end subroutine RmTimesXhy

   ! A.xhy.B

   function RmTimes_xhy(A, B) result(C)
      type(Rmat) :: C
      type(Rmat), intent(in) :: A, B

      C%a => null()
      C%temporary => null()
      call Clean(C)
      call TimesXhy(C, A, B)
      call SetTemp(C)

   end function RmTimes_xhy


   ! Rmat times transpose of a Rmat

   subroutine RmTimesXyh(C, A, B)
      type(Rmat), intent(inout) :: C
      type(Rmat), intent(in)  :: A, B

      integer :: k, m, n


      m = A%nrow
      k = A%ncol
      n = B%nrow
      if (k /= B%ncol) then
         call MatranError('RmTimesXyh in RmatProduct: &
                          &Dimensions incompatible for multiplication.')
      end if

      call GuardTemp(A)
      call GuardTemp(B)

      call ReshapeAry(C, m, n)

      if (k==0 .or. m==0 .or. n==0) then

         ! Null matrix.  Set the tag and get out.

         C%tag = 'GE'
      else
         if (A%tag=='LT' .and. B%tag=='UT' .or. &
             A%tag=='UT' .and. B%tag=='Lt' ) then
            C%tag = A%tag
         else
            C%tag = 'GE'
         end if

         if (B%tag == 'UT' .and. n==k) then

            C%a(1:m,1:n) = A%a(1:m,1:n)
#ifdef dbl
            call dtrmm('R', 'U', 'T', 'N', m, n, ONE, &
                                         B%a, B%narow, &
                                         C%a, A%narow)
#endif
#ifdef sngl
            call strmm('R', 'U', 'T', 'N', m, n, ONE, &
                                         B%a, B%narow, &
                                         C%a, A%narow)
#endif

         else if (B%tag == 'LT' .and. n==k) then
            C%a(1:m,1:n) = A%a(1:m,1:n)
#ifdef dbl
            call dtrmm('R', 'L', 'T', 'N', m, n, ONE, &
                                         B%a, B%narow, &
                                         C%a, A%narow)
#endif
#ifdef sngl
            call strmm('R', 'L', 'T', 'N', m, n, ONE, &
                                         B%a, B%narow, &
                                         C%a, A%narow)
#endif
         else

#ifdef dbl
            call dgemm('N', 'T', m, n, k, ONE, A%a, A%narow, B%a, B%narow, &
                                                       ZERO,  C%a, C%narow)
#endif
#ifdef sngl
            call sgemm('N', 'T', m, n, k, ONE, A%a, A%narow, B%a, B%narow, &
                                                       ZERO,  C%a, C%narow)
#endif
         end if

      end if

      call CleanTemp(A)
      call CleanTemp(B)

   end subroutine RmTimesXyh

   ! A.xhy.B

   function RmTimes_xyh(A, B) result(C)
      type(Rmat) :: C
      type(Rmat), intent(in) :: A, B

      C%a => null()
      C%temporary => null()
      call Clean(C)
      call TimesXyh(C, A, B)
      call SetTemp(C)

   end function RmTimes_xyh

   ! Cross product matrix of a Rmat

   subroutine RmTimesXhx(C, A)
      type(Rmat), intent(inout) :: C
      type(Rmat), intent(in)  :: A

      integer :: i, j, k, n

      n = A%ncol
      k = A%nrow

      call GuardTemp(A)
      call ReshapeAry(C, n, n)

      if (n==0 .or. k==0) then

         ! Null matrix.  Set the tag and get out.

         C%tag = 'GE'

      else

         C%tag = 'HP'

#ifdef dbl
         call dsyrk('U', 'T',  n, k, ONE, A%a, A%narow, &
                                  ZERO, C%a, C%narow)
#endif
#ifdef sngl
         call ssyrk('U', 'T',  n, k, ONE, A%a, A%narow, &
                                  ZERO, C%a, C%narow)
#endif
         do i=2, n
            do j=1,i-1
               C%a(i,j) = C%a(j,i)
            end do
         end do

      end if

      call CleanTemp(A)

   end subroutine RmTimesXhx

   ! .xhx.A

   function RmTimes_xhx(A) result(C)
      type(Rmat) :: C
      type(Rmat), intent(in) :: A

      C%a => null()
      C%temporary => null()
      call Clean(C)
      call TimesXhx(C, A)
      call Settemp(C)
   end function RmTimes_xhx


   ! Exterior product of a Rmat

   subroutine RmTimesXxh(C, A)
      type(Rmat), intent(inout) :: C
      type(Rmat), intent(in)  :: A

      integer :: i, j, k, n

      n = A%nrow
      k = A%ncol

      call GuardTemp(A)
      call ReshapeAry(C, n, n)

      if (n==0 .or. k==0) then

         ! Null matrix.  Set the tag and get out.

         C%tag = 'GE'

      else

         C%tag = 'HP'
#ifdef dbl
         call dsyrk('U', 'N',  n, k, ONE, A%a, A%narow, &
                                  ZERO, C%a, C%narow)
#endif
#ifdef sngl
         call ssyrk('U', 'N',  n, k, ONE, A%a, A%narow, &
                                  ZERO, C%a, C%narow)
#endif
         do i=2, n
            do j=1,i-1
               C%a(i,j) = C%a(j,i)
            end do
         end do

      end if

      call CleanTemp(A)

   end subroutine RmTimesXxh

   ! .xxh.A

   function RmTimes_xxh(A) result(C)
      type(Rmat) :: C
      type(Rmat), intent(in) :: A

      C%a => null()
      C%temporary => null()
      call Clean(C)
      call TimesXxh(C, A)
      call SetTemp(C)

   end function RmTimes_xxh


end module RmatProduct_m
