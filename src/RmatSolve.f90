module RmatSolve_m

use MatranUtil_m
use Rmat_m
use RmatLudpp_m
use RmatChol_m
use RmatPivot_m

implicit none

#ifdef OVERVIEW

RmatSolve implements the products of a matrix with its inverse.  As is
true of all Matran operations, the product comes in two forms: an
explicit call to a subroutine and a function implementing a binary or
unary operator.  The actual computation takes place in the former,
which is called by the latter.  They also differ in the sizes of the
arrays concerned.  if the result is an mxn matrix, the subroutine will
try to fit it into the current array C, reshaping it only if
necessary.  The functional form always returns the result in an mxn
array.


The routines are

   Solve
      subroutine RmSolveDiv(C, A, s): C = A/s
      subroutine RmSolveXiy(C, A, B): C = A^-1*B
      subroutine RmSolveXihy(C, A, B): C = A^-H*B
      subroutine RmSolveXyi(C, A, B): C = A*B^-1
      subroutine RmSolveXyih(C, A, B) C = A*B^-h

   /
      functionRmSolve_div(A,S): C = A/s
   .xiy.
      function RmSolve_xiy(A, B) result(C): C = A^-1*B
   .xihy.
      function RmSolve_xihy(A, B) result(C): C = A^-H*B
   .xyi.
      function RmSolve_xyi(A, B) result(C): C = A*B^-1
   .xyih.
      function RmSolve_xyih(A, B) result(C): C = A*B^-H

In the names of the procedures

   x is the first argument
   y is the second argument
   i is the inverse operator
   h is the conjugate transpose operator

Thus xyih stands for x times y inverse conjugate transpose.

The structure of the subroutines is essentially the same.

   Check the dimensions of the matrices.

   Branch on the type of A (LT, UT, GE, HE, HP).

      For GE, HE, and HP obtain an appropriate deomposition.

      Check for singularity.

      Compute C.

      deallocate storage.

      Set the type of C.

   return

Since the subroutines differ only in minor respects (most particularly
in the routines called to compute C), detailed comments are provided
only for RmSolveXiy.

Author: Pete Stewart
Oct  1 2003

#endif


real(wp), parameter :: ONE = 1

   interface SolveDiv
      module procedure RmSolveDiv
   end interface

   interface operator(/)
      module procedure RmSolve_div
   end interface

   interface SolveXiy
      module procedure RmSolveXiy
   end interface

   interface operator(.xiy.)
      module procedure RmSolve_xiy
   end interface

   interface SolveXihy
      module procedure RmSolveXihy
   end interface

   interface operator(.xihy.)
      module procedure RmSolve_xihy
   end interface

   interface SolveXyi
      module procedure RmSolveXyi
   end interface

   interface operator(.xyi.)
      module procedure RmSolve_xyi
   end interface

   interface SolveXyih
      module procedure RmSolveXyih
   end interface

   interface operator(.xyih.)
      module procedure RmSolve_xyih
   end interface



contains

   !  Computes C = A/s.

   subroutine RmSolveDiv(C, A, s)
      type(Rmat), intent(out) :: C
      type(Rmat), intent(in) :: A
      real(wp), intent(in) :: s

      real(wp) :: t

      if (s .eq. 0) &
         call MatranError('RmSolveDiv in RmatSolve:&
                            & Attempt to divide by zero.')
      call GuardTemp(A)

      t = 1/s
      C = A
      C%a = t*C%a

      C%tag = A%tag
      if (s < 0 .and. A%tag=='PO')  C%tag = 'HE'

      call CleanTemp(A)

   end subroutine RmSolveDiv

   function RmSolve_div(A, s) result(C)
      type(Rmat) :: C
      type(Rmat), intent(in) :: A
      real(wp), intent(in) :: s

      C%a => null()
      C%temporary => null()
      call Clean(C)
      call SolveDiv(C, A, s)
      call SetTemp(C)

   end function RmSolve_div

   !  Computes C = A^-1*B for various types of Rmats

   subroutine RmSolveXiy(C, A, B, lu, chl)
      type(Rmat), intent(inout) :: C
      type(Rmat), intent(in) :: A
      type(Rmat), intent(in) :: B
      type(RmatLudpp), intent(inout), optional, target :: lu
      type(RmatChol), intent(inout), optional, target :: chl

      integer :: i
      type(RmatLudpp), target :: lua
      type(RmatLudpp), pointer :: lud
      type(RmatChol), target :: chla
      type(Rmat), pointer :: R

      call GuardTemp(A)
      call GuardTemp(B)

      ! Check Dimensions.

      if (A%nrow==0 .or. A%ncol==0)&
         call MatranError('RmSolveXiy in RmatSolve:&
                            & Cannot invert a null matrix.')

      if (A%nrow/=A%ncol .or. A%nrow/=B%nrow)&
         call MatranError('RmSolveXiy in RmatSolve:&
                            & Inconsistent dimensions.')

      ! Actual calculations are between C and A.

      C = B

      ! Branch on the type of A.

      select case(A%tag)

      case('LT')  ! A is lower triangular

         ! Check for singularity.

         do i=1, A%nrow
            if (A%a(i,i) == 0)&
               call MatranError('RmSolveXiy in RmatSolve:&
                            & Lower triangular A is singular.')
         end do

         ! Compute the product.

#ifdef dbl
         call dtrsm('L', 'L', 'N', 'N', C%nrow, C%ncol, ONE,&
                     A%a, A%narow, C%a, C%narow)
#endif
#ifdef sngl
         call strsm('L', 'L', 'N', 'N', C%nrow, C%ncol, ONE,&
                     A%a, A%narow, C%a, C%narow)
#endif

         ! Set type of C.

         if (B%tag == 'LT') then
            C%tag = 'LT'
         else
            C%tag = 'GE'
         end if


      case('UT') ! A is upper triangular.

         ! Check for singularity.

         do i=1, A%nrow
            if (A%a(i,i) == 0)&
               call MatranError('RmSolveXiy in RmatSolve:&
                            & Upper triangular A is singular.')
         end do

         ! Compute the product.

#ifdef dbl
         call dtrsm('L', 'U', 'N', 'N', C%nrow, C%ncol, ONE,&
                     A%a, A%narow, C%a, C%narow)
#endif
#ifdef sngl
         call strsm('L', 'U', 'N', 'N', C%nrow, C%ncol, ONE,&
                     A%a, A%narow, C%a, C%narow)
#endif

         ! Set type of C.

         if (B%tag == 'UT') then
            C%tag = 'UT'
         else
            C%tag = 'GE'
         end if


      case('GE', 'HE') ! A is a general or a symmetric matrix.

         ! To save computatins when a system with the same
         ! matrix must be repeatedly solved, the user is given
         ! the option of furnishing a pivoted LU decomposition.
         ! If the decomposition is initialized, it is used
         ! in the computation.  If it is not, and Ludpp is
         ! computed and returned to the caller with C.  If
         ! no LU decomposition is present, one is computed
         ! and destroyed on return.  

         if (present(lu)) then
            if (.not.lu%companion) call Ludpp(lu, A)
            lud => lu
         else
            call Ludpp(lua, A)
            lud => lua
         end if

         ! At this point lud points to an Ludpp of A.

         ! Check for singularity of U.

         do i=1, lud%U%nrow
            if (lud%U%a(i,i) == 0)&
               call MatranError('RmSolveXiy in RmatSolve:&
                            & U-factor of A is singular.')
         end do

         ! Compute the product.

         call PivotRow(C, lud%pvt, lud%npvt)
#ifdef dbl
         call dtrsm('L', 'L', 'N', 'N', C%nrow, C%ncol, ONE, &
                     lud%L%a, lud%L%narow, C%a, C%narow)
         call dtrsm('L', 'U', 'N', 'N', C%nrow, C%ncol, ONE, &
                     lud%U%a, lud%U%narow, C%a, C%narow)
#endif
#ifdef sngl
         call strsm('L', 'L', 'N', 'N', C%nrow, C%ncol, ONE, &
                     lud%L%a, lud%L%narow, C%a, C%narow)
         call strsm('L', 'U', 'N', 'N', C%nrow, C%ncol, ONE, &
                     lud%U%a, lud%U%narow, C%a, C%narow)
#endif

         ! Set the type of C.

         C%tag = 'GE'

         ! If necessary, deallocate.

         if (.not.present(lu)) then
            call Clean(lud)
         end if
  

      case('HP') ! Positive definite matrix

         ! To save computatins when a system with the same
         ! matrix must be repeatedly solved, the user is given
         ! the option of furnishing a Cholesky decomposition chl.
         ! If chl%companion is .true., it is used.
         ! in the computation.  If it is not, an Cholesky factor is
         ! computed and returned to the caller with C.  If
         ! no Cholesky factor is present, one is computed
         ! and destroyed on return.  

         if (present(chl)) then
            if (.not.chl%companion) call Chol(chl, A)
            R => chl%R
         else
            call Chol(chla, A)
            R => chla%R
         end if

         ! Check for singularity of R.

         do i=1, R%nrow
            if (R%a(i,i) == 0)&
               call MatranError('RmSolveXiy in RmatSolve:&
                            & Cholesky-factor of A is singular.')
         end do


         ! Solve the system.

#ifdef dbl
         call dtrsm('L', 'U', 'T', 'N', C%nrow, C%ncol, ONE, &
                     R%a, R%narow, C%a, C%narow)
         call dtrsm('L', 'U', 'N', 'N', C%nrow, C%ncol, ONE, &
                     R%a, R%narow, C%a, C%narow)
#endif
#ifdef sngl
         call strsm('L', 'U', 'T', 'N', C%nrow, C%ncol, ONE, &
                     R%a, R%narow, C%a, C%narow)
         call strsm('L', 'U', 'N', 'N', C%nrow, C%ncol, ONE, &
                     R%a, R%narow, C%a, C%narow)
#endif

         ! Set the type of C.

         C%tag = 'GE'
         ! If necessary, deallocate.

         if (.not.present(chl)) then
            call Clean(chla)
         end if

      case default

         call Matranerror('RmSolveXiy in RmatSolve: Illegal tag.')

      end select

      ! Deallocate storage of temporaries.

      call CleanTemp(A)
      call CleanTemp(B)

   end subroutine RmSolveXiy

!  Computes C = A.xiy.B.

   function RmSolve_xiy(A, B) result(C)
      type(Rmat) :: C
      type(Rmat), intent(in) :: A, B

      C%a => null()
      C%temporary => null()
      call Clean(C)
      call SolveXiy(C, A, B)
      call SetTemp(C)
   end function RmSolve_xiy

! Xihy

   subroutine RmSolveXihy(C, A, B, lu)
      type(Rmat), intent(inout) :: C
      type(Rmat), intent(in) :: A
      type(Rmat), intent(in) :: B
      type(RmatLudpp), intent(inout), optional, target :: lu

      integer :: i
      type(RmatLudpp), target :: lua
      type(RmatLudpp), pointer :: lud


      if (A%nrow==0 .or. A%ncol==0)&
         call MatranError('RmSolveXihy in RmatSolve:&
                            & Cannot invert a null matrix.')

      if (A%nrow/=A%ncol .or. A%nrow/=B%nrow)&
         call MatranError('RmSolveXihy in RmatSolve:&
                            & Inconsistent dimensions.')

      call GuardTemp(A)
      call GuardTemp(B)

      C = B

      select case(A%tag)

      case('LT')

         do i=1, A%nrow
            if (A%a(i,i) == 0)&
               call MatranError('RmSolveXihy in RmatSolve:&
                            & Lower triangular A is singular.')
         end do

#ifdef dbl
         call dtrsm('L', 'L', 'T', 'N', C%nrow, C%ncol, ONE,&
                     A%a, A%narow, C%a, C%narow)
#endif
#ifdef sngl
         call strsm('L', 'L', 'T', 'N', C%nrow, C%ncol, ONE,&
                     A%a, A%narow, C%a, C%narow)
#endif


         if (B%tag == 'UT') then
            C%tag = 'UT'
         else
            C%tag = 'GE'
         end if


      case('UT')

         do i=1, A%nrow
            if (A%a(i,i) == 0)&
               call MatranError('RmSolveXihy in RmatSolve:&
                            & Upper triangular A is singular.')
         end do

#ifdef dbl
         call dtrsm('L', 'U', 'T', 'N', C%nrow, C%ncol, ONE,&
                     A%a, A%narow, C%a, C%narow)
#endif
#ifdef sngl
         call strsm('L', 'U', 'T', 'N', C%nrow, C%ncol, ONE,&
                     A%a, A%narow, C%a, C%narow)
#endif

         if (B%tag == 'LT') then
            C%tag = 'LT'
         else
            C%tag = 'GE'
         end if


      case('GE')

         if (present(lu)) then
            if (.not.lu%companion) call Ludpp(lu, A)
            lud => lu
         else
            call Ludpp(lua, A)
            lud => lua
         end if

         do i=1, lud%U%nrow
            if (lud%U%a(i,i) == 0)&
               call MatranError('RmSolveXihy in RmatSolve:&
                            & U-factor of A is singular.')
         end do

#ifdef dbl
         call dtrsm('L', 'U', 'T', 'N', C%nrow, C%ncol, ONE, &
                     lud%U%a, lud%U%narow, C%a, C%narow)
         call dtrsm('L', 'L', 'T', 'N', C%nrow, C%ncol, ONE, &
                     lud%L%a, lud%L%narow, C%a, C%narow)
#endif
#ifdef sngl
         call strsm('L', 'U', 'T', 'N', C%nrow, C%ncol, ONE, &
                     lud%U%a, lud%U%narow, C%a, C%narow)
         call strsm('L', 'L', 'T', 'N', C%nrow, C%ncol, ONE, &
                     lud%L%a, lud%L%narow, C%a, C%narow)
#endif

         call PivotInvRow(C, lud%pvt, lud%npvt)

         C%tag = 'GE'

         if (.not.present(lu)) then
            call Clean(lud)
         end if
  

      case default
         call MatranError('RmSolveXihy in RmatSolve:&
                            & Operator .xihy. not implemented for&
                            & Rmats of type HE or HP. Use .xiy..')

      end select

      ! Deallocate temporaries.

      call CleanTemp(A)
      call CleanTemp(B)


   end subroutine RmSolveXihy
   
   ! xihy
   
   function RmSolve_xihy(A, B) result(C)
      type(Rmat) :: C
      type(Rmat), intent(in) :: A, B

      C%a => null()
      C%temporary => null()
      call Clean(C)
      call SolveXihy(C, A, B)
      call SetTemp(C)
   end function RmSolve_xihy

   ! Xyi

   subroutine RmSolveXyi(C, A, B, lu, chl)
      type(Rmat), intent(inout) :: C
      type(Rmat), intent(in) :: A
      type(Rmat), intent(in) :: B
      type(RmatLudpp), intent(inout), optional, target :: lu
      type(RmatChol), intent(inout), optional, target :: chl

      integer :: i
      type(RmatLudpp), target :: lub
      type(RmatLudpp), pointer :: lud
      type(RmatChol), target :: chlb
      type(Rmat), pointer :: R

      if (B%nrow==0 .or. B%ncol==0)&
         call MatranError('RmSolveXyi in RmatSolve:&
                            & Cannot invert a null matrix.')

      if (B%nrow/=B%ncol .or. A%ncol/=B%nrow)&
         call MatranError('RmSolveXyi in RmatSolve:&
                            & Inconsistent dimensions.')

      call GuardTemp(A)
      call GuardTemp(B)

      C = A

      select case(B%tag)

      case('LT')

         do i=1, B%nrow
            if (B%a(i,i) == 0)&
               call MatranError('RmSolveXyi in RmatSolve:&
                            & Lower triangularB is singular.')
         end do

#ifdef dbl
         call dtrsm('R', 'L', 'N', 'N', C%nrow, C%ncol, ONE,&
                     B%a, B%narow, C%a, C%narow)
#endif
#ifdef sngl
         call strsm('R', 'L', 'N', 'N', C%nrow, C%ncol, ONE,&
                     B%a, B%narow, C%a, C%narow)
#endif

         if (A%tag == 'LT') then
            C%tag = 'LT'
         else
            C%tag = 'GE'
         end if


      case('UT')

         do i=1, B%nrow
            if (B%a(i,i) == 0)&
               call MatranError('RmSolveXyi in RmatSolve:&
                            & Lower triangular B is singular.')
         end do

#ifdef dbl
         call dtrsm('R', 'U', 'N', 'N', C%nrow, C%ncol, ONE,&
                     B%a, B%narow, C%a, C%narow)
#endif
#ifdef sngl
         call strsm('R', 'U', 'N', 'N', C%nrow, C%ncol, ONE,&
                     B%a, B%narow, C%a, C%narow)
#endif

         if (A%tag == 'UT') then
            C%tag = 'UT'
         else
            C%tag = 'GE'
         end if

      case('GE', 'HE')

         if (present(lu)) then
            if (.not.lu%companion) call Ludpp(lu, B)
            lud => lu
         else
            call Ludpp(lub, B)
            lud => lub
         end if

         do i=1, lud%U%nrow
            if (lud%U%a(i,i) == 0)&
               call MatranError('RmSolveXiy in RmatSolve:&
                            & U-factor of B is singular.')
         end do

#ifdef dbl
         call dtrsm('R', 'U', 'N', 'N', C%nrow, C%ncol, ONE, &
                     lud%U%a, lud%U%narow, C%a, C%narow)
         call dtrsm('R', 'L', 'N', 'N', C%nrow, C%ncol, ONE, &
                     lud%L%a, lud%L%narow, C%a, C%narow)
#endif
#ifdef sngl
         call strsm('R', 'U', 'N', 'N', C%nrow, C%ncol, ONE, &
                     lud%U%a, lud%U%narow, C%a, C%narow)
         call strsm('R', 'L', 'N', 'N', C%nrow, C%ncol, ONE, &
                     lud%L%a, lud%L%narow, C%a, C%narow)
#endif

         call PivotInvCol(C, lud%pvt, lud%npvt)

         C%tag = 'GE'

         if (.not.present(lu)) then
            call Clean(lud)
         end if
   

      case('HP')

         if (present(chl)) then
            if (.not.chl%companion) call Chol(chl, B)
            R => chl%R
         else
            call Chol(chlb, B)
            R => chlb%R
         end if

         do i=1, R%nrow
            if (R%a(i,i) == 0)&
               call MatranError('RmSolveXiy in RmatSolve:&
                            & Cholesky-factor of B is singular.')
         end do

#ifdef dbl
         call dtrsm('R', 'U', 'N', 'N', C%nrow, C%ncol, ONE, &
                     R%a, R%narow, C%a, C%narow)
         call dtrsm('R', 'U', 'T', 'N', C%nrow, C%ncol, ONE, &
                     R%a, R%narow, C%a, C%narow)
#endif
#ifdef sngl
         call strsm('R', 'U', 'N', 'N', C%nrow, C%ncol, ONE, &
                     R%a, R%narow, C%a, C%narow)
         call strsm('R', 'U', 'T', 'N', C%nrow, C%ncol, ONE, &
                     R%a, R%narow, C%a, C%narow)
#endif

         C%tag = 'GE'

         if (.not.present(chl)) then
            call Clean(chlb)
         end if
  
      case default

         call Matranerror('RmSolveXyi in RmatSolve: Illegal tag.')

      end select

      ! Deallocate temp vars

      call CleanTemp(A)
      call CleanTemp(B)

   end subroutine RmSolveXyi

   ! xyi

   function RmSolve_xyi(A, B) result(C)
      type(Rmat) :: C
      type(Rmat), intent(in) :: A, B

      C%a => null()
      C%temporary => null()
      call Clean(C)
      call SolveXyi(C, A, B)
      call SetTemp(C)

   end function RmSolve_xyi

   ! Xyih

   subroutine RmSolveXyih(C, A, B, lu)
      type(Rmat), intent(inout) :: C
      type(Rmat), intent(in) :: A
      type(Rmat), intent(in) :: B
      type(RmatLudpp), intent(inout), optional, target :: lu

      integer :: i
      type(RmatLudpp), target :: lub
      type(RmatLudpp), pointer :: lud

      if (B%nrow==0 .or. B%ncol==0)&
         call MatranError('RmSolveXyih in RmatSolve:&
                            & Cannot invert a null matrix.')

      if (B%nrow/=B%ncol .or. A%ncol/=B%nrow)&
         call MatranError('RmSolveXyi in RmatSolve:&
                            & Inconsistent dimensions.')

      call GuardTemp(A)
      call GuardTemp(B)
      C = A

      select case(B%tag)

      case('LT')

         do i=1, B%nrow
            if (B%a(i,i) == 0)&
               call MatranError('RmSolveXyi in RmatSolve:&
                            & Lower triangular B is singular.')
         end do

#ifdef dbl
         call dtrsm('R', 'L', 'T', 'N', C%nrow, C%ncol, ONE,&
                     B%a, B%narow, C%a, C%narow)
#endif
#ifdef sngl
         call strsm('R', 'L', 'T', 'N', C%nrow, C%ncol, ONE,&
                     B%a, B%narow, C%a, C%narow)
#endif

         if (A%tag == 'LT') then
            C%tag = 'LT'
         else
            C%tag = 'GE'
         end if

      case('UT')

         do i=1, B%nrow
            if (B%a(i,i) == 0)&
               call MatranError('RmSolveXyih in RmatSolve:&
                            & Upper triangular B is singular.')
         end do

         
#ifdef dbl
         call dtrsm('R', 'U', 'T', 'N', C%nrow, C%ncol, ONE,&
                     B%a, B%narow, C%a, C%narow)
#endif
#ifdef sngl
         call strsm('R', 'U', 'T', 'N', C%nrow, C%ncol, ONE,&
                     B%a, B%narow, C%a, C%narow)
#endif

         if (A%tag == 'UT') then
            C%tag = 'UT'
         else
            C%tag = 'GE'
         end if

      case('GE')

         if (present(lu)) then
            if (.not.lu%companion) call Ludpp(lu, B)
            lud => lu
         else
            call Ludpp(lub, B)
            lud => lub
         end if

         do i=1, lud%U%nrow
            if (lud%U%a(i,i) == 0)&
               call MatranError('RmSolveXiyh in RmatSolve:&
                            & U-factor of B is singular.')
         end do

         call PivotCol(C, lud%pvt, lud%npvt)
#ifdef dbl
         call dtrsm('R', 'L', 'T', 'N', C%nrow, C%ncol, ONE, &
                     lud%L%a, lud%L%narow, C%a, C%narow)
         call dtrsm('R', 'U', 'T', 'N', C%nrow, C%ncol, ONE, &
                     lud%U%a, lud%U%narow, C%a, C%narow)
#endif
#ifdef sngl
         call strsm('R', 'L', 'T', 'N', C%nrow, C%ncol, ONE, &
                     lud%L%a, lud%L%narow, C%a, C%narow)
         call strsm('R', 'U', 'T', 'N', C%nrow, C%ncol, ONE, &
                     lud%U%a, lud%U%narow, C%a, C%narow)
#endif

         C%tag = 'GE'

         if (.not.present(lu)) then
            call Clean(lud)
         end if
  
      case default
         call MatranError('RmSolveXyih in RmatSolve:&
                            & Operator .xyih. not implemented for&
                            & Rmats of type HE or PO.  Use .xyi.')         

      end select

      ! Deallocate temp vars

      call CleanTemp(A)
      call CleanTemp(B)

   end subroutine RmSolveXyih

   function RmSolve_xyih(A, B) result(C)
      type(Rmat) :: C
      type(Rmat), intent(in) :: A, B

      C%a => null()
      C%temporary => null()
      call Clean(C)
      call SolveXyih(C, A, B)
      call SetTemp(C)

   end function RmSolve_xyih



end module RmatSolve_m
