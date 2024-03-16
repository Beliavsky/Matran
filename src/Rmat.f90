module Rmat_m
use MatranUtil_m
implicit none

#ifdef OVERVIEW

This module implements a struture Dmat representing a real
nrow x ncol matrix contained in an narow x nacol array.  Indexing for
the array starts at (1,1).  If A is a Dmat, we distinguish between the
object A, the matrix A, and the array A.


Components

   real(wp), pointer :: ary(:,:) = null()
      The array containing the matrix.

   integer :: nrow = 0
      The number of rows in the matrix.

   integer :: ncol = 0
      The number of columns in the matrix.

   integer :: narow = 0
      The number of rows in the array.

   integer :: nacol = 0
      The number of columns in the array.

   character(2) :: tag = 'GE'
      The type of the matrix.  There are five types.

         GE -- General matrix
         LT -- Lower triangular matrix
         UT -- Upper triangular matrix
         HE -- Hermitian matrix
         HP -- Hermitian positive (semi) definite matrix

      The package contains no packed forms.  For an example, a lower
      triangular matrix is stored with explicit zeros in its lower
      half.  It is the responsibility of the user to establish and
      maintain the pattern of such a matrix.

   logical :: adjustable = .true.
      If adjustable is .true. the interface is permitted to adjust the
      shape of the array to accommodate computed results.  For example,
      if we write C = A + B and one of the dimensions of the *matrix*
      A and B are greater than the corresponding dimensions of the
      *array* C, then the function implementing + (in RmatSum) will
      reallocate C%array to contain the sum, but only if adjustable is
      .true.
   
   integer, pointer: temporary => 0
      If temporary is larger than one, the Ddiag is an intermediate
      result in the value of an expression or the result of a function.
      Its array must be deallocated by the last routine that uses it.

Support

   subroutine RmSetTemp(A) (instance of SetTemp)

      Makes a Rmat a temporary


   subroutine RmGuardTemp(A) (instance of GuardTemp)

      Protects a temporary from premature deallocation

   subroutine RmCleanTemp(A) (instance of CleanTemp)

      Deallocates the array of a temporary Rmat, if necessary.

   subroutine RmReshapAry(A, m, n) (instance of ReshapeAry)

      Finds enough storage to hold an mxn matrix.


   subroutine RmEqualsRm(A, B)

      Overloads = to assign B to A.

   subroutine RmEqualsAry(A, ary)

      Overloads = to copy ary to A%a

   function RmFromAry(ary) result(C)

      Overloads .rm. to produce C = ary

   subroutine RmEqualsRowCol(A, dims)

      Overloads = to initialize A to a zero dims(1) x dims(2) matrx
      in a dims(3) x dims(4) array.

   function RmFromRowCol(dims) result(C)

      Overloads .rm. to produce C = dims

   subroutine RmEqualsRs(A, s)

      Overloads = to create a 1x1 matrix with element s.

   function RmFromRs(s) result(C)

      Overloads .rm. to produce C = s.

Author: Pete Stewart
Aug 20 2003

#endif

   type Rmat
      real(wp), pointer &               ! The matrix array
              :: a(:,:) => null()       !
      integer :: nrow = 0               ! Number of rows in the matrix
      integer :: ncol = 0               ! Number of columns in the matrix
      integer :: narow = 0              ! Number of rows in the array
      integer :: nacol = 0              ! Number of columns in the array
      character(2) &                    ! Type of matrix
              :: tag = 'GE'             !
      logical :: adjustable =.true.     ! Adjustable array
      integer,pointer &                 ! Intermediate value
              :: temporary => null()      
   end type Rmat

   interface assignment (=)
      module procedure RmEqualsRm, RmEqualsAry, RmEqualsRowCol, RmEqualsRs
   end interface

   interface operator (.rm.)
      module procedure RmFromRs, RmFromAry, RmFromRowCol
   end interface

   interface ReshapeAry
      module procedure RmReshapeAry
   end interface

   interface Clean
      module procedure RmClean
   end interface

   interface CleanTemp
      module procedure RmCleanTemp
   end interface

   interface GuardTemp
      module procedure RmGuardTemp
   end interface

   interface SetTemp
      module procedure RmSetTemp
   end interface



contains

   ! Sanitizer for Rmat.

   subroutine RmClean(A)

      type(Rmat), intent(inout) :: A

      if (associated(A%a)) deallocate(A%a)
      A%nrow = 0
      A%ncol = 0
      A%narow = 0
      A%nacol = 0
      A%tag = 'GE'
      A%adjustable = .true.

      if(associated(A%a)) &
           deallocate(A%a)
      if (associated(A%temporary)) &
           deallocate(A%temporary)

   end subroutine RmClean

   ! If A is temporary and A%temporary>1, decrement A%temporary.  If
   ! A%temporary == 1 after decrememtation, deallocate A%a and
   ! A%temporary

   subroutine RmCleanTemp(A)
      type(Rmat) :: A
      integer, pointer :: t
      real(wp), pointer :: s(:,:)

      if(associated(A%temporary)) then
         t=>A%temporary
         if(t > 1)  t = t - 1

         if (t == 1) then
            s=>A%a
            deallocate(s)
            deallocate(t)
         end if
      end if
   end subroutine RmCleanTemp
   


   ! Increase A%temporary if A is temporary.
   
   subroutine RmGuardTemp(A)
      type(Rmat) :: A
      integer, pointer :: t
      
      if(associated(A%temporary)) then
         t=>A%temporary
         t = t + 1
      end if
   end subroutine RmGuardTemp


   ! Set A to be a temporary variable.

   subroutine RmSetTemp(A)
      type(Rmat), intent(inout) :: A

      if(.not. associated(A%temporary)) &
           allocate( A%temporary)
      A%temporary = 1

   end subroutine RmSetTemp


   ! RmEqualsRm overloads = to initialize Rmat A to Rmat B.

   subroutine RmEqualsRm(A, B)
      type(Rmat), intent(inout) :: A
      type(Rmat), intent(in) :: B

      call GuardTemp(B)

      call ReshapeAry(A, B%nrow, B%ncol)

      A%a(1:A%nrow, 1:A%ncol) = B%a(1:B%nrow,1:B%ncol)
      A%tag = B%tag
      call CleanTemp(B)

   end subroutine RmEqualsRm


   ! RmEqualsAry overloads = to set Rmat A%a to ary.

   subroutine RmEqualsAry(A, ary)

      type(Rmat), intent(inout) :: A
      real(wp), intent(in) :: ary(:,:)

      integer :: shp(2)
  
      shp = shape(ary)

      call ReshapeAry(A, shp(1), shp(2))

      A%a(1:A%nrow, 1:A%ncol) = ary


   end subroutine RmEqualsAry

   ! RmFromAry overloads .rm. to produce C = ary.

   function RmFromAry(ary) result(C)
      type(Rmat) :: C
      real(wp), intent(in) :: ary(:,:)


      C%a => null()
      C%temporary => null()
      call Clean(C)

      C = ary
      call SetTemp(C)

   end function RmFromAry
   
   ! RmEqualsRowCol overloads = to create
   ! an m x n zero matrix in a ma x na array, as
   ! specified by the array dims.  If size(dims) = 2,
   ! m = ma = dims(1) and n = na = dims(2).  If size(dims) = 4,
   ! m = dims(1), n = dims(2), ma = dims(3), na = dims(4).

   subroutine RmEqualsRowCol(A, dims)

      type(Rmat), intent(inout) :: A
      integer, intent(in) :: dims(:)

      integer sz

      sz = size(dims)

      if (sz/=2 .and. sz/=4) then
         call MatranError("RmEqualsRowCol in Rmat: Wrong number of &
                          &inputs.")
      end if

      if (sz==2) then
         if (dims(1)<0 .or. dims(2)<0) &
            call MatranError("RmEqualsRowCol in Rmat: Inconsistent &
                               &dimensions.")
         call ReshapeAry(A, dims(1), dims(2))
      else
         if (dims(1)<0 .or. dims(2)<0 .or. &
            dims(1)>dims(3) .or. dims(2)>dims(4)) &
            call MatranError("RmRowCol in Rmat: Inconsistent &
                             &dimensions.")
         A%narow = dims(3)
         A%nacol = dims(4)
         if (associated(A%a)) deallocate(A%a)
         allocate(A%a(A%narow, A%nacol))
         A%nrow = dims(1)
         A%ncol = dims(2)
         A%a = 0.0
      end if
      A%tag = 'GE'

   end subroutine RmEqualsRowCol

   function RmFromRowCol(dims) result(C)
      type(Rmat) :: C
      integer, intent(in) :: dims(:)

      C%a => null()
      C%temporary => null()
      call Clean(C)

      C = dims
      call SetTemp(C)

   end function RmFromRowCol

   ! RmEqualsRs overloads = to create a 1x1 Rmat whose element
   ! is the scalar s.

   subroutine RmEqualsRs(A, s)

      type(Rmat), intent(inout) :: A
      real(wp), intent(in) :: s
      
      Call ReshapeAry(A, 1, 1)
      A%a(1,1) = s
      A%tag = 'GE'

   end subroutine RmEqualsRs

   ! RmFromRs overloads .rm. to produce a 1x1 Rmat whose element
   ! is the scalar s.

   function RmFromRs(s) result(C)
      type(Rmat) :: C
      real(wp), intent(in) :: s

      C%a => null()
      C%temporary => null()
      call Clean(C)

      C = s

      call SetTemp(C)

   end function RmFromRs

   ! RmReshapeAry reshapes the array of a Rmat.

   subroutine RmReshapeAry(A, m, n)
      type(Rmat) :: A
      integer :: m, n

      ! RmReshapAry reinitilizes A to be and  mxn zero matrix
      ! in an mxn array.  Type and adjust are left unchanged.

      if (.not.associated(A%a)) then

         ! A is clean, allocate storage for ary.

         allocate(A%a(m, n))
         A%narow = m
         A%nacol = n
      else
         if (m>A%narow .or. n>A%nacol) then

            ! Current Rmat A does not have enough storage.
            ! Allocate if adjustable.  Otherwise throw error message.

            if (A%adjustable) then
               deallocate(A%a)
               allocate(A%a(m, n))
               A%narow = m
               A%nacol = n
            else
               call MatranError("RmReshapeAry in Rmat: Cannot change ary &
                            &in an unadjustable Rmat.")
            end if
         end if
      end if
      A%nrow = m
      A%ncol = n
      A%a = 0.0

   end subroutine RmReshapeAry

end module Rmat_m
