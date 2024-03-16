module MatranUtil_m
implicit none

#ifdef Overview

MatranUtil_m is the root the MatWrap module.  It contains global
constants and subprograms.  Currently it contains a parameter to
define the kind for working precision, error handlers, and routines to
reallocate storage for arrays.

Authors: CheRung Lee and Pete Stewart
Aug 21 2003

#endif

   ! Kind for double precision.

#ifdef sngl
   integer, parameter :: wp = kind(1.0e0)
#endif

#ifdef dbl
   integer, parameter :: wp = kind(1.0d0)
#endif

   ! ReshapeAry is a generic function widely used in Matran.
   ! Specifically, given one or two dimensions, depending on the kind
   ! of array, ReshapAry determines if the the array can hold a
   ! subarray of those dimensions.  If so it returns.  Otherwise it
   ! reallocates the array to be of the size of the input dimensions.
   ! In both cases the array is set to zero.

   interface ReshapeAry
      module procedure ReshapeAryI1, ReshapeAryI2
      module procedure ReshapeAryR1, ReshapeAryR2
      module procedure ReshapeAryC1, ReshapeAryC2
   end interface


contains

   ! Since Matran uses BLAS and LAPACK routines, it has
   ! two subroutines to print error messages.  The first
   ! prints error messages that result from Matran's own
   ! error checking.  The second prints appropriate an
   ! error message, when an LAPACK routine gives an error
   ! return.  Since LAPACK communicates its status on completion
   ! by the parameter info, the integer info is also printed
   ! out.  Placing an optional info parameter in the
   ! calling sequence of the Matran program that calls
   ! will cause Matran to suppress the error message and
   ! return the value of info from the LAPACK routine.

   ! Matran error handler

   subroutine MatranError(ErrorMessage)
      character(*), intent(in) :: ErrorMessage

      print *, ErrorMessage
      stop

   end subroutine MatranError

   ! Support error handler.

   subroutine SupportError(ErrorMessage, info)
      character(*), intent(in) :: ErrorMessage
      integer, intent(in) :: info

      print *, ErrorMessage, info
      stop

   end subroutine SupportError


   ! The following six subroutines implement ReshapeAry for
   ! linear and rectangular integer, real, and complex arrays.

   ! Reshape a linear integer array

   subroutine ReshapeAryI1(Ary, n)
      integer, pointer    :: Ary(:)
      integer, intent(in) :: n

      if (associated(Ary)) then
         if (size(Ary)<n) then
            deallocate(Ary)
            allocate(Ary(n))
         end if
      else
         allocate(Ary(n))
      end if
      Ary = 0      
   end subroutine ReshapeAryI1

   ! Reshape a rectangular integer array

   subroutine ReshapeAryI2(Ary, m, n)
      integer, pointer    :: Ary(:,:)
      integer, intent(in) :: m, n

      integer :: shp(2)
  
      if (associated(Ary)) then
         shp = shape(Ary)
         if (m>shp(1) .or. n>shp(2)) then
            deallocate(Ary)
            allocate(Ary(m, n))
         end if
      else
         allocate(Ary(m, n))
      end if
      Ary = 0
   end subroutine ReshapeAryI2

   ! Reshape a linear real arrray

   subroutine ReshapeAryR1(Ary, n)
      real(wp), pointer   :: Ary(:)
      integer, intent(in) :: n

      if (associated(Ary)) then
         if (size(Ary)<n) then
            deallocate(Ary)
            allocate(Ary(n))
         end if
      else
         allocate(Ary(n))
      end if
      Ary = 0.0
   end subroutine ReshapeAryR1

   ! Reshape a rectangular double-precision arrray

   subroutine ReshapeAryR2(Ary, m, n)
      real(wp), pointer   :: Ary(:,:)
      integer, intent(in) :: m, n

      integer :: shp(2)
  
      if (associated(Ary)) then
         shp = shape(Ary)
         if (m>shp(1) .or. n>shp(2)) then
            deallocate(Ary)
            allocate(Ary(m, n))
         end if
      else
         allocate(Ary(m, n))
      end if
      Ary = 0.0
   end subroutine ReshapeAryR2

   ! Reshape a linear complex double-precision arrray

   subroutine ReshapeAryC1(Ary, n)
      complex(wp), pointer :: Ary(:)
      integer, intent(in)   :: n

      if (associated(Ary)) then
         if (size(Ary)<n) then
            deallocate(Ary)
            allocate(Ary(n))
         end if
      else
         allocate(Ary(n))
      end if
      Ary = 0.0
   end subroutine ReshapeAryC1

   ! Reshape a rectangular complex double-precision arrray

   subroutine ReshapeAryC2(Ary, m, n)
      complex(wp), pointer  :: Ary(:, :)
      integer               :: m, n
      integer               :: shp(2)

      if (associated(Ary)) then
         shp = shape(Ary)
         if (m>shp(1) .or. n>shp(2)) then
            deallocate(Ary)
            allocate(Ary(m, n))
         end if
      else
         allocate(Ary(m, n))
      end if
      Ary = 0.0
   end subroutine ReshapeAryC2

end module MatranUtil_m

