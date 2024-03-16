module Rdiag_m
use MatranUtil_m
use Rmat_m
implicit none

#ifdef OVERVIEW

This module implements a type Rmat that represents a diagonal
matrix of order order contained in a one-dimensional array of
length na.

Components

   real(wp), pointer & :: a(:) => null()
      The array containing the diagonal matrix

   integer :: order
      The order of the matrix

   integer :: na
      The length of ary

   logical :: adjustable = .true.  
      If adjustable is .true. a program is permitted to adjust the
      shape of the array to accomodate computed results.  For example,
      if we write F = D + E the corresponding order of the diagonal
      matrices A and B are greater than the order of the *array* C,
      then the function implementing + (in RdiagSum) will reallocate
      C%array to contain the sum, but only if adjustable is .true.

   integer :: temporary = 0
      If temporary is larger than one, the Rdiag is an intermediate 
      result in the value of an expression or the rusult of a function.
      Its array must be deallocated by the last routine that uses it.

Support

   subroutine RdClean(D)  (Instance of Clean)

      Deallocates D%a, if necessary, and reinitializes the
      compontents of A.

   subroutine RdSetTemp(D) (Instance of SetTemp)

      Makes a Rdiag a temporary.

   subroutine RdGuardTemp(D)  (Instance of GuardTemp)

      Protects a temporary from premature deallocation

   subroutine RdCleanTemp(D)  (Instance of CleanTemp)

      Deallocates the array of a temporary Rdiag, if necessary.

   subroutine RdReshapeAry(D, n) (Instance of ReshapeAry)

      If necessary, allocates enough storage for D%a to hold
       a Rdiag of order n and initilizes D%a to zero.

   subroutine RdEqualsRd(D, E)

      Overloads = to assign the Rdiag E to Rdiag D.

   subroutine RdEqualsAry(D, a)

      Overloads = to create in Rdiag B a Rdiag whose ary contents
         are those of the array a.

   subroutine RmEqualsRd(A, D)

      Overloads = to create a diagonal Rmat A whose diagonal is D.

   function RmFromRd(D) result(A)

      Overloads .rm. to produce A = D

   function RdFromAry(ary) result(D)

      RdFromAry overloads .rd. to produce D = ary.

   subroutine RdEqualsOrder(D, dim)

      Overloads = to create a zero Didiag of order dim(1) in
                  an array of length dim(2).

   function RdFromOrder(order) result(D)

      RdFromOrder overloads .rd. to produce D = dim.

   subroutine RdEqualsDs(D, s)

      Overloads = to create a Rdiag of order 1 whose element
      is the scalar s.

   function RdFromDs(s) result(D)

      Overloads .rd. to produce a Rdiag of order 1 whose element
      is the scalar s.


Author: Pete Stewart
Apr 30 2003

#endif      

   type Rdiag
      real(wp), pointer &             ! The matrix array
           :: a(:) => null()
      integer :: order = 0            ! The order of the matrix
      integer :: na = 0               ! The length of the array
      logical :: adjustable = .true.  ! Adjustable array
      integer, pointer&               ! Intermediate value
           :: temporary => null()     
   end type Rdiag

   interface Clean
      module procedure RdClean
   end interface

   interface ReshapeAry
      module procedure RdReshapeAry
   end interface

   interface CleanTemp
      module procedure RdCleanTemp
   end interface

   interface GuardTemp
      module procedure RdGuardTemp
   end interface

   interface SetTemp
      module procedure RdSetTemp
   end interface


   interface assignment (=)
      module procedure RdEqualsRd, RdEqualsAry, RmEqualsRd,&
                       RdEqualsOrder, RdEqualsDs
   end interface

   interface operator (.rd.)
      module procedure RdFromAry, RdFromOrder, RdFromDs
   end interface 

   interface operator (.rm.)
      module procedure RmFromRd
   end interface

contains
   
   ! Sanitizer for Rdiag, implementing Clean

   subroutine RdClean(D)

      ! DcClean deallocates the storage for D%a and reinitializaes
      ! the Rdiag to default tolerences.

      type(Rdiag), intent(inout) :: D

      if (associated(D%a)) then
         deallocate(D%a)
      end if
      D%a => null()
      D%order = 0
      D%na = 0
      D%adjustable = .true.
      if (associated(D%temporary)) deallocate(D%temporary)
      D%temporary => null()
   end subroutine RdClean

   ! Reshape array for Rdiag.  Implements ReshapeAry

   subroutine RdReshapeAry(D, order)

      ! RdReshapeAry enlarges D%a, if necessary, to hold a
      ! Rdiag of size order.  Will throw error return if
      ! input size of D%a is insufficient and D is not
      ! adjustable.

      type(Rdiag) :: D
      integer     :: order

      if (order < 0)&
            call MatranError("RdReshapeAry in Rdiag: Negative order.")

      if (.not.associated(D%a)) then
         if (D%adjustable) then
            allocate(D%a(order))
            D%na = order
         else
            call MatranError("RdReshapeAry in Rdiag: cannot change&
                              & array in an unadjustable Rdiag.")
         end if
      else
         if (order > D%na) then
            if (D%adjustable) then
               deallocate(D%a)
               allocate(D%a(order))
               D%na = order
            else
               call MatranError("RdReshapeAry in Rdiag: cannot change&
                              & arrayy in an unadjustable Rdiag.")
            end if
         end if
      end if

      D%order = order
      D%a = 0.0

   end subroutine RdReshapeAry

   ! Deallocates the array of a temporary Rdiag.

   subroutine RdCleanTemp(D)
      type(Rdiag) :: D
      integer, pointer :: t=>null()
      real(wp), pointer :: a(:)=>null()

      if(associated(D%temporary)) then
         t=>D%temporary

         if (t > 1)&
            t = t - 1
         
         if (t == 1) then
            a => D%a
            deallocate(a)
            deallocate(t)
         end if
      end if
   end subroutine RdCleanTemp


   ! Increase the D%temporary if D is temporary

   subroutine RdGuardTemp(D)
      type(Rdiag) :: D
      integer, pointer :: t=>null()
      
      if(associated(D%temporary)) then
         t=>D%temporary
         t = t + 1
      end if
   end subroutine RdGuardTemp


   ! Set D to be a temp variable.

   subroutine RdSetTemp(D)
      type(Rdiag) :: D

      integer, pointer :: t=>null()

      if(.not. associated(D%temporary)) &
           allocate( D%temporary)
      t => D%temporary
      t = 1

   end subroutine RdSetTemp

   ! Implements Rdiag = Rdiag.

   subroutine RdEqualsRd(D, E)

      ! RdEqualsRd sets D to a Rdiag whose diagonal is the
      ! same as the diagonal of E.

      type(Rdiag), intent(inout) :: D
      type(Rdiag), intent(in) :: E

      call GuardTemp(E)
      call ReshapeAry(D, E%order)

      D%a = 0
      D%a(1:D%order) = E%a(1:E%order)

      call CleanTemp(E)
   end subroutine RdEqualsRd
      
   !  Implements Rdiag = array.

   subroutine RdEqualsAry(D, a)

      ! Rdiag sets D to a Rdiag whose ary is equal to a.

      type(Rdiag), intent(inout) :: D
      real(wp), intent(in) :: a(:)

      integer :: order


      order = size(a)
      call ReshapeAry(D, order)

      D%a = 0
      D%a(1:order) = a
      

   end subroutine RdEqualsAry

   ! RdFromAry overloads .rd. to produce D = ary

   function RdFromAry(ary) result(D)
      type(Rdiag) :: D
      real(wp), intent(in) :: ary(:)

      D%a => null()
      D%temporary => null()
      call Clean(D)

      D = ary
      call SetTemp(D)

   end function RdFromAry


   ! Implements Rmat = Rdiag.

   subroutine RmEqualsRd(A, D)

      ! RmEqualsRd sets a to a diagonal Rmat whose diagonal
      ! is equal to the diagoal of D.

      type(Rmat), intent(inout) :: A
      type(Rdiag), intent(in) :: D

      integer i, n

      n = D%order

      call GuardTemp(D)
      call ReshapeAry(A, n, n)

      do i=1,n
         A%a(i,i) = D%a(i)
      end do

      call CleanTemp(D)
   end subroutine RmEqualsRd

   ! RmFromRd overloads .rm. to produce A = D

   function RmFromRd(D) result(A)
      type(Rmat) :: A
      type(Rdiag), intent(in) :: D

      A%a => null()
      A%temporary => null()
      call Clean(A)

      A = D

      call SetTemp(A)

   end function RmFromRd

   ! Implements Rmat = order.

   subroutine RdEqualsOrder(D, dims)

      ! RdEqualsOrder creates a zero Rdiag whose size is
      ! equal to order.

      type(Rdiag), intent(inout) :: D
      integer, intent(in) :: dims(:)

      integer :: sz

      sz = size(dims)

      if (sz/=1 .and. sz/=2)&
         call MatranError("RdEqualsOrder in Rdiag: Wrong number of &
                           &inputs.")

      if (sz == 1) then
         if (dims(1) < 0)&
              call MatranError("RdEqualsOrder in Rdiag: Negative order.")
         call ReshapeAry(D, dims(1))
      else
         if (dims(1)<0 .or. dims(1) > dims(2))&
              call MatranError("RdEqualsOrder in Rdiag: Inconsistent &
                   &Dimensions.")
         D%na = dims(2)
         if (associated(D%a)) deallocate(D%a)
         allocate(D%a(D%na))
         D%order = dims(1)
         D%a = 0
         call ReshapeAry(D, dims(1))
      end if

   end subroutine RdEqualsOrder

   ! RdFromOrder overloads .rd. to produce D = order

   function RdFromOrder(order) result(D)
      type(Rdiag) :: D
      integer, intent(in) :: order(:)

      D%a => null()
      D%temporary => null()
      call Clean(D)

      D = order
      call SetTemp(D)

   end function RdFromOrder

   ! RdEqualsDs overloads = to create a Rdiag of order 1 whose element
   ! is the scalar s.

   subroutine RdEqualsDs(D, s)

      type(Rdiag), intent(inout) :: D
      real(wp), intent(in) :: s

      Call ReshapeAry(D, 1)
      D%a(1) = s

   end subroutine RdEqualsDs

   ! RdFromDs overloads .rd. to produce a Rdiag of order 1 whose element
   ! is the scalar s.

   function RdFromDs(s) result(D)
      type(Rdiag) :: D
      real(wp), intent(in) :: s

      D%a => null()
      D%temporary => null()
      call Clean(D)

      D = s
      call SetTemp(D)

   end function RdFromDs

end module Rdiag_m
