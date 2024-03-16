module RmatBorder_m

use MatranUtil_m
use Rmat_m

implicit none

#ifdef OVERVIEW

RmatBorder enlarges a Rmat by adding other Rmats into
its border. The borders of a matrix are specified by their 
directions: East, West, North, South, NorthEast, NorthWest, 
SouthEast and SouthWest of the matrix.

   BorderE
     subroutine RmBorderE(A, B): A = [A, B]     
   BorderS
     subroutine RmBorderS(A, B): A = [A; B]
   BorderW
     subroutine RmBorderW(A, B): A = [B, A]     
   BorderN
     subroutine RmBorderN(A, B): A = [B; A]
   BorderSE
     subroutine RmBorderSE(A, S, E, SE): A = [A, E; S, SE]
   BorderSW
     subroutine RmBorderSE(A, S, W, SW): A = [W, A; SW, S]
   BorderNE
     subroutine RmBorderSE(A, N, E, NE): A = [N, NE; A, E]
   BorderNW
     subroutine RmBorderSE(A, N, W, NW): A = [NW, N; W, A]

These routines make use of a generic subroutine

   Enlarge(A, m, n, pos)

That increases the size of its array without destroying the matrix.

Author: Che-Rung Lee
Jun 17 2003

#endif

   interface Enlarge
      module procedure RmEnlarge
   end interface

   interface BorderE
      module procedure RmBorderE
   end interface

   interface BorderS
      module procedure RmBorderS
   end interface

   interface BorderW
      module procedure RmBorderW
   end interface

   interface BorderN
      module procedure RmBorderN
   end interface

   interface BorderNE
      module procedure RmBorderNE
   end interface

   interface BorderSE
      module procedure RmBorderSE
   end interface

   interface BorderSW
      module procedure RmBorderSW
   end interface

   interface BorderNW
      module procedure RmBorderNW
   end interface

contains

   ! Enlarge a matrix A without destroying its content

   subroutine RmEnlarge(A, m, n, pos)
      type(Rmat)        :: A
      integer           :: m, n   ! the dim of A to be enlarged
      character(2)      :: pos    ! the directions of A to be enlarged

      real(wp), pointer :: T(:,:) ! temp var to hold the data of A
      integer           :: ma, na

      if (.not.associated(A%a)) then
         call ReshapeAry(A, m, n)
      else
         if (m<A%nrow .or. n<A%ncol) &
              call MatranError("RmEnlarge in RmatBorder: array &
                                &cannot fit into the specified dims.")

         if (m>A%narow .or. n>A%nacol) then
            
            ! A does not have enough space.
            ! Allocate new space if possible.

            ma = A%nrow
            na = A%ncol

            if (A%adjustable) then
               allocate(T(m,n))
               T = 0.0

               ! copy A's data to T.

               select case (pos)
               case('SE')
                  T(1:ma, 1:na) = A%a(1:ma, 1:na)
               case('SW')
                  T(1:ma, n-na+1:n) = A%a(1:ma, 1:na)
               case('NE')
                  T(m-ma+1:m, 1:na) = A%a(1:ma, 1:na)
               case('NW')
                  T(m-ma+1:m, n-na+1:n) = A%a(1:ma, 1:na)
               end select

               deallocate(A%a)
               A%a => T
               A%narow = m
               A%nacol = n
            else
               call MatranError("RmEnlarge in RmatBorder: Cannot &
                            &change array in an unadjustable Rmat.")
            end if
         else 

            ! A has enough space. Move A's data to the location
            ! specified by pos.

            select case (pos)
            case('SE')
            case('SW')
               A%a(1:ma,na+1:n) = A%a(1:ma, 1:na)
            case('NE')
               A%a(m-ma+1:m, 1:n) = A%a(1:ma, 1:na)
            case('NW')
               A%a(m-ma+1:m, n-na+1:n) = A%a(1:ma, 1:na)
            end select
         end if

         A%nrow = m
         A%ncol = n
      end if
   end subroutine RmEnlarge

   !  A = [A, B]
   
   subroutine RmBorderE(A, B)
      type(Rmat), intent(inout) :: A
      type(Rmat), intent(in)    :: B
      integer :: m, n, nn

      call GuardTemp(B)

      if (A%nrow /= B%nrow) &
         call MatranError('RmBorderE in RmatBorder: &
                               &Incompatible dimensions.')

      m = A%nrow 
      n = A%ncol
      nn = A%ncol + B%ncol

      call RmEnlarge(A, m, nn, 'SE')
      A%a(1:m, n+1:nn) = B%a(1:m,1:B%ncol)

      ! Adjust A's type
      if (A%tag /= 'UT') A%tag = 'GE'

      call CleanTemp(B)
   end subroutine RmBorderE

   ! A = [B, A]
    
   subroutine RmBorderW(A, B)
      type(Rmat), intent(inout) :: A
      type(Rmat), intent(in)    :: B

      integer           :: n, m, nn
      
      call GuardTemp(B)

      if (A%nrow /= B%nrow) &
         call MatranError('RmBorderW in RmatBorder: &
                               &Incompatible dimensions.')

      m = A%nrow
      n = A%ncol
      nn = A%ncol+ B%ncol
   
      call RmEnlarge(A, m, nn, 'SW')
      A%a(1:m, 1:B%ncol) = B%a(1:m, 1:B%ncol)

      ! Adjust the type.
 
      A%tag = 'GE'
      if(B%tag == 'UT') A%tag = 'UT'

      call CleanTemp(B)
   end subroutine RmBorderW


   !  A = [A; B]
   
   subroutine RmBorderS(A, B)
      type(Rmat), intent(inout) :: A
      type(Rmat), intent(in)    :: B

      integer           :: m, n, mm
      
      call GuardTemp(B)

      if (A%ncol /= B%ncol) &
         call MatranError('RmBorderS in RmatBorder: &
                               &Incompatible dimensions.')

      m = A%nrow
      mm = A%nrow + B%nrow
      n = A%ncol

      call RmEnlarge(A, mm, n, 'SE')
      A%a(m+1:mm, 1:n) = B%a(1:B%nrow, 1:n)

      ! Adjust the type.
      if (A%tag /= 'LT') A%tag = 'GE'

      call CleanTemp(B)
   end subroutine RmBorderS


   ! A = [B; A]
    
   subroutine RmBorderN(A, B)
      type(Rmat), intent(inout) :: A
      type(Rmat), intent(in)    :: B
      integer           :: n, m
      
      call GuardTemp(B)

      if (A%ncol /= B%ncol) &
         call MatranError('RmBorderN in RmatBorder: &
                               &Incompatible dimensions.')
      m = A%nrow + B%nrow
      n = A%ncol

      call RmEnlarge(A, m, n, 'NE')
      A%a(1:B%nrow, 1:n) = B%a(1:B%nrow, 1:n)

      ! Adjust the type.

      A%tag = 'GE'
      if (B%tag == 'LT') A%tag = 'LT'

      call CleanTemp(B)
   end subroutine RmBorderN

   ! A = [A E; S SE]

   subroutine RmBorderSE(A, S, E, SE)
      type(Rmat), intent(inout) :: A
      type(Rmat), intent(in)    :: S
      type(Rmat), intent(in)    :: E
      type(Rmat), intent(in)    :: SE
      integer :: m1, n1, m2, n2, m3, n3

      call GuardTemp(S)
      call GuardTemp(E)
      call GuardTemp(SE)

      m1 = A%nrow
      n1 = A%ncol
      m2 = SE%nrow
      n2 = SE%ncol
      
      if (S%ncol/=n1 .or. S%nrow/=m2 .or. &
          E%ncol/=n2 .or. E%nrow/=m1) &
         call MatranError('RmBorderSE in RmatBorder: &
                               &Incompatible dimensions.')


      m3  = m1 + m2
      n3  = n1 + n2

      call RmEnlarge(A, m3, n3, 'SE')
      A%a(1:m1, n1+1:n3) = E%a(1:m1, 1:n2)
      A%a(m1+1:m3, 1:n1)  = S%a(1:m2, 1:n1)
      A%a(m1+1:m3, n1+1:n3) = SE%a(1:m2, 1: n2)
      
      ! Adjust the type.
      A%tag = 'GE'

      call CleanTemp(S)
      call CleanTemp(E)
      call CleanTemp(SE)

   end subroutine RmBorderSE


   ! A = [W A; SW S]

   subroutine RmBorderSW(A, S, W, SW)
      type(Rmat), intent(inout) :: A
      type(Rmat), intent(in)    :: S
      type(Rmat), intent(in)    :: W
      type(Rmat), intent(in)    :: SW
      integer :: m3, n3, m1, n1, m2, n2

      call GuardTemp(S)
      call GuardTemp(W)
      call GuardTemp(SW)

      m1 = A%nrow
      n1 = A%ncol
      m2 = SW%nrow
      n2 = SW%ncol
      
      if (S%ncol/=n1 .or. S%nrow/=m2 .or. &
          W%ncol/=n2 .or. W%nrow/=m1) &
         call MatranError('RmBorderSW in RmatBorder: &
                               &Incompatible dimensions.')      
      m3  = m1 + m2
      n3  = n1 + n2

      call RmEnlarge(A, m3, n3, 'SW')
      A%a(1:m1, 1:n2) = W%a(1:m1, 1:n2)
      A%a(m1+1:m3, n2+1:n3) = S%a(1:m2, 1:n1)
      A%a(m1+1:m3, 1:n2) = SW%a(1:m2, 1:n2)
      
      ! Adjust the type.
      A%tag = 'GE'

      call CleanTemp(S)
      call CleanTemp(W)
      call CleanTemp(SW)
   end subroutine RmBorderSW

   ! A = [N NE; A E]

   subroutine RmBorderNE(A, N, E, NE)
      type(Rmat), intent(inout) :: A
      type(Rmat), intent(in)    :: N
      type(Rmat), intent(in)    :: E
      type(Rmat), intent(in)    :: NE
      integer :: m3, n3, m1, n1, m2, n2

      call GuardTemp(N)
      call GuardTemp(E)
      call GuardTemp(NE)

      m1 = A%nrow
      n1 = A%ncol
      m2 = NE%nrow
      n2 = NE%ncol
      
      if (N%ncol/=n1 .or. N%nrow/=m2 .or. &
          E%ncol/=n2 .or. E%nrow/=m1) &
         call MatranError('RmBorderNE in RmatBorder: &
                               &Incompatible dimensions.')

      m3  = m1 + m2
      n3  = n1 + n2

      call RmEnlarge(A, m3, n3, 'NE')
      A%a(1:m2, 1:n1) = N%a(1:m2, 1:n1)
      A%a(m2+1:m3, n1+1:n3) = E%a(1:m1, 1:n2)
      A%a(1:m2, n1+1:n3) = NE%a(1:m2, 1:n2)
      
      ! Adjust the type.
      A%tag = 'GE'


      call CleanTemp(N)
      call CleanTemp(E)
      call CleanTemp(NE)
   end subroutine RmBorderNE


   ! A = [W A; SW S]

   subroutine RmBorderNW(A, N, W, NW)
      type(Rmat), intent(inout) :: A
      type(Rmat), intent(in)    :: N
      type(Rmat), intent(in)    :: W
      type(Rmat), intent(in)    :: NW
      integer :: m3, n3, m1, n1, m2, n2

      call GuardTemp(N)
      call GuardTemp(W)
      call GuardTemp(NW)

      m1 = A%nrow
      n1 = A%ncol
      m2 = NW%nrow
      n2 = NW%ncol
      
      if (N%ncol/=n1 .or. N%nrow/=m2 .or. &
          W%ncol/=n2 .or. W%nrow/=m1) &
         call MatranError('RmBorderNW in RmatBorder: &
                               &Incompatible dimensions.')
      m3 = m1 + m2
      n3 = n1 + n2

      call RmEnlarge(A, m3, n3, 'NW')
      A%a(m2+1:m3, 1:n2) = W%a(1:m1, 1:n2)
      A%a(1:m2, n2+1:n3) = N%a(1:m2, 1:n1)
      A%a(1:m2, 1:n2) = NW%a(1:m2, 1:n2)
      
      ! Adjust the type.
      A%tag = 'GE'

      call CleanTemp(N)
      call CleanTemp(W)
      call CleanTemp(NW)
   end subroutine RmBorderNW

end module RmatBorder_m
