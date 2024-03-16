module RmatRand_m
use MatranUtil_m
use Rmat_m
implicit none

#ifdef OVERVIEW

RmatRand_m produces matrices with random entries with distribution X,
where

   X = U: Uniformally distributed on [0,1]
   X = N: Normally distrubed (0,1)

The subroutine form has two forms.

   call RandX(A, m)
      On return A is an mxm random Rmat.

   call RandX(A, m, n)
      On return A is an mxn random Rmat.

There are also functional forms.

   RrandX(m)
      Returns an mxm random matrix.

   RrandX(m,n)
      Returns an mxn random matrix.

Author: Che Rung Lee
Jun 17 2003

#endif

   interface RandU
      module procedure RmRandU
   end interface

   interface RandN
      module procedure RmRandN
   end interface

   interface RrandU
      module procedure RrandU
   end interface

   interface RrandN
      module procedure RrandN
   end interface

contains
      
   ! Sets A to a uniformly distributed random matrix.

   subroutine RmRandU(A, m, n)
      type(Rmat), intent(inout) :: A
      integer, intent(in) :: m
      integer, optional, intent(in) :: n

      integer :: mm, nn

      mm = m
      if(present(n)) then
         nn = n
      else
         nn = m
      end if

      if (mm < 0 .or. nn < 0)&
         call MatranError("RmRandU in RmatRand: Negative dimension.")


      call ReshapeAry(A, mm, nn)
      call Random_NUmber(A%a(1:mm, 1:nn))

      A%tag = 'GE'
         
   end subroutine RmRandU

   ! Returns a uniformally distributed random matrix

   function RrandU(m, n) result(A)

      type(Rmat) :: A
      integer :: m
      integer, optional :: n
      integer nn

      A%a => null()
      A%temporary => null()
      call Clean(A)

      if (present(n)) then
         nn = n
      else
         nn = m
      end if

      call RandU(A, m, nn)

      call SetTemp(A)

   end function RrandU

   ! Sets A to a normally distributed random matrix.

   subroutine RmRandN(A, m, n)
      type(Rmat), intent(inout) :: A
      integer, intent(in) :: m
      integer, optional, intent(in) :: n

      integer :: i, j, mm, nn

      mm = m
      if(present(n)) then
         nn = n
      else
         nn = m
      end if

      if (mm < 0 .or. nn < 0)&
         call MatranError("RmRandN in RmatRand: Negative dimension.")


      call ReshapeAry(A, mm, nn)

      do i=1, A%nrow
         do j = 1, A%ncol
            A%a(i,j) = RmNormal()
         end do
      end do

      A%tag = 'GE'
         
   end subroutine RmRandN

   ! Returns a normally distributed random matrix.

   function RrandN(m, n) result(A)

      type(Rmat) :: A
      integer :: m
      integer, optional :: n
      integer nn

      A%a => null()
      A%temporary => null()
      call Clean(A)

      if (present(n)) then
         nn = n
      else
         nn = m
      end if

      call RandN(A, m, nn)

      call SetTemp(A)

   end function RrandN



   ! Generates a psuedo-random normal number.
   ! Algorithm due to Leva [TOMS 18 (1992) 454-455].

   function RmNormal() result(N)
      real(wp):: S= 0.449871, T = -0.386595
      real(wp):: A= 0.19600, B = 0.25472
      real(wp):: R1= 0.27597, R2 = 0.27846

      real U, V, X, Y, Q, N

      do 
         call Random_Number(U)
         call Random_Number(V)
         V = 1.7156 * (V - 0.5)
         X  = U -  S
         Y  = ABS(V) - T
         Q  = X**2 + Y*(A*Y - B*X)
         if (Q .LT. R1) exit
         if(.NOT. (Q .GT. R2) .AND. .NOT. &
              (V**2 .GT. -4.0*LOG(U)*U**2)) exit
      end do
      N = V/U
      
   end function RmNormal

      

end module RmatRand_m
