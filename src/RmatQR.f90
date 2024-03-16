module RmatQR_m

use MatranUtil_m
use Rmat_m

implicit none

#ifdef OVERVIEW

Let A be an mxn matrix with m>=n.  Then there is an mxm orthogonal
matrix Q such that

         (R)
   A = Q*( )   (*)
         (0)

where R is an  nxn upper triangular matrix.  If Q = (Q1, Q2), then 

   A = Q1*R    (**)

If m<n, then the factorization takes the form

   A = QR      (***)

where R is an mxn upper triangular matrrix.

The module RmatQR_m uses the lapack routine DGEQRF to compute one of
the three above factorizations.  The result is returned in a defined
type RmatQR, defined below.  The generic subroutine for computing the
QR decomposition is

   QR(qr, A, fullq, mywork).

The arguments are explained below.

Author: Che Rung Lee
Jun 17 2003

#endif

   type RmatQR
      type(Rmat) :: Q                ! The Q-factor
      type(Rmat) :: R                ! The R-factor
      logical :: companion = .false. ! True if the decomposition is
                                     ! associated with a Rmat of interest
   end type RmatQR

   interface QR
      module procedure RmQR
   end interface

   interface Clean
      module procedure RmCleanQR
   end interface

contains

   ! Clean for type RmatQR

   subroutine RmCleanQR(qr)
      type(RmatQR), intent(inout) :: qr

      call Clean(qr%Q)
      call Clean(qr%R)

      qr%companion = .false.
   end subroutine RmCleanQR


   
   subroutine RmQR(qr, A, fullq, mywork)
      type(RmatQR), intent(out), target   :: QR
         ! The QR decmposition of A

      type(Rmat), intent(in) :: A
         ! The Rmat whose QR decomposition is to be computed

      logical, intent(in), optional :: fullq
         ! This controls whether a full or truncated Q factor
         ! is computed.  If fullq is absent or present and false,
         ! the factorization (**) or (***) is computed depending
         ! on the dimensions of A.  If fullq is present and true
         ! the factorization (*) or (**) is computed depending
         ! on the dimension of A.

      real(wp), pointer, optional:: mywork(:)
         ! Optional work array for dgeqrf.  If mywork is present
         ! it will be used, possibly after a reallocation to
         ! obtain enough memory.  Mywork is not deallocated
         ! by RmQR.

      integer :: icase, info, j, lwork, m, mm, n, nn, min_mn

      real(wp) :: twork(1), v(A%nrow), u(A%nrow), tau(min(A%nrow, A%ncol))
 
      real(wp), pointer :: work(:)=>null()

      logical :: usemywork, full

      real(wp), parameter :: ZERO = 0, ONE = 1

      type(Rmat), pointer ::  Q, R

      call GuardTemp(A)

      m = A%nrow
      n = A%ncol
      min_mn = min(m,n)

      ! Check if a full Q is wanted.

      full = .false.       
      if(present(fullq)) full = fullq

      ! Get working storage for dgeqrf.
#ifdef dbl
      call dgeqrf(m, n, A%a, A%narow, tau, twork, -1, info)
#endif
#ifdef sngl
      call sgeqrf(m, n, A%a, A%narow, tau, twork, -1, info)
#endif

      lwork = twork(1)

      if(present(mywork)) then
         usemywork = .true.
         call ReshapeAry(mywork, lwork)
         work => mywork
       else
          usemywork = .false.
          call ReshapeAry(work, lwork)
      end if

      v(1) = 1
      ! Decide the size of Q and R and copy A
      ! Case 1:  m>n and full --> Q(m,m), R(n,n), A->Q
      ! Case 2:  m>n but not full -> Q(m,n), R(n,n), A->Q
      ! Case 3:  m<=n -> Q(m,m), R(m,n), A->R

      Q => QR%Q
      R => QR%R
       
      if(full .or. m<n) then
         call ReshapeAry(Q, m, m)
         if(m>n) then
            icase = 1
            call ReshapeAry(R, n, n)
            Q%a(1:m,1:n) = A%a(1:m, 1:n)
         else
            icase = 3
            call ReshapeAry(R, m, n)
           R%a(1:m,1:n) = A%a(1:m, 1:n)
         end if
      else
         icase = 2
         call ReshapeAry(Q, m, n)
         call ReshapeAry(R, n, n)
         Q%a(1:m,1:n) = A%a(1:m, 1:n)
      end if

      R%tag = 'UT'
      Q%tag = 'GE'

  
      ! Compute theQR decomposition.

      select case (icase)

      case(1)  
                 
         ! result is in Q 
         ! Q is full (mxm)
         ! m > n
         
#ifdef dbl
         call dgeqrf(m, n, Q%a, Q%narow, tau, work, lwork, info)
#endif
#ifdef sngl
         call sgeqrf(m, n, Q%a, Q%narow, tau, work, lwork, info)
#endif

         do j = n+1, m
            Q%a(j,j) = 1
         end do      

         do j = n, 1, -1
            mm = m - j + 1

            R%a(1:j,j) = Q%a(1:j,j)
            v(2:mm)    = Q%a((j+1):m,j)
            Q%a(1:m,j) = 0
            Q%a(j,j)   = 1

            ! u(1:(m-j+1)) = Q(j:m, j:m)' * A(j:m,j)

#ifdef dbl
            call dgemv('T', mm, mm, ONE, Q%a(j:m, j:m), Q%narow-j+1, &
                        v(1:mm), 1, ZERO, u(1:mm), 1)
#endif
#ifdef sngl
            call sgemv('T', mm, mm, ONE, Q%a(j:m, j:m), Q%narow-j+1, &
                        v(1:mm), 1, ZERO, u(1:mm), 1)
#endif

            ! Q(j:m, j:m) += -tau(j) * A(j:m,j) * u(1:mm)

#ifdef dbl
            call dger(mm, mm,  -tau(j), v(1:mm), 1, &
                      u(1:mm), 1, Q%a(j:m, j:m), Q%narow-j+1)
#endif
#ifdef sngl
            call sger(mm, mm,  -tau(j), v(1:mm), 1, &
                      u(1:mm), 1, Q%a(j:m, j:m), Q%narow-j+1)
#endif

         end do

      
      case(2)  
             
         ! result is in Q
         ! Q is not full (mxn)
         ! m > n

#ifdef dbl
         call dgeqrf(m, n, Q%a, Q%narow, tau, work, lwork, info)
#endif
#ifdef sngl
         call sgeqrf(m, n, Q%a, Q%narow, tau, work, lwork, info)
#endif
          
         do j = n, 1, -1
         
            mm = m - j + 1
            nn = n - j + 1

            R%a(1:j,j) = Q%a(1:j,j)
            v(2:mm)    = Q%a((j+1):m,j)
            Q%a(1:m,j) = 0
            Q%a(j,j)   = 1
            

            ! u(1:(n-j+1)) = Q(j:m, j:n)' * A(j:m,j)

#ifdef dbl
            call dgemv('T', mm, nn, ONE, Q%a(j:m, j:n), Q%narow-j+1 , &
                        v(1:mm), 1, ZERO, u(1:nn), 1)
#endif
#ifdef sngl
            call sgemv('T', mm, nn, ONE, Q%a(j:m, j:n), Q%narow-j+1 , &
                        v(1:mm), 1, ZERO, u(1:nn), 1)
#endif

            ! Q(j:m, j:n) += -tau(j) * A(j:m,j) * u(1:(n-j+1))

#ifdef dbl
            call dger(mm, nn,  -tau(j), v(1:mm), 1, &
                      u(1:nn), 1, Q%a(j:m, j:n), Q%narow-j+1 )
#endif
#ifdef sngl
            call sger(mm, nn,  -tau(j), v(1:mm), 1, &
                      u(1:nn), 1, Q%a(j:m, j:n), Q%narow-j+1 )
#endif

         end do


      case(3)     
            
         ! result is in R
         ! m < n

#ifdef dbl
         call dgeqrf(m, n, R%a, R%narow, tau, work, lwork, info)
#endif
#ifdef sngl
         call sgeqrf(m, n, R%a, R%narow, tau, work, lwork, info)
#endif

         do j = m, 1, -1
         
            mm = m - j + 1

            v(2:mm) = R%a((j+1):m,j)
            R%a((j+1):m,j) = 0
            Q%a(j,j) = 1

            
            ! u(1:mm) = Q(j:m, j:m)' * A(j:m,j)

#ifdef dbl
            call dgemv('T', mm, mm, ONE, Q%a(j:m, j:m),  Q%narow-j+1, &
                        v(1:mm), 1, ZERO, u(1:mm), 1)
#endif
#ifdef sngl
            call sgemv('T', mm, mm, ONE, Q%a(j:m, j:m),  Q%narow-j+1, &
                        v(1:mm), 1, ZERO, u(1:mm), 1)
#endif

            ! Q(j:m, j:m) += -tau(j) * A(j:m,j) * u(1:mm)

#ifdef dbl
            call dger(mm, mm, -tau(j), v(1:mm), 1, &
                      u(1:mm), 1, Q%a(j:m, j:m),  Q%narow-j+1)
#endif
#ifdef sngl
            call sger(mm, mm, -tau(j), v(1:mm), 1, &
                      u(1:mm), 1, Q%a(j:m, j:m),  Q%narow-j+1)
#endif

         end do

      end select

      QR%companion = .true.

      ! Clean up.

      call CleanTemp(A)

      if(.not. usemywork) deallocate(work)

   end subroutine RmQR

end module RmatQR_m
