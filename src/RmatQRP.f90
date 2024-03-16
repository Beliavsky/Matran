module RmatQRP_m

use MatranUtil_m
use Rmat_m

implicit none

#ifdef OVERVIEW
Let A be an mxn matrix with m>=n.  Then there is an mxm orthogonal
matrix Q and a permutation matrix P such that

          (R)
   AP = Q*( )   (*)
          (0)

where R is an  nxn upper triangular matrix.  If Q = (Q1, Q2), then 

   AP = Q1*R    (**)

If m<n, then the factorization takes the form

   AP = QR      (***)

where R is an mxn upper triangular matrrix.

The module RmatQR_m uses the lapack routine DGEQRF to compute one of
the three above factorizations.  The result is returned in a defined
type RmatQR, defined below.  The generic subroutine for computing the
QR decomposition is 

   QR(qr, A, fullq, firstcols, mywork).

The arguments are explained below.

Author: Che Rung Lee
Jun 17 2003

#endif

   ! QR decomposition with pivoting

   type RmatQRP
      type(Rmat)       :: Q        ! The Q-factor
      type(Rmat)       :: R        ! The R-factor
      integer, pointer &           ! The pivot array
               :: pvt(:)=>null()   !
      integer          :: npvt=0   ! the number of pivots
      logical :: companion=.false. ! True if the decomposition is
                                   ! associated with a Rmat of iterest
   end type RmatQRp

   ! The permutation P in (*), (**), and (**) is represented by
   ! interchages specified in the array pvt.  Specifically, for
   ! any row vector x, x*P can be computed as follows.
   !
   !   do i=1,npvt
   !      t = x(i)
   !      x(i) = x(pvt(i))
   !      x(pvt(i)) = t
   !   end do
   !
   !   For more see the Pivot suite.

   interface QRP
      module procedure RmQRP
   end interface

   interface Clean
      module procedure RmCleanQRP
   end interface

contains

   ! Clean for type RmatQRP

   subroutine RmCleanQRP(qr)
      type(RmatQRP), intent(inout) :: qr

      call Clean(qr%Q)
      call Clean(qr%R)
      if (associated(qr%pvt)) deallocate(qr%pvt)

      qr%companion = .false.
   end subroutine RmCleanQRP



   subroutine RmQRP(QR, A, fullq, firstcols, mywork)

      type(RmatQRP), intent(out),target  :: QR
         ! The QRP decmposition of A

      type(Rmat), intent(in) :: A
         ! The Rmat whose QRP decomposition is to be computed

      logical, intent(in), optional :: fullq
         ! This controls wether a full or truncated Q factor
         ! is computed.  If fullq is absent or present and false,
         ! the factorization (**) or (***) is computed depending
         ! on the dimensions of A.  If fullq is present and true
         ! the factorization (*) or (**) is computed depending
         ! on the dimension of A.


      integer, intent(in), optional :: firstcols(:)
         ! If present, QRP will swap all columns j for which
         ! firstcols(j)/=0 to the beginning of the matrix and
         ! freeze them there during the pivoting process.

      real(wp), pointer, optional :: mywork(:)
         ! Optional work array for dgeqrf.  If mywork is present
         ! it will be used, possibly after a reallocation to
         ! obtain enough memory.  Mywork is not deallocated
         ! by RmQRP.


      real(wp), parameter :: ZERO = 0, ONE = 1

      integer:: icase, info, lwork, m, n, min_mn

      real(wp):: twork(1)

      type(Rmat), pointer ::  Q, R

      real(wp), pointer:: work(:)

      logical :: full, usemywork

      integer :: i, j, mm, nn

      integer :: colord(A%ncol)

      real(wp) :: v(A%nrow), u(A%nrow), tau(min(A%ncol, A%nrow))


      call GuardTemp(A)

      m = A%nrow
      n = A%ncol

      min_mn = min(m,n)

      ! Check if a full Q is wanted.

      full = .false.       
      if (present(fullq)) full = fullq

      ! Set up the pivot array.

      colord = 0
      if (present(firstcols)) then
         do i=1, min(n, size(firstcols))
            colord(i) = firstcols(i)
         end do
      end if

      ! Get working storage for dgeqp3.

#ifdef dbl
      call dgeqp3(m, n, A%a, A%narow, colord, tau, twork, -1, info)
#endif
#ifdef sngl
      call sgeqp3(m, n, A%a, A%narow, colord, tau, twork, -1, info)
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
      ! Case 3:  m<n -> Q(m,m), R(m,n), A->R

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

      ! Compute the QRP decomposition.

      select case (icase)

      case(1)  
                 
         ! result is in Q 
         ! Q is full (mxm)
         ! m > n

#ifdef dbl
         call dgeqp3(m, n, Q%a, Q%narow, colord, tau, work, lwork, info)
#endif
#ifdef sngl
         call sgeqp3(m, n, Q%a, Q%narow, colord, tau, work, lwork, info)
#endif

         do j = n+1, m
            Q%a(j,j) = 1
         end do
      
         do j = n, 1, -1
         
            mm = m - j + 1

            R%a(1:j,j) = Q%a(1:j,j)
            v(2:mm)      = Q%a((j+1):m,j)
            Q%a(1:m,j) = 0
            Q%a(j,j)   = 1
            
            ! u(1:(m-j+1)) = Q(j:m, j:m)' * A(j:m,j)

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
            call dger(mm, mm,  -tau(j), v(1:mm), 1, &
                      u(1:mm), 1, Q%a(j:m, j:m),  Q%narow-j+1)
#endif
#ifdef sngl
            call sger(mm, mm,  -tau(j), v(1:mm), 1, &
                      u(1:mm), 1, Q%a(j:m, j:m),  Q%narow-j+1)
#endif

         end do

      
      case(2)  
             
         ! result is in Q
         ! Q is not full (mxn)
         ! m > n

#ifdef dbl
         call dgeqp3(m, n, Q%a, Q%narow, colord, tau, work, lwork, info)
#endif
#ifdef sngl
         call sgeqp3(m, n, Q%a, Q%narow, colord, tau, work, lwork, info)
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
            call dgemv('T', mm, nn, ONE, Q%a(j:m, j:n),  Q%narow-j+1, &
                        v(1:mm), 1, ZERO, u(1:nn), 1)
#endif
#ifdef sngl
            call sgemv('T', mm, nn, ONE, Q%a(j:m, j:n),  Q%narow-j+1, &
                        v(1:mm), 1, ZERO, u(1:nn), 1)
#endif

            ! Q(j:m, j:n) += -tau(j) * A(j:m,j) * u(1:(n-j+1))

#ifdef dbl
            call dger(mm, nn,  -tau(j), v(1:mm), 1, &
                      u(1:nn), 1, Q%a(j:m, j:n),  Q%narow-j+1)
#endif
#ifdef sngl
            call sger(mm, nn,  -tau(j), v(1:mm), 1, &
                      u(1:nn), 1, Q%a(j:m, j:n),  Q%narow-j+1)
#endif

         end do


      case(3)     
            
         ! result is in R
         ! m < n

#ifdef dbl
         call dgeqp3(m, n, R%a, R%narow, colord, tau, work, lwork, info)
#endif
#ifdef sngl
         call sgeqp3(m, n, R%a, R%narow, colord, tau, work, lwork, info)
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


      ! Change the permutation array to an exchange array

      QR%npvt = min_mn

      call ReshapeAry(QR%pvt,min_mn)

      do j=1,min_mn
         QR%pvt(j) = colord(j)
         i = j
         do while (colord(i) /= j)
            i = colord(i)
         end do
         colord(i) = colord(j)
      end do

      ! Clean up.

      call CleanTemp(A)
      if(.not. usemywork) deallocate(work)
      QR%companion = .true.

   end subroutine RmQRP

end module RmatQRP_m
