program RmatEigDebug

use MatranUtil_m
use RmatEig_m
use RmatProduct_m
use RmatNorm_m
use RmatSum_m
use Rdiag_m
use RdiagProduct_m
use RmatSolve_m

implicit none

   type(Rmat)        :: A, B, C, D
   type(RmatEig)     :: E
   type(Rdiag)       :: lam
   complex(wp)       :: tmp
   real(wp)          :: sum
   real(wp), pointer :: wv(:,:), mywork(:)
   integer           :: tn, i, j, k, n


   print *, ' '
   print *, "Enter test number"
   read *, tn
   print *, 'Test number ', tn


   select case (tn)

   case(1) 
      print *, "Enter matrix size "
      read *, n
      print *, 'n = ', n

      A = (/n, n/)
      call Random_Number(A%a)
      A%a = A%a - 0.5D0

      call Eig(E, A, .true.)

      ! compute AX - XEig

      sum = 0.0D0
      do i = 1,n
         do j = 1, n
            tmp = (0.0D0, 0.0D0)
            do k = 1, n
               tmp = tmp + A%a(i,k) * E%X(k,j)
            end do 
            tmp = tmp - E%X(i,j)*E%D(j)
            sum = sum + tmp*conjg(tmp)
         end do
      end do

      print *, sqrt(sum)/normF(A)
      call clean(E)

   case(2)
      print *, "Enter matrix size "
      read *, n
      print *, 'n = ', n

      A = (/n, n/)
      call Random_Number(A%a)
      A%a = A%a - 0.5D0

      call eig(E, A, .false., .true.)

      ! compute X'A - EigX'

      sum = 0.0D0
      do i = 1, n
         do j = 1, n
            tmp = (0.0D0, 0.0D0)
            do k = 1, n
               tmp = tmp + A%a(k,j) * conjg(E%Y(k,i))
            end do 
            tmp = tmp - conjg(E%Y(j,i))*E%D(i)
            sum = sum + tmp*conjg(tmp)
         end do
      end do

      print *, sqrt(sum)/normF(A)
      call clean(E)

   case(3) 
     print *, "Enter matrix size "
      read *, n
      print *, 'n = ', n

      A = (/n, n/)
      call Random_Number(A%a)
      A%a = A%a - 0.5D0

      call Eig(E, A, .true., wrv=wv, mywork=mywork )

      ! compute AX - XEig

      sum = 0.0D0
      do i = 1,n
         do j = 1, n
            tmp = (0.0D0, 0.0D0)
            do k = 1, n
               tmp = tmp + A%a(i,k) * E%X(k,j)
            end do 
            tmp = tmp - E%X(i,j)*E%D(j)
            sum = sum + tmp*conjg(tmp)
         end do
      end do

      print *, sqrt(sum)/normF(A)
      call clean(E)

   case(4)
      print *, "Enter matrix size "
      read *, n
      print *, 'n = ', n

      A = (/n, n/)
      call Random_Number(A%a)
      A%a = A%a - 0.5D0

      call Eig(E, A, .false., .true., wlv=wv, mywork=mywork )

      ! compute X'A - EigX'

      sum = 0.0D0
      do i = 1, n
         do j = 1, n
            tmp = (0.0D0, 0.0D0)
            do k = 1, n
               tmp = tmp + A%a(k,j) * conjg(E%Y(k,i))
            end do 
            tmp = tmp - conjg(E%Y(j,i))*E%D(i)
            sum = sum + tmp*conjg(tmp)
         end do
      end do

      print *, sqrt(sum)/normF(A)
      call clean(E)

   case(5)

      C = (/5,5/)
      call random_number(C%a)
      lam = (/5/)
      do i=1,5
         lam%a(i) = i
      end do
      A = C.xiy.(lam*C)

      call Eig(E, A)
      print *, E%D
      call clean(E)

   case default

      print *, 'No such test.'

   end select

end program RmatEigDebug
