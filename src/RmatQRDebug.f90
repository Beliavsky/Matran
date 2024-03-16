program RmatQRDebug

use MatranUtil_m
use RmatQR_m
use RmatProduct_m
use RmatNorm_m
use RmatSum_m
use RmatPrint_m
use RmatSubmatrix_m
use RmatEye_m
implicit none

   type(Rmat)   :: A, B
   type(RmatQr) :: dqr
   integer tn, i, j, m, n
   real(wp), pointer :: mywork(:)=>null()

   print *, ' '
   print *, "Enter test number"
   read *, tn
   print *, 'Test number ', tn


   select case (tn)

   case(1) 
      A = (/6,4/)
      do i=1,6
         do j=1,4
            A%a(i,j) = i+j
         end do
      end do

      call qr(dqr, A, .true.)
      B = col(dqr%Q, 1, 4)*dqr%R
      print *, NormF(Reye(dqr%Q%ncol) - .xhx.dqr%Q)
      print *, NormF(A - B)/normF(A)
      call clean(dqr)

   case(2)

      A = (/6,4/)
      do i=1,6
         do j=1,4
            A%a(i,j) = i+j
         end do
      end do

      call qr(dqr, A)
      B = dqr%Q*dqr%R
      print *, NormF(Reye(dqr%Q%ncol) - .xhx.dqr%Q)
      print *, NormF(A - B)/normF(A)
      call clean(dqr)

   case(3)

      A = (/4,6/)
      do i=1,4
         do j=1,6
            A%a(i,j) = i+j
         end do
      end do

      call qr(dqr, A)

      B = dqr%Q*dqr%R
      print *, NormF(Reye(dqr%Q%ncol) - .xhx.dqr%Q)
      print *, NormF(A - B)/normF(A)
      call clean(dqr)

  case(4)
      print *, 'Enter m and n'
      read *, m, n
      print *, 'm = ', m, 'n = ', n 

      A = (/m,n/)

      call Random_Number(A%a)

      call QR(dqr, A, .true., mywork)

      B = col(dqr%Q, 1, dqr%R%nrow)*dqr%R
      print *, NormF(Reye(dqr%Q%ncol) - .xhx.dqr%Q)
      print *, NormF(A - B)/normF(A)
      call clean(dqr)

   case default

      print *, 'No such test.'

  end select

end program RmatQrDebug
