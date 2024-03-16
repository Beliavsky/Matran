program RmatSpecDebug

use MatranUtil_m
use Rmat_m
use RmatSpec_m
use RmatProduct_m
use RmatNorm_m
use RmatSum_m
use RdiagProduct_m
use RmatTranspose_m

implicit none


   type(Rmat)     :: A, B, C
   type(RmatSpec) :: E
   integer        :: i, n, tn
   real(wp), pointer :: mywork(:)=>null()

   print *, ' '
   print *, "Enter test number"
   read *, tn
   print *, 'Test number ', tn


   select case (tn)

   case(1)

      print *, 'Enter n.'
      read *, n
      print *, 'n = ', n

      A = (/n, n/)

      call Random_Number(A%a)
   
      A = A + .ctp.A

      A%tag = 'HE'

      call Spec(E, A, .true.)


      print *, NormF(A*E%V - E%V*E%D)/NormF(A)

      call Clean(E)

   case(2)

      print *, 'Enter n.'
      read *, n
      print *, 'n = ', n

      A = (/n, n/)

      call Random_Number(A%a)
   
      A = A + .ctp.A

      A%tag = 'HE'

      call Spec(E, A, .true., mywork=mywork)

      print *, NormF(A*E%V - E%V*E%D)/NormF(A)

      call Clean(E)

   case(3)

      print *, 'Enter n.'
      read *, n
      print *, 'n = ', n

      A = (/n, n/)

      call Random_Number(A%a)
   
      B = .xhx.A

      B%tag = 'HE'

      call Spec(E, .xhx.A, .true., mywork=mywork)

      print *, NormF(B*E%V - E%V*E%D)/NormF(B)

      call Clean(E)

   case default

      print *, 'No such test.'

   end select

end program RmatSpecDebug
