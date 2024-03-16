program RmatRealSchurDebug

use MatranUtil_m
use RmatRealSchur_m
use RmatProduct_m
use RmatNorm_m
use RmatSum_m
use RmatPrint_m
use RmatEye_m
use RmatRand_m

implicit none

   type(Rmat)        :: A
   type(RmatRealSchur)   :: S
   integer           :: tn, n, i, i1, i2
   real(wp), pointer :: work(:)=>null()

   print *, ' '
   print *, "Enter test number"
   read *, tn, n
   print *, 'Test number ', tn
   print *, 'n = ', n
   

   A = RrandN(n)


   select case (tn)

   case(1)

      call Schur(S, A, .true., mywork=work)


      print *, normF((A*S%U)-(S%U*S%T))/normF(A), associated(work)
      print *, NormF(Reye(n) - .xhx.S%U)

      call clean(S)

   case(2)

      call Schur(S, A,.true.)
      print *, normF((A*S%U)-(S%U*S%T))/normF(A), associated(work)
      print *, NormF(Reye(n) - .xhx.S%U)

      print *, S%D

      print *, "Enter i1,i2"
      read *, i1,i2
      print *, 'i1 = ', i1, 'i2 = ', i2


      call ReorderSchur(S, i1, i2)
      print *, normF((A*S%U)-(S%U*S%T))/normF(A), associated(work)
      print *, NormF(Reye(n) - .xhx.S%U)
      print *, i1, i2
      print *, S%D

      call clean(S)

   case default

      print *, 'No such test.'

   end select

end program RmatRealSchurDebug
