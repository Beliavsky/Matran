program RmatSvdDebug

use MatranUtil_m
use Rmat_m
use RdiagProduct_m
use RmatSvd_m
use RmatProduct_m
use RmatNorm_m
use RmatSum_m
use RmatRand_m
use RmatSubmatrix_m
use RmatEye_m

implicit none

   type(Rmat)    :: A
   type(RmatSvd) :: S
   integer       :: m, n
   real(wp), pointer :: mywork(:)=>null()
   
   print *, 'Enter m and n'
   read *, m, n
   print *, 'm = ', m, 'n = ', n

   A = RrandN(m,n)

   call svd(S, A, .true., .true.)

   print *, normF(A - S%U*S%D*.ctp.S%V)/normF(A)
   if (m>n) then
      print *, NormF(Reye(n) - .xhx.S%U)
      print *, NormF(Reye(n) - .xhx.S%V)
   else
      print *, NormF(Reye(m) - .xhx.S%U)
      print *, NormF(Reye(m) - .xhx.S%V)
   end if

   call svd(S, A, .true., .true., .true.)

   if (m > n) then
!      print *, NormF(A - sbm(S%U,1,m,1,n)*S%D*.ctp.S%V)
      print *, NormF(A - .rm.S%U%a(1:m,1:n)*S%D*.ctp.S%V)
   else
      print *, NormF(A - S%U*S%D*.ctp.sbm(S%V,1,n,1,m))
   end if
   print *, NormF(Reye(m) - .xhx.S%U)
   print *, NormF(Reye(n) - .xhx.S%V)

   call svd(S, A, .true., .true., mywork = mywork)
   print *, normF(A - S%U*S%D*.ctp.S%V)/normF(A)
   if (m>n) then
      print *, NormF(Reye(n) - .xhx.S%U)
      print *, NormF(Reye(n) - .xhx.S%V)
   else
      print *, NormF(Reye(m) - .xhx.S%U)
      print *, NormF(Reye(m) - .xhx.S%V)
   end if

   call svd(S, A, .true., .true., mywork = mywork)
   print *, normF(A - S%U*S%D*.ctp.S%V)/normF(A)
   if (m>n) then
      print *, NormF(Reye(n) - .xhx.S%U)
      print *, NormF(Reye(n) - .xhx.S%V)
   else
      print *, NormF(Reye(m) - .xhx.S%U)
      print *, NormF(Reye(m) - .xhx.S%V)
   end if


end program RmatSvdDebug
