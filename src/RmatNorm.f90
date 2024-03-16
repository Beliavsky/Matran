module RmatNorm_m
use MatranUtil_m
use Rmat_m
implicit none

#ifdef OVERVIEW

RmatNorm_m contains functions for computing various norms of a Rmat A.
Specifically,

   Normf(A)   = the Frobenius norm of A = sqrt(sum_{ij} a_{ij}^2)
   Norm1(A)   = the 1-norm of A         = max_j sum_i abs(a_{ij})
   Norminf(A) = the infinity norm of A  = max_i sum_j abs(a_{ij})

Note that the 2-norm (aka the spectral norm) is not included in this
suite, since it requires the computation of a spectral decomposition.
To get the 2-norm, use Norm2_m

Author: Che Rung Lee, Pete Stewart
Aug 21 2003

#endif

   interface Norm
      module procedure RmNormf
   end interface Norm

   interface Normf
      module procedure RmNormf
   end interface Normf

   interface Norm1
      module procedure RmNorm1
   end interface

   interface Norminf
      module procedure RmNorminf
   end interface

contains

   ! Compute the Frobenius norm of matrix A

   function RmNormf(A)
      
      real(wp) :: RmNormf
      type (Rmat), intent(in) :: A

      real(wp) :: amax, s
      integer :: i, j

      call GuardTemp(A)

      RmNormf = 0.0

      if (A%ncol==0 .or. A%nrow==0) return

      amax = 0.0

      do j=1,A%ncol
         do i=1,A%nrow
            amax = max(amax, abs(A%a(i,j)))
         end do
      end do

      if (amax == 0.0) goto 99999  ! Return

      s = 1.0/amax

      do j=1,A%ncol
         do i=1,A%nrow
            RmNormf = RmNormf + (s*A%a(i,j))**2
         end do
      end do

      RmNormf = sqrt(RmNormf)/s
      
99999 call CleanTemp(A)

   end function RmNormf


   ! Compute 1 norm of matrix A

   function RmNorm1(A) result(norm)
      real(wp)  :: norm
      type(Rmat):: A
      real(wp)  :: amax
      integer   :: j

      call GuardTemp(A)

      norm = 0.0
      do j = 1,A%ncol
         amax = sum(abs(A%a(1:A%nrow,j)))
         if(amax>norm) norm = amax
      end do

      call CleanTemp(A)
   end function RmNorm1


   ! Compute the infinity norm of a matrx

   function RmNormInf(A) result(norm)
      real(wp)  :: norm
      type(Rmat):: A

      real(wp) :: amax
      integer :: i

      call GuardTemp(A)
      norm = 0.0
      do i = 1,A%nrow
         amax = sum(abs(A%a(i,1:A%ncol)))
         if(amax>norm) norm = amax
      end do

      call CleanTemp(A)
   end function RmNormInf

end module RmatNorm_m

