module MatranRealCore_m

   ! Root module

   use MatranUtil_m

   ! The two matrix objects

   use Rmat_m; use Rdiag_m

   ! Matrix operations

   use RmatTranspose_m; use RmatSum_m;  use RmatProduct_m
   use RmatSolve_m    ; use RmatJoin_m; use RmatBorder_m
   use RmatSubmatrix_m

   use RdiagSum_m; use RdiagProduct_m; use RdiagSolve_m

   ! Matrix Miscelania

   use RdiagDiag_m; use RmatEye_m;   use RmatNorm_m;
   use RmatPivot_m; use RmatPrint_m; use RmatRand_m

   ! Decompositions

   use RmatLudpp_m; use RmatChol_m

end module MatranRealCore_m
