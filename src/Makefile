.SUFFIXES : .o .f90

FC = gfortran
PREC = dbl
LIB = liblapack.a librefblas.a
#LIB = lapack.a blas.a

.f95.o :
	$(FC) -c -cpp -D$(PREC) $<

objs = \
Rdiag.o \
RdiagDiag.o \
RdiagJoin.o \
RdiagProduct.o \
RdiagSolve.o \
RdiagSum.o \
Rmat.o \
RmatBorder.o \
RmatChol.o \
RmatEig.o \
RmatEye.o \
RmatInv.o \
RmatJoin.o \
RmatLudpp.o \
RmatNorm.o \
RmatNorm2.o \
RmatPivot.o \
RmatPrint.o \
RmatProduct.o \
RmatQR.o \
RmatQRP.o \
RmatRand.o \
RmatRealSchur.o \
RmatSVD.o \
RmatSolve.o \
RmatSpec.o \
RmatSubmatrix.o \
RmatSum.o \
RmatTranspose.o \
MatranRealCore.o \
MatranUtil.o



all: $(objs)

clean :
	rm -f *.o *.mod  Matran.a core

# Make chsuf changes all files of the form <file>.s1 to <file>.s2

chsuf :
	ls -d *.s1 | sed 's/\(.*\)\.s1$$/mv & \1.s2/' | sh

arch :
	rm -f Matran.a; ar -r Matran.a *.o

MatranUtilDebug.o : MatranUtil.o
matranutildebug : MatranUtilDebug.o
	$(FC) MatranUtil.o MatranUtilDebug.o

RmatPrint.o: Rmat.o MatranUtil.o
RmatPrintDebug.o: RmatPrint.o Rmat.o
rmatprintdebug: RmatPrintDebug.o
	$(FC) MatranUtil.o Rmat.o RmatPrint.o RmatPrintDebug.o 

Rmat.o : MatranUtil.o
RmatDebug.o : Rmat.o RmatPrint.o MatranUtil.o
rmatdebug : RmatDebug.o
	$(FC) MatranUtil.o Rmat.o RmatPrint.o RmatDebug.o 


RmatSum.o : Rmat.o MatranUtil.o
RmatSumDebug.o : RmatSum.o RmatPrint.o Rmat.o
rmatsumdebug : RmatSumDebug.o
	$(FC) MatranUtil.o Rmat.o RmatSum.o RmatSumDebug.o RmatPrint.o 


RmatProduct.o : Rmat.o MatranUtil.o
RmatProductDebug.o : RmatProduct.o RmatPrint.o
rmatproductdebug : RmatProductDebug.o
	$(FC) MatranUtil.o Rmat.o RmatProduct.o RmatProductDebug.o\
              RmatPrint.o $(LIB) 

RmatPivot.o : Rmat.o MatranUtil.o
RmatPivotDebug.o : RmatPivot.o Rmat.o RmatPivot.o RmatPrint.o
rmatpivotdebug: RmatPivotDebug.o
	$(FC) MatranUtil.o Rmat.o RmatPivot.o RmatPivotDebug.o\
              RmatPrint.o $(LIB)

RmatNorm.o : Rmat.o  MatranUtil.o
RmatNormDebug.o : RmatNorm.o Rmat.o
rmatnormdebug: RmatNormDebug.o
	$(FC) MatranUtil.o Rmat.o RmatNorm.o RmatNormDebug.o 

RmatLudpp.o : Rmat.o
RmatLudppDebug.o :: RmatProduct.o RmatSum.o RmatNorm.o RmatPivot.o\
                   RmatLudpp.o RmatPrint.o
rmatludppdebug: RmatLudppDebug.o
	$(FC) MatranUtil.o Rmat.o RmatProduct.o RmatSum.o RmatNorm.o\
              RmatPivot.o RmatLudpp.o RmatPrint.o\
              RmatLudppDebug.o $(LIB) 

RmatChol.o : Rmat.o MatranUtil.o
RmatCholDebug.o : Rmat.o RmatProduct.o RmatSum.o RmatNorm.o RmatChol.o\
                 RmatPrint.o
rmatcholdebug: RmatCholDebug.o
	$(FC) MatranUtil.o Rmat.o RmatProduct.o RmatSum.o RmatNorm.o\
              RmatChol.o RmatCholDebug.o RmatPrint.o $(LIB) 

RmatSolve.o : RmatChol.o RmatLudpp.o RmatPivot.o
RmatSolveDebug.o: RmatProduct.o RmatSum.o RmatNorm.o RmatSolve.o RmatPivot.o\
                  Rmat.o RmatPrint.o MatranUtil.o
rmatsolvedebug: RmatSolveDebug.o
	$(FC) MatranUtil.o Rmat.o RmatProduct.o RmatSum.o RmatNorm.o\
              RmatChol.o RmatLudpp.o RmatSolve.o RmatSolveDebug.o\
              RmatPivot.o RmatPrint.o $(LIB) 

RmatTranspose.o : Rmat.o
RmatTransposeDebug.o: RmatTranspose.o Rmat.o RmatPrint.o
rmattransposedebug: RmatTransposeDebug.o
	$(FC) MatranUtil.o Rmat.o RmatTranspose.o RmatTransposeDebug.o\
              RmatPrint.o 

RmatJoin.o: Rmat.o MatranUtil.o
RmatJoinDebug.o: RmatPrint.o RmatJoin.o RmatSum.o RmatNorm.o Rmat.o \
                 MatranUtil.o
rmatjoindebug: RmatJoinDebug.o
	$(FC) MatranUtil.o Rmat.o RmatPrint.o RmatJoin.o \
              RmatSum.o RmatNorm.o RmatJoinDebug.o

RmatSubmatrix.o:      Rmat.o MatranUtil.o
RmatSubmatrixDebug.o: RmatSubmatrix.o Rmat.o RmatPrint.o MatranUtil.o
rmatsubmatrixdebug:   RmatSubmatrixDebug.o
	$(FC) MatranUtil.o Rmat.o RmatPrint.o  \
              RmatSubmatrixDebug.o RmatSubmatrix.o

RmatBorder.o:      Rmat.o MatranUtil.o
RmatBorderDebug.o: RmatSubmatrix.o Rmat.o RmatJoin.o RmatSum.o RmatNorm.o \
                   RmatBorder.o RmatPrint.o MatranUtil.o
rmatborderdebug: RmatBorderDebug.o
	$(FC) MatranUtil.o Rmat.o RmatBorder.o RmatJoin.o \
              RmatBorderDebug.o RmatSubmatrix.o RmatSum.o RmatNorm.o \
              RmatPrint.o


RmatQR.o :     Rmat.o MatranUtil.o
RmatQRDebug.o : RmatProduct.o RmatNorm.o RmatSum.o RmatQR.o RmatPrint.o \
                RmatSubmatrix.o RmatEye.o
rmatqrdebug:    RmatQRDebug.o 
	$(FC)  MatranUtil.o Rmat.o RmatQR.o RmatProduct.o RmatSum.o \
               RmatNorm.o RmatQRDebug.o RmatPrint.o RmatSubmatrix.o \
               RmatEye.o $(LIB)

RmatQRP.o :      Rmat.o MatranUtil.o
RmatQRPDebug.o : RmatProduct.o RmatNorm.o RmatSum.o RmatQRP.o RmatPrint.o \
                 RmatSubmatrix.o RmatEye.o RmatPivot.o MatranUtil.o
rmatqrpdebug:    RmatQRPDebug.o 
	$(FC)    MatranUtil.o Rmat.o RmatQRP.o RmatProduct.o RmatSum.o \
                 RmatNorm.o RmatQRPDebug.o RmatPrint.o RmatSubmatrix.o \
                 RmatEye.o RmatPivot.o $(LIB)

RmatEig.o :     Rmat.o MatranUtil.o
RmatEigDebug.o : RmatProduct.o RmatNorm.o RmatSum.o RmatEig.o Rdiag.o \
                 RdiagProduct.o RmatSolve.o MatranUtil.o
rmateigdebug:    RmatEigDebug.o 
	$(FC)   MatranUtil.o Rmat.o RmatEig.o RmatProduct.o RmatSum.o \
                RmatNorm.o RmatEigDebug.o Rdiag.o RdiagProduct.o \
                RmatSolve.o RmatChol.o RmatLudpp.o RmatPivot.o \
                $(LIB)

RmatRealSchur.o :     Rmat.o MatranUtil.o
RmatRealSchurDebug.o : RmatProduct.o RmatNorm.o RmatSum.o RmatRealSchur.o \
                       RmatEye.o RmatSolve.o RmatPrint.o RmatRand.o
rmatrealschurdebug:    RmatRealSchurDebug.o  MatranUtil.o
	$(FC)   MatranUtil.o Rmat.o RmatRealSchur.o RmatProduct.o RmatSum.o \
                RmatNorm.o RmatRealSchurDebug.o RmatPrint.o RmatRand.o\
                RmatEye.o $(LIB)

RmatSpec.o :    Rmat.o Rdiag.o MatranUtil.o
RmatSpecDebug.o: RmatProduct.o RmatNorm.o RmatSum.o RmatSpec.o Rdiag.o \
                 Rmat.o RmatTranspose.o RdiagProduct.o MatranUtil.o
rmatspecdebug:  RmatSpecDebug.o
	$(FC)   MatranUtil.o Rmat.o RmatSpec.o RmatProduct.o RmatSum.o \
                RmatNorm.o RmatSpecDebug.o Rdiag.o RdiagProduct.o \
                RmatTranspose.o $(LIB)

RmatEye.o:  Rmat.o
RmatEyeDebug.o: Rmat.o RmatEye.o RmatPrint.o
rmateyedebug: RmatEyeDebug.o
	$(FC)  MatranUtil.o Rmat.o RmatEye.o RmatPrint.o RmatEyeDebug.o

RmatRand.o	: Rmat.o MatranUtil.o
RmatRandDebug.o	: Rmat.o RmatRand.o RmatPrint.o MatranUtil.o
rmatranddebug	: RmatRandDebug.o
	$(FC) MatranUtil.o Rmat.o RmatPrint.o RmatRand.o RmatRandDebug.o


RmatNorm2.o:      Rmat.o RmatProduct.o RmatSpec.o RmatNorm.o MatranUtil.o
RmatNorm2Debug.o: Rmat.o Rdiag.o RmatQR.o RmatTranspose.o RmatProduct.o \
                  RdiagProduct.o RmatNorm2.o MatranUtil.o
rmatnorm2debug:   RmatNorm2Debug.o
	$(FC)     RmatNorm2.o Rmat.o Rdiag.o RmatQR.o RmatTranspose.o \
                  RmatProduct.o RdiagProduct.o MatranUtil.o RmatNorm2Debug.o\
                  RmatSpec.o $(LIB)

RmatInv.o      : Rmat.o RmatChol.o RmatLudpp.o MatranUtil.o
RmatInvDebug.o : RmatProduct.o RmatSum.o RmatNorm.o RmatInv.o RmatEye.o \
	         RmatChol.o RmatLudpp.o MatranUtil.o Rmat.o MatranUtil.o
rmatinvdebug   : RmatInvDebug.o
	$(FC)  MatranUtil.o Rmat.o RmatInv.o RmatProduct.o RmatInvDebug.o \
	       RmatSum.o RmatNorm.o RmatEye.o RmatChol.o RmatLudpp.o \
	       $(LIB)

Rdiag.o: MatranUtil.o Rmat.o
RdiagDebug.o: Rdiag.o Rmat.o RmatPrint.o
rdiagdebug: RdiagDebug.o
	$(FC)   MatranUtil.o RdiagDebug.o Rdiag.o Rmat.o RmatPrint.o 


RdiagSum: MatranUtil.o Rdiag.o Rmat.o
RdiagSumDebug.o: Rdiag.o RdiagSum.o Rmat.o RmatPrint.o MatranUtil.o
rdiagsumdebug:   RdiagSumDebug.o
	$(FC)   MatranUtil.o RdiagSumDebug.o Rdiag.o Rmat.o \
                RmatPrint.o RdiagSum.o 

RdiagProduct.o:      Rdiag.o Rmat.o MatranUtil.o
RdiagProductDebug.o: Rdiag.o RdiagProduct.o Rmat.o RmatPrint.o \
                     MatranUtil.o
rdiagproductdebug:   RdiagProductDebug.o
	$(FC)        MatranUtil.o RdiagProductDebug.o Rdiag.o Rmat.o \
                     RmatPrint.o RdiagProduct.o $(LIB)

RdiagSolve.o:      Rdiag.o Rmat.o MatranUtil.o
RdiagSolveDebug.o: Rdiag.o RdiagSolve.o Rmat.o RmatPrint.o MatranUtil.o
rdiagsolvedebug:   RdiagSolveDebug.o
	$(FC)      MatranUtil.o RdiagSolveDebug.o RdiagSolve.o \
                   Rdiag.o Rmat.o RmatPrint.o $(LIB)

RdiagDiag.o:      Rdiag.o Rmat.o MatranUtil.o
RdiagDiagDebug.o: Rdiag.o Rmat.o RdiagDiag.o MatranUtil.o
rdiagdiagdebug:   RdiagDiagDebug.o
	$(FC)     MatranUtil.o Rdiag.o Rmat.o RdiagDiag.o RdiagDiagDebug.o 

RdiagJoin.o:      Rdiag.o Rmat.o MatranUtil.o
RdiagJoinDebug.o: Rdiag.o RdiagJoin.o Rmat.o RmatPrint.o MatranUtil.o
rdiagjoindebug:   RdiagJoinDebug.o
	$(FC)     MatranUtil.o Rdiag.o Rmat.o RdiagJoin.o RdiagJoinDebug.o \
                  RmatPrint.o

RmatSVD.o:      Rmat.o RmatTranspose.o MatranUtil.o
RmatSVDDebug.o:  Rmat.o Rdiag.o RmatTranspose.o RmatSVD.o RdiagProduct.o \
                 RmatProduct.o RmatSum.o RmatNorm.o RmatRand.o \
                 RmatSubmatrix.o RmatEye.o MatranUtil.o
rmatsvddebug:    RmatSVDDebug.o 
	$(FC)   MatranUtil.o Rmat.o Rdiag.o RmatSVD.o RmatProduct.o RmatSum.o \
                RmatTranspose.o RmatNorm.o RdiagProduct.o RmatSVDDebug.o \
                RmatRand.o RmatSubmatrix.o RmatEye.o $(LIB)

MatranRealCore.o: MatranUtil.o Rmat.o Rdiag.o RmatTranspose.o \
                  RmatSum.o RmatProduct.o RmatSolve.o RmatJoin.o \
                  RmatBorder.o RmatSubmatrix.o RdiagSum.o RdiagProduct.o \
                  RdiagSolve.o RdiagDiag.o RmatEye.o RmatNorm.o \
                  RmatPivot.o RmatPrint.o RmatRand.o RmatLudpp.o \
                  RmatChol.o
matranrealcore: MatranRealCore.o
	

