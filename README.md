# Matran
Matrix wrapper for Fortran 95 by G. W. Stewart, copied from https://www.cs.umd.edu/~stewart/matran/Matran.html, which has the following description:

## Introduction
Matran (pronounced MAY-tran) is a Fortran 95 wrapper that implements matrix operations and computes matrix decompositions using Lapack and the Blas. Although Matran is not based on a formally defined matrix language, it provides the flavor and convenience of coding in matrix oriented systems like Matlab, Octave, etc. By using routines from Lapack and the Blas, Matran allows the user to obtain the computational benefits of these packages with minimal fuss and bother.
Detailed information about Matran may be found in the Matran Writeup (ps, pdf). Here we give a general overview of Matran, its organization, and its capabilities.

## Organization
Matran consists of a number of Fortran 95 modules, which fall into five categories.

### The Matran utility module
This module contains global parameters, error handlers, and utility routines.

### Matrix types
A Matran matrix type is essentially a storage class for matrices. Currently implemented are the Rmat and the Rdiag. The Rmat represents matrices that can be conveniently stored in a rectangular array (whose dimensions may be larger than that of the matrix). The Rdiag represents diagonal matrices whose diagonal elements are stored in a linear array. Matran can be compiled to give a single or double precision package (but not both at the same time).
The Rmat has a special tag field that specifies if the matrix is general, triangular, Hermitian, or positive (semi) definite. Matran uses this field to the most efficient way to perform operations with the matrices. There are no packed representations. For example, the zero elements of a triangular matrix are explicitly stored in the Rmat array.

The Rdiag implements a diagonal matrix with its diagonal elements stored in a linear array.

Although they are not currently implemented, Matran will eventually have complex types, Cmat and Cdiag, corresponding to the real types.

### Matrix operations
Matran overloads old operators and defines new operators to compute transposes (.trp.A, .ctp.A), sums (A+B, -A), and products (A*B)of matrices. In addition it defines operators to compute quantities like A-1B (A.xiy.B). There are also operators for concatenating matrices and extracting submatrices.
These operations are organized into collections of modules called suites. For example, there are two modules in the product suite that compute matrix products. The module RmatProduct_m implements products between Rmats, taking due account of the tag components mentioned above. The module RdiagProduct implements products between Rdiags and between Rdiags and Rmats. When the type Cmat is introduced, the module CmatProduct_m will be responsible for implementing products between any combination of Rmats, Rdiags, and Cmats.

### Matrix miscelanea
Matran provides a number of useful matrix functions (e.g., norms) and constructors (e.g., identity matrices). It also provides a module to pretty print a matrix.

### Decompositions
Matwrap provides types and constructors for the following decompositions: LU with partial pivoting, QR, QR with column pivoting, Cholesky, spectral, SVD, real Schur, and the eigendecompositon of a general matrix.
A recurring problem in matrix oriented languages is how to prevent the recomputation of matrix decompositions in loops. Matran provides features that help the programmer circumvent this problem.

<b>Openness</b><br>
Modules in Fortran 95 cannot see private parts of other modules. For this reason Matran's modules have no private parts at all. Although this goes against the grain of object oriented programming, it actually makes a good deal of sense when it comes to matrix computations. The number of uses of matrices (which have been around for 150 years) is far greater than the number of methods that can be folded into an object oriented package. The present arrangements allow the user to extend the package efficiently by operating directly on the arrays containing the matrices. For example, all the powerful array operations of Fortran 95 are directly available to the Matran programmer. The dark side of this flexibility is that the programmer must exercise discipline in enforcing a number of Matran conventions.

<b>Memory management</b><br>
Matran handles all memory allocation automatically. It is conservative and attempts to use available storage wherever possible. Unfortunately, the storage of matrix types and decompositions declared in subprograms as well as temporaries passed to subprograms are not automatically deallocated when the program returns from the subprogram, and the coder must explicitly deallocate the storage. However, Matran provides special generic subroutines to make this task easy.
Matran allows considerable control over the creation of temporary matrix objects by shadowing matrix operators with subroutine forms. For example, the statement `C = A + B`
requires a temporary Rmat to hold the intermediate result A + B. On the other hand the statement `call Plus(C, A, B)` accomplishes the same thing with no temporaries.
Matran performs most of its computations by calling appropriate Lapack and Blas routines, some of which require the user to furnish working storage. Ordinarily, this storage is silently allocated and deallocated by Matran. But through optional arguments, the coder can furnish the storage explicitly, thus reducing calls to the allocator.

## Current Status
Matran is at a point where it can perform useful computations (see the last section of the Matran Writeup (ps, pdf for an example), which is why it is being distributed. However, it can certainly be improved, and a second purpose of the distribution is to get feedback before Matran gets so large it is hard to change. Please send comments and suggestions to me at stewart@cs.umd.edu. There are three areas in which Matran will likely change:
The addition of types and operations for complex matrices.
The expansion of the matrix miscelanea class of modules.
The replacement of pointer arrays with allocatable arrays. This will happen when the allocatable array extensions to Fortran 95 are widely and correctly implemented.
