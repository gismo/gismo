/** @file linearAlgebra_example.cpp

    @brief Tutorial on matrix operations and linear algebra

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris
*/

#include <iostream>

#include <gismo.h>

using namespace gismo;

int main(int argc, char**argv)
{
    gsCmdLine cmd("Tutorial on matrix operations and linear algebra.");
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

#ifdef EIGEN_VECTORIZE
    gsDebug << "Vectorization is enabled in Eigen.\n";
#endif
#if EIGEN_HAS_RVALUE_REFERENCES
    gsDebug << "Eigen has rvalue references.\n";
#endif

    gsInfo << "G+Smo uses Eigen v"<< EIGEN_WORLD_VERSION<<"."
           <<EIGEN_MAJOR_VERSION<<"."<<EIGEN_MINOR_VERSION
           <<" for matrix computations.\n";

    // A matrix with entries of type real_t, and allocated size 3x3
    gsMatrix<real_t> A (3,3);
    // The comman initializer lets us fill the matrix. Note that the
    // matrix must have the correct size for this to work
    A << 2,2,3,  4,5,6,  7,8,10;
    A(0,0) -= 1 ;

    // If the type of the entries of the matrix is not given, the
    // default type is real_t (e.g. double)
    gsMatrix<> E (3,1);
    gsVector<> c (3);
    E << 2,2,3 ;

    // Even if two matrices do not have the same size we can assign one
    // to the other and the result will be two identical matrices
    c = E ;
    // Similary we can assign a column or any other expression to a
    // matrix variable
    c = A.col(1) ;

    // The two extra arguments are the number of rows and columns,
    // statically known at compile time
    gsMatrix<real_t,2,2> F ;
    F << 2,2,3, 4 ;
    gsVector<> w (2);
    w= F.row(1);

    gsInfo << "vector c:\n" << c <<"\n"<< E << "\n";

    gsInfo << "vector as diagonal:\n" << gsMatrix<>( c.asDiagonal() ) << "\n";

    gsInfo << "E.sum():\n" << E.sum() << "\n";

    gsInfo << "dyn: " << Dynamic << "\n";

    gsVector<> b (3);
    b << 3, 3, 4;

    gsInfo << "Here is the matrix A:\n" << A << "\n";
    gsInfo << "Here is the vector b:\n" << b << "\n";

    gsInfo << "Size of A: " << A.rows() << " x " << A.cols()  << "\n";
    gsInfo << "Determinant of A: "<< A.determinant()  << "\n";
    gsInfo << "Transpose of A:\n"<< A.transpose()  << "\n";

    // Note that A.transpose() does not alter A, but only returns an
    // expression of its transpose.
    // We can transpose A in place as follows:
    A.transposeInPlace();
    A.transposeInPlace();//second transposition results in the original A matrix

    // Note that A.inverse() does not alter A. To save the inverse in place of A use
    // A.inverseInPlace()
    gsInfo << "Inverse of A:\n"<< A.inverse()  << "\n";

    // We can initialize a matrix using other matrices by "<<"
    // Note that AAA must have the right size!
    gsMatrix<> AAA(3,12);
    AAA << A, A.transpose(), A.adjugate(), A;
    gsInfo << "A block matrix containing [A, A.tranpose(), A.adjugate(), A] :\n"<< AAA  << "\n";
    AAA.blockTransposeInPlace(3);
    gsInfo << "Block-wise transposition of the above:\n"<< AAA  << "\n";
    AAA.blockTransposeInPlace(6);
    gsInfo << "Block-wise transposition of the above seen as 3x6 blocks:\n"<< AAA  << "\n";

    gsVector<index_t> perm(3);
    perm << 2,1,0;
    gsInfo << "Here is the row permutation ("<<perm.transpose()<<") of matrix A:\n"
         << perm.asPermutation() * A << "\n";
    gsInfo << "Here is the column permutation ("<<perm.transpose()<<") of matrix A:\n"
         <<  A * perm.asPermutation() << "\n";

    gsInfo << "Here is the matrix A:\n" << A << "\n";

    gsVector<> x;
    // Computes a factorization (LU,QR) and solves Ax=b for the unknown x using
    // this factorization
    x= A.partialPivLu().solve(b);
    //x= A.fullPivLu().solve(b);
    //x= A.colPivHouseholderQr().solve(b);
    gsInfo << "The solution of Ax=b is:\n" << x << "\n";
    gsInfo << "Verification, A*x is:\n" << A*x << "\n";

    gsInfo << "The dot product x.b is : " <<  x.transpose()* b<< "\n"; //x.dot(b)
    gsInfo << "The dot product x.b is : " <<  x.dot( b )<< "\n"; //x.dot(b)

    gsInfo << "The product x*bt is : \n" << x *  b.transpose() << "\n";

    gsMatrix<> M  = x *  b.transpose() ;

    gsMatrix<real_t,3,3> W;

    W << 2,2,3,  4,5,6,  7,8,10;
    gsMatrix<> R2 = W * W ; //x * b.transpose() ;

    gsInfo << "Block of A of size (2,2), starting at (1,1):\n"<< A.block<2,2>(1,1) << "\n";
    // if the block size is not known at compile time:  A.block(1,1,2,2) << "\n";

    gsInfo << "Reverse matrix:\n"<< A.colwise().reverse() << "\n";

    gsSparseMatrix<> B(3,3);
    B.insert(0,0) = 1 ;
    B.insert(1,1) = 2 ;
    B.insert(2,2) = 3 ;
    B(1,1) += 3 ;

    gsInfo << "Here is a sparse matrix B:\n" << B<< " and B(1,1) is "<< B.coeffRef(1,1) << "\n";
    gsInfo << "Matrix B has "<<B.nonZeros()  << " non-zero elements"<< "\n";
    gsInfo << "Here is the product A*B:\n" << A*B << "\n";

    gsVector3d<> v1;
    v1 << 1,2,3;
    gsVector3d<> v2;
    v2 << 3,2,1;

    gsInfo << " dot product: "<< v1.dot(v2) << "\n";
    gsInfo << " cross product: "<< v1.cross(v2) << "\n";

    gsInfo << " dot product of matrix columns: "<< A.col(0).adjoint() * A.col(1) << "\n";
    gsInfo << " Another way: converts 1x1 matrix to value: "<< (A.col(0).transpose() * A.col(1) ).value() << "\n";

    gsMatrix<> r;
    A.firstMinor(0, 0, r);
    gsInfo << "Here are some minors of A:\n" << r  << "\n";
    A.firstMinor(1, 2, r);
    gsInfo << r  << "\n";
    A.firstMinor(2, 0, r);
    gsInfo << r  << "\n";
    A.firstMinor(2, 2, r);
    gsInfo << r  << "\n";

    r.setZero(2,2);
    gsInfo <<"Set matrix to zero setZero():\n"<< r <<"\n";
    r.setOnes();
    gsInfo <<"Set matrix to all ones setOnes():\n"<< r <<"\n";
    r.setConstant(3);
    gsInfo <<"Set matrix to all a constant setConstant(3):\n"<< r <<"\n";
    r.setRandom(2,2); // SLE_11_SP4
    gsInfo <<"Set matrix to random entires setRandom():\n"<< r <<"\n";

#ifndef GISMO_WITH_GMP // eigenvalues will not work for rational arithmetic types

    gsInfo << " Eigenvalues of non-symmetric matrix: "<< A.eigenvalues().transpose() << "\n";
    gsInfo << " Eigenvectors of non-symmetric matrix: \n"
         << gsMatrix<>::EigenSolver(A).eigenvectors() << "\n";

    gsInfo << " Eigenvalues of symmetric matrix (A's lower triangular part): "
         << A.selfadjointView<Lower>().eigenvalues().transpose()  << "\n";

    gsInfo << " Eigenvalues of symmetric matrix (A's upper triangular part): "
         << A.selfadjointView<Upper>().eigenvalues().transpose()  << "\n";

#endif

    return 0;

}


