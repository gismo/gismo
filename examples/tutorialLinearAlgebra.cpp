/** @file tutorialLinearAlgebra.cpp

    @brief Tutorial on matrix operations and linear algebra

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris
*/

#include <iostream>

#include <gismo.h>


using namespace std;
using namespace gismo;

int main()
{
#ifdef EIGEN_VECTORIZE
    cout << "Vectorization is enabled in Eigen."<< endl;
#endif

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

    cout << "vector c:\n" << c <<"\n"<< E << endl;
  
    cout << "vector as diagonal:\n" << gsMatrix<>( c.asDiagonal() ) << endl;

    cout << "E.sum():\n" << E.sum() << endl;

    cout << "dyn: " << Dynamic << endl;

    gsVector<> b (3);
    b << 3, 3, 4;

    cout << "Here is the matrix A:\n" << A << endl;
    cout << "Here is the vector b:\n" << b << endl;
 
    cout << "Size of A: " << A.rows() << " x " << A.cols()  << endl;
    cout << "Determinant of A: "<< A.determinant()  << endl;
    cout << "Transpose of A:\n"<< A.transpose()  << endl;

    // Note that A.transpose() does not alter A, but only returns an
    // expression of its transpose.
    // We can transpose A in place as follows:
    A.transposeInPlace();
    A.transposeInPlace();//second transposition results in the original A matrix

    // Note that A.inverse() does not alter A. To save the inverse in place of A use
    // A.inverseInPlace()
    cout << "Inverse of A:\n"<< A.inverse()  << endl;

    // We can initialize a matrix using other matrices by "<<"
    // Note that AAA must have the right size!
    gsMatrix<> AAA(3,12);
    AAA << A, A.transpose(), A.adjugate(), A;
    cout << "A block matrix containing [A, A.tranpose(), A.adjugate(), A] :\n"<< AAA  << endl;
    AAA.blockTransposeInPlace(3);
    cout << "Block-wise transposition of the above:\n"<< AAA  << endl;
    AAA.blockTransposeInPlace(6);
    cout << "Block-wise transposition of the above seen as 3x6 blocks:\n"<< AAA  << endl;  

    gsVector<index_t> perm(3);
    perm << 2,1,0;
    cout << "Here is the row permutation ("<<perm.transpose()<<") of matrix A:\n" 
         << perm.asPermutation() * A << endl;
    cout << "Here is the column permutation ("<<perm.transpose()<<") of matrix A:\n" 
         <<  A * perm.asPermutation() << endl;
  
    gsVector<> x;
    // Computes QR factorization and solved Ax=b for the unknown x using
    // this factorization
    x= A.colPivHouseholderQr().solve(b);
    cout << "The solution of Ax=b is:\n" << x << endl;
    cout << "Verification, A*x is:\n" << A*x << endl;

    cout << "The dot product x.b is : " <<  x.transpose()* b<< endl; //x.dot(b)
    cout << "The dot product x.b is : " <<  x.dot( b )<< endl; //x.dot(b)

    cout << "The product x*bt is : \n" << x *  b.transpose() << endl;

    gsMatrix<> M  = x *  b.transpose() ;

    gsMatrix<real_t,3,3> W;

    W << 2,2,3,  4,5,6,  7,8,10;
    gsMatrix<> R2 = W * W ; //x * b.transpose() ;

    cout << "Block of A of size (2,2), starting at (1,1):\n"<< A.block<2,2>(1,1) << endl;
    // if the block size is not known at compile time:  A.block(1,1,2,2) << endl;

    cout << "Reverse matrix:\n"<< A.colwise().reverse() << endl;

    gsSparseMatrix<> B(3,3);
    B.insert(0,0) = 1 ;
    B.insert(1,1) = 2 ;
    B.insert(2,2) = 3 ;

    B(1,1) += 3 ;

    cout << "Here is a sparse matrix B:\n" << B<< " and B(1,1) is "<< B.coeffRef(1,1) << endl;
    cout << "Matrix B has "<<B.nonZeros()  << " non-zero elements"<< endl;
    cout << "Here is the product A*B:\n" << A*B << endl;

    gsVector3d<> v1;
    v1 << 1,2,3;
    gsVector3d<> v2;
    v2 << 3,2,1;

    cout << " dot product: "<< v1.dot(v2) << endl;
    cout << " cross product: "<< v1.cross(v2) << endl;
 
    cout << " dot product of matrix columns: "<< A.col(0).adjoint() * A.col(1) << endl;
    cout << " Another way: converts 1x1 matrix to value: "<< (A.col(0).transpose() * A.col(1) ).value() << endl;

    gsMatrix<> r;
    A.firstMinor(0, 0, r);
    cout << "Here are some minors of A:\n" << r  << endl;
    A.firstMinor(1, 2, r);
    cout << r  << endl;
    A.firstMinor(2, 0, r);
    cout << r  << endl;
    A.firstMinor(2, 2, r);
    cout << r  << endl;

    r.setRandom(2,2);
    cout <<"Set matrix to zero setZero():\n"<< r <<"\n";
    r.setOnes();
    cout <<"Set matrix to all ones setOnes():\n"<< r <<"\n";
    r.setConstant(3);
    cout <<"Set matrix to all a constant setConstant(3):\n"<< r <<"\n";
    r.setRandom();
    cout <<"Set matrix to random entires setRandom():\n"<< r <<"\n";

    cout << " Eigenvalues of non-symmetric matrix: "<< A.eigenvalues().transpose() << endl;
    cout << " Eigenvectors of non-symmetric matrix: \n"
         << Eigen::EigenSolver<gsMatrix<>::Base>(A).eigenvectors() << endl;

    cout << " Eigenvalues of symmetric matrix (A's lower triangular part): "
         << A.selfadjointView<Lower>().eigenvalues().transpose()  << endl;

    cout << " Eigenvalues of symmetric matrix (A's upper triangular part): "
         << A.selfadjointView<Upper>().eigenvalues().transpose()  << endl;

    return 0;

}


