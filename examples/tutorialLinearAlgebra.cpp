#include <iostream>

#include <gismo.h>


using namespace std;
using namespace gismo;

int main()
{
#ifdef EIGEN_VECTORIZE
  cout << "Vectorization is enabled in Eigen."<< endl;
#endif

  gsMatrix<> A (3,3);
  A << 2,2,3,  4,5,6,  7,8,10;
  A(0,0) -= 1 ;

  gsMatrix<> E (3,1);
  gsVector<> c (3);
  E << 2,2,3 ;

  c = E ;
  c = A.col(1) ;

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
  A.transposeInPlace();
  A.transposeInPlace();
  cout << "Inverse of A:\n"<< A.inverse()  << endl;
  
  gsVector<> x;
  x= A.colPivHouseholderQr().solve(b);

  cout << "The solution of Ax=b is:\n" << x << endl;

  cout << "Verification, A*x is:\n" << A*x << endl;


  cout << "The dot product x.b is : " <<  x.transpose()* b<< endl; //x.dot(b)
  cout << "The dot product x.b is : " <<  x.dot( b )<< endl; //x.dot(b)

  cout << "The product x*bt is : \n" << x *  b.transpose() << endl;

  gsMatrix<> M  = x *  b.transpose() ;

  //gsMatrix<real_t,Dynamic,Dynamic> W;
  gsMatrix<real_t,3,3> W;

  W << 2,2,3,  4,5,6,  7,8,10;
  gsMatrix<> R2= W * W ; //x *  b.transpose() ;

  //cout << "The cross product x times b is: " << x.cross(b) << endl;


  cout << "Block of A of size (2,2), starting at (1,1):\n"<< A.block<2,2>(1,1) << endl;
  // OR:  A.block(1,1,2,2) << endl;

  cout << "Reverse matrix:\n"<< A.colwise().reverse() << endl;

  gsSparseMatrix<> B(3,3);
  B.insert(0,0) = 1 ;
  B.insert(1,1) = 2 ;
  B.insert(2,2) = 3 ;

  B(1,1) += 3 ;

  cout << "Here is a sparse matrix B:\n" << B<< " and B(1,1) is "<< B.coeffRef(1,1) << endl;
  cout << "Matrix B has "<<B.nonZeros()  << " non-zero elements"<< endl;
  cout << "Here is the product A*B:\n" << A*B << endl;


  cout << "--- Check bug ---" << endl;
  gsSparseMatrix<> I(9,9);
  for (index_t i=0; i<9; ++i )
      I(i,i)=1;
  gsMatrix<> r(9,1);
  r.setZero();
  gsSparseSolver<>::CGDiagonal solver;
  gsMatrix<> rr(9,1);
  rr =  solver.compute(I).solve( r );
  cout << " sol: "<< rr.transpose()  << endl;


  gsVector3d<> v1;
  v1 << 1,2,3;
  gsVector3d<> v2;
  v2 << 3,2,1;

  cout << " dot product: "<< v1.dot(v2) << endl;
  cout << " cross product: "<< v1.cross(v2) << endl;
 
  cout << " dot product of matrix columns: "<< A.col(0).adjoint() * A.col(1) << endl;
  cout << " Another way: converts 1x1 matrix to value: "<< (A.col(0).transpose() * A.col(1) ).value() << endl;

  A.firstMinor(0, 0, r);
  cout << "Here are some minors of A:\n" << r  << endl;
  A.firstMinor(1, 2, r);
  cout << r  << endl;
  A.firstMinor(2, 0, r);
  cout << r  << endl;
  A.firstMinor(2, 2, r);
  cout << r  << endl;


  cout << " Eigenvalues of non-symmetric matrix: "<< A.eigenvalues().transpose() << endl;
  cout << " Eigenvectors of non-symmetric matrix: \n"
       << Eigen::EigenSolver<gsMatrix<>::Base>(A).eigenvectors() << endl;

  cout << " Eigenvalues of symmetric matrix (A's lower triangular part): "
       << A.selfadjointView<Lower>().eigenvalues().transpose()  << endl;

  cout << " Eigenvalues of symmetric matrix (A's upper triangular part): "
       << A.selfadjointView<Upper>().eigenvalues().transpose()  << endl;

  return 0;

}


// setZero(), setOnes(), setIdentity(), setConstant(value), setRandom()

// Note. 3D: Eigen::Matrix<ArrayXd, Dynamic, Dynamic>
