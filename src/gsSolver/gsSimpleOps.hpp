/** @file gsSimpleOps.hpp

    @brief Collection of some simple preconditioners.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): C. Hofreither
*/
#include <gsSolver/gsSimpleOps.h>
#include <gsCore/gsLinearAlgebra.h>

namespace gismo
{

template<typename T>
void dampedRichardsonSweep(const gsSparseMatrix<T>& A, gsMatrix<T>& x, const gsMatrix<T>& f, T tau)
{
    GISMO_ASSERT( A.rows() == x.rows() && x.rows() == f.rows() && A.cols() == A.rows() && x.cols() == f.cols(),
        "Dimensions do not match.");

    gsMatrix<T> temp = f - A * x;
    x += tau * temp;
}

template<typename T>
void jacobiSweep(const gsSparseMatrix<T>& A, gsMatrix<T>& x, const gsMatrix<T>& f)
{
    GISMO_ASSERT( A.rows() == x.rows() && x.rows() == f.rows() && A.cols() == A.rows() && x.cols() == f.cols(),
        "Dimensions do not match.");
    
    GISMO_ASSERT( f.cols() == 1, "This operator is only implemented for a single right-hand side." );

    gsMatrix<T> temp = f - A * x;
    temp.array() /= A.diagonal().array();
    x += temp;
}

template<typename T>
void dampedJacobiSweep(const gsSparseMatrix<T>& A, gsMatrix<T>& x, const gsMatrix<T>& f, T tau)
{
    GISMO_ASSERT( A.rows() == x.rows() && x.rows() == f.rows() && A.cols() == A.rows() && x.cols() == f.cols(),
        "Dimensions do not match.");
    
    GISMO_ASSERT( f.cols() == 1, "This operator is only implemented for a single right-hand side." );

    gsMatrix<T> temp = f - A * x;
    temp.array() /= A.diagonal().array();
    x += tau * temp;
}

template<typename T>
void gaussSeidelSweep(const gsSparseMatrix<T>& A, gsMatrix<T>& x, const gsMatrix<T>& f)
{
    GISMO_ASSERT( A.rows() == x.rows() && x.rows() == f.rows() && A.cols() == A.rows() && x.cols() == f.cols(),
        "Dimensions do not match.");
    
    GISMO_ASSERT( f.cols() == 1, "This operator is only implemented for a single right-hand side." );

    // A is supposed to be symmetric, so it doesn't matter if it's stored in row- or column-major order
    for (int i = 0; i < A.outerSize(); ++i)
    {
        T diag = 0;
        T sum  = 0;

        for (typename gsSparseMatrix<T>::InnerIterator it(A,i); it; ++it)
        {
            sum += it.value() * x( it.index() );        // compute A.x
            if (it.index() == i)
                diag = it.value();
        }

        x(i) += (f(i) - sum) / diag;
    }
}

template<typename T>
void reverseGaussSeidelSweep(const gsSparseMatrix<T>& A, gsMatrix<T>& x, const gsMatrix<T>& f)
{
    GISMO_ASSERT( A.rows() == x.rows() && x.rows() == f.rows() && A.cols() == A.rows() && x.cols() == f.cols(),
        "Dimensions do not match.");
    
    GISMO_ASSERT( f.cols() == 1, "This operator is only implemented for a single right-hand side." );

    // A is supposed to be symmetric, so it doesn't matter if it's stored in row- or column-major order
    for (int i = A.outerSize() - 1; i >= 0; --i)
    {
        T diag = 0;
        T sum = 0;

        for (typename gsSparseMatrix<T>::InnerIterator it(A,i); it; ++it)
        {
            sum += it.value() * x( it.index() );        // compute A.x
            if (it.index() == i)
                diag = it.value();
        }

        x(i) += (f(i) - sum) / diag;
    }
}

//Assumes A is symmetric (not needed)!
template<typename T>
void gaussSeidelSingleBlock(const gsSparseMatrix<T>& A, gsMatrix<T>& x, const gsMatrix<T>& f, gsVector<index_t>& DoFs)
{
    GISMO_ASSERT( A.rows() == x.rows() && x.rows() == f.rows() && A.cols() == A.rows() && x.cols() == f.cols(),
        "Dimensions do not match.");
    
    GISMO_ASSERT( f.cols() == 1, "This operator is only implemented for a single right-hand side." );    
    
    //Sorting from lowest to highest
    DoFs.sortByColumn(0);
    const index_t size = DoFs.rows();

    GISMO_ASSERT ( DoFs(size-1)< A.cols(), "The given DoF is higher than the size of the matrix");

    gsMatrix<T> Dblock(size, size);
    gsMatrix<T> residual(size, 1);
    for (int i = 0; i< size; i++)
    {
        //Symmetry is assumed here!
        residual(i,0) = f(DoFs(i),0) - (A.innerVector(DoFs(i)).transpose()*x).value();
        for (int j = 0; j< size; j++)
            Dblock(i,j) = A.coeff(DoFs(i), DoFs(j));
    }
    //Multiply residual with the inverse of the diagonal block
    residual = Dblock.inverse()*residual;

    //Update
    for (int i = 0; i< size; i++)
       x(DoFs(i),0) += residual(i,0);
}

} // namespace gismo

