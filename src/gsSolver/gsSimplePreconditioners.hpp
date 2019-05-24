/** @file gsSimplePreconditioners.hpp

    @brief Collection of some simple preconditioners.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): C. Hofreither
*/

namespace gismo
{

namespace internal
{

template<typename T>
void gaussSeidelSweep(const gsSparseMatrix<T> & A, gsMatrix<T>& x, const gsMatrix<T>& f)
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
void reverseGaussSeidelSweep(const gsSparseMatrix<T> & A, gsMatrix<T>& x, const gsMatrix<T>& f)
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

} // namespace internal

} // namespace gismo
