/** @file gsKroneckerOp.hpp

    @brief Provides a linear operator representing the Kronecker product of linear operators

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): C. Hofreither, S. Takacs
*/

namespace gismo
{

/// @cond
template <typename T>
void gsKroneckerOp<T>::apply(const std::vector<typename gsLinearOperator<T>::Ptr> & ops, const gsMatrix<T> & input, gsMatrix<T> & x)
{
    GISMO_ASSERT( !ops.empty(), "Zero-term Kronecker product" );
    const index_t nrOps = ops.size();

    if (nrOps == 1)        // deal with single-operator case efficiently
    {
        ops[0]->apply(input, x);
        return;
    }

    //index_t rows = 1; // we don't compute rows as we do not need it
    index_t sz = 1;

    for (index_t i = 0; i < nrOps; ++i)
    {
        //rows *= ops[i]->rows();
        sz *= ops[i]->cols();
    }

    GISMO_ASSERT (sz == input.rows(), "The input matrix has wrong size.");
    const index_t n = input.cols();

    // Note: algorithm relies on col-major matrices
    gsMatrix<T, Dynamic, Dynamic, ColMajor> q0, q1;
    gsMatrix<T> temp;

    // size: sz x n
    q0 = input;

    for (index_t i = nrOps - 1; i >= 0; --i)
    {

        const index_t cols_i = ops[i]->cols();
        const index_t rows_i = ops[i]->rows();
        const index_t r_i  = sz / cols_i;

        // Re-order right-hand sides
        q0.resize(cols_i, n * r_i);

        // Apply operator
        ops[i]->apply(q0, temp);
        GISMO_ASSERT (temp.rows() == rows_i && temp.cols() == n * r_i, "The linear operator returned a matrix with unexpected size.");

        // Transpose solution component-wise
        if (n == 1)
            q1 = temp.transpose();
        else
        {
            q1.resize(r_i, n * rows_i);
            for (index_t k = 0; k != n; ++k)
                q1.middleCols(k*rows_i, rows_i) = temp.middleCols(k*r_i, r_i).transpose();
        }

        q1.swap( q0 ); // move solution as next right-hand side
        sz = ( sz / cols_i) * rows_i; // update now the dimensionality such that in the end sz = rows
        //sz % cols_i == 0, since sz *= ops[i]->cols();
    }

    //GISMO_ASSERT (sz == rows, "Internal error."); // see one above

    q0.resize(sz, n);
    x.swap( q0 );
}
/// @endcond

template <typename T>
void gsKroneckerOp<T>::apply(const gsMatrix<T> & input, gsMatrix<T> & x) const
{
    apply(m_ops, input, x);
}

template <typename T>
index_t gsKroneckerOp<T>::rows() const
{
    index_t rows = 1;
    for (size_t i = 0; i < m_ops.size(); ++i)
        rows *= m_ops[i]->rows();
    return rows;
}

template <typename T>
index_t gsKroneckerOp<T>::cols() const
{
    index_t cols = 1;
    for (size_t i = 0; i < m_ops.size(); ++i)
        cols *= m_ops[i]->cols();
    return cols;
}


}
