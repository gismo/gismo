#include <gsSolver/gsKronecker.h>

namespace gismo
{

void applyKronecker(const std::vector< gsLinearOperator* > & ops, const gsMatrix<>& x, gsMatrix<>& result)
{
    if (ops.size() == 1)        // deal with single-operator case efficiently
    {
        ops[0]->apply(x, result);
        return;
    }

    int sz = 1;

    for (unsigned i = 0; i < ops.size(); ++i)
    {
        GISMO_ASSERT(ops[i]->cols() == ops[i]->rows(), "Kronecker product only implemented for square operators");

        sz *= ops[i]->cols();
    }

    GISMO_ASSERT (sz == x.rows(), "Wrong size for input matrix");
    const index_t n = x.cols();

    // Note: algorithm relies on col-major matrices    
    gsMatrix<real_t, Dynamic, Dynamic, ColMajor> q0, q1;
    gsMatrix<real_t> temp1, temp2;

    // size: sz x n
    q0 = x;

    for (int i = (int)ops.size() - 1; i >= 0; --i)
    {
        // Re-order right-hand sides
        const int sz_i = ops[i]->cols();
        const int r_i  = sz / sz_i;
        q0.resize(sz_i, n * r_i);

        // Transpose solution component-wise
        q1.resize(r_i, n * sz_i);

        if (n == 1)     // optimization: save one matrix copy if input has only one column
        {
            ops[i]->apply(q0, temp2);
            q1 = temp2.transpose();
        }
        else
        {
            for (index_t k = 0; k != n; ++k)
            {
                temp1 = q0.middleCols(k*r_i, r_i);      //size: sz_i * r_i (coeffs for k-th rhs vector)
                ops[i]->apply(temp1, temp2);
                q1.middleCols(k*sz_i, sz_i) = temp2.transpose();
            }
        }

        q1.swap( q0 ); // move solution as next right-hand side
    }

    q0.resize(sz, n);
    result.swap( q0 );
}



void gsKroneckerProduct::calcSize()
{
    m_size = 1;

    for (unsigned i = 0; i < m_ops.size(); ++i)
    {
        GISMO_ASSERT(m_ops[i]->cols() == m_ops[i]->rows(), "Kronecker product only implemented for square operators");

        m_size *= m_ops[i]->cols();
    }
}


}

