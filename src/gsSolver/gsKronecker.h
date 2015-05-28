/** @file gsKronecker.h

    @brief Provides functions and classes for working with Kronecker products of matrices and operators.

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): C. Hofreither
*/
#pragma once

#include <vector>
#include <gsCore/gsExport.h>
#include <gsSolver/gsLinearOperator.h>

namespace gismo
{

/**
  Compute the application of the Kronecker product
    kron(ops[0], ops[1], ..., ops[n-1]) * x
  and store it in \a result without computing the large Kronecker product matrix itself.
*/
void GISMO_EXPORT applyKronecker(const std::vector< gsLinearOperator* > & ops, const gsMatrix<>& x, gsMatrix<>& result);



/// Class for representing a Kronecker product of linear operators
class GISMO_EXPORT gsKroneckerProduct : public gsLinearOperator
{
public:
    /// Kronecker product of a given list of operators. Takes ownership of the operators.
    gsKroneckerProduct(const std::vector< gsLinearOperator* >& ops)
        : m_ops(ops)
    {
        GISMO_ASSERT( !m_ops.empty(), "Zero-term Kronecker product" );
        calcSize();
    }

    /// Convenience constructor for Kronecker product of two linear operators
    gsKroneckerProduct(gsLinearOperator * op1, gsLinearOperator * op2)
    {
        m_ops.resize(2);
        m_ops[0] = op1;
        m_ops[1] = op2;
        calcSize();
    }

    /// Convenience constructor for Kronecker product of three linear operators
    gsKroneckerProduct(gsLinearOperator * op1, gsLinearOperator * op2, gsLinearOperator * op3)
    {
        m_ops.resize(3);
        m_ops[0] = op1;
        m_ops[1] = op2;
        m_ops[2] = op3;
        calcSize();
    }

    virtual ~gsKroneckerProduct()
    {
        freeAll(m_ops);
    }

    virtual void apply(const gsMatrix<real_t> & input, gsMatrix<real_t> & result) const
    {
        applyKronecker(m_ops, input, result);
    }

    virtual index_t rows() const    { return m_size; }
    virtual index_t cols() const    { return m_size; }

private:
    // disable copy constructor and assignment operator since for now we have no good
    // way of cloning the linear operators in m_ops
    gsKroneckerProduct(const gsKroneckerProduct& other);
    gsKroneckerProduct& operator=(const gsKroneckerProduct& other);

private:
    void calcSize();

    std::vector< gsLinearOperator* > m_ops;
    int m_size;
};


}

