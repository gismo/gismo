/** @file gsKroneckerOp.h

    @brief Provides a linear operator representing the Kronecker product of linear operators

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): C. Hofreither, S. Takacs
*/
#pragma once

#include <gsCore/gsExport.h>
#include <gsSolver/gsLinearOperator.h>

namespace gismo
{

/// @brief Class for representing a Kronecker product of operators of type \a gsLinearOperator
///
/// \ingroup Solver
template <class T>
class gsKroneckerOp : public gsLinearOperator<T>
{
    typedef typename gsLinearOperator<T>::Ptr BasePtr;
public:

    /// Shared pointer for gsKroneckerOp
    typedef memory::shared_ptr<gsKroneckerOp> Ptr;

    /// Unique pointer for gsKroneckerOp
    typedef memory::unique_ptr<gsKroneckerOp> uPtr;

    /// Kronecker product of a given list of operators
    gsKroneckerOp(const std::vector< BasePtr >& ops)
        : m_ops(ops)
    {
        GISMO_ASSERT( !m_ops.empty(), "Zero-term Kronecker product" );
        calcSize();
    }

    /// Convenience constructor for Kronecker product of two linear operators
    gsKroneckerOp(const BasePtr & op1, const BasePtr & op2)
        : m_ops(2)
    {
        m_ops[0] = op1;
        m_ops[1] = op2;
        calcSize();
    }

    /// Convenience constructor for Kronecker product of three linear operators
    gsKroneckerOp(const BasePtr & op1, const BasePtr & op2, const BasePtr & op3)
        : m_ops(3)
    {
        m_ops[0] = op1;
        m_ops[1] = op2;
        m_ops[2] = op3;
        calcSize();
    }

    /// Make command returning a smart pointer
    static uPtr make(const std::vector< BasePtr >& ops)
    { return memory::make_unique( new gsKroneckerOp(ops) ); }

    /// Make command returning a smart pointer
    /// Convenience function for Kronecker product of two linear operators
    static uPtr make(const BasePtr & op1, const BasePtr & op2)
    { return memory::make_unique( new gsKroneckerOp(op1,op2) ); }

    /// Make command returning a smart pointer
    /// Convenience function for Kronecker product of three linear operators
    static uPtr make(const BasePtr & op1, const BasePtr & op2, const BasePtr & op3)
    { return memory::make_unique( new gsKroneckerOp(op1,op2,op3) ); }

    virtual void apply(const gsMatrix<T> & input, gsMatrix<T> & result) const;
    virtual index_t rows() const    { return m_rows; }
    virtual index_t cols() const    { return m_cols; }

    static void applyKronecker(const std::vector<BasePtr> &, const gsMatrix<T> & input, gsMatrix<T> & result);

private:
    void calcSize();

    std::vector<BasePtr> m_ops;
    index_t m_rows, m_cols;
};

}

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsKronecker.hpp)
#endif
