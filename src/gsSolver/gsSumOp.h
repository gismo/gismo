/** @file gsSumOp.h

    @brief Provides the sum of \a gsLinearOperator s as a \a gsLinearOperator.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): S. Takacs
*/
#pragma once

#include <gsSolver/gsLinearOperator.h>

namespace gismo
{

/// @brief Class for representing the sum of objects of type \a gsLinearOperator as \a gsLinearOperator
///
/// @ingroup Solver
template<typename T>
class gsSumOp GISMO_FINAL : public gsLinearOperator<T>
{
    typedef typename gsLinearOperator<T>::Ptr BasePtr;
public:

    /// Shared pointer for gsSumOp
    typedef memory::shared_ptr<gsSumOp> Ptr;

    /// Unique pointer for gsSumOp
    typedef memory::unique_ptr<gsSumOp> uPtr;

    /// Empty constructor. To be filled with addOperator()
    gsSumOp() : m_ops(0) {}

    /// Constructor taking a vector of Linear Operators
    gsSumOp(std::vector<BasePtr> ops)
        : m_ops(give(ops))
    {
#ifndef NDEBUG
        const size_t sz = m_ops.size();
        for (size_t i=0; i<sz; ++i)
        {
            GISMO_ASSERT ( m_ops[0]->rows() == m_ops[i]->rows() && m_ops[0]->cols() == m_ops[i]->cols(),
                           "Dimensions of the operators do not fit." );
        }
#endif
    }

    /// Constructor taking two Linear Operators
    gsSumOp(BasePtr op0, BasePtr op1)
        : m_ops(2)
    {
        GISMO_ASSERT ( op0->rows() == op1->rows() && op0->cols() == op1->cols(), "Dimensions of the operators do not fit." );
        m_ops[0] = give(op0); m_ops[1] = give(op1);
    }

    /// Constructor taking three Linear Operators
    gsSumOp(BasePtr op0, BasePtr op1, BasePtr op2 )
        : m_ops(3)
    {
        GISMO_ASSERT ( op0->rows() == op1->rows() && op0->cols() == op1->cols()
                        && op0->rows() == op2->rows() && op0->cols() == op2->cols(), "Dimensions of the operators do not fit." );
        m_ops[0] = give(op0); m_ops[1] = give(op1); m_ops[2] = give(op2);
    }

    /// Make command returning a smart pointer
    static uPtr make()
    { return uPtr( new gsSumOp() ); }

    /// Make command returning a smart pointer
    static uPtr make(std::vector<BasePtr> ops)
    { return uPtr( new gsSumOp(give(ops)) ); }

    /// Make command returning a smart pointer
    static uPtr make(BasePtr op0, BasePtr op1)
    { return uPtr( new gsSumOp(give(op0),give(op1)) ); }

    /// Make command returning a smart pointer
    static uPtr make(BasePtr op0, BasePtr op1, BasePtr op2)
    { return uPtr( new gsSumOp(give(op0),give(op1),give(op2)) ); }

    /// Add another operator
    void addOperator(BasePtr op)
    {
        GISMO_ASSERT ( m_ops.empty() || ( op->rows() == m_ops[0]->rows() && op->cols() == m_ops[0]->cols() ),
                       "Dimensions of the operators do not fit." );
        m_ops.push_back(give(op));
    }

    void apply(const gsMatrix<T> & input, gsMatrix<T> & x) const
    {
        GISMO_ASSERT ( !m_ops.empty(), "gsSumOp::apply does not work for 0 operators." );

        // Here, we could make a permanently allocated vector
        gsMatrix<T> tmp;
        const size_t sz = m_ops.size();

        m_ops[0]->apply(input,x);
        for (size_t i=1; i<sz; ++i)
        {
            m_ops[i]->apply(input,tmp);
            x += tmp;
        }
    }

    index_t rows() const
    {
        GISMO_ASSERT( !m_ops.empty(), "gsSumOp::rows does not work for 0 operators." );
        return m_ops[0]->rows();
    }

    index_t cols() const
    {
        GISMO_ASSERT( !m_ops.empty(), "gsSumOp::cols does not work for 0 operators." );
        return m_ops[0]->cols();
    }

    /// Return a vector of shared pointers to all operators
    const std::vector<BasePtr>& getOps() const { return m_ops; }

private:
    std::vector<BasePtr> m_ops;

};

}
