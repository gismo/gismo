/** @file gsSumOfOperatorsOp.h

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
class gsSumOfOperatorsOp GISMO_FINAL : public gsLinearOperator<T>
{
    typedef typename gsLinearOperator<T>::Ptr BasePtr;
public:

    /// Shared pointer for gsSumOfOperatorsOp
    typedef memory::shared_ptr<gsSumOfOperatorsOp> Ptr;

    /// Unique pointer for gsSumOfOperatorsOp
    typedef memory::unique_ptr<gsSumOfOperatorsOp> uPtr;

    /// Empty constructor. To be filled with addOperator()
    gsSumOfOperatorsOp() : m_sz(0), m_ops(0) {}

    /// Constructor taking a vector of Linear Operators
    gsSumOfOperatorsOp(const std::vector<BasePtr>& ops)
        : m_sz(ops.size()), m_ops(ops)
    {
#ifndef NDEBUG
        for (index_t i=0; i<m_sz; ++i)
        {
            GISMO_ASSERT ( m_ops[0]->rows() == m_ops[i]->rows() && m_ops[0]->cols() == m_ops[i]->cols(),
                           "Dimensions of the operators do not fit." );
        }
#endif
    }

    /// Constructor taking two Linear Operators
    gsSumOfOperatorsOp(const BasePtr& op0, const BasePtr& op1)
        : m_sz(2), m_ops(2)
    {
        GISMO_ASSERT ( op0->rows() == op1->rows() && op0->cols() == op1->cols(), "Dimensions of the operators do not fit." );
        m_ops[0] = op0; m_ops[1] = op1;
    }

    /// Constructor taking three Linear Operators
    gsSumOfOperatorsOp(const BasePtr& op0, const BasePtr& op1, const BasePtr& op2 )
        : m_sz(3), m_ops(3)
    {
        GISMO_ASSERT ( op0->rows() == op1->rows() && op0->cols() == op1->cols()
                        && op0->rows() == op2->rows() && op0->cols() == op2->cols(), "Dimensions of the operators do not fit." );
        m_ops[0] = op0; m_ops[1] = op1; m_ops[2] = op2;
    }

    /// Make command returning a smart pointer
    static uPtr make()
    { return memory::make_unique( new gsSumOfOperatorsOp() ); }

    /// Make command returning a smart pointer
    static uPtr make( const std::vector<BasePtr>& ops )
    { return memory::make_unique( new gsSumOfOperatorsOp( ops ) ); }

    /// Make command returning a smart pointer
    static uPtr make( const BasePtr& op0, const BasePtr& op1 )
    { return memory::make_unique( new gsSumOfOperatorsOp( op0, op1 ) ); }

    /// Make command returning a smart pointer
    static uPtr make( const BasePtr& op0, const BasePtr& op1, const BasePtr& op2 )
    { return memory::make_unique( new gsSumOfOperatorsOp( op0, op1, op2 ) ); }

    /// Add another operator
    void addOperator( const BasePtr& op )
    {
        GISMO_ASSERT ( m_sz == 0 || ( op->rows() == m_ops[0]->rows() && op->cols() == m_ops[0]->cols() ),
                       "Dimensions of the operators do not fit." );
        ++m_sz;
        m_ops.push_back( op );
    }

    void apply(const gsMatrix<T> & input, gsMatrix<T> & result) const
    {
        GISMO_ASSERT ( m_sz>0, "gsSumOfOperatorsOp::apply does not work for 0 operators." );

        // Here, we could make a permanently allocated vector
        gsMatrix<T> tmp;

        m_ops[0]->apply(input,result);
        for (index_t i=1; i<m_sz; ++i)
        {
            m_ops[i]->apply(input,tmp);
            result += tmp;
        }
    }

    index_t rows() const
    {
        GISMO_ASSERT( m_sz>0, "gsSumOfOperatorsOp::rows does not work for 0 operators." );
        return m_ops[0]->rows();
    }

    index_t cols() const
    {
        GISMO_ASSERT( m_sz>0, "gsSumOfOperatorsOp::cols does not work for 0 operators." );
        return m_ops[0]->cols();
    }

    /// Return a vector of shared pointers to all operators
    const std::vector<BasePtr>& getOps() const { return m_ops; }

private:
    index_t m_sz;
    std::vector<BasePtr> m_ops;

};

}
