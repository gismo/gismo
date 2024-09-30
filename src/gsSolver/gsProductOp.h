/** @file gsProductOp.h

    @brief A class representing the product of \a gsLinearOperator s.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): C. Hofer
    Created on: 2016-08-30
*/


#pragma once

#include <gsSolver/gsLinearOperator.h>

namespace gismo {
/**
@brief Class for representing the product of objects of type \a gsLinearOperator as \a gsLinearOperator

For given operators \f$ A_1, A_2, ..., A_N\f$, it implements \f$ A_N,
..., A_1\f$, so \f$ A_1 \f$ is applied first, and \f$ A_N \f$ applied
last. The operator has dimensions A_N.rows() x A_1.cols().

For composition of preconditioners, cf. also \a gsCompositePrecOp

@ingroup Solver
    */
template<typename T>
class gsProductOp GISMO_FINAL : public gsLinearOperator<T>
{
    typedef typename gsLinearOperator<T>::Ptr BasePtr;
public:

    /// Shared pointer for gsProductOp
    typedef memory::shared_ptr<gsProductOp> Ptr;

    /// Unique pointer for gsProductOp
    typedef memory::unique_ptr<gsProductOp> uPtr;

    /// Empty constructor. To be filled with addOperator()
    gsProductOp() : m_ops() {}

    /// Constructor taking a vector of Linear Operators
    gsProductOp(std::vector<BasePtr> ops) : m_ops(give(ops))
    {
#ifndef NDEBUG
        const index_t sz = m_ops.size();
        for (index_t i=0; i<sz-1; ++i)
        {
            GISMO_ASSERT ( m_ops[i]->rows() == m_ops[i+1]->cols(),
                           "Dimensions of the operators do not fit." );
        }
#endif
    }

    /// Constructor taking two Linear Operators
    gsProductOp(BasePtr op0, BasePtr op1) : m_ops(2)
    {
        GISMO_ASSERT ( op0->rows() == op1->cols(),
                       "Dimensions of the operators do not fit." );
        m_ops[0]=give(op0); m_ops[1]=give(op1);
    }

    /// Constructor taking three Linear Operators
    gsProductOp(BasePtr op0, BasePtr op1, BasePtr op2) : m_ops(3)
    {
        GISMO_ASSERT ( op0->rows() == op1->cols() && op1->rows() == op2->cols(),
                       "Dimensions of the operators do not fit." );
        m_ops[0]=give(op0); m_ops[1]=give(op1); m_ops[2]=give(op2);
    }

    /// Make command returning a smart pointer
    static uPtr make()
    { return uPtr( new gsProductOp() ); }

    /// Make command returning a smart pointer
    static uPtr make(std::vector<BasePtr> ops)
    { return uPtr( new gsProductOp(give(ops)) ); }

    /// Make command returning a smart pointer
    static uPtr make(BasePtr op0, BasePtr op1)
    { return uPtr( new gsProductOp(give(op0),give(op1)) ); }

    /// Make command returning a smart pointer
    static uPtr make(BasePtr op0, BasePtr op1, BasePtr op2)
    { return uPtr( new gsProductOp(give(op0),give(op1),give(op2)) ); }

    /// Add another operator at the end
    void addOperator( BasePtr op )
    {
        GISMO_ASSERT ( m_ops.empty() || m_ops.back()->rows() == op->cols(),
                       "Dimensions of the operators do not fit." );
        m_ops.push_back(give(op));
    }

    void apply(const gsMatrix<T> & input, gsMatrix<T> & x) const
    {
        // The product of 0 operators is the identity
        if ( m_ops.empty() ) { x = input; return; }

        // Here, we could make a permanently allocated vector
        gsMatrix<T> temp;
        const size_t sz = m_ops.size();

        m_ops[0]->apply(input,x);
        for(size_t i=1; i<sz;++i)
        {
            temp.swap(x);
            m_ops[i]->apply(temp,x);
        }
    }

    index_t rows() const {
        GISMO_ASSERT( !m_ops.empty(), "gsProductOp::rows does not work for 0 operators.");
        return m_ops.back()->rows();

    }

    index_t cols() const {
        GISMO_ASSERT( !m_ops.empty(), "gsProductOp::cols does not work for 0 operators.");
        return m_ops.front()->cols();
    }

    /// Return a vector of shared pointers to all operators
    const std::vector<BasePtr>& getOps() const { return m_ops; }

private:
    std::vector<BasePtr> m_ops;
};

} // namespace gismo
