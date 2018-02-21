/** @file gsProductOfOperatorsOp.h

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

/// @brief Class for representing the product of objects of type \a gsLinearOperator as \a gsLinearOperator
///
/// Note, for given operators \$f A_1, A_2, ..., A_N\$f, it implements
/// \$f A_N \cdots A_1$, so \$ A_1 \$f is applied first, and \$ A_N \$f applied last.
///
/// For composition of preconditioners, cf. also \a gsCompositionOfPreconditionersOp
///
/// @ingroup Solver
template<typename T>
class gsProductOfOperatorsOp GISMO_FINAL : public gsLinearOperator<T>
{
    typedef typename gsLinearOperator<T>::Ptr BasePtr;
public:

    /// Shared pointer for gsProductOfOperatorsOp
    typedef memory::shared_ptr<gsProductOfOperatorsOp> Ptr;

    /// Unique pointer for gsProductOfOperatorsOp
    typedef memory::unique_ptr<gsProductOfOperatorsOp> uPtr;

    /// Empty constructor. To be filled with addOperator()
    gsProductOfOperatorsOp() : m_ops() {}

    /// Constructor taking a vector of Linear Operators
    gsProductOfOperatorsOp(const std::vector<BasePtr>& ops) : m_ops(ops)
    {
#ifndef NDEBUG
        for (size_t i=0; i<m_ops.size()-1; ++i)
        {
            GISMO_ASSERT ( m_ops[i]->rows() == m_ops[i+1]->cols(),
                           "Dimensions of the operators do not fit." );
        }
#endif
    }

    /// Constructor taking two Linear Operators
    gsProductOfOperatorsOp(const BasePtr & op0, const BasePtr & op1 ) : m_ops(2)
    {
        GISMO_ASSERT ( op0->rows() == op1->cols(),
                       "Dimensions of the operators do not fit." );
        m_ops[0]=op0; m_ops[1]=op1;
    }

    /// Constructor taking three Linear Operators
    gsProductOfOperatorsOp(const BasePtr & op0, const BasePtr & op1, const BasePtr & op2) : m_ops(3)
    {
        GISMO_ASSERT ( op0->rows() == op1->cols() && op1->rows() == op2->cols(),
                       "Dimensions of the operators do not fit." );
        m_ops[0]=op0; m_ops[1]=op1; m_ops[2]=op2;
    }

    /// Make command returning a smart pointer
    static uPtr make()
    { return memory::make_unique( new gsProductOfOperatorsOp() ); }

    /// Make command returning a smart pointer
    static uPtr make(const std::vector<BasePtr>& ops)
    { return memory::make_unique( new gsProductOfOperatorsOp(ops) ); }

    /// Make command returning a smart pointer
    static uPtr make(const BasePtr & op0, const BasePtr & op1)
    { return memory::make_unique( new gsProductOfOperatorsOp(op0,op1) ); }

    /// Make command returning a smart pointer
    static uPtr make(const BasePtr & op0, const BasePtr & op1, const BasePtr & op2)
    { return memory::make_unique( new gsProductOfOperatorsOp(op0,op1,op2) ); }

    /// Add another operator at the end
    void addOperator( const BasePtr& op )
    {
        GISMO_ASSERT ( m_ops.size() == 0 || m_ops.back()->cols() == op->rows(),
                       "Dimensions of the operators do not fit." );
        m_ops.push_back( op );
    }

    void apply(const gsMatrix<T> & input, gsMatrix<T> & x) const
    {
        // The product of 0 operators is the identity
        if ( m_ops.size() == 0 ) { x = input; return; }

        // Here, we could make a permanently allocated vector
        gsMatrix<T> temp;

        m_ops[0]->apply(input,x);
        for(size_t i=1; i<m_ops.size();++i)
        {
            temp.swap(x);
            m_ops[i]->apply(temp,x);
        }
    }

    index_t rows() const {
        GISMO_ASSERT( m_ops.size()>0, "gsProductOfOperatorsOp::rows does not work for 0 operators.");
        return m_ops.back()->rows();

    }

    index_t cols() const {
        GISMO_ASSERT( m_ops.size()>0, "gsProductOfOperatorsOp::cols does not work for 0 operators.");
        return m_ops.front()->cols();
    }

    /// Return a vector of shared pointers to all operators
    const std::vector<BasePtr>& getOps() const { return m_ops; }

private:
    std::vector<BasePtr> m_ops;
};

} // namespace gismo
