/** @file gsCompositePrecOp.h

    @brief Allows to represent the composition of preconditioners.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): S. Takacs
*/

#pragma once

#include<gsCore/gsLinearAlgebra.h>
#include<gsSolver/gsLinearOperator.h>

namespace gismo
{

/// @brief This class represents the composition of preconditioners of type \a gsPreconditionerOp.
///
/// If the preconditioners have the iteration matrices \f$ I-P_iA \f$, the overall preconditioner
/// has the iteration matrix \f$ (I-P_nA)\cdots (I-P_1A) \f$
///
/// This should not be confused with \a gsProductOp, which would yield
///
/// \f$ I - P_n \cdots P_1 A \f$.
///
/// @ingroup Solver
template <class T>
class gsCompositePrecOp GISMO_FINAL : public gsPreconditionerOp<T>
{
    typedef typename gsPreconditionerOp<T>::Ptr BasePtr;
public:

    /// Shared pointer for gsCompositePrecOp
    typedef memory::shared_ptr<gsCompositePrecOp> Ptr;

    /// Unique pointer for gsCompositePrecOp
    typedef memory::unique_ptr<gsCompositePrecOp> uPtr;

    /// Empty constructor. To be filled with addOperator()
    gsCompositePrecOp() : m_ops() {}

    /// Constructor taking a vector of preconditioners
    gsCompositePrecOp(std::vector< BasePtr > ops)
        : m_ops(give(ops)) {}

    /// Convenience constructor taking two preconditioners
    gsCompositePrecOp(BasePtr op0, BasePtr op1)
        : m_ops(2)
    {
        m_ops[0] = give(op0); m_ops[1] = give(op1);
    }

    /// Convenience constructor taking three preconditioners
    gsCompositePrecOp(BasePtr op0, BasePtr op1, BasePtr op2)
        : m_ops(3)
    {
        m_ops[0] = give(op0); m_ops[1] = give(op1); m_ops[2] = give(op2);
    }

    /// Make command returning a smart pointer
    static uPtr make(std::vector< BasePtr > ops)
    {
        return uPtr( new gsCompositePrecOp(give(ops)) );
    }

    /// Make command returning a smart pointer
    static uPtr make(BasePtr op0, BasePtr op1)
    {
        return uPtr( new gsCompositePrecOp(give(op0),give(op1)) );
    }

    /// Make command returning a smart pointer
    static uPtr make(BasePtr op0, BasePtr op1, BasePtr op2)
    {
        return uPtr( new gsCompositePrecOp(give(op0),give(op1),give(op2)) );
    }

    /// Add another operator at the end
    void addOperator(BasePtr op)
    {
        m_ops.push_back(give(op));
    }

    /// Apply the smoother for the equation Ax=rhs and update the current iterate x.
    virtual void step(const gsMatrix<T>& rhs, gsMatrix<T>& x) const
    {
        const size_t sz = m_ops.size();
        for ( size_t i=0; i<sz; ++i )
            m_ops[i]->step(rhs,x);
    }

    /// Apply the transposed smoother for the equation Ax=rhs and update the current iterate x.
    virtual void stepT(const gsMatrix<T>& rhs, gsMatrix<T>& x) const
    {
        const index_t sz = m_ops.size();
        for ( index_t i=sz-1; i>=0; --i )
            m_ops[i]->stepT(rhs,x);
    }

    typename gsLinearOperator<T>::Ptr underlyingOp() const
    { GISMO_NO_IMPLEMENTATION }

    index_t rows() const
    {
        GISMO_ASSERT( !m_ops.empty(), "gsCompositePrecOp::rows does not work for 0 operators.");
        return m_ops[0]->rows();
    }

    index_t cols() const
    {
        GISMO_ASSERT( !m_ops.empty(), "gsCompositePrecOp::rows does not work for 0 operators.");
        return m_ops[0]->cols();
    }

private:
    std::vector< BasePtr > m_ops;
};


} // namespace gismo
