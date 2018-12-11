/** @file gsAdditiveOp.h

    @brief Allows to set up additive Schwarz type preconditioners

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): S. Takacs
*/

#pragma once

#include <gsSolver/gsPreconditioner.h>
#include <gsSolver/gsMatrixOp.h>
#include <gsCore/gsMultiPatch.h>

namespace gismo
{

/// @brief Generic preconditioner which applies an arbitrary linear operator to the residual.
///
/// \code{.cpp}
///    gsLinearOperator<>::Ptr pc = gsAdditiveOp<>::make( transfers, ops );
/// \endcode
///
/// is equivalent to
///
/// \code{.cpp}
///    gsSumOp<>::Ptr pc = gsSumOp<>::make();
///
///    for (index_t i=0; i<transfers.size(); ++i)
///       s->addOperator(
///          gsProductOp<>::make(
///             makeMatrixOp(transfers[i].transpose()),
///             ops[i],
///             makeMatrixOp(transfers[i])
///          )
///       );
/// \endcode
///
/// but much faster.
///
/// @ingroup Solvers

template<class T>
class gsAdditiveOp : public gsLinearOperator<T>
{
    typedef typename gsLinearOperator<T>::Ptr   OpPtr;
    typedef std::vector<OpPtr>   OpContainer;
    typedef gsSparseMatrix<T,RowMajor>   Transfer;
    typedef std::vector<Transfer>   TransferContainer;

public:

    /// Shared pointer
    typedef memory::shared_ptr<gsAdditiveOp> Ptr;

    /// Unique pointer
    typedef memory::unique_ptr<gsAdditiveOp> uPtr;

    /// Default Constructor
    gsAdditiveOp() : m_transfers(), m_ops() {}

    /// Constructor
    gsAdditiveOp(TransferContainer transfers, OpContainer ops)
    : m_transfers(give(transfers)), m_ops(give(ops))
    {
#ifndef NDEBUG
        GISMO_ASSERT( m_transfers.size() == m_ops.size(), "Sizes do not agree" );
        const size_t sz = m_transfers.size();
        for (size_t i=0; i<sz; ++i)
        {
            GISMO_ASSERT ( m_transfers[i].rows()==m_transfers[0].rows() )
                       && m_transfers[i].cols() == m_op[i]->rows()
                       && m_op[i]->cols() == m_op[i]->rows(),
                       "Dimensions of the operators do not fit." );
        }
#endif        
    }

    /// Make function
    static uPtr make(TransferContainer transfers, OpContainer ops)
    { return uPtr( new gsAdditiveOp( give(transfers), give(ops) ) ); }

    /// Make function
    static uPtr make(std::pair<TransferContainer,OpContainer> transfer_and_ops)
    { return uPtr( new gsAdditiveOp( give(transfer_and_ops.first), give(transfer_and_ops.second) ) ); }
    
    /// Add another operator
    void addOperator(Transfer transfer, OpPtr op)
    {
        m_transfers.push_back(give(transfer));
        m_ops.push_back(give(op));
        GISMO_ASSERT ( transfer.rows()==m_transfers[0].rows()
                       && transfer.cols() == op->rows()
                       && op->cols() == op->rows(),
                       "Dimensions of the operators do not fit." );
    }

    void apply(const gsMatrix<T>& input, gsMatrix<T>& x) const;

    index_t rows() const
    {
        GISMO_ASSERT( !m_transfers.empty(), "gsAdditiveOp::rows does not work for 0 operators." );
        return m_transfers[0].rows();
    }

    index_t cols() const
    {
        GISMO_ASSERT( !m_transfers.empty(), "gsAdditiveOp::cols does not work for 0 operators." );
        return m_transfers[0].cols();
    }

protected:
    TransferContainer m_transfers;
    OpContainer m_ops;

};

} // namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsAdditiveOp.hpp)
#endif
