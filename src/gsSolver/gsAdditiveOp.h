/** @file gsAdditiveOp.h

    @brief Allows to set up additive Schwarz type preconditioners

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

/// @brief Generic preconditioner which applies an arbitrary linear operator to the residual.
///
/// This preconditioner realizes \f$ \sum_{i=1}^n T_i A_i T_i^T \f$, where the
/// \f$ T_i \f$ are the transfer matrices and the \f$ A_i \f$ are linear operators
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
class gsAdditiveOp GISMO_FINAL : public gsLinearOperator<T>
{
    typedef typename gsLinearOperator<T>::Ptr   OpPtr;
    typedef std::vector<OpPtr>                  OpContainer;
    typedef gsSparseMatrix<T,RowMajor>          Transfer;
    typedef memory::shared_ptr<Transfer>        TransferPtr;
    typedef std::vector<Transfer>               TransferContainer;
    typedef std::vector<TransferPtr>            TransferPtrContainer;

public:

    /// Shared pointer
    typedef memory::shared_ptr<gsAdditiveOp> Ptr;

    /// Unique pointer
    typedef memory::unique_ptr<gsAdditiveOp> uPtr;

    /// Default Constructor
    gsAdditiveOp() : m_transfers(), m_ops() {}

    /// @brief Constructor
    ///
    /// The operator realizes \f$ \sum_{i=1}^n T_i A_i T_i^T \f$
    ///
    /// @param transfers  transfer matrices \f$ T_i \f$
    /// @param ops        local operators \f$ A_i \f$
    gsAdditiveOp(TransferContainer transfers, OpContainer ops)
    : m_transfers(), m_ops(give(ops))
    {
        const size_t sz = transfers.size();
        m_transfers.reserve(sz);
        for (size_t i=0; i<sz; ++i)
            m_transfers.push_back( transfers[i].moveToPtr() );
#ifndef NDEBUG
        GISMO_ASSERT( m_transfers.size() == m_ops.size(), "Sizes do not agree" );
        for (size_t i=0; i<sz; ++i)
        {
            GISMO_ASSERT ( m_transfers[i]->rows()==m_transfers[0]->rows()
                       && m_transfers[i]->cols() == m_ops[i]->rows()
                       && m_ops[i]->cols() == m_ops[i]->rows(),
                       "Dimensions of the operators do not fit." );
        }
#endif
    }
    
    /// @brief Constructor
    ///
    /// The operator realizes \f$ \sum_{i=1}^n T_i A_i T_i^T \f$
    ///
    /// @param transfers  transfer matrices \f$ T_i \f$
    /// @param ops        local operators \f$ A_i \f$
    gsAdditiveOp(TransferPtrContainer transfers, OpContainer ops)
    : m_transfers(give(transfers)), m_ops(give(ops))
    {
#ifndef NDEBUG
        GISMO_ASSERT( m_transfers.size() == m_ops.size(), "Sizes do not agree" );
        const size_t sz = m_transfers.size();
        for (size_t i=0; i<sz; ++i)
        {
            GISMO_ASSERT ( m_transfers[i]->rows()==m_transfers[0]->rows()
                       && m_transfers[i]->cols() == m_ops[i]->rows()
                       && m_ops[i]->cols() == m_ops[i]->rows(),
                       "Dimensions of the operators do not fit." );
        }
#endif
    }

    /// Make function
    ///
    /// This function allows to obtain an empty instance
    static uPtr make()
    { return uPtr( new gsAdditiveOp() ); }

    /// Make function
    ///
    /// The operator realizes \f$ \sum_{i=1}^n T_i A_i T_i^T \f$
    ///
    /// @param transfers  transfer matrices \f$ T_i \f$
    /// @param ops        local operators \f$ A_i \f$
    static uPtr make(TransferContainer transfers, OpContainer ops)
    { return uPtr( new gsAdditiveOp( give(transfers), give(ops) ) ); }

    /// Make function
    ///
    /// The operator realizes \f$ \sum_{i=1}^n T_i A_i T_i^T \f$
    ///
    /// @param transfers  transfer matrices \f$ T_i \f$
    /// @param ops        local operators \f$ A_i \f$
    static uPtr make(TransferPtrContainer transfers, OpContainer ops)
    { return uPtr( new gsAdditiveOp( give(transfers), give(ops) ) ); }


    /// Add another entry to the sum
    ///
    /// @param transfer   the additional transfer matrix \f$ T_i \f$
    /// @param op         the additional operator \f$ A_i \f$
    void addOperator(Transfer transfer, OpPtr op)
    {
        m_transfers.push_back(transfer.moveToPtr());
        m_ops.push_back(give(op));
        GISMO_ASSERT ( m_transfers.back()->rows()==m_transfers[0]->rows()
                       && m_transfers.back()->cols() == m_ops.back()->rows()
                       && m_ops.back()->cols() == m_ops.back()->rows(),
                       "Dimensions of the operators do not fit." );
    }
    
    /// Add another entry to the sum
    ///
    /// @param transfer   the additional transfer matrix \f$ T_i \f$
    /// @param op         the additional operator \f$ A_i \f$
    void addOperator(TransferPtr transfer, OpPtr op)
    {
        m_transfers.push_back(give(transfer));
        m_ops.push_back(give(op));
        GISMO_ASSERT ( m_transfers.back()->rows()==m_transfers[0]->rows()
                       && m_transfers.back()->cols() == m_ops.back()->rows()
                       && m_ops.back()->cols() == m_ops.back()->rows(),
                       "Dimensions of the operators do not fit." );
    }

    void apply(const gsMatrix<T>& input, gsMatrix<T>& x) const;

    index_t rows() const
    {
        GISMO_ASSERT( !m_transfers.empty(), "gsAdditiveOp::rows does not work for 0 operators." );
        return m_transfers[0]->rows();
    }

    index_t cols() const
    {
        GISMO_ASSERT( !m_transfers.empty(), "gsAdditiveOp::cols does not work for 0 operators." );
        return m_transfers[0]->rows();
    }

protected:
    TransferPtrContainer m_transfers;   ///< Transfer matrices
    OpContainer m_ops;                  ///< Operators to be applied in the subspaces

};

} // namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsAdditiveOp.hpp)
#endif
