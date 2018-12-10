/** @file gsAdditivePrecOp.h

    @brief Allows to set up additive Schwarz type smoothers

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
///    gsPreconditionerOp<>::Ptr pc = gsAdditivePrecOp<>::make( A, transfers, ops );
/// \endcode
///
/// is equivalent to
///
/// \code{.cpp}
///    gsSumOp<>::Ptr s = gsSumOp<>::make();
///
///    for (index_t i=0; i<transfers.size(); ++i)
///       s->addOperator(
///          gsProductOp<>::make(
///             makeMatrixOp(transfers[i].transpose()),
///             ops[i],
///             makeMatrixOp(transfers[i])
///          )
///       );
///     gsPreconditionerOp<>::Ptr pc = gsPreconditionerFromOp<>::make( A, s );
/// \endcode
///
/// but much faster.
///
/// @ingroup Solvers

template<class T>
class gsAdditivePrecOp : public gsPreconditionerOp<T>
{
    typedef gsPreconditionerOp<T> Base;
    typedef typename gsLinearOperator<T>::Ptr OpPtr;
public:

    /// Shared pointer
    typedef memory::shared_ptr<gsAdditivePrecOp> Ptr;

    /// Unique pointer
    typedef memory::unique_ptr<gsAdditivePrecOp> uPtr;

    /// Type of container of transfer matrices
    typedef std::vector< gsSparseMatrix<T,RowMajor> > TransferContainer;

    /// Type of container of linear operators
    typedef std::vector<OpPtr> OpContainer;

    /// Default Constructor
    gsAdditivePrecOp() : m_A(), m_transfers(), m_ops(), m_damping(1) {}

    /// Constructor
    gsAdditivePrecOp(const OpPtr& A,
                       TransferContainer transfers,
                       OpContainer ops,
                       T damping = 1.0)
    : m_A(A), m_transfers(give(transfers)), m_ops(give(ops)), m_damping(damping)
    {
        GISMO_ASSERT( m_ops.size() == m_transfers.size(), "Sizes do not agree" );
    }

    /// Make function
    static uPtr make(const OpPtr& A,
                     TransferContainer transfers,
                     OpContainer ops,
                     T damping = 1.0)
    { return uPtr( new gsAdditivePrecOp( A, give(transfers), give(ops), damping ) ); }

    /// Make function
    static uPtr make(const OpPtr& A,
                     std::pair<
                         TransferContainer,
                         OpContainer
                     > transfer_and_ops,
                     T damping = 1.0)
    { return uPtr( new gsAdditivePrecOp( A, give(transfer_and_ops.first), give(transfer_and_ops.second), damping ) ); }

    virtual void step(const gsMatrix<T>& f, gsMatrix<T>& x) const;

    index_t rows() const { return m_A->rows();}
    index_t cols() const { return m_A->cols();}
    OpPtr underlyingOp() const { return m_A; }

protected:
    OpPtr m_A;
    TransferContainer m_transfers;
    OpContainer m_ops;
    T m_damping;
    mutable gsMatrix<T> m_res, m_res_local, m_corr_local;

};

} // namespace gismo
