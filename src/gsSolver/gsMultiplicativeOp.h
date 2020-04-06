/** @file gsMultiplicativeOp.h

    @brief Allows to set up multiplicative Schwarz type preconditioners

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): S. Takacs, J. Sogn
*/

#pragma once

#include <gsSolver/gsPreconditioner.h>
#include <gsSolver/gsMatrixOp.h>

namespace gismo
{

/// @brief Generic preconditioner which applies an arbitrary linear operator to the residual.
///
/// @ingroup Solvers

template<class T>
class gsMultiplicativeOp GISMO_FINAL : public gsPreconditionerOp<T>
{
    typedef typename gsLinearOperator<T>::Ptr   OpPtr;
    typedef std::vector<OpPtr>   OpContainer;
    typedef gsSparseMatrix<T,RowMajor>   Transfer;
    typedef std::vector<Transfer>   TransferContainer;
    typedef gsVector<size_t>  TransferVector;
    typedef std::vector<TransferVector > TransferVectorContainer;

public:

    /// Shared pointer
    typedef memory::shared_ptr<gsMultiplicativeOp> Ptr;

    /// Unique pointer
    typedef memory::unique_ptr<gsMultiplicativeOp> uPtr;

    /// Default Constructor
    //gsMultiplicativeOp() : m_transfers(), m_ops() {}

    /// @brief Constructor
    ///
    /// The operator realizes \f$ \prod_{i=1}^n (I - T_i A_i T_i^T) \f$
    ///
    /// @param underlying  underlying matrix \f$ A \f$
    /// @param transfers   transfer matrices \f$ T_i \f$
    /// @param ops         local operators \f$ A_i \f$
    /// @note: This takes A by reference
    gsMultiplicativeOp(memory::shared_ptr< gsSparseMatrix<T> > underlying, TransferContainer transfers, OpContainer ops)
        : m_underlying(give(underlying)), m_transfers(give(transfers)), m_ops(give(ops)), m_optReady(false)
    {
#ifndef NDEBUG
        GISMO_ASSERT( m_underlying->rows() == m_underlying->cols(), "Dimensions of the operators do not fit." );
        GISMO_ASSERT( m_transfers.size() == m_ops.size(), "Sizes do not agree" );
        const size_t sz = m_transfers.size();
        for (size_t i=0; i<sz; ++i)
        {
            GISMO_ASSERT ( m_transfers[i].rows()==m_underlying->rows()
                       && m_transfers[i].cols() == m_ops[i]->rows()
                       && m_ops[i]->cols() == m_ops[i]->rows(),
                       "Dimensions of the operators do not fit."<<m_transfers[i].rows()<<"=="<<m_underlying->rows()<<
                       "&&"<<m_transfers[i].cols()<<"=="<< m_ops[i]->rows()<<
	               "&&"<<m_ops[i]->cols()<<"=="<<m_ops[i]->rows()<<"; i="<<i<<"; sz="<<sz );
        }
#endif
    }

    /// Make function
    ///
    /// The operator realizes \f$ \sum_{i=1}^n T_i A_i T_i^T \f$
    ///
    /// @param underlying  underlying matrix \f$ A \f$
    /// @param transfers  transfer matrices \f$ T_i \f$
    /// @param ops        local operators \f$ A_i \f$
    /// @note: This takes A by reference
    static uPtr make(memory::shared_ptr< gsSparseMatrix<T> > underlying, TransferContainer transfers, OpContainer ops)
    { return uPtr( new gsMultiplicativeOp( give(underlying), give(transfers), give(ops) ) ); }

    /// Make function
    ///
    /// The operator realizes \f$ \sum_{i=1}^n T_i A_i T_i^T \f$
    ///
    /// @param underlying  underlying matrix \f$ A \f$
    /// @note: This takes A by reference
    static uPtr make(memory::shared_ptr< gsSparseMatrix<T> > underlying)
    { return uPtr( new gsMultiplicativeOp( give(underlying), TransferContainer(), OpContainer() ) ); }

    /// Add another entry to the sum
    ///
    /// @param transfer   the additional transfer matrix \f$ T_i \f$
    /// @param op         the additional operator \f$ A_i \f$
    void addOperator(Transfer transfer, OpPtr op)
    {
        m_transfers.push_back(transfer); //TODO: should give
        m_ops.push_back(op); //TODO: should give
        GISMO_ASSERT ( transfer.rows()==m_underlying->rows()
                       && transfer.cols() == op->rows()
                       && op->cols() == op->rows(),
                       "Dimensions of the operators do not fit."<<transfer.rows()<<"=="<<m_underlying->rows()<<"&&"<<transfer.cols()<<"=="<<op->rows()<<"&&"<<op->cols()<<"=="<<op->rows() );
    }

    void step(const gsMatrix<T>& rhs, gsMatrix<T>& x) const;

    void stepT(const gsMatrix<T>& rhs, gsMatrix<T>& x) const;

private:
    
    void stepOpt(const gsMatrix<T>& rhs, gsMatrix<T>& x) const;

    void stepTOpt(const gsMatrix<T>& rhs, gsMatrix<T>& x) const;

    OpPtr underlyingOp() const
    {
        return makeMatrixOp(m_underlying);
    }

    index_t rows() const
    {
        return m_underlying->rows();
    }

    index_t cols() const
    {
        return m_underlying->cols();
    }
    
public:

    ///Initialize an optimization of the step() and stepT() function
    ///by assuming that the transformation matrix have the form (N,n)
    ///where each column has one non-zero entry with the value 1.
    ///Then createsa data-structure: std::vector< gsVector<index_t>>
    ///where number of rows in the vector is \a n and the i'th entry
    ///is the row index of the k'th transformation matrix.
    ///Need a column-major underlying matrix.
    void initOptimize()
    {
        //TODO: Check if it is already colmajor, then this can be skipped.
        m_underlyingColMajorMat = *m_underlying;
        
        const size_t sz = m_transfers.size();
        
        for (size_t k = 0; k < sz; ++k)
        {
            size_t localSz = m_transfers[k].cols();             
            TransferVector transVec;
            transVec.setZero(localSz);
            for (int  j=0; j < m_transfers[k].outerSize(); ++j)
            {
                for (typename gsSparseMatrix<T,RowMajor>::InnerIterator it(m_transfers[k],j); it; ++it) 
                {
                    transVec(it.col()) = it.row();
                }
            }
            m_transfersVec.push_back(transVec);
        }

        m_optReady = true;
    }

protected:
    memory::shared_ptr< gsSparseMatrix<T> > m_underlying;    
    TransferContainer m_transfers;   ///< Transfer matrices
    OpContainer m_ops;               ///< Operators to be applied in the subspaces

    //Used for the optimization
    bool m_optReady; ///< Optimization switch
    TransferVectorContainer m_transfersVec; ///< Different format of m_transfers for optimization
    gsSparseMatrix<T,ColMajor> m_underlyingColMajorMat; ///< Underlying sparse matrix stored as colmajor

};

} // namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsMultiplicativeOp.hpp)
#endif
