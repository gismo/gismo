/** @file gsBlockOp.hpp

    @brief Simple class create a block preconditioner structure.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): J. Sogn
*/

namespace gismo
{

template<typename T>
gsBlockOp<T>::gsBlockOp(index_t nRows, index_t nCols)
{
    m_blockPrec.resize(nRows, nCols);
    m_blockTargetPositions.setZero(nRows);
    m_blockInputPositions.setZero(nCols);
    // Fill up all block entries with null pointers.
    for (index_t i = 0; i < nRows; ++i)
        for (index_t j = 0; j < nCols; ++j)
            m_blockPrec(i,j).reset();
}

template<typename T>
void gsBlockOp<T>::addOperator(index_t row, index_t col, const BasePtr& op)
{
    GISMO_ASSERT( row >= 0 && row < m_blockPrec.rows(), "The given row is not feasible." );
    GISMO_ASSERT( col >= 0 && col < m_blockPrec.cols(), "The given column is not feasible." );
    GISMO_ASSERT( op->rows() == m_blockTargetPositions[row] || m_blockTargetPositions[row] == 0,
                  "The size of the given preconditioner does not fit to the other preconditioners in the same row." );
    GISMO_ASSERT( op->cols() == m_blockInputPositions[col] || m_blockInputPositions[col] == 0,
                  "The size of the given preconditioner does not fit to the other preconditioners in the same column." );

    m_blockPrec(row, col) = op;
    m_blockTargetPositions[row] = op->rows();
    m_blockInputPositions[col] = op->cols();
}


template<typename T>
void gsBlockOp<T>::apply(const gsMatrix<T> & input, gsMatrix<T> & result) const
{
    result.setZero(m_blockTargetPositions.sum(), input.cols());
    gsVector<index_t> singleCol(1);
    singleCol <<  input.cols();
    typename gsMatrix<T>::BlockView resultBlocks = result.blockView(m_blockTargetPositions, singleCol);

    for (index_t i = 0; i < m_blockPrec.rows() ; ++i)
    {
        index_t inputIndex = 0;
        for (index_t j = 0; j < m_blockPrec.cols(); ++j)
        {
            if (!m_blockPrec(i,j))// if the block is a null pointer
            {
                inputIndex += m_blockInputPositions(j);
                continue;
            }

            gsMatrix<T> tmp_result;
            m_blockPrec(i,j)->apply(input.block(inputIndex,0,m_blockInputPositions(j),input.cols()),tmp_result);
            resultBlocks(i) += tmp_result;
            inputIndex += m_blockInputPositions(j);
        }
    }
}

}
