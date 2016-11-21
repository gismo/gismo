/** @file gsBlockOp.hpp

    @brief Simple class create a block preconditioner structure.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): J. Sogn
*/
#pragma once

#include <gsSolver/gsBlockOp.h>

namespace gismo
{

template<typename T>
gsBlockOp<T>::gsBlockOp(index_t nRows, index_t nCols)
{
    blockPrec.resize(nRows, nCols);
    blockTargetPositions.setZero(nRows);
    blockInputPositions.setZero(nCols);
    // Fill up all block entries with null pointers.
    for (index_t i = 0; i < nRows; ++i)
        for (index_t j = 0; j < nCols; ++j)
            blockPrec(i,j).reset();
}

template<typename T>
void gsBlockOp<T>::addOperator(index_t row, index_t col, const BasePtr& op)
{
    GISMO_ASSERT( row >= 0 && row < blockPrec.rows(), "The given row is not feasible." );
    GISMO_ASSERT( col >= 0 && col < blockPrec.cols(), "The given column is not feasible." );
    GISMO_ASSERT( op->rows() == blockTargetPositions[row] || blockTargetPositions[row] == 0,
                  "The size of the given preconditioner does not fit to the other preconditioners in the same row." );
    GISMO_ASSERT( op->cols() == blockInputPositions[col] || blockInputPositions[col] == 0,
                  "The size of the given preconditioner does not fit to the other preconditioners in the same column." );
    
    blockPrec(row, col) = op;
    blockTargetPositions[row] = op->rows();
    blockInputPositions[col] = op->cols();
}


template<typename T>
void gsBlockOp<T>::apply(const gsMatrix<T> & input, gsMatrix<T> & result) const
{
    result.setZero(blockTargetPositions.sum(), input.cols());
    gsVector<index_t> singleCol(1);
    singleCol <<  input.cols();
    typename gsMatrix<T>::BlockView resultBlocks = result.blockView(blockTargetPositions, singleCol);

    for (index_t i = 0; i < blockPrec.rows() ; ++i)
    {
        index_t inputIndex = 0;
        for (index_t j = 0; j < blockPrec.cols(); ++j)
        {
            if (!blockPrec(i,j))// if the block is a null pointer
            {
                inputIndex += blockInputPositions(j);
                continue;
            }

            gsMatrix<T> tmp_result;
            blockPrec(i,j)->apply(input.block(inputIndex,0,blockInputPositions(j),input.cols()),tmp_result);
            resultBlocks(i) += tmp_result;
            inputIndex += blockInputPositions(j);
        }
    }
}

template<typename T>
const typename gsBlockOp<T>::BasePtr & gsBlockOp<T>::getOperator(index_t row, index_t col) const
{
    GISMO_ASSERT((bool)blockPrec(row, col)!=0, "No linear operator exists in this block");
    return blockPrec(row,col);
}

}

