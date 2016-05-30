/** @file gsBlockPreconditioner.cpp

    @brief Simple class create a block preconditioner structure.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): J. Sogn
*/
#include <gsSolver/gsBlockPreconditioner.h>

namespace gismo
{


gsBlockPreconditioner::gsBlockPreconditioner(index_t nRows, index_t nCols)
{
    blockPrec.resize(nRows, nCols);
    blockTargetPositions.setZero(nCols);
    blockInputPositions.setZero(nRows);
    // Fill up all block entries with null pointers.
    for (index_t i = 0; i < nRows; ++i)
        for (index_t j = 0; j < nCols; ++j)
            blockPrec(i,j).reset();
}


void gsBlockPreconditioner::addPreconditioner(const BasePtr& prec, index_t row, index_t col)
{
    GISMO_ASSERT( row >= 0 && row < blockPrec.rows(), "The given row is not feasible." );
    GISMO_ASSERT( col >= 0 && col < blockPrec.cols(), "The given column is not feasible." );
    GISMO_ASSERT( prec->rows() == blockTargetPositions[row] || blockTargetPositions[row] == 0,
                  "The size of the given preconditioner does not fit to the other preconditioners in the same row." );
    GISMO_ASSERT( prec->cols() == blockTargetPositions[col] || blockTargetPositions[col] == 0,
                  "The size of the given preconditioner does not fit to the other preconditioners in the same column." );
    
    blockPrec(row, col) = prec;
    blockTargetPositions[row] = prec->rows();
    blockInputPositions[col] = prec->cols();
}


void gsBlockPreconditioner::apply(const gsMatrix<real_t> & input, gsMatrix<real_t> & result) const
{
    result.setZero(blockTargetPositions.sum(), input.cols());
    gsVector<index_t> singleCol(1);
    singleCol << 1;
    gsMatrix<real_t>::BlockView resultBlocks= result.blockView(blockTargetPositions, singleCol);

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

            gsMatrix<real_t> tmp_result;
            blockPrec(i,j)->apply(input.block(inputIndex,0,blockInputPositions(j),input.cols()),tmp_result);
            resultBlocks(i) += tmp_result;
            inputIndex += blockInputPositions(j);
        }
    }
}

}

