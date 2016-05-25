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
    blockPrec(row, col) = prec;
    blockTargetPositions[row] = prec->rows();
    blockInputPositions[col] = prec->cols();
    if (!consistencyCheck())
        GISMO_ERROR("Block preconditioners do not have correct dimension");
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


bool gsBlockPreconditioner::consistencyCheck()
{
    for (index_t i = 0; i < blockPrec.rows(); ++i)
    {
        if (!blockPrec(i,0))// if the block is a null pointer
            continue;

        index_t r = blockPrec(i,0)->rows();

        for (index_t j = 1; j < blockPrec.cols(); ++j)
        {
            if (!blockPrec(i,j))// if the block is a null pointer
                continue;

            if (blockPrec(i,j)->rows() != r)
            {
                gsWarn << "Block preconditioners do not have correct dimension" << std::endl;
                return false;
            }
        }
    }
    for (index_t j = 0; j < blockPrec.cols(); ++j)
    {
        if (!blockPrec(0,j))// if the block is a null pointer
            continue;

        index_t c = blockPrec(0,j)->cols();

        for (index_t i = 1; i < blockPrec.rows(); ++i)
        {
            if (!blockPrec(i,j))// if the block is a null pointer
                continue;

            if (blockPrec(i,j)->rows() != c)
            {
                gsWarn << "Block preconditioners do not have correct dimension" << std::endl;
                return false;
            }
        }
    }

    return true;
}


}

