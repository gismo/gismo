/** @file gsBlockPreconditioner.h

    @brief Simple class create a block preconditioner structure.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): J. Sogn
*/
#pragma once

#include <gsCore/gsLinearAlgebra.h>
#include <gsSolver/gsPreconditioner.h>

namespace gismo
{

/** \brief Simple class create a block preconditioner structure.
 *
 * Let \f$C\f$ be a preconditioner for the system of equations \f$ A\mathbf{x} =  \mathbf{f}\f$.
 * We instead wish to solve the preconditioned system \f$ CA\mathbf{x} =  C\mathbf{f}\f$.
 *
 * This class allows \f$C\f$ to be a block structure of preconditioners i.e \f$C\f$:
 * \f[
     \left( \begin{array}{cccc}
     C_{00} & C_{01} & \ldots & C_{0n}  \\
     C_{10} & C_{11} & \ldots & C_{1n}  \\
     \vdots & \vdots & \ddots & \vdots  \\
     C_{n0} & C_{n1} & \ldots & C_{nn}
     \end{array}
     \right)
     \f]
 * Where \f$C_{ij}\f$ are pre-defined preconditioners which all have the apply(input, dest) function defined.
 */
class gsBlockPreconditioner : public gsPreconditioner
{
public:
    gsBlockPreconditioner(index_t nRows, index_t nCols)
    {
        blockPrec.resize(nRows, nCols);
        blockTargetPositions.setZero(nCols);
        blockInputPositions.setZero(nRows);
        gsPreconditioner * null_ptr = NULL;
        // Fill up all block entries with null pointers.
        for (index_t i = 0; i< nRows; ++i)
            for (index_t j = 0; j< nCols; ++j)
                blockPrec(i,j) = null_ptr;
    }

    /**
     * @brief Add a preconditioner \f$C_{ij}\f$ to the block structure
     * @param prec Pointer the preconditioner
     * @param row row position in the block preconditioner
     * @param col column position in the block preconditioner
     */
    void addPreconditioner(gsPreconditioner * prec, index_t row, index_t col)
    {
        blockPrec(row, col) = prec;
        blockTargetPositions[row] = prec->rows();
        blockInputPositions[col] = prec->cols();
        if (!ConsistencyCheck())
            GISMO_ERROR("Block preconditioners do not have correct dimension");
    }

    /**
     * @brief Apply the correct segment of the input vector on the preconditioners in the block structure and store the result.
     * @param input  Input vector
     * @param result Result vector
     */
    void apply(const gsMatrix<real_t> & input, gsMatrix<real_t> & result) const
    {
        result.setZero(blockTargetPositions.sum(), input.cols());
        gsVector<index_t> singleCol(1);
        singleCol << 1;
        gsMatrix<real_t>::BlockView resultBlocks= result.blockView(blockTargetPositions, singleCol);

        for (index_t i = 0; i< blockPrec.rows() ; ++i)
        {
            index_t inputIndex = 0;
            for (index_t j = 0; j< blockPrec.cols(); ++j)
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

    index_t rows() const {return blockTargetPositions.sum();}
    index_t cols() const {return blockInputPositions.sum() ;}

private:

    /**
     * @brief Loops through every preconditioner and checks that the dimensions are correct
     */
    bool ConsistencyCheck()
    {
        for (index_t i = 0; i< blockPrec.rows(); ++i)
        {
            if (!blockPrec(i,0))// if the block is a null pointer
                    continue;

            index_t r = blockPrec(i,0)->rows();

            for (index_t j = 1; j< blockPrec.cols(); ++j)
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
        for (index_t j = 0; j< blockPrec.cols(); ++j)
        {
            if (!blockPrec(0,j))// if the block is a null pointer
                    continue;

            index_t c = blockPrec(0,j)->cols();

            for (index_t i = 1; i< blockPrec.rows(); ++i)
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

    Eigen::Array<gsPreconditioner *, Dynamic, Dynamic> blockPrec;

    //Contains the size of the target vector for each block
    gsVector<index_t> blockTargetPositions;
    //Contains the size of the input vector for each block
    gsVector<index_t> blockInputPositions;

};

} // namespace gismo
