/** @file gsBlockOp.h

    @brief Simple class create a block preconditioner structure.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): J. Sogn
*/
#pragma once

#include <gsCore/gsExport.h>
#include <gsCore/gsLinearAlgebra.h>
#include <gsSolver/gsLinearOperator.h>

namespace gismo
{

/** \brief Simple class create a block operator structure.
 *
 * This class represents a linear operator \f$C\f$ having block structure:
 * \f[
 *   C =
         \left( \begin{array}{cccc}
         C_{00} & C_{01} & \ldots & C_{0m}  \\
         C_{10} & C_{11} & \ldots & C_{1m}  \\
         \vdots & \vdots & \ddots & \vdots  \\
         C_{n0} & C_{n1} & \ldots & C_{nm}
         \end{array}
         \right),
   \f]
 * where \f$C_{ij}\f$ are themselves gsLinearOperators.
 *
 * The number of blocks (m and n) are specified in the constructor. The blocks \f$C_{ij}\f$ are
 * defined using addOperator(i,j,...). Unspecified blocks are considered to be 0.
 *
 * \ingroup Solver
 */
template<class T>
class gsBlockOp GISMO_FINAL : public gsLinearOperator<T>
{
public:

    /// Shared pointer for gsBlockOp
    typedef memory::shared_ptr< gsBlockOp<T> > Ptr;

    /// Unique pointer for gsBlockOp
    typedef memory::unique_ptr< gsBlockOp<T> > uPtr;

    /// Base class
    typedef memory::shared_ptr< gsLinearOperator<T> > BasePtr;

    /// Constructor. Takes the number of blocks (nRows, nCols). Provide the contents of the blocks with addOperator
    gsBlockOp(index_t nRows, index_t nCols);

    /// Make function returning a smart pointer
    static uPtr make(index_t nRows, index_t nCols)
    { return memory::make_unique( new gsBlockOp(nRows,nCols) ); }

    /**
     * @brief Add a preconditioner \f$C_{ij}\f$ to the block structure
     * @param row row position in the block operator
     * @param col column position in the block operator
     * @param op shared pointer to the operator
     */
    void addOperator(index_t row, index_t col, const BasePtr& op);

    /**
    * @brief Returns the pointer to a linear operator of a specific block (if existent)
    * @param row row position in the block operator
    * @param col column position in the block operator
    *
    * @note  The result can be a null-pointer
    */
    const BasePtr & getOperator(index_t row, index_t col) const {
        return m_blockPrec(row,col);
    }

    /**
     * @brief Apply the correct segment of the input vector on the preconditioners in the block structure and store the result.
     * @param input  Input vector
     * @param result Result vector
     */
    void apply(const gsMatrix<T> & input, gsMatrix<T> & result) const;

    /// Number of row blocks
    index_t rowBlocks() const {return m_blockPrec.rows();}
    /// Number of col blocks
    index_t colBlocks() const {return m_blockPrec.cols();}

    index_t rows() const {return m_blockTargetPositions.sum();}
    index_t cols() const {return m_blockInputPositions.sum() ;}

private:

    Eigen::Array<BasePtr, Dynamic, Dynamic> m_blockPrec;

    //Contains the size of the target vector for each block
    gsVector<index_t> m_blockTargetPositions;
    //Contains the size of the input vector for each block
    gsVector<index_t> m_blockInputPositions;

};

} // namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsBlockOp.hpp)
#endif
