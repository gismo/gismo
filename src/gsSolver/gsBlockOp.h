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

/** \brief Simple class to create a block operator structure.
 *
 * This class represents a linear operator \f$C\f$ having block structure:
 * \f[
 *   C =
         \begin{pmatrix}
         C_{00} & C_{01} & \ldots & C_{0m}  \\
         C_{10} & C_{11} & \ldots & C_{1m}  \\
         \vdots & \vdots & \ddots & \vdots  \\
         C_{n0} & C_{n1} & \ldots & C_{nm}
         \end{pmatrix},
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

/** \brief Represents the solution of a block system
 *
 * Represents the inverse of
 * \f[
 *   M =
         \begin{pmatrix}
         A & C^T \\
         C & 0
         \end{pmatrix},
   \f]
 * where the user provides \f$ A^{-1} \f$ as linear operator representing the
 * inverse of a symmetric and positive definite matrix A and \f$ C \f$ as a
 * matrix.
 *
 * \ingroup Solver
 */
template<class T=real_t> //TODO: forward decs
class gsBlockSolverOp GISMO_FINAL : public gsLinearOperator<T>
{
public:

    /// Shared pointer for gsBlockSolverOp
    typedef memory::shared_ptr< gsBlockSolverOp<T> > Ptr;

    /// Unique pointer for gsBlockSolverOp
    typedef memory::unique_ptr< gsBlockSolverOp<T> > uPtr;

    /// Base class
    typedef memory::shared_ptr< gsLinearOperator<T> > BasePtr;

    /// Constructor. Takes the number of blocks (nRows, nCols). Provide the contents of the blocks with addOperator
    gsBlockSolverOp(BasePtr Ainv, gsSparseMatrix<T> C)
        : m_Ainv(give(Ainv)), m_C(give(C))
    {
        GISMO_ASSERT(m_Ainv->rows()==m_Ainv->cols(), "Ainv is assumed to be symmetric.");
        GISMO_ASSERT(m_Ainv->cols()==m_C.cols(), "Dimensions do not agree: "<<m_Ainv->cols()<<"!="<<m_C.cols());

        gsMatrix<T> Ct = m_C.transpose();
        gsMatrix<T> res;
        m_Ainv->apply(Ct, res);
        gsMatrix<T> S = Ct.transpose() * res;
        m_Sinv = makeCholeskySolver(S);
    }

    /// Make function returning a smart pointer
    static uPtr make(BasePtr Ainv, gsSparseMatrix<T> C)
    { return memory::make_unique( new gsBlockSolverOp(give(Ainv), give(C)) ); }

    void apply(const gsMatrix<T> & input, gsMatrix<T> & result) const
    {
        result.setZero(input.rows(), input.cols());

        gsMatrix<T> q;
        m_Ainv->apply(input.topRows(m_Ainv->rows()), q);

        gsMatrix<T> z = m_C * q - input.bottomRows(m_C.rows());
        gsMatrix<T> Sinvz;
        m_Sinv->apply(z, Sinvz);

        gsMatrix<T> AinvCSinvz;
        m_Ainv->apply(m_C.transpose() * Sinvz, AinvCSinvz);

        result.topRows(m_Ainv->rows()) = q - AinvCSinvz;
        result.bottomRows(m_C.rows()) = Sinvz;
    }

    index_t rows() const { return m_Ainv->rows()+m_C.rows(); }
    index_t cols() const { return m_Ainv->rows()+m_C.rows(); }

private:

    BasePtr m_Ainv, m_Sinv;
    gsSparseMatrix<T> m_C;

};



} // namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsBlockOp.hpp)
#endif
