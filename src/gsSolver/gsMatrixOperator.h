/** @file gsMatrixOp.h

    @brief Simple adapter classes to use matrices or linear solvers as gsLinearOperators

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): J. Sogn, C. Hofreither
*/
#pragma once

#include <gsCore/gsLinearAlgebra.h>
#include <gsSolver/gsLinearOperator.h>

namespace gismo
{

/// @brief Simple adapter class to use a matrix (or matrix-like object) as a linear operator. Needed for the iterative method classes.
///
/// \ingroup Solver

template <class MatrixType>
class gsMatrixOp : public gsLinearOperator
{
public:

    /// Shared pointer for gsMatrixOp
    typedef memory::shared_ptr<gsMatrixOp> Ptr;

    /// Unique pointer for gsMatrixOp   
    typedef typename memory::unique<gsMatrixOp>::ptr uPtr;
    
    
    gsMatrixOp(const MatrixType& mat, bool sym=false)
        : m_mat(mat), m_symmetric(sym) {}
        
    static Ptr make(const MatrixType& mat, bool sym=false) { return shared( new gsMatrixOp(mat,sym) ); }

    void apply(const gsMatrix<real_t> & input, gsMatrix<real_t> & x) const
    {
        if (m_symmetric)
            x.noalias() = m_mat.template selfadjointView<Lower>() * input;
        else
            x.noalias() = m_mat * input;
    }

    index_t rows() const {return m_mat.rows();}

    index_t cols() const {return m_mat.cols();}

    ///Returns the matrix
    const MatrixType& matrix() const { return m_mat; }

private:
    const MatrixType& m_mat;
    bool m_symmetric;
};

/** This essentially just calls the gsMatrixOp constructor, but the
 * use of a template functions allows us to let the compiler do type inference,
 * so we don't need to type out the matrix type explicitly.
 */
template <class MatrixType>
typename gsMatrixOp<MatrixType>::Ptr makeMatrixOp(const MatrixType& mat, bool sym=false)
{
    return shared(new gsMatrixOp<MatrixType>(mat, sym));
}


/// @brief Simple adapter class to use the transpose of a matrix (or matrix-like object) as
/// a linear operator. This should, of course, be done without transposing the matrix itself.
///
/// \ingroup Solver

template <class MatrixType>
class gsTransposedMatrixOp : public gsLinearOperator
{
public:

    /// Shared pointer for gsTransposedMatrixOp
    typedef memory::shared_ptr<gsTransposedMatrixOp> Ptr;

    /// Unique pointer for gsTransposedMatrixOp   
    typedef typename memory::unique<gsTransposedMatrixOp>::ptr uPtr;
    
    gsTransposedMatrixOp(const MatrixType& mat)
        : m_mat(mat) {}

    static Ptr make(const MatrixType& mat) { return shared( new gsTransposedMatrixOp(mat) ); }
        
    void apply(const gsMatrix<real_t> & input, gsMatrix<real_t> & x) const
    {
        x.noalias() = m_mat.transpose() * input;
    }

    index_t rows() const {return m_mat.cols();}

    index_t cols() const {return m_mat.rows();}
    
private:
    const MatrixType& m_mat;
};

/** This essentially just calls the gsTransposedMatrixOp constructor, but the
 * use of a template functions allows us to let the compiler do type inference,
 * so we don't need to type out the matrix type explicitly.
 */
template <class MatrixType>
typename gsTransposedMatrixOp<MatrixType>::Ptr makeTransposedMatrixOp(const MatrixType& mat)
{
    return shared(new gsTransposedMatrixOp<MatrixType>(mat));
}


/// @brief Simple adapter class to use an Eigen solver (having a compute() and a solve() method) as a linear operator.
///
/// \ingroup Solver

template <class SolverType>
class gsSolverOp : public gsLinearOperator
{
public:

    /// Shared pointer for gsBasis
    typedef memory::shared_ptr<gsSolverOp> Ptr;

    /// Unique pointer for gsBasis   
    typedef typename memory::unique<gsSolverOp>::ptr uPtr;
    
    template <class MatrixType>
    gsSolverOp(const MatrixType& mat)
    {
        GISMO_ASSERT(mat.rows() == mat.cols(), "Need square matrix");
        m_size = mat.rows();

        m_solver.compute(mat);
    }

    template <class MatrixType>
    static Ptr make(const MatrixType& mat) { return shared( new gsSolverOp(mat) ); }    
    
    void apply(const gsMatrix<real_t> & input, gsMatrix<real_t> & x) const
    {
        x = m_solver.solve(input);
    }

    index_t rows() const { return m_size; }

    index_t cols() const { return m_size; }

    /// Access the solver class
    SolverType& solver()                { return m_solver; }

    /// Access the solver class
    const SolverType& solver() const    { return m_solver; }

private:
    SolverType m_solver;
    index_t m_size;
};


/// @brief Convenience function to create an LU solver with partial pivoting (for dense matrices) as a gsLinearOperator.
template <class T, int _Rows, int _Cols, int _Opt>
typename gsSolverOp< Eigen::PartialPivLU< Eigen::Matrix<T, _Rows, _Cols, _Opt> > >::Ptr  makePartialPivLUSolver(const gsMatrix<T, _Rows, _Cols, _Opt> & mat)
{
    return shared( new gsSolverOp< Eigen::PartialPivLU< Eigen::Matrix<T, _Rows, _Cols, _Opt> > >(mat) );
}

/// @brief Convenience function to create an LU solver with full pivoting (for dense matrices) as a gsLinearOperator.
template <class T, int _Rows, int _Cols, int _Opt>
typename gsSolverOp< Eigen::FullPivLU< Eigen::Matrix<T, _Rows, _Cols, _Opt> > >::Ptr  makeFullPivLUSolver(const gsMatrix<T, _Rows, _Cols, _Opt> & mat)
{
    return shared( new gsSolverOp< Eigen::FullPivLU< Eigen::Matrix<T, _Rows, _Cols, _Opt> > >(mat) );
}


/// @brief Convenience function to create a Cholesky (LDL^T) solver (for dense matrices) as a gsLinearOperator.
///
/// @note Works only on symmetric (stored in lower half) and positive (semi-)definite matrices.
template <class T, int _Rows, int _Cols, int _Opt>
typename gsSolverOp< Eigen::LDLT< Eigen::Matrix<T, _Rows, _Cols, _Opt> > >::Ptr  makeCholeskySolver(const gsMatrix<T, _Rows, _Cols, _Opt> & mat)
{
    return shared( new gsSolverOp< Eigen::LDLT< Eigen::Matrix<T, _Rows, _Cols, _Opt> > >(mat) );
}


/// @brief Convenience function to create a sparse LU solver as a gsLinearOperator.
///
/// @note This uses the default COLAMD column ordering.
template <typename T, int _Opt, typename _Index>
typename gsSolverOp< typename gsSparseSolver<T>::LU >::Ptr  makeSparseLUSolver(const gsSparseMatrix<T,_Opt,_Index> & mat)
{
    return shared( new gsSolverOp< typename gsSparseSolver<T>::LU >(mat) );
}


/// @brief Convenience function to create a sparse Cholesky (simplicial LDL^T) solver as a gsLinearOperator.
///
/// @note Works only on sparse, symmetric (stored in lower half) and positive definite matrices.
template <typename T, int _Opt, typename _Index>
typename gsSolverOp< typename gsSparseSolver<T>::SimplicialLDLT >::Ptr  makeSparseCholeskySolver(const gsSparseMatrix<T,_Opt,_Index> & mat)
{
    return shared( new gsSolverOp<typename  gsSparseSolver<T>::SimplicialLDLT >(mat) );
}


} // namespace gismo
