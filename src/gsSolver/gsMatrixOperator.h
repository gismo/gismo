/** @file gsMatrixOperator.h

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
class gsMatrixOperator : public gsLinearOperator
{
public:

    /// Shared pointer for gsMatrixOperator
    typedef memory::shared_ptr<gsMatrixOperator> Ptr;

    /// Unique pointer for gsMatrixOperator   
    typedef typename memory::unique<gsMatrixOperator>::ptr uPtr;
    
    
    gsMatrixOperator(const MatrixType& mat, bool sym=false)
        : m_mat(mat), m_symmetric(sym) {}
        
    static Ptr make(const MatrixType& mat, bool sym=false) { return shared( new gsMatrixOperator(mat,sym) ); }

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

/** This essentially just calls the gsMatrixOperator constructor, but the
 * use of a template functions allows us to let the compiler do type inference,
 * so we don't need to type out the matrix type explicitly.
 */
template <class MatrixType>
typename gsMatrixOperator<MatrixType>::Ptr makeMatrixOperator(const MatrixType& mat, bool sym=false)
{
    return shared(new gsMatrixOperator<MatrixType>(mat, sym));
}


/// @brief Simple adapter class to use the transpose of a matrix (or matrix-like object) as
/// a linear operator. This should, of course, be done without transposing the matrix itself.
///
/// \ingroup Solver

template <class MatrixType>
class gsTransposedMatrixOperator : public gsLinearOperator
{
public:

    /// Shared pointer for gsTransposedMatrixOperator
    typedef memory::shared_ptr<gsTransposedMatrixOperator> Ptr;

    /// Unique pointer for gsTransposedMatrixOperator   
    typedef typename memory::unique<gsTransposedMatrixOperator>::ptr uPtr;
    
    gsTransposedMatrixOperator(const MatrixType& mat)
        : m_mat(mat) {}

    static Ptr make(const MatrixType& mat) { return shared( new gsTransposedMatrixOperator(mat) ); }
        
    void apply(const gsMatrix<real_t> & input, gsMatrix<real_t> & x) const
    {
        x.noalias() = m_mat.transpose() * input;
    }

    index_t rows() const {return m_mat.cols();}

    index_t cols() const {return m_mat.rows();}
    
private:
    const MatrixType& m_mat;
};

/** This essentially just calls the gsTransposedMatrixOperator constructor, but the
 * use of a template functions allows us to let the compiler do type inference,
 * so we don't need to type out the matrix type explicitly.
 */
template <class MatrixType>
typename gsTransposedMatrixOperator<MatrixType>::Ptr makeTransposedMatrixOperator(const MatrixType& mat)
{
    return shared(new gsTransposedMatrixOperator<MatrixType>(mat));
}


/// @brief Simple adapter class to use an Eigen solver (having a compute() and a solve() method) as a linear operator.
///
/// \ingroup Solver

template <class SolverType>
class gsSolverOperator : public gsLinearOperator
{
public:

    /// Shared pointer for gsBasis
    typedef memory::shared_ptr<gsSolverOperator> Ptr;

    /// Unique pointer for gsBasis   
    typedef typename memory::unique<gsSolverOperator>::ptr uPtr;
    
    template <class MatrixType>
    gsSolverOperator(const MatrixType& mat)
    {
        GISMO_ASSERT(mat.rows() == mat.cols(), "Need square matrix");
        m_size = mat.rows();

        m_solver.compute(mat);
    }

    template <class MatrixType>
    static Ptr make(const MatrixType& mat) { return shared( new gsSolverOperator(mat) ); }    
    
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
typename gsSolverOperator< Eigen::PartialPivLU< Eigen::Matrix<T, _Rows, _Cols, _Opt> > >::Ptr  makePartialPivLUSolver(const gsMatrix<T, _Rows, _Cols, _Opt> & mat)
{
    return shared( new gsSolverOperator< Eigen::PartialPivLU< Eigen::Matrix<T, _Rows, _Cols, _Opt> > >(mat) );
}

/// @brief Convenience function to create an LU solver with full pivoting (for dense matrices) as a gsLinearOperator.
template <class T, int _Rows, int _Cols, int _Opt>
typename gsSolverOperator< Eigen::FullPivLU< Eigen::Matrix<T, _Rows, _Cols, _Opt> > >::Ptr  makeFullPivLUSolver(const gsMatrix<T, _Rows, _Cols, _Opt> & mat)
{
    return shared( new gsSolverOperator< Eigen::FullPivLU< Eigen::Matrix<T, _Rows, _Cols, _Opt> > >(mat) );
}


/// @brief Convenience function to create a Cholesky (LDL^T) solver (for dense matrices) as a gsLinearOperator.
///
/// @note Works only on symmetric (stored in lower half) and positive (semi-)definite matrices.
template <class T, int _Rows, int _Cols, int _Opt>
typename gsSolverOperator< Eigen::LDLT< Eigen::Matrix<T, _Rows, _Cols, _Opt> > >::Ptr  makeCholeskySolver(const gsMatrix<T, _Rows, _Cols, _Opt> & mat)
{
    return shared( new gsSolverOperator< Eigen::LDLT< Eigen::Matrix<T, _Rows, _Cols, _Opt> > >(mat) );
}


/// @brief Convenience function to create a sparse LU solver as a gsLinearOperator.
///
/// @note This uses the default COLAMD column ordering.
template <typename T, int _Opt, typename _Index>
typename gsSolverOperator< typename gsSparseSolver<T>::LU >::Ptr  makeSparseLUSolver(const gsSparseMatrix<T,_Opt,_Index> & mat)
{
    return shared( new gsSolverOperator< typename gsSparseSolver<T>::LU >(mat) );
}


/// @brief Convenience function to create a sparse Cholesky (simplicial LDL^T) solver as a gsLinearOperator.
///
/// @note Works only on sparse, symmetric (stored in lower half) and positive definite matrices.
template <typename T, int _Opt, typename _Index>
typename gsSolverOperator< typename gsSparseSolver<T>::SimplicialLDLT >::Ptr  makeSparseCholeskySolver(const gsSparseMatrix<T,_Opt,_Index> & mat)
{
    return shared( new gsSolverOperator<typename  gsSparseSolver<T>::SimplicialLDLT >(mat) );
}


} // namespace gismo
