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

/**
  * @brief Simple adapter class to use a matrix (or matrix-like object) as a linear operator. Needed for the iterative method classes.
  *
  * \ingroup Solver
  */
  
template <class Derived>
class gsMatrixOp : public gsLinearOperator<typename Derived::Scalar>
{
public:
    typedef typename Derived::Scalar T;
    
    /// Shared pointer for gsMatrixOp
    typedef typename memory::shared<gsMatrixOp>::ptr Ptr;

    /// Unique pointer for gsMatrixOp   
    typedef typename memory::unique<gsMatrixOp>::ptr uPtr;
    
    // The matrix type
    typedef typename Eigen::EigenBase<Derived> MatrixType;

    /// Shared pointer to the matrix type
    typedef typename memory::shared<MatrixType>::ptr MatrixPtr;
    
    /// @brief Constructor taking a reference
    /// @note This does not copy the matrix. Make sure that the matrix is not deleted too early or provide a shared pointer.
    gsMatrixOp(const MatrixType& mat, bool sym=false) : m_mat(memory::make_shared_not_owned(&mat.derived())), m_symmetric(sym)    {}
    
    /// Constructor taking a shared pointer
    gsMatrixOp(const MatrixPtr& mat, bool sym=false) : m_mat(mat->derived()), m_symmetric(sym)    {}

    /// @brief Make function returning a smart pointer
    /// @note This does not copy the matrix. Make sure that the matrix is not deleted too early or provide a shared pointer.
    static Ptr make(const MatrixType& mat, bool sym=false)
        { return memory::make_shared( new gsMatrixOp(mat,sym) ); }

    /// Make function returning a smart pointer
    static Ptr make(const MatrixPtr& mat, bool sym=false)
        { return memory::make_shared( new gsMatrixOp(mat,sym) ); }

    void apply(const gsMatrix<T> & input, gsMatrix<T> & x) const
    {
        if (m_symmetric)
            x.noalias() = (*m_mat).template selfadjointView<Lower>() * input;
        else
            x.noalias() = (*m_mat) * input;
    }

    index_t rows() const {return m_mat->rows();}

    index_t cols() const {return m_mat->cols();}

    ///Returns the matrix
    const MatrixType& matrix() const { return *m_mat; }

private:
    const typename memory::shared<Derived>::ptr m_mat;
    bool m_symmetric;
};

/** @brief This essentially just calls the gsMatrixOp constructor, but the use of a template functions allows us to let the compiler
  * do type inference, so we don't need to type out the matrix type explicitly.
  * 
  * @note If a matrix is provided, only a reference is stored. Make sure that the matrix is not deleted too early or provide a shared
  * pointer.
  */
template <class Derived>
typename gsMatrixOp<Derived>::Ptr makeMatrixOp(const Eigen::EigenBase<Derived>& mat, bool sym=false)
{
    return memory::make_shared(new gsMatrixOp<Derived>(mat, sym));
}

/** @brief This essentially just calls the gsMatrixOp constructor, but the use of a template functions allows us to let the compiler
  * do type inference, so we don't need to type out the matrix type explicitly.
  */
template <class Derived>
typename gsMatrixOp<Derived>::Ptr makeMatrixOp(const memory::shared< Eigen::EigenBase<Derived> >& mat, bool sym=false)
{
    return memory::make_shared(new gsMatrixOp<Derived>(mat, sym));
}



/** @brief Simple adapter class to use an Eigen solver (having a compute() and a solve() method) as a linear operator.
  *
  * \ingroup Solver
  */
template <class SolverType>
class gsSolverOp : public gsLinearOperator<typename SolverType::Scalar>
{
public:
    typedef typename SolverType::Scalar T;

    typedef typename SolverType::MatrixType MatrixType;
    
    /// Shared pointer for gsSolverOp
    typedef typename memory::shared<gsSolverOp>::ptr Ptr;

    /// Unique pointer for gsSolverOp   
    typedef typename memory::unique<gsSolverOp>::ptr uPtr;
    
    
    /// Constructor taking a matrix
    gsSolverOp(const MatrixType& mat)
    {
        GISMO_ASSERT(mat.rows() == mat.cols(), "Need square matrix");
        m_size = mat.rows();

        m_solver.compute(mat);
    }

    /// Constructor taking a shared pointer
    gsSolverOp(const typename memory::shared<MatrixType>::ptr& mat)
    {
        GISMO_ASSERT(mat->rows() == mat->cols(), "Need square matrix");
        m_size = mat->rows();

        m_solver.compute(*mat);
    }
    
    /// Make function taking a matrix OR a shared pointer
    static Ptr make(const MatrixType& mat) { return memory::make_shared( new gsSolverOp(mat) ); }    
    
    void apply(const gsMatrix<T> & input, gsMatrix<T> & x) const
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
    return memory::make_shared( new gsSolverOp< Eigen::PartialPivLU< Eigen::Matrix<T, _Rows, _Cols, _Opt> > >(mat) );
}

/// @brief Convenience function to create an LU solver with partial pivoting (for dense matrices) as a gsLinearOperator taking a shared pointer.
template <class T, int _Rows, int _Cols, int _Opt>
typename gsSolverOp< Eigen::PartialPivLU< Eigen::Matrix<T, _Rows, _Cols, _Opt> > >::Ptr  makePartialPivLUSolver(const typename memory::shared< gsMatrix<T, _Rows, _Cols, _Opt> >::ptr & mat)
{
    return memory::make_shared( new gsSolverOp< Eigen::PartialPivLU< Eigen::Matrix<T, _Rows, _Cols, _Opt> > >(mat) );
}


/// @brief Convenience function to create an LU solver with full pivoting (for dense matrices) as a gsLinearOperator.
template <class T, int _Rows, int _Cols, int _Opt>
typename gsSolverOp< Eigen::FullPivLU< Eigen::Matrix<T, _Rows, _Cols, _Opt> > >::Ptr  makeFullPivLUSolver(const gsMatrix<T, _Rows, _Cols, _Opt> & mat)
{
    return memory::make_shared( new gsSolverOp< Eigen::FullPivLU< Eigen::Matrix<T, _Rows, _Cols, _Opt> > >(mat) );
}

/// @brief Convenience function to create an LU solver with full pivoting (for dense matrices) as a gsLinearOperator taking a shared pointer.
template <class T, int _Rows, int _Cols, int _Opt>
typename gsSolverOp< Eigen::FullPivLU< Eigen::Matrix<T, _Rows, _Cols, _Opt> > >::Ptr  makeFullPivLUSolver(const typename memory::shared< gsMatrix<T, _Rows, _Cols, _Opt> >::ptr & mat)
{
    return memory::make_shared( new gsSolverOp< Eigen::FullPivLU< Eigen::Matrix<T, _Rows, _Cols, _Opt> > >(mat) );
}


/// @brief Convenience function to create a Cholesky (LDL^T) solver (for dense matrices) as a gsLinearOperator.
///
/// @note Works only on symmetric (stored in lower half) and positive (semi-)definite matrices.
template <class T, int _Rows, int _Cols, int _Opt>
typename gsSolverOp< Eigen::LDLT< Eigen::Matrix<T, _Rows, _Cols, _Opt> > >::Ptr  makeCholeskySolver(const gsMatrix<T, _Rows, _Cols, _Opt> & mat)
{
    return memory::make_shared( new gsSolverOp< Eigen::LDLT< Eigen::Matrix<T, _Rows, _Cols, _Opt> > >(mat) );
}

/// @brief Convenience function to create a Cholesky (LDL^T) solver (for dense matrices) as a gsLinearOperator taking a shared pointer.
///
/// @note Works only on symmetric (stored in lower half) and positive (semi-)definite matrices.
template <class T, int _Rows, int _Cols, int _Opt>
typename gsSolverOp< Eigen::LDLT< Eigen::Matrix<T, _Rows, _Cols, _Opt> > >::Ptr  makeCholeskySolver(const typename memory::shared< gsMatrix<T, _Rows, _Cols, _Opt> >::ptr & mat)
{
    return memory::make_shared( new gsSolverOp< Eigen::LDLT< Eigen::Matrix<T, _Rows, _Cols, _Opt> > >(mat) );
}


/// @brief Convenience function to create a sparse LU solver as a gsLinearOperator.
///
/// @note This uses the default COLAMD column ordering.
template <typename T, int _Opt, typename _Index>
typename gsSolverOp< typename gsSparseSolver<T>::LU >::Ptr  makeSparseLUSolver(const gsSparseMatrix<T,_Opt,_Index> & mat)
{
    return memory::make_shared( new gsSolverOp< typename gsSparseSolver<T>::LU >(mat) );
}

/// @brief Convenience function to create a sparse LU solver as a gsLinearOperator taking a shared pointer.
///
/// @note This uses the default COLAMD column ordering.
template <typename T, int _Opt, typename _Index>
typename gsSolverOp< typename gsSparseSolver<T>::LU >::Ptr  makeSparseLUSolver(const typename memory::shared< gsSparseMatrix<T,_Opt,_Index> >::ptr & mat)
{
    return memory::make_shared( new gsSolverOp< typename gsSparseSolver<T>::LU >(mat) );
}


/// @brief Convenience function to create a sparse Cholesky (simplicial LDL^T) solver as a gsLinearOperator.
///
/// @note Works only on sparse, symmetric (stored in lower half) and positive definite matrices.
template <typename T, int _Opt, typename _Index>
typename gsSolverOp< typename gsSparseSolver<T>::SimplicialLDLT >::Ptr  makeSparseCholeskySolver(const gsSparseMatrix<T,_Opt,_Index> & mat)
{
    return memory::make_shared( new gsSolverOp<typename  gsSparseSolver<T>::SimplicialLDLT >(mat) );
}

/// @brief Convenience function to create a sparse Cholesky (simplicial LDL^T) solver as a gsLinearOperator.
///
/// @note Works only on sparse, symmetric (stored in lower half) and positive definite matrices taking a shared pointer.
template <typename T, int _Opt, typename _Index>
typename gsSolverOp< typename gsSparseSolver<T>::SimplicialLDLT >::Ptr  makeSparseCholeskySolver(const typename memory::shared< gsSparseMatrix<T,_Opt,_Index> >::ptr & mat)
{
    return memory::make_shared( new gsSolverOp<typename  gsSparseSolver<T>::SimplicialLDLT >(mat) );
}


} // namespace gismo
