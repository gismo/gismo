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
  
template <class MatrixType>
class gsMatrixOp : public gsLinearOperator<typename MatrixType::Scalar>
{
public:
    typedef typename MatrixType::Scalar T;
    
    /// Shared pointer for gsMatrixOp
    typedef memory::shared_ptr<gsMatrixOp> Ptr;

    /// Unique pointer for gsMatrixOp   
    typedef typename memory::unique<gsMatrixOp>::ptr uPtr;
    
    /// Shared pointer to the matrix type
    typedef memory::shared_ptr<MatrixType> MatrixPtr;
    
    /// @brief Constructor taking a reference
    /// @note This does not copy the matrix. Make sure that the matrix is not deleted too early or provide a shared pointer.
    gsMatrixOp(const MatrixType& mat, bool sym=false)
    : m_mat(shared_not_owned(&mat)), m_symmetric(sym)
    {}
    
    /// Constructor taking a shared pointer
    gsMatrixOp(const MatrixPtr& mat, bool sym=false)
        : m_mat(mat), m_symmetric(sym) {}

    /// @brief Make function returning a smart pointer
    /// @note This does not copy the matrix. Make sure that the matrix is not deleted too early or provide a shared pointer.    
    static Ptr make(const MatrixType& mat, bool sym=false) { return shared( new gsMatrixOp(mat,sym) ); }

    /// Make function returning a smart pointer
    static Ptr make(const MatrixPtr& mat, bool sym=false) { return shared( new gsMatrixOp(mat,sym) ); }

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
    const MatrixPtr m_mat;
    bool m_symmetric;
};

/** @brief This essentially just calls the gsMatrixOp constructor, but the use of a template functions allows us to let the compiler
  * do type inference, so we don't need to type out the matrix type explicitly.
  * 
  * @note If a matrix is provided, only a reference is stored. Make sure that the matrix is not deleted too early or provide a shared
  * pointer.
  */
template <class MatrixType>
typename gsMatrixOp<MatrixType>::Ptr makeMatrixOp(const MatrixType& mat, bool sym=false)
{
    return memory::make_shared(new gsMatrixOp<MatrixType>(mat, sym));
}

/** @brief This essentially just calls the gsMatrixOp constructor, but the use of a template functions allows us to let the compiler
  * do type inference, so we don't need to type out the matrix type explicitly.
  */
template <class MatrixType>
typename gsMatrixOp<MatrixType>::Ptr makeMatrixOp(const memory::shared_ptr<MatrixType>& mat, bool sym=false)
{
    return memory::make_shared(new gsMatrixOp<MatrixType>(mat, sym));
}



/**
  * @brief Simple adapter class to use the transpose of a matrix as a linear operator.
  * This should, of course, be done without transposing the matrix itself.
  *
  * \ingroup Solver
  */
  
template <class MatrixType>
class gsTransposedMatrixOp : public gsLinearOperator<typename MatrixType::Scalar>
{
public:
    typedef typename MatrixType::Scalar T;
    
    /// Shared pointer for gsTransposedMatrixOp
    typedef memory::shared_ptr<gsTransposedMatrixOp> Ptr;

    /// Unique pointer for gsTransposedMatrixOp   
    typedef typename memory::unique<gsTransposedMatrixOp>::ptr uPtr;
    
    /// Shared pointer to the matrix type
    typedef memory::shared_ptr<MatrixType> MatrixPtr;
    
    /// @brief Constructor taking a reference
    /// @note This does not copy the matrix. Make sure that the matrix is not deleted too early or provide a shared pointer.
    gsTransposedMatrixOp(const MatrixType& mat)
    : m_mat(shared_not_owned(&mat))
    {}
    
    /// Constructor taking a shared pointer
    gsTransposedMatrixOp(const MatrixPtr& mat)
        : m_mat(mat) {}

    /// @brief Make function returning a smart pointer
    /// @note This does not copy the matrix. Make sure that the matrix is not deleted too early or provide a shared pointer.    
    static Ptr make(const MatrixType& mat) { return memory::make_shared( new gsTransposedMatrixOp(mat) ); }

    /// Make function returning a smart pointer
    static Ptr make(const MatrixPtr& mat) { return memory::make_shared( new gsTransposedMatrixOp(mat) ); }


    void apply(const gsMatrix<T> & input, gsMatrix<T> & x) const
    {
        x.noalias() = (*m_mat).transpose() * input;
    }

    index_t rows() const {return m_mat->cols();}

    index_t cols() const {return m_mat->rows();}

private:
    const MatrixPtr m_mat;
};


/** @brief This essentially just calls the gsTransposedMatrixOp constructor, but the use of a template functions allows us to let the
  * compiler do type inference, so we don't need to type out the matrix type explicitly.
  * 
  * @note If a matrix is provided, only a reference is stored. Make sure that the matrix is not deleted too early or provide a shared
  * pointer.
  */
template <class MatrixType>
typename gsTransposedMatrixOp<MatrixType>::Ptr makeTransposedMatrixOp(const MatrixType& mat)
{
    return memory::make_shared(new gsTransposedMatrixOp<MatrixType>(mat));
}

/** @brief This essentially just calls the gsTransposedMatrixOp constructor, but the use of a template functions allows us to let the
  * compiler do type inference, so we don't need to type out the matrix type explicitly.
  * 
  */
template <class MatrixType>
typename gsTransposedMatrixOp<MatrixType>::Ptr makeTransposedMatrixOp(const memory::shared_ptr<MatrixType>& mat)
{
    return memory::make_shared(new gsTransposedMatrixOp<MatrixType>(mat));
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
    
    /// Shared pointer for gsSolverOp
    typedef memory::shared_ptr<gsSolverOp> Ptr;

    /// Unique pointer for gsSolverOp   
    typedef typename memory::unique<gsSolverOp>::ptr uPtr;
    
    
    /// Constructor taking a matrix
    template <class MatrixType>
    gsSolverOp(const MatrixType& mat)
    {
        GISMO_ASSERT(mat.rows() == mat.cols(), "Need square matrix");
        m_size = mat.rows();

        m_solver.compute(mat);
    }

    /// Constructor taking a shared pointer
    template <class MatrixType>
    gsSolverOp(const memory::shared_ptr<MatrixType>& mat)
    {
        GISMO_ASSERT(mat->rows() == mat->cols(), "Need square matrix");
        m_size = mat->rows();

        m_solver.compute(*mat);
    }
    
    /// Make function taking a matrix OR a shared pointer
    template <class MatrixType>
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
typename gsSolverOp< Eigen::PartialPivLU< Eigen::Matrix<T, _Rows, _Cols, _Opt> > >::Ptr  makePartialPivLUSolver(const memory::shared_ptr< gsMatrix<T, _Rows, _Cols, _Opt> > & mat)
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
typename gsSolverOp< Eigen::FullPivLU< Eigen::Matrix<T, _Rows, _Cols, _Opt> > >::Ptr  makeFullPivLUSolver(const memory::shared_ptr< gsMatrix<T, _Rows, _Cols, _Opt> > & mat)
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
typename gsSolverOp< Eigen::LDLT< Eigen::Matrix<T, _Rows, _Cols, _Opt> > >::Ptr  makeCholeskySolver(const memory::shared_ptr< gsMatrix<T, _Rows, _Cols, _Opt> > & mat)
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
typename gsSolverOp< typename gsSparseSolver<T>::LU >::Ptr  makeSparseLUSolver(const memory::shared_ptr< gsSparseMatrix<T,_Opt,_Index> > & mat)
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
typename gsSolverOp< typename gsSparseSolver<T>::SimplicialLDLT >::Ptr  makeSparseCholeskySolver(const memory::shared_ptr< gsSparseMatrix<T,_Opt,_Index> > & mat)
{
    return memory::make_shared( new gsSolverOp<typename  gsSparseSolver<T>::SimplicialLDLT >(mat) );
}


} // namespace gismo
