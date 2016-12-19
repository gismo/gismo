/** @file gsSimpleOps.h

    @brief Collection of some simple preconditioners.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): J. Sogn, C. Hofreither, A. Mantzaflaris
*/
#pragma once

#include <gsCore/gsLinearAlgebra.h>
#include <gsSolver/gsLinearOperator.h>

namespace gismo
{
    
/// @brief Update \a x with a Richardson sweep
/// \ingroup Solver
GISMO_EXPORT void dampedRichardsonSweep(const gsSparseMatrix<real_t>& A, gsMatrix<real_t>& x, const gsMatrix<real_t>& f, real_t tau = (real_t)(1.));

/// @brief Update \a x with a Jacobi sweep
/// \ingroup Solver
GISMO_EXPORT void JacobiSweep(const gsSparseMatrix<real_t>& A, gsMatrix<real_t>& x, const gsMatrix<real_t>& f);

/// @brief Update \a x with a damped Jacobi sweep
/// \ingroup Solver
GISMO_EXPORT void dampedJacobiSweep(const gsSparseMatrix<real_t>& A, gsMatrix<real_t>& x, const gsMatrix<real_t>& f, real_t tau = (real_t)(0.5));

/// @brief Update \a x with a forward Gauss-Seidel sweep
/// \ingroup Solver
GISMO_EXPORT void gaussSeidelSweep(const gsSparseMatrix<real_t>& A, gsMatrix<real_t>& x, const gsMatrix<real_t>& f);

/// @brief Update \a x with a backward Gauss-Seidel sweep
/// \ingroup Solver
GISMO_EXPORT void reverseGaussSeidelSweep(const gsSparseMatrix<real_t>& A, gsMatrix<real_t>& x, const gsMatrix<real_t>& f);

/// @brief Preforms a block Gauss-Seidel on the degrees of freedom in DoFs.
/// \inrgoup Solver
GISMO_EXPORT void gaussSeidelSingleBlock(const gsSparseMatrix<real_t>& A, gsMatrix<real_t>& x, const gsMatrix<real_t>& f, gsVector<index_t>& DoFs);


/// @brief Richardson preconditioner
///
/// \ingroup Solver
template <typename MatrixType>
class gsRichardsonOp : public gsSteppableOperator<typename MatrixType::Scalar>
{
    typedef typename memory::shared_ptr<MatrixType> MatrixPtr;
    typedef typename MatrixType::Nested              NestedMatrix;

public:
    typedef typename MatrixType::Scalar T;
    
    /// Shared pointer for gsRichardsonOp
    typedef typename memory::shared_ptr< gsRichardsonOp > Ptr;

    /// Unique pointer for gsRichardsonOp   
    typedef memory::unique_ptr< gsRichardsonOp > uPtr;    

    /// @brief Constructor with given matrix
    explicit gsRichardsonOp(const MatrixType& _mat)
    : m_mat(), m_expr(_mat.derived()), m_tau(1) {}
    
    /// @brief Constructor with shared pointer to matrix
    explicit gsRichardsonOp(const MatrixPtr& _mat)
    : m_mat(_mat), m_expr(m_mat->derived()), m_tau(1) { }
    
    static Ptr make(const MatrixType& _mat) 
    { return memory::make_shared( new gsRichardsonOp(_mat) ); }

    static Ptr make(const MatrixPtr& _mat) 
    { return memory::make_shared( new gsRichardsonOp(_mat) ); }

    void step(const gsMatrix<T> & rhs, gsMatrix<T> & x) const
    {
        GISMO_ASSERT( m_expr.rows() == rhs.rows() && m_expr.cols() == m_expr.rows(),
                      "Dimensions do not match.");

        x += m_tau * ( rhs - m_expr * x );
    }

    // We use our own apply implementation as we can save one multiplication. This is important if the number
    // of sweeps is 1: Then we can save *all* multiplications.
    void apply(const gsMatrix<T> & input, gsMatrix<T> & x) const
    {
        GISMO_ASSERT( m_expr.rows() == input.rows() && m_expr.cols() == m_expr.rows(),
                      "Dimensions do not match.");

        // For the first sweep, we do not need to multiply with the matrix
        x.noalias() = m_tau * input;
        
        for (index_t k = 1; k < m_num_of_sweeps; ++k)
            x += m_tau * ( input - m_expr * x );
    }

    index_t rows() const {return m_expr.rows();}
    index_t cols() const {return m_expr.cols();}

    /// Set scaling parameter
    void setScaling(const T tau) { m_tau = tau; }

    ///Returns the matrix
    const MatrixType & matrix() const { return m_mat; }

private:
    const MatrixPtr m_mat;  ///< Shared pointer to matrix (if needed)
    NestedMatrix    m_expr; ///< Nested Eigen expression

    using gsLinearOperator<T>::m_num_of_sweeps;
    T m_tau;
};

/**
   \brief Returns a shared pointer to a Richardson operator referring on \a mat
*/
template <class Derived>
typename gsRichardsonOp<Derived>::Ptr makeRichardsonOp(const Eigen::EigenBase<Derived>& mat)
{ return gsRichardsonOp<Derived>::make(mat.derived()); }

/// @brief Jacobi preconditioner
///
/// Requires a positive definite matrix.
///
/// \ingroup Solver
template <typename MatrixType>
class gsJacobiOp : public gsSteppableOperator<typename MatrixType::Scalar>
{
    typedef typename memory::shared_ptr<MatrixType> MatrixPtr;
    typedef typename MatrixType::Nested              NestedMatrix;

public:
    typedef typename MatrixType::Scalar T;
    
    /// Shared pointer for gsJacobiOp
    typedef typename memory::shared_ptr< gsJacobiOp > Ptr;

    /// Unique pointer for gsJacobiOp   
    typedef memory::unique_ptr< gsJacobiOp > uPtr;    

    /// @brief Constructor with given matrix
    explicit gsJacobiOp(const MatrixType& _mat)
    : m_mat(), m_expr(_mat.derived()), m_tau(1) {}
    
    /// @brief Constructor with shared pointer to matrix
    explicit gsJacobiOp(const MatrixPtr& _mat)
    : m_mat(_mat), m_expr(m_mat->derived()), m_tau(1) { }
    
    static Ptr make(const MatrixType& _mat) 
    { return memory::make_shared( new gsJacobiOp(_mat) ); }

    static Ptr make(const MatrixPtr& _mat) 
    { return memory::make_shared( new gsJacobiOp(_mat) ); }
  
    void step(const gsMatrix<T> & rhs, gsMatrix<T> & x) const
    {
        GISMO_ASSERT( m_expr.rows() == rhs.rows() && m_expr.cols() == m_expr.rows(),
                      "Dimensions do not match.");
        
        GISMO_ASSERT( rhs.cols() == 1, "This operator is only implemented for a single right-hand side." );

        x.array() += m_tau * ( rhs - m_expr * x ).array() / m_expr.diagonal().array();
    }
  
    // We use our own apply implementation as we can save one multiplication. This is important if the number
    // of sweeps is 1: Then we can save *all* multiplications.
    void apply(const gsMatrix<T> & input, gsMatrix<T> & x) const
    {
        GISMO_ASSERT( m_expr.rows() == input.rows() && m_expr.cols() == m_expr.rows(),
                      "Dimensions do not match.");
        
        GISMO_ASSERT( input.cols() == 1, "This operator is only implemented for a single right-hand side." );

        // For the first sweep, we do not need to multiply with the matrix
        x.array() = m_tau * input.array() / m_expr.diagonal().array();
        
        for (index_t k = 1; k < m_num_of_sweeps; ++k)
            x.array() += m_tau * ( input - m_expr * x ).array() / m_expr.diagonal().array();
    }

    index_t rows() const {return m_expr.rows();}
    index_t cols() const {return m_expr.cols();}

    /// Set scaling parameter
    void setScaling(const T tau) { m_tau = tau; }

    ///Returns the matrix
    NestedMatrix matrix() const { return m_expr; }

private:
    const MatrixPtr m_mat;  ///< Shared pointer to matrix (if needed)
    NestedMatrix    m_expr; ///< Nested Eigen expression
    using gsLinearOperator<T>::m_num_of_sweeps;
    T m_tau;
};

/**
   \brief Returns a shared pointer to a Jacobi operator referring on \a mat
*/
template <class Derived>
typename gsJacobiOp<Derived>::Ptr makeJacobiOp(const Eigen::EigenBase<Derived>& mat)
{ return gsJacobiOp<Derived>::make(mat.derived()); }


/// @brief Gauss-Seidel preconditioner
///
/// Requires a positive definite matrix.
///
/// \ingroup Solver
template <typename MatrixType>
class gsGaussSeidelOp : public gsSteppableOperator<typename MatrixType::Scalar>
{
    typedef typename memory::shared_ptr<MatrixType> MatrixPtr;
    typedef typename MatrixType::Nested              NestedMatrix;

public:
    typedef typename MatrixType::Scalar T;

    /// Shared pointer for gsGaussSeidelOp
    typedef typename memory::shared_ptr< gsGaussSeidelOp > Ptr;

    /// Unique pointer for gsGaussSeidelOp   
    typedef memory::unique_ptr< gsGaussSeidelOp > uPtr;   
    
    /// @brief Constructor with given matrix
    explicit gsGaussSeidelOp(const MatrixType& _mat)
    : m_mat(), m_expr(_mat.derived()) {}
    
    /// @brief Constructor with shared pointer to matrix
    explicit gsGaussSeidelOp(const MatrixPtr& _mat)
    : m_mat(_mat), m_expr(m_mat->derived()) { }
    
    static Ptr make(const MatrixType& _mat) 
    { return memory::make_shared( new gsGaussSeidelOp(_mat) ); }

    static Ptr make(const MatrixPtr& _mat) 
    { return memory::make_shared( new gsGaussSeidelOp(_mat) ); }

    void step(const gsMatrix<T> & rhs, gsMatrix<T> & x) const
    {
        gaussSeidelSweep(m_expr,x,rhs);
    }

    index_t rows() const {return m_expr.rows();}
    index_t cols() const {return m_expr.cols();}

    ///Returns the matrix
    const MatrixType & matrix() const { return m_mat; }
    
private:
    const MatrixPtr m_mat;  ///< Shared pointer to matrix (if needed)
    NestedMatrix    m_expr; ///< Nested Eigen expression
};

/**
   \brief Returns a shared pointer to a Gauss-Seidel operator referring on \a mat
*/
template <class Derived>
typename gsGaussSeidelOp<Derived>::Ptr makeGaussSeidelOp(const Eigen::EigenBase<Derived>& mat)
{ return gsGaussSeidelOp<Derived>::make(mat.derived()); }

/// @brief Symmetric Gauss-Seidel preconditioner
///
/// Requires a positive definite matrix. Does first
/// one forward Gauss-Seidel sweep then one backward
/// Gauss-Seidel sweep.
///
/// \ingroup Solver
template <typename MatrixType>
class gsSymmetricGaussSeidelOp : public gsSteppableOperator<typename MatrixType::Scalar>
{
    typedef typename memory::shared_ptr<MatrixType> MatrixPtr;
    typedef typename MatrixType::Nested          NestedMatrix;

public:
    typedef typename MatrixType::Scalar T;
    
    /// Shared pointer for gsSymmetricGaussSeidelOp
    typedef typename memory::shared_ptr< gsSymmetricGaussSeidelOp > Ptr;

    /// Unique pointer for gsSymmetricGaussSeidelOp   
    typedef memory::unique_ptr< gsSymmetricGaussSeidelOp > uPtr; 
    
    /// @brief Constructor with given matrix
    explicit gsSymmetricGaussSeidelOp(const MatrixType& _mat)
    : m_mat(), m_expr(_mat.derived()) {}

    /// @brief Constructor with shared pointer to matrix
    explicit gsSymmetricGaussSeidelOp(const MatrixPtr& _mat)
    : m_mat(_mat), m_expr(m_mat->derived()) {}
    
    static Ptr make(const MatrixType& _mat) 
    { return memory::make_shared( new gsSymmetricGaussSeidelOp(_mat) ); }

    static Ptr make(const MatrixPtr& _mat) 
    { return memory::make_shared( new gsSymmetricGaussSeidelOp(_mat) ); }

    void step(const gsMatrix<T> & rhs, gsMatrix<T> & x) const
    {
        gaussSeidelSweep(m_expr,x,rhs);
        //x.array() *= m_expr.diagonal().array();
        reverseGaussSeidelSweep(m_expr,x,rhs);
    }

    index_t rows() const {return m_expr.rows();}

    index_t cols() const {return m_expr.cols();}


    ///Returns the matrix
    NestedMatrix matrix() const { return m_expr; }

private:
    const MatrixPtr m_mat;  ///< Shared pointer to matrix (if needed)
    NestedMatrix    m_expr; ///< Nested Eigen expression

};

/**
   \brief Returns a shared pointer to a Symmetric Gauss-Seidel operator referring on \a mat
*/
template <class Derived>
typename gsSymmetricGaussSeidelOp<Derived>::Ptr
makeSymmetricGaussSeidelOp(const Eigen::EigenBase<Derived>& mat)
{ return gsSymmetricGaussSeidelOp<Derived>::make(mat.derived()); }


} // namespace gismo
