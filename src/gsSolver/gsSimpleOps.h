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
template<typename T>
GISMO_EXPORT void dampedRichardsonSweep(const gsSparseMatrix<T>& A, gsMatrix<T>& x, const gsMatrix<T>& f, T tau = (T)(1.));

/// @brief Update \a x with a Jacobi sweep
/// \ingroup Solver
template<typename T>
GISMO_EXPORT void jacobiSweep(const gsSparseMatrix<T>& A, gsMatrix<T>& x, const gsMatrix<T>& f);

/// @brief Update \a x with a damped Jacobi sweep
/// \ingroup Solver
template<typename T>
GISMO_EXPORT void dampedJacobiSweep(const gsSparseMatrix<T>& A, gsMatrix<T>& x, const gsMatrix<T>& f, T tau = (T)(0.5));

/// @brief Update \a x with a forward Gauss-Seidel sweep
/// \ingroup Solver
template<typename T>
GISMO_EXPORT void gaussSeidelSweep(const gsSparseMatrix<T>& A, gsMatrix<T>& x, const gsMatrix<T>& f);

/// @brief Update \a x with a backward Gauss-Seidel sweep
/// \ingroup Solver
template<typename T>
GISMO_EXPORT void reverseGaussSeidelSweep(const gsSparseMatrix<T>& A, gsMatrix<T>& x, const gsMatrix<T>& f);

/// @brief Preforms a block Gauss-Seidel on the degrees of freedom in DoFs.
/// \ingroup Solver
template<typename T>
GISMO_EXPORT void gaussSeidelSingleBlock(const gsSparseMatrix<T>& A, gsMatrix<T>& x, const gsMatrix<T>& f, gsVector<index_t>& DoFs);


/// @brief Richardson preconditioner
///
/// \ingroup Solver
template <typename MatrixType>
class gsRichardsonOp : public gsSteppableOperator<typename MatrixType::Scalar>
{
    typedef typename memory::shared_ptr<MatrixType>  MatrixPtr;
    typedef typename MatrixType::Nested              NestedMatrix;

public:
    /// Scalar type
    typedef typename MatrixType::Scalar T;

    /// Shared pointer for gsRichardsonOp
    typedef memory::shared_ptr< gsRichardsonOp > Ptr;

    /// Unique pointer for gsRichardsonOp
    typedef memory::unique_ptr< gsRichardsonOp > uPtr;
 
    /// Base class
    typedef gsSteppableOperator<T> Base;

    /// @brief Constructor with given matrix
    explicit gsRichardsonOp(const MatrixType& _mat)
    : m_mat(), m_expr(_mat.derived()), m_tau(1) {}

    /// @brief Constructor with shared pointer to matrix
    explicit gsRichardsonOp(const MatrixPtr& _mat)
    : m_mat(_mat), m_expr(m_mat->derived()), m_tau(1) { }

    static uPtr make(const MatrixType& _mat) 
    { return memory::make_unique( new gsRichardsonOp(_mat) ); }

    static uPtr make(const MatrixPtr& _mat) 
    { return memory::make_unique( new gsRichardsonOp(_mat) ); }

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
    void setScaling(const T tau) { m_tau = tau;  }

    /// Get scaling parameter
    void getScaling()            { return m_tau; }

    /// Get the default options as gsOptionList object
    static gsOptionList defaultOptions()
    {
        gsOptionList opt = Base::defaultOptions();
        opt.addReal( "Scaling", "Scaling parameter of the Richardson iteration", 1 );
        return opt;
    }

    /// Set options based on a gsOptionList object
    virtual void setOptions(const gsOptionList & opt)
    {
        Base::setOptions(opt);
        m_tau = opt.askReal( "Scaling", m_tau );
    }

    /// Returns the matrix
    const MatrixType & matrix() const { return m_mat; }

private:
    const MatrixPtr m_mat;  ///< Shared pointer to matrix (if needed)
    NestedMatrix    m_expr; ///< Nested Eigen expression

    using Base::m_num_of_sweeps;
    T m_tau;
};

/**
   \brief Returns a smart pointer to a Richardson operator referring on \a mat
*/
template <class Derived>
typename gsRichardsonOp<Derived>::uPtr makeRichardsonOp(const Eigen::EigenBase<Derived>& mat)
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
    typedef typename MatrixType::Nested             NestedMatrix;

public:
    /// Scalar type
    typedef typename MatrixType::Scalar T;

    /// Shared pointer for gsJacobiOp
    typedef memory::shared_ptr< gsJacobiOp > Ptr;

    /// Unique pointer for gsJacobiOp   
    typedef memory::unique_ptr< gsJacobiOp > uPtr;

    /// Base class
    typedef gsSteppableOperator<T> Base;

    /// @brief Constructor with given matrix
    explicit gsJacobiOp(const MatrixType& _mat)
    : m_mat(), m_expr(_mat.derived()), m_tau(1) {}

    /// @brief Constructor with shared pointer to matrix
    explicit gsJacobiOp(const MatrixPtr& _mat)
    : m_mat(_mat), m_expr(m_mat->derived()), m_tau(1) { }

    static uPtr make(const MatrixType& _mat) 
    { return memory::make_unique( new gsJacobiOp(_mat) ); }

    static uPtr make(const MatrixPtr& _mat) 
    { return memory::make_unique( new gsJacobiOp(_mat) ); }
 
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
    void setScaling(const T tau) { m_tau = tau;  }

    /// Get scaling parameter
    void getScaling()            { return m_tau; }

    /// Get the default options as gsOptionList object
    static gsOptionList defaultOptions()
    {
        gsOptionList opt = Base::defaultOptions();
        opt.addReal( "Scaling", "Scaling parameter of the Jacobi iteration", 1 );
        return opt;
    }

    /// Set options based on a gsOptionList object
    virtual void setOptions(const gsOptionList & opt)
    {
        Base::setOptions(opt);
        m_tau = opt.askReal( "Scaling", m_tau );
    }

    /// Returns the matrix
    NestedMatrix matrix() const { return m_expr; }

private:
    const MatrixPtr m_mat;  ///< Shared pointer to matrix (if needed)
    NestedMatrix    m_expr; ///< Nested Eigen expression
    using Base::m_num_of_sweeps;
    T m_tau;
};


/**
   \brief Returns a smart pointer to a Jacobi operator referring on \a mat
*/
template <class Derived>
typename gsJacobiOp<Derived>::uPtr makeJacobiOp(const Eigen::EigenBase<Derived>& mat)
{ return gsJacobiOp<Derived>::make(mat.derived()); }


namespace gsGaussSeidel
{
    enum ordering
    {
        forward = 0,
        reverse = 1,
        symmetric = 2
    };
}

/// @brief Gauss-Seidel preconditioner
///
/// Requires a positive definite matrix.
///
/// \ingroup Solver
template <typename MatrixType, gsGaussSeidel::ordering ordering = gsGaussSeidel::forward>
class gsGaussSeidelOp : public gsSteppableOperator<typename MatrixType::Scalar>
{
    typedef typename memory::shared_ptr<MatrixType> MatrixPtr;
    typedef typename MatrixType::Nested             NestedMatrix;

public:
    /// Scalar type
    typedef typename MatrixType::Scalar T;

    /// Shared pointer for gsGaussSeidelOp
    typedef memory::shared_ptr< gsGaussSeidelOp > Ptr;

    /// Unique pointer for gsGaussSeidelOp
    typedef memory::unique_ptr< gsGaussSeidelOp > uPtr;

    /// Base class
    typedef gsSteppableOperator<T> Base;

    /// @brief Constructor with given matrix
    explicit gsGaussSeidelOp(const MatrixType& _mat)
    : m_mat(), m_expr(_mat.derived()) {}

    /// @brief Constructor with shared pointer to matrix
    explicit gsGaussSeidelOp(const MatrixPtr& _mat)
    : m_mat(_mat), m_expr(m_mat->derived()) { }

    static uPtr make(const MatrixType& _mat) 
    { return memory::make_unique( new gsGaussSeidelOp(_mat) ); }

    static uPtr make(const MatrixPtr& _mat) 
    { return memory::make_unique( new gsGaussSeidelOp(_mat) ); }

    void step(const gsMatrix<T> & rhs, gsMatrix<T> & x) const
    {
        if (ordering == gsGaussSeidel::forward )
            gaussSeidelSweep<T>(m_expr,x,rhs);
        if (ordering == gsGaussSeidel::reverse )
            reverseGaussSeidelSweep<T>(m_expr,x,rhs);
        if (ordering == gsGaussSeidel::symmetric )
        {
            gaussSeidelSweep<T>(m_expr,x,rhs);
            reverseGaussSeidelSweep<T>(m_expr,x,rhs);
        }
    }

    index_t rows() const {return m_expr.rows();}
    index_t cols() const {return m_expr.cols();}

    /// Returns the matrix
    const MatrixType & matrix() const { return m_mat; }

private:
    const MatrixPtr m_mat;  ///< Shared pointer to matrix (if needed)
    NestedMatrix    m_expr; ///< Nested Eigen expression
};

/**
   \brief Returns a smart pointer to a Gauss-Seidel operator referring on \a mat
*/
template <class Derived>
typename gsGaussSeidelOp<Derived>::uPtr makeGaussSeidelOp(const Eigen::EigenBase<Derived>& mat)
{ return gsGaussSeidelOp<Derived>::make(mat.derived()); }

/**
   \brief Returns a smart pointer to a reverse Gauss-Seidel operator referring on \a mat
*/
template <class Derived>
typename gsGaussSeidelOp<Derived,gsGaussSeidel::reverse>::uPtr makeReverseGaussSeidelOp(const Eigen::EigenBase<Derived>& mat)
{ return gsGaussSeidelOp<Derived,gsGaussSeidel::reverse>::make(mat.derived()); }

/**
   \brief Returns a smart pointer to a symmetric Gauss-Seidel operator referring on \a mat
*/
template <class Derived>
typename gsGaussSeidelOp<Derived,gsGaussSeidel::symmetric>::uPtr makeSymmetricGaussSeidelOp(const Eigen::EigenBase<Derived>& mat)
{ return gsGaussSeidelOp<Derived,gsGaussSeidel::symmetric>::make(mat.derived()); }

} // namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsSimpleOps.hpp)
#endif
