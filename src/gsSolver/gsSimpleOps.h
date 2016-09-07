/** @file gsSimpleOps.h

    @brief Colloction of some simple preconditioners.

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
    
/// Update \a x with a Richardson sweep
GISMO_EXPORT void dampedRichardsonSweep(const Eigen::SparseMatrix<real_t>& A, gsMatrix<real_t>& x, const gsMatrix<real_t>& f, real_t tau = (real_t)(1.));

/// Update \a x with a Jacobi sweep
GISMO_EXPORT void JacobiSweep(const Eigen::SparseMatrix<real_t>& A, gsMatrix<real_t>& x, const gsMatrix<real_t>& f);

/// Update \a x with a damped Jacobi sweep
GISMO_EXPORT void dampedJacobiSweep(const Eigen::SparseMatrix<real_t>& A, gsMatrix<real_t>& x, const gsMatrix<real_t>& f, real_t tau = (real_t)(0.5));

/// Update \a x with a forward Gauss-Seidel sweep
GISMO_EXPORT void gaussSeidelSweep(const Eigen::SparseMatrix<real_t>& A, gsMatrix<real_t>& x, const gsMatrix<real_t>& f);

/// Update \a x with a backward Gauss-Seidel sweep
GISMO_EXPORT void reverseGaussSeidelSweep(const Eigen::SparseMatrix<real_t>& A, gsMatrix<real_t>& x, const gsMatrix<real_t>& f);

/// Preforms a block Gauss-Seidel on the degrees of freedom in DoFs.
GISMO_EXPORT void gaussSeidelSingleBlock(const Eigen::SparseMatrix<real_t>& A, gsMatrix<real_t>& x, const gsMatrix<real_t>& f, gsVector<index_t>& DoFs);


/// @brief Richardson preconditioner
///
template <typename MatrixType, int UpLo = Eigen::Lower>
class gsRichardsonOp : public gsLinearOperator
{
public:

    /// Shared pointer for gsRichardsonOp
    typedef memory::shared_ptr< gsRichardsonOp > Ptr;

    /// Unique pointer for gsRichardsonOp   
    typedef typename memory::unique< gsRichardsonOp >::ptr uPtr;    
    
    /// @brief Contructor with given matrix
    gsRichardsonOp(const MatrixType& _mat, real_t tau = 1.)
        : m_mat(_mat), m_numOfSweeps(1), m_tau(tau) {}
        
    static Ptr make(const MatrixType& _mat, real_t tau = 1.) { return shared( new gsRichardsonOp(_mat,tau) ); }

    void apply(const gsMatrix<real_t> & input, gsMatrix<real_t> & x) const
    {
        assert( m_mat.rows() == input.rows() && m_mat.cols() == m_mat.rows() && input.cols() == 1);

        // For the first sweep, we do not need to multiply with the matrix
        x = input;
        x *= m_tau;
        
        for (index_t k = 1; k < m_numOfSweeps; ++k)
        {
            gsMatrix<real_t> temp = input - m_mat * x;
            x += m_tau * temp;
        }
    }

    index_t rows() const {return m_mat.rows();}
    index_t cols() const {return m_mat.cols();}

    /// Set number of sweeps.
    void setNumOfSweeps(index_t n) {
        GISMO_ASSERT ( n > 0, "Number of sweeps needs to be positive. ");
        m_numOfSweeps=n;
    }

    ///Returns the matrix
    MatrixType matrix() const { return m_mat; }

private:
    MatrixType m_mat;
    index_t m_numOfSweeps;
    real_t m_tau;
};

/// @brief Jacobi preconditioner
///
/// Requires a positive definite matrix.
template <typename MatrixType, int UpLo = Eigen::Lower>
class gsJacobiOp : public gsLinearOperator
{
public:

    /// Shared pointer for gsJacobiOp
    typedef memory::shared_ptr< gsJacobiOp > Ptr;

    /// Unique pointer for gsJacobiOp   
    typedef typename memory::unique< gsJacobiOp >::ptr uPtr;    

    /// @brief Contructor with given matrix
    gsJacobiOp(const MatrixType& _mat, real_t tau = 1.)
        : m_mat(_mat), m_numOfSweeps(1), m_tau(tau) {}
        
    static Ptr make(const MatrixType& _mat, real_t tau = 1.) { return shared( new gsJacobiOp(_mat,tau) ); }

    void apply(const gsMatrix<real_t> & input, gsMatrix<real_t> & x) const
    {
        assert( m_mat.rows() == input.rows() && m_mat.cols() == m_mat.rows() && input.cols() == 1);

        // For the first sweep, we do not need to multiply with the matrix
        x = input;
        x.array() /= m_mat.diagonal().array();
        x *= m_tau;
        
        for (index_t k = 1; k < m_numOfSweeps; ++k)
        {
            gsMatrix<real_t> temp = input - m_mat * x;
            temp.array() /= m_mat.diagonal().array();
            x += m_tau * temp;
        }
    }

    index_t rows() const {return m_mat.rows();}
    index_t cols() const {return m_mat.cols();}

    /// Set number of sweeps.
    void setNumOfSweeps(index_t n) {
        GISMO_ASSERT ( n > 0, "Number of sweeps needs to be positive. ");
        m_numOfSweeps=n;
    }

    ///Returns the matrix
    MatrixType matrix() const { return m_mat; }

private:
    MatrixType m_mat;
    index_t m_numOfSweeps;
    real_t m_tau;
};

/// @brief Gauss-Seidel preconditioner
///
/// Requires a positive definite matrix.
template <typename MatrixType, int UpLo = Eigen::Lower>
class gsGaussSeidelOp : public gsLinearOperator
{
public:

    /// Shared pointer for gsGaussSeidelOp
    typedef memory::shared_ptr< gsGaussSeidelOp > Ptr;

    /// Unique pointer for gsGaussSeidelOp   
    typedef typename memory::unique< gsGaussSeidelOp >::ptr uPtr;   
    
    /// @brief Contructor with given matrix
    gsGaussSeidelOp(const MatrixType& _mat)
        : m_mat(_mat), m_numOfSweeps(1) {}
        
    static Ptr make(const MatrixType& _mat) { return shared( new gsGaussSeidelOp(_mat) ); }

    void apply(const gsMatrix<real_t> & input, gsMatrix<real_t> & x) const
    {
        x.setZero(rows(), input.cols());

        for (index_t k = 0; k < m_numOfSweeps; ++k)
        {
            gaussSeidelSweep(m_mat,x,input);
        }
    }

    index_t rows() const {return m_mat.rows();}
    index_t cols() const {return m_mat.cols();}

    /// Set number of sweeps of to symmetric Gauss-Seidel perform (default is 1).
    void setNumOfSweeps(index_t n) {
        GISMO_ASSERT ( n > 0, "Number of sweeps needs to be positive. ");
        m_numOfSweeps=n;
    }

    ///Returns the matrix
    MatrixType matrix() const { return m_mat; }

private:
    MatrixType m_mat;
    index_t m_numOfSweeps;
};


/// @brief Symmetric Gauss-Seidel preconditioner
///
/// Requires a positive definite matrix. Does first
/// one forward Gauss-Seidel sweep then one backward
/// Gauss-Seidel sweep.
template <typename MatrixType, int UpLo = Eigen::Lower>
class gsSymmetricGaussSeidelOp : public gsLinearOperator
{
public:

    /// Shared pointer for gsSymmetricGaussSeidelOp
    typedef memory::shared_ptr< gsSymmetricGaussSeidelOp > Ptr;

    /// Unique pointer for gsSymmetricGaussSeidelOp   
    typedef typename memory::unique< gsSymmetricGaussSeidelOp >::ptr uPtr; 
    
    /// @brief Contructor with given matrix
    gsSymmetricGaussSeidelOp(const MatrixType& _mat, index_t numOfSweeps = 1)
        : m_mat(_mat), m_numOfSweeps(numOfSweeps) {}
        
    static Ptr make(const MatrixType& _mat, index_t numOfSweeps = 1) { return shared( new gsSymmetricGaussSeidelOp(_mat,numOfSweeps) ); }

    void apply(const gsMatrix<real_t> & input, gsMatrix<real_t> & x) const
    {
        x.setZero(rows(), input.cols());

        for (index_t k = 0; k < m_numOfSweeps; ++k)
        {
            gaussSeidelSweep(m_mat,x,input);
            //x.array() *= m_mat.diagonal().array();
            reverseGaussSeidelSweep(m_mat,x,input);
        }
    }

    index_t rows() const {return m_mat.rows();}

    index_t cols() const {return m_mat.cols();}

    /// Set number of sweeps of to symmetric Gauss-Seidel perform (default is 1).
    void setNumOfSweeps(index_t n)    { m_numOfSweeps= n; }

    ///Returns the matrix
    MatrixType matrix() const { return m_mat; }

private:
    MatrixType m_mat;
    index_t m_numOfSweeps;
};


} // namespace gismo
