/** @file gsSimplePreconditioners.h

    @brief Colloction of some simple preconditioners.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): J. Sogn, C. Hofreither
*/
#pragma once

#include <gsCore/gsLinearAlgebra.h>
#include <gsCore/gsMultiPatch.h>
#include <gsCore/gsMultiBasis.h>
#include <gsSolver/gsLinearOperator.h>
#include <gsAssembler/gsGenericAssembler.h>

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
class gsRichardsonPreconditioner : public gsLinearOperator
{
public:

    /// @brief Contructor with given matrix
    gsRichardsonPreconditioner(const MatrixType& _mat, real_t tau = 1.)
        : m_mat(_mat), m_numOfSweeps(1), m_tau(tau) {}

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
class gsJacobiPreconditioner : public gsLinearOperator
{
public:

    /// @brief Contructor with given matrix
    gsJacobiPreconditioner(const MatrixType& _mat, real_t tau = 1.)
        : m_mat(_mat), m_numOfSweeps(1), m_tau(tau) {}

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
class gsGaussSeidelPreconditioner : public gsLinearOperator
{
public:

    /// @brief Contructor with given matrix
    gsGaussSeidelPreconditioner(const MatrixType& _mat)
        : m_mat(_mat), m_numOfSweeps(1) {}

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
class gsSymmetricGaussSeidelPreconditioner : public gsLinearOperator
{
public:

    /// @brief Contructor with given matrix
    gsSymmetricGaussSeidelPreconditioner(const MatrixType& _mat, index_t numOfSweeps = 1)
        : m_mat(_mat), m_numOfSweeps(numOfSweeps) {}

    /// @brief Contructor with build the mass matrix from \a patches and \a basis
    //TODO: is this really what a "simple" preconditioner should do?
    gsSymmetricGaussSeidelPreconditioner(const gsMultiPatch<real_t> patches, gsMultiBasis<real_t> basis, index_t numOfSweeps = 1)
        : m_numOfSweeps(numOfSweeps)
    {
        //Assemble the mass matrix for the pressure space
        gsGenericAssembler<real_t> massConst(patches, basis);
        const gsSparseMatrix<> & massMatrixBtmp = massConst.assembleMass();

        //Get full matrix (not just lower triangular)
        massMatrixBtmp.cols();
        m_mat = massConst.fullMatrix();
    }

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
