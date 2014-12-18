/** @file gsPreconditioner.h

    @brief Simple abstract class for (preconditioner) operators.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): J. Sogn
*/
#pragma once

#include <gsCore/gsLinearAlgebra.h>
#include <gsSolver/gsSmoother.h>

namespace gismo
{

/// @brief Simple abstract class for (preconditioner) operators.
///
/// Simple abstract class for (preconditioner) operators.
/// The derived classes have to contain the functions: apply(), cols(), and rows().
class gsPreconditioner
{
public:

    gsPreconditioner()
    {
    }

    /**
     * @brief apply the operator on the input vector and store the result in x
     * @param input Input vector
     * @param x     result vector
     */
    virtual void apply(const gsMatrix<real_t> & input, gsMatrix<real_t> & x) const = 0;

    ///Returns the number of rows in the preconditioner
    virtual index_t rows() const = 0;

    ///Returns the number of columns in the preconditioner
    virtual index_t cols() const = 0;

}; // gsPreconditioner



/// @brief Symmetric Gauss-Seidel preconditioner
///
/// Requires a positive definite matrix. Does first
/// one forward Gauss-Seidel sweep then one backward
/// Gauss-Seidel sweep.
template <typename MatrixType, int UpLo = Eigen::Lower>
class gsSymmetricGaussSeidelPreconditioner : public gsPreconditioner
{
public:

    gsSymmetricGaussSeidelPreconditioner(const MatrixType& _mat, index_t numOfSweeps = 1)
        : m_mat(_mat), m_numOfSweeps(numOfSweeps) {}

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
    const MatrixType& m_mat; 
    index_t m_numOfSweeps;
};

/// @brief Identity preconditioner ("do nothing"), must be square!
class gsIdentityPreconditioner : public gsPreconditioner
{
public:

    gsIdentityPreconditioner(index_t dim) : m_dim(dim) {}


    void apply(const gsMatrix<real_t> & input, gsMatrix<real_t> & x) const
    {
        x = input;
    }

    index_t rows() const {return m_dim;}

    index_t cols() const {return m_dim;}

private:
    const index_t m_dim;
};


} // namespace gismo
