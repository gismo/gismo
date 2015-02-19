/** @file gsSimplePreconditioners.h

    @brief Colloction of some simple preconditioners.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): J. Sogn
*/
#pragma once

#include <gsCore/gsLinearAlgebra.h>
#include <gsCore/gsMultiPatch.h>
#include <gsCore/gsMultiBasis.h>
#include <gsSolver/gsLinearOperator.h>
#include <gsAssembler/gsGenericAssembler.h>
#include <gsSolver/gsSmoother.h>

namespace gismo
{

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

/// @brief Identity preconditioner ("do nothing"), must be square!
class gsIdentityPreconditioner : public gsLinearOperator
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
