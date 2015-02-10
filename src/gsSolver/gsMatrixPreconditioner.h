/** @file gsMatrixPreconditioner.h

    @brief Simple class to give a matrix preconditioner properties

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): J. Sogn
*/
#pragma once

#include <gsCore/gsLinearAlgebra.h>
#include <gsSolver/gsLinearOperator.h>

namespace gismo
{

/// @brief Simple adapter class to use a matrix (or matrix-like object) as a linear operator. Needed for the iterative method classes.
///
/// \ingroup Solver

template <typename MatrixType>
class gsMatrixPreconditioner : public gsLinearOperator
{
public:

    gsMatrixPreconditioner(const MatrixType& _matPre, bool sym=false)
        : matPre(_matPre), m_symmetric(sym)
    {
    }

    /// Destructor
    ~gsMatrixPreconditioner()
    {
    }

    void apply(const gsMatrix<real_t> & input, gsMatrix<real_t> & x) const
    {
        if (m_symmetric)
            x.noalias() = matPre.template selfadjointView<Lower>() * input;
        else
            x.noalias() = matPre * input;
    }

    index_t rows() const {return matPre.rows();}

    index_t cols() const {return matPre.cols();}

    ///Returns the matrix
    const MatrixType& matrix() const { return matPre; }

private:
    const MatrixType& matPre;
    bool m_symmetric;
};

/** This essentially just calls the gsMatrixPreconditioner constructor, but the
 * use of a template functions allows us to let the compiler do type inference,
 * so we don't need to type out the matrix type explicitly.
 */
template <typename MatrixType>
gsMatrixPreconditioner<MatrixType>* makeMatrixOperator(const MatrixType& mat, bool sym=false)
{
    return new gsMatrixPreconditioner<MatrixType>(mat, sym);
}

} // namespace gismo
