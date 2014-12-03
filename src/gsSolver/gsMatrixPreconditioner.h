/** @file gsMatrixPreconditioner.h

    @brief Simple class to give a matrix preconditior properties

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): J. Sogn
*/
#pragma once

#include <gsCore/gsLinearAlgebra.h>
#include <gsSolver/gsPreconditioner.h>

namespace gismo
{

///Simple class to give a matrix (or sparse matrix) preconditior properties which is needed for the interative method classes.

template <typename MatrixType, int UpLo = Eigen::Lower>
class gsMatrixPreconditioner : public gsPreconditioner
{
public:

    gsMatrixPreconditioner(const MatrixType& _matPre) : matPre(_matPre)
    {
    }

    void apply(const gsMatrix<real_t> & input, gsMatrix<real_t> & x) const
    {
        x.noalias() = matPre.template selfadjointView<UpLo>() * input;
    }

    index_t rows() const {return matPre.rows();}

    index_t cols() const {return matPre.cols();}

    ///Returns the matrix
    MatrixType matrix() const { return matPre; }

private:
    const MatrixType& matPre;
};

} // namespace gismo
