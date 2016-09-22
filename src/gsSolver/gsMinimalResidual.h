/** @file gsMinimalResidual.h

    @brief Preconditioned iterative solver using the minimal residual method.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): J. Sogn
*/
#pragma once

#include <gsSolver/gsIterativeSolver.h>

namespace gismo
{

class GISMO_EXPORT gsMinimalResidual : public gsIterativeSolver
{
public:
    typedef gsMatrix<real_t>    VectorType;

    /// Constructor for general linear operator
    gsMinimalResidual( const gsLinearOperator<>::Ptr& mat, index_t max_iter=1000, real_t tol=1e-10 )
        : gsIterativeSolver(mat, max_iter, tol) {}

    /// Constructor for sparse matrix
    ///
    /// @note: This does not copy the matrix. So, make sure that the matrix is not deleted.
    template<int _Options, typename _Index>
    gsMinimalResidual( const gsSparseMatrix<real_t, _Options, _Index > & mat, index_t max_iter=1000, real_t tol=1e-10 )
        : gsIterativeSolver(makeMatrixOp(mat, true), max_iter, tol) {}

    /// Constructor for dense matrix
    ///
    /// @note: This does not copy the matrix. So, make sure that the matrix is not deleted.
    template<int _Rows, int _Cols, int _Options>
    gsMinimalResidual( const gsMatrix<real_t, _Rows, _Cols, _Options> & mat, index_t max_iter=1000, real_t tol=1e-10 )
        : gsIterativeSolver(makeMatrixOp(mat, true), max_iter, tol) {}

    bool initIteration( const VectorType& rhs, VectorType& x, const gsLinearOperator<>& precond );
    bool step( VectorType& x, const gsLinearOperator<>& precond );


private:
    using gsIterativeSolver::m_mat;
    using gsIterativeSolver::m_max_iters;
    using gsIterativeSolver::m_tol;
    using gsIterativeSolver::m_num_iter;
    using gsIterativeSolver::m_initial_error;
    using gsIterativeSolver::m_error;

    
    gsMatrix<real_t> vPrew, v, vNew, wPrew, w, wNew,zNew, z, xPrew, m_rhs, residual, tmp, tmp2;
    real_t eta, gammaPrew, gamma, gammaNew, sPrew, s, sNew, cPrew, c, cNew;
};

} // namespace gismo

