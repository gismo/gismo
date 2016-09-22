/** @file gsGMRes.h

    @brief Preconditioned iterative solver using the generalized minimal residual method.

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

class GISMO_EXPORT gsGMRes: public gsIterativeSolver
{
public:
    typedef gsMatrix<real_t>    VectorType;

    ///Contructor for general linear operator
    gsGMRes( const gsLinearOperator<>::Ptr& mat, index_t max_iter=1000, real_t tol=1e-10 )
        : gsIterativeSolver(mat, max_iter, tol) {}

    ///Contructor for sparse matrix
    ///
    /// @note: This does not copy the matrix. So, make sure that the matrix is not deleted.
    template<int _Options, typename _Index>
    gsGMRes( const gsSparseMatrix<real_t, _Options, _Index > & mat, index_t max_iter=1000, real_t tol=1e-10 )
        : gsIterativeSolver(makeMatrixOp(mat, true), max_iter, tol) {}

    ///Contructor for dense matrix
    ///
    /// @note: This does not copy the matrix. So, make sure that the matrix is not deleted.
    template<int _Rows, int _Cols, int _Options>
    gsGMRes( const gsMatrix<real_t, _Rows, _Cols, _Options> & mat, index_t max_iter=1000, real_t tol=1e-10 )
        : gsIterativeSolver(makeMatrixOp(mat, true), max_iter, tol) {}

    void initIteration( const VectorType& rhs, const VectorType& x0, const gsLinearOperator<>& precond );

    void solve( const VectorType& rhs, VectorType& x, const gsLinearOperator<>& precond );

    /// Solve system without preconditioner
    void solve( const VectorType& rhs, VectorType& x )
    {
        gsIdentityOp<> preConId(m_mat->rows());
        solve(rhs, x, preConId);
    }

    bool step( VectorType& x, const gsLinearOperator<>& precond );

private:

    /// Solves the Upper triangular system Ry = gg
    /// and stores the solution in the private member y.
    void solveUpperTriangular(const gsMatrix<real_t> & R, const gsMatrix<real_t> & gg)
    {
       y = R.triangularView<Eigen::Upper>().solve(gg);
    }

private:
    using gsIterativeSolver::m_mat;
    using gsIterativeSolver::m_error;
    using gsIterativeSolver::m_max_iters;
    using gsIterativeSolver::m_num_iter;
    using gsIterativeSolver::m_tol;

    gsMatrix<real_t> xInit, tmp, g, g_tmp, h_tmp, y, w;
    gsMatrix<real_t> m_rhs, residual;
    gsMatrix<real_t> H_prew, H, Omega, Omega_prew, Omega_tmp, Omega_prew_tmp;
    std::vector<gsMatrix<real_t> > v;
    real_t residualNorm2, threshold, rhsNorm2, beta;
};

} // namespace gismo
