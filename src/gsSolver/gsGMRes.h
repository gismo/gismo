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

    /// Contructor. See gsIterativeSolver to find out what you can pass for mat.
    template< typename OperatorType >
    gsGMRes( const OperatorType& mat, index_t max_iters=1000, real_t tol=1e-10 )
        : gsIterativeSolver(mat, max_iters, tol) {}

    bool initIteration( const VectorType& rhs, VectorType& x, const gsLinearOperator<>& precond );
    bool step( VectorType& x, const gsLinearOperator<>& precond );
    void finalizeIteration( const VectorType& rhs, VectorType& x );

private:

    /// Solves the Upper triangular system Ry = gg
    /// and stores the solution in the private member y.
    void solveUpperTriangular(const gsMatrix<real_t> & R, const gsMatrix<real_t> & gg)
    {
       y = R.triangularView<Eigen::Upper>().solve(gg);
    }

private:
    using gsIterativeSolver::m_mat;
    using gsIterativeSolver::m_max_iters;
    using gsIterativeSolver::m_tol;
    using gsIterativeSolver::m_num_iter;
    using gsIterativeSolver::m_initial_error;
    using gsIterativeSolver::m_error;


    gsMatrix<real_t> xInit, tmp, g, g_tmp, h_tmp, y, w;
    gsMatrix<real_t> m_rhs, residual;
    gsMatrix<real_t> H_prew, H, Omega, Omega_prew, Omega_tmp, Omega_prew_tmp;
    std::vector<gsMatrix<real_t> > v;
    real_t beta;
};

} // namespace gismo
