/** @file gsGMRes.h

    @brief Preconditioned iterative solver using the generalized minimal residual method.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): J. Sogn
*/
#pragma once

#include <gismo.h>
#include <gsSolver/gsIterativeSolver.h>

namespace gismo
{

class GISMO_EXPORT gsGMRes: public gsIterativeSolver
{
public:
    typedef gsMatrix<real_t>                VectorType;

    ///Contructor for general linear operator
    gsGMRes(const gsLinearOperator& _mat, index_t _maxIt=1000, real_t _tol=1e-10)
        : gsIterativeSolver(_mat, _maxIt, _tol) {}

    ///Contructor for sparse matrix
    template<typename T, int _Options, typename _Index>
    gsGMRes(const gsSparseMatrix<T, _Options, _Index > & _mat, index_t _maxIt=1000, real_t _tol=1e-10)
        : gsIterativeSolver(makeMatrixOperator(_mat, true), _maxIt, _tol) {}

    ///Contructor for dense matrix
    template<class T, int _Rows, int _Cols, int _Options>
    gsGMRes(const gsMatrix<T, _Rows, _Cols, _Options> & _mat, index_t _maxIt=1000, real_t _tol=1e-10)
        : gsIterativeSolver(makeMatrixOperator(_mat, true), _maxIt, _tol) {}

    void initIteration( const VectorType& rhs, const VectorType& x0, const gsLinearOperator& precond);

    void solve(const VectorType& rhs, VectorType& x, const gsLinearOperator& precond);

    /// Solve system without preconditioner
    void solve(const VectorType& rhs, VectorType& x)
    {
        gsIdentityPreconditioner preConId(m_mat.rows());
        solve(rhs, x, preConId);
    }

    bool step( VectorType& x, const gsLinearOperator& precond );

private:

    /// Solves the Upper triangular system Ry = g
    /// and stores the solution in the private member y.
    void solveUpperTriangular(gsMatrix<> R, gsMatrix<> g)
    {
       y = R.triangularView<Eigen::Upper>().solve(g);
    }

private:
    using gsIterativeSolver::m_mat;
    using gsIterativeSolver::m_error;
    using gsIterativeSolver::m_maxIters;
    using gsIterativeSolver::m_numIter;
    using gsIterativeSolver::m_tol;

    gsMatrix<real_t> xInit, tmp, g, g_tmp, h_tmp, y, w;
    gsMatrix<real_t> m_rhs, residual;
    gsMatrix<real_t> H_prew, H, Omega, Omega_prew, Omega_tmp, Omega_prew_tmp;
    std::vector<gsMatrix<real_t> > v;
    real_t residualNorm2, threshold, rhsNorm2, beta;
};

} // namespace gismo

//////////////////////////////////////////////////
//////////////////////////////////////////////////


#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsGMRes.cpp)
#endif
