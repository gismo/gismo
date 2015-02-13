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
    typedef gsMatrix<real_t>                VectorType;

    /// Constructor for general linear operator
    gsMinimalResidual(const gsLinearOperator& _mat, index_t _maxIt=1000, real_t _tol=1e-10)
        : gsIterativeSolver(_mat, _maxIt, _tol) {}

    /// Constructor for sparse matrix
    template<typename T, int _Options, typename _Index>
    gsMinimalResidual(const gsSparseMatrix<T, _Options, _Index > & _mat, index_t _maxIt=1000, real_t _tol=1e-10)
        : gsIterativeSolver(makeMatrixOperator(_mat, true), _maxIt, _tol) {}

    /// Constructor for dense matrix
    template<class T, int _Rows, int _Cols, int _Options>
    gsMinimalResidual(const gsMatrix<T, _Rows, _Cols, _Options> & _mat, index_t _maxIt=1000, real_t _tol=1e-10)
        : gsIterativeSolver(makeMatrixOperator(_mat, true), _maxIt, _tol) {}

    void initIteration( const VectorType& rhs, const VectorType& x0, const gsLinearOperator& precond);

    void solve(const VectorType& rhs, VectorType& x, const gsLinearOperator& precond)
        {
            initIteration(rhs, x, precond);

            while(m_numIter < m_maxIters)
            {
                if (step(x, precond))
                    break;
                m_numIter++;
            }
            m_error = std::sqrt(residualNorm2 / rhsNorm2);
        }

    bool step( VectorType& x, const gsLinearOperator& precond );


private:
    using gsIterativeSolver::m_mat;
    using gsIterativeSolver::m_error;
    using gsIterativeSolver::m_maxIters;
    using gsIterativeSolver::m_numIter;
    using gsIterativeSolver::m_tol;

    gsMatrix<real_t> vPrew, v, vNew, wPrew, w, wNew,zNew, z,xPrew, m_rhs, residual, tmp, tmp2;
    real_t residualNorm2, threshold, rhsNorm2;
    real_t eta,gammaPrew,gamma,gammaNew,sPrew,s,sNew,cPrew,c,cNew;
};

} // namespace gismo
