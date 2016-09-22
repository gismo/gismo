/** @file gsConjugateGradient.h

    @brief Conjugate gradient solver

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): C. Hofreither
*/

#pragma once

#include <gsSolver/gsIterativeSolver.h>

namespace gismo
{

/** The conjugate gradient implementation from Eigen, adapted to allow for more
 *  general preconditioners and better iteration control. Also capable of using
 *  a gsLinearOperator as matrix.
 *
 *
 *  Only implemented for single right hand side!
 */

class GISMO_EXPORT gsConjugateGradient : public gsIterativeSolver
{
public:
    typedef gsMatrix<real_t>    VectorType;

    /// Constructor for general linear operator
    gsConjugateGradient( const gsLinearOperator<>::Ptr& mat, int max_iter=1000, real_t tol=1e-10, bool calcEigenval=false )
        : gsIterativeSolver(mat, max_iter, tol), m_calcEigenvals(calcEigenval), m_eigsAreCalculated(false) {}

    /// Constructor for sparse matrix
    ///
    /// @note: This does not copy the matrix. So, make sure that the matrix is not deleted.
    template<int _Options, typename _Index>
    gsConjugateGradient( const gsSparseMatrix<real_t, _Options, _Index > & mat, index_t max_iter=1000, real_t tol=1e-10, bool calcEigenval=false )
        : gsIterativeSolver(makeMatrixOp(mat, true), max_iter, tol), m_calcEigenvals(calcEigenval), m_eigsAreCalculated(false)  {}

    /// Constructor for dense matrix
    ///
    /// @note: This does not copy the matrix. So, make sure that the matrix is not deleted.
    template<int _Rows, int _Cols, int _Options>
    gsConjugateGradient( const gsMatrix<real_t, _Rows, _Cols, _Options> & mat, index_t max_iter=1000, real_t tol=1e-10, bool calcEigenval=false )
        : gsIterativeSolver(makeMatrixOp(mat, true), max_iter, tol), m_calcEigenvals(calcEigenval), m_eigsAreCalculated(false)  {}

    void initIteration( const VectorType& rhs, VectorType& x0, const gsLinearOperator<>& precond );

    void solve( const VectorType& rhs, VectorType& x, const gsLinearOperator<>& precond )
        {
            initIteration(rhs, x, precond);

            while(m_num_iter < m_max_iters)
            {
                m_num_iter++;
                if (step(x, precond))
                    break;
            }
            m_error = math::sqrt(residualNorm2 / rhsNorm2);

        }

    /// Solve system without preconditioner
    void solve( const VectorType& rhs, VectorType& x )
    {
        gsIdentityOp<> preConId(m_mat->rows());
        solve(rhs, x, preConId);
    }

    bool step( VectorType& x, const gsLinearOperator<>& precond );

    /// @brief specify if you want to store data for eigenvalue estimation
    /// @param flag true stores the coefficients of the lancos matrix, false not.
    void setCalcEigenvalues( bool flag )     { m_calcEigenvals = flag ;}

    /// @brief returns the condition number of the (preconditioned) system matrix
    real_t getConditionNumber();

    /// @brief returns the eigenvalues of the Lanczos matrix
    void getEigenvalues( gsMatrix<real_t>& eigs );

private:
    using gsIterativeSolver::m_mat;
    using gsIterativeSolver::m_error;
    using gsIterativeSolver::m_max_iters;
    using gsIterativeSolver::m_num_iter;
    using gsIterativeSolver::m_tol;

    VectorType z, tmp, tmp2, p;
    VectorType residual;
    real_t absNew, residualNorm2, threshold, rhsNorm2;

    bool m_calcEigenvals;
    bool m_eigsAreCalculated;

    std::vector<real_t> delta, gamma;
};

} // namespace gismo

