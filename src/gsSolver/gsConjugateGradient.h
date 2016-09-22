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
 */

class GISMO_EXPORT gsConjugateGradient : public gsIterativeSolver
{
public:
    typedef gsMatrix<real_t>    VectorType;

    /// Contructor. See gsIterativeSolver for details.
    template< typename OperatorType >
    gsConjugateGradient( const OperatorType& mat, const gsLinearOperator<>::Ptr& precond, index_t max_iters=1000, real_t tol=1e-10, bool calcEigenval=false )
        : gsIterativeSolver(mat, precond, max_iters, tol), m_calcEigenvals(calcEigenval), m_eigsAreCalculated(false) {}

    /// Contructor. See gsIterativeSolver for details.
    template< typename OperatorType >
    gsConjugateGradient( const OperatorType& mat, index_t max_iters=1000, real_t tol=1e-10, bool calcEigenval=false )
        : gsIterativeSolver(mat, max_iters, tol), m_calcEigenvals(calcEigenval), m_eigsAreCalculated(false) {}

    bool initIteration( const VectorType& rhs, VectorType& x );
    bool step( VectorType& x );

    /// @brief specify if you want to store data for eigenvalue estimation
    /// @param flag true stores the coefficients of the lancos matrix, false not.
    void setCalcEigenvalues( bool flag )     { m_calcEigenvals = flag ;}

    /// @brief returns the condition number of the (preconditioned) system matrix
    real_t getConditionNumber();

    /// @brief returns the eigenvalues of the Lanczos matrix
    void getEigenvalues( gsMatrix<real_t>& eigs );

private:
    using gsIterativeSolver::m_mat;
    using gsIterativeSolver::m_precond;
    using gsIterativeSolver::m_max_iters;
    using gsIterativeSolver::m_tol;
    using gsIterativeSolver::m_num_iter;
    using gsIterativeSolver::m_initial_error;
    using gsIterativeSolver::m_error;


    VectorType z, tmp, tmp2, p;
    VectorType residual;
    real_t absNew;

    bool m_calcEigenvals;
    bool m_eigsAreCalculated;

    std::vector<real_t> delta, gamma;
};

} // namespace gismo

