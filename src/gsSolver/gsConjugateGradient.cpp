/** @file gsConjugateGradient.cpp

    @brief Conjugate gradient solver

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): C. Hofreither
*/
#include <gsSolver/gsConjugateGradient.h>
#include <gsSolver/gsSolverUtils.h>

namespace gismo
{

bool gsConjugateGradient::initIteration( const gsConjugateGradient::VectorType& rhs, gsConjugateGradient::VectorType& x )
{
    if (Base::initIteration(rhs,x))
        return true;
    
    int n = m_mat->cols();
    int m = 1;                                                          // == rhs.cols();
    m_tmp.resize(n,m);
    m_update.resize(n,m);

    m_mat->apply(x,m_tmp);                                              // apply the system matrix
    m_res = rhs - m_tmp;                                                // initial residual

    m_precond->apply(m_res,m_update);                                   // initial search direction

    m_abs_new = m_res.col(0).dot(m_update.col(0)); // the square of the absolute value of r scaled by invM

    if (m_calcEigenvals)
    {
        delta.clear();
        delta.resize(1,0);
        delta.reserve(m_max_iters);

        gamma.clear();
        gamma.reserve(m_max_iters);

        m_eigsAreCalculated = true;
    }
    return false;
}


bool gsConjugateGradient::step( gsConjugateGradient::VectorType& x )
{
    m_mat->apply(m_update,m_tmp);                                      // apply system matrix

    real_t alpha = m_abs_new / m_update.col(0).dot(m_tmp.col(0));      // the amount we travel on dir
    if (m_calcEigenvals)
        delta.back()+=(1./alpha);

    x += alpha * m_update;                                             // update solution
    m_res -= alpha * m_tmp;                                            // update residual

    m_error = m_res.norm() / m_rhs_norm;
    if (m_error < m_tol)
        return true;

    m_precond->apply(m_res, m_tmp);                                    // approximately solve for "A tmp = residual"

    real_t abs_old = m_abs_new;

    m_abs_new = m_res.col(0).dot(m_tmp.col(0));   // update the absolute value of r
    real_t beta = m_abs_new / abs_old;                                 // calculate the Gram-Schmidt value used to create the new search direction
    m_update = m_tmp + beta * m_update;                                // update search direction

    if (m_calcEigenvals)
    {
        gamma.push_back(-math::sqrt(beta)/alpha);
        delta.push_back(beta/alpha);
    }
    return false;
}

real_t gsConjugateGradient::getConditionNumber()
{
    GISMO_ASSERT(m_eigsAreCalculated,"No data for eigenvalues was collected, call setCalcEigenvalues(true) and solve with an arbitrary right hand side");
    gsLanczosMatrix<real_t> L(gamma,delta);

    return L.maxEigenvalue()/L.minEigenvalue();
}

void gsConjugateGradient::getEigenvalues(gsMatrix<real_t>& eigs )
{
   GISMO_ASSERT(m_eigsAreCalculated,"No data for eigenvalues was collected, call setCalcEigenvalues(true) and solve with an arbitrary right hand side");

   gsLanczosMatrix<real_t> LM(gamma,delta);
   gsSparseMatrix<real_t> L;
   LM.matrixForm(L);
   // there is probably a better option...
   Eigen::SelfAdjointEigenSolver<Eigen::Matrix<real_t,Dynamic,Dynamic> > eigensolver(L.toDense());

   eigs = eigensolver.eigenvalues();
}


}


