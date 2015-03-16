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

void gsConjugateGradient::initIteration(const gsConjugateGradient::VectorType& rhs, gsConjugateGradient::VectorType& x0, const gsLinearOperator& precond)
{
    GISMO_ASSERT(rhs.cols()== 1, "Implemented only for single column right hand side matrix");

    int n = m_mat.cols();
    int m = 1; // == rhs.cols();
    z.resize(n,m);
    tmp.resize(n,m);
    tmp2.resize(n,m);
    p.resize(n,m);

    m_mat.apply(x0,tmp2);  //apply the system matrix
    residual = rhs - tmp2; //initial residual

    precond.apply(residual, p);      //initial search direction

    absNew = Eigen::numext::real(residual.col(0).dot(p.col(0)));  // the square of the absolute value of r scaled by invM
    rhsNorm2 = rhs.squaredNorm();
    residualNorm2 = 0;
    threshold = m_tol*m_tol*rhsNorm2;
    m_numIter = 0;

    if(m_calcEigenvals)
    {
        delta.clear();
        delta.resize(1,0);
        delta.reserve(m_maxIters);

        gamma.clear();
        gamma.reserve(m_maxIters);
    }
}


bool gsConjugateGradient::step( gsConjugateGradient::VectorType& x, const gsLinearOperator& precond )
{
    m_mat.apply(p,tmp); //apply system matrix

    real_t alpha = absNew / p.col(0).dot(tmp.col(0));   // the amount we travel on dir
    if(m_calcEigenvals)
        delta.back()+=(1./alpha);

    x += alpha * p;                       // update solution
    residual -= alpha * tmp;              // update residual

    residualNorm2 = residual.squaredNorm();
    if(residualNorm2 < threshold)
        return true;

    precond.apply(residual, z);          // approximately solve for "A z = residual"

    real_t absOld = absNew;

    absNew = Eigen::numext::real(residual.col(0).dot(z.col(0)));     // update the absolute value of r
    real_t beta = absNew / absOld;            // calculate the Gram-Schmidt value used to create the new search direction
    p = z + beta * p;                             // update search direction

    if(m_calcEigenvals)
    {
        gamma.push_back(-std::sqrt(beta)/alpha);
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
    //there is probably a better option...
   Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigensolver(L.toDense());

   eigs = eigensolver.eigenvalues();
}


}


