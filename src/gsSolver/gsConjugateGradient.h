/** @file gsConjugateGradient.hpp

    @brief Conjugate gradient implementation from Eigen

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

class gsConjugateGradient : public gsIterativeSolver
{
public:
    typedef gsMatrix<real_t>    VectorType;

    /// Constructor for general linear operator
    gsConjugateGradient(const gsLinearOperator& _mat, int _maxIt=1000, real_t _tol=1e-10)
        : gsIterativeSolver(_mat, _maxIt, _tol){}

    /// Constructor for sparse matrix
    template<typename T, int _Options, typename _Index>
    gsConjugateGradient(const gsSparseMatrix<T, _Options, _Index > & _mat, index_t _maxIt=1000, real_t _tol=1e-10)
        : gsIterativeSolver(makeMatrixOperator(_mat.template selfadjointView<Lower>()), _maxIt, _tol) {}

    /// Constructor for dense matrix
    template<class T, int _Rows, int _Cols, int _Options>
    gsConjugateGradient(const gsMatrix<T, _Rows, _Cols, _Options> & _mat, index_t _maxIt=1000, real_t _tol=1e-10)
        : gsIterativeSolver(makeMatrixOperator(_mat.template selfadjointView<Lower>()), _maxIt, _tol) {}

    void initIteration(const VectorType& rhs, VectorType& x0, const gsLinearOperator& precond)
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

    }

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

    bool step( VectorType& x, const gsLinearOperator& precond )
        {

            m_mat.apply(p,tmp); //apply system matrix

            real_t alpha = absNew / p.col(0).dot(tmp.col(0));   // the amount we travel on dir
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
            return false;
        }

private:
    using gsIterativeSolver::m_mat;
    using gsIterativeSolver::m_error;
    using gsIterativeSolver::m_maxIters;
    using gsIterativeSolver::m_numIter;
    using gsIterativeSolver::m_tol;

    VectorType z, tmp, tmp2, p;
    VectorType residual;
    real_t absNew, residualNorm2, threshold, rhsNorm2;
};

} // namespace gismo
