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

class gsMinimalResidual : public gsIterativeSolver
{
public:
    typedef gsMatrix<real_t>                VectorType;

    ///Contructor for general linear operator
    gsMinimalResidual(const gsLinearOperator& _mat, index_t _maxIt=1000, real_t _tol=1e-10)
        : gsIterativeSolver(_mat, _maxIt, _tol) {}

    ///Contructor for sparse matrix
    template<typename T, int _Options, typename _Index>
    gsMinimalResidual(const gsSparseMatrix<T, _Options, _Index > & _mat, index_t _maxIt=1000, real_t _tol=1e-10)
        : gsIterativeSolver(_mat, _maxIt, _tol) {}

    ///Contructor for dense matrix
    template<class T, int _Rows, int _Cols, int _Options>
    gsMinimalResidual(const gsMatrix<T, _Rows, _Cols, _Options> & _mat, index_t _maxIt=1000, real_t _tol=1e-10)
        : gsIterativeSolver(_mat, _maxIt, _tol) {}

    void initIteration( const VectorType& rhs, const VectorType& x0, const gsLinearOperator& precond)
    {
        GISMO_ASSERT(rhs.cols()== 1, "Implemented only for single columns right hand side matrix");

        int n = m_mat.cols();
        int m = 1;//rhs.cols();
        m_rhs = rhs;
        xPrew = x0;
        vPrew.setZero(n,m); vNew.setZero(n,m);
        wPrew.setZero(n,m); w.setZero(n,m); wNew.setZero(n,m);
        tmp2.setZero(n,1);

        m_mat.apply(x0,tmp2);
        v = m_rhs - tmp2;

        precond.apply(v, z);

        gammaPrew = 1; gamma = math::sqrt(z.col(0).dot(v.col(0))); gammaNew = 1;
        eta = gamma;
        sPrew = 0; s = 0; sNew = 0;
        cPrew = 1; c = 1; cNew = 1;
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
            z /= gamma;
            m_mat.apply(z,tmp);

            real_t delta = z.col(0).dot(tmp.col(0));
            vNew = tmp - (delta/gamma)*v - (gamma/gammaPrew)*vPrew;
            precond.apply(vNew, zNew);
            gammaNew = math::sqrt(zNew.col(0).dot(vNew.col(0)));
            real_t a0 = c*delta - cPrew*s*gamma;
            real_t a1 = math::sqrt(a0*a0 + gammaNew*gammaNew);
            real_t a2 = s*delta + cPrew*c*gamma;
            real_t a3 = sPrew*gamma;
            cNew = a0/a1;
            sNew = gammaNew/a1;
            wNew = (z - a3*wPrew - a2*w)/a1;
            x = xPrew + cNew*eta*wNew;
            eta = -sNew*eta;

            //Test for convergence
            m_mat.apply(x,tmp2);
            residual = m_rhs - tmp2;
            residualNorm2 = residual.squaredNorm();
            if(residualNorm2 < threshold)
                return true;

            //Update variables
            vPrew = v; v = vNew;
            wPrew = w; w = wNew;
            z = zNew;
            xPrew = x;
            gammaPrew = gamma; gamma = gammaNew;
            sPrew = s; s = sNew;
            cPrew = c; c = cNew;
            return false;
        }


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
