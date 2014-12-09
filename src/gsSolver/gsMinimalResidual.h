/** @file gsMinimalResidual.h

    @brief Preconditioned interative solver using the minimal residual method.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): J. Sogn
*/
#pragma once

#include <gsCore/gsLinearAlgebra.h>
#include <gsSolver/gsMatrixPreconditioner.h>

namespace gismo
{
//in case of MatrixType is a child of gsPreconditioner, it must have the typedefs Scalar and RealScalar


// if MatrixType is a gsPreconditioner, then do its apply method
void applyMatrix(const gsPreconditioner& mat,const gsMatrix<real_t>& x, gsMatrix<real_t>& res)
{
    mat.apply(x,res);
}

//general case for matrix multiplication
template< typename MatrixType, int UpLo>
void applyMatrix(const MatrixType& mat, const gsMatrix<real_t>& x, gsMatrix<real_t>& res)
{
     res.noalias() = mat.template selfadjointView<UpLo>() * x;
}


template <typename MatrixType, int UpLo = Eigen::Lower>
class gsMinimalResidual
{
public:
    typedef typename MatrixType::Scalar Scalar;
    typedef typename MatrixType::RealScalar RealScalar;

    gsMinimalResidual(const MatrixType& _mat, int _maxIt=1000, RealScalar _tol=1e-10)
        : mat(_mat), tol(_tol), maxIters(_maxIt), numIter(0)
    {
    }

    template<typename Rhs, typename Dest, typename Preconditioner>
    void initIteration( const Rhs& rhs, const Dest& x0, const Preconditioner& precond)
    {
        GISMO_ASSERT(rhs.cols()== 1, "Implemented only for single collum right hand side matrix");

        int n = mat.cols();
        int m = 1;//rhs.cols();
        m_rhs = rhs;
        xPrew = x0;
        vPrew.setZero(n,m); vNew.setZero(n,m);
        wPrew.setZero(n,m); w.setZero(n,m); wNew.setZero(n,m);
        tmp2.setZero(n,1);

        applyMatrix(mat,x0,tmp2);
        v = m_rhs - tmp2;

        precond.apply(v, z);

        gammaPrew = 1; gamma = math::sqrt(z.col(0).dot(v.col(0))); gammaNew = 1;
        eta = gamma;
        sPrew = 0; s = 0; sNew = 0;
        cPrew = 1; c = 1; cNew = 1;
        rhsNorm2 = rhs.squaredNorm();
        residualNorm2 = 0;
        threshold = tol*tol*rhsNorm2;
        numIter = 0;
    }

    template<typename Rhs, typename Dest, typename Preconditioner>
    void solve(const Rhs& rhs, Dest& x, const Preconditioner& precond)
        {
            initIteration(rhs, x, precond);

            while(numIter < maxIters)
            {
                if (step(x, precond))
                    break;
                numIter++;
            }
            m_error = std::sqrt(residualNorm2 / rhsNorm2);
        }

    template<typename Dest, typename Preconditioner>
    bool step( Dest& x, const Preconditioner& precond )
        {
            z /= gamma;
            applyMatrix(mat,z,tmp);

            RealScalar delta = z.col(0).dot(tmp.col(0));
            vNew = tmp - (delta/gamma)*v - (gamma/gammaPrew)*vPrew;
            precond.apply(vNew, zNew);
            gammaNew = math::sqrt(zNew.col(0).dot(vNew.col(0)));
            RealScalar a0 = c*delta - cPrew*s*gamma;
            RealScalar a1 = math::sqrt(a0*a0 + gammaNew*gammaNew);
            RealScalar a2 = s*delta + cPrew*c*gamma;
            RealScalar a3 = sPrew*gamma;
            cNew = a0/a1;
            sNew = gammaNew/a1;
            wNew = (z - a3*wPrew - a2*w)/a1;
            x = xPrew + cNew*eta*wNew;
            eta = -sNew*eta;

            //Test for convergence
             applyMatrix(mat,x,tmp2);
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

    int iterations() const { return numIter; }
    RealScalar error() const { return m_error; }

private:
    const MatrixType& mat;

    RealScalar tol;
    int maxIters;
    int numIter;
    gsMatrix<real_t> vPrew, v, vNew, wPrew, w, wNew,zNew, z,xPrew, m_rhs, residual, tmp, tmp2;
    RealScalar residualNorm2, threshold, rhsNorm2, m_error;
    RealScalar eta,gammaPrew,gamma,gammaNew,sPrew,s,sNew,cPrew,c,cNew;
};

} // namespace gismo
