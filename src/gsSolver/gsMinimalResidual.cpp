/** @file gsMinimalResidual.cpp

    @brief Preconditioned iterative solver using the minimal residual method.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): J. Sogn
*/
#include <gsSolver/gsMinimalResidual.h>

namespace gismo
{

bool gsMinimalResidual::initIteration( const gsMinimalResidual::VectorType& rhs, gsMinimalResidual::VectorType& x )
{
    if (Base::initIteration(rhs,x))
        return true;
    
    int n = m_mat->cols();
    int m = 1;//rhs.cols();
    m_rhs = rhs;

    vPrew.setZero(n,m); vNew.setZero(n,m);
    wPrew.setZero(n,m); w.setZero(n,m); wNew.setZero(n,m);
    tmp2.setZero(n,1);

    m_mat->apply(x,tmp2);
    v = m_rhs - tmp2;

    m_precond->apply(v, z);

    gammaPrew = 1; gamma = math::sqrt(z.col(0).dot(v.col(0))); gammaNew = 1;
    eta = gamma;
    sPrew = 0; s = 0; sNew = 0;
    cPrew = 1; c = 1; cNew = 1;
    
    return false;
}


bool gsMinimalResidual::step( gsMinimalResidual::VectorType& x )
{
    z /= gamma;
    m_mat->apply(z,tmp);

    real_t delta = z.col(0).dot(tmp.col(0));
    vNew = tmp - (delta/gamma)*v - (gamma/gammaPrew)*vPrew;
    m_precond->apply(vNew, zNew);
    gammaNew = math::sqrt(zNew.col(0).dot(vNew.col(0)));
    real_t a0 = c*delta - cPrew*s*gamma;
    real_t a1 = math::sqrt(a0*a0 + gammaNew*gammaNew);
    real_t a2 = s*delta + cPrew*c*gamma;
    real_t a3 = sPrew*gamma;
    cNew = a0/a1;
    sNew = gammaNew/a1;
    wNew = (z - a3*wPrew - a2*w)/a1;
    x += cNew*eta*wNew;
    eta = -sNew*eta;

    //Test for convergence
    m_mat->apply(x,tmp2);
    residual = m_rhs - tmp2;
    m_error = residual.norm() / m_rhs_norm;
    if (m_error < m_tol)
        return true;

    //Update variables
    vPrew = v; v = vNew;
    wPrew = w; w = wNew;
    z = zNew;
    gammaPrew = gamma; gamma = gammaNew;
    sPrew = s; s = sNew;
    cPrew = c; c = cNew;
    return false;
}

}

