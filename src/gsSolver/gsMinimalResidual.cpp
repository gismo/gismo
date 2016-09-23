/** @file gsMinimalResidual.cpp

    @brief Preconditioned iterative solver using the minimal residual method.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): J. Sogn, S. Takacs
*/
#include <gsSolver/gsMinimalResidual.h>

namespace gismo
{

bool gsMinimalResidual::initIteration( const gsMinimalResidual::VectorType& rhs, gsMinimalResidual::VectorType& x )
{
    if (Base::initIteration(rhs,x))
        return true;
    
    int n = m_mat->cols();
    int m = 1; // = rhs.cols();

    vPrew.setZero(n,m); vNew.setZero(n,m);
    wPrew.setZero(n,m); w.setZero(n,m); wNew.setZero(n,m);
    AwPrew.setZero(n,m); Aw.setZero(n,m); AwNew.setZero(n,m);

    m_mat->apply(x,neg_residual);
    neg_residual -= rhs;

    v = -neg_residual;
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
    m_mat->apply(z,Az);

    real_t delta = z.col(0).dot(Az.col(0));
    vNew = Az - (delta/gamma)*v - (gamma/gammaPrew)*vPrew;
    m_precond->apply(vNew, zNew);
    gammaNew = math::sqrt(zNew.col(0).dot(vNew.col(0)));
    real_t a0 = c*delta - cPrew*s*gamma;
    real_t a1 = math::sqrt(a0*a0 + gammaNew*gammaNew);
    real_t a2 = s*delta + cPrew*c*gamma;
    real_t a3 = sPrew*gamma;
    cNew = a0/a1;
    sNew = gammaNew/a1;
    wNew = (z - a3*wPrew - a2*w)/a1;
    AwNew = (Az - a3*AwPrew - a2*Aw)/a1;
    x += cNew*eta*wNew;
    neg_residual += cNew*eta*AwNew;
    eta = -sNew*eta;

    //Test for convergence
    m_error = neg_residual.norm() / m_rhs_norm;
    if (m_error < m_tol)
        return true;

    //Update variables
    vPrew.swap(v); v.swap(vNew);     // for us the same as: vPrew = v; v = vNew;
    wPrew.swap(w); w.swap(wNew);     // for us the same as: wPrew = w; w = wNew;
    AwPrew.swap(Aw); Aw.swap(AwNew); // for us the same as: AwPrew = Aw; Aw = AwNew;
    z.swap(zNew);                    // for us the same as: z = zNew;
    gammaPrew = gamma; gamma = gammaNew;
    sPrew = s; s = sNew;
    cPrew = c; c = cNew;
    return false;
}

}

