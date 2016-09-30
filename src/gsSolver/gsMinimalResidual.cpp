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

    vPrev.setZero(n,m); vNew.setZero(n,m);
    wPrev.setZero(n,m); w.setZero(n,m); wNew.setZero(n,m);
    if (!m_inexact_residual)
        { AwPrev.setZero(n,m); Aw.setZero(n,m); AwNew.setZero(n,m); }

    m_mat->apply(x,negResidual);
    negResidual -= rhs;

    m_error = negResidual.norm() / m_rhs_norm;
    if (m_error < m_tol)
        return true;

    v = -negResidual;
    m_precond->apply(v, z);

    gammaPrev = 1; gamma = math::sqrt(z.col(0).dot(v.col(0))); gammaNew = 1;
    eta = gamma;
    sPrev = 0; s = 0; sNew = 0;
    cPrev = 1; c = 1; cNew = 1;

    return false;
}


bool gsMinimalResidual::step( gsMinimalResidual::VectorType& x )
{
    z /= gamma;
    m_mat->apply(z,Az);

    real_t delta = z.col(0).dot(Az.col(0));
    vNew = Az - (delta/gamma)*v - (gamma/gammaPrev)*vPrev;
    m_precond->apply(vNew, zNew);
    gammaNew = math::sqrt(zNew.col(0).dot(vNew.col(0)));
    const real_t a0 = c*delta - cPrev*s*gamma;
    const real_t a1 = math::sqrt(a0*a0 + gammaNew*gammaNew);
    const real_t a2 = s*delta + cPrev*c*gamma;
    const real_t a3 = sPrev*gamma;
    cNew = a0/a1;
    sNew = gammaNew/a1;
    wNew = (z - a3*wPrev - a2*w)/a1;
    if (!m_inexact_residual)
        AwNew = (Az - a3*AwPrev - a2*Aw)/a1;
    x += cNew*eta*wNew;
    if (!m_inexact_residual)
        negResidual += cNew*eta*AwNew;

    if (m_inexact_residual)
        m_error *= std::fabs(sNew); // see https://eigen.tuxfamily.org/dox-devel/unsupported/MINRES_8h_source.html
    else
        m_error = negResidual.norm() / m_rhs_norm;

    eta = -sNew*eta;
    
    // Test for convergence
    if (m_error < m_tol)
        return true;

    //Update variables
    vPrev.swap(v); v.swap(vNew);     // for us the same as: vPrev = v; v = vNew;
    wPrev.swap(w); w.swap(wNew);     // for us the same as: wPrev = w; w = wNew;
    if (!m_inexact_residual)
        { AwPrev.swap(Aw); Aw.swap(AwNew); } // for us the same as: AwPrev = Aw; Aw = AwNew;
    z.swap(zNew);                    // for us the same as: z = zNew;
    gammaPrev = gamma; gamma = gammaNew;
    sPrev = s; s = sNew;
    cPrev = c; c = cNew;
    return false;
}


void gsMinimalResidual::finalizeIteration( gsMinimalResidual::VectorType& x )
{
    GISMO_UNUSED(x);
    // cleanup temporaries
    negResidual.clear();
    vPrev.clear(); v.clear(); vNew.clear();
    wPrev.clear(); w.clear(); wNew.clear();
    AwPrev.clear(); Aw.clear(); AwNew.clear();
    zNew.clear(); z.clear(); Az.clear();
}



}

