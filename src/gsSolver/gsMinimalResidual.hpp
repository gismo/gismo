/** @file gsMinimalResidual.hpp

    @brief Preconditioned iterative solver using the minimal residual method.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): J. Sogn, S. Takacs
*/

namespace gismo
{

template<class T>
bool gsMinimalResidual<T>::initIteration( const typename gsMinimalResidual<T>::VectorType& rhs, typename gsMinimalResidual<T>::VectorType& x )
{
    Base::initIteration(rhs,x);
    //if (Base::initIteration(rhs,x)) return true; // z will not be initialized!

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

template<class T>
bool gsMinimalResidual<T>::step( typename gsMinimalResidual<T>::VectorType& x )
{
    z /= gamma;
    m_mat->apply(z,Az);

    T delta = z.col(0).dot(Az.col(0));
    vNew = Az - (delta/gamma)*v - (gamma/gammaPrev)*vPrev;
    m_precond->apply(vNew, zNew);
    gammaNew = math::sqrt(zNew.col(0).dot(vNew.col(0)));
    const T a0 = c*delta - cPrev*s*gamma;
    const T a1 = math::sqrt(a0*a0 + gammaNew*gammaNew);
    const T a2 = s*delta + cPrev*c*gamma;
    const T a3 = sPrev*gamma;
    cNew = a0/a1;
    sNew = gammaNew/a1;
    wNew = (z - a3*wPrev - a2*w)/a1;
    if (!m_inexact_residual)
        AwNew = (Az - a3*AwPrev - a2*Aw)/a1;
    x += cNew*eta*wNew;
    if (!m_inexact_residual)
        negResidual += cNew*eta*AwNew;

    if (m_inexact_residual)
        m_error *= math::abs(sNew); // see https://eigen.tuxfamily.org/dox-devel/unsupported/MINRES_8h_source.html
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

template<class T>
void gsMinimalResidual<T>::finalizeIteration( typename gsMinimalResidual<T>::VectorType& )
{
    // cleanup temporaries
    negResidual.clear();
    vPrev.clear(); v.clear(); vNew.clear();
    wPrev.clear(); w.clear(); wNew.clear();
    AwPrev.clear(); Aw.clear(); AwNew.clear();
    zNew.clear(); z.clear(); Az.clear();
}



}
