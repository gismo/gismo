/** @file gsGradientMethod.hpp

    @brief Gradient iteration

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): S. Takacs
*/

#include <gsSolver/gsGradientMethod.h>

namespace gismo
{

template<class T>
bool gsGradientMethod<T>::initIteration( const typename gsGradientMethod<T>::VectorType& rhs,
                                            typename gsGradientMethod<T>::VectorType& x )
{
    if (Base::initIteration(rhs,x))
        return true;

    m_mat->apply(x,m_tmp);
    m_res = rhs - m_tmp;

    m_error = m_res.norm() / m_rhs_norm;
    return m_error < m_tol;

}

template<class T>
bool gsGradientMethod<T>::step( typename gsGradientMethod<T>::VectorType& x )
{
    m_precond->apply(m_res,m_update);
    m_mat->apply(m_update,m_tmp);
    x += m_damping * m_update;
    m_res -= m_damping * m_tmp;
    m_error = m_res.norm() / m_rhs_norm;
    return m_error < m_tol;
}


} // end namespace gismo
