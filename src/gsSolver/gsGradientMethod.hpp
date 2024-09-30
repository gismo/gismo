/** @file gsGradientMethod.hpp

    @brief Gradient iteration

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): S. Takacs
*/

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

    T step_size;
    if (m_adapt_step_size)
        step_size = m_tmp.col(0).dot(m_res.col(0)) / m_tmp.col(0).dot(m_tmp.col(0));
    else
        step_size = m_step_size;

    x += step_size * m_update;
    m_res -= step_size * m_tmp;
    m_error = m_res.norm() / m_rhs_norm;
    return m_error < m_tol;
}


} // end namespace gismo
