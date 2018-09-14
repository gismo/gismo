/** @file gsField.hpp

    @brief Provides implementation of the Field class.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): J. Vogl
*/

#pragma once

#include <gsAssembler/gsExprHelper.h>
#include <gsAssembler/gsExprEvaluator.h>
#include <gsCore/gsDomainIterator.h>

namespace gismo
{


/// Computes the L2-distance between the two fields, on the physical domain
template <class T>
T gsField<T>::distanceL2(gsField<T> const & field, int numEvals) const
{
    const gsMultiPatch<T> & mp = this->patches();
    gsMultiBasis<T> mb(mp);
    gsExprEvaluator<T> ev;
    ev.setIntegrationElements(mb);
    typename gsExprEvaluator<T>::geometryMap G = ev.getMap(mp);
    typename gsExprEvaluator<T>::variable f1   = // getCoeff
        (m_parametric ? ev.getVariable(*m_fields) : ev.getVariable(*m_fields, G) );
    typename gsExprEvaluator<T>::variable f2   =
        (field.m_parametric ? ev.getVariable(*field.m_fields) : ev.getVariable(*field.m_fields, G) );
    return math::sqrt( ev.integral((f1 - f2).sqNorm() * meas(G)) );
}

/// Computes the L2-distance between the field and a function \a func on the physical domain
template <class T>
T gsField<T>::distanceL2(gsFunction<T> const & func,
                         bool isFunc_param,
                         int numEvals) const
{
    gsMultiBasis<T> mb(this->patches());
    return distanceL2(func, mb, isFunc_param, numEvals);
}

/// Computes the L2-distance between the field and a function \a
/// func on the physical domain, using mesh from B
template <class T>
T gsField<T>::distanceL2(gsFunction<T> const & func,
                         gsMultiBasis<T> const & B,
                         bool isFunc_param,
                         int numEvals) const
{
    GISMO_UNUSED(numEvals);// todo: subdivided quadrature elements
    const gsMultiPatch<T> & mp = this->patches();
    gsExprEvaluator<T> ev;
    ev.setIntegrationElements(B);
    typename gsExprEvaluator<T>::geometryMap G = ev.getMap(mp);
    typename gsExprEvaluator<T>::variable f1   =
        (m_parametric ? ev.getVariable(*m_fields) : ev.getVariable(*m_fields, G) );
    typename gsExprEvaluator<T>::variable f2   =
        (isFunc_param ? ev.getVariable(func) : ev.getVariable(func, G) );
    gsDebugVar(m_parametric);
    gsDebugVar(isFunc_param);
    return math::sqrt( ev.integral((f1 - f2).sqNorm() * meas(G)) );
}

/// Computes the H1-distance between the field and a function \a
/// func on the physical domain
template <class T>
T gsField<T>::distanceH1(gsFunction<T> const & func,
                         bool isFunc_param,
                         int) const
{
    gsMultiBasis<T> mb(this->patches());
    return distanceH1(func, mb, isFunc_param);
}

/// Computes the H1-distance between the field and a function \a
/// func on the physical domain, using mesh from B
template <class T>
T gsField<T>::distanceH1(gsFunction<T> const & func,
                         gsMultiBasis<T> const & B,
                         bool isFunc_param,
                         int) const
{
    const gsMultiPatch<T> & mp = this->patches();
    gsExprEvaluator<T> ev;
    ev.setIntegrationElements(B);
    typename gsExprEvaluator<T>::geometryMap G = ev.getMap(mp);
    typename gsExprEvaluator<T>::variable f1   =
        (m_parametric ? ev.getVariable(*m_fields) : ev.getVariable(*m_fields, G) );
    typename gsExprEvaluator<T>::variable f2   =
        (isFunc_param ? ev.getVariable(func) : ev.getVariable(func, G) );

    if (m_parametric && isFunc_param)
        return math::sqrt(ev.integral( ( igrad(f1,G) - igrad(f2,G)).sqNorm()*meas(G) ) );
    if (m_parametric)
        return math::sqrt(ev.integral( ( igrad(f1,G) - igrad(f2)).sqNorm()*meas(G) ) );
    if (isFunc_param)
        return math::sqrt(ev.integral( ( igrad(f1) - igrad(f2,G)).sqNorm()*meas(G) ) );
    return math::sqrt(ev.integral( ( igrad(f1) - igrad(f2)).sqNorm()*meas(G) ) );
}

/// Computes the DG-distance between the field and a function \a
/// func on the physical domain
template <class T>
T gsField<T>::distanceDG(gsFunction<T>
                         const & func,
                         bool isFunc_param,
                         int) const
{
    if (m_parametric) // isogeometric field
        return 0; //igaFieldDGDistance(*this, func, isFunc_param);
    else
    {
        gsWarn << "DG norm not implemented.\n";
        return -1;
    }
}

} // namespace gismo
