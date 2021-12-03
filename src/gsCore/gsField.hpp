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

template <class T>
T gsField<T>::distanceL2(gsFunctionSet<T> const & func,
                         gsMultiBasis<T> const & B,
                         bool isFunc_param,
                         int numEvals) const
{
    GISMO_UNUSED(numEvals);// todo: subdivided quadrature elements
    gsExprEvaluator<T> ev;
    ev.setIntegrationElements(B);
    auto G = ev.getMap(this->patches());
    auto f1a = ev.getVariable(*m_fields);
    auto f1b = ev.getVariable(*m_fields, G);
    auto f2a = ev.getVariable(func);
    auto f2b = ev.getVariable(func, G);
    if (m_parametric && isFunc_param) return math::sqrt( ev.integral((f1a-f2a).sqNorm() * meas(G)) );
    if (m_parametric) return math::sqrt( ev.integral((f1a-f2b).sqNorm() * meas(G)) );
    if (isFunc_param) return math::sqrt( ev.integral((f1b-f2a).sqNorm() * meas(G)) );
    return math::sqrt( ev.integral((f1b-f2b).sqNorm() * meas(G)) );
}

template <class T>
T gsField<T>::distanceH1(gsFunctionSet<T> const & func,
                         gsMultiBasis<T> const & B,
                         bool isFunc_param,
                         int) const
{
    gsExprEvaluator<T> ev;
    ev.setIntegrationElements(B);
    auto G = ev.getMap(this->patches());
    auto f1a = ev.getVariable(*m_fields);
    auto f1b = ev.getVariable(*m_fields, G);
    auto f2a = ev.getVariable(func);
    auto f2b = ev.getVariable(func, G);
    if (m_parametric && isFunc_param)
        return math::sqrt(ev.integral( ( igrad(f1a,G) - igrad(f2a,G)).sqNorm()*meas(G) ) );
    if (m_parametric)
        return math::sqrt(ev.integral( ( igrad(f1a,G) - igrad(f2b)).sqNorm()*meas(G) ) );
    if (isFunc_param)
        return math::sqrt(ev.integral( ( igrad(f1b) - igrad(f2a,G)).sqNorm()*meas(G) ) );
    return math::sqrt(ev.integral( ( igrad(f1b) - igrad(f2b)).sqNorm()*meas(G) ) );
}

template <class T>
T gsField<T>::distanceH2(gsFunctionSet<T> const & func,
                         bool isFunc_param) const
{
    const gsMultiPatch<T> & mp = this->patches();
    gsMultiBasis<T> mb;
    if (const gsMultiPatch<T>* imp = dynamic_cast<const gsMultiPatch<T>*>(m_fields.get()))
        mb = gsMultiBasis<T>(*imp);
    else
        mb = gsMultiBasis<T>(mp);

    gsExprEvaluator<T> ev;
    ev.setIntegrationElements(mb);
    auto G = ev.getMap(this->patches());
    auto f1a = ev.getVariable(*m_fields);
    auto f1b = ev.getVariable(*m_fields, G);
    auto f2a = ev.getVariable(func);
    auto f2b = ev.getVariable(func, G);

    if (m_parametric && isFunc_param)
        return math::sqrt(ev.integral( ( ihess(f1a,G) - ihess(f2a,G)).sqNorm()*meas(G) ) );
    if (m_parametric)
        return math::sqrt(ev.integral( (ihess(f1a,G) - ihess(f2b)).sqNorm()*meas(G) ) );
    if (isFunc_param)
        return math::sqrt(ev.integral( ( ihess(f1b) - ihess(f2a,G)).sqNorm()*meas(G) ) );
    return math::sqrt(ev.integral( ( ihess(f1b) - ihess(f2b)).sqNorm()*meas(G) ) );
}

template <class T>
T gsField<T>::distanceL2(gsField<T> const & field, int numEvals) const
{
    return distanceL2(*field.m_fields, field.m_parametric, numEvals);
}

template <class T>
T gsField<T>::distanceL2(gsFunctionSet<T> const & func,
                         bool isFunc_param,
                         int numEvals) const
{
    if (const gsMultiPatch<T>* mp = dynamic_cast<const gsMultiPatch<T>*>(m_fields.get()))
        return distanceL2(func, gsMultiBasis<T>(*mp), isFunc_param, numEvals);
    gsMultiBasis<T> mb(this->patches());
    return distanceL2(func, mb, isFunc_param, numEvals);
}

template <class T>
T gsField<T>::distanceH1(gsFunctionSet<T> const & func,
                         bool isFunc_param,
                         int) const
{
    if (const gsMultiPatch<T>* mp = dynamic_cast<const gsMultiPatch<T>*>(m_fields.get()))
        return distanceH1(func, gsMultiBasis<T>(*mp), isFunc_param);
    gsMultiBasis<T> mb(this->patches());
    return distanceH1(func, mb, isFunc_param);
}

template <class T>
T gsField<T>::distanceDG(gsFunctionSet<T>
                         const & ,
                         bool ,
                         int) const
{
    GISMO_NO_IMPLEMENTATION
    /*
    if (m_parametric) // isogeometric field
        return 0; //igaFieldDGDistance(*this, func, isFunc_param);
    else
    {
        gsWarn << "DG norm not implemented.\n";
        return -1;
    }
    */
}

} // namespace gismo
