/** @file gsFieldDistance.hpp

    @brief Implementation of the gsField distance functions.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): J. Vogl
*/

#pragma once

#include <gsAssembler/gsFieldDistance.h>
#include <gsExpressions/gsExprHelper.h>
#include <gsAssembler/gsExprEvaluator.h>
#include <gsDomain/gsDomainIterator.h>

namespace gismo
{

template<class T>
T distanceL2(const gsField<T> & field, const gsFunctionSet<T> & func,
             const gsMultiBasis<T> & B, bool isFunc_param)
{
    gsExprEvaluator<T> ev;
    ev.setIntegrationDomain(B.domain());
    auto G = ev.getMap(field.patches());
    auto f1a = ev.getVariable(field.fields());
    auto f1b = ev.getVariable(field.fields(), G);
    auto f2a = ev.getVariable(func);
    auto f2b = ev.getVariable(func, G);
    const bool param = field.isParametric();
    if (param && isFunc_param) return math::sqrt( ev.integral((f1a-f2a).sqNorm() * meas(G)) );
    if (param)                 return math::sqrt( ev.integral((f1a-f2b).sqNorm() * meas(G)) );
    if (isFunc_param)          return math::sqrt( ev.integral((f1b-f2a).sqNorm() * meas(G)) );
    return math::sqrt( ev.integral((f1b-f2b).sqNorm() * meas(G)) );
}

template<class T>
T distanceH1(const gsField<T> & field, const gsFunctionSet<T> & func,
             const gsMultiBasis<T> & B, bool isFunc_param)
{
    gsExprEvaluator<T> ev;
    ev.setIntegrationDomain(B.domain());
    auto G = ev.getMap(field.patches());
    auto f1a = ev.getVariable(field.fields());
    auto f1b = ev.getVariable(field.fields(), G);
    auto f2a = ev.getVariable(func);
    auto f2b = ev.getVariable(func, G);
    const bool param = field.isParametric();
    if (param && isFunc_param)
        return math::sqrt(ev.integral( ( igrad(f1a,G) - igrad(f2a,G)).sqNorm()*meas(G) ) );
    if (param)
        return math::sqrt(ev.integral( ( igrad(f1a,G) - igrad(f2b)).sqNorm()*meas(G) ) );
    if (isFunc_param)
        return math::sqrt(ev.integral( ( igrad(f1b) - igrad(f2a,G)).sqNorm()*meas(G) ) );
    return math::sqrt(ev.integral( ( igrad(f1b) - igrad(f2b)).sqNorm()*meas(G) ) );
}

template<class T>
T distanceH2(const gsField<T> & field, const gsFunctionSet<T> & func,
             bool isFunc_param)
{
    const gsMultiPatch<T> & mp = field.patches();
    gsMultiBasis<T> mb;
    if (const gsMultiPatch<T>* imp = dynamic_cast<const gsMultiPatch<T>*>(&field.fields()))
        mb = gsMultiBasis<T>(*imp);
    else
        mb = gsMultiBasis<T>(mp);

    gsExprEvaluator<T> ev;
    ev.setIntegrationDomain(mb.domain());
    auto G = ev.getMap(field.patches());
    auto f1a = ev.getVariable(field.fields());
    auto f1b = ev.getVariable(field.fields(), G);
    auto f2a = ev.getVariable(func);
    auto f2b = ev.getVariable(func, G);
    const bool param = field.isParametric();
    if (param && isFunc_param)
        return math::sqrt(ev.integral( ( ihess(f1a,G) - ihess(f2a,G)).sqNorm()*meas(G) ) );
    if (param)
        return math::sqrt(ev.integral( (ihess(f1a,G) - ihess(f2b)).sqNorm()*meas(G) ) );
    if (isFunc_param)
        return math::sqrt(ev.integral( ( ihess(f1b) - ihess(f2a,G)).sqNorm()*meas(G) ) );
    return math::sqrt(ev.integral( ( ihess(f1b) - ihess(f2b)).sqNorm()*meas(G) ) );
}

template<class T>
T distanceL2(const gsField<T> & field, const gsField<T> & other)
{
    return distanceL2(field, other.fields(), other.isParametric());
}

template<class T>
T distanceL2(const gsField<T> & field, const gsFunctionSet<T> & func,
             bool isFunc_param)
{
    if (const gsMultiPatch<T>* mp = dynamic_cast<const gsMultiPatch<T>*>(&field.fields()))
        return distanceL2(field, func, gsMultiBasis<T>(*mp), isFunc_param);
    gsMultiBasis<T> mb(field.patches());
    return distanceL2(field, func, mb, isFunc_param);
}

template<class T>
T distanceH1(const gsField<T> & field, const gsFunctionSet<T> & func,
             bool isFunc_param)
{
    if (const gsMultiPatch<T>* mp = dynamic_cast<const gsMultiPatch<T>*>(&field.fields()))
        return distanceH1(field, func, gsMultiBasis<T>(*mp), isFunc_param);
    gsMultiBasis<T> mb(field.patches());
    return distanceH1(field, func, mb, isFunc_param);
}

template<class T>
T distanceDG(const gsField<T> &, const gsFunctionSet<T> &, bool)
{
    GISMO_ERROR("distanceDG is not implemented");
}

} // namespace gismo
