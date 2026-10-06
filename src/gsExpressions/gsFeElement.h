/** @file gsFeElement.h

    @brief Defines an element as an expression

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris
               H.M. Verhelst
*/

#pragma once

#include <gsExpressions/gsGeometryMap.h>

namespace gismo
{
namespace expr
{

/**
 * @brief Expression for a finite element, collecting relevant expressions
 * @ingroup Expressions
 * @tparam T The expression type
 */
template<class T>
class gsFeElement
{
    friend class cdiam_expr<T>;
    const gsExprHelper<T> * m_exprdata;

    gsFeElement(const gsFeElement &);
public:
    typedef T Scalar;

    gsFeElement(const gsExprHelper<T> & eh) : m_exprdata(&eh) { }

    bool isValid() const { return nullptr!=m_exprdata; }

    const gsVector<T> & weights() const {return m_exprdata->weights();}

    template<class E> inline
    integral_expr<E> integral(const _expr<E>& ff) const
    { return integral_expr<E>(*this,ff); }

    typedef integral_expr<T> AreaRetType;
    AreaRetType area() const
    { return integral(_expr<T,true>(1)); }

    typedef integral_expr<meas_expr<T> > PHAreaRetType;
    /// The diameter of the element on the physical space
    PHAreaRetType area(const gsGeometryMap<Scalar> & _G) const
    { return integral(meas_expr<T>(_G)); }

    typedef pow_expr<integral_expr<T> > DiamRetType;
    /// The diameter of the element (on parameter space)
    DiamRetType diam() const //-> int(1)^(1/d)
    { return pow(integral(_expr<T,true>(1)),(T)(1)/(T)(2)); }

    typedef pow_expr<integral_expr<meas_expr<T> > > PHDiamRetType;
    /// The diameter of the element on the physical space
    PHDiamRetType diam(const gsGeometryMap<Scalar> & _G) const
    { return pow(integral(meas_expr<T>(_G)),(T)(1)/(T)(2)); }

    //auto points() const {return point_expr<T>(*this);}
    //index_t dim() { return di->

    void print(std::ostream &os) const { os << "e"; }

    void parse(gsExprHelper<T> & evList) const
    {
        GISMO_ERROR("Call desired member of element expression instead.");
    }
};

}// namespace expr
}// namespace gismo