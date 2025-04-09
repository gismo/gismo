/** @file integral_expr.h

    @brief Defines an expression as an integral over an element

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris
               H.M. Verhelst
*/

#pragma once

#include <gsExpressions/gsFeElement.h>

namespace gismo
{
namespace expr
{

/**
 * @brief Expression for the element diameter
 * @ingroup Expressions
 * @tparam E The expression type
 */
template<class E>
class integral_expr : public _expr<integral_expr<E> >
{
public:
    //typedef typename E::Scalar Scalar;
    typedef real_t Scalar;
    mutable Scalar m_val;
private:
    const gsFeElement<Scalar> & _e; ///<Reference to the element
    typename _expr<E>::Nested_t _ff;
public:
    enum {Space= 0, ScalarValued= 1, ColBlocks = 0};

    integral_expr(const gsFeElement<Scalar> & el, const _expr<E> & u)
    : m_val(-1), _e(el), _ff(u) { }

    const Scalar & eval(const index_t k) const
    {
        GISMO_UNUSED(k);
        GISMO_ENSURE(_e.isValid(), "Element is valid within integrals only.");
        // if (0==k)
        {
            const Scalar * w = _e.weights().data();
            m_val = (*w) * _ff.val().eval(0);
            for (index_t j = 1; j != _e.weights().rows(); ++j)
                m_val += (*(++w)) * _ff.val().eval(j);
        }
        return m_val;
    }

    inline const integral_expr<E> & val() const { return *this; }
    inline index_t rows() const { return 0; }
    inline index_t cols() const { return 0; }
    void parse(gsExprHelper<Scalar> & evList) const
    {
        _ff.parse(evList);
    }

    const gsFeSpace<Scalar> & rowVar() const { return gsNullExpr<Scalar>::get(); }
    const gsFeSpace<Scalar> & colVar() const { return gsNullExpr<Scalar>::get(); }

    void print(std::ostream &os) const
    {
        os << "integral(";
        _ff.print(os);
        os <<")";
    }
};

}// namespace expr
}// namespace gismo