/** @file onormal_expr.h

    @brief Defines the onormal expression

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris
               H.M. Verhelst
*/

#pragma once

namespace gismo
{
namespace expr
{

/**
 * @brief Expression for the outer pointing normal of a geometry map. This
 *        expression is valid only at the boundaries of a geometric patch
 * @ingroup Expressions
 * @tparam T The type of the expression
 */
template<class T>
class onormal_expr : public _expr<onormal_expr<T> >
{
    typename gsGeometryMap<T>::Nested_t _G;

public:
    typedef T Scalar;
    enum {Space = 0, ScalarValued= 0, ColBlocks= 0};

    explicit onormal_expr(const gsGeometryMap<T> & G) : _G(G) { }

    auto eval(const index_t k) const -> decltype(_G.data().outNormals.col(k))
    { return _G.data().outNormals.col(k); }

    index_t rows() const { return  _G.source().targetDim(); }
    index_t cols() const { return 1; }

    const gsFeSpace<T> & rowVar() const {return gsNullExpr<T>::get();}
    const gsFeSpace<T> & colVar() const {return gsNullExpr<T>::get();}

    void parse(gsExprHelper<Scalar> & evList) const
    {
        evList.add(_G);
        _G.data().flags |= NEED_OUTER_NORMAL;
    }

    void print(std::ostream &os) const { os << "nv("; _G.print(os); os <<")"; }
};

/**
 * @brief Expression for the outer pointing normal of a geometry map. This
 *        expression is valid only at the boundaries of a geometric patch
 * @ingroup Expressions
 * @tparam T The type of the expression
 */
template<class T> EIGEN_STRONG_INLINE
onormal_expr<T> nv(const gsGeometryMap<T> & u) { return onormal_expr<T>(u); }

}// namespace expr
}// namespace gismo