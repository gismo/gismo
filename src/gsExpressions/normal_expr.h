/** @file normal_expr.h

    @brief Defines the normal expression

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
 * @brief Expression for the out of plane surface normal of a geometry map
 * @ingroup Expressions
 * @tparam T The type of the expression
 */
template<class T>
class normal_expr : public _expr<normal_expr<T> >
{
    typename gsGeometryMap<T>::Nested_t _G;

public:
    typedef T Scalar;
    enum {Space = 0, ScalarValued= 0, ColBlocks= 0};

    normal_expr(const gsGeometryMap<T> & G) : _G(G)
    {
        GISMO_ENSURE( _G.source().domainDim()+1 == _G.source().targetDim(), "Surface normal requires codimension 1");
    }

    auto eval(const index_t k) const -> decltype(_G.data().normals.col(k))
    { return _G.data().normals.col(k); }

    index_t rows() const { return _G.source().targetDim(); }
    index_t cols() const { return 1; }

    const gsFeSpace<T> & rowVar() const {return gsNullExpr<T>::get();}
    const gsFeSpace<T> & colVar() const {return gsNullExpr<T>::get();}

    void parse(gsExprHelper<Scalar> & evList) const
    {
        evList.add(_G);
        _G.data().flags |= NEED_NORMAL;
    }

    void print(std::ostream &os) const { os << "sn("; _G.print(os); os <<")"; }
};

/**
 * @brief Expression for the out of plane surface normal of a geometry map
 * @ingroup Expressions
 * @tparam T The type of the expression
 */
template<class T> EIGEN_STRONG_INLINE
normal_expr<T> sn(const gsGeometryMap<T> & u) { return normal_expr<T>(u); }

}// namespace expr
}// namespace gismo