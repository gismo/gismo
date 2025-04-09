/** @file meas_expr.h

    @brief Defines the measure expression

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
 * @brief Expression for the measure of a geometry map
 * @ingroup Expressions
 * @tparam T The expression type
 */
template<class T>
class meas_expr : public _expr<meas_expr<T> >
{
    typename gsGeometryMap<T>::Nested_t _G;

public:
    enum {Space = 0, ScalarValued = 1, ColBlocks = 0};

    typedef T Scalar;

    meas_expr(const gsGeometryMap<T> & G) : _G(G) { }

    T eval(const index_t k) const
    {
        return _G.data().measures.at(k);
    }

    index_t rows() const { return 0; }
    index_t cols() const { return 0; }
    void parse(gsExprHelper<Scalar> & evList) const
    {
        evList.add(_G);
        _G.data().flags |= NEED_MEASURE;
    }

    const gsFeSpace<T> & rowVar() const { return gsNullExpr<T>::get(); }
    const gsFeSpace<T> & colVar() const { return gsNullExpr<T>::get(); }

    void print(std::ostream &os) const { os << "meas("; _G.print(os); os <<")"; }
};

/**
 * @brief Function to create a measure expression
 * @ingroup Expressions
 * @param G The geometry map
 */
template<class T> EIGEN_STRONG_INLINE
meas_expr<T> meas(const gsGeometryMap<T> & G) { return meas_expr<T>(G); }

}// namespace expr
}// namespace gismo