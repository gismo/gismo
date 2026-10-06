/** @file nabla2_expr.h

    @brief Defines the nabla2 expression

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
 * @brief Expression for the nabla2 (\f$\nabla^2\f$) of a finite element variable
 *        see also https://en.wikipedia.org/wiki/Del
 * @note  Transposed pure second derivatives are returned as a matrix
 * @ingroup Expressions
 * @tparam T The type of the expression
 */
template<class T>
class nabla2_expr : public _expr<nabla2_expr<T> >
{
    typename gsFeVariable<T>::Nested_t u;

public:
    typedef T Scalar;
    enum{Space = 1};

    /* // todo
       nabla2_expr(const gsGeometryMap<T> & G)
       : m_data(G.data()) { }
    */

    nabla2_expr(const gsFeVariable<T> & _u)
    : u(_u)
    { }

    MatExprType eval(const index_t k) const
    {
        // numActive x parDim
        return u.data().values[2]
            .reShapeCol(k, u.data().values[2].rows()/u.cSize(), u.cSize() )
            .topRows(u.parDim()).transpose();
    }

    index_t rows() const { return u.rows();   }
    index_t cols() const { return u.parDim(); }

    void parse(gsExprHelper<Scalar> & evList) const
    {
        evList.add(u);
        u.data().flags |= NEED_DERIV2;
    }

    const gsFeSpace<T> & rowVar() const { return u.rowVar(); }
    const gsFeSpace<T> & colVar() const
    {return gsNullExpr<T>::get();}
};

/**
 * @brief Expression for the nabla2 (\f$\nabla^2\f$) of a finite element variable
 * @ingroup Expressions
 * @param u The expression
 */
template<class T>
nabla2_expr<T> nabla2(const gsFeVariable<T> & u) { return nabla2_expr<T>(u); }
// #define lapl(x) nabla2(x).sum() // assume tarDim==1

}// namespace expr
}// namespace gismo