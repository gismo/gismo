/** @file col_expr.h

    @brief Defines the column expression

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
 * @brief Expression for the column of a matrix
 * @ingroup Expressions
 * @tparam E The expression type
 */
template<class E>
class col_expr : public _expr<col_expr<E> >
{
    typename E::Nested_t _c;
    const index_t _i;
public:
    typedef typename E::Scalar Scalar;
    typedef const col_expr<E> Nested_t;

    enum { Space = E::Space, ScalarValued = 0, ColBlocks = 0 };

    col_expr(const E & c, const index_t i) : _c(c), _i(i) { }

public:

    //ConstColXpr
    inline MatExprType eval(const index_t k) const { return _c.eval(k).col(_i); }

    index_t rows() const { return _c.rows(); }
    index_t cols() const { return 1; }
    void parse(gsExprHelper<Scalar> & evList) const { _c.parse(evList); }

    const gsFeSpace<Scalar> & rowVar() const { return _c.rowVar(); }
    const gsFeSpace<Scalar> & colVar() const { return _c.colVar(); }

    void print(std::ostream &os) const { os<<_c<<"["<<_i<<"]"; }
};

}// namespace expr
}// namespace gismo