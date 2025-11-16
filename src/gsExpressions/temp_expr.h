/** @file temp_expr.h

    @brief Defines the temp expression

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
 * @brief Expression for an evaluation of the (sub-)expression in temporary memory
 * @ingroup Expressions
 * @tparam E The type of the expression
 */
template<class E>
class temp_expr : public _expr<temp_expr<E> >
{
    typename E::Nested_t _u;
    typedef typename E::Scalar Scalar;
    mutable gsMatrix<Scalar> tmp;

public:
    temp_expr(_expr<E> const& u)
    : _u(u) { }

public:
    enum {Space = E::Space, ScalarValued = E::ScalarValued,
        ColBlocks = E::ColBlocks};

    // template<bool S  = ColBlocks>
    // typename util::enable_if<S,MatExprType>::type
    const gsMatrix<Scalar> & eval(const index_t k) const
    {
        tmp = _u.eval(k);
        return tmp;
    }

    index_t rows() const { return _u.rows(); }
    index_t cols() const { return _u.cols(); }
    void parse(gsExprHelper<Scalar> & evList) const { _u.parse(evList); }
    const gsFeSpace<Scalar> & rowVar() const { return _u.rowVar(); }
    const gsFeSpace<Scalar> & colVar() const { return _u.colVar(); }

    void print(std::ostream &os) const { _u.print(os); }
};

}// namespace expr
}// namespace gismo