/** @file symm_expr.h

    @brief Defines the symm expression

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
 * @brief Expression for the symmetrization operation
 *        Evaluates
 *        \f[
 *          \text{symm}(u) = (u * u^T)
 *        \f]
 * @ingroup Expressions
 * @tparam E The type of the expression
 */
template <typename E>
class symm_expr : public _expr<symm_expr<E> >
{
    typename E::Nested_t _u;

    mutable gsMatrix<typename E::Scalar> tmp;
public:
    typedef typename E::Scalar Scalar;

    enum { Space = (0==E::Space ? 0 : E::Space), ScalarValued= E::ScalarValued, ColBlocks= E::ColBlocks };

    symm_expr(_expr<E> const& u)
    : _u(u) { }

    MatExprType eval(const index_t k) const
    {
        //const MatExprType tmp = _u.eval(k);
        tmp = _u.eval(k);
        // todo: avoid temporary or not ?
        return tmp * tmp.transpose();
    }

    index_t rows() const { return _u.rows(); }
    index_t cols() const { return _u.rows(); }

    void parse(gsExprHelper<Scalar> & evList) const
    { _u.parse(evList); }

    const gsFeSpace<Scalar> & rowVar() const { return _u.rowVar(); }
    const gsFeSpace<Scalar> & colVar() const { return _u.rowVar(); }

    void print(std::ostream &os) const { os << "symm("; _u.print(os); os <<")"; }
};

}// namespace expr
}// namespace gismo