/** @file ppart_expr.h

    @brief Defines the positive part

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
 * @brief Expression for the component-wise positive part
 *        Evaluates
 *       \f[
 *          \text{ppart}(u) = \max(u, 0)
 *       \f]
 * @ingroup Expressions
 * @tparam E The type of the expression
 */
template<class E>
class ppart_expr : public _expr<ppart_expr<E> >
{
public:
    typedef typename E::Scalar Scalar;
    enum {ScalarValued = E::ScalarValued, Space = E::Space, ColBlocks= E::ColBlocks};
private:
    typename E::Nested_t _u;
    mutable gsMatrix<Scalar> res;
public:

    ppart_expr(_expr<E> const& u) : _u(u) { }

    const gsMatrix<Scalar> & eval(index_t k) const
    {
        res = _u.eval(k).cwiseMax(0.0); // component-wise maximum with zero
        return res;
    }


    index_t rows() const { return _u.rows(); }
    index_t cols() const { return _u.cols(); }

    void parse(gsExprHelper<Scalar> & el) const
    { _u.parse(el); }

    const gsFeSpace<Scalar> & rowVar() const {return _u.rowVar();}
    const gsFeSpace<Scalar> & colVar() const {return _u.colVar();}

    void print(std::ostream &os) const { os<<"posPart("; _u.print(os); os <<")"; }
};


}// namespace expr
}// namespace gismo