/** @file ppartval_expr.h

    @brief Defines the positive part (scalar)

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
class ppartval_expr : public _expr<ppartval_expr<E> >
{
  typename E::Nested_t _u;
 public:
  typedef typename E::Scalar Scalar;
  enum {ScalarValued = 1, Space = 0, ColBlocks= 0};
  mutable Scalar res;
 public:

  ppartval_expr(_expr<E> const& u) : _u(u) { }

  Scalar & eval(index_t k) const
  {
    res = std::max(0.0,_u.eval(k));
    return res; // component-wise maximum with zero
  }

  index_t rows() const { return 0; }
  index_t cols() const { return 0; }

  void parse(gsExprHelper<Scalar> & evList) const
  { _u.parse(evList); }

  const gsFeSpace<Scalar> & rowVar() const {return gsNullExpr<Scalar>::get();}
  const gsFeSpace<Scalar> & colVar() const {return gsNullExpr<Scalar>::get();}

  void print(std::ostream &os) const { os<<"posPart("; _u.print(os); os <<")"; }
};



}// namespace expr
}// namespace gismo