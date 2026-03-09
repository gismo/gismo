/** @file exp_expr.h

    @brief Defines the exp expression

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
 * @brief Expression for the exponentiation of a given expression
 * @ingroup Expressions
 * @tparam E The expression type
 */
template<class E>
class exp_expr : public _expr<exp_expr<E> >
{
  typename E::Nested_t _u;
 public:
  typedef typename E::Scalar Scalar;
  enum {ScalarValued = 1, Space = E::Space, ColBlocks= 0};

  exp_expr(_expr<E> const& u) : _u(u) { }

  Scalar eval(const index_t k) const
  {
    const Scalar v = _u.val().eval(k);
    return math::exp(v);
  }

  static index_t rows() { return 0; }
  static index_t cols() { return 0; }

  void parse(gsExprHelper<Scalar> & el) const
  { _u.parse(el); }

  static bool isScalar() { return true; }

  const gsFeSpace<Scalar> & rowVar() const {return gsNullExpr<Scalar>::get();}
  const gsFeSpace<Scalar> & colVar() const {return gsNullExpr<Scalar>::get();}

  void print(std::ostream &os) const { os<<"exp("; _u.print(os); os <<")"; }
};

}// namespace expr
}// namespace gismo