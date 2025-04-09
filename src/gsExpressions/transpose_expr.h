/** @file transpose_expr.h

    @brief Defines the transpose expression

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
 * @brief Expression for the transpose of an expression
 * @ingroup Expressions
 * @tparam E The type of the expression
 * @tparam cw The direction of the transpose (column-wise or row-wise)
 */
template<class E, bool cw>
class transpose_expr : public _expr<transpose_expr<E,cw> >
{
    typename E::Nested_t _u;

public:

    typedef typename E::Scalar Scalar;

    transpose_expr(_expr<E> const& u)
    : _u(u) { }

public:
    enum {ColBlocks = E::ColBlocks, ScalarValued=E::ScalarValued};
    enum {Space = cw?E::Space:(E::Space==1?2:(E::Space==2?1:E::Space))};

    mutable Temporary_t res;
    const Temporary_t & eval(const index_t k) const
    {
        if (E::ColBlocks)
            res = _u.eval(k).blockTranspose( _u.cardinality() );
        else
            res = _u.eval(k).transpose();
        return res;
    }

    index_t rows() const { return _u.cols(); }

    index_t cols() const { return _u.rows(); }

    void parse(gsExprHelper<Scalar> & evList) const
    { _u.parse(evList); }

    const gsFeSpace<Scalar> & rowVar() const { return cw?_u.rowVar():_u.colVar(); }
    const gsFeSpace<Scalar> & colVar() const { return cw?_u.colVar():_u.rowVar(); }

    index_t cardinality_impl() const { return _u.cardinality_impl(); }

    void print(std::ostream &os) const { os<<"("; _u.print(os); os <<")\u1D40"; }
private:
/*
  template<class U> EIGEN_STRONG_INLINE MatExprType
  eval_impl(const U k, typename util::enable_if<1==ColBlocks,U>::type* = nullptr)
  { return _u.eval(k).blockTranspose(_u.cols()/_u.rows()); }

  template<class U> EIGEN_STRONG_INLINE MatExprType
  eval_impl(const U k, typename util::enable_if<0==ColBlocks,U>::type* = nullptr)
  { return _u.eval(k).transpose(); }
*/
};

}// namespace expr
}// namespace gismo