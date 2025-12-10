/** @file abs_expr.h

    @brief Defines the abs expression

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
    \brief Expression for the absolute value of an expression
    \ingroup Expressions
    \tparam E The expression type
*/
template<class E>
class abs_expr  : public _expr<abs_expr<E> >
{
    typename E::Nested_t _u;

public:
    typedef typename E::Scalar Scalar;
    explicit abs_expr(_expr<E> const& u) : _u(u) { }

public:
    enum {Space= 0, ScalarValued= 1, ColBlocks= 0};

    Scalar eval(const index_t k) const { return abs_expr::eval_impl(_u,k); }

    index_t rows() const { return _u.rows(); }
    index_t cols() const { return _u.cols(); }
    void parse(gsExprHelper<Scalar> & evList) const
    { _u.parse(evList); }

    static bool isScalar() { return true; }

    const gsFeSpace<Scalar> & rowVar() const {return gsNullExpr<Scalar>::get();}
    const gsFeSpace<Scalar> & colVar() const {return gsNullExpr<Scalar>::get();}

    void print(std::ostream &os) const { _u.print(os); }

    // Math functions. eg
    // sqrt_mexpr<T> sqrt() { return sqrt_mexpr<T>(*this); }
private:
    template<class U> static inline
    typename util::enable_if<U::ScalarValued,Scalar>::type
    eval_impl(const U & u, const index_t k) {return math::abs(u.eval(k)); }
    template<class U> static inline
    typename util::enable_if<!U::ScalarValued,gsMatrix<Scalar> >::type
    eval_impl(const U & u, const index_t k) { return u.eval(k).cwiseAbs(); }
};

/**
 * @brief Returns the absolute value of an expression
 * @param u The expression
 * @ingroup Expressions
 */
template<class E> EIGEN_STRONG_INLINE
abs_expr<E> abs(const E & u) { return abs_expr<E>(u); }

}// namespace expr
}// namespace gismo