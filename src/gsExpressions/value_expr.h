/** @file value_expr.h

    @brief Defines the value expression

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
 * @brief Expression adaptor for a scalar-valued expression
 * @ingroup Expressions
 * @tparam E The type of the expression
 */
template<class E>
class value_expr  : public _expr<value_expr<E> >
{
    typename E::Nested_t _u;

public:
    typedef typename E::Scalar Scalar;
    value_expr(_expr<E> const& u) : _u(u)
    {
        // rows/cols not known at construction time
        //GISMO_ASSERT(u.rows()*u.cols()<=1, "Expression\n"<<u<<"is not a scalar.");
    }

public:
    enum {Space= 0, ScalarValued= 1, ColBlocks= 0};

    Scalar eval(const index_t k) const { return eval_impl(_u,k); }

    // enables expr.val().val()
    inline value_expr<E> val() const { return *this; }
    index_t rows() const { return 0; }
    index_t cols() const { return 0; }
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
    eval_impl(const U & u, const index_t k) { return u.eval(k); }

    template<class U> static inline
    typename util::enable_if<!U::ScalarValued,Scalar>::type
    eval_impl(const U & u, const index_t k) { return u.eval(k).value(); }
};

}// namespace expr
}// namespace gismo