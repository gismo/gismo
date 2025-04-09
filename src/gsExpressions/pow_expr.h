/** @file pow_expr.h

    @brief Defines the pow expression

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
/**
 * @brief Represents an expression for computing the power of a base expression raised to a given exponent.
 * @tparam E The type of the base expression.
 * @ingroup Expressions
 */
{
template<class E>
class pow_expr : public _expr<pow_expr<E> >
{
    typename E::Nested_t _u;

public:
    typedef typename E::Scalar Scalar;
    enum {ScalarValued = 1, Space = E::Space, ColBlocks= E::ColBlocks};

    Scalar _q;// power

    pow_expr(_expr<E> const& u, Scalar q) : _u(u), _q(q) { }

    Scalar eval(const index_t k) const
    {
        const Scalar v = _u.val().eval(k);
        return math::pow(v,_q);
    }

    static index_t rows() { return 0; }
    static index_t cols() { return 0; }

    void parse(gsExprHelper<Scalar> & el) const
    { _u.parse(el); }

    static bool isScalar() { return true; }

    const gsFeSpace<Scalar> & rowVar() const {return gsNullExpr<Scalar>::get();}
    const gsFeSpace<Scalar> & colVar() const {return gsNullExpr<Scalar>::get();}

    void print(std::ostream &os) const { os<<"pow("; _u.print(os); os <<")"; }
};

/**
 * @brief Creates a power expression by raising the given expression to the specified power.
 *
 * @tparam E The type of the expression being raised to a power.
 * @param u The input expression to be raised to the power.
 * @param q The exponent to which the expression is raised.
 * @return A pow_expr object representing the result of raising the input expression to the power q.
 */
template<class E> pow_expr<E>
pow(_expr<E> const& u, real_t q) { return pow_expr<E>(u,q); }

}// namespace expr
}// namespace gismo