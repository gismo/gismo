/** @file sign_expr.h

    @brief Defines the sign expression

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
 * @brief Expression for the sign of a finite element variable
 * @ingroup Expressions
 * @tparam E The type of the expression
 */
template<class E>
class sign_expr : public _expr<sign_expr<E> >
{
    typename E::Nested_t _u;
    typename E::Scalar _tol;
public:
    typedef typename E::Scalar Scalar;
    enum {ScalarValued = 1, Space = E::Space, ColBlocks= 0};

    sign_expr(_expr<E> const& u, Scalar tolerance = 0.0) : _u(u),_tol(tolerance){
        GISMO_ASSERT( _tol >= 0, "Tolerance for sign_expr should be a positive number.");
    }

    Scalar eval(const index_t k) const
    {
        const Scalar v = _u.val().eval(k);
        return ( v>_tol ? Scalar(1) : ( v<-_tol ? Scalar(-1) : Scalar(0) ) );
    }

    static index_t rows() { return 0; }
    static index_t cols() { return 0; }

    void parse(gsExprHelper<Scalar> & el) const
    { _u.parse(el); }

    static bool isScalar() { return true; }

    const gsFeSpace<Scalar> & rowVar() const {return gsNullExpr<Scalar>::get();}
    const gsFeSpace<Scalar> & colVar() const {return gsNullExpr<Scalar>::get();}

    void print(std::ostream &os) const { os<<"sgn("; _u.print(os); os <<")"; }
};

}// namespace expr
}// namespace gismo