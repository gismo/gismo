/** @file gsVector_expr.h

    @brief Defines vector expressions

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

    #define GISMO_EXPR_VECTOR_EXPRESSION(name, mname, isSv)                 \
    template<class E> class name##_##expr  : public _expr<name##_##expr<E> > { \
        typename E::Nested_t _u;                                        \
    public:                                                             \
    typedef typename E::Scalar Scalar;                                  \
    enum {Space= E::Space, ScalarValued= isSv, ColBlocks= E::ColBlocks}; \
    name##_##expr(_expr<E> const& u) : _u(u) { }                        \
    mutable Temporary_t tmp;                                            \
    const Temporary_t & eval(const index_t k) const {                   \
        tmp = _u.eval(k).mname(); return tmp; }                         \
    index_t rows() const { return isSv ? 0 : _u.rows(); }               \
    index_t cols() const { return isSv ? 0 : _u.cols(); }               \
    void parse(gsExprHelper<Scalar> & evList) const { _u.parse(evList); } \
    const gsFeSpace<Scalar> & rowVar() const {return gsNullExpr<Scalar>::get();} \
    const gsFeSpace<Scalar> & colVar() const {return gsNullExpr<Scalar>::get();} \
    void print(std::ostream &os) const                                  \
        { os << #name <<"("; _u.print(os); os <<")"; }                  \
    };

/// Eucledian Norm
GISMO_EXPR_VECTOR_EXPRESSION(norm,norm,1);
/// Squared Eucledian Norm
GISMO_EXPR_VECTOR_EXPRESSION(sqNorm,squaredNorm,1);
/// Normalization of a vector to unit measure
GISMO_EXPR_VECTOR_EXPRESSION(normalized,normalized,0);
/// Inverse of a matrix expression
GISMO_EXPR_VECTOR_EXPRESSION(inv,cramerInverse,0);
// GISMO_EXPR_VECTOR_EXPRESSION(cwSqr,array().square,0)
// GISMO_EXPR_VECTOR_EXPRESSION(sum,array().sum,1)
// GISMO_EXPR_VECTOR_EXPRESSION(sqrt,array().sqrt,0)
//GISMO_EXPR_VECTOR_EXPRESSION(abs,array().abs,0)

//Determinant
GISMO_EXPR_VECTOR_EXPRESSION(det,determinant,1);

//GISMO_EXPR_VECTOR_EXPRESSION(replicate,replicate,0);

#undef GISMO_EXPR_VECTOR_EXPRESSION

}// namespace expr
}// namespace gismo