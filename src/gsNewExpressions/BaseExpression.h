/** @file BaseObject.h
@brief

This file is part of the G+Smo library.

This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this
file, You can obtain one at http://mozilla.org/MPL/2.0/.

Author(s): H.M.Verhelst
*/

#pragma once

#include <gsMatrix/gsMatrix.h>

namespace gismo
{
namespace Expr
{
// IsConstant: Flag that indicates if the expression is constant, e.g., its derivatives are zero
// Space: Flag that indicates whether the expression is a space
template <typename E>
class BaseExpression
{
protected:
    BaseExpression() : deriv_(ExpressionTraits<E>::deriv) {};
    BaseExpression(const BaseExpression&) : deriv_(ExpressionTraits<E>::deriv) {};

    mutable short_t deriv_; // Used to store the required derivative order

public:

    typedef typename ExpressionTraits<E>::Scalar Scalar;
    static constexpr size_t order = ExpressionTraits<E>::order;
    static constexpr size_t space = ExpressionTraits<E>::space;
    static constexpr size_t deriv = ExpressionTraits<E>::deriv;
    static constexpr bool isConstant = ExpressionTraits<E>::isConstant;

    size_t domainDim() const { return static_cast<const E&>(*this).domainDim(); }

    const std::array<size_t,order> & sizes() const { return static_cast<const E&>(*this).sizes(); }

    void print(std::ostream & os) const { static_cast<const E&>(*this).print(os); }

    gsMatrix<Scalar> eval(const index_t k) const { return static_cast<const E&>(*this).eval(k); }

    void parse(gismo::ExpressionHelper<Scalar> &) const
    {
        GISMO_NO_IMPLEMENTATION;
        // static_cast<E const&>(*this).parse(helper);
    }

    const SpaceObject<Scalar, space, order> & rowVar() const
    {
        GISMO_NO_IMPLEMENTATION;
    //     return static_cast<E const&>(*this).rowVar();
    }

    const SpaceObject<Scalar, space, order> & colVar() const
    {
        GISMO_NO_IMPLEMENTATION;
    //     return static_cast<E const&>(*this).colVar();
    }

public:

    void setDerivative(short_t d) const
    {
        deriv_ = math::max(deriv_, d);
    }

    short_t getDerivative() const { return deriv_; }


    TransposeExpression<E> tr() const
    {
        return TransposeExpression<E>(static_cast<E const&>(*this));
    }

    ArrayExpression<E> array() const
    {
        return ArrayExpression<E>(static_cast<E const&>(*this));
    }


    // Overload conversions
    operator E&()             { return static_cast<      E&>(*this); }
    operator E const&() const { return static_cast<const E&>(*this); }

    E const & derived() const { return static_cast<const E&>(*this); }
};

/// Stream operator for expressions
template <typename E>
std::ostream &operator<<(std::ostream &os, const BaseExpression<E> & b)
{b.print(os); return os; }

}//namespace Expr
}//namespace gismo