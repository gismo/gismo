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
// Space: Flag that indicates whether the expression is a Space
template <typename E>
class BaseExpression
{
protected:
    BaseExpression() : Deriv_(ExpressionTraits<E>::Deriv) {};
    BaseExpression(const BaseExpression&) : Deriv_(ExpressionTraits<E>::Deriv) {};

    mutable short_t Deriv_; // Used to store the required derivative order

public:

    typedef typename ExpressionTraits<E>::Scalar Scalar;
    static constexpr size_t Order = ExpressionTraits<E>::Order;
    static constexpr size_t Space = ExpressionTraits<E>::Space;
    static constexpr size_t Deriv = ExpressionTraits<E>::Deriv;
    static constexpr bool IsConstant = ExpressionTraits<E>::IsConstant;

    // Getters of static members
    static size_t order() {return Order;};
    static size_t space() {return Space;};
    static size_t deriv() {return Deriv;};
    static size_t isConstant() {return IsConstant;};

    size_t domainDim() const { return static_cast<const E&>(*this).domainDim(); }

    const std::array<size_t,Order> & sizes() const { return static_cast<const E&>(*this).sizes(); }

    void print(std::ostream & os) const { static_cast<const E&>(*this).print(os); }

    ExpressionValue<Scalar> eval(const index_t k) const { return static_cast<const E&>(*this).eval(k); }

    bool isZero() const { return false; }

    void parse(gismo::ExpressionHelper<Scalar> &) const
    {
        GISMO_NO_IMPLEMENTATION;
        // static_cast<E const&>(*this).parse(helper);
    }

    const SpaceObject<Scalar, Space, Order> & test () const { return Space==SpaceType::Both ? static_cast<E const&>(*this).test() : NullObject<Scalar,Space,Order>::get(); }
    const SpaceObject<Scalar, Space, Order> & trial() const { return Space==SpaceType::Both ? static_cast<E const&>(*this).trial() : NullObject<Scalar,Space,Order>::get(); }

public:

    void setDerivative(short_t d) const
    {
        Deriv_ = math::max(Deriv_, d);
    }

    short_t getDerivative() const { return Deriv_; }


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

template <typename E, typename SpaceObj>
VariationObject<BaseExpression<E>, SpaceObj> variation(const BaseExpression<E> & sol,const SpaceObj & space)
{
    GISMO_ERROR("Not implemented.");
}

}//namespace Expr
}//namespace gismo