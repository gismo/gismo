/** @file BaseExpression.h
    @brief Base class for all expressions in the new expression system

    This file defines the CRTP base class for all expressions. All expression
    types (variables, operators, spaces, etc.) inherit from BaseExpression.

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
    // Forward declarations
    template <class T, enum SpaceType _Space, size_t _Order> class NullObject;
    template <class T, enum SpaceType _Space, size_t _Order> class SpaceObject;

/**
 * @brief Base class for all expressions using CRTP (Curiously Recurring Template Pattern)
 * 
 * This is the foundation of the expression system. All expressions inherit from this class
 * and provide their own implementation of key methods like parse(), eval(), etc.
 * 
 * @tparam E The derived expression type (CRTP parameter)
 * 
 * Key concepts:
 * - Order: The polynomial order of the expression (0 for constants, 1 for linear, etc.)
 * - Space: Indicates if expression belongs to a test/trial space (Test, Trial, or None)
 * - Deriv_: Tracks required derivative order for evaluation
 * 
 * The expression system uses a two-phase evaluation:
 * 1. parse() - Set up required data flags and propagate derivative requirements
 * 2. eval() - Evaluate the expression at quadrature points
 */
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
    static constexpr SpaceType Space = ExpressionTraits<E>::Space;
    static constexpr size_t Deriv = ExpressionTraits<E>::Deriv;
    static constexpr bool IsConstant = ExpressionTraits<E>::IsConstant;

    // Getters of static members
    static size_t order() {return Order;};
    static SpaceType space() {return Space;};
    static size_t deriv() {return Deriv;};
    static size_t isConstant() {return IsConstant;};

    size_t domainDim() const { return static_cast<const E&>(*this).domainDim(); }

    const std::array<size_t,Order> & sizes() const { return static_cast<const E&>(*this).sizes(); }

    void print(std::ostream & os) const { static_cast<const E&>(*this).print(os); }

    ExpressionResult<Scalar> eval(const index_t k) const { return static_cast<const E&>(*this).eval(k); }

    bool isZero() const { return false; }

    void parse(gismo::ExpressionHelper<Scalar> & helper) const
    {
        static_cast<E const&>(*this).parse(helper);
    }

    // Note: test() and trial() are NOT defined in BaseExpression
    // Each expression class implements them with appropriate return types
    // This avoids incomplete type issues with decltype and allows flexible return types

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

    // Component access for vector expressions (Order 1)
    template<typename Expr = E>
    typename std::enable_if<Expr::Order == 1, ComponentExpression<E, 1>>::type
    operator[](index_t i) const
    {
        return ComponentExpression<E, 1>(static_cast<const E&>(*this), i);
    }

    // Component access for matrix expressions (Order 2)
    template<typename Expr = E>
    typename std::enable_if<Expr::Order == 2, ComponentExpression<E, 2>>::type
    operator()(index_t i, index_t j) const
    {
        return ComponentExpression<E, 2>(static_cast<const E&>(*this), i, j);
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