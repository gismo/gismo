/** @file ComponentExpression.h

    @brief Expression for accessing components of vector/matrix expressions

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M.Verhelst
*/

#pragma once

namespace gismo
{
namespace Expr
{

// Forward declaration
template <typename E, size_t _Order> class ComponentExpression;

// ExpressionTraits for ComponentExpression - always returns a scalar (Order 0)
template <typename E, size_t _Order>
struct ExpressionTraits<ComponentExpression<E, _Order>>
{
    typedef E ExprType;
    
    typedef typename ExpressionTraits<E>::Scalar Scalar;
    static constexpr size_t Order = 0; // Component access always returns scalar
    static constexpr SpaceType Space = ExpressionTraits<E>::Space;
    static constexpr size_t Deriv = ExpressionTraits<E>::Deriv;
    static constexpr bool IsConstant = ExpressionTraits<E>::IsConstant;
};

/**
 * \brief Expression for accessing a component of a vector (_Order=1) or matrix (_Order=2)
 * \tparam E The expression type
 * \tparam _Order The order of the parent expression (1 for vector, 2 for matrix)
 */
template <typename E, size_t _Order>
class ComponentExpression : public BaseExpression<ComponentExpression<E, _Order>>
{
    using Base = BaseExpression<ComponentExpression<E, _Order>>;
    typedef typename Base::Scalar T;

public:
    using Base::Order;
    using Base::Space;
    using Base::Deriv;
    using Base::IsConstant;
    using Base::Scalar;

    typedef typename ExpressionTraits<ComponentExpression<E, _Order>>::ExprType ExprType;

    const ExprType& expr() const { return static_cast<const ExprType&>(expr_); }

    // For vector component access [i]
    ComponentExpression(const ExprType& expr, index_t i)
    : BaseExpression<ComponentExpression<E, _Order>>(), expr_(expr), i_(i), j_(-1)
    {
        static_assert(_Order == 1, "Single index constructor only for Order 1 expressions");
    }

    // For matrix component access (i,j)
    ComponentExpression(const ExprType& expr, index_t i, index_t j)
    : BaseExpression<ComponentExpression<E, _Order>>(), expr_(expr), i_(i), j_(j)
    {
        static_assert(_Order == 2, "Double index constructor only for Order 2 expressions");
    }

    void parse(gismo::ExpressionHelper<T> & helper) const
    {
        expr_.parse(helper);
        expr_.setDerivative(Deriv);
    }

    void print(std::ostream & os) const
    {
        if constexpr (_Order == 1)
            os << "(" << expr_ << ")[" << i_ << "]";
        else
            os << "(" << expr_ << ")(" << i_ << "," << j_ << ")";
    }

    std::array<size_t, 0> sizes() const 
    { 
        return std::array<size_t, 0>(); // Scalar result
    }

    size_t domainDim() const { return expr_.domainDim(); }

    // Default implementation: extract component from evaluated expression
    ExpressionResult<T> eval(const index_t k) const
    {
        auto result = expr_.eval(k);
        
        if constexpr (_Order == 1)
        {
            // Vector component
            GISMO_ASSERT(i_ >= 0 && i_ < result.rows(), "Component index out of bounds");
            return result(i_, 0);
        }
        else // _Order == 2
        {
            // Matrix component
            GISMO_ASSERT(i_ >= 0 && i_ < result.rows(), "Row index out of bounds");
            GISMO_ASSERT(j_ >= 0 && j_ < result.cols(), "Col index out of bounds");
            return result(i_, j_);
        }
    }

    // Forward test/trial to parent expression
    auto test() const -> decltype(std::declval<const ComponentExpression&>().expr().test())
    {
        return expr().test();
    }

    auto trial() const -> decltype(std::declval<const ComponentExpression&>().expr().trial())
    {
        return expr().trial();
    }

protected:
    const ExprType expr_;
    const index_t i_;
    const index_t j_; // Only used for Order 2
};

} // namespace Expr
} // namespace gismo
