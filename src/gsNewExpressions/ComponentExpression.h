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
        // Set derivative order BEFORE parsing so that parse() knows what flags to set
        expr_.setDerivative(Deriv);
        expr_.parse(helper);
    }

    void print(std::ostream & os) const
    {
        print_impl(os, SizeTag<_Order>());
    }

private:
    void print_impl(std::ostream & os, SizeTag<1>) const
    {
        os << "(" << expr_ << ")[" << i_ << "]";
    }
    void print_impl(std::ostream & os, SizeTag<2>) const
    {
        os << "(" << expr_ << ")(" << i_ << "," << j_ << ")";
    }

public:
    std::array<size_t, 0> sizes() const 
    { 
        return std::array<size_t, 0>(); // Scalar result
    }

    size_t domainDim() const { return expr_.domainDim(); }

    // Default implementation: extract component from evaluated expression
    ExpressionResult<T> eval(const index_t k) const
    {
        return eval_impl(k, SizeTag<_Order>());
    }

private:
    ExpressionResult<T> eval_impl(const index_t k, SizeTag<1>) const
    {
        auto result = expr_.eval(k);
        // Vector component - extract single value
        const auto& mat = result();
        GISMO_ASSERT(i_ >= 0 && i_ < mat.rows(), "Component index out of bounds");
        gsMatrix<T> component(1, 1);
        component(0, 0) = mat(i_, 0);
        return ExpressionResult<T>(component);
    }
    ExpressionResult<T> eval_impl(const index_t k, SizeTag<2>) const
    {
        auto result = expr_.eval(k);
        // Matrix component - extract single value
        const auto& mat = result();
        GISMO_ASSERT(i_ >= 0 && i_ < mat.rows(), "Row index out of bounds");
        GISMO_ASSERT(j_ >= 0 && j_ < mat.cols(), "Col index out of bounds");
        gsMatrix<T> component(1, 1);
        component(0, 0) = mat(i_, j_);
        return ExpressionResult<T>(component);
    }

public:
    // Forward test/trial to parent expression
    auto test() const -> decltype(std::declval<const ComponentExpression&>().expr().test())
    {
        return expr().test();
    }

    auto trial() const -> decltype(std::declval<const ComponentExpression&>().expr().trial())
    {
        return expr().trial();
    }
    
    // Forward data() to parent expression
    auto data() const -> decltype(std::declval<const ComponentExpression&>().expr().data())
    {
        return expr().data();
    }
    
    // Get component index (for use by derivative expressions)
    index_t component() const { return i_; }

protected:
    const ExprType expr_;
    const index_t i_;
    const index_t j_; // Only used for Order 2
};

} // namespace Expr
} // namespace gismo
