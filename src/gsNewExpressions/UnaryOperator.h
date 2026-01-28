/** @file UnaryOperator.h

    @brief Base class for unary operators

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

// ExpressionTraits for UnaryOperator - forwards template parameters
template <typename Operator>
struct ExpressionTraits<UnaryOperator<Operator>>
{
    typedef typename ExpressionTraits<Operator>::ExprType ExprType;

    typedef typename ExpressionTraits<Operator>::Scalar Scalar;
    static constexpr size_t Order = ExpressionTraits<Operator>::Order;
    static constexpr SpaceType Space = ExpressionTraits<Operator>::Space;
    static constexpr size_t Deriv = ExpressionTraits<Operator>::Deriv;
    static constexpr bool IsConstant = ExpressionTraits<Operator>::IsConstant;
};

/**
 * \brief Base class for unary operators
 * \tparam Expr The operand expression type
 */
template <typename Operator>
class UnaryOperator : public BaseExpression<Operator>
{
    using Base = BaseExpression<Operator>;
    typedef typename Base::Scalar T;

public:
    // Expose the static traits publicly so they can be accessed as LhsExpr::Order etc.
    using Base::Order;
    using Base::Space;
    using Base::Deriv;
    using Base::IsConstant;
    using Base::Scalar;

    typedef typename ExpressionTraits<Operator>::ExprType ExprType;

    const ExprType& expr() const { return static_cast<const ExprType&>(expr_); }

    void parse(gismo::ExpressionHelper<T> & helper) const
    {
        expr_.parse(helper);
        expr_.setDerivative(Deriv);
    }

    void print(std::ostream & os) const
    {
        static_cast<const Operator&>(*this).print(os);
    }

    // Note: sizes() and domainDim() are NOT defined here - derived classes must implement them
    // This allows each derived class to define its own order and sizes behavior

    template<class Op = Operator>
    auto test() const -> decltype(std::declval<const UnaryOperator&>().expr().test()) 
    { 
        return expr().test();  
    } 
    template<class Op = Operator>
    auto trial() const -> decltype(std::declval<const UnaryOperator&>().expr().trial()) 
    { 
        return expr().trial(); 
    }


protected:
    UnaryOperator(const ExprType& expr)
    : BaseExpression<Operator>(), expr_(expr)
    {
    }

    const ExprType expr_;
};

}//namespace Expr
}//namespace gismo