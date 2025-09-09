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
    static constexpr size_t order = ExpressionTraits<Operator>::order;
    static constexpr size_t space = ExpressionTraits<Operator>::space;
    static constexpr size_t deriv = ExpressionTraits<Operator>::deriv;
    static constexpr bool isConstant = ExpressionTraits<Operator>::isConstant;
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
    using Base::deriv;
public:
    typedef typename ExpressionTraits<Operator>::ExprType ExprType;

    const ExprType& expr() const { return static_cast<const ExprType&>(expr_); }

    void parse(gismo::ExpressionHelper<T> & helper) const
    {
        expr_.parse(helper);
        expr_.setDerivative(deriv);
    }

    void print(std::ostream & os) const
    {
        static_cast<const Operator&>(*this).print(os);
    }

    // Note: sizes() and domainDim() are NOT defined here - derived classes must implement them
    // This allows each derived class to define its own order and sizes behavior

protected:
    UnaryOperator(const ExprType& expr)
    : BaseExpression<Operator>(), expr_(expr)
    {
    }

    const ExprType expr_;
};

}//namespace Expr
}//namespace gismo