/** @file TransposeExpression.h

@brief

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

template <typename E>
struct ExpressionTraits<TransposeExpression<E>>
{
    typedef E ExprType; // Needed for UnaryOperator

    typedef typename ExpressionTraits<E>::Scalar Scalar;
    static constexpr size_t Order = ExpressionTraits<E>::Order;
    static constexpr SpaceType Space = ExpressionTraits<E>::Space;
    static constexpr size_t Deriv = ExpressionTraits<E>::Deriv;
    static constexpr bool IsConstant = ExpressionTraits<E>::IsConstant;
};

template<typename E>
class TransposeExpression : public UnaryOperator< TransposeExpression<E> >
{
    static_assert(ExpressionTraits<E>::Order == 1 || ExpressionTraits<E>::Order == 2,
                    "TransposeExpression: Only vector (order 1) or matrix (order 2) expressions can be transposed.");

    using Base = UnaryOperator<TransposeExpression<E>>;

protected:
    mutable std::array<size_t, Base::Order> sizes_;

public:
    TransposeExpression(const E& expr)
    :
    Base(expr)
    {
        initializeSizes(expr, SizeTag<Base::Order>());
    }

private:
    void initializeSizes(const E& expr, SizeTag<1>)
    {
        // Vector transpose: sizes remain the same
        sizes_ = expr.sizes();
    }

    void initializeSizes(const E& expr, SizeTag<2>)
    {
        // Matrix transpose: swap dimensions
        auto orig_sizes = expr.sizes();
        sizes_[0] = orig_sizes[1];
        sizes_[1] = orig_sizes[0];
    }

public:

    const std::array<size_t, Base::Order> & sizes() const
    {
        return sizes_;
    }

    size_t domainDim() const
    {
        return expr_.domainDim();
    }

    ExpressionResult<typename Base::Scalar> eval(const index_t k) const
    {
        ExpressionResult<typename Base::Scalar> val = expr_.eval(k);
        ExpressionResult<typename Base::Scalar> result(val.rowCardinality(), val.colCardinality());
        
        for (index_t i = 0; i < result.rowCardinality(); ++i)
        {
            for (index_t j = 0; j < result.colCardinality(); ++j)
            {
                result(i, j) = val(i, j).transpose();
            }
        }
        
        return result;
    }

    void parse(ExpressionHelper<typename Base::Scalar> & helper) const
    {
        // Set derivative order first so child expression knows what to request
        expr_.setDerivative(Base::Deriv);
        expr_.parse(helper);
    }

    void print(std::ostream & os) const
    {
        os<<"("<<expr_<<")\u1D40";
    }

protected:
    using Base::expr_;

};

// Factory function for easier creation
template <typename E>
TransposeExpression<E> transpose(const E& expr)
{
    return TransposeExpression<E>(expr);
}

}//namespace Expr
}//namespace gismo