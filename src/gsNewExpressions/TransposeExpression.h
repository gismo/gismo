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
    static constexpr size_t order = ExpressionTraits<E>::order;
    static constexpr size_t space = ExpressionTraits<E>::space;
    static constexpr size_t deriv = ExpressionTraits<E>::deriv;
    static constexpr bool isConstant = ExpressionTraits<E>::isConstant;
};

template<typename E>
class TransposeExpression : public UnaryOperator< TransposeExpression<E> >
{
    static_assert(ExpressionTraits<E>::order == 1 || ExpressionTraits<E>::order == 2,
                    "TransposeExpression: Only vector (order 1) or matrix (order 2) expressions can be transposed.");

    using Base = UnaryOperator<TransposeExpression<E>>;

protected:
    mutable std::array<size_t, Base::order> sizes_;

public:
    TransposeExpression(const E& expr)
    :
    Base(expr)
    {
        initializeSizes(expr, std::integral_constant<size_t, Base::order>{});
    }

private:
    void initializeSizes(const E& expr, std::integral_constant<size_t, 1>)
    {
        // Vector transpose: sizes remain the same
        sizes_ = expr.sizes();
    }

    void initializeSizes(const E& expr, std::integral_constant<size_t, 2>)
    {
        // Matrix transpose: swap dimensions
        auto orig_sizes = expr.sizes();
        sizes_[0] = orig_sizes[1];
        sizes_[1] = orig_sizes[0];
    }

public:

    const std::array<size_t, Base::order> & sizes() const
    {
        return sizes_;
    }

    size_t domainDim() const
    {
        return expr_.domainDim();
    }

    gsMatrix<typename Base::Scalar> eval(const index_t k) const
    {
        return expr_.eval(k).transpose();
    }

    void parse(ExpressionHelper<typename Base::Scalar> & helper) const
    {
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