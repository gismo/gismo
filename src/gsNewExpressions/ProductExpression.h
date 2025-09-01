/** @file ProductExpression.h

    @brief Product expression class

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

// --- ExpressionTraits specializations for ProductExpression ---
// Primary template: Catches all unsupported combinations with a compile-time error
template <typename LhsExpr, typename RhsExpr, size_t LhsOrder, size_t RhsOrder, size_t LhsSpace, size_t RhsSpace>
struct ExpressionTraits<ProductExpression<LhsExpr, RhsExpr, LhsOrder, RhsOrder, LhsSpace, RhsSpace>>
{
    typedef LhsExpr LhsType;
    typedef RhsExpr RhsType;

    typedef typename ExpressionTraits<LhsExpr>::Scalar Scalar;
    static constexpr size_t order = 0; // Default/fallback, should not be used
    static constexpr size_t space = Space::None; // TODO: Define appropriate space logic
    static constexpr size_t deriv = ExpressionTraits<RhsExpr>::deriv; //TODO // Conservative approach
    static constexpr bool isConstant = ExpressionTraits<LhsExpr>::isConstant && ExpressionTraits<RhsExpr>::isConstant;
};

template <typename LhsExpr, typename RhsExpr, size_t LhsOrder, size_t RhsOrder, size_t LhsSpace, size_t RhsSpace>
class ProductExpression
{
    static_assert(std::is_same<LhsExpr, void>::value,
                  "ProductExpression: Unsupported tensor order combination for product.");
};

// --- Partial specialization 1: Scalar multiplication
template <typename LhsExpr, typename RhsExpr, size_t RhsOrder, size_t LhsSpace, size_t RhsSpace>
struct ExpressionTraits<ProductExpression<LhsExpr, RhsExpr, 0, RhsOrder, LhsSpace, RhsSpace>>
{
    typedef LhsExpr LhsType;
    typedef RhsExpr RhsType;

    typedef typename ExpressionTraits<LhsExpr>::Scalar Scalar;
    static constexpr size_t order = RhsOrder;
    static constexpr size_t space = ExpressionTraits<RhsExpr>::space; // TODO: Define appropriate space logic
    static constexpr size_t deriv = ExpressionTraits<RhsExpr>::deriv; // Conservative approach
    static constexpr bool isConstant = ExpressionTraits<LhsExpr>::isConstant && ExpressionTraits<RhsExpr>::isConstant;
};

template <typename LhsExpr, typename RhsExpr, size_t RhsOrder, size_t LhsSpace, size_t RhsSpace>
class ProductExpression<LhsExpr, RhsExpr, 0, RhsOrder, LhsSpace, RhsSpace>
 : public BinaryOperator<ProductExpression<LhsExpr, RhsExpr, 0, RhsOrder, LhsSpace, RhsSpace>>
{
    using Base = BinaryOperator<ProductExpression<LhsExpr, RhsExpr, 0, RhsOrder, LhsSpace, RhsSpace>>;

public:
    ProductExpression(const LhsExpr& lhs, const RhsExpr& rhs)
        : Base(lhs, rhs)
    {
    }

    const std::array<size_t, Base::order> & sizes() const
    {
        return this->rhs_expr_.sizes();
    }

    size_t domainDim() const
    {
        return this->rhs_expr_.domainDim();
    }

    gsMatrix<typename Base::Scalar> eval(const index_t k) const
    {
        gsMatrix<typename Base::Scalar> lhs_val = this->lhs_expr_.eval(k);
        gsMatrix<typename Base::Scalar> rhs_val = this->rhs_expr_.eval(k);
        rhs_val.array() *= lhs_val.value(); // Element-wise multiplication
        return rhs_val; // Return the modified rhs_val
    }


    void print(std::ostream & os) const override
    {
        os<<this->lhs_expr_<<"*"<<this->rhs_expr_;
    }

private:
    using Base::lhs_expr_;
    using Base::rhs_expr_;
};

// --- Partial specialization 2: Multiplication of equal order (matrix only)
template <typename LhsExpr, typename RhsExpr, size_t LhsSpace, size_t RhsSpace>
struct ExpressionTraits<ProductExpression<LhsExpr, RhsExpr, 2, 2, LhsSpace, RhsSpace>>
{
    typedef LhsExpr LhsType;
    typedef RhsExpr RhsType;

    typedef typename ExpressionTraits<LhsExpr>::Scalar Scalar;
    static constexpr size_t order = 2;
    static constexpr size_t space = ExpressionTraits<RhsExpr>::space; // TODO: Define appropriate space logic
    static constexpr size_t deriv = ExpressionTraits<RhsExpr>::deriv; // Conservative approach
    static constexpr bool isConstant = ExpressionTraits<LhsExpr>::isConstant && ExpressionTraits<RhsExpr>::isConstant;
};

template <typename LhsExpr, typename RhsExpr, size_t LhsSpace, size_t RhsSpace>
class ProductExpression<LhsExpr, RhsExpr, 2, 2, LhsSpace, RhsSpace>
 : public BinaryOperator<ProductExpression<LhsExpr, RhsExpr, 2, 2, LhsSpace, RhsSpace>>
{
    using Base = BinaryOperator<ProductExpression<LhsExpr, RhsExpr, 2, 2, LhsSpace, RhsSpace>>;
    using Base::order;

public:
    using Scalar = typename Base::Scalar;

protected:

    std::array<size_t,order> sizes_;

public:
    ProductExpression(const LhsExpr& lhs, const RhsExpr& rhs)
        : Base(lhs, rhs)
    {
        GISMO_ENSURE(lhs.sizes()[1]==rhs.sizes()[0],
                     "Size mismatch in * operator.\n"<<
                     "lhs ("<<lhs.sizes()[0]<<" x "<<lhs.sizes()[1]<<") = "<<lhs<<"\n"<<
                     "rhs ("<<rhs.sizes()[0]<<" x "<<rhs.sizes()[1]<<") = "<<rhs<<")\n");
        sizes_[0] = lhs.sizes()[0];
        sizes_[1] = rhs.sizes()[1];
    }

    const std::array<size_t, order> & sizes() const { return sizes_; }

    size_t domainDim() const
    {
        gsWarn<<"Correct?\n";
        return this->lhs_expr_.domainDim();
    }

    gsMatrix<Scalar> eval(const index_t k) const
    {
        gsMatrix<Scalar> lhs_val = this->lhs_expr_.eval(k);
        gsMatrix<Scalar> rhs_val = this->rhs_expr_.eval(k);
        return lhs_val * rhs_val;
    }

    void print(std::ostream & os) const override
    {
        os<<this->lhs_expr_<<"*"<<this->rhs_expr_;
    }

private:
    using Base::lhs_expr_;
    using Base::rhs_expr_;
};

// --- ExpressionTraits for matrix-vector product (2,1) -> vector (1)
template <typename LhsExpr, typename RhsExpr, size_t LhsSpace, size_t RhsSpace>
struct ExpressionTraits<ProductExpression<LhsExpr, RhsExpr, 2, 1, LhsSpace, RhsSpace>>
{
    using Scalar = typename ExpressionTraits<LhsExpr>::Scalar;
    static constexpr size_t order = 1;
    static constexpr size_t space = ExpressionTraits<RhsExpr>::space; // TODO: Define appropriate space logic
    static constexpr size_t deriv = ExpressionTraits<RhsExpr>::deriv; //TODO
    static constexpr bool isConstant = ExpressionTraits<LhsExpr>::isConstant && ExpressionTraits<RhsExpr>::isConstant;
};

// --- Partial specialization 3: Multiplication of matrix with vector
template <typename LhsExpr, typename RhsExpr, size_t LhsSpace, size_t RhsSpace>
class ProductExpression<LhsExpr, RhsExpr, 2, 1, LhsSpace, RhsSpace>
 : public BinaryOperator<ProductExpression<LhsExpr, RhsExpr, 2, 1, LhsSpace, RhsSpace>>
{
    using Base = BinaryOperator<ProductExpression<LhsExpr, RhsExpr, 2, 1, LhsSpace, RhsSpace>>;
    using Base::order;

public:
    using Scalar = typename Base::Scalar;

protected:

    std::array<size_t,order> sizes_;

public:
    ProductExpression(const LhsExpr& lhs, const RhsExpr& rhs)
        : Base(lhs, rhs)
    {
        GISMO_ENSURE(lhs.sizes()[1]==rhs.sizes()[0],
                     "Size mismatch in * operator.\n"<<
                     "lhs ("<<lhs.sizes()[0]<<" x "<<lhs.sizes()[1]<<") = "<<lhs<<"\n"<<
                     "rhs ("<<rhs.sizes()[0]<<" x "<<rhs.sizes()[1]<<") = "<<rhs<<")\n");
        sizes_[0] = lhs.sizes()[0];
    }

    const std::array<size_t, order> & sizes() const { return sizes_; }

    size_t domainDim() const
    {
        gsWarn<<"Correct?\n";
        return this->lhs_expr_.domainDim();
    }


    gsMatrix<Scalar> eval(const index_t k) const
    {
        gsMatrix<Scalar> lhs_val = this->lhs_expr_.eval(k);
        gsMatrix<Scalar> rhs_val = this->rhs_expr_.eval(k);
        return lhs_val * rhs_val;
    }

    void print(std::ostream & os) const override
    {
        os<<this->lhs_expr_<<"*"<<this->rhs_expr_;
    }

private:
    using Base::lhs_expr_;
    using Base::rhs_expr_;
};

// // Generic operator* to create ProductExpression instances
// template <typename LhsExpr, typename RhsExpr>
// ProductExpression<LhsExpr, RhsExpr,LhsExpr::order,RhsExpr::order>
// operator*(const LhsExpr& lhs, const RhsExpr& rhs)
// {
//     return ProductExpression<LhsExpr, RhsExpr,LhsExpr::order,RhsExpr::order>(lhs, rhs);
// }

// Specialization for scalar multiplication
// Specialization: Scalar × Scalar
template <typename LhsExpr, typename RhsExpr>
typename std::enable_if<
    ExpressionTraits<LhsExpr>::order == 0 &&
    ExpressionTraits<RhsExpr>::order == 0,
    ProductExpression<LhsExpr, RhsExpr, 0, 0, ExpressionTraits<LhsExpr>::space, ExpressionTraits<RhsExpr>::space>
>::type
operator*(const LhsExpr& lhs, const RhsExpr& rhs)
{
    return ProductExpression<LhsExpr, RhsExpr, 0, 0, ExpressionTraits<LhsExpr>::space, ExpressionTraits<RhsExpr>::space>(lhs, rhs);
}

// Scalar multiplication operators for primitive types
// Note: These must be carefully written to avoid recursion with ConstantObject

// Disabled to prevent recursion with expression objects
// template <typename Scalar>
// ProductExpression<ConstantObject<Scalar,0>, ConstantObject<Scalar,0>, 0, 0, 0, 0>
// operator*(const Scalar lhs, const Scalar rhs)
// {
//     return ProductExpression<ConstantObject<Scalar,0>, ConstantObject<Scalar,0>, 0, 0, 0, 0>(
//         ConstantObject<Scalar,0>(lhs),
//         ConstantObject<Scalar,0>(rhs)
//     );
// }

// Temporarily disabled to avoid recursion - these need better SFINAE conditions
/*
template <typename RhsExpr>
typename std::enable_if<
    !std::is_same<RhsExpr, ConstantObject<typename ExpressionTraits<RhsExpr>::Scalar,0>>::value,
    decltype(ConstantObject<typename ExpressionTraits<RhsExpr>::Scalar,0>(std::declval<typename ExpressionTraits<RhsExpr>::Scalar>()) * std::declval<RhsExpr>())
>::type
operator*(const typename ExpressionTraits<RhsExpr>::Scalar lhs, const RhsExpr& rhs)
{
    return ConstantObject<typename ExpressionTraits<RhsExpr>::Scalar,0>(lhs) * rhs;
}

template <typename LhsExpr>
typename std::enable_if<
    !std::is_same<LhsExpr, ConstantObject<typename ExpressionTraits<LhsExpr>::Scalar,0>>::value,
    decltype(std::declval<LhsExpr>() * ConstantObject<typename ExpressionTraits<LhsExpr>::Scalar,0>(std::declval<typename ExpressionTraits<LhsExpr>::Scalar>()))
>::type
operator*(const LhsExpr& lhs, const typename ExpressionTraits<LhsExpr>::Scalar rhs)
{
    return lhs * ConstantObject<typename ExpressionTraits<LhsExpr>::Scalar,0>(rhs);
}
*/

// Specialization: Scalar × Non-scalar
template <typename LhsExpr, typename RhsExpr>
typename std::enable_if<
    ExpressionTraits<LhsExpr>::order == 0 &&
    ExpressionTraits<RhsExpr>::order != 0,
    ProductExpression<LhsExpr, RhsExpr, 0, ExpressionTraits<RhsExpr>::order, ExpressionTraits<LhsExpr>::space, ExpressionTraits<RhsExpr>::space>
>::type
operator*(const LhsExpr& lhs, const RhsExpr& rhs)
{
    return ProductExpression<LhsExpr, RhsExpr, 0, ExpressionTraits<RhsExpr>::order, ExpressionTraits<LhsExpr>::space, ExpressionTraits<RhsExpr>::space>(lhs, rhs);
}

template <typename RhsExpr>
typename std::enable_if<
    !std::is_same<RhsExpr, ConstantObject<typename ExpressionTraits<RhsExpr>::Scalar,0>>::value,
    ProductExpression<ConstantObject<typename ExpressionTraits<RhsExpr>::Scalar,0>, RhsExpr, 0, ExpressionTraits<RhsExpr>::order, Space::None, ExpressionTraits<RhsExpr>::space>
>::type
operator*(const typename ExpressionTraits<RhsExpr>::Scalar lhs, const RhsExpr& rhs)
{
    return ProductExpression<ConstantObject<typename ExpressionTraits<RhsExpr>::Scalar,0>, RhsExpr, 0, ExpressionTraits<RhsExpr>::order, Space::None, ExpressionTraits<RhsExpr>::space>(ConstantObject<typename ExpressionTraits<RhsExpr>::Scalar,0>(lhs), rhs);
}

template <typename LhsExpr>
typename std::enable_if<
    !std::is_same<LhsExpr, ConstantObject<typename ExpressionTraits<LhsExpr>::Scalar,0>>::value,
    ProductExpression<LhsExpr, ConstantObject<typename ExpressionTraits<LhsExpr>::Scalar,0>, ExpressionTraits<LhsExpr>::order, 0, ExpressionTraits<LhsExpr>::space, Space::None>
>::type
operator*(const LhsExpr& lhs, const typename ExpressionTraits<LhsExpr>::Scalar rhs)
{
    return ProductExpression<LhsExpr, ConstantObject<typename ExpressionTraits<LhsExpr>::Scalar,0>, ExpressionTraits<LhsExpr>::order, 0, ExpressionTraits<LhsExpr>::space, Space::None>(lhs, ConstantObject<typename ExpressionTraits<LhsExpr>::Scalar,0>(rhs));
}

// Specialization: Expression × Scalar
template <typename LhsExpr, typename RhsExpr>
typename std::enable_if<
    ExpressionTraits<LhsExpr>::order != 0 &&
    ExpressionTraits<RhsExpr>::order == 0,
    ProductExpression<RhsExpr, LhsExpr, 0, ExpressionTraits<LhsExpr>::order, ExpressionTraits<RhsExpr>::space, ExpressionTraits<LhsExpr>::space>
>::type
operator*(const LhsExpr& lhs, const RhsExpr& rhs)
{
    return ProductExpression<RhsExpr, LhsExpr, 0, ExpressionTraits<LhsExpr>::order, ExpressionTraits<RhsExpr>::space, ExpressionTraits<LhsExpr>::space>(rhs, lhs);
}

// Specialization for matrix-vector product
template <typename LhsExpr, typename RhsExpr>
typename std::enable_if<
    ExpressionTraits<LhsExpr>::order == 2 &&
    ExpressionTraits<RhsExpr>::order == 1,
    ProductExpression<LhsExpr, RhsExpr, 2, 1, ExpressionTraits<LhsExpr>::space, ExpressionTraits<RhsExpr>::space>
>::type
operator*(const LhsExpr& mat, const RhsExpr& vec)
{
    return ProductExpression<LhsExpr, RhsExpr, 2, 1, ExpressionTraits<LhsExpr>::space, ExpressionTraits<RhsExpr>::space>(mat, vec);
}

template <typename LhsExpr, typename RhsExpr>
typename std::enable_if<
    ExpressionTraits<LhsExpr>::order == 1 &&
    ExpressionTraits<RhsExpr>::order == 2,
    TransposeExpression<ProductExpression<RhsExpr, LhsExpr, 2, 1, ExpressionTraits<RhsExpr>::space, ExpressionTraits<LhsExpr>::space>>
>::type
operator*(const TransposeExpression<LhsExpr>& vec, const RhsExpr& mat)
{
    return TransposeExpression<ProductExpression<RhsExpr, LhsExpr, 2, 1, ExpressionTraits<RhsExpr>::space, ExpressionTraits<LhsExpr>::space>>(mat, vec.expr());
}

// Matrix-matrix product
template <typename LhsExpr, typename RhsExpr>
typename std::enable_if<
    ExpressionTraits<LhsExpr>::order == 2 &&
    ExpressionTraits<RhsExpr>::order == 2,
    ProductExpression<LhsExpr, RhsExpr, 2, 2, ExpressionTraits<LhsExpr>::space, ExpressionTraits<RhsExpr>::space>
>::type
operator*(const LhsExpr& lhs, const RhsExpr& rhs)
{
    return ProductExpression<LhsExpr, RhsExpr, 2, 2, ExpressionTraits<LhsExpr>::space, ExpressionTraits<RhsExpr>::space>(lhs, rhs);
}

// // Specialization for transpose of a vector with a vector (more specific than generic)
// template <typename LhsExpr, typename RhsExpr>
// auto operator*(const TransposeExpression<LhsExpr>& lhs, const RhsExpr& rhs)
//     -> ProductExpression<TransposeExpression<LhsExpr>, RhsExpr>
// {
//     GISMO_ERROR("Multiplication of a transposed vector and a vector is not defined.");
//     return ProductExpression<TransposeExpression<LhsExpr>, RhsExpr>(lhs, rhs);
// }

// // Specialization for vector with transpose of a vector (more specific than generic)
// template <typename LhsExpr, typename RhsExpr>
// auto operator*(const LhsExpr& lhs, const TransposeExpression<RhsExpr>& rhs)
//     -> ProductExpression<LhsExpr, TransposeExpression<RhsExpr>>
// {
//     GISMO_ERROR("Multiplication of a vector and a transposed vector is not defined.");
//     return ProductExpression<LhsExpr, TransposeExpression<RhsExpr>>(lhs, rhs);
// }

}//namespace Expr
}//namespace gismo