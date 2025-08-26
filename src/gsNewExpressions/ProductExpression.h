/** @file BaseObject.h

    @brief

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M.Verhelst
*/

#pragma once

// Forward declarations
namespace gismo { namespace Expr {
template <typename E> class TransposeExpression;
} }

namespace gismo
{
namespace Expr
{

// --- ProductExpression using Partial Specialization (Redesigned) ---

// Primary template: Catches all unsupported combinations with a compile-time error
template <typename LhsExpr, typename RhsExpr, typename Enable = void>
class ProductExpression
{
    static_assert(std::is_same<LhsExpr, void>::value,
                  "ProductExpression: Unsupported tensor order combination for product.");
};

// --- Partial specialization 1: Scalar multiplication
template <typename LhsExpr, typename RhsExpr>
struct ExpressionTraits<ProductExpression<LhsExpr, RhsExpr,
    typename std::enable_if<(ExpressionTraits<LhsExpr>::order==0 &&
                             ExpressionTraits<LhsExpr>::space==Space::None &&
                             ExpressionTraits<RhsExpr>::space==Space::None)>::type>>
{
    typedef typename ExpressionTraits<RhsExpr>::Scalar Scalar;

    static constexpr size_t order = ExpressionTraits<RhsExpr>::order;
    static constexpr size_t space = ExpressionTraits<RhsExpr>::space;
    // CORRECT??
    static constexpr size_t deriv = ExpressionTraits<RhsExpr>::deriv;
    static constexpr bool isConstant = ExpressionTraits<LhsExpr>::isConstant &&
                                      ExpressionTraits<RhsExpr>::isConstant;
};

template <typename LhsExpr, typename RhsExpr>
class ProductExpression<LhsExpr, RhsExpr,
    typename std::enable_if<(ExpressionTraits<LhsExpr>::order==0 &&
                             ExpressionTraits<LhsExpr>::space==Space::None &&
                             ExpressionTraits<RhsExpr>::space==Space::None)>::type
> : public BaseExpression<ProductExpression<LhsExpr, RhsExpr>>
{
public:
// Use ExpressionTraits for this specific ProductExpression
    typedef typename ExpressionTraits<ProductExpression<LhsExpr, RhsExpr>>::Scalar Scalar;
    static constexpr size_t order = ExpressionTraits<ProductExpression<LhsExpr, RhsExpr>>::order;
    static constexpr size_t space = ExpressionTraits<ProductExpression<LhsExpr, RhsExpr>>::space;
    static constexpr size_t deriv = ExpressionTraits<ProductExpression<LhsExpr, RhsExpr>>::deriv;
    static constexpr bool isConstant = ExpressionTraits<ProductExpression<LhsExpr, RhsExpr>>::isConstant;

    const std::array<size_t, order> & sizes() const { return rhs_expr_.sizes(); }

    size_t domainDim() const
    {
        gsWarn<<"Correct?\n";
        return lhs_expr_.domainDim();
    }

    const LhsExpr& lhs() const {return lhs_expr_;}
    const RhsExpr& rhs() const {return rhs_expr_;}

public:
    ProductExpression(const LhsExpr& lhs, const RhsExpr& rhs)
        : lhs_expr_(lhs),
          rhs_expr_(rhs)
    {
    }

    gsMatrix<Scalar> eval(const index_t k) const
    {
        gsMatrix<Scalar> lhs_val = lhs_expr_.eval(k);
        gsMatrix<Scalar> rhs_val = rhs_expr_.eval(k);
        rhs_val.array() *= lhs_val.value(); // Element-wise multiplication
        return rhs_val; // Return the modified rhs_val
    }

    void parse(gismo::ExpressionHelper<Scalar> & helper) const
    {
        lhs_expr_.parse(helper);
        rhs_expr_.parse(helper);
    }

    // TODO
    // const SpaceObject<Scalar, space, order> & rowVar() const
    // {
    //     return lhs_expr_.rowVar(); // Use left operand's row variable
    // }

    // const SpaceObject<Scalar, space, order> & colVar() const
    // {
    //     return lhs_expr_.colVar(); // Use left operand's column variable
    // }

    void print(std::ostream & os) const
    {
        os<<lhs_expr_<<"*"<<rhs_expr_;
    }

private:
    const LhsExpr& lhs_expr_;
    const RhsExpr& rhs_expr_;
};

// --- Partial specialization 2: Multiplication of equal order (matrix only)
template <typename LhsExpr, typename RhsExpr>
struct ExpressionTraits<ProductExpression<LhsExpr, RhsExpr,
    typename std::enable_if<(ExpressionTraits<LhsExpr>::order==2 &&
                             ExpressionTraits<LhsExpr>::space==Space::None &&
                             ExpressionTraits<RhsExpr>::order==2 &&
                             ExpressionTraits<RhsExpr>::space==Space::None)>::type>>
{
    typedef typename ExpressionTraits<RhsExpr>::Scalar Scalar;

    static constexpr size_t order = ExpressionTraits<RhsExpr>::order;
    static constexpr size_t space = ExpressionTraits<RhsExpr>::space;
    // CORRECT??
    static constexpr size_t deriv = ExpressionTraits<RhsExpr>::deriv;
    static constexpr bool isConstant = ExpressionTraits<LhsExpr>::isConstant &&
                                      ExpressionTraits<RhsExpr>::isConstant;
};
template <typename LhsExpr, typename RhsExpr>
class ProductExpression<LhsExpr, RhsExpr,
    typename std::enable_if<(ExpressionTraits<LhsExpr>::order==2 &&
                             ExpressionTraits<LhsExpr>::space==Space::None &&
                             ExpressionTraits<RhsExpr>::order==2 &&
                             ExpressionTraits<RhsExpr>::space==Space::None)>::type
> : public BaseExpression<ProductExpression<LhsExpr, RhsExpr>>
{
public:
// Scalar and Order are from ExpressionTraits
    typedef typename ExpressionTraits<RhsExpr>::Scalar Scalar;
    static constexpr size_t order = 2;
    static constexpr size_t space = ExpressionTraits<RhsExpr>::space;
    // TODO:
    static constexpr size_t deriv = ExpressionTraits<RhsExpr>::deriv;
    static constexpr bool isConstant = ExpressionTraits<LhsExpr>::isConstant && ExpressionTraits<RhsExpr>::isConstant;

    const std::array<size_t, order> & sizes() const { return sizes_; }

    size_t domainDim() const
    {
        gsWarn<<"Correct?\n";
        return lhs_expr_.domainDim();
    }

    const LhsExpr& lhs() const {return lhs_expr_;}
    const RhsExpr& rhs() const {return rhs_expr_;}

protected:

    std::array<size_t,order> sizes_;

public:
    ProductExpression(const LhsExpr& lhs, const RhsExpr& rhs)
        : lhs_expr_(lhs),
          rhs_expr_(rhs)
    {
        GISMO_ENSURE(lhs.sizes()[1]==rhs.sizes()[0],
                     "Size mismatch in * operator.\n"<<
                     "lhs ("<<lhs.sizes()[0]<<" x "<<lhs.sizes()[1]<<") = "<<lhs<<"\n"<<
                     "rhs ("<<rhs.sizes()[0]<<" x "<<rhs.sizes()[1]<<") = "<<rhs<<")\n");
        sizes_[0] = lhs.sizes()[0];
        sizes_[1] = rhs.sizes()[1];
    }

    gsMatrix<Scalar> eval(const index_t k) const
    {
        gsMatrix<Scalar> lhs_val = lhs_expr_.eval(k);
        gsMatrix<Scalar> rhs_val = rhs_expr_.eval(k);
        return lhs_val * rhs_val;
    }

    void parse(gismo::ExpressionHelper<Scalar> & helper) const
    {
        lhs_expr_.parse(helper);
        rhs_expr_.parse(helper);
    }

    const SpaceObject<Scalar, space, order> & rowVar() const
    {
        return lhs_expr_.rowVar(); // Use left operand's row variable
    }

    const SpaceObject<Scalar, space, order> & colVar() const
    {
        return lhs_expr_.colVar(); // Use left operand's column variable
    }

    void print(std::ostream & os) const
    {
        os<<lhs_expr_<<"*"<<rhs_expr_;
    }

private:
    const LhsExpr& lhs_expr_;
    const RhsExpr& rhs_expr_;
};

// --- Partial specialization 2: Multiplication of matrix and vector
template <typename LhsExpr, typename RhsExpr>
class ProductExpression<LhsExpr, RhsExpr,
    typename std::enable_if<(ExpressionTraits<LhsExpr>::order==2 &&
                             ExpressionTraits<LhsExpr>::space==Space::None &&
                             ExpressionTraits<RhsExpr>::order==1 &&
                             ExpressionTraits<RhsExpr>::space==Space::None)>::type
> : public BaseExpression<ProductExpression<LhsExpr, RhsExpr>>
{

public:
// Scalar and Order are from ExpressionTraits
    typedef typename ExpressionTraits<RhsExpr>::Scalar Scalar;
    static constexpr size_t order = 1;
    static constexpr size_t space = ExpressionTraits<RhsExpr>::space;
    // TODO:
    // static constexpr size_t deriv = ExpressionTraits<RhsExpr>::deriv;
    static constexpr bool isConstant = ExpressionTraits<LhsExpr>::isConstant && ExpressionTraits<RhsExpr>::isConstant;

    const std::array<size_t, order> & sizes() const { return sizes_; }

    size_t domainDim() const
    {
        gsWarn<<"Correct?\n";
        return lhs_expr_.domainDim();
    }

    const LhsExpr& lhs() const {return lhs_expr_;}
    const RhsExpr& rhs() const {return rhs_expr_;}

protected:

    std::array<size_t,order> sizes_;

public:
    ProductExpression(const LhsExpr& lhs, const RhsExpr& rhs)
        : lhs_expr_(lhs),
          rhs_expr_(rhs)
    {
        GISMO_ENSURE(lhs.sizes()[1]==rhs.sizes()[0],
                     "Size mismatch in * operator.\n"<<
                     "lhs ("<<lhs.sizes()[0]<<" x "<<lhs.sizes()[1]<<") = "<<lhs<<"\n"<<
                     "rhs ("<<rhs.sizes()[0]<<" x "<<rhs.sizes()[1]<<") = "<<rhs<<")\n");
        sizes_[0] = lhs.sizes()[0];
    }

    gsMatrix<Scalar> eval(const index_t k) const
    {
        gsMatrix<Scalar> lhs_val = lhs_expr_.eval(k);
        gsMatrix<Scalar> rhs_val = rhs_expr_.eval(k);
        return lhs_val * rhs_val;
    }

    void parse(gismo::ExpressionHelper<Scalar> & helper) const
    {
        lhs_expr_.parse(helper);
        rhs_expr_.parse(helper);
    }

    const SpaceObject<Scalar, space, order> & rowVar() const
    {
        return lhs_expr_.rowVar(); // Use left operand's row variable
    }

    const SpaceObject<Scalar, space, order> & colVar() const
    {
        return lhs_expr_.colVar(); // Use left operand's column variable
    }

    void print(std::ostream & os) const
    {
        os<<lhs_expr_<<"*"<<rhs_expr_;
    }

private:
    const LhsExpr& lhs_expr_;
    const RhsExpr& rhs_expr_;
};


// Generic operator* to create ProductExpression instances
template <typename LhsExpr, typename RhsExpr>
ProductExpression<LhsExpr, RhsExpr>
operator*(const LhsExpr& lhs, const RhsExpr& rhs)
{
    return ProductExpression<LhsExpr, RhsExpr>(lhs, rhs);
}

// Specialization for scalar multiplication
// Specialization: Expression × Scalar -> Scalar × Expression
template <typename LhsExpr, typename RhsExpr>
typename std::enable_if<
    (ExpressionTraits<RhsExpr>::order == 0 && ExpressionTraits<RhsExpr>::space == Space::None &&
     !(ExpressionTraits<LhsExpr>::order == 0 && ExpressionTraits<LhsExpr>::space == Space::None)),
    ProductExpression<RhsExpr, LhsExpr>  // Note: swapped order
>::type
operator*(const LhsExpr& lhs, const RhsExpr& rhs)
{
    return ProductExpression<RhsExpr, LhsExpr>(rhs, lhs);  // Swap arguments
}

// Specialization for matrix-vector product
template <typename LhsExpr, typename RhsExpr>
typename std::enable_if<
    (ExpressionTraits<LhsExpr>::order == 1 && ExpressionTraits<RhsExpr>::order == 2 &&
     ExpressionTraits<LhsExpr>::space == Space::None && ExpressionTraits<RhsExpr>::space == Space::None),
    ProductExpression<RhsExpr,LhsExpr>
>::type
operator*(const LhsExpr& lhs, const RhsExpr& rhs)
{
    return TransposeExpression<ProductExpression<RhsExpr, LhsExpr>>(rhs, lhs);
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