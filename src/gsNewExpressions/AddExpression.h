/** @file BaseObject.h

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

template <typename LhsExpr, typename RhsExpr, typename Enable>
struct ExpressionTraits<AddExpression<LhsExpr, RhsExpr, Enable>>
{
    // Static assertions to ensure compatibility
    static_assert(ExpressionTraits<LhsExpr>::order == ExpressionTraits<RhsExpr>::order,
                  "AddExpression requires same order operands");
    static_assert(ExpressionTraits<LhsExpr>::space == ExpressionTraits<RhsExpr>::space,
                  "AddExpression requires same space operands");

    typedef typename ExpressionTraits<LhsExpr>::Scalar Scalar;

    static constexpr size_t order = ExpressionTraits<LhsExpr>::order;
    static constexpr size_t space = ExpressionTraits<LhsExpr>::space;
    static constexpr size_t deriv = ExpressionTraits<LhsExpr>::deriv;
    static constexpr bool isConstant = ExpressionTraits<LhsExpr>::isConstant &&
                                      ExpressionTraits<RhsExpr>::isConstant;
};

// --- AddExpression using Partial Specialization (Redesigned) ---

// Primary template: Catches all unsupported combinations with a compile-time error
template <typename LhsExpr, typename RhsExpr, typename Enable = void>
class AddExpression
{
    static_assert(std::is_same<LhsExpr, void>::value,
                  "AddExpression: Unsupported tensor order combination for addition.");
};

// --- Partial Specialization 1: Addition of two expressions of the SAME ORDER (X + X) ---
template <typename LhsExpr, typename RhsExpr>
class AddExpression<LhsExpr, RhsExpr,
    typename std::enable_if<(ExpressionTraits<LhsExpr>::order) == (ExpressionTraits<RhsExpr>::order) &&
                            (ExpressionTraits<LhsExpr>::space) == (ExpressionTraits<RhsExpr>::space)>::type
> : public BaseExpression<AddExpression<LhsExpr, RhsExpr>>
{
public:
// Scalar and Order are from ExpressionTraits
    typedef typename ExpressionTraits<LhsExpr>::Scalar Scalar;
    static constexpr size_t order = ExpressionTraits<LhsExpr>::order;
    static constexpr size_t space = ExpressionTraits<LhsExpr>::space;
    // TODO:
    // static constexpr size_t deriv = ExpressionTraits<LhsExpr>::deriv;
    static constexpr bool isConstant = ExpressionTraits<LhsExpr>::isConstant && ExpressionTraits<RhsExpr>::isConstant;

    const std::array<size_t, order> & sizes() const
    {
        return lhs_expr_.sizes(); // Use left operand's sizes
    }

    size_t domainDim() const
    {
        return lhs_expr_.domainDim();
    }

    const LhsExpr& lhs() const {return lhs_expr_;}
    const RhsExpr& rhs() const {return rhs_expr_;}

public:
    AddExpression(const LhsExpr& lhs, const RhsExpr& rhs)
        : lhs_expr_(lhs),
          rhs_expr_(rhs)
    {
    }

    gsMatrix<Scalar> eval(const index_t k) const
    {
        gsMatrix<Scalar> lhs_val = lhs_expr_.eval(k);
        gsMatrix<Scalar> rhs_val = rhs_expr_.eval(k);
        lhs_val.array() += rhs_val.array(); // Element-wise addition
        return lhs_val; // Return the modified lhs_val
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
        os<<lhs_expr_<<"+"<<rhs_expr_;
    }

private:
    const LhsExpr& lhs_expr_;
    const RhsExpr& rhs_expr_;
};

// Generic operator+ to create AddExpression instances
template <typename LhsExpr, typename RhsExpr>
AddExpression<LhsExpr, RhsExpr>
operator+(const LhsExpr& lhs, const RhsExpr& rhs)
{
    return AddExpression<LhsExpr, RhsExpr>(lhs, rhs);
}

// Specialization for transpose of a vector with a vector (more specific than generic)
template <typename LhsExpr, typename RhsExpr>
auto operator+(const TransposeExpression<LhsExpr>& lhs, const RhsExpr& rhs)
    -> AddExpression<TransposeExpression<LhsExpr>, RhsExpr>
{
    GISMO_ERROR("Addition of a transposed vector and a vector is not defined.");
    return AddExpression<TransposeExpression<LhsExpr>, RhsExpr>(lhs, rhs);
}

// Specialization for vector with transpose of a vector (more specific than generic)
template <typename LhsExpr, typename RhsExpr>
auto operator+(const LhsExpr& lhs, const TransposeExpression<RhsExpr>& rhs)
    -> AddExpression<LhsExpr, TransposeExpression<RhsExpr>>
{
    GISMO_ERROR("Addition of a vector and a transposed vector is not defined.");
    return AddExpression<LhsExpr, TransposeExpression<RhsExpr>>(lhs, rhs);
}

// // --- Partial Specialization 1: Addition of two expressions of the SAME ORDER (X + X) ---
// template <typename LhsExpr, typename RhsExpr>
// class AddExpression<ArrayExpression<LhsExpr>, RhsExpr,
//     typename std::enable_if<0 == (ExpressionTraits<RhsExpr>::order)>::type // Simplified condition
// > : public BaseObject<ExpressionTraits<LhsExpr>::Scalar,
//                           ExpressionTraits<LhsExpr>::order,
//                           ExpressionTraits<LhsExpr>::isConstant && ExpressionTraits<RhsExpr>::isConstant,
//                           0> // Use LhsExpr's Scalar and Order directly
// {
// public:
// // Scalar and Order are from ExpressionTraits
//     typedef typename ExpressionTraits<LhsExpr>::Scalar Scalar;
//     static constexpr size_t order = ExpressionTraits<LhsExpr>::order;

// public:
//     AddExpression(const LhsExpr& lhs, const RhsExpr& rhs)
//         : BaseObject<typename ExpressionTraits<LhsExpr>::Scalar, order, ExpressionTraits<LhsExpr>::isConstant && ExpressionTraits<RhsExpr>::isConstant, 0>(rhs.sizes()), // Pass RhsExpr as Derived to BaseObject
//           lhs_expr_(lhs),
//           rhs_expr_(rhs)
//     {
//     }

//     gsMatrix<Scalar> eval(const index_t k) const
//     {
//         gsMatrix<Scalar> lhs_val = lhs_expr_.eval(k);
//         gsMatrix<Scalar> rhs_val = rhs_expr_.eval(k);
//         lhs_val.array() += rhs_val; // Element-wise addition
//         return lhs_val; // Return the modified lhs_val
//     }

// private:
//     const LhsExpr& lhs_expr_;
//     const RhsExpr& rhs_expr_;
// };



// // --- Partial Specialization 2: Scalar (Order 0) + Higher Order (Order N > 0) ---
// template <typename LhsExpr, typename RhsExpr>
// class AddExpression<LhsExpr, RhsExpr,
//     typename std::enable_if<(ExpressionTraits<LhsExpr>::order == 0) && (ExpressionTraits<RhsExpr>::order > 0)>::type // Simplified condition
// > : public BaseObject<typename ExpressionTraits<LhsExpr>::Scalar, ExpressionTraits<RhsExpr>::order> // Base on RhsExpr for Scalar and Order
// {
// public:
//     typedef typename ExpressionTraits<RhsExpr>::Scalar Scalar;
//     static constexpr size_t order = ExpressionTraits<RhsExpr>::order;

// public:
//     AddExpression(const LhsExpr& lhs, const RhsExpr& rhs)
//         : BaseObject<typename ExpressionTraits<LhsExpr>::Scalar, order>(rhs.sizes()), // Pass RhsExpr as Derived to BaseObject
//           lhs_expr_(lhs),
//           rhs_expr_(rhs)
//     {}

//     gsMatrix<Scalar> eval(const index_t k) const
//     {
//         gsMatrix<Scalar> lhs_val = lhs_expr_.eval(k);
//         gsMatrix<Scalar> rhs_val = rhs_expr_.eval(k);
//         rhs_val.array() += lhs_val.value(); // Add scalar to each element of lhs_val
//         return rhs_val;
//     }

// private:
//     const LhsExpr& lhs_expr_;
//     const RhsExpr& rhs_expr_;
// };


// // Scalar + Scalar
// AddExpression<ScalarExpression<real_t>, ScalarExpression<real_t>> operator+(const ScalarExpression<real_t>& lhs, const ScalarExpression<real_t>& rhs)
// {
//     return AddExpression<ScalarExpression<real_t>, ScalarExpression<real_t>>(lhs, rhs);
// }

// // Scalar + Vector
// AddExpression<ScalarExpression<real_t>, VectorExpression<real_t>> operator+(const ScalarExpression<real_t>& lhs, const VectorExpression<real_t>& rhs)
// {
//     return AddExpression<ScalarExpression<real_t>, VectorExpression<real_t>>(lhs, rhs);
// }

// // Vector + Scalar
// AddExpression<ScalarExpression<real_t>, VectorExpression<real_t>> operator+(const VectorExpression<real_t>& lhs, const ScalarExpression<real_t>& rhs)
// {
//     return AddExpression<ScalarExpression<real_t>, VectorExpression<real_t>>(rhs, lhs);
// }

// // Vector + Vector
// AddExpression<VectorExpression<real_t>, VectorExpression<real_t>> operator+(const VectorExpression<real_t>& lhs, const VectorExpression<real_t>& rhs)
// {
//     return AddExpression<VectorExpression<real_t>, VectorExpression<real_t>>(lhs, rhs);
// }

// // Scalar + Matrix
// AddExpression<ScalarExpression<real_t>, MatrixExpression<real_t>> operator+(const ScalarExpression<real_t>& lhs, const MatrixExpression<real_t>& rhs)
// {
//     return AddExpression<ScalarExpression<real_t>, MatrixExpression<real_t>>(lhs, rhs);
// }

// // Matrix + Scalar
// AddExpression<ScalarExpression<real_t>, MatrixExpression<real_t>> operator+(const MatrixExpression<real_t>& lhs, const ScalarExpression<real_t>& rhs)
// {
//     return AddExpression<ScalarExpression<real_t>, MatrixExpression<real_t>>(rhs, lhs);
// }

}//namespace Expr
}//namespace gismo