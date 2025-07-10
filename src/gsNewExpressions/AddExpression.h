/** @file BaseExpression.h

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
    typename std::enable_if<(LhsExpr::order) == (RhsExpr::order)>::type // Simplified condition
> : public BaseExpression<LhsExpr,
                          typename LhsExpr::Scalar,
                          LhsExpr::order,
                          LhsExpr::isConstant && RhsExpr::isConstant> // Use LhsExpr's Scalar and Order directly
{
public:
// Scalar and Order are directly from LhsExpr/RhsExpr
    typedef typename LhsExpr::Scalar Scalar;
    static constexpr int order = LhsExpr::order;

public:
    AddExpression(const LhsExpr& lhs, const RhsExpr& rhs)
        : BaseExpression<RhsExpr, typename LhsExpr::Scalar, order>(rhs.sizes()), // Pass RhsExpr as Derived to BaseExpression
          lhs_expr_(lhs),
          rhs_expr_(rhs)
    {}

    gsMatrix<Scalar> eval(const index_t k) const
    {
        gsMatrix<Scalar> lhs_val = lhs_expr_.eval(k);
        gsMatrix<Scalar> rhs_val = rhs_expr_.eval(k);
        lhs_val.array() += rhs_val.array(); // Element-wise addition
        return lhs_val; // Return the modified lhs_val
    }

private:
    const LhsExpr& lhs_expr_;
    const RhsExpr& rhs_expr_;
};



// // --- Partial Specialization 2: Scalar (Order 0) + Higher Order (Order N > 0) ---
// template <typename LhsExpr, typename RhsExpr>
// class AddExpression<LhsExpr, RhsExpr,
//     typename std::enable_if<(LhsExpr::order == 0) && (RhsExpr::order > 0)>::type // Simplified condition
// > : public BaseExpression<RhsExpr, typename LhsExpr::Scalar, RhsExpr::order> // Base on RhsExpr for Scalar and Order
// {
// public:
//     typedef typename RhsExpr::Scalar Scalar;
//     static constexpr int order = RhsExpr::order;

// public:
//     AddExpression(const LhsExpr& lhs, const RhsExpr& rhs)
//         : BaseExpression<RhsExpr, typename LhsExpr::Scalar, order>(rhs.sizes()), // Pass RhsExpr as Derived to BaseExpression
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


// Scalar + Scalar
AddExpression<ScalarExpression<real_t>, ScalarExpression<real_t>> operator+(const ScalarExpression<real_t>& lhs, const ScalarExpression<real_t>& rhs)
{
    return AddExpression<ScalarExpression<real_t>, ScalarExpression<real_t>>(lhs, rhs);
}

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

// Vector + Vector
AddExpression<VectorExpression<real_t>, VectorExpression<real_t>> operator+(const VectorExpression<real_t>& lhs, const VectorExpression<real_t>& rhs)
{
    return AddExpression<VectorExpression<real_t>, VectorExpression<real_t>>(lhs, rhs);
}

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