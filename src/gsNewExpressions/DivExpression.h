/** @file DivExpression.h

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

template <typename E, typename Enable>
struct ExpressionTraits<DivExpression<E, Enable>>
{
    // Static assertions to ensure compatibility
    // static_assert(ExpressionTraits<E>::order == 0,
    //               "DivExpression requires a scalar (order 0) operand");
    // static_assert(ExpressionTraits<E>::space != Space::None,
    //               "DivExpression requires the operand to be defined in a space");

    typedef typename ExpressionTraits<E>::Scalar Scalar;
    static constexpr size_t order = ExpressionTraits<E>::order-1; // Divergence results in a vector (order 1)
    static constexpr size_t space = ExpressionTraits<E>::space;
    static constexpr size_t deriv = ExpressionTraits<E>::deriv + 1; // Increment derivative order
    static constexpr bool isConstant = ExpressionTraits<E>::isConstant;
};

// --- DivExpression using Partial Specialization (Redesigned) ---
// Primary template: Catches all unsupported combinations with a compile-time error
template <typename E, typename Enable = void>
class DivExpression : public BaseExpression<DivExpression<E>>
{
    static_assert(std::is_same<Enable, void>::value,
                  "TODO");
};

// --- Partial Specialization 1: Divergence of a constant expression ---
template <typename E>
class DivExpression<E,
    typename std::enable_if<(ExpressionTraits<E>::isConstant)>::type
> : public BaseExpression<DivExpression<E, typename std::enable_if<(ExpressionTraits<E>::isConstant)>::type>>
{
    using Base = BaseExpression<DivExpression<E, typename std::enable_if<(ExpressionTraits<E>::isConstant)>::type>>;

public:
    typedef typename ExpressionTraits<DivExpression<E>>::Scalar Scalar;
    static constexpr size_t order = ExpressionTraits<DivExpression<E>>::order;
    static constexpr size_t space = ExpressionTraits<DivExpression<E>>::space;
    static constexpr size_t deriv = ExpressionTraits<DivExpression<E>>::deriv;
    static constexpr bool isConstant = ExpressionTraits<DivExpression<E>>::isConstant;

    const std::array<size_t, order> & sizes() const
    {
        return sizes_;
    }

    size_t domainDim() const
    {
        return expr_.domainDim();
    }

private:
    std::array<size_t,order> sizes_;
    mutable gsMatrix<Scalar> tmp;
    const E& expr_;

public:
    DivExpression(const E& expr)
    :
    BaseExpression<DivExpression<E>>(),
    expr_(expr)
    {
        // Assumes symmetry!
        for (short_t d=0; d!=E::order; d++)
        {
            GISMO_ENSURE(expr_.sizes()[d] == expr_.sizes()[0]);
            GISMO_ENSURE(expr_.sizes()[d] == expr_.domainDim());
        }

        for (short_t d=1; d!=E::order; d++)
            sizes_[d-1] = expr_.sizes()[d];
    }

    gsMatrix<Scalar> eval(const index_t k) const
    {
        GISMO_UNUSED(k);
        tmp.resize(expr_.source().domainDim(),1);
        tmp.setZero();
        return tmp;
    }

    void parse(gismo::ExpressionHelper<Scalar> & helper) const
    {
    }

};

// --- Partial Specialization 2: divient of a variable object ---
template <typename E>
class DivExpression<E,
    typename std::enable_if<(ExpressionTraits<E>::isConstant==false)>::type
> : public BaseExpression<DivExpression<E, typename std::enable_if<(ExpressionTraits<E>::isConstant==false)>::type>>
{
    using Base = BaseExpression<DivExpression<E, typename std::enable_if<(ExpressionTraits<E>::isConstant==false)>::type>>;
public:
    typedef typename ExpressionTraits<DivExpression<E>>::Scalar Scalar;
    static constexpr size_t order = ExpressionTraits<DivExpression<E>>::order;
    static constexpr size_t space = ExpressionTraits<DivExpression<E>>::space;
    static constexpr size_t deriv = ExpressionTraits<DivExpression<E>>::deriv;
    static constexpr bool isConstant = ExpressionTraits<DivExpression<E>>::isConstant;

    const std::array<size_t, order> & sizes() const
    {
        return sizes_;
    }

    size_t domainDim() const
    {
        return expr_.domainDim();
    }

private:
    std::array<size_t,order> sizes_;
    mutable gsMatrix<Scalar> tmp;
    E & expr_;

public:
    DivExpression(E & expr)
    :
    BaseExpression<DivExpression<E>>(),
    expr_(expr)
    {
        // Assumes symmetry!
        for (short_t d=0; d!=E::order; d++)
        {
            GISMO_ENSURE(expr_.sizes()[d] == expr_.sizes()[0]);
            GISMO_ENSURE(expr_.sizes()[d] == expr_.domainDim());
        }

        for (short_t d=1; d!=E::order; d++)
            sizes_[d-1] = expr_.sizes()[d];
    }

    gsMatrix<Scalar> eval(const index_t k) const
    {
        GISMO_UNUSED(k);
        tmp.resize(expr_.domainDim(),1);
        tmp.setZero();
        return tmp;
    }

    void parse(gismo::ExpressionHelper<Scalar> & helper) const
    {
        expr_.parse(helper);
        expr_.setDerivative(deriv);
    }

    void print(std::ostream & os) const
    {
        os<<"\u2207\u2022("<<expr_<<")";
    }

};

// Generic factory function for easy creation
template <typename E>
DivExpression<E> div(E& expr)
{
    return DivExpression<E>(expr);
}

// Partial specialization for addition
template <typename LhsExpr, typename RhsExpr>
auto div(const AddExpression<LhsExpr,RhsExpr> expr)
-> AddExpression<DivExpression<LhsExpr>, DivExpression<RhsExpr>>
{
    return AddExpression<DivExpression<LhsExpr>, DivExpression<RhsExpr>>(div(expr.lhs()),div(expr.rhs()));
}

// --- Partial Specialization 2: Divergence of a VariableObject ---



// // --- Partial Specialization 1: Addition of two expressions of the SAME ORDER (X + X) ---
// template <typename LhsExpr, typename RhsExpr>
// class AddExpression<ArrayExpression<LhsExpr>, RhsExpr,
//     typename std::enable_if<0 == (RhsExpr::order)>::type // Simplified condition
// > : public BaseObject<LhsExpr,
//                           typename LhsExpr::Scalar,
//                           LhsExpr::order,
//                           LhsExpr::isConstant && RhsExpr::isConstant,
//                           0> // Use LhsExpr's Scalar and Order directly
// {
// public:
// // Scalar and Order are directly from LhsExpr/RhsExpr
//     typedef typename LhsExpr::Scalar Scalar;
//     static constexpr int order = LhsExpr::order;

// public:
//     AddExpression(const LhsExpr& lhs, const RhsExpr& rhs)
//         : BaseObject<RhsExpr, typename LhsExpr::Scalar, order, LhsExpr::isConstant && RhsExpr::isConstant, 0>(rhs.sizes()), // Pass RhsExpr as Derived to BaseObject
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
//     typename std::enable_if<(LhsExpr::order == 0) && (RhsExpr::order > 0)>::type // Simplified condition
// > : public BaseObject<RhsExpr, typename LhsExpr::Scalar, RhsExpr::order> // Base on RhsExpr for Scalar and Order
// {
// public:
//     typedef typename RhsExpr::Scalar Scalar;
//     static constexpr int order = RhsExpr::order;

// public:
//     AddExpression(const LhsExpr& lhs, const RhsExpr& rhs)
//         : BaseObject<RhsExpr, typename LhsExpr::Scalar, order>(rhs.sizes()), // Pass RhsExpr as Derived to BaseObject
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