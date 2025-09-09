/** @file AddExpression.h

    @brief Addition expression class

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


// --- Partial Specialization: Addition of two expressions of the SAME ORDER (X + X) ---
// ExpressionTraits for same-order addition (N,N) -> N
template <typename LhsExpr, typename RhsExpr, size_t Order, size_t LhsSpace, size_t RhsSpace>
struct ExpressionTraits<AddExpression<LhsExpr, RhsExpr, Order, Order, LhsSpace, RhsSpace>>
{
    // Static assertions to ensure compatibility
    static_assert(ExpressionTraits<LhsExpr>::order == ExpressionTraits<RhsExpr>::order,
                  "AddExpression requires same order operands");
    static_assert(ExpressionTraits<LhsExpr>::space == ExpressionTraits<RhsExpr>::space,
                  "AddExpression requires same space operands");

    typedef LhsExpr LhsType;
    typedef RhsExpr RhsType;

    typedef typename ExpressionTraits<LhsExpr>::Scalar Scalar;
    static constexpr size_t order = Order;
    static constexpr size_t space = ExpressionTraits<LhsExpr>::space;
    static constexpr size_t deriv = ExpressionTraits<RhsExpr>::deriv; //TODO
    static constexpr bool isConstant = ExpressionTraits<LhsExpr>::isConstant && ExpressionTraits<RhsExpr>::isConstant;
};

template <typename LhsExpr, typename RhsExpr, size_t Order, size_t LhsSpace, size_t RhsSpace>
class AddExpression<LhsExpr, RhsExpr, Order, Order, LhsSpace, RhsSpace> : public BinaryOperator<AddExpression<LhsExpr, RhsExpr, Order, Order, LhsSpace, RhsSpace>>
{
    using Base = BinaryOperator<AddExpression<LhsExpr, RhsExpr, Order, Order, LhsSpace, RhsSpace>>;
public:
    AddExpression(const LhsExpr& lhs, const RhsExpr& rhs)
        : Base(lhs, rhs)
    {
        for (size_t d=0; d!=Base::order; ++d)
        {
            GISMO_ENSURE(this->lhs_expr_.sizes()[d] == this->rhs_expr_.sizes()[d],"AddExpression requires same sizes in each dimension, but got "+std::to_string(this->lhs_expr_.sizes()[d])+" and "+std::to_string(this->rhs_expr_.sizes()[d])+" in dimension "+std::to_string(d));
        }
        // Copy sizes from left operand to base class sizes_
        for (size_t d=0; d!=Base::order; ++d)
        {
            this->sizes_[d] = this->lhs_expr_.sizes()[d];
        }
    }

    size_t domainDim() const
    {
        return this->lhs_expr_.domainDim(); // Use left operand's domain dimension
    }

    gsMatrix<typename Base::Scalar> eval(const index_t k) const
    {
        gsMatrix<typename Base::Scalar> lhs_val = this->lhs_expr_.eval(k);
        gsMatrix<typename Base::Scalar> rhs_val = this->rhs_expr_.eval(k);
        lhs_val.array() += rhs_val.array(); // Element-wise addition
        return lhs_val; // Return the modified lhs_val
    }

    void print(std::ostream & os) const
    {
        // gsDebug<<"AddExpression print called\n";
        // this->lhs_expr_.print(os);
        // gsDebug<<this->lhs_expr_<<"\n";
        // gsDebug<<this->rhs_expr_<<"\n";
        os << this->lhs_expr_ << "+" << this->rhs_expr_;
    }

private:
};

// Partial Specialization: Addition of scalar (0) and higher order (N > 0) is now handled by operator+ SFINAE functions

// TODO: Temporarily disabled due to ambiguity with (0,0) case - need proper SFINAE
/*
// --- Partial Specialization: Addition of scalar (0) and higher order (N > 0) ---
template <typename LhsExpr, typename RhsExpr, size_t RhsOrder, size_t LhsSpace, size_t RhsSpace>
class AddExpression<LhsExpr, RhsExpr, 0, RhsOrder, LhsSpace, RhsSpace>
 : public BinaryOperator<LhsExpr, RhsExpr, 0, RhsOrder, LhsSpace, RhsSpace>
{
    using Base = BinaryOperator<LhsExpr, RhsExpr, 0, RhsOrder, LhsSpace, RhsSpace>;
public:
    typedef typename ExpressionTraits<AddExpression<LhsExpr, RhsExpr, 0, RhsOrder, LhsSpace, RhsSpace>>::Scalar Scalar;
    static constexpr size_t order = ExpressionTraits<AddExpression<LhsExpr, RhsExpr, 0, RhsOrder, LhsSpace, RhsSpace>>::order;
    static constexpr size_t space = ExpressionTraits<AddExpression<LhsExpr, RhsExpr, 0, RhsOrder, LhsSpace, RhsSpace>>::space;
    static constexpr size_t deriv = ExpressionTraits<AddExpression<LhsExpr, RhsExpr, 0, RhsOrder, LhsSpace, RhsSpace>>::deriv;
    static constexpr bool isConstant = ExpressionTraits<AddExpression<LhsExpr, RhsExpr, 0, RhsOrder, LhsSpace, RhsSpace>>::isConstant;

public:
    AddExpression(const LhsExpr& lhs, const RhsExpr& rhs)
        : BinaryOperator<LhsExpr, RhsExpr, 0, RhsOrder, LhsSpace, RhsSpace>(lhs, rhs)
    {
    }

    gsMatrix<Scalar> eval(const index_t k) const
    {
        gsMatrix<Scalar> lhs_val = this->lhs_expr_.eval(k);
        gsMatrix<Scalar> rhs_val = this->rhs_expr_.eval(k);
        rhs_val.array() += lhs_val(0,0); // Add scalar to each element
        return rhs_val;
    }

    void print(std::ostream & os) const
    {
        os << this->lhs_expr_ << "+" << this->rhs_expr_;
    }

    const std::array<size_t, order> & sizes() const
    {
        return this->rhs_expr_.sizes();
    }

    size_t domainDim() const
    {
        return this->rhs_expr_.domainDim();
    }

private:
};

// --- Partial Specialization: Addition of higher order (N > 0) and scalar (0) ---
template <typename LhsExpr, typename RhsExpr, size_t LhsOrder, size_t LhsSpace, size_t RhsSpace>
class AddExpression<LhsExpr, RhsExpr, LhsOrder, 0, LhsSpace, RhsSpace>
 : public BinaryOperator<
     typename ExpressionTraits<AddExpression<LhsExpr, RhsExpr, LhsOrder, 0, LhsSpace, RhsSpace>>::Scalar,
     ExpressionTraits<AddExpression<LhsExpr, RhsExpr, LhsOrder, 0, LhsSpace, RhsSpace>>::order,
     ExpressionTraits<AddExpression<LhsExpr, RhsExpr, LhsOrder, 0, LhsSpace, RhsSpace>>::space,
     LhsExpr, RhsExpr, LhsOrder, 0, LhsSpace, RhsSpace>
{
    using Base = BinaryOperator<
        typename ExpressionTraits<AddExpression<LhsExpr, RhsExpr, LhsOrder, 0, LhsSpace, RhsSpace>>::Scalar,
        ExpressionTraits<AddExpression<LhsExpr, RhsExpr, LhsOrder, 0, LhsSpace, RhsSpace>>::order,
        ExpressionTraits<AddExpression<LhsExpr, RhsExpr, LhsOrder, 0, LhsSpace, RhsSpace>>::space,
        LhsExpr, RhsExpr, LhsOrder, 0, LhsSpace, RhsSpace>;
public:
    typedef typename ExpressionTraits<AddExpression<LhsExpr, RhsExpr, LhsOrder, 0, LhsSpace, RhsSpace>>::Scalar Scalar;
    static constexpr size_t order = ExpressionTraits<AddExpression<LhsExpr, RhsExpr, LhsOrder, 0, LhsSpace, RhsSpace>>::order;
    static constexpr size_t space = ExpressionTraits<AddExpression<LhsExpr, RhsExpr, LhsOrder, 0, LhsSpace, RhsSpace>>::space;
    static constexpr size_t deriv = ExpressionTraits<AddExpression<LhsExpr, RhsExpr, LhsOrder, 0, LhsSpace, RhsSpace>>::deriv;
    static constexpr bool isConstant = ExpressionTraits<AddExpression<LhsExpr, RhsExpr, LhsOrder, 0, LhsSpace, RhsSpace>>::isConstant;

public:
    AddExpression(const LhsExpr& lhs, const RhsExpr& rhs)
        : BinaryOperator<
            typename ExpressionTraits<AddExpression<LhsExpr, RhsExpr, LhsOrder, 0, LhsSpace, RhsSpace>>::Scalar,
            ExpressionTraits<AddExpression<LhsExpr, RhsExpr, LhsOrder, 0, LhsSpace, RhsSpace>>::order,
            ExpressionTraits<AddExpression<LhsExpr, RhsExpr, LhsOrder, 0, LhsSpace, RhsSpace>>::space,
            LhsExpr, RhsExpr, LhsOrder, 0, LhsSpace, RhsSpace>(lhs, rhs)
    {
    }

    gsMatrix<Scalar> eval(const index_t k) const
    {
        gsMatrix<Scalar> lhs_val = this->lhs_expr_.eval(k);
        gsMatrix<Scalar> rhs_val = this->rhs_expr_.eval(k);
        lhs_val.array() += rhs_val(0,0); // Add scalar to each element
        return lhs_val;
    }

    void print(std::ostream & os) const
    {
        os << this->lhs_expr_ << "+" << this->rhs_expr_;
    }

    const std::array<size_t, order> & sizes() const
    {
        return this->lhs_expr_.sizes();
    }

    size_t domainDim() const
    {
        return this->lhs_expr_.domainDim();
    }

private:
};
*/// Generic operator+ to create AddExpression instances using SFINAE
template <typename LhsExpr, typename RhsExpr>
typename std::enable_if<
    ExpressionTraits<LhsExpr>::order == ExpressionTraits<RhsExpr>::order &&
    ExpressionTraits<LhsExpr>::space == ExpressionTraits<RhsExpr>::space,
    AddExpression<LhsExpr, RhsExpr, ExpressionTraits<LhsExpr>::order, ExpressionTraits<RhsExpr>::order, ExpressionTraits<LhsExpr>::space, ExpressionTraits<RhsExpr>::space>
>::type
operator+(const BaseExpression<LhsExpr>& lhs, const BaseExpression<RhsExpr>& rhs)
{
    return AddExpression<LhsExpr, RhsExpr, ExpressionTraits<LhsExpr>::order, ExpressionTraits<RhsExpr>::order, ExpressionTraits<LhsExpr>::space, ExpressionTraits<RhsExpr>::space>(lhs, rhs);
}

// Specialization for transpose expressions (these should trigger compile errors)
template <typename LhsExpr, typename RhsExpr>
typename std::enable_if<
    ExpressionTraits<LhsExpr>::order == ExpressionTraits<RhsExpr>::order &&
    ExpressionTraits<LhsExpr>::space == ExpressionTraits<RhsExpr>::space,
    AddExpression<TransposeExpression<LhsExpr>, RhsExpr, ExpressionTraits<TransposeExpression<LhsExpr>>::order, ExpressionTraits<RhsExpr>::order, ExpressionTraits<TransposeExpression<LhsExpr>>::space, ExpressionTraits<RhsExpr>::space>
>::type
operator+(const TransposeExpression<LhsExpr>& lhs, const BaseExpression<RhsExpr>& rhs)
{
    GISMO_ERROR("Addition of a transposed vector and a vector is not defined.");
    return AddExpression<TransposeExpression<LhsExpr>, RhsExpr, ExpressionTraits<TransposeExpression<LhsExpr>>::order, ExpressionTraits<RhsExpr>::order, ExpressionTraits<TransposeExpression<LhsExpr>>::space, ExpressionTraits<RhsExpr>::space>(lhs, rhs);
}

template <typename LhsExpr, typename RhsExpr>
typename std::enable_if<
    ExpressionTraits<LhsExpr>::order == ExpressionTraits<RhsExpr>::order &&
    ExpressionTraits<LhsExpr>::space == ExpressionTraits<RhsExpr>::space,
    AddExpression<LhsExpr, TransposeExpression<RhsExpr>, ExpressionTraits<LhsExpr>::order, ExpressionTraits<TransposeExpression<RhsExpr>>::order, ExpressionTraits<LhsExpr>::space, ExpressionTraits<TransposeExpression<RhsExpr>>::space>
>::type
operator+(const BaseExpression<LhsExpr>& lhs, const TransposeExpression<RhsExpr>& rhs)
{
    GISMO_ERROR("Addition of a vector and a transposed vector is not defined.");
    return AddExpression<LhsExpr, TransposeExpression<RhsExpr>, ExpressionTraits<LhsExpr>::order, ExpressionTraits<TransposeExpression<RhsExpr>>::order, ExpressionTraits<LhsExpr>::space, ExpressionTraits<TransposeExpression<RhsExpr>>::space>(lhs, rhs);
}

// Specialization for scalars
template <typename LhsExpr>
typename std::enable_if<
    !std::is_same<LhsExpr, ConstantObject<typename ExpressionTraits<LhsExpr>::Scalar,0>>::value,
    AddExpression<LhsExpr, ConstantObject<typename ExpressionTraits<LhsExpr>::Scalar,0>, ExpressionTraits<LhsExpr>::order, 0, ExpressionTraits<LhsExpr>::space, Space::None>
>::type
operator+(const BaseExpression<LhsExpr>& lhs, const typename ExpressionTraits<LhsExpr>::Scalar rhs)
{
    return AddExpression<LhsExpr, ConstantObject<typename ExpressionTraits<LhsExpr>::Scalar,0>, ExpressionTraits<LhsExpr>::order, 0, ExpressionTraits<LhsExpr>::space, Space::None>(lhs, ConstantObject<typename ExpressionTraits<LhsExpr>::Scalar,0>(rhs));
}

template <typename RhsExpr>
typename std::enable_if<
    !std::is_same<RhsExpr, ConstantObject<typename ExpressionTraits<RhsExpr>::Scalar,0>>::value,
    AddExpression<ConstantObject<typename ExpressionTraits<RhsExpr>::Scalar,0>, RhsExpr, 0, ExpressionTraits<RhsExpr>::order, Space::None, ExpressionTraits<RhsExpr>::space>
>::type
operator+(const typename ExpressionTraits<RhsExpr>::Scalar lhs, const BaseExpression<RhsExpr>& rhs)
{
    return AddExpression<ConstantObject<typename ExpressionTraits<RhsExpr>::Scalar,0>, RhsExpr, 0, ExpressionTraits<RhsExpr>::order, Space::None, ExpressionTraits<RhsExpr>::space>(ConstantObject<typename ExpressionTraits<RhsExpr>::Scalar,0>(lhs), rhs);
}

// Specialization for matrices
template <typename LhsExpr>
typename std::enable_if<
    !std::is_same<LhsExpr, ConstantObject<typename ExpressionTraits<LhsExpr>::Scalar,ExpressionTraits<LhsExpr>::order>>::value,
    AddExpression<LhsExpr, ConstantObject<typename ExpressionTraits<LhsExpr>::Scalar,ExpressionTraits<LhsExpr>::order>, ExpressionTraits<LhsExpr>::order, ExpressionTraits<LhsExpr>::order, ExpressionTraits<LhsExpr>::space, ExpressionTraits<LhsExpr>::space>
>::type
operator+(const BaseExpression<LhsExpr>& lhs, const gsMatrix<typename ExpressionTraits<LhsExpr>::Scalar>& rhs)
{
    return AddExpression<LhsExpr, ConstantObject<typename ExpressionTraits<LhsExpr>::Scalar,ExpressionTraits<LhsExpr>::order>, ExpressionTraits<LhsExpr>::order, ExpressionTraits<LhsExpr>::order, ExpressionTraits<LhsExpr>::space, ExpressionTraits<LhsExpr>::space>(lhs, ConstantObject<typename ExpressionTraits<LhsExpr>::Scalar,ExpressionTraits<LhsExpr>::order>(rhs));
}

template <typename RhsExpr>
typename std::enable_if<
    !std::is_same<RhsExpr, ConstantObject<typename ExpressionTraits<RhsExpr>::Scalar,ExpressionTraits<RhsExpr>::order>>::value,
    AddExpression<ConstantObject<typename ExpressionTraits<RhsExpr>::Scalar,ExpressionTraits<RhsExpr>::order>, RhsExpr, ExpressionTraits<RhsExpr>::order, ExpressionTraits<RhsExpr>::order, ExpressionTraits<RhsExpr>::space, ExpressionTraits<RhsExpr>::space>
>::type
operator+(const gsMatrix<typename ExpressionTraits<RhsExpr>::Scalar>& lhs, const BaseExpression<RhsExpr>& rhs)
{
    return AddExpression<ConstantObject<typename ExpressionTraits<RhsExpr>::Scalar,ExpressionTraits<RhsExpr>::order>, RhsExpr, ExpressionTraits<RhsExpr>::order, ExpressionTraits<RhsExpr>::order, ExpressionTraits<RhsExpr>::space, ExpressionTraits<RhsExpr>::space>(ConstantObject<typename ExpressionTraits<RhsExpr>::Scalar,ExpressionTraits<RhsExpr>::order>(lhs), rhs);
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