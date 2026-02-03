/** @file ProductExpression.h

    @brief Product expression class for element-wise multiplication

    Implements element-wise (Hadamard) product of expressions. The product
    operation has special meaning for space types:
    - Test × Trial: Creates a bilinear form (e.g., mass matrix u*v)
    - None × None: Element-wise multiplication of fields
    - Test × None or Trial × None: Linear form

    Space type rules:
    - If either operand is Test or Trial, result inherits that space
    - Test × Trial → checks for consistency, used in bilinear forms
    - None × None → None

    Order behavior:
    - Product of Order 0 × Order 0 → Order 0 (scalar × scalar → scalar)
    - Product of Order 1 × Order 0 → Order 1 (vector × scalar → vector)
    - Product of Order 0 × Order 1 → Order 1 (scalar × vector → vector)
    - Order 1 × Order 1 requires inner product (use dot() instead)

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

// --- Generic ExpressionTraits for ProductExpression ---
// Single unified traits template that handles all valid combinations
template <typename _LhsExpr, typename _RhsExpr, size_t _LhsOrder, size_t _RhsOrder, enum SpaceType _LhsSpace, enum SpaceType _RhsSpace>
struct ExpressionTraits<ProductExpression<_LhsExpr, _RhsExpr, _LhsOrder, _RhsOrder, _LhsSpace, _RhsSpace>>
{
    typedef _LhsExpr LhsType;
    typedef _RhsExpr RhsType;

    typedef typename ExpressionTraits<_LhsExpr>::Scalar Scalar;
    
    // Product order logic: depends on the tensor operation type
    // For scalar * anything: keep the other order
    // For matrix-vector: 2*1 -> 1, 1*2 -> 1 (dot products)
    // For matrix-matrix: 2*2 -> 2
    static constexpr size_t Order = (_LhsOrder == 0) ? _RhsOrder : 
                                    (_RhsOrder == 0) ? _LhsOrder :
                                    ((_LhsOrder == 2 && _RhsOrder == 1) ? 1 :
                                    ((_LhsOrder == 1 && _RhsOrder == 2) ? 1 :
                                    ((_LhsOrder == 2 && _RhsOrder == 2) ? 2 : 0)));
    
    // Space logic: prefer non-None space; if both different non-None, use Both
    static constexpr SpaceType Space = (_LhsSpace != SpaceType::None && _RhsSpace != SpaceType::None && _LhsSpace != _RhsSpace) 
                                   ? SpaceType::Both
                                   : ((_LhsSpace != SpaceType::None) ? _LhsSpace : _RhsSpace);
    
    static constexpr size_t Deriv = (ExpressionTraits<_LhsExpr>::Deriv > ExpressionTraits<_RhsExpr>::Deriv)
                                   ? ExpressionTraits<_LhsExpr>::Deriv : ExpressionTraits<_RhsExpr>::Deriv;
    static constexpr bool IsConstant = ExpressionTraits<_LhsExpr>::IsConstant && ExpressionTraits<_RhsExpr>::IsConstant;
};

// --- Unified ProductExpression class with enable_if specializations for eval() ---
template <typename _LhsExpr, typename _RhsExpr, size_t _LhsOrder, size_t _RhsOrder, enum SpaceType _LhsSpace, enum SpaceType _RhsSpace>
class ProductExpression : public BinaryOperator<ProductExpression<_LhsExpr, _RhsExpr, _LhsOrder, _RhsOrder, _LhsSpace, _RhsSpace>>
{
    using Base = BinaryOperator<ProductExpression<_LhsExpr, _RhsExpr, _LhsOrder, _RhsOrder, _LhsSpace, _RhsSpace>>;

    static_assert(_LhsSpace == _RhsSpace || 
                  _LhsSpace == SpaceType::None || 
                  _RhsSpace == SpaceType::None ||
                  (_LhsSpace == SpaceType::Test && _RhsSpace == SpaceType::Trial) ||
                  (_LhsSpace == SpaceType::Trial && _RhsSpace == SpaceType::Test),
                  "ProductExpression requires compatible spaces: same space, one None, or Test×Trial combination.");
    
    // Additional check: Vector dot product (Order 1 * Order 1) with Test and/or Trial spaces is not allowed
    // because the standard * operator for vectors would compute dot product (collapse to scalar)
    // and lose the bilinear form structure. Use inner() explicitly for bilinear forms.
    static_assert(!(_LhsOrder == 1 && _RhsOrder == 1 && 
                    ((_LhsSpace != SpaceType::None) || (_RhsSpace != SpaceType::None))),
                  "Vector-vector product (e.g., grad(u) * grad(v)) with Test/Trial spaces is ambiguous and not allowed. "
                  "Use inner(lhs, rhs) to create a bilinear form, or explicitly use .transpose() for outer products.");

public:
    ProductExpression(const _LhsExpr& lhs, const _RhsExpr& rhs)
        : Base(lhs, rhs)
    {
        // Space consistency checks - conditions are compile-time constants, 
        // compiler will optimize away dead branches
        if (_LhsSpace == _RhsSpace && _LhsSpace == SpaceType::Trial)
            GISMO_ASSERT(lhs.trial()==rhs.trial(),"ProductExpression requires both operands to have the same trial space.");
        if (_LhsSpace == _RhsSpace && _LhsSpace == SpaceType::Test)
            GISMO_ASSERT(lhs.test()==rhs.test(),"ProductExpression requires both operands to have the same test space.");
        if (_LhsSpace == SpaceType::Both && _RhsSpace == SpaceType::Both)
        {
            GISMO_ASSERT(lhs.test()==rhs.test(),"ProductExpression requires both operands to have the same test space.");
            GISMO_ASSERT(lhs.trial()==rhs.trial(),"ProductExpression requires both operands to have the same trial space.");
        }

        // Set sizes based on the product operation type
        if (_LhsOrder == 0 && _RhsOrder > 0) {
            // Scalar * Expression: inherit sizes from RHS
            for (size_t d = 0; d != _RhsOrder; ++d) {
                this->sizes_[d] = this->rhs_expr_.sizes()[d];
            }
        } else if (_LhsOrder > 0 && _RhsOrder == 0) {
            // Expression * Scalar: inherit sizes from LHS
            for (size_t d = 0; d != _LhsOrder; ++d) {
                this->sizes_[d] = this->lhs_expr_.sizes()[d];
            }
        } else if (_LhsOrder == 2 && _RhsOrder == 2) {
            // Matrix * Matrix: (m,n) * (n,p) -> (m,p)
            GISMO_ENSURE(this->lhs_expr_.sizes()[1] == this->rhs_expr_.sizes()[0], 
                        "ProductExpression: Size mismatch in matrix-matrix product");
            this->sizes_[0] = this->lhs_expr_.sizes()[0];
            this->sizes_[1] = this->rhs_expr_.sizes()[1];
        } else if (_LhsOrder == 2 && _RhsOrder == 1) {
            // Matrix * Vector: (m,n) * (n,) -> (m,)
            GISMO_ENSURE(this->lhs_expr_.sizes()[1] == this->rhs_expr_.sizes()[0], 
                        "ProductExpression: Size mismatch in matrix-vector product");
            this->sizes_[0] = this->lhs_expr_.sizes()[0];
        } else if (_LhsOrder == 1 && _RhsOrder == 2) {
            // Vector * Matrix: handled specially (transpose involved)
            GISMO_ENSURE(this->lhs_expr_.sizes()[0] == this->rhs_expr_.sizes()[0], 
                        "ProductExpression: Size mismatch in vector-matrix product");
            this->sizes_[0] = this->rhs_expr_.sizes()[1];
        } else if (_LhsOrder == 0 && _RhsOrder == 0) {
            // Scalar * Scalar: no sizes needed
        } else if (_LhsOrder == 1 && _RhsOrder == 1) {
            // Vector . Vector (dot product): becomes scalar, no sizes
        }
    }

    void print(std::ostream & os) const
    {
        os << "(" << this->lhs() << " * " << this->rhs() << ")";
    }

    size_t domainDim() const
    {
        return (_LhsOrder >= _RhsOrder) ? this->lhs_expr_.domainDim() : this->rhs_expr_.domainDim();
    }

    // --- eval() specialization 0: Fallback/default for unhandled combinations ---
    // This handles edge cases not covered by the specific specializations
    template <size_t LO = _LhsOrder, size_t RO = _RhsOrder, size_t LS = _LhsSpace, size_t RS = _RhsSpace>
    typename std::enable_if<!((LO == 0 && RO > 0) || (LO > 0 && RO == 0)),
                            ExpressionResult<typename Base::Scalar>>::type
    eval(const index_t k) const
    {
        return this->lhs_expr_.eval(k) * this->rhs_expr_.eval(k);
    }

    // --- eval() specialization 2: Scalar * Expression (0, N, *, None) ---
    template <size_t LO = _LhsOrder, size_t RO = _RhsOrder>
    typename std::enable_if<(LO == 0) && (RO > 0),
                            ExpressionResult<typename Base::Scalar>>::type
    eval(const index_t k) const
    {
        return ScalarExpressionResult<typename Base::Scalar>(this->lhs_expr_.eval(k)) * this->rhs_expr_.eval(k);
    }

    // --- eval() specialization 4: Expression * Scalar (N, 0, Space, *) ---
    template <size_t LO = _LhsOrder, size_t RO = _RhsOrder>
    typename std::enable_if<(LO > 0) && (RO == 0),
                            ExpressionResult<typename Base::Scalar>>::type
    eval(const index_t k) const
    {
        return this->lhs_expr_.eval(k) * ScalarExpressionResult<typename Base::Scalar>(this->rhs_expr_.eval(k));
    }

private:
    // Helper to get trial space - returns actual trial if available, otherwise from operand's trial()
    template <SpaceType S1, SpaceType S2>
    auto get_trial_impl(TrueTag /* lhs has trial */, FalseTag) const 
        -> decltype(this->lhs_expr_.trial())
    {
        return this->lhs_expr_.trial();
    }
    
    template <SpaceType S1, SpaceType S2>
    auto get_trial_impl(FalseTag, TrueTag /* rhs has trial */) const 
        -> decltype(this->rhs_expr_.trial())
    {
        return this->rhs_expr_.trial();
    }
    
    template <SpaceType S1, SpaceType S2>
    auto get_trial_impl(FalseTag, FalseTag) const 
        -> decltype(this->lhs_expr_.trial())
    {
        // Neither has trial, but we still need to return something - delegate to lhs which will return NullObject
        return this->lhs_expr_.trial();
    }

    // Helper to get test space - returns actual test if available, otherwise from operand's test()
    template <SpaceType S1, SpaceType S2>
    auto get_test_impl(TrueTag /* lhs has test */, FalseTag) const 
        -> decltype(this->lhs_expr_.test())
    {
        return this->lhs_expr_.test();
    }
    
    template <SpaceType S1, SpaceType S2>
    auto get_test_impl(FalseTag, TrueTag /* rhs has test */) const 
        -> decltype(this->rhs_expr_.test())
    {
        return this->rhs_expr_.test();
    }
    
    template <SpaceType S1, SpaceType S2>
    auto get_test_impl(FalseTag, FalseTag) const 
        -> decltype(this->lhs_expr_.test())
    {
        // Neither has test, but we still need to return something - delegate to lhs which will return NullObject
        return this->lhs_expr_.test();
    }

public:
    auto trial() const -> decltype(this->get_trial_impl<_LhsSpace, _RhsSpace>(
        typename BoolToTag<_LhsSpace == SpaceType::Trial || _LhsSpace == SpaceType::Both>::type(),
        typename BoolToTag<_RhsSpace == SpaceType::Trial || _RhsSpace == SpaceType::Both>::type()))
    {
        return this->get_trial_impl<_LhsSpace, _RhsSpace>(
            typename BoolToTag<_LhsSpace == SpaceType::Trial || _LhsSpace == SpaceType::Both>::type(),
            typename BoolToTag<_RhsSpace == SpaceType::Trial || _RhsSpace == SpaceType::Both>::type());
    }

    auto test() const -> decltype(this->get_test_impl<_LhsSpace, _RhsSpace>(
        typename BoolToTag<_LhsSpace == SpaceType::Test || _LhsSpace == SpaceType::Both>::type(),
        typename BoolToTag<_RhsSpace == SpaceType::Test || _RhsSpace == SpaceType::Both>::type()))
    {
        return this->get_test_impl<_LhsSpace, _RhsSpace>(
            typename BoolToTag<_LhsSpace == SpaceType::Test || _LhsSpace == SpaceType::Both>::type(),
            typename BoolToTag<_RhsSpace == SpaceType::Test || _RhsSpace == SpaceType::Both>::type());
    }

};

// Generic operator* to create ProductExpression instances using SFINAE
template <typename _LhsExpr, typename _RhsExpr>
typename std::enable_if<
    _LhsExpr::Order == 0 && _RhsExpr::Order == 0,
    ProductExpression<_LhsExpr, _RhsExpr, 0, 0, _LhsExpr::Space, _RhsExpr::Space>
>::type
operator*(const BaseExpression<_LhsExpr>& lhs, const BaseExpression<_RhsExpr>& rhs)
{
    return ProductExpression<_LhsExpr, _RhsExpr, 0, 0, _LhsExpr::Space, _RhsExpr::Space>(lhs, rhs);
}

// Specialization: Scalar Expression × Non-scalar Expression
template <typename _LhsExpr, typename _RhsExpr>
typename std::enable_if<
    _LhsExpr::Order == 0 && _RhsExpr::Order != 0,
    ProductExpression<_LhsExpr, _RhsExpr, 0, _RhsExpr::Order, _LhsExpr::Space, _RhsExpr::Space>
>::type
operator*(const BaseExpression<_LhsExpr>& lhs, const BaseExpression<_RhsExpr>& rhs)
{
    return ProductExpression<_LhsExpr, _RhsExpr, 0, _RhsExpr::Order, _LhsExpr::Space, _RhsExpr::Space>(lhs, rhs);
}

// Specialization: Expression × Scalar Expression
template <typename _LhsExpr, typename _RhsExpr>
typename std::enable_if<
    _LhsExpr::Order != 0 && _RhsExpr::Order == 0,
    ProductExpression<_RhsExpr, _LhsExpr, 0, _LhsExpr::Order, _RhsExpr::Space, _LhsExpr::Space>
>::type
operator*(const BaseExpression<_LhsExpr>& lhs, const BaseExpression<_RhsExpr>& rhs)
{
    return ProductExpression<_RhsExpr, _LhsExpr, 0, _LhsExpr::Order, _RhsExpr::Space, _LhsExpr::Space>(rhs, lhs);
}

// Specialization: Expression × Scalar primitive
template <typename _LhsExpr>
ProductExpression<ConstantObject<typename ExpressionTraits<_LhsExpr>::Scalar,0>, _LhsExpr, 0, _LhsExpr::Order, SpaceType::None, _LhsExpr::Space>
operator*(const BaseExpression<_LhsExpr>& lhs, const typename ExpressionTraits<_LhsExpr>::Scalar rhs)
{
    return ProductExpression<ConstantObject<typename ExpressionTraits<_LhsExpr>::Scalar,0>, _LhsExpr, 0, _LhsExpr::Order, SpaceType::None, _LhsExpr::Space>(ConstantObject<typename ExpressionTraits<_LhsExpr>::Scalar,0>(rhs), lhs);
}

// Specialization: Scalar primitive × Expression
template <typename _RhsExpr>
ProductExpression<ConstantObject<typename ExpressionTraits<_RhsExpr>::Scalar,0>, _RhsExpr, 0, _RhsExpr::Order, SpaceType::None, _RhsExpr::Space>
operator*(const typename ExpressionTraits<_RhsExpr>::Scalar lhs, const BaseExpression<_RhsExpr>& rhs)
{
    return ProductExpression<ConstantObject<typename ExpressionTraits<_RhsExpr>::Scalar,0>, _RhsExpr, 0, _RhsExpr::Order, SpaceType::None, _RhsExpr::Space>(ConstantObject<typename ExpressionTraits<_RhsExpr>::Scalar,0>(lhs), rhs);
}

// Specialization for vector-vector dot product (Order 1 * Order 1 = Order 0)
template <typename _LhsExpr, typename _RhsExpr>
typename std::enable_if<
    _LhsExpr::Order == 1 && _RhsExpr::Order == 1,
    ProductExpression<_LhsExpr, _RhsExpr, 1, 1, _LhsExpr::Space, _RhsExpr::Space>
>::type
operator*(const BaseExpression<_LhsExpr>& lhs, const BaseExpression<_RhsExpr>& rhs)
{
    return ProductExpression<_LhsExpr, _RhsExpr, 1, 1, _LhsExpr::Space, _RhsExpr::Space>(lhs, rhs);
}

// Specialization for matrix-vector product
template <typename _LhsExpr, typename _RhsExpr>
typename std::enable_if<
    _LhsExpr::Order == 2 && _RhsExpr::Order == 1,
    ProductExpression<_LhsExpr, _RhsExpr, 2, 1, _LhsExpr::Space, _RhsExpr::Space>
>::type
operator*(const BaseExpression<_LhsExpr>& lhs, const BaseExpression<_RhsExpr>& rhs)
{
    return ProductExpression<_LhsExpr, _RhsExpr, 2, 1, _LhsExpr::Space, _RhsExpr::Space>(lhs, rhs);
}

// Specialization for vector-matrix product (row vector × matrix = row vector)
template <typename _LhsExpr, typename _RhsExpr>
typename std::enable_if<
    _LhsExpr::Order == 1 && _RhsExpr::Order == 2,
    ProductExpression<_LhsExpr, _RhsExpr, 1, 2, _LhsExpr::Space, _RhsExpr::Space>
>::type
operator*(const BaseExpression<_LhsExpr>& lhs, const BaseExpression<_RhsExpr>& rhs)
{
    return ProductExpression<_LhsExpr, _RhsExpr, 1, 2, _LhsExpr::Space, _RhsExpr::Space>(lhs, rhs);
}

// Specialization for matrix-vector product (vector primitive)
template <typename _LhsExpr, int _Rows>
typename std::enable_if<
    _LhsExpr::Order == 2,
    ProductExpression<_LhsExpr, ConstantObject<typename ExpressionTraits<_LhsExpr>::Scalar,1>, 2, 1, _LhsExpr::Space, SpaceType::None>
>::type
operator*(const BaseExpression<_LhsExpr>& lhs, const gsVector<typename ExpressionTraits<_LhsExpr>::Scalar, _Rows>& rhs)
{
    return ProductExpression<_LhsExpr, ConstantObject<typename ExpressionTraits<_LhsExpr>::Scalar,1>, 2, 1, _LhsExpr::Space, SpaceType::None>(lhs, ConstantObject<typename ExpressionTraits<_LhsExpr>::Scalar,1>(rhs));
}

// Specialization for matrix-vector product (matrix primitive)
template <typename _RhsExpr, int _Rows, int _Cols>
typename std::enable_if<
    _RhsExpr::Order == 1,
    ProductExpression<ConstantObject<typename ExpressionTraits<_RhsExpr>::Scalar,2>, _RhsExpr, 2, 1, SpaceType::None, _RhsExpr::Space>
>::type
operator*(const gsMatrix<typename ExpressionTraits<_RhsExpr>::Scalar, _Rows, _Cols>& lhs, const BaseExpression<_RhsExpr>& rhs)
{
    return ProductExpression<ConstantObject<typename ExpressionTraits<_RhsExpr>::Scalar,2>, _RhsExpr, 2, 1, SpaceType::None, _RhsExpr::Space>(ConstantObject<typename ExpressionTraits<_RhsExpr>::Scalar,2>(lhs), rhs);
}

// Specialization for vector-matrix product
template <typename _LhsExpr, typename _RhsExpr>
typename std::enable_if<
    _LhsExpr::Order == 1 && _RhsExpr::Order == 2,
    TransposeExpression<ProductExpression<_RhsExpr, _LhsExpr, 2, 1, _RhsExpr::Space, _LhsExpr::Space>>
>::type
operator*(const TransposeExpression<_LhsExpr>& lhs, const BaseExpression<_RhsExpr>& rhs)
{
    auto product = ProductExpression<_RhsExpr, _LhsExpr, 2, 1, _RhsExpr::Space, _LhsExpr::Space>(rhs, lhs.expr());
    return TransposeExpression<ProductExpression<_RhsExpr, _LhsExpr, 2, 1, _RhsExpr::Space, _LhsExpr::Space>>(product);

}

// Specialization for vector-matrix product (matrix primitive)
template <typename _LhsExpr, int _Rows, int _Cols>
typename std::enable_if<
    _LhsExpr::Order == 1,
    ProductExpression<ConstantObject<typename ExpressionTraits<_LhsExpr>::Scalar,2>, _LhsExpr, 2, 1, SpaceType::None, _LhsExpr::Space>
>::type
operator*(const TransposeExpression<_LhsExpr>& lhs, const gsMatrix<typename ExpressionTraits<_LhsExpr>::Scalar, _Rows, _Cols>& rhs)
{
    return ProductExpression<ConstantObject<typename ExpressionTraits<_LhsExpr>::Scalar,2>, _LhsExpr, 2, 1, SpaceType::None, _LhsExpr::Space>(ConstantObject<typename ExpressionTraits<_LhsExpr>::Scalar,2>(rhs), lhs.expr());
}

// Specialization for vector-matrix product (vector primitive)
template <typename _RhsExpr, int _Rows>
typename std::enable_if<
    _RhsExpr::Order == 2,
    ProductExpression<_RhsExpr, TransposeExpression<ConstantObject<typename ExpressionTraits<_RhsExpr>::Scalar,1>>, 2, 1, _RhsExpr::Space, SpaceType::None>
>::type
operator*(const gsEigen::Transpose<gsVector<typename ExpressionTraits<_RhsExpr>::Scalar, _Rows>>& lhs, const BaseExpression<_RhsExpr>& rhs)
{
    return ProductExpression<_RhsExpr, TransposeExpression<ConstantObject<typename ExpressionTraits<_RhsExpr>::Scalar,1>>, 2, 1, _RhsExpr::Space, SpaceType::None>(rhs, TransposeExpression<ConstantObject<typename ExpressionTraits<_RhsExpr>::Scalar,1>>(ConstantObject<typename ExpressionTraits<_RhsExpr>::Scalar,1>(lhs)));
}   

// Specialization for transpose-vector dot product: (vector)^T * vector = scalar
template <typename _LhsExpr, typename _RhsExpr>
typename std::enable_if<
    _LhsExpr::Order == 1 && _RhsExpr::Order == 1,
    InnerProductExpression<_LhsExpr, _RhsExpr, 1, 1, _LhsExpr::Space, _RhsExpr::Space>
>::type
operator*(const TransposeExpression<_LhsExpr>& lhs, const BaseExpression<_RhsExpr>& rhs)
{
    return InnerProductExpression<_LhsExpr, _RhsExpr, 1, 1, _LhsExpr::Space, _RhsExpr::Space>(lhs.expr(), rhs);
}

// Specialization for vector-matrix dot product: vector * (vector)^T = matrix
template <typename _LhsExpr, typename _RhsExpr>
typename std::enable_if<
    _LhsExpr::Order == 1 && _RhsExpr::Order == 1,
    OuterProductExpression<_LhsExpr, _RhsExpr, 2, 2, _LhsExpr::Space, _RhsExpr::Space>
>::type
operator*(const BaseExpression<_LhsExpr>& lhs, const TransposeExpression<_RhsExpr>&rhs)
{
    return OuterProductExpression<_LhsExpr, _RhsExpr, 2, 2, _LhsExpr::Space, _RhsExpr::Space>(lhs, rhs.expr());
}

// Matrix-matrix product
template <typename _LhsExpr, typename _RhsExpr>
typename std::enable_if<
    _LhsExpr::Order == 2 && _RhsExpr::Order == 2,
    ProductExpression<_LhsExpr, _RhsExpr, 2, 2, _LhsExpr::Space, _RhsExpr::Space>
>::type
operator*(const BaseExpression<_LhsExpr>& lhs, const BaseExpression<_RhsExpr>& rhs)
{
    return ProductExpression<_LhsExpr, _RhsExpr, 2, 2, _LhsExpr::Space, _RhsExpr::Space>(lhs, rhs);
}

// Specialization for matrix-matrix product (first argument primitive)
template <typename _RhsExpr, int _Rows, int _Cols>
typename std::enable_if<
    _RhsExpr::Order == 2,
    ProductExpression<ConstantObject<typename ExpressionTraits<_RhsExpr>::Scalar,2>, _RhsExpr, 2, 2, SpaceType::None, _RhsExpr::Space>
>::type
operator*(const gsMatrix<typename ExpressionTraits<_RhsExpr>::Scalar, _Rows, _Cols>& lhs, const BaseExpression<_RhsExpr>& rhs)
{
    return ProductExpression<ConstantObject<typename ExpressionTraits<_RhsExpr>::Scalar,2>, _RhsExpr, 2, 2, SpaceType::None, _RhsExpr::Space>(ConstantObject<typename ExpressionTraits<_RhsExpr>::Scalar,2>(lhs), rhs);
}

// Specialization for matrix-matrix product (second argument primitive)
template <typename _LhsExpr, int _Rows, int _Cols>
typename std::enable_if<
    _LhsExpr::Order == 2,
    ProductExpression<_LhsExpr, ConstantObject<typename ExpressionTraits<_LhsExpr>::Scalar,2>, 2, 2, _LhsExpr::Space, SpaceType::None>
>::type
operator*(const BaseExpression<_LhsExpr>& lhs, const  gsMatrix<typename ExpressionTraits<_LhsExpr>::Scalar, _Rows, _Cols>& rhs)
{
    return ProductExpression<_LhsExpr, ConstantObject<typename ExpressionTraits<_LhsExpr>::Scalar,2>, 2, 2, _LhsExpr::Space, SpaceType::None>(lhs, ConstantObject<typename ExpressionTraits<_LhsExpr>::Scalar,2>(rhs));
}

// Specialization for NullObject
template <typename _T, enum SpaceType _LhsSpace, size_t _LhsOrder, typename _RhsExpr>
auto operator*(const NullObject<_T,_LhsSpace,_LhsOrder>& /* lhs */, const BaseExpression<_RhsExpr>& /* rhs */)
-> NullObject<_T,SpaceType::None,0>
{
    return NullObject<_T,SpaceType::None,0>::get();
}

template <typename _LhsExpr, typename _T, enum SpaceType _RhsSpace, size_t _RhsOrder>
auto operator*(const BaseExpression<_LhsExpr>& /* lhs */, const NullObject<_T,_RhsSpace,_RhsOrder>& /* rhs */)
-> NullObject<_T,SpaceType::None,0>
{
    return NullObject<_T,SpaceType::None,0>::get();
}

// Specialization: operator-
template <typename _Expr>
auto operator-(const BaseExpression<_Expr>& expr)
-> decltype(-1.0*expr)
{
    return -1.0*expr;
}

// Variation of ProductExpression where LHS is a SolutionObject and RhsExpr is NOT a SolutionObject
template <typename _LhsExpr, typename _RhsExpr, typename _SpaceObject>
typename std::enable_if<
    !std::is_base_of<SolutionObject<typename _LhsExpr::Scalar,_LhsExpr::Space,_LhsExpr::Order>, _RhsExpr>::value,
    decltype(variation(std::declval<_RhsExpr>(), std::declval<_SpaceObject>()))
>::type
variation(const ProductExpression<SolutionObject<typename _LhsExpr::Scalar,_LhsExpr::Space,_LhsExpr::Order>,_RhsExpr,_LhsExpr::Order,_RhsExpr::Order,SolutionObject<typename _LhsExpr::Scalar,_LhsExpr::Space,_LhsExpr::Order>::Order,_RhsExpr::Space> & expr,
          const _SpaceObject & space)
{
    return expr.lhs()*variation(expr.rhs(), space);
}

// Variation of ProductExpression where RHS is a SolutionObject and LhsExpr is NOT a SolutionObject
template <typename _LhsExpr, typename _RhsExpr, typename _SpaceObject>
typename std::enable_if<
    !std::is_base_of<SolutionObject<typename _RhsExpr::Scalar,_RhsExpr::Space,_RhsExpr::Order>, _LhsExpr>::value,
    decltype(variation(std::declval<_LhsExpr>(), std::declval<_SpaceObject>()))
>::type
variation(const ProductExpression<_LhsExpr,SolutionObject<typename _RhsExpr::Scalar,_RhsExpr::Space,_RhsExpr::Order>,_LhsExpr::Order,_RhsExpr::Order,_LhsExpr::Space,SolutionObject<typename _RhsExpr::Scalar,_RhsExpr::Space,_RhsExpr::Order>::Order> & expr,
          const _SpaceObject & space)
{
    return variation(expr.lhs(), space)*expr.rhs();
}

template<typename _LhsExpr, typename _RhsExpr, typename _SpaceObject>
auto variation(const ProductExpression<_LhsExpr,_RhsExpr,_LhsExpr::Order,_RhsExpr::Order,_LhsExpr::Space,_RhsExpr::Space> & expr,
               const _SpaceObject & space)
-> decltype(variation(expr.lhs(), space)*expr.rhs() + expr.lhs()*variation(expr.rhs(), space))
{
    return variation(expr.lhs(), space)*expr.rhs() + expr.lhs()*variation(expr.rhs(), space);
}

}//namespace Expr
}//namespace gismo
