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

// --- Generic ExpressionTraits for ProductExpression ---
// Single unified traits template that handles all valid combinations
template <typename _LhsExpr, typename _RhsExpr, size_t _LhsOrder, size_t _RhsOrder, size_t _LhsSpace, size_t _RhsSpace>
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
    static constexpr size_t Space = (_LhsSpace != SpaceType::None && _RhsSpace != SpaceType::None && _LhsSpace != _RhsSpace) 
                                   ? static_cast<size_t>(SpaceType::Both)
                                   : ((_LhsSpace != SpaceType::None) ? _LhsSpace : _RhsSpace);
    
    static constexpr size_t Deriv = (ExpressionTraits<_LhsExpr>::Deriv > ExpressionTraits<_RhsExpr>::Deriv)
                                   ? ExpressionTraits<_LhsExpr>::Deriv : ExpressionTraits<_RhsExpr>::Deriv;
    static constexpr bool IsConstant = ExpressionTraits<_LhsExpr>::IsConstant && ExpressionTraits<_RhsExpr>::IsConstant;
};

// --- Unified ProductExpression class with enable_if specializations for eval() ---
template <typename _LhsExpr, typename _RhsExpr, size_t _LhsOrder, size_t _RhsOrder, size_t _LhsSpace, size_t _RhsSpace>
class ProductExpression : public BinaryOperator<ProductExpression<_LhsExpr, _RhsExpr, _LhsOrder, _RhsOrder, _LhsSpace, _RhsSpace>>
{
    using Base = BinaryOperator<ProductExpression<_LhsExpr, _RhsExpr, _LhsOrder, _RhsOrder, _LhsSpace, _RhsSpace>>;

public:
    ProductExpression(const _LhsExpr& lhs, const _RhsExpr& rhs)
        : Base(lhs, rhs)
    {
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

    size_t domainDim() const
    {
        return (_LhsOrder >= _RhsOrder) ? this->lhs_expr_.domainDim() : this->rhs_expr_.domainDim();
    }

    // --- eval() specialization 0: Fallback/default for unhandled combinations ---
    // This handles edge cases not covered by the specific specializations
    template <size_t LO = _LhsOrder, size_t RO = _RhsOrder, size_t LS = _LhsSpace, size_t RS = _RhsSpace>
    typename std::enable_if<!((LO == 0 && RO == 0 && LS == SpaceType::None && RS == SpaceType::None) ||
                              (LO == 0 && RO > 0 && RS == SpaceType::None) ||
                              (LO == 0 && RO > 0 && RS != SpaceType::None) ||
                              (LO > 0 && RO == 0) ||
                              (LO == 2 && RO == 2 && LS == SpaceType::None && RS == SpaceType::None) ||
                              (LO == 2 && RO == 2 && LS != SpaceType::None && RS != SpaceType::None && LS == RS) ||
                              (LO == 2 && RO == 1 && LS == SpaceType::None && RS == SpaceType::None) ||
                              (LO == 2 && RO == 1 && LS != SpaceType::None && RS != SpaceType::None && LS == RS)),
                            ExpressionValue<typename Base::Scalar>>::type
    eval(const index_t k) const
    {
        gsDebug<<"ProductExpression: eval fallback (handling edge case) at k="<<k<<"\n";
        ExpressionValue<typename Base::Scalar> lhs_val = this->lhs_expr_.eval(k);
        ExpressionValue<typename Base::Scalar> rhs_val = this->rhs_expr_.eval(k);
        
        // For (0,0) with different spaces or other edge cases, treat as scalar multiplication
        if (lhs_val.rowCardinality() == 1 && lhs_val.colCardinality() == 1 &&
            rhs_val.rowCardinality() == 1 && rhs_val.colCardinality() == 1)
        {
            ExpressionValue<typename Base::Scalar> result(1, 1);
            result(0, 0) = lhs_val(0, 0) * rhs_val(0, 0);
            return result;
        }
        
        // General case: elementwise multiplication with broadcasting if needed
        ExpressionValue<typename Base::Scalar> result(rhs_val.rowCardinality(), rhs_val.colCardinality());
        
        if (lhs_val.rowCardinality() == 1 && lhs_val.colCardinality() == 1)
        {
            // LHS is scalar, broadcast to RHS size
            typename Base::Scalar scalar_value = lhs_val(0, 0)(0, 0);
            for (index_t i = 0; i < result.rowCardinality(); ++i)
            {
                for (index_t j = 0; j < result.colCardinality(); ++j)
                {
                    result(i, j) = rhs_val(i, j).array() * scalar_value;
                }
            }
        }
        else if (rhs_val.rowCardinality() == 1 && rhs_val.colCardinality() == 1)
        {
            // RHS is scalar, broadcast to LHS size
            typename Base::Scalar scalar_value = rhs_val(0, 0)(0, 0);
            result.resize(lhs_val.rowCardinality(), lhs_val.colCardinality());
            for (index_t i = 0; i < result.rowCardinality(); ++i)
            {
                for (index_t j = 0; j < result.colCardinality(); ++j)
                {
                    result(i, j) = lhs_val(i, j).array() * scalar_value;
                }
            }
        }
        else
        {
            // Both non-scalar: elementwise multiplication
            result.resize(lhs_val.rowCardinality(), lhs_val.colCardinality());
            for (index_t i = 0; i < result.rowCardinality(); ++i)
            {
                for (index_t j = 0; j < result.colCardinality(); ++j)
                {
                    result(i, j) = lhs_val(i, j) * rhs_val(i, j);
                }
            }
        }
        
        return result;
    }

    // --- eval() specialization 1: Scalar (0,0) * Scalar (0,0) ---
    template <size_t LO = _LhsOrder, size_t RO = _RhsOrder, size_t LS = _LhsSpace, size_t RS = _RhsSpace>
    typename std::enable_if<LO == 0 && RO == 0 && LS == SpaceType::None && RS == SpaceType::None,
                            ExpressionValue<typename Base::Scalar>>::type
    eval(const index_t k) const
    {
        gsDebug<<"ProductExpression: eval Scalar * Scalar at k="<<k<<"\n";
        ExpressionValue<typename Base::Scalar> lhs_val = this->lhs_expr_.eval(k);
        ExpressionValue<typename Base::Scalar> rhs_val = this->rhs_expr_.eval(k);
        
        ExpressionValue<typename Base::Scalar> result(1, 1);
        result(0, 0) = lhs_val(0, 0) * rhs_val(0, 0);
        return result;
    }

    // --- eval() specialization 2: Scalar * Expression (0, N, *, None) ---
    template <size_t LO = _LhsOrder, size_t RO = _RhsOrder, size_t LS = _LhsSpace, size_t RS = _RhsSpace>
    typename std::enable_if<(LO == 0) && (RO > 0) && (RS == SpaceType::None),
                            ExpressionValue<typename Base::Scalar>>::type
    eval(const index_t k) const
    {
        gsDebug<<"ProductExpression: eval Scalar * Expression (None space) at k="<<k<<"\n";
        ExpressionValue<typename Base::Scalar> lhs_val = this->lhs_expr_.eval(k);
        ExpressionValue<typename Base::Scalar> rhs_val = this->rhs_expr_.eval(k);
        
        typename Base::Scalar scalar_value = lhs_val(0, 0)(0, 0);
        ExpressionValue<typename Base::Scalar> result(1, 1);
        result(0, 0) = rhs_val(0, 0).array() * scalar_value;
        return result;
    }

    // --- eval() specialization 3: Scalar * Expression (0, N, *, Space) ---
    template <size_t LO = _LhsOrder, size_t RO = _RhsOrder, size_t LS = _LhsSpace, size_t RS = _RhsSpace>
    typename std::enable_if<(LO == 0) && (RO > 0) && (RS != SpaceType::None),
                            ExpressionValue<typename Base::Scalar>>::type
    eval(const index_t k) const
    {
        gsDebug<<"ProductExpression: eval Scalar * Expression (space-dependent) at k="<<k<<"\n";
        ExpressionValue<typename Base::Scalar> lhs_val = this->lhs_expr_.eval(k);
        ExpressionValue<typename Base::Scalar> rhs_val = this->rhs_expr_.eval(k);
        
        typename Base::Scalar scalar_value = lhs_val(0, 0)(0, 0);
        ExpressionValue<typename Base::Scalar> result(rhs_val.rowCardinality(), rhs_val.colCardinality());
        
        for (index_t i = 0; i < result.rowCardinality(); ++i)
        {
            for (index_t j = 0; j < result.colCardinality(); ++j)
            {
                result(i, j) = rhs_val(i, j).array() * scalar_value;
            }
        }
        
        return result;
    }

    // --- eval() specialization 4: Expression * Scalar (N, 0, Space, *) ---
    template <size_t LO = _LhsOrder, size_t RO = _RhsOrder, size_t LS = _LhsSpace, size_t RS = _RhsSpace>
    typename std::enable_if<(LO > 0) && (RO == 0),
                            ExpressionValue<typename Base::Scalar>>::type
    eval(const index_t k) const
    {
        gsDebug<<"ProductExpression: eval Expression * Scalar at k="<<k<<"\n";
        ExpressionValue<typename Base::Scalar> lhs_val = this->lhs_expr_.eval(k);
        ExpressionValue<typename Base::Scalar> rhs_val = this->rhs_expr_.eval(k);
        
        typename Base::Scalar scalar_value = rhs_val(0, 0)(0, 0);
        ExpressionValue<typename Base::Scalar> result(lhs_val.rowCardinality(), lhs_val.colCardinality());
        
        for (index_t i = 0; i < result.rowCardinality(); ++i)
        {
            for (index_t j = 0; j < result.colCardinality(); ++j)
            {
                result(i, j) = lhs_val(i, j).array() * scalar_value;
            }
        }
        
        return result;
    }

    // --- eval() specialization 5: Matrix-Matrix product (2,2) None x None ---
    template <size_t LO = _LhsOrder, size_t RO = _RhsOrder, size_t LS = _LhsSpace, size_t RS = _RhsSpace>
    typename std::enable_if<LO == 2 && RO == 2 && LS == SpaceType::None && RS == SpaceType::None,
                            ExpressionValue<typename Base::Scalar>>::type
    eval(const index_t k) const
    {
        gsDebug<<"ProductExpression: eval Matrix * Matrix (None, None) at k="<<k<<"\n";
        ExpressionValue<typename Base::Scalar> lhs_val = this->lhs_expr_.eval(k);
        ExpressionValue<typename Base::Scalar> rhs_val = this->rhs_expr_.eval(k);
        
        ExpressionValue<typename Base::Scalar> result(1, 1);
        result(0, 0) = lhs_val(0, 0) * rhs_val(0, 0);
        return result;
    }

    // --- eval() specialization 6: Matrix-Matrix product (2,2) Same Space x Same Space ---
    template <size_t LO = _LhsOrder, size_t RO = _RhsOrder, size_t LS = _LhsSpace, size_t RS = _RhsSpace>
    typename std::enable_if<LO == 2 && RO == 2 && LS != SpaceType::None && RS != SpaceType::None && LS == RS,
                            ExpressionValue<typename Base::Scalar>>::type
    eval(const index_t k) const
    {
        gsDebug<<"ProductExpression: eval Matrix * Matrix (same space) at k="<<k<<"\n";
        ExpressionValue<typename Base::Scalar> lhs_val = this->lhs_expr_.eval(k);
        ExpressionValue<typename Base::Scalar> rhs_val = this->rhs_expr_.eval(k);
        
        GISMO_ENSURE(lhs_val.rowCardinality() == rhs_val.rowCardinality() &&
                    lhs_val.colCardinality() == rhs_val.colCardinality(),
                    "ProductExpression: Cardinality mismatch in matrix-matrix product");
        
        ExpressionValue<typename Base::Scalar> result(lhs_val.rowCardinality(), lhs_val.colCardinality());
        
        for (index_t i = 0; i < result.rowCardinality(); ++i)
        {
            for (index_t j = 0; j < result.colCardinality(); ++j)
            {
                result(i, j) = lhs_val(i, j) * rhs_val(i, j);
            }
        }
        
        return result;
    }

    // --- eval() specialization 7: Matrix-Vector product (2,1) None x None ---
    template <size_t LO = _LhsOrder, size_t RO = _RhsOrder, size_t LS = _LhsSpace, size_t RS = _RhsSpace>
    typename std::enable_if<LO == 2 && RO == 1 && LS == SpaceType::None && RS == SpaceType::None,
                            ExpressionValue<typename Base::Scalar>>::type
    eval(const index_t k) const
    {
        gsDebug<<"ProductExpression: eval Matrix * Vector (None, None) at k="<<k<<"\n";
        ExpressionValue<typename Base::Scalar> lhs_val = this->lhs_expr_.eval(k);
        ExpressionValue<typename Base::Scalar> rhs_val = this->rhs_expr_.eval(k);
        
        ExpressionValue<typename Base::Scalar> result(1, 1);
        result(0, 0) = lhs_val(0, 0) * rhs_val(0, 0);
        return result;
    }

    // --- eval() specialization 8: Matrix-Vector product (2,1) Same Space x Same Space ---
    template <size_t LO = _LhsOrder, size_t RO = _RhsOrder, size_t LS = _LhsSpace, size_t RS = _RhsSpace>
    typename std::enable_if<LO == 2 && RO == 1 && LS != SpaceType::None && RS != SpaceType::None && LS == RS,
                            ExpressionValue<typename Base::Scalar>>::type
    eval(const index_t k) const
    {
        gsDebug<<"ProductExpression: eval Matrix * Vector (same space) at k="<<k<<"\n";
        ExpressionValue<typename Base::Scalar> lhs_val = this->lhs_expr_.eval(k);
        ExpressionValue<typename Base::Scalar> rhs_val = this->rhs_expr_.eval(k);
        
        GISMO_ENSURE(lhs_val.rowCardinality() == rhs_val.rowCardinality() &&
                    lhs_val.colCardinality() == rhs_val.colCardinality(),
                    "ProductExpression: Cardinality mismatch in matrix-vector product");
        
        ExpressionValue<typename Base::Scalar> result(lhs_val.rowCardinality(), lhs_val.colCardinality());
        
        for (index_t i = 0; i < result.rowCardinality(); ++i)
        {
            for (index_t j = 0; j < result.colCardinality(); ++j)
            {
                result(i, j) = lhs_val(i, j) * rhs_val(i, j);
            }
        }
        
        return result;
    }

    void print(std::ostream & os) const
    {
        os << this->lhs_expr_ << "*" << this->rhs_expr_;
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
ProductExpression<ConstantObject<typename ExpressionTraits<_LhsExpr>::Scalar,0>, _LhsExpr, 0, _LhsExpr::Order, 0, _LhsExpr::Space>
operator*(const BaseExpression<_LhsExpr>& lhs, const typename ExpressionTraits<_LhsExpr>::Scalar rhs)
{
    return ProductExpression<ConstantObject<typename ExpressionTraits<_LhsExpr>::Scalar,0>, _LhsExpr, 0, _LhsExpr::Order, 0, _LhsExpr::Space>(ConstantObject<typename ExpressionTraits<_LhsExpr>::Scalar,0>(rhs), lhs);
}

// Specialization: Scalar primitive × Expression
template <typename _RhsExpr>
ProductExpression<ConstantObject<typename ExpressionTraits<_RhsExpr>::Scalar,0>, _RhsExpr, 0, _RhsExpr::Order, 0, _RhsExpr::Space>
operator*(const typename ExpressionTraits<_RhsExpr>::Scalar lhs, const BaseExpression<_RhsExpr>& rhs)
{
    return ProductExpression<ConstantObject<typename ExpressionTraits<_RhsExpr>::Scalar,0>, _RhsExpr, 0, _RhsExpr::Order, 0, _RhsExpr::Space>(ConstantObject<typename ExpressionTraits<_RhsExpr>::Scalar,0>(lhs), rhs);
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
template <typename _T, size_t _LhsSpace, size_t _LhsOrder, typename _RhsExpr>
auto operator*(const NullObject<_T,_LhsSpace,_LhsOrder>& /* lhs */, const BaseExpression<_RhsExpr>& /* rhs */)
-> NullObject<_T,SpaceType::None,0>
{
    return NullObject<_T,SpaceType::None,0>::get();
}

template <typename _LhsExpr, typename _T, size_t _RhsSpace, size_t _RhsOrder>
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
