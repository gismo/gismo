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


// --- Generic ExpressionTraits for AddExpression ---
template <typename _LhsExpr, typename _RhsExpr, size_t _LhsOrder, size_t _RhsOrder, size_t _LhsSpace, size_t _RhsSpace>
struct ExpressionTraits<AddExpression<_LhsExpr, _RhsExpr, _LhsOrder, _RhsOrder, _LhsSpace, _RhsSpace>>
{
    typedef _LhsExpr LhsType;
    typedef _RhsExpr RhsType;

    typedef typename ExpressionTraits<_LhsExpr>::Scalar Scalar;
    static constexpr size_t Order = (_LhsOrder > _RhsOrder) ? _LhsOrder : _RhsOrder; // Max order
    static constexpr size_t Space = (_LhsSpace != SpaceType::None && _RhsSpace != SpaceType::None && _LhsSpace != _RhsSpace) 
                                   ? static_cast<size_t>(SpaceType::Both)
                                   : ((_LhsSpace != SpaceType::None) ? _LhsSpace : _RhsSpace); // Both if different non-None spaces
    static constexpr size_t Deriv = (ExpressionTraits<_LhsExpr>::Deriv > ExpressionTraits<_RhsExpr>::Deriv)
                                   ? ExpressionTraits<_LhsExpr>::Deriv : ExpressionTraits<_RhsExpr>::Deriv;
    static constexpr bool IsConstant = ExpressionTraits<_LhsExpr>::IsConstant && ExpressionTraits<_RhsExpr>::IsConstant;
};

// --- Unified AddExpression class with enable_if specializations for eval() ---
template <typename _LhsExpr, typename _RhsExpr, size_t _LhsOrder, size_t _RhsOrder, size_t _LhsSpace, size_t _RhsSpace>
class AddExpression : public BinaryOperator<AddExpression<_LhsExpr, _RhsExpr, _LhsOrder, _RhsOrder, _LhsSpace, _RhsSpace>>
{
    using Base = BinaryOperator<AddExpression<_LhsExpr, _RhsExpr, _LhsOrder, _RhsOrder, _LhsSpace, _RhsSpace>>;

private:
    // Helper to emit warning for Test + Trial case
    template <size_t LS = _LhsSpace, size_t RS = _RhsSpace>
    typename std::enable_if<(LS == SpaceType::Test && RS == SpaceType::Trial) || 
                            (LS == SpaceType::Trial && RS == SpaceType::Test)>::type
    check_bilinear_warning() const
    {
        gsWarn << "AddExpression: Test + Trial spaces broadcast to bilinear form. Consider using specialized bilinear operators.\n";
    }

    template <size_t LS = _LhsSpace, size_t RS = _RhsSpace>
    typename std::enable_if<!(( LS == SpaceType::Test && RS == SpaceType::Trial) || 
                             (LS == SpaceType::Trial && RS == SpaceType::Test))>::type
    check_bilinear_warning() const
    {
        // No warning for other cases
    }

public:
    AddExpression(const _LhsExpr& lhs, const _RhsExpr& rhs)
        : Base(lhs, rhs)
    {
        check_bilinear_warning();
        
        if (_LhsOrder == _RhsOrder && _LhsOrder > 0) {
            for (size_t d = 0; d != _LhsOrder; ++d) {
                GISMO_ENSURE(this->lhs_expr_.sizes()[d] == this->rhs_expr_.sizes()[d],
                    "AddExpression requires same sizes");
                this->sizes_[d] = this->lhs_expr_.sizes()[d];
            }
        } else if (_LhsOrder > _RhsOrder) {
            for (size_t d = 0; d != _LhsOrder; ++d) {
                this->sizes_[d] = this->lhs_expr_.sizes()[d];
            }
        } else if (_RhsOrder > _LhsOrder) {
            for (size_t d = 0; d != _RhsOrder; ++d) {
                this->sizes_[d] = this->rhs_expr_.sizes()[d];
            }
        }
    }

    size_t domainDim() const
    {
        return (_LhsOrder >= _RhsOrder) ? this->lhs_expr_.domainDim() : this->rhs_expr_.domainDim();
    }

    // --- eval() specialization 1: None + None ---
    template <size_t LS = _LhsSpace, size_t RS = _RhsSpace>
    typename std::enable_if<LS == SpaceType::None && RS == SpaceType::None,
                            ExpressionValue<typename Base::Scalar>>::type
    eval(const index_t k) const
    {
        gsDebug<<"AddExpression: eval None + None at k="<<k<<"\n";
        ExpressionValue<typename Base::Scalar> lhs_val = this->lhs_expr_.eval(k);
        ExpressionValue<typename Base::Scalar> rhs_val = this->rhs_expr_.eval(k);
        
        // Both are scalar (1,1) - no loop needed
        ExpressionValue<typename Base::Scalar> result(1, 1);
        const auto& lhs_mat = lhs_val(0, 0);
        const auto& rhs_mat = rhs_val(0, 0);
        
        if (lhs_mat.rows() == rhs_mat.rows() && lhs_mat.cols() == rhs_mat.cols())
        {
            result(0, 0) = lhs_mat + rhs_mat;
        }
        else if (lhs_mat.rows() == 1 && lhs_mat.cols() == 1)
        {
            result(0, 0) = rhs_mat.array() + lhs_mat(0, 0);
        }
        else if (rhs_mat.rows() == 1 && rhs_mat.cols() == 1)
        {
            result(0, 0) = lhs_mat.array() + rhs_mat(0, 0);
        }
        else
        {
            GISMO_ERROR("AddExpression: Incompatible matrix dimensions. LHS=(" 
                       << lhs_mat.rows() << "," << lhs_mat.cols() << "), RHS=(" 
                       << rhs_mat.rows() << "," << rhs_mat.cols() << ")");
        }
        return result;
    }

    // --- eval() specialization 2: None + Space ---
    template <size_t LS = _LhsSpace, size_t RS = _RhsSpace>
    typename std::enable_if<LS == SpaceType::None && RS != SpaceType::None,
                            ExpressionValue<typename Base::Scalar>>::type
    eval(const index_t k) const
    {
        gsDebug<<"AddExpression: eval None + Space at k="<<k<<"\n";
        ExpressionValue<typename Base::Scalar> lhs_val = this->lhs_expr_.eval(k);
        ExpressionValue<typename Base::Scalar> rhs_val = this->rhs_expr_.eval(k);
        
        // LHS is scalar (1,1), RHS has space dependencies
        ExpressionValue<typename Base::Scalar> result(rhs_val.rowCardinality(), rhs_val.colCardinality());
        const auto& lhs_mat = lhs_val(0, 0);
        
        for (index_t i = 0; i < result.rowCardinality(); ++i)
        {
            for (index_t j = 0; j < result.colCardinality(); ++j)
            {
                const auto& rhs_mat = rhs_val(i, j);
                
                if (lhs_mat.rows() == rhs_mat.rows() && lhs_mat.cols() == rhs_mat.cols())
                {
                    result(i, j) = lhs_mat + rhs_mat;
                }
                else if (lhs_mat.rows() == 1 && lhs_mat.cols() == 1)
                {
                    result(i, j) = rhs_mat.array() + lhs_mat(0, 0);
                }
                else if (rhs_mat.rows() == 1 && rhs_mat.cols() == 1)
                {
                    result(i, j) = lhs_mat.array() + rhs_mat(0, 0);
                }
                else
                {
                    GISMO_ERROR("AddExpression: Incompatible matrix dimensions at (" << i << "," << j << ")");
                }
            }
        }
        return result;
    }

    // --- eval() specialization 3: Space + None ---
    template <size_t LS = _LhsSpace, size_t RS = _RhsSpace>
    typename std::enable_if<LS != SpaceType::None && RS == SpaceType::None,
                            ExpressionValue<typename Base::Scalar>>::type
    eval(const index_t k) const
    {
        gsDebug<<"AddExpression: eval Space + None at k="<<k<<"\n";
        ExpressionValue<typename Base::Scalar> lhs_val = this->lhs_expr_.eval(k);
        ExpressionValue<typename Base::Scalar> rhs_val = this->rhs_expr_.eval(k);
        
        // LHS has space dependencies, RHS is scalar (1,1)
        ExpressionValue<typename Base::Scalar> result(lhs_val.rowCardinality(), lhs_val.colCardinality());
        const auto& rhs_mat = rhs_val(0, 0);
        
        for (index_t i = 0; i < result.rowCardinality(); ++i)
        {
            for (index_t j = 0; j < result.colCardinality(); ++j)
            {
                const auto& lhs_mat = lhs_val(i, j);
                
                if (lhs_mat.rows() == rhs_mat.rows() && lhs_mat.cols() == rhs_mat.cols())
                {
                    result(i, j) = lhs_mat + rhs_mat;
                }
                else if (lhs_mat.rows() == 1 && lhs_mat.cols() == 1)
                {
                    result(i, j) = rhs_mat.array() + lhs_mat(0, 0);
                }
                else if (rhs_mat.rows() == 1 && rhs_mat.cols() == 1)
                {
                    result(i, j) = lhs_mat.array() + rhs_mat(0, 0);
                }
                else
                {
                    GISMO_ERROR("AddExpression: Incompatible matrix dimensions at (" << i << "," << j << ")");
                }
            }
        }
        return result;
    }

    // --- eval() specialization 4: Same Space + Same Space ---
    template <size_t LS = _LhsSpace, size_t RS = _RhsSpace>
    typename std::enable_if<LS != SpaceType::None && RS != SpaceType::None && LS == RS,
                            ExpressionValue<typename Base::Scalar>>::type
    eval(const index_t k) const
    {
        gsDebug<<"AddExpression: eval Same Space + Same Space at k="<<k<<"\n";
        ExpressionValue<typename Base::Scalar> lhs_val = this->lhs_expr_.eval(k);
        ExpressionValue<typename Base::Scalar> rhs_val = this->rhs_expr_.eval(k);
        
        // Both have same non-None space: require exact cardinality match
        GISMO_ENSURE(lhs_val.rowCardinality() == rhs_val.rowCardinality() &&
                    lhs_val.colCardinality() == rhs_val.colCardinality(),
                    "AddExpression: Cardinality mismatch. LHS=(" << lhs_val.rowCardinality() << "," << lhs_val.colCardinality() <<
                    "), RHS=(" << rhs_val.rowCardinality() << "," << rhs_val.colCardinality() << ")");
        
        ExpressionValue<typename Base::Scalar> result(lhs_val.rowCardinality(), lhs_val.colCardinality());
        
        for (index_t i = 0; i < result.rowCardinality(); ++i)
        {
            for (index_t j = 0; j < result.colCardinality(); ++j)
            {
                const auto& lhs_mat = lhs_val(i, j);
                const auto& rhs_mat = rhs_val(i, j);
                
                if (lhs_mat.rows() == rhs_mat.rows() && lhs_mat.cols() == rhs_mat.cols())
                {
                    result(i, j) = lhs_mat + rhs_mat;
                }
                else if (lhs_mat.rows() == 1 && lhs_mat.cols() == 1)
                {
                    result(i, j) = rhs_mat.array() + lhs_mat(0, 0);
                }
                else if (rhs_mat.rows() == 1 && rhs_mat.cols() == 1)
                {
                    result(i, j) = lhs_mat.array() + rhs_mat(0, 0);
                }
                else
                {
                    GISMO_ERROR("AddExpression: Incompatible matrix dimensions at (" << i << "," << j << ")");
                }
            }
        }
        
        return result;
    }

    // --- eval() specialization 5: Test + Trial ---
    template <size_t LS = _LhsSpace, size_t RS = _RhsSpace>
    typename std::enable_if<LS == SpaceType::Test && RS == SpaceType::Trial,
                            ExpressionValue<typename Base::Scalar>>::type
    eval(const index_t k) const
    {
        // TODO: Check more efficient implementation (assuming constant size of matrix per test/trial combination).
        // TODO: Consider re-using this for specialization 6 with swapped roles.
        gsDebug<<"AddExpression: eval Test + Trial at k="<<k<<"\n";
        ExpressionValue<typename Base::Scalar> lhs_val = this->lhs_expr_.eval(k);
        ExpressionValue<typename Base::Scalar> rhs_val = this->rhs_expr_.eval(k);
        
        // Bilinear broadcast: (N,1) + (1,M) -> (N,M)
        index_t result_rows = std::max(lhs_val.rowCardinality(), rhs_val.rowCardinality());
        index_t result_cols = std::max(lhs_val.colCardinality(), rhs_val.colCardinality());
        ExpressionValue<typename Base::Scalar> result(result_rows, result_cols);
        
        for (index_t i = 0; i < result_rows; ++i)
        {
            for (index_t j = 0; j < result_cols; ++j)
            {
                index_t lhs_i = (lhs_val.rowCardinality() > 1) ? i : 0;
                index_t lhs_j = (lhs_val.colCardinality() > 1) ? j : 0;
                index_t rhs_i = (rhs_val.rowCardinality() > 1) ? i : 0;
                index_t rhs_j = (rhs_val.colCardinality() > 1) ? j : 0;
                
                const auto& lhs_mat = lhs_val(lhs_i, lhs_j);
                const auto& rhs_mat = rhs_val(rhs_i, rhs_j);
                
                if (lhs_mat.rows() == rhs_mat.rows() && lhs_mat.cols() == rhs_mat.cols())
                {
                    result(i, j) = lhs_mat + rhs_mat;
                }
                else if (lhs_mat.rows() == 1 && lhs_mat.cols() == 1)
                {
                    result(i, j) = rhs_mat.array() + lhs_mat(0, 0);
                }
                else if (rhs_mat.rows() == 1 && rhs_mat.cols() == 1)
                {
                    result(i, j) = lhs_mat.array() + rhs_mat(0, 0);
                }
                else
                {
                    GISMO_ERROR("AddExpression: Incompatible matrix dimensions at (" << i << "," << j << "). LHS=(" 
                               << lhs_mat.rows() << "," << lhs_mat.cols() << "), RHS=(" 
                               << rhs_mat.rows() << "," << rhs_mat.cols() << ")");
                }
            }
        }
        
        return result;
    }

    // --- eval() specialization 6: Trial + Test ---
    template <size_t LS = _LhsSpace, size_t RS = _RhsSpace>
    typename std::enable_if<LS == SpaceType::Trial && RS == SpaceType::Test,
                            ExpressionValue<typename Base::Scalar>>::type
    eval(const index_t k) const
    {
        gsDebug<<"AddExpression: eval Trial + Test at k="<<k<<"\n";
        ExpressionValue<typename Base::Scalar> lhs_val = this->lhs_expr_.eval(k);
        ExpressionValue<typename Base::Scalar> rhs_val = this->rhs_expr_.eval(k);
        
        // Bilinear broadcast: (1,M) + (N,1) -> (N,M)
        index_t result_rows = std::max(lhs_val.rowCardinality(), rhs_val.rowCardinality());
        index_t result_cols = std::max(lhs_val.colCardinality(), rhs_val.colCardinality());
        ExpressionValue<typename Base::Scalar> result(result_rows, result_cols);
        
        for (index_t i = 0; i < result_rows; ++i)
        {
            for (index_t j = 0; j < result_cols; ++j)
            {
                index_t lhs_i = (lhs_val.rowCardinality() > 1) ? i : 0;
                index_t lhs_j = (lhs_val.colCardinality() > 1) ? j : 0;
                index_t rhs_i = (rhs_val.rowCardinality() > 1) ? i : 0;
                index_t rhs_j = (rhs_val.colCardinality() > 1) ? j : 0;
                
                const auto& lhs_mat = lhs_val(lhs_i, lhs_j);
                const auto& rhs_mat = rhs_val(rhs_i, rhs_j);
                
                if (lhs_mat.rows() == rhs_mat.rows() && lhs_mat.cols() == rhs_mat.cols())
                {
                    result(i, j) = lhs_mat + rhs_mat;
                }
                else if (lhs_mat.rows() == 1 && lhs_mat.cols() == 1)
                {
                    result(i, j) = rhs_mat.array() + lhs_mat(0, 0);
                }
                else if (rhs_mat.rows() == 1 && rhs_mat.cols() == 1)
                {
                    result(i, j) = lhs_mat.array() + rhs_mat(0, 0);
                }
                else
                {
                    GISMO_ERROR("AddExpression: Incompatible matrix dimensions at (" << i << "," << j << "). LHS=(" 
                               << lhs_mat.rows() << "," << lhs_mat.cols() << "), RHS=(" 
                               << rhs_mat.rows() << "," << rhs_mat.cols() << ")");
                }
            }
        }
        
        return result;
    }

    void print(std::ostream & os) const
    {
        os << this->lhs_expr_ << "+" << this->rhs_expr_;
    }

    std::string label() const { return this->lhs_expr_.label() + " + " + this->rhs_expr_.label(); }
    
};

// Generic operator+ to create AddExpression instances using SFINAE
template <typename _LhsExpr, typename _RhsExpr>
typename std::enable_if<
    _LhsExpr::Order == _RhsExpr::Order,
    AddExpression<_LhsExpr, _RhsExpr, _LhsExpr::Order, _RhsExpr::Order, _LhsExpr::Space, _RhsExpr::Space>
>::type
operator+(const BaseExpression<_LhsExpr>& lhs, const BaseExpression<_RhsExpr>& rhs)
{
    return AddExpression<_LhsExpr, _RhsExpr, _LhsExpr::Order, _RhsExpr::Order, _LhsExpr::Space, _RhsExpr::Space>(lhs, rhs);
}

// TODO: Add of NullObject

// Specialization for transpose expressions (these should trigger compile errors)
template <typename _LhsExpr, typename _RhsExpr>
typename std::enable_if<
    _LhsExpr::Order == _RhsExpr::Order &&
    _LhsExpr::Space == _RhsExpr::Space,
    AddExpression<TransposeExpression<_LhsExpr>, _RhsExpr, ExpressionTraits<TransposeExpression<_LhsExpr>>::Order, _RhsExpr::Order, ExpressionTraits<TransposeExpression<_LhsExpr>>::Space, _RhsExpr::Space>
>::type
operator+(const TransposeExpression<_LhsExpr>& lhs, const BaseExpression<_RhsExpr>& rhs)
{
    GISMO_ERROR("Addition of a transposed vector and a vector is not defined.");
    return AddExpression<TransposeExpression<_LhsExpr>, _RhsExpr, ExpressionTraits<TransposeExpression<_LhsExpr>>::Order, _RhsExpr::Order, ExpressionTraits<TransposeExpression<_LhsExpr>>::Space, _RhsExpr::Space>(lhs, rhs);
}

template <typename _LhsExpr, typename _RhsExpr>
typename std::enable_if<
    _LhsExpr::Order == _RhsExpr::Order &&
    _LhsExpr::Space == _RhsExpr::Space,
    AddExpression<_LhsExpr, TransposeExpression<_RhsExpr>, _LhsExpr::Order, ExpressionTraits<TransposeExpression<_RhsExpr>>::Order, _LhsExpr::Space, ExpressionTraits<TransposeExpression<_RhsExpr>>::Space>
>::type
operator+(const BaseExpression<_LhsExpr>& lhs, const TransposeExpression<_RhsExpr>& rhs)
{
    GISMO_ERROR("Addition of a vector and a transposed vector is not defined.");
    return AddExpression<_LhsExpr, TransposeExpression<_RhsExpr>, _LhsExpr::Order, ExpressionTraits<TransposeExpression<_RhsExpr>>::Order, _LhsExpr::Space, ExpressionTraits<TransposeExpression<_RhsExpr>>::Space>(lhs, rhs);
}

// Specialization for scalars
template <typename _LhsExpr>
typename std::enable_if<
    !std::is_same<_LhsExpr, ConstantObject<typename ExpressionTraits<_LhsExpr>::Scalar,0>>::value,
    AddExpression<_LhsExpr, ConstantObject<typename ExpressionTraits<_LhsExpr>::Scalar,0>, _LhsExpr::Order, 0, _LhsExpr::Space, SpaceType::None>
>::type
operator+(const BaseExpression<_LhsExpr>& lhs, const typename ExpressionTraits<_LhsExpr>::Scalar rhs)
{
    return AddExpression<_LhsExpr, ConstantObject<typename ExpressionTraits<_LhsExpr>::Scalar,0>, _LhsExpr::Order, 0, _LhsExpr::Space, SpaceType::None>(lhs, ConstantObject<typename ExpressionTraits<_LhsExpr>::Scalar,0>(rhs));
}

template <typename _RhsExpr>
typename std::enable_if<
    !std::is_same<_RhsExpr, ConstantObject<typename ExpressionTraits<_RhsExpr>::Scalar,0>>::value,
    AddExpression<ConstantObject<typename ExpressionTraits<_RhsExpr>::Scalar,0>, _RhsExpr, 0, _RhsExpr::Order, SpaceType::None, _RhsExpr::Space>
>::type
operator+(const typename ExpressionTraits<_RhsExpr>::Scalar lhs, const BaseExpression<_RhsExpr>& rhs)
{
    return AddExpression<ConstantObject<typename ExpressionTraits<_RhsExpr>::Scalar,0>, _RhsExpr, 0, _RhsExpr::Order, SpaceType::None, _RhsExpr::Space>(ConstantObject<typename ExpressionTraits<_RhsExpr>::Scalar,0>(lhs), rhs);
}

// Specialization for vectors
template <typename _LhsExpr, int _Rows>
typename std::enable_if<
    !std::is_same<_LhsExpr, ConstantObject<typename ExpressionTraits<_LhsExpr>::Scalar,_LhsExpr::Order>>::value,
    AddExpression<_LhsExpr, ConstantObject<typename ExpressionTraits<_LhsExpr>::Scalar,_LhsExpr::Order>, _LhsExpr::Order, _LhsExpr::Order, _LhsExpr::Space, _LhsExpr::Space>
>::type
operator+(const BaseExpression<_LhsExpr>& lhs, const gsVector<typename ExpressionTraits<_LhsExpr>::Scalar,_Rows>& rhs)
{
    return AddExpression<_LhsExpr, ConstantObject<typename ExpressionTraits<_LhsExpr>::Scalar,_LhsExpr::Order>, _LhsExpr::Order, _LhsExpr::Order, _LhsExpr::Space, _LhsExpr::Space>(lhs, ConstantObject<typename ExpressionTraits<_LhsExpr>::Scalar,1>(rhs));
}

template <typename _RhsExpr, int _Rows>
typename std::enable_if<
    !std::is_same<_RhsExpr, ConstantObject<typename ExpressionTraits<_RhsExpr>::Scalar,_RhsExpr::Order>>::value,
    AddExpression<ConstantObject<typename ExpressionTraits<_RhsExpr>::Scalar,_RhsExpr::Order>, _RhsExpr, _RhsExpr::Order, _RhsExpr::Order, _RhsExpr::Space, _RhsExpr::Space>
>::type
operator+(const gsVector<typename ExpressionTraits<_RhsExpr>::Scalar,_Rows>& lhs, const BaseExpression<_RhsExpr>& rhs)
{
    return AddExpression<ConstantObject<typename ExpressionTraits<_RhsExpr>::Scalar,_RhsExpr::Order>, _RhsExpr, _RhsExpr::Order, _RhsExpr::Order, _RhsExpr::Space, _RhsExpr::Space>(ConstantObject<typename ExpressionTraits<_RhsExpr>::Scalar,1>(lhs), rhs);
}

// Specialization for matrices
template <typename _LhsExpr, int _Rows, int _Cols>
typename std::enable_if<
    !std::is_same<_LhsExpr, ConstantObject<typename ExpressionTraits<_LhsExpr>::Scalar,_LhsExpr::Order>>::value,
    AddExpression<_LhsExpr, ConstantObject<typename ExpressionTraits<_LhsExpr>::Scalar,_LhsExpr::Order>, _LhsExpr::Order, _LhsExpr::Order, _LhsExpr::Space, _LhsExpr::Space>
>::type
operator+(const BaseExpression<_LhsExpr>& lhs, const gsMatrix<typename ExpressionTraits<_LhsExpr>::Scalar,_Rows,_Cols>& rhs)
{
    return AddExpression<_LhsExpr, ConstantObject<typename ExpressionTraits<_LhsExpr>::Scalar,_LhsExpr::Order>, _LhsExpr::Order, _LhsExpr::Order, _LhsExpr::Space, _LhsExpr::Space>(lhs, ConstantObject<typename ExpressionTraits<_LhsExpr>::Scalar,2>(rhs));
}

template <typename _RhsExpr, int _Rows, int _Cols>
typename std::enable_if<
    !std::is_same<_RhsExpr, ConstantObject<typename ExpressionTraits<_RhsExpr>::Scalar,_RhsExpr::Order>>::value,
    AddExpression<ConstantObject<typename ExpressionTraits<_RhsExpr>::Scalar,_RhsExpr::Order>, _RhsExpr, _RhsExpr::Order, _RhsExpr::Order, _RhsExpr::Space, _RhsExpr::Space>
>::type
operator+(const gsMatrix<typename ExpressionTraits<_RhsExpr>::Scalar,_Rows,_Cols>& lhs, const BaseExpression<_RhsExpr>& rhs)
{
    return AddExpression<ConstantObject<typename ExpressionTraits<_RhsExpr>::Scalar,_RhsExpr::Order>, _RhsExpr, _RhsExpr::Order, _RhsExpr::Order, _RhsExpr::Space, _RhsExpr::Space>(ConstantObject<typename ExpressionTraits<_RhsExpr>::Scalar,2>(lhs), rhs);
}

template <typename _LhsExpr, typename _RhsExpr>
auto operator+(const NullObject<typename _LhsExpr::Scalar, _LhsExpr::Space, _LhsExpr::Order>& /* lhs */, const BaseExpression<_RhsExpr>& rhs)
-> BaseExpression<_RhsExpr>
{
    return rhs;
}

template <typename _LhsExpr, typename _RhsExpr>
auto operator+(const BaseExpression<_LhsExpr>& lhs, const NullObject<typename _RhsExpr::Scalar, _RhsExpr::Space, _RhsExpr::Order>& /* rhs */)
-> BaseExpression<_LhsExpr>
{
    return lhs;
}

// TODO: avoid ambiguity with (0,0) case
// template <typename _LhsExpr, typename _RhsExpr, typename _SpaceObject>
// SpaceObject<typename LhsExpr::Scalar,SpaceType::None,0>
// variation(const AddExpression<LhsExpr,RhsExpr,LhsExpr::Order,RhsExpr::Order,LhsExpr::Space,RhsExpr::Space> & /* expr */,
//           const _SpaceObject & /* space */)
// {
//     gsInfo<<"ZERO VARIATION ADD\n";
//     return NullObject<typename LhsExpr::Scalar,SpaceType::None,0>::get();
// }

// Solution + rhs expression

// Variation of AddExpression where LHS is a SolutionObject and RhsExpr is NOT a SolutionObject
template <typename _LhsExpr, typename _RhsExpr, typename _SpaceObject>
typename std::enable_if<
    !std::is_base_of<SolutionObject<typename _LhsExpr::Scalar,_LhsExpr::Space,_LhsExpr::Order>, _RhsExpr>::value,
    decltype(variation(std::declval<_RhsExpr>(), std::declval<_SpaceObject>()))
>::type
variation(const AddExpression<SolutionObject<typename _LhsExpr::Scalar,_LhsExpr::Space,_LhsExpr::Order>,_RhsExpr,_LhsExpr::Order,_RhsExpr::Order,SolutionObject<typename _LhsExpr::Scalar,_LhsExpr::Space,_LhsExpr::Order>::Order,_RhsExpr::Space> & expr,
          const _SpaceObject & space)
{
    return variation(expr.rhs(), space);
}

// Variation of AddExpression where RHS is a SolutionObject and LhsExpr is NOT a SolutionObject
template <typename _LhsExpr, typename _RhsExpr, typename _SpaceObject>
typename std::enable_if<
    !std::is_base_of<SolutionObject<typename _RhsExpr::Scalar,_RhsExpr::Space,_RhsExpr::Order>, _LhsExpr>::value,
    decltype(variation(std::declval<_LhsExpr>(), std::declval<_SpaceObject>()))
>::type
variation(const AddExpression<_LhsExpr,SolutionObject<typename _RhsExpr::Scalar,_RhsExpr::Space,_RhsExpr::Order>,_LhsExpr::Order,_RhsExpr::Order,_LhsExpr::Space,SolutionObject<typename _RhsExpr::Scalar,_RhsExpr::Space,_RhsExpr::Order>::Order> & expr,
          const _SpaceObject & space)
{
    return variation(expr.lhs(), space);
}

template<typename _LhsExpr, typename _RhsExpr, typename _SpaceObject>
auto variation(const AddExpression<_LhsExpr,_RhsExpr,_LhsExpr::Order,_RhsExpr::Order,_LhsExpr::Space,_RhsExpr::Space> & expr,
               const _SpaceObject & space)
-> decltype(variation(expr.lhs(), space) + variation(expr.rhs(), space))
{
    return variation(expr.lhs(), space) + variation(expr.rhs(), space);
}

}//namespace Expr
}//namespace gismo