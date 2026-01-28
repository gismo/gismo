/** @file InnerProductExpression.h

    @brief Inner product expression class

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
// --- InnerProductExpression using Partial Specialization (Redesigned) ---

// --- Generic ExpressionTraits for InnerProductExpression ---
// Inner product always results in a scalar (Order=0) regardless of operand orders/spaces
template <typename _LhsExpr, typename _RhsExpr, size_t _LhsOrder, size_t _RhsOrder, enum SpaceType _LhsSpace, enum SpaceType _RhsSpace>
struct ExpressionTraits<InnerProductExpression<_LhsExpr, _RhsExpr, _LhsOrder, _RhsOrder, _LhsSpace, _RhsSpace>>
{
    typedef _LhsExpr LhsType;
    typedef _RhsExpr RhsType;
    typedef typename ExpressionTraits<_LhsExpr>::Scalar Scalar;
    static constexpr size_t Order = 0;  // Inner product always produces a scalar
    
    // Space logic: if both different non-None, use Both; otherwise prefer non-None space
    static constexpr SpaceType Space = (_LhsSpace != SpaceType::None && _RhsSpace != SpaceType::None && _LhsSpace != _RhsSpace) 
                                   ? SpaceType::Both
                                   : ((_LhsSpace != SpaceType::None) ? _LhsSpace : _RhsSpace);
    
    static constexpr size_t Deriv = 0;
    static constexpr bool IsConstant = ExpressionTraits<_LhsExpr>::IsConstant && ExpressionTraits<_RhsExpr>::IsConstant;
};

// --- Unified InnerProductExpression (All Order/Space combinations via enable_if) ---
template <typename _LhsExpr, typename _RhsExpr, size_t _LhsOrder, size_t _RhsOrder, enum SpaceType _LhsSpace, enum SpaceType _RhsSpace>
class InnerProductExpression
 : public BinaryOperator<InnerProductExpression<_LhsExpr, _RhsExpr, _LhsOrder, _RhsOrder, _LhsSpace, _RhsSpace>>
{
    using Base = BinaryOperator<InnerProductExpression<_LhsExpr, _RhsExpr, _LhsOrder, _RhsOrder, _LhsSpace, _RhsSpace>>;
    
    static_assert(_LhsSpace == _RhsSpace || 
                  (_LhsSpace == SpaceType::Test && _RhsSpace == SpaceType::Trial) ||
                  (_LhsSpace == SpaceType::Trial && _RhsSpace == SpaceType::Test),
                  "InnerProductExpression requires same space or Test×Trial combination.");

public:
    InnerProductExpression(const _LhsExpr& lhs, const _RhsExpr& rhs) : Base(lhs, rhs)
    {
        // For same-space combinations, check they're the same space ID
        if (_LhsSpace == _RhsSpace)
        {
            if (_LhsSpace == SpaceType::Trial)
                GISMO_ASSERT(lhs.trial()==rhs.trial(),"InnerProductExpression requires both operands to have the same trial space. Probably lhs and rhs are different spaces (e.g., with different ID).");
            else if (_LhsSpace == SpaceType::Test)
                GISMO_ASSERT(lhs.test()==rhs.test(),"InnerProductExpression requires both operands to have the same test space. Probably lhs and rhs are different spaces (e.g., with different ID).");
        }
    }
    size_t domainDim() const { return this->lhs_expr_.domainDim(); }

    // Inner product always produces a scalar (Order=0), but test/trial return the underlying space objects
    // The returned SpaceObject should have the order of the base space (e.g., 0 for scalar spaces)
    auto test() const -> decltype(std::declval<_LhsExpr>().test())
    {
        if (_LhsSpace == SpaceType::Test)
            return this->lhs_expr_.test();
        else
            return this->rhs_expr_.test();
    }

    auto trial() const -> decltype(std::declval<_RhsExpr>().trial())
    {
        if (_LhsSpace == SpaceType::Trial)
            return this->lhs_expr_.trial();
        else
            return this->rhs_expr_.trial();
    }

    // --- eval() specialization 1: Vector dot (Order=1), Space=None ---
    template <size_t LO = _LhsOrder, size_t RO = _RhsOrder, size_t LS = _LhsSpace, size_t RS = _RhsSpace>
    typename std::enable_if<LO == 1 && RO == 1 && LS == SpaceType::None && RS == SpaceType::None,
                            ExpressionResult<typename Base::Scalar>>::type
    eval(const index_t k) const
    {
        ExpressionResult<typename Base::Scalar> lhs_eval = this->lhs_expr_.eval(k);
        ExpressionResult<typename Base::Scalar> rhs_eval = this->rhs_expr_.eval(k);

        GISMO_ASSERT(lhs_eval.rowCardinality() == 1 && lhs_eval.colCardinality() == 1,
                    "InnerProductExpression (1,1,None,None): Expected scalar cardinality");
        GISMO_ASSERT(rhs_eval.rowCardinality() == 1 && rhs_eval.colCardinality() == 1,
                    "InnerProductExpression (1,1,None,None): Expected scalar cardinality");
        
        ExpressionResult<typename Base::Scalar> result_val(1, 1);
        const auto& lhs_mat = lhs_eval(0, 0);
        const auto& rhs_mat = rhs_eval(0, 0);
        typename Base::Scalar dot_result = 0;
        for (index_t r = 0; r < lhs_mat.rows(); ++r)
            dot_result += lhs_mat(r, 0) * rhs_mat(r, 0);
        gsMatrix<typename Base::Scalar> res(1, 1);
        res(0, 0) = dot_result;
        result_val(0, 0) = res;
        return result_val;
    }

    // --- eval() specialization 2: Vector dot (Order=1), Space=Test ---
    template <size_t LO = _LhsOrder, size_t RO = _RhsOrder, size_t LS = _LhsSpace, size_t RS = _RhsSpace>
    typename std::enable_if<(LO == 1) && (RO == 1) && (LS == SpaceType::Test) && (RS == SpaceType::Test),
                            ExpressionResult<typename Base::Scalar>>::type
    eval(const index_t k) const
    {
        ExpressionResult<typename Base::Scalar> lhs_eval = this->lhs_expr_.eval(k);
        ExpressionResult<typename Base::Scalar> rhs_eval = this->rhs_expr_.eval(k);

        GISMO_ENSURE(lhs_eval.rowCardinality() == rhs_eval.rowCardinality() &&
                    lhs_eval.colCardinality() == rhs_eval.colCardinality(),
                    "InnerProductExpression (1,1,Test,Test): Cardinality mismatch");
        
        ExpressionResult<typename Base::Scalar> result_val(lhs_eval.rowCardinality(), lhs_eval.colCardinality());

        for (index_t i = 0; i < result_val.rowCardinality(); ++i)
        {
            for (index_t j = 0; j < result_val.colCardinality(); ++j)
            {
                const auto& lhs_mat = lhs_eval(i, j);
                const auto& rhs_mat = rhs_eval(i, j);
                typename Base::Scalar dot_result = 0;
                for (index_t r = 0; r < lhs_mat.rows(); ++r)
                    dot_result += lhs_mat(r, 0) * rhs_mat(r, 0);
                gsMatrix<typename Base::Scalar> res(1, 1);
                res(0, 0) = dot_result;
                result_val(i, j) = res;
            }
        }
        return result_val;
    }

    // --- eval() specialization 3: Vector dot (Order=1), Space=Trial ---
    template <size_t LO = _LhsOrder, size_t RO = _RhsOrder, size_t LS = _LhsSpace, size_t RS = _RhsSpace>
    typename std::enable_if<LO == 1 && RO == 1 && LS == SpaceType::Trial && RS == SpaceType::Trial,
                            ExpressionResult<typename Base::Scalar>>::type
    eval(const index_t k) const
    {
        ExpressionResult<typename Base::Scalar> lhs_eval = this->lhs_expr_.eval(k);
        ExpressionResult<typename Base::Scalar> rhs_eval = this->rhs_expr_.eval(k);

        GISMO_ENSURE(lhs_eval.rowCardinality() == rhs_eval.rowCardinality() &&
                    lhs_eval.colCardinality() == rhs_eval.colCardinality(),
                    "InnerProductExpression (1,1,Trial,Trial): Cardinality mismatch");
        
        ExpressionResult<typename Base::Scalar> result_val(lhs_eval.rowCardinality(), lhs_eval.colCardinality());

        for (index_t i = 0; i < result_val.rowCardinality(); ++i)
        {
            for (index_t j = 0; j < result_val.colCardinality(); ++j)
            {
                const auto& lhs_mat = lhs_eval(i, j);
                const auto& rhs_mat = rhs_eval(i, j);
                typename Base::Scalar dot_result = 0;
                for (index_t r = 0; r < lhs_mat.rows(); ++r)
                    dot_result += lhs_mat(r, 0) * rhs_mat(r, 0);
                gsMatrix<typename Base::Scalar> res(1, 1);
                res(0, 0) = dot_result;
                result_val(i, j) = res;
            }
        }
        return result_val;
    }

    // --- eval() specialization 4: Vector dot (Order=1), Space=Test×Trial ---
    template <size_t LO = _LhsOrder, size_t RO = _RhsOrder, size_t LS = _LhsSpace, size_t RS = _RhsSpace>
    typename std::enable_if<LO == 1 && RO == 1 && LS == SpaceType::Test && RS == SpaceType::Trial,
                            ExpressionResult<typename Base::Scalar>>::type
    eval(const index_t k) const
    {
        ExpressionResult<typename Base::Scalar> lhs_eval = this->lhs_expr_.eval(k);
        ExpressionResult<typename Base::Scalar> rhs_eval = this->rhs_expr_.eval(k);

        // Test has row cardinality, Trial has col cardinality
        ExpressionResult<typename Base::Scalar> result_val(lhs_eval.rowCardinality(), rhs_eval.colCardinality());

        for (index_t i = 0; i < result_val.rowCardinality(); ++i)
        {
            for (index_t j = 0; j < result_val.colCardinality(); ++j)
            {
                const auto& lhs_mat = lhs_eval(i, 0);
                const auto& rhs_mat = rhs_eval(0, j);
                typename Base::Scalar dot_result = 0;
                for (index_t r = 0; r < lhs_mat.rows(); ++r)
                    dot_result += lhs_mat(r, 0) * rhs_mat(r, 0);
                gsMatrix<typename Base::Scalar> res(1, 1);
                res(0, 0) = dot_result;
                result_val(i, j) = res;
            }
        }
        return result_val;
    }

    // --- eval() specialization 5: Vector dot (Order=1), Space=Trial×Test ---
    template <size_t LO = _LhsOrder, size_t RO = _RhsOrder, size_t LS = _LhsSpace, size_t RS = _RhsSpace>
    typename std::enable_if<LO == 1 && RO == 1 && LS == SpaceType::Trial && RS == SpaceType::Test,
                            ExpressionResult<typename Base::Scalar>>::type
    eval(const index_t k) const
    {
        ExpressionResult<typename Base::Scalar> lhs_eval = this->lhs_expr_.eval(k);
        ExpressionResult<typename Base::Scalar> rhs_eval = this->rhs_expr_.eval(k);

        // Trial has col cardinality, Test has row cardinality
        ExpressionResult<typename Base::Scalar> result_val(rhs_eval.rowCardinality(), lhs_eval.colCardinality());

        for (index_t i = 0; i < result_val.rowCardinality(); ++i)
        {
            for (index_t j = 0; j < result_val.colCardinality(); ++j)
            {
                const auto& lhs_mat = lhs_eval(0, j);
                const auto& rhs_mat = rhs_eval(i, 0);
                typename Base::Scalar dot_result = 0;
                for (index_t r = 0; r < lhs_mat.rows(); ++r)
                    dot_result += lhs_mat(r, 0) * rhs_mat(r, 0);
                gsMatrix<typename Base::Scalar> res(1, 1);
                res(0, 0) = dot_result;
                result_val(i, j) = res;
            }
        }
        return result_val;
    }

    // --- eval() specialization 6: Matrix ddot (Order=2), Space=None ---
    template <size_t LO = _LhsOrder, size_t RO = _RhsOrder, size_t LS = _LhsSpace, size_t RS = _RhsSpace>
    typename std::enable_if<LO == 2 && RO == 2 && LS == SpaceType::None && RS == SpaceType::None,
                            ExpressionResult<typename Base::Scalar>>::type
    eval(const index_t k) const
    {
        ExpressionResult<typename Base::Scalar> lhs_val = this->lhs_expr_.eval(k);
        ExpressionResult<typename Base::Scalar> rhs_val = this->rhs_expr_.eval(k);
        
        GISMO_ASSERT(lhs_val.rowCardinality() == 1 && lhs_val.colCardinality() == 1,
                    "InnerProductExpression (2,2,None,None): Expected scalar cardinality");
        GISMO_ASSERT(rhs_val.rowCardinality() == 1 && rhs_val.colCardinality() == 1,
                    "InnerProductExpression (2,2,None,None): Expected scalar cardinality");
        
        ExpressionResult<typename Base::Scalar> result(1, 1);
        const auto& lhs_mat = lhs_val(0, 0);
        const auto& rhs_mat = rhs_val(0, 0);
        result(0, 0) = (lhs_mat.array()*rhs_mat.array()).sum();
        return result;
    }

    // --- eval() specialization 7: Matrix ddot (Order=2), Space=Test ---
    template <size_t LO = _LhsOrder, size_t RO = _RhsOrder, size_t LS = _LhsSpace, size_t RS = _RhsSpace>
    typename std::enable_if<LO == 2 && RO == 2 && LS == SpaceType::Test && RS == SpaceType::Test,
                            ExpressionResult<typename Base::Scalar>>::type
    eval(const index_t k) const
    {
        ExpressionResult<typename Base::Scalar> lhs_val = this->lhs_expr_.eval(k);
        ExpressionResult<typename Base::Scalar> rhs_val = this->rhs_expr_.eval(k);
        
        GISMO_ENSURE(lhs_val.rowCardinality() == rhs_val.rowCardinality() &&
                    lhs_val.colCardinality() == rhs_val.colCardinality(),
                    "InnerProductExpression (2,2,Test,Test): Cardinality mismatch");
        
        ExpressionResult<typename Base::Scalar> result(lhs_val.rowCardinality(), lhs_val.colCardinality());
        
        for (index_t i = 0; i < result.rowCardinality(); ++i)
            for (index_t j = 0; j < result.colCardinality(); ++j)
            {
                const auto& lhs_mat = lhs_val(i, j);
                const auto& rhs_mat = rhs_val(i, j);
                result(i, j) = (lhs_mat.array()*rhs_mat.array()).sum();
            }
        
        return result;
    }

    // --- eval() specialization 8: Matrix ddot (Order=2), Space=Trial ---
    template <size_t LO = _LhsOrder, size_t RO = _RhsOrder, size_t LS = _LhsSpace, size_t RS = _RhsSpace>
    typename std::enable_if<LO == 2 && RO == 2 && LS == SpaceType::Trial && RS == SpaceType::Trial,
                            ExpressionResult<typename Base::Scalar>>::type
    eval(const index_t k) const
    {
        ExpressionResult<typename Base::Scalar> lhs_val = this->lhs_expr_.eval(k);
        ExpressionResult<typename Base::Scalar> rhs_val = this->rhs_expr_.eval(k);
        
        GISMO_ENSURE(lhs_val.rowCardinality() == rhs_val.rowCardinality() &&
                    lhs_val.colCardinality() == rhs_val.colCardinality(),
                    "InnerProductExpression (2,2,Trial,Trial): Cardinality mismatch");
        
        ExpressionResult<typename Base::Scalar> result(lhs_val.rowCardinality(), lhs_val.colCardinality());
        
        for (index_t i = 0; i < result.rowCardinality(); ++i)
            for (index_t j = 0; j < result.colCardinality(); ++j)
            {
                const auto& lhs_mat = lhs_val(i, j);
                const auto& rhs_mat = rhs_val(i, j);
                result(i, j) = (lhs_mat.array()*rhs_mat.array()).sum();
            }
        
        return result;
    }

    // --- eval() specialization 9: Matrix ddot (Order=2), Space=Test×Trial ---
    template <size_t LO = _LhsOrder, size_t RO = _RhsOrder, size_t LS = _LhsSpace, size_t RS = _RhsSpace>
    typename std::enable_if<LO == 2 && RO == 2 && LS == SpaceType::Test && RS == SpaceType::Trial,
                            ExpressionResult<typename Base::Scalar>>::type
    eval(const index_t k) const
    {
        ExpressionResult<typename Base::Scalar> lhs_val = this->lhs_expr_.eval(k);
        ExpressionResult<typename Base::Scalar> rhs_val = this->rhs_expr_.eval(k);
        
        // Test has row cardinality, Trial has col cardinality
        ExpressionResult<typename Base::Scalar> result(lhs_val.rowCardinality(), rhs_val.colCardinality());
        
        for (index_t i = 0; i < result.rowCardinality(); ++i)
            for (index_t j = 0; j < result.colCardinality(); ++j)
            {
                const auto& lhs_mat = lhs_val(i, 0);
                const auto& rhs_mat = rhs_val(0, j);
                typename Base::Scalar sum_val = (lhs_mat.array()*rhs_mat.array()).sum();
                result(i, j).resize(1, 1);
                result(i, j)(0, 0) = sum_val;
            }
        
        return result;
    }

    // --- eval() specialization 10: Matrix ddot (Order=2), Space=Trial×Test ---
    template <size_t LO = _LhsOrder, size_t RO = _RhsOrder, size_t LS = _LhsSpace, size_t RS = _RhsSpace>
    typename std::enable_if<LO == 2 && RO == 2 && LS == SpaceType::Trial && RS == SpaceType::Test,
                            ExpressionResult<typename Base::Scalar>>::type
    eval(const index_t k) const
    {
        ExpressionResult<typename Base::Scalar> lhs_val = this->lhs_expr_.eval(k);
        ExpressionResult<typename Base::Scalar> rhs_val = this->rhs_expr_.eval(k);
        
        // Trial has col cardinality, Test has row cardinality
        ExpressionResult<typename Base::Scalar> result(rhs_val.rowCardinality(), lhs_val.colCardinality());
        
        for (index_t i = 0; i < result.rowCardinality(); ++i)
            for (index_t j = 0; j < result.colCardinality(); ++j)
            {
                const auto& lhs_mat = lhs_val(0, j);
                const auto& rhs_mat = rhs_val(i, 0);
                result(i, j) = (lhs_mat.array()*rhs_mat.array()).sum();
            }
        
        return result;
    }

    void print(std::ostream & os) const 
    { 
        if (_LhsOrder == 1 && _RhsOrder == 1)
            os<<this->lhs_expr_<<"⋅"<<this->rhs_expr_;
        else
            os<<this->lhs_expr_<<":"<<this->rhs_expr_;
    }
};

// Generic dot operator to create InnerProductExpression instances using SFINAE
// Vector-vector inner product (None space)
template <typename _LhsExpr, typename _RhsExpr>
typename std::enable_if<
    ((_LhsExpr::Order == 1) &&
    (_RhsExpr::Order == 1)) &&
    ((_LhsExpr::Space == SpaceType::None) && (_RhsExpr::Space == SpaceType::None)),
    InnerProductExpression<_LhsExpr, _RhsExpr, 1, 1, SpaceType::None, SpaceType::None>
>::type
dot(const BaseExpression<_LhsExpr>& lhs, const BaseExpression<_RhsExpr>& rhs)
{
    return InnerProductExpression<_LhsExpr, _RhsExpr, 1, 1, SpaceType::None, SpaceType::None>(lhs, rhs);
}

// Vector-vector inner product (Test space)
template <typename _LhsExpr, typename _RhsExpr>
typename std::enable_if<
    ((_LhsExpr::Order == 1) &&
    (_RhsExpr::Order == 1)) &&
    ((_LhsExpr::Space == SpaceType::Test) && (_RhsExpr::Space == SpaceType::Test)),
    InnerProductExpression<_LhsExpr, _RhsExpr, 1, 1, SpaceType::Test, SpaceType::Test>
>::type
dot(const BaseExpression<_LhsExpr>& lhs, const BaseExpression<_RhsExpr>& rhs)
{
    return InnerProductExpression<_LhsExpr, _RhsExpr, 1, 1, SpaceType::Test, SpaceType::Test>(lhs, rhs);
}

// Vector-vector inner product (Trial space)
template <typename _LhsExpr, typename _RhsExpr>
typename std::enable_if<
    ((_LhsExpr::Order == 1) &&
    (_RhsExpr::Order == 1)) &&
    ((_LhsExpr::Space == SpaceType::Trial) && (_RhsExpr::Space == SpaceType::Trial)),
    InnerProductExpression<_LhsExpr, _RhsExpr, 1, 1, SpaceType::Trial, SpaceType::Trial>
>::type
dot(const BaseExpression<_LhsExpr>& lhs, const BaseExpression<_RhsExpr>& rhs)
{
    return InnerProductExpression<_LhsExpr, _RhsExpr, 1, 1, SpaceType::Trial, SpaceType::Trial>(lhs, rhs);
}

// Vector-vector inner product (Test × Trial)
template <typename _LhsExpr, typename _RhsExpr>
typename std::enable_if<
    ((_LhsExpr::Order == 1) &&
    (_RhsExpr::Order == 1)) &&
    ((_LhsExpr::Space == SpaceType::Test) && (_RhsExpr::Space == SpaceType::Trial)),
    InnerProductExpression<_LhsExpr, _RhsExpr, 1, 1, SpaceType::Test, SpaceType::Trial>
>::type
dot(const BaseExpression<_LhsExpr>& lhs, const BaseExpression<_RhsExpr>& rhs)
{
    return InnerProductExpression<_LhsExpr, _RhsExpr, 1, 1, SpaceType::Test, SpaceType::Trial>(lhs, rhs);
}

// Vector-vector inner product (Trial × Test)
template <typename _LhsExpr, typename _RhsExpr>
typename std::enable_if<
    ((_LhsExpr::Order == 1) &&
    (_RhsExpr::Order == 1)) &&
    ((_LhsExpr::Space == SpaceType::Trial) && (_RhsExpr::Space == SpaceType::Test)),
    InnerProductExpression<_LhsExpr, _RhsExpr, 1, 1, SpaceType::Trial, SpaceType::Test>
>::type
dot(const BaseExpression<_LhsExpr>& lhs, const BaseExpression<_RhsExpr>& rhs)
{
    return InnerProductExpression<_LhsExpr, _RhsExpr, 1, 1, SpaceType::Trial, SpaceType::Test>(lhs, rhs);
}

// `inner` alias for `dot` (vector-vector)
template <typename _LhsExpr, typename _RhsExpr>
auto inner(const BaseExpression<_LhsExpr>& lhs, const BaseExpression<_RhsExpr>& rhs)
-> typename std::enable_if<
    (_LhsExpr::Order == 1) && (_RhsExpr::Order == 1),
    decltype(dot(lhs, rhs))
>::type
{
    return dot(lhs, rhs);
}

// Matrix-matrix inner product (None space)
template <typename _LhsExpr, typename _RhsExpr>
typename std::enable_if<
    (_LhsExpr::Order == 2) &&
    (_RhsExpr::Order == 2) &&
    ((_LhsExpr::Space == SpaceType::None) && (_RhsExpr::Space == SpaceType::None)),
    InnerProductExpression<_LhsExpr, _RhsExpr, 2, 2, SpaceType::None, SpaceType::None>
>::type
ddot(const BaseExpression<_LhsExpr>& lhs, const BaseExpression<_RhsExpr>& rhs)
{
    return InnerProductExpression<_LhsExpr, _RhsExpr, 2, 2, SpaceType::None, SpaceType::None>(lhs, rhs);
}

// Matrix-matrix inner product (Test space)
template <typename _LhsExpr, typename _RhsExpr>
typename std::enable_if<
    (_LhsExpr::Order == 2) &&
    (_RhsExpr::Order == 2) &&
    ((_LhsExpr::Space == SpaceType::Test) && (_RhsExpr::Space == SpaceType::Test)),
    InnerProductExpression<_LhsExpr, _RhsExpr, 2, 2, SpaceType::Test, SpaceType::Test>
>::type
ddot(const BaseExpression<_LhsExpr>& lhs, const BaseExpression<_RhsExpr>& rhs)
{
    return InnerProductExpression<_LhsExpr, _RhsExpr, 2, 2, SpaceType::Test, SpaceType::Test>(lhs, rhs);
}

// Matrix-matrix inner product (Trial space)
template <typename _LhsExpr, typename _RhsExpr>
typename std::enable_if<
    (_LhsExpr::Order == 2) &&
    (_RhsExpr::Order == 2) &&
    ((_LhsExpr::Space == SpaceType::Trial) && (_RhsExpr::Space == SpaceType::Trial)),
    InnerProductExpression<_LhsExpr, _RhsExpr, 2, 2, SpaceType::Trial, SpaceType::Trial>
>::type
ddot(const BaseExpression<_LhsExpr>& lhs, const BaseExpression<_RhsExpr>& rhs)
{
    return InnerProductExpression<_LhsExpr, _RhsExpr, 2, 2, SpaceType::Trial, SpaceType::Trial>(lhs, rhs);
}

// Matrix-matrix inner product (Test × Trial)
template <typename _LhsExpr, typename _RhsExpr>
typename std::enable_if<
    (_LhsExpr::Order == 2) &&
    (_RhsExpr::Order == 2) &&
    ((_LhsExpr::Space == SpaceType::Test) && (_RhsExpr::Space == SpaceType::Trial)),
    InnerProductExpression<_LhsExpr, _RhsExpr, 2, 2, SpaceType::Test, SpaceType::Trial>
>::type
ddot(const BaseExpression<_LhsExpr>& lhs, const BaseExpression<_RhsExpr>& rhs)
{
    return InnerProductExpression<_LhsExpr, _RhsExpr, 2, 2, SpaceType::Test, SpaceType::Trial>(lhs, rhs);
}

// Matrix-matrix inner product (Trial × Test)
template <typename _LhsExpr, typename _RhsExpr>
typename std::enable_if<
    (_LhsExpr::Order == 2) &&
    (_RhsExpr::Order == 2) &&
    ((_LhsExpr::Space == SpaceType::Trial) && (_RhsExpr::Space == SpaceType::Test)),
    InnerProductExpression<_LhsExpr, _RhsExpr, 2, 2, SpaceType::Trial, SpaceType::Test>
>::type
ddot(const BaseExpression<_LhsExpr>& lhs, const BaseExpression<_RhsExpr>& rhs)
{
    return InnerProductExpression<_LhsExpr, _RhsExpr, 2, 2, SpaceType::Trial, SpaceType::Test>(lhs, rhs);
}

// `inner` alias for `ddot` (matrix-matrix)
template <typename _LhsExpr, typename _RhsExpr>
auto inner(const BaseExpression<_LhsExpr>& lhs, const BaseExpression<_RhsExpr>& rhs)
-> typename std::enable_if<
    (_LhsExpr::Order == 2) && (_RhsExpr::Order == 2),
    decltype(ddot(lhs, rhs))
>::type
{
    return ddot(lhs, rhs);
}

// Specialization for NullObject 
template <typename _T, enum SpaceType _LhsSpace, size_t _LhsOrder, typename _RhsExpr>
auto dot(const NullObject<_T, _LhsSpace, _LhsOrder>& /* lhs */, const BaseExpression<_RhsExpr>& /* rhs */)
-> NullObject<_T, SpaceType::None, 0>
{
    return NullObject<_T, SpaceType::None, 0>();
}

template <typename _LhsExpr, typename _T, enum SpaceType _RhsSpace, size_t _RhsOrder>
auto dot(const BaseExpression<_LhsExpr>& /* lhs */, const NullObject<_T, _RhsSpace, _RhsOrder>& /* rhs */)
-> NullObject<_T, SpaceType::None, 0>
{
    return NullObject<_T, SpaceType::None, 0>();
}

// TODO: Implement products also for primitive types

// Specialized dot product for vector · gradient (directional derivative)
// This computes (A·∇)B where A is a vector and ∇B is the gradient matrix
// Mathematically: A^T * ∇B = ∇B^T*A, but we need the result as a vector (not transpose)
template <typename _LhsExpr, typename _RhsExpr>
auto dot(const BaseExpression<_LhsExpr>& lhs, const GradExpression<_RhsExpr, _RhsExpr::Order, _RhsExpr::Space, _RhsExpr::IsConstant>& rhs)
-> decltype(transpose(rhs) * lhs)
{
    return transpose(rhs) * lhs;
}

// Variation of InnerProductExpression where LHS is a SolutionObject and RhsExpr is NOT a SolutionObject
template <typename _LhsExpr, typename _RhsExpr, typename _SpaceObject>
typename std::enable_if<
    !std::is_base_of<SolutionObject<typename _LhsExpr::Scalar,_LhsExpr::Space,_LhsExpr::Order>, _RhsExpr>::value,
    decltype(variation(std::declval<_RhsExpr>(), std::declval<_SpaceObject>()))
>::type
variation(const InnerProductExpression<SolutionObject<typename _LhsExpr::Scalar,_LhsExpr::Space,_LhsExpr::Order>,_RhsExpr,_LhsExpr::Order,_RhsExpr::Order,SolutionObject<typename _LhsExpr::Scalar,_LhsExpr::Space,_LhsExpr::Order>::Order,_RhsExpr::Space> & expr,
          const _SpaceObject & space)
{
    return inner(expr.lhs(), variation(expr.rhs(), space));
}

// Variation of InnerProductExpression where RHS is a SolutionObject and LhsExpr is NOT a SolutionObject
template <typename _LhsExpr, typename _RhsExpr, typename _SpaceObject>
typename std::enable_if<
    !std::is_base_of<SolutionObject<typename _RhsExpr::Scalar,_RhsExpr::Space,_RhsExpr::Order>, _LhsExpr>::value,
    decltype(variation(std::declval<_LhsExpr>(), std::declval<_SpaceObject>()))
>::type
variation(const InnerProductExpression<_LhsExpr,SolutionObject<typename _RhsExpr::Scalar,_RhsExpr::Space,_RhsExpr::Order>,_LhsExpr::Order,_RhsExpr::Order,_LhsExpr::Space,SolutionObject<typename _RhsExpr::Scalar,_RhsExpr::Space,_RhsExpr::Order>::Order> & expr,
          const _SpaceObject & space)
{
    return inner(variation(expr.lhs(), space), expr.rhs());
}

template<typename _LhsExpr, typename _RhsExpr, typename _SpaceObject>
auto variation(const InnerProductExpression<_LhsExpr,_RhsExpr,_LhsExpr::Order,_RhsExpr::Order,_LhsExpr::Space,_RhsExpr::Space> & expr,
               const _SpaceObject & space)
-> decltype(inner(variation(expr.lhs(), space), expr.rhs()) + inner(expr.lhs(), variation(expr.rhs(), space)))
{
    return inner(variation(expr.lhs(), space), expr.rhs()) + inner(expr.lhs(), variation(expr.rhs(), space));
}

}//namespace Expr
}//namespace gismo