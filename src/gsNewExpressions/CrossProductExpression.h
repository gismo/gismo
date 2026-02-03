/** @file CrossProductExpression.h

    @brief Cross product expression class

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

// --- Expression Traits for CrossProductExpression ---
template <typename _LhsExpr, typename _RhsExpr, size_t _LhsOrder, size_t _RhsOrder, enum SpaceType _LhsSpace, enum SpaceType _RhsSpace>
struct ExpressionTraits<CrossProductExpression<_LhsExpr, _RhsExpr, _LhsOrder, _RhsOrder, _LhsSpace, _RhsSpace>>
{
    typedef _LhsExpr LhsType;
    typedef _RhsExpr RhsType;
    typedef typename ExpressionTraits<_LhsExpr>::Scalar Scalar;
    static constexpr size_t Order = ExpressionTraits<_LhsExpr>::Order;
    static constexpr SpaceType Space = (_LhsSpace != SpaceType::None && _RhsSpace != SpaceType::None && _LhsSpace != _RhsSpace) 
                                   ? SpaceType::Both
                                   : ((_LhsSpace != SpaceType::None) ? _LhsSpace : _RhsSpace);
    static constexpr size_t Deriv = (ExpressionTraits<_LhsExpr>::Deriv > ExpressionTraits<_RhsExpr>::Deriv) ? 
                                     ExpressionTraits<_LhsExpr>::Deriv + 1 : ExpressionTraits<_RhsExpr>::Deriv + 1;
    static constexpr bool IsConstant = ExpressionTraits<_LhsExpr>::IsConstant && ExpressionTraits<_RhsExpr>::IsConstant;
};

// --- CrossProductExpression with unified class using enable_if on eval() ---

// Unified primary template (handles all Space combinations)
template <typename _LhsExpr, typename _RhsExpr, size_t _LhsOrder, size_t _RhsOrder, enum SpaceType _LhsSpace, enum SpaceType _RhsSpace>
class CrossProductExpression
 : public BinaryOperator<CrossProductExpression<_LhsExpr, _RhsExpr, _LhsOrder, _RhsOrder, _LhsSpace, _RhsSpace>>
{
    using Base = BinaryOperator<CrossProductExpression<_LhsExpr, _RhsExpr, _LhsOrder, _RhsOrder, _LhsSpace, _RhsSpace>>;

    static_assert(_LhsSpace == _RhsSpace || 
                  _LhsSpace == SpaceType::None || 
                  _RhsSpace == SpaceType::None ||
                  (_LhsSpace == SpaceType::Test && _RhsSpace == SpaceType::Trial) ||
                  (_LhsSpace == SpaceType::Trial && _RhsSpace == SpaceType::Test),
                  "CrossProductExpression requires compatible spaces: same space, one None, or Test×Trial combination.");
public:
    CrossProductExpression(const _LhsExpr& lhs, const _RhsExpr& rhs) : Base(lhs, rhs)
    {
        // Space consistency checks (runtime to check space ID)
        // If two spaces of the same type are used, they must be identical
        if (_LhsSpace==_RhsSpace)
        {
            if (_LhsSpace == SpaceType::Trial)
                GISMO_ASSERT(lhs.trial()==rhs.trial(),"CrossProductExpression requires both operands to have the same trial space. Probably lhs and rhs are different spaces (e.g., with different ID).");
            else if (_LhsSpace == SpaceType::Test)
                GISMO_ASSERT(lhs.test()==rhs.test(),"CrossProductExpression requires both operands to have the same test space. Probably lhs and rhs are different spaces (e.g., with different ID).");
        }
        else if (_LhsSpace==SpaceType::Both && _RhsSpace==SpaceType::Both)
        {
            GISMO_ASSERT(lhs.test()==rhs.test(),"CrossProductExpression requires both operands to have the same test space. Probably lhs and rhs are different spaces (e.g., with different ID).");
            GISMO_ASSERT(lhs.trial()==rhs.trial(),"CrossProductExpression requires both operands to have the same trial space. Probably lhs and rhs are different spaces (e.g., with different ID).");
        }

        GISMO_ENSURE(lhs.sizes()[0]==3,"lhs must be a vector of size 3");
        GISMO_ENSURE(rhs.sizes()[0]==3,"rhs must be a vector of size 3");
        this->sizes_[0] = 3;
    }

    size_t domainDim() const { return this->lhs_expr_.domainDim(); }

    // --- eval() specialization 1: Space=None ---
    template <size_t LS = _LhsSpace, size_t RS = _RhsSpace>
    typename std::enable_if<LS == SpaceType::None && RS == SpaceType::None,
                            ExpressionResult<typename Base::Scalar>>::type
    eval(const index_t k) const
    {
        ExpressionResult<typename Base::Scalar> lhs_val = this->lhs_expr_.eval(k);
        ExpressionResult<typename Base::Scalar> rhs_val = this->rhs_expr_.eval(k);

        GISMO_ASSERT(lhs_val.rowCardinality() == 1 && lhs_val.colCardinality() == 1,
                    "CrossProductExpression (None,None): Expected scalar cardinality");
        GISMO_ASSERT(rhs_val.rowCardinality() == 1 && rhs_val.colCardinality() == 1,
                    "CrossProductExpression (None,None): Expected scalar cardinality");

        ExpressionResult<typename Base::Scalar> result(1, 1);
        const auto& a = lhs_val(0, 0);
        const auto& b = rhs_val(0, 0);
        gsMatrix<typename Base::Scalar> cross_result(3, 1);
        cross_result(0, 0) = a(1, 0) * b(2, 0) - a(2, 0) * b(1, 0);
        cross_result(1, 0) = a(2, 0) * b(0, 0) - a(0, 0) * b(2, 0);
        cross_result(2, 0) = a(0, 0) * b(1, 0) - a(1, 0) * b(0, 0);
        result(0, 0) = cross_result;
        return result;
    }

    // --- eval() specialization 2: Space=Test ---
    template <size_t LS = _LhsSpace, size_t RS = _RhsSpace>
    typename std::enable_if<LS == SpaceType::Test && RS == SpaceType::Test,
                            ExpressionResult<typename Base::Scalar>>::type
    eval(const index_t k) const
    {
        ExpressionResult<typename Base::Scalar> lhs_val = this->lhs_expr_.eval(k);
        ExpressionResult<typename Base::Scalar> rhs_val = this->rhs_expr_.eval(k);

        GISMO_ENSURE(lhs_val.rowCardinality() == rhs_val.rowCardinality() &&
                    lhs_val.colCardinality() == rhs_val.colCardinality(),
                    "CrossProductExpression (Test,Test): Cardinality mismatch");

        ExpressionResult<typename Base::Scalar> result(lhs_val.rowCardinality(), lhs_val.colCardinality());

        for (index_t i = 0; i < result.rowCardinality(); ++i)
            for (index_t j = 0; j < result.colCardinality(); ++j)
            {
                const auto& a = lhs_val(i, j);
                const auto& b = rhs_val(i, j);
                gsMatrix<typename Base::Scalar> cross_result(3, 1);
                cross_result(0, 0) = a(1, 0) * b(2, 0) - a(2, 0) * b(1, 0);
                cross_result(1, 0) = a(2, 0) * b(0, 0) - a(0, 0) * b(2, 0);
                cross_result(2, 0) = a(0, 0) * b(1, 0) - a(1, 0) * b(0, 0);
                result(i, j) = cross_result;
            }

        return result;
    }

    // --- eval() specialization 3: Space=Trial ---
    template <size_t LS = _LhsSpace, size_t RS = _RhsSpace>
    typename std::enable_if<LS == SpaceType::Trial && RS == SpaceType::Trial,
                            ExpressionResult<typename Base::Scalar>>::type
    eval(const index_t k) const
    {
        ExpressionResult<typename Base::Scalar> lhs_val = this->lhs_expr_.eval(k);
        ExpressionResult<typename Base::Scalar> rhs_val = this->rhs_expr_.eval(k);

        GISMO_ENSURE(lhs_val.rowCardinality() == rhs_val.rowCardinality() &&
                    lhs_val.colCardinality() == rhs_val.colCardinality(),
                    "CrossProductExpression (Trial,Trial): Cardinality mismatch");

        ExpressionResult<typename Base::Scalar> result(lhs_val.rowCardinality(), lhs_val.colCardinality());

        for (index_t i = 0; i < result.rowCardinality(); ++i)
            for (index_t j = 0; j < result.colCardinality(); ++j)
            {
                const auto& a = lhs_val(i, j);
                const auto& b = rhs_val(i, j);
                gsMatrix<typename Base::Scalar> cross_result(3, 1);
                cross_result(0, 0) = a(1, 0) * b(2, 0) - a(2, 0) * b(1, 0);
                cross_result(1, 0) = a(2, 0) * b(0, 0) - a(0, 0) * b(2, 0);
                cross_result(2, 0) = a(0, 0) * b(1, 0) - a(1, 0) * b(0, 0);
                result(i, j) = cross_result;
            }

        return result;
    }

    void print(std::ostream & os) const { os<<this->lhs_expr_<<"×"<<this->rhs_expr_; }

    const SpaceObject<typename Base::Scalar, SpaceType::Trial, Base::Order>& trial() const
    {
        // Return trial from whichever operand has it, or NullObject if neither
        if (_LhsSpace==SpaceType::Trial || _LhsSpace==SpaceType::Both)
            return static_cast<const SpaceObject<typename Base::Scalar, SpaceType::Trial, Base::Order>&>(this->lhs_expr_.trial());
        else
            return static_cast<const SpaceObject<typename Base::Scalar, SpaceType::Trial, Base::Order>&>(this->rhs_expr_.trial());
    }

    const SpaceObject<typename Base::Scalar, SpaceType::Test, Base::Order>& test() const
    {
        // Return test from whichever operand has it, or NullObject if neither
        if (_LhsSpace==SpaceType::Test || _LhsSpace==SpaceType::Both)
            return static_cast<const SpaceObject<typename Base::Scalar, SpaceType::Test, Base::Order>&>(this->lhs_expr_.test());
        else
            return static_cast<const SpaceObject<typename Base::Scalar, SpaceType::Test, Base::Order>&>(this->rhs_expr_.test());
    }
};

// Generic cross operator to create CrossProductExpression instances using SFINAE
// None-space variant
template <typename _LeftExpr, typename _RightExpr>
typename std::enable_if<
    (ExpressionTraits<_LeftExpr>::Order == 1) && (ExpressionTraits<_LeftExpr>::Space == SpaceType::None) &&
    (ExpressionTraits<_RightExpr>::Order == 1) && (ExpressionTraits<_RightExpr>::Space == SpaceType::None),
    CrossProductExpression<_LeftExpr, _RightExpr, 1, 1, SpaceType::None, SpaceType::None>
>::type
cross(const BaseExpression<_LeftExpr>& lhs, const BaseExpression<_RightExpr>& rhs)
{
    return CrossProductExpression<_LeftExpr, _RightExpr, 1, 1, SpaceType::None, SpaceType::None>(lhs, rhs);
}

// Test-space variant
template <typename _LeftExpr, typename _RightExpr>
typename std::enable_if<
    (ExpressionTraits<_LeftExpr>::Order == 1) && (ExpressionTraits<_LeftExpr>::Space == SpaceType::Test) &&
    (ExpressionTraits<_RightExpr>::Order == 1) && (ExpressionTraits<_RightExpr>::Space == SpaceType::Test),
    CrossProductExpression<_LeftExpr, _RightExpr, 1, 1, SpaceType::Test, SpaceType::Test>
>::type
cross(const BaseExpression<_LeftExpr>& lhs, const BaseExpression<_RightExpr>& rhs)
{
    return CrossProductExpression<_LeftExpr, _RightExpr, 1, 1, SpaceType::Test, SpaceType::Test>(lhs, rhs);
}

// Trial-space variant
template <typename _LeftExpr, typename _RightExpr>
typename std::enable_if<
    (ExpressionTraits<_LeftExpr>::Order == 1) && (ExpressionTraits<_LeftExpr>::Space == SpaceType::Trial) &&
    (ExpressionTraits<_RightExpr>::Order == 1) && (ExpressionTraits<_RightExpr>::Space == SpaceType::Trial),
    CrossProductExpression<_LeftExpr, _RightExpr, 1, 1, SpaceType::Trial, SpaceType::Trial>
>::type
cross(const BaseExpression<_LeftExpr>& lhs, const BaseExpression<_RightExpr>& rhs)
{
    return CrossProductExpression<_LeftExpr, _RightExpr, 1, 1, SpaceType::Trial, SpaceType::Trial>(lhs, rhs);
}

// Specializations for transpose expressions (these should trigger compile errors)
template <typename _LeftExpr, typename _RightExpr>
void cross(const TransposeExpression<_LeftExpr>& lhs, const _RightExpr& rhs)
{
    GISMO_ERROR("Cross product is not defined for vectors with different transposition");
}

template <typename _LeftExpr, typename _RightExpr>
void cross(const _LeftExpr& lhs, const TransposeExpression<_RightExpr>& rhs)
{
    GISMO_ERROR("Cross product is not defined for vectors with different transposition");
}

template <typename _LeftExpr, typename _RightExpr>
auto cross(const TransposeExpression<_LeftExpr>& lhs, const TransposeExpression<_RightExpr>& rhs)
    -> decltype(cross(lhs.expr(), rhs.expr()).transpose())
{
    return cross(lhs.expr(), rhs.expr()).transpose();
}

// Specialization for NullObject 
template <typename _T, enum SpaceType _LhsSpace, size_t _LhsOrder, typename _RhsExpr>
auto cross(const NullObject<_T, _LhsSpace, _LhsOrder>& /* lhs */, const BaseExpression<_RhsExpr>& /* rhs */)
-> NullObject<_T, SpaceType::None, 0>
{
    return NullObject<_T, SpaceType::None, 0>();
}

template <typename _LhsExpr, typename _T, enum SpaceType _RhsSpace, size_t _RhsOrder>
auto cross(const BaseExpression<_LhsExpr>& /* lhs */, const NullObject<_T, _RhsSpace, _RhsOrder>& /* rhs */)
-> NullObject<_T, SpaceType::None, 0>
{
    return NullObject<_T, SpaceType::None, 0>();
}

// Specialization for Expression x Vector primitive
template <typename _LhsExpr, int _Rows>
typename std::enable_if<
    (ExpressionTraits<_LhsExpr>::Order == 1) && (ExpressionTraits<_LhsExpr>::Space == SpaceType::None) &&
    (_Rows == 3),
    CrossProductExpression<_LhsExpr, ConstantObject< typename ExpressionTraits<_LhsExpr>::Scalar, 1>, 1, 1, SpaceType::None, SpaceType::None>
>::type
cross(const BaseExpression<_LhsExpr>& lhs, const gsVector<typename ExpressionTraits<_LhsExpr>::Scalar, _Rows>& rhs)
{
    ConstantObject<typename ExpressionTraits<_LhsExpr>::Scalar, 1> rhs_expr(rhs);
    return CrossProductExpression<_LhsExpr, ConstantObject<typename ExpressionTraits<_LhsExpr>::Scalar, 1>, 1, 1, SpaceType::None, SpaceType::None>(lhs, rhs_expr);
}

// Specialization for Vector primitive x Expression
template <typename _RhsExpr, int _Rows>
typename std::enable_if<
    (ExpressionTraits<_RhsExpr>::Order == 1) && (ExpressionTraits<_RhsExpr>::Space == SpaceType::None) &&
    (_Rows == 3),
    CrossProductExpression<ConstantObject< typename ExpressionTraits<_RhsExpr>::Scalar, 1>, _RhsExpr, 1, 1, SpaceType::None, SpaceType::None>
>::type
cross(const gsVector<typename ExpressionTraits<_RhsExpr>::Scalar, _Rows>& lhs, const BaseExpression<_RhsExpr>& rhs)
{
    ConstantObject<typename ExpressionTraits<_RhsExpr>::Scalar, 1> lhs_expr(lhs);
    return CrossProductExpression<ConstantObject<typename ExpressionTraits<_RhsExpr>::Scalar, 1>, _RhsExpr, 1, 1, SpaceType::None, SpaceType::None>(lhs_expr, rhs);
}

// TODO: Implement variations
// Variation of CrossProductExpression where LHS is a SolutionObject and RhsExpr is NOT a SolutionObject
template <typename _LhsExpr, typename _RhsExpr, typename _SpaceObject>
typename std::enable_if<
    !std::is_base_of<SolutionObject<typename _LhsExpr::Scalar,_LhsExpr::Space,_LhsExpr::Order>, _RhsExpr>::value,
    decltype(variation(std::declval<_LhsExpr>(), std::declval<_SpaceObject>()))
>::type
variation(const CrossProductExpression<SolutionObject<typename _LhsExpr::Scalar,_LhsExpr::Space,_LhsExpr::Order>,_RhsExpr,_LhsExpr::Order,_RhsExpr::Order,SolutionObject<typename _LhsExpr::Scalar,_LhsExpr::Space,_LhsExpr::Order>::Order,_RhsExpr::Space> & expr,
          const _SpaceObject & space)
{
    return cross(variation(expr.lhs(), space), expr.rhs());
}

// Variation of CrossProductExpression where RHS is a SolutionObject and LhsExpr is NOT a SolutionObject
template <typename _LhsExpr, typename _RhsExpr, typename _SpaceObject>
typename std::enable_if<
    !std::is_base_of<SolutionObject<typename _RhsExpr::Scalar,_RhsExpr::Space,_RhsExpr::Order>, _LhsExpr>::value,
    decltype(variation(std::declval<_LhsExpr>(), std::declval<_SpaceObject>()))
>::type
variation(const CrossProductExpression<_LhsExpr,SolutionObject<typename _RhsExpr::Scalar,_RhsExpr::Space,_RhsExpr::Order>,_LhsExpr::Order,_RhsExpr::Order,_LhsExpr::Space,SolutionObject<typename _RhsExpr::Scalar,_RhsExpr::Space,_RhsExpr::Order>::Order> & expr,
          const _SpaceObject & space)
{
    return cross(expr.lhs(), variation(expr.rhs(), space));
}

template<typename _LhsExpr, typename _RhsExpr, typename _SpaceObject>
auto variation(const CrossProductExpression<_LhsExpr,_RhsExpr,_LhsExpr::Order,_RhsExpr::Order,_LhsExpr::Space,_RhsExpr::Space> & expr,
               const _SpaceObject & space)
-> decltype(cross(variation(expr.lhs(), space), expr.rhs()) + cross(expr.lhs(), variation(expr.rhs(), space)))
{
    return cross(variation(expr.lhs(), space), expr.rhs()) + cross(expr.lhs(), variation(expr.rhs(), space));
}


}//namespace Expr
}//namespace gismo