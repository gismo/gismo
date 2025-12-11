/** @file DivisionExpression.h

    @brief Division expression class

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

// --- DivisionExpression with unified class using enable_if on eval() ---


// ExpressionTraits for DivisionExpression
template <typename _LhsExpr, typename _RhsExpr, size_t _LhsOrder, size_t _RhsOrder, size_t _LhsSpace, size_t _RhsSpace>
struct ExpressionTraits<DivisionExpression<_LhsExpr, _RhsExpr, _LhsOrder, _RhsOrder, _LhsSpace, _RhsSpace>>
{
    typedef _LhsExpr LhsType;
    typedef _RhsExpr RhsType;
    typedef typename ExpressionTraits<_LhsExpr>::Scalar Scalar;
    static constexpr size_t Order = _LhsOrder;
    // Space logic: prefer non-None space; if both different non-None, use Both
    static constexpr size_t Space = (_LhsSpace != SpaceType::None && _RhsSpace != SpaceType::None && _LhsSpace != _RhsSpace) 
                                   ? static_cast<size_t>(SpaceType::Both)
                                   : ((_LhsSpace != SpaceType::None) ? _LhsSpace : _RhsSpace);
    
    static constexpr size_t Deriv = (ExpressionTraits<_LhsExpr>::Deriv > ExpressionTraits<_RhsExpr>::Deriv) ? ExpressionTraits<_LhsExpr>::Deriv + 1 : ExpressionTraits<_RhsExpr>::Deriv + 1;
    static constexpr bool IsConstant = ExpressionTraits<_LhsExpr>::IsConstant && ExpressionTraits<_RhsExpr>::IsConstant;
};

// Unified primary template (handles all Space combinations)
template <typename _LhsExpr, typename _RhsExpr, size_t _LhsOrder, size_t _RhsOrder, size_t _LhsSpace, size_t _RhsSpace>
class DivisionExpression
 : public BinaryOperator<DivisionExpression<_LhsExpr, _RhsExpr, _LhsOrder, _RhsOrder, _LhsSpace, _RhsSpace>>
{
    using Base = BinaryOperator<DivisionExpression<_LhsExpr, _RhsExpr, _LhsOrder, _RhsOrder, _LhsSpace, _RhsSpace>>;
public:
    DivisionExpression(const _LhsExpr& lhs, const _RhsExpr& rhs) : Base(lhs, rhs)
    {
        // Division: result has same size as left operand
        for (size_t d = 0; d < Base::Order; ++d) {
            this->sizes_[d] = lhs.sizes()[d];
        }
    }

    size_t domainDim() const { return this->lhs_expr_.domainDim(); }

    // --- eval() specialization 1: LhsSpace==None ---
    template <size_t LS = _LhsSpace, size_t RS = _RhsSpace>
    typename std::enable_if<LS == SpaceType::None && RS == SpaceType::None,
                            ExpressionValue<typename Base::Scalar>>::type
    eval(const index_t k) const
    {
        ExpressionValue<typename Base::Scalar> lhs_val = this->lhs_expr_.eval(k);
        ExpressionValue<typename Base::Scalar> rhs_val = this->rhs_expr_.eval(k);
        
        // Both have Space=None, cardinality is (1,1) - no loop needed
        typename Base::Scalar divisor = rhs_val(0, 0)(0, 0);
        ExpressionValue<typename Base::Scalar> result(1, 1);
        result(0, 0) = lhs_val(0, 0).array() / divisor;
        return result;
    }

    // --- eval() specialization 2: LhsSpace!=None ---
    template <size_t LS = _LhsSpace, size_t RS = _RhsSpace>
    typename std::enable_if<LS != SpaceType::None && RS == SpaceType::None,
                            ExpressionValue<typename Base::Scalar>>::type
    eval(const index_t k) const
    {
        ExpressionValue<typename Base::Scalar> lhs_val = this->lhs_expr_.eval(k);
        ExpressionValue<typename Base::Scalar> rhs_val = this->rhs_expr_.eval(k);
        
        // RHS must be scalar (Space=None), extract divisor
        typename Base::Scalar divisor = rhs_val(0, 0)(0, 0);
        
        // LHS has space dependency: loop over all basis functions (no runtime branching on cardinality)
        ExpressionValue<typename Base::Scalar> result(lhs_val.rowCardinality(), lhs_val.colCardinality());
        for (index_t i = 0; i < result.rowCardinality(); ++i)
        {
            for (index_t j = 0; j < result.colCardinality(); ++j)
            {
                result(i, j) = lhs_val(i, j).array() / divisor;
            }
        }
        
        return result;
    }

    // --- eval() specialization 3: RHS Space!=None ---
    template <size_t LS = _LhsSpace, size_t RS = _RhsSpace>
    typename std::enable_if<RS != SpaceType::None && LS == SpaceType::None,
                            ExpressionValue<typename Base::Scalar>>::type
    eval(const index_t k) const
    {
        ExpressionValue<typename Base::Scalar> lhs_val = this->lhs_expr_.eval(k);
        ExpressionValue<typename Base::Scalar> rhs_val = this->rhs_expr_.eval(k);

        // LHS can be scalar or vector/matrix, but space-independent (Space=None)
        typename Base::Scalar numerator = lhs_val(0, 0);

        // RHS has space dependency: loop over all basis functions (no runtime branching on cardinality)
        ExpressionValue<typename Base::Scalar> result(rhs_val.rowCardinality(), rhs_val.colCardinality());
        for (index_t i = 0; i < result.rowCardinality(); ++i)
        {
            for (index_t j = 0; j < result.colCardinality(); ++j)
            {
                typename Base::Scalar divisor = rhs_val(i, j)(0, 0);
                result(i, j) = numerator.array() / divisor;
            }
        }
        
        return result;
    }

    // --- eval() specialization 4: Both LHS and RHS Space!=None ---
    template <size_t LS = _LhsSpace, size_t RS = _RhsSpace>
    typename std::enable_if<LS != SpaceType::None && RS != SpaceType::None,
                            ExpressionValue<typename Base::Scalar>>::type
    eval(const index_t k) const
    {
        ExpressionValue<typename Base::Scalar> lhs_val = this->lhs_expr_.eval(k);
        ExpressionValue<typename Base::Scalar> rhs_val = this->rhs_expr_.eval(k);

        GISMO_ENSURE(lhs_val.rowCardinality() == rhs_val.rowCardinality() &&
                    lhs_val.colCardinality() == rhs_val.colCardinality(),
                    "DivisionExpression (Both,Both): Cardinality mismatch");

        // Both LHS and RHS have space dependency: loop over all basis functions (no runtime branching on cardinality)
        ExpressionValue<typename Base::Scalar> result(lhs_val.rowCardinality(), lhs_val.colCardinality());
        for (index_t i = 0; i < result.rowCardinality(); ++i)
        {
            for (index_t j = 0; j < result.colCardinality(); ++j)
            {
                typename Base::Scalar divisor = rhs_val(i, j)(0, 0);
                result(i, j) = lhs_val(i, j).array() / divisor;
            }
        }
        
        return result;
    }

    void print(std::ostream & os) const { os<<this->lhs_expr_<<"/"<<this->rhs_expr_; }

private:
    using Base::lhs_expr_;
    using Base::rhs_expr_;
};

// Generic operator/ to create DivisionExpression instances using SFINAE
template <typename _LhsExpr, typename _RhsExpr>
typename std::enable_if<_RhsExpr::Order == 0,
                        DivisionExpression<_LhsExpr, _RhsExpr, _LhsExpr::Order, 0, _LhsExpr::Space, _RhsExpr::Space>>::type
operator/(const BaseExpression<_LhsExpr>& lhs, const BaseExpression<_RhsExpr>& rhs)
{
    return DivisionExpression<_LhsExpr, _RhsExpr, _LhsExpr::Order, 0, _LhsExpr::Space, _RhsExpr::Space>(lhs, rhs);
}

// Specialization for NullObject
template <typename _T, size_t _LhsSpace, size_t _LhsOrder, typename _RhsExpr>
auto operator/(const NullObject<_T, _LhsSpace, _LhsOrder>& /* lhs */, const BaseExpression<_RhsExpr>& /* rhs */)
-> NullObject<_T,SpaceType::None,0>
{
    return NullObject<_T,SpaceType::None,0>::get();
}

template <typename _LhsExpr, typename _T, size_t _RhsSpace, size_t _RhsOrder>
auto operator/(const BaseExpression<_LhsExpr>& /* lhs */, const NullObject<_T, _RhsSpace, _RhsOrder>& /* rhs */)
-> NullObject<typename _LhsExpr::Scalar,SpaceType::None,0>
{
    static_assert(_RhsOrder != 0, "Division by NullObject of order 0 is undefined.");
    return NullObject<typename _LhsExpr::Scalar,SpaceType::None,0>::get();
}

// Specialization for Expression / Scalar primitive
template <typename _LhsExpr>
auto operator/(const BaseExpression<_LhsExpr>& lhs, const typename _LhsExpr::Scalar & scalar)
-> DivisionExpression<_LhsExpr, ConstantObject<typename _LhsExpr::Scalar,0>, _LhsExpr::Order, 0, _LhsExpr::Space, SpaceType::None>
{
    return DivisionExpression<_LhsExpr, ConstantObject<typename _LhsExpr::Scalar,0>, _LhsExpr::Order, 0, _LhsExpr::Space, SpaceType::None>(lhs, ConstantObject<typename _LhsExpr::Scalar,0>(scalar));
}

// Specialization for Scalar primitive / Expression
template <typename _RhsExpr>
auto operator/(const typename _RhsExpr::Scalar & scalar, const BaseExpression<_RhsExpr>& rhs)
-> DivisionExpression<ConstantObject<typename _RhsExpr::Scalar,0>, _RhsExpr, 0, _RhsExpr::Order, SpaceType::None, _RhsExpr::Space>
{
    return DivisionExpression<ConstantObject<typename _RhsExpr::Scalar,0>, _RhsExpr, 0, _RhsExpr::Order, SpaceType::None, _RhsExpr::Space>(ConstantObject<typename _RhsExpr::Scalar,0>(scalar), rhs);
}

// // Specialization for Scalar primitive / Scalar primitive
// template <typename T>
// auto operator/(const T & lhs, const T & rhs)
// -> DivisionExpression<ConstantObject<T,0>, ConstantObject<T,0>, 0, 0, SpaceType::None, SpaceType::None>
// {
//     return DivisionExpression<ConstantObject<T,0>, ConstantObject<T,0>, 0, 0, SpaceType::None, SpaceType::None>(ConstantObject<T,0>(lhs), ConstantObject<T,0>(rhs));
// }

// Variation of DivisionExpression where LHS is a SolutionObject and RhsExpr is NOT a SolutionObject
template <typename _LhsExpr, typename _RhsExpr, typename _SpaceObject>
typename std::enable_if<
    !std::is_base_of<SolutionObject<typename _LhsExpr::Scalar,_LhsExpr::Space,_LhsExpr::Order>, _RhsExpr>::value,
    decltype(variation(std::declval<_RhsExpr>(), std::declval<_SpaceObject>()))
>::type
variation(const DivisionExpression<SolutionObject<typename _LhsExpr::Scalar,_LhsExpr::Space,_LhsExpr::Order>,_RhsExpr,_LhsExpr::Order,_RhsExpr::Order,SolutionObject<typename _LhsExpr::Scalar,_LhsExpr::Space,_LhsExpr::Order>::Order,_RhsExpr::Space> & expr,
          const _SpaceObject & space)
{
    return (variation(expr.lhs(), space)) / expr.rhs();
}

// Variation of DivisionExpression where RHS is a SolutionObject and LhsExpr is NOT a SolutionObject
template <typename _LhsExpr, typename _RhsExpr, typename _SpaceObject>
typename std::enable_if<
    !std::is_base_of<SolutionObject<typename _RhsExpr::Scalar,_RhsExpr::Space,_RhsExpr::Order>, _LhsExpr>::value,
    decltype(variation(std::declval<_LhsExpr>(), std::declval<_SpaceObject>()))
>::type
variation(const DivisionExpression<_LhsExpr,SolutionObject<typename _RhsExpr::Scalar,_RhsExpr::Space,_RhsExpr::Order>,_LhsExpr::Order,_RhsExpr::Order,_LhsExpr::Space,SolutionObject<typename _RhsExpr::Scalar,_RhsExpr::Space,_RhsExpr::Order>::Order> & expr,
          const _SpaceObject & space)
{
    return (- expr.lhs() * variation(expr.rhs(), space)) / (expr.rhs() * expr.rhs());
}

template<typename _LhsExpr, typename _RhsExpr, typename _SpaceObject>
auto variation(const DivisionExpression<_LhsExpr,_RhsExpr,_LhsExpr::Order,_RhsExpr::Order,_LhsExpr::Space,_RhsExpr::Space> & expr,
               const _SpaceObject & space)
-> decltype((variation(expr.lhs(), space) * expr.rhs() - expr.lhs() * variation(expr.rhs(), space)) / (expr.rhs() * expr.rhs()))
{
    return (variation(expr.lhs(), space) * expr.rhs() - expr.lhs() * variation(expr.rhs(), space)) / (expr.rhs() * expr.rhs());
}

}//namespace Expr
}//namespace gismo