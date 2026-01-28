/** @file SubtractExpression.h

    @brief Subtraction expression class

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

// --- Generic ExpressionTraits for SubtractExpression ---
template <typename _LhsExpr, typename _RhsExpr, size_t _LhsOrder, size_t _RhsOrder, enum SpaceType _LhsSpace, enum SpaceType _RhsSpace>
struct ExpressionTraits<SubtractExpression<_LhsExpr, _RhsExpr, _LhsOrder, _RhsOrder, _LhsSpace, _RhsSpace>>
{
    typedef _LhsExpr LhsType;
    typedef _RhsExpr RhsType;
    typedef typename ExpressionTraits<_LhsExpr>::Scalar Scalar;
    static constexpr size_t Order = (_LhsOrder > _RhsOrder) ? _LhsOrder : _RhsOrder;
    static constexpr SpaceType Space = (_LhsSpace != SpaceType::None && _RhsSpace != SpaceType::None && _LhsSpace != _RhsSpace) 
                                   ? SpaceType::Both
                                   : ((_LhsSpace != SpaceType::None) ? _LhsSpace : _RhsSpace);
    static constexpr size_t Deriv = (ExpressionTraits<_LhsExpr>::Deriv > ExpressionTraits<_RhsExpr>::Deriv)
                                   ? ExpressionTraits<_LhsExpr>::Deriv : ExpressionTraits<_RhsExpr>::Deriv;
    static constexpr bool IsConstant = ExpressionTraits<_LhsExpr>::IsConstant && ExpressionTraits<_RhsExpr>::IsConstant;
};

// --- Unified SubtractExpression class with enable_if specializations for eval() ---
template <typename _LhsExpr, typename _RhsExpr, size_t _LhsOrder, size_t _RhsOrder, enum SpaceType _LhsSpace, enum SpaceType _RhsSpace>
class SubtractExpression : public BinaryOperator<SubtractExpression<_LhsExpr, _RhsExpr, _LhsOrder, _RhsOrder, _LhsSpace, _RhsSpace>>
{
    using Base = BinaryOperator<SubtractExpression<_LhsExpr, _RhsExpr, _LhsOrder, _RhsOrder, _LhsSpace, _RhsSpace>>;

    static_assert(_LhsSpace == _RhsSpace, "SubtractExpression requires both operands to have the same space.");
    static_assert(_LhsOrder == _RhsOrder, "SubtractExpression requires both operands to have the same order.");

public:
    SubtractExpression(const _LhsExpr& lhs, const _RhsExpr& rhs) : Base(lhs, rhs)
    {
        GISMO_ASSERT(_LhsSpace!=SpaceType::Trial || lhs.trial()==rhs.trial(),"SubtractExpression requires both operands to have the same trial space. Probably lhs and rhs are different spaces (e.g., with different ID).");
        GISMO_ASSERT(_LhsSpace!=SpaceType::Test  || lhs.test()==rhs.test(),"SubtractExpression requires both operands to have the same test space. Probably lhs and rhs are different spaces (e.g., with different ID).");
        
        if (_LhsOrder == _RhsOrder && _LhsOrder > 0) {
            for (size_t d = 0; d != _LhsOrder; ++d) {
                GISMO_ENSURE(this->lhs_expr_.sizes()[d] == this->rhs_expr_.sizes()[d], "SubtractExpression requires same sizes");
                this->sizes_[d] = this->lhs_expr_.sizes()[d];
            }
        } else if (_LhsOrder > _RhsOrder) {
            for (size_t d = 0; d != _LhsOrder; ++d) this->sizes_[d] = this->lhs_expr_.sizes()[d];
        } else if (_RhsOrder > _LhsOrder) {
            for (size_t d = 0; d != _RhsOrder; ++d) this->sizes_[d] = this->rhs_expr_.sizes()[d];
        }
    }

    size_t domainDim() const { return (_LhsOrder >= _RhsOrder) ? this->lhs_expr_.domainDim() : this->rhs_expr_.domainDim(); }

    ExpressionResult<typename Base::Scalar> eval(const index_t k) const
    {
        return this->lhs_expr_.eval(k) - this->rhs_expr_.eval(k);
    }

    void print(std::ostream & os) const { os<<this->lhs_expr_<<"-"<<this->rhs_expr_; }

    // We can take the test space of the LHS expression since both sides must match
    auto test() const -> decltype(std::declval<_LhsExpr>().test())
    {
        return this->lhs_expr_.test();
    }
    // We can take the trial space of the LHS expression since both sides must match
    auto trial() const -> decltype(std::declval<_LhsExpr>().trial())
    {
        return this->lhs_expr_.trial();
    }
};

// Generic operator- to create SubtractExpression instances
template <typename _LhsExpr, typename _RhsExpr>
typename std::enable_if<
    _LhsExpr::Order == _RhsExpr::Order &&
    _LhsExpr::Space == _RhsExpr::Space,
    SubtractExpression<_LhsExpr, _RhsExpr, _LhsExpr::Order, _RhsExpr::Order, _LhsExpr::Space, _RhsExpr::Space>
>::type
operator-(const BaseExpression<_LhsExpr>& lhs, const BaseExpression<_RhsExpr>& rhs)
{
    return SubtractExpression<_LhsExpr, _RhsExpr, _LhsExpr::Order, _RhsExpr::Order, _LhsExpr::Space, _RhsExpr::Space>(lhs, rhs);
}

// Specialization for scalars
template <typename _LhsExpr>
typename std::enable_if<
    !std::is_same<_LhsExpr, ConstantObject<typename ExpressionTraits<_LhsExpr>::Scalar,0>>::value,
    SubtractExpression<_LhsExpr, ConstantObject<typename ExpressionTraits<_LhsExpr>::Scalar,0>, _LhsExpr::Order, 0, _LhsExpr::Space, SpaceType::None>
>::type
operator-(const BaseExpression<_LhsExpr>& lhs, const typename ExpressionTraits<_LhsExpr>::Scalar rhs)
{
    return SubtractExpression<_LhsExpr, ConstantObject<typename ExpressionTraits<_LhsExpr>::Scalar,0>, _LhsExpr::Order, 0, _LhsExpr::Space, SpaceType::None>(lhs, ConstantObject<typename ExpressionTraits<_LhsExpr>::Scalar,0>(rhs));
}

template <typename _RhsExpr>
typename std::enable_if<
    !std::is_same<_RhsExpr, ConstantObject<typename ExpressionTraits<_RhsExpr>::Scalar,0>>::value,
    SubtractExpression<ConstantObject<typename ExpressionTraits<_RhsExpr>::Scalar,0>, _RhsExpr, 0, _RhsExpr::Order, SpaceType::None, _RhsExpr::Space>
>::type
operator-(const typename ExpressionTraits<_RhsExpr>::Scalar lhs, const BaseExpression<_RhsExpr>& rhs)
{
    return SubtractExpression<ConstantObject<typename ExpressionTraits<_RhsExpr>::Scalar,0>, _RhsExpr, 0, _RhsExpr::Order, SpaceType::None, _RhsExpr::Space>(ConstantObject<typename ExpressionTraits<_RhsExpr>::Scalar,0>(lhs), rhs);
}

// Specialization for vectors
template <typename _LhsExpr, int _Rows>
typename std::enable_if<
    !std::is_same<_LhsExpr, ConstantObject<typename ExpressionTraits<_LhsExpr>::Scalar,_LhsExpr::Order>>::value,
    SubtractExpression<_LhsExpr, ConstantObject<typename ExpressionTraits<_LhsExpr>::Scalar,_LhsExpr::Order>, _LhsExpr::Order, _LhsExpr::Order, _LhsExpr::Space, _LhsExpr::Space>
>::type
operator-(const BaseExpression<_LhsExpr>& lhs, const gsVector<typename ExpressionTraits<_LhsExpr>::Scalar,_Rows>& rhs)
{
    return SubtractExpression<_LhsExpr, ConstantObject<typename ExpressionTraits<_LhsExpr>::Scalar,_LhsExpr::Order>, _LhsExpr::Order, _LhsExpr::Order, _LhsExpr::Space, _LhsExpr::Space>(lhs, ConstantObject<typename ExpressionTraits<_LhsExpr>::Scalar,1>(rhs));
}

template <typename _RhsExpr, int _Rows>
typename std::enable_if<
    !std::is_same<_RhsExpr, ConstantObject<typename ExpressionTraits<_RhsExpr>::Scalar,_RhsExpr::Order>>::value,
    SubtractExpression<ConstantObject<typename ExpressionTraits<_RhsExpr>::Scalar,_RhsExpr::Order>, _RhsExpr, _RhsExpr::Order, _RhsExpr::Order, _RhsExpr::Space, _RhsExpr::Space>
>::type
operator-(const gsVector<typename ExpressionTraits<_RhsExpr>::Scalar,_Rows>& lhs, const BaseExpression<_RhsExpr>& rhs)
{
    return SubtractExpression<ConstantObject<typename ExpressionTraits<_RhsExpr>::Scalar,_RhsExpr::Order>, _RhsExpr, _RhsExpr::Order, _RhsExpr::Order, _RhsExpr::Space, _RhsExpr::Space>(ConstantObject<typename ExpressionTraits<_RhsExpr>::Scalar,1>(lhs), rhs);
}

// Specialization for matrices
template <typename _LhsExpr, int _Rows, int _Cols>
typename std::enable_if<
    !std::is_same<_LhsExpr, ConstantObject<typename ExpressionTraits<_LhsExpr>::Scalar,_LhsExpr::Order>>::value,
    SubtractExpression<_LhsExpr, ConstantObject<typename ExpressionTraits<_LhsExpr>::Scalar,_LhsExpr::Order>, _LhsExpr::Order, _LhsExpr::Order, _LhsExpr::Space, _LhsExpr::Space>
>::type
operator-(const BaseExpression<_LhsExpr>& lhs, const gsMatrix<typename ExpressionTraits<_LhsExpr>::Scalar,_Rows,_Cols>& rhs)
{
    return SubtractExpression<_LhsExpr, ConstantObject<typename ExpressionTraits<_LhsExpr>::Scalar,_LhsExpr::Order>, _LhsExpr::Order, _LhsExpr::Order, _LhsExpr::Space, _LhsExpr::Space>(lhs, ConstantObject<typename ExpressionTraits<_LhsExpr>::Scalar,2>(rhs));
}

template <typename _RhsExpr, int _Rows, int _Cols>
typename std::enable_if<
    !std::is_same<_RhsExpr, ConstantObject<typename ExpressionTraits<_RhsExpr>::Scalar,_RhsExpr::Order>>::value,
    SubtractExpression<ConstantObject<typename ExpressionTraits<_RhsExpr>::Scalar,_RhsExpr::Order>, _RhsExpr, _RhsExpr::Order, _RhsExpr::Order, _RhsExpr::Space, _RhsExpr::Space>
>::type
operator-(const gsMatrix<typename ExpressionTraits<_RhsExpr>::Scalar,_Rows,_Cols>& lhs, const BaseExpression<_RhsExpr>& rhs)
{
    return SubtractExpression<ConstantObject<typename ExpressionTraits<_RhsExpr>::Scalar,_RhsExpr::Order>, _RhsExpr, _RhsExpr::Order, _RhsExpr::Order, _RhsExpr::Space, _RhsExpr::Space>(ConstantObject<typename ExpressionTraits<_RhsExpr>::Scalar,2>(lhs), rhs);
}

template <typename _T, enum SpaceType _LhsSpace, size_t _LhsOrder, typename _RhsExpr>
auto operator-(const NullObject<_T, _LhsSpace, _LhsOrder>& /* lhs */, const BaseExpression<_RhsExpr>& rhs)
-> decltype(-rhs)
{
    return -rhs;
}

template <typename _LhsExpr, typename _T, enum SpaceType _RhsSpace, size_t _RhsOrder>
auto operator-(const BaseExpression<_LhsExpr>& lhs, const NullObject<_T, _RhsSpace, _RhsOrder>& /* rhs */)
-> BaseExpression<_LhsExpr>
{
    return lhs;
}


// TODO: avoid ambiguity with (0,0) case
// template <typename _LhsExpr, typename _RhsExpr, typename _SpaceObject>
// SpaceObject<typename LhsExpr::Scalar,SpaceType::None,0>
// variation(const SubtractExpression<LhsExpr,RhsExpr,LhsExpr::Order,RhsExpr::Order,LhsExpr::Space,RhsExpr::Space> & /* expr */,
//           const _SpaceObject & /* space */)
// {
//     gsInfo<<"ZERO VARIATION ADD\n";
//     return NullObject<typename LhsExpr::Scalar,SpaceType::None,0>::get();
// }

// Variation of SubtractExpression where LHS is a SolutionObject and RhsExpr is NOT a SolutionObject
template <typename _LhsExpr, typename _RhsExpr, typename _SpaceObject>
typename std::enable_if<
    !std::is_base_of<SolutionObject<typename _LhsExpr::Scalar,_LhsExpr::Space,_LhsExpr::Order>, _RhsExpr>::value,
    decltype(variation(std::declval<_RhsExpr>(), std::declval<_SpaceObject>()))
>::type
variation(const SubtractExpression<SolutionObject<typename _LhsExpr::Scalar,_LhsExpr::Space,_LhsExpr::Order>,_RhsExpr,_LhsExpr::Order,_RhsExpr::Order,SolutionObject<typename _LhsExpr::Scalar,_LhsExpr::Space,_LhsExpr::Order>::Order,_RhsExpr::Space> & expr,
          const _SpaceObject & space)
{
    return -variation(expr.rhs(), space);
}

// Variation of SubtractExpression where RHS is a SolutionObject and LhsExpr is NOT a SolutionObject
template <typename _LhsExpr, typename _RhsExpr, typename _SpaceObject>
typename std::enable_if<
    !std::is_base_of<SolutionObject<typename _RhsExpr::Scalar,_RhsExpr::Space,_RhsExpr::Order>, _LhsExpr>::value,
    decltype(variation(std::declval<_LhsExpr>(), std::declval<_SpaceObject>()))
>::type
variation(const SubtractExpression<_LhsExpr,SolutionObject<typename _RhsExpr::Scalar,_RhsExpr::Space,_RhsExpr::Order>,_LhsExpr::Order,_RhsExpr::Order,_LhsExpr::Space,SolutionObject<typename _RhsExpr::Scalar,_RhsExpr::Space,_RhsExpr::Order>::Order> & expr,
          const _SpaceObject & space)
{
    return variation(expr.lhs(), space);
}

template<typename _LhsExpr, typename _RhsExpr, typename _SpaceObject>
auto variation(const SubtractExpression<_LhsExpr,_RhsExpr,_LhsExpr::Order,_RhsExpr::Order,_LhsExpr::Space,_RhsExpr::Space> & expr,
               const _SpaceObject & space)
-> decltype(variation(expr.lhs(), space) - variation(expr.rhs(), space))
{
    return variation(expr.lhs(), space) - variation(expr.rhs(), space);
}

}//namespace Expr
}//namespace gismo