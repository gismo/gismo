/** @file OuterProductExpression.h

    @brief Outer product expression class

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

// --- Expression Traits for OuterProductExpression ---
template <typename _LhsExpr, typename _RhsExpr, size_t _LhsOrder, size_t _RhsOrder, size_t _LhsSpace, size_t _RhsSpace>
struct ExpressionTraits<OuterProductExpression<_LhsExpr, _RhsExpr, _LhsOrder, _RhsOrder, _LhsSpace, _RhsSpace>>
{
    typedef _LhsExpr LhsType;
    typedef _RhsExpr RhsType;
    typedef typename ExpressionTraits<_LhsExpr>::Scalar Scalar;
    static constexpr size_t Order = ExpressionTraits<_LhsExpr>::Order + 1;
    static constexpr size_t Space = (SpaceType::None == _LhsSpace && SpaceType::None == _RhsSpace) ? SpaceType::None :
                                    (SpaceType::None != _LhsSpace) ? _LhsSpace : _RhsSpace;
    static constexpr size_t Deriv = (ExpressionTraits<_LhsExpr>::Deriv > ExpressionTraits<_RhsExpr>::Deriv) ? 
                                     ExpressionTraits<_LhsExpr>::Deriv : ExpressionTraits<_RhsExpr>::Deriv;
    static constexpr bool IsConstant = ExpressionTraits<_LhsExpr>::IsConstant && ExpressionTraits<_RhsExpr>::IsConstant;
};

// --- OuterProductExpression with unified class using enable_if on eval() ---

// Unified primary template (handles all Space combinations)
template <typename _LhsExpr, typename _RhsExpr, size_t _LhsOrder, size_t _RhsOrder, size_t _LhsSpace, size_t _RhsSpace>
class OuterProductExpression
 : public BinaryOperator<OuterProductExpression<_LhsExpr, _RhsExpr, _LhsOrder, _RhsOrder, _LhsSpace, _RhsSpace>>
{
    using Base = BinaryOperator<OuterProductExpression<_LhsExpr, _RhsExpr, _LhsOrder, _RhsOrder, _LhsSpace, _RhsSpace>>;
public:
    using Scalar = typename Base::Scalar;
    OuterProductExpression(const _LhsExpr& lhs, const _RhsExpr& rhs) : Base(lhs, rhs)
    {
        this->sizes_[0] = lhs.sizes()[0];
        this->sizes_[1] = rhs.sizes()[0];
    }

    size_t domainDim() const { return this->lhs_expr_.domainDim(); }

    // --- eval() specialization 1: Space=None ---
    template <size_t LS = _LhsSpace, size_t RS = _RhsSpace>
    typename std::enable_if<LS == SpaceType::None && RS == SpaceType::None,
                            ExpressionValue<Scalar>>::type
    eval(const index_t k) const
    {
        ExpressionValue<Scalar> lhs_val = this->lhs_expr_.eval(k);
        ExpressionValue<Scalar> rhs_val = this->rhs_expr_.eval(k);

        GISMO_ASSERT(lhs_val.rowCardinality() == 1 && lhs_val.colCardinality() == 1,
                    "OuterProductExpression (None,None): Expected scalar cardinality");
        GISMO_ASSERT(rhs_val.rowCardinality() == 1 && rhs_val.colCardinality() == 1,
                    "OuterProductExpression (None,None): Expected scalar cardinality");
        
        ExpressionValue<Scalar> result(1, 1);
        gsMatrix<Scalar> lhs_mat = lhs_val(0, 0);
        gsMatrix<Scalar> rhs_mat = rhs_val(0, 0);
        GISMO_ASSERT(lhs_mat.rows()==1 || lhs_mat.cols()==1,"lhs should be a vector");
        GISMO_ASSERT(rhs_mat.rows()==1 || rhs_mat.cols()==1,"rhs should be a vector");
        lhs_mat.resize(lhs_mat.size(),1);
        rhs_mat.resize(rhs_mat.size(),1);
        result(0, 0) = lhs_mat.operator*(rhs_mat.transpose());
        return result;
    }

    // --- eval() specialization 2: Space=Test ---
    template <size_t LS = _LhsSpace, size_t RS = _RhsSpace>
    typename std::enable_if<LS == SpaceType::Test && RS == SpaceType::Test,
                            ExpressionValue<Scalar>>::type
    eval(const index_t k) const
    {
        ExpressionValue<Scalar> lhs_val = this->lhs_expr_.eval(k);
        ExpressionValue<Scalar> rhs_val = this->rhs_expr_.eval(k);

        GISMO_ENSURE(lhs_val.rowCardinality() == rhs_val.rowCardinality() &&
                    lhs_val.colCardinality() == rhs_val.colCardinality(),
                    "OuterProductExpression (Test,Test): Cardinality mismatch");
        
        ExpressionValue<Scalar> result(lhs_val.rowCardinality(), lhs_val.colCardinality());

        for (index_t i = 0; i < result.rowCardinality(); ++i)
            for (index_t j = 0; j < result.colCardinality(); ++j)
            {
                gsMatrix<Scalar> lhs_mat = lhs_val(i, j);
                gsMatrix<Scalar> rhs_mat = rhs_val(i, j);
                GISMO_ASSERT(lhs_mat.rows()==1 || lhs_mat.cols()==1,"lhs should be a vector");
                GISMO_ASSERT(rhs_mat.rows()==1 || rhs_mat.cols()==1,"rhs should be a vector");
                lhs_mat.resize(lhs_mat.size(),1);
                rhs_mat.resize(rhs_mat.size(),1);
                result(i, j) = lhs_mat.operator*(rhs_mat.transpose());
            }

        return result;
    }

    // --- eval() specialization 3: Space=Trial ---
    template <size_t LS = _LhsSpace, size_t RS = _RhsSpace>
    typename std::enable_if<LS == SpaceType::Trial && RS == SpaceType::Trial,
                            ExpressionValue<Scalar>>::type
    eval(const index_t k) const
    {
        ExpressionValue<Scalar> lhs_val = this->lhs_expr_.eval(k);
        ExpressionValue<Scalar> rhs_val = this->rhs_expr_.eval(k);

        GISMO_ENSURE(lhs_val.rowCardinality() == rhs_val.rowCardinality() &&
                    lhs_val.colCardinality() == rhs_val.colCardinality(),
                    "OuterProductExpression (Trial,Trial): Cardinality mismatch");
        
        ExpressionValue<Scalar> result(lhs_val.rowCardinality(), lhs_val.colCardinality());

        for (index_t i = 0; i < result.rowCardinality(); ++i)
            for (index_t j = 0; j < result.colCardinality(); ++j)
            {
                gsMatrix<Scalar> lhs_mat = lhs_val(i, j);
                gsMatrix<Scalar> rhs_mat = rhs_val(i, j);
                GISMO_ASSERT(lhs_mat.rows()==1 || lhs_mat.cols()==1,"lhs should be a vector");
                GISMO_ASSERT(rhs_mat.rows()==1 || rhs_mat.cols()==1,"rhs should be a vector");
                lhs_mat.resize(lhs_mat.size(),1);
                rhs_mat.resize(rhs_mat.size(),1);
                result(i, j) = lhs_mat.operator*(rhs_mat.transpose());
            }

        return result;
    }

    void print(std::ostream & os) const { os<<this->lhs_expr_<<"\u2297"<<this->rhs_expr_; }
};

// --- Function overloads

// Generic outer for None-space expressions (order 1 x order 1)
template<typename _LhsExpr, typename _RhsExpr>
typename std::enable_if<(_LhsExpr::Order == 1) &&
                        (_RhsExpr::Order == 1) &&
                        (_LhsExpr::Space == SpaceType::None) &&
                        (_RhsExpr::Space == SpaceType::None),
                        OuterProductExpression<_LhsExpr, _RhsExpr, 1, 1, SpaceType::None, SpaceType::None>>::type
outer(const BaseExpression<_LhsExpr>& lhs, const BaseExpression<_RhsExpr>& rhs)
{
    return OuterProductExpression<_LhsExpr, _RhsExpr, 1, 1, SpaceType::None, SpaceType::None>(lhs, rhs);
}

// Generic outer for space expressions (Test x Test)
template<typename _LhsExpr, typename _RhsExpr>
typename std::enable_if<(_LhsExpr::Order == 1) &&
                        (_RhsExpr::Order == 1) &&
                        (_LhsExpr::Space == SpaceType::Test) &&
                        (_RhsExpr::Space == SpaceType::Test),
                        OuterProductExpression<_LhsExpr, _RhsExpr, 1, 1, SpaceType::Test, SpaceType::Test>>::type
outer(const BaseExpression<_LhsExpr>& lhs, const BaseExpression<_RhsExpr>& rhs)
{
    return OuterProductExpression<_LhsExpr, _RhsExpr, 1, 1, SpaceType::Test, SpaceType::Test>(lhs, rhs);
}

// Generic outer for space expressions (Trial x Trial)
template<typename _LhsExpr, typename _RhsExpr>
typename std::enable_if<(_LhsExpr::Order == 1) &&
                        (_RhsExpr::Order == 1) &&
                        (_LhsExpr::Space == SpaceType::Trial) &&
                        (_RhsExpr::Space == SpaceType::Trial),
                        OuterProductExpression<_LhsExpr, _RhsExpr, 1, 1, SpaceType::Trial, SpaceType::Trial>>::type
outer(const BaseExpression<_LhsExpr>& lhs, const BaseExpression<_RhsExpr>& rhs)
{
    return OuterProductExpression<_LhsExpr, _RhsExpr, 1, 1, SpaceType::Trial, SpaceType::Trial>(lhs, rhs);
}

// TODO: Implement products also for primitive types


// Specialization for NullObject 
template <typename _T, size_t _LhsSpace, size_t _LhsOrder, typename _RhsExpr>
auto outer(const NullObject<_T, _LhsSpace, _LhsOrder>& /* lhs */, const BaseExpression<_RhsExpr>& /* rhs */)
-> NullObject<_T, SpaceType::None, 0>
{
    return NullObject<_T, SpaceType::None, 0>();
}

template <typename _LhsExpr, typename _T, size_t _RhsSpace, size_t _RhsOrder>
auto outer(const BaseExpression<_LhsExpr>& /* lhs */, const NullObject<_T, _RhsSpace, _RhsOrder>& /* rhs */)
-> NullObject<_T, SpaceType::None, 0>
{
    return NullObject<_T, SpaceType::None, 0>();
}

// Specialization for Expression * Vector primitive type
template <typename _LhsExpr, int _Rows>
auto outer(const BaseExpression<_LhsExpr>& lhs, const gsVector<typename _LhsExpr::Scalar, _Rows>& rhs)
-> OuterProductExpression<_LhsExpr, ConstantObject<typename _LhsExpr::Scalar,1>, _LhsExpr::Order, 1, _LhsExpr::Space, SpaceType::None>
{
    return OuterProductExpression<_LhsExpr, ConstantObject<typename _LhsExpr::Scalar,1>, _LhsExpr::Order, 1, _LhsExpr::Space, SpaceType::None>(lhs, ConstantObject<typename _LhsExpr::Scalar,1>(rhs));
}

// Specialization for Vector primitive type * Expression
template <typename _RhsExpr, int _Rows>
auto outer(const gsVector<typename _RhsExpr::Scalar, _Rows>& lhs, const BaseExpression<_RhsExpr>& rhs)
-> OuterProductExpression<ConstantObject<typename _RhsExpr::Scalar,1>, _RhsExpr, 1, _RhsExpr::Order, SpaceType::None, _RhsExpr::Space>
{
    return OuterProductExpression<ConstantObject<typename _RhsExpr::Scalar,1>, _RhsExpr, 1, _RhsExpr::Order, SpaceType::None, _RhsExpr::Space>(ConstantObject<typename _RhsExpr::Scalar,1>(lhs), rhs);
}

// Variation of OuterProductExpression where LHS is a SolutionObject and RhsExpr is NOT a SolutionObject
template <typename _LhsExpr, typename _RhsExpr, typename _SpaceObject>
typename std::enable_if<
    !std::is_base_of<SolutionObject<typename _LhsExpr::Scalar,_LhsExpr::Space,_LhsExpr::Order>, _RhsExpr>::value,
    decltype(variation(std::declval<_RhsExpr>(), std::declval<_SpaceObject>()))
>::type
variation(const OuterProductExpression<SolutionObject<typename _LhsExpr::Scalar,_LhsExpr::Space,_LhsExpr::Order>,_RhsExpr,_LhsExpr::Order,_RhsExpr::Order,SolutionObject<typename _LhsExpr::Scalar,_LhsExpr::Space,_LhsExpr::Order>::Order,_RhsExpr::Space> & expr,
          const _SpaceObject & space)
{
    return outer(expr.lhs(), variation(expr.rhs(), space));
}

// Variation of OuterProductExpression where RHS is a SolutionObject and LhsExpr is NOT a SolutionObject
template <typename _LhsExpr, typename _RhsExpr, typename _SpaceObject>
typename std::enable_if<
    !std::is_base_of<SolutionObject<typename _RhsExpr::Scalar,_RhsExpr::Space,_RhsExpr::Order>, _LhsExpr>::value,
    decltype(variation(std::declval<_LhsExpr>(), std::declval<_SpaceObject>()))
>::type
variation(const OuterProductExpression<_LhsExpr,SolutionObject<typename _RhsExpr::Scalar,_RhsExpr::Space,_RhsExpr::Order>,_LhsExpr::Order,_RhsExpr::Order,_LhsExpr::Space,SolutionObject<typename _RhsExpr::Scalar,_RhsExpr::Space,_RhsExpr::Order>::Order> & expr,
          const _SpaceObject & space)
{
    return outer(variation(expr.lhs(), space), expr.rhs());
}

template<typename _LhsExpr, typename _RhsExpr, typename _SpaceObject>
auto variation(const OuterProductExpression<_LhsExpr,_RhsExpr,_LhsExpr::Order,_RhsExpr::Order,_LhsExpr::Space,_RhsExpr::Space> & expr,
               const _SpaceObject & space)
-> decltype(outer(variation(expr.lhs(), space), expr.rhs()) + outer(expr.lhs(), variation(expr.rhs(), space)))
{
    return outer(variation(expr.lhs(), space), expr.rhs()) + outer(expr.lhs(), variation(expr.rhs(), space));
}

}//namespace Expr
}//namespace gismo