/** @file BaseObject.h

    @brief

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M.Verhelst
*/

#pragma once

#include <gsCore/gsLinearAlgebra.h>

#include <gsNewExpressions/ExpressionUtils.h>
#include <gsNewExpressions/ExpressionValue.h>
#include <gsNewExpressions/ExpressionForwardDeclarations.h>

#include <gsNewExpressions/BaseExpression.h>
#include <gsNewExpressions/UnaryOperator.h>
#include <gsNewExpressions/BinaryOperator.h>

#include <gsNewExpressions/BaseObject.h>

#include <gsNewExpressions/NullObject.h>
#include <gsNewExpressions/ConstantObject.h>
#include <gsNewExpressions/VariableObject.h>
#include <gsNewExpressions/SpaceObject.h>
#include <gsNewExpressions/VariationObject.h>
#include <gsNewExpressions/SolutionObject.h>

#include <gsNewExpressions/AddExpression.h>
#include <gsNewExpressions/SubtractExpression.h>

#include <gsNewExpressions/ProductExpression.h>
#include <gsNewExpressions/InnerProductExpression.h>
#include <gsNewExpressions/OuterProductExpression.h>
#include <gsNewExpressions/CrossProductExpression.h>

#include <gsNewExpressions/DivisionExpression.h>

#include <gsNewExpressions/GradExpression.h>
#include <gsNewExpressions/DivExpression.h>
#include <gsNewExpressions/CurlExpression.h>
#include <gsNewExpressions/LaplExpression.h>

#include <gsNewExpressions/TransposeExpression.h>

namespace gismo
{
namespace Expr
{

/*
  Traits class for expressions
*/
template <typename E>
struct ExpressionTraits
{
public:
    typedef real_t Scalar;//todo
    static constexpr size_t Order = 0; // Invalid order by default
    static constexpr size_t Space = SpaceType::None;
    static constexpr size_t Deriv = 0;
    static constexpr bool IsConstant = false;
};

// Specialization for const types - delegate to non-const version
template <typename E>
struct ExpressionTraits<const E>
{
    typedef typename ExpressionTraits<E>::Scalar Scalar;
    static constexpr size_t Order = ExpressionTraits<E>::Order;
    static constexpr size_t Space = ExpressionTraits<E>::Space;
    static constexpr size_t Deriv = ExpressionTraits<E>::Deriv;
    static constexpr bool IsConstant = ExpressionTraits<E>::IsConstant;
};



// // Specialization for BaseObject to avoid circular dependency
// template <class T, size_t _Order, bool _IsConstant, short_t _Space>
// struct ExpressionTraits<BaseObject<T, _Order, _IsConstant, _Space>>
// {
// public:
//     typedef T Scalar;
//     typedef const BaseObject<T, _Order, _IsConstant, _Space> Nested_t;

//     static constexpr size_t Order = _Order;
//     static constexpr size_t Space = _Space;
// };

/**
 * The new expressions module has the following key components:
 * - BaseObject: Provides the base declaration of an expression
 * - {Scalar,Vector,Matrix}Expression: represent a {Scalar,Vector,Matrix} defined by constant types (e.g., T, gsVector<T>, gsMatrix<T>)
 * - {Scalar,Vector}VariableObject: represent a VariableObject defined as a {Scalar,Vector} based on a gsFunction --> gsFeVariableObject
 * - {Scalar,Vector}SpaceObject: represent a space defined as a {Scalar,Vector} based on a gsFunction --> gsFeSpace
 *
 *
 * Design questions:
 * - Would there be a way to avoid {Scalar,Vector,Matrix} classes and have just one?
 * - Should we implement a parser or evaluation class keeping track of the derivatives requested per function? Can we keep a maxDeriv member keeping track of this?
 * -
 */

    // template <class T>
    // ScalarExpression<T,>






}//namespace Expr
}//namespace gismo