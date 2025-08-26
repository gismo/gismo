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
#include <gsNewExpressions/ExpressionForwardDeclarations.h>

#include <gsNewExpressions/BaseExpression.h>


#include <gsNewExpressions/BaseObject.h>

#include <gsNewExpressions/NullObject.h>
#include <gsNewExpressions/ConstantObject.h>
#include <gsNewExpressions/VariableObject.h>
#include <gsNewExpressions/SpaceObject.h>

#include <gsNewExpressions/AddExpression.h>

#include <gsNewExpressions/ProductExpression.h>
#include <gsNewExpressions/DotProductExpression.h>
#include <gsNewExpressions/CrossProductExpression.h>

#include <gsNewExpressions/GradExpression.h>
#include <gsNewExpressions/DivExpression.h>

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
    static constexpr size_t order = 0; // Invalid order by default
    static constexpr size_t space = Space::None;
    static constexpr size_t deriv = 0;
    static constexpr bool isConstant = false;
};

// // Specialization for BaseObject to avoid circular dependency
// template <class T, size_t _order, bool _isConstant, short_t _space>
// struct ExpressionTraits<BaseObject<T, _order, _isConstant, _space>>
// {
// public:
//     typedef T Scalar;
//     typedef const BaseObject<T, _order, _isConstant, _space> Nested_t;

//     static constexpr size_t order = _order;
//     static constexpr size_t space = _space;
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