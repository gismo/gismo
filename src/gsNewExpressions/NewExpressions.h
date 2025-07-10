/** @file BaseExpression.h

    @brief

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M.Verhelst
*/

#pragma once

#include <gsCore/gsLinearAlgebra.h>
#include <gsNewExpressions/BaseExpression.h>
#include <gsNewExpressions/ScalarExpression.h>
#include <gsNewExpressions/VectorExpression.h>
#include <gsNewExpressions/MatrixExpression.h>
#include <gsNewExpressions/ArrayExpression.h>
#include <gsNewExpressions/AddExpression.h>

namespace gismo
{
namespace Expr
{

/**
 * The new expressions module has the following key components:
 * - BaseExpression: Provides the base declaration of an expression
 * - {Scalar,Vector,Matrix}Expression: represent a {Scalar,Vector,Matrix} defined by constant types (e.g., T, gsVector<T>, gsMatrix<T>)
 * - {Scalar,Vector}VariableExpression: represent a variable defined as a {Scalar,Vector} based on a gsFunction --> gsFeVariable
 * - {Scalar,Vector}SpaceExpression: represent a space defined as a {Scalar,Vector} based on a gsFunction --> gsFeSpace
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