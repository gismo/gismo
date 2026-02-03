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

#include <gsNewExpressions/ExpressionForwardDeclarations.h>
#include <gsNewExpressions/ExpressionUtils.h>
#include <gsNewExpressions/ExpressionResult.h>
#include <gsNewExpressions/ExpressionTraits.h>

#include <gsNewExpressions/BaseExpression.h>
#include <gsNewExpressions/UnaryOperator.h>
#include <gsNewExpressions/BinaryOperator.h>

#include <gsNewExpressions/BaseObject.h>

#include <gsNewExpressions/NullObject.h>
#include <gsNewExpressions/ConstantObject.h>
#include <gsNewExpressions/VariableObject.h>
#include <gsNewExpressions/SpaceObject.h>
#include <gsNewExpressions/ComponentSpaceObject.h>
#include <gsNewExpressions/VariationObject.h>
#include <gsNewExpressions/SolutionObject.h>

#include <gsNewExpressions/AddExpression.h>
#include <gsNewExpressions/SubtractExpression.h>

#include <gsNewExpressions/ProductExpression.h>
#include <gsNewExpressions/InnerProductExpression.h>
#include <gsNewExpressions/OuterProductExpression.h>
#include <gsNewExpressions/CrossProductExpression.h>

#include <gsNewExpressions/DivisionExpression.h>

#include <gsNewExpressions/TransposeExpression.h>
#include <gsNewExpressions/ComponentExpression.h>

#include <gsNewExpressions/GradExpression.h>
#include <gsNewExpressions/DivExpression.h>
#include <gsNewExpressions/CurlExpression.h>
#include <gsNewExpressions/LaplExpression.h>

#include <gsNewExpressions/GeometryMap.h>
#include <gsNewExpressions/MeasureExpression.h>
#include <gsNewExpressions/NormalExpressions.h>
#include <gsNewExpressions/JacobianExpressions.h>
#include <gsNewExpressions/ExpressionVisitor.h>

namespace gismo
{
namespace Expr
{

}//namespace Expr
}//namespace gismo