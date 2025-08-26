/** @file ExpressionForwardDeclarations.h

    @brief Forward declarations for all expression classes

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M.Verhelst
*/

#pragma once

#include <gsCore/gsConfig.h>

namespace gismo
{
    // Forward declaration for ExpressionHelper
    template <class T> class ExpressionHelper;

namespace Expr
{
    // Forward declarations for all expression classes

    ////////////////////////////////////////////////////////////////
    // BASE EXPRESSIONS
    ////////////////////////////////////////////////////////////////

    template <typename E>
    struct ExpressionTraits;

    template <typename E>
    class BaseExpression;

    ////////////////////////////////////////////////////////////////
    // OBJECTS
    ////////////////////////////////////////////////////////////////

    // Base Object
    template <class T, size_t _order, bool _isConstant, size_t _space>
    class BaseObject;

    // Null Object
    template <class T, size_t _order>
    class NullObject;

    // Constant Objects
    template <class T, size_t _order>
    class ConstantObject;

    // VariableObject Objects
    template <class T, size_t _order, bool _isConstant>
    class VariableObject;

    // Space Objects
    template <class T, size_t _space, size_t _order>
    class SpaceObject;


    ////////////////////////////////////////////////////////////////
    // PRODUCT EXPRESSIONS
    ////////////////////////////////////////////////////////////////

    template <typename LhsExpr, typename RhsExpr, typename Enable>
    class ProductExpression;

    template <typename LhsExpr, typename RhsExpr, typename Enable>
    class DotProductExpression;

    template <typename LhsExpr, typename RhsExpr, typename Enable>
    class CrossProductExpression;

    ////////////////////////////////////////////////////////////////
    // ADDITION EXPRESSIONS
    ////////////////////////////////////////////////////////////////

    template <typename LhsExpr, typename RhsExpr, typename Enable>
    class AddExpression;

    ////////////////////////////////////////////////////////////////
    // SUBTRACTION EXPRESSIONS
    ////////////////////////////////////////////////////////////////

    template <typename LhsExpr, typename RhsExpr, typename Enable>
    class SubExpression;

    ////////////////////////////////////////////////////////////////
    // DIFFERENTIAL OPERATORS
    ////////////////////////////////////////////////////////////////

    template <typename E, typename Enable>
    class GradExpression;

    template <typename E, typename Enable>
    class DivExpression;

    ////////////////////////////////////////////////////////////////
    // OTHER OPERATORS
    ////////////////////////////////////////////////////////////////

    template <typename E>
    class TransposeExpression;

    template <typename E>
    class ArrayExpression;


}//namespace Expr
}//namespace gismo
