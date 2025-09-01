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
    // Binary Operators
    ////////////////////////////////////////////////////////////////

    template <typename E>
    class BinaryOperator;

    template <typename LhsExpr, typename RhsExpr, size_t LhsOrder, size_t RhsOrder, size_t LhsSpace, size_t RhsSpace>
    class ProductExpression;

    template <typename LhsExpr, typename RhsExpr, size_t LhsOrder, size_t RhsOrder, size_t LhsSpace, size_t RhsSpace>
    class InnerProductExpression;

    template <typename LhsExpr, typename RhsExpr, size_t LhsOrder, size_t RhsOrder, size_t LhsSpace, size_t RhsSpace>
    class OuterProductExpression;

    template <typename LhsExpr, typename RhsExpr, size_t LhsOrder, size_t RhsOrder, size_t LhsSpace, size_t RhsSpace>
    class CrossProductExpression;

    template <typename LhsExpr, typename RhsExpr, size_t LhsOrder, size_t RhsOrder, size_t LhsSpace, size_t RhsSpace>
    class DivisionExpression;

    template <typename LhsExpr, typename RhsExpr, size_t LhsOrder, size_t RhsOrder, size_t LhsSpace, size_t RhsSpace>
    class AddExpression;

    template <typename LhsExpr, typename RhsExpr, size_t LhsOrder, size_t RhsOrder, size_t LhsSpace, size_t RhsSpace>
    class SubtractExpression;

    ////////////////////////////////////////////////////////////////
    // Unary OPERATORS
    ////////////////////////////////////////////////////////////////

    template <typename E>
    class UnaryOperator;

    template <typename E, size_t Order, size_t Space, size_t IsConstant>
    class GradExpression;

    template <typename E, size_t Order, size_t Space, size_t IsConstant>
    class DivExpression;

    template <typename E, size_t Order, size_t Space, size_t IsConstant>
    class CurlExpression;

    template <typename E, size_t Order, size_t Space, size_t IsConstant>
    class LaplExpression;

    ////////////////////////////////////////////////////////////////
    // OTHER OPERATORS
    ////////////////////////////////////////////////////////////////

    template <typename E>
    class TransposeExpression;

    template <typename E>
    class ArrayExpression;


}//namespace Expr
}//namespace gismo
