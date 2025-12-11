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
    template <typename E>
    class BaseObject;

    // Null Object
    template <class T, size_t _Space, size_t _Order>
    class NullObject;

    // Constant Objects
    template <class T, size_t _Order>
    class ConstantObject;

    // VariableObject Objects
    template <class T, size_t _Order, bool _IsConstant>
    class VariableObject;

    // Space Objects
    template <class T, size_t _Space, size_t _Order>
    class SpaceObject;
    template <class E, class SpaceType> 
    class VariationObject;

    // Solution Objects
    template <class T, size_t _Space, size_t _Order>
    class SolutionObject;

    ////////////////////////////////////////////////////////////////
    // Binary Operators
    ////////////////////////////////////////////////////////////////

    template <typename E>
    class BinaryOperator;

    template <typename _LhsExpr, typename _RhsExpr, size_t LhsOrder, size_t RhsOrder, size_t _LhsSpace, size_t _RhsSpace>
    class ProductExpression;

    template <typename _LhsExpr, typename _RhsExpr, size_t LhsOrder, size_t RhsOrder, size_t _LhsSpace, size_t _RhsSpace>
    class InnerProductExpression;

    template <typename _LhsExpr, typename _RhsExpr, size_t LhsOrder, size_t RhsOrder, size_t _LhsSpace, size_t _RhsSpace>
    class OuterProductExpression;

    template <typename _LhsExpr, typename _RhsExpr, size_t LhsOrder, size_t RhsOrder, size_t _LhsSpace, size_t _RhsSpace>
    class CrossProductExpression;

    template <typename _LhsExpr, typename _RhsExpr, size_t LhsOrder, size_t RhsOrder, size_t _LhsSpace, size_t _RhsSpace>
    class DivisionExpression;

    template <typename _LhsExpr, typename _RhsExpr, size_t LhsOrder, size_t RhsOrder, size_t _LhsSpace, size_t _RhsSpace>
    class AddExpression;

    template <typename _LhsExpr, typename _RhsExpr, size_t LhsOrder, size_t RhsOrder, size_t _LhsSpace, size_t _RhsSpace>
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

    template <typename _E, size_t _Order, size_t _Space, size_t _IsConstant>
    class CurlExpression;

    template <typename _E, size_t _Order, size_t _Space, size_t _IsConstant>
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
