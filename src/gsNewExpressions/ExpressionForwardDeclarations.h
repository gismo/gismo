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

namespace Expr
{
    enum SpaceType
    {
        None = 0,
        Test = 1,
        Trial = 2,
        Both = 3
    };

    // Tag dispatch types for C++11 compatibility
    // Use these instead of std::integral_constant for cleaner code
    struct TrueTag {};
    struct FalseTag {};

    // Helper to convert bool to tag type
    template <bool B> struct BoolToTag { typedef FalseTag type; };
    template <> struct BoolToTag<true> { typedef TrueTag type; };

    // Size tag for dispatching on tensor order
    template <size_t N> struct SizeTag { static const size_t value = N; };
}

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

    // Constant Objects
    template <class T, size_t _Order>
    class ConstantObject;

    // VariableObject Objects
    template <class T, size_t _Order, bool _IsConstant>
    class VariableObject;

    // Space Objects
    template <class T, enum SpaceType _Space, size_t _Order>
    class SpaceObject;
    // Null Object
    template <class T, enum SpaceType _Space, size_t _Order>
    class NullObject;
    template <class E, class SpaceType> 
    class VariationObject;

    // Solution Objects
    template <class T, enum SpaceType _Space, size_t _Order>
    class SolutionObject;

    // Component Space Objects
    template <class T, enum SpaceType _Space>
    class ComponentSpaceObject;

    ////////////////////////////////////////////////////////////////
    // Binary Operators
    ////////////////////////////////////////////////////////////////

    template <typename E>
    class BinaryOperator;

    template <typename _LhsExpr, typename _RhsExpr, size_t LhsOrder, size_t RhsOrder, enum SpaceType _LhsSpace, enum SpaceType _RhsSpace>
    class ProductExpression;

    template <typename _LhsExpr, typename _RhsExpr, size_t LhsOrder, size_t RhsOrder, enum SpaceType _LhsSpace, enum SpaceType _RhsSpace>
    class InnerProductExpression;

    template <typename _LhsExpr, typename _RhsExpr, size_t LhsOrder, size_t RhsOrder, enum SpaceType _LhsSpace, enum SpaceType _RhsSpace>
    class OuterProductExpression;

    template <typename _LhsExpr, typename _RhsExpr, size_t LhsOrder, size_t RhsOrder, enum SpaceType _LhsSpace, enum SpaceType _RhsSpace>
    class CrossProductExpression;

    template <typename _LhsExpr, typename _RhsExpr, size_t LhsOrder, size_t RhsOrder, enum SpaceType _LhsSpace, enum SpaceType _RhsSpace>
    class DivisionExpression;

    template <typename _LhsExpr, typename _RhsExpr, size_t LhsOrder, size_t RhsOrder, enum SpaceType _LhsSpace, enum SpaceType _RhsSpace>
    class AddExpression;

    template <typename _LhsExpr, typename _RhsExpr, size_t LhsOrder, size_t RhsOrder, enum SpaceType _LhsSpace, enum SpaceType _RhsSpace>
    class SubtractExpression;

    ////////////////////////////////////////////////////////////////
    // Unary OPERATORS
    ////////////////////////////////////////////////////////////////

    template <typename E>
    class UnaryOperator;

    template <typename E, size_t Order, enum SpaceType _Space, size_t IsConstant>
    class GradExpression;

    template <typename E, size_t Order, enum SpaceType _Space, size_t IsConstant>
    class DivExpression;

    template <typename _E, size_t _Order, enum SpaceType _Space, size_t _IsConstant>
    class CurlExpression;

    template <typename _E, size_t _Order, enum SpaceType _Space, size_t _IsConstant>
    class LaplExpression;

    ////////////////////////////////////////////////////////////////
    // OTHER OPERATORS
    ////////////////////////////////////////////////////////////////

    template <typename E>
    class TransposeExpression;

    template <typename E>
    class ArrayExpression;

    template <typename E, size_t _Order>
    class ComponentExpression;

    template <class T>
    class GeometryMap;

    template <class T>
    class MeasureExpression;

    template <class T>
    class SurfaceNormalExpression;

    template <class T>
    class BoundaryNormalExpression;

    template <class T>
    class BoundaryTangentExpression;

    template <class T>
    class JacobianExpression;

    template <class T>
    class InverseJacobianExpression;

}//namespace Expr
}//namespace gismo
