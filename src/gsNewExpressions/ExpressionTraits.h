/** @file ExpressionTraits.h

    @brief Traits class for expression types

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

/*
  @brief Traits class for expression types
  
  This class provides compile-time information about an expression type,
  including its scalar type (T), tensor order, space type, and other properties.
  
  Every expression type must specialize this template to define:
  - Scalar: The scalar/element type (real_t, float, complex<real_t>, etc.)
           Unary and binary operators extract this from their operands
  - Order: The tensor order (0=scalar, 1=vector, 2=matrix, etc.)
  - Space: The space type (None, Test, Trial, Both)
  - Deriv: The maximum derivative order required
  - IsConstant: Whether the expression is spatially constant
  
  Example specializations:
  - SpaceObject<T, Space, Order> defines Scalar = T
  - GradExpression<E, ...> defines Scalar = ExpressionTraits<E>::Scalar
  - ProductExpression<L, R, ...> defines Scalar = ExpressionTraits<L>::Scalar
  
  This design allows the framework to work with any scalar type:
  - double (default real_t)
  - float
  - complex<double>
  - Any type supporting the required operations (+, *, norm(), etc.)
*/
template <typename E>
struct ExpressionTraits
{
public:
    // No default Scalar type - every expression must provide one
    // For base expressions (SpaceObject, ConstantObject):
    //   typedef T Scalar;
    // For derived expressions (grad, div, products, etc.):
    //   typedef typename ExpressionTraits<operand>::Scalar Scalar;
    
    static constexpr size_t Order = 0; // Invalid order by default
    static constexpr SpaceType Space = SpaceType::None;
    static constexpr size_t Deriv = 0;
    static constexpr bool IsConstant = false;
};

// Specialization for const types - delegate to non-const version
template <typename E>
struct ExpressionTraits<const E>
{
    typedef typename ExpressionTraits<E>::Scalar Scalar;
    static constexpr size_t Order = ExpressionTraits<E>::Order;
    static constexpr SpaceType Space = ExpressionTraits<E>::Space;
    static constexpr size_t Deriv = ExpressionTraits<E>::Deriv;
    static constexpr bool IsConstant = ExpressionTraits<E>::IsConstant;
};

}//namespace Expr
}//namespace gismo
