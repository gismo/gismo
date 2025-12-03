/**
   Automatic differentiation data type for C++, depends on the Eigen
   linear algebra library.

   Copyright (c) 2012 by Wenzel Jakob. Based on code by Jon Kaldor
   and Eitan Grinspun.

   Modifications for G+Smo, Angelos Mantzaflaris, 2015
   Additional modifications for multi-mode support, H.M. Verhelst, 2025

   This library is free software; you can redistribute it and/or
   modify it under the terms of the GNU Lesser General Public
   License as published by the Free Software Foundation; either
   version 2.1 of the License, or (at your option) any later version.

   This library is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   Lesser General Public License for more details.
*/

#pragma once

#include "gsAutoDiff2_fwd.h"
#include <cmath>
#include <sstream>

// Include forward-mode autodiff types
// Note: Do not include eigen.hpp files here as they pull in Eigen/Core too early
// They will be included from gsAutoDiffEigen.h after Eigen plugins are set up
#include <autodiff/forward/dual.hpp>
// we need real_t and util::to_string behaviour
#include <gsCore/gsForwardDeclarations.h>
#include <type_traits>
#include <utility>

// Reverse-mode autodiff is available via <autodiff/reverse/var.hpp>
// For now, VarAdaptor is not included by default due to namespace compatibility issues
// Users can include gsVarAdaptor.h directly if needed

// Short-hand type aliases for autodiff types
using autodiff::detail::Dual;

namespace autodiff {
   // Helper to extract the base scalar type from autodiff types
   template<class T>
   struct autodiff_real_helper {
      typedef T type;
   };

   // Specialization for Dual numbers - recursively extract the value type
   template<typename T, typename G>
   struct autodiff_real_helper<autodiff::detail::Dual<T, G>> {
      typedef typename autodiff_real_helper<T>::type type;
   };

   template<class T>
   using autodiff_real = typename autodiff_real_helper<T>::type;
}

namespace gismo {
namespace math {

// Bring in std math functions by default
using std::floor;
using std::ceil;
using std::round;
using std::trunc;
using std::abs;
using std::log10;
using std::isinf;
using std::isnan;
using std::isfinite;} // namespace math

} // namespace gismo

// Include implementations for each autodiff type
// This follows the same pattern as exprtk_codi_adaptor.hpp

#define AUTODIFF_TYPE dual
#include "gsAutoDiff2_impl.h"
#undef AUTODIFF_TYPE

// Support for reverse-mode autodiff::var
// NOTE: Reverse-mode AD (var) cannot be precompiled into the library due to namespace
// compatibility issues between autodiff (expects Eigen::) and GISMO (uses gsEigen::).
// 
// To use var in your code:
// 1. Include <autodiff/reverse/var.hpp> directly in your source file
// 2. Use var for scalar computations and element-wise matrix operations
// 3. For gradients, use autodiff::derivatives() with autodiff::reverse::detail::wrt()
//
// See reverseAutoDiff.cpp for complete working examples of:
//   - Optimization with many parameters (var's strength)
//   - Element-wise operations on GISMO matrix data
//   - Gradient computation in O(1) passes regardless of parameter count

// Bring gismo::math functions into global scope for Eigen ADL
// This makes ceil, floor, isinf, isnan, isfinite visible to Eigen
using gismo::math::ceil;
using gismo::math::floor;
using gismo::math::isinf;
using gismo::math::isnan;
using gismo::math::isfinite;

#undef Eigen

// Avoid duplicate definitions on multiple includes
#ifndef GISMO_AUTODIFF_TO_STRING_OVERLOADS
#define GISMO_AUTODIFF_TO_STRING_OVERLOADS
namespace gismo {
namespace util {

template <class C, class = void>
struct has_val : std::false_type {};
template <class C>
struct has_val<C, std::void_t<decltype(std::declval<const C&>().val)>> : std::true_type {};

// Convert autodiff expression-like types to a `real_t` and stringify that.
// Disable for types that are already `real_t` to avoid recursion.
template <typename E, std::enable_if_t<has_val<E>::value && !std::is_same<std::decay_t<E>, real_t>::value, int> = 0>
inline std::string to_string(const E &value)
{
   return util::to_string(static_cast<real_t>(value.val));
}

// ExprPtr special-case (reverse-mode expression pointers)
template <typename T>
inline std::string to_string(const autodiff::reverse::detail::ExprPtr<T> &expr)
{
   if (expr) return util::to_string(static_cast<real_t>(expr->val));
   return std::string("null");
}

} // namespace util
} // namespace gismo
#endif

