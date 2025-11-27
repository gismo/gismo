/**
   Support for reverse-mode automatic differentiation using autodiff::var

   This header provides integration between autodiff's reverse-mode AD (var)
   and G+Smo's matrix types. It should be included AFTER gismo.h in user code.

   Copyright H.M. Verhelst, 2025

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

// Map Eigen namespace to gsEigen for autodiff compatibility
#define Eigen gsEigen

// Reverse-mode autodiff types with Eigen integration
#include <autodiff/reverse/var.hpp>
#include <autodiff/reverse/var/eigen.hpp>  // Eigen integration for var

// Undefine the mapping
#undef Eigen

namespace gismo {

/// Helper to extract the value from a var
inline double val(const autodiff::var& v) {
    return static_cast<double>(v);
}

} // namespace gismo

// Note: gsMatrix<var> and gsVector<var> are now fully supported!
// 
// Usage example:
//   #include <gismo.h>
//   #include <gsAutoDiff/gsAutoDiffVar.h>
//   
//   gsMatrix<var> params(n, 1);
//   // ... initialize params ...
//   var objective = /* compute from params */;
//   auto grad = autodiff::gradient(objective, params);
//
// Key features:
// - Full gsMatrix<var> support via namespace mapping
// - autodiff::gradient() works directly with gsMatrix
// - autodiff::hessian() for second derivatives
// - Efficient: O(1) passes for gradient regardless of parameter count
//
// Current limitations:
// - Some advanced operations (cramerInverse) may not work
// - Use standard matrix operations: +, -, *, transpose, norm, etc.
