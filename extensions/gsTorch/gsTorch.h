/** @file gsTorch.h

    @brief Main include file of gsTorch extension
    
    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): M. Moller
*/

#pragma once

#include <gsTorch/gsTensor.h>
#include <gsTorch/gsTensorOptions.h>
#include <gsTorch/gsTensorWrapper.h>

namespace gismo {
  
    namespace math {

    // LibTorch types cannot be be used as native types in G+Smo and
    // are therefore not included in gsCore/gsMath.h
    // However, to make math functions accessible under the namespace
    // gismo::math::func we import them here.
    using torch::abs;
    using torch::acos;
    using torch::asin;
    using torch::atan2;
    using torch::atan;
    using torch::ceil;
    using torch::cos;
    using torch::cosh;
    using torch::exp;
    using torch::floor;
    //using torch::frexp;
    //using torch::ldexp;
    using torch::log10;
    using torch::log;
    using torch::max;
    using torch::min;
    using torch::pow;
    using torch::sin;
    using torch::sinh;
    using torch::sqrt;
    using torch::tan;
    using torch::tanh;
    using torch::real;
    using torch::imag;
    using torch::conj;

    // dummy
    inline torch::Tensor frexp(const torch::Tensor & a, int* b) { return a; }
    inline torch::Tensor ldexp(const torch::Tensor & a, int* b) { return a; }
    
  } // end namespace math
} // end namespace gismo
