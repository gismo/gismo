/**
   Forward declarations for automatic differentiation types

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

// Define Eigen macro before including autodiff headers
#define Eigen gsEigen

namespace autodiff {
namespace detail {
   template<typename T, typename G> struct Dual;
}
}

namespace gismo {

namespace math {
namespace constant_autodiff {
   static const double e       =  2.71828182845904523536028747135266249775724709369996;
   static const double pi      =  3.14159265358979323846264338327950288419716939937510;
   static const double pi_2    =  1.57079632679489661923132169163975144209858469968755;
   static const double pi_4    =  0.78539816339744830961566084581987572104929234984378;
   static const double pi_180  =  0.01745329251994329576923690768488612713442871888542;
   static const double _1_pi   =  0.31830988618379067153776752674502872406891929148091;
   static const double _2_pi   =  0.63661977236758134307553505349005744813783858296183;
   static const double _180_pi = 57.29577951308232087679815481410517033240547246656443;
   static const double log2    =  0.69314718055994530941723212145817656807550013436026;
   static const double sqrt2   =  1.41421356237309504880168872420969807856967187537695;
} // namespace constant_autodiff
} // namespace math
} // namespace gismo
