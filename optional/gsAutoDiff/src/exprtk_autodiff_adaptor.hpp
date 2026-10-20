/** @file exprtk_autodiff_adaptor.hpp

    @brief Provides an exprtk adaptor for autodiff
    arithmetic types

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H. Verhelst
*/

#pragma once

#include <string>
#include <gsAutoDiff/gsAutoDiff.h>

#include <gsAutoDiff/exprtk_autodiff_forward.hpp>
#include "exprtk.hpp"

namespace exprtk
{
namespace details
{

namespace constant_autodiff
{
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

namespace autodiff_helper
{
// True if every component of the forward-mode number x is zero
template <typename U>
inline bool is_zero(const U& x)
{
    if constexpr (autodiff::detail::isDual<U>)
        return is_zero(x.val) && is_zero(x.grad);
    else
        return x == U(0);
}

// True if x is a forward-mode number whose derivative parts are all zero,
// i.e. a constant. Reverse-mode variables are never reported as constant.
template <typename U>
inline bool is_constant(const U& x)
{
    if constexpr (autodiff::detail::isDual<U>)
        return is_zero(x.grad) && is_constant(x.val);
    else
        return std::is_arithmetic<U>::value;
}

// x^c for a forward-mode number x and a real constant c, differentiated as
// c * x^(c-1) * x'. Unlike pow(dual, dual), it does not evaluate log(x) or
// c/x, and unlike autodiff's pow(dual, real), whose value is x^(c-1) * x,
// it gives 0^0 = 1. The exponents 0 and 1 are exact, so that x^0 and x^1
// have finite derivatives of every order at x = 0.
template <typename U, typename R>
inline U pow_constant(const U& x, const R& c)
{
    if constexpr (autodiff::detail::isDual<U>)
    {
        if (c == R(0)) return U(1);
        if (c == R(1)) return x;
        U result;
        result.val  = pow_constant(x.val, c);
        result.grad = c * pow_constant(x.val, c - R(1)) * x.grad;
        return result;
    }
    else
    {
        using std::pow;
        return pow(x, c);
    }
}
} // namespace autodiff_helper
}
}

#define AUTODIFF_TYPE autodiff_dual_t
#include "exprtk_autodiff_adaptor.h"
#undef AUTODIFF_TYPE

#define AUTODIFF_TYPE autodiff_dual2nd_t
#include "exprtk_autodiff_adaptor.h"
#undef AUTODIFF_TYPE

#define AUTODIFF_TYPE autodiff_var_t  
#include "exprtk_autodiff_adaptor.h"
#undef AUTODIFF_TYPE
