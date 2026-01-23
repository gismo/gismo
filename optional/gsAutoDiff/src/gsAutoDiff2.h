/** @file gsAutoDiff2.h

    @brief Automatic differentiation data type for C++

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst
*/

#pragma once

#include <gsCore/gsConfig.h>
#include <gsCore/gsForwardDeclarations.h>

// Suppress warnings from autodiff library
#if defined(__GNUC__) || defined(__clang__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wparentheses"
#pragma GCC diagnostic ignored "-Wunused-parameter"
#endif

// Include autodiff compatibility layer
#include <gsAutoDiff/gsAutoDiffTraits.h>
#include <gsAutoDiff/gsAutoDiffUtils.h>

#if defined(__GNUC__) || defined(__clang__)
#pragma GCC diagnostic pop
#endif

namespace gismo {
    // Forward AD using autodiff::detail::Dual
    using dual_t = autodiff::detail::Dual<GISMO_COEFF_TYPE, GISMO_COEFF_TYPE>;
    using autodiff_dual_t = dual_t; // For exprtk
    
    // Reverse AD using autodiff::var
    using var_t = autodiff::var;
}

// Provide ostream printing for autodiff scalar types so expression printing works
inline std::ostream & operator<<(std::ostream &os, const autodiff::detail::Dual<GISMO_COEFF_TYPE, GISMO_COEFF_TYPE> &d)
{
    using autodiff::val;
    os << val(d);
    return os;
}

inline std::ostream & operator<<(std::ostream &os, const autodiff::var &v)
{
    using autodiff::val;
    os << val(v);
    return os;
}