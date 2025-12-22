/** @file exprtk_autodiff_forward.hpp

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

// Note: gsAutoDiffEigen.h should already be included by client code BEFORE this header
// to ensure Eigen NumTraits specializations are defined before Eigen is used.
#include <gsAutoDiff/gsAutoDiff2.h>

typedef GISMO_COEFF_TYPE                                            autodiff_real_t;
typedef autodiff::detail::Dual<autodiff_real_t,autodiff_real_t>     autodiff_dual_t;
typedef autodiff::reverse::detail::Variable<autodiff_real_t>        autodiff_var_t;

#define AUTODIFF_TYPE autodiff_dual_t
#include <gsAutoDiff/exprtk_autodiff_forward.h>
#undef AUTODIFF_TYPE

#define AUTODIFF_TYPE autodiff_var_t
#include <gsAutoDiff/exprtk_autodiff_forward.h>
#undef AUTODIFF_TYPE
