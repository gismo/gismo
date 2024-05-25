/** @file gsOptim.h

    @brief Provides declaration of an optimization problem.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris
*/

#include <gsCore/gsLinearAlgebra.h>
#include <gsOptimizer/gsOptProblem.h>
#include <gsOptimizer/gsOptimizer.h>

#   define Eigen gsEigen
#define OPTIM_ENABLE_EIGEN_WRAPPERS
#include <optim/header_only_version/optim.hpp>
#   undef Eigen

using namespace optim;

#pragma once

namespace gismo
{


/*
 * This file can potentially wrap the root finding solvers of the optim library, see
 * https://optimlib.readthedocs.io/en/latest/api/root_finding_algo_index.html
 */

} // end namespace gismo

// // note: statically compiled in header-only mode
// #ifndef GISMO_BUILD_LIB
// #include GISMO_HPP_HEADER(gsOptim.hpp)
// #endif
