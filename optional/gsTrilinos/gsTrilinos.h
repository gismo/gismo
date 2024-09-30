/** @file gsTrilinos.h

    @brief Headers for trilinos-based objects and functionality

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris

    Links:
    http://p4est.github.io/papers/BangerthBursteddeHeisterEtAl11.pdf
*/

#pragma once

namespace gismo
{

/** @namespace gismo::trilinos

    @brief
    This namespace contains wrappers for the Trilinos library
*/
namespace trilinos { }

} // namespace gismo

#include <gsTrilinos/SparseMatrix.h>
#include <gsTrilinos/Vector.h>
#include <gsTrilinos/Operator.h>
#include <gsTrilinos/gsTrilinosSolvers.h>
#include <gsTrilinos/gsTrilinosNonLinear.h>
#include <gsTrilinos/gsTrilinosEigenProblem.h>
