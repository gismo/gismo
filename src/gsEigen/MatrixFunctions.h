/** @file MatrixFunctions.h

    @brief Provides matrix functions (exponential, logarithm, square root,
    powers, and the trigonometric/hyperbolic family) as members of
    gsEigen::MatrixBase.

    Entry points, called on any gsMatrix-derived expression: \c .exp()
    (returns a \c gsEigen::MatrixExponentialReturnValue<Derived>), \c .log(),
    \c .sqrt(), \c .pow(p), \c .cos(), \c .sin(), \c .cosh(), \c .sinh(), and
    \c .matrixFunction(f) for a user-supplied scalar function \c f.

    G+Smo compiles Eigen with the token \c Eigen renamed to \c gsEigen
    (see gsCore/gsLinearAlgebra.h), and gsMatrix / gsSparseMatrix derive
    from gsEigen types. Include this header instead of
    <unsupported/Eigen/MatrixFunctions>: a direct include outside the rename
    would place the module in namespace ::Eigen, apart from the gsEigen types
    it has to operate on.

    \warning If a translation unit has already included
    <unsupported/Eigen/MatrixFunctions> directly, that header's include guard
    turns this wrapper into a no-op and the module's symbols stay in ::Eigen.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst
*/

#pragma once

#include <gsCore/gsLinearAlgebra.h>

#define Eigen gsEigen
#define eigen_assert( cond ) GISMO_ASSERT( cond, "" )
#include <unsupported/Eigen/MatrixFunctions>
#undef eigen_assert
#undef Eigen
