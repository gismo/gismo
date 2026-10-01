/** @file IterativeSolvers.h

    @brief Provides iterative sparse solvers beyond Eigen's core module, in
    namespace gsEigen.

    Entry points, all templated on a gsEigen-compatible sparse matrix type
    (e.g. \c gsSparseMatrix<T>::Base): \c gsEigen::GMRES, \c DGMRES,
    \c MINRES, \c IDRS, \c BiCGSTABL, \c IDRSTABL, and the
    \c gsEigen::IncompleteLU preconditioner.

    G+Smo compiles Eigen with the token \c Eigen renamed to \c gsEigen
    (see gsCore/gsLinearAlgebra.h), and gsMatrix / gsSparseMatrix derive
    from gsEigen types. Include this header instead of
    <unsupported/Eigen/IterativeSolvers>: a direct include outside the rename
    would place the module in namespace ::Eigen, apart from the gsEigen types
    it has to operate on.

    \warning If a translation unit has already included
    <unsupported/Eigen/IterativeSolvers> directly, that header's include guard
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
#include <unsupported/Eigen/IterativeSolvers>
#undef eigen_assert
#undef Eigen
