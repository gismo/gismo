/** @file gsTrilinosNonLinear.h

    @brief Wrappers for Trilinos NOX/LOCA

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris
*/

#pragma once

#include <gsTrilinos/Vector.h>

namespace gismo
{

namespace trilinos
{

namespace solver
{


struct NonLinearPrivate;

/**

 */
class GISMO_EXPORT NonLinear
{
public:

    // https://trilinos.org/docs/dev/packages/nox/doc/html/nox_epetra_tutorial.html
    NonLinear(const SparseMatrix & A );

    ~NonLinear();

    const Vector & solve( const Vector & b );

    void getSolution(gsVector<> & sol, const int rank = 0) const;

protected:
        virtual void solveProblem() = 0;

protected:
        NonLinearPrivate * my;
};

}// namespace solver
}// namespace trilinos
}// namespace gismo
