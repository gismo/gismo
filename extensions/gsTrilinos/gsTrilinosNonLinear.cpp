/** @file gsTrilinosSolvers.cpp

    @brief Wrappers for Trilinos NOX/LOCA

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris
*/

#include <gsTrilinos/gsTrilinosNonLinear.h>
#include "gsTrilinosHeaders.h"

#include "NOX.H"
#include "NOX_Epetra.H"

namespace gismo
{

namespace trilinos
{

namespace solver
{

struct NonLinearPrivate
{
    //Epetra_LinearProblem Problem;

    Vector solution;
};


NonLinear::NonLinear(const SparseMatrix & A)
: my(new NonLinearPrivate)
{

    
}

NonLinear::~NonLinear()
{
    delete my;
}

const Vector & NonLinear::solve( const Vector & b )
{

    
    return my->solution;
}

void NonLinear::getSolution( gsVector<> & sol, const int rank) const
{
    my->solution.copyTo(sol,rank);
}

};// namespace solver
};// namespace trilinos
};// namespace gismo

