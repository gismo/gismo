/** @file gsTrilinosHeaders.h

    @brief Header files from Trilinos

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris
*/

#pragma once

#include <gsCore/gsConfig.h>
#include <gsCore/gsTemplateTools.h>

// This defines useful macros like HAVE_MPI, which is defined if and
// only if Epetra was built with MPI enabled.
#include <Epetra_config.h>

#ifdef HAVE_MPI
// Your code is an existing MPI code, so it presumably includes mpi.h directly.
//#  include <mpi.h>

#ifndef GISMO_WITH_MPI
#  ifdef _MSC_VER
// MSVC and GCC >= 4.4.7
#    pragma message ("warning: MPI is disabled.")
#  else
// GCC
#    warning "MPI is disabled."
#  endif
#endif

// Epetra's wrapper for MPI_Comm.  This header file only exists if
// Epetra was built with MPI enabled.
#  include <Epetra_MpiComm.h>
#else
// Wrapper for a "communicator" containing only one process.  This
// header file always exists, whether or not Epetra was built with MPI
// enabled.
#  include <Epetra_SerialComm.h>
#endif // HAVE_MPI

//#include "Epetra_Version.h" // note: contains non-inlined code

#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_FECrsMatrix.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_LinearProblem.h"

#include <Epetra_Export.h>
#include <Epetra_Import.h>

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include <Tpetra_Operator.hpp>
#include <Tpetra_MultiVector_decl.hpp>
#include <Tpetra_CrsMatrix_decl.hpp>
//#include <Tpetra_RowMatrix_decl.hpp>
