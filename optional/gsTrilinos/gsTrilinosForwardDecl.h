/** @file gsTrilinosForwardDecl.h

    @brief Forward declarations for Trilinos extensions

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris
*/

class Epetra_Vector;
class Epetra_MultiVector;
class Epetra_CrsMatrix;
class Epetra_FECrsMatrix;
class Epetra_BlockMap;
class Epetra_Operator;
class Epetra_RowMatrix;

namespace Teuchos
{
template<typename O> class RCP;
}
