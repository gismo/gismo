/** @file gsEigenDeclarations.h

    @brief Extra forward declarations related to the Eigen library

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): A. Mantzaflaris
*/

namespace Eigen 
{

template<typename MatrixType,int RowFactor> class BlockDiag;
template<typename MatrixType,int RowFactor> class BlockTranspose;

namespace internal 
{
template<typename MatrixType> struct adjugate_impl;
}

template<typename XprType, typename IndicesType> class RowSelection;

}
