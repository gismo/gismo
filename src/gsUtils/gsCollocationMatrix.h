/** @file gsCollocationMatrix.h

    @brief Provides declaration of routines to compute the collocation matrix.

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): A. Mantzaflaris
*/

#pragma once

#include <gsCore/gsForwardDeclarations.h>

namespace gismo
{

///  Computes the collocation matrix of basis on the points u
template<class T>
void gsCollocationMatrix_into (gsBasis<T> const& basis, gsMatrix<T> const& u, 
                       gsSparseMatrix<T> & res);


///  Computes the collocation matrix of basis on the points u
template<class T>
gsSparseMatrix<T> * gsCollocationMatrix (gsBasis<T> const& basis, gsMatrix<T> const& u);


} // namespace gismo


#ifndef GISMO_HEADERS_ONLY
#include GISMO_HPP_HEADER(gsCollocationMatrix.hpp)
#endif
