/** @file gsLinearAlgebra.h

    @brief This is the main header file that collects wrappers of Eigen for linear algebra.

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): A. Mantzaflaris
*/


# pragma once

#include <gsCore/gsForwardDeclarations.h>

#include <gsCore/gsMath.h>

// Eigen linear algebra library (http://eigen.tuxfamily.org)
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>

#ifdef GISMO_WITH_SUPERLU
#include <Eigen/SuperLUSupport>
#endif

namespace gismo
{

using Eigen::internal::cast; // from Core/MathFunctions.h

//Global variables related to gsMatrix ( see also
//gsForwardDeclarations,h)
using Eigen::Dynamic ;//=-1

using Eigen::Lower;
using Eigen::Upper;

// Values for matrix align options
using Eigen::RowMajor;//=0
using Eigen::ColMajor;//=1
using Eigen::AutoAlign;//

template<class T, int _Rows, int _Cols> class gsAsMatrix ;
template<class T, int _Rows, int _Cols> class gsAsConstMatrix ;

template<class T, int _Rows> class gsAsVector ;
template<class T, int _Rows> class gsAsConstVector ;

// helper template for changing the dimension of a matrix
template <int Dim, int Change>
struct ChangeDim
{
    enum { D = Dim + Change };
};
template <int Change>
struct ChangeDim<Dynamic, Change>
{
    enum { D = Dynamic };
};


}; // namespace gismo


#include <gsMatrix/gsMatrixBlockView.h>
#include <gsMatrix/gsMatrix.h>
#include <gsMatrix/gsVector.h>
#include <gsMatrix/gsAsMatrix.h>
#include <gsMatrix/gsSparseMatrix.h>
#include <gsMatrix/gsSparseVector.h>



