
# pragma once

#include <gsCore/gsForwardDeclarations.h>

#include <gsCore/gsMath.h>

// Eigen linear algebra library (http://eigen.tuxfamily.org)
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>

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



