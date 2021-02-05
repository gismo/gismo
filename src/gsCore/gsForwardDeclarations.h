/** @file gsForwardDeclarations.h

    @brief Provides forward declarations of types and structs.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris
*/

#pragma once

// STD includes
#include <vector>
#include <iterator>
#include <set>
#include <map>
#include <stack>
#include <algorithm>
#include <functional>
#include <limits>

#include <gsCore/gsConfig.h>
#include <gsCore/gsDebug.h>
//#include <gsCore/gsExport.h>  // included by gsMemory.h
#include <gsCore/gsMemory.h>
#include <gsUtils/gsUtils.h>

namespace gismo
{

/**
    \brief Enumeration flags which define the needs of an evaluation object.

    \enum gsNeedEnum

    \ingroup enums
*/
enum gsNeedEnum
{
    NEED_VALUE             = 1U << 0, ///< Value of the object
    NEED_DERIV             = 1U << 1, ///< Gradient of the object
    NEED_GRAD              = NEED_DERIV,
    NEED_JACOBIAN          = 1U << 2, ///< Jacobian of the object
    NEED_MEASURE           = 1U << 3, ///< The density of the measure pull back
    NEED_GRAD_TRANSFORM    = 1U << 4, ///< Gradient transformation matrix
    NEED_DIV               = 1U << 5, ///< Div operator
    NEED_CURL              = 1U << 6, ///< Curl operator
    NEED_DERIV2            = 1U << 7, ///< Second derivatives
    NEED_2ND_DER           = NEED_DERIV2,
    NEED_HESSIAN           = 1U << 8, ///< Hessian matrix
    NEED_LAPLACIAN         = 1U << 9, ///< Laplacian
    NEED_ACTIVE            = 1U <<10, ///< Active function ids
    NEED_NORMAL            = 1U <<11, ///< Normal vector of the object
    NEED_OUTER_NORMAL      = 1U <<12, ///< Outward normal on the boundary

    SAME_ELEMENT           = 1U <<15  ///< Enable optimizations based on the assumption that all evaluation points are in the same bezier domain
};

/** @name ForwardDeclarations
 *  These are forward declarations of classes.
 */

///@{

// exclude from Doxygen processing
/// @cond

//Core objects
template <class T=real_t>                class gsBasis;
template <class T=real_t>                class gsGeometry;
template <class T=real_t>                class gsGeometrySlice;
template <short_t d, class T=real_t>     class gsGenericGeometry;
template <class T=real_t>                class gsConstantBasis;
template <class T=real_t>                class gsBasisFun;

class  gsBoxTopology;
class  boxSide;
struct patchSide;
struct boxCorner;
struct patchCorner;
struct boxComponent;
struct patchComponent;
struct boundaryInterface;

template <class T=real_t>                class gsCurve;
template <class T=real_t>                class gsSurface;
template <class T=real_t>                class gsVolume;
template <class T=real_t>                class gsBulk;

template <class T=real_t>                class gsDomainIterator;

template <class T = real_t, int D=-1>    class gsTensorDomainIterator;

template <class T, int D=-1, class uiter=typename std::vector<T>::const_iterator>
                                         class gsTensorDomainBoundaryIterator;

template <class T=real_t>                class gsDomain;
template <class T=real_t>                class gsFunctionSet;
template <class T=real_t>                class gsFunction;
template <class T=real_t>                class gsFuncCoordinate;
template <class T=real_t>                class gsFuncData;
template <class T=real_t>                class gsMapData;
template <class T=real_t>                class gsFunctionExpr;
template <class T=real_t>                class gsPiecewiseFunction;
template <class T=real_t>                class gsConstantFunction;
template <class T=real_t>                class gsAffineFunction;
template <class T=real_t>                class gsMultiPatch;

// Bases
template <class basis_t >                class gsRationalBasis;
template <short_t d, class T=real_t>     class gsTensorBasis;
template <short_t d, class T=real_t>     class gsHTensorBasis;

template <class T=real_t>                class gsKnotVector;
//template <class T=real_t>              class gsCompactKnotVector;
template <class T=real_t>                class gsBSplineBasis;
template <class T=real_t>                class gsNurbsBasis;
template <short_t d, class T=real_t>     class gsTensorBSplineBasis;
template <short_t d, class T=real_t>     class gsTensorNurbsBasis;
template <short_t d, class T=real_t>     struct gsBSplineTraits;

template <short_t d, class T=real_t>     class gsCompositeIncrSmoothnessBasis;
template <short_t d, class T=real_t>     class gsCompositeGeom;

template <class T=real_t>                class gsBernsteinBasis;
template <short_t d, class T=real_t>     class gsTensorBernsteinBasis;

//template <class T=real_t>              class gsHKnotVector;
template <short_t d, class T=real_t>     class gsHBSplineBasis;
template <short_t d, class T=real_t>     class gsTHBSplineBasis;
template <short_t d, class T=real_t>     class gsTHBSpline;

// Geometries
template <class T=real_t>                class gsBSpline;
template <class T=real_t>                class gsNurbs;
template <class T=real_t>                class gsBezier;
template <short_t d, class T=real_t>     class gsTensorBSpline;
template <short_t d, class T=real_t>     class gsTensorNurbs;
template <short_t d, class T=real_t>     class gsTensorBezier;
template <short_t d, class T=real_t>     class gsHBSpline;
template <class T=real_t>                class gsTrimSurface;

// Quadrature rules
template <class T=real_t>                class gsQuadRule;
template <class T=real_t>                class gsGaussRule;
template <class T=real_t>                class gsGalerkinMethod;

// Domains
// template <class T=real_t>             class gsTensorDomain;
template <short_t d, class T=real_t>     class gsHFitting;

template <class Z, int mode, short_t d=-1,
         bool = //std::is_integral<Z>::value>
         std::numeric_limits<Z>::is_integer && mode!=3>
                                         class gsGridIterator;


// Pde
template <class T=real_t>                class gsPde;
template <class T=real_t>                class gsPoissonPde;
template <class T=real_t>                class gsConvDiffRePde;

template <class T=real_t>                class gsAssembler;
template <class T=real_t>                class gsStokesAssembler;
template <class T=real_t>                class gsGenericAssembler;
template <class T=real_t>                class gsPoissonAssembler;
template <class T=real_t>                class gsCDRAssembler;
template <class T=real_t>                class gsSolverUtils;
template <class T=real_t, bool symm=false>  class gsSparseSystem;

template <class T=real_t>                class gsExprAssembler;
template <class T=real_t>                class gsExprEvaluator;

// More
template <class T=real_t>                class gsCurveLoop;
template <class T=real_t>                class gsPlanarDomain;
template <class T=real_t>                class gsField;
template <class T=real_t>                class gsMesh;
template <class T=real_t>                class gsHeMesh;

template <int d, class T=real_t>         class gsLineSegment;

template <class T=real_t>                class gsFileData;
class gsFileManager;

template <class T=real_t>                class gsSolid;
template <class T=real_t>                class gsSolidVertex;
template <class T=real_t>                class gsSolidHeVertex;
template <class T=real_t>                class gsVolumeBlock;
template <class T=real_t>                class gsSolidHalfEdge;
template <class T=real_t>                class gsSolidHalfFace;
template <class T=real_t>                class gsTriMeshToSolid;

template <class T=real_t>                class gsOptParameterization;
template <class T=real_t>                class gsQualityMeasure;
template <class T=real_t>                class gsCuttingLoop;
template <class T=real_t>                class gsInterpOption;
template <class T=real_t>                class gsMVInterpolation;
template <class T=real_t>                class gsVolumeSegment;
template <class T=real_t>                class gsCompositeTopology;
template <class T=real_t>                class gsBasisEvaluator;
template <class T=real_t>                class gsMultiBasis;

template <class T=real_t>                class gsBemLaplace;
template <class T=real_t>                class gsBemSolution;

template <class T=real_t>                class gsBoundaryConditions;


template <class T=real_t>                class gsVertex;
template <class T=real_t>                class gsCell;
template <class T=real_t>                class gsFace;
template <class T=real_t>                class gsEdge;

template <class T=real_t>                class gsHeVertex;
template <class T=real_t>                class gsHalfFace;
template <class T=real_t>                class gsHalfEdge;
template <class T=real_t>                class gsMeshElement;
template <class T=real_t>                class gsCurvatureSmoothing;
template <class T=real_t>                class gsFitting;


template <class T=real_t>                class gsCurveFitting;

template <class T=real_t>               struct gsNurbsCreator;

template <class T=real_t>               struct gsFieldCreator;

class gsOptionList;

template<class T = real_t, int _Rows=-1, int _Cols=-1,
         int _Options  = 0|((_Rows==1 && _Cols!=1)?0x1:0)> class gsMatrix;
template<class T = real_t, int _Rows=-1, int _Options = 0> class gsVector;

template<class T= real_t, int _Rows=-1, int _Cols=-1> class gsAsConstMatrix;
template<class T= real_t, int _Rows=-1, int _Cols=-1> class gsAsMatrix;

template<class T= real_t, int _Rows=-1> class gsAsVector;
template<class T= real_t, int _Rows=-1> class gsAsConstVector;

template<class T = real_t> class gsVector3d;

template<typename T=real_t, int _Options=0, typename _Index = index_t>
class gsSparseMatrix;

template<typename T=real_t, int _Options=0, typename _Index = index_t>
class gsSparseVector;

template <class T=real_t>                class gsSparseEntries;

// gsSolver

template <class T=real_t>                class gsLinearOperator;
template <class T=real_t>                class gsScaledOp;
template <class T=real_t>                class gsIdentityOp;

template <class T=real_t>                class gsPreconditionerOp;
template <class T=real_t>                class gsPreconditionerFromOp;

template <class T=real_t>                class gsAdditiveOp;
template <class T=real_t>                class gsSumOp;
template <class T=real_t>                class gsProductOp;
template <class T=real_t>                class gsCompositePrecOp;
template <class T=real_t>                class gsKroneckerOp;
template <class T=real_t>                class gsBlockOp;
template <class T=real_t>                class gsPatchPreconditionersCreator;

// gsMultiGrid

template <class T=real_t>                class gsMultiGridOp;
template <class T=real_t>                class gsGridHierarchy;

// gsIeti

template <class T=real_t>                class gsIetiMapper;
template <class T=real_t>                class gsIetiSystem;
template <class T=real_t>                class gsPrimalSystem;
template <class T=real_t>                class gsScaledDirichletPrec;

/// @endcond

///@}

} // end namespace gismo
