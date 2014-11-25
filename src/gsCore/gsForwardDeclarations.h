/** @file gsForwardDeclarations.h

    @brief Provides forward declarations of types and structs.
    
    This file is part of the G+Smo library.
    
    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris
*/

#pragma once

#include <gsCore/gsConfig.h>
#include <gsCore/gsDebug.h>
#include <gsCore/gsMemory.h>


// STD includes
#include <vector>
#include <set>
#include <map>
#include <stack>
#include <algorithm>

namespace gismo 
{

/** 
    \brief Enumeration flags which define the needs of an evaluation object.

    \enum gsNeedEnum
*/
enum gsNeedEnum
{
    NEED_VALUE             = 1U << 0, /// Value of the object

    NEED_NORMAL            = 1U <<16, /// Normal vector of the object
    NEED_OUTER_NORMAL      = 1U <<17, /// Outward normal on the boundary

    NEED_GRAD              = 1U << 1, /// Gradient of the object
    NEED_JACOBIAN          = 1U << 2, /// Jacobian of the object

    NEED_MEASURE           = 1U << 3, /// The density of the measure pull back
    NEED_GRAD_TRANSFORM    = 1U << 4, /// Gradient transformatin matrix
    NEED_DIV               = 1U << 5, /// Div operator
    NEED_CURL              = 1U << 6, /// Curl operator 

    NEED_2ND_DER           = 1U << 7, /// Second derivatives
    NEED_HESSIAN           = 1U << 8, /// Hessian matrix

    NEED_LAPLACIAN         = 1U << 9  /// Laplacian
};

/** @name ForwardDeclarations
 *  These are forward declarations of classes.
 */

///@{ 

// exclude from Doxygen processing
/// @cond

//Core objects
template< class T = real_t>  class gsBasis;
template< class T = real_t>  class gsGeometry;
template<class Basis_t> class gsGenericGeometry;
template< class T = real_t> class gsGeometryEvaluator;

class gsBoxTopology;

template< class T = real_t>  class gsCurve;
template< class T = real_t>  class gsSurface;
template< class T = real_t>  class gsVolume;
template< class T = real_t>  class gsBulk;

template <class T = real_t> class gsDomainIterator;

template<class T = real_t, int D=-1> class gsTensorDomainIterator;

template<class T, int D=-1, typename uiter = typename std::vector<T>::const_iterator>
class gsTensorDomainBoundaryIterator;

template< class T = real_t>  class gsDomain;
template< class T = real_t>  class gsFunction;
template< class T = real_t>  class gsFunctionExpr;
template< class T = real_t>  class gsConstantFunction;
template< class T = real_t>  class gsMultiPatch;
template< class T = real_t>  class gsBVProblem;

// Bases
template< class T = real_t>  class gsKnotVector;
template< class T = real_t>  class gsCompactKnotVector;
template< class T = real_t>  class gsHKnotVector;

template< class T = real_t, class KnotVectorType=gsKnotVector<T> > class gsBSplineBasis;
template< class T = real_t, class KnotVectorType=gsKnotVector<T> > class gsNurbsBasis;
template< class T = real_t>  class gsBernsteinBasis;
//template< class T = real_t>  class gsRatBernsteinBasis;

template<unsigned d, class T = real_t>  class gsCompositeBasis;
template<unsigned d, class T = real_t>  class gsCompositeGeom;

//template<unsigned d, class basis_t > class gsTensorBasis;
template<unsigned d, class T = real_t, class KnotVectorType=gsKnotVector<T> > class gsTensorBSplineBasis;
template<unsigned d, class T = real_t, class KnotVectorType=gsKnotVector<T> > class gsTensorNurbsBasis;
template<unsigned d, class T = real_t> class gsTensorBernsteinBasis;


template<class basis_t > class gsRationalBasis;
template<unsigned d, class T = real_t>  class gsHTensorBasis;
template<unsigned d, class T = real_t>  class gsHBSplineBasis;
template<unsigned d, class T = real_t>  class gsTHBSplineBasis;
template<unsigned d, class T = real_t>  class gsTHBSpline;

// Geometries
template< class T = real_t, class KnotVectorType=gsKnotVector<T> > class gsBSpline;
template< class T = real_t, class KnotVectorType=gsKnotVector<T> >  class gsNurbs;
template< class T = real_t>  class gsBezier;
template<unsigned d, class T = real_t, class KnotVectorType=gsKnotVector<T> > class gsTensorBSpline;
template<unsigned d, class T = real_t, class KnotVectorType = gsKnotVector<T> > class gsTensorNurbs;
template<unsigned d, class T = real_t> class gsTensorBezier;
template<unsigned d, class T = real_t> class gsHBSpline;
template< class T = real_t>  class gsTrimSurface;

// Quadrature rules
template<class T = real_t> class gsQuadRule;
template<class T = real_t> class gsGaussRule;

template<class T = real_t> class gsGalerkinMethod;

// Domains
// template< class T = real_t>  class gsTensorDomain;
template< class T = real_t>   class gsHFitting;

// Pde
template< class T = real_t>  class gsPde;
template< class T = real_t>  class gsPoissonPde;
template< class T = real_t>  class gsStokesAssembler;
template< class T = real_t>  class gsPoissonAssembler;
template< class T = real_t>  class gsSolverUtils;
template< class T = real_t>  class gsInterpolationAssembler;
template< class T = real_t>  class gsSparseSystem;
template< class T = real_t>  class gsConvDiffRePde;

// More
template< class T = real_t>  class gsCurveLoop;
template< class T = real_t>  class gsPlanarDomain;
template< class T = real_t>  class gsField;
template< class T = real_t>  class gsMesh;
template< class T = real_t>  class gsHeMesh;


template< class T = real_t>  class gsSolid;
template< class T = real_t>  class gsSolidVertex;
template< class T = real_t>  class gsSolidHeVertex;
template< class T = real_t>  class gsVolumeBlock;
template <class T = real_t>  class gsSolidHalfEdge;
template <class T = real_t>  class gsSolidHalfFace;
template< class T = real_t>  class gsTriMeshToSolid;

template <class T=real_t>    class gsOptParameterization;
template <class T=real_t>    class gsQualityMeasure;
template <class T=real_t>    class gsCuttingLoop;
template <class T=real_t>    class gsInterpOption;
template <class T=real_t>    class gsMVInterpolation;
template <class T=real_t>    class gsVolumeSegment;
template <class T=real_t>    class gsCompositeTopology;
template <class T=real_t>    class gsBasisEvaluator;
template <class T=real_t>    class gsMFunctionExpr;
template <class T=real_t>    class gsMultiBasis;

template <class T=real_t>    class gsBemLaplace;
template <class T=real_t>    class gsBemSolution;

template <class T=real_t>    class gsBoundaryConditions;


template <class T=real_t >  class gsVertex;
template <class T=real_t >  class gsCell;
template <class T=real_t >  class gsFace;
template <class T=real_t >  class gsEdge;

template <class T=real_t >  class gsHeVertex;
template <class T=real_t >  class gsHalfFace;
template <class T=real_t >  class gsHalfEdge;
template <class T=real_t >  class gsMeshElement;
template <class T=real_t >  class gsCurvatureSmoothing;
template <class T=real_t >  class gsFitting;


template< class T = real_t>  class gsCurveFitting;

template<class T = real_t> struct gsNurbsCreator;

template<class T = real_t> struct gsFieldCreator;


template<class T = real_t, int _Rows=-1, int _Cols=-1, 
         int _Options  = 0 | ( (_Rows==1 && _Cols!=1) ? 0x1 : 0 ) > class gsMatrix;

template<class T= real_t, int _Rows=-1, int _Cols=-1> class gsAsConstMatrix;
template<class T= real_t, int _Rows=-1, int _Cols=-1> class gsAsMatrix;

template<class T = real_t, int _Rows=-1> class gsVector;

template<class T = real_t> class gsVector3d;

template<typename T=real_t, int _Options=0, typename _Index = index_t>
class gsSparseMatrix;

template<typename T=real_t, int _Options=0, typename _Index = index_t>
class gsSparseVector;

/// @endcond

///@}

} // end namespace gismo
