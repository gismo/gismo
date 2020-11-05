/** @file gismo.h

    @brief Main header to be included by clients using the G+Smo library.

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): A. Mantzaflaris
*/

#ifndef __GISMO_H__
#define __GISMO_H__

/** @namespace gismo

    @brief
    The \gismo namespace, containing all definitions for the library.

    \note This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/
namespace gismo
{ 

/** @namespace gismo::internal

    @brief
    This namespace contains functionalities that is internal to the library.
*/
namespace internal 
{ }

}


/* ----------- Core ----------- */
// #include <gsCore/gsConfig.h>
// #include <gsCore/gsDebug.h>
// #include <gsCore/gsExport.h>
// #include <gsCore/gsMemory.h>
#include <gsCore/gsForwardDeclarations.h>
//#include <gsCore/gsJITCompiler.h>

#include <gsCore/gsLinearAlgebra.h>
#include <gsCore/gsFunctionSet.h>
#include <gsCore/gsFuncData.h>
#include <gsCore/gsFunction.h>
#include <gsCore/gsPiecewiseFunction.h>
#include <gsCore/gsBoundary.h>

#include <gsCore/gsGeometry.h>
#include <gsCore/gsGeometrySlice.h>
#include <gsCore/gsCurve.h>
#include <gsCore/gsSurface.h>
#include <gsCore/gsVolume.h>
#include <gsCore/gsBulk.h>
#include <gsCore/gsGenericGeometry.h>

#include <gsCore/gsConstantFunction.h>
#include <gsCore/gsAffineFunction.h>
#include <gsCore/gsFunctionExpr.h>
#include <gsCore/gsBasisFun.h>

#include <gsCore/gsBoxTopology.h>
#include <gsCore/gsMultiPatch.h>
#include <gsCore/gsField.h>

#include <gsCore/gsBasis.h>

#include <gsCore/gsFieldCreator.h>

#include <gsCore/gsDomainIterator.h>

// #include <gsCore/gsTemplateTools.h> // included by gsForwardDeclarations -> gsMemory

// Tensors
#include <gsTensor/gsTensorDomainIterator.h>
#include <gsTensor/gsTensorDomainBoundaryIterator.h>
#include <gsTensor/gsGridIterator.h>
#include <gsTensor/gsGenericTensorBasis.h>

/* ----------- Nurbs ----------- */
#include <gsNurbs/gsKnotVector.h>
#include <gsNurbs/gsBSplineBasis.h>
#include <gsNurbs/gsBSpline.h>
#include <gsNurbs/gsNurbsBasis.h>
#include <gsNurbs/gsNurbs.h>
#include <gsNurbs/gsTensorBSplineBasis.h>
#include <gsNurbs/gsTensorBSpline.h>
#include <gsNurbs/gsTensorNurbsBasis.h>
#include <gsNurbs/gsTensorNurbs.h>
#include <gsNurbs/gsNurbsCreator.h>

/* ----------- HSplines ----------- */
#include <gsHSplines/gsHBSplineBasis.h>
#include <gsHSplines/gsHBSpline.h>
#include <gsHSplines/gsTHBSplineBasis.h>
#include <gsHSplines/gsTHBSpline.h>
#include <gsHSplines/gsHFitting.h>

/* ----------- Modeling ----------- */
#include <gsModeling/gsTrimSurface.h>
#include <gsModeling/gsCurveLoop.h>
#include <gsModeling/gsPlanarDomain.h>
#include <gsModeling/gsSolid.h> 
#include <gsUtils/gsMesh/gsMesh.h>
#include <gsUtils/gsMesh/gsHalfEdgeMesh.h>
#include <gsModeling/gsTriMeshToSolid.h>
//#include <gsSegment/gsVolumeSegment.h> 
#include <gsModeling/gsFitting.h>
#include <gsModeling/gsCurveFitting.h>

#include <gsModeling/gsSpringPatch.h>
#include <gsModeling/gsCoonsPatch.h>
#include <gsModeling/gsCrossApPatch.h>

#include <gsModeling/gsLineSegment.h>
#include <gsModeling/gsParametrization.h>

/* ----------- Pde ----------- */
#include <gsPde/gsBoundaryConditions.h>
#include <gsPde/gsConvDiffRePde.h>
#include <gsPde/gsEulerBernoulliBeamPde.h>
#include <gsPde/gsPoissonPde.h>
#include <gsPde/gsStokesPde.h>
//#include <gsPde/gsNewtonIterator.h>

/* ----------- MultiGrid ----------- */
#include <gsMultiGrid/gsMultiGrid.h>
#include <gsMultiGrid/gsGridHierarchy.h>

/* ----------- Quadrature ----------- */
#include <gsAssembler/gsQuadRule.h>
#include <gsAssembler/gsQuadrature.h>

/* ----------- Assembler ----------- */
#include <gsAssembler/gsAssembler.h>
#include <gsAssembler/gsGenericAssembler.h>
#include <gsAssembler/gsPoissonAssembler.h>
#include <gsAssembler/gsCDRAssembler.h>
#include <gsAssembler/gsHeatEquation.h>

#include <gsAssembler/gsExprHelper.h>
#include <gsAssembler/gsExprAssembler.h>
#include <gsAssembler/gsExprEvaluator.h>

/* ----------- Solver ----------- */
#include <gsSolver/gsLinearOperator.h>
#include <gsSolver/gsMinimalResidual.h>
#include <gsSolver/gsGMRes.h>
#include <gsSolver/gsGradientMethod.h>
#include <gsSolver/gsConjugateGradient.h>
#include <gsSolver/gsPreconditioner.h>
#include <gsSolver/gsAdditiveOp.h>
#include <gsSolver/gsBlockOp.h>
#include <gsSolver/gsCompositePrecOp.h>
#include <gsSolver/gsProductOp.h>
#include <gsSolver/gsSimplePreconditioners.h>
#include <gsSolver/gsSumOp.h>
#include <gsSolver/gsKroneckerOp.h>
#include <gsSolver/gsPatchPreconditionersCreator.h>
#include <gsSolver/gsLanczosMatrix.h>

/* ----------- IO ----------- */
#include <gsIO/gsOptionList.h>
#include <gsIO/gsCmdLine.h>
#include <gsIO/gsFileData.h>
#include <gsIO/gsFileManager.h>
#include <gsIO/gsWriteParaview.h>
#include <gsIO/gsParaviewCollection.h>
#include <gsIO/gsReadFile.h>
#include <gsUtils/gsPointGrid.h>
#include <gsIO/gsXmlUtils.h>

/* ----------- MPI ----------- */
#include <gsMpi/gsMpi.h>

/* ----------- Utilities ----------- */
//#include <gsUtils/gsUtils.h> - in gsForwardDeclarations.h
#include <gsUtils/gsStopwatch.h>
#include <gsUtils/gsFunctionWithDerivatives.h>

/* ----------- Extension ----------- */
#ifdef GISMO_WITH_ADIFF
#include <gsAutoDiff.h>
#endif

/*
#if defined(gismo_EXPORTS) || defined(gismo_dev_EXPORTS)
#  ifdef _MSC_VER
// MSVC and GCC >= 4.4.7
#    pragma message ("warning: The gismo.h header is for clients using the library.")
#  else
// GCC
#    warning "The gismo.h header is for clients using the library."
#  endif
#endif
*/

#endif // __GISMO_H__
