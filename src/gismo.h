/** @file gismo.h

    @brief Main header to be included by clients using the G+Smo library.

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): A. Mantzaflaris
*/


#pragma once

/** @namespace gismo

    @brief
    The gismo namespace, containing all definitions for the Gismo library.
*/
namespace gismo
{ }

#include <gsCore/gsConfig.h>
#include <gsCore/gsDebug.h>
#include <gsCore/gsForwardDeclarations.h>

// Core
#include <gsCore/gsLinearAlgebra.h>
#include <gsCore/gsFunctionExpr.h>
#include <gsCore/gsMFunctionExpr.h>
#include <gsCore/gsConstantFunction.h>
#include <gsCore/gsMultiPatch.h>
#include <gsCore/gsField.h>

// Domain iterators
#include <gsCore/gsDomainIterator.h>
#include <gsTensor/gsTensorDomainIterator.h>
#include <gsTensor/gsTensorDomainBoundaryIterator.h>

// Nurbs
#include <gsNurbs/gsKnotVector.h>
#include <gsNurbs/gsCompactKnotVector.h>
#include <gsNurbs/gsBSplineBasis.h>
#include <gsNurbs/gsBSpline.h>
#include <gsNurbs/gsNurbsBasis.h>
#include <gsNurbs/gsNurbs.h>
#include <gsNurbs/gsTensorBSplineBasis.h>
#include <gsNurbs/gsTensorBSpline.h>
#include <gsNurbs/gsTensorNurbsBasis.h>
#include <gsNurbs/gsTensorNurbs.h>
#include <gsNurbs/gsNurbsCreator.h>

// Thbs
#include <gsThbs/gsHBSplineBasis.h>
#include <gsThbs/gsHBSpline.h>
#include <gsThbs/gsTHBSplineBasis.h>
#include <gsThbs/gsTHBSpline.h>
#include <gsThbs/gsHFitting.h>

// Modeling
#include <gsModeling/gsTrimSurface.h>
#include <gsModeling/gsCurveLoop.h>
#include <gsModeling/gsPlanarDomain.h>
#include <gsModeling/gsSolid.h> 
#include <gsUtils/gsMesh/gsMesh.h>
#include <gsModeling/gsTriMeshToSolid.h>
//#include <gsSegment/gsVolumeSegment.h> 
#include <gsModeling/gsFitting.h>
#include <gsModeling/gsCurveFitting.h>

// Pde
#include <gsPde/gsConvDiffRePde.h>
#include <gsPde/gsEulerBernoulliBeamPde.h>
#include <gsPde/gsPoissonPde.h>
#include <gsPde/gsStokesPde.h>
#include <gsPde/gsBVProblem.h>
// Norms
#include <gsAssembler/gsNormL2.h>
#include <gsAssembler/gsNormH1.h>

// Quadrature
#include <gsAssembler/gsQuadRule.h>
#include <gsAssembler/gsGaussRule.h>

// Assembler
#include <gsAssembler/gsAssemblerBase.h>
#include <gsAssembler/gsGenericAssembler.h>
#include <gsAssembler/gsPoissonAssembler.h>

// Solver
#include <gsSolver/gsMinimalResidual.h>
#include <gsSolver/gsPreconditioner.h>

// IO
#include <gsIO/gsCmdLine.h>
#include <gsIO/gsWriteParaview.h>
#include <gsIO/gsReadFile.h>
#include <gsUtils/gsPointGrid.h>

// Utilities
#include <gsCore/gsMemory.h>
#include <gsUtils/gsNorms.h>
#include <gsUtils/gsStopwatch.h>
#include <gsCore/gsFieldCreator.h>
#include <gsUtils/gsInterpolate.h>
#include <gsUtils/gsCollocationMatrix.h>

