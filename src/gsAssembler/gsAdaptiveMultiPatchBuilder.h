/** @file gsAdaptiveMultiPatchBuilder.h

    @brief Provides generic routines for adaptive refinement.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): M. BAHARI
*/

#ifndef GS_ADAPTIVE_MULTIPATCH_BUILDER_H
#define GS_ADAPTIVE_MULTIPATCH_BUILDER_H

#include <fstream>  // For file operations

using namespace gismo;

class GISMO_EXPORT gsAdaptiveMultiPatchBuilder
{
public:

    /** @brief gsAdaptiveMultiPatchBuilder: Main constructor of the r-refinement class
    * @param mapping initial geometry mapping
    * @param numRefine number of uniform refinement steps to perform on the basis before solving
    * @param maxIter maximum number of iterations for the Picard loop
    * @param IntensityMAE intensity of the density function for the Monge-Ampere problem
    * @param numReduce number of degree reduction steps to perform on the basis before solving
    * @param exactGeo if true, boundaries are not adapted (only interfaces) and the MMPDE solver is used even for a single patch; collocation then restores the exact boundary control points
    * @param extraRefCmp number of extra refinement steps for the composition mapping
    */
    // Constructor for one patch compoosition mapping
    gsAdaptiveMultiPatchBuilder(const gsMultiPatch<> mapping,
                                index_t numRefine   = 0,
                                index_t maxIter     = 30,
                                double IntensityMAE = 9.0,
                                index_t numReduce   = 0,
                                index_t numElevate  = 0,
                                bool exactGeo       = false);
    
    // ... optimal Monge-Ampere (or moving mesh) mapping in square to itself
    mutable gsMultiPatch<> MAmapping;

    // ... containes error values over elements, which considered piecewise density function vector
    mutable gismo::gsMatrix<> errorVector;    

    // m_maxIter: max iterations, in moving mesh we want to change max iteration since we start with adaptive mapping
    index_t m_maxIter;

    // degees of freedom used in the computation
    int DoFs;
public:
    // Method to perform uniform refinement on the basis
    void uniformRefine(const index_t numRefine = 1);

    // Method to build a density function from analytic form: we project first f o F into a spline space (avoid composing three functions)
    gsMultiPatch<> buildAnalyticDensity(const gsFunctionExpr<> &f) const;

    // Build and return a density as a MultiPatch object from marked elements using local h-refinement strategies
    gsMultiPatch<> buildDensity(const gsMultiBasis<> Hbasis, const  std::vector<bool> elMarked, const index_t setRhogrid = 0, const  index_t setRhoZero = 0) const;

    //-----------------------------------------
    //  functions to build mapping from density
    //-----------------------------------------
    // Method to build a multipatch Monge-Ampere mapping: tolMAE is tolerance in Picard iterations
    void buildMultiPatch(const gsMultiPatch<> &density, const double tolMAE = 1e-8) const;

    //---------------------------------------------------------------------------------
    //  functions to project the composition of Initial mapping and moving mesh mapping
    //---------------------------------------------------------------------------------
    // Method to build a multipatch adaptive mapping by projection the composition of geometry maps :: L2-projection
    gsMultiPatch<> buildCompMultiPatch(const int quadValue = 1, const bool& sepBoundary = false) const;

    // Method to build a multipatch adaptive mapping by projection the composition of geometry maps :: fitting (penalized least sqaure)
    gsMultiPatch<> buildFitCompMultiPatch(const int numElData = 50, const real_t lambda = 0, const bool& sepboundary = false) const;

    // computes the projection of a composition and return a MultiPatch object :: Collocation
    gsMultiPatch<> buildColCompMultiPatch() const;
    
    //----------------------------------------
    // Useful functions for moving mesh
    //----------------------------------------
    // ... assemble the mass matrix for a given basis in one dimension
    gsSparseMatrix<> assembleMass(const gsBasis<>& basis) const;

    // ... assemble the stiffness matrix for a given basis in one dimension
    gsSparseMatrix<> assembleStiffness(const gsBasis<>& basis) const;
    
    // ... extract boundary condition for each direction 
    gsBoundaryConditions<> boundaryConditionsForDirection( const gsBoundaryConditions<>& bc, index_t direction ) const;

    // ... apply dirichlet by elimination
    void eliminateDirichlet1D(const gsBoundaryConditions<>& bc, const gsOptionList& opt, gsSparseMatrix<> & result) const;

    // ... correct the boundary constrol points only in two dimensions
    void CorrectBoundary(gsMultiPatch<>& Psi, const index_t& patchNumber, const index_t& patch_cmp, const gsMatrix<>& xsoly0, const gsMatrix<>& xsoly1, const gsMatrix<>& x0soly, const gsMatrix<>& x1soly, const bool& corners = false) const;

    // Project control points following normal direction at the boundaries for square domain for moving mesh mapping
    void setBoundaryControlPointsAlongNormal(gsMultiPatch<>& Psi) const;

    // Method to build a inverse multipatch adaptive mapping by projection the composition of geometry maps : fitting
    //gsMultiPatch<> buildInverseMultiPatch(const gsMultiPatch<> lastMAEmapping, const int numElData = 50, const real_t lambda = 0., const bool UpdateInTime = true) const;

    // Method to find the span of a knot vector
    index_t find_span(const gsKnotVector<double>& knots, const index_t& degree, const double& x) const;

    // Method to compute the basis functions and their derivatives
    void basis_functions(const gsKnotVector<double>& knots, const index_t& degree, const double& x, index_t& span,
                     gsVector<double>& d0,
                     gsVector<double>& d1) const;

    // Method to compute the right-hand side vector for composition of mpLeft and MAE mapping assembly in 1D
    void assemble_rhsvector_1d(const gsBasis<>& basis, const index_t& p1, const gsKnotVector<double>& knots_1, const double& valpt, const index_t& dir_valpt,  const index_t& comp_nb,
                           gsVector<double>& rhs) const;
    // Method to compute the right-hand side vector for adaptive multi-patch assembly in 2D
    void assemble_rhsvector_2d(const index_t& p1, const index_t& p2,
                           const gsKnotVector<double>& knots_1, const gsKnotVector<double>& knots_2,
                           const gsMatrix<double>& vector_mp, const gsMatrix<double>& vector_un,
                           gsMatrix<double>& rhs) const;
    // Method to compute the right-hand side vector for adaptive multi-patch assembly in 3D
    void assemble_rhsvector_3d(const index_t& p1, const index_t& p2, const index_t& p3,
                           const gsKnotVector<double>& knots_1, const gsKnotVector<double>& knots_2, const gsKnotVector<double>& knots_3,
                           const gsMatrix<double>& vector_mp, const gsMatrix<double>& vector_un,
                           gsMatrix<double>& rhs) const;

private:
    // Fast diagonalization solver with Dirichlet conditions on the shared basis
    gsPatchPreconditionersCreator<double>::Poisson_FastDiag dirichletPoissonSolver() const;
    // Multipatch version of buildMultiPatch (one mapping of the unit square per patch)
    void buildMultiPatchMMPDE(const gsMultiPatch<> &density, const double tolMAE) const;

    // 1D problem (w(phi) phi')' = 0 along one side of the unit square; returns the tangential coefficients of the edge mapping
    gsMatrix<> solveEdgeMapping(const gsFunction<>& rho, const index_t side) const;

    gsMultiBasis<double> m_basis;
    //... Identity mapping in square to itself
    gsMultiPatch<double> identity_mp; 
    gsMultiPatch<double> initial_mapping;
    double m_IntensityMAE;
    bool m_exactGeo;
    gsFunctionExpr<double> neumann_id;
    gsBoundaryConditions<double> bc_mae;
public:
    // Public members for the mapping basis and Poisson solvers
    gsMultiBasis<double> mapping_basis;
    gsPatchPreconditionersCreator<double>::Poisson_FastDiag Poisson;
    // Dirichlet counterpart of Poisson (used by the multipatch solver)
    gsPatchPreconditionersCreator<double>::Poisson_FastDiag PoissonDir;
};

#endif // GS_ADAPTIVE_MULTIPATCH_BUILDER_H