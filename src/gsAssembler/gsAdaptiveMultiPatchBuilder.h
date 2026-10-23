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
                                index_t maxIter     = 100,
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
    
    // m_mae: tolerance for the Picard loop in the moving mesh solver
    double tolMAE = 1e-8;

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
    // Reuse MAmapping when available; perform at most m_maxIter Picard updates.
    void buildMultiPatch(const gsMultiPatch<> &density) const;

    //---------------------------------------------------------------------------------
    //  functions to project the composition of Initial mapping and moving mesh mapping
    //---------------------------------------------------------------------------------
    // computes the projection of a composition and return a MultiPatch object :: Collocation
    gsMultiPatch<> buildColCompMultiPatch() const;

    // Minimize squared value error + gradientWeight * squared reference-gradient error.
    // Keep collocated boundary traces and enforce C0 continuity at matching interfaces.
    gsMultiPatch<> buildH1CompMultiPatch(const real_t gradientWeight = 1., const real_t quadValue = 2.) const;
    
    //----------------------------------------
    // Useful functions for moving mesh
    //----------------------------------------
    // ... assemble the mass matrix for a given basis in one dimension
    gsSparseMatrix<> assembleMass(const gsBasis<>& basis) const;

    // ... assemble the stiffness matrix for a given basis in one dimension
    gsSparseMatrix<> assembleStiffness(const gsBasis<>& basis) const;

    // ... apply dirichlet by elimination
    void eliminateDirichlet1D(const gsBoundaryConditions<>& bc, const gsOptionList& opt, gsSparseMatrix<> & result) const;

    // Project control points following normal direction at the boundaries for square domain for moving mesh mapping
    void setBoundaryControlPointsAlongNormal(gsMultiPatch<>& Psi) const;

    
    // Transfer a scalar density using the shared mapping/density basis.
    void assemble_rhsvector_2d(const gsBasis<>& basis,
                           const gsMatrix<double>& vector_mp, const gsMatrix<double>& vector_un,
                           gsMatrix<double>& rhs) const;
    void assemble_rhsvector_3d(const gsBasis<>& basis,
                           const gsMatrix<double>& vector_mp, const gsMatrix<double>& vector_un,
                           gsMatrix<double>& rhs) const;

private:
    // Fast diagonalization solver with Dirichlet conditions on the shared basis
    gsPatchPreconditionersCreator<double>::Poisson_FastDiag dirichletPoissonSolver() const;
    // Multipatch version of buildMultiPatch (one mapping of the unit square per patch)
    void buildMultiPatchMMPDE(const gsMultiPatch<> &density) const;
    // Hierarchical 3D sweep: shared edges, shared faces, then volume interiors.
    void buildMultiPatchMMPDE3D(const gsMultiPatch<> &density) const;

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
