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
    * @param mapping initial geometry mapping, with arbitrary finite parameter boxes
    * @param numRefine number of uniform refinement steps to perform on the basis before solving
    * @param maxIter maximum number of iterations for the Picard loop
    * @param IntensityMAE intensity of the density function for the Monge-Ampere problem
    * @param numReduce number of degree reduction steps to perform on the basis before solving
    * @param numElevate number of degree elevation steps for the moving mesh basis, used when numReduce is zero
    * @param exactGeo if true, boundaries are not adapted (only interfaces) and the MMPDE solver is used even for a single patch; collocation then restores the exact boundary control points
    */
    // Constructor for one patch compoosition mapping
    gsAdaptiveMultiPatchBuilder(const gsMultiPatch<> mapping,
                                index_t numRefine   = 0,
                                index_t maxIter     = 100,
                                real_t IntensityMAE = 9.0,
                                index_t numReduce   = 0,
                                index_t numElevate  = 0,
                                bool exactGeo       = false);

    // Optimal Monge-Ampere (or moving mesh) mapping [0,1]^d -> [0,1]^d.
    mutable gsMultiPatch<> MAmapping;

    // A_p : [0,1]^d -> the original parameter box of patch p.
    const gsMultiPatch<> & parameterConversion() const { return m_parameterConversion; }

    // F_p o A_p, represented exactly by reparameterizing a copy of the input.
    // Use this geometry with MAmapping for direct composition and evaluation.
    const gsMultiPatch<> & normalizedGeometry() const { return initial_mapping; }

    // ... containes error values over elements, which considered piecewise density function vector
    mutable gismo::gsMatrix<> errorVector;

    // m_maxIter: max iterations, in moving mesh we want to change max iteration since we start with adaptive mapping
    index_t m_maxIter;

    // m_mae: tolerance for the Picard loop in the moving mesh solver
    real_t tolMAE = 1e-8;

    // degees of freedom used in the computation
    int DoFs;
public:
    // Method to perform uniform refinement on the basis
    void uniformRefine(const index_t numRefine = 1);

    // Project f o F_p o A_p into the unit-domain moving mesh spline space.
    gsMultiPatch<> buildAnalyticDensity(const gsFunctionExpr<> &f) const;

    // Build a density from marked elements of Hbasis, defined on unit parameter boxes.
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
                           const gsMatrix<>& vector_mp, const gsMatrix<>& vector_un,
                           gsMatrix<>& rhs) const;
    void assemble_rhsvector_3d(const gsBasis<>& basis,
                           const gsMatrix<>& vector_mp, const gsMatrix<>& vector_un,
                           gsMatrix<>& rhs) const;

private:
    // Fast diagonalization solver with Dirichlet conditions on the shared basis
    gsPatchPreconditionersCreator<>::Poisson_FastDiag dirichletPoissonSolver() const;
    // Multipatch version of buildMultiPatch (one mapping of the unit square per patch)
    void buildMultiPatchMMPDE(const gsMultiPatch<> &density) const;
    // Hierarchical 3D sweep: shared edges, shared faces, then volume interiors.
    void buildMultiPatchMMPDE3D(const gsMultiPatch<> &density) const;

    gsMultiBasis<> m_basis;
    //... Identity mapping in square to itself
    gsMultiPatch<> identity_mp;
    gsMultiPatch<> m_parameterConversion;
    // The normalized input geometry F_p o A_p; coefficients and weights are preserved.
    gsMultiPatch<> initial_mapping;
    real_t m_IntensityMAE;
    bool m_exactGeo;
    gsFunctionExpr<> neumann_id;
    gsBoundaryConditions<> bc_mae;
public:
    // Composition bases on [0,1]^d and Poisson solvers.
    gsMultiBasis<> mapping_basis;
    gsPatchPreconditionersCreator<>::Poisson_FastDiag Poisson;
    // Dirichlet counterpart of Poisson (used by the multipatch solver)
    gsPatchPreconditionersCreator<>::Poisson_FastDiag PoissonDir;
};

#endif // GS_ADAPTIVE_MULTIPATCH_BUILDER_H
