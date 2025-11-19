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
    */
    // Constructor
    gsAdaptiveMultiPatchBuilder(const gsMultiPatch<> mapping,
                                index_t numRefine   = 0,
                                index_t maxIter     = 30,
                                double IntensityMAE = 9.0,
                                index_t numReduce   = 0);
    //... Square mapping
    gsMultiPatch<> mp; 
    // ... optimal Monge-Ampere mapping in unit-square to itself
    mutable gsMultiPatch<> MAmapping;
    // ... density coefs
    mutable gismo::gsMatrix<> errorVector;    
    // m_maxIter: max iterations, in moving mesh we want to change max iteration since we start with adaptive mapping
    index_t m_maxIter;
    // degees of freedom used in the computation
    int DoFs;
public:

    //... uniform refinement
    void uniformRefine(const index_t numRefine = 1);

    // Project control points following  normal direction at the boundaries for square domain 
    void NormalProjectPts(gsMultiPatch<>& Psi) const;

    // Method to build a density function from analytic form
    gsMultiPatch<> buildAnalyticDensity(const gsFunctionExpr<> &f) const;

    // Build and return a density as a MultiPatch object from marked elements using local h-refinement strategies
    gsMultiPatch<> buildDensity(const gsMultiBasis<> Givbasis, const  std::vector<bool> elMarked, const  index_t setRhoZero = 0) const;
    
    //-----------------------------------------
    //  functions to build mapping from density
    //-----------------------------------------
    // Method to build a multipatch Monge-Ampere mapping
    void buildMultiPatch(const gsMultiPatch<> &density, const double tolMAE = 1e-5) const;

    // Method to build a multipatch adaptive mapping by projection the composition of geometry maps : L2-projection
    gsMultiPatch<> buildCompMultiPatch(const gsMultiBasis<> Cbasis, const int quadValue = 1) const;

    // Method to build a multipatch adaptive mapping by projection the composition of geometry maps : fitting
    gsMultiPatch<> buildFitCompMultiPatch(const gsMultiBasis<> Cbasis, const int numElData = 10, const real_t lambda = 0.) const;

    // computes the projection of a composition and return a MultiPatch object :: Collocation
    gsMultiPatch<> buildColCompMultiPatch(const gsMultiBasis<> Cbasis) const;
    
    //----------------------------------------
    // Useful functions for time moving meshes
    //----------------------------------------
    // Method to find the span of a knot vector
    index_t find_span(const gsKnotVector<double>& knots, const index_t& degree, const double& x) const;

    // Method to compute the basis functions and their derivatives
    void basis_functions(const gsKnotVector<double>& knots, const index_t& degree, const double& x, index_t& span,
                     gsVector<double>& d0,
                     gsVector<double>& d1) const;
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
    gsMultiBasis<double> m_basis;
    gsMultiPatch<double> m_mapping;
    double m_IntensityMAE;
    gsBoundaryConditions<> bc_mae;
public:
    gsPatchPreconditionersCreator<double>::Poisson_FastDiag Poisson;
};

#endif // GS_ADAPTIVE_MULTIPATCH_BUILDER_H