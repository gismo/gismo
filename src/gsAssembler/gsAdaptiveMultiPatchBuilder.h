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
    // Constructor
    gsAdaptiveMultiPatchBuilder(const gsMultiBasis<double> basis,
                                const gsMultiPatch<> mapping,
                                index_t maxIter     = 30,
                                double IntensityMAE = 9.0,
                                index_t numReduce = 0);

    // Method to project normal control points
    void ProjectionNormalCPoints(gsMultiPatch<>& Psi, int boxMaxNumber = 1) const;

    // Method to build a density function
    gsMultiPatch<> buildAnalyticDensity(const gsFunctionExpr<> &f) const;

    // Method to build a density function from a given solution vector
    gsMultiPatch<> buildDensity(const std::vector<double> &elwiseERROR, const double eps = 0.1, index_t circleN = 0, bool maxminVar = true) const;

    // Build and return a density as a MultiPatch object from solution vector using local h-refinement strategies
    gsMultiPatch<> buildStrategyDensity(const std::vector<double> &elwiseERROR, const double MarkPercentage = 0.9) const;

    // Method to build a multipatch adaptive mapping
    gsMultiPatch<> buildMultiPatch(const gsMultiPatch<> &density) const;

    // Method to build a multipatch adaptive mapping TODO: CAN BE OPTIMIZED @BAHARI
    gsMultiPatch<> buildMovingMultiPatch(const gsMultiPatch<> &density, gsMultiPatch<> lsPsi, int Niter = 0) const;

    // Method to build a multipatch adaptive mapping by projection the composition of geometry maps
    gsMultiPatch<> buildCompMultiPatch(gsMultiPatch<> Psitp, double quadValue = 1.) const;

    // Method to find the span of a knot vector
    index_t find_span(const gsKnotVector<double>& knots, const index_t& degree, const double& x) const;

    // Method to compute the basis functions and their derivatives
    void basis_functions(const gsKnotVector<double>& knots, const index_t& degree, const double& x, index_t& span,
                     gsVector<double>& d0,
                     gsVector<double>& d1) const;

    // Method to compute the right-hand side vector for the adaptive multi-patch assembly
    void assemble_rhsvector_ad(const index_t& p1, const index_t& p2,
                           const gsKnotVector<double>& knots_1, const gsKnotVector<double>& knots_2,
                           const gsMatrix<double>& vector_mp, const gsMatrix<double>& vector_cp,
                           const gsMatrix<double>& vector_un, gsMatrix<double>& rhs) const;

private:
    gsMultiBasis<double> n_basis;
    gsMultiBasis<double> m_basis;
    gsMultiPatch<double> m_mapping;
    gsMultiPatch<double> mp;
    index_t m_maxIter;
    double m_IntensityMAE;
    gsBoundaryConditions<> bc_mae;
public:
    gsPatchPreconditionersCreator<double>::Poisson_FastDiag Poisson;
};

#endif // GS_ADAPTIVE_MULTIPATCH_BUILDER_H
