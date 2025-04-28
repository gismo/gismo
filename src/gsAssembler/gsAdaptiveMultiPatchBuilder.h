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
                                const index_t numElevate,
                                index_t maxIter     = 30,
                                double IntensityMAE = 9.0);

    // Method to project normal control points
    void ProjectionNormalCPoints(gsMultiPatch<>& Psi, int boxMaxNumber = 1) const;

    // Method to build a density function
    gsMultiPatch<> buildAnalyticDensity(const gsFunctionExpr<> &f) const;

    // Method to build a density function from a given solution vector
    gsMultiPatch<> buildDensity(const std::vector<double> &elwiseERROR, const double eps = 0.1, index_t circleN = 0) const;

    // Method to build a multipatch adaptive mapping
    gsMultiPatch<> buildMultiPatch(const gsMultiPatch<> &density, bool composition=true) const;

private:
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
