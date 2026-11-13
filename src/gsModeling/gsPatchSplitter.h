/** @file gsPatchSplitter.h
    @brief Automated Quadrilateral Domain Decomposition Engine for Planar Domains in G+Smo.

    Decomposes arbitrary planar domains with outer boundaries and interior holes into
    a conforming, non-degenerate multi-patch quadrilateral B-spline geometry (gsMultiPatch)
    certified with strictly positive Jacobian determinant (det(J) > 0).

    Features:
    - Strict invariance of outer boundaries and holes (exact CAD B-spline geometry preserved).
    - Conforming topological quad ring (annulus) around every interior hole.
    - Dual-method decomposition: Boundary-Aligned Harmonic Cross-Field and Medial Axis Transform (MAT).
    - Topological preset template decomposition fallback.
    - Pareto-optimal automated winner selection minimizing patch count and distortion energy.
    - Winslow Barrier functional optimization for parameterization injectivity and regularity.
    - Direct downstream compatibility with gsAsG1GeometryGenerator and gsASG1LocalFitter.

    This file is part of the G+Smo library.
    Author(s): F. Hasanova, Antigravity
*/

#pragma once

#include <gismo.h>
#include <gsModeling/gsPlanarDomain.h>
#include <gsCore/gsMultiPatch.h>
#include <gsModeling/gsAsG1GeometryGenerator.h>

#include <vector>
#include <string>
#include <memory>

namespace gismo {

template <class T = real_t>
class gsPatchSplitter {
public:
    enum class Strategy {
        AutomaticBest,         // Runs both Cross-Field and MAT; picks minimal patches & distortion
        HarmonicCrossField,    // Cross-field / frame-field separatrix tracing only
        MedialAxisSkeleton,    // Voronoi MAT chordal projection only
        TemplateDecomposition  // Topological preset templates
    };

    struct SplittingCandidate {
        gsMultiPatch<T> multiPatch;
        int numPatches = 0;
        T minJacobianDet = 0;
        T distortionEnergy = 1e30;
        bool isValid = false;
        std::string methodUsed;
    };

public:
    gsPatchSplitter() = default;
    explicit gsPatchSplitter(const gsPlanarDomain<T>& domain);

    // Ingestion
    void setDomain(const gsPlanarDomain<T>& domain);
    void loadDomainFromXml(const std::string& xml_file);

    // Configuration
    void setStrategy(Strategy strategy);
    void setDegree(int degree = 3);
    void setWeights(T weightPatches = 1.0, T weightDistortion = 0.05);

    // Execution
    bool computeSplitting();

    // Diagnostics & Validation
    bool isConforming() const;
    bool checkJacobians(T* min_det = nullptr) const;
    bool checkPatchCornersConvex(T* min_tri_area2 = nullptr) const;
    bool checkNoSelfIntersections() const;
    bool verifyBoundaryMatch(T tol = 1e-4, T* maxDistOut = nullptr) const;
    static bool verifyBoundaryMatch(const gsPlanarDomain<T>& dom, const gsMultiPatch<T>& mp, T tol = 1e-4, T* maxDistOut = nullptr);
    static bool validateMultiPatch(const gsPlanarDomain<T>& dom, const gsMultiPatch<T>& mp, T bdrTol = 1e-5, std::string* reason = nullptr);
    std::vector<int> getVertexValences() const;

    // Export & Downstream Interface
    const gsMultiPatch<T>& getMultiPatch() const;
    gsMultiPatch<T> getBilinearSkeleton() const;
    void configureAsG1Generator(gsAsG1GeometryGenerator<T>& generator) const;

    void exportParaview(const std::string& filename_prefix, unsigned npts = 900) const;

    // Additional Accessors
    const SplittingCandidate& bestCandidate() const { return m_bestCandidate; }
    Strategy strategy() const { return m_strategy; }
    int degree() const { return m_degree; }
    const gsPlanarDomain<T>& domain() const { return m_domain; }

private:
    // Splitting Core
    SplittingCandidate computeCrossFieldSplitting();
    SplittingCandidate computeMedialAxisSplitting();
    SplittingCandidate computeTemplateSplitting();
    SplittingCandidate computeHoleFreeStructuralSplitting();
    SplittingCandidate computeHoleFreeTriangulationSplitting();

    // Optimization & Post-Processing
    bool optimizeWinslowBarrier(gsMultiPatch<T>& mp, T minDetThreshold = 1e-7);
    void ensureExactBoundaryConformity(gsMultiPatch<T>& mp);
    T evaluateDistortionEnergy(const gsMultiPatch<T>& mp, T* minDetOut = nullptr) const;

private:
    gsPlanarDomain<T>   m_domain;
    Strategy            m_strategy = Strategy::AutomaticBest;
    int                 m_degree = 3;
    T                   m_weightPatches = 1.0;
    T                   m_weightDistortion = 0.05;

    SplittingCandidate  m_bestCandidate;
};

} // namespace gismo

#ifndef GISMO_BUILD_LIB
#include <gsModeling/gsPatchSplitter.hpp>
#endif

