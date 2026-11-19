/** @file gsPatchSplitter.h
    @brief Automated Quadrilateral Domain Decomposition & Polychord Splitting Engine for G+Smo.

    Converts either:
    - Mode A: An existing gsMultiPatch<T> (with degenerate ~180 deg corners or elongated aspect ratios)
    - Mode B: A planar B-Rep domain gsPlanarDomain<T> (outer loop + optional interior hole loops)

    into a strictly conforming, analysis-suitable quadrilateral gsMultiPatch<T> suitable
    for C^1/C^2 Argyris spline discretizations and AS-G1 reparametrization.

    Key Properties:
    - Strict 1-to-1 interface conformity (zero hanging nodes / T-junctions).
    - Exact CAD boundary preservation (outer and hole curves split exactly via knot insertion).
    - Bounded corner tangent angles (rejecting collinear ~180 deg and acute sliver corners).
    - 4-RoSy Harmonic Cross-Field + Separatrix Simplification & Structural/MAT decomposition.
    - Global Polychord integer optimization for near-square aspect ratios with minimal patches.
    - Winslow barrier + gsSpringPatch / gsBarrierPatch bijective interior parameterization.
    - Direct ArgyrisSplines formatting without reparametrization (degree elevation, minimum knot spans, orientation unification).

    This file is part of the G+Smo library.
    Author(s): F. Hasanova, Antigravity
*/

#pragma once

#include <gismo.h>
#include <gsIO/gsOptionList.h>
#include <gsModeling/gsPlanarDomain.h>
#include <gsCore/gsMultiPatch.h>
#include <gsModeling/gsAsG1GeometryGenerator.h>

#include <vector>
#include <string>
#include <memory>

namespace gismo {

template <class T = real_t>
class GISMO_EXPORT gsPatchSplitter {
public:
    enum class Strategy {
        AutomaticBest,         // Runs Cross-Field, Structural/MAT, and Templates; picks best valid layout
        HarmonicCrossField,    // 4-RoSy harmonic cross-field + separatrix simplification
        MedialAxisSkeleton,    // Concave-cut, structural & hole-cell decomposition
        TemplateDecomposition, // Topological preset templates & conforming triangulation fallback
        PolychordOnly          // Polychord aspect-ratio splitting only (for existing gsMultiPatch)
    };

    struct SplittingCandidate {
        gsMultiPatch<T> multiPatch;
        int numPatches = 0;
        T minJacobianDet = 0;
        T distortionEnergy = 1e30;
        T minCornerAngleDeg = 0;
        T maxCornerAngleDeg = 180;
        T maxAspectRatio = 1;
        bool isValid = false;
        std::string methodUsed;
    };

    struct QualityReport {
        int nPatches = 0;
        int nInterfaces = 0;
        int nBoundarySides = 0;
        int nExtraordinaryVertices = 0;
        std::vector<int> intValences;
        T minDetJ = 0;
        T maxDetJ = 0;
        T minAspectRatio = 1;
        T maxAspectRatio = 1;
        T minAngleDeg = 0;
        T maxAngleDeg = 0;
        bool hasHangingNodes = false;
        bool isConforming = false;
        bool isValid = false;
    };

public:
    gsPatchSplitter();

    /// Mode B: Construct from a planar domain (outer boundary loop + optional hole loops)
    explicit gsPatchSplitter(const gsPlanarDomain<T>& domain);

    /// Mode A: Construct from an existing multi-patch (re-blocking / corner & polychord splitting)
    explicit gsPatchSplitter(const gsMultiPatch<T>& inputMp);

    // Ingestion
    void setDomain(const gsPlanarDomain<T>& domain);
    void setMultiPatch(const gsMultiPatch<T>& inputMp);
    void loadDomainFromXml(const std::string& xml_file);

    // Options & Configuration
    gsOptionList& options() { return m_options; }
    const gsOptionList& options() const { return m_options; }
    static gsOptionList defaultOptions();

    void setStrategy(Strategy strategy);
    void setDegree(int degree = 3);
    void setWeights(T weightPatches = 2.0, T weightDistortion = 0.01);

    // Execution
    /// Full 4-stage pipeline: coarse conforming decomposition -> global polychord optimization
    /// -> relaxation/parameterization -> ArgyrisSplines formatting (degree & min knot spans).
    gsMultiPatch<T> compute();

    /// Legacy coarse-layout pipeline (preserves coarse 1-patch convex quads when called directly).
    bool computeSplitting();

    // Diagnostics & Validation
    bool isConforming() const;
    bool checkJacobians(T* min_det = nullptr) const;
    bool checkPatchCornersConvex(T* min_tri_area2 = nullptr) const;
    bool checkNoSelfIntersections() const;
    bool verifyBoundaryMatch(T tol = 1e-4, T* maxDistOut = nullptr) const;
    static bool verifyBoundaryMatch(const gsPlanarDomain<T>& dom, const gsMultiPatch<T>& mp, T tol = 1e-4, T* maxDistOut = nullptr);
    static bool validateMultiPatch(const gsPlanarDomain<T>& dom, const gsMultiPatch<T>& mp, T bdrTol = 1e-5, std::string* reason = nullptr,
                                   T minAngleDeg = 12.0, T maxAngleDeg = 155.0);
    std::vector<int> getVertexValences() const;

    QualityReport analyzeQuality(const gsMultiPatch<T>& result) const;
    void printQualityReport(const gsMultiPatch<T>& result) const;

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
    SplittingCandidate computeMultiHoleStructuralTemplates();
    SplittingCandidate computeHoleFreeTriangulationSplitting();
    SplittingCandidate computeMultiHoleTriangulationSplitting();

    // Mode A Helpers
    SplittingCandidate computeModeACoarseLayout();
    static bool extractPlanarDomainFromMultiPatch(const gsMultiPatch<T>& mp, gsPlanarDomain<T>& outDom);
    bool splitDegenerateSinglePatchOGrid(const gsMultiPatch<T>& inMp, gsMultiPatch<T>& outMp) const;

    // Stage 2 & 3 Helpers
    bool optimizePolychordSplits(gsMultiPatch<T>& mp) const;
    bool reducePolychordStrips(gsMultiPatch<T>& mp);
    bool harmonizeMultiPatchKnots(gsMultiPatch<T>& mp) const;
    bool relaxAndParameterize(gsMultiPatch<T>& mp, bool runBarrier = true);
    bool formatForArgyris(gsMultiPatch<T>& mp) const;

    // Optimization & Post-Processing
    bool optimizeWinslowBarrier(gsMultiPatch<T>& mp, T minDetThreshold = 1e-7) const;
    void ensureExactBoundaryConformity(gsMultiPatch<T>& mp);
    T evaluateDistortionEnergy(const gsMultiPatch<T>& mp, T* minDetOut = nullptr,
                               T* minAngleOut = nullptr, T* maxAngleOut = nullptr) const;
    T getDomainMinCornerAngleDeg() const;

private:
    gsPlanarDomain<T>   m_domain;
    gsMultiPatch<T>     m_inputMp;
    bool                m_hasInputMp = false;
    bool                m_hasDomain = false;
    mutable T           m_cachedDomainMinCornerAngle = -1;

    gsOptionList        m_options;
    Strategy            m_strategy = Strategy::AutomaticBest;
    int                 m_degree = 3;
    T                   m_weightPatches = 2.0;
    T                   m_weightDistortion = 0.01;

    SplittingCandidate  m_bestCandidate;
};

} // namespace gismo

#ifndef GISMO_BUILD_LIB
#include <gsModeling/gsPatchSplitter.hpp>
#endif
