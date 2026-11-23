/** @file gsPatchSplitter.hpp
    @brief Implementation of gsPatchSplitter<T>: Planar Quadrilateral Domain Decomposition.

    Decomposes 2D planar domains (gsPlanarDomain) into conforming multi-patch
    tensor B-spline geometries (gsMultiPatch) suitable for isogeometric analysis
    and downstream AS-G1 reparametrization (gsASG1LocalFitter / gsAsG1GeometryGenerator).

    Key Properties & Implemented Pipeline:
    1. Exact CAD boundary preservation: outer boundary curves (and hole curves)
       are retained verbatim or subdivided via exact B-spline splitting (splitAt)
       and knot harmonization without perturbing boundary control points.
    2. Hole-free decomposition: evaluates structural quad layouts (convex 1-patch
       preservation, O-grid, ribbon, reflex L-cut, star fans/kites, comb spine-teeth)
       and triangulation-based 1-to-3 Catmull-Clark quadification fallback, selecting
       the lowest-cost valid candidate by patch count and distortion energy.
    3. Domains with holes: quad-ring and block-ring decomposition templates for
       single-hole and supported multi-hole topologies.
    4. Parameterization & validation: Coons blending with interior-only, barrier-guarded
       smoothing; checks conformity, exact boundary match, per-patch nonzero
       single-sign Jacobian determinant, strictly convex patch corners, and absence
       of self-intersections or patch overlaps.

    This file is part of the G+Smo library.
    Author(s): F. Hasanova, Antigravity
*/

#pragma once

#include <gsModeling/gsPatchSplitter.h>
#include <gsModeling/gsPlanarDomainGenerator.h>
#include <gsModeling/gsPatchSplitterGeom.hpp>
#include <gsModeling/gsPatchSplitterTriMesh.hpp>
#include <gsModeling/gsPatchSplitterTemplates.hpp>
#include <gsModeling/gsPatchSplitterMAT.hpp>
#include <gsModeling/gsPatchSplitterPolychord.hpp>
#include <gsModeling/gsPatchSplitterCrossField.hpp>
#include <gsModeling/gsPatchSplitterMotorcycle.hpp>
#include <gsModeling/gsPatchSplitterMerged.hpp>
#include <gsModeling/gsASG1LocalFitter.h>
#include <gsModeling/gsCoonsPatch.h>
#include <gsNurbs/gsBSpline.h>
#include <gsNurbs/gsTensorBSpline.h>
#include <gsNurbs/gsKnotVector.h>
#include <gsIO/gsReadFile.h>
#include <gsIO/gsWriteParaview.h>
#include <gsMatrix/gsSparseMatrix.h>
#include <gsMatrix/gsSparseSolver.h>

#include <vector>
#include <cmath>
#include <algorithm>
#include <map>
#include <set>
#include <iostream>
#include <iomanip>
#include <limits>
#include <memory>

namespace gismo {

// ============================================================================
// Constructor & Ingestion
// ============================================================================

template <class T>
gsOptionList gsPatchSplitter<T>::defaultOptions()
{
    gsOptionList opt;
    opt.addInt("Strategy", "0: Auto/Hybrid, 1: CrossField, 2: ConcaveCut/MAT, 3: PolychordOnly", 0);
    opt.addReal("TargetAspectRatio", "Target patch aspect ratio Lu/Lv", 1.0);
    opt.addReal("MaxAspectRatio", "Preferred maximum patch aspect ratio", 2.0);
    opt.addReal("PatchCountWeight", "Penalty weight per patch in polychord optimization", 0.18);
    opt.addReal("AngleThresholdDeg", "Boundary corner angle above which corner is split", 135.0);
    opt.addReal("MinCornerAngleDeg", "Minimum admissible patch corner tangent angle (deg)", 12.0);
    opt.addReal("MaxCornerAngleDeg", "Maximum admissible patch corner tangent angle (deg)", 155.0);
    opt.addReal("StripCollapseRatio", "Relative distance to snap parallel separatrices", 0.25);
    opt.addReal("SmoothJointTolDeg", "Boundary angle deviation from 180 deg below which CAD curve joints are treated as smooth G1", 22.0);
    opt.addReal("TMeshARWeight", "Aspect-ratio regularization weight in T-mesh integer quantization (multigrid-friendly)", 0.35);
    opt.addInt("Degree", "Target B-spline degree (>= 3 for ArgyrisSplines)", 3);
    opt.addInt("Regularity", "Target spline regularity r (1 <= r <= p - 2)", 1);
    opt.addInt("MinKnotSpans", "Minimum knot spans per patch direction (>= 3 for ArgyrisSplines)", 3);
    opt.addInt("MaxPatches", "Maximum admissible patch count for macro-patch decomposition", 140);
    opt.addSwitch("EnablePolychordSplitting", "Run global polychord aspect-ratio optimization in compute()", true);
    opt.addSwitch("ApplyASG1Fitting", "Run gsASG1LocalFitter in formatForArgyris()", false);
    opt.addInt("ParameterizationMethod", "0: Coons + Spring + Winslow, 1: Coons + gsBarrierPatch (if enabled)", 0);
    return opt;
}

template <class T>
gsPatchSplitter<T>::gsPatchSplitter()
    : m_options(defaultOptions())
{
}

template <class T>
gsPatchSplitter<T>::gsPatchSplitter(const gsPlanarDomain<T>& domain)
    : m_domain(domain),
      m_hasDomain(true),
      m_options(defaultOptions())
{
}

template <class T>
gsPatchSplitter<T>::gsPatchSplitter(const gsMultiPatch<T>& inputMp)
    : m_inputMp(inputMp),
      m_hasInputMp(true),
      m_options(defaultOptions())
{
}

template <class T>
void gsPatchSplitter<T>::setDomain(const gsPlanarDomain<T>& domain)
{
    m_domain = domain;
    m_hasDomain = true;
    m_hasInputMp = false;
    m_didDomainReblock = false;
    m_cachedDomainMinCornerAngle = -1;
    m_bestCandidate = SplittingCandidate();
}

template <class T>
void gsPatchSplitter<T>::setMultiPatch(const gsMultiPatch<T>& inputMp)
{
    m_inputMp = inputMp;
    m_hasInputMp = true;
    m_didDomainReblock = false;
    m_bestCandidate = SplittingCandidate();
}

template <class T>
void gsPatchSplitter<T>::loadDomainFromXml(const std::string& xml_file)
{
    gsReadFile<T> rf(xml_file);
    memory::unique_ptr<gsPlanarDomain<T>> pDom = rf;
    if (pDom) {
        m_domain = *pDom;
        m_hasDomain = true;
        m_hasInputMp = false;
        m_didDomainReblock = false;
        m_cachedDomainMinCornerAngle = -1;
        m_bestCandidate = SplittingCandidate();
    } else {
        GISMO_ERROR("XML file " << xml_file << " does not contain a gsPlanarDomain.");
    }
}

template <class T>
void gsPatchSplitter<T>::setStrategy(Strategy strategy)
{
    m_strategy = strategy;
    m_options.setInt("Strategy", static_cast<int>(strategy));
}

template <class T>
void gsPatchSplitter<T>::setDegree(int degree)
{
    m_degree = std::max(1, degree);
    m_options.setInt("Degree", m_degree);
}

template <class T>
void gsPatchSplitter<T>::setWeights(T weightPatches, T weightDistortion)
{
    m_weightPatches = weightPatches;
    m_weightDistortion = weightDistortion;
}

template <class T>
T gsPatchSplitter<T>::getDomainMinCornerAngleDeg() const
{
    if (m_cachedDomainMinCornerAngle < 0 && m_hasDomain) {
        m_cachedDomainMinCornerAngle = internal_splitter::computeDomainMinCornerAngleDeg(m_domain);
    }
    return (m_cachedDomainMinCornerAngle >= 0) ? m_cachedDomainMinCornerAngle : T(90);
}

// ============================================================================
// Approach 1: Boundary-Aligned 4-RoSy Harmonic Cross-Field Decomposition
// ============================================================================

template <class T>
typename gsPatchSplitter<T>::SplittingCandidate
gsPatchSplitter<T>::computeCrossFieldSplitting()
{
    using namespace internal_splitter;
    SplittingCandidate bestCand;
    bestCand.methodUsed = "HarmonicCrossField";
    bestCand.isValid = false;

    if (m_domain.numLoops() == 0) return bestCand;

    const gsMatrix<T> bbox = m_domain.boundingBox();
    const T diag = std::hypot(bbox(0, 1) - bbox(0, 0), bbox(1, 1) - bbox(1, 0));
    if (diag <= T(1e-12)) return bestCand;

    const T domMinAng = getDomainMinCornerAngleDeg();
    const T minAngTol = std::min(static_cast<T>(m_options.getReal("MinCornerAngleDeg")), T(0.85) * domMinAng);
    const T effMinAng = (m_domain.numHoles() > 0)
                      ? std::min(minAngTol, std::max(T(4.0), minAngTol - T(4.0)))
                      : std::min(minAngTol, std::max(T(4.0), minAngTol - T(3.5)));
    const T maxAngTol = static_cast<T>(m_options.getReal("MaxCornerAngleDeg"));
    const T effMaxAng = std::max(maxAngTol, T(164.9));
    const int maxPatchesCap = std::max(1, static_cast<int>(m_options.getInt("MaxPatches")));
    const T stripRatio = static_cast<T>(m_options.getReal("StripCollapseRatio"));
    const T angThresh  = static_cast<T>(m_options.getReal("AngleThresholdDeg"));
    const T smoothJointTol = static_cast<T>(m_options.getReal("SmoothJointTolDeg"));
    const T tmeshARWeight  = static_cast<T>(m_options.getReal("TMeshARWeight"));

    // Stage 0: Virtual G1-Consolidated CAD Boundary Loops
    std::vector<VirtualCADLoop<T>> vloops = buildVirtualCADLoops(m_domain, m_degree, smoothJointTol);

    // Stage 1: Solve Curvature-Regularized Robin MBO 4-RoSy Cross-Field ONCE on the background FEM mesh
    const T targetH = diag / ((m_domain.numHoles() >= 4) ? T(26) : T(22));
    TriMeshData<T> fineMesh;
    if (!buildConstrainedDomainTriangulation(m_domain, m_degree, targetH, true, 0, fineMesh))
        return bestCand;
    applyDomainLawsonFlips(fineMesh);
    refineInteriorTriangles(fineMesh, targetH, 2);

    CrossFieldMeshAccel<T> accel;
    accel.build(fineMesh);
    CrossFieldSolution<T> cf = solveRegularized4RoSy(fineMesh, vloops, &accel);

    const bool cfDbg = (std::getenv("GISMO_CF_DEBUG") != nullptr);
    if (cfDbg) {
        int totPinned = 0;
        for (const auto& vl : vloops) totPinned += static_cast<int>(vl.corners.size());
        std::cerr << "[CF_DBG] Stage0/1: loops=" << vloops.size() << " pinnedCorners=" << totPinned
                  << " robinSings=" << cf.singularities.size() << "\n";
    }

    auto countNegDense = [&](const gsMultiPatch<T>& mpCheck) -> int {
        int cnt = 0;
        for (size_t p = 0; p < mpCheck.nPatches(); ++p)
            if (evalPatchMinDetDense(mpCheck.patch(p), 49) <= T(1e-8)) ++cnt;
        return cnt;
    };
    auto repairMicroFolds = [&](gsMultiPatch<T>& mpFix) {
        for (size_t p = 0; p < mpFix.nPatches(); ++p) {
            T d0 = evalPatchMinDetDense(mpFix.patch(p), 49);
            if (d0 <= T(1e-8)) {
                if (auto* tb = dynamic_cast<gsTensorBSpline<2, T>*>(&mpFix.patch(p))) {
                    gsTensorBSpline<2, T> bk(*tb);
                    quickHarmonicPatchRepair(*tb, T(1e-7));
                    T d1 = evalPatchMinDetDense(*tb, 49);
                    if (d1 < d0) { *tb = bk; d1 = d0; }
                    if (d1 <= T(1e-8)) {
                        try {
                            gsMultiPatch<T> bndC;
                            bndC.addPatch(tb->boundary(boundary::south));
                            bndC.addPatch(tb->boundary(boundary::east));
                            bndC.addPatch(tb->boundary(boundary::north));
                            bndC.addPatch(tb->boundary(boundary::west));
                            gsSpringPatch<T> sp(bndC);
                            const gsGeometry<T>& res = sp.compute();
                            if (const auto* sTb = dynamic_cast<const gsTensorBSpline<2, T>*>(&res)) {
                                if (evalPatchMinDetDense(*sTb, 49) > d1)
                                    *tb = *sTb;
                            }
                        } catch (...) {}
                    }
                }
            }
        }
    };
    auto scoreCand = [&](const SplittingCandidate& c) -> T {
        std::vector<std::vector<patchCorner>> evs;
        c.multiPatch.topology().getEVs(evs, false);
        T maxER = evalMultiPatchMaxEdgeRatio(c.multiPatch);
        T ratioPen = (maxER > T(5.0)) ? T(1.5) * (maxER - T(5.0)) : T(0);
        T angPen = ((c.maxCornerAngleDeg > T(146.0)) ? T(0.35) * (c.maxCornerAngleDeg - T(146.0)) : T(0))
                 + ((c.minCornerAngleDeg < T(22.0))  ? T(0.60) * (T(22.0) - c.minCornerAngleDeg) : T(0));
        T detPen = (c.minJacobianDet < T(1e-4)) ? T(6.0) : T(0);
        T globalBonus = (c.numPatches <= 32) ? -T(2.5) : T(0);
        return m_weightPatches * static_cast<T>(c.numPatches)
             + m_weightDistortion * c.distortionEnergy
             + T(0.25) * static_cast<T>(evs.size())
             + ratioPen + angPen + detPen + globalBonus;
    };

    const int compactTarget = (m_domain.numHoles() == 0) ? 6
                            : (m_domain.numHoles() == 1 ? 6 : std::max(12, 5 * m_domain.numHoles() + 4));
    auto evaluateAndRegister = [&](gsMultiPatch<T>& candMp, const std::string& label) -> bool {
        if (candMp.nPatches() == 0 || static_cast<int>(candMp.nPatches()) > maxPatchesCap + 80) {
            if (cfDbg) std::cerr << "[CF_DBG] " << label << " P=" << candMp.nPatches() << " skip: empty or > headroom\n";
            return false;
        }
        // Skip expensive relaxation if we already have a much smaller well-conditioned candidate
        if (bestCand.isValid && static_cast<int>(candMp.nPatches()) > bestCand.numPatches + 12 &&
            bestCand.minJacobianDet > T(2e-4) && bestCand.maxCornerAngleDeg <= T(150.0))
            return false;

        gsASG1LocalFitter<T>::unifyPatchOrientations(candMp);
        candMp.computeTopology(T(1e-4), false, false);
        if (2 * candMp.nInterfaces() + candMp.nBoundary() != 4 * candMp.nPatches()) {
            if (cfDbg) std::cerr << "[CF_DBG] " << label << " P=" << candMp.nPatches() << " FAIL: non-conforming IF=" << candMp.nInterfaces() << " B=" << candMp.nBoundary() << "\n";
            return false;
        }

        auto checkConforming = [](gsMultiPatch<T>& m) -> bool {
            m.computeTopology(T(1e-4), false, false);
            return (2 * m.nInterfaces() + m.nBoundary() == 4 * m.nPatches());
        };

        repairMicroFolds(candMp);

        // 1. Pre-reduction interior-vertex smoothing + bad-corner angle relaxation (on clean low-knot curves)
        if (candMp.nPatches() >= 2 && candMp.nPatches() <= 160) {
            gsMultiPatch<T> preV(candMp);
            int negBefore = countNegDense(candMp);
            relaxInteriorVertices(candMp, 5);
            relaxBadCornerAngles(candMp, effMinAng, effMaxAng, 2);
            repairMicroFolds(candMp);
            if (!checkConforming(candMp) || countNegDense(candMp) > negBefore)
                candMp = give(preV);
        }

        // 2. Single-pass Polychord dual-sheet strip collapse (NO repeated harmonizeMultiPatchKnots blow-up!)
        if (candMp.nPatches() > 2) {
            int negBefore = countNegDense(candMp);
            gsMultiPatch<T> redMp(candMp);
            if (reducePolychordStrips(redMp) && redMp.nPatches() > 0 && redMp.nPatches() < candMp.nPatches()) {
                redMp.computeTopology(T(1e-4), false, false);
                repairMicroFolds(redMp);
                if (checkConforming(redMp) && countNegDense(redMp) <= negBefore)
                    candMp = give(redMp);
            }
        }

        if (static_cast<int>(candMp.nPatches()) > maxPatchesCap) {
            if (cfDbg) std::cerr << "[CF_DBG] " << label << " P=" << candMp.nPatches() << " skip: > cap(" << maxPatchesCap << ") after strip reduction\n";
            return false;
        }

        // 3. Post-reduction interior-vertex smoothing, corner-angle rescue, and Winslow + coupled relaxation
        if (candMp.nPatches() >= 2) {
            gsMultiPatch<T> preV(candMp);
            int negBefore = countNegDense(candMp);
            relaxInteriorVertices(candMp, 5);
            relaxBadCornerAngles(candMp, effMinAng, effMaxAng, 3);
            repairMicroFolds(candMp);
            if (!checkConforming(candMp) || countNegDense(candMp) > negBefore)
                candMp = give(preV);
        }
        {
            gsMultiPatch<T> preWin(candMp);
            int negBefore = countNegDense(candMp);
            optimizeWinslowBarrier(candMp);
            repairMicroFolds(candMp);
            if (!checkConforming(candMp) || countNegDense(candMp) > negBefore)
                candMp = give(preWin);
        }
        if (candMp.nPatches() >= 2 && candMp.nPatches() <= 140) {
            gsMultiPatch<T> preCpl(candMp);
            int negBefore = countNegDense(candMp);
            optimizeCoupledInterfacesAndInteriors(candMp, T(1e-7), 22, true);
            repairMicroFolds(candMp);
            if (!checkConforming(candMp) || countNegDense(candMp) > negBefore)
                candMp = give(preCpl);
        }
        if (candMp.nPatches() >= 2) {
            gsMultiPatch<T> preAng(candMp);
            int negBefore = countNegDense(candMp);
            relaxBadCornerAngles(candMp, effMinAng, effMaxAng, 3);
            repairMicroFolds(candMp);
            if (!checkConforming(candMp) || countNegDense(candMp) > negBefore)
                candMp = give(preAng);
        }

        gsASG1LocalFitter<T>::unifyPatchOrientations(candMp);
        candMp.computeTopology(T(1e-4), false, false);
        repairMicroFolds(candMp);

        // 4. Symmetric 3-span interior-DOF rescue for strongly curved low-patch layouts (e.g. amoeba P=5/24)
        //    where 1-span cubic interiors (only 2x2 interior coefs) lack freedom to untangle curved lobes!
        if (countNegDense(candMp) > 0 && candMp.nPatches() <= 64) {
            gsMultiPatch<T> refMp(candMp);
            const int negPre = countNegDense(refMp);
            for (size_t p = 0; p < refMp.nPatches(); ++p) {
                if (auto* tb = dynamic_cast<gsTensorBSpline<2, T>*>(&refMp.patch(p))) {
                    safeNormalizeKnotVector01(tb->knots(0));
                    safeNormalizeKnotVector01(tb->knots(1));
                    for (int d = 0; d < 2; ++d) {
                        if (tb->knots(d).unique().size() < 4) {
                            tb->insertKnot(T(1) / T(3), d, 1);
                            tb->insertKnot(T(2) / T(3), d, 1);
                        }
                    }
                }
            }
            harmonizeMultiPatchKnots(refMp);
            repairMicroFolds(refMp);
            optimizeWinslowBarrier(refMp);
            repairMicroFolds(refMp);
            if (refMp.nPatches() >= 2) {
                optimizeCoupledInterfacesAndInteriors(refMp, T(1e-7), 18, true);
                repairMicroFolds(refMp);
            }
            gsASG1LocalFitter<T>::unifyPatchOrientations(refMp);
            if (checkConforming(refMp) && countNegDense(refMp) < negPre)
                candMp = give(refMp);
        }

        if (static_cast<int>(candMp.nPatches()) > maxPatchesCap)
            return false;

        std::string vErr;
        if (!validateMultiPatch(m_domain, candMp, T(1e-4), &vErr, effMinAng, effMaxAng)) {
            if (cfDbg) {
                if (vErr.find("boundary") != std::string::npos) {
                    T bDist = 0;
                    verifyBoundaryMatch(m_domain, candMp, T(1e-4), &bDist);
                    std::cerr << "[CF_DBG] " << label << " P=" << candMp.nPatches()
                              << " FAIL validate: " << vErr << " maxBdrDist=" << bDist
                              << " IF=" << candMp.nInterfaces() << " B=" << candMp.nBoundary() << "\n";
                } else {
                    std::cerr << "[CF_DBG] " << label << " P=" << candMp.nPatches() << " FAIL validate: " << vErr << "\n";
                }
            }
            return false;
        }
        int nNeg = countNegDense(candMp);
        if (nNeg > 0) {
            if (cfDbg) std::cerr << "[CF_DBG] " << label << " P=" << candMp.nPatches() << " FAIL denseNeg=" << nNeg << "\n";
            return false;
        }

        T minDet = 0, minA = 0, maxA = 180;
        T energy = evaluateDistortionEnergy(candMp, &minDet, &minA, &maxA);
        if (minDet <= T(1e-8) || maxA > effMaxAng) {
            if (cfDbg) std::cerr << "[CF_DBG] " << label << " P=" << candMp.nPatches() << " FAIL post minDet=" << minDet << " maxA=" << maxA << "\n";
            return false;
        }
        T maxER = evalMultiPatchMaxEdgeRatio(candMp);
        if (cfDbg) std::cerr << "[CF_DBG] " << label << " P=" << candMp.nPatches() << " VALID! minDet=" << minDet << " ang=[" << minA << "," << maxA << "] maxER=" << maxER << "\n";

        SplittingCandidate cand;
        cand.multiPatch = give(candMp);
        cand.numPatches = static_cast<int>(cand.multiPatch.nPatches());
        cand.minJacobianDet = minDet;
        cand.distortionEnergy = energy;
        cand.minCornerAngleDeg = minA;
        cand.maxCornerAngleDeg = maxA;
        cand.methodUsed = label.empty() ? "HarmonicCrossField" : label;
        cand.isValid = true;

        if (!bestCand.isValid || scoreCand(cand) < scoreCand(bestCand))
            bestCand = give(cand);

        T bestER = evalMultiPatchMaxEdgeRatio(bestCand.multiPatch);
        return (bestCand.isValid && bestCand.numPatches <= compactTarget &&
                bestCand.maxCornerAngleDeg <= T(146.0) && bestCand.minCornerAngleDeg >= T(22.0) &&
                bestER <= T(5.5));
    };

    // Stage 2A: Poisson-Torsion Curved 5-Patch O-Grid for smooth 0-corner hole-free domains (e.g. circle, amoeba)
    if (m_domain.numHoles() == 0 && !vloops.empty() && vloops[0].corners.empty()) {
        gsMultiPatch<T> ogridMp;
        bool builtO = buildHarmonicCurvedOGrid(m_domain, vloops, fineMesh, accel, m_degree, ogridMp);
        if (cfDbg) std::cerr << "[CF_DBG] Stage2A PoissonOGrid: ok=" << builtO << " P=" << ogridMp.nPatches() << "\n";
        if (builtO) {
            if (evaluateAndRegister(ogridMp, "HarmonicCrossField"))
                return bestCand;
        }
    }

    // Stage 2B + Stage 3 (Primary Engine 1): Simultaneous Motorcycle T-Mesh + Native Integer Quantization
    for (int reflexMode : {0, 1}) {
        TMeshGraph<T> tmesh;
        bool builtT = buildSimultaneousMotorcycleTMesh(
            m_domain, vloops, fineMesh, cf, accel, m_degree, reflexMode, tmesh);
        if (cfDbg) {
            std::cerr << "[CF_DBG] Stage2B Motorcycle(rm=" << reflexMode << "): ok=" << builtT
                      << " F=" << tmesh.faces.size() << " E=" << tmesh.segs.size() << "\n";
        }
        if (!builtT) continue;
        bool quantOk = quantizeTMeshMinimalPatches(tmesh, tmeshARWeight, maxPatchesCap);
        if (cfDbg) {
            int estP = 0;
            for (const auto& rf : tmesh.faces) {
                int nu = 0, nv = 0;
                for (int he : rf.sideHEs[0]) nu += tmesh.segs[tmesh.halfEdges[he].segIdx].quantLen;
                for (int he : rf.sideHEs[1]) nv += tmesh.segs[tmesh.halfEdges[he].segIdx].quantLen;
                estP += nu * nv;
            }
            std::cerr << "[CF_DBG] Stage3 Quantizer(rm=" << reflexMode << "): ok=" << quantOk << " estP=" << estP << "\n";
        }
        if (!quantOk) continue;

        gsMultiPatch<T> motoMp;
        if (emitQuantizedTMeshPatches(tmesh, vloops, m_degree, motoMp)) {
            if (evaluateAndRegister(motoMp, "HarmonicCrossField"))
                return bestCand;
        }
    }

    // Engine 1b (Legacy fallback): MBO Curved-Streamline Separatrix Arrangement + Corridor Merging + Local Even-Gon Splitting
    {
        gsMultiPatch<T> sepMp;
        std::string label;
        if (buildCrossFieldSeparatrixFromSolved(
                m_domain, fineMesh, cf, accel, m_degree, stripRatio, angThresh, maxPatchesCap, sepMp, label)) {
            if (evaluateAndRegister(sepMp, label))
                return bestCand;
        }
    }

    // Engine 2: MBO Cross-Field-Aligned Coarse Blossom-Quad + Dual-Path Zipper + Curved Cross-Field Interior Edges
    runCrossFieldBlossomQuadEngine(
        m_domain, fineMesh, cf, accel, m_degree, minAngTol, maxAngTol, maxPatchesCap, evaluateAndRegister);

    return bestCand;
}

// ============================================================================
// Hole-Free Structural & Triangulation-Quad Decomposition Helpers
// ============================================================================

template <class T>
typename gsPatchSplitter<T>::SplittingCandidate
gsPatchSplitter<T>::computeHoleFreeStructuralSplitting()
{
    using namespace internal_splitter;

    SplittingCandidate bestCand;
    bestCand.methodUsed = "MedialAxisSkeleton";
    bestCand.isValid = false;

    const auto& outer = m_domain.outer();
    const int K = outer.size();
    if (K == 0) return bestCand;

    std::vector<CurveSeg<T>> origSegs;
    origSegs.reserve(K);
    for (int i = 0; i < K; ++i)
        origSegs.emplace_back(outer.curve(i), m_degree);

    std::vector<gsPoint2D<T>> outerPoly = sampleCurveLoopPoly(outer, 24);
    gsPoint2D<T> centroid = polygonCentroid(outerPoly);

    // Analyze outer vertices & sharp/reflex angles
    std::vector<gsPoint2D<T>> verts(K);
    std::vector<bool> isSharp(K, false), isReflex(K, false);
    int nSharp = 0, nReflex = 0;
    for (int i = 0; i < K; ++i) {
        verts[i] = origSegs[i].p0();
        const auto& prevBsp = *origSegs[(i + K - 1) % K].bsp;
        const auto& nextBsp = *origSegs[i].bsp;
        auto din  = evalCurveDeriv(prevBsp, T(1));
        auto dout = evalCurveDeriv(nextBsp, T(0));
        T n1 = std::hypot(din.x, din.y), n2 = std::hypot(dout.x, dout.y);
        if (n1 > 1e-14 && n2 > 1e-14) {
            T cr = (din.x * dout.y - din.y * dout.x) / (n1 * n2);
            T dt = (din.x * dout.x + din.y * dout.y) / (n1 * n2);
            T ang = std::atan2(cr, dt); // CCW: >0 convex turn, <0 reflex turn
            if (std::abs(ang) > T(0.10)) { // > ~5.7 deg
                isSharp[i] = true;
                ++nSharp;
                if (ang < -T(0.10)) {
                    isReflex[i] = true;
                    ++nReflex;
                }
            }
        }
    }

    auto tryFinalize = [&](gsMultiPatch<T>& mp, const std::string& label, SplittingCandidate& out) -> bool {
        if (mp.nPatches() == 0) return false;
        optimizeWinslowBarrier(mp);
        mp.computeTopology(T(1e-5), false, false);
        const T domMinAng = getDomainMinCornerAngleDeg();
        const T minAngTol = std::min(static_cast<T>(m_options.getReal("MinCornerAngleDeg")), T(0.85) * domMinAng);
        const T maxAngTol = static_cast<T>(m_options.getReal("MaxCornerAngleDeg"));
        if (!validateMultiPatch(m_domain, mp, T(1e-5), nullptr, minAngTol, maxAngTol)) return false;
        T minDet = 0, minA = 0, maxA = 180;
        T energy = evaluateDistortionEnergy(mp, &minDet, &minA, &maxA);
        out.multiPatch = mp;
        out.numPatches = static_cast<int>(mp.nPatches());
        out.minJacobianDet = minDet;
        out.distortionEnergy = energy;
        out.minCornerAngleDeg = minA;
        out.maxCornerAngleDeg = maxA;
        out.methodUsed = label;
        out.isValid = true;
        return true;
    };

    auto registerIfBetter = [&](const SplittingCandidate& cand) {
        if (!cand.isValid) return;
        if (!bestCand.isValid) { bestCand = cand; return; }
        T qNew = m_weightPatches * cand.numPatches + m_weightDistortion * cand.distortionEnergy;
        T qOld = m_weightPatches * bestCand.numPatches + m_weightDistortion * bestCand.distortionEnergy;
        if (qNew < qOld) bestCand = cand;
    };

    // ------------------------------------------------------------------------
    // 1. Single-patch check: "convexish" 4-curve hole-free quad -> 1 patch
    //    (A quad is convexish iff no corner lies in the convex hull of the other 3,
    //     which is equivalent to isStrictlyConvexQuad(c4); also allow convex domain
    //     corner angles up to ~174 deg with Winslow/refined-CP untangling for S-ducts.)
    // ------------------------------------------------------------------------
    if (K == 4) {
        gsPoint2D<T> c4[4] = {verts[0], verts[1], verts[2], verts[3]};
        if (isStrictlyConvexQuad(c4)) {
            const T radToDeg = T(180) / T(3.14159265358979323846);
            T domMinC = T(180), domMaxC = T(0);
            for (int i = 0; i < 4; ++i) {
                auto dPrev = internal_splitter::evalCurveDeriv(*origSegs[(i + 3) % 4].bsp, T(1));
                auto dCurr = internal_splitter::evalCurveDeriv(*origSegs[i].bsp, T(0));
                T nP = std::hypot(dPrev.x, dPrev.y), nC = std::hypot(dCurr.x, dCurr.y);
                if (nP > T(1e-12) && nC > T(1e-12)) {
                    T ux = -dPrev.x / nP, uy = -dPrev.y / nP;
                    T vx = dCurr.x / nC, vy = dCurr.y / nC;
                    T ang = std::atan2(vx * uy - vy * ux, ux * vx + uy * vy) * radToDeg;
                    if (ang <= T(0)) ang += T(360);
                    domMinC = std::min(domMinC, ang);
                    domMaxC = std::max(domMaxC, ang);
                }
            }
            if (domMinC > T(3) && domMaxC < T(174)) {
                const T localMinAng = std::min(static_cast<T>(m_options.getReal("MinCornerAngleDeg")), T(0.85) * domMinC);
                const T localMaxAng = std::max(static_cast<T>(m_options.getReal("MaxCornerAngleDeg")), std::min(T(175), domMaxC + T(4)));
                for (int rLevel = 0; rLevel <= 2; ++rLevel) {
                    gsMultiPatch<T> mp;
                    mp.addPatch(typename gsGeometry<T>::uPtr(
                        makeCoonsFromCCWCycle(origSegs[0], origSegs[1], origSegs[2], origSegs[3], m_degree)));
                    if (rLevel > 0) {
                        if (auto* tb = dynamic_cast<gsTensorBSpline<2, T>*>(&mp.patch(0)))
                            tb->uniformRefine(rLevel);
                    }
                    optimizeWinslowBarrier(mp);
                    mp.computeTopology(T(1e-5), false, false);
                    if (validateMultiPatch(m_domain, mp, T(1e-5), nullptr, localMinAng, localMaxAng)) {
                        T minDet = 0, minA = 0, maxA = 180;
                        T energy = evaluateDistortionEnergy(mp, &minDet, &minA, &maxA);
                        SplittingCandidate cand;
                        cand.multiPatch = give(mp);
                        cand.numPatches = 1;
                        cand.minJacobianDet = minDet;
                        cand.distortionEnergy = energy;
                        cand.minCornerAngleDeg = minA;
                        cand.maxCornerAngleDeg = maxA;
                        cand.methodUsed = "SingleConvexQuad";
                        cand.isValid = true;
                        return cand;
                    }
                }
            }
        }
    }

    // ------------------------------------------------------------------------
    // 2. O-Grid (5 patches): ideal for 1-corner / 2-corner shapes & 4-station shapes
    // ------------------------------------------------------------------------
    {
        std::vector<CurveSeg<T>> bEdges;
        if (K == 1) {
            auto sHalf = origSegs[0].splitArc(T(0.5));
            auto sQ1   = sHalf.first.splitArc(T(0.5));
            auto sQ2   = sHalf.second.splitArc(T(0.5));
            bEdges = {sQ1.first, sQ1.second, sQ2.first, sQ2.second};
        } else if (K == 2) {
            auto s0 = origSegs[0].splitArc(T(0.5));
            auto s1 = origSegs[1].splitArc(T(0.5));
            bEdges = {s0.first, s0.second, s1.first, s1.second};
        } else if (K == 4) {
            bEdges = origSegs;
        }

        if (bEdges.size() == 4) {
            gsPoint2D<T> P[4];
            for (int k = 0; k < 4; ++k) P[k] = bEdges[k].p0();
            gsPoint2D<T> pAvg((P[0].x + P[1].x + P[2].x + P[3].x) * T(0.25),
                              (P[0].y + P[1].y + P[2].y + P[3].y) * T(0.25));

            std::vector<gsPoint2D<T>> centers = {centroid, pAvg};
            // Also test shifting center away from any reflex station (e.g. heart cleft)
            for (int k = 0; k < 4; ++k) {
                gsPoint2D<T> opp = P[(k + 2) % 4];
                centers.emplace_back(T(0.75) * pAvg.x + T(0.25) * opp.x,
                                     T(0.75) * pAvg.y + T(0.25) * opp.y);
            }

            T alphas[3] = {T(0.45), T(0.35), T(0.55)};
            for (const auto& ctr : centers) {
                if (!pointInPolyStrict(ctr, outerPoly, T(1e-5))) continue;
                for (T alpha : alphas) {
                    gsPoint2D<T> Q[4];
                    bool qInside = true;
                    for (int k = 0; k < 4; ++k) {
                        Q[k] = gsPoint2D<T>(ctr.x + alpha * (P[k].x - ctr.x),
                                            ctr.y + alpha * (P[k].y - ctr.y));
                        if (!pointInPolyStrict(Q[k], outerPoly, T(1e-5))) { qInside = false; break; }
                    }
                    if (!qInside || !isStrictlyConvexQuad(Q)) continue;

                    CurveSeg<T> qEdge[4], spoke[4];
                    for (int k = 0; k < 4; ++k) {
                        qEdge[k] = makeStraightSeg(Q[k], Q[(k + 1) % 4], m_degree);
                        spoke[k] = makeStraightSeg(Q[k], P[k], m_degree);
                    }

                    gsMultiPatch<T> mp;
                    mp.addPatch(typename gsGeometry<T>::uPtr(
                        makeCoonsFromCCWCycle(qEdge[0], qEdge[1], qEdge[2], qEdge[3], m_degree)));
                    for (int k = 0; k < 4; ++k) {
                        int kn = (k + 1) % 4;
                        mp.addPatch(typename gsGeometry<T>::uPtr(
                            makeCoonsFromCCWCycle(qEdge[k], spoke[kn], bEdges[k].reversed(), spoke[k].reversed(), m_degree)));
                    }

                    SplittingCandidate cand;
                    if (tryFinalize(mp, "MedialOGrid", cand)) {
                        registerIfBetter(cand);
                    }
                }
            }
        }
    }

    // ------------------------------------------------------------------------
    // 3. Generalized Ribbon / Spine Decomposition (for 2-curve crescents/lenses)
    // ------------------------------------------------------------------------
    if (K == 2) {
        int mSegsList[2] = {3, 4};
        T tipShifts[2] = {T(0.18), T(0.28)};
        for (int mSegs : mSegsList) {
            std::vector<CurveSeg<T>> botSegs(mSegs), topSegs(mSegs);
            CurveSeg<T> rem0 = origSegs[0];
            CurveSeg<T> rem1 = origSegs[1].reversed(); // Both run V0 -> V1
            for (int s = 0; s < mSegs - 1; ++s) {
                T frac = T(1) / T(mSegs - s);
                auto sp0 = rem0.splitArc(frac);
                botSegs[s] = sp0.first; rem0 = sp0.second;
                auto sp1 = rem1.splitArc(frac);
                topSegs[s] = sp1.first; rem1 = sp1.second;
            }
            botSegs[mSegs - 1] = rem0;
            topSegs[mSegs - 1] = rem1;

            int nMid = mSegs - 1;
            for (T delta : tipShifts) {
                std::vector<gsPoint2D<T>> M(nMid);
                for (int i = 0; i < nMid; ++i) {
                    auto pb = botSegs[i].p1();
                    auto pt = topSegs[i].p1();
                    M[i] = gsPoint2D<T>(T(0.5) * (pb.x + pt.x), T(0.5) * (pb.y + pt.y));
                }
                if (nMid >= 2) {
                    gsPoint2D<T> M0_orig = M[0], Mlast_orig = M[nMid - 1];
                    M[0]        = gsPoint2D<T>((1 - delta) * M0_orig.x + delta * M[1].x,
                                               (1 - delta) * M0_orig.y + delta * M[1].y);
                    M[nMid - 1] = gsPoint2D<T>((1 - delta) * Mlast_orig.x + delta * M[nMid - 2].x,
                                               (1 - delta) * Mlast_orig.y + delta * M[nMid - 2].y);
                }

                std::vector<CurveSeg<T>> bSpoke(nMid), tSpoke(nMid), spine(nMid - 1);
                for (int i = 0; i < nMid; ++i) {
                    bSpoke[i] = makeStraightSeg(botSegs[i].p1(), M[i], m_degree);
                    tSpoke[i] = makeStraightSeg(M[i], topSegs[i].p1(), m_degree);
                }
                for (int i = 0; i < nMid - 1; ++i)
                    spine[i] = makeStraightSeg(M[i], M[i + 1], m_degree);

                gsMultiPatch<T> mp;
                // Left tip kite: (V0, B1, M0, T1)
                mp.addPatch(typename gsGeometry<T>::uPtr(
                    makeCoonsFromCCWCycle(botSegs[0], bSpoke[0], tSpoke[0], topSegs[0].reversed(), m_degree)));
                // Middle pairs
                for (int i = 0; i < nMid - 1; ++i) {
                    mp.addPatch(typename gsGeometry<T>::uPtr(
                        makeCoonsFromCCWCycle(botSegs[i + 1], bSpoke[i + 1], spine[i].reversed(), bSpoke[i].reversed(), m_degree)));
                    mp.addPatch(typename gsGeometry<T>::uPtr(
                        makeCoonsFromCCWCycle(spine[i], tSpoke[i + 1], topSegs[i + 1].reversed(), tSpoke[i].reversed(), m_degree)));
                }
                // Right tip kite: (M_last, B_last, V1, T_last)
                mp.addPatch(typename gsGeometry<T>::uPtr(
                    makeCoonsFromCCWCycle(bSpoke[nMid - 1].reversed(), botSegs[mSegs - 1],
                                          topSegs[mSegs - 1].reversed(), tSpoke[nMid - 1].reversed(), m_degree)));

                SplittingCandidate cand;
                if (tryFinalize(mp, "MedialRibbon", cand)) {
                    registerIfBetter(cand);
                }
            }
        }
    }

    // ------------------------------------------------------------------------
    // 3b. Curved Channel / Serpentine Ribbon Decomposition (for K == 4 ducts)
    //     Splits along opposite boundary curves into M stations, bisecting any
    //     obtuse (> 140 deg) corners cleanly with Catmull-Clark edge halving.
    // ------------------------------------------------------------------------
    if (K == 4) {
        std::vector<gsPoint2D<T>> densePoly = sampleCurveLoopPoly(outer, 96);
        auto segStrictlyInsideDense = [&](const gsPoint2D<T>& a, const gsPoint2D<T>& b) -> bool {
            const int np = static_cast<int>(densePoly.size());
            const T abx = b.x - a.x, aby = b.y - a.y;
            for (int i = 0; i < np; ++i) {
                const auto& c = densePoly[i];
                const auto& d = densePoly[(i + 1) % np];
                const T cdx = d.x - c.x, cdy = d.y - c.y;
                const T denom = abx * cdy - aby * cdx;
                if (std::abs(denom) <= T(1e-14)) continue;
                const T t = ((c.x - a.x) * cdy - (c.y - a.y) * cdx) / denom;
                const T u = ((c.x - a.x) * aby - (c.y - a.y) * abx) / denom;
                if (t > T(0.015) && t < T(0.985) && u >= -T(1e-6) && u <= T(1) + T(1e-6))
                    return false;
            }
            const T sampleT[5] = {T(0.1), T(0.25), T(0.5), T(0.75), T(0.9)};
            for (T t : sampleT) {
                gsPoint2D<T> mid((T(1) - t) * a.x + t * b.x, (T(1) - t) * a.y + t * b.y);
                if (!pointInPolyStrict(mid, densePoly, T(1e-6))) return false;
            }
            return true;
        };

        auto emitTri1to3 = [&](gsMultiPatch<T>& mp,
                               const gsPoint2D<T>& A, const gsPoint2D<T>& B, const gsPoint2D<T>& C,
                               const CurveSeg<T>& eAB0, const CurveSeg<T>& eAB1, const gsPoint2D<T>& mAB,
                               const CurveSeg<T>& eBC0, const CurveSeg<T>& eBC1, const gsPoint2D<T>& mBC,
                               const CurveSeg<T>& eCA0, const CurveSeg<T>& eCA1, const gsPoint2D<T>& mCA) -> bool
        {
            gsPoint2D<T> O_mid((mAB.x + mBC.x + mCA.x) / T(3), (mAB.y + mBC.y + mCA.y) / T(3));
            gsPoint2D<T> O_tri((A.x + B.x + C.x) / T(3), (A.y + B.y + C.y) / T(3));
            auto check3 = [&](const gsPoint2D<T>& Oc) -> bool {
                gsPoint2D<T> q0[4] = {A, mAB, Oc, mCA};
                gsPoint2D<T> q1[4] = {B, mBC, Oc, mAB};
                gsPoint2D<T> q2[4] = {C, mCA, Oc, mBC};
                return isStrictlyConvexQuad(q0) && isStrictlyConvexQuad(q1) && isStrictlyConvexQuad(q2);
            };
            gsPoint2D<T> O = O_mid;
            if (!check3(O)) {
                if (check3(O_tri)) O = O_tri;
                else return false;
            }
            auto s0 = makeStraightSeg(O, mAB, m_degree);
            auto s1 = makeStraightSeg(O, mBC, m_degree);
            auto s2 = makeStraightSeg(O, mCA, m_degree);
            mp.addPatch(typename gsGeometry<T>::uPtr(makeCoonsFromCCWCycle(eAB0, s0.reversed(), s2, eCA1, m_degree)));
            mp.addPatch(typename gsGeometry<T>::uPtr(makeCoonsFromCCWCycle(eBC0, s1.reversed(), s0, eAB1, m_degree)));
            mp.addPatch(typename gsGeometry<T>::uPtr(makeCoonsFromCCWCycle(eCA0, s2.reversed(), s1, eBC1, m_degree)));
            return true;
        };

        for (int base = 0; base < 2; ++base) {
            CurveSeg<T> botSeg   = origSegs[base];
            CurveSeg<T> rightSeg = origSegs[(base + 1) % 4];
            CurveSeg<T> topSeg   = origSegs[(base + 2) % 4].reversed(); // L[M] -> R[M]
            CurveSeg<T> leftSeg  = origSegs[(base + 3) % 4].reversed(); // L[0] -> L[M]

            const int mVals[4] = {4, 6, 8, 10};
            for (int mIdx = 0; mIdx < 4; ++mIdx) {
                const int M = mVals[mIdx];
                std::vector<CurveSeg<T>> lSegs(M), rSegs(M);
                {
                    CurveSeg<T> remL = leftSeg, remR = rightSeg;
                    for (int j = 0; j < M - 1; ++j) {
                        T frac = T(1) / static_cast<T>(M - j);
                        auto spL = remL.splitArc(frac);
                        auto spR = remR.splitArc(frac);
                        lSegs[j] = spL.first;  remL = spL.second;
                        rSegs[j] = spR.first;  remR = spR.second;
                    }
                    lSegs[M - 1] = remL;
                    rSegs[M - 1] = remR;
                }

                std::vector<gsPoint2D<T>> Lpts(M + 1), Rpts(M + 1);
                std::vector<CurveSeg<T>> hSegs(M + 1);
                Lpts[0] = botSeg.p0(); Rpts[0] = botSeg.p1(); hSegs[0] = botSeg;
                Lpts[M] = topSeg.p0(); Rpts[M] = topSeg.p1(); hSegs[M] = topSeg;

                bool cutsValid = true;
                for (int j = 1; j < M; ++j) {
                    Lpts[j] = lSegs[j].p0();
                    Rpts[j] = rSegs[j].p0();
                    if (!segStrictlyInsideDense(Lpts[j], Rpts[j])) {
                        cutsValid = false;
                        break;
                    }
                    hSegs[j] = makeStraightSeg(Lpts[j], Rpts[j], m_degree);
                }
                if (!cutsValid) continue;

                // Option A: Direct M-quad channel strip
                {
                    gsMultiPatch<T> mpA;
                    for (int j = 0; j < M; ++j) {
                        mpA.addPatch(typename gsGeometry<T>::uPtr(
                            makeCoonsFromCCWCycle(hSegs[j], rSegs[j], hSegs[j + 1].reversed(), lSegs[j].reversed(), m_degree)));
                    }
                    SplittingCandidate candA;
                    if (tryFinalize(mpA, "CurvedChannelStrip", candA))
                        registerIfBetter(candA);
                }

                // Option B: Conforming diagonal split of obtuse (> 140 deg) end/strip quads + Catmull-Clark midpoint quadification
                std::vector<CurveSeg<T>> h0(M + 1), h1(M + 1), l0(M), l1(M), r0(M), r1(M);
                std::vector<gsPoint2D<T>> hm(M + 1), lm(M), rm(M);
                for (int j = 0; j <= M; ++j) {
                    auto sp = hSegs[j].splitArc(T(0.5));
                    h0[j] = sp.first; h1[j] = sp.second; hm[j] = h0[j].p1();
                }
                for (int j = 0; j < M; ++j) {
                    auto spL = lSegs[j].splitArc(T(0.5));
                    auto spR = rSegs[j].splitArc(T(0.5));
                    l0[j] = spL.first; l1[j] = spL.second; lm[j] = l0[j].p1();
                    r0[j] = spR.first; r1[j] = spR.second; rm[j] = r0[j].p1();
                }

                gsMultiPatch<T> mpB;
                bool okB = true;
                for (int j = 0; j < M && okB; ++j) {
                    std::unique_ptr<gsTensorBSpline<2, T>> testP(
                        makeCoonsFromCCWCycle(hSegs[j], rSegs[j], hSegs[j + 1].reversed(), lSegs[j].reversed(), m_degree));
                    T angs[4], minA = 0, maxA = 180;
                    bool degen = false;
                    evalPatchCornerAnglesDeg(*testP, angs, minA, maxA, degen);

                    gsPoint2D<T> P[4] = {Lpts[j], Rpts[j], Rpts[j + 1], Lpts[j + 1]};
                    int bestObtuse = -1;
                    T worstA = T(140);
                    for (int c = 0; c < 4; ++c) {
                        if (angs[c] > worstA) { worstA = angs[c]; bestObtuse = c; }
                    }

                    if (bestObtuse >= 0 && !degen) {
                        if (bestObtuse == 1 || bestObtuse == 3) {
                            // Diagonal P1 -> P3
                            if (!segStrictlyInsideDense(P[1], P[3])) { okB = false; break; }
                            auto dSeg = makeStraightSeg(P[1], P[3], m_degree);
                            auto dSp  = dSeg.splitArc(T(0.5));
                            gsPoint2D<T> dm = dSp.first.p1();
                            if (!emitTri1to3(mpB, P[0], P[1], P[3],
                                             h0[j], h1[j], hm[j],
                                             dSp.first, dSp.second, dm,
                                             l1[j].reversed(), l0[j].reversed(), lm[j]) ||
                                !emitTri1to3(mpB, P[1], P[2], P[3],
                                             r0[j], r1[j], rm[j],
                                             h1[j + 1].reversed(), h0[j + 1].reversed(), hm[j + 1],
                                             dSp.second.reversed(), dSp.first.reversed(), dm)) {
                                okB = false;
                            }
                        } else {
                            // Diagonal P0 -> P2
                            if (!segStrictlyInsideDense(P[0], P[2])) { okB = false; break; }
                            auto dSeg = makeStraightSeg(P[0], P[2], m_degree);
                            auto dSp  = dSeg.splitArc(T(0.5));
                            gsPoint2D<T> dm = dSp.first.p1();
                            if (!emitTri1to3(mpB, P[0], P[1], P[2],
                                             h0[j], h1[j], hm[j],
                                             r0[j], r1[j], rm[j],
                                             dSp.second.reversed(), dSp.first.reversed(), dm) ||
                                !emitTri1to3(mpB, P[0], P[2], P[3],
                                             dSp.first, dSp.second, dm,
                                             h1[j + 1].reversed(), h0[j + 1].reversed(), hm[j + 1],
                                             l1[j].reversed(), l0[j].reversed(), lm[j])) {
                                okB = false;
                            }
                        }
                    } else {
                        gsPoint2D<T> qCorners[4] = {Lpts[j], Rpts[j], Rpts[j + 1], Lpts[j + 1]};
                        if (!isStrictlyConvexQuad(qCorners)) { okB = false; break; }
                        gsPoint2D<T> O((hm[j].x + rm[j].x + hm[j + 1].x + lm[j].x) * T(0.25),
                                       (hm[j].y + rm[j].y + hm[j + 1].y + lm[j].y) * T(0.25));
                        auto sBot   = makeStraightSeg(O, hm[j], m_degree);
                        auto sRight = makeStraightSeg(O, rm[j], m_degree);
                        auto sTop   = makeStraightSeg(O, hm[j + 1], m_degree);
                        auto sLeft  = makeStraightSeg(O, lm[j], m_degree);
                        mpB.addPatch(typename gsGeometry<T>::uPtr(
                            makeCoonsFromCCWCycle(h0[j], sBot.reversed(), sLeft, l0[j].reversed(), m_degree)));
                        mpB.addPatch(typename gsGeometry<T>::uPtr(
                            makeCoonsFromCCWCycle(h1[j], r0[j], sRight.reversed(), sBot, m_degree)));
                        mpB.addPatch(typename gsGeometry<T>::uPtr(
                            makeCoonsFromCCWCycle(sRight, r1[j], h1[j + 1].reversed(), sTop.reversed(), m_degree)));
                        mpB.addPatch(typename gsGeometry<T>::uPtr(
                            makeCoonsFromCCWCycle(sLeft.reversed(), sTop, h0[j + 1].reversed(), l1[j].reversed(), m_degree)));
                    }
                }
                if (okB && mpB.nPatches() > 0) {
                    SplittingCandidate candB;
                    if (tryFinalize(mpB, "CurvedChannelSplit", candB))
                        registerIfBetter(candB);
                }
            }
        }
    }

    // ------------------------------------------------------------------------
    // 3c. Single-Reflex Dart / Boomerang Decomposition (K == 4, nReflex == 1)
    //     Splits edges incident to opposite corner opp=(r+2)%4 at frac t,
    //     creating Quad(opp, mA, Vr, mB) + 2 Triangles, then 1-step Catmull-Clark -> 10 quads.
    // ------------------------------------------------------------------------
    if (K == 4 && nReflex == 1) {
        int r = -1;
        for (int i = 0; i < 4; ++i) if (isReflex[i]) { r = i; break; }
        if (r >= 0) {
            int opp = (r + 2) % 4;
            int v1  = (opp + 1) % 4;
            int v3  = (r + 1) % 4;
            const CurveSeg<T>& segOA = origSegs[opp];     // opp -> v1
            const CurveSeg<T>& seg1r = origSegs[v1];      // v1  -> Vr
            const CurveSeg<T>& segr3 = origSegs[r];       // Vr  -> v3
            const CurveSeg<T>& seg3O = origSegs[v3];      // v3  -> opp

            std::vector<gsPoint2D<T>> densePoly = sampleCurveLoopPoly(outer, 48);
            const T fracs[4] = {T(0.50), T(0.35), T(0.45), T(0.25)};
            for (T frac : fracs) {
                auto spA = segOA.splitArc(frac);          // opp -> mA -> v1
                auto spB = seg3O.splitArc(T(1) - frac);   // v3  -> mB -> opp
                gsPoint2D<T> mA = spA.first.p1();
                gsPoint2D<T> mB = spB.first.p1();
                gsPoint2D<T> Vr = verts[r];
                gsPoint2D<T> Vopp = verts[opp];
                gsPoint2D<T> V1 = verts[v1];
                gsPoint2D<T> V3 = verts[v3];

                if (!isSegmentInsidePoly(mA, Vr, densePoly) || !isSegmentInsidePoly(mB, Vr, densePoly)) continue;

                auto cA = makeStraightSeg(mA, Vr, m_degree);
                auto cB = makeStraightSeg(mB, Vr, m_degree);
                auto spEA0 = spA.first.splitArc(T(0.5));
                auto spEA1 = spA.second.splitArc(T(0.5));
                auto sp1r  = seg1r.splitArc(T(0.5));
                auto spr3  = segr3.splitArc(T(0.5));
                auto spEB0 = spB.first.splitArc(T(0.5));
                auto spEB1 = spB.second.splitArc(T(0.5));
                auto spCA  = cA.splitArc(T(0.5));
                auto spCB  = cB.splitArc(T(0.5));

                auto m_eA0 = spEA0.first.p1(), m_eA1 = spEA1.first.p1();
                auto m_1r  = sp1r.first.p1(),  m_r3  = spr3.first.p1();
                auto m_eB0 = spEB0.first.p1(), m_eB1 = spEB1.first.p1();
                auto m_cA  = spCA.first.p1(),  m_cB  = spCB.first.p1();

                auto addTri1to3 = [&](gsMultiPatch<T>& mpOut,
                                      const gsPoint2D<T>& A, const gsPoint2D<T>& B, const gsPoint2D<T>& C,
                                      const CurveSeg<T>& eAB0, const CurveSeg<T>& eAB1, const gsPoint2D<T>& mAB,
                                      const CurveSeg<T>& eBC0, const CurveSeg<T>& eBC1, const gsPoint2D<T>& mBC,
                                      const CurveSeg<T>& eCA0, const CurveSeg<T>& eCA1, const gsPoint2D<T>& mCA) -> bool
                {
                    gsPoint2D<T> O;
                    if (!internal_splitter::pickOptimalTri1to3Center(A, B, C, mAB, mBC, mCA, O))
                        return false;
                    auto s0 = makeStraightSeg(O, mAB, m_degree);
                    auto s1 = makeStraightSeg(O, mBC, m_degree);
                    auto s2 = makeStraightSeg(O, mCA, m_degree);
                    mpOut.addPatch(typename gsGeometry<T>::uPtr(makeCoonsFromCCWCycle(eAB0, s0.reversed(), s2, eCA1, m_degree)));
                    mpOut.addPatch(typename gsGeometry<T>::uPtr(makeCoonsFromCCWCycle(eBC0, s1.reversed(), s0, eAB1, m_degree)));
                    mpOut.addPatch(typename gsGeometry<T>::uPtr(makeCoonsFromCCWCycle(eCA0, s2.reversed(), s1, eBC1, m_degree)));
                    return true;
                };

                // First try 4-triangle split (diagonal Vopp -> Vr) into 12 quads
                if (isSegmentInsidePoly(Vopp, Vr, densePoly)) {
                    auto cDiag = makeStraightSeg(Vopp, Vr, m_degree);
                    auto spDiag = cDiag.splitArc(T(0.5));
                    auto m_diag = spDiag.first.p1();
                    gsMultiPatch<T> mp12;
                    bool ok4 =
                        addTri1to3(mp12, Vopp, mA, Vr,
                                   spEA0.first, spEA0.second, m_eA0,
                                   spCA.first, spCA.second, m_cA,
                                   spDiag.second.reversed(), spDiag.first.reversed(), m_diag) &&
                        addTri1to3(mp12, Vopp, Vr, mB,
                                   spDiag.first, spDiag.second, m_diag,
                                   spCB.second.reversed(), spCB.first.reversed(), m_cB,
                                   spEB1.first, spEB1.second, m_eB1) &&
                        addTri1to3(mp12, mA, V1, Vr,
                                   spEA1.first, spEA1.second, m_eA1,
                                   sp1r.first, sp1r.second, m_1r,
                                   spCA.second.reversed(), spCA.first.reversed(), m_cA) &&
                        addTri1to3(mp12, Vr, V3, mB,
                                   spr3.first, spr3.second, m_r3,
                                   spEB0.first, spEB0.second, m_eB0,
                                   spCB.first, spCB.second, m_cB);
                    if (ok4) {
                        SplittingCandidate cand12;
                        if (tryFinalize(mp12, "DartSplit12", cand12))
                            registerIfBetter(cand12);
                    }
                }

                gsPoint2D<T> qCorners[4] = {Vopp, mA, Vr, mB};
                if (!isStrictlyConvexQuad(qCorners)) continue;

                gsMultiPatch<T> mpD;
                gsPoint2D<T> Oq((m_eA0.x + m_cA.x + m_cB.x + m_eB1.x) * T(0.25),
                                (m_eA0.y + m_cA.y + m_cB.y + m_eB1.y) * T(0.25));
                auto qS0 = makeStraightSeg(Oq, m_eA0, m_degree);
                auto qS1 = makeStraightSeg(Oq, m_cA, m_degree);
                auto qS2 = makeStraightSeg(Oq, m_cB, m_degree);
                auto qS3 = makeStraightSeg(Oq, m_eB1, m_degree);
                mpD.addPatch(typename gsGeometry<T>::uPtr(
                    makeCoonsFromCCWCycle(spEA0.first, qS0.reversed(), qS3, spEB1.second, m_degree)));
                mpD.addPatch(typename gsGeometry<T>::uPtr(
                    makeCoonsFromCCWCycle(spEA0.second, spCA.first, qS1.reversed(), qS0, m_degree)));
                mpD.addPatch(typename gsGeometry<T>::uPtr(
                    makeCoonsFromCCWCycle(qS1, spCA.second, spCB.second.reversed(), qS2.reversed(), m_degree)));
                mpD.addPatch(typename gsGeometry<T>::uPtr(
                    makeCoonsFromCCWCycle(qS3.reversed(), qS2, spCB.first.reversed(), spEB1.first, m_degree)));

                bool okTriA = addTri1to3(mpD, mA, V1, Vr,
                                         spEA1.first, spEA1.second, m_eA1,
                                         sp1r.first, sp1r.second, m_1r,
                                         spCA.second.reversed(), spCA.first.reversed(), m_cA);
                bool okTriB = addTri1to3(mpD, Vr, V3, mB,
                                         spr3.first, spr3.second, m_r3,
                                         spEB0.first, spEB0.second, m_eB0,
                                         spCB.first, spCB.second, m_cB);
                if (okTriA && okTriB) {
                    SplittingCandidate candD;
                    if (tryFinalize(mpD, "DartSplit10", candD))
                        registerIfBetter(candD);
                }
            }
        }
    }

    // ------------------------------------------------------------------------
    // 4. L-Shape Single-Reflex Cut (3 patches, generalized to any reflex index)
    // ------------------------------------------------------------------------
    if (K == 6 && nReflex == 1) {
        int r = -1;
        for (int i = 0; i < K; ++i) if (isReflex[i]) { r = i; break; }
        if (r >= 0) {
            int eA = (r + 2) % 6;
            int eB = (r + 3) % 6;
            auto spA = origSegs[eA].splitArc(T(0.5));
            auto spB = origSegs[eB].splitArc(T(0.5));
            gsPoint2D<T> Vr = verts[r];
            auto cutA = makeStraightSeg(spA.first.p1(), Vr, m_degree);
            auto cutB = makeStraightSeg(spB.first.p1(), Vr, m_degree);

            gsMultiPatch<T> mp;
            // Corner block: (V_{eB+0.5}, V_{r+4}, V_r, V_{eA+0.5}) -> wait trace edges cleanly:
            // spB.first ends at mB; spB.second starts at mB and ends at V_{eB+1}=V_{r+4}
            mp.addPatch(typename gsGeometry<T>::uPtr(
                makeCoonsFromCCWCycle(spB.second, origSegs[(r + 4) % 6], origSegs[(r + 5) % 6], cutB.reversed(), m_degree)));
            // Second arm block: (mA, cutA, cutB_rev, spB.first) -> wait mB -> mA:
            mp.addPatch(typename gsGeometry<T>::uPtr(
                makeCoonsFromCCWCycle(spB.first, cutB, cutA.reversed(), spA.second, m_degree)));
            // Wait check cycle order: spA.second (mA->V_{eA+1}=mB? No, eB=eA+1! So spA.second ends at V_{eB}=start of spB.first!):
            // Corner at mB: spB.first (V_{eB} -> mB), cutB (mB -> Vr), cutA.reversed() (Vr -> mA), spA.second (mA -> V_{eB})! That IS CCW!
            // Third block: (V_r, origSegs[r], origSegs[(r+1)%6], spA.first, cutA):
            mp.addPatch(typename gsGeometry<T>::uPtr(
                makeCoonsFromCCWCycle(origSegs[r], origSegs[(r + 1) % 6], spA.first, cutA, m_degree)));

            SplittingCandidate cand;
            if (tryFinalize(mp, "ReflexCut3Patch", cand)) {
                registerIfBetter(cand);
            }
        }
    }

    // ------------------------------------------------------------------------
    // 5a. Alternating Star Kites (K even >= 6, K/2 patches, spokes to reflex vertices)
    // ------------------------------------------------------------------------
    if (K >= 6 && (K % 2 == 0) && nReflex == K / 2 && pointInPolyStrict(centroid, outerPoly, T(1e-5))) {
        int startIdx = isReflex[0] ? 0 : 1;
        bool altOk = true;
        for (int i = 0; i < K; ++i) {
            bool expectRef = ((i % 2) == startIdx);
            if (isReflex[i] != expectRef) { altOk = false; break; }
        }
        if (altOk) {
            std::vector<CurveSeg<T>> spokes(K / 2);
            for (int j = 0; j < K / 2; ++j) {
                int vIdx = (startIdx + 2 * j) % K;
                spokes[j] = makeStraightSeg(centroid, verts[vIdx], m_degree);
            }
            gsMultiPatch<T> mp;
            for (int j = 0; j < K / 2; ++j) {
                int jn = (j + 1) % (K / 2);
                int e0 = (startIdx + 2 * j) % K;
                int e1 = (e0 + 1) % K;
                mp.addPatch(typename gsGeometry<T>::uPtr(
                    makeCoonsFromCCWCycle(spokes[j], origSegs[e0], origSegs[e1], spokes[jn].reversed(), m_degree)));
            }
            SplittingCandidate cand;
            if (tryFinalize(mp, "AlternatingStarKites", cand)) {
                registerIfBetter(cand);
            }
        }
    }

    // ------------------------------------------------------------------------
    // 5b. Centroid-to-Midpoint Star Fan (K patches, for 0-reflex star domains)
    // ------------------------------------------------------------------------
    if (K >= 3 && nReflex == 0 && pointInPolyStrict(centroid, outerPoly, T(1e-5))) {
        std::vector<CurveSeg<T>> h0(K), h1(K), spokes(K);
        for (int i = 0; i < K; ++i) {
            auto sp = origSegs[i].splitArc(T(0.5));
            h0[i] = sp.first;
            h1[i] = sp.second;
            spokes[i] = makeStraightSeg(centroid, h0[i].p1(), m_degree);
        }
        gsMultiPatch<T> mp;
        for (int i = 0; i < K; ++i) {
            int iprev = (i + K - 1) % K;
            mp.addPatch(typename gsGeometry<T>::uPtr(
                makeCoonsFromCCWCycle(spokes[iprev], h1[iprev], h0[i], spokes[i].reversed(), m_degree)));
        }
        SplittingCandidate cand;
        if (tryFinalize(mp, "CentroidStarFan", cand)) {
            registerIfBetter(cand);
        }
    }

    // ------------------------------------------------------------------------
    // 5c. Generalized Reflex-Sector Star Kite Decomposition
    //     Connects an interior center C to all reflex vertices; sectors of
    //     length s=1,2 become single kite quads, and longer sectors subdivide
    //     via interior spokes to boundary edge midpoints (e.g. 7-star -> 4 quads).
    // ------------------------------------------------------------------------
    if (nReflex >= 2) {
        std::vector<int> rIdx;
        for (int i = 0; i < K; ++i) {
            if (isReflex[i]) rIdx.push_back(i);
        }
        gsPoint2D<T> rAvg(0, 0);
        for (int ri : rIdx) {
            rAvg.x += verts[ri].x;
            rAvg.y += verts[ri].y;
        }
        rAvg.x /= static_cast<T>(rIdx.size());
        rAvg.y /= static_cast<T>(rIdx.size());

        std::vector<gsPoint2D<T>> starCenters = {rAvg, centroid};
        for (const auto& ctr : starCenters) {
            if (!pointInPolyStrict(ctr, outerPoly, T(1e-5))) continue;

            gsMultiPatch<T> mp;
            bool okStar = true;
            std::vector<CurveSeg<T>> rSpokes(rIdx.size());
            for (size_t a = 0; a < rIdx.size(); ++a) {
                if (!isSegmentInsidePoly(ctr, verts[rIdx[a]], outerPoly)) {
                    okStar = false;
                    break;
                }
                rSpokes[a] = makeStraightSeg(ctr, verts[rIdx[a]], m_degree);
            }
            if (!okStar) continue;

            for (size_t a = 0; a < rIdx.size() && okStar; ++a) {
                size_t b = (a + 1) % rIdx.size();
                int rA = rIdx[a];
                int rB = rIdx[b];
                int span = (rB - rA + K) % K;
                if (span == 0) span = K;

                if (span == 1) {
                    auto sp = origSegs[rA].splitArc(T(0.5));
                    mp.addPatch(typename gsGeometry<T>::uPtr(
                        makeCoonsFromCCWCycle(rSpokes[a], sp.first, sp.second, rSpokes[b].reversed(), m_degree)));
                } else if (span == 2) {
                    mp.addPatch(typename gsGeometry<T>::uPtr(
                        makeCoonsFromCCWCycle(rSpokes[a], origSegs[rA], origSegs[(rA + 1) % K], rSpokes[b].reversed(), m_degree)));
                } else {
                    // Walk sector in 2-edge steps; if 3 edges remain at end, split middle edge.
                    int pos = 0;
                    CurveSeg<T> curSpoke = rSpokes[a];
                    while (span - pos > 3) {
                        int e0 = (rA + pos) % K;
                        int e1 = (rA + pos + 1) % K;
                        int vNext = (rA + pos + 2) % K;
                        if (!isSegmentInsidePoly(ctr, verts[vNext], outerPoly)) {
                            okStar = false;
                            break;
                        }
                        CurveSeg<T> nextSpoke = makeStraightSeg(ctr, verts[vNext], m_degree);
                        mp.addPatch(typename gsGeometry<T>::uPtr(
                            makeCoonsFromCCWCycle(curSpoke, origSegs[e0], origSegs[e1], nextSpoke.reversed(), m_degree)));
                        curSpoke = nextSpoke;
                        pos += 2;
                    }
                    if (!okStar) break;
                    int rem = span - pos;
                    if (rem == 2) {
                        int e0 = (rA + pos) % K;
                        int e1 = (rA + pos + 1) % K;
                        mp.addPatch(typename gsGeometry<T>::uPtr(
                            makeCoonsFromCCWCycle(curSpoke, origSegs[e0], origSegs[e1], rSpokes[b].reversed(), m_degree)));
                    } else if (rem == 3) {
                        int e0 = (rA + pos) % K;
                        int e1 = (rA + pos + 1) % K;
                        int e2 = (rA + pos + 2) % K;
                        auto spMid = origSegs[e1].splitArc(T(0.5));
                        gsPoint2D<T> mPt = spMid.first.p1();
                        if (!isSegmentInsidePoly(ctr, mPt, outerPoly)) {
                            okStar = false;
                            break;
                        }
                        CurveSeg<T> mSpoke = makeStraightSeg(ctr, mPt, m_degree);
                        mp.addPatch(typename gsGeometry<T>::uPtr(
                            makeCoonsFromCCWCycle(curSpoke, origSegs[e0], spMid.first, mSpoke.reversed(), m_degree)));
                        mp.addPatch(typename gsGeometry<T>::uPtr(
                            makeCoonsFromCCWCycle(mSpoke, spMid.second, origSegs[e2], rSpokes[b].reversed(), m_degree)));
                    } else {
                        okStar = false;
                    }
                }
            }
            if (!okStar) continue;

            SplittingCandidate cand;
            if (tryFinalize(mp, "ReflexStarKites", cand)) {
                registerIfBetter(cand);
            }
        }
    }

    // ------------------------------------------------------------------------
    // 6. Comb / Fingered-Channel Spine Decomposition
    //    Detects T-tooth comb domains (K = 4*T, T >= 2) and tries both:
    //    (a) 3T - 3 conforming patches (diagonal end teeth -> 9 patches for T=4)
    //    (b) 3T - 1 conforming patches (orthogonal grid -> 11 patches for T=4)
    // ------------------------------------------------------------------------
    if (K >= 8 && (K % 4 == 0) && nReflex == 2 * (K / 4 - 1)) {
        const int Tcnt = K / 4;
        for (int kBot = 0; kBot < K; ++kBot) {
            auto idxOf = [&](int off) -> int { return (kBot + off) % K; };
            bool patOk = !isReflex[idxOf(0)] && !isReflex[idxOf(1)] &&
                         !isReflex[idxOf(2)] && !isReflex[idxOf(3)];
            for (int m = 0; m < Tcnt - 1 && patOk; ++m) {
                int baseOff = 4 + 4 * m;
                if (!isReflex[idxOf(baseOff)] || !isReflex[idxOf(baseOff + 1)] ||
                    isReflex[idxOf(baseOff + 2)] || isReflex[idxOf(baseOff + 3)]) {
                    patOk = false;
                }
            }
            if (!patOk) continue;

            // Collect valley floor points ordered right-to-left (R[t], L[t] for t = 0..Tcnt-2)
            // and split parameters along bottom curve eBot = origSegs[kBot] (ordered ascending left-to-right).
            const CurveSeg<T>& eBot = origSegs[idxOf(0)];
            const CurveSeg<T>& eR   = origSegs[idxOf(1)];
            const CurveSeg<T>& eL   = origSegs[idxOf(K - 1)];

            std::vector<gsPoint2D<T>> rPts(Tcnt - 1), lPts(Tcnt - 1);
            std::vector<T> rParams(Tcnt - 1), lParams(Tcnt - 1);
            for (int t = 0; t < Tcnt - 1; ++t) {
                rPts[t] = verts[idxOf(4 + 4 * t)];
                lPts[t] = verts[idxOf(5 + 4 * t)];
                rParams[t] = eBot.closestParam(rPts[t]);
                lParams[t] = eBot.closestParam(lPts[t]);
            }

            // Build ordered valley vertices left-to-right: L_{T-2}, R_{T-2}, ..., L_0, R_0
            std::vector<gsPoint2D<T>> topNodes;
            std::vector<T> botSplitU;
            for (int t = Tcnt - 2; t >= 0; --t) {
                topNodes.push_back(lPts[t]);
                botSplitU.push_back(lParams[t]);
                topNodes.push_back(rPts[t]);
                botSplitU.push_back(rParams[t]);
            }

            bool monoU = !botSplitU.empty() && botSplitU.front() > T(0.02) && botSplitU.back() < T(0.98);
            for (size_t i = 1; i < botSplitU.size() && monoU; ++i) {
                if (botSplitU[i] <= botSplitU[i - 1] + T(0.02)) monoU = false;
            }
            if (!monoU) continue;

            auto splitMulti = [&](const CurveSeg<T>& crv, const std::vector<T>& us) -> std::vector<CurveSeg<T>> {
                std::vector<CurveSeg<T>> segs;
                CurveSeg<T> rem = crv;
                T prevU = T(0);
                for (T u : us) {
                    T localU = (u - prevU) / std::max(T(1e-9), T(1) - prevU);
                    auto sp = rem.splitAtParam(localU);
                    segs.push_back(sp.first);
                    rem = sp.second;
                    prevU = u;
                }
                segs.push_back(rem);
                return segs;
            };

            // (a) Diagonal-end comb decomposition: 3*Tcnt - 3 patches (no side splits)
            {
                std::vector<T> innerU;
                for (size_t i = 1; i + 1 < botSplitU.size(); ++i)
                    innerU.push_back(botSplitU[i]);
                auto botSegs9 = splitMulti(eBot, innerU); // size = 2*Tcnt - 3

                bool diagOk = true;
                gsPoint2D<T> vBotL = verts[idxOf(0)];
                gsPoint2D<T> vBotR = verts[idxOf(1)];
                if (!isSegmentInsidePoly(vBotL, topNodes.front(), outerPoly) ||
                    !isSegmentInsidePoly(topNodes.back(), vBotR, outerPoly)) {
                    diagOk = false;
                }
                for (size_t i = 1; i + 1 < topNodes.size() && diagOk; ++i) {
                    if (!isSegmentInsidePoly(botSegs9[i - 1].p1(), topNodes[i], outerPoly))
                        diagOk = false;
                }

                if (diagOk) {
                    gsMultiPatch<T> mp;
                    CurveSeg<T> diagL = makeStraightSeg(vBotL, topNodes.front(), m_degree);
                    CurveSeg<T> diagR = makeStraightSeg(topNodes.back(), vBotR, m_degree);

                    // Leftmost tooth (t = Tcnt - 1) and rightmost tooth (t = 0)
                    mp.addPatch(typename gsGeometry<T>::uPtr(
                        makeCoonsFromCCWCycle(diagL, origSegs[idxOf(K - 3)], origSegs[idxOf(K - 2)], eL, m_degree)));
                    mp.addPatch(typename gsGeometry<T>::uPtr(
                        makeCoonsFromCCWCycle(diagR, eR, origSegs[idxOf(2)], origSegs[idxOf(3)], m_degree)));

                    // Tooth bases left-to-right between topNodes[2*j-1] (= R_{...}) and topNodes[2*j] (= L_{...})
                    int nMidTop = static_cast<int>(topNodes.size()) - 1; // 2*Tcnt - 3
                    std::vector<CurveSeg<T>> spineTop(nMidTop);
                    for (int m = 0; m < nMidTop; ++m) {
                        if (m % 2 == 0) {
                            // Valley floor right-to-left index: gapIdx = (Tcnt - 2) - m/2
                            int gapIdx = (Tcnt - 2) - m / 2;
                            spineTop[m] = origSegs[idxOf(4 + 4 * gapIdx)].reversed();
                        } else {
                            // Middle tooth base across inner tooth t = (Tcnt - 1) - (m + 1)/2
                            int toothIdx = (Tcnt - 1) - (m + 1) / 2;
                            spineTop[m] = makeStraightSeg(topNodes[m], topNodes[m + 1], m_degree);
                            int eDropR = idxOf(3 + 4 * toothIdx);
                            int eTop   = idxOf(2 + 4 * toothIdx);
                            int eRiseL = idxOf(5 + 4 * (toothIdx - 1));
                            mp.addPatch(typename gsGeometry<T>::uPtr(
                                makeCoonsFromCCWCycle(spineTop[m], origSegs[eRiseL], origSegs[eTop], origSegs[eDropR], m_degree)));
                        }
                    }

                    // Spine quads left-to-right (m = 0 .. nMidTop - 1)
                    std::vector<CurveSeg<T>> pillars(nMidTop + 1);
                    pillars[0] = diagL.reversed();
                    for (int m = 1; m < nMidTop; ++m)
                        pillars[m] = makeStraightSeg(botSegs9[m - 1].p1(), topNodes[m], m_degree);
                    pillars[nMidTop] = diagR;

                    for (int m = 0; m < nMidTop; ++m) {
                        mp.addPatch(typename gsGeometry<T>::uPtr(
                            makeCoonsFromCCWCycle(botSegs9[m], pillars[m + 1], spineTop[m].reversed(), pillars[m].reversed(), m_degree)));
                    }

                    SplittingCandidate cand;
                    if (tryFinalize(mp, "CombDiagonal", cand)) {
                        registerIfBetter(cand);
                    }
                }
            }

            // (b) Orthogonal-spine comb decomposition: 3*Tcnt - 1 patches
            {
                T uL = eL.closestParam(topNodes.front());
                T uR = eR.closestParam(topNodes.back());
                if (uL > T(0.05) && uL < T(0.95) && uR > T(0.05) && uR < T(0.95)) {
                    auto spL = eL.splitAtParam(uL); // spL.first: v(K-1)->midL, spL.second: midL->v(0)
                    auto spR = eR.splitAtParam(uR); // spR.first: v(1)->midR, spR.second: midR->v(2)
                    auto botSegs11 = splitMulti(eBot, botSplitU); // size = 2*Tcnt - 1

                    gsPoint2D<T> midL = spL.first.p1();
                    gsPoint2D<T> midR = spR.first.p1();

                    std::vector<gsPoint2D<T>> allTopNodes;
                    allTopNodes.push_back(midL);
                    allTopNodes.insert(allTopNodes.end(), topNodes.begin(), topNodes.end());
                    allTopNodes.push_back(midR);

                    int nCols = 2 * Tcnt - 1;
                    bool gridOk = isSegmentInsidePoly(midL, topNodes.front(), outerPoly) &&
                                  isSegmentInsidePoly(topNodes.back(), midR, outerPoly);
                    for (int m = 1; m < nCols && gridOk; ++m) {
                        if (!isSegmentInsidePoly(botSegs11[m - 1].p1(), allTopNodes[m], outerPoly))
                            gridOk = false;
                    }

                    if (gridOk) {
                        gsMultiPatch<T> mp;
                        std::vector<CurveSeg<T>> spineTop(nCols);
                        for (int m = 0; m < nCols; ++m) {
                            if (m % 2 == 1) {
                                int gapIdx = (Tcnt - 2) - (m - 1) / 2;
                                spineTop[m] = origSegs[idxOf(4 + 4 * gapIdx)].reversed();
                            } else {
                                spineTop[m] = makeStraightSeg(allTopNodes[m], allTopNodes[m + 1], m_degree);
                                int toothIdx = (Tcnt - 1) - m / 2;
                                if (toothIdx == Tcnt - 1) {
                                    mp.addPatch(typename gsGeometry<T>::uPtr(
                                        makeCoonsFromCCWCycle(spineTop[m], origSegs[idxOf(K - 3)], origSegs[idxOf(K - 2)], spL.first, m_degree)));
                                } else if (toothIdx == 0) {
                                    mp.addPatch(typename gsGeometry<T>::uPtr(
                                        makeCoonsFromCCWCycle(spineTop[m], spR.second, origSegs[idxOf(2)], origSegs[idxOf(3)], m_degree)));
                                } else {
                                    int eDropR = idxOf(3 + 4 * toothIdx);
                                    int eTop   = idxOf(2 + 4 * toothIdx);
                                    int eRiseL = idxOf(5 + 4 * (toothIdx - 1));
                                    mp.addPatch(typename gsGeometry<T>::uPtr(
                                        makeCoonsFromCCWCycle(spineTop[m], origSegs[eRiseL], origSegs[eTop], origSegs[eDropR], m_degree)));
                                }
                            }
                        }

                        std::vector<CurveSeg<T>> pillars(nCols + 1);
                        pillars[0] = spL.second.reversed();
                        for (int m = 1; m < nCols; ++m)
                            pillars[m] = makeStraightSeg(botSegs11[m - 1].p1(), allTopNodes[m], m_degree);
                        pillars[nCols] = spR.first;

                        for (int m = 0; m < nCols; ++m) {
                            mp.addPatch(typename gsGeometry<T>::uPtr(
                                makeCoonsFromCCWCycle(botSegs11[m], pillars[m + 1], spineTop[m].reversed(), pillars[m].reversed(), m_degree)));
                        }

                        SplittingCandidate cand;
                        if (tryFinalize(mp, "CombOrthogonal", cand)) {
                            registerIfBetter(cand);
                        }
                    }
                }
            }
            break;
        }
    }

    if (!bestCand.isValid || nSharp <= 2 || bestCand.numPatches > 10) {
        const T domMinAng = getDomainMinCornerAngleDeg();
        const T minAngTol = std::min(static_cast<T>(m_options.getReal("MinCornerAngleDeg")), T(0.85) * domMinAng);
        const T maxAngTol = static_cast<T>(m_options.getReal("MaxCornerAngleDeg"));
        if (!bestCand.isValid || nSharp <= 2) {
            gsMultiPatch<T> offsetMp;
            if (internal_splitter::buildSmoothOffsetOGrid(m_domain, m_degree, offsetMp, minAngTol, std::max(maxAngTol, T(164.9)))) {
                SplittingCandidate candO;
                if (tryFinalize(offsetMp, "SmoothOffsetOGrid", candO))
                    registerIfBetter(candO);
            }
        }
        if (!bestCand.isValid || bestCand.numPatches > 10) {
            gsMultiPatch<T> recMp;
            std::string recLabel;
            if (internal_splitter::computeRecursiveCellSplitting(
                    m_domain, m_degree, m_weightPatches, m_weightDistortion,
                    minAngTol, std::max(maxAngTol, T(164.9)), recMp, recLabel)) {
                SplittingCandidate candR;
                if (tryFinalize(recMp, recLabel, candR))
                    registerIfBetter(candR);
            }
        }
    }

    return bestCand;
}

template <class T>
typename gsPatchSplitter<T>::SplittingCandidate
gsPatchSplitter<T>::computeHoleFreeTriangulationSplitting()
{
    using namespace internal_splitter;

    SplittingCandidate bestCand;
    bestCand.methodUsed = "HarmonicCrossField";
    bestCand.isValid = false;

    const auto& outer = m_domain.outer();
    const int K = outer.size();
    if (K == 0) return bestCand;

    std::vector<gsPoint2D<T>> outerPoly = sampleCurveLoopPoly(outer, 24);
    gsPoint2D<T> centroid = polygonCentroid(outerPoly);

    auto candCost = [&](const SplittingCandidate& c) -> T {
        return m_weightPatches * T(c.numPatches) + m_weightDistortion * c.distortionEnergy;
    };
    auto updateBest = [&](const SplittingCandidate& c) {
        if (!c.isValid) return;
        if (!bestCand.isValid || candCost(c) < candCost(bestCand))
            bestCand = c;
    };

    auto tryFinalize = [&](gsMultiPatch<T>& mp, SplittingCandidate& out) -> bool {
        if (mp.nPatches() == 0) return false;
        optimizeWinslowBarrier(mp);
        mp.computeTopology(T(1e-5), false, false);
        const T domMinAng = getDomainMinCornerAngleDeg();
        const T minAngTol = std::min(static_cast<T>(m_options.getReal("MinCornerAngleDeg")), T(0.85) * domMinAng);
        const T maxAngTol = static_cast<T>(m_options.getReal("MaxCornerAngleDeg"));
        if (!validateMultiPatch(m_domain, mp, T(1e-5), nullptr, minAngTol, maxAngTol))
            return false;
        T minDet = 0, minA = 0, maxA = 180;
        T energy = evaluateDistortionEnergy(mp, &minDet, &minA, &maxA);
        out.multiPatch = mp;
        out.numPatches = static_cast<int>(mp.nPatches());
        out.minJacobianDet = minDet;
        out.distortionEnergy = energy;
        out.minCornerAngleDeg = minA;
        out.maxCornerAngleDeg = maxA;
        out.methodUsed = "HarmonicCrossField";
        out.isValid = true;
        return true;
    };

    struct Tri { int v[3]; int e[3]; };  // edges: v0->v1, v1->v2, v2->v0
    struct Quad { int v[4]; int e[4]; }; // edges: v0->v1, v1->v2, v2->v3, v3->v0
    struct EdgeRec {
        int a, b;
        CurveSeg<T> seg; // oriented a -> b
        CurveSeg<T> h0, h1;
        gsPoint2D<T> mid;
        bool splitDone = false;
        bool isBoundary = false;
    };
    struct MixedMesh {
        std::vector<gsPoint2D<T>> pts;
        std::vector<bool> isSmooth;
        std::vector<EdgeRec> edges;
        std::vector<Tri> tris;
        std::vector<Quad> quads;
    };

    auto ensureSplit = [&](MixedMesh& mesh, int eIdx) {
        auto& er = mesh.edges[eIdx];
        if (!er.splitDone) {
            auto sp = er.seg.splitArc(T(0.5));
            er.h0 = sp.first;
            er.h1 = sp.second;
            er.mid = er.h0.p1();
            er.splitDone = true;
        }
    };

    auto getSegDir = [&](const MixedMesh& mesh, int eIdx, int u, int v) -> CurveSeg<T> {
        const auto& er = mesh.edges[eIdx];
        return (er.a == u && er.b == v) ? er.seg : er.seg.reversed();
    };

    auto getHalvesDir = [&](const MixedMesh& mesh, int eIdx, int u, int v,
                            CurveSeg<T>& out0, CurveSeg<T>& out1, gsPoint2D<T>& outMid) {
        const auto& er = mesh.edges[eIdx];
        outMid = er.mid;
        if (er.a == u && er.b == v) {
            out0 = er.h0;
            out1 = er.h1;
        } else {
            out0 = er.h1.reversed();
            out1 = er.h0.reversed();
        }
    };

    // Shared quadifier for any conforming Tri+Quad mixed mesh given edge-split flags.
    // Every triangle must have all 3 edges split; every quad must have opposite pairs matched.
    auto emitQuadPatches = [&](MixedMesh& mesh,
                               const std::vector<bool>& needSplit,
                               SplittingCandidate& candOut) -> bool
    {
        for (size_t i = 0; i < mesh.edges.size(); ++i) {
            if (needSplit[i]) ensureSplit(mesh, static_cast<int>(i));
        }

        gsMultiPatch<T> mp;

        // 1. Emit unmatched triangles (1-to-3 Catmull-Clark split)
        for (const auto& t : mesh.tris) {
            if (!needSplit[t.e[0]] || !needSplit[t.e[1]] || !needSplit[t.e[2]]) return false;
            int v0 = t.v[0], v1 = t.v[1], v2 = t.v[2];
            CurveSeg<T> e01_a, e01_b, e12_a, e12_b, e20_a, e20_b;
            gsPoint2D<T> m01, m12, m20;
            getHalvesDir(mesh, t.e[0], v0, v1, e01_a, e01_b, m01);
            getHalvesDir(mesh, t.e[1], v1, v2, e12_a, e12_b, m12);
            getHalvesDir(mesh, t.e[2], v2, v0, e20_a, e20_b, m20);

            gsPoint2D<T> O;
            if (!internal_splitter::pickOptimalTri1to3Center(mesh.pts[v0], mesh.pts[v1], mesh.pts[v2],
                                                             m01, m12, m20, O))
                return false;

            auto s01 = makeStraightSeg(O, m01, m_degree);
            auto s12 = makeStraightSeg(O, m12, m_degree);
            auto s20 = makeStraightSeg(O, m20, m_degree);
            mp.addPatch(typename gsGeometry<T>::uPtr(
                makeCoonsFromCCWCycle(e01_a, s01.reversed(), s20, e20_b, m_degree)));
            mp.addPatch(typename gsGeometry<T>::uPtr(
                makeCoonsFromCCWCycle(e12_a, s12.reversed(), s01, e01_b, m_degree)));
            mp.addPatch(typename gsGeometry<T>::uPtr(
                makeCoonsFromCCWCycle(e20_a, s20.reversed(), s12, e12_b, m_degree)));
        }

        // 2. Emit quads (0-split, 2-split along one opposite pair, or 4-split)
        for (const auto& q : mesh.quads) {
            bool sp0 = needSplit[q.e[0]], sp1 = needSplit[q.e[1]];
            bool sp2 = needSplit[q.e[2]], sp3 = needSplit[q.e[3]];
            if (sp0 != sp2 || sp1 != sp3) return false;
            int v0 = q.v[0], v1 = q.v[1], v2 = q.v[2], v3 = q.v[3];

            auto hasAnyBoundaryEdge = [&](int ea, int eb, int ec, int ed) -> bool {
                return mesh.edges[ea].isBoundary || mesh.edges[eb].isBoundary ||
                       mesh.edges[ec].isBoundary || mesh.edges[ed].isBoundary;
            };
            bool anyBnd = hasAnyBoundaryEdge(q.e[0], q.e[1], q.e[2], q.e[3]);

            if (!sp0 && !sp1) {
                gsPoint2D<T> c[4] = {mesh.pts[v0], mesh.pts[v1], mesh.pts[v2], mesh.pts[v3]};
                if (!anyBnd && !isStrictlyConvexQuad(c)) return false;
                mp.addPatch(typename gsGeometry<T>::uPtr(
                    makeCoonsFromCCWCycle(getSegDir(mesh, q.e[0], v0, v1),
                                          getSegDir(mesh, q.e[1], v1, v2),
                                          getSegDir(mesh, q.e[2], v2, v3),
                                          getSegDir(mesh, q.e[3], v3, v0),
                                          m_degree)));
            } else if (sp0 && !sp1) {
                CurveSeg<T> e0a, e0b, e2a, e2b;
                gsPoint2D<T> m01, m23;
                getHalvesDir(mesh, q.e[0], v0, v1, e0a, e0b, m01);
                getHalvesDir(mesh, q.e[2], v2, v3, e2a, e2b, m23);
                gsPoint2D<T> qa[4] = {mesh.pts[v0], m01, m23, mesh.pts[v3]};
                gsPoint2D<T> qb[4] = {m01, mesh.pts[v1], mesh.pts[v2], m23};
                if (!anyBnd && (!isStrictlyConvexQuad(qa) || !isStrictlyConvexQuad(qb))) return false;
                auto sMid = makeStraightSeg(m01, m23, m_degree);
                mp.addPatch(typename gsGeometry<T>::uPtr(
                    makeCoonsFromCCWCycle(e0a, sMid, e2b, getSegDir(mesh, q.e[3], v3, v0), m_degree)));
                mp.addPatch(typename gsGeometry<T>::uPtr(
                    makeCoonsFromCCWCycle(e0b, getSegDir(mesh, q.e[1], v1, v2), e2a, sMid.reversed(), m_degree)));
            } else if (!sp0 && sp1) {
                CurveSeg<T> e1a, e1b, e3a, e3b;
                gsPoint2D<T> m12, m30;
                getHalvesDir(mesh, q.e[1], v1, v2, e1a, e1b, m12);
                getHalvesDir(mesh, q.e[3], v3, v0, e3a, e3b, m30);
                gsPoint2D<T> qa[4] = {mesh.pts[v0], mesh.pts[v1], m12, m30};
                gsPoint2D<T> qb[4] = {m30, m12, mesh.pts[v2], mesh.pts[v3]};
                if (!anyBnd && (!isStrictlyConvexQuad(qa) || !isStrictlyConvexQuad(qb))) return false;
                auto sMid = makeStraightSeg(m12, m30, m_degree);
                mp.addPatch(typename gsGeometry<T>::uPtr(
                    makeCoonsFromCCWCycle(getSegDir(mesh, q.e[0], v0, v1), e1a, sMid, e3b, m_degree)));
                mp.addPatch(typename gsGeometry<T>::uPtr(
                    makeCoonsFromCCWCycle(sMid.reversed(), e1b, getSegDir(mesh, q.e[2], v2, v3), e3a, m_degree)));
            } else {
                CurveSeg<T> e0a, e0b, e1a, e1b, e2a, e2b, e3a, e3b;
                gsPoint2D<T> m01, m12, m23, m30;
                getHalvesDir(mesh, q.e[0], v0, v1, e0a, e0b, m01);
                getHalvesDir(mesh, q.e[1], v1, v2, e1a, e1b, m12);
                getHalvesDir(mesh, q.e[2], v2, v3, e2a, e2b, m23);
                getHalvesDir(mesh, q.e[3], v3, v0, e3a, e3b, m30);

                gsPoint2D<T> O_mid((m01.x + m12.x + m23.x + m30.x) / T(4),
                                   (m01.y + m12.y + m23.y + m30.y) / T(4));
                gsPoint2D<T> O_v((mesh.pts[v0].x + mesh.pts[v1].x + mesh.pts[v2].x + mesh.pts[v3].x) / T(4),
                                 (mesh.pts[v0].y + mesh.pts[v1].y + mesh.pts[v2].y + mesh.pts[v3].y) / T(4));
                auto check4 = [&](const gsPoint2D<T>& Oc) -> bool {
                    gsPoint2D<T> q0[4] = {mesh.pts[v0], m01, Oc, m30};
                    gsPoint2D<T> q1[4] = {mesh.pts[v1], m12, Oc, m01};
                    gsPoint2D<T> q2[4] = {mesh.pts[v2], m23, Oc, m12};
                    gsPoint2D<T> q3[4] = {mesh.pts[v3], m30, Oc, m23};
                    return isStrictlyConvexQuad(q0) && isStrictlyConvexQuad(q1)
                        && isStrictlyConvexQuad(q2) && isStrictlyConvexQuad(q3);
                };
                gsPoint2D<T> O = O_mid;
                if (!check4(O)) {
                    if (check4(O_v)) {
                        O = O_v;
                    } else {
                        bool found = false;
                        for (T w = T(0.1); w <= T(0.91); w += T(0.1)) {
                            gsPoint2D<T> Otry((1 - w) * O_mid.x + w * O_v.x,
                                              (1 - w) * O_mid.y + w * O_v.y);
                            if (check4(Otry)) { O = Otry; found = true; break; }
                        }
                        if (!found) return false;
                    }
                }
                auto s01 = makeStraightSeg(O, m01, m_degree);
                auto s12 = makeStraightSeg(O, m12, m_degree);
                auto s23 = makeStraightSeg(O, m23, m_degree);
                auto s30 = makeStraightSeg(O, m30, m_degree);
                mp.addPatch(typename gsGeometry<T>::uPtr(
                    makeCoonsFromCCWCycle(e0a, s01.reversed(), s30, e3b, m_degree)));
                mp.addPatch(typename gsGeometry<T>::uPtr(
                    makeCoonsFromCCWCycle(e1a, s12.reversed(), s01, e0b, m_degree)));
                mp.addPatch(typename gsGeometry<T>::uPtr(
                    makeCoonsFromCCWCycle(e2a, s23.reversed(), s12, e1b, m_degree)));
                mp.addPatch(typename gsGeometry<T>::uPtr(
                    makeCoonsFromCCWCycle(e3a, s30.reversed(), s23, e2b, m_degree)));
            }
        }
        return tryFinalize(mp, candOut);
    };

    // Opposite-edge parity propagation across quads starting from unmatched triangles.
    auto propagateOppositeSplits = [&](const MixedMesh& mesh, std::vector<bool>& needSplit) {
        needSplit.assign(mesh.edges.size(), false);
        for (const auto& t : mesh.tris) {
            needSplit[t.e[0]] = true;
            needSplit[t.e[1]] = true;
            needSplit[t.e[2]] = true;
        }
        bool changed = true;
        while (changed) {
            changed = false;
            for (const auto& q : mesh.quads) {
                if (needSplit[q.e[0]] != needSplit[q.e[2]]) {
                    needSplit[q.e[0]] = needSplit[q.e[2]] = true;
                    changed = true;
                }
                if (needSplit[q.e[1]] != needSplit[q.e[3]]) {
                    needSplit[q.e[1]] = needSplit[q.e[3]] = true;
                    changed = true;
                }
            }
        }
    };

    // Ensure no smooth (G1 ~180 deg) boundary vertex has triangle degree == 1 by bisecting
    // the opposite interior edge (Rivara 1-edge bisection preserves conforming mesh).
    auto fixSmoothDegreeOneVertices = [&](MixedMesh& mesh) {
        for (int iter = 0; iter < 8; ++iter) {
            std::vector<int> vDeg(mesh.pts.size(), 0);
            for (const auto& t : mesh.tris) {
                ++vDeg[t.v[0]]; ++vDeg[t.v[1]]; ++vDeg[t.v[2]];
            }
            int badTri = -1, badCorner = -1;
            for (int ti = 0; ti < static_cast<int>(mesh.tris.size()) && badTri < 0; ++ti) {
                for (int k = 0; k < 3; ++k) {
                    int vi = mesh.tris[ti].v[k];
                    if (vDeg[vi] != 1) continue;
                    int vn = mesh.tris[ti].v[(k + 1) % 3];
                    int vp = mesh.tris[ti].v[(k + 2) % 3];
                    T dx1 = mesh.pts[vn].x - mesh.pts[vi].x, dy1 = mesh.pts[vn].y - mesh.pts[vi].y;
                    T dx2 = mesh.pts[vp].x - mesh.pts[vi].x, dy2 = mesh.pts[vp].y - mesh.pts[vi].y;
                    T l1 = std::hypot(dx1, dy1), l2 = std::hypot(dx2, dy2);
                    T cosA = (l1 > T(1e-14) && l2 > T(1e-14)) ? (dx1 * dx2 + dy1 * dy2) / (l1 * l2) : T(1);
                    bool obtuseCorner = (cosA < T(-0.766)); // > 140 deg
                    if (mesh.isSmooth[vi] || obtuseCorner) {
                        int oppEdge = mesh.tris[ti].e[(k + 1) % 3];
                        if (!mesh.edges[oppEdge].isBoundary) {
                            badTri = ti;
                            badCorner = k;
                            break;
                        }
                    }
                }
            }
            if (badTri < 0) break;

            Tri t0 = mesh.tris[badTri];
            int v = t0.v[badCorner];
            int n = t0.v[(badCorner + 1) % 3];
            int p = t0.v[(badCorner + 2) % 3];
            int e_vn = t0.e[badCorner];
            int e_np = t0.e[(badCorner + 1) % 3];
            int e_pv = t0.e[(badCorner + 2) % 3];

            // Find neighbor triangle across e_np
            int mateTri = -1, mateEdgePos = -1;
            for (int ti = 0; ti < static_cast<int>(mesh.tris.size()); ++ti) {
                if (ti == badTri) continue;
                for (int k = 0; k < 3; ++k) {
                    if (mesh.tris[ti].e[k] == e_np) {
                        mateTri = ti;
                        mateEdgePos = k;
                        break;
                    }
                }
                if (mateTri >= 0) break;
            }
            if (mateTri < 0) break;

            Tri t1 = mesh.tris[mateTri];
            int w = t1.v[(mateEdgePos + 2) % 3];
            int e_pw = t1.e[(mateEdgePos + 1) % 3]; // edge from n->w or p->w depending on winding
            int e_wn = t1.e[(mateEdgePos + 2) % 3];
            // Wait: in t1, e[mateEdgePos] goes p -> n (opposite to t0's n -> p),
            // so t1.v[mateEdgePos] == p, t1.v[(mateEdgePos+1)%3] == n, w == t1.v[(mateEdgePos+2)%3].
            if (t1.v[mateEdgePos] != p || t1.v[(mateEdgePos + 1) % 3] != n) {
                // Adjust edge names to match (p, n, w) ordering
                e_pw = t1.e[(mateEdgePos + 1) % 3];
                e_wn = t1.e[(mateEdgePos + 2) % 3];
            } else {
                e_pw = t1.e[(mateEdgePos + 2) % 3];
                e_wn = t1.e[(mateEdgePos + 1) % 3];
            }

            ensureSplit(mesh, e_np);
            int mIdx = static_cast<int>(mesh.pts.size());
            mesh.pts.push_back(mesh.edges[e_np].mid);
            mesh.isSmooth.push_back(false);

            CurveSeg<T> s_nm, s_mp;
            gsPoint2D<T> midTmp;
            getHalvesDir(mesh, e_np, n, p, s_nm, s_mp, midTmp);

            auto addEdge = [&](int a, int b, const CurveSeg<T>& s) -> int {
                int id = static_cast<int>(mesh.edges.size());
                mesh.edges.push_back({a, b, s, CurveSeg<T>(), CurveSeg<T>(), gsPoint2D<T>(), false, false});
                return id;
            };
            int e_nm = addEdge(n, mIdx, s_nm);
            int e_mp = addEdge(mIdx, p, s_mp);
            int e_vm = addEdge(v, mIdx, makeStraightSeg(mesh.pts[v], mesh.pts[mIdx], m_degree));
            int e_wm = addEdge(w, mIdx, makeStraightSeg(mesh.pts[w], mesh.pts[mIdx], m_degree));

            // Replace badTri and mateTri cleanly
            mesh.tris[badTri]  = Tri{{v, n, mIdx}, {e_vn, e_nm, e_vm}};
            mesh.tris[mateTri] = Tri{{p, mIdx, w}, {e_mp, e_wm, e_pw}};
            mesh.tris.push_back(Tri{{p, v, mIdx}, {e_pv, e_vm, e_mp}});
            mesh.tris.push_back(Tri{{mIdx, n, w}, {e_nm, e_wn, e_wm}});
        }
    };

    // Lawson Delaunay edge flips on interior straight edges inside outerPoly.
    auto applyLawsonDelaunayFlips = [&](MixedMesh& mesh) {
        int nT = static_cast<int>(mesh.tris.size());
        int nE = static_cast<int>(mesh.edges.size());
        if (nT < 2 || nE == 0) return;
        for (int pass = 0; pass < nT * 2; ++pass) {
            struct EdgeAdj { int t0 = -1, k0 = -1, t1 = -1, k1 = -1; };
            std::vector<EdgeAdj> eAdj(nE);
            std::vector<int> vertDeg(mesh.pts.size(), 0);
            for (int ti = 0; ti < nT; ++ti) {
                for (int k = 0; k < 3; ++k) {
                    ++vertDeg[mesh.tris[ti].v[k]];
                    int eid = mesh.tris[ti].e[k];
                    if (eid >= 0 && eid < nE) {
                        if (eAdj[eid].t0 < 0) { eAdj[eid].t0 = ti; eAdj[eid].k0 = k; }
                        else if (eAdj[eid].t1 < 0) { eAdj[eid].t1 = ti; eAdj[eid].k1 = k; }
                    }
                }
            }
            bool flipped = false;
            for (int eIdx = 0; eIdx < nE; ++eIdx) {
                if (mesh.edges[eIdx].isBoundary) continue;
                int tA = eAdj[eIdx].t0, kA = eAdj[eIdx].k0;
                int tB = eAdj[eIdx].t1, kB = eAdj[eIdx].k1;
                if (tA < 0 || tB < 0) continue;
                Tri A = mesh.tris[tA], B = mesh.tris[tB];
                int a = A.v[kA], b = A.v[(kA + 1) % 3], c = A.v[(kA + 2) % 3];
                int d = B.v[(kB + 2) % 3];
                gsPoint2D<T> qa[4] = {mesh.pts[c], mesh.pts[a], mesh.pts[d], mesh.pts[b]};
                if (!isStrictlyConvexQuad(qa)) continue;
                if (!isSegmentInsidePoly(mesh.pts[c], mesh.pts[d], outerPoly)) continue;

                T curMin0, curMax0, curMin1, curMax1;
                internal_splitter::triMinMaxAnglesDeg(mesh.pts[a], mesh.pts[b], mesh.pts[c], curMin0, curMax0);
                internal_splitter::triMinMaxAnglesDeg(mesh.pts[b], mesh.pts[a], mesh.pts[d], curMin1, curMax1);
                T newMin0, newMax0, newMin1, newMax1;
                internal_splitter::triMinMaxAnglesDeg(mesh.pts[c], mesh.pts[a], mesh.pts[d], newMin0, newMax0);
                internal_splitter::triMinMaxAnglesDeg(mesh.pts[c], mesh.pts[d], mesh.pts[b], newMin1, newMax1);
                T curMin = std::min(curMin0, curMin1), newMin = std::min(newMin0, newMin1);
                T curMax = std::max(curMax0, curMax1), newMax = std::max(newMax0, newMax1);
                if (newMin <= curMin + T(1e-3) && !(std::abs(newMin - curMin) <= T(1e-3) && newMax < curMax - T(1e-2)))
                    continue;

                // Protect smooth boundary vertices from dropping to deg == 1 unless needed
                if ((mesh.isSmooth[a] && vertDeg[a] <= 2) || (mesh.isSmooth[b] && vertDeg[b] <= 2)) continue;

                // Find edge indices around A and B matching c->a, a->d, d->b, b->c
                int e_bc = A.e[(kA + 1) % 3];
                int e_ca = A.e[(kA + 2) % 3];
                int e_ad = -1, e_db = -1;
                for (int k = 0; k < 3; ++k) {
                    int u = B.v[k], v = B.v[(k + 1) % 3];
                    if ((u == a && v == d) || (u == d && v == a)) e_ad = B.e[k];
                    if ((u == d && v == b) || (u == b && v == d)) e_db = B.e[k];
                }
                if (e_ad < 0 || e_db < 0) continue;

                mesh.edges[eIdx] = {c, d, makeStraightSeg(mesh.pts[c], mesh.pts[d], m_degree),
                                    CurveSeg<T>(), CurveSeg<T>(), gsPoint2D<T>(), false, false};
                mesh.tris[tA] = Tri{{c, a, d}, {e_ca, e_ad, eIdx}};
                mesh.tris[tB] = Tri{{c, d, b}, {eIdx, e_db, e_bc}};
                flipped = true;
                break;
            }
            if (!flipped) break;
        }
    };

    // Given a valid conforming Tri mesh (quads empty), test:
    //   1) direct 1-to-3 quadification
    //   2) triangle joining (bounded backtracking + greedy) -> OddStripFinish AND CatmullClarkMixedFinish
    auto evaluateTriangleJoining = [&](const MixedMesh& triMesh) {
        // Direct 1-to-3 tri split baseline
        {
            MixedMesh m1 = triMesh;
            std::vector<bool> allTrue(m1.edges.size(), true);
            SplittingCandidate c1;
            if (emitQuadPatches(m1, allTrue, c1)) updateBest(c1);
        }

        const int nT = static_cast<int>(triMesh.tris.size());
        if (nT < 2) return;

        const T radToDeg = T(180) / T(3.14159265358979323846);
        const T domMinA = getDomainMinCornerAngleDeg();
        const T qMinTol = std::min(static_cast<T>(m_options.getReal("MinCornerAngleDeg")), T(0.85) * domMinA);
        const T qMaxTol = static_cast<T>(m_options.getReal("MaxCornerAngleDeg"));

        struct PairEdge {
            int t0, t1;
            Quad q;
            T score;
        };

        struct EdgeAdj { int t0 = -1, k0 = -1, t1 = -1, k1 = -1; };
        std::vector<EdgeAdj> eAdj(triMesh.edges.size());
        for (int ti = 0; ti < nT; ++ti) {
            for (int k = 0; k < 3; ++k) {
                int eid = triMesh.tris[ti].e[k];
                if (eid >= 0 && eid < static_cast<int>(triMesh.edges.size())) {
                    if (eAdj[eid].t0 < 0) { eAdj[eid].t0 = ti; eAdj[eid].k0 = k; }
                    else if (eAdj[eid].t1 < 0) { eAdj[eid].t1 = ti; eAdj[eid].k1 = k; }
                }
            }
        }

        std::vector<PairEdge> pairEdges;
        for (int eIdx = 0; eIdx < static_cast<int>(triMesh.edges.size()); ++eIdx) {
            if (triMesh.edges[eIdx].isBoundary) continue;
            int tA = eAdj[eIdx].t0, kA = eAdj[eIdx].k0;
            int tB = eAdj[eIdx].t1, kB = eAdj[eIdx].k1;
            if (tA < 0 || tB < 0) continue;
            Tri A = triMesh.tris[tA], B = triMesh.tris[tB];
            int a = A.v[kA], b = A.v[(kA + 1) % 3], c = A.v[(kA + 2) % 3];
            int d = B.v[(kB + 2) % 3];
            gsPoint2D<T> qc[4] = {triMesh.pts[c], triMesh.pts[a], triMesh.pts[d], triMesh.pts[b]};
            T minArea2 = 0;
            if (!isStrictlyConvexQuad(qc, &minArea2)) continue;

            int e_bc = A.e[(kA + 1) % 3];
            int e_ca = A.e[(kA + 2) % 3];
            int e_ad = -1, e_db = -1;
            for (int k = 0; k < 3; ++k) {
                int u = B.v[k], v = B.v[(k + 1) % 3];
                if ((u == a && v == d) || (u == d && v == a)) e_ad = B.e[k];
                if ((u == d && v == b) || (u == b && v == d)) e_db = B.e[k];
            }
            if (e_ad < 0 || e_db < 0) continue;

            int qv[4] = {c, a, d, b};
            int qe[4] = {e_ca, e_ad, e_db, e_bc};
            auto edgeDir = [&](int ei, int fromV, int toV, bool atEnd) -> gsPoint2D<T> {
                const auto& er = triMesh.edges[ei];
                if (er.isBoundary && er.seg.bsp) {
                    bool fwd = (er.a == fromV);
                    T u = (atEnd == fwd) ? T(1) : T(0);
                    gsPoint2D<T> dv = evalCurveDeriv(*er.seg.bsp, u);
                    T sgn = ((fwd && !atEnd) || (!fwd && atEnd)) ? T(1) : T(-1);
                    T nm = std::hypot(dv.x, dv.y);
                    if (nm > T(1e-14)) return gsPoint2D<T>(sgn * dv.x / nm, sgn * dv.y / nm);
                }
                T chSign = atEnd ? T(-1) : T(1);
                T dx = chSign * (triMesh.pts[toV].x - triMesh.pts[fromV].x);
                T dy = chSign * (triMesh.pts[toV].y - triMesh.pts[fromV].y);
                T nm = std::hypot(dx, dy);
                return (nm > T(1e-14)) ? gsPoint2D<T>(dx / nm, dy / nm) : gsPoint2D<T>(1, 0);
            };
            bool angOk = true;
            T minAng = T(360), maxAng = T(0), maxDev90 = T(0);
            T lenMin = T(1e30), lenMax = T(0);
            for (int k = 0; k < 4; ++k) {
                int kn = (k + 1) % 4, kp = (k + 3) % 4;
                gsPoint2D<T> tOut = edgeDir(qe[k],  qv[k], qv[kn], false);
                gsPoint2D<T> tBack = edgeDir(qe[kp], qv[kp], qv[k], true);
                T cr = tOut.x * tBack.y - tOut.y * tBack.x;
                T dt = tOut.x * tBack.x + tOut.y * tBack.y;
                T ang = std::atan2(cr, dt) * radToDeg;
                if (ang <= T(0)) ang += T(360);
                if (ang < qMinTol - T(1.5) || ang > qMaxTol) { angOk = false; break; }
                minAng = std::min(minAng, ang);
                maxAng = std::max(maxAng, ang);
                maxDev90 = std::max(maxDev90, std::abs(ang - T(90)));
                T eLen = triMesh.pts[qv[k]].dist(triMesh.pts[qv[kn]]);
                lenMin = std::min(lenMin, eLen);
                lenMax = std::max(lenMax, eLen);
            }
            if (!angOk) continue;

            T effMinQ = std::min(minAng, T(180) - maxAng);
            T arRatio = (lenMax > T(1e-14)) ? (lenMin / lenMax) : T(1);
            T sc = effMinQ + T(2.25) * arRatio - T(0.15) * maxDev90;
            pairEdges.push_back({tA, tB, Quad{{c, a, d, b}, {e_ca, e_ad, e_db, e_bc}}, sc});
        }
        if (pairEdges.empty()) return;

        std::sort(pairEdges.begin(), pairEdges.end(),
                  [](const PairEdge& p, const PairEdge& q) { return p.score > q.score; });

        std::vector<int> baseDeg(triMesh.pts.size(), 0);
        for (const auto& t : triMesh.tris) {
            ++baseDeg[t.v[0]]; ++baseDeg[t.v[1]]; ++baseDeg[t.v[2]];
        }

        // Collect distinct valid matchings (greedy + bounded backtracking)
        std::vector<std::vector<int>> candidateMatchings;
        auto pushUnique = [&](const std::vector<int>& m) {
            std::vector<int> ms = m;
            std::sort(ms.begin(), ms.end());
            for (const auto& prev : candidateMatchings)
                if (prev == ms) return;
            candidateMatchings.push_back(ms);
        };

        // 1. Greedy matching
        {
            std::vector<bool> used(nT, false);
            std::vector<int> deg = baseDeg;
            std::vector<int> gMatch;
            for (int pi = 0; pi < static_cast<int>(pairEdges.size()); ++pi) {
                const auto& pe = pairEdges[pi];
                if (used[pe.t0] || used[pe.t1]) continue;
                int a = pe.q.v[1], b = pe.q.v[3];
                if ((triMesh.isSmooth[a] && deg[a] <= 2) || (triMesh.isSmooth[b] && deg[b] <= 2))
                    continue;
                used[pe.t0] = used[pe.t1] = true;
                --deg[a]; --deg[b];
                gMatch.push_back(pi);
            }
            if (!gMatch.empty()) pushUnique(gMatch);
        }

        // 2. Bounded backtracking matching up to nT <= 26 to maximize pairs and quality
        if (nT <= 26) {
            std::vector<std::vector<int>> triPairs(nT);
            for (int pi = 0; pi < static_cast<int>(pairEdges.size()); ++pi) {
                triPairs[pairEdges[pi].t0].push_back(pi);
                triPairs[pairEdges[pi].t1].push_back(pi);
            }
            std::vector<bool> used(nT, false);
            std::vector<int> deg = baseDeg;
            std::vector<int> cur, bestMaxCard, bestMaxScore;
            int maxCard = 0;
            T maxScoreVal = -T(1);
            int nodes = 0;

            std::function<void(int, T)> dfs = [&](int ti, T curScore) {
                if (++nodes > 4000) return;
                while (ti < nT && used[ti]) ++ti;
                if (ti == nT) {
                    if (static_cast<int>(cur.size()) > maxCard) {
                        maxCard = static_cast<int>(cur.size());
                        bestMaxCard = cur;
                    }
                    T totalVal = T(cur.size()) * T(10) + curScore;
                    if (totalVal > maxScoreVal) {
                        maxScoreVal = totalVal;
                        bestMaxScore = cur;
                    }
                    return;
                }
                int remTris = 0;
                for (int j = ti; j < nT; ++j) if (!used[j]) ++remTris;
                if (static_cast<int>(cur.size()) + remTris / 2 < maxCard) return;

                // Try matching ti with each partner
                for (int pi : triPairs[ti]) {
                    const auto& pe = pairEdges[pi];
                    int other = (pe.t0 == ti) ? pe.t1 : pe.t0;
                    if (used[other]) continue;
                    int a = pe.q.v[1], b = pe.q.v[3];
                    if ((triMesh.isSmooth[a] && deg[a] <= 2) || (triMesh.isSmooth[b] && deg[b] <= 2))
                        continue;
                    used[ti] = used[other] = true;
                    --deg[a]; --deg[b];
                    cur.push_back(pi);
                    dfs(ti + 1, curScore + pe.score);
                    cur.pop_back();
                    ++deg[a]; ++deg[b];
                    used[ti] = used[other] = false;
                    if (nodes > 4000) return;
                }
                // Option: leave ti unmatched
                used[ti] = true;
                dfs(ti + 1, curScore);
                used[ti] = false;
            };
            dfs(0, T(0));
            if (!bestMaxCard.empty()) pushUnique(bestMaxCard);
            if (!bestMaxScore.empty()) pushUnique(bestMaxScore);
        }

        // Evaluate OddStripFinish and CatmullClarkMixedFinish on each matching
        for (const auto& m : candidateMatchings) {
            MixedMesh joined;
            joined.pts = triMesh.pts;
            joined.isSmooth = triMesh.isSmooth;
            joined.edges = triMesh.edges;
            std::vector<bool> usedTri(nT, false);
            for (int pi : m) {
                usedTri[pairEdges[pi].t0] = true;
                usedTri[pairEdges[pi].t1] = true;
                joined.quads.push_back(pairEdges[pi].q);
            }
            for (int ti = 0; ti < nT; ++ti) {
                if (!usedTri[ti]) joined.tris.push_back(triMesh.tris[ti]);
            }

            // Finish A: Opposite-edge parity propagation (minimal quad split)
            {
                MixedMesh jA = joined;
                std::vector<bool> needSplit;
                propagateOppositeSplits(jA, needSplit);
                SplittingCandidate cA;
                if (emitQuadPatches(jA, needSplit, cA)) updateBest(cA);
            }
            // Finish B: Uniform 1-edge Catmull-Clark mixed split fallback
            {
                MixedMesh jB = joined;
                std::vector<bool> allTrue(jB.edges.size(), true);
                SplittingCandidate cB;
                if (emitQuadPatches(jB, allTrue, cB)) updateBest(cB);
            }
        }
    };

    // Build initial curved boundary segs (subdividing K<=2 loops to >= 4 segments)
    std::vector<CurveSeg<T>> baseSegs;
    if (K == 1) {
        auto sH = CurveSeg<T>(outer.curve(0), m_degree).splitArc(T(0.5));
        auto s1 = sH.first.splitArc(T(0.5));
        auto s2 = sH.second.splitArc(T(0.5));
        baseSegs = {s1.first, s1.second, s2.first, s2.second};
    } else if (K == 2) {
        auto s0 = CurveSeg<T>(outer.curve(0), m_degree).splitArc(T(0.5));
        auto s1 = CurveSeg<T>(outer.curve(1), m_degree).splitArc(T(0.5));
        baseSegs = {s0.first, s0.second, s1.first, s1.second};
    } else {
        baseSegs.reserve(K);
        for (int i = 0; i < K; ++i) baseSegs.emplace_back(outer.curve(i), m_degree);
    }

    auto computeSmoothFlags = [&](const std::vector<CurveSeg<T>>& segs) -> std::vector<bool> {
        int M = static_cast<int>(segs.size());
        std::vector<bool> flags(M, false);
        const T radToDeg = T(180) / T(3.14159265358979323846);
        for (int i = 0; i < M; ++i) {
            const auto& prevBsp = *segs[(i + M - 1) % M].bsp;
            const auto& nextBsp = *segs[i].bsp;
            auto din = evalCurveDeriv(prevBsp, T(1));
            auto dout = evalCurveDeriv(nextBsp, T(0));
            T n1 = std::hypot(din.x, din.y), n2 = std::hypot(dout.x, dout.y);
            T cr = (n1 > 1e-14 && n2 > 1e-14) ? (din.x * dout.y - din.y * dout.x) / (n1 * n2) : T(0);
            T dt = (n1 > 1e-14 && n2 > 1e-14) ? (din.x * dout.x + din.y * dout.y) / (n1 * n2) : T(1);
            T turnRad = std::atan2(cr, dt);
            T intAngDeg = T(180) - turnRad * radToDeg;
            // Mark smooth (~180 deg) OR obtuse (> 140 deg) vertices so Rivara bisection splits them
            flags[i] = (std::abs(turnRad) <= T(0.10)) || (intAngDeg > T(140));
        }
        return flags;
    };

    // Curved-safe ear clipping supporting sharp-first protection and multiple scoring modes
    auto buildEarTriangulation = [&](int mode, MixedMesh& outMesh) -> bool {
        struct RingNode {
            int ptIdx;
            CurveSeg<T> outSeg; // curve along active boundary leaving this node
        };
        std::vector<RingNode> ring;
        std::vector<gsPoint2D<T>> pts;
        std::vector<bool> isSmooth = computeSmoothFlags(baseSegs);
        std::vector<EdgeRec> edges;
        std::map<std::pair<int, int>, int> edgeMap;

        for (size_t i = 0; i < baseSegs.size(); ++i) {
            int idx = static_cast<int>(i);
            pts.push_back(baseSegs[i].p0());
            ring.push_back({idx, baseSegs[i]});
        }
        for (size_t i = 0; i < baseSegs.size(); ++i) {
            int a = static_cast<int>(i);
            int b = static_cast<int>((i + 1) % baseSegs.size());
            int id = static_cast<int>(edges.size());
            edges.push_back({a, b, baseSegs[i], CurveSeg<T>(), CurveSeg<T>(), gsPoint2D<T>(), false, true});
            edgeMap[{std::min(a, b), std::max(a, b)}] = id;
        }

        auto getOrAddDiag = [&](int a, int b) -> int {
            int lo = std::min(a, b), hi = std::max(a, b);
            auto it = edgeMap.find({lo, hi});
            if (it != edgeMap.end()) return it->second;
            int id = static_cast<int>(edges.size());
            edges.push_back({lo, hi, makeStraightSeg(pts[lo], pts[hi], m_degree),
                             CurveSeg<T>(), CurveSeg<T>(), gsPoint2D<T>(), false, false});
            edgeMap[{lo, hi}] = id;
            return id;
        };

        std::vector<int> triCount(pts.size(), 0);
        std::vector<Tri> tris;
        int subDivs = 0;

        auto pointInTri = [&](const gsPoint2D<T>& A, const gsPoint2D<T>& B,
                              const gsPoint2D<T>& C, const gsPoint2D<T>& P) -> bool {
            return cross2(A, B, P) >= -1e-9 && cross2(B, C, P) >= -1e-9 && cross2(C, A, P) >= -1e-9;
        };

        while (ring.size() > 3) {
            int m = static_cast<int>(ring.size());
            int bestIdx = -1;
            T bestScore = -T(1e30);

            // Two passes: pass 0 avoids leaving any smooth vertex with final deg == 1
            for (int pass = 0; pass < 2 && bestIdx < 0; ++pass) {
                for (int idx = 0; idx < m; ++idx) {
                    int iPrev = (idx + m - 1) % m;
                    int iNext = (idx + 1) % m;
                    int ia = ring[iPrev].ptIdx;
                    int ib = ring[idx].ptIdx;
                    int ic = ring[iNext].ptIdx;

                    if (pass == 0 && isSmooth[ib] && triCount[ib] == 0) continue;

                    const auto& A = pts[ia];
                    const auto& B = pts[ib];
                    const auto& C = pts[ic];
                    T area2 = cross2(A, B, C);
                    if (area2 <= T(1e-8)) continue;

                    if (!isSegmentInsidePoly(A, C, outerPoly)) continue;

                    bool clear = true;
                    for (int k = 0; k < m && clear; ++k) {
                        int ik = ring[k].ptIdx;
                        if (ik == ia || ik == ib || ik == ic) continue;
                        if (pointInTri(A, B, C, pts[ik])) clear = false;
                    }
                    // Also check interior samples of other active ring curves against curved bulge
                    for (int k = 0; k < m && clear; ++k) {
                        if (k == iPrev || k == idx) continue;
                        const auto& bsp = *ring[k].outSeg.bsp;
                        gsMatrix<T> ev;
                        gsMatrix<T> par(1, 3);
                        par(0, 0) = T(0.25); par(0, 1) = T(0.5); par(0, 2) = T(0.75);
                        bsp.eval_into(par, ev);
                        for (int s = 0; s < 3; ++s) {
                            if (pointInTri(A, B, C, gsPoint2D<T>(ev(0, s), ev(1, s)))) {
                                clear = false;
                                break;
                            }
                        }
                    }
                    if (!clear) continue;

                    T sc = T(0);
                    T dAC2 = std::pow(A.x - C.x, 2) + std::pow(A.y - C.y, 2);
                    if (mode == 0) {
                        T l2 = std::pow(A.x - B.x, 2) + std::pow(A.y - B.y, 2)
                             + std::pow(B.x - C.x, 2) + std::pow(B.y - C.y, 2) + dAC2;
                        sc = area2 / (l2 + T(1e-14));
                    } else if (mode == 1) {
                        sc = T(1) / (dAC2 + T(1e-14));
                    } else {
                        T lAB = std::hypot(A.x - B.x, A.y - B.y);
                        T lBC = std::hypot(B.x - C.x, B.y - C.y);
                        sc = area2 / ((lAB * lBC) * (dAC2 + T(1e-14)) + T(1e-14));
                    }
                    if (isSmooth[ib]) sc *= T(0.2);
                    if (sc > bestScore) {
                        bestScore = sc;
                        bestIdx = idx;
                    }
                }
            }

            if (bestIdx < 0) {
                // Subdivide active ring curve with largest chord sagitta (or longest curved edge)
                if (subDivs >= 24) return false;
                int splitK = -1;
                T maxSag = T(0);
                int longestK = -1;
                T maxLen = T(0);
                for (int k = 0; k < m; ++k) {
                    int kNext = (k + 1) % m;
                    auto sp = ring[k].outSeg.splitArc(T(0.5));
                    gsPoint2D<T> mid = sp.first.p1();
                    const auto& P = pts[ring[k].ptIdx];
                    const auto& Q = pts[ring[kNext].ptIdx];
                    T chord = std::hypot(Q.x - P.x, Q.y - P.y);
                    if (chord > maxLen) { maxLen = chord; longestK = k; }
                    T sag = (chord > 1e-12) ? std::abs(cross2(P, Q, mid)) / (chord * chord + T(1e-12)) : T(0);
                    if (sag > maxSag) {
                        maxSag = sag;
                        splitK = k;
                    }
                }
                if (splitK < 0 || maxSag < T(1e-4)) splitK = longestK;
                if (splitK < 0) return false;

                int kNext = (splitK + 1) % m;
                int u = ring[splitK].ptIdx;
                int v = ring[kNext].ptIdx;
                int eOld = edgeMap[{std::min(u, v), std::max(u, v)}];
                auto sp = ring[splitK].outSeg.splitArc(T(0.5));
                gsPoint2D<T> midPt = sp.first.p1();

                int newIdx = static_cast<int>(pts.size());
                pts.push_back(midPt);
                isSmooth.push_back(true);
                triCount.push_back(0);

                // Keep eOld as u -> midPt, add midPt -> v
                bool oldWasBdr = edges[eOld].isBoundary;
                edges[eOld] = {u, newIdx, sp.first, CurveSeg<T>(), CurveSeg<T>(), gsPoint2D<T>(), false, oldWasBdr};
                edgeMap.erase({std::min(u, v), std::max(u, v)});
                edgeMap[{std::min(u, newIdx), std::max(u, newIdx)}] = eOld;
                int eSecond = static_cast<int>(edges.size());
                edges.push_back({newIdx, v, sp.second, CurveSeg<T>(), CurveSeg<T>(), gsPoint2D<T>(), false, oldWasBdr});
                edgeMap[{std::min(newIdx, v), std::max(newIdx, v)}] = eSecond;

                ring[splitK].outSeg = sp.first;
                ring.insert(ring.begin() + splitK + 1, {newIdx, sp.second});
                ++subDivs;
                continue;
            }

            int iPrev = (bestIdx + m - 1) % m;
            int iNext = (bestIdx + 1) % m;
            int ia = ring[iPrev].ptIdx;
            int ib = ring[bestIdx].ptIdx;
            int ic = ring[iNext].ptIdx;

            int e0 = edgeMap[{std::min(ia, ib), std::max(ia, ib)}];
            int e1 = edgeMap[{std::min(ib, ic), std::max(ib, ic)}];
            int e2 = getOrAddDiag(ia, ic);
            tris.push_back({{ia, ib, ic}, {e0, e1, e2}});
            ++triCount[ia];
            ++triCount[ib];
            ++triCount[ic];

            const auto& er2 = edges[e2];
            ring[iPrev].outSeg = (er2.a == ia && er2.b == ic) ? er2.seg : er2.seg.reversed();
            ring.erase(ring.begin() + bestIdx);
        }

        if (ring.size() != 3) return false;
        int ia = ring[0].ptIdx, ib = ring[1].ptIdx, ic = ring[2].ptIdx;
        if (cross2(pts[ia], pts[ib], pts[ic]) <= T(1e-8)) return false;
        int e0 = getOrAddDiag(ia, ib);
        int e1 = getOrAddDiag(ib, ic);
        int e2 = getOrAddDiag(ic, ia);
        tris.push_back({{ia, ib, ic}, {e0, e1, e2}});

        outMesh.pts = std::move(pts);
        outMesh.isSmooth = std::move(isSmooth);
        outMesh.edges = std::move(edges);
        outMesh.tris = std::move(tris);
        outMesh.quads.clear();
        return true;
    };

    // Special K==3 reflex triangle -> 2 triangles by cutting reflex corner to opposite midpoint
    if (K == 3) {
        std::vector<gsPoint2D<T>> polyPts = {baseSegs[0].p0(), baseSegs[1].p0(), baseSegs[2].p0()};
        for (int r = 0; r < 3; ++r) {
            int p = (r + 2) % 3, opp = (r + 1) % 3;
            if (cross2(polyPts[p], polyPts[r], polyPts[opp]) < -1e-8) {
                auto spOpp = baseSegs[opp].splitArc(T(0.5));
                MixedMesh mesh;
                mesh.pts = {polyPts[0], polyPts[1], polyPts[2], spOpp.first.p1()};
                mesh.isSmooth = computeSmoothFlags(baseSegs);
                mesh.isSmooth.push_back(false);
                auto addEdge = [&](int a, int b, const CurveSeg<T>& s, bool isBdr) -> int {
                    int id = static_cast<int>(mesh.edges.size());
                    mesh.edges.push_back({a, b, s, CurveSeg<T>(), CurveSeg<T>(), gsPoint2D<T>(), false, isBdr});
                    return id;
                };
                int e_r_opp = addEdge(r, opp, baseSegs[r], true);
                int e_opp_m = addEdge(opp, 3, spOpp.first, true);
                int e_m_r   = addEdge(3, r, makeStraightSeg(mesh.pts[3], mesh.pts[r], m_degree), false);
                int e_m_p   = addEdge(3, p, spOpp.second, true);
                int e_p_r   = addEdge(p, r, baseSegs[p], true);
                mesh.tris = {
                    {{r, opp, 3}, {e_r_opp, e_opp_m, e_m_r}},
                    {{r, 3, p},   {e_m_r,   e_m_p,   e_p_r}}
                };
                evaluateTriangleJoining(mesh);
                break;
            }
        }
    }

    // Strategy T1: General curved ear-clipping (3 priority orders) + Rivara smooth fix + Lawson flips + triangle joining
    for (int mode = 0; mode < 3; ++mode) {
        MixedMesh mesh;
        if (!buildEarTriangulation(mode, mesh)) continue;
        fixSmoothDegreeOneVertices(mesh);
        evaluateTriangleJoining(mesh);

        MixedMesh flipped = mesh;
        applyLawsonDelaunayFlips(flipped);
        evaluateTriangleJoining(flipped);
    }

    // Strategy T2: Centroid-Star Triangulation + triangle joining
    if (pointInPolyStrict(centroid, outerPoly, T(1e-5))) {
        int M = static_cast<int>(baseSegs.size());
        std::vector<gsPoint2D<T>> pts(M + 1);
        bool allVisible = true;
        for (int i = 0; i < M; ++i) {
            pts[i] = baseSegs[i].p0();
            if (!isSegmentInsidePoly(centroid, pts[i], outerPoly)) allVisible = false;
        }
        pts[M] = centroid;

        if (allVisible && M >= 3) {
            MixedMesh starMesh;
            starMesh.pts = pts;
            starMesh.isSmooth = computeSmoothFlags(baseSegs);
            starMesh.isSmooth.push_back(false);
            for (int i = 0; i < M; ++i) {
                starMesh.edges.push_back({i, (i + 1) % M, baseSegs[i],
                                          CurveSeg<T>(), CurveSeg<T>(), gsPoint2D<T>(), false, true});
            }
            for (int i = 0; i < M; ++i) {
                starMesh.edges.push_back({M, i, makeStraightSeg(centroid, pts[i], m_degree),
                                          CurveSeg<T>(), CurveSeg<T>(), gsPoint2D<T>(), false, false});
            }
            for (int i = 0; i < M; ++i) {
                int j = (i + 1) % M;
                starMesh.tris.push_back({{M, i, j}, {M + i, i, M + j}});
            }
            evaluateTriangleJoining(starMesh);
        }
    }

    return bestCand;
}

template <class T>
typename gsPatchSplitter<T>::SplittingCandidate
gsPatchSplitter<T>::computeMultiHoleTriangulationSplitting(bool coarseOnly)
{
    using namespace internal_splitter;

    SplittingCandidate bestCand;
    bestCand.methodUsed = "MultiHoleTriangulation";
    bestCand.isValid = false;

    if (m_domain.numLoops() == 0) return bestCand;
    const gsMatrix<T> bbox = m_domain.boundingBox();
    const T diag = std::hypot(bbox(0, 1) - bbox(0, 0), bbox(1, 1) - bbox(1, 0));
    const T domMinAng = getDomainMinCornerAngleDeg();
    const T minAngTol = std::min(static_cast<T>(m_options.getReal("MinCornerAngleDeg")), T(0.85) * domMinAng);
    const T maxAngTol = static_cast<T>(m_options.getReal("MaxCornerAngleDeg"));

    auto updateBest = [&](const SplittingCandidate& c) {
        if (!c.isValid) return;
        if (!bestCand.isValid) { bestCand = c; return; }
        T qNew = m_weightPatches * T(c.numPatches) + m_weightDistortion * c.distortionEnergy;
        T qOld = m_weightPatches * T(bestCand.numPatches) + m_weightDistortion * bestCand.distortionEnergy;
        if (qNew < qOld) bestCand = c;
    };

    int globalCoupledBudget = coarseOnly ? 3 : 16;
    std::function<void(TriMeshData<T>&, bool, bool)> quadifyTriMesh;
    quadifyTriMesh = [&](TriMeshData<T>& mesh, bool coarseRelaxed, bool allowParityRepair) {
        const int nT = static_cast<int>(mesh.tris.size());
        if (nT == 0) return;

        for (auto& er : mesh.edges) {
            if (!er.seg.bsp)
                er.seg = makeStraightSeg(mesh.pts[er.a], mesh.pts[er.b], m_degree);
            if (!er.splitDone) {
                auto sp = er.seg.splitArc(T(0.5));
                er.h0 = sp.first;
                er.h1 = sp.second;
                er.mid = sp.first.p1();
                er.splitDone = true;
                er.splitFrac = T(0.5);
            }
        }

        const T effMinValAng = coarseRelaxed ? std::min(minAngTol, std::max(T(4), minAngTol - T(3.5))) : minAngTol;
        const T effMaxValAng = coarseRelaxed ? std::max(maxAngTol, T(164.2))     : maxAngTol;
        const T effMinPairAng = std::min(minAngTol, std::max(T(4), minAngTol - T(4.0)));
        const T effMaxPairAng = std::max(maxAngTol, T(164.2));

        auto emitFromMixed = [&](const std::vector<typename TriMeshData<T>::Tri>& tris,
                                 const std::vector<typename TriMeshData<T>::Quad>& quads,
                                 const std::vector<bool>& needSplit,
                                 std::vector<int>* outFailedQuads = nullptr) -> bool {
            std::vector<typename TriMeshData<T>::EdgeRec> edges = mesh.edges;

            auto resplitEdge = [&](int eIdx, T paramU) {
                auto& er = edges[eIdx];
                T u = std::max(T(0.10), std::min(T(0.90), paramU));
                auto sp = er.seg.splitAtParam(u);
                er.h0 = sp.first;
                er.h1 = sp.second;
                er.mid = sp.first.p1();
                er.splitDone = true;
                er.splitFrac = u;
            };
            auto getSegDir = [&](int eIdx, int u, int v) -> CurveSeg<T> {
                const auto& er = edges[eIdx];
                return (er.a == u && er.b == v) ? er.seg : er.seg.reversed();
            };
            auto getHalvesDir = [&](int eIdx, int u, int v, CurveSeg<T>& out0, CurveSeg<T>& out1, gsPoint2D<T>& outMid) {
                const auto& er = edges[eIdx];
                outMid = er.mid;
                if (er.a == u && er.b == v) { out0 = er.h0; out1 = er.h1; }
                else                         { out0 = er.h1.reversed(); out1 = er.h0.reversed(); }
            };

            for (const auto& t : tris)
                if (!needSplit[t.e[0]] || !needSplit[t.e[1]] || !needSplit[t.e[2]]) return false;
            for (const auto& q : quads)
                if (needSplit[q.e[0]] != needSplit[q.e[2]] || needSplit[q.e[1]] != needSplit[q.e[3]]) return false;

            auto chooseQuad4Center = [&](const typename TriMeshData<T>::Quad& q,
                                         const CurveSeg<T>& e0a, const CurveSeg<T>& e0b, const gsPoint2D<T>& m01,
                                         const CurveSeg<T>& e1a, const CurveSeg<T>& e1b, const gsPoint2D<T>& m12,
                                         const CurveSeg<T>& e2a, const CurveSeg<T>& e2b, const gsPoint2D<T>& m23,
                                         const CurveSeg<T>& e3a, const CurveSeg<T>& e3b, const gsPoint2D<T>& m30,
                                         gsPoint2D<T>& outO, T& outEffMin) -> bool {
                int v0 = q.v[0], v1 = q.v[1], v2 = q.v[2], v3 = q.v[3];
                gsPoint2D<T> O_m((m01.x + m12.x + m23.x + m30.x) / T(4), (m01.y + m12.y + m23.y + m30.y) / T(4));
                gsPoint2D<T> O_v((mesh.pts[v0].x + mesh.pts[v1].x + mesh.pts[v2].x + mesh.pts[v3].x) / T(4),
                                 (mesh.pts[v0].y + mesh.pts[v1].y + mesh.pts[v2].y + mesh.pts[v3].y) / T(4));
                std::vector<gsPoint2D<T>> oCands = {O_m, O_v};
                for (T w : {T(0.25), T(0.5), T(0.75)})
                    oCands.emplace_back((T(1) - w) * O_m.x + w * O_v.x, (T(1) - w) * O_m.y + w * O_v.y);
                for (const auto& vp : {mesh.pts[v0], mesh.pts[v1], mesh.pts[v2], mesh.pts[v3]})
                    oCands.emplace_back(T(0.75) * O_m.x + T(0.25) * vp.x, T(0.75) * O_m.y + T(0.25) * vp.y);

                bool found = false;
                T bestEff = -T(1e9);
                for (const auto& Oc : oCands) {
                    gsPoint2D<T> q0[4] = {mesh.pts[v0], m01, Oc, m30};
                    gsPoint2D<T> q1[4] = {mesh.pts[v1], m12, Oc, m01};
                    gsPoint2D<T> q2[4] = {mesh.pts[v2], m23, Oc, m12};
                    gsPoint2D<T> q3[4] = {mesh.pts[v3], m30, Oc, m23};
                    if (!isStrictlyConvexQuad(q0) || !isStrictlyConvexQuad(q1) ||
                        !isStrictlyConvexQuad(q2) || !isStrictlyConvexQuad(q3))
                        continue;
                    auto s01 = makeStraightSeg(Oc, m01, m_degree);
                    auto s12 = makeStraightSeg(Oc, m12, m_degree);
                    auto s23 = makeStraightSeg(Oc, m23, m_degree);
                    auto s30 = makeStraightSeg(Oc, m30, m_degree);
                    T mnA, mxA, eff0, eff1, eff2, eff3;
                    eval4SegCornerAnglesDeg(e0a, s01.reversed(), s30, e3b, mnA, mxA, eff0);
                    eval4SegCornerAnglesDeg(e1a, s12.reversed(), s01, e0b, mnA, mxA, eff1);
                    eval4SegCornerAnglesDeg(e2a, s23.reversed(), s12, e1b, mnA, mxA, eff2);
                    eval4SegCornerAnglesDeg(e3a, s30.reversed(), s23, e2b, mnA, mxA, eff3);
                    T effMin = std::min(std::min(eff0, eff1), std::min(eff2, eff3));
                    if (!found || effMin > bestEff) {
                        bestEff = effMin;
                        outO = Oc;
                        outEffMin = effMin;
                        found = true;
                    }
                }
                return found;
            };

            auto evalFaceSplitQuality = [&](bool isTri, int fIdx, T& minAngOut, T& maxAngOut, T& effMinOut) -> bool {
                minAngOut = T(360); maxAngOut = T(0); effMinOut = -T(1);
                if (isTri) {
                    const auto& t = tris[fIdx];
                    int v0 = t.v[0], v1 = t.v[1], v2 = t.v[2];
                    CurveSeg<T> e01a, e01b, e12a, e12b, e20a, e20b;
                    gsPoint2D<T> m01, m12, m20, O;
                    getHalvesDir(t.e[0], v0, v1, e01a, e01b, m01);
                    getHalvesDir(t.e[1], v1, v2, e12a, e12b, m12);
                    getHalvesDir(t.e[2], v2, v0, e20a, e20b, m20);
                    if (!internal_splitter::pickOptimalTri1to3CenterCurved(
                            mesh.pts[v0], mesh.pts[v1], mesh.pts[v2],
                            e01a, e01b, e12a, e12b, e20a, e20b,
                            edges[t.e[0]].isBoundary, edges[t.e[1]].isBoundary, edges[t.e[2]].isBoundary, O))
                        return false;
                    auto s01 = makeStraightSeg(O, m01, m_degree);
                    auto s12 = makeStraightSeg(O, m12, m_degree);
                    auto s20 = makeStraightSeg(O, m20, m_degree);
                    T mn, mx, ef;
                    eval4SegCornerAnglesDeg(e01a, s01.reversed(), s20, e20b, mn, mx, ef);
                    minAngOut = std::min(minAngOut, mn); maxAngOut = std::max(maxAngOut, mx);
                    eval4SegCornerAnglesDeg(e12a, s12.reversed(), s01, e01b, mn, mx, ef);
                    minAngOut = std::min(minAngOut, mn); maxAngOut = std::max(maxAngOut, mx);
                    eval4SegCornerAnglesDeg(e20a, s20.reversed(), s12, e12b, mn, mx, ef);
                    minAngOut = std::min(minAngOut, mn); maxAngOut = std::max(maxAngOut, mx);
                    effMinOut = std::min(minAngOut, T(180) - maxAngOut);
                    return (minAngOut >= effMinValAng && maxAngOut <= effMaxValAng);
                }
                const auto& q = quads[fIdx];
                bool sp0 = needSplit[q.e[0]], sp1 = needSplit[q.e[1]];
                int v0 = q.v[0], v1 = q.v[1], v2 = q.v[2], v3 = q.v[3];
                if (!sp0 && !sp1) {
                    gsPoint2D<T> c[4] = {mesh.pts[v0], mesh.pts[v1], mesh.pts[v2], mesh.pts[v3]};
                    if (!isStrictlyConvexQuad(c)) return false;
                    eval4SegCornerAnglesDeg(getSegDir(q.e[0], v0, v1), getSegDir(q.e[1], v1, v2),
                                            getSegDir(q.e[2], v2, v3), getSegDir(q.e[3], v3, v0),
                                            minAngOut, maxAngOut, effMinOut);
                    return (minAngOut >= effMinValAng && maxAngOut <= effMaxValAng);
                } else if (sp0 && !sp1) {
                    CurveSeg<T> e0a, e0b, e2a, e2b; gsPoint2D<T> m01, m23;
                    getHalvesDir(q.e[0], v0, v1, e0a, e0b, m01);
                    getHalvesDir(q.e[2], v2, v3, e2a, e2b, m23);
                    gsPoint2D<T> qa[4] = {mesh.pts[v0], m01, m23, mesh.pts[v3]};
                    gsPoint2D<T> qb[4] = {m01, mesh.pts[v1], mesh.pts[v2], m23};
                    if (!isStrictlyConvexQuad(qa) || !isStrictlyConvexQuad(qb)) return false;
                    auto sMid = makeStraightSeg(m01, m23, m_degree);
                    T mn, mx, ef;
                    eval4SegCornerAnglesDeg(e0a, sMid, e2b, getSegDir(q.e[3], v3, v0), mn, mx, ef);
                    minAngOut = std::min(minAngOut, mn); maxAngOut = std::max(maxAngOut, mx);
                    eval4SegCornerAnglesDeg(e0b, getSegDir(q.e[1], v1, v2), e2a, sMid.reversed(), mn, mx, ef);
                    minAngOut = std::min(minAngOut, mn); maxAngOut = std::max(maxAngOut, mx);
                    effMinOut = std::min(minAngOut, T(180) - maxAngOut);
                    return (minAngOut >= effMinValAng && maxAngOut <= effMaxValAng);
                } else if (!sp0 && sp1) {
                    CurveSeg<T> e1a, e1b, e3a, e3b; gsPoint2D<T> m12, m30;
                    getHalvesDir(q.e[1], v1, v2, e1a, e1b, m12);
                    getHalvesDir(q.e[3], v3, v0, e3a, e3b, m30);
                    gsPoint2D<T> qa[4] = {mesh.pts[v0], mesh.pts[v1], m12, m30};
                    gsPoint2D<T> qb[4] = {m30, m12, mesh.pts[v2], mesh.pts[v3]};
                    if (!isStrictlyConvexQuad(qa) || !isStrictlyConvexQuad(qb)) return false;
                    auto sMid = makeStraightSeg(m12, m30, m_degree);
                    T mn, mx, ef;
                    eval4SegCornerAnglesDeg(getSegDir(q.e[0], v0, v1), e1a, sMid, e3b, mn, mx, ef);
                    minAngOut = std::min(minAngOut, mn); maxAngOut = std::max(maxAngOut, mx);
                    eval4SegCornerAnglesDeg(sMid.reversed(), e1b, getSegDir(q.e[2], v2, v3), e3a, mn, mx, ef);
                    minAngOut = std::min(minAngOut, mn); maxAngOut = std::max(maxAngOut, mx);
                    effMinOut = std::min(minAngOut, T(180) - maxAngOut);
                    return (minAngOut >= effMinValAng && maxAngOut <= effMaxValAng);
                } else {
                    CurveSeg<T> e0a, e0b, e1a, e1b, e2a, e2b, e3a, e3b; gsPoint2D<T> m01, m12, m23, m30, O;
                    getHalvesDir(q.e[0], v0, v1, e0a, e0b, m01);
                    getHalvesDir(q.e[1], v1, v2, e1a, e1b, m12);
                    getHalvesDir(q.e[2], v2, v3, e2a, e2b, m23);
                    getHalvesDir(q.e[3], v3, v0, e3a, e3b, m30);
                    T dummyEff = 0;
                    if (!chooseQuad4Center(q, e0a, e0b, m01, e1a, e1b, m12, e2a, e2b, m23, e3a, e3b, m30, O, dummyEff))
                        return false;
                    auto s01 = makeStraightSeg(O, m01, m_degree);
                    auto s12 = makeStraightSeg(O, m12, m_degree);
                    auto s23 = makeStraightSeg(O, m23, m_degree);
                    auto s30 = makeStraightSeg(O, m30, m_degree);
                    T mn, mx, ef;
                    eval4SegCornerAnglesDeg(e0a, s01.reversed(), s30, e3b, mn, mx, ef);
                    minAngOut = std::min(minAngOut, mn); maxAngOut = std::max(maxAngOut, mx);
                    eval4SegCornerAnglesDeg(e1a, s12.reversed(), s01, e0b, mn, mx, ef);
                    minAngOut = std::min(minAngOut, mn); maxAngOut = std::max(maxAngOut, mx);
                    eval4SegCornerAnglesDeg(e2a, s23.reversed(), s12, e1b, mn, mx, ef);
                    minAngOut = std::min(minAngOut, mn); maxAngOut = std::max(maxAngOut, mx);
                    eval4SegCornerAnglesDeg(e3a, s30.reversed(), s23, e2b, mn, mx, ef);
                    minAngOut = std::min(minAngOut, mn); maxAngOut = std::max(maxAngOut, mx);
                    effMinOut = std::min(minAngOut, T(180) - maxAngOut);
                    return (minAngOut >= effMinValAng && maxAngOut <= effMaxValAng);
                }
            };

            // Build Coons multi-patch from current `edges`, recording per-face/per-edge patch ownership.
            // When `useCurvedSpokes` is true, internal 1-to-3 / 1-to-2 / 1-to-4 spokes use symmetric
            // cubic Hermite curves departing interior centers symmetrically and arriving orthogonally at boundaries.
            auto tryBuildFromCurrentEdges = [&](std::vector<int>* outBadSplitEdges = nullptr,
                                                std::vector<int>* outBadIntEdges = nullptr,
                                                std::vector<int>* outBadQuads = nullptr,
                                                bool useCurvedSpokes = false) -> bool {
                gsMultiPatch<T> mp;
                std::vector<std::vector<int>> patchEdges;
                std::vector<int> patchQuadIdx;
                auto recordEdges = [&](const std::vector<int>& eList, int qIdx) {
                    if (outBadSplitEdges)
                        for (int eid : eList)
                            if (needSplit[eid]) outBadSplitEdges->push_back(eid);
                    if (outBadIntEdges)
                        for (int eid : eList)
                            if (!edges[eid].isBoundary && !needSplit[eid]) outBadIntEdges->push_back(eid);
                    if (outBadQuads && qIdx >= 0)
                        outBadQuads->push_back(qIdx);
                };
                auto inwardNormalAtMid = [&](int eIdx, const gsPoint2D<T>& mPt, const gsPoint2D<T>& refIn) -> gsPoint2D<T> {
                    gsPoint2D<T> toIn(refIn.x - mPt.x, refIn.y - mPt.y);
                    T ln = std::hypot(toIn.x, toIn.y);
                    gsPoint2D<T> fallback = (ln > T(1e-12)) ? gsPoint2D<T>(toIn.x / ln, toIn.y / ln) : gsPoint2D<T>(1, 0);
                    if (!edges[eIdx].isBoundary || !edges[eIdx].h0.bsp) return fallback;
                    gsPoint2D<T> td = evalCurveDeriv(*edges[eIdx].h0.bsp, T(1));
                    T tn = std::hypot(td.x, td.y);
                    if (tn <= T(1e-12)) return fallback;
                    gsPoint2D<T> n1(-td.y / tn, td.x / tn), n2(td.y / tn, -td.x / tn);
                    return (n1.x * fallback.x + n1.y * fallback.y >= n2.x * fallback.x + n2.y * fallback.y) ? n1 : n2;
                };
                auto buildSpoke = [&](const gsPoint2D<T>& O, const gsPoint2D<T>& mPt, int eIdx,
                                      const gsPoint2D<T>& tO) -> CurveSeg<T> {
                    if (!useCurvedSpokes) return makeStraightSeg(O, mPt, m_degree);
                    gsPoint2D<T> nM = inwardNormalAtMid(eIdx, mPt, O);
                    return internal_splitter::makeSymmetricHermiteSpoke(O, mPt, tO, nM, mesh.loopPolys, m_degree);
                };
                for (const auto& t : tris) {
                    int v0 = t.v[0], v1 = t.v[1], v2 = t.v[2];
                    CurveSeg<T> e01_a, e01_b, e12_a, e12_b, e20_a, e20_b;
                    gsPoint2D<T> m01, m12, m20;
                    getHalvesDir(t.e[0], v0, v1, e01_a, e01_b, m01);
                    getHalvesDir(t.e[1], v1, v2, e12_a, e12_b, m12);
                    getHalvesDir(t.e[2], v2, v0, e20_a, e20_b, m20);

                    gsPoint2D<T> O;
                    if (!internal_splitter::pickOptimalTri1to3CenterCurved(
                            mesh.pts[v0], mesh.pts[v1], mesh.pts[v2],
                            e01_a, e01_b, e12_a, e12_b, e20_a, e20_b,
                            edges[t.e[0]].isBoundary, edges[t.e[1]].isBoundary, edges[t.e[2]].isBoundary, O)) {
                        recordEdges({t.e[0], t.e[1], t.e[2]}, -1);
                        return false;
                    }
                    // Symmetric 120-deg star frame at O aligned to (m01-O, m12-O, m20-O)
                    const T twoPi3 = T(2.09439510239319549);
                    T a0 = std::atan2(m01.y - O.y, m01.x - O.x);
                    T a1 = std::atan2(m12.y - O.y, m12.x - O.x);
                    T a2 = std::atan2(m20.y - O.y, m20.x - O.x);
                    T sx = std::cos(a0) + std::cos(a1 - twoPi3) + std::cos(a2 - T(2) * twoPi3);
                    T sy = std::sin(a0) + std::sin(a1 - twoPi3) + std::sin(a2 - T(2) * twoPi3);
                    T baseA = (std::hypot(sx, sy) > T(1e-9)) ? std::atan2(sy, sx) : a0;
                    gsPoint2D<T> d0(std::cos(baseA), std::sin(baseA));
                    gsPoint2D<T> d1(std::cos(baseA + twoPi3), std::sin(baseA + twoPi3));
                    gsPoint2D<T> d2(std::cos(baseA + T(2) * twoPi3), std::sin(baseA + T(2) * twoPi3));

                    auto s01 = buildSpoke(O, m01, t.e[0], d0);
                    auto s12 = buildSpoke(O, m12, t.e[1], d1);
                    auto s20 = buildSpoke(O, m20, t.e[2], d2);
                    mp.addPatch(typename gsGeometry<T>::uPtr(makeCoonsFromCCWCycle(e01_a, s01.reversed(), s20, e20_b, m_degree)));
                    patchEdges.push_back({t.e[0], t.e[2]}); patchQuadIdx.push_back(-1);
                    mp.addPatch(typename gsGeometry<T>::uPtr(makeCoonsFromCCWCycle(e12_a, s12.reversed(), s01, e01_b, m_degree)));
                    patchEdges.push_back({t.e[1], t.e[0]}); patchQuadIdx.push_back(-1);
                    mp.addPatch(typename gsGeometry<T>::uPtr(makeCoonsFromCCWCycle(e20_a, s20.reversed(), s12, e12_b, m_degree)));
                    patchEdges.push_back({t.e[2], t.e[1]}); patchQuadIdx.push_back(-1);
                }
                for (int qi = 0; qi < static_cast<int>(quads.size()); ++qi) {
                    const auto& q = quads[qi];
                    bool sp0 = needSplit[q.e[0]], sp1 = needSplit[q.e[1]];
                    int v0 = q.v[0], v1 = q.v[1], v2 = q.v[2], v3 = q.v[3];
                    auto pushQBad = [&]() {
                        recordEdges({q.e[0], q.e[1], q.e[2], q.e[3]}, qi);
                    };
                    auto hasAnyBoundaryEdge = [&](int ea, int eb, int ec, int ed) -> bool {
                        return mesh.edges[ea].isBoundary || mesh.edges[eb].isBoundary ||
                               mesh.edges[ec].isBoundary || mesh.edges[ed].isBoundary;
                    };
                    bool anyBnd = hasAnyBoundaryEdge(q.e[0], q.e[1], q.e[2], q.e[3]);
                    if (!sp0 && !sp1) {
                        gsPoint2D<T> c[4] = {mesh.pts[v0], mesh.pts[v1], mesh.pts[v2], mesh.pts[v3]};
                        if (!anyBnd && !isStrictlyConvexQuad(c)) { pushQBad(); return false; }
                        mp.addPatch(typename gsGeometry<T>::uPtr(
                            makeCoonsFromCCWCycle(getSegDir(q.e[0], v0, v1), getSegDir(q.e[1], v1, v2),
                                                  getSegDir(q.e[2], v2, v3), getSegDir(q.e[3], v3, v0), m_degree)));
                        patchEdges.push_back({q.e[0], q.e[1], q.e[2], q.e[3]}); patchQuadIdx.push_back(qi);
                    } else if (sp0 && !sp1) {
                        CurveSeg<T> e0a, e0b, e2a, e2b; gsPoint2D<T> m01, m23;
                        getHalvesDir(q.e[0], v0, v1, e0a, e0b, m01);
                        getHalvesDir(q.e[2], v2, v3, e2a, e2b, m23);
                        gsPoint2D<T> qa[4] = {mesh.pts[v0], m01, m23, mesh.pts[v3]};
                        gsPoint2D<T> qb[4] = {m01, mesh.pts[v1], mesh.pts[v2], m23};
                        if (!anyBnd && (!isStrictlyConvexQuad(qa) || !isStrictlyConvexQuad(qb))) { pushQBad(); return false; }
                        auto sMid = useCurvedSpokes
                            ? internal_splitter::makeSymmetricHermiteSpoke(
                                  m01, m23, inwardNormalAtMid(q.e[0], m01, m23),
                                  inwardNormalAtMid(q.e[2], m23, m01), mesh.loopPolys, m_degree)
                            : makeStraightSeg(m01, m23, m_degree);
                        mp.addPatch(typename gsGeometry<T>::uPtr(makeCoonsFromCCWCycle(e0a, sMid, e2b, getSegDir(q.e[3], v3, v0), m_degree)));
                        patchEdges.push_back({q.e[0], q.e[2], q.e[3]}); patchQuadIdx.push_back(qi);
                        mp.addPatch(typename gsGeometry<T>::uPtr(makeCoonsFromCCWCycle(e0b, getSegDir(q.e[1], v1, v2), e2a, sMid.reversed(), m_degree)));
                        patchEdges.push_back({q.e[0], q.e[1], q.e[2]}); patchQuadIdx.push_back(qi);
                    } else if (!sp0 && sp1) {
                        CurveSeg<T> e1a, e1b, e3a, e3b; gsPoint2D<T> m12, m30;
                        getHalvesDir(q.e[1], v1, v2, e1a, e1b, m12);
                        getHalvesDir(q.e[3], v3, v0, e3a, e3b, m30);
                        gsPoint2D<T> qa[4] = {mesh.pts[v0], mesh.pts[v1], m12, m30};
                        gsPoint2D<T> qb[4] = {m30, m12, mesh.pts[v2], mesh.pts[v3]};
                        if (!anyBnd && (!isStrictlyConvexQuad(qa) || !isStrictlyConvexQuad(qb))) { pushQBad(); return false; }
                        auto sMid = useCurvedSpokes
                            ? internal_splitter::makeSymmetricHermiteSpoke(
                                  m12, m30, inwardNormalAtMid(q.e[1], m12, m30),
                                  inwardNormalAtMid(q.e[3], m30, m12), mesh.loopPolys, m_degree)
                            : makeStraightSeg(m12, m30, m_degree);
                        mp.addPatch(typename gsGeometry<T>::uPtr(makeCoonsFromCCWCycle(getSegDir(q.e[0], v0, v1), e1a, sMid, e3b, m_degree)));
                        patchEdges.push_back({q.e[0], q.e[1], q.e[3]}); patchQuadIdx.push_back(qi);
                        mp.addPatch(typename gsGeometry<T>::uPtr(makeCoonsFromCCWCycle(sMid.reversed(), e1b, getSegDir(q.e[2], v2, v3), e3a, m_degree)));
                        patchEdges.push_back({q.e[1], q.e[2], q.e[3]}); patchQuadIdx.push_back(qi);
                    } else {
                        CurveSeg<T> e0a, e0b, e1a, e1b, e2a, e2b, e3a, e3b; gsPoint2D<T> m01, m12, m23, m30;
                        getHalvesDir(q.e[0], v0, v1, e0a, e0b, m01);
                        getHalvesDir(q.e[1], v1, v2, e1a, e1b, m12);
                        getHalvesDir(q.e[2], v2, v3, e2a, e2b, m23);
                        getHalvesDir(q.e[3], v3, v0, e3a, e3b, m30);
                        gsPoint2D<T> O;
                        T dummyEff = 0;
                        if (!chooseQuad4Center(q, e0a, e0b, m01, e1a, e1b, m12, e2a, e2b, m23, e3a, e3b, m30, O, dummyEff)) {
                            pushQBad();
                            return false;
                        }
                        const T pi2 = T(1.57079632679489662);
                        T b0 = std::atan2(m01.y - O.y, m01.x - O.x);
                        T b1 = std::atan2(m12.y - O.y, m12.x - O.x);
                        T b2 = std::atan2(m23.y - O.y, m23.x - O.x);
                        T b3 = std::atan2(m30.y - O.y, m30.x - O.x);
                        T qx = std::cos(b0) + std::cos(b1 - pi2) + std::cos(b2 - T(2) * pi2) + std::cos(b3 - T(3) * pi2);
                        T qy = std::sin(b0) + std::sin(b1 - pi2) + std::sin(b2 - T(2) * pi2) + std::sin(b3 - T(3) * pi2);
                        T qA = (std::hypot(qx, qy) > T(1e-9)) ? std::atan2(qy, qx) : b0;
                        auto s01 = buildSpoke(O, m01, q.e[0], gsPoint2D<T>(std::cos(qA), std::sin(qA)));
                        auto s12 = buildSpoke(O, m12, q.e[1], gsPoint2D<T>(std::cos(qA + pi2), std::sin(qA + pi2)));
                        auto s23 = buildSpoke(O, m23, q.e[2], gsPoint2D<T>(std::cos(qA + T(2) * pi2), std::sin(qA + T(2) * pi2)));
                        auto s30 = buildSpoke(O, m30, q.e[3], gsPoint2D<T>(std::cos(qA + T(3) * pi2), std::sin(qA + T(3) * pi2)));
                        mp.addPatch(typename gsGeometry<T>::uPtr(makeCoonsFromCCWCycle(e0a, s01.reversed(), s30, e3b, m_degree)));
                        patchEdges.push_back({q.e[0], q.e[3]}); patchQuadIdx.push_back(qi);
                        mp.addPatch(typename gsGeometry<T>::uPtr(makeCoonsFromCCWCycle(e1a, s12.reversed(), s01, e0b, m_degree)));
                        patchEdges.push_back({q.e[1], q.e[0]}); patchQuadIdx.push_back(qi);
                        mp.addPatch(typename gsGeometry<T>::uPtr(makeCoonsFromCCWCycle(e2a, s23.reversed(), s12, e1b, m_degree)));
                        patchEdges.push_back({q.e[2], q.e[1]}); patchQuadIdx.push_back(qi);
                        mp.addPatch(typename gsGeometry<T>::uPtr(makeCoonsFromCCWCycle(e3a, s30.reversed(), s23, e2b, m_degree)));
                        patchEdges.push_back({q.e[3], q.e[2]}); patchQuadIdx.push_back(qi);
                    }
                }
                if (mp.nPatches() == 0) return false;
                bool angPreOk = true;
                int badAngCount = 0;
                for (size_t p = 0; p < mp.nPatches(); ++p) {
                    T angs[4], mnA = 0, mxA = 180;
                    bool degen = false;
                    evalPatchCornerAnglesDeg(mp.patch(p), angs, mnA, mxA, degen);
                    if (degen || mnA < effMinValAng || mxA > effMaxValAng) {
                        angPreOk = false;
                        ++badAngCount;
                        recordEdges(patchEdges[p], patchQuadIdx[p]);
                    }
                }
                // If coarse layout has only mild angle/fold issues, try coupled curved-interface relaxation before rejecting!
                bool triedCoupled = false;
                const int maxFixableBad = std::min<int>(4, std::max<int>(1, static_cast<int>(mp.nPatches()) / 6));
                if (!angPreOk && globalCoupledBudget > 0 && mp.nPatches() >= 2 && mp.nPatches() <= 80 && badAngCount <= maxFixableBad) {
                    --globalCoupledBudget;
                    gsMultiPatch<T> cMp(mp);
                    cMp.computeTopology(T(1e-5), false, false);
                    if (2 * cMp.nInterfaces() + cMp.nBoundary() == 4 * cMp.nPatches() &&
                        optimizeCoupledInterfacesAndInteriors(cMp, T(1e-7), 14, false)) {
                        bool cAngOk = true;
                        for (size_t p = 0; p < cMp.nPatches(); ++p) {
                            T angs[4], mnA = 0, mxA = 180; bool degen = false;
                            evalPatchCornerAnglesDeg(cMp.patch(p), angs, mnA, mxA, degen);
                            if (degen || mnA < effMinValAng || mxA > effMaxValAng) { cAngOk = false; break; }
                        }
                        if (cAngOk) {
                            mp = give(cMp);
                            angPreOk = true;
                            triedCoupled = true;
                        }
                    }
                }
                if (!angPreOk) return false;
                gsMultiPatch<T> preWinslow = mp;
                if (!triedCoupled) optimizeWinslowBarrier(mp);
                mp.computeTopology(T(1e-5), false, false);
                std::string valErr;
                bool valOk = validateMultiPatch(m_domain, mp, T(1e-4), &valErr, effMinValAng, effMaxValAng);
                if (!valOk) {
                    preWinslow.computeTopology(T(1e-5), false, false);
                    if (validateMultiPatch(m_domain, preWinslow, T(1e-4), &valErr, effMinValAng, effMaxValAng)) {
                        mp = give(preWinslow);
                        valOk = true;
                    }
                }
                if (!valOk && !triedCoupled && globalCoupledBudget > 0 && mp.nPatches() >= 2 && mp.nPatches() <= 80) {
                    --globalCoupledBudget;
                    gsMultiPatch<T> cMp(preWinslow);
                    cMp.computeTopology(T(1e-5), false, false);
                    if (2 * cMp.nInterfaces() + cMp.nBoundary() == 4 * cMp.nPatches() &&
                        optimizeCoupledInterfacesAndInteriors(cMp, T(1e-7), 14, false)) {
                        cMp.computeTopology(T(1e-5), false, false);
                        if (validateMultiPatch(m_domain, cMp, T(1e-4), &valErr, effMinValAng, effMaxValAng)) {
                            mp = give(cMp);
                            valOk = true;
                        }
                    }
                }
                if (valOk) {
                    gsASG1LocalFitter<T>::unifyPatchOrientations(mp);
                    mp.computeTopology(T(1e-5), false, false);
                    for (size_t p = 0; p < mp.nPatches(); ++p) {
                        T d49 = internal_splitter::evalPatchMinDetDense(mp.patch(p), 49);
                        if (d49 <= T(1e-8)) {
                            if (auto* tb = dynamic_cast<gsTensorBSpline<2, T>*>(&mp.patch(p))) {
                                gsMatrix<T> saveC = tb->coefs();
                                internal_splitter::quickHarmonicPatchRepair(*tb, T(1e-7));
                                d49 = internal_splitter::evalPatchMinDetDense(*tb, 49);
                                if (d49 <= T(1e-8)) tb->coefs() = saveC;
                            }
                        }
                        if (d49 <= T(1e-8)) {
                            recordEdges(patchEdges[p], patchQuadIdx[p]);
                            valOk = false;
                        }
                    }
                }
                if (!valOk) {
                    for (size_t p = 0; p < mp.nPatches(); ++p) {
                        T angs[4], mnA = 0, mxA = 180; bool degen = false;
                        evalPatchCornerAnglesDeg(mp.patch(p), angs, mnA, mxA, degen);
                        bool singleSign = isPatchDetSingleSign(mp.patch(p), 17, T(1e-7));
                        if (degen || mnA < effMinValAng || mxA > effMaxValAng || !singleSign) {
                            recordEdges(patchEdges[p], patchQuadIdx[p]);
                        }
                    }
                    return false;
                }
                T minDet = 0, minA = 0, maxA = 180;
                T energy = evaluateDistortionEnergy(mp, &minDet, &minA, &maxA);
                SplittingCandidate cand;
                cand.multiPatch = give(mp);
                cand.numPatches = static_cast<int>(cand.multiPatch.nPatches());
                cand.minJacobianDet = minDet;
                cand.distortionEnergy = energy;
                cand.minCornerAngleDeg = minA;
                cand.maxCornerAngleDeg = maxA;
                cand.methodUsed = "MultiHoleTriangulation";
                cand.isValid = true;
                updateBest(cand);
                return true;
            };

            // Step 1: Fast baseline with untouched straight interior chords and 0.5 midpoints
            std::vector<int> badEdges, badIntEdges, badQuads;
            if (tryBuildFromCurrentEdges(&badEdges, &badIntEdges, &badQuads)) return true;

            if (badEdges.empty() && badIntEdges.empty()) {
                if (outFailedQuads) {
                    std::sort(badQuads.begin(), badQuads.end());
                    badQuads.erase(std::unique(badQuads.begin(), badQuads.end()), badQuads.end());
                    *outFailedQuads = std::move(badQuads);
                }
                return false;
            }

            std::sort(badEdges.begin(), badEdges.end());
            badEdges.erase(std::unique(badEdges.begin(), badEdges.end()), badEdges.end());
            if (badEdges.size() > 24) badEdges.resize(24);
            std::sort(badIntEdges.begin(), badIntEdges.end());
            badIntEdges.erase(std::unique(badIntEdges.begin(), badIntEdges.end()), badIntEdges.end());
            if (badIntEdges.size() > 24) badIntEdges.resize(24);

            // Step 2: Adaptive curved-interior-chord + split-point rescue ONLY on edges bordering failing patches
            struct IncFace { bool isTri; int idx; };
            std::vector<std::vector<IncFace>> edgeFaces(edges.size());
            for (int ti = 0; ti < static_cast<int>(tris.size()); ++ti)
                for (int k = 0; k < 3; ++k)
                    edgeFaces[tris[ti].e[k]].push_back({true, ti});
            for (int qi = 0; qi < static_cast<int>(quads.size()); ++qi)
                for (int k = 0; k < 4; ++k)
                    edgeFaces[quads[qi].e[k]].push_back({false, qi});

            // Step 2a: Bend failing unsplit interior chords into orthogonal-arrival cubic Beziers
            auto optimizeCurvedInteriorEdge = [&](int eIdx) -> bool {
                const auto& inc = edgeFaces[eIdx];
                if (inc.empty() || edges[eIdx].isBoundary || needSplit[eIdx]) return false;
                int va = edges[eIdx].a, vb = edges[eIdx].b;
                const gsPoint2D<T>& A = mesh.pts[va];
                const gsPoint2D<T>& B = mesh.pts[vb];
                T L = A.dist(B);
                if (L <= T(1e-9)) return false;

                gsPoint2D<T> dAB((B.x - A.x) / L, (B.y - A.y) / L);
                gsPoint2D<T> perp(-dAB.y, dAB.x);

                auto inwardNormalAt = [&](int vIdx, const gsPoint2D<T>& toOther) -> gsPoint2D<T> {
                    if (!mesh.isBoundaryPt[vIdx]) return toOther;
                    gsPoint2D<T> tAvg(mesh.bndTanIn[vIdx].x + mesh.bndTanOut[vIdx].x,
                                      mesh.bndTanIn[vIdx].y + mesh.bndTanOut[vIdx].y);
                    T tn = std::hypot(tAvg.x, tAvg.y);
                    if (tn <= T(1e-12)) return toOther;
                    gsPoint2D<T> n1(-tAvg.y / tn, tAvg.x / tn);
                    gsPoint2D<T> n2(tAvg.y / tn, -tAvg.x / tn);
                    gsPoint2D<T> nUse = (n1.x * toOther.x + n1.y * toOther.y >= n2.x * toOther.x + n2.y * toOther.y) ? n1 : n2;
                    gsPoint2D<T> probe(mesh.pts[vIdx].x + T(0.05) * L * nUse.x,
                                       mesh.pts[vIdx].y + T(0.05) * L * nUse.y);
                    if (!isPointStrictlyInsideDomain(probe, mesh.loopPolys, T(1e-6)))
                        return toOther;
                    return nUse;
                };

                gsPoint2D<T> nA = inwardNormalAt(va, dAB);
                gsPoint2D<T> nB = inwardNormalAt(vb, gsPoint2D<T>(-dAB.x, -dAB.y));

                auto evalInc = [&](bool checkDet) -> std::pair<bool, T> {
                    T minAll = T(1e9);
                    bool okAll = true;
                    for (const auto& f : inc) {
                        T mnA, mxA, efA;
                        bool okA = evalFaceSplitQuality(f.isTri, f.idx, mnA, mxA, efA);
                        okAll = okAll && okA;
                        minAll = std::min(minAll, okA ? efA : (efA - T(200)));
                        if (checkDet && !f.isTri) {
                            const auto& q = quads[f.idx];
                            if (!needSplit[q.e[0]] && !needSplit[q.e[1]]) {
                                std::unique_ptr<gsTensorBSpline<2, T>> qp(
                                    makeCoonsFromCCWCycle(getSegDir(q.e[0], q.v[0], q.v[1]),
                                                          getSegDir(q.e[1], q.v[1], q.v[2]),
                                                          getSegDir(q.e[2], q.v[2], q.v[3]),
                                                          getSegDir(q.e[3], q.v[3], q.v[0]), m_degree));
                                if (!qp || !isPatchDetSingleSign(*qp, 11, T(1e-8))) {
                                    okAll = false;
                                    minAll = std::min(minAll, T(-300));
                                }
                            }
                        }
                    }
                    return {okAll, minAll};
                };

                auto baseRes = evalInc(true);
                CurveSeg<T> oldSeg = edges[eIdx].seg;
                CurveSeg<T> bestSeg = oldSeg;
                bool bestOk = baseRes.first;
                T bestScore = baseRes.second;
                bool improved = false;

                std::vector<CurveSeg<T>> segCands;
                const T weights[5][2] = {{T(0.45), T(0.45)}, {T(0.75), T(0.75)}, {T(0.65), T(0.0)}, {T(0.0), T(0.65)}, {T(0.95), T(0.95)}};
                for (int wi = 0; wi < 5; ++wi) {
                    gsPoint2D<T> tA((T(1) - weights[wi][0]) * dAB.x + weights[wi][0] * nA.x,
                                    (T(1) - weights[wi][0]) * dAB.y + weights[wi][0] * nA.y);
                    gsPoint2D<T> tB((T(1) - weights[wi][1]) * (-dAB.x) + weights[wi][1] * nB.x,
                                    (T(1) - weights[wi][1]) * (-dAB.y) + weights[wi][1] * nB.y);
                    T na = std::hypot(tA.x, tA.y), nb = std::hypot(tB.x, tB.y);
                    if (na <= T(1e-12) || nb <= T(1e-12)) continue;
                    tA.x /= na; tA.y /= na; tB.x /= nb; tB.y /= nb;
                    gsPoint2D<T> p1(A.x + (L / T(3)) * tA.x, A.y + (L / T(3)) * tA.y);
                    gsPoint2D<T> p2(B.x + (L / T(3)) * tB.x, B.y + (L / T(3)) * tB.y);
                    segCands.push_back(makeCubicBezierInteriorSeg(A, p1, p2, B, m_degree));
                }
                for (T bow : {T(-0.16), T(0.16), T(-0.28), T(0.28)}) {
                    gsPoint2D<T> p1(A.x + (L / T(3)) * dAB.x + bow * L * perp.x,
                                    A.y + (L / T(3)) * dAB.y + bow * L * perp.y);
                    gsPoint2D<T> p2(A.x + (T(2) * L / T(3)) * dAB.x + bow * L * perp.x,
                                    A.y + (T(2) * L / T(3)) * dAB.y + bow * L * perp.y);
                    segCands.push_back(makeCubicBezierInteriorSeg(A, p1, p2, B, m_degree));
                }

                for (const auto& cs : segCands) {
                    bool insideOk = true;
                    for (int s = 1; s < 8; ++s) {
                        gsPoint2D<T> pt = evalCurvePoint(*cs.bsp, static_cast<T>(s) / T(8));
                        if (!isPointStrictlyInsideDomain(pt, mesh.loopPolys, T(1e-5))) { insideOk = false; break; }
                    }
                    if (!insideOk) continue;
                    edges[eIdx].seg = cs;
                    auto res = evalInc(true);
                    bool accept = (!bestOk && res.first) ||
                                  (res.first == bestOk && res.second > bestScore + T(1.0));
                    if (accept) {
                        bestOk = res.first;
                        bestScore = res.second;
                        bestSeg = cs;
                        improved = true;
                    }
                }
                edges[eIdx].seg = improved ? bestSeg : oldSeg;
                return improved;
            };

            auto optimizeIncidentSplitEdge = [&](int eIdx) -> bool {
                const auto& inc = edgeFaces[eIdx];
                if (inc.empty()) return false;
                T curMinAll = T(1e9);
                bool curOkAll = true;
                for (const auto& f : inc) {
                    T mnA, mxA, efA;
                    bool okA = evalFaceSplitQuality(f.isTri, f.idx, mnA, mxA, efA);
                    curOkAll = curOkAll && okA;
                    curMinAll = std::min(curMinAll, okA ? efA : (efA - T(200)));
                }

                T oldFrac = edges[eIdx].splitFrac;
                auto oldH0 = edges[eIdx].h0;
                auto oldH1 = edges[eIdx].h1;
                auto oldMid = edges[eIdx].mid;

                std::vector<T> tCands = {arcLengthSplitParam(*edges[eIdx].seg.bsp, T(0.5))};
                for (const auto& f : inc) {
                    if (f.isTri) {
                        const auto& t = tris[f.idx];
                        for (int k = 0; k < 3; ++k)
                            if (t.e[k] == eIdx)
                                tCands.push_back(edges[eIdx].seg.closestParam(mesh.pts[t.v[(k + 2) % 3]]));
                    } else {
                        const auto& q = quads[f.idx];
                        for (int k = 0; k < 4; ++k)
                            if (q.e[k] == eIdx)
                                tCands.push_back(edges[eIdx].seg.closestParam(edges[q.e[(k + 2) % 4]].mid));
                    }
                }
                for (T tc : {T(0.14), T(0.22), T(0.30), T(0.38), T(0.62), T(0.70), T(0.78), T(0.86)})
                    tCands.push_back(tc);

                T bestScore = curMinAll;
                bool bestOk = curOkAll;
                T bestFrac = oldFrac;
                CurveSeg<T> bestH0 = oldH0, bestH1 = oldH1;
                gsPoint2D<T> bestMid = oldMid;
                bool improved = false;
                for (T tTry : tCands) {
                    T uClamped = std::max(T(0.10), std::min(T(0.90), tTry));
                    resplitEdge(eIdx, uClamped);
                    T minAll = T(1e9);
                    bool okAll = true;
                    for (const auto& f : inc) {
                        T mnA, mxA, efA;
                        bool okA = evalFaceSplitQuality(f.isTri, f.idx, mnA, mxA, efA);
                        okAll = okAll && okA;
                        minAll = std::min(minAll, okA ? efA : (efA - T(200)));
                    }
                    // On curved boundary edges, also verify the actual Coons Jacobian of incident faces doesn't flip
                    if (okAll && edges[eIdx].isBoundary) {
                        for (const auto& f : inc) {
                            if (f.isTri) {
                                const auto& t = tris[f.idx];
                                int v0 = t.v[0], v1 = t.v[1], v2 = t.v[2];
                                CurveSeg<T> e01a, e01b, e12a, e12b, e20a, e20b; gsPoint2D<T> m01, m12, m20, O;
                                getHalvesDir(t.e[0], v0, v1, e01a, e01b, m01);
                                getHalvesDir(t.e[1], v1, v2, e12a, e12b, m12);
                                getHalvesDir(t.e[2], v2, v0, e20a, e20b, m20);
                                if (!internal_splitter::pickOptimalTri1to3CenterCurved(
                                        mesh.pts[v0], mesh.pts[v1], mesh.pts[v2],
                                        e01a, e01b, e12a, e12b, e20a, e20b,
                                        edges[t.e[0]].isBoundary, edges[t.e[1]].isBoundary, edges[t.e[2]].isBoundary, O)) {
                                    okAll = false; break;
                                }
                                auto s01 = makeStraightSeg(O, m01, m_degree), s12 = makeStraightSeg(O, m12, m_degree), s20 = makeStraightSeg(O, m20, m_degree);
                                std::unique_ptr<gsTensorBSpline<2, T>> p0(makeCoonsFromCCWCycle(e01a, s01.reversed(), s20, e20b, m_degree));
                                std::unique_ptr<gsTensorBSpline<2, T>> p1(makeCoonsFromCCWCycle(e12a, s12.reversed(), s01, e01b, m_degree));
                                std::unique_ptr<gsTensorBSpline<2, T>> p2(makeCoonsFromCCWCycle(e20a, s20.reversed(), s12, e12b, m_degree));
                                if (!isPatchDetSingleSign(*p0, 11, T(1e-8)) || !isPatchDetSingleSign(*p1, 11, T(1e-8)) || !isPatchDetSingleSign(*p2, 11, T(1e-8))) {
                                    okAll = false; break;
                                }
                            }
                        }
                    }
                    bool accept = (!bestOk && okAll) ||
                                  (okAll == bestOk && minAll > bestScore + T(1.0));
                    if (accept) {
                        bestOk = okAll;
                        bestScore = minAll;
                        bestFrac = uClamped;
                        bestH0 = edges[eIdx].h0;
                        bestH1 = edges[eIdx].h1;
                        bestMid = edges[eIdx].mid;
                        improved = true;
                    }
                }
                // Phase 2: For INTERIOR split edges, also free `er.mid` in 2D off the straight chord A -> B!
                if (!edges[eIdx].isBoundary) {
                    const gsPoint2D<T>& A = mesh.pts[edges[eIdx].a];
                    const gsPoint2D<T>& B = mesh.pts[edges[eIdx].b];
                    T L = A.dist(B);
                    if (L > T(1e-9)) {
                        gsPoint2D<T> perp(-(B.y - A.y) / L, (B.x - A.x) / L);
                        gsPoint2D<T> baseM = bestMid;
                        for (T off : {T(-0.20), T(-0.10), T(0.10), T(0.20)}) {
                            gsPoint2D<T> candM(baseM.x + off * L * perp.x, baseM.y + off * L * perp.y);
                            if (!isPointStrictlyInsideDomain(candM, mesh.loopPolys, T(1e-5)) ||
                                !isSegmentInsideDomain(A, candM, mesh.loopPolys) ||
                                !isSegmentInsideDomain(candM, B, mesh.loopPolys))
                                continue;
                            edges[eIdx].mid = candM;
                            edges[eIdx].h0 = makeStraightSeg(A, candM, m_degree);
                            edges[eIdx].h1 = makeStraightSeg(candM, B, m_degree);
                            T minAll = T(1e9);
                            bool okAll = true;
                            for (const auto& f : inc) {
                                T mnA, mxA, efA;
                                bool okA = evalFaceSplitQuality(f.isTri, f.idx, mnA, mxA, efA);
                                okAll = okAll && okA;
                                minAll = std::min(minAll, okA ? efA : (efA - T(200)));
                            }
                            bool accept = (!bestOk && okAll) ||
                                          (okAll == bestOk && minAll > bestScore + T(1.0));
                            if (accept) {
                                bestOk = okAll;
                                bestScore = minAll;
                                bestH0 = edges[eIdx].h0;
                                bestH1 = edges[eIdx].h1;
                                bestMid = candM;
                                improved = true;
                            }
                        }
                    }
                }
                if (improved) {
                    edges[eIdx].splitFrac = bestFrac;
                    edges[eIdx].h0 = bestH0;
                    edges[eIdx].h1 = bestH1;
                    edges[eIdx].mid = bestMid;
                    return true;
                }
                edges[eIdx].splitFrac = oldFrac;
                edges[eIdx].h0 = oldH0;
                edges[eIdx].h1 = oldH1;
                edges[eIdx].mid = oldMid;
                return false;
            };

            bool anyMoved = false;
            for (int eIdx : badIntEdges) {
                if (optimizeCurvedInteriorEdge(eIdx)) anyMoved = true;
            }
            for (int eIdx : badEdges) {
                if (optimizeIncidentSplitEdge(eIdx)) anyMoved = true;
            }
            if (anyMoved) {
                if (tryBuildFromCurrentEdges(nullptr, nullptr, outFailedQuads, false))
                    return true;
            }
            // Phase 2: Try symmetric curved Hermite spokes (120 deg -> 90 deg at boundaries)
            if ((!tris.empty() || !badEdges.empty()) && (tris.size() + quads.size() <= 90)) {
                if (tryBuildFromCurrentEdges(nullptr, nullptr, outFailedQuads, true))
                    return true;
            }
            if (outFailedQuads && outFailedQuads->empty()) {
                std::sort(badQuads.begin(), badQuads.end());
                badQuads.erase(std::unique(badQuads.begin(), badQuads.end()), badQuads.end());
                *outFailedQuads = std::move(badQuads);
            }
            return false;
        };

        // 1. Curved-aware max-min-angle triangle pairing + Edmonds Blossom matching + relaxed tier + local swap
        struct PairEdge { int t0, t1; typename TriMeshData<T>::Quad q; T effMin; T score; bool relaxed; };
        std::vector<PairEdge> pairEdges;
        const T relMinPairAng = std::max(T(3), effMinPairAng - T(3.0));
        const T relMaxPairAng = std::min(T(172), effMaxPairAng + T(7.5));
        for (int eIdx = 0; eIdx < static_cast<int>(mesh.edges.size()); ++eIdx) {
            if (mesh.edges[eIdx].isBoundary) continue;
            int tA = -1, kA = -1, tB = -1, kB = -1;
            for (int ti = 0; ti < nT; ++ti) {
                for (int k = 0; k < 3; ++k) {
                    if (mesh.tris[ti].e[k] == eIdx) {
                        if (tA < 0) { tA = ti; kA = k; } else { tB = ti; kB = k; }
                    }
                }
            }
            if (tA < 0 || tB < 0) continue;
            const auto& A = mesh.tris[tA];
            const auto& B = mesh.tris[tB];
            int a = A.v[kA], b = A.v[(kA + 1) % 3], c = A.v[(kA + 2) % 3];
            int d = B.v[(kB + 2) % 3];
            if (B.v[kB] != b || B.v[(kB + 1) % 3] != a) continue;
            int qv[4] = {c, a, d, b};
            int e_bc = A.e[(kA + 1) % 3], e_ca = A.e[(kA + 2) % 3];
            int e_ad = B.e[(kB + 1) % 3], e_db = B.e[(kB + 2) % 3];
            int qe[4] = {e_ca, e_ad, e_db, e_bc};
            T effMin = 0, sc = 0;
            if (evalMeshQuadUnion(mesh, qv, qe, effMinPairAng, effMaxPairAng, effMin, sc, false)) {
                pairEdges.push_back({tA, tB, {{c, a, d, b}, {e_ca, e_ad, e_db, e_bc}}, effMin, sc, false});
            } else if (evalMeshQuadUnion(mesh, qv, qe, relMinPairAng, relMaxPairAng, effMin, sc, true)) {
                pairEdges.push_back({tA, tB, {{c, a, d, b}, {e_ca, e_ad, e_db, e_bc}}, effMin, sc - T(30), true});
            }
        }
        std::sort(pairEdges.begin(), pairEdges.end(),
                  [](const PairEdge& p, const PairEdge& q) { return p.score > q.score; });

        std::vector<int> initBaseDeg(mesh.pts.size(), 0);
        for (const auto& t : mesh.tris) { ++initBaseDeg[t.v[0]]; ++initBaseDeg[t.v[1]]; ++initBaseDeg[t.v[2]]; }
        std::vector<std::vector<std::pair<int, int>>> triAdjStrict(nT), triAdjAll(nT);
        for (int pi = 0; pi < static_cast<int>(pairEdges.size()); ++pi) {
            const auto& pe = pairEdges[pi];
            triAdjAll[pe.t0].push_back({pe.t1, pi});
            triAdjAll[pe.t1].push_back({pe.t0, pi});
            if (!pe.relaxed) {
                triAdjStrict[pe.t0].push_back({pe.t1, pi});
                triAdjStrict[pe.t1].push_back({pe.t0, pi});
            }
        }
        std::vector<T> triMinAng(nT, T(60));
        for (int ti = 0; ti < nT; ++ti) {
            T mxA = T(0);
            triMinMaxAnglesDeg(mesh.pts[mesh.tris[ti].v[0]],
                               mesh.pts[mesh.tris[ti].v[1]],
                               mesh.pts[mesh.tris[ti].v[2]], triMinAng[ti], mxA);
        }

        auto runMatchingStrategy = [&](T minEffThresh, bool allowRelaxed,
                                       std::vector<typename TriMeshData<T>::Quad>& joinedQuadsOut,
                                       std::vector<typename TriMeshData<T>::Tri>& remTrisOut,
                                       std::vector<int>* remTriIdxOut = nullptr) {
            std::vector<int> baseDeg = initBaseDeg;
            std::vector<int> matchedPair(nT, -1), partnerTri(nT, -1);

            for (int tier = 0; tier < 2; ++tier) {
                T thresh = (tier == 0) ? minEffThresh : T(0);
                if (tier == 1 && minEffThresh <= T(0)) break;
                for (int pi = 0; pi < static_cast<int>(pairEdges.size()); ++pi) {
                    const auto& pe = pairEdges[pi];
                    if (pe.relaxed && !allowRelaxed) continue;
                    if (pe.effMin < thresh) continue;
                    if (matchedPair[pe.t0] >= 0 || matchedPair[pe.t1] >= 0) continue;
                    int a = pe.q.v[1], b = pe.q.v[3];
                    if (baseDeg[a] <= internal_splitter::minRequiredBoundaryTriValence(mesh, a) ||
                        baseDeg[b] <= internal_splitter::minRequiredBoundaryTriValence(mesh, b))
                        continue;
                    matchedPair[pe.t0] = pi; matchedPair[pe.t1] = pi;
                    partnerTri[pe.t0] = pe.t1; partnerTri[pe.t1] = pe.t0;
                    --baseDeg[a]; --baseDeg[b];
                }
            }

            // Edmonds' Blossom augmenting paths (handles odd cycles around holes and odd-valence vertices)
            auto runBlossomPass = [&](const std::vector<std::vector<std::pair<int, int>>>& useAdj) {
                std::vector<int> orderUnmatched;
                for (int ti = 0; ti < nT; ++ti)
                    if (matchedPair[ti] < 0) orderUnmatched.push_back(ti);
                std::sort(orderUnmatched.begin(), orderUnmatched.end(),
                          [&](int a, int b) { return triMinAng[a] < triMinAng[b]; });

                for (int startTi : orderUnmatched) {
                    if (matchedPair[startTi] >= 0) continue;
                    internal_splitter::augmentBlossomMatchingFrom(
                        startTi, nT, useAdj, partnerTri, matchedPair,
                        [&](const std::vector<int>& addPairs, const std::vector<int>& remPairs) -> bool {
                            std::map<int, int> delta;
                            for (int pi : remPairs) {
                                ++delta[pairEdges[pi].q.v[1]];
                                ++delta[pairEdges[pi].q.v[3]];
                            }
                            for (int pi : addPairs) {
                                --delta[pairEdges[pi].q.v[1]];
                                --delta[pairEdges[pi].q.v[3]];
                            }
                            for (const auto& kv : delta) {
                                if (baseDeg[kv.first] + kv.second < internal_splitter::minRequiredBoundaryTriValence(mesh, kv.first))
                                    return false;
                            }
                            for (const auto& kv : delta)
                                baseDeg[kv.first] += kv.second;
                            return true;
                        });
                }
            };
            runBlossomPass(triAdjStrict);
            if (allowRelaxed)
                runBlossomPass(triAdjAll);

            // Local 1-edge shift improvement: replace (t0, t1) with higher-angle (t1, freeT2)
            for (int sweep = 0; sweep < 2; ++sweep) {
                bool changed = false;
                for (int pi = 0; pi < static_cast<int>(pairEdges.size()); ++pi) {
                    const auto& newPe = pairEdges[pi];
                    if (newPe.relaxed && !allowRelaxed) continue;
                    int tA = newPe.t0, tB = newPe.t1;
                    bool aFree = (matchedPair[tA] < 0), bFree = (matchedPair[tB] < 0);
                    if (aFree == bFree) continue;
                    int freeT = aFree ? tA : tB;
                    int occT  = aFree ? tB : tA;
                    int oldPi = matchedPair[occT];
                    const auto& oldPe = pairEdges[oldPi];
                    if (newPe.effMin <= oldPe.effMin + T(3.0)) continue;
                    int occMate = partnerTri[occT];
                    int ca = newPe.q.v[1], cb = newPe.q.v[3];
                    bool sameVerts = (ca == oldPe.q.v[1] && cb == oldPe.q.v[3]) || (ca == oldPe.q.v[3] && cb == oldPe.q.v[1]);
                    if (!sameVerts && (baseDeg[ca] <= internal_splitter::minRequiredBoundaryTriValence(mesh, ca) ||
                                       baseDeg[cb] <= internal_splitter::minRequiredBoundaryTriValence(mesh, cb)))
                        continue;
                    ++baseDeg[oldPe.q.v[1]]; ++baseDeg[oldPe.q.v[3]];
                    matchedPair[occMate] = -1; partnerTri[occMate] = -1;
                    matchedPair[freeT] = pi; matchedPair[occT] = pi;
                    partnerTri[freeT] = occT; partnerTri[occT] = freeT;
                    --baseDeg[ca]; --baseDeg[cb];
                    changed = true;
                }
                if (!changed) break;
            }

            joinedQuadsOut.clear();
            remTrisOut.clear();
            if (remTriIdxOut) remTriIdxOut->clear();
            for (int ti = 0; ti < nT; ++ti) {
                if (matchedPair[ti] >= 0) {
                    if (partnerTri[ti] > ti) joinedQuadsOut.push_back(pairEdges[matchedPair[ti]].q);
                } else {
                    remTrisOut.push_back(mesh.tris[ti]);
                    if (remTriIdxOut) remTriIdxOut->push_back(ti);
                }
            }
        };

        auto propagateAndEmit = [&](const std::vector<typename TriMeshData<T>::Quad>& joinedQuads,
                                    const std::vector<typename TriMeshData<T>::Tri>& remTris) -> bool {
            if (joinedQuads.empty()) return false;
            auto propagateMask = [&](std::vector<bool>& mask) {
                bool changed = true;
                while (changed) {
                    changed = false;
                    for (const auto& q : joinedQuads) {
                        if (mask[q.e[0]] != mask[q.e[2]]) {
                            mask[q.e[0]] = mask[q.e[2]] = true;
                            changed = true;
                        }
                        if (mask[q.e[1]] != mask[q.e[3]]) {
                            mask[q.e[1]] = mask[q.e[3]] = true;
                            changed = true;
                        }
                    }
                }
            };

            // Dual-Path Unmatched-Triangle Zipper:
            // When remTris.size() is even (>= 2), pair unmatched triangles via vertex-disjoint dual
            // paths through interior quad edges ONLY. Each path turns 2 tris + m quads into
            // 2 + 2m (opposite sides) or 2 + 3m (adjacent sides) strictly convex conforming quads
            // WITHOUT splitting any boundary edge or propagating needSplit across the rest of the mesh!
            auto tryZipUnmatchedTriPairs = [&]() -> bool {
                const int nR = static_cast<int>(remTris.size());
                const int nQ = static_cast<int>(joinedQuads.size());
                if (nR < 2 || (nR % 2) != 0 || nQ == 0) return false;

                // Map each interior mesh edge to the two incident mixed faces:
                // faceId in [0, nQ-1] = quad, [nQ, nQ+nR-1] = remTri
                const int nF = nQ + nR;
                struct FaceSide { int fId; int side; };
                std::vector<std::vector<FaceSide>> edgeInc(mesh.edges.size());
                for (int qi = 0; qi < nQ; ++qi)
                    for (int s = 0; s < 4; ++s)
                        edgeInc[joinedQuads[qi].e[s]].push_back({qi, s});
                for (int ri = 0; ri < nR; ++ri)
                    for (int s = 0; s < 3; ++s)
                        edgeInc[remTris[ri].e[s]].push_back({nQ + ri, s});

                struct DualAdj { int toF; int eId; int mySide; int toSide; };
                std::vector<std::vector<DualAdj>> dualAdj(nF);
                for (int eid = 0; eid < static_cast<int>(mesh.edges.size()); ++eid) {
                    if (mesh.edges[eid].isBoundary) continue;
                    if (edgeInc[eid].size() != 2) continue;
                    auto a = edgeInc[eid][0], b = edgeInc[eid][1];
                    dualAdj[a.fId].push_back({b.fId, eid, a.side, b.side});
                    dualAdj[b.fId].push_back({a.fId, eid, b.side, a.side});
                }

                std::vector<bool> faceUsed(nF, false);
                struct ZipQuadSpec { CurveSeg<T> e[4]; gsPoint2D<T> c[4]; };
                std::vector<ZipQuadSpec> extraQuads;

                auto getMeshSegDir = [&](int eIdx, int u, int v) -> CurveSeg<T> {
                    const auto& er = mesh.edges[eIdx];
                    return (er.a == u && er.b == v) ? er.seg : er.seg.reversed();
                };

                int zippedTris = 0;
                for (int startRi = 0; startRi < nR; ++startRi) {
                    int sF = nQ + startRi;
                    if (faceUsed[sF]) continue;

                    // BFS shortest path through unused quads to another unused remTri
                    std::vector<int> prevF(nF, -1), prevEdge(nF, -1), prevSideFrom(nF, -1), prevSideTo(nF, -1);
                    std::queue<int> bfs;
                    prevF[sF] = sF;
                    bfs.push(sF);
                    int targetF = -1;
                    while (!bfs.empty() && targetF < 0) {
                        int u = bfs.front(); bfs.pop();
                        for (const auto& nb : dualAdj[u]) {
                            int v = nb.toF;
                            if (faceUsed[v] || prevF[v] >= 0) continue;
                            // Skip direct tri-tri edge if union isn't convex; routing through >=1 quad always works
                            if (u == sF && v >= nQ) {
                                const auto& A = remTris[u - nQ];
                                const auto& B = remTris[v - nQ];
                                int kA = nb.mySide, kB = nb.toSide;
                                gsPoint2D<T> qc[4] = {
                                    mesh.pts[A.v[(kA + 2) % 3]], mesh.pts[A.v[kA]],
                                    mesh.pts[B.v[(kB + 2) % 3]], mesh.pts[A.v[(kA + 1) % 3]]
                                };
                                if (!isStrictlyConvexQuad(qc)) continue;
                            }
                            prevF[v] = u;
                            prevEdge[v] = nb.eId;
                            prevSideFrom[v] = nb.mySide;
                            prevSideTo[v] = nb.toSide;
                            if (v >= nQ) { targetF = v; break; }
                            bfs.push(v);
                        }
                    }
                    if (targetF < 0) return false;

                    // Reconstruct path faces[0..L] (faces[0]=sF, faces[L]=targetF) and edges[0..L-1]
                    std::vector<int> pathF, pathE, sideOut, sideIn;
                    for (int cur = targetF; cur != sF; cur = prevF[cur]) {
                        pathF.push_back(cur);
                        pathE.push_back(prevEdge[cur]);
                        sideOut.push_back(prevSideFrom[cur]);
                        sideIn.push_back(prevSideTo[cur]);
                    }
                    pathF.push_back(sF);
                    std::reverse(pathF.begin(), pathF.end());
                    std::reverse(pathE.begin(), pathE.end());
                    std::reverse(sideOut.begin(), sideOut.end());
                    std::reverse(sideIn.begin(), sideIn.end());

                    const int L = static_cast<int>(pathE.size());
                    // Direct convex tri-tri merge (L == 1)
                    if (L == 1) {
                        const auto& A = remTris[pathF[0] - nQ];
                        const auto& B = remTris[pathF[1] - nQ];
                        int kA = sideOut[0], kB = sideIn[0];
                        int a = A.v[kA], b = A.v[(kA + 1) % 3], c = A.v[(kA + 2) % 3], d = B.v[(kB + 2) % 3];
                        ZipQuadSpec zq;
                        zq.e[0] = getMeshSegDir(A.e[(kA + 2) % 3], c, a);
                        zq.e[1] = getMeshSegDir(B.e[(kB + 1) % 3], a, d);
                        zq.e[2] = getMeshSegDir(B.e[(kB + 2) % 3], d, b);
                        zq.e[3] = getMeshSegDir(A.e[(kA + 1) % 3], b, c);
                        zq.c[0] = mesh.pts[c]; zq.c[1] = mesh.pts[a]; zq.c[2] = mesh.pts[d]; zq.c[3] = mesh.pts[b];
                        if (!isStrictlyConvexQuad(zq.c)) return false;
                        extraQuads.push_back(zq);
                        faceUsed[pathF[0]] = faceUsed[pathF[1]] = true;
                        zippedTris += 2;
                        continue;
                    }

                    // Path length L >= 2: place 1 interior split point M[i] on each pathE[i] (i = 0..L-1)
                    std::vector<gsPoint2D<T>> M(L);
                    auto quadCentroid = [&](int qIdx) -> gsPoint2D<T> {
                        const auto& q = joinedQuads[qIdx];
                        return gsPoint2D<T>(
                            T(0.25) * (mesh.pts[q.v[0]].x + mesh.pts[q.v[1]].x + mesh.pts[q.v[2]].x + mesh.pts[q.v[3]].x),
                            T(0.25) * (mesh.pts[q.v[0]].y + mesh.pts[q.v[1]].y + mesh.pts[q.v[2]].y + mesh.pts[q.v[3]].y));
                    };

                    // Pick M[0] shifted slightly into first quad pathF[1] so tri pathF[0] becomes strictly convex quad
                    auto pickEndShift = [&](int triF, int triSide, int qF, gsPoint2D<T>& outM) -> bool {
                        const auto& t = remTris[triF - nQ];
                        int u = t.v[triSide], v = t.v[(triSide + 1) % 3], opp = t.v[(triSide + 2) % 3];
                        gsPoint2D<T> mid(T(0.5) * (mesh.pts[u].x + mesh.pts[v].x), T(0.5) * (mesh.pts[u].y + mesh.pts[v].y));
                        gsPoint2D<T> ctr = quadCentroid(qF);
                        const auto& q = joinedQuads[qF];
                        gsPoint2D<T> qp[4] = {mesh.pts[q.v[0]], mesh.pts[q.v[1]], mesh.pts[q.v[2]], mesh.pts[q.v[3]]};
                        T bestSc = -T(1);
                        bool ok = false;
                        for (T al : {T(0.22), T(0.30), T(0.15), T(0.36)}) {
                            gsPoint2D<T> cand(mid.x + al * (ctr.x - mid.x), mid.y + al * (ctr.y - mid.y));
                            if (!isPointStrictlyInsideDomain(cand, mesh.loopPolys, T(1e-5))) continue;
                            // Must stay strictly inside q
                            bool inQ = true;
                            for (int k = 0; k < 4; ++k)
                                if (cross2(qp[k], qp[(k + 1) % 4], cand) <= T(1e-7)) { inQ = false; break; }
                            if (!inQ) continue;
                            // CCW quad for tri: (opp, u, cand, v)
                            gsPoint2D<T> tq[4] = {mesh.pts[opp], mesh.pts[u], cand, mesh.pts[v]};
                            if (!isStrictlyConvexQuad(tq)) continue;
                            if (!isSegmentInsideDomain(mesh.pts[u], cand, mesh.loopPolys) ||
                                !isSegmentInsideDomain(cand, mesh.pts[v], mesh.loopPolys))
                                continue;
                            CurveSeg<T> e0 = getMeshSegDir(t.e[(triSide + 2) % 3], opp, u);
                            CurveSeg<T> e1 = makeStraightSeg(mesh.pts[u], cand, m_degree);
                            CurveSeg<T> e2 = makeStraightSeg(cand, mesh.pts[v], m_degree);
                            CurveSeg<T> e3 = getMeshSegDir(t.e[(triSide + 1) % 3], v, opp);
                            T mnA, mxA, efA;
                            eval4SegCornerAnglesDeg(e0, e1, e2, e3, mnA, mxA, efA);
                            if (mnA < effMinValAng || mxA > effMaxValAng) continue;
                            if (!ok || efA > bestSc) { bestSc = efA; outM = cand; ok = true; }
                        }
                        return ok;
                    };

                    if (!pickEndShift(pathF[0], sideOut[0], pathF[1], M[0])) return false;
                    if (!pickEndShift(pathF[L], sideIn[L - 1], pathF[L - 1], M[L - 1])) return false;
                    for (int i = 1; i < L - 1; ++i) {
                        const auto& er = mesh.edges[pathE[i]];
                        M[i] = gsPoint2D<T>(T(0.5) * (mesh.pts[er.a].x + mesh.pts[er.b].x),
                                            T(0.5) * (mesh.pts[er.a].y + mesh.pts[er.b].y));
                    }

                    std::vector<ZipQuadSpec> pathQuads;
                    // End triangle 0 -> 1 quad: (opp, u, M[0], v)
                    {
                        const auto& t = remTris[pathF[0] - nQ];
                        int s = sideOut[0];
                        int u = t.v[s], v = t.v[(s + 1) % 3], opp = t.v[(s + 2) % 3];
                        ZipQuadSpec zq;
                        zq.e[0] = getMeshSegDir(t.e[(s + 2) % 3], opp, u);
                        zq.e[1] = makeStraightSeg(mesh.pts[u], M[0], m_degree);
                        zq.e[2] = makeStraightSeg(M[0], mesh.pts[v], m_degree);
                        zq.e[3] = getMeshSegDir(t.e[(s + 1) % 3], v, opp);
                        zq.c[0] = mesh.pts[opp]; zq.c[1] = mesh.pts[u]; zq.c[2] = M[0]; zq.c[3] = mesh.pts[v];
                        pathQuads.push_back(zq);
                    }
                    // End triangle L -> 1 quad: (opp, u, M[L-1], v)
                    {
                        const auto& t = remTris[pathF[L] - nQ];
                        int s = sideIn[L - 1];
                        int u = t.v[s], v = t.v[(s + 1) % 3], opp = t.v[(s + 2) % 3];
                        ZipQuadSpec zq;
                        zq.e[0] = getMeshSegDir(t.e[(s + 2) % 3], opp, u);
                        zq.e[1] = makeStraightSeg(mesh.pts[u], M[L - 1], m_degree);
                        zq.e[2] = makeStraightSeg(M[L - 1], mesh.pts[v], m_degree);
                        zq.e[3] = getMeshSegDir(t.e[(s + 1) % 3], v, opp);
                        zq.c[0] = mesh.pts[opp]; zq.c[1] = mesh.pts[u]; zq.c[2] = M[L - 1]; zq.c[3] = mesh.pts[v];
                        pathQuads.push_back(zq);
                    }

                    // Intermediate quads pathF[1..L-1]: each has 2 split sides (sA = sideIn[i-1] with M[i-1], sB = sideOut[i] with M[i])
                    bool pathOk = true;
                    for (int i = 1; i < L && pathOk; ++i) {
                        const auto& q = joinedQuads[pathF[i]];
                        int sA = sideIn[i - 1], sB = sideOut[i];
                        const gsPoint2D<T>& mA = M[i - 1];
                        const gsPoint2D<T>& mB = M[i];
                        if (sA == sB) { pathOk = false; break; }
                        int diff = (sB - sA + 4) % 4;
                        if (diff == 2) {
                            // Opposite sides: rotate so sA is side 0 (v0->v1), sB is side 2 (v2->v3)
                            int r = sA;
                            int v0 = q.v[r], v1 = q.v[(r + 1) % 4], v2 = q.v[(r + 2) % 4], v3 = q.v[(r + 3) % 4];
                            int e1 = q.e[(r + 1) % 4], e3 = q.e[(r + 3) % 4];
                            if (!isSegmentInsideDomain(mA, mB, mesh.loopPolys)) { pathOk = false; break; }
                            CurveSeg<T> sMid = makeStraightSeg(mA, mB, m_degree);
                            ZipQuadSpec qa, qb;
                            qa.e[0] = makeStraightSeg(mesh.pts[v0], mA, m_degree);
                            qa.e[1] = sMid;
                            qa.e[2] = makeStraightSeg(mB, mesh.pts[v3], m_degree);
                            qa.e[3] = getMeshSegDir(e3, v3, v0);
                            qa.c[0] = mesh.pts[v0]; qa.c[1] = mA; qa.c[2] = mB; qa.c[3] = mesh.pts[v3];

                            qb.e[0] = makeStraightSeg(mA, mesh.pts[v1], m_degree);
                            qb.e[1] = getMeshSegDir(e1, v1, v2);
                            qb.e[2] = makeStraightSeg(mesh.pts[v2], mB, m_degree);
                            qb.e[3] = sMid.reversed();
                            qb.c[0] = mA; qb.c[1] = mesh.pts[v1]; qb.c[2] = mesh.pts[v2]; qb.c[3] = mB;

                            if (!isStrictlyConvexQuad(qa.c) || !isStrictlyConvexQuad(qb.c)) { pathOk = false; break; }
                            pathQuads.push_back(qa);
                            pathQuads.push_back(qb);
                        } else {
                            // Adjacent sides: rotate so split sides are side 0 (v0->v1, m0) and side 1 (v1->v2, m1),
                            // untouched sides are side 2 (v2->v3) and side 3 (v3->v0), opposite corner is v3!
                            int r = (diff == 1) ? sA : sB;
                            const gsPoint2D<T>& m0 = (diff == 1) ? mA : mB;
                            const gsPoint2D<T>& m1 = (diff == 1) ? mB : mA;
                            int v0 = q.v[r], v1 = q.v[(r + 1) % 4], v2 = q.v[(r + 2) % 4], v3 = q.v[(r + 3) % 4];
                            int e2 = q.e[(r + 2) % 4], e3 = q.e[(r + 3) % 4];
                            gsPoint2D<T> ctr = quadCentroid(pathF[i]);
                            gsPoint2D<T> oCands[4] = {
                                gsPoint2D<T>((m0.x + m1.x + mesh.pts[v3].x) / T(3), (m0.y + m1.y + mesh.pts[v3].y) / T(3)),
                                ctr,
                                gsPoint2D<T>(T(0.5) * (ctr.x + mesh.pts[v3].x), T(0.5) * (ctr.y + mesh.pts[v3].y)),
                                gsPoint2D<T>(T(0.25) * (m0.x + m1.x + mesh.pts[v0].x + mesh.pts[v2].x),
                                             T(0.25) * (m0.y + m1.y + mesh.pts[v0].y + mesh.pts[v2].y))
                            };
                            gsPoint2D<T> bestO = ctr;
                            T bestEf = -T(1);
                            bool foundO = false;
                            for (int ci = 0; ci < 4; ++ci) {
                                const auto& O = oCands[ci];
                                if (!isPointStrictlyInsideDomain(O, mesh.loopPolys, T(1e-5))) continue;
                                gsPoint2D<T> q0[4] = {mesh.pts[v0], m0, O, mesh.pts[v3]};
                                gsPoint2D<T> q1[4] = {m0, mesh.pts[v1], m1, O};
                                gsPoint2D<T> q2[4] = {O, m1, mesh.pts[v2], mesh.pts[v3]};
                                if (!isStrictlyConvexQuad(q0) || !isStrictlyConvexQuad(q1) || !isStrictlyConvexQuad(q2)) continue;
                                if (!isSegmentInsideDomain(m0, O, mesh.loopPolys) ||
                                    !isSegmentInsideDomain(m1, O, mesh.loopPolys) ||
                                    !isSegmentInsideDomain(mesh.pts[v3], O, mesh.loopPolys))
                                    continue;
                                CurveSeg<T> sm0 = makeStraightSeg(m0, O, m_degree);
                                CurveSeg<T> sm1 = makeStraightSeg(O, m1, m_degree);
                                CurveSeg<T> sv3 = makeStraightSeg(O, mesh.pts[v3], m_degree);
                                T mn0, mx0, ef0, mn1, mx1, ef1, mn2, mx2, ef2;
                                eval4SegCornerAnglesDeg(makeStraightSeg(mesh.pts[v0], m0, m_degree), sm0, sv3, getMeshSegDir(e3, v3, v0), mn0, mx0, ef0);
                                eval4SegCornerAnglesDeg(makeStraightSeg(m0, mesh.pts[v1], m_degree), makeStraightSeg(mesh.pts[v1], m1, m_degree), sm1.reversed(), sm0.reversed(), mn1, mx1, ef1);
                                eval4SegCornerAnglesDeg(sm1, makeStraightSeg(m1, mesh.pts[v2], m_degree), getMeshSegDir(e2, v2, v3), sv3.reversed(), mn2, mx2, ef2);
                                T efMin = std::min(ef0, std::min(ef1, ef2));
                                if (!foundO || efMin > bestEf) { bestEf = efMin; bestO = O; foundO = true; }
                            }
                            if (!foundO) { pathOk = false; break; }
                            CurveSeg<T> sm0 = makeStraightSeg(m0, bestO, m_degree);
                            CurveSeg<T> sm1 = makeStraightSeg(bestO, m1, m_degree);
                            CurveSeg<T> sv3 = makeStraightSeg(bestO, mesh.pts[v3], m_degree);
                            ZipQuadSpec q0, q1, q2;
                            q0.e[0] = makeStraightSeg(mesh.pts[v0], m0, m_degree); q0.e[1] = sm0; q0.e[2] = sv3; q0.e[3] = getMeshSegDir(e3, v3, v0);
                            q0.c[0] = mesh.pts[v0]; q0.c[1] = m0; q0.c[2] = bestO; q0.c[3] = mesh.pts[v3];
                            q1.e[0] = makeStraightSeg(m0, mesh.pts[v1], m_degree); q1.e[1] = makeStraightSeg(mesh.pts[v1], m1, m_degree); q1.e[2] = sm1.reversed(); q1.e[3] = sm0.reversed();
                            q1.c[0] = m0; q1.c[1] = mesh.pts[v1]; q1.c[2] = m1; q1.c[3] = bestO;
                            q2.e[0] = sm1; q2.e[1] = makeStraightSeg(m1, mesh.pts[v2], m_degree); q2.e[2] = getMeshSegDir(e2, v2, v3); q2.e[3] = sv3.reversed();
                            q2.c[0] = bestO; q2.c[1] = m1; q2.c[2] = mesh.pts[v2]; q2.c[3] = mesh.pts[v3];
                            pathQuads.push_back(q0);
                            pathQuads.push_back(q1);
                            pathQuads.push_back(q2);
                        }
                    }
                    if (!pathOk) return false;
                    for (int f : pathF) faceUsed[f] = true;
                    for (auto& pq : pathQuads) extraQuads.push_back(std::move(pq));
                    zippedTris += 2;
                }
                if (zippedTris != nR) return false;

                // Assemble pure-quad multipatch (untouched quads + zipped path quads)
                gsMultiPatch<T> mp;
                for (int qi = 0; qi < nQ; ++qi) {
                    if (faceUsed[qi]) continue;
                    const auto& q = joinedQuads[qi];
                    gsPoint2D<T> c[4] = {mesh.pts[q.v[0]], mesh.pts[q.v[1]], mesh.pts[q.v[2]], mesh.pts[q.v[3]]};
                    bool anyBnd = mesh.edges[q.e[0]].isBoundary || mesh.edges[q.e[1]].isBoundary ||
                                  mesh.edges[q.e[2]].isBoundary || mesh.edges[q.e[3]].isBoundary;
                    if (!anyBnd && !isStrictlyConvexQuad(c)) return false;
                    mp.addPatch(typename gsGeometry<T>::uPtr(
                        makeCoonsFromCCWCycle(getMeshSegDir(q.e[0], q.v[0], q.v[1]),
                                              getMeshSegDir(q.e[1], q.v[1], q.v[2]),
                                              getMeshSegDir(q.e[2], q.v[2], q.v[3]),
                                              getMeshSegDir(q.e[3], q.v[3], q.v[0]), m_degree)));
                }
                for (const auto& zq : extraQuads) {
                    mp.addPatch(typename gsGeometry<T>::uPtr(
                        makeCoonsFromCCWCycle(zq.e[0], zq.e[1], zq.e[2], zq.e[3], m_degree)));
                }
                if (mp.nPatches() == 0) return false;
                for (size_t p = 0; p < mp.nPatches(); ++p) {
                    T angs[4], mnA = 0, mxA = 180; bool degen = false;
                    evalPatchCornerAnglesDeg(mp.patch(p), angs, mnA, mxA, degen);
                    if (degen || mnA < effMinValAng || mxA > effMaxValAng) return false;
                }
                gsMultiPatch<T> preWinslow = mp;
                optimizeWinslowBarrier(mp);
                mp.computeTopology(T(1e-5), false, false);
                std::string valErr;
                bool valOk = validateMultiPatch(m_domain, mp, T(1e-4), &valErr, effMinValAng, effMaxValAng);
                if (!valOk) {
                    preWinslow.computeTopology(T(1e-5), false, false);
                    if (validateMultiPatch(m_domain, preWinslow, T(1e-4), &valErr, effMinValAng, effMaxValAng)) {
                        mp = give(preWinslow);
                        valOk = true;
                    }
                }
                if (!valOk) return false;
                gsASG1LocalFitter<T>::unifyPatchOrientations(mp);
                mp.computeTopology(T(1e-5), false, false);
                for (size_t p = 0; p < mp.nPatches(); ++p) {
                    T d49 = internal_splitter::evalPatchMinDetDense(mp.patch(p), 49);
                    if (d49 <= T(1e-8)) {
                        if (auto* tb = dynamic_cast<gsTensorBSpline<2, T>*>(&mp.patch(p))) {
                            internal_splitter::quickHarmonicPatchRepair(*tb, T(1e-7));
                            d49 = internal_splitter::evalPatchMinDetDense(*tb, 49);
                        }
                    }
                    if (d49 <= T(1e-8)) return false;
                }

                T minDet = 0, minA = 0, maxA = 180;
                T energy = evaluateDistortionEnergy(mp, &minDet, &minA, &maxA);
                SplittingCandidate cand;
                cand.multiPatch = give(mp);
                cand.numPatches = static_cast<int>(cand.multiPatch.nPatches());
                cand.minJacobianDet = minDet;
                cand.distortionEnergy = energy;
                cand.minCornerAngleDeg = minA;
                cand.maxCornerAngleDeg = maxA;
                cand.methodUsed = "MultiHoleTriangulation";
                cand.isValid = true;
                updateBest(cand);
                return true;
            };

            int prevP = bestCand.isValid ? bestCand.numPatches : 999999;
            if (!remTris.empty() && tryZipUnmatchedTriPairs()) {
                return bestCand.isValid && bestCand.numPatches < prevP;
            }

            std::vector<bool> needSplit(mesh.edges.size(), false);
            for (const auto& t : remTris)
                needSplit[t.e[0]] = needSplit[t.e[1]] = needSplit[t.e[2]] = true;
            propagateMask(needSplit);

            std::vector<int> failedQuads;
            if (emitFromMixed(remTris, joinedQuads, needSplit, &failedQuads)) {
                return bestCand.isValid && bestCand.numPatches < prevP;
            }

            // Selective single-strip split rescue: if only a few quads failed, split 1 polychord strip
            // through a failing quad instead of quadrupling every quad via allTrue!
            if (!failedQuads.empty()) {
                int triedStrips = 0;
                for (int qi : failedQuads) {
                    if (triedStrips >= 4) break;
                    if (qi < 0 || qi >= static_cast<int>(joinedQuads.size())) continue;
                    const auto& q = joinedQuads[qi];
                    for (int dir = 0; dir < 2; ++dir) {
                        int ea = q.e[dir], eb = q.e[dir + 2];
                        if (needSplit[ea] && needSplit[eb]) continue;
                        std::vector<bool> stripMask = needSplit;
                        stripMask[ea] = stripMask[eb] = true;
                        propagateMask(stripMask);
                        bool isAllTrue = true;
                        for (size_t ei = 0; ei < stripMask.size(); ++ei)
                            if (!stripMask[ei]) { isAllTrue = false; break; }
                        if (isAllTrue) continue;
                        ++triedStrips;
                        if (emitFromMixed(remTris, joinedQuads, stripMask))
                            return bestCand.isValid && bestCand.numPatches < prevP;
                    }
                }
            }

            // Last resort on this matching: global 1-to-4 / 1-to-3 split ONLY when partial/strip failed
            std::vector<bool> allTrue(mesh.edges.size(), true);
            emitFromMixed(remTris, joinedQuads, allTrue);
            return bestCand.isValid && bestCand.numPatches < prevP;
        };

        std::vector<typename TriMeshData<T>::Quad> jQuads0, jQuads1, jQuadsR;
        std::vector<typename TriMeshData<T>::Tri> rem0, rem1, remR;
        std::vector<int> remIdx0;
        if (!pairEdges.empty()) {
            runMatchingStrategy(T(0), false, jQuads0, rem0, &remIdx0);
            propagateAndEmit(jQuads0, rem0);

            // If strict matching left unpaired tris (which trigger global needSplit cascades),
            // try relaxed-angle Blossom matching + curved chords to reach 0 unpaired tris!
            if (!rem0.empty()) {
                runMatchingStrategy(T(0), true, jQuadsR, remR);
                if (remR.size() < rem0.size())
                    propagateAndEmit(jQuadsR, remR);
            }

            T highThresh = std::max(effMinPairAng + T(8), T(22));
            runMatchingStrategy(highThresh, false, jQuads1, rem1);
            if (jQuads1.size() != jQuads0.size() || rem1.size() != rem0.size())
                propagateAndEmit(jQuads1, rem1);

            // Parity & unmatched-pair repair:
            // (A) If odd nT leaves unmatched tris, bisect 1 boundary edge of worst unmatched tri.
            // (B) If even nT still leaves 2..10 unmatched tris (and no zero-cascade solution yet),
            //     bisect a shared interior edge between two adjacent unmatched tris OR bisect 2 boundary edges.
            const bool wantRepair = allowParityRepair && nT <= 220 && !rem0.empty() &&
                                    (!bestCand.isValid || bestCand.numPatches > nT);
            if (wantRepair && (nT % 2) == 1) {
                int candTri = -1, candEdge = -1, candK = -1;
                T bestLen = -T(1);
                for (size_t ri = 0; ri < remIdx0.size(); ++ri) {
                    int ti = remIdx0[ri];
                    const auto& rt = mesh.tris[ti];
                    for (int k = 0; k < 3; ++k) {
                        int eid = rt.e[k];
                        if (!mesh.edges[eid].isBoundary) continue;
                        T el = mesh.pts[rt.v[k]].dist(mesh.pts[rt.v[(k + 1) % 3]]);
                        if (el > bestLen) { bestLen = el; candEdge = eid; candK = k; candTri = ti; }
                    }
                }
                if (candTri >= 0 && candEdge >= 0) {
                    const auto bt = mesh.tris[candTri];
                    int u = bt.v[candK], v = bt.v[(candK + 1) % 3], opp = bt.v[(candK + 2) % 3];
                    int e_vo = bt.e[(candK + 1) % 3], e_ou = bt.e[(candK + 2) % 3];
                    auto oldEr = mesh.edges[candEdge];
                    T tMid = arcLengthSplitParam(*oldEr.seg.bsp, T(0.5));
                    T tPrj = std::max(T(0.22), std::min(T(0.78), oldEr.seg.closestParam(mesh.pts[opp])));
                    for (T tTry : {tPrj, tMid}) {
                        TriMeshData<T> rMesh = mesh;
                        auto sp = oldEr.seg.splitAtParam(tTry);
                        gsPoint2D<T> sPt = sp.first.p1();
                        if (cross2(rMesh.pts[u], sPt, rMesh.pts[opp]) <= T(1e-9) ||
                            cross2(sPt, rMesh.pts[v], rMesh.pts[opp]) <= T(1e-9))
                            continue;
                        if (!isSegmentInsideDomain(rMesh.pts[opp], sPt, rMesh.loopPolys)) continue;
                        int mIdx = static_cast<int>(rMesh.pts.size());
                        rMesh.pts.push_back(sPt);
                        rMesh.isBoundaryPt.push_back(true);
                        rMesh.isSmoothBnd.push_back(true);
                        rMesh.isSharpCorner.push_back(false);
                        rMesh.cornerAngleDeg.push_back(T(180));
                        auto din = evalCurveDeriv(*sp.first.bsp, T(1));
                        T dn = std::hypot(din.x, din.y);
                        gsPoint2D<T> tDir = (dn > T(1e-14)) ? gsPoint2D<T>(din.x / dn, din.y / dn) : gsPoint2D<T>(1, 0);
                        rMesh.bndTanIn.push_back(tDir); rMesh.bndTanOut.push_back(tDir);
                        int oldA = oldEr.a, oldB = oldEr.b;
                        rMesh.edges[candEdge] = {oldA, mIdx, sp.first, CurveSeg<T>(), CurveSeg<T>(), gsPoint2D<T>(), false, true, oldEr.loopIdx};
                        int e_mb = static_cast<int>(rMesh.edges.size());
                        rMesh.edges.push_back({mIdx, oldB, sp.second, CurveSeg<T>(), CurveSeg<T>(), gsPoint2D<T>(), false, true, oldEr.loopIdx});
                        int e_om = static_cast<int>(rMesh.edges.size());
                        rMesh.edges.push_back({opp, mIdx, makeStraightSeg(rMesh.pts[opp], sPt, m_degree),
                                               CurveSeg<T>(), CurveSeg<T>(), gsPoint2D<T>(), false, false, -1});
                        auto halfFor = [&](int ep) { return (ep == oldA) ? candEdge : e_mb; };
                        rMesh.tris[candTri] = {{u, mIdx, opp}, {halfFor(u), e_om, e_ou}};
                        rMesh.tris.push_back({{mIdx, v, opp}, {halfFor(v), e_vo, e_om}});
                        applyDomainLawsonFlips(rMesh);
                        flipForBetterQuadUnions(rMesh, effMinPairAng, effMaxPairAng);
                        quadifyTriMesh(rMesh, coarseRelaxed, false);
                        break;
                    }
                }
            } else if (wantRepair && (nT % 2) == 0 && remIdx0.size() >= 2 && remIdx0.size() <= 10) {
                TriMeshData<T> rMesh = mesh;
                bool repaired = false;
                // Subcase 2a: split shared interior edge between two adjacent unmatched triangles
                std::set<int> remSet(remIdx0.begin(), remIdx0.end());
                for (int ti : remIdx0) {
                    if (repaired) break;
                    const auto& rt = rMesh.tris[ti];
                    for (int k = 0; k < 3; ++k) {
                        int eid = rt.e[k];
                        if (rMesh.edges[eid].isBoundary) continue;
                        int mate = -1;
                        for (int tj : remIdx0) {
                            if (tj == ti) continue;
                            if (rMesh.tris[tj].e[0] == eid || rMesh.tris[tj].e[1] == eid || rMesh.tris[tj].e[2] == eid) {
                                mate = tj; break;
                            }
                        }
                        if (mate >= 0) {
                            int ord[3] = {k, (k + 1) % 3, (k + 2) % 3};
                            if (trySplitTriangleEdges(rMesh, ti, ord, m_degree, T(4.0))) {
                                repaired = true;
                                break;
                            }
                        }
                    }
                }
                // Subcase 2b: split 1 edge on each of 2 unmatched tris (keeps nT even)
                if (!repaired) {
                    int nDone = 0;
                    for (int ti : remIdx0) {
                        if (nDone >= 2 || ti >= static_cast<int>(rMesh.tris.size())) break;
                        const auto& bt = rMesh.tris[ti];
                        int ord[3] = {0, 1, 2};
                        T lens[3];
                        for (int k = 0; k < 3; ++k)
                            lens[k] = rMesh.pts[bt.v[k]].dist(rMesh.pts[bt.v[(k + 1) % 3]]);
                        std::sort(ord, ord + 3, [&](int a, int b) {
                            bool ba = rMesh.edges[bt.e[a]].isBoundary, bb = rMesh.edges[bt.e[b]].isBoundary;
                            if (ba != bb) return ba > bb;
                            return lens[a] > lens[b];
                        });
                        if (trySplitTriangleEdges(rMesh, ti, ord, m_degree, T(4.0)))
                            ++nDone;
                    }
                    repaired = (nDone > 0 && (rMesh.tris.size() % 2) == 0);
                }
                if (repaired) {
                    applyDomainLawsonFlips(rMesh);
                    flipForBetterQuadUnions(rMesh, effMinPairAng, effMaxPairAng);
                    quadifyTriMesh(rMesh, coarseRelaxed, false);
                }
            }
        }
        if (!bestCand.isValid)
        {
            std::vector<bool> allTrue(mesh.edges.size(), true);
            emitFromMixed(mesh.tris, {}, allTrue);
        }
    };

    const int earlyCap = std::max(29, 7 * m_domain.numHoles() + 14);
    const T flipMinAng = std::min(minAngTol, std::max(T(4), minAngTol - T(4.0)));
    const T flipMaxAng = std::max(maxAngTol, T(164.2));
    TriMeshData<T> savedCoarseCdt;
    bool hasSavedCoarse = false;

    // Coarse-first passes: true CDT + feature projections + quad-union co-optimization + targeted sliver repair
    const int maxCoarsePasses = coarseOnly ? 2 : 5;
    for (int coarsePass = 0; coarsePass < maxCoarsePasses && (!bestCand.isValid || (coarsePass <= 1 && bestCand.numPatches > earlyCap)); ++coarsePass) {
        const bool cRef = (coarsePass > 0);
        const bool cRefineHoles = (coarsePass >= 3);
        const T cH = (coarsePass <= 1) ? (diag * T(0.65))
                   : ((coarsePass == 2) ? (diag * T(0.45))
                   : ((coarsePass == 3) ? (diag * T(0.32)) : (diag * T(0.25))));
        const T cBend = (coarsePass <= 1) ? T(1.10) : T(0.65);
        const int cNeedleRounds = (coarsePass >= 3) ? 2 : 1;

        for (int featPass = 0; featPass < 2; ++featPass) {
            const bool useFeat = (featPass == 1);
            if (useFeat && bestCand.isValid && (coarsePass > 1 || bestCand.numPatches <= earlyCap)) break;

            for (int mode = 0; mode < 3 && (!bestCand.isValid || (coarsePass <= 1 && bestCand.numPatches > earlyCap)); ++mode) {
                TriMeshData<T> coarse;
                if (!buildConstrainedDomainTriangulation(m_domain, m_degree, cH, cRef, mode, coarse, cBend, cRefineHoles, useFeat, true))
                    continue;
                if (coarsePass == 0 && !useFeat)
                    quadifyTriMesh(coarse, true, true);

                TriMeshData<T> cdt = coarse;
                applyDomainLawsonFlips(cdt);
                fixUndervalentAndReflexCorners(cdt, m_degree, 12);
                applyDomainLawsonFlips(cdt);
                if (!hasSavedCoarse && coarsePass == 0 && !useFeat) {
                    savedCoarseCdt = cdt;
                    hasSavedCoarse = true;
                }
                if (!bestCand.isValid || (coarsePass <= 1 && cdt.tris.size() <= 75 && bestCand.numPatches > earlyCap))
                    quadifyTriMesh(cdt, true, true);

                if (!bestCand.isValid || (coarsePass <= 1 && cdt.tris.size() <= 90 && bestCand.numPatches > earlyCap)) {
                    TriMeshData<T> qOpt = cdt;
                    flipForBetterQuadUnions(qOpt, flipMinAng, flipMaxAng);
                    quadifyTriMesh(qOpt, true, true);
                }

                if (!bestCand.isValid) {
                    TriMeshData<T> mild = cdt;
                    refineNeedleCircumcenters(mild, cH * T(0.55), cNeedleRounds, T(16.0), true, 16);
                    bisectObtuseTriangles(mild, m_degree, T(146), 28);
                    applyDomainLawsonFlips(mild);
                    flipForBetterQuadUnions(mild, flipMinAng, flipMaxAng);
                    quadifyTriMesh(mild, true, true);
                }
            }
        }
    }

    // Defect-local adaptive CDT repair on the coarse mesh BEFORE jumping to global diag/8 or diag/10!
    if (!coarseOnly && !bestCand.isValid && hasSavedCoarse) {
        TriMeshData<T> locMesh = savedCoarseCdt;
        for (int locRound = 0; locRound < 10 && !bestCand.isValid && locMesh.tris.size() <= 180; ++locRound) {
            int worstTi = -1;
            T worstScore = T(-1);
            for (int ti = 0; ti < static_cast<int>(locMesh.tris.size()); ++ti) {
                const auto& t = locMesh.tris[ti];
                T mnA, mxA;
                triMinMaxAnglesDeg(locMesh.pts[t.v[0]], locMesh.pts[t.v[1]], locMesh.pts[t.v[2]], mnA, mxA);
                T sc = std::max(T(0), T(180) - mxA < mnA ? (T(150) - (T(180) - mxA)) : (T(30) - mnA));
                if (mxA > T(135)) sc += (mxA - T(135)) * T(2);
                if (mnA < T(18))  sc += (T(18) - mnA) * T(2.5);
                if (sc > worstScore) { worstScore = sc; worstTi = ti; }
            }
            if (worstTi < 0) break;
            const auto bt = locMesh.tris[worstTi];
            int ord[3] = {0, 1, 2};
            T lens[3];
            for (int k = 0; k < 3; ++k)
                lens[k] = locMesh.pts[bt.v[k]].dist(locMesh.pts[bt.v[(k + 1) % 3]]);
            std::sort(ord, ord + 3, [&](int a, int b) { return lens[a] > lens[b]; });
            if (!trySplitTriangleEdges(locMesh, worstTi, ord, m_degree, T(3.5)))
                break;
            applyDomainLawsonFlips(locMesh);
            flipForBetterQuadUnions(locMesh, flipMinAng, flipMaxAng);
            quadifyTriMesh(locMesh, true, true);
        }
    }

    for (int refinePass = 0; !coarseOnly && refinePass < 2 && !bestCand.isValid; ++refinePass) {
        const bool doRef = (refinePass == 1);
        const T hStep = doRef ? (diag / T(6)) : (diag / T(8));
        for (int mode = 0; mode < 3 && (!bestCand.isValid || !doRef); ++mode) {
            if (mode == 2 && bestCand.isValid) continue;
            TriMeshData<T> coarse;
            if (!buildConstrainedDomainTriangulation(m_domain, m_degree, hStep, doRef, mode, coarse, T(0.44), true, true, true))
                continue;
            TriMeshData<T> cdt = coarse;
            applyDomainLawsonFlips(cdt);
            quadifyTriMesh(cdt, false, true);

            if (!bestCand.isValid) {
                TriMeshData<T> qOpt = cdt;
                flipForBetterQuadUnions(qOpt, flipMinAng, flipMaxAng);
                quadifyTriMesh(qOpt, false, true);
            }

            // Try obtuse/needle bisection WITHOUT global area-based centroid 1-to-3 flooding first
            if (!bestCand.isValid) {
                TriMeshData<T> localRef = cdt;
                bisectObtuseTriangles(localRef, m_degree, T(146), 28);
                applyDomainLawsonFlips(localRef);
                flipForBetterQuadUnions(localRef, flipMinAng, flipMaxAng);
                quadifyTriMesh(localRef, false, true);
            }

            if (!bestCand.isValid) {
                TriMeshData<T> steiner = cdt;
                refineInteriorTriangles(steiner, hStep, 2);
                bisectObtuseTriangles(steiner, m_degree, T(146), 28);
                applyDomainLawsonFlips(steiner);
                flipForBetterQuadUnions(steiner, flipMinAng, flipMaxAng);
                quadifyTriMesh(steiner, false, true);
            }
        }
    }

    if (!coarseOnly && !bestCand.isValid) {
        TriMeshData<T> fine;
        if (buildConstrainedDomainTriangulation(m_domain, m_degree, diag / T(10), true, 0, fine, T(0.44), true, true, true)) {
            applyDomainLawsonFlips(fine);
            TriMeshData<T> fineNoCent = fine;
            bisectObtuseTriangles(fineNoCent, m_degree, T(146), 28);
            applyDomainLawsonFlips(fineNoCent);
            flipForBetterQuadUnions(fineNoCent, flipMinAng, flipMaxAng);
            quadifyTriMesh(fineNoCent, false, true);
            if (!bestCand.isValid) {
                refineInteriorTriangles(fine, diag / T(10), 2);
                bisectObtuseTriangles(fine, m_degree, T(146), 28);
                applyDomainLawsonFlips(fine);
                flipForBetterQuadUnions(fine, flipMinAng, flipMaxAng);
                quadifyTriMesh(fine, false, false);
            }
        }
    }

    // Fallback: Ruppert-style strict-interior circumcenter refinement for acute needles + coarse diag/4 bisect
    if (!coarseOnly && !bestCand.isValid) {
        for (int fPass = 0; fPass < 8 && !bestCand.isValid; ++fPass) {
            T fH = (fPass == 0 || fPass == 4 || fPass == 6 || fPass == 7) ? (diag / T(10))
                 : ((fPass == 1) ? (diag / T(5))
                 : ((fPass == 2 || fPass == 5) ? (diag / T(4)) : (diag / T(7))));
            bool fRef = (fPass != 2 && fPass != 5);
            int fMode = (fPass == 1) ? 1 : 0;
            int nRounds = (fPass == 2) ? 2 : ((fPass >= 3) ? 3 : 2);
            bool runNeedle = (fPass != 2);
            T needleAng = (fPass >= 3) ? T(18) : T(16);
            T maxQ = (fPass == 2) ? T(146) : ((fPass >= 3) ? T(142) : T(144));
            int maxSp = (fPass == 2) ? 28 : ((fPass >= 3) ? 56 : 48);
            T minQ = (fPass >= 3) ? T(16.5) : T(14.5);
            bool allowCentroid = (fPass >= 5);
            bool fRelax = (fPass >= 4);
            T fBend = (fPass == 7) ? T(0.30) : T(0.44);
            TriMeshData<T> needle;
            if (!buildConstrainedDomainTriangulation(m_domain, m_degree, fH, fRef, fMode, needle, fBend, true, false, false) &&
                !buildConstrainedDomainTriangulation(m_domain, m_degree, fH, fRef, fMode, needle, fBend, true, true, false))
                continue;
            applyDomainLawsonFlips(needle);
            // Try needle/obtuse repair without global centroid flooding first
            {
                TriMeshData<T> needleLean = needle;
                if (runNeedle)
                    refineNeedleCircumcenters(needleLean, fH, nRounds, needleAng, allowCentroid);
                bisectObtuseTriangles(needleLean, m_degree, maxQ, maxSp, minQ);
                applyDomainLawsonFlips(needleLean);
                flipForBetterQuadUnions(needleLean, flipMinAng, flipMaxAng);
                quadifyTriMesh(needleLean, fRelax, true);
            }
            if (!bestCand.isValid) {
                refineInteriorTriangles(needle, fH, nRounds);
                if (runNeedle)
                    refineNeedleCircumcenters(needle, fH, nRounds, needleAng, allowCentroid);
                bisectObtuseTriangles(needle, m_degree, maxQ, maxSp, minQ);
                applyDomainLawsonFlips(needle);
                flipForBetterQuadUnions(needle, flipMinAng, flipMaxAng);
                quadifyTriMesh(needle, fRelax, false);
            }
        }
    }

    if (bestCand.isValid && bestCand.numPatches > 10) {
        gsMultiPatch<T> redMp = bestCand.multiPatch;
        reducePolychordStrips(redMp);
        if (redMp.nPatches() > 0 && redMp.nPatches() < bestCand.multiPatch.nPatches()) {
            redMp.computeTopology(T(1e-5), false, false);
            optimizeWinslowBarrier(redMp);
            gsASG1LocalFitter<T>::unifyPatchOrientations(redMp);
            std::string redErr;
            if (validateMultiPatch(m_domain, redMp, T(1e-4), &redErr,
                                   std::min(minAngTol, std::max(T(4), minAngTol - T(3.5))), std::max(maxAngTol, T(164.2)))) {
                bool denseOk = true;
                for (size_t p = 0; p < redMp.nPatches(); ++p) {
                    if (internal_splitter::evalPatchMinDetDense(redMp.patch(p), 49) <= T(1e-8)) {
                        denseOk = false;
                        break;
                    }
                }
                if (denseOk) {
                    T minDet = 0, minA = 0, maxA = 180;
                    T energy = evaluateDistortionEnergy(redMp, &minDet, &minA, &maxA);
                    bestCand.multiPatch = give(redMp);
                    bestCand.numPatches = static_cast<int>(bestCand.multiPatch.nPatches());
                    bestCand.minJacobianDet = minDet;
                    bestCand.distortionEnergy = energy;
                    bestCand.minCornerAngleDeg = minA;
                    bestCand.maxCornerAngleDeg = maxA;
                }
            }
        }
    }

    return bestCand;
}

// ============================================================================
// Approach 2: Medial Axis Transform (MAT) Skeleton Decomposition
// ============================================================================

template <class T>
typename gsPatchSplitter<T>::SplittingCandidate
gsPatchSplitter<T>::computeMedialAxisSplitting()
{
    using namespace internal_splitter;
    SplittingCandidate bestCand;
    bestCand.methodUsed = "MedialAxisTransform";
    bestCand.isValid = false;

    if (m_domain.numLoops() == 0) return bestCand;

    const T domMinAng = getDomainMinCornerAngleDeg();
    const T minAngTol = std::min(static_cast<T>(m_options.getReal("MinCornerAngleDeg")), T(0.85) * domMinAng);
    const T minHoleAngTol = (m_domain.numHoles() > 0)
                          ? std::min(minAngTol, std::max(T(4.0), minAngTol - T(4.0)))
                          : std::min(minAngTol, std::max(T(4.0), minAngTol - T(3.5)));
    const T maxAngTol = static_cast<T>(m_options.getReal("MaxCornerAngleDeg"));
    const T effMaxAng = std::max(maxAngTol, T(164.9));
    const int maxPatchesCap = std::max(1, static_cast<int>(m_options.getInt("MaxPatches")));

    struct MATConfig { T hFrac; T pruneDeg; T juncMerge; };
    const MATConfig configs[6] = {
        {T(1.0) / T(14.0), T(54.0), T(0.68)},
        {T(1.0) / T(20.0), T(48.0), T(0.55)},
        {T(1.0) / T(12.0), T(58.0), T(0.75)},
        {T(1.0) / T(16.0), T(50.0), T(0.60)},
        {T(1.0) / T(22.0), T(45.0), T(0.45)},
        {T(1.0) / T(28.0), T(38.0), T(0.40)}
    };

    for (int ci = 0; ci < 6; ++ci) {
        gsMultiPatch<T> matMp;
        if (!gsPatchSplitterMAT<T>::buildMATMultiPatch(
                m_domain, m_degree, matMp, configs[ci].hFrac, configs[ci].pruneDeg, configs[ci].juncMerge))
            continue;
        if (matMp.nPatches() == 0) continue;

        gsASG1LocalFitter<T>::unifyPatchOrientations(matMp);
        matMp.computeTopology(T(1e-4), false, false);

        auto countNegDense = [&](const gsMultiPatch<T>& mpCheck) -> int {
            int cnt = 0;
            for (size_t p = 0; p < mpCheck.nPatches(); ++p)
                if (evalPatchMinDetDense(mpCheck.patch(p), 49) <= T(1e-8)) ++cnt;
            return cnt;
        };
        auto repairMicroFolds = [&](gsMultiPatch<T>& mpFix) {
            for (size_t p = 0; p < mpFix.nPatches(); ++p) {
                T d0 = evalPatchMinDetDense(mpFix.patch(p), 49);
                if (d0 <= T(1e-8)) {
                    if (auto* tb = dynamic_cast<gsTensorBSpline<2, T>*>(&mpFix.patch(p))) {
                        gsTensorBSpline<2, T> bk(*tb);
                        quickHarmonicPatchRepair(*tb, T(1e-7));
                        if (evalPatchMinDetDense(*tb, 49) < d0) *tb = bk;
                    }
                }
            }
        };

        repairMicroFolds(matMp);

        // 1. Untangle interior CPs and curve shared MAT interfaces + slide boundary footpoints first
        //    (monotonic dense-injectivity rollback guards so coarse 17x17 grids never fold a valid multi-span patch!)
        {
            gsMultiPatch<T> preWin(matMp);
            int negBefore = countNegDense(matMp);
            optimizeWinslowBarrier(matMp);
            repairMicroFolds(matMp);
            if (countNegDense(matMp) > negBefore)
                matMp = give(preWin);
        }
        if (matMp.nPatches() >= 2) {
            gsMultiPatch<T> preCpl(matMp);
            int negBefore = countNegDense(matMp);
            optimizeCoupledInterfacesAndInteriors(matMp, T(1e-7), 24, true);
            repairMicroFolds(matMp);
            if (countNegDense(matMp) > negBefore)
                matMp = give(preCpl);
        }

        // 2. Collapse redundant 2-wide medial centerlines & short slivers (aspect-ratio guarded, up to 2 passes)
        for (int rPass = 0; rPass < 2 && matMp.nPatches() > 2; ++rPass) {
            int negBefore = countNegDense(matMp);
            gsMultiPatch<T> redMp(matMp);
            if (reducePolychordStrips(redMp) && redMp.nPatches() > 0 && redMp.nPatches() < matMp.nPatches()) {
                redMp.computeTopology(T(1e-4), false, false);
                repairMicroFolds(redMp);
                {
                    gsMultiPatch<T> preRedOpt(redMp);
                    int negRed = countNegDense(redMp);
                    optimizeWinslowBarrier(redMp);
                    if (redMp.nPatches() >= 2)
                        optimizeCoupledInterfacesAndInteriors(redMp, T(1e-7), 18, true);
                    repairMicroFolds(redMp);
                    if (countNegDense(redMp) > negRed)
                        redMp = give(preRedOpt);
                }
                if (countNegDense(redMp) <= negBefore)
                    matMp = give(redMp);
                else
                    break;
            } else {
                break;
            }
        }

        // 3. Final dense certification
        gsASG1LocalFitter<T>::unifyPatchOrientations(matMp);
        matMp.computeTopology(T(1e-4), false, false);
        repairMicroFolds(matMp);
        if (static_cast<int>(matMp.nPatches()) > maxPatchesCap)
            continue;

        std::string vErr;
        if (validateMultiPatch(m_domain, matMp, T(1e-4), &vErr, minHoleAngTol, effMaxAng)) {
            bool denseOk = true;
            for (size_t p = 0; p < matMp.nPatches(); ++p) {
                if (evalPatchMinDetDense(matMp.patch(p), 49) <= T(1e-8)) {
                    denseOk = false;
                    break;
                }
            }
            if (!denseOk) continue;

            SplittingCandidate cand;
            cand.numPatches = static_cast<int>(matMp.nPatches());
            cand.distortionEnergy = evaluateDistortionEnergy(matMp, &cand.minJacobianDet,
                                                             &cand.minCornerAngleDeg, &cand.maxCornerAngleDeg);
            cand.isValid = (cand.minJacobianDet > T(1e-8) && cand.maxCornerAngleDeg <= effMaxAng &&
                            cand.numPatches <= maxPatchesCap);
            if (cand.isValid) {
                cand.multiPatch = give(matMp);
                cand.methodUsed = "MedialAxisTransform";
                if (!bestCand.isValid ||
                    (m_weightPatches * cand.numPatches + m_weightDistortion * cand.distortionEnergy <
                     m_weightPatches * bestCand.numPatches + m_weightDistortion * bestCand.distortionEnergy)) {
                    bestCand = give(cand);
                }
                if (bestCand.isValid && bestCand.numPatches <= std::max(12, 6 * m_domain.numHoles() + 4))
                    break;
            }
        }
    }
    return bestCand;
}

// ============================================================================
// Approach 4: Field-Guided Medial Macro-Pairing ("Merged" 5-Module Algorithm)
// ============================================================================

template <class T>
typename gsPatchSplitter<T>::SplittingCandidate
gsPatchSplitter<T>::computeFieldGuidedMedialSplitting()
{
    using namespace internal_splitter;
    SplittingCandidate bestCand;
    bestCand.methodUsed = "FieldGuidedMedialPairing";
    bestCand.isValid = false;

    if (m_domain.numLoops() == 0) return bestCand;

    const bool dbg = (std::getenv("GISMO_MERGED_DEBUG") != nullptr);
    const T domMinAng = getDomainMinCornerAngleDeg();
    const T minAngTol = std::min(static_cast<T>(m_options.getReal("MinCornerAngleDeg")), T(0.85) * domMinAng);
    const T effMinAng = (m_domain.numHoles() > 0)
                      ? std::min(minAngTol, std::max(T(4.0), minAngTol - T(4.0)))
                      : std::min(minAngTol, std::max(T(4.0), minAngTol - T(3.5)));
    const T maxAngTol = static_cast<T>(m_options.getReal("MaxCornerAngleDeg"));
    const T effMaxAng = std::max(maxAngTol, T(164.9));
    const int maxPatchesCap = std::max(1, static_cast<int>(m_options.getInt("MaxPatches")));
    const T smoothJointTol = static_cast<T>(m_options.getReal("SmoothJointTolDeg"));

    auto countNegDense = [&](const gsMultiPatch<T>& mpCheck) -> int {
        int cnt = 0;
        for (size_t p = 0; p < mpCheck.nPatches(); ++p)
            if (!isPatchDetSingleSign(mpCheck.patch(p), 21, T(1e-7)) ||
                evalPatchMinDetDense(mpCheck.patch(p), 49) <= T(1e-8))
                ++cnt;
        return cnt;
    };
    auto repairMicroFolds = [&](gsMultiPatch<T>& mpFix) {
        for (size_t p = 0; p < mpFix.nPatches(); ++p) {
            T d0 = evalPatchMinDetDense(mpFix.patch(p), 49);
            if (d0 <= T(1e-8)) {
                if (auto* tb = dynamic_cast<gsTensorBSpline<2, T>*>(&mpFix.patch(p))) {
                    gsTensorBSpline<2, T> bk(*tb);
                    quickHarmonicPatchRepair(*tb, T(1e-7));
                    T d1 = evalPatchMinDetDense(*tb, 49);
                    if (d1 < d0) { *tb = bk; }
                }
            }
        }
    };
    auto checkConforming = [](gsMultiPatch<T>& m, size_t maxAllowedBnd = std::numeric_limits<size_t>::max()) -> bool {
        m.computeTopology(T(1e-4), false, false);
        if (2 * m.nInterfaces() + m.nBoundary() != 4 * m.nPatches()) return false;
        if (m.nBoundary() > maxAllowedBnd) return false;
        return true;
    };

    struct MergedCfg { T hFrac; T thetaPrune; T jMerge; T alphaMAT; };
    const MergedCfg configs[4] = {
        {T(1.0) / T(24.0), T(45.0), T(0.32), T(2.0)},
        {T(1.0) / T(30.0), T(38.0), T(0.28), T(1.2)},
        {T(1.0) / T(18.0), T(52.0), T(0.45), T(2.8)},
        {T(1.0) / T(36.0), T(34.0), T(0.55), T(1.5)}
    };

    for (int ci = 0; ci < 4; ++ci) {
        gsMultiPatch<T> mp;
        if (!buildFieldGuidedMedialMultiPatch(
                m_domain, m_degree, smoothJointTol, mp,
                configs[ci].hFrac, configs[ci].thetaPrune, configs[ci].jMerge, configs[ci].alphaMAT))
            continue;
        if (mp.nPatches() == 0 || static_cast<int>(mp.nPatches()) > maxPatchesCap + 40)
            continue;

        gsASG1LocalFitter<T>::unifyPatchOrientations(mp);
        if (!checkConforming(mp)) {
            if (dbg) std::cerr << "[MERGED] cfg=" << ci << " P=" << mp.nPatches() << " FAIL non-conforming\n";
            continue;
        }
        harmonizeMultiPatchKnots(mp);
        const size_t initBndSides = mp.nBoundary();
        repairMicroFolds(mp);

        // 1. Interior vertex smoothing + corner-angle relaxation
        if (mp.nPatches() >= 2 && mp.nPatches() <= 160) {
            gsMultiPatch<T> preV(mp);
            int negB = countNegDense(mp);
            size_t bndB = mp.nBoundary();
            relaxInteriorVertices(mp, 5);
            relaxBadCornerAngles(mp, effMinAng, effMaxAng, 3);
            repairMicroFolds(mp);
            if (!checkConforming(mp, bndB) || countNegDense(mp) > negB)
                mp = give(preV);
        }

        // 2. Winslow barrier + Coupled moving-interface Winslow + CAD boundary split sliding
        {
            gsMultiPatch<T> preW(mp);
            int negB = countNegDense(mp);
            size_t bndB = mp.nBoundary();
            optimizeWinslowBarrier(mp);
            repairMicroFolds(mp);
            if (!checkConforming(mp, bndB) || countNegDense(mp) > negB)
                mp = give(preW);
        }
        if (mp.nPatches() >= 2 && mp.nPatches() <= 140) {
            gsMultiPatch<T> preC(mp);
            int negB = countNegDense(mp);
            size_t bndB = mp.nBoundary();
            bool hasSkewed = false;
            for (size_t p = 0; p < mp.nPatches(); ++p) {
                T angs[4], mnA = 0, mxA = 180; bool deg = false;
                evalPatchCornerAnglesDeg(mp.patch(p), angs, mnA, mxA, deg);
                if (deg || mnA < T(30.0) || mxA > T(140.0)) { hasSkewed = true; break; }
            }
            optimizeCoupledInterfacesAndInteriors(mp, T(1e-7), hasSkewed ? 36 : 24, true);
            repairMicroFolds(mp);
            if (!checkConforming(mp, bndB) || countNegDense(mp) > negB)
                mp = give(preC);
        }

        // 2b. Polychord strip collapse (up to 2 passes) with knot harmonization & post-relaxation
        for (int rPass = 0; rPass < 2 && mp.nPatches() > 2; ++rPass) {
            int negB = countNegDense(mp);
            size_t bndB = mp.nBoundary();
            gsMultiPatch<T> redMp(mp);
            if (reducePolychordStrips(redMp) && redMp.nPatches() > 0 && redMp.nPatches() < mp.nPatches()) {
                redMp.computeTopology(T(1e-4), false, false);
                harmonizeMultiPatchKnots(redMp);
                repairMicroFolds(redMp);
                {
                    gsMultiPatch<T> preR(redMp);
                    int negR = countNegDense(redMp);
                    size_t bndR = redMp.nBoundary();
                    optimizeWinslowBarrier(redMp);
                    if (redMp.nPatches() >= 2 && redMp.nPatches() <= 140)
                        optimizeCoupledInterfacesAndInteriors(redMp, T(1e-7), 16, true);
                    repairMicroFolds(redMp);
                    if (!checkConforming(redMp, bndR) || countNegDense(redMp) > negR)
                        redMp = give(preR);
                }
                if (checkConforming(redMp, bndB) && countNegDense(redMp) <= negB)
                    mp = give(redMp);
                else
                    break;
            } else {
                break;
            }
        }

        // 2c. Final corner-angle & interior-vertex polishing + Stage 4 Topological O-Ring Pruning & Absorption!
        if (mp.nPatches() >= 2 && mp.nPatches() <= 160) {
            gsMultiPatch<T> preP(mp);
            int negB = countNegDense(mp);
            size_t bndB = mp.nBoundary();
            relaxInteriorVertices(mp, 5);
            relaxBadCornerAngles(mp, effMinAng, effMaxAng, 6);
            repairMicroFolds(mp);
            if (!checkConforming(mp, bndB) || countNegDense(mp) > negB)
                mp = give(preP);
        }
        if (mp.nPatches() >= 5) {
            gsMultiPatch<T> prePrune(mp);
            int negB = countNegDense(mp);
            int nPruned = pruneParasiticORings(m_domain, mp, m_degree, smoothJointTol, effMinAng, effMaxAng, dbg);
            if (nPruned > 0) {
                size_t prunedBnd = mp.nBoundary();
                harmonizeMultiPatchKnots(mp);
                repairMicroFolds(mp);
                if (mp.nPatches() >= 2 && mp.nPatches() <= 140)
                    optimizeCoupledInterfacesAndInteriors(mp, T(1e-7), 18, true);
                relaxInteriorVertices(mp, 4);
                relaxBadCornerAngles(mp, effMinAng, effMaxAng, 5);
                repairMicroFolds(mp);
                if (!checkConforming(mp, prunedBnd) || countNegDense(mp) > negB ||
                    !verifyBoundaryMatch(m_domain, mp, T(1e-4))) {
                    mp = give(prePrune);
                }
            }
        }
        (void)initBndSides;

        // 2d. Targeted Local 1-to-5 O-Grid / Sharp-Kite Rescue ONLY on the 1..14 individual stubborn patches
        //     that remain folded (matching validateMultiPatch's exact 1e-7 threshold!) or have out-of-bounds angles.
        const std::vector<std::vector<gsPoint2D<T>>> emptyPolys;
        for (int oRound = 0; oRound < 2; ++oRound) {
            auto collectBad = [&](const gsMultiPatch<T>& m, std::vector<size_t>& badOut, int& nNegOut, int& nAngOut) {
                badOut.clear(); nNegOut = 0; nAngOut = 0;
                for (size_t p = 0; p < m.nPatches(); ++p) {
                    bool isNeg = (!isPatchDetSingleSign(m.patch(p), 21, T(1e-7)) ||
                                  evalPatchMinDetDense(m.patch(p), 49) <= T(1e-8));
                    T angs[4], mnA = 0, mxA = 180; bool deg = false;
                    evalPatchCornerAnglesDeg(m.patch(p), angs, mnA, mxA, deg);
                    bool isBadA = (deg || mxA > effMaxAng || mnA < effMinAng);
                    if (isNeg) ++nNegOut;
                    if (isBadA) ++nAngOut;
                    if (isNeg || isBadA) badOut.push_back(p);
                }
            };

            std::vector<size_t> badIds;
            int oldNeg = 0, oldAng = 0;
            collectBad(mp, badIds, oldNeg, oldAng);
            if (badIds.empty() || badIds.size() > 14) break;
            if (static_cast<int>(mp.nPatches() + 5 * badIds.size()) > maxPatchesCap) break;

            std::set<size_t> badSet(badIds.begin(), badIds.end());
            gsMultiPatch<T> trialMp;
            int nReplaced = 0;
            for (size_t p = 0; p < mp.nPatches(); ++p) {
                if (!badSet.count(p)) {
                    trialMp.addPatch(mp.patch(p).clone());
                    continue;
                }
                typename gsGeometry<T>::uPtr sS = mp.patch(p).boundary(boundary::south);
                typename gsGeometry<T>::uPtr sE = mp.patch(p).boundary(boundary::east);
                typename gsGeometry<T>::uPtr sN = mp.patch(p).boundary(boundary::north);
                typename gsGeometry<T>::uPtr sW = mp.patch(p).boundary(boundary::west);
                auto* bS = dynamic_cast<gsBSpline<T>*>(sS.get());
                auto* bE = dynamic_cast<gsBSpline<T>*>(sE.get());
                auto* bN = dynamic_cast<gsBSpline<T>*>(sN.get());
                auto* bW = dynamic_cast<gsBSpline<T>*>(sW.get());
                bool didSub = false;
                if (bS && bE && bN && bW) {
                    std::array<CurveSeg<T>, 4> b4 = {
                        CurveSeg<T>(*bS), CurveSeg<T>(*bE),
                        CurveSeg<T>(*bN).reversed(), CurveSeg<T>(*bW).reversed()
                    };
                    std::vector<gsPoint2D<T>> poly;
                    poly.reserve(64);
                    for (int k = 0; k < 4; ++k)
                        for (int s = 0; s < 16; ++s)
                            poly.push_back(evalCurvePoint(*b4[k].bsp, static_cast<T>(s) / T(16)));
                    T area2 = T(0);
                    for (size_t i = 0; i < poly.size(); ++i)
                        area2 += poly[i].x * poly[(i + 1) % poly.size()].y - poly[(i + 1) % poly.size()].x * poly[i].y;
                    if (area2 < T(0)) {
                        b4 = {b4[3].reversed(), b4[2].reversed(), b4[1].reversed(), b4[0].reversed()};
                    }
                    // Use exact patch boundary `capPoly` (emptyPolys avoids false rejections from coarse CAD polygonization!)
                    std::vector<std::array<CurveSeg<T>, 4>> ogQuads;
                    bool ogRes = emitRobustOGrid5(b4, m_degree, emptyPolys, ogQuads, true, std::max(T(4.5), effMinAng));
                    if (ogRes) {
                        std::vector<std::unique_ptr<gsTensorBSpline<2, T>>> built;
                        bool allOk = !ogQuads.empty();
                        for (const auto& oq : ogQuads) {
                            std::unique_ptr<gsTensorBSpline<2, T>> np(
                                makeCoonsFromCCWCycle(oq[0], oq[1], oq[2], oq[3], m_degree));
                            if (!np) { allOk = false; break; }
                            if (evalPatchMinDetDense(*np, 25) <= T(1e-7))
                                quickHarmonicPatchRepair(*np, T(1e-7));
                            if ((!isPatchDetSingleSign(*np, 21, T(1e-7)) || evalPatchMinDetDense(*np, 25) <= T(1e-8)) &&
                                np->knots(0).unique().size() < 3 && np->knots(1).unique().size() < 3) {
                                np->insertKnot(T(0.5), 0, 1);
                                np->insertKnot(T(0.5), 1, 1);
                                quickHarmonicPatchRepair(*np, T(1e-7));
                                if (!isPatchDetSingleSign(*np, 21, T(1e-7)) || evalPatchMinDetDense(*np, 25) <= T(1e-8)) {
                                    np->insertKnot(T(0.25), 0, 1);
                                    np->insertKnot(T(0.75), 0, 1);
                                    np->insertKnot(T(0.25), 1, 1);
                                    np->insertKnot(T(0.75), 1, 1);
                                    quickHarmonicPatchRepair(*np, T(1e-7));
                                }
                                if (!isPatchDetSingleSign(*np, 21, T(1e-7)) || evalPatchMinDetDense(*np, 25) <= T(1e-8)) {
                                    gsMultiPatch<T> tmpO;
                                    tmpO.addPatch(np->clone());
                                    if (optimizeWinslowBarrier(tmpO, T(1e-7))) {
                                        if (auto* fx = dynamic_cast<gsTensorBSpline<2, T>*>(&tmpO.patch(0)))
                                            *np = *fx;
                                    }
                                }
                            }
                            if (!isPatchDetSingleSign(*np, 21, T(1e-7)) || evalPatchMinDetDense(*np, 25) <= T(1e-8)) {
                                allOk = false; break;
                            }
                            built.push_back(std::move(np));
                        }
                        if (allOk) {
                            for (auto& np : built)
                                trialMp.addPatch(typename gsGeometry<T>::uPtr(np.release()));
                            didSub = true;
                            ++nReplaced;
                        } else if (dbg) {
                            std::cerr << "[MERGED] 2d p=" << p << " emitOk(n=" << ogQuads.size() << ") but Coons allOk=0\n";
                        }
                    } else if (dbg) {
                        T angs[4], mnA = 0, mxA = 180; bool deg = false;
                        evalPatchCornerAnglesDeg(mp.patch(p), angs, mnA, mxA, deg);
                        std::cerr << "[MERGED] 2d p=" << p << " emitRobustOGrid5=0 minDet=" << evalPatchMinDetDense(mp.patch(p), 49)
                                  << " angs=[" << angs[0] << "," << angs[1] << "," << angs[2] << "," << angs[3] << "] area2=" << area2 << "\n";
                    }
                }
                if (!didSub)
                    trialMp.addPatch(mp.patch(p).clone());
            }
            if (nReplaced == 0) {
                if (dbg) std::cerr << "[MERGED] 2d r=" << oRound << " bad=" << badIds.size() << " (neg=" << oldNeg << " ang=" << oldAng << ") replaced=0\n";
                break;
            }

            size_t preBnd2d = mp.nBoundary();
            gsASG1LocalFitter<T>::unifyPatchOrientations(trialMp);
            if (!checkConforming(trialMp, preBnd2d)) break;
            harmonizeMultiPatchKnots(trialMp);
            repairMicroFolds(trialMp);
            {
                gsMultiPatch<T> preO(trialMp);
                int nPre = countNegDense(trialMp);
                optimizeWinslowBarrier(trialMp);
                relaxInteriorVertices(trialMp, 5);
                relaxBadCornerAngles(trialMp, effMinAng, effMaxAng, 8);
                repairMicroFolds(trialMp);
                if (!checkConforming(trialMp, preBnd2d) || countNegDense(trialMp) > nPre)
                    trialMp = give(preO);
            }
            std::vector<size_t> newBad;
            int newNeg = 0, newAng = 0;
            collectBad(trialMp, newBad, newNeg, newAng);
            if (dbg) {
                std::cerr << "[MERGED] 2d r=" << oRound << " bad=" << badIds.size()
                          << " (" << oldNeg << "," << oldAng << ") -> repl=" << nReplaced
                          << " new=(" << newNeg << "," << newAng << ") P=" << trialMp.nPatches() << "\n";
            }
            if (checkConforming(trialMp, preBnd2d) &&
                (newNeg < oldNeg || (newNeg == oldNeg && newAng < oldAng))) {
                mp = give(trialMp);
                if (newNeg == 0 && newAng == 0) break;
            } else {
                break;
            }
        }

        // 3. Optional interior DOF rescue for strongly curved macro-patches (up to 140 patches)
        if (countNegDense(mp) > 0 && mp.nPatches() <= 140) {
            gsMultiPatch<T> refMp(mp);
            int negPre = countNegDense(refMp);
            size_t bndPre3 = mp.nBoundary();
            for (size_t p = 0; p < refMp.nPatches(); ++p) {
                if (auto* tb = dynamic_cast<gsTensorBSpline<2, T>*>(&refMp.patch(p))) {
                    safeNormalizeKnotVector01(tb->knots(0));
                    safeNormalizeKnotVector01(tb->knots(1));
                    for (int d = 0; d < 2; ++d) {
                        if (tb->knots(d).unique().size() < 4) {
                            tb->insertKnot(T(1) / T(3), d, 1);
                            tb->insertKnot(T(2) / T(3), d, 1);
                        }
                    }
                }
            }
            harmonizeMultiPatchKnots(refMp);
            repairMicroFolds(refMp);
            optimizeWinslowBarrier(refMp);
            repairMicroFolds(refMp);
            if (refMp.nPatches() >= 2) {
                optimizeCoupledInterfacesAndInteriors(refMp, T(1e-7), 18, true);
                relaxInteriorVertices(refMp, 5);
                relaxBadCornerAngles(refMp, effMinAng, effMaxAng, 10);
                repairMicroFolds(refMp);
            }
            gsASG1LocalFitter<T>::unifyPatchOrientations(refMp);
            if (checkConforming(refMp, bndPre3) && countNegDense(refMp) < negPre)
                mp = give(refMp);
        }

        gsASG1LocalFitter<T>::unifyPatchOrientations(mp);
        mp.computeTopology(T(1e-4), false, false);
        repairMicroFolds(mp);
        if (static_cast<int>(mp.nPatches()) > maxPatchesCap) continue;

        std::string vErr;
        if (!validateMultiPatch(m_domain, mp, T(1e-4), &vErr, effMinAng, effMaxAng)) {
            // Targeted extra corner-angle polish if all patches are dense-positive and only a corner angle missed!
            if (countNegDense(mp) == 0 && vErr.find("corner tangent angle") != std::string::npos &&
                mp.nPatches() >= 2 && mp.nPatches() <= 160) {
                gsMultiPatch<T> polMp(mp);
                relaxInteriorVertices(polMp, 6);
                relaxBadCornerAngles(polMp, effMinAng, effMaxAng, 10);
                repairMicroFolds(polMp);
                if (checkConforming(polMp) && countNegDense(polMp) == 0) {
                    std::string vErr2;
                    if (validateMultiPatch(m_domain, polMp, T(1e-4), &vErr2, effMinAng, effMaxAng)) {
                        mp = give(polMp);
                        vErr.clear();
                    } else {
                        vErr = vErr2;
                    }
                }
            }
            if (!vErr.empty()) {
                if (dbg) std::cerr << "[MERGED] cfg=" << ci << " P=" << mp.nPatches() << " FAIL validate: " << vErr << "\n";
                continue;
            }
        }
        if (countNegDense(mp) > 0) {
            if (dbg) std::cerr << "[MERGED] cfg=" << ci << " P=" << mp.nPatches() << " FAIL denseNeg=" << countNegDense(mp) << "\n";
            continue;
        }

        T minDet = 0, minA = 0, maxA = 180;
        T energy = evaluateDistortionEnergy(mp, &minDet, &minA, &maxA);
        if (minDet <= T(1e-8) || maxA > effMaxAng) continue;

        SplittingCandidate cand;
        cand.multiPatch = give(mp);
        cand.numPatches = static_cast<int>(cand.multiPatch.nPatches());
        cand.minJacobianDet = minDet;
        cand.distortionEnergy = energy;
        cand.minCornerAngleDeg = minA;
        cand.maxCornerAngleDeg = maxA;
        cand.methodUsed = "FieldGuidedMedialPairing";
        cand.isValid = true;

        if (dbg) {
            std::cerr << "[MERGED] cfg=" << ci << " VALID! P=" << cand.numPatches
                      << " minDet=" << minDet << " ang=[" << minA << "," << maxA << "]\n";
        }

        if (!bestCand.isValid ||
            (m_weightPatches * cand.numPatches + m_weightDistortion * cand.distortionEnergy <
             m_weightPatches * bestCand.numPatches + m_weightDistortion * bestCand.distortionEnergy)) {
            bestCand = give(cand);
        }

        // Early exit if already compact and well-conditioned
        int compactCap = std::max(10, 5 * m_domain.numHoles() + 6);
        if (bestCand.isValid && bestCand.numPatches <= compactCap &&
            bestCand.maxCornerAngleDeg <= T(148.0) && bestCand.minCornerAngleDeg >= T(20.0))
            break;
    }

    return bestCand;
}

template <class T>
typename gsPatchSplitter<T>::SplittingCandidate
gsPatchSplitter<T>::computeMultiHoleStructuralTemplates()
{
    SplittingCandidate candidate;
    candidate.methodUsed = "TemplateDecomposition";
    candidate.isValid = false;

    if (m_domain.numLoops() <= 1) {
        return candidate;
    }

    const auto& outer = m_domain.outer();
    int nHoles = m_domain.numHoles();
    int nOuterCurves = outer.size();

    const T domMinAng = getDomainMinCornerAngleDeg();
    const T minAngTol = std::min(static_cast<T>(m_options.getReal("MinCornerAngleDeg")), T(0.85) * domMinAng);
    const T minHoleAngTol = std::min(minAngTol, T(4.0));
    const T maxAngTol = static_cast<T>(m_options.getReal("MaxCornerAngleDeg"));
    auto tryFinalizeHole = [&](SplittingCandidate& cand) -> bool {
        const bool dbg = (std::getenv("GISMO_SPLITTER_DEBUG") != nullptr);
        if (cand.multiPatch.nPatches() == 0) { if (dbg) std::cerr << "[DBG] empty multiPatch\n"; return false; }
        cand.multiPatch.computeTopology(T(1e-4), false, false);
        optimizeWinslowBarrier(cand.multiPatch);
        gsASG1LocalFitter<T>::unifyPatchOrientations(cand.multiPatch);
        std::string vErr;
        if (!validateMultiPatch(m_domain, cand.multiPatch, T(1e-4), &vErr, minHoleAngTol, std::max(maxAngTol, T(164.9)))) {
            if (dbg) {
                T md=0, ma=0, xa=180;
                evaluateDistortionEnergy(cand.multiPatch, &md, &ma, &xa);
                std::cerr << "[DBG] validateMultiPatch failed (P=" << cand.multiPatch.nPatches()
                          << ", IF=" << cand.multiPatch.topology().nInterfaces()
                          << ", BND=" << cand.multiPatch.topology().nBoundary()
                          << ", domMin=" << domMinAng << ", tol=[" << minHoleAngTol << "," << std::max(maxAngTol, T(164.9)) << "]"
                          << ", minDet=" << md << ", ang=[" << ma << "," << xa << "]): " << vErr << "\n";
            }
            cand.multiPatch.clear();
            cand.isValid = false;
            return false;
        }
        T minPos20 = std::numeric_limits<T>::max();
        for (size_t p = 0; p < cand.multiPatch.nPatches(); ++p) {
            T pMin = 0, pMax = 0, pAbs = 0, pEn = 0;
            internal_splitter::evalPatchJacobianAndWinslow(cand.multiPatch.patch(p), 20, pMin, pMax, pAbs, pEn);
            minPos20 = std::min(minPos20, pMin);
        }
        cand.numPatches = static_cast<int>(cand.multiPatch.nPatches());
        cand.distortionEnergy = evaluateDistortionEnergy(cand.multiPatch, &cand.minJacobianDet,
                                                         &cand.minCornerAngleDeg, &cand.maxCornerAngleDeg);
        cand.isValid = (minPos20 > T(1e-8) && cand.minJacobianDet > T(1e-8) && cand.maxCornerAngleDeg < std::max(maxAngTol, T(164.9)));
        if (!cand.isValid) {
            if (dbg) std::cerr << "[DBG] post-check failed (P=" << cand.numPatches
                               << " minPos20=" << minPos20 << " minDet=" << cand.minJacobianDet
                               << " ang=[" << cand.minCornerAngleDeg << "," << cand.maxCornerAngleDeg << "])\n";
            cand.multiPatch.clear();
            return false;
        }
        return true;
    };

    using internal_splitter::CurveSeg;
    std::vector<CurveSeg<T>> oCCW = internal_splitter::orientLoopSegments(outer, m_degree, true);

    auto addBlockCCW = [&](const CurveSeg<T>& e01, const CurveSeg<T>& e12,
                           const CurveSeg<T>& e23, const CurveSeg<T>& e30,
                           const std::vector<const gsCurveLoop<T>*>& blockHoles) -> bool {
        return internal_splitter::build4SidedMultiHoleGrid(
            e01, e12, e23.reversed(), e30.reversed(), blockHoles, candidate.multiPatch, m_degree);
    };

    // ------------------------------------------------------------------------
    // Case 2: 4 Outer Curves, Any Number of Holes (1, 2, 3, ..., N)
    // ------------------------------------------------------------------------
    if (nOuterCurves == 4) {
        if (nHoles == 1 && m_domain.loop(1).size() == 4) {
            internal_splitter::buildQuadRingByProximity(outer, m_domain.loop(1), candidate.multiPatch, m_degree);
            if (tryFinalizeHole(candidate)) return candidate;
        }
        std::vector<const gsCurveLoop<T>*> allHoles;
        for (int h = 0; h < nHoles; ++h) allHoles.push_back(&m_domain.loop(h + 1));
        for (int rot = 0; rot < 2; ++rot) {
            candidate.multiPatch.clear();
            if (addBlockCCW(oCCW[rot], oCCW[(rot+1)%4], oCCW[(rot+2)%4], oCCW[(rot+3)%4], allHoles)) {
                if (tryFinalizeHole(candidate)) return candidate;
            }
        }
        if (nHoles >= 2) {
            candidate.multiPatch.clear();
            if (internal_splitter::build4SidedInnerAxisGrid(oCCW, allHoles, candidate.multiPatch, m_degree)) {
                if (tryFinalizeHole(candidate)) return candidate;
            }
        }
        candidate.multiPatch.clear();
    }

    // Generalized rotation-invariant recursive cell decomposition (handles arbitrary convex/non-convex polygons,
    // rotated/bent/swirled shapes, and 1..N holes without axis-aligned hardcoded templates).
    // Includes Phase 3 per-cell hybrid leaf fallback: keeps sibling 4-patch hole rings and meshes
    // only a failing leaf cell with a conforming zero-cut-split sub-domain triangulation.
    {
        int fallbackCalls = 0;
        auto leafHybridFallback = [&](const std::vector<CurveSeg<T>>& leafCCW,
                                      const std::vector<bool>& leafIsCut,
                                      const gsCurveLoop<T>* optHole,
                                      std::vector<std::unique_ptr<gsGeometry<T>>>& outLeafP) -> bool {
            if (leafCCW.size() < 3 || ++fallbackCalls > 4) return false;
            std::vector<gsCurve<T>*> outerCrvs;
            outerCrvs.reserve(leafCCW.size());
            for (const auto& cs : leafCCW) {
                if (!cs.bsp) {
                    for (auto* p : outerCrvs) delete p;
                    return false;
                }
                outerCrvs.push_back(cs.bsp->clone().release());
            }
            gsCurveLoop<T>* outerLp = new gsCurveLoop<T>(outerCrvs);
            std::vector<gsCurveLoop<T>*> loops;
            loops.push_back(outerLp);
            if (optHole)
                loops.push_back(optHole->clone().release());
            gsPlanarDomain<T> subDom(loops);

            gsPatchSplitter<T> subSplit(subDom);
            subSplit.setDegree(m_degree);
            subSplit.options().setReal("MinCornerAngleDeg", m_options.getReal("MinCornerAngleDeg"));
            subSplit.options().setReal("MaxCornerAngleDeg", std::max(maxAngTol, T(164.9)));
            SplittingCandidate subCand = subSplit.computeMultiHoleTriangulationSplitting(true);
            if (!subCand.isValid || subCand.multiPatch.nPatches() == 0 || subCand.numPatches > 48)
                return false;

            // Verify 100% conformity along every `isCutEdge`: no cut edge may be split into multiple patch sides!
            gsMultiPatch<T> subMp = subCand.multiPatch;
            subMp.computeTopology(T(1e-4), false, false);
            for (size_t ci = 0; ci < leafCCW.size(); ++ci) {
                if (!leafIsCut[ci]) continue;
                gsPoint2D<T> pA = leafCCW[ci].p0(), pB = leafCCW[ci].p1();
                gsMatrix<T> cpr = leafCCW[ci].bsp->parameterRange();
                gsPoint2D<T> pM = internal_splitter::evalCurvePoint(*leafCCW[ci].bsp, T(0.5) * (cpr(0, 0) + cpr(0, 1)));
                bool exactSideFound = false;
                for (auto bit = subMp.bBegin(); bit != subMp.bEnd(); ++bit) {
                    auto bc = subMp.patch(bit->patch).boundary(bit->side());
                    gsMatrix<T> pr = bc->parameterRange(), u(1, 3), ev;
                    u << pr(0, 0), T(0.5) * (pr(0, 0) + pr(0, 1)), pr(0, 1);
                    bc->eval_into(u, ev);
                    gsPoint2D<T> q0(ev(0, 0), ev(1, 0)), qm(ev(0, 1), ev(1, 1)), q1(ev(0, 2), ev(1, 2));
                    bool endMatch = (q0.dist(pA) < T(1e-4) && q1.dist(pB) < T(1e-4)) ||
                                    (q0.dist(pB) < T(1e-4) && q1.dist(pA) < T(1e-4));
                    if (endMatch && qm.dist(pM) < T(2e-4)) {
                        exactSideFound = true;
                        break;
                    }
                }
                if (!exactSideFound) return false;
            }
            outLeafP.clear();
            outLeafP.reserve(subMp.nPatches());
            for (size_t p = 0; p < subMp.nPatches(); ++p)
                outLeafP.emplace_back(subMp.patch(p).clone().release());
            return true;
        };

        (void)leafHybridFallback;
        gsMultiPatch<T> recMp;
        std::string recLabel;
        if (internal_splitter::computeRecursiveCellSplitting<T>(
                m_domain, m_degree, m_weightPatches, m_weightDistortion,
                minHoleAngTol, std::max(maxAngTol, T(164.9)), recMp, recLabel, nullptr)) {
            candidate.multiPatch = give(recMp);
            candidate.methodUsed = recLabel.empty() ? std::string("RecursiveCellTemplate") : recLabel;
            if (tryFinalizeHole(candidate)) return candidate;
            candidate.multiPatch.clear();
        }
    }

    return candidate;
}

// ============================================================================
// Approach 3: Topological Preset Template Decomposition
// ============================================================================

template <class T>
typename gsPatchSplitter<T>::SplittingCandidate
gsPatchSplitter<T>::computeTemplateSplitting()
{
    SplittingCandidate candidate;
    candidate.methodUsed = "TemplateDecomposition";
    candidate.isValid = false;

    if (m_domain.numLoops() == 0) {
        return candidate;
    }

    if (m_domain.numHoles() == 0) {
        candidate = computeHoleFreeStructuralSplitting();
    } else {
        candidate = computeMultiHoleStructuralTemplates();
    }
    if (candidate.isValid) {
        if (candidate.methodUsed.empty() || candidate.methodUsed == "MedialAxisSkeleton")
            candidate.methodUsed = "TemplateDecomposition";
    }
    return candidate;
}

// ============================================================================
// Mode A Helpers (Existing gsMultiPatch<T> -> Coarse Conforming Layout)
// ============================================================================

template <class T>
bool gsPatchSplitter<T>::extractPlanarDomainFromMultiPatch(const gsMultiPatch<T>& mp, gsPlanarDomain<T>& outDom)
{
    using namespace internal_splitter;
    if (mp.nPatches() == 0) return false;

    gsMultiPatch<T> mpCopy(mp);
    mpCopy.computeTopology(T(1e-4), false, false);

    struct BndEdge {
        std::unique_ptr<gsCurve<T>> crv;
        gsPoint2D<T> p0, p1;
        bool used = false;
    };
    std::vector<BndEdge> edges;
    for (auto it = mpCopy.bBegin(); it != mpCopy.bEnd(); ++it) {
        auto g = mpCopy.patch(it->patch).boundary(it->side());
        std::unique_ptr<gsCurve<T>> crv(dynamic_cast<gsCurve<T>*>(g.release()));
        if (!crv) return false;
        gsPoint2D<T> a = startPoint(*crv);
        gsPoint2D<T> b = endPoint(*crv);
        edges.push_back({std::move(crv), a, b, false});
    }
    if (edges.empty()) return false;

    gsMatrix<T> bbox;
    mpCopy.boundingBox(bbox);
    const T diag = std::hypot(bbox(0, 1) - bbox(0, 0), bbox(1, 1) - bbox(1, 0));
    const T tol = std::max(T(1e-6), T(1e-4) * diag);

    std::vector<std::unique_ptr<gsCurveLoop<T>>> loops;
    std::vector<T> loopAreas;

    for (size_t s = 0; s < edges.size(); ++s) {
        if (edges[s].used) continue;
        edges[s].used = true;
        std::vector<gsCurve<T>*> crvs;
        crvs.push_back(edges[s].crv->clone().release());
        gsPoint2D<T> head = edges[s].p0;
        gsPoint2D<T> tail = edges[s].p1;

        while (tail.dist(head) > tol) {
            int best = -1;
            bool rev = false;
            T bestD = tol * T(10);
            for (size_t k = 0; k < edges.size(); ++k) {
                if (edges[k].used) continue;
                T d0 = tail.dist(edges[k].p0);
                T d1 = tail.dist(edges[k].p1);
                if (d0 < bestD) { bestD = d0; best = static_cast<int>(k); rev = false; }
                if (d1 < bestD) { bestD = d1; best = static_cast<int>(k); rev = true; }
            }
            if (best < 0) break;
            edges[best].used = true;
            std::unique_ptr<gsCurve<T>> c(edges[best].crv->clone().release());
            if (rev) {
                safeReverseCurve(*c);
                tail = edges[best].p0;
            } else {
                tail = edges[best].p1;
            }
            crvs.push_back(c.release());
        }
        auto lp = std::unique_ptr<gsCurveLoop<T>>(new gsCurveLoop<T>(crvs));
        auto poly = sampleCurveLoopPoly(*lp, 20);
        T area2 = T(0);
        for (size_t i = 0; i < poly.size(); ++i) {
            const auto& p = poly[i];
            const auto& q = poly[(i + 1) % poly.size()];
            area2 += p.x * q.y - q.x * p.y;
        }
        loopAreas.push_back(std::abs(area2) * T(0.5));
        loops.push_back(std::move(lp));
    }

    if (loops.empty()) return false;
    size_t outerIdx = 0;
    for (size_t i = 1; i < loops.size(); ++i)
        if (loopAreas[i] > loopAreas[outerIdx]) outerIdx = i;

    outDom = gsPlanarDomain<T>(loops[outerIdx].release());
    for (size_t i = 0; i < loops.size(); ++i) {
        if (i == outerIdx) continue;
        outDom.insertHole(loops[i].release());
    }
    return true;
}

template <class T>
bool gsPatchSplitter<T>::splitDegenerateSinglePatchOGrid(const gsMultiPatch<T>& inMp, gsMultiPatch<T>& outMp) const
{
    using namespace internal_splitter;
    if (inMp.nPatches() != 1) return false;

    const auto& p0 = inMp.patch(0);
    std::unique_ptr<gsCurve<T>> cS(dynamic_cast<gsCurve<T>*>(p0.boundary(boundary::south).release()));
    std::unique_ptr<gsCurve<T>> cE(dynamic_cast<gsCurve<T>*>(p0.boundary(boundary::east).release()));
    std::unique_ptr<gsCurve<T>> cN(dynamic_cast<gsCurve<T>*>(p0.boundary(boundary::north).release()));
    std::unique_ptr<gsCurve<T>> cW(dynamic_cast<gsCurve<T>*>(p0.boundary(boundary::west).release()));
    if (!cS || !cE || !cN || !cW) return false;
    safeReverseCurve(*cN);
    safeReverseCurve(*cW);

    int deg = std::max(m_degree, m_options.getInt("Degree"));
    std::vector<CurveSeg<T>> bEdges = {
        CurveSeg<T>(*cS, deg),
        CurveSeg<T>(*cE, deg),
        CurveSeg<T>(*cN, deg),
        CurveSeg<T>(*cW, deg)
    };

    // Ensure CCW ordering of bEdges
    gsPoint2D<T> P[4];
    for (int k = 0; k < 4; ++k) P[k] = bEdges[k].p0();
    std::vector<gsPoint2D<T>> poly;
    for (int k = 0; k < 4; ++k) {
        for (int s = 0; s < 16; ++s) {
            T u = static_cast<T>(s) / T(16);
            poly.push_back(evalCurvePoint(*bEdges[k].bsp, u));
        }
    }
    T area2 = T(0);
    for (size_t i = 0; i < poly.size(); ++i) {
        const auto& a = poly[i];
        const auto& b = poly[(i + 1) % poly.size()];
        area2 += a.x * b.y - b.x * a.y;
    }
    if (area2 < T(0)) {
        std::reverse(bEdges.begin(), bEdges.end());
        for (auto& e : bEdges) e = e.reversed();
        for (int k = 0; k < 4; ++k) P[k] = bEdges[k].p0();
        std::reverse(poly.begin(), poly.end());
    }

    gsPoint2D<T> ctr = polygonCentroid(poly);
    if (!pointInPolyStrict(ctr, poly, T(1e-5))) {
        ctr = gsPoint2D<T>((P[0].x + P[1].x + P[2].x + P[3].x) * T(0.25),
                           (P[0].y + P[1].y + P[2].y + P[3].y) * T(0.25));
    }

    T bestEnergy = std::numeric_limits<T>::max();
    gsMultiPatch<T> bestMp;
    const T alphas[3] = {T(0.45), T(0.35), T(0.55)};
    const T maxAngTol = static_cast<T>(m_options.getReal("MaxCornerAngleDeg"));

    for (T alpha : alphas) {
        gsPoint2D<T> Q[4];
        bool ok = true;
        for (int k = 0; k < 4; ++k) {
            Q[k] = gsPoint2D<T>(ctr.x + alpha * (P[k].x - ctr.x),
                                ctr.y + alpha * (P[k].y - ctr.y));
            if (!pointInPolyStrict(Q[k], poly, T(1e-5))) { ok = false; break; }
        }
        if (!ok || !isStrictlyConvexQuad(Q)) continue;

        CurveSeg<T> qEdge[4], spoke[4];
        for (int k = 0; k < 4; ++k) {
            qEdge[k] = makeStraightSeg(Q[k], Q[(k + 1) % 4], deg);
            spoke[k] = makeStraightSeg(Q[k], P[k], deg);
        }

        gsMultiPatch<T> candMp;
        candMp.addPatch(typename gsGeometry<T>::uPtr(
            makeCoonsFromCCWCycle(qEdge[0], qEdge[1], qEdge[2], qEdge[3], deg)));
        for (int k = 0; k < 4; ++k) {
            int kn = (k + 1) % 4;
            candMp.addPatch(typename gsGeometry<T>::uPtr(
                makeCoonsFromCCWCycle(qEdge[k], spoke[kn], bEdges[k].reversed(), spoke[k].reversed(), deg)));
        }

        gsASG1LocalFitter<T>::unifyPatchOrientations(candMp);
        candMp.computeTopology(T(1e-5), false, false);
        optimizeWinslowBarrier(candMp);
        candMp.computeTopology(T(1e-5), false, false);
        T minDet = 0, minA = 0, maxA = 180;
        T energy = evaluateDistortionEnergy(candMp, &minDet, &minA, &maxA);
        if (minDet > T(1e-7) && maxA <= maxAngTol && energy < bestEnergy) {
            bestEnergy = energy;
            bestMp = give(candMp);
        }
    }

    if (bestMp.nPatches() == 5) {
        outMp = give(bestMp);
        return true;
    }
    return false;
}

template <class T>
typename gsPatchSplitter<T>::SplittingCandidate
gsPatchSplitter<T>::computeModeACoarseLayout()
{
    using namespace internal_splitter;
    SplittingCandidate cand;
    cand.methodUsed = "ModeAPolychord";
    if (m_inputMp.nPatches() == 0) return cand;

    gsMultiPatch<T> mp(m_inputMp);
    gsASG1LocalFitter<T>::unifyPatchOrientations(mp);
    mp.computeTopology(T(1e-4), false, false);

    const T maxAngTol = static_cast<T>(m_options.getReal("MaxCornerAngleDeg"));
    const T minAngTol = static_cast<T>(m_options.getReal("MinCornerAngleDeg"));
    bool hasDegenerateCorners = false;
    bool hasReflexCorners = false;
    bool hasBadJacobian = false;

    for (size_t p = 0; p < mp.nPatches(); ++p) {
        T pMinDet = 0, pMaxDet = 0, pMinAbs = 0, pEn = 0;
        evalPatchJacobianAndWinslow(mp.patch(p), 17, pMinDet, pMaxDet, pMinAbs, pEn);
        if (pMinDet <= T(1e-7) && pMaxDet >= -T(1e-7))
            hasBadJacobian = true;

        T angs[4], minA = 0, maxA = 180;
        bool degen = false;
        evalPatchCornerAnglesDeg(mp.patch(p), angs, minA, maxA, degen);
        if (degen || maxA > maxAngTol || minA < minAngTol)
            hasDegenerateCorners = true;
        if (maxA >= T(200))
            hasReflexCorners = true;
    }

    // 1. Single patch with ~180 deg smooth corners (e.g. 1-patch disk/ellipse) -> 5-patch O-grid
    if (mp.nPatches() == 1 && (hasDegenerateCorners || hasBadJacobian) && !hasReflexCorners) {
        gsMultiPatch<T> ogridMp;
        if (splitDegenerateSinglePatchOGrid(mp, ogridMp)) {
            cand.multiPatch = give(ogridMp);
            cand.numPatches = static_cast<int>(cand.multiPatch.nPatches());
            cand.distortionEnergy = evaluateDistortionEnergy(cand.multiPatch, &cand.minJacobianDet,
                                                             &cand.minCornerAngleDeg, &cand.maxCornerAngleDeg);
            cand.methodUsed = "ModeAButterflyOGrid";
            cand.isValid = (cand.minJacobianDet > T(1e-8));
            return cand;
        }
    }

    // 2. If degenerate/reflex/invalid corners remain, extract planar domain and re-split via Mode B
    if (hasDegenerateCorners || hasReflexCorners || hasBadJacobian) {
        gsPlanarDomain<T> extractedDom;
        if (extractPlanarDomainFromMultiPatch(mp, extractedDom)) {
            m_domain = extractedDom;
            m_hasDomain = true;
            m_didDomainReblock = true;
            if (computeSplitting()) {
                return m_bestCandidate;
            }
        }
    }

    // 3. Otherwise keep existing conforming coarse topology for polychord optimization
    cand.multiPatch = give(mp);
    cand.numPatches = static_cast<int>(cand.multiPatch.nPatches());
    cand.distortionEnergy = evaluateDistortionEnergy(cand.multiPatch, &cand.minJacobianDet,
                                                     &cand.minCornerAngleDeg, &cand.maxCornerAngleDeg);
    cand.isValid = (cand.numPatches > 0 && cand.minJacobianDet > T(1e-8));
    return cand;
}

// ============================================================================
// Automated Winner-Selection Metric, Full 4-Stage Pipeline & Execution
// ============================================================================

template <class T>
gsMultiPatch<T> gsPatchSplitter<T>::compute()
{
    m_degree = std::max(1, m_options.getInt("Degree"));
    int stratOpt = m_options.getInt("Strategy");
    if (stratOpt >= 0 && stratOpt <= 4)
        m_strategy = static_cast<Strategy>(stratOpt);

    // Stage 1: Coarse Conforming Quad Layout (Mode A or Mode B)
    if (m_hasInputMp && m_inputMp.nPatches() > 0 && (!m_hasDomain || m_strategy == Strategy::PolychordOnly)) {
        m_bestCandidate = computeModeACoarseLayout();
    } else if (m_hasDomain && m_domain.numLoops() > 0) {
        computeSplitting();
    } else if (m_hasInputMp && m_inputMp.nPatches() > 0) {
        m_bestCandidate = computeModeACoarseLayout();
    } else {
        m_bestCandidate = SplittingCandidate();
        return gsMultiPatch<T>();
    }

    if (m_bestCandidate.multiPatch.nPatches() == 0)
        return m_bestCandidate.multiPatch;

    gsMultiPatch<T> workMp = m_bestCandidate.multiPatch;

    // Stage 2a: Dual-Sheet / Polychord Strip Reduction (Mode B or Mode A->B reblocking, aiming for minimal patch count)
    // Guarded by dense 49x49 Jacobian check so Stage 2a never degrades a certified Stage 1 layout.
    const bool isDomainModeB = m_hasDomain && (!m_hasInputMp || m_didDomainReblock);
    if (isDomainModeB && workMp.nPatches() > 2) {
        bool stage1DenseOk = true;
        for (size_t p = 0; p < workMp.nPatches(); ++p) {
            if (internal_splitter::evalPatchMinDetDense(workMp.patch(p), 49) <= T(1e-8)) {
                stage1DenseOk = false; break;
            }
        }
        gsMultiPatch<T> preStripMp;
        if (stage1DenseOk) preStripMp = workMp;
        reducePolychordStrips(workMp);
        if (stage1DenseOk) {
            bool brokeStrip = false;
            for (size_t p = 0; p < workMp.nPatches(); ++p) {
                if (internal_splitter::evalPatchMinDetDense(workMp.patch(p), 49) <= T(1e-8)) {
                    brokeStrip = true; break;
                }
            }
            if (brokeStrip) workMp = give(preStripMp);
        }
    }

    // Stage 2b: Global Polychord Aspect-Ratio Optimization (zero hanging nodes)
    // Skip re-splitting Mode B multipatches when already a minimal single quad or >= 4 patches
    // (preserves minimal 4-patch annuli and 5-patch butterfly O-grids).
    const bool keepModeBCoarse =
        isDomainModeB &&
        ((m_domain.numHoles() == 0 && m_domain.outer().size() == 4 && workMp.nPatches() == 1) ||
         workMp.nPatches() >= 4);
    if (m_options.getSwitch("EnablePolychordSplitting") && !keepModeBCoarse) {
        bool preSplitDenseOk = true;
        for (size_t p = 0; p < workMp.nPatches(); ++p) {
            T mn20 = 0, mx20 = 0, ma20 = 0, en20 = 0;
            internal_splitter::evalPatchJacobianAndWinslow(workMp.patch(p), 20, mn20, mx20, ma20, en20);
            if (internal_splitter::evalPatchMinDetDense(workMp.patch(p), 49) <= T(1e-8) || mn20 <= T(1e-8)) {
                preSplitDenseOk = false; break;
            }
        }
        gsMultiPatch<T> preSplitMp;
        if (preSplitDenseOk) preSplitMp = workMp;
        optimizePolychordSplits(workMp);
        if (preSplitDenseOk) {
            bool brokeSplit = false;
            for (size_t p = 0; p < workMp.nPatches(); ++p) {
                T mn20 = 0, mx20 = 0, ma20 = 0, en20 = 0;
                internal_splitter::evalPatchJacobianAndWinslow(workMp.patch(p), 20, mn20, mx20, ma20, en20);
                if (internal_splitter::evalPatchMinDetDense(workMp.patch(p), 49) <= T(1e-8) || mn20 <= T(1e-8)) {
                    brokeSplit = true; break;
                }
            }
            if (brokeSplit) workMp = give(preSplitMp);
        }
    }

    // Stage 3: Knot Harmonization, Relaxation & Bijective Parameterization
    relaxAndParameterize(workMp, true);

    // Stage 4: ArgyrisSplines Post-Processing & Verification
    formatForArgyris(workMp);

    // Update cached candidate metrics
    m_bestCandidate.multiPatch = workMp;
    m_bestCandidate.numPatches = static_cast<int>(workMp.nPatches());
    m_bestCandidate.distortionEnergy = evaluateDistortionEnergy(
        workMp, &m_bestCandidate.minJacobianDet,
        &m_bestCandidate.minCornerAngleDeg, &m_bestCandidate.maxCornerAngleDeg);
    m_bestCandidate.isValid = (m_bestCandidate.numPatches > 0 && m_bestCandidate.minJacobianDet > T(1e-9));

    return m_bestCandidate.multiPatch;
}

template <class T>
bool gsPatchSplitter<T>::computeSplitting()
{
    if (m_hasInputMp && !m_hasDomain) {
        m_bestCandidate = computeModeACoarseLayout();
        if (m_bestCandidate.multiPatch.nPatches() > 0)
            ensureExactBoundaryConformity(m_bestCandidate.multiPatch);
        return m_bestCandidate.isValid;
    }

    if (m_strategy == Strategy::PolychordOnly && m_hasInputMp) {
        m_bestCandidate = computeModeACoarseLayout();
        if (m_bestCandidate.isValid) {
            optimizePolychordSplits(m_bestCandidate.multiPatch);
            relaxAndParameterize(m_bestCandidate.multiPatch, true);
            m_bestCandidate.numPatches = static_cast<int>(m_bestCandidate.multiPatch.nPatches());
        }
        return m_bestCandidate.isValid;
    }

    if (m_strategy == Strategy::FieldGuidedMedialPairing) {
        m_bestCandidate = computeFieldGuidedMedialSplitting();
    } else if (m_strategy == Strategy::HarmonicCrossField) {
        m_bestCandidate = computeCrossFieldSplitting();
        if (!m_bestCandidate.isValid) {
            m_bestCandidate = computeMedialAxisSplitting();
        }
    } else if (m_strategy == Strategy::MedialAxisSkeleton) {
        m_bestCandidate = computeMedialAxisSplitting();
        if (!m_bestCandidate.isValid) {
            m_bestCandidate = computeTemplateSplitting();
        }
    } else if (m_strategy == Strategy::TemplateDecomposition) {
        m_bestCandidate = computeTemplateSplitting();
    } else {
        // Strategy::AutomaticBest:
        // Staged cascade: fast structural/recursive templates first with early-exit on compact layouts,
        // falling back to MAT / Cross-Field / Merged / CDT triangulation only if needed.
        SplittingCandidate cfCand;
        SplittingCandidate matCand;
        SplittingCandidate mergedCand;
        SplittingCandidate tmplCand;

        const int nHoles = m_domain.numHoles();

        // Fast-path: convexish 4-corner hole-free quad -> keep as 1 patch cleanly
        if (nHoles == 0 && m_domain.outer().size() == 4) {
            matCand = computeMedialAxisSplitting();
            if (matCand.isValid && matCand.numPatches == 1) {
                m_bestCandidate = matCand;
                return true;
            }
        }

        if (nHoles > 0) {
            mergedCand = computeFieldGuidedMedialSplitting();
            cfCand = computeCrossFieldSplitting();
            matCand = computeMedialAxisSplitting();
        } else {
            tmplCand = computeHoleFreeStructuralSplitting();
            if (!matCand.isValid)
                matCand = computeMedialAxisSplitting();
            bool haveCompact =
                (tmplCand.isValid && tmplCand.numPatches <= 6 && tmplCand.maxCornerAngleDeg <= T(150.0)) ||
                (matCand.isValid  && matCand.numPatches <= 6  && matCand.maxCornerAngleDeg <= T(150.0));
            if (!haveCompact) {
                mergedCand = computeFieldGuidedMedialSplitting();
                cfCand = computeCrossFieldSplitting();
            }
        }

        const int maxPatchesCap = std::max(1, static_cast<int>(m_options.getInt("MaxPatches")));
        std::vector<SplittingCandidate> validPool;
        if (mergedCand.isValid && mergedCand.numPatches > 0 && mergedCand.numPatches <= maxPatchesCap) validPool.push_back(mergedCand);
        if (cfCand.isValid && cfCand.numPatches > 0 && cfCand.numPatches <= maxPatchesCap) validPool.push_back(cfCand);
        if (matCand.isValid && matCand.numPatches > 0 && matCand.numPatches <= maxPatchesCap) validPool.push_back(matCand);
        if (tmplCand.isValid && tmplCand.numPatches > 0 && tmplCand.numPatches <= maxPatchesCap) validPool.push_back(tmplCand);

        if (validPool.empty()) {
            m_bestCandidate = mergedCand.numPatches > 0 ? mergedCand
                            : (matCand.numPatches > 0 ? matCand : (cfCand.numPatches > 0 ? cfCand : tmplCand));
            return false;
        }

        const bool preferSingleQuad = (m_domain.numHoles() == 0 && m_domain.outer().size() == 4);
        auto computeQ = [this](const SplittingCandidate& c) -> T {
            std::vector<std::vector<patchCorner>> evs;
            c.multiPatch.topology().getEVs(evs, false);
            int nEV = static_cast<int>(evs.size());
            T angPen = (c.maxCornerAngleDeg > T(145.0)) ? T(0.1) * (c.maxCornerAngleDeg - T(145.0)) : T(0);
            return m_weightPatches * c.numPatches
                 + m_weightDistortion * c.distortionEnergy
                 + T(0.25) * nEV
                 + angPen;
        };

        std::sort(validPool.begin(), validPool.end(),
                  [&](const SplittingCandidate& a, const SplittingCandidate& b) {
                      if (preferSingleQuad && (a.numPatches == 1) != (b.numPatches == 1))
                          return a.numPatches == 1;
                      return computeQ(a) < computeQ(b);
                  });

        m_bestCandidate = validPool.front();
    }

    if (m_bestCandidate.multiPatch.nPatches() > 0) {
        ensureExactBoundaryConformity(m_bestCandidate.multiPatch);
    }

    return m_bestCandidate.isValid;
}

// ============================================================================
// Optimization & Quality Evaluation (Winslow Barrier Functional)
// ============================================================================

template <class T>
T gsPatchSplitter<T>::evaluateDistortionEnergy(const gsMultiPatch<T>& mp, T* minDetOut,
                                               T* minAngleOut, T* maxAngleOut) const
{
    if (mp.nPatches() == 0) {
        if (minDetOut) *minDetOut = T(0);
        if (minAngleOut) *minAngleOut = T(0);
        if (maxAngleOut) *maxAngleOut = T(180);
        return T(1e30);
    }

    T totalEnergy = 0.0;
    T globalMinAbsDet = std::numeric_limits<T>::max();
    T globalMinAng = T(360);
    T globalMaxAng = T(0);
    bool allSingleSign = true;

    for (size_t p = 0; p < mp.nPatches(); ++p) {
        T pMinDet = 0, pMaxDet = 0, pMinAbs = 0, pEnergy = 0;
        internal_splitter::evalPatchJacobianAndWinslow(mp.patch(p), 17, pMinDet, pMaxDet, pMinAbs, pEnergy);

        bool pSingleSign = (pMinDet > T(1e-10)) || (pMaxDet < -T(1e-10));
        if (!pSingleSign) {
            allSingleSign = false;
            pEnergy += T(1e6);
        }
        if (pMinAbs < globalMinAbsDet) globalMinAbsDet = pMinAbs;

        T angs[4], pMinA = 0, pMaxA = 180;
        bool degen = false;
        internal_splitter::evalPatchCornerAnglesDeg(mp.patch(p), angs, pMinA, pMaxA, degen);
        globalMinAng = std::min(globalMinAng, pMinA);
        globalMaxAng = std::max(globalMaxAng, pMaxA);

        T angPen = T(0);
        for (int k = 0; k < 4; ++k) {
            T a = angs[k];
            if (a < T(35))  angPen += ((T(35) - a) / T(10)) * ((T(35) - a) / T(10));
            if (a > T(140)) angPen += ((a - T(140)) / T(10)) * ((a - T(140)) / T(10));
        }
        if (degen) angPen += T(1e4);

        totalEnergy += pEnergy + T(5.0) * angPen;
    }

    if (minDetOut) *minDetOut = allSingleSign ? globalMinAbsDet : -globalMinAbsDet;
    if (minAngleOut) *minAngleOut = globalMinAng;
    if (maxAngleOut) *maxAngleOut = globalMaxAng;
    return totalEnergy;
}

template <class T>
bool gsPatchSplitter<T>::optimizeWinslowBarrier(gsMultiPatch<T>& mp, T minDetThreshold) const
{
    if (mp.nPatches() == 0) return false;

    bool allOk = true;
    for (size_t p = 0; p < mp.nPatches(); ++p) {
        auto* tp = dynamic_cast<gsTensorBSpline<2, T>*>(&mp.patch(p));
        if (!tp) continue;

        T minDet = 0, maxDet = 0, minAbs = 0, curEnergy = 0;
        internal_splitter::evalPatchJacobianAndWinslow(*tp, 21, minDet, maxDet, minAbs, curEnergy);

        // If patch is uniformly negative everywhere, flip parameter orientation without altering geometry
        if (maxDet < -minDetThreshold) {
            internal_splitter::safeReverseTensorBSpline(*tp, 0);
            internal_splitter::evalPatchJacobianAndWinslow(*tp, 21, minDet, maxDet, minAbs, curEnergy);
        }

        const int n_u = tp->basis().size(0);
        const int n_v = tp->basis().size(1);
        gsMatrix<T>& coefs = tp->coefs();
        bool valid = (minDet > minDetThreshold) || (maxDet < -minDetThreshold);
        if (n_u <= 2 || n_v <= 2) {
            if (!valid) allOk = false;
            continue;
        }
        // Skip expensive BB optimization on patches that are already injective and near-isometric
        if (valid && minDet > minDetThreshold && curEnergy <= T(2.45)) {
            continue;
        }

        const int Mopt = std::max(21, std::min(96, std::max(n_u, n_v) + 4));
        internal_splitter::TensorGridEval1D<T> tgOpt;
        if (!tgOpt.init(*tp, Mopt)) {
            if (!valid) allOk = false;
            continue;
        }

        const gsMatrix<T> origCoefs = coefs;
        gsMatrix<T> harmCoefs;
        internal_splitter::computeHarmonicInteriorCoefs(*tp, harmCoefs);

        // Build corner-bilinear interior net for emergency recovery
        gsMatrix<T> bilCoefs = origCoefs;
        {
            gsPoint2D<T> P00(origCoefs(0, 0), origCoefs(0, 1));
            gsPoint2D<T> P10(origCoefs(n_u - 1, 0), origCoefs(n_u - 1, 1));
            gsPoint2D<T> P01(origCoefs((n_v - 1) * n_u, 0), origCoefs((n_v - 1) * n_u, 1));
            gsPoint2D<T> P11(origCoefs(n_u - 1 + (n_v - 1) * n_u, 0), origCoefs(n_u - 1 + (n_v - 1) * n_u, 1));
            gsMatrix<T> gr_u = tp->basis().knots(0).greville();
            gsMatrix<T> gr_v = tp->basis().knots(1).greville();
            for (int j = 1; j < n_v - 1; ++j) {
                T eta = gr_v(0, j);
                for (int i = 1; i < n_u - 1; ++i) {
                    T xi = gr_u(0, i);
                    int idx = i + j * n_u;
                    bilCoefs(idx, 0) = (1 - xi) * (1 - eta) * P00.x + xi * (1 - eta) * P10.x +
                                       (1 - xi) * eta * P01.x + xi * eta * P11.x;
                    bilCoefs(idx, 1) = (1 - xi) * (1 - eta) * P00.y + xi * (1 - eta) * P10.y +
                                       (1 - xi) * eta * P01.y + xi * eta * P11.y;
                }
            }
        }

        // Fast candidate selection on Mopt x Mopt: origCoefs vs Greville-Harmonic blends (and bilinear recovery if still folded)
        gsMatrix<T> bestCoefs = origCoefs, Cx, Cy;
        tgOpt.coefsToMatrices(origCoefs, Cx, Cy);
        T bestMin = 0, bestMax = 0, bestAbs = 0, bestEn = 0;
        tgOpt.evalWinslow(Cx, Cy, bestMin, bestMax, bestAbs, bestEn);
        bool bestValid = (bestMin > minDetThreshold);

        auto tryBlend = [&](const gsMatrix<T>& A, const gsMatrix<T>& B, T alpha) {
            gsMatrix<T> cand = (T(1) - alpha) * A + alpha * B;
            tgOpt.coefsToMatrices(cand, Cx, Cy);
            T md = 0, mxd = 0, ma = 0, en = 0;
            tgOpt.evalWinslow(Cx, Cy, md, mxd, ma, en);
            bool cValid = (md > minDetThreshold);
            if ((!bestValid && cValid) ||
                (!bestValid && !cValid && md > bestMin + T(1e-12)) ||
                (bestValid && cValid && en < bestEn)) {
                bestValid = cValid;
                bestMin = md;
                bestMax = mxd;
                bestEn = en;
                bestCoefs = std::move(cand);
            }
        };

        for (T alpha : {T(0.4), T(0.75), T(1.0)})
            tryBlend(origCoefs, harmCoefs, alpha);
        if (!bestValid) {
            for (T beta : {T(0.35), T(0.7), T(1.0)})
                tryBlend(harmCoefs, bilCoefs, beta);
        }

        // Analytic-gradient Barzilai-Borwein + Armijo line-search minimizing regularized Winslow energy
        const T tikW = (std::max(n_u, n_v) > 20) ? T(1e-5) : T(0);
        auto evalWinObjGrad = [&](const gsMatrix<T>& x, T epsVal, T& objOut, T& mdOut, T& mxOut, gsMatrix<T>& gOut) -> void {
            tgOpt.coefsToMatrices(x, Cx, Cy);
            tgOpt.evalRegularizedWinslowAndGrad(Cx, Cy, epsVal, objOut, mdOut, mxOut, gOut);
            if (tikW > T(0)) {
                gsMatrix<T> diff = x - harmCoefs;
                T sumSq = T(0);
                for (int j = 1; j < n_v - 1; ++j) {
                    for (int i = 1; i < n_u - 1; ++i) {
                        int idx = i + j * n_u;
                        T dx = diff(idx, 0), dy = diff(idx, 1);
                        sumSq += dx * dx + dy * dy;
                        gOut(idx, 0) += T(2) * tikW * dx;
                        gOut(idx, 1) += T(2) * tikW * dy;
                    }
                }
                objOut += tikW * sumSq;
            }
        };

        T curMin = bestMin, curMax = bestMax;
        T eps = bestValid
            ? std::max(T(1e-7), T(1e-4) * std::max(T(0), curMax))
            : std::max(T(1.2) * std::max(T(0), -curMin), T(0.03) * std::max(std::abs(curMax), T(1e-6)));

        gsMatrix<T> xCur = bestCoefs, gCur, xPrev, gPrev;
        T objCur = 0;
        evalWinObjGrad(xCur, eps, objCur, curMin, curMax, gCur);

        // Estimate characteristic patch scale for safe initial step normalization
        T bboxDiag = std::hypot(origCoefs.col(0).maxCoeff() - origCoefs.col(0).minCoeff(),
                                origCoefs.col(1).maxCoeff() - origCoefs.col(1).minCoeff());
        const T maxStepDisp = std::max(T(1e-5), T(0.15) * bboxDiag / std::max(T(2), static_cast<T>(std::max(n_u, n_v))));
        const int maxIters = bestValid ? 8 : 24;
        bool hasPrev = false;

        for (int it = 0; it < maxIters; ++it) {
            T gNorm = gCur.norm();
            if (gNorm <= T(1e-12)) break;

            T stepScale = maxStepDisp / std::max(T(1e-14), gCur.cwiseAbs().maxCoeff());
            if (hasPrev) {
                gsMatrix<T> sDiff = xCur - xPrev;
                gsMatrix<T> yDiff = gCur - gPrev;
                T yy = (yDiff.array() * yDiff.array()).sum();
                T sy = (sDiff.array() * yDiff.array()).sum();
                if (yy > T(1e-24) && sy > T(1e-24)) {
                    T bbStep = sy / yy;
                    T bbDisp = bbStep * gCur.cwiseAbs().maxCoeff();
                    if (bbDisp > maxStepDisp)
                        bbStep = stepScale;
                    stepScale = std::max(T(1e-4) * stepScale, bbStep);
                }
            }

            T slope = -(gCur.array() * gCur.array()).sum();
            bool accepted = false;
            gsMatrix<T> xNext, gNext;
            T objNext = 0, mdNext = 0, mxdNext = 0;

            T alpha = stepScale;
            for (int ls = 0; ls < 7; ++ls) {
                xNext = xCur - alpha * gCur;
                evalWinObjGrad(xNext, eps, objNext, mdNext, mxdNext, gNext);
                bool okBarrier = (!bestValid) || (mdNext > minDetThreshold);
                if (okBarrier && objNext <= objCur + T(1e-4) * alpha * slope) {
                    accepted = true;
                    break;
                }
                alpha *= T(0.5);
            }
            if (!accepted) break;

            xPrev = std::move(xCur);
            gPrev = std::move(gCur);
            xCur = std::move(xNext);
            gCur = std::move(gNext);
            curMin = mdNext;
            curMax = mxdNext;
            objCur = objNext;
            hasPrev = true;

            if (curMin > minDetThreshold) {
                bestCoefs = xCur;
                if (!bestValid) {
                    bestValid = true;
                    eps = std::max(T(1e-7), T(1e-4) * std::max(T(0), curMax));
                    evalWinObjGrad(xCur, eps, objCur, curMin, curMax, gCur);
                    hasPrev = false;
                }
            } else if (!bestValid) {
                if (curMin > bestMin + T(1e-12)) {
                    bestMin = curMin;
                    bestCoefs = xCur;
                }
                T newEps = std::max(T(1.1) * std::max(T(0), -curMin), T(0.25) * eps);
                if (newEps < T(0.7) * eps) {
                    eps = newEps;
                    evalWinObjGrad(xCur, eps, objCur, curMin, curMax, gCur);
                    hasPrev = false;
                }
            }
        }

        // Certify on 21x21; keep optimized coefs if valid (or if orig wasn't valid and we improved minDet)
        coefs = bestCoefs;
        internal_splitter::evalPatchJacobianAndWinslow(*tp, 21, minDet, maxDet, minAbs, curEnergy);
        valid = (minDet > minDetThreshold) || (maxDet < -minDetThreshold);
        if (!valid) {
            // Fall back to origCoefs only if orig was actually better on 21x21
            tp->setCoefs(origCoefs);
            T oMin = 0, oMax = 0, oAbs = 0, oEn = 0;
            internal_splitter::evalPatchJacobianAndWinslow(*tp, 21, oMin, oMax, oAbs, oEn);
            bool oValid = (oMin > minDetThreshold) || (oMax < -minDetThreshold);
            if (!oValid && minDet > oMin) {
                tp->setCoefs(bestCoefs);
            } else {
                valid = oValid;
            }
        }

        if (!valid) allOk = false;
    }

    mp.computeTopology(T(1e-5));
    return allOk;
}

template <class T>
bool gsPatchSplitter<T>::optimizeCoupledInterfacesAndInteriors(
    gsMultiPatch<T>& mp,
    T minDetThreshold,
    int maxIters,
    bool slideBoundarySplits) const
{
    using namespace internal_splitter;
    const int nP = static_cast<int>(mp.nPatches());
    if (nP < 2) return optimizeWinslowBarrier(mp, minDetThreshold);

    gsMultiPatch<T> origMp(mp);
    gsASG1LocalFitter<T>::unifyPatchOrientations(mp);
    harmonizeMultiPatchKnots(mp);
    mp.computeTopology(T(1e-4), false, false);

    // Collect TensorBSpline pointers and dimensions
    std::vector<gsTensorBSpline<2, T>*> tbs(nP, nullptr);
    std::vector<int> nu(nP, 0), nv(nP, 0), pOff(nP + 1, 0);
    for (int p = 0; p < nP; ++p) {
        tbs[p] = dynamic_cast<gsTensorBSpline<2, T>*>(&mp.patch(p));
        if (!tbs[p]) { mp = give(origMp); return false; }
        nu[p] = tbs[p]->basis().size(0);
        nv[p] = tbs[p]->basis().size(1);
        pOff[p + 1] = pOff[p] + nu[p] * nv[p];
    }
    const int totalLocal = pOff[nP];
    if (totalLocal == 0) { mp = give(origMp); return false; }

    // Helper: local CP index along side s (west=0, east=1, south=2, north=3) at slot k
    auto sideCount = [&](int p, int s) -> int { return (s < 2) ? nv[p] : nu[p]; };
    auto sideCP = [&](int p, int s, int k) -> int {
        if (s == 0) return pOff[p] + 0 + k * nu[p];
        if (s == 1) return pOff[p] + (nu[p] - 1) + k * nu[p];
        if (s == 2) return pOff[p] + k + 0 * nu[p];
        return pOff[p] + k + (nv[p] - 1) * nu[p];
    };
    auto cpPos = [&](int localIdx) -> gsPoint2D<T> {
        int p = static_cast<int>(std::upper_bound(pOff.begin(), pOff.end(), localIdx) - pOff.begin()) - 1;
        int loc = localIdx - pOff[p];
        return gsPoint2D<T>(tbs[p]->coefs()(loc, 0), tbs[p]->coefs()(loc, 1));
    };

    // Union-Find linking shared interface CPs across patches for exact C0 conformity
    std::vector<int> uf(totalLocal);
    for (int i = 0; i < totalLocal; ++i) uf[i] = i;
    std::function<int(int)> ufFind = [&](int x) -> int {
        return (uf[x] == x) ? x : (uf[x] = ufFind(uf[x]));
    };
    auto ufUnion = [&](int a, int b) {
        a = ufFind(a); b = ufFind(b);
        if (a != b) uf[b] = a;
    };

    std::vector<std::vector<int>> ifcChains; // ordered global local-rep indices along each interface
    for (auto it = mp.iBegin(); it != mp.iEnd(); ++it) {
        const auto& ps1 = it->first();
        const auto& ps2 = it->second();
        int p1 = static_cast<int>(ps1.patch), s1 = ps1.side().m_index - 1;
        int p2 = static_cast<int>(ps2.patch), s2 = ps2.side().m_index - 1;
        int n1 = sideCount(p1, s1), n2 = sideCount(p2, s2);
        if (n1 != n2 || n1 < 2) continue;
        gsPoint2D<T> a0 = cpPos(sideCP(p1, s1, 0)), a1 = cpPos(sideCP(p1, s1, n1 - 1));
        gsPoint2D<T> b0 = cpPos(sideCP(p2, s2, 0)), b1 = cpPos(sideCP(p2, s2, n2 - 1));
        bool flipped = (a0.dist(b1) + a1.dist(b0)) < (a0.dist(b0) + a1.dist(b1));
        if ((flipped ? (a0.dist(b1) + a1.dist(b0)) : (a0.dist(b0) + a1.dist(b1))) > T(1e-3))
            continue;
        std::vector<int> chain(n1);
        for (int k = 0; k < n1; ++k) {
            int idx1 = sideCP(p1, s1, k);
            int idx2 = sideCP(p2, s2, flipped ? (n1 - 1 - k) : k);
            ufUnion(idx1, idx2);
            chain[k] = idx1;
        }
        ifcChains.push_back(std::move(chain));
    }

    // Compact global DOF ids and mark all domain-boundary CPs as fixed
    std::vector<int> rootToDof(totalLocal, -1);
    int nDof = 0;
    for (int i = 0; i < totalLocal; ++i) {
        int r = ufFind(i);
        if (rootToDof[r] < 0) rootToDof[r] = nDof++;
    }
    std::vector<int> localToDof(totalLocal);
    std::vector<bool> dofFixed(nDof, false);
    gsMatrix<T> xDof(nDof, 2);
    xDof.setZero();
    std::vector<int> dofCount(nDof, 0);
    for (int p = 0; p < nP; ++p) {
        const auto& cf = tbs[p]->coefs();
        for (int loc = 0; loc < nu[p] * nv[p]; ++loc) {
            int d = rootToDof[ufFind(pOff[p] + loc)];
            localToDof[pOff[p] + loc] = d;
            xDof(d, 0) += cf(loc, 0);
            xDof(d, 1) += cf(loc, 1);
            ++dofCount[d];
        }
        for (int s = 0; s < 4; ++s) {
            boxSide bs = static_cast<boxSide>(s + 1);
            if (!mp.isInterface(patchSide(p, bs))) {
                int ns = sideCount(p, s);
                for (int k = 0; k < ns; ++k)
                    dofFixed[localToDof[sideCP(p, s, k)]] = true;
            }
        }
    }
    for (int d = 0; d < nDof; ++d) {
        if (dofCount[d] > 0) { xDof(d, 0) /= dofCount[d]; xDof(d, 1) /= dofCount[d]; }
    }

    // Convert interface chains into global DOF id sequences for smoothness regularization
    std::vector<std::vector<int>> ifcDofChains;
    ifcDofChains.reserve(ifcChains.size());
    for (const auto& ch : ifcChains) {
        std::vector<int> dch(ch.size());
        for (size_t k = 0; k < ch.size(); ++k) dch[k] = localToDof[ch[k]];
        ifcDofChains.push_back(std::move(dch));
    }

    // Precompute 1D evaluation matrices and initial patch area targets (prevents scale-collapse!)
    std::vector<TensorGridEval1D<T>> tgs(nP);
    std::vector<T> patchAreaRef(nP, T(1));
    T sumPosArea = T(0);
    int nPosArea = 0;
    for (int p = 0; p < nP; ++p) {
        int Mopt = std::max(17, std::min(33, std::max(nu[p], nv[p]) + 4));
        if (!tgs[p].init(*tbs[p], Mopt)) { mp = give(origMp); return false; }
        T pMn = 0, pMx = 0, pMa = 0, pEn = 0;
        evalPatchJacobianAndWinslow(*tbs[p], 13, pMn, pMx, pMa, pEn);
        patchAreaRef[p] = std::max(pMa, T(0));
        if (pMa > T(1e-12)) { sumPosArea += pMa; ++nPosArea; }
    }
    const T meanPatchArea = (nPosArea > 0) ? (sumPosArea / static_cast<T>(nPosArea)) : T(1);
    const T minTargetArea = std::max(T(1e-8), T(0.20) * meanPatchArea);
    for (int p = 0; p < nP; ++p)
        patchAreaRef[p] = std::max(patchAreaRef[p], minTargetArea);

    // Precompute initial chord lengths squared per interface chain for scale-invariant fairness
    std::vector<T> ifcChordSq(ifcDofChains.size(), T(1));
    for (size_t ci = 0; ci < ifcDofChains.size(); ++ci) {
        const auto& ch = ifcDofChains[ci];
        if (ch.size() >= 2) {
            T dx = xDof(ch.back(), 0) - xDof(ch.front(), 0);
            T dy = xDof(ch.back(), 1) - xDof(ch.front(), 1);
            ifcChordSq[ci] = std::max(dx * dx + dy * dy, T(1e-12));
        }
    }

    auto scatterToPatches = [&](const gsMatrix<T>& x) {
        for (int p = 0; p < nP; ++p) {
            gsMatrix<T>& cf = tbs[p]->coefs();
            for (int loc = 0; loc < nu[p] * nv[p]; ++loc) {
                int d = localToDof[pOff[p] + loc];
                cf(loc, 0) = x(d, 0);
                cf(loc, 1) = x(d, 1);
            }
        }
    };

    const T lamArea  = T(0.35);
    const T lamBend  = T(0.05);
    const T lamOrtho = T(0.22);
    auto evalGlobalObjGrad = [&](const gsMatrix<T>& x, T epsVal, T& objOut, T& gMinDet, T& gMaxDet, gsMatrix<T>& gOut) {
        objOut = T(0);
        gMinDet = std::numeric_limits<T>::max();
        gMaxDet = -std::numeric_limits<T>::max();
        gOut.setZero(nDof, 2);

        gsMatrix<T> Cx, Cy, pGrad, pCoefs;
        for (int p = 0; p < nP; ++p) {
            pCoefs.resize(nu[p] * nv[p], 2);
            for (int loc = 0; loc < nu[p] * nv[p]; ++loc) {
                int d = localToDof[pOff[p] + loc];
                pCoefs(loc, 0) = x(d, 0);
                pCoefs(loc, 1) = x(d, 1);
            }
            tgs[p].coefsToMatrices(pCoefs, Cx, Cy);
            T pObj = 0, pMn = 0, pMx = 0;
            tgs[p].evalRegularizedWinslowAndGradAll(Cx, Cy, epsVal, pObj, pMn, pMx, pGrad, patchAreaRef[p], lamArea, lamOrtho);
            objOut += pObj;
            if (pMn < gMinDet) gMinDet = pMn;
            if (pMx > gMaxDet) gMaxDet = pMx;
            for (int loc = 0; loc < nu[p] * nv[p]; ++loc) {
                int d = localToDof[pOff[p] + loc];
                if (!dofFixed[d]) {
                    gOut(d, 0) += pGrad(loc, 0);
                    gOut(d, 1) += pGrad(loc, 1);
                }
            }
        }

        // Scale-invariant discrete second-difference fairness regularizer along every curved interior interface
        for (size_t ci = 0; ci < ifcDofChains.size(); ++ci) {
            const auto& ch = ifcDofChains[ci];
            const int nc = static_cast<int>(ch.size());
            if (nc < 3) continue;
            T wB = lamBend / (ifcChordSq[ci] * static_cast<T>(nc - 2));
            for (int k = 1; k + 1 < nc; ++k) {
                int dm = ch[k - 1], d0 = ch[k], dp = ch[k + 1];
                T rx = x(dm, 0) - T(2) * x(d0, 0) + x(dp, 0);
                T ry = x(dm, 1) - T(2) * x(d0, 1) + x(dp, 1);
                objOut += wB * (rx * rx + ry * ry);
                if (!dofFixed[dm]) { gOut(dm, 0) += T(2) * wB * rx; gOut(dm, 1) += T(2) * wB * ry; }
                if (!dofFixed[d0]) { gOut(d0, 0) -= T(4) * wB * rx; gOut(d0, 1) -= T(4) * wB * ry; }
                if (!dofFixed[dp]) { gOut(dp, 0) += T(2) * wB * rx; gOut(dp, 1) += T(2) * wB * ry; }
            }
        }
    };

    T curObj = 0, curMin = 0, curMax = 0;
    gsMatrix<T> gCur, xPrev, gPrev;
    evalGlobalObjGrad(xDof, T(1e-6), curObj, curMin, curMax, gCur);

    bool initValid = (curMin > minDetThreshold);
    T eps = initValid
        ? std::max(T(1e-7), T(1e-4) * std::max(T(0), curMax))
        : std::max(T(1.2) * std::max(T(0), -curMin), T(0.03) * std::max(std::abs(curMax), T(1e-6)));
    evalGlobalObjGrad(xDof, eps, curObj, curMin, curMax, gCur);

    gsMatrix<T> xCur = xDof, bestX = xDof;
    bool bestValid = initValid;
    T bestMin = curMin, bestObj = curObj;

    T bboxDiag = std::hypot(xDof.col(0).maxCoeff() - xDof.col(0).minCoeff(),
                            xDof.col(1).maxCoeff() - xDof.col(1).minCoeff());
    const T maxDisp = std::max(T(1e-5), T(0.08) * bboxDiag / std::max(T(3), std::sqrt(static_cast<T>(nP))));
    bool hasPrev = false;

    for (int it = 0; it < maxIters; ++it) {
        T gMax = gCur.cwiseAbs().maxCoeff();
        if (gMax <= T(1e-12)) break;
        T stepScale = maxDisp / gMax;
        if (hasPrev) {
            gsMatrix<T> sDiff = xCur - xPrev;
            gsMatrix<T> yDiff = gCur - gPrev;
            T yy = (yDiff.array() * yDiff.array()).sum();
            T sy = (sDiff.array() * yDiff.array()).sum();
            if (yy > T(1e-24) && sy > T(1e-24)) {
                T bbStep = sy / yy;
                if (bbStep * gMax <= maxDisp)
                    stepScale = std::max(T(1e-4) * stepScale, bbStep);
            }
        }

        T slope = -(gCur.array() * gCur.array()).sum();
        bool accepted = false;
        gsMatrix<T> xNext, gNext;
        T objNext = 0, mdNext = 0, mxNext = 0;
        T alpha = stepScale;
        for (int ls = 0; ls < 7; ++ls) {
            xNext = xCur - alpha * gCur;
            evalGlobalObjGrad(xNext, eps, objNext, mdNext, mxNext, gNext);
            bool okBarrier = (!bestValid) || (mdNext > minDetThreshold);
            if (okBarrier && objNext <= curObj + T(1e-4) * alpha * slope) {
                accepted = true;
                break;
            }
            alpha *= T(0.5);
        }
        if (!accepted) break;

        xPrev = std::move(xCur); gPrev = std::move(gCur);
        xCur = std::move(xNext); gCur = std::move(gNext);
        curObj = objNext; curMin = mdNext; curMax = mxNext;
        hasPrev = true;

        if (curMin > minDetThreshold) {
            if (!bestValid || curObj < bestObj) {
                bestX = xCur;
                bestObj = curObj;
                bestMin = curMin;
            }
            if (!bestValid) {
                bestValid = true;
                eps = std::max(T(1e-7), T(1e-4) * std::max(T(0), curMax));
                evalGlobalObjGrad(xCur, eps, curObj, curMin, curMax, gCur);
                hasPrev = false;
            }
        } else if (!bestValid) {
            if (curMin > bestMin + T(1e-12)) {
                bestMin = curMin;
                bestX = xCur;
            }
            T newEps = std::max(T(1.1) * std::max(T(0), -curMin), T(0.25) * eps);
            if (newEps < T(0.7) * eps) {
                eps = newEps;
                evalGlobalObjGrad(xCur, eps, curObj, curMin, curMax, gCur);
                hasPrev = false;
            }
        }
    }

    scatterToPatches(bestX);

    // Phase 1B: Exact CAD Boundary Split-Parameter Sliding (1D u_b along smooth CAD curves)
    // Slides smooth 2-patch boundary split points along the exact joined CAD B-spline via splitAt,
    // prioritized worst-corner-skew-first across the entire mesh (unblocked beyond 16 vertices!).
    if (slideBoundarySplits && m_hasDomain && m_domain.numLoops() > 0 && bestValid) {
        const T minAngTol = std::min(static_cast<T>(m_options.getReal("MinCornerAngleDeg")), T(0.85) * getDomainMinCornerAngleDeg());
        const T maxAngTol = std::max(static_cast<T>(m_options.getReal("MaxCornerAngleDeg")), T(164.2));
        auto extractSide = [&](const gsTensorBSpline<2, T>& tb, int sIdx) -> CurveSeg<T> {
            typename gsGeometry<T>::uPtr b = tb.boundary(static_cast<boxSide>(sIdx + 1));
            auto* bsp = dynamic_cast<gsBSpline<T>*>(b.get());
            return bsp ? CurveSeg<T>(*bsp) : CurveSeg<T>();
        };
        for (int sweep = 0; sweep < 2; ++sweep) {
            mp.computeTopology(T(1e-4), false, false);
            bool anySlid = false;
            struct SlideIfc { int p1, s1, p2, s2; T skew; };
            std::vector<SlideIfc> ifcList;
            for (auto it = mp.iBegin(); it != mp.iEnd(); ++it) {
                int p1 = static_cast<int>(it->first().patch), s1 = it->first().side().m_index - 1;
                int p2 = static_cast<int>(it->second().patch), s2 = it->second().side().m_index - 1;
                if (p1 == p2) continue;
                T a1[4], mn1 = 90, mx1 = 90, a2[4], mn2 = 90, mx2 = 90; bool d1 = false, d2 = false;
                evalPatchCornerAnglesDeg(mp.patch(p1), a1, mn1, mx1, d1);
                evalPatchCornerAnglesDeg(mp.patch(p2), a2, mn2, mx2, d2);
                T sk = std::max(std::max(T(90) - mn1, mx1 - T(90)), std::max(T(90) - mn2, mx2 - T(90)));
                ifcList.push_back({p1, s1, p2, s2, sk});
            }
            std::sort(ifcList.begin(), ifcList.end(), [](const SlideIfc& a, const SlideIfc& b) { return a.skew > b.skew; });

            int triedVerts = 0;
            for (const auto& item : ifcList) {
                if (triedVerts >= 32) break;
                if (item.skew < T(10.0)) continue;
                int p1 = item.p1, s1 = item.s1, p2 = item.p2, s2 = item.s2;
                auto* tb1 = dynamic_cast<gsTensorBSpline<2, T>*>(&mp.patch(p1));
                auto* tb2 = dynamic_cast<gsTensorBSpline<2, T>*>(&mp.patch(p2));
                if (!tb1 || !tb2) continue;

                // Transverse sides of p1 and p2 meeting at the endpoints of this shared interface
                int tr1[2] = {(s1 < 2) ? 2 : 0, (s1 < 2) ? 3 : 1};
                int tr2[2] = {(s2 < 2) ? 2 : 0, (s2 < 2) ? 3 : 1};

                for (int eEnd = 0; eEnd < 2 && triedVerts < 32; ++eEnd) {
                    int bS1 = tr1[eEnd];
                    if (mp.isInterface(patchSide(p1, static_cast<boxSide>(bS1 + 1)))) continue;
                    CurveSeg<T> segB1 = extractSide(*tb1, bS1);
                    CurveSeg<T> segIfc = extractSide(*tb1, s1);
                    if (!segB1.bsp || !segIfc.bsp) continue;
                    // Identify the corner vertex V shared by side s1 and boundary side bS1
                    gsPoint2D<T> V = (segB1.p0().dist(segIfc.p0()) < T(1e-4) || segB1.p0().dist(segIfc.p1()) < T(1e-4))
                                   ? segB1.p0() : segB1.p1();
                    // Find matching boundary side bS2 on p2 touching V
                    int bS2 = -1;
                    CurveSeg<T> segB2;
                    for (int c2 = 0; c2 < 2; ++c2) {
                        int candS2 = tr2[c2];
                        if (mp.isInterface(patchSide(p2, static_cast<boxSide>(candS2 + 1)))) continue;
                        CurveSeg<T> cs2 = extractSide(*tb2, candS2);
                        if (cs2.bsp && (cs2.p0().dist(V) < T(1e-4) || cs2.p1().dist(V) < T(1e-4))) {
                            bS2 = candS2; segB2 = cs2; break;
                        }
                    }
                    if (bS2 < 0) continue;

                    // Orient segB1 -> V -> segB2 continuously along the CAD boundary
                    CurveSeg<T> b1ToV = (segB1.p1().dist(V) < segB1.p0().dist(V)) ? segB1 : segB1.reversed();
                    CurveSeg<T> b2FromV = (segB2.p0().dist(V) < segB2.p1().dist(V)) ? segB2 : segB2.reversed();
                    gsPoint2D<T> tIn = evalCurveDeriv(*b1ToV.bsp, T(1));
                    gsPoint2D<T> tOut = evalCurveDeriv(*b2FromV.bsp, T(0));
                    T nIn = std::hypot(tIn.x, tIn.y), nOut = std::hypot(tOut.x, tOut.y);
                    if (nIn <= T(1e-12) || nOut <= T(1e-12)) continue;
                    // Only slide smooth split vertices (cos > 0.96 => collinear within ~16 deg, never sharp corners!)
                    if ((tIn.x * tOut.x + tIn.y * tOut.y) / (nIn * nOut) < T(0.96)) continue;

                    // Ensure ONLY p1 and p2 touch V (clean valence-2 boundary split vertex)
                    int incP = 0;
                    for (int q = 0; q < nP; ++q) {
                        const auto& qcf = mp.patch(q).coefs();
                        int qnu = dynamic_cast<const gsTensorBSpline<2, T>*>(&mp.patch(q))->basis().size(0);
                        int qnv = dynamic_cast<const gsTensorBSpline<2, T>*>(&mp.patch(q))->basis().size(1);
                        int cIdx[4] = {0, qnu - 1, qnu * qnv - 1, (qnv - 1) * qnu};
                        for (int ck : cIdx) {
                            if (std::hypot(qcf(ck, 0) - V.x, qcf(ck, 1) - V.y) < T(1e-4)) { ++incP; break; }
                        }
                    }
                    if (incP != 2) continue;

                    auto joinedPtr = joinCurveSegsExact(std::vector<CurveSeg<T>>{b1ToV, b2FromV}, m_degree);
                    if (!joinedPtr) continue;
                    CurveSeg<T> joinedArc(*joinedPtr);
                    T u0 = joinedArc.closestParam(V);
                    if (u0 < T(0.15) || u0 > T(0.85)) continue;
                    ++triedVerts;

                    // Measure current combined Winslow energy on p1 + p2
                    T mn1 = 0, mx1 = 0, ma1 = 0, en1 = 0, mn2 = 0, mx2 = 0, ma2 = 0, en2 = 0;
                    evalPatchJacobianAndWinslow(*tb1, 21, mn1, mx1, ma1, en1);
                    evalPatchJacobianAndWinslow(*tb2, 21, mn2, mx2, ma2, en2);
                    T baseEn = en1 + en2;

                    // Extract the remaining 3 untouched sides of p1 and p2
                    auto getOrientPatch = [&](const gsTensorBSpline<2, T>& tb, int ifcSide, int bndSide,
                                              CurveSeg<T>& outOppIfc, CurveSeg<T>& outOppBnd) {
                        int oppIfc = (ifcSide ^ 1);
                        int oppBnd = (bndSide ^ 1);
                        outOppIfc = extractSide(tb, oppIfc);
                        outOppBnd = extractSide(tb, oppBnd);
                    };
                    CurveSeg<T> p1OppIfc, p1OppBnd, p2OppIfc, p2OppBnd;
                    getOrientPatch(*tb1, s1, bS1, p1OppIfc, p1OppBnd);
                    getOrientPatch(*tb2, s2, bS2, p2OppIfc, p2OppBnd);
                    // Other endpoint W of the shared interface (opposite V)
                    gsPoint2D<T> W = (segIfc.p0().dist(V) < segIfc.p1().dist(V)) ? segIfc.p1() : segIfc.p0();
                    CurveSeg<T> ifcWtoV = (segIfc.p0().dist(W) < segIfc.p1().dist(W)) ? segIfc : segIfc.reversed();

                    T bestSlideEn = baseEn;
                    std::unique_ptr<gsTensorBSpline<2, T>> bestP1, bestP2;
                    for (T du : {T(-0.16), T(-0.08), T(0.08), T(0.16)}) {
                        T uTry = std::max(T(0.12), std::min(T(0.88), u0 + du));
                        if (std::abs(uTry - u0) < T(0.02)) continue;
                        auto spArc = joinedArc.splitAtParam(uTry);
                        gsPoint2D<T> Vnew = spArc.first.p1();

                        // Smoothly shift shared interface endpoint V -> Vnew using linear ramp of CPs (preserves curvature!)
                        gsBSpline<T> newIfcBsp(*ifcWtoV.bsp);
                        gsMatrix<T>& icf = newIfcBsp.coefs();
                        const int nic = icf.rows();
                        T dx = Vnew.x - V.x, dy = Vnew.y - V.y;
                        for (int k = 1; k < nic; ++k) {
                            T ramp = static_cast<T>(k) / static_cast<T>(nic - 1);
                            icf(k, 0) += ramp * dx;
                            icf(k, 1) += ramp * dy;
                        }
                        CurveSeg<T> newIfc(newIfcBsp); // W -> Vnew

                        // Rebuild p1 and p2 as conforming Coons patches from their 4 exact sides
                        auto buildRebuiltPatch = [&](const CurveSeg<T>& bSide, const CurveSeg<T>& oppIfc,
                                                     const CurveSeg<T>& oppBnd, const CurveSeg<T>& ifcSeg) -> std::unique_ptr<gsTensorBSpline<2, T>> {
                            // Chain 4 segments around the patch loop starting from bSide
                            CurveSeg<T> c0 = bSide;
                            gsPoint2D<T> cur = c0.p1();
                            auto pickNext = [&](const CurveSeg<T>& a, const CurveSeg<T>& b, const CurveSeg<T>& c,
                                                bool& uA, bool& uB, bool& uC, CurveSeg<T>& out) -> bool {
                                const CurveSeg<T>* arr[3] = {&a, &b, &c};
                                bool* used[3] = {&uA, &uB, &uC};
                                for (int k = 0; k < 3; ++k) {
                                    if (*used[k]) continue;
                                    if (arr[k]->p0().dist(cur) < T(1e-3)) { *used[k] = true; out = *arr[k]; cur = out.p1(); return true; }
                                    if (arr[k]->p1().dist(cur) < T(1e-3)) { *used[k] = true; out = arr[k]->reversed(); cur = out.p1(); return true; }
                                }
                                return false;
                            };
                            bool uA = false, uB = false, uC = false;
                            CurveSeg<T> c1, c2, c3;
                            if (!pickNext(ifcSeg, oppBnd, oppIfc, uA, uB, uC, c1) ||
                                !pickNext(ifcSeg, oppBnd, oppIfc, uA, uB, uC, c2) ||
                                !pickNext(ifcSeg, oppBnd, oppIfc, uA, uB, uC, c3) ||
                                cur.dist(c0.p0()) > T(1e-3))
                                return nullptr;
                            // Ensure CCW orientation
                            gsPoint2D<T> cc[4] = {c0.p0(), c1.p0(), c2.p0(), c3.p0()};
                            T ar2 = cross2(cc[0], cc[1], cc[2]) + cross2(cc[0], cc[2], cc[3]);
                            gsTensorBSpline<2, T>* raw = (ar2 >= T(0))
                                ? makeCoonsFromCCWCycle(c0, c1, c2, c3, m_degree)
                                : makeCoonsFromCCWCycle(c3.reversed(), c2.reversed(), c1.reversed(), c0.reversed(), m_degree);
                            if (!raw) return nullptr;
                            std::unique_ptr<gsTensorBSpline<2, T>> ptr(raw);
                            T mdSlide = evalPatchMinDetDense(*ptr, 49);
                            if (mdSlide <= minDetThreshold) {
                                quickHarmonicPatchRepair(*ptr, std::max(minDetThreshold, T(1e-7)));
                                mdSlide = evalPatchMinDetDense(*ptr, 49);
                            }
                            if (mdSlide <= minDetThreshold && mdSlide > -T(0.05)) {
                                gsMultiPatch<T> tmpOne;
                                tmpOne.addPatch(ptr->clone());
                                if (optimizeWinslowBarrier(tmpOne, std::max(minDetThreshold, T(1e-7)))) {
                                    if (auto* fx = dynamic_cast<gsTensorBSpline<2, T>*>(&tmpOne.patch(0)))
                                        *ptr = *fx;
                                }
                            }
                            return ptr;
                        };

                        auto np1 = buildRebuiltPatch(spArc.first, p1OppIfc, p1OppBnd, newIfc);
                        auto np2 = buildRebuiltPatch(spArc.second, p2OppIfc, p2OppBnd, newIfc);
                        if (!np1 || !np2) continue;

                        T d49_1 = evalPatchMinDetDense(*np1, 49);
                        T d49_2 = evalPatchMinDetDense(*np2, 49);
                        if (d49_1 <= minDetThreshold || d49_2 <= minDetThreshold) continue;

                        T angs[4], mnA1 = 0, mxA1 = 180, mnA2 = 0, mxA2 = 180; bool deg1 = false, deg2 = false;
                        evalPatchCornerAnglesDeg(*np1, angs, mnA1, mxA1, deg1);
                        evalPatchCornerAnglesDeg(*np2, angs, mnA2, mxA2, deg2);
                        if (deg1 || deg2 || mnA1 < minAngTol || mxA1 > maxAngTol || mnA2 < minAngTol || mxA2 > maxAngTol)
                            continue;

                        T tmn1, tmx1, tma1, ten1, tmn2, tmx2, tma2, ten2;
                        evalPatchJacobianAndWinslow(*np1, 21, tmn1, tmx1, tma1, ten1);
                        evalPatchJacobianAndWinslow(*np2, 21, tmn2, tmx2, tma2, ten2);
                        if (ten1 + ten2 < bestSlideEn * T(0.97)) {
                            bestSlideEn = ten1 + ten2;
                            bestP1 = std::move(np1);
                            bestP2 = std::move(np2);
                        }
                    }
                    if (bestP1 && bestP2) {
                        *tb1 = *bestP1;
                        *tb2 = *bestP2;
                        anySlid = true;
                        break;
                    }
                }
            }
            if (!anySlid) break;
        }
        harmonizeMultiPatchKnots(mp);
    }

    gsASG1LocalFitter<T>::unifyPatchOrientations(mp);
    mp.computeTopology(T(1e-4), false, false);

    // Fast exit if rescue optimization failed to achieve positive Jacobian on optimization grid
    if (!bestValid && !initValid) {
        mp = give(origMp);
        return false;
    }

    // Final safety rollback: never degrade a previously dense-injective multipatch or crush minDet / angles
    auto denseMinOf = [](const gsMultiPatch<T>& m) -> T {
        T gMin = std::numeric_limits<T>::max();
        for (size_t p = 0; p < m.nPatches(); ++p)
            gMin = std::min(gMin, evalPatchMinDetDense(m.patch(p), 49));
        return gMin;
    };
    auto angleBoundsOf = [](const gsMultiPatch<T>& m, T& minA, T& maxA) {
        minA = T(180); maxA = T(0);
        for (size_t p = 0; p < m.nPatches(); ++p) {
            T angs[4], mn = 0, mx = 180; bool deg = false;
            evalPatchCornerAnglesDeg(m.patch(p), angs, mn, mx, deg);
            if (deg) { minA = T(0); maxA = T(180); return; }
            minA = std::min(minA, mn);
            maxA = std::max(maxA, mx);
        }
    };
    T dOrig = denseMinOf(origMp);
    T dNew = denseMinOf(mp);
    bool topoValid = (2 * mp.nInterfaces() + mp.nBoundary() == 4 * mp.nPatches());
    bool reject = !topoValid ||
                  (dOrig > minDetThreshold && dNew < T(0.92) * dOrig) ||
                  (dOrig <= minDetThreshold && dNew < dOrig);
    if (!reject && dOrig > minDetThreshold) {
        T oMinA = 180, oMaxA = 0, nMinA = 180, nMaxA = 0;
        angleBoundsOf(origMp, oMinA, oMaxA);
        angleBoundsOf(mp, nMinA, nMaxA);
        if ((oMinA >= T(15) && nMinA < T(0.90) * oMinA) ||
            (oMaxA <= T(165) && nMaxA > std::max(T(165), oMaxA + T(2)))) {
            reject = true;
        }
    }
    if (reject) {
        mp = give(origMp);
        return (dOrig > minDetThreshold);
    }
    if (m_hasDomain && !verifyBoundaryMatch(m_domain, mp, T(1e-4))) {
        mp = give(origMp);
        return (dOrig > minDetThreshold);
    }
    return (dNew > minDetThreshold);
}

template <class T>
void gsPatchSplitter<T>::ensureExactBoundaryConformity(gsMultiPatch<T>& mp)
{
    if (mp.nPatches() > 0) {
        mp.computeTopology(T(1e-5));
    }
}

// ============================================================================
// Diagnostics & Validation
// ============================================================================

template <class T>
bool gsPatchSplitter<T>::isConforming() const
{
    const auto& mp = m_bestCandidate.multiPatch;
    if (mp.nPatches() == 0) return false;
    for (size_t p = 0; p < mp.nPatches(); ++p) {
        if (mp.patch(p).domainDim() != 2) return false;
    }
    gsMatrix<T> uMat(1, 1), pt1, pt2;
    for (const auto& ifc : mp.interfaces()) {
        auto b1 = mp.patch(ifc.first().patch).boundary(ifc.first().side());
        auto b2 = mp.patch(ifc.second().patch).boundary(ifc.second().side());
        gsMatrix<T> pr1 = b1->parameterRange(), pr2 = b2->parameterRange();
        uMat(0, 0) = T(0.5) * (pr1(0, 0) + pr1(0, 1));
        b1->eval_into(uMat, pt1);
        uMat(0, 0) = T(0.5) * (pr2(0, 0) + pr2(0, 1));
        b2->eval_into(uMat, pt2);
        if (std::hypot(pt1(0, 0) - pt2(0, 0), pt1(1, 0) - pt2(1, 0)) > T(1e-4))
            return false;
    }
    return true;
}

template <class T>
bool gsPatchSplitter<T>::checkJacobians(T* min_det) const
{
    const auto& mp = m_bestCandidate.multiPatch;
    if (mp.nPatches() == 0) {
        if (min_det) *min_det = T(0);
        return false;
    }
    T mAbs = std::numeric_limits<T>::max();
    bool ok = true;
    for (size_t p = 0; p < mp.nPatches(); ++p) {
        T pAbs = 0;
        if (!internal_splitter::isPatchDetSingleSign(mp.patch(p), 21, T(1e-7), &pAbs))
            ok = false;
        if (pAbs < mAbs) mAbs = pAbs;
    }
    if (min_det) *min_det = ok ? mAbs : -mAbs;
    return ok;
}

template <class T>
bool gsPatchSplitter<T>::checkPatchCornersConvex(T* min_tri_area2) const
{
    const auto& mp = m_bestCandidate.multiPatch;
    if (mp.nPatches() == 0) return false;
    T worst = std::numeric_limits<T>::max();
    bool ok = true;
    gsMatrix<T> u(2, 1), val;
    const T uv[4][2] = {{0, 0}, {1, 0}, {1, 1}, {0, 1}};
    for (size_t p = 0; p < mp.nPatches(); ++p) {
        gsPoint2D<T> c[4];
        for (int k = 0; k < 4; ++k) {
            u(0, 0) = uv[k][0]; u(1, 0) = uv[k][1];
            mp.patch(p).eval_into(u, val);
            c[k] = gsPoint2D<T>(val(0, 0), val(1, 0));
        }
        T triA = 0;
        if (!internal_splitter::isStrictlyConvexQuad(c, &triA)) ok = false;
        if (triA < worst) worst = triA;
    }
    if (min_tri_area2) *min_tri_area2 = worst;
    return ok;
}

template <class T>
bool gsPatchSplitter<T>::checkNoSelfIntersections() const
{
    const auto& mp = m_bestCandidate.multiPatch;
    size_t nP = mp.nPatches();
    if (nP == 0) return false;

    struct PatchAABB { T xMin, xMax, yMin, yMax; };
    std::vector<PatchAABB> boxes(nP);
    std::vector<std::vector<gsPoint2D<T>>> polys(nP);

    int globalNSide = 24;
    for (size_t p = 0; p < nP; ++p) {
        if (const auto* tb = dynamic_cast<const gsTensorBSpline<2, T>*>(&mp.patch(p)))
            globalNSide = std::max(globalNSide, std::min(64, 2 * std::max<int>(tb->basis().component(0).size(), tb->basis().component(1).size())));
    }

    for (size_t p = 0; p < nP; ++p) {
        polys[p] = internal_splitter::samplePatchBoundaryPoly(mp.patch(p), globalNSide);
        T xMin = std::numeric_limits<T>::max(), xMax = -std::numeric_limits<T>::max();
        T yMin = std::numeric_limits<T>::max(), yMax = -std::numeric_limits<T>::max();
        for (const auto& pt : polys[p]) {
            xMin = std::min(xMin, pt.x); xMax = std::max(xMax, pt.x);
            yMin = std::min(yMin, pt.y); yMax = std::max(yMax, pt.y);
        }
        boxes[p] = {xMin - T(1e-7), xMax + T(1e-7), yMin - T(1e-7), yMax + T(1e-7)};

        int m = static_cast<int>(polys[p].size());
        for (int i = 0; i < m; ++i) {
            for (int j = i + 2; j < m; ++j) {
                if (i == 0 && j == m - 1) continue;
                if (internal_splitter::segProperIntersect(polys[p][i], polys[p][(i + 1) % m],
                                                          polys[p][j], polys[p][(j + 1) % m]))
                    return false;
            }
        }
    }

    gsMatrix<T> u(2, 9), val;
    const T pts[3] = {T(0.25), T(0.5), T(0.75)};
    int col = 0;
    for (int a = 0; a < 3; ++a) {
        for (int b = 0; b < 3; ++b) {
            u(0, col) = pts[a]; u(1, col) = pts[b];
            ++col;
        }
    }
    std::vector<std::vector<gsPoint2D<T>>> innerPts(nP, std::vector<gsPoint2D<T>>(9));
    for (size_t p = 0; p < nP; ++p) {
        mp.patch(p).eval_into(u, val);
        for (int k = 0; k < 9; ++k) {
            innerPts[p][k] = gsPoint2D<T>(val(0, k), val(1, k));
        }
    }

    for (size_t p = 0; p < nP; ++p) {
        for (size_t q = p + 1; q < nP; ++q) {
            if (boxes[p].xMax < boxes[q].xMin || boxes[p].xMin > boxes[q].xMax ||
                boxes[p].yMax < boxes[q].yMin || boxes[p].yMin > boxes[q].yMax)
                continue;

            int mp1 = static_cast<int>(polys[p].size());
            int mq1 = static_cast<int>(polys[q].size());
            bool sharedSide[4][4] = {{false}};
            for (int sp = 0; sp < 4; ++sp) {
                const auto& p0 = polys[p][sp * globalNSide];
                const auto& pm = polys[p][sp * globalNSide + globalNSide / 2];
                const auto& p1 = polys[p][((sp + 1) % 4) * globalNSide];
                for (int sq = 0; sq < 4; ++sq) {
                    const auto& q0 = polys[q][sq * globalNSide];
                    const auto& qm = polys[q][sq * globalNSide + globalNSide / 2];
                    const auto& q1 = polys[q][((sq + 1) % 4) * globalNSide];
                    if (((p0.dist(q0) < T(1e-4) && p1.dist(q1) < T(1e-4)) ||
                         (p0.dist(q1) < T(1e-4) && p1.dist(q0) < T(1e-4))) &&
                        pm.dist(qm) < T(1e-4)) {
                        sharedSide[sp][sq] = true;
                    }
                }
            }
            for (int i = 0; i < mp1; ++i) {
                const int sp = i / globalNSide;
                const auto& A = polys[p][i];
                const auto& B = polys[p][(i + 1) % mp1];
                for (int j = 0; j < mq1; ++j) {
                    if (sharedSide[sp][j / globalNSide]) continue;
                    const auto& C = polys[q][j];
                    const auto& D = polys[q][(j + 1) % mq1];
                    if (A.dist(C) < T(1e-5) || A.dist(D) < T(1e-5) ||
                        B.dist(C) < T(1e-5) || B.dist(D) < T(1e-5))
                        continue;
                    if (internal_splitter::segProperIntersect(A, B, C, D))
                        return false;
                }
            }
            for (int k = 0; k < 9; ++k) {
                if (internal_splitter::pointInPolyStrict(innerPts[p][k], polys[q]))
                    return false;
            }
        }
    }
    return true;
}

template <class T>
bool gsPatchSplitter<T>::validateMultiPatch(
    const gsPlanarDomain<T>& dom,
    const gsMultiPatch<T>& mp,
    T bdrTol,
    std::string* reason,
    T minAngleDeg,
    T maxAngleDeg)
{
    if (mp.nPatches() == 0) {
        if (reason) *reason = "empty multipatch";
        return false;
    }

    const bool hasTopo = (mp.nInterfaces() > 0 || mp.nBoundary() > 0) &&
                         (2 * mp.nInterfaces() + mp.nBoundary() == 4 * mp.nPatches());
    gsMultiPatch<T> mpTopo;
    const gsMultiPatch<T>& topoMp = hasTopo ? mp : mpTopo;
    if (!hasTopo) {
        mpTopo = mp;
        mpTopo.computeTopology(std::max(T(1e-4), bdrTol), false, false);
        if (2 * mpTopo.nInterfaces() + mpTopo.nBoundary() != 4 * mpTopo.nPatches()) {
            if (reason) *reason = "non-conforming topology (hanging nodes or unmatched sides)";
            return false;
        }
    }

    T bdrErr = 0;
    if (!verifyBoundaryMatch(dom, mp, bdrTol, &bdrErr)) {
        if (std::getenv("GISMO_SPLITTER_DEBUG")) {
            std::cerr << "[DBG verifyBoundaryMatch FAIL] bdrErr=" << bdrErr << " tol=" << bdrTol << "\n";
            for (size_t p = 0; p < topoMp.nPatches(); ++p) {
                for (int s = 1; s <= 4; ++s) {
                    boxSide bs = static_cast<boxSide>(s);
                    if (!topoMp.isInterface(patchSide(p, bs))) {
                        auto bc = topoMp.patch(p).boundary(bs);
                        gsMatrix<T> pr = bc->parameterRange(), u(1, 3), val;
                        u << pr(0,0), T(0.5)*(pr(0,0)+pr(0,1)), pr(0,1);
                        bc->eval_into(u, val);
                        std::cerr << "  bnd p=" << p << " s=" << s << ": (" << val(0,0) << "," << val(1,0) << ")->(" << val(0,1) << "," << val(1,1) << ")->(" << val(0,2) << "," << val(1,2) << ")\n";
                    }
                }
            }
        }
        if (reason) *reason = "boundary mismatch";
        return false;
    }

    gsMatrix<T> u(2, 1), val;
    const T uv[4][2] = {{0, 0}, {1, 0}, {1, 1}, {0, 1}};
    size_t nP = mp.nPatches();
    std::vector<std::vector<gsPoint2D<T>>> polys(nP);

    for (size_t p = 0; p < nP; ++p) {
        if (!internal_splitter::isPatchDetSingleSign(mp.patch(p), 21, T(1e-7))) {
            if (reason) *reason = "detJ zero or sign flip within patch";
            return false;
        }
        gsPoint2D<T> c[4];
        for (int k = 0; k < 4; ++k) {
            u(0, 0) = uv[k][0]; u(1, 0) = uv[k][1];
            mp.patch(p).eval_into(u, val);
            c[k] = gsPoint2D<T>(val(0, 0), val(1, 0));
        }
        for (int k = 0; k < 4; ++k) {
            int kn = (k + 1) % 4;
            if (std::hypot(c[kn].x - c[k].x, c[kn].y - c[k].y) < T(1e-8)) {
                if (reason) *reason = "non-convex or degenerate quad corners";
                return false;
            }
        }
        T angs[4], minA = 0, maxA = 180;
        bool degen = false;
        internal_splitter::evalPatchCornerAnglesDeg(mp.patch(p), angs, minA, maxA, degen);
        if (degen || minA < minAngleDeg || maxA > maxAngleDeg) {
            if (reason) {
                std::ostringstream ss;
                ss << "corner tangent angle out of bounds on patch " << p 
                   << " (degen=" << degen << " minA=" << minA << " maxA=" << maxA 
                   << " bounds=[" << minAngleDeg << "," << maxAngleDeg << "])";
                *reason = ss.str();
            }
            return false;
        }
        int globalNSide = 24;
        for (size_t p2 = 0; p2 < nP; ++p2) {
            if (const auto* tb2 = dynamic_cast<const gsTensorBSpline<2, T>*>(&mp.patch(p2)))
                globalNSide = std::max(globalNSide, std::min(64, 2 * std::max<int>(tb2->basis().component(0).size(), tb2->basis().component(1).size())));
        }
        polys[p] = internal_splitter::samplePatchBoundaryPoly(mp.patch(p), globalNSide);
        // Only apply straight-chord corner convexity to patches with near-linear boundaries
        // (curved spokes / curved hole arcs can validly curve past straight-chord collinearity).
        T maxEdgeDev = 0;
        {
            const int perSide = static_cast<int>(polys[p].size()) / 4;
            for (int k = 0; k < 4 && perSide > 1; ++k) {
                const gsPoint2D<T>& A = c[k];
                const gsPoint2D<T>& B = c[(k + 1) % 4];
                T len = std::hypot(B.x - A.x, B.y - A.y);
                if (len <= T(1e-12)) continue;
                for (int s = 1; s < perSide; ++s) {
                    const gsPoint2D<T>& q = polys[p][k * perSide + s];
                    T dist = std::abs((B.x - A.x) * (q.y - A.y) - (B.y - A.y) * (q.x - A.x)) / len;
                    maxEdgeDev = std::max(maxEdgeDev, dist);
                }
            }
        }
        if (maxEdgeDev < T(1e-4) && !internal_splitter::isStrictlyConvexQuad(c)) {
            if (reason) *reason = "non-convex or degenerate quad corners";
            return false;
        }
        int m = static_cast<int>(polys[p].size());
        for (int i = 0; i < m; ++i) {
            for (int j = i + 2; j < m; ++j) {
                if (i == 0 && j == m - 1) continue;
                if (internal_splitter::segProperIntersect(polys[p][i], polys[p][(i + 1) % m],
                                                          polys[p][j], polys[p][(j + 1) % m])) {
                    if (reason) *reason = "patch boundary self-intersection";
                    return false;
                }
            }
        }
    }

    struct PatchAABB { T xMin, xMax, yMin, yMax; };
    std::vector<PatchAABB> boxes(nP);
    for (size_t p = 0; p < nP; ++p) {
        T xMin = std::numeric_limits<T>::max(), xMax = -std::numeric_limits<T>::max();
        T yMin = std::numeric_limits<T>::max(), yMax = -std::numeric_limits<T>::max();
        for (const auto& pt : polys[p]) {
            xMin = std::min(xMin, pt.x); xMax = std::max(xMax, pt.x);
            yMin = std::min(yMin, pt.y); yMax = std::max(yMax, pt.y);
        }
        boxes[p] = {xMin - T(1e-7), xMax + T(1e-7), yMin - T(1e-7), yMax + T(1e-7)};
    }

    gsMatrix<T> uInner(2, 9), valInner;
    const T pts[3] = {T(0.25), T(0.5), T(0.75)};
    int colIn = 0;
    for (int a = 0; a < 3; ++a) {
        for (int b = 0; b < 3; ++b) {
            uInner(0, colIn) = pts[a]; uInner(1, colIn) = pts[b];
            ++colIn;
        }
    }
    std::vector<std::vector<gsPoint2D<T>>> innerPts(nP, std::vector<gsPoint2D<T>>(9));
    for (size_t p = 0; p < nP; ++p) {
        mp.patch(p).eval_into(uInner, valInner);
        for (int k = 0; k < 9; ++k) {
            innerPts[p][k] = gsPoint2D<T>(valInner(0, k), valInner(1, k));
        }
    }

    for (size_t p = 0; p < nP; ++p) {
        for (size_t q = p + 1; q < nP; ++q) {
            if (boxes[p].xMax < boxes[q].xMin || boxes[p].xMin > boxes[q].xMax ||
                boxes[p].yMax < boxes[q].yMin || boxes[p].yMin > boxes[q].yMax)
                continue;

            int mp1 = static_cast<int>(polys[p].size());
            int mq1 = static_cast<int>(polys[q].size());
            const int perP = mp1 / 4;
            const int perQ = mq1 / 4;
            bool sharedSide[4][4] = {{false}};
            if (perP > 0 && perQ > 0) {
                for (int sp = 0; sp < 4; ++sp) {
                    const auto& p0 = polys[p][sp * perP];
                    const auto& pm = polys[p][sp * perP + perP / 2];
                    const auto& p1 = polys[p][((sp + 1) % 4) * perP];
                    for (int sq = 0; sq < 4; ++sq) {
                        const auto& q0 = polys[q][sq * perQ];
                        const auto& qm = polys[q][sq * perQ + perQ / 2];
                        const auto& q1 = polys[q][((sq + 1) % 4) * perQ];
                        if (((p0.dist(q0) < T(1e-4) && p1.dist(q1) < T(1e-4)) ||
                             (p0.dist(q1) < T(1e-4) && p1.dist(q0) < T(1e-4))) &&
                            pm.dist(qm) < T(1e-4)) {
                            sharedSide[sp][sq] = true;
                        }
                    }
                }
            }
            for (int i = 0; i < mp1; ++i) {
                const int sp = (perP > 0) ? (i / perP) : 0;
                const auto& A = polys[p][i];
                const auto& B = polys[p][(i + 1) % mp1];
                for (int j = 0; j < mq1; ++j) {
                    if (perQ > 0 && sharedSide[sp][j / perQ]) continue;
                    const auto& C = polys[q][j];
                    const auto& D = polys[q][(j + 1) % mq1];
                    if (A.dist(C) < T(1e-5) || A.dist(D) < T(1e-5) ||
                        B.dist(C) < T(1e-5) || B.dist(D) < T(1e-5)) {
                        continue;
                    }
                    if (internal_splitter::segProperIntersect(A, B, C, D)) {
                        if (std::getenv("GISMO_MERGED_DEBUG") != nullptr) {
                            std::cerr << "  [DBG cross-isect] p=" << p << " (s=" << sp << ") q=" << q
                                      << " (s=" << (perQ > 0 ? j / perQ : 0) << ") A=(" << A.x << "," << A.y
                                      << ") B=(" << B.x << "," << B.y << ") C=(" << C.x << "," << C.y
                                      << ") D=(" << D.x << "," << D.y << ")\n";
                        }
                        if (reason) *reason = "cross-patch boundary intersection";
                        return false;
                    }
                }
            }
            for (int k = 0; k < 9; ++k) {
                if (internal_splitter::pointInPolyStrict(innerPts[p][k], polys[q])) {
                    if (reason) *reason = "patch overlap";
                    return false;
                }
            }
        }
    }

    return true;
}

template <class T>
bool gsPatchSplitter<T>::verifyBoundaryMatch(T tol, T* maxDistOut) const
{
    return verifyBoundaryMatch(m_domain, m_bestCandidate.multiPatch, tol, maxDistOut);
}

template <class T>
bool gsPatchSplitter<T>::verifyBoundaryMatch(
    const gsPlanarDomain<T>& dom,
    const gsMultiPatch<T>& mp,
    T tol,
    T* maxDistOut)
{
    if (mp.nPatches() == 0) return false;

    // 1. Collect all original domain curves
    std::vector<const gsGeometry<T>*> origCurves;
    const auto& outer = dom.outer();
    for (int i = 0; i < outer.size(); ++i) {
        origCurves.push_back(&outer.curve(i));
    }
    for (int h = 0; h < dom.numHoles(); ++h) {
        const auto& hole = dom.loop(h + 1);
        for (int j = 0; j < hole.size(); ++j) {
            origCurves.push_back(&hole.curve(j));
        }
    }

    // 2. Collect all MultiPatch exterior boundary curves
    std::vector<typename gsGeometry<T>::uPtr> mpBdrCurves;
    for (size_t p = 0; p < mp.nPatches(); ++p) {
        for (int s = 1; s <= 4; ++s) {
            boxSide side = static_cast<boxSide>(s);
            if (!mp.isInterface(patchSide(p, side))) {
                mpBdrCurves.push_back(mp.patch(p).boundary(side));
            }
        }
    }

    // Pre-sample curves once for fast closest-bracket lookup, then refine continuously via golden-section search
    struct SampledCurve {
        const gsGeometry<T>* crv;
        T u0, u1;
        std::vector<T> u, x, y;
        T xMin, xMax, yMin, yMax;
    };
    auto buildSampled = [](const gsGeometry<T>* crv, int nSeg = 40) -> SampledCurve {
        SampledCurve sc;
        sc.crv = crv;
        gsMatrix<T> pr = crv->parameterRange();
        sc.u0 = pr(0, 0); sc.u1 = pr(0, 1);
        sc.u.resize(nSeg + 1);
        sc.x.resize(nSeg + 1);
        sc.y.resize(nSeg + 1);
        gsMatrix<T> uMat(1, nSeg + 1), val;
        for (int i = 0; i <= nSeg; ++i) {
            T ui = sc.u0 + (sc.u1 - sc.u0) * (static_cast<T>(i) / static_cast<T>(nSeg));
            sc.u[i] = ui;
            uMat(0, i) = ui;
        }
        crv->eval_into(uMat, val);
        sc.xMin = std::numeric_limits<T>::max();
        sc.xMax = -std::numeric_limits<T>::max();
        sc.yMin = std::numeric_limits<T>::max();
        sc.yMax = -std::numeric_limits<T>::max();
        for (int i = 0; i <= nSeg; ++i) {
            sc.x[i] = val(0, i);
            sc.y[i] = val(1, i);
            sc.xMin = std::min(sc.xMin, sc.x[i]);
            sc.xMax = std::max(sc.xMax, sc.x[i]);
            sc.yMin = std::min(sc.yMin, sc.y[i]);
            sc.yMax = std::max(sc.yMax, sc.y[i]);
        }
        return sc;
    };

    std::vector<SampledCurve> origSampled, mpSampled;
    origSampled.reserve(origCurves.size());
    for (const auto* c : origCurves) origSampled.push_back(buildSampled(c));
    mpSampled.reserve(mpBdrCurves.size());
    for (const auto& c : mpBdrCurves) mpSampled.push_back(buildSampled(c.get()));

    auto aabbDistSq = [](T px, T py, const SampledCurve& sc) -> T {
        T dx = (px < sc.xMin) ? (sc.xMin - px) : ((px > sc.xMax) ? (px - sc.xMax) : T(0));
        T dy = (py < sc.yMin) ? (sc.yMin - py) : ((py > sc.yMax) ? (py - sc.yMax) : T(0));
        return dx * dx + dy * dy;
    };

    auto pointToCurvesDist = [&](T px, T py, const std::vector<SampledCurve>& targetCurves) -> T {
        T bestPolyD2 = std::numeric_limits<T>::max();
        int bestC = -1;
        int bestSeg = -1;

        for (size_t c = 0; c < targetCurves.size(); ++c) {
            const auto& sc = targetCurves[c];
            if (aabbDistSq(px, py, sc) >= bestPolyD2) continue;

            int nSeg = static_cast<int>(sc.u.size()) - 1;
            for (int i = 0; i < nSeg; ++i) {
                T ax = sc.x[i], ay = sc.y[i];
                T dx = sc.x[i + 1] - ax, dy = sc.y[i + 1] - ay;
                T L2 = dx * dx + dy * dy;
                T t = (L2 > T(1e-24)) ? std::max(T(0), std::min(T(1), ((px - ax) * dx + (py - ay) * dy) / L2)) : T(0);
                T qx = ax + t * dx - px, qy = ay + t * dy - py;
                T d2 = qx * qx + qy * qy;
                if (d2 < bestPolyD2) {
                    bestPolyD2 = d2;
                    bestC = static_cast<int>(c);
                    bestSeg = i;
                }
            }
        }

        if (bestC < 0) return std::numeric_limits<T>::max();

        // Refine on best curve via golden-section search
        const auto& sc = targetCurves[bestC];
        int nSeg = static_cast<int>(sc.u.size()) - 1;
        T du = (sc.u1 - sc.u0) / static_cast<T>(nSeg);
        T a = std::max(sc.u0, sc.u[bestSeg] - T(0.5) * du);
        T b = std::min(sc.u1, sc.u[bestSeg + 1] + T(0.5) * du);
        const T phi = T(0.6180339887498949);
        gsMatrix<T> uMat(1, 1), val;
        auto evalD2 = [&](T u) -> T {
            uMat(0, 0) = u;
            sc.crv->eval_into(uMat, val);
            T dx = val(0, 0) - px, dy = val(1, 0) - py;
            return dx * dx + dy * dy;
        };
        T c1 = b - phi * (b - a), c2 = a + phi * (b - a);
        T f1 = evalD2(c1), f2 = evalD2(c2);
        for (int iter = 0; iter < 18; ++iter) {
            if (f1 < f2) {
                b = c2; c2 = c1; f2 = f1;
                c1 = b - phi * (b - a); f1 = evalD2(c1);
            } else {
                a = c1; c1 = c2; f1 = f2;
                c2 = a + phi * (b - a); f2 = evalD2(c2);
            }
        }
        T bestFinal = std::min(bestPolyD2, std::min(f1, f2));
        return std::sqrt(bestFinal);
    };

    T maxDist = 0;

    // Forward check: every original boundary curve point must match MultiPatch boundary
    for (const auto& oSc : origSampled) {
        for (size_t i = 0; i < oSc.u.size(); ++i) {
            T d = pointToCurvesDist(oSc.x[i], oSc.y[i], mpSampled);
            if (d > maxDist) maxDist = d;
        }
    }

    // Backward check: every MultiPatch boundary point must match original domain curves
    for (const auto& bSc : mpSampled) {
        for (size_t i = 0; i < bSc.u.size(); ++i) {
            T d = pointToCurvesDist(bSc.x[i], bSc.y[i], origSampled);
            if (d > maxDist) maxDist = d;
        }
    }

    if (maxDistOut) *maxDistOut = maxDist;
    return (maxDist <= tol);
}

template <class T>
std::vector<int> gsPatchSplitter<T>::getVertexValences() const
{
    return analyzeQuality(m_bestCandidate.multiPatch).intValences;
}

// ============================================================================
// Export & Downstream Interface
// ============================================================================

template <class T>
const gsMultiPatch<T>& gsPatchSplitter<T>::getMultiPatch() const
{
    return m_bestCandidate.multiPatch;
}

template <class T>
gsMultiPatch<T> gsPatchSplitter<T>::getBilinearSkeleton() const
{
    gsMultiPatch<T> skeleton;
    const auto& mp = m_bestCandidate.multiPatch;

    gsKnotVector<T> kv_u(0.0, 1.0, 0, 2, 1, 1);
    gsKnotVector<T> kv_v(0.0, 1.0, 0, 2, 1, 1);

    for (size_t p = 0; p < mp.nPatches(); ++p) {
        const auto* tp = dynamic_cast<const gsTensorBSpline<2, T>*>(&mp.patch(p));
        if (!tp) continue;
        int nu = tp->basis().size(0);
        int nv = tp->basis().size(1);
        const auto& c = tp->coefs();

        gsMatrix<T> corner_coefs(4, 2);
        corner_coefs.row(0) = c.row(0);
        corner_coefs.row(1) = c.row(nu - 1);
        corner_coefs.row(2) = c.row((nv - 1) * nu);
        corner_coefs.row(3) = c.row(nu - 1 + (nv - 1) * nu);

        skeleton.addPatch(typename gsGeometry<T>::uPtr(
            new gsTensorBSpline<2, T>(kv_u, kv_v, corner_coefs)));
    }

    skeleton.computeTopology();
    return skeleton;
}

template <class T>
void gsPatchSplitter<T>::configureAsG1Generator(gsAsG1GeometryGenerator<T>& generator) const
{
    const auto& mp = m_bestCandidate.multiPatch;
    if (mp.nPatches() == 0) return;

    int V = static_cast<int>(mp.nPatches());
    generator.setValence(V);
    generator.setDegree(m_degree);

    std::vector<gsVector<T>> corners;
    for (size_t p = 0; p < mp.nPatches(); ++p) {
        const auto* tp = dynamic_cast<const gsTensorBSpline<2, T>*>(&mp.patch(p));
        if (!tp) continue;
        int nu = tp->basis().size(0);
        const auto& c = tp->coefs();
        gsVector<T> pt(2);
        pt << c(nu - 1, 0), c(nu - 1, 1);
        corners.push_back(pt);
    }

    generator.setBoundaryCorners(corners);
    generator.setUseVirtualParallelogram(true);
}

template <class T>
bool gsPatchSplitter<T>::postProcessSplit(gsMultiPatch<T>& mp)
{
    if (mp.nPatches() == 0) return false;
    bool changed = false;

    auto multiPatchMinDetDense = [](const gsMultiPatch<T>& m) -> T {
        T gMin = std::numeric_limits<T>::max();
        for (size_t p = 0; p < m.nPatches(); ++p)
            gMin = std::min(gMin, internal_splitter::evalPatchMinDetDense(m.patch(p), 49));
        return gMin;
    };

    if (m_hasDomain && !m_hasInputMp && mp.nPatches() > 2) {
        gsMultiPatch<T> cand(mp);
        if (reducePolychordStrips(cand) && cand.nPatches() < mp.nPatches()) {
            gsASG1LocalFitter<T>::unifyPatchOrientations(cand);
            cand.computeTopology(T(1e-4), false, false);
            for (size_t p = 0; p < cand.nPatches(); ++p) {
                if (internal_splitter::evalPatchMinDetDense(cand.patch(p), 49) <= T(1e-8)) {
                    if (auto* tb = dynamic_cast<gsTensorBSpline<2, T>*>(&cand.patch(p)))
                        internal_splitter::quickHarmonicPatchRepair(*tb, T(1e-7));
                }
            }
            T minDet = multiPatchMinDetDense(cand);
            T maxBdrDist = 0;
            bool bdrOk = verifyBoundaryMatch(m_domain, cand, T(1e-4), &maxBdrDist);
            if (minDet > T(1e-8) && bdrOk) {
                mp = give(cand);
                changed = true;
            }
        }
    }

    gsMultiPatch<T> preOpt(mp);
    T preMinDet = multiPatchMinDetDense(preOpt);

    // Phase 1: Jointly optimize curved interior interfaces, 2D interior vertices, and CAD boundary splits
    if (mp.nPatches() >= 2) {
        gsMultiPatch<T> cCand(mp);
        if (optimizeCoupledInterfacesAndInteriors(cCand, T(1e-7), 24, true)) {
            T cMinDet = multiPatchMinDetDense(cCand);
            T maxBdrDist = 0;
            bool bdrOk = (!m_hasDomain) || verifyBoundaryMatch(m_domain, cCand, T(1e-4), &maxBdrDist);
            if (bdrOk && ((preMinDet > T(1e-8) && cMinDet >= T(0.92) * preMinDet) ||
                          (preMinDet <= T(1e-8) && cMinDet > preMinDet))) {
                mp = give(cCand);
                changed = true;
            }
        }
    }

    // Repair any patch with a hidden micro-fold on the dense 49x49 grid
    for (size_t p = 0; p < mp.nPatches(); ++p) {
        T pMin = internal_splitter::evalPatchMinDetDense(mp.patch(p), 49);
        if (pMin <= T(1e-7)) {
            if (auto* tb = dynamic_cast<gsTensorBSpline<2, T>*>(&mp.patch(p))) {
                gsMatrix<T> saveC = tb->coefs();
                internal_splitter::quickHarmonicPatchRepair(*tb, T(1e-7));
                T pAfter = internal_splitter::evalPatchMinDetDense(*tb, 49);
                if (pAfter < pMin)
                    tb->coefs() = saveC;
                else if (pAfter > pMin)
                    changed = true;
            }
        }
    }

    gsASG1LocalFitter<T>::unifyPatchOrientations(mp);
    mp.computeTopology(T(1e-4), false, false);
    T postMinDet = multiPatchMinDetDense(mp);
    if (preMinDet > T(1e-8) && postMinDet <= T(1e-8)) {
        mp = give(preOpt);
    }
    return changed;
}

template <class T>
void gsPatchSplitter<T>::exportParaviewMultiPatch(const gsMultiPatch<T>& mp,
                                                  const std::string& filename_prefix,
                                                  unsigned samplesPerEdge,
                                                  bool exportKnots)
{
    const size_t nP = mp.nPatches();
    if (nP == 0) return;

    gsMultiPatch<T> workMp(mp);
    workMp.computeTopology(T(1e-4), false, false);

    // Map (patch, side) -> true if exterior/hole domain boundary
    std::set<std::pair<int, int>> domainBoundarySides;
    for (auto it = workMp.bBegin(); it != workMp.bEnd(); ++it) {
        domainBoundarySides.insert({static_cast<int>(it->patch), static_cast<int>(it->side().index())});
    }

    // Choose smooth adaptive sampling per edge
    unsigned baseSegs = 24;
    if (samplesPerEdge > 0) {
        baseSegs = std::max(8u, std::min(64u, samplesPerEdge));
    } else if (nP > 150) {
        baseSegs = 12;
    } else if (nP > 64) {
        baseSegs = 16;
    } else if (nP > 24) {
        baseSegs = 20;
    }

    std::vector<unsigned> patchNu(nP, baseSegs), patchNv(nP, baseSegs);
    if (samplesPerEdge == 0) {
        auto edgeNeedsBoost = [](const gsGeometry<T>& patch, bool alongU, T fixedParam) -> bool {
            gsMatrix<T> ab = patch.parameterRange();
            T u0 = ab(0, 0), u1 = ab(0, 1), v0 = ab(1, 0), v1 = ab(1, 1);
            const int M = 9;
            gsMatrix<T> uv(2, M);
            for (int k = 0; k < M; ++k) {
                T s = static_cast<T>(k) / static_cast<T>(M - 1);
                if (alongU) {
                    uv(0, k) = u0 + s * (u1 - u0);
                    uv(1, k) = v0 + fixedParam * (v1 - v0);
                } else {
                    uv(0, k) = u0 + fixedParam * (u1 - u0);
                    uv(1, k) = v0 + s * (v1 - v0);
                }
            }
            gsMatrix<T> pts = patch.eval(uv);
            T x0 = pts(0, 0), y0 = pts(1, 0);
            T x1 = pts(0, M - 1), y1 = pts(1, M - 1);
            T chord = std::hypot(x1 - x0, y1 - y0);
            if (chord <= T(1e-12)) return false;
            T maxDev = T(0);
            for (int k = 1; k < M - 1; ++k) {
                T dx = pts(0, k) - x0, dy = pts(1, k) - y0;
                T dist = std::abs(dx * (y1 - y0) - dy * (x1 - x0)) / chord;
                if (dist > maxDev) maxDev = dist;
            }
            return (maxDev > T(4e-3) * chord);
        };
        for (size_t p = 0; p < nP; ++p) {
            const auto& patch = workMp.patch(p);
            if (edgeNeedsBoost(patch, true, T(0)) || edgeNeedsBoost(patch, true, T(1)))
                patchNu[p] = std::max(patchNu[p], 20u);
            if (edgeNeedsBoost(patch, false, T(0)) || edgeNeedsBoost(patch, false, T(1)))
                patchNv[p] = std::max(patchNv[p], 20u);
        }
    }

    // 1. Build merged VTU UnstructuredGrid of smooth quads with PatchID, JacobianDet, ScaledJacobian
    std::vector<T> ptCoords;
    std::vector<int> quadConn;
    std::vector<int> cellPatchId;
    std::vector<T> cellDetJ;
    std::vector<T> cellScaledJ;

    // 2. Build merged VTP PolyData of smooth patch boundary polylines
    std::vector<T> bndCoords;
    std::vector<std::vector<int>> bndLines;
    std::vector<int> bndPatchId;
    std::vector<int> bndIsDomain;

    size_t estPts = 0, estCells = 0;
    for (size_t p = 0; p < nP; ++p) {
        estPts += (patchNu[p] + 1) * (patchNv[p] + 1);
        estCells += patchNu[p] * patchNv[p];
    }
    ptCoords.reserve(estPts * 3);
    quadConn.reserve(estCells * 4);
    cellPatchId.reserve(estCells);
    cellDetJ.reserve(estCells);
    cellScaledJ.reserve(estCells);

    for (size_t p = 0; p < nP; ++p) {
        const auto& patch = workMp.patch(p);
        const unsigned nu = patchNu[p];
        const unsigned nv = patchNv[p];
        const unsigned nx = nu + 1;
        const unsigned ny = nv + 1;

        gsMatrix<T> ab = patch.parameterRange();
        T u0 = ab(0, 0), u1 = ab(0, 1);
        T v0 = ab(1, 0), v1 = ab(1, 1);

        gsMatrix<T> uv(2, nx * ny);
        for (unsigned j = 0; j < ny; ++j) {
            T v = v0 + (v1 - v0) * static_cast<T>(j) / static_cast<T>(nv);
            for (unsigned i = 0; i < nx; ++i) {
                T u = u0 + (u1 - u0) * static_cast<T>(i) / static_cast<T>(nu);
                uv(0, i + j * nx) = u;
                uv(1, i + j * nx) = v;
            }
        }
        gsMatrix<T> evalPts = patch.eval(uv);
        const int gDim = static_cast<int>(evalPts.rows());
        const int basePt = static_cast<int>(ptCoords.size() / 3);
        for (index_t k = 0; k < evalPts.cols(); ++k) {
            ptCoords.push_back(evalPts(0, k));
            ptCoords.push_back(gDim > 1 ? evalPts(1, k) : T(0));
            ptCoords.push_back(gDim > 2 ? evalPts(2, k) : T(0));
        }

        // Evaluate Jacobians at cell midpoints
        gsMatrix<T> uvMid(2, nu * nv);
        for (unsigned j = 0; j < nv; ++j) {
            T vm = v0 + (v1 - v0) * (static_cast<T>(j) + T(0.5)) / static_cast<T>(nv);
            for (unsigned i = 0; i < nu; ++i) {
                T um = u0 + (u1 - u0) * (static_cast<T>(i) + T(0.5)) / static_cast<T>(nu);
                uvMid(0, i + j * nu) = um;
                uvMid(1, i + j * nu) = vm;
            }
        }
        gsMatrix<T> jacAll = patch.jacobian(uvMid);

        for (unsigned j = 0; j < nv; ++j) {
            for (unsigned i = 0; i < nu; ++i) {
                int c00 = basePt + static_cast<int>(i + j * nx);
                int c10 = basePt + static_cast<int>((i + 1) + j * nx);
                int c11 = basePt + static_cast<int>((i + 1) + (j + 1) * nx);
                int c01 = basePt + static_cast<int>(i + (j + 1) * nx);
                quadConn.push_back(c00);
                quadConn.push_back(c10);
                quadConn.push_back(c11);
                quadConn.push_back(c01);

                const index_t cid = static_cast<index_t>(i + j * nu);
                T j00 = jacAll(0, 2 * cid);
                T j01 = jacAll(0, 2 * cid + 1);
                T j10 = jacAll(1, 2 * cid);
                T j11 = jacAll(1, 2 * cid + 1);
                T det = j00 * j11 - j10 * j01;
                T n1 = std::hypot(j00, j10);
                T n2 = std::hypot(j01, j11);
                T sj = (n1 * n2 > T(1e-18)) ? (det / (n1 * n2)) : T(0);

                cellPatchId.push_back(static_cast<int>(p));
                cellDetJ.push_back(det);
                cellScaledJ.push_back(sj);
            }
        }

        // Extract the 4 exact boundary polylines using the identical grid points
        auto pushBndPolyline = [&](const std::vector<index_t>& localIdx, boxSide side) {
            std::vector<int> line;
            line.reserve(localIdx.size());
            for (index_t k : localIdx) {
                int idx = static_cast<int>(bndCoords.size() / 3);
                bndCoords.push_back(evalPts(0, k));
                bndCoords.push_back(gDim > 1 ? evalPts(1, k) : T(0));
                bndCoords.push_back(gDim > 2 ? evalPts(2, k) : T(0));
                line.push_back(idx);
            }
            bndLines.push_back(std::move(line));
            bndPatchId.push_back(static_cast<int>(p));
            int isDom = domainBoundarySides.count({static_cast<int>(p), static_cast<int>(side.index())}) ? 1 : 0;
            bndIsDomain.push_back(isDom);
        };

        std::vector<index_t> southIdx(nx), northIdx(nx), westIdx(ny), eastIdx(ny);
        for (unsigned i = 0; i < nx; ++i) {
            southIdx[i] = static_cast<index_t>(i);
            northIdx[i] = static_cast<index_t>(i + nv * nx);
        }
        for (unsigned j = 0; j < ny; ++j) {
            westIdx[j] = static_cast<index_t>(j * nx);
            eastIdx[j] = static_cast<index_t>(nu + j * nx);
        }
        pushBndPolyline(southIdx, boundary::south);
        pushBndPolyline(eastIdx,  boundary::east);
        pushBndPolyline(northIdx, boundary::north);
        pushBndPolyline(westIdx,  boundary::west);
    }

    const std::string dir  = gsFileManager::getPath(filename_prefix);
    const std::string stem = gsFileManager::getBasename(filename_prefix);
    if (!dir.empty()) gsFileManager::mkdir(dir);

    const std::string vtuPath = filename_prefix + "_patches.vtu";
    const std::string vtpPath = filename_prefix + "_boundaries.vtp";
    const std::string pvdPath = filename_prefix + ".pvd";

    // Write <prefix>_patches.vtu
    {
        const size_t numPts = ptCoords.size() / 3;
        const size_t numCells = cellPatchId.size();
        std::ofstream f(vtuPath.c_str());
        f << std::fixed << std::setprecision(12);
        f << "<?xml version=\"1.0\"?>\n"
          << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n"
          << "  <UnstructuredGrid>\n"
          << "    <Piece NumberOfPoints=\"" << numPts << "\" NumberOfCells=\"" << numCells << "\">\n"
          << "      <Points>\n"
          << "        <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";
        for (size_t i = 0; i < numPts; ++i)
            f << ptCoords[3 * i] << " " << ptCoords[3 * i + 1] << " " << ptCoords[3 * i + 2] << "\n";
        f << "        </DataArray>\n"
          << "      </Points>\n"
          << "      <Cells>\n"
          << "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
        for (size_t c = 0; c < numCells; ++c)
            f << quadConn[4 * c] << " " << quadConn[4 * c + 1] << " "
              << quadConn[4 * c + 2] << " " << quadConn[4 * c + 3] << "\n";
        f << "        </DataArray>\n"
          << "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
        for (size_t c = 1; c <= numCells; ++c)
            f << (4 * c) << ((c % 16 == 0) ? "\n" : " ");
        f << "\n        </DataArray>\n"
          << "        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
        for (size_t c = 0; c < numCells; ++c)
            f << "9" << (((c + 1) % 32 == 0) ? "\n" : " ");
        f << "\n        </DataArray>\n"
          << "      </Cells>\n"
          << "      <CellData Scalars=\"PatchID\">\n"
          << "        <DataArray type=\"Int32\" Name=\"PatchID\" format=\"ascii\">\n";
        for (size_t c = 0; c < numCells; ++c)
            f << cellPatchId[c] << (((c + 1) % 32 == 0) ? "\n" : " ");
        f << "\n        </DataArray>\n"
          << "        <DataArray type=\"Float64\" Name=\"JacobianDet\" format=\"ascii\">\n";
        for (size_t c = 0; c < numCells; ++c)
            f << cellDetJ[c] << (((c + 1) % 8 == 0) ? "\n" : " ");
        f << "\n        </DataArray>\n"
          << "        <DataArray type=\"Float64\" Name=\"ScaledJacobian\" format=\"ascii\">\n";
        for (size_t c = 0; c < numCells; ++c)
            f << cellScaledJ[c] << (((c + 1) % 8 == 0) ? "\n" : " ");
        f << "\n        </DataArray>\n"
          << "      </CellData>\n"
          << "    </Piece>\n"
          << "  </UnstructuredGrid>\n"
          << "</VTKFile>\n";
    }

    // Write <prefix>_boundaries.vtp
    {
        const size_t numPts = bndCoords.size() / 3;
        const size_t numLines = bndLines.size();
        std::ofstream f(vtpPath.c_str());
        f << std::fixed << std::setprecision(12);
        f << "<?xml version=\"1.0\"?>\n"
          << "<VTKFile type=\"PolyData\" version=\"0.1\" byte_order=\"LittleEndian\">\n"
          << "  <PolyData>\n"
          << "    <Piece NumberOfPoints=\"" << numPts << "\" NumberOfVerts=\"0\" NumberOfLines=\""
          << numLines << "\" NumberOfStrips=\"0\" NumberOfPolys=\"0\">\n"
          << "      <Points>\n"
          << "        <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";
        for (size_t i = 0; i < numPts; ++i)
            f << bndCoords[3 * i] << " " << bndCoords[3 * i + 1] << " " << bndCoords[3 * i + 2] << "\n";
        f << "        </DataArray>\n"
          << "      </Points>\n"
          << "      <Lines>\n"
          << "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
        for (const auto& line : bndLines) {
            for (int v : line) f << v << " ";
            f << "\n";
        }
        f << "        </DataArray>\n"
          << "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
        size_t off = 0;
        for (size_t i = 0; i < numLines; ++i) {
            off += bndLines[i].size();
            f << off << (((i + 1) % 16 == 0) ? "\n" : " ");
        }
        f << "\n        </DataArray>\n"
          << "      </Lines>\n"
          << "      <CellData Scalars=\"IsDomainBoundary\">\n"
          << "        <DataArray type=\"Int32\" Name=\"IsDomainBoundary\" format=\"ascii\">\n";
        for (size_t i = 0; i < numLines; ++i)
            f << bndIsDomain[i] << (((i + 1) % 32 == 0) ? "\n" : " ");
        f << "\n        </DataArray>\n"
          << "        <DataArray type=\"Int32\" Name=\"PatchID\" format=\"ascii\">\n";
        for (size_t i = 0; i < numLines; ++i)
            f << bndPatchId[i] << (((i + 1) % 32 == 0) ? "\n" : " ");
        f << "\n        </DataArray>\n"
          << "      </CellData>\n"
          << "    </Piece>\n"
          << "  </PolyData>\n"
          << "</VTKFile>\n";
    }

    // Optional standalone IGA knot mesh (<prefix>_knots.vtp) — not added to .pvd by default to keep ParaView view uncluttered
    if (exportKnots) {
        std::vector<T> kCoords;
        std::vector<std::vector<int>> kLines;
        for (size_t p = 0; p < nP; ++p) {
            const auto* tb = dynamic_cast<const gsTensorBSpline<2, T>*>(&workMp.patch(p));
            if (!tb) continue;
            auto uKnots = tb->knots(0).unique();
            auto vKnots = tb->knots(1).unique();
            const unsigned nSamples = 16;
            T u0 = uKnots.front(), u1 = uKnots.back();
            T v0 = vKnots.front(), v1 = vKnots.back();
            for (size_t iu = 1; iu + 1 < uKnots.size(); ++iu) {
                gsMatrix<T> uv(2, nSamples + 1);
                for (unsigned k = 0; k <= nSamples; ++k) {
                    uv(0, k) = uKnots[iu];
                    uv(1, k) = v0 + (v1 - v0) * static_cast<T>(k) / static_cast<T>(nSamples);
                }
                gsMatrix<T> ev = tb->eval(uv);
                std::vector<int> ln;
                ln.reserve(nSamples + 1);
                for (unsigned k = 0; k <= nSamples; ++k) {
                    int idx = static_cast<int>(kCoords.size() / 3);
                    kCoords.push_back(ev(0, k));
                    kCoords.push_back(ev.rows() > 1 ? ev(1, k) : T(0));
                    kCoords.push_back(ev.rows() > 2 ? ev(2, k) : T(0));
                    ln.push_back(idx);
                }
                kLines.push_back(std::move(ln));
            }
            for (size_t iv = 1; iv + 1 < vKnots.size(); ++iv) {
                gsMatrix<T> uv(2, nSamples + 1);
                for (unsigned k = 0; k <= nSamples; ++k) {
                    uv(0, k) = u0 + (u1 - u0) * static_cast<T>(k) / static_cast<T>(nSamples);
                    uv(1, k) = vKnots[iv];
                }
                gsMatrix<T> ev = tb->eval(uv);
                std::vector<int> ln;
                ln.reserve(nSamples + 1);
                for (unsigned k = 0; k <= nSamples; ++k) {
                    int idx = static_cast<int>(kCoords.size() / 3);
                    kCoords.push_back(ev(0, k));
                    kCoords.push_back(ev.rows() > 1 ? ev(1, k) : T(0));
                    kCoords.push_back(ev.rows() > 2 ? ev(2, k) : T(0));
                    ln.push_back(idx);
                }
                kLines.push_back(std::move(ln));
            }
        }
        const std::string knotsPath = filename_prefix + "_knots.vtp";
        const size_t numPts = kCoords.size() / 3;
        const size_t numLines = kLines.size();
        std::ofstream f(knotsPath.c_str());
        f << std::fixed << std::setprecision(12);
        f << "<?xml version=\"1.0\"?>\n"
          << "<VTKFile type=\"PolyData\" version=\"0.1\" byte_order=\"LittleEndian\">\n"
          << "  <PolyData>\n"
          << "    <Piece NumberOfPoints=\"" << numPts << "\" NumberOfVerts=\"0\" NumberOfLines=\""
          << numLines << "\" NumberOfStrips=\"0\" NumberOfPolys=\"0\">\n"
          << "      <Points>\n"
          << "        <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";
        for (size_t i = 0; i < numPts; ++i)
            f << kCoords[3 * i] << " " << kCoords[3 * i + 1] << " " << kCoords[3 * i + 2] << "\n";
        f << "        </DataArray>\n"
          << "      </Points>\n"
          << "      <Lines>\n"
          << "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
        for (const auto& line : kLines) {
            for (int v : line) f << v << " ";
            f << "\n";
        }
        f << "        </DataArray>\n"
          << "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
        size_t off = 0;
        for (size_t i = 0; i < numLines; ++i) {
            off += kLines[i].size();
            f << off << (((i + 1) % 16 == 0) ? "\n" : " ");
        }
        f << "\n        </DataArray>\n"
          << "      </Lines>\n"
          << "    </Piece>\n"
          << "  </PolyData>\n"
          << "</VTKFile>\n";
    }

    // Write clean ParaView Collection <prefix>.pvd with explicit timestep, group, name, part
    {
        std::ofstream pvd(pvdPath.c_str());
        pvd << "<?xml version=\"1.0\"?>\n"
            << "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\">\n"
            << "  <Collection>\n"
            << "    <DataSet timestep=\"0\" group=\"patches\" name=\"Patches\" part=\"0\" file=\"" << stem << "_patches.vtu\"/>\n"
            << "    <DataSet timestep=\"0\" group=\"patch_boundaries\" name=\"PatchBoundaries\" part=\"1\" file=\"" << stem << "_boundaries.vtp\"/>\n"
            << "  </Collection>\n"
            << "</VTKFile>\n";
    }
}

template <class T>
void gsPatchSplitter<T>::exportParaview(const std::string& filename_prefix, unsigned npts) const
{
    if (m_bestCandidate.multiPatch.nPatches() > 0) {
        unsigned edgeSamples = (npts == 900) ? 0u : static_cast<unsigned>(std::max(8.0, std::round(std::sqrt(static_cast<double>(npts)))));
        exportParaviewMultiPatch(m_bestCandidate.multiPatch, filename_prefix, edgeSamples, false);
    }
}

} // namespace gismo
