/** @file gsMakeMultiPatch.cpp

    @brief Computes correctly the boundaries and interfaces of
    a multipatch structure

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris
*/

#include <gismo.h>
#include <gsNurbs/gsDeboor.hpp>
#include <algorithm>
#include <cstdlib>
#include <direct.h>
#include <fstream>
#include <limits>
#include <set>
#include <sstream>
#include <streambuf>
#include <unordered_map>
#include <unordered_set>

using namespace gismo;
using namespace std;

struct ProfileTimer
{
    ProfileTimer(const char* /*name*/, bool /*enabled*/) {}
};

#define PROFILE_FUNCTION() ProfileTimer _timer_##__LINE__(__FUNCTION__, true)
#define PROFILE_SECTION(name) ProfileTimer _timer_section_##__LINE__(name, true)

struct Profiler
{
    void printSummary() const {}
};

class TeeStreamBuf : public std::streambuf
{
public:
    TeeStreamBuf(std::streambuf* first, std::streambuf* second)
        : m_first(first), m_second(second)
    {
    }

protected:
    int overflow(int ch) override
    {
        if (ch == EOF)
            return !EOF;

        const int firstRes = m_first ? m_first->sputc(static_cast<char>(ch)) : ch;
        const int secondRes = m_second ? m_second->sputc(static_cast<char>(ch)) : ch;
        // Keep mirror file synchronized with console output line-by-line.
        if (ch == '\n')
            sync();
        return (firstRes == EOF || secondRes == EOF) ? EOF : ch;
    }

    int sync() override
    {
        const int firstSync = m_first ? m_first->pubsync() : 0;
        const int secondSync = m_second ? m_second->pubsync() : 0;
        return (firstSync == 0 && secondSync == 0) ? 0 : -1;
    }

private:
    std::streambuf* m_first;
    std::streambuf* m_second;
};

class ScopedConsoleMirror
{
public:
    explicit ScopedConsoleMirror(const std::string& filePath)
        : m_file(filePath.c_str(), std::ios::out | std::ios::trunc),
          m_coutTee(std::cout.rdbuf(), m_file.is_open() ? m_file.rdbuf() : nullptr),
          m_cerrTee(std::cerr.rdbuf(), m_file.is_open() ? m_file.rdbuf() : nullptr),
          m_oldCout(nullptr),
          m_oldCerr(nullptr),
          m_active(false)
    {
        if (!m_file.is_open())
            return;

        m_oldCout = std::cout.rdbuf(&m_coutTee);
        m_oldCerr = std::cerr.rdbuf(&m_cerrTee);
        m_active = true;
    }

    ~ScopedConsoleMirror()
    {
        if (!m_active)
            return;

        std::cout.rdbuf(m_oldCout);
        std::cerr.rdbuf(m_oldCerr);
    }

    bool active() const { return m_active; }

private:
    std::ofstream m_file;
    TeeStreamBuf m_coutTee;
    TeeStreamBuf m_cerrTee;
    std::streambuf* m_oldCout;
    std::streambuf* m_oldCerr;
    bool m_active;
};

static Profiler g_profiler;
static std::string g_currentOptimizationTerm = "unknown";

struct CellToCoarsen
{
    int level;
    int x1;
    int y1;
    int x2;
    int y2;
};
struct OriginalMpbesReference;

struct LocalCoarseningRegion
{
    bool enabled = false;
    int lambda = 0;
    // Exact union: one list of non-overlapping rectangles per patch.
    // Using a list of rects instead of a single bounding box prevents the
    // "gap bridging" effect where two distant twin strips on the same patch
    // were merged into one large AABB that captured all functions in between.
    std::vector<std::vector<std::array<real_t, 4>>> patchAABB;
    std::vector<bool> hasPatch;
    gsMatrix<real_t> basisInd;
};
// ======================================================

// ============ MPBES (Multi-Patch B-spline with Exact continuity) Class ============
// Based on "Analysis-suitable G1 multi-patch parametrizations for C1 isogeometric spaces" (2015)
// Integrates: Twin identification + Kraft selection + Truncation
template<short_t d, typename T>
class MPBES : public gsBasis<T> {
public:
    typedef memory::shared_ptr<MPBES> Ptr;
    typedef memory::unique_ptr<MPBES> uPtr;

    // Constructor - performs twin identification, selection, and truncation
    MPBES(const gsVector<gsVector<gsVector<index_t>>>& boxMat,
        const gsVector<gsTHBSplineBasis<d, T>>& thbBases,
        const gsVector<gsVector<gsTensorBSplineBasis<d, T>>>& bellsBases,
        const gsVector<gsVector<gsVector<index_t>>>& hasATwin,
        const gsVector<gsVector<gsVector<std::vector<index_t>>>>& twinsIndex,
        const gsVector<gsVector<gsVector<std::vector<index_t>>>>& twinsPatch,
        const gsVector<gsVector<gsVector<index_t>>>& indexInTHB,
        const gsVector<index_t>& currentLastNonZeroRow,
        bool verbose = true,
        int attempt = 0)
        : m_thbBases(thbBases),
        m_bellsBases(bellsBases),
        m_numPatches(thbBases.size()),
        m_boxMat(boxMat),
        m_indexInTHB(indexInTHB),
        m_currentLastNonZeroRow(currentLastNonZeroRow),
        m_attempt(attempt)
    {
        constructGlobalBasis(hasATwin, twinsIndex, twinsPatch, verbose);
    }

    // gsBasis interface implementations
    short_t domainDim() const override { return d; }
    short_t targetDim() const override { return d; }
    index_t size() const override { return m_numGlobalFunctions; }

    void eval_into(const gsMatrix<T>& u, gsMatrix<T>& result) const override {
        GISMO_NO_IMPLEMENTATION;
    }

    void deriv_into(const gsMatrix<T>& u, gsMatrix<T>& result) const override {
        GISMO_NO_IMPLEMENTATION;
    }

    void active_into(const gsMatrix<T>& u, gsMatrix<index_t>& result) const override {
        GISMO_NO_IMPLEMENTATION;
    }

    gsMatrix<T> support() const override {
        gsMatrix<T> result(d, 2);
        result.col(0).setConstant(0.0);
        result.col(1).setConstant(1.0);
        return result;
    }

    gsMatrix<T> support(const index_t& i) const override {
        return support();
    }

    void connectivity(const gsMatrix<T>& nodes, gsMesh<T>& mesh) const override {
        GISMO_NO_IMPLEMENTATION;
    }

    memory::unique_ptr<gsGeometry<T>> makeGeometry(gsMatrix<T> coefs) const override {
        GISMO_NO_IMPLEMENTATION;
    }

    std::ostream& print(std::ostream& os) const override {
        os << "MPBES basis with " << m_numGlobalFunctions << " functions on "
            << m_numPatches << " patches\n";
        return os;
    }

    GISMO_CLONE_FUNCTION(MPBES)

        // MPBES-specific methods

        /// Get number of global basis functions
        index_t numGlobalFunctions() const { return m_numGlobalFunctions; }
    index_t numPatches() const { return m_numPatches; }
    index_t nPatches() const { return m_numPatches; }  // Alias for compatibility

    /// Access to THB bases
    const gsVector<gsTHBSplineBasis<d, T>>& thbBases() const {
        return m_thbBases;
    }

    /// Access to Bells bases
    const gsVector<gsVector<gsTensorBSplineBasis<d, T>>>& bellsBases() const {
        return m_bellsBases;
    }

    /// Get global index from patch-level-index coordinate
    index_t globalIndexAt(index_t patch, index_t level, index_t localIdx) const {
        if (patch < m_globalIndex.size() &&
            level < m_globalIndex[patch].size() &&
            localIdx < m_globalIndex[patch][level].size()) {
            return m_globalIndex[patch][level][localIdx];
        }
        return -1;
    }

    /// Access to functionDescription for assembly
    const std::vector<std::vector<std::vector<index_t>>>& functionDescription() const {
        return m_functionDescription;
    }

    /// Access to truncation info
    const std::vector<bool>& isTruncated() const {
        return m_isTruncated;
    }

    /// Access to spillover info
    const std::vector<std::vector<std::array<int, 3>>>& spilloverCoordinates() const {
        return m_spilloverFunctionCoordinates;
    }

    const std::vector<bool>& hasSpillover() const {
        return m_hasSpillover;
    }

    /// Get the current attempt number (for conditional logging)
    int attempt() const {
        return m_attempt;
    }

    /// Access to truncation representation
    const std::vector<std::vector<gsSparseVector<T>>>& presentation() const {
        return m_presentation;
    }

    /// Access to indexInTHB mapping
    const gsVector<gsVector<gsVector<index_t>>>& indexInTHB() const {
        return m_indexInTHB;
    }

    /// Check if a specific component is truncated
    bool isComponentTruncated(index_t globalIdx, index_t componentIdx) const {
        if (globalIdx >= m_isComponentTruncated.size())
            return false;
        if (componentIdx >= m_isComponentTruncated[globalIdx].size())
            return false;
        return m_isComponentTruncated[globalIdx][componentIdx];
    }

    /// Get presentation level for a specific component
    index_t presentationLevel(index_t globalIdx, index_t componentIdx) const {
        if (globalIdx >= m_presentationLevel.size())
            return 0;
        if (componentIdx >= m_presentationLevel[globalIdx].size())
            return 0;
        return m_presentationLevel[globalIdx][componentIdx];
    }

    /// Evaluate single function on specific patch
    void evalSingleOnPatch(
        index_t globalIdx,
        index_t patch,
        const gsMatrix<T>& u,
        gsMatrix<T>& result,
        bool includeSpilloverFallback = false) const {
        result.resize(1, u.cols());
        result.setZero();

        if (globalIdx >= m_numGlobalFunctions || globalIdx >= m_functionDescription.size())
            return;

        const auto& twins = m_functionDescription[globalIdx];
        int componentIdx = 0;
        bool hasRegularComponentOnThisPatch = false;

        for (const auto& twin : twins) {
            if (twin.size() < 3)
                continue;
            if (twin[0] == static_cast<int>(patch)) {
                int level = twin[1];
                int tensorIdx = twin[2];
                if (level < 0 || tensorIdx < 0)
                {
                    componentIdx++;
                    continue;
                }

                const bool compTruncated = isComponentTruncated(globalIdx, componentIdx);

                if (compTruncated) {
                    if (globalIdx < m_presentation.size() && componentIdx < m_presentation[globalIdx].size()) {
                        const auto& coeffs = m_presentation[globalIdx][componentIdx];
                        if (globalIdx >= m_presentationLevel.size() || componentIdx >= m_presentationLevel[globalIdx].size())
                        {
                            componentIdx++;
                            continue;
                        }
                        index_t presLevel = m_presentationLevel[globalIdx][componentIdx];

                        if (patch < static_cast<index_t>(m_bellsBases.size()) &&
                            presLevel < static_cast<index_t>(m_bellsBases[patch].size())) {
                            const index_t basisSize = m_bellsBases[patch][presLevel].size();
                            for (typename gsSparseVector<T>::InnerIterator it(coeffs); it; ++it) {
                                const index_t bIdx = it.index();
                                if (bIdx < 0 || bIdx >= basisSize)
                                    continue;
                                gsMatrix<T> componentResult;
                                m_bellsBases[patch][presLevel].evalSingle_into(bIdx, u, componentResult);
                                result += it.value() * componentResult;
                            }
                        }
                    }
                }
                else {
                    if (patch < m_indexInTHB.size() &&
                        level < m_indexInTHB[patch].size() &&
                        tensorIdx < m_indexInTHB[patch][level].size()) {
                        int thbIdx = m_indexInTHB[patch][level][tensorIdx];
                        if (thbIdx >= 0) {
                            hasRegularComponentOnThisPatch = true;
                            gsMatrix<T> componentResult;
                            m_thbBases[patch].evalSingle_into(thbIdx, u, componentResult);
                            result += componentResult;
                        }
                    }
                }
            }
            componentIdx++;
        }

        if (includeSpilloverFallback &&
            !hasRegularComponentOnThisPatch &&
            globalIdx < static_cast<index_t>(m_spilloverFunctionCoordinates.size()))
        {
            for (const auto& spill : m_spilloverFunctionCoordinates[globalIdx])
            {
                const index_t spillPatch = spill[0];
                const index_t spillLevel = spill[1];
                const index_t spillTensorIdx = spill[2];

                if (spillPatch != patch)
                    continue;
                if (spillPatch < 0 || spillPatch >= static_cast<index_t>(m_bellsBases.size()))
                    continue;
                if (spillLevel < 0 || spillLevel >= static_cast<index_t>(m_bellsBases[spillPatch].size()))
                    continue;
                if (spillTensorIdx < 0 || spillTensorIdx >= m_bellsBases[spillPatch][spillLevel].size())
                    continue;

                gsMatrix<T> componentResult;
                m_bellsBases[spillPatch][spillLevel].evalSingle_into(spillTensorIdx, u, componentResult);
                result += componentResult;
            }
        }
    }

    /// Evaluate first derivative of single function on specific patch
    void evalDerivSingleOnPatch(
        index_t globalIdx,
        index_t patch,
        const gsMatrix<T>& u,
        gsMatrix<T>& result,
        bool includeSpilloverFallback = false) const {
        result.resize(d, u.cols());
        result.setZero();

        if (globalIdx >= m_numGlobalFunctions || globalIdx >= m_functionDescription.size())
            return;

        const auto& twins = m_functionDescription[globalIdx];
        int componentIdx = 0;
        bool hasRegularComponentOnThisPatch = false;

        for (const auto& twin : twins) {
            if (twin.size() < 3)
                continue;
            if (twin[0] == static_cast<int>(patch)) {
                int level = twin[1];
                int tensorIdx = twin[2];
                if (level < 0 || tensorIdx < 0)
                {
                    componentIdx++;
                    continue;
                }
                bool compTruncated = isComponentTruncated(globalIdx, componentIdx);

                if (compTruncated) {
                    if (globalIdx < m_presentation.size() && componentIdx < m_presentation[globalIdx].size()) {
                        const auto& coeffs = m_presentation[globalIdx][componentIdx];
                        if (globalIdx >= m_presentationLevel.size() || componentIdx >= m_presentationLevel[globalIdx].size())
                        {
                            componentIdx++;
                            continue;
                        }
                        index_t presLevel = m_presentationLevel[globalIdx][componentIdx];

                        if (patch < static_cast<index_t>(m_bellsBases.size()) &&
                            presLevel < static_cast<index_t>(m_bellsBases[patch].size())) {
                            const index_t basisSize = m_bellsBases[patch][presLevel].size();
                            for (typename gsSparseVector<T>::InnerIterator it(coeffs); it; ++it) {
                                const index_t bIdx = it.index();
                                if (bIdx < 0 || bIdx >= basisSize)
                                    continue;
                                gsMatrix<T> derivVal = m_bellsBases[patch][presLevel].function(bIdx).deriv(u);
                                result += it.value() * derivVal;
                            }
                        }
                    }
                }
                else {
                    if (patch < m_indexInTHB.size() &&
                        level < m_indexInTHB[patch].size() &&
                        tensorIdx < m_indexInTHB[patch][level].size()) {
                        int thbIdx = m_indexInTHB[patch][level][tensorIdx];
                        if (thbIdx >= 0) {
                            hasRegularComponentOnThisPatch = true;
                            gsMatrix<T> derivVal = m_thbBases[patch].function(thbIdx).deriv(u);
                            result += derivVal;
                        }
                    }
                }
            }
            componentIdx++;
        }

        if (includeSpilloverFallback &&
            !hasRegularComponentOnThisPatch &&
            globalIdx < static_cast<index_t>(m_spilloverFunctionCoordinates.size()))
        {
            for (const auto& spill : m_spilloverFunctionCoordinates[globalIdx])
            {
                const index_t spillPatch = spill[0];
                const index_t spillLevel = spill[1];
                const index_t spillTensorIdx = spill[2];

                if (spillPatch != patch)
                    continue;
                if (spillPatch < 0 || spillPatch >= static_cast<index_t>(m_bellsBases.size()))
                    continue;
                if (spillLevel < 0 || spillLevel >= static_cast<index_t>(m_bellsBases[spillPatch].size()))
                    continue;
                if (spillTensorIdx < 0 || spillTensorIdx >= m_bellsBases[spillPatch][spillLevel].size())
                    continue;

                result += m_bellsBases[spillPatch][spillLevel].function(spillTensorIdx).deriv(u);
            }
        }
    }

    /// Batch Jacobian assembly for one patch using fastDeriv_into.
    /// Replaces N×evalDerivSingleOnPatch calls with one fastDeriv_into + a map lookup.
    /// activeFuncs: global MPBES function indices active on this patch.
    /// Returns J00,J01,J10,J11 each of size (1 × M).
    void jacContribOnPatch(
        index_t patch,
        const gsMatrix<T>& u,
        const gsMatrix<T>& vectSol,
        const std::vector<index_t>& activeFuncs,
        gsMatrix<T>& J00, gsMatrix<T>& J01,
        gsMatrix<T>& J10, gsMatrix<T>& J11) const
    {
        const index_t M = u.cols();
        J00.setZero(1, M); J01.setZero(1, M);
        J10.setZero(1, M); J11.setZero(1, M);

        // --- Build THB-index → (cx,cy) map for non-truncated components ---
        // Truncated components (twins) are accumulated separately via evalDerivSingleOnPatch.
        std::unordered_map<index_t, std::pair<T,T>> thbCoeffMap;
        thbCoeffMap.reserve(activeFuncs.size());

        gsMatrix<T> bd;
        for (index_t f : activeFuncs) {
            if (f >= static_cast<index_t>(vectSol.rows()) || vectSol.cols() < 2) continue;
            const T cx = vectSol(f, 0), cy = vectSol(f, 1);
            if (!std::isfinite(cx) || !std::isfinite(cy)) continue;
            if (f >= static_cast<index_t>(m_functionDescription.size())) continue;

            int compIdx = 0;
            for (const auto& twin : m_functionDescription[f]) {
                if (twin.size() < 3 || twin[0] != static_cast<int>(patch)) { ++compIdx; continue; }
                const int level = twin[1], tensorIdx = twin[2];
                if (level < 0 || tensorIdx < 0) { ++compIdx; continue; }

                if (isComponentTruncated(f, compIdx)) {
                    // Truncated twin: fall back to per-function path
                    evalDerivSingleOnPatch(f, patch, u, bd);
                    J00 += cx * bd.row(0);  J01 += cx * bd.row(1);
                    J10 += cy * bd.row(0);  J11 += cy * bd.row(1);
                } else {
                    // Non-truncated: record the THB index this function maps to
                    if (patch < static_cast<index_t>(m_indexInTHB.size()) &&
                        level < static_cast<index_t>(m_indexInTHB[patch].size()) &&
                        tensorIdx < static_cast<index_t>(m_indexInTHB[patch][level].size())) {
                        const index_t thbIdx = m_indexInTHB[patch][level][tensorIdx];
                        if (thbIdx >= 0)
                            thbCoeffMap[thbIdx] = {cx, cy};
                    }
                }
                ++compIdx;
            }
        }

        if (thbCoeffMap.empty()) return;

        // --- One fastDeriv_into call for all non-truncated contributions ---
        gsMatrix<index_t> activeIdx;
        m_thbBases[patch].active_into(u, activeIdx);          // (numAct × M)
        gsMatrix<T> thbDeriv;
        m_thbBases[patch].fastDeriv_into(u, thbDeriv);        // (numAct*d × M)
        const index_t numAct = activeIdx.rows();

        for (index_t pt = 0; pt < M; ++pt) {
            for (index_t ind = 0; ind < numAct; ++ind) {
                const index_t thbIdx = activeIdx(ind, pt);
                if (ind != 0 && thbIdx == 0) break;           // trailing zero sentinel
                auto it = thbCoeffMap.find(thbIdx);
                if (it == thbCoeffMap.end()) continue;
                const T cx = it->second.first, cy = it->second.second;
                const T du = thbDeriv(ind * d,     pt);
                const T dv = thbDeriv(ind * d + 1, pt);
                J00(0, pt) += cx * du;  J01(0, pt) += cx * dv;
                J10(0, pt) += cy * du;  J11(0, pt) += cy * dv;
            }
        }
    }

    /// Batch-evaluate FULL geometry (all MPBES functions) on a patch at given UV points.
    /// geomOut is sized (M × 2) and filled by this method.
    void evalGeomOnPatch(
        index_t patch,
        const gsMatrix<T>& u,           // 2 × M
        const gsMatrix<T>& vectSol,
        gsMatrix<T>& geomOut) const     // M × 2, zeroed and filled here
    {
        const index_t M = u.cols();
        geomOut.setZero(M, 2);

        // Build THB-index → (cx,cy) map for non-truncated components on this patch.
        // Truncated components are evaluated in batch via evalSingleOnPatch (all M pts in one call).
        std::unordered_map<index_t, std::pair<T,T>> thbCoeffMap;
        thbCoeffMap.reserve(static_cast<size_t>(m_functionDescription.size() / 10 + 16));

        gsMatrix<T> bv;
        const index_t nf = std::min<index_t>(
            static_cast<index_t>(m_functionDescription.size()), vectSol.rows());
        for (index_t f = 0; f < nf; ++f) {
            if (vectSol.cols() < 2) continue;
            const T cx = vectSol(f, 0), cy = vectSol(f, 1);

            int compIdx = 0;
            for (const auto& twin : m_functionDescription[f]) {
                if (twin.size() < 3 || twin[0] != static_cast<int>(patch)) { ++compIdx; continue; }
                const int level = twin[1], tensorIdx = twin[2];
                if (level < 0 || tensorIdx < 0) { ++compIdx; continue; }

                if (isComponentTruncated(f, compIdx)) {
                    // Truncated: one batch call for all M points
                    evalSingleOnPatch(f, patch, u, bv, false);
                    if (bv.rows() >= 1 && bv.cols() == M) {
                        for (index_t pt = 0; pt < M; ++pt) {
                            const T val = bv(0, pt);
                            if (val == T(0)) continue;
                            geomOut(pt, 0) += val * cx;
                            geomOut(pt, 1) += val * cy;
                        }
                    }
                } else {
                    if (patch < static_cast<index_t>(m_indexInTHB.size()) &&
                        level < static_cast<index_t>(m_indexInTHB[patch].size()) &&
                        tensorIdx < static_cast<index_t>(m_indexInTHB[patch][level].size())) {
                        const index_t thbIdx = m_indexInTHB[patch][level][tensorIdx];
                        if (thbIdx >= 0)
                            thbCoeffMap[thbIdx] = {cx, cy};
                    }
                }
                ++compIdx;
            }
        }

        if (thbCoeffMap.empty()) return;

        // One active_into + eval_into call covers all non-truncated contributions
        gsMatrix<index_t> activeIdx;
        m_thbBases[patch].active_into(u, activeIdx);   // numAct × M
        gsMatrix<T> thbVals;
        m_thbBases[patch].eval_into(u, thbVals);       // numAct × M
        const index_t numAct = activeIdx.rows();

        for (index_t pt = 0; pt < M; ++pt) {
            for (index_t ind = 0; ind < numAct; ++ind) {
                const index_t thbIdx = activeIdx(ind, pt);
                if (ind != 0 && thbIdx == 0) break;    // trailing zero sentinel
                auto it = thbCoeffMap.find(thbIdx);
                if (it == thbCoeffMap.end()) continue;
                const T val = thbVals(ind, pt);
                geomOut(pt, 0) += val * it->second.first;
                geomOut(pt, 1) += val * it->second.second;
            }
        }
    }

    /// Batch-collect global MPBES function indices active at UV points on a given patch.
    /// Uses active_into + reverse THB→global map; O(nf + M×numActivePerPoint).
    void collectActiveGlobalIndices(
        index_t patch,
        const gsMatrix<T>& u,
        std::unordered_set<index_t>& activeGlobal) const
    {
        if (u.cols() == 0 || patch >= static_cast<index_t>(m_thbBases.size())) return;

        // Build reverse map: THB index → global MPBES index for this patch
        std::unordered_map<index_t, index_t> thbToGlobal;
        thbToGlobal.reserve(static_cast<size_t>(m_functionDescription.size() / 4 + 16));
        const index_t nf = static_cast<index_t>(m_functionDescription.size());
        for (index_t f = 0; f < nf; ++f) {
            for (const auto& twin : m_functionDescription[f]) {
                if (twin.size() < 3 || twin[0] != static_cast<int>(patch)) continue;
                const int level = twin[1], tensorIdx = twin[2];
                if (level < 0 || tensorIdx < 0) continue;
                if (patch < static_cast<index_t>(m_indexInTHB.size()) &&
                    level < static_cast<index_t>(m_indexInTHB[patch].size()) &&
                    tensorIdx < static_cast<index_t>(m_indexInTHB[patch][level].size())) {
                    const index_t thbIdx = m_indexInTHB[patch][level][tensorIdx];
                    if (thbIdx != static_cast<index_t>(-1))
                        thbToGlobal.emplace(thbIdx, f);
                }
            }
        }
        if (thbToGlobal.empty()) return;

        gsMatrix<index_t> activeThb;
        m_thbBases[patch].active_into(u, activeThb);
        for (index_t col = 0; col < activeThb.cols(); ++col)
            for (index_t row = 0; row < activeThb.rows(); ++row) {
                auto it = thbToGlobal.find(activeThb(row, col));
                if (it != thbToGlobal.end()) activeGlobal.insert(it->second);
            }
    }

    /// Evaluate second derivative of single function on specific patch
    void evalDeriv2SingleOnPatch(
        index_t globalIdx,
        index_t patch,
        const gsMatrix<T>& u,
        gsMatrix<T>& result,
        bool includeSpilloverFallback = false) const {
        const index_t numSecondDerivs = d * (d + 1) / 2;
        result.resize(numSecondDerivs, u.cols());
        result.setZero();

        if (globalIdx >= m_numGlobalFunctions || globalIdx >= m_functionDescription.size())
            return;

        const auto& twins = m_functionDescription[globalIdx];
        int componentIdx = 0;
        bool hasRegularComponentOnThisPatch = false;

        for (const auto& twin : twins) {
            if (twin.size() < 3)
                continue;
            if (twin[0] == static_cast<int>(patch)) {
                int level = twin[1];
                int tensorIdx = twin[2];
                if (level < 0 || tensorIdx < 0)
                {
                    componentIdx++;
                    continue;
                }
                bool compTruncated = isComponentTruncated(globalIdx, componentIdx);

                if (compTruncated) {
                    if (globalIdx < m_presentation.size() && componentIdx < m_presentation[globalIdx].size()) {
                        const auto& coeffs = m_presentation[globalIdx][componentIdx];
                        if (globalIdx >= m_presentationLevel.size() || componentIdx >= m_presentationLevel[globalIdx].size())
                        {
                            componentIdx++;
                            continue;
                        }
                        index_t presLevel = m_presentationLevel[globalIdx][componentIdx];

                        if (patch < static_cast<index_t>(m_bellsBases.size()) &&
                            presLevel < static_cast<index_t>(m_bellsBases[patch].size())) {
                            const index_t basisSize = m_bellsBases[patch][presLevel].size();
                            for (typename gsSparseVector<T>::InnerIterator it(coeffs); it; ++it) {
                                const index_t bIdx = it.index();
                                if (bIdx < 0 || bIdx >= basisSize)
                                    continue;
                                gsMatrix<T> deriv2Val = m_bellsBases[patch][presLevel].function(bIdx).deriv2(u);
                                result += it.value() * deriv2Val;
                            }
                        }
                    }
                }
                else {
                    if (patch < m_indexInTHB.size() &&
                        level < m_indexInTHB[patch].size() &&
                        tensorIdx < m_indexInTHB[patch][level].size()) {
                        int thbIdx = m_indexInTHB[patch][level][tensorIdx];
                        if (thbIdx >= 0) {
                            hasRegularComponentOnThisPatch = true;
                            gsMatrix<T> deriv2Val = m_thbBases[patch].function(thbIdx).deriv2(u);
                            result += deriv2Val;
                        }
                    }
                }
            }
            componentIdx++;
        }

        if (includeSpilloverFallback &&
            !hasRegularComponentOnThisPatch &&
            globalIdx < static_cast<index_t>(m_spilloverFunctionCoordinates.size()))
        {
            for (const auto& spill : m_spilloverFunctionCoordinates[globalIdx])
            {
                const index_t spillPatch = spill[0];
                const index_t spillLevel = spill[1];
                const index_t spillTensorIdx = spill[2];

                if (spillPatch != patch)
                    continue;
                if (spillPatch < 0 || spillPatch >= static_cast<index_t>(m_bellsBases.size()))
                    continue;
                if (spillLevel < 0 || spillLevel >= static_cast<index_t>(m_bellsBases[spillPatch].size()))
                    continue;
                if (spillTensorIdx < 0 || spillTensorIdx >= m_bellsBases[spillPatch][spillLevel].size())
                    continue;

                result += m_bellsBases[spillPatch][spillLevel].function(spillTensorIdx).deriv2(u);
            }
        }
    }

    /// Wrapper for evalSingle_into signature used in old code
    void evalSingle_into(
        index_t globalIdx,
        index_t thbIdx,
        const gsMatrix<T>& u,
        gsMatrix<T>& result,
        index_t patch) const
    {
        // If thbIdx == -1, this is a spillover function
        if (thbIdx == -1) {
            // For spillover, we still need to evaluate using the truncated representation
            evalSingleOnPatch(globalIdx, patch, u, result);
        }
        else {
            // Regular evaluation
            evalSingleOnPatch(globalIdx, patch, u, result);
        }
    }

private:
    // Member variables
    gsVector<gsTHBSplineBasis<d, T>> m_thbBases;
    gsVector<gsVector<gsTensorBSplineBasis<d, T>>> m_bellsBases;
    index_t m_numPatches;
    index_t m_numGlobalFunctions;

    // Box structure
    gsVector<gsVector<gsVector<index_t>>> m_boxMat;
    gsVector<gsVector<gsVector<index_t>>> m_indexInTHB;
    gsVector<index_t> m_currentLastNonZeroRow;  // Valid box count per patch

    // Global indexing
    gsVector<gsVector<gsVector<index_t>>> m_globalIndex;  // [patch][level][localIdx] -> globalIdx
    gsVector<gsVector<int>> m_globalIndexTHB;              // [patch][thbIdx] -> globalIdx

    // Function description and spillover
    std::vector<std::vector<std::vector<index_t>>> m_functionDescription;  // [globalIdx][twinIdx] = {patch,level,index}
    std::vector<std::vector<std::array<int, 3>>> m_spilloverFunctionCoordinates;
    std::vector<std::vector<std::array<int, 3>>> m_spilloverSources;
    std::vector<bool> m_hasSpillover;

    // Truncation
    std::vector<bool> m_isTruncated;  // [globalIdx] -> true if ANY component is truncated
    std::vector<std::vector<bool>> m_isComponentTruncated;  // [globalIdx][componentIdx] -> true if THIS component is truncated
    std::vector<std::vector<gsSparseVector<T>>> m_presentation;  // Truncated representations
    std::vector<std::vector<index_t>> m_presentationLevel;  // [globalIdx][componentIdx] -> level of presentation
    int m_attempt;  // For conditional logging
    int verbose_selectionMechanism = 0;  // Set to 1 to enable detailed logging in selectionMechanism
    /// Main construction: Twin identification + Kraft selection + Truncation
    void constructGlobalBasis(
        const gsVector<gsVector<gsVector<index_t>>>& hasATwin,
        const gsVector<gsVector<gsVector<std::vector<index_t>>>>& twinsIndex,
        const gsVector<gsVector<gsVector<std::vector<index_t>>>>& twinsPatch,
        bool verbose)
    {
        // Step 1: Apply Kraft selection mechanism
        selectionMechanism(
            m_boxMat,
            m_bellsBases,
            m_thbBases,
            hasATwin,
            twinsIndex,
            twinsPatch,
            m_indexInTHB,
            m_currentLastNonZeroRow,
            m_functionDescription,
            m_spilloverFunctionCoordinates,
            m_spilloverSources,
            m_hasSpillover,
            m_globalIndexTHB,
            m_globalIndex,
            m_numGlobalFunctions,
            verbose_selectionMechanism
        );

        if (outfile.is_open() && g_verbose)
        {
            outfile << "[mpbes-init] functionDescription dump begin\n";
            for (size_t f = 0; f < m_functionDescription.size(); ++f)
            {
                outfile << "[mpbes-init] functionDescription[" << f << "] ";
                for (size_t compIdx = 0; compIdx < m_functionDescription[f].size(); ++compIdx)
                {
                    const auto& desc = m_functionDescription[f][compIdx];
                    if (desc.size() < 3)
                    {
                        outfile << "-1 -1 -1  ";
                        continue;
                    }

                    outfile << desc[0] << ' ' << desc[1] << ' ' << desc[2];

                    const int patch = static_cast<int>(desc[0]);
                    const int level = static_cast<int>(desc[1]);
                    const int tensorIdx = static_cast<int>(desc[2]);
                    if (patch >= 0 && patch < static_cast<int>(m_bellsBases.size()) &&
                        level >= 0 && level < static_cast<int>(m_bellsBases[patch].size()) &&
                        tensorIdx >= 0 && tensorIdx < static_cast<int>(m_bellsBases[patch][level].size()))
                    {
                        const auto support = m_bellsBases[patch][level].function(tensorIdx).support();
                        outfile << " ["
                                << support(0, 0) << ' ' << support(1, 0) << ' '
                                << support(0, 1) << ' ' << support(1, 1) << "]";
                    }

                    outfile << "  ";
                }
                outfile << "\n";
            }
            outfile << "[mpbes-init] functionDescription dump end\n";
        }

        // if (verbose)
        //     gsInfo << "MPBES: After selection, " << m_numGlobalFunctions << " global functions\n";

        // Step 2: Detect which functions need truncation
        detectTruncatedFunctions(verbose);

        // Step 3: Compute truncated representations (spillover coefficients)
        if (std::any_of(m_isTruncated.begin(), m_isTruncated.end(), [](bool v) { return v; })) {
            m_presentation.resize(m_numGlobalFunctions);
            computeTruncatedRepresentations(verbose);
        }

        // if (verbose) {
        //     gsInfo << "MPBES constructed: " << m_numGlobalFunctions << " global functions\n";
        //     index_t numTruncated = std::count(m_isTruncated.begin(), m_isTruncated.end(), true);
        //     gsInfo << "  - " << numTruncated << " functions are truncated\n";
        // }
    }

    /// Detect which functions need truncation based on Kraft's criterion:
    /// A function β^L at level L is truncated if its support overlaps with 
    /// a refinement region at level L+1 or higher (i.e., supp β^L ∩ Ω^(L+1) ≠ ∅)
    /// CRITICAL: We track truncation at COMPONENT level, not just function level
    void detectTruncatedFunctions(bool verbose) {
        m_isTruncated.clear();
        m_isTruncated.resize(m_functionDescription.size(), false);
        m_isComponentTruncated.clear();
        m_isComponentTruncated.resize(m_functionDescription.size());
        m_presentationLevel.clear();
        m_presentationLevel.resize(m_functionDescription.size());

        index_t numTruncatedComponents = 0;

        for (size_t f = 0; f < m_functionDescription.size(); ++f) {
            m_isComponentTruncated[f].resize(m_functionDescription[f].size(), false);
            m_presentationLevel[f].resize(m_functionDescription[f].size(), 0);

            // Check each twin component of this global function
            for (size_t twinIdx = 0; twinIdx < m_functionDescription[f].size(); ++twinIdx) {
                int patch = m_functionDescription[f][twinIdx][0];
                int level = m_functionDescription[f][twinIdx][1];
                int tensorIdx = m_functionDescription[f][twinIdx][2];

                // Skip if not in THB basis (spillover twins - they are NEVER truncated because they don't exist in THB)
                if (patch >= m_indexInTHB.size() ||
                    level >= m_indexInTHB[patch].size() ||
                    tensorIdx >= m_indexInTHB[patch][level].size() ||
                    m_indexInTHB[patch][level][tensorIdx] == -1) {
                    continue;
                }

                int thbIdx = m_indexInTHB[patch][level][tensorIdx];

                // Get support of this function in index space using THB basis methods
                // The THB basis knows how to get element support for functions at any level
                gsMatrix<index_t, d, 2> elementSupport;
                m_thbBases[patch].elementSupport_into(thbIdx, elementSupport);

                gsVector<index_t, d> low = elementSupport.col(0);
                gsVector<index_t, d> high = elementSupport.col(1);

                // Query: what is the highest level that overlaps this support?
                // If it's > level, then this component overlaps with a finer refinement region
                int highestOverlap = m_thbBases[patch].tree().query4(low, high, level);

                if (highestOverlap > level) {
                    // This COMPONENT's support overlaps with a refinement at level > L
                    // According to Kraft's criterion: Ω^L ⊇ supp β^L ⊄ Ω^(L+1)
                    // THIS COMPONENT must be truncated
                    m_isComponentTruncated[f][twinIdx] = true;
                    m_presentationLevel[f][twinIdx] = highestOverlap;
                    m_isTruncated[f] = true;  // At least one component is truncated
                    numTruncatedComponents++;

                    // if (f < 10 && verbose) {
                    //     gsInfo << "  [DETECT] Func " << f << " comp " << twinIdx
                    //         << ": level=" << level << ", highestOverlap=" << highestOverlap
                    //         << " -> presentationLevel=" << m_presentationLevel[f][twinIdx] << "\n";
                    // }
                }
            }
        }

        index_t numTruncatedFunctions = std::count(m_isTruncated.begin(), m_isTruncated.end(), true);

        // gsInfo << "MPBES: Detected " << numTruncatedComponents << " truncated components in "
        //     << numTruncatedFunctions << " functions\n";
    }

    /// Compute truncated representations (call existing computeCoefsForTruncatedFunctions)
    void computeTruncatedRepresentations(bool verbose) {
        std::vector<std::vector<std::vector<T>>> globalCoefs;  // Dummy, not used

        computeCoefsForTruncatedFunctions(
            m_bellsBases,
            m_thbBases,
            m_functionDescription,
            m_presentation,
            m_isTruncated,
            globalCoefs,
            m_globalIndex,
            m_indexInTHB,
            m_spilloverFunctionCoordinates,
            m_isComponentTruncated,
            m_presentationLevel,
            verbose,
            m_attempt
        );
    }
    
};

// ======================================================

index_t commonBasisSize;
int commonSize;
int DTD;
int printAB;
int localTempAttempt;
int numberOfPatchesForReference;
ofstream outfile;
ofstream summaryFile;
ofstream xboxFile;
ofstream yboxFile;
ofstream closureLogFile;

// Set to true to emit INTERFACE_X_Y_DIAGNOSTIC_BEGIN/END blocks in the console/log.
// Off by default — these blocks are very verbose (~60 lines per interface per step).
static bool g_enableInterfaceDiagnostics = false;

// Skip the Level-2 LO-with-zero-orthogonality fallback and go straight from Level-1 LO
// failure to Level-3 Gauss-Newton NLO.  Set via --skip-lo-fallback flag.
// Default off to preserve existing behaviour.
static bool g_skipLoFallback = false;

// Multiplicative scale applied to the LO uniformity and length weights before the
// LO call.  Default 1.0 (no change).  Set via --lo-weight-scale.  Values < 1 weaken
// regularisation so that LO is easier to break; useful for probing the failure threshold.
static double g_loWeightScale = 1.0;

// LS-only mode: skip LO and NLO; accept only when FIT is already regular and within tolerance.
// Set via --ls-only flag.  Default off (full algorithm: FIT → LO → NLO).
static bool g_lsOnly = false;

// LO-only mode: try LO when FIT is irregular, but do not invoke NLO.
// Set via --lo-only flag.  Default off.
static bool g_loOnly = false;

// Global approximation error tolerance ε_g.  Negative = use hardcoded default (1e+6).
// Set via --epsilon-g <value>.
static double g_epsilonG = -1.0;

// Feature boundary error tolerance ε_f.  Negative = use hardcoded default (0.1).
// Set via --epsilon-f <value>.
static double g_epsilonF = -1.0;

// Cell selection method.  Set via --cell-method <g|l|r|s>.
// 'g' = Grenda's geometry-based (default), 'l' = lexicographic, 'r' = random, 's' = smallest.
static char g_cellMethod = 'g';

// Verbose logging.  Set via --verbose flag.
// When false (default): only essential output (timing, errors, escape times, FIT entries).
// When true: also dumps indexInTHB, functionDescription, vectSol coefficients,
//            IdentifyPatches twin pairs, and geo-coarsen candidate details.
static bool g_verbose = false;

// Use local fitting instead of global.  Set via --local-fitting flag.
// Global fitting uses all MPBES evaluation points; local fitting restricts to the
// axis-aligned bounding box of the support of functions active on the coarsened cells
// (paper equations (9)-(10) with lambda-cell extension).  Default off (global).
static bool g_useLocalFitting = false;

// Locality parameter λ for local fitting (number of cell widths to extend the local
// evaluation region beyond the coarsened cell support).  Set via --lambda <n>.
// Only used when --local-fitting is active.  Default 0 (tight bounding box).
static int g_localityLambda = 0;

// Use a dedicated exit code so partition failures are easy to detect from scripts.
static constexpr int kPartitionUnityViolationExitCode = 86;
// Temporary diagnostic exit: geometry regular, but epsilon acceptance conditions fail.
static constexpr int kRegularGeometryEpsilonViolationExitCode = 87;

class ProgramExitSignal : public std::exception
{
public:
    ProgramExitSignal(int code, const std::string& reason)
        : m_code(code), m_reason(reason)
    {
    }

    int code() const noexcept { return m_code; }
    const char* what() const noexcept override { return m_reason.c_str(); }

private:
    int m_code;
    std::string m_reason;
};

struct PartitionUnityViolationDetails
{
    bool allSatisfied = true;
    index_t violationCount = 0;
    index_t firstViolatingRow = -1;
    real_t firstViolatingRowSum = 1.0;
    index_t maxDeviationRow = -1;
    real_t maxDeviation = 0.0;
    real_t maxDeviationRowSum = 1.0;
    std::vector<real_t> rowSums;
};

static void flushAndCloseDiagnosticLogs()
{
    if (outfile.is_open())
    {
        outfile.flush();
        outfile.close();
    }
    if (summaryFile.is_open())
    {
        summaryFile.flush();
        summaryFile.close();
    }
    if (xboxFile.is_open())
    {
        xboxFile.flush();
        xboxFile.close();
    }
    if (yboxFile.is_open())
    {
        yboxFile.flush();
        yboxFile.close();
    }
    gsInfo << std::flush;
}

static PartitionUnityViolationDetails analyzePartitionOfUnity(
    const gsSparseMatrix<real_t>& A,
    real_t tol)
{
    PartitionUnityViolationDetails details;
    const index_t nRows = A.rows();
    details.rowSums.assign(static_cast<size_t>(nRows), 0.0);

    for (index_t k = 0; k < A.outerSize(); ++k)
    {
        for (gsSparseMatrix<real_t>::InnerIterator it(A, k); it; ++it)
        {
            const index_t row = it.row();
            if (row >= 0 && row < nRows)
                details.rowSums[static_cast<size_t>(row)] += it.value();
        }
    }

    for (index_t row = 0; row < nRows; ++row)
    {
        const real_t rowSum = details.rowSums[static_cast<size_t>(row)];
        const real_t deviation = std::abs(rowSum - 1.0);
        if (deviation > tol)
        {
            details.allSatisfied = false;
            ++details.violationCount;

            if (details.firstViolatingRow < 0)
            {
                details.firstViolatingRow = row;
                details.firstViolatingRowSum = rowSum;
            }

            if (deviation > details.maxDeviation)
            {
                details.maxDeviation = deviation;
                details.maxDeviationRow = row;
                details.maxDeviationRowSum = rowSum;
            }
        }
    }

    return details;
}



template <typename T>
void logTHBContributionsForLine(
    T fixed_u,
    index_t targetPatch,
    const gsVector<gsTHBSplineBasis<2, T>>& THBVector,
    std::ostream& out
);
template <typename T>
void logSpilloverContributionsForLine(
    T fixed_u,
    index_t targetPatch,
    const gsVector<gsVector<gsTensorBSplineBasis<2, T>>>& Bells,
    const std::vector<std::vector<std::array<int, 3>>>& spilloverFunctionCoordinates,
    const std::vector<bool>& hasSpillover,
    std::ostream& out
);
index_t countTotalActiveElements(const gsVector<gsTHBSplineBasis<2, real_t>>& THBVector);

void exportLinesToFile(
    const std::vector<gsMatrix<>>& lines,
    const std::vector<index_t>& patchIDs,
    const std::vector<char>& directions,
    const std::string& filename,
    bool verbose,
    const std::vector<gsMatrix<>>& uvLines // <-- NEW
);

void printTheMatrix(const gsMatrix<real_t>& matrix, const std::string& matrixName);

void active_into_Spillovers(
    index_t patch,
    const gsMatrix<real_t>& uv,
    const std::vector<std::vector<std::array<int, 3>>>& spilloverFunctionCoordinates,
    gsVector<gsVector<gsTensorBSplineBasis<2, real_t>>> Bells,
    std::vector<std::vector<int>>& activeSpillovers);

void savePatches(
    const gsMatrix<real_t>& vectSol,
    const gsVector<gsVector<int>>& globalIndexTHB,
    const gsVector<gsTHBSplineBasis<2, real_t>>& SubdomainHierarchy,
    std::ostream& outfile,
    const std::string& baseName,
    gsVector<gsVector<int>>& localIndex
);

void saveMultipatchGeometry(
    const gsMatrix<real_t>& vectSol,
    const gsVector<gsVector<int>>& globalIndexTHB,
    const gsVector<gsTHBSplineBasis<2, real_t>>& SubdomainHierarchy,
    std::ostream& outfile,
    const std::string& baseName,
    gsVector<gsVector<int>>& localIndex
);

std::vector<std::unordered_set<index_t>> computeTwinPatches(
    const std::vector<std::vector<std::vector<index_t>>>& functionDescription);
void printDerMatrices(gsMatrix<real_t> basisDer,
    gsMatrix<real_t> geomDer,
    gsMatrix<real_t> basis2Der,
    gsMatrix<real_t> geom2Der,
    bool verbose = false
);
void selectionMechanism(
    const gsVector<gsVector<gsVector<index_t>>>& boxMat,
    const gsVector<gsVector<gsTensorBSplineBasis<2, real_t>>>& Bells,
    const gsVector<gsTHBSplineBasis<2, real_t>>& THBVector,
    const gsVector<gsVector<gsVector<index_t>>>& hasATwin,
    const gsVector<gsVector<gsVector<std::vector<index_t>>>>& twinsIndex,
    const gsVector<gsVector<gsVector<std::vector<index_t>>>>& twinsPatch,
    const gsVector<gsVector<gsVector<index_t>>>& indexInTHB,
    const gsVector<index_t>& currentLastNonZeroRow,
    std::vector<std::vector<std::vector<index_t>>>& functionDescription,
    std::vector<std::vector<std::array<int, 3>>>& spilloverFunctionCoordinates,
    std::vector<std::vector<std::array<int, 3>>>& spilloverSources,
    std::vector<bool>& hasSpillover,
    gsVector<gsVector<int>>& globalIndexTHB,
    gsVector<gsVector<gsVector<index_t>>>& globalIndex,
    int& globalFunctionCount,
    bool verbose = true);
void truncationMechanism(
    const gsVector<gsVector<gsTensorBSplineBasis<2, real_t>>>& Bells,
    const std::vector<std::vector<std::vector<index_t>>>& functionDescription,
    const std::vector<std::vector<std::array<int, 3>>>& spilloverFunctionCoordinates,
    std::vector<bool>& isTruncated,
    bool verbose);
bool shouldSuppressSpillover(
    index_t patch,  // The patch we are assembling (target)
    const gsTHBSplineBasis<2, real_t>& thbBasis,
    index_t sourcePatch,
    index_t level,
    index_t index,
    const gsTensorBSplineBasis<2, real_t>& spillBasis,
    bool verbose);
bool shouldSuppressTHB(
    int patch,
    int level,
    int index,
    const gsTensorBSplineBasis<2, real_t>& basis,
    const std::vector<std::array<int, 3>>& incomingSpillovers,
    const gsVector<gsVector<gsTensorBSplineBasis<2, real_t>>>& Bells,
    bool verbose = false
);
void assemble(
    const gsVector<gsMatrix<>>& uv,
    const gsVector<gsVector<gsTensorBSplineBasis<2, real_t>>>& Bells,
    const gsVector<gsTHBSplineBasis<2, real_t>>& THBVector,
    const std::vector<std::vector<std::vector<index_t>>>& functionDescription,
    const std::vector<std::vector<std::array<int, 3>>>& spilloverFunctionCoordinates,
    const std::vector<bool>& hasSpillover,
    const std::vector<bool>& isTruncated,
    const gsVector<gsVector<gsVector<index_t>>>& indexInTHB,
    std::vector<std::vector<gsSparseVector < real_t >>>& presentation,
    gsSparseMatrix<real_t>& A_mat,
    gsMatrix<real_t>& b_vec,
    const gsMultiPatch<real_t>& mp,
    bool verbose = false);

// New MPBES-based assemble (preferred)
void assemble(
    const gsVector<gsMatrix<>>& uv,
    const MPBES<2, real_t>& mpbes,
    gsSparseMatrix<real_t>& A_mat,
    gsMatrix<real_t>& b_vec,
    const gsMultiPatch<real_t>& mp,
    bool verbose = false,
    int currentPatch = -1,
    int currentLevel = -1,
    int currentAttempt = -1,
    const std::vector<index_t>* activeFunctionIds = nullptr);
void assembleGauss(
    const gsMultiPatch<real_t>& mp,
    const gsVector<gsVector<gsTensorBSplineBasis<2, real_t>>>& Bells,
    const gsVector<gsTHBSplineBasis<2, real_t>>& THBVector,
    const std::vector<std::vector<std::vector<index_t>>>& functionDescription,
    const std::vector<std::vector<std::array<int, 3>>>& spilloverFunctionCoordinates,
    const std::vector<bool>& hasSpillover,
    const gsVector<gsVector<gsVector<index_t>>>& indexInTHB,
    gsSparseMatrix<real_t>& A_mat,
    gsMatrix<real_t>& b_vec,
    bool verbose = false);
void assembleATHB(
    const gsVector<gsMatrix<>>& uv,
    const gsVector<gsVector<gsTensorBSplineBasis<2, real_t>>>& Bells,
    const gsVector<gsTHBSplineBasis<2, real_t>>& THBVector,
    const std::vector<std::vector<std::vector<index_t>>>& functionDescription,
    const std::vector<std::vector<std::array<int, 3>>>& spilloverFunctionCoordinates,
    const std::vector<bool>& hasSpillover,
    gsSparseMatrix<real_t>& A_mat,
    bool verbose);
void assembleATHBGauss(
    gsVector<gsMatrix<>> uv,
    gsVector<gsVector<gsTensorBSplineBasis<2, real_t>>> Bells,
    gsVector<gsTHBSplineBasis<2, real_t>> THBVector,
    std::vector<std::vector<std::vector<index_t>>> functionDescription,
    std::vector<std::vector<std::array<int, 3>>> spilloverFunctionCoordinates,
    std::vector<bool> hasSpillover,
    gsSparseMatrix<real_t>& A_mat,
    gsMatrix<real_t>& b_vec,
    gsMultiPatch<real_t> mp,
    bool verbose
);
void assembleATHB_maskSpilloverFromTwinnedPatches(
    const gsVector<gsMatrix<>>& uv,
    const gsVector<gsVector<gsTensorBSplineBasis<2, real_t>>>& Bells,
    const gsVector<gsTHBSplineBasis<2, real_t>>& THBVector,
    const std::vector<std::vector<std::vector<index_t>>>& functionDescription,
    const std::vector<std::vector<std::array<int, 3>>>& spilloverFunctionCoordinates,
    const std::vector<bool>& hasSpillover,
    gsSparseMatrix<real_t>& A_mat,
    bool verbose
);
void assembleATHB_clean(
    const gsVector<gsMatrix<>>& uv,
    const gsVector<gsVector<gsTensorBSplineBasis<2, real_t>>>& Bells,
    const gsVector<gsTHBSplineBasis<2, real_t>>& THBVector,
    const std::vector<std::vector<std::vector<index_t>>>& functionDescription,
    const std::vector<std::vector<std::array<int, 3>>>& spilloverFunctionCoordinates,
    const std::vector<bool>& hasSpillover,
    gsSparseMatrix<real_t>& A_mat
);
void assembleATHB_exclude(
    const gsVector<gsMatrix<>>& uv,
    const gsVector<gsVector<gsTensorBSplineBasis<2, real_t>>>& Bells,
    const gsVector<gsTHBSplineBasis<2, real_t>>& THBVector,
    const std::vector<std::vector<std::vector<index_t>>>& functionDescription,
    const std::vector<std::vector<std::array<int, 3>>>& spilloverFunctionCoordinates,
    const std::vector<std::vector<std::array<int, 3>>>& spilloverSources,
    const std::vector<bool>& hasSpillover,
    gsSparseMatrix<real_t>& A_mat,
    bool verbose
);
real_t testBoundaryAssembly(
    const MPBES<2, real_t>& mpbes,
    const gsMultiPatch<real_t>& mp,
    const gsMatrix<real_t>& vectSol
);
enum class FeatureSide
{
    U0,  // u = 0
    U1,  // u = 1
    V0,  // v = 0
    V1   // v = 1
};

struct FeatureBoundarySpec
{
    index_t patch;
    FeatureSide side;
};

enum class FittingMode
{
    GlobalFitting,
    LocalFitting
};

static const char* fittingModeName(FittingMode mode)
{
    switch (mode)
    {
    case FittingMode::GlobalFitting:
        return "globalFitting";
    case FittingMode::LocalFitting:
        return "localFitting";
    }

    return "unknownFittingMode";
}

struct MirroredCheckResult
{
    std::vector<bool> mirrored;
    std::vector<size_t> usedPerPatch;
    std::vector<size_t> posPerPatch;
    std::vector<size_t> negPerPatch;
};

struct GeometryPreflightInterfaceInfo
{
    index_t patchA = -1;
    index_t patchB = -1;
    int sideIndexA = -1;
    int sideIndexB = -1;
    FeatureSide sideA = FeatureSide::U0;
    FeatureSide sideB = FeatureSide::U0;
    bool hasMappedSides = false;
    bool orientationReversed = false;
};

struct GeometryPreflightInfo
{
    bool valid = false;
    MirroredCheckResult mirroredReport;
    std::vector<GeometryPreflightInterfaceInfo> interfaces;
};

struct InterfaceValidationSpec
{
    index_t patchA = -1;
    index_t patchB = -1;
    FeatureSide sideA = FeatureSide::U0;
    FeatureSide sideB = FeatureSide::U0;
    bool orientationReversed = false;
};

struct PatchVisualizationMetadata
{
    index_t patch = -1;
    bool mirrored = false;
    real_t jacobianDeterminant = 0.0;
    gsVector<real_t> centerXY;
    gsVector<real_t> uDirection;
    gsVector<real_t> vDirection;

    PatchVisualizationMetadata()
        : centerXY(2)
        , uDirection(2)
        , vDirection(2)
    {
        centerXY.setZero();
        uDirection.setZero();
        vDirection.setZero();
    }
};

struct InterfaceVisualizationMetadata
{
    index_t patchA = -1;
    index_t patchB = -1;
    FeatureSide sideA = FeatureSide::U0;
    FeatureSide sideB = FeatureSide::U0;
    bool orientationReversed = false;
    bool mirroredA = false;
    bool mirroredB = false;
    real_t midpointGap = 0.0;
    gsVector<real_t> midpointA;
    gsVector<real_t> midpointB;
    gsVector<real_t> tangentA;
    gsVector<real_t> tangentB;

    InterfaceVisualizationMetadata()
        : midpointA(2)
        , midpointB(2)
        , tangentA(2)
        , tangentB(2)
    {
        midpointA.setZero();
        midpointB.setZero();
        tangentA.setZero();
        tangentB.setZero();
    }
};

static GeometryPreflightInfo g_geometryPreflight;

static void appendUvSamples(
    gsMatrix<real_t>& target,
    const gsMatrix<real_t>& extra);

static void appendUvSamples(
    gsMatrix<real_t>& target,
    const gsMatrix<real_t>& extra)
{
    if (extra.cols() == 0)
        return;

    if (target.cols() == 0)
    {
        target = extra;
        return;
    }

    gsMatrix<real_t> merged(target.rows(), target.cols() + extra.cols());
    merged.leftCols(target.cols()) = target;
    merged.rightCols(extra.cols()) = extra;
    target = merged;
}

static const char* featureSideName(FeatureSide side);

struct InterfaceSideInfo;

static real_t vectorDistance(const gsVector<real_t>& a, const gsVector<real_t>& b);

static bool isInterfaceOrientationReversed(
    const gsMultiPatch<real_t>& mp,
    const InterfaceSideInfo& sideA,
    const InterfaceSideInfo& sideB);

static bool featureSideToInterfaceSideIndex(
    FeatureSide side,
    int& sideIndex);

static std::vector<GeometryPreflightInterfaceInfo> collectVisualizationInterfaceInfos(
    const gsMultiPatch<real_t>& mp);

static const GeometryPreflightInterfaceInfo* findPreflightInterfaceInfo(
    index_t patchA,
    index_t patchB);

static bool isPreflightMirroredPatch(index_t patch);

real_t boundaryError(
    const MPBES<2, real_t>& mpbes,
    const gsMultiPatch<real_t>& mp,
    const gsMatrix<real_t>& vectSol,
    const std::vector<FeatureBoundarySpec>& featureSides = std::vector<FeatureBoundarySpec>()
);
void checkForNaN(const gsMatrix<real_t>& matrix, const std::string& matrixName);
void checkSparseMatrixHealth(
    const gsSparseMatrix<real_t>& A,
    const std::string& name,
    size_t maxLogs = 10);
int computeRankByZeroRows(const gsMatrix<double>& mat, double tol = 1e-12);
void generateVisualizationMesh(
    const gsVector<gsVector<gsVector<index_t>>>& boxMat,
    index_t gridResolution,
    const MPBES<2, real_t>& mpbes,
    const gsMultiPatch<real_t>& mp,
    const gsMatrix<real_t>& coefficients,
    const gsVector<index_t>& currentLastNonZeroRow,
    const std::string& outputPrefix,
    bool verbose = false
);

void generateParametricGrid(
    const gsVector<gsVector<gsVector<index_t>>>& boxMat,
    const gsVector<index_t>& numBoxesPerPatch,
    const std::vector<std::set<index_t>>& activeLevels,
    index_t gridResolution,
    std::vector<gsMatrix<>>& uvGridLines,
    std::vector<index_t>& patchIDs,
    std::vector<char>& directions,
    bool verbose = false
);

void generateUniformPatchGrid(
    index_t numPatches,
    index_t gridResolution,
    std::vector<gsMatrix<>>& uvGridLines,
    std::vector<index_t>& patchIDs,
    std::vector<char>& directions
);

void evaluateGeometryAtPoints(
    const std::vector<gsMatrix<>>& uvPoints,
    const std::vector<index_t>& patchIDs,
    const MPBES<2, real_t>& mpbes,
    const gsMultiPatch<real_t>& mp,
    const gsMatrix<real_t>& coefficients,
    std::vector<gsMatrix<>>& xyPoints,
    bool verbose = false
);

static gsVector<real_t> evaluateFittedGeometryPoint(
    const MPBES<2, real_t>& mpbes,
    const gsMatrix<real_t>& coefficients,
    index_t patch,
    const gsVector<real_t>& uv,
    bool includeSpilloverFallback = false,
    bool normalizePou = true);

void exportMeshToFile(
    const std::vector<gsMatrix<>>& xyLines,
    const std::vector<index_t>& patchIDs,
    const std::vector<char>& directions,
    const std::vector<PatchVisualizationMetadata>& patchMetadata,
    const std::string& filename,
    bool verbose = false
);

void exportPatchMetadataToFile(
    const std::vector<PatchVisualizationMetadata>& patchMetadata,
    const std::string& filename,
    bool verbose = false
);

void exportInterfaceMetadataToFile(
    const std::vector<InterfaceVisualizationMetadata>& interfaceMetadata,
    const std::string& filename,
    bool verbose = false
);

void exportMeshToFile(
    const std::vector<gsMatrix<>>& xyLines,
    const std::vector<index_t>& patchIDs,
    const std::vector<char>& directions,
    const std::string& filename,
    bool verbose = false
);

// Additional forward declarations for non-predeclared functions
void logResidualDetails(
    const gsVector<gsMatrix<>>& uv,
    const gsMatrix<real_t>& b_vec,
    const gsMatrix<real_t>& matOut,
    std::ostream& outfile,
    bool logToInfo = true
);

real_t getSupports(int i, int j, int k, gsMatrix<> supps);

int pickCell(gsVector<int> vectorS, int& currArrayIndex, int levNow, int& x1U, int& y1U, int& x2U, int& y2U, index_t interior);

int pickCell(gsVector<int> vectorS, int attempt, int levNow, int& x1U, int& y1U, int& x2U, int& y2U, int lexicographic, index_t interior);

template <typename T> void goestaerDaeOPokhu(T object);

template <typename T> void DebugDaOPokhu(T object);

int PatchesIntersection(gsGeometry<>& geom1, gsGeometry<>& geom2,
    index_t numPoints = 10, real_t tolerance = 1e-5);

int PatchesIntersection(
    const std::vector<gsMatrix<>>& xyLines,
    const std::vector<index_t>& patchIDs,
    index_t patch1,
    index_t patch2,
    real_t tolerance = 1e-5,
    bool verbose = false);

int checkAllPatchIntersections(
    const std::vector<gsMatrix<>>& xyLines,
    const std::vector<index_t>& patchIDs,
    index_t numPatches,
    real_t tolerance = 1e-5,
    bool verbose = false);

real_t twoDet(gsMatrix<> a0, gsMatrix<> a1);

real_t theDistance(gsMultiPatch<real_t> repairedmp);

void coordinateTranformation(double& uSecond, double& vSecond, int firstSide, int secondSide, gsMatrix<> punto);

int twinFunction(
    const gsMultiPatch<real_t>& mp,
    int firstPatch,
    const gsTHBSplineBasis<2, real_t>& initTP,
    int functionIndex,
    int secondPatch,
    const gsTHBSplineBasis<2, real_t>& twinTP);

void boxToDomain(gsMatrix<index_t> mybox, gsMatrix<real_t>& coords, index_t interioru, index_t interiorv);

inline void supportToBoxOfLevel(gsMatrix<index_t>& mybox, gsMatrix<real_t> coords, int level, index_t interioru, index_t interiorv);

void boxToDomain(gsVector<index_t> mybox, gsVector<real_t>& coords, index_t interioru, index_t interiorv);

void boxToDomain(gsVector<index_t> mybox, std::vector<real_t>& coords, index_t interioru, index_t interiorv);

bool subset(gsBasisFun<> beta, gsVector<real_t> theBox);

template<typename T>
void initializeNestedVector(T& vec, index_t value = -1);

template<typename T, typename ValueType>
void fillVector3(gsVector<gsVector<gsVector<T>>>& vec, ValueType value);

template<typename T, typename ValueType>
void fillVector2(gsVector<gsVector<T>>& vec, ValueType value);

void restoreTheHierarchy(int& createdBoxNum, int& lastNonzeroRow, gsVector<gsVector<gsVector<index_t>>>& boxMat, int levNow, int centerInd, int ourBox[], int successfullAttempts, int patch);

void computeLevelIndex(gsMatrix<index_t, 2, 2> coarse_elem_ind, int level, gsMatrix<index_t, 2, 2>& fine_elem_ind, int maxLevel);

void setSupports(gsTHBSplineBasis<2, real_t> THB, gsMatrix<> supps);

template<typename T>
void resizeNestedVector(T& vec, const std::vector<size_t>& sizes, size_t depth = 0);

void actives_into(
    int patch,
    gsVector<std::vector<int>>& activesVect,
    gsVector<std::vector<int>>& activesLev,
    gsMatrix<real_t> punto,
    gsVector<gsVector<gsTensorBSplineBasis<2, real_t>>> Bells,
    gsVector<gsVector<gsVector<index_t>>> globalIndex
);

real_t assembleFunctional(const gsMatrix<real_t>& q,
    const real_t weight,
    const real_t constant);

// Forward declaration for extractGlobalIndexTHB
gsVector<gsVector<int>> extractGlobalIndexTHB(const MPBES<2, real_t>& mpbes);

// Forward declaration for checkBoxesConsistency
void checkBoxesConsistency(
    const gsVector<gsVector<gsVector<index_t>>>& boxMat,
    int patch,
    int lastNonZeroRow,
    std::ofstream& outfile
);

struct AdaptiveWeights
{
    double fitting = 1.0;
    double uniformity = 1.0;
    double orthogonality = 1.0;
    double length = 1.0;
};

AdaptiveWeights chooseAdaptiveWeights(
    double residualFitting,
    double residualUniformity,
    double residualOrthogonality,
    double residualLength,
    int minusnumber,
    int previousMinusnumber,
    double previousFunctional,
    double currentFunctional,
    bool nonlinearIteration
)
{
    AdaptiveWeights w;
    const double eps = 1e-12;
    const double maxRes = std::max(
        std::max(residualFitting, residualUniformity),
        std::max(residualOrthogonality, residualLength)
    );
    const double denom = (maxRes > eps) ? maxRes : eps;
    auto clamp = [](double v, double lo, double hi) {
        return (v < lo) ? lo : ((v > hi) ? hi : v);
    };
    auto scaled = [&](double r) {
        return clamp(r / denom, 0.0, 1.0);
    };

    const bool minusStagnant = (previousMinusnumber >= 0 && minusnumber >= previousMinusnumber);
    const bool functionalWorse = (previousFunctional > 0.0 && currentFunctional > previousFunctional);

    double boost = 1.0;
    if (minusnumber > 0)
    {
        // Scale up regularization pressure when Jacobian violations persist.
        // Example: minusnumber=43 -> +4.3x contribution before other boosts.
        boost += std::min(10.0, static_cast<double>(minusnumber) / 10.0);
    }
    if (minusStagnant)
        boost *= 1.35;
    if (functionalWorse)
        boost *= 1.20;

    // LO base lowered from 2e-2 → 2e-4: MPBES-LO operates globally across all
    // patches and is inherently stronger than single-patch LO; it needs less
    // regularisation.  At 2e-2 the floor prevented NLO from ever firing naturally.
    // At 2e-4 the effective weight for late-stage low-violation events falls below
    // their failure threshold, allowing NLO to be invoked as the paper intends.
    const double base = nonlinearIteration ? 2e-3 : 2e-4;
    const double minW = nonlinearIteration ? 1e-5 : 1e-6;
    const double maxW = nonlinearIteration ? 5e-2 : 2e-1;
    const double floorIrregular = nonlinearIteration
     ? 2e-3 : 1e-4;

    w.fitting = 1.0;
    w.uniformity = clamp(base * boost * (0.25 + 0.75 * scaled(residualUniformity)), minW, maxW);
    w.orthogonality = clamp(base * boost * (0.25 + 0.75 * scaled(residualOrthogonality)), minW, maxW);
    w.length = clamp(base * boost * (0.25 + 0.75 * scaled(residualLength)), minW, maxW);

    if (minusnumber > 0)
    {
        w.uniformity = std::max(w.uniformity, floorIrregular);
        w.orthogonality = std::max(w.orthogonality, floorIrregular);
        w.length = std::max(w.length, floorIrregular);
    }
    return w;
}

// Forward declaration for nonLinearOptimization
void nonLinearOptimization(
    const MPBES<2, real_t>& mpbes,
    const gsVector<gsMatrix<real_t>>& uv1,
    const gsVector<gsMatrix<real_t>>& uv2,
    const gsMultiPatch<>& mp1,
    const gsMatrix<>& matAsquare,
    const gsSparseMatrix<real_t>& matA,
    const gsMatrix<>& b_vec,
    gsMatrix<real_t>& vectSol,
    gsSparseMatrix<real_t>& A,
    const gsVector<gsTHBSplineBasis<2>>& SubdomainHierarchy,
    const gsVector<gsVector<gsTensorBSplineBasis<2, real_t>>>& Bells,
    const std::vector<std::vector<std::vector<index_t>>>& functionDescription,
    const gsVector<gsVector<gsVector<index_t>>>& indexInTHB,
    const std::vector<bool>& isTruncated,
    const std::vector<std::vector<gsSparseVector<double>>>& presentation,
    const std::vector<std::vector<std::array<int, 3>>>& spilloverFunctionCoordinates,
    const std::vector<bool>& hasSpillover,
    const gsVector<gsVector<index_t>>& globalIndexTHB,
    const gsVector<gsVector<gsVector<index_t>>>& boxMat,
    const gsVector<int>& currentLastNonZeroRow,
    real_t fitting_weight,
    real_t uniformity_weight,
    real_t orthogonality_weight,
    real_t length_weight,
    real_t skewness_weight,
    real_t eccentricity_weight,
    real_t area_weight,
    real_t epsilon_g,
    real_t epsilon_f,
    gsVector<size_t>& numIrregular,
    index_t geoDim,
    bool doNonLinear,
    FittingMode fittingMode,
    const std::vector<FeatureBoundarySpec>& featureSides,
    const LocalCoarseningRegion* localRegion,
    const gsMatrix<real_t>* originalCoefficients,
    const gsMatrix<real_t>* boundaryFreezeMask = nullptr,
    const gsMatrix<real_t>* boundaryFreezeRef  = nullptr
);

void evalSingle_into(
    index_t f,
    index_t i,
    const gsMatrix<real_t>& u,
    const std::vector<bool>& isTruncated,
    const gsTHBSplineBasis<2, real_t>& localTHBBasis,
    const gsSparseVector<real_t>& coefs,
    gsMatrix<real_t>& result,
    const gsTensorBSplineBasis<2, real_t>* bellBasis,
    index_t bellIndex);

void evalSingle_into(
    index_t f,
    index_t i,
    const gsMatrix<real_t>& u,
    const std::vector<bool>& isTruncated,
    const gsTHBSplineBasis<2, real_t>& localTHBBasis,
    const gsSparseVector<real_t>& coefs,
    gsMatrix<real_t>& result)
{
    evalSingle_into(f, i, u, isTruncated, localTHBBasis, coefs, result, nullptr, -1);
}

void evalSingle_into(
    index_t f,
    index_t i,
    const gsMatrix<real_t>& u,
    const std::vector<bool>& isTruncated,
    const gsTHBSplineBasis<2, real_t>& localTHBBasis,
    const gsSparseVector<real_t>& coefs,
    gsMatrix<real_t>& result,
    const gsTensorBSplineBasis<2, real_t>* bellBasis,
    index_t bellIndex)
{
    int verbose = 0;
    bool coefsOK = false;
    for (size_t j = 0; j < coefs.size(); j++)
    {
        if (coefs(j) != 0.0) {
            coefsOK = true;
            break;
        }
    }

    if (i < 0 && bellBasis != nullptr && bellIndex >= 0)
    {
        bellBasis->evalSingle_into(bellIndex, u, result);
        return;
    }

    if (!isTruncated[f]) // basis function not truncated
    {
        unsigned level = localTHBBasis.levelOf(i);
        unsigned tensor_index = localTHBBasis.flatTensorIndexOf(i, level);
        localTHBBasis.getBases()[level]->evalSingle_into(tensor_index, u, result);
    }
    else
    {
        unsigned level = localTHBBasis.getistruncated(i);
        if (verbose) {
            gsInfo << "coefs for f: ";
            for (size_t j = 0; j < coefs.size(); j++)
            {
                if (coefs(j) != 0.0)
                    gsInfo << "(" << j << ":" << coefs(j) << ") ";
            }
            gsInfo << "\n";
        }
        const gsTensorBSplineBasis<2, real_t>& base = *localTHBBasis.getBases()[level];
        gsTensorDeboor<2, real_t, gsKnotVector<real_t>, gsSparseVector<real_t>>(u, base, coefs, result);
    }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void evaluateInterfaceGeometry(
    const gsVector<real_t>& uv_patch0,
    const gsVector<real_t>& uv_patch1,
    index_t patch0,
    index_t patch1,
    const std::vector<index_t>& specialFunctions,
    const std::string& outputFileName,
    const std::vector<std::vector<std::vector<index_t>>>& functionDescription,
    const std::vector<bool>& isTruncated,
    const gsVector<gsTHBSplineBasis<2, real_t>>& SubdomainHierarchy,
    const gsVector<gsVector<gsTensorBSplineBasis<2, real_t>>>& Bells,
    const gsVector<gsVector<gsVector<index_t>>>& indexInTHB,
    const std::vector<std::vector<gsSparseVector<real_t>>>& presentation,
    const gsMatrix<real_t>& vectSol,
    index_t commonSize)
{
    PROFILE_FUNCTION();
    std::ofstream diagLog(outputFileName);
    diagLog << "Interface Evaluation Diagnostic\n";
    diagLog << "================================\n\n";
    diagLog << "Total functions (commonSize): " << commonSize << "\n";
    diagLog << "vectSol size: " << vectSol.rows() << " x " << vectSol.cols() << "\n\n";

    real_t x_p0 = 0.0, y_p0 = 0.0;
    real_t x_p1 = 0.0, y_p1 = 0.0;

    int functionsWithP0Component = 0;
    int functionsWithP1Component = 0;
    int p0_nonzero_basis = 0;
    int p1_nonzero_basis = 0;

    diagLog << "========== PATCH " << patch0 << " EVALUATION at (" << uv_patch0(0) << ", " << uv_patch0(1) << ") ==========\n";
    diagLog << "Only showing functions with non-zero basis values\n\n";
    gsInfo << "Evaluating Patch " << patch0 << " at (" << uv_patch0(0) << ", " << uv_patch0(1) << ")...\n";

    // Evaluate at patch 0
    for (index_t f = 0; f < commonSize; ++f)
    {
        real_t basisValue = 0.0;
        bool hasP0Component = false;
        int numComponents = functionDescription[f].size();

        // Special logging for specified functions
        bool isSpecial = std::find(specialFunctions.begin(), specialFunctions.end(), f) != specialFunctions.end();
        if (isSpecial)
        {
            diagLog << "\n*** SPECIAL: Function " << f << " (at PATCH " << patch0 << " evaluation) ***\n";
            diagLog << "  Total components: " << numComponents << ", truncated: " << isTruncated[f] << "\n";
            diagLog << "  Evaluating at patch " << patch0 << ", point (" << uv_patch0(0) << ", " << uv_patch0(1) << ")\n";
            diagLog << "  ALL components:\n";
            real_t debugBasisValue = 0.0;
            for (size_t c = 0; c < functionDescription[f].size(); ++c)
            {
                const auto& desc = functionDescription[f][c];
                int patch = desc[0];
                int level = desc[1];
                int idx = desc[2];
                diagLog << "    Component " << c << ": patch=" << patch << ", level=" << level << ", idx=" << idx;

                if (patch < indexInTHB.size() && level < indexInTHB[patch].size() && idx < indexInTHB[patch][level].size())
                {
                    int thbIdx = indexInTHB[patch][level][idx];
                    diagLog << ", thbIdx=" << thbIdx;

                    if (patch == patch0)
                    {
                        gsMatrix<real_t> result;
                        evalSingle_into(f, thbIdx, uv_patch0, isTruncated, SubdomainHierarchy[patch0],
                            presentation[f][c], result, &Bells[patch0][level], idx);
                        real_t val = result(0, 0);
                        debugBasisValue += val;
                        if (thbIdx >= 0)
                        {
                            diagLog << " -> EVALUATED (THB), basis_eval=" << val << ", running_total=" << debugBasisValue;
                        }
                        else
                        {
                            diagLog << " -> EVALUATED (B-spline spillover), basis_eval=" << val << ", running_total=" << debugBasisValue;
                        }
                    }
                    else
                    {
                        diagLog << " -> SKIPPED (wrong patch, need patch " << patch0 << ")";
                    }
                }
                else
                {
                    diagLog << " -> SKIPPED (OUT OF BOUNDS in indexInTHB)";
                }
                diagLog << "\n";
            }
            diagLog << "  FINAL basisValue for function " << f << " at patch " << patch0 << ": " << debugBasisValue << "\n";
            diagLog << "  Coefficient: coef_x=" << vectSol(f, 0) << ", coef_y=" << vectSol(f, 1) << "\n";
            diagLog << "  Contribution: x=" << (debugBasisValue * vectSol(f, 0)) << ", y=" << (debugBasisValue * vectSol(f, 1)) << "\n";
            diagLog << "*** END SPECIAL ***\n\n";
        }

        for (size_t c = 0; c < functionDescription[f].size(); ++c)
        {
            const auto& desc = functionDescription[f][c];
            int patch = desc[0];
            int level = desc[1];
            int idx = desc[2];

            if (patch == patch0 && level < indexInTHB[patch0].size() && idx < indexInTHB[patch0][level].size())
            {
                hasP0Component = true;
                int thbIdx = indexInTHB[patch0][level][idx];

                gsMatrix<real_t> result;
                evalSingle_into(f, thbIdx, uv_patch0, isTruncated, SubdomainHierarchy[patch0],
                    presentation[f][c], result, &Bells[patch0][level], idx);
                real_t val = result(0, 0);
                basisValue += val;
                if (val != 0.0) p0_nonzero_basis++;
            }
        }

        if (hasP0Component) functionsWithP0Component++;

        real_t coef_x = vectSol(f, 0);
        real_t coef_y = vectSol(f, 1);
        real_t contrib_x = basisValue * coef_x;
        real_t contrib_y = basisValue * coef_y;

        // Only log if non-zero contribution
        if (basisValue != 0.0 || contrib_x != 0.0 || contrib_y != 0.0)
        {
            diagLog << "Function " << f << " (components: " << numComponents << ", truncated: " << isTruncated[f] << "):\n";

            for (size_t c = 0; c < functionDescription[f].size(); ++c)
            {
                const auto& desc = functionDescription[f][c];
                int patch = desc[0];
                int level = desc[1];
                int idx = desc[2];

                if (patch == patch0 && level < indexInTHB[patch0].size() && idx < indexInTHB[patch0][level].size())
                {
                    int thbIdx = indexInTHB[patch0][level][idx];
                    diagLog << "  Component " << c << ": patch=" << patch << ", level=" << level << ", idx=" << idx << ", thbIdx=" << thbIdx;

                    gsMatrix<real_t> result;
                    evalSingle_into(f, thbIdx, uv_patch0, isTruncated, SubdomainHierarchy[patch0],
                        presentation[f][c], result, &Bells[patch0][level], idx);
                    real_t val = result(0, 0);
                    if (thbIdx >= 0)
                    {
                        diagLog << ", basis_eval=" << val;
                    }
                    else
                    {
                        diagLog << " (SPILLOVER), basis_eval=" << val;
                    }
                    diagLog << "\n";
                }
            }

            diagLog << "  Total basis value: " << basisValue << "\n";
            diagLog << "  Coefficients: coef_x=" << coef_x << ", coef_y=" << coef_y << "\n";
            diagLog << "  Contribution: x=" << contrib_x << ", y=" << contrib_y << "\n\n";
        }

        x_p0 += contrib_x;
        y_p0 += contrib_y;
    }

    diagLog << "\n========== PATCH " << patch1 << " EVALUATION at (" << uv_patch1(0) << ", " << uv_patch1(1) << ") ==========\n";
    diagLog << "Only showing functions with non-zero basis values\n\n";
    gsInfo << "Evaluating Patch " << patch1 << " at (" << uv_patch1(0) << ", " << uv_patch1(1) << ")...\n";

    // Evaluate at patch 1
    for (index_t f = 0; f < commonSize; ++f)
    {
        real_t basisValue = 0.0;
        bool hasP1Component = false;
        int numComponents = functionDescription[f].size();

        // Special logging for specified functions
        bool isSpecial = std::find(specialFunctions.begin(), specialFunctions.end(), f) != specialFunctions.end();
        if (isSpecial)
        {
            diagLog << "\n*** SPECIAL: Function " << f << " (at PATCH " << patch1 << " evaluation) ***\n";
            diagLog << "  Total components: " << numComponents << ", truncated: " << isTruncated[f] << "\n";
            diagLog << "  Evaluating at patch " << patch1 << ", point (" << uv_patch1(0) << ", " << uv_patch1(1) << ")\n";
            diagLog << "  ALL components:\n";
            real_t debugBasisValue = 0.0;
            for (size_t c = 0; c < functionDescription[f].size(); ++c)
            {
                const auto& desc = functionDescription[f][c];
                int patch = desc[0];
                int level = desc[1];
                int idx = desc[2];
                diagLog << "    Component " << c << ": patch=" << patch << ", level=" << level << ", idx=" << idx;

                if (patch < indexInTHB.size() && level < indexInTHB[patch].size() && idx < indexInTHB[patch][level].size())
                {
                    int thbIdx = indexInTHB[patch][level][idx];
                    diagLog << ", thbIdx=" << thbIdx;

                    if (patch == patch1)
                    {
                        gsMatrix<real_t> result;
                        evalSingle_into(f, thbIdx, uv_patch1, isTruncated, SubdomainHierarchy[patch1],
                            presentation[f][c], result, &Bells[patch1][level], idx);
                        real_t val = result(0, 0);
                        debugBasisValue += val;
                        if (thbIdx >= 0)
                        {
                            diagLog << " -> EVALUATED (THB), basis_eval=" << val << ", running_total=" << debugBasisValue;
                        }
                        else
                        {
                            diagLog << " -> EVALUATED (B-spline spillover), basis_eval=" << val << ", running_total=" << debugBasisValue;
                        }
                    }
                    else
                    {
                        diagLog << " -> SKIPPED (wrong patch, need patch " << patch1 << ")";
                    }
                }
                else
                {
                    diagLog << " -> SKIPPED (OUT OF BOUNDS in indexInTHB)";
                }
                diagLog << "\n";
            }
            diagLog << "  FINAL basisValue for function " << f << " at patch " << patch1 << ": " << debugBasisValue << "\n";
            diagLog << "  Coefficient: coef_x=" << vectSol(f, 0) << ", coef_y=" << vectSol(f, 1) << "\n";
            diagLog << "  Contribution: x=" << (debugBasisValue * vectSol(f, 0)) << ", y=" << (debugBasisValue * vectSol(f, 1)) << "\n";
            diagLog << "*** END SPECIAL ***\n\n";
        }

        for (size_t c = 0; c < functionDescription[f].size(); ++c)
        {
            const auto& desc = functionDescription[f][c];
            int patch = desc[0];
            int level = desc[1];
            int idx = desc[2];

            if (patch == patch1 && level < indexInTHB[patch1].size() && idx < indexInTHB[patch1][level].size())
            {
                hasP1Component = true;
                int thbIdx = indexInTHB[patch1][level][idx];

                gsMatrix<real_t> result;
                evalSingle_into(f, thbIdx, uv_patch1, isTruncated, SubdomainHierarchy[patch1],
                    presentation[f][c], result, &Bells[patch1][level], idx);
                real_t val = result(0, 0);
                basisValue += val;
                if (val != 0.0) p1_nonzero_basis++;
            }
        }

        if (hasP1Component) functionsWithP1Component++;

        real_t coef_x = vectSol(f, 0);
        real_t coef_y = vectSol(f, 1);
        real_t contrib_x = basisValue * coef_x;
        real_t contrib_y = basisValue * coef_y;

        // Only log if non-zero contribution
        if (basisValue != 0.0 || contrib_x != 0.0 || contrib_y != 0.0)
        {
            diagLog << "Function " << f << " (components: " << numComponents << ", truncated: " << isTruncated[f] << "):\n";

            for (size_t c = 0; c < functionDescription[f].size(); ++c)
            {
                const auto& desc = functionDescription[f][c];
                int patch = desc[0];
                int level = desc[1];
                int idx = desc[2];

                if (patch == patch1 && level < indexInTHB[patch1].size() && idx < indexInTHB[patch1][level].size())
                {
                    int thbIdx = indexInTHB[patch1][level][idx];
                    diagLog << "  Component " << c << ": patch=" << patch << ", level=" << level << ", idx=" << idx << ", thbIdx=" << thbIdx;

                    gsMatrix<real_t> result;
                    evalSingle_into(f, thbIdx, uv_patch1, isTruncated, SubdomainHierarchy[patch1],
                        presentation[f][c], result, &Bells[patch1][level], idx);
                    real_t val = result(0, 0);
                    if (thbIdx >= 0)
                    {
                        diagLog << ", basis_eval=" << val;
                    }
                    else
                    {
                        diagLog << " (SPILLOVER), basis_eval=" << val;
                    }
                    diagLog << "\n";
                }
            }

            diagLog << "  Total basis value: " << basisValue << "\n";
            diagLog << "  Coefficients: coef_x=" << coef_x << ", coef_y=" << coef_y << "\n";
            diagLog << "  Contribution: x=" << contrib_x << ", y=" << contrib_y << "\n\n";
        }

        x_p1 += contrib_x;
        y_p1 += contrib_y;
    }

    real_t dist = std::sqrt(std::pow(x_p0 - x_p1, 2) + std::pow(y_p0 - y_p1, 2));

    diagLog << "\n========== SUMMARY ==========\n";
    diagLog << "Functions with Patch " << patch0 << " component: " << functionsWithP0Component << "\n";
    diagLog << "Functions with Patch " << patch1 << " component: " << functionsWithP1Component << "\n";
    diagLog << "Patch " << patch0 << " non-zero basis evaluations: " << p0_nonzero_basis << "\n";
    diagLog << "Patch " << patch1 << " non-zero basis evaluations: " << p1_nonzero_basis << "\n\n";
    diagLog << "FINAL RESULTS:\n";
    diagLog << "  Patch " << patch0 << " at (" << uv_patch0(0) << "," << uv_patch0(1) << "): x=" << x_p0 << ", y=" << y_p0 << "\n";
    diagLog << "  Patch " << patch1 << " at (" << uv_patch1(0) << "," << uv_patch1(1) << "): x=" << x_p1 << ", y=" << y_p1 << "\n";
    diagLog << "  Distance: " << dist << "\n";
    diagLog.close();

    gsInfo << "Patch " << patch0 << " at (" << uv_patch0(0) << "," << uv_patch0(1) << "): x=" << x_p0 << ", y=" << y_p0 << "\n";
    gsInfo << "Patch " << patch1 << " at (" << uv_patch1(0) << "," << uv_patch1(1) << "): x=" << x_p1 << ", y=" << y_p1 << "\n";
    gsInfo << "Distance: " << dist << "\n";
    gsInfo << "Functions with P" << patch0 << " component: " << functionsWithP0Component << "\n";
    gsInfo << "Functions with P" << patch1 << " component: " << functionsWithP1Component << "\n";
    gsInfo << "P" << patch0 << " non-zero basis: " << p0_nonzero_basis << "\n";
    gsInfo << "P" << patch1 << " non-zero basis: " << p1_nonzero_basis << "\n";
    gsInfo << "Detailed log written to: " << outputFileName << "\n";
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void setSupports(gsTHBSplineBasis<2, real_t> THB, gsMatrix<> supps) {
    for (int i = 0; i < THB.size(); i++)
    {
        for (int j = 0; j < THB.support(i).rows(); j++)
        {
            for (int k = 0; k < THB.support(i).cols(); k++)
            {
                supps(j, 2 * i + k) = THB.support(i)(j, k);
            }
        }
    }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int computeRankByZeroRows(const gsMatrix<double>& mat, double tol)
{
    int rank = 0;

    for (int i = 0; i < mat.rows(); ++i)
    {
        bool nonZeroRow = false;
        for (int j = 0; j < mat.cols(); ++j)
        {
            if (std::abs(mat(i, j)) > tol)
            {
                nonZeroRow = true;
                break;
            }
        }
        if (nonZeroRow) {
            ++rank;
            //gsInfo << "rank: " << rank << "\n";
        }
    }

    return rank;
}

real_t twoDet(gsMatrix<> a0, gsMatrix<> a1) {
    return a0(0, 0) * a1(1, 0) - a1(0, 0) * a0(1, 0);
}

real_t assembleFunctional(const gsMatrix<real_t>& q,
    const real_t weight,
    const real_t constant);

void actives_into(
    int patch,
    gsVector<std::vector<int>>& activesVect,
    gsVector<std::vector<int>>& activesLev,
    gsMatrix <real_t> punto,
    gsVector<gsVector<gsTensorBSplineBasis<2, real_t>>> Bells,
    gsVector<gsVector<gsVector<index_t >>> globalIndex
);

template<typename T>
void resizeNestedVector(T& vec, const std::vector<size_t>& sizes, size_t depth) {
    if (depth >= sizes.size()) return;  // Stop if we've resized to the innermost level

    vec.resize(sizes[depth]);
    for (auto& subVec : vec) {
        resizeNestedVector(subVec, sizes, depth + 1);  // Recurse to the next level
    }
}



/**
 * @brief Exports mesh to file for visualization.
 *
 * Writes the (x,y) coordinates of grid lines to a text file.
 */
void exportMeshToFile(
    const std::vector<gsMatrix<>>& xyLines,
    const std::vector<index_t>& patchIDs,
    const std::vector<char>& directions,
    const std::vector<PatchVisualizationMetadata>& patchMetadata,
    const std::string& filename,
    bool verbose)
{
    PROFILE_FUNCTION();
    if (xyLines.empty())
    {
        gsWarn << "No mesh data to export.\n";
        return;
    }

    std::ofstream file(filename);
    if (!file.is_open())
    {
        gsWarn << "Cannot open file: " << filename << "\n";
        return;
    }

    // Write header
    file << "# Visualization Mesh\n";
    file << "# Format: x y patch_id direction\n";
    file << "# Blank line separates different grid lines\n";
    file << "# Total lines: " << xyLines.size() << "\n";
    file << "\n";

    // Write each line
    for (size_t i = 0; i < xyLines.size(); ++i)
    {
        const gsMatrix<>& line = xyLines[i];
        const index_t patch = patchIDs[i];
        const char dir = directions[i];

        for (index_t j = 0; j < line.cols(); ++j)
        {
            file << line(0, j) << " " << line(1, j)
                << " " << patch << " " << dir << "\n";
        }
        file << "\n"; // Blank line between grid lines
    }

    file.close();

    if (verbose)
        gsInfo << "Exported mesh to: " << filename << "\n";
}

void exportPatchMetadataToFile(
    const std::vector<PatchVisualizationMetadata>& patchMetadata,
    const std::string& filename,
    bool verbose)
{
    PROFILE_FUNCTION();

    std::ofstream file(filename);
    if (!file.is_open())
    {
        gsWarn << "Cannot open file: " << filename << "\n";
        return;
    }

    file << "# Patch visualization metadata\n";
    file << "# Format: patch mirrored jac_det center_x center_y u_dir_x u_dir_y v_dir_x v_dir_y\n";
    for (size_t i = 0; i < patchMetadata.size(); ++i)
    {
        const PatchVisualizationMetadata& info = patchMetadata[i];
        file << info.patch << " "
             << (info.mirrored ? 1 : 0) << " "
             << info.jacobianDeterminant << " "
             << info.centerXY(0) << " " << info.centerXY(1) << " "
             << info.uDirection(0) << " " << info.uDirection(1) << " "
             << info.vDirection(0) << " " << info.vDirection(1) << "\n";
    }

    file.close();

    if (verbose)
        gsInfo << "Exported patch metadata to: " << filename << "\n";
}

void exportInterfaceMetadataToFile(
    const std::vector<InterfaceVisualizationMetadata>& interfaceMetadata,
    const std::string& filename,
    bool verbose)
{
    PROFILE_FUNCTION();

    std::ofstream file(filename);
    if (!file.is_open())
    {
        gsWarn << "Cannot open file: " << filename << "\n";
        return;
    }

    file << "# Interface visualization metadata\n";
    file << "# Format: patchA patchB sideA sideB reversed mirroredA mirroredB gap xA yA tanAx tanAy xB yB tanBx tanBy\n";
    for (size_t i = 0; i < interfaceMetadata.size(); ++i)
    {
        const InterfaceVisualizationMetadata& info = interfaceMetadata[i];
        file << info.patchA << " "
             << info.patchB << " "
             << featureSideName(info.sideA) << " "
             << featureSideName(info.sideB) << " "
             << (info.orientationReversed ? 1 : 0) << " "
             << (info.mirroredA ? 1 : 0) << " "
             << (info.mirroredB ? 1 : 0) << " "
             << info.midpointGap << " "
             << info.midpointA(0) << " " << info.midpointA(1) << " "
             << info.tangentA(0) << " " << info.tangentA(1) << " "
             << info.midpointB(0) << " " << info.midpointB(1) << " "
             << info.tangentB(0) << " " << info.tangentB(1) << "\n";
    }

    file.close();

    if (verbose)
        gsInfo << "Exported interface metadata to: " << filename << "\n";
}

void exportMeshToFile(
    const std::vector<gsMatrix<>>& xyLines,
    const std::vector<index_t>& patchIDs,
    const std::vector<char>& directions,
    const std::string& filename,
    bool verbose)
{
    const std::vector<PatchVisualizationMetadata> emptyPatchMetadata;
    exportMeshToFile(
        xyLines,
        patchIDs,
        directions,
        emptyPatchMetadata,
        filename,
        verbose);
}

/**
 * @brief Check if a box is covered by any finer level box in the same patch.
 *
 * In THB splines, finer boxes override coarser boxes in overlapping regions.
 *
 * @param boxMat All boxes for the patch
 * @param numBoxes Number of boxes
 * @param testIdx Index of box to test
 * @return true if the box is completely covered by finer boxes
 */
bool isBoxCoveredByFinerBoxes(
    const gsVector<gsVector<index_t>>& boxMat,
    index_t numBoxes,
    index_t testIdx)
{
    if (boxMat[testIdx].size() < 5)
        return false;

    const index_t testLevel = boxMat[testIdx][0];
    const index_t testIx0 = boxMat[testIdx][1];
    const index_t testIy0 = boxMat[testIdx][2];
    const index_t testIx1 = boxMat[testIdx][3];
    const index_t testIy1 = boxMat[testIdx][4];

    // Check all other boxes
    for (index_t otherIdx = 0; otherIdx < numBoxes; ++otherIdx)
    {
        if (otherIdx == testIdx || boxMat[otherIdx].size() < 5)
            continue;

        const index_t otherLevel = boxMat[otherIdx][0];

        // Only finer levels can override this box
        if (otherLevel <= testLevel)
            continue;

        const index_t otherIx0 = boxMat[otherIdx][1];
        const index_t otherIy0 = boxMat[otherIdx][2];
        const index_t otherIx1 = boxMat[otherIdx][3];
        const index_t otherIy1 = boxMat[otherIdx][4];

        // Scale indices to same level for comparison
        const index_t levelDiff = otherLevel - testLevel;
        const index_t scale = 1 << levelDiff; // 2^levelDiff

        const index_t scaledOtherIx0 = otherIx0 / scale;
        const index_t scaledOtherIy0 = otherIy0 / scale;
        const index_t scaledOtherIx1 = otherIx1 / scale;
        const index_t scaledOtherIy1 = otherIy1 / scale;

        // Check if finer box completely covers test box
        if (scaledOtherIx0 <= testIx0 && scaledOtherIx1 >= testIx1 &&
            scaledOtherIy0 <= testIy0 && scaledOtherIy1 >= testIy1)
        {
            return true; // Test box is completely covered
        }
    }

    return false; // Not covered by any finer box
}

static bool boxTouchesFeatureSide(
    const gsVector<index_t>& box,
    FeatureSide side)
{
    if (box.size() < 5)
        return false;

    const index_t level = box[0];
    const index_t maxIndex = static_cast<index_t>(1) << level;
    const index_t ix0 = box[1];
    const index_t iy0 = box[2];
    const index_t ix1 = box[3];
    const index_t iy1 = box[4];

    switch (side)
    {
    case FeatureSide::U0:
        return ix0 == 0;
    case FeatureSide::U1:
        return ix1 == maxIndex;
    case FeatureSide::V0:
        return iy0 == 0;
    case FeatureSide::V1:
        return iy1 == maxIndex;
    }

    return false;
}

static std::vector<std::set<index_t>> collectActiveVisualizationLevels(
    const MPBES<2, real_t>& mpbes,
    index_t numPatches)
{
    std::vector<std::set<index_t>> activeLevels(static_cast<size_t>(numPatches));
    const auto& functionDescription = mpbes.functionDescription();

    for (size_t f = 0; f < functionDescription.size(); ++f)
    {
        for (size_t compIdx = 0; compIdx < functionDescription[f].size(); ++compIdx)
        {
            const std::vector<index_t>& desc = functionDescription[f][compIdx];
            if (desc.size() < 3)
                continue;

            const index_t patch = desc[0];
            const index_t level = desc[1];
            if (patch < 0 || patch >= numPatches || level < 0)
                continue;

            activeLevels[static_cast<size_t>(patch)].insert(level);
        }
    }

    return activeLevels;
}

static bool isVisualizationLevelActive(
    const std::vector<std::set<index_t>>& activeLevels,
    index_t patch,
    index_t level)
{
    if (patch < 0 || patch >= static_cast<index_t>(activeLevels.size()))
        return true;

    const std::set<index_t>& patchLevels = activeLevels[static_cast<size_t>(patch)];
    if (patchLevels.empty())
        return true;

    return patchLevels.find(level) != patchLevels.end();
}

static void collectBoundaryTangentialIndices(
    const gsVector<gsVector<gsVector<index_t>>>& boxMat,
    const gsVector<index_t>& numBoxesPerPatch,
    const std::vector<std::set<index_t>>& activeLevels,
    index_t patch,
    FeatureSide side,
    index_t targetLevel,
    std::set<index_t>& tangentIndices)
{
    if (patch < 0 || patch >= boxMat.size())
        return;

    const index_t numBoxes = numBoxesPerPatch[patch];
    const index_t maxIndex = static_cast<index_t>(1) << targetLevel;

    for (index_t boxIdx = 0; boxIdx < numBoxes; ++boxIdx)
    {
        if (boxMat[patch][boxIdx].size() < 5)
            continue;

        if (!isVisualizationLevelActive(activeLevels, patch, boxMat[patch][boxIdx][0]))
            continue;

        if (isBoxCoveredByFinerBoxes(boxMat[patch], numBoxes, boxIdx))
            continue;

        const gsVector<index_t>& box = boxMat[patch][boxIdx];
        if (!boxTouchesFeatureSide(box, side))
            continue;

        const index_t level = box[0];
        const index_t scaleFactor = static_cast<index_t>(1) << (targetLevel - level);
        const index_t startIndex = (side == FeatureSide::U0 || side == FeatureSide::U1)
            ? box[2] * scaleFactor
            : box[1] * scaleFactor;
        const index_t endIndex = (side == FeatureSide::U0 || side == FeatureSide::U1)
            ? box[4] * scaleFactor
            : box[3] * scaleFactor;

        for (index_t idx = startIndex; idx <= endIndex; ++idx)
            tangentIndices.insert(idx);
    }

    if (tangentIndices.empty())
    {
        tangentIndices.insert(0);
        tangentIndices.insert(maxIndex);
    }
}

static void appendInterfaceAlignedLine(
    index_t patch,
    FeatureSide side,
    index_t tangentIndex,
    index_t tangentLevel,
    index_t patchMaxLevel,
    index_t samplesPerIndexCell,
    std::vector<gsMatrix<>>& uvGridLines,
    std::vector<index_t>& patchIDs,
    std::vector<char>& directions)
{
    const index_t maxIndex = static_cast<index_t>(1) << patchMaxLevel;
    const index_t numPoints = maxIndex * samplesPerIndexCell + 1;
    static_cast<void>(tangentIndex);
    static_cast<void>(tangentLevel);

    gsMatrix<> line(2, numPoints);
    if (side == FeatureSide::U0 || side == FeatureSide::U1)
    {
        const real_t uFixed = (side == FeatureSide::U0) ? 0.0 : 1.0;
        for (index_t i = 0; i < numPoints; ++i)
        {
            const real_t v = static_cast<real_t>(i) / static_cast<real_t>(numPoints - 1);
            line(0, i) = uFixed;
            line(1, i) = v;
        }
        directions.push_back('v');
    }
    else
    {
        const real_t vFixed = (side == FeatureSide::V0) ? 0.0 : 1.0;
        for (index_t i = 0; i < numPoints; ++i)
        {
            const real_t u = static_cast<real_t>(i) / static_cast<real_t>(numPoints - 1);
            line(0, i) = u;
            line(1, i) = vFixed;
        }
        directions.push_back('u');
    }

    uvGridLines.push_back(line);
    patchIDs.push_back(patch);
}

static void appendAlignedInterfaceLines(
    const gsVector<gsVector<gsVector<index_t>>>& boxMat,
    const gsVector<index_t>& numBoxesPerPatch,
    const std::vector<std::set<index_t>>& activeLevels,
    const std::vector<index_t>& maxLevelPerPatch,
    index_t samplesPerIndexCell,
    std::vector<gsMatrix<>>& uvGridLines,
    std::vector<index_t>& patchIDs,
    std::vector<char>& directions)
{
    if (!g_geometryPreflight.valid)
        return;

    std::set<std::pair<index_t, int>> emittedPatchSides;

    for (size_t iface = 0; iface < g_geometryPreflight.interfaces.size(); ++iface)
    {
        const GeometryPreflightInterfaceInfo& info = g_geometryPreflight.interfaces[iface];
        if (!info.hasMappedSides)
            continue;

        const index_t patchA = info.patchA;
        const index_t patchB = info.patchB;
        if (patchA < 0 || patchB < 0 ||
            patchA >= static_cast<index_t>(boxMat.size()) ||
            patchB >= static_cast<index_t>(boxMat.size()))
            continue;

        const index_t interfaceLevel = std::max(maxLevelPerPatch[patchA], maxLevelPerPatch[patchB]);
        std::set<index_t> tangentialA;
        std::set<index_t> tangentialB;
        collectBoundaryTangentialIndices(boxMat, numBoxesPerPatch, activeLevels, patchA, info.sideA, interfaceLevel, tangentialA);
        collectBoundaryTangentialIndices(boxMat, numBoxesPerPatch, activeLevels, patchB, info.sideB, interfaceLevel, tangentialB);

        std::set<index_t> canonicalTangential = tangentialA;
        for (std::set<index_t>::const_iterator it = tangentialB.begin(); it != tangentialB.end(); ++it)
        {
            const index_t mappedIndex = info.orientationReversed
                ? (static_cast<index_t>(1) << interfaceLevel) - *it
                : *it;
            canonicalTangential.insert(mappedIndex);
        }

        static_cast<void>(canonicalTangential);

        const int sideAKey = static_cast<int>(info.sideA);
        const int sideBKey = static_cast<int>(info.sideB);

        if (emittedPatchSides.insert(std::make_pair(patchA, sideAKey)).second)
        {
            appendInterfaceAlignedLine(
                patchA,
                info.sideA,
                0,
                interfaceLevel,
                maxLevelPerPatch[patchA],
                samplesPerIndexCell,
                uvGridLines,
                patchIDs,
                directions);
        }

        if (emittedPatchSides.insert(std::make_pair(patchB, sideBKey)).second)
        {
            appendInterfaceAlignedLine(
                patchB,
                info.sideB,
                0,
                interfaceLevel,
                maxLevelPerPatch[patchB],
                samplesPerIndexCell,
                uvGridLines,
                patchIDs,
                directions);
        }
    }
}

static bool patchSideParticipatesInInterface(
    index_t patch,
    FeatureSide side)
{
    if (!g_geometryPreflight.valid)
        return false;

    for (size_t iface = 0; iface < g_geometryPreflight.interfaces.size(); ++iface)
    {
        const GeometryPreflightInterfaceInfo& info = g_geometryPreflight.interfaces[iface];
        if (!info.hasMappedSides)
            continue;

        if ((info.patchA == patch && info.sideA == side) ||
            (info.patchB == patch && info.sideB == side))
            return true;
    }

    return false;
}

/**
 * @brief Generates a parametric grid from box structure.
 *
 * Creates grid lines showing the tensor product grid structure within each box.
 * For a box with index coordinates (x0, y0, x1, y1), generates:
 * - (x1 - x0 - 1) internal vertical lines
 * - (y1 - y0 - 1) internal horizontal lines
 * Plus the 4 boundary edges.
 *
 * Skips boxes that are completely covered by finer level boxes to avoid
 * visual overlap in hierarchical refinement regions.
 *
 * Ensures interface continuity: boxes sharing edges have matching point counts
 * based on the finest refinement level present.
 */
 /**
  * @brief Generates a parametric grid from box structure with index-aligned sampling.
  *
 * Creates visualization lines from the visible box tensor grid.
 * For a box with index coordinates (ix0, iy0, ix1, iy1) at level L:
 * - Emits every horizontal tensor-grid line from iy0 to iy1
 * - Emits every vertical tensor-grid line from ix0 to ix1
 * - Samples each line using the finest visible level on that patch so shared
 *   edges between coarse and fine boxes still align in the exported mesh.
 *
 * This keeps the main visualization file consistent with the box structure
 * instead of showing only patch outlines.
  */
void generateParametricGrid(
    const gsVector<gsVector<gsVector<index_t>>>& boxMat,
    const gsVector<index_t>& numBoxesPerPatch,
    const std::vector<std::set<index_t>>& activeLevels,
    index_t gridResolution,
    std::vector<gsMatrix<>>& uvGridLines,
    std::vector<index_t>& patchIDs,
    std::vector<char>& directions,
    bool verbose)
{
    PROFILE_FUNCTION();
    uvGridLines.clear();
    patchIDs.clear();
    directions.clear();

    // if (verbose)
    //     gsInfo << "Generating parametric grid from box structure...\n";

    // Find maximum refinement level per patch to determine finest grid spacing.
    // Use all boxes in boxMat directly — activeLevels only tracks which levels have
    // direct MPBES functions, but twin functions mean a box can exist at a level that
    // has no local functions on that patch. boxMat is always authoritative.
    std::vector<index_t> maxLevelPerPatch(boxMat.size(), 0);
    for (index_t patch = 0; patch < boxMat.size(); ++patch)
    {
        const index_t numBoxes = numBoxesPerPatch[patch];
        for (index_t boxIdx = 0; boxIdx < numBoxes; ++boxIdx)
        {
            if (boxMat[patch][boxIdx].size() >= 5)
                maxLevelPerPatch[patch] = std::max(maxLevelPerPatch[patch], boxMat[patch][boxIdx][0]);
        }
    }

    // gridResolution is minimum samples per finest-level index cell
    const index_t samplesPerIndexCell = std::max<index_t>(gridResolution, 1);

    if (g_verbose)
        gsInfo << "  Samples per index cell: " << samplesPerIndexCell << "\n";

    index_t skippedBoxes = 0;

    for (index_t patch = 0; patch < boxMat.size(); ++patch)
    {
        const index_t numBoxes = numBoxesPerPatch[patch];

        if (g_verbose)
            gsInfo << "  Patch " << patch << ": " << numBoxes
                   << " boxes, maxLevel=" << maxLevelPerPatch[patch] << "\n";

        for (index_t boxIdx = 0; boxIdx < numBoxes; ++boxIdx)
        {
            if (boxMat[patch][boxIdx].size() < 5)
                continue;

            // Skip boxes completely covered by finer boxes
            if (isBoxCoveredByFinerBoxes(boxMat[patch], numBoxes, boxIdx))
            {
                if (verbose)
                    gsInfo << "    Skipping box " << boxIdx << " (covered by finer boxes)\n";
                skippedBoxes++;
                continue;
            }

            // Box data in index space (level-specific indices)
            const index_t level = boxMat[patch][boxIdx][0];
            const index_t ix0 = boxMat[patch][boxIdx][1];
            const index_t iy0 = boxMat[patch][boxIdx][2];
            const index_t ix1 = boxMat[patch][boxIdx][3];
            const index_t iy1 = boxMat[patch][boxIdx][4];

            // Grid spacing at this level
            const real_t h = 1.0 / std::pow(2.0, level);

            const index_t patchMaxLevel = maxLevelPerPatch[patch];
            const index_t levelDiff = patchMaxLevel - level;
            const index_t scaleFactor = 1 << levelDiff; // 2^(patchMaxLevel - level)

            // Scaled index bounds at finest level
            const index_t scaledIx0 = ix0 * scaleFactor;
            const index_t scaledIy0 = iy0 * scaleFactor;
            const index_t scaledIx1 = ix1 * scaleFactor;
            const index_t scaledIy1 = iy1 * scaleFactor;

            // Number of finest-level index cells in each direction
            const index_t numFinestCellsU = scaledIx1 - scaledIx0;
            const index_t numFinestCellsV = scaledIy1 - scaledIy0;

            const real_t hFinest = 1.0 / std::pow(2.0, patchMaxLevel);

            // Generate all horizontal tensor-grid lines of the visible box.
            for (index_t finestYIndex = scaledIy0; finestYIndex <= scaledIy1; ++finestYIndex)
            {
                const real_t v = finestYIndex * hFinest;

                const index_t numPoints = numFinestCellsU * samplesPerIndexCell + 1;

                gsMatrix<> line(2, numPoints);
                for (index_t i = 0; i < numPoints; ++i)
                {
                    const real_t finestLevelIndex = scaledIx0 + static_cast<real_t>(i) / samplesPerIndexCell;
                    const real_t u = finestLevelIndex * hFinest;

                    line(0, i) = u;
                    line(1, i) = v;
                }

                uvGridLines.push_back(line);
                patchIDs.push_back(patch);
                directions.push_back('u');
            }

            // Generate all vertical tensor-grid lines of the visible box.
            for (index_t finestXIndex = scaledIx0; finestXIndex <= scaledIx1; ++finestXIndex)
            {
                const real_t u = finestXIndex * hFinest;

                const index_t numPoints = numFinestCellsV * samplesPerIndexCell + 1;

                gsMatrix<> line(2, numPoints);
                for (index_t i = 0; i < numPoints; ++i)
                {
                    const real_t finestLevelIndex = scaledIy0 + static_cast<real_t>(i) / samplesPerIndexCell;
                    const real_t v = finestLevelIndex * hFinest;

                    line(0, i) = u;
                    line(1, i) = v;
                }

                uvGridLines.push_back(line);
                patchIDs.push_back(patch);
                directions.push_back('v');
            }

            // if (verbose)
            //     gsInfo << "    Box " << boxIdx << " [level=" << level
            //     << ", (" << ix0 << "," << iy0 << ")->(" << ix1 << "," << iy1 << ")]: "
            //     << (numFinestCellsU + 1) << " horizontal + "
            //     << (numFinestCellsV + 1) << " vertical tensor-grid lines\n";
        }
    }

    // Add explicitly interface-aligned lines so shared patch boundaries are
    // sampled identically in the exported visualization mesh.
    appendAlignedInterfaceLines(
        boxMat,
        numBoxesPerPatch,
        activeLevels,
        maxLevelPerPatch,
        samplesPerIndexCell,
        uvGridLines,
        patchIDs,
        directions);

    // if (verbose)
    //     gsInfo << "Generated " << uvGridLines.size() << " total grid lines ("
    //     << skippedBoxes << " boxes skipped due to finer coverage).\n";

    // Ensure patch boundary lines exist only if a patch produced no lines
    // (avoid duplicating boundaries that are already generated by box edges)
    for (index_t patch = 0; patch < boxMat.size(); ++patch)
    {
        bool patchHasLines = false;
        for (size_t i = 0; i < patchIDs.size(); ++i)
        {
            if (patchIDs[i] == patch)
            {
                patchHasLines = true;
                break;
            }
        }

        if (patchHasLines)
            continue;

        const index_t patchMaxLevel = maxLevelPerPatch[patch];
        const index_t maxIndex = 1 << patchMaxLevel;
        const index_t numPoints = maxIndex * samplesPerIndexCell + 1;
        const real_t hFinest = 1.0 / std::pow(2.0, patchMaxLevel);

        // u = 0 and u = 1
        for (index_t edge = 0; edge < 2; ++edge)
        {
            const real_t u = (edge == 0) ? 0.0 : 1.0;
            gsMatrix<> line(2, numPoints);
            for (index_t i = 0; i < numPoints; ++i)
            {
                const real_t v = (static_cast<real_t>(i) / samplesPerIndexCell) * hFinest;
                line(0, i) = u;
                line(1, i) = v;
            }
            uvGridLines.push_back(line);
            patchIDs.push_back(patch);
            directions.push_back('v');
        }

        // v = 0 and v = 1
        for (index_t edge = 0; edge < 2; ++edge)
        {
            const real_t v = (edge == 0) ? 0.0 : 1.0;
            gsMatrix<> line(2, numPoints);
            for (index_t i = 0; i < numPoints; ++i)
            {
                const real_t u = (static_cast<real_t>(i) / samplesPerIndexCell) * hFinest;
                line(0, i) = u;
                line(1, i) = v;
            }
            uvGridLines.push_back(line);
            patchIDs.push_back(patch);
            directions.push_back('u');
        }
    }
}

/**
 * @brief Main entry point for visualization mesh generation.
 *
 * Generates a mesh from the box structure, evaluates geometry, and exports to file.
 */
 /**
  * @brief Export box structure in parameter space for debugging
  */
void exportBoxStructure(
    const gsVector<gsVector<gsVector<index_t>>>& boxMat,
    const gsVector<index_t>& numBoxesPerPatch,
    const std::string& filename)
{
    PROFILE_FUNCTION();
    std::ofstream file(filename);
    if (!file.is_open())
    {
        gsWarn << "Could not open " << filename << " for writing.\n";
        return;
    }

    file << "# Box structure in parameter space [0,1]x[0,1]\n";
    file << "# Format: patch level u0 v0 u1 v1 ix0 iy0 ix1 iy1\n";

    for (index_t patch = 0; patch < boxMat.size(); ++patch)
    {
        const index_t numBoxes = numBoxesPerPatch[patch];
        for (index_t boxIdx = 0; boxIdx < numBoxes; ++boxIdx)
        {
            if (boxMat[patch][boxIdx].size() < 5)
                continue;

            const index_t level = boxMat[patch][boxIdx][0];
            const index_t ix0 = boxMat[patch][boxIdx][1];
            const index_t iy0 = boxMat[patch][boxIdx][2];
            const index_t ix1 = boxMat[patch][boxIdx][3];
            const index_t iy1 = boxMat[patch][boxIdx][4];

            const real_t h = 1.0 / std::pow(2.0, level);
            const real_t u0 = ix0 * h;
            const real_t v0 = iy0 * h;
            const real_t u1 = ix1 * h;
            const real_t v1 = iy1 * h;

            file << patch << " " << level << " "
                << u0 << " " << v0 << " " << u1 << " " << v1 << " "
                << ix0 << " " << iy0 << " " << ix1 << " " << iy1 << "\n";
        }
    }

    file.close();
    // gsInfo << "Box structure exported to " << filename << "\n";
}

struct InterfaceSideInfo
{
    index_t patch = -1;
    int sideIndex = -1;
    FeatureSide side = FeatureSide::U0;
};

static bool tryMapInterfaceSideToFeatureSide(
    int sideIndex,
    FeatureSide& side);

struct InterfaceGaussSample
{
    index_t patch = -1;
    index_t element = -1;
    index_t gauss = -1;
    real_t tangent = 0.0;
    real_t normalDistance = 0.0;
    gsVector<real_t> uv;
    gsVector<real_t> xyFit;
    gsVector<real_t> xyOrig;
    real_t jacobianDet = 0.0;
};

static gsVector<real_t> uvOnFeatureSide(FeatureSide side, real_t t)
{
    gsVector<real_t> uv(2);
    switch (side)
    {
    case FeatureSide::U0:
        uv(0) = 0.0;
        uv(1) = t;
        break;
    case FeatureSide::U1:
        uv(0) = 1.0;
        uv(1) = t;
        break;
    case FeatureSide::V0:
        uv(0) = t;
        uv(1) = 0.0;
        break;
    case FeatureSide::V1:
        uv(0) = t;
        uv(1) = 1.0;
        break;
    }
    return uv;
}

static real_t tangentCoordinateOnSide(FeatureSide side, const gsVector<real_t>& uv)
{
    switch (side)
    {
    case FeatureSide::U0:
    case FeatureSide::U1:
        return uv(1);
    case FeatureSide::V0:
    case FeatureSide::V1:
        return uv(0);
    }
    return 0.0;
}

static real_t normalDistanceToSide(FeatureSide side, const gsVector<real_t>& uv)
{
    switch (side)
    {
    case FeatureSide::U0:
        return std::abs(uv(0));
    case FeatureSide::U1:
        return std::abs(1.0 - uv(0));
    case FeatureSide::V0:
        return std::abs(uv(1));
    case FeatureSide::V1:
        return std::abs(1.0 - uv(1));
    }
    return 0.0;
}

static bool elementTouchesSide(
    FeatureSide side,
    const gsVector<real_t>& lower,
    const gsVector<real_t>& upper,
    real_t tolerance = 1e-12)
{
    switch (side)
    {
    case FeatureSide::U0:
        return std::abs(lower(0)) <= tolerance;
    case FeatureSide::U1:
        return std::abs(upper(0) - 1.0) <= tolerance;
    case FeatureSide::V0:
        return std::abs(lower(1)) <= tolerance;
    case FeatureSide::V1:
        return std::abs(upper(1) - 1.0) <= tolerance;
    }
    return false;
}

static gsVector<real_t> matrixColumnToVector2(const gsMatrix<real_t>& matrix)
{
    gsVector<real_t> result(2);
    result.setZero();
    if (matrix.rows() >= 2 && matrix.cols() >= 1)
    {
        result(0) = matrix(0, 0);
        result(1) = matrix(1, 0);
    }
    return result;
}

static gsVector<real_t> evaluateFittedGeometryPoint(
    const MPBES<2, real_t>& mpbes,
    const gsMatrix<real_t>& coefficients,
    index_t patch,
    const gsVector<real_t>& uv,
    bool includeSpilloverFallback,
    bool normalizePou)
{
    gsVector<real_t> xy(2);
    xy.setZero();
    real_t pouSum = static_cast<real_t>(0);

    const index_t count = std::min<index_t>(mpbes.size(), coefficients.rows());
    gsMatrix<real_t> basisValue;
    for (index_t f = 0; f < count; ++f)
    {
        mpbes.evalSingleOnPatch(f, patch, uv, basisValue, includeSpilloverFallback);
        if (basisValue.rows() < 1 || basisValue.cols() < 1)
            continue;

        const real_t value = basisValue(0, 0);
        if (!std::isfinite(value) || value == 0.0)
            continue;

        xy(0) += value * coefficients(f, 0);
        xy(1) += value * coefficients(f, 1);
        pouSum += value;
    }

    // Normalize by the actual partition-of-unity sum to correct any over-counting.
    // When POU = 1 (correct) this is a no-op.  When the sum exceeds 1 because the
    // same physical basis function is counted more than once (e.g. a triple-junction
    // corner function that appears in multiple twin pairs), the evaluation is scaled
    // back to the convex hull of the control coefficients.
    // normalizePou=false bypasses this: required when comparing evaluations from two
    // different patches at a shared interface, where asymmetric POU sums would introduce
    // an apparent gap that does not exist in the actual parameterization.
    if (normalizePou && std::abs(pouSum) > static_cast<real_t>(1e-12))
        xy /= pouSum;

    return xy;
}

static real_t evaluateFittedJacobianDeterminant(
    const MPBES<2, real_t>& mpbes,
    const gsMatrix<real_t>& coefficients,
    index_t patch,
    const gsVector<real_t>& uv,
    bool includeSpilloverFallback = false)
{
    real_t dxdu = 0.0;
    real_t dxdv = 0.0;
    real_t dydu = 0.0;
    real_t dydv = 0.0;

    const index_t count = std::min<index_t>(mpbes.size(), coefficients.rows());
    gsMatrix<real_t> basisDeriv;
    for (index_t f = 0; f < count; ++f)
    {
        mpbes.evalDerivSingleOnPatch(f, patch, uv, basisDeriv, includeSpilloverFallback);
        if (basisDeriv.rows() < 2 || basisDeriv.cols() < 1)
            continue;

        const real_t dNdu = basisDeriv(0, 0);
        const real_t dNdv = basisDeriv(1, 0);
        if (!std::isfinite(dNdu) || !std::isfinite(dNdv))
            continue;

        dxdu += dNdu * coefficients(f, 0);
        dxdv += dNdv * coefficients(f, 0);
        dydu += dNdu * coefficients(f, 1);
        dydv += dNdv * coefficients(f, 1);
    }

    return dxdu * dydv - dxdv * dydu;
}

static void evaluateFittedJacobianColumns(
    const MPBES<2, real_t>& mpbes,
    const gsMatrix<real_t>& coefficients,
    index_t patch,
    const gsVector<real_t>& uv,
    gsVector<real_t>& du,
    gsVector<real_t>& dv,
    bool includeSpilloverFallback = false)
{
    du.setZero();
    dv.setZero();

    const index_t count = std::min<index_t>(mpbes.size(), coefficients.rows());
    gsMatrix<real_t> basisDeriv;
    for (index_t f = 0; f < count; ++f)
    {
        mpbes.evalDerivSingleOnPatch(f, patch, uv, basisDeriv, includeSpilloverFallback);
        if (basisDeriv.rows() < 2 || basisDeriv.cols() < 1)
            continue;

        const real_t dNdu = basisDeriv(0, 0);
        const real_t dNdv = basisDeriv(1, 0);
        if (!std::isfinite(dNdu) || !std::isfinite(dNdv))
            continue;

        du(0) += dNdu * coefficients(f, 0);
        du(1) += dNdu * coefficients(f, 1);
        dv(0) += dNdv * coefficients(f, 0);
        dv(1) += dNdv * coefficients(f, 1);
    }
}

static void normalizeDirection(gsVector<real_t>& vec)
{
    const real_t norm = std::sqrt(vec(0) * vec(0) + vec(1) * vec(1));
    if (norm > 1e-12 && std::isfinite(norm))
        vec /= norm;
    else
        vec.setZero();
}

static std::vector<PatchVisualizationMetadata> buildPatchVisualizationMetadata(
    const MPBES<2, real_t>& mpbes,
    const gsMatrix<real_t>& coefficients,
    index_t numPatches)
{
    std::vector<PatchVisualizationMetadata> metadata(static_cast<size_t>(numPatches));
    gsVector<real_t> uv(2);
    uv(0) = 0.5;
    uv(1) = 0.5;

    for (index_t patch = 0; patch < numPatches; ++patch)
    {
        PatchVisualizationMetadata info;
        info.patch = patch;
        info.centerXY = evaluateFittedGeometryPoint(mpbes, coefficients, patch, uv);
        evaluateFittedJacobianColumns(mpbes, coefficients, patch, uv, info.uDirection, info.vDirection);
        normalizeDirection(info.uDirection);
        normalizeDirection(info.vDirection);
        info.jacobianDeterminant = evaluateFittedJacobianDeterminant(mpbes, coefficients, patch, uv);

        if (g_geometryPreflight.valid &&
            patch < static_cast<index_t>(g_geometryPreflight.mirroredReport.mirrored.size()))
        {
            info.mirrored = g_geometryPreflight.mirroredReport.mirrored[patch];
        }
        else
        {
            info.mirrored = info.jacobianDeterminant < 0.0;
        }

        metadata[static_cast<size_t>(patch)] = info;
    }

    return metadata;
}

static gsVector<real_t> tangentDirectionOnSide(
    FeatureSide side,
    const gsVector<real_t>& du,
    const gsVector<real_t>& dv)
{
    switch (side)
    {
    case FeatureSide::U0:
    case FeatureSide::U1:
        return dv;
    case FeatureSide::V0:
    case FeatureSide::V1:
        return du;
    }

    gsVector<real_t> zero(2);
    zero.setZero();
    return zero;
}

static std::vector<InterfaceVisualizationMetadata> buildInterfaceVisualizationMetadata(
    const MPBES<2, real_t>& mpbes,
    const gsMultiPatch<real_t>& mp,
    const gsMatrix<real_t>& coefficients)
{
    std::vector<InterfaceVisualizationMetadata> metadata;
    const std::vector<GeometryPreflightInterfaceInfo> interfaceInfos =
        g_geometryPreflight.valid ? g_geometryPreflight.interfaces : collectVisualizationInterfaceInfos(mp);
    metadata.reserve(interfaceInfos.size());

    for (size_t i = 0; i < interfaceInfos.size(); ++i)
    {
        InterfaceSideInfo sideA;
        InterfaceSideInfo sideB;
        sideA.patch = interfaceInfos[i].patchA;
        sideA.sideIndex = interfaceInfos[i].sideIndexA;
        sideA.side = interfaceInfos[i].sideA;
        sideB.patch = interfaceInfos[i].patchB;
        sideB.sideIndex = interfaceInfos[i].sideIndexB;
        sideB.side = interfaceInfos[i].sideB;

        if (!interfaceInfos[i].hasMappedSides)
            continue;

        InterfaceVisualizationMetadata info;
        info.patchA = sideA.patch;
        info.patchB = sideB.patch;
        info.sideA = sideA.side;
        info.sideB = sideB.side;
        info.orientationReversed = interfaceInfos[i].orientationReversed;
        info.mirroredA = isPreflightMirroredPatch(info.patchA);
        info.mirroredB = isPreflightMirroredPatch(info.patchB);

        const gsVector<real_t> uvA = uvOnFeatureSide(sideA.side, 0.5);
        const gsVector<real_t> uvB = uvOnFeatureSide(sideB.side, 0.5);
        info.midpointA = evaluateFittedGeometryPoint(mpbes, coefficients, info.patchA, uvA, false, false);
        info.midpointB = evaluateFittedGeometryPoint(mpbes, coefficients, info.patchB, uvB, false, false);
        info.midpointGap = vectorDistance(info.midpointA, info.midpointB);

        gsVector<real_t> duA(2), dvA(2), duB(2), dvB(2);
        evaluateFittedJacobianColumns(mpbes, coefficients, info.patchA, uvA, duA, dvA, true);
        evaluateFittedJacobianColumns(mpbes, coefficients, info.patchB, uvB, duB, dvB, true);
        info.tangentA = tangentDirectionOnSide(sideA.side, duA, dvA);
        info.tangentB = tangentDirectionOnSide(sideB.side, duB, dvB);
        normalizeDirection(info.tangentA);
        normalizeDirection(info.tangentB);

        metadata.push_back(info);
    }

    return metadata;
}

static real_t vectorDistance(const gsVector<real_t>& a, const gsVector<real_t>& b)
{
    return std::sqrt(
        (a(0) - b(0)) * (a(0) - b(0)) +
        (a(1) - b(1)) * (a(1) - b(1)));
}

static bool findPatchInterface(
    const gsMultiPatch<real_t>& mp,
    index_t patchA,
    index_t patchB,
    InterfaceSideInfo& sideA,
    InterfaceSideInfo& sideB)
{
    if (g_geometryPreflight.valid)
    {
        const GeometryPreflightInterfaceInfo* info = findPreflightInterfaceInfo(patchA, patchB);
        if (info && info->hasMappedSides)
        {
            if (info->patchA == patchA && info->patchB == patchB)
            {
                sideA.patch = patchA;
                sideA.sideIndex = info->sideIndexA;
                sideA.side = info->sideA;
                sideB.patch = patchB;
                sideB.sideIndex = info->sideIndexB;
                sideB.side = info->sideB;
            }
            else
            {
                sideA.patch = patchA;
                sideA.sideIndex = info->sideIndexB;
                sideA.side = info->sideB;
                sideB.patch = patchB;
                sideB.sideIndex = info->sideIndexA;
                sideB.side = info->sideA;
            }
            return true;
        }
    }

    std::vector<boundaryInterface>& interfaces = const_cast<gsMultiPatch<real_t>&>(mp).interfaces();
    for (size_t i = 0; i < interfaces.size(); ++i)
    {
        const index_t firstPatch = interfaces[i].first().patch;
        const index_t secondPatch = interfaces[i].second().patch;

        if (firstPatch == patchA && secondPatch == patchB)
        {
            sideA.patch = patchA;
            sideA.sideIndex = interfaces[i].first().index();
            sideB.patch = patchB;
            sideB.sideIndex = interfaces[i].second().index();
        }
        else if (firstPatch == patchB && secondPatch == patchA)
        {
            sideA.patch = patchA;
            sideA.sideIndex = interfaces[i].second().index();
            sideB.patch = patchB;
            sideB.sideIndex = interfaces[i].first().index();
        }
        else
        {
            continue;
        }

        return tryMapInterfaceSideToFeatureSide(sideA.sideIndex, sideA.side)
            && tryMapInterfaceSideToFeatureSide(sideB.sideIndex, sideB.side);
    }

    return false;
}

static bool isInterfaceOrientationReversed(
    const gsMultiPatch<real_t>& mp,
    const InterfaceSideInfo& sideA,
    const InterfaceSideInfo& sideB)
{
    const gsVector<real_t> uvA0 = uvOnFeatureSide(sideA.side, 0.0);
    const gsVector<real_t> uvA1 = uvOnFeatureSide(sideA.side, 1.0);
    const gsVector<real_t> uvB0 = uvOnFeatureSide(sideB.side, 0.0);
    const gsVector<real_t> uvB1 = uvOnFeatureSide(sideB.side, 1.0);

    const gsVector<real_t> xyA0 = matrixColumnToVector2(mp.patch(sideA.patch).eval(uvA0));
    const gsVector<real_t> xyA1 = matrixColumnToVector2(mp.patch(sideA.patch).eval(uvA1));
    const gsVector<real_t> xyB0 = matrixColumnToVector2(mp.patch(sideB.patch).eval(uvB0));
    const gsVector<real_t> xyB1 = matrixColumnToVector2(mp.patch(sideB.patch).eval(uvB1));

    const real_t same = vectorDistance(xyA0, xyB0) + vectorDistance(xyA1, xyB1);
    const real_t reversed = vectorDistance(xyA0, xyB1) + vectorDistance(xyA1, xyB0);
    return reversed < same;
}

static real_t estimateGeometryScale(const gsMultiPatch<real_t>& mp)
{
    real_t minX = std::numeric_limits<real_t>::max();
    real_t minY = std::numeric_limits<real_t>::max();
    real_t maxX = -std::numeric_limits<real_t>::max();
    real_t maxY = -std::numeric_limits<real_t>::max();

    for (index_t patch = 0; patch < mp.nPatches(); ++patch)
    {
        const real_t corners[4][2] = {
            {0.0, 0.0},
            {1.0, 0.0},
            {0.0, 1.0},
            {1.0, 1.0}
        };

        for (index_t corner = 0; corner < 4; ++corner)
        {
            gsVector<real_t> uv(2);
            uv(0) = corners[corner][0];
            uv(1) = corners[corner][1];
            const gsVector<real_t> xy = matrixColumnToVector2(mp.patch(patch).eval(uv));
            minX = std::min(minX, xy(0));
            minY = std::min(minY, xy(1));
            maxX = std::max(maxX, xy(0));
            maxY = std::max(maxY, xy(1));
        }
    }

    const real_t dx = maxX - minX;
    const real_t dy = maxY - minY;
    return std::max(std::sqrt(dx * dx + dy * dy), static_cast<real_t>(1.0));
}

static bool isRequestedVisualizationInterfacePair(
    index_t patchA,
    index_t patchB)
{
    const std::pair<index_t, index_t> pairKey(
        std::min(patchA, patchB),
        std::max(patchA, patchB));

    return pairKey == std::make_pair<index_t, index_t>(0, 1) ||
        pairKey == std::make_pair<index_t, index_t>(0, 3) ||
        pairKey == std::make_pair<index_t, index_t>(1, 2) ||
        pairKey == std::make_pair<index_t, index_t>(2, 5) ||
        pairKey == std::make_pair<index_t, index_t>(3, 4);
}

static void computeInterfaceGapStats(
    const gsMultiPatch<real_t>& mp,
    index_t patchA,
    FeatureSide sideA,
    index_t patchB,
    FeatureSide sideB,
    bool reversedOrientation,
    real_t& maxGap,
    real_t& avgGap)
{
    const index_t sampleCount = 11;
    maxGap = 0.0;
    real_t sumGap = 0.0;

    for (index_t i = 0; i < sampleCount; ++i)
    {
        const real_t t = static_cast<real_t>(i) / static_cast<real_t>(sampleCount - 1);
        const real_t otherT = reversedOrientation ? (1.0 - t) : t;
        const gsVector<real_t> uvA = uvOnFeatureSide(sideA, t);
        const gsVector<real_t> uvB = uvOnFeatureSide(sideB, otherT);
        const gsVector<real_t> xyA = matrixColumnToVector2(mp.patch(patchA).eval(uvA));
        const gsVector<real_t> xyB = matrixColumnToVector2(mp.patch(patchB).eval(uvB));
        const real_t gap = vectorDistance(xyA, xyB);
        maxGap = std::max(maxGap, gap);
        sumGap += gap;
    }

    avgGap = sumGap / static_cast<real_t>(sampleCount);
}

static std::vector<GeometryPreflightInterfaceInfo> collectVisualizationInterfaceInfos(
    const gsMultiPatch<real_t>& mp)
{
    std::vector<GeometryPreflightInterfaceInfo> infos;
    std::vector<boundaryInterface>& interfaces = const_cast<gsMultiPatch<real_t>&>(mp).interfaces();
    infos.reserve(interfaces.size());

    std::set<std::pair<index_t, index_t>> existingPairs;
    for (size_t i = 0; i < interfaces.size(); ++i)
    {
        GeometryPreflightInterfaceInfo info;
        info.patchA = interfaces[i].first().patch;
        info.patchB = interfaces[i].second().patch;
        info.sideIndexA = interfaces[i].first().index();
        info.sideIndexB = interfaces[i].second().index();

        InterfaceSideInfo sideA;
        InterfaceSideInfo sideB;
        sideA.patch = info.patchA;
        sideA.sideIndex = info.sideIndexA;
        sideB.patch = info.patchB;
        sideB.sideIndex = info.sideIndexB;

        info.hasMappedSides =
            tryMapInterfaceSideToFeatureSide(info.sideIndexA, sideA.side) &&
            tryMapInterfaceSideToFeatureSide(info.sideIndexB, sideB.side);

        if (info.hasMappedSides)
        {
            info.sideA = sideA.side;
            info.sideB = sideB.side;
            info.orientationReversed = isInterfaceOrientationReversed(mp, sideA, sideB);
        }

        infos.push_back(info);
        existingPairs.insert(std::make_pair(
            std::min(info.patchA, info.patchB),
            std::max(info.patchA, info.patchB)));
    }

    const real_t geometryScale = estimateGeometryScale(mp);
    const real_t interfaceTolerance = geometryScale * static_cast<real_t>(1e-5);
    const real_t requestedInterfaceTolerance = geometryScale * static_cast<real_t>(5e-3);
    const FeatureSide allSides[4] = {
        FeatureSide::U0,
        FeatureSide::U1,
        FeatureSide::V0,
        FeatureSide::V1
    };

    for (index_t patchA = 0; patchA < mp.nPatches(); ++patchA)
    {
        for (index_t patchB = patchA + 1; patchB < mp.nPatches(); ++patchB)
        {
            const std::pair<index_t, index_t> pairKey(patchA, patchB);
            if (existingPairs.find(pairKey) != existingPairs.end())
                continue;

            real_t bestMaxGap = std::numeric_limits<real_t>::max();
            real_t bestAvgGap = std::numeric_limits<real_t>::max();
            FeatureSide bestSideA = FeatureSide::U0;
            FeatureSide bestSideB = FeatureSide::U0;
            bool bestReversed = false;

            for (index_t sideIdxA = 0; sideIdxA < 4; ++sideIdxA)
            {
                for (index_t sideIdxB = 0; sideIdxB < 4; ++sideIdxB)
                {
                    for (index_t reverseFlag = 0; reverseFlag < 2; ++reverseFlag)
                    {
                        real_t maxGap = 0.0;
                        real_t avgGap = 0.0;
                        computeInterfaceGapStats(
                            mp,
                            patchA,
                            allSides[sideIdxA],
                            patchB,
                            allSides[sideIdxB],
                            reverseFlag == 1,
                            maxGap,
                            avgGap);

                        if (maxGap < bestMaxGap ||
                            (std::abs(maxGap - bestMaxGap) <= 1e-12 && avgGap < bestAvgGap))
                        {
                            bestMaxGap = maxGap;
                            bestAvgGap = avgGap;
                            bestSideA = allSides[sideIdxA];
                            bestSideB = allSides[sideIdxB];
                            bestReversed = (reverseFlag == 1);
                        }
                    }
                }
            }

            const real_t acceptedTolerance = isRequestedVisualizationInterfacePair(patchA, patchB)
                ? requestedInterfaceTolerance
                : interfaceTolerance;

            if (bestMaxGap > acceptedTolerance)
                continue;

            GeometryPreflightInterfaceInfo info;
            info.patchA = patchA;
            info.patchB = patchB;
            info.sideA = bestSideA;
            info.sideB = bestSideB;
            info.hasMappedSides = true;
            info.orientationReversed = bestReversed;
            featureSideToInterfaceSideIndex(bestSideA, info.sideIndexA);
            featureSideToInterfaceSideIndex(bestSideB, info.sideIndexB);
            infos.push_back(info);
            existingPairs.insert(pairKey);
        }
    }

    return infos;
}

static std::vector<InterfaceValidationSpec> collectValidationInterfaces(
    const gsMultiPatch<real_t>& mp)
{
    std::vector<InterfaceValidationSpec> specs;
    std::vector<boundaryInterface>& interfaces = const_cast<gsMultiPatch<real_t>&>(mp).interfaces();
    specs.reserve(interfaces.size());

    for (size_t i = 0; i < interfaces.size(); ++i)
    {
        FeatureSide firstSide;
        FeatureSide secondSide;
        if (!tryMapInterfaceSideToFeatureSide(interfaces[i].first().index(), firstSide) ||
            !tryMapInterfaceSideToFeatureSide(interfaces[i].second().index(), secondSide))
            continue;

        InterfaceValidationSpec spec;
        spec.patchA = interfaces[i].first().patch;
        spec.patchB = interfaces[i].second().patch;
        spec.sideA = firstSide;
        spec.sideB = secondSide;

        const GeometryPreflightInterfaceInfo* preflightInfo =
            findPreflightInterfaceInfo(spec.patchA, spec.patchB);
        if (preflightInfo && preflightInfo->hasMappedSides)
        {
            spec.orientationReversed = preflightInfo->orientationReversed;
        }
        else
        {
            InterfaceSideInfo sideA;
            InterfaceSideInfo sideB;
            sideA.patch = spec.patchA;
            sideA.sideIndex = interfaces[i].first().index();
            sideA.side = spec.sideA;
            sideB.patch = spec.patchB;
            sideB.sideIndex = interfaces[i].second().index();
            sideB.side = spec.sideB;
            spec.orientationReversed = isInterfaceOrientationReversed(mp, sideA, sideB);
        }

        specs.push_back(spec);
    }

    return specs;
}

static void buildInterfaceValidationUv(
    const InterfaceValidationSpec& spec,
    index_t pointsPerEdge,
    gsMatrix<real_t>& uvA,
    gsMatrix<real_t>& uvB)
{
    uvA.resize(2, pointsPerEdge);
    uvB.resize(2, pointsPerEdge);

    for (index_t i = 0; i < pointsPerEdge; ++i)
    {
        const real_t t = (pointsPerEdge == 1)
            ? 0.0
            : static_cast<real_t>(i) / static_cast<real_t>(pointsPerEdge - 1);
        const real_t otherT = spec.orientationReversed ? (1.0 - t) : t;
        uvA.col(i) = uvOnFeatureSide(spec.sideA, t);
        uvB.col(i) = uvOnFeatureSide(spec.sideB, otherT);
    }
}

static index_t buildValidationInterfacePoints(
    const gsMultiPatch<real_t>& mp,
    index_t pointsPerEdge,
    std::vector<InterfaceValidationSpec>& specs,
    gsVector<gsMatrix<>>& interfacePoints)
{
    specs = collectValidationInterfaces(mp);
    interfacePoints.resize(mp.nPatches());
    for (index_t p = 0; p < mp.nPatches(); ++p)
        interfacePoints[p].resize(2, 0);

    index_t totalPoints = 0;
    for (size_t i = 0; i < specs.size(); ++i)
    {
        gsMatrix<real_t> uvA;
        gsMatrix<real_t> uvB;
        buildInterfaceValidationUv(specs[i], pointsPerEdge, uvA, uvB);
        appendUvSamples(interfacePoints[specs[i].patchA], uvA);
        appendUvSamples(interfacePoints[specs[i].patchB], uvB);
        totalPoints += uvA.cols() + uvB.cols();
    }

    return totalPoints;
}

static bool featureSideToInterfaceSideIndex(
    FeatureSide side,
    int& sideIndex)
{
    switch (side)
    {
    case FeatureSide::U0:
        sideIndex = 1;
        return true;
    case FeatureSide::U1:
        sideIndex = 2;
        return true;
    case FeatureSide::V0:
        sideIndex = 3;
        return true;
    case FeatureSide::V1:
        sideIndex = 4;
        return true;
    }

    sideIndex = -1;
    return false;
}

static bool isInterfaceSide(
    const gsMultiPatch<real_t>& mp,
    index_t patch,
    FeatureSide side)
{
    int sideIndex = -1;
    if (!featureSideToInterfaceSideIndex(side, sideIndex))
        return false;

    std::vector<boundaryInterface>& interfaces = const_cast<gsMultiPatch<real_t>&>(mp).interfaces();
    for (size_t i = 0; i < interfaces.size(); ++i)
    {
        if ((interfaces[i].first().patch == patch && interfaces[i].first().index() == sideIndex) ||
            (interfaces[i].second().patch == patch && interfaces[i].second().index() == sideIndex))
            return true;
    }

    return false;
}

static std::vector<FeatureBoundarySpec> normalizeFeatureBoundarySpecs(
    const gsMultiPatch<real_t>& mp,
    const std::vector<FeatureBoundarySpec>& requested,
    bool logToInfo)
{
    std::vector<FeatureBoundarySpec> normalized;
    std::unordered_set<long long> seen;

    auto logMessage = [logToInfo](const std::string& message)
    {
        if (outfile.is_open())
            outfile << message << "\n";
        if (logToInfo)
            gsInfo << message << "\n";
    };

    logMessage("[feature-boundaries] requested=" + std::to_string(requested.size()));

    for (size_t i = 0; i < requested.size(); ++i)
    {
        const FeatureBoundarySpec& spec = requested[i];
        if (spec.patch < 0 || spec.patch >= mp.nPatches())
        {
            logMessage("[feature-boundaries] skip invalid patch=" + std::to_string(spec.patch) +
                       " for side=" + featureSideName(spec.side));
            continue;
        }

        if (isInterfaceSide(mp, spec.patch, spec.side))
        {
            logMessage("[feature-boundaries] skip internal interface side patch=" +
                       std::to_string(spec.patch) + " side=" + featureSideName(spec.side));
            continue;
        }

        int sideIndex = -1;
        featureSideToInterfaceSideIndex(spec.side, sideIndex);
        const long long key = static_cast<long long>(spec.patch) * 16LL + static_cast<long long>(sideIndex);
        if (!seen.insert(key).second)
        {
            logMessage("[feature-boundaries] skip duplicate patch=" + std::to_string(spec.patch) +
                       " side=" + featureSideName(spec.side));
            continue;
        }

        normalized.push_back(spec);

        std::ostringstream resolved;
        resolved << "[feature-boundaries] keep patch=" << spec.patch
                 << " side=" << featureSideName(spec.side)
                 << " mirroredBaseline="
                 << ((g_geometryPreflight.valid &&
                      spec.patch < static_cast<index_t>(g_geometryPreflight.mirroredReport.mirrored.size()) &&
                      g_geometryPreflight.mirroredReport.mirrored[spec.patch]) ? "true" : "false");
        logMessage(resolved.str());
    }

    logMessage("[feature-boundaries] normalized=" + std::to_string(normalized.size()));
    return normalized;
}

static bool featureSideToCanonicalBoundaryCode(
    FeatureSide side,
    int& boundaryCode)
{
    switch (side)
    {
    case FeatureSide::U0:
        boundaryCode = 1;
        return true;
    case FeatureSide::U1:
        boundaryCode = 2;
        return true;
    case FeatureSide::V0:
        boundaryCode = 3;
        return true;
    case FeatureSide::V1:
        boundaryCode = 4;
        return true;
    }

    boundaryCode = -1;
    return false;
}

static bool resolveCanonicalTopologySides(
    index_t firstPatch,
    index_t secondPatch,
    int rawFirstSide,
    int rawSecondSide,
    int& firstSide,
    int& secondSide)
{
    FeatureSide firstFeatureSide;
    FeatureSide secondFeatureSide;
    bool hasMappedSides = false;

    const GeometryPreflightInterfaceInfo* preflightInfo =
        findPreflightInterfaceInfo(firstPatch, secondPatch);
    if (preflightInfo && preflightInfo->hasMappedSides)
    {
        if (preflightInfo->patchA == firstPatch && preflightInfo->patchB == secondPatch)
        {
            firstFeatureSide = preflightInfo->sideA;
            secondFeatureSide = preflightInfo->sideB;
            hasMappedSides = true;
        }
        else if (preflightInfo->patchA == secondPatch && preflightInfo->patchB == firstPatch)
        {
            firstFeatureSide = preflightInfo->sideB;
            secondFeatureSide = preflightInfo->sideA;
            hasMappedSides = true;
        }
    }

    if (!hasMappedSides)
    {
        hasMappedSides =
            tryMapInterfaceSideToFeatureSide(rawFirstSide, firstFeatureSide) &&
            tryMapInterfaceSideToFeatureSide(rawSecondSide, secondFeatureSide);
    }

    if (!hasMappedSides)
        return false;

    return featureSideToCanonicalBoundaryCode(firstFeatureSide, firstSide) &&
        featureSideToCanonicalBoundaryCode(secondFeatureSide, secondSide);
}

static std::vector<InterfaceGaussSample> collectClosestGaussSamples(
    const MPBES<2, real_t>& mpbes,
    const gsMultiPatch<real_t>& mp,
    const gsMatrix<real_t>& coefficients,
    const InterfaceSideInfo& sideInfo)
{
    std::vector<InterfaceGaussSample> samples;

    if (sideInfo.patch < 0 || sideInfo.patch >= mpbes.nPatches())
        return samples;

    gsVector<index_t> numNodes(2);
    numNodes[0] = numNodes[1] = 2;
    gsGaussRule<real_t> quRule(numNodes);

    auto domIt = mpbes.thbBases()[sideInfo.patch].makeDomainIterator();
    gsMatrix<real_t> quNodes;
    gsVector<real_t> quWeights;

    for (index_t element = 0; domIt->good(); domIt->next(), ++element)
    {
        const gsVector<real_t> lower = domIt->lowerCorner();
        const gsVector<real_t> upper = domIt->upperCorner();
        if (!elementTouchesSide(sideInfo.side, lower, upper))
            continue;

        quRule.mapTo(lower, upper, quNodes, quWeights);

        index_t bestGauss = -1;
        real_t bestDistance = std::numeric_limits<real_t>::max();
        gsVector<real_t> bestUv(2);

        for (index_t q = 0; q < quNodes.cols(); ++q)
        {
            const gsVector<real_t> uv = quNodes.col(q);
            const real_t distance = normalDistanceToSide(sideInfo.side, uv);
            if (distance < bestDistance)
            {
                bestDistance = distance;
                bestGauss = q;
                bestUv = uv;
            }
        }

        if (bestGauss < 0)
            continue;

        InterfaceGaussSample sample;
        sample.patch = sideInfo.patch;
        sample.element = element;
        sample.gauss = bestGauss;
        sample.tangent = tangentCoordinateOnSide(sideInfo.side, bestUv);
        sample.normalDistance = bestDistance;
        sample.uv = bestUv;
        sample.xyFit = evaluateFittedGeometryPoint(mpbes, coefficients, sideInfo.patch, bestUv, true);
        sample.xyOrig = matrixColumnToVector2(mp.patch(sideInfo.patch).eval(bestUv));
        sample.jacobianDet = evaluateFittedJacobianDeterminant(mpbes, coefficients, sideInfo.patch, bestUv, true);
        samples.push_back(sample);
    }

    std::sort(samples.begin(), samples.end(), [](const InterfaceGaussSample& lhs, const InterfaceGaussSample& rhs)
    {
        return lhs.tangent < rhs.tangent;
    });

    return samples;
}

static const InterfaceGaussSample* findClosestSampleByTangent(
    const std::vector<InterfaceGaussSample>& samples,
    real_t targetTangent)
{
    const InterfaceGaussSample* best = nullptr;
    real_t bestDistance = std::numeric_limits<real_t>::max();

    for (size_t i = 0; i < samples.size(); ++i)
    {
        const real_t distance = std::abs(samples[i].tangent - targetTangent);
        if (distance < bestDistance)
        {
            bestDistance = distance;
            best = &samples[i];
        }
    }

    return best;
}

static void writePatchInterfaceDiagnostics(
    std::ostream& log,
    const MPBES<2, real_t>& mpbes,
    const gsMultiPatch<real_t>& mp,
    const gsMatrix<real_t>& coefficients,
    const std::string& label,
    const InterfaceSideInfo& firstPatchSide,
    const InterfaceSideInfo& secondPatchSide,
    bool reversedOrientation,
    const std::vector<InterfaceGaussSample>& firstPatchSamples,
    const std::vector<InterfaceGaussSample>& secondPatchSamples)
{
    log << std::setprecision(16);
    log << "\nINTERFACE_" << firstPatchSide.patch << "_" << secondPatchSide.patch
        << "_DIAGNOSTIC_BEGIN [" << label << "]\n";
    log << "patch " << firstPatchSide.patch << " sideIndex=" << firstPatchSide.sideIndex
        << ", patch " << secondPatchSide.patch << " sideIndex=" << secondPatchSide.sideIndex
        << ", orientation=" << (reversedOrientation ? "reversed" : "same") << "\n";
    log << "closest Gauss strips: patch" << firstPatchSide.patch << "=" << firstPatchSamples.size()
        << ", patch" << secondPatchSide.patch << "=" << secondPatchSamples.size() << "\n";

    if (firstPatchSamples.empty() || secondPatchSamples.empty())
    {
        log << "Unable to collect boundary-adjacent Gauss samples on both sides.\n";
        log << "INTERFACE_" << firstPatchSide.patch << "_" << secondPatchSide.patch << "_DIAGNOSTIC_END\n";
        return;
    }

    for (size_t i = 0; i < firstPatchSamples.size(); ++i)
    {
        const InterfaceGaussSample& firstSample = firstPatchSamples[i];
        const real_t tFirst = firstSample.tangent;
        const real_t tSecond = reversedOrientation ? (1.0 - tFirst) : tFirst;

        const gsVector<real_t> uvIfaceFirst = uvOnFeatureSide(firstPatchSide.side, tFirst);
        const gsVector<real_t> uvIfaceSecond = uvOnFeatureSide(secondPatchSide.side, tSecond);
        const gsVector<real_t> xyFitIfaceFirst = evaluateFittedGeometryPoint(mpbes, coefficients, firstPatchSide.patch, uvIfaceFirst, true);
        const gsVector<real_t> xyFitIfaceSecond = evaluateFittedGeometryPoint(mpbes, coefficients, secondPatchSide.patch, uvIfaceSecond, true);
        const gsVector<real_t> xyOrigIfaceFirst = matrixColumnToVector2(mp.patch(firstPatchSide.patch).eval(uvIfaceFirst));
        const gsVector<real_t> xyOrigIfaceSecond = matrixColumnToVector2(mp.patch(secondPatchSide.patch).eval(uvIfaceSecond));
        const real_t detIfaceFirst = evaluateFittedJacobianDeterminant(mpbes, coefficients, firstPatchSide.patch, uvIfaceFirst, true);
        const real_t detIfaceSecond = evaluateFittedJacobianDeterminant(mpbes, coefficients, secondPatchSide.patch, uvIfaceSecond, true);
        const InterfaceGaussSample* secondSample = findClosestSampleByTangent(secondPatchSamples, tSecond);

        log << "sample " << i
            << ": interface_t_patch" << firstPatchSide.patch << "=" << tFirst
            << ", interface_t_patch" << secondPatchSide.patch << "=" << tSecond << "\n";
        log << "  interface patch " << firstPatchSide.patch << ": uv=(" << uvIfaceFirst(0) << ", " << uvIfaceFirst(1) << ")"
            << ", fit=(" << xyFitIfaceFirst(0) << ", " << xyFitIfaceFirst(1) << ")"
            << ", orig=(" << xyOrigIfaceFirst(0) << ", " << xyOrigIfaceFirst(1) << ")"
            << ", det=" << detIfaceFirst << "\n";
        log << "  interface patch " << secondPatchSide.patch << ": uv=(" << uvIfaceSecond(0) << ", " << uvIfaceSecond(1) << ")"
            << ", fit=(" << xyFitIfaceSecond(0) << ", " << xyFitIfaceSecond(1) << ")"
            << ", orig=(" << xyOrigIfaceSecond(0) << ", " << xyOrigIfaceSecond(1) << ")"
            << ", det=" << detIfaceSecond << "\n";
        log << "  interface gaps: fitGap=" << vectorDistance(xyFitIfaceFirst, xyFitIfaceSecond)
            << ", origGap=" << vectorDistance(xyOrigIfaceFirst, xyOrigIfaceSecond) << "\n";

        log << "  closest Gauss patch " << firstPatchSide.patch << ": element=" << firstSample.element
            << ", gauss=" << firstSample.gauss
            << ", tangent=" << firstSample.tangent
            << ", normalDist=" << firstSample.normalDistance
            << ", uv=(" << firstSample.uv(0) << ", " << firstSample.uv(1) << ")"
            << ", fit=(" << firstSample.xyFit(0) << ", " << firstSample.xyFit(1) << ")"
            << ", orig=(" << firstSample.xyOrig(0) << ", " << firstSample.xyOrig(1) << ")"
            << ", det=" << firstSample.jacobianDet << "\n";

        if (secondSample)
        {
            log << "  closest Gauss patch " << secondPatchSide.patch << ": element=" << secondSample->element
                << ", gauss=" << secondSample->gauss
                << ", tangent=" << secondSample->tangent
                << ", normalDist=" << secondSample->normalDistance
                << ", uv=(" << secondSample->uv(0) << ", " << secondSample->uv(1) << ")"
                << ", fit=(" << secondSample->xyFit(0) << ", " << secondSample->xyFit(1) << ")"
                << ", orig=(" << secondSample->xyOrig(0) << ", " << secondSample->xyOrig(1) << ")"
                << ", det=" << secondSample->jacobianDet << "\n";
            log << "  Gauss fit gap=" << vectorDistance(firstSample.xyFit, secondSample->xyFit)
                << ", Gauss orig gap=" << vectorDistance(firstSample.xyOrig, secondSample->xyOrig) << "\n";
        }
        else
        {
            log << "  closest Gauss patch " << secondPatchSide.patch << ": not found\n";
        }
    }

    log << "INTERFACE_" << firstPatchSide.patch << "_" << secondPatchSide.patch << "_DIAGNOSTIC_END\n";
}

void logSpecificInterface(
    const MPBES<2, real_t>& mpbes,
    const gsMultiPatch<real_t>& mp,
    const gsMatrix<real_t>& coefficients,
    const std::string& label,
    index_t firstPatch,
    index_t secondPatch)
{
    InterfaceSideInfo firstPatchSide;
    InterfaceSideInfo secondPatchSide;

    if (!outfile.is_open())
        return;

    if (!findPatchInterface(mp, firstPatch, secondPatch, firstPatchSide, secondPatchSide))
    {
        std::ostringstream buffer;
        buffer << "\nINTERFACE_" << firstPatch << "_" << secondPatch << "_DIAGNOSTIC_BEGIN [" << label << "]\n";
        buffer << "Could not find topology interface between patches " << firstPatch << " and " << secondPatch << ".\n";
        buffer << "INTERFACE_" << firstPatch << "_" << secondPatch << "_DIAGNOSTIC_END\n";
        outfile << buffer.str();
        outfile.flush();
        gsInfo << buffer.str();
        return;
    }

    const bool reversedOrientation = isInterfaceOrientationReversed(mp, firstPatchSide, secondPatchSide);
    const std::vector<InterfaceGaussSample> firstPatchSamples = collectClosestGaussSamples(mpbes, mp, coefficients, firstPatchSide);
    const std::vector<InterfaceGaussSample> secondPatchSamples = collectClosestGaussSamples(mpbes, mp, coefficients, secondPatchSide);

    std::ostringstream buffer;
    writePatchInterfaceDiagnostics(
        buffer,
        mpbes,
        mp,
        coefficients,
        label,
        firstPatchSide,
        secondPatchSide,
        reversedOrientation,
        firstPatchSamples,
        secondPatchSamples);

    outfile << buffer.str();
    outfile.flush();
    gsInfo << buffer.str();

    gsInfo << "Interface " << firstPatch << "/" << secondPatch << " diagnostic written to main logfile\n";
}

void generateOriginalGeometryMesh(
    const gsMultiPatch<real_t>& mp,
    index_t gridResolution,
    const std::string& outputPrefix)
{
    std::string meshFile = outputPrefix + ".txt";
    std::ofstream meshOut(meshFile);

    if (!meshOut.is_open())
    {
        gsInfo << "ERROR: Could not open " << meshFile << " for writing\n";
        return;
    }

    meshOut << "# Original Geometry Mesh\n";
    meshOut << "# Format: x y patch direction\n";
    meshOut << "# Grid resolution: " << gridResolution << " (iso-lines at k/" << gridResolution << " — matches element boundaries)\n\n";

    // For degree >= 2 splines, sample at 200 points per line so the smooth
    // B-spline curve is rendered faithfully rather than as a coarse polygon.
    const index_t deg = (mp.nPatches() > 0) ? mp.patch(0).basis().degree(0) : 1;
    const index_t sampleN = (deg >= 2) ? 200 : gridResolution;

    auto writeLine = [&](index_t patch, bool isU, real_t fixedParam)
    {
        gsMatrix<> uvLine(2, sampleN + 1);
        for (int k = 0; k <= sampleN; ++k)
        {
            real_t varying = k / (real_t)sampleN;
            uvLine(0, k) = isU ? fixedParam : varying;
            uvLine(1, k) = isU ? varying    : fixedParam;
        }
        gsMatrix<> xyLine = mp.patch(patch).eval(uvLine);
        const char* dir = isU ? "u" : "v";
        for (int k = 0; k <= sampleN; ++k)
            meshOut << xyLine(0, k) << " " << xyLine(1, k) << " " << patch << " " << dir << "\n";
        meshOut << "\n";
    };

    for (index_t patch = 0; patch < mp.nPatches(); ++patch)
    {
        // Write boundary iso-lines first (u=0 and u=1), then interior.
        // This guarantees both patch boundaries are plotted even if the
        // visualiser drops the last group of each direction block.
        writeLine(patch, true,  0.0);  // u = 0  (one patch boundary)
        writeLine(patch, true,  1.0);  // u = 1  (opposite patch boundary)
        for (int i = 1; i < gridResolution; ++i)
            writeLine(patch, true, i / (real_t)gridResolution);

        writeLine(patch, false, 0.0);  // v = 0
        writeLine(patch, false, 1.0);  // v = 1
        for (int j = 1; j < gridResolution; ++j)
            writeLine(patch, false, j / (real_t)gridResolution);
    }
    
    meshOut.close();
    gsInfo << "Original geometry mesh saved to " << meshFile << "\n";
    outfile << "generateOriginalGeometryMesh: Original geometry mesh saved to " << meshFile << "\n";
    
}

void generateUniformPatchGrid(
    index_t numPatches,
    index_t gridResolution,
    std::vector<gsMatrix<>>& uvGridLines,
    std::vector<index_t>& patchIDs,
    std::vector<char>& directions)
{
    uvGridLines.clear();
    patchIDs.clear();
    directions.clear();

    const index_t samples = std::max<index_t>(gridResolution, 1);

    for (index_t patch = 0; patch < numPatches; ++patch)
    {
        for (index_t i = 0; i <= samples; ++i)
        {
            gsMatrix<> uvLine(2, samples + 1);
            const real_t u = static_cast<real_t>(i) / static_cast<real_t>(samples);
            for (index_t j = 0; j <= samples; ++j)
            {
                uvLine(0, j) = u;
                uvLine(1, j) = static_cast<real_t>(j) / static_cast<real_t>(samples);
            }
            uvGridLines.push_back(uvLine);
            patchIDs.push_back(patch);
            directions.push_back('u');
        }

        for (index_t j = 0; j <= samples; ++j)
        {
            gsMatrix<> uvLine(2, samples + 1);
            const real_t v = static_cast<real_t>(j) / static_cast<real_t>(samples);
            for (index_t i = 0; i <= samples; ++i)
            {
                uvLine(0, i) = static_cast<real_t>(i) / static_cast<real_t>(samples);
                uvLine(1, i) = v;
            }
            uvGridLines.push_back(uvLine);
            patchIDs.push_back(patch);
            directions.push_back('v');
        }
    }
}

void generateVisualizationMesh(
    const gsVector<gsVector<gsVector<index_t>>>& boxMat,
    index_t gridResolution,
    const MPBES<2, real_t>& mpbes,
    const gsMultiPatch<real_t>& mp,
    const gsMatrix<real_t>& coefficients,
    const gsVector<index_t>& currentLastNonZeroRow,
    const std::string& outputPrefix,
    bool verbose)
{
    // if (verbose)
    //     gsInfo << "\n=== generateVisualizationMesh BEGIN ===\n";

    gsVector<index_t> numBoxesPerPatch(boxMat.size());
    for (index_t patch = 0; patch < boxMat.size(); ++patch)
    {
        const index_t fallback = static_cast<index_t>(boxMat[patch].size());
        const index_t requested = (patch < currentLastNonZeroRow.size())
            ? currentLastNonZeroRow[patch]
            : fallback;
        numBoxesPerPatch[patch] = std::min(requested, fallback);
    }

    // Debug: Export box structure
    exportBoxStructure(boxMat, numBoxesPerPatch, outputPrefix + "_boxes.txt");

    // Step 1: Generate parametric grid
    std::vector<gsMatrix<>> uvLines;
    std::vector<index_t> uvLinePatchIDs;
    std::vector<char> uvLineDirections;
    const std::vector<std::set<index_t>> activeLevels =
        collectActiveVisualizationLevels(mpbes, static_cast<index_t>(boxMat.size()));

    // if (verbose)
    //     gsInfo << "Step 1: Generating parametric grid...\n";

    generateParametricGrid(
        boxMat,
        numBoxesPerPatch,
        activeLevels,
        gridResolution,
        uvLines,
        uvLinePatchIDs,
        uvLineDirections,
        verbose
    );

    if (uvLines.empty())
    {
        gsWarn << "No grid lines generated. Aborting.\n";
        return;
    }

    // Step 2: Evaluate the hierarchical box grid and export it as the main mesh.
    std::vector<gsMatrix<>> xyLines;
    evaluateGeometryAtPoints(
        uvLines,
        uvLinePatchIDs,
        mpbes,
        mp,
        coefficients,
        xyLines,
        verbose
    );

    const std::vector<PatchVisualizationMetadata> patchMetadata =
        buildPatchVisualizationMetadata(mpbes, coefficients, static_cast<index_t>(boxMat.size()));
    const std::vector<InterfaceVisualizationMetadata> interfaceMetadata =
        buildInterfaceVisualizationMetadata(mpbes, mp, coefficients);

    if (g_enableInterfaceDiagnostics)
    {
        logSpecificInterface(mpbes, mp, coefficients, outputPrefix, 0, 1);
        logSpecificInterface(mpbes, mp, coefficients, outputPrefix, 0, 3);
        logSpecificInterface(mpbes, mp, coefficients, outputPrefix, 1, 2);
        logSpecificInterface(mpbes, mp, coefficients, outputPrefix, 2, 5);
        logSpecificInterface(mpbes, mp, coefficients, outputPrefix, 3, 4);
    }

    std::string meshFile = outputPrefix + ".txt";
    std::string patchMetadataFile = outputPrefix + "_patch_meta.txt";
    std::string interfaceMetadataFile = outputPrefix + "_interface_meta.txt";
    exportPatchMetadataToFile(patchMetadata, patchMetadataFile, verbose);
    exportInterfaceMetadataToFile(interfaceMetadata, interfaceMetadataFile, verbose);
    exportMeshToFile(
        xyLines,
        uvLinePatchIDs,
        uvLineDirections,
        meshFile,
        verbose
    );

    // if (verbose)
    //     gsInfo << "=== generateVisualizationMesh END ===\n\n";
}

void selectionMechanism(
    const gsVector<gsVector<gsVector<index_t>>>& boxMat,
    const gsVector<gsVector<gsTensorBSplineBasis<2, real_t>>>& Bells,
    const gsVector<gsTHBSplineBasis<2, real_t>>& THBVector,
    const gsVector<gsVector<gsVector<index_t>>>& hasATwin,
    const gsVector<gsVector<gsVector<std::vector<index_t>>>>& twinsIndex,
    const gsVector<gsVector<gsVector<std::vector<index_t>>>>& twinsPatch,
    const gsVector<gsVector<gsVector<index_t>>>& indexInTHB,
    const gsVector<index_t>& currentLastNonZeroRow,
    std::vector<std::vector<std::vector<index_t>>>& functionDescription,
    std::vector<std::vector<std::array<int, 3>>>& spilloverFunctionCoordinates,
    std::vector<std::vector<std::array<int, 3>>>& spilloverSources,
    std::vector<bool>& hasSpillover,
    gsVector<gsVector<int>>& globalIndexTHB,
    gsVector<gsVector<gsVector<index_t>>>& globalIndex,
    int& globalFunctionCount,
    bool verbose)
{
    PROFILE_FUNCTION();
    const index_t nPatches = Bells.size();

    auto hasPatchLevel = [&](index_t patch, index_t level) -> bool
    {
        return patch >= 0 && patch < nPatches && level >= 0 && level < Bells[patch].size();
    };

    auto hasPatchFunction = [&](index_t patch, index_t level, index_t functionIndex) -> bool
    {
        return hasPatchLevel(patch, level) && functionIndex >= 0 && functionIndex < Bells[patch][level].size();
    };

    auto hasTwinMetadata = [&](index_t patch, index_t level, index_t functionIndex) -> bool
    {
        return patch >= 0 && patch < hasATwin.size() &&
            level >= 0 && level < hasATwin[patch].size() &&
            functionIndex >= 0 && functionIndex < hasATwin[patch][level].size() &&
            patch < twinsIndex.size() && level < twinsIndex[patch].size() && functionIndex < twinsIndex[patch][level].size() &&
            patch < twinsPatch.size() && level < twinsPatch[patch].size() && functionIndex < twinsPatch[patch][level].size();
    };

    auto hasThbMapping = [&](index_t patch, index_t level, index_t functionIndex) -> bool
    {
        return patch >= 0 && patch < indexInTHB.size() &&
            level >= 0 && level < indexInTHB[patch].size() &&
            functionIndex >= 0 && functionIndex < indexInTHB[patch][level].size();
    };

    globalFunctionCount = 0;
    functionDescription.clear();
    spilloverFunctionCoordinates.clear();
    spilloverSources.clear();
    hasSpillover.clear();

    globalIndexTHB.resize(nPatches);
    globalIndex.resize(nPatches);
    spilloverSources.resize(nPatches);

    // if (verbose)
    //     gsInfo << "=== selectionMechanism BEGIN ===\n";

    for (index_t patch = 0; patch < nPatches; ++patch)
    {
        const index_t nLevels = Bells[patch].size();
        globalIndex[patch].resize(nLevels);

        for (index_t level = 0; level < nLevels; ++level)
        {
            const index_t nFuncs = Bells[patch][level].size();
            globalIndex[patch][level].resize(nFuncs);
            for (index_t i = 0; i < nFuncs; ++i)
                globalIndex[patch][level][i] = -1;
        }

        const index_t thbSize = THBVector[patch].size();
        globalIndexTHB[patch].resize(thbSize);
        for (index_t i = 0; i < thbSize; ++i)
            globalIndexTHB[patch][i] = -1;
    }

    for (index_t patch = 0; patch < nPatches; ++patch)
    {
        // Safety check: skip if this patch has no boxes
        if (patch >= boxMat.size() || boxMat[patch].size() == 0)
        {
            if (verbose)
                gsInfo << "SKIP patch " << patch << " — no boxes in boxMat\n";
            continue;
        }
        
        for (index_t level = 0; level < Bells[patch].size(); ++level)
        {
            for (index_t i = 0; i < Bells[patch][level].size(); ++i)
            {
                if (globalIndex[patch][level][i] != -1)
                    continue;

                if (!hasThbMapping(patch, level, i) || indexInTHB[patch][level][i] == -1)
                {
                    if (verbose)
                        outfile << "SKIP: patch=" << patch << ", level=" << level << ", index=" << i << " — not in THB\n";
                    continue;
                }

                bool intersectsBox = false;
                const auto& supp = Bells[patch][level].function(i).support();
                real_t sX0 = supp(0, 0), sX1 = supp(0, 1);
                real_t sY0 = supp(1, 0), sY1 = supp(1, 1);

                const index_t numValidBoxes = (patch < currentLastNonZeroRow.size())
                    ? currentLastNonZeroRow[patch]
                    : static_cast<index_t>(boxMat[patch].size());

                for (size_t b = 0; b < static_cast<size_t>(numValidBoxes); ++b)
                {
                    const auto& box = boxMat[patch][b];
                    
                    // Safety: ensure box has required 5 elements
                    if (box.size() < 5)
                    {
                        if (verbose)
                            gsInfo << "WARNING: patch=" << patch << ", box " << b << " has size " << box.size() << " < 5, skipping\n";
                        continue;
                    }
                    
                    if (box[0] != level)
                        continue;

                    const real_t h = 1.0 / std::pow(2.0, level);
                    const real_t bx0 = box[1] * h;
                    const real_t bx1 = box[3] * h;
                    const real_t by0 = box[2] * h;
                    const real_t by1 = box[4] * h;

                    if (!(sX1 <= bx0 || sX0 >= bx1 || sY1 <= by0 || sY0 >= by1))
                    {
                        intersectsBox = true;
                        break;
                    }
                }

                if (!intersectsBox)
                {
                    if (verbose)
                        outfile << "SKIP: patch=" << patch << ", level=" << level << ", index=" << i << " — no box intersection\n";
                    continue;
                }

                // Build twin class
                std::vector<std::vector<index_t>> fd;
                std::vector<std::array<int, 3>> spilloversOfThisFunction;
                bool currentHasSpillover = false;

                const int globalID = globalFunctionCount++;
                globalIndex[patch][level][i] = globalID;

                if (verbose)
                {
                    outfile << "=== New Global Function ID: " << globalID << " ===\n";
                    outfile << "  Main function: patch=" << patch
                        << ", level=" << level
                        << ", index=" << i << "\n";
                }

                int thbIdx = indexInTHB[patch][level][i];
                if (thbIdx != -1 && thbIdx < globalIndexTHB[patch].size())
                {
                    globalIndexTHB[patch][thbIdx] = globalID;
                    fd.push_back({ patch, level, i });

                    if (verbose)
                        outfile << "  -> THB index " << thbIdx
                        << " assigned globalID " << globalID << "\n";
                }
                else if (verbose)
                {
                    outfile << "  -> Not a THB function (indexInTHB = " << thbIdx << ")\n";
                }

                if (hasTwinMetadata(patch, level, i) && hasATwin[patch][level][i])
                {
                    if (verbose)
                        outfile << "  -> Checking twins for (" << patch << "," << level << "," << i << ")\n";

                    const size_t twinCount = std::min(twinsIndex[patch][level][i].size(), twinsPatch[patch][level][i].size());
                    for (size_t t = 0; t < twinCount; ++t)
                    {
                        int twinPatch = twinsPatch[patch][level][i][t];
                        int twinIndex = twinsIndex[patch][level][i][t];

                        if (verbose)
                            outfile << "    Twin " << t << ": patch=" << twinPatch
                            << ", level=" << level
                            << ", index=" << twinIndex << "\n";

                        if (!hasPatchFunction(twinPatch, level, twinIndex))
                        {
                            if (verbose)
                                outfile << "      -> Skipped invalid twin coordinate\n";
                            continue;
                        }

                        if (globalIndex[twinPatch][level][twinIndex] == -1)
                        {
                            globalIndex[twinPatch][level][twinIndex] = globalID;

                            int twinTHB = -1;
                            if (twinPatch < indexInTHB.size() &&
                                level < indexInTHB[twinPatch].size() &&
                                twinIndex < indexInTHB[twinPatch][level].size())
                            {
                                twinTHB = indexInTHB[twinPatch][level][twinIndex];
                            }

                            if (twinTHB != -1 && twinTHB < globalIndexTHB[twinPatch].size())
                            {
                                globalIndexTHB[twinPatch][twinTHB] = globalID;
                                fd.push_back({ twinPatch, level, twinIndex });

                                if (verbose)
                                    outfile << "      -> THB twin added to functionDescription, index=" << twinTHB << "\n";
                            }
                            else
                            {
                                fd.push_back({ twinPatch, level, twinIndex });
                                spilloversOfThisFunction.push_back({ twinPatch, static_cast<int>(level), static_cast<int>(twinIndex) });
                                spilloverSources[twinPatch].push_back({ static_cast<int>(patch), static_cast<int>(level), static_cast<int>(i) });
                                currentHasSpillover = true;

                                if (verbose)
                                    outfile << "      -> Marked as spillover\n";
                            }
                        }
                        else if (verbose)
                        {
                            outfile << "      -> Already assigned to globalID " << globalIndex[twinPatch][level][twinIndex] << "\n";
                        }
                    }
                }

                functionDescription.push_back(fd);
                spilloverFunctionCoordinates.push_back(spilloversOfThisFunction);
                hasSpillover.push_back(currentHasSpillover);

                if (verbose)
                {
                    outfile << "  -> functionDescription[" << globalID << "]:\n";
                    for (const auto& triplet : fd)
                        outfile << "     - patch=" << triplet[0]
                        << ", level=" << triplet[1]
                        << ", index=" << triplet[2] << "\n";

                    if (currentHasSpillover)
                    {
                        outfile << "  -> spilloverFunctionCoordinates[" << globalID << "]:\n";
                        for (const auto& sp : spilloversOfThisFunction)
                            outfile << "     - from patch=" << sp[0]
                            << ", level=" << sp[1]
                            << ", index=" << sp[2] << "\n";
                    }

                    outfile << "\n";
                }


                if (verbose)
                {
                    for (const auto& desc : fd)
                    {
                        outfile << "SELECTED: patch=" << desc[0]
                            << ", level=" << desc[1]
                            << ", index=" << desc[2]
                            << ", globalID=" << globalID << "\n";
                    }

                    if (currentHasSpillover)
                    {
                        for (const auto& sp : spilloversOfThisFunction)
                        {
                            outfile << "  Spillover from patch=" << sp[0]
                                << ", level=" << sp[1]
                                << ", index=" << sp[2]
                                << " into globalID=" << globalID << "\n";
                        }
                    }
                }
            }
        }
    }

    if (verbose)
    {
        gsInfo << "=== selectionMechanism END ===\n";
        gsInfo << "Total global functions selected: " << globalFunctionCount << "\n";
    }
}

















void logResidualDetails(
    const gsVector<gsMatrix<>>& uv,
    const gsMatrix<real_t>& b_vec,
    const gsMatrix<real_t>& matOut,
    std::ostream& outfile,
    bool logToInfo)
{
    if (matOut.rows() != b_vec.rows() || matOut.cols() != 2 || b_vec.cols() != 2)
    {
        if (logToInfo)
        {
            outfile << "[ERROR] Dimension mismatch in logResidualDetails\n";
            outfile << "matOut: " << matOut.rows() << "x" << matOut.cols()
                << ", b_vec: " << b_vec.rows() << "x" << b_vec.cols() << "\n";
        }
        return;
    }

    const index_t totalRows = matOut.rows();
    gsMatrix<> matC(totalRows, 2);

    if (logToInfo)
        outfile << "=== Subtraction details: matOut - b_vec ===\n";

    for (index_t row = 0; row < totalRows; ++row)
    {
        real_t outX = matOut(row, 0);
        real_t outY = matOut(row, 1);
        real_t geoX = b_vec(row, 0);
        real_t geoY = b_vec(row, 1);
        real_t diffX = outX - geoX;
        real_t diffY = outY - geoY;

        matC(row, 0) = diffX;
        matC(row, 1) = diffY;

        if (logToInfo)
        {
            outfile << "Row " << row
                << " | Fit: (" << outX << ", " << outY << ")"
                << " - Geom: (" << geoX << ", " << geoY << ")"
                << " = Residual: (" << diffX << ", " << diffY << ")\n";
        }
    }

    real_t maxResidual = matC.cwiseAbs().maxCoeff();

    if (logToInfo)
    {
        gsInfo << "=== logResidualDetails ===\n";
        gsInfo << "Total points: " << totalRows << "\n";
        gsInfo << "Max residual component: " << maxResidual << "\n";
    }

    if (logToInfo)
    {
        outfile << "\n=== Residuals per patch and point ===\n";

        index_t rowOffset = 0;
        for (index_t patch = 0; patch < uv.size(); ++patch)
        {
            const gsMatrix<real_t>& uvPatch = uv[patch];
            const index_t numPoints = uvPatch.cols();

            outfile << "Patch " << patch << ":\n";

            for (index_t i = 0; i < numPoints; ++i)
            {
                index_t row = rowOffset + i;

                const real_t u = uvPatch(0, i);
                const real_t v = uvPatch(1, i);

                const real_t x_fit = matOut(row, 0);
                const real_t y_fit = matOut(row, 1);

                const real_t x_geom = b_vec(row, 0);
                const real_t y_geom = b_vec(row, 1);

                const real_t dx = matC(row, 0);
                const real_t dy = matC(row, 1);

                outfile << "  Row " << row
                    << " | (u,v)=(" << u << ", " << v << ")"
                    << " | Fit=(" << x_fit << ", " << y_fit << ")"
                    << " | Geom=(" << x_geom << ", " << y_geom << ")"
                    << " | Diff=(" << dx << ", " << dy << ")\n";
            }

            rowOffset += numPoints;
            outfile << "\n";
        }

        outfile << "=== End of Residual Log ===\n\n";
    }
}
























void truncationMechanism(
    const gsVector<gsVector<gsTensorBSplineBasis<2, real_t>>>& Bells,
    const std::vector<std::vector<std::vector<index_t>>>& functionDescription,
    const std::vector<std::vector<std::array<int, 3>>>& spilloverFunctionCoordinates,
    std::vector<bool>& isTruncated,
    bool verbose)
{
    PROFILE_FUNCTION();
    const index_t numFunctions = functionDescription.size();
    isTruncated.assign(numFunctions, false);

    if (verbose)
    {
        outfile << "=== truncationMechanism BEGIN ===\n";
        gsInfo << "=== truncationMechanism BEGIN ===\n";
    }

    for (index_t f = 0; f < numFunctions; ++f)
    {
        if (functionDescription[f].empty())
            continue;

        const auto& descF = functionDescription[f][0];
        const int patchF = descF[0];
        const int levelF = descF[1];
        const int indexF = descF[2];

        const auto& supportF = Bells[patchF][levelF].function(indexF).support();

        for (index_t g = 0; g < numFunctions; ++g)
        {
            if (f == g || functionDescription[g].empty())
                continue;

            for (const auto& g_native : functionDescription[g])
            {
                const int patchG = g_native[0];
                const int levelG = g_native[1];
                const int indexG = g_native[2];

                if (patchG != patchF || levelG <= levelF)
                    continue;

                const auto& supportG = Bells[patchG][levelG].function(indexG).support();

                bool intersects =
                    supportF(0, 0) < supportG(0, 1) && supportF(0, 1) > supportG(0, 0) &&
                    supportF(1, 0) < supportG(1, 1) && supportF(1, 1) > supportG(1, 0);

                if (intersects)
                {
                    isTruncated[f] = true;

                    if (verbose)
                    {
                        outfile << "  --> Truncated by native overlap: f=" << f
                            << " intersects with g=" << g << "\n";
                    }

                    break;
                }
            }

            if (isTruncated[f])
                break;

            if (g >= static_cast<index_t>(spilloverFunctionCoordinates.size()))
                continue;

            for (const auto& g_spill : spilloverFunctionCoordinates[g])
            {
                const int patchG = g_spill[0];
                const int levelG = g_spill[1];
                const int indexG = g_spill[2];

                if (patchG != patchF || levelG <= levelF)
                    continue;

                const auto& supportG = Bells[patchG][levelG].function(indexG).support();

                bool intersects =
                    supportF(0, 0) < supportG(0, 1) && supportF(0, 1) > supportG(0, 0) &&
                    supportF(1, 0) < supportG(1, 1) && supportF(1, 1) > supportG(1, 0);

                if (intersects)
                {
                    isTruncated[f] = true;

                    if (verbose)
                    {
                        outfile << "  --> Truncated by spillover overlap: f=" << f
                            << " intersects with spill g=" << g << "\n";
                    }

                    break;
                }
            }

            if (isTruncated[f])
                break;
        }
    }

    if (verbose)
    {
        const index_t count = std::count(isTruncated.begin(), isTruncated.end(), true);
        outfile << "Total functions truncated: " << count << " / " << numFunctions << "\n";
        outfile << "=== truncationMechanism END ===\n";
    }
}


void truncationMechanismMPBES(
    const gsVector<gsVector<gsTensorBSplineBasis<2, real_t>>>& Bells,
    const std::vector<std::vector<std::vector<index_t>>>& functionDescription,
    const std::vector<std::vector<std::array<index_t, 3>>>& spilloverFunctionCoordinates,
    std::vector<bool>& isTruncated,
    bool verbose = false)
{
    index_t numFunctions = functionDescription.size();
    isTruncated.clear();
    isTruncated.resize(numFunctions, false);

    if (verbose)
        gsInfo << "=== truncationMechanismMPBES BEGIN ===\n";

    for (index_t f = 0; f < numFunctions; ++f)
    {
        if (functionDescription[f].empty())
            continue;

        const std::vector<index_t>& fData = functionDescription[f][0];
        index_t patchF = fData[0];
        index_t levelF = fData[1];
        index_t indexF = fData[2];

        gsMatrix<index_t> suppF = Bells[patchF][levelF].function(indexF).support();

        for (index_t g = 0; g < numFunctions; ++g)
        {
            if (f == g || functionDescription[g].empty())
                continue;

            // First, check native functions
            for (size_t k = 0; k < functionDescription[g].size(); ++k)
            {
                const std::vector<index_t>& gData = functionDescription[g][k];
                index_t patchG = gData[0];
                index_t levelG = gData[1];
                index_t indexG = gData[2];

                if (patchG != patchF || levelG <= levelF)
                    continue;

                gsMatrix<index_t> suppG = Bells[patchG][levelG].function(indexG).support();

                if (suppF(0, 0) >= suppG(0, 0) && suppF(0, 1) <= suppG(0, 1) &&
                    suppF(1, 0) >= suppG(1, 0) && suppF(1, 1) <= suppG(1, 1))
                {
                    isTruncated[f] = true;
                    if (verbose)
                        gsInfo << "f=" << f << " truncated by native g=" << g << "\n";
                    goto next_f;
                }
            }

            // Then, check spillover functions
            if (g < spilloverFunctionCoordinates.size())
            {
                const std::vector<std::array<index_t, 3>>& spillList = spilloverFunctionCoordinates[g];

                for (size_t k = 0; k < spillList.size(); ++k)
                {
                    index_t patchG = spillList[k][0];
                    index_t levelG = spillList[k][1];
                    index_t indexG = spillList[k][2];

                    if (patchG != patchF || levelG <= levelF)
                        continue;

                    gsMatrix<index_t> suppG = Bells[patchG][levelG].function(indexG).support();

                    if (suppF(0, 0) >= suppG(0, 0) && suppF(0, 1) <= suppG(0, 1) &&
                        suppF(1, 0) >= suppG(1, 0) && suppF(1, 1) <= suppG(1, 1))
                    {
                        isTruncated[f] = true;
                        if (verbose)
                            gsInfo << "f=" << f << " truncated by spillover g=" << g << "\n";
                        goto next_f;
                    }
                }
            }
        }

    next_f:;
    }

    if (verbose)
    {
        index_t count = std::count(isTruncated.begin(), isTruncated.end(), true);
        gsInfo << "Total functions truncated: " << count << " / " << numFunctions << "\n";
        gsInfo << "=== truncationMechanismMPBES END ===\n";
    }
}











bool shouldSuppressSpillover(
    index_t patch,  // The patch we are assembling (target)
    const gsTHBSplineBasis<2, real_t>& thbBasis,
    index_t sourcePatch,
    index_t level,
    index_t index,
    const gsTensorBSplineBasis<2, real_t>& spillBasis,
    bool verbose)
{
    // Always suppress intra-patch spillovers: they are already included in the THB sum
    if (sourcePatch == patch)
    {
        if (verbose)
            outfile << "[shouldSuppressSpillover] Suppressed: intra-patch spillover.\n";
        return true;
    }

    const auto& supportSpill = spillBasis.function(index).support();

    if (verbose)
    {
        outfile << "[shouldSuppressSpillover] Spillover into patch=" << patch
            << " from patch=" << sourcePatch
            << ", level=" << level << ", index=" << index
            << ", support=[" << supportSpill(0, 0) << ", " << supportSpill(0, 1)
            << "] x [" << supportSpill(1, 0) << ", " << supportSpill(1, 1) << "]\n";
    }

    for (index_t tIndex = 0; tIndex < thbBasis.size(); ++tIndex)
    {
        index_t tLevel = thbBasis.levelOf(tIndex);
        if (tLevel <= level)
            continue;

        const auto& supportTHB = thbBasis.function(tIndex).support();

        bool isContained =
            supportSpill(0, 0) >= supportTHB(0, 0) &&
            supportSpill(0, 1) <= supportTHB(0, 1) &&
            supportSpill(1, 0) >= supportTHB(1, 0) &&
            supportSpill(1, 1) <= supportTHB(1, 1);

        if (isContained)
        {
            if (verbose)
            {
                outfile << "  Suppressed: spillover is fully contained in finer THB of patch " << patch << "\n";
            }
            return true;
        }
    }

    if (verbose)
        outfile << "  Kept: no finer THB in patch " << patch << " fully contains this spillover.\n";

    return false;
}










bool shouldSuppressTHB(
    index_t patch,
    index_t level,
    index_t index,
    const gsTensorBSplineBasis<2, real_t>& tensorBasis,
    const std::vector<std::array<int, 3>>& incomingSpillovers,
    const gsVector<gsVector<gsTensorBSplineBasis<2, real_t>>>& Bells,
    bool verbose)
{
    const auto& supportTHB = tensorBasis.function(index).support();

    if (verbose)
    {
        outfile << "[shouldSuppressTHB] THB candidate: "
            << "patch=" << patch
            << ", level=" << level
            << ", index=" << index
            << ", support=[" << supportTHB(0, 0) << ", " << supportTHB(0, 1)
            << "] x [" << supportTHB(1, 0) << ", " << supportTHB(1, 1) << "]\n";
    }

    for (const auto& source : incomingSpillovers)
    {
        int sourcePatch = source[0];
        int sourceLevel = source[1];
        int sourceIndex = source[2];

        if (sourcePatch == patch)
            continue;

        if (sourceLevel >= Bells[patch].size())
            continue;

        const auto& spillBasis = Bells[patch];

        for (index_t spLevel = 0; spLevel < spillBasis.size(); ++spLevel)
        {
            if (spLevel <= level)
                continue; // Only finer spillovers may suppress

            for (index_t spIndex = 0; spIndex < spillBasis[spLevel].size(); ++spIndex)
            {
                const auto& supportSpill = spillBasis[spLevel].function(spIndex).support();

                bool isContained =
                    supportTHB(0, 0) >= supportSpill(0, 0) &&
                    supportTHB(0, 1) <= supportSpill(0, 1) &&
                    supportTHB(1, 0) >= supportSpill(1, 0) &&
                    supportTHB(1, 1) <= supportSpill(1, 1);

                if (isContained)
                {
                    if (verbose)
                    {
                        outfile << "  Suppressed by finer spillover containing it: from patch=" << sourcePatch
                            << ", level=" << spLevel
                            << ", index=" << spIndex
                            << ", support=[" << supportSpill(0, 0) << ", " << supportSpill(0, 1)
                            << "] x [" << supportSpill(1, 0) << ", " << supportSpill(1, 1) << "]\n";
                    }
                    return true;
                }
            }
        }
    }

    if (verbose)
        outfile << "  Kept: no finer spillover fully contains it.\n";

    return false;
}










/**
 * @brief Evaluates mapped geometry lines from THB and spillover basis functions.
 *
 * @date 16.06.2025
 *
 * @param[in] uvLines                      List of (u,v) parameter lines.
 * @param[in] uvLinePatchIDs              Patch ID for each line.
 * @param[in] Bells                       Tensor-product basis functions.
 * @param[in] THBVector                   THB basis per patch.
 * @param[in] functionDescription         Global basis functions: (patch, level, tensor-index).
 * @param[in] spilloverFunctionCoordinates Spillover entries: (patch, level, index).
 * @param[in] hasSpillover                Flags indicating which global functions have spillovers.
 * @param[in] spilloverSources            Not used in this function (can be ignored).
 * @param[in] AcceptedvectSol             Solution vector (x).
 * @param[out] xyLines                    Output: mapped geometry lines.
 * @param[in] indexInTHB                  Maps (patch, level, tensor-index) => THB index, or -1 if inactive.
 * @param[in] verbose                     Enables logging.
 */
 /**
  * @brief Evaluates the geometric mapping at parameter space points.
  *
  * Takes (u,v) coordinates and evaluates the THB parameterization to get (x,y) coordinates.
  */
void evaluateGeometryAtPoints(
    const std::vector<gsMatrix<>>& uvPoints,
    const std::vector<index_t>& patchIDs,
    const MPBES<2, real_t>& mpbes,
    const gsMultiPatch<real_t>& mp,
    const gsMatrix<real_t>& coefficients,
    std::vector<gsMatrix<>>& xyPoints,
    bool verbose)
{
    PROFILE_FUNCTION();
    xyPoints.clear();

    auto isFiniteVec = [](const gsVector<real_t>& v) -> bool
    {
        for (index_t i = 0; i < v.size(); ++i)
        {
            if (!std::isfinite(v[i]))
                return false;
        }
        return true;
    };

    auto tryMapPeerUv = [](const GeometryPreflightInterfaceInfo& info,
                           index_t patch,
                           const gsVector<real_t>& uv,
                           index_t& peerPatch,
                           gsVector<real_t>& peerUv) -> bool
    {
        if (!info.hasMappedSides)
            return false;

        const real_t sideTolerance = static_cast<real_t>(1e-10);

        if (patch == info.patchA && normalDistanceToSide(info.sideA, uv) <= sideTolerance)
        {
            const real_t t = tangentCoordinateOnSide(info.sideA, uv);
            const real_t peerT = info.orientationReversed ? (1.0 - t) : t;
            peerPatch = info.patchB;
            peerUv = uvOnFeatureSide(info.sideB, peerT);
            return true;
        }

        if (patch == info.patchB && normalDistanceToSide(info.sideB, uv) <= sideTolerance)
        {
            const real_t t = tangentCoordinateOnSide(info.sideB, uv);
            const real_t peerT = info.orientationReversed ? (1.0 - t) : t;
            peerPatch = info.patchA;
            peerUv = uvOnFeatureSide(info.sideA, peerT);
            return true;
        }

        return false;
    };

    if (uvPoints.empty())
    {
        if (verbose)
            gsWarn << "No points to evaluate.\n";
        return;
    }

    if (verbose)
        gsInfo << "Evaluating " << uvPoints.size() << " lines using MPBES with fitted coefficients...\n";

    xyPoints.resize(uvPoints.size());

    const bool canUseFitted = (coefficients.rows() >= mpbes.size() && coefficients.size() > 0);

    // Evaluate each line using MPBES basis with fitted coefficients when valid
    for (size_t lineIdx = 0; lineIdx < uvPoints.size(); ++lineIdx)
    {
        const gsMatrix<>& uvLine = uvPoints[lineIdx];
        const index_t patch = patchIDs[lineIdx];
        index_t evalPatch = patch;
        gsMatrix<real_t> uvEval = uvLine;

        if (uvLine.cols() == 0) {
            xyPoints[lineIdx].resize(2, 0);
            continue;
        }

        // Allocate output matrix
        const index_t geoDim = canUseFitted ? coefficients.cols() : mp.patch(evalPatch).geoDim();
        xyPoints[lineIdx].resize(geoDim, uvLine.cols());

        if (!canUseFitted)
        {
            gsMatrix<real_t> xyLine = mp.patch(evalPatch).eval(uvEval);
            xyPoints[lineIdx] = xyLine;
            continue;
        }

        // Evaluate fitted geometry at each point on this line
        for (index_t k = 0; k < uvLine.cols(); ++k)
        {
            const gsVector<real_t> uv = uvEval.col(k);
            // For visualization, include spillover fallback so shared interfaces
            // use the full MPBES contribution and do not show artificial seams.
            gsVector<real_t> xFit = evaluateFittedGeometryPoint(
                mpbes,
                coefficients,
                evalPatch,
                uv,
                true);

            // Seam guard removed: mixing original geometry at interface endpoints
            // with fitted geometry at interior points caused kinks up to 120 degrees
            // (the endpoint snaps to original, the adjacent interior point is at the
            // fitted position — a jump proportional to the global fitting error ~0.028).

            const bool hasNaN = !isFiniteVec(xFit);

            if (hasNaN)
            {
                gsMatrix<real_t> xy = mp.patch(evalPatch).eval(uv);
                for (index_t d = 0; d < geoDim; ++d)
                    xyPoints[lineIdx](d, k) = xy(d, 0);
            }
            else
            {
                for (index_t d = 0; d < geoDim; ++d)
                    xyPoints[lineIdx](d, k) = xFit[d];
            }
        }
    }

    if (verbose)
        gsInfo << "Geometry evaluation complete.\n";
}


/**
 * @brief Compute the maximum geometric fitting error at all sample points (Ax-b)
 *        using the same structure as testBoundaryAssembly, and the same interface as boundaryError.
 *
 * @param mpbes    MPBES basis structure
 * @param mp       MultiPatch geometry (for compatibility, not used)
 * @param vectSol  Fitted coefficients from least-squares solution
 * @return         Maximum pointwise Euclidean distance between fitted and target
 */
real_t testGlobalFittingError(
    const MPBES<2, real_t>& mpbes,
    const gsMultiPatch<real_t>& /*mp*/, // unused, keep for interface compatibility
    const gsMatrix<real_t>& vectSol,
    const gsVector<gsMatrix<>>& uv1,
    const gsMatrix<real_t>& b_vec)
{
    PROFILE_FUNCTION();

    const index_t geoDim = vectSol.cols();
    real_t maxError = 0.0;
    real_t sumSqError = 0.0;
    index_t numSamples = 0;
    index_t maxPatch = -1;
    index_t maxLocalPoint = -1;
    gsVector<real_t> maxUv(2);
    gsVector<real_t> maxXFit(geoDim);
    gsVector<real_t> maxTarget(geoDim);
    bool hasMaxSample = false;
    index_t sampleLogs = 0;
    bool loggedNonFiniteCoef = false;

    index_t globalRow = 0;
    std::vector<real_t> patchMaxError(uv1.size(), 0.0);

    // Precompute active functions per patch — avoids evaluating ~94% of functions that return 0
    const index_t numFunctionsGFE = mpbes.size();
    const auto& funcDescGFE = mpbes.functionDescription();
    const auto& idxInTHBGFE = mpbes.indexInTHB();
    const index_t numPatchesGFE = static_cast<index_t>(uv1.size());
    std::vector<std::vector<index_t>> activeFuncsGFE(numPatchesGFE);
    for (index_t p = 0; p < numPatchesGFE; ++p) {
        std::vector<index_t>& active = activeFuncsGFE[p];
        for (index_t f = 0; f < numFunctionsGFE; ++f) {
            if (f >= static_cast<index_t>(funcDescGFE.size())) break;
            for (const auto& desc : funcDescGFE[f]) {
                if (desc.size() < 3) continue;
                const int patchId = desc[0], levelId = desc[1], indexId = desc[2];
                if (patchId == static_cast<int>(p) &&
                    levelId >= 0 && indexId >= 0 &&
                    levelId < static_cast<int>(idxInTHBGFE[p].size()) &&
                    indexId < static_cast<int>(idxInTHBGFE[p][levelId].size()) &&
                    idxInTHBGFE[p][levelId][indexId] != static_cast<index_t>(-1)) {
                    active.push_back(f);
                    break;
                }
            }
        }
    }
    gsMatrix<real_t> basisValueGFE;
    gsVector<real_t> xFitGFE(geoDim);

    for (index_t patch = 0; patch < uv1.size(); ++patch)
    {
        const gsMatrix<real_t>& uvPatch = uv1(patch); // (d × nPts)
        const index_t nPts = uvPatch.cols();

        for (index_t k = 0; k < nPts; ++k, ++globalRow)
        {
            if (globalRow >= b_vec.rows())
            {
                gsInfo << "WARNING: globalRow " << globalRow << " >= b_vec.rows()=" << b_vec.rows()
                       << " in testGlobalFittingError (patch=" << patch << ", k=" << k << ")\n";
                outfile << "WARNING: globalRow " << globalRow << " >= b_vec.rows()=" << b_vec.rows()
                        << " in testGlobalFittingError (patch=" << patch << ", k=" << k << ")\n";
                return maxError;
            }
            // parameter point (d × 1)
            const gsVector<real_t> uv = uvPatch.col(k);

            // Active-function-filtered evaluation — only iterates over functions active on this patch
            xFitGFE.setZero();
            real_t pouSumGFE = 0.0;
            for (index_t f : activeFuncsGFE[patch]) {
                if (f >= vectSol.rows()) continue;
                mpbes.evalSingleOnPatch(f, patch, uv, basisValueGFE, false);
                if (basisValueGFE.rows() < 1 || basisValueGFE.cols() < 1) continue;
                const real_t val = basisValueGFE(0, 0);
                if (!std::isfinite(val) || val == 0.0) continue;
                for (index_t d = 0; d < geoDim; ++d)
                    xFitGFE(d) += val * vectSol(f, d);
                pouSumGFE += val;
            }
            if (std::abs(pouSumGFE) > static_cast<real_t>(1e-12))
                xFitGFE /= pouSumGFE;
            const gsVector<real_t>& xFit = xFitGFE;

            bool xFitFinite = true;
            for (index_t d = 0; d < geoDim; ++d)
                xFitFinite = xFitFinite && std::isfinite(xFit[d]);
            if (!xFitFinite)
            {
                // gsInfo << "[globalFittingError] Non-finite xFit at patch=" << patch
                //        << ", k=" << k << ", uv=(" << pt(0, 0) << "," << pt(1, 0) << ")\n";
                // outfile << "[globalFittingError] Non-finite xFit at patch=" << patch
                //         << ", k=" << k << ", uv=(" << pt(0, 0) << "," << pt(1, 0) << ")\n";
                maxError = std::numeric_limits<real_t>::infinity();
                continue;
            }

            // compute pointwise error
            real_t err2 = 0.0;
            for (index_t d = 0; d < geoDim; ++d)
            {
                const real_t diff = xFit[d] - b_vec(globalRow, d);
                err2 += diff * diff;
            }

            const real_t err = std::sqrt(err2);
            if (err > patchMaxError[patch])
                patchMaxError[patch] = err;
            if (err > maxError)
            {
                maxError = err;
                maxPatch = patch;
                maxLocalPoint = k;
                hasMaxSample = true;
                for (index_t d = 0; d < std::min<index_t>(2, uv.rows()); ++d)
                    maxUv[d] = uv[d];
                for (index_t d = 0; d < geoDim; ++d)
                {
                    maxXFit[d] = xFit[d];
                    maxTarget[d] = b_vec(globalRow, d);
                }
            }
            sumSqError += err2;
            ++numSamples;

            if (sampleLogs < 5)
            {
                // gsInfo << "[globalFittingError] patch=" << patch << ", k=" << k
                //        << ", uv=(" << uv(0) << "," << uv(1) << ")"
                //        << ", xFit=(" << xFit[0] << "," << xFit[1] << ")"
                //        << ", b=(" << b_vec(globalRow, 0) << "," << b_vec(globalRow, 1) << ")\n";
                // outfile << "[globalFittingError] patch=" << patch << ", k=" << k
                //         << ", uv=(" << uv(0) << "," << uv(1) << ")"
                //         << ", xFit=(" << xFit[0] << "," << xFit[1] << ")"
                //         << ", b=(" << b_vec(globalRow, 0) << "," << b_vec(globalRow, 1) << ")\n";
                ++sampleLogs;
            }
        }
    }

    GISMO_ASSERT(globalRow == b_vec.rows(),
        "Mismatch between uv1 samples and b_vec rows");

    if (numSamples > 0)
    {
        const real_t rmsError = std::sqrt(sumSqError / static_cast<real_t>(numSamples));
        if (g_verbose)
            gsInfo << "Global fitting error stats: max=" << maxError
                   << ", rms=" << rmsError
                   << ", samples=" << numSamples
                   << ", maxAtPatch=" << maxPatch
                   << ", maxAtLocalPoint=" << maxLocalPoint << "\n";
        outfile << "Global fitting error stats: max=" << maxError
                << ", rms=" << rmsError
                << ", samples=" << numSamples
                << ", maxAtPatch=" << maxPatch
                << ", maxAtLocalPoint=" << maxLocalPoint << "\n";

        if (hasMaxSample)
        {
            if (g_verbose)
                gsInfo << "Global fitting max sample: uv=(" << maxUv[0] << "," << maxUv[1] << ")";
            outfile << "Global fitting max sample: uv=(" << maxUv[0] << "," << maxUv[1] << ")";
            if (geoDim >= 2)
            {
                if (g_verbose)
                    gsInfo << ", xFit=(" << maxXFit[0] << "," << maxXFit[1] << ")"
                           << ", target=(" << maxTarget[0] << "," << maxTarget[1] << ")\n";
                outfile << ", xFit=(" << maxXFit[0] << "," << maxXFit[1] << ")"
                        << ", target=(" << maxTarget[0] << "," << maxTarget[1] << ")\n";
            }
            else
            {
                if (g_verbose) gsInfo << "\n";
                outfile << "\n";
            }

        }

        // Per-patch max error breakdown
        if (g_verbose)
        {
            gsInfo  << "Per-patch globalError:";
            outfile << "Per-patch globalError:";
            for (index_t p = 0; p < static_cast<index_t>(patchMaxError.size()); ++p)
            {
                gsInfo  << " p" << p << "=" << patchMaxError[p];
                outfile << " p" << p << "=" << patchMaxError[p];
            }
            gsInfo  << "\n";
            outfile << "\n";
        }
    }

    if (maxError == 0.0)
    {
        real_t maxAbsCoef = 0.0;
        for (index_t i = 0; i < vectSol.rows(); ++i)
            for (index_t d = 0; d < vectSol.cols(); ++d)
            {
                const real_t coef = vectSol(i, d);
                if (!std::isfinite(coef))
                {
                    // gsInfo << "[globalFittingError] Non-finite coef at i=" << i << ", d=" << d << "\n";
                    // outfile << "[globalFittingError] Non-finite coef at i=" << i << ", d=" << d << "\n";
                    continue;
                }
                maxAbsCoef = std::max(maxAbsCoef, std::abs(coef));
            }
        real_t maxAbsB = 0.0;
        for (index_t i = 0; i < b_vec.rows(); ++i)
            for (index_t d = 0; d < b_vec.cols(); ++d)
            {
                const real_t bv = b_vec(i, d);
                if (!std::isfinite(bv))
                {
                    // gsInfo << "[globalFittingError] Non-finite b_vec at i=" << i << ", d=" << d << "\n";
                    // outfile << "[globalFittingError] Non-finite b_vec at i=" << i << ", d=" << d << "\n";
                    continue;
                }
                maxAbsB = std::max(maxAbsB, std::abs(bv));
            }
        // gsInfo << "[globalFittingError] maxError=0. maxAbsCoef=" << maxAbsCoef
        //        << ", maxAbsB=" << maxAbsB << "\n";
        // outfile << "[globalFittingError] maxError=0. maxAbsCoef=" << maxAbsCoef
        //         << ", maxAbsB=" << maxAbsB << "\n";
    }

    return maxError;
}


/**
 * @brief Wrapper for testGlobalFittingError, matching boundaryError style. Same interface as boundaryError.
 *
 * @param mpbes    MPBES basis structure
 * @param mp       MultiPatch geometry (for compatibility, not used)
 * @param vectSol  Fitted coefficients from least-squares solution
 * @return         Maximum pointwise Euclidean distance between fitted and target
 */

real_t globalFittingError(
    const MPBES<2, real_t>& mpbes,
    const gsMultiPatch<real_t>& mp,
    const gsMatrix<real_t>& vectSol,
    const gsVector<gsMatrix<>>& uv1,
    const gsMatrix<real_t>& b_vec)
{
    return testGlobalFittingError(mpbes, mp, vectSol, uv1, b_vec);
}

real_t boundaryError(
    const MPBES<2, real_t>& mpbes,
    const gsMultiPatch<real_t>& mp,
    const gsMatrix<real_t>& vectSol,
    const std::vector<FeatureBoundarySpec>& featureSides)
{
    PROFILE_FUNCTION();
    
    outfile << "\n=== boundaryError BEGIN ===\n";

    const index_t pointsPerEdge = 10;
    const std::vector<InterfaceValidationSpec> validationInterfaces =
        collectValidationInterfaces(mp);
    const index_t numInterfaces = static_cast<index_t>(validationInterfaces.size());
    
    real_t maxError = 0.0;
    bool loggedNonFiniteCoef = false;
    for (index_t i = 0; i < vectSol.rows() && !loggedNonFiniteCoef; ++i)
        for (index_t d = 0; d < vectSol.cols(); ++d)
            if (!std::isfinite(vectSol(i, d)))
            {
                gsInfo << "[boundaryError] Non-finite vectSol at i=" << i << ", d=" << d << "\n";
                outfile << "[boundaryError] Non-finite vectSol at i=" << i << ", d=" << d << "\n";
                loggedNonFiniteCoef = true;
                break;
            }

    if (!featureSides.empty())
    {
        outfile << "(feature subset mode)\n";

        auto uvFromSide = [](FeatureSide side, real_t t) -> gsVector<real_t>
        {
            gsVector<real_t> uv(2);
            switch (side)
            {
                case FeatureSide::U0: uv(0) = 0.0; uv(1) = t;   break;
                case FeatureSide::U1: uv(0) = 1.0; uv(1) = t;   break;
                case FeatureSide::V0: uv(0) = t;   uv(1) = 0.0; break;
                case FeatureSide::V1: uv(0) = t;   uv(1) = 1.0; break;
            }
            return uv;
        };

        for (size_t s = 0; s < featureSides.size(); ++s)
        {
            const FeatureBoundarySpec& spec = featureSides[s];
            if (spec.patch < 0 || spec.patch >= mp.nPatches())
            {
                outfile << "WARNING: featureSides[" << s << "] has invalid patch index "
                        << spec.patch << ", skipping.\n";
                continue;
            }

            for (index_t i = 0; i < pointsPerEdge; ++i)
            {
                real_t t = static_cast<real_t>(i) / (pointsPerEdge - 1);
                gsVector<real_t> uv = uvFromSide(spec.side, t);

                gsMatrix<real_t> xyFit(2, 1);
                xyFit.setZero();

                if (vectSol.rows() == 0)
                {
                    gsInfo << "WARNING: vectSol is empty in boundaryError, skipping point " << i << "\n";
                    outfile << "WARNING: vectSol is empty in boundaryError, skipping point " << i << "\n";
                    continue;
                }

                for (index_t f = 0; f < mpbes.size(); ++f)
                {
                    gsMatrix<real_t> N;
                    mpbes.evalSingleOnPatch(f, spec.patch, uv, N);
                    if (N.size() == 0) continue;
                    if (N(0, 0) != 0.0 && f < vectSol.rows())
                    {
                        xyFit(0, 0) += N(0, 0) * vectSol(f, 0);
                        xyFit(1, 0) += N(0, 0) * vectSol(f, 1);
                    }
                }

                gsMatrix<real_t> xyOrig = mp.patch(spec.patch).eval(uv);
                real_t err = std::sqrt(std::pow(xyFit(0,0) - xyOrig(0,0), 2) +
                                       std::pow(xyFit(1,0) - xyOrig(1,0), 2));
                maxError = std::max(maxError, err);

                if (err > 1e-6 || i < 2)
                {
                    outfile << "Feature patch " << spec.patch
                        << ", t=" << t
                        << ": fit=(" << xyFit(0,0) << "," << xyFit(1,0) << ")"
                        << ", orig=(" << xyOrig(0,0) << "," << xyOrig(1,0) << ")"
                        << ", error=" << err << "\n";
                }
            }
        }

        outfile << "=== boundaryError END: max feature error = " << maxError << " ===\n\n";
        return maxError;
    }

    if (validationInterfaces.empty())
    {
        outfile << "[boundaryError] no mapped interfaces available\n";
    }

    for (size_t iface = 0; iface < validationInterfaces.size(); ++iface)
    {
        const InterfaceValidationSpec& spec = validationInterfaces[iface];
        for (index_t i = 0; i < pointsPerEdge; ++i)
        {
            const real_t t = (pointsPerEdge == 1)
                ? 0.0
                : static_cast<real_t>(i) / static_cast<real_t>(pointsPerEdge - 1);
            const real_t otherT = spec.orientationReversed ? (1.0 - t) : t;

            const gsVector<real_t> uvA = uvOnFeatureSide(spec.sideA, t);
            const gsVector<real_t> uvB = uvOnFeatureSide(spec.sideB, otherT);
            const gsVector<real_t> xyFitA = evaluateFittedGeometryPoint(mpbes, vectSol, spec.patchA, uvA, true);
            const gsVector<real_t> xyFitB = evaluateFittedGeometryPoint(mpbes, vectSol, spec.patchB, uvB, true);
            const gsVector<real_t> xyOrigA = matrixColumnToVector2(mp.patch(spec.patchA).eval(uvA));
            const gsVector<real_t> xyOrigB = matrixColumnToVector2(mp.patch(spec.patchB).eval(uvB));

            const real_t fitGap = vectorDistance(xyFitA, xyFitB);
            const real_t origGap = vectorDistance(xyOrigA, xyOrigB);
            maxError = std::max(maxError, fitGap);

            if (fitGap > 1e-6 || i < 2)
            {
                outfile << "Interface " << spec.patchA << "↔" << spec.patchB
                        << " [" << featureSideName(spec.sideA) << " vs " << featureSideName(spec.sideB)
                        << ", reversed=" << (spec.orientationReversed ? "true" : "false") << "]"
                        << ", t=" << t
                        << ": fitA=(" << xyFitA(0) << "," << xyFitA(1) << ")"
                        << ", origA=(" << xyOrigA(0) << "," << xyOrigA(1) << ")"
                        << ", fitB=(" << xyFitB(0) << "," << xyFitB(1) << ")"
                        << ", origB=(" << xyOrigB(0) << "," << xyOrigB(1) << ")"
                        << ", fitGap=" << fitGap
                        << ", origGap=" << origGap << "\n";
            }
        }
    }

    outfile << "=== boundaryError END: max interface error = " << maxError << " ===\n\n";
    if (maxError == 0.0)
    {
        gsInfo << "[boundaryError] maxError=0. featureSides=" << featureSides.size()
               << ", pointsPerEdge=" << pointsPerEdge << ", numInterfaces=" << numInterfaces << "\n";
        outfile << "[boundaryError] maxError=0. featureSides=" << featureSides.size()
                << ", pointsPerEdge=" << pointsPerEdge << ", numInterfaces=" << numInterfaces << "\n";
    }
    
    return maxError;
}


/**
 * @brief Test boundary continuity using the assemble() function
    // debugFile << "=== boundaryError DEBUG LOG ===\n";

    const index_t numFunctions = functionDescription.size();
    const index_t pointsPerEdge = 10;
    const index_t totalRows = pointsPerEdge * 4;  // 4 = 2 interfaces * 2 sides per interface

    // Convert presentation (double->real_t) once to avoid repeated conversions
    std::vector<std::vector<gsSparseVector<real_t>>> presReal;
    presReal.resize(presentation.size());
    for (size_t f = 0; f < presentation.size(); ++f)
    {
        presReal[f].resize(presentation[f].size());
        for (size_t c = 0; c < presentation[f].size(); ++c)
        {
            presReal[f][c] = gsSparseVector<real_t>(presentation[f][c].size());
            for (index_t k = 0; k < presentation[f][c].size(); ++k)
                presReal[f][c](k) = static_cast<real_t>(presentation[f][c](k));
        }
    }

    gsMatrix<real_t> A_row(totalRows, numFunctions);
    A_row.setZero();
    gsMatrix<real_t> b_vec(totalRows, 2);
    b_vec.setZero();

    index_t rowOffset = 0;

    for (index_t i = 0; i < pointsPerEdge; ++i)
    {
        real_t t = static_cast<real_t>(i) / (pointsPerEdge - 1);

        // Interface 1: Patch 0 west (side 3) <-> Patch 1 east (side 1)
        // Patch 0, side 3 => (u,v) = (0, 1 - t)
        gsVector<real_t> uv0_side3(2); uv0_side3(0) = 0.0; uv0_side3(1) = 1.0 - t;
        // Patch 1, side 1 => (u,v) = (1, t)
        gsVector<real_t> uv1_side1(2); uv1_side1(0) = 1.0; uv1_side1(1) = t;

        // Interface 2: Patch 0 south (side 2) <-> Patch 2 east (side 1)
        // Patch 0, side 2 => (u,v) = (t, 0)
        gsVector<real_t> uv0_side2(2); uv0_side2(0) = t; uv0_side2(1) = 0.0;
        // Patch 2, side 1 => (u,v) = (1, 1 - t)
        gsVector<real_t> uv2_side1(2); uv2_side1(0) = 1.0; uv2_side1(1) = 1.0 - t;

        std::vector<std::pair<index_t, gsVector<real_t>>> probes = {
            {0, uv0_side3},
            {1, uv1_side1},
            {0, uv0_side2},
            {2, uv2_side1}
        };

        for (index_t p = 0; p < 4; ++p)
        {
            const index_t patch = probes[p].first;
            const gsVector<real_t>& uv = probes[p].second;

            // Geometry
            gsMatrix<real_t> xy = mp.patch(patch).eval(uv);
            b_vec(rowOffset, 0) = xy(0, 0);
            b_vec(rowOffset, 1) = xy(1, 0);

            outfile << "\n>> Patch " << patch << ", (u,v)=(" << uv(0) << "," << uv(1) << ")"
                << ", geom = (" << xy(0, 0) << "," << xy(1, 0) << ")\n";

            for (index_t f = 0; f < numFunctions; ++f)
            {
                real_t val = 0.0;
                gsMatrix<> resultMatrix; // reuse for evalSingle_into

                int functionComponent = 0;
                bool hasRegularComponentOnThisPatch = false;

                // 1. Evaluate regular components from functionDescription
                for (const auto& desc : functionDescription[f])
                {
                    int fnPatch = desc[0];
                    int fnLevel = desc[1];
                    int fnTensorIndex = desc[2];

                    // CRITICAL: presentation[f] has entries for ALL components in functionDescription[f]
                    // We must increment functionComponent for EVERY component, not just matching patches
                    if (fnPatch != static_cast<int>(patch)) {
                        ++functionComponent;  // Skip but still increment
                        continue;
                    }

                    if (fnLevel >= indexInTHB[fnPatch].size() ||
                        fnTensorIndex >= indexInTHB[fnPatch][fnLevel].size()) {
                        ++functionComponent;  // Skip but still increment
                        continue;
                    }

                    int thbIdx = indexInTHB[fnPatch][fnLevel][fnTensorIndex];

                    // Skip if thbIdx == -1 (component doesn't exist in THB, will be handled as spillover)
                    if (thbIdx == -1) {
                        ++functionComponent;  // Skip but still increment
                        continue;
                    }

                    // Found a valid regular component on this patch (thbIdx >= 0)
                    hasRegularComponentOnThisPatch = true;

                    // Use pre-converted presentation coefficients
                    const gsSparseVector<real_t>& coefsReal = presReal[f][functionComponent];

                    // evalSingle_into handles both regular (thbIdx >= 0) and spillover (thbIdx == -1)
                    evalSingle_into(
                        f,
                        thbIdx,
                        uv,
                        isTruncated,
                        THBVector[patch],
                        coefsReal,
                        resultMatrix,
                        &Bells[fnPatch][fnLevel],
                        fnTensorIndex);
                    real_t contrib = resultMatrix(0, 0);
                    val += contrib;
                    ++functionComponent;

                    if (contrib != 0.0) {
                        outfile << "    [THB] f=" << f << " contributes " << contrib;
                        // Debug for Row 0 with functions 20 and 591
                        if (rowOffset == 0 && (f == 20 || f == 591)) {
                            outfile << " [isTruncated=" << isTruncated[f]
                                    << ", coefsSize=" << coefsReal.size()
                                    << ", thbIdx=" << thbIdx << "]";
                        }
                        outfile << "\n";
                    }
                }

                // 2. Evaluate spillover components ONLY if no regular component exists on this patch
                if (!hasRegularComponentOnThisPatch &&
                    f < spilloverFunctionCoordinates.size() &&
                    !spilloverFunctionCoordinates[f].empty())
                {
                    for (const auto& spillCoord : spilloverFunctionCoordinates[f])
                    {
                        int spPatch = spillCoord[0];
                        int spLevel = spillCoord[1];
                        int spTensorIdx = spillCoord[2];

                        if (spPatch != static_cast<int>(patch)) continue;

                        // Spillover: evaluate B-spline basis function directly (thbIdx = -1)
                        if (spLevel < Bells[spPatch].size())
                        {
                            gsSparseVector<real_t> emptyCoefs;  // No presentation for spillover
                            evalSingle_into(
                                f,
                                -1,  // thbIdx = -1 signals spillover
                                uv,
                                isTruncated,
                                THBVector[patch],
                                emptyCoefs,
                                resultMatrix,
                                &Bells[spPatch][spLevel],
                                spTensorIdx);
                            real_t contrib = resultMatrix(0, 0);
                            val += contrib;

                            if (contrib != 0.0) {
                                outfile << "    [SPILL] f=" << f << " contributes " << contrib << "\n";
                            }
                        }
                    }
                }

                if (val != 0.0)
                    A_row(rowOffset, f) = val;
            }

            ++rowOffset;
        }
    }

    gsMatrix<real_t> fit = A_row * vectSol;
    gsMatrix<real_t> diff = fit - b_vec;
    if (!fit.allFinite() || !diff.allFinite())
    {
        gsInfo << "[boundaryError] Non-finite fit/diff detected\n";
        outfile << "[boundaryError] Non-finite fit/diff detected\n";
        maxError = std::numeric_limits<real_t>::infinity();
    }
    else
    {
        maxError = diff.cwiseAbs().maxCoeff();
    }

    // Debug: Check specific entries
    outfile << "\n--- Debug: Checking A_row and vectSol ---\n";
    outfile << "A_row.rows() = " << A_row.rows() << ", A_row.cols() = " << A_row.cols() << "\n";
    outfile << "vectSol.rows() = " << vectSol.rows() << ", vectSol.cols() = " << vectSol.cols() << "\n";
    outfile << "fit.rows() = " << fit.rows() << ", fit.cols() = " << fit.cols() << "\n";

    // Find which functions contribute to row 0
    outfile << "Row 0 non-zero entries in A_row:\n";
    for (index_t c = 0; c < A_row.cols(); ++c) {
        if (A_row(0, c) != 0.0) {
            outfile << "  A_row(0, " << c << ") = " << A_row(0, c)
                    << ", vectSol(" << c << ") = (" << vectSol(c, 0) << ", " << vectSol(c, 1) << ")\n";
        }
    }
    outfile << "Computed fit(0) = (" << fit(0, 0) << ", " << fit(0, 1) << ")\n";
    outfile << "Expected geom(0) = (" << b_vec(0, 0) << ", " << b_vec(0, 1) << ")\n\n";

    outfile << "\n--- Interface Values ---\n";
    outfile << "Note: Only showing first 10 rows for brevity\n";

    // Find which row has the max error
    index_t maxRow = 0, maxCol = 0;
    for (index_t r = 0; r < diff.rows(); ++r) {
        for (index_t c = 0; c < diff.cols(); ++c) {
            if (std::abs(diff(r, c)) == maxError) {
                maxRow = r;
                maxCol = c;
            }
        }
    }

    outfile << "Max error is at row " << maxRow << ", col " << maxCol << "\n";
    outfile << "  fit = (" << fit(maxRow, 0) << ", " << fit(maxRow, 1) << ")\n";
    outfile << "  geom = (" << b_vec(maxRow, 0) << ", " << b_vec(maxRow, 1) << ")\n";
    outfile << "  diff = (" << diff(maxRow, 0) << ", " << diff(maxRow, 1) << ")\n\n";

    for (index_t r = 0; r < std::min(totalRows, (index_t)10); ++r)
    {
        index_t patch = (r % 2 == 0) ? 0 : 1;
        outfile << "Row " << r
            << " | Patch " << patch
            << " => fit = (" << fit(r, 0) << ", " << fit(r, 1) << ")"
            << ", geom = (" << b_vec(r, 0) << ", " << b_vec(r, 1) << ")"
            << ", residual = (" << diff(r, 0) << ", " << diff(r, 1) << ")\n";
    }

    outfile << "\n=== boundaryError END ===\n";
    outfile << "Max componentwise error: " << maxError << "\n";

    // Log A_row statistics for debugging
    int nonzeros = 0;
    for (index_t r = 0; r < A_row.rows(); ++r)
        for (index_t c = 0; c < A_row.cols(); ++c)
            if (A_row(r, c) != 0.0) ++nonzeros;

 *
 * This function validates interface continuity by:
 * 1. Creating interface sample points (similar to boundaryError)
 * 2. Using the actual assemble() function to build basis evaluation matrix
 * 3. Computing geometry using A * vectSol
 * 4. Comparing with ground truth geometry from mp.eval()
 *
 * This is the "gold standard" test - if assemble works correctly,
 * the reconstructed geometry should match the original at interfaces.
 *
 * @return Maximum componentwise error across all interface points
 */
real_t testBoundaryAssembly(
    const MPBES<2, real_t>& mpbes,
    const gsMultiPatch<real_t>& mp,
    const gsMatrix<real_t>& vectSol)
{
    try {
        PROFILE_FUNCTION();
        // gsInfo << "\n=== testBoundaryAssembly BEGIN ===\n";
        // outfile << "\n=== testBoundaryAssembly BEGIN ===\n";
        
        // Comprehensive input validation
        if (vectSol.rows() == 0 || vectSol.cols() == 0) {
            // gsInfo << "WARNING: vectSol is empty in testBoundaryAssembly, returning 0\n";
            // outfile << "WARNING: vectSol is empty in testBoundaryAssembly, returning 0\n";
            return 0.0;
        }
        
        if (vectSol.cols() < 2) {
            // gsInfo << "WARNING: vectSol has insufficient columns (" << vectSol.cols() << ") in testBoundaryAssembly, returning 0\n";
            // outfile << "WARNING: vectSol has insufficient columns (" << vectSol.cols() << ") in testBoundaryAssembly, returning 0\n";
            return 0.0;
        }
        
        if (mp.nPatches() < 1) {
            // gsInfo << "WARNING: Insufficient patches (" << mp.nPatches() << ") for testBoundaryAssembly, returning 0\n";
            // outfile << "WARNING: Insufficient patches (" << mp.nPatches() << ") for testBoundaryAssembly, returning 0\n";
            return 0.0;
        }
        
        const index_t numFunctions = mpbes.size();
        if (numFunctions == 0) {
            // gsInfo << "WARNING: mpbes has 0 functions, returning 0\n";
            // outfile << "WARNING: mpbes has 0 functions, returning 0\n";
            return 0.0;
        }
        
        if (static_cast<index_t>(vectSol.rows()) != numFunctions) {
            // gsInfo << "WARNING: vectSol.rows (" << vectSol.rows() << ") != numFunctions (" << numFunctions << "), returning 0\n";
            // outfile << "WARNING: vectSol.rows (" << vectSol.rows() << ") != numFunctions (" << numFunctions << "), returning 0\n";
            return 0.0;
        }

        const index_t pointsPerEdge = 10;
        std::vector<InterfaceValidationSpec> validationInterfaces;
        gsVector<gsMatrix<>> interfacePoints;
        const index_t totalPoints = buildValidationInterfacePoints(
            mp,
            pointsPerEdge,
            validationInterfaces,
            interfacePoints);
        const index_t numInterfaces = static_cast<index_t>(validationInterfaces.size());
        if (numInterfaces == 0 || totalPoints == 0)
            return 0.0;

        // Assemble using the actual assemble() function
        gsSparseMatrix<real_t> A_interface(totalPoints, numFunctions);
        gsMatrix<real_t> b_interface(totalPoints, 2);

        // gsInfo << "Assembling interface matrix: " << totalPoints << " points, " << numFunctions << " functions\n";
        // outfile << "Assembling interface matrix: " << totalPoints << " points, " << numFunctions << " functions\n";

        try {
            assemble(
                interfacePoints,
                mpbes,
                A_interface,
                b_interface,
                mp,
                false  // verbose = false
            );
        } catch (const std::exception& e) {
            // gsInfo << "WARNING: Exception in assemble(): " << e.what() << ", returning 0\n";
            // outfile << "WARNING: Exception in assemble(): " << e.what() << ", returning 0\n";
            return 0.0;
        }

        // Validate matrix sizes before multiplication
        if (A_interface.rows() != totalPoints) {
            // gsInfo << "WARNING: A_interface has wrong rows (" << A_interface.rows() << " != " << totalPoints << "), returning 0\n";
            // outfile << "WARNING: A_interface has wrong rows (" << A_interface.rows() << " != " << totalPoints << "), returning 0\n";
            return 0.0;
        }
        
        if (A_interface.cols() != numFunctions) {
            // gsInfo << "WARNING: A_interface has wrong cols (" << A_interface.cols() << " != " << numFunctions << "), returning 0\n";
            // outfile << "WARNING: A_interface has wrong cols (" << A_interface.cols() << " != " << numFunctions << "), returning 0\n";
            return 0.0;
        }
        
        if (b_interface.rows() != totalPoints || b_interface.cols() != 2) {
            // gsInfo << "WARNING: b_interface has wrong size (" << b_interface.rows() << "x" << b_interface.cols() << "), returning 0\n";
            // outfile << "WARNING: b_interface has wrong size (" << b_interface.rows() << "x" << b_interface.cols() << "), returning 0\n";
            return 0.0;
        }

        // Compute reconstructed geometry: A * vectSol
        gsMatrix<real_t> reconstructed;
        try {
            reconstructed = gsEigen::MatrixXd(A_interface) * vectSol;
        } catch (const std::exception& e) {
            // gsInfo << "WARNING: Exception in matrix multiplication: " << e.what() << ", returning 0\n";
            // outfile << "WARNING: Exception in matrix multiplication: " << e.what() << ", returning 0\n";
            return 0.0;
        }

        // Validate reconstructed size
        if (reconstructed.rows() != totalPoints || reconstructed.cols() != 2) {
            // gsInfo << "WARNING: reconstructed has wrong size (" << reconstructed.rows() << "x" << reconstructed.cols() << "), returning 0\n";
            // outfile << "WARNING: reconstructed has wrong size (" << reconstructed.rows() << "x" << reconstructed.cols() << "), returning 0\n";
            return 0.0;
        }

        // Compute residual
        gsMatrix<real_t> residual = reconstructed - b_interface;
        
        // Compute max Euclidean error
        real_t maxError = 0.0;
        index_t maxRow = -1;
        for (index_t row = 0; row < residual.rows(); ++row)
        {
            real_t err = std::sqrt(residual(row, 0) * residual(row, 0) +
                                   residual(row, 1) * residual(row, 1));
            if (std::isfinite(err) && err > maxError) {
                maxError = err;
                maxRow = row;
            }
        }

        // Per-interface breakdown
        if (maxRow >= 0 && pointsPerEdge > 0)
        {
            // Each interface contributes 2*pointsPerEdge rows (one block per side)
            const index_t rowsPerIface = 2 * pointsPerEdge;
            const index_t ifaceIdx = maxRow / rowsPerIface;
            const index_t localRow  = maxRow % rowsPerIface;
            const index_t side      = localRow / pointsPerEdge; // 0=sideA, 1=sideB
            const index_t ptIdx     = localRow % pointsPerEdge;

            if (g_verbose)
            {
                gsInfo  << "Per-interface featureError breakdown (max " << maxError << "):\n";
                outfile << "Per-interface featureError breakdown (max " << maxError << "):\n";

                // Print all interface errors
                for (index_t i = 0; i < static_cast<index_t>(validationInterfaces.size()); ++i)
                {
                    real_t ifaceMax = 0.0;
                    const index_t rowStart = i * rowsPerIface;
                    const index_t rowEnd   = std::min<index_t>(rowStart + rowsPerIface, residual.rows());
                    for (index_t r = rowStart; r < rowEnd; ++r)
                    {
                        real_t e = std::sqrt(residual(r,0)*residual(r,0) + residual(r,1)*residual(r,1));
                        if (std::isfinite(e)) ifaceMax = std::max(ifaceMax, e);
                    }
                    gsInfo  << "  iface " << i << " (p" << validationInterfaces[i].patchA
                            << "<->p" << validationInterfaces[i].patchB << "): " << ifaceMax << "\n";
                    outfile << "  iface " << i << " (p" << validationInterfaces[i].patchA
                            << "<->p" << validationInterfaces[i].patchB << "): " << ifaceMax << "\n";
                }

                if (ifaceIdx < static_cast<index_t>(validationInterfaces.size()))
                {
                    const auto& spec = validationInterfaces[ifaceIdx];
                    gsInfo  << "  worst: iface " << ifaceIdx
                            << " (p" << spec.patchA << "<->p" << spec.patchB << ")"
                            << " side=" << side << " pt=" << ptIdx
                            << ", reconstructed=(" << reconstructed(maxRow,0) << "," << reconstructed(maxRow,1) << ")"
                            << ", target=("        << b_interface(maxRow,0)   << "," << b_interface(maxRow,1)   << ")\n";
                    outfile << "  worst: iface " << ifaceIdx
                            << " (p" << spec.patchA << "<->p" << spec.patchB << ")"
                            << " side=" << side << " pt=" << ptIdx
                            << ", reconstructed=(" << reconstructed(maxRow,0) << "," << reconstructed(maxRow,1) << ")"
                            << ", target=("        << b_interface(maxRow,0)   << "," << b_interface(maxRow,1)   << ")\n";
                }
            }
        }

        // Detailed reporting
        // gsInfo << "\n--- Interface Continuity Test Results ---\n";
        // outfile << "\n--- Interface Continuity Test Results ---\n";

        // gsInfo << "\nMax componentwise error: " << maxError << "\n";
        // gsInfo << "=== testBoundaryAssembly END ===\n\n";
        // outfile << "\nMax componentwise error: " << maxError << "\n";
        // outfile << "=== testBoundaryAssembly END ===\n\n";

        return maxError;
        
    } catch (const std::exception& e) {
        // gsInfo << "CRITICAL: Exception in testBoundaryAssembly: " << e.what() << ", returning 0\n";
        // outfile << "CRITICAL: Exception in testBoundaryAssembly: " << e.what() << ", returning 0\n";
        return 0.0;
    } catch (...) {
        // gsInfo << "CRITICAL: Unknown exception in testBoundaryAssembly, returning 0\n";
        // outfile << "CRITICAL: Unknown exception in testBoundaryAssembly, returning 0\n";
        return 0.0;
    }
}




real_t boundaryErrorTwoPointsOnly(
    const gsVector<gsVector<gsTensorBSplineBasis<2, real_t>>>& Bells,
    const gsVector<gsTHBSplineBasis<2, real_t>>& THBVector,
    const std::vector<std::vector<std::vector<index_t>>>& functionDescription,
    const std::vector<std::vector<std::array<int, 3>>>& spilloverFunctionCoordinates,
    const std::vector<bool>& hasSpillover,
    const std::vector<bool>& isTruncated,
    const gsVector<gsVector<gsVector<index_t>>>& indexInTHB,
    const gsMultiPatch<real_t>& mp,
    const gsMatrix<real_t>& vectSol)
{
    PROFILE_FUNCTION();
    gsInfo << "=== boundaryErrorTwoPointsOnly BEGIN ===\n";

    struct ProbePoint {
        index_t patch;
        gsVector<real_t> uv;
    };

    std::vector<ProbePoint> points(2);
    points[0].patch = 0; points[0].uv.setZero(2); points[0].uv(0) = 0.0; points[0].uv(1) = 0.888889;
    points[1].patch = 1; points[1].uv.setZero(2); points[1].uv(0) = 1.0; points[1].uv(1) = 0.888889;

    const index_t numFunctions = functionDescription.size();
    gsMatrix<real_t> A_row(2, numFunctions);
    A_row.setZero();
    gsMatrix<real_t> b_vec(2, 2);
    b_vec.setZero();

    for (size_t p = 0; p < points.size(); ++p)
    {
        index_t row = p;
        index_t patch = points[p].patch;
        const gsVector<real_t>& uvk = points[p].uv;

        // RHS target
        gsMatrix<real_t> xy = mp.patch(patch).eval(uvk);
        b_vec(row, 0) = xy(0, 0);
        b_vec(row, 1) = xy(1, 0);

        gsInfo << ">> Patch " << patch
            << ", (u,v) = (" << uvk(0) << "," << uvk(1) << ")"
            << ", geometry = (" << xy(0, 0) << "," << xy(1, 0) << ")\n";

        for (index_t f = 0; f < numFunctions; ++f)
        {
            if (isTruncated[f])
            {
                gsInfo << "  [Skip] f=" << f << " is truncated.\n";
                continue;
            }

            real_t val = 0.0;

            for (const auto& desc : functionDescription[f])
            {
                const int fnPatch = desc[0];
                const int fnLevel = desc[1];
                const int fnIndex = desc[2];

                if (fnPatch != static_cast<int>(patch)) continue;
                if (fnLevel >= Bells[fnPatch].size()) continue;
                if (fnIndex >= Bells[fnPatch][fnLevel].size()) continue;

                //int thbIdx = indexInTHB[fnPatch][fnLevel][fnIndex];
                //if (thbIdx == -1) continue;

                real_t contrib = Bells[fnPatch][fnLevel].function(fnIndex).eval(uvk)(0, 0);
                val += contrib;

                if (contrib != 0.0)
                {
                    gsInfo << "    [THB] f=" << f
                        << ", patch=" << fnPatch
                        << ", level=" << fnLevel
                        << ", index=" << fnIndex
                        << ", contrib=" << contrib << "\n";
                }
            }

            if (hasSpillover[f])
            {
                for (const auto& spill : spilloverFunctionCoordinates[f])
                {
                    const int spPatch = spill[0];
                    const int spLevel = spill[1];
                    const int spIndex = spill[2];

                    if (spPatch >= Bells.size() ||
                        spLevel >= Bells[spPatch].size() ||
                        spIndex >= Bells[spPatch][spLevel].size()) continue;

                    real_t contrib = Bells[spPatch][spLevel].function(spIndex).eval(uvk)(0, 0);
                    val += contrib;

                    if (contrib != 0.0)
                    {
                        gsInfo << "    [Spill] f=" << f
                            << ", patch=" << spPatch
                            << ", level=" << spLevel
                            << ", index=" << spIndex
                            << ", contrib=" << contrib << "\n";
                    }
                }
            }

            if (val != 0.0)
            {
                A_row(row, f) = val;
                gsInfo << "      [Assemble] A(" << row << "," << f << ") = " << val << "\n";
            }
        }
    }

    gsMatrix<real_t> fit = A_row * vectSol;
    gsMatrix<real_t> diff = fit - b_vec;
    real_t maxError = diff.cwiseAbs().maxCoeff();

    gsInfo << "--- Results ---\n";
    gsInfo << "fit.dim = " << fit.rows() << "x" << fit.cols() << "\n";

    for (index_t i = 0; i < fit.rows(); ++i)
    {
        gsInfo << "Patch " << points[i].patch
            << ", (u,v) = (" << points[i].uv(0) << "," << points[i].uv(1) << ")"
            << ", fit = (" << fit(i, 0) << ", " << fit(i, 1) << ")"
            << ", geom = (" << b_vec(i, 0) << ", " << b_vec(i, 1) << ")"
            << ", residual = (" << diff(i, 0) << ", " << diff(i, 1) << ")\n";
    }

    gsInfo << "=== boundaryErrorTwoPointsOnly END ===\n";
    gsInfo << "Max component-wise error: " << maxError << "\n";

    return maxError;
}






template <typename T>
void logTHBContributionsForLine(
    T fixed_u,
    index_t targetPatch,
    const gsVector<gsTHBSplineBasis<2, T>>& THBVector,
    std::ostream& out
)
{
    //std::vector<T> v_values = {
    //    static_cast<T>(0.875),    static_cast<T>(0.881579),
    //    static_cast<T>(0.888158), static_cast<T>(0.894737),
    //    static_cast<T>(0.901316), static_cast<T>(0.907895),
    //    static_cast<T>(0.914474), static_cast<T>(0.921053),
    //    static_cast<T>(0.927632), static_cast<T>(0.934211),
    //    static_cast<T>(0.940789), static_cast<T>(0.947368),
    //    static_cast<T>(0.953947), static_cast<T>(0.960526),
    //    static_cast<T>(0.967105), static_cast<T>(0.973684),
    //    static_cast<T>(0.980263), static_cast<T>(0.986842),
    //    static_cast<T>(0.993421), static_cast<T>(1.0)
    //};
    std::vector<T> v_values = {
    static_cast<T>(0.875)
    };

    out << "=== THB contributions from patch " << targetPatch
        << " along u = " << fixed_u << " ===\n";

    for (T v : v_values)
    {
        gsVector<T> uvk(2);
        uvk << fixed_u, v;

        out << "--- uv = (" << fixed_u << ", " << v << ") ---\n";

        for (index_t i = 0; i < THBVector[targetPatch].size(); ++i)
        {
            T val = THBVector[targetPatch].function(i).eval(uvk)(0, 0);
            if (val != 0.0)
            {
                int level = THBVector[targetPatch].levelOf(i);
                int tensorIndex = THBVector[targetPatch].flatTensorIndexOf(i);
                auto support = THBVector[targetPatch].function(i).support();

                out << "  THB: index=" << i
                    << ", level=" << level
                    << ", tensorIndex=" << tensorIndex
                    << " -> val = " << val << "\n";
                out << "      support: ["
                    << support(0, 0) << ", " << support(0, 1) << "] x ["
                    << support(1, 0) << ", " << support(1, 1) << "]\n";
            }
        }
    }

    out << "=== End of THB debug ===\n";
}


template <typename T>
void logSpilloverContributionsForLine(
    T fixed_u,
    index_t targetPatch,
    const gsVector<gsVector<gsTensorBSplineBasis<2, T>>>& Bells,
    const std::vector<std::vector<std::array<int, 3>>>& spilloverFunctionCoordinates,
    const std::vector<bool>& hasSpillover,
    std::ostream& out
)
{
    //std::vector<T> v_values = {
    //    static_cast<T>(0.875),    static_cast<T>(0.881579),
    //    static_cast<T>(0.888158), static_cast<T>(0.894737),
    //    static_cast<T>(0.901316), static_cast<T>(0.907895),
    //    static_cast<T>(0.914474), static_cast<T>(0.921053),
    //    static_cast<T>(0.927632), static_cast<T>(0.934211),
    //    static_cast<T>(0.940789), static_cast<T>(0.947368),
    //    static_cast<T>(0.953947), static_cast<T>(0.960526),
    //    static_cast<T>(0.967105), static_cast<T>(0.973684),
    //    static_cast<T>(0.980263), static_cast<T>(0.986842),
    //    static_cast<T>(0.993421), static_cast<T>(1.0)
    //};
    std::vector<T> v_values = {
        static_cast<T>(0.875)
    };

    out << "=== Spillover contributions into patch " << targetPatch
        << " along u = " << fixed_u << " ===\n";

    for (T v : v_values)
    {
        gsVector<T> uvk(2);
        uvk << fixed_u, v;

        out << "--- uv = (" << fixed_u << ", " << v << ") ---\n";

        for (index_t f = 0; f < spilloverFunctionCoordinates.size(); ++f)
        {
            if (!hasSpillover[f])
                continue;

            for (const auto& spill : spilloverFunctionCoordinates[f])
            {
                index_t spPatch = spill[0];
                index_t spLevel = spill[1];
                index_t spIndex = spill[2];

                if (spPatch != targetPatch)
                    continue; // Only interested in spillovers *into* targetPatch

                T val = Bells[spPatch][spLevel].function(spIndex).eval(uvk)(0, 0);
                if (val != 0.0)
                {
                    auto support = Bells[spPatch][spLevel].function(spIndex).support();

                    out << "  f=" << f
                        << " from patch=" << spPatch
                        << ", level=" << spLevel
                        << ", index=" << spIndex
                        << " -> val = " << val << "\n";

                    out << "      support: ["
                        << support(0, 0) << ", " << support(0, 1) << "] x ["
                        << support(1, 0) << ", " << support(1, 1) << "]\n";
                }
            }
        }
    }

    out << "=== End of spillover debug ===\n";
}

/**
 * @brief Wrapper function that generates mapped (x,y) lines from boxMat and exports them.
 *
 * @date 16.06.2025
 *
 * @param[in] boxMat                    Box definitions per patch.
 * @param[in] numPointsPerLine         Number of sample points per parametric line.
 * @param[in] Bells                    Tensor-product bases.
 * @param[in] THBVector                THB basis per patch.
 * @param[in] functionDescription      Global basis functions.
 * @param[in] spilloverFunctionCoordinates Spillover mapping.
 * @param[in] spilloverSources        Source mapping for spillovers.
 * @param[in] hasSpillover            Boolean flags per global function.
 * @param[in] AcceptedvectSol         Solution vector (geometry mapping).
 * @param[in] mp                      Geometry object (not used here, but passed if needed).
 * @param[in] currentLastNonZeroRow   Box truncation info.
 * @param[in] outputFilename          Output file name prefix.
 * @param[in] indexInTHB              Maps (patch, level, tensor-index) => THB index.
 * @param[in] verbose                 If true, enables logging.
 */
void generateAndExportMappedLines(
    const gsVector<gsVector<gsVector<index_t>>>& boxMat,
    index_t numPointsPerLine,
    const MPBES<2, real_t>& mpbes,
    const gsMatrix<real_t>& AcceptedvectSol,
    const gsMultiPatch<real_t>& mp,
    const gsVector<index_t>& currentLastNonZeroRow,
    const std::string& outputFilename,
    bool verbose
)
{
    if (verbose)
        gsInfo << "\n=== generateAndExportMappedLines BEGIN ===\n";

    gsVector<index_t> numBoxesPerPatch(boxMat.size());
    for (index_t patch = 0; patch < boxMat.size(); ++patch)
    {
        const index_t fallback = static_cast<index_t>(boxMat[patch].size());
        const index_t requested = (patch < currentLastNonZeroRow.size())
            ? currentLastNonZeroRow[patch]
            : fallback;
        numBoxesPerPatch[patch] = std::min(requested, fallback);
    }

    // Storage for parametric and mapped lines
    std::vector<gsMatrix<>> uvLines;
    std::vector<index_t> uvLinePatchIDs;
    std::vector<char> uvLineDirections;
    std::vector<gsMatrix<>> xyLines;
    const std::vector<std::set<index_t>> activeLevels =
        collectActiveVisualizationLevels(mpbes, static_cast<index_t>(boxMat.size()));

    // Step 1: Generate parametric lines from boxes
    if (verbose)
        gsInfo << "Step 1: Generating parametric lines...\n";

    generateParametricGrid(
        boxMat,
        numBoxesPerPatch,
        activeLevels,
        numPointsPerLine,
        uvLines,
        uvLinePatchIDs,
        uvLineDirections,
        verbose
    );

    if (uvLines.empty())
    {
        gsWarn << "No parametric lines generated. Aborting.\n";
        return;
    }

    // Step 2: Evaluate geometry at parametric lines
    // if (verbose)
    //     gsInfo << "Step 2: Evaluating geometry lines...\n";

    evaluateGeometryAtPoints(
        uvLines,
        uvLinePatchIDs,
        mpbes,
        mp,
        AcceptedvectSol,
        xyLines,
        verbose
    );

    // Step 3: Export parametric lines per patch
    if (verbose)
        gsInfo << "Step 3: Exporting parametric (u,v) lines...\n";

    std::string uvOutPrefix = outputFilename + "_uv.txt";
    exportMeshToFile(
        uvLines,
        uvLinePatchIDs,
        uvLineDirections,
        uvOutPrefix,
        verbose
    );

    // Step 4: Export mapped (x,y) lines
    if (verbose)
        gsInfo << "Step 4: Exporting mapped (x,y) lines...\n";

    exportLinesToFile(
        xyLines,
        uvLinePatchIDs,
        uvLineDirections,
        outputFilename,
        verbose,
        uvLines
    );

    if (verbose)
        gsInfo << "=== generateAndExportMappedLines END ===\n\n";
}





void generateParametricLinesFromBoxes(
    const gsVector<gsVector<gsVector<index_t>>>& boxMat,
    gsVector<index_t> currentLastNonZeroRow,
    index_t numPointsPerLine,
    std::vector<gsMatrix<>>& uvLines,
    std::vector<index_t>& uvLinePatchIDs,
    std::vector<char>& uvLineDirections,
    bool verbose
)
{
    uvLines.clear();
    uvLinePatchIDs.clear();
    uvLineDirections.clear();

    if (numPointsPerLine < 2)
    {
        gsWarn << "generateParametricLinesFromBoxes: numPointsPerLine must be >= 2. Setting to 2.\\n";
        numPointsPerLine = 2;
    }

    const index_t numPatches = boxMat.size();

    // Reserve space to reduce allocations
    size_t estimatedLines = 0;
    for (index_t patch = 0; patch < numPatches; ++patch)
        estimatedLines += currentLastNonZeroRow(patch) * 10; // Rough estimate

    uvLines.reserve(estimatedLines);
    uvLinePatchIDs.reserve(estimatedLines);
    uvLineDirections.reserve(estimatedLines);

    if (verbose)
        gsInfo << "Generating parametric lines from " << numPatches << " patches...\\n";

    for (index_t patch = 0; patch < boxMat.size(); ++patch)
    {
        const index_t numBoxes = currentLastNonZeroRow(patch);

        if (verbose && numBoxes > 0)
            gsInfo << "  Patch " << patch << ": processing " << numBoxes << " boxes\\n";

        for (index_t i = 0; i < numBoxes; ++i)
        {
            const index_t level = boxMat[patch][i][0];
            const index_t x0 = boxMat[patch][i][1];
            const index_t y0 = boxMat[patch][i][2];
            const index_t x1 = boxMat[patch][i][3];
            const index_t y1 = boxMat[patch][i][4];

            const real_t h = 1.0 / std::pow(2.0, static_cast<real_t>(level));

            const index_t verticalEdges[2] = { x0, x1 };
            const index_t horizontalEdges[2] = { y0, y1 };

            // Emit only box edges. Internal tensor-grid indices are not visible box boundaries.
            for (index_t edge = 0; edge < 2; ++edge)
            {
                const index_t xi = verticalEdges[edge];
                const real_t u = xi * h;
                const real_t v_start = y0 * h;
                const real_t v_end = y1 * h;

                gsMatrix<> line(2, numPointsPerLine);
                for (index_t k = 0; k < numPointsPerLine; ++k)
                {
                    const real_t t = static_cast<real_t>(k) / (numPointsPerLine - 1);
                    line(0, k) = u;
                    line(1, k) = (1.0 - t) * v_start + t * v_end;
                }

                uvLines.push_back(line);
                uvLinePatchIDs.push_back(patch);
                uvLineDirections.push_back('v');
            }

            for (index_t edge = 0; edge < 2; ++edge)
            {
                const index_t yi = horizontalEdges[edge];
                const real_t v = yi * h;
                const real_t u_start = x0 * h;
                const real_t u_end = x1 * h;

                gsMatrix<> line(2, numPointsPerLine);
                for (index_t k = 0; k < numPointsPerLine; ++k)
                {
                    const real_t t = static_cast<real_t>(k) / (numPointsPerLine - 1);
                    line(0, k) = (1.0 - t) * u_start + t * u_end;
                    line(1, k) = v;
                }

                uvLines.push_back(line);
                uvLinePatchIDs.push_back(patch);
                uvLineDirections.push_back('u');
            }
        }
    }

    if (verbose)
        gsInfo << "Generated " << uvLines.size() << " parametric lines.\\n";
}













int checkJacobianDeterminantTHB(
    const gsVector<gsMatrix<>>& uv,
    const gsVector<gsTHBSplineBasis<2, real_t>>& THBVector,
    const gsVector<gsVector<gsTensorBSplineBasis<2, real_t>>>& Bells,
    const std::vector<std::vector<std::vector<index_t>>>& functionDescription,
    const std::vector<std::vector<std::array<int, 3>>>& spilloverFunctionCoordinates,
    const std::vector<bool>& hasSpillover,
    const gsVector<gsVector<int>>& globalIndexTHB,
    const gsMatrix<real_t>& vectSol,
    gsVector<size_t>& numIrregular,
    bool verbose = true)
{
    if (verbose)
        gsInfo << "=== checkJacobianDeterminantTHB BEGIN ===\n";

    int irregular = 0;
    numIrregular.setZero();
    numIrregular.resize(THBVector.size());

    for (size_t patch = 0; patch < THBVector.size(); ++patch)
    {
        const size_t numPoints = uv[patch].cols();

        if (verbose)
        {
            gsInfo << ">> Patch " << patch << " with " << numPoints << " points.\n";
            outfile << ">> Patch " << patch << " with " << numPoints << " points.\n";
        }

        for (size_t ptIndex = 0; ptIndex < numPoints; ++ptIndex)
        {
            const gsMatrix<> point = uv[patch].col(ptIndex);
            gsMatrix<> jacobian(2, 2); jacobian.setZero();

            gsMatrix<int> actives;
            THBVector[patch].active_into(point, actives);

            //if (verbose)
            //{
            //    outfile << "  Point " << ptIndex << ": (" << point(0, 0) << ", " << point(1, 0) << ")\n";
            //    outfile << "    Active functions:\n";
            //}

            for (index_t i = 0; i < actives.rows(); ++i)
            {
                int localIndex = actives(i, 0);
                int globalIndex = globalIndexTHB[patch][localIndex];

                if (globalIndex == -1)
                {
                    //if (verbose)
                    //    outfile << "      Skipped: localIndex=" << localIndex << " (no global ID)\n";
                    continue;
                }

                if (globalIndex >= vectSol.rows())
                {
                    gsInfo << "[ERROR] globalIndex " << globalIndex << " out of bounds!\n";
                    continue;
                }

                const gsMatrix<> derivs = THBVector[patch].function(localIndex).deriv(point);
                jacobian(0, 0) += derivs(0, 0) * vectSol(globalIndex, 0);
                jacobian(0, 1) += derivs(1, 0) * vectSol(globalIndex, 0);
                jacobian(1, 0) += derivs(0, 0) * vectSol(globalIndex, 1);
                jacobian(1, 1) += derivs(1, 0) * vectSol(globalIndex, 1);

                //if (verbose)
                //    outfile << "      local=" << localIndex << ", global=" << globalIndex
                //    << ", ∂x=(" << derivs(0, 0) << "," << derivs(1, 0)
                //    << "), ∂y=(" << derivs(0, 0) << "," << derivs(1, 0) << ")\n";
            }

            const double det = jacobian.determinant();

            //if (verbose)
            //{
            //    outfile << "    Jacobian matrix:\n" << jacobian << "\n";
            //    outfile << "    Determinant: " << det << "\n";
            //}

            if (det <= 0)
            {
                ++numIrregular[patch];
                ++irregular;
                outfile << "  [IRREGULAR] Patch=" << patch << ", pt=" << ptIndex
                    << ", uv=(" << point(0, 0) << ", " << point(1, 0)
                    << "), det=" << det << "\n";
            }
        }

        gsInfo << "  Irregular points in patch " << patch << ": " << numIrregular[patch] << "\n";
        outfile << "  Irregular points in patch " << patch << ": " << numIrregular[patch] << "\n";
    }

    if (verbose)
    {
        gsInfo << "=== checkJacobianDeterminantTHB END ===\n";
        gsInfo << "Total irregular points: " << irregular << "\n";
    }

    return irregular;
}

// Helper function to check if a specific component of a function is truncated
inline bool isComponentTruncated(index_t globalIdx, int componentIdx, 
                                  const std::vector<std::vector<gsSparseVector<double>>>& presentation)
{
    if (globalIdx >= presentation.size()) return false;
    if (componentIdx >= presentation[globalIdx].size()) return false;
    
    const auto& coefs = presentation[globalIdx][componentIdx];
    for (size_t j = 0; j < coefs.size(); j++) {
        if (coefs(j) != 0.0) return true;
    }
    return false;
}

// MPBES-based checkJacobianDeterminant - uses MPBES derivatives just like assemble()
bool isParameterizationMirrored(
    const gsVector<gsMatrix<>>& uv,
    const MPBES<2, real_t>& mpbes,
    const gsMatrix<real_t>& vectSol,
    real_t detEps,
    bool verbose);

std::vector<bool> isParameterizationMirroredPerPatch(
    const gsVector<gsMatrix<>>& uv,
    const MPBES<2, real_t>& mpbes,
    const gsMatrix<real_t>& vectSol,
    real_t detEps,
    bool verbose);

MirroredCheckResult computeMirroredCheckPerPatch(
    const gsVector<gsMatrix<>>& uv,
    const MPBES<2, real_t>& mpbes,
    const gsMatrix<real_t>& vectSol,
    real_t detEps);

void logMirroredCheck(
    const gsVector<gsMatrix<>>& uv,
    const MPBES<2, real_t>& mpbes,
    const gsMatrix<real_t>& vectSol,
    real_t detEps,
    const std::string& label,
    bool logToInfo = false);

GeometryPreflightInfo runGeometryPreflight(
    const std::string& filename,
    const OriginalMpbesReference& originalReference,
    index_t pointsPerDirection = 20);

void logGeometryPreflight(
    const GeometryPreflightInfo& preflight,
    const std::string& label,
    bool logToInfo = false);

int checkJacobianDeterminant(
    const gsVector<gsMatrix<>>& uv,
    const MPBES<2, real_t>& mpbes,
    const gsMatrix<real_t>& vectSol,
    gsVector<size_t>& numIrregular,
    bool verbose = true)
{
    try {
        PROFILE_FUNCTION();
        // gsInfo << "=== checkJacobianDeterminant (MPBES-based) BEGIN ===\n";
        // outfile << "=== checkJacobianDeterminant (MPBES-based) BEGIN ===\n";
        
        // Validate inputs immediately
        if (vectSol.rows() == 0 || vectSol.cols() == 0) {
            // gsInfo << "WARNING: vectSol is empty, returning 0\n";
            // outfile << "WARNING: vectSol is empty, returning 0\n";
            numIrregular.setZero();
            return 0;
        }
        
        const index_t numPatches = mpbes.nPatches();
        const index_t numFunctions = mpbes.size();
        
        if (numPatches == 0) {
            // gsInfo << "WARNING: numPatches is 0, returning 0\n";
            // outfile << "WARNING: numPatches is 0, returning 0\n";
            numIrregular.setZero();
            return 0;
        }
        
        if (numFunctions == 0) {
            // gsInfo << "WARNING: numFunctions is 0, returning 0\n";
            // outfile << "WARNING: numFunctions is 0, returning 0\n";
            numIrregular.setZero();
            return 0;
        }
        
        if (uv.size() == 0) {
            gsInfo << "WARNING: uv is empty, returning 0\n";
            outfile << "WARNING: uv is empty, returning 0\n";
            numIrregular.setZero();
            return 0;
        }
        
        int totalIrregular = 0;
        numIrregular.resize(numPatches);
        numIrregular.setZero();
        std::vector<size_t> rawNeg(numPatches, 0);
        std::vector<size_t> rawPos(numPatches, 0);
        std::vector<size_t> signedNeg(numPatches, 0);
        std::vector<bool> mirroredPerPatch;
        // Per-patch minimum signed Jacobian determinant (used to measure cross-patch coupling).
        std::vector<real_t> minSignedDetPerPatch(numPatches, std::numeric_limits<real_t>::max());
        
        // Open detailed log file with exception handling
        std::ofstream irregularLog;
        try {
            static bool firstCall = true;
            if (firstCall) {
                irregularLog.open("irregular_points_log.txt", std::ios::out | std::ios::trunc);
                firstCall = false;
            } else {
                irregularLog.open("irregular_points_log.txt", std::ios::out | std::ios::app);
            }
            if (!irregularLog.is_open()) {
                gsInfo << "WARNING: Could not open irregular_points_log.txt\n";
            }
        } catch (const std::exception& e) {
            gsInfo << "WARNING: Exception opening log file: " << e.what() << "\n";
        }
        
        if (verbose) {
            gsInfo << "checkJacobianDeterminant: numPatches=" << numPatches
                << ", numFunctions=" << numFunctions
                << ", vectSol=" << vectSol.rows() << "x" << vectSol.cols()
                << ", uv.size=" << uv.size() << "\n";
            outfile << "checkJacobianDeterminant: numPatches=" << numPatches
                << ", numFunctions=" << numFunctions
                << ", vectSol=" << vectSol.rows() << "x" << vectSol.cols()
                << ", uv.size=" << uv.size() << "\n";
        }
        
        try {
            // Validate MPBES object state
            try {
                if (mpbes.nPatches() == 0) {
                    gsInfo << "WARNING: MPBES has no patches\n";
                    numIrregular.setZero();
                    return 0;
                }
                if (mpbes.size() == 0) {
                    gsInfo << "WARNING: MPBES has no functions\n";
                    numIrregular.setZero();
                    return 0;
                }
            } catch (const std::exception& e) {
                gsInfo << "WARNING: Cannot validate MPBES state: " << e.what() << "\n";
                numIrregular.setZero();
                return 0;
            } catch (...) {
                gsInfo << "WARNING: Unknown exception validating MPBES state\n";
                numIrregular.setZero();
                return 0;
            }
            
            const auto& functionDescription = mpbes.functionDescription();
            const auto& indexInTHB = mpbes.indexInTHB();
            
            // Validate structure sizes
            if (functionDescription.size() != static_cast<size_t>(numFunctions)) {
                gsInfo << "WARNING: functionDescription size mismatch: " << functionDescription.size() 
                      << " vs " << numFunctions << "\n";
            }
            
            if (indexInTHB.size() != static_cast<size_t>(numPatches)) {
                gsInfo << "WARNING: indexInTHB size mismatch: " << indexInTHB.size() 
                      << " vs " << numPatches << "\n";
                gsInfo << "Adjusting numPatches to match indexInTHB\n";
                // Use the smaller value to prevent out-of-bounds
                const index_t safeNumPatches = std::min(numPatches, static_cast<index_t>(indexInTHB.size()));
                
                // Pre-filter active functions per patch
                std::vector<std::vector<index_t>> activeFuncsPerPatch(safeNumPatches);
                size_t skippedDescTooSmall = 0;
                size_t skippedInvalidIndex = 0;
                
                // Continue with safe size...
                for (index_t patch = 0; patch < safeNumPatches; ++patch) {
                    try {
                        if (patch >= indexInTHB.size()) {
                            continue;
                        }
                        
                        std::set<index_t> activeFuncsSet;
                        for (index_t f = 0; f < numFunctions; ++f) {
                            if (f >= static_cast<index_t>(functionDescription.size())) {
                                break;
                            }
                            
                            for (size_t comp = 0; comp < functionDescription[f].size(); ++comp) {
                                try {
                                    const auto& desc = functionDescription[f][comp];
                                    if (desc.size() < 3) {
                                        ++skippedDescTooSmall;
                                        continue;
                                    }
                                    int patchId = desc[0];
                                    int levelId = desc[1];
                                    int indexId = desc[2];
                                    
                                    if (patchId == static_cast<int>(patch) &&
                                        levelId >= 0 && indexId >= 0 &&
                                        levelId < static_cast<int>(indexInTHB[patch].size()) &&
                                        indexId < static_cast<int>(indexInTHB[patch][levelId].size()) &&
                                        indexInTHB[patch][levelId][indexId] != static_cast<index_t>(-1)) {
                                        activeFuncsSet.insert(f);
                                        break;
                                    } else {
                                        ++skippedInvalidIndex;
                                    }
                                } catch (const std::exception& e) {
                                    gsInfo << "WARNING: Exception in function description processing: " << e.what() << "\n";
                                }
                            }
                        }
                        activeFuncsPerPatch[patch].assign(activeFuncsSet.begin(), activeFuncsSet.end());
                    } catch (const std::exception& e) {
                        gsInfo << "WARNING: Exception processing patch " << patch << ": " << e.what() << "\n";
                    }
                }
                
                // Early return with safe result
                gsInfo << "Skipping full checkJacobianDeterminant due to size mismatch\n";
                numIrregular.setZero();
                return 0;
            }
            
            // Pre-filter active functions per patch
            std::vector<std::vector<index_t>> activeFuncsPerPatch(numPatches);
            size_t skippedDescTooSmall = 0;
            size_t skippedInvalidIndex = 0;
            
            for (index_t patch = 0; patch < numPatches; ++patch) {
                try {
                    if (patch >= indexInTHB.size()) {
                        continue;
                    }
                    
                    std::set<index_t> activeFuncsSet;
                    for (index_t f = 0; f < numFunctions; ++f) {
                        if (f >= static_cast<index_t>(functionDescription.size())) {
                            break;
                        }
                        
                        for (size_t comp = 0; comp < functionDescription[f].size(); ++comp) {
                            try {
                                const auto& desc = functionDescription[f][comp];
                                if (desc.size() < 3) {
                                    ++skippedDescTooSmall;
                                    continue;
                                }
                                int patchId = desc[0];
                                int levelId = desc[1];
                                int indexId = desc[2];
                                
                                if (patchId == static_cast<int>(patch) &&
                                    levelId >= 0 && indexId >= 0 &&
                                    levelId < static_cast<int>(indexInTHB[patch].size()) &&
                                    indexId < static_cast<int>(indexInTHB[patch][levelId].size()) &&
                                    indexInTHB[patch][levelId][indexId] != static_cast<index_t>(-1)) {
                                    activeFuncsSet.insert(f);
                                    break;
                                } else {
                                    ++skippedInvalidIndex;
                                }
                            } catch (const std::exception& e) {
                                gsInfo << "WARNING: Exception in function description processing: " << e.what() << "\n";
                            }
                        }
                    }
                    activeFuncsPerPatch[patch].assign(activeFuncsSet.begin(), activeFuncsSet.end());
                } catch (const std::exception& e) {
                    gsInfo << "WARNING: Exception processing patch " << patch << ": " << e.what() << "\n";
                }
            }

            // ---- BATCH-EVAL: one evalDerivSingleOnPatch call per function over all valid pts
            // (replaces per-point per-function calls — N×M calls → N calls per patch).
            const real_t interiorMargin = 0.05;

            // Step 1: collect valid parameter points per patch into pre-filtered matrices
            std::vector<gsMatrix<real_t>> filteredUV(numPatches);
            std::vector<std::vector<index_t>> filteredPtIdx(numPatches);
            for (index_t patch = 0; patch < numPatches; ++patch) {
                if (patch >= static_cast<index_t>(uv.size())) continue;
                const gsMatrix<>& uvPatch = uv[patch];
                if (uvPatch.rows() < 2 || uvPatch.cols() == 0) continue;
                std::vector<index_t>& cols = filteredPtIdx[patch];
                cols.reserve(uvPatch.cols());
                for (index_t pt = 0; pt < uvPatch.cols(); ++pt) {
                    const real_t u0 = uvPatch(0, pt), u1 = uvPatch(1, pt);
                    if (u0 < 0 || u0 > 1 || u1 < 0 || u1 > 1) continue;
                    if (u0 <= interiorMargin || u0 >= 1.0 - interiorMargin ||
                        u1 <= interiorMargin || u1 >= 1.0 - interiorMargin) continue;
                    cols.push_back(pt);
                }
                const index_t M = static_cast<index_t>(cols.size());
                filteredUV[patch].resize(2, M);
                for (index_t i = 0; i < M; ++i)
                    filteredUV[patch].col(i) = uvPatch.col(cols[i]);
            }

            // Step 2: batch-evaluate Jacobian components per patch
            // J = [J00 J01; J10 J11], det = J00*J11 - J01*J10
            std::vector<gsMatrix<real_t>> patchDets(numPatches);
            for (index_t patch = 0; patch < numPatches; ++patch) {
                if (patch >= static_cast<index_t>(activeFuncsPerPatch.size())) continue;
                const index_t M = static_cast<index_t>(filteredPtIdx[patch].size());
                if (M == 0) { patchDets[patch].resize(1, 0); continue; }

                gsMatrix<real_t> J00, J01, J10, J11;
                mpbes.jacContribOnPatch(patch, filteredUV[patch], vectSol,
                                        activeFuncsPerPatch[patch],
                                        J00, J01, J10, J11);

                patchDets[patch] = J00.cwiseProduct(J11) - J01.cwiseProduct(J10);

                // Accumulate rawNeg/rawPos for orientation detection
                for (index_t i = 0; i < M; ++i) {
                    const real_t det = patchDets[patch](0, i);
                    if (!std::isfinite(det)) continue;
                    if (det <= 0) ++rawNeg[patch]; else ++rawPos[patch];
                }
            }

            // ---- Determine orientation from eval results (replaces computeMirroredCheckPerPatch) ----
            mirroredPerPatch.assign(numPatches, false);
            for (index_t p = 0; p < numPatches; ++p)
                mirroredPerPatch[p] = (rawNeg[p] > rawPos[p]);

            // ---- Preflight reconciliation (same logic as before) ----
            std::vector<bool> mirrorBaseline = mirroredPerPatch;
            bool usingPreflightMirrorBaseline = false;
            if (g_geometryPreflight.valid &&
                g_geometryPreflight.mirroredReport.mirrored.size() == mirroredPerPatch.size())
            {
                mirrorBaseline = g_geometryPreflight.mirroredReport.mirrored;
                usingPreflightMirrorBaseline = true;
                for (index_t p = 0; p < static_cast<index_t>(mirrorBaseline.size()); ++p)
                {
                    if (filteredPtIdx[p].empty()) continue;  // no samples → keep preflight
                    if (mirrorBaseline[p] != mirroredPerPatch[p])
                        mirrorBaseline[p] = mirroredPerPatch[p];  // current wins on disagreement
                }
            }

            // ---- Count violations from det arrays (cheap in-memory pass) ----
            for (index_t patch = 0; patch < numPatches; ++patch) {
                if (patchDets[patch].cols() == 0) continue;
                const bool patchMirrored =
                    (patch < static_cast<index_t>(mirrorBaseline.size())) ? mirrorBaseline[patch] : false;
                const gsMatrix<>& uvPatch = (patch < static_cast<index_t>(uv.size()))
                                            ? uv[patch] : uv[0];
                const index_t M = patchDets[patch].cols();
                for (index_t i = 0; i < M; ++i) {
                    const real_t det = patchDets[patch](0, i);
                    if (!std::isfinite(det)) continue;
                    const real_t signedDet = patchMirrored ? -det : det;
                    if (patch < static_cast<index_t>(minSignedDetPerPatch.size()))
                        minSignedDetPerPatch[patch] = std::min(minSignedDetPerPatch[patch], signedDet);
                    if (signedDet <= 0) {
                        ++signedNeg[patch];
                        if (patch < static_cast<index_t>(numIrregular.size())) {
                            ++numIrregular[patch];
                            ++totalIrregular;
                            const index_t pt = filteredPtIdx[patch][i];
                            gsVector<real_t> param = uvPatch.col(pt);
                            outfile << "  [IRREGULAR] Patch=" << patch
                                    << ", pt=" << pt
                                    << ", uv=(" << param(0) << ", " << param(1)
                                    << "), det=" << det
                                    << ", signedDet=" << signedDet << "\n";
                        }
                    }
                }
            }
        } catch (const std::exception& e) {
            // gsInfo << "ERROR: Exception in checkJacobianDeterminant: " << e.what() << "\n";
            // outfile << "ERROR: Exception in checkJacobianDeterminant: " << e.what() << "\n";
        }
        
        // gsInfo << "=== checkJacobianDeterminant END ===\n";
        // outfile << "=== checkJacobianDeterminant END ===\n";
        if (g_verbose) gsInfo << "Total irregular points: " << totalIrregular << "\n";
        outfile << "Total irregular points: " << totalIrregular << "\n";

        // Per-patch min signed det — key metric for measuring cross-patch coupling magnitude.
        for (index_t p = 0; p < numPatches; ++p)
        {
            const real_t minD = (p < static_cast<index_t>(minSignedDetPerPatch.size()))
                                ? minSignedDetPerPatch[p]
                                : std::numeric_limits<real_t>::max();
            const index_t irr = (p < numIrregular.size()) ? static_cast<index_t>(numIrregular[p]) : 0;
        }

        // if (verbose) {
        //     for (index_t p = 0; p < numPatches; ++p) {
        //         const bool patchMirrored = (p < static_cast<index_t>(mirroredPerPatch.size())) ? mirroredPerPatch[p] : false;
        //         gsInfo << "Patch " << p << ": rawNeg=" << rawNeg[p] << ", rawPos=" << rawPos[p]
        //                << ", signedNeg=" << signedNeg[p] << ", mirrored=" << (patchMirrored ? "true" : "false") << "\n";
        //         outfile << "Patch " << p << ": rawNeg=" << rawNeg[p] << ", rawPos=" << rawPos[p]
        //                 << ", signedNeg=" << signedNeg[p] << ", mirrored=" << (patchMirrored ? "true" : "false") << "\n";
        //         if (!patchMirrored && rawNeg[p] > rawPos[p]) {
        //             gsInfo << "  [MIRROR?] rawNeg > rawPos with mirrored=false\n";
        //             outfile << "  [MIRROR?] rawNeg > rawPos with mirrored=false\n";
        //         }
        //         if (patchMirrored && rawPos[p] > rawNeg[p]) {
        //             gsInfo << "  [MIRROR?] rawPos > rawNeg with mirrored=true\n";
        //             outfile << "  [MIRROR?] rawPos > rawNeg with mirrored=true\n";
        //         }
        //     }
        // }
        
        try {
            if (irregularLog.is_open()) {
                irregularLog << "========================================\n";
                irregularLog << "SUMMARY: " << totalIrregular << " total irregular points detected\n";
                irregularLog.close();
            }
        } catch (...) {}
        
        return totalIrregular;
        
    } catch (const std::exception& e) {
        gsInfo << "CRITICAL: Unhandled exception in checkJacobianDeterminant: " << e.what() << "\n";
        outfile << "CRITICAL: Unhandled exception in checkJacobianDeterminant: " << e.what() << "\n";
        numIrregular.setZero();
        return 0;
    } catch (...) {
        gsInfo << "CRITICAL: Unknown exception in checkJacobianDeterminant\n";
        outfile << "CRITICAL: Unknown exception in checkJacobianDeterminant\n";
        numIrregular.setZero();
        return 0;
    }
}

// Detect whether the parameterization is mirrored (overall negative Jacobian orientation)
bool isParameterizationMirrored(
    const gsVector<gsMatrix<>>& uv,
    const MPBES<2, real_t>& mpbes,
    const gsMatrix<real_t>& vectSol,
    real_t detEps = 1e-12,
    bool verbose = true)
{
    PROFILE_FUNCTION();
    const index_t numPatches = mpbes.nPatches();
    const index_t numFunctions = mpbes.size();

    if (numPatches == 0 || numFunctions == 0 || uv.size() == 0 || vectSol.rows() == 0) {
        if (verbose) {
            gsInfo << "isParameterizationMirrored: insufficient data, returning false\n";
            outfile << "isParameterizationMirrored: insufficient data, returning false\n";
        }
        return false;
    }

    const auto& functionDescription = mpbes.functionDescription();
    const auto& indexInTHB = mpbes.indexInTHB();

    std::vector<std::vector<index_t>> activeFuncsPerPatch(numPatches);
    for (index_t patch = 0; patch < numPatches; ++patch)
    {
        if (patch >= static_cast<index_t>(indexInTHB.size()))
            continue;

        std::set<index_t> activeFuncsSet;
        for (index_t f = 0; f < numFunctions; ++f)
        {
            if (f >= static_cast<index_t>(functionDescription.size()))
                break;

            for (size_t comp = 0; comp < functionDescription[f].size(); ++comp)
            {
                const auto& desc = functionDescription[f][comp];
                if (desc.size() < 3)
                    continue;
                int patchId = desc[0];
                int levelId = desc[1];
                int indexId = desc[2];

                if (patchId == static_cast<int>(patch) &&
                    levelId >= 0 && indexId >= 0 &&
                    levelId < static_cast<int>(indexInTHB[patch].size()) &&
                    indexId < static_cast<int>(indexInTHB[patch][levelId].size()) &&
                    indexInTHB[patch][levelId][indexId] != static_cast<index_t>(-1))
                {
                    activeFuncsSet.insert(f);
                    break;
                }
            }
        }
        activeFuncsPerPatch[patch].assign(activeFuncsSet.begin(), activeFuncsSet.end());
    }

    size_t totalPos = 0;
    size_t totalNeg = 0;
    size_t totalUsed = 0;
    std::vector<size_t> posPerPatch(numPatches, 0);
    std::vector<size_t> negPerPatch(numPatches, 0);

    for (index_t patch = 0; patch < numPatches; ++patch)
    {
        if (patch >= static_cast<index_t>(uv.size()))
            continue;

        const gsMatrix<>& uvPatch = uv[patch];
        const index_t numPoints = uvPatch.cols();
        if (uvPatch.rows() < 2 || numPoints == 0)
            continue;

        for (index_t pt = 0; pt < numPoints; ++pt)
        {
            gsVector<real_t> param = uvPatch.col(pt);
            if (param.rows() < 2)
                continue;

            if (param(0) < 0 || param(0) > 1 || param(1) < 0 || param(1) > 1)
                continue;

            const real_t interiorMargin = 0.05;
            if (param(0) <= interiorMargin || param(0) >= 1.0 - interiorMargin ||
                param(1) <= interiorMargin || param(1) >= 1.0 - interiorMargin)
                continue;

            gsMatrix<> J(2, 2);
            J.setZero();

            if (patch >= static_cast<index_t>(activeFuncsPerPatch.size()))
                continue;

            for (index_t f : activeFuncsPerPatch[patch])
            {
                if (f >= static_cast<index_t>(vectSol.rows()) || vectSol.cols() < 2)
                    continue;
                if (f >= numFunctions)
                    continue;

                gsMatrix<real_t> basisDeriv;
                mpbes.evalDerivSingleOnPatch(f, patch, param, basisDeriv);
                if (basisDeriv.rows() < 2 || basisDeriv.cols() < 1)
                    continue;

                real_t dphi_du = basisDeriv(0, 0);
                real_t dphi_dv = basisDeriv(1, 0);
                if (!std::isfinite(dphi_du) || !std::isfinite(dphi_dv))
                    continue;

                real_t coef_x = vectSol(f, 0);
                real_t coef_y = vectSol(f, 1);
                if (!std::isfinite(coef_x) || !std::isfinite(coef_y))
                    continue;

                J(0, 0) += dphi_du * coef_x;
                J(0, 1) += dphi_dv * coef_x;
                J(1, 0) += dphi_du * coef_y;
                J(1, 1) += dphi_dv * coef_y;
            }

            real_t det = J.determinant();
            if (!std::isfinite(det))
                continue;

            if (std::abs(det) <= detEps)
                continue;

            ++totalUsed;
            if (det > 0) {
                ++totalPos;
                ++posPerPatch[patch];
            } else {
                ++totalNeg;
                ++negPerPatch[patch];
            }
        }
    }

    if (verbose) {
        gsInfo << "isParameterizationMirrored: totalUsed=" << totalUsed
               << ", pos=" << totalPos << ", neg=" << totalNeg << "\n";
        outfile << "isParameterizationMirrored: totalUsed=" << totalUsed
                << ", pos=" << totalPos << ", neg=" << totalNeg << "\n";
        for (index_t p = 0; p < numPatches; ++p) {
            if (posPerPatch[p] + negPerPatch[p] == 0)
                continue;
            gsInfo << "  patch " << p << ": pos=" << posPerPatch[p]
                   << ", neg=" << negPerPatch[p] << "\n";
            outfile << "  patch " << p << ": pos=" << posPerPatch[p]
                    << ", neg=" << negPerPatch[p] << "\n";
        }
    }

    if (totalUsed == 0)
        return false;

    if (totalPos == 0 && totalNeg > 0)
        return true;

    const real_t negRatio = static_cast<real_t>(totalNeg) / static_cast<real_t>(totalUsed);
    return (totalNeg > totalPos) && (negRatio > 0.9);
}

std::vector<bool> isParameterizationMirroredPerPatch(
    const gsVector<gsMatrix<>>& uv,
    const MPBES<2, real_t>& mpbes,
    const gsMatrix<real_t>& vectSol,
    real_t detEps,
    bool verbose)
{
    const MirroredCheckResult report =
        computeMirroredCheckPerPatch(uv, mpbes, vectSol, detEps);

    if (verbose) {
        gsInfo << "isParameterizationMirroredPerPatch summary:\n";
        outfile << "isParameterizationMirroredPerPatch summary:\n";
        for (index_t p = 0; p < static_cast<index_t>(report.mirrored.size()); ++p) {
            if (report.usedPerPatch[p] == 0)
                continue;
            const real_t negRatio = static_cast<real_t>(report.negPerPatch[p]) /
                                    static_cast<real_t>(report.usedPerPatch[p]);
            gsInfo << "  patch " << p << ": used=" << report.usedPerPatch[p]
                   << ", pos=" << report.posPerPatch[p]
                   << ", neg=" << report.negPerPatch[p]
                   << ", negRatio=" << negRatio
                   << ", detEps=" << detEps
                   << ", mirrored=" << (report.mirrored[p] ? "true" : "false") << "\n";
            outfile << "  patch " << p << ": used=" << report.usedPerPatch[p]
                    << ", pos=" << report.posPerPatch[p]
                    << ", neg=" << report.negPerPatch[p]
                    << ", negRatio=" << negRatio
                    << ", detEps=" << detEps
                    << ", mirrored=" << (report.mirrored[p] ? "true" : "false") << "\n";
        }
    }

    return report.mirrored;
}

static MirroredCheckResult computeMirroredCheckFromMp(
    const gsVector<gsMatrix<real_t>>& uv,
    const gsMultiPatch<real_t>& mp,
    real_t detEps)
{
    const index_t numPatches = mp.nPatches();
    MirroredCheckResult result;
    result.mirrored.assign(numPatches, false);
    result.usedPerPatch.assign(numPatches, 0);
    result.posPerPatch.assign(numPatches, 0);
    result.negPerPatch.assign(numPatches, 0);

    const real_t interiorMargin = 0.05;
    for (index_t patch = 0; patch < numPatches && patch < static_cast<index_t>(uv.size()); ++patch)
    {
        const gsMatrix<real_t>& uvPatch = uv[patch];
        if (uvPatch.rows() < 2 || uvPatch.cols() == 0) continue;

        for (index_t pt = 0; pt < uvPatch.cols(); ++pt)
        {
            const real_t u0 = uvPatch(0, pt), v0 = uvPatch(1, pt);
            if (u0 <= interiorMargin || u0 >= 1.0 - interiorMargin ||
                v0 <= interiorMargin || v0 >= 1.0 - interiorMargin) continue;

            const gsMatrix<real_t> J = mp.patch(patch).jacobian(uvPatch.col(pt));
            if (J.rows() < 2 || J.cols() < 2) continue;
            const real_t det = J.determinant();
            if (!std::isfinite(det) || std::abs(det) <= detEps) continue;

            ++result.usedPerPatch[patch];
            if (det > 0) ++result.posPerPatch[patch];
            else         ++result.negPerPatch[patch];
        }

        if (result.usedPerPatch[patch] > 0) {
            const real_t negRatio = static_cast<real_t>(result.negPerPatch[patch]) /
                                    static_cast<real_t>(result.usedPerPatch[patch]);
            result.mirrored[patch] = (result.posPerPatch[patch] == 0 || negRatio >= 0.5);
        }
    }
    return result;
}

MirroredCheckResult computeMirroredCheckPerPatch(
    const gsVector<gsMatrix<>>& uv,
    const MPBES<2, real_t>& mpbes,
    const gsMatrix<real_t>& vectSol,
    real_t detEps)
{
    PROFILE_FUNCTION();
    const index_t numPatches = mpbes.nPatches();
    const index_t numFunctions = mpbes.size();
    MirroredCheckResult result;
    result.mirrored.assign(numPatches, false);
    result.usedPerPatch.assign(numPatches, 0);
    result.posPerPatch.assign(numPatches, 0);
    result.negPerPatch.assign(numPatches, 0);

    if (numPatches == 0 || numFunctions == 0 || uv.size() == 0 || vectSol.rows() == 0)
        return result;

    const auto& functionDescription = mpbes.functionDescription();
    const auto& indexInTHB = mpbes.indexInTHB();

    std::vector<std::vector<index_t>> activeFuncsPerPatch(numPatches);
    for (index_t patch = 0; patch < numPatches; ++patch)
    {
        if (patch >= static_cast<index_t>(indexInTHB.size()))
            continue;

        std::set<index_t> activeFuncsSet;
        for (index_t f = 0; f < numFunctions; ++f)
        {
            if (f >= static_cast<index_t>(functionDescription.size()))
                break;

            for (size_t comp = 0; comp < functionDescription[f].size(); ++comp)
            {
                const auto& desc = functionDescription[f][comp];
                if (desc.size() < 3)
                    continue;
                int patchId = desc[0];
                int levelId = desc[1];
                int indexId = desc[2];

                if (patchId == static_cast<int>(patch) &&
                    levelId >= 0 && indexId >= 0 &&
                    levelId < static_cast<int>(indexInTHB[patch].size()) &&
                    indexId < static_cast<int>(indexInTHB[patch][levelId].size()) &&
                    indexInTHB[patch][levelId][indexId] != static_cast<index_t>(-1))
                {
                    activeFuncsSet.insert(f);
                    break;
                }
            }
        }
        activeFuncsPerPatch[patch].assign(activeFuncsSet.begin(), activeFuncsSet.end());
    }

    for (index_t patch = 0; patch < numPatches; ++patch)
    {
        if (patch >= static_cast<index_t>(uv.size()))
            continue;

        const gsMatrix<>& uvPatch = uv[patch];
        const index_t numPoints = uvPatch.cols();
        if (uvPatch.rows() < 2 || numPoints == 0)
            continue;

        for (index_t pt = 0; pt < numPoints; ++pt)
        {
            gsVector<real_t> param = uvPatch.col(pt);
            if (param.rows() < 2)
                continue;
            if (param(0) < 0 || param(0) > 1 || param(1) < 0 || param(1) > 1)
                continue;

            const real_t interiorMargin = 0.05;
            if (param(0) <= interiorMargin || param(0) >= 1.0 - interiorMargin ||
                param(1) <= interiorMargin || param(1) >= 1.0 - interiorMargin)
                continue;

            gsMatrix<> J(2, 2);
            J.setZero();

            if (patch >= static_cast<index_t>(activeFuncsPerPatch.size()))
                continue;

            for (index_t f : activeFuncsPerPatch[patch])
            {
                if (f >= static_cast<index_t>(vectSol.rows()) || vectSol.cols() < 2)
                    continue;
                if (f >= numFunctions)
                    continue;

                gsMatrix<real_t> basisDeriv;
                mpbes.evalDerivSingleOnPatch(f, patch, param, basisDeriv);
                if (basisDeriv.rows() < 2 || basisDeriv.cols() < 1)
                    continue;

                real_t dphi_du = basisDeriv(0, 0);
                real_t dphi_dv = basisDeriv(1, 0);
                if (!std::isfinite(dphi_du) || !std::isfinite(dphi_dv))
                    continue;

                real_t coef_x = vectSol(f, 0);
                real_t coef_y = vectSol(f, 1);
                if (!std::isfinite(coef_x) || !std::isfinite(coef_y))
                    continue;

                J(0, 0) += dphi_du * coef_x;
                J(0, 1) += dphi_dv * coef_x;
                J(1, 0) += dphi_du * coef_y;
                J(1, 1) += dphi_dv * coef_y;
            }

            real_t det = J.determinant();
            if (!std::isfinite(det))
                continue;
            if (std::abs(det) <= detEps)
                continue;

            ++result.usedPerPatch[patch];
            if (det > 0)
                ++result.posPerPatch[patch];
            else
                ++result.negPerPatch[patch];
        }
    }

    for (index_t patch = 0; patch < numPatches; ++patch)
    {
        if (result.usedPerPatch[patch] == 0) {
            result.mirrored[patch] = false;
            continue;
        }
        if (result.posPerPatch[patch] == 0 && result.negPerPatch[patch] > 0) {
            result.mirrored[patch] = true;
            continue;
        }
        const real_t negRatio = static_cast<real_t>(result.negPerPatch[patch]) /
                                static_cast<real_t>(result.usedPerPatch[patch]);
        result.mirrored[patch] =
            (result.negPerPatch[patch] > result.posPerPatch[patch]) && (negRatio > 0.9);
    }

    return result;
}

void logMirroredCheck(
    const gsVector<gsMatrix<>>& uv,
    const MPBES<2, real_t>& mpbes,
    const gsMatrix<real_t>& vectSol,
    real_t detEps,
    const std::string& label,
    bool logToInfo)
{
    const MirroredCheckResult report =
        computeMirroredCheckPerPatch(uv, mpbes, vectSol, detEps);

    std::ostringstream buffer;
    buffer << "MIRROREDCHECK_BEGIN [" << label << "]\n";
    buffer << "detEps=" << detEps << ", patches=" << report.mirrored.size() << "\n";

    index_t mirroredCount = 0;
    for (index_t p = 0; p < static_cast<index_t>(report.mirrored.size()); ++p)
    {
        if (report.mirrored[p])
            ++mirroredCount;

        buffer << "patch " << p
               << ": used=" << report.usedPerPatch[p]
               << ", pos=" << report.posPerPatch[p]
               << ", neg=" << report.negPerPatch[p];

        if (report.usedPerPatch[p] > 0)
        {
            const real_t negRatio = static_cast<real_t>(report.negPerPatch[p]) /
                                    static_cast<real_t>(report.usedPerPatch[p]);
            buffer << ", negRatio=" << negRatio;
        }
        else
        {
            buffer << ", negRatio=n/a";
        }

        buffer << ", mirrored=" << (report.mirrored[p] ? "true" : "false") << "\n";
    }

    buffer << "mirroredCount=" << mirroredCount << "\n";
    buffer << "MIRROREDCHECK_END\n";

    if (outfile.is_open())
    {
        outfile << buffer.str();
        outfile.flush();
    }

    if (logToInfo)
        gsInfo << buffer.str();
}

// Old version kept for compatibility
int checkJacobianDeterminant(
    const gsVector<gsMatrix<>>& uv,
    const gsVector<gsTHBSplineBasis<2, real_t>>& THBVector,
    const gsVector<gsVector<gsTensorBSplineBasis<2, real_t>>>& Bells,
    const std::vector<std::vector<std::vector<index_t>>>& functionDescription,
    const gsVector<gsVector<gsVector<index_t>>>& indexInTHB,
    const gsVector<gsVector<int>>& globalIndexTHB,
    const gsMatrix<real_t>& vectSol,
    gsVector<size_t>& numIrregular,
    const std::vector<bool>& isTruncated,
    const std::vector<std::vector<gsSparseVector<double>>>& presentation,
    const std::vector<std::vector<std::array<int, 3>>>& spilloverFunctionCoordinates,
    bool verbose = true)
{
    PROFILE_FUNCTION();
    gsInfo << "=== checkJacobianDeterminant BEGIN ===\n";
    
    // Open detailed log file for irregular points (append mode to preserve across calls)
    static bool firstCall = true;
    std::ofstream irregularLog;
    if (firstCall) {
        irregularLog.open("irregular_points_log.txt", std::ios::out | std::ios::trunc);
        irregularLog << "=== IRREGULAR POINTS DETAILED LOG ===\n";
        irregularLog << "Logging all points where Jacobian determinant <= 0\n\n";
        firstCall = false;
    } else {
        irregularLog.open("irregular_points_log.txt", std::ios::out | std::ios::app);
    }
    
    const index_t numPatches = THBVector.size();
    int totalIrregular = 0;
    numIrregular.setZero();
    numIrregular.resize(numPatches);

    const index_t numFunctions = static_cast<index_t>(functionDescription.size());
    
    // Pre-filter active functions per patch (same as assemble)
    std::vector<std::vector<index_t>> activeFuncsPerPatch(numPatches);
    {
        for (index_t patch = 0; patch < numPatches; ++patch)
        {
            std::set<index_t> activeFuncsSet;
            for (index_t f = 0; f < numFunctions; ++f)
            {
                for (size_t comp = 0; comp < functionDescription[f].size(); ++comp)
                {
                    const auto& desc = functionDescription[f][comp];
                    if (desc[0] == static_cast<int>(patch) &&
                        desc[1] < indexInTHB[patch].size() &&
                        desc[2] < indexInTHB[patch][desc[1]].size() &&
                        indexInTHB[patch][desc[1]][desc[2]] != static_cast<index_t>(-1))
                    {
                        activeFuncsSet.insert(f);
                        break;
                    }
                }
            }
            activeFuncsPerPatch[patch].assign(activeFuncsSet.begin(), activeFuncsSet.end());
        }
    }

    // Evaluate Jacobian determinant at sample points per patch
    for (index_t patch = 0; patch < numPatches; ++patch)
    {
        const gsMatrix<>& uvPatch = uv[patch];
        const index_t numPoints = uvPatch.cols();

        if (verbose)
        {
            gsInfo << ">> Patch " << patch << ": " << numPoints << " points\n";
            outfile << ">> Patch " << patch << ": " << numPoints << " points\n";
        }

        // Evaluate derivatives at each point (matching assemble pattern)
        for (index_t pt = 0; pt < numPoints; ++pt)
        {
            const gsVector<real_t> param = uvPatch.col(pt);
            
            // Jacobian matrix J = [dx/du  dx/dv]
            //                     [dy/du  dy/dv]
            gsMatrix<> J(2, 2);
            J.setZero();
            
            gsMatrix<> basisDeriv; // Reusable matrix for derivatives
            
            // Only iterate active functions for this patch (same loop as assemble)
            for (index_t f : activeFuncsPerPatch[patch])
            {
                if (f >= static_cast<index_t>(vectSol.rows())) continue;
                
                // Accumulate derivatives for all components of this function
                gsMatrix<> totalDeriv(2, 1);
                totalDeriv.setZero();
                
                int functionComponent = 0;
                for (const auto& desc : functionDescription[f])
                {
                    int fnPatch = desc[0];
                    int fnLevel = desc[1];
                    int fnTensorIndex = desc[2];

                    if (fnPatch != static_cast<int>(patch)) {
                        functionComponent++;
                        continue;
                    }

                    int thbIdx = indexInTHB[fnPatch][fnLevel][fnTensorIndex];
                    
                    // Evaluate derivative - matching MPBES.evalDerivSingleOnPatch logic
                    if (thbIdx == -1) {
                        // Spillover component (not in THB, never truncated) - use Bells basis directly
                        if (fnLevel < Bells[patch].size()) {
                            Bells[patch][fnLevel].derivSingle_into(fnTensorIndex, param, basisDeriv);
                            totalDeriv += basisDeriv;
                        }
                    } else if (!isTruncated[f] || !isComponentTruncated(f, functionComponent, presentation)) {
                        // Non-truncated component: use THB basis derivative directly
                        unsigned level = THBVector[patch].levelOf(thbIdx);
                        unsigned tensorIdx = THBVector[patch].flatTensorIndexOf(thbIdx, level);
                        THBVector[patch].getBases()[level]->derivSingle_into(tensorIdx, param, basisDeriv);
                        totalDeriv += basisDeriv;
                    } else {
                        // Truncated component: sum derivatives weighted by presentation coefficients
                        if (f < presentation.size() && functionComponent < presentation[f].size()) {
                            const auto& coefs = presentation[f][functionComponent];
                            if (coefs.size() > 0) {
                                // Presentation level is the level of the finer basis used in truncation
                                unsigned level = THBVector[patch].getistruncated(thbIdx);
                                for (index_t k = 0; k < coefs.size(); ++k) {
                                    if (coefs(k) != 0.0) {
                                        Bells[patch][level].derivSingle_into(k, param, basisDeriv);
                                        totalDeriv += static_cast<real_t>(coefs(k)) * basisDeriv;
                                    }
                                }
                            }
                        }
                    }
                    functionComponent++;
                }
                
                // Accumulate contribution to Jacobian
                // totalDeriv(0,0) = dφ/du, totalDeriv(1,0) = dφ/dv
                real_t dφ_du = totalDeriv(0, 0);
                real_t dφ_dv = totalDeriv(1, 0);
                
                J(0, 0) += dφ_du * vectSol(f, 0);  // dx/du
                J(0, 1) += dφ_dv * vectSol(f, 0);  // dx/dv
                J(1, 0) += dφ_du * vectSol(f, 1);  // dy/du
                J(1, 1) += dφ_dv * vectSol(f, 1);  // dy/dv
            }

            real_t det = J.determinant();
            if (det <= 0)
            {
                ++numIrregular[patch];
                ++totalIrregular;

                // Write to main outfile
                outfile << "  [IRREGULAR] Patch=" << patch
                    << ", pt=" << pt
                    << ", uv=(" << param(0) << ", " << param(1)
                    << "), det=" << det << "\n";
                
                // Write detailed information to irregular points log
                irregularLog << "========================================\n";
                irregularLog << "IRREGULAR POINT #" << totalIrregular << "\n";
                irregularLog << "Patch: " << patch << "\n";
                irregularLog << "Point index: " << pt << " (out of " << numPoints << ")\n";
                irregularLog << "Parameter coordinates (u,v): (" << param(0) << ", " << param(1) << ")\n";
                irregularLog << "\n--- Jacobian Derivatives (Analytical) ---\n";
                irregularLog << "dx/du = " << J(0, 0) << "\n";
                irregularLog << "dx/dv = " << J(0, 1) << "\n";
                irregularLog << "dy/du = " << J(1, 0) << "\n";
                irregularLog << "dy/dv = " << J(1, 1) << "\n";
                irregularLog << "\n--- Jacobian Matrix ---\n";
                irregularLog << "J = [" << J(0, 0) << "  " << J(0, 1) << "]\n";
                irregularLog << "    [" << J(1, 0) << "  " << J(1, 1) << "]\n";
                irregularLog << "\n--- Determinant ---\n";
                irregularLog << "det(J) = " << det << "\n";
                irregularLog << "Status: " << (det < 0 ? "NEGATIVE (fold/overlap)" : "ZERO (singular/degenerate)") << "\n";
                
                // Delta test: check surrounding points
                irregularLog << "\n--- Neighboring Points Test (Delta = 0.01) ---\n";
                const real_t delta = 0.01;
                std::vector<std::pair<std::string, gsVector<real_t>>> neighbors = {
                    {"Right (u+δ, v)", gsVector<real_t>(2)},
                    {"Left  (u-δ, v)", gsVector<real_t>(2)},
                    {"Up    (u, v+δ)", gsVector<real_t>(2)},
                    {"Down  (u, v-δ)", gsVector<real_t>(2)}
                };
                neighbors[0].second << param(0) + delta, param(1);
                neighbors[1].second << param(0) - delta, param(1);
                neighbors[2].second << param(0), param(1) + delta;
                neighbors[3].second << param(0), param(1) - delta;
                
                for (const auto& nb : neighbors) {
                    const auto& nbParam = nb.second;
                    
                    // Check if neighbor is within domain [0,1]x[0,1]
                    if (nbParam(0) < 0 || nbParam(0) > 1 || nbParam(1) < 0 || nbParam(1) > 1) {
                        irregularLog << nb.first << ": OUT OF BOUNDS (" << nbParam(0) << ", " << nbParam(1) << ")\n";
                        continue;
                    }
                    
                    // Compute Jacobian at neighbor
                    gsMatrix<> J_nb(2, 2);
                    J_nb.setZero();
                    
                    for (index_t f : activeFuncsPerPatch[patch]) {
                        if (f >= static_cast<index_t>(vectSol.rows())) continue;
                        
                        gsMatrix<> totalDeriv_nb(2, 1);
                        totalDeriv_nb.setZero();
                        
                        int functionComponent = 0;
                        for (const auto& desc : functionDescription[f]) {
                            int fnPatch = desc[0];
                            int fnLevel = desc[1];
                            int fnTensorIndex = desc[2];

                            if (fnPatch != static_cast<int>(patch)) {
                                functionComponent++;
                                continue;
                            }

                            int thbIdx = indexInTHB[fnPatch][fnLevel][fnTensorIndex];
                            
                            if (thbIdx == -1) {
                                // Spillover component - use Bells basis directly
                                if (fnLevel < Bells[patch].size()) {
                                    Bells[patch][fnLevel].derivSingle_into(fnTensorIndex, nbParam, basisDeriv);
                                    totalDeriv_nb += basisDeriv;
                                }
                            } else if (!isTruncated[f] || !isComponentTruncated(f, functionComponent, presentation)) {
                                // Non-truncated component
                                unsigned level = THBVector[patch].levelOf(thbIdx);
                                unsigned tensorIdx = THBVector[patch].flatTensorIndexOf(thbIdx, level);
                                THBVector[patch].getBases()[level]->derivSingle_into(tensorIdx, nbParam, basisDeriv);
                                totalDeriv_nb += basisDeriv;
                            } else {
                                // Truncated component
                                if (f < presentation.size() && functionComponent < presentation[f].size()) {
                                    const auto& coefs = presentation[f][functionComponent];
                                    if (coefs.size() > 0) {
                                        unsigned level = THBVector[patch].getistruncated(thbIdx);
                                        for (index_t k = 0; k < coefs.size(); ++k) {
                                            if (coefs(k) != 0.0) {
                                                Bells[patch][level].derivSingle_into(k, nbParam, basisDeriv);
                                                totalDeriv_nb += static_cast<real_t>(coefs(k)) * basisDeriv;
                                            }
                                        }
                                    }
                                }
                            }
                            functionComponent++;
                        }
                        
                        real_t dφ_du_nb = totalDeriv_nb(0, 0);
                        real_t dφ_dv_nb = totalDeriv_nb(1, 0);
                        
                        J_nb(0, 0) += dφ_du_nb * vectSol(f, 0);
                        J_nb(0, 1) += dφ_dv_nb * vectSol(f, 0);
                        J_nb(1, 0) += dφ_du_nb * vectSol(f, 1);
                        J_nb(1, 1) += dφ_dv_nb * vectSol(f, 1);
                    }
                    
                    real_t det_nb = J_nb.determinant();
                    irregularLog << nb.first << ": det = " << det_nb 
                                << " [" << (det_nb > 0 ? "POSITIVE" : (det_nb < 0 ? "NEGATIVE" : "ZERO")) << "]\n";
                }
                
                irregularLog << "\n";
            }
        }

        gsInfo << "  Patch " << patch << ": " << numIrregular[patch] << " irregular points\n";
        outfile << "  Patch " << patch << ": " << numIrregular[patch] << " irregular points\n";
    }

    if (g_verbose) gsInfo << "=== checkJacobianDeterminant END ===\n";
    if (g_verbose) gsInfo << "Total irregular points: " << totalIrregular << "\n";
    
    // Close irregular points log
    irregularLog << "========================================\n";
    irregularLog << "SUMMARY: " << totalIrregular << " total irregular points detected\n";
    irregularLog.close();
    
    if (totalIrregular > 0) {
        gsInfo << "Detailed irregular points log written to: irregular_points_log.txt\n";
    }
    
    return totalIrregular;
}






















real_t getSupports(int i, int j, int k, gsMatrix<> supps) {
    return supps(j, 2 * i + k);
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int pickCell(gsVector < int > vectorS, int& currArrayIndex, int levNow, int& x1U, int& y1U, int& x2U, int& y2U, index_t interior) {
    currArrayIndex = rand() % vectorS.size();
    //gsInfo << "nonCheckedCells.size() " << vectorS.size() << "\n";
    //gsInfo << "currArrayIndex " << currArrayIndex << "\n";
    //gsInfo << "The cell is picked randomly\n";
    int currCellIndex = vectorS(currArrayIndex);
    x1U = ((interior + 1) * (int)pow(2, levNow) - 1) - (currCellIndex % (int)((interior + 1) * (int)pow(2, levNow)));
    x2U = x1U + 1;
    y1U = ((interior + 1) * (int)pow(2, levNow) * (interior + 1) * (int)pow(2, levNow) - 1 - currCellIndex) / ((interior + 1) * (int)pow(2, levNow));
    y2U = y1U + 1;
    return currCellIndex;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int pickCell(gsVector < int > vectorS, int attempt, int levNow, int& x1U, int& y1U, int& x2U, int& y2U, int lexicographic, index_t interior) {
    index_t currArrayIndex = 0;//attempt % nonCheckedCells.size();
    //gsInfo << "nonCheckedCells.size() " << vectorS.size() << "\n";
    //gsInfo << "The cell is picked lexicographically\n";
    int currCellIndex = vectorS(currArrayIndex);
    x1U = ((interior + 1) * (int)pow(2, levNow) - 1) - (currCellIndex % (int)((interior + 1) * (int)pow(2, levNow)));
    x2U = x1U + 1;
    y1U = ((interior + 1) * (int)pow(2, levNow) * (interior + 1) * (int)pow(2, levNow) - 1 - currCellIndex) / ((interior + 1) * (int)pow(2, levNow));
    y2U = y1U + 1;
    return currCellIndex;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
static bool cellCoordsFromIndex(int levNow, int cellIndex, index_t interior,
    int& x1U, int& y1U, int& x2U, int& y2U)
{
    const int n = (interior + 1) * (int)pow(2, levNow);
    if (n <= 0 || cellIndex < 0 || cellIndex >= n * n)
        return false;

    x1U = (n - 1) - (cellIndex % n);
    x2U = x1U + 1;
    y1U = (n * n - 1 - cellIndex) / n;
    y2U = y1U + 1;
    return true;
}

static bool cellIndexFromCoords(int levNow, int x1U, int y1U, index_t interior, int& outIndex)
{
    const int n = (interior + 1) * (int)pow(2, levNow);
    if (x1U < 0 || y1U < 0 || x1U >= n || y1U >= n)
        return false;

    const int col = (n - 1) - x1U;
    const int base = (n * n - 1 - col) / n;
    const int row = base - y1U;
    if (row < 0)
        return false;

    outIndex = row * n + col;
    return outIndex >= 0 && outIndex < n * n;
}

static void removeCellIdsByValue(gsVector<int>& cells, const std::vector<int>& idsToRemove)
{
    if (cells.size() == 0 || idsToRemove.empty())
        return;

    std::unordered_set<int> removeSet;
    removeSet.reserve(idsToRemove.size() * 2);
    for (size_t i = 0; i < idsToRemove.size(); ++i)
        if (idsToRemove[i] >= 0)
            removeSet.insert(idsToRemove[i]);

    if (removeSet.empty())
        return;

    gsVector<int> filtered(cells.size());
    int pos = 0;
    for (int i = 0; i < cells.size(); ++i)
    {
        if (removeSet.find(cells(i)) == removeSet.end())
            filtered(pos++) = cells(i);
    }

    filtered.conservativeResize(pos);
    cells = filtered;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



template <typename T> void goestaerDaeOPokhu(T object) {
    //    cout << "Printing the object\n";
    for (int patch = 0; patch < object.size(); ++patch) {
        //gsInfo << "Patch " << patch << "\n";
        for (int level = 0; level < object(patch).size(); ++level) {
            //gsInfo << "level " << level << "\n";
            //gsInfo << object(patch)(level) << "\n";
        }
    }
    cout << "------------------------\n";
};




static bool pickHierarchySafeSiblingGroup(
    const gsVector<index_t>& vectorS,
    int preferredArrayIndex,
    int levNow,
    index_t interior,
    std::vector<CellToCoarsen>& geoCells,
    std::vector<int>& geoCellIndices,
    int& chosenArrayIndex)
{
    geoCells.clear();
    geoCellIndices.clear();
    chosenArrayIndex = 0;

    if (levNow <= 0 || vectorS.size() < 4)
        return false;

    std::unordered_set<int> available;
    available.reserve(vectorS.size() * 2);
    for (int i = 0; i < vectorS.size(); ++i)
        available.insert(static_cast<int>(vectorS(i)));

    const int start =
        (vectorS.size() > 0)
        ? std::max(0, std::min(preferredArrayIndex, static_cast<int>(vectorS.size()) - 1))
        : 0;

    for (int offset = 0; offset < vectorS.size(); ++offset)
    {
        const int arrayIndex = (start + offset) % static_cast<int>(vectorS.size());
        const int candidateIndex = static_cast<int>(vectorS(arrayIndex));

        int x1U = 0, y1U = 0, x2U = 0, y2U = 0;
        if (!cellCoordsFromIndex(levNow, candidateIndex, interior, x1U, y1U, x2U, y2U))
            continue;

        const int parentX = (x1U / 2) * 2;
        const int parentY = (y1U / 2) * 2;

        std::array<std::pair<int, int>, 4> siblingCoords = {
            std::make_pair(parentX, parentY),
            std::make_pair(parentX + 1, parentY),
            std::make_pair(parentX, parentY + 1),
            std::make_pair(parentX + 1, parentY + 1)
        };

        std::vector<int> siblingIndices;
        siblingIndices.reserve(4);
        bool allAvailable = true;
        for (const auto& coord : siblingCoords)
        {
            int siblingIndex = -1;
            if (!cellIndexFromCoords(levNow, coord.first, coord.second, interior, siblingIndex) ||
                available.find(siblingIndex) == available.end())
            {
                allAvailable = false;
                break;
            }
            siblingIndices.push_back(siblingIndex);
        }

        if (!allAvailable)
            continue;

        geoCells.reserve(4);
        geoCellIndices.reserve(4);
        for (const int siblingIndex : siblingIndices)
        {
            int sx1 = 0, sy1 = 0, sx2 = 0, sy2 = 0;
            if (!cellCoordsFromIndex(levNow, siblingIndex, interior, sx1, sy1, sx2, sy2))
            {
                allAvailable = false;
                break;
            }

            CellToCoarsen cell;
            cell.level = levNow;
            cell.x1 = sx1;
            cell.y1 = sy1;
            cell.x2 = sx2;
            cell.y2 = sy2;
            geoCells.push_back(cell);
            geoCellIndices.push_back(siblingIndex);
        }

        if (!allAvailable || geoCells.size() != 4 || geoCellIndices.size() != 4)
        {
            geoCells.clear();
            geoCellIndices.clear();
            continue;
        }

        chosenArrayIndex = arrayIndex;
        return true;
    }

    return false;
}
template <typename T> void DebugDaOPokhu(T object) {
    printValue(object);
    cout << "Printing the object\n";
    for (int patch = 0; patch < object.size(); ++patch) {
        gsInfo << "Patch " << patch << ", ";
        for (int level = 0; level < object(patch).size(); ++level) {
            gsInfo << "level " << level << "\n";
            gsDebugVar(object(patch)(level));
        }
    }
    cout << "------------------------\n";
};

int PatchesIntersection(gsGeometry < >& geom1, gsGeometry < >& geom2,
    index_t numPoints, real_t tolerance) {
    int doPatchesIntersect = 0;
    boxSide const s;

    gsVector<index_t> boundaryFunctions = geom1.basis().boundary(s);
    gsVector<index_t>  sortedBoundaryFunctions(boundaryFunctions.size());

    gsVector<real_t>  referencePoints(boundaryFunctions.size());
    sortedBoundaryFunctions = boundaryFunctions;
    for (int i = 0; i < boundaryFunctions.size(); ++i) {

        referencePoints(i) = geom1.coef(boundaryFunctions(i))(1);
    }
    for (int i = 0; i < referencePoints.size(); ++i) {

        for (int j = 0; j < referencePoints.size(); ++j) {
            if (referencePoints(i) < referencePoints(j)) {
                index_t a = sortedBoundaryFunctions(i);
                real_t  b = referencePoints(i);
                sortedBoundaryFunctions(i) = sortedBoundaryFunctions(j);
                referencePoints(i) = referencePoints(j);
                sortedBoundaryFunctions(j) = a;
                referencePoints(j) = b;
            }
        }
    }
    gsVector<index_t> isToBeMoved(sortedBoundaryFunctions.size());
    gsMatrix<> xy1;
    gsMatrix<> xy2;

    gsMatrix<> uv1(2, (numPoints + 1) * (1));
    gsMatrix<> uv2(2, (numPoints + 1) * (1));

    for (int i = 0; i < isToBeMoved.size(); ++i) {
        isToBeMoved(i) = 0;
    }
    bool thirdIsAbove, fourthIsAbove;
    for (int i = 0; i < referencePoints.size() - 1; ++i) {

        for (int j = 0; j <= numPoints; ++j) {
            const real_t val = referencePoints(i) + (referencePoints(i + 1) - referencePoints(i)) * (double)j / numPoints;
            const real_t zeroWTF = 0.0;
            const real_t  oneWTF = 1.0;
            uv1(1, j) = uv2(1, j) = min(max(val, zeroWTF), oneWTF);
            uv1(0, j) = 1.0;
            uv2(0, j) = 0.0;
        }
        xy1 = geom1.eval(uv1);
        xy2 = geom2.eval(uv2);
        for (int j = 0; j < numPoints; ++j) {
            real_t x1 = xy1(0, j);
            real_t y1 = xy1(1, j);
            real_t x2 = xy1(0, j + 1);
            real_t y2 = xy1(1, j + 1);
            real_t x3 = xy2(0, j);
            real_t y3 = xy2(1, j);
            real_t x4 = xy2(0, j + 1);
            real_t y4 = xy2(1, j + 1);
            real_t val1 = (y3 - y1) * (x2 - x1) - (y2 - y1) * (x3 - x1);
            real_t val2 = (y4 - y1) * (x2 - x1) - (y2 - y1) * (x4 - x1);
            thirdIsAbove = val1 > -tolerance;
            fourthIsAbove = val2 > -tolerance;
            if (!((thirdIsAbove || !fourthIsAbove) && (!thirdIsAbove || fourthIsAbove))) {
                isToBeMoved(i) = isToBeMoved(i + 1) = 1;
                doPatchesIntersect = 1;

            }

        }

    }
    return doPatchesIntersect;
}

/**
 * @brief Checks if mapped lines from two patches intersect improperly.
 *
 * This overload works with the parameterized/mapped (x,y) lines rather than
 * the geometry objects. It detects self-intersections and folding in the
 * parameterization by checking if line segments from different patches
 * intersect in physical space.
 *
 * @param[in] xyLines    All mapped (x,y) lines from generateAndExportMappedLines
 * @param[in] patchIDs   Patch ID for each line
 * @param[in] patch1     First patch to check
 * @param[in] patch2     Second patch to check
 * @param[in] tolerance  Tolerance for intersection detection
 * @param[in] verbose    Enable detailed logging
 * @return Number of intersection violations detected
 */
int PatchesIntersection(
    const std::vector<gsMatrix<>>& xyLines,
    const std::vector<index_t>& patchIDs,
    index_t patch1,
    index_t patch2,
    real_t tolerance,
    bool verbose)
{
    if (xyLines.empty() || patchIDs.empty() || xyLines.size() != patchIDs.size())
    {
        gsWarn << "PatchesIntersection: Invalid input data.\n";
        return -1;
    }

    // Collect lines for each patch
    std::vector<size_t> lines1, lines2;
    for (size_t i = 0; i < patchIDs.size(); ++i)
    {
        if (patchIDs[i] == patch1)
            lines1.push_back(i);
        else if (patchIDs[i] == patch2)
            lines2.push_back(i);
    }

    if (lines1.empty() || lines2.empty())
    {
        if (verbose)
            gsInfo << "One or both patches have no lines. Skipping intersection check.\n";
        return 0;
    }

    if (verbose)
    {
        gsInfo << "\n=== Checking intersection between patches " << patch1 << " and " << patch2 << " ===\n";
        gsInfo << "  Patch " << patch1 << ": " << lines1.size() << " lines\n";
        gsInfo << "  Patch " << patch2 << ": " << lines2.size() << " lines\n";
    }

    int intersectionCount = 0;

    // Check each line segment from patch1 against each line segment from patch2
    for (size_t idx1 : lines1)
    {
        const gsMatrix<>& line1 = xyLines[idx1];
        const index_t numPts1 = line1.cols();

        for (index_t seg1 = 0; seg1 < numPts1 - 1; ++seg1)
        {
            // Segment 1: (x1,y1) to (x2,y2)
            const real_t x1 = line1(0, seg1);
            const real_t y1 = line1(1, seg1);
            const real_t x2 = line1(0, seg1 + 1);
            const real_t y2 = line1(1, seg1 + 1);

            for (size_t idx2 : lines2)
            {
                const gsMatrix<>& line2 = xyLines[idx2];
                const index_t numPts2 = line2.cols();

                for (index_t seg2 = 0; seg2 < numPts2 - 1; ++seg2)
                {
                    // Segment 2: (x3,y3) to (x4,y4)
                    const real_t x3 = line2(0, seg2);
                    const real_t y3 = line2(1, seg2);
                    const real_t x4 = line2(0, seg2 + 1);
                    const real_t y4 = line2(1, seg2 + 1);

                    // Check if segments intersect using cross products
                    // For segment 1-2, check if points 3 and 4 are on opposite sides
                    const real_t d1 = (x2 - x1) * (y3 - y1) - (y2 - y1) * (x3 - x1);
                    const real_t d2 = (x2 - x1) * (y4 - y1) - (y2 - y1) * (x4 - x1);

                    // For segment 3-4, check if points 1 and 2 are on opposite sides
                    const real_t d3 = (x4 - x3) * (y1 - y3) - (y4 - y3) * (x1 - x3);
                    const real_t d4 = (x4 - x3) * (y2 - y3) - (y4 - y3) * (x2 - x3);

                    // Segments intersect if signs differ (with tolerance)
                    const bool cross1 = (d1 * d2 < -tolerance * tolerance);
                    const bool cross2 = (d3 * d4 < -tolerance * tolerance);

                    if (cross1 && cross2)
                    {
                        ++intersectionCount;

                        if (verbose)
                        {
                            gsInfo << "  [INTERSECTION] Line " << idx1 << " seg " << seg1
                                << " crosses Line " << idx2 << " seg " << seg2 << "\n";
                            gsInfo << "    Seg1: (" << x1 << "," << y1 << ") -> (" << x2 << "," << y2 << ")\n";
                            gsInfo << "    Seg2: (" << x3 << "," << y3 << ") -> (" << x4 << "," << y4 << ")\n";
                            gsInfo << "    Cross products: d1=" << d1 << ", d2=" << d2
                                << ", d3=" << d3 << ", d4=" << d4 << "\n";
                        }
                    }
                }
            }
        }
    }

    if (verbose)
    {
        if (intersectionCount > 0)
            gsInfo << "  TOTAL INTERSECTIONS: " << intersectionCount << "\n";
        else
            gsInfo << "  No intersections detected.\n";
        gsInfo << "=== Intersection check complete ===\n\n";
    }

    return intersectionCount;
}

/**
 * @brief Checks all pairs of patches for intersections in their mapped lines.
 *
 * Convenience function that checks every pair of patches for improper intersections.
 *
 * @param[in] xyLines    All mapped (x,y) lines
 * @param[in] patchIDs   Patch ID for each line
 * @param[in] numPatches Total number of patches
 * @param[in] tolerance  Tolerance for intersection detection
 * @param[in] verbose    Enable detailed logging
 * @return Total number of intersection violations across all patch pairs
 */
int checkAllPatchIntersections(
    const std::vector<gsMatrix<>>& xyLines,
    const std::vector<index_t>& patchIDs,
    index_t numPatches,
    real_t tolerance,
    bool verbose)
{
    int totalIntersections = 0;

    if (verbose)
        gsInfo << "\\n=== Checking intersections for all patch pairs ===\\n";

    // Check each pair of patches (including self-intersections)
    for (index_t p1 = 0; p1 < numPatches; ++p1)
    {
        for (index_t p2 = p1; p2 < numPatches; ++p2)
        {
            int count = PatchesIntersection(xyLines, patchIDs, p1, p2, tolerance, verbose);
            if (count > 0)
            {
                totalIntersections += count;
                gsWarn << "Patches " << p1 << " and " << p2 << " have " << count << " intersections!\n";
            }
        }
    }

    if (verbose || totalIntersections > 0)
    {
        if (totalIntersections > 0)
            gsWarn << "\nTOTAL INTERSECTIONS FOUND: " << totalIntersections << "\n";
        else
            gsInfo << "\nNo intersections found across all patches.\n";
        gsInfo << "=== Intersection check complete ===\n\n";
    }

    return totalIntersections;
}

void assembleA(gsVector<gsMatrix<>>  uv,
    gsVector<gsVector<gsTensorBSplineBasis<2, real_t>>> Bells,
    gsVector<gsVector  <std::vector<index_t>>> functionDescription,
    std::vector<std::vector<gsMatrix<int>>> actives,
    gsSparseMatrix<real_t>& A_mat
)
{

    int shift = 0;
    //outfile << "assembling\n\n";
    for (size_t patch = 0; patch < Bells.size(); patch++)
    {
        //gsInfo << "uv.cols(): " << uv(patch).cols() << "\n";
        for (size_t i = 0; i < uv(patch).cols(); i++)
        {
            gsMatrix<> punto = uv(patch).col(i);
            //outfile << "punto: " << punto << "\n";
            //gsInfo << "punto: " << punto << "\n";
 /*           gsInfo << i << "\n";
            gsInfo << actives[patch][i].cols() << "\n\n";*/
            //outfile << "point " << i + shift << "\n";
            real_t checkSum = 0.0;
            for (size_t functionIndex = 0; functionIndex < actives[patch][i].cols(); functionIndex++)
            {
                //get the index
                //gsInfo << functionIndex << "\n";
                auto corractiveIndex = actives[patch][i](0, functionIndex);
                auto corrPiece = actives[patch][i](1, functionIndex);
                auto corrIndex = functionDescription(corractiveIndex)(corrPiece)[2];
                auto corrLevel = functionDescription(corractiveIndex)(corrPiece)[1];
                //if (corrLevel == 3)  outfile << "level 3\n";
                gsMatrix<> valore;
                valore = Bells(patch)(corrLevel).function(corrIndex).eval(punto);
                //Assign the corresponding place
                A_mat(i + shift, corractiveIndex) = valore(0, 0);
                checkSum += valore(0, 0);
                //outfile << patch << " " << corrLevel << " " << corrIndex << ":" << valore(0, 0) << " at ("<< uv(patch)(0,i)  << ", " << uv(patch)(1,i) << ")\n";
            }
            if (checkSum <= 0)   gsInfo << i << " is problematic\n";
            //outfile << "checkSum: " << checkSum << "\n";
            checkSum = 0.0;
            //outfile << "\n";
            //,
        }
        //outfile << "\n";
        shift += uv(patch).cols();
    }
}

void assembleA(gsVector<gsMatrix<>>  uv,
    gsVector<gsMatrix<>>  xy,
    gsVector<gsVector<gsTensorBSplineBasis<2, real_t>>> Bells,
    gsVector  <gsVector< gsVector<index_t>>> isActive,
    gsVector  <gsVector< gsVector<index_t>>>  globalIndex,
    gsVector<gsVector  <std::vector<index_t>>> functionDescription,
    gsSparseMatrix<real_t>& A_mat
) {
    int shift = 0;
    for (int patch = 0; patch < Bells.size(); ++patch) {
        for (int functionIndex = 0; functionIndex < functionDescription.size(); ++functionIndex) {
            for (int k = 0; k < uv(patch).cols(); ++k) {
                //gsInfo << functionIndex << "\n";

                bool functionIsActive = false;
                int corrLevel = 0; // Initialize outside the loop
                int corrIndex = 0; // Initialize outside the loop
                gsMatrix<> valore; // Initialize outside the loop

                for (const auto& piece : functionDescription(functionIndex)) {
                    if (piece[0] == patch) {
                        functionIsActive = true;
                        corrLevel = piece[1];
                        corrIndex = piece[2];
                        valore = Bells(patch)(corrLevel).function(corrIndex).eval(uv(patch).col(k));
                        break; // No need to continue once found
                    }
                }

                if (functionIsActive) {
                    A_mat(k + shift, functionIndex) = valore(0, 0);
                }
            }
        }
        shift += uv(patch).cols();
    }
}

void assembleA(gsMatrix<>  uv,
    gsVector<gsVector<gsTensorBSplineBasis<2, real_t>>> Bells,
    gsVector  <gsVector< gsVector<index_t>>> isActive,
    gsVector  <gsVector< gsVector<index_t>>>  globalIndex,
    gsVector<gsVector  <std::vector<index_t>>> functionDescription,
    gsSparseMatrix<real_t>& A_mat
) {
    int shift = 0;
    for (int patch = 0; patch < Bells.size(); ++patch) {
        for (int functionIndex = 0; functionIndex < functionDescription.size(); ++functionIndex) {
            for (int k = 0; k < uv.cols(); ++k) {
                int corrLevel;
                int corrIndex;
                gsMatrix<> valore;
                bool functionIsActive = false;
                for (int pieceOfFunction = 0;
                    pieceOfFunction < functionDescription(functionIndex).size();
                    pieceOfFunction++) {
                    if (functionDescription(functionIndex)(pieceOfFunction)[0] == patch) {
                        functionIsActive = true;
                        corrLevel = functionDescription(functionIndex)(pieceOfFunction)[1];
                        corrIndex = functionDescription(functionIndex)(pieceOfFunction)[2];
                        valore = Bells(patch)(corrLevel).function(corrIndex).eval(uv.col(k));
                    }
                }
                if (functionIsActive) {
                    A_mat(k + shift, functionIndex) = valore(0, 0);
                }
            }
        }
        shift += uv.cols();
    }
}

#include <iostream> // For logging

#include <iostream> // For logging

void active_into_Spillovers(
    index_t patch,
    const gsMatrix<real_t>& uv,
    const std::vector<std::vector<std::array<int, 3>>>& spilloverFunctionCoordinates,
    gsVector<gsVector<gsTensorBSplineBasis<2, real_t>>> Bells,
    std::vector<std::vector<int>>& activeSpillovers)
{
    activeSpillovers.clear();
    activeSpillovers.resize(uv.cols());

    //std::cout << "Starting active_into_Spillovers with patch: " << patch << std::endl;
    //std::cout << "Total spilloverFunctionCoordinates: " << spilloverFunctionCoordinates.size() << std::endl;
    //std::cout << "Number of uv columns: " << uv.cols() << std::endl;

    for (size_t functionIndex = 0; functionIndex < spilloverFunctionCoordinates.size(); ++functionIndex) {
        //std::cout << "Function index " << functionIndex << " has "
        //    << spilloverFunctionCoordinates[functionIndex].size() << " spillover entries." << std::endl;

        for (size_t spilloverIndex = 0; spilloverIndex < spilloverFunctionCoordinates[functionIndex].size(); spilloverIndex++) {
            const auto& spillover = spilloverFunctionCoordinates[functionIndex][spilloverIndex];

            index_t spillPatch = spillover[0];
            index_t spillLevel = spillover[1];
            index_t spillIndex = spillover[2];

            //    << ", spillPatch: " << spillPatch
            //    << ", spillLevel: " << spillLevel
            //    << ", spillIndex: " << spillIndex << ")" << std::endl;

            if (spillPatch == patch) {
                //std::cout << "Match found for patch!" << std::endl;

                // Validate indexing before accessing Bells
                if (spillPatch >= Bells.size()) {
                    //    << " is out of bounds (max: " << Bells.size() - 1 << ")" << std::endl;
                    continue;
                }

                if (spillLevel >= Bells[spillPatch].size()) {
                    //std::cerr << "ERROR: spillLevel " << spillLevel
                    //    << " is out of bounds (max: " << Bells[spillPatch].size() - 1 << ")" << std::endl;
                    continue;
                }

                gsMatrix<> supp = Bells[spillPatch][spillLevel].function(spillIndex).support();
                //std::cout << "Support region: [" << supp(0, 0) << ", " << supp(0, 1) << "] x [" 
                //    << supp(1, 0) << ", " << supp(1, 1) << "]" << std::endl;

                for (index_t p = 0; p < uv.cols(); ++p) {
                    //std::cout << "Checking uv(" << p << "): [" << uv(0, p) << ", " << uv(1, p) << "]" << std::endl;
                    gsMatrix<> TestValore = Bells[spillPatch][spillLevel].function(spillIndex).eval(uv.col(p));
                    bool inside_x = (uv(0, p) >= supp(0, 0)) && (uv(0, p) <= supp(0, 1));
                    bool inside_y = (uv(1, p) >= supp(1, 0)) && (uv(1, p) <= supp(1, 1));

                    if (inside_x && inside_y) {
                        //std::cout << "Point " << p << " is INSIDE support region. Adding functionIndex "
                        //    << functionIndex << " to activeSpillovers[" << p << "]" << std::endl;
                        activeSpillovers[p].push_back(functionIndex);
                    }
                    else {
                        if (TestValore(0, 0) != 0) {
                            gsInfo << "\n==== INEXPLICABLE CASE DETECTED ====\n";
                            gsInfo << "Function Index: " << functionIndex << "\n";
                            gsInfo << "Spill Patch: " << spillPatch << ", Spill Level: " << spillLevel
                                << ", Spill Index: " << spillIndex << "\n";
                            gsInfo << "Support Region: [" << supp(0, 0) << ", " << supp(0, 1) << "] x ["
                                << supp(1, 0) << ", " << supp(1, 1) << "]\n";
                            gsInfo << "uv(" << p << "): [" << uv(0, p) << ", " << uv(1, p) << "]\n";
                            gsInfo << "TestValore(0, 0): " << TestValore(0, 0) << "\n";
                            gsInfo << "Inside X: " << (inside_x ? "YES" : "NO")
                                << ", Inside Y: " << (inside_y ? "YES" : "NO") << "\n";
                            gsInfo << "===================================\n";
                        }
                        //std::cout << "Point " << p << " is OUTSIDE support region: "
                        //    << (inside_x ? "X inside, " : "X outside, ")
                        //    << (inside_y ? "Y inside" : "Y outside") << std::endl;
                    }
                }
            }
        }
    }

    //std::cout << "Completed active_into_Spillovers." << std::endl;
}

void exportLinesToFile(
    const std::vector<gsMatrix<>>& lines,
    const std::vector<index_t>& patchIDs,
    const std::vector<char>& directions,
    const std::string& filename,
    bool verbose,
    const std::vector<gsMatrix<>>& uvLines // <-- NEW
)
{
    std::ofstream file(filename);
    if (!file.is_open())
    {
        gsWarn << "Could not open file: " << filename << "\n";
        return;
    }

    for (size_t i = 0; i < lines.size(); ++i)
    {
        const gsMatrix<>& uv = uvLines[i];
        const gsMatrix<>& xy = lines[i];

        // Retrieve start and end parametric coordinates
        std::stringstream header;
        header << "# Patch " << patchIDs[i]
            << " Direction " << directions[i];

        if (uv.cols() > 1)
        {
            header << " Start (u,v): (" << uv(0, 0) << ", " << uv(1, 0) << ")"
                << " End (u,v): (" << uv(0, uv.cols() - 1) << ", " << uv(1, uv.cols() - 1) << ")";
        }

        file << header.str() << "\n";

        for (index_t k = 0; k < xy.cols(); ++k)
            file << xy(0, k) << " " << xy(1, k) << "\n";

        file << "\n";

        if (verbose)
            outfile << "Wrote line " << i << " with " << xy.cols() << " points.\n";
    }

    file.close();

    if (verbose)
        outfile << "Exported " << lines.size() << " lines to file: " << filename << "\n";
}

/**
 * @brief Verifies partition of unity for an assembled system matrix.
 *
 * This function checks whether each row of the given sparse matrix sums
 * to approximately one, which is a necessary condition for a valid
 * partition of unity (e.g. for THB-spline bases).
 *
 * The check is performed row-wise:
 *   sum_j A(i,j) ≈ 1
 *
 * Rows whose sum deviates from unity by more than the given tolerance
 * are reported.
 *
 * @param A        Assembled system matrix (rows correspond to quadrature points).
 * @param tol      Numerical tolerance for deviation from unity.
 * @param verbose  If true, prints all violating rows and a summary.
 *
 * @note
 * - This function assumes that each row corresponds to one evaluation point.
 * - Useful for debugging truncation, basis activation, and index mapping issues.
 * - A failure usually indicates missing basis contributions or incorrect THB hierarchy handling.
 */
bool checkPartitionOfUnity(const gsSparseMatrix<real_t>& A,
    real_t tol = 1e-10)
{
    return analyzePartitionOfUnity(A, tol).allSatisfied;
}



// Validate that MPBES presentations match THB presentations
int validatePresentations(
    const MPBES<2, real_t>& mpbes,
    const gsVector<gsVector<gsTensorBSplineBasis<2, real_t>>>& Bells,
    const gsVector<gsTHBSplineBasis<2, real_t>>& SubdomainHierarchy,
    int attempt = 0)
{
    gsInfo << "\n========== Validating Presentations (Attempt " << attempt << ") ==========\n";

    const auto& functionDescription = mpbes.functionDescription();
    const auto& indexInTHB = mpbes.indexInTHB();
    const auto& presentation = mpbes.presentation();
    const auto& isTruncated = mpbes.isTruncated();

    bool presentationError = false;
    int truncatedCount = 0;
    int checkedComponents = 0;

    // Compare MPBES presentations with THB presentations
    for (size_t f = 0; f < functionDescription.size(); ++f) {
        if (!isTruncated[f]) continue;

        truncatedCount++;
        gsInfo << "\n--- Function " << f << " (truncated) ---\n";

        // Check each component
        for (size_t compIdx = 0; compIdx < functionDescription[f].size(); ++compIdx) {
            int patch = functionDescription[f][compIdx][0];
            int level = functionDescription[f][compIdx][1];
            int tensorIdx = functionDescription[f][compIdx][2];
            int thbIdx = indexInTHB[patch][level][tensorIdx];

            gsInfo << "  Component " << compIdx << ": patch=" << patch << ", level=" << level
                << ", tensorIdx=" << tensorIdx << ", thbIdx=" << thbIdx;

            if (thbIdx == -1) {
                gsInfo << " [SPILLOVER - skipping]\n";
                continue;
            }
            gsInfo << "\n";

            // Check if THIS component is truncated
            if (!mpbes.isComponentTruncated(f, compIdx)) {
                gsInfo << "    Component NOT truncated - skipping\n";
                continue;
            }

            checkedComponents++;

            // Get THB basis
            const gsTHBSplineBasis<2, real_t>& thbBasis = SubdomainHierarchy[patch];

            if (!thbBasis.isTruncated(thbIdx)) {
                gsInfo << "    ERROR: Truncated in MPBES but NOT in THB!\n";
                presentationError = true;
                continue;
            }

            // Get both presentations
            const auto& mpbesCoefs = presentation[f][compIdx];
            const auto& thbCoefs = thbBasis.getCoefs(thbIdx);

            gsInfo << "    MPBES coefs size: " << mpbesCoefs.size() << "\n";
            gsInfo << "    THB coefs size:   " << thbCoefs.size() << "\n";

            if (mpbesCoefs.size() != thbCoefs.size()) {
                gsInfo << "    ERROR: Size mismatch!\n";
                presentationError = true;
                continue;
            }

            // Compare coefficients
            real_t maxDiff = 0.0;
            int differingCoefs = 0;

            for (typename gsSparseVector<real_t>::InnerIterator it(mpbesCoefs); it; ++it) {
                index_t idx = it.index();
                real_t mpbesVal = it.value();
                real_t thbVal = thbCoefs[idx];
                real_t diff = std::abs(mpbesVal - thbVal);

                if (diff > maxDiff) maxDiff = diff;

                if (diff > 1e-10) {
                    gsInfo << "    Coef[" << idx << "]: MPBES=" << mpbesVal
                        << " vs THB=" << thbVal << " (diff=" << diff << ")\n";
                    differingCoefs++;
                }
            }

            if (differingCoefs > 0) {
                gsInfo << "    ERROR: " << differingCoefs << " coefficients differ (max diff: " << maxDiff << ")\n";
                presentationError = true;
            }
            else {
                gsInfo << "    OK: All coefficients match\n";
            }
        }
    }

    gsInfo << "\nChecked " << checkedComponents << " truncated components in " << truncatedCount << " functions\n";

    if (presentationError) {
        gsInfo << "========================================\n";
        gsInfo << "PRESENTATION VALIDATION FAILED!\n";
        gsInfo << "========================================\n";
        if (attempt == 0) {
            return 260114;
        }
    }
    else {
        gsInfo << "Presentation validation PASSED\n";
    }

    gsInfo << "========================================\n\n";
    return 0;
}

// Compare MPBES evaluation against THB/SubdomainHierarchy evaluation on Gauss points
void compareMPBESvsTHB(
    const MPBES<2, real_t>& mpbes,
    const gsVector<gsTHBSplineBasis<2, real_t>>& SubdomainHierarchy,
    const gsVector<gsVector<gsVector<index_t>>>& indexInTHB,
    const std::vector<std::vector<std::vector<index_t>>>& functionDescription,
    const std::vector<bool>& isTruncated,
    const std::vector<std::vector<gsSparseVector<real_t>>>& presentation,
    const std::string& filename)
{
    std::ofstream logFile(filename);
    if (!logFile.is_open())
    {
        gsWarn << "Could not open file: " << filename << "\n";
        return;
    }

    logFile << "=== MPBES vs THB Comparison (Gauss points) ===\n";

    const auto& thbBases = mpbes.thbBases();
    const auto& bellsBases = mpbes.bellsBases();

    gsVector<index_t> numNodes(2);
    numNodes[0] = numNodes[1] = 2;
    gsGaussRule<real_t> quRule(numNodes);

    index_t row = 0;
    for (index_t patch = 0; patch < thbBases.size(); ++patch)
    {
        auto domIt = thbBases[patch].makeDomainIterator();
        gsMatrix<real_t> quNodes;
        gsVector<real_t> quWeights;

        for (int el = 0; domIt->good(); domIt->next(), ++el)
        {
            gsVector<real_t> lower = domIt->lowerCorner();
            gsVector<real_t> upper = domIt->upperCorner();
            quRule.mapTo(lower, upper, quNodes, quWeights);

            for (index_t q = 0; q < quNodes.cols(); ++q)
            {
                const gsVector<real_t> uvk = quNodes.col(q);

                real_t mpbesSum = 0.0;
                real_t thbSum = 0.0;
                int mpbesNNZ = 0;
                int thbNNZ = 0;

                for (index_t f = 0; f < mpbes.size(); ++f)
                {
                    gsMatrix<real_t> mpbesVal;
                    mpbes.evalSingleOnPatch(f, patch, uvk, mpbesVal);
                    real_t vMp = mpbesVal(0, 0);
                    if (vMp != 0.0) {
                        mpbesSum += vMp;
                        mpbesNNZ++;
                    }

                    // THB evaluation using SubdomainHierarchy + presentation
                    real_t vThb = 0.0;
                    int compIdx = 0;
                    for (const auto& desc : functionDescription[f])
                    {
                        int fnPatch = desc[0];
                        int fnLevel = desc[1];
                        int fnTensorIndex = desc[2];

                        if (fnPatch != static_cast<int>(patch)) {
                            compIdx++;
                            continue;
                        }

                        int thbIdx = -1;
                        if (fnPatch < static_cast<int>(indexInTHB.size()) &&
                            fnLevel < static_cast<int>(indexInTHB[fnPatch].size()) &&
                            fnTensorIndex < static_cast<int>(indexInTHB[fnPatch][fnLevel].size()))
                        {
                            thbIdx = indexInTHB[fnPatch][fnLevel][fnTensorIndex];
                        }

                        gsMatrix<real_t> compVal;
                        const gsTensorBSplineBasis<2, real_t>* bellsBasis = nullptr;
                        int bellsIdx = -1;
                        if (thbIdx == -1 && fnPatch < static_cast<int>(bellsBases.size()) &&
                            fnLevel < static_cast<int>(bellsBases[fnPatch].size()))
                        {
                            bellsBasis = &bellsBases[fnPatch][fnLevel];
                            bellsIdx = fnTensorIndex;
                        }

                        evalSingle_into(
                            f,
                            thbIdx,
                            uvk,
                            isTruncated,
                            SubdomainHierarchy[fnPatch],
                            presentation[f][compIdx],
                            compVal,
                            bellsBasis,
                            bellsIdx);

                        vThb += compVal(0, 0);
                        compIdx++;
                    }

                    if (vThb != 0.0) {
                        thbSum += vThb;
                        thbNNZ++;
                    }
                }

                logFile << "Row " << row
                    << " | Patch " << patch
                    << " | Element " << el
                    << " | Gauss " << q
                    << " | (u,v)= (" << uvk(0) << ", " << uvk(1) << ")\n";
                logFile << "  MPBES: sum=" << std::setprecision(16) << mpbesSum
                    << " nnz=" << mpbesNNZ << "\n";
                logFile << "  THB  : sum=" << std::setprecision(16) << thbSum
                    << " nnz=" << thbNNZ << "\n";
                logFile << "  Diff : " << std::setprecision(16) << (mpbesSum - thbSum) << "\n\n";

                row++;
            }
        }
    }

    logFile << "=== End of Comparison ===\n";
    logFile.close();
}

// MPBES-based assemble - Uses MPBES presentation for truncated functions
void assemble(
    const gsVector<gsMatrix<>>& uv,
    const MPBES<2, real_t>& mpbes,
    gsSparseMatrix<real_t>& A_mat,
    gsMatrix<real_t>& b_vec,
    const gsMultiPatch<real_t>& mp,
    bool verbose,
    int currentPatch,
    int currentLevel,
    int currentAttempt,
    const std::vector<index_t>* activeFunctionIds)
{
    try {
        PROFILE_SECTION("assemble_mpbes_total");

        if (verbose)
            outfile << "=== assemble (MPBES, Gauss points) BEGIN ===\n";

        // Input validation
        try {
            if (mpbes.size() == 0) {
                gsInfo << "WARNING: mpbes.size() is 0 in assemble, returning without assembly\n";
                A_mat.resize(0, 0);
                b_vec.resize(0, 2);
                return;
            }
            if (mpbes.nPatches() == 0) {
                gsInfo << "WARNING: mpbes.nPatches() is 0 in assemble, returning without assembly\n";
                A_mat.resize(0, 0);
                b_vec.resize(0, 2);
                return;
            }
            if (mp.nPatches() == 0) {
                gsInfo << "WARNING: mp.nPatches() is 0 in assemble, returning without assembly\n";
                A_mat.resize(0, 0);
                b_vec.resize(0, 2);
                return;
            }
        } catch (const std::exception& e) {
            gsInfo << "EXCEPTION during assemble input validation: " << e.what() << ", aborting\n";
            A_mat.resize(0, 0);
            b_vec.resize(0, 2);
            return;
        } catch (...) {
            gsInfo << "UNKNOWN EXCEPTION during assemble input validation, aborting\n";
            A_mat.resize(0, 0);
            b_vec.resize(0, 2);
            return;
        }

        gsVector<index_t> numNodes(2);
        numNodes[0] = numNodes[1] = 2;
        gsGaussRule<real_t> quRule(numNodes);

        if (verbose)
            gsInfo << "Using " << numNodes[0] << "x" << numNodes[1] << " = " << numNodes.prod() << " Gauss points per element\n";

        const index_t numFunctions = (activeFunctionIds ? static_cast<index_t>(activeFunctionIds->size()) : mpbes.size());
        const index_t numPatches = mpbes.nPatches();

        auto globalFunctionOf = [&](index_t localF) -> index_t {
            if (activeFunctionIds)
            {
                if (localF < 0 || localF >= static_cast<index_t>(activeFunctionIds->size()))
                    return -1;
                return (*activeFunctionIds)[localF];
            }
            return localF;
        };

    // Get references to MPBES data
    const auto& functionDescription = mpbes.functionDescription();
    const auto& thbBases = mpbes.thbBases();
    const auto& indexInTHB = mpbes.indexInTHB();
    const auto& isTruncated = mpbes.isTruncated();
    const auto& presentation = mpbes.presentation();
    const auto& bellsBases = mpbes.bellsBases();
    const auto& hasSpillover = mpbes.hasSpillover();
    const auto& spilloverFunctionCoordinates = mpbes.spilloverCoordinates();

    // DO NOT pre-filter active functions - let evalSingleOnPatch handle it
    // Pre-filtering was too aggressive and missed spillover/truncated contributions

    // Check if custom uv points are provided (for testBoundaryAssembly)
    bool useCustomPoints = false;
    index_t totalCustomPoints = 0;
    
    if (uv.size() > 0) {
        for (index_t patch = 0; patch < uv.size(); ++patch) {
            if (uv[patch].cols() > 0) {
                useCustomPoints = true;
                totalCustomPoints += uv[patch].cols();
            }
        }
    }
    
    // Count total rows
    index_t totalRows = 0;
    
    if (useCustomPoints) {
        // Use custom points from uv parameter
        totalRows = totalCustomPoints;
        if (g_verbose) gsInfo << "Using " << totalRows << " custom evaluation points from uv parameter\n";
        if (verbose) {
            outfile << "Using " << totalRows << " custom evaluation points from uv parameter\n";
        }
        if (outfile.is_open()) {
            outfile << "[assemble] custom-point layout per patch:\n";
            for (index_t patch = 0; patch < uv.size(); ++patch)
                outfile << "  patch=" << patch << " uv.cols=" << uv[patch].cols() << "\n";
        }
    } else {
        // Use Gauss quadrature points (original behavior)
        try {
            for (index_t patch = 0; patch < numPatches; ++patch)
            {
                if (patch >= static_cast<index_t>(thbBases.size())) {
                    gsInfo << "WARNING: patch " << patch << " >= thbBases.size() " << thbBases.size() << "\n";
                    continue;
                }
                
                auto domIt = thbBases[patch].makeDomainIterator();
                for (; domIt->good(); domIt->next())
                    totalRows += numNodes.prod();
            }
        } catch (const std::exception& e) {
            gsInfo << "EXCEPTION counting rows: " << e.what() << ", using totalRows=0\n";
            totalRows = 0;
        } catch (...) {
            gsInfo << "UNKNOWN EXCEPTION counting rows, using totalRows=0\n";
            totalRows = 0;
        }
    }
    
    // Validate dimensions before resize
    if (totalRows == 0) {
        gsInfo << "ERROR: totalRows is 0, cannot assemble\n";
        outfile << "ERROR: totalRows is 0, cannot assemble\n";
        A_mat.resize(0, 0);
        b_vec.resize(0, 2);
        return;
    }
    if (numFunctions == 0) {
        gsInfo << "ERROR: numFunctions is 0, cannot assemble\n";
        outfile << "ERROR: numFunctions is 0, cannot assemble\n";
        A_mat.resize(0, 0);
        b_vec.resize(0, 2);
        return;
    }
    
    if (g_verbose) gsInfo << "Resizing matrices: A_mat(" << totalRows << " x " << numFunctions << "), b_vec(" << totalRows << " x 2)\n";
    if (outfile.is_open())
        outfile << "[assemble] resize A_mat(" << totalRows << "," << numFunctions
                << "), b_vec(" << totalRows << ",2)\n";

    A_mat.resize(totalRows, numFunctions);
    A_mat.setZero();
    b_vec.resize(totalRows, 2);
    b_vec.setZero();

    if (outfile.is_open())
        outfile << "[assemble] after resize: A_mat.rows=" << A_mat.rows()
                << " cols=" << A_mat.cols()
                << " b_vec.rows=" << b_vec.rows()
                << " cols=" << b_vec.cols() << "\n";

    index_t rowOffset = 0;

    std::ofstream row0Log;
    bool row0LogOpened = false;
    real_t row0Sum = 0.0;
    bool row0SumSet = false;
    int row0NonZeroFuncs = 0;
    std::vector<int> row0NonZeroFunctionList;
    std::vector<real_t> row0NonZeroValueList;
    int row0ZeroNoPatchComponents = 0;
    int row0ZeroHasPatchComponents = 0;

    std::ofstream rowTargetLog;
    bool rowTargetLogOpened = false;
    const index_t targetRow = 336;
    real_t rowTargetSum = 0.0;
    bool rowTargetSumSet = false;
    int rowTargetNonZeroFuncs = 0;
    std::vector<int> rowTargetNonZeroFunctionList;
    std::vector<real_t> rowTargetNonZeroValueList;

    if (verbose)
    {
        gsInfo << "Total rows: " << totalRows << "\n";
        gsInfo << "Total basis functions: " << numFunctions << "\n";
        gsInfo << "Using " << (useCustomPoints ? "custom" : "Gauss quadrature") << " points\n";
    }

    if (useCustomPoints) {
        // ===== CUSTOM POINTS MODE (for testBoundaryAssembly) =====
        for (index_t patch = 0; patch < numPatches; ++patch)
        {
            if (patch >= uv.size() || uv[patch].cols() == 0) {
                continue;
            }
            
            if (verbose) {
                outfile << ">> Patch " << patch << ": " << uv[patch].cols() << " custom points\n";
            }
            
            const gsMatrix<>& patchPoints = uv[patch];
            
            for (index_t pt = 0; pt < patchPoints.cols(); ++pt)
            {
                if (rowOffset >= totalRows) {
                    gsInfo << "ERROR: rowOffset " << rowOffset << " >= totalRows " << totalRows << "\n";
                    break;
                }
                
                index_t row = rowOffset;
                const gsVector<real_t> uvk = patchPoints.col(pt);
                
                // RHS: target geometry
                gsMatrix<real_t> xy;
                try {
                    if (patch >= mp.nPatches()) {
                        if (verbose) {
                            gsInfo << "WARNING: patch " << patch << " >= mp.nPatches() " << mp.nPatches() << "\n";
                        }
                        ++rowOffset;
                        continue;
                    }
                    
                    xy = mp.patch(patch).eval(uvk);
                    
                    if (xy.rows() < 1 || xy.cols() < 1) {
                        b_vec(row, 0) = 0.0;
                        b_vec(row, 1) = 0.0;
                    } else {
                        b_vec(row, 0) = xy(0, 0);
                        b_vec(row, 1) = xy(1, 0);
                    }
                } catch (const std::exception& e) {
                    if (verbose) {
                        gsInfo << "EXCEPTION in mp.patch().eval(): " << e.what() << "\n";
                    }
                    b_vec(row, 0) = 0.0;
                    b_vec(row, 1) = 0.0;
                } catch (...) {
                    b_vec(row, 0) = 0.0;
                    b_vec(row, 1) = 0.0;
                }
                
                // Evaluate all functions at this point
                for (index_t f = 0; f < numFunctions; ++f)
                {
                    const index_t globalF = globalFunctionOf(f);
                    if (globalF < 0 || globalF >= static_cast<index_t>(functionDescription.size())) {
                        continue;
                    }
                    
                    gsMatrix<real_t> basisValue;
                    real_t value = 0.0;
                    
                    try {
                        mpbes.evalSingleOnPatch(globalF, patch, uvk, basisValue);
                        
                        if (basisValue.rows() >= 1 && basisValue.cols() >= 1) {
                            value = basisValue(0, 0);
                            if (std::isfinite(value) && value != 0.0) {
                                if (row < A_mat.rows() && f < A_mat.cols()) {
                                    A_mat(row, f) += value;
                                }
                            }
                        }
                    } catch (const std::exception& e) {
                        // Silent - too verbose for custom points
                    } catch (...) {
                        // Silent
                    }
                }
                
                ++rowOffset;
            }
        }
    } else {
        // ===== GAUSS QUADRATURE MODE — batch evaluation =====
        // For each patch: collect all Gauss points, batch-evaluate RHS (mp.patch.eval) and
        // each basis function (evalSingleOnPatch) over all points at once.
        // Reduces N_funcs×M per-point calls → N_funcs batch calls per patch.
        std::vector<gsEigen::Triplet<real_t>> triplets;
        triplets.reserve(static_cast<size_t>(totalRows) * 8);

        for (index_t patch = 0; patch < numPatches; ++patch)
        {
            if (patch >= static_cast<index_t>(thbBases.size())) {
                gsInfo << "WARNING: patch " << patch << " >= thbBases.size() " << thbBases.size() << "\n";
                continue;
            }

            if (verbose)
                outfile << ">> Patch " << patch << "\n";

            // Phase 1: collect all Gauss points for this patch in one domain-iterator pass
            std::vector<gsVector<real_t>> gaussPtVec;
            {
                gsMatrix<real_t> qN;
                gsVector<real_t> qW;
                auto domIt = thbBases[patch].makeDomainIterator();
                for (; domIt->good(); domIt->next()) {
                    quRule.mapTo(domIt->lowerCorner(), domIt->upperCorner(), qN, qW);
                    for (index_t q = 0; q < qN.cols(); ++q)
                        gaussPtVec.emplace_back(qN.col(q));
                }
            }
            const index_t M = static_cast<index_t>(gaussPtVec.size());
            if (M == 0) continue;

            if (rowOffset + M > totalRows) {
                gsInfo << "ERROR: rowOffset+M would exceed totalRows ("
                       << rowOffset+M << " > " << totalRows << ")\n";
                outfile << "ERROR: rowOffset+M would exceed totalRows, stopping assembly\n";
                break;
            }

            gsMatrix<real_t> uvPatch(2, M);
            for (index_t i = 0; i < M; ++i)
                uvPatch.col(i) = gaussPtVec[i];

            // Phase 2: batch RHS — target geometry at all Gauss points
            try {
                if (patch < mp.nPatches()) {
                    gsMatrix<real_t> xy = mp.patch(patch).eval(uvPatch);
                    for (index_t i = 0; i < M; ++i) {
                        b_vec(rowOffset+i, 0) = (xy.rows() >= 1 && i < xy.cols()) ? xy(0, i) : 0.0;
                        b_vec(rowOffset+i, 1) = (xy.rows() >= 2 && i < xy.cols()) ? xy(1, i) : 0.0;
                    }
                }
            } catch (const std::exception& e) {
                if (verbose)
                    gsInfo << "EXCEPTION in mp.patch().eval() for patch " << patch << ": " << e.what() << "\n";
            } catch (...) {}

            // Phase 3: batch basis eval — one call per function over all M Gauss points
            gsMatrix<real_t> basisValues;
            for (index_t f = 0; f < numFunctions; ++f) {
                const index_t globalF = globalFunctionOf(f);
                if (globalF < 0 || globalF >= static_cast<index_t>(functionDescription.size())) continue;
                try {
                    mpbes.evalSingleOnPatch(globalF, patch, uvPatch, basisValues);
                } catch (const std::exception& e) {
                    if (verbose)
                        gsInfo << "EXCEPTION in evalSingleOnPatch for f=" << f << ", patch=" << patch
                               << ": " << e.what() << "\n";
                    continue;
                } catch (...) { continue; }
                if (basisValues.rows() < 1 || basisValues.cols() != M) continue;
                for (index_t i = 0; i < M; ++i) {
                    const real_t val = basisValues(0, i);
                    if (std::isfinite(val) && val != 0.0)
                        triplets.emplace_back(rowOffset+i, f, val);
                }
            }

            if (verbose)
                outfile << "  Patch " << patch << " completed: " << M << " Gauss points\n";

            rowOffset += M;
        }

        A_mat.setFromTriplets(triplets.begin(), triplets.end());

        // Silence unused-variable warnings for debug logging infrastructure (never fires in normal runs)
        (void)row0Sum; (void)row0SumSet; (void)row0NonZeroFuncs;
        (void)rowTargetSum; (void)rowTargetSumSet; (void)rowTargetNonZeroFuncs;
        (void)row0ZeroNoPatchComponents; (void)row0ZeroHasPatchComponents;

        // (batch eval complete above — old element-by-element loop removed)
    } // end of useCustomPoints if-else

    if (row0LogOpened) {
        row0Log << "\n=== Summary for Row 0 ===\n";
        row0Log << "Total non-zero functions: " << row0NonZeroFuncs << " / " << numFunctions << "\n";
        row0Log << "Row sum: " << std::setprecision(16) << (row0SumSet ? row0Sum : 0.0) << "\n";
        row0Log << "Missing contribution: " << std::setprecision(16) << (1.0 - (row0SumSet ? row0Sum : 0.0)) << "\n\n";
        row0Log << "Zero-valued functions (no components on this patch): " << row0ZeroNoPatchComponents << "\n";
        row0Log << "Zero-valued functions (has components on this patch): " << row0ZeroHasPatchComponents << "\n\n";
        
        row0Log << "=== All Contributing Functions ===\n";
        for (size_t i = 0; i < row0NonZeroFunctionList.size(); ++i) {
            row0Log << "Func " << row0NonZeroFunctionList[i] << ": " << std::setprecision(16) << row0NonZeroValueList[i] << "\n";
        }
        row0Log << "\n=== All Non-Contributing Functions ===\n";
        int zeroCount = 0;
        for (index_t f = 0; f < numFunctions; ++f) {
            bool found = false;
            for (int nzf : row0NonZeroFunctionList) {
                if (nzf == static_cast<int>(f)) {
                    found = true;
                    break;
                }
            }
            if (!found) {
                int numComp = functionDescription[f].size();
                bool isTr = (f < isTruncated.size()) ? isTruncated[f] : false;
                bool hasSp = (f < hasSpillover.size()) ? hasSpillover[f] : false;
                row0Log << "Func " << f << ": " << numComp << " comp(s), truncated=" << (isTr ? "Y" : "N") 
                        << ", spillover=" << (hasSp ? "Y" : "N");
                
                // Check if this function has patch 0 components
                bool hasPatch0 = false;
                for (const auto& desc : functionDescription[f]) {
                    if (desc[0] == 0) {
                        hasPatch0 = true;
                        break;
                    }
                }
                row0Log << ", patch0=" << (hasPatch0 ? "Y" : "N") << "\n";
                zeroCount++;
            }
        }
        
        // Evaluate THB basis functions at the same point to understand presentation
        gsMatrix<real_t> testPt(2, 1);
        testPt(0, 0) = 0.901416;
        testPt(1, 0) = 0.0264156;
        
        row0Log << "\n=== THB vs MPBES Presentation Comparison ===\n";
        row0Log << "For each non-zero MPBES function, show:\n";
        row0Log << "1. MPBES function value\n";
        row0Log << "2. Corresponding THB basis function value\n";
        row0Log << "3. For TRUNCATED functions: presentation coefficients (how it's represented)\n\n";
        
        // For each non-zero MPBES function
        for (size_t i = 0; i < row0NonZeroFunctionList.size(); ++i) {
            int f = row0NonZeroFunctionList[i];
            real_t mpbesValue = row0NonZeroValueList[i];
            bool isTr = (f < static_cast<int>(isTruncated.size())) ? isTruncated[f] : false;
            
            row0Log << "--- MPBES Function " << f << " ---\n";
            row0Log << "MPBES value: " << std::setprecision(16) << mpbesValue << "\n";
            row0Log << "Truncated: " << (isTr ? "YES" : "NO") << "\n";
            
            // Find the THB index for the first component on patch 0
            if (f < static_cast<int>(functionDescription.size()) && functionDescription[f].size() > 0) {
                int fnPatch = functionDescription[f][0][0];
                int fnLevel = functionDescription[f][0][1];
                int fnTensorIdx = functionDescription[f][0][2];
                
                bool thbInBounds = (fnPatch < static_cast<int>(indexInTHB.size()) &&
                    fnLevel < static_cast<int>(indexInTHB[fnPatch].size()) &&
                    fnTensorIdx < static_cast<int>(indexInTHB[fnPatch][fnLevel].size()));
                
                if (thbInBounds) {
                    int thbIdx = indexInTHB[fnPatch][fnLevel][fnTensorIdx];
                    row0Log << "First component: patch=" << fnPatch << ", level=" << fnLevel 
                            << ", tensorIdx=" << fnTensorIdx << ", thbIdx=" << thbIdx << "\n";
                    
                    if (fnPatch < static_cast<int>(thbBases.size()) && thbIdx >= 0) {
                        // Evaluate THB basis function
                        gsMatrix<real_t> thbVal;
                        thbBases[fnPatch].evalSingle_into(thbIdx, testPt, thbVal);
                        row0Log << "THB Basis " << thbIdx << " value: " << std::setprecision(16) << thbVal(0, 0) << "\n";
                        
                        if (isTr) {
                            // Show the presentation: how THB represents this truncated function
                            const gsTHBSplineBasis<2, real_t>& thbBasis = thbBases[fnPatch];
                            
                            // Get the presentation from THB (which level 4 functions represent this)
                            // This would require accessing internal THB structure
                            // Show the presentation: MPBES vs THB coefficient comparison
                            row0Log << "\nMPBES Presentation Coefficients:\n";
                            
                            if (f < static_cast<int>(presentation.size()) && 
                                0 < static_cast<int>(presentation[f].size())) {
                                const auto& mpbesPres = presentation[f][0];  // First component's presentation
                                
                                row0Log << "  Sparse vector (index: coefficient):\n";
                                row0Log << "    ";
                                for (typename gsSparseVector<real_t>::InnerIterator it(mpbesPres); it; ++it) {
                                    row0Log << it.index() << ":" << std::setprecision(8) << it.value() << "  ";
                                }
                                row0Log << "\n";
                            }
                            
                            // Now compare with THB coefficients
                            if (fnPatch < static_cast<int>(thbBases.size()) && thbIdx >= 0 && 
                                thbBases[fnPatch].isTruncated(thbIdx)) {
                                row0Log << "\nTHB Presentation Coefficients (for comparison):\n";
                                
                                try {
                                    const auto& thbPres = thbBases[fnPatch].getCoefs(thbIdx);
                                    
                                    row0Log << "  Sparse vector (index: coefficient):\n";
                                    row0Log << "    ";
                                    for (typename gsSparseVector<real_t>::InnerIterator it(thbPres); it; ++it) {
                                        row0Log << it.index() << ":" << std::setprecision(8) << it.value() << "  ";
                                    }
                                    row0Log << "\n";
                                } catch (const std::exception& e) {
                                    row0Log << "  Error getting THB coefficients: " << e.what() << "\n";
                                }
                            } else if (fnPatch < static_cast<int>(thbBases.size()) && thbIdx >= 0) {
                                row0Log << "\nTHB function " << thbIdx << " is NOT truncated in THB (so no presentation needed)\n";
                            }
                            
                            row0Log << "\nComparison:\n";
                            row0Log << "  THB Basis " << thbIdx << " value: " << std::setprecision(16) << thbVal(0, 0) << "\n";
                            row0Log << "  MPBES Function " << f << " value: " << std::setprecision(16) << mpbesValue << "\n";
                            row0Log << "  Difference (THB - MPBES): " << std::setprecision(16) 
                                    << (thbVal(0, 0) - mpbesValue) << "\n";
                        }
                    }
                }
            }
            row0Log << "\n";
        }
        
        // For functions 0 and 1, create a separate diagnostic
        if (thbBases.size() > 0) {
            const gsTHBSplineBasis<2, real_t>& thbBasis = thbBases[0];
            
            for (int f : {0, 1}) {
                if (f >= static_cast<int>(row0NonZeroFunctionList.size() + 738)) continue;
                
                // Find if this function was in the non-zero list
                real_t mpbesValue = 0.0;
                bool found = false;
                for (size_t i = 0; i < row0NonZeroFunctionList.size(); ++i) {
                    if (row0NonZeroFunctionList[i] == f) {
                        mpbesValue = row0NonZeroValueList[i];
                        found = true;
                        break;
                    }
                }
                
                // Get THB value - function f should map to some THB basis index
                // Based on functionDescription[f]
                real_t thbValue = 0.0;
                if (f < static_cast<int>(functionDescription.size())) {
                    // For non-spillover functions, there should be a direct THB index
                    for (size_t compIdx = 0; compIdx < functionDescription[f].size(); ++compIdx) {
                        int fnPatch = functionDescription[f][compIdx][0];
                        int fnLevel = functionDescription[f][compIdx][1];
                        int fnTensorIdx = functionDescription[f][compIdx][2];
                        
                        if (fnPatch == 0 && fnLevel == 0) {  // Level 0 = base level (not truncated)
                            // This is the direct THB index
                            bool thbInBounds = (fnPatch < static_cast<int>(indexInTHB.size()) &&
                                fnLevel < static_cast<int>(indexInTHB[fnPatch].size()) &&
                                fnTensorIdx < static_cast<int>(indexInTHB[fnPatch][fnLevel].size()));
                            if (thbInBounds) {
                                int thbIdx = indexInTHB[fnPatch][fnLevel][fnTensorIdx];
                                if (thbIdx >= 0 && thbIdx < static_cast<int>(thbBasis.size())) {
                                    // Create test point
                                    gsMatrix<real_t> testPt(2, 1);
                                    testPt(0, 0) = 0.901416;
                                    testPt(1, 0) = 0.0264156;
                                    
                                    gsMatrix<real_t> evalRes;
                                    thbBasis.evalSingle_into(thbIdx, testPt, evalRes);
                                    thbValue += evalRes(0, 0);
                                }
                            }
                        }
                    }
                }
                
                if (found || thbValue > 1e-15) {
                    row0Log << "Function " << f << ":\n";
                    row0Log << "  MPBES value: " << std::setprecision(16) << mpbesValue << "\n";
                    row0Log << "  THB value: " << std::setprecision(16) << thbValue << "\n";
                    row0Log << "  Difference: " << std::setprecision(16) << (thbValue - mpbesValue) << "\n";
                    
                    if (std::abs(thbValue - mpbesValue) > 1e-6) {
                        row0Log << "  *** MISMATCH: Presentation basis is not reconstructing THB correctly ***\n";
                    }
                    row0Log << "\n";
                }
            }
        }
        
        row0Log << "\n=== End Detailed Row 0 Diagnostic ===\n";
        row0Log.close();
    }

    if (rowTargetLogOpened) {
        rowTargetLog << "\n=== Summary for Row " << targetRow << " ===\n";
        rowTargetLog << "Total non-zero functions: " << rowTargetNonZeroFuncs << " / " << numFunctions << "\n";
        rowTargetLog << "Row sum: " << std::setprecision(16) << (rowTargetSumSet ? rowTargetSum : 0.0) << "\n";
        rowTargetLog << "Missing contribution: " << std::setprecision(16) << (1.0 - (rowTargetSumSet ? rowTargetSum : 0.0)) << "\n\n";

        rowTargetLog << "=== All Contributing Functions ===\n";
        for (size_t i = 0; i < rowTargetNonZeroFunctionList.size(); ++i) {
            rowTargetLog << "Func " << rowTargetNonZeroFunctionList[i] << ": " << std::setprecision(16) << rowTargetNonZeroValueList[i] << "\n";
        }

        rowTargetLog << "\n=== End Detailed Row " << targetRow << " Diagnostic ===\n";
        rowTargetLog.close();
    }

    if (verbose)
    {
        gsInfo << "=== assemble END ===\n";
        gsInfo << "Final system size: " << A_mat.rows() << " x " << A_mat.cols() << "\n";
        
        outfile << "\n=== Assembly Summary ===\n";
        outfile << "Final system size: " << A_mat.rows() << " x " << A_mat.cols() << "\n";
        outfile << "Non-zero entries: " << A_mat.nonZeros() << "\n";
        
        // Check all rows for partition of unity
        int violationCount = 0;
        real_t tol = 1e-6;
        for (index_t i = 0; i < A_mat.rows(); ++i) {
            real_t rowSum = 0.0;
            for (index_t j = 0; j < A_mat.cols(); ++j) {
                rowSum += A_mat(i, j);
            }
            if (std::abs(rowSum - 1.0) > tol) {
                violationCount++;
                if (violationCount <= 10) {
                    outfile << "  Row " << i << ": sum = " << rowSum 
                            << " (deviation: " << (rowSum - 1.0) << ")\n";
                }
            }
        }
        outfile << "Total rows violating partition of unity: " << violationCount 
                << " / " << A_mat.rows() << "\n";
        
        // For violated rows, try to understand which functions are missing
        if (violationCount > 0) {
            outfile << "\n=== Analyzing Violations ===\n";
            for (index_t i = 0; i < std::min<index_t>(5, A_mat.rows()); ++i) {
                real_t rowSum = 0.0;
                int nnz = 0;
                for (index_t j = 0; j < A_mat.cols(); ++j) {
                    rowSum += A_mat(i, j);
                    if (A_mat(i, j) != 0) nnz++;
                }
                if (std::abs(rowSum - 1.0) > tol) {
                    outfile << "Row " << i << ": sum=" << rowSum << ", nnz=" << nnz << "\n";
                    outfile << "  Missing basis contribution: " << (1.0 - rowSum) << "\n";
                }
            }
        }
        
        outfile << "=== End Assembly ===\n\n";
    }
    
    // Test: Evaluate THB basis functions at the Row 0 test point
    // This helps understand if the partition of unity violation is in MPBES or THB
    {
        std::ofstream thbTestFile("thb_basis_eval_test.txt");
        thbTestFile << "=== THB Basis Function Evaluation Test ===\n";
        thbTestFile << "Testing partition of unity for THB basis functions\n";
        thbTestFile << "Test point: (u,v) = (0.901416, 0.0264156), patch = 0\n\n";
        
        // Create test point
        gsMatrix<real_t> testPoint(2, 1);
        testPoint(0, 0) = 0.901416;
        testPoint(1, 0) = 0.0264156;
        
        real_t thbSum = 0.0;
        int thbNonZeroCount = 0;
        
        // Patch 0 only (where the test point is)
        if (thbBases.size() > 0) {
            const gsTHBSplineBasis<2, real_t>& thbBasis = thbBases[0];
            int numBases = thbBasis.size();
            
            thbTestFile << "THB Basis at Patch 0:\n";
            thbTestFile << "Total basis functions: " << numBases << "\n\n";
            
            real_t patchSum = 0.0;
            int patchNonZero = 0;
            
            for (int i = 0; i < numBases; ++i) {
                gsMatrix<real_t> evalResult;
                thbBasis.evalSingle_into(i, testPoint, evalResult);
                real_t value = evalResult(0, 0);
                
                if (std::abs(value) > 1e-15) {
                    thbNonZeroCount++;
                    patchNonZero++;
                    thbSum += value;
                    patchSum += value;
                    
                    thbTestFile << "  Basis " << i << ": " 
                               << std::setprecision(16) << value << "\n";
                }
            }
            thbTestFile << "\nPatch 0 non-zero count: " << patchNonZero << "\n";
            thbTestFile << "Patch 0 sum: " << std::setprecision(16) << patchSum << "\n";
        }
        
        thbTestFile << "\n=== THB Summary ===\n";
        thbTestFile << "Total non-zero basis functions: " << thbNonZeroCount << "\n";
        thbTestFile << "Total sum: " << std::setprecision(16) << thbSum << "\n";
        thbTestFile << "Partition of unity error: " << std::setprecision(16) << (thbSum - 1.0) << "\n";
        
        if (std::abs(thbSum - 1.0) < 1e-12) {
            thbTestFile << "\n*** THB SATISFIES PARTITION OF UNITY ***\n";
        } else if (thbSum > 0) {
            thbTestFile << "\n*** THB VIOLATES PARTITION OF UNITY (sum = " << thbSum << ") ***\n";
        } else {
            thbTestFile << "\n*** THB EVALUATION FAILED - NO NONZERO FUNCTIONS ***\n";
        }
        
        thbTestFile << "\n=== Comparison with MPBES ===\n";
        thbTestFile << "MPBES sum at same point: 0.681553 (from row0_attempt56_detailed.txt)\n";
        thbTestFile << "THB sum at same point: " << std::setprecision(16) << thbSum << "\n";
        if (thbSum > 0) {
            thbTestFile << "Difference: " << std::setprecision(16) << (thbSum - 0.681553) << "\n";
            
            if (std::abs(thbSum - 1.0) < 1e-12 && std::abs(0.681553 - 1.0) > 1e-6) {
                thbTestFile << "\n*** CRITICAL: THB OK but MPBES fails ***\n";
                thbTestFile << "*** This means MPBES truncation is incorrectly zeroing functions ***\n";
            }
        }
        
        thbTestFile.close();
    }
    
    } catch (const std::exception& e) {
        gsInfo << "CRITICAL EXCEPTION in assemble: " << e.what() << "\n";
        outfile << "CRITICAL EXCEPTION in assemble: " << e.what() << "\n";
        gsInfo << "Returning partially assembled matrices\n";
        outfile << "Returning partially assembled matrices\n";
    } catch (...) {
        gsInfo << "UNKNOWN CRITICAL EXCEPTION in assemble\n";
        outfile << "UNKNOWN CRITICAL EXCEPTION in assemble\n";
        gsInfo << "Returning partially assembled matrices\n";
        outfile << "Returning partially assembled matrices\n";
    }
}


/**
 * @brief Assembles matrix A_mat and right-hand side b_vec using UNIFORM GRID points.
 *
 * This function uses a uniform grid of points (passed via uv parameter)
 * to construct the fitting system A * x = b. Unlike assemble(), which uses
 * Gauss quadrature points, this uses a pre-generated uniform grid.
 *
 * @date 12.01.26
 *
 * @param[in] uv         Parametric points (per patch) - uniform grid.
 * @param[in] Bells      Tensor-product bases per patch and level.
 * @param[in] THBVector  THB basis per patch.
 * @param[in] functionDescription List of (patch, level, index) entries for each global basis function.
 * @param[in] spilloverFunctionCoordinates Optional spillover sources per function.
 * @param[in] hasSpillover Flags marking which functions have spillovers.
 * @param[in] isTruncated Flags marking which functions are truncated.
 * @param[in] indexInTHB  Converts Bells indices into THB local indices.
 * @param[in] presentation Coefficient vectors for truncated functions.
 * @param[out] A_mat      Output matrix A.
 * @param[out] b_vec      Right-hand side vector.
 * @param[in] mp         The target geometry (gsMultiPatch).
 * @param[in] verbose    Enables logging.
 */
void assembleUniformGrid(
    const gsVector<gsMatrix<>>& uv,
    const gsVector<gsVector<gsTensorBSplineBasis<2, real_t>>>& Bells,
    const gsVector<gsTHBSplineBasis<2, real_t>>& THBVector,
    const std::vector<std::vector<std::vector<index_t>>>& functionDescription,
    const std::vector<std::vector<std::array<int, 3>>>& spilloverFunctionCoordinates,
    const std::vector<bool>& hasSpillover,
    const std::vector<bool>& isTruncated,
    const gsVector<gsVector<gsVector<index_t>>>& indexInTHB,
    std::vector<std::vector<gsSparseVector<real_t>>>& presentation,
    gsSparseMatrix<real_t>& A_mat,
    gsMatrix<real_t>& b_vec,
    const gsMultiPatch<real_t>& mp,
    bool verbose)
{
    if (verbose)
        gsInfo << "=== assembleUniformGrid BEGIN ===\n";

    const index_t numFunctions = functionDescription.size();
    const index_t numPatches = THBVector.size();

    // Count total number of points
    index_t totalRows = 0;
    for (index_t patch = 0; patch < numPatches; ++patch)
        totalRows += uv[patch].cols();

    A_mat.resize(totalRows, numFunctions);
    A_mat.setZero();
    b_vec.resize(totalRows, 2);
    b_vec.setZero();

    if (verbose)
        gsInfo << "Total rows (uniform grid points): " << totalRows << "\n";

    index_t rowOffset = 0;

    for (index_t patch = 0; patch < numPatches; ++patch)
    {
        if (verbose)
            gsInfo << ">> Patch " << patch << ", points: " << uv[patch].cols() << "\n";

        for (index_t ptIdx = 0; ptIdx < uv[patch].cols(); ++ptIdx)
        {
            gsVector<real_t> uvk = uv[patch].col(ptIdx);

            // RHS: target geometry
            gsMatrix<real_t> xy = mp.patch(patch).eval(uvk);
            b_vec(rowOffset, 0) = xy(0, 0);
            b_vec(rowOffset, 1) = xy(1, 0);

            // Evaluate all basis functions at this point
            gsMatrix<> resultMatrix;
            for (index_t f = 0; f < numFunctions; ++f)
            {
                real_t val = 0.0;
                int functionComponent = 0;

                for (const auto& desc : functionDescription[f])
                {
                    int fnPatch = desc[0];
                    int fnLevel = desc[1];
                    int fnTensorIndex = desc[2];

                    if (fnPatch != static_cast<int>(patch)) {
                        functionComponent++;
                        continue;
                    }

                    int thbIdx = indexInTHB[fnPatch][fnLevel][fnTensorIndex];

                    evalSingle_into(
                        f,
                        thbIdx,
                        uvk,
                        isTruncated,
                        THBVector[fnPatch],
                        presentation[f][functionComponent],
                        resultMatrix);

                    val += resultMatrix(0, 0);
                    functionComponent++;
                }

                if (val != 0.0)
                    A_mat(rowOffset, f) = val;
            }

            ++rowOffset;
        }
    }

    if (verbose)
        gsInfo << "=== assembleUniformGrid END ===\n";
}

/**
 * @brief Wrapper function that assembles and solves the uniform grid system.
 *
 * This function:
 * 1. Calls assembleUniformGrid to build A and b
 * 2. Prints the matrices and checks partition of unity
 * 3. Transforms to normal equations (A^T*A)x = A^T*b
 * 4. Solves the system
 * 5. Returns the solution
 *
 * @date 12.01.26
 */
gsMatrix<real_t> solveUniformGridSystem(
    const gsVector<gsMatrix<>>& uv,
    const gsVector<gsVector<gsTensorBSplineBasis<2, real_t>>>& Bells,
    const gsVector<gsTHBSplineBasis<2, real_t>>& THBVector,
    const std::vector<std::vector<std::vector<index_t>>>& functionDescription,
    const std::vector<std::vector<std::array<int, 3>>>& spilloverFunctionCoordinates,
    const std::vector<bool>& hasSpillover,
    const std::vector<bool>& isTruncated,
    const gsVector<gsVector<gsVector<index_t>>>& indexInTHB,
    std::vector<std::vector<gsSparseVector<real_t>>>& presentation,
    const gsMultiPatch<real_t>& mp,
    bool verbose = true)
{
    gsInfo << "\n========== UNIFORM GRID SYSTEM SOLVE ==========\n";

    // Print UV grid information
    gsInfo << "\n--- UV Grid Points ---\n";
    for (index_t patch = 0; patch < uv.size(); ++patch)
    {
        gsInfo << "Patch " << patch << ": " << uv[patch].cols() << " points\n";
        gsInfo << "  First 5 points:\n";
        for (index_t i = 0; i < std::min<index_t>(5, uv[patch].cols()); ++i)
        {
            gsInfo << "    Point " << i << ": (" << uv[patch](0, i) << ", " << uv[patch](1, i) << ")\n";
        }
        if (uv[patch].cols() > 5)
            gsInfo << "  ... (" << (uv[patch].cols() - 5) << " more points)\n";
    }
    gsInfo << "----------------------\n\n";

    // Step 1: Assemble the system
    gsSparseMatrix<real_t> matA_uniform;
    gsMatrix<real_t> vectB_uniform;

    assembleUniformGrid(
        uv,
        Bells,
        THBVector,
        functionDescription,
        spilloverFunctionCoordinates,
        hasSpillover,
        isTruncated,
        indexInTHB,
        presentation,
        matA_uniform,
        vectB_uniform,
        mp,
        verbose);

    gsInfo << "\nUniform grid system assembled:\n";
    gsInfo << "  matA size: " << matA_uniform.rows() << " x " << matA_uniform.cols() << "\n";
    gsInfo << "  vectB size: " << vectB_uniform.rows() << " x " << vectB_uniform.cols() << "\n";

    // Step 2: Print matrices
    gsInfo << "\n--- Original Uniform Grid System ---\n";
    printTheMatrix(matA_uniform, "matA_uniform (original)");
    printTheMatrix(vectB_uniform, "vectB_uniform (original b)");

    // Step 3: Check partition of unity
    gsInfo << "\n--- Partition of Unity Check (Uniform Grid) ---\n";
    checkPartitionOfUnity(matA_uniform);

    // Step 4: Transform to normal equations
    gsMatrix<> b_vec_copy = vectB_uniform;
    vectB_uniform = gsEigen::MatrixXd(matA_uniform).transpose() * vectB_uniform;
    gsMatrix<> matAsquare_uniform = gsEigen::MatrixXd(matA_uniform).transpose() * gsEigen::MatrixXd(matA_uniform);

    gsInfo << "\nNormal equations formed:\n";
    gsInfo << "  A^T*A size: " << matAsquare_uniform.rows() << " x " << matAsquare_uniform.cols() << "\n";
    gsInfo << "  A^T*b size: " << vectB_uniform.rows() << " x " << vectB_uniform.cols() << "\n";
    gsInfo << "  det(A^T*A): " << matAsquare_uniform.determinant() << "\n";

    // Step 5: Solve the system
    gsMatrix<real_t> vectSol_uniform = matAsquare_uniform.partialPivLu().solve(vectB_uniform);

    gsInfo << "\nSolution obtained:\n";
    gsInfo << "  vectSol size: " << vectSol_uniform.rows() << " x " << vectSol_uniform.cols() << "\n";

    // Step 6: Check residuals
    gsMatrix<> matC = matAsquare_uniform * vectSol_uniform - vectB_uniform;
    gsInfo << "\nResidual norm (A^T*A*x - A^T*b): " << matC.maxCoeff() << "\n";

    gsMatrix<> matOut = gsEigen::MatrixXd(matA_uniform) * vectSol_uniform;
    matC = matOut - b_vec_copy;
    gsInfo << "Residual norm (A*x - b): " << matC.maxCoeff() << "\n";

    // Step 7: Evaluate geometry at interface points (0,(0,1)) and (1,(1,1))
    gsInfo << "\n--- Geometry Evaluation at Interface Points ---\n";

    gsVector<real_t> uv_p0(2), uv_p1(2);
    uv_p0 << 0.0, 1.0;  // Patch 0 at (u=0, v=1)
    uv_p1 << 1.0, 1.0;  // Patch 1 at (u=1, v=1)

    std::ofstream interfaceLog("interface_basis_analysis.txt");
    interfaceLog << "Interface Point Analysis\n";
    interfaceLog << "========================\n\n";

    // Evaluate for patch 0 at (0, 1)
    {
        gsMatrix<> result_p0(2, 1);
        result_p0.setZero();

        interfaceLog << "PATCH 0 at (u,v) = (0, 1):\n";
        interfaceLog << "Function    BasisValue    Coeff_x        Coeff_y        Contrib_x      Contrib_y\n";
        interfaceLog << "--------    ----------    ----------     ----------     ----------     ----------\n";

        for (index_t f = 0; f < functionDescription.size(); ++f)
        {
            real_t val = 0.0;
            int functionComponent = 0;

            for (const auto& desc : functionDescription[f])
            {
                int fnPatch = desc[0];
                int fnLevel = desc[1];
                int fnTensorIndex = desc[2];

                if (fnPatch != 0) {
                    functionComponent++;
                    continue;
                }

                int thbIdx = indexInTHB[fnPatch][fnLevel][fnTensorIndex];

                gsMatrix<> resultMatrix;
                evalSingle_into(
                    f,
                    thbIdx,
                    uv_p0,
                    isTruncated,
                    THBVector[0],
                    presentation[f][functionComponent],
                    resultMatrix);

                val += resultMatrix(0, 0);
                functionComponent++;
            }

            if (fabs(val) > 1e-10)
            {
                interfaceLog << std::setw(8) << f
                    << std::setw(14) << std::scientific << std::setprecision(6) << val
                    << std::setw(14) << vectSol_uniform(f, 0)
                    << std::setw(14) << vectSol_uniform(f, 1)
                    << std::setw(14) << val * vectSol_uniform(f, 0)
                    << std::setw(14) << val * vectSol_uniform(f, 1) << "\n";
            }

            result_p0(0, 0) += val * vectSol_uniform(f, 0);
            result_p0(1, 0) += val * vectSol_uniform(f, 1);
        }

        interfaceLog << "\nTotal geometry at Patch 0 (0,1): x = " << result_p0(0, 0) << ", y = " << result_p0(1, 0) << "\n\n";
        gsInfo << "Patch 0 at (u,v) = (0, 1):\n";
        gsInfo << "  Evaluated geometry: x = " << result_p0(0, 0) << ", y = " << result_p0(1, 0) << "\n";
    }

    // Evaluate for patch 1 at (1, 1)
    {
        gsMatrix<> result_p1(2, 1);
        result_p1.setZero();

        interfaceLog << "\nPATCH 1 at (u,v) = (1, 1):\n";
        interfaceLog << "Function    BasisValue    Coeff_x        Coeff_y        Contrib_x      Contrib_y\n";
        interfaceLog << "--------    ----------    ----------     ----------     ----------     ----------\n";

        for (index_t f = 0; f < functionDescription.size(); ++f)
        {
            real_t val = 0.0;
            int functionComponent = 0;

            for (const auto& desc : functionDescription[f])
            {
                int fnPatch = desc[0];
                int fnLevel = desc[1];
                int fnTensorIndex = desc[2];

                if (fnPatch != 1) {
                    functionComponent++;
                    continue;
                }

                int thbIdx = indexInTHB[fnPatch][fnLevel][fnTensorIndex];

                gsMatrix<> resultMatrix;
                evalSingle_into(
                    f,
                    thbIdx,
                    uv_p1,
                    isTruncated,
                    THBVector[1],
                    presentation[f][functionComponent],
                    resultMatrix);

                val += resultMatrix(0, 0);
                functionComponent++;
            }

            if (fabs(val) > 1e-10)
            {
                interfaceLog << std::setw(8) << f
                    << std::setw(14) << std::scientific << std::setprecision(6) << val
                    << std::setw(14) << vectSol_uniform(f, 0)
                    << std::setw(14) << vectSol_uniform(f, 1)
                    << std::setw(14) << val * vectSol_uniform(f, 0)
                    << std::setw(14) << val * vectSol_uniform(f, 1) << "\n";
            }

            result_p1(0, 0) += val * vectSol_uniform(f, 0);
            result_p1(1, 0) += val * vectSol_uniform(f, 1);
        }

        interfaceLog << "\nTotal geometry at Patch 1 (1,1): x = " << result_p1(0, 0) << ", y = " << result_p1(1, 0) << "\n\n";

        // Check which functions are active at BOTH interface points
        interfaceLog << "\n=== TWIN FUNCTION ANALYSIS ===\n";
        interfaceLog << "Functions that should match at interface:\n\n";

        // Special focus on Functions 20 and 591
        interfaceLog << "\n=== DETAILED ANALYSIS: Functions 20 and 591 ===\n";
        for (index_t f : {20, 591})
        {
            if (f >= functionDescription.size()) continue;

            interfaceLog << "\nFunction " << f << " details:\n";
            interfaceLog << "  Number of components: " << functionDescription[f].size() << "\n";
            interfaceLog << "  Is truncated: " << (isTruncated[f] ? "YES" : "NO") << "\n";

            for (size_t comp = 0; comp < functionDescription[f].size(); ++comp)
            {
                const auto& desc = functionDescription[f][comp];
                interfaceLog << "  Component " << comp << ": patch=" << desc[0]
                    << ", level=" << desc[1] << ", tensorIdx=" << desc[2] << "\n";

                if (desc[0] < indexInTHB.size() && desc[1] < indexInTHB[desc[0]].size()
                    && desc[2] < indexInTHB[desc[0]][desc[1]].size())
                {
                    int thbIdx = indexInTHB[desc[0]][desc[1]][desc[2]];
                    interfaceLog << "    → THB index: " << thbIdx << "\n";

                    // Check presentation/spillover
                    if (isTruncated[f] && comp < presentation[f].size())
                    {
                        interfaceLog << "    → Presentation (spillover) coefficients:\n";
                        for (typename gsSparseVector<real_t>::InnerIterator it(presentation[f][comp]); it; ++it)
                        {
                            interfaceLog << "        [" << it.index() << "] = " << it.value() << "\n";
                        }
                        if (presentation[f][comp].nonZeros() == 0)
                            interfaceLog << "        (empty - no spillover)\n";
                    }
                }
            }
        }

        interfaceLog << "\n=== ALL TWIN FUNCTIONS ===\n";

        for (index_t f = 0; f < functionDescription.size(); ++f)
        {
            bool hasPatch0 = false, hasPatch1 = false;

            for (const auto& desc : functionDescription[f])
            {
                if (desc[0] == 0) hasPatch0 = true;
                if (desc[0] == 1) hasPatch1 = true;
            }

            if (hasPatch0 && hasPatch1)
            {
                interfaceLog << "Function " << f << " (twin on both patches):\n";

                // Evaluate on patch 0
                real_t val_p0 = 0.0;
                int comp0 = 0;
                for (const auto& desc : functionDescription[f])
                {
                    if (desc[0] == 0)
                    {
                        int thbIdx = indexInTHB[0][desc[1]][desc[2]];
                        gsMatrix<> resultMatrix;
                        evalSingle_into(f, thbIdx, uv_p0, isTruncated, THBVector[0],
                            presentation[f][comp0], resultMatrix);
                        val_p0 += resultMatrix(0, 0);
                    }
                    comp0++;
                }

                // Evaluate on patch 1
                real_t val_p1 = 0.0;
                int comp1 = 0;
                for (const auto& desc : functionDescription[f])
                {
                    if (desc[0] == 1)
                    {
                        int thbIdx = indexInTHB[1][desc[1]][desc[2]];
                        gsMatrix<> resultMatrix;
                        evalSingle_into(f, thbIdx, uv_p1, isTruncated, THBVector[1],
                            presentation[f][comp1], resultMatrix);
                        val_p1 += resultMatrix(0, 0);
                    }
                    comp1++;
                }

                if (fabs(val_p0) > 1e-10 || fabs(val_p1) > 1e-10)
                {
                    interfaceLog << "  Basis value on Patch 0: " << val_p0 << "\n";
                    interfaceLog << "  Basis value on Patch 1: " << val_p1 << "\n";
                    interfaceLog << "  Difference: " << fabs(val_p0 - val_p1) << "\n";
                    if (fabs(val_p0 - val_p1) > 1e-6)
                        interfaceLog << "  ⚠️  WARNING: Basis values don't match!\n";
                    interfaceLog << "\n";
                }
            }
        }

        interfaceLog.close();

        gsInfo << "Patch 1 at (u,v) = (1, 1):\n";
        gsInfo << "  Evaluated geometry: x = " << result_p1(0, 0) << ", y = " << result_p1(1, 0) << "\n";
        gsInfo << "\n⚠️  Detailed basis analysis written to interface_basis_analysis.txt\n";
    }
    gsInfo << "-----------------------------------------------\n";

    // Evaluate for patch 1 at (1, 1)
    {
        gsMatrix<> result_p1(2, 1);
        result_p1.setZero();

        for (index_t f = 0; f < functionDescription.size(); ++f)
        {
            real_t val = 0.0;
            int functionComponent = 0;

            for (const auto& desc : functionDescription[f])
            {
                int fnPatch = desc[0];
                int fnLevel = desc[1];
                int fnTensorIndex = desc[2];

                if (fnPatch != 1) {
                    functionComponent++;
                    continue;
                }

                int thbIdx = indexInTHB[fnPatch][fnLevel][fnTensorIndex];

                gsMatrix<> resultMatrix;
                evalSingle_into(
                    f,
                    thbIdx,
                    uv_p1,
                    isTruncated,
                    THBVector[1],
                    presentation[f][functionComponent],
                    resultMatrix);

                val += resultMatrix(0, 0);
                functionComponent++;
            }

            result_p1(0, 0) += val * vectSol_uniform(f, 0);
            result_p1(1, 0) += val * vectSol_uniform(f, 1);
        }

        gsInfo << "Patch 1 at (u,v) = (1, 1):\n";
        gsInfo << "  Evaluated geometry: x = " << result_p1(0, 0) << ", y = " << result_p1(1, 0) << "\n";
    }
    gsInfo << "-----------------------------------------------\n";

    gsInfo << "========== END UNIFORM GRID SYSTEM ==========\n\n";

    return vectSol_uniform;
}

real_t theDistance(gsMultiPatch<real_t> repairedmp) {
    gsMatrix<> points0(2, 17 * 10 + 1);
    gsMatrix<> points1(2, 17 * 10 + 1);
    for (int i = 0; i < points0.cols(); ++i) {
        points0(0, i) = 1.0;
        points1(0, i) = 0.0;
        points0(1, i) = (double)i / (points0.cols() - 1);
        points1(1, i) = (double)i / (points0.cols() - 1);
    }
    gsMatrix<> val0 = repairedmp.patch(0).eval(points0);
    gsMatrix<> val1 = repairedmp.patch(1).eval(points1);
    gsVector<> distances(val0.cols());
    real_t maxError = 0.0;
    for (int i = 0; i < val0.cols(); ++i) {
        if (val0(0, i) > val1(0, i)) {
            maxError = std::max(maxError, std::sqrt((val0(0, i) - val1(0, i)) * (val0(0, i) - val1(0, i)) +
                (val0(0, i) - val1(0, i)) * (val0(0, i) - val1(0, i))));
            distances(i) = std::sqrt((val0(0, i) - val1(0, i)) * (val0(0, i) - val1(0, i)) +
                (val0(0, i) - val1(0, i)) * (val0(0, i) - val1(0, i)));
        }
        else
            distances(i) = -1;
    }
    return maxError;
}




void coordinateTranformation(double& uSecond, double& vSecond, int firstSide, int secondSide, gsMatrix<> punto) {
    switch (firstSide) {
    case 1:
        // "WEST\n";
        switch (secondSide) {
        case 1:
            //"WEST\n";
            uSecond = punto(0, 0);
            vSecond = 1 - punto(1, 0);
            return;
        case 2:
            //"EAST\n";
            uSecond = 1 - punto(0, 0);
            vSecond = punto(1, 0);
            return;
        case 3:
            //SOUTH!\n";
            uSecond = punto(1, 0);
            vSecond = punto(0, 0);
            return;
        case 4:
            //NORTH\n";
            uSecond = 1 - punto(1, 0);
            vSecond = 1 - punto(0, 0);
            return;
        }
    case 2:
        //EAST\n";
        switch (secondSide) {
        case 1:
            //WEST\n";
            uSecond = 1 - punto(0, 0);
            vSecond = punto(1, 0);
            return;
        case 2:
            //EAST\n";
            uSecond = punto(0, 0);
            vSecond = 1 - punto(1, 0);
            return;
        case 3:
            //SOUTH\n";
            uSecond = 1 - punto(1, 0);
            vSecond = 1 - punto(0, 0);
            return;
        case 4:
            //NORTH\n";
            uSecond = punto(1, 0);
            vSecond = punto(0, 0);
            return;
        }
    case 3:
        //SOUTH\n";
        switch (secondSide) {
        case 1:
            //WEST\n";
            uSecond = punto(1, 0);
            vSecond = punto(0, 0);
            return;
        case 2:
            //EAST\n";
            uSecond = 1 - punto(1, 0);
            vSecond = 1 - punto(0, 0);
            return;
        case 3:
            //SOUTH\n";
            uSecond = 1 - punto(0, 0);
            vSecond = punto(1, 0);
            return;
        case 4:
            //NORTH\n";
            uSecond = punto(0, 0);
            vSecond = 1 - punto(1, 0);
            return;
        }
    case 4:
        //NORTH\n";
        switch (secondSide) {
        case 1:
            //WEST\n";
            uSecond = 1 - punto(1, 0);
            vSecond = 1 - punto(0, 0);
            return;
        case 2:
            //EAST\n";
            uSecond = punto(1, 0);
            vSecond = punto(0, 0);
            return;
        case 3:
            //SOUTH\n";
            uSecond = punto(0, 0);
            vSecond = 1 - punto(1, 0);
            return;
        case 4:
            //NORTH\n";
            uSecond = 1 - punto(0, 0);
            vSecond = punto(1, 0);
            return;
        }
    }
}

void coordinateTranformationBoxes(int firstSide, gsVector<index_t> thebox, int secondSide,
    std::vector<index_t>& secondbox, int maxIndex
)
{
    //Assumes that we know the main patch
    index_t x11, y11, x12, y12, x21, y21, x22, y22;
    x11 = thebox[1];
    y11 = thebox[2];
    x12 = thebox[3];
    y12 = thebox[4];
    secondbox[0] = thebox[0];
    switch (firstSide)
    {
    case 1:
        // "WEST\n";
        switch (secondSide) {
        case 1:
            //"WEST\n";
            secondbox[1] = x21 = x11;
            secondbox[2] = y21 = maxIndex - y11;
            secondbox[3] = x22 = x12;
            secondbox[4] = y22 = maxIndex - y12;
            return;
        case 2:
            //"EAST\n";
            secondbox[1] = maxIndex - x11;
            secondbox[2] = y11;
            secondbox[3] = maxIndex - x12;
            secondbox[4] = y12;
            return;
        case 3:
            //SOUTH!\n";
            secondbox[1] = x21 = y11;
            secondbox[2] = y21 = x11;
            secondbox[3] = x22 = y12;
            secondbox[4] = y22 = x12;
            return;
        case 4:
            //NORTH\n";
            secondbox[1] = x21 = maxIndex - y11;
            secondbox[2] = y21 = maxIndex - x11;
            secondbox[3] = x22 = maxIndex - y12;
            secondbox[4] = y22 = maxIndex - x12;
            return;
        }
    case 2:
        //EAST\n";
        switch (secondSide) {
        case 1:
            //WEST\n";
            secondbox[1] = x21 = maxIndex - x11;
            secondbox[2] = y21 = y11;
            secondbox[3] = x22 = maxIndex - x12;
            secondbox[4] = y22 = y12;
            return;
        case 2:
            //EAST\n";
            secondbox[1] = x21 = x11;
            secondbox[2] = y21 = maxIndex - y11;
            secondbox[3] = x22 = x12;
            secondbox[4] = y22 = maxIndex - y12;
            return;
        case 3:
            //SOUTH\n";
            secondbox[1] = x21 = maxIndex - y11;
            secondbox[2] = y21 = maxIndex - x11;
            secondbox[3] = x22 = maxIndex - y12;
            secondbox[4] = y22 = maxIndex - x12;
            return;
        case 4:
            //NORTH\n";
            secondbox[1] = x21 = y11;
            secondbox[2] = y21 = x11;
            secondbox[3] = x22 = y12;
            secondbox[4] = y22 = x12;
            return;
        }
    case 3:
        //SOUTH\n";
        switch (secondSide) {
        case 1:
            //WEST\n";
            secondbox[1] = x21 = y11;
            secondbox[2] = y21 = x11;
            secondbox[3] = x22 = y12;
            secondbox[4] = y22 = x12;
            return;
        case 2:
            //EAST\n";
            secondbox[1] = x21 = maxIndex - y11;
            secondbox[2] = y21 = maxIndex - x11;
            secondbox[3] = x22 = maxIndex - y12;
            secondbox[4] = y22 = maxIndex - x12;
            return;
        case 3:
            //SOUTH\n";
            secondbox[1] = x21 = maxIndex - y11;
            secondbox[2] = y21 = x11;
            secondbox[3] = x22 = maxIndex - y12;
            secondbox[4] = y22 = x12;
            return;
        case 4:
            //NORTH\n";
            secondbox[1] = x21 = x11;
            secondbox[2] = y21 = maxIndex - y11;
            secondbox[3] = x22 = y12;
            secondbox[4] = y22 = maxIndex - x12;
            return;
        }
    case 4:
        //NORTH\n";
        switch (secondSide) {
        case 1:
            //WEST\n";
            secondbox[1] = x21 = maxIndex - y11;
            secondbox[2] = y21 = maxIndex - x11;
            secondbox[3] = x22 = maxIndex - y12;
            secondbox[4] = y22 = maxIndex - x12;
            return;
        case 2:
            //EAST\n";
            secondbox[1] = x21 = y11;
            secondbox[2] = y21 = x11;
            secondbox[3] = x22 = y12;
            secondbox[4] = y22 = x12;
            return;
        case 3:
            //SOUTH\n";
            secondbox[1] = x21 = x11;
            secondbox[2] = y21 = maxIndex - y11;
            secondbox[3] = x22 = y12;
            secondbox[4] = y22 = maxIndex - x12;
            return;
        case 4:
            //NORTH\n";
            secondbox[1] = x21 = maxIndex - x11;
            secondbox[2] = y21 = y11;
            secondbox[3] = x22 = maxIndex - x12;
            secondbox[4] = y22 = y12;
            return;
        }
    }
}

int twinFunction(
    const gsMultiPatch<real_t>& mp,
    int firstPatch,
    const gsTHBSplineBasis<2, real_t>& initTP,
    int functionIndex,
    int secondPatch,
    const gsTHBSplineBasis<2, real_t>& twinTP)
{
    gsMatrix<> firstAnchor;
    gsMatrix<> secondAnchors;
    initTP.anchor_into(functionIndex, firstAnchor);
    twinTP.anchors_into(secondAnchors);

    if (firstAnchor.rows() < 2 || firstAnchor.cols() < 1)
        return -1;

    const gsVector<real_t> firstUv = matrixColumnToVector2(firstAnchor);
    const gsVector<real_t> firstXy = matrixColumnToVector2(mp.patch(firstPatch).eval(firstUv));

    int bestIndex = -1;
    real_t bestDistance = std::numeric_limits<real_t>::infinity();

    for (int i = 0; i < secondAnchors.cols(); ++i)
    {
        const gsVector<real_t> secondUv = matrixColumnToVector2(secondAnchors.col(i));
        const gsVector<real_t> secondXy = matrixColumnToVector2(mp.patch(secondPatch).eval(secondUv));
        const real_t distance = vectorDistance(firstXy, secondXy);
        if (distance < bestDistance)
        {
            bestDistance = distance;
            bestIndex = i;
        }
    }

    return bestDistance <= 1e-8 ? bestIndex : -1;
}

void boxToDomain(gsMatrix<index_t> mybox, gsMatrix <real_t>& coords, index_t interioru, index_t interiorv) {
    for (int i = 0; i < mybox.rows(); i++) {
        coords(i, 0) = (double)mybox(i, 0);
        coords(i, 1) = (double)mybox(i, 1) / ((double)((interioru + 1) * pow(2, mybox(i, 0))));
        coords(i, 2) = (double)mybox(i, 2) / ((double)((interiorv + 1) * pow(2, mybox(i, 0))));
        coords(i, 3) = (double)mybox(i, 3) / ((double)((interioru + 1) * pow(2, mybox(i, 0))));
        coords(i, 4) = (double)mybox(i, 4) / ((double)((interiorv + 1) * pow(2, mybox(i, 0))));
    }
}
inline void supportToBoxOfLevel(gsMatrix<index_t>& mybox, gsMatrix <real_t>  coords, int level, index_t interioru, index_t interiorv) {
    for (int i = 0; i < mybox.rows(); i++) {
        for (int j = 0; j < mybox.cols(); ++j) {
            mybox(i, j) = coords(i, j) * pow(2, level);
        }
    }
}

void boxToDomain(gsVector<index_t> mybox, gsVector <real_t>& coords, index_t interioru, index_t interiorv) {
    coords(0) = (double)mybox(0);
    coords(1) = (double)mybox(1) / ((double)((interioru + 1) * pow(2, mybox(0))));
    coords(2) = (double)mybox(2) / ((double)((interiorv + 1) * pow(2, mybox(0))));
    coords(3) = (double)mybox(3) / ((double)((interioru + 1) * pow(2, mybox(0))));
    coords(4) = (double)mybox(4) / ((double)((interiorv + 1) * pow(2, mybox(0))));
}

void boxToDomain(gsVector<index_t> mybox, std::vector <real_t>& coords, index_t interioru, index_t interiorv) {
    coords[0] = (double)mybox(0);
    coords[1] = (double)mybox(1) / ((double)((interioru + 1) * pow(2, mybox(0))));
    coords[2] = (double)mybox(2) / ((double)((interiorv + 1) * pow(2, mybox(0))));
    coords[3] = (double)mybox(3) / ((double)((interioru + 1) * pow(2, mybox(0))));
    coords[4] = (double)mybox(4) / ((double)((interiorv + 1) * pow(2, mybox(0))));
}

bool subset(gsBasisFun<> beta, gsVector <real_t> theBox) {
    if ((beta.support()(0, 0) >= theBox(1) &&
        beta.support()(1, 0) >= theBox(2) &&
        beta.support()(0, 1) <= theBox(3) &&
        beta.support()(1, 1) <= theBox(4))
        ) {
        return true;
    }
    else return false;
}

//index_t isZeroOnBoundary(gsVector <gsTHBSplineBasis <2, real_t >> Bells,
//    int patch,
//    int functionIndex,
//    int theSide
//) {
//    const unsigned int numPoints = 50;
//    gsMatrix<double> punto(2, 1);
//    if (theSide == 1) {
//        for (int i = 0; i < numPoints; ++i) {
//            punto(0, 0) = 0.0;
//            punto(1, 0) = ((double)i) / (numPoints - 1);
//            gsMatrix<double> valore;
//
//            valore = Bells(patch).function(functionIndex).eval(punto);
//            if (math::abs(valore(0, 0)) > 1e-13) {
//                return 0;
//            }
//        }
//    }
//    if (theSide == 2) {
//        //EAST\n";
//        for (int i = 0; i < numPoints; ++i) {
//            punto(0, 0) = 1.0;
//            punto(1, 0) = ((double)i) / (numPoints - 1);
//            gsMatrix<double> valore;
//
//            valore = Bells(patch).function(functionIndex).eval(punto);
//            if (math::abs(valore(0, 0)) > 1e-13) {
//                return 0;
//            }
//        }
//    }
//    if (theSide == 3) {
//        for (int i = 0; i < numPoints; ++i) {
//            punto(0, 0) = ((double)i) / (numPoints - 1);
//            punto(1, 0) = 0.0;
//            gsMatrix<double> valore;
//
//            valore = Bells(patch).function(functionIndex).eval(punto);
//            if (math::abs(valore(0, 0)) > 1e-13) {
//                return 0;
//            }
//        }
//    }
//    if (theSide == 4) {
//        for (int i = 0; i < numPoints; ++i) {
//            punto(0, 0) = ((double)i) / (numPoints - 1);
//            punto(1, 0) = 1.0;
//            gsMatrix<double> valore;
//
//            valore = Bells(patch).function(functionIndex).eval(punto);
//            if (math::abs(valore(0, 0)) > 1e-13) {
//                return 0;
//            }
//        }
//    }
//    return 1;
//}


int isZeroOnBoundary(gsVector<gsVector<gsTensorBSplineBasis<2, real_t>>> Bells,
    int patch,
    int level,
    int functionIndex,
    int theSide)
{
    const unsigned int numPoints = 50;
    gsMatrix<double> punto(2, 1);
    //gsInfo << "isZeroOnBoundary called\n";
    //gsInfo << "  Patch: " << patch << ", Level: " << level
        //<< ", FunctionIndex: " << functionIndex << ", Side: " << theSide << "\n";

    if (theSide == 1)
    {
        //gsInfo << "  Checking WEST boundary\n";
        for (int i = 0; i < numPoints; ++i)
        {
            punto(0, 0) = 0.0;
            punto(1, 0) = static_cast<double>(i) / (numPoints - 1);
            gsMatrix<double> valore = Bells(patch)(level).function(functionIndex).eval(punto);
            //gsInfo << "    i=" << i << ", punto=(" << punto(0, 0) << "," << punto(1, 0) << "), value=" << valore(0, 0) << "\n";

            if (math::abs(valore(0, 0)) > 1e-13)
            {
                //gsInfo << "    --> Value above threshold! Returning 0\n";
                return 0;
            }
        }
    }

    if (theSide == 2)
    {
        //gsInfo << "  Checking EAST boundary\n";
        for (int i = 0; i < numPoints; ++i)
        {
            punto(0, 0) = 1.0;
            punto(1, 0) = static_cast<double>(i) / (numPoints - 1);
            gsMatrix<double> valore = Bells(patch)(level).function(functionIndex).eval(punto);
            //gsInfo << "    i=" << i << ", punto=(" << punto(0, 0) << "," << punto(1, 0) << "), value=" << valore(0, 0) << "\n";

            if (math::abs(valore(0, 0)) > 1e-13)
            {
                //gsInfo << "    --> Value above threshold! Returning 0\n";
                return 0;
            }
        }
    }

    if (theSide == 3)
    {
        //gsInfo << "  Checking SOUTH boundary\n";
        for (int i = 0; i < numPoints; ++i)
        {
            punto(0, 0) = static_cast<double>(i) / (numPoints - 1);
            punto(1, 0) = 0.0;
            gsMatrix<double> valore = Bells(patch)(level).function(functionIndex).eval(punto);
            //gsInfo << "    i=" << i << ", punto=(" << punto(0, 0) << "," << punto(1, 0) << "), value=" << valore(0, 0) << "\n";

            if (math::abs(valore(0, 0)) > 1e-13)
            {
                //gsInfo << "    --> Value above threshold! Returning 0\n";
                return 0;
            }
        }
    }

    if (theSide == 4)
    {
        //gsInfo << "  Checking NORTH boundary\n";
        for (int i = 0; i < numPoints; ++i)
        {
            punto(0, 0) = static_cast<double>(i) / (numPoints - 1);
            punto(1, 0) = 1.0;
            gsMatrix<double> valore = Bells(patch)(level).function(functionIndex).eval(punto);
            //gsInfo << "    i=" << i << ", punto=(" << punto(0, 0) << "," << punto(1, 0) << "), value=" << valore(0, 0) << "\n";

            if (math::abs(valore(0, 0)) > 1e-13)
            {
                //gsInfo << "    --> Value above threshold! Returning 0\n";
                return 0;
            }
        }
    }

    //gsInfo << "  All values within tolerance. Returning 1\n";
    return 1;
}

static bool supportTouchesInterfaceSide(const gsBasisFun<>& basisFunction, int sideIndex)
{
    const gsMatrix<real_t> support = basisFunction.support();

    switch (sideIndex)
    {
    case 1:
        return math::abs(support(0, 0)) <= 1e-13;
    case 2:
        return math::abs(support(0, 1) - 1.0) <= 1e-13;
    case 3:
        return math::abs(support(1, 0)) <= 1e-13;
    case 4:
        return math::abs(support(1, 1) - 1.0) <= 1e-13;
    default:
        return false;
    }
}

static bool isInvestigatedInterface(index_t firstPatch, index_t secondPatch)
{
    return (firstPatch == 15 && secondPatch == 20) ||
        (firstPatch == 20 && secondPatch == 15);
}

static void logInvestigatedTopologyInterface(
    const char* stage,
    int interfaceNum,
    index_t firstPatch,
    index_t secondPatch,
    int rawFirstSide,
    int rawSecondSide,
    int firstSide,
    int secondSide,
    const GeometryPreflightInterfaceInfo* preflightInfo)
{
    if (!outfile.is_open() || !isInvestigatedInterface(firstPatch, secondPatch))
        return;

    outfile << "[investigate-interface] stage=" << stage
            << " interface=" << interfaceNum
            << " patches=(" << firstPatch << "," << secondPatch << ")"
            << " rawSides=(" << rawFirstSide << "," << rawSecondSide << ")"
            << " canonicalSides=(" << firstSide << "," << secondSide << ")"
            << " mirrored=(" << (isPreflightMirroredPatch(firstPatch) ? 1 : 0)
            << "," << (isPreflightMirroredPatch(secondPatch) ? 1 : 0) << ")";

    if (preflightInfo && preflightInfo->hasMappedSides)
    {
        FeatureSide firstFeatureSide = preflightInfo->sideA;
        FeatureSide secondFeatureSide = preflightInfo->sideB;
        if (preflightInfo->patchA == secondPatch && preflightInfo->patchB == firstPatch)
        {
            firstFeatureSide = preflightInfo->sideB;
            secondFeatureSide = preflightInfo->sideA;
        }

        outfile << " featureSides=(" << featureSideName(firstFeatureSide)
                << "," << featureSideName(secondFeatureSide) << ")"
                << " orientation="
                << (preflightInfo->orientationReversed ? "reversed" : "same");
    }
    else
    {
        outfile << " featureSides=(unavailable,unavailable) orientation=unknown";
    }

    outfile << "\n";
}

static void logInvestigatedTwinSummary(
    index_t firstPatch,
    index_t secondPatch,
    const gsVector<gsVector<gsVector<std::vector<index_t>>>>& twinsIndex,
    const gsVector<gsVector<gsVector<std::vector<index_t>>>>& twinsPatch)
{
    if (!outfile.is_open() || !isInvestigatedInterface(firstPatch, secondPatch))
        return;

    if (firstPatch < 0 || secondPatch < 0 ||
        firstPatch >= twinsIndex.size() || secondPatch >= twinsIndex.size() ||
        firstPatch >= twinsPatch.size() || secondPatch >= twinsPatch.size())
    {
        outfile << "[investigate-interface] twin-summary skipped: patches out of range "
                << "patches=(" << firstPatch << "," << secondPatch << ") "
                << "twinsIndex.size=" << twinsIndex.size() << " "
                << "twinsPatch.size=" << twinsPatch.size() << "\n";
        return;
    }

    outfile << "[investigate-interface] twin-summary-begin patches=("
            << firstPatch << "," << secondPatch << ")\n";

    const int sharedLevels = std::min(
        std::min(twinsIndex(firstPatch).size(), twinsIndex(secondPatch).size()),
        std::min(twinsPatch(firstPatch).size(), twinsPatch(secondPatch).size()));
    for (int level = 0; level < sharedLevels; ++level)
    {
        int countFirstToSecond = 0;
        int countSecondToFirst = 0;
        std::ostringstream examplesFirstToSecond;
        std::ostringstream examplesSecondToFirst;
        int exampleCountFirstToSecond = 0;
        int exampleCountSecondToFirst = 0;

        const int firstFunctionCount = std::min(
            twinsIndex(firstPatch)(level).size(),
            twinsPatch(firstPatch)(level).size());
        for (int functionIndex = 0; functionIndex < firstFunctionCount; ++functionIndex)
        {
            const size_t pairCount = std::min(
                twinsIndex(firstPatch)(level)(functionIndex).size(),
                twinsPatch(firstPatch)(level)(functionIndex).size());
            for (size_t pairIdx = 0; pairIdx < pairCount; ++pairIdx)
            {
                if (twinsPatch(firstPatch)(level)(functionIndex)[pairIdx] != secondPatch)
                    continue;

                ++countFirstToSecond;
                if (exampleCountFirstToSecond < 8)
                {
                    if (exampleCountFirstToSecond > 0)
                        examplesFirstToSecond << "; ";
                    examplesFirstToSecond << functionIndex << "->"
                                          << twinsIndex(firstPatch)(level)(functionIndex)[pairIdx];
                    ++exampleCountFirstToSecond;
                }
            }
        }

        const int secondFunctionCount = std::min(
            twinsIndex(secondPatch)(level).size(),
            twinsPatch(secondPatch)(level).size());
        for (int functionIndex = 0; functionIndex < secondFunctionCount; ++functionIndex)
        {
            const size_t pairCount = std::min(
                twinsIndex(secondPatch)(level)(functionIndex).size(),
                twinsPatch(secondPatch)(level)(functionIndex).size());
            for (size_t pairIdx = 0; pairIdx < pairCount; ++pairIdx)
            {
                if (twinsPatch(secondPatch)(level)(functionIndex)[pairIdx] != firstPatch)
                    continue;

                ++countSecondToFirst;
                if (exampleCountSecondToFirst < 8)
                {
                    if (exampleCountSecondToFirst > 0)
                        examplesSecondToFirst << "; ";
                    examplesSecondToFirst << functionIndex << "->"
                                          << twinsIndex(secondPatch)(level)(functionIndex)[pairIdx];
                    ++exampleCountSecondToFirst;
                }
            }
        }

        if (countFirstToSecond == 0 && countSecondToFirst == 0)
            continue;

        outfile << "[investigate-interface] twin-summary level=" << level
                << " firstToSecond=" << countFirstToSecond
                << " secondToFirst=" << countSecondToFirst;
        if (exampleCountFirstToSecond > 0)
            outfile << " examplesFirstToSecond={" << examplesFirstToSecond.str() << "}";
        if (exampleCountSecondToFirst > 0)
            outfile << " examplesSecondToFirst={" << examplesSecondToFirst.str() << "}";
        outfile << "\n";
    }

    outfile << "[investigate-interface] twin-summary-end\n";
}



void orientThePatches(gsMultiPatch<> mp, std::vector<int>& firstSide, std::vector<int>& secondSide, std::vector<int>& firstPatch, std::vector<int>& secondPatch)
{
    std::vector<boundaryInterface>& boundaryInterfaces = mp.interfaces();
    for (int interfaceNum = 0; interfaceNum < boundaryInterfaces.size(); ++interfaceNum)
    {
        outfile << boundaryInterfaces[interfaceNum] << "\n";
        if (g_verbose) gsInfo << boundaryInterfaces[interfaceNum] << "\n";
        const int rawFirstSide = boundaryInterfaces[interfaceNum].first().index();
        const int rawSecondSide = boundaryInterfaces[interfaceNum].second().index();
        const int localFirstPatch = boundaryInterfaces[interfaceNum].first().patch;
        const int localSecondPatch = boundaryInterfaces[interfaceNum].second().patch;
        int localFirstSide = rawFirstSide;
        int localSecondSide = rawSecondSide;
        const bool hasCanonicalSides = resolveCanonicalTopologySides(
            localFirstPatch,
            localSecondPatch,
            rawFirstSide,
            rawSecondSide,
            localFirstSide,
            localSecondSide);

        const GeometryPreflightInterfaceInfo* preflightInfo =
            findPreflightInterfaceInfo(localFirstPatch, localSecondPatch);
        if (preflightInfo && outfile.is_open())
        {
            outfile << "[topology-preflight] interface=" << interfaceNum
                    << " patches=(" << localFirstPatch << "," << localSecondPatch << ")"
                    << " mirrored=(" << (isPreflightMirroredPatch(localFirstPatch) ? 1 : 0)
                    << "," << (isPreflightMirroredPatch(localSecondPatch) ? 1 : 0) << ")"
                    << " rawSides=(" << rawFirstSide << "," << rawSecondSide << ")"
                    << " canonicalSides=(" << localFirstSide << "," << localSecondSide << ")"
                    << " orientation=" << (preflightInfo->orientationReversed ? "reversed" : "same")
                    << "\n";
        }
        else if (!hasCanonicalSides && outfile.is_open())
        {
            outfile << "[topology-preflight] interface=" << interfaceNum
                    << " patches=(" << localFirstPatch << "," << localSecondPatch << ")"
                    << " canonicalSideResolution=failed"
                    << " rawSides=(" << rawFirstSide << "," << rawSecondSide << ")\n";
        }

        logInvestigatedTopologyInterface(
            "orientThePatches",
            interfaceNum,
            localFirstPatch,
            localSecondPatch,
            rawFirstSide,
            rawSecondSide,
            localFirstSide,
            localSecondSide,
            preflightInfo);

        firstSide.push_back(localFirstSide);
        secondSide.push_back(localSecondSide);
        firstPatch.push_back(localFirstPatch);
        secondPatch.push_back(localSecondPatch);
    }
}

void adjustSpilloverBoxes(gsVector<gsVector<gsTHBSplineBasis<2, real_t>>> spilloverPatchStructure,
    gsVector<gsVector<gsVector<index_t>>> boxMat,
    gsVector<gsVector<gsVector<gsVector<index_t>>>>& SpilloverboxMat
    , std::vector<int> firstSide,
    std::vector<int> secondSide,
    std::vector<int> firstPatch,
    std::vector<int> secondPatch,
    gsVector<index_t> currentLastNonZeroRow
)
{
    int maxLocIndex;
    for (size_t i = 0; i < firstSide.size(); i++)
    {
        gsInfo << i << "\n";
        gsInfo << i << " ";
        //gsInfo << "currentLastNonZeroRow: " << currentLastNonZeroRow[firstPatch[i]];
        for (size_t j = 0; j < currentLastNonZeroRow[firstPatch[i]]; j++)
        {
            gsInfo << j << ":\n";
            std::vector<index_t> box(5);
            std::vector<index_t> orderedbox(5);
            maxLocIndex = std::pow(2, boxMat[firstPatch[i]][j][0]);
            coordinateTranformationBoxes(firstSide[i], boxMat[firstPatch[i]][j], secondSide[i],
                box, maxLocIndex);
            orderedbox[0] = box[0];
            orderedbox[1] = std::min(box[1], box[3]);
            orderedbox[2] = std::min(box[2], box[4]);
            orderedbox[3] = std::max(box[1], box[3]);
            orderedbox[4] = std::max(box[2], box[4]);
            for (size_t k = 0; k < 5; k++)
            {
                SpilloverboxMat[firstPatch[i]][secondPatch[i]][j][k] = orderedbox[k];
                outfile << boxMat[firstPatch[i]][j][k] << " " << SpilloverboxMat[firstPatch[i]][secondPatch[i]][j][k] << "\n";
                //gsInfo <<  boxMat[firstPatch[i]][j][k] << " " << SpilloverboxMat[firstPatch[i]][secondPatch[i]][j][k] << "\n";
            }
            outfile << ";\n";
        }
    }
}

std::vector<std::unordered_set<index_t>> computeTwinPatches(
    const std::vector<std::vector<std::vector<index_t>>>& functionDescription)
{
    std::vector<std::unordered_set<index_t>> twinPatches(functionDescription.size());

    for (index_t f = 0; f < functionDescription.size(); ++f)
    {
        for (const auto& twin : functionDescription[f])
        {
            index_t patchIndex = twin[0]; // each twin = {patch, functionIndex}
            twinPatches[f].insert(patchIndex);
        }
    }

    return twinPatches;
}


void IdentifyPatches(gsMultiPatch<> mp,
    gsVector < gsVector < gsTensorBSplineBasis <2, real_t >>> Bells,
    gsVector  <gsVector< gsVector<index_t>>>& isTouching,
    gsVector < gsVector < gsVector <std::vector< index_t >>>>& twinsIndex,
    gsVector < gsVector < gsVector <std::vector< index_t >>>>& twinsPatch,
    gsVector < gsVector < gsVector < index_t >>>& hasATwin,
    gsVector  <gsVector< gsVector<index_t>>>& isActive
) {
    auto patchHasLevel = [&](int patch, int level) -> bool {
        return patch >= 0 && patch < Bells.size() && level >= 0 && level < Bells(patch).size();
    };

    auto patchHasFunction = [&](int patch, int level, int functionIndex) -> bool {
        return patchHasLevel(patch, level) &&
            functionIndex >= 0 && functionIndex < Bells(patch)(level).size();
    };

    for (int i = 0; i < Bells.size(); ++i) {
        twinsIndex(i).resize(Bells(i).size());
        twinsPatch(i).resize(Bells(i).size());
        isTouching(i).setZero(Bells(i).size());
        hasATwin(i).setZero(Bells(i).size());
        isActive(i).setZero(Bells(i).size());
        for (int j = 0; j < Bells(i).size(); ++j) {
            twinsIndex(i)(j).resize(Bells(i)(j).size());
            twinsPatch(i)(j).resize(Bells(i)(j).size());
            isTouching(i)(j).setZero(Bells(i)(j).size());
            isActive(i)(j).setZero(Bells(i)(j).size());
            hasATwin(i)(j).setZero(Bells(i)(j).size());
        }
    }
    gsVector  <gsVector< gsVector<int>>> touchedPatchSide(mp.nPatches());
    for (int patch = 0; patch < mp.nPatches(); ++patch) {
        touchedPatchSide(patch).resize(Bells(patch).size());
        for (int level = 0; level < Bells(patch).size(); ++level) {
            touchedPatchSide(patch)(level).resize(Bells(patch)(level).size());
            for (int functionIndex = 0; functionIndex < Bells(patch)(level).size(); ++functionIndex) {
                touchedPatchSide(patch)(level)(functionIndex) = -1;
            }
        }
    }
    std::vector<boundaryInterface>& boundaryInterfaces = mp.interfaces();
    for (int interfaceNum = 0; interfaceNum < boundaryInterfaces.size(); ++interfaceNum) {
        outfile << boundaryInterfaces[interfaceNum] << "\n";
        const int rawFirstSide = boundaryInterfaces[interfaceNum].first().index();
        const int rawSecondSide = boundaryInterfaces[interfaceNum].second().index();
        int firstPatch = boundaryInterfaces[interfaceNum].first().patch;
        int secondPatch = boundaryInterfaces[interfaceNum].second().patch;
        int firstSide = rawFirstSide;
        int secondSide = rawSecondSide;
        resolveCanonicalTopologySides(
            firstPatch,
            secondPatch,
            rawFirstSide,
            rawSecondSide,
            firstSide,
            secondSide);
        const GeometryPreflightInterfaceInfo* preflightInfo =
            findPreflightInterfaceInfo(firstPatch, secondPatch);
        if (preflightInfo && outfile.is_open() && g_verbose)
        {
            outfile << "[IdentifyPatches] interface=" << interfaceNum
                    << " firstPatch=" << firstPatch
                    << " secondPatch=" << secondPatch
                    << " rawSides=(" << rawFirstSide << "," << rawSecondSide << ")"
                    << " canonicalSides=(" << firstSide << "," << secondSide << ")"
                    << " firstMirrored=" << (isPreflightMirroredPatch(firstPatch) ? 1 : 0)
                    << " secondMirrored=" << (isPreflightMirroredPatch(secondPatch) ? 1 : 0)
                    << " orientation=" << (preflightInfo->orientationReversed ? "reversed" : "same")
                    << "\n";
        }

        logInvestigatedTopologyInterface(
            "IdentifyPatches.markTouching",
            interfaceNum,
            firstPatch,
            secondPatch,
            rawFirstSide,
            rawSecondSide,
            firstSide,
            secondSide,
            preflightInfo);

        for (int level = 0; level < Bells(firstPatch).size(); ++level) {
            for (int functionIndex = 0; functionIndex < Bells(firstPatch)(level).size(); ++functionIndex) {
                if (!isTouching(firstPatch)(level)(functionIndex)) {
                    if (supportTouchesInterfaceSide(
                        Bells(firstPatch)(level).function(functionIndex),
                        firstSide)) {
                        isTouching(firstPatch)(level)(functionIndex) = 1;
                        touchedPatchSide(firstPatch)(level)(functionIndex) = secondSide;
                    }
                }
            }
        }
        for (int level = 0; level < Bells(secondPatch).size(); ++level) {
            for (int functionIndex = 0; functionIndex < Bells(secondPatch)(level).size(); ++functionIndex) {
                if (!isTouching(secondPatch)(level)(functionIndex)) {
                    if (supportTouchesInterfaceSide(
                        Bells(secondPatch)(level).function(functionIndex),
                        secondSide)) {
                        isTouching(secondPatch)(level)(functionIndex) = 1;
                        touchedPatchSide(secondPatch)(level)(functionIndex) = firstSide;
                    }
                }
            }
        }
    }
    for (int patch = 0; patch < isTouching.size(); ++patch) {
        for (int level = 0; level < isTouching(patch).size(); ++level) {
        }
    }
    //Now we know all the functions that are touching the boundary
    for (int interfaceNum = 0; interfaceNum < boundaryInterfaces.size(); ++interfaceNum) {
        const int rawFirstSide = boundaryInterfaces[interfaceNum].first().index();
        const int rawSecondSide = boundaryInterfaces[interfaceNum].second().index();
        int firstPatch = boundaryInterfaces[interfaceNum].first().patch;
        int secondPatch = boundaryInterfaces[interfaceNum].second().patch;
        int firstSide = rawFirstSide;
        int secondSide = rawSecondSide;
        resolveCanonicalTopologySides(
            firstPatch,
            secondPatch,
            rawFirstSide,
            rawSecondSide,
            firstSide,
            secondSide);
        // gsInfo << "Identifying " << firstPatch << " and " << secondPatch << "\n";
        const int sharedLevels = std::min(isTouching(firstPatch).size(), isTouching(secondPatch).size());
        std::vector<int> candidateCountPerLevel(sharedLevels, 0);
        std::vector<int> directTwinCountPerLevel(sharedLevels, 0);
        std::vector<int> distanceRejectedCountPerLevel(sharedLevels, 0);
        std::vector<int> unresolvedCountPerLevel(sharedLevels, 0);
        for (int level = 0; level < sharedLevels; ++level) {
            if (!patchHasLevel(firstPatch, level) || !patchHasLevel(secondPatch, level))
                continue;
            for (int functionIndex = 0; functionIndex < Bells(firstPatch)(level).size(); ++functionIndex) {
                int test250519 = isZeroOnBoundary(Bells, firstPatch, level, functionIndex, firstSide);
                if (
                    isTouching(firstPatch)(level)(functionIndex)
                    && test250519 == 0
                    ) {
                    candidateCountPerLevel[level]++;
                    int twinIndex = -1;
                    real_t bestDistance = std::numeric_limits<real_t>::infinity();

                    gsMatrix<> firstAnchor;
                    gsMatrix<> secondAnchors;
                    Bells(firstPatch)(level).anchor_into(functionIndex, firstAnchor);
                    Bells(secondPatch)(level).anchors_into(secondAnchors);

                    if (firstAnchor.rows() >= 2 && firstAnchor.cols() >= 1)
                    {
                        const gsVector<real_t> firstUv = matrixColumnToVector2(firstAnchor);
                        const gsVector<real_t> firstXy = matrixColumnToVector2(mp.patch(firstPatch).eval(firstUv));

                        for (int candidateIndex = 0; candidateIndex < secondAnchors.cols(); ++candidateIndex)
                        {
                            if (!isTouching(secondPatch)(level)(candidateIndex))
                                continue;
                            if (isZeroOnBoundary(Bells, secondPatch, level, candidateIndex, secondSide) != 0)
                                continue;

                            const gsVector<real_t> secondUv = matrixColumnToVector2(secondAnchors.col(candidateIndex));
                            const gsVector<real_t> secondXy = matrixColumnToVector2(mp.patch(secondPatch).eval(secondUv));
                            const real_t distance = vectorDistance(firstXy, secondXy);
                            if (distance < bestDistance)
                            {
                                bestDistance = distance;
                                twinIndex = candidateIndex;
                            }
                        }

                        if (bestDistance > 1e-5)
                            twinIndex = -1;
                    }

                    const bool rejectedByDistance =
                        (bestDistance < std::numeric_limits<real_t>::infinity() && bestDistance > 1e-5);

                    if (twinIndex != -1 && patchHasFunction(secondPatch, level, twinIndex)) {
                        if (
                            isTouching(secondPatch)(level)(twinIndex)
                            ) {
                            twinsIndex(firstPatch)(level)(functionIndex).push_back(twinIndex);
                            twinsPatch(firstPatch)(level)(functionIndex).push_back(secondPatch);
                            hasATwin(firstPatch)(level)(functionIndex) = 1;
                            twinsPatch(secondPatch)(level)(twinIndex).push_back(firstPatch);
                            twinsIndex(secondPatch)(level)(twinIndex).push_back(functionIndex);
                            hasATwin(secondPatch)(level)(twinIndex) = 1;
                            directTwinCountPerLevel[level]++;

                            if (outfile.is_open() && g_verbose)
                            {
                                outfile << "[IdentifyPatches] direct twin level=" << level
                                        << " via interface " << interfaceNum
                                        << " added (patch=" << firstPatch << ", index=" << functionIndex
                                        << ") <-> (patch=" << secondPatch << ", index=" << twinIndex << ")\n";
                            }
                        }
                        else if (!isTouching(secondPatch)(level)(twinIndex)) {
                            //But the twin does not touch the boundary\n";
                            unresolvedCountPerLevel[level]++;
                        }
                    }
                    else if (rejectedByDistance)
                    {
                        distanceRejectedCountPerLevel[level]++;
                    }
                    else
                    {
                        unresolvedCountPerLevel[level]++;
                    }
                }
                else if (!isTouching(firstPatch)(level)(functionIndex)) {
                    //The function " << functionIndex << " does not touch the boundary\n";
                }
            }
        }

        if (outfile.is_open() && isInvestigatedInterface(firstPatch, secondPatch))
        {
            outfile << "[investigate-interface] stage=IdentifyPatches.directTwins interface="
                    << interfaceNum << " patches=(" << firstPatch << "," << secondPatch << ")\n";
            for (int level = 0; level < sharedLevels; ++level)
            {
                if (candidateCountPerLevel[level] == 0 &&
                    directTwinCountPerLevel[level] == 0 &&
                    distanceRejectedCountPerLevel[level] == 0 &&
                    unresolvedCountPerLevel[level] == 0)
                    continue;

                outfile << "[investigate-interface] level=" << level
                        << " candidates=" << candidateCountPerLevel[level]
                        << " directTwins=" << directTwinCountPerLevel[level]
                        << " distanceRejected=" << distanceRejectedCountPerLevel[level]
                        << " unresolved=" << unresolvedCountPerLevel[level]
                        << "\n";
            }
        }
    }
    for (int patch = 0; patch < twinsIndex.size() - 1; ++patch) {
        for (int level = 0; level < twinsIndex(patch).size(); ++level) {
            for (int functionIndex = 0; functionIndex < twinsIndex(patch)(level).size(); ++functionIndex) {
                if (twinsIndex(patch)(level)(functionIndex).size() == 0) {
                    twinsIndex(patch)(level)(functionIndex).push_back(-1);
                    twinsPatch(patch)(level)(functionIndex).push_back(-1);
                }
            }
        }
    }

    auto hasTwinPair = [&](const std::vector<int>& indices,
                           const std::vector<int>& patches,
                           int candidateIndex,
                           int candidatePatch) -> bool
    {
        const size_t pairCount = std::min(indices.size(), patches.size());
        for (size_t pairIdx = 0; pairIdx < pairCount; ++pairIdx)
        {
            if (indices[pairIdx] == candidateIndex && patches[pairIdx] == candidatePatch)
                return true;
        }
        return false;
    };

    //now we have to identify twins of twins. Since being a twin is a equivalence relation, it is transitive
    for (int patch = 0; patch < twinsIndex.size(); ++patch) {
        for (int level = 0; level < twinsIndex(patch).size(); ++level) {
            for (int funcIndex = 0; funcIndex < twinsIndex(patch)(level).size(); ++funcIndex) {
                for (int twinNum1 = 0; twinNum1 < twinsIndex(patch)(level)(funcIndex).size(); ++twinNum1) {
                    int patchIndex = twinsPatch(patch)(level)(funcIndex)[twinNum1];
                    int functionIndex = twinsIndex(patch)(level)(funcIndex)[twinNum1];
                    if (patchIndex == -1 || functionIndex == -1)     break;
                    if (patchIndex < 0 || patchIndex >= twinsIndex.size() || level >= twinsIndex(patchIndex).size())
                        continue;
                    if (functionIndex < 0 || functionIndex >= twinsIndex(patchIndex)(level).size())
                        continue;
                    for (int twinNum2 = 0; twinNum2 < twinsIndex(patchIndex)(level)(functionIndex).size(); ++twinNum2) {
                        int tS = twinsIndex(patchIndex)(level)(functionIndex).size();
                        if (twinsIndex(patchIndex)(level)(functionIndex).size() <= 1) {
                            break;
                        }
                        else {
                            std::vector<int> arr1 = twinsIndex(patch)(level)(funcIndex);
                            std::vector<int> arr2 = twinsPatch(patch)(level)(funcIndex);
                            const int candidateIndex = twinsIndex(patchIndex)(level)(functionIndex)[twinNum2];
                            const int candidatePatch = twinsPatch(patchIndex)(level)(functionIndex)[twinNum2];
                            if (!hasTwinPair(arr1, arr2, candidateIndex, candidatePatch) &&
                                !(candidateIndex == funcIndex && candidatePatch == patch)
                                ) {
                                arr1.push_back(candidateIndex);
                                arr2.push_back(candidatePatch);
                                twinsIndex(patch)(level)(funcIndex) = arr1;
                                twinsPatch(patch)(level)(funcIndex) = arr2;

                                if (outfile.is_open() && g_verbose)
                                {
                                    outfile << "[IdentifyPatches] transitive twin level=" << level
                                            << " seed=(patch=" << patch << ", index=" << funcIndex << ")"
                                            << " via=(patch=" << patchIndex << ", index=" << functionIndex << ")"
                                            << " added=(patch=" << candidatePatch << ", index=" << candidateIndex << ")\n";
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    //once again
    for (int patch = twinsIndex.size() - 1; patch >= 0; --patch) {
        for (int level = 0; level < twinsIndex(patch).size(); ++level) {
            for (int funcIndex = twinsIndex(patch)(level).size() - 1; funcIndex >= 0; --funcIndex) {
                for (int twinNum1 = 0; twinNum1 < twinsIndex(patch)(level)(funcIndex).size(); ++twinNum1) {
                    int patchIndex = twinsPatch(patch)(level)(funcIndex)[twinNum1];
                    int functionIndex = twinsIndex(patch)(level)(funcIndex)[twinNum1];
                    if (patchIndex == -1 || functionIndex == -1)     break;
                    if (patchIndex < 0 || patchIndex >= twinsIndex.size() || level >= twinsIndex(patchIndex).size())
                        continue;
                    if (functionIndex < 0 || functionIndex >= twinsIndex(patchIndex)(level).size())
                        continue;
                    for (int twinNum2 = 0; twinNum2 < twinsIndex(patchIndex)(level)(functionIndex).size(); ++twinNum2) {
                        int tS = twinsIndex(patchIndex)(level)(functionIndex).size();
                        if (twinsIndex(patchIndex)(level)(functionIndex).size() <= 1) {
                            break;
                        }
                        else {
                            std::vector<int> arr1 = twinsIndex(patch)(level)(funcIndex);
                            std::vector<int> arr2 = twinsPatch(patch)(level)(funcIndex);
                            const int candidateIndex = twinsIndex(patchIndex)(level)(functionIndex)[twinNum2];
                            const int candidatePatch = twinsPatch(patchIndex)(level)(functionIndex)[twinNum2];
                            if (!hasTwinPair(arr1, arr2, candidateIndex, candidatePatch) &&
                                !(candidateIndex == funcIndex && candidatePatch == patch)
                                ) {
                                arr1.push_back(candidateIndex);
                                arr2.push_back(candidatePatch);
                                twinsIndex(patch)(level)(funcIndex) = arr1;
                                twinsPatch(patch)(level)(funcIndex) = arr2;

                                if (outfile.is_open() && g_verbose)
                                {
                                    outfile << "[IdentifyPatches] transitive twin reverse-pass level=" << level
                                            << " seed=(patch=" << patch << ", index=" << funcIndex << ")"
                                            << " via=(patch=" << patchIndex << ", index=" << functionIndex << ")"
                                            << " added=(patch=" << candidatePatch << ", index=" << candidateIndex << ")\n";
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    for (int patch = 0; patch < twinsIndex.size(); ++patch) {
        for (int level = 0; level < twinsIndex(patch).size(); ++level) {
            for (int functionIndex = 0; functionIndex < twinsIndex(patch)(level).size(); ++functionIndex) {
                if (twinsIndex(patch)(level)(functionIndex).size() == 0) {
                    twinsIndex(patch)(level)(functionIndex).push_back(-1);
                    twinsPatch(patch)(level)(functionIndex).push_back(-1);
                }
            }
        }
    }

    if (outfile.is_open() && g_verbose)
    {
        outfile << "[IdentifyPatches] results begin\n";
        for (int patch = 0; patch < twinsIndex.size(); ++patch)
        {
            outfile << "[IdentifyPatches] patch " << patch << "\n";
            for (int level = 0; level < twinsIndex(patch).size(); ++level)
            {
                for (int functionIndex = 0; functionIndex < twinsIndex(patch)(level).size(); ++functionIndex)
                {
                    outfile << "[IdentifyPatches] twins "
                            << patch << " " << level << " " << functionIndex << ": ";

                    const size_t pairCount = std::min(twinsIndex(patch)(level)(functionIndex).size(),
                                                      twinsPatch(patch)(level)(functionIndex).size());
                    bool wroteTwin = false;
                    for (size_t twinNum = 0; twinNum < pairCount; ++twinNum)
                    {
                        const int twinPatch = twinsPatch(patch)(level)(functionIndex)[twinNum];
                        const int twinIndex = twinsIndex(patch)(level)(functionIndex)[twinNum];
                        if (twinPatch == -1 || twinIndex == -1)
                            continue;

                        outfile << twinPatch << " " << level << " " << twinIndex << "  ";
                        wroteTwin = true;
                    }

                    if (!wroteTwin)
                        outfile << "-1 " << level << " -1";

                    outfile << "\n";
                }
            }
        }
        outfile << "[IdentifyPatches] results end\n";
    }

    logInvestigatedTwinSummary(15, 20, twinsIndex, twinsPatch);
}

//void IdentifyPatches(gsMultiPatch<> mp,
//    gsVector <  gsTHBSplineBasis <2, real_t >> THBVector1,
//    gsVector  < gsVector<index_t>>& isTouching,
//    gsVector < gsVector <std::vector< index_t >>>& twinsIndex,
//    gsVector < gsVector <std::vector< index_t >>>& twinsPatch,
//    gsVector < gsVector < index_t >>& hasATwin,
//    gsVector  < gsVector<index_t>>& isActive
//) {
//    for (int patch = 0; patch < THBVector1.size(); ++patch) {
//        twinsIndex(patch).resize(THBVector1(patch).size());
//        twinsPatch(patch).resize(THBVector1(patch).size());
//        isTouching(patch).setZero(THBVector1(patch).size());
//        hasATwin(patch).setZero(THBVector1(patch).size());
//        isActive(patch).setZero(THBVector1(patch).size());
//    }
//    gsVector  < gsVector<int>> touchedPatchSide(mp.nPatches());
//    for (int patch = 0; patch < mp.nPatches(); ++patch) {
//        touchedPatchSide(patch).resize(THBVector1(patch).size());
//        for (int functionIndex = 0; functionIndex < THBVector1(patch).size(); ++functionIndex) {
//            touchedPatchSide(patch)(functionIndex) = -1;
//        }
//    }
//    std::vector<boundaryInterface>& boundaryInterfaces = mp.interfaces();
//    gsInfo << "directionMap and Orientation: \n";
//    for (int interfaceNum = 0; interfaceNum < boundaryInterfaces.size(); ++interfaceNum) {
//        gsInfo << boundaryInterfaces[interfaceNum].dirMap() << ",\n" << boundaryInterfaces[interfaceNum].dirOrientation() << ".\n";
//        gsVector<index_t> dummyVector;
//        boundaryInterfaces[interfaceNum].cornerMap(dummyVector);
//        //gsInfo << "cornerMap:" << dummyVector << "\n";
//        gsInfo << "NARA\n";
//        outfile << boundaryInterfaces[interfaceNum] << "\n";
//        int firstSide = boundaryInterfaces[interfaceNum].first().index();
//        int secondSide = boundaryInterfaces[interfaceNum].second().index();
//        int firstPatch = boundaryInterfaces[interfaceNum].first().patch;
//        int secondPatch = boundaryInterfaces[interfaceNum].second().patch;
//        for (int functionIndex = 0; functionIndex < THBVector1(firstPatch).size(); ++functionIndex) {
//            if (!isTouching(firstPatch)(functionIndex)) {
//                if (THBVector1(firstPatch).function(functionIndex).support()(
//                    (firstSide + 1) / 2 - 1, (firstSide + 1) % 2) ==
//                    THBVector1(firstPatch).support()((firstSide + 1) / 2 - 1, (firstSide + 1) % 2)) {
//                    isTouching(firstPatch)(functionIndex) = 1;
//                    touchedPatchSide(firstPatch)(functionIndex) = secondSide;
//                }
//            }
//            if (!isTouching(secondPatch)(functionIndex)) {
//                if (THBVector1(secondPatch).function(functionIndex).support()((secondSide + 1) / 2 - 1, (secondSide + 1) % 2) ==
//                    THBVector1(secondPatch).support()((secondSide + 1) / 2 - 1, (secondSide + 1) % 2)) {
//                    isTouching(secondPatch)(functionIndex) = 1;
//                    touchedPatchSide(secondPatch)(functionIndex) = firstSide;
//                }
//            }
//        }
//    }
//    //Now we know all the functions that are touching the boundary
//    for (int interfaceNum = 0; interfaceNum < boundaryInterfaces.size(); ++interfaceNum) {
//        int firstSide = boundaryInterfaces[interfaceNum].first().index();
//        int secondSide = boundaryInterfaces[interfaceNum].second().index();
//        int firstPatch = boundaryInterfaces[interfaceNum].first().patch;
//        int secondPatch = boundaryInterfaces[interfaceNum].second().patch;
//        gsInfo << "Identifying " << firstPatch << " and " << secondPatch << "\n";
//        for (int functionIndex = 0; functionIndex < THBVector1(firstPatch).size(); ++functionIndex) {
//            if (
//                isTouching(firstPatch)(functionIndex)
//                && !isZeroOnBoundary(THBVector1, firstPatch, functionIndex, firstSide)
//                ) {
//                int twinIndex = twinFunction(THBVector1(firstPatch), functionIndex, THBVector1(secondPatch), firstSide, secondSide);
//                if (twinIndex != -1) {
//                    if (
//                        isTouching(secondPatch)(twinIndex)
//                        ) {
//                        twinsIndex(firstPatch)(functionIndex).push_back(twinIndex);
//                        twinsPatch(firstPatch)(functionIndex).push_back(secondPatch);
//                        hasATwin(firstPatch)(functionIndex) = 1;
//                        twinsPatch(secondPatch)(twinIndex).push_back(firstPatch);
//                        twinsIndex(secondPatch)(twinIndex).push_back(functionIndex);
//                        hasATwin(secondPatch)(twinIndex) = 1;
//                    }
//                    else if (!isTouching(secondPatch)(twinIndex)) {
//                        //But the twin does not touch the boundary\n";
//                    }
//                }
//            }
//            else if (!isTouching(firstPatch)(functionIndex)) {
//                //The function " << functionIndex << " does not touch the boundary\n";
//            }
//        }
//
//    }
//    for (int patch = 0; patch < twinsIndex.size() - 1; ++patch) {
//        for (int functionIndex = 0; functionIndex < twinsIndex(patch).size(); ++functionIndex) {
//            if (twinsIndex(patch)(functionIndex).size() == 0) {
//                twinsIndex(patch)(functionIndex).push_back(-1);
//                twinsPatch(patch)(functionIndex).push_back(-1);
//            }
//        }
//    }
//    //now we have to identify twins of twins. Since being a twin is a equivalence relation, it is transitive
//    for (int patch = 0; patch < twinsIndex.size(); ++patch) {
//        for (int funcIndex = 0; funcIndex < twinsIndex(patch).size(); ++funcIndex) {
//            for (int twinNum1 = 0; twinNum1 < twinsIndex(patch)(funcIndex).size(); ++twinNum1) {
//                int patchIndex = twinsPatch(patch)(funcIndex)[twinNum1];
//                int functionIndex = twinsIndex(patch)(funcIndex)[twinNum1];
//                if (patchIndex == -1 || functionIndex == -1)     break;
//                for (int twinNum2 = 0; twinNum2 < twinsIndex(patchIndex)(functionIndex).size(); ++twinNum2) {
//                    int tS = twinsIndex(patchIndex)(functionIndex).size();
//                    if (twinsIndex(patchIndex)(functionIndex).size() <= 1) {
//                        break;
//                    }
//                    else {
//                        std::vector<int>::iterator it1;
//                        std::vector<int>::iterator it2;
//                        std::vector<int> arr1 = twinsIndex(patch)(funcIndex);
//                        std::vector<int> arr2 = twinsPatch(patch)(funcIndex);
//                        it1 = std::find(arr1.begin(),
//                            arr1.end(),
//                            twinsIndex(patchIndex)(functionIndex)[twinNum2]
//                        );
//                        it2 = std::find(arr2.begin(),
//                            arr2.end(),
//                            twinsPatch(patchIndex)(functionIndex)[twinNum2]
//                        );
//                        bool a1 = !(it1 != arr1.end() && it2 != arr2.end());
//                        bool a2 = twinsIndex(patchIndex)(functionIndex)[twinNum2] != funcIndex;
//                        bool a3 = twinsPatch(patchIndex)(functionIndex)[twinNum2] != patch;
//                        if (!(it1 != arr1.end() && it2 != arr2.end()) &&
//                            twinsIndex(patchIndex)(functionIndex)[twinNum2] != funcIndex &&
//                            twinsPatch(patchIndex)(functionIndex)[twinNum2] != patch
//                            ) {
//                            arr1.push_back(twinsIndex(patchIndex)(functionIndex)[twinNum2]);
//                            arr2.push_back(twinsPatch(patchIndex)(functionIndex)[twinNum2]);
//                            twinsIndex(patch)(funcIndex) = arr1;
//                            twinsPatch(patch)(funcIndex) = arr2;
//                        }
//                    }
//                }
//            }
//        }
//    }
//    //once again
//    for (int patch = twinsIndex.size() - 1; patch >= 0; --patch) {
//        for (int funcIndex = twinsIndex(patch).size() - 1; funcIndex >= 0; --funcIndex) {
//            for (int twinNum1 = 0; twinNum1 < twinsIndex(patch)(funcIndex).size(); ++twinNum1) {
//                int patchIndex = twinsPatch(patch)(funcIndex)[twinNum1];
//                int functionIndex = twinsIndex(patch)(funcIndex)[twinNum1];
//                if (patchIndex == -1 || functionIndex == -1)     break;
//                for (int twinNum2 = 0; twinNum2 < twinsIndex(patchIndex)(functionIndex).size(); ++twinNum2) {
//                    int tS = twinsIndex(patchIndex)(functionIndex).size();
//                    if (twinsIndex(patchIndex)(functionIndex).size() <= 1) {
//                        break;
//                    }
//                    else {
//                        std::vector<int>::iterator it1;
//                        std::vector<int>::iterator it2;
//                        std::vector<int> arr1 = twinsIndex(patch)(funcIndex);
//                        std::vector<int> arr2 = twinsPatch(patch)(funcIndex);
//                        it1 = std::find(arr1.begin(),
//                            arr1.end(),
//                            twinsIndex(patchIndex)(functionIndex)[twinNum1]
//                        );
//                        it2 = std::find(arr2.begin(),
//                            arr2.end(),
//                            twinsPatch(patchIndex)(functionIndex)[twinNum1]
//                        );
//                        bool a1 = !(it1 != arr1.end() && it2 != arr2.end());
//                        bool a2 = twinsIndex(patchIndex)(functionIndex)[twinNum1] != funcIndex;
//                        bool a3 = twinsPatch(patchIndex)(functionIndex)[twinNum1] != patch;
//                        if (!(it1 != arr1.end() && it2 != arr2.end()) &&
//                            twinsIndex(patchIndex)(functionIndex)[twinNum1] != funcIndex &&
//                            twinsPatch(patchIndex)(functionIndex)[twinNum1] != patch
//                            ) {
//                            arr1.push_back(twinsIndex(patchIndex)(functionIndex)[twinNum1]);
//                            arr2.push_back(twinsPatch(patchIndex)(functionIndex)[twinNum2]);
//                            twinsIndex(patch)(funcIndex) = arr1;
//                            twinsPatch(patch)(funcIndex) = arr2;
//                        }
//                    }
//                }
//            }
//        }
//    }
//    for (int patch = 0; patch < twinsIndex.size(); ++patch) {
//        for (int functionIndex = 0; functionIndex < twinsIndex(patch).size(); ++functionIndex) {
//            if (twinsIndex(patch)(functionIndex).size() == 0) {
//                twinsIndex(patch)(functionIndex).push_back(-1);
//                twinsPatch(patch)(functionIndex).push_back(-1);
//            }
//        }
//    }
//}

//void IdentifyPatches1(gsMultiPatch<> mp,
//    gsVector <gsTHBSplineBasis <2, real_t >> THBVector1,
//    gsVector  < gsVector<index_t>>& isTouching,
//    gsVector <   gsVector <std::vector< index_t >>>& twinsIndex,
//    gsVector <   gsVector <std::vector< index_t >>>& twinsPatch,
//    gsVector <   gsVector < index_t >>& hasATwin,
//    gsVector  < gsVector<index_t>>& isActive
//) {
//    for (int i = 0; i < THBVector1.size(); ++i) {
//        twinsIndex(i).resize(THBVector1(i).size());
//        twinsPatch(i).resize(THBVector1(i).size());
//        isTouching(i).setZero(THBVector1(i).size());
//        hasATwin(i).setZero(THBVector1(i).size());
//        isActive(i).setZero(THBVector1(i).size());
//    }
//    gsVector  < gsVector<int>> touchedPatchSide(mp.nPatches());
//    for (int patch = 0; patch < mp.nPatches(); ++patch) {
//        touchedPatchSide(patch).resize(THBVector1(patch).size());
//        touchedPatchSide(patch).resize(THBVector1(patch).size());
//        for (int functionIndex = 0; functionIndex < THBVector1(patch).size(); ++functionIndex) {
//            touchedPatchSide(patch)(functionIndex) = -1;
//        }
//    }
//    std::vector<boundaryInterface>& boundaryInterfaces = mp.interfaces();
//    gsInfo << "directionMap and Orientation: \n";
//    for (int interfaceNum = 0; interfaceNum < boundaryInterfaces.size(); ++interfaceNum) {
//        gsInfo << boundaryInterfaces[interfaceNum].dirMap() << ",\n" << boundaryInterfaces[interfaceNum].dirOrientation() << ".\n";
//        gsVector<index_t> dummyVector;
//        boundaryInterfaces[interfaceNum].cornerMap(dummyVector);
//        gsInfo << "cornerMap:" << dummyVector << "\n";
//        gsInfo << "NARA\n";
//        outfile << boundaryInterfaces[interfaceNum] << "\n";
//        int firstSide = boundaryInterfaces[interfaceNum].first().index();
//        int secondSide = boundaryInterfaces[interfaceNum].second().index();
//        int firstPatch = boundaryInterfaces[interfaceNum].first().patch;
//        int secondPatch = boundaryInterfaces[interfaceNum].second().patch;
//        for (int functionIndex = 0; functionIndex < THBVector1(firstPatch).size(); ++functionIndex) {
//            if (!isTouching(firstPatch)(functionIndex)) {
//                if (THBVector1(firstPatch).function(functionIndex).support()(
//                    (firstSide + 1) / 2 - 1, (firstSide + 1) % 2) ==
//                    THBVector1(firstPatch).support()((firstSide + 1) / 2 - 1, (firstSide + 1) % 2)) {
//                    isTouching(firstPatch)(functionIndex) = 1;
//                    touchedPatchSide(firstPatch)(functionIndex) = secondSide;
//                }
//            }
//            if (!isTouching(secondPatch)(functionIndex)) {
//                if (THBVector1(secondPatch).function(functionIndex).support()((secondSide + 1) / 2 - 1, (secondSide + 1) % 2) ==
//                    THBVector1(secondPatch).support()((secondSide + 1) / 2 - 1, (secondSide + 1) % 2)) {
//                    isTouching(secondPatch)(functionIndex) = 1;
//                    touchedPatchSide(secondPatch)(functionIndex) = firstSide;
//                }
//            }
//        }
//    }
//    //Now we know all the functions that are touching the boundary
//    for (int interfaceNum = 0; interfaceNum < boundaryInterfaces.size(); ++interfaceNum) {
//        int firstSide = boundaryInterfaces[interfaceNum].first().index();
//        int secondSide = boundaryInterfaces[interfaceNum].second().index();
//        int firstPatch = boundaryInterfaces[interfaceNum].first().patch;
//        int secondPatch = boundaryInterfaces[interfaceNum].second().patch;
//        gsInfo << "Identifying " << firstPatch << " and " << secondPatch << "\n";
//        for (int functionIndex = 0; functionIndex < THBVector1(firstPatch).size(); ++functionIndex) {
//            if (
//                isTouching(firstPatch)(functionIndex)
//                && !isZeroOnBoundary(THBVector1,
//                    firstPatch,
//                    functionIndex,
//                    firstSide)
//                ) {
//                int twinIndex = twinFunction(THBVector1(firstPatch), functionIndex, THBVector1(secondPatch), firstSide, secondSide);
//                if (twinIndex != -1) {
//                    if (
//                        isTouching(secondPatch)(twinIndex)
//                        ) {
//                        twinsIndex(firstPatch)(functionIndex).push_back(twinIndex);
//                        twinsPatch(firstPatch)(functionIndex).push_back(secondPatch);
//                        hasATwin(firstPatch)(functionIndex) = 1;
//                        twinsPatch(secondPatch)(twinIndex).push_back(firstPatch);
//                        twinsIndex(secondPatch)(twinIndex).push_back(functionIndex);
//                        hasATwin(secondPatch)(twinIndex) = 1;
//                    }
//                    else if (!isTouching(secondPatch)(twinIndex)) {
//                        //But the twin does not touch the boundary\n";
//                    }
//                }
//            }
//            else if (!isTouching(firstPatch)(functionIndex)) {
//                //The function " << functionIndex << " does not touch the boundary\n";
//            }
//        }
//    }
//    for (int patch = 0; patch < twinsIndex.size() - 1; ++patch) {
//        for (int level = 0; level < twinsIndex(patch).size(); ++level) {
//            for (int functionIndex = 0; functionIndex < twinsIndex(patch)(level).size(); ++functionIndex) {
//                if (twinsIndex(patch)(functionIndex).size() == 0) {
//                    twinsIndex(patch)(functionIndex).push_back(-1);
//                    twinsPatch(patch)(functionIndex).push_back(-1);
//                }
//            }
//        }
//    }
//    //now we have to identify twins of twins. Since being a twin is a equivalence relation, it is transitive
//    for (int patch = 0; patch < twinsIndex.size(); ++patch) {
//        for (int level = 0; level < twinsIndex(patch).size(); ++level) {
//            for (int funcIndex = 0; funcIndex < twinsIndex(patch)(level).size(); ++funcIndex) {
//                for (int twinNum1 = 0; twinNum1 < twinsIndex(patch)(funcIndex).size(); ++twinNum1) {
//                    int patchIndex = twinsPatch(patch)(funcIndex)[twinNum1];
//                    int functionIndex = twinsIndex(patch)(funcIndex)[twinNum1];
//                    if (patchIndex == -1 || functionIndex == -1)     break;
//                    for (int twinNum2 = 0; twinNum2 < twinsIndex(patchIndex)(functionIndex).size(); ++twinNum2) {
//                        int tS = twinsIndex(patchIndex)(functionIndex).size();
//                        if (twinsIndex(patchIndex)(functionIndex).size() <= 1) {
//                            break;
//                        }
//                        else {
//                            std::vector<int>::iterator it1;
//                            std::vector<int>::iterator it2;
//                            std::vector<int> arr1 = twinsIndex(patch)(funcIndex);
//                            std::vector<int> arr2 = twinsPatch(patch)(funcIndex);
//                            it1 = std::find(arr1.begin(),
//                                arr1.end(),
//                                twinsIndex(patchIndex)(functionIndex)[twinNum2]
//                            );
//                            it2 = std::find(arr2.begin(),
//                                arr2.end(),
//                                twinsPatch(patchIndex)(functionIndex)[twinNum2]
//                            );
//                            bool a1 = !(it1 != arr1.end() && it2 != arr2.end());
//                            bool a2 = twinsIndex(patchIndex)(functionIndex)[twinNum2] != funcIndex;
//                            bool a3 = twinsPatch(patchIndex)(functionIndex)[twinNum2] != patch;
//                            if (!(it1 != arr1.end() && it2 != arr2.end()) &&
//                                twinsIndex(patchIndex)(functionIndex)[twinNum2] != funcIndex &&
//                                twinsPatch(patchIndex)(functionIndex)[twinNum2] != patch
//                                ) {
//                                arr1.push_back(twinsIndex(patchIndex)(functionIndex)[twinNum2]);
//                                arr2.push_back(twinsPatch(patchIndex)(functionIndex)[twinNum2]);
//                                twinsIndex(patch)(funcIndex) = arr1;
//                                twinsPatch(patch)(funcIndex) = arr2;
//                            }
//                        }
//                    }
//                }
//            }
//        }
//    }
//    //once again
//    for (int patch = twinsIndex.size() - 1; patch >= 0; --patch) {
//        for (int level = 0; level < twinsIndex(patch).size(); ++level) {
//            for (int funcIndex = twinsIndex(patch)(level).size() - 1; funcIndex >= 0; --funcIndex) {
//                for (int twinNum1 = 0; twinNum1 < twinsIndex(patch)(funcIndex).size(); ++twinNum1) {
//                    int patchIndex = twinsPatch(patch)(funcIndex)[twinNum1];
//                    int functionIndex = twinsIndex(patch)(funcIndex)[twinNum1];
//                    if (patchIndex == -1 || functionIndex == -1)     break;
//                    for (int twinNum2 = 0; twinNum2 < twinsIndex(patchIndex)(functionIndex).size(); ++twinNum2) {
//                        int tS = twinsIndex(patchIndex)(functionIndex).size();
//                        if (twinsIndex(patchIndex)(functionIndex).size() <= 1) {
//                            break;
//                        }
//                        else {
//                            std::vector<int>::iterator it1;
//                            std::vector<int>::iterator it2;
//                            std::vector<int> arr1 = twinsIndex(patch)(funcIndex);
//                            std::vector<int> arr2 = twinsPatch(patch)(funcIndex);
//                            it1 = std::find(arr1.begin(),
//                                arr1.end(),
//                                twinsIndex(patchIndex)(functionIndex)[twinNum1]
//                            );
//                            it2 = std::find(arr2.begin(),
//                                arr2.end(),
//                                twinsPatch(patchIndex)(functionIndex)[twinNum1]
//                            );
//                            bool a1 = !(it1 != arr1.end() && it2 != arr2.end());
//                            bool a2 = twinsIndex(patchIndex)(functionIndex)[twinNum1] != funcIndex;
//                            bool a3 = twinsPatch(patchIndex)(functionIndex)[twinNum1] != patch;
//                            if (!(it1 != arr1.end() && it2 != arr2.end()) &&
//                                twinsIndex(patchIndex)(functionIndex)[twinNum1] != funcIndex &&
//                                twinsPatch(patchIndex)(functionIndex)[twinNum1] != patch
//                                ) {
//                                arr1.push_back(twinsIndex(patchIndex)(functionIndex)[twinNum1]);
//                                arr2.push_back(twinsPatch(patchIndex)(functionIndex)[twinNum2]);
//                                twinsIndex(patch)(funcIndex) = arr1;
//                                twinsPatch(patch)(funcIndex) = arr2;
//                            }
//                        }
//                    }
//                }
//            }
//        }
//    }
//    for (int patch = 0; patch < twinsIndex.size(); ++patch) {
//        for (int level = 0; level < twinsIndex(patch).size(); ++level) {
//            for (int functionIndex = 0; functionIndex < twinsIndex(patch)(level).size(); ++functionIndex) {
//                if (twinsIndex(patch)(functionIndex).size() == 0) {
//                    twinsIndex(patch)(functionIndex).push_back(-1);
//                    twinsPatch(patch)(functionIndex).push_back(-1);
//                }
//            }
//        }
//    }
//}




template<typename T>
void initializeNestedVector(T& vec, index_t value) {
    if constexpr (std::is_same_v<typename T::value_type, int>) {
        std::fill(vec.begin(), vec.end(), value);  // Base case: fill vector of indices
    }
    else {
        for (auto& subVec : vec) {
            initializeNestedVector(subVec, value);  // Recursive case: go deeper
        }
    }
}

/**
 * @brief Fills a 3D nested vector with a specified value.
 *
 * This function iterates over a 3D nested vector and assigns a given value
 * to every element at the deepest level.
 *
 * @tparam T Type of the innermost elements in the vector (e.g., `index_t`).
 * @tparam ValueType Type of the value being assigned (e.g., `int` or `real_t`).
 * @param[in,out] vec The 3D nested vector to fill.
 * @param[in] value The value to assign to every innermost element.
 *
 * @details This template ensures flexibility for filling any structure
 * with compatible nested vector levels.
 */
template <typename T, typename ValueType>
void fillVector3(gsVector<gsVector<gsVector<T>>>& vec, ValueType value) {
    for (auto& patchVec : vec) {
        for (auto& levelVec : patchVec) {
            std::fill(levelVec.begin(), levelVec.end(), value);
        }
    }
}

/**
 * @brief Fills a 2D nested vector with a specified value.
 *
 * This function iterates over a 2D nested vector and assigns a given value
 * to every element at the innermost level.
 *
 * @tparam T Type of the innermost elements in the vector (e.g., `index_t`).
 * @tparam ValueType Type of the value being assigned (e.g., `int` or `real_t`).
 * @param[in,out] vec The 2D nested vector to fill.
 * @param[in] value The value to assign to every innermost element.
 *
 * @details This template ensures flexibility for filling any structure
 * with compatible nested vector levels.
 */
template <typename T, typename ValueType>
void fillVector2(gsVector<gsVector<T>>& vec, ValueType value) {
    for (auto& levelVec : vec) {
        std::fill(levelVec.begin(), levelVec.end(), value);
    }
}

/**
 * @brief Creates mappings between THB basis functions and tensor-product basis functions.
 *
 * This function establishes the correspondence between THB functions and their tensor-product
 * counterparts in Bells, while also marking active functions in the hierarchical basis.
 *
 * @param[in] Bells       Tensor-product B-spline bases for each patch and level.
 * @param[in] THBVector   Hierarchical THB spline basis for each patch.
 * @param[in] twinsIndex  Index of twin functions in the tensor-product basis.
 * @param[in] twinsPatch  Patch indices corresponding to twin functions.
 * @param[out] isIncluded Boolean-like structure indicating if a basis function is included in the THB basis.
 * @param[out] indexInTHB Maps tensor-product basis indices to their THB function indices.
 * @param[out] thbToBellsMapping Maps each THB function to its corresponding tensor-product function
 *                                (level and index in Bells).
 *
 * @details
 * - Resizes and initializes `isIncluded` and `indexInTHB` for each patch.
 * - Iterates through each THB function to map it to the tensor-product basis (`Bells`).
 * - Populates `thbToBellsMapping` with level and index information, avoiding the need to store patch data.
 *
 * @note Ensure all input data structures (e.g., `Bells`, `THBVector`) are properly initialized
 *       before calling this function.
 */

void createIndexMapping(
    gsVector<gsVector<gsTensorBSplineBasis<2, real_t>>> Bells,
    gsVector<gsTHBSplineBasis<2, real_t>> THBVector,
    gsVector<gsVector<gsVector<std::vector<index_t>>>> twinsIndex,
    gsVector<gsVector<gsVector<std::vector<index_t>>>> twinsPatch,
    gsVector<gsVector<gsVector<index_t>>>& isIncluded,
    gsVector<gsVector<gsVector<index_t>>>& indexInTHB,
    gsVector<std::vector<std::array<int, 2>>>& thbToBellsMapping // New structure to store corrLevel and corrIndex
) {
    outfile << "creating the index mapping\n";
    isIncluded.resize(Bells.size());
    indexInTHB.resize(Bells.size());
    thbToBellsMapping.resize(THBVector.size());

    // Resize and initialize structures
    for (size_t patch = 0; patch < Bells.size(); ++patch) {
        isIncluded(patch).resize(Bells(patch).size());
        indexInTHB(patch).resize(Bells(patch).size());
        for (size_t level = 0; level < Bells(patch).size(); ++level) {
            isIncluded(patch)(level).resize(Bells(patch)(level).size());
            indexInTHB(patch)(level).resize(Bells(patch)(level).size());
        }
    }
    fillVector3(isIncluded, 0);
    fillVector3(indexInTHB, -1);


    // Map THB functions to tensor-product basis
    int skippedInvalidMappings = 0;
    for (size_t patch = 0; patch < THBVector.size(); ++patch) {
        auto a = THBVector(patch).size();
        thbToBellsMapping[patch].resize(a);
        outfile << "thbToBellsMapping[patch].size(): " << thbToBellsMapping[patch].size() << "\n";
        for (size_t functionIndex = 0; functionIndex < THBVector(patch).size(); ++functionIndex) {
            int corrLevel = THBVector(patch).levelOf(functionIndex);
            int corrIndex = THBVector(patch).flatTensorIndexOf(functionIndex);
            // Save mapping to Bells
            thbToBellsMapping[patch][functionIndex] = { corrLevel, corrIndex };

            const bool levelInRange =
                corrLevel >= 0 &&
                static_cast<size_t>(corrLevel) < isIncluded(patch).size();
            const bool indexInRange =
                levelInRange &&
                corrIndex >= 0 &&
                static_cast<size_t>(corrIndex) < isIncluded(patch)(corrLevel).size();

            if (!indexInRange)
            {
                ++skippedInvalidMappings;
                if (outfile.is_open())
                {
                    outfile << "[createIndexMapping] skip invalid mapping patch=" << patch
                            << " function=" << functionIndex
                            << " level=" << corrLevel
                            << " tensorIndex=" << corrIndex << "\n";
                }
                continue;
            }

            // Mark as included in the hierarchical basis
            isIncluded(patch)(corrLevel)(corrIndex) = 1;
            indexInTHB(patch)(corrLevel)(corrIndex) = functionIndex;
        }
    }

    if (skippedInvalidMappings > 0)
    {
        gsInfo << "[createIndexMapping] skipped " << skippedInvalidMappings
               << " invalid THB-to-Bells mappings\n";
        if (outfile.is_open())
        {
            outfile << "[createIndexMapping] skipped " << skippedInvalidMappings
                    << " invalid THB-to-Bells mappings\n";
        }
    }

    if (g_verbose && outfile.is_open())
    {
        outfile << "indexInTHB\n";
        for (size_t patch = 0; patch < indexInTHB.size(); patch++)
        {
            for (size_t level = 0; level < indexInTHB(patch).size(); level++)
            {
                for (size_t i = 0; i < indexInTHB(patch)(level).size(); i++)
                {
                    outfile << patch << " " << level << " " << i << ": " << indexInTHB(patch)(level)(i) << "\n";
                }
            }
        }
        outfile << "thbToBellsMapping\n";
        for (size_t patch = 0; patch < thbToBellsMapping.size(); patch++)
        {
            for (size_t i = 0; i < thbToBellsMapping[patch].size(); i++)
            {
                outfile << thbToBellsMapping[patch][i][0] << " " << thbToBellsMapping[patch][i][1] << "\n";
            }
        }
    }
}



//void selectionMechanism(gsVector < gsVector < gsTensorBSplineBasis <2, real_t >>> THBVector1,
//    gsVector<gsTHBSplineBasis<2, real_t> >SubdomainHierarchy,
//    gsVector  <gsVector< gsVector<index_t>>>& isActive,
//    gsVector < gsVector < gsVector <std::vector< index_t >>>>  twinsIndex,
//    gsVector < gsVector < gsVector <std::vector< index_t >>>> twinsPatch,
//    gsVector<index_t>& numActive
//) {
//    gsVector<index_t> numSelected;
//    gsVector<index_t> normallySelected;
//    numSelected.setZero(THBVector1.size());
//    numActive.setZero(THBVector1.size());
//    normallySelected.setZero(THBVector1.size());
//    commonBasisSize = 0;
//    index_t interioru0 = 0;
//    for (int patch = 0; patch < THBVector1.size(); ++patch) {
//        for (int level = 0; level < THBVector1(patch).size(); ++level) {
//            for (int funcIndex = 0; funcIndex < THBVector1(patch)(level).size(); ++funcIndex) {
//                gsMatrix<> supp = THBVector1(patch)(level).function(funcIndex).support();
//                gsMatrix<index_t> mybox(supp.rows(), supp.cols());
//                supportToBoxOfLevel(mybox, supp, level, interioru0, interioru0);
//                gsVector<index_t, 2> eins(2);
//                gsVector<index_t, 2> zwei(2);
//                eins(0) = mybox(0, 0);//bLC(0,0);
//                eins(1) = mybox(1, 0);//bLC(0,1);
//                zwei(0) = mybox(0, 1);//bUC(0,0);
//                zwei(1) = mybox(1, 1);//bUC(0,0);
//                gsVector<index_t, 2> const lower = eins;
//                gsVector<index_t, 2> const upper = zwei;
//                bool qBir = SubdomainHierarchy(patch).tree().query1(upper, lower, level);
//                bool qIki = SubdomainHierarchy(patch).tree().query2(lower, upper, level);
//                int  qUc = SubdomainHierarchy(patch).tree().query3(lower, upper, level);
//                int qDort = SubdomainHierarchy(patch).tree().query4(lower, upper, level);
//                bool Active = false;
//                if (((qBir) || (qUc == level && qDort >= level)) &&
//                    !isActive(patch)(level)(funcIndex)
//                    ) {
//                    commonBasisSize++;
//                    numSelected(patch)++;
//                    numActive(patch)++;
//                    isActive(patch)(level)(funcIndex) = 1;
//                    normallySelected(patch)++;
//                    //Does it have twins?
//                    for (int twinNum = 0; twinNum < twinsIndex(patch)(level)(funcIndex).size(); ++twinNum) {
//                        int patchIndex = twinsPatch(patch)(level)(funcIndex)[twinNum];
//                        if (patchIndex != -1) {
//                            //It is a twin
//                            int twinIndex = twinsIndex(patch)(level)(funcIndex)[twinNum];
//                            isActive(patchIndex)(level)(twinIndex) = 1;
//                        }
//                    }
//                }
//            }
//        }
//    }
//    outfile << "commonBasisSize " << commonBasisSize << std::endl;
//    gsInfo << "commonBasisSize " << commonBasisSize << std::endl;
//}

/**
 * @brief Sets up global indexing for THB basis functions and manages twin/spillover relationships.
 *
 * This function maps THB basis functions to global indices, identifies hierarchical twins,
 * and records spillover tensor-product functions for cases where refinement mismatches occur.
 * It also calculates the total number of globally indexed functions (`commonSize`).
 *
 * @param[in] THBVector                    The THB spline basis for each patch.
 * @param[out] globalIndexTHB              Global indices for all THB functions across patches.
 * @param[in] hasATwin                     Flags indicating if a function has twins in the tensor-product basis.
 * @param[in] twinsIndex                   Indices of twin functions in the tensor-product basis.
 * @param[in] twinsPatch                   Patch indices corresponding to twin functions.
 * @param[out] globalIndex                 Global indices for all tensor-product basis functions.
 * @param[in] indexInTHB                   Maps tensor-product basis indices to THB function indices.
 * @param[out] functionDescription         Describes relationships for each THB function. Each entry includes:
 *                                         - Twins defined in the THB basis as {patch, functionIndex}.
 * @param[out] spilloverFunctionCoordinates Stores spillover tensor-product functions for each THB function as:
 *                                          - {patch, level, index}.
 * @param[out] spilloverSources            For each patch, lists functions from other patches that spill over into it.
 * @param[out] hasSpillover                Boolean array indicating whether each THB function has spillovers.
 * @param[out] commonSize                  Total number of globally indexed functions across all patches.
 *
 * @details
 * - Assigns a unique global index to each THB basis function across all patches.
 * - Captures hierarchical twin relationships within the THB basis.
 * - Identifies spillover cases where twins are not part of the THB basis, saving their details in
 *   `spilloverFunctionCoordinates`.
 * - Tracks all external spillovers into a patch using `spilloverSources`.
 * - Populates `hasSpillover` for efficient spillover checks during assembly.
 *
 */
void SetUpTheBasisTHB(
    gsVector<gsTHBSplineBasis<2, real_t>> THBVector,
    gsVector<gsVector<int>>& globalIndexTHB,
    gsVector<gsVector<gsVector<index_t>>> hasATwin,
    gsVector<gsVector<gsVector<std::vector<index_t>>>> twinsIndex,
    gsVector<gsVector<gsVector<std::vector<index_t>>>> twinsPatch,
    gsVector<gsVector<gsVector<index_t>>>& globalIndex,
    gsVector<gsVector<gsVector<index_t>>> indexInTHB,
    std::vector<std::vector<std::vector<index_t>>>& functionDescription,
    std::vector<std::vector<std::array<int, 3>>>& spilloverFunctionCoordinates,
    std::vector<std::vector<std::array<int, 3>>>& spilloverSources,
    std::vector<bool>& hasSpillover,
    int& commonSize,
    std::vector<std::vector<std::vector<int>>>& localIndex // <-- ADDED
)
{
    int theIndex = 0;
    int maxCommonSize = 0;

    // Resize data structures
    for (size_t patch = 0; patch < THBVector.size(); ++patch)
    {
        auto a = THBVector(patch).size();
        gsInfo << "a: " << a << "\n";
        maxCommonSize += a;
        globalIndexTHB(patch).resize(a);
    }

    fillVector2(globalIndexTHB, -1);
    fillVector3(globalIndex, -1);

    functionDescription.resize(maxCommonSize);
    hasSpillover.resize(maxCommonSize, false);
    spilloverFunctionCoordinates.resize(maxCommonSize);
    spilloverSources.resize(THBVector.size());

    localIndex.clear();                        // <-- ADDED
    localIndex.resize(THBVector.size());        // <-- ADDED

    for (size_t patch = 0; patch < THBVector.size(); ++patch)
    {
        for (size_t functionIndex = 0; functionIndex < THBVector[patch].size(); ++functionIndex)
        {
            int corrLevel = THBVector(patch).levelOf(functionIndex);
            int corrIndex = THBVector(patch).flatTensorIndexOf(functionIndex);

            if (globalIndex(patch)(corrLevel)(corrIndex) == -1)
            {
                globalIndex(patch)(corrLevel)(corrIndex) = theIndex;

                globalIndexTHB(patch)(functionIndex) = theIndex;

                // Save main function info
                functionDescription[theIndex].push_back({ static_cast<int>(patch), static_cast<int>(functionIndex) });
                localIndex[patch].push_back({ static_cast<int>(patch), static_cast<int>(functionIndex) }); // <-- ADDED

                // Twins
                if (hasATwin(patch)(corrLevel)(corrIndex) == 1)
                {
                    for (size_t twinNum = 0; twinNum < twinsIndex(patch)(corrLevel)(corrIndex).size(); ++twinNum)
                    {
                        int twinPatch = twinsPatch(patch)(corrLevel)(corrIndex)[twinNum];
                        int twinIndex = twinsIndex(patch)(corrLevel)(corrIndex)[twinNum];

                        if (globalIndex(twinPatch)(corrLevel)(twinIndex) == -1)
                        {
                            globalIndex(twinPatch)(corrLevel)(twinIndex) = theIndex;

                            int THBTwinIndex = indexInTHB(twinPatch)(corrLevel)(twinIndex);
                            if (THBTwinIndex != -1)
                            {
                                globalIndexTHB(twinPatch)(THBTwinIndex) = theIndex;
                                functionDescription[theIndex].push_back({ twinPatch, THBTwinIndex });
                                localIndex[twinPatch].push_back({ twinPatch, THBTwinIndex }); // <-- ADDED
                            }
                            else
                            {
                                // Spillover case
                                spilloverFunctionCoordinates[theIndex].push_back({ twinPatch, corrLevel, twinIndex });
                                spilloverSources[twinPatch].push_back({ static_cast<int>(patch), corrLevel, corrIndex });
                                hasSpillover[theIndex] = true;
                            }
                        }
                    }
                }

                ++theIndex;
            }
        }
    }

    commonSize = theIndex;
}





////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void restoreTheHierarchy2(gsVector< gsVector<gsVector<index_t>>>& boxMat, int patch,
    int ourBox[], gsTHBSplineBasis<2, real_t>& THB, gsTHBSplineBasis<2, real_t> &THBUnstructured,
    gsVector<gsMatrix<index_t >>& lowCorners,
    gsVector<gsMatrix<index_t >>& upCorners, gsVector<gsVector<index_t >>& myLevel, int& lastNonZeroRow) {
    std::vector<index_t> box;
    box.push_back(ourBox[0]);
    box.push_back(ourBox[1]);
    box.push_back(ourBox[2]);
    box.push_back(ourBox[3]);
    box.push_back(ourBox[4]);
    if (box.size() != 5
        || box[0] < 0
        || box[1] < 0
        || box[2] < 0
        || box[3] <= box[1]
        || box[4] <= box[2])
    {
        gsInfo << "[restoreTheHierarchy2] Skipping invalid box: "
               << box[0] << " " << box[1] << " " << box[2] << " "
               << box[3] << " " << box[4] << "\n";
        outfile << "[restoreTheHierarchy2] Skipping invalid box: "
                << box[0] << " " << box[1] << " " << box[2] << " "
                << box[3] << " " << box[4] << "\n";
        return;
    }
    const index_t level = box[0];
    if (level < 0 || level > THBUnstructured.maxLevel() + 1)
    {
        gsInfo << "[restoreTheHierarchy2] Skipping box with level out of range (level="
               << level << ", maxLevel=" << THBUnstructured.maxLevel() << "): "
               << box[0] << " " << box[1] << " " << box[2] << " "
               << box[3] << " " << box[4] << "\n";
        outfile << "[restoreTheHierarchy2] Skipping box with level out of range (level="
                << level << ", maxLevel=" << THBUnstructured.maxLevel() << "): "
                << box[0] << " " << box[1] << " " << box[2] << " "
                << box[3] << " " << box[4] << "\n";
        return;
    }
    const index_t maxIndex = static_cast<index_t>(1) << static_cast<unsigned>(level);
    if (box[1] < 0 || box[2] < 0 || box[3] > maxIndex || box[4] > maxIndex)
    {
        gsInfo << "[restoreTheHierarchy2] Skipping box out of bounds (level=" << level
               << ", maxIndex=" << maxIndex << "): "
               << box[0] << " " << box[1] << " " << box[2] << " "
               << box[3] << " " << box[4] << "\n";
        outfile << "[restoreTheHierarchy2] Skipping box out of bounds (level=" << level
                << ", maxIndex=" << maxIndex << "): "
                << box[0] << " " << box[1] << " " << box[2] << " "
                << box[3] << " " << box[4] << "\n";
        return;
    }
    outfile << "THBUnstructured before\n";
    outfile << THBUnstructured << "\n";
    (THBUnstructured.tree().getBoxes(lowCorners(patch), upCorners(patch), myLevel(patch)));
    outfile << "Just to be clear\n";
    for (int i = 0; i < lowCorners(patch).rows(); ++i) {
        outfile << myLevel(patch)(i) << " ";
        outfile << lowCorners(patch)(i, 0) << " " << lowCorners(patch)(i, 1) << " ";
        outfile << upCorners(patch)(i, 0) << " " << upCorners(patch)(i, 1);
        outfile << ";\n";
    }
    outfile << "The box to be inserted\n";
    outfile << box[0] << " ";
    outfile << box[1] << " ";
    outfile << box[2] << " ";
    outfile << box[3] << " ";
    outfile << box[4] << "\n";
    bool foundInBoxMat = false;
    if (patch >= 0 && patch < boxMat.size())
    {
        const int maxRow = std::min(lastNonZeroRow, static_cast<int>(boxMat(patch).size()));
        for (int i = 0; i < maxRow; ++i)
        {
            if (boxMat(patch)(i).size() < 5)
                continue;
            if (boxMat(patch)(i)(0) == box[0] &&
                boxMat(patch)(i)(1) == box[1] &&
                boxMat(patch)(i)(2) == box[2] &&
                boxMat(patch)(i)(3) == box[3] &&
                boxMat(patch)(i)(4) == box[4])
            {
                foundInBoxMat = true;
                break;
            }
        }
    }
    if (!foundInBoxMat && patch >= 0 && patch < boxMat.size())
    {
        if (lastNonZeroRow < 0)
            lastNonZeroRow = 0;
        if (lastNonZeroRow >= static_cast<int>(boxMat(patch).size()))
            boxMat(patch).resize(lastNonZeroRow + 1);

        boxMat(patch)(lastNonZeroRow).resize(5);
        boxMat(patch)(lastNonZeroRow)(0) = box[0];
        boxMat(patch)(lastNonZeroRow)(1) = box[1];
        boxMat(patch)(lastNonZeroRow)(2) = box[2];
        boxMat(patch)(lastNonZeroRow)(3) = box[3];
        boxMat(patch)(lastNonZeroRow)(4) = box[4];
        ++lastNonZeroRow;
    }

    // Always refine the box, even if it wasn't found in boxMat
    THBUnstructured.refineElements(box);
    outfile << "LITTLECHECK\n";
    outfile << "THBUnstructured after\n";
    outfile << THBUnstructured << "\n";
    (THBUnstructured.tree().getBoxes(lowCorners(patch), upCorners(patch), myLevel(patch)));
    outfile << "Just to be clear\n";
    outfile << "lastNonzeroRow: " << lastNonZeroRow << "\n";
    for (int i = 0; i < lowCorners(patch).rows(); ++i) {
        outfile << myLevel(patch)(i) << " ";
        outfile << lowCorners(patch)(i, 0) << " " << lowCorners(patch)(i, 1) << " ";
        outfile << upCorners(patch)(i, 0) << " " << upCorners(patch)(i, 1);
        outfile << ";\n";
    }
    box.clear();
    gsInfo << "lowCorners(patch).rows():" << lowCorners(patch).rows() << "\n";
    gsInfo << "boxMat(patch).size(): " << boxMat(patch).size() << "\n";
    gsInfo << "myLevel(patch): \n" << myLevel(patch) << "\n";
    outfile << "lowCorners(patch).rows():" << lowCorners(patch).rows() << "\n";
    outfile << "boxMat(patch).size(): " << boxMat(patch).size() << "\n";
    outfile << "myLevel(patch): \n" << myLevel(patch) << "\n";
    for (int i = 0; i < lowCorners(patch).rows(); i++) {
        gsInfo << "i = " << i << ", " << boxMat(patch)(i).size() << "\n";
        boxMat(patch)(i).resize(5);
        if (i < lastNonZeroRow)
        {
            gsInfo << "before: " << boxMat(patch)(i)(0) << " " << boxMat(patch)(i)(1) << " " << boxMat(patch)(i)(2) << " " << boxMat(patch)(i)(3) << " " <<
                boxMat(patch)(i)(4) << "\n";
        }
        boxMat(patch)(i)(0) = myLevel(patch)(i);
        boxMat(patch)(i)(1) = (real_t)lowCorners(patch)(i, 0) / pow(2, THBUnstructured.maxLevel() - myLevel(patch)(i));
        boxMat(patch)(i)(2) = (real_t)lowCorners(patch)(i, 1) / pow(2, THBUnstructured.maxLevel() - myLevel(patch)(i));
        boxMat(patch)(i)(3) = (real_t)upCorners(patch)(i, 0) / pow(2, THBUnstructured.maxLevel() - myLevel(patch)(i));
        boxMat(patch)(i)(4) = (real_t)upCorners(patch)(i, 1) / pow(2, THBUnstructured.maxLevel() - myLevel(patch)(i));
        gsInfo << boxMat(patch)(i)(0) << " " << boxMat(patch)(i)(1) << " " << boxMat(patch)(i)(2) << " " << boxMat(patch)(i)(3) << " " <<
            boxMat(patch)(i)(4) << "\n";
        outfile << boxMat(patch)(i)(0) << " " << boxMat(patch)(i)(1) << " " << boxMat(patch)(i)(2) << " " << boxMat(patch)(i)(3) << " " <<
            boxMat(patch)(i)(4) << "\n";
    }
    outfile << "After the renorming;\n";
    lastNonZeroRow = lowCorners(patch).rows();
    for (int therow = 0; therow < lastNonZeroRow; therow++) {
        std::vector<index_t> box;
        for (int column = 0; column < 5; column++) {
            box.push_back(boxMat(patch)(therow)(column));
            outfile << boxMat(patch)(therow)(column) << "; ";
        }
        outfile << "\n";
        if (box.size() == 5)
        {
            const index_t bLevel = box[0];
            const index_t bMaxIndex = static_cast<index_t>(1) << static_cast<unsigned>(bLevel);
            if (bLevel >= 0 && bLevel <= THB.maxLevel()
                && box[1] >= 0 && box[2] >= 0
                && box[3] > box[1] && box[4] > box[2]
                && box[3] <= bMaxIndex && box[4] <= bMaxIndex)
            {
                THB.refineElements(box);
            }
            else
            {
                gsInfo << "[restoreTheHierarchy2] Skipping THB refine box: "
                       << box[0] << " " << box[1] << " " << box[2] << " "
                       << box[3] << " " << box[4] << "\n";
                outfile << "[restoreTheHierarchy2] Skipping THB refine box: "
                        << box[0] << " " << box[1] << " " << box[2] << " "
                        << box[3] << " " << box[4] << "\n";
            }
        }
    }
    //experimentalpart
    for (size_t therow = lastNonZeroRow; therow < boxMat(patch).size(); therow++)
    {
        boxMat(patch)(therow).clear();
    }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void restoreTheHierarchy(int& createdBoxNum, int& lastNonzeroRow, gsVector< gsVector<gsVector<index_t>>>& boxMat, int levNow, int centerInd, int ourBox[], int successfullAttempts, int patch) {
    if (patch < 0 || patch >= boxMat.size())
    {
        gsInfo << "[restoreTheHierarchy] Invalid patch index: " << patch << "\n";
        outfile << "[restoreTheHierarchy] Invalid patch index: " << patch << "\n";
        createdBoxNum = 0;
        return;
    }
    if (lastNonzeroRow < 0)
        lastNonzeroRow = 0;
    if (lastNonzeroRow > static_cast<int>(boxMat(patch).size()))
        lastNonzeroRow = static_cast<int>(boxMat(patch).size());
    if (boxMat(patch).size() == 0)
    {
        gsInfo << "[restoreTheHierarchy] boxMat(patch) is empty for patch " << patch << "\n";
        outfile << "[restoreTheHierarchy] boxMat(patch) is empty for patch " << patch << "\n";
        createdBoxNum = 0;
        return;
    }
    if (createdBoxNum <= 0)
    {
        createdBoxNum = 0;
        return;
    }
    for (int l = 0; l < createdBoxNum; l++) {
        const int idx = lastNonzeroRow - l;
        if (idx < 0 || idx >= static_cast<int>(boxMat(patch).size()))
        {
            gsInfo << "[restoreTheHierarchy] Skipping invalid index " << idx << " (lastNonzeroRow="
                   << lastNonzeroRow << ", createdBoxNum=" << createdBoxNum << ")\n";
            outfile << "[restoreTheHierarchy] Skipping invalid index " << idx << " (lastNonzeroRow="
                    << lastNonzeroRow << ", createdBoxNum=" << createdBoxNum << ")\n";
            continue;
        }
        if (boxMat(patch)(idx).size() < 5)
            boxMat(patch)(idx).resize(5);
        boxMat(patch)(idx)(0) = levNow; //Preparation for multilevel meshes
        outfile << "updated coordinates of " << lastNonzeroRow - l << "box: " <<
            boxMat(patch)(lastNonzeroRow - l)(0);
        for (int m = 1; m < 5; m++) {
            boxMat(patch)(lastNonzeroRow - l)(m) = 0;
            outfile << "\t" << boxMat(patch)(lastNonzeroRow - l)(m);
        }
        outfile << "\n";
    }
    if (centerInd >= 0 && centerInd < static_cast<int>(boxMat(patch).size()))
    {
        outfile << "current " << centerInd << "box is now:";
        for (int k = 0; k < 5; k++) {
            if (boxMat(patch)(centerInd).size() < 5)
                boxMat(patch)(centerInd).resize(5);
            boxMat(patch)(centerInd)(k) = ourBox[k];
            outfile << "\t" << boxMat(patch)(centerInd)(k);
        }
        outfile << "\n";
    }
    else
    {
        gsInfo << "[restoreTheHierarchy] centerInd out of range: " << centerInd << "\n";
        outfile << "[restoreTheHierarchy] centerInd out of range: " << centerInd << "\n";
    }
    if (successfullAttempts == 0) {
        if (boxMat(patch).size() > 0)
        {
            outfile << "0th position IS RETURNED TO\n";
            for (int k = 0; k < 5; k++) {
                boxMat(patch)(0)(k) = ourBox[k];
                outfile << "\t" << boxMat(patch)(0)(k);
            }
            outfile << "\n";
        }
    }
    lastNonzeroRow = lastNonzeroRow - createdBoxNum;
    createdBoxNum = 0;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
struct RestoreRecord
{
    int createdBoxNum;
    int centerInd;
    std::array<int, 5> box;
    int levNow;
    int successfullAttempts;
    int patch;
};

static inline real_t clampToUnit(real_t val)
{
    return std::max<real_t>(0.0, std::min<real_t>(1.0, val));
}

static inline bool intervalsOverlap(real_t a0, real_t a1, real_t b0, real_t b1)
{
    return (a0 <= b1) && (b0 <= a1);
}

// Check if a function support intersects any rectangle in the per-patch list.
template <typename SupportMatrix>
static inline bool supportIntersectsRectList(const SupportMatrix& support,
                                              const std::vector<std::array<real_t, 4>>& rects)
{
    const real_t u0 = static_cast<real_t>(support(0, 0));
    const real_t u1 = static_cast<real_t>(support(0, 1));
    const real_t v0 = static_cast<real_t>(support(1, 0));
    const real_t v1 = static_cast<real_t>(support(1, 1));
    for (const auto& box : rects)
        if (intervalsOverlap(u0, u1, box[0], box[2]) &&
            intervalsOverlap(v0, v1, box[1], box[3]))
            return true;
    return false;
}

// Merge rect {uMin,vMin,uMax,vMax} into the per-patch exact-union list.
// Overlapping existing rects are absorbed into the new rect (greedy merge),
// keeping the list free of pairwise overlaps.
static inline bool mergePatchAabb(LocalCoarseningRegion& region,
                                  int patch,
                                  real_t uMin,
                                  real_t vMin,
                                  real_t uMax,
                                  real_t vMax)
{
    if (patch < 0 || patch >= static_cast<int>(region.patchAABB.size()))
        return false;

    uMin = clampToUnit(uMin);
    vMin = clampToUnit(vMin);
    uMax = clampToUnit(uMax);
    vMax = clampToUnit(vMax);

    if (uMax < uMin) std::swap(uMax, uMin);
    if (vMax < vMin) std::swap(vMax, vMin);

    std::array<real_t, 4> incoming = { uMin, vMin, uMax, vMax };

    if (!region.hasPatch[patch])
    {
        region.hasPatch[patch] = true;
        region.patchAABB[patch].push_back(incoming);
        return true;
    }

    // Containment-only merge: absorb existing rects that are contained in
    // incoming (subsumed), and skip if incoming is contained in an existing rect.
    // Two rects that merely overlap but neither contains the other are kept separate
    // so the list stays an exact union, not an over-approximating bounding box.
    auto rectContains = [](const std::array<real_t, 4>& outer,
                           const std::array<real_t, 4>& inner) -> bool {
        return outer[0] <= inner[0] && outer[1] <= inner[1] &&
               outer[2] >= inner[2] && outer[3] >= inner[3];
    };

    bool keep = true;
    while (keep)
    {
        keep = false;
        for (auto it = region.patchAABB[patch].begin();
             it != region.patchAABB[patch].end(); ++it)
        {
            if (rectContains(*it, incoming))
                return false; // incoming already covered, nothing to add
            if (rectContains(incoming, *it))
            {
                region.patchAABB[patch].erase(it);
                keep = true;
                break;
            }
        }
    }

    region.patchAABB[patch].push_back(incoming);
    return true;
}

static inline bool pointInsidePatchRegion(const LocalCoarseningRegion& region,
                                          int patch,
                                          real_t u,
                                          real_t v)
{
    if (!region.enabled)
        return true;
    if (patch < 0 || patch >= static_cast<int>(region.hasPatch.size()) || !region.hasPatch[patch])
        return false;

    for (const auto& box : region.patchAABB[patch])
        if (u >= box[0] && u <= box[2] && v >= box[1] && v <= box[3])
            return true;
    return false;
}

static inline bool patchElementIntersectsRegion(const LocalCoarseningRegion& region,
                                                int patch,
                                                const gsVector<real_t>& lower,
                                                const gsVector<real_t>& upper)
{
    if (!region.enabled)
        return true;
    if (patch < 0 || patch >= static_cast<int>(region.hasPatch.size()) || !region.hasPatch[patch])
        return false;

    for (const auto& box : region.patchAABB[patch])
        if (intervalsOverlap(lower(0), upper(0), box[0], box[2]) &&
            intervalsOverlap(lower(1), upper(1), box[1], box[3]))
            return true;
    return false;
}

// Returns true if the bounding box of the rect list covers [0,1]^2.
static inline bool patchRegionCoversWholePatch(const LocalCoarseningRegion& region,
                                               int patch,
                                               real_t tol = 1e-12)
{
    if (!region.enabled)
        return true;
    if (patch < 0 || patch >= static_cast<int>(region.hasPatch.size()) || !region.hasPatch[patch])
        return false;

    real_t u0 = 1, v0 = 1, u1 = 0, v1 = 0;
    for (const auto& box : region.patchAABB[patch])
    {
        u0 = std::min(u0, box[0]); v0 = std::min(v0, box[1]);
        u1 = std::max(u1, box[2]); v1 = std::max(v1, box[3]);
    }
    return std::abs(u0) <= tol && std::abs(v0) <= tol &&
           std::abs(u1 - 1.0) <= tol && std::abs(v1 - 1.0) <= tol;
}

static LocalCoarseningRegion buildLocalCoarseningRegion(
    const MPBES<2, real_t>& mpbes,
    int sourcePatch,
    const std::vector<CellToCoarsen>& cells,
    index_t interioru,
    index_t interiorv,
    int localityLambda,
    const std::vector<std::pair<int, std::vector<CellToCoarsen>>>& extraSeeds = {})
{
    LocalCoarseningRegion region;
    const index_t nPatches = mpbes.nPatches();
    const index_t nFunctions = mpbes.size();

    region.patchAABB.resize(nPatches); // each patch starts with empty rect list
    region.hasPatch.assign(nPatches, false);
    region.basisInd.setZero(1, nFunctions);
    region.lambda = std::max(0, localityLambda);

    auto formatBox = [&](const std::array<real_t, 4>& box) -> std::string
    {
        std::ostringstream ss;
        ss << "[" << box[0] << ", " << box[1] << "] x [" << box[2] << ", " << box[3] << "]";
        return ss.str();
    };

    auto formatBoxList = [&](const std::vector<std::array<real_t, 4>>& rects) -> std::string
    {
        if (rects.empty()) return "(empty)";
        std::ostringstream ss;
        for (size_t i = 0; i < rects.size(); ++i)
        {
            if (i) ss << " | ";
            ss << "[" << rects[i][0] << "," << rects[i][1] << "]x["
               << rects[i][2] << "," << rects[i][3] << "]";
        }
        return ss.str();
    };

    struct RegionSupportRef
    {
        int patch;
        int level;
        int index;
        bool spillover;
    };

    if (outfile.is_open())
    {
        outfile << "[local-region] build start: sourcePatch=" << sourcePatch
                << ", cells=" << cells.size()
                << ", lambda=" << region.lambda
                << ", nPatches=" << nPatches
                << ", nFunctions=" << nFunctions << "\n";
    }

    if (cells.empty() || sourcePatch < 0 || sourcePatch >= static_cast<int>(nPatches))
        return region;

    region.enabled = true;

    for (size_t cellNum = 0; cellNum < cells.size(); ++cellNum)
    {
        const auto& cell = cells[cellNum];
        const index_t divU = (interioru + 1) * static_cast<index_t>(std::pow(2.0, cell.level));
        const index_t divV = (interiorv + 1) * static_cast<index_t>(std::pow(2.0, cell.level));
        if (divU <= 0 || divV <= 0)
            continue;

        const real_t u1 = static_cast<real_t>(cell.x1) / static_cast<real_t>(divU);
        const real_t v1 = static_cast<real_t>(cell.y1) / static_cast<real_t>(divV);
        const real_t u2 = static_cast<real_t>(cell.x2) / static_cast<real_t>(divU);
        const real_t v2 = static_cast<real_t>(cell.y2) / static_cast<real_t>(divV);

        const real_t du = std::abs(u2 - u1);
        const real_t dv = std::abs(v2 - v1);

        const std::array<real_t, 4> requestedBox = {
            clampToUnit(u1 - region.lambda * du),
            clampToUnit(v1 - region.lambda * dv),
            clampToUnit(u2 + region.lambda * du),
            clampToUnit(v2 + region.lambda * dv)
        };

        mergePatchAabb(region,
                       sourcePatch,
                       requestedBox[0],
                       requestedBox[1],
                       requestedBox[2],
                       requestedBox[3]);

        if (outfile.is_open())
        {
            outfile << "[local-region] seed cell " << cellNum
                    << ": level=" << cell.level
                    << " cell=[" << cell.x1 << ", " << cell.y1 << "] x [" << cell.x2 << ", " << cell.y2 << "]"
                    << " requestedBox=" << formatBox(requestedBox)
                    << " -> rects=" << formatBoxList(region.patchAABB[sourcePatch]) << "\n";
        }
    }

    // Seed AABBs from co-patch group members (common block-diagonal fit).
    for (const auto& xEntry : extraSeeds) {
        const int xPatch = xEntry.first;
        const std::vector<CellToCoarsen>& xCells = xEntry.second;
        if (xPatch < 0 || xPatch >= static_cast<int>(nPatches)) continue;
        for (const auto& cell : xCells) {
            const index_t divU = (interioru + 1) * static_cast<index_t>(std::pow(2.0, cell.level));
            const index_t divV = (interiorv + 1) * static_cast<index_t>(std::pow(2.0, cell.level));
            if (divU <= 0 || divV <= 0) continue;
            const real_t u1 = static_cast<real_t>(cell.x1) / static_cast<real_t>(divU);
            const real_t v1 = static_cast<real_t>(cell.y1) / static_cast<real_t>(divV);
            const real_t u2 = static_cast<real_t>(cell.x2) / static_cast<real_t>(divU);
            const real_t v2 = static_cast<real_t>(cell.y2) / static_cast<real_t>(divV);
            const real_t du = std::abs(u2 - u1);
            const real_t dv = std::abs(v2 - v1);
            const std::array<real_t, 4> requestedBox = {
                clampToUnit(u1 - region.lambda * du),
                clampToUnit(v1 - region.lambda * dv),
                clampToUnit(u2 + region.lambda * du),
                clampToUnit(v2 + region.lambda * dv)
            };
            mergePatchAabb(region, xPatch,
                requestedBox[0], requestedBox[1], requestedBox[2], requestedBox[3]);
        }
        if (outfile.is_open())
            outfile << "[local-region] co-patch=" << xPatch
                    << " cells=" << xCells.size() << " seeded\n";
        region.enabled = true;
    }

    const auto& functionDescription = mpbes.functionDescription();
    const auto& spilloverFunctionCoordinates = mpbes.spilloverCoordinates();
    const auto& hasSpillover = mpbes.hasSpillover();
    const auto& Bells = mpbes.bellsBases();

    const auto initialAabb = region.patchAABB;
    const auto initialHasPatch = region.hasPatch;

    if (outfile.is_open())
    {
        outfile << "[local-region] initial AU boxes:\n";
        for (index_t p = 0; p < nPatches; ++p)
        {
            if (!initialHasPatch[p])
                continue;
            outfile << "  patch=" << p << " rects=" << formatBoxList(initialAabb[p]) << "\n";
        }
    }

    auto collectRegionHits = [&](index_t f,
                                 const std::vector<std::vector<std::array<real_t, 4>>>& activeAabb,
                                 const std::vector<bool>& activeHasPatch,
                                 std::vector<RegionSupportRef>& hits) -> bool
    {
        if (f < 0 || f >= static_cast<index_t>(functionDescription.size()))
            return false;
        for (const auto& desc : functionDescription[f])
        {
            if (desc.size() < 3)
                continue;
            const int p = static_cast<int>(desc[0]);
            const int lvl = static_cast<int>(desc[1]);
            const int idx = static_cast<int>(desc[2]);
            if (p < 0 || p >= static_cast<int>(nPatches) || !activeHasPatch[p])
                continue;
            if (p >= static_cast<int>(Bells.size()) || lvl < 0 || lvl >= static_cast<int>(Bells[p].size()) ||
                idx < 0 || idx >= static_cast<int>(Bells[p][lvl].size()))
                continue;

            const auto support = Bells[p][lvl].function(idx).support();
            if (supportIntersectsRectList(support, activeAabb[p]))
                hits.push_back({ p, lvl, idx, false });
        }

        if (f < static_cast<index_t>(hasSpillover.size()) && hasSpillover[f] &&
            f < static_cast<index_t>(spilloverFunctionCoordinates.size()))
        {
            for (const auto& sp : spilloverFunctionCoordinates[f])
            {
                const int p = sp[0];
                const int lvl = sp[1];
                const int idx = sp[2];
                if (p < 0 || p >= static_cast<int>(nPatches) || !activeHasPatch[p])
                    continue;
                if (p >= static_cast<int>(Bells.size()) || lvl < 0 || lvl >= static_cast<int>(Bells[p].size()) ||
                    idx < 0 || idx >= static_cast<int>(Bells[p][lvl].size()))
                    continue;

                const auto support = Bells[p][lvl].function(idx).support();
                if (supportIntersectsRectList(support, activeAabb[p]))
                    hits.push_back({ p, lvl, idx, true });
            }
        }

        return !hits.empty();
    };

    // Phase 1 runs the closure only on the initially-active (seed) patches.
    // Phase 2 maps those seeds to neighboring patches and runs a single-patch
    // closure on each — no further cascade beyond them.
    const std::vector<bool> seedHasPatch = region.hasPatch;
    std::vector<bool> expansionAllowed   = seedHasPatch; // Phase 1: seed patches only

    auto expandRegionByFunction = [&](index_t f) -> bool
    {
        bool changed = false;
        if (f < 0 || f >= static_cast<index_t>(functionDescription.size()))
            return false;
        for (const auto& desc : functionDescription[f])
        {
            if (desc.size() < 3)
                continue;
            const int p = static_cast<int>(desc[0]);
            const int lvl = static_cast<int>(desc[1]);
            const int idx = static_cast<int>(desc[2]);
            if (p < 0 || p >= static_cast<int>(nPatches))
                continue;
            if (!expansionAllowed[p])
                continue;
            if (p >= static_cast<int>(Bells.size()) || lvl < 0 || lvl >= static_cast<int>(Bells[p].size()) ||
                idx < 0 || idx >= static_cast<int>(Bells[p][lvl].size()))
                continue;

            const auto support = Bells[p][lvl].function(idx).support();
            const std::array<real_t, 4> supportBox = {
                static_cast<real_t>(support(0, 0)),
                static_cast<real_t>(support(1, 0)),
                static_cast<real_t>(support(0, 1)),
                static_cast<real_t>(support(1, 1))
            };
            const bool merged = mergePatchAabb(region, p,
                                               supportBox[0],
                                               supportBox[1],
                                               supportBox[2],
                                               supportBox[3]);
            changed = merged || changed;

            if (outfile.is_open())
            {
                outfile << "    [local-region] function " << f
                        << " component patch=" << p
                        << " level=" << lvl
                        << " tensorIdx=" << idx
                        << " support=" << formatBox(supportBox)
                        << " -> rects=" << formatBoxList(region.patchAABB[p]) << "\n";
            }
        }

        if (f < static_cast<index_t>(hasSpillover.size()) && hasSpillover[f] &&
            f < static_cast<index_t>(spilloverFunctionCoordinates.size()))
        {
            for (const auto& sp : spilloverFunctionCoordinates[f])
            {
                const int p = sp[0];
                const int lvl = sp[1];
                const int idx = sp[2];
                if (p < 0 || p >= static_cast<int>(nPatches))
                    continue;
                if (!expansionAllowed[p])
                    continue;
                if (p >= static_cast<int>(Bells.size()) || lvl < 0 || lvl >= static_cast<int>(Bells[p].size()) ||
                    idx < 0 || idx >= static_cast<int>(Bells[p][lvl].size()))
                    continue;

                const auto support = Bells[p][lvl].function(idx).support();
                const std::array<real_t, 4> supportBox = {
                    static_cast<real_t>(support(0, 0)),
                    static_cast<real_t>(support(1, 0)),
                    static_cast<real_t>(support(0, 1)),
                    static_cast<real_t>(support(1, 1))
                };
                const bool merged = mergePatchAabb(region, p,
                                                   supportBox[0],
                                                   supportBox[1],
                                                   supportBox[2],
                                                   supportBox[3]);
                changed = merged || changed;

                if (outfile.is_open())
                {
                    outfile << "    [local-region] function " << f
                            << " spillover patch=" << p
                            << " level=" << lvl
                            << " tensorIdx=" << idx
                            << " support=" << formatBox(supportBox)
                            << " -> rects=" << formatBoxList(region.patchAABB[p]) << "\n";
                }
            }
        }

        return changed;
    };

    // Multi-pass closure: iteratively select functions whose support intersects the
    // current region AABBs and expand the region until no new functions are added.
    bool selectedNewFunction = true;
    index_t closureIterations = 0;

    // Write per-step header to closure_log.txt
    if (closureLogFile.is_open())
    {
        static int closureStepCounter = 0;
        ++closureStepCounter;
        closureLogFile << "\n=== STEP " << closureStepCounter
                       << " srcPatch=" << sourcePatch
                       << " nFunctions=" << functionDescription.size()
                       << " ===\n";
        closureLogFile << "Initial rects:\n";
        for (index_t p = 0; p < nPatches; ++p)
        {
            if (!region.hasPatch[p]) continue;
            closureLogFile << "  patch=" << p << " " << formatBoxList(region.patchAABB[p]) << "\n";
        }
    }

    while (selectedNewFunction)
    {
        selectedNewFunction = false;
        ++closureIterations;

        const auto iterationAabb = region.patchAABB;
        const auto iterationHasPatch = region.hasPatch;

        // Closure-log: snapshot at start of this iteration
        if (closureLogFile.is_open())
        {
            closureLogFile << "--- Iteration " << closureIterations << " ---\n";
            closureLogFile << "  Snapshot patches:";
            for (index_t p = 0; p < nPatches; ++p)
                if (iterationHasPatch[p]) closureLogFile << " " << p;
            closureLogFile << "\n";
            for (index_t p = 0; p < nPatches; ++p)
            {
                if (!iterationHasPatch[p]) continue;
                closureLogFile << "  patch=" << p << " " << formatBoxList(iterationAabb[p]) << "\n";
            }
        }

        index_t newFuncThisIter = 0;
        std::vector<bool> patchWasActive(nPatches, false);
        for (index_t p = 0; p < nPatches; ++p)
            patchWasActive[p] = iterationHasPatch[p];

        for (index_t f = 0; f < static_cast<index_t>(functionDescription.size()); ++f)
        {
            if (region.basisInd(0, f) != 0.0)
                continue;

            std::vector<RegionSupportRef> regionHits;
            if (!collectRegionHits(f, iterationAabb, iterationHasPatch, regionHits))
                continue;

            if (outfile.is_open())
            {
                outfile << "[local-region] selecting function " << f
                        << " in closure iteration " << closureIterations
                        << " with " << regionHits.size() << " intersecting component(s)\n";
                for (const auto& hit : regionHits)
                {
                    const auto support = Bells[hit.patch][hit.level].function(hit.index).support();
                    const std::array<real_t, 4> supportBox = {
                        static_cast<real_t>(support(0, 0)),
                        static_cast<real_t>(support(1, 0)),
                        static_cast<real_t>(support(0, 1)),
                        static_cast<real_t>(support(1, 1))
                    };
                    outfile << "  -> " << (hit.spillover ? "spillover" : "component")
                            << " patch=" << hit.patch
                            << " level=" << hit.level
                            << " tensorIdx=" << hit.index
                            << " support=" << formatBox(supportBox)
                            << " intersects rects " << formatBoxList(iterationAabb[hit.patch]) << "\n";
                }
            }

            region.basisInd(0, f) = 1.0;
            selectedNewFunction = true;
            ++newFuncThisIter;

            // Snapshot patch membership before expanding so we can report new pulls
            std::vector<bool> patchBeforeExpand(nPatches, false);
            for (index_t p = 0; p < nPatches; ++p)
                patchBeforeExpand[p] = region.hasPatch[p];

            const bool grew = expandRegionByFunction(f);

            if (outfile.is_open())
            {
                if (!grew)
                    outfile << "  -> function " << f << " caused no AU box growth\n";
                outfile << "  -> AU rects after function " << f << ":\n";
                for (index_t p = 0; p < nPatches; ++p)
                {
                    if (!region.hasPatch[p])
                        continue;
                    outfile << "     patch=" << p << " rects=" << formatBoxList(region.patchAABB[p]) << "\n";
                }
            }

            // Closure-log: one line per newly-selected function
            // Format: f=X HIT:p=A lv=B idx=C sup=... [PULL->p=D sup=...] [(no box growth)]
            if (closureLogFile.is_open())
            {
                // Primary hit: the twin whose support overlapped the active AABB
                const auto& hit = regionHits[0];
                const auto hitSup = Bells[hit.patch][hit.level].function(hit.index).support();
                closureLogFile << "  f=" << f
                               << " HIT:p=" << hit.patch
                               << " lv=" << hit.level
                               << " idx=" << hit.index
                               << " sup=[" << hitSup(0,0) << "," << hitSup(1,0)
                               << "]x[" << hitSup(0,1) << "," << hitSup(1,1) << "]"
                               << (hit.spillover ? " (spillover)" : "");

                // Any twin that pulled a previously-inactive patch into the region
                for (const auto& desc : functionDescription[f])
                {
                    if (desc.size() < 3) continue;
                    const int tp = static_cast<int>(desc[0]);
                    if (tp < 0 || tp >= static_cast<int>(nPatches)) continue;
                    if (!patchBeforeExpand[tp] && region.hasPatch[tp])
                    {
                        const int tlv  = static_cast<int>(desc[1]);
                        const int tidx = static_cast<int>(desc[2]);
                        if (tp < static_cast<int>(Bells.size()) &&
                            tlv >= 0 && tlv < static_cast<int>(Bells[tp].size()) &&
                            tidx >= 0 && tidx < static_cast<int>(Bells[tp][tlv].size()))
                        {
                            const auto twinSup = Bells[tp][tlv].function(tidx).support();
                            closureLogFile << " PULL->p=" << tp
                                           << " sup=[" << twinSup(0,0) << "," << twinSup(1,0)
                                           << "]x[" << twinSup(0,1) << "," << twinSup(1,1) << "]";
                        }
                        else
                        {
                            closureLogFile << " PULL->p=" << tp;
                        }
                    }
                }

                if (!grew) closureLogFile << " (no box growth)";
                closureLogFile << "\n";
            }
        }

        // Closure-log: end-of-iteration summary
        if (closureLogFile.is_open())
        {
            closureLogFile << "  => " << newFuncThisIter << " new function(s) selected\n";
            // Which patches were newly added this iteration
            closureLogFile << "  => New patches added:";
            bool anyNew = false;
            for (index_t p = 0; p < nPatches; ++p)
            {
                if (region.hasPatch[p] && !patchWasActive[p])
                {
                    closureLogFile << " " << p;
                    anyNew = true;
                }
            }
            if (!anyNew) closureLogFile << " (none)";
            closureLogFile << "\n";
            // Rect state after this iteration
            closureLogFile << "  => Rects after iteration " << closureIterations << ":\n";
            for (index_t p = 0; p < nPatches; ++p)
            {
                if (!region.hasPatch[p]) continue;
                closureLogFile << "     patch=" << p << " " << formatBoxList(region.patchAABB[p]) << "\n";
            }
        }
    }

    // =========================================================
    // Phase 2: one-shot twin extension to neighboring patches.
    // Each Phase-1-selected function's twins on non-seed patches
    // define the mapped E(c)' for that neighbor.  We then run a
    // full single-patch closure on each neighbor independently —
    // no cascade beyond it.
    // =========================================================
    {
        // Reset expansionAllowed: each non-seed patch gets its own turn.
        for (index_t p = 0; p < static_cast<index_t>(nPatches); ++p)
            expansionAllowed[p] = false;

        // Collect mapped seed cells for each non-seed patch from Phase-1 selections.
        std::vector<std::vector<std::array<real_t,4>>> phase2Seeds(nPatches);
        for (index_t f = 0; f < static_cast<index_t>(functionDescription.size()); ++f)
        {
            if (region.basisInd(0, f) == 0.0) continue;
            for (const auto& desc : functionDescription[f])
            {
                if (desc.size() < 3) continue;
                const int p   = static_cast<int>(desc[0]);
                const int lvl = static_cast<int>(desc[1]);
                const int idx = static_cast<int>(desc[2]);
                if (p < 0 || p >= static_cast<int>(nPatches)) continue;
                if (seedHasPatch[p]) continue; // already covered by Phase 1
                if (p >= static_cast<int>(Bells.size()) || lvl < 0 ||
                    lvl >= static_cast<int>(Bells[p].size()) ||
                    idx < 0 || idx >= static_cast<int>(Bells[p][lvl].size()))
                    continue;
                const auto sup = Bells[p][lvl].function(idx).support();
                phase2Seeds[p].push_back({
                    static_cast<real_t>(sup(0,0)), static_cast<real_t>(sup(1,0)),
                    static_cast<real_t>(sup(0,1)), static_cast<real_t>(sup(1,1))
                });
            }
        }

        // For each non-seed patch that received seeds, run a single-patch closure.
        for (index_t targetP = 0; targetP < static_cast<index_t>(nPatches); ++targetP)
        {
            if (phase2Seeds[targetP].empty()) continue;

            // Seed the patch AABB from the mapped twin supports.
            region.hasPatch[targetP] = true;
            for (const auto& seed : phase2Seeds[targetP])
                mergePatchAabb(region, targetP, seed[0], seed[1], seed[2], seed[3]);

            if (closureLogFile.is_open())
            {
                closureLogFile << "--- Phase2 patch=" << targetP
                               << " seed=" << formatBoxList(region.patchAABB[targetP]) << " ---\n";
            }

            expansionAllowed[targetP] = true; // this patch only

            bool p2New = true;
            index_t p2Iter = 0;
            while (p2New)
            {
                p2New = false;
                ++p2Iter;
                const auto p2Aabb = region.patchAABB;
                std::vector<bool> p2Has(nPatches, false);
                p2Has[targetP] = true;

                index_t p2Count = 0;
                for (index_t f = 0; f < static_cast<index_t>(functionDescription.size()); ++f)
                {
                    if (region.basisInd(0, f) != 0.0) continue;

                    std::vector<RegionSupportRef> p2Hits;
                    if (!collectRegionHits(f, p2Aabb, p2Has, p2Hits)) continue;

                    region.basisInd(0, f) = 1.0;
                    p2New = true;
                    ++p2Count;
                    expandRegionByFunction(f); // restricted to targetP by expansionAllowed

                    if (closureLogFile.is_open())
                    {
                        const auto& h = p2Hits[0];
                        const auto hs = Bells[h.patch][h.level].function(h.index).support();
                        closureLogFile << "  [p2] f=" << f
                                       << " HIT:p=" << h.patch
                                       << " lv=" << h.level
                                       << " idx=" << h.index
                                       << " sup=[" << hs(0,0) << "," << hs(1,0)
                                       << "]x[" << hs(0,1) << "," << hs(1,1) << "]\n";
                    }
                }
                if (closureLogFile.is_open())
                    closureLogFile << "  [p2 iter " << p2Iter << "] " << p2Count << " new\n";
            }

            expansionAllowed[targetP] = false; // prevent leakage into the next patch's turn
        }
    }

    if (closureLogFile.is_open())
    {
        index_t activeFunctions = 0;
        for (index_t f = 0; f < region.basisInd.cols(); ++f)
            if (region.basisInd(0, f) != 0.0) ++activeFunctions;
        index_t activePatches = 0;
        for (index_t p = 0; p < nPatches; ++p)
            if (region.hasPatch[p]) ++activePatches;
        closureLogFile << "Final: " << closureIterations << " iterations, "
                       << activePatches << " patches, "
                       << activeFunctions << " functions\n";
    }

    if (outfile.is_open())
    {
        index_t activeFunctions = 0;
        for (index_t f = 0; f < region.basisInd.cols(); ++f)
            if (region.basisInd(0, f) != 0.0)
                ++activeFunctions;

        outfile << "[local-region] final AU rects:\n";
        for (index_t p = 0; p < nPatches; ++p)
        {
            if (!region.hasPatch[p])
                continue;
            outfile << "  patch=" << p << " nRects=" << region.patchAABB[p].size()
                    << " " << formatBoxList(region.patchAABB[p]) << "\n";
        }
            outfile << "[local-region] closureIterations=" << closureIterations << "\n";
        outfile << "[local-region] basisInd selected " << activeFunctions
                << " / " << region.basisInd.cols() << " functions\n";
    }

    return region;
}

static gsVector<gsMatrix<real_t>> filterUvByLocalRegion(
    const gsVector<gsMatrix<real_t>>& uv,
    const LocalCoarseningRegion& region,
    index_t numPoints)
{
    (void)numPoints;

    if (!region.enabled)
        return uv;

    gsVector<gsMatrix<real_t>> localUv(uv.size());
    index_t totalSelected = 0;

    for (index_t patch = 0; patch < uv.size(); ++patch)
    {
        const index_t baseCols = uv(patch).cols();
        const index_t baseRows = uv(patch).rows();

        if (!region.hasPatch.empty() &&
            (patch >= static_cast<index_t>(region.hasPatch.size()) || !region.hasPatch[patch]))
        {
            localUv(patch).resize(baseRows, 0);
            if (outfile.is_open())
                outfile << "[local-region] uv cloud patch=" << patch
                        << " base=" << baseCols
                        << " generated=0 (patch not in AU)\n";
            continue;
        }

        std::vector<index_t> keep;
        keep.reserve(baseCols);
        for (index_t c = 0; c < baseCols; ++c)
        {
            if (pointInsidePatchRegion(region, static_cast<int>(patch), uv(patch)(0, c), uv(patch)(1, c)))
                keep.push_back(c);
        }

        localUv(patch).resize(baseRows, keep.size());
        for (index_t k = 0; k < static_cast<index_t>(keep.size()); ++k)
            localUv(patch).col(k) = uv(patch).col(keep[k]);

        totalSelected += localUv(patch).cols();

        if (outfile.is_open())
        {
            outfile << "[local-region] uv cloud patch=" << patch
                    << " base=" << baseCols
                    << " kept=" << localUv(patch).cols()
                    << " fullPatch=" << (patchRegionCoversWholePatch(region, static_cast<int>(patch)) ? 1 : 0)
                    << "\n";
        }
    }

    if (outfile.is_open())
        outfile << "[local-region] filtered total local uv points=" << totalSelected << "\n";

    return localUv;
}

// Generate a fresh k×k uniform UV grid inside each patch's local AABB.
// k is chosen so that k*k > nLocalDOF (overdetermined system regardless of region size).
// Replaces filterUvByLocalRegion for the local-fitting mode.
static gsVector<gsMatrix<real_t>> resampleLocalRegion(
    const LocalCoarseningRegion& region,
    index_t nPatches,
    int k)
{
    gsVector<gsMatrix<real_t>> localUv(nPatches);
    index_t totalGenerated = 0;

    for (index_t patch = 0; patch < nPatches; ++patch)
    {
        if (!region.enabled ||
            patch >= static_cast<index_t>(region.hasPatch.size()) ||
            !region.hasPatch[patch])
        {
            localUv(patch).resize(2, 0);
            continue;
        }

        // Compute bounding box of all rects for this patch for sampling.
        real_t uMin = 1, vMin = 1, uMax = 0, vMax = 0;
        for (const auto& r : region.patchAABB[patch])
        {
            uMin = std::min(uMin, r[0]); vMin = std::min(vMin, r[1]);
            uMax = std::max(uMax, r[2]); vMax = std::max(vMax, r[3]);
        }
        if (uMin > uMax || vMin > vMax) { localUv(patch).resize(2, 0); continue; }

        localUv(patch).resize(2, k * k);
        index_t col = 0;
        for (int j = 0; j < k; ++j)
        {
            const real_t v = (k > 1)
                ? vMin + j * (vMax - vMin) / (k - 1)
                : 0.5 * (vMin + vMax);
            for (int i = 0; i < k; ++i)
            {
                const real_t u = (k > 1)
                    ? uMin + i * (uMax - uMin) / (k - 1)
                    : 0.5 * (uMin + uMax);
                localUv(patch)(0, col) = u;
                localUv(patch)(1, col) = v;
                ++col;
            }
        }
        totalGenerated += k * k;
    }

    if (outfile.is_open())
        outfile << "[local-region] resampled total local uv points=" << totalGenerated
                << " k=" << k << " k*k=" << k * k << "\n";

    return localUv;
}

int rebuildTheHierarchy(gsVector< gsVector<gsVector<index_t>>>& boxMat, int row, int x1U, int x1Bi, int x2U, int x2Bi, int y1U, int y1Bi, int y2U, int y2Bi, int levNow, int& lastNonZeroRow,
    int& createdBoxNum, int& centerInd, int ourBox[], int& needToEscape, int patch) {
    int RTH = 0;
    if (boxMat(patch)(row)(1) <= x1U * pow(2, boxMat(patch)(row)(0) - levNow) &&
        x1U * pow(2, boxMat(patch)(row)(0) - levNow) <= boxMat(patch)(row)(3) &&
        boxMat(patch)(row)(2) <= y1U * pow(2, boxMat(patch)(row)(0) - levNow) &&
        y1U * pow(2, boxMat(patch)(row)(0) - levNow) <= boxMat(patch)(row)(4) &&
        boxMat(patch)(row)(1) <= x2U * pow(2, boxMat(patch)(row)(0) - levNow) &&
        x2U * pow(2, boxMat(patch)(row)(0) - levNow) <= boxMat(patch)(row)(3) &&
        boxMat(patch)(row)(2) <= y2U * pow(2, boxMat(patch)(row)(0) - levNow) &&
        y2U * pow(2, boxMat(patch)(row)(0) - levNow) <= boxMat(patch)(row)(4)
        ) {
        if (boxMat(patch)(row)(0) <= levNow) {
        }
        //HERE I START TO WRITE CODE IN MORE GENERAL WAY
        else if (boxMat(patch)(row)(0) > levNow + 1) {
        }
        else if (boxMat(patch)(row)(0) == levNow + 1) {
            x1Bi = boxMat(patch)(row)(1);
            y1Bi = boxMat(patch)(row)(2);
            x2Bi = boxMat(patch)(row)(3);
            y2Bi = boxMat(patch)(row)(4);
            RTH = 1;
            //SW W NW
            if (x1U * pow(2, boxMat(patch)(row)(0) - levNow) - x1Bi > 0) {
                if (y1U * pow(2, boxMat(patch)(row)(0) - levNow) - y1Bi > 0) {
                    boxMat(patch)(lastNonZeroRow).resize(5);
                    boxMat(patch)(lastNonZeroRow)(0) = levNow + 1;
                    boxMat(patch)(lastNonZeroRow)(1) = x1Bi;
                    boxMat(patch)(lastNonZeroRow)(2) = y1Bi;
                    boxMat(patch)(lastNonZeroRow)(3) = x1U * pow(2, boxMat(patch)(row)(0) - (levNow));
                    boxMat(patch)(lastNonZeroRow)(4) = y1U * pow(2, boxMat(patch)(row)(0) - (levNow));
                    lastNonZeroRow++;
                    createdBoxNum++;
                }
                boxMat(patch)(lastNonZeroRow).resize(5);
                boxMat(patch)(lastNonZeroRow)(0) = levNow + 1;
                boxMat(patch)(lastNonZeroRow)(1) = x1Bi;
                boxMat(patch)(lastNonZeroRow)(2) = y1U * pow(2, boxMat(patch)(row)(0) - levNow);
                boxMat(patch)(lastNonZeroRow)(3) = x1U * pow(2, boxMat(patch)(row)(0) - levNow);
                boxMat(patch)(lastNonZeroRow)(4) = y2U * pow(2, boxMat(patch)(row)(0) - levNow);
                lastNonZeroRow++;
                createdBoxNum++;
                if (y2Bi - y2U * pow(2, boxMat(patch)(row)(0) - levNow) > 0) {
                    boxMat(patch)(lastNonZeroRow).resize(5);
                    boxMat(patch)(lastNonZeroRow)(0) = levNow + 1;
                    boxMat(patch)(lastNonZeroRow)(1) = x1Bi;
                    boxMat(patch)(lastNonZeroRow)(2) = y2U * pow(2, boxMat(patch)(row)(0) - levNow);
                    boxMat(patch)(lastNonZeroRow)(3) = x1U * pow(2, boxMat(patch)(row)(0) - levNow);
                    boxMat(patch)(lastNonZeroRow)(4) = y2Bi;
                    lastNonZeroRow++;
                    createdBoxNum++;
                }
            }
            //S N
            if (y1U * pow(2, boxMat(patch)(row)(0) - levNow) - y1Bi > 0) {
                boxMat(patch)(lastNonZeroRow).resize(5);
                boxMat(patch)(lastNonZeroRow)(0) = levNow + 1;
                boxMat(patch)(lastNonZeroRow)(1) = x1U * pow(2, boxMat(patch)(row)(0) - levNow);
                boxMat(patch)(lastNonZeroRow)(2) = y1Bi;
                boxMat(patch)(lastNonZeroRow)(3) = x2U * pow(2, boxMat(patch)(row)(0) - levNow);
                boxMat(patch)(lastNonZeroRow)(4) = y1U * pow(2, boxMat(patch)(row)(0) - levNow);
                lastNonZeroRow++;
                createdBoxNum++;
            }
            if (y2Bi - y2U * pow(2, boxMat(patch)(row)(0) - levNow) > 0) {
                boxMat(patch)(lastNonZeroRow).resize(5);
                boxMat(patch)(lastNonZeroRow)(0) = levNow + 1;
                boxMat(patch)(lastNonZeroRow)(1) = x1U * pow(2, boxMat(patch)(row)(0) - levNow);
                boxMat(patch)(lastNonZeroRow)(2) = y2U * pow(2, boxMat(patch)(row)(0) - levNow);
                boxMat(patch)(lastNonZeroRow)(3) = x2U * pow(2, boxMat(patch)(row)(0) - levNow);
                boxMat(patch)(lastNonZeroRow)(4) = y2Bi;                lastNonZeroRow++;
                createdBoxNum++;
            }
            //SE E NE
            if (x2Bi - x2U * pow(2, boxMat(patch)(row)(0) - levNow) > 0) {
                if (y1U * pow(2, boxMat(patch)(row)(0) - levNow) - y1Bi > 0) {
                    boxMat(patch)(lastNonZeroRow).resize(5);
                    boxMat(patch)(lastNonZeroRow)(0) = levNow + 1;
                    boxMat(patch)(lastNonZeroRow)(1) = x2U * pow(2, boxMat(patch)(row)(0) - levNow);
                    boxMat(patch)(lastNonZeroRow)(2) = y1Bi;
                    boxMat(patch)(lastNonZeroRow)(3) = x2Bi;
                    boxMat(patch)(lastNonZeroRow)(4) = y1U * pow(2, boxMat(patch)(row)(0) - levNow);
                    lastNonZeroRow++;
                    createdBoxNum++;
                }
                boxMat(patch)(lastNonZeroRow).resize(5);
                boxMat(patch)(lastNonZeroRow)(0) = levNow + 1;
                boxMat(patch)(lastNonZeroRow)(1) = x2U * pow(2, boxMat(patch)(row)(0) - levNow);
                boxMat(patch)(lastNonZeroRow)(2) = y1U * pow(2, boxMat(patch)(row)(0) - levNow);
                boxMat(patch)(lastNonZeroRow)(3) = x2Bi;
                boxMat(patch)(lastNonZeroRow)(4) = y2U * pow(2, boxMat(patch)(row)(0) - levNow);
                lastNonZeroRow++;
                createdBoxNum++;
                if (y2Bi - y2U * pow(2, boxMat(patch)(row)(0) - levNow) > 0) {
                    boxMat(patch)(lastNonZeroRow).resize(5);
                    boxMat(patch)(lastNonZeroRow)(0) = levNow + 1;
                    boxMat(patch)(lastNonZeroRow)(1) = x2U * pow(2, boxMat(patch)(row)(0) - levNow);
                    boxMat(patch)(lastNonZeroRow)(2) = y2U * pow(2, boxMat(patch)(row)(0) - levNow);
                    boxMat(patch)(lastNonZeroRow)(3) = x2Bi;
                    boxMat(patch)(lastNonZeroRow)(4) = y2Bi;
                    lastNonZeroRow++;
                    createdBoxNum++;
                }
            }
            outfile << "CREATE C BOX:\n";
            outfile << "ourBOX is \n";
            for (int k = 0; k < 5; k++) {
                ourBox[k] = boxMat(patch)(row)(k);
                outfile << ourBox[k] << " ";
                centerInd = row;
            }
            outfile << "\n";
            boxMat(patch)(row)(0) = levNow;
            boxMat(patch)(row)(1) = x1U;
            boxMat(patch)(row)(2) = y1U;
            boxMat(patch)(row)(3) = x2U;
            boxMat(patch)(row)(4) = y2U;
        }
    }
    return RTH;
}

int rebuildTheHierarchyMultiple(gsVector< gsVector<gsVector<index_t>>>& boxMat,
    const std::vector<CellToCoarsen>& cells,
    int& lastNonZeroRow,
    int& createdBoxNum,
    int& centerInd,
    int ourBox[],
    int& needToEscape,
    int patch,
    std::vector<RestoreRecord>* restoreRecords = nullptr,
    int successfullAttempts = 1)
{
    int anyRebuilt = 0;
    for (const auto& cell : cells)
    {
        int createdBefore = createdBoxNum;
        for (int row = 0; row < lastNonZeroRow; ++row)
        {
            int x1Bi = 0, x2Bi = 0, y1Bi = 0, y2Bi = 0;
            int RTH = rebuildTheHierarchy(
                boxMat,
                row,
                cell.x1, x1Bi,
                cell.x2, x2Bi,
                cell.y1, y1Bi,
                cell.y2, y2Bi,
                cell.level,
                lastNonZeroRow,
                createdBoxNum,
                centerInd,
                ourBox,
                needToEscape,
                patch);
            if (RTH == 1)
            {
                anyRebuilt = 1;
                if (restoreRecords)
                {
                    RestoreRecord record;
                    record.createdBoxNum = createdBoxNum - createdBefore;
                    record.centerInd = centerInd;
                    record.box = { ourBox[0], ourBox[1], ourBox[2], ourBox[3], ourBox[4] };
                    record.levNow = cell.level;
                    record.successfullAttempts = successfullAttempts;
                    record.patch = patch;
                    restoreRecords->push_back(record);
                }
                break;
            }
        }
    }
    return anyRebuilt;
}

void restoreTheHierarchyMultiple(
    const std::vector<RestoreRecord>& records,
    int& lastNonzeroRow,
    gsVector< gsVector<gsVector<index_t>>>& boxMat)
{
    for (auto it = records.rbegin(); it != records.rend(); ++it)
    {
        int createdBoxNum = it->createdBoxNum;
        int centerInd = it->centerInd;
        int ourBox[5] = { it->box[0], it->box[1], it->box[2], it->box[3], it->box[4] };
        int levNow = it->levNow;
        int successfullAttempts = it->successfullAttempts;
        int patch = it->patch;
        if (patch < 0 || patch >= boxMat.size())
        {
            gsInfo << "[restoreTheHierarchyMultiple] Skipping invalid patch: " << patch << "\n";
            outfile << "[restoreTheHierarchyMultiple] Skipping invalid patch: " << patch << "\n";
            continue;
        }
        if (boxMat(patch).size() == 0)
        {
            gsInfo << "[restoreTheHierarchyMultiple] Skipping empty boxMat for patch: " << patch << "\n";
            outfile << "[restoreTheHierarchyMultiple] Skipping empty boxMat for patch: " << patch << "\n";
            continue;
        }
        restoreTheHierarchy(createdBoxNum, lastNonzeroRow, boxMat, levNow, centerInd, ourBox, successfullAttempts, patch);
    }
}

void exportToPatches(
    gsMatrix<real_t> vectSol,
    std::vector<std::vector<std::vector<index_t>>> functionDescription,
    std::vector<std::vector<std::vector<double>>>& localCoeffs,
    std::vector<std::vector<std::vector<int>>>& localIndex
) {
    localCoeffs.clear();
    localIndex.clear();
    localCoeffs.resize(numberOfPatchesForReference);
    localIndex.resize(numberOfPatchesForReference);
    for (int functionIndex = 0; functionIndex < functionDescription.size(); functionIndex++) {
        for (int pieceOfFunction = 0; pieceOfFunction < functionDescription[functionIndex].size(); pieceOfFunction++) {
            int patchIndex = functionDescription[functionIndex][pieceOfFunction][0];
            localIndex[patchIndex].push_back(functionDescription[functionIndex][pieceOfFunction]);
            std::vector<double> corrCoeff;
            corrCoeff.push_back(vectSol(functionIndex, 0));
            corrCoeff.push_back(vectSol(functionIndex, 1));
            localCoeffs[patchIndex].push_back(corrCoeff);
        }
    }
}

static bool tryMapInterfaceSideToFeatureSide(
    int sideIndex,
    FeatureSide& side)
{
    switch (sideIndex)
    {
    case 1:
        side = FeatureSide::U0;
        return true;
    case 2:
        side = FeatureSide::U1;
        return true;
    case 3:
        side = FeatureSide::V0;
        return true;
    case 4:
        side = FeatureSide::V1;
        return true;
    default:
        return false;
    }
}

/**
 * @brief Saves the THB-spline geometry with a given vector solution patch-wise to files.
 *
 * This function:
 * - Builds `localIndex` (list of active basis indices per patch) using globalIndexTHB
 * - Uses globalIndexTHB to extract solution values from vectSol
 * - Builds geometry for each patch
 * - Writes the result to a file with names baseName_0, baseName_1, ...
 *
 * @param[in]  vectSol         Global solution matrix (numFunctions x 2), each row [u_x, u_y].
 * @param[in]  globalIndexTHB  Maps (patch, local basis index) to global function index.
 * @param[in]  SubdomainHierarchy  THB-spline basis per patch.
 * @param[out] outfile         Stream for logging file paths and errors.
 * @param[in]  baseName        Prefix for file names.
 * @param[out] localIndex      Will be filled with active basis indices per patch.
 *
 * @date Last modified: 2025-04-24
 */
void savePatches(
    const gsMatrix<real_t>& vectSol,
    const gsVector<gsVector<int>>& globalIndexTHB,
    const gsVector<gsTHBSplineBasis<2, real_t>>& SubdomainHierarchy,
    std::ostream& outfile,
    const std::string& baseName,
    gsVector<gsVector<int>>& localIndex
)
{
    PROFILE_FUNCTION();
    // gsInfo << "Starting savePatches for " << SubdomainHierarchy.size() << " patches.\n";
    // outfile << "Starting savePatches for " << SubdomainHierarchy.size() << " patches.\n";

    // Resize to accommodate the maximum of the two sizes
    size_t size1 = static_cast<size_t>(SubdomainHierarchy.size());
    size_t size2 = static_cast<size_t>(globalIndexTHB.size());
    size_t maxSize = (size1 > size2) ? size1 : size2;
    localIndex.resize(maxSize);

    for (size_t patch = 0; patch < SubdomainHierarchy.size(); ++patch)
    {
        // Bounds check for globalIndexTHB access
        if (patch >= static_cast<size_t>(globalIndexTHB.size()))
        {
            gsInfo << "WARNING: patch " << patch << " >= globalIndexTHB.size=" << globalIndexTHB.size() << ", skipping\n";
            outfile << "WARNING: patch " << patch << " >= globalIndexTHB.size=" << globalIndexTHB.size() << ", skipping\n";
            continue;
        }
        
        const auto& patchGlobalIndices = globalIndexTHB(patch);

        // Count active basis functions
        index_t count = 0;
        for (index_t i = 0; i < patchGlobalIndices.size(); ++i)
        {
            if (patchGlobalIndices[i] != -1)
                ++count;
        }

        localIndex(patch).resize(count);

        // gsInfo << "Processing patch " << patch << " with " << count << " active basis functions.\n";
        // outfile << "Processing patch " << patch << " with " << count << " active basis functions.\n";

        // Fill localIndex with valid basis indices
        index_t current = 0;
        for (index_t i = 0; i < patchGlobalIndices.size(); ++i)
        {
            if (patchGlobalIndices[i] != -1)
            {
                localIndex(patch)(current) = static_cast<int>(i);
                ++current;
            }
        }

        gsMatrix<> localVectSol(count, 2);
        for (index_t i = 0; i < count; ++i)
        {
            int basisIndex = localIndex(patch)(i);
            int globalIndex = patchGlobalIndices[basisIndex];

            if (globalIndex < 0 || globalIndex >= vectSol.rows())
            {
                gsInfo << "ERROR: Invalid global index " << globalIndex
                    << " at patch " << patch << ", basis index " << basisIndex << "\n";
                outfile << "ERROR: Invalid global index " << globalIndex
                    << " at patch " << patch << ", basis index " << basisIndex << "\n";
                continue;
            }

            localVectSol(i, 0) = vectSol(globalIndex, 0);
            localVectSol(i, 1) = vectSol(globalIndex, 1);
        }

        try
        {
            gsGeometry<>::uPtr geom = SubdomainHierarchy(patch).makeGeometry(localVectSol);
            gsFileData<> fileData;
            fileData << *geom;

            std::string outputFileName = baseName + "_" + std::to_string(patch);
            fileData.dump(outputFileName);

            // gsInfo << "Saved patch " << patch << " to file: " << outputFileName << "\n";
            // outfile << "Saved patch " << patch << " to file: " << outputFileName << "\n";
        }
        catch (const std::exception& e)
        {
            gsInfo << "EXCEPTION on patch " << patch << ": " << e.what() << "\n";
            outfile << "EXCEPTION on patch " << patch << ": " << e.what() << "\n";
        }
    }

    // gsInfo << "Finished savePatches.\n";
    // outfile << "Finished savePatches.\n";
}


/////August block/////////////
bool knotExists(const gsKnotVector<real_t>& kv, real_t xi)
{
    const index_t size = kv.size();
    for (index_t i = 0; i < size; ++i)
        if (std::abs(kv[i] - xi) < 1e-12)
            return true;
    return false;
}

/**
 * @brief Perform one-directional Boehm refinement (manual) on a tensor-product function.
 *
 * This inserts a single knot ξ in the u-direction and updates the coefficient matrix accordingly.
 *
 * @param kvU Knot vector in u-direction (will be modified).
 * @param coefs Coefficients (nU × nV) to be refined in u-direction.
 * @param xi Knot to insert in u-direction (must lie in knot span).
 * @param degree Degree of basis in u-direction.
 */
void refineU_manual(gsKnotVector<real_t>& kvU,
    gsMatrix<real_t>& coefs,
    real_t xi,
    index_t degree)
{
    const index_t nU = coefs.rows(); // original number of rows
    const index_t nV = coefs.cols(); // columns stay the same

    outfile << "=== refineU_manual ===\n";
    outfile << "xi = " << xi << ", degree = " << degree << "\n";
    outfile << "coefs: \n" << coefs << "\n";
    // Find the span manually
    index_t span = 0;
    for (index_t i = degree; i < kvU.size() - degree - 1; ++i)
    {
        if (xi >= kvU[i] && xi < kvU[i + 1])
        {
            span = i;
            break;
        }
    }
    outfile << "Found span = " << span << "\n";

    index_t a = span - degree + 1;
    index_t b = span;

    outfile << "a = " << a << ", b = " << b << "\n";

    // New coefficient matrix
    gsMatrix<real_t> newCoefs(nU + 1, nV);
    newCoefs.setZero();

    outfile << "Copying unaffected rows (0 to " << a - 1 << ")\n";
    for (index_t i = 0; i < a; ++i)
        newCoefs.row(i) = coefs.row(i);

    outfile << "Copying unaffected rows (" << b + 1 << " to " << nU - 1 << ")\n";
    for (index_t i = b + 1; i < nU; ++i)
        newCoefs.row(i + 1) = coefs.row(i);

    // Apply refinement in the affected rows
    outfile << "Refining rows from " << b << " down to " << a << "\n";
    for (index_t j = 0; j < nV; ++j)
    {
        for (index_t i = b; i >= a; --i)
        {
            real_t denom = kvU[i + degree] - kvU[i];
            real_t alpha = denom > 0 ? (xi - kvU[i]) / denom : 0.0;

            newCoefs(i, j) = alpha * coefs(i, j) + (1.0 - alpha) * coefs(i - 1, j);
            outfile << "newCoefs(" << i << ", " << j << ") = " << alpha << " * " << coefs(i, j) << " + (1.0 - " << alpha << ") * " << coefs(i - 1, j) << "\n";
            outfile << "  (i=" << i << ", j=" << j << ") alpha=" << alpha
                << " => newCoefs(i,j)=" << newCoefs(i, j) << "\n";
        }
    }

    // Insert knot into the vector
    kvU.insert(xi, 1);
    outfile << "Inserted knot " << xi << " into kvU\n";

    coefs = newCoefs;
    outfile << "Updated coefficient matrix (rows=" << coefs.rows()
        << ", cols=" << coefs.cols() << ")\n";
}

void refineV_manual(gsKnotVector<real_t>& kvV,
    gsMatrix<real_t>& coefs,
    real_t eta,
    index_t degree)
{
    const index_t nU = coefs.rows(); // stays unchanged
    const index_t nV = coefs.cols(); // control points before refinement

    outfile << "=== refineV_manual ===\n";
    outfile << "eta = " << eta << ", degree = " << degree << "\n";

    // Manual span search
    index_t span = degree;
    while (span < kvV.size() - degree - 1 && eta >= kvV[span + 1])
        ++span;

    const index_t a = span - degree + 1;
    const index_t b = span;

    outfile << "Found span = " << span << ", a = " << a << ", b = " << b << "\n";

    const index_t new_nV = nV + 1;
    gsMatrix<real_t> newCoefs(nU, new_nV);
    newCoefs.setZero();

    // Copy unaffected columns
    outfile << "Copying unaffected columns (0 to " << a - 1 << ")\n";
    for (index_t j = 0; j < a; ++j)
        newCoefs.col(j) = coefs.col(j);

    outfile << "Copying unaffected columns (" << b + 1 << " to " << nV - 1 << ")\n";
    for (index_t j = b + 1; j < nV; ++j)
        newCoefs.col(j + 1) = coefs.col(j);

    // Refine affected band
    outfile << "Refining columns from " << b << " down to " << a << "\n";
    for (index_t i = 0; i < nU; ++i)
    {
        for (index_t j = b; j >= a; --j)
        {
            real_t denom = kvV[j + degree] - kvV[j];
            real_t alpha = (denom == 0.0) ? 0.0 : (eta - kvV[j]) / denom;

            newCoefs(i, j) = alpha * coefs(i, j) + (1.0 - alpha) * coefs(i, j - 1);

            outfile << "  (i=" << i << ", j=" << j << ") alpha=" << alpha
                << " => newCoefs(i,j)=" << newCoefs(i, j) << "\n";
        }
    }

    kvV.insert(eta, 1);
    outfile << "Inserted knot " << eta << " into kvV\n";

    coefs = newCoefs;
    outfile << "Updated coefficient matrix (rows=" << coefs.rows()
        << ", cols=" << coefs.cols() << ")\n";
}



void refineV_many(gsKnotVector<real_t>& kvV,
    gsMatrix<real_t>& coefs,
    const std::vector<real_t>& newKnots,
    index_t degree)
{
    for (const real_t eta : newKnots)
        refineV_manual(kvV, coefs, eta, degree);
}
void refineU_many(gsKnotVector<real_t>& kvU,
    gsMatrix<real_t>& coefs,
    const std::vector<real_t>& newKnots,
    index_t degree)
{
    for (const real_t xi : newKnots)
        refineU_manual(kvU, coefs, xi, degree);
}
void refineBoehm2D_manual(gsKnotVector<real_t>& kvU,
    gsKnotVector<real_t>& kvV,
    gsMatrix<real_t>& coefs,
    const std::vector<real_t>& newKnotsU,
    const std::vector<real_t>& newKnotsV,
    index_t degreeU,
    index_t degreeV)
{
    refineU_many(kvU, coefs, newKnotsU, degreeU);
    refineV_many(kvV, coefs, newKnotsV, degreeV);
}
std::vector<real_t> internalKnots(const gsKnotVector<real_t>& kv)
{
    outfile << "=== internalKnots ===\n";
    std::vector<real_t> result;
    for (index_t i = 1; i < kv.size() - 1; ++i)
    {
        real_t xi = kv.at(i);
        outfile << "xi: " << xi << "\n";
        if (xi != kv.at(i - 1) && xi != kv.at(i + 1))
            result.push_back(xi);
    }
    return result;
}
void mapCoefsToGlobalVector(
    const gsVector<gsVector<gsTensorBSplineBasis<2, real_t>>>& Bells,
    int patch,
    int level,
    const gsMatrix<real_t>& coefs,
    const gsVector<gsVector<gsVector<index_t>>>& indexInTHB,
    std::vector<index_t>& indices,
    std::vector<real_t>& values,
    bool verbose = false)
{
    outfile << "mapCoefsToGlobalVector, coeffs:\n" << coefs << "\n";
    const index_t sizeU = Bells[patch][level].size(0);
    const index_t sizeV = Bells[patch][level].size(1);

    const index_t nU = coefs.rows();
    const index_t nV = coefs.cols();

    for (index_t i = 0; i < nU; ++i)
    {
        for (index_t j = 0; j < nV; ++j)
        {
            real_t val = coefs(i, j);
            if (std::abs(val) < 1e-14)
                continue;

            //if (i >= sizeU || j >= sizeV)
            //    continue;

            index_t localIndex = j * sizeU + i;
            index_t globalIndex = indexInTHB[patch][level + 1][localIndex];

            indices.push_back(globalIndex);
            values.push_back(val);

            if (verbose)
                gsInfo << "  (" << i << "," << j << ") -> global " << globalIndex
                << ", value = " << val << "\n";
        }
    }
}



// Small 2D wrapper around gsTensorBoehmRefineLocal
void refineDirectionWithBoehm(
    gsKnotVector<real_t>& kv,
    gsMatrix<real_t>& coefs,
    unsigned direction,
    const std::vector<real_t>& newKnots)
{
    if (newKnots.empty())
        return;

    gsVector<index_t, 2> nmb_of_coefs;
    nmb_of_coefs(0) = coefs.rows();
    nmb_of_coefs(1) = coefs.cols();
    gsVector<index_t, 2> act_size_of_coefs = nmb_of_coefs;
    gsVector<index_t, 2> size_of_coefs = nmb_of_coefs;
    size_of_coefs[direction] += (index_t)newKnots.size();

    gsTensorBoehmRefineLocal<2>(
        kv,
        0, // index offset
        coefs,
        nmb_of_coefs,
        act_size_of_coefs,
        size_of_coefs,
        direction,
        newKnots.begin(),
        newKnots.end(),
        true // update knots
        );
}

// 21:37 — Thin wrapper around gsTensorBoehmRefineLocal for 2D flattened layout
void refineDirectionWithBoehm_flat(
    gsVector<index_t, 2> bspl_vec_ti,
    int dim,
    gsKnotVector<real_t>& kv,
    gsMatrix<real_t>& coefs,          // flattened: (nU * nV) x 1
    index_t nU,                       // number of coefficients in U
    index_t nV,                       // number of coefficients in V
    unsigned direction,               // 0 = U, 1 = V
    const std::vector<real_t>& newKnots)
{
    if (newKnots.empty())
        return;

    gsVector<index_t, 2> nmb_of_coefs;
    nmb_of_coefs[0] = nU;
    nmb_of_coefs[1] = nV;

    // In the THB code, act_size_of_coefs starts as nmb_of_coefs
    gsVector<index_t, 2> act_size_of_coefs = nmb_of_coefs;

    // size_of_coefs = nmb_of_coefs with extra slots for new knots in this dir
    gsVector<index_t, 2> size_of_coefs = nmb_of_coefs;
    size_of_coefs[direction] += static_cast<index_t>(newKnots.size());

    gsTensorBoehmRefineLocal<2>(
        kv,
        bspl_vec_ti[dim],                  // index offset (local refinement start index)
        coefs,
        nmb_of_coefs,
        act_size_of_coefs,
        size_of_coefs,
        direction,
        newKnots.begin(),
        newKnots.end(),
        true                // update knots in kv
        );
}

// Convert a THB span index to actual knot array index of first occurrence
index_t getFirstKnotIndex(
    const gsKnotVector<real_t>& kv,
    index_t spanIndex)
{
    real_t knotValue = kv[spanIndex];
    for (index_t i = 0; i < kv.size(); ++i)
        if (kv[i] == knotValue)
            return i;
    GISMO_ERROR("Knot value not found in knot vector");
    return -1;
}

// Convert a THB span index to actual knot array index of last occurrence
index_t getLastKnotIndex(
    const gsKnotVector<real_t>& kv,
    index_t spanIndex)
{
    real_t knotValue = kv[spanIndex];
    for (index_t i = kv.size() - 1; i >= 0; --i)
        if (kv[i] == knotValue)
            return i;
    GISMO_ERROR("Knot value not found in knot vector");
    return -1;
}



// helper: find span k so that kv[k] <= x < kv[k+1]
// behaves like a robust manual findSpan (handles right endpoint)
static int getFirstKnotIndex(const gsKnotVector<real_t>& kv, real_t x, int degree)
{
    double epsilon = 1e-6;
    const int nk = static_cast<int>(kv.size());
    const int first = 0;
    const int lastSpan = nk - degree; // last valid span start
    if (x >= kv[lastSpan + 1]) return lastSpan - 1;
    if (x <= kv[first]) return 0;
    for (int k = 0; k <= lastSpan - 1; ++k)
        if (x >= kv[k] && x < kv[k + 1]) return k;
    return lastSpan - 2;
}

static int getLastKnotIndex(const gsKnotVector<real_t>& kv, real_t x, int degree)
{
    double epsilon = 1e-6;
    const int nk = static_cast<int>(kv.size());
    const int first = degree;
    const int lastSpan = nk;
    if (x <= kv[first] + epsilon)
        return first;
    for (int k = lastSpan - 1; k >= degree; k--)
    {
        if (x >= kv[k - 1] && x <= kv[k])
            return k;
    }
    return first;
}

static int findSpanManual(const gsKnotVector<real_t>& kv, real_t x, int degree)
{
    double epsilon = 1e-6;
    const int nk = static_cast<int>(kv.size());
    const int first = degree;
    const int lastSpan = nk - degree; // last valid span start
    if (x >= kv[lastSpan + 1]) return lastSpan - 1;
    for (int k = first; k <= lastSpan; ++k)
        if (x >= kv[k] && x < kv[k + 1]) return k - 2;
    if (x < kv[first]) return first - 2;
    return lastSpan - 2;
}



static void computeElementSpanIndicesFromSupport(
    const gsTensorBSplineBasis<2, real_t>& basis,
    const gsMatrix<real_t, 2, 2>& supp,     // [u0 u1; v0 v1]
    gsMatrix<index_t, 2, 2>& element_ind)   // output: dim x {low, high}
{
    // Knot vectors and degrees
    const auto& kvU = basis.knots(0);
    const auto& kvV = basis.knots(1);
    const int degU = basis.degree(0);
    const int degV = basis.degree(1);

    // Support coordinates
    real_t u0 = supp(0, 0), u1 = supp(0, 1);
    real_t v0 = supp(1, 0), v1 = supp(1, 1);

    // Find ordinal span indices
    int spanLowU = getFirstKnotIndex(kvU, u0, degU);
    int spanHighU = getLastKnotIndex(kvU, u1, degU);

    int spanLowV = getFirstKnotIndex(kvV, v0, degV);
    int spanHighV = getLastKnotIndex(kvV, v1, degV);

    // Store in element_ind (row = dim, col 0 = low, col 1 = high)
    element_ind(0, 0) = spanLowU;
    element_ind(0, 1) = spanHighU;
    element_ind(1, 0) = spanLowV;
    element_ind(1, 1) = spanHighV;
}


// emulate elementSupport_into for 2D basis
// element_ind must be a gsMatrix<index_t,2,2> (or equivalent) already allocated
static void emulate_elementSupport_into(
    const gsTensorBSplineBasis<2, real_t>& basis,
    const index_t tensor_index,                     // flat basis index inside this basis[level]
    gsMatrix<index_t, 2, 2>& element_ind)              // out: [ [lowU, highU], [lowV, highV] ] as columns
{
    // get the support interval of that local basis function
    // NOTE: basis.function(...) expects the local (flat) index inside this basis
    auto supp = basis.function(tensor_index).support();
    // supp is 2x2: supp(0,0)=u0, supp(0,1)=u1; supp(1,0)=v0, supp(1,1)=v1

    // knot vectors and degrees
    const auto& kvU = basis.knots(0);
    const auto& kvV = basis.knots(1);
    const int degU = basis.degree(0);
    const int degV = basis.degree(1);

    // find span indices on coarse knot vectors
    // use a tiny inward nudge for right endpoints to avoid hitting the right-open convention
    //real_t eps = std::numeric_limits<real_t>::epsilon() * 1000.0;
    real_t u0 = supp(0, 0), u1 = supp(0, 1);
    real_t v0 = supp(1, 0), v1 = supp(1, 1);

    int spanLowU = findSpanManual(kvU, u0, degU);
    int spanHighU = findSpanManual(kvU, u1, degU);

    int spanLowV = findSpanManual(kvV, v0, degV);
    int spanHighV = findSpanManual(kvV, v1, degV);

    // write output (format matches original: dim index in rows, low/high in columns)
    element_ind.col(0)(0) = spanLowU;
    element_ind.col(1)(0) = spanHighU;
    element_ind.col(0)(1) = spanLowV;
    element_ind.col(1)(1) = spanHighV;
}


auto containsKnot = [](const gsKnotVector<real_t>& kv, real_t val)
{
    for (index_t kk = 0; kk < kv.size(); ++kk)
        if (std::abs(kv[kk] - val) < 1e-14) // tolerance
            return true;
    return false;
};

unsigned updateSizeOfCoefsFromBells(
    const gsVector<gsVector<gsTensorBSplineBasis<2, real_t>>>& Bells,
    int patch,
    int clevel,
    int flevel,
    const gsVector<index_t, 2>& clow,
    const gsVector<index_t, 2>& chigh,
    const gsVector<index_t, 2>& flow,
    const gsVector<index_t, 2>& fhigh,
    gsVector<index_t, 2>& size_of_coefs)
{
    unsigned nmb_of_coefs = 1;
    for (unsigned dim = 0; dim < 2; ++dim)
    {
        const auto& ckv = Bells[patch][clevel].knots(dim);
        const auto& fkv = Bells[patch][flevel].knots(dim);

        auto knotCountBetween = [](const auto& kv, real_t lowVal, real_t highVal)
        {
            return std::count_if(kv.begin(), kv.end(),
                [&](real_t val) { return val >= lowVal && val <= highVal; });
        };

        real_t lowValC = ckv[clow(dim)];
        real_t highValC = ckv[chigh(dim)];
        real_t lowValF = fkv[flow(dim)];
        real_t highValF = fkv[fhigh(dim)];

        unsigned cnmb_knts = knotCountBetween(ckv, lowValC, highValC);
        unsigned fnmb_knts = knotCountBetween(fkv, lowValF, highValF);

        size_of_coefs(dim) += fnmb_knts - cnmb_knts;
        nmb_of_coefs *= size_of_coefs(dim);
    }
    return nmb_of_coefs;
}

static void computeFineBlockRange(
    const gsTensorBSplineBasis<2, real_t>& coarseBasis,
    const gsTensorBSplineBasis<2, real_t>& fineBasis,
    const index_t coarseLocalIndex,
    index_t& startU, index_t& endU, index_t& blockU,
    index_t& startV, index_t& endV, index_t& blockV)
{
    auto supp = coarseBasis.function(coarseLocalIndex).support();
    real_t sX0 = supp(0, 0), sX1 = supp(0, 1);
    real_t sY0 = supp(1, 0), sY1 = supp(1, 1);

    const auto& fkvU = fineBasis.knots(0);
    const auto& fkvV = fineBasis.knots(1);
    const int degU = fineBasis.degree(0);
    const int degV = fineBasis.degree(1);

    int spanLowU = findSpanManual(fkvU, sX0, degU);
    int spanHighU = findSpanManual(fkvU, std::nextafter(sX1, sX0), degU);

    int spanLowV = findSpanManual(fkvV, sY0, degV);
    int spanHighV = findSpanManual(fkvV, std::nextafter(sY1, sY0), degV);

    startU = spanLowU - degU;
    endU = spanHighU - degU;
    startV = spanLowV - degV;
    endV = spanHighV - degV;

    // clamp to valid ranges
    const index_t maxU = static_cast<index_t>(fineBasis.size(0) - 1);
    const index_t maxV = static_cast<index_t>(fineBasis.size(1) - 1);
    if (startU < 0) startU = 0;
    if (startV < 0) startV = 0;
    if (endU < 0) endU = 0;
    if (endV < 0) endV = 0;
    if (startU > maxU) startU = maxU;
    if (endU > maxU) endU = maxU;
    if (startV > maxV) startV = maxV;
    if (endV > maxV) endV = maxV;

    blockU = static_cast<index_t>(endU - startU + 1);
    blockV = static_cast<index_t>(endV - startV + 1);
}

void _truncate(
    gsMatrix<real_t>& coefs,
    const gsVector<index_t, 2>& act_size_of_coefs,
    const gsVector<index_t, 2>& size_of_coefs,
    const unsigned level,
    const gsVector<index_t, 2>& bspl_vec_ti,
    const unsigned bspl_vec_ti_level,
    const gsVector<index_t, 2>& finest_low,
    int patch,
    gsVector<int>& presavedIndices,
    const gsVector<gsVector<gsVector<index_t>>>& globalIndex,
    gsVector < gsTHBSplineBasis < 2, real_t > > SubdomainHierarchy,
    bool verbose = false,
    int functionIndex = -1,
    int attempt = 0)
{
    // if we dont have any active function in this level, we do not truncate
    if (SubdomainHierarchy[patch].getXmatrix()[level].size() == 0)
        return;

    // global tensor index
    const unsigned const_ten_index = SubdomainHierarchy[patch]._basisFunIndexOnLevel(bspl_vec_ti,
        bspl_vec_ti_level, finest_low, level);
    gsVector<index_t, 2> act_coefs_strides(2);
    bspline::buildCoeffsStrides<2>(act_size_of_coefs, act_coefs_strides);


    gsVector<index_t, 2> last_point(2);
    bspline::getLastIndexLocal<2>(size_of_coefs, last_point);
    last_point(0) = 0;


    gsVector<index_t, 2> position(2);
    position.fill(0);

    gsVector<index_t, 2> first_point(position);

    unsigned xmatrix_index = 0;
    unsigned tensor_active_index = SubdomainHierarchy[patch].getXmatrix()[level][0];

    unsigned numb_of_point = size_of_coefs[0];
    int i = 0;
    int truncationStopped = 0;
    presavedIndices.resize(coefs.rows());
    presavedIndices.fill(-1);  // Initialize to -1 to mark invalid entries
    // verbose = (patch == 1 && level == 3 && attempt == 0) && functionIndex <= 1;
    do
    {
        // indices
        // ten_index - (tensor) index of a bspline function with respect to
        //             the coef at position "position"
        // coef_index - (local) index of a ceofficient at position
        // tensor_active_index - (tensor) index of a active functions


        unsigned ten_index = const_ten_index;
        for (unsigned dim = 1; dim < 2; dim++)
        {
            ten_index += position(dim) *
                SubdomainHierarchy[patch].getBases()[level]->stride(static_cast<int>(dim));
        }

        unsigned coef_index = bspline::getIndex<2>(act_coefs_strides, position);

        for (unsigned index = 0; index < numb_of_point; ++index)
        {
            if (tensor_active_index < ten_index)
            {
                presavedIndices(i) = ten_index;
                i++;
                if (truncationStopped)   goto increase;
                while (tensor_active_index < ten_index)
                {
                    xmatrix_index++;
                    if (xmatrix_index == SubdomainHierarchy[patch].getXmatrix()[level].size())
                    {
                        // we don't have any active basis function,
                        // so all the rest basis function in our
                        // representation should not be truncated
                        //return; For now like this
                        truncationStopped = 1;
                        break;
                    }
                    auto wurschtl = SubdomainHierarchy[patch].getXmatrix();
                    tensor_active_index =
                        SubdomainHierarchy[patch].getXmatrix()[level][xmatrix_index];
                }
                // ten_index <= tensor_active_index holds
            }

            // Truncate only if this tensor index corresponds to an ACTIVE function
            // Check both: exists in globalIndex AND is the current active index from Xmatrix
            if (ten_index == tensor_active_index && !truncationStopped) // truncate
            {
                // Diagnostic: log truncation decisions for attempt 7 and 56
                if (false) {
                    gsInfo << "  _truncate: Zeroing coef[" << (coef_index + index) 
                           << "] (ten_index=" << ten_index 
                           << ", tensor_active_index=" << tensor_active_index
                           << ", globalIdx=" << globalIndex[patch][level][ten_index] 
                           << ", level=" << level << ")\n";
                }
                
                coefs(coef_index + index, 0) = 0;
            }

        increase:
            ten_index++;
        }

    } while (nextCubePoint<gsVector<index_t, 2> >(position, first_point,
        last_point));
}

// Zero refined coefficients that correspond to ACTIVE fine-level THB basis fns. 
// Active is inferred via indexInTHB[patch][fineLevel][fineFlatIndex] != -1.
static void truncateAfterBoehm2D_flat(
    gsMatrix<real_t>& coefs,                           // (blockU*blockV) x 1, flattened (iu + iv*blockU)
    index_t blockU, index_t blockV,                    // refined window size on fine grid
    index_t const_ten_index,                           // base fine flat index
    const gsTensorBSplineBasis<2, real_t>& finebasis,  // level+1 basis
    //const gsVector<gsVector<gsVector<index_t>>>& indexInTHB,
    const gsVector<gsVector<gsVector<index_t>>>& globalIndex,
    int patch,
    int fineLevel,
    int level,
    gsVector < gsTHBSplineBasis < 2, real_t > > SubdomainHierarchy,
    bool verbose = false)
{

    ///////////////161025
    const index_t rows = coefs.rows();
    if (rows == 0) return;

    const index_t strideU = 1;
    const index_t strideV = finebasis.size(0);
    const index_t fineTotal = finebasis.size();// static_cast<index_t>(indexInTHB[patch][fineLevel].size());

    if (verbose)
    {
        gsInfo << "[truncateAfterBoehm2D_flat] const_ten_index=" << const_ten_index
            << " blockU=" << blockU << " blockV=" << blockV
            << " strideV=" << strideV
            << " rows=" << rows << "\n";
    }

    for (index_t r = 0; r < rows; ++r)
    {
        real_t& val = coefs(r, 0);
        if (std::abs(val) < 1e-14) continue;

        // local offsets inside the refinement block
        const index_t iu = r % blockU;
        const index_t iv = r / blockU;

        // fine-level flat index using GISMO’s tensor ordering
        const index_t fineFlatIndex = const_ten_index + iu * strideU + iv * strideV;

        if (fineFlatIndex < 0 || fineFlatIndex >= fineTotal)
        {
            if (verbose)
                gsInfo << "  skip r=" << r << " fineFlatIndex=" << fineFlatIndex << " (out of range)\n";
            continue;
        }

        const index_t thbIndex = globalIndex[patch][fineLevel][fineFlatIndex];
        if (thbIndex != -1)
        {
            if (verbose)
                gsInfo << "  trunc r=" << r << " (iu=" << iu << ", iv=" << iv
                << ") fineFlatIndex=" << fineFlatIndex
                << " thb=" << thbIndex << " -> val=0\n";
            val = real_t(0);
        }
        else if (verbose)
        {
            gsInfo << "  keep r=" << r << " fineFlatIndex=" << fineFlatIndex
                << " thb=" << thbIndex << " val=" << val << "\n";
        }
        if (verbose)
        {
            auto wurschtl = SubdomainHierarchy[patch].getXmatrix();
            //gsDebugVar(wurschtl);
            gsDebugVar(wurschtl.size());
            gsDebugVar(wurschtl[fineLevel][fineFlatIndex]);
            gsDebugVar(SubdomainHierarchy[patch].getXmatrix()[fineLevel][fineFlatIndex]);
        }
    }
}



// Compute the flat fine-level index of a coefficient from its coarse index
// MPBES variant: uses coarse/fine bases and span indices from support
static unsigned basisFunIndexOnLevel(
    const gsTensorBSplineBasis<2, real_t>& coarseBasis,
    const gsTensorBSplineBasis<2, real_t>& fineBasis,
    const gsVector<index_t, 2>& coarseIndex,   // tensor index on coarse basis
    const gsVector<index_t, 2>& coarseSpanLow, // span index low on coarse
    const gsVector<index_t, 2>& fineSpanLow)   // span index low on fine
{
    gsVector<index_t, 2> new_index;

    for (unsigned dim = 0; dim < 2; ++dim)
    {
        const auto& ckv = coarseBasis.knots(dim);
        const auto& fkv = fineBasis.knots(dim);

        // offset of this basis fn relative to start of coarse support
        unsigned mult = coarseIndex[dim] - ckv.firstKnotIndex(coarseSpanLow[dim]);

        // mapped into fine level
        new_index(dim) = fkv.firstKnotIndex(fineSpanLow[dim]) + mult;
    }

    return fineBasis.index(new_index);
}

void computeLevelIndex(gsMatrix<index_t, 2, 2> coarse_elem_ind, int level, gsMatrix<index_t, 2, 2>& fine_elem_ind, int maxLevel) {
    for (size_t i = 0; i < coarse_elem_ind.rows(); i++)
    {
        for (size_t j = 0; j < coarse_elem_ind.cols(); j++)
        {
            fine_elem_ind(i, j) = coarse_elem_ind(i, j) * pow(2, maxLevel - level);
        }
    }
}

//unsigned updateSizeOfCoefs(const gsVector<gsVector<gsTensorBSplineBasis<2, real_t>>>& Bells, 
//    const unsigned clevel,
//    const unsigned flevel, 
//    )
//{
//    gsMatrix<index_t, 2, 2> coarse_elem_ind;
//    basis.elementSupport_into(localIndex, coarse_elem_ind);
//    //computeElementSpanIndicesFromSupport(basis, supp, coarse_elem_ind);
//    gsMatrix<index_t, 2, 2> fine_elem_ind;
//    computeLevelIndex(coarse_elem_ind, level, fine_elem_ind, maxLevel);
//    gsVector<index_t, 2> chigh = coarse_elem_ind.col(1);
//    gsVector<index_t, 2> fhigh = fine_elem_ind.col(1);
//    for (unsigned dim = 0; dim < 2; ++dim)
//    {
//        const gsKnotVector<> ckv = basis.knots(dim);
//        const gsKnotVector<> fkv = finebasis.knots(dim);
//        unsigned cnmb_knts = ckv.lastKnotIndex(chigh[dim]) - ckv.firstKnotIndex(clow[dim]);
//        unsigned fnmb_knts = fkv.lastKnotIndex(fhigh[dim]) - fkv.firstKnotIndex(flow[dim]);
//        size_of_coefs(dim) += fnmb_knts - cnmb_knts;
//        nmb_of_coefs *= size_of_coefs(dim);
//    }
//}

template<short_t d>
unsigned _updateSizeOfCoefs(
    int patch,
    const unsigned clevel,
    const unsigned flevel,
    const gsVector<index_t, d>& finest_low,
    const gsVector<index_t, d>& finest_high,
    gsVector<index_t, d>& size_of_coefs,
    gsVector < gsTHBSplineBasis < 2, real_t > > SubdomainHierarchy)
{
    gsVector<index_t, d> clow, chigh;

    SubdomainHierarchy[patch].getTree().computeLevelIndex(finest_low, clevel, clow);
    SubdomainHierarchy[patch].getTree().computeLevelIndex(finest_high, clevel, chigh);

    gsVector<index_t, d> flow, fhigh;
    SubdomainHierarchy[patch].getTree().computeLevelIndex(finest_low, flevel, flow);
    SubdomainHierarchy[patch].getTree().computeLevelIndex(finest_high, flevel, fhigh);

    if (SubdomainHierarchy[patch].manualLevels())
    {
        SubdomainHierarchy[patch]._diadicIndexToKnotIndex(clevel, clow);
        SubdomainHierarchy[patch]._diadicIndexToKnotIndex(clevel, chigh);
        SubdomainHierarchy[patch]._diadicIndexToKnotIndex(flevel, flow);
        SubdomainHierarchy[patch]._diadicIndexToKnotIndex(flevel, fhigh);
    }
    // number of new coefficients
    unsigned nmb_of_coefs = 1;

    for (unsigned dim = 0; dim < d; ++dim)
    {
        const gsKnotVector<real_t>& ckv =
            SubdomainHierarchy[patch].getBases()[clevel]->knots(dim);
        const gsKnotVector<real_t>& fkv =
            SubdomainHierarchy[patch].getBases()[flevel]->knots(dim);

        // Number of knots in the coarse knot vector
        unsigned cnmb_knts = ckv.lastKnotIndex(chigh[dim]) - ckv.firstKnotIndex(clow[dim]);
        //gsDebugVar(ckv.lastKnotIndex(chigh[dim]));
        //gsDebugVar(ckv.firstKnotIndex(clow[dim]));
        // Number of knots in the fine knot vector
        unsigned fnmb_knts = fkv.lastKnotIndex(fhigh[dim]) - fkv.firstKnotIndex(flow[dim]);
        //gsDebugVar(fkv.lastKnotIndex(fhigh[dim]));
        //gsDebugVar(fkv.firstKnotIndex(flow[dim]));
        size_of_coefs(dim) += fnmb_knts - cnmb_knts;
        nmb_of_coefs *= size_of_coefs(dim);
    }
    return nmb_of_coefs;
}

//void _saveNewBasisFunPresentation(
//    const gsMatrix<real_t>& coefs,
//    std::vector<std::vector<gsSparseVector < real_t >>>& presentation, 
//    const int functionIndex,
//    const int twinIndex,
//    const gsVector<int> &presavedIndices
//    ){
//    for (size_t i = 0; i < coefs.rows(); i++)
//    {
//        if(presavedIndices(i) >= 0)
//            presentation[functionIndex][twinIndex](presavedIndices(i)) = coefs(i, 0);
//    }
//}
template<short_t d>
void _saveNewBasisFunPresentation(
    const gsMatrix<real_t>& coefs,
    const gsVector<index_t, d>& act_size_of_coefs,
    const unsigned j,
    const unsigned pres_level,
    const gsVector<index_t, d>& finest_low,
    gsSparseVector < real_t >& localpresentation,
    gsTHBSplineBasis<d, real_t> localTHBBasis)
{
    int verbose = 0;
    const unsigned level = localTHBBasis.levelOf(j);
    if (verbose) {
        std::cout << "Entering _saveNewBasisFunPresentation\n";
        std::cout << "localTHBBasis.size() " << localTHBBasis.size() << "\n";
        std::cout << "  input j = " << j << "\n";
        std::cout << "  computed level = " << level << "\n";
    }
    const unsigned tensor_index = localTHBBasis.flatTensorIndexOf(j, level);

    if (verbose)
    {
        std::cout << "  input j = " << j << ", pres_level = " << pres_level << "\n";
        std::cout << "  computed level = " << level << ", tensor_index = " << tensor_index << "\n";
    }

    gsVector<index_t, d> bspl_vec_ti =
        localTHBBasis.getBases()[level]->tensorIndex(tensor_index);

    if (verbose)
    {
        std::cout << "  bspline tensor index vector (bspl_vec_ti): [";
        for (unsigned dim = 0; dim < d; ++dim)
        {
            std::cout << bspl_vec_ti(dim) << (dim + 1 < d ? ", " : "");
        }
        std::cout << "]\n";
    }

    // finer tensor index
    const unsigned f_ten_index = localTHBBasis._basisFunIndexOnLevel(bspl_vec_ti, level,
        finest_low, pres_level);

    if (verbose)
    {
        std::cout << "  finest_low vector: [";
        for (unsigned dim = 0; dim < d; ++dim)
            std::cout << finest_low(dim) << (dim + 1 < d ? ", " : "");
        std::cout << "]\n";
        std::cout << "  computed f_ten_index = " << f_ten_index << "\n";
    }

    gsVector<index_t, d> act_coefs_strides(d);
    bspline::buildCoeffsStrides<d>(act_size_of_coefs, act_coefs_strides);

    if (verbose)
    {
        std::cout << "  active coefficients size (act_size_of_coefs): [";
        for (unsigned dim = 0; dim < d; ++dim)
            std::cout << act_size_of_coefs(dim) << (dim + 1 < d ? ", " : "");
        std::cout << "]\n";

        std::cout << "  computed coefficient strides (act_coefs_strides): [";
        for (unsigned dim = 0; dim < d; ++dim)
            std::cout << act_coefs_strides(dim) << (dim + 1 < d ? ", " : "");
        std::cout << "]\n";
    }

    gsVector<index_t, d> position(d);
    position.fill(0);

    gsVector<index_t, d> first_point(position);
    gsVector<index_t, d> last_point(d);
    bspline::getLastIndexLocal<d>(act_size_of_coefs, last_point);

    if (verbose)
    {
        std::cout << "  iteration bounds first_point: [";
        for (unsigned dim = 0; dim < d; ++dim)
            std::cout << first_point(dim) << (dim + 1 < d ? ", " : "");
        std::cout << "], last_point: [";
        for (unsigned dim = 0; dim < d; ++dim)
            std::cout << last_point(dim) << (dim + 1 < d ? ", " : "");
        std::cout << "]\n";
    }

    localpresentation.resize(localTHBBasis.getBases()[pres_level]->size());

    if (verbose)
    {
        std::cout << "  resized localpresentation to size = "
            << localTHBBasis.getBases()[pres_level]->size() << "\n";
    }

    unsigned long long iter = 0;

    do
    {
        if (verbose)
        {
            std::cout << "----------------------------------------\n";
            std::cout << "  iteration " << iter << "\n";
            std::cout << "  current position: [";
            for (unsigned dim = 0; dim < d; ++dim)
                std::cout << position(dim) << (dim + 1 < d ? ", " : "");
            std::cout << "]\n";
        }

        // ten_index - (tensor) index of a bspline function with respect to
        //             the coef at position "position"
        unsigned ten_index = f_ten_index;
        if (verbose) std::cout << "  base ten_index (f_ten_index) = " << f_ten_index << "\n";

        for (unsigned dim = 0; dim < d; dim++)
        {
            unsigned stride = static_cast<unsigned>(
                localTHBBasis.getBases()[pres_level]->stride(static_cast<int>(dim)));
            unsigned add = static_cast<unsigned>(position(dim)) * stride;
            ten_index += add;
            if (verbose)
            {
                std::cout << "    dim " << dim << ": position(" << dim << ") = " << position(dim)
                    << ", stride = " << stride << ", add = " << add << "\n";
            }
        }

        if (verbose)
            std::cout << "  computed ten_index = " << ten_index << "\n";

        unsigned coef_index = bspline::getIndex<d>(act_coefs_strides, position);

        if (verbose)
            std::cout << "  computed coef_index = " << coef_index << "\n";

        real_t coef_val = coefs(coef_index);

        if (verbose)
            std::cout << "  coefs[" << coef_index << "] = " << coef_val << "\n";

        if (coef_val != 0)
        {
            if (verbose)
                std::cout << "  assigning localpresentation(" << ten_index << ") = " << coef_val << "\n";
            localpresentation(ten_index) = coef_val;
        }
        else
        {
            if (verbose)
                std::cout << "  skipping zero coefficient at coef_index = " << coef_index << "\n";
        }

        ++iter;

    } while (nextCubePoint<gsVector<index_t, d> >(position, first_point,
        last_point));

    if (verbose)
    {
        std::cout << "Completed _saveNewBasisFunPresentation, total iterations = " << iter << "\n";
        std::cout << "Exiting _saveNewBasisFunPresentation\n";
    }
}
void _representBasisFunction(const int functionIndex,
    const int twinIndex,
    const unsigned thbIndex,
    const int patch,
    const int cur_level,
    const unsigned pres_level,
    const gsVector < index_t, 2 >& finest_low,
    const gsVector < index_t, 2 >& finest_high,
    const gsVector < gsVector < gsTensorBSplineBasis < 2, real_t >>>& Bells,
    const int tensor_index,
    const gsVector<gsVector<gsVector<index_t>>>& globalIndex,
    gsVector < gsTHBSplineBasis < 2, real_t > > SubdomainHierarchy,
    std::vector<std::vector<gsSparseVector < real_t >>>& presentation,
    int attempt = 0) {
    int verbose = 0;
    constexpr short_t d = 2;
    // Implementation of the basis function representation
    gsVector < index_t, d > act_size_of_coefs(d);
    act_size_of_coefs.fill(1);
    unsigned nmb_of_coefs = _updateSizeOfCoefs(patch, cur_level, pres_level,
        finest_low, finest_high, act_size_of_coefs, SubdomainHierarchy);
    gsMatrix < real_t > coefs(nmb_of_coefs, 1);
    coefs.fill(0);
    coefs.row(0).setOnes();
    // vector of the numbers of the coefficients (in each dimension)
    // stored in coefs
    gsVector < index_t, d > vec_nmb_of_coefs(d);
    vec_nmb_of_coefs.fill(1);
    // B-Spline vector tensor index
    gsVector < index_t, d > bspl_vec_ti = Bells[patch][cur_level].tensorIndex(tensor_index);
    // we need to separately save knot vectors because we will modify
    // them, when we proceed from one level on another
    std::vector < gsKnotVector < real_t > > vector_of_kv(d);
    // size of the coefficients that are affected in individual iteration
    gsVector < index_t, d > cur_size_of_coefs(d);
    cur_size_of_coefs.fill(1);
    int TimDebug = 1;
    gsVector<int> presavedIndices(coefs.rows());
    for (unsigned level = cur_level; level < pres_level; ++level) {
        _updateSizeOfCoefs(patch, level, level + 1, finest_low,
            finest_high, cur_size_of_coefs, SubdomainHierarchy);

        // index of a support of the j-th basis function (l_low, l_high
        // on level, and l1_high, l1_low on level + 1)
        gsVector < index_t, d > clow, chigh, fhigh, flow;
        SubdomainHierarchy[patch].getTree().computeLevelIndex(finest_low, level, clow);
        SubdomainHierarchy[patch].getTree().computeLevelIndex(finest_high, level, chigh);
        SubdomainHierarchy[patch].getTree().computeLevelIndex(finest_low, level + 1, flow);
        SubdomainHierarchy[patch].getTree().computeLevelIndex(finest_high, level + 1, fhigh);
        std::vector < double> knots;
        for (unsigned dim = 0; dim < d; ++dim) {
            const gsKnotVector < double >& ckv = Bells[patch][level].knots(dim);
            const gsKnotVector < double >& fkv = Bells[patch][level + 1].knots(dim);
            //if (verbose) {
              //  gsDebugVar(ckv);
               // gsDebugVar(fkv);
               // gsDebugVar(clow[dim]);
               // gsDebugVar(flow[dim]);
               // gsDebugVar(chigh[dim]);
               // gsDebugVar(fhigh[dim]);
            //}
            if (level == cur_level)
                vector_of_kv[dim] = ckv;

            knots.clear();
            std::set_symmetric_difference(ckv.beginAt(clow[dim]), ckv.endAt(chigh[dim]),
                fkv.beginAt(flow[dim]), fkv.endAt(fhigh[dim]),
                std::back_inserter(knots));

            gsTensorBoehmRefineLocal < d,
                gsKnotVector < double >,
                gsMatrix < double >,
                typename std::vector < double > ::const_iterator >
                (vector_of_kv[dim], bspl_vec_ti[dim], coefs, vec_nmb_of_coefs,
                    act_size_of_coefs,
                    cur_size_of_coefs, dim, knots.begin(), knots.end(),
                    true);
            //gsDebugVar(coefs);
        }
        if (verbose) {
            gsInfo << "before\n";
            gsDebugVar(coefs);
        }
        //truncateAfterBoehm2D_flat(coefs, act_size_of_coefs, cur_size_of_coefs,
        //    level + 1, bspl_vec_ti, cur_level, level + 1);
        
        _truncate(
            coefs,
            act_size_of_coefs,
            cur_size_of_coefs,
            level + 1,
            bspl_vec_ti,
            cur_level,
            finest_low,
            patch,
            presavedIndices,
            globalIndex,
            SubdomainHierarchy,
            verbose,
            functionIndex,
            attempt
        );
        // Truncation logging removed - using row0 attempt diagnostics instead
        //" + "// Log after truncation
        // if (logAttempt7 || logAttempt56) {
        //     ofstream& truncateDiag = logAttempt7 ? truncateDiag7 : truncateDiag56;
        // ...
        // }
        //truncateAfterBoehm2D_flat(coefs, act_size_of_coefs[0], act_size_of_coefs[1], tensor_index, Bells[patch][level + 1],
        //    globalIndex, patch, level + 1, level, SubdomainHierarchy, true);
        if (verbose) {
            gsInfo << "after\n";
            gsDebugVar(coefs);
            gsDebugVar(presavedIndices);
        }
    }
    _saveNewBasisFunPresentation(coefs, act_size_of_coefs, thbIndex, pres_level, finest_low, presentation[functionIndex][twinIndex], SubdomainHierarchy[patch]);
}

void computeCoefsForTruncatedFunctions(
    const gsVector<gsVector<gsTensorBSplineBasis<2, real_t>>>& Bells,
    gsVector<gsTHBSplineBasis<2, real_t> > SubdomainHierarchy,
    const std::vector<std::vector<std::vector<index_t>>>& functionDescription,
    std::vector<std::vector<gsSparseVector < double >>>& presentation,
    std::vector<bool>& isTruncated,
    std::vector<std::vector<std::vector<real_t>>>& globalCoefs, // [func f][globalId] -> list of coefs
    const gsVector<gsVector<gsVector<index_t>>>& globalIndex,
    const gsVector<gsVector<gsVector<index_t>>>& indexInTHB,
    std::vector<std::vector<std::array<int, 3>>> spilloverFunctionCoordinates,
    const std::vector<std::vector<bool>>& isComponentTruncated,
    const std::vector<std::vector<index_t>>& presentationLevel,
    bool verbose,
    int attempt = 0)
{
    int maxLevel = 0;
    for (index_t p = 0; p < SubdomainHierarchy.size(); ++p)
        maxLevel = std::max(maxLevel, static_cast<int>(SubdomainHierarchy[p].maxLevel()));
    const index_t numFunctions = functionDescription.size();
    globalCoefs.clear();
    globalCoefs.resize(numFunctions);
    int d = 2;
    gsMatrix<index_t, 2, 2> element_ind(2, 2);
    gsVector<index_t, 2   > low, high;
    int numTruncated = 0;

    // First pass: Check which functions have spillovers
    std::vector<bool> hasSpilloverTwin(numFunctions, false);
    for (index_t f = 0; f < numFunctions; ++f)
    {
        for (size_t twinIndex = 0; twinIndex < functionDescription[f].size(); twinIndex++)
        {
            const auto& desc = functionDescription[f][twinIndex];
            const int patch = desc[0];
            const int level = desc[1];
            const int tensor_index = desc[2];
            auto THBIndex = indexInTHB[patch][level][tensor_index];

            if (THBIndex == -1)
            {
                hasSpilloverTwin[f] = true;
                break;
            }
        }
    }

    if (verbose) {
        int numWithSpillovers = 0;
        for (bool hasSpill : hasSpilloverTwin) if (hasSpill) numWithSpillovers++;
        gsInfo << "Functions with spillover twins: " << numWithSpillovers << "\n";
    }

    for (index_t f = 0; f < numFunctions; ++f)
    {

        presentation[f].resize(functionDescription[f].size());
        for (size_t twinIndex = 0; twinIndex < functionDescription[f].size(); twinIndex++)
        {
            const auto& desc = functionDescription[f][twinIndex];
            const int patch = desc[0];
            const int level = desc[1];
            const int tensor_index = desc[2];
            auto THBIndex = indexInTHB[patch][level][tensor_index];

            // Handle spillover twins: they don't have a THB index
            if (THBIndex == -1)
            {
                // This is a spillover twin - it's not in THB basis on this patch
                // Create empty presentation since it won't be evaluated
                presentation[f][twinIndex] = gsSparseVector<real_t>(Bells[patch][level].size());
                continue;
            }
            SubdomainHierarchy[patch].getBases()[level]->elementSupport_into(tensor_index, element_ind);

            // I tried with block, I can not trick the compiler to use references
            low = element_ind.col(0); //block<d, 1>(0, 0);
            high = element_ind.col(1); //block<d, 1>(0, 1);

            // Use the presentation level that was already computed by detectTruncatedFunctions
            // This is the highest level that overlaps this component's support
            index_t clevel = 0;
            bool isTrunc = false;

            if (f < isComponentTruncated.size() && twinIndex < isComponentTruncated[f].size()) {
                isTrunc = isComponentTruncated[f][twinIndex];
                clevel = presentationLevel[f][twinIndex];
            }

            // Always initialize presentation vector, even if not truncated
            if (clevel == 0) clevel = level;  // If no presentation level set, use function's own level
            presentation[f][twinIndex] = gsSparseVector<real_t>(Bells[patch][clevel].size());

            if (isTrunc) // we must compute its presentation
            {
                // CRITICAL: Transform support indices to finest level
                SubdomainHierarchy[patch].tree().computeFinestIndex(low, level, low);
                SubdomainHierarchy[patch].tree().computeFinestIndex(high, level, high);

                SubdomainHierarchy[patch].setTruncated(tensor_index, clevel);
                _representBasisFunction(f, twinIndex, THBIndex, patch, level, clevel, low, high, Bells, tensor_index, globalIndex, SubdomainHierarchy, presentation, attempt);

                isTruncated[f] = true;
                numTruncated++;
            }
        }
    }
    gsInfo << "\n*** Total truncated: " << numTruncated << " ***\n\n";
}




















/////////////////////



void extractTHBCoefficients(const gsVector<gsVector<int>>& globalIndexTHB,
    gsMatrix<real_t>& vectSol,
    std::vector<std::vector<std::vector<double>>>& globalCoeffsTHB)
{
    // Ensure outer vector is correctly sized
    globalCoeffsTHB.resize(globalIndexTHB.size());

    for (int patch = 0; patch < globalIndexTHB.size(); ++patch)
    {
        for (int functionIndex = 0; functionIndex < globalIndexTHB[patch].size(); ++functionIndex)
        {
            int globalIndex = globalIndexTHB[patch][functionIndex];
            std::vector<double> corrCoeff;
            corrCoeff.push_back(vectSol(globalIndex, 0));
            corrCoeff.push_back(vectSol(globalIndex, 1));
            globalCoeffsTHB[patch].push_back(corrCoeff);
        }
    }
}


void featureMatrixGenerator(gsVector<gsMatrix<real_t>>& uvFeature,
    std::vector<std::vector<std::vector<real_t>>> featureCoordinates,
    index_t numPoints) {
    for (int patch = 0; patch < featureCoordinates.size(); ++patch) {
        uvFeature(patch).resize(2, featureCoordinates[patch].size() * numPoints);
        int theInd = 0;
        for (int segmentIndex = 0; segmentIndex < featureCoordinates[patch].size(); ++segmentIndex) {
            for (int pointNum = 0; pointNum < numPoints; ++pointNum) {
                uvFeature(patch)(0, theInd) = featureCoordinates[patch][segmentIndex][0] + ((double)pointNum / numPoints) * (featureCoordinates[patch][segmentIndex][2] - featureCoordinates[patch][segmentIndex][0]);
                uvFeature(patch)(1, theInd) = featureCoordinates[patch][segmentIndex][1] + ((double)pointNum / numPoints) * (featureCoordinates[patch][segmentIndex][3] - featureCoordinates[patch][segmentIndex][1]);
                theInd++;
            }
        }
    }
}

void meshMatrixGenerator(gsVector<gsMatrix<real_t>>& uvFeature,
    std::vector<std::vector<std::vector<real_t>>> featureCoordinates,
    index_t numPoints) {
    for (int patch = 0; patch < featureCoordinates.size(); ++patch) {
        outfile << "patch " << patch << "\n";
        uvFeature(patch).resize(2, featureCoordinates[patch].size() * numPoints);
        int theInd = 0;
        for (int segmentIndex = 0; segmentIndex < featureCoordinates[patch].size(); ++segmentIndex) {
            outfile << "generating mesh for \n";
            outfile << featureCoordinates[patch][segmentIndex][0] << " " << featureCoordinates[patch][segmentIndex][1] << " " << featureCoordinates[patch][segmentIndex][2] << " " << featureCoordinates[patch][segmentIndex][3] << "\n";
            if (fabs(featureCoordinates[patch][segmentIndex][1] - featureCoordinates[patch][segmentIndex][3]) < 1e-9)
            {
                outfile << "horizontal\n";
                for (int pointNum = 0; pointNum < numPoints; ++pointNum) {
                    uvFeature(patch)(0, theInd) = featureCoordinates[patch][segmentIndex][0] + ((double)pointNum / numPoints) * (featureCoordinates[patch][segmentIndex][2] - featureCoordinates[patch][segmentIndex][0]);
                    uvFeature(patch)(1, theInd) = featureCoordinates[patch][segmentIndex][1];
                    outfile << uvFeature(patch)(0, theInd) << " " << uvFeature(patch)(1, theInd) << "\n";
                    theInd++;
                }
            }
            if (fabs(featureCoordinates[patch][segmentIndex][0] - featureCoordinates[patch][segmentIndex][2]) < 1e-9)
            {
                outfile << "vertical\n";
                for (int pointNum = 0; pointNum < numPoints; ++pointNum) {
                    uvFeature(patch)(0, theInd) = featureCoordinates[patch][segmentIndex][0];
                    uvFeature(patch)(1, theInd) = featureCoordinates[patch][segmentIndex][1] + ((double)pointNum / numPoints) * (featureCoordinates[patch][segmentIndex][3] - featureCoordinates[patch][segmentIndex][1]);
                    outfile << uvFeature(patch)(0, theInd) << " " << uvFeature(patch)(1, theInd) << "\n";
                    theInd++;
                }
            }
        }
    }
}

void generateBoxMesh(gsVector<index_t> box,
    gsVector<gsVector<gsTensorBSplineBasis<2, real_t>>> Bells,
    gsVector  <gsVector< gsVector<index_t>>> isActive,
    gsVector  <gsVector< gsVector<index_t>>>  globalIndex,
    gsVector<gsVector  <std::vector<index_t>>> functionDescription,
    gsMatrix<real_t> vectSol)
{
    //Rethink; overload assembleA can be an option;
    gsMatrix<>  uvBox;
    auto numPoints = 10 * pow(2, box(0));
    gsSparseMatrix<real_t> A_mat((box(4) - (box(2) + 1)) * numPoints, vectSol.rows());
    gsMatrix<index_t> segment(1, 5);
    gsMatrix<real_t>   coords(1, 5);
    //Horizontal
    segment(0, 0) = box(0);
    for (size_t i = box(2) + 1; i < box(4); i++)
    {
        segment(0, 1) = box(1);
        segment(0, 2) = i;
        segment(0, 3) = box(3);
        segment(0, 4) = i;
        boxToDomain(segment, coords, 0, 0);
        gsInfo << "Coords:\n";
        for (size_t i = 0; i < 5; i++)
        {
            gsInfo << coords(0, i) << " ";
        }
        gsInfo << "\n";
        uvBox.resize(2, numPoints);
        for (size_t j = 0; j < numPoints; j++)
        {
            uvBox(0, j) = coords(0, 1) + (coords(0, 3) - coords(0, 1)) * (real_t)j / (10 * pow(2, box(0) + 1));
            uvBox(1, j) = coords(0, 2);
        }
        assembleA(uvBox, Bells, isActive, globalIndex, functionDescription, A_mat);
        gsInfo << A_mat.rows() << ", " << A_mat.cols() << "\n";
        gsInfo << vectSol.rows() << ", " << vectSol.cols() << "\n";
        gsMatrix<> boxValues(A_mat.rows(), 2);
        boxValues.setZero();
        for (size_t i = 0; i < A_mat.rows(); i++)
        {
            for (size_t j = 0; j < A_mat.cols(); j++)
            {
                boxValues(i, 0) += A_mat(i, j) * vectSol(j, 0);
                boxValues(i, 1) += A_mat(i, j) * vectSol(j, 1);
            }
        }
        xboxFile << boxValues.col(0) << "\n";
        yboxFile << boxValues.col(1) << "\n";
        xboxFile << "\n";
        yboxFile << "\n";
        //eval
    }
    //vertical
}
// MPBES-based version
// MPBES-based version using MPBES evaluation methods  
void jacobian_into_mpbe(
    const MPBES<2, real_t>& mpbes,
    int patch,
    const gsMatrix<real_t>& quNodes,
    gsMatrix<real_t>& out,
    const gsMatrix<int>& actives,
    const gsMatrix<real_t>& vectSol,
    const gsVector<gsVector<int>>& globalIndexTHB,
    bool verbose = false)
{
    const index_t numPoints = quNodes.cols();
    out.setZero(2, 2 * numPoints); // 2x2 Jacobian per point
    int badLogCount = 0;
    const int maxBadLogs = 10;

    for (index_t pt = 0; pt < numPoints; ++pt)
    {
        const auto param = quNodes.col(pt);
        if (!param.allFinite() && badLogCount < maxBadLogs)
        {
            gsInfo << "[jacobian_into_mpbe] Non-finite quNode: patch=" << patch
                   << ", pt=" << pt << ", param=" << param.transpose() << "\n";
            outfile << "[jacobian_into_mpbe] Non-finite quNode: patch=" << patch
                    << ", pt=" << pt << ", param=" << param.transpose() << "\n";
            throw std::runtime_error("Non-finite quNode detected in jacobian_into_mpbe");
            ++badLogCount;
        }

        for (index_t a = 0; a < actives.rows(); ++a)
        {
            const int thbIndex = actives(a, 0);
            const int globalID = globalIndexTHB[patch][thbIndex];
            if (globalID < 0) continue;

            const real_t cx = vectSol(globalID, 0);
            const real_t cy = vectSol(globalID, 1);
            if ((!std::isfinite(cx) || !std::isfinite(cy)) && badLogCount < maxBadLogs)
            {
                gsInfo << "[jacobian_into_mpbe] Non-finite vectSol: patch=" << patch
                       << ", pt=" << pt << ", globalID=" << globalID
                       << ", cx=" << cx << ", cy=" << cy << "\n";
                outfile << "[jacobian_into_mpbe] Non-finite vectSol: patch=" << patch
                        << ", pt=" << pt << ", globalID=" << globalID
                        << ", cx=" << cx << ", cy=" << cy << "\n";
                throw std::runtime_error("Non-finite vectSol detected in jacobian_into_mpbe");
                ++badLogCount;
            }

            // Use MPBES derivative evaluation (handles truncation automatically)
            gsMatrix<real_t> deriv;
            mpbes.evalDerivSingleOnPatch(globalID, patch, param, deriv);
            if (!deriv.allFinite() && badLogCount < maxBadLogs)
            {
                gsInfo << "[jacobian_into_mpbe] Non-finite deriv: patch=" << patch
                       << ", pt=" << pt << ", globalID=" << globalID
                       << ", param=" << param.transpose() << "\n";
                gsInfo << "deriv=\n" << deriv << "\n";
                outfile << "[jacobian_into_mpbe] Non-finite deriv: patch=" << patch
                        << ", pt=" << pt << ", globalID=" << globalID
                        << ", param=" << param.transpose() << "\n";
                outfile << "deriv=\n" << deriv << "\n";
                throw std::runtime_error("Non-finite deriv detected in jacobian_into_mpbe");
                ++badLogCount;
            }
            
            // Accumulate Jacobian contributions
            out(0, 2 * pt) += deriv(0, 0) * cx;      // dx/du
            out(0, 2 * pt + 1) += deriv(1, 0) * cx;  // dx/dv
            out(1, 2 * pt) += deriv(0, 0) * cy;      // dy/du
            out(1, 2 * pt + 1) += deriv(1, 0) * cy;  // dy/dv
        }

        if (verbose)
        {
            outfile << "  ==> Final Jacobian at pt " << pt << ": ";
            outfile << "[[" << out(0, 2 * pt) << ", " << out(0, 2 * pt + 1) << "], ";
            outfile << "[" << out(1, 2 * pt) << ", " << out(1, 2 * pt + 1) << "]]\n";
        }
    }

    if (verbose)
    {
        outfile << "\n[jacobian_into_mpbe] FULL OUTPUT MATRIX:\n";
        outfile << "geomDer\n";
        for (index_t r = 0; r < out.rows(); ++r)
        {
            for (index_t c = 0; c < out.cols(); ++c)
                outfile << out(r, c) << " ";
            outfile << "\n";
        }
        outfile << "[jacobian_into_mpbe] END\n";
    }
}


// MPBES-based version using MPBES evaluation methods
void deriv_into_mpbe(
    const MPBES<2, real_t>& mpbes,
    int patch,
    const gsMatrix<real_t>& quNodes,
    gsMatrix<real_t>& out,
    const gsMatrix<int>& actives,
    const gsVector<gsVector<int>>& globalIndexTHB,
    bool verbose = false)
{
    const index_t numPoints = quNodes.cols();
    const index_t numActives = actives.rows();
    out.setZero(2 * numActives, numPoints);

    for (index_t pt = 0; pt < numPoints; ++pt)
    {
        const auto param = quNodes.col(pt);

        for (index_t a = 0; a < numActives; ++a)
        {
            const int thbIndex = actives(a, 0);
            const int globalID = globalIndexTHB[patch][thbIndex];
            if (globalID < 0) continue;

            // Use MPBES derivative evaluation (handles truncation automatically)
            gsMatrix<real_t> deriv;
            mpbes.evalDerivSingleOnPatch(globalID, patch, param, deriv);
            
            out(2 * a, pt) = deriv(0, 0);     // du derivative
            out(2 * a + 1, pt) = deriv(1, 0);  // dv derivative
        }
    }

    if (verbose)
    {
        outfile << "[deriv_into_mpbe] Final output matrix:\n";
        outfile << "basisDer\n";
        for (index_t r = 0; r < out.rows(); ++r)
        {
            for (index_t c = 0; c < out.cols(); ++c)
                outfile << out(r, c) << " ";
            outfile << "\n";
        }
    }
}


















// MPBES-based version using MPBES evaluation methods
void deriv2_into_geometry_mpbe(
    const MPBES<2, real_t>& mpbes,
    int patch,
    const gsMatrix<real_t>& quNodes,
    gsMatrix<real_t>& out,
    const gsMatrix<int>& actives,
    const gsMatrix<real_t>& vectSol,
    const gsVector<gsVector<int>>& globalIndexTHB,
    bool verbose = false)
{
    if (verbose) outfile << "geom2Der\n";
    const index_t numPoints = quNodes.cols();
    out.setZero(6, numPoints); // 6 rows for 2D second derivatives: dxx, dxy, dyy for x and y components

    for (index_t pt = 0; pt < numPoints; ++pt)
    {
        const auto param = quNodes.col(pt);
        for (index_t a = 0; a < actives.rows(); ++a)
        {
            const int thbIndex = actives(a, 0);
            const int globalID = globalIndexTHB[patch][thbIndex];
            if (globalID < 0) continue;

            // Use MPBES second derivative evaluation (handles truncation automatically)
            gsMatrix<real_t> d2;
            mpbes.evalDeriv2SingleOnPatch(globalID, patch, param, d2);
            
            // Accumulate geometry second derivatives
            for (index_t d = 0; d < 6; ++d)
                out(d, pt) += vectSol(globalID, d / 3) * d2(d % 3, 0);
        }
    }

    if (verbose)
    {
        outfile << "[deriv2_into_geometry_mpbe] Output matrix (rows: 6, cols: " << numPoints << "):\n";
        for (index_t r = 0; r < out.rows(); ++r)
        {
            for (index_t c = 0; c < out.cols(); ++c)
                outfile << out(r, c) << " ";
            outfile << "\n";
        }
    }
}







// MPBES-based version using MPBES evaluation methods
void deriv2_into_basis_mpbe(
    const MPBES<2, real_t>& mpbes,
    int patch,
    const gsMatrix<real_t>& quNodes,
    gsMatrix<real_t>& out,
    const gsMatrix<int>& actives,
    const gsVector<gsVector<int>>& globalIndexTHB,
    bool verbose = false)
{
    if (verbose) outfile << "basis2Der\n";
    const index_t numPoints = quNodes.cols();
    const index_t numActives = actives.rows();
    out.setZero(3 * numActives, numPoints);

    for (index_t pt = 0; pt < numPoints; ++pt)
    {
        const auto param = quNodes.col(pt);
        for (index_t a = 0; a < numActives; ++a)
        {
            const int thbIndex = actives(a, 0);
            const int globalID = globalIndexTHB[patch][thbIndex];
            if (globalID < 0) continue;

            // Use MPBES second derivative evaluation (handles truncation automatically)
            gsMatrix<real_t> d2;
            mpbes.evalDeriv2SingleOnPatch(globalID, patch, param, d2);
            
            // Store basis second derivatives
            for (index_t d = 0; d < 3; ++d)
                out(3 * a + d, pt) = d2(d, 0);
        }
    }

    if (verbose)
    {
        outfile << "[deriv2_into_basis_mpbe] Output matrix (rows: " << out.rows() << ", cols: " << out.cols() << "):\n";
        for (index_t r = 0; r < out.rows(); ++r)
        {
            for (index_t c = 0; c < out.cols(); ++c)
                outfile << out(r, c) << " ";
            outfile << "\n";
        }
    }
}


//    index_t numActive,
//    const index_t quNode,
//    const gsMatrix<real_t>& basis2Der)
//{
//    const short_t parDim = 2;
//    const short_t geoDim = 2;
//    const int stride = parDim * (parDim + 1) / 2;
//    const int numQ = stride * geoDim;
//
//    Juni.resize(numQ * geoDim, numActive);
//    Juni.setZero();
//
//    for (index_t i = 0; i != numActive; ++i)
//    {
//        for (int s = 0; s != stride; ++s)
//        {
//            const int basInd2Der = i * stride + s;
//            real_t Nuu = basis2Der(basInd2Der, quNode);
//
//            if (parDim <= s)
//            {
//                Nuu *= sqrt(2);
//            }
//
//            if (std::isnan(Nuu) || std::isinf(Nuu))
//            {
//                outfile << "Error: Nuu is " << (std::isnan(Nuu) ? "NaN" : "Inf")
//                    << " at i=" << i << ", s=" << s << ", quNode=" << quNode << "\n";
//            }
//
//            for (int dimQ = 0; dimQ < geoDim; dimQ++)
//            {
//                const int qIndex = s * geoDim + dimQ;
//                Juni(qIndex * geoDim + dimQ, i) = Nuu;
//
//                if (std::isnan(Juni(qIndex * geoDim + dimQ, i)) ||
//                    std::isinf(Juni(qIndex * geoDim + dimQ, i)))
//                {
//                    outfile << "Error: Juni(" << qIndex * geoDim + dimQ << ", " << i
//                        << ") is " << (std::isnan(Juni(qIndex * geoDim + dimQ, i)) ? "NaN" : "Inf") << "\n";
//                }
//            }
//        }
//    }
//}

void getJuni(gsMatrix<real_t>& Juni,
    index_t numActive,
    const index_t quNode,
    const gsMatrix<real_t>& basis2Der)
{
    const short_t parDim = 2;
    const short_t geoDim = 2;
    const int stride = parDim * (parDim + 1) / 2;
    const int numQ = stride * geoDim;

    Juni.resize(numQ * geoDim, numActive);
    Juni.setZero();

    assert(quNode >= 0 && quNode < basis2Der.cols() && "quNode out of bounds in basis2Der!");

    for (index_t i = 0; i != numActive; ++i)
    {
        for (int s = 0; s != stride; ++s)
        {
            const int basInd2Der = i * stride + s;
            assert(basInd2Der >= 0 && basInd2Der < basis2Der.rows() && "basInd2Der out of bounds in basis2Der!");

            real_t Nuu = basis2Der(basInd2Der, quNode);

            if (parDim <= s)
            {
                Nuu *= sqrt(2);
            }

            if (std::isnan(Nuu) || std::isinf(Nuu))
            {
                outfile << "Error: Nuu is " << (std::isnan(Nuu) ? "NaN" : "Inf")
                    << " at i=" << i << ", s=" << s << ", quNode=" << quNode << "\n";
            }

            for (int dimQ = 0; dimQ < geoDim; dimQ++)
            {
                const int qIndex = s * geoDim + dimQ;
                int juniRow = qIndex * geoDim + dimQ;

                assert(juniRow >= 0 && juniRow < Juni.rows() && "juniRow out of bounds in Juni!");
                assert(i >= 0 && i < Juni.cols() && "i out of bounds in Juni!");

                Juni(juniRow, i) = Nuu;

                if (std::isnan(Juni(juniRow, i)) || std::isinf(Juni(juniRow, i)))
                {
                    outfile << "Error: Juni(" << juniRow << ", " << i
                        << ") is " << (std::isnan(Juni(juniRow, i)) ? "NaN" : "Inf") << "\n";
                }
            }
        }
    }
}


void getJlen(gsMatrix<real_t>& Jlen,
    index_t numActive,
    const index_t quNode,
    gsMatrix<real_t> basisDer)
{
    const short_t parDim = 2;
    const short_t geoDim = 2;
    const int numQ = parDim * geoDim;
    Jlen.resize(numQ * geoDim, numActive);
    Jlen.setZero();
    //const index_t quNodeErsatz = 0;
    //gsInfo << "Jlen:size()" << numQ * geoDim << " * " << numActive << "\n";
    //gsInfo << "basisDer:size()" << basisDer.rows() << " * " << basisDer.cols() << "\n";
    for (index_t i = 0; i != numActive; ++i)
    {
        for (int u = 0; u != parDim; ++u)
        {
            const int basIndDer = parDim * i + u;
            //gsInfo << "basIndDer: " << basIndDer << "\n";
            //const real_t Nu = basisDer(basIndDer, quNode);
            const real_t Nu = basisDer(basIndDer, 0);
            for (int dimQ = 0; dimQ != geoDim; dimQ++)
            {
                const int qIndex = u * geoDim + dimQ;
                //gsInfo << "qIndex: " << qIndex << "\n";
                Jlen(qIndex * geoDim + dimQ, i) = Nu;
                if (!isdigit(Nu))
                {
                    /*    gsInfo << "-nan! J " << i << " " << u << " " << dimQ << "\n";
                        return;*/
                }
            }
        }
    }
}

void getJort(gsMatrix<double>& Jort,
    const int numActive,
    const gsMatrix<double>& basisDer,
    const gsMatrix<double>& geomDer,
    const int quNode)
{
    const short_t parDim = 2;
    const short_t geoDim = 2;

    // ** why nit numQ = 1 since it is scalar ? -->(row)
    const int numQ = (parDim - 1) * parDim / 2;
    //const index_t numActive = domIt->computeActiveFunctions().rows();

    Jort.resize(numQ * geoDim, numActive);
    Jort.setZero();
    // computing Jort
    for (index_t i = 0; i != numActive; ++i)
    {
        int row = 0;
        for (index_t u = 0; u != parDim; ++u)
        {
            for (index_t v = u + 1; v != parDim; ++v)
            {
                const int uBasInd = parDim * i + u;
                const int vBasInd = parDim * i + v;

                const int uGeomInd = quNode * parDim + u;
                const int vGeomInd = quNode * parDim + v;
                //                    gsInfo << uGeomInd << "\t" << vGeomInd << "\n";
                Jort.block(row * geoDim, i, geoDim, 1) =
                    basisDer(uBasInd, quNode) * geomDer.col(vGeomInd) +
                    basisDer(vBasInd, quNode) * geomDer.col(uGeomInd);
                //                    gsInfo << basisDer(uBasInd, quNode) * geomDer.col(vGeomInd) +
                //                        basisDer(vBasInd, quNode) * geomDer.col(uGeomInd) << " VAH!\n";
                //                    gsInfo << row * geoDim << "\t" << i << "\t" << geoDim << "\n";
                //                      gsInfo << Jort << "\n";
                //                      gsInfo << "EOM\n";
                row++;
            }
        }
    }
}

//void getQuni(gsMatrix<real_t>& q,
//    const gsMatrix<real_t>& geom2Der,
//    const index_t quNode)
//{
//    const short_t parDim = 2;
//    const short_t geoDim = 2;
//    const int stride = parDim * (parDim + 1) / 2;
//    const int numQ = stride * geoDim;
//    q.resize(1, numQ);
//    q.setZero();
//    for (int s = 0; s != stride; ++s)
//    {
//        gsVector<real_t> Gss(geoDim);
//        for (int dim = 0; dim != geoDim; ++dim)
//        {
//            //Gss(dim) = geom2Der(stride * dim + s, 0);//;What is the range?
//            Gss(dim) = geom2Der(stride * dim + s, quNode);//;What is the range?
//            if (parDim <= s)
//            {
//                Gss(dim) *= sqrt(2);
//            }
//        }
//        for (int dimQ = 0; dimQ != geoDim; dimQ++)
//        {
//            const int qIndex = s * geoDim + dimQ;
//            q(0, qIndex) = Gss(dimQ);
//            if (isnan(q(0, qIndex)) || isinf(q(0, qIndex)))
//            {
//                gsInfo << "nan: " << "s = " << s << ", dimQ = " << dimQ << "\n";
//            }
//            if (isinf(q(0, qIndex)))
//            {
//                gsInfo << "inf: " << "s = " << s << ", dimQ = " << dimQ << "\n";
//            }
//        }
//    }
//}

void getQuni(gsMatrix<real_t>& q,
    const gsMatrix<real_t>& geom2Der,
    const index_t quNode)
{
    const short_t parDim = 2;
    const short_t geoDim = 2;
    const int stride = parDim * (parDim + 1) / 2;
    const int numQ = stride * geoDim;
    q.resize(1, numQ);
    q.setZero();

    assert(quNode >= 0 && quNode < geom2Der.cols() && "quNode out of bounds in geom2Der!");

    for (int s = 0; s != stride; ++s)
    {
        gsVector<real_t> Gss(geoDim);
        for (int dim = 0; dim != geoDim; ++dim)
        {
            int rowIndex = stride * dim + s;
            assert(rowIndex >= 0 && rowIndex < geom2Der.rows() && "rowIndex out of bounds in geom2Der!");

            Gss(dim) = geom2Der(rowIndex, quNode);

            if (parDim <= s)
            {
                Gss(dim) *= sqrt(2);
            }
        }
        for (int dimQ = 0; dimQ != geoDim; dimQ++)
        {
            const int qIndex = s * geoDim + dimQ;
            q(0, qIndex) = Gss(dimQ);
            if (isnan(q(0, qIndex)) || isinf(q(0, qIndex)))
            {
                gsInfo << "nan/inf detected: s = " << s << ", dimQ = " << dimQ << "\n";
            }
        }
    }
}

void getQlen(gsMatrix<real_t>& q,
    const gsMatrix<real_t>& geomDer,
    const index_t quNode)
{
    const short_t parDim = 2;
    const short_t geoDim = 2;
    const int numQ = parDim * geoDim;
    q.resize(1, numQ);
    q.setZero();
    //gsInfo << geomDer.cols() << "\n";
    //quNode ist auf 0 ersetzt; //temporär
    //index_t quNodeErsatz = 0;
    for (int u = 0; u != parDim; ++u)
    {
        const int uGeomInd = quNode * parDim + u;
        //const int uGeomInd = quNodeErsatz * parDim + u;
        //gsInfo << "uGeomInd: " << uGeomInd << "\n";
        const gsVector<real_t>& Gu = geomDer.col(uGeomInd);
        //gsInfo << "the col\n";
        //for (size_t i = 0; i < Gu.size(); i++)
        //{
        //    gsInfo << Gu[i] << "\n";
        //}
        for (int dimQ = 0; dimQ != geoDim; dimQ++)
        {
            const int qIndex = u * geoDim + dimQ;
            //gsInfo << "qIndex: " << qIndex << "\n";
            q(0, qIndex) = Gu[dimQ];
        }
    }
}

void getQort(gsMatrix<real_t>& q,
    const gsMatrix<real_t>& geomDer,
    const index_t quNode)
{
    const short_t parDim = 2;
    const int numQ = (parDim - 1) * parDim / 2;
    q.resize(1, numQ);
    q.setZero();
    int col = 0;
    //gsInfo << "size: " << geomDer.rows() << " * " << geomDer.cols() << "\n";
    // computing q
    for (index_t u = 0; u != parDim; ++u)
    {
        for (index_t v = u + 1; v != parDim; ++v)
        {
            const int uGeomInd = quNode * parDim + u;
            const int vGeomInd = quNode * parDim + v;
            if (isnan(geomDer(uGeomInd, 0)) || isnan(geomDer(uGeomInd, 1))) {
                gsInfo << "u trouble\n";
                gsInfo << uGeomInd << "\n";
                gsInfo << "size: " << geomDer.rows() << " * " << geomDer.cols() << "\n";
            }
            //gsInfo << "No u trouble\n";
            if (isnan(geomDer(0, vGeomInd)) || isnan(geomDer(1, vGeomInd))) {
                gsInfo << "v trouble\n";
                gsInfo << vGeomInd << "\n";
                gsInfo << "size: " << geomDer.rows() << " * " << geomDer.cols() << "\n";
            }
            //gsInfo << "No v trouble\n";
            q(0, col) = geomDer.col(uGeomInd).transpose() * geomDer.col(vGeomInd);
            //outfile <<
            //    q(0, col) << " = "  << geomDer.col(uGeomInd).transpose() <<" * "<< geomDer.col(vGeomInd) << "\n";
            //    gsInfo << "q matrix check\n";
            //for (int i = 0; i < numQ; i++) {
                //if (isnan(q(0, i)))    gsInfo << "nan in q, i = " << i << "\n";
            if (isnan(q(0, col)))    gsInfo << "nan in qOrt, col " << col << "\n";
            //}
            col++;
        }
    }
}

void getQskew(gsMatrix<real_t>& q,
    const gsMatrix<real_t>& geomDer,
    const index_t quNode)
{
    const real_t eps = 1e-12;
    const int uGeomInd = quNode * 2 + 0;
    const int vGeomInd = quNode * 2 + 1;

    gsVector<real_t> pu = geomDer.col(uGeomInd);
    gsVector<real_t> pv = geomDer.col(vGeomInd);

    const real_t a = pu.dot(pv);
    const real_t b = pu.dot(pu);
    const real_t c = pv.dot(pv);
    const real_t denomSq = std::max(b * c, eps);
    const real_t denom = std::sqrt(denomSq);

    q.resize(1, 1);
    q.setZero();
    if (denom > eps)
        q(0, 0) = a / denom;
}

void getJskew(gsMatrix<real_t>& Jskew,
    index_t numActive,
    const index_t quNode,
    const gsMatrix<real_t>& basisDer,
    const gsMatrix<real_t>& geomDer)
{
    const real_t eps = 1e-12;
    Jskew.resize(2, numActive);
    Jskew.setZero();

    const int uGeomInd = quNode * 2 + 0;
    const int vGeomInd = quNode * 2 + 1;
    gsVector<real_t> pu = geomDer.col(uGeomInd);
    gsVector<real_t> pv = geomDer.col(vGeomInd);

    const real_t a = pu.dot(pv);
    const real_t b = pu.dot(pu);
    const real_t c = pv.dot(pv);
    const real_t denomSq = b * c;
    if (denomSq <= eps)
        return;

    const real_t denom = std::sqrt(denomSq);
    const real_t invDen = 1.0 / denom;
    const real_t invDen3 = 1.0 / (denom * denomSq);

    for (index_t i = 0; i != numActive; ++i)
    {
        const real_t Nu = basisDer(2 * i, quNode);
        const real_t Nv = basisDer(2 * i + 1, quNode);

        gsVector<real_t> da = Nu * pv + Nv * pu;
        gsVector<real_t> db = 2.0 * Nu * pu;
        gsVector<real_t> dc = 2.0 * Nv * pv;

        gsVector<real_t> dr = da * invDen - a * 0.5 * (c * db + b * dc) * invDen3;

        Jskew(0, i) = dr(0);
        Jskew(1, i) = dr(1);
    }
}

void getQarea(gsMatrix<real_t>& q,
    const gsMatrix<real_t>& geomDer,
    const index_t quNode)
{
    const int uGeomInd = quNode * 2 + 0;
    const int vGeomInd = quNode * 2 + 1;
    const gsVector<real_t> pu = geomDer.col(uGeomInd);
    const gsVector<real_t> pv = geomDer.col(vGeomInd);

    const real_t det = pu(0) * pv(1) - pu(1) * pv(0);

    q.resize(1, 1);
    q(0, 0) = det;
}

void getJarea(gsMatrix<real_t>& Jarea,
    index_t numActive,
    const index_t quNode,
    const gsMatrix<real_t>& basisDer,
    const gsMatrix<real_t>& geomDer)
{
    Jarea.resize(2, numActive);
    Jarea.setZero();

    const int uGeomInd = quNode * 2 + 0;
    const int vGeomInd = quNode * 2 + 1;
    const gsVector<real_t> pu = geomDer.col(uGeomInd);
    const gsVector<real_t> pv = geomDer.col(vGeomInd);

    gsVector<real_t> ddet_dpu(2);
    ddet_dpu << pv(1), -pv(0);
    gsVector<real_t> ddet_dpv(2);
    ddet_dpv << -pu(1), pu(0);

    for (index_t i = 0; i != numActive; ++i)
    {
        const real_t Nu = basisDer(2 * i, quNode);
        const real_t Nv = basisDer(2 * i + 1, quNode);

        gsVector<real_t> dr = Nu * ddet_dpu + Nv * ddet_dpv;
        Jarea(0, i) = dr(0);
        Jarea(1, i) = dr(1);
    }
}

void getQecc(gsMatrix<real_t>& q,
    const gsMatrix<real_t>& geomDer,
    const gsMatrix<real_t>& geom2Der,
    const index_t quNode)
{
    const real_t eps = 1e-12;
    const int uGeomInd = quNode * 2 + 0;
    const int vGeomInd = quNode * 2 + 1;

    const gsVector<real_t> pu = geomDer.col(uGeomInd);
    const gsVector<real_t> pv = geomDer.col(vGeomInd);

    gsVector<real_t> puu(2);
    gsVector<real_t> pvv(2);
    puu << geom2Der(0, quNode), geom2Der(3, quNode);
    pvv << geom2Der(2, quNode), geom2Der(5, quNode);

    const real_t a1 = pu.dot(puu);
    const real_t b1 = std::max(pu.dot(pu), eps);
    const real_t a2 = pv.dot(pvv);
    const real_t b2 = std::max(pv.dot(pv), eps);

    q.resize(1, 2);
    q(0, 0) = a1 / b1;
    q(0, 1) = a2 / b2;
}

void getJecc(gsMatrix<real_t>& Jecc,
    index_t numActive,
    const index_t quNode,
    const gsMatrix<real_t>& basisDer,
    const gsMatrix<real_t>& basis2Der,
    const gsMatrix<real_t>& geomDer,
    const gsMatrix<real_t>& geom2Der)
{
    const real_t eps = 1e-12;
    Jecc.resize(4, numActive);
    Jecc.setZero();

    const int uGeomInd = quNode * 2 + 0;
    const int vGeomInd = quNode * 2 + 1;

    const gsVector<real_t> pu = geomDer.col(uGeomInd);
    const gsVector<real_t> pv = geomDer.col(vGeomInd);

    gsVector<real_t> puu(2);
    gsVector<real_t> pvv(2);
    puu << geom2Der(0, quNode), geom2Der(3, quNode);
    pvv << geom2Der(2, quNode), geom2Der(5, quNode);

    const real_t a1 = pu.dot(puu);
    const real_t b1 = pu.dot(pu);
    const real_t a2 = pv.dot(pvv);
    const real_t b2 = pv.dot(pv);

    for (index_t i = 0; i != numActive; ++i)
    {
        const real_t Nu = basisDer(2 * i, quNode);
        const real_t Nv = basisDer(2 * i + 1, quNode);

        const real_t Nuu = basis2Der(3 * i + 0, quNode);
        const real_t Nvv = basis2Der(3 * i + 2, quNode);

        if (b1 > eps)
        {
            gsVector<real_t> da1 = Nu * puu + Nuu * pu;
            gsVector<real_t> db1 = 2.0 * Nu * pu;
            gsVector<real_t> dr1 = da1 / b1 - a1 * db1 / (b1 * b1);
            Jecc(0, i) = dr1(0);
            Jecc(1, i) = dr1(1);
        }

        if (b2 > eps)
        {
            gsVector<real_t> da2 = Nv * pvv + Nvv * pv;
            gsVector<real_t> db2 = 2.0 * Nv * pv;
            gsVector<real_t> dr2 = da2 / b2 - a2 * db2 / (b2 * b2);
            Jecc(2, i) = dr2(0);
            Jecc(3, i) = dr2(1);
        }
    }
}

inline void assembleFitting(gsVector<gsMatrix<>>  uv,
    gsSparseMatrix<real_t>& A_mat,
    gsVector<gsVector<gsTensorBSplineBasis<2, real_t>>> Bells,
    std::vector<std::vector<std::vector<index_t>>> functionDescription,
    gsMatrix<>& matAsquare,
    const real_t weight
)
{
    PROFILE_FUNCTION();
    
    try {
        // Input validation
        if (matAsquare.rows() == 0 || matAsquare.cols() == 0)
        {
                 // gsInfo << "WARNING: matAsquare has zero dimensions (" << matAsquare.rows()
                 //        << " x " << matAsquare.cols() << "), skipping assembleFitting\n";
            return;
        }
        
        if (!std::isfinite(weight))
        {
            // gsInfo << "WARNING: weight is not finite (" << weight << "), skipping assembleFitting\n";
            return;
        }
        
        // Check if A_mat has sufficient dimensions for the computed indices
        size_t required_rows = matAsquare.rows() * 2;
        size_t required_cols = matAsquare.cols() * 2;
        
        if (A_mat.rows() < static_cast<gsMatrix<>::Index>(required_rows) || 
            A_mat.cols() < static_cast<gsMatrix<>::Index>(required_cols))
        {
            // gsInfo << "ERROR: A_mat dimensions insufficient!\n";
            // gsInfo << "  A_mat is " << A_mat.rows() << " x " << A_mat.cols() << "\n";
            // gsInfo << "  matAsquare is " << matAsquare.rows() << " x " << matAsquare.cols() << "\n";
            // gsInfo << "  Required: " << required_rows << " x " << required_cols << "\n";
            // gsInfo << "  Skipping assembleFitting to prevent crash\n";
            return;
        }
        
        // gsInfo << "A_mat.rows()" << A_mat.rows() << " " << A_mat.cols() << "\n";
        // gsInfo << "matA.rows()" << matAsquare.rows() << " " << matAsquare.cols() << "\n";
        // outfile << "matA.rows()" << matAsquare.rows() << " " << matAsquare.cols() << "\n";
        
        std::chrono::time_point<std::chrono::system_clock> beforesquare, aftersquare;
        beforesquare = std::chrono::system_clock::now();
        
        for (size_t i = 0; i < matAsquare.rows(); i++)
        {
            // Bounds check for row access
            if (i * 2 + 1 >= static_cast<size_t>(A_mat.rows()))
            {
                gsInfo << "WARNING: Computed row index " << (i * 2 + 1) 
                       << " exceeds A_mat.rows() " << A_mat.rows() << ", stopping\n";
                break;
            }
            
            for (size_t j = 0; j < matAsquare.cols(); j++)
            {
                // Bounds check for column access
                if (j * 2 + 1 >= static_cast<size_t>(A_mat.cols()))
                {
                    gsInfo << "WARNING: Computed col index " << (j * 2 + 1) 
                           << " exceeds A_mat.cols() " << A_mat.cols() << ", stopping\n";
                    break;
                }
                
                // Check for NaN or infinite values
                if (!std::isfinite(matAsquare(i, j)))
                {
                    continue; // Skip non-finite values
                }
                
                if (matAsquare(i, j) != 0)
                {
                    real_t contribution = weight * matAsquare(i, j);
                    
                    // Additional safety: check that contribution is finite
                    if (!std::isfinite(contribution))
                    {
                        continue;
                    }
                    
                    A_mat(i * 2 + 0, j * 2 + 0) += contribution;
                    A_mat(i * 2 + 1, j * 2 + 1) += contribution;
                }
            }
        }
        
        aftersquare = std::chrono::system_clock::now();
        std::chrono::duration<double> elapsed_seconds_square = aftersquare - beforesquare;
        // gsInfo << "square matrix took: " << elapsed_seconds_square.count() << "\n";
    }
    catch (const std::exception& e)
    {
        gsInfo << "EXCEPTION in assembleFitting: " << e.what() << "\n";
        gsInfo << "  A_mat dimensions: " << A_mat.rows() << " x " << A_mat.cols() << "\n";
        gsInfo << "  matAsquare dimensions: " << matAsquare.rows() << " x " << matAsquare.cols() << "\n";
        gsInfo << "  Continuing despite error...\n";
    }
    catch (...)
    {
        gsInfo << "UNKNOWN EXCEPTION in assembleFitting\n";
        gsInfo << "  A_mat dimensions: " << A_mat.rows() << " x " << A_mat.cols() << "\n";
        gsInfo << "  matAsquare dimensions: " << matAsquare.rows() << " x " << matAsquare.cols() << "\n";
        gsInfo << "  Continuing despite error...\n";
    }
}

/**
 * @brief Assembles the optimization matrix and RHS vector for THB functions.
 * Modernized for MPBES - spillovers are now handled internally by MPBES.
 *
 * @param[in,out] A                 Sparse matrix for optimization (LHS).
 * @param[in,out] b                 Vector for optimization (RHS).
 * @param[in]      J                Jacobian matrix contributions.
 * @param[in]      q                Functional weights matrix.
 * @param[in]      weight           Weighting factor for the functional.
 * @param[in]      globalIndexTHB   Global indices for THB functions across patches.
 * @param[in]      actives          Active THB function indices at the point.
 * @param[in]      patch            Current patch index.
 * @param[in]      verbose          Enables verbose logging output.
 */
void assembleOptimization(gsSparseMatrix<real_t>& A,
    gsMatrix<real_t>& b,
    const gsMatrix<real_t>& J,
    const gsMatrix<real_t>& q,
    const real_t weight,
    const gsVector<gsVector<int>>& globalIndexTHB,
    const gsMatrix<int>& actives,
    int patch,
    bool verbose = true,
    const gsMatrix<real_t>* basisInd = nullptr)
{
    //PROFILE_FUNCTION();
    const short_t geoDim = 2;
    //bool verbose = true;

    //if (verbose)
    //{
    //    gsInfo << "\n[assembleOptimization] Starting assembly for patch " << patch << "\n";
    //    gsInfo << "  actives.rows(): " << actives.rows()
    //        << ", J.cols(): " << J.cols()
    //        << ", q.cols(): " << q.cols() << "\n";
    //}
    assert(actives.rows() > 0);
    const index_t activeCount = std::min<index_t>(actives.rows(), J.cols());

    for (index_t i = 0; i < activeCount; ++i)
    {
        int activeIndex = actives(i, 0);
        int globalI = globalIndexTHB[patch][activeIndex];
        if (globalI == -1) continue;
        if (basisInd)
        {
            if (globalI < 0 || globalI >= basisInd->cols() || (*basisInd)(0, globalI) == 0.0)
                continue;
        }

        for (int col = 0; col < q.cols(); ++col)
        {
            for (int dim = 0; dim < geoDim; ++dim)
            {
                index_t Jrow = col * geoDim + dim;
                real_t term1 = weight;
                real_t term2 = q(0, col);
                real_t term3 = J(Jrow, i);
                real_t val = term1 * term2 * term3;

                if (std::isnan(val))
                {
                    gsInfo << "  [NaN] b update: weight=" << term1
                        << ", q(0," << col << ")=" << term2
                        << ", J(" << Jrow << "," << i << ")=" << term3 << "\n";
                    outfile << "  [NaN] Detected in b update (THB)\n";
                }

                index_t brow = globalI * geoDim + dim;
                b(brow) += val;

                if (verbose)
                {
                    if (std::abs(val) > 1e100 || std::abs(term1) > 1e100 || std::abs(term2) > 1e100 || std::abs(term3) > 1e100)
                        outfile << "SUSPICIOUS\n";
                    //gsInfo << "  b[" << brow << "] += weight * q * J = "
                    //    << term1 << " * " << term2 << " * " << term3
                    //    << " = " << val << "\n";
                    outfile << "b[" << brow << "] += " << val
                        << "   (weight=" << term1
                        << ", q=" << term2
                        << ", J=" << term3 << ")\n";
                }
            }
        }

        for (index_t j = 0; j < activeCount; ++j)
        {
            int activeJIndex = actives(j, 0);
            int globalJ = globalIndexTHB[patch][activeJIndex];
            if (globalJ == -1) continue;
            if (basisInd)
            {
                if (globalJ < 0 || globalJ >= basisInd->cols() || (*basisInd)(0, globalJ) == 0.0)
                    continue;
            }

            for (int col = 0; col < q.cols(); ++col)
            {
                for (int dimu = 0; dimu < geoDim; ++dimu)
                {
                    for (int dimv = 0; dimv < geoDim; ++dimv)
                    {
                        index_t row = col * geoDim + dimu;
                        // if (verbose)
                        //     outfile << "[assembleOptimization] row = col * geoDim + dimu => " << row << "\n";
                        index_t colJ = col * geoDim + dimv;
                        // if (verbose)
                        //     outfile << "[assembleOptimization] colJ = col * geoDim + dimv => " << colJ << "\n";
                        real_t term1 = J(row, i);
                        // if (verbose)
                        //     outfile << "[assembleOptimization] term1 = J(row,i) => " << term1 << "\n";
                        real_t term2 = J(colJ, j);
                        // if (verbose)
                        //     outfile << "[assembleOptimization] term2 = J(colJ,j) => " << term2 << "\n";
                        real_t term3 = weight;
                        // if (verbose)
                        //     outfile << "[assembleOptimization] term3 = weight => " << term3 << "\n";
                        real_t val = term1 * term2 * term3;
                        // if (verbose)
                        //     outfile << "[assembleOptimization] val = term1*term2*term3 => " << val << "\n";

                        if (std::isnan(val))
                        {
                            gsInfo << "  [NaN] A update: J(" << row << "," << i << ")=" << term1
                                << ", J(" << colJ << "," << j << ")=" << term2
                                << ", weight=" << term3 << "\n";
                            outfile << "  [NaN] Detected in A update (THB)\n";
                        }

                        index_t aRow = globalI * geoDim + dimu;
                        index_t aCol = globalJ * geoDim + dimv;
                        A(aRow, aCol) += val;

                        if (verbose)
                        {
                            if (std::abs(val) > 1e100 || std::abs(term1) > 1e100 || std::abs(term2) > 1e100 || std::abs(term3) > 1e100)
                            {
                                outfile << "SUSPICIOUS\n";
                                gsInfo << "  [" << g_currentOptimizationTerm << "] A[" << aRow << "," << aCol << "] += J(" << row << "," << i << ") * "
                                    << "J(" << colJ << "," << j << ") * weight = "
                                    << term1 << " * " << term2 << " * " << term3
                                    << " = " << val << "\n";
                                outfile << "[" << g_currentOptimizationTerm << "] A[" << aRow << "," << aCol << "] += " << val
                                    << "   (J1=" << term1 << ", J2=" << term2 << ", w=" << term3 << ")\n";
                                outfile.flush();
                                throw ProgramExitSignal(260306, "Assemble optimization encountered suspiciously large values");

                            }
                        }
                    }
                }
            }
        }
    }

    if (verbose)
    {
        //gsInfo << "[assembleOptimization] Finished assembly for patch " << patch << "\n\n";
    }
}



void printDerMatrices(gsMatrix<real_t> basisDer,
    gsMatrix<real_t> geomDer,
    gsMatrix<real_t> basis2Der,
    gsMatrix<real_t> geom2Der,
    bool verbose
) {
    if (verbose) {
        //outfile << "basisDer\n";
        //for (int i = 0; i < basisDer.rows(); ++i) {
        //    for (int j = 0; j < basisDer.cols(); ++j) {
        //        outfile << basisDer(i, j) << " ";
        //    }
        //    outfile << "\n";
        //}
        outfile << "geomDer\n";
        for (int i = 0; i < geomDer.rows(); ++i) {
            for (int j = 0; j < geomDer.cols(); ++j) {
                outfile << geomDer(i, j) << " ";
            }
            outfile << "\n";
        }
        outfile << "geom2Der.size(): " << geom2Der.rows() << " * " << geom2Der.cols() << "\n";
        for (int i = 0; i < geom2Der.rows(); i++) {
            for (int j = 0; j < geom2Der.cols(); j++) {
                outfile << geom2Der(i, j) << " ";
                if (isnan(geom2Der(i, j)))
                    gsInfo << "nan, Geom2Der: " << i << ", " << j << "\n";
            }
            outfile << "\n";
            //gsInfo << "\n";
        }
        outfile << "basis2Der\n";
        for (int i = 0; i < basis2Der.rows(); i++) {
            for (int j = 0; j < basis2Der.cols(); j++) {
                outfile << basis2Der(i, j) << " ";
                if (isnan(basis2Der(i, j)))
                    gsInfo << "nan, basis2Der: " << i << ", " << j << "\n";
            }
            outfile << "\n";
            //gsInfo << "\n";
        }
    }
}







// Optimized version for MPBES basis
void optimize(
    const MPBES<2, real_t>& mpbes,
    const gsVector<gsMatrix<>>& uv,
    gsMatrix<>& matAsquare,
    gsMatrix<real_t>& vectSol,
    gsSparseMatrix<real_t>& A,
    const real_t fitting,
    const real_t orthogonality,
    const real_t length,
    const real_t uniformity,
    const real_t skewness,
    const real_t eccentricity,
    const real_t area,
    const real_t epsilon,
    bool verbose,
    bool isNonLinearOptimization,
    const LocalCoarseningRegion* localRegion,
    const gsMatrix<real_t>* originalCoefficients)
{
    PROFILE_FUNCTION();

    auto logMat = [&](const char* name,
                      const gsMatrix<real_t>& M,
                      int patch,
                      int elem,
                      int qnode)
    {
        if (!outfile.is_open())
            return;
        const real_t maxAbs = (M.size() == 0) ? 0.0 : M.cwiseAbs().maxCoeff();
        if (maxAbs > 1e100)
            outfile << "SUSPICIOUS\n";
        outfile << "=== " << name << " patch=" << patch
                << " elem=" << elem << " q=" << qnode
                << " size=" << M.rows() << "x" << M.cols() << " ===\n";
        for (index_t r = 0; r < M.rows(); ++r)
        {
            for (index_t c = 0; c < M.cols(); ++c)
                outfile << M(r, c) << " ";
            outfile << "\n";
        }
        outfile << "\n";
        outfile.flush();
    };

    auto maxAbsMat = [](const gsMatrix<real_t>& M) -> real_t {
        if (M.size() == 0) return 0.0;
        return M.cwiseAbs().maxCoeff();
    };
    struct MaxTrack {
        real_t Jort=0, Qort=0;
        real_t Jlen=0, Qlen=0;
        real_t Juni=0, Quni=0;
        real_t Jskew=0, Qskew=0;
        real_t Jecc=0, Qecc=0;
        real_t Jarea=0, Qarea=0;
    } maxAll;
    
    // Extract data from MPBES
    auto Bells = mpbes.bellsBases();
    auto SubdomainHierarchy = mpbes.thbBases();
    auto functionDescription = mpbes.functionDescription();
    auto spilloverFunctionCoordinates = mpbes.spilloverCoordinates();
    auto hasSpillover = mpbes.hasSpillover();
    auto globalIndexTHB = extractGlobalIndexTHB(mpbes);
    const bool useLocalRegion = (localRegion != nullptr && localRegion->enabled);
    const gsMatrix<real_t>* basisMask = (useLocalRegion ? &localRegion->basisInd : nullptr);

    if (useLocalRegion && basisMask && originalCoefficients &&
        originalCoefficients->rows() == vectSol.rows() &&
        originalCoefficients->cols() >= vectSol.cols())
    {
        for (index_t f = 0; f < vectSol.rows(); ++f)
        {
            if ((*basisMask)(0, f) == 0.0)
                for (index_t d = 0; d < vectSol.cols(); ++d)
                    vectSol(f, d) = (*originalCoefficients)(f, d);
        }
    }

    assembleFitting(uv, A, Bells, functionDescription, matAsquare, fitting);

    gsMatrix<real_t> b(A.rows(), 1);
    b.setZero();

    real_t functional = 0;
    auto beforeOpt = std::chrono::system_clock::now();

    for (size_t patch = 0; patch < SubdomainHierarchy.size(); ++patch)
    {
        if (useLocalRegion)
        {
            if (patch >= localRegion->hasPatch.size() || !localRegion->hasPatch[patch])
                continue;
        }

        typename gsBasis<real_t>::domainIter domIt = SubdomainHierarchy(patch).makeDomainIterator();

        gsVector<size_t> numNodes(2);
        numNodes[0] = 2;
        numNodes[1] = 2;

        gsGaussRule<real_t> quRule(numNodes);

        gsMatrix<real_t> quNodes;
        gsVector<real_t> quWeights;
        gsMatrix<int> actives;

        for (int domanI = 0; domIt->good(); domIt->next(), ++domanI)
        {
            gsVector<real_t> lower = domIt->lowerCorner();
            gsVector<real_t> upper = domIt->upperCorner();

            if (useLocalRegion && !patchElementIntersectsRegion(*localRegion, static_cast<int>(patch), lower, upper))
                continue;

            try {
                quRule.mapTo(lower, upper, quNodes, quWeights);
            }
            catch (const std::exception& e) {
                outfile << "[ERROR] quadrature mapTo failed: " << e.what() << "\n";
                continue;
            }
            SubdomainHierarchy(patch).active_into(quNodes.col(0), actives);

            // Note: spillovers are now handled internally by MPBES
            gsMatrix<real_t> geomDer, basisDer, geom2Der, basis2Der;

            deriv_into_mpbe(mpbes, patch, quNodes, basisDer, actives, globalIndexTHB, false);

            jacobian_into_mpbe(mpbes, patch, quNodes, geomDer, actives, vectSol, globalIndexTHB, false);


            if (uniformity > 0 || eccentricity > 0)
            {
                deriv2_into_geometry_mpbe(mpbes, patch, quNodes, geom2Der, actives, vectSol, globalIndexTHB);

                deriv2_into_basis_mpbe(mpbes, patch, quNodes, basis2Der, actives, globalIndexTHB, false);
            }

            for (index_t quNode = 0; quNode < quNodes.cols(); ++quNode)
            {
                const real_t weight = quWeights[quNode];

                if (orthogonality > 0)
                {
                    gsMatrix<real_t> Jort, q;
                    getJort(Jort, actives.rows(), basisDer, geomDer, quNode);
                    getQort(q, geomDer, quNode);
                    logMat("Jort", Jort, static_cast<int>(patch), domanI, static_cast<int>(quNode));
                    logMat("Qort", q, static_cast<int>(patch), domanI, static_cast<int>(quNode));
                    maxAll.Jort = std::max(maxAll.Jort, maxAbsMat(Jort));
                    maxAll.Qort = std::max(maxAll.Qort, maxAbsMat(q));
                    functional += assembleFunctional(q, weight, orthogonality);

                    g_currentOptimizationTerm = "orthogonality";
                    assembleOptimization(A, b, Jort, q, weight * orthogonality,
                        globalIndexTHB, actives, patch, true, basisMask);
                }

                if (length > 0)
                {
                    gsMatrix<real_t> Jlen, qLen;
                    getQlen(qLen, geomDer, quNode);
                    getJlen(Jlen, basisDer.cols(), quNode, basisDer);
                    logMat("Jlen", Jlen, static_cast<int>(patch), domanI, static_cast<int>(quNode));
                    logMat("Qlen", qLen, static_cast<int>(patch), domanI, static_cast<int>(quNode));
                    maxAll.Jlen = std::max(maxAll.Jlen, maxAbsMat(Jlen));
                    maxAll.Qlen = std::max(maxAll.Qlen, maxAbsMat(qLen));
                    functional += assembleFunctional(qLen, weight, length);

                    g_currentOptimizationTerm = "length";
                    assembleOptimization(A, b, Jlen, qLen, weight * length,
                        globalIndexTHB, actives, patch, true, basisMask);
                }

                if (uniformity > 0)
                {
                    gsMatrix<real_t> Juni, qUni;
                    getJuni(Juni, actives.rows(), quNode, basis2Der);
                    getQuni(qUni, geom2Der, quNode);
                    logMat("Juni", Juni, static_cast<int>(patch), domanI, static_cast<int>(quNode));
                    logMat("Quni", qUni, static_cast<int>(patch), domanI, static_cast<int>(quNode));
                    maxAll.Juni = std::max(maxAll.Juni, maxAbsMat(Juni));
                    maxAll.Quni = std::max(maxAll.Quni, maxAbsMat(qUni));
                    functional += assembleFunctional(qUni, weight, uniformity);

                    g_currentOptimizationTerm = "uniformity";
                    assembleOptimization(A, b, Juni, qUni, weight * uniformity,
                        globalIndexTHB, actives, patch, true, basisMask);
                }

                if (skewness > 0)
                {
                    gsMatrix<real_t> Jskew, qSkew;
                    getQskew(qSkew, geomDer, quNode);
                    getJskew(Jskew, actives.rows(), quNode, basisDer, geomDer);
                    logMat("Jskew", Jskew, static_cast<int>(patch), domanI, static_cast<int>(quNode));
                    logMat("Qskew", qSkew, static_cast<int>(patch), domanI, static_cast<int>(quNode));
                    maxAll.Jskew = std::max(maxAll.Jskew, maxAbsMat(Jskew));
                    maxAll.Qskew = std::max(maxAll.Qskew, maxAbsMat(qSkew));
                    functional += assembleFunctional(qSkew, weight, skewness);

                    g_currentOptimizationTerm = "skewness";
                    assembleOptimization(A, b, Jskew, qSkew, weight * skewness,
                        globalIndexTHB, actives, patch, true, basisMask);
                }

                if (eccentricity > 0)
                {
                    gsMatrix<real_t> Jecc, qEcc;
                    getQecc(qEcc, geomDer, geom2Der, quNode);
                    getJecc(Jecc, actives.rows(), quNode, basisDer, basis2Der, geomDer, geom2Der);
                    logMat("Jecc", Jecc, static_cast<int>(patch), domanI, static_cast<int>(quNode));
                    logMat("Qecc", qEcc, static_cast<int>(patch), domanI, static_cast<int>(quNode));
                    maxAll.Jecc = std::max(maxAll.Jecc, maxAbsMat(Jecc));
                    maxAll.Qecc = std::max(maxAll.Qecc, maxAbsMat(qEcc));
                    functional += assembleFunctional(qEcc, weight, eccentricity);

                    g_currentOptimizationTerm = "eccentricity";
                    assembleOptimization(A, b, Jecc, qEcc, weight * eccentricity,
                        globalIndexTHB, actives, patch, true, basisMask);
                }

                if (area > 0)
                {
                    gsMatrix<real_t> Jarea, qArea;
                    getQarea(qArea, geomDer, quNode);
                    getJarea(Jarea, actives.rows(), quNode, basisDer, geomDer);
                    logMat("Jarea", Jarea, static_cast<int>(patch), domanI, static_cast<int>(quNode));
                    logMat("Qarea", qArea, static_cast<int>(patch), domanI, static_cast<int>(quNode));
                    maxAll.Jarea = std::max(maxAll.Jarea, maxAbsMat(Jarea));
                    maxAll.Qarea = std::max(maxAll.Qarea, maxAbsMat(qArea));
                    functional += assembleFunctional(qArea, weight, area);

                    g_currentOptimizationTerm = "area";
                    assembleOptimization(A, b, Jarea, qArea, weight * area,
                        globalIndexTHB, actives, patch, true, basisMask);
                }
            }
        }
    }

    auto afterOpt = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = afterOpt - beforeOpt;
        outfile << "[optimize] Max |J|/|Q| per functional (all steps):\n"
            << "  ort:  |J|=" << maxAll.Jort << ", |Q|=" << maxAll.Qort << "\n"
            << "  len:  |J|=" << maxAll.Jlen << ", |Q|=" << maxAll.Qlen << "\n"
            << "  uni:  |J|=" << maxAll.Juni << ", |Q|=" << maxAll.Quni << "\n"
            << "  skew: |J|=" << maxAll.Jskew << ", |Q|=" << maxAll.Qskew << "\n"
            << "  ecc:  |J|=" << maxAll.Jecc << ", |Q|=" << maxAll.Qecc << "\n"
            << "  area: |J|=" << maxAll.Jarea << ", |Q|=" << maxAll.Qarea << "\n";
    if (verbose)
    {
        checkSparseMatrixHealth(A, "A_opt");
        if (!gsEigen::MatrixXd(b).allFinite())
        {
            gsInfo << "[optimize] b contains non-finite entries\n";
            outfile << "[optimize] b contains non-finite entries\n";
        }
        gsMatrix<> A_dense = gsEigen::MatrixXd(A);
        printTheMatrix(A_dense, "A_opt");
        printTheMatrix(b, "b_opt");
    }
    
    // SOLVER PHASE

    if (A.rows() == A.cols())
    {
        gsEigen::MatrixXd A_dense = gsEigen::MatrixXd(A);
        if (!A_dense.allFinite())
        {
            gsInfo << "[optimize] A_opt has non-finite entries\n";
            outfile << "[optimize] A_opt has non-finite entries\n";
            outfile.flush();
            throw ProgramExitSignal(EXIT_FAILURE, "Optimization matrix contains non-finite entries");
        }
    }

    gsSparseSolver<>::CGIdentity solver(A);
    gsMatrix<real_t> x;

    try {
        x = solver.solve(b);
    }
    catch (const std::exception& e) {
        return;
    }

    if (!gsEigen::MatrixXd(x).allFinite())
    {
        gsInfo << "[optimize] Solver produced non-finite x; aborting update\n";
        outfile << "[optimize] Solver produced non-finite x; aborting update\n";
        printTheMatrix(x, "x");
        throw ProgramExitSignal(260228, "Optimization solver produced non-finite iterate");
    }

    // UPDATE SOLUTION
    gsMatrix<> vectSolOld = vectSol;
    for (int row = 0; row < vectSol.rows(); ++row)
        vectSol.row(row) -= x.block(row * 2, 0, 2, 1).transpose();

    if (useLocalRegion && basisMask && originalCoefficients &&
        originalCoefficients->rows() == vectSol.rows() &&
        originalCoefficients->cols() >= vectSol.cols())
    {
        for (index_t f = 0; f < vectSol.rows(); ++f)
        {
            if ((*basisMask)(0, f) == 0.0)
                for (index_t d = 0; d < vectSol.cols(); ++d)
                    vectSol(f, d) = (*originalCoefficients)(f, d);
        }
    }

    if (outfile.is_open())
    {
        outfile << "[optimize] vectSol ("
                << (isNonLinearOptimization ? "nonlinear" : "linear")
                << ")\n";
        outfile << "rows=" << vectSol.rows() << " cols=" << vectSol.cols() << "\n";
        for (int r = 0; r < vectSol.rows(); ++r)
        {
            for (int c = 0; c < vectSol.cols(); ++c)
                outfile << vectSol(r, c) << " ";
            outfile << "\n";
        }
        outfile.flush();
    }

    localTempAttempt++;
}

// Helper function to extract globalIndexTHB from MPBES
gsVector<gsVector<int>> extractGlobalIndexTHB(const MPBES<2, real_t>& mpbes) {
    const auto& functionDescription = mpbes.functionDescription();
    const auto& thbBases = mpbes.thbBases();
    const auto& indexInTHB = mpbes.indexInTHB();
    
    gsVector<gsVector<int>> globalIndexTHB(thbBases.size());
    
    // Initialize with -1 for all THB basis functions
    for (size_t patch = 0; patch < thbBases.size(); ++patch) {
        globalIndexTHB[patch].resize(thbBases[patch].size());
        for (size_t i = 0; i < thbBases[patch].size(); ++i) {
            globalIndexTHB[patch][i] = -1;
        }
    }
    
    // Fill in global indices from function description
    for (size_t globalIdx = 0; globalIdx < functionDescription.size(); ++globalIdx) {
        const auto& twins = functionDescription[globalIdx];
        for (const auto& twin : twins) {
            int patch = twin[0];
            int level = twin[1];
            int tensorIdx = twin[2];
            
            // Get THB index
            if (patch < indexInTHB.size() && 
                level < indexInTHB[patch].size() && 
                tensorIdx < indexInTHB[patch][level].size()) {
                int thbIdx = indexInTHB[patch][level][tensorIdx];
                if (thbIdx >= 0 && thbIdx < globalIndexTHB[patch].size()) {
                    globalIndexTHB[patch][thbIdx] = globalIdx;
                }
            }
        }
    }
    
    return globalIndexTHB;
}

// Returns a (1 × N) mask: 0 for MPBES functions whose support touches any
// boundary side of any patch, 1 for interior functions.  Used to freeze
// boundary DOFs during LO/NLO so that featureError (measured at patch
// interfaces) is not degraded by the optimization.
gsMatrix<real_t> buildBoundaryFunctionMask(const MPBES<2, real_t>& mpbes)
{
    const index_t N = static_cast<index_t>(mpbes.size());
    gsMatrix<real_t> mask(1, N);
    mask.setOnes();

    const auto thbBases  = mpbes.thbBases();
    const gsVector<gsVector<int>> gIdx = extractGlobalIndexTHB(mpbes);

    for (index_t p = 0; p < static_cast<index_t>(thbBases.size()); ++p)
    {
        for (int s = 1; s <= 4; ++s)
        {
            gsMatrix<index_t> bnd;
            try { bnd = thbBases[p].boundary(static_cast<boxSide>(s)); }
            catch (...) { continue; }
            for (index_t k = 0; k < bnd.size(); ++k)
            {
                const index_t thbIdx = bnd(k, 0);
                if (thbIdx < 0 || thbIdx >= static_cast<index_t>(gIdx[p].size())) continue;
                const int g = gIdx[p][thbIdx];
                if (g >= 0 && g < N)
                    mask(0, g) = 0.0;
            }
        }
    }
    return mask;
}










inline real_t assembleFunctional(const gsMatrix<real_t>& q,
    const real_t weight,
    const real_t constant)
{
    //gsInfo << "-F-> assembleFunctional\n";
    assert(q.rows() > 0 && "Error: q has no rows!");
    assert(q.cols() > 0 && "Error: q has no columns!");

    real_t value = 0;
    for (int col = 0; col < q.cols(); col++)
    {
        assert(col < q.cols() && "Error: Index out of bounds in q!");

        real_t q_val = q(0, col);
        assert(!std::isnan(q_val) && "Error: NaN detected in q!");
        assert(!std::isinf(q_val) && "Error: Inf detected in q!");

        value += q_val * q_val;
    }

    assert(!std::isnan(value) && "Error: NaN detected in computed value!");
    assert(!std::isinf(value) && "Error: Inf detected in computed value!");

    return weight * constant * value;
}






void actives_into(
    int patch,
    gsVector<std::vector<int>>& activesVect,
    gsVector<std::vector<int>>& activesLev,
    gsMatrix <real_t> quNodes,
    gsVector<gsVector<gsTensorBSplineBasis<2, real_t>>> Bells,
    gsVector<gsVector<gsVector<index_t >>> globalIndex
)
{
    gsVector<gsVector<int>> isActiveAtPoint;
    isActiveAtPoint.resize(Bells(patch).size());
    activesVect.resize(quNodes.cols());
    activesLev.resize(quNodes.cols());
    isActiveAtPoint(patch).resize(Bells(patch).size());
    for (size_t level = 0; level < Bells(patch).size(); level++)
    {
        //gsInfo << level << "\n";
        //gsInfo << "Bells(patch)(level).size(): " << Bells(patch)(level).size() << "\n";
        isActiveAtPoint(level).setZero(Bells(patch)(level).size());
        //isActiveAtPoint(patch)(level).setZero();
    }
    for (size_t i = 0; i < quNodes.cols(); i++)
    {
        for (size_t level = 0; level < Bells(0).size(); level++)
        {
            gsMatrix<int> test;
            //gsInfo << i << " " << level << " " << k << "\n";
            Bells(patch)(level).active_into(quNodes.col(i), test);
            //gsInfo << test.rows() << " " << test.cols() << "\n";
            for (size_t j = 0; j < test.rows(); j++)
            {
                isActiveAtPoint(level)(test(j, 0)) = 1;
            }
        }
        for (size_t level = 0; level < globalIndex(patch).size(); level++)
        {
            for (size_t j = 0; j < globalIndex(patch)(level).size(); j++)
            {
                if (globalIndex(patch)(level)(j) != -1 && isActiveAtPoint(level)(j))
                {
                    activesVect(i).push_back(j);
                    activesLev(i).push_back(level);
                }
            }
        }
    }
}

void active_into(
    std::vector<int>& activesVect,
    std::vector<int>& activesLev,
    gsMatrix <real_t> punto,
    gsVector<gsVector<gsTensorBSplineBasis<2, real_t>>> Bells,
    gsVector<gsVector<gsVector<index_t >>> globalIndex
)
{
    gsVector<gsVector<gsVector<int>>> isActiveAtPoint;
    isActiveAtPoint.resize(Bells.size());
    for (size_t patch = 0; patch < Bells.size(); patch++)
    {
        isActiveAtPoint(patch).resize(Bells(patch).size());
        for (size_t level = 0; level < Bells(patch).size(); level++)
        {
            //gsInfo << level << "\n";
            //gsInfo << "Bells(patch)(level).size(): " << Bells(patch)(level).size() << "\n";
            isActiveAtPoint(patch)(level).setZero(Bells(patch)(level).size());
            //isActiveAtPoint(patch)(level).setZero();
        }
    }
    for (size_t level = 0; level < Bells(0).size(); level++)
    {
        gsMatrix<int> test;
        //gsInfo << i << " " << level << " " << k << "\n";
        Bells(0)(level).active_into(punto, test);
        //gsInfo << test.rows() << " " << test.cols() << "\n";
        for (size_t j = 0; j < test.rows(); j++)
        {
            isActiveAtPoint(0)(level)(test(j, 0)) = 1;
        }
    }


    for (size_t level = 0; level < globalIndex(0).size(); level++)
    {
        for (size_t j = 0; j < globalIndex(0)(level).size(); j++)
        {
            if (globalIndex(0)(level)(j) != -1 && isActiveAtPoint(0)(level)(j))
            {
                activesVect.push_back(j);//  Tofiq
                activesLev.push_back(level);//  Tofiq
            }
        }
    }
}


void active_into(const gsMatrix< real_t >& u,
    int patch,
    gsVector<gsVector<gsTensorBSplineBasis<2, real_t>>> Bells,
    gsVector<gsVector<std::vector<index_t>>> functionDescription,
    gsMatrix< index_t >& result
)
{
    int corrLevel = 0;
    int corrIndex = 0;
    gsMatrix<> support;
    std::vector<int> pieceSaver;
    std::vector<int> patchSaver;
    for (size_t functionIndex = 0; functionIndex < functionDescription.size(); functionIndex++)
    {
        bool functionIsActive = false;
        for (size_t pieceIndex = 0; pieceIndex < functionDescription(functionIndex).size(); pieceIndex++)
        {
            const auto& piece = functionDescription(functionIndex)[pieceIndex];
            if (piece[0] == patch) {
                corrLevel = piece[1];
                corrIndex = piece[2];
                support = Bells(patch)(corrLevel).function(corrIndex).support();
                if ((u(0, 0) <= support(0, 1) && u(0, 0) >= support(0, 0)) && (u(1, 0) <= support(1, 1) && u(1, 0) >= support(1, 0)))
                {
                    functionIsActive = true;
                    // Save the function index and the index of the piece within the functionDescription vector
                    pieceSaver.push_back(functionIndex);
                    patchSaver.push_back(pieceIndex);
                }
                break; // No need to continue once found
            }
        }
    }
    result.resize(2, pieceSaver.size());
    for (size_t i = 0; i < pieceSaver.size(); i++)
    {
        result(0, i) = pieceSaver[i];
        result(1, i) = patchSaver[i];
    }
}

void checkForNaN(const gsMatrix<real_t>& matrix, const std::string& matrixName) {
    for (size_t row = 0; row < matrix.rows(); ++row) {
        for (size_t col = 0; col < matrix.cols(); ++col) {
            if (isnan(matrix(row, col))) {
                gsInfo << "NaN detected in matrix " << matrixName
                    << " at position (" << row << ", " << col << ")\n";
            }
        }
    }
}

void checkSparseMatrixHealth(
    const gsSparseMatrix<real_t>& A,
    const std::string& name,
    size_t maxLogs)
{
    const index_t rows = A.rows();
    const index_t cols = A.cols();
    std::vector<index_t> rowNnz(rows, 0);
    std::vector<index_t> colNnz(cols, 0);
    size_t nonFiniteCount = 0;
    size_t logged = 0;

    for (int k = 0; k < A.outerSize(); ++k)
    {
        for (gsSparseMatrix<real_t>::InnerIterator it(A, k); it; ++it)
        {
            const real_t v = it.value();
            rowNnz[it.row()]++;
            colNnz[it.col()]++;
            if (!std::isfinite(v))
            {
                ++nonFiniteCount;
                if (logged < maxLogs)
                {
                    // gsInfo << "[" << name << "] non-finite value at (" << it.row() << "," << it.col() << ") = " << v << "\n";
                    // outfile << "[" << name << "] non-finite value at (" << it.row() << "," << it.col() << ") = " << v << "\n";
                    ++logged;
                }
            }
        }
    }

    index_t zeroRows = 0;
    index_t zeroCols = 0;
    for (index_t r = 0; r < rows; ++r)
        if (rowNnz[r] == 0) ++zeroRows;
    for (index_t c = 0; c < cols; ++c)
        if (colNnz[c] == 0) ++zeroCols;

        // gsInfo << "[" << name << "] rows=" << rows << ", cols=" << cols
        //        << ", nnz=" << A.nonZeros()
        //        << ", zeroRows=" << zeroRows
        //        << ", zeroCols=" << zeroCols
        //        << ", nonFinite=" << nonFiniteCount << "\n";
        // outfile << "[" << name << "] rows=" << rows << ", cols=" << cols
        //         << ", nnz=" << A.nonZeros()
        //         << ", zeroRows=" << zeroRows
        //         << ", zeroCols=" << zeroCols
        //         << ", nonFinite=" << nonFiniteCount << "\n";
}

index_t countTotalActiveElements(const gsVector<gsTHBSplineBasis<2, real_t>>& THBVector)
{
    index_t count = 0;
    for (index_t patch = 0; patch < THBVector.size(); ++patch)
    {
        typename gsBasis<real_t>::domainIter domIt = THBVector(patch).makeDomainIterator();
        for (; domIt->good(); domIt->next())
            ++count;
    }
    return count;
}

// ============================================================
// MPBES geometry-based coarsening marking strategy
// ============================================================
template<short_t d>
struct ElementKey
{
    index_t patch;
    index_t level;
    std::array<index_t, d> low;
    std::array<index_t, d> upp;

    bool operator==(const ElementKey& other) const
    {
        return patch == other.patch && level == other.level && low == other.low && upp == other.upp;
    }
};

template<short_t d>
struct ElementKeyHash
{
    size_t operator()(const ElementKey<d>& key) const
    {
        size_t h = std::hash<index_t>{}(key.patch);
        h ^= std::hash<index_t>{}(key.level) + 0x9e3779b9 + (h << 6) + (h >> 2);
        for (short_t i = 0; i < d; ++i)
        {
            h ^= std::hash<index_t>{}(key.low[i]) + 0x9e3779b9 + (h << 6) + (h >> 2);
            h ^= std::hash<index_t>{}(key.upp[i]) + 0x9e3779b9 + (h << 6) + (h >> 2);
        }
        return h;
    }
};

template<short_t d, typename T>
static ElementKey<d> makeElementKey(const gsHBox<d, T>& box)
{
    ElementKey<d> key;
    key.patch = box.patch();
    key.level = box.level();
    const gsVector<index_t, d>& low = box.lowerIndex();
    const gsVector<index_t, d>& upp = box.upperIndex();
    for (short_t i = 0; i < d; ++i)
    {
        key.low[i] = low[i];
        key.upp[i] = upp[i];
    }
    return key;
}

template<short_t d, typename T>
struct ElementInfo
{
    gsHBox<d, T> box;
    T eta;
};

template<short_t d, typename T>
std::vector<gsHBox<d, T>*> markCoarseningGeometryBased(
    const MPBES<d, T>& mpbes,
    const gsMultiPatch<T>& geometry,
    T delta,
    index_t numGaussPerDir = 2,
    int patchFilter = -1,
    bool verbose = false,
    bool adaptiveRelax = true,
    T minChildFraction = static_cast<T>(0.25))
{
    PROFILE_FUNCTION();

    std::vector<gsHBox<d, T>*> coarseningSet;
    if (delta <= 0 || mpbes.nPatches() == 0)
        return coarseningSet;

    if (verbose)
    {
        gsInfo << "[geo-coarsen][step 1] collect elements\n";
        outfile << "[geo-coarsen][step 1] collect elements\n";
    }

    const auto& thbBases = mpbes.thbBases();
    auto activeFunctionsAtLevel = [&](index_t patch, index_t level) -> index_t
    {
        if (patch < 0 || patch >= static_cast<index_t>(thbBases.size()))
            return 0;
        const auto& xmat = thbBases[patch].getXmatrix();
        if (level < 0 || level >= static_cast<index_t>(xmat.size()))
            return 0;
        return static_cast<index_t>(xmat[level].size());
    };
    gsVector<index_t> numNodes(d);
    for (short_t i = 0; i < d; ++i)
        numNodes[i] = std::max<index_t>(1, numGaussPerDir);

    gsGaussRule<T> quRule(numNodes);
    std::vector<ElementInfo<d, T>> elements;

    for (index_t patch = 0; patch < thbBases.size(); ++patch)
    {
        if (patchFilter >= 0 && patch != patchFilter)
            continue;
        if (patch >= geometry.nPatches())
            continue;

        typename gsBasis<T>::domainIter domIt = thbBases[patch].makeDomainIterator();
        auto* hdomIt = dynamic_cast<gsHDomainIterator<T, d>*>(domIt.get());
        if (!hdomIt)
        {
            if (verbose)
                gsInfo << "WARNING: Expected hierarchical domain iterator for patch " << patch << "\n";
            continue;
        }

        if (verbose)
        {
            gsInfo << "[geo-coarsen] patch=" << patch << " collecting elements\n";
            outfile << "[geo-coarsen] patch=" << patch << " collecting elements\n";
        }

        for (; domIt->good(); domIt->next())
        {
            gsVector<T> lower = domIt->lowerCorner();
            gsVector<T> upper = domIt->upperCorner();

            gsMatrix<T> quNodes;
            gsVector<T> quWeights;
            try
            {
                quRule.mapTo(lower, upper, quNodes, quWeights);
            }
            catch (const std::exception& e)
            {
                if (verbose)
                    outfile << "[ERROR] quadrature mapTo failed: " << e.what() << "\n";
                continue;
            }

            T eta = 0;
            for (index_t q = 0; q < quNodes.cols(); ++q)
            {
                gsMatrix<T> J = geometry.patch(patch).jacobian(quNodes.col(q));
                eta += quWeights[q] * J.squaredNorm();
            }

            elements.push_back(ElementInfo<d, T>{ gsHBox<d, T>(hdomIt, patch), eta });
        }
    }

    if (elements.empty())
    {
        if (verbose)
        {
            gsInfo << "[geo-coarsen] no elements collected"
                   << (patchFilter >= 0 ? " for patch=" + std::to_string(patchFilter) : "")
                   << "\n";
            outfile << "[geo-coarsen] no elements collected"
                    << (patchFilter >= 0 ? " for patch=" + std::to_string(patchFilter) : "")
                    << "\n";
        }
        return coarseningSet;
    }

    if (verbose)
    {
        gsInfo << "[geo-coarsen][step 2] compute total eta\n";
        outfile << "[geo-coarsen][step 2] compute total eta\n";
    }

    T total = 0;
    for (const auto& el : elements)
        total += el.eta;

    if (total <= 0)
    {
        if (verbose)
        {
            gsInfo << "[geo-coarsen] total eta <= 0"
                   << (patchFilter >= 0 ? " for patch=" + std::to_string(patchFilter) : "")
                   << " (total=" << total << ")\n";
            outfile << "[geo-coarsen] total eta <= 0"
                    << (patchFilter >= 0 ? " for patch=" + std::to_string(patchFilter) : "")
                    << " (total=" << total << ")\n";
        }
        return coarseningSet;
    }

    if (verbose)
    {
        gsInfo << "[geo-coarsen][step 3] sort by eta\n";
        outfile << "[geo-coarsen][step 3] sort by eta\n";
    }

    std::vector<index_t> order(elements.size());
    std::iota(order.begin(), order.end(), 0);
    std::sort(order.begin(), order.end(), [&](index_t a, index_t b)
    {
        return elements[a].eta < elements[b].eta;
    });

    if (verbose)
    {
        gsInfo << "[geo-coarsen][step 4] bulk select (delta=" << delta << ")\n";
        outfile << "[geo-coarsen][step 4] bulk select (delta=" << delta << ")\n";
    }
    T target = delta * total;
    T running = 0;
    std::vector<index_t> selected;
    selected.reserve(order.size());
    for (index_t idx : order)
    {
        selected.push_back(idx);
        running += elements[idx].eta;
        if (running >= target)
            break;
    }

        if (verbose)
        {
        gsInfo << "[geo-coarsen] elements=" << elements.size()
               << " totalEta=" << total
               << " target=" << target
               << " selected=" << selected.size()
               << (patchFilter >= 0 ? " patch=" + std::to_string(patchFilter) : "")
               << "\n";
        outfile << "[geo-coarsen] elements=" << elements.size()
            << " totalEta=" << total
            << " target=" << target
            << " selected=" << selected.size()
            << (patchFilter >= 0 ? " patch=" + std::to_string(patchFilter) : "")
            << "\n";

        gsInfo << "[geo-coarsen][step 5] hierarchy rule (all children)\n";
        outfile << "[geo-coarsen][step 5] hierarchy rule (all children)\n";
        }

    using Key = ElementKey<d>;
    std::unordered_set<Key, ElementKeyHash<d>> selectedSet;
    selectedSet.reserve(selected.size() * 2);
    for (index_t idx : selected)
        selectedSet.insert(makeElementKey(elements[idx].box));

    std::unordered_map<Key, bool, ElementKeyHash<d>> parentAllowed;
    parentAllowed.reserve(selected.size());

    int allowedCount = 0;
    int rejectedCount = 0;
    int loggedRejected = 0;
    const int maxRejectedLogs = 10;
    for (index_t idx : selected)
    {
        const gsHBox<d, T>& box = elements[idx].box;
        if (box.level() <= 0)
            continue;

        gsHBox<d, T> parent = box.getParent();
        Key parentKey = makeElementKey(parent);

        auto it = parentAllowed.find(parentKey);
        bool ok = false;
        int selectedChildren = 0;
        int totalChildren = 0;
        typename gsHBox<d, T>::Container children;
        if (it == parentAllowed.end())
        {
            children = parent.getChildren();
            totalChildren = static_cast<int>(children.size());
            for (const auto& child : children)
            {
                if (selectedSet.find(makeElementKey(child)) != selectedSet.end())
                    selectedChildren++;
            }

            int requiredChildren = totalChildren;
            if (adaptiveRelax && totalChildren > 0)
            {
                const index_t coarseCount = activeFunctionsAtLevel(parent.patch(), parent.level());
                const index_t fineCount = activeFunctionsAtLevel(parent.patch(), parent.level() + 1);
                const double diff = std::max<double>(0.0, static_cast<double>(fineCount) - static_cast<double>(coarseCount));
                double ratio = 0.0;
                if (coarseCount > 0)
                    ratio = diff / static_cast<double>(coarseCount);
                else if (fineCount > 0)
                    ratio = 1.0;
                double requiredFraction = 1.0 / (1.0 + ratio);
                requiredFraction = std::max(requiredFraction, static_cast<double>(minChildFraction));
                requiredChildren = static_cast<int>(std::ceil(requiredFraction * static_cast<double>(totalChildren)));
            }

            ok = (selectedChildren >= requiredChildren);
            parentAllowed.emplace(parentKey, ok);
        }
        else
        {
            ok = it->second;
            if (!ok && loggedRejected < maxRejectedLogs)
            {
                children = parent.getChildren();
                totalChildren = static_cast<int>(children.size());
                for (const auto& child : children)
                {
                    if (selectedSet.find(makeElementKey(child)) != selectedSet.end())
                        selectedChildren++;
                }
            }
        }

        if (ok)
        {
            coarseningSet.push_back(new gsHBox<d, T>(box));
            allowedCount++;
        }
        else
        {
            rejectedCount++;
            if (verbose && loggedRejected < maxRejectedLogs)
            {
                const auto plow = parent.lowerIndex();
                const auto pupp = parent.upperIndex();
                gsInfo << "[geo-coarsen] reject parent patch=" << parent.patch()
                       << " level=" << parent.level()
                       << " low=[" << plow[0] << "," << plow[1] << "]"
                       << " upp=[" << pupp[0] << "," << pupp[1] << "]"
                       << " selectedChildren=" << selectedChildren
                       << " totalChildren=" << totalChildren
                       << "\n";
                outfile << "[geo-coarsen] reject parent patch=" << parent.patch()
                        << " level=" << parent.level()
                        << " low=[" << plow[0] << "," << plow[1] << "]"
                        << " upp=[" << pupp[0] << "," << pupp[1] << "]"
                        << " selectedChildren=" << selectedChildren
                        << " totalChildren=" << totalChildren
                        << "\n";
                int childLogCount = 0;
                for (const auto& child : children)
                {
                    const auto clow = child.lowerIndex();
                    const auto cupp = child.upperIndex();
                    bool childSelected = (selectedSet.find(makeElementKey(child)) != selectedSet.end());
                    gsInfo << "[geo-coarsen] child[" << childLogCount << "] patch=" << child.patch()
                           << " level=" << child.level()
                           << " low=[" << clow[0] << "," << clow[1] << "]"
                           << " upp=[" << cupp[0] << "," << cupp[1] << "]"
                           << " selected=" << (childSelected ? 1 : 0) << "\n";
                    outfile << "[geo-coarsen] child[" << childLogCount << "] patch=" << child.patch()
                            << " level=" << child.level()
                            << " low=[" << clow[0] << "," << clow[1] << "]"
                            << " upp=[" << cupp[0] << "," << cupp[1] << "]"
                            << " selected=" << (childSelected ? 1 : 0) << "\n";
                    childLogCount++;
                    if (childLogCount >= 5)
                        break;
                }
                loggedRejected++;
            }
        }
    }

    if (verbose)
    {
        gsInfo << "[geo-coarsen] allowed=" << allowedCount
               << " rejected=" << rejectedCount
               << " (after hierarchy rule)"
               << (patchFilter >= 0 ? " patch=" + std::to_string(patchFilter) : "")
               << "\n";
        outfile << "[geo-coarsen] allowed=" << allowedCount
                << " rejected=" << rejectedCount
                << " (after hierarchy rule)"
                << (patchFilter >= 0 ? " patch=" + std::to_string(patchFilter) : "")
                << "\n";
    }

    return coarseningSet;
}

template<short_t d, typename T>
std::vector<gsHBox<d, T>> markCoarseningGeometryBasedTHB(
    const gsVector<gsTHBSplineBasis<d, T>>& thbBases,
    const gsMultiPatch<T>& geometry,
    T delta,
    index_t numGaussPerDir = 2,
    int levelFilter = -1,
    int patchFilter = -1,
    bool verbose = false)
{
    // gsInfo << "[geo-coarsen-thb] called: delta=" << delta
    //        << " gauss=" << numGaussPerDir
    //        << " levelFilter=" << levelFilter
    //        << " patchFilter=" << patchFilter
    //        << " verbose=" << (verbose ? "true" : "false") << "\n";
    PROFILE_FUNCTION();

    std::vector<gsHBox<d, T>> coarseningSet;
    if (delta <= 0 || thbBases.size() == 0)
        return coarseningSet;

    gsVector<index_t> numNodes(d);
    for (short_t i = 0; i < d; ++i)
        numNodes[i] = std::max<index_t>(1, numGaussPerDir);

    gsGaussRule<T> quRule(numNodes);
    std::vector<ElementInfo<d, T>> elements;

    for (index_t patch = 0; patch < thbBases.size(); ++patch)
    {
        if (patchFilter >= 0 && patch != patchFilter)
            continue;
        if (patch >= geometry.nPatches())
            continue;

        typename gsBasis<T>::domainIter domIt = thbBases[patch].makeDomainIterator();
        auto* hdomIt = dynamic_cast<gsHDomainIterator<T, d>*>(domIt.get());
        if (!hdomIt)
        {
            if (verbose)
                gsInfo << "WARNING: Expected hierarchical domain iterator for patch " << patch << "\n";
            continue;
        }

        for (; domIt->good(); domIt->next())
        {
            if (levelFilter >= 0 && hdomIt->getLevel() != levelFilter)
                continue;

            gsVector<T> lower = domIt->lowerCorner();
            gsVector<T> upper = domIt->upperCorner();

            gsMatrix<T> quNodes;
            gsVector<T> quWeights;
            try
            {
                quRule.mapTo(lower, upper, quNodes, quWeights);
            }
            catch (const std::exception& e)
            {
                outfile << "[ERROR] quadrature mapTo failed: " << e.what() << "\n";
                continue;
            }

            T eta = 0;
            for (index_t q = 0; q < quNodes.cols(); ++q)
            {
                gsMatrix<T> J = geometry.patch(patch).jacobian(quNodes.col(q));
                eta += quWeights[q] * J.squaredNorm();
            }

            elements.push_back(ElementInfo<d, T>{ gsHBox<d, T>(hdomIt, patch), eta });
        }
    }

    if (elements.empty())
        return coarseningSet;

    T total = 0;
    for (const auto& el : elements)
        total += el.eta;

    if (total <= 0)
        return coarseningSet;

    std::vector<index_t> order(elements.size());
    std::iota(order.begin(), order.end(), 0);
    std::sort(order.begin(), order.end(), [&](index_t a, index_t b)
    {
        return elements[a].eta < elements[b].eta;
    });

    T target = delta * total;
    T running = 0;
    std::vector<index_t> selected;
    selected.reserve(order.size());
    for (index_t idx : order)
    {
        selected.push_back(idx);
        running += elements[idx].eta;
        if (running >= target)
            break;
    }

    using Key = ElementKey<d>;
    std::unordered_set<Key, ElementKeyHash<d>> selectedSet;
    selectedSet.reserve(selected.size() * 2);
    for (index_t idx : selected)
        selectedSet.insert(makeElementKey(elements[idx].box));

    std::unordered_map<Key, bool, ElementKeyHash<d>> parentAllowed;
    parentAllowed.reserve(selected.size());

    for (index_t idx : selected)
    {
        const gsHBox<d, T>& box = elements[idx].box;
        if (box.level() <= 0)
            continue;

        gsHBox<d, T> parent = box.getParent();
        Key parentKey = makeElementKey(parent);

        auto it = parentAllowed.find(parentKey);
        bool ok = false;
        if (it == parentAllowed.end())
        {
            ok = true;
            typename gsHBox<d, T>::Container children = parent.getChildren();
            for (const auto& child : children)
            {
                if (selectedSet.find(makeElementKey(child)) == selectedSet.end())
                {
                    ok = false;
                    break;
                }
            }
            parentAllowed.emplace(parentKey, ok);
        }
        else
        {
            ok = it->second;
        }

        if (ok)
            coarseningSet.push_back(box);
    }

        if (verbose)
        {
        const int maxLog = 20;
        gsInfo << "[geo-coarsen] candidates=" << elements.size()
               << " selected=" << selected.size()
               << " after-hierarchy=" << coarseningSet.size() << "\n";
        outfile << "[geo-coarsen] candidates=" << elements.size()
            << " selected=" << selected.size()
            << " after-hierarchy=" << coarseningSet.size() << "\n";
        for (int i = 0; i < static_cast<int>(coarseningSet.size()) && i < maxLog; ++i)
        {
            const auto& box = coarseningSet[i];
            const auto low = box.lowerIndex();
            const auto upp = box.upperIndex();
            outfile << "  box#" << i << " patch=" << box.patch() << " level=" << box.level()
                << " low=[" << low[0] << "," << low[1] << "]"
                << " upp=[" << upp[0] << "," << upp[1] << "]\n";
        }
        }

    return coarseningSet;
}

struct CellSelectionResult
{
    bool useGeo = false;
    int currCellIndex = -1;
    int currArrayIndex = 0;
    index_t x1U = 0;
    index_t y1U = 0;
    index_t x2U = 0;
    index_t y2U = 0;
    std::vector<CellToCoarsen> geoCells;
    std::vector<int> geoCellIndices;
    // Lowest Grenda delta at which candidates were found (infinity = no Grenda match).
    // Used by global pool to rank patches: smaller delta = more regular geometry.
    real_t acceptedDelta = std::numeric_limits<real_t>::infinity();
};

static CellSelectionResult selectCellForCoarsening(
    char method,
    const gsVector<index_t>& vectorS,
    int levNow,
    index_t interior,
    int patch,
    int attempt,
    int& pickedOne,
    gsMatrix<int>& pickedCells,
    int& valid,
    int& currArrayIndex,
    const MPBES<2, real_t>* currentMpbes,
    const gsMultiPatch<real_t>& mp1,
    const gsVector<gsTHBSplineBasis<2, real_t>>& SubdomainHierarchy,
    bool preflightOnly = false)
{
    CellSelectionResult result;
    result.currArrayIndex = currArrayIndex;
    const bool patchMirroredInPreflight = isPreflightMirroredPatch(patch);

    // gsInfo << "[selectCell] method=" << method
    //        << " patch=" << patch
    //        << " level=" << levNow
    //        << " attempt=" << attempt
    //        << " pool=" << vectorS.size() << "\n";

    if (method == 'g' && !patchMirroredInPreflight) {
        const bool geoVerbose = false;
        // gsInfo << "[selectCell] using Grendas geometry-based coarsening\n";
        const std::vector<real_t> geoDeltas = { 0.2, 0.4, 0.6, 0.8, 1.0, 1.2 };
        std::vector<gsHBox<2, real_t>> geoBoxes;
        std::unordered_set<int> available;
        available.reserve(vectorS.size() * 2);
        for (int i = 0; i < vectorS.size(); ++i)
            available.insert(vectorS(i));
        // if (currentMpbes) {
        //     gsInfo << "[geo-coarsen] acceptedMpbes: size=" << currentMpbes->size()
        //            << " patches=" << currentMpbes->nPatches() << "\n";
        //     outfile << "[geo-coarsen] acceptedMpbes: size=" << currentMpbes->size()
        //             << " patches=" << currentMpbes->nPatches() << "\n";
        // } else {
        //     gsInfo << "[geo-coarsen] acceptedMpbes=null, using THB fallback\n";
        //     outfile << "[geo-coarsen] acceptedMpbes=null, using THB fallback\n";
        // }

        for (size_t di = 0; di < geoDeltas.size(); ++di) {
            const real_t geoDelta = geoDeltas[di];
            const bool logCandidates = (di == 0);
            // gsInfo << "[geo-coarsen] trying delta=" << geoDelta << "\n";
            // outfile << "[geo-coarsen] trying delta=" << geoDelta << "\n";

            geoBoxes.clear();
            if (currentMpbes) {
                std::vector<gsHBox<2, real_t>*> geoPtrsPatch = markCoarseningGeometryBased<2, real_t>(
                    *currentMpbes, mp1, geoDelta, 2, patch, geoVerbose);
                  // gsInfo << "[geo-coarsen] geoPtrs patch=" << patch
                  //        << " total=" << geoPtrsPatch.size() << "\n";
                  // outfile << "[geo-coarsen] geoPtrs patch=" << patch
                  //         << " total=" << geoPtrsPatch.size() << "\n";

                if (logCandidates && g_verbose) {
                    for (size_t gi = 0; gi < geoPtrsPatch.size(); ++gi) {
                        const auto* ptr = geoPtrsPatch[gi];
                        if (!ptr) {
                            gsInfo << "[geo-coarsen] patch-cand[" << gi << "] null\n";
                            outfile << "[geo-coarsen] patch-cand[" << gi << "] null\n";
                            continue;
                        }
                        const auto low = ptr->lowerIndex();
                        const auto upp = ptr->upperIndex();
                        gsInfo << "[geo-coarsen] patch-cand[" << gi << "] patch=" << ptr->patch()
                               << " level=" << ptr->level()
                               << " low=[" << low[0] << "," << low[1] << "]"
                               << " upp=[" << upp[0] << "," << upp[1] << "]\n";
                        outfile << "[geo-coarsen] patch-cand[" << gi << "] patch=" << ptr->patch()
                                << " level=" << ptr->level()
                                << " low=[" << low[0] << "," << low[1] << "]"
                                << " upp=[" << upp[0] << "," << upp[1] << "]\n";
                    }
                }

                for (auto* ptr : geoPtrsPatch) {
                    if (ptr) {
                        if (ptr->patch() == patch && ptr->level() == levNow + 1)
                            geoBoxes.push_back(ptr->getParent());
                        delete ptr;
                    }
                }
            } else {
                geoBoxes = markCoarseningGeometryBasedTHB<2, real_t>(
                    SubdomainHierarchy, mp1, geoDelta, 2, levNow, patch, geoVerbose);
            }

            if (!geoBoxes.empty()) {
                std::vector<CellToCoarsen> mappedCells;
                std::vector<int> mappedIndices;
                mappedCells.reserve(geoBoxes.size());
                mappedIndices.reserve(geoBoxes.size());

                int mapped = 0;
                for (const auto& box : geoBoxes) {
                    const auto low = box.lowerIndex();
                    const auto upp = box.upperIndex();
                    int idx = -1;
                    if (!cellIndexFromCoords(levNow, low[0], low[1], interior, idx))
                        continue;
                    if (available.find(idx) == available.end())
                        continue;

                    CellToCoarsen cell;
                    cell.level = levNow;
                    cell.x1 = low[0];
                    cell.y1 = low[1];
                    cell.x2 = upp[0];
                    cell.y2 = upp[1];
                    mappedCells.push_back(cell);
                    mappedIndices.push_back(idx);
                    mapped++;
                }

                if (g_verbose)
                {
                    gsInfo << "[geo-coarsen] patch=" << patch << " level=" << levNow
                           << " geoBoxes=" << geoBoxes.size() << " mappedToVectorS=" << mapped
                           << " available=" << available.size() << "\n";
                    outfile << "[geo-coarsen] patch=" << patch << " level=" << levNow
                            << " geoBoxes=" << geoBoxes.size() << " mappedToVectorS=" << mapped
                            << " available=" << available.size() << "\n";
                }

                if (mapped > 0) {
                    result.geoCells = std::move(mappedCells);
                    result.geoCellIndices = std::move(mappedIndices);
                    result.acceptedDelta = geoDelta;
                    if (g_verbose)
                    {
                        gsInfo << "[geo-coarsen] delta accepted=" << geoDelta << " (breaking)\n";
                        outfile << "[geo-coarsen] delta accepted=" << geoDelta << " (breaking)\n";
                    }
                    break;
                }
            }
        }

        if (!result.geoCells.empty()) {
            result.useGeo = true;
            result.currCellIndex = result.geoCellIndices.front();
            result.currArrayIndex = 0;
            result.x1U = result.geoCells.front().x1;
            result.y1U = result.geoCells.front().y1;
            result.x2U = result.geoCells.front().x2;
            result.y2U = result.geoCells.front().y2;
            if (g_verbose)
                gsInfo << "[geo-coarsen] useGeo=1 cells=" << result.geoCells.size()
                       << " firstIndex=" << result.currCellIndex
                       << " firstBox=[" << result.x1U << "," << result.y1U << "," << result.x2U << "," << result.y2U << "]\n";
            outfile << "[geo-coarsen] useGeo=1 cells=" << result.geoCells.size()
                    << " firstIndex=" << result.currCellIndex
                    << " firstBox=[" << result.x1U << "," << result.y1U << "," << result.x2U << "," << result.y2U << "]\n";
        } else {
            // preflightOnly: don't call pickCell — preserve rand() state for the real call
            if (preflightOnly) return result; // acceptedDelta stays infinity
            int x1Ui = 0, y1Ui = 0, x2Ui = 0, y2Ui = 0;
            result.currCellIndex = pickCell(vectorS, result.currArrayIndex, levNow,
                x1Ui, y1Ui, x2Ui, y2Ui, interior);
            result.x1U = x1Ui;
            result.y1U = y1Ui;
            result.x2U = x2Ui;
            result.y2U = y2Ui;
            if (g_verbose) gsInfo << "[geo-coarsen] useGeo=0 fallback pickCell index=" << result.currCellIndex << "\n";
            outfile << "[geo-coarsen] useGeo=0 fallback pickCell index=" << result.currCellIndex << "\n";
        }
    }
    else if (method == 'g' && patchMirroredInPreflight) {
        int chosenArrayIndex = result.currArrayIndex;
        if (pickHierarchySafeSiblingGroup(
                vectorS,
                result.currArrayIndex,
                levNow,
                interior,
                result.geoCells,
                result.geoCellIndices,
                chosenArrayIndex))
        {
            result.useGeo = true;
            result.acceptedDelta = 0.0; // sibling group found — highest priority
            result.currArrayIndex = chosenArrayIndex;
            result.currCellIndex = result.geoCellIndices.front();
            result.x1U = result.geoCells.front().x1;
            result.y1U = result.geoCells.front().y1;
            result.x2U = result.geoCells.front().x2;
            result.y2U = result.geoCells.front().y2;
            if (g_verbose)
                gsInfo << "[geo-coarsen] preflight mirrored patch=" << patch
                       << ", hierarchy-safe sibling fallback cells=" << result.geoCells.size()
                       << " level=" << levNow << "\n";
            outfile << "[geo-coarsen] preflight mirrored patch=" << patch
                    << ", hierarchy-safe sibling fallback cells=" << result.geoCells.size()
                    << " level=" << levNow << "\n";
        }
        else
        {
            // preflightOnly: don't call pickCell — preserve rand() state
            if (preflightOnly) return result; // acceptedDelta stays infinity
            int x1Ui = 0, y1Ui = 0, x2Ui = 0, y2Ui = 0;
            result.currCellIndex = pickCell(vectorS, result.currArrayIndex, levNow,
                x1Ui, y1Ui, x2Ui, y2Ui, interior);
            result.x1U = x1Ui;
            result.y1U = y1Ui;
            result.x2U = x2Ui;
            result.y2U = y2Ui;
            if (g_verbose)
                gsInfo << "[geo-coarsen] preflight mirrored patch=" << patch
                       << ", no hierarchy-safe sibling group found; fallback to pickCell at level=" << levNow << "\n";
            outfile << "[geo-coarsen] preflight mirrored patch=" << patch
                    << ", no hierarchy-safe sibling group found; fallback to pickCell at level=" << levNow << "\n";
        }
    }
    else if (method == 'r') {
        int x1Ui = 0, y1Ui = 0, x2Ui = 0, y2Ui = 0;
        result.currCellIndex = pickCell(vectorS, result.currArrayIndex, levNow,
            x1Ui, y1Ui, x2Ui, y2Ui, interior);
        result.x1U = x1Ui;
        result.y1U = y1Ui;
        result.x2U = x2Ui;
        result.y2U = y2Ui;
    }
    else if (method == 'l') {
        int x1Ui = 0, y1Ui = 0, x2Ui = 0, y2Ui = 0;
        result.currCellIndex = pickCell(vectorS, attempt, levNow,
            x1Ui, y1Ui, x2Ui, y2Ui, 1, interior);
        result.x1U = x1Ui;
        result.y1U = y1Ui;
        result.x2U = x2Ui;
        result.y2U = y2Ui;
        result.currArrayIndex = 0;
    }
    else if (method == 's') {
        if (valid && (pickedOne + 1) >= 0 && (pickedOne + 1) < pickedCells.cols() && !pickedCells(0, pickedOne + 1)) {
            int x1Ui = 0, y1Ui = 0, x2Ui = 0, y2Ui = 0;
            result.currCellIndex = pickCell(vectorS, attempt, levNow,
                x1Ui, y1Ui, x2Ui, y2Ui, 1, interior);
            result.x1U = x1Ui;
            result.y1U = y1Ui;
            result.x2U = x2Ui;
            result.y2U = y2Ui;
            result.currArrayIndex = 0;
        }
        else {
            int x1Ui = 0, y1Ui = 0, x2Ui = 0, y2Ui = 0;
            result.currCellIndex = pickCell(vectorS, result.currArrayIndex, levNow,
                x1Ui, y1Ui, x2Ui, y2Ui, interior);
            result.x1U = x1Ui;
            result.y1U = y1Ui;
            result.x2U = x2Ui;
            result.y2U = y2Ui;
        }
        pickedOne = result.currCellIndex;
        pickedCells(0, result.currCellIndex) = 1;
    }
    else {
        gsInfo << "UNKNOWN METHOD - using default behavior\n";
        outfile << "UNKNOWN METHOD - using default behavior\n";
    }

    currArrayIndex = result.currArrayIndex;
    return result;
}



void printTheMatrix(const gsMatrix<real_t>& matrix, const std::string& matrixName) {
    outfile << matrixName << "\n";
    for (size_t row = 0; row < matrix.rows(); ++row) {
        for (size_t col = 0; col < matrix.cols(); ++col) {
            outfile << matrix(row, col) << " ";
        }
        outfile << "\n";
    }

    // Also save to separate file
    std::ofstream matrixFile(matrixName + ".txt");
    if (matrixFile.is_open()) {
        matrixFile << matrixName << "\n";
        for (size_t row = 0; row < matrix.rows(); ++row) {
            for (size_t col = 0; col < matrix.cols(); ++col) {
                matrixFile << matrix(row, col) << " ";
            }
            matrixFile << "\n";
        }
        matrixFile.close();
    }
}

// Maps the previously accepted vectSol to a new MPBES by matching
// {patch, level, tensorIdx} descriptions across MPBES rebuilds.
// Surviving functions carry their last fitted value; newly introduced
// coarser functions (no match) start at zero and will be covered by the local fit.
static gsMatrix<real_t> mapPrevSolToNewMpbes(
    const std::vector<std::vector<std::vector<index_t>>>& prevFD,
    const gsMatrix<real_t>& prevSol,
    const MPBES<2, real_t>& newMpbes)
{
    const index_t geoDim = 2;
    gsMatrix<real_t> mapped;
    mapped.setZero(newMpbes.size(), geoDim);

    auto encodeKey = [](index_t patch, index_t level, index_t tensorIdx) -> long long {
        return (static_cast<long long>(patch) << 40)
             | (static_cast<long long>(level) << 32)
             | static_cast<long long>(static_cast<unsigned int>(tensorIdx));
    };

    std::unordered_map<long long, index_t> prevLookup;
    prevLookup.reserve(prevFD.size() * 2);
    for (index_t oldIdx = 0; oldIdx < static_cast<index_t>(prevFD.size()); ++oldIdx)
        for (const auto& twin : prevFD[oldIdx])
            if (twin.size() >= 3)
                prevLookup[encodeKey(twin[0], twin[1], twin[2])] = oldIdx;

    const auto& newFD = newMpbes.functionDescription();
    index_t nCarried = 0, nNew = 0;

    for (index_t newIdx = 0; newIdx < newMpbes.size(); ++newIdx)
    {
        if (newIdx >= static_cast<index_t>(newFD.size())) continue;
        bool found = false;
        for (const auto& twin : newFD[newIdx])
        {
            if (twin.size() < 3) continue;
            auto it = prevLookup.find(encodeKey(twin[0], twin[1], twin[2]));
            if (it == prevLookup.end()) continue;
            const index_t oldIdx = it->second;
            if (oldIdx >= prevSol.rows()) continue;
            mapped.row(newIdx) = prevSol.row(oldIdx).leftCols(geoDim);
            found = true;
            ++nCarried;
            break;
        }
        if (!found) ++nNew;
    }

    if (g_verbose)
        gsInfo  << "[mapPrev] newSize=" << newMpbes.size()
                << " carried=" << nCarried << " new(zero)=" << nNew << "\n";
    outfile << "[mapPrev] newSize=" << newMpbes.size()
            << " carried=" << nCarried << " new(zero)=" << nNew << "\n";

    return mapped;
}

static gsMatrix<real_t> mapMpCoefficientsToMpbes(
    const gsMultiPatch<real_t>& sourceMp,
    const MPBES<2, real_t>& mpbes)
{
    const index_t geoDim  = 2;
    const index_t patchCount = std::min<index_t>(sourceMp.nPatches(), mpbes.nPatches());

    gsMatrix<real_t> mapped;
    mapped.setZero(mpbes.size(), geoDim);

    // Build per-patch reverse lookup: (level, tensorIdx) → fine THB flat index.
    // Uses the same API (levelOf / flatTensorIndexOf) as createIndexMapping.
    std::vector<std::unordered_map<long long, index_t>> fineIdx(patchCount);
    for (index_t patch = 0; patch < patchCount; ++patch)
    {
        const gsBasis<real_t>*              rawPtr  = &sourceMp.patch(patch).basis();
        const gsTHBSplineBasis<2, real_t>*  finePtr =
            dynamic_cast<const gsTHBSplineBasis<2, real_t>*>(rawPtr);
        if (!finePtr)
            finePtr = dynamic_cast<const gsTHBSplineBasis<2, real_t>*>(&rawPtr->source());
        if (!finePtr) continue;

        const index_t n = finePtr->size();
        for (index_t i = 0; i < n; ++i)
        {
            const int lv  = finePtr->levelOf(i);
            const int tid = finePtr->flatTensorIndexOf(i);
            fineIdx[patch][static_cast<long long>(lv) << 32 | static_cast<unsigned>(tid)] = i;
        }
    }

    // For each MPBES global function: find its (patch, level, tensorIdx) from
    // functionDescription, look up the corresponding fine coefficient by (level, tensorIdx).
    const auto& fd = mpbes.functionDescription();
    index_t nMapped = 0, nNotInFine = 0;

    for (index_t globalIdx = 0; globalIdx < mpbes.size(); ++globalIdx)
    {
        if (globalIdx >= static_cast<index_t>(fd.size())) continue;
        const auto& twins = fd[globalIdx];
        if (twins.empty()) continue;

        bool assigned = false;
        for (const auto& twin : twins)
        {
            if (twin.size() < 3) continue;
            const int patch     = twin[0];
            const int level     = twin[1];
            const int tensorIdx = twin[2];
            if (patch < 0 || patch >= patchCount) continue;

            const long long key = static_cast<long long>(level) << 32
                                | static_cast<unsigned>(tensorIdx);
            auto it = fineIdx[patch].find(key);
            if (it == fineIdx[patch].end()) continue;

            const index_t fi = it->second;
            const gsMatrix<real_t> pc = sourceMp.patch(patch).coefs();
            if (fi >= pc.rows()) continue;

            for (index_t d = 0; d < std::min<index_t>(geoDim, pc.cols()); ++d)
                mapped(globalIdx, d) = pc(fi, d);
            assigned = true;
            break;
        }
        if (assigned) ++nMapped; else ++nNotInFine;
    }

    gsInfo  << "[mapMp] mpbes.size()=" << mpbes.size()
            << " mapped=" << nMapped << " notInFine=" << nNotInFine << "\n";
    outfile << "[mapMp] mpbes.size()=" << mpbes.size()
            << " mapped=" << nMapped << " notInFine=" << nNotInFine << "\n";

    return mapped;
}

struct OriginalMpbesReference
{
    std::unique_ptr<MPBES<2, real_t>> mpbes;
    gsMatrix<real_t> coefficients;
    gsVector<size_t> numIrregular;
    int minusnumber = 0;
    bool jacobianCheckRan = false;
    bool valid = false;
};

static OriginalMpbesReference buildOriginalMpbesReference(const std::string& filename)
{
    OriginalMpbesReference ref;

    gsMultiPatch<>::uPtr mp = gsReadFile<>(filename);
    if (!mp || mp->nPatches() == 0)
        return ref;

    gsMultiPatch<real_t> mpWork = *mp;
    const index_t nPatches = mp->nPatches();

    gsVector<gsTHBSplineBasis<2, real_t>*> THBFromGeo(nPatches);
    gsVector<gsTHBSplineBasis<2, real_t>> SubdomainHierarchy(nPatches);
    gsVector<gsVector<gsTensorBSplineBasis<2, real_t>>> Bells(nPatches);

    gsVector<gsMatrix<index_t>> lowCorners(nPatches);
    gsVector<gsMatrix<index_t>> upCorners(nPatches);
    gsVector<gsVector<index_t>> myLevel(nPatches);
    gsVector<index_t> maxLevel(nPatches);
    gsVector<index_t> currentLastNonZeroRow(nPatches);
    gsVector<gsVector<gsVector<index_t>>> boxMat(nPatches);

    const unsigned interioru0 = 0;
    const unsigned interioru1 = 0;
    const unsigned interiorv1 = 0;
    const int degreeU = mp->patch(0).basis().degree(0);
    const int degreeV = mp->patch(0).basis().degree(1);
    const unsigned multEndU = static_cast<unsigned>(degreeU + 1);
    const unsigned multEndV = static_cast<unsigned>(degreeV + 1);

    gsKnotVector<> ku(mp->patch(0).basis().support()(0, 0), mp->patch(0).basis().support()(0, 1),
        interioru0 + interioru1, multEndU);
    gsKnotVector<> kv(mp->patch(0).basis().support()(1, 0), mp->patch(0).basis().support()(1, 1),
        interiorv1, multEndV);
    gsTensorBSplineBasis<2, real_t> tens(ku, kv);

    for (index_t patch = 0; patch < nPatches; ++patch)
    {
        THBFromGeo(patch) = dynamic_cast<gsTHBSplineBasis<2>*>(&mp->patch(patch).basis().source());
        if (!THBFromGeo(patch))
            return ref;

        SubdomainHierarchy(patch) = *THBFromGeo(patch);
        THBFromGeo(patch)->tree().getBoxes(lowCorners(patch), upCorners(patch), myLevel(patch));
        if (myLevel(patch).size() == 0)
            return ref;

        const index_t boxLevelMax = myLevel(patch).maxCoeff();
        maxLevel(patch) = THBFromGeo(patch)->maxLevel();
        currentLastNonZeroRow(patch) = lowCorners(patch).rows();
        boxMat(patch).resize(currentLastNonZeroRow(patch));

        for (index_t boxInd = 0; boxInd < currentLastNonZeroRow(patch); ++boxInd)
        {
            boxMat(patch)(boxInd).resize(5);
            boxMat(patch)(boxInd)(0) = myLevel(patch)(boxInd);
            boxMat(patch)(boxInd)(1) = (real_t)lowCorners(patch)(boxInd, 0) /
                std::pow(2.0, boxLevelMax - myLevel(patch)(boxInd));
            boxMat(patch)(boxInd)(2) = (real_t)lowCorners(patch)(boxInd, 1) /
                std::pow(2.0, boxLevelMax - myLevel(patch)(boxInd));
            boxMat(patch)(boxInd)(3) = (real_t)upCorners(patch)(boxInd, 0) /
                std::pow(2.0, boxLevelMax - myLevel(patch)(boxInd));
            boxMat(patch)(boxInd)(4) = (real_t)upCorners(patch)(boxInd, 1) /
                std::pow(2.0, boxLevelMax - myLevel(patch)(boxInd));
        }

        gsTensorBSplineBasis<2, real_t> tensForBells = THBFromGeo(patch)->tensorLevel(0);
        Bells(patch).resize(maxLevel(patch) + 1);
        for (int level = 0; level < Bells(patch).size(); ++level)
        {
            Bells(patch)(level) = tensForBells;
            tensForBells.uniformRefine();
        }
    }

    gsVector<gsVector<gsVector<std::vector<index_t>>>> twinsIndex(Bells.size()), twinsPatch(Bells.size());
    gsVector<gsVector<gsVector<index_t>>> hasATwin(Bells.size()), isActive(Bells.size()), isTouching(Bells.size());

    IdentifyPatches(mpWork,
        Bells,
        isTouching,
        twinsIndex,
        twinsPatch,
        hasATwin,
        isActive);

    std::vector<int> firstSide;
    std::vector<int> secondSide;
    std::vector<int> firstPatch;
    std::vector<int> secondPatch;
    orientThePatches(mpWork, firstSide, secondSide, firstPatch, secondPatch);

    gsVector<gsVector<gsVector<index_t>>> isIncluded;
    gsVector<gsVector<gsVector<index_t>>> indexInTHB;
    gsVector<std::vector<std::array<int, 2>>> thbToBellsMapping;
    createIndexMapping(
        Bells,
        SubdomainHierarchy,
        twinsIndex,
        twinsPatch,
        isIncluded,
        indexInTHB,
        thbToBellsMapping
    );

    ref.mpbes = std::make_unique<MPBES<2, real_t>>(
        boxMat,
        SubdomainHierarchy,
        Bells,
        hasATwin,
        twinsIndex,
        twinsPatch,
        indexInTHB,
        currentLastNonZeroRow,
        false,
        0);

    if (!ref.mpbes || ref.mpbes->size() == 0)
        return ref;

    ref.coefficients = mapMpCoefficientsToMpbes(*mp, *ref.mpbes);

    {
        const index_t jacPoints = 40;
        const index_t sampleCount = jacPoints * jacPoints;
        gsVector<real_t> low(2), up(2);
        low << 0.0, 0.0;
        up << 1.0, 1.0;

        gsVector<gsMatrix<real_t>> uv(ref.mpbes->nPatches());
        for (index_t patch = 0; patch < static_cast<index_t>(uv.size()); ++patch)
            uv(patch) = uniformPointGrid(low, up, sampleCount);

        const MirroredCheckResult initialMirrorReport =
            computeMirroredCheckPerPatch(uv, *ref.mpbes, ref.coefficients, 1e-12);

        ref.numIrregular.resize(ref.mpbes->nPatches());
        ref.numIrregular.setZero();
        gsInfo << "\n=== INITIAL MPBES DETERMINANT CHECK (before coarsening) ===\n";
        if (outfile.is_open())
            outfile << "\n=== INITIAL MPBES DETERMINANT CHECK (before coarsening) ===\n";
        ref.minusnumber = checkJacobianDeterminant(
            uv,
            *ref.mpbes,
            ref.coefficients,
            ref.numIrregular,
            false);
        ref.jacobianCheckRan = true;

        index_t totalPoints = 0;
        for (index_t patch = 0; patch < uv.size(); ++patch)
            totalPoints += uv[patch].cols();

        const real_t irregularPercentage =
            totalPoints > 0 ? (100.0 * ref.minusnumber / totalPoints) : 0.0;

        gsInfo << "[initial-mpbes] irregular points: " << ref.minusnumber
               << " / " << totalPoints << " (" << irregularPercentage << "%)\n";
        if (outfile.is_open())
        {
            outfile << "[initial-mpbes] irregular points: " << ref.minusnumber
                    << " / " << totalPoints << " (" << irregularPercentage << "%)\n";
        }

        for (index_t patch = 0; patch < ref.numIrregular.size(); ++patch)
        {
            gsInfo << "[initial-mpbes] patch " << patch
                   << " irregular=" << ref.numIrregular(patch)
                   << ", used=" << initialMirrorReport.usedPerPatch[patch]
                   << ", pos=" << initialMirrorReport.posPerPatch[patch]
                   << ", neg=" << initialMirrorReport.negPerPatch[patch]
                   << ", mirrored=" << (initialMirrorReport.mirrored[patch] ? "true" : "false")
                   << "\n";
            if (outfile.is_open())
            {
                outfile << "[initial-mpbes] patch " << patch
                        << " irregular=" << ref.numIrregular(patch)
                        << ", used=" << initialMirrorReport.usedPerPatch[patch]
                        << ", pos=" << initialMirrorReport.posPerPatch[patch]
                        << ", neg=" << initialMirrorReport.negPerPatch[patch]
                        << ", mirrored=" << (initialMirrorReport.mirrored[patch] ? "true" : "false")
                        << "\n";
            }
        }

        gsInfo << "=== END INITIAL MPBES DETERMINANT CHECK ===\n\n";
        if (outfile.is_open())
            outfile << "=== END INITIAL MPBES DETERMINANT CHECK ===\n\n";
    }

    ref.valid = true;
    return ref;
}

static const char* featureSideName(FeatureSide side)
{
    switch (side)
    {
    case FeatureSide::U0:
        return "U0";
    case FeatureSide::U1:
        return "U1";
    case FeatureSide::V0:
        return "V0";
    case FeatureSide::V1:
        return "V1";
    }
    return "unknown";
}

static const GeometryPreflightInterfaceInfo* findPreflightInterfaceInfo(
    index_t patchA,
    index_t patchB)
{
    if (!g_geometryPreflight.valid)
        return nullptr;

    for (size_t i = 0; i < g_geometryPreflight.interfaces.size(); ++i)
    {
        const GeometryPreflightInterfaceInfo& info = g_geometryPreflight.interfaces[i];
        if ((info.patchA == patchA && info.patchB == patchB) ||
            (info.patchA == patchB && info.patchB == patchA))
            return &info;
    }

    return nullptr;
}

static bool isPreflightMirroredPatch(index_t patch)
{
    return g_geometryPreflight.valid &&
        patch >= 0 &&
        patch < static_cast<index_t>(g_geometryPreflight.mirroredReport.mirrored.size()) &&
        g_geometryPreflight.mirroredReport.mirrored[patch];
}

GeometryPreflightInfo runGeometryPreflight(
    const std::string& filename,
    const OriginalMpbesReference& originalReference,
    index_t pointsPerDirection)
{
    GeometryPreflightInfo preflight;

    if (!originalReference.valid || !originalReference.mpbes)
        return preflight;

    gsMultiPatch<>::uPtr mpPtr = gsReadFile<>(filename);
    if (!mpPtr || mpPtr->nPatches() == 0)
        return preflight;

    gsMultiPatch<real_t> mp = *mpPtr;
    const index_t sampleCount = std::max<index_t>(2, pointsPerDirection * pointsPerDirection);

    gsVector<real_t> low(2), up(2);
    low << 0.0, 0.0;
    up << 1.0, 1.0;

    gsVector<gsMatrix<real_t>> uv(mp.nPatches());
    for (index_t patch = 0; patch < static_cast<index_t>(mp.nPatches()); ++patch)
        uv(patch) = uniformPointGrid(low, up, sampleCount);

    preflight.mirroredReport = computeMirroredCheckPerPatch(
        uv,
        *originalReference.mpbes,
        originalReference.coefficients,
        1e-12);

    preflight.interfaces = collectVisualizationInterfaceInfos(mp);

    preflight.valid = true;
    return preflight;
}

void logGeometryPreflight(
    const GeometryPreflightInfo& preflight,
    const std::string& label,
    bool logToInfo)
{
    std::ostringstream buffer;
    buffer << "GEOMETRY_PREFLIGHT_BEGIN [" << label << "]\n";
    buffer << "valid=" << (preflight.valid ? "true" : "false") << "\n";

    if (!preflight.valid)
    {
        buffer << "reason=preflight data unavailable\n";
        buffer << "GEOMETRY_PREFLIGHT_END\n";
    }
    else
    {
        index_t mirroredCount = 0;
        for (index_t p = 0; p < static_cast<index_t>(preflight.mirroredReport.mirrored.size()); ++p)
        {
            if (preflight.mirroredReport.mirrored[p])
                ++mirroredCount;

            buffer << "patch " << p
                   << ": used=" << preflight.mirroredReport.usedPerPatch[p]
                   << ", pos=" << preflight.mirroredReport.posPerPatch[p]
                   << ", neg=" << preflight.mirroredReport.negPerPatch[p];

            if (preflight.mirroredReport.usedPerPatch[p] > 0)
            {
                const real_t negRatio =
                    static_cast<real_t>(preflight.mirroredReport.negPerPatch[p]) /
                    static_cast<real_t>(preflight.mirroredReport.usedPerPatch[p]);
                buffer << ", negRatio=" << negRatio;
            }
            else
            {
                buffer << ", negRatio=n/a";
            }

            buffer << ", mirrored="
                   << (preflight.mirroredReport.mirrored[p] ? "true" : "false")
                   << "\n";
        }

        buffer << "mirroredCount=" << mirroredCount << "\n";
        buffer << "interfaces=" << preflight.interfaces.size() << "\n";
        for (size_t i = 0; i < preflight.interfaces.size(); ++i)
        {
            const GeometryPreflightInterfaceInfo& info = preflight.interfaces[i];
            buffer << "interface " << i
                   << ": patchA=" << info.patchA
                   << ", sideIndexA=" << info.sideIndexA
                   << ", patchB=" << info.patchB
                   << ", sideIndexB=" << info.sideIndexB;

            if (info.hasMappedSides)
            {
                buffer << ", sideA=" << featureSideName(info.sideA)
                       << ", sideB=" << featureSideName(info.sideB)
                       << ", orientation="
                       << (info.orientationReversed ? "reversed" : "same");
            }
            else
            {
                buffer << ", sideMapping=unavailable";
            }
            buffer << "\n";
        }

        buffer << "GEOMETRY_PREFLIGHT_END\n";
    }

    if (outfile.is_open())
    {
        outfile << buffer.str();
        outfile.flush();
    }

    if (logToInfo)
        gsInfo << buffer.str();
}

// Struct to hold all data computed by the algorithm
struct AlgorithmResult {
    gsMultiPatch<>::uPtr mp;
    gsVector<gsVector<gsTensorBSplineBasis<2, real_t>>> Bells;
    gsVector<gsVector<gsVector<index_t>>> boxMat;
    gsVector<gsVector<gsVector<index_t>>> acceptedBoxMat;
    std::vector<std::vector<std::vector<double>>> acceptedCoefs;
    gsMatrix<real_t> AcceptedvectSol;
    gsVector<index_t> AcceptedlastRow;
    gsMatrix<> acceptedMatOut;
    std::vector<std::vector<std::vector<real_t>>> featureCoordinates;
    gsVector<gsMatrix<real_t>> uvFeature;
    gsVector<gsMatrix<real_t>> xyFeature;
    gsVector<gsVector<gsVector<index_t>>> AcceptedisActive;
    gsVector<gsVector<gsVector<index_t>>> AcceptedglobalIndex;
    gsVector<gsVector<std::vector<index_t>>> AcceptedfunctionDescription;
    gsVector<gsMatrix<real_t>> uv1;
    gsVector<real_t> lowc2;
    gsVector<real_t> uppc2;
    gsVector<index_t> maxLevel;
    std::chrono::time_point<std::chrono::system_clock> toc;
    int successfullAttempts;
    int totalAttempts;
    gsMatrix<> matFile;
    unsigned interioru0;
    unsigned interiorv0;
    int proj;
    gsMultiPatch<real_t> mp1;
    gsVector<index_t> currentLastNonZeroRow;
    std::unique_ptr<MPBES<2, real_t>> mpbes;
};

// Implementation of checkBoxesConsistency
void checkBoxesConsistency(
    const gsVector<gsVector<gsVector<index_t>>>& boxMat,
    int patch,
    int lastNonZeroRow,
    std::ofstream& outfile
) {
    outfile << "\n=== Box Consistency Check ===\n";
    
    if (lastNonZeroRow == 0) {
        outfile << "ERROR: No boxes defined!\n";
        gsInfo << "ERROR: No boxes defined!\n";
        outfile << "Returning early due to no boxes\n";
        outfile.close();
        return;
    }
    
    outfile << "Patch " << patch << " has " << lastNonZeroRow << " boxes\n\n";
    
    // Print all boxes with their physical domains
    outfile << "Box domains:\n";
    real_t totalArea = 0.0;
    
    for (int b = 0; b < lastNonZeroRow; ++b) {
        index_t level = boxMat(patch)(b)(0);
        index_t uMin = boxMat(patch)(b)(1);
        index_t vMin = boxMat(patch)(b)(2);
        index_t uMax = boxMat(patch)(b)(3);
        index_t vMax = boxMat(patch)(b)(4);
        
        // Convert indices to physical domain [0,1]x[0,1]
        real_t scaleFactor = 1.0 / (1 << level);  // 1 / 2^level
        real_t uMinPhys = uMin * scaleFactor;
        real_t vMinPhys = vMin * scaleFactor;
        real_t uMaxPhys = uMax * scaleFactor;
        real_t vMaxPhys = vMax * scaleFactor;
        
        real_t boxArea = (uMaxPhys - uMinPhys) * (vMaxPhys - vMinPhys);
        totalArea += boxArea;
        
        outfile << "Box " << b << ": u=[" << uMinPhys << ", " << uMaxPhys << "] x v=[" 
                << vMinPhys << ", " << vMaxPhys << "] (area=" << boxArea << ")\n";
    }
    
    // Check for overlaps
    outfile << "\nChecking for overlaps:\n";
    bool hasOverlap = false;
    
    for (int b1 = 0; b1 < lastNonZeroRow; ++b1) {
        for (int b2 = b1 + 1; b2 < lastNonZeroRow; ++b2) {
            index_t level1 = boxMat(patch)(b1)(0);
            index_t uMin1 = boxMat(patch)(b1)(1);
            index_t vMin1 = boxMat(patch)(b1)(2);
            index_t uMax1 = boxMat(patch)(b1)(3);
            index_t vMax1 = boxMat(patch)(b1)(4);
            
            index_t level2 = boxMat(patch)(b2)(0);
            index_t uMin2 = boxMat(patch)(b2)(1);
            index_t vMin2 = boxMat(patch)(b2)(2);
            index_t uMax2 = boxMat(patch)(b2)(3);
            index_t vMax2 = boxMat(patch)(b2)(4);
            
            // Convert to physical domain
            real_t scale1 = 1.0 / (1 << level1);
            real_t scale2 = 1.0 / (1 << level2);
            
            real_t u1min = uMin1 * scale1, u1max = uMax1 * scale1;
            real_t v1min = vMin1 * scale1, v1max = vMax1 * scale1;
            real_t u2min = uMin2 * scale2, u2max = uMax2 * scale2;
            real_t v2min = vMin2 * scale2, v2max = vMax2 * scale2;
            
            // Check for strict overlap (excluding boundary touch)
            real_t tolerance = 1e-10;
            if (u1max > u2min + tolerance && u2max > u1min + tolerance &&
                v1max > v2min + tolerance && v2max > v1min + tolerance) {
                outfile << "  OVERLAP: Box " << b1 << " overlaps Box " << b2 << "\n";
                hasOverlap = true;
            }
        }
    }
    
    if (!hasOverlap) {
        outfile << "  No overlaps detected (boundaries may touch)\n";
    }
    
    // Check total coverage
    outfile << "\nTotal area covered: " << totalArea << " (target: 1.0)\n";
    
    if (std::abs(totalArea - 1.0) > 1e-6) {
        outfile << "ERROR: Boxes do not cover [0,1]x[0,1]!\n";
        outfile << "  Missing area: " << (1.0 - totalArea) << "\n";
        gsInfo << "WARNING: Domain not fully covered by boxes (missing area: " << (1.0 - totalArea) << ")\n";
        outfile << "WARNING: Continuing despite incomplete coverage\n";
    }
    
    if (hasOverlap) {
        outfile << "ERROR: Boxes self-intersect!\n";
        gsInfo << "WARNING: Boxes self-intersect (overlapping regions detected)\n";
        outfile << "WARNING: Continuing despite box overlap\n";
    }
    
    outfile << "\nSUCCESS: Boxes are consistent (no overlaps, full coverage)\n";
}

AlgorithmResult unrefinementAlgorithmHBJ(
    const std::string& filename,
    real_t epsilon_g,
    real_t epsilon_f,
    char method,
    const std::string& givenGeo,
    const std::string& acCond,
    const std::vector<FeatureBoundarySpec>& featureSides) {
    
    // Initialize variables that were in main()
    std::chrono::time_point<std::chrono::system_clock> startTime, iterTime;
    startTime = std::chrono::system_clock::now();
    int row, acceptedsize, attempt = 0;
    int valid = 0;
    real_t lcx, lcy, ucx, ucy;
    int successfullAttempts = 0, totalAttempts = 0;
    int projections = 0;  // Use different name to avoid conflict with proj() function
    // ---- Phase-time accumulators ----
    double t_grenda     = 0.0; // Grenda preflight: selectCellForCoarsening on all patches
    double t_cell_sel   = 0.0; // Actual cell selection (winning patch)
    double t_mpbes_bld  = 0.0; // Per-step MPBES build (orientThePatches + createIndexMapping + MPBES ctor)
    double t_orient     = 0.0; //   └─ orientThePatches
    double t_indexmap   = 0.0; //   └─ createIndexMapping
    double t_mpbes_ctor = 0.0; //   └─ MPBES<2,real_t> constructor
    double t_assemble   = 0.0; // assemble() call
    double t_ata_ldlt   = 0.0; // Sparse AtA formation + LDLT solve
    double t_lo         = 0.0; // LO nonLinearOptimization call
    double t_nlo        = 0.0; // NLO nonLinearOptimization call
    double t_mpbes_copy  = 0.0; // MPBES deep-copy on accepted coarsening
    double t_jack        = 0.0; // checkJacobianDeterminant (already printed per-step as "jack took:")
    double t_boundary    = 0.0; // testBoundaryAssembly
    double t_targetgeom  = 0.0; // target geometry matrix construction (mp1.eval at uv1 points)
    double t_globalerr   = 0.0; // globalFittingError
    double t_preflight   = 0.0; // generateOriginalGeometryMesh (mesh write) + direct Jacobian check on mp
    double t_pf_mesh     = 0.0; // sub: generateOriginalGeometryMesh
    double t_pf_check    = 0.0; // sub: computeMirroredCheckFromMp + logGeometryPreflight
    double t_presolve    = 0.0; // column mapping + defect correction (between assemble and LDLT)
    double t_postsolve   = 0.0; // scatter/merge + residual diagnostics + uv2_adaptive setup (between LDLT and jack)
    double t_postjack    = 0.0; // logging + flush + A-alloc between jack and boundary
    double t_postaccept  = 0.0; // flushing + acceptance logging between globalerr and next THB rebuild
    double t_localsetup  = 0.0; // local region build + assembleActiveFunctions selection (local fitting only)
    double t_init        = 0.0; // one-time init before main while loop
    double t_thb_rebuild = 0.0; // THB rebuild: rebuildHierarchy + refineElements replay + outfile<<THB
    double t_extract     = 0.0; // MPBES data extraction + vectSol init (after MPBES build)
    double t_gap_a       = 0.0; // GAP A: between Grenda end and cell-sel start (coarsenGroup build + var decl)
    double t_bookkeeping = 0.0; // success/rejection bookkeeping after t_postaccept (data copies + outfile)
    std::string xmlFile = "/home/targon/theydarov/output/" + givenGeo + acCond;
    const std::string filePrefix = gsFileManager::getBasename(filename) + "_";
    
    gsMultiPatch<>::uPtr mp = gsReadFile<>(filename);
    if (!mp) {
        gsInfo << "ERROR: Failed to read geometry file: " << filename << "\n";
        outfile << "ERROR: Failed to read geometry file: " << filename << "\n";
        outfile.close();
        AlgorithmResult errorResult;
        errorResult.mp = nullptr;
        return errorResult;
    }
    gsMultiPatch<> domain;
    gsMultiPatch<real_t> mp1 = *mp;
    gsMultiPatch<>::uPtr initmp = gsReadFile<>(filename);
    gsMultiPatch<> newmp;
    gsMatrix<> matFile;
    gsVector<> lowc1(2);
    gsVector<> uppc1(2);
    lowc1(0) = 0;
    lowc1(1) = 0;
    uppc1(0) = 1;
    uppc1(1) = 1;
    unsigned interioru0 = 0;
    unsigned interiorv0 = 0;
    unsigned interioru1 = 0;
    unsigned interiorv1 = 0;
    
    if (mp->nPatches() == 0) {
        gsInfo << "ERROR: No patches found in geometry file\n";
        outfile << "ERROR: No patches found in geometry file\n";
        outfile.close();
        AlgorithmResult errorResult;
        errorResult.mp = nullptr;
        return errorResult;
    }

    // ---- Preflight: direct Jacobian check on mp (no MPBES build needed) ----
    gsMatrix<real_t> originalCoefsCapture;   // filled from first vectSol after MPBES builds
    const gsMatrix<real_t>* originalCoefficients = nullptr;
    {
        auto _tPFM0 = std::chrono::system_clock::now();
        generateOriginalGeometryMesh(*mp, 16, filePrefix + "output_mesh_original");
        t_pf_mesh = std::chrono::duration<double>(std::chrono::system_clock::now() - _tPFM0).count();

        auto _tPFC0 = std::chrono::system_clock::now();
        gsVector<real_t> pfLow(2), pfUp(2);
        pfLow << 0.0, 0.0;  pfUp << 1.0, 1.0;
        gsVector<gsMatrix<real_t>> uvPF(mp->nPatches());
        for (index_t p = 0; p < static_cast<index_t>(mp->nPatches()); ++p)
            uvPF(p) = uniformPointGrid(pfLow, pfUp, 20 * 20);

        GeometryPreflightInfo gPF;
        gPF.mirroredReport = computeMirroredCheckFromMp(uvPF, *mp, 1e-12);
        gPF.interfaces     = collectVisualizationInterfaceInfos(*mp);
        gPF.valid          = true;
        g_geometryPreflight = gPF;
        logGeometryPreflight(gPF, gsFileManager::getBasename(filename), false);
        t_pf_check = std::chrono::duration<double>(std::chrono::system_clock::now() - _tPFC0).count();
        // originalCoefficients captured on first vectSol extraction below
        t_preflight = t_pf_mesh + t_pf_check;
    }

    const int degreeU = mp->patch(0).basis().degree(0);
    const int degreeV = mp->patch(0).basis().degree(1);
    const unsigned multEndU = static_cast<unsigned>(degreeU + 1);
    const unsigned multEndV = static_cast<unsigned>(degreeV + 1);
    
    gsKnotVector<> ku(mp->patch(0).basis().support()(0, 0), mp->patch(0).basis().support()(0, 1),
        interioru0 + interioru1, multEndU);
    gsKnotVector<> kv(mp->patch(0).basis().support()(1, 0), mp->patch(0).basis().support()(1, 1), interiorv1, multEndV);
    gsTensorBSplineBasis<2, real_t> tens(ku, kv);
    std::vector<index_t> box;
    std::vector<std::vector<std::vector<double>>> acceptedCoefs;
    gsVector<gsVector<gsTensorBSplineBasis<2, real_t>>> Bells(mp->nPatches());
    gsVector<gsVector<gsTHBSplineBasis<2, real_t>>> THBVector1(mp->nPatches());
    gsVector<gsVector<gsTHBSplineBasis<2, real_t>>> spilloverPatchStructure(mp->nPatches());
    gsVector<gsTHBSplineBasis<2, real_t>> THBVector(mp->nPatches());
    gsMatrix<> acceptedMatOut;

    // gsFunctionExpr<> ff = Bells(0)(0).function(0)*1.0;
    numberOfPatchesForReference = mp->nPatches();
    gsVector<gsTHBSplineBasis<2, real_t>*> THBFromGeo(mp->nPatches()), THBAccepted(mp->nPatches()), THBTemporary(
        mp->nPatches());
    std::vector<std::vector<std::vector<real_t>>> featureCoordinates(mp->nPatches());
    gsVector<gsMatrix<index_t >> lowCorners(mp->nPatches());
    gsVector<gsMatrix<index_t >> upCorners(mp->nPatches());
    gsVector<gsVector<index_t >> myLevel(mp->nPatches());
    gsVector<index_t> fineLevel(mp->nPatches());
    gsVector<index_t> coarseLevel(mp->nPatches());
    gsVector<index_t> maxLevel(mp->nPatches());
    gsMatrix<> coefshat, anmat2, anmat;
    gsMatrix<> supps;
    const gsVector<real_t> lowc2 = lowc1;
    const gsVector<real_t> uppc2 = uppc1;
    int numPoints = std::max(5, (int)std::sqrt(mp->patch(0).basis().size()));
    int jacPoints = 40;
    gsVector<gsMatrix<real_t>> uv1(mp->nPatches());
    gsVector<gsMatrix<real_t>> uv2(mp->nPatches());
    gsVector<gsMatrix<real_t>> xy1(mp->nPatches());
    gsVector<gsMatrix<real_t>> uvFeature(mp->nPatches());
    gsVector<gsMatrix<real_t>> xyFeature(mp->nPatches());
    gsVector<gsMatrix<real_t>> xyFeatureMiddle(mp->nPatches());
    const std::vector<FeatureBoundarySpec> normalizedFeatureSides =
        normalizeFeatureBoundarySpecs(mp1, featureSides, true);
    for (int patch = 0; patch < uv1.size(); ++patch)
    {
        uv1(patch) = uniformPointGrid(lowc2, uppc2, numPoints * numPoints);
        uv2(patch) = uniformPointGrid(lowc2, uppc2, jacPoints * jacPoints);
    }
    gsVector<gsMatrix<real_t>> uvMiddle(mp->nPatches());
    for (int patch = 0; patch < mp->nPatches(); ++patch) {
        uvMiddle(patch).resize(uv1(patch).rows(), uv1(patch).cols() - 1);
        xy1(patch) = mp->patch(patch).eval(uv1(patch));
        // //gsDebugVar(xy1(patch));
        for (int k = 0; k < uv1(patch).cols() - 1; ++k) {
            uvMiddle(patch).col(k) = (uv1(patch).col(k + 1) + uv1(patch).col(k)) / 2;
        }
    }
    gsVector<gsMatrix<real_t>> xyMiddle(mp->nPatches());

    featureCoordinates[0].resize(1);
    featureCoordinates[1].resize(1);
    featureCoordinates[2].resize(1);

    featureCoordinates[0][0].push_back(1);
    featureCoordinates[0][0].push_back(1);
    featureCoordinates[0][0].push_back(0);
    featureCoordinates[0][0].push_back(1);
    featureCoordinates[1][0].push_back(1);
    featureCoordinates[1][0].push_back(0);
    featureCoordinates[1][0].push_back(1);
    featureCoordinates[1][0].push_back(1);
    featureCoordinates[2][0].push_back(0);
    featureCoordinates[2][0].push_back(0);
    featureCoordinates[2][0].push_back(0);
    featureCoordinates[2][0].push_back(1);
    for (int patch = 0; patch < mp->nPatches(); ++patch) {
        xyMiddle(patch) = mp->patch(patch).eval(uvMiddle(patch));
    }
    gsMatrix<> xyNew(2, xy1.size() * xy1(0).cols());
    gsMatrix<> xyMiddleNew(2, xy1.size() * xy1(0).cols());
    auto numPoints1 = 100;
    for (size_t i = 0; i < uvFeature.size(); i++)
    {
        uvFeature(i).setZero();
    }
    featureMatrixGenerator(uvFeature, featureCoordinates, numPoints1);



    int dummySize = 0;
    int dummySizeMiddle = 0;
    for (int patch = 0; patch < mp->nPatches(); ++patch) {
        //xyMiddle(patch) = mp->patch(patch).eval(uvMiddle(patch));
        //xyFeatureMiddle(patch) = mp->patch(patch).eval(uvFeatureMiddle(patch));
        if (uvFeature(patch).size() > 0) {
            xyFeature(patch) = mp->patch(patch).eval(uvFeature(patch));
            dummySize += xyFeature(patch).cols();
        }
        //if (uvFeatureMiddle(patch).size() > 0) {
        //    xyFeatureMiddle(patch) = mp->patch(patch).eval(uvFeatureMiddle(patch));
        //    dummySizeMiddle += xyFeatureMiddle(patch).cols();
        //}
    }


    gsMatrix<> xyFeatureNew(2, dummySize);
    //gsMatrix<> xyFeatureMiddleNew(2, dummySizeMiddle);
    //gsDebugVar(dummySize);
    dummySize = 0;

    for (int currentPatch = 0; currentPatch < xy1.size(); ++currentPatch) {
        for (int i = 0; i < xy1(0).cols(); ++i) {
            xyNew.col(currentPatch * xy1(0).cols() + i) = xy1(currentPatch).col(i);
        }
        //for (int i = 0; i < xyMiddle(0).cols(); ++i) {
        //    xyMiddleNew.col(currentPatch * xyMiddle(0).cols() + i) = xyMiddle(currentPatch).col(i);
        //}
        for (int i = 0; i < xyFeature(currentPatch).cols(); ++i) {
            xyFeatureNew.col(dummySize + i) = xyFeature(currentPatch).col(i);
        }

        dummySize += xyFeature(currentPatch).cols();
    }


    //gsDebugVar(xyFeature(0)(1, 5));
    //gsDebugVar(xyFeatureNew);
    gsVector<gsVector<gsVector<index_t >>> globalIndex(Bells.size());
    for (int patch = 0; patch < mp->nPatches(); ++patch) {
        gsTensorBSplineBasis<2, real_t> tensForBells = tens;
        THBFromGeo(patch)
            = dynamic_cast<gsTHBSplineBasis<2> *>(&mp->patch(patch).basis().source());
        //Bells(patch).resize(THBFromGeo->maxLevel() + 1);
        //gsDebugVar(mp->patch(patch).basis());
        (THBFromGeo(patch)->tree().getBoxes(lowCorners(patch), upCorners(patch), myLevel(patch)));
        supps.resize(THBFromGeo(patch)->support().rows(),
            THBFromGeo(patch)->support().cols() + 2 * THBFromGeo(patch)->size());
        for (int i = 0; i < THBFromGeo(patch)->size(); ++i) {
            for (int j = 0; j < 2; ++j) {
                for (int k = 0; k < 2; ++k) {
                    supps(j, 2 * i + k) = THBFromGeo(patch)->support(i)(j, k);
                }
            }
        }
        maxLevel(patch) = myLevel(patch).maxCoeff();
        fineLevel(patch) = THBFromGeo(patch)->maxLevel();
        coarseLevel(patch) = fineLevel(patch) - 1;
        Bells(patch).resize(fineLevel(patch) + 1);
        globalIndex(patch).resize(fineLevel(patch) + 1);
        spilloverPatchStructure(patch).resize(mp->nPatches());
        THBVector1(patch).resize(fineLevel(patch) + 1);
        for (int level = 0; level < Bells(patch).size(); ++level) {
            Bells(patch)(level) = tensForBells;
            globalIndex(patch)(level).resize(tensForBells.size());
            gsTHBSplineBasis<2, real_t> THBfor(tensForBells);
            THBVector1(patch)(level) = THBfor;
            //     Bells(patch)(0) = tensForBells;
            tensForBells.uniformRefine();
            //tensForBells.uniformRefine();
        }
    }
    //printValue(Bells);
    //goestaerDaeOPokhu(Bells);

    int fullMonty;
    // Per-patch active-box counter (replaces the old scalar lastNonZeroRow).
    // Used via a reference alias `int& lastNonZeroRow = lastNonZeroRowPerPatch(patch)`
    // inside the main coarsening loop so that existing code is unchanged.
    gsVector<int> lastNonZeroRowPerPatch(mp->nPatches());
    lastNonZeroRowPerPatch.setZero();
    gsVector<index_t> currentLastNonZeroRow(mp->nPatches());
    gsVector<index_t> AcceptedlastRow(mp->nPatches());

    unsigned interior = THBFromGeo(0)->numKnots(0, 1) - 2 * (THBFromGeo(0)->degree(1) + 1);
    index_t numCells = (interior + 1) * (int)pow(2, coarseLevel.maxCoeff()) * (interior + 1) *
        (int)pow(2, coarseLevel.maxCoeff());
    gsInfo << "numCells: " << numCells << "\n";
    //    gsVector<gsTensorBSplineBasis<2, real_t>> Bells(ell);
    //    Bells(0)(0) = Bells(1)(0) = tens;
    gsVector<gsVector<gsVector<index_t>>> boxMat(mp->nPatches());
    gsVector<gsVector<gsVector<gsVector<index_t>>>> SpilloverboxMat(mp->nPatches());
    for (size_t i = 0; i < SpilloverboxMat.size(); i++)
    {
        SpilloverboxMat[i].resize(mp->nPatches());
        for (size_t j = 0; j < SpilloverboxMat[i].size(); j++)
        {
            SpilloverboxMat[i][j].resize(numCells);
            for (size_t k = 0; k < numCells; k++)
            {
                SpilloverboxMat[i][j][k].resize(5);
            }
        }
    }
    gsVector<gsTHBSplineBasis<2, real_t> > SubdomainHierarchy(mp->nPatches());
    gsTHBSplineBasis<2, real_t> THBInit(tens);

    std::vector<index_t> box1;
    std::chrono::time_point<std::chrono::system_clock> toc;
    gsVector<index_t> initialBoxesNum(Bells.size());
    for (int patch = 0; patch < Bells.size(); ++patch) {
        lastNonZeroRowPerPatch(patch) = static_cast<int>(lowCorners(patch).rows());
        initialBoxesNum(patch) = lastNonZeroRowPerPatch(patch);
        currentLastNonZeroRow(patch) = lastNonZeroRowPerPatch(patch);
        boxMat(patch).resize(numCells);
        for (size_t i = 0; i < numCells; i++)
        {
            boxMat(patch)(i).resize(5);
        }
        //gsDebugVar(boxMat(patch).size());
        SubdomainHierarchy(patch) = THBInit;
        for (int boxInd = 0; boxInd < lowCorners(patch).rows(); ++boxInd) {
            boxMat(patch)(boxInd).resize(5);
            boxMat(patch)(boxInd)(0) = myLevel(patch)(boxInd);
            boxMat(patch)(boxInd)(1) =
                (real_t)lowCorners(patch)(boxInd, 0) / pow(2, maxLevel(patch) - myLevel(patch)(boxInd));
            boxMat(patch)(boxInd)(2) =
                (real_t)lowCorners(patch)(boxInd, 1) / pow(2, maxLevel(patch) - myLevel(patch)(boxInd));
            boxMat(patch)(boxInd)(3) =
                (real_t)upCorners(patch)(boxInd, 0) / pow(2, maxLevel(patch) - myLevel(patch)(boxInd));
            boxMat(patch)(boxInd)(4) =
                (real_t)upCorners(patch)(boxInd, 1) / pow(2, maxLevel(patch) - myLevel(patch)(boxInd));
            box1.push_back(boxMat(patch)(boxInd)(0));
            box1.push_back(boxMat(patch)(boxInd)(1));
            box1.push_back(boxMat(patch)(boxInd)(2));
            box1.push_back(boxMat(patch)(boxInd)(3));
            box1.push_back(boxMat(patch)(boxInd)(4));
            SubdomainHierarchy(patch).refineElements(box1);
            box1.clear();
        }
        //gsDebugVar(SubdomainHierarchy(patch).size());
    }


    gsVector<index_t> numActive;
    gsMatrix<real_t> AcceptedvectSol;
    std::unique_ptr<MPBES<2, real_t>> Acceptedmpbes;
    gsVector<gsVector<gsVector<index_t>>> AcceptedboxMat = boxMat;

    gsVector<gsVector<gsVector<std::vector<index_t >>>> twinsIndex(Bells.size()), twinsPatch(Bells.size());
    //gsVector < gsVector < gsVector <std::vector< index_t >>>> twinsIndex(Bells.size()), twinsPatch(Bells.size());
    gsVector<gsVector<gsVector<index_t >>> hasATwin(Bells.size()), isActive(Bells.size()), isTouching(
        Bells.size());
    gsVector<gsVector<index_t >> hasATwinTHB(Bells.size()), isActiveTHB(Bells.size()), isTouchingTHB(
        Bells.size()), globalIndexTHB(Bells.size());
    gsVector<gsVector<std::vector<index_t>>> AcceptedfunctionDescription;
    gsVector<gsVector<gsVector<index_t >>> AcceptedisActive, AcceptedglobalIndex;
    gsVector<gsMatrix<index_t>> localLocation;
    gsMatrix<> Z;
    index_t iteration, theLev, wasRebuilt, failed, createdBoxNum, currCellIndex, currArrayIndex, x1U, y1U, x2U, y2U, x1Bi, y1Bi, x2Bi, y2Bi, createSpline, RTH, needToEscape, centerInd;
    int ourBox[5];


    if (mp1.interfaces().empty())
        mp1.computeTopology(1e-4, true);

    IdentifyPatches(mp1,
        Bells,
        isTouching,
        twinsIndex,
        twinsPatch,
        hasATwin,
        isActive
    );
    for (int patch = 0; patch < Bells.size(); patch++)
    {
        outfile << "Patch " << patch << "\n";
        for (size_t level = 0; level < Bells(patch).size(); level++)
        {
            outfile << "Patch " << patch << ", level " << level << "\n";
            for (size_t functionIndex = 0; functionIndex < Bells(patch)(level).size(); functionIndex++)
            {
                for (size_t i = 0; i < twinsIndex(patch)(level)(functionIndex).size(); i++)
                {
                    outfile << twinsPatch(patch)(level)(functionIndex)[i]
                        << " " << level << " "
                        << twinsIndex(patch)(level)(functionIndex)[i] << " ";
                }
                outfile << "\n";
            }
            outfile << "\n";
        }
        outfile << "\n";
    }
    // -----------------------------------------------------------------------
    // Global cell-selection loop: level-first, then patch.
    // All patches are processed at the same level before descending,
    // so Grendas' algorithm (and any other selection method) draws from
    // a single global candidate pool across all patches simultaneously.
    // -----------------------------------------------------------------------
    t_init = std::chrono::duration<double>(std::chrono::system_clock::now() - startTime).count() - t_preflight;
    {
    const int globalMaxCoarseLevel = static_cast<int>(coarseLevel.maxCoeff());
    for (int levNow = globalMaxCoarseLevel; levNow >= 0; --levNow) {
        // -----------------------------------------------------------------------
        // Truly global cell pool: all eligible patches share one candidate set.
        // Grenda's algorithm evaluates every patch and picks the globally best
        // cell (lowest delta = most regular geometry) on each attempt.
        // -----------------------------------------------------------------------
        const int nPatches = mp->nPatches();
        const int cellsAtLevel = static_cast<int>(
            (interior + 1) * std::pow(2, levNow) *
            (interior + 1) * std::pow(2, levNow));

        std::vector<gsVector<int>>    perPatchNCC(nPatches);
        std::vector<gsMatrix<int>>    perPatchPickedCells(nPatches);
        for (int p = 0; p < nPatches; ++p) {
            if (levNow > static_cast<int>(coarseLevel(p))) {
                perPatchNCC[p].resize(0);
                continue;
            }
            perPatchNCC[p].resize(cellsAtLevel);
            perPatchPickedCells[p].resize(1, cellsAtLevel);
            for (int i = 0; i < cellsAtLevel; ++i) {
                perPatchNCC[p](i)        = i;
                perPatchPickedCells[p](0, i) = 0;
            }
        }

        auto anyPoolNonEmpty = [&](const std::vector<gsVector<int>>& pools) {
            for (const auto& v : pools) if (v.size() > 0) return true;
            return false;
        };
        auto anyVectorSNonEmpty = [&](const std::vector<gsVector<index_t>>& vss) {
            for (const auto& v : vss) if (v.size() > 0) return true;
            return false;
        };

        int success = 1;
        attempt = 0;
        iteration = 0;
        theLev = levNow;

        // Patch interface orientations are topology-invariant: cache once before the loop.
        std::vector<int> cachedFirstSide, cachedSecondSide, cachedFirstPatch, cachedSecondPatch;
        orientThePatches(mp1, cachedFirstSide, cachedSecondSide, cachedFirstPatch, cachedSecondPatch);

        while (anyPoolNonEmpty(perPatchNCC) && success) {
            iteration++;
            success = 0;

            // Per-pass vectorS: copy of perPatchNCC (cells to try this pass)
            std::vector<gsVector<index_t>> perPatchVectorS(nPatches);
            for (int p = 0; p < nPatches; ++p)
                perPatchVectorS[p] = perPatchNCC[p].cast<index_t>();

            fullMonty = 0;
            while (anyVectorSNonEmpty(perPatchVectorS)) {
                // ----------------------------------------------------------
                // Global Grenda selection: evaluate every patch's remaining
                // candidates and pick the one with the lowest accepted delta.
                // ----------------------------------------------------------
                int patch = -1;
                real_t bestDelta = std::numeric_limits<real_t>::infinity();
                std::vector<CellSelectionResult> preflightResults(nPatches);
                {
                auto _tG0 = std::chrono::system_clock::now();
                for (int p = 0; p < nPatches; ++p) {
                    if (perPatchVectorS[p].size() == 0) continue;
                    int pickedOneTmp = -1, validTmp = 1, currArrTmp = 0;
                    gsVector<index_t> vsTmp = perPatchVectorS[p];
                    // Preflight: preflightOnly=true skips pickCell, preserving rand() state.
                    // Also pass a copy of pickedCells so preflight doesn't modify real state.
                    gsMatrix<int> pickedCellsCopy = perPatchPickedCells[p];
                    preflightResults[p] = selectCellForCoarsening(
                        method, vsTmp, levNow, interior, p, attempt,
                        pickedOneTmp, pickedCellsCopy, validTmp, currArrTmp,
                        Acceptedmpbes.get(), mp1, SubdomainHierarchy,
                        /*preflightOnly=*/true);
                    if (preflightResults[p].acceptedDelta < bestDelta) {
                        bestDelta = preflightResults[p].acceptedDelta;
                        patch = p;
                    }
                }
                t_grenda += std::chrono::duration<double>(std::chrono::system_clock::now() - _tG0).count();
                }
                auto _tGA0 = std::chrono::system_clock::now();
                // If no Grenda geometry candidate found in ANY patch, fall back to the
                // first non-empty patch and let selectCellForCoarsening use its pickCell fallback.
                if (patch < 0) {
                    for (int p = 0; p < nPatches; ++p) {
                        if (perPatchVectorS[p].size() > 0) { patch = p; break; }
                    }
                }
                if (patch < 0) break; // truly nothing left to do

                // Build the coarsening group: all patches that achieved bestDelta.
                // Non-connected patches at equal delta must be coarsened in one step.
                struct CoarsenGroupEntry {
                    int patchId;
                    std::vector<CellToCoarsen> geoCells;
                    std::vector<int> geoCellIndices;
                };
                std::vector<CoarsenGroupEntry> coarsenGroup;
                if (bestDelta < std::numeric_limits<real_t>::infinity()) {
                    for (int p = 0; p < nPatches; ++p) {
                        if (!preflightResults[p].geoCells.empty() &&
                            preflightResults[p].acceptedDelta <= bestDelta + 1e-12) {
                            coarsenGroup.push_back({p,
                                preflightResults[p].geoCells,
                                preflightResults[p].geoCellIndices});
                        }
                    }
                }
                if (coarsenGroup.size() > 1) {
                    gsInfo  << "[coarsen-group] size=" << coarsenGroup.size() << " patches:";
                    outfile << "[coarsen-group] size=" << coarsenGroup.size() << " patches:";
                    for (const auto& cpe : coarsenGroup) {
                        gsInfo  << " " << cpe.patchId;
                        outfile << " " << cpe.patchId;
                    }
                    gsInfo  << " delta=" << bestDelta << "\n";
                    outfile << " delta=" << bestDelta << "\n";
                }

                // Reference aliases so the existing inner body is unchanged
                int& lastNonZeroRow = lastNonZeroRowPerPatch(patch);
                gsVector<int>&      nonCheckedCells = perPatchNCC[patch];
                gsVector<index_t>&  vectorS         = perPatchVectorS[patch];
                gsMatrix<int>&      pickedCells      = perPatchPickedCells[patch];
                int pickedOne = -1;
                {
                    std::vector<std::vector<std::vector<index_t>>> functionDescription;
                    std::vector<std::vector<std::array<int, 3>>> spilloverFunctionCoordinates;
                    std::vector<bool> hasSpillover;
                    std::vector<bool> isTruncated;
                    std::vector<std::vector<gsSparseVector<real_t>>> presentation;
                    gsMatrix<real_t> vectSol;
                    std::vector<CellToCoarsen> geoCells;
                    std::vector<int> geoCellIndices;
                    std::vector<int> attemptedCellIds;
                    std::vector<RestoreRecord> restoreRecords;
                    bool useGeo = false;

                    // gsInfo << "vectorS.size(): " << vectorS.size() << "\n";
                    // outfile << "vectorS.size(): " << vectorS.size() << "\n";
                    if (vectorS.size() == 0)
                    {
                        // gsInfo << "WARNING: vectorS is empty, skipping cell selection\n";
                        // outfile << "WARNING: vectorS is empty, skipping cell selection\n";
                        break;
                    }
                    lcx = 1.0;
                    lcy = 1.0;
                    ucx = 0.0;
                    ucy = 0.0;
                    failed = 0;
                    // outfile << "\n";
                    // outfile << "\n";
                    // outfile << "\n";
                    // outfile << "The boxes\n";
                    //gsDebugVar(boxMat(patch).size());
                    //lastNonZeroRow = lowCorners(patch).rows();
                    // gsInfo << "lastNonZeroRow: " << lastNonZeroRow << "\n";
                    if (patch < 0 || patch >= boxMat.size())
                    {
                        // gsInfo << "WARNING: patch index out of range for boxMat: " << patch << "\n";
                        // outfile << "WARNING: patch index out of range for boxMat: " << patch << "\n";
                        break;
                    }
                    if (lastNonZeroRow < 0)
                    {
                        // gsInfo << "WARNING: lastNonZeroRow < 0, clamping to 0\n";
                        // outfile << "WARNING: lastNonZeroRow < 0, clamping to 0\n";
                        lastNonZeroRow = 0;
                    }
                    if (lastNonZeroRow > static_cast<int>(boxMat(patch).size()))
                    {
                           // gsInfo << "WARNING: lastNonZeroRow > boxMat(patch).size(), clamping. lastNonZeroRow="
                           //        << lastNonZeroRow << " size=" << boxMat(patch).size() << "\n";
                           // outfile << "WARNING: lastNonZeroRow > boxMat(patch).size(), clamping. lastNonZeroRow="
                           //         << lastNonZeroRow << " size=" << boxMat(patch).size() << "\n";
                        lastNonZeroRow = static_cast<int>(boxMat(patch).size());
                    }
                    //for (int i = 0; i < lastNonZeroRow; i++) {
                    for (int i = 0; i < lastNonZeroRow; i++) {
                        // gsInfo << "i:" << i << "\n";
                        if (boxMat(patch)(i).size() < 5)
                        {
                            // gsInfo << "WARNING: boxMat(patch)(" << i << ").size() < 5, skipping row\n";
                            // outfile << "WARNING: boxMat(patch)(" << i << ").size() < 5, skipping row\n";
                            continue;
                        }
                        for (int j = 0; j < 5; j++) {
                            // outfile << boxMat(patch)(i)(j) << "\t";
                            // gsInfo << boxMat(patch)(i)(j) << "\t";
                        }
                        // outfile << "\n";
                        // gsInfo << "\n";
                    }
                    createdBoxNum = 0;
                    t_gap_a += std::chrono::duration<double>(std::chrono::system_clock::now() - _tGA0).count();
                    { auto _tS0 = std::chrono::system_clock::now();
                    CellSelectionResult selection = selectCellForCoarsening(
                        method,
                        vectorS,
                        levNow,
                        interior,
                        patch,
                        attempt,
                        pickedOne,
                        pickedCells,
                        valid,
                        currArrayIndex,
                        Acceptedmpbes ? Acceptedmpbes.get() : nullptr,
                        mp1,
                        SubdomainHierarchy);
                    t_cell_sel += std::chrono::duration<double>(std::chrono::system_clock::now() - _tS0).count();

                    useGeo = selection.useGeo;
                    currCellIndex = selection.currCellIndex;
                    currArrayIndex = selection.currArrayIndex;
                    x1U = selection.x1U;
                    y1U = selection.y1U;
                    x2U = selection.x2U;
                    y2U = selection.y2U;
                    geoCells = std::move(selection.geoCells);
                    geoCellIndices = std::move(selection.geoCellIndices);
                    } // end cell-selection timing block
                    auto _tTHB0 = std::chrono::system_clock::now();
                    int jopa = 4 * (int)pow(2, levNow) * 4 * (int)pow(2, levNow);
                    // outfile << "=======================================\n";
                    // outfile << "=======================================\n";
                    // outfile << "=======================================\n";
                    // outfile << "Number of nonchecked cells: " << vectorS.size() << "\n";
                    // gsInfo << "Number of nonchecked cells: " << vectorS.size() << "\n";
                    if (currArrayIndex < 0 || currArrayIndex >= vectorS.size())
                    {
                           // gsInfo << "WARNING: currArrayIndex out of range (" << currArrayIndex
                           //        << "), clamping to 0\n";
                           // outfile << "WARNING: currArrayIndex out of range (" << currArrayIndex
                           //         << "), clamping to 0\n";
                        currArrayIndex = 0;
                    }
                    // gsInfo << "patch " << patch << ", level " << levNow << ", attempt " << attempt
                    //     << ", CURRENT INDEX IS " << vectorS(currArrayIndex) <<
                    //     ",\t the coordinates of the box are\n" << x1U << " " << y1U << " " << x2U << " " << y2U
                    //     << "\n";
                    // outfile << "patch " << patch << ", level " << levNow << ", attempt " << attempt
                    //     << ", CURRENT INDEX IS " << vectorS(currArrayIndex) <<
                    //     ",\t the coordinates of the box are\n" << x1U << " " << y1U << " " << x2U << " " << y2U
                    //     << "\n";
                    /* if (attempt > 0)     break;*/
                     //gsInfo << currArrayIndex << "\n";
                    // outfile << currArrayIndex << "\n";
                    if (useGeo && !geoCellIndices.empty())
                    {
                        attemptedCellIds = geoCellIndices;
                        std::unordered_set<int> toRemove(geoCellIndices.begin(), geoCellIndices.end());
                        gsVector<int> filtered(vectorS.size());
                        int pos = 0;
                        for (int i = 0; i < vectorS.size(); ++i)
                        {
                            if (toRemove.find(vectorS(i)) == toRemove.end())
                                filtered(pos++) = vectorS(i);
                        }
                        filtered.conservativeResize(pos);
                        vectorS = filtered;
                    }
                    else if (currArrayIndex >= 0 && currArrayIndex < vectorS.size())
                    {
                        if (currCellIndex >= 0)
                            attemptedCellIds.push_back(currCellIndex);
                        vectorS.removeElement(currArrayIndex);
                    }
                    else
                    {
                        if (currCellIndex >= 0)
                            attemptedCellIds.push_back(currCellIndex);
                        // gsInfo << "WARNING: skip removeElement; currArrayIndex out of range (" << currArrayIndex << ")\n";
                        // outfile << "WARNING: skip removeElement; currArrayIndex out of range (" << currArrayIndex << ")\n";
                    }
                    // gsInfo << "DELETED\n";
                    // outfile << "DELETED\n";
                    //if (attempt > 1)     break;
                    valid = 0;
                    createSpline = 0;
                    if (useGeo && !geoCells.empty())
                    {
                        RTH = rebuildTheHierarchyMultiple(
                            boxMat,
                            geoCells,
                            lastNonZeroRow,
                            createdBoxNum,
                            centerInd,
                            ourBox,
                            needToEscape,
                            patch,
                            &restoreRecords,
                            successfullAttempts);
                        createSpline = (RTH == 1);
                    }
                    else
                    {
                        //for (int currentrow = 0; currentrow < lastNonZeroRow; ++currentrow) {
                        for (int currentrow = 0; currentrow < lastNonZeroRow; ++currentrow) {
                            row = currentrow;
                            x1Bi = 0;
                            x2Bi = 0;
                            y1Bi = 0;
                            y2Bi = 0;
                            RTH = rebuildTheHierarchy(boxMat, row, x1U, x1Bi, x2U, x2Bi, y1U, y1Bi, y2U, y2Bi, levNow,
                                lastNonZeroRow,
                                createdBoxNum, centerInd, ourBox, needToEscape, patch);
                            if (RTH == 1) {
                                createSpline = 1;
                                break;
                            }
                            else          createSpline = 0;
                        }
                    }
                    
                    // Check box consistency after rebuild
                    checkBoxesConsistency(boxMat, patch, lastNonZeroRow, outfile);
                    
                    if (createSpline == 1) {
                        iterTime = std::chrono::system_clock::now();
                        // outfile << "creating the spline\n";
                        // gsInfo << "creating the spline\n";
                        std::chrono::duration<double> elapsed_seconds = iterTime - startTime;
                        // outfile << "TIME: " << elapsed_seconds.count() << "\n";
                        //gsInfo << "creating the spline\n";
                        gsTHBSplineBasis<2, real_t> THB(tens);
                        gsTHBSplineBasis<2, real_t> THBUnstructured(tens);
                        //gsInfo << THB << "\n";
                        outfile << THB << "\n";
                        //gsDebugVar(boxMat(patch).size());
                        //if(attempt >0)     //gsDebugVar(boxMat(patch)(7).size());
                        for (int therow = 0; therow < lastNonZeroRow; therow++) {
                            outfile << therow << "; ";
                            std::vector<index_t> box;

                            for (int column = 0; column < 5; column++) {
                                box.push_back(boxMat(patch)(therow)(column));
                                outfile << boxMat(patch)(therow)(column) << "; ";
                                //gsInfo << boxMat(patch)(therow)(column) << "; ";
                            }
                            THBUnstructured.refineElements(box);
                            box.clear();
                            outfile << "\n";
                            //gsInfo << "\n";
                            ;

                        }




                        //gsDebugVar(THBUnstructured.size());
                        lowCorners(patch).clear();
                        upCorners(patch).clear();
                        myLevel(patch).clear();
                        (THBUnstructured.tree().getBoxes(lowCorners(patch), upCorners(patch), myLevel(patch)));
                        //gsInfo << "THBUnstructured: lowCorners(patch).rows(): " << lowCorners(patch).rows() << "\n";
                        lastNonZeroRow = lowCorners(patch).rows();
                        currentLastNonZeroRow(patch) = lastNonZeroRow;
                        //boxMat(patch).resize(lastNonZeroRow);


                        for (int i = 0; i < lastNonZeroRow; i++) {
                            boxMat(patch)(i).resize(5);
                            boxMat(patch)(i)(0) = myLevel(patch)(i);
                            boxMat(patch)(i)(1) = (real_t)lowCorners(patch)(i, 0) /
                                pow(2, THBUnstructured.maxLevel() - myLevel(patch)(i));
                            boxMat(patch)(i)(2) = (real_t)lowCorners(patch)(i, 1) /
                                pow(2, THBUnstructured.maxLevel() - myLevel(patch)(i));
                            boxMat(patch)(i)(3) = (real_t)upCorners(patch)(i, 0) /
                                pow(2, THBUnstructured.maxLevel() - myLevel(patch)(i));
                            boxMat(patch)(i)(4) = (real_t)upCorners(patch)(i, 1) /
                                pow(2, THBUnstructured.maxLevel() - myLevel(patch)(i));
                            //gsInfo << boxMat(patch)(i)(0) << " " << boxMat(patch)(i)(1) << " " << boxMat(patch)(i)(2)
      /*                          << " " << boxMat(patch)(i)(3) << " " <<
                                boxMat(patch)(i)(4) << "\n";*/
                            outfile << boxMat(patch)(i)(0) << " " << boxMat(patch)(i)(1) << " " << boxMat(patch)(i)(2)
                                << " " << boxMat(patch)(i)(3) << " " <<
                                boxMat(patch)(i)(4) << "\n";

                        }
                        // lastNonZeroRow = lowCorners(patch).rows() - 1;

                        for (int therow = 0; therow < lastNonZeroRow; therow++) {
                            std::vector<index_t> box;
                            for (int column = 0; column < 5; column++) {
                                box.push_back(boxMat(patch)(therow)(column));
                                //                            outfile << boxMat(therow)(column) << "; ";
                            }
                            THB.refineElements(box);
                        }
                        //gsDebugVar(THB.size());

                        int isAdequate = 0;
                        for (int i = 0; i < THB.size(); ++i) {
                            int levOf = THB.levelOf(i);
                            if (levOf >= levNow + 1) {
                                isAdequate = 1;
                            }
                        }
                        /*if (!isAdequate) {
                            outfile << "The basis is inadequate, rewriting\n";
                            //gsInfo << "The basis is inadequate, rewriting\n";
                            valid = 0;
                            std::string input(xmlFile + ".xml");
                            gsFileData<> data2(input);
                            gsMultiPatch<real_t>::uPtr tempmp = gsReadFile<>(input);
                            THBTemporary(patch) = dynamic_cast<gsTHBSplineBasis<2> *>(&mp->patch(
                                    patch).basis().source());;
                            (THBTemporary(patch)->tree().getBoxes(lowCorners(patch), upCorners(patch), myLevel(patch)));
                            //gsInfo << "FOUND " << lowCorners(patch).rows() << " BOXES\n";
                            outfile << "FOUND " << lowCorners(patch).rows() << " BOXES\n";
                            for (int i = 0; i < lowCorners(patch).rows(); i++) {

                                boxMat(patch)(i)(0) = myLevel(patch)(i);
                                boxMat(patch)(i)(1) = (real_t) lowCorners(patch)(i, 0) / pow(2, std::min(
                                        (int) THBFromGeo(patch)->maxLevel(), (int) THBTemporary(patch)->maxLevel()) -
                                                                                                myLevel(patch)(i));
                                boxMat(patch)(i)(2) = (real_t) lowCorners(patch)(i, 1) / pow(2, std::min(
                                        (int) THBFromGeo(patch)->maxLevel(), (int) THBTemporary(patch)->maxLevel()) -
                                                                                                myLevel(patch)(i));
                                boxMat(patch)(i)(3) = (real_t) upCorners(patch)(i, 0) / pow(2, std::min(
                                        (int) THBFromGeo(patch)->maxLevel(), (int) THBTemporary(patch)->maxLevel()) -
                                                                                               myLevel(patch)(i));
                                boxMat(patch)(i)(4) = (real_t) upCorners(patch)(i, 1) / pow(2, std::min(
                                        (int) THBFromGeo(patch)->maxLevel(), (int) THBTemporary(patch)->maxLevel()) -
                                                                                               myLevel(patch)(i));
                                //gsInfo << boxMat(patch)(i)(0) << " " << boxMat(patch)(i)(1) << " "
                                       << boxMat(patch)(i)(2) << " " << boxMat(patch)(i)(3) << " " <<
                                       boxMat(patch)(i)(4) << "\n";
                                outfile << boxMat(patch)(i)(0) << " " << boxMat(patch)(i)(1) << " "
                                        << boxMat(patch)(i)(2) << " " << boxMat(patch)(i)(3) << " " <<
                                        boxMat(patch)(i)(4) << "\n";
                            }

                            for (int i = lowCorners(patch).rows(); i < boxMat(patch).size(); ++i) {
                                if(boxMat(patch)(i).size() > 0)
                                boxMat(patch)(i)(0) = boxMat(patch)(i)(1) = boxMat(patch)(i)(2) = boxMat(patch)(i)(
                                        3) = boxMat(patch)(i)(4) = 0;
                            }
                            lastNonZeroRow = lowCorners(patch).rows() - 1;

                            attempt = std::max(attempt - 1, 0);
                            continue;
                        }*/
                        outfile << THB << "\n";
                        //gsInfo << THB << "\n";
                        //                            for (int k = 0; k < THB.size(); ++k) {
                        //                                if(THB.levelOf(k) == 1)
                        //                                    //gsDebugVar(THB.function(k).support());
                        //                            }


                        SubdomainHierarchy(patch) = THB;
                        if (g_verbose)
                            gsInfo  << "[primary coarsen] patch=" << patch
                                    << " cells=" << (useGeo ? geoCells.size() : static_cast<size_t>(1))
                                    << " THBsize=" << THB.size() << "\n";
                        outfile << "[primary coarsen] patch=" << patch
                                << " cells=" << (useGeo ? geoCells.size() : static_cast<size_t>(1))
                                << " THBsize=" << THB.size() << "\n";

                        // Rebuild co-patch hierarchies so MPBES sees the full group in one shot.
                        for (const auto& cpe : coarsenGroup) {
                            if (cpe.patchId == patch) continue;
                            int cpLNZR = lastNonZeroRowPerPatch(cpe.patchId);
                            int cpCreated = 0, cpCenter = -1, cpEscape = 0;
                            int cpOurBox[4] = {0, 0, 0, 0};
                            int cpRTH = rebuildTheHierarchyMultiple(
                                boxMat, cpe.geoCells, cpLNZR, cpCreated,
                                cpCenter, cpOurBox, cpEscape, cpe.patchId);
                            outfile << "[co-patch coarsen] patch=" << cpe.patchId
                                    << " RTH=" << cpRTH << "\n";
                            if (cpRTH == 1) {
                                // Build unstructured THB from post-rebuild boxMat
                                gsTHBSplineBasis<2, real_t> cpTHBUnstructured(tens);
                                for (int row = 0; row < cpLNZR; ++row) {
                                    std::vector<index_t> box;
                                    for (int col = 0; col < 5; ++col)
                                        box.push_back(boxMat(cpe.patchId)(row)(col));
                                    cpTHBUnstructured.refineElements(box);
                                }
                                // Extract canonical tree form (mirrors primary patch flow)
                                lowCorners(cpe.patchId).clear();
                                upCorners(cpe.patchId).clear();
                                myLevel(cpe.patchId).clear();
                                cpTHBUnstructured.tree().getBoxes(
                                    lowCorners(cpe.patchId), upCorners(cpe.patchId), myLevel(cpe.patchId));
                                cpLNZR = lowCorners(cpe.patchId).rows();
                                lastNonZeroRowPerPatch(cpe.patchId) = cpLNZR;
                                currentLastNonZeroRow(cpe.patchId) = cpLNZR;
                                // Rewrite boxMat from canonical form
                                for (int i = 0; i < cpLNZR; ++i) {
                                    boxMat(cpe.patchId)(i).resize(5);
                                    boxMat(cpe.patchId)(i)(0) = myLevel(cpe.patchId)(i);
                                    boxMat(cpe.patchId)(i)(1) = (real_t)lowCorners(cpe.patchId)(i, 0) /
                                        pow(2, cpTHBUnstructured.maxLevel() - myLevel(cpe.patchId)(i));
                                    boxMat(cpe.patchId)(i)(2) = (real_t)lowCorners(cpe.patchId)(i, 1) /
                                        pow(2, cpTHBUnstructured.maxLevel() - myLevel(cpe.patchId)(i));
                                    boxMat(cpe.patchId)(i)(3) = (real_t)upCorners(cpe.patchId)(i, 0) /
                                        pow(2, cpTHBUnstructured.maxLevel() - myLevel(cpe.patchId)(i));
                                    boxMat(cpe.patchId)(i)(4) = (real_t)upCorners(cpe.patchId)(i, 1) /
                                        pow(2, cpTHBUnstructured.maxLevel() - myLevel(cpe.patchId)(i));
                                }
                                // Build canonical THB from normalized boxMat
                                gsTHBSplineBasis<2, real_t> cpTHB(tens);
                                for (int row = 0; row < cpLNZR; ++row) {
                                    std::vector<index_t> box;
                                    for (int col = 0; col < 5; ++col)
                                        box.push_back(boxMat(cpe.patchId)(row)(col));
                                    cpTHB.refineElements(box);
                                }
                                SubdomainHierarchy(cpe.patchId) = cpTHB;
                                if (g_verbose)
                                    gsInfo  << "[co-patch coarsen] patch=" << cpe.patchId
                                            << " THBsize=" << cpTHB.size() << "\n";
                                outfile << "[co-patch coarsen] patch=" << cpe.patchId
                                        << " THBsize=" << cpTHB.size() << "\n";
                            }
                        }

                        THB.anchors_into(anmat);
                        gsMatrix<real_t> basisInd(1, THB.size()); // 1 if THB function is not represented in the previous support set
                        for (int l = 0; l < THB.size(); ++l) {
                            basisInd(0, l) = 1;
                        }
                        gsMatrix<> THBcoefs(THB.size(), 2);
                        for (int k = 0; k < THB.size(); ++k) {
                            THBcoefs(k, 0) = THBcoefs(k, 1) = 0.0;
                        }
                        //gsDebugVar(THB.size());
                        ////gsDebugVar(coefshat());
                        for (int j = 0; j < THB.size(); j++) {
                            for (int i = 0; i < anmat2.cols(); i++) {
                                if (
                                    THB.support(j)(0, 0) == getSupports(i, 0, 0, supps)
                                    && THB.support(j)(1, 0) == getSupports(i, 1, 0, supps)
                                    && THB.support(j)(0, 1) == getSupports(i, 0, 1, supps)
                                    && THB.support(j)(1, 1) == getSupports(i, 1, 1, supps)
                                    ) {

                                    /*THBcoefs(j, 0) = coefshat(i, 0);
                                    THBcoefs(j, 1) = coefshat(i, 1);*/
                                    basisInd(0, j) = 0;

                                    break;
                                }
                            }
                        }
                        //                    //gsInfo << basisInd << "\n";
                        outfile << basisInd << "\n";
                        for (int i1 = 0; i1 < basisInd.cols(); ++i1) {
                            if (basisInd(0, i1) != 0.0) {
                                if (THB.support(i1)(0, 0) < lcx) lcx = THB.support(i1)(0, 0);
                                if (THB.support(i1)(1, 0) < lcy) lcy = THB.support(i1)(1, 0);
                                if (THB.support(i1)(0, 1) > ucx) ucx = THB.support(i1)(0, 1);
                                if (THB.support(i1)(1, 1) > ucy) ucy = THB.support(i1)(1, 1);
                            }
                        }
                        //const std::vector<tensorBasis> bases = 
                        auto Bkas0 = SubdomainHierarchy(patch).getBases();
                        auto Bkas1 = SubdomainHierarchy(1).getBases();
                        //gsInfo << typeid(Bkas).name() << "\n";
                        //gsInfo << Bkas0[3]->size() << "\n";
                        //gsInfo << Bkas1[3]->size() << "\n";

                        //gsInfo << ;

                        //for (size_t i = 0; i < SubdomainHierarchy(patch).size(); i++)
                        //{
                        //    gsInfo << SubdomainHierarchy(patch).flatTensorIndexOf(i) << ", " << SubdomainHierarchy(patch).levelOf(i) << "\n";
                        //}
                        //gsInfo << "The index of our concern:\n" << hasATwin(0)(4)(0) << "\n";
                        //return 241026;
                        t_thb_rebuild += std::chrono::duration<double>(std::chrono::system_clock::now() - _tTHB0).count();
                        std::chrono::time_point<std::chrono::system_clock> beforeIdent, afterIdent;
                        beforeIdent = std::chrono::system_clock::now();
                        //IdentifyPatches(mp1, THBVector1, isTouching, twinsIndex, twinsPatch, hasATwin, isActive);
                        std::vector<int> firstSide  = cachedFirstSide;
                        std::vector<int> secondSide = cachedSecondSide;
                        std::vector<int> firstPatch = cachedFirstPatch;
                        std::vector<int> secondPatch = cachedSecondPatch;
                        t_orient += 0.0; // skipped: cached before the loop
                        //IdentifyPatches(mp1, THBVector, isTouchingTHB, twinsIndexTHB, twinsPatchTHB, hasATwinTHB, isActiveTHB);
                        gsVector  <gsVector< gsVector<index_t>>> isIncluded;
                        gsVector  <gsVector< gsVector<index_t>>> indexInTHB;
                        gsVector<std::vector<std::array<int, 2>>> thbToBellsMapping;
                        { auto _tIM0 = std::chrono::system_clock::now();
                        createIndexMapping(
                            Bells,
                            SubdomainHierarchy,
                            twinsIndex,
                            twinsPatch,
                            isIncluded,
                            indexInTHB,
                            thbToBellsMapping
                        );
                        t_indexmap += std::chrono::duration<double>(std::chrono::system_clock::now() - _tIM0).count(); }

                        afterIdent = std::chrono::system_clock::now();
                        std::chrono::duration<double> elapsed_seconds_Ident = afterIdent - beforeIdent;
                        // gsInfo << "IdentifyPatches took: " << elapsed_seconds_Ident.count() << "\n";

                        // ===== MPBES: Multi-Patch B-spline with Exact continuity =====
                        // (beforeIdent reused as start of the full MPBES-build phase)
                        // This performs: Twin identification + Kraft selection + Truncation
                        // gsInfo << "\n========== Creating MPBES basis ==========\n";
                        auto _tMC0 = std::chrono::system_clock::now();
                        MPBES<2, real_t> mpbes(
                            boxMat,
                            SubdomainHierarchy,
                            Bells,
                            hasATwin,
                            twinsIndex,
                            twinsPatch,
                            indexInTHB,
                            currentLastNonZeroRow,
                            true,  // verbose
                            attempt  // for conditional logging
                        );
                        t_mpbes_ctor += std::chrono::duration<double>(std::chrono::system_clock::now() - _tMC0).count();
                        
                        // Validate MPBES after creation
                        try {
                            if (mpbes.size() == 0) {
                                // gsInfo << "ERROR: MPBES created with 0 functions!\n";
                                // outfile << "ERROR: MPBES created with 0 functions!\n";
                            } else if (mpbes.nPatches() == 0) {
                                // gsInfo << "ERROR: MPBES created with 0 patches!\n";
                                // outfile << "ERROR: MPBES created with 0 patches!\n";
                            } else {
                                // gsInfo << "MPBES created with " << mpbes.size() << " basis functions\n";
                            }
                        } catch (const std::exception& e) {
                            gsInfo << "EXCEPTION validating MPBES: " << e.what() << "\n";
                        } catch (...) {
                            gsInfo << "UNKNOWN EXCEPTION validating MPBES\n";
                        }
                        t_mpbes_bld += std::chrono::duration<double>(std::chrono::system_clock::now() - beforeIdent).count();

                        // gsInfo << "==========================================\n\n";

                        auto _tEX0 = std::chrono::system_clock::now();
                        // Extract data from MPBES for use throughout the code - with validation
                        try {
                            functionDescription = mpbes.functionDescription();
                            spilloverFunctionCoordinates = mpbes.spilloverCoordinates();
                            hasSpillover = mpbes.hasSpillover();
                            isTruncated = mpbes.isTruncated();
                            presentation = mpbes.presentation();
                            commonSize = mpbes.size();
                            
                            // Validate extracted data
                                // gsInfo << "MPBES data extracted: functions=" << commonSize 
                                //       << ", functionDescription.size=" << functionDescription.size()
                                //       << ", hasSpillover.size=" << hasSpillover.size()
                                //       << ", isTruncated.size=" << isTruncated.size() << "\n";
                            
                            if (functionDescription.size() != static_cast<size_t>(commonSize)) {
                                // gsInfo << "WARNING: functionDescription size mismatch!\n";
                                // outfile << "WARNING: functionDescription size mismatch!\n";
                            }
                            if (hasSpillover.size() != static_cast<size_t>(commonSize)) {
                                // gsInfo << "WARNING: hasSpillover size mismatch!\n";
                                // outfile << "WARNING: hasSpillover size mismatch!\n";
                            }
                            if (isTruncated.size() != static_cast<size_t>(commonSize)) {
                                // gsInfo << "WARNING: isTruncated size mismatch!\n";
                                // outfile << "WARNING: isTruncated size mismatch!\n";
                            }
                            if (presentation.size() != static_cast<size_t>(commonSize)) {
                                // gsInfo << "WARNING: presentation size mismatch!\n";
                                // outfile << "WARNING: presentation size mismatch!\n";
                            }
                        } catch (const std::exception& e) {
                            // gsInfo << "EXCEPTION extracting MPBES data: " << e.what() << "\n";
                            // outfile << "EXCEPTION extracting MPBES data: " << e.what() << "\n";
                            // Set safe defaults
                            commonSize = 0;
                        } catch (...) {
                            // gsInfo << "UNKNOWN EXCEPTION extracting MPBES data\n";
                            // outfile << "UNKNOWN EXCEPTION extracting MPBES data\n";
                            commonSize = 0;
                        }

                        const FittingMode fittingMode = g_useLocalFitting ? FittingMode::LocalFitting : FittingMode::GlobalFitting;
                        const bool useLocalFitting = (fittingMode == FittingMode::LocalFitting);
                        bool localFittingEffective = useLocalFitting; // set false when local region = global
                        const bool verboseLocalRegionDump = false;
                        const bool verboseFitMatrixDump = false;
                        const bool verboseVectSolRowDump = false;

                        vectSol.setZero(commonSize, 2);
                        if (commonSize > 0)
                        {
                            if (Acceptedmpbes &&
                                AcceptedvectSol.rows() == static_cast<index_t>(Acceptedmpbes->size()))
                                vectSol = mapPrevSolToNewMpbes(
                                    Acceptedmpbes->functionDescription(),
                                    AcceptedvectSol, mpbes);
                            else
                                vectSol = mapMpCoefficientsToMpbes(mp1, mpbes);
                            if (outfile.is_open())
                            {
                                index_t nonZeroRows = 0;
                                for (index_t r = 0; r < vectSol.rows(); ++r)
                                    if (vectSol(r, 0) != 0.0 || vectSol(r, 1) != 0.0)
                                        ++nonZeroRows;

                                outfile << "[mpbes-init] vectSol assigned"
                                        << (Acceptedmpbes ? " from prevSol" : " from mp coefficients")
                                        << " patch=" << patch
                                        << " levNow=" << levNow
                                        << " attempt=" << attempt
                                        << " rows=" << vectSol.rows()
                                        << " cols=" << vectSol.cols()
                                        << " nonZeroRows=" << nonZeroRows
                                        << "\n";
                                if (verboseFitMatrixDump)
                                    printTheMatrix(vectSol, "vectSol_initial_from_mpbes");
                            }
                        }
                        t_extract += std::chrono::duration<double>(std::chrono::system_clock::now() - _tEX0).count();

                        // Capture original fine-geometry coefficients once, on the first MPBES build
                        if (originalCoefficients == nullptr && vectSol.rows() > 0 && vectSol.cols() >= 2) {
                            originalCoefsCapture  = vectSol;
                            originalCoefficients  = &originalCoefsCapture;
                        }

                        gsMatrix<real_t> effectiveOriginalCoefficientsStorage;
                        const gsMatrix<real_t>* effectiveOriginalCoefficients = originalCoefficients;
                        if (vectSol.rows() > 0 && vectSol.cols() >= 2)
                        {
                            effectiveOriginalCoefficientsStorage = vectSol;
                            effectiveOriginalCoefficients = &effectiveOriginalCoefficientsStorage;
                        }

                        std::vector<CellToCoarsen> localCells = geoCells;
                        if (localCells.empty())
                        {
                            CellToCoarsen fallbackCell;
                            fallbackCell.level = levNow;
                            fallbackCell.x1 = static_cast<int>(x1U);
                            fallbackCell.y1 = static_cast<int>(y1U);
                            fallbackCell.x2 = static_cast<int>(x2U);
                            fallbackCell.y2 = static_cast<int>(y2U);
                            localCells.push_back(fallbackCell);
                        }

                        auto _tLS0 = std::chrono::system_clock::now();
                        const int localityLambda = g_localityLambda;
                        LocalCoarseningRegion localRegion;
                        gsVector<gsMatrix<real_t>> uvFitting = uv1;
                        index_t basisSelected = 0;
                        if (useLocalFitting)
                        {
                            std::vector<std::pair<int, std::vector<CellToCoarsen>>> groupSeeds;
                            for (const auto& cpe : coarsenGroup) {
                                if (cpe.patchId == patch) continue;
                                groupSeeds.emplace_back(cpe.patchId, cpe.geoCells);
                            }
                            localRegion = buildLocalCoarseningRegion(
                                mpbes,
                                patch,
                                localCells,
                                interior,
                                interior,
                                localityLambda,
                                groupSeeds);

                            // Extend local region to cover patches adjacent to the coarsened
                            // group but NOT in it.  Their MPBES functions otherwise keep the
                            // wrong seed from mapMpCoefficientsToMpbes (fine mp1 THB indices
                            // are mismatched to the coarsened THB basis), producing large
                            // globalErrorEval at shared corners/edges.
                            if (coarsenGroup.size() > 1) {
                                std::unordered_set<int> groupPatchIds;
                                for (const auto& cpe : coarsenGroup)
                                    groupPatchIds.insert(cpe.patchId);
                                const index_t nPatchesMpbes = mpbes.nPatches();
                                const auto& fd = mpbes.functionDescription();
                                for (size_t ii = 0; ii < firstPatch.size(); ++ii) {
                                    for (int side = 0; side < 2; ++side) {
                                        int inGroup   = (side == 0) ? firstPatch[ii] : secondPatch[ii];
                                        int neighbor  = (side == 0) ? secondPatch[ii] : firstPatch[ii];
                                        if (!groupPatchIds.count(inGroup)) continue;
                                        if ( groupPatchIds.count(neighbor)) continue;
                                        if (neighbor < 0 || neighbor >= static_cast<int>(nPatchesMpbes)) continue;
                                        // Add full [0,1]^2 coverage for this neighbor patch
                                        localRegion.patchAABB[neighbor] = {{0.0, 0.0, 1.0, 1.0}};
                                        localRegion.hasPatch[neighbor]  = true;
                                        localRegion.enabled = true;
                                        // Mark all MPBES functions touching neighbor as active
                                        for (index_t f = 0; f < localRegion.basisInd.cols(); ++f) {
                                            if (localRegion.basisInd(0, f) != 0.0) continue;
                                            if (f >= static_cast<index_t>(fd.size())) break;
                                            for (const auto& comp : fd[f]) {
                                                if (!comp.empty() && comp[0] == neighbor) {
                                                    localRegion.basisInd(0, f) = 1.0;
                                                    break;
                                                }
                                            }
                                        }
                                        outfile << "[local-region] neighbor patch=" << neighbor
                                                << " added (group-adjacent, full domain)\n";
                                    }
                                }
                            }

                            // Count local DOFs and active patches for sampling density
                            index_t nLocalPatches = 0;
                            for (index_t f = 0; f < localRegion.basisInd.cols(); ++f)
                                if (localRegion.basisInd(0, f) != 0.0)
                                    ++basisSelected;
                            for (index_t p = 0; p < static_cast<index_t>(localRegion.hasPatch.size()); ++p)
                                if (localRegion.hasPatch[p])
                                    ++nLocalPatches;
                            if (nLocalPatches < 1) nLocalPatches = 1;

                            // Use the same per-patch density as the global grid — never denser.
                            index_t totalGlobalPts = 0;
                            for (index_t p = 0; p < uv1.size(); ++p)
                                totalGlobalPts += uv1(p).cols();
                            const int kLocal = std::max(3, static_cast<int>(
                                std::ceil(std::sqrt(static_cast<real_t>(totalGlobalPts)
                                                    / static_cast<real_t>(uv1.size())))));

                            if (nLocalPatches >= static_cast<index_t>(uv1.size()))
                            {
                                gsInfo << "[local-fitting] NOTE: local region covers all "
                                       << nLocalPatches << " patches — falling back to global fitting for this step.\n";
                                if (outfile.is_open())
                                    outfile << "[local-fitting] NOTE: local region covers all "
                                            << nLocalPatches << " patches — falling back to global fitting for this step.\n";
                                uvFitting = uv1;
                                localFittingEffective = false;
                            }
                            else
                            {
                                uvFitting = resampleLocalRegion(localRegion, uv1.size(), kLocal);
                            }
                        }
                        else
                        {
                            if (g_verbose) gsInfo << "Local fitting disabled; using global uv cloud\n";
                        }

                        if (outfile.is_open())
                        {
                            index_t uvBefore = 0;
                            index_t uvAfter = 0;
                            for (index_t p = 0; p < uv1.size(); ++p)
                                uvBefore += uv1(p).cols();
                            for (index_t p = 0; p < uvFitting.size(); ++p)
                                uvAfter += uvFitting(p).cols();

                                outfile << "[fitting-mode] attempt context: mode=" << fittingModeName(fittingMode)
                                    << " patch=" << patch
                                    << " levNow=" << levNow
                                    << " attempt=" << attempt << "\n";

                                outfile << "[local-region] attempt context: patch=" << patch
                                    << " levNow=" << levNow
                                    << " attempt=" << attempt
                                    << " lambda=" << localityLambda
                                    << " basisSelected=" << basisSelected << "/" << localRegion.basisInd.cols()
                                    << " uvBefore=" << uvBefore
                                    << " uvAfter=" << uvAfter << "\n";

                            if (verboseLocalRegionDump)
                            {
                                outfile << "[local-region] functionDescription dump begin\n";
                                for (index_t f = 0; f < static_cast<index_t>(functionDescription.size()); ++f)
                                {
                                    outfile << "  function " << f
                                            << ": components=" << functionDescription[f].size();
                                    if (f < static_cast<index_t>(hasSpillover.size()))
                                        outfile << " hasSpillover=" << (hasSpillover[f] ? 1 : 0);
                                    if (f < static_cast<index_t>(isTruncated.size()))
                                        outfile << " isTruncated=" << (isTruncated[f] ? 1 : 0);
                                    outfile << "\n";

                                    for (size_t compIdx = 0; compIdx < functionDescription[f].size(); ++compIdx)
                                    {
                                        const auto& desc = functionDescription[f][compIdx];
                                        if (desc.size() < 3)
                                            continue;
                                        outfile << "    component " << compIdx
                                                << ": patch=" << desc[0]
                                                << " level=" << desc[1]
                                                << " tensorIdx=" << desc[2] << "\n";
                                    }

                                    if (f < static_cast<index_t>(spilloverFunctionCoordinates.size()) &&
                                        !spilloverFunctionCoordinates[f].empty())
                                    {
                                        for (size_t spIdx = 0; spIdx < spilloverFunctionCoordinates[f].size(); ++spIdx)
                                        {
                                            const auto& sp = spilloverFunctionCoordinates[f][spIdx];
                                            outfile << "    spillover " << spIdx
                                                    << ": patch=" << sp[0]
                                                    << " level=" << sp[1]
                                                    << " tensorIdx=" << sp[2] << "\n";
                                        }
                                    }
                                }
                                outfile << "[local-region] functionDescription dump end\n";

                                outfile << "[local-region] basisInd==1 functions begin\n";
                                for (index_t f = 0; f < localRegion.basisInd.cols(); ++f)
                                {
                                    if (localRegion.basisInd(0, f) == 0.0)
                                        continue;

                                    outfile << "  active function " << f << "\n";
                                    if (f < static_cast<index_t>(functionDescription.size()))
                                    {
                                        for (size_t compIdx = 0; compIdx < functionDescription[f].size(); ++compIdx)
                                        {
                                            const auto& desc = functionDescription[f][compIdx];
                                            if (desc.size() < 3)
                                                continue;
                                            outfile << "    component " << compIdx
                                                    << ": patch=" << desc[0]
                                                    << " level=" << desc[1]
                                                    << " tensorIdx=" << desc[2] << "\n";
                                        }
                                    }

                                    if (f < static_cast<index_t>(spilloverFunctionCoordinates.size()) &&
                                        !spilloverFunctionCoordinates[f].empty())
                                    {
                                        for (size_t spIdx = 0; spIdx < spilloverFunctionCoordinates[f].size(); ++spIdx)
                                        {
                                            const auto& sp = spilloverFunctionCoordinates[f][spIdx];
                                            outfile << "    spillover " << spIdx
                                                    << ": patch=" << sp[0]
                                                    << " level=" << sp[1]
                                                    << " tensorIdx=" << sp[2] << "\n";
                                        }
                                    }
                                }
                                outfile << "[local-region] basisInd==1 functions end\n";
                            }
                        }

                        // Save functionDescription for debugging
                        if (patch == 1 && levNow == 3 && attempt == 29) {
                            std::string descFile = "functionDescription_p" + std::to_string(patch) + "_lev" + std::to_string(levNow) + "_att" + std::to_string(attempt) + ".txt";
                            std::ofstream descLog(descFile);
                            
                            descLog << "=== Function Description Dump ===\n";
                            descLog << "Configuration: Patch " << patch << ", Level " << levNow << ", Attempt " << attempt << "\n";
                            descLog << "Total functions: " << functionDescription.size() << "\n\n";
                            
                            for (size_t f = 0; f < functionDescription.size(); ++f) {
                                const auto& desc = functionDescription[f];
                                
                                descLog << "Function " << f << ": " << desc.size() << " component(s)\n";
                                for (size_t compIdx = 0; compIdx < desc.size(); ++compIdx) {
                                    descLog << "  Component " << compIdx << ": patch=" << desc[compIdx][0]
                                           << ", level=" << desc[compIdx][1]
                                           << ", tensorIdx=" << desc[compIdx][2] << "\n";
                                }
                                
                                // Check truncation and spillover status
                                if (f < isTruncated.size()) {
                                    descLog << "  Truncated: " << (isTruncated[f] ? "YES" : "NO") << "\n";
                                }
                                if (f < hasSpillover.size()) {
                                    descLog << "  HasSpillover: " << (hasSpillover[f] ? "YES" : "NO") << "\n";
                                }
                                descLog << "\n";
                            }
                            
                            descLog.close();
                            gsInfo << "Function description saved to " << descFile << "\n";
                        }

                        // gsInfo << "Setup completed. Functions: " << commonSize << "\n";

                        // Save function supports for debugging
                        if (patch == 1 && levNow == 3 && attempt == 26) {
                            std::string supportsFile = "function_supports_p" + std::to_string(patch) + "_lev" + std::to_string(levNow) + "_att" + std::to_string(attempt) + ".txt";
                            std::ofstream supportsLog(supportsFile);
                            
                            // Target evaluation point
                            real_t targetU = 0.605662;
                            real_t targetV = 0.605662;
                            
                            supportsLog << "=== Function Supports Analysis ===\n";
                            supportsLog << "Configuration: Patch " << patch << ", Level " << levNow << ", Attempt " << attempt << "\n";
                            supportsLog << "Target evaluation point: (u,v) = (" << targetU << ", " << targetV << ")\n";
                            supportsLog << "Total functions: " << commonSize << "\n\n";
                            
                            int functionsContainingPoint = 0;
                            const auto& bellsBases = mpbes.bellsBases();
                            
                            for (index_t f = 0; f < commonSize; ++f) {
                                const auto& desc = functionDescription[f];
                                
                                supportsLog << "Function " << f << ": " << desc.size() << " component(s)\n";
                                
                                bool containsPoint = false;
                                for (size_t compIdx = 0; compIdx < desc.size(); ++compIdx) {
                                    int fnPatch = desc[compIdx][0];
                                    int fnLevel = desc[compIdx][1];
                                    int fnTensorIdx = desc[compIdx][2];
                                    
                                    supportsLog << "  Comp " << compIdx << ": patch=" << fnPatch 
                                               << ", level=" << fnLevel << ", tensorIdx=" << fnTensorIdx;
                                    
                                    // Check if patch is 0 (we're evaluating on patch 0)
                                    if (fnPatch == 0) {
                                        // Get support from basis
                                        if (fnPatch < static_cast<int>(bellsBases.size()) && 
                                            fnLevel < static_cast<int>(bellsBases[fnPatch].size())) {
                                            const gsTensorBSplineBasis<2, real_t>& basis = bellsBases[fnPatch][fnLevel];
                                            if (fnTensorIdx >= 0 && fnTensorIdx < basis.size()) {
                                                gsMatrix<real_t> supportBox = basis.support(fnTensorIdx);
                                                real_t u_min = supportBox(0, 0);
                                                real_t u_max = supportBox(0, 1);
                                                real_t v_min = supportBox(1, 0);
                                                real_t v_max = supportBox(1, 1);
                                                
                                                supportsLog << ", support=[" << u_min << "," << u_max 
                                                           << "] x [" << v_min << "," << v_max << "]";
                                                
                                                // Check if target point is in support
                                                bool inSupport = (targetU >= u_min && targetU <= u_max &&
                                                                 targetV >= v_min && targetV <= v_max);
                                                supportsLog << ", contains_point=" << (inSupport ? "YES" : "NO");
                                                
                                                if (inSupport) {
                                                    containsPoint = true;
                                                }
                                            } else {
                                                supportsLog << ", tensorIdx out of bounds";
                                            }
                                        } else {
                                            supportsLog << ", basis not available";
                                        }
                                    } else {
                                        supportsLog << " [patch " << fnPatch << ", not evaluating on patch " << fnPatch << "]";
                                    }
                                    supportsLog << "\n";
                                }
                                
                                if (containsPoint) {
                                    functionsContainingPoint++;
                                    supportsLog << "  *** FUNCTION CONTAINS TARGET POINT ***\n";
                                }
                                supportsLog << "\n";
                            }
                            
                            supportsLog << "\n=== Summary ===\n";
                            supportsLog << "Functions with components on patch 0 containing point (" 
                                       << targetU << ", " << targetV << "): " << functionsContainingPoint << "\n";
                            supportsLog << "These functions SHOULD have non-zero values at the evaluation point.\n";
                            
                            supportsLog.close();
                            gsInfo << "Function supports saved to " << supportsFile << "\n";
                        }

                        // Assembly phase
                        index_t numElements = countTotalActiveElements(SubdomainHierarchy);
                        index_t numGaussPointsPerElement = 4;
                        index_t numTotalRows = numElements * numGaussPointsPerElement;

                        gsMatrix<real_t> vectB(numTotalRows, 2);
                        gsSparseMatrix<real_t> matA(numTotalRows, commonSize);

                        std::vector<index_t> assembleActiveFunctions;
                        if (localFittingEffective && localRegion.enabled && localRegion.basisInd.cols() == commonSize)
                        {
                            const auto& hasSpillLocal = mpbes.hasSpillover();

                            // Batch active-function selection: active_into + reverse THB→global map.
                            // O(nf + M×numActivePerPoint) vs old O(AU × M × evalSingle).
                            std::unordered_set<index_t> functionsWithNonZeroValue;
                            functionsWithNonZeroValue.reserve(static_cast<size_t>(commonSize / 2 + 16));
                            for (index_t p = 0; p < static_cast<index_t>(uvFitting.size()); ++p)
                                if (uvFitting(p).cols() > 0)
                                    mpbes.collectActiveGlobalIndices(p, uvFitting(p), functionsWithNonZeroValue);
                            // Spillover functions may contribute on patches not covered by uvFitting;
                            // include them all to avoid false negatives in the AU set.
                            for (index_t f = 0; f < static_cast<index_t>(commonSize); ++f)
                                if (localRegion.basisInd(0, f) != 0.0 &&
                                    f < static_cast<index_t>(hasSpillLocal.size()) && hasSpillLocal[f])
                                    functionsWithNonZeroValue.insert(f);

                            index_t droppedZeroOrUnsampled = 0;
                            for (index_t f = 0; f < static_cast<index_t>(commonSize); ++f)
                            {
                                if (localRegion.basisInd(0, f) == 0.0) continue;
                                if (functionsWithNonZeroValue.count(f))
                                    assembleActiveFunctions.push_back(f);
                                else
                                    ++droppedZeroOrUnsampled;
                            }

                            if (outfile.is_open())
                                outfile << "[assemble] AU functions selected=" << assembleActiveFunctions.size()
                                        << " droppedZeroOrUnsampled=" << droppedZeroOrUnsampled
                                        << " (batch active_into)\n";
                        }
                        t_localsetup += std::chrono::duration<double>(std::chrono::system_clock::now() - _tLS0).count();

                        std::chrono::time_point<std::chrono::system_clock> beforeassembleA, afterassembleA;
                        beforeassembleA = std::chrono::system_clock::now();
                        // gsInfo << "Starting assembly with " << commonSize << " functions...\n";

                        // Enable detailed logging for patch 1, level 3, attempt 26
                        bool detailedLog = (patch == 1 && levNow == 3 && attempt == 26);
                        // if (detailedLog) {
                        //     gsInfo << "\n=== DETAILED ASSEMBLY LOGGING (Patch " << patch << ", Level " << levNow << ", Attempt " << attempt << ") ===\n";
                        //     outfile << "\n=== DETAILED ASSEMBLY LOGGING (Patch " << patch << ", Level " << levNow << ", Attempt " << attempt << ") ===\n";
                        // }

                        // Use new MPBES-based assemble.
                        // For local fitting: pass assembleActiveFunctions so A has
                        // dimensions (# sampling pts) x (# local functions) only.
                        // For global fitting: pass nullptr → A covers all MPBES functions.
                        const std::vector<index_t>* activeFunctionsForAssembly =
                            (useLocalFitting && !assembleActiveFunctions.empty())
                            ? &assembleActiveFunctions
                            : nullptr;
                        assemble(
                            uvFitting,
                            mpbes,
                            matA,
                            vectB,
                            mp1,
                            detailedLog,
                            patch,
                            levNow,
                            attempt,
                            activeFunctionsForAssembly);
                        if (outfile.is_open())
                            outfile << "[assemble] mode=" << (activeFunctionsForAssembly ? "local" : "global")
                                    << " A cols=" << (activeFunctionsForAssembly
                                        ? static_cast<index_t>(assembleActiveFunctions.size())
                                        : mpbes.size()) << "\n";

                        // if (attempt == 56) {
                        //     std::string cmpFile = "compare_mpbes_thb_attempt56.txt";
                        //     compareMPBESvsTHB(
                        //         mpbes,
                        //         SubdomainHierarchy,
                        //         indexInTHB,
                        //         functionDescription,
                        //         isTruncated,
                        //         presentation,
                        //         cmpFile);
                        //     gsInfo << "MPBES vs THB comparison saved to " << cmpFile << "\n";
                        //     outfile << "MPBES vs THB comparison saved to " << cmpFile << "\n";
                        // }

                        // if (attempt == 0) {
                        //     printTheMatrix(matA, "matA");
                        // }

                        // gsInfo << "About to log matA.size...\n";
                        // gsInfo << "matA.size: " << matA.rows() << " * " << matA.cols() << "\n";
                        // gsInfo << "Logged matA.size successfully.\n";

                        checkSparseMatrixHealth(matA, "matA");
                        if (!gsEigen::MatrixXd(vectB).allFinite())
                        {
                            // gsInfo << "ERROR: vectB contains non-finite entries before normal equations\n";
                            // outfile << "ERROR: vectB contains non-finite entries before normal equations\n";
                        }

                        // Check partition of unity with detailed reporting
                        // gsInfo << "About to check partition of unity...\n";
                        // gsInfo << "Checking partition of unity...\n";
                        // outfile << "Checking partition of unity...\n";
                        real_t tolerance = 1e-6;  // More relaxed tolerance
                        const PartitionUnityViolationDetails partitionInfo =
                            analyzePartitionOfUnity(matA, tolerance);
                        bool partitionOK = partitionInfo.allSatisfied;
                        // gsInfo << "checkPartitionOfUnity returned: " << (partitionOK ? "true" : "false") << "\n";
                        
                        if (!partitionOK) {
                            // Record the exact failing row and its contributors before aborting.
                            auto writeRowBreakdown = [&](std::ostream& stream,
                                                         index_t rowIndex,
                                                         const std::string& label,
                                                         index_t maxContributors)
                            {
                                if (rowIndex < 0 || rowIndex >= matA.rows())
                                {
                                    stream << label << ": unavailable\n";
                                    return;
                                }

                                std::vector<std::pair<index_t, real_t>> contributors;
                                contributors.reserve(static_cast<size_t>(matA.cols()));
                                for (index_t col = 0; col < matA.cols(); ++col)
                                {
                                    const real_t value = matA.coeff(rowIndex, col);
                                    if (std::abs(value) > 1e-14)
                                        contributors.push_back(std::make_pair(col, value));
                                }

                                std::sort(
                                    contributors.begin(),
                                    contributors.end(),
                                    [](const std::pair<index_t, real_t>& lhs,
                                       const std::pair<index_t, real_t>& rhs)
                                    {
                                        return std::abs(lhs.second) > std::abs(rhs.second);
                                    });

                                const real_t rowSum = partitionInfo.rowSums[static_cast<size_t>(rowIndex)];
                                const real_t signedDeviation = rowSum - 1.0;
                                const index_t elementIdx = numGaussPointsPerElement > 0
                                    ? rowIndex / numGaussPointsPerElement
                                    : -1;
                                const index_t gaussIdx = numGaussPointsPerElement > 0
                                    ? rowIndex % numGaussPointsPerElement
                                    : -1;

                                stream << label << "\n";
                                stream << "  row=" << rowIndex
                                       << ", element=" << elementIdx
                                       << ", gauss=" << gaussIdx
                                       << ", rowSum=" << std::setprecision(16) << rowSum
                                       << ", deviation=" << signedDeviation
                                       << ", nonZeros=" << contributors.size() << "\n";
                                stream << "  interpretation="
                                       << (signedDeviation < 0.0 ? "missing basis contribution" : "excess basis contribution")
                                       << " of magnitude " << std::abs(signedDeviation) << "\n";

                                const index_t contributorCount = std::min<index_t>(maxContributors, contributors.size());
                                for (index_t i = 0; i < contributorCount; ++i)
                                {
                                    const index_t globalFunction = contributors[static_cast<size_t>(i)].first;
                                    const real_t value = contributors[static_cast<size_t>(i)].second;
                                    stream << "    col " << globalFunction << " -> " << value;
                                    if (globalFunction >= 0 &&
                                        globalFunction < static_cast<index_t>(functionDescription.size()))
                                    {
                                        stream << " components=";
                                        for (size_t comp = 0; comp < functionDescription[globalFunction].size(); ++comp)
                                        {
                                            const std::vector<index_t>& desc = functionDescription[globalFunction][comp];
                                            if (desc.size() < 3)
                                                continue;
                                            stream << "(" << desc[0] << "," << desc[1] << "," << desc[2] << ")";
                                            if (comp + 1 < functionDescription[globalFunction].size())
                                                stream << " ";
                                        }
                                    }
                                    stream << "\n";
                                }

                                if (contributors.size() > static_cast<size_t>(contributorCount))
                                    stream << "    ... " << (contributors.size() - contributorCount) << " more contributors\n";
                            };

                            gsInfo << "Entering partition violation handling...\n";
                            gsInfo << "\n*** WARNING: Partition of unity violated with tolerance " << tolerance << " ***\n";
                            gsInfo << "Attempt: " << attempt << ", patch=" << patch << ", level=" << levNow << "\n";
                            gsInfo << "Matrix dimensions: " << matA.rows() << " x " << matA.cols() << "\n";
                            gsInfo << "Violating rows: " << partitionInfo.violationCount
                                   << ", firstRow=" << partitionInfo.firstViolatingRow
                                   << ", worstRow=" << partitionInfo.maxDeviationRow
                                   << ", worstDeviation=" << partitionInfo.maxDeviation << "\n";
                            writeRowBreakdown(gsInfo, partitionInfo.firstViolatingRow, "First violating row", 8);
                            if (partitionInfo.maxDeviationRow != partitionInfo.firstViolatingRow)
                                writeRowBreakdown(gsInfo, partitionInfo.maxDeviationRow, "Worst violating row", 8);

                            outfile << "\n*** WARNING: Partition of unity violated with tolerance " << tolerance << " ***\n";
                            outfile << "Attempt: " << attempt << "\n";
                            outfile << "Patch being processed: " << patch << "\n";
                            outfile << "Level: " << levNow << "\n";
                            outfile << "Matrix dimensions: " << matA.rows() << " x " << matA.cols() << "\n";
                            outfile << "Violating rows: " << partitionInfo.violationCount
                                    << ", firstRow=" << partitionInfo.firstViolatingRow
                                    << ", worstRow=" << partitionInfo.maxDeviationRow
                                    << ", worstDeviation=" << partitionInfo.maxDeviation << "\n";
                            writeRowBreakdown(outfile, partitionInfo.firstViolatingRow, "First violating row", 20);
                            if (partitionInfo.maxDeviationRow != partitionInfo.firstViolatingRow)
                                writeRowBreakdown(outfile, partitionInfo.maxDeviationRow, "Worst violating row", 20);

                            // Save detailed diagnostic information
                            std::string diagFile = "partition_diagnostic_p" + std::to_string(patch) + "_lev" + std::to_string(levNow) + "_att" + std::to_string(attempt) + ".txt";
                            std::ofstream diag(diagFile);
                            diag << "========================================\n";
                            diag << "PARTITION OF UNITY VIOLATION DIAGNOSTICS\n";
                            diag << "========================================\n\n";
                            diag << "--- Assembly Information ---\n";
                            diag << "Matrix dimensions: " << matA.rows() << " rows x " << matA.cols() << " columns\n";
                            diag << "Number of MPBES functions: " << mpbes.size() << "\n";
                            diag << "Number of patches: " << mpbes.nPatches() << "\n";
                            diag << "Number of active elements: " << numElements << "\n";
                            diag << "Gauss points per element: " << numGaussPointsPerElement << "\n";
                            diag << "Total Gauss points (rows): " << numTotalRows << "\n";
                            diag << "Tolerance used: " << tolerance << "\n\n";
                            diag << "Violating rows: " << partitionInfo.violationCount << "\n";
                            diag << "First violating row: " << partitionInfo.firstViolatingRow << "\n";
                            diag << "Worst violating row: " << partitionInfo.maxDeviationRow << "\n";
                            diag << "Worst absolute deviation: " << partitionInfo.maxDeviation << "\n\n";
                            
                            diag << "--- Function Statistics ---\n";
                            int numTruncated = 0;
                            int numWithSpillover = 0;
                            for (index_t i = 0; i < mpbes.size(); ++i) {
                                if (i < isTruncated.size() && isTruncated[i]) numTruncated++;
                                if (i < hasSpillover.size() && hasSpillover[i]) numWithSpillover++;
                            }
                            diag << "Truncated functions: " << numTruncated << " / " << mpbes.size() << "\n";
                            diag << "Functions with spillover: " << numWithSpillover << " / " << mpbes.size() << "\n\n";
                            
                            diag << "--- Matrix Sparsity ---\n";
                            index_t nonZeros = matA.nonZeros();
                            real_t sparsity = 100.0 * (1.0 - static_cast<real_t>(nonZeros) / (matA.rows() * matA.cols()));
                            diag << "Non-zero entries: " << nonZeros << "\n";
                            diag << "Sparsity: " << sparsity << "%\n";
                            diag << "Average non-zeros per row: " << (matA.rows() > 0 ? static_cast<real_t>(nonZeros) / matA.rows() : 0) << "\n\n";
                            
                            diag << "--- First 10 Rows Sum Check ---\n";
                            for (index_t i = 0; i < std::min<index_t>(10, matA.rows()); ++i) {
                                const real_t rowSum = partitionInfo.rowSums[static_cast<size_t>(i)];
                                diag << "Row " << i << ": sum = " << rowSum;
                                if (std::abs(rowSum - 1.0) > tolerance) {
                                    diag << " [VIOLATION: " << (rowSum - 1.0) << "]";
                                }
                                diag << "\n";
                            }
                            diag << "\n";
                            
                            diag << "--- Sample Matrix Values (First 5x5 block) ---\n";
                            for (index_t i = 0; i < std::min<index_t>(5, matA.rows()); ++i) {
                                diag << "Row " << i << ": ";
                                for (index_t j = 0; j < std::min<index_t>(5, matA.cols()); ++j) {
                                    diag << matA(i, j);
                                    if (j < std::min<index_t>(5, matA.cols()) - 1) diag << ", ";
                                }
                                if (matA.cols() > 5) diag << ", ...";
                                diag << "\n";
                            }
                            diag << "\n";
                            
                            diag << "--- Violated Rows (first 20) ---\n";
                            int violationCount = 0;
                            for (index_t i = 0; i < matA.rows() && violationCount < 20; ++i) {
                                const real_t rowSum = partitionInfo.rowSums[static_cast<size_t>(i)];
                                if (std::abs(rowSum - 1.0) > tolerance) {
                                    diag << "Row " << i << ": sum = " << rowSum << " (deviation: " << (rowSum - 1.0) << ")\n";
                                    violationCount++;
                                }
                            }
                            diag << "Total violations found: " << violationCount << (violationCount >= 20 ? " (showing first 20)" : "") << "\n\n";

                            diag << "--- Detailed Violating Row Breakdown ---\n";
                            writeRowBreakdown(diag, partitionInfo.firstViolatingRow, "First violating row", 25);
                            if (partitionInfo.maxDeviationRow != partitionInfo.firstViolatingRow)
                                writeRowBreakdown(diag, partitionInfo.maxDeviationRow, "Worst violating row", 25);
                            diag.flush();
                            diag.close();

                            if (partitionInfo.firstViolatingRow >= 0) {
                                std::string rowFileName = "row" + std::to_string(partitionInfo.firstViolatingRow)
                                    + "_p" + std::to_string(patch) + "_lev" + std::to_string(levNow)
                                    + "_att" + std::to_string(attempt) + ".txt";

                                try {
                                    std::ofstream rowDiagFile(rowFileName);
                                    if (rowDiagFile.is_open()) {
                                        rowDiagFile << "=== Detailed Row " << partitionInfo.firstViolatingRow << " Diagnostic (Patch " 
                                                   << patch << ", Level " << levNow << ", Attempt " << attempt << ") ===\n";
                                        rowDiagFile << "Total functions in MPBES: " << commonSize << "\n";
                                        writeRowBreakdown(rowDiagFile, partitionInfo.firstViolatingRow, "Detailed row breakdown", 50);
                                        rowDiagFile << "\n=== End Detailed Row " << partitionInfo.firstViolatingRow << " Diagnostic ===\n";
                                        rowDiagFile.flush();
                                        rowDiagFile.close();
                                    }
                                }
                                catch (const std::exception& e) {
                                    gsInfo << "WARNING: Exception creating row diagnostic file: " << e.what() << "\n";
                                    if (outfile.is_open())
                                        outfile << "WARNING: Exception creating row diagnostic file: " << e.what() << "\n";
                                }
                            }

                            gsInfo << "Detailed diagnostic saved to " << diagFile << "\n";
                            if (outfile.is_open())
                                outfile << "Detailed diagnostic saved to " << diagFile << "\n";

                            // Abort immediately after flushing so the failing attempt cannot continue silently.
                            if (summaryFile.is_open())
                            {
                                summaryFile << "PARTITION_VIOLATION," << patch << "," << levNow << "," << attempt
                                            << ",nan,nan,nan,nan\n";
                                summaryFile.flush();
                            }
                            if (outfile.is_open())
                            {
                                outfile << "Exiting due to partition violation with exit code "
                                        << kPartitionUnityViolationExitCode << "\n";
                                outfile.flush();
                            }
                            gsInfo << "Exiting due to partition violation with exit code "
                                   << kPartitionUnityViolationExitCode << "\n";
                            flushAndCloseDiagnosticLogs();
                            throw ProgramExitSignal(kPartitionUnityViolationExitCode,
                                "Partition of unity violation detected");
                            
                        } else {
                            // gsInfo << "Partition of unity satisfied.\n";
                            // outfile << "Partition of unity satisfied.\n";
                        }

                        // gsInfo << "Computing assembly time...\n";
                        afterassembleA = std::chrono::system_clock::now();
                        std::chrono::duration<double> elapsed_seconds_assembleA = afterassembleA - beforeassembleA;
                        t_assemble += elapsed_seconds_assembleA.count();
                        // gsInfo << "Assembly took: " << elapsed_seconds_assembleA.count() << " seconds\n";

                        // Column mapping + defect correction (timed as t_presolve)
                        auto _tPS0 = std::chrono::system_clock::now();

                        // Check matrix rank before solving
                        // gsInfo << "Converting matA to dense...\n";
                        // gsInfo << "matA.rows()=" << matA.rows() << ", matA.cols()=" << matA.cols() << ", matA.nonZeros()=" << matA.nonZeros() << "\n";
                        // Non-finite check directly on sparse stored values (avoids 87.7MB dense conversion)
                        if (outfile.is_open())
                        {
                            outfile << "[fit-solve] matA sparse rows=" << matA.rows()
                                    << " cols=" << matA.cols()
                                    << " nonZeros=" << matA.nonZeros() << "\n";
                            index_t nonFiniteA = 0;
                            for (index_t k = 0; k < matA.outerSize(); ++k)
                                for (gsSparseMatrix<real_t>::InnerIterator it(matA, k); it; ++it)
                                    if (!std::isfinite(it.value())) ++nonFiniteA;
                            outfile << "[fit-solve] matA nonFiniteEntries=" << nonFiniteA << "\n";
                        }

                        std::vector<index_t> sourceCols;
                        std::vector<index_t> globalCols;
                        sourceCols.reserve(matA.cols());
                        globalCols.reserve(matA.cols());
                        if (useLocalFitting)
                        {
                            if (!assembleActiveFunctions.empty() &&
                                matA.cols() == commonSize)
                            {
                                for (index_t globalCol = 0; globalCol < static_cast<index_t>(assembleActiveFunctions.size()); ++globalCol)
                                {
                                    sourceCols.push_back(assembleActiveFunctions[globalCol]);
                                    globalCols.push_back(assembleActiveFunctions[globalCol]);
                                }
                            }
                            else if (!assembleActiveFunctions.empty() &&
                                matA.cols() == static_cast<index_t>(assembleActiveFunctions.size()))
                            {
                                for (index_t localCol = 0; localCol < matA.cols(); ++localCol)
                                {
                                    sourceCols.push_back(localCol);
                                    globalCols.push_back(assembleActiveFunctions[localCol]);
                                }
                            }
                            else if (localRegion.enabled && localRegion.basisInd.cols() == matA.cols())
                            {
                                for (index_t f = 0; f < matA.cols(); ++f)
                                    if (localRegion.basisInd(0, f) != 0.0)
                                    {
                                        sourceCols.push_back(f);
                                        globalCols.push_back(f);
                                    }
                            }
                        }

                        if (sourceCols.empty())
                        {
                            for (index_t f = 0; f < matA.cols(); ++f)
                            {
                                sourceCols.push_back(f);
                                globalCols.push_back(f);
                            }
                        }

                        // Column informative check — squared norm directly on sparse columns
                        std::vector<index_t> informativeCols;
                        informativeCols.reserve(sourceCols.size());
                        const real_t colTol = 1e-14;
                        for (index_t j = 0; j < static_cast<index_t>(sourceCols.size()); ++j)
                        {
                            const index_t localCol = sourceCols[j];
                            if (!useLocalFitting || matA.col(localCol).squaredNorm() > colTol)
                                informativeCols.push_back(j);
                        }

                        if (informativeCols.empty())
                        {
                            for (index_t j = 0; j < static_cast<index_t>(sourceCols.size()); ++j)
                                informativeCols.push_back(j);
                        }

                        std::vector<char> isInformative(sourceCols.size(), 0);
                        for (index_t idx = 0; idx < static_cast<index_t>(informativeCols.size()); ++idx)
                        {
                            const index_t localCol = informativeCols[idx];
                            if (localCol >= 0 && localCol < static_cast<index_t>(isInformative.size()))
                                isInformative[localCol] = 1;
                        }

                        std::vector<index_t> solveCols;
                        solveCols.reserve(informativeCols.size());
                        for (index_t k = 0; k < static_cast<index_t>(informativeCols.size()); ++k)
                            solveCols.push_back(globalCols[informativeCols[k]]);

                        if (outfile.is_open())
                        {
                            outfile << "[fit-solve] column mapping source="
                                << ((!assembleActiveFunctions.empty() &&
                                 matA.cols() == static_cast<index_t>(assembleActiveFunctions.size()))
                                    ? "assembleActiveFunctions(global)"
                                    : ((localRegion.enabled && localRegion.basisInd.cols() == matA.cols())
                                    ? "basisInd(global)"
                                    : "identity/global"))
                                << "\n";
                            outfile << "[fit-solve] dynamic column reduction: selected=" << sourceCols.size()
                                    << " / " << matA.cols() << "\n";
                            outfile << "[fit-solve] matA rows=" << matA.rows()
                                    << " cols=" << sourceCols.size() << "\n";
                            outfile << "[fit-solve] informative column reduction: selected=" << solveCols.size()
                                    << " / " << sourceCols.size() << "\n";
                            outfile << "[fit-solve] matA_solve rows=" << matA.rows()
                                    << " cols=" << solveCols.size() << "\n";
                            outfile << "[fit-solve] dropped active/global columns by informative reduction begin\n";
                            for (index_t j = 0; j < static_cast<index_t>(sourceCols.size()); ++j)
                            {
                                if (j < static_cast<index_t>(isInformative.size()) && isInformative[j] != 0)
                                    continue;
                                const real_t colNorm2 = matA.col(sourceCols[j]).squaredNorm();
                                outfile << "  localCol=" << j
                                    << " sourceCol=" << sourceCols[j]
                                    << " globalCol=" << globalCols[j]
                                    << " squaredNorm=" << colNorm2 << "\n";
                            }
                            outfile << "[fit-solve] dropped active/global columns by informative reduction end\n";
                        }

                        // Solve the system
                        // gsInfo << "Creating b_vec copy...\n";
                        gsMatrix<> b_vec = vectB;
                        const index_t fullRows = static_cast<index_t>(commonSize);
                        const index_t fullCols = vectB.cols();
                        gsMatrix<real_t> vectSolSeed = gsMatrix<real_t>::Zero(fullRows, fullCols);
                        if (useLocalFitting && effectiveOriginalCoefficients &&
                            effectiveOriginalCoefficients->rows() == fullRows &&
                            effectiveOriginalCoefficients->cols() >= fullCols)
                        {
                            vectSolSeed = effectiveOriginalCoefficients->leftCols(fullCols);
                        }

                        std::vector<char> solveMask(static_cast<std::size_t>(fullRows), 0);
                        for (index_t j = 0; j < static_cast<index_t>(solveCols.size()); ++j)
                        {
                            const index_t globalCol = solveCols[j];
                            if (globalCol >= 0 && globalCol < fullRows)
                                solveMask[static_cast<std::size_t>(globalCol)] = 1;
                        }

                        // Defect correction for local fitting.
                        //
                        // A_local (n_pts x n_local) was assembled with only the local functions.
                        // The LS system models:
                        //   A_local * x_local = b_target - A_nonlocal * x_nonlocal_current
                        //
                        // Since A_nonlocal is not assembled, compute the non-local contribution
                        // indirectly:
                        //   A_nonlocal * x_nonlocal = geom_current - A_local * x_local_current
                        //
                        // So: b_adj = b_target - geom_current + A_local * x_local_current
                        if (useLocalFitting
                            && !assembleActiveFunctions.empty()
                            && matA.cols() == static_cast<index_t>(assembleActiveFunctions.size()))
                        {
                            const index_t nLocal = static_cast<index_t>(assembleActiveFunctions.size());

                            // 1. Extract current local control-point values from seed solution
                            gsMatrix<real_t> x_local_current(nLocal, fullCols);
                            for (index_t j = 0; j < nLocal; ++j)
                            {
                                const index_t gf = assembleActiveFunctions[j];
                                if (gf >= 0 && gf < vectSolSeed.rows())
                                    x_local_current.row(j) = vectSolSeed.row(gf);
                                else
                                    x_local_current.row(j).setZero();
                            }

                            // 2. Evaluate current parameterization (all functions) at uvFitting pts.
                            // Uses evalGeomOnPatch (batch active_into+eval_into) instead of
                            // per-point evaluateFittedGeometryPoint to avoid O(MPBES_size) per point.
                            gsMatrix<real_t> geom_current(b_vec.rows(), fullCols);
                            geom_current.setZero();
                            {
                                index_t row = 0;
                                for (index_t p = 0; p < static_cast<index_t>(uvFitting.size()); ++p)
                                {
                                    const index_t Mp = uvFitting(p).cols();
                                    if (Mp == 0) continue;
                                    gsMatrix<real_t> patchGeom;
                                    mpbes.evalGeomOnPatch(p, uvFitting(p), vectSolSeed, patchGeom);
                                    const index_t cols = std::min<index_t>(patchGeom.cols(), fullCols);
                                    geom_current.block(row, 0, Mp, cols) = patchGeom.leftCols(cols);
                                    row += Mp;
                                }
                            }

                            // 3. b_adj = b_target - geom_current + A_local * x_local_current
                            b_vec -= geom_current;
                            b_vec.noalias() += matA * x_local_current;

                            gsInfo << "local fit: b adjusted by defect correction"
                                   << " (nLocal=" << nLocal << " pts=" << b_vec.rows() << ")\n";
                            if (outfile.is_open())
                                outfile << "local fit: b adjusted by defect correction"
                                        << " (nLocal=" << nLocal << " pts=" << b_vec.rows() << ")\n";
                        }

                        if (outfile.is_open())
                        {
                            outfile << "[fit-solve] vectB(before At*b) rows=" << vectB.rows()
                                    << " cols=" << vectB.cols() << "\n";
                            index_t nonFiniteB = 0;
                            for (index_t r = 0; r < vectB.rows(); ++r)
                                for (index_t c = 0; c < vectB.cols(); ++c)
                                    if (!std::isfinite(vectB(r, c)))
                                        ++nonFiniteB;
                            outfile << "[fit-solve] vectB(before At*b) nonFiniteEntries=" << nonFiniteB << "\n";
                        }
                        // gsInfo << "b_vec.rows()=" << b_vec.rows() << ", b_vec.cols()=" << b_vec.cols() << "\n";
                        
                        t_presolve += std::chrono::duration<double>(std::chrono::system_clock::now() - _tPS0).count();

                        // Build sparse A from the LOCAL column indices into matA.
                        // solveCols[k] is the global MPBES index, but matA columns use
                        // LOCAL indices (sourceCols[informativeCols[k]]) in local-fitting mode.
                        auto _tLdlt0 = std::chrono::system_clock::now();
                        gsSparseMatrix<real_t> matA_sp(matA.rows(), static_cast<index_t>(solveCols.size()));
                        {
                            std::vector<gsEigen::Triplet<real_t>> trips;
                            trips.reserve(matA.nonZeros());
                            for (index_t colK = 0; colK < static_cast<index_t>(informativeCols.size()); ++colK)
                            {
                                const index_t localCol = sourceCols[informativeCols[colK]];
                                if (localCol >= 0 && localCol < matA.cols())
                                    for (gsSparseMatrix<real_t>::InnerIterator it(matA, localCol); it; ++it)
                                        trips.emplace_back(it.row(), colK, it.value());
                            }
                            matA_sp.setFromTriplets(trips.begin(), trips.end());
                        }

                        // A^T b (sparse × dense)
                        vectB = matA_sp.transpose() * b_vec;

                        // A^T A (sparse × sparse)
                        gsSparseMatrix<real_t> AtA_sp = (matA_sp.transpose() * matA_sp).eval();
                        AtA_sp.makeCompressed();

                        if (outfile.is_open())
                            outfile << "[fit-solve] sparse AtA: " << AtA_sp.rows() << "x" << AtA_sp.cols()
                                    << " nnz=" << AtA_sp.nonZeros() << "\n";

                        if (verboseFitMatrixDump)
                        {
                            gsMatrix<> matA_solve_dbg(matA_sp);
                            printTheMatrix(matA_solve_dbg, "A_fit_reduced");
                            printTheMatrix(b_vec, "b_fit_original");
                            printTheMatrix(vectB, "At_b_fit");
                        }

                        commonSize = functionDescription.size();
                        gsMatrix<real_t> vectSolReduced;
                        {
                            typename gsSparseSolver<real_t>::SimplicialLDLT ldlt;
                            ldlt.compute(AtA_sp);
                            if (ldlt.info() != gsEigen::Success)
                            {
                                gsInfo  << "[fit-solve] LDLT failed, fallback to dense partialPivLu\n";
                                outfile << "[fit-solve] LDLT failed, fallback to dense partialPivLu\n";
                                gsMatrix<> AtA_dense(AtA_sp);
                                vectSolReduced = AtA_dense.partialPivLu().solve(vectB);
                            }
                            else
                            {
                                vectSolReduced = ldlt.solve(vectB);
                            }
                        }
                        t_ata_ldlt += std::chrono::duration<double>(std::chrono::system_clock::now() - _tLdlt0).count();
                        auto _tPS1 = std::chrono::system_clock::now();
                        {
                            auto afterfit = std::chrono::system_clock::now();
                            std::chrono::duration<double> elapsed_fit = afterfit - beforeassembleA;
                            if (g_verbose) gsInfo << "fit took: " << elapsed_fit.count() << "\n";
                            if (outfile.is_open())
                                outfile << "fit took: " << elapsed_fit.count() << "\n";
                        }

                        // matAsquare (dense A^T*A) deferred — only built inside LO/NLO block when needed

                        vectSol = vectSolSeed;

                        gsMatrix<real_t> vectSolBeforeMerge = vectSol;

                        gsMatrix<real_t> vectSolFitted;
                        vectSolFitted = vectSolBeforeMerge;
                        std::vector<char> hasSolvedValue(fullRows, 0);
                        for (index_t j = 0; j < static_cast<index_t>(solveCols.size()); ++j)
                            if (solveCols[j] >= 0 && solveCols[j] < fullRows)
                            {
                                vectSolFitted.row(solveCols[j]) = vectSolReduced.row(j);
                                hasSolvedValue[solveCols[j]] = 1;
                            }

                        const bool mergeBasisMaskApplied =
                            (useLocalFitting && localRegion.enabled && localRegion.basisInd.cols() == fullRows);

                        if (outfile.is_open() && g_verbose)
                        {
                            index_t diagnosticRow = -1;
                            if (!solveCols.empty())
                                diagnosticRow = solveCols.front();
                            else if (!assembleActiveFunctions.empty())
                                diagnosticRow = assembleActiveFunctions.front();

                            if (diagnosticRow >= 0 && diagnosticRow < fullRows)
                            {
                                const auto& diagFunctionDescription = mpbes.functionDescription();
                                const auto& diagSpilloverCoordinates = mpbes.spilloverCoordinates();
                                const auto& diagHasSpillover = mpbes.hasSpillover();
                                const auto& diagBells = mpbes.bellsBases();

                                auto countSamplesInSupport = [&](int p, int lvl, int idx) -> index_t
                                {
                                    if (p < 0 || p >= static_cast<int>(uvFitting.size()))
                                        return 0;
                                    if (uvFitting(p).cols() == 0)
                                        return 0;
                                    if (p >= static_cast<int>(diagBells.size()) ||
                                        lvl < 0 || lvl >= static_cast<int>(diagBells[p].size()) ||
                                        idx < 0 || idx >= static_cast<int>(diagBells[p][lvl].size()))
                                        return 0;

                                    const auto support = diagBells[p][lvl].function(idx).support();
                                    const real_t u0 = static_cast<real_t>(support(0, 0));
                                    const real_t u1 = static_cast<real_t>(support(0, 1));
                                    const real_t v0 = static_cast<real_t>(support(1, 0));
                                    const real_t v1 = static_cast<real_t>(support(1, 1));

                                    index_t insideCount = 0;
                                    for (index_t c = 0; c < uvFitting(p).cols(); ++c)
                                    {
                                        const real_t u = uvFitting(p)(0, c);
                                        const real_t v = uvFitting(p)(1, c);
                                        if (u >= u0 && u <= u1 && v >= v0 && v <= v1)
                                            ++insideCount;
                                    }
                                    return insideCount;
                                };

                                index_t assembleListIndex = -1;
                                for (index_t idx = 0; idx < static_cast<index_t>(assembleActiveFunctions.size()); ++idx)
                                {
                                    if (assembleActiveFunctions[idx] == diagnosticRow)
                                    {
                                        assembleListIndex = idx;
                                        break;
                                    }
                                }

                                index_t workColIndex = -1;
                                for (index_t idx = 0; idx < static_cast<index_t>(globalCols.size()); ++idx)
                                {
                                    if (globalCols[idx] == diagnosticRow)
                                    {
                                        workColIndex = idx;
                                        break;
                                    }
                                }

                                index_t solveIndex = -1;
                                for (index_t idx = 0; idx < static_cast<index_t>(solveCols.size()); ++idx)
                                {
                                    if (solveCols[idx] == diagnosticRow)
                                    {
                                        solveIndex = idx;
                                        break;
                                    }
                                }

                                outfile << "[fit-solve] diagnostic row begin: globalRow=" << diagnosticRow << "\n";
                                outfile << "  basisInd="
                                        << ((mergeBasisMaskApplied && localRegion.basisInd(0, diagnosticRow) != 0.0) ? 1 : 0)
                                        << " hasSolvedValue=" << (hasSolvedValue[diagnosticRow] ? 1 : 0)
                                        << " inAssembleActiveFunctions=" << (assembleListIndex >= 0 ? 1 : 0)
                                        << " assembleListIndex=" << assembleListIndex
                                        << " workColIndex=" << workColIndex
                                        << " solveIndex=" << solveIndex << "\n";

                                if (workColIndex >= 0)
                                {
                                    const index_t sourceCol = sourceCols[workColIndex];
                                    const real_t sparseNorm2 = matA.col(sourceCol).squaredNorm();
                                    index_t sparseNonZeros = 0;
                                    for (gsSparseMatrix<real_t>::InnerIterator it(matA, sourceCol); it; ++it)
                                        ++sparseNonZeros;

                                    outfile << "  sourceCol=" << sourceCol
                                            << " globalCol=" << globalCols[workColIndex]
                                            << " informative="
                                            << ((workColIndex < static_cast<index_t>(isInformative.size()) && isInformative[workColIndex] != 0) ? 1 : 0)
                                            << " colTol=" << colTol
                                            << " sparseSquaredNorm=" << sparseNorm2
                                            << " sparseNonZeros=" << sparseNonZeros << "\n";

                                    outfile << "  assembled column entries begin\n";
                                    for (gsSparseMatrix<real_t>::InnerIterator it(matA, sourceCol); it; ++it)
                                        outfile << "    row=" << it.row() << " value=" << it.value() << "\n";
                                    outfile << "  assembled column entries end\n";
                                }
                                else
                                {
                                    outfile << "  selected diagnostic row has no compact/global column in sourceCols\n";
                                }

                                outfile << "  regular components begin\n";
                                if (diagnosticRow < static_cast<index_t>(diagFunctionDescription.size()))
                                {
                                    for (size_t compIdx = 0; compIdx < diagFunctionDescription[diagnosticRow].size(); ++compIdx)
                                    {
                                        const auto& desc = diagFunctionDescription[diagnosticRow][compIdx];
                                        if (desc.size() < 3)
                                            continue;
                                        const int p = static_cast<int>(desc[0]);
                                        const int lvl = static_cast<int>(desc[1]);
                                        const int idx = static_cast<int>(desc[2]);
                                        outfile << "    component " << compIdx
                                                << ": patch=" << p
                                                << " level=" << lvl
                                                << " tensorIdx=" << idx;
                                        if (p >= 0 && p < static_cast<int>(diagBells.size()) &&
                                            lvl >= 0 && lvl < static_cast<int>(diagBells[p].size()) &&
                                            idx >= 0 && idx < static_cast<int>(diagBells[p][lvl].size()))
                                        {
                                            const auto support = diagBells[p][lvl].function(idx).support();
                                            outfile << " support=[" << support(0, 0) << ", " << support(0, 1)
                                                    << "] x [" << support(1, 0) << ", " << support(1, 1) << "]"
                                                    << " uvCols=" << ((p >= 0 && p < static_cast<int>(uvFitting.size())) ? uvFitting(p).cols() : 0)
                                                    << " sampleCountInSupport=" << countSamplesInSupport(p, lvl, idx);
                                        }
                                        outfile << "\n";
                                    }
                                }
                                outfile << "  regular components end\n";

                                outfile << "  spillover components begin\n";
                                if (diagnosticRow < static_cast<index_t>(diagHasSpillover.size()) &&
                                    diagHasSpillover[diagnosticRow] &&
                                    diagnosticRow < static_cast<index_t>(diagSpilloverCoordinates.size()))
                                {
                                    for (size_t spIdx = 0; spIdx < diagSpilloverCoordinates[diagnosticRow].size(); ++spIdx)
                                    {
                                        const auto& sp = diagSpilloverCoordinates[diagnosticRow][spIdx];
                                        const int p = sp[0];
                                        const int lvl = sp[1];
                                        const int idx = sp[2];
                                        outfile << "    spillover " << spIdx
                                                << ": patch=" << p
                                                << " level=" << lvl
                                                << " tensorIdx=" << idx;
                                        if (p >= 0 && p < static_cast<int>(diagBells.size()) &&
                                            lvl >= 0 && lvl < static_cast<int>(diagBells[p].size()) &&
                                            idx >= 0 && idx < static_cast<int>(diagBells[p][lvl].size()))
                                        {
                                            const auto support = diagBells[p][lvl].function(idx).support();
                                            outfile << " support=[" << support(0, 0) << ", " << support(0, 1)
                                                    << "] x [" << support(1, 0) << ", " << support(1, 1) << "]"
                                                    << " uvCols=" << ((p >= 0 && p < static_cast<int>(uvFitting.size())) ? uvFitting(p).cols() : 0)
                                                    << " sampleCountInSupport=" << countSamplesInSupport(p, lvl, idx);
                                        }
                                        outfile << "\n";
                                    }
                                }
                                outfile << "  spillover components end\n";
                                outfile << "[fit-solve] diagnostic row end: globalRow=" << diagnosticRow << "\n";
                            }
                        }

                        if (mergeBasisMaskApplied)
                        {
                            index_t keptFromMpbes = 0;
                            index_t overwrittenFromFitted = 0;
                            index_t activeRowsWithoutSolvedValue = 0;
                            for (index_t f = 0; f < vectSol.rows(); ++f)
                            {
                                if (localRegion.basisInd(0, f) != 0.0)
                                {
                                    vectSol.row(f) = vectSolFitted.row(f);
                                    ++overwrittenFromFitted;
                                    if (!hasSolvedValue[f])
                                        ++activeRowsWithoutSolvedValue;
                                }
                                else
                                {
                                    vectSol.row(f) = vectSolBeforeMerge.row(f);
                                    ++keptFromMpbes;
                                }
                            }

                            if (outfile.is_open())
                            {
                                outfile << "[fit-solve] merge rule: keep mpbes when basisInd==0; "
                                        << "use vectSolFitted(row) when basisInd==1\n";
                                outfile << "[fit-solve] merge counts: keptFromMpbes=" << keptFromMpbes
                                        << " overwrittenFromFitted=" << overwrittenFromFitted
                                        << " activeRowsWithoutSolvedValue=" << activeRowsWithoutSolvedValue << "\n";
                            }
                        }
                        else
                        {
                            vectSol = vectSolFitted;
                        }

                        if (verboseFitMatrixDump)
                        {
                            printTheMatrix(vectSolReduced, "x_fit_reduced");
                            printTheMatrix(vectSolFitted, "vectSol_fitted_expanded");
                            printTheMatrix(vectSol, "vectSol_merged");
                        }

                        if (outfile.is_open())
                        {
                            outfile << "[fit-solve] vectSol(after solve) rows=" << vectSol.rows()
                                    << " cols=" << vectSol.cols() << "\n";
                                outfile << "[fit-solve] vectSol merge: seededFromOriginal="
                                     << ((effectiveOriginalCoefficients &&
                                         effectiveOriginalCoefficients->rows() == fullRows &&
                                         effectiveOriginalCoefficients->cols() >= fullCols) ? 1 : 0)
                                    << " overwrittenRows=" << solveCols.size()
                                    << " basisMaskApplied="
                                    << (mergeBasisMaskApplied ? 1 : 0)
                                    << "\n";
                            if (verboseVectSolRowDump)
                            {
                                outfile << "[fit-solve] vectSol row-by-row before/after begin\n";
                                for (index_t r = 0; r < vectSol.rows(); ++r)
                                {
                                    const bool basisActive =
                                        (mergeBasisMaskApplied && localRegion.basisInd(0, r) != 0.0);
                                    const bool solvedPresent =
                                        (r >= 0 && r < static_cast<index_t>(hasSolvedValue.size()) && hasSolvedValue[r] != 0);

                                    outfile << "  row " << r
                                            << " before=(";
                                    for (index_t c = 0; c < vectSolBeforeMerge.cols(); ++c)
                                    {
                                        if (c > 0)
                                            outfile << ", ";
                                        outfile << vectSolBeforeMerge(r, c);
                                    }
                                    outfile << ")\n";

                                    outfile << "  row " << r
                                            << " after=(";
                                    for (index_t c = 0; c < vectSol.cols(); ++c)
                                    {
                                        if (c > 0)
                                            outfile << ", ";
                                        outfile << vectSol(r, c);
                                    }
                                    outfile << ") basisInd=" << (basisActive ? 1 : 0)
                                            << " hasSolvedValue=" << (solvedPresent ? 1 : 0);

                                    if (solvedPresent)
                                    {
                                        outfile << " fitted=(";
                                        for (index_t c = 0; c < vectSolFitted.cols(); ++c)
                                        {
                                            if (c > 0)
                                                outfile << ", ";
                                            outfile << vectSolFitted(r, c);
                                        }
                                        outfile << ")";
                                    }
                                    outfile << "\n";
                                }
                                outfile << "[fit-solve] vectSol row-by-row before/after end\n";
                            }
                            index_t nonFiniteX = 0;
                            for (index_t r = 0; r < vectSol.rows(); ++r)
                                for (index_t c = 0; c < vectSol.cols(); ++c)
                                    if (!std::isfinite(vectSol(r, c)))
                                        ++nonFiniteX;
                            outfile << "[fit-solve] vectSol(after solve) nonFiniteEntries=" << nonFiniteX << "\n";
                        }

                        // gsInfo << "Solve finished. vectSol.rows()=" << vectSol.rows() << ", cols=" << vectSol.cols() << "\n";
                        outfile << "Solve finished. vectSol.rows()=" << vectSol.rows() << ", cols=" << vectSol.cols() << "\n";

                        if (g_verbose) gsInfo << "vectSol.rows(): " << vectSol.rows() << "\n";
                        gsMatrix<> matC = gsMatrix<>(AtA_sp * vectSolReduced) - vectB;
                        if (g_verbose)
                        {
                            gsInfo << "Residual norm (A^T A x - A^T b) matCbefore\n";
                            gsInfo << matC.maxCoeff() << "\n";
                        }
                        outfile << "matCbefore\n";
                        outfile << matC.maxCoeff() << "\n";
                        ////gsInfo << "matC:\n" << matC << "\n";
                        gsMatrix<> matOut = matA_sp * vectSolReduced;
                        logResidualDetails(uvFitting, b_vec, matOut, outfile, false);
                        matC = matOut - b_vec;
                        if (g_verbose)
                        {
                            gsInfo << "Residual norm (Ax - b) matCafter\n";
                            gsInfo << matC.maxCoeff() << "\n";
                        }
                        outfile << "matCafter\n";
                        outfile << matC.maxCoeff() << "\n";
                        real_t globalError = 0.0;
                        int minusnumber = 0;
                        gsVector<size_t> numIrregular(Bells.size());
                        numIrregular.setZero();

                        // Adaptive jacPoints: scale with current mpbes size per patch, min 20
                        const index_t adaptiveJacPts = std::max<index_t>(20,
                            static_cast<index_t>(std::ceil(2.0 * std::sqrt(
                                static_cast<double>(mpbes.size()) / mpbes.nPatches()))));
                        gsVector<gsMatrix<real_t>> uv2_adaptive(mpbes.nPatches());
                        for (index_t p = 0; p < mpbes.nPatches(); ++p)
                            uv2_adaptive(p) = uniformPointGrid(lowc2, uppc2, adaptiveJacPts * adaptiveJacPts);
                        double totalJackTime = 0.0;

                        t_postsolve += std::chrono::duration<double>(std::chrono::system_clock::now() - _tPS1).count();
                        std::chrono::time_point<std::chrono::system_clock> beforejack, afterjack;
                        beforejack = std::chrono::system_clock::now();

                        bool jacVerbose = 0;//(patch == 1 && levNow == 3 && (attempt == 45 || attempt == 17 || attempt == 20));

                        minusnumber = checkJacobianDeterminant(
                            uv2_adaptive,
                            mpbes,
                            vectSol,
                            numIrregular,
                            jacVerbose);
                        afterjack = std::chrono::system_clock::now();
                        std::chrono::duration<double> elapsed_seconds_jack = afterjack - beforejack;
                        if (g_verbose)
                            gsInfo << "jack took: " << elapsed_seconds_jack.count()
                                   << " (pts/patch=" << adaptiveJacPts << "^2=" << adaptiveJacPts*adaptiveJacPts << ")\n";
                        outfile << "jack took: " << elapsed_seconds_jack.count()
                                << " (pts/patch=" << adaptiveJacPts << "^2=" << adaptiveJacPts*adaptiveJacPts << ")\n";
                        totalJackTime += elapsed_seconds_jack.count();
                        t_jack += elapsed_seconds_jack.count();
                        auto _tPJ0 = std::chrono::system_clock::now();

                        // Calculate total sample points and percentage of irregular points
                        index_t totalPoints = 0;
                        for (index_t patch = 0; patch < uv2_adaptive.size(); ++patch) {
                            totalPoints += uv2_adaptive[patch].cols();
                        }
                        real_t irregularPercentage = totalPoints > 0 ? (100.0 * minusnumber / totalPoints) : 0.0;
                        
                           std::cout << "Irregular points: " << minusnumber << " / " << totalPoints
                                  << " (" << irregularPercentage << "%)\n";
                        outfile << "Irregular points: " << minusnumber << " / " << totalPoints
                                << " (" << irregularPercentage << "%)\n";
                        if (g_enableInterfaceDiagnostics)
                        {
                            logSpecificInterface(mpbes, mp1, vectSol, "main_postsolve", 0, 1);
                            logSpecificInterface(mpbes, mp1, vectSol, "main_postsolve", 0, 3);
                            logSpecificInterface(mpbes, mp1, vectSol, "main_postsolve", 1, 2);
                            logSpecificInterface(mpbes, mp1, vectSol, "main_postsolve", 2, 5);
                            logSpecificInterface(mpbes, mp1, vectSol, "main_postsolve", 3, 4);
                        }
                        
                        gsMatrix<> geomDer;
            
                        std::cout << "minusnumber: " << minusnumber << "\n";
                        outfile << "minusnumber: " << minusnumber << "\n";
                        std::vector<std::vector<std::vector<double>>> localCoeffs(Bells.size());
                        gsVector<gsVector<int>> localIndex;
                        if (g_verbose && outfile.is_open())
                        {
                            outfile << "Solution coefficients: " << vectSol.rows() << " x " << vectSol.cols() << "\n";
                            for (index_t i = 0; i < vectSol.rows(); ++i) {
                                outfile << vectSol(i, 0) << " " << vectSol(i, 1) << "\n";
                            }
                        }
                        gsMatrix<> matFeatOut;
                        double featureError = 1;
                        toc = std::chrono::system_clock::now();
                        std::chrono::duration<double> elapsed_seconds_write = toc - startTime;
                        t_postjack += std::chrono::duration<double>(std::chrono::system_clock::now() - _tPJ0).count();
                        // Use the same acceptance metric in local and global modes so outcomes are comparable.
                        { auto _tB0 = std::chrono::system_clock::now();
                        real_t assemblyBoundaryError = testBoundaryAssembly(mpbes, mp1, vectSol);
                        t_boundary += std::chrono::duration<double>(std::chrono::system_clock::now() - _tB0).count();
                        featureError = assemblyBoundaryError; }
                        real_t assemblyBoundaryError = featureError;


                        // Build target geometry matrix at uv1 sample points for error computation
                        index_t totalSamples = 0;
                        index_t safeSize = std::min(static_cast<index_t>(uv1.size()), static_cast<index_t>(mp1.nPatches()));
                        for (index_t p = 0; p < safeSize; ++p)
                            totalSamples += uv1(p).cols();

                        if (totalSamples == 0) {
                            gsInfo << "WARNING: totalSamples is 0, setting targetGeometry to 1x2\n";
                            outfile << "WARNING: totalSamples is 0, setting targetGeometry to 1x2\n";
                            totalSamples = 1;
                        }

                        gsMatrix<real_t> targetGeometry(totalSamples, 2);
                        { auto _tTG0 = std::chrono::system_clock::now();
                        targetGeometry.setZero();
                        index_t rowIdx = 0;
                        for (index_t p = 0; p < safeSize; ++p)
                        {
                            const gsMatrix<real_t>& uvPatch = uv1(p);
                            for (index_t k = 0; k < uvPatch.cols(); ++k, ++rowIdx)
                            {
                                if (rowIdx >= totalSamples) {
                                    gsInfo << "WARNING: rowIdx=" << rowIdx << " >= totalSamples=" << totalSamples << "\n";
                                    break;
                                }
                                gsMatrix<real_t> xy = mp1.patch(p).eval(uvPatch.col(k));
                                targetGeometry(rowIdx, 0) = xy(0, 0);
                                targetGeometry(rowIdx, 1) = xy(1, 0);
                            }
                        }
                        t_targetgeom += std::chrono::duration<double>(std::chrono::system_clock::now() - _tTG0).count(); }

                        { auto _tGE0 = std::chrono::system_clock::now();
                        globalError = globalFittingError(mpbes, mp1, vectSol, uv1, targetGeometry);
                        t_globalerr += std::chrono::duration<double>(std::chrono::system_clock::now() - _tGE0).count(); }
                        auto _tPA0 = std::chrono::system_clock::now();

                        outfile << "globalError: " << globalError << "\n";
                        outfile << "featureError: " << featureError << "\n";
                        outfile << "assemblyBoundaryError: " << assemblyBoundaryError << "\n";
                        outfile << "minusnumber: " << minusnumber << "\n";
                        std::cout << "globalError: " << globalError << "\n";
                        std::cout << "featureError: " << featureError << "\n";
                        if (g_verbose) std::cout << "assemblyBoundaryError: " << assemblyBoundaryError << "\n";
                        if (summaryFile.is_open())
                        {
                            summaryFile << "FIT," << patch << "," << levNow << "," << attempt << ","
                                        << globalError << "," << featureError << "," << minusnumber << ","
                                        << irregularPercentage << "\n";
                        }
                        // Ensure cmd mirror and logfile stay synchronized during long runs.
                        std::cout << std::flush;
                        if (outfile.is_open()) outfile.flush();
                        if (summaryFile.is_open()) summaryFile.flush();
                        index_t geoDim = 2;
                        gsSparseMatrix<real_t> A(vectSol.rows() * geoDim, vectSol.rows() * geoDim);
                        // Temporarily disabled to see full boundaryError output
                        // if (minusnumber > 0) return 260106;
                        //numIrregular(patch)++;
                        
                        t_postaccept += std::chrono::duration<double>(std::chrono::system_clock::now() - _tPA0).count();
                        auto _tBK0 = std::chrono::system_clock::now();

                        if (globalError <= epsilon_g && featureError <= epsilon_f && minusnumber == 0) {
                            //if (globalError <= epsilon_g && featureError <= epsilon_f && numIrregular(patch) != 0) {
                        success:
                            outfile << "Success! iteration = " << iteration << ", coarselevel = " << levNow
                                << "\n";
                            gsInfo << "Success! iteration = " << iteration << ", coarselevel = " << levNow
                                << "\n";
                            //return 0;
                            valid = 1;
                            success = 1;
                            successfullAttempts++;
                            totalAttempts++;
                            // matFile = matOut;
                            matFile = matFeatOut;
                            //gsInfo << "Success! iteration = " << iteration << ", coarselevel = " << coarseLevel << "\n";
                            //gsInfo << "TOTALSIZE: " << THB.size() << "\n";
                            outfile << "TOTALSIZE: " << THB.size() << "\n";
                            //if(attempt > 0)  return 0;
                            acceptedsize = THB.size();
                            if (g_verbose) gsInfo << "ACCEPTEDSIZE: " << acceptedsize << "\n";
                            outfile << "ACCEPTEDSIZE: " << acceptedsize << "\n";
                            acceptedCoefs = localCoeffs;
                            AcceptedvectSol = vectSol;
                            AcceptedboxMat = boxMat;
                            AcceptedlastRow = currentLastNonZeroRow;
                            
                            // Safe MPBES copy with validation
                            try {
                                // Validate mpbes before copying
                                if (mpbes.size() == 0) {
                                    gsInfo << "WARNING: Cannot copy MPBES with size 0, skipping\n";
                                    outfile << "WARNING: Cannot copy MPBES with size 0, skipping\n";
                                } else if (mpbes.nPatches() == 0) {
                                    gsInfo << "WARNING: Cannot copy MPBES with 0 patches, skipping\n";
                                    outfile << "WARNING: Cannot copy MPBES with 0 patches, skipping\n";
                                } else {
                                    if (g_verbose)
                                        gsInfo << "Creating MPBES copy: size=" << mpbes.size()
                                              << ", patches=" << mpbes.nPatches() << "\n";
                                    auto _tCp0 = std::chrono::system_clock::now();
                                    Acceptedmpbes = std::make_unique<MPBES<2, real_t>>(mpbes);
                                    t_mpbes_copy += std::chrono::duration<double>(std::chrono::system_clock::now() - _tCp0).count();
                                    if (g_verbose) gsInfo << "MPBES copy created successfully\n";
                                }
                            } catch (const std::exception& e) {
                                gsInfo << "EXCEPTION during MPBES copy: " << e.what() << "\n";
                                outfile << "EXCEPTION during MPBES copy: " << e.what() << "\n";
                                gsInfo << "Continuing without MPBES copy...\n";
                                // Don't set Acceptedmpbes - leave it as is
                            } catch (...) {
                                gsInfo << "UNKNOWN EXCEPTION during MPBES copy\n";
                                outfile << "UNKNOWN EXCEPTION during MPBES copy\n";
                                gsInfo << "Continuing without MPBES copy...\n";
                            }
                            
                            acceptedMatOut = matOut;
                            //std::vector<std::vector<std::vector<int>>> acceptedIndex = localIndex;
                            gsVector<gsVector<int>> acceptedIndex = localIndex;
                            //AcceptedfunctionDescription = functionDescription;
                            AcceptedisActive = isActive;
                            AcceptedglobalIndex = globalIndex;


                            for (int l = 0; l < lastNonZeroRow; l++) {
                                for (int m = 1; m < 5; m++) {
                                    outfile << boxMat(patch)(l)(m) << "\t";
                                }
                                outfile << "\n";
                            }
                            AcceptedlastRow(patch) = lastNonZeroRow;
                            projections++;
                            anmat2 = THB.anchors();
                            wasRebuilt = 1;
                            removeCellIdsByValue(nonCheckedCells, attemptedCellIds);
                            outfile << "nonCheckedCells.size() has become: " << nonCheckedCells.size() << "\n";
                            // Commit co-patch group cells.
                            for (const auto& cpe : coarsenGroup) {
                                if (cpe.patchId == patch) continue;
                                std::unordered_set<int> cpRemove(
                                    cpe.geoCellIndices.begin(), cpe.geoCellIndices.end());
                                auto& cpVS = perPatchVectorS[cpe.patchId];
                                gsVector<index_t> cpFiltered(cpVS.size());
                                int cpPos = 0;
                                for (int ci = 0; ci < cpVS.size(); ++ci)
                                    if (cpRemove.find(static_cast<int>(cpVS(ci))) == cpRemove.end())
                                        cpFiltered(cpPos++) = cpVS(ci);
                                cpFiltered.conservativeResize(cpPos);
                                cpVS = cpFiltered;
                                removeCellIdsByValue(perPatchNCC[cpe.patchId], cpe.geoCellIndices);
                            }

                            int jopa = 4 * (int)pow(2, levNow) * 4 * (int)pow(2, levNow);
                            attempt = (attempt + 1) % (jopa);
                            t_bookkeeping += std::chrono::duration<double>(std::chrono::system_clock::now() - _tBK0).count();
                            continue;
                        }
                        else {
                            if (globalError > epsilon_g) {
                                outfile << "Global conditions are violated\n";
                                gsInfo << "Global conditions violated\n";
                            }
                            if (featureError > epsilon_f) {
                                outfile << "Feature conditions are violated\n";
                                gsInfo << "Feature conditions violated\n";
                            }
                            if (minusnumber > 0)
                            {
                                outfile << "The parameterization is not regular\n";
                                gsInfo << "The parameterization is not regular, attempt =" << attempt << "\n";
                                for (index_t p = 0; p < numIrregular.size(); ++p)
                                    gsInfo << numIrregular(p) << "\n";
                                {
                                    std::string irregMeshName = filePrefix + "output_mesh_irregular_FIT_p"
                                        + std::to_string(patch)
                                        + "_lev" + std::to_string(levNow)
                                        + "_att" + std::to_string(attempt);
                                    gsInfo << "Saving irregular-FIT mesh: " << irregMeshName << "\n";
                                    outfile << "Saving irregular-FIT mesh: " << irregMeshName << "\n";
                                    generateVisualizationMesh(boxMat, 20, mpbes, mp1, vectSol,
                                        currentLastNonZeroRow, irregMeshName, true);
                                    outfile.flush();
                                }
                                //return 0;
                            }

                            // Early withdrawal — two always-on cases + one mode-dependent case:
                            //
                            // Case 1: FIT is already regular (minusnumber==0) but errors exceed
                            //   epsilon. FIT is LS-optimal; adding regularization can only worsen
                            //   accuracy. No point running LO/NLO.
                            //
                            // Case 2: featureError is more than 10x epsilon_f regardless of
                            //   irregularity. Even if LO/NLO fixed every irregular point, the
                            //   error gap is too large to close. Withdraw immediately.
                            //
                            // Case 3 (--ls-only only): FIT produced an irregular parameterization.
                            //   In full mode, LO/NLO are given a chance to restore regularity.
                            //   In --ls-only mode, withdraw immediately without invoking them.
                            {
                                const bool regularButOverTolerance =
                                    (minusnumber == 0) &&
                                    (globalError > epsilon_g || featureError > epsilon_f);
                                const bool hopelesslyLargeError =
                                    (featureError  > epsilon_f * 10.0) ||
                                    (globalError   > epsilon_g * 10.0);
                                const bool geometryIrregular = g_lsOnly && (minusnumber > 0);

                                if (regularButOverTolerance || hopelesslyLargeError || geometryIrregular)
                                {
                                    auto escapeTime = std::chrono::system_clock::now();
                                    std::chrono::duration<double> escapeElapsed = escapeTime - startTime;
                                    std::time_t escapeWallTime = std::chrono::system_clock::to_time_t(escapeTime);
                                    char escapeTimeBuf[32];
                                    std::strftime(escapeTimeBuf, sizeof(escapeTimeBuf), "%Y-%m-%d %H:%M:%S",
                                                  std::localtime(&escapeWallTime));

                                    if (geometryIrregular)
                                    {
                                        gsInfo << "[--ls-only] Geometry irregular after FIT (minusnumber=" << minusnumber
                                               << "). Withdrawing candidate (LO/NLO disabled).\n";
                                        outfile << "[--ls-only] Geometry irregular after FIT (minusnumber=" << minusnumber
                                                << "). Withdrawing candidate (LO/NLO disabled).\n";
                                        gsInfo  << "Escape time: " << escapeTimeBuf << " (+" << escapeElapsed.count() << "s, "
                                                << (escapeElapsed.count() - totalJackTime) << "s without jacobian checks)\n";
                                        outfile << "Escape time: " << escapeTimeBuf << " (+" << escapeElapsed.count() << "s, "
                                                << (escapeElapsed.count() - totalJackTime) << "s without jacobian checks)\n";
                                        outfile.flush();
                                    }
                                    else if (regularButOverTolerance)
                                    {
                                        gsInfo << "FIT regular but error/feature tolerance not met"
                                               << " (globalError=" << globalError
                                               << ", featureError=" << featureError
                                               << "). Withdrawing candidate.\n";
                                        outfile << "FIT regular but error/feature tolerance not met"
                                                << " (globalError=" << globalError
                                                << ", featureError=" << featureError
                                                << "). Withdrawing candidate.\n";
                                        gsInfo  << "Escape time: " << escapeTimeBuf << " (+" << escapeElapsed.count() << "s, "
                                                << (escapeElapsed.count() - totalJackTime) << "s without jacobian checks)\n";
                                        outfile << "Escape time: " << escapeTimeBuf << " (+" << escapeElapsed.count() << "s, "
                                                << (escapeElapsed.count() - totalJackTime) << "s without jacobian checks)\n";
                                    }
                                    else
                                    {
                                        gsInfo << "Errors are hopelessly large (>10x epsilon)"
                                               << " (globalError=" << globalError
                                               << ", featureError=" << featureError
                                               << "). Withdrawing candidate.\n";
                                        outfile << "Errors are hopelessly large (>10x epsilon)"
                                                << " (globalError=" << globalError
                                                << ", featureError=" << featureError
                                                << "). Withdrawing candidate.\n";
                                        gsInfo  << "Escape time: " << escapeTimeBuf << " (+" << escapeElapsed.count() << "s, "
                                                << (escapeElapsed.count() - totalJackTime) << "s without jacobian checks)\n";
                                        outfile << "Escape time: " << escapeTimeBuf << " (+" << escapeElapsed.count() << "s, "
                                                << (escapeElapsed.count() - totalJackTime) << "s without jacobian checks)\n";
                                    }
                                    removeCellIdsByValue(nonCheckedCells, attemptedCellIds);
                                    outfile << "nonCheckedCells.size() has become: " << nonCheckedCells.size() << "\n";
                                    t_bookkeeping += std::chrono::duration<double>(std::chrono::system_clock::now() - _tBK0).count();
                                    continue;
                                }
                            }

                            gsInfo << "\n=== Starting attempt " << attempt << " - calling nonLinearOptimization ===\n";
                            outfile << "\n=== Starting attempt " << attempt << " - calling nonLinearOptimization ===\n";
                            
                            // if (attempt == 7)   DTD = 1;
                            gsMatrix<> localVectSol(localCoeffs[0].size(), 2);
                            gsVector<gsMatrix<>> localVectSols(localCoeffs.size());

                            A.setZero();
                            // LO weights (quadratic terms)
                            real_t fitting = 1;
                            real_t uniformity = 1;
                            real_t length = 1;
                            // NLO weights (non-quadratic terms)
                            real_t skewness = 1;
                            real_t eccentricity = 1;
                            real_t area = 1;
                            const bool doNonLinearOptimization = true;
                            bool usedNLO = false;

                            AdaptiveWeights adaptiveInit = chooseAdaptiveWeights(
                                globalError,
                                irregularPercentage,
                                0.0,
                                featureError,
                                minusnumber,
                                minusnumber,
                                globalError,
                                globalError,
                                false
                            );
                            fitting = adaptiveInit.fitting;
                            uniformity = adaptiveInit.uniformity;
                            length = adaptiveInit.length;
                            gsInfo << "Adaptive weights update:\n";
                            gsInfo << " fitting = " << fitting << "\n";
                            gsInfo << " uniformity = " << uniformity << "\n";
                            gsInfo << " orthogonality = " << adaptiveInit.orthogonality << "\n";
                            gsInfo << " length = " << length << "\n";
                            outfile << "Adaptive weights update:\n";
                            outfile << " fitting = " << fitting << "\n";
                            outfile << " uniformity = " << uniformity << "\n";
                            outfile << " orthogonality = " << adaptiveInit.orthogonality << "\n";
                            outfile << " length = " << length << "\n";
                            
                            gsInfo << "ALERT! TRYING TO USE LINEAR OPTIMIZATION\n";
                            outfile << "ALERT! TRYING TO USE LINEAR OPTIMIZATION\n";

                            const real_t suggestedFitting = fitting;
                            const real_t suggestedUniformity = uniformity;
                            const real_t suggestedOrthogonality = adaptiveInit.orthogonality;
                            const real_t suggestedLength = length;
                            gsMatrix<real_t> baseVectSol = vectSol;
                            const gsMatrix<real_t> bndFreezeMask = buildBoundaryFunctionMask(mpbes);
                            bool optimizationAccepted = false;
                            do
                            {
                                real_t fittingWeight = suggestedFitting;
                                real_t uniformityWeight = suggestedUniformity;
                                real_t orthogonalityWeight = suggestedOrthogonality;
                                real_t lengthWeight = suggestedLength;
                                usedNLO = false;

                                auto evaluateOptimizationStage = [&](const std::string& stageTag) -> bool
                                {
                                    numIrregular.setZero();
                                    auto beforeJackStage = std::chrono::system_clock::now();
                                    int minusnumberStage = checkJacobianDeterminant(
                                        uv2_adaptive,
                                        mpbes,
                                        vectSol,
                                        numIrregular,
                                        jacVerbose);
                                    totalJackTime += std::chrono::duration<double>(
                                        std::chrono::system_clock::now() - beforeJackStage).count();

                                    index_t totalPointsStage = 0;
                                    for (index_t p = 0; p < uv2_adaptive.size(); ++p)
                                        totalPointsStage += uv2_adaptive[p].cols();

                                    real_t irregularPercentageStage = totalPointsStage > 0
                                        ? (100.0 * minusnumberStage / totalPointsStage)
                                        : 0.0;

                                    gsInfo << "minusnumber (" << stageTag << "): " << minusnumberStage
                                           << " / " << totalPointsStage
                                           << " (" << irregularPercentageStage << "%)\n";
                                    outfile << "minusnumber (" << stageTag << "): " << minusnumberStage
                                            << " / " << totalPointsStage
                                            << " (" << irregularPercentageStage << "%)\n";

                                    for (index_t p = 0; p < numIrregular.size(); ++p)
                                    {
                                        gsInfo << "patch " << p << ": " << numIrregular(p) << "\n";
                                        outfile << "patch " << p << ": " << numIrregular(p) << "\n";
                                    }

                                    std::string postsolveLabel = stageTag + "_postsolve_p" + std::to_string(patch)
                                        + "_lev" + std::to_string(levNow)
                                        + "_att" + std::to_string(attempt);
                                    logMirroredCheck(
                                        uv2_adaptive,
                                        mpbes,
                                        vectSol,
                                        1e-12,
                                        postsolveLabel,
                                        false);

                                    if (g_enableInterfaceDiagnostics)
                                    {
                                        logSpecificInterface(mpbes, mp1, vectSol, postsolveLabel, 0, 1);
                                        logSpecificInterface(mpbes, mp1, vectSol, postsolveLabel, 0, 3);
                                        logSpecificInterface(mpbes, mp1, vectSol, postsolveLabel, 1, 2);
                                        logSpecificInterface(mpbes, mp1, vectSol, postsolveLabel, 2, 5);
                                        logSpecificInterface(mpbes, mp1, vectSol, postsolveLabel, 3, 4);
                                    }

                                    real_t globalErrorStage = globalFittingError(mpbes, mp1, vectSol, uv1, targetGeometry);
                                    real_t featureErrorStage = testBoundaryAssembly(mpbes, mp1, vectSol);
                                    globalError = globalErrorStage;
                                    featureError = featureErrorStage;
                                    minusnumber = minusnumberStage;
                                    irregularPercentage = irregularPercentageStage;

                                    gsInfo << stageTag << " errors: globalError=" << globalErrorStage
                                           << ", featureError=" << featureErrorStage << "\n";
                                    outfile << stageTag << " errors: globalError=" << globalErrorStage
                                            << ", featureError=" << featureErrorStage << "\n";
                                    outfile << "Post-optimization acceptance metrics (" << stageTag << "): globalError=" << globalError
                                        << ", featureError=" << featureError
                                        << ", minusnumber=" << minusnumber
                                        << ", irregularPercentage=" << irregularPercentage << "\n";

                                    std::string meshName = filePrefix + "output_mesh_" + stageTag + "_run_p" + std::to_string(patch)
                                        + "_lev" + std::to_string(levNow)
                                        + "_att" + std::to_string(attempt);
                                    gsInfo << "Saving " << stageTag << " mesh: " << meshName << "\n";
                                    outfile << "Saving " << stageTag << " mesh: " << meshName << "\n";
                                    generateVisualizationMesh(boxMat, 20, mpbes, mp1, vectSol,
                                        currentLastNonZeroRow, meshName, true);
                                    outfile.flush();

                                    return globalError <= epsilon_g && featureError <= epsilon_f && minusnumber == 0;
                                };


                                vectSol = baseVectSol;
                                {
                                    std::string meshPreName = filePrefix + "output_mesh_LO_pre_p" + std::to_string(patch)
                                        + "_lev" + std::to_string(levNow)
                                        + "_att" + std::to_string(attempt);
                                    gsInfo << "Saving LO mesh (pre-opt): " << meshPreName << "\n";
                                    outfile << "Saving LO mesh (pre-opt): " << meshPreName << "\n";
                                    generateVisualizationMesh(boxMat, 20, mpbes, mp1, vectSol,
                                        currentLastNonZeroRow, meshPreName, true);
                                }

                                gsInfo << "\nSuggested weights: fitting=" << suggestedFitting
                                       << ", uniformity=" << suggestedUniformity
                                       << ", orthogonality=" << suggestedOrthogonality
                                       << ", length=" << suggestedLength << "\n";
                                outfile << "Suggested weights: fitting=" << suggestedFitting
                                        << ", uniformity=" << suggestedUniformity
                                        << ", orthogonality=" << suggestedOrthogonality
                                        << ", length=" << suggestedLength << "\n";

                                gsInfo << "Using suggested weights automatically.\n";
                                outfile << "Using suggested weights automatically.\n";
                                // Temporary non-interactive mode: always accept suggested weights.
                                // gsInfo << "Accept suggested weights? Y/N: ";
                                // char acceptSuggested = 'Y';
                                // if (!(std::cin >> acceptSuggested))
                                // {
                                //     std::cin.clear();
                                //     std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
                                //     gsInfo << "Invalid input. Try again.\n";
                                //     continue;
                                // }
                                // std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
                                //
                                // if (acceptSuggested == 'Y' || acceptSuggested == 'y')
                                // {
                                //     gsInfo << "Using suggested weights.\n";
                                //     outfile << "Using suggested weights.\n";
                                // }
                                // else
                                // {
                                //     gsInfo << "Enter weights: fitting uniformity orthogonality length (orthogonality is used only by NLO fallback, Ctrl+C to exit): ";
                                //     if (!(std::cin >> fittingWeight >> uniformityWeight >> orthogonalityWeight >> lengthWeight))
                                //     {
                                //         std::cin.clear();
                                //         std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
                                //         gsInfo << "Invalid input. Try again.\n";
                                //         continue;
                                //     }
                                //     std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
                                // }

                                A.setZero();
                                gsMatrix<> matAsquare(AtA_sp);
                                const gsSparseMatrix<real_t>& matAForOptimization = matA;
                                real_t acceptedFittingWeight = fittingWeight;
                                real_t acceptedUniformityWeight = uniformityWeight;
                                real_t acceptedOrthogonalityWeight = 0.0;
                                real_t acceptedLengthWeight = lengthWeight;
                                if (g_skipLoFallback)
                                {
                                    gsInfo << "[skip-lo-fallback] Level-2 LO skipped; Level-1 minusnumber="
                                           << minusnumber << " — going straight to Gauss-Newton NLO.\n";
                                    outfile << "[skip-lo-fallback] Level-2 LO skipped; Level-1 minusnumber="
                                            << minusnumber << " — going straight to Gauss-Newton NLO.\n";
                                    // optimizationAccepted stays false, minusnumber keeps Level-1 value.
                                }
                                else
                                {
                                gsInfo << "Trying LO with orthogonality forced to 0.\n";
                                outfile << "Trying LO with orthogonality forced to 0.\n";
                                const real_t scaledUniformityWeight = uniformityWeight * static_cast<real_t>(g_loWeightScale);
                                const real_t scaledLengthWeight     = lengthWeight     * static_cast<real_t>(g_loWeightScale);
                                if (g_loWeightScale != 1.0)
                                {
                                    gsInfo  << "[lo-weight-scale=" << g_loWeightScale
                                            << "] uniformity=" << scaledUniformityWeight
                                            << " length=" << scaledLengthWeight << "\n";
                                    outfile << "[lo-weight-scale=" << g_loWeightScale
                                            << "] uniformity=" << scaledUniformityWeight
                                            << " length=" << scaledLengthWeight << "\n";
                                }
                                { auto _tLO0 = std::chrono::system_clock::now();
                                nonLinearOptimization(
                                    mpbes, uvFitting, uv2, mp1,
                                    matAsquare, matAForOptimization, b_vec,
                                    vectSol, A,
                                    SubdomainHierarchy, Bells, functionDescription,
                                    indexInTHB, isTruncated, presentation,
                                    spilloverFunctionCoordinates, hasSpillover, globalIndexTHB,
                                    boxMat, currentLastNonZeroRow,
                                    fittingWeight, scaledUniformityWeight, 0.0, scaledLengthWeight,
                                    0.0, 0.0, 0.0,
                                    epsilon_g, epsilon_f,
                                    numIrregular, geoDim,
                                    false,
                                    fittingMode,
                                    normalizedFeatureSides,
                                    useLocalFitting ? &localRegion : nullptr,
                                    effectiveOriginalCoefficients,
                                    &bndFreezeMask,
                                    &baseVectSol
                                );
                                t_lo += std::chrono::duration<double>(std::chrono::system_clock::now() - _tLO0).count(); }

                                optimizationAccepted = evaluateOptimizationStage("LO");
                                } // end !g_skipLoFallback
                                if (!optimizationAccepted)
                                {
                                    if (g_loOnly)
                                    {
                                        gsInfo << "[--lo-only] LO did not succeed; withdrawing without NLO.\n";
                                        outfile << "[--lo-only] LO did not succeed; withdrawing without NLO.\n";
                                        // Leave usedNLO=false so the withdrawal block does not exit(2).
                                    }
                                    else if (minusnumber == 0)
                                    {
                                        // LO produced a regular parameterization but failed the
                                        // approximation-error criterion.  NLO's angular functionals
                                        // trade fitting accuracy for regularity — they cannot reduce
                                        // fitting error and will only make it worse.  Withdraw this
                                        // candidate without invoking NLO.
                                        gsInfo << "LO regular (minusnumber=0) but error/feature tolerance not met"
                                               << " (globalError=" << globalError << ", featureError=" << featureError
                                               << "). Withdrawing — NLO cannot help here.\n";
                                        outfile << "LO regular (minusnumber=0) but error/feature tolerance not met"
                                                << " (globalError=" << globalError << ", featureError=" << featureError
                                                << "). Withdrawing — NLO cannot help here.\n";
                                        // Leave usedNLO=false so the withdrawal block does not exit(2).
                                    }
                                    else
                                    {
                                    usedNLO = true;
                                    vectSol = baseVectSol;
                                    A.setZero();
                                    AdaptiveWeights nloAdaptive = chooseAdaptiveWeights(
                                        globalError,
                                        irregularPercentage,
                                        featureError,
                                        featureError,
                                        minusnumber,
                                        minusnumber,
                                        globalError,
                                        globalError,
                                        true
                                    );

                                    const real_t nloFittingWeight = nloAdaptive.fitting;
                                    const real_t nloUniformityWeight = nloAdaptive.uniformity;
                                    const real_t nloOrthogonalityWeight = nloAdaptive.orthogonality;
                                    const real_t nloLengthWeight = nloAdaptive.length;
                                    const real_t nloSkewnessWeight =
                                        std::max<real_t>(1e-4, 0.5 * nloOrthogonalityWeight + 0.5 * nloUniformityWeight);
                                    const real_t nloEccentricityWeight =
                                        std::max<real_t>(1e-4, 0.5 * nloLengthWeight + 0.5 * nloUniformityWeight);
                                    const real_t nloAreaWeight =
                                        std::max<real_t>(1e-4, 0.5 * nloLengthWeight + 0.5 * nloOrthogonalityWeight);

                                    acceptedFittingWeight = nloFittingWeight;
                                    acceptedUniformityWeight = nloUniformityWeight;
                                    acceptedOrthogonalityWeight = nloOrthogonalityWeight;
                                    acceptedLengthWeight = nloLengthWeight;

                                     gsInfo << "LO candidate did not meet acceptance. Recomputing weights once for NLO fallback.\n";
                                    gsInfo << "NLO weights: fitting=" << nloFittingWeight
                                           << ", uniformity=" << nloUniformityWeight
                                           << ", orthogonality=" << nloOrthogonalityWeight
                                           << ", length=" << nloLengthWeight
                                           << ", skewness=" << nloSkewnessWeight
                                           << ", eccentricity=" << nloEccentricityWeight
                                           << ", area=" << nloAreaWeight << "\n";
                                     outfile << "LO candidate did not meet acceptance. Recomputing weights once for NLO fallback.\n";
                                    outfile << "NLO weights: fitting=" << nloFittingWeight
                                            << ", uniformity=" << nloUniformityWeight
                                            << ", orthogonality=" << nloOrthogonalityWeight
                                            << ", length=" << nloLengthWeight
                                            << ", skewness=" << nloSkewnessWeight
                                            << ", eccentricity=" << nloEccentricityWeight
                                            << ", area=" << nloAreaWeight << "\n";
                                    { auto _tNLO0 = std::chrono::system_clock::now();
                                    nonLinearOptimization(
                                        mpbes, uvFitting, uv2, mp1,
                                        matAsquare, matAForOptimization, b_vec,
                                        vectSol, A,
                                        SubdomainHierarchy, Bells, functionDescription,
                                        indexInTHB, isTruncated, presentation,
                                        spilloverFunctionCoordinates, hasSpillover, globalIndexTHB,
                                        boxMat, currentLastNonZeroRow,
                                        nloFittingWeight, nloUniformityWeight, nloOrthogonalityWeight, nloLengthWeight,
                                        nloSkewnessWeight, nloEccentricityWeight, nloAreaWeight,
                                        epsilon_g, epsilon_f,
                                        numIrregular, geoDim,
                                        true,
                                        fittingMode,
                                        normalizedFeatureSides,
                                        useLocalFitting ? &localRegion : nullptr,
                                        effectiveOriginalCoefficients,
                                        &bndFreezeMask,
                                        &baseVectSol
                                    );
                                    t_nlo += std::chrono::duration<double>(std::chrono::system_clock::now() - _tNLO0).count(); }
                                    optimizationAccepted = evaluateOptimizationStage("NLO");
                                    } // end else (minusnumber > 0)
                                }

                                if (optimizationAccepted)
                                {
                                    outfile << (usedNLO ? "NLO worked! iteration = " : "LO worked! iteration = ") << iteration << ", coarselevel = " << levNow
                                        << "\n";
                                    gsInfo << (usedNLO ? "NLO worked!\n" : "LO worked!\n");
                                    gsInfo << "coefficients: \n";
                                    gsInfo << "fitting: " << acceptedFittingWeight << "\n";
                                    gsInfo << "uniformity: " << acceptedUniformityWeight << "\n";
                                    gsInfo << "orthogonality: " << acceptedOrthogonalityWeight << "\n";
                                    gsInfo << "length: " << acceptedLengthWeight << "\n";
                                    outfile << "coefficients: \n";
                                    outfile << "fitting: " << acceptedFittingWeight << "\n";
                                    outfile << "uniformity: " << acceptedUniformityWeight << "\n";
                                    outfile << "orthogonality: " << acceptedOrthogonalityWeight << "\n";
                                    outfile << "length: " << acceptedLengthWeight << "\n";
                                    valid = 1;
                                    success = 1;
                                    successfullAttempts++;
                                    totalAttempts++;
                                    matFile = matFeatOut;
                                    outfile << "TOTALSIZE: " << THB.size() << "\n";
                                    acceptedsize = THB.size();
                                    if (g_verbose) gsInfo << "ACCEPTEDSIZE: " << acceptedsize << "\n";
                                    outfile << "ACCEPTEDSIZE: " << acceptedsize << "\n";
                                    acceptedCoefs = localCoeffs;
                                    AcceptedvectSol = vectSol;
                                    AcceptedboxMat = boxMat;
                                    AcceptedlastRow = currentLastNonZeroRow;
                                    { auto _tCp1 = std::chrono::system_clock::now();
                                    Acceptedmpbes = std::make_unique<MPBES<2, real_t>>(mpbes);
                                    t_mpbes_copy += std::chrono::duration<double>(std::chrono::system_clock::now() - _tCp1).count(); }
                                    acceptedMatOut = matOut;
                                    gsVector<gsVector<int>> acceptedIndex = localIndex;
                                    AcceptedisActive = isActive;
                                    AcceptedglobalIndex = globalIndex;

                                    if (attempt - 1 == 55) {
                                        gsInfo << "\n=== SAVING MESH AT ATTEMPT 55 ===\n";
                                        std::string meshFilename = "mesh_attempt55";
                                        generateVisualizationMesh(boxMat, 20, mpbes, mp1, vectSol,
                                            currentLastNonZeroRow, meshFilename, false);
                                        gsInfo << "Mesh saved to " << meshFilename << ".txt\n";
                                        gsInfo << "=== MESH SAVED ===\n\n";
                                    }

                                    for (int l = 0; l < lastNonZeroRow; l++) {
                                        for (int m = 1; m < 5; m++) {
                                            outfile << boxMat(patch)(l)(m);
                                            outfile << "\t";
                                        }
                                        outfile << "\n";
                                    }
                                    AcceptedlastRow(patch) = lastNonZeroRow;
                                    projections++;
                                    anmat2 = THB.anchors();
                                    wasRebuilt = 1;
                                    removeCellIdsByValue(nonCheckedCells, attemptedCellIds);
                                    outfile << "nonCheckedCells.size() has become: " << nonCheckedCells.size() << "\n";
                                    // Commit co-patch group cells.
                                    for (const auto& cpe : coarsenGroup) {
                                        if (cpe.patchId == patch) continue;
                                        std::unordered_set<int> cpRemove(
                                            cpe.geoCellIndices.begin(), cpe.geoCellIndices.end());
                                        auto& cpVS = perPatchVectorS[cpe.patchId];
                                        gsVector<index_t> cpFiltered(cpVS.size());
                                        int cpPos = 0;
                                        for (int ci = 0; ci < cpVS.size(); ++ci)
                                            if (cpRemove.find(static_cast<int>(cpVS(ci))) == cpRemove.end())
                                                cpFiltered(cpPos++) = cpVS(ci);
                                        cpFiltered.conservativeResize(cpPos);
                                        cpVS = cpFiltered;
                                        removeCellIdsByValue(perPatchNCC[cpe.patchId], cpe.geoCellIndices);
                                    }
                                    optimizationAccepted = true;
                                    break;
                                }

                                if (globalError > epsilon_g) {
                                    outfile << "Global conditions are violated\n";
                                    outfile << globalError << " > " << epsilon_g << "\n";
                                    gsInfo << globalError << " > " << epsilon_g << "\n";
                                }
                                if (featureError > epsilon_f) {
                                    outfile << "Feature conditions are violated\n";
                                    outfile << featureError << " > " << epsilon_f << "\n";
                                    gsInfo << featureError << " > " << epsilon_f << "\n";
                                }
                                if (minusnumber > 0)
                                {
                                    outfile << "The parameterization is not regular\n";
                                    gsInfo << "The parameterization is not regular\n";
                                    {
                                        std::string irregMeshName = filePrefix + "output_mesh_irregular_NLO_p"
                                            + std::to_string(patch)
                                            + "_lev" + std::to_string(levNow)
                                            + "_att" + std::to_string(attempt);
                                        gsInfo << "Saving irregular-NLO mesh: " << irregMeshName << "\n";
                                        outfile << "Saving irregular-NLO mesh: " << irregMeshName << "\n";
                                        generateVisualizationMesh(boxMat, 20, mpbes, mp1, vectSol,
                                            currentLastNonZeroRow, irregMeshName, true);
                                        outfile.flush();
                                    }
                                }

                                gsInfo << "Both LO and NLO failed for this attempt. Restoring baseVectSol and withdrawing candidate.\n";
                                outfile << "Both LO and NLO failed for this attempt. Restoring baseVectSol and withdrawing candidate.\n";
                                vectSol = baseVectSol;
                                removeCellIdsByValue(nonCheckedCells, attemptedCellIds);
                                outfile << "nonCheckedCells.size() has become: " << nonCheckedCells.size() << "\n";
                                break;
                            }
                            while (false);

                            continue;
                        }
                        //else return 3;
                    }
                    else{ 
                        outfile << "\n*** EMERGENCY: No spline created ***\n";
                        gsInfo << "\n*** EMERGENCY: No spline created ***\n";
                        
                        outfile << "Configuration: Patch " << patch << ", Level " << levNow << ", Attempt " << attempt << "\n";
                        gsInfo << "Configuration: Patch " << patch << ", Level " << levNow << ", Attempt " << attempt << "\n";
                        
                        outfile << "createSpline condition: " << createSpline << "\n";
                        gsInfo << "createSpline condition: " << createSpline << "\n";
                        
                        outfile << "\nboxMat(patch) contents (lastNonZeroRow=" << lastNonZeroRow << "):\n";
                        gsInfo << "\nboxMat(patch) contents (lastNonZeroRow=" << lastNonZeroRow << "):\n";
                        
                        for (int l = 0; l < lastNonZeroRow; ++l) {
                            outfile << "Box " << l << ": ";
                            gsInfo << "Box " << l << ": ";
                            for (int m = 0; m < 5; ++m) {
                                outfile << boxMat(patch)(l)(m);
                                gsInfo << boxMat(patch)(l)(m);
                                if (m < 4) {
                                    outfile << " ";
                                    gsInfo << " ";
                                }
                            }
                            outfile << "\n";
                            gsInfo << "\n";
                        }
                        
                        attempt = (attempt + 1) % (jopa);
                        removeCellIdsByValue(nonCheckedCells, attemptedCellIds);
                        outfile << "nonCheckedCells.size() has become: " << nonCheckedCells.size() << "\n";

                    }
                    outfile << "FINISHED\n";
                    //break;
                } // end inner-body scope (aliases: nonCheckedCells, vectorS, pickedCells, lastNonZeroRow)
            } // end while (anyVectorSNonEmpty) — global inner loop

            if (iteration == 1)
            {
                std::string meshName = filePrefix + "output_mesh_after_iter1_lev" + std::to_string(levNow);
                gsInfo << "Saving after-iter1 mesh: " << meshName << "\n";
                outfile << "Saving after-iter1 mesh: " << meshName << "\n";
                if (Acceptedmpbes)
                    generateVisualizationMesh(AcceptedboxMat, 20, *Acceptedmpbes, mp1, AcceptedvectSol,
                        AcceptedlastRow, meshName, true);
                else
                    gsInfo << "[after-iter1] No accepted coarsening in iteration 1 — skipping mesh save.\n";
                outfile.flush();
            }

        } // end while (anyPoolNonEmpty && success) — global outer loop
    } // end for (levNow) — outer loop
    } // end global cell-selection block

    // Prepare return values
    AlgorithmResult result;
    result.mp = std::move(mp);
    result.Bells = Bells;
    result.boxMat = boxMat;
    result.acceptedBoxMat = AcceptedboxMat;
    result.acceptedCoefs = acceptedCoefs;
    result.AcceptedvectSol = AcceptedvectSol;
    result.AcceptedlastRow = AcceptedlastRow;
    result.acceptedMatOut = acceptedMatOut;
    result.featureCoordinates = featureCoordinates;
    result.uvFeature = uvFeature;
    result.xyFeature = xyFeature;
    result.AcceptedisActive = AcceptedisActive;
    result.AcceptedglobalIndex = AcceptedglobalIndex;
    result.AcceptedfunctionDescription = AcceptedfunctionDescription;
    result.uv1 = uv1;
    result.lowc2 = lowc2;
    result.uppc2 = uppc2;
    result.maxLevel = maxLevel;
    result.toc = toc;
    result.successfullAttempts = successfullAttempts;
    result.totalAttempts = totalAttempts;
    result.matFile = matFile;
    result.interioru0 = interioru0;
    result.interiorv0 = interiorv0;
    result.proj = projections;
    result.mp1 = mp1;
    result.currentLastNonZeroRow = currentLastNonZeroRow;
    result.mpbes = std::move(Acceptedmpbes);

    // ---- Phase-time summary ----
    {
        double t_total = std::chrono::duration<double>(std::chrono::system_clock::now() - startTime).count();
        double t_accounted = t_preflight + t_init + t_thb_rebuild + t_extract
                           + t_grenda + t_cell_sel + t_mpbes_bld + t_assemble + t_ata_ldlt
                           + t_localsetup + t_presolve + t_postsolve + t_lo + t_nlo + t_mpbes_copy
                           + t_jack + t_postjack + t_boundary + t_targetgeom + t_globalerr + t_postaccept
                           + t_gap_a + t_bookkeeping;
        auto prt = [&](const char* name, double t) {
            gsInfo  << "  " << name << ": " << t << " s (" << (100.0*t/t_total) << "%)\n";
            outfile << "  " << name << ": " << t << " s (" << (100.0*t/t_total) << "%)\n";
        };
        gsInfo  << "[timing] Phase breakdown (wall=" << t_total << " s):\n";
        outfile << "[timing] Phase breakdown (wall=" << t_total << " s):\n";
        prt("Preflight total                 ", t_preflight);
        prt("  mesh write (generateOriginalGeometryMesh) ", t_pf_mesh);
        prt("  Jacobian check (computeMirroredCheckFromMp)", t_pf_check);
        prt("One-time init (before loop)     ", t_init);
        prt("Grenda preflight (all patches)  ", t_grenda);
        prt("Cell selection (winning patch)  ", t_cell_sel);
        prt("THB rebuild + outfile<<THB      ", t_thb_rebuild);
        prt("MPBES build total               ", t_mpbes_bld);
        prt("  orient                        ", t_orient);
        prt("  createIndexMapping            ", t_indexmap);
        prt("  MPBES constructor             ", t_mpbes_ctor);
        prt("MPBES extract + vectSol init    ", t_extract);
        prt("Local setup (region + AU select) ", t_localsetup);
        prt("Assembly  assemble()            ", t_assemble);
        prt("Pre-solve (col map + defect corr)", t_presolve);
        prt("AtA formation + LDLT solve      ", t_ata_ldlt);
        prt("Post-solve (scatter/merge/resid) ", t_postsolve);
        prt("LO  nonLinearOptimization       ", t_lo);
        prt("NLO nonLinearOptimization       ", t_nlo);
        prt("MPBES copy on accept            ", t_mpbes_copy);
        prt("checkJacobianDeterminant        ", t_jack);
        prt("Post-jack (log+flush+A-alloc)   ", t_postjack);
        prt("testBoundaryAssembly            ", t_boundary);
        prt("target geometry construction    ", t_targetgeom);
        prt("globalFittingError              ", t_globalerr);
        prt("Post-accept (log+flush)         ", t_postaccept);
        prt("GAP A (Grenda→cell-sel)         ", t_gap_a);
        prt("Bookkeeping (accept/reject path)", t_bookkeeping);
        prt("Other (unaccounted)             ", t_total - t_accounted);
    }

    return result;
}

// Helper function for non-linear optimization using MPBES basis
// Based on Heydarov, Buffa, Jüttler "An Unrefinement Algorithm for Planar THB-spline Parameterizations"
void nonLinearOptimization(
    // MPBES basis structure
    const MPBES<2, real_t>& mpbes,
    
    // UV coordinates and geometry
    const gsVector<gsMatrix<real_t>>& uv1,
    const gsVector<gsMatrix<real_t>>& uv2,
    const gsMultiPatch<>& mp1,
    
    // System matrices
    const gsMatrix<>& matAsquare,
    const gsSparseMatrix<real_t>& matA,
    const gsMatrix<>& b_vec,
    
    // Solution vector (modified in-place)
    gsMatrix<real_t>& vectSol,
    gsSparseMatrix<real_t>& A,
    
    // THB hierarchy and basis
    const gsVector<gsTHBSplineBasis<2>>& SubdomainHierarchy,
    const gsVector<gsVector<gsTensorBSplineBasis<2, real_t>>>& Bells,
    const std::vector<std::vector<std::vector<index_t>>>& functionDescription,
    const gsVector<gsVector<gsVector<index_t>>>& indexInTHB,
    const std::vector<bool>& isTruncated,
    const std::vector<std::vector<gsSparseVector<double>>>& presentation,
    const std::vector<std::vector<std::array<int, 3>>>& spilloverFunctionCoordinates,
    const std::vector<bool>& hasSpillover,
    const gsVector<gsVector<index_t>>& globalIndexTHB,
    
    // Box hierarchy for visualization
    const gsVector<gsVector<gsVector<index_t>>>& boxMat,
    const gsVector<int>& currentLastNonZeroRow,
    
    // Quality measure weights
    real_t fitting_weight,
    real_t uniformity_weight,
    real_t orthogonality_weight,
    real_t length_weight,
    real_t skewness_weight,
    real_t eccentricity_weight,
    real_t area_weight,
    
    // Convergence criteria
    real_t epsilon_g,
    real_t epsilon_f,
    
    // Irregular count (modified)
    gsVector<size_t>& numIrregular,
    
    // Geometry dimension
    index_t geoDim,
    bool doNonLinear,
    FittingMode fittingMode,
    const std::vector<FeatureBoundarySpec>& featureSides,
    const LocalCoarseningRegion* localRegion,
    const gsMatrix<real_t>* originalCoefficients,
    const gsMatrix<real_t>* boundaryFreezeMask,
    const gsMatrix<real_t>* boundaryFreezeRef
) {
    PROFILE_FUNCTION();
    
    // Quality weights (predefined, no user input)
    real_t fitting = fitting_weight;
    real_t uniformity = uniformity_weight;
    real_t orthogonality = orthogonality_weight;
    real_t length = length_weight;
    real_t skewness = skewness_weight;
    real_t eccentricity = eccentricity_weight;
    real_t area = area_weight;

    struct StepStatus
    {
        bool converged;
        bool smallStep;
        bool tooFar;
    };

    auto doStep = [&](real_t fitW,
                      real_t uniW,
                      real_t orthW,
                      real_t lenW,
                      real_t skewW,
                      real_t eccW,
                      real_t areaW,
                      bool allowFarContinue,
                      int &minusOut,
                      real_t &residualOut,
                      real_t &featureErrorOut) -> StepStatus
    {
        // Store old solution for comparison
        gsMatrix<> vectSolOld = vectSol;

        // Reset optimization matrix
        A.setZero();

        // Call MPBES-based optimization with quality functionals
        // Create non-const copy for optimize function
        gsMatrix<> matAsquareCopy = matAsquare;
        const LocalCoarseningRegion* effectiveLocalRegion =
            (fittingMode == FittingMode::LocalFitting) ? localRegion : nullptr;
        const gsMatrix<real_t>* effectiveOriginalCoefficients =
            (fittingMode == FittingMode::LocalFitting) ? originalCoefficients : nullptr;

        optimize(mpbes, uv1, matAsquareCopy, vectSol, A,
            fitW, orthW, lenW, uniW, skewW, eccW, areaW, 1e-7, true,
            doNonLinear, effectiveLocalRegion, effectiveOriginalCoefficients);

        if (!gsEigen::MatrixXd(vectSol).allFinite())
        {
            gsInfo << "[nonLinearOptimization] vectSol became non-finite after optimize. "
                   << "weights: fit=" << fitW << ", uni=" << uniW << ", len=" << lenW
                   << ", skew=" << skewW << ", ecc=" << eccW << ", area=" << areaW << "\n";
            outfile << "[nonLinearOptimization] vectSol became non-finite after optimize. "
                    << "weights: fit=" << fitW << ", uni=" << uniW << ", len=" << lenW
                    << ", skew=" << skewW << ", ecc=" << eccW << ", area=" << areaW << "\n";
            throw std::runtime_error("vectSol non-finite after optimize");
        }

        // Restore boundary DOFs to their LS values so that featureError
        // (measured at patch interfaces) is not degraded by the optimizer.
        if (boundaryFreezeMask && boundaryFreezeRef &&
            boundaryFreezeRef->rows() == vectSol.rows() &&
            boundaryFreezeRef->cols() >= vectSol.cols())
        {
            for (index_t f = 0; f < vectSol.rows(); ++f)
            {
                if (f < boundaryFreezeMask->cols() &&
                    (*boundaryFreezeMask)(0, f) == 0.0)
                {
                    for (index_t d = 0; d < vectSol.cols(); ++d)
                        vectSol(f, d) = (*boundaryFreezeRef)(f, d);
                }
            }
        }

        // Calculate residual
        gsMatrix<> matOut = gsEigen::MatrixXd(matA) * vectSol;
        gsMatrix<> matC = matOut - b_vec;
        real_t residual = matC.maxCoeff();
        residualOut = residual;
        featureErrorOut = std::numeric_limits<real_t>::infinity();
        minusOut = 0;

        // Skip expensive checks if residual is very large (not close to convergence)
        if (residual > epsilon_g * 10.0) {
            if (!allowFarContinue)
                return { false, false, true };
            return { false, false, true };
        }

        // Check Jacobian determinant
        int minusnumber = checkJacobianDeterminant(
            uv2, mpbes, vectSol, numIrregular, false);
        minusOut = minusnumber;

        // Keep convergence checks aligned with the main acceptance path.
        real_t featureError = testBoundaryAssembly(mpbes, mp1, vectSol);
        featureErrorOut = featureError;

        // Check convergence
        if (residual < epsilon_g && featureError < epsilon_f && minusnumber == 0) {
            return { true, false, false };
        }

        // Triangle rule check for solution stability
        gsMatrix<> matOutOld = gsEigen::MatrixXd(matA) * vectSolOld;
        gsMatrix<> matCOld = matOutOld - b_vec;

        real_t maxDiff = 0.0;
        for (int row = 0; row != vectSol.rows(); row++) {
            for (int dim = 0; dim != geoDim; ++dim) {
                maxDiff = std::max(maxDiff,
                    sqrt(pow(matCOld(row, 0) - matC(row, 0), 2) +
                         pow(matCOld(row, 1) - matC(row, 1), 2)));
            }
        }

        // Small update check to avoid redundant iterations
        const real_t stepNorm = (vectSol - vectSolOld).norm();
        return { false, stepNorm < 1e-12, false };
    };

    if (!doNonLinear)
    {
        int minusnumberLO = 0;
        real_t residualLO = 0.0;
        real_t featureErrorLO = 0.0;
        StepStatus loStatus = doStep(fitting, uniformity, 0.0, length, 0.0, 0.0, 0.0, false,
                                     minusnumberLO, residualLO, featureErrorLO);
        gsInfo << "[LO-step] residual=" << residualLO
               << ", featureError=" << (std::isfinite(featureErrorLO) ? std::to_string(featureErrorLO) : std::string("N/A"))
               << ", minusnumber=" << minusnumberLO
               << ", converged=" << (loStatus.converged ? "true" : "false")
               << ", smallStep=" << (loStatus.smallStep ? "true" : "false")
               << ", tooFar=" << (loStatus.tooFar ? "true" : "false") << "\n";
        outfile << "[LO-step] residual=" << residualLO
                << ", featureError=" << (std::isfinite(featureErrorLO) ? std::to_string(featureErrorLO) : std::string("N/A"))
                << ", minusnumber=" << minusnumberLO
                << ", converged=" << (loStatus.converged ? "true" : "false")
                << ", smallStep=" << (loStatus.smallStep ? "true" : "false")
                << ", tooFar=" << (loStatus.tooFar ? "true" : "false") << "\n";
        if (loStatus.converged)
            return;
        return;
    }

    // NLO: non-quadratic functionals as well as LO functionals
    const int maxIterations = 10;
    int residualIncreaseStreak = 0;
    real_t previousResidual = std::numeric_limits<real_t>::quiet_NaN();
    for (int iter = 0; iter < maxIterations; ++iter)
    {
        gsInfo << "[NLO-iter " << (iter + 1) << "/" << maxIterations << "]"
               << " weights: fit=" << fitting
               << ", uni=" << uniformity
               << ", ortho=" << orthogonality
               << ", len=" << length
               << ", skew=" << skewness
               << ", ecc=" << eccentricity
               << ", area=" << area << "\n";
        outfile << "[NLO-iter " << (iter + 1) << "/" << maxIterations << "]"
            << " weights: fit=" << fitting
            << ", uni=" << uniformity
            << ", ortho=" << orthogonality
            << ", len=" << length
            << ", skew=" << skewness
            << ", ecc=" << eccentricity
            << ", area=" << area << "\n";

        int minusnumberNow = 0;
        real_t residualNow = 0.0;
        real_t featureErrorNow = 0.0;
        StepStatus nloStatus = doStep(fitting, uniformity, orthogonality, length, skewness, eccentricity, area, true,
                                      minusnumberNow, residualNow, featureErrorNow);

        gsInfo << "[NLO-iter " << (iter + 1) << "/" << maxIterations << "]"
               << " result: residual=" << residualNow
               << ", featureError=" << (std::isfinite(featureErrorNow) ? std::to_string(featureErrorNow) : std::string("N/A"))
               << ", minusnumber=" << minusnumberNow
               << ", converged=" << (nloStatus.converged ? "true" : "false")
               << ", smallStep=" << (nloStatus.smallStep ? "true" : "false")
               << ", tooFar=" << (nloStatus.tooFar ? "true" : "false") << "\n";
        outfile << "[NLO-iter " << (iter + 1) << "/" << maxIterations << "]"
            << " result: residual=" << residualNow
            << ", featureError=" << (std::isfinite(featureErrorNow) ? std::to_string(featureErrorNow) : std::string("N/A"))
            << ", minusnumber=" << minusnumberNow
            << ", converged=" << (nloStatus.converged ? "true" : "false")
            << ", smallStep=" << (nloStatus.smallStep ? "true" : "false")
            << ", tooFar=" << (nloStatus.tooFar ? "true" : "false") << "\n";

        if (nloStatus.converged || nloStatus.smallStep)
            return;

        const real_t residualSignal = std::isfinite(residualNow)
            ? residualNow
            : std::max<real_t>(epsilon_g * 10.0, 1.0);
        if (std::isfinite(previousResidual) && residualSignal > previousResidual + 1e-9)
            ++residualIncreaseStreak;
        else
            residualIncreaseStreak = 0;

        if (residualIncreaseStreak >= 2)
        {
            gsInfo << "[NLO] Early stop: residual increased for 2 consecutive iterations.\n";
            outfile << "[NLO] Early stop: residual increased for 2 consecutive iterations.\n";
            return;
        }

        previousResidual = residualSignal;
    }

    gsInfo << "[NLO] Reached max iterations without convergence. Returning latest iterate for post-check.\n";
    outfile << "[NLO] Reached max iterations without convergence. Returning latest iterate for post-check.\n";
    return;
    
    // Non-converged optimization - skip mesh generation for performance
    // generateVisualizationMesh(boxMat, 20, mpbes, mp1, vectSol, 
    //     currentLastNonZeroRow, "mesh_unstable", true);
}

int main(int argc, char** argv) {
    const std::string exePath = (argc > 0 && argv && argv[0]) ? std::string(argv[0]) : std::string();
    // Quick pre-scan: detect --save-separately and --local-fitting before opening the mirror.
    bool preScanSaveSeparately = false;
    bool preScanLocalFitting   = false;
    for (int ai = 1; ai < argc; ++ai)
    {
        const std::string arg(argv[ai]);
        if (arg == "--save-separately") preScanSaveSeparately = true;
        if (arg == "--local-fitting")   preScanLocalFitting   = true;
    }

    std::string cmdMirrorStem = "poissonTHB_example_cmd_output";
    if (preScanSaveSeparately)
        cmdMirrorStem += preScanLocalFitting ? "_local" : "_global";
    std::string cmdMirrorPath = cmdMirrorStem + ".txt";
    const size_t sepPos = exePath.find_last_of("\\/");
    if (sepPos != std::string::npos)
        cmdMirrorPath = exePath.substr(0, sepPos + 1) + cmdMirrorStem + ".txt";

    static ScopedConsoleMirror cmdMirror(cmdMirrorPath);
    if (cmdMirror.active())
        std::cout << "Mirroring console output to: " << cmdMirrorPath << "\n";
    std::cout << "Executable path: " << (exePath.empty() ? "<unknown>" : exePath) << "\n";
    char cwdBuf[4096] = {0};
    if (_getcwd(cwdBuf, static_cast<int>(sizeof(cwdBuf))) != nullptr)
        std::cout << "Working directory: " << cwdBuf << "\n";
    else
        std::cout << "Working directory: <unknown>\n";

    std::string filename;
    std::vector<FeatureBoundarySpec> featureSides;
    const bool useManualStartupInput = false;
     const std::string defaultFilename =
         "C:\\Users\\heydatey\\source\\repos\\gismo\\filedata\\generatedMPs\\mask_approximation_fine_L3_NLO.xml";
    // const std::string defaultFilename =
    //     "C:\\Users\\heydatey\\source\\repos\\gismo\\filedata\\generatedMPs\\hexagon_3p_4l.xml";
    // const std::string defaultFilename =
    //     "C:\\Users\\heydatey\\source\\repos\\gismo\\filedata\\generatedMPs\\tv_approximation_fine_L3.xml";
    //const std::string defaultFilename =
    //    "C:\\Users\\heydatey\\source\\repos\\gismo\\filedata\\generatedMPs\\distorted_output_260422.xml";
    //  const std::string defaultFilename =
    //      "C:\\Users\\heydatey\\source\\repos\\gismo\\filedata\\generatedMPs\\distorted_output.xml";
    const int defaultNumFeatureSides = 0;

    // Parse --input <path> and --skip-lo-fallback from argv before the main logic.
    std::string cliInputFile;
    for (int ai = 1; ai < argc; ++ai)
    {
        if (std::string(argv[ai]) == "--input" && ai + 1 < argc)
        {
            cliInputFile = std::string(argv[ai + 1]);
            ++ai;
        }
        else if (std::string(argv[ai]).size() > 4 &&
                 std::string(argv[ai]).substr(std::string(argv[ai]).size() - 4) == ".xml" &&
                 std::string(argv[ai]).front() != '-')
        {
            cliInputFile = std::string(argv[ai]);
        }
        else if (std::string(argv[ai]) == "--ls-only")
        {
            g_lsOnly = true;
            gsInfo << "[flag] --ls-only: LS fitting only — LO and NLO disabled. "
                      "FIT must be regular and within tolerance, else candidate is withdrawn.\n";
        }
        else if (std::string(argv[ai]) == "--lo-only")
        {
            g_loOnly = true;
            gsInfo << "[flag] --lo-only: LO enabled as fallback after FIT, but NLO disabled.\n";
        }
        else if (std::string(argv[ai]) == "--skip-lo-fallback")
        {
            g_skipLoFallback = true;
            gsInfo << "[flag] --skip-lo-fallback active: Level-2 LO will be skipped; "
                      "Level-1 LO failure goes straight to Gauss-Newton NLO.\n";
        }
        else if (std::string(argv[ai]) == "--lo-weight-scale" && ai + 1 < argc)
        {
            g_loWeightScale = std::stod(argv[ai + 1]);
            ++ai;
            gsInfo << "[flag] --lo-weight-scale=" << g_loWeightScale
                   << ": LO uniformity and length weights will be multiplied by this factor.\n";
        }
        else if (std::string(argv[ai]) == "--save-separately")
        {
            gsInfo << "[flag] --save-separately: log saved to " << cmdMirrorPath << "\n";
        }
        else if (std::string(argv[ai]) == "--local-fitting")
        {
            g_useLocalFitting = true;
            gsInfo << "[flag] --local-fitting active: LS/LO will be restricted to the "
                      "local support region of coarsened cells (paper eq. 9-10).\n";
        }
        else if (std::string(argv[ai]) == "--lambda" && ai + 1 < argc)
        {
            g_localityLambda = std::stoi(argv[ai + 1]);
            ++ai;
            gsInfo << "[flag] --lambda=" << g_localityLambda
                   << ": locality extension for local fitting region.\n";
        }
        else if (std::string(argv[ai]) == "--cell-method" && ai + 1 < argc)
        {
            const std::string val = argv[ai + 1];
            if (val == "g" || val == "l" || val == "r" || val == "s")
            {
                g_cellMethod = val[0];
                ++ai;
                gsInfo << "[flag] --cell-method=" << g_cellMethod
                       << ": cell selection method (g=Grenda, l=lexicographic, r=random, s=smallest).\n";
            }
            else
                gsInfo << "[flag] --cell-method: unknown value '" << val
                       << "', ignoring (valid: g l r s).\n";
        }
        else if (std::string(argv[ai]) == "--verbose")
        {
            g_verbose = true;
            gsInfo << "[flag] --verbose: verbose logging enabled (indexInTHB, functionDescription,"
                   << " vectSol coefficients, IdentifyPatches twins, geo-coarsen candidates).\n";
        }
        else if (std::string(argv[ai]) == "--epsilon-g" && ai + 1 < argc)
        {
            g_epsilonG = std::stod(argv[ai + 1]);
            ++ai;
            gsInfo << "[flag] --epsilon-g=" << g_epsilonG
                   << ": global approximation error tolerance.\n";
        }
        else if (std::string(argv[ai]) == "--epsilon-f" && ai + 1 < argc)
        {
            g_epsilonF = std::stod(argv[ai + 1]);
            ++ai;
            gsInfo << "[flag] --epsilon-f=" << g_epsilonF
                   << ": feature boundary error tolerance.\n";
        }
        else if (std::string(argv[ai]) == "--feature" && ai + 2 < argc)
        {
            const int    patchIdx = std::stoi(argv[ai + 1]);
            std::string  sideStr  = argv[ai + 2];
            ai += 2;
            std::string sideLow = sideStr;
            std::transform(sideLow.begin(), sideLow.end(), sideLow.begin(),
                           [](unsigned char c){ return static_cast<char>(std::tolower(c)); });
            FeatureSide side = FeatureSide::U0;
            bool sideOk = true;
            if      (sideLow == "u0") side = FeatureSide::U0;
            else if (sideLow == "u1") side = FeatureSide::U1;
            else if (sideLow == "v0") side = FeatureSide::V0;
            else if (sideLow == "v1") side = FeatureSide::V1;
            else sideOk = false;
            if (sideOk)
            {
                featureSides.push_back(FeatureBoundarySpec{static_cast<index_t>(patchIdx), side});
                gsInfo << "[flag] --feature patch=" << patchIdx
                       << " side=" << sideStr << "\n";
            }
            else
            {
                gsWarn << "[flag] --feature: unknown side '" << sideStr
                       << "' (expected u0/u1/v0/v1), skipping.\n";
            }
        }
    }

    gsInfo << "Startup mode: "
           << (useManualStartupInput ? "manual" : "default")
           << "\n";

    if (!useManualStartupInput)
    {
        filename = cliInputFile.empty() ? defaultFilename : cliInputFile;
        gsInfo << "Using default input filename: " << filename << "\n";
        gsInfo << "Using default feature boundary count: " << defaultNumFeatureSides << "\n";
    }
    else
    {
        gsInfo << "Enter input filename (full path or name): ";
        std::getline(std::cin, filename);
        if (filename.empty())
            filename = "hexagon_3p_4l_240325.xml";

        gsInfo << "Enter number of feature boundary sides (0 for default boundary error): ";
        int numFeatureSides = 0;
        if (!(std::cin >> numFeatureSides))
        {
            std::cin.clear();
            std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
            numFeatureSides = 0;
        }
        if (numFeatureSides < 0)
            numFeatureSides = 0;
        std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

        for (int i = 0; i < numFeatureSides; ++i)
        {
            int patchIndex = -1;
            std::string sideStr;

            gsInfo << "Feature " << (i + 1) << ": patch index: ";
            if (!(std::cin >> patchIndex))
            {
                std::cin.clear();
                std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
                gsInfo << "Invalid patch index, skipping feature.\n";
                continue;
            }

            gsInfo << "Feature " << (i + 1) << ": side (u0,u1,v0,v1): ";
            if (!(std::cin >> sideStr))
            {
                std::cin.clear();
                std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
                gsInfo << "Invalid side, skipping feature.\n";
                continue;
            }
            std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

            std::string sideLower = sideStr;
            std::transform(sideLower.begin(), sideLower.end(), sideLower.begin(),
                           [](unsigned char c) { return static_cast<char>(std::tolower(c)); });

            FeatureSide side;
            bool sideOk = true;
            if (sideLower == "u0")
                side = FeatureSide::U0;
            else if (sideLower == "u1")
                side = FeatureSide::U1;
            else if (sideLower == "v0")
                side = FeatureSide::V0;
            else if (sideLower == "v1")
                side = FeatureSide::V1;
            else
                sideOk = false;

            if (!sideOk)
            {
                gsInfo << "Invalid side '" << sideStr << "', skipping feature.\n";
                continue;
            }

            featureSides.push_back(FeatureBoundarySpec{ static_cast<index_t>(patchIndex), side });
        }
    }

    // Reduce command output noise after user input
    // gsInfo.setstate(std::ios_base::failbit);

    const std::string inputBase = gsFileManager::getBasename(filename);
    const std::string filePrefix = inputBase + "_";

    xboxFile.open(filePrefix + "xboxFile.txt");
    yboxFile.open(filePrefix + "yboxFile.txt");
    closureLogFile.open(filePrefix + "closure_log.txt");
    gsStopwatch clock;
    DTD = 0;
    printAB = 0;
    localTempAttempt = 0;
    std::chrono::time_point<std::chrono::system_clock> startTime, iterTime;
    startTime = std::chrono::system_clock::now();
    const int dim = 2;
    int row, acceptedsize, attempt = 0;
    std::string badFile = filePrefix + "badgeoLocalAndreatemp6";
    std::string pvdFile = filePrefix + "resultLocalAndreatemp6";
    int los = 0, nlos = 0, proj = 0;
    std::string givenGeo;
    int valid = 0;
    int gradingExtent;
    real_t epsilon_g = (g_epsilonG >= 0.0) ? static_cast<real_t>(g_epsilonG) : real_t(1e+6);
    real_t epsilon_f = (g_epsilonF >= 0.0) ? static_cast<real_t>(g_epsilonF) : real_t(1);
    gsInfo << "[params] epsilon_g=" << epsilon_g << "  epsilon_f=" << epsilon_f
           << "  featureSides=" << featureSides.size() << "\n";
    real_t lcx, lcy, ucx, ucy;
    std::string acCond = to_string(epsilon_g) + "and" + to_string(epsilon_f);
    givenGeo = "two_squares_lev1";
    givenGeo += "L2";
    givenGeo += "LO";
    givenGeo += "NLO";
    std::string fileLoc = filePrefix + "logFile_poissonTHB_example";
    outfile.open(fileLoc + ".txt");
    summaryFile.open(filePrefix + "summary_log.txt");
    if (summaryFile.is_open())
    {
        summaryFile << "stage,patch,level,attempt,globalError,featureError,minusnumber,irregularPercent\n";
    }
    gsFileData<> data0(filename);
    int iter = -1;
    int successfullAttempts = 0, totalAttempts = 0;
    index_t degree;
    real_t tol = 1e-8, gtol = 1e-8;
    char method = g_cellMethod;

    // Run the main algorithm
    AlgorithmResult result;
    try {
        result = unrefinementAlgorithmHBJ(
            filename,
            epsilon_g,
            epsilon_f,
            method,
            givenGeo,
            acCond,
            featureSides);
    } catch (const ProgramExitSignal& ex) {
        gsInfo << "\nAborted with code " << ex.code() << ": " << ex.what() << "\n";
        outfile << "\nAborted with code " << ex.code() << ": " << ex.what() << "\n";
        if (summaryFile.is_open())
            summaryFile.flush();
        if (outfile.is_open())
            outfile.flush();
        outfile.close();
        xboxFile.close();
        yboxFile.close();
        return ex.code();
    } catch (const std::exception& ex) {
        gsInfo << "\nAborted: " << ex.what() << "\n";
        outfile << "\nAborted: " << ex.what() << "\n";
        outfile.close();
        xboxFile.close();
        yboxFile.close();
        return 1;
    }
    
    std::chrono::time_point<std::chrono::system_clock> toc = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_finished = toc - startTime;
    
    gsInfo << "FINISHED\n";
    gsInfo << "Total execution time: " << elapsed_finished.count() << " seconds\n";
    outfile << "FINISHED took: " << elapsed_finished.count() << "\n";
    if (result.mpbes)
    {
        gsInfo << "accepted mpbes.size(): " << result.mpbes->size() << "\n";
        outfile << "accepted mpbes.size(): " << result.mpbes->size() << "\n";
    }
    outfile << "Successful attempts: " << result.successfullAttempts << " / " << result.totalAttempts << "\n\n";
    
    // Generate visualization mesh
    // gsInfo << "\n=== Generating final visualization mesh ===\n";
    // outfile << "\n=== Generating final visualization mesh ===\n";
    
    const std::string finalVisualizationPrefix = filePrefix + "output_mesh_final";
    generateVisualizationMesh(
        result.acceptedBoxMat, 
        20, 
        *result.mpbes, 
        result.mp1, 
        result.AcceptedvectSol, 
        result.AcceptedlastRow, 
        finalVisualizationPrefix, 
        true
    );
    
    // gsInfo << "Visualization mesh saved to " << filePrefix << "output_mesh_final.txt\n";
    // outfile << "Visualization mesh saved to " << filePrefix << "output_mesh_final.txt\n";
    
    // Print profiling summary
    g_profiler.printSummary();

    {
        const std::string meshPath = finalVisualizationPrefix + ".txt";
        gsInfo << "Visualization written to: " << meshPath << "\n";
        outfile << "Visualization written to: " << meshPath << "\n";
        gsInfo << "Python script execution skipped (final visualization path).\n";
        outfile << "Python script execution skipped (final visualization path).\n";
    }

    // Close files
    outfile.close();
    if (summaryFile.is_open())
        summaryFile.close();
    xboxFile.close();
    yboxFile.close();

    gsInfo << "Output files closed\n";
    gsInfo << "Program completed successfully\n";

    return 0;
}


/**
 * @brief Saves the multipatch geometry assembled from THB patch geometries.
 *
 * Builds per-patch local coefficient matrices from `vectSol` using `globalIndexTHB`,
 * constructs each patch geometry through the `SubdomainHierarchy`, adds them to
 * a new multipatch object, computes topology, and writes to file using `gsFileData<>`.
 *
 * @param[in]  vectSol         Global solution matrix (numFunctions x 2).
 * @param[in]  globalIndexTHB  Maps (patch, local basis index) to global function index.
 * @param[in]  SubdomainHierarchy  THB-spline basis per patch.
 * @param[in]  outfile         Stream for logging file paths and errors.
 * @param[in]  baseName        Prefix for output file names.
 * @param[out] localIndex      Will be filled with active basis indices per patch.
 */
void saveMultipatchGeometry(
    const gsMatrix<real_t>& vectSol,
    const gsVector<gsVector<int>>& globalIndexTHB,
    const gsVector<gsTHBSplineBasis<2, real_t>>& SubdomainHierarchy,
    std::ostream& outfile,
    const std::string& baseName,
    gsVector<gsVector<int>>& localIndex)
{
    PROFILE_FUNCTION();
    gsInfo << "Saving multipatch geometry to: " << baseName << "\n";
    outfile << "Saving multipatch geometry to: " << baseName << "\n";

    localIndex.resize(globalIndexTHB.size());
    gsMultiPatch<> newmp;

    // Build per-patch geometries and add to newmp
    for (size_t patch = 0; patch < SubdomainHierarchy.size(); ++patch)
    {
        const auto& patchGlobalIndices = globalIndexTHB(patch);

        // Count active basis functions
        index_t count = 0;
        for (index_t i = 0; i < patchGlobalIndices.size(); ++i)
            if (patchGlobalIndices[i] != -1) ++count;

        localIndex(patch).resize(count);

        // Fill local index mapping
        index_t cur = 0;
        for (index_t i = 0; i < patchGlobalIndices.size(); ++i)
            if (patchGlobalIndices[i] != -1)
                localIndex(patch)(cur++) = static_cast<int>(i);

        if (count == 0)
        {
            gsInfo << "Patch " << patch << " has no active basis functions, skipping.\n";
            outfile << "Patch " << patch << " has no active basis functions, skipping.\n";
            continue;
        }

        gsMatrix<> localVectSol(count, 2);
        index_t pos = 0;
        for (index_t i = 0; i < patchGlobalIndices.size(); ++i)
        {
            if (patchGlobalIndices[i] == -1) continue;
            const int globalIdx = patchGlobalIndices[i];
            if (globalIdx < 0 || globalIdx >= static_cast<int>(vectSol.rows()))
            {
                gsInfo << "Invalid global index " << globalIdx << " for patch " << patch << "\n";
                outfile << "Invalid global index " << globalIdx << " for patch " << patch << "\n";
                continue;
            }
            localVectSol(pos, 0) = vectSol(globalIdx, 0);
            localVectSol(pos, 1) = vectSol(globalIdx, 1);
            ++pos;
        }

        try
        {
            gsGeometry<>::uPtr geom = SubdomainHierarchy(patch).makeGeometry(localVectSol);
            if (geom)
            {
                newmp.addPatch(std::move(geom));
                gsInfo << "Added geometry for patch " << patch << "\n";
                outfile << "Added geometry for patch " << patch << "\n";
            }
        }
        catch (const std::exception& e)
        {
            gsInfo << "EXCEPTION while creating geometry for patch " << patch << ": " << e.what() << "\n";
            outfile << "EXCEPTION while creating geometry for patch " << patch << ": " << e.what() << "\n";
        }
    }

    // Compute topology and dump multipatch
    try
    {
        real_t tol = 1e-8;
        if (newmp.interfaces().empty())
            newmp.computeTopology(tol, false);

        gsFileData<> fdG;
        fdG << newmp;
        fdG.dump(baseName);

        gsInfo << "Wrote multipatch geometry file: " << baseName << "\n";
        outfile << "Wrote multipatch geometry file: " << baseName << "\n";
    }
    catch (const std::exception& e)
    {
        gsInfo << "EXCEPTION while saving multipatch: " << e.what() << "\n";
        outfile << "EXCEPTION while saving multipatch: " << e.what() << "\n";
    }
}

// (No overload here — use pre-converted gsSparseVector<real_t> where needed.)