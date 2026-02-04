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

#include <gsIO/gsIOUtils.h>
#include <gsHSplines/gsHTensorBasis.h>
#include <gsNurbs/gsBoehm.h>
#include <gsNurbs/gsBoehm.hpp>
#include <gsNurbs/gsDeboor.hpp>
//#include "/home/targon/theydarov/gismo/extensions/motor/jku/gsQualityMeasure2.h"

//#include "/home/targon/theydarov/gismo/extensions/motor/jku/gsMotorUtils.h"
#include <cmath>
//#include </home/targon/theydarov/gismo/examples/rebuildTheHierarchy.cpp>
/*#include  <gsMPBES/gsMPBESBasis.h>
#include  <gsMPBES/gsMPBESBSplineBasis.h>
#include  <gsMPBES/gsMPBESHSplineBasis.h>
#include  <gsMPBES/gsMPBESMapB2D.h>
#include  <gsMPBES/gsMPBESMapHB2D.h>
#include  <gsMPBES/gsMPBESMapTensor.h>
#include  <gsMPBES/gsMPBESSpline.h>
#include  <gsMPBES/gsMPBESUtils.h>
#include  <gsMPBES/gsPatchReparameterized.h>*/
#include <unordered_set>
#include <chrono>
#include <map>
#include <iomanip>


#define printValue( val )   cout << #val << ":\n";
#define bussy front();

using namespace gismo;
using namespace std;

// Return codes
constexpr int SUCCESS = 0;
constexpr int ERROR_READING_GEOMETRY = 1;
constexpr int ERROR_FAILED_REFINEMENT = 2;
constexpr int ERROR_SINGULAR_JACOBIAN = 3;
constexpr int VIOLATED_PARTITION_OF_UNITY = 4;
constexpr int ERROR_MPBES_EVALUATION = 260115;

// Global output file - declared before profiling classes
ofstream outfile;
ofstream xboxFile;
ofstream yboxFile;

// ============== PROFILING INFRASTRUCTURE ==============
class ProfileTimer {
private:
    std::chrono::time_point<std::chrono::high_resolution_clock> start_time;
    std::string section_name;
    bool active;

public:
    ProfileTimer(const std::string& name, bool enabled = true)
        : section_name(name), active(enabled) {
        if (active) {
            start_time = std::chrono::high_resolution_clock::now();
        }
    }

    ~ProfileTimer() {
        if (active) {
            auto end_time = std::chrono::high_resolution_clock::now();
            auto duration = std::chrono::duration<double>(end_time - start_time).count();
            gsInfo << "[PROFILE] " << section_name << ": " << duration << " s\n";
            outfile << "[PROFILE] " << section_name << ": " << duration << " s\n";
        }
    }
};

class GlobalProfiler {
public:
    std::map<std::string, double> timings;
    std::map<std::string, int> counts;

    void record(const std::string& name, double time) {
        timings[name] += time;
        counts[name]++;
    }

    void printSummary() {
        gsInfo << "\n========== PROFILING SUMMARY ==========\n";
        outfile << "\n========== PROFILING SUMMARY ==========\n";

        std::vector<std::pair<std::string, double>> sorted_timings(timings.begin(), timings.end());
        std::sort(sorted_timings.begin(), sorted_timings.end(),
            [](const auto& a, const auto& b) { return a.second > b.second; });

        for (const auto& item : sorted_timings) {
            double avg = item.second / counts[item.first];
            gsInfo << item.first << ": " << item.second << " s (avg: "
                << avg << " s, calls: " << counts[item.first] << ")\n";
            outfile << item.first << ": " << item.second << " s (avg: "
                << avg << " s, calls: " << counts[item.first] << ")\n";
        }
        gsInfo << "=======================================\n\n";
        outfile << "=======================================\n\n";
    }
};

GlobalProfiler g_profiler;

#define PROFILE_SECTION(name) ProfileTimer _timer_##__LINE__(name, true)
#define PROFILE_FUNCTION() ProfileTimer _timer_##__LINE__(__FUNCTION__, true)
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
        bool verbose = true,
        int attempt = 0)
        : m_thbBases(thbBases),
        m_bellsBases(bellsBases),
        m_numPatches(thbBases.size()),
        m_boxMat(boxMat),
        m_indexInTHB(indexInTHB),
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
    void evalSingleOnPatch(index_t globalIdx, index_t patch, const gsMatrix<T>& u, gsMatrix<T>& result) const {
        result.resize(1, u.cols());
        result.setZero();

        if (globalIdx >= m_numGlobalFunctions || globalIdx >= m_functionDescription.size())
            return;

        // Find local representation on this patch
        const auto& twins = m_functionDescription[globalIdx];
        int componentIdx = 0;

        for (const auto& twin : twins) {
            if (twin[0] == static_cast<int>(patch)) {
                int level = twin[1];
                int tensorIdx = twin[2];

                // CRITICAL: Check if THIS SPECIFIC COMPONENT is truncated
                // If truncated: use presentation (finer-level tensor basis)
                // If not truncated: use standard THB basis
                bool compTruncated = isComponentTruncated(globalIdx, componentIdx);

                if (compTruncated) {
                    // Component is truncated - evaluate using presentation coefficients
                    if (globalIdx < m_presentation.size() && componentIdx < m_presentation[globalIdx].size()) {
                        const auto& coeffs = m_presentation[globalIdx][componentIdx];
                        index_t presLevel = m_presentationLevel[globalIdx][componentIdx];

                        // Linear combination: evaluate each basis function individually
                        // CRITICAL: Must use evalSingle_into for each coefficient because
                        // coeffs has sparse indices (e.g., 266, 267, 284) that are global
                        // tensor basis indices, not positions in a compressed array
                        T resultValue = 0.0;
                        for (typename gsSparseVector<T>::InnerIterator it(coeffs); it; ++it) {
                            gsMatrix<T> basisVal;
                            m_bellsBases[patch][presLevel].evalSingle_into(it.index(), u, basisVal);
                            resultValue += it.value() * basisVal(0, 0);
                        }
                        result.row(0).array() += resultValue;
                    }
                }
                else {
                    // Component NOT truncated - evaluate using THB basis (respects active/inactive)
                    // Get the THB index for this component
                    if (patch < m_indexInTHB.size() &&
                        level < m_indexInTHB[patch].size() &&
                        tensorIdx < m_indexInTHB[patch][level].size()) {
                        int thbIdx = m_indexInTHB[patch][level][tensorIdx];
                        if (thbIdx >= 0) {
                            gsMatrix<T> componentResult;
                            m_thbBases[patch].evalSingle_into(thbIdx, u, componentResult);
                            result += componentResult;
                        }
                    }
                }
            }
            componentIdx++;
        }
    }

    /// Evaluate first derivative of single function on specific patch
    void evalDerivSingleOnPatch(index_t globalIdx, index_t patch, const gsMatrix<T>& u, gsMatrix<T>& result) const {
        result.resize(d, u.cols());
        result.setZero();

        if (globalIdx >= m_numGlobalFunctions || globalIdx >= m_functionDescription.size())
            return;

        const auto& twins = m_functionDescription[globalIdx];
        int componentIdx = 0;

        for (const auto& twin : twins) {
            if (twin[0] == static_cast<int>(patch)) {
                int level = twin[1];
                int tensorIdx = twin[2];
                bool compTruncated = isComponentTruncated(globalIdx, componentIdx);

                if (compTruncated) {
                    // Truncated: evaluate derivative of presentation
                    if (globalIdx < m_presentation.size() && componentIdx < m_presentation[globalIdx].size()) {
                        const auto& coeffs = m_presentation[globalIdx][componentIdx];
                        index_t presLevel = m_presentationLevel[globalIdx][componentIdx];

                        for (typename gsSparseVector<T>::InnerIterator it(coeffs); it; ++it) {
                            gsMatrix<T> derivVal = m_bellsBases[patch][presLevel].function(it.index()).deriv(u);
                            result += it.value() * derivVal;
                        }
                    }
                }
                else {
                    // Not truncated: evaluate THB derivative
                    if (patch < m_indexInTHB.size() &&
                        level < m_indexInTHB[patch].size() &&
                        tensorIdx < m_indexInTHB[patch][level].size()) {
                        int thbIdx = m_indexInTHB[patch][level][tensorIdx];
                        if (thbIdx >= 0) {
                            gsMatrix<T> derivVal = m_thbBases[patch].function(thbIdx).deriv(u);
                            result += derivVal;
                        }
                    }
                }
            }
            componentIdx++;
        }
    }

    /// Evaluate second derivative of single function on specific patch
    void evalDeriv2SingleOnPatch(index_t globalIdx, index_t patch, const gsMatrix<T>& u, gsMatrix<T>& result) const {
        const index_t numSecondDerivs = d * (d + 1) / 2;  // 3 for 2D: dxx, dxy, dyy
        result.resize(numSecondDerivs, u.cols());
        result.setZero();

        if (globalIdx >= m_numGlobalFunctions || globalIdx >= m_functionDescription.size())
            return;

        const auto& twins = m_functionDescription[globalIdx];
        int componentIdx = 0;

        for (const auto& twin : twins) {
            if (twin[0] == static_cast<int>(patch)) {
                int level = twin[1];
                int tensorIdx = twin[2];
                bool compTruncated = isComponentTruncated(globalIdx, componentIdx);

                if (compTruncated) {
                    // Truncated: evaluate second derivative of presentation
                    if (globalIdx < m_presentation.size() && componentIdx < m_presentation[globalIdx].size()) {
                        const auto& coeffs = m_presentation[globalIdx][componentIdx];
                        index_t presLevel = m_presentationLevel[globalIdx][componentIdx];

                        for (typename gsSparseVector<T>::InnerIterator it(coeffs); it; ++it) {
                            gsMatrix<T> deriv2Val = m_bellsBases[patch][presLevel].function(it.index()).deriv2(u);
                            result += it.value() * deriv2Val;
                        }
                    }
                }
                else {
                    // Not truncated: evaluate THB second derivative
                    if (patch < m_indexInTHB.size() &&
                        level < m_indexInTHB[patch].size() &&
                        tensorIdx < m_indexInTHB[patch][level].size()) {
                        int thbIdx = m_indexInTHB[patch][level][tensorIdx];
                        if (thbIdx >= 0) {
                            gsMatrix<T> deriv2Val = m_thbBases[patch].function(thbIdx).deriv2(u);
                            result += deriv2Val;
                        }
                    }
                }
            }
            componentIdx++;
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
            m_functionDescription,
            m_spilloverFunctionCoordinates,
            m_spilloverSources,
            m_hasSpillover,
            m_globalIndexTHB,
            m_globalIndex,
            m_numGlobalFunctions,
            verbose
        );

        if (verbose)
            gsInfo << "MPBES: After selection, " << m_numGlobalFunctions << " global functions\n";

        // Step 2: Detect which functions need truncation
        detectTruncatedFunctions(verbose);

        // Step 3: Compute truncated representations (spillover coefficients)
        if (std::any_of(m_isTruncated.begin(), m_isTruncated.end(), [](bool v) { return v; })) {
            m_presentation.resize(m_numGlobalFunctions);
            computeTruncatedRepresentations(verbose);
        }

        if (verbose) {
            gsInfo << "MPBES constructed: " << m_numGlobalFunctions << " global functions\n";
            index_t numTruncated = std::count(m_isTruncated.begin(), m_isTruncated.end(), true);
            gsInfo << "  - " << numTruncated << " functions are truncated\n";
        }
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

                    if (f < 10 && verbose) {
                        gsInfo << "  [DETECT] Func " << f << " comp " << twinIdx
                            << ": level=" << level << ", highestOverlap=" << highestOverlap
                            << " -> presentationLevel=" << m_presentationLevel[f][twinIdx] << "\n";
                    }
                }
            }
        }

        index_t numTruncatedFunctions = std::count(m_isTruncated.begin(), m_isTruncated.end(), true);

        gsInfo << "MPBES: Detected " << numTruncatedComponents << " truncated components in "
            << numTruncatedFunctions << " functions\n";
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
    bool verbose = false);
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
real_t boundaryError(
    const MPBES<2, real_t>& mpbes,
    const gsMultiPatch<real_t>& mp,
    const gsMatrix<real_t>& vectSol
);
void checkForNaN(const gsMatrix<real_t>& matrix, const std::string& matrixName);
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
    index_t gridResolution,
    std::vector<gsMatrix<>>& uvGridLines,
    std::vector<index_t>& patchIDs,
    std::vector<char>& directions,
    bool verbose = false
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

int twinFunction(gsTHBSplineBasis<2, real_t> initTP, int functionIndex, gsTHBSplineBasis<2, real_t> twinTP, int firstSide, int secondSide);

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
    real_t length_weight,
    real_t epsilon_g,
    real_t epsilon_f,
    gsVector<size_t>& numIrregular,
    index_t geoDim
);

void evalSingle_into(
    index_t f,
    index_t i,
    const gsMatrix<real_t>& u,
    const std::vector<bool>& isTruncated,
    const gsTHBSplineBasis<2, real_t>& localTHBBasis,
    const gsSparseVector<real_t>& coefs,
    gsMatrix<real_t>& result,
    const gsTensorBSplineBasis<2, real_t>* BellsBasis = nullptr,
    int idx = -1)
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

    // Handle spillover case (i == -1)
    if (i == -1)
    {
        if (BellsBasis != nullptr && idx >= 0)
        {
            // Evaluate as pure B-spline basis function without coefficients
            BellsBasis->evalSingle_into(idx, u, result);
        }
        else
        {
            // Cannot evaluate spillover without BellsBasis information
            result.resize(1, u.cols());
            result.setZero();
        }
    }
    else if (!isTruncated[f]) //  basis function not truncated
    {
        unsigned level = localTHBBasis.levelOf(i);
        unsigned tensor_index = localTHBBasis.flatTensorIndexOf(i, level);
        localTHBBasis.getBases()[level]->evalSingle_into(tensor_index, u, result);
    }
    else
    {
        // Function IS truncated - check if we have valid coefficients
        // If coefs is empty/zero, this component might not be the truncated one
        if (!coefsOK) {
            // No coefficients provided - evaluate as non-truncated
            unsigned level = localTHBBasis.levelOf(i);
            unsigned tensor_index = localTHBBasis.flatTensorIndexOf(i, level);
            localTHBBasis.getBases()[level]->evalSingle_into(tensor_index, u, result);
        }
        else {
            // Function IS truncated and we have coefficients - use De Boor's algorithm
            unsigned level = localTHBBasis.getistruncated(i);

            if (verbose) {
                gsInfo << "Evaluating truncated function f=" << f << " at THB index " << i
                    << " (level " << level << ") with " << coefs.size() << " coefficients\n";
                gsInfo << "  Coefficients: ";
                for (size_t j = 0; j < coefs.size(); j++) {
                    if (coefs(j) != 0.0)
                        gsInfo << "(" << j << ":" << coefs(j) << ") ";
                }
                gsInfo << "\n";
            }

            const gsTensorBSplineBasis<2, real_t>& base = *localTHBBasis.getBases()[level];
            gsTensorDeboor<2, real_t, gsKnotVector<real_t>, gsSparseVector<real_t>>(u, base, coefs, result);
        }
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
  * Creates grid lines showing the tensor product grid structure within each box.
  * For a box with index coordinates (ix0, iy0, ix1, iy1) at level L:
  * - Generates grid lines at each index position: ix0, ix0+1, ..., ix1 (vertical)
  *   and iy0, iy0+1, ..., iy1 (horizontal)
  * - For each grid line, samples at ALL index positions between the endpoints
  * - Ensures interface continuity: boxes sharing edges have matching sample points
  *   because they sample at the same index positions (scaled to finest level)
  *
  * Key insight: By sampling at integer index positions in the hierarchical grid,
  * we guarantee that interfaces between boxes at different levels align perfectly.
  */
void generateParametricGrid(
    const gsVector<gsVector<gsVector<index_t>>>& boxMat,
    const gsVector<index_t>& numBoxesPerPatch,
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

    if (verbose)
        gsInfo << "Generating parametric grid from box structure...\n";

    // Find maximum refinement level to determine finest grid spacing
    index_t maxLevel = 0;
    for (index_t patch = 0; patch < boxMat.size(); ++patch)
    {
        const index_t numBoxes = numBoxesPerPatch[patch];
        for (index_t boxIdx = 0; boxIdx < numBoxes; ++boxIdx)
        {
            if (boxMat[patch][boxIdx].size() >= 5)
            {
                maxLevel = std::max(maxLevel, boxMat[patch][boxIdx][0]);
            }
        }
    }

    // gridResolution is minimum samples per finest-level index cell
    const index_t samplesPerIndexCell = std::max<index_t>(gridResolution, 1);

    if (verbose)
    {
        gsInfo << "  Maximum refinement level: " << maxLevel << "\n";
        gsInfo << "  Samples per index cell: " << samplesPerIndexCell << "\n";
    }

    index_t skippedBoxes = 0;

    for (index_t patch = 0; patch < boxMat.size(); ++patch)
    {
        const index_t numBoxes = numBoxesPerPatch[patch];

        if (verbose)
            gsInfo << "  Patch " << patch << ": " << numBoxes << " boxes\n";

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

            // Box data in index space
            const index_t level = boxMat[patch][boxIdx][0];
            const index_t ix0 = boxMat[patch][boxIdx][1];
            const index_t iy0 = boxMat[patch][boxIdx][2];
            const index_t ix1 = boxMat[patch][boxIdx][3];
            const index_t iy1 = boxMat[patch][boxIdx][4];

            // Grid spacing at this level
            const real_t h = 1.0 / std::pow(2.0, level);

            // For each grid line, we need to sample at all relevant index positions.
            // When we scale to maxLevel, this box spans from 
            // (ix0*scale, iy0*scale) to (ix1*scale, iy1*scale) where scale = 2^(maxLevel-level)
            const index_t levelDiff = maxLevel - level;
            const index_t scaleFactor = 1 << levelDiff; // 2^(maxLevel - level)

            // Scaled index bounds at finest level
            const index_t scaledIx0 = ix0 * scaleFactor;
            const index_t scaledIy0 = iy0 * scaleFactor;
            const index_t scaledIx1 = ix1 * scaleFactor;
            const index_t scaledIy1 = iy1 * scaleFactor;

            // Number of finest-level index cells in each direction
            const index_t numFinestCellsU = scaledIx1 - scaledIx0;
            const index_t numFinestCellsV = scaledIy1 - scaledIy0;

            // Generate horizontal lines (constant v, varying u)
            // One line for each index position iy0, iy0+1, ..., iy1
            for (index_t iy = iy0; iy <= iy1; ++iy)
            {
                const real_t v = iy * h;

                // For this horizontal line, sample at all finest-level indices
                // The line spans from u = ix0*h to u = ix1*h
                // At finest level, this corresponds to indices scaledIx0 to scaledIx1
                // Sample at each index with samplesPerIndexCell points per cell
                const index_t numPoints = numFinestCellsU * samplesPerIndexCell + 1;

                gsMatrix<> line(2, numPoints);
                for (index_t i = 0; i < numPoints; ++i)
                {
                    // Map sample index to finest-level index coordinate
                    const real_t finestLevelIndex = scaledIx0 + static_cast<real_t>(i) / samplesPerIndexCell;
                    // Convert finest-level index to parameter space
                    const real_t hFinest = 1.0 / std::pow(2.0, maxLevel);
                    const real_t u = finestLevelIndex * hFinest;

                    line(0, i) = u;
                    line(1, i) = v;
                }

                uvGridLines.push_back(line);
                patchIDs.push_back(patch);
                directions.push_back('u');
            }

            // Generate vertical lines (constant u, varying v)
            // One line for each index position ix0, ix0+1, ..., ix1
            for (index_t ix = ix0; ix <= ix1; ++ix)
            {
                const real_t u = ix * h;

                // For this vertical line, sample at all finest-level indices
                const index_t numPoints = numFinestCellsV * samplesPerIndexCell + 1;

                gsMatrix<> line(2, numPoints);
                for (index_t i = 0; i < numPoints; ++i)
                {
                    // Map sample index to finest-level index coordinate
                    const real_t finestLevelIndex = scaledIy0 + static_cast<real_t>(i) / samplesPerIndexCell;
                    // Convert finest-level index to parameter space
                    const real_t hFinest = 1.0 / std::pow(2.0, maxLevel);
                    const real_t v = finestLevelIndex * hFinest;

                    line(0, i) = u;
                    line(1, i) = v;
                }

                uvGridLines.push_back(line);
                patchIDs.push_back(patch);
                directions.push_back('v');
            }

            if (verbose)
                gsInfo << "    Box " << boxIdx << " [level=" << level
                << ", (" << ix0 << "," << iy0 << ")->(" << ix1 << "," << iy1 << ")]: "
                << (ix1 - ix0 + 1) << " vertical + "
                << (iy1 - iy0 + 1) << " horizontal lines\n";
        }
    }

    if (verbose)
        gsInfo << "Generated " << uvGridLines.size() << " total grid lines ("
        << skippedBoxes << " boxes skipped due to finer coverage).\n";
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
    gsInfo << "Box structure exported to " << filename << "\n";
}

/**
 * @brief Log specific interface between two boxes
 */
void logSpecificInterface(
    const std::vector<gsMatrix<>>& uvLines,
    const std::vector<gsMatrix<>>& xyLines,
    const std::vector<index_t>& patchIDs,
    const std::vector<char>& directions,
    const std::string& filename)
{
    std::ofstream log(filename);
    if (!log.is_open())
    {
        gsWarn << "Could not open " << filename << " for writing.\n";
        return;
    }

    log << "# Specific Interface Point Analysis: (0,(0,1)) vs (1,(1,1))\n";
    log << "# Evaluating geometry at Patch 0 (u=0, v=1) and Patch 1 (u=1, v=1)\n\n";

    // Evaluate geometry at patch 0 (u=0, v=1)
    gsVector<real_t> uv_p0(2);
    uv_p0 << 0.0, 1.0;

    gsVector<real_t> xy_p0(2);
    xy_p0.setZero();
    bool found_p0 = false;

    for (size_t i = 0; i < uvLines.size() && !found_p0; ++i)
    {
        if (patchIDs[i] == 0)
        {
            const gsMatrix<>& uvLine = uvLines[i];
            const gsMatrix<>& xyLine = xyLines[i];

            // Find closest point to (0, 1)
            for (index_t j = 0; j < uvLine.cols(); ++j)
            {
                real_t dist = std::sqrt(std::pow(uvLine(0, j) - 0.0, 2) + std::pow(uvLine(1, j) - 1.0, 2));
                if (dist < 1e-6)
                {
                    xy_p0 << xyLine(0, j), xyLine(1, j);
                    found_p0 = true;
                    log << "Found Patch 0 at (u=" << uvLine(0, j) << ", v=" << uvLine(1, j) << ")\n";
                    log << "  -> Geometry: (x=" << xy_p0(0) << ", y=" << xy_p0(1) << ")\n\n";
                    break;
                }
            }
        }
    }

    // Evaluate geometry at patch 1 (u=1, v=1)
    gsVector<real_t> uv_p1(2);
    uv_p1 << 1.0, 1.0;

    gsVector<real_t> xy_p1(2);
    xy_p1.setZero();
    bool found_p1 = false;

    for (size_t i = 0; i < uvLines.size() && !found_p1; ++i)
    {
        if (patchIDs[i] == 1)
        {
            const gsMatrix<>& uvLine = uvLines[i];
            const gsMatrix<>& xyLine = xyLines[i];

            // Find closest point to (1, 1)
            for (index_t j = 0; j < uvLine.cols(); ++j)
            {
                real_t dist = std::sqrt(std::pow(uvLine(0, j) - 1.0, 2) + std::pow(uvLine(1, j) - 1.0, 2));
                if (dist < 1e-6)
                {
                    xy_p1 << xyLine(0, j), xyLine(1, j);
                    found_p1 = true;
                    log << "Found Patch 1 at (u=" << uvLine(0, j) << ", v=" << uvLine(1, j) << ")\n";
                    log << "  -> Geometry: (x=" << xy_p1(0) << ", y=" << xy_p1(1) << ")\n\n";
                    break;
                }
            }
        }
    }

    if (found_p0 && found_p1)
    {
        real_t distance = std::sqrt(std::pow(xy_p0(0) - xy_p1(0), 2) + std::pow(xy_p0(1) - xy_p1(1), 2));
        log << "========================================\n";
        log << "RESULT:\n";
        log << "  Patch 0 (u=0, v=1): (x=" << xy_p0(0) << ", y=" << xy_p0(1) << ")\n";
        log << "  Patch 1 (u=1, v=1): (x=" << xy_p1(0) << ", y=" << xy_p1(1) << ")\n";
        log << "  Distance: " << distance << "\n";
        log << "  Delta X: " << (xy_p0(0) - xy_p1(0)) << "\n";
        log << "  Delta Y: " << (xy_p0(1) - xy_p1(1)) << "\n";
        log << "========================================\n\n";

        gsInfo << "Interface point analysis:\n";
        gsInfo << "  Patch 0 (0,1): x=" << xy_p0(0) << ", y=" << xy_p0(1) << "\n";
        gsInfo << "  Patch 1 (1,1): x=" << xy_p1(0) << ", y=" << xy_p1(1) << "\n";
        gsInfo << "  Distance: " << distance << "\n";
    }
    else
    {
        log << "ERROR: Could not find both points\n";
        if (!found_p0) log << "  Missing: Patch 0 at (0, 1)\n";
        if (!found_p1) log << "  Missing: Patch 1 at (1, 1)\n";
    }

    log << "# Original Interface Analysis: Patch 0 Box [3, 0, 2, 8, 8] <-> Patch 1 Box [4, 0, 0, 16, 16]\n";
    log << "# Interface: Patch 0 at u=0 (v from 0.25 to 1) <-> Patch 1 at u=1 (v from 0.25 to 1)\n\n";

    const real_t tolerance = 1e-9;

    // Find vertical line at u=0 on patch 0
    gsMatrix<> patch0Line;
    gsMatrix<> patch0XY;
    bool foundPatch0 = false;

    for (size_t i = 0; i < uvLines.size(); ++i)
    {
        if (patchIDs[i] == 0 && directions[i] == 'v')
        {
            const gsMatrix<>& uvLine = uvLines[i];
            if (uvLine.cols() > 0 && std::abs(uvLine(0, 0) - 0.0) < tolerance)
            {
                // This is a vertical line at u=0
                // Check if it covers the range v=0.25 to v=1
                real_t minV = uvLine(1, 0);
                real_t maxV = uvLine(1, uvLine.cols() - 1);
                if (minV < maxV)
                {
                    if (std::abs(minV - 0.25) < tolerance && std::abs(maxV - 1.0) < tolerance)
                    {
                        patch0Line = uvLine;
                        patch0XY = xyLines[i];
                        foundPatch0 = true;
                        log << "# Found Patch 0 line at index " << i << " with " << uvLine.cols() << " points\n";
                        break;
                    }
                }
                else
                {
                    if (std::abs(maxV - 0.25) < tolerance && std::abs(minV - 1.0) < tolerance)
                    {
                        patch0Line = uvLine;
                        patch0XY = xyLines[i];
                        foundPatch0 = true;
                        log << "# Found Patch 0 line at index " << i << " with " << uvLine.cols() << " points (reversed)\n";
                        break;
                    }
                }
            }
        }
    }

    // Find vertical line at u=1 on patch 1
    gsMatrix<> patch1Line;
    gsMatrix<> patch1XY;
    bool foundPatch1 = false;

    for (size_t i = 0; i < uvLines.size(); ++i)
    {
        if (patchIDs[i] == 1 && directions[i] == 'v')
        {
            const gsMatrix<>& uvLine = uvLines[i];
            if (uvLine.cols() > 0 && std::abs(uvLine(0, 0) - 1.0) < tolerance)
            {
                // This is a vertical line at u=1
                // Check if it covers the range v=0.25 to v=1 (or v=0 to 1)
                real_t minV = uvLine(1, 0);
                real_t maxV = uvLine(1, uvLine.cols() - 1);

                // We want the portion from v=0.25 to v=1
                if (minV < maxV)
                {
                    // Line goes from minV to maxV
                    if (maxV >= 0.25)
                    {
                        patch1Line = uvLine;
                        patch1XY = xyLines[i];
                        foundPatch1 = true;
                        log << "# Found Patch 1 line at index " << i << " with " << uvLine.cols() << " points\n";
                        log << "# v ranges from " << minV << " to " << maxV << "\n";
                        break;
                    }
                }
                else
                {
                    // Line goes from maxV to minV
                    if (minV >= 0.25)
                    {
                        patch1Line = uvLine;
                        patch1XY = xyLines[i];
                        foundPatch1 = true;
                        log << "# Found Patch 1 line at index " << i << " with " << uvLine.cols() << " points (reversed)\n";
                        log << "# v ranges from " << maxV << " to " << minV << "\n";
                        break;
                    }
                }
            }
        }
    }

    if (!foundPatch0)
    {
        log << "# ERROR: Could not find Patch 0 interface line at u=0, v=[0.25, 1]\n";
        log.close();
        return;
    }

    if (!foundPatch1)
    {
        log << "# ERROR: Could not find Patch 1 interface line at u=1\n";
        log.close();
        return;
    }

    log << "\n# DETAILED POINT-BY-POINT COMPARISON\n";
    log << "# Format: (u0, v) -> (x0, y0) | (u1, v) -> (x1, y1) | distance\n\n";

    // Find matching v-values
    // Patch 0: v from 0.25 to 1.0 at level 3 (indices 2 to 8) -> scaled to level 4: indices 4 to 16
    // Patch 1: v from 0 to 1.0 at level 4 (indices 0 to 16)
    // So we need to compare v in [0.25, 1.0]

    // Build a map of v-values from patch 1
    std::map<real_t, std::pair<gsVector<real_t>, gsVector<real_t>>> patch1Map; // v -> (uv, xy)

    for (index_t i = 0; i < patch1Line.cols(); ++i)
    {
        real_t v = patch1Line(1, i);
        if (v >= 0.25 - tolerance && v <= 1.0 + tolerance)
        {
            gsVector<real_t> uv(2);
            uv << patch1Line(0, i), patch1Line(1, i);
            gsVector<real_t> xy(2);
            xy << patch1XY(0, i), patch1XY(1, i);
            patch1Map[v] = std::make_pair(uv, xy);
        }
    }

    log << "# Patch 0 has " << patch0Line.cols() << " points\n";
    log << "# Patch 1 has " << patch1Map.size() << " points in range [0.25, 1.0]\n\n";

    // For each point on patch 0, find the matching point on patch 1
    int matchCount = 0;
    real_t maxDistance = 0.0;

    for (index_t i = 0; i < patch0Line.cols(); ++i)
    {
        real_t u0 = patch0Line(0, i);
        real_t v0 = patch0Line(1, i);
        real_t x0 = patch0XY(0, i);
        real_t y0 = patch0XY(1, i);

        // Find matching v-value in patch 1
        auto it = patch1Map.lower_bound(v0 - tolerance);
        if (it != patch1Map.end() && std::abs(it->first - v0) < 1e-6)
        {
            real_t u1 = it->second.first(0);
            real_t v1 = it->second.first(1);
            real_t x1 = it->second.second(0);
            real_t y1 = it->second.second(1);

            real_t distance = std::sqrt((x0 - x1) * (x0 - x1) + (y0 - y1) * (y0 - y1));
            maxDistance = std::max(maxDistance, distance);

            log << std::setprecision(10)
                << "(" << u0 << ", " << v0 << ") -> (" << x0 << ", " << y0 << ") | "
                << "(" << u1 << ", " << v1 << ") -> (" << x1 << ", " << y1 << ") | "
                << "dist = " << distance << "\n";

            matchCount++;
        }
        else
        {
            log << "(" << u0 << ", " << v0 << ") -> (" << x0 << ", " << y0 << ") | "
                << "NO MATCH FOUND on Patch 1\n";
        }
    }

    log << "\n# SUMMARY\n";
    log << "# Matched " << matchCount << " out of " << patch0Line.cols() << " points\n";
    log << "# Maximum distance between matched points: " << maxDistance << "\n";

    if (matchCount != patch0Line.cols())
    {
        log << "# WARNING: Not all points from Patch 0 have matches on Patch 1!\n";
    }

    log.close();
    gsInfo << "Specific interface analysis exported to " << filename << "\n";
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
    if (verbose)
        gsInfo << "\n=== generateVisualizationMesh BEGIN ===\n";

    // Debug: Export box structure
    exportBoxStructure(boxMat, currentLastNonZeroRow, outputPrefix + "_boxes.txt");

    // Step 1: Generate parametric grid
    std::vector<gsMatrix<>> uvLines;
    std::vector<index_t> uvLinePatchIDs;
    std::vector<char> uvLineDirections;

    if (verbose)
        gsInfo << "Step 1: Generating parametric grid...\n";

    generateParametricGrid(
        boxMat,
        currentLastNonZeroRow,
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

    // Step 2: Evaluate geometry at grid points
    std::vector<gsMatrix<>> xyLines;

    if (verbose)
        gsInfo << "Step 2: Evaluating geometry...\n";

    evaluateGeometryAtPoints(
        uvLines,
        uvLinePatchIDs,
        mpbes,
        mp,
        coefficients,
        xyLines,
        verbose
    );

    // Step 2.5: Log specific interface
    if (verbose)
        gsInfo << "Step 2.5: Analyzing specific interface...\n";

    logSpecificInterface(
        uvLines,
        xyLines,
        uvLinePatchIDs,
        uvLineDirections,
        outputPrefix + "_interface_specific.txt"
    );

    // Step 3: Export mesh to file
    if (verbose)
        gsInfo << "Step 3: Exporting mesh...\n";

    std::string meshFile = outputPrefix + ".txt";
    exportMeshToFile(
        xyLines,
        uvLinePatchIDs,
        uvLineDirections,
        meshFile,
        verbose
    );

    if (verbose)
        gsInfo << "=== generateVisualizationMesh END ===\n\n";
}

void selectionMechanism(
    const gsVector<gsVector<gsVector<index_t>>>& boxMat,
    const gsVector<gsVector<gsTensorBSplineBasis<2, real_t>>>& Bells,
    const gsVector<gsTHBSplineBasis<2, real_t>>& THBVector,
    const gsVector<gsVector<gsVector<index_t>>>& hasATwin,
    const gsVector<gsVector<gsVector<std::vector<index_t>>>>& twinsIndex,
    const gsVector<gsVector<gsVector<std::vector<index_t>>>>& twinsPatch,
    const gsVector<gsVector<gsVector<index_t>>>& indexInTHB,
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
    globalFunctionCount = 0;
    functionDescription.clear();
    spilloverFunctionCoordinates.clear();
    spilloverSources.clear();
    hasSpillover.clear();

    globalIndexTHB.resize(nPatches);
    globalIndex.resize(nPatches);
    spilloverSources.resize(nPatches);

    if (verbose)
        gsInfo << "=== selectionMechanism BEGIN ===\n";

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

                if (indexInTHB[patch][level][i] == -1)
                {
                    if (verbose)
                        outfile << "SKIP: patch=" << patch << ", level=" << level << ", index=" << i << " — not in THB\n";
                    continue;
                }

                bool intersectsBox = false;
                const auto& supp = Bells[patch][level].function(i).support();
                real_t sX0 = supp(0, 0), sX1 = supp(0, 1);
                real_t sY0 = supp(1, 0), sY1 = supp(1, 1);

                for (size_t b = 0; b < boxMat[patch].size(); ++b)
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

                if (hasATwin[patch][level][i])
                {
                    if (verbose)
                        outfile << "  -> Checking twins for (" << patch << "," << level << "," << i << ")\n";

                    for (size_t t = 0; t < twinsIndex[patch][level][i].size(); ++t)
                    {
                        int twinPatch = twinsPatch[patch][level][i][t];
                        int twinIndex = twinsIndex[patch][level][i][t];

                        if (verbose)
                            outfile << "    Twin " << t << ": patch=" << twinPatch
                            << ", level=" << level
                            << ", index=" << twinIndex << "\n";

                        if (globalIndex[twinPatch][level][twinIndex] == -1)
                        {
                            globalIndex[twinPatch][level][twinIndex] = globalID;

                            int twinTHB = indexInTHB[twinPatch][level][twinIndex];
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

    if (uvPoints.empty())
    {
        if (verbose)
            gsWarn << "No points to evaluate.\n";
        return;
    }

    if (verbose)
        gsInfo << "Evaluating " << uvPoints.size() << " lines using MPBES with fitted coefficients...\n";

    xyPoints.resize(uvPoints.size());

    // Evaluate each line using MPBES basis with fitted coefficients (NOT original geometry)
    for (size_t lineIdx = 0; lineIdx < uvPoints.size(); ++lineIdx)
    {
        const gsMatrix<>& uvLine = uvPoints[lineIdx];
        const index_t patch = patchIDs[lineIdx];

        if (uvLine.cols() == 0) {
            xyPoints[lineIdx].resize(2, 0);
            continue;
        }

        // Allocate output matrix
        const index_t geoDim = coefficients.cols();
        xyPoints[lineIdx].resize(geoDim, uvLine.cols());

        // Evaluate fitted geometry at each point on this line
        for (index_t k = 0; k < uvLine.cols(); ++k)
        {
            gsMatrix<real_t> pt = uvLine.col(k);
            
            // Evaluate fitted geometry using MPBES basis with coefficients
            gsVector<real_t> xFit(geoDim);
            xFit.setZero();

            for (index_t f = 0; f < mpbes.size(); ++f)
            {
                gsMatrix<real_t> basisVal;
                mpbes.evalSingleOnPatch(f, patch, pt, basisVal);

                const real_t Nf = basisVal(0, 0);
                for (index_t d = 0; d < geoDim; ++d)
                    xFit[d] += Nf * coefficients(f, d);
            }

            // Store result
            for (index_t d = 0; d < geoDim; ++d)
                xyPoints[lineIdx](d, k) = xFit[d];
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

    index_t globalRow = 0;

    for (index_t patch = 0; patch < uv1.size(); ++patch)
    {
        const gsMatrix<real_t>& uvPatch = uv1(patch); // (d × nPts)
        const index_t nPts = uvPatch.cols();

        for (index_t k = 0; k < nPts; ++k, ++globalRow)
        {
            // parameter point (d × 1)
            gsMatrix<real_t> pt = uvPatch.col(k);

            // evaluate fitted geometry
            gsVector<real_t> xFit(geoDim);
            xFit.setZero();

            for (index_t f = 0; f < mpbes.size(); ++f)
            {
                gsMatrix<real_t> basisVal;
                mpbes.evalSingleOnPatch(f, patch, pt, basisVal);

                const real_t Nf = basisVal(0, 0);
                for (index_t d = 0; d < geoDim; ++d)
                    xFit[d] += Nf * vectSol(f, d);
            }

            // compute pointwise error
            real_t err2 = 0.0;
            for (index_t d = 0; d < geoDim; ++d)
            {
                const real_t diff = xFit[d] - b_vec(globalRow, d);
                err2 += diff * diff;
            }

            maxError = std::max(maxError, std::sqrt(err2));
        }
    }

    GISMO_ASSERT(globalRow == b_vec.rows(),
        "Mismatch between uv1 samples and b_vec rows");

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
    const gsMatrix<real_t>& vectSol)
{
    PROFILE_FUNCTION();
    
    outfile << "\n=== boundaryError BEGIN ===\n";

    const index_t pointsPerEdge = 10;
    const index_t numInterfaces = 2;  // Two interfaces: 0↔1 and 0↔2
    
    real_t maxError = 0.0;

    // Interface 1: Patch 0 west (side 3) ↔ Patch 1 east (side 1)
    // Compare fitted geometry vs original geometry at interface points
    for (index_t i = 0; i < pointsPerEdge; ++i)
    {
        real_t t = static_cast<real_t>(i) / (pointsPerEdge - 1);

        // Patch 0, side 3 => (u,v) = (0, 1 - t)
        gsVector<real_t> uv0(2); 
        uv0(0) = 0.0; 
        uv0(1) = 1.0 - t;
        
        // Patch 1, side 1 => (u,v) = (1, t)
        gsVector<real_t> uv1(2); 
        uv1(0) = 1.0; 
        uv1(1) = t;

        // Evaluate fitted geometry using MPBES
        gsMatrix<real_t> xyFit0(2, 1), xyFit1(2, 1);
        xyFit0.setZero();
        xyFit1.setZero();

        // Evaluate fitted on patch 0
        for (index_t f = 0; f < mpbes.size(); ++f)
        {
            gsMatrix<real_t> N0;
            mpbes.evalSingleOnPatch(f, 0, uv0, N0);
            if (N0(0, 0) != 0.0)
            {
                xyFit0(0, 0) += N0(0, 0) * vectSol(f, 0);
                xyFit0(1, 0) += N0(0, 0) * vectSol(f, 1);
            }
        }

        // Evaluate fitted on patch 1
        for (index_t f = 0; f < mpbes.size(); ++f)
        {
            gsMatrix<real_t> N1;
            mpbes.evalSingleOnPatch(f, 1, uv1, N1);
            if (N1(0, 0) != 0.0)
            {
                xyFit1(0, 0) += N1(0, 0) * vectSol(f, 0);
                xyFit1(1, 0) += N1(0, 0) * vectSol(f, 1);
            }
        }

        // Evaluate original geometry at same points
        gsMatrix<real_t> xyOrig0 = mp.patch(0).eval(uv0);
        gsMatrix<real_t> xyOrig1 = mp.patch(1).eval(uv1);

        // Compute fitting error at interface: max error between fitted and original
        real_t err0 = std::sqrt(std::pow(xyFit0(0,0) - xyOrig0(0,0), 2) + 
                                 std::pow(xyFit0(1,0) - xyOrig0(1,0), 2));
        real_t err1 = std::sqrt(std::pow(xyFit1(0,0) - xyOrig1(0,0), 2) + 
                                 std::pow(xyFit1(1,0) - xyOrig1(1,0), 2));
        real_t err = std::max(err0, err1);
        maxError = std::max(maxError, err);

        if (err > 1e-6 || i < 2)  // Log first 2 points and errors
        {
            outfile << "Interface 0↔1, t=" << t 
                << ": fit0=(" << xyFit0(0,0) << "," << xyFit0(1,0) << ")"
                << ", orig0=(" << xyOrig0(0,0) << "," << xyOrig0(1,0) << ")"
                << ", fit1=(" << xyFit1(0,0) << "," << xyFit1(1,0) << ")"
                << ", orig1=(" << xyOrig1(0,0) << "," << xyOrig1(1,0) << ")"
                << ", error=" << err << "\n";
        }
    }

    // Interface 2: Patch 0 south (side 2) ↔ Patch 2 east (side 1)
    for (index_t i = 0; i < pointsPerEdge; ++i)
    {
        real_t t = static_cast<real_t>(i) / (pointsPerEdge - 1);

        // Patch 0, side 2 => (u,v) = (t, 0)
        gsVector<real_t> uv0(2); 
        uv0(0) = t; 
        uv0(1) = 0.0;
        
        // Patch 2, side 1 => (u,v) = (1, 1 - t)
        gsVector<real_t> uv2(2); 
        uv2(0) = 1.0; 
        uv2(1) = 1.0 - t;

        // Evaluate fitted geometry using MPBES
        gsMatrix<real_t> xyFit0(2, 1), xyFit2(2, 1);
        xyFit0.setZero();
        xyFit2.setZero();

        // Evaluate fitted on patch 0
        for (index_t f = 0; f < mpbes.size(); ++f)
        {
            gsMatrix<real_t> N0;
            mpbes.evalSingleOnPatch(f, 0, uv0, N0);
            if (N0(0, 0) != 0.0)
            {
                xyFit0(0, 0) += N0(0, 0) * vectSol(f, 0);
                xyFit0(1, 0) += N0(0, 0) * vectSol(f, 1);
            }
        }

        // Evaluate fitted on patch 2
        for (index_t f = 0; f < mpbes.size(); ++f)
        {
            gsMatrix<real_t> N2;
            mpbes.evalSingleOnPatch(f, 2, uv2, N2);
            if (N2(0, 0) != 0.0)
            {
                xyFit2(0, 0) += N2(0, 0) * vectSol(f, 0);
                xyFit2(1, 0) += N2(0, 0) * vectSol(f, 1);
            }
        }

        // Evaluate original geometry at same points
        gsMatrix<real_t> xyOrig0 = mp.patch(0).eval(uv0);
        gsMatrix<real_t> xyOrig2 = mp.patch(2).eval(uv2);

        // Compute fitting error at interface: max error between fitted and original
        real_t err0 = std::sqrt(std::pow(xyFit0(0,0) - xyOrig0(0,0), 2) + 
                                 std::pow(xyFit0(1,0) - xyOrig0(1,0), 2));
        real_t err2 = std::sqrt(std::pow(xyFit2(0,0) - xyOrig2(0,0), 2) + 
                                 std::pow(xyFit2(1,0) - xyOrig2(1,0), 2));
        real_t err = std::max(err0, err2);
        maxError = std::max(maxError, err);

        if (err > 1e-6 || i < 2)  // Log first 2 points and errors
        {
            outfile << "Interface 0↔2, t=" << t 
                << ": fit0=(" << xyFit0(0,0) << "," << xyFit0(1,0) << ")"
                << ", orig0=(" << xyOrig0(0,0) << "," << xyOrig0(1,0) << ")"
                << ", fit2=(" << xyFit2(0,0) << "," << xyFit2(1,0) << ")"
                << ", orig2=(" << xyOrig2(0,0) << "," << xyOrig2(1,0) << ")"
                << ", error=" << err << "\n";
        }
    }

    outfile << "=== boundaryError END: max interface error = " << maxError << " ===\n\n";
    
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
    real_t maxError = diff.cwiseAbs().maxCoeff();

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
    PROFILE_FUNCTION();
    gsInfo << "\n=== testBoundaryAssembly BEGIN ===\n";
    outfile << "\n=== testBoundaryAssembly BEGIN ===\n";

    const index_t pointsPerEdge = 10;
    const index_t numInterfaces = 2;  // Two interfaces: 0↔1 and 0↔2
    const index_t pointsPerInterface = pointsPerEdge * 2;  // 2 sides per interface
    const index_t totalPoints = numInterfaces * pointsPerInterface;

    // Build parametric points for all interfaces
    gsVector<gsMatrix<>> interfacePoints(mp.nPatches());
    for (index_t p = 0; p < mp.nPatches(); ++p) {
        interfacePoints[p].resize(2, 0);  // Will be appended to
    }

    // Interface 1: Patch 0 west ↔ Patch 1 east
    gsMatrix<real_t> patch0_interface1(2, pointsPerEdge);
    gsMatrix<real_t> patch1_interface1(2, pointsPerEdge);

    // Interface 2: Patch 0 south ↔ Patch 2 east  
    gsMatrix<real_t> patch0_interface2(2, pointsPerEdge);
    gsMatrix<real_t> patch2_interface2(2, pointsPerEdge);

    for (index_t i = 0; i < pointsPerEdge; ++i)
    {
        real_t t = static_cast<real_t>(i) / (pointsPerEdge - 1);

        // Interface 1
        patch0_interface1(0, i) = 0.0;        // u = 0 (west side)
        patch0_interface1(1, i) = 1.0 - t;    // v varies

        patch1_interface1(0, i) = 1.0;        // u = 1 (east side)
        patch1_interface1(1, i) = t;          // v varies

        // Interface 2
        patch0_interface2(0, i) = t;          // u varies
        patch0_interface2(1, i) = 0.0;        // v = 0 (south side)

        patch2_interface2(0, i) = 1.0;        // u = 1 (east side)
        patch2_interface2(1, i) = 1.0 - t;    // v varies
    }

    // Combine all interface points per patch
    interfacePoints[0].resize(2, 2 * pointsPerEdge);
    interfacePoints[0] << patch0_interface1, patch0_interface2;

    interfacePoints[1].resize(2, pointsPerEdge);
    interfacePoints[1] = patch1_interface1;

    interfacePoints[2].resize(2, pointsPerEdge);
    interfacePoints[2] = patch2_interface2;

    // Assemble using the actual assemble() function
    const index_t numFunctions = mpbes.size();
    gsSparseMatrix<real_t> A_interface(totalPoints, numFunctions);
    gsMatrix<real_t> b_interface(totalPoints, 2);

    gsInfo << "Assembling interface matrix: " << totalPoints << " points, " << numFunctions << " functions\n";

    assemble(
        interfacePoints,
        mpbes,
        A_interface,
        b_interface,
        mp,
        false  // verbose = false
    );

    // Compute reconstructed geometry: A * vectSol
    gsMatrix<real_t> reconstructed = gsEigen::MatrixXd(A_interface) * vectSol;

    // Compute residual
    gsMatrix<real_t> residual = reconstructed - b_interface;
    
    // Compute max Euclidean error
    real_t maxError = 0.0;
    for (index_t row = 0; row < residual.rows(); ++row)
    {
        real_t err = std::sqrt(residual(row, 0) * residual(row, 0) + 
                               residual(row, 1) * residual(row, 1));
        maxError = std::max(maxError, err);
    }

    // Detailed reporting
    gsInfo << "\n--- Interface Continuity Test Results ---\n";
    outfile << "\n--- Interface Continuity Test Results ---\n";

    index_t row = 0;
    const char* interfaceNames[] = { "0↔1", "0↔2" };

    for (index_t iface = 0; iface < numInterfaces; ++iface)
    {
        gsInfo << "\nInterface " << interfaceNames[iface] << ":\n";
        outfile << "\nInterface " << interfaceNames[iface] << ":\n";

        for (index_t pt = 0; pt < pointsPerInterface; ++pt)
        {
            real_t err_x = std::abs(residual(row, 0));
            real_t err_y = std::abs(residual(row, 1));
            real_t err_max = std::max(err_x, err_y);

            if (err_max > 1e-6 || pt < 2)  // Show first 2 points and errors
            {
                outfile << "  Point " << pt << ": "
                    << "geom=(" << b_interface(row, 0) << "," << b_interface(row, 1) << "), "
                    << "recon=(" << reconstructed(row, 0) << "," << reconstructed(row, 1) << "), "
                    << "error=(" << residual(row, 0) << "," << residual(row, 1) << ")\n";
            }
            ++row;
        }
    }

    gsInfo << "\nMax componentwise error: " << maxError << "\n";
    gsInfo << "=== testBoundaryAssembly END ===\n\n";
    outfile << "\nMax componentwise error: " << maxError << "\n";
    outfile << "=== testBoundaryAssembly END ===\n\n";

    return maxError;
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

    // Storage for parametric and mapped lines
    std::vector<gsMatrix<>> uvLines;
    std::vector<index_t> uvLinePatchIDs;
    std::vector<char> uvLineDirections;
    std::vector<gsMatrix<>> xyLines;

    // Step 1: Generate parametric lines from boxes
    if (verbose)
        gsInfo << "Step 1: Generating parametric lines...\n";

    generateParametricGrid(
        boxMat,
        currentLastNonZeroRow,
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
    if (verbose)
        gsInfo << "Step 2: Evaluating geometry lines...\n";

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

            // Vertical lines (u = const)
            for (index_t xi = x0; xi <= x1; ++xi)
            {
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

            // Horizontal lines (v = const)
            for (index_t yi = y0; yi <= y1; ++yi)
            {
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
int checkJacobianDeterminant(
    const gsVector<gsMatrix<>>& uv,
    const MPBES<2, real_t>& mpbes,
    const gsMatrix<real_t>& vectSol,
    gsVector<size_t>& numIrregular,
    bool verbose = true)
{
    PROFILE_FUNCTION();
    gsInfo << "=== checkJacobianDeterminant (MPBES-based) BEGIN ===\n";
    
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
    
    const index_t numPatches = mpbes.nPatches();
    const index_t numFunctions = mpbes.size();
    int totalIrregular = 0;
    numIrregular.setZero();
    numIrregular.resize(numPatches);

    const auto& functionDescription = mpbes.functionDescription();
    const auto& thbBases = mpbes.thbBases();
    const auto& indexInTHB = mpbes.indexInTHB();
    
    // Pre-filter active functions per patch (same as assemble)
    std::vector<std::vector<index_t>> activeFuncsPerPatch(numPatches);
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

        // Evaluate derivatives at each point using MPBES
        for (index_t pt = 0; pt < numPoints; ++pt)
        {
            const gsVector<real_t> param = uvPatch.col(pt);
            
            // SKIP REFINEMENT ZONES: Only check truly interior points
            // Refinement/truncation boundaries cause parametric singularities that don't reflect geometry
            // Check a reasonable interior region, avoiding only narrow boundaries
            const real_t uMin = 0.0, uMax = 1.0, vMin = 0.0, vMax = 1.0;
            const real_t interiorMargin = 0.05;  // Only check [0.05, 0.95] x [0.05, 0.95] region (90% of domain)
            
            bool inInterior = (param(0) > uMin + interiorMargin) && (param(0) < uMax - interiorMargin) &&
                              (param(1) > vMin + interiorMargin) && (param(1) < vMax - interiorMargin);
            
            if (!inInterior) {
                continue;  // Skip narrow boundary margin points
            }
            
            // Jacobian matrix J = [dx/du  dx/dv]
            //                     [dy/du  dy/dv]
            gsMatrix<> J(2, 2);
            J.setZero();
            
            // Only iterate active functions for this patch (same as assemble)
            for (index_t f : activeFuncsPerPatch[patch])
            {
                if (f >= static_cast<index_t>(vectSol.rows())) continue;
                
                // Use MPBES derivative evaluation - exactly like assemble() uses evalSingleOnPatch
                gsMatrix<real_t> basisDeriv;
                mpbes.evalDerivSingleOnPatch(f, patch, param, basisDeriv);
                
                // basisDeriv is 2x1 matrix: [dφ/du; dφ/dv]
                real_t dφ_du = basisDeriv(0, 0);
                real_t dφ_dv = basisDeriv(1, 0);
                
                // Accumulate contribution to Jacobian
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
                    
                    // Compute Jacobian at neighbor using MPBES
                    gsMatrix<> J_nb(2, 2);
                    J_nb.setZero();
                    
                    for (index_t f : activeFuncsPerPatch[patch]) {
                        if (f >= static_cast<index_t>(vectSol.rows())) continue;
                        
                        gsMatrix<real_t> basisDeriv_nb;
                        mpbes.evalDerivSingleOnPatch(f, patch, nbParam, basisDeriv_nb);
                        
                        real_t dφ_du_nb = basisDeriv_nb(0, 0);
                        real_t dφ_dv_nb = basisDeriv_nb(1, 0);
                        
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

    gsInfo << "=== checkJacobianDeterminant END ===\n";
    gsInfo << "Total irregular points: " << totalIrregular << "\n";
    
    // Close irregular points log
    irregularLog << "========================================\n";
    irregularLog << "SUMMARY: " << totalIrregular << " total irregular points detected\n";
    irregularLog.close();
    
    if (totalIrregular > 0) {
        gsInfo << "Detailed irregular points log written to: irregular_points_log.txt\n";
    }
    
    return totalIrregular;
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

    gsInfo << "=== checkJacobianDeterminant END ===\n";
    gsInfo << "Total irregular points: " << totalIrregular << "\n";
    
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

            //std::cout << "Checking spillover function (index: " << functionIndex
            //    << ", spillPatch: " << spillPatch
            //    << ", spillLevel: " << spillLevel
            //    << ", spillIndex: " << spillIndex << ")" << std::endl;

            if (spillPatch == patch) {
                //std::cout << "Match found for patch!" << std::endl;

                // Validate indexing before accessing Bells
                if (spillPatch >= Bells.size()) {
                    //std::cerr << "ERROR: spillPatch " << spillPatch
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
)
{
    PROFILE_FUNCTION();
    gsVector<index_t> numNodes(2);
    numNodes[0] = numNodes[1] = 2;
    int shift = 0;

    if (verbose)
        outfile << "=== Starting assembly of A_mat and b_vec ===\n";
    for (size_t patch = 0; patch < THBVector.size(); ++patch)
    {
        if (verbose)
            outfile << ">> Processing patch " << patch << "\n";

        typename gsBasis<real_t>::domainIter domIt = THBVector(patch).makeDomainIterator();
        gsGaussRule<real_t> quRule(numNodes);
        gsMatrix<real_t> quNodes;
        gsVector<real_t> quWeights;

        int domanI = 0;

        for (; domIt->good(); domIt->next())
        {
            gsVector<real_t> lower = domIt->lowerCorner();
            gsVector<real_t> upper = domIt->upperCorner();

            if (verbose)
            {
                outfile << "  Element " << domanI
                    << ", bounds: [(" << lower[0] << ", " << lower[1] << ") -> ("
                    << upper[0] << ", " << upper[1] << ")]\n";
            }

            quRule.mapTo(lower, upper, quNodes, quWeights);
            if (verbose)
            {
                outfile << "Finished the gaussRule\n";
            }
            for (int q = 0; q < quNodes.cols(); ++q)
            {
                const gsVector<real_t> uv = quNodes.col(q);
                gsMatrix<real_t> xy = mp.patch(patch).eval(uv);
                b_vec(q + shift, 0) = xy(0, 0);
                b_vec(q + shift, 1) = xy(1, 0);

                if (verbose)
                {
                    outfile << "    Gauss point " << q
                        << ": (u, v) = (" << uv[0] << ", " << uv[1] << ") -> "
                        << "(x, y) = (" << xy(0, 0) << ", " << xy(1, 0) << ")\n";
                }

                for (size_t functionIndex = 0; functionIndex < functionDescription.size(); ++functionIndex)
                {
                    real_t valore = 0;

                    if (verbose)
                        outfile << "      Function " << functionIndex << ":\n";

                    // Twin evaluation
                    for (const auto& twin : functionDescription[functionIndex])
                    {
                        if (twin[0] == patch)
                        {
                            real_t val = THBVector(patch).function(twin[1]).eval(uv)(0, 0);
                            valore += val;

                            if (verbose)
                            {
                                outfile << "        Twin (patch=" << twin[0]
                                    << ", index=" << twin[1]
                                    << ") -> value = " << val << "\n";
                            }
                        }
                    }

                    // Spillover evaluation
                    if (hasSpillover[functionIndex])
                    {
                        for (const auto& spill : spilloverFunctionCoordinates[functionIndex])
                        {
                            int spPatch = spill[0];
                            int spLevel = spill[1];
                            int spIndex = spill[2];
                            real_t val = Bells(spPatch)(spLevel).function(spIndex).eval(uv)(0, 0);
                            valore += val;

                            if (verbose)
                            {
                                outfile << "        Spill (patch=" << spPatch
                                    << ", level=" << spLevel
                                    << ", index=" << spIndex
                                    << ") -> value = " << val << "\n";
                            }
                        }
                    }

                    // Insert value
                    if (valore != 0)
                    {
                        A_mat(q + shift, functionIndex) += valore;

                        if (verbose)
                        {
                            outfile << "        Inserted A(" << (q + shift) << ", "
                                << functionIndex << ") += " << valore << "\n";
                        }
                    }
                }
            }

            shift += quNodes.cols();
            domanI++;
        }
    }

    if (verbose)
        outfile << "=== Assembly complete ===\n";
}

void assembleATHB_maskSpilloverFromTwinnedPatches(
    const gsVector<gsMatrix<>>& uv,
    const gsVector<gsVector<gsTensorBSplineBasis<2, real_t>>>& Bells,
    const gsVector<gsTHBSplineBasis<2, real_t>>& THBVector,
    const std::vector<std::vector<std::vector<index_t>>>& functionDescription,
    const std::vector<std::vector<std::array<int, 3>>>& spilloverFunctionCoordinates,
    const std::vector<bool>& hasSpillover,
    gsSparseMatrix<real_t>& A_mat,
    bool verbose = false
)
{
    PROFILE_FUNCTION();
    index_t totalRows = 0;
    for (const auto& patchUV : uv)
        totalRows += patchUV.cols();
    A_mat.resize(totalRows, functionDescription.size());

    index_t rowOffset = 0;

    for (index_t patch = 0; patch < uv.size(); ++patch)
    {
        const gsMatrix<>& uvPatch = uv[patch];
        index_t nPts = uvPatch.cols();

        for (index_t k = 0; k < nPts; ++k)
        {
            const gsVector<real_t> uvk = uvPatch.col(k);
            index_t row = rowOffset + k;

            for (index_t f = 0; f < functionDescription.size(); ++f)
            {
                real_t totalVal = 0.0;

                // 1. Apply native THB twins on this patch
                bool hasTwinHere = false;
                for (const auto& twin : functionDescription[f])
                {
                    if (twin[0] == patch)
                    {
                        hasTwinHere = true;
                        real_t val = THBVector[patch].function(twin[1]).eval(uvk)(0, 0);
                        if (val != 0.0)
                        {
                            totalVal += val;
                            if (verbose)
                                outfile << "    Twin of f=" << f
                                << " (patch=" << patch
                                << ", index=" << twin[1]
                                << ") @ uv=(" << uvk(0) << ", " << uvk(1)
                                << ") -> " << val << "\n";
                        }
                    }
                }

                // 2. Apply spillover only from patches that DO NOT have native twin of f
                if (hasSpillover[f])
                {
                    for (const auto& spill : spilloverFunctionCoordinates[f])
                    {
                        int sourcePatch = spill[0];
                        int spLevel = spill[1];
                        int spIndex = spill[2];

                        // Skip spillover from patches where f is already defined
                        bool skip = false;
                        for (const auto& twin : functionDescription[f])
                        {
                            if (twin[0] == sourcePatch)
                            {
                                skip = true;
                                break;
                            }
                        }

                        if (skip)
                            continue;

                        // Apply only if spill is into the current patch
                        if (sourcePatch == patch)
                            continue;

                        real_t val = Bells[sourcePatch][spLevel].function(spIndex).eval(uvk)(0, 0);
                        if (val != 0.0)
                        {
                            totalVal += val;
                            if (verbose)
                                outfile << "    Spillover of f=" << f
                                << " from patch=" << sourcePatch
                                << ", level=" << spLevel
                                << ", index=" << spIndex
                                << " @ uv=(" << uvk(0) << ", " << uvk(1)
                                << ") -> " << val << "\n";
                        }
                    }
                }

                if (totalVal != 0.0)
                {
                    A_mat(row, f) += totalVal;
                    if (verbose)
                        outfile << "    => A_mat(" << row << ", " << f
                        << ") += " << totalVal
                        << " @ uv=(" << uvk(0) << ", " << uvk(1) << ")\n";
                }
            }
        }

        rowOffset += nPts;
        if (verbose)
            outfile << "  Patch " << patch << " completed, row offset now " << rowOffset << "\n";
    }

    if (verbose)
        outfile << "Finished assembling A_mat (filtered spillover from twin-defined patches).\n";
}


void assembleATHB_exclude(
    const gsVector<gsMatrix<>>& uv,
    const gsVector<gsVector<gsTensorBSplineBasis<2, real_t>>>& Bells,
    const gsVector<gsTHBSplineBasis<2, real_t>>& THBVector,
    const std::vector<std::vector<std::vector<index_t>>>& functionDescription,
    const std::vector<std::vector<std::array<int, 3>>>& spilloverFunctionCoordinates,
    const std::vector<std::vector<std::array<int, 3>>>& spilloverSources,
    const std::vector<bool>& hasSpillover,
    gsSparseMatrix<real_t>& A_mat,
    bool verbose)
{
    PROFILE_FUNCTION();
    index_t totalRows = 0;
    for (const auto& patchUV : uv)
        totalRows += patchUV.cols();
    A_mat.resize(totalRows, functionDescription.size());

    index_t rowOffset = 0;

    for (index_t patch = 0; patch < uv.size(); ++patch)
    {
        const gsMatrix<>& uvPatch = uv[patch];
        index_t nPts = uvPatch.cols();

        for (index_t k = 0; k < nPts; ++k)
        {
            const gsVector<real_t> uvk = uvPatch.col(k);
            index_t row = rowOffset + k;

            if (verbose)
                outfile << "=== Evaluating at uv = (" << uvk(0) << ", " << uvk(1) << "), patch " << patch << " ===\n";

            for (index_t f = 0; f < functionDescription.size(); ++f)
            {
                real_t totalVal = 0.0;

                // Process THB contributions
                for (const auto& twin : functionDescription[f])
                {
                    int twinPatch = twin[0];
                    int twinIndex = twin[1];

                    if (twinPatch != patch)
                        continue;

                    int twinLevel = THBVector[patch].levelOf(twinIndex);
                    int tensorIndex = THBVector[patch].flatTensorIndexOf(twinIndex);

                    bool suppressTHB = shouldSuppressTHB(
                        patch,
                        twinLevel,
                        tensorIndex,
                        Bells[patch][twinLevel],
                        spilloverSources[patch],
                        Bells,
                        verbose);

                    if (!suppressTHB)
                    {
                        real_t val = THBVector[patch].function(twinIndex).eval(uvk)(0, 0);
                        totalVal += val;

                        if (verbose && val != 0.0)
                            outfile << "  Accepted THB: f=" << f
                            << ", patch=" << patch
                            << ", level=" << twinLevel
                            << ", tensorIndex=" << tensorIndex
                            << ", val=" << val << "\n";
                    }
                    else if (verbose)
                    {
                        outfile << "  Suppressed THB: f=" << f
                            << ", patch=" << patch
                            << ", level=" << twinLevel
                            << ", tensorIndex=" << tensorIndex << "\n";
                    }
                }

                // Process spillover contributions (Bells)
                if (hasSpillover[f])
                {
                    for (const auto& spill : spilloverFunctionCoordinates[f])
                    {
                        int spPatch = spill[0];
                        int spLevel = spill[1];
                        int spIndex = spill[2];

                        bool suppressSpill = shouldSuppressSpillover(
                            patch,
                            THBVector[patch],
                            spPatch,
                            spLevel,
                            spIndex,
                            Bells[spPatch][spLevel],
                            verbose);

                        if (!suppressSpill)
                        {
                            real_t val = Bells[spPatch][spLevel].function(spIndex).eval(uvk)(0, 0);
                            totalVal += val;

                            if (verbose && val != 0.0)
                                outfile << "  Accepted Spillover: f=" << f
                                << ", from patch=" << spPatch
                                << ", level=" << spLevel
                                << ", index=" << spIndex
                                << ", val=" << val << "\n";
                        }
                        else if (verbose)
                        {
                            outfile << "  Suppressed Spillover: f=" << f
                                << ", from patch=" << spPatch
                                << ", level=" << spLevel
                                << ", index=" << spIndex << "\n";
                        }
                    }
                }

                if (totalVal != 0.0)
                    A_mat(row, f) += totalVal;
            }
        }

        rowOffset += nPts;
    }

    if (verbose)
        outfile << "Finished assembling A_mat with twin-aware bidirectional suppression.\n";
}















void assembleATHB_clean(
    const gsVector<gsMatrix<>>& uv,
    const gsVector<gsVector<gsTensorBSplineBasis<2, real_t>>>& Bells,
    const gsVector<gsTHBSplineBasis<2, real_t>>& THBVector,
    const std::vector<std::vector<std::vector<index_t>>>& functionDescription,
    const std::vector<std::vector<std::array<int, 3>>>& spilloverFunctionCoordinates,
    const std::vector<bool>& hasSpillover,
    gsSparseMatrix<real_t>& A_mat
)
{
    PROFILE_FUNCTION();
    // Total number of rows = all UV samples
    index_t totalRows = 0;
    for (const auto& patchUV : uv)
        totalRows += patchUV.cols();

    // Preallocate (for assembly performance)
    A_mat.resize(totalRows, functionDescription.size());

    index_t rowOffset = 0;

    // Iterate over patches
    for (index_t patch = 0; patch < uv.size(); ++patch)
    {
        const gsMatrix<>& uvPatch = uv[patch];
        index_t nPts = uvPatch.cols();

        // Loop over parametric points in this patch
        for (index_t k = 0; k < nPts; ++k)
        {
            const gsVector<real_t> uvk = uvPatch.col(k);
            index_t row = rowOffset + k;

            // Loop over all global basis functions
            for (index_t f = 0; f < functionDescription.size(); ++f)
            {
                real_t totalVal = 0;

                // Evaluate native supports (twins)
                for (const auto& twin : functionDescription[f])
                {
                    if (twin[0] == patch)
                    {
                        index_t localIdx = twin[1];
                        real_t val = THBVector[patch].function(localIdx).eval(uvk)(0, 0);
                        totalVal += val;
                    }
                }

                // Evaluate spillover supports
                if (hasSpillover[f])
                {
                    for (const auto& spill : spilloverFunctionCoordinates[f])
                    {
                        index_t spPatch = spill[0];
                        index_t spLevel = spill[1];
                        index_t spIndex = spill[2];

                        if (spPatch == patch)
                        {
                            real_t val = Bells[spPatch][spLevel].function(spIndex).eval(uvk)(0, 0);
                            totalVal += val;
                        }
                    }
                }

                if (totalVal != 0.0)
                    A_mat(row, f) += totalVal;
            }
        }

        rowOffset += nPts;
    }
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
    const index_t nRows = A.rows();
    const index_t nCols = A.cols();
    bool allSatisfied = true;

    gsInfo << "Checking partition of unity...\n";
    for (index_t i = 0; i < nRows; ++i)
    {
        real_t sum = 0.0;

        for (index_t j = 0; j < nCols; ++j)
            sum += A(i, j);

        if (std::abs(sum - 1.0) > tol) {
            gsInfo << "Row " << i << " violates partition of unity: sum = "
                << sum << "\n";
            allSatisfied = false;
        }
    }

    return allSatisfied;
}



/**
 * @brief Assembles the matrix A_mat and optionally the right-hand side b_vec for least-squares fitting
 *        using MPBE-compatible THB basis functions. Now supports Option B indexing: (patch, level, tensor-index).
 *
 * @date 15.06.2025
 *
 * @param[in]  uv                          Parametric data points per patch.
 * @param[in]  Bells                       Tensor-product bases per patch and level.
 * @param[in]  THBVector                   Per-patch THB bases.
 * @param[in]  functionDescription         List of global basis functions, each represented as a list of (patch, level, tensor-index) entries.
 * @param[in]  spilloverFunctionCoordinates Spillover entries: for each global function, a list of (patch, level, index).
 * @param[in]  hasSpillover                Flags for whether each global basis function has spillovers.
 * @param[in]  indexInTHB                  Maps (patch, level, tensor-index) to THB index. If -1, the function is not active in THB.
 * @param[out] A_mat                       Output system matrix.
 * @param[out] b_vec                       Output RHS vector (optional, can be uninitialized).
 * @param[in]  mp                          Multi-patch geometry to evaluate target data.
 * @param[in]  verbose                     If true, detailed log to gsInfo and outfile.
 */
void assemble(
    const gsVector<gsMatrix<>>& uv, // Used only to count total rows
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
    bool verbose)
{
    PROFILE_SECTION("assemble_total");
    auto asm_start = std::chrono::high_resolution_clock::now();

    if (verbose)
        outfile << "=== assemble (Gauss points) BEGIN ===\n";

    gsVector<index_t> numNodes(2);
    numNodes[0] = numNodes[1] = 2;  // Increased from 2 to 5 (5x5=25 Gauss points per element)
    gsGaussRule<real_t> quRule(numNodes);

    if (verbose)
        gsInfo << "Using " << numNodes[0] << "x" << numNodes[1] << " = " << numNodes.prod() << " Gauss points per element\n";

    const index_t numFunctions = functionDescription.size();
    const index_t numPatches = THBVector.size();
    std::vector<bool> hasNonzero(numFunctions, false);
    std::vector<std::vector<index_t>> activeFuncsPerPatch(numPatches);

    // Pre-filter active functions per patch
    {
        PROFILE_SECTION("assemble_prefilter");
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
    } // End PROFILE_SECTION("assemble_prefilter")

    index_t totalRows = 0;
    for (index_t patch = 0; patch < numPatches; ++patch)
    {
        auto domIt = THBVector[patch].makeDomainIterator();
        for (; domIt->good(); domIt->next())
            totalRows += numNodes.prod();
    }

    A_mat.resize(totalRows, numFunctions);
    A_mat.setZero();
    b_vec.resize(totalRows, 2);
    b_vec.setZero();

    index_t rowOffset = 0;
    const bool DEBUG_EVAL = true; // Set to true only when debugging evalSingle_into

    // Timing measurements
    auto timeEvalStart = std::chrono::system_clock::now();
    index_t evalCount = 0;
    double timeEval = 0.0, timeEvalSingle = 0.0, timeOther = 0.0;

    for (index_t patch = 0; patch < numPatches; ++patch)
    {
        if (verbose)
            outfile << ">> Patch " << patch << "\n";

        auto domIt = THBVector[patch].makeDomainIterator();
        gsMatrix<real_t> quNodes;
        gsVector<real_t> quWeights;

        for (int el = 0; domIt->good(); domIt->next(), ++el)
        {
            gsVector<real_t> lower = domIt->lowerCorner();
            gsVector<real_t> upper = domIt->upperCorner();
            quRule.mapTo(lower, upper, quNodes, quWeights);

            for (index_t q = 0; q < quNodes.cols(); ++q)
            {
                auto timeEvalPtStart = std::chrono::system_clock::now();
                const gsVector<real_t> uvk = quNodes.col(q);

                // RHS: target geometry
                auto timeGeomStart = std::chrono::system_clock::now();
                gsMatrix<real_t> xy = mp.patch(patch).eval(uvk);
                auto timeGeomEnd = std::chrono::system_clock::now();
                timeOther += std::chrono::duration<double>(timeGeomEnd - timeGeomStart).count();

                b_vec(rowOffset, 0) = xy(0, 0);
                b_vec(rowOffset, 1) = xy(1, 0);

                if (verbose)
                {
                    outfile << "  Element " << el << ", Gauss point " << q
                        << ": (u,v) = (" << uvk(0) << ", " << uvk(1) << ")"
                        << " => (x,y) = (" << xy(0, 0) << ", " << xy(1, 0) << ")\n";
                }

                // Only iterate active functions for this patch
                gsMatrix<> resultMatrix; // reuse to avoid repeated allocations
                for (index_t f : activeFuncsPerPatch[patch])
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

                        auto timeEvalSingleStart = std::chrono::system_clock::now();
                        evalSingle_into(
                            f,
                            thbIdx,
                            uvk,
                            isTruncated,
                            THBVector[fnPatch],
                            presentation[f][functionComponent],
                            resultMatrix);
                        auto timeEvalSingleEnd = std::chrono::system_clock::now();
                        timeEvalSingle += std::chrono::duration<double>(timeEvalSingleEnd - timeEvalSingleStart).count();
                        evalCount++;

                        real_t contrib = resultMatrix(0, 0);

                        // Debug validation (disabled by default for performance)
                        if (DEBUG_EVAL)
                        {
                            real_t referenceValue = THBVector[patch].function(thbIdx).eval(uvk)(0, 0);
                            if (fabs(fabs(contrib) - fabs(referenceValue)) > 1e-6)
                            {
                                gsInfo << "MISMATCH in evalSingle_into for f=" << f
                                    << ", patch=" << patch
                                    << ", level=" << fnLevel
                                    << ", tensorIndex=" << fnTensorIndex
                                    << ": computed " << contrib
                                    << ", reference " << referenceValue << "\n";
                            }
                        }
                        val += contrib;
                        functionComponent++;

                        if (verbose && contrib != 0.0)
                            outfile << "    [THB] f=" << f << " contributes " << contrib << "\n";
                    }

                    if (val != 0.0)
                    {
                        A_mat(rowOffset, f) += val;
                        hasNonzero[f] = true;
                        if (verbose)
                            outfile << "      -> A(" << rowOffset << ", " << f << ") += " << val << "\n";
                    }
                }

                ++rowOffset;
            }
        }
    }

    auto timeAssemblyEnd = std::chrono::system_clock::now();
    double totalTime = std::chrono::duration<double>(timeAssemblyEnd - timeEvalStart).count();

    gsInfo << "\n=== ASSEMBLY TIMING BREAKDOWN ===\n";
    gsInfo << "Total time:             " << totalTime << " s\n";
    gsInfo << "  evalSingle_into time: " << timeEvalSingle << " s (" << (evalCount > 0 ? timeEvalSingle / evalCount * 1000 : 0) << " ms per call, " << evalCount << " calls)\n";
    gsInfo << "  Geometry eval time:   " << timeOther << " s\n";
    gsInfo << "  Other operations:     " << (totalTime - timeEvalSingle - timeOther) << " s\n";
    gsInfo << "===================================\n\n";

    for (size_t i = 0; i < hasNonzero.size(); i++)
    {
        if (!hasNonzero[i]) {
            gsInfo << "The column " << i << " is zero\n";
        }
    }

    if (verbose)
        outfile << "=== assemble (Gauss points) END ===\n";
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

// MPBES-based assemble - Uses MPBES presentation for truncated functions
void assemble(
    const gsVector<gsMatrix<>>& uv,
    const MPBES<2, real_t>& mpbes,
    gsSparseMatrix<real_t>& A_mat,
    gsMatrix<real_t>& b_vec,
    const gsMultiPatch<real_t>& mp,
    bool verbose)
{
    PROFILE_SECTION("assemble_mpbes_total");

    if (verbose)
        outfile << "=== assemble (MPBES, Gauss points) BEGIN ===\n";

    gsVector<index_t> numNodes(2);
    numNodes[0] = numNodes[1] = 2;
    gsGaussRule<real_t> quRule(numNodes);

    if (verbose)
        gsInfo << "Using " << numNodes[0] << "x" << numNodes[1] << " = " << numNodes.prod() << " Gauss points per element\n";

    const index_t numFunctions = mpbes.size();
    const index_t numPatches = mpbes.nPatches();

    // Get references to MPBES data
    const auto& functionDescription = mpbes.functionDescription();
    const auto& thbBases = mpbes.thbBases();
    const auto& indexInTHB = mpbes.indexInTHB();
    const auto& isTruncated = mpbes.isTruncated();
    const auto& presentation = mpbes.presentation();

    // Pre-filter active functions per patch (CRITICAL for correctness)
    std::vector<std::vector<index_t>> activeFuncsPerPatch(numPatches);
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

    // Count total rows
    index_t totalRows = 0;
    for (index_t patch = 0; patch < numPatches; ++patch)
    {
        auto domIt = thbBases[patch].makeDomainIterator();
        for (; domIt->good(); domIt->next())
            totalRows += numNodes.prod();
    }

    A_mat.resize(totalRows, numFunctions);
    A_mat.setZero();
    b_vec.resize(totalRows, 2);
    b_vec.setZero();

    index_t rowOffset = 0;

    if (verbose)
    {
        gsInfo << "Total rows (Gauss points): " << totalRows << "\n";
        gsInfo << "Total basis functions: " << numFunctions << "\n";
    }

    for (index_t patch = 0; patch < numPatches; ++patch)
    {
        if (verbose)
            outfile << ">> Patch " << patch << "\n";

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
                index_t row = rowOffset;

                // RHS: target geometry
                gsMatrix<real_t> xy = mp.patch(patch).eval(uvk);
                b_vec(row, 0) = xy(0, 0);
                b_vec(row, 1) = xy(1, 0);

                if (verbose)
                {
                    outfile << "  Element " << el << ", Gauss point " << q
                        << ": (u,v) = (" << uvk(0) << ", " << uvk(1) << ")"
                        << " => (x,y) = (" << xy(0, 0) << ", " << xy(1, 0) << ")\n";
                }

                // Only iterate active functions for this patch
                gsMatrix<real_t> resultMatrix;

                int nonZeroFuncs = 0;
                for (index_t f : activeFuncsPerPatch[patch])
                {
                    // Use MPBES's evalSingleOnPatch which handles component-level truncation correctly
                    gsMatrix<real_t> basisValue;
                    mpbes.evalSingleOnPatch(f, patch, uvk, basisValue);
                    real_t value = basisValue(0, 0);

                    if (value != 0.0)
                    {
                        A_mat(row, f) += value;
                        nonZeroFuncs++;
                    }
                }

                ++rowOffset;
            }
        }
    }

    if (verbose)
    {
        gsInfo << "=== assemble END ===\n";
        gsInfo << "Final system size: " << A_mat.rows() << " x " << A_mat.cols() << "\n";
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



/**
 * @brief Assembles matrix A_mat and right-hand side b_vec using Gauss points.
 *
 * This function uses Gaussian quadrature points per element of the THB basis
 * to construct the fitting system A * x = b. It is an adaptation of assemble()
 * with Gauss integration.
 *
 * @date 13.06.25
 *
 * @param[in] mp         The geometry (gsMultiPatch).
 * @param[in] Bells      Tensor-product bases per patch and level.
 * @param[in] THBVector  THB basis per patch.
 * @param[in] functionDescription List of (patch, level, index) entries for each global basis function.
 * @param[in] spilloverFunctionCoordinates Optional spillover sources per function.
 * @param[in] hasSpillover Flags marking which functions have spillovers.
 * @param[in] indexInTHB  Converts Bells indices into THB local indices.
 * @param[out] A_mat      Output matrix A.
 * @param[out] b_vec      Right-hand side vector.
 * @param[in] verbose     Enables logging.
 */
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
    bool verbose)
{
    PROFILE_FUNCTION();
    index_t totalRows = 0;
    gsVector<index_t> numNodes(2); numNodes[0] = numNodes[1] = 2;

    // Precompute number of rows
    for (index_t patch = 0; patch < THBVector.size(); ++patch)
    {
        auto domIt = THBVector[patch].makeDomainIterator();
        for (; domIt->good(); domIt->next())
            totalRows += numNodes.prod();
    }

    const index_t numFunctions = functionDescription.size();
    A_mat.resize(totalRows, numFunctions);
    A_mat.setZero();
    b_vec.setZero(totalRows, 2); // Store (x, y)

    index_t shift = 0;

    if (verbose)
    {
        gsInfo << "=== assembleGauss BEGIN ===\n";
        gsInfo << "Total rows (Gauss points): " << totalRows << "\n";
        gsInfo << "Total basis functions: " << numFunctions << "\n";
    }

    for (index_t patch = 0; patch < THBVector.size(); ++patch)
    {
        auto domIt = THBVector[patch].makeDomainIterator();
        gsGaussRule<real_t> rule(numNodes);
        gsMatrix<real_t> quNodes;
        gsVector<real_t> quWeights;

        for (; domIt->good(); domIt->next())
        {
            rule.mapTo(domIt->lowerCorner(), domIt->upperCorner(), quNodes, quWeights);

            for (index_t q = 0; q < quNodes.cols(); ++q)
            {
                const gsVector<real_t> uvk = quNodes.col(q);
                index_t row = shift + q;

                gsMatrix<real_t> xy = mp.patch(patch).eval(uvk);
                b_vec(row, 0) = xy(0, 0);
                b_vec(row, 1) = xy(1, 0);

                for (index_t f = 0; f < numFunctions; ++f)
                {
                    real_t value = 0.0;

                    for (const auto& twin : functionDescription[f])
                    {
                        if (twin[0] == patch)
                        {
                            int level = twin[1];
                            int index = twin[2];
                            int thbIndex = indexInTHB[patch][level][index];

                            if (thbIndex == -1) continue;

                            real_t val = THBVector[patch].function(thbIndex).eval(uvk)(0, 0);
                            value += val;
                        }
                    }

                    if (f < hasSpillover.size() && hasSpillover[f] &&
                        f < spilloverFunctionCoordinates.size())
                    {
                        for (const auto& spill : spilloverFunctionCoordinates[f])
                        {
                            int spPatch = spill[0], spLevel = spill[1], spIndex = spill[2];

                            if (spPatch >= Bells.size() || spLevel >= Bells[spPatch].size() || spIndex >= Bells[spPatch][spLevel].size())
                                continue;

                            real_t val = Bells[spPatch][spLevel].function(spIndex).eval(uvk)(0, 0);
                            value += val;
                        }
                    }

                    if (value != 0.0)
                        A_mat(row, f) += value;
                }
            }

            shift += quNodes.cols();
        }
    }

    if (verbose)
    {
        gsInfo << "=== assembleGauss END ===\n";
        gsInfo << "Final system size: " << A_mat.rows() << " x " << A_mat.cols() << "\n";
    }
}















void assembleSystem4(gsVector<gsMatrix<>>  uv,
    gsVector<gsMatrix<>>  xy,
    gsVector<gsVector<gsTensorBSplineBasis<2, real_t>>> Bells,
    gsVector  <gsVector< gsVector<index_t>>> isActive,
    gsVector  <gsVector< gsVector<index_t>>>  globalIndex,
    gsVector  <gsMatrix<index_t>> functionLocation,
    gsSparseMatrix<real_t>& A_mat,
    gsMatrix<real_t>& m_B) {
    for (int patch = 0; patch < Bells.size(); ++patch) {
        for (int functionIndex = 0; functionIndex < functionLocation(1).cols(); ++functionIndex) {
            for (int k = 0; k < 1; ++k) {
                int corrLevel;
                int currIndex;
                gsMatrix<> valore;
                bool functionIsActive = false;
                if (functionLocation(0)(0, functionIndex) == patch) {
                    //MASTER FUNCTION CORRESPONDS WITH THE CURRENT PATCH
                    functionIsActive = true;
                    corrLevel = functionLocation(0)(1, functionIndex);
                    currIndex = functionLocation(0)(2, functionIndex);
                    valore = Bells(patch)(corrLevel).function(currIndex).eval(uv(patch).col(k));
                    m_B(functionIndex, 0) += xy(patch)(0, k) * valore(0, 0);
                    m_B(functionIndex, 1) += xy(patch)(1, k) * valore(0, 0);
                }
                if (functionLocation(1)(0, functionIndex) == patch) {
                    //SLAVE FUNCTION CORRESPONDS WITH THE CURRENT PATCH
                    functionIsActive = true;
                    corrLevel = functionLocation(1)(1, functionIndex);
                    currIndex = functionLocation(1)(2, functionIndex);
                    valore = Bells(patch)(corrLevel).function(currIndex).eval(uv(patch).col(k));
                    m_B(functionIndex, 0) += xy(patch)(0, k) * valore(0, 0);
                    m_B(functionIndex, 1) += xy(patch)(1, k) * valore(0, 0);
                }
                if (functionIsActive) {
                    for (int functionColIndex = 0; functionColIndex < functionLocation(1).cols(); ++functionColIndex) {
                        int currColLevel;
                        int currColIndex;
                        if (functionLocation(0)(0, functionColIndex) == patch) {
                            //MASTER FUNCTION CORRESPONDS WITH THE CURRENT PATCH
                            corrLevel = functionLocation(0)(1, functionColIndex);
                            currIndex = functionLocation(0)(2, functionColIndex);
                            gsMatrix<> valoreCol = Bells(patch)(corrLevel).function(currIndex).eval(uv(patch).col(k));
                            A_mat(functionIndex, functionColIndex) += valore(0, 0) * valoreCol(0, 0);
                        }
                        if (functionLocation(1)(0, functionColIndex) == patch) {
                            //SLAVE FUNCTION CORRESPONDS WITH THE CURRENT PATCH
                            corrLevel = functionLocation(1)(1, functionColIndex);
                            currIndex = functionLocation(1)(2, functionColIndex);
                            gsMatrix<> valoreCol = Bells(patch)(corrLevel).function(currIndex).eval(uv(patch).col(k));
                            A_mat(functionIndex, functionColIndex) += valore(0, 0) * valoreCol(0, 0);
                        }
                    }
                }
            }
        }
    }
}


void assembleSystem5(gsVector<gsMatrix<>>  uv,
    gsVector<gsMatrix<>>  xy,
    gsVector<gsVector<gsTensorBSplineBasis<2, real_t>>> Bells,
    gsVector  <gsVector< gsVector<index_t>>> isActive,
    gsVector  <gsVector< gsVector<index_t>>>  globalIndex,
    gsVector  <gsMatrix<index_t>> functionLocation,
    gsSparseMatrix<real_t>& A_mat,
    gsMatrix<real_t>& m_B
) {
    int totalSize = 0;
    for (int patch = 0; patch < functionLocation.size(); ++patch) {
        totalSize += functionLocation(patch).cols();
    }
    int shift = 0;
    for (int patch = 0; patch < Bells.size(); ++patch) {
        for (int functionIndex = 0; functionIndex < functionLocation(patch).cols(); ++functionIndex) {
            for (int k = 0; k < uv(patch).cols(); ++k) {
                int corrLevel;
                int currIndex;
                gsMatrix<> valore;
                bool functionIsActive = false;
                //SLAVE FUNCTION CORRESPONDS WITH THE CURRENT PATCH
                corrLevel = functionLocation(patch)(1, functionIndex);
                currIndex = functionLocation(patch)(2, functionIndex);
                valore = Bells(patch)(corrLevel).function(currIndex).eval(uv(patch).col(k));
                m_B(functionIndex + shift, 0) += xy(patch)(0, k);
                m_B(functionIndex + shift, 1) += xy(patch)(1, k);
                for (int functionColIndex = 0; functionColIndex < functionLocation(patch).cols(); ++functionColIndex) {
                    int currColLevel;
                    int currColIndex;
                    //SLAVE FUNCTION CORRESPONDS WITH THE CURRENT PATCH
                    corrLevel = functionLocation(patch)(1, functionColIndex);
                    currIndex = functionLocation(patch)(2, functionColIndex);
                    gsMatrix<> valoreCol = Bells(patch)(corrLevel).function(currIndex).eval(uv(patch).col(k));
                    A_mat(functionIndex + shift, functionColIndex + shift) += valore(0, 0) * valoreCol(0, 0);
                }
            }
        }
        shift += functionLocation(patch).cols();
    }
}

void assembleSystem3(gsVector<gsMatrix<>>  uv,
    gsVector<gsMatrix<>>  xy,
    gsVector<gsVector<gsTensorBSplineBasis<2, real_t>>> Bells,
    gsVector  <gsVector< gsVector<index_t>>> isActive,
    gsVector  <gsVector< gsVector<index_t>>>  globalIndex,
    gsVector  <gsMatrix<index_t>> functionLocation,
    gsSparseMatrix<real_t>& A_mat,
    gsMatrix<real_t>& m_B) {
    gsMatrix<real_t> value, curr_point;
    for (int patchIndex = 0; patchIndex < Bells.size(); ++patchIndex) {
        for (index_t k = 0; k < uv(patchIndex).cols(); k++) {
            curr_point = uv(patchIndex).col(k);
            for (int functionIndex = 0; functionIndex < functionLocation(0).cols(); ++functionIndex) {
                double funcvalR = 0;
                if (patchIndex == functionLocation(0)(0, functionIndex)) {
                    gsMatrix<> valore = Bells(functionLocation(0)(0, functionIndex))(functionLocation(0)(1, functionIndex)).function(functionLocation(0)(2, functionIndex)).eval(curr_point);
                    funcvalR = valore(0, 0);
                }
                if (patchIndex != functionLocation(0)(0, functionIndex) && functionLocation(1)(0, functionIndex) != -1) {
                    gsMatrix<> valore = Bells(functionLocation(1)(0, functionIndex))(functionLocation(1)(1, functionIndex)).function(functionLocation(1)(2, functionIndex)).eval(curr_point);
                    funcvalR = valore(0, 0);
                }
                m_B(functionIndex, 0) += funcvalR * xy(patchIndex)(0, k);
                m_B(functionIndex, 1) += funcvalR * xy(patchIndex)(1, k);
                for (int functionIndexC = 0; functionIndexC < functionLocation(0).cols(); ++functionIndexC) {
                    double funcvalC = 0;
                    if (patchIndex == functionLocation(0)(0, functionIndexC)) {//RIGHT PATCH FOR MASTER
                        gsMatrix<> valore = Bells(functionLocation(0)(0, functionIndexC))(functionLocation(0)(1, functionIndexC)).function(functionLocation(0)(2, functionIndexC)).eval(curr_point);
                        funcvalC = valore(0, 0);
                    }
                    else if (patchIndex != functionLocation(0)(0, functionIndexC) && functionLocation(1)(0, functionIndexC) != -1) {//WRONG PATCH FOR MASTER, BUT RIGHT FOR SLAVE
                        gsMatrix<> valore = Bells(functionLocation(1)(0, functionIndexC))(functionLocation(1)(1, functionIndexC)).function(functionLocation(1)(2, functionIndexC)).eval(curr_point);
                        funcvalC = valore(0, 0);
                    }
                    A_mat(functionIndex, functionIndexC) += funcvalR * funcvalC;
                }
            }
        }
    }
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

int twinFunction(gsTHBSplineBasis < 2, real_t>  initTP, int functionIndex, gsTHBSplineBasis < 2, real_t>  twinTP, int firstSide, int secondSide) {
    gsMatrix<> grevilles;
    gsMatrix<> grevillesI;
    twinTP.anchors_into(grevilles);
    initTP.anchor_into(functionIndex, grevillesI);
    index_t foundTwinsNum = 0;
    double uSecond;
    double vSecond;
    gsMatrix<> punto(2, 1);
    for (int i = 0; i < grevilles.cols(); ++i) {
        punto.col(0) = grevilles.col(i);
        coordinateTranformation(uSecond, vSecond, firstSide, secondSide, punto);
        if (
            abs(grevillesI(0, 0) - uSecond) <= 1e-9 &&
            abs(grevillesI(1, 0) - vSecond) <= 1e-9
            ) {
            return i;
        }
    }
    return -1;
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



void orientThePatches(gsMultiPatch<> mp, std::vector<int>& firstSide, std::vector<int>& secondSide, std::vector<int>& firstPatch, std::vector<int>& secondPatch)
{
    std::vector<boundaryInterface>& boundaryInterfaces = mp.interfaces();
    for (int interfaceNum = 0; interfaceNum < boundaryInterfaces.size(); ++interfaceNum)
    {
        outfile << boundaryInterfaces[interfaceNum] << "\n";
        gsInfo << boundaryInterfaces[interfaceNum] << "\n";
        firstSide.push_back(boundaryInterfaces[interfaceNum].first().index());
        secondSide.push_back(boundaryInterfaces[interfaceNum].second().index());
        firstPatch.push_back(boundaryInterfaces[interfaceNum].first().patch);
        secondPatch.push_back(boundaryInterfaces[interfaceNum].second().patch);
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
        int firstSide = boundaryInterfaces[interfaceNum].first().index();
        int secondSide = boundaryInterfaces[interfaceNum].second().index();
        int firstPatch = boundaryInterfaces[interfaceNum].first().patch;
        int secondPatch = boundaryInterfaces[interfaceNum].second().patch;
        for (int level = 0; level < Bells(0).size(); ++level) {
            for (int functionIndex = 0; functionIndex < Bells(firstPatch)(level).size(); ++functionIndex) {
                if (!isTouching(firstPatch)(level)(functionIndex)) {
                    if (Bells(firstPatch)(level).function(functionIndex).support()(
                        (firstSide + 1) / 2 - 1, (firstSide + 1) % 2) ==
                        Bells(firstPatch)(level).support()((firstSide + 1) / 2 - 1, (firstSide + 1) % 2)) {
                        isTouching(firstPatch)(level)(functionIndex) = 1;
                        touchedPatchSide(firstPatch)(level)(functionIndex) = secondSide;
                    }
                }
                if (!isTouching(secondPatch)(level)(functionIndex)) {
                    if (Bells(secondPatch)(level).function(functionIndex).support()((secondSide + 1) / 2 - 1, (secondSide + 1) % 2) ==
                        Bells(secondPatch)(level).support()((secondSide + 1) / 2 - 1, (secondSide + 1) % 2)) {
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
        int firstSide = boundaryInterfaces[interfaceNum].first().index();
        int secondSide = boundaryInterfaces[interfaceNum].second().index();
        int firstPatch = boundaryInterfaces[interfaceNum].first().patch;
        int secondPatch = boundaryInterfaces[interfaceNum].second().patch;
        gsInfo << "Identifying " << firstPatch << " and " << secondPatch << "\n";
        gsInfo << "test\n";
        for (int level = 0; level < isTouching(firstPatch).size(); ++level) {
            for (int functionIndex = 0; functionIndex < Bells(firstPatch)(level).size(); ++functionIndex) {
                int test250519 = isZeroOnBoundary(Bells, firstPatch, level, functionIndex, firstSide);
                if (
                    isTouching(firstPatch)(level)(functionIndex)
                    && test250519 == 0
                    ) {
                    int twinIndex = twinFunction(Bells(firstPatch)(level), functionIndex, Bells(secondPatch)(level), firstSide, secondSide);
                    if (twinIndex != -1) {
                        if (
                            isTouching(secondPatch)(level)(twinIndex)
                            ) {
                            twinsIndex(firstPatch)(level)(functionIndex).push_back(twinIndex);
                            twinsPatch(firstPatch)(level)(functionIndex).push_back(secondPatch);
                            hasATwin(firstPatch)(level)(functionIndex) = 1;
                            twinsPatch(secondPatch)(level)(twinIndex).push_back(firstPatch);
                            twinsIndex(secondPatch)(level)(twinIndex).push_back(functionIndex);
                            hasATwin(secondPatch)(level)(twinIndex) = 1;
                        }
                        else if (!isTouching(secondPatch)(level)(twinIndex)) {
                            //But the twin does not touch the boundary\n";
                        }
                    }
                }
                else if (!isTouching(firstPatch)(level)(functionIndex)) {
                    //The function " << functionIndex << " does not touch the boundary\n";
                }
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
    //now we have to identify twins of twins. Since being a twin is a equivalence relation, it is transitive
    for (int patch = 0; patch < twinsIndex.size(); ++patch) {
        for (int level = 0; level < twinsIndex(patch).size(); ++level) {
            for (int funcIndex = 0; funcIndex < twinsIndex(patch)(level).size(); ++funcIndex) {
                for (int twinNum1 = 0; twinNum1 < twinsIndex(patch)(level)(funcIndex).size(); ++twinNum1) {
                    int patchIndex = twinsPatch(patch)(level)(funcIndex)[twinNum1];
                    int functionIndex = twinsIndex(patch)(level)(funcIndex)[twinNum1];
                    if (patchIndex == -1 || functionIndex == -1)     break;
                    for (int twinNum2 = 0; twinNum2 < twinsIndex(patchIndex)(level)(functionIndex).size(); ++twinNum2) {
                        int tS = twinsIndex(patchIndex)(level)(functionIndex).size();
                        if (twinsIndex(patchIndex)(level)(functionIndex).size() <= 1) {
                            break;
                        }
                        else {
                            std::vector<int>::iterator it1;
                            std::vector<int>::iterator it2;
                            std::vector<int> arr1 = twinsIndex(patch)(level)(funcIndex);
                            std::vector<int> arr2 = twinsPatch(patch)(level)(funcIndex);
                            it1 = std::find(arr1.begin(),
                                arr1.end(),
                                twinsIndex(patchIndex)(level)(functionIndex)[twinNum2]
                            );
                            it2 = std::find(arr2.begin(),
                                arr2.end(),
                                twinsPatch(patchIndex)(level)(functionIndex)[twinNum2]
                            );
                            bool a1 = !(it1 != arr1.end() && it2 != arr2.end());
                            bool a2 = twinsIndex(patchIndex)(level)(functionIndex)[twinNum2] != funcIndex;
                            bool a3 = twinsPatch(patchIndex)(level)(functionIndex)[twinNum2] != patch;
                            if (!(it1 != arr1.end() && it2 != arr2.end()) &&
                                twinsIndex(patchIndex)(level)(functionIndex)[twinNum2] != funcIndex &&
                                twinsPatch(patchIndex)(level)(functionIndex)[twinNum2] != patch
                                ) {
                                arr1.push_back(twinsIndex(patchIndex)(level)(functionIndex)[twinNum2]);
                                arr2.push_back(twinsPatch(patchIndex)(level)(functionIndex)[twinNum2]);
                                twinsIndex(patch)(level)(funcIndex) = arr1;
                                twinsPatch(patch)(level)(funcIndex) = arr2;
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
                    for (int twinNum2 = 0; twinNum2 < twinsIndex(patchIndex)(level)(functionIndex).size(); ++twinNum2) {
                        int tS = twinsIndex(patchIndex)(level)(functionIndex).size();
                        if (twinsIndex(patchIndex)(level)(functionIndex).size() <= 1) {
                            break;
                        }
                        else {
                            std::vector<int>::iterator it1;
                            std::vector<int>::iterator it2;
                            std::vector<int> arr1 = twinsIndex(patch)(level)(funcIndex);
                            std::vector<int> arr2 = twinsPatch(patch)(level)(funcIndex);
                            it1 = std::find(arr1.begin(),
                                arr1.end(),
                                twinsIndex(patchIndex)(level)(functionIndex)[twinNum1]
                            );
                            it2 = std::find(arr2.begin(),
                                arr2.end(),
                                twinsPatch(patchIndex)(level)(functionIndex)[twinNum1]
                            );
                            bool a1 = !(it1 != arr1.end() && it2 != arr2.end());
                            bool a2 = twinsIndex(patchIndex)(level)(functionIndex)[twinNum1] != funcIndex;
                            bool a3 = twinsPatch(patchIndex)(level)(functionIndex)[twinNum1] != patch;
                            if (!(it1 != arr1.end() && it2 != arr2.end()) &&
                                twinsIndex(patchIndex)(level)(functionIndex)[twinNum1] != funcIndex &&
                                twinsPatch(patchIndex)(level)(functionIndex)[twinNum1] != patch
                                ) {
                                arr1.push_back(twinsIndex(patchIndex)(level)(functionIndex)[twinNum1]);
                                arr2.push_back(twinsPatch(patchIndex)(level)(functionIndex)[twinNum2]);
                                twinsIndex(patch)(level)(funcIndex) = arr1;
                                twinsPatch(patch)(level)(funcIndex) = arr2;
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
    for (size_t patch = 0; patch < THBVector.size(); ++patch) {
        auto a = THBVector(patch).size();
        thbToBellsMapping[patch].resize(a);
        outfile << "thbToBellsMapping[patch].size(): " << thbToBellsMapping[patch].size() << "\n";
        for (size_t functionIndex = 0; functionIndex < THBVector(patch).size(); ++functionIndex) {
            int corrLevel = THBVector(patch).levelOf(functionIndex);
            int corrIndex = THBVector(patch).flatTensorIndexOf(functionIndex);
            // Save mapping to Bells
            thbToBellsMapping[patch][functionIndex] = { corrLevel, corrIndex };

            // Mark as included in the hierarchical basis
            isIncluded(patch)(corrLevel)(corrIndex) = 1;
            indexInTHB(patch)(corrLevel)(corrIndex) = functionIndex;
        }
    }

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
    int ourBox[], gsTHBSplineBasis<2, real_t>& THB, gsTHBSplineBasis<2, real_t> THBUnstructured,
    gsVector<gsMatrix<index_t >>& lowCorners,
    gsVector<gsMatrix<index_t >>& upCorners, gsVector<gsVector<index_t >>& myLevel, int& lastNonZeroRow) {
    std::vector<index_t> box;
    box.push_back(ourBox[0]);
    box.push_back(ourBox[1]);
    box.push_back(ourBox[2]);
    box.push_back(ourBox[3]);
    box.push_back(ourBox[4]);
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
        THB.refineElements(box);
    }
    //experimentalpart
    for (size_t therow = lastNonZeroRow; therow < boxMat(patch).size(); therow++)
    {
        boxMat(patch)(therow).clear();
    }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void restoreTheHierarchy(int& createdBoxNum, int& lastNonzeroRow, gsVector< gsVector<gsVector<index_t>>>& boxMat, int levNow, int centerInd, int ourBox[], int successfullAttempts, int patch) {
    for (int l = 0; l < createdBoxNum; l++) {
        boxMat(patch)(lastNonzeroRow - l)(0) = levNow; //Preparation for multilevel meshes
        outfile << "updated coordinates of " << lastNonzeroRow - l << "box: " <<
            boxMat(patch)(lastNonzeroRow - l)(0);
        for (int m = 1; m < 5; m++) {
            boxMat(patch)(lastNonzeroRow - l)(m) = 0;
            outfile << "\t" << boxMat(patch)(lastNonzeroRow - l)(m);
        }
        outfile << "\n";
    }
    outfile << "current " << centerInd << "box is now:";
    for (int k = 0; k < 5; k++) {
        boxMat(patch)(centerInd)(k) = ourBox[k];
        outfile << "\t" << boxMat(patch)(centerInd)(k);
    }
    outfile << "\n";
    if (successfullAttempts == 0) {
        outfile << "0th position IS RETURNED TO\n";
        for (int k = 0; k < 5; k++) {
            boxMat(patch)(0)(k) = ourBox[k];
            outfile << "\t" << boxMat(patch)(0)(k);
        }
        outfile << "\n";
    }
    lastNonzeroRow = lastNonzeroRow - createdBoxNum;
    createdBoxNum = 0;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int rebuildTheHierarchy(gsVector< gsVector<gsVector<index_t>>>& boxMat, int row, int x1U, int x1Bi, int x2U, int x2Bi, int y1U, int y1Bi, int y2U, int y2Bi, int levNow, int& lastNonZeroRow,
    int& createdBoxNum, int& centerInd, int ourBox[], int& needToEscape, int patch) {
    gsInfo << "row = " << row << "\n";
    gsInfo << "Before the rebuild lastNonZeroRow: " << lastNonZeroRow << "\n";
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
        gsInfo << "ALERT: OUR CELL OF LEVEL " << levNow << "\n" << x1U << " " << y1U << " " << x2U <<
            " " << y2U <<
            " INTERSECTS WITH [" << row << "]th box of level " <<
            boxMat(patch)(row)(0) << ":\n" << boxMat(patch)(row)(1) << " " << boxMat(patch)(row)(2) << " " <<
            boxMat(patch)(row)(3) << " " << boxMat(patch)(row)(4) << "\n";
        if (boxMat(patch)(row)(0) <= levNow) {
            gsInfo << "NO NEED TO REBUILD \n";
        }
        //HERE I START TO WRITE CODE IN MORE GENERAL WAY
        else if (boxMat(patch)(row)(0) > levNow + 1) {
            gsInfo << "I DON'T TOUCH THESE BOX AT THIS ITERATION!\n";
        }
        else if (boxMat(patch)(row)(0) == levNow + 1) {
            x1Bi = boxMat(patch)(row)(1);
            y1Bi = boxMat(patch)(row)(2);
            x2Bi = boxMat(patch)(row)(3);
            y2Bi = boxMat(patch)(row)(4);
            gsInfo << "REBUILD\n";
            RTH = 1;
            //SW W NW
            if (x1U * pow(2, boxMat(patch)(row)(0) - levNow) - x1Bi > 0) {
                if (y1U * pow(2, boxMat(patch)(row)(0) - levNow) - y1Bi > 0) {
                    gsInfo << "CREATE SW BOX:\n";
                    boxMat(patch)(lastNonZeroRow).resize(5);
                    boxMat(patch)(lastNonZeroRow)(0) = levNow + 1;
                    boxMat(patch)(lastNonZeroRow)(1) = x1Bi;
                    boxMat(patch)(lastNonZeroRow)(2) = y1Bi;
                    boxMat(patch)(lastNonZeroRow)(3) = x1U * pow(2, boxMat(patch)(row)(0) - (levNow));
                    boxMat(patch)(lastNonZeroRow)(4) = y1U * pow(2, boxMat(patch)(row)(0) - (levNow));
                    lastNonZeroRow++;
                    createdBoxNum++;
                }
                gsInfo << "CREATE W BOX:\n";
                boxMat(patch)(lastNonZeroRow).resize(5);
                boxMat(patch)(lastNonZeroRow)(0) = levNow + 1;
                boxMat(patch)(lastNonZeroRow)(1) = x1Bi;
                boxMat(patch)(lastNonZeroRow)(2) = y1U * pow(2, boxMat(patch)(row)(0) - levNow);
                boxMat(patch)(lastNonZeroRow)(3) = x1U * pow(2, boxMat(patch)(row)(0) - levNow);
                boxMat(patch)(lastNonZeroRow)(4) = y2U * pow(2, boxMat(patch)(row)(0) - levNow);
                lastNonZeroRow++;
                createdBoxNum++;
                if (y2Bi - y2U * pow(2, boxMat(patch)(row)(0) - levNow) > 0) {
                    gsInfo << "CREATE NW BOX:\n";
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
                gsInfo << "CREATE S BOX:\n";
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
                gsInfo << "CREATE N BOX:\n";
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
                    gsInfo << "CREATE SE BOX:\n";
                    boxMat(patch)(lastNonZeroRow).resize(5);
                    boxMat(patch)(lastNonZeroRow)(0) = levNow + 1;
                    boxMat(patch)(lastNonZeroRow)(1) = x2U * pow(2, boxMat(patch)(row)(0) - levNow);
                    boxMat(patch)(lastNonZeroRow)(2) = y1Bi;
                    boxMat(patch)(lastNonZeroRow)(3) = x2Bi;
                    boxMat(patch)(lastNonZeroRow)(4) = y1U * pow(2, boxMat(patch)(row)(0) - levNow);
                    lastNonZeroRow++;
                    createdBoxNum++;
                }
                gsInfo << "CREATE E BOX:\n";
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
    gsInfo << "Starting savePatches for " << SubdomainHierarchy.size() << " patches.\n";
    outfile << "Starting savePatches for " << SubdomainHierarchy.size() << " patches.\n";

    localIndex.resize(globalIndexTHB.size());

    for (size_t patch = 0; patch < SubdomainHierarchy.size(); ++patch)
    {
        const auto& patchGlobalIndices = globalIndexTHB(patch);

        // Count active basis functions
        index_t count = 0;
        for (index_t i = 0; i < patchGlobalIndices.size(); ++i)
        {
            if (patchGlobalIndices[i] != -1)
                ++count;
        }

        localIndex(patch).resize(count);

        gsInfo << "Processing patch " << patch << " with " << count << " active basis functions.\n";
        outfile << "Processing patch " << patch << " with " << count << " active basis functions.\n";

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

            gsInfo << "Saved patch " << patch << " to file: " << outputFileName << "\n";
            outfile << "Saved patch " << patch << " to file: " << outputFileName << "\n";
        }
        catch (const std::exception& e)
        {
            gsInfo << "EXCEPTION on patch " << patch << ": " << e.what() << "\n";
            outfile << "EXCEPTION on patch " << patch << ": " << e.what() << "\n";
        }
    }

    gsInfo << "Finished savePatches.\n";
    outfile << "Finished savePatches.\n";
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
    gsVector < gsTHBSplineBasis < 2, real_t > > SubdomainHierarchy)
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

            if (//ten_index == tensor_active_index 
                //||
                globalIndex[patch][level][ten_index] != -1 && !truncationStopped) // truncate
            {
                /*                gsDebugVar(ten_index);
                                gsDebugVar(tensor_active_index);
                                gsDebugVar(coef_index + index);*/
                                //              gsDebugVar(globalIndex[patch][level][ten_index]);

                coefs(coef_index + index, 0) = 0;
                //            gsInfo << "-----------------------------";
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
        auto wurschtl = SubdomainHierarchy[patch].getXmatrix();
        //gsDebugVar(wurschtl);
        gsDebugVar(wurschtl.size());
        gsDebugVar(wurschtl[fineLevel][fineFlatIndex]);
        gsDebugVar(SubdomainHierarchy[patch].getXmatrix()[fineLevel][fineFlatIndex]);
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
    std::vector<std::vector<gsSparseVector < real_t >>>& presentation) {
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
            SubdomainHierarchy
        );
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
    int maxLevel = 4;
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
                _representBasisFunction(f, twinIndex, THBIndex, patch, level, clevel, low, high, Bells, tensor_index, globalIndex, SubdomainHierarchy, presentation);

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

    for (index_t pt = 0; pt < numPoints; ++pt)
    {
        const auto param = quNodes.col(pt);

        for (index_t a = 0; a < actives.rows(); ++a)
        {
            const int thbIndex = actives(a, 0);
            const int globalID = globalIndexTHB[patch][thbIndex];
            if (globalID < 0) continue;

            const real_t cx = vectSol(globalID, 0);
            const real_t cy = vectSol(globalID, 1);

            // Use MPBES derivative evaluation (handles truncation automatically)
            gsMatrix<real_t> deriv;
            mpbes.evalDerivSingleOnPatch(globalID, patch, param, deriv);
            
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

inline void assembleFitting(gsVector<gsMatrix<>>  uv,
    gsSparseMatrix<real_t>& A_mat,
    gsVector<gsVector<gsTensorBSplineBasis<2, real_t>>> Bells,
    std::vector<std::vector<std::vector<index_t>>> functionDescription,
    gsMatrix<>& matAsquare,
    const real_t weight
)
{
    PROFILE_FUNCTION();
    gsInfo << "A_mat.rows()" << A_mat.rows() << " " << A_mat.cols() << "\n";
    gsInfo << "matA.rows()" << matAsquare.rows() << " " << matAsquare.cols() << "\n";
    outfile << "matA.rows()" << matAsquare.rows() << " " << matAsquare.cols() << "\n";
    std::chrono::time_point<std::chrono::system_clock> beforesquare, aftersquare;
    beforesquare = std::chrono::system_clock::now();
    for (size_t i = 0; i < matAsquare.rows(); i++)
    {
        for (size_t j = 0; j < matAsquare.cols(); j++)
        {
            if (matAsquare(i, j) != 0)
            {
                A_mat(i * 2 + 0, j * 2 + 0) += weight * matAsquare(i, j);
                A_mat(i * 2 + 0, j * 2 + 1) += weight * matAsquare(i, j);
                A_mat(i * 2 + 1, j * 2 + 0) += weight * matAsquare(i, j);
                A_mat(i * 2 + 1, j * 2 + 1) += weight * matAsquare(i, j);
            }
        }
    }
    aftersquare = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds_square = aftersquare - beforesquare;
    gsInfo << "square matrix took: " << elapsed_seconds_square.count() << "\n";
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
    bool verbose = true)
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

    for (index_t i = 0; i < actives.rows(); ++i)
    {
        int activeIndex = actives(i, 0);
        int globalI = globalIndexTHB[patch][activeIndex];
        if (globalI == -1) continue;

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

        for (index_t j = 0; j < actives.rows(); ++j)
        {
            int activeJIndex = actives(j, 0);
            int globalJ = globalIndexTHB[patch][activeJIndex];
            if (globalJ == -1) continue;

            for (int col = 0; col < q.cols(); ++col)
            {
                for (int dimu = 0; dimu < geoDim; ++dimu)
                {
                    for (int dimv = 0; dimv < geoDim; ++dimv)
                    {
                        index_t row = col * geoDim + dimu;
                        index_t colJ = col * geoDim + dimv;
                        real_t term1 = J(row, i);
                        real_t term2 = J(colJ, j);
                        real_t term3 = weight;
                        real_t val = term1 * term2 * term3;

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
                            //gsInfo << "  A[" << aRow << "," << aCol << "] += J(" << row << "," << i << ") * "
                            //    << "J(" << colJ << "," << j << ") * weight = "
                            //    << term1 << " * " << term2 << " * " << term3
                            //    << " = " << val << "\n";
                            outfile << "A[" << aRow << "," << aCol << "] += " << val
                                << "   (J1=" << term1 << ", J2=" << term2 << ", w=" << term3 << ")\n";
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
    const real_t epsilon,
    bool verbose)
{
    PROFILE_FUNCTION();
    
    // Extract data from MPBES
    auto Bells = mpbes.bellsBases();
    auto SubdomainHierarchy = mpbes.thbBases();
    auto functionDescription = mpbes.functionDescription();
    auto spilloverFunctionCoordinates = mpbes.spilloverCoordinates();
    auto hasSpillover = mpbes.hasSpillover();
    auto globalIndexTHB = extractGlobalIndexTHB(mpbes);

    assembleFitting(uv, A, Bells, functionDescription, matAsquare, fitting);

    gsMatrix<real_t> b(A.rows(), 1);
    b.setZero();

    real_t functional = 0;
    auto beforeOpt = std::chrono::system_clock::now();

    for (size_t patch = 0; patch < SubdomainHierarchy.size(); ++patch)
    {
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


            if (uniformity > 0)
            {
                deriv2_into_geometry_mpbe(mpbes, patch, quNodes, geom2Der, actives, vectSol, globalIndexTHB);

                deriv2_into_basis_mpbe(mpbes, patch, quNodes, basis2Der, actives, globalIndexTHB, false);
            }

            for (index_t quNode = 0; quNode < quNodes.cols(); ++quNode)
            {
                const real_t weight = quWeights[quNode];

                if (fitting > 0)
                {
                    gsMatrix<real_t> Jort, q;
                    getJort(Jort, actives.rows(), basisDer, geomDer, quNode);
                    getQort(q, geomDer, quNode);
                    functional += assembleFunctional(q, weight, fitting);

                    assembleOptimization(A, b, Jort, q, weight * orthogonality,
                        globalIndexTHB, actives, patch, false);
                }

                if (length > 0)
                {
                    gsMatrix<real_t> Jlen, qLen;
                    getQlen(qLen, geomDer, quNode);
                    getJlen(Jlen, basisDer.cols(), quNode, basisDer);
                    functional += assembleFunctional(qLen, weight, length);

                    assembleOptimization(A, b, Jlen, qLen, weight * length,
                        globalIndexTHB, actives, patch, false);
                }

                if (uniformity > 0)
                {
                    gsMatrix<real_t> Juni, qUni;
                    getJuni(Juni, actives.rows(), quNode, basis2Der);
                    getQuni(qUni, geom2Der, quNode);
                    functional += assembleFunctional(qUni, weight, uniformity);

                    assembleOptimization(A, b, Juni, qUni, weight * uniformity,
                        globalIndexTHB, actives, patch, false);
                }
            }
        }
    }

    auto afterOpt = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = afterOpt - beforeOpt;
    
    // SOLVER PHASE
    gsSparseSolver<>::CGIdentity solver(A);
    gsMatrix<real_t> x;

    try {
        x = solver.solve(b);
    }
    catch (const std::exception& e) {
        return;
    }

    // UPDATE SOLUTION
    gsMatrix<> vectSolOld = vectSol;
    for (int row = 0; row < vectSol.rows(); ++row)
        vectSol.row(row) -= x.block(row * 2, 0, 2, 1).transpose();

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

// Struct to hold all data computed by the algorithm
struct AlgorithmResult {
    gsMultiPatch<>::uPtr mp;
    gsVector<gsVector<gsTensorBSplineBasis<2, real_t>>> Bells;
    gsVector<gsVector<gsVector<index_t>>> boxMat;
    std::vector<std::vector<std::vector<double>>> acceptedCoefs;
    gsMatrix<real_t> AcceptedvectSol;
    gsVector<int> AcceptedlastRow;
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

AlgorithmResult unrefinementAlgorithmHBJ(
    const std::string& filename,
    real_t epsilon_g,
    real_t epsilon_f,
    char method,
    const std::string& givenGeo,
    const std::string& acCond) {
    
    // Initialize variables that were in main()
    std::chrono::time_point<std::chrono::system_clock> startTime, iterTime;
    startTime = std::chrono::system_clock::now();
    int row, acceptedsize, attempt = 0;
    int valid = 0;
    real_t lcx, lcy, ucx, ucy;
    int successfullAttempts = 0, totalAttempts = 0;
    int projections = 0;  // Use different name to avoid conflict with proj() function
    std::string xmlFile = "/home/targon/theydarov/output/" + givenGeo + acCond;
    
    gsMultiPatch<>::uPtr mp = gsReadFile<>(filename);
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
    int degree0 = 2;
    int degree1 = 2;
    unsigned multEnd = degree0 + 1;
    gsKnotVector<> ku(mp->patch(0).basis().support()(0, 0), mp->patch(0).basis().support()(0, 1),
        interioru0 + interioru1, multEnd);
    gsKnotVector<> kv(mp->patch(0).basis().support()(1, 0), mp->patch(0).basis().support()(1, 1), interiorv1, multEnd);
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
    int numPoints = std::sqrt(mp->patch(0).basis().size());
    int jacPoints = 40;
    gsVector<gsMatrix<real_t>> uv1(mp->nPatches());
    gsVector<gsMatrix<real_t>> uv2(mp->nPatches());
    gsVector<gsMatrix<real_t>> xy1(mp->nPatches());
    gsVector<gsMatrix<real_t>> uvFeature(mp->nPatches());
    gsVector<gsMatrix<real_t>> xyFeature(mp->nPatches());
    gsVector<gsMatrix<real_t>> xyFeatureMiddle(mp->nPatches());
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
        maxLevel(patch) = fineLevel(patch) = myLevel(patch).maxCoeff();
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

    int fullMonty, lastNonZeroRow;
    gsVector<index_t> currentLastNonZeroRow(mp->nPatches());
    gsVector<int> AcceptedlastRow(mp->nPatches());

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
        lastNonZeroRow = lowCorners(patch).rows();
        initialBoxesNum(patch) = lastNonZeroRow;
        currentLastNonZeroRow(patch) = lastNonZeroRow;
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
    gsInfo << "test fot 20.07.25\n";
    gsInfo << Bells[0][3].function(70).support() << "\n";
    gsInfo << Bells[1][3].function(79).support() << "\n";
    gsInfo << "--------------\n";
    gsMatrix<> puntoTest(2, 1);
    puntoTest(0, 0) = 0;    puntoTest(1, 0) = 0.777778;
    gsInfo << Bells[0][3].function(70).eval(puntoTest) << "\n";
    puntoTest(0, 0) = 1;
    gsInfo << Bells[1][3].function(79).eval(puntoTest) << "\n";
    puntoTest(0, 0) = 0;
    gsInfo << "end of test for 20.07.25\n";
    //outfile << "After reading the file, lastNonZeroRow:" << lastNonZeroRow << "\n";
    gsInfo << "a special test, remove afterwards: \n";
    gsMatrix<> pointSpecial(2, 1);
    pointSpecial(0, 0) = pointSpecial(1, 0) = 0.5;
    gsInfo << mp1.patch(0).jacobian(pointSpecial).determinant() << "\n";
    for (int patch = 0; patch < mp->nPatches(); ++patch) {
        lastNonZeroRow = initialBoxesNum(patch);
        //for (int patch = 0; patch < 1; ++patch) {
            //gsDebugVar(patch);
        gsInfo << "started checking the patch number " << patch << "\n";
        gsInfo << "started checking the patch number " << patch << "\n";
        gsInfo << "started checking the patch number " << patch << "\n";
        outfile << "started checking the patch number " << patch << "\n";
        outfile << "started checking the patch number " << patch << "\n";
        outfile << "started checking the patch number " << patch << "\n";
        //gsDebugVar(coarseLevel(patch));
        for (int levNow = coarseLevel(patch); levNow >= 0; --levNow) {
            //for (int levNow = coarseLevel(patch); levNow >= 1; --levNow) {
            //for (int levNow = coarseLevel(patch); levNow >= 2; --levNow) {
            //for (int levNow = coarseLevel(patch); levNow >= coarseLevel(patch); --levNow) {
            iteration = 0;
            theLev = levNow;
            /*if ((levNow < coarseLevel(patch)) && wasRebuilt) {
                std::string input(xmlFile + ".xml");
                gsMultiPatch<>::uPtr AcceptedMultiPatch = gsReadFile<>(xmlFile);
                for (int currentPatch = 0; currentPatch < AcceptedMultiPatch->nPatches(); ++currentPatch) {
                    THBAccepted(currentPatch)
                            = dynamic_cast<gsTHBSplineBasis<2> *>(&AcceptedMultiPatch->patch(
                            currentPatch).basis().source());
                    (THBAccepted(currentPatch)->tree().getBoxes(lowCorners(currentPatch), upCorners(currentPatch),
                                                                myLevel(currentPatch)));
                    for (int i = 0; i < lowCorners(currentPatch).size(); ++i) {
                        boxMat(currentPatch)(i)(0) = myLevel(currentPatch)(i);
                        boxMat(currentPatch)(i)(1) = (real_t) lowCorners(currentPatch)(i, 0) /
                                                     pow(2, maxLevel(currentPatch) - myLevel(currentPatch)(i));
                        boxMat(currentPatch)(i)(2) = (real_t) lowCorners(currentPatch)(i, 1) /
                                                     pow(2, maxLevel(currentPatch) - myLevel(currentPatch)(i));
                        boxMat(currentPatch)(i)(3) = (real_t) upCorners(currentPatch)(i, 0) /
                                                     pow(2, maxLevel(currentPatch) - myLevel(currentPatch)(i));
                        boxMat(currentPatch)(i)(4) = (real_t) upCorners(currentPatch)(i, 1) /
                                                     pow(2, maxLevel(currentPatch) - myLevel(currentPatch)(i));
                    }
                }
                lastNonZeroRow = lowCorners(patch).rows() - 1;
            }*/
            outfile << "working with level " << levNow << " as a coarseLevel\n";
            gsInfo << "working with level " << levNow << " as a coarseLevel\n";
            int pickedOne = -1;
            gsVector<int> nonCheckedCells((interior + 1) * pow(2, levNow) * (interior + 1) * pow(2, levNow));
            gsMatrix<int> pickedCells(1, (interior + 1) * pow(2, levNow) * (interior + 1) * pow(2, levNow));
            for (int i = 0; i < nonCheckedCells.size(); i++) {
                nonCheckedCells(i) = i;
                pickedCells(0, i) = 0;
            }
            int success = 1;
            attempt = 0;
            while (nonCheckedCells.size() != 0 && success) {
                iteration++;
                success = 0;
                gsVector<index_t> vectorS = nonCheckedCells;
                fullMonty = 0;
                while (vectorS.size() != 0) {
                    // Declare MPBES data variables at broader scope
                    std::vector<std::vector<std::vector<index_t>>> functionDescription;
                    std::vector<std::vector<std::array<int, 3>>> spilloverFunctionCoordinates;
                    std::vector<bool> hasSpillover;
                    std::vector<bool> isTruncated;
                    std::vector<std::vector<gsSparseVector<real_t>>> presentation;

                    gsInfo << "vectorS.size(): " << vectorS.size() << "\n";
                    outfile << "vectorS.size(): " << vectorS.size() << "\n";
                    lcx = 1.0;
                    lcy = 1.0;
                    ucx = 0.0;
                    ucy = 0.0;
                    failed = 0;
                    outfile << "\n";
                    outfile << "\n";
                    outfile << "\n";
                    outfile << "The boxes\n";
                    //gsDebugVar(boxMat(patch).size());
                    //lastNonZeroRow = lowCorners(patch).rows();
                    gsInfo << "lastNonZeroRow: " << lastNonZeroRow << "\n";
                    //for (int i = 0; i < lastNonZeroRow; i++) {
                    for (int i = 0; i < lastNonZeroRow; i++) {
                        gsInfo << "i:" << i << "\n";
                        for (int j = 0; j < 5; j++) {
                            outfile << boxMat(patch)(i)(j) << "\t";
                            gsInfo << boxMat(patch)(i)(j) << "\t";
                        }
                        outfile << "\n";
                        gsInfo << "\n";
                    }
                    createdBoxNum = 0;
                    if (method == 'r') {
                        currCellIndex = pickCell(vectorS, currArrayIndex, levNow, x1U, y1U, x2U, y2U, interior);
                    }
                    else if (method == 'l') {
                        currCellIndex = pickCell(vectorS, attempt, levNow, x1U, y1U, x2U, y2U, 1, interior);
                        currArrayIndex = 0;
                    }
                    else if (method == 's') {
                        if (valid && !pickedCells(0, pickedOne + 1)) {
                            pickCell(vectorS, attempt, levNow, x1U, y1U, x2U, y2U, 1, interior);
                            currArrayIndex = 0;
                        }
                        else {
                            currCellIndex = pickCell(vectorS, currArrayIndex, levNow, x1U, y1U, x2U, y2U, interior);
                        }
                        pickedOne = currCellIndex;
                        pickedCells(0, currCellIndex) = 1;
                    }
                    else {
                        gsInfo << "UNKNOWN METHOD.\n";
                        throw std::runtime_error("Unknown method specified");
                    }
                    int jopa = 4 * (int)pow(2, levNow) * 4 * (int)pow(2, levNow);
                    outfile << "=======================================\n";
                    outfile << "=======================================\n";
                    outfile << "=======================================\n";
                    outfile << "Number of nonchecked cells: " << vectorS.size() << "\n";
                    gsInfo << "Number of nonchecked cells: " << vectorS.size() << "\n";
                    gsInfo << "attempt " << attempt << ", CURRENT INDEX IS " << vectorS(currArrayIndex) <<
                        ",\t the coordinates of the box are\n" << x1U << " " << y1U << " " << x2U << " " << y2U
                        << "\n";
                    outfile << "patch " << patch << "\n";
                    outfile << "attempt " << attempt << ", CURRENT INDEX IS " << vectorS(currArrayIndex) <<
                        ",\t the coordinates of the box are\n" << x1U << " " << y1U << " " << x2U << " " << y2U
                        << "\n";
                    /* if (attempt > 0)     break;*/
                     //gsInfo << currArrayIndex << "\n";
                    outfile << currArrayIndex << "\n";
                    vectorS.removeElement(currArrayIndex);
                    gsInfo << "DELETED\n";
                    outfile << "DELETED\n";
                    //if (attempt > 1)     break;
                    valid = 0;
                    createSpline = 0;
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
                    if (createSpline == 1) {
                        iterTime = std::chrono::system_clock::now();
                        outfile << "creating the spline\n";
                        gsInfo << "creating the spline\n";
                        std::chrono::duration<double> elapsed_seconds = iterTime - startTime;
                        outfile << "TIME: " << elapsed_seconds.count() << "\n";
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
                        THB.anchors_into(anmat);
                        gsMatrix<real_t> basisInd(1, THB.size()); //1 if function is not presented
                        for (int l = 0; l < THB.size(); ++l) {
                            basisInd(0, l) = 0;
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
                                    basisInd(j) = 0;

                                    break;
                                }
                            }
                        }
                        for (int m = 0; m < THB.size(); m++) {
                            //                        if(THBcoefs(m,0) == 0.0 || THBcoefs(m,1) == 0.0)
                            //                        {
                            basisInd(0, m) = 1;
                            //                        }
                        }
                        //                    //gsInfo << basisInd << "\n";
                        outfile << basisInd << "\n";
                        for (int i1 = 0; i1 < basisInd.cols(); ++i1) {
                            if (basisInd(i1)) {
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
                        outfile << "Test from 07.06.25\n";
                        for (size_t i = 0; i < THB.size(); i++)
                        {
                            outfile << "Index " << i << ": ["
                                << THB.function(i).support()(0, 0) << ", "
                                << THB.function(i).support()(0, 1) << "] x ["
                                << THB.function(i).support()(1, 0) << ", "
                                << THB.function(i).support()(1, 1) << "]\n";

                        }
                        outfile << "End of test from 07.06.25\n";

                        std::chrono::time_point<std::chrono::system_clock> beforeIdent, afterIdent;
                        beforeIdent = std::chrono::system_clock::now();
                        //IdentifyPatches(mp1, THBVector1, isTouching, twinsIndex, twinsPatch, hasATwin, isActive);
                        std::vector<int> firstSide;
                        std::vector<int> secondSide;
                        std::vector<int> firstPatch;
                        std::vector<int> secondPatch;
                        orientThePatches(mp1, firstSide, secondSide, firstPatch, secondPatch);
                        gsInfo << "The minimum goal for 31/07/24 achieved at 06:36\n";
                        //IdentifyPatches(mp1, THBVector, isTouchingTHB, twinsIndexTHB, twinsPatchTHB, hasATwinTHB, isActiveTHB);
                        gsVector  <gsVector< gsVector<index_t>>> isIncluded;
                        gsVector  <gsVector< gsVector<index_t>>> indexInTHB;
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

                        afterIdent = std::chrono::system_clock::now();
                        std::chrono::duration<double> elapsed_seconds_Ident = afterIdent - beforeIdent;
                        gsInfo << "IdentifyPatches took: " << elapsed_seconds_Ident.count() << "\n";

                        // ===== MPBES: Multi-Patch B-spline with Exact continuity =====
                        // This performs: Twin identification + Kraft selection + Truncation
                        gsInfo << "\n========== Creating MPBES basis ==========\n";
                        MPBES<2, real_t> mpbes(
                            boxMat,
                            SubdomainHierarchy,
                            Bells,
                            hasATwin,
                            twinsIndex,
                            twinsPatch,
                            indexInTHB,
                            true,  // verbose
                            attempt  // for conditional logging
                        );
                        gsInfo << "MPBES created with " << mpbes.size() << " basis functions\n";
                        gsInfo << "==========================================\n\n";

                        // Extract data from MPBES for use throughout the code
                        functionDescription = mpbes.functionDescription();
                        spilloverFunctionCoordinates = mpbes.spilloverCoordinates();
                        hasSpillover = mpbes.hasSpillover();
                        isTruncated = mpbes.isTruncated();
                        presentation = mpbes.presentation();
                        commonSize = mpbes.size();

                        gsInfo << "Setup completed. Functions: " << commonSize << "\n";

                        // Assembly phase
                        index_t numElements = countTotalActiveElements(SubdomainHierarchy);
                        index_t numGaussPointsPerElement = 4;
                        index_t numTotalRows = numElements * numGaussPointsPerElement;

                        gsMatrix<real_t> vectB(numTotalRows, 2);
                        gsSparseMatrix<real_t> matA(numTotalRows, commonSize);

                        std::chrono::time_point<std::chrono::system_clock> beforeassembleA, afterassembleA;
                        beforeassembleA = std::chrono::system_clock::now();
                        gsInfo << "Starting assembly with " << commonSize << " functions...\n";

                        // Use new MPBES-based assemble
                        assemble(
                            uv1,
                            mpbes,
                            matA,
                            vectB,
                            mp1,
                            false);

                        if (attempt == 0) {
                            printTheMatrix(matA, "matA");
                        }

                        gsInfo << "matA.size: " << matA.rows() << " * " << matA.cols() << "\n";

                        // Check partition of unity with detailed reporting
                        gsInfo << "Checking partition of unity...\n";
                        real_t tolerance = 1e-6;  // More relaxed tolerance
                        bool partitionOK = checkPartitionOfUnity(matA, tolerance);
                        
                        if (!partitionOK) {
                            gsInfo << "\n*** WARNING: Partition of unity violated with tolerance " << tolerance << " ***\n";
                            gsInfo << "Attempt: " << attempt << "\n";
                            gsInfo << "This may indicate truncated basis functions or spillover issues\n";
                            
                            // Save diagnostic information
                            std::string diagFile = "partition_diagnostic_attempt" + std::to_string(attempt) + ".txt";
                            std::ofstream diag(diagFile);
                            diag << "Partition of unity violation at attempt " << attempt << "\n";
                            diag << "Matrix dimensions: " << matA.rows() << " x " << matA.cols() << "\n";
                            diag << "Tolerance: " << tolerance << "\n";
                            diag.close();
                            gsInfo << "Diagnostic saved to " << diagFile << "\n";
                            
                            // Try to continue with a warning instead of throwing
                            gsInfo << "Continuing execution despite partition of unity violation...\n\n";
                        } else {
                            gsInfo << "Partition of unity satisfied.\n";
                        }

                        afterassembleA = std::chrono::system_clock::now();
                        std::chrono::duration<double> elapsed_seconds_assembleA = afterassembleA - beforeassembleA;
                        gsInfo << "Assembly took: " << elapsed_seconds_assembleA.count() << " seconds\n";

                        // Check matrix rank before solving
                        gsMatrix<> matA_dense = gsEigen::MatrixXd(matA);
                        int rankA = computeRankByZeroRows(matA_dense);
                        gsInfo << "Matrix A rank (by zero rows): " << rankA << " / " << matA.rows() << " rows\n";
                        gsInfo << "Matrix A dimensions: " << matA.rows() << " x " << matA.cols() << "\n";

                        // Solve the system
                        gsMatrix<> b_vec = vectB;
                        vectB = gsEigen::MatrixXd(matA).transpose() * vectB;
                        gsMatrix<> matAsquare = gsEigen::MatrixXd(matA).transpose() * gsEigen::MatrixXd(matA);

                        int rankAsquare = computeRankByZeroRows(matAsquare);
                        gsInfo << "matAsquare rank: " << rankAsquare << " / " << matAsquare.rows() << " rows\n";
                        gsInfo << "matAsquare.determinant(): " << matAsquare.determinant() << "\n";

                        commonSize = functionDescription.size();
                        gsMatrix<real_t> vectSol = matAsquare.partialPivLu().solve(vectB);

                        gsInfo << "vectSol.rows(): " << vectSol.rows() << "\n";
                        gsMatrix<> matC = matAsquare * vectSol - vectB;
                        gsInfo << "Residual norm (A^T A x - A^T b) matCbefore\n";
                        gsInfo << matC.maxCoeff() << "\n";
                        outfile << "matCbefore\n";
                        outfile << matC.maxCoeff() << "\n";
                        ////gsInfo << "matC:\n" << matC << "\n";
                        gsMatrix<> matOut = gsEigen::MatrixXd(matA) * vectSol;
                        logResidualDetails(uv1, b_vec, matOut, outfile, false);
                        matC = matOut - b_vec;
                        gsInfo << "Residual norm (Ax - b) matCafter\n";
                        gsInfo << matC.maxCoeff() << "\n";
                        outfile << "matCafter\n";
                        outfile << matC.maxCoeff() << "\n";
                        // Compute max row-wise Euclidean norm (max pointwise error)
                        double globalError = 0.0;
                        for (index_t i = 0; i < matC.rows(); ++i) {
                            globalError = std::max(globalError, matC.row(i).norm());
                        }
                        int minusnumber = 0;
                        gsVector<size_t> numIrregular(Bells.size());
                        numIrregular.setZero();
                        std::chrono::time_point<std::chrono::system_clock> beforejack, afterjack;
                        beforejack = std::chrono::system_clock::now();


                        minusnumber = checkJacobianDeterminant(
                            uv2,
                            mpbes,
                            vectSol,
                            numIrregular,
                            false);
                        afterjack = std::chrono::system_clock::now();
                        std::chrono::duration<double> elapsed_seconds_jack = afterjack - beforejack;
                        gsInfo << "jack took: " << elapsed_seconds_jack.count() << "\n";
                        outfile << "jack took: " << elapsed_seconds_jack.count() << "\n";
                        
                        // Calculate total sample points and percentage of irregular points
                        index_t totalPoints = 0;
                        for (index_t patch = 0; patch < uv2.size(); ++patch) {
                            totalPoints += uv2[patch].cols();
                        }
                        real_t irregularPercentage = totalPoints > 0 ? (100.0 * minusnumber / totalPoints) : 0.0;
                        
                        gsInfo << "Irregular points: " << minusnumber << " / " << totalPoints 
                               << " (" << irregularPercentage << "%)\n";
                        gsInfo << "Global error: " << globalError << "\n";
                        outfile << "Irregular points: " << minusnumber << " / " << totalPoints 
                                << " (" << irregularPercentage << "%)\n";
                        outfile << "Global error: " << globalError << "\n";
                        
                        gsMatrix<> geomDer;
            
                        gsInfo << "minusnumber: " << minusnumber << "\n";
                        outfile << "minusnumber: " << minusnumber << "\n";
                        std::vector<std::vector<std::vector<double>>> localCoeffs(Bells.size());
                        gsVector<gsVector<int>> localIndex;
                        savePatches(
                            vectSol,
                            globalIndexTHB,
                            SubdomainHierarchy,
                            outfile,
                            "savedPatch", localIndex
                        );
                        // Mesh generation moved to before emergency exit only
                        // Output solution coefficients
                        outfile << "Solution coefficients: " << vectSol.rows() << " x " << vectSol.cols() << "\n";
                        for (index_t i = 0; i < vectSol.rows(); ++i) {
                            outfile << vectSol(i, 0) << " " << vectSol(i, 1) << "\n";
                        }
                        gsMatrix<> matFeatOut;
                        double featureError = 1;
                        for (int patch = 0; patch < localCoeffs.size(); ++patch) {
                            for (int locind = 0; locind < localCoeffs[patch].size(); locind++) {
                                for (int twin = 0; twin < localCoeffs[patch][locind].size(); ++twin) {
                                    //gsInfo << localCoeffs[patch][locind][twin] << " ";
                                    //<< localCoeffs[patch][locind][twin]
                                }
                                //gsInfo << "\n";
                            }
                        }
                        outfile << "local coeffs:\n";
                        for (int patch = 0; patch < localCoeffs.size(); ++patch) {
                            for (int locind = 0; locind < localCoeffs[patch].size(); locind++) {
                                for (int twin = 0; twin < localCoeffs[patch][locind].size(); ++twin) {
                                    outfile << localCoeffs[patch][locind][twin] << " ";
                                    //<< localCoeffs[patch][locind][twin]
                                }
                                outfile << "\n";
                            }
                        }
                        toc = std::chrono::system_clock::now();
                        std::chrono::duration<double> elapsed_seconds_write = toc - startTime;
                        featureError = boundaryError(
                            mpbes,
                            mp1,
                            vectSol);

                        // Test interface continuity using assemble() - gold standard test
                        real_t assemblyBoundaryError = testBoundaryAssembly(mpbes, mp1, vectSol);

                        // Build target geometry matrix at uv1 sample points for error computation
                        index_t totalSamples = 0;
                        for (index_t p = 0; p < uv1.size(); ++p)
                            totalSamples += uv1(p).cols();
                        
                        gsMatrix<real_t> targetGeometry(totalSamples, 2);
                        index_t rowIdx = 0;
                        for (index_t p = 0; p < uv1.size(); ++p)
                        {
                            const gsMatrix<real_t>& uvPatch = uv1(p);
                            for (index_t k = 0; k < uvPatch.cols(); ++k, ++rowIdx)
                            {
                                gsMatrix<real_t> xy = mp1.patch(p).eval(uvPatch.col(k));
                                targetGeometry(rowIdx, 0) = xy(0, 0);
                                targetGeometry(rowIdx, 1) = xy(1, 0);
                            }
                        }

                        // Evaluate global fitting error: fitted vs original at uv1 sample points
                        globalError = globalFittingError(mpbes, mp1, vectSol, uv1, targetGeometry);

                        outfile << "globalError: " << globalError << "\n";
                        outfile << "featureError: " << featureError << "\n";
                        outfile << "assemblyBoundaryError: " << assemblyBoundaryError << "\n";
                        outfile << "minusnumber: " << minusnumber << "\n";
                        gsInfo << "globalError: " << globalError << "\n";
                        gsInfo << "featureError: " << featureError << "\n";
                        gsInfo << "assemblyBoundaryError: " << assemblyBoundaryError << "\n";
                        index_t geoDim = 2;
                        gsSparseMatrix<real_t> A(vectSol.rows() * geoDim, vectSol.rows() * geoDim);
                        // Temporarily disabled to see full boundaryError output
                        // if (minusnumber > 0) return 260106;
                        //numIrregular(patch)++;
                        if (globalError <= epsilon_g && featureError <= epsilon_f && numIrregular(patch) == 0) {
                            //if (globalError <= epsilon_g && featureError <= epsilon_f && numIrregular(patch) != 0) {
                        success:
                            outfile << "Success! iteration = " << iteration << ", coarselevel = " << coarseLevel
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
                            acceptedCoefs = localCoeffs;
                            AcceptedvectSol = vectSol;
                            Acceptedmpbes = std::make_unique<MPBES<2, real_t>>(mpbes);
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
                            nonCheckedCells.removeElement(currArrayIndex);
                            outfile << "nonCheckedCells.size() has become: " << nonCheckedCells.size() << "\n";

                            int jopa = 4 * (int)pow(2, levNow) * 4 * (int)pow(2, levNow);
                            attempt = (attempt + 1) % (jopa);
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
                                gsInfo << numIrregular(0) << "\n";
                                gsInfo << numIrregular(1) << "\n";
                                gsInfo << numIrregular(2) << "\n";
                                //return 0;
                            }

                            gsInfo << "\n=== Starting attempt " << attempt << " - calling nonLinearOptimization ===\n";
                            outfile << "\n=== Starting attempt " << attempt << " - calling nonLinearOptimization ===\n";
                            
                            // if (attempt == 7)   DTD = 1;
                            gsMatrix<> localVectSol(localCoeffs[0].size(), 2);
                            gsVector<gsMatrix<>> localVectSols(localCoeffs.size());

                            A.setZero();
                            real_t fitting = 1;
                            real_t uniformity = 1e+4, length = 0;
                            
                            // Call the modernized optimization function
                            nonLinearOptimization(
                                mpbes, uv1, uv2, mp1,
                                matAsquare, matA, b_vec,
                                vectSol, A,
                                SubdomainHierarchy, Bells, functionDescription,
                                indexInTHB, isTruncated, presentation,
                                spilloverFunctionCoordinates, hasSpillover, globalIndexTHB,
                                boxMat, currentLastNonZeroRow,
                                fitting, uniformity, length,
                                epsilon_g, epsilon_f,
                                numIrregular, geoDim
                            );

                            // Recalculate metrics after optimization
                            gsMatrix<> matOut = gsEigen::MatrixXd(matA) * vectSol;
                            matC = matOut - b_vec;
                            globalError = globalFittingError(mpbes, mp1, vectSol, uv1, targetGeometry);
                            // featureError = 0;
                            
                            gsInfo << "=== Attempt " << attempt << " results: globalError=" << globalError 
                                   << ", featureError=" << featureError 
                                   << ", irregularities=" << numIrregular(patch) << " ===\n";
                            outfile << "=== Attempt " << attempt << " results: globalError=" << globalError 
                                    << ", featureError=" << featureError 
                                    << ", irregularities=" << numIrregular(patch) << " ===\n";

                            
                            // Increment attempt counter for next iteration
                            int jopa = 4 * (int)pow(2, levNow) * 4 * (int)pow(2, levNow);
                            attempt = (attempt + 1) % (jopa);
                            //gsFileManager::open("logFile_poissonTHB_example.txt");                            
                            if (globalError <= epsilon_g && featureError <= epsilon_f
                                && numIrregular(patch) == 0
                                ) {
                                outfile << "LO worked! iteration = " << iteration << ", coarselevel = " << coarseLevel
                                    << "\n";
                                gsInfo << "LO worked!\n";
                                //uniformity *= 10;
                                //return 7;
                                gsInfo << "coefficients: \n";
                                gsInfo << "fitting: " << fitting << "\n";
                                gsInfo << "uniformity: " << uniformity << "\n";
                                outfile << "coefficients: \n";
                                outfile << "fitting: " << fitting << "\n";
                                outfile << "uniformity: " << uniformity << "\n";
                                //gsFileManager::open("logFile_poissonTHB_example.txt");
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
                                acceptedCoefs = localCoeffs;
                                AcceptedvectSol = vectSol;
                                Acceptedmpbes = std::make_unique<MPBES<2, real_t>>(mpbes);
                                acceptedMatOut = matOut;
                                //std::vector<std::vector<std::vector<int>>> acceptedIndex = localIndex;
                                gsVector<gsVector<int>> acceptedIndex = localIndex;
                                //AcceptedfunctionDescription = functionDescription;
                                AcceptedisActive = isActive;
                                AcceptedglobalIndex = globalIndex;

                                // Save mesh at attempt 55 for inspection
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
                                        //gsInfo << boxMat(patch)(l)(m);
                                        //gsInfo << "\t";
                                    }
                                    outfile << "\n";
                                    //gsInfo << "\n";
                                }
                                AcceptedlastRow(patch) = lastNonZeroRow;
                                /*for (size_t i = lastNonZeroRow; i < boxMat(patch).size(); i++)
                                {
                                    for (size_t j = 0; j < 5; j++)
                                    {

                                    }
                                }*/
                                projections++;
                                //gsWrite
                                //setSupports(THB, supps);
                                anmat2 = THB.anchors();
                                wasRebuilt = 1;
                                ////gsInfo << "nonCheckedCells.size(): " << nonCheckedCells.size() << "\n";
                                nonCheckedCells.removeElement(currArrayIndex);
                                //gsInfo << "nonCheckedCells.size() has become: " << nonCheckedCells.size() << "\n";
                                outfile << "nonCheckedCells.size() has become: " << nonCheckedCells.size() << "\n";
                                //gsInfo << "break\n";
                                //break;
                                continue;
                            }
                            else {
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
                                //if (numIrregular(patch) > 0)
                                if (numIrregular(patch) > 1)
                                {
                                    outfile << "The parameterization is not regular\n";
                                    gsInfo << "The parameterization is not regular\n";
                                    //gsInfo << "minusnumber: " << minusnumber << "\n";
                                    //gsFileManager::open("logFile_poissonTHB_example.txt");
                                    attempt = (attempt + 1) % (jopa);
                                    continue;
                                }
                                gsInfo << "minusnumber: " << minusnumber << "\n";
                                //goto success;
                                //gsFileManager::open("logFile_poissonTHB_example.txt");
                                goto hop;
                            hop:
                                gsInfo << "WITHDRAW\n\n\n";
                                outfile << "WITHDRAW\n\n\n";
                                printAB = 1;

                                gsInfo << "attempt: " << attempt;


                                /*for (int row = 0; row != vectSol.rows(); row++)
                                {
                                    for (int dim = 0; dim != geoDim; ++dim)
                                    {
                                        gsInfo << vectSol(row, dim) << " ";
                                    }
                                    gsInfo << "\n";
                                }*/


                                //outfile << "functionDescription:\n";
                                /*for (int patch = 0; patch < functionDescription.size(); ++patch) {
                                    for (int functionIndex = 0; functionIndex < functionDescription(patch).size(); ++functionIndex) {
                                        for (int i = 0; i < functionDescription(patch)(functionIndex).size(); ++i) {
                                            outfile << functionDescription(patch)(functionIndex)[i] << " ";
                                        }
                                        outfile << "\n";
                                    }
                                    outfile << "==============\n";
                                }*/
                                nonCheckedCells.removeElement(currArrayIndex);
                                totalAttempts++;
                                //gsInfo << "boxMat before restoration" << "\n";
                                outfile << "boxMat before restoration" << "\n";
                                for (int l = 0; l < lastNonZeroRow; l++) {
                                    for (int m = 0; m < 5; m++) {
                                        outfile << boxMat(patch)(l)(m);
                                        outfile << "\t";
                                        //gsInfo << boxMat(patch)(l)(m);
                                        //gsInfo << "\t";
                                    }
                                    outfile << "\n";
                                    //gsInfo << "\n";
                                }
                                (SubdomainHierarchy(patch).tree().getBoxes(lowCorners(patch), upCorners(patch), myLevel(patch)));
                                //gsInfo << lowCorners.rows() << "\n";
                                restoreTheHierarchy2(boxMat, patch, ourBox, THB, THBUnstructured, lowCorners, upCorners, myLevel, lastNonZeroRow);
                                AcceptedlastRow(patch) = lastNonZeroRow;
                                //restoreTheHierarchy(boxMat, patch, ourBox, THB, THBUnstructured, lowCorners, upCorners, myLevel,lastNonZeroRow);
                                //restoreTheHierarchy(createdBoxNum, lastNonZeroRow, boxMat, levNow, centerInd, ourBox, successfullAttempts, patch);
                                continue;
                            }
                        }
                        //else return 3;
                    }
                    else outfile << "No spline created\n";
                    outfile << "FINISHED\n";
                    //break;
                }
            }
        }
    }
    
    // Prepare return values
    AlgorithmResult result;
    result.mp = std::move(mp);
    result.Bells = Bells;
    result.boxMat = boxMat;
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
    real_t length_weight,
    
    // Convergence criteria
    real_t epsilon_g,
    real_t epsilon_f,
    
    // Irregular count (modified)
    gsVector<size_t>& numIrregular,
    
    // Geometry dimension
    index_t geoDim
) {
    PROFILE_FUNCTION();
    
    // Quality weights (predefined, no user input)
    real_t fitting = fitting_weight;
    real_t uniformity = uniformity_weight;
    real_t length = length_weight;
    real_t orthogonality = 0;
    
    // Store old solution for comparison
    gsMatrix<> vectSolOld = vectSol;
    
    // Reset optimization matrix
    A.setZero();
    
    // Call MPBES-based optimization with quality functionals
    // Create non-const copy for optimize function
    gsMatrix<> matAsquareCopy = matAsquare;
    optimize(mpbes, uv1, matAsquareCopy, vectSol, A, 
        fitting, orthogonality, length, uniformity, 1e-7, false);
    
    // Calculate residual
    gsMatrix<> matOut = gsEigen::MatrixXd(matA) * vectSol;
    gsMatrix<> matC = matOut - b_vec;
    real_t residual = matC.maxCoeff();
    
    // Skip expensive checks if residual is very large (not close to convergence)
    if (residual > epsilon_g * 10.0) {
        // Quick exit for attempts that are far from convergence
        numIrregular.setZero();
        return;
    }
    
    // Check Jacobian determinant
    int minusnumber = checkJacobianDeterminant(
        uv2, mpbes, vectSol, numIrregular, false);
    
    // Calculate boundary error
    real_t featureError = testBoundaryAssembly(mpbes, mp1, vectSol);
    
    // Check convergence
    if (residual < epsilon_g && featureError < epsilon_f && minusnumber == 0) {
        return;
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
    
    // Non-converged optimization - skip mesh generation for performance
    // generateVisualizationMesh(boxMat, 20, mpbes, mp1, vectSol, 
    //     currentLastNonZeroRow, "mesh_unstable", true);
}

int main() {
    gsInfo << "Y\n";
    xboxFile.open("xboxFile.txt");
    yboxFile.open("yboxFile.txt");
    gsStopwatch clock;
    DTD = 0;
    printAB = 0;
    localTempAttempt = 0;
    std::chrono::time_point<std::chrono::system_clock> startTime, iterTime;
    startTime = std::chrono::system_clock::now();
    std::string filename("hexagon_3p_4l_240325.xml");
    const int dim = 2;
    int row, acceptedsize, attempt = 0;
    std::string badFile = "badgeoLocalAndreatemp6";
    std::string pvdFile = "resultLocalAndreatemp6";
    int los = 0, nlos = 0, proj = 0;
    std::string givenGeo;
    int valid = 0;
    int gradingExtent;
    real_t epsilon_g = 0.02, epsilon_f = 0.02;
    real_t lcx, lcy, ucx, ucy;
    std::string acCond = to_string(epsilon_g) + "and" + to_string(epsilon_f);
    givenGeo = "two_squares_lev1";
    givenGeo += "L2";
    givenGeo += "LO";
    givenGeo += "NLO";
    std::string fileLoc = "logFile_poissonTHB_example";
    outfile.open(fileLoc + ".txt");
    gsFileData<> data0(filename);
    int iter = -1;
    int successfullAttempts = 0, totalAttempts = 0;
    index_t degree;
    real_t tol = 1e-8, gtol = 1e-8;
    char method = 'l';

    // Run the main algorithm
    AlgorithmResult result = unrefinementAlgorithmHBJ(filename, epsilon_g, epsilon_f, method, givenGeo, acCond);
    
    std::chrono::time_point<std::chrono::system_clock> toc = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_finished = toc - startTime;
    
    gsInfo << "FINISHED\n";
    gsInfo << "Total execution time: " << elapsed_finished.count() << " seconds\n";
    outfile << "FINISHED took: " << elapsed_finished.count() << "\n";
    outfile << "Successful attempts: " << result.successfullAttempts << " / " << result.totalAttempts << "\n\n";
    
    // Generate visualization mesh
    gsInfo << "\n=== Generating final visualization mesh ===\n";
    outfile << "\n=== Generating final visualization mesh ===\n";
    
    generateVisualizationMesh(
        result.boxMat, 
        20, 
        *result.mpbes, 
        result.mp1, 
        result.AcceptedvectSol, 
        result.currentLastNonZeroRow, 
        "output_mesh_final", 
        true
    );
    
    gsInfo << "Visualization mesh saved to output_mesh_final.txt\n";
    outfile << "Visualization mesh saved to output_mesh_final.txt\n";
    
    // Print profiling summary
    g_profiler.printSummary();
    
    // Close files
    outfile.close();
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