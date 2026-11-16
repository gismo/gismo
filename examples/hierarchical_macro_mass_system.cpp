/** @file hierarchical_macro_mass_system.cpp

    @brief Builds local mass moment systems for each atomic cell of a
    hierarchical domain partition and estimates the numerical rank of the
    local product space T_r^M = span{B_i * B_j|_{D_r}}.

    For every gsHTree leaf (= one atomic cell = one macro region):
      - Computes reference moments  b_{ij} = integral_{D_r} B_i * B_j
      - Builds the sample matrix    S_{aq} = B_i(x_q) * B_j(x_q)
      - Estimates the numerical rank of S via truncated SVD
      - Reports estimated free quadrature points ~ ceil(rank / 3)

    The ratio sum(ceil(rank_r / 3)) / sum(Gauss pts) is the first viability
    test: if it is << 1, a moment-fitted macro rule can beat standard Gauss.

    Usage:
        ./hierarchical_macro_mass_system [--plot] [--rankSome] [-c N]

    Options:
        --rankSome   Compute actual SVD rank for one representative cell per level
        -c N         Compute actual SVD rank for cell N specifically

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst
*/

#include <gismo.h>
#include <unordered_map>
#include <set>

using namespace gismo;

// ============================================================
// Atomic cell data structures (same as hierarchical_domain_partition.cpp)
// ============================================================

struct ElementBox
{
    gsVector<real_t> lower;
    gsVector<real_t> upper;
};

struct AtomicCell
{
    index_t            level;
    gsVector<index_t>  lowerIndex;
    gsVector<index_t>  upperIndex;
    gsVector<real_t>   lowerParam;
    gsVector<real_t>   upperParam;
    std::vector<ElementBox> elements;
};

/**
 * @brief Collect all gsHTree leaves as AtomicCells.
 * @note Requires !thb.manualLevels().
 */
template<short_t d, class T>
std::vector<AtomicCell> collectAtomicCells(const gsHTensorBasis<d,T> & thb)
{
    GISMO_ASSERT( !thb.manualLevels(),
        "collectAtomicCells: manualLevels bases require a dyadic-to-knot-index "
        "conversion not implemented here." );

    std::vector<AtomicCell> cells;
    auto leaf = thb.tree().beginLeafIterator();

    for (; leaf.good(); leaf.next())
    {
        AtomicCell cell;
        cell.level      = leaf.level();
        cell.lowerIndex = leaf.lowerCorner();
        cell.upperIndex = leaf.upperCorner();

        cell.lowerParam.resize(d);
        cell.upperParam.resize(d);

        for (short_t dd = 0; dd < d; ++dd)
        {
            const gsKnotVector<T> & kv =
                thb.tensorLevel(cell.level).component(dd).knots();
            cell.lowerParam[dd] = kv( cell.lowerIndex[dd] );
            cell.upperParam[dd] = kv( cell.upperIndex[dd] );
        }

        gsVector<index_t> alpha(d);
        alpha = cell.lowerIndex;

        bool more = true;
        while (more)
        {
            ElementBox elem;
            elem.lower.resize(d);
            elem.upper.resize(d);
            for (short_t dd = 0; dd < d; ++dd)
            {
                const gsKnotVector<T> & kv =
                    thb.tensorLevel(cell.level).component(dd).knots();
                elem.lower[dd] = kv( alpha[dd]     );
                elem.upper[dd] = kv( alpha[dd] + 1 );
            }
            cell.elements.push_back(elem);

            more = false;
            for (short_t dd = d - 1; dd >= 0; --dd)
            {
                if (alpha[dd] < cell.upperIndex[dd] - 1)
                {
                    ++alpha[dd];
                    more = true;
                    break;
                }
                alpha[dd] = cell.lowerIndex[dd];
            }
        }

        cells.push_back(cell);
    }

    return cells;
}

// ============================================================
// Local mass moment system
// ============================================================

struct LocalMassSystem
{
    index_t macroId;
    index_t level;
    index_t numElements;

    // Ordered list of (unique) product pairs (i,j) with i <= j
    std::vector< std::pair<index_t,index_t> > pairs;

    // Reference mass moments b[row] = integral_{D_r} B_i * B_j
    gsVector<real_t> b;

    // Rank upper bounds of the local product space T_r^M = span{B_i * B_j|_{D_r}}.
    //
    // oldRankBound      = min(numPairs, (n_x + 2p)(n_y + 2p))
    //   — dimension of the degree-2p *max-continuity* (C^{2p-1}) tensor spline
    //   space.  This over-estimates the target: product B_i*B_j has continuity
    //   min(p-1, p-1) = p-1, not 2p-1.
    //
    // correctedRankBound = min(numPairs, ((p+1)n_x + p)((p+1)n_y + p))
    //   — dimension of the degree-2p, C^{p-1} tensor spline space.  In 1D with
    //   n knot spans the dimension is (2p+1) + (n-1)*(p+1) = (p+1)*n + p.
    //   This is the correct conservative bound for the product space.
    //
    // actualRank = -1 means "not computed" (requires SVD, only on request).
    index_t oldRankBound;
    index_t correctedRankBound;
    index_t actualRank;
};

/**
 * @brief Build the local mass moment system for one atomic cell.
 *
 * Uses `gsFunctionSet::compute` with `SAME_ELEMENT` per element so that the
 * HB-tree level-detection runs only once per element (not once per Gauss point),
 * active indices are cached as a single column, and pair-row lookups happen only
 * once per element.  Pair storage uses `std::unordered_map` for O(1) access.
 *
 * @tparam d        domain dimension
 * @tparam T        scalar type
 * @param cell      the atomic cell
 * @param thb       the hierarchical basis
 * @param macroId   index of this cell in the global cell list
 * @param overQuad  extra quadrature points per direction above degree
 */
template<short_t d, class T>
LocalMassSystem buildLocalMassSystem(
    const AtomicCell         & cell,
    const gsHTensorBasis<d,T>& thb,
    index_t                    macroId,
    index_t                    overQuad = 4)
{
    LocalMassSystem sys;
    sys.macroId     = macroId;
    sys.level       = cell.level;
    sys.numElements = (index_t)cell.elements.size();

    // High-order Gauss-Legendre rule at the leaf's level
    gsOptionList opts;
    opts.addInt   ("quRule", "Quadrature rule", gsQuadrature::GaussLegendre);
    opts.addReal  ("quA",    "Points: quA*deg + quB", T(1));
    opts.addInt   ("quB",    "Points: quA*deg + quB", (int)overQuad);
    opts.addSwitch("overInt","Over-integration", false);
    typename gsQuadRule<T>::uPtr rule = gsQuadrature::getPtr(thb.tensorLevel(cell.level), opts);

    // O2: flat 64-bit hash for (i,j) pairs — O(1) average vs std::map O(log N)
    struct PairHash {
        size_t operator()(std::pair<index_t,index_t> p) const noexcept {
            return std::hash<size_t>()(
                (size_t(uint32_t(p.first)) << 32) | size_t(uint32_t(p.second)));
        }
    };
    std::unordered_map<std::pair<index_t,index_t>, index_t, PairHash> rowOfPair;
    rowOfPair.reserve(256);
    std::vector<T> bVec;

    // O1+O3: For an atomic cell at level l all points lie in Ω^l \ Ω^{l+1}.
    // The truncation correction of any THB function vanishes there (level-(l+1)
    // functions are zero outside Ω^{l+1}), so THB values = tensor-level-l values.
    // Using the plain tensor basis avoids all HB-tree overhead.
    // With SAME_ELEMENT the active lookup and level detection run only once
    // per element (for col 0); the allActives() column is reused for all q.
    // Note: allActives() returns tensor-level-l indices — sufficient as unique
    //       pair identifiers within one cell (all actives share the same level).
    const gsTensorBSplineBasis<d,T>& tensorBasis =
        dynamic_cast<const gsTensorBSplineBasis<d,T>&>(thb.tensorLevel(cell.level));

    gsFuncData<T> fdata;
    fdata.flags = NEED_VALUE | NEED_ACTIVE | SAME_ELEMENT;

    // Pre-allocate once per cell (avoids per-element heap allocation)
    gsMatrix<T> nodes;
    gsVector<T> weights;

    for (const ElementBox & e : cell.elements)
    {
        rule->mapTo(e.lower, e.upper, nodes, weights);
        tensorBasis.compute(nodes, fdata);

        const index_t na = fdata.allActives().rows();
        const index_t nq = (index_t)weights.size();

        // Discover pairs once (col 0), cache row indices
        std::vector<index_t> rowCache;
        rowCache.reserve((na * (na + 1)) / 2);

        for (index_t a = 0; a < na; ++a)
        {
            const index_t bi = fdata.allActives()(a, 0);
            for (index_t c = a; c < na; ++c)
            {
                const index_t bj  = fdata.allActives()(c, 0);
                const index_t ii  = std::min(bi, bj);
                const index_t jj  = std::max(bi, bj);
                auto          key = std::make_pair(ii, jj);

                auto it = rowOfPair.find(key);
                index_t row;
                if (it == rowOfPair.end())
                {
                    row = (index_t)rowOfPair.size();
                    rowOfPair[key] = row;
                    bVec.push_back(T(0));
                }
                else
                    row = it->second;

                rowCache.push_back(row);
            }
        }

        // Accumulate bVec for all q — no map access
        for (index_t q = 0; q < nq; ++q)
        {
            const T w = weights[q];
            index_t k = 0;
            for (index_t a = 0; a < na; ++a)
                for (index_t c = a; c < na; ++c, ++k)
                    bVec[rowCache[k]] +=
                        w * fdata.allValues()(a, q) * fdata.allValues()(c, q);
        }
    }

    // Store pairs and b
    const index_t numPairs = (index_t)rowOfPair.size();
    sys.pairs.resize(numPairs);
    for (const auto & kv : rowOfPair)
        sys.pairs[kv.second] = kv.first;

    sys.b.resize(numPairs);
    for (index_t r = 0; r < numPairs; ++r)
        sys.b[r] = bVec[r];

    // Rank upper bounds
    const index_t p   = (index_t)thb.degree(0);
    const index_t n_x = cell.upperIndex[0] - cell.lowerIndex[0];
    const index_t n_y = cell.upperIndex[1] - cell.lowerIndex[1];

    sys.oldRankBound = std::min(numPairs, (n_x + 2*p) * (n_y + 2*p));

    // Degree-2p, C^{p-1}: 1D dim = (p+1)*n + p
    const index_t dim1d_x = (p + 1) * n_x + p;
    const index_t dim1d_y = (p + 1) * n_y + p;
    sys.correctedRankBound = std::min(numPairs, dim1d_x * dim1d_y);

    sys.actualRank = -1; // computed on request

    return sys;
}

// ============================================================
// Actual numerical rank via sample matrix SVD (optional, slow for large cells)
// ============================================================

/**
 * @brief Compute the actual numerical rank of the local product space T_r^M.
 *
 * Builds S(numPairs × nSamples) where each column is [B_i(x_q)*B_j(x_q)]
 * at a Gauss sample point.  Subsamples columns uniformly when totalQ > maxCols
 * to keep the SVD tractable.  Uses BDCSVD with threshold 1e-10 * s_max.
 *
 * @param maxCols  cap on number of sample columns (default 5000)
 */
template<short_t d, class T>
index_t computeSampleRank(
    const AtomicCell         & cell,
    const gsHTensorBasis<d,T>& thb,
    const LocalMassSystem    & sys,
    index_t                    overQuad = 4,
    index_t                    maxCols  = 5000)
{
    const index_t numPairs = (index_t)sys.pairs.size();
    if (numPairs == 0) return 0;

    // Build an inverse map: pair (i,j) -> row index
    struct PairHash {
        size_t operator()(std::pair<index_t,index_t> p) const noexcept {
            return std::hash<size_t>()(
                (size_t(uint32_t(p.first)) << 32) | size_t(uint32_t(p.second)));
        }
    };
    std::unordered_map<std::pair<index_t,index_t>, index_t, PairHash> rowOfPair;
    rowOfPair.reserve((size_t)numPairs * 2);
    for (index_t r = 0; r < numPairs; ++r)
        rowOfPair[sys.pairs[r]] = r;

    // Gather all Gauss columns (or a uniform subsample)
    gsOptionList opts;
    opts.addInt   ("quRule", "Quadrature rule", gsQuadrature::GaussLegendre);
    opts.addReal  ("quA",    "Points: quA*deg + quB", T(1));
    opts.addInt   ("quB",    "Points: quA*deg + quB", (int)overQuad);
    opts.addSwitch("overInt","Over-integration", false);
    typename gsQuadRule<T>::uPtr rule =
        gsQuadrature::getPtr(thb.tensorLevel(cell.level), opts);

    // Count total Gauss points first
    index_t totalQ = 0;
    for (const ElementBox & e : cell.elements)
    {
        gsMatrix<T> nodes; gsVector<T> weights;
        rule->mapTo(e.lower, e.upper, nodes, weights);
        totalQ += (index_t)weights.size();
    }

    const index_t nSamples = std::min(totalQ, maxCols);
    const index_t stride   = (nSamples < totalQ) ?
                              (totalQ / nSamples) : 1;

    // Use plain Eigen matrix — gsMatrix<T> incompatible with gsEigen::BDCSVD
    // (missing NestedExpression typedef in G+Smo's gsMatrix wrapper).
    typedef gsEigen::Matrix<T, gsEigen::Dynamic, gsEigen::Dynamic> EigenMat;
    EigenMat S = EigenMat::Zero(numPairs, nSamples);

    const gsTensorBSplineBasis<d,T>& tensorBasis =
        dynamic_cast<const gsTensorBSplineBasis<d,T>&>(
            thb.tensorLevel(cell.level));

    gsFuncData<T> fdata;
    fdata.flags = NEED_VALUE | NEED_ACTIVE | SAME_ELEMENT;

    gsMatrix<T> nodes; gsVector<T> weights;
    index_t globalQ = 0;  // index over all Gauss points
    index_t colIdx  = 0;  // index into S columns

    for (const ElementBox & e : cell.elements)
    {
        rule->mapTo(e.lower, e.upper, nodes, weights);
        tensorBasis.compute(nodes, fdata);

        const index_t na = fdata.allActives().rows();
        const index_t nq = (index_t)weights.size();

        for (index_t q = 0; q < nq; ++q, ++globalQ)
        {
            if (colIdx >= nSamples) break;
            if ((globalQ % stride) != 0) continue;

            for (index_t a = 0; a < na; ++a)
            {
                const index_t bi = fdata.allActives()(a, 0);
                for (index_t c = a; c < na; ++c)
                {
                    const index_t bj  = fdata.allActives()(c, 0);
                    const index_t ii  = std::min(bi, bj);
                    const index_t jj  = std::max(bi, bj);
                    auto          it  = rowOfPair.find({ii, jj});
                    if (it != rowOfPair.end())
                        S(it->second, colIdx) +=
                            fdata.allValues()(a, q) * fdata.allValues()(c, q);
                }
            }
            ++colIdx;
        }
        if (colIdx >= nSamples) break;
    }

    // Trim S to actual columns filled
    S.conservativeResize(numPairs, colIdx);

    // BDCSVD rank
    gsEigen::BDCSVD<EigenMat> svd(S);
    const gsEigen::Matrix<T, gsEigen::Dynamic, 1>& sv = svd.singularValues();
    if (sv.size() == 0) return 0;

    const T threshold = T(1e-10) * sv[0];
    index_t rank = 0;
    for (index_t k = 0; k < sv.size(); ++k)
        if (sv[k] > threshold) ++rank;

    return rank;
}

// ============================================================
// main
// ============================================================

int main(int argc, char * argv[])
{
    bool   plot     = false;
    bool   rankSome = false;
    index_t overQ   = 4;   // quB: extra points per direction above degree
    index_t rankCell = -1; // specific cell for SVD (-1 = none)

    gsCmdLine cmd("Local mass moment system diagnostics for hierarchical domains.");
    cmd.addSwitch("plot",     "Write ParaView mesh file",                    plot);
    cmd.addSwitch("rankSome", "Compute actual SVD rank for one cell/level",  rankSome);
    cmd.addInt   ("o", "overQ",    "Extra Gauss points per direction",       overQ);
    cmd.addInt   ("c", "rankCell", "Compute actual SVD rank for cell N",     rankCell);
    try { cmd.getValues(argc, argv); } catch (int rv) { return rv; }

    // ------------------------------------------------------------------
    // Build THB-spline basis with L-shaped local refinement
    // ------------------------------------------------------------------
    gsKnotVector<> kv1(0, 1.0, 0, 3, 1);
    kv1.uniformRefine();
    kv1.uniformRefine();
    gsKnotVector<> kv2 = kv1;

    gsTensorBSplineBasis<2, real_t> tbsb(kv1, kv2);
    gsTHBSplineBasis<2, real_t>     thb(tbsb);

    {
        thb.uniformRefine();
        for (index_t lvl = 1; lvl <= 5; ++lvl)
        {
            const index_t N = 4 * (1 << lvl);
            std::vector<index_t> boxes;
            boxes.insert(boxes.end(), {lvl, 0, 0, N, 4});
            boxes.insert(boxes.end(), {lvl, 0, 0, 4, N});
            thb.refineElements(boxes);
        }
    }

    const index_t deg = thb.degree(0); // assumed isotropic
    gsInfo << "Basis: degree=" << deg
           << ", maxLevel=" << thb.maxLevel()
           << ", DOF=" << thb.size() << "\n\n";

    if (plot)
    {
        gsMesh<> mesh(thb);
        gsWriteParaview(mesh, "hierarchical_macro_mesh");
    }

    // ------------------------------------------------------------------
    // Step 1: collect atomic cells
    // ------------------------------------------------------------------
    const std::vector<AtomicCell> cells = collectAtomicCells<2, real_t>(thb);
    gsInfo << "Atomic cells: " << cells.size() << "\n\n";

    // ------------------------------------------------------------------
    // Step 2: build local mass systems and print per-cell diagnostics
    // ------------------------------------------------------------------

    // Determine which cells get actual SVD rank
    std::set<index_t> svdCells;
    if (rankCell >= 0 && rankCell < (index_t)cells.size())
        svdCells.insert(rankCell);

    const index_t gaussPtsPerElem = (deg + 1) * (deg + 1); // standard Gauss

    gsInfo << std::setw(6)  << "macro"
           << std::setw(6)  << "level"
           << std::setw(8)  << "elems"
           << std::setw(8)  << "pairs"
           << std::setw(9)  << "old_ub"
           << std::setw(12) << "corrected"
           << std::setw(9)  << "actual"
           << std::setw(10) << "ceil(r/3)"
           << std::setw(12) << "Gauss pts"
           << std::setw(8)  << "ratio"
           << "\n";
    gsInfo << std::string(88, '-') << "\n";

    index_t sumOld       = 0;
    index_t sumCorrected = 0;
    index_t sumActual    = 0;
    index_t sumFree      = 0;
    index_t sumGauss     = 0;
    index_t sumElems     = 0;
    index_t numActual    = 0;

    // If --rankSome: pick the smallest cell (by numPairs) per level
    // We need to build systems first, then compute SVD for selected cells.
    std::vector<LocalMassSystem> systems;
    systems.reserve(cells.size());

    for (index_t r = 0; r < (index_t)cells.size(); ++r)
        systems.push_back(buildLocalMassSystem<2, real_t>(cells[r], thb, r, overQ));

    // Auto-select one cell per level for --rankSome
    if (rankSome)
    {
        std::map<index_t, index_t> bestPerLevel; // level -> cell with min pairs
        for (index_t r = 0; r < (index_t)systems.size(); ++r)
        {
            const index_t lv = systems[r].level;
            auto it = bestPerLevel.find(lv);
            if (it == bestPerLevel.end() ||
                (index_t)systems[r].pairs.size() <
                (index_t)systems[it->second].pairs.size())
                bestPerLevel[lv] = r;
        }
        for (const auto & kv : bestPerLevel)
            svdCells.insert(kv.second);
    }

    // Compute actual ranks for selected cells
    for (index_t r : svdCells)
    {
        gsInfo << "  Computing SVD rank for cell " << r
               << " (pairs=" << systems[r].pairs.size() << ") ...\n";
        systems[r].actualRank =
            computeSampleRank<2, real_t>(cells[r], thb, systems[r], overQ);
    }
    if (!svdCells.empty()) gsInfo << "\n";

    for (index_t r = 0; r < (index_t)cells.size(); ++r)
    {
        const LocalMassSystem & sys = systems[r];

        // Use correctedRankBound for free-point estimate; fall back to actual rank if available
        const index_t rankForFree = (sys.actualRank >= 0) ?
                                     sys.actualRank : sys.correctedRankBound;
        const index_t nFree  = (rankForFree + 2) / 3; // ceil(rank/3)
        const index_t nGauss = sys.numElements * gaussPtsPerElem;

        gsInfo << std::setw(6)  << r
               << std::setw(6)  << sys.level
               << std::setw(8)  << sys.numElements
               << std::setw(8)  << (index_t)sys.pairs.size()
               << std::setw(9)  << sys.oldRankBound
               << std::setw(12) << sys.correctedRankBound;

        if (sys.actualRank >= 0)
            gsInfo << std::setw(9) << sys.actualRank;
        else
            gsInfo << std::setw(9) << "---";

        gsInfo << std::setw(10) << nFree
               << std::setw(12) << nGauss
               << std::setw(8)  << std::fixed << std::setprecision(3)
               << (nGauss > 0 ? (real_t)nFree / (real_t)nGauss : 0.0)
               << "\n";

        sumOld       += sys.oldRankBound;
        sumCorrected += sys.correctedRankBound;
        sumFree      += nFree;
        sumGauss     += nGauss;
        sumElems     += sys.numElements;
        if (sys.actualRank >= 0) { sumActual += sys.actualRank; ++numActual; }
    }

    // ------------------------------------------------------------------
    // Global summary
    // ------------------------------------------------------------------
    gsInfo << "\n";
    gsInfo << "Global summary\n";
    gsInfo << "==============\n";
    gsInfo << "  Total atomic cells               : " << cells.size()  << "\n";
    gsInfo << "  Total elements                   : " << sumElems      << "\n";
    gsInfo << "  Sum of old rank bounds (n+2p)^2  : " << sumOld        << "\n";
    gsInfo << "  Sum of corrected bounds ((p+1)n+p)^2: " << sumCorrected << "\n";
    if (numActual > 0)
        gsInfo << "  Sum of actual ranks (" << numActual << " cells)      : " << sumActual << "\n";
    gsInfo << "  Sum of ceil(rank/3) [corrected]  : " << sumFree       << "\n";
    gsInfo << "  Total Gauss pts (p+1)^2 * elems  : " << sumGauss      << "\n";
    gsInfo << "  Ratio free/Gauss [old bound]     : "
           << std::fixed << std::setprecision(4)
           << (sumGauss > 0 ? (real_t)(sumOld+2)/3 / (real_t)sumGauss : 0.0)
           << "\n";
    gsInfo << "  Ratio free/Gauss [corrected]     : "
           << std::fixed << std::setprecision(4)
           << (sumGauss > 0 ? (real_t)sumFree / (real_t)sumGauss : 0.0)
           << "\n";
    gsInfo << "  (ceil(r/3) column uses corrected bound; falls back to actual rank if computed)\n\n";

    const real_t ratio = (sumGauss > 0) ? (real_t)sumFree / (real_t)sumGauss : 1.0;
    if (ratio < 1.0)
        gsInfo << ">>> Viability: PROMISING — estimated free pts < Gauss pts.\n";
    else
        gsInfo << ">>> Viability: NOT PROMISING — estimated free pts >= Gauss pts.\n";

    return 0;
}
