/** @file hierarchical_partition_study.cpp

    @brief Partition study for gsHPatchRule: how many quadrature nodes does the
    patch rule need when the active region of each level of a THB-spline mesh
    is cut into axis-aligned boxes in different ways?

    The active region of every hierarchical level is partitioned into boxes in
    three ways ("candidates"):
      - kdtree   : the leaves of the hierarchical kd-tree, as gsHPatchRule uses
                   them;
      - pairwise : repeated merging of same-level boxes whose union is a box;
      - sweep    : a canonical slab decomposition followed by a merge, with the
                   axis order chosen per level.

    For every (form, candidate) the driver reports the number of boxes, the
    exact number of nodes of the patch rule on those boxes, a lower bound that
    holds for any exact partition, and the partition time. The node count comes
    from a model that mirrors gsHPatchRule::_build; the model is validated
    against the real rule on the kdtree candidate and then reused unchanged for
    the other two. Forms are the ones of hierarchical_assembly_hpatch: mass
    (target degree 2p, regularity r), stiffness (2p, r-1) and biharmonic
    (2p, r-2, only for r >= 2).

    Background. The node count of a box expands into terms proportional to its
    volume, its cut area and the number of boxes, so the candidates trade the
    number of boxes against the length of the cuts. The minimum-box partition
    of a 2D rectilinear region is polynomial (Eppstein, arXiv:0908.3916); in 3D
    it is reported NP-complete (Dielissen and Kaldewaij, Inf. Process. Lett.
    38, 1991). The candidates are therefore heuristics, compared by measured
    node count. Every partition is checked to be exact (no overlap, no gap, no
    cell outside the region); the driver aborts otherwise.

    Output: a console table (lines starting with "row ") and, with --csv, a CSV
    file with the header
    form,dim,p,r,levels,pattern,admissible,jump,nElem,maxLevel,candidate,
    axisOrder,nBoxes,nNodes,nNodesLowerBound,nNodesRule,nFallbackRule,
    modelMatchesRule,tPartition,reps,threads,host,cpu,baseRefs,band
    (one line, no spaces), and one row per (form, candidate) in the order
    kdtree, pairwise, sweep. nNodesRule, nFallbackRule and modelMatchesRule are
    only filled for kdtree, axisOrder only for sweep. tPartition is in seconds.

    Usage:
        ./hierarchical_partition_study [-P deg] [-R reg] [-l levels] [-d dim]
            [--meshPattern lshape|center|diagonal] [--admissible]
            [--admissibleJump m] [--shuffle seed] [--baseRefs k] [--band scaled|thin] [--reps n] [--csv file]

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst
*/

#include <gismo.h>
#include <gsAssembler/gsHPatchRule.h>
#include <gsHSplines/gsHElement.h>
#include <gsHSplines/gsHElementHelper.h>
#include <gsParallel/gsOpenMP.h>
#include <gsCore/gsSysInfo.h>

#include <algorithm>
#include <array>
#include <fstream>
#include <map>
#include <numeric>
#include <random>
#include <sstream>
#include <unistd.h>

using namespace gismo;

namespace
{

enum MeshPattern { PatternLShape, PatternCenter, PatternDiagonal };

// Width of the refinement bands: BandScaled keeps the refined region's physical
// size for every --baseRefs >= 1 (the centre box is shifted to an even origin,
// which moves it by one level-targetLevel element at --baseRefs 1); BandThin
// keeps the bands 4 target-level elements wide.
enum BandVariant { BandScaled, BandThin };

// ============================================================
// Mesh construction (same meshes as hierarchical_assembly_hpatch)
// ============================================================

bool parseMeshPattern(const std::string & pattern, MeshPattern & result)
{
    if (pattern == "lshape")   { result = PatternLShape;   return true; }
    if (pattern == "center")   { result = PatternCenter;   return true; }
    if (pattern == "diagonal") { result = PatternDiagonal; return true; }
    return false;
}

bool parseBandVariant(const std::string & band, BandVariant & result)
{
    if (band == "scaled")
    {
        result = BandScaled;
        return true;
    }
    if (band == "thin")
    {
        result = BandThin;
        return true;
    }
    return false;
}

const char * bandName(BandVariant band)
{
    return band == BandScaled ? "scaled" : "thin";
}

// Appends one refinement box [targetLevel, low(d), upp(d)] to the flat encoding
// read by gsHTensorBasis::refineElements (stride 2d+1).
template<short_t d>
void appendTargetBox(const typename gsHTensorBasis<d,real_t>::point & low,
                     const typename gsHTensorBasis<d,real_t>::point & upp,
                     index_t targetLevel,
                     std::vector<index_t> & boxes)
{
    boxes.push_back(targetLevel);
    for (short_t dd = 0; dd < d; ++dd)
        boxes.push_back(low[dd]);
    for (short_t dd = 0; dd < d; ++dd)
        boxes.push_back(upp[dd]);
}

// Refinement boxes of one level for the chosen pattern, in level-targetLevel
// element indices; the pattern spans the first 2^baseRefs * 2^targetLevel of the
// 2^(baseRefs+1+targetLevel) elements per direction.
template<short_t d>
std::vector<index_t> makePatternBoxes(index_t targetLevel, MeshPattern pattern,
                                      index_t baseRefs, BandVariant band)
{
    // K base elements per direction; N is the pattern extent in level-targetLevel
    // element indices. Scaled bands are K wide (fixed physical size), thin
    // bands stay 4 elements wide. baseRefs >= 1 makes K, N and both widths even,
    // so every box corner is a multiple of 2 and aligns with parent cells; the
    // centre origin is rounded down to an even index for the same reason.
    const index_t K = index_t(1) << baseRefs;
    const index_t N = K * (index_t(1) << targetLevel);
    const index_t width = (band == BandScaled ? K : std::min<index_t>(4, N));
    std::vector<index_t> boxes;
    typename gsHTensorBasis<d,real_t>::point low, upp;

    switch (pattern)
    {
        case PatternLShape:
            boxes.reserve(d * (2 * d + 1));
            for (short_t DD = 0; DD < d; ++DD)
            {
                low.setZero();
                for (short_t dd = 0; dd < d; ++dd)
                    upp[dd] = (dd == DD ? N : width);
                appendTargetBox<d>(low, upp, targetLevel, boxes);
            }
            break;

        case PatternCenter:
            boxes.reserve(2 * d + 1);
            for (short_t dd = 0; dd < d; ++dd)
            {
                low[dd] = ((N - width) / 4) * 2;
                upp[dd] = low[dd] + width;
            }
            appendTargetBox<d>(low, upp, targetLevel, boxes);
            break;

        case PatternDiagonal:
            boxes.reserve((N / width) * (2 * d + 1));
            for (index_t start = 0; start < N; start += width)
            {
                for (short_t dd = 0; dd < d; ++dd)
                {
                    low[dd] = start;
                    upp[dd] = start + width;
                }
                appendTargetBox<d>(low, upp, targetLevel, boxes);
            }
            break;
    }

    return boxes;
}

// Marks the elements of the pattern boxes (parent-level elements) so that the
// admissible closure can be applied to them.
template<short_t d, class T>
typename gsHElementHelper<d,T>::HElementContainer
makeMarkedElementsFromBoxes(const std::vector<index_t> & boxes,
                            index_t targetLevel,
                            const gsHElementHelper<d,T> & helper)
{
    typedef gsHElementHelper<d,T> helper_t;
    typedef typename helper_t::point_t point_t;

    GISMO_ENSURE(targetLevel > 0, "Target refinement level must be positive.");
    GISMO_ENSURE(boxes.size() % (2 * d + 1) == 0, "Invalid refinement box encoding.");

    typename helper_t::HElementContainer marked;
    const index_t parentLevel = targetLevel - 1;
    const size_t  stride = 2 * d + 1;

    for (size_t i = 0; i < boxes.size() / stride; ++i)
    {
        const size_t offset = i * stride;
        point_t low, upp;
        for (short_t dd = 0; dd < d; ++dd)
        {
            GISMO_ENSURE(boxes[offset + 1 + dd] % 2 == 0 &&
                         boxes[offset + 1 + d + dd] % 2 == 0,
                         "Pattern boxes must align with parent elements.");
            low[dd] = boxes[offset + 1 + dd] / 2;
            upp[dd] = boxes[offset + 1 + d + dd] / 2;
        }

        typename helper_t::HElementContainer local = helper.explode(low, upp, parentLevel);
        marked.insert(local.begin(), local.end());
    }

    return marked;
}

// Permutes whole refinement tuples (stride 2d+1). The kd-tree is built by
// inserting the tuples one by one, so the insertion order is the "history" of
// the tree; shuffling it probes whether a partition depends on that history.
// A null generator leaves the order unchanged.
template<short_t d>
void shuffleTuples(std::vector<index_t> & boxes, std::mt19937 * rng)
{
    if (rng == nullptr) return;
    const size_t stride = 2 * d + 1;
    std::vector<size_t> perm(boxes.size() / stride);
    std::iota(perm.begin(), perm.end(), size_t(0));
    std::shuffle(perm.begin(), perm.end(), *rng);

    std::vector<index_t> shuffled;
    shuffled.reserve(boxes.size());
    for (size_t k = 0; k < perm.size(); ++k)
        shuffled.insert(shuffled.end(),
                        boxes.begin() + perm[k] * stride,
                        boxes.begin() + (perm[k] + 1) * stride);
    boxes.swap(shuffled);
}

template<short_t d, class T>
void refineLevel(gsTHBSplineBasis<d,T> & thb,
                 MeshPattern pattern,
                 index_t baseRefs,
                 BandVariant band,
                 index_t targetLevel,
                 bool admissible,
                 index_t admissibleJump,
                 std::mt19937 * rng)
{
    std::vector<index_t> boxes = makePatternBoxes<d>(targetLevel, pattern, baseRefs, band);
    if (!admissible)
    {
        shuffleTuples<d>(boxes, rng);
        thb.refineElements(boxes);
        return;
    }

    gsHElementHelper<d,T> helper(thb);
    typename gsHElementHelper<d,T>::HElementContainer marked =
        makeMarkedElementsFromBoxes<d,T>(boxes, targetLevel, helper);

    marked = helper.markAdmissible(marked, admissibleJump);
    std::vector<index_t> refBoxes = helper.toRefBoxes(marked, false);
    shuffleTuples<d>(refBoxes, rng);
    thb.refineElements(refBoxes);
}

// ============================================================
// Boxes and partitions
// ============================================================

/// Axis-aligned box of one level, in level-local unique-knot indices:
/// the box covers the elements lo[dd] <= i_dd < up[dd] of that level.
template<short_t d>
struct Box
{
    index_t lo[d];
    index_t up[d];
};

/// A partition of the active region: one vector of boxes per level
/// (empty for a level without active cells).
template<short_t d>
using Partition = std::vector<std::vector<Box<d>>>;

template<short_t d>
bool boxLess(const Box<d> & x, const Box<d> & y)
{
    for (short_t dd = 0; dd < d; ++dd)
        if (x.lo[dd] != y.lo[dd]) return x.lo[dd] < y.lo[dd];
    for (short_t dd = 0; dd < d; ++dd)
        if (x.up[dd] != y.up[dd]) return x.up[dd] < y.up[dd];
    return false;
}

char axisLetter(short_t dd) { return "xyzw"[dd]; }

/// Candidate `kdtree`: the leaves of the kd-tree, collected per level in
/// leaf-iterator order (the order in which gsHPatchRule::_build visits them).
/// Complexity: O(number of leaves).
template<short_t d, class T>
Partition<d> collectKdLeaves(const gsHTensorBasis<d,T> & basis)
{
    Partition<d> parts(basis.maxLevel() + 1);
    auto leaf = basis.tree().beginLeafIterator();
    for (; leaf.good(); leaf.next())
    {
        Box<d> b;
        for (short_t dd = 0; dd < d; ++dd)
        {
            b.lo[dd] = leaf.lowerCorner()[dd];
            b.up[dd] = leaf.upperCorner()[dd];
        }
        parts[leaf.level()].push_back(b);
    }
    return parts;
}

template<short_t d>
long long countBoxes(const Partition<d> & parts)
{
    long long n = 0;
    for (size_t l = 0; l < parts.size(); ++l)
        n += static_cast<long long>(parts[l].size());
    return n;
}

// ============================================================
// Node model
// ============================================================

/// Number of nodes of element-wise Gauss quadrature on a box, with basis degree
/// p+1 points per element and direction: prod_d (up-lo)(p_d+1). This is the
/// nGaussTotal of gsHPatchRule::_build, which the fallback fill reproduces.
template<short_t d, class T>
long long gaussCountOfBox(const gsTHBSplineBasis<d,T> & thb,
                          const index_t * lo, const index_t * up)
{
    long long g = 1;
    for (short_t dd = 0; dd < d; ++dd)
        g *= static_cast<long long>(up[dd] - lo[dd]) * (thb.degree(dd) + 1);
    return g;
}

/// Node count of gsHPatchRule (volume rule, fixDir = -1) on one box at level
/// `lev`; `countFallback` is set when _build would fill the box with
/// element-wise Gauss because the patch rule has no fewer nodes.
///
/// Mirror of gsHPatchRule::_build. In direction dd the target knot vector has
/// multiplicity tDeg+1 at both ends of the box and, at an interior unique knot
/// of multiplicity mu_j in the level knot vector, min(tDeg - tReg - 1 + mu_j,
/// tDeg + 1). Its spline space has dimension N0 = (sum of multiplicities) -
/// (tDeg+1). _ensureEven1D leaves an even N0 alone and inserts one knot into an
/// odd one, so N = N0 + (N0 % 2) and the 1D rule has N/2 = ceil(N0/2) nodes.
/// The tensor rule has the product of these; it is replaced by Gauss when it
/// is not smaller than the Gauss count (the decision is per box).
///
/// Facts used by the validation: N0 <= (2p+1) n_d when tDeg = 2p, so
/// ceil(N0/2) <= (p+1) n_d with equality only for n_d = 1. The count fallback
/// can therefore only fire on boxes one element wide in every direction, and
/// there both counts coincide. The Newton iteration of _compute1D can also
/// fail and trigger the fallback; the model assumes it converges.
///
/// Complexity: O(sum_d (up_d - lo_d)) multiplicity lookups.
template<short_t d, class T>
long long nodesOfBox(const gsTHBSplineBasis<d,T> & thb, index_t lev,
                     const index_t * lo, const index_t * up,
                     index_t tDeg, index_t tReg, bool & countFallback)
{
    long long candidate = 1;
    for (short_t dd = 0; dd < d; ++dd)
    {
        const gsKnotVector<T> & kv = thb.tensorLevel(lev).component(dd).knots();
        long long sumMult = 0;
        for (index_t i = lo[dd]; i <= up[dd]; ++i)
        {
            if (i == lo[dd] || i == up[dd])
                sumMult += tDeg + 1;
            else
            {
                const index_t muj = kv.multiplicity(kv(i));
                sumMult += math::min(tDeg - tReg - 1 + muj, tDeg + 1);
            }
        }
        const long long N0 = sumMult - (tDeg + 1);
        candidate *= (N0 + 1) / 2;
    }

    const long long gauss = gaussCountOfBox<d,T>(thb, lo, up);
    countFallback = (candidate >= gauss);
    return countFallback ? gauss : candidate;
}

template<short_t d, class T>
long long nodesOfBoxes(const gsTHBSplineBasis<d,T> & thb, index_t lev,
                       const std::vector<Box<d>> & boxes,
                       index_t tDeg, index_t tReg, long long * nCountFallback = nullptr)
{
    long long n = 0;
    for (size_t k = 0; k < boxes.size(); ++k)
    {
        bool fb = false;
        n += nodesOfBox<d,T>(thb, lev, boxes[k].lo, boxes[k].up, tDeg, tReg, fb);
        if (nCountFallback != nullptr && fb) ++(*nCountFallback);
    }
    return n;
}

template<short_t d, class T>
long long nodesOfPartition(const gsTHBSplineBasis<d,T> & thb, const Partition<d> & parts,
                           index_t tDeg, index_t tReg)
{
    long long n = 0;
    for (size_t l = 0; l < parts.size(); ++l)
        n += nodesOfBoxes<d,T>(thb, static_cast<index_t>(l), parts[l], tDeg, tReg);
    return n;
}

/// Lower bound on the node count of ANY exact partition of the active region,
/// ceil(nElem * a^d / 2^d) with a = tDeg - tReg. It does not depend on the
/// candidate; the three rows of a form carry the same value.
///
/// Proof. In a direction with n elements in the box, every interior
/// multiplicity is min(a-1+mu_j, tDeg+1) >= a (mu_j >= 1, and a <= tDeg+1
/// because tReg >= -1), so N0 >= (tDeg+1) + (n-1) a >= n a and the 1D node
/// count ceil(N0/2) is at least n a / 2. A box with E elements thus has at
/// least E (a/2)^d patch nodes. The Gauss count per element, (p+1)^d, is also
/// at least (a/2)^d because a <= 2p+1, so the fallback does not undercut it.
/// Summing over the boxes of an exact partition gives nElem (a/2)^d, and the
/// ceiling holds since node counts are integers. Exact integer arithmetic.
long long nodesLowerBound(long long nElem, index_t d, index_t tDeg, index_t tReg)
{
    const long long a = tDeg - tReg;
    long long pw = 1, twoPw = 1;
    for (index_t dd = 0; dd < d; ++dd) { pw *= a; twoPw *= 2; }
    return (nElem * pw + twoPw - 1) / twoPw;
}

// ============================================================
// Candidate pairwise
// ============================================================

/// True when x and y have the same extents in every axis except `a`.
template<short_t d>
bool sameCrossSection(const Box<d> & x, const Box<d> & y, short_t a)
{
    for (short_t dd = 0; dd < d; ++dd)
        if (dd != a && (x.lo[dd] != y.lo[dd] || x.up[dd] != y.up[dd]))
            return false;
    return true;
}

/// Orders boxes by their extents in all axes except `a`, then by lo[a]; boxes
/// of equal cross-section are then adjacent and sorted along `a`.
template<short_t d>
struct CrossSectionLess
{
    short_t a;
    bool operator()(const Box<d> & x, const Box<d> & y) const
    {
        for (short_t dd = 0; dd < d; ++dd)
        {
            if (dd == a) continue;
            if (x.lo[dd] != y.lo[dd]) return x.lo[dd] < y.lo[dd];
            if (x.up[dd] != y.up[dd]) return x.up[dd] < y.up[dd];
        }
        return x.lo[a] < y.lo[a];
    }
};

/// Candidate `pairwise` on the boxes of one level: repeatedly merges boxes
/// whose union is a box, i.e. that agree in all axes but one and touch along
/// it, until a full sweep over the axes merges nothing.
///
/// Complexity: one pass over an axis sorts the L boxes, O(L log L), and merges
/// maximal runs in one scan. Every productive pass removes at least one box, so
/// there are at most L passes (O(L^2 log L) worst case; usually a few).
///
/// Node count: a merge keeps the cross-section and gives, along the merge
/// axis, N0_merged = N0_A + N0_B + mult_shared - (tDeg+1) <= N0_A + N0_B, so
/// nodes(merged) <= nodes(A) + nodes(B) and the merged partition never has
/// more nodes than its input.
template<short_t d>
void mergePairwise(std::vector<Box<d>> & boxes)
{
    bool changed = true;
    while (changed)
    {
        changed = false;
        for (short_t a = 0; a < d; ++a)
        {
            CrossSectionLess<d> cmp;
            cmp.a = a;
            std::sort(boxes.begin(), boxes.end(), cmp);

            std::vector<Box<d>> merged;
            merged.reserve(boxes.size());
            for (size_t k = 0; k < boxes.size(); ++k)
            {
                if (!merged.empty() && sameCrossSection<d>(merged.back(), boxes[k], a)
                    && merged.back().up[a] == boxes[k].lo[a])
                {
                    merged.back().up[a] = boxes[k].up[a];
                    changed = true;
                }
                else
                    merged.push_back(boxes[k]);
            }
            boxes.swap(merged);
        }
    }
}

template<short_t d>
Partition<d> pairwisePartition(const Partition<d> & kd)
{
    Partition<d> out = kd;
    for (size_t l = 0; l < out.size(); ++l)
        mergePairwise<d>(out[l]);
    return out;
}

// ============================================================
// Candidate sweep
// ============================================================

/// Canonical slab decomposition of the union of the disjoint boxes B for the
/// axis order sigma, from depth k on.
///
/// At depth k the cut set C holds every face coordinate of B along a =
/// sigma[k]. Each slab [c0,c1] between consecutive cuts is decomposed by the
/// recursive call (axes sigma[k+1..]) on the boxes spanning the slab, clipped
/// to it; a piece is glued to the piece of the previous slab that has the same
/// extents in the remaining axes and touches it. The last axis merges touching
/// intervals directly.
///
/// Why this is canonical. C contains every face coordinate of the kd leaves,
/// including interior faces that depend on the insertion history. At a cut
/// where the cross-section of the region is the same on both sides the
/// recursive results on both sides are identical, so every piece continues
/// across the cut and the cut leaves no trace. Every other cut is a face of
/// the region itself. The output therefore depends on the region and sigma,
/// not on the kd-tree. Merging happens inner axes first (sigma[d-1] up to
/// sigma[0]); a second pass would find nothing, because pieces touching along
/// sigma[k] with equal other extents are already glued inside the slabs.
///
/// Complexity. A depth has at most 2L-1 slabs for L boxes and the pieces before
/// gluing number at most (2L)^(d-1) L = O(L^d), so the result can contain more
/// boxes than L (a plus-shaped region made of 2 kd boxes gives 3 sweep boxes).
/// The naive slab filter costs O(|B|) per slab, and ordered maps add a log
/// factor: O(L^d log L) per axis order, O(d! L^d log L) over all orders. In
/// practice the pieces are far fewer than the bound. Ordered containers only,
/// so the result does not depend on hash order.
template<short_t d>
std::vector<Box<d>> sweepRec(const std::vector<Box<d>> & B,
                             const std::array<short_t,d> & sigma, short_t k)
{
    const short_t a = sigma[k];
    std::vector<Box<d>> result;

    if (k == d - 1)
    {
        // All boxes share their extents in the other axes (fixed by the outer
        // slabs): merge touching intervals along a.
        std::vector<Box<d>> sorted(B);
        std::sort(sorted.begin(), sorted.end(),
                  [a](const Box<d> & x, const Box<d> & y) { return x.lo[a] < y.lo[a]; });
        for (size_t i = 0; i < sorted.size(); ++i)
        {
            if (!result.empty() && result.back().up[a] == sorted[i].lo[a])
                result.back().up[a] = sorted[i].up[a];
            else
                result.push_back(sorted[i]);
        }
        return result;
    }

    std::vector<index_t> cuts;
    cuts.reserve(2 * B.size());
    for (size_t i = 0; i < B.size(); ++i)
    {
        cuts.push_back(B[i].lo[a]);
        cuts.push_back(B[i].up[a]);
    }
    std::sort(cuts.begin(), cuts.end());
    cuts.erase(std::unique(cuts.begin(), cuts.end()), cuts.end());

    typedef std::vector<index_t> Key;
    std::map<Key, size_t> open;
    for (size_t c = 0; c + 1 < cuts.size(); ++c)
    {
        const index_t c0 = cuts[c], c1 = cuts[c + 1];
        std::vector<Box<d>> S;
        for (size_t i = 0; i < B.size(); ++i)
            if (B[i].lo[a] <= c0 && B[i].up[a] >= c1)
            {
                Box<d> clipped = B[i];
                clipped.lo[a] = c0;
                clipped.up[a] = c1;
                S.push_back(clipped);
            }
        if (S.empty())
        {
            open.clear();
            continue;
        }

        const std::vector<Box<d>> P = sweepRec<d>(S, sigma, k + 1);
        std::map<Key, size_t> newOpen;
        for (size_t i = 0; i < P.size(); ++i)
        {
            Key key;
            for (short_t j = k + 1; j < d; ++j)
            {
                key.push_back(P[i].lo[sigma[j]]);
                key.push_back(P[i].up[sigma[j]]);
            }
            const std::map<Key, size_t>::const_iterator it = open.find(key);
            if (it != open.end() && result[it->second].up[a] == c0)
            {
                result[it->second].up[a] = c1;
                newOpen[key] = it->second;
            }
            else
            {
                result.push_back(P[i]);
                newOpen[key] = result.size() - 1;
            }
        }
        open.swap(newOpen);
    }
    return result;
}

/// Candidate `sweep`: per level, runs sweepRec for every axis order (starting
/// from 0,1,..,d-1 in lexicographic order) and keeps the order with the fewest
/// nodes under the node model, ties to the lexicographically smallest order.
/// The levels are disjoint regions partitioned independently, so the per-level
/// choice is never worse than a single order for all levels. `axisOrder` lists
/// the chosen orders from level 0 upward joined by '|' ('-' for an empty level).
/// Not guaranteed to have no more nodes than kdtree: it refines first and then
/// merges greedily.
/// Complexity: d! orders times the cost of sweepRec, plus the node model.
template<short_t d, class T>
Partition<d> sweepPartition(const gsTHBSplineBasis<d,T> & thb, const Partition<d> & kd,
                            index_t tDeg, index_t tReg, std::string & axisOrder)
{
    Partition<d> out(kd.size());
    axisOrder.clear();
    for (size_t l = 0; l < kd.size(); ++l)
    {
        if (l > 0) axisOrder += '|';
        if (kd[l].empty())
        {
            axisOrder += '-';
            continue;
        }

        std::array<short_t,d> sigma, best;
        for (short_t dd = 0; dd < d; ++dd) sigma[dd] = dd;
        best = sigma;
        long long bestNodes = -1;
        do
        {
            std::vector<Box<d>> pieces = sweepRec<d>(kd[l], sigma, 0);
            const long long n = nodesOfBoxes<d,T>(thb, static_cast<index_t>(l), pieces, tDeg, tReg);
            if (bestNodes < 0 || n < bestNodes)
            {
                bestNodes = n;
                best = sigma;
                std::sort(pieces.begin(), pieces.end(), boxLess<d>);
                out[l].swap(pieces);
            }
        } while (std::next_permutation(sigma.begin(), sigma.end()));

        for (short_t dd = 0; dd < d; ++dd) axisOrder += axisLetter(best[dd]);
    }
    return out;
}

// ============================================================
// Exactness check
// ============================================================

/// Sorted linear indices of all level-`lev` elements covered by `boxes`.
/// Cell (i_0..i_{d-1}) maps to sum_dd i_dd * prod_{e<dd} n_e with n_e the
/// number of elements of the level in direction e.
template<short_t d, class T>
std::vector<long long> coveredCells(const gsTHBSplineBasis<d,T> & thb, index_t lev,
                                    const std::vector<Box<d>> & boxes)
{
    long long stride[d];
    long long s = 1;
    for (short_t dd = 0; dd < d; ++dd)
    {
        stride[dd] = s;
        s *= thb.tensorLevel(lev).component(dd).knots().numElements();
    }

    std::vector<long long> cells;
    for (size_t k = 0; k < boxes.size(); ++k)
    {
        index_t idx[d];
        for (short_t dd = 0; dd < d; ++dd) idx[dd] = boxes[k].lo[dd];
        bool more = true;
        while (more)
        {
            long long lin = 0;
            for (short_t dd = 0; dd < d; ++dd) lin += idx[dd] * stride[dd];
            cells.push_back(lin);

            short_t dd = 0;
            for (; dd < d; ++dd)
            {
                if (++idx[dd] < boxes[k].up[dd]) break;
                idx[dd] = boxes[k].lo[dd];
            }
            more = (dd < d);
        }
    }
    std::sort(cells.begin(), cells.end());
    return cells;
}

/// Aborts unless `cand` covers exactly the cells of the kd leaves `kdCells`
/// (per level, sorted): a duplicate is an overlap, a missing index a gap and an
/// extra index a cell outside the region.
template<short_t d, class T>
void ensureExact(const gsTHBSplineBasis<d,T> & thb, const Partition<d> & cand,
                 const std::vector<std::vector<long long>> & kdCells,
                 const std::string & name)
{
    GISMO_ENSURE(cand.size() == kdCells.size(),
                 "partition " << name << ": wrong number of levels");
    for (size_t l = 0; l < cand.size(); ++l)
    {
        const std::vector<long long> cells = coveredCells<d,T>(thb, static_cast<index_t>(l), cand[l]);
        const size_t m = std::min(cells.size(), kdCells[l].size());
        for (size_t i = 0; i < m; ++i)
            GISMO_ENSURE(cells[i] == kdCells[l][i],
                         "partition " << name << ", level " << l << ": cell index "
                         << cells[i] << " differs from the region's " << kdCells[l][i]
                         << " at sorted position " << i << " (overlap, gap or foreign cell)");
        GISMO_ENSURE(cells.size() == kdCells[l].size(),
                     "partition " << name << ", level " << l << ": covers " << cells.size()
                     << " cells, the region has " << kdCells[l].size()
                     << " (first differing position " << m << ")");
    }
}

// ============================================================
// Model validation against the real rule (kdtree only)
// ============================================================

struct RuleCheck
{
    long long nNodesRule     = 0;  ///< nodes of the real rule over all elements
    long long nFallbackRule  = 0;  ///< cells with usedFallback (count + Newton)
    long long fModel         = 0;  ///< leaves the model puts in count fallback
    long long nNewton        = 0;  ///< leaves where the real rule fell back unpredicted
    long long nMismatch      = 0;  ///< leaves matching neither case
    bool      matches        = false;
};

/// Builds the real gsHPatchRule and compares it leaf by leaf with the model.
/// A leaf matches when the real node count equals the model's. When the model
/// predicted the patch rule (strictly fewer nodes than Gauss) and the real
/// count is the Gauss count, the Newton iteration failed there: a Newton
/// fallback, which a model cannot predict. Any other leaf is a mismatch.
/// `matches` requires no mismatch and nFallbackRule == fModel + nNewton, so
/// with nNewton == 0 it implies that the modelled total equals nNodesRule.
template<short_t d, class T>
RuleCheck validateModel(const gsTHBSplineBasis<d,T> & thb, const Partition<d> & kd,
                        index_t tDeg, index_t tReg)
{
    RuleCheck rc;
    gsHPatchRule<T> rule(*thb.domain(), tDeg, tReg, false);

    gsMatrix<T> nodes;
    gsVector<T> weights;
    for (auto & elem : thb.domain()->allElements())
    {
        rule.mapTo(elem.lowerCorner(), elem.upperCorner(), nodes, weights);
        rc.nNodesRule += weights.size();
    }
    rc.nFallbackRule = rule.numFallbackCells();

    long long sumLeaves = 0;
    gsVector<T> lower(d), upper(d);
    for (size_t l = 0; l < kd.size(); ++l)
    {
        for (size_t k = 0; k < kd[l].size(); ++k)
        {
            const Box<d> & b = kd[l][k];
            for (short_t dd = 0; dd < d; ++dd)
            {
                const gsKnotVector<T> & kv =
                    thb.tensorLevel(static_cast<index_t>(l)).component(dd).knots();
                lower[dd] = kv(b.lo[dd]);
                upper[dd] = kv(b.up[dd]);
            }
            rule.mapTo(lower, upper, nodes, weights);
            const long long real = weights.size();
            sumLeaves += real;

            bool fb = false;
            const long long model =
                nodesOfBox<d,T>(thb, static_cast<index_t>(l), b.lo, b.up, tDeg, tReg, fb);
            const long long gauss = gaussCountOfBox<d,T>(thb, b.lo, b.up);
            if (fb) ++rc.fModel;
            if (real == model) continue;
            if (!fb && real == gauss) ++rc.nNewton;
            else ++rc.nMismatch;
        }
    }
    GISMO_ENSURE(sumLeaves == rc.nNodesRule,
                 "per-leaf node counts (" << sumLeaves << ") do not sum to the element-loop count ("
                 << rc.nNodesRule << ")");
    rc.matches = (rc.nMismatch == 0 && rc.nFallbackRule == rc.fModel + rc.nNewton);
    return rc;
}

// ============================================================
// Reporting
// ============================================================

struct Row
{
    std::string form, candidate, axisOrder;
    long long   nBoxes = 0, nNodes = 0, lowerBound = 0;
    bool        hasRule = false;
    RuleCheck   rule;
    real_t      tPartition = 0;
};

struct FormSpec
{
    std::string name;
    index_t     tDeg, tReg;
};

std::string sanitize(std::string s)
{
    for (size_t i = 0; i < s.size(); ++i)
        if (s[i] == ',' || s[i] == '\n' || s[i] == '\r') s[i] = ';';
    return s;
}

std::string sci(real_t v)
{
    std::ostringstream os;
    os << std::scientific << std::setprecision(6) << v;
    return os.str();
}

/// Smallest wall time of `reps` runs of f.
template<class F>
real_t minTime(index_t reps, F f)
{
    real_t best = -1;
    for (index_t r = 0; r < reps; ++r)
    {
        gsStopwatch sw;
        f();
        const real_t t = sw.stop();
        if (best < 0 || t < best) best = t;
    }
    return best;
}

// ============================================================
// Templated driver
// ============================================================

template<short_t d, class T>
int run(index_t numLevels, index_t p, index_t r, MeshPattern pattern, BandVariant band,
        const std::string & patternName, bool admissible, index_t admissibleJump,
        index_t shuffleSeed, index_t reps, index_t baseRefs,
        const std::string & csvFile)
{
    // Basis: same coarse tensor basis and refinement sequence as
    // hierarchical_assembly_hpatch, so the meshes are the benchmark's meshes.
    gsKnotVector<T> kv(0, 1.0, 0, p + 1, 1);
    for (index_t i = 0; i < baseRefs; ++i)
        kv.uniformRefine();
    std::vector<gsKnotVector<T>> kvs;
    for (short_t dd = 0; dd < d; ++dd) kvs.push_back(kv);
    gsTensorBSplineBasis<d,T> tbsb(kvs);
    gsTHBSplineBasis<d,T>     thb(tbsb);

    // The generator exists only with --shuffle; it perturbs the order in which
    // refinement boxes enter the kd-tree.
    std::mt19937 rngStorage(shuffleSeed >= 0 ? static_cast<unsigned>(shuffleSeed) : 0u);
    std::mt19937 * rng = shuffleSeed >= 0 ? &rngStorage : nullptr;

    thb.uniformRefine();
    for (index_t lvl = 1; lvl <= numLevels; ++lvl)
        refineLevel<d,T>(thb, pattern, baseRefs, band, lvl, admissible, admissibleJump, rng);

    GISMO_ENSURE(!thb.manualLevels(), "manualLevels bases are not supported.");

    const long long nElem = thb.numElements();
    const Partition<d> kdRef = collectKdLeaves<d,T>(thb);

    // The kd leaves are the reference region: disjoint cells whose number must
    // equal the number of active elements. Their sorted cell indices are what
    // every candidate has to reproduce.
    std::vector<std::vector<long long>> kdCells(kdRef.size());
    long long totalCells = 0;
    for (size_t l = 0; l < kdRef.size(); ++l)
    {
        kdCells[l] = coveredCells<d,T>(thb, static_cast<index_t>(l), kdRef[l]);
        for (size_t i = 1; i < kdCells[l].size(); ++i)
            GISMO_ENSURE(kdCells[l][i] != kdCells[l][i - 1],
                         "kd leaves overlap at level " << l << ", cell index " << kdCells[l][i]);
        totalCells += static_cast<long long>(kdCells[l].size());
    }
    GISMO_ENSURE(totalCells == nElem,
                 "kd leaves cover " << totalCells << " cells, the basis has " << nElem << " elements");

    // Header block.
    gsInfo << "dim=" << d << "  p=" << p << "  r=" << r << "  levels=" << numLevels
           << "  pattern=" << patternName
           << "  admissible=" << (admissible ? 1 : 0) << "  jump=" << admissibleJump
           << "  shuffle=" << shuffleSeed << "  reps=" << reps << "\n";
    gsInfo << "nElem=" << nElem << "  maxLevel=" << thb.maxLevel()
           << "  DOF=" << thb.size() << "\n";
    gsInfo << "kd leaves per level:";
    for (size_t l = 0; l < kdRef.size(); ++l) gsInfo << " " << kdRef[l].size();
    gsInfo << "\n\n";

    // Forms: target degree 2p, target regularity as in the benchmark driver.
    std::vector<FormSpec> forms;
    forms.push_back({"mass", 2 * p, r});
    forms.push_back({"stiffness", 2 * p, r - 1});
    if (r >= 2)
        forms.push_back({"biharmonic", 2 * p, r - 2});

    std::vector<Row> rows;
    for (size_t f = 0; f < forms.size(); ++f)
    {
        const FormSpec & fs = forms[f];
        const long long lb = nodesLowerBound(nElem, d, fs.tDeg, fs.tReg);

        // kdtree: leaf collection only (the model and the real rule are not timed).
        Partition<d> kd;
        Row rowKd;
        rowKd.tPartition = minTime(reps, [&]() { kd = collectKdLeaves<d,T>(thb); });
        ensureExact<d,T>(thb, kd, kdCells, "kdtree");
        rowKd.form = fs.name; rowKd.candidate = "kdtree";
        rowKd.nBoxes = countBoxes<d>(kd);
        rowKd.nNodes = nodesOfPartition<d,T>(thb, kd, fs.tDeg, fs.tReg);
        rowKd.lowerBound = lb;
        rowKd.hasRule = true;
        rowKd.rule = validateModel<d,T>(thb, kd, fs.tDeg, fs.tReg);

        // pairwise: the copy of the kd boxes is part of the timed region.
        Partition<d> pw;
        Row rowPw;
        rowPw.tPartition = minTime(reps, [&]() { pw = pairwisePartition<d>(kd); });
        ensureExact<d,T>(thb, pw, kdCells, "pairwise");
        rowPw.form = fs.name; rowPw.candidate = "pairwise";
        rowPw.nBoxes = countBoxes<d>(pw);
        rowPw.nNodes = nodesOfPartition<d,T>(thb, pw, fs.tDeg, fs.tReg);
        rowPw.lowerBound = lb;
        GISMO_ENSURE(rowPw.nNodes <= rowKd.nNodes,
                     "pairwise has more nodes (" << rowPw.nNodes << ") than kdtree ("
                     << rowKd.nNodes << ") for form " << fs.name);

        // sweep: all axis orders on all levels, including the model evaluations
        // that select the order.
        Partition<d> sw;
        Row rowSw;
        rowSw.tPartition = minTime(reps, [&]() {
            sw = sweepPartition<d,T>(thb, kd, fs.tDeg, fs.tReg, rowSw.axisOrder); });
        ensureExact<d,T>(thb, sw, kdCells, "sweep");
        rowSw.form = fs.name; rowSw.candidate = "sweep";
        rowSw.nBoxes = countBoxes<d>(sw);
        rowSw.nNodes = nodesOfPartition<d,T>(thb, sw, fs.tDeg, fs.tReg);
        rowSw.lowerBound = lb;

        const Row * three[3] = {&rowKd, &rowPw, &rowSw};
        for (index_t c = 0; c < 3; ++c)
        {
            const Row & w = *three[c];
            gsInfo << "row " << w.form << " " << w.candidate
                   << "  axisOrder=" << (w.axisOrder.empty() ? "-" : w.axisOrder)
                   << "  nBoxes=" << w.nBoxes << "  nNodes=" << w.nNodes
                   << "  nNodesLowerBound=" << w.lowerBound
                   << "  tPartition=" << sci(w.tPartition);
            if (w.hasRule)
                gsInfo << "  nNodesRule=" << w.rule.nNodesRule
                       << "  nFallbackRule=" << w.rule.nFallbackRule
                       << "  F_model=" << w.rule.fModel
                       << "  nNewton=" << w.rule.nNewton
                       << "  mismatchLeaves=" << w.rule.nMismatch
                       << "  modelMatchesRule=" << (w.rule.matches ? 1 : 0);
            gsInfo << "\n";
            rows.push_back(w);
        }
        gsInfo << "sweep<=kdtree: " << (rowSw.nNodes <= rowKd.nNodes ? 1 : 0) << " " << fs.name << "\n";
        gsInfo << "pairwise<=kdtree: " << (rowPw.nNodes <= rowKd.nNodes ? 1 : 0) << " " << fs.name << "\n";
    }

    gsInfo << "partition check: OK\n";

    // CSV: one line per (form, candidate); free-text fields have commas and
    // newlines replaced so that every row has exactly 24 fields.
    if (!csvFile.empty())
    {
        char hostBuf[256] = {0};
        gethostname(hostBuf, sizeof(hostBuf) - 1);
        const std::string host = sanitize(hostBuf);
        const std::string cpu  = sanitize(gsSysInfo::getCpuInfo());

        std::ofstream csv(csvFile.c_str());
        GISMO_ENSURE(csv.good(), "cannot open CSV file '" << csvFile << "'");
        csv << "form,dim,p,r,levels,pattern,admissible,jump,nElem,maxLevel,candidate,"
               "axisOrder,nBoxes,nNodes,nNodesLowerBound,nNodesRule,nFallbackRule,"
               "modelMatchesRule,tPartition,reps,threads,host,cpu,baseRefs,band\n";
        for (size_t i = 0; i < rows.size(); ++i)
        {
            const Row & w = rows[i];
            csv << w.form << "," << d << "," << p << "," << r << "," << numLevels << ","
                << patternName << "," << (admissible ? 1 : 0) << "," << admissibleJump << ","
                << nElem << "," << thb.maxLevel() << "," << w.candidate << ","
                << w.axisOrder << "," << w.nBoxes << "," << w.nNodes << "," << w.lowerBound << ",";
            if (w.hasRule)
                csv << w.rule.nNodesRule << "," << w.rule.nFallbackRule << ","
                    << (w.rule.matches ? 1 : 0) << ",";
            else
                csv << ",,,";
            csv << sci(w.tPartition) << "," << reps << "," << omp_get_max_threads() << ","
                << host << "," << cpu << "," << baseRefs << "," << bandName(band) << "\n";
        }
        gsInfo << "wrote " << csvFile << "\n";
    }
    return EXIT_SUCCESS;
}

} // namespace

// ============================================================
// main
// ============================================================

int main(int argc, char * argv[])
{
    index_t p              = 2;
    index_t r              = 1;
    index_t numLevels      = 3;
    index_t dim            = 2;
    index_t admissibleJump = 2;
    index_t shuffleSeed    = -1;
    index_t reps           = 3;
    index_t baseRefs       = 2;
    bool    admissible     = false;
    std::string meshPattern = "lshape";
    std::string band = "scaled";
    std::string csvFile;

    gsCmdLine cmd("Partition study: node count of gsHPatchRule for three box partitions "
                  "(kdtree, pairwise, sweep) of the active region of each level.");
    cmd.addInt   ("P", "deg",    "Polynomial degree",           p);
    cmd.addInt   ("R", "reg",    "Regularity",                  r);
    cmd.addInt   ("l", "levels", "Number of refinement levels", numLevels);
    cmd.addInt   ("d", "dim",    "Problem dimension (1..4)",    dim);
    cmd.addInt   ("",  "admissibleJump",
                  "Jump parameter m for admissible refinement", admissibleJump);
    cmd.addString("",  "meshPattern",
                  "Mesh pattern: lshape, center, diagonal", meshPattern);
    cmd.addSwitch("",  "admissible",
                  "Apply admissible closure to the marked refinement region", admissible);
    cmd.addInt   ("",  "shuffle",
                  "Seed for shuffling the order of refinement boxes (-1: no shuffle)", shuffleSeed);
    cmd.addInt   ("",  "reps",
                  "Repetitions of each timed partition; the minimum is reported (>= 1)", reps);
    cmd.addString("",  "band",
                  "Refinement band variant: scaled (band width grows with --baseRefs) or thin (4 elements wide)", band);
    cmd.addInt   ("",  "baseRefs",
                  "Uniform refinements of the base knot vector "
                  "(2: 4 elements per direction before the level-0 refinement)", baseRefs);
    cmd.addString("",  "csv",
                  "CSV output file (header and one row per form and candidate); empty: none", csvFile);
    try { cmd.getValues(argc, argv); } catch (int rv) { return rv; }

    if (r < 0)
    {
        gsInfo << "Error: regularity must be >= 0.\n";
        return 2;
    }
    if (baseRefs < 1)
    {
        gsInfo << "Error: baseRefs must be >= 1 (pattern boxes must align with parent elements).\n";
        return 2;
    }
    if (r >= p)
    {
        gsInfo << "Error: regularity must be < degree.\n";
        return EXIT_FAILURE;
    }
    if (admissibleJump < 1)
    {
        gsInfo << "Error: admissibleJump must be >= 1.\n";
        return EXIT_FAILURE;
    }
    if (reps < 1)
    {
        gsInfo << "Error: reps must be >= 1.\n";
        return EXIT_FAILURE;
    }

    BandVariant bandVariant = BandScaled;
    if (!parseBandVariant(band, bandVariant))
    {
        gsInfo << "Error: unknown band '" << band
               << "'. Expected one of: scaled, thin.\n";
        return 2;
    }

    MeshPattern pattern = PatternLShape;
    if (!parseMeshPattern(meshPattern, pattern))
    {
        gsInfo << "Error: unknown meshPattern '" << meshPattern
               << "'. Expected one of: lshape, center, diagonal.\n";
        return EXIT_FAILURE;
    }

    if      (dim == 1) return run<1, real_t>(numLevels, p, r, pattern, bandVariant, meshPattern, admissible, admissibleJump, shuffleSeed, reps, baseRefs, csvFile);
    else if (dim == 2) return run<2, real_t>(numLevels, p, r, pattern, bandVariant, meshPattern, admissible, admissibleJump, shuffleSeed, reps, baseRefs, csvFile);
    else if (dim == 3) return run<3, real_t>(numLevels, p, r, pattern, bandVariant, meshPattern, admissible, admissibleJump, shuffleSeed, reps, baseRefs, csvFile);
    else if (dim == 4) return run<4, real_t>(numLevels, p, r, pattern, bandVariant, meshPattern, admissible, admissibleJump, shuffleSeed, reps, baseRefs, csvFile);

    gsInfo << "Error: dimension must be 1..4.\n";
    return EXIT_FAILURE;
}
