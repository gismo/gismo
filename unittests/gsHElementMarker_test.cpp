/** @file gsHElementMarker_test.cpp

    @brief Tests for gsHElementMarker:
    - markCrs dispatches on the "CoarsenRule" option, not "RefineRule";
    - exact DoF counts after refining one cell of a THB basis with the
      refinement extension (isotropic and anisotropic degrees);
    - admissibility, partition of unity and non-negativity under random
      multi-round refinement and refine/coarsen steps;
    - refinement sets that did not come from the last markRef() (hand-edited,
      stale, built with markAdmissible) are refined to admissible meshes;
    - Admissible = false with Extension = true still grows the space;
    - refined regions survive a combined refine/coarsen round.
    - the "CoarsenGroupRule" option (any child / all children / summed) selects
      different sibling groups for coarsening, and is rejected without GARU.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): Testing
**/

#include "gismo_unittest.h"

#include <gsHSplines/gsHElementMarker.h>
#include <gsHSplines/gsTHBSplineBasis.h>

#include <cstdint>
#include <functional>
#include <map>
#include <set>
#include <sstream>

using namespace gismo;

namespace
{

// std::set<element_t,...> with a custom comparator may not define
// operator==, so compare in terms of the comparator instead of relying on
// std::set::operator==.
template <short_t d, class T>
bool sameSet(const std::set<gsHElement<d,T>, typename gsHElement<d,T>::Compare> & a,
             const std::set<gsHElement<d,T>, typename gsHElement<d,T>::Compare> & b)
{
    typename gsHElement<d,T>::Compare cmp;
    if (a.size() != b.size()) return false;
    auto ia = a.begin(); auto ib = b.begin();
    for (; ia != a.end(); ++ia, ++ib)
        if (cmp(*ia,*ib) || cmp(*ib,*ia)) return false;
    return true;
}

// A THB-spline basis over gsKnotVector<real_t>(0,1,3,3) in both directions,
// refined once so that level-1 elements exist.
gsTHBSplineBasis<2,real_t> makeFixtureBasis()
{
    gsKnotVector<real_t> kv(0, 1, 3, 3);
    gsTensorBSplineBasis<2,real_t> tbasis(kv, kv);
    gsTHBSplineBasis<2,real_t> basis(tbasis);
    basis.refineElements({{1, 0, 0, 2, 2}}); // level 1, box [0,2] x [0,2]
    return basis;
}

typedef gsHElementMarker<2,real_t> Marker;

// Level of the leaf element the iterator points at.
index_t leafLevel(const gsDomainIterator<real_t> & it)
{
    return static_cast<const gsHDomainIterator<real_t,2>&>(it).getLevel();
}

// Class-2 T-admissibility: on every leaf, the truncated functions acting on
// the leaf span at most two consecutive levels (gsHElementHelper::isAdmissible).
bool admissibleClass2(const gsTHBSplineBasis<2,real_t> & b, std::string * msg)
{
    const gsHElementHelper<2,real_t> h(b);
    if (h.isAdmissible(2))
        return true;

    const gsHElementHelper<2,real_t>::HElementContainer bad = h.getNonAdmissibleElements(2);
    const gsHElement<2,real_t> & e = *bad.begin();
    std::ostringstream os;
    os << "leaf level " << e.level() << " from (" << e.lowerCorner().transpose() << ") to ("
       << e.upperCorner().transpose() << ") spans more than two active levels";
    if (msg) *msg = os.str();
    gsTestInfo << os.str() << "\n";
    return false;
}

// Partition of unity and non-negativity, sampled at the centre of every leaf
// and on a uniform 6x6 grid including the boundary.
void checkSpace(const gsTHBSplineBasis<2,real_t> & b, real_t & pouErr, real_t & minVal)
{
    std::vector<real_t> px, py;
    for (auto it = b.domain()->beginAll(); it != b.domain()->endAll(); ++it)
    {
        const gsMatrix<real_t> lo = it.lowerCorner(), up = it.upperCorner();
        px.push_back(0.5*(lo(0,0) + up(0,0)));
        py.push_back(0.5*(lo(1,0) + up(1,0)));
    }
    for (int i = 0; i <= 5; ++i)
        for (int j = 0; j <= 5; ++j)
        {
            px.push_back(i/5.0);
            py.push_back(j/5.0);
        }
    gsMatrix<real_t> pts(2, px.size());
    for (size_t j = 0; j != px.size(); ++j)
    {
        pts(0,j) = px[j];
        pts(1,j) = py[j];
    }
    gsMatrix<real_t> v = b.eval(pts);
    pouErr = 0;
    for (index_t j = 0; j < v.cols(); ++j)
        pouErr = (std::max)(pouErr, math::abs(v.col(j).sum() - 1));
    minVal = v.minCoeff();
}

// Portable deterministic generator with values in [0,1).
struct Lcg
{
    uint32_t s;
    explicit Lcg(uint32_t seed) : s(seed*2654435761u + 1u) {}
    real_t operator()()
    {
        s = 1664525u*s + 1013904223u;
        return (s >> 8) * (1.0/16777216.0);
    }
};

gsTHBSplineBasis<2,real_t> makeUniformTHB(int p, int nel)
{
    gsKnotVector<real_t> kv(0, 1, nel-1, p+1);
    gsTensorBSplineBasis<2,real_t> tb(kv, kv);
    return gsTHBSplineBasis<2,real_t>(tb);
}

// err = 1 on the leaves selected by pred(level, ix, iy), where (ix,iy) is the
// integer index of the leaf on its own level (nel spans on level 0), else 0.
// Returns the number of selected leaves through nSel.
// nx and ny are the numbers of spans on level 0 in the two directions.
std::vector<real_t> errorsMarking(const gsTHBSplineBasis<2,real_t> & b, int nx, int ny,
                                  const std::function<bool(index_t,index_t,index_t)> & pred,
                                  index_t & nSel)
{
    std::vector<real_t> err(b.numElements(), 0.0);
    nSel = 0;
    for (auto it = b.domain()->beginAll(); it != b.domain()->endAll(); ++it)
    {
        const index_t lv = leafLevel(*it.get());
        const real_t scale = real_t(index_t(1) << lv);
        const gsMatrix<real_t> lo = it.lowerCorner();
        const index_t ix = static_cast<index_t>(std::floor(lo(0,0)*nx*scale + 0.5));
        const index_t iy = static_cast<index_t>(std::floor(lo(1,0)*ny*scale + 0.5));
        if (pred(lv, ix, iy))
        {
            err[it.id()] = 1.0;
            ++nSel;
        }
    }
    return err;
}

// Square variant: nel spans on level 0 in both directions.
std::vector<real_t> errorsMarking(const gsTHBSplineBasis<2,real_t> & b, int nel,
                                  const std::function<bool(index_t,index_t,index_t)> & pred,
                                  index_t & nSel)
{
    return errorsMarking(b, nel, nel, pred, nSel);
}

gsOptionList markerOptions(int maxLevel)
{
    gsOptionList o = Marker::defaultOptions();
    o.setInt("RefineRule", 1);
    o.setReal("RefineParam", 0.5);
    o.setInt("MaxLevel", maxLevel);
    o.setInt("Jump", 2);
    o.setSwitch("Extension", true);
    o.setSwitch("Admissible", true);
    return o;
}

// One refinement round with a fresh marker on the current basis; returns the
// number of refinement boxes.
index_t refineOnce(gsTHBSplineBasis<2,real_t> & b, const std::vector<real_t> & err,
                   const gsOptionList & opts)
{
    Marker m(b, opts);
    m.setErrors(err);
    const Marker::HElementContainer marked = m.markRef();
    const std::vector<index_t> boxes = m.toRefBoxes(marked);
    b.refineElements(boxes);
    return static_cast<index_t>(boxes.size() / 5);
}

// Random marking of the leaves below maxLevel with probability q; at least
// one eligible leaf is marked. In single-level mode only leaves of one
// randomly drawn level are eligible.
std::vector<real_t> randomMarking(const gsTHBSplineBasis<2,real_t> & b, Lcg & rng,
                                  bool single, int maxLevel, real_t q)
{
    index_t maxPresent = 0;
    for (auto it = b.domain()->beginAll(); it != b.domain()->endAll(); ++it)
        maxPresent = (std::max)(maxPresent, leafLevel(*it.get()));
    const index_t top = (std::min)(maxPresent, index_t(maxLevel-1));
    const index_t selL = static_cast<index_t>(rng() * (top + 1));

    std::vector<real_t> err(b.numElements(), 0.0);
    index_t first = -1, nMarked = 0;
    for (auto it = b.domain()->beginAll(); it != b.domain()->endAll(); ++it)
    {
        const index_t lv = leafLevel(*it.get());
        if (lv >= maxLevel || (single && lv != selL))
            continue;
        if (first < 0) first = it.id();
        if (rng() < q)
        {
            err[it.id()] = 1.0;
            ++nMarked;
        }
    }
    if (nMarked == 0 && first >= 0)
        err[first] = 1.0;
    return err;
}

// After refining by `boxes` and unrefining, every cell of a refinement box
// (target level L) and every level L-1 cell the box snaps to must still lie
// on a level >= L element. Returns the number of violating boxes.
int refinedRegionLost(const gsTHBSplineBasis<2,real_t> & b, const std::vector<index_t> & boxes)
{
    int bad = 0;
    for (size_t k = 0; k + 5 <= boxes.size(); k += 5)
    {
        const index_t L = boxes[k];
        if (L < 1) continue;
        std::vector<real_t> u[2], us[2];
        for (int i = 0; i < 2; ++i)
        {
            u[i]  = b.tensorLevel(L).knots(i).unique();
            us[i] = b.tensorLevel(L-1).knots(i).unique();
        }
        bool lost = false;
        for (index_t ix = boxes[k+1]; ix < boxes[k+3] && !lost; ++ix)
            for (index_t iy = boxes[k+2]; iy < boxes[k+4] && !lost; ++iy)
            {
                gsMatrix<real_t> x(2,1);
                x(0,0) = 0.5*(u[0][ix] + u[0][ix+1]);
                x(1,0) = 0.5*(u[1][iy] + u[1][iy+1]);
                lost = b.getLevelAtPoint(x) < L;
            }
        for (index_t ix = boxes[k+1]/2; ix <= (boxes[k+3]+1)/2 - 1 && !lost; ++ix)
            for (index_t iy = boxes[k+2]/2; iy <= (boxes[k+4]+1)/2 - 1 && !lost; ++iy)
            {
                gsMatrix<real_t> x(2,1);
                x(0,0) = 0.5*(us[0][ix] + us[0][ix+1]);
                x(1,0) = 0.5*(us[1][iy] + us[1][iy+1]);
                lost = b.getLevelAtPoint(x) < L;
            }
        if (lost) ++bad;
    }
    return bad;
}

// 4x4(x4) cells of degree 2 with three level-0 cells refined to level 1:
// group A = cell (0,0[,0]), B = (2,0[,0]) and C = (0,2[,0]).
gsTHBSplineBasis<2,real_t> makeGroupBasis2D()
{
    gsKnotVector<real_t> kv(0, 1, 3, 3);
    gsTensorBSplineBasis<2,real_t> tb(kv, kv);
    gsTHBSplineBasis<2,real_t> b(tb);
    b.refineElements({1,0,0,2,2, 1,4,0,6,2, 1,0,4,2,6});
    return b;
}

gsTHBSplineBasis<3,real_t> makeGroupBasis3D()
{
    gsKnotVector<real_t> kv(0, 1, 3, 3);
    gsTensorBSplineBasis<3,real_t> tb(kv, kv, kv);
    gsTHBSplineBasis<3,real_t> b(tb);
    b.refineElements({1,0,0,0,2,2,2, 1,4,0,0,6,2,2, 1,0,4,0,2,6,2});
    return b;
}

// Element errors with the maximum 1 on the level-0 cell (3,..,3), 0.5 on every
// other level-0 leaf (never a coarsening candidate) and, for the three sibling
// groups of makeGroupBasis{2,3}D:
//   A: d=2 children (bx,by) = (0,0):0.005 (1,0):0.03 (0,1):0.04 (1,1):0.2;
//      d=3 child (0,0,0):0.005, the other seven 0.2
//   B: all children 0.008
//   C: all children 0.003
template <short_t d>
std::vector<real_t> groupErrors(const gsTHBSplineBasis<d,real_t> & b)
{
    const index_t nel = 4;
    std::vector<real_t> err(b.numElements(), 0.0);
    for (auto it = b.domain()->beginAll(); it != b.domain()->endAll(); ++it)
    {
        const index_t lv = static_cast<const gsHDomainIterator<real_t,d> *>(it.get())->getLevel();
        const gsMatrix<real_t> lo = it.lowerCorner();
        index_t idx[d], cell[d], bit[d];
        bool allMax = true, groupA = true, groupB = true, groupC = true;
        for (short_t i = 0; i != d; ++i)
        {
            idx[i]  = static_cast<index_t>(std::floor(lo(i,0)*nel*real_t(index_t(1) << lv) + 0.5));
            cell[i] = (lv == 0 ? idx[i] : idx[i] >> 1);
            bit[i]  = (lv == 0 ? 0 : idx[i] & 1);
            allMax  = allMax && cell[i] == 3;
            groupA  = groupA && cell[i] == 0;
            groupB  = groupB && cell[i] == (i == 0 ? 2 : 0);
            groupC  = groupC && cell[i] == (i == 1 ? 2 : 0);
        }
        if (lv == 0)
            err[it.id()] = (allMax ? 1.0 : 0.5);
        else if (groupA)
        {
            const bool zero = (bit[0] + bit[1] + (d == 3 ? bit[d-1] : 0) == 0);
            real_t e = 0.2;
            if (zero)
                e = 0.005;
            else if (d == 2 && bit[0] == 1 && bit[1] == 0)
                e = 0.03;
            else if (d == 2 && bit[0] == 0 && bit[1] == 1)
                e = 0.04;
            err[it.id()] = e;
        }
        else if (groupB)
            err[it.id()] = 0.008;
        else if (groupC)
            err[it.id()] = 0.003;
    }
    return err;
}

// GARU coarsening options of the group-rule fixtures: threshold 0.01 times the
// maximal error 1, no admissibility filtering, so that only the group statistic
// decides which groups are selected.
gsOptionList groupRuleOptions(const gsOptionList & base, index_t groupRule)
{
    gsOptionList o = base;
    o.setInt("CoarsenRule", 1);
    o.setReal("CoarsenParam", 0.01);
    o.setSwitch("Admissible", false);
    o.setSwitch("Extension", true);
    o.setInt("MaxLevel", 3);
    o.setInt("CoarsenGroupRule", groupRule);
    return o;
}

// Coarsening boxes of a marked set (toCrsBoxes chunks of 2d+1 entries) with
// the number of times each occurs in the raw output.
template <short_t d>
std::map<std::vector<index_t>, int> crsBoxCounts(const gsHElementMarker<d,real_t> & m,
                                                 const typename gsHElementMarker<d,real_t>::HElementContainer & marked)
{
    const std::vector<index_t> raw = m.toCrsBoxes(marked);
    std::map<std::vector<index_t>, int> counts;
    for (size_t k = 0; k + (2*d+1) <= raw.size(); k += 2*d+1)
        ++counts[std::vector<index_t>(raw.begin() + k, raw.begin() + k + 2*d+1)];
    return counts;
}

template <short_t d>
std::set<std::vector<index_t> > keysOf(const std::map<std::vector<index_t>, int> & m)
{
    std::set<std::vector<index_t> > s;
    for (std::map<std::vector<index_t>, int>::const_iterator it = m.begin(); it != m.end(); ++it)
        s.insert(it->first);
    return s;
}

// Coarsening box {0, lo..., hi...} (level-0 index units) of the level-1
// children of the level-0 cell c.
template <short_t d>
std::vector<index_t> groupBox(const index_t (&c)[d])
{
    std::vector<index_t> box(1, 0);
    for (short_t i = 0; i != d; ++i) box.push_back(c[i]);
    for (short_t i = 0; i != d; ++i) box.push_back(c[i] + 1);
    return box;
}

// Expected box -> multiplicity map for the groups selected by `rule` (see the
// table at the group-rule tests): rule 0 keeps one child of A and every child
// of B and C, rule 1 keeps B and C, rule 2 keeps C, all with 2^d children.
template <short_t d>
std::map<std::vector<index_t>, int> expectedGroupBoxes(index_t rule)
{
    index_t a[d], b[d], c[d];
    for (short_t i = 0; i != d; ++i) { a[i] = 0; b[i] = (i == 0 ? 2 : 0); c[i] = (i == 1 ? 2 : 0); }
    const int nChild = 1 << d;
    std::map<std::vector<index_t>, int> m;
    if (rule == 0) m[groupBox<d>(a)] = 1;
    if (rule <= 1) m[groupBox<d>(b)] = nChild;
    m[groupBox<d>(c)] = nChild;
    return m;
}

// Constructs a marker on `b`, sets `err` and runs markCrs().
template <short_t d>
void runCrs(const gsTHBSplineBasis<d,real_t> & b, const gsOptionList & opts, const std::vector<real_t> & err)
{
    gsHElementMarker<d,real_t> m(b, opts);
    m.setErrors(err);
    m.markCrs();
}

// 4x4 cells of degree 2 with the level-0 cell (0,0) refined to level 1 and its
// level-1 child (0,0) refined further to level 2: the level-1 sibling group of
// cell (0,0) has only three of its four children as leaves.
gsTHBSplineBasis<2,real_t> makeIncompleteGroupBasis2D()
{
    gsKnotVector<real_t> kv(0, 1, 3, 3);
    gsTensorBSplineBasis<2,real_t> tb(kv, kv);
    gsTHBSplineBasis<2,real_t> b(tb);
    b.refineElements({1,0,0,2,2, 2,0,0,2,2});
    return b;
}

// Errors for makeIncompleteGroupBasis2D: maximum 1 on the level-0 cell (3,3),
// 0.5 on the other level-0 and on the level-2 leaves, 0.003 on the three
// level-1 leaves (all far below CoarsenParam*max = 0.01).
std::vector<real_t> incompleteGroupErrors(const gsTHBSplineBasis<2,real_t> & b)
{
    std::vector<real_t> err(b.numElements(), 0.5);
    for (auto it = b.domain()->beginAll(); it != b.domain()->endAll(); ++it)
    {
        const index_t lv = static_cast<const gsHDomainIterator<real_t,2> *>(it.get())->getLevel();
        const gsMatrix<real_t> lo = it.lowerCorner();
        if (lv == 0 && std::floor(lo(0,0)*4 + 0.5) == 3 && std::floor(lo(1,0)*4 + 0.5) == 3)
            err[it.id()] = 1.0;
        else if (lv == 1)
            err[it.id()] = 0.003;
    }
    return err;
}

} // anonymous namespace

SUITE(gsHElementMarker_test)
{

    TEST(MarkCrs_DispatchesOnCoarsenRule_NotRefineRule)
    {
        gsTHBSplineBasis<2,real_t> basis = makeFixtureBasis();

        // Strictly increasing, non-uniform error spread so the threshold,
        // percentage and fraction rules genuinely disagree.
        const size_t nel = basis.numElements();
        std::vector<real_t> err(nel);
        for (size_t i = 0; i != nel; i++)
            err[i] = static_cast<real_t>((i+1)*(i+1));

        // With this fixture (19 elements: 15 at level 0, the 4 refined
        // sub-elements at level 1 carrying the highest errors, since the
        // level-1 ids are the last -- and thus highest-error -- in ascending
        // order), the coarsening candidates are ONLY the 4 level-1 elements
        // (level==0 is always skipped). GARU/PUCA scan ascending order and
        // only reach them once the threshold/percentage covers most of the
        // 15 level-0 elements first; 0.85 was picked (probed empirically) to
        // put GARU and PUCA at genuinely different, nonzero counts.
        const real_t coarsenParam = 0.85;

        gsOptionList optsA = gsHElementMarker<2,real_t>::defaultOptions();
        optsA.setSwitch("Admissible", false); // isolate rule dispatch from _markCrs_admissible
        optsA.setInt("RefineRule", 1);
        optsA.setInt("CoarsenRule", 2);
        optsA.setReal("CoarsenParam", coarsenParam);
        gsHElementMarker<2,real_t> markerA(basis, optsA);
        markerA.setErrors(err);
        auto A = markerA.markCrs(); // empty refined => sibling filters skipped

        gsOptionList optsB = optsA;
        optsB.setInt("RefineRule", 2); // differs from A in RefineRule only
        gsHElementMarker<2,real_t> markerB(basis, optsB);
        markerB.setErrors(err);
        auto B = markerB.markCrs();

        gsOptionList optsC = optsA;
        optsC.setInt("CoarsenRule", 1); // differs from A in CoarsenRule only
        gsHElementMarker<2,real_t> markerC(basis, optsC);
        markerC.setErrors(err);
        auto C = markerC.markCrs();

        gsTestInfo << "sizes: A=" << A.size() << " B=" << B.size() << " C=" << C.size() << "\n";

        // Guard: the fixture must actually distinguish rule 1 from rule 2,
        // otherwise the A==B / A!=C checks below would be vacuous.
        const bool bEqualsC = sameSet(B, C);
        CHECK(bEqualsC == false);

        // markCrs follows CoarsenRule (A and B share CoarsenRule=2, differ
        // only in RefineRule, which markCrs must ignore).
        const bool aEqualsB = sameSet(A, B);
        CHECK(aEqualsB);

        // markCrs does not follow RefineRule: A and C share RefineRule=1 and
        // differ in CoarsenRule (2 vs 1), so the coarsening sets differ.
        const bool aEqualsC = sameSet(A, C);
        CHECK(!aEqualsC);

        // Additional configuration: CoarsenRule = 3 (BULK) should likewise
        // differ from CoarsenRule = 1 (C), if the fixture separates them.
        gsOptionList optsD = optsA;
        optsD.setInt("CoarsenRule", 3);
        gsHElementMarker<2,real_t> markerD(basis, optsD);
        markerD.setErrors(err);
        auto D = markerD.markCrs();
        gsTestInfo << "sizes: D(BULK)=" << D.size() << "\n";
        const bool dEqualsC = sameSet(D, C);
        if (!dEqualsC)
            CHECK(!dEqualsC);
    }

    // Refining one isolated level-0 element of an 8x8 uniform THB basis adds
    // THB functions at every degree. The refined region is the marked cell
    // plus a one-cell ring for p = 2..4 (floor(p/2) spans, boxes snapping to
    // whole cells), clamped at the boundary. For p = 2 (100 functions) the
    // counts follow from [2,5]^2 interior, [0,2]x[2,5] edge, [0,2]^2 corner in
    // cell units: level-1 functions have support [max(0,i-2),min(16,i+1)] in
    // fine spans, level-0 ones [max(0,j-2),j+1] in cells. Interior: 4x4 fine
    // functions added, 1x1 coarse removed (+16-1). Edge: +16, coarse 2x1
    // removed (-2). Corner: +16, coarse 2x2 removed (-4). The p = 3 and p = 4
    // counts were obtained by running the refinement. The extension is
    // symmetric, so the upper edge and corner give the lower edge and corner
    // counts and exercise the upper clamp.
    TEST(Extension_IsolatedElementChangesSpace)
    {
        const int nel = 8;
        const index_t cell[5][2] = {{3,3},{0,3},{0,0},{7,3},{7,7}};
        const char * name[5] = {"interior", "edge", "corner", "upper edge", "upper corner"};
        const index_t before[3] = {100, 121, 144};
        const index_t expected[3][5] = {{115, 114, 112, 114, 112},
                                        {130, 133, 133, 133, 133},
                                        {148, 152, 156, 152, 156}};

        for (int p = 1; p <= 6; ++p)
            for (int c = 0; c < 5; ++c)
            {
                gsTHBSplineBasis<2,real_t> b = makeUniformTHB(p, nel);
                const index_t sizeBefore = b.size(), elBefore = b.numElements();
                index_t nSel = 0;
                const index_t cx = cell[c][0], cy = cell[c][1];
                std::vector<real_t> err = errorsMarking(b, nel,
                    [=](index_t l, index_t ix, index_t iy) { return l == 0 && ix == cx && iy == cy; }, nSel);
                CHECK_EQUAL(1, nSel);
                if (nSel != 1) continue;

                refineOnce(b, err, markerOptions(3));

                gsTestInfo << "p=" << p << " " << name[c] << ": " << sizeBefore << " -> " << b.size() << "\n";
                CHECK(b.numElements() > elBefore);
                CHECK(b.size() > sizeBefore);
                if (p >= 2 && p <= 4)
                {
                    CHECK_EQUAL(before[p-2], sizeBefore);
                    CHECK_EQUAL(expected[p-2][c], b.size());
                }
            }
    }

    // Degrees (3,1) on 8x6 cells: floor(p/2) is 1 in x and 0 in y, so the
    // extension adds a one-cell ring in x only; a mix-up of the per-direction
    // degree changes the counts. Counts for the interior cell (3,2) and the
    // upper corner (7,5) were obtained by running the refinement.
    TEST(Extension_AnisotropicSingleCellCounts)
    {
        const int nx = 8, ny = 6;
        const index_t cell[2][2] = {{3,2},{7,5}};
        const index_t expected[2] = {80, 83};
        for (int c = 0; c < 2; ++c)
        {
            gsKnotVector<real_t> kx(0, 1, nx-1, 4), ky(0, 1, ny-1, 2);
            gsTensorBSplineBasis<2,real_t> tb(kx, ky);
            gsTHBSplineBasis<2,real_t> b(tb);
            const index_t sizeBefore = b.size();
            const index_t cx = cell[c][0], cy = cell[c][1];
            index_t nSel = 0;
            std::vector<real_t> err = errorsMarking(b, nx, ny,
                [=](index_t l, index_t ix, index_t iy) { return l == 0 && ix == cx && iy == cy; }, nSel);
            CHECK_EQUAL(1, nSel);

            refineOnce(b, err, markerOptions(3));

            gsTestInfo << "cell (" << cx << "," << cy << "): " << sizeBefore << " -> " << b.size() << "\n";
            CHECK_EQUAL(77, sizeBefore);
            CHECK(b.size() > sizeBefore);
            CHECK_EQUAL(expected[c], b.size());
            std::string msg;
            CHECK(admissibleClass2(b, &msg));
        }
    }

    // Random multi-round refinement (single-level and mixed-level marking)
    // keeps partition of unity, non-negativity and class-2 admissibility.
    // Seeds 5, 6 and 10 at p = 4 are cases where an extension that ignores the
    // closure produces non-admissible meshes.
    TEST(Extension_RandomMarkingKeepsSpaceValid)
    {
        const int nel = 8, maxLevel = 3, rounds = 5;
        const int seedList[] = {0, 1, 5, 6, 10};
        const real_t q = 0.1;
        for (int p = 1; p <= 4; ++p)
            for (int mode = 0; mode < 2; ++mode)
                for (int seed : seedList)
                {
                    Lcg rng(seed + 100*p + 1000*mode);
                    gsTHBSplineBasis<2,real_t> b = makeUniformTHB(p, nel);
                    for (int r = 0; r < rounds; ++r)
                    {
                        refineOnce(b, randomMarking(b, rng, mode == 0, maxLevel, q), markerOptions(maxLevel));

                        std::string msg;
                        real_t pouErr, minVal;
                        checkSpace(b, pouErr, minVal);
                        const bool adm = admissibleClass2(b, &msg);
                        const bool ok = adm && pouErr < 1e-12 && minVal > -1e-14;
                        CHECK(ok);
                        if (!ok)
                        {
                            gsTestInfo << "FAIL p=" << p << (mode == 0 ? " single" : " mixed")
                                       << " seed=" << seed << " round=" << r << " pouErr=" << pouErr
                                       << " minVal=" << minVal << " " << msg << "\n";
                            break;
                        }
                    }
                }
    }

    // Degree 4 on a 4x4 mesh: a three-round marking sequence that reaches
    // level 2 stays admissible.
    TEST(Extension_MinimalDegree4CaseIsAdmissible)
    {
        const int nel = 4, p = 4;
        gsKnotVector<real_t> kv(0, 1, nel-1, p+1);
        gsTensorBSplineBasis<2,real_t> tb(kv, kv);
        gsTHBSplineBasis<2,real_t> b(tb);
        const gsOptionList opts = markerOptions(2);

        typedef std::pair<index_t,index_t> Idx;
        const std::vector<std::vector<std::pair<index_t,Idx> > > rounds = {
            { {0, Idx(2,2)} },
            { {0, Idx(3,0)} },
            { {1, Idx(6,0)}, {1, Idx(7,1)} } };

        for (size_t r = 0; r != rounds.size(); ++r)
        {
            index_t nSel = 0;
            const std::vector<std::pair<index_t,Idx> > & marks = rounds[r];
            std::vector<real_t> err = errorsMarking(b, nel,
                [&marks](index_t l, index_t ix, index_t iy)
                {
                    for (size_t k = 0; k != marks.size(); ++k)
                        if (marks[k].first == l && marks[k].second.first == ix && marks[k].second.second == iy)
                            return true;
                    return false;
                }, nSel);
            CHECK_EQUAL(static_cast<index_t>(marks.size()), nSel);
            refineOnce(b, err, opts);
        }

        bool hasLevel2 = false;
        for (auto it = b.domain()->beginAll(); it != b.domain()->endAll(); ++it)
            hasLevel2 = hasLevel2 || leafLevel(*it.get()) == 2;
        CHECK(hasLevel2);

        std::string msg;
        real_t pouErr, minVal;
        checkSpace(b, pouErr, minVal);
        CHECK(admissibleClass2(b, &msg));
        CHECK(pouErr < 1e-12);
        CHECK(minVal > -1e-14);
    }

    // Combined refine/coarsen rounds on one pre-update mesh: the refined
    // regions survive the unrefinement (checked by sampling the level at
    // points), and the mesh stays admissible with partition of unity.
    TEST(Extension_RefineCoarsenRoundTrip)
    {
        const int nel = 8, maxLevel = 3, seeds = 2;
        int coarsenedRounds = 0;
        for (int p = 2; p <= 4; ++p)
            for (int seed = 0; seed < seeds; ++seed)
            {
                Lcg rng(seed + 100*p + 7);
                gsTHBSplineBasis<2,real_t> b = makeUniformTHB(p, nel);
                for (int r = 0; r < 3; ++r)
                    refineOnce(b, randomMarking(b, rng, false, maxLevel, 0.1), markerOptions(maxLevel));

                for (int r = 0; r < 2; ++r)
                {
                    // 1 on refine marks, 0 on coarsening candidates, 0.3 neutral.
                    std::vector<real_t> err(b.numElements(), 0.3);
                    index_t nRef = 0;
                    for (auto it = b.domain()->beginAll(); it != b.domain()->endAll(); ++it)
                    {
                        const index_t lv = leafLevel(*it.get());
                        if (lv < maxLevel && rng() < 0.05)
                        {
                            err[it.id()] = 1.0;
                            ++nRef;
                        }
                    }
                    if (nRef == 0)
                        err[b.domain()->beginAll().id()] = 1.0;
                    for (auto it = b.domain()->beginAll(); it != b.domain()->endAll(); ++it)
                        if (leafLevel(*it.get()) >= 1 && err[it.id()] < 0.5 && rng() < 0.5)
                            err[it.id()] = 0.0;

                    Marker m(b, markerOptions(maxLevel));
                    m.setErrors(err);
                    const Marker::HElementContainer markedRef = m.markRef();
                    const std::vector<index_t> boxes = m.toRefBoxes(markedRef);
                    const Marker::HElementContainer markedCrs = m.markCrs(markedRef);
                    const std::vector<index_t> crsBoxes = m.toCrsBoxes(markedCrs);

                    if (!boxes.empty()) b.refineElements(boxes);
                    const size_t nRefined = b.numElements();
                    if (!crsBoxes.empty()) b.unrefineElements(crsBoxes);
                    if (!markedCrs.empty() && b.numElements() < nRefined)
                        ++coarsenedRounds;

                    CHECK_EQUAL(0, refinedRegionLost(b, boxes));

                    std::string msg;
                    real_t pouErr, minVal;
                    checkSpace(b, pouErr, minVal);
                    const bool ok = admissibleClass2(b, &msg) && pouErr < 1e-12 && minVal > -1e-14;
                    CHECK(ok);
                    if (!ok)
                    {
                        gsTestInfo << "FAIL p=" << p << " seed=" << seed << " round=" << r
                                   << " pouErr=" << pouErr << " minVal=" << minVal << " " << msg << "\n";
                        break;
                    }
                }
            }
        gsTestInfo << "rounds that coarsened: " << coarsenedRounds << "\n";
        CHECK(coarsenedRounds > 0);
    }

    // Refinement sets that the last markRef() did not return, on multi-level
    // meshes where the admissible closure adds elements: a hand-edited set
    // (markRef result plus one leaf; the raw marked leaves), a stale set
    // (result of an earlier markRef() after setErrors() changed), and a set
    // built with gsHElementHelper::markAdmissible. toRefBoxes must recompute
    // the closure, so the refined mesh is admissible, the space grows and the
    // partition of unity holds.
    TEST(Extension_ForeignSetsStayAdmissible)
    {
        const int nel = 8, maxLevel = 3;
        const int seedList[] = {0, 1, 5, 6, 10};
        const real_t q = 0.1;
        int rawDiffers = 0, staleDiffers = 0, callerDiffers = 0, plainViolations = 0;

        for (int p = 2; p <= 4; ++p)
            for (int seed : seedList)
                for (int rounds = 2; rounds <= 3; ++rounds)
                {
                    Lcg rng(seed + 100*p + 31 + 1000*rounds);
                    gsTHBSplineBasis<2,real_t> b = makeUniformTHB(p, nel);
                    for (int r = 0; r < rounds; ++r)
                        refineOnce(b, randomMarking(b, rng, false, maxLevel, q), markerOptions(maxLevel));

                    const gsHElementHelper<2,real_t> h(b);
                    const Marker::HElementContainer all = h.toElements();
                    const std::vector<real_t> err1 = randomMarking(b, rng, false, maxLevel, q);
                    const std::vector<real_t> err2 = randomMarking(b, rng, false, maxLevel, q);

                    Marker::HElementContainer raw;
                    for (auto it = b.domain()->beginAll(); it != b.domain()->endAll(); ++it)
                        if (err1[it.id()] == 1.0)
                            raw.insert(h.toElement(it.lowerCorner(), it.upperCorner(), leafLevel(*it.get())));

                    Marker m(b, markerOptions(maxLevel));
                    m.setErrors(err1);
                    const Marker::HElementContainer R = m.markRef();

                    Marker::HElementContainer hand = R;
                    std::vector<gsHElement<2,real_t> > cand;
                    for (const gsHElement<2,real_t> & e : all)
                        if (e.level() < maxLevel && R.count(e) == 0)
                            cand.push_back(e);
                    if (!cand.empty())
                        hand.insert(cand[static_cast<size_t>(rng() * cand.size())]);

                    Marker ms(b, markerOptions(maxLevel));
                    ms.setErrors(err1);
                    const Marker::HElementContainer R1 = ms.markRef();
                    ms.setErrors(err2);
                    const Marker::HElementContainer R2 = ms.markRef();

                    const Marker::HElementContainer caller = h.markAdmissible(raw, 2);
                    Marker mc(b, markerOptions(maxLevel));

                    if (!sameSet(raw, R)) ++rawDiffers;
                    if (!sameSet(R1, R2)) ++staleDiffers;
                    if (!sameSet(caller, raw)) ++callerDiffers;

                    struct Scenario { const char * name; const Marker::HElementContainer * set; Marker * marker; };
                    const Scenario scenarios[] = { {"hand-edited", &hand, &m}, {"raw marked", &raw, &m},
                                                   {"stale", &R1, &ms}, {"caller-built", &caller, &mc} };
                    for (const Scenario & sc : scenarios)
                    {
                        gsTHBSplineBasis<2,real_t> bm = b;
                        bm.refineElements(sc.marker->toRefBoxes(*sc.set));

                        std::string msg;
                        real_t pouErr, minVal;
                        checkSpace(bm, pouErr, minVal);
                        const bool ok = admissibleClass2(bm, &msg) && bm.size() > b.size()
                                        && pouErr < 1e-12 && minVal > -1e-14;
                        CHECK(ok);
                        if (!ok)
                            gsTestInfo << "FAIL " << sc.name << " p=" << p << " seed=" << seed
                                       << " rounds=" << rounds << " pouErr=" << pouErr
                                       << " minVal=" << minVal << " " << msg << "\n";

                        // Plain extension of every element, without the closure.
                        gsTHBSplineBasis<2,real_t> bp = b;
                        bp.refineElements(h.toRefBoxes(*sc.set, true));
                        if (!gsHElementHelper<2,real_t>(bp).isAdmissible(2))
                            ++plainViolations;
                    }
                }

        gsTestInfo << "raw!=R: " << rawDiffers << ", stale R1!=R2: " << staleDiffers
                   << ", caller!=seeds: " << callerDiffers << ", plain-extension violations: "
                   << plainViolations << "\n";
        CHECK(rawDiffers > 0);
        CHECK(staleDiffers > 0);
        CHECK(callerDiffers > 0);
        CHECK(plainViolations > 0);
    }

    // Admissible = false with Extension = true makes no admissibility claim,
    // but every refinement round must still change the (nested) THB space and
    // keep partition of unity and non-negativity.
    TEST(Extension_NonAdmissibleModeGrowsSpace)
    {
        const int nel = 8, maxLevel = 3, rounds = 4;
        const int seedList[] = {0, 1, 5, 6, 10};
        const real_t q = 0.1;
        gsOptionList opts = markerOptions(maxLevel);
        opts.setSwitch("Admissible", false);
        for (int p = 1; p <= 4; ++p)
            for (int mode = 0; mode < 2; ++mode)
                for (int seed : seedList)
                {
                    Lcg rng(seed + 100*p + 1000*mode + 7);
                    gsTHBSplineBasis<2,real_t> b = makeUniformTHB(p, nel);
                    for (int r = 0; r < rounds; ++r)
                    {
                        const index_t sizeBefore = b.size();
                        refineOnce(b, randomMarking(b, rng, mode == 0, maxLevel, q), opts);

                        real_t pouErr, minVal;
                        checkSpace(b, pouErr, minVal);
                        const bool ok = b.size() > sizeBefore && pouErr < 1e-12 && minVal > -1e-14;
                        CHECK(ok);
                        if (!ok)
                        {
                            gsTestInfo << "FAIL p=" << p << (mode == 0 ? " single" : " mixed")
                                       << " seed=" << seed << " round=" << r << " size "
                                       << sizeBefore << " -> " << b.size() << " pouErr=" << pouErr
                                       << " minVal=" << minVal << "\n";
                            break;
                        }
                    }
                }
    }

    // Group rules on 4x4 cells (p=2) with sibling groups A, B, C refined to
    // level 1 and CoarsenParam*max = 0.01 (see groupErrors for the errors).
    //   group | min    max    sqrt(sum err^2) | rule 0  rule 1  rule 2
    //   A     | 0.005  0.2    0.206 (d=3 0.53)| yes     no      no
    //   B     | 0.008  0.008  0.016 (d=3 0.023)| yes    yes     no
    //   C     | 0.003  0.003  0.006 (d=3 0.0085)| yes   yes     yes
    // Rule 0 returns the sub-threshold children themselves (A: 1, B and C:
    // 2^d each), rules 1 and 2 all 2^d children of every passing group.
    // The call is made with empty and non-empty `refined`.
    TEST(CoarsenGroupRule_ThreeRules2D)
    {
        gsTHBSplineBasis<2,real_t> basis = makeGroupBasis2D();
        CHECK_EQUAL(25u, basis.numElements());
        const std::vector<real_t> err = groupErrors<2>(basis);

        for (index_t rule = 0; rule != 3; ++rule)
        {
            gsOptionList opts = groupRuleOptions(Marker::defaultOptions(), rule);
            Marker m(basis, opts);
            m.setErrors(err);

            const Marker::HElementContainer res = m.markCrs();
            const std::map<std::vector<index_t>, int> got = crsBoxCounts<2>(m, res);
            const std::map<std::vector<index_t>, int> want = expectedGroupBoxes<2>(rule);
            CHECK(keysOf<2>(got) == keysOf<2>(want));
            CHECK(got == want);
            CHECK_EQUAL(rule == 0 ? 9u : (rule == 1 ? 8u : 4u), res.size());

            // refinement of the level-0 cell (3,3) only (error 1 > 0.9)
            gsOptionList refOpts = opts;
            refOpts.setInt("RefineRule", 1);
            refOpts.setReal("RefineParam", 0.9);
            Marker mr(basis, refOpts);
            mr.setErrors(err);
            const Marker::HElementContainer R = mr.markRef();
            CHECK_EQUAL(1u, R.size());
            const Marker::HElementContainer resR = mr.markCrs(R);
            const std::map<std::vector<index_t>, int> gotR = crsBoxCounts<2>(mr, resR);
            CHECK(keysOf<2>(gotR) == keysOf<2>(want));
            CHECK(gotR == want);
            CHECK_EQUAL(res.size(), resR.size());

            if (rule == 0)
            {
                // A contributes only its 0.005 child, the one at the origin
                index_t nA = 0;
                for (Marker::HElementContainer::const_iterator it = res.begin(); it != res.end(); ++it)
                    if (it->lowerCorner()(0) < 2 && it->lowerCorner()(1) < 2)
                    {
                        ++nA;
                        CHECK_EQUAL(1u, it->level());
                        CHECK_EQUAL(0, it->lowerCorner()(0));
                        CHECK_EQUAL(0, it->lowerCorner()(1));
                    }
                CHECK_EQUAL(1, nA);
            }
        }
    }

    // Same groups in 3D (8 children each); group A has one child of 0.005 and
    // seven of 0.2.
    TEST(CoarsenGroupRule_ThreeRules3D)
    {
        typedef gsHElementMarker<3,real_t> Marker3;
        gsTHBSplineBasis<3,real_t> basis = makeGroupBasis3D();
        CHECK_EQUAL(85u, basis.numElements());
        const std::vector<real_t> err = groupErrors<3>(basis);

        for (index_t rule = 0; rule != 3; ++rule)
        {
            gsOptionList opts = groupRuleOptions(Marker3::defaultOptions(), rule);
            Marker3 m(basis, opts);
            m.setErrors(err);

            const Marker3::HElementContainer res = m.markCrs();
            const std::map<std::vector<index_t>, int> got = crsBoxCounts<3>(m, res);
            const std::map<std::vector<index_t>, int> want = expectedGroupBoxes<3>(rule);
            CHECK(keysOf<3>(got) == keysOf<3>(want));
            CHECK(got == want);
            CHECK_EQUAL(rule == 0 ? 17u : (rule == 1 ? 16u : 8u), res.size());

            if (rule == 0)
            {
                index_t nA = 0;
                for (Marker3::HElementContainer::const_iterator it = res.begin(); it != res.end(); ++it)
                    if (it->lowerCorner()(0) < 2 && it->lowerCorner()(1) < 2 && it->lowerCorner()(2) < 2)
                    {
                        ++nA;
                        CHECK_EQUAL(1u, it->level());
                        CHECK_EQUAL(0, it->lowerCorner()(0));
                        CHECK_EQUAL(0, it->lowerCorner()(1));
                        CHECK_EQUAL(0, it->lowerCorner()(2));
                    }
                CHECK_EQUAL(1, nA);
            }
        }
    }

    // Without the option, with the default value and with an explicit 0 the
    // marker selects as in rule 0 (groups A, B, C; 9 elements on the 2D fixture).
    TEST(CoarsenGroupRule_DefaultIsAnyChild)
    {
        gsTHBSplineBasis<2,real_t> basis = makeGroupBasis2D();
        const std::vector<real_t> err = groupErrors<2>(basis);
        const std::map<std::vector<index_t>, int> want = expectedGroupBoxes<2>(0);

        // default options, CoarsenGroupRule untouched
        gsOptionList def = Marker::defaultOptions();
        def.setInt("CoarsenRule", 1);
        def.setReal("CoarsenParam", 0.01);
        def.setSwitch("Admissible", false);
        CHECK_EQUAL(0, def.askInt("CoarsenGroupRule", -1));

        // hand-built list without the key
        gsOptionList bare;
        bare.addInt("CoarsenRule", "", 1);
        bare.addReal("CoarsenParam", "", 0.01);
        bare.addSwitch("Admissible", "", false);

        // explicit 0
        gsOptionList zero = def;
        zero.setInt("CoarsenGroupRule", 0);

        const gsOptionList * lists[3] = {&def, &bare, &zero};
        for (int k = 0; k != 3; ++k)
        {
            Marker m(basis, *lists[k]);
            m.setErrors(err);
            const Marker::HElementContainer res = m.markCrs();
            CHECK(crsBoxCounts<2>(m, res) == want);
            CHECK_EQUAL(9u, res.size());
        }
    }

    // A nonzero group rule is only defined for GARU (CoarsenRule 1).
    TEST(CoarsenGroupRule_RequiresGaru)
    {
        gsTHBSplineBasis<2,real_t> basis = makeGroupBasis2D();
        const std::vector<real_t> err = groupErrors<2>(basis);
        gsOptionList o = groupRuleOptions(Marker::defaultOptions(), 0);

        o.setInt("CoarsenRule", 3); o.setInt("CoarsenGroupRule", 2);
        CHECK_THROW(runCrs<2>(basis, o, err), std::runtime_error);
        o.setInt("CoarsenRule", 2); o.setInt("CoarsenGroupRule", 2);
        CHECK_THROW(runCrs<2>(basis, o, err), std::runtime_error);
        o.setInt("CoarsenRule", 3); o.setInt("CoarsenGroupRule", 1);
        CHECK_THROW(runCrs<2>(basis, o, err), std::runtime_error);

        o.setInt("CoarsenRule", 3); o.setInt("CoarsenGroupRule", 0);
        runCrs<2>(basis, o, err);
        o.setInt("CoarsenRule", 1); o.setInt("CoarsenGroupRule", 2);
        runCrs<2>(basis, o, err);
    }


    // One level-1 sibling group has a child refined further, so only three of
    // its four children are leaves; their errors (0.003) are far below the
    // threshold 0.01. Rule 0 selects those leaves individually; rules 1 and 2
    // require all four children to be leaves and select nothing.
    TEST(CoarsenGroupRule_IncompleteGroupExcluded2D)
    {
        gsTHBSplineBasis<2,real_t> basis = makeIncompleteGroupBasis2D();
        CHECK_EQUAL(22u, basis.numElements());
        const std::vector<real_t> err = incompleteGroupErrors(basis);

        for (index_t rule = 0; rule != 3; ++rule)
        {
            gsOptionList opts = groupRuleOptions(Marker::defaultOptions(), rule);
            Marker m(basis, opts);
            m.setErrors(err);
            const Marker::HElementContainer res = m.markCrs();
            if (rule == 0)
            {
                CHECK_EQUAL(3u, res.size());
                std::set<std::vector<index_t> > corners;
                for (Marker::HElementContainer::const_iterator it = res.begin(); it != res.end(); ++it)
                {
                    CHECK_EQUAL(1u, it->level());
                    corners.insert({it->lowerCorner()(0), it->lowerCorner()(1)});
                }
                const std::set<std::vector<index_t> > want = {{1,0}, {0,1}, {1,1}};
                CHECK(corners == want);
            }
            else
                CHECK_EQUAL(0u, res.size());
        }
    }

    // CoarsenGroupRule outside {0,1,2} is rejected even for GARU.
    TEST(CoarsenGroupRule_RejectsOutOfRange)
    {
        gsTHBSplineBasis<2,real_t> basis = makeGroupBasis2D();
        const std::vector<real_t> err = groupErrors<2>(basis);
        gsOptionList o = groupRuleOptions(Marker::defaultOptions(), 0);

        o.setInt("CoarsenGroupRule", 5);
        CHECK_THROW(runCrs<2>(basis, o, err), std::runtime_error);
        o.setInt("CoarsenGroupRule", -1);
        CHECK_THROW(runCrs<2>(basis, o, err), std::runtime_error);
    }
}
