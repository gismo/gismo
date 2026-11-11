/** @file gsHTree_test.cpp

    @brief Tests for gsHTree and the hierarchical (THB) bases built on it:
    leaf view, side queries, merge, polylines, index-level scaling, XML
    round trips, refinement/unrefinement, conversion and interface repair.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H. Verhelst
**/

#include "gismo_unittest.h"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <limits>
#include <sstream>
#include <vector>

using namespace gismo;

namespace
{
    typedef gsHTree<2>        tree2;
    typedef tree2::point      point2;

    /// A leaf (or box) with its level, ordered for order-independent comparison.
    struct Cell
    {
        index_t level;
        std::array<index_t, 2> lo;
        std::array<index_t, 2> up;

        bool operator==(const Cell & o) const
        { return level == o.level && lo == o.lo && up == o.up; }

        bool operator<(const Cell & o) const
        {
            if (level != o.level) return level < o.level;
            if (lo != o.lo) return lo < o.lo;
            return up < o.up;
        }
    };

    /// Sorted leaves in level-local indices, read through the leaf iterator.
    std::vector<Cell> leafList(const tree2 & t)
    {
        std::vector<Cell> r;
        for (tree2::const_literator it = t.beginLeafIterator(); it.good(); it.next())
        {
            Cell c;
            c.level = it.level();
            const point2 lo = it.lowerCorner(), up = it.upperCorner();
            for (size_t k = 0; k < 2; ++k)
            {
                c.lo[k] = lo[k];
                c.up[k] = up[k];
            }
            r.push_back(c);
        }
        std::sort(r.begin(), r.end());
        return r;
    }

    /// Sorted getBoxes() output (corners in the index space of getMaxInsLevel()).
    std::vector<Cell> boxSignatures(const tree2 & t)
    {
        gsMatrix<index_t> b1, b2;
        gsVector<index_t> lvl;
        t.getBoxes(b1, b2, lvl);
        std::vector<Cell> r;
        for (index_t i = 0; i < lvl.size(); ++i)
        {
            Cell c;
            c.level = lvl[i];
            for (size_t k = 0; k < 2; ++k)
            {
                c.lo[k] = b1(i, k);
                c.up[k] = b2(i, k);
            }
            r.push_back(c);
        }
        std::sort(r.begin(), r.end());
        return r;
    }

    /// Level of every unit cell of the level-L grid, row by row. Two trees
    /// describe the same hierarchical mesh iff these maps agree for an L not
    /// below either maximal level; the kd-splits, and hence getBoxes(), may
    /// differ.
    std::vector<int> levelMap(const tree2 & t, unsigned L)
    {
        const unsigned s = t.getIndexLevel() - L;
        const index_t nx = t.upperCorner()[0] >> s, ny = t.upperCorner()[1] >> s;
        std::vector<int> m;
        m.reserve(static_cast<size_t>(nx * ny));
        point2 lo, hi;
        for (index_t j = 0; j < ny; ++j)
            for (index_t i = 0; i < nx; ++i)
            {
                lo << i, j;
                hi << i + 1, j + 1;
                m.push_back(t.query4(lo, hi, static_cast<int>(L)));
            }
        return m;
    }

    /// Level along side s, sampled on the cells of the level-F grid (F not
    /// below getMaxInsLevel()). Cells not covered by any side box keep -1.
    std::vector<int> sideProfile(const tree2 & t, boxSide s, unsigned F)
    {
        gsMatrix<index_t> b1, b2;
        gsVector<index_t> lvl;
        t.getBoxesOnSide(s, b1, b2, lvl);
        const unsigned M = t.getMaxInsLevel();
        const short_t q = 1 - s.direction();
        const index_t n = t.upperCorner()[q] >> (t.getIndexLevel() - F);
        std::vector<int> p(static_cast<size_t>(n), -1);
        for (index_t i = 0; i < lvl.size(); ++i)
        {
            const index_t lo = b1(i, q) << (F - M), hi = b2(i, q) << (F - M);
            for (index_t c = lo; c < hi && c < n; ++c)
                p[static_cast<size_t>(c)] = static_cast<int>(lvl[i]);
        }
        return p;
    }

    bool fullyCovered(const std::vector<int> & p)
    { return std::find(p.begin(), p.end(), -1) == p.end(); }

    /// Sets level `lev` on cells [a,b) of a side profile.
    void setRange(std::vector<int> & p, size_t a, size_t b, int lev)
    { for (size_t c = a; c < b; ++c) p[c] = lev; }

    /// 2 x n matrix of deterministic points in (0,1)^2.
    gsMatrix<real_t> testPoints(index_t n)
    {
        gsMatrix<real_t> pts(2, n);
        for (index_t k = 0; k < n; ++k)
        {
            const real_t a = 0.6180339887498949 * (k + 1);
            const real_t b = 0.4142135623730951 * (k + 1);
            pts(0, k) = a - std::floor(a);
            pts(1, k) = b - std::floor(b);
        }
        return pts;
    }

    /// Degree 2 tensor basis with 8 x 8 elements.
    gsTensorBSplineBasis<2, real_t> tensor8x8()
    {
        gsKnotVector<real_t> kv(0, 1, 7, 3);
        return gsTensorBSplineBasis<2, real_t>(kv, kv);
    }

    /// Non-trivial planar control points, one row per function.
    gsMatrix<real_t> geomCoefs(index_t n)
    {
        gsMatrix<real_t> c(n, 2);
        for (index_t i = 0; i < n; ++i)
        {
            c(i, 0) = std::sin(0.7 * i) + i * 0.01;
            c(i, 1) = std::cos(1.3 * i);
        }
        return c;
    }

    /// The 8 x 8 THB basis refined by the given boxes (level-local indices).
    gsTHBSplineBasis<2> thb8x8(const std::vector<index_t> & boxes)
    {
        gsTHBSplineBasis<2> b(tensor8x8());
        if (!boxes.empty())
            b.refineElements(boxes);
        return b;
    }

    real_t maxDiff(const gsMatrix<real_t> & a, const gsMatrix<real_t> & b)
    {
        if (a.rows() != b.rows() || a.cols() != b.cols())
            return std::numeric_limits<real_t>::max();
        return (a - b).cwiseAbs().maxCoeff();
    }

    /// Tree of a basis that is a gsHTensorBasis<2> (e.g. a patch basis of a gsMultiPatch or gsMultiBasis).
    const tree2 & treeOf(const gsBasis<real_t> & b)
    { return static_cast<const gsHTensorBasis<2, real_t> &>(b).tree(); }

    std::string tempXml(const std::string & tag)
    {
        static int anchor = 0;
        std::ostringstream os;
        os << gsFileManager::getTempPath() << "/gsHTree_test_" << tag << "_"
           << reinterpret_cast<std::size_t>(&anchor) << "_"
           << std::chrono::high_resolution_clock::now().time_since_epoch().count()
           << ".xml";
        return os.str();
    }

    /// Tree with a level-1 box [2,2]-[6,6] and the nested level-2 box [6,6]-[10,10].
    tree2 nestedTree()
    {
        point2 upp, lo, hi;
        upp << 8, 8;
        tree2 t(upp);
        lo << 2, 2; hi << 6, 6;
        t.insertBox(lo, hi, 1);
        lo << 6, 6; hi << 10, 10;
        t.insertBox(lo, hi, 2);
        return t;
    }

    /// Every leaf of t: level-local corners agree with the global corners of
    /// the raw data and lie inside the grid of the leaf's level.
    template<class It>
    bool leafViewConsistent(const tree2 & t, It it)
    {
        bool ok = true;
        point2 g;
        for (; it.good(); it.next())
        {
            const unsigned lvl = static_cast<unsigned>(it.level());
            t.local2globalIndex(it.lowerCorner(), lvl, g);
            ok = ok && (g == it.data().lowerCorner());
            t.local2globalIndex(it.upperCorner(), lvl, g);
            ok = ok && (g == it.data().upperCorner());
            ok = ok && (it.level() == it.data().level());
            for (short_t k = 0; k < 2; ++k)
            {
                ok = ok && it.lowerCorner()[k] >= 0;
                ok = ok && it.upperCorner()[k] <= (index_t(8) << lvl);
            }
        }
        return ok;
    }
}

SUITE(gsHTree_test)
{
    TEST(LeafViewMatchesGlobalData)
    {
        tree2 t = nestedTree();
        const tree2 & ct = t;

        CHECK(leafViewConsistent(t, t.beginLeafIterator()));
        CHECK(leafViewConsistent(ct, ct.beginLeafIterator()));

        // The level-2 leaves cover exactly the box [6,6]-[10,10].
        index_t lo[2] = { std::numeric_limits<index_t>::max(), std::numeric_limits<index_t>::max() };
        index_t up[2] = { 0, 0 };
        index_t area = 0;
        for (tree2::const_literator it = ct.beginLeafIterator(); it.good(); it.next())
        {
            if (it.level() != 2) continue;
            const point2 l = it.lowerCorner(), u = it.upperCorner();
            for (short_t k = 0; k < 2; ++k)
            {
                lo[k] = (std::min)(lo[k], l[k]);
                up[k] = (std::max)(up[k], u[k]);
            }
            area += (u[0] - l[0]) * (u[1] - l[1]);
        }
        CHECK_EQUAL(6, lo[0]);
        CHECK_EQUAL(6, lo[1]);
        CHECK_EQUAL(10, up[0]);
        CHECK_EQUAL(10, up[1]);
        CHECK_EQUAL(16, area);
    }

    TEST(GetBoxesOnSideMatchesBruteForce)
    {
        tree2 t = nestedTree();
        point2 lo, hi;
        lo << 0, 0; hi << 4, 16;
        t.insertBox(lo, hi, 1);
        lo << 24, 0; hi << 32, 8;
        t.insertBox(lo, hi, 2);

        gsMatrix<index_t> a1, a2;
        gsVector<index_t> al;
        t.getBoxes(a1, a2, al);
        const unsigned shift = t.getIndexLevel() - t.getMaxInsLevel();

        const boundary::side sides[4] = { boundary::west, boundary::east, boundary::south, boundary::north };
        for (int si = 0; si < 4; ++si)
        {
            const boundary::side s = sides[si];
            const short_t q = (s - 1) / 2;
            gsMatrix<index_t> b1, b2;
            gsVector<index_t> bl;
            t.getBoxesOnSide(s, b1, b2, bl);

            std::vector<index_t> rows;
            for (index_t i = 0; i < al.size(); ++i)
            {
                const bool touches = (s % 2 == 1) ? (a1(i, q) == 0)
                    : ((a2(i, q) << shift) == t.upperCorner()[q]);
                if (touches) rows.push_back(i);
            }
            CHECK(!rows.empty());
            CHECK_EQUAL(rows.size(), static_cast<size_t>(bl.size()));
            if (rows.size() == static_cast<size_t>(bl.size()))
                for (size_t r = 0; r < rows.size(); ++r)
                {
                    CHECK_EQUAL(al[rows[r]], bl[r]);
                    for (short_t k = 0; k < 2; ++k)
                    {
                        CHECK_EQUAL(a1(rows[r], k), b1(r, k));
                        CHECK_EQUAL(a2(rows[r], k), b2(r, k));
                    }
                }
        }

        // Independent expectation on the 32 cells of the level-2 grid.
        const unsigned F = 2;
        std::vector<int> west(32, 1), east(32, 0), south(32, 0), north(32, 0);
        setRange(east, 0, 8, 2);
        setRange(south, 0, 8, 1);
        setRange(south, 24, 32, 2);
        setRange(north, 0, 8, 1);
        const std::vector<int> pw = sideProfile(t, boundary::west, F);
        const std::vector<int> pe = sideProfile(t, boundary::east, F);
        const std::vector<int> ps = sideProfile(t, boundary::south, F);
        const std::vector<int> pn = sideProfile(t, boundary::north, F);
        CHECK(fullyCovered(pw) && fullyCovered(pe) && fullyCovered(ps) && fullyCovered(pn));
        CHECK(west == pw);
        CHECK(east == pe);
        CHECK(south == ps);
        CHECK(north == pn);
    }

    TEST(MergeEqualsBothRefinements)
    {
        point2 upp, lo, hi;
        upp << 8, 8;
        tree2 A(upp), B(upp), C(upp);
        lo << 2, 2; hi << 6, 6;
        A.insertBox(lo, hi, 1);
        C.insertBox(lo, hi, 1);
        lo << 16, 16; hi << 24, 24;
        B.insertBox(lo, hi, 2);
        C.insertBox(lo, hi, 2);
        lo << 4, 4; hi << 8, 8;
        B.insertBox(lo, hi, 1);
        C.insertBox(lo, hi, 1);

        const std::vector<int> mapA = levelMap(A, 2), mapB = levelMap(B, 2);
        const tree2 r = tree2::merge(A, B);

        CHECK_EQUAL(C.getMaxInsLevel(), r.getMaxInsLevel());
        CHECK(levelMap(C, 2) == levelMap(r, 2));
        CHECK(mapA == levelMap(A, 2));
        CHECK(mapB == levelMap(B, 2));
    }

    TEST(PolylinesOneBox)
    {
        typedef std::vector<index_t> seg;
        point2 upp, lo, hi;
        upp << 8, 8;
        tree2 t(upp);

        // a segment is {x1,y1,x2,y2}
        const index_t full[4][4] = { {0,0,0,8}, {0,0,8,0}, {0,8,8,8}, {8,0,8,8} };
        std::vector<seg> expected(4);
        for (int i = 0; i < 4; ++i) expected[i].assign(full[i], full[i] + 4);

        std::vector<std::vector<std::vector<std::vector<index_t> > > > res = t.getPolylines();
        CHECK_EQUAL(1u, res.size());
        if (res.size() == 1 && res[0].size() == 1)
        {
            CHECK_EQUAL(4u, res[0][0].size());
            std::vector<seg> s = res[0][0];
            std::sort(s.begin(), s.end());
            CHECK(expected == s);
        }
        else
            CHECK(false);

        lo << 2, 2; hi << 6, 6;
        t.insertBox(lo, hi, 1);
        const index_t box[4][4] = { {2,2,2,6}, {2,2,6,2}, {2,6,6,6}, {6,2,6,6} };
        for (int i = 0; i < 4; ++i) expected[i].assign(box[i], box[i] + 4);

        res = t.getPolylines();
        CHECK_EQUAL(2u, res.size());
        if (res.size() == 2 && res[1].size() == 1)
        {
            CHECK_EQUAL(4u, res[1][0].size());
            std::vector<seg> s = res[1][0];
            std::sort(s.begin(), s.end());
            CHECK(expected == s);
        }
        else
            CHECK(false);
    }

    TEST(MultiplyDivideByTwo)
    {
        tree2 t = nestedTree();
        const tree2 copy(t);
        const point2 upp0 = copy.upperCorner();
        const std::vector<Cell> ref = leafList(copy);

        t.multiplyByTwo();
        CHECK(t.upperCorner() == 2 * upp0);
        CHECK_EQUAL(copy.getMaxInsLevel(), t.getMaxInsLevel());
        std::vector<Cell> doubled = ref;
        for (size_t i = 0; i < doubled.size(); ++i)
            for (size_t k = 0; k < 2; ++k)
            {
                doubled[i].lo[k] *= 2;
                doubled[i].up[k] *= 2;
            }
        CHECK(doubled == leafList(t));

        t.divideByTwo();
        CHECK(upp0 == t.upperCorner());
        CHECK(ref == leafList(t));
    }

    TEST(MovedFromTreeIsWalkable)
    {
        tree2 a = nestedTree();
        const size_t nLeaves = leafList(a).size();
        CHECK(nLeaves > 1);

        tree2 b(give(a));
        CHECK_EQUAL(nLeaves, leafList(b).size());

        // the moved-from tree has no leaves but can be walked, re-assigned and destroyed
        CHECK_EQUAL(0u, leafList(a).size());
        a = nestedTree();
        CHECK_EQUAL(nLeaves, leafList(a).size());

        tree2 c = nestedTree();
        tree2 d;
        d = give(c);
        CHECK_EQUAL(0u, leafList(c).size());
        CHECK_EQUAL(nLeaves, leafList(d).size());
        c = d;
        CHECK_EQUAL(nLeaves, leafList(c).size());

        gsTHBSplineBasis<2> b1 = thb8x8({1, 2, 2, 6, 6});
        const index_t n1 = b1.size();
        gsTHBSplineBasis<2> b2;
        b2 = give(b1);
        CHECK_EQUAL(n1, b2.size());
        CHECK_EQUAL(0u, leafList(b1.tree()).size());
        b1 = thb8x8({});
        CHECK_EQUAL(tensor8x8().size(), b1.size());
    }

    TEST(UniformRefineMultiLevel)
    {
        gsTHBSplineBasis<2> b = thb8x8({1, 2, 2, 6, 6, 2, 6, 6, 10, 10});
        const std::vector<Cell> before = boxSignatures(b.tree());
        const size_t nel = b.numElements();
        const index_t maxLevel = b.maxLevel();
        CHECK(before.size() > 1);

        b.uniformRefine();
        const std::vector<Cell> after = boxSignatures(b.tree());
        CHECK_EQUAL(maxLevel, b.maxLevel());
        CHECK_EQUAL(4 * nel, b.numElements());
        CHECK_EQUAL(before.size(), after.size());
        std::vector<Cell> doubled = before;
        for (size_t i = 0; i < doubled.size(); ++i)
            for (size_t k = 0; k < 2; ++k)
            {
                doubled[i].lo[k] *= 2;
                doubled[i].up[k] *= 2;
            }
        std::sort(doubled.begin(), doubled.end());
        CHECK(doubled == after);

        b.uniformCoarsen();
        CHECK(before == boxSignatures(b.tree()));
        CHECK_EQUAL(nel, b.numElements());
        CHECK_EQUAL(maxLevel, b.maxLevel());
    }

    TEST(XmlRoundTrip)
    {
        const gsMatrix<real_t> pts = testPoints(100);

        // (a) three levels
        {
            gsTHBSplineBasis<2> b = thb8x8({1, 2, 2, 6, 6, 2, 6, 6, 10, 10});
            gsTHBSpline<2> g(b, geomCoefs(b.size()));
            const std::string path = tempXml("a");
            gsFileData<> fw;
            fw.setFloatPrecision(17);
            fw << g;
            fw.dump(path);
            gsFileData<> fr(path);
            gsTHBSpline<2>::uPtr r = fr.getFirst< gsTHBSpline<2> >();
            std::remove(path.c_str());

            CHECK(r.get() != NULL);
            if (r.get() != NULL)
            {
                const unsigned L = g.basis().maxLevel();
                CHECK_EQUAL(g.basis().maxLevel(), r->basis().maxLevel());
                CHECK_EQUAL(g.basis().size(), r->basis().size());
                CHECK(levelMap(g.basis().tree(), L) == levelMap(r->basis().tree(), L));
                CHECK(maxDiff(g.eval(pts), r->eval(pts)) <= 1e-14);
                // getBoxes() of a re-inserted tree need not split identically
                CHECK(boxSignatures(g.basis().tree()) == boxSignatures(r->basis().tree()));
            }
        }

        // (b) the written boxes carry level-local indices
        {
            gsTHBSplineBasis<2> b = thb8x8({1, 2, 2, 6, 6});
            internal::gsXmlTree data;
            internal::gsXmlNode * nd = internal::gsXml< gsTHBSplineBasis<2, real_t> >::put(b, data);
            CHECK(nd != NULL);
            index_t lo[2] = { std::numeric_limits<index_t>::max(), std::numeric_limits<index_t>::max() };
            index_t up[2] = { 0, 0 };
            index_t area = 0;
            int nBoxes = 0;
            bool inRange = true;
            for (internal::gsXmlNode * bx = nd->first_node("box"); bx; bx = bx->next_sibling("box"))
            {
                const int L = atoi(bx->first_attribute("level")->value());
                std::istringstream ss(bx->value());
                index_t c[4];
                for (int k = 0; k < 4; ++k) ss >> c[k];
                inRange = inRange && c[0] >= 0 && c[0] < c[2] && c[2] <= (index_t(8) << L)
                                  && c[1] >= 0 && c[1] < c[3] && c[3] <= (index_t(8) << L);
                CHECK_EQUAL(1, L);
                lo[0] = (std::min)(lo[0], c[0]); lo[1] = (std::min)(lo[1], c[1]);
                up[0] = (std::max)(up[0], c[2]); up[1] = (std::max)(up[1], c[3]);
                area += (c[2] - c[0]) * (c[3] - c[1]);
                ++nBoxes;
            }
            CHECK(nBoxes > 0);
            CHECK(inRange);
            CHECK_EQUAL(2, lo[0]);
            CHECK_EQUAL(2, lo[1]);
            CHECK_EQUAL(6, up[0]);
            CHECK_EQUAL(6, up[1]);
            CHECK_EQUAL(16, area);
        }

        // (c) a data-file geometry with geoDim 3
        {
            gsFileData<> fd("surfaces/thbs_face_3levels.xml");
            gsTHBSpline<2>::uPtr g = fd.getFirst< gsTHBSpline<2> >();
            CHECK(g.get() != NULL);
            if (g.get() != NULL)
            {
                const std::string path = tempXml("c");
                gsFileData<> fw;
                fw.setFloatPrecision(17);
                fw << *g;
                fw.dump(path);
                gsFileData<> fr(path);
                gsTHBSpline<2>::uPtr r = fr.getFirst< gsTHBSpline<2> >();
                std::remove(path.c_str());

                CHECK(r.get() != NULL);
                if (r.get() != NULL)
                {
                    gsMatrix<real_t> supp = g->support();
                    gsMatrix<real_t> p = testPoints(100);
                    for (index_t k = 0; k < p.cols(); ++k)
                        for (index_t i = 0; i < 2; ++i)
                            p(i, k) = supp(i, 0) + (supp(i, 1) - supp(i, 0)) * p(i, k);

                    const unsigned L = g->basis().maxLevel();
                    CHECK_EQUAL(3, g->geoDim());
                    CHECK_EQUAL(g->basis().maxLevel(), r->basis().maxLevel());
                    CHECK_EQUAL(g->basis().size(), r->basis().size());
                    CHECK(levelMap(g->basis().tree(), L) == levelMap(r->basis().tree(), L));
                    CHECK(maxDiff(g->eval(p), r->eval(p)) <= 1e-14);
                    CHECK(boxSignatures(g->basis().tree()) == boxSignatures(r->basis().tree()));
                }
            }
        }
    }

    TEST(RefineUnrefineIdentity)
    {
        const gsMatrix<real_t> pts = testPoints(100);
        const gsTHBSplineBasis<2> fresh = thb8x8({});
        const gsMatrix<real_t> freshEval = fresh.eval(pts);

        // (a) a single refinement is undone completely
        {
            gsTHBSplineBasis<2> b = thb8x8({1, 2, 2, 6, 6});
            CHECK_EQUAL(1, b.maxLevel());
            b.unrefineElements({0, 1, 1, 3, 3});
            CHECK_EQUAL(0, b.maxLevel());
            CHECK_EQUAL(0u, b.tree().getMaxInsLevel());
            CHECK_EQUAL(1, b.tree().leafSize());
            CHECK_EQUAL(fresh.size(), b.size());
            CHECK(maxDiff(freshEval, b.eval(pts)) <= 1e-14);
        }

        // (b) nested refinements are undone level by level
        {
            const gsTHBSplineBasis<2> ref1 = thb8x8({1, 2, 2, 6, 6});
            const gsMatrix<real_t> ref1Eval = ref1.eval(pts);
            gsTHBSplineBasis<2> b = thb8x8({1, 2, 2, 6, 6});
            b.refineElements({2, 6, 6, 10, 10});
            CHECK_EQUAL(2, b.maxLevel());
            b.unrefineElements({1, 3, 3, 5, 5});
            CHECK_EQUAL(1, b.maxLevel());
            CHECK(levelMap(ref1.tree(), 1) == levelMap(b.tree(), 1));
            CHECK_EQUAL(ref1.size(), b.size());
            CHECK(maxDiff(ref1Eval, b.eval(pts)) <= 1e-14);

            b.unrefineElements({0, 1, 1, 3, 3});
            CHECK_EQUAL(0, b.maxLevel());
            CHECK_EQUAL(fresh.size(), b.size());
            CHECK(maxDiff(freshEval, b.eval(pts)) <= 1e-14);
        }
    }

    TEST(UniformRefineCoarsenWithCoefsKeepGeometry)
    {
        const gsMatrix<real_t> pts = testPoints(100);
        gsTHBSplineBasis<2> b = thb8x8({1, 2, 2, 6, 6});
        const gsTHBSpline<2> G(b, geomCoefs(b.size()));
        const gsMatrix<real_t> ref = G.eval(pts);
        const std::vector<Cell> sig = boxSignatures(G.basis().tree());

        gsTHBSpline<2> copy(G);
        copy.basis().uniformRefine_withCoefs(copy.coefs());
        CHECK_EQUAL(copy.basis().size(), copy.coefs().rows());
        CHECK(copy.basis().size() > G.basis().size());
        CHECK(maxDiff(ref, copy.eval(pts)) <= 1e-10);

        copy.basis().uniformCoarsen_withCoefs(copy.coefs());
        CHECK_EQUAL(G.basis().size(), copy.coefs().rows());
        CHECK(maxDiff(ref, copy.eval(pts)) <= 1e-10);
        CHECK(sig == boxSignatures(copy.basis().tree()));
    }

    TEST(ConvertToBSpline)
    {
        const gsMatrix<real_t> pts = testPoints(100);
        gsTHBSplineBasis<2> b = thb8x8({1, 2, 2, 6, 6});
        const gsTHBSpline<2> G(b, geomCoefs(b.size()));
        const gsMatrix<real_t> ref = G.eval(pts);

        gsTHBSpline<2> copy(G);
        gsTensorBSpline<2, real_t> tb;
        copy.convertToBSpline(tb);
        CHECK(tb.basis().size() > G.basis().size());
        CHECK(maxDiff(ref, tb.eval(pts)) <= 1e-10);
    }

    TEST(RepairInterfaceTwoPatches)
    {
        // (a) two 4 x 4 patches sharing the side x = 1
        {
            gsKnotVector<real_t> kv(0, 1, 3, 3);
            gsTHBSplineBasis<2> basis((gsTensorBSplineBasis<2, real_t>(kv, kv)));
            gsMatrix<real_t> c0 = basis.anchors().transpose();
            gsMatrix<real_t> c1 = c0;
            c1.col(0).array() += 1;

            gsMultiPatch<real_t> mp;
            mp.addPatch(gsTHBSpline<2>(basis, c0));
            mp.addPatch(gsTHBSpline<2>(basis, c1));
            mp.computeTopology();
            CHECK_EQUAL(1u, mp.interfaces().size());
            const boundaryInterface bi = mp.interfaces().front();

            mp.patch(0).basis().refineElements_withCoefs(mp.patch(0).coefs(), {1, 6, 0, 8, 4});

            const unsigned F = 1;

            // multi-basis path, on bases that still have the mismatch
            gsMultiBasis<real_t> mb(mp);
            CHECK(mb.repairInterface(bi));
            CHECK(!mb.repairInterface(bi));
            CHECK(sideProfile(treeOf(mb.basis(bi.first().patch)), bi.first().side(), F)
                  == sideProfile(treeOf(mb.basis(bi.second().patch)), bi.second().side(), F));

            // multi-patch path
            CHECK(mp.repairInterface(bi));
            CHECK(!mp.repairInterface(bi));
            const tree2 & t1 = treeOf(mp.patch(bi.first().patch).basis());
            const tree2 & t2 = treeOf(mp.patch(bi.second().patch).basis());
            const std::vector<int> p1 = sideProfile(t1, bi.first().side(), F);
            const std::vector<int> p2 = sideProfile(t2, bi.second().side(), F);
            CHECK(fullyCovered(p1));
            CHECK(p1 == p2);
            CHECK(std::find(p1.begin(), p1.end(), 1) != p1.end());
        }

        // (b) three patches of an L-shape
        {
            gsMultiPatch<real_t> mp;
            gsReadFile<>("planar/lshape2d_3patches_thb.xml", mp);
            CHECK_EQUAL(3u, mp.nPatches());
            mp.patch(0).basis().refineElements_withCoefs(mp.patch(0).coefs(), {2, 0, 0, 1, 1});

            const std::vector<boundaryInterface> & ifs = mp.interfaces();
            const boundaryInterface * west0 = NULL;
            for (size_t i = 0; i < ifs.size(); ++i)
                if ((ifs[i].first().patch == 0 && ifs[i].first().side() == boundary::west)
                    || (ifs[i].second().patch == 0 && ifs[i].second().side() == boundary::west))
                    west0 = &ifs[i];
            CHECK(west0 != NULL);
            if (west0 != NULL)
            {
                const boundaryInterface it = *west0;
                CHECK(mp.repairInterface(it));
                CHECK(!mp.repairInterface(it));
            }

            mp.repairInterfaces();
            for (size_t i = 0; i < ifs.size(); ++i)
            {
                const boundaryInterface bi = ifs[i];
                CHECK(!mp.repairInterface(bi));
                const tree2 & t1 = treeOf(mp.patch(bi.first().patch).basis());
                const tree2 & t2 = treeOf(mp.patch(bi.second().patch).basis());
                const unsigned F = (std::max)(t1.getMaxInsLevel(), t2.getMaxInsLevel());
                const std::vector<int> p1 = sideProfile(t1, bi.first().side(), F);
                const std::vector<int> p2 = sideProfile(t2, bi.second().side(), F);
                CHECK(fullyCovered(p1));
                CHECK(p1 == p2);
            }
        }
    }
}
