/** @file gsDofMapper_test.cpp

    @brief Tests for gsDofMapper

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst, M.M. Hodzelmans
*/

#include "gismo_unittest.h"

#include <algorithm>
#include <cstring>
#include <functional>
#include <map>
#include <new>
#include <sstream>

namespace
{
typedef std::vector<index_t> vec;
typedef std::vector<std::pair<index_t,index_t> > pairs;

/*
  Fixture M2: two patches (4 and 5 local dofs), two components. Both
  components couple dof 3 of patch 0 with dof 0 of patch 1; component 1
  additionally couples dofs 1 and 2 of patch 1. Dof 0 of patch 0 (component 0)
  and dof 4 of patch 1 (component 1) are eliminated.

  Numbering after finalize() at shift 0. A shift s adds s to every global
  index; stored values (tags, anyPreImages positions) and patch-local
  indices do not shift.

  comp | patch | local dof -> global index           | kind
  -----+-------+-------------------------------------+------------------------------
   0   |   0   | 0->13, 1->0, 2->1, 3->6             | 0 eliminated; 3 coupled
   0   |   1   | 0->6, 1->2, 2->3, 3->4, 4->5        | 0 coupled
   1   |   0   | 0->7, 1->8, 2->9, 3->11             | 3 coupled
   1   |   1   | 0->11, 1->12, 2->12, 3->10, 4->14   | 0, 1, 2 coupled; 4 eliminated

  Free indices are 0-12 (component 0 owns 0-6, component 1 owns 7-12);
  eliminated indices are 13 (component 0) and 14 (component 1). The coupled
  band is {6} for component 0 and {11,12} for component 1.
  freeSize() = 13, boundarySize() = 2, size() = 15, coupledSize() = 3,
  firstIndex(0) = 0, firstIndex(1) = 7.
*/
gsDofMapper makeM2()
{
    gsVector<index_t> sz(2); sz << 4, 5;
    gsDofMapper m(sz, 2);
    m.matchDof(0,3, 1,0, 0);
    m.matchDof(0,3, 1,0, 1);
    m.matchDof(1,1, 1,2, 1);
    m.eliminateDof(0,0, 0);
    m.eliminateDof(4,1, 1);
    m.finalize();
    return m;
}

vec toVec(const gsVector<index_t> & v) { return vec(v.data(), v.data() + v.size()); }
} // namespace

namespace {

typedef std::pair<index_t,index_t> dofPair;

gsMatrix<index_t> idx(std::initializer_list<index_t> l)
{
    gsMatrix<index_t> r(static_cast<index_t>(l.size()), 1);
    index_t i = 0;
    for (index_t v : l) r(i++, 0) = v;
    return r;
}

std::string join(const std::vector<dofPair> & v)
{
    std::ostringstream os;
    for (const dofPair & p : v)
        os << "(" << p.first << "," << p.second << ")";
    return os.str();
}

// "off=offset(k)... ps=patchSize(k,0)...;patchSize(k,1)..."
std::string layout(const gsDofMapper & m)
{
    std::ostringstream os;
    os << "off=";
    for (size_t k = 0; k != m.numPatches(); ++k)
        os << (k ? " " : "") << m.offset(k);
    os << " ps=";
    for (index_t c = 0; c != m.numComponents(); ++c)
        for (size_t k = 0; k != m.numPatches(); ++k)
            os << (k ? " " : (c ? ";" : "")) << m.patchSize(k, c);
    return os.str();
}

// index(i,k,c) of every local dof, as "c0k0|c0k1;c1k0|c1k1"
std::string indices(const gsDofMapper & m)
{
    std::ostringstream os;
    for (index_t c = 0; c != m.numComponents(); ++c)
        for (size_t k = 0; k != m.numPatches(); ++k)
        {
            if (c || k) os << (k ? "|" : ";");
            for (index_t i = 0; i != static_cast<index_t>(m.patchSize(k, c)); ++i)
                os << (i ? "," : "") << m.index(i, k, c);
        }
    return os.str();
}

// Counters, layout, index/bindex/coupling of every local dof, and preimages
std::string state(const gsDofMapper & m)
{
    std::ostringstream os;
    m.print(os);
    os << layout(m) << " map=" << m.mapSize() << "\n";
    for (index_t c = 0; c != m.numComponents(); ++c)
    {
        os << "c" << c << ":" << m.size(c) << "," << m.freeSize(c) << "," << m.firstIndex(c);
        for (size_t k = 0; k != m.numPatches(); ++k)
            for (index_t i = 0; i != static_cast<index_t>(m.patchSize(k, c)); ++i)
            {
                os << " " << m.index(i, k, c);
                if (m.is_boundary(i, k, c)) os << "b" << m.bindex(i, k, c);
                if (m.is_coupled(i, k, c)) os << "c";
            }
        os << "\n";
    }
    std::vector<dofPair> pre;
    for (index_t gl = 0; gl != m.size(); ++gl)
    {
        m.preImage(gl, pre);
        os << join(pre) << " ";
    }
    return os.str();
}

// indexOnPatch(gl,k) per patch, '-' where gl is not on patch k
std::string onPatches(const gsDofMapper & m, index_t gl)
{
    std::ostringstream os;
    for (size_t k = 0; k != m.numPatches(); ++k)
    {
        index_t loc = -1;
        os << (k ? " " : "") << (m.indexOnPatch(gl, k, loc) ? std::to_string(loc) : std::string("-"));
    }
    return os.str();
}

// Matching and elimination on the component with the given index
typedef std::function<void(gsDofMapper &, index_t)> compOps;

// Builds gsDofMapper(sizes) with ops[c] on component c and compares each
// component with the scalar mapper of column c with ops[c]. Its free dofs
// come after those of the preceding components; its eliminated dofs come
// after all free dofs and the eliminated dofs of the preceding components.
gsDofMapper checkAgainstScalars(const gsMatrix<index_t> & sizes, const std::vector<compOps> & ops)
{
    const index_t nc = sizes.cols(), np = sizes.rows();
    gsDofMapper m(sizes);
    for (index_t c = 0; c != nc; ++c)
        ops[c](m, c);
    m.finalize();

    std::vector<gsDofMapper> s;
    index_t F = 0, E = 0, cpld = 0;
    size_t map = 0;
    for (index_t c = 0; c != nc; ++c)
    {
        const gsVector<index_t> col = sizes.col(c);
        s.push_back(gsDofMapper(col));
        ops[c](s.back(), 0);
        s.back().finalize();
        F += s.back().freeSize();
        E += s.back().boundarySize();
        cpld += s.back().coupledSize();
        map += s.back().mapSize();
    }
    std::ostringstream expected, actual;
    expected << "free=" << F << " elim=" << E << " size=" << F + E << " cpld=" << cpld << " map=" << map;
    actual << "free=" << m.freeSize() << " elim=" << m.boundarySize() << " size=" << m.size()
           << " cpld=" << m.coupledSize() << " map=" << m.mapSize();
    std::vector<dofPair> pre;
    index_t Fb = 0, Eb = 0, n = 0;
    for (index_t c = 0; c != nc; ++c)
    {
        const gsDofMapper & sc = s[c];
        const gsVector<index_t> asVec = m.asVector(c);
        expected << "\nc" << c << " size=" << sc.size() << " free=" << sc.freeSize() << " total=" << sc.totalSize();
        actual   << "\nc" << c << " size=" << m.size(c) << " free=" << m.freeSize(c) << " total=" << m.totalSize(c);
        for (index_t k = 0, j = 0; k != np; ++k)
        {
            expected << "\n k" << k << " ps=" << sc.patchSize(k);
            actual   << "\n k" << k << " ps=" << m.patchSize(k, c);
            if (sc.patchSize(k) != m.patchSize(k, c))
                continue;
            for (index_t i = 0; i != static_cast<index_t>(sc.patchSize(k)); ++i, ++j, ++n)
            {
                const index_t g = sc.index(i, k), gl = m.index(i, k, c);
                const index_t e = sc.is_free_index(g) ? Fb + g : F + Eb + g - sc.freeSize();
                sc.preImage(g, pre);
                expected << " " << e << "c" << c << "/" << e << "/" << e << join(pre) << onPatches(sc, g);
                m.preImage(gl, pre);
                actual << " " << gl << "c" << m.componentOf(gl) << "/" << m.mapIndex(n) << "/" << asVec[j]
                       << join(pre) << onPatches(m, gl);
            }
        }
        Fb += sc.freeSize();
        Eb += sc.boundarySize();
    }
    CHECK_EQUAL(expected.str(), actual.str());
    return m;
}

// Raviart-Thomas type pair of degree 2 on nPatches unit squares side by side,
// nx x ny elements each: component c has degree 3 in direction c, 2 in the
// other. Component 0 is eliminated on the west side of patch 0, component 1
// on the south side of every patch.
struct rtFixture
{
    gsMatrix<index_t> sizes;
    std::vector<compOps> ops;

    rtFixture(unsigned nx, unsigned ny, index_t nPatches) : sizes(nPatches, 2)
    {
        const std::vector<boundaryInterface> ifs =
            gsNurbsCreator<real_t>::BSplineSquareGrid(nPatches, 1, 1.0).interfaces();
        for (index_t c = 0; c != 2; ++c)
        {
            const gsTensorBSplineBasis<2,real_t> bc(gsKnotVector<real_t>(0, 1, nx - 1, c == 0 ? 4 : 3),
                                                    gsKnotVector<real_t>(0, 1, ny - 1, c == 1 ? 4 : 3));
            sizes.col(c).setConstant(bc.size());
            ops.push_back([bc, ifs, c, nPatches](gsDofMapper & m, index_t comp) {
                gsMatrix<index_t> b1, b2;
                for (const boundaryInterface & bi : ifs)
                {
                    bc.matchWith(bi, bc, b1, b2);
                    m.matchDofs(bi.first().patch, b1, bi.second().patch, b2, comp);
                }
                if (0 == c)
                    m.markBoundary(0, bc.boundary(boundary::west), comp);
                else
                    for (index_t k = 0; k != nPatches; ++k)
                        m.markBoundary(k, bc.boundary(boundary::south), comp);
            });
        }
    }
};

} // anonymous namespace

SUITE(gsDofMapper_test)
{
    TEST(DofMapperT1_shift0)
    {
        gsDofMapper m = makeM2();

        const vec as0 = {13,0,1,6, 6,2,3,4,5};
        const vec as1 = {7,8,9,11, 11,12,12,10,14};
        CHECK(toVec(m.asVector(0)) == as0);
        CHECK(toVec(m.asVector(1)) == as1);
        CHECK_EQUAL(13, m.freeSize());
        CHECK_EQUAL(15, m.size());
        CHECK_EQUAL(2, m.boundarySize());

        const index_t comp[15] = {0,0,0,0,0,0,0, 1,1,1,1,1,1, 0,1};
        bool componentsOk = true;
        for (index_t g = 0; g < 15; ++g)
        {
            const index_t c = m.componentOf(g);
            CHECK_EQUAL(comp[g], c);
            componentsOk = componentsOk && (c == comp[g]);
        }

        const vec free0 = {1,2,3};
        const vec free1 = {0,1,2,3,4};
        const vec bnd0 = {0};
        CHECK(toVec(m.findFree(0)) == free0);
        CHECK(toVec(m.findFree(1)) == free1);
        CHECK(toVec(m.findBoundary(0)) == bnd0);
        CHECK_EQUAL(0, (index_t)m.findBoundary(1).size());

        // preImage(), anyPreImage() and indexOnPatch() search the component reported by componentOf()
        if (!componentsOk) return;

        const pairs pre[15] = {
            {{0,1}}, {{0,2}}, {{1,1}}, {{1,2}}, {{1,3}}, {{1,4}}, {{0,3},{1,0}},
            {{0,0}}, {{0,1}}, {{0,2}}, {{1,3}}, {{0,3},{1,0}}, {{1,1},{1,2}},
            {{0,0}}, {{1,4}} };
        for (index_t g = 0; g < 15; ++g)
        {
            pairs got;
            m.preImage(g, got);
            CHECK(got == pre[g]);
        }

        const std::pair<index_t,index_t> any0(0,1), any6(0,3), any12(1,1), any13(0,0), any14(1,4);
        CHECK(m.anyPreImage(0) == any0);
        CHECK(m.anyPreImage(6) == any6);
        CHECK(m.anyPreImage(12) == any12);
        CHECK(m.anyPreImage(13) == any13);
        CHECK(m.anyPreImage(14) == any14);

        index_t l = -1;
        CHECK(m.indexOnPatch(6, 0, l));
        CHECK_EQUAL(3, l);
        CHECK(m.indexOnPatch(6, 1, l));
        CHECK_EQUAL(0, l);
        CHECK(m.indexOnPatch(12, 1, l));
        CHECK_EQUAL(1, l);
        CHECK(m.indexOnPatch(13, 0, l));
        CHECK_EQUAL(0, l);
        CHECK(m.indexOnPatch(14, 1, l));
        CHECK_EQUAL(4, l);
        CHECK(!m.indexOnPatch(0, 1, l));
    }

    TEST(DofMapperT1_shift13)
    {
        gsDofMapper m = makeM2();
        m.setShift(13);

        const vec as0 = {26,13,14,19, 19,15,16,17,18};
        const vec as1 = {20,21,22,24, 24,25,25,23,27};
        CHECK(toVec(m.asVector(0)) == as0);
        CHECK(toVec(m.asVector(1)) == as1);
        CHECK_EQUAL(13, m.freeSize());
        CHECK_EQUAL(15, m.size());
        CHECK_EQUAL(2, m.boundarySize());

        const index_t comp[15] = {0,0,0,0,0,0,0, 1,1,1,1,1,1, 0,1};
        bool componentsOk = true;
        for (index_t g = 0; g < 15; ++g)
        {
            const index_t c = m.componentOf(13 + g);
            CHECK_EQUAL(comp[g], c);
            componentsOk = componentsOk && (c == comp[g]);
        }

        const vec free0 = {1,2,3};
        const vec free1 = {0,1,2,3,4};
        const vec bnd0 = {0};
        CHECK(toVec(m.findFree(0)) == free0);
        CHECK(toVec(m.findFree(1)) == free1);
        CHECK(toVec(m.findBoundary(0)) == bnd0);
        CHECK_EQUAL(0, (index_t)m.findBoundary(1).size());

        // preImage(), anyPreImage() and indexOnPatch() search the component reported by componentOf()
        if (!componentsOk) return;

        const pairs pre[15] = {
            {{0,1}}, {{0,2}}, {{1,1}}, {{1,2}}, {{1,3}}, {{1,4}}, {{0,3},{1,0}},
            {{0,0}}, {{0,1}}, {{0,2}}, {{1,3}}, {{0,3},{1,0}}, {{1,1},{1,2}},
            {{0,0}}, {{1,4}} };
        for (index_t g = 0; g < 15; ++g)
        {
            pairs got;
            m.preImage(13 + g, got);
            CHECK(got == pre[g]);
        }

        const std::pair<index_t,index_t> any0(0,1), any6(0,3), any12(1,1), any13(0,0), any14(1,4);
        CHECK(m.anyPreImage(13 + 0) == any0);
        CHECK(m.anyPreImage(13 + 6) == any6);
        CHECK(m.anyPreImage(13 + 12) == any12);
        CHECK(m.anyPreImage(13 + 13) == any13);
        CHECK(m.anyPreImage(13 + 14) == any14);

        index_t l = -1;
        CHECK(m.indexOnPatch(13 + 6, 0, l));
        CHECK_EQUAL(3, l);
        CHECK(m.indexOnPatch(13 + 6, 1, l));
        CHECK_EQUAL(0, l);
        CHECK(m.indexOnPatch(13 + 12, 1, l));
        CHECK_EQUAL(1, l);
        CHECK(m.indexOnPatch(13 + 13, 0, l));
        CHECK_EQUAL(0, l);
        CHECK(m.indexOnPatch(13 + 14, 1, l));
        CHECK_EQUAL(4, l);
        CHECK(!m.indexOnPatch(13 + 0, 1, l));
    }

    TEST(DofMapperT2_shift0)
    {
        gsDofMapper m = makeM2();

        const vec cpl0 = {3};
        const vec cpl1 = {0};
        const vec fu0 = {1,2};
        const vec fu1 = {1,2,3,4};
        CHECK(toVec(m.findCoupled(0)) == cpl0);
        CHECK(toVec(m.findCoupled(1)) == cpl1);
        CHECK(toVec(m.findFreeUncoupled(0)) == fu0);
        CHECK(toVec(m.findFreeUncoupled(1)) == fu1);

        const bool coupled[15] = {0,0,0,0,0,0,1, 0,0,0,0,1,1, 0,0};
        for (index_t g = 0; g < 15; ++g)
            CHECK_EQUAL(coupled[g], m.is_coupled_index(g));

        CHECK(m.is_coupled(3,0,0));
        CHECK(m.is_coupled(0,1,0));
        CHECK(m.is_coupled(3,0,1));
        CHECK(m.is_coupled(0,1,1));
        CHECK(m.is_coupled(1,1,1));
        CHECK(m.is_coupled(2,1,1));
        CHECK(!m.is_coupled(1,0,0));
        CHECK(!m.is_coupled(3,1,1));
        CHECK(!m.is_coupled(0,0,0));
        CHECK(!m.is_coupled(4,1,1));

        gsDofMapper t = makeM2();
        t.markCoupledAsTagged();
        const bool tagged[15] = {0,0,0,0,0,0,1, 0,0,0,0,1,1, 0,0};
        for (index_t g = 0; g < 15; ++g)
            CHECK_EQUAL(tagged[g], t.is_tagged_index(g));
        CHECK_EQUAL(3, t.taggedSize());
        CHECK_EQUAL(3, t.coupledSize());
        CHECK(t.is_tagged(3,0,0));
        CHECK(t.is_tagged(0,1,0));
        CHECK(t.is_tagged(3,0,1));
        CHECK(t.is_tagged(0,1,1));
        CHECK(t.is_tagged(1,1,1));
        CHECK(t.is_tagged(2,1,1));
        CHECK(!t.is_tagged(3,1,1));
        CHECK(!t.is_tagged(0,0,0));
        CHECK(!t.is_tagged(4,1,1));
    }

    TEST(DofMapperT2_shift13)
    {
        gsDofMapper m = makeM2();
        m.setShift(13);

        const vec cpl0 = {3};
        const vec cpl1 = {0};
        const vec fu0 = {1,2};
        const vec fu1 = {1,2,3,4};
        CHECK(toVec(m.findCoupled(0)) == cpl0);
        CHECK(toVec(m.findCoupled(1)) == cpl1);
        CHECK(toVec(m.findFreeUncoupled(0)) == fu0);
        CHECK(toVec(m.findFreeUncoupled(1)) == fu1);

        const bool coupled[15] = {0,0,0,0,0,0,1, 0,0,0,0,1,1, 0,0};
        for (index_t g = 0; g < 15; ++g)
            CHECK_EQUAL(coupled[g], m.is_coupled_index(13 + g));

        CHECK(m.is_coupled(3,0,0));
        CHECK(m.is_coupled(0,1,0));
        CHECK(m.is_coupled(3,0,1));
        CHECK(m.is_coupled(0,1,1));
        CHECK(m.is_coupled(1,1,1));
        CHECK(m.is_coupled(2,1,1));
        CHECK(!m.is_coupled(1,0,0));
        CHECK(!m.is_coupled(3,1,1));
        CHECK(!m.is_coupled(0,0,0));
        CHECK(!m.is_coupled(4,1,1));
    }

    TEST(DofMapperT2_oneComponent)
    {
        gsVector<index_t> s1(1); s1 << 3;
        gsDofMapper m(s1);
        m.markCoupled(2,0);
        m.finalize();
        CHECK(m.is_coupled(2,0,0));
        m.markCoupledAsTagged();
        CHECK(m.is_tagged(2,0,0));
        CHECK(!m.is_tagged(0,0,0));
        CHECK(!m.is_tagged(1,0,0));
        CHECK_EQUAL(1, m.taggedSize());
    }

    TEST(DofMapperT3_taggedShift)
    {
        // Tag first, then shift
        gsDofMapper m = makeM2();
        m.markTagged(1,1,1);
        m.setShift(13);
        CHECK(m.is_tagged(1,1,1));
        CHECK(m.is_tagged(2,1,1));
        CHECK(m.is_tagged_index(13 + 12));
        CHECK(!m.is_tagged(3,1,1));
        CHECK_EQUAL(0, m.tindex(1,1,1));
        CHECK_EQUAL(1, m.taggedSize());

        // Shift first, then tag
        gsDofMapper f = makeM2();
        f.setShift(13);
        f.markTagged(1,1,1);
        CHECK_EQUAL(1, f.taggedSize());
        CHECK(f.is_tagged_index(13 + 12));
        CHECK(f.is_tagged(1,1,1));

        // Shift first, then tag all coupled dofs
        gsDofMapper g = makeM2();
        g.setShift(13);
        g.markCoupledAsTagged();
        CHECK_EQUAL(3, g.taggedSize());
        CHECK(g.is_tagged_index(13 + 6));
        CHECK(g.is_tagged_index(13 + 11));
        CHECK(g.is_tagged_index(13 + 12));
        CHECK(g.is_tagged(3,0,0));
        CHECK(g.is_tagged(0,1,0));
        CHECK(g.is_tagged(3,0,1));
        CHECK(g.is_tagged(0,1,1));
        CHECK(g.is_tagged(1,1,1));
        CHECK(g.is_tagged(2,1,1));
        CHECK(!g.is_tagged(3,1,1));
        CHECK(!g.is_tagged(1,0,0));
    }

    TEST(DofMapperT4_shift0)
    {
        gsDofMapper m = makeM2();
        m.markTagged(1,0,0);
        m.markTagged(2,1,0);
        m.markCoupledAsTagged();

        // Stored tags are {0,3,6,11,12}
        const vec t0 = {1,3};
        const vec t1 = {0,2};
        CHECK(toVec(m.findTagged(0)) == t0);
        CHECK(toVec(m.findTagged(1)) == t1);
    }

    TEST(DofMapperT4_shift13)
    {
        gsDofMapper m = makeM2();
        m.markTagged(1,0,0);
        m.markTagged(2,1,0);
        m.markCoupledAsTagged();
        m.setShift(13);

        // Stored tags are {0,3,6,11,12}
        const vec t0 = {1,3};
        const vec t1 = {0,2};
        CHECK(toVec(m.findTagged(0)) == t0);
        CHECK(toVec(m.findTagged(1)) == t1);
    }

    TEST(DofMapperT5_shift0)
    {
        // Same matches as M2 but nothing eliminated: 8 + 7 dofs
        gsVector<index_t> sz(2); sz << 4, 5;
        gsDofMapper noElim(sz, 2);
        noElim.matchDof(0,3, 1,0, 0);
        noElim.matchDof(0,3, 1,0, 1);
        noElim.matchDof(1,1, 1,2, 1);
        noElim.finalize();
        CHECK_EQUAL(15, (index_t)noElim.anyPreImages(0).size());

        gsDofMapper m = makeM2();

        const pairs want0 = {{0,1},{0,2},{1,1},{1,2},{1,3},{1,4},{0,3},
                             {-1,0},{-1,0},{-1,0},{-1,0},{-1,0},{-1,0},
                             {0,0},{-1,0}};
        const pairs want1 = {{-1,0},{-1,0},{-1,0},{-1,0},{-1,0},{-1,0},{-1,0},
                             {0,0},{0,1},{0,2},{1,3},{0,3},{1,1},
                             {-1,0},{1,4}};
        CHECK(m.anyPreImages(0) == want0);
        CHECK(m.anyPreImages(1) == want1);
    }

    TEST(DofMapperT5_shift13)
    {
        gsVector<index_t> sz(2); sz << 4, 5;
        gsDofMapper noElim(sz, 2);
        noElim.matchDof(0,3, 1,0, 0);
        noElim.matchDof(0,3, 1,0, 1);
        noElim.matchDof(1,1, 1,2, 1);
        noElim.finalize();
        noElim.setShift(13);
        CHECK_EQUAL(15, (index_t)noElim.anyPreImages(0).size());

        gsDofMapper m = makeM2();
        m.setShift(13);

        const pairs want0 = {{0,1},{0,2},{1,1},{1,2},{1,3},{1,4},{0,3},
                             {-1,0},{-1,0},{-1,0},{-1,0},{-1,0},{-1,0},
                             {0,0},{-1,0}};
        const pairs want1 = {{-1,0},{-1,0},{-1,0},{-1,0},{-1,0},{-1,0},{-1,0},
                             {0,0},{0,1},{0,2},{1,3},{0,3},{1,1},
                             {-1,0},{1,4}};
        CHECK(m.anyPreImages(0) == want0);
        CHECK(m.anyPreImages(1) == want1);
    }

    TEST(DofMapperT6_firstIndex)
    {
        gsDofMapper m = makeM2();
        CHECK_EQUAL(0, m.firstIndex(0));
        CHECK_EQUAL(7, m.firstIndex(1));
    }

    TEST(DofMapperT7_shift0)
    {
        gsDofMapper m = makeM2();

        // Keys are global indices, values are local dofs of the given patch
        typedef std::map<index_t,index_t> imap;
        imap inv0 = m.inverseOnPatch(0);
        imap want0;
        want0[0] = 1; want0[1] = 2; want0[6] = 3; want0[7] = 0;
        want0[8] = 1; want0[9] = 2; want0[11] = 3; want0[13] = 0;
        CHECK(inv0 == want0);

        imap inv1 = m.inverseOnPatch(1);
        imap want1;
        want1[2] = 1; want1[3] = 2; want1[4] = 3; want1[5] = 4; want1[6] = 0;
        want1[10] = 3; want1[11] = 0; want1[14] = 4;
        CHECK_EQUAL(9, (index_t)inv1.size());
        // Key 12 is shared by local dofs 1 and 2; which one is reported is unspecified
        CHECK_EQUAL(1, (index_t)inv1.count(12));
        CHECK(inv1[12] == 1 || inv1[12] == 2);
        inv1.erase(12);
        CHECK(inv1 == want1);
    }

    // Eliminated local dofs, counted once per patch-local occurrence, in all components
    TEST(DofMapperT14_boundarySizeWithDuplicates)
    {
        CHECK_EQUAL(2, makeM2().boundarySizeWithDuplicates());

        // An eliminated dof shared by two patches counts twice
        gsVector<index_t> sz(2); sz << 4, 5;
        gsDofMapper m(sz, 2);
        m.matchDof(0,3, 1,0, 1);
        m.eliminateDof(3,0, 1);
        m.finalize();
        CHECK_EQUAL(1, m.boundarySize());
        CHECK_EQUAL(2, m.boundarySizeWithDuplicates());
    }

    // Tags are final dof numbers, which exist only after finalize()
    TEST(DofMapperT15_markTaggedNeedsFinalize)
    {
        gsVector<index_t> sz(2); sz << 4, 5;
        gsDofMapper m(sz, 2);
        CHECK_THROW(m.markTagged(1,0,0), std::runtime_error);
    }

    TEST(DofMapperT8_setIdentityReuse)
    {
        gsDofMapper m = makeM2();
        m.markTagged(1,0,0);
        m.setIdentity(3,20,2);

        CHECK_EQUAL(0, m.taggedSize());
        CHECK_EQUAL(2, (index_t)m.numComponents());
        CHECK_EQUAL(3, (index_t)m.numPatches());
        CHECK_EQUAL(0, (index_t)m.offset(0));
        CHECK_EQUAL(0, (index_t)m.offset(1));
        CHECK_EQUAL(0, (index_t)m.offset(2));
        CHECK_EQUAL(40, (index_t)m.mapSize());
        CHECK_EQUAL(20, (index_t)m.totalSize(0));
        CHECK_EQUAL(20, (index_t)m.totalSize(1));
        const vec zeros(20, 0);
        CHECK(toVec(m.asVector(0)) == zeros);
        CHECK(toVec(m.asVector(1)) == zeros);

        m.finalize();
        CHECK_EQUAL(40, m.size());
        CHECK_EQUAL(40, m.freeSize());
        const vec id0 = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19};
        const vec id1 = {20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39};
        CHECK(toVec(m.asVector(0)) == id0);
        CHECK(toVec(m.asVector(1)) == id1);
    }

    TEST(DofMapperT12_inverseAsVectorOtherComp)
    {
        // Two patches (2 and 3 dofs), two components, no matching or
        // elimination: component 0 owns global indices 0-4 and component 1
        // owns 5-9. Entries of inverseAsVector(c) outside component c's
        // index block have no preimage and are -1.
        gsVector<index_t> sz(2); sz << 2, 3;
        gsDofMapper m(sz, 2);
        m.finalize();

        CHECK_EQUAL(10, m.size());
        CHECK_EQUAL(0, m.index(0,0,0));
        CHECK_EQUAL(4, m.index(2,1,0));
        CHECK_EQUAL(5, m.index(0,0,1));
        CHECK_EQUAL(9, m.index(2,1,1));

        const vec want0 = {0,1,2,3,4,-1,-1,-1,-1,-1};
        const vec want1 = {-1,-1,-1,-1,-1,0,1,2,3,4};
        const gsVector<index_t> v0 = m.inverseAsVector(0);
        const gsVector<index_t> v1 = m.inverseAsVector(1);
        CHECK(toVec(v0) == want0);
        CHECK(toVec(v1) == want1);
        CHECK_EQUAL(5, (index_t)std::count(v0.data(), v0.data()+v0.size(), -1));
        CHECK_EQUAL(5, (index_t)std::count(v1.data(), v1.data()+v1.size(), -1));
    }

    TEST(DofMapperA1_inverseOnPatchShift13)
    {
        gsDofMapper m = makeM2();
        m.setShift(13);

        // Keys are shifted global indices, values are local dofs of the given patch
        typedef std::map<index_t,index_t> imap;
        imap inv0 = m.inverseOnPatch(0);
        imap want0;
        want0[13] = 1; want0[14] = 2; want0[19] = 3; want0[20] = 0;
        want0[21] = 1; want0[22] = 2; want0[24] = 3; want0[26] = 0;
        CHECK(inv0 == want0);

        imap inv1 = m.inverseOnPatch(1);
        imap want1;
        want1[15] = 1; want1[16] = 2; want1[17] = 3; want1[18] = 4; want1[19] = 0;
        want1[23] = 3; want1[24] = 0; want1[27] = 4;
        CHECK_EQUAL(9, (index_t)inv1.size());
        // Key 25 is shared by local dofs 1 and 2; which one is reported is unspecified
        CHECK_EQUAL(1, (index_t)inv1.count(25));
        CHECK(inv1.count(25) == 1 && (inv1[25] == 1 || inv1[25] == 2));
        inv1.erase(25);
        CHECK(inv1 == want1);
    }

#ifndef NDEBUG
    // GISMO_ASSERT throws std::logic_error only when NDEBUG is not defined
    TEST(DofMapperA6_isFreeIndexRange)
    {
        // Valid global indices of M2 at shift 13 are 13..27; free ones are 13..25
        gsDofMapper m = makeM2();
        m.setShift(13);
        CHECK(m.is_free_index(13));
        CHECK(m.is_free_index(25));
        CHECK(!m.is_free_index(26));
        CHECK(!m.is_free_index(27));
        CHECK_THROW(m.is_free_index(12), std::logic_error);
        CHECK_THROW(m.is_free_index(28), std::logic_error);
    }

    TEST(DofMapperA6_isBoundaryIndexRange)
    {
        // Valid global indices of M2 at shift 13 are 13..27; eliminated ones are 26 and 27
        gsDofMapper m = makeM2();
        m.setShift(13);
        CHECK(!m.is_boundary_index(13));
        CHECK(!m.is_boundary_index(25));
        CHECK(m.is_boundary_index(26));
        CHECK(m.is_boundary_index(27));
        CHECK_THROW(m.is_boundary_index(12), std::logic_error);
        CHECK_THROW(m.is_boundary_index(28), std::logic_error);
    }

    TEST(DofMapperA6_globalToBindexRange)
    {
        // Eliminated indices 26 and 27 are boundary dofs 0 and 1; 12 and 28 are out of range
        gsDofMapper m = makeM2();
        m.setShift(13);
        CHECK_EQUAL(0, m.global_to_bindex(26));
        CHECK_EQUAL(1, m.global_to_bindex(27));
        CHECK_THROW(m.global_to_bindex(12), std::logic_error);
        CHECK_THROW(m.global_to_bindex(28), std::logic_error);
    }

    TEST(DofMapperA6_componentOfRange)
    {
        // Component 0 owns free 13..19 and eliminated 26; component 1 owns free 20..25 and eliminated 27
        gsDofMapper m = makeM2();
        m.setShift(13);
        CHECK_EQUAL(0, m.componentOf(13));
        CHECK_EQUAL(1, m.componentOf(25));
        CHECK_EQUAL(0, m.componentOf(26));
        CHECK_EQUAL(1, m.componentOf(27));
        CHECK_THROW(m.componentOf(12), std::logic_error);
        CHECK_THROW(m.componentOf(28), std::logic_error);
    }

    TEST(DofMapperA6_emptyMapperRejectsAll)
    {
        // A finalized default-constructed mapper has size() == 0, so it has no valid global index
        gsDofMapper m;
        m.finalize();
        CHECK_EQUAL(0, m.size());
        CHECK_THROW(m.is_free_index(0), std::logic_error);
        CHECK_THROW(m.is_boundary_index(0), std::logic_error);
        CHECK_THROW(m.global_to_bindex(0), std::logic_error);
    }
#endif

    TEST(DofMapperR1_raggedOracle)
    {
        // Two patches, three components with different sizes per patch:
        //
        //   size(k,c)    c=0  c=1  c=2
        //   patch 0       2    5    1
        //   patch 1       3    1    6
        //   totalSize     5    6    7     mapSize 18
        //   offset(1,c)   2    5    1     offset(0,c) = 0
        //
        // Component c occupies the flat block that starts at 0, 5, 11.
        gsMatrix<index_t> M(2,3);
        M << 2,5,1,
             3,1,6;
        const index_t blockStart[3] = {0, 5, 11};
        gsDofMapper m(M);

        CHECK_EQUAL(2, (index_t)m.numPatches());
        CHECK_EQUAL(18, (index_t)m.mapSize());
        CHECK_EQUAL(5, (index_t)m.totalSize(0));
        CHECK_EQUAL(6, (index_t)m.totalSize(1));
        CHECK_EQUAL(7, (index_t)m.totalSize(2));
        for (index_t c = 0; c < 3; ++c)
        {
            CHECK_EQUAL(M(0,c), (index_t)m.patchSize(0,c));
            CHECK_EQUAL(M(1,c), (index_t)m.patchSize(1,c));
            CHECK_EQUAL(0, (index_t)m.offset(0,c));
            CHECK_EQUAL(M(0,c), (index_t)m.offset(1,c));
        }

        // Without conditions the numbering is the concatenation of the blocks
        gsDofMapper id(M);
        id.finalize();
        const vec want0 = {0,1,2,3,4};
        const vec want1 = {5,6,7,8,9,10};
        const vec want2 = {11,12,13,14,15,16,17};
        CHECK(toVec(id.asVector(0)) == want0);
        CHECK(toVec(id.asVector(1)) == want1);
        CHECK(toVec(id.asVector(2)) == want2);
        for (index_t c = 0; c < 3; ++c)
            for (index_t k = 0; k < 2; ++k)
                for (index_t i = 0; i < M(k,c); ++i)
                {
                    CHECK_EQUAL(blockStart[c] + (index_t)id.offset(k,c) + i, id.index(i,k,c));
                    CHECK_EQUAL(c, id.componentOf(id.index(i,k,c)));
                }

        // Eliminating dofs 0 and 4 of component 1 on patch 0 sends them behind
        // all free dofs: component 1 keeps 4 free dofs numbered 5..8 after
        // component 0, and component 2 follows with 9..15
        gsMatrix<index_t> b(2,1);
        b << 0, 4;
        m.markBoundary(0, b, 1);
        m.finalize();

        const vec bnd0 = {0,1,2,3,4};
        const vec bnd1 = {16,5,6,7,17,8};
        const vec bnd2 = {9,10,11,12,13,14,15};
        CHECK(toVec(m.asVector(0)) == bnd0);
        CHECK(toVec(m.asVector(1)) == bnd1);
        CHECK(toVec(m.asVector(2)) == bnd2);
        CHECK_EQUAL(18, m.size());
        CHECK_EQUAL(16, m.freeSize());
        CHECK_EQUAL(2, m.boundarySize());
        CHECK_EQUAL(5, m.freeSize(0));
        CHECK_EQUAL(4, m.freeSize(1));
        CHECK_EQUAL(7, m.freeSize(2));
        CHECK_EQUAL(6, m.size(1));
        CHECK_EQUAL(0, m.bindex(0,0,1));
        CHECK_EQUAL(1, m.bindex(4,0,1));

        vec all;
        for (index_t c = 0; c < 3; ++c)
        {
            const vec v = toVec(m.asVector(c));
            all.insert(all.end(), v.begin(), v.end());
        }
        vec sorted = all;
        std::sort(sorted.begin(), sorted.end());
        vec range(18);
        for (index_t n = 0; n < 18; ++n) range[n] = n;
        CHECK(sorted == range);

        for (index_t c = 0; c < 3; ++c)
            for (index_t k = 0; k < 2; ++k)
                for (index_t i = 0; i < M(k,c); ++i)
                    CHECK_EQUAL(c, m.componentOf(m.index(i,k,c)));

        for (index_t n = 0; n < 18; ++n)
            CHECK_EQUAL(all[n], m.mapIndex(n));
    }

    TEST(DofMapperR2_matrixCtorMatchesVectorCtor)
    {
        // Two patches of 3 and 2 dofs, three components. The matrix with the
        // size vector replicated per component describes the same mapper.
        gsVector<index_t> v(2);
        v << 3, 2;
        gsMatrix<index_t> M = v.replicate(1,3);
        gsDofMapper a(M), b(v, 3);

        gsMatrix<index_t> bd(2,1);
        bd << 0, 2;
        gsDofMapper * both[2] = {&a, &b};
        for (index_t j = 0; j < 2; ++j)
        {
            both[j]->markBoundary(0, bd, 1);
            both[j]->matchDof(0,1, 1,0, 2);
            both[j]->eliminateDof(1,1, 0);
            both[j]->finalize();
        }

        CHECK_EQUAL(14, a.size());
        CHECK_EQUAL(a.size(), b.size());
        CHECK_EQUAL(a.freeSize(), b.freeSize());
        CHECK_EQUAL(a.boundarySize(), b.boundarySize());
        CHECK_EQUAL(a.coupledSize(), b.coupledSize());
        CHECK_EQUAL((index_t)a.mapSize(), (index_t)b.mapSize());
        CHECK_EQUAL((index_t)a.numPatches(), (index_t)b.numPatches());
        for (index_t c = 0; c < 3; ++c)
        {
            CHECK(toVec(a.asVector(c)) == toVec(b.asVector(c)));
            for (index_t k = 0; k < 2; ++k)
            {
                CHECK_EQUAL((index_t)a.patchSize(k,c), (index_t)b.patchSize(k,c));
                CHECK_EQUAL((index_t)a.offset(k,c), (index_t)b.offset(k,c));
            }
        }
    }

    TEST(DofMapperR3_minusOneComponent)
    {
        // Component -1 acts on the same local dof in every component; on a
        // ragged mapper it is only allowed for patches of equal size in all
        // components. Patch 0 has 3 dofs in both components, patch 1 has 2 and 4.
        gsMatrix<index_t> R(2,2);
        R << 3,3,
             2,4;
        {
            gsDofMapper r(R);
            CHECK_THROW(r.eliminateDof(0, 1, -1), std::runtime_error);
        }
        {
            gsDofMapper r(R);
            CHECK_THROW(r.matchDof(0,2, 1,0, -1), std::runtime_error);
        }
        {
            // Component 0 holds patch 0 (dof 0 eliminated, dofs 1 and 2 coupled
            // into one) and patch 1 (dofs 0,1 free): free 0,1 and coupled 2,
            // eliminated 8. Component 1 follows: free 3..6, coupled 7,
            // eliminated 9.
            gsDofMapper r(R);
            r.eliminateDof(0, 0, -1);
            r.matchDof(0,1, 0,2, -1);
            r.finalize();
            const vec want0 = {8,2,2,0,1};
            const vec want1 = {9,7,7,3,4,5,6};
            CHECK(toVec(r.asVector(0)) == want0);
            CHECK(toVec(r.asVector(1)) == want1);
            CHECK_EQUAL(10, r.size());
            CHECK_EQUAL(8, r.freeSize());
            CHECK_EQUAL(2, r.boundarySize());
        }

        // Uniform mapper: component -1 reaches every component
        gsVector<index_t> v(2);
        v << 3, 2;
        gsDofMapper m(v, 2);
        m.eliminateDof(1,0, -1);
        m.matchDof(0,2, 1,0, -1);
        m.finalize();
        for (index_t c = 0; c < 2; ++c)
        {
            CHECK(m.is_boundary(1,0,c));
            CHECK_EQUAL(m.index(2,0,c), m.index(0,1,c));
        }
        CHECK_EQUAL(2, m.boundarySize());
        CHECK_EQUAL(6, m.freeSize());
        CHECK_EQUAL(8, m.size());
    }

    TEST(DofMapperR4_defaultCtor)
    {
        // Construct in a poisoned buffer so that a member the default
        // constructor leaves uninitialised reads back as 0xAB bytes
        alignas(gsDofMapper) unsigned char buf[sizeof(gsDofMapper)];
        std::memset(buf, 0xAB, sizeof buf);
        gsDofMapper * m = new (buf) gsDofMapper();
        CHECK_EQUAL(1, (index_t)m->numPatches());
        CHECK_EQUAL(0, (index_t)m->offset(0));
        CHECK_EQUAL(0, (index_t)m->patchSize(0));
        CHECK_EQUAL(0, m->firstIndex(0));
        CHECK_EQUAL(0, m->numComponents());
        m->~gsDofMapper();
    }

    // Equal columns give the mapper of gsDofMapper(column, nComp), also after
    // coupling and elimination
    TEST(matrix_equal_columns_as_vector)
    {
        gsVector<index_t> sizes(3);
        sizes << 3, 0, 4;
        for (index_t nc = 1; nc <= 3; nc += 2)
        {
            gsDofMapper a(sizes, nc), b(gsMatrix<index_t>(sizes.replicate(1, nc)));
            for (gsDofMapper * m : {&a, &b})
            {
                m->matchDof(0, 2, 2, 0, -1);
                m->matchDof(0, 1, 2, 1, nc - 1);
                m->matchDof(2, 1, 2, 3, 0);
                m->eliminateDof(0, 0, -1);
                m->markBoundary(2, idx({0, 3}), nc - 1);
                m->markCoupled(2, 2, 0);
                m->finalize();
            }
            CHECK_EQUAL(state(a), state(b));
        }
    }

    // Each component numbers its dofs like a scalar mapper of its column.
    // Component 1 has an empty patch, component 3 is left untouched.
    TEST(matrix_components_as_scalar_mappers)
    {
        gsMatrix<index_t> sizes(3, 4);
        sizes << 5, 4, 3, 2,
                 7, 0, 3, 2,
                 2, 6, 3, 0;
        // Array of mapper-operations to be applied to each component
        std::vector<compOps> ops;
        ops.push_back([](gsDofMapper & m, index_t c) {
            m.matchDofs(0, idx({1, 3}), 1, idx({0, 2}), c);
            m.eliminateDof(0, 0, c);
            m.matchDof(0, 0, 1, 1, c);
            m.markCoupled(1, 2, c);
            m.matchDof(1, 6, 2, 0, c); });
        ops.push_back([](gsDofMapper & m, index_t c) {
            m.matchDof(0, 3, 2, 0, c);
            m.eliminateDof(5, 2, c);
            m.matchDof(0, 2, 2, 5, c);
            m.markBoundary(2, idx({1, 2}), c); });
        ops.push_back([](gsDofMapper & m, index_t c) {
            for (index_t k = 0; k != 3; ++k)
                m.eliminateDof(1, k, c);
            m.matchDof(0, 2, 1, 0, c);
            m.matchDof(1, 2, 2, 0, c); });
        ops.push_back([](gsDofMapper &, index_t) { });
        const gsDofMapper m = checkAgainstScalars(sizes, ops);
        CHECK_EQUAL("off=0 5 12 ps=5 7 2;4 0 6;3 3 3;2 2 0", layout(m));
    }

    // Raviart-Thomas type pair on two patches of 4x2 elements: 28 and 30 dofs
    // per patch. Compared component-wise to scalar mappers.
    TEST(matrix_raviart_thomas_bases)
    {
        const rtFixture rt(4, 2, 2);
        const gsDofMapper m = checkAgainstScalars(rt.sizes, rt.ops);
        CHECK(m.coupledSize() > 0);
    }

    // The per-patch sizes decide the layout, not the component totals (58 each)
    TEST(matrix_patch_sizes_differ_with_equal_totals)
    {
        gsMatrix<index_t> sizes(2, 2);
        sizes << 28, 31,
                 30, 27;
        gsDofMapper m(sizes);
        m.finalize();
        CHECK_EQUAL("off=0 28 ps=28 30;31 27", layout(m));

        // Flat local dof 29 is on patch 1 for component 0, on patch 0 for component 1
        std::vector<dofPair> pre;
        m.preImage(29, pre);
        CHECK_EQUAL("(1,1)", join(pre));
        m.preImage(58 + 29, pre);
        CHECK_EQUAL("(0,29)", join(pre));
        CHECK(dofPair(1, 0) == m.anyPreImage(58 + 31));
        CHECK_EQUAL("29 -", onPatches(m, 58 + 29));
        CHECK_EQUAL("- 26", onPatches(m, 58 + 57));
    }

    // comp = -1 throws on unequal components and leaves the mapper unchanged;
    // debug builds reject a local dof beyond its patch
    TEST(matrix_unequal_components_reject_shared_calls)
    {
        gsMatrix<index_t> sizes(2, 2);
        sizes << 4, 2,
                 4, 6;
        gsDofMapper m(sizes);
        CHECK_THROW(m.matchDof(0, 3, 1, 0, -1), std::runtime_error);
        CHECK_THROW(m.eliminateDof(0, 0, -1), std::runtime_error);
    #ifndef NDEBUG
        CHECK_THROW(m.eliminateDof(6, 1, 1), std::logic_error);
        CHECK_THROW(m.matchDof(1, 6, 0, 0, 1), std::logic_error);
        // Within the component, beyond the patch
        CHECK_THROW(m.eliminateDof(3, 0, 1), std::logic_error);
        CHECK_THROW(m.matchDof(0, 2, 1, 0, 1), std::logic_error);
    #endif
        m.finalize();
        CHECK_EQUAL("off=0 4 ps=4 4;2 6 0,1,2,3|4,5,6,7;8,9|10,11,12,13,14,15", layout(m) + " " + indices(m));
    }

    TEST(matrix_invalid_sizes_throw)
    {
        const index_t big = std::numeric_limits<index_t>::max() - 1;
        gsMatrix<index_t> noRows(0, 2), noCols(2, 0), negative(2, 2), compTotal(2, 1), total(1, 2);
        negative << 3, 4, -1, 2;
        compTotal << big, 1;
        total << big, 1;
        CHECK_THROW(gsDofMapper m(noRows), std::runtime_error);
        CHECK_THROW(gsDofMapper m(noCols), std::runtime_error);
        CHECK_THROW(gsDofMapper m(negative), std::runtime_error);
        CHECK_THROW(gsDofMapper m(compTotal), std::runtime_error);
        CHECK_THROW(gsDofMapper m(total), std::runtime_error);
    }

    // The matrix constructor is explicit; Eigen expressions still select the
    // vector constructor (one component, one patch per entry)
    static_assert(!std::is_convertible<gsMatrix<index_t>, gsDofMapper>::value,
                  "gsDofMapper(const gsMatrix<index_t> &) must be explicit");
    TEST(matrix_constructor_overloads)
    {
        gsVector<index_t> v(2);
        v << 3, 4;
        gsMatrix<index_t> M(2, 2);
        M << 3, 5,
             4, 6;
        const gsDofMapper cst(gsVector<index_t>::Constant(2, 5)), plus(v.array() + 1), col(M.col(1)),
                          sum(v + v), row(v.transpose()), mat(M);
        for (const gsDofMapper * m : {&cst, &plus, &col, &sum, &row})
        {
            CHECK_EQUAL(1, m->numComponents());
            CHECK_EQUAL(2u, m->numPatches());
        }
        CHECK_EQUAL(2, mat.numComponents());
        CHECK_EQUAL(2u, mat.numPatches());
    }

    // setIdentity() on a moved-from mapper, which has no offsets left
    TEST(identity_after_move)
    {
        gsVector<index_t> sizes(2);
        sizes << 5, 7;
        gsDofMapper m(sizes);
        const gsDofMapper moved(std::move(m));
        m.setIdentity(2, 4, 1);
        m.finalize();
        CHECK_EQUAL("off=0 0 ps=0 4 |0,1,2,3", layout(m) + " " + indices(m));
    }

    // After setIdentity() the dofs are addressed through any patch, as callers
    // with mapped bases do
    TEST(identity_then_writes)
    {
        gsDofMapper m;
        m.setIdentity(3, 6, 2);
        m.markBoundary(1, idx({4}), 0);
        m.matchDof(0, 1, 2, 3, 1);
        m.finalize();
        CHECK_EQUAL("off=0 0 0 ps=0 0 6;0 0 6 ||0,1,2,3,10,4;||5,9,6,9,7,8", layout(m) + " " + indices(m));
    }

    // Patch-wise queries of a ragged, shifted mapper agree with the per-dof
    // queries; the find* functions report component 0
    TEST(matrix_patch_queries_with_shift)
    {
        gsMatrix<index_t> sizes(2, 2);
        sizes << 3, 5,
                 4, 2;
        gsDofMapper m(sizes);
        m.eliminateDof(0, 0, 0);
        m.matchDof(0, 2, 1, 0, 0);
        m.matchDof(0, 4, 1, 0, 1);
        m.eliminateDof(1, 1, 1);
        m.finalize();
        m.setShift(7);
        m.markTagged(1, 0, 0);
        m.markTagged(1, 1, 1);
        CHECK_EQUAL(2, m.taggedSize());
        CHECK(m.is_tagged(1, 1, 1));
    #ifndef NDEBUG
        CHECK_THROW(m.patchSize(0, 2), std::logic_error);
    #endif

        for (index_t k = 0; k != 2; ++k)
        {
            const std::map<index_t,index_t> inv = m.inverseOnPatch(k);
            size_t n = 0;
            for (index_t c = 0; c != 2; ++c)
                for (index_t i = 0; i != static_cast<index_t>(m.patchSize(k, c)); ++i, ++n)
                {
                    const index_t gl = m.index(i, k, c);
                    index_t loc = -1;
                    CHECK(m.indexOnPatch(gl, k, loc) && loc == i);
                    CHECK_EQUAL(c, m.componentOf(gl));
                    CHECK(inv.count(gl) && inv.at(gl) == i);
                }
            CHECK_EQUAL(n, inv.size());

            vec fr, bd, cp, fu, tg;
            for (index_t i = 0; i != static_cast<index_t>(m.patchSize(k)); ++i)
            {
                (m.is_free(i, k) ? fr : bd).push_back(i);
                if (m.is_coupled(i, k)) cp.push_back(i);
                else if (m.is_free(i, k)) fu.push_back(i);
                if (m.is_tagged(i, k)) tg.push_back(i);
            }
            CHECK(fr == toVec(m.findFree(k)));
            CHECK(bd == toVec(m.findBoundary(k)));
            CHECK(cp == toVec(m.findCoupled(k)));
            CHECK(fu == toVec(m.findFreeUncoupled(k)));
            CHECK(tg == toVec(m.findTagged(k)));
        }
        // Patch 0, component 0: dof 0 eliminated, dof 2 coupled with patch 1, dof 1 tagged
        CHECK(vec({1, 2}) == toVec(m.findFree(0)));
        CHECK(vec({0}) == toVec(m.findBoundary(0)));
        CHECK(vec({2}) == toVec(m.findCoupled(0, 1)));
        CHECK(vec({1}) == toVec(m.findFreeUncoupled(0)));
        CHECK(vec({1}) == toVec(m.findTagged(0)));
    }
}
