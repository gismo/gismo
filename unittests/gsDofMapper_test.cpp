/** @file gsDofMapper_test.cpp

    @brief Tests for gsDofMapper

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst
*/

#include "gismo_unittest.h"

#include <cstring>
#include <map>
#include <new>

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
  index; stored values (tags, anyPreImages positions, inverseOnPatch keys)
  and patch-local indices do not shift.

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
        const vec tagged = {6,11,12};
        CHECK(t.getTagged() == tagged);
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
        const vec tagM = {12};
        CHECK(m.is_tagged(1,1,1));
        CHECK(m.is_tagged(2,1,1));
        CHECK(m.is_tagged_index(13 + 12));
        CHECK(!m.is_tagged(3,1,1));
        CHECK_EQUAL(0, m.tindex(1,1,1));
        CHECK(m.getTagged() == tagM);

        // Shift first, then tag
        gsDofMapper f = makeM2();
        f.setShift(13);
        f.markTagged(1,1,1);
        const vec tagF = {12};
        CHECK(f.getTagged() == tagF);
        CHECK(f.is_tagged(1,1,1));

        // Shift first, then tag all coupled dofs
        gsDofMapper g = makeM2();
        g.setShift(13);
        g.markCoupledAsTagged();
        const vec tagG = {6,11,12};
        CHECK(g.getTagged() == tagG);
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

    TEST(DofMapperT7_shift13)
    {
        gsDofMapper m = makeM2();
        m.setShift(13);

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

    TEST(DofMapperT11_defaultCtorBoundaryShift)
    {
        // Construct in a poisoned buffer so that a member the default
        // constructor leaves uninitialised reads back as 0xAB bytes
        alignas(gsDofMapper) unsigned char buf[sizeof(gsDofMapper)];
        std::memset(buf, 0xAB, sizeof buf);
        gsDofMapper * m = new (buf) gsDofMapper();
        m->finalize();
        CHECK_EQUAL(0, m->global_to_bindex(0));
        m->~gsDofMapper();
    }
}
