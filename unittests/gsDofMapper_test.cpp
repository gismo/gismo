/** @file gsDofMapper_test.cpp

    @brief Tests for gsDofMapper

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst
*/

#include "gismo_unittest.h"

#include <algorithm>
#include <cstring>
#include <iostream>
#include <map>
#include <new>
#include <sstream>

SUITE(gsDofMapper_test)
{
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

    /*
      Mutes std::cout and std::cerr while it is alive; UnitTest++ reports
      failures through fprintf(stderr), so a failing check inside the block stays
      visible. GISMO_ASSERT and GISMO_ENSURE
      print their message to std::cerr before they throw, and
      UnitTest::deactivate_output() only mutes std::cout, so both are needed
      to keep the log of a passing run free of the expected messages.
    */
    struct MuteOutput
    {
        MuteOutput() : m_cerr(std::cerr.rdbuf(m_sink.rdbuf()))
        { UnitTest::deactivate_output(); }
        ~MuteOutput()
        { UnitTest::reactivate_output(); std::cerr.rdbuf(m_cerr); }
    private:
        std::ostringstream m_sink;
        std::streambuf * m_cerr;
    };
    } // namespace

    // Verifies the global numbering of the M2 fixture at shift 0 and every
    // query that reads it: dense index tables, sizes, component of an index,
    // free/boundary lists, preimages and patch-local lookup. All expected
    // values are read off the M2 table above.
    TEST(numbering_shift0)
    {
        gsDofMapper m = makeM2();

        // The two rows per component of the M2 table, and the size counters stated below it
        const vec as0 = {13,0,1,6, 6,2,3,4,5};
        const vec as1 = {7,8,9,11, 11,12,12,10,14};
        CHECK(toVec(m.asVector(0)) == as0);
        CHECK(toVec(m.asVector(1)) == as1);
        CHECK_EQUAL(13, m.freeSize());
        CHECK_EQUAL(15, m.size());
        CHECK_EQUAL(2, m.boundarySize());

        // Component owning each global index: component 0 owns 0-6 and 13, component 1 owns 7-12 and 14
        const index_t comp[15] = {0,0,0,0,0,0,0, 1,1,1,1,1,1, 0,1};
        bool componentsOk = true;
        for (index_t g = 0; g < 15; ++g)
        {
            const index_t c = m.componentOf(g);
            CHECK_EQUAL(comp[g], c);
            componentsOk = componentsOk && (c == comp[g]);
        }

        // Free and eliminated patch-local dofs of component 0: patch 0 has dof 0 eliminated, patch 1 has none
        const vec free0 = {1,2,3};
        const vec free1 = {0,1,2,3,4};
        const vec bnd0 = {0};
        CHECK(toVec(m.findFree(0)) == free0);
        CHECK(toVec(m.findFree(1)) == free1);
        CHECK(toVec(m.findBoundary(0)) == bnd0);
        CHECK_EQUAL(0, (index_t)m.findBoundary(1).size());

        // preImage(), anyPreImage() and indexOnPatch() search the component reported by componentOf()
        if (!componentsOk) return;

        // Inverse of the M2 table: for each global index the (patch, local dof) pairs mapped to it
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

        // anyPreImage() returns the first pair of the preImage() list above
        const std::pair<index_t,index_t> any0(0,1), any6(0,3), any12(1,1), any13(0,0), any14(1,4);
        CHECK(m.anyPreImage(0) == any0);
        CHECK(m.anyPreImage(6) == any6);
        CHECK(m.anyPreImage(12) == any12);
        CHECK(m.anyPreImage(13) == any13);
        CHECK(m.anyPreImage(14) == any14);

        // indexOnPatch(g, k, l) is true exactly when the preImage() list of g contains a pair on patch k
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

    // Same queries as numbering_shift0 after setShift(13): global indices
    // move by 13, patch-local indices and the sizes stay. Expected values are
    // the M2 table with 13 added to every global index.
    TEST(numbering_shift13)
    {
        gsDofMapper m = makeM2();
        m.setShift(13);

        // The M2 table rows plus 13; sizes are shift independent
        const vec as0 = {26,13,14,19, 19,15,16,17,18};
        const vec as1 = {20,21,22,24, 24,25,25,23,27};
        CHECK(toVec(m.asVector(0)) == as0);
        CHECK(toVec(m.asVector(1)) == as1);
        CHECK_EQUAL(13, m.freeSize());
        CHECK_EQUAL(15, m.size());
        CHECK_EQUAL(2, m.boundarySize());

        // Component owning each shifted index 13 + g: the same as for g at shift 0
        const index_t comp[15] = {0,0,0,0,0,0,0, 1,1,1,1,1,1, 0,1};
        bool componentsOk = true;
        for (index_t g = 0; g < 15; ++g)
        {
            const index_t c = m.componentOf(13 + g);
            CHECK_EQUAL(comp[g], c);
            componentsOk = componentsOk && (c == comp[g]);
        }

        // Patch-local lists do not shift
        const vec free0 = {1,2,3};
        const vec free1 = {0,1,2,3,4};
        const vec bnd0 = {0};
        CHECK(toVec(m.findFree(0)) == free0);
        CHECK(toVec(m.findFree(1)) == free1);
        CHECK(toVec(m.findBoundary(0)) == bnd0);
        CHECK_EQUAL(0, (index_t)m.findBoundary(1).size());

        // preImage(), anyPreImage() and indexOnPatch() search the component reported by componentOf()
        if (!componentsOk) return;

        // Preimages of the shifted index 13 + g are those of g at shift 0
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

        // Representative preimages as at shift 0
        const std::pair<index_t,index_t> any0(0,1), any6(0,3), any12(1,1), any13(0,0), any14(1,4);
        CHECK(m.anyPreImage(13 + 0) == any0);
        CHECK(m.anyPreImage(13 + 6) == any6);
        CHECK(m.anyPreImage(13 + 12) == any12);
        CHECK(m.anyPreImage(13 + 13) == any13);
        CHECK(m.anyPreImage(13 + 14) == any14);

        // Patch-local lookup as at shift 0, with shifted arguments
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

    // Verifies the coupled-dof queries of M2 at shift 0: the coupled and
    // free-uncoupled lists per patch, the per-index and per-dof coupled flags,
    // and that markCoupledAsTagged() tags exactly the coupled dofs. Expected
    // values follow from the coupled band in the M2 table above ({6} for
    // component 0, {11,12} for component 1).
    TEST(coupled_dofs_shift0)
    {
        gsDofMapper m = makeM2();

        // Lists refer to component 0 only: patch 0 couples local dof 3 (index 6), patch 1 local dof 0
        // (index 6). Free uncoupled are the free local dofs minus the coupled one.
        const vec cpl0 = {3};
        const vec cpl1 = {0};
        const vec fu0 = {1,2};
        const vec fu1 = {1,2,3,4};
        CHECK(toVec(m.findCoupled(0)) == cpl0);
        CHECK(toVec(m.findCoupled(1)) == cpl1);
        CHECK(toVec(m.findFreeUncoupled(0)) == fu0);
        CHECK(toVec(m.findFreeUncoupled(1)) == fu1);

        // findCoupled(k, j) keeps the coupled dofs of patch k (component 0) whose global index also occurs
        // on patch j. Index 6 is the only coupled one and sits at local dof 3 of patch 0 and local
        // dof 0 of patch 1, so (0,1) gives {3} (local indices of patch k) and (1,0) gives {0}.
        const vec cpl01 = {3};
        const vec cpl10 = {0};
        CHECK(toVec(m.findCoupled(0, 1)) == cpl01);
        CHECK(toVec(m.findCoupled(1, 0)) == cpl10);

        // Coupled band of both components by global index: 6 for component 0, 11 and 12 for component 1
        const bool coupled[15] = {0,0,0,0,0,0,1, 0,0,0,0,1,1, 0,0};
        for (index_t g = 0; g < 15; ++g)
            CHECK_EQUAL(coupled[g], m.is_coupled_index(g));

        // Per-dof flags: the dofs matched in makeM2() are coupled, an unmatched or eliminated one is not
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

        // Tagging the coupled dofs tags exactly the coupled band, so the tagged flags equal the coupled flags
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

    // The coupled-dof queries of coupled_dofs_shift0 after setShift(13). The
    // patch-local lists are shift independent; the per-index flags are those
    // of M2 at index 13 + g.
    TEST(coupled_dofs_shift13)
    {
        gsDofMapper m = makeM2();
        m.setShift(13);

        // Patch-local lists as at shift 0
        const vec cpl0 = {3};
        const vec cpl1 = {0};
        const vec fu0 = {1,2};
        const vec fu1 = {1,2,3,4};
        CHECK(toVec(m.findCoupled(0)) == cpl0);
        CHECK(toVec(m.findCoupled(1)) == cpl1);
        CHECK(toVec(m.findFreeUncoupled(0)) == fu0);
        CHECK(toVec(m.findFreeUncoupled(1)) == fu1);

        // Derivation as at shift 0: patch-local results do not depend on the shift
        const vec cpl01 = {3};
        const vec cpl10 = {0};
        CHECK(toVec(m.findCoupled(0, 1)) == cpl01);
        CHECK(toVec(m.findCoupled(1, 0)) == cpl10);

        // Coupled flag of the shifted index 13 + g equals that of g at shift 0
        const bool coupled[15] = {0,0,0,0,0,0,1, 0,0,0,0,1,1, 0,0};
        for (index_t g = 0; g < 15; ++g)
            CHECK_EQUAL(coupled[g], m.is_coupled_index(13 + g));

        // Per-dof flags are patch-local and as at shift 0
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

    // Verifies that a one-component mapper with a single coupled dof reports
    // it as coupled and that markCoupledAsTagged() tags only that dof. The
    // expectations follow from the setup: only dof 2 of the one patch is marked.
    TEST(coupled_dofs_one_component)
    {
        // One patch with 3 dofs, dof 2 coupled
        gsVector<index_t> s1(1); s1 << 3;
        gsDofMapper m(s1);
        m.markCoupled(2,0);
        m.finalize();
        CHECK(m.is_coupled(2,0,0));

        // Tagging the coupled dofs tags dof 2 and no other
        m.markCoupledAsTagged();
        CHECK(m.is_tagged(2,0,0));
        CHECK(!m.is_tagged(0,0,0));
        CHECK(!m.is_tagged(1,0,0));
        CHECK_EQUAL(1, m.taggedSize());
    }

    // Verifies that a tag, which is stored as an unshifted value, follows the
    // shift whether it is set before or after setShift(13). Expected indices
    // are the M2 table values plus 13: dof 1 of patch 1 in component 1 sits at
    // index 12 (so 25), and the coupled band of M2 is {6,11,12} (so 19, 24, 25).
    TEST(tag_follows_shift)
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

    // Verifies findTagged() at shift 0: it lists the patch-local dofs of
    // component 0 whose index is tagged. The expected lists are read off the
    // M2 table above.
    TEST(find_tagged_shift0)
    {
        // Tag local dof 1 of patch 0 (index 0), local dof 2 of patch 1 (index 3) and the coupled band
        gsDofMapper m = makeM2();
        m.markTagged(1,0,0);
        m.markTagged(2,1,0);
        m.markCoupledAsTagged();

        // Stored tags are {0,3,6,11,12}. In component 0 (indices 0-6) these are 0, 3 and 6: patch 0 has
        // index 0 at local dof 1 and index 6 at local dof 3, patch 1 has index 3 at local dof 2 and
        // index 6 at local dof 0.
        const vec t0 = {1,3};
        const vec t1 = {0,2};
        CHECK(toVec(m.findTagged(0)) == t0);
        CHECK(toVec(m.findTagged(1)) == t1);
    }

    // Same as find_tagged_shift0 with the shift applied after tagging: tags
    // are stored unshifted, so the patch-local lists do not change.
    TEST(find_tagged_shift13)
    {
        // Same tags as at shift 0, then shift
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

    // Verifies anyPreImages(c): one (patch, local dof) per global index, with
    // {-1,0} for indices owned by another component, at shift 0. Expected
    // pairs are the first preimage of each index in the M2 table.
    TEST(any_preimages_shift0)
    {
        // Same matches as M2 but nothing eliminated: 8 + 7 dofs
        gsVector<index_t> sz(2); sz << 4, 5;
        gsDofMapper noElim(sz, 2);
        noElim.matchDof(0,3, 1,0, 0);
        noElim.matchDof(0,3, 1,0, 1);
        noElim.matchDof(1,1, 1,2, 1);
        noElim.finalize();
        CHECK_EQUAL(15, (index_t)noElim.anyPreImages(0).size());

        // Component 0 fills positions 0-6 and 13; component 1 fills 7-12 and 14; the rest is {-1,0}
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

    // Same as any_preimages_shift0 at shift 13: the dense vectors are still
    // indexed by gl - firstIndex(0), so the contents are those of shift 0.
    TEST(any_preimages_shift13)
    {
        // Mapper without eliminated dofs, shifted
        gsVector<index_t> sz(2); sz << 4, 5;
        gsDofMapper noElim(sz, 2);
        noElim.matchDof(0,3, 1,0, 0);
        noElim.matchDof(0,3, 1,0, 1);
        noElim.matchDof(1,1, 1,2, 1);
        noElim.finalize();
        noElim.setShift(13);
        CHECK_EQUAL(15, (index_t)noElim.anyPreImages(0).size());

        // Shifted M2: position gl - firstIndex(0) holds the preimage of gl, as in the shift-0 table
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

    // Verifies firstIndex(c): the first global index of component c. The
    // M2 table gives 0 for component 0 and 7 for component 1.
    TEST(first_index)
    {
        // Component 0 owns indices 0-6, so component 1 starts at 7
        gsDofMapper m = makeM2();
        CHECK_EQUAL(0, m.firstIndex(0));
        CHECK_EQUAL(7, m.firstIndex(1));
    }

    // Verifies inverseOnPatch(k) at shift 0: it maps each global index used by
    // patch k to a local dof of patch k, over all components. Expected maps are
    // the M2 table restricted to one patch.
    TEST(inverse_on_patch_shift0)
    {
        gsDofMapper m = makeM2();

        // Patch 0: the table rows of patch 0 for both components, inverted. Keys are global
        // indices, values are local dofs of the given patch
        typedef std::map<index_t,index_t> imap;
        imap inv0 = m.inverseOnPatch(0);
        imap want0;
        want0[0] = 1; want0[1] = 2; want0[6] = 3; want0[7] = 0;
        want0[8] = 1; want0[9] = 2; want0[11] = 3; want0[13] = 0;
        CHECK(inv0 == want0);

        // Patch 1: nine distinct indices because local dofs 1 and 2 of component 1 share index 12
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
    TEST(boundary_size_with_duplicates)
    {
        // M2 eliminates one dof in each component, none shared between patches: 2 occurrences
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
    TEST(mark_tagged_needs_finalize)
    {
        // Not finalized: markTagged() has no final numbering to store, which is a GISMO_ENSURE
        // failure in every build
        gsVector<index_t> sz(2); sz << 4, 5;
        gsDofMapper m(sz, 2);
        {
            MuteOutput mute;
            CHECK_THROW(m.markTagged(1,0,0), std::runtime_error);
        }
    }

    // Verifies that setIdentity() on a mapper that was already in use resets
    // it: tags are cleared, the sizes follow the new arguments and finalize()
    // numbers the dofs consecutively per component. Expected values follow
    // from setIdentity(3,20,2): 2 components of 20 dofs each, 3 patches.
    TEST(set_identity_reuse)
    {
        // Reuse a finalized mapper that holds a tag
        gsDofMapper m = makeM2();
        m.markTagged(1,0,0);
        m.setIdentity(3,20,2);

        // Tags are cleared; 2 components x 20 dofs = 40 map entries, 20 per component
        CHECK_EQUAL(0, m.taggedSize());
        CHECK_EQUAL(2, (index_t)m.numComponents());
        CHECK_EQUAL(3, (index_t)m.numPatches());
        CHECK_EQUAL(40, (index_t)m.mapSize());
        CHECK_EQUAL(20, (index_t)m.totalSize(0));
        CHECK_EQUAL(20, (index_t)m.totalSize(1));
        const vec zeros(20, 0);
        CHECK(toVec(m.asVector(0)) == zeros);
        CHECK(toVec(m.asVector(1)) == zeros);

        // After finalize() the identity numbering is 0-19 for component 0 and 20-39 for component 1
        m.finalize();
        CHECK_EQUAL(40, m.size());
        CHECK_EQUAL(40, m.freeSize());
        const vec id0 = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19};
        const vec id1 = {20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39};
        CHECK(toVec(m.asVector(0)) == id0);
        CHECK(toVec(m.asVector(1)) == id1);
    }

    // Verifies that the default constructor initialises the boundary shift.
    // The only CHECK sits behind #ifdef NDEBUG, so this test verifies the
    // initialisation only in builds without assertions, such as CI Release
    // builds. With assertions live the default-constructed boundary shift
    // cannot be observed through the public API: initPatchDofs() and
    // setIdentity() both reset it, and global_to_bindex(0) on the empty range
    // asserts before it reads it.
    TEST(default_ctor_boundary_shift)
    {
        // Construct in a poisoned buffer so that a member the default
        // constructor leaves uninitialised reads back as 0xAB bytes
        alignas(gsDofMapper) unsigned char buf[sizeof(gsDofMapper)];
        std::memset(buf, 0xAB, sizeof buf);
        gsDofMapper * m = new (buf) gsDofMapper();
        m->finalize();
#ifdef NDEBUG
        // Index 0 is outside the empty range [0,0) and asserts when NDEBUG is not defined;
        // with a zero boundary shift it maps to boundary dof 0
        CHECK_EQUAL(0, m->global_to_bindex(0));
#endif
        m->~gsDofMapper();
    }

    // Verifies inverseAsVector(c) for a mapper without matching or
    // elimination: each component owns a contiguous index block and the
    // entries outside it are -1. Expected values follow from the layout
    // described in the setup comment below.
    TEST(inverse_as_vector_other_component)
    {
        // Two patches (2 and 3 dofs), two components, no matching or
        // elimination: component 0 owns global indices 0-4 and component 1
        // owns 5-9. Entries of inverseAsVector(c) outside component c's
        // index block have no preimage and are -1.
        gsVector<index_t> sz(2); sz << 2, 3;
        gsDofMapper m(sz, 2);
        m.finalize();

        // Indices are consecutive over patches within a component: patch 1 starts after the 2 dofs of patch 0
        CHECK_EQUAL(10, m.size());
        CHECK_EQUAL(0, m.index(0,0,0));
        CHECK_EQUAL(4, m.index(2,1,0));
        CHECK_EQUAL(5, m.index(0,0,1));
        CHECK_EQUAL(9, m.index(2,1,1));

        // Position g of inverseAsVector(c) holds the position of g in the component-c block, or -1
        const vec want0 = {0,1,2,3,4,-1,-1,-1,-1,-1};
        const vec want1 = {-1,-1,-1,-1,-1,0,1,2,3,4};
        const gsVector<index_t> v0 = m.inverseAsVector(0);
        const gsVector<index_t> v1 = m.inverseAsVector(1);
        CHECK(toVec(v0) == want0);
        CHECK(toVec(v1) == want1);
        CHECK_EQUAL(5, (index_t)std::count(v0.data(), v0.data()+v0.size(), -1));
        CHECK_EQUAL(5, (index_t)std::count(v1.data(), v1.data()+v1.size(), -1));
    }

    // Verifies permuteFreeDofs() on a single-component mapper at shift 13: the
    // free indices are reversed (p maps k to 3-k) and the tag follows its dof.
    // Expected indices are 13 + p(k) for local dof k = 0..3.
    TEST(permute_single_component)
    {
        // One patch of 4 free dofs, local dof 1 tagged, shifted by 13
        gsVector<index_t> sz(1); sz << 4;
        gsDofMapper m(sz);
        m.finalize();
        m.markTagged(1,0);
        m.setShift(13);

        // Reverse the numbering of the free dofs
        gsVector<index_t> p(4); p << 3, 2, 1, 0;
        m.permuteFreeDofs(p, 0);

        // Local dof k moves to free index p(k), reported with the shift 13 added; the tag stays on local dof 1
        CHECK_EQUAL(16, m.index(0,0));
        CHECK_EQUAL(15, m.index(1,0));
        CHECK_EQUAL(14, m.index(2,0));
        CHECK_EQUAL(13, m.index(3,0));
        CHECK(m.is_tagged(1,0));
        CHECK_EQUAL(1, m.taggedSize());
    }

    // Verifies that permuteFreeDofs() keeps a tagged eliminated dof tagged
    // while it renumbers the free dofs. Eliminated indices lie past the free
    // block and are not permuted, so the tag of local dof 3 stays at index 3
    // and the tag of local dof 0 moves with its dof from index 0 to p(0) = 2.
    TEST(permute_keeps_eliminated_tagged)
    {
        // One patch of 4 dofs, local dof 3 eliminated: free indices 0-2, eliminated index 3.
        // Tag the free local dof 0 and the eliminated local dof 3.
        gsVector<index_t> sz(1); sz << 4;
        gsDofMapper m(sz);
        m.eliminateDof(3,0);
        m.finalize();
        m.markTagged(0,0);
        m.markTagged(3,0);
        CHECK_EQUAL(3, m.freeSize());
        CHECK_EQUAL(2, m.taggedSize());

        // Reverse the three free indices
        gsVector<index_t> p(3); p << 2, 1, 0;
        m.permuteFreeDofs(p, 0);

        // Free dofs are renumbered, the eliminated one is not
        CHECK_EQUAL(2, m.index(0,0));
        CHECK_EQUAL(0, m.index(2,0));
        CHECK_EQUAL(3, m.index(3,0));

        // Both tags survive: the free one at its new index 2, the eliminated one at index 3
        CHECK_EQUAL(2, m.taggedSize());
        CHECK(m.is_tagged(0,0));
        CHECK(m.is_tagged_index(2));
        CHECK(!m.is_tagged_index(0));
        CHECK(m.is_tagged(3,0));
        CHECK(m.is_tagged_index(3));
    }

    // Verifies that permuteFreeDofs() rejects a mapper with more than one
    // component. This is a GISMO_ENSURE and therefore throws in every build.
    TEST(permute_rejects_multi_component)
    {
        // Two components over one patch of 3 dofs
        gsVector<index_t> sz(1); sz << 3;
        gsDofMapper m(sz, 2);
        m.finalize();

        gsVector<index_t> p(3); p << 2, 1, 0;
        {
            MuteOutput mute;
            CHECK_THROW(m.permuteFreeDofs(p, 0), std::exception);
        }
    }

    // Verifies inverseOnPatch(k) of M2 at shift 13: the keys are the shift-0
    // keys of inverse_on_patch_shift0 plus 13, the values are unchanged local dofs.
    TEST(inverse_on_patch_shift13)
    {
        gsDofMapper m = makeM2();
        m.setShift(13);

        // Patch 0: the T7 map with 13 added to every key
        // Keys are shifted global indices, values are local dofs of the given patch
        typedef std::map<index_t,index_t> imap;
        imap inv0 = m.inverseOnPatch(0);
        imap want0;
        want0[13] = 1; want0[14] = 2; want0[19] = 3; want0[20] = 0;
        want0[21] = 1; want0[22] = 2; want0[24] = 3; want0[26] = 0;
        CHECK(inv0 == want0);

        // Patch 1: the T7 map with 13 added to every key; shared key 12 becomes 25
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

    // The range checks below are GISMO_ASSERTs: they throw std::logic_error
    // only when NDEBUG is not defined. The in-range checks run in every build.

    // Verifies is_free_index() at shift 13: valid global indices of M2 are
    // 13..27 and the free ones 13..25; 12 and 28 are outside the range.
    TEST(is_free_index_range)
    {
        // Valid global indices of M2 at shift 13 are 13..27; free ones are 13..25
        gsDofMapper m = makeM2();
        m.setShift(13);
        CHECK(m.is_free_index(13));
        CHECK(m.is_free_index(25));
        CHECK(!m.is_free_index(26));
        CHECK(!m.is_free_index(27));
        {
            MuteOutput mute;
            CHECK_THROW_IN_DEBUG(m.is_free_index(12), std::logic_error);
            CHECK_THROW_IN_DEBUG(m.is_free_index(28), std::logic_error);
        }
    }

    // Verifies is_boundary_index() at shift 13: the eliminated indices of M2
    // are 26 and 27; 12 and 28 are outside the range.
    TEST(is_boundary_index_range)
    {
        // Valid global indices of M2 at shift 13 are 13..27; eliminated ones are 26 and 27
        gsDofMapper m = makeM2();
        m.setShift(13);
        CHECK(!m.is_boundary_index(13));
        CHECK(!m.is_boundary_index(25));
        CHECK(m.is_boundary_index(26));
        CHECK(m.is_boundary_index(27));
        {
            MuteOutput mute;
            CHECK_THROW_IN_DEBUG(m.is_boundary_index(12), std::logic_error);
            CHECK_THROW_IN_DEBUG(m.is_boundary_index(28), std::logic_error);
        }
    }

    // Verifies global_to_bindex() at shift 13: the two eliminated indices of
    // M2 are boundary dofs 0 and 1; 12 and 28 are outside the range.
    TEST(global_to_bindex_range)
    {
        // Eliminated indices 26 and 27 are boundary dofs 0 and 1; 12 and 28 are out of range
        gsDofMapper m = makeM2();
        m.setShift(13);
        CHECK_EQUAL(0, m.global_to_bindex(26));
        CHECK_EQUAL(1, m.global_to_bindex(27));
        {
            MuteOutput mute;
            CHECK_THROW_IN_DEBUG(m.global_to_bindex(12), std::logic_error);
            CHECK_THROW_IN_DEBUG(m.global_to_bindex(28), std::logic_error);
        }
    }

    // Verifies componentOf() at shift 13: component ownership of the shifted
    // indices follows the M2 table; 12 and 28 are outside the range.
    TEST(component_of_range)
    {
        // Component 0 owns free 13..19 and eliminated 26; component 1 owns free 20..25 and eliminated 27
        gsDofMapper m = makeM2();
        m.setShift(13);
        CHECK_EQUAL(0, m.componentOf(13));
        CHECK_EQUAL(1, m.componentOf(25));
        CHECK_EQUAL(0, m.componentOf(26));
        CHECK_EQUAL(1, m.componentOf(27));
        {
            MuteOutput mute;
            CHECK_THROW_IN_DEBUG(m.componentOf(12), std::logic_error);
            CHECK_THROW_IN_DEBUG(m.componentOf(28), std::logic_error);
        }
    }

    // Verifies that a finalized default-constructed mapper has no valid
    // global index: its size is 0, so every index query is out of range.
    TEST(empty_mapper_rejects_all)
    {
        // A finalized default-constructed mapper has size() == 0, so it has no valid global index
        gsDofMapper m;
        m.finalize();
        CHECK_EQUAL(0, m.size());
        {
            MuteOutput mute;
            CHECK_THROW_IN_DEBUG(m.is_free_index(0), std::logic_error);
            CHECK_THROW_IN_DEBUG(m.is_boundary_index(0), std::logic_error);
            CHECK_THROW_IN_DEBUG(m.global_to_bindex(0), std::logic_error);
        }
    }

    // Verifies that permuteFreeDofs() rejects component 1 of a single-component
    // mapper. This is a GISMO_ENSURE and therefore throws in every build.
    TEST(permute_rejects_component_one)
    {
        // One component, so only component 0 exists
        gsVector<index_t> sz(1); sz << 4;
        gsDofMapper m(sz);
        m.finalize();

        gsVector<index_t> p(4); p << 3, 2, 1, 0;
        {
            MuteOutput mute;
            CHECK_THROW(m.permuteFreeDofs(p, 1), std::runtime_error);
        }
    }
}
