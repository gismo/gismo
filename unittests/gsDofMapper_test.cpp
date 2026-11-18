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
#include <iostream>
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

SUITE(gsDofMapper_test)
{
    // Helpers for the tests on ragged mappers (per-patch, per-component sizes)
    namespace
    {
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
        const index_t nP = static_cast<index_t>(m.numPatches());
        std::ostringstream os;
        os << "off=";
        for (index_t k = 0; k != nP; ++k)
            os << (k ? " " : "") << m.offset(k);
        os << " ps=";
        for (index_t c = 0; c != m.numComponents(); ++c)
            for (index_t k = 0; k != nP; ++k)
                os << (k ? " " : (c ? ";" : "")) << m.patchSize(k, c);
        return os.str();
    }

    // index(i,k,c) of every local dof, as "c0k0|c0k1;c1k0|c1k1"
    std::string indices(const gsDofMapper & m)
    {
        const index_t nP = static_cast<index_t>(m.numPatches());
        std::ostringstream os;
        for (index_t c = 0; c != m.numComponents(); ++c)
            for (index_t k = 0; k != nP; ++k)
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
        const index_t nP = static_cast<index_t>(m.numPatches());
        std::ostringstream os;
        m.print(os);
        os << layout(m) << " map=" << m.mapSize() << "\n";
        for (index_t c = 0; c != m.numComponents(); ++c)
        {
            os << "c" << c << ":" << m.size(c) << "," << m.freeSize(c) << "," << m.firstIndex(c);
            for (index_t k = 0; k != nP; ++k)
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
        const index_t nP = static_cast<index_t>(m.numPatches());
        std::ostringstream os;
        for (index_t k = 0; k != nP; ++k)
        {
            index_t loc = -1;
            os << (k ? " " : "") << (m.indexOnPatch(gl, k, loc) ? std::to_string(loc) : std::string("-"));
        }
        return os.str();
    }

    // Matching and elimination applied to one component of a mapper (the index is the component)
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

        rtFixture(index_t nx, index_t ny, index_t nPatches) : sizes(nPatches, 2)
        {
            const std::vector<boundaryInterface> ifs =
                gsNurbsCreator<real_t>::BSplineSquareGrid(nPatches, 1, 1.0).interfaces();
            for (index_t c = 0; c != 2; ++c)
            {
                const gsTensorBSplineBasis<2,real_t> bc(gsKnotVector<real_t>(0, 1, static_cast<unsigned>(nx - 1), c == 0 ? 4 : 3),
                                                        gsKnotVector<real_t>(0, 1, static_cast<unsigned>(ny - 1), c == 1 ? 4 : 3));
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

        // The same lists for component 1 (argument c = 1; the j argument stays -1 for "any patch").
        // Patch 0 has 4 local dofs with indices 7,8,9,11: local dof 3 (index 11) is coupled, local dofs 0-2
        // (indices 7-9) are free uncoupled. Patch 1 has 5 local dofs with indices 11,12,12,10,14: local
        // dofs 0,1,2 (indices 11,12,12) are coupled, local dof 3 (index 10) is free uncoupled and local
        // dof 4 (index 14) is eliminated; patch 0 has nothing eliminated in component 1.
        const vec cpl0c1 = {3};
        const vec cpl1c1 = {0,1,2};
        const vec fu0c1 = {0,1,2};
        const vec fu1c1 = {3};
        const vec bnd1c1 = {4};
        CHECK(toVec(m.findCoupled(0, -1, 1)) == cpl0c1);
        CHECK(toVec(m.findCoupled(1, -1, 1)) == cpl1c1);
        CHECK(toVec(m.findFreeUncoupled(0, 1)) == fu0c1);
        CHECK(toVec(m.findFreeUncoupled(1, 1)) == fu1c1);
        CHECK(toVec(m.findBoundary(1, 1)) == bnd1c1);
        CHECK_EQUAL(0, (index_t)m.findBoundary(0, 1).size());

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

    // Verifies that boundarySizeWithDuplicates() counts every patch-local
    // occurrence of an eliminated dof, in all components, while
    // boundarySize() counts each eliminated dof once. Expected values are
    // counted by hand from the M2 table and from a two-patch mapper whose
    // single eliminated dof is shared by both patches.
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

    // Verifies that markTagged() on a mapper that is not finalized throws
    // std::runtime_error. Tags are final dof numbers, which exist only after
    // finalize(), so the call has nothing to store.
    TEST(mark_tagged_needs_finalize)
    {
        // Not finalized: markTagged() has no final numbering to store, which is a GISMO_ENSURE
        // failure in every build
        gsVector<index_t> sz(2); sz << 4, 5;
        gsDofMapper m(sz, 2);
        CHECK_THROW(m.markTagged(1,0,0), std::runtime_error);
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
        // Every patch starts at 0, so only the last patch carries the 20 dofs
        for (index_t k = 0; k != 3; ++k)
        {
            CHECK_EQUAL(0, (index_t)m.offset(k));
            CHECK_EQUAL(k == 2 ? 20 : 0, (index_t)m.patchSize(k, 0));
            CHECK_EQUAL(k == 2 ? 20 : 0, (index_t)m.patchSize(k, 1));
        }
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

    // Verifies inverseOnPatch(k) of M2 at shift 13: the keys are the shift-0
    // keys of inverse_on_patch_shift0 plus 13, the values are unchanged local dofs.
    TEST(inverse_on_patch_shift13)
    {
        gsDofMapper m = makeM2();
        m.setShift(13);

        // Patch 0: the map of inverse_on_patch_shift0 with 13 added to every key
        // Keys are shifted global indices, values are local dofs of the given patch
        typedef std::map<index_t,index_t> imap;
        imap inv0 = m.inverseOnPatch(0);
        imap want0;
        want0[13] = 1; want0[14] = 2; want0[19] = 3; want0[20] = 0;
        want0[21] = 1; want0[22] = 2; want0[24] = 3; want0[26] = 0;
        CHECK(inv0 == want0);

        // Patch 1: the map of inverse_on_patch_shift0 with 13 added to every key; shared key 12 becomes 25
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
        CHECK_THROW_IN_DEBUG(m.is_free_index(12), std::logic_error);
        CHECK_THROW_IN_DEBUG(m.is_free_index(28), std::logic_error);
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
        CHECK_THROW_IN_DEBUG(m.is_boundary_index(12), std::logic_error);
        CHECK_THROW_IN_DEBUG(m.is_boundary_index(28), std::logic_error);
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
        CHECK_THROW_IN_DEBUG(m.global_to_bindex(12), std::logic_error);
        CHECK_THROW_IN_DEBUG(m.global_to_bindex(28), std::logic_error);
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
        CHECK_THROW_IN_DEBUG(m.componentOf(12), std::logic_error);
        CHECK_THROW_IN_DEBUG(m.componentOf(28), std::logic_error);
    }

    // Verifies that a finalized default-constructed mapper has no valid
    // global index: its size is 0, so every index query is out of range.
    TEST(empty_mapper_rejects_all)
    {
        // A finalized default-constructed mapper has size() == 0, so it has no valid global index
        gsDofMapper m;
        m.finalize();
        CHECK_EQUAL(0, m.size());
        CHECK_THROW_IN_DEBUG(m.is_free_index(0), std::logic_error);
        CHECK_THROW_IN_DEBUG(m.is_boundary_index(0), std::logic_error);
        CHECK_THROW_IN_DEBUG(m.global_to_bindex(0), std::logic_error);
    }

    // Verifies the layout and numbering of a mapper whose components have
    // different sizes on each patch (size matrix constructor): sizes, offsets,
    // the identity numbering without conditions, and the numbering after two
    // dofs of component 1 are eliminated. The oracle is the hand-counted table
    // below; no second gsDofMapper is involved.
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
        // The column sums are the totals 2+3, 5+1, 1+6; the mapSize 18 is
        // their sum 5+6+7 (equal to the sum of all six entries); the blocks
        // start at 0, 5 and 5+6 = 11.
        gsMatrix<index_t> M(2,3);
        M << 2,5,1,
             3,1,6;
        const index_t blockStart[3] = {0, 5, 11};
        gsDofMapper m(M);

        // Sizes and offsets before any condition, as in the table above
        CHECK_EQUAL(2, (index_t)m.numPatches());
        CHECK_EQUAL(18, (index_t)m.mapSize());
        CHECK_EQUAL(5, (index_t)m.totalSize(0));
        CHECK_EQUAL(6, (index_t)m.totalSize(1));
        CHECK_EQUAL(7, (index_t)m.totalSize(2));
        // patchSize(k,c) is the matrix entry; patch 0 starts at offset 0 and patch 1 right after it,
        // so offset(1,c) equals the size of patch 0 in component c
        for (index_t c = 0; c < 3; ++c)
        {
            CHECK_EQUAL(M(0,c), (index_t)m.patchSize(0,c));
            CHECK_EQUAL(M(1,c), (index_t)m.patchSize(1,c));
            CHECK_EQUAL(0, (index_t)m.offset(0,c));
            CHECK_EQUAL(M(0,c), (index_t)m.offset(1,c));
        }

        // Without conditions the numbering is the concatenation of the blocks: component c holds
        // the consecutive indices blockStart[c] ... blockStart[c] + totalSize(c) - 1
        gsDofMapper id(M);
        id.finalize();
        const vec want0 = {0,1,2,3,4};
        const vec want1 = {5,6,7,8,9,10};
        const vec want2 = {11,12,13,14,15,16,17};
        CHECK(toVec(id.asVector(0)) == want0);
        CHECK(toVec(id.asVector(1)) == want1);
        CHECK(toVec(id.asVector(2)) == want2);
        // Local dof i of patch k in component c has index blockStart[c] + offset(k,c) + i, and that
        // index belongs to component c
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
        // Free dofs: 18 - 2 = 16, indices 0..15; the eliminated dofs get 16 and 17. Component 1
        // on patch 0 is {16,5,6,7,17}, followed by its single dof 8 on patch 1.
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
        // Counters: size() still counts the eliminated dofs (18); freeSize() = 18 - 2 = 16; per component
        // 5 (nothing eliminated), 6 - 2 = 4 and 7 (nothing eliminated), which add up to 16
        CHECK_EQUAL(18, m.size());
        CHECK_EQUAL(16, m.freeSize());
        CHECK_EQUAL(2, m.boundarySize());
        CHECK_EQUAL(5, m.freeSize(0));
        CHECK_EQUAL(4, m.freeSize(1));
        CHECK_EQUAL(7, m.freeSize(2));
        CHECK_EQUAL(6, m.size(1));
        // The eliminated dofs are numbered 0, 1, ... among themselves: dof 0 is boundary dof 0, dof 4 is 1
        CHECK_EQUAL(0, m.bindex(0,0,1));
        CHECK_EQUAL(1, m.bindex(4,0,1));

        // Without matching, every dof keeps a distinct index: the concatenated numbering of the three
        // components is a permutation of 0..17
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

        // The index of every local dof, free or eliminated, lies in the range owned by its own component
        for (index_t c = 0; c < 3; ++c)
            for (index_t k = 0; k < 2; ++k)
                for (index_t i = 0; i < M(k,c); ++i)
                    CHECK_EQUAL(c, m.componentOf(m.index(i,k,c)));

        // mapIndex(n) is the n-th entry of the flat map, which is asVector(0), asVector(1), asVector(2)
        // one after the other, i.e. the vector "all" built above
        for (index_t n = 0; n < 18; ++n)
            CHECK_EQUAL(all[n], m.mapIndex(n));
    }

    // Verifies that the size matrix with three equal columns builds the
    // mapper of the size vector with three components. The mapper from the
    // vector is checked directly against a hand count of the indices; the one
    // from the matrix is pinned to it by the comparison, which alone would not
    // detect a fault in the common construction of the layout.
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

        // Hand count. Before conditions 3 + 2 = 5 dofs per component, 15 in all. Component 2 identifies
        // dof 1 of patch 0 with dof 0 of patch 1, which leaves 14 indices. Eliminated dofs keep their
        // index, so the eliminations do not change the size: 3 of them (dof 1 of patch 1 in
        // component 0, dofs 0 and 2 of patch 0 in component 1), hence 14 - 3 = 11 free dofs.
        CHECK_EQUAL(14, a.size());
        CHECK_EQUAL(11, a.freeSize());
        CHECK_EQUAL(3, a.boundarySize());
        CHECK_EQUAL(1, a.coupledSize());
        // Numbering: free dofs first, component by component, free uncoupled before coupled, then the
        // eliminated dofs in component order.
        //   component 0: free 0,1,2 (patch 0) and 3 (patch 1, dof 0); eliminated dof is 11
        //   component 1: free 4 (patch 0, dof 1), 5 and 6 (patch 1); eliminated dofs are 12 and 13
        //   component 2: free uncoupled 7 (patch 0, dof 0), 8 (patch 0, dof 2), 9 (patch 1, dof 1);
        //                the identified pair of dofs shares the last free index 10
        const vec want0 = {0,1,2,3,11};
        const vec want1 = {12,4,13,5,6};
        const vec want2 = {7,10,8,10,9};
        CHECK(toVec(a.asVector(0)) == want0);
        CHECK(toVec(a.asVector(1)) == want1);
        CHECK(toVec(a.asVector(2)) == want2);

        // Both constructors agree on all counters, layouts and index vectors
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

    // Verifies the component argument -1 ("every component"): it throws
    // std::runtime_error when a patch it touches has different sizes in the
    // components, and on patches of equal size it acts on all components.
    // The numbers are hand-counted from the size matrices below.
    TEST(DofMapperR3_minusOneComponent)
    {
        // Component -1 acts on the same local dof in every component; on a
        // ragged mapper it is only allowed for patches of equal size in all
        // components. Patch 0 has 3 dofs in both components, patch 1 has 2 and 4.
        gsMatrix<index_t> R(2,2);
        R << 3,3,
             2,4;
        // Eliminating dof 0 of patch 1 in all components: patch 1 has 2 dofs in component 0 and
        // 4 in component 1, so "the same local dof" is not defined there
        {
            gsDofMapper r(R);
            CHECK_THROW(r.eliminateDof(0, 1, -1), std::runtime_error);
        }
        // Matching dof 2 of patch 0 with dof 0 of patch 1 in all components: the same reason, patch 1
        // is involved
        {
            gsDofMapper r(R);
            CHECK_THROW(r.matchDof(0,2, 1,0, -1), std::runtime_error);
        }
        {
            // Patch 0 is allowed since both of its components have 3 dofs. Of the 3 + 2 + 3 + 4 = 12
            // dofs the identification of dofs 1 and 2 removes one per component, leaving 10 indices,
            // and two dofs (dof 0 of patch 0 in each component) are eliminated, leaving 8 free.
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
        // Dof 1 of patch 0 is eliminated and dof 2 of patch 0 is identified with dof 0 of patch 1,
        // in both components
        for (index_t c = 0; c < 2; ++c)
        {
            CHECK(m.is_boundary(1,0,c));
            CHECK_EQUAL(m.index(2,0,c), m.index(0,1,c));
        }
        // 2 components x (3 + 2) = 10 dofs, one identification per component leaves 8 indices;
        // one elimination per component gives 2 eliminated and 8 - 2 = 6 free ones
        CHECK_EQUAL(2, m.boundarySize());
        CHECK_EQUAL(6, m.freeSize());
        CHECK_EQUAL(8, m.size());
    }

    // Verifies that a default-constructed gsDofMapper has the documented
    // default layout: one patch, no component, and offsets {0,0} (see
    // m_offset in gsDofMapper.h), so that the offset and size of patch 0 are 0.
    TEST(DofMapperR4_defaultCtor)
    {
        // Construct in a poisoned buffer so that a member the default
        // constructor leaves uninitialised reads back as 0xAB bytes. The
        // poison is effective only at -O0 or with -fno-lifetime-dse: GCC may
        // drop stores to an object's storage made before its constructor runs.
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

    // Verifies that a size matrix with nc equal columns builds the mapper of
    // gsDofMapper(column, nc), also after coupling and elimination. The vector
    // mapper is checked directly against a hand count of the sizes; the matrix
    // mapper is pinned to it by the comparison (which alone would not detect a
    // fault shared by the constructors).
    TEST(matrix_equal_columns_as_vector)
    {
        // Three patches with 3, 0 and 4 dofs (the empty patch in the middle), called a0..a2 on patch 0
        // and b0..b3 on patch 2. The operations mix the argument -1 (every component), nc - 1 (the last
        // component) and 0 (the first). For nc = 1 all three name the same component, so every
        // operation lands on one mapper; for nc = 3 they are distinct, with component 1 touched only
        // by the -1 operations. That is why nc runs over 1 and 3.
        gsVector<index_t> sizes(3);
        sizes << 3, 0, 4;
        for (index_t nc = 1; nc <= 3; nc += 2)
        {
            gsDofMapper a(sizes, nc), b(gsMatrix<index_t>(sizes.replicate(1, nc)));
            for (gsDofMapper * m : {&a, &b})
            {
                m->matchDof(0, 2, 2, 0, -1);        // a2 ~ b0 in every component
                m->matchDof(0, 1, 2, 1, nc - 1);    // a1 ~ b1 in the last component
                m->matchDof(2, 1, 2, 3, 0);         // b1 ~ b3 in component 0
                m->eliminateDof(0, 0, -1);          // a0 eliminated in every component
                m->markBoundary(2, idx({0, 3}), nc - 1); // b0 and b3 eliminated in the last component
                m->markCoupled(2, 2, 0);            // b2 coupled in component 0
                m->finalize();
            }

            // Hand count of the indices per component (7 dofs before identification).
            // nc = 1 (all operations on component 0): classes {a2,b0}, {a1,b1,b3}, {a0}, {b2}, so size 4;
            // {a2,b0} and {a1,b1,b3} contain an eliminated dof and {a0} is eliminated: 3 eliminated,
            // free 4 - 3 = 1 (b2).
            // nc = 3: component 0 has {a2,b0}, {b1,b3}, {a1}, {b2}, {a0}: size 5 and only a0 eliminated,
            // free 4. Component 1 only has a2 ~ b0 and a0 eliminated: size 7 - 1 = 6, free 5.
            // Component 2 has {a2,b0} (eliminated through b0), {a1,b1}, {a0}, {b2}, {b3}: size
            // 7 - 2 = 5; eliminated {a2,b0}, a0 and b3, so free 2.
            const vec wantSize = nc == 1 ? vec{4} : vec{5, 6, 5};
            const vec wantFree = nc == 1 ? vec{1} : vec{4, 5, 2};
            for (index_t c = 0; c != nc; ++c)
            {
                CHECK_EQUAL(wantSize[c], a.size(c));
                CHECK_EQUAL(wantFree[c], a.freeSize(c));
                // The operation with argument -1 acts on every component
                CHECK_EQUAL(a.index(2, 0, c), a.index(0, 2, c));
                CHECK(a.is_boundary(0, 0, c));
            }
            CHECK_EQUAL(nc == 1 ? 4 : 16, a.size());
            CHECK_EQUAL(nc == 1 ? 1 : 11, a.freeSize());

            // The matrix constructor gives the same mapper in every observable respect
            CHECK_EQUAL(state(a), state(b));
        }
    }

    // Verifies that each component of a mapper built from a size matrix is
    // numbered like the scalar mapper of its column when the same operations
    // are applied (checkAgainstScalars compares counters, indices, preimages
    // and the flat map against the scalar mappers). Component 1 and component 3
    // each have an empty patch, and component 3 is also left untouched.
    TEST(matrix_components_as_scalar_mappers)
    {
        // Rows are the 3 patches, columns the 4 components
        gsMatrix<index_t> sizes(3, 4);
        sizes << 5, 4, 3, 2,
                 7, 0, 3, 2,
                 2, 6, 3, 0;
        // Array of mapper-operations to be applied to each component
        std::vector<compOps> ops;
        // Component 0 (patches of 5, 7, 2 dofs): two dofs of patch 0 matched with two of patch 1, an
        // elimination, a further match, a lone coupled mark and a match of patch 1 with patch 2
        ops.push_back([](gsDofMapper & m, index_t c) {
            m.matchDofs(0, idx({1, 3}), 1, idx({0, 2}), c);
            m.eliminateDof(0, 0, c);
            m.matchDof(0, 0, 1, 1, c);
            m.markCoupled(1, 2, c);
            m.matchDof(1, 6, 2, 0, c); });
        // Component 1 (4, 0, 6 dofs; patch 1 is empty): matches that skip the empty patch, an
        // elimination of a matched dof and a boundary mark on two dofs of patch 2
        ops.push_back([](gsDofMapper & m, index_t c) {
            m.matchDof(0, 3, 2, 0, c);
            m.eliminateDof(5, 2, c);
            m.matchDof(0, 2, 2, 5, c);
            m.markBoundary(2, idx({1, 2}), c); });
        // Component 2 (3 dofs on each patch): dof 1 eliminated on every patch, and a chain of matches
        // across the three patches
        ops.push_back([](gsDofMapper & m, index_t c) {
            for (index_t k = 0; k != 3; ++k)
                m.eliminateDof(1, k, c);
            m.matchDof(0, 2, 1, 0, c);
            m.matchDof(1, 2, 2, 0, c); });
        // Component 3 (2, 2, 0 dofs): no operation, an untouched component must stay the identity
        ops.push_back([](gsDofMapper &, index_t) { });
        const gsDofMapper m = checkAgainstScalars(sizes, ops);
        // Patch offsets are the running sums of the component 0 column: 0, 5, 5 + 7 = 12; patchSize(k,c)
        // lists the columns of the matrix one after the other, each separated by ';'
        CHECK_EQUAL("off=0 5 12 ps=5 7 2;4 0 6;3 3 3;2 2 0", layout(m));
    }

    // Verifies the per-component numbering for two bases of Raviart-Thomas
    // type on two patches of 4x2 elements: 28 and 30 dofs per patch (in the
    // two components respectively, on each patch)
    // (7 x 4 and 6 x 5 tensor-product functions), matched across the patch
    // interface and partly eliminated, compared with scalar mappers.
    TEST(matrix_raviart_thomas_bases)
    {
        const rtFixture rt(4, 2, 2);
        const gsDofMapper m = checkAgainstScalars(rt.sizes, rt.ops);
        // The interface of the two patches carries matched dofs, which are coupled; without them the
        // comparison would not exercise any matching. checkAgainstScalars establishes that each
        // component is numbered like the scalar mapper of its own basis column.
        CHECK(m.coupledSize() > 0);
    }

    // Verifies that the per-patch sizes decide the layout, not the component
    // totals (58 in both components).
    TEST(matrix_patch_sizes_differ_with_equal_totals)
    {
        // Component 0 has 28 and 30 dofs on the two patches, component 1 has 31 and 27: both total 58
        gsMatrix<index_t> sizes(2, 2);
        sizes << 28, 31,
                 30, 27;
        gsDofMapper m(sizes);
        m.finalize();
        CHECK_EQUAL("off=0 28 ps=28 30;31 27", layout(m));

        // Without conditions index i is the flat position. Component 0 owns 0..57: patch 0 holds 0..27
        // and patch 1 holds 28..57, so index 29 is dof 29 - 28 = 1 of patch 1. Component 1 owns 58..115:
        // patch 0 holds 58..88 (31 dofs) and patch 1 holds 89..115 (27 dofs). With equal splits of 29 per
        // patch these results would differ.
        std::vector<dofPair> pre;
        m.preImage(29, pre);
        CHECK_EQUAL("(1,1)", join(pre));
        // 58 + 29 is dof 29 of patch 0 in component 1
        m.preImage(58 + 29, pre);
        CHECK_EQUAL("(0,29)", join(pre));
        // 58 + 31 = 89 is the first index of patch 1 in component 1
        CHECK(dofPair(1, 0) == m.anyPreImage(58 + 31));
        // "29 -": local dof 29 on patch 0 and not on patch 1. "- 26": 58 + 57 = 115 is the last index
        // of component 1, local dof 115 - 89 = 26 of patch 1 and not on patch 0
        CHECK_EQUAL("29 -", onPatches(m, 58 + 29));
        CHECK_EQUAL("- 26", onPatches(m, 58 + 57));
    }

    // Verifies that the argument -1 (every component) throws on patches whose
    // sizes differ between components and leaves the mapper unchanged, and
    // that debug builds reject a local dof beyond its patch.
    TEST(matrix_unequal_components_reject_shared_calls)
    {
        // Patch 0 has 4 dofs in component 0 and 2 in component 1; patch 1 has 4 and 6. Both patches
        // differ between the components
        gsMatrix<index_t> sizes(2, 2);
        sizes << 4, 2,
                 4, 6;
        gsDofMapper m(sizes);
        // "The same local dof of every component" does not exist on these patches: GISMO_ENSURE
        // throws in every build
        CHECK_THROW(m.matchDof(0, 3, 1, 0, -1), std::runtime_error);
        CHECK_THROW(m.eliminateDof(0, 0, -1), std::runtime_error);
        // Range checks of the local dof are GISMO_ASSERTs, which exist only when NDEBUG is not
        // defined: they throw std::logic_error in debug builds and compile away otherwise
    #ifndef NDEBUG
        // Dof 6 of patch 1 in component 1 (6 dofs, valid 0..5)
        CHECK_THROW(m.eliminateDof(6, 1, 1), std::logic_error);
        CHECK_THROW(m.matchDof(1, 6, 0, 0, 1), std::logic_error);
        // Within the component, beyond the patch: patch 0 has only 2 dofs in component 1, but the
        // component holds 8
        CHECK_THROW(m.eliminateDof(3, 0, 1), std::logic_error);
        CHECK_THROW(m.matchDof(0, 2, 1, 0, 1), std::logic_error);
    #endif
        m.finalize();
        // Unchanged mapper: no identification and no elimination happened, so the indices are the
        // consecutive flat positions. Component 0 holds 0..3 (patch 0) and 4..7 (patch 1), component 1
        // holds 8,9 (patch 0) and 10..15 (patch 1)
        CHECK_EQUAL("off=0 4 ps=4 4;2 6 0,1,2,3|4,5,6,7;8,9|10,11,12,13,14,15", layout(m) + " " + indices(m));
    }

    // Verifies that invalid size matrices throw std::runtime_error from
    // the constructor.
    TEST(matrix_invalid_sizes_throw)
    {
        // The total number of dofs may be at most max() - 1: global indices span all components and
        // the coupling ids count from 1. big is that bound, so one more dof must be rejected
        const index_t big = std::numeric_limits<index_t>::max() - 1;
        gsMatrix<index_t> noRows(0, 2), noCols(2, 0), negative(2, 2), compTotal(2, 1), total(1, 2);
        negative << 3, 4, -1, 2;
        compTotal << big, 1;
        total << big, 1;
        // No patches
        CHECK_THROW(gsDofMapper m(noRows), std::runtime_error);
        // No components
        CHECK_THROW(gsDofMapper m(noCols), std::runtime_error);
        // A negative entry
        CHECK_THROW(gsDofMapper m(negative), std::runtime_error);
        // One component whose two patches sum to big + 1
        CHECK_THROW(gsDofMapper m(compTotal), std::runtime_error);
        // Two components (one patch) whose sizes sum to big + 1 over the whole mapper
        CHECK_THROW(gsDofMapper m(total), std::runtime_error);
    }

    // The matrix constructor is explicit; Eigen expressions still select the
    // vector constructor (one component, one patch per entry)
    static_assert(!std::is_convertible<gsMatrix<index_t>, gsDofMapper>::value,
                  "gsDofMapper(const gsMatrix<index_t> &) must be explicit");

    // Verifies which constructor each argument selects: a vector, or any Eigen
    // expression of vector shape, gives one component with one patch per
    // entry; a gsMatrix gives one patch per row and one component per column.
    TEST(matrix_constructor_overloads)
    {
        gsVector<index_t> v(2);
        v << 3, 4;
        gsMatrix<index_t> M(2, 2);
        M << 3, 5,
             4, 6;
        // Five vector-shaped arguments (a constant expression, an array expression, a column of a
        // matrix, a sum and a row vector), each of 2 entries, and one 2 x 2 matrix
        const gsDofMapper cst(gsVector<index_t>::Constant(2, 5)), plus(v.array() + 1), col(M.col(1)),
                          sum(v + v), row(v.transpose()), mat(M);
        // The vector-shaped ones have 1 component and 2 patches (2 entries); 2u is unsigned because
        // numPatches() returns size_t
        for (const gsDofMapper * m : {&cst, &plus, &col, &sum, &row})
        {
            CHECK_EQUAL(1, m->numComponents());
            CHECK_EQUAL(2u, m->numPatches());
        }
        // The 2 x 2 matrix has 2 rows (patches) and 2 columns (components)
        CHECK_EQUAL(2, mat.numComponents());
        CHECK_EQUAL(2u, mat.numPatches());
    }

    // Verifies that setIdentity() works on a moved-from mapper, which has no
    // offsets left.
    TEST(identity_after_move)
    {
        gsVector<index_t> sizes(2);
        sizes << 5, 7;
        gsDofMapper m(sizes);
        const gsDofMapper moved(give(m));
        m.setIdentity(2, 4, 1);
        m.finalize();
        // setIdentity(2, 4, 1) is 2 patches, 4 dofs, 1 component: all offsets are 0 and the last
        // patch holds all 4 dofs (patch 0 holds none), numbered 0..3
        CHECK_EQUAL("off=0 0 ps=0 4 |0,1,2,3", layout(m) + " " + indices(m));
    }

    // Verifies that after setIdentity() the dofs are addressed through any
    // patch, as callers with mapped bases do: an operation on patch 1 or patch 0
    // reaches the flat dofs of the one patch that holds them.
    TEST(identity_then_writes)
    {
        gsDofMapper m;
        m.setIdentity(3, 6, 2);
        m.markBoundary(1, idx({4}), 0);
        m.matchDof(0, 1, 2, 3, 1);
        m.finalize();
        // setIdentity(3, 6, 2): 3 patches, 6 dofs, 2 components, all of them held by patch 2.
        // Component 0: flat dof 4 is eliminated, the others take free indices 0,1,2,3,4 (dofs 0,1,2,3,5)
        // and dof 4 follows after all 10 free dofs of the mapper, index 10.
        // Component 1: flat dofs 1 and 3 are identified; its 5 classes are numbered from 5, free
        // uncoupled first (dofs 0,2,4,5 get 5,6,7,8) and the coupled class last, 9.
        CHECK_EQUAL("off=0 0 0 ps=0 0 6;0 0 6 ||0,1,2,3,10,4;||5,9,6,9,7,8", layout(m) + " " + indices(m));
    }

    // Verifies that the patch-wise queries (inverseOnPatch, indexOnPatch,
    // componentOf and the find* lists) of a ragged, shifted mapper agree with
    // the per-dof queries, in both components, and checks hand-derived lists.
    TEST(matrix_patch_queries_with_shift)
    {
        // Component 0 has 3 and 4 dofs on the patches, component 1 has 5 and 2.
        // Component 0: dof 0 of patch 0 eliminated, dof 2 of patch 0 identified with dof 0 of patch 1.
        // Component 1: dof 4 of patch 0 identified with dof 0 of patch 1, dof 1 of patch 1 eliminated.
        // The shift is set after finalize(); the tags follow, dof 1 of patch 0 in component 0 and
        // dof 1 of patch 1 (an eliminated dof) in component 1.
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
        // Patch 0 has no component 2: the range check is a GISMO_ASSERT, present in debug builds only
    #ifndef NDEBUG
        CHECK_THROW(m.patchSize(0, 2), std::logic_error);
    #endif

        for (index_t k = 0; k != 2; ++k)
        {
            // Every local dof of every component appears in inverseOnPatch(k) under its global index,
            // indexOnPatch returns its local index and componentOf the component it was taken from;
            // n counts the dofs of the patch over both components
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

            // Reference lists by brute force over the per-dof predicates, for each component, compared
            // with the find* functions of that component: free, boundary (not free), coupled,
            // free-uncoupled (free and not coupled) and tagged
            for (index_t c = 0; c != 2; ++c)
            {
                vec fr, bd, cp, fu, tg;
                for (index_t i = 0; i != static_cast<index_t>(m.patchSize(k, c)); ++i)
                {
                    (m.is_free(i, k, c) ? fr : bd).push_back(i);
                    if (m.is_coupled(i, k, c)) cp.push_back(i);
                    else if (m.is_free(i, k, c)) fu.push_back(i);
                    if (m.is_tagged(i, k, c)) tg.push_back(i);
                }
                CHECK(fr == toVec(m.findFree(k, c)));
                CHECK(bd == toVec(m.findBoundary(k, c)));
                CHECK(cp == toVec(m.findCoupled(k, -1, c)));
                CHECK(fu == toVec(m.findFreeUncoupled(k, c)));
                CHECK(tg == toVec(m.findTagged(k, c)));
            }
        }
        // Patch 0, component 0: dof 0 eliminated, dof 2 coupled with patch 1, dof 1 tagged
        CHECK(vec({1, 2}) == toVec(m.findFree(0)));
        CHECK(vec({0}) == toVec(m.findBoundary(0)));
        CHECK(vec({2}) == toVec(m.findCoupled(0, 1)));
        CHECK(vec({1}) == toVec(m.findFreeUncoupled(0)));
        CHECK(vec({1}) == toVec(m.findTagged(0)));
        // Component 1: patch 0 (5 dofs) has dof 4 coupled with patch 1 and dofs 0-3 free uncoupled;
        // patch 1 (2 dofs) has dof 0 coupled and dof 1 eliminated and tagged
        CHECK(vec({0, 1, 2, 3, 4}) == toVec(m.findFree(0, 1)));
        CHECK(vec({1}) == toVec(m.findBoundary(1, 1)));
        CHECK(vec({4}) == toVec(m.findCoupled(0, 1, 1)));
        CHECK(vec({0, 1, 2, 3}) == toVec(m.findFreeUncoupled(0, 1)));
        CHECK(vec({1}) == toVec(m.findTagged(1, 1)));
    }
}
