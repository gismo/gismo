/** @file gsDofMapper_test.cpp

    @brief Tests for gismo::gsDofMapper with ragged (per-component,
    per-patch varying) numbers of dofs.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s):

**/

#include "gismo_unittest.h"
#include <gsCore/gsDofMapper.h>
#include <set>

using namespace gismo;

namespace {

// The primary ragged fixture: 3 components over 3 patches, with
// pairwise-different, non-multiple totals (15/20/23) and per-patch sizes
// that also differ per component within every patch. This is deliberately
// NOT symmetric: a fixture built by replicating one basis over several
// components cannot discriminate a ragged offset/stride bug because every
// component then shares identical sizes and any error divides out.
std::vector<gsVector<index_t> > raggedSizes()
{
    std::vector<gsVector<index_t> > sz(3);
    sz[0].resize(3); sz[0] << 5, 7, 3;   // total 15
    sz[1].resize(3); sz[1] << 2, 11, 7;  // total 20
    sz[2].resize(3); sz[2] << 9, 4, 10;  // total 23
    return sz;
}

gsDofMapper raggedMapper()
{
    gsDofMapper m(raggedSizes());
    m.finalize();
    return m;
}

// Closed-form global index for a mapper built from raggedSizes() with no
// coupling and no elimination, derived from gsDofMapper::finalize() /
// finalizeComp() (gsDofMapper.cpp): finalizeComp(c) starts curFreeDof at
// flatOffset(c) = sum_{c'<c} totalSize(c') (cumulative prefix already in
// m_numFreeDofs[c] at that point, and nothing eliminated), then assigns
// dofs in flat (patch,local) order, i.e. flatOffset(c) + patchOffset(c,k) + i
// with patchOffset(c,k) = sum_{k'<k} patchSize(k',c).
index_t closedFormIndex(const std::vector<gsVector<index_t> > & sz,
                        index_t i, index_t k, index_t c)
{
    index_t flatOffset = 0;
    for (index_t cp = 0; cp < c; ++cp)
        flatOffset += sz[cp].sum();

    index_t patchOffset = 0;
    for (index_t kp = 0; kp < k; ++kp)
        patchOffset += sz[c][kp];

    return flatOffset + patchOffset + i;
}

} // anonymous namespace

SUITE(gsDofMapper_test)
{

// 1. index(i,k,c) is injective over the whole (c,k,i) space, its image is
//    exactly [0,size()), and it agrees with the closed form.
TEST(ragged_index_injective_and_onto)
{
    const std::vector<gsVector<index_t> > sz = raggedSizes();
    gsDofMapper m = raggedMapper();

    std::vector<int> seen(m.size(), 0);
    index_t nVisited = 0;

    for (index_t c = 0; c < 3; ++c)
        for (index_t k = 0; k < 3; ++k)
            for (index_t i = 0; i < sz[c][k]; ++i)
            {
                const index_t gl = m.index(i, k, c);
                CHECK(gl >= 0 && gl < static_cast<index_t>(seen.size()));
                ++seen[gl];
                ++nVisited;

                CHECK_EQUAL(closedFormIndex(sz, i, k, c), gl);
            }

    CHECK_EQUAL(58, nVisited);
    CHECK_EQUAL(58, m.size());
    for (size_t j = 0; j != seen.size(); ++j)
        CHECK_EQUAL(1, seen[j]);
}

// 2. componentOf(index(i,k,c)) == c for every (i,k,c) of the ragged fixture.
TEST(ragged_componentOf_roundtrip)
{
    const std::vector<gsVector<index_t> > sz = raggedSizes();
    gsDofMapper m = raggedMapper();

    for (index_t c = 0; c < 3; ++c)
        for (index_t k = 0; k < 3; ++k)
            for (index_t i = 0; i < sz[c][k]; ++i)
                CHECK_EQUAL(c, m.componentOf(m.index(i, k, c)));
}

// 3. preImage(index(i,k,c)) is exactly the singleton (k,i): the size()==1
//    half is load-bearing, since "result contains (k,i)" alone would also
//    pass on a preImage() that returns every pre-image of the mapper.
TEST(ragged_preImage)
{
    const std::vector<gsVector<index_t> > sz = raggedSizes();
    gsDofMapper m = raggedMapper();

    std::vector<std::pair<index_t,index_t> > result;
    for (index_t c = 0; c < 3; ++c)
        for (index_t k = 0; k < 3; ++k)
            for (index_t i = 0; i < sz[c][k]; ++i)
            {
                m.preImage(m.index(i, k, c), result);
                CHECK_EQUAL(1u, (unsigned)result.size());
                CHECK(result[0] == std::make_pair(k, i));
            }

    // The last local dof of the last patch of every component: this is
    // where an off-by-one in a ragged offset table hides.
    m.preImage(m.index(2, 2, 0), result);
    CHECK_EQUAL(1u, (unsigned)result.size());
    CHECK(result[0] == std::make_pair(index_t(2), index_t(2)));

    m.preImage(m.index(6, 2, 1), result);
    CHECK_EQUAL(1u, (unsigned)result.size());
    CHECK(result[0] == std::make_pair(index_t(2), index_t(6)));

    m.preImage(m.index(9, 2, 2), result);
    CHECK_EQUAL(1u, (unsigned)result.size());
    CHECK(result[0] == std::make_pair(index_t(2), index_t(9)));
}

// 4. patchSize(k,c), totalSize(c), numPatches(), componentsSize() and
//    offset(k,c) reproduce the fixture table exactly, including patch 2.
TEST(ragged_patchSize)
{
    gsDofMapper m = raggedMapper();

    CHECK_EQUAL(3u, (unsigned)m.numPatches());
    CHECK_EQUAL(3u, (unsigned)m.componentsSize());

    const index_t expectedSize[3][3] = { {5,7,3}, {2,11,7}, {9,4,10} };
    const index_t expectedTotal[3]   = {15,20,23};
    const size_t  expectedOffset[3][3] = { {0,5,12}, {0,2,13}, {0,9,13} };

    for (index_t c = 0; c < 3; ++c)
    {
        CHECK_EQUAL(expectedTotal[c], (index_t)m.totalSize(c));
        for (index_t k = 0; k < 3; ++k)
        {
            CHECK_EQUAL((size_t)expectedSize[c][k], m.patchSize(k, c));
            CHECK_EQUAL(expectedOffset[c][k], m.offset(k, c));
        }
    }
}

// 5. mapIndex(flatOffset(c) + j) == asVector(c)[j] over the full sweep of
//    every component. The sweep must reach j >= totalSize(0)==15 within
//    component 1 and cover all of component 2: with this fixture
//    flatOffset(1) == 15 == totalSize(0), so a stale n/patchDofSizes.front()
//    style implementation still agrees with the correct answer for all of
//    component 1's j < 15; divergence begins only past that point.
TEST(ragged_mapIndex_roundtrip)
{
    gsDofMapper m = raggedMapper();
    const index_t flatOffset[3] = {0, 15, 35};
    const index_t total[3]      = {15, 20, 23};

    for (index_t c = 0; c < 3; ++c)
    {
        const gsVector<index_t> vec = m.asVector(c);
        for (index_t j = 0; j < total[c]; ++j)
            CHECK_EQUAL(vec[j], m.mapIndex(flatOffset[c] + j));
    }
}

// 6. mapSize() sums the per-component totals, and the untouched ragged
//    mapper is a permutation (nothing eliminated).
TEST(ragged_mapSize_and_isPermutation)
{
    const std::vector<gsVector<index_t> > sz = raggedSizes();
    gsDofMapper m = raggedMapper();

    index_t sum = 0;
    for (index_t c = 0; c < 3; ++c)
        sum += sz[c].sum();

    CHECK_EQUAL(sum, (index_t)m.mapSize());
    CHECK_EQUAL(58u, (unsigned)m.mapSize());
    CHECK(m.isPermutation());
}

// 7. markBoundary on one (patch,component) pair must not leak into other
//    components sharing the same patch index, even though their patch
//    sizes differ (patch 1: 7 dofs in component 0, 11 in component 1).
TEST(ragged_markBoundary_component_isolation)
{
    gsDofMapper m(raggedSizes());

    gsMatrix<index_t> dofs(3, 1);
    dofs << 0, 3, 10;              // local dof 10 only exists on component 1's patch 1
    m.markBoundary(1, dofs, 1);
    m.finalize();

    CHECK_EQUAL(15, m.freeSize(0));
    CHECK_EQUAL(23, m.freeSize(2));
    CHECK_EQUAL(17, m.freeSize(1)); // 20 - 3
    CHECK_EQUAL(3,  m.boundarySize());

    std::set<index_t> bidx;
    const index_t eliminated[3] = {0, 3, 10};
    for (int t = 0; t < 3; ++t)
    {
        const index_t b = m.bindex(eliminated[t], 1, 1);
        CHECK(b >= 0 && b < m.boundarySize());
        bidx.insert(b);
    }
    CHECK_EQUAL(3u, (unsigned)bidx.size());

    const std::vector<gsVector<index_t> > sz = raggedSizes();

    // Component 1: eliminated exactly at patch 1, local dofs {0,3,10}.
    for (index_t k = 0; k < 3; ++k)
        for (index_t i = 0; i < sz[1][k]; ++i)
        {
            const bool expected = (k == 1) && (i == 0 || i == 3 || i == 10);
            CHECK_EQUAL(expected, m.is_boundary(i, k, 1));
        }

    // Components 0 and 2: nothing eliminated anywhere.
    for (index_t c = 0; c < 3; c += 2)
        for (index_t k = 0; k < 3; ++k)
            for (index_t i = 0; i < sz[c][k]; ++i)
                CHECK(!m.is_boundary(i, k, c));
}

// 8. Backward-compatibility gate: the uniform constructor (one basis
//    replicated over nComp components) must keep behaving exactly as
//    before the ragged rewrite of gsDofMapper's internals.
TEST(uniform_backward_compatibility)
{
    gsVector<index_t> sz(3);
    sz << 4, 6, 5;
    gsDofMapper m(sz, 3);
    m.finalize();

    const index_t off[3] = {0, 4, 10};
    const index_t compTotal = 15;

    CHECK_EQUAL(3u, (unsigned)m.numPatches());
    CHECK_EQUAL(3u, (unsigned)m.componentsSize());
    CHECK_EQUAL(45, m.size());
    CHECK_EQUAL(45, m.freeSize());
    CHECK_EQUAL(0, m.boundarySize());
    CHECK_EQUAL(0, m.coupledSize());
    CHECK_EQUAL(45u, (unsigned)m.mapSize());
    CHECK(m.isPermutation());

    for (index_t c = 0; c < 3; ++c)
    {
        for (index_t k = 0; k < 3; ++k)
        {
            CHECK_EQUAL((size_t)sz[k], m.patchSize(k, c));
            CHECK_EQUAL((size_t)off[k], m.offset(k, c));
        }

        const gsVector<index_t> vec = m.asVector(c);
        for (index_t j = 0; j < compTotal; ++j)
        {
            CHECK_EQUAL(vec[j], m.mapIndex(compTotal * c + j));
            for (index_t k = 0; k < 3; ++k)
                for (index_t i = 0; i < sz[k]; ++i)
                    if (off[k] + i == j)
                        CHECK_EQUAL(compTotal * c + off[k] + i, m.index(i, k, c));
        }
    }

    gsDofMapper d;
    CHECK_EQUAL(1u, (unsigned)d.numPatches());
}

// 9. Regression test for a heap-buffer overflow that ragged components
// expose in indexOnPatch(gl,k,local): scanning component-0's offset
// table (m_offset[0]) for every
// component, so for component 2 -- the longest component, total 23 vs
// component 0's 15 -- the scan ran past the end of m_dofs[0] whenever gl
// belonged to component 2 (Valgrind: "Invalid read of size 4 ... 4 bytes
// after a block of size 40"). Exercising every patch of the longest
// component pins this: it fails (or crashes) on the pre-fix code and
// passes cleanly once indexOnPatch uses componentOf(gl) to pick the
// right offset table (gsDofMapper.cpp:645-646).
TEST(ragged_indexOnPatch_longest_component_regression)
{
    const std::vector<gsVector<index_t> > sz = raggedSizes();
    gsDofMapper m = raggedMapper();

    const index_t longest = 2; // component 2, total 23, the largest of the three
    for (index_t k = 0; k < 3; ++k)
        for (index_t i = 0; i < sz[longest][k]; ++i)
        {
            const index_t gl = m.index(i, k, longest);

            index_t local = -1;
            CHECK(m.indexOnPatch(gl, k, local));
            CHECK_EQUAL(i, local);

            for (index_t kOther = 0; kOther < 3; ++kOther)
                if (kOther != k)
                    CHECK(!m.indexOnPatch(gl, kOther, local));
        }
}

// 10. anyPreImages(comp)'s contract is size() entries (NOT
//     totalSize(comp)), one per global dof, filled with (patch,local) for
//     dofs owned by comp and left at the sentinel (-1,0) otherwise. Also
//     exercises the one real caller pattern, gsMultiPatch::toMesh()
//     (gsCore/gsMultiPatch.hpp:670), which indexes anyPreImages()[j] for
//     every j < mapper.size() on a single-component mapper.
TEST(ragged_anyPreImages_contract_and_caller_pattern)
{
    const std::vector<gsVector<index_t> > sz = raggedSizes();
    gsDofMapper m = raggedMapper();
    const index_t total[3] = {15, 20, 23};

    for (index_t c = 0; c < 3; ++c)
    {
        const std::vector<std::pair<index_t,index_t> > pi = m.anyPreImages(c);
        CHECK_EQUAL(58u, (unsigned)pi.size());

        index_t filled = 0;
        for (size_t j = 0; j != pi.size(); ++j)
        {
            if (pi[j].first != -1)
                ++filled;
            else
                CHECK_EQUAL(0, pi[j].second);
        }
        CHECK_EQUAL(total[c], filled);
    }

    // Single-component caller pattern: anyPreImages()[j] must be in range
    // and agree with anyPreImage(j) for every global dof j < mapper.size().
    gsVector<index_t> single(2);
    single << 6, 9;
    gsDofMapper single_m(single);
    single_m.finalize();

    const std::vector<std::pair<index_t,index_t> > pi = single_m.anyPreImages();
    CHECK_EQUAL((size_t)single_m.size(), pi.size());
    for (index_t j = 0; j < single_m.size(); ++j)
    {
        CHECK(pi[j].first >= 0 && pi[j].first < (index_t)single_m.numPatches());
        CHECK(pi[j] == single_m.anyPreImage(j));
    }
}

// 11. findTagged has no prior behaviour to regress against, so its oracle is
//     built from first principles: tag a known set with markTagged (NOT
//     markCoupledAsTagged,
//     which tags a set disjoint from the coupled dofs -- a pre-existing,
//     deliberately untouched defect that would give a wrong expectation
//     here), then assert findTagged returns exactly that set per (k,c).
//     Ordering trap: markTagged goes through index(), so it must be called
//     AFTER finalize().
TEST(ragged_findTagged_explicit_markTagged)
{
    gsDofMapper m = raggedMapper(); // already finalized

    m.markTagged(0, 1, 1);
    m.markTagged(3, 1, 1);
    m.markTagged(10, 1, 1);  // component 1's patch 1 has 11 dofs; the
                             // ragged discriminator (component 0's patch 1
                             // has only 7)
    m.markTagged(2, 2, 2);

    const gsVector<index_t> tagged11 = m.findTagged(1, 1);
    CHECK_EQUAL(3, tagged11.size());
    CHECK_EQUAL(0,  tagged11[0]);
    CHECK_EQUAL(3,  tagged11[1]);
    CHECK_EQUAL(10, tagged11[2]);

    const gsVector<index_t> tagged22 = m.findTagged(2, 2);
    CHECK_EQUAL(1, tagged22.size());
    CHECK_EQUAL(2, tagged22[0]);

    // Every other (k,c) pair must report no tagged dofs at all.
    for (index_t c = 0; c < 3; ++c)
        for (index_t k = 0; k < 3; ++k)
            if (!(k == 1 && c == 1) && !(k == 2 && c == 2))
                CHECK_EQUAL(0, m.findTagged(k, c).size());
}

}
