/** @file gsFeSpaceMapper_test.cpp

    @brief Regression tests for gsFeSpace::setupMapper (B0): an explicitly
    installed gsDofMapper must survive a gsExprAssembler::initSystem() cycle
    unchanged, instead of being silently discarded and rebuilt by
    gsFeSpaceData::init().

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s):

**/

#include "gismo_unittest.h"
#include <gsAssembler/gsExprAssembler.h>
#include <gsAssembler/gsDofMapperCreator.h>

using namespace gismo;

namespace {

// 2 unit squares side by side, one interface (patch0 east <-> patch1 west).
gsMultiBasis<real_t> twoPatchBasis()
{
    gsMultiPatch<real_t> mp = gsNurbsCreator<real_t>::BSplineSquareGrid(2, 1, 1.0);
    gsMultiBasis<real_t> mb(mp);
    mb.degreeElevate(1);
    mb.uniformRefine(1);
    return mb;
}

// Ragged mapper: 3 components over the same 2 patches, with pairwise-different,
// non-multiple per-patch sizes (component totals 34 / 40 / 42), so that
// mapSize() != dim*source().size() and every patchSize(k,c) differs from
// every other -- a per-(k,c) mismatch cannot hide behind an accidentally
// symmetric fixture.
gsDofMapper raggedMapper()
{
    std::vector<gsVector<index_t> > sizes(3);
    sizes[0].resize(2); sizes[0] << 20, 14;
    sizes[1].resize(2); sizes[1] <<  9, 31;
    sizes[2].resize(2); sizes[2] << 25, 17;
    gsDofMapper mapper(sizes);
    mapper.finalize();
    return mapper;
}

} // anonymous namespace

SUITE(gsFeSpaceMapper_test)
{

TEST(ragged_mapper_survives_initSystem)
{
    gsMultiBasis<real_t> mb = twoPatchBasis();
    gsDofMapper ragged = raggedMapper();

    // Fixture self-check: a broken Phase A must show up here, not in the
    // assembler round-trip below.
    CHECK_EQUAL(3, ragged.numComponents());
    CHECK_EQUAL(size_t(2), ragged.numPatches());
    for (index_t c = 0; c < 3; ++c)
        for (index_t k = 0; k < 2; ++k)
        {
            gsVector<index_t> expected(2);
            expected << (c==0 ? 20 : (c==1 ? 9 : 25)), (c==0 ? 14 : (c==1 ? 31 : 17));
            CHECK_EQUAL(static_cast<size_t>(expected[k]), ragged.patchSize(k,c));
        }
    CHECK_EQUAL(0, ragged.boundarySize());

    const size_t mapSizeBefore   = ragged.mapSize();
    const index_t sizeBefore     = ragged.size();
    const index_t freeSizeBefore = ragged.freeSize();
    index_t freeSizeCompBefore[3];
    size_t  patchSizeBefore[2][3];
    for (index_t c = 0; c < 3; ++c)
    {
        freeSizeCompBefore[c] = ragged.freeSize(c);
        for (index_t k = 0; k < 2; ++k)
            patchSizeBefore[k][c] = ragged.patchSize(k,c);
    }

    gsExprAssembler<real_t> A(1, 1);
    A.setIntegrationElements(mb);
    // Space dim must match the mapper's component count (3), or setupMapper's
    // GISMO_ENSURE rejects the install outright -- see
    // mismatched_component_mapper_rejected below.
    auto u = A.getSpace(mb, 3);
    u.setupMapper(ragged);

    A.initSystem();

    // The B0 gap this task pins: gsFeSpaceData::init() always builds with
    // conforming=false and has NO size check for a ragged mapper (Phase A
    // ragged components are not yet coverable by the uniform size-identity
    // check in setupMapper), so an installed ragged mapper is currently
    // accepted with no validation at all. The correct check
    // (mapSize()==Sum_c source(c).size()) is not expressible until B1
    // introduces per-component sources; until then this test documents the
    // "installed survives, unchecked" contract, not a validated one.
    CHECK(u.mapper().isFinalized());
    CHECK_EQUAL(3, (index_t)u.mapper().componentsSize());
    CHECK_EQUAL(mapSizeBefore, u.mapper().mapSize());
    CHECK_EQUAL(sizeBefore, u.mapper().size());
    CHECK_EQUAL(freeSizeBefore, u.mapper().freeSize());
    for (index_t c = 0; c < 3; ++c)
    {
        CHECK_EQUAL(freeSizeCompBefore[c], u.mapper().freeSize(c));
        for (index_t k = 0; k < 2; ++k)
            CHECK_EQUAL(patchSizeBefore[k][c], u.mapper().patchSize(k,c));
    }

    // The guard must hold across repeated resetDimensions() calls, not just
    // the first.
    A.initSystem();

    CHECK(u.mapper().isFinalized());
    CHECK_EQUAL(3, (index_t)u.mapper().componentsSize());
    CHECK_EQUAL(mapSizeBefore, u.mapper().mapSize());
    CHECK_EQUAL(sizeBefore, u.mapper().size());
    CHECK_EQUAL(freeSizeBefore, u.mapper().freeSize());
    for (index_t c = 0; c < 3; ++c)
    {
        CHECK_EQUAL(freeSizeCompBefore[c], u.mapper().freeSize(c));
        for (index_t k = 0; k < 2; ++k)
            CHECK_EQUAL(patchSizeBefore[k][c], u.mapper().patchSize(k,c));
    }
}

TEST(uniform_space_mapper_built_by_init)
{
    gsMultiBasis<real_t> mb = twoPatchBasis();
    gsExprAssembler<real_t> A(1, 1);
    A.setIntegrationElements(mb);
    auto u = A.getSpace(mb, 2);

    // Neither setup() nor setupMapper() is called: the "explicitly installed"
    // flag must default to false so init() still builds the mapper normally.
    A.initSystem();

    index_t total = 0;
    for (size_t k = 0; k < mb.nBases(); ++k)
        total += mb.basis(k).size();

    CHECK(u.mapper().isFinalized());
    CHECK_EQUAL(static_cast<size_t>(2*total), u.mapper().mapSize());
    CHECK_EQUAL(2, (index_t)u.mapper().numComponents());
    CHECK_EQUAL((index_t)mb.nBases(), (index_t)u.mapper().numPatches());
    for (index_t c = 0; c < 2; ++c)
        for (size_t k = 0; k < mb.nBases(); ++k)
            CHECK_EQUAL(static_cast<size_t>(mb.basis(k).size()), u.mapper().patchSize(k,c));
    CHECK_EQUAL(2*total, A.numDofs());
}

TEST(uniform_installed_mapper_unchanged)
{
    gsMultiBasis<real_t> mb = twoPatchBasis();
    gsBoundaryConditions<real_t> bcs; // empty: no eliminated dof anywhere

    // Reference path: mapper built the normal way, via setup().
    gsExprAssembler<real_t> Aref(1, 1);
    Aref.setIntegrationElements(mb);
    auto uref = Aref.getSpace(mb, 2);
    uref.setup(bcs, dirichlet::homogeneous, /*_icont=*/0);
    Aref.initSystem();
    gsDofMapper expected = uref.mapper();

    // Fixture precondition and the reason this test can fail at all:
    // _icont==0 makes setup() build with conforming=true (matches the
    // interface), whereas gsFeSpaceData::init() always builds with
    // conforming=false. An accidental rebuild collapses coupledSize() to 0
    // and changes every asVector(c), so this is the discriminating oracle,
    // not a vacuous re-read of componentsSize().
    CHECK(expected.coupledSize() > 0);

    // Install path: an independent copy of the same mapper, installed
    // explicitly on a fresh space.
    gsExprAssembler<real_t> A(1, 1);
    A.setIntegrationElements(mb);
    auto u = A.getSpace(mb, 2);
    gsDofMapper installed = expected;
    u.setupMapper(installed);
    A.initSystem();

    for (index_t c = 0; c < 2; ++c)
    {
        gsVector<index_t> vExpected = expected.asVector(c);
        gsVector<index_t> vGot      = u.mapper().asVector(c);
        CHECK_EQUAL(vExpected.size(), vGot.size());
        for (index_t i = 0; i < vExpected.size(); ++i)
            CHECK_EQUAL(vExpected[i], vGot[i]);
    }
    CHECK_EQUAL(expected.freeSize(), u.mapper().freeSize());
    CHECK_EQUAL(expected.boundarySize(), u.mapper().boundarySize());
    CHECK_EQUAL(expected.coupledSize(), u.mapper().coupledSize());
    CHECK_EQUAL(expected.taggedSize(), u.mapper().taggedSize());
    CHECK_EQUAL(expected.size(), u.mapper().size());
    CHECK_EQUAL(expected.mapSize(), u.mapper().mapSize());

    CHECK_EQUAL(Aref.numDofs(), A.numDofs());
    CHECK_EQUAL(Aref.numTestDofs(), A.numTestDofs());
}

// setupMapper's component-count check is a GISMO_ENSURE, not a GISMO_ASSERT,
// so it survives NDEBUG: a mapper whose numComponents() disagrees with the
// space's dim() is rejected in every build configuration rather than being
// silently rebuilt by init().
TEST(mismatched_component_mapper_rejected)
{
    gsMultiBasis<real_t> mb = twoPatchBasis();
    gsDofMapper ragged = raggedMapper(); // 3 components

    gsExprAssembler<real_t> A(1, 1);
    A.setIntegrationElements(mb);
    auto u = A.getSpace(mb, /*dim=*/2); // dim disagrees with ragged's 3 components

    CHECK_THROW(u.setupMapper(ragged), std::exception);
}

} // SUITE(gsFeSpaceMapper_test)
