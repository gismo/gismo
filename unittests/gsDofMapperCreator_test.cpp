/** @file gsDofMapperCreator_test.cpp

    @brief Tests for gismo::createMapper (gsAssembler/gsDofMapperCreator.h)

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s):

**/

#include "gismo_unittest.h"
#include <gsAssembler/gsDofMapperCreator.h>
#include <gsMSplines/gsMappedBasis.h>

using namespace gismo;

namespace {

// 2 unit squares side by side, one interface (patch0 east <-> patch1 west).
// Degree elevated + refined so that the interface carries several dofs.
gsMultiBasis<real_t> twoPatchBasis()
{
    gsMultiPatch<real_t> mp = gsNurbsCreator<real_t>::BSplineSquareGrid(2, 1, 1.0);
    gsMultiBasis<real_t> mb(mp);
    mb.degreeElevate(1);   // bilinear -> biquadratic
    mb.uniformRefine(3);
    return mb;
}

// Replays the conforming loop of createMapper: number of identified dof pairs
// over all (non-contact) interfaces of the given topology.
index_t countMatchedPairs(const gsMultiBasis<real_t> & mb, const gsBoxTopology & topology)
{
    index_t nMatched = 0;
    gsMatrix<index_t> b1, b2;
    for (gsBoxTopology::const_iiterator it = topology.iBegin(); it != topology.iEnd(); ++it)
    {
        if (it->type() == interaction::contact) continue;
        const gsBasis<real_t> & basis1 = mb.basis(it->first().patch);
        const gsBasis<real_t> & basis2 = mb.basis(it->second().patch);
        basis1.matchWith(*it, basis2, b1, b2);
        nMatched += b1.rows();
    }
    return nMatched;
}

} // anonymous namespace


SUITE(gsDofMapperCreator_test)
{

// 1. Single gsBasis, several components: plain identity mapper.
TEST(single_basis)
{
    gsMultiBasis<real_t> mb = twoPatchBasis();
    const gsBasis<real_t> & b = mb.basis(0);

    gsDofMapper m = createMapper(b, 2);
    m.finalize();

    CHECK_EQUAL(1u, (unsigned)m.numPatches());
    CHECK_EQUAL(static_cast<index_t>(b.size()*2), m.size());
    CHECK_EQUAL(m.size(), m.freeSize());
    CHECK_EQUAL(0, m.boundarySize());
    CHECK_EQUAL(0, m.coupledSize());
    CHECK(m.isPermutation());

    // Parity with the surviving gsDofMapper constructor
    gsVector<index_t> sz(1);
    sz[0] = b.size();
    gsDofMapper ref(sz, 2);
    ref.finalize();

    CHECK_EQUAL(ref.size()        , m.size()        );
    CHECK_EQUAL(ref.freeSize()    , m.freeSize()    );
    CHECK_EQUAL(ref.boundarySize(), m.boundarySize());
    CHECK_EQUAL((unsigned)ref.numPatches(), (unsigned)m.numPatches());
}

// 2. Multipatch, conforming == false: nothing is identified.
TEST(multipatch_nonconforming)
{
    gsMultiBasis<real_t> mb = twoPatchBasis();

    gsDofMapper m = createMapper(mb, 1, /*conforming=*/false, /*finalize=*/true);

    CHECK_EQUAL(static_cast<index_t>(mb.totalSize()), m.size());
    CHECK_EQUAL(0, m.coupledSize());
}

// 3. Multipatch, conforming == true: exactly the matchWith pairs are identified.
TEST(multipatch_conforming)
{
    gsMultiBasis<real_t> mb = twoPatchBasis();

    const index_t nMatched = countMatchedPairs(mb, mb.topology());
    CHECK(nMatched > 0);

    gsDofMapper m = createMapper(mb, 1, /*conforming=*/true, /*finalize=*/true);

    CHECK_EQUAL(static_cast<index_t>(mb.totalSize()) - nMatched, m.size());
    CHECK(m.size() < static_cast<index_t>(mb.totalSize()));
    CHECK(m.coupledSize() > 0);
}

// 4. The finalize flag (the semantics that differed between the removed
//    gsDofMapper constructors and the removed gsMultiBasis::getMapper).
TEST(finalize_flag)
{
    gsMultiBasis<real_t> mb = twoPatchBasis();

    gsDofMapper mUnfinal = createMapper(mb, 1, true, /*finalize=*/false);
    gsDofMapper mFinal   = createMapper(mb, 1, true, /*finalize=*/true );

    CHECK(!mUnfinal.isFinalized());
    CHECK( mFinal  .isFinalized());

    mUnfinal.finalize();
    CHECK(mUnfinal.isFinalized());
    CHECK_EQUAL(mFinal.size()    , mUnfinal.size()    );
    CHECK_EQUAL(mFinal.freeSize(), mUnfinal.freeSize());
}

// 5. Dirichlet boundary conditions are eliminated.
TEST(dirichlet_elimination)
{
    gsMultiBasis<real_t> mb = twoPatchBasis();

    gsFunctionExpr<real_t> g("0", 2);
    gsBoundaryConditions<real_t> bc;
    bc.addCondition(0, boundary::west, condition_type::dirichlet, &g);

    gsDofMapper m = createMapper(mb, bc, 1, 0, /*conforming=*/true, /*finalize=*/true);

    const index_t nb = mb.basis(0).boundary(boundary::west).rows();
    CHECK(nb > 0);
    CHECK_EQUAL(nb, m.boundarySize());
    CHECK_EQUAL(m.size() - m.boundarySize(), m.freeSize());
}

// 6. Contact interfaces must NOT be glued by the conforming loop.
TEST(contact_interface_skipped)
{
    gsMultiBasis<real_t> mb = twoPatchBasis();
    CHECK_EQUAL(1u, (unsigned)mb.topology().nInterfaces());

    const boundaryInterface bi0 = *mb.topology().iBegin();

    // (a) contact interface: nothing may be identified
    gsBoxTopology topoContact(2, static_cast<index_t>(mb.nBases()));
    boundaryInterface biContact = bi0;
    biContact.setAsContact();
    topoContact.addInterface(biContact);
    topoContact.addAutoBoundaries();

    gsDofMapper mContact = createMapper(mb, topoContact, 1, /*conforming=*/true, /*finalize=*/true);
    CHECK_EQUAL(static_cast<index_t>(mb.totalSize()), mContact.size());
    CHECK_EQUAL(0, mContact.coupledSize());

    // (b) control: the very same topology with a conforming interface DOES shrink
    gsBoxTopology topoConf(2, static_cast<index_t>(mb.nBases()));
    topoConf.addInterface(bi0);
    topoConf.addAutoBoundaries();

    gsDofMapper mConf = createMapper(mb, topoConf, 1, /*conforming=*/true, /*finalize=*/true);
    CHECK(mConf.size() < static_cast<index_t>(mb.totalSize()));
    CHECK_EQUAL(static_cast<index_t>(mb.totalSize()) - countMatchedPairs(mb, topoConf),
                mConf.size());
    CHECK(mConf.coupledSize() > 0);
}

// 7. The strategy-enum overload agrees with the boolean form, and drops the
//    boundary conditions for any non-elimination strategy.
TEST(strategy_overload_parity)
{
    gsMultiBasis<real_t> mb = twoPatchBasis();

    gsFunctionExpr<real_t> g("0", 2);
    gsBoundaryConditions<real_t> bc;
    bc.addCondition(0, boundary::west, condition_type::dirichlet, &g);

    gsDofMapper mElim = createMapper(mb, bc, dirichlet::elimination, iFace::glue, 1, 0, true);
    gsDofMapper mBool = createMapper(mb, bc, 1, 0, /*conforming=*/true, /*finalize=*/true);

    CHECK_EQUAL(mBool.size()        , mElim.size()        );
    CHECK_EQUAL(mBool.freeSize()    , mElim.freeSize()    );
    CHECK_EQUAL(mBool.boundarySize(), mElim.boundarySize());
    CHECK(mElim.boundarySize() > 0); // the bc is really there in this branch

    // Non-elimination strategy: the bc is ignored entirely
    gsDofMapper mNitsche = createMapper(mb, bc, dirichlet::nitsche, iFace::glue, 1, 0, true);
    gsDofMapper mFree    = createMapper(mb, 1, /*conforming=*/true, /*finalize=*/true);

    CHECK_EQUAL(0, mNitsche.boundarySize());
    CHECK_EQUAL(mFree.size()    , mNitsche.size()    );
    CHECK_EQUAL(mFree.freeSize(), mNitsche.freeSize());
}

// 8. unk == -1: all unknown indices are affected by Dirichlet conditions.
TEST(unknown_minus_one)
{
    gsMultiBasis<real_t> mb = twoPatchBasis();
    const index_t nbWest = mb.basis(0).boundary(boundary::west).rows();

    gsFunctionExpr<real_t> g("0", 2);
    gsBoundaryConditions<real_t> bc;

    // Set Dirichlet on patch 0, west boundary for unk=1 (not 0)
    bc.addCondition(0, boundary::west, condition_type::dirichlet, &g, 1);

    // unk == 0: must skip the condition (unknown 0 is not matched)
    gsDofMapper m0 = createMapper(mb, bc, 1, 0, false, true);
    CHECK_EQUAL(0, m0.boundarySize());

    // Add a second condition for unk=3
    bc.addCondition(0, boundary::east, condition_type::dirichlet, &g, 3);
    const index_t nbEast = mb.basis(0).boundary(boundary::east).rows();

    // unk == 1: only eliminates unk=1 (west)
    gsDofMapper m1 = createMapper(mb, bc, 1, 1, false, true);
    CHECK_EQUAL(nbWest, m1.boundarySize());

    // unk == -1: eliminates ALL unknowns (both west and east)
    gsDofMapper mAll = createMapper(mb, bc, 1, -1, false, true);
    CHECK_EQUAL(nbWest + nbEast, mAll.boundarySize());
}

// 9. Multipatch with multiple components: conforming matching per component.
TEST(multipatch_multicomponent_conforming)
{
    gsMultiBasis<real_t> mb = twoPatchBasis();
    const index_t nMatchedSingle = countMatchedPairs(mb, mb.topology());
    CHECK(nMatchedSingle > 0);

    gsDofMapper m = createMapper(mb, 3, true, true);

    CHECK_EQUAL(3 * (static_cast<index_t>(mb.totalSize()) - nMatchedSingle), m.size());
    CHECK_EQUAL(3u, (unsigned)m.componentsSize());
    CHECK(m.coupledSize() > 0);
    // The coupling is replicated per component
    CHECK_EQUAL(3 * nMatchedSingle, m.coupledSize());
}

// 10. gsMappedBasis with Dirichlet BCs exercises the setIdentity path.
TEST(mapped_basis_with_bcs)
{
    auto geom = gsNurbsCreator<real_t>::BSplineSquare(2);
    gsMultiBasis<real_t> mb(geom->basis());
    mb.basis(0).uniformRefine(2);
    const index_t sz = mb.basis(0).size();

    // Create a trivial (identity) mapped basis
    gsSparseMatrix<real_t> ident(sz, sz);
    ident.setIdentity();
    gsMappedBasis<2, real_t> mapB(mb, ident);

    gsFunctionExpr<real_t> g("0", 2);
    gsBoundaryConditions<real_t> bc;
    bc.addCondition(0, boundary::west, condition_type::dirichlet, &g);

    gsDofMapper m = createMapper(mapB, bc, 1, 0, false, true);

    CHECK_EQUAL(static_cast<index_t>(sz), m.size());
    const index_t nb = mb.basis(0).boundary(boundary::west).rows();
    CHECK_EQUAL(nb, m.boundarySize());
    CHECK_EQUAL(m.size() - m.boundarySize(), m.freeSize());
}

// Tests 12-16 (ragged createMapper overload coverage) are inserted here,
// between tests 10 and 11, so they sit beside the fixture helpers above;
// numbering continues after test 11 below.
// 12. Ragged primary overload (vector<const gsFunctionSet<T>*>): two
//     components with genuinely different basis sizes must produce a mapper
//     whose per-patch sizes differ per component, and whose per-component
//     eliminated count matches that component's own basis boundary size.
TEST(ragged_two_components_different_sizes)
{
    gsMultiBasis<real_t> mb0 = twoPatchBasis();
    gsMultiBasis<real_t> mb1 = twoPatchBasis();
    mb1.basis(0).uniformRefine(); // component 1's patch 0 only: now ragged

    CHECK(mb0.basis(0).size() != mb1.basis(0).size());

    gsFunctionExpr<real_t> g("0", 2);
    gsBoundaryConditions<real_t> bc;
    bc.addCondition(0, boundary::west, condition_type::dirichlet, &g, 0, false, 0); // component 0
    bc.addCondition(1, boundary::east, condition_type::dirichlet, &g, 0, false, 1); // component 1

    std::vector<const gsFunctionSet<real_t>*> basesPerComp(2);
    basesPerComp[0] = &mb0;
    basesPerComp[1] = &mb1;

    gsDofMapper m = createMapper(basesPerComp, mb0.topology(), bc, /*unk=*/-1,
                                 /*conforming=*/false, /*finalize=*/true);

    CHECK(m.patchSize(0, 0) != m.patchSize(0, 1));
    CHECK_EQUAL((size_t)mb0.basis(0).size(), m.patchSize(0, 0));
    CHECK_EQUAL((size_t)mb1.basis(0).size(), m.patchSize(0, 1));

    // Per-component eliminated count (valid because Dirichlet elimination
    // leaves the mapper with actual boundary dofs, i.e. !allFree()).
    const index_t nbWest0 = mb0.basis(0).boundary(boundary::west).rows();
    const index_t nbEast1 = mb1.basis(1).boundary(boundary::east).rows();
    CHECK(!m.allFree());
    CHECK_EQUAL(nbWest0, m.size(0) - m.freeSize(0));
    CHECK_EQUAL(nbEast1, m.size(1) - m.freeSize(1));
}

// 13. The vector<gsMultiBasis<T>> convenience overload agrees with the
//     primary pointer-based overload on the same inputs, WITH conforming=true.
//     The overload's only distinguishing job is deriving the topology from
//     basesPerComp.front().topology() (gsDofMapperCreator.hpp), and topology
//     is read only inside createMapper's `if (conforming)` branch, so a
//     conforming=false comparison (as in the earlier draft of this test)
//     cannot exercise that branch at all. The input here is uniform (same
//     basis in both components, no uniformRefine) so matchDofs's equal-size
//     assert does not fire.
TEST(ragged_multibasis_overload_agrees_with_primary)
{
    gsMultiBasis<real_t> mb0 = twoPatchBasis();
    gsMultiBasis<real_t> mb1 = twoPatchBasis();

    gsFunctionExpr<real_t> g("0", 2);
    gsBoundaryConditions<real_t> bc;
    bc.addCondition(0, boundary::west, condition_type::dirichlet, &g, 0, false, -1);

    std::vector<const gsFunctionSet<real_t>*> ptrs(2);
    ptrs[0] = &mb0;
    ptrs[1] = &mb1;
    gsDofMapper mPrimary = createMapper(ptrs, mb0.topology(), bc, /*unk=*/-1,
                                        /*conforming=*/true, /*finalize=*/true);

    std::vector<gsMultiBasis<real_t> > vecMb(2);
    vecMb[0] = mb0;
    vecMb[1] = mb1;
    gsDofMapper mConv = createMapper(vecMb, bc, /*unk=*/-1,
                                     /*conforming=*/true, /*finalize=*/true);

    for (index_t c = 0; c < 2; ++c)
    {
        const gsVector<index_t> a = mPrimary.asVector(c);
        const gsVector<index_t> b = mConv.asVector(c);
        CHECK(a.size() > 0);
        CHECK_EQUAL(a.size(), b.size());
        CHECK_ARRAY_EQUAL(a.data(), b.data(), a.size());
    }
    CHECK_EQUAL(mPrimary.freeSize()    , mConv.freeSize()    );
    CHECK_EQUAL(mPrimary.boundarySize(), mConv.boundarySize());
    CHECK_EQUAL(mPrimary.size()        , mConv.size()        );
}

// 14. A uniform input (the same basis replicated per component) run through
//     the new ragged overload must reproduce the OLD single-basis overload
//     exactly: the strongest available check that the new numbering path
//     did not quietly diverge from the surviving one.
TEST(ragged_uniform_input_matches_old_overload)
{
    gsMultiBasis<real_t> mb = twoPatchBasis();

    std::vector<const gsFunctionSet<real_t>*> basesPerComp(2, &mb);
    gsBoundaryConditions<real_t> bc; // empty

    gsDofMapper mRagged = createMapper(basesPerComp, mb.topology(), bc, /*unk=*/-1,
                                       /*conforming=*/true, /*finalize=*/true);
    gsDofMapper mOld     = createMapper(mb, /*nComp=*/2, /*conforming=*/true,
                                       /*finalize=*/true);

    for (index_t c = 0; c < 2; ++c)
    {
        const gsVector<index_t> a = mRagged.asVector(c);
        const gsVector<index_t> b = mOld.asVector(c);
        CHECK(a.size() > 0);
        CHECK_EQUAL(a.size(), b.size());
        CHECK_ARRAY_EQUAL(a.data(), b.data(), a.size());
    }
    CHECK_EQUAL(mOld.freeSize()    , mRagged.freeSize()    );
    CHECK_EQUAL(mOld.boundarySize(), mRagged.boundarySize());
    CHECK_EQUAL(mOld.size()        , mRagged.size()        );
}

// 15. matchDofs asserts equal-size index blocks (gsDofMapper.cpp:173): a
//     genuinely ragged interface (one component's patch refined, the other
//     not) must be built with conforming=false, or finalize() aborts.
TEST(ragged_matchDofs_requires_conforming_false)
{
    gsMultiBasis<real_t> mb0 = twoPatchBasis();
    gsMultiBasis<real_t> mb1 = twoPatchBasis();
    mb1.basis(0).uniformRefine(); // patch 0 of component 1 no longer matches
                                   // patch 1 of component 1 across the interface

    std::vector<const gsFunctionSet<real_t>*> basesPerComp(2);
    basesPerComp[0] = &mb0;
    basesPerComp[1] = &mb1;
    gsBoundaryConditions<real_t> bc;

    CHECK_THROW_IN_DEBUG(createMapper(basesPerComp, mb0.topology(), bc, /*unk=*/-1,
                             /*conforming=*/true, /*finalize=*/true),
               std::logic_error);

    // The same fixture with conforming=false does not throw.
    gsDofMapper m = createMapper(basesPerComp, mb0.topology(), bc, /*unk=*/-1,
                                 /*conforming=*/false, /*finalize=*/true);
    CHECK(m.isFinalized());
}

// 16. Out-of-range BC component indices are silently dropped by the ragged
//     overloads (documented, current behaviour -- not "fixed" here).
TEST(ragged_out_of_range_bc_component_dropped)
{
    gsMultiBasis<real_t> mb0 = twoPatchBasis();
    gsMultiBasis<real_t> mb1 = twoPatchBasis();

    gsFunctionExpr<real_t> g("0", 2);
    gsBoundaryConditions<real_t> bc;
    // Only 2 components (indices 0,1) are built below; component 5 never matches.
    bc.addCondition(0, boundary::west, condition_type::dirichlet, &g, 0, false, 5);

    std::vector<const gsFunctionSet<real_t>*> basesPerComp(2);
    basesPerComp[0] = &mb0;
    basesPerComp[1] = &mb1;

    gsDofMapper m = createMapper(basesPerComp, mb0.topology(), bc, /*unk=*/-1,
                                 /*conforming=*/false, /*finalize=*/true);

    CHECK_EQUAL(0, m.boundarySize());

    // Positive control: the identical condition with an in-range component
    // (0) IS applied, so the "dropped" assertion above is not vacuously true
    // on a createMapper that drops every boundary condition.
    const index_t nbWest = mb0.basis(0).boundary(boundary::west).rows();
    gsBoundaryConditions<real_t> bcInRange;
    bcInRange.addCondition(0, boundary::west, condition_type::dirichlet, &g, 0, false, 0);
    gsDofMapper mInRange = createMapper(basesPerComp, mb0.topology(), bcInRange, /*unk=*/-1,
                                        /*conforming=*/false, /*finalize=*/true);
    CHECK_EQUAL(nbWest, mInRange.boundarySize());
}

// 11. Primary 7-arg overload called directly with all arguments.
TEST(primary_7arg_overload)
{
    gsMultiBasis<real_t> mb = twoPatchBasis();

    gsFunctionExpr<real_t> g("0", 2);
    gsBoundaryConditions<real_t> bc;
    bc.addCondition(0, boundary::west, condition_type::dirichlet, &g);

    gsDofMapper m = createMapper(mb, mb.topology(), bc, /*nComp=*/1, /*unk=*/0,
                                 /*conforming=*/true, /*finalize=*/true);

    const index_t nMatched = countMatchedPairs(mb, mb.topology());
    const index_t nb       = mb.basis(0).boundary(boundary::west).rows();
    CHECK(nMatched > 0);
    CHECK(nb       > 0);
    CHECK_EQUAL(static_cast<index_t>(mb.totalSize()) - nMatched, m.size());
    CHECK_EQUAL(nb, m.boundarySize());
    CHECK_EQUAL(m.size() - m.boundarySize(), m.freeSize());
    CHECK(m.coupledSize() > 0);
}

// Regression guard for the defect this refactor fixed.
//
// gsMultiBasis::getMapper(bool conforming, index_t nComp, ...) used to discard its
// `conforming` argument and always build gsDofMapper(*this, topology(), nComp), i.e.
// always glued. So a caller asking for iFace::dg silently got conforming interfaces.
// createMapper must honour conforming = (is == iFace::glue) in BOTH branches.
//
// Without this test the suite cannot distinguish the fix from the bug: every other
// strategy-overload assertion uses iFace::glue, where both behaviours agree.
TEST(strategy_overload_honours_interface_strategy)
{
    gsMultiBasis<real_t> mb = twoPatchBasis();

    gsFunctionExpr<real_t> g("0", 2);
    gsBoundaryConditions<real_t> bc;
    bc.addCondition(0, boundary::west, condition_type::dirichlet, &g);

    // Non-elimination branch (the one that carried the bug): dg must NOT glue.
    gsDofMapper dgFree   = createMapper(mb, bc, dirichlet::nitsche, iFace::dg  , 1, 0, true);
    gsDofMapper glueFree = createMapper(mb, bc, dirichlet::nitsche, iFace::glue, 1, 0, true);

    CHECK_EQUAL(static_cast<index_t>(mb.totalSize()), dgFree.size());
    CHECK(glueFree.size() < dgFree.size());          // glue really does identify dofs
    CHECK_EQUAL(0, dgFree.coupledSize());            // and dg really does not

    // Elimination branch: dg must not glue either, but the Dirichlet dofs still go.
    gsDofMapper dgElim   = createMapper(mb, bc, dirichlet::elimination, iFace::dg  , 1, 0, true);
    gsDofMapper glueElim = createMapper(mb, bc, dirichlet::elimination, iFace::glue, 1, 0, true);

    CHECK_EQUAL(glueElim.boundarySize(), dgElim.boundarySize());
    CHECK(glueElim.size() < dgElim.size());
    CHECK_EQUAL(0, dgElim.coupledSize());
}

}
