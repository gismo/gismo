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

// GISMO_ENSURE reports its reason on std::cerr and throws a bare
// std::runtime_error, so a test that cares why a call was rejected reads
// the stream.
class CerrCapture
{
public:
    CerrCapture() : m_old(std::cerr.rdbuf(m_buf.rdbuf())) { }
    ~CerrCapture() { std::cerr.rdbuf(m_old); }
    bool contains(const std::string & s) const
    { return std::string::npos != m_buf.str().find(s); }
private:
    std::ostringstream m_buf;
    std::streambuf * m_old;
};

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


// --- per-component fixtures ------------------------------------------------
//
// The two components of a 2D Raviart-Thomas pair, degrees and regularities
// transposed between them:
//
//     component 0 = S^{3,2}_{2,1}      component 1 = S^{2,3}_{1,2}
//
// on an e0 x e1 mesh per patch.  With p-r == 1 in every direction the
// component dimensions per patch are (e0+3)(e1+2) and (e0+2)(e1+3): equal on
// an isotropic mesh (4x4: 42 and 42), unequal on an anisotropic one (4x2: 28
// and 30).  Their traces on an east/west side differ in either case (e1+2
// against e1+3), so the interface matching tells the components apart even
// where their sizes coincide.

gsTensorBSplineBasis<2,real_t> * rtComponentBasis(short_t p0, short_t p1,
                                                  index_t e0, index_t e1)
{
    // Maximal regularity (interior multiplicity p-r == 1) in both directions.
    gsKnotVector<real_t> kv0(0.0, 1.0, e0-1, p0+1, 1, p0);
    gsKnotVector<real_t> kv1(0.0, 1.0, e1-1, p1+1, 1, p1);
    return new gsTensorBSplineBasis<2,real_t>(kv0, kv1);
}

// One RT component on two unit squares side by side (patch 0 east == patch 1
// west), the same component basis on both patches.
gsMultiBasis<real_t> rtComponent(index_t c, index_t e0, index_t e1)
{
    gsMultiPatch<real_t> mp = gsNurbsCreator<real_t>::BSplineSquareGrid(2, 1, 1.0);
    const short_t p0 = (0 == c ? 3 : 2), p1 = (0 == c ? 2 : 3);
    gsMultiBasis<real_t>::BasisContainer bases;
    bases.push_back(rtComponentBasis(p0, p1, e0, e1));
    bases.push_back(rtComponentBasis(p0, p1, e0, e1));
    return gsMultiBasis<real_t>(bases, mp.topology());
}

std::vector<gsMultiBasis<real_t> > rtPair(index_t e0, index_t e1)
{
    std::vector<gsMultiBasis<real_t> > r;
    r.push_back(rtComponent(0, e0, e1));
    r.push_back(rtComponent(1, e0, e1));
    return r;
}

std::vector<const gsFunctionSet<real_t>*> pointers(const std::vector<gsMultiBasis<real_t> > & v)
{
    std::vector<const gsFunctionSet<real_t>*> r;
    for (size_t c = 0; c != v.size(); ++c)
        r.push_back(&v[c]);
    return r;
}

// One condition of every kind the creator applies, each on its own side or
// corner of the two-patch fixture.  A condition is added when \a keep says
// its component selection applies, and is then given the component \a as
// maps it to; this builds both the per-component input (every selection
// kept as it is) and, for the oracle, the conditions one component sees on
// its own (selections of that component and -1 kept, as component -1).
struct Selection
{
    // component selections of the Dirichlet(west), Dirichlet(north),
    // Clamped, Collapsed, Coupled and corner conditions, in that order;
    // omit leaves the condition out
    index_t sel[6];
};
const index_t omit = -99;

gsBoundaryConditions<real_t> allKinds(const Selection & s, index_t onlyComp = -2)
{
    static gsFunctionExpr<real_t> g("0", 2);
    gsBoundaryConditions<real_t> bc;
    // onlyComp == -2: keep every selection as it is.  Otherwise keep the
    // conditions that apply to component onlyComp, as broadcasts.
    struct Keep
    {
        static bool applies(index_t sel, index_t only)
        { return omit != sel && (-2 == only || -1 == sel || sel == only); }
        static index_t as(index_t sel, index_t only) { return -2 == only ? sel : -1; }
    };
    if (Keep::applies(s.sel[0], onlyComp))
        bc.addCondition(0, boundary::west , condition_type::dirichlet, &g, 0, false, Keep::as(s.sel[0], onlyComp));
    if (Keep::applies(s.sel[1], onlyComp))
        bc.addCondition(1, boundary::north, condition_type::dirichlet, &g, 0, false, Keep::as(s.sel[1], onlyComp));
    if (Keep::applies(s.sel[2], onlyComp))
        bc.addCondition(1, boundary::east , condition_type::clamped  , &g, 0, false, Keep::as(s.sel[2], onlyComp));
    if (Keep::applies(s.sel[3], onlyComp))
        bc.addCondition(1, boundary::south, condition_type::collapsed, &g, 0, false, Keep::as(s.sel[3], onlyComp));
    if (Keep::applies(s.sel[4], onlyComp))
        bc.addCoupled(0, boundary::south, 0, boundary::north, 2, 0, Keep::as(s.sel[4], onlyComp));
    if (Keep::applies(s.sel[5], onlyComp))
        bc.addCornerValue(boundary::southeast, 0.0, 0, 0, Keep::as(s.sel[5], onlyComp));
    return bc;
}

// Checks component c of the per-component mapper m against the one-component
// mapper ref built from component c's basis alone.  The numbering of a
// component depends only on the conditions applied to it, so, shifted into
// that component's free and eliminated blocks, the two must agree dof for
// dof -- which holds only if every interface, boundary and corner of
// component c was taken from component c's own basis.
void checkComponentAgainst(const gsDofMapper & m, index_t c, const gsDofMapper & ref)
{
    index_t freeBefore = 0, elimBefore = 0;
    for (index_t d = 0; d != c; ++d)
    {
        freeBefore += m.freeSize(d);
        elimBefore += m.size(d) - m.freeSize(d);
    }
    CHECK_EQUAL(ref.size()    , m.size(c));
    CHECK_EQUAL(ref.freeSize(), m.freeSize(c));
    CHECK_EQUAL(ref.numPatches(), m.numPatches());
    for (size_t k = 0; k != ref.numPatches(); ++k)
    {
        CHECK_EQUAL(ref.patchSize(k), m.patchSize(k, c));
        for (size_t i = 0; i != ref.patchSize(k); ++i)
        {
            const index_t r = ref.index(i, k);
            const index_t expected = ref.is_free_index(r)
                ? r + freeBefore
                : r - ref.freeSize() + m.freeSize() + elimBefore;
            CHECK_EQUAL(expected, m.index(i, k, c));
        }
    }
}

// Dof-for-dof equality of two finalized mappers, declared flag included.
void checkSameMapper(const gsDofMapper & a, const gsDofMapper & b)
{
    CHECK_EQUAL(a.numComponents(), b.numComponents());
    CHECK_EQUAL(a.numPatches()   , b.numPatches());
    CHECK_EQUAL(a.size()         , b.size());
    CHECK_EQUAL(a.freeSize()     , b.freeSize());
    CHECK_EQUAL(a.boundarySize() , b.boundarySize());
    CHECK_EQUAL(a.coupledSize()  , b.coupledSize());
    CHECK_EQUAL(a.hasDistinctComponentSpaces(), b.hasDistinctComponentSpaces());
    CHECK(a.layout() == b.layout());
    for (index_t c = 0; c != a.numComponents(); ++c)
        for (size_t k = 0; k != a.numPatches(); ++k)
        {
            CHECK_EQUAL(a.patchSize(k, c), b.patchSize(k, c));
            for (size_t i = 0; i != a.patchSize(k, c); ++i)
                CHECK_EQUAL(a.index(i, k, c), b.index(i, k, c));
        }
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

// 8b. Corner conditions with unknown == -1 are a wildcard: they must be
// eliminated regardless of the requested unk, while a corner condition with
// an explicit, different unknown must still be skipped.
TEST(corner_wildcard_unknown)
{
    gsMultiBasis<real_t> mb = twoPatchBasis();
    gsBoundaryConditions<real_t> bc;

    // Wildcard corner (unknown == -1): matches any requested unk.
    bc.addCornerValue(boundary::southwest, 0.5, 0, -1);
    // Explicit unknown == 1 corner: must be skipped when unk == 0.
    bc.addCornerValue(boundary::northwest, 0.5, 0, 1);

    gsDofMapper m0 = createMapper(mb, bc, 1, 0, false, true);
    CHECK_EQUAL(1, m0.boundarySize()); // only the wildcard corner

    gsDofMapper m1 = createMapper(mb, bc, 1, 1, false, true);
    CHECK_EQUAL(2, m1.boundarySize()); // wildcard + explicit unk==1 corner
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


// =========================================================================
// Per-component bases
// =========================================================================

// The RT pair on an ISOTROPIC mesh: both components have 42 dofs per patch,
// yet the mapper is declared distinct and rejected by the uniform
// evaluator, and each component is numbered from its own basis -- its own
// interface trace (6 against 7 dofs) included.
TEST(per_component_rt_isotropic)
{
    const std::vector<gsMultiBasis<real_t> > rt = rtPair(4, 4);
    CHECK_EQUAL(42, rt[0].basis(0).size());
    CHECK_EQUAL(42, rt[1].basis(0).size());

    const gsDofMapper m = createMapper(rt, gsBoundaryConditions<real_t>(), 0, true, true);

    CHECK(m.hasDistinctComponentSpaces());
    CHECK(!m.hasUniformComponents());
    CHECK_EQUAL(2, m.numComponents());
    CHECK_EQUAL(countMatchedPairs(rt[0], rt[0].topology()), 6);
    CHECK_EQUAL(countMatchedPairs(rt[1], rt[1].topology()), 7);
    CHECK_EQUAL(2*42 - 6, m.size(0));
    CHECK_EQUAL(2*42 - 7, m.size(1));
    for (index_t c = 0; c != 2; ++c)
        checkComponentAgainst(m, c, createMapper(rt[c], 1, true, true));
}

// The RT pair on an ANISOTROPIC mesh: 28 and 30 dofs per patch.
TEST(per_component_rt_anisotropic)
{
    const std::vector<gsMultiBasis<real_t> > rt = rtPair(4, 2);
    CHECK_EQUAL(28, rt[0].basis(0).size());
    CHECK_EQUAL(30, rt[1].basis(0).size());

    const gsDofMapper m = createMapper(pointers(rt), rt[0].topology(),
                                       gsBoundaryConditions<real_t>(), 0, true, true);

    CHECK(m.hasDistinctComponentSpaces());
    CHECK(!m.hasUniformComponents());
    CHECK_EQUAL(28u, m.patchSize(1, 0));
    CHECK_EQUAL(30u, m.patchSize(1, 1));
    for (index_t c = 0; c != 2; ++c)
        checkComponentAgainst(m, c, createMapper(rt[c], 1, true, true));
}

// Every condition kind, with component-specific and broadcast selections,
// against unequal and against equal-sized component bases.  Each component
// must see exactly the conditions selecting it or -1, applied to its own
// basis.
TEST(per_component_every_condition_kind)
{
    const Selection selections[] =
    {
        {{ 1, -1,  0,  1,  0,  1}},   // mixed specific and broadcast
        {{-1, -1, -1, -1, -1, -1}},   // everything broadcast
        {{ 0,  1,  1,  0,  1,  0}},   // everything specific, swapped
    };
    const index_t meshes[2][2] = {{4, 4}, {4, 2}};

    for (index_t mesh = 0; mesh != 2; ++mesh)
    {
        const std::vector<gsMultiBasis<real_t> > rt = rtPair(meshes[mesh][0], meshes[mesh][1]);
        for (size_t s = 0; s != sizeof(selections)/sizeof(selections[0]); ++s)
        {
            const gsDofMapper m = createMapper(rt, allKinds(selections[s]), 0, true, true);
            CHECK(m.hasDistinctComponentSpaces());
            CHECK(m.boundarySize() > 0);
            CHECK(m.coupledSize() > 0);
            for (index_t c = 0; c != 2; ++c)
                checkComponentAgainst(m, c,
                    createMapper(rt[c], allKinds(selections[s], c), 1, 0, true, true));
        }
    }
}

// One function-set object for every component is delegated to the
// single-basis creator: the same mapper, not declared distinct.
TEST(per_component_same_object_delegates)
{
    gsMultiBasis<real_t> mb = twoPatchBasis();
    std::vector<const gsFunctionSet<real_t>*> same(3, &mb);

    gsFunctionExpr<real_t> g("0", 2);
    gsBoundaryConditions<real_t> bc;
    bc.addCondition(0, boundary::west, condition_type::dirichlet, &g, 0, false, 1);
    bc.addCondition(1, boundary::east, condition_type::dirichlet, &g, 0, false, -1);
    bc.addCornerValue(boundary::northwest, 0.0, 0, 0, 2);

    for (index_t fin = 0; fin != 2; ++fin)
    {
        gsDofMapper a = createMapper(same, mb.topology(), bc, 0, true, 1 == fin);
        gsDofMapper b = createMapper(mb  , mb.topology(), bc, 3, 0, true, 1 == fin);
        CHECK_EQUAL(1 == fin, a.isFinalized());
        if (0 == fin) { a.finalize(); b.finalize(); }
        CHECK(!a.hasDistinctComponentSpaces());
        CHECK(a.hasUniformComponents());
        checkSameMapper(a, b);
    }

    // one component is the single-basis case too
    gsBoundaryConditions<real_t> bc1;
    bc1.addCondition(0, boundary::west, condition_type::dirichlet, &g, 0, false, -1);
    bc1.addCornerValue(boundary::northeast, 0.0, 1, 0, 0);
    std::vector<const gsFunctionSet<real_t>*> one(1, &mb);
    checkSameMapper(createMapper(one, mb.topology(), bc1, 0, true, true),
                    createMapper(mb , mb.topology(), bc1, 1, 0, true, true));
}

// Distinct but equal function-set objects are distinct component spaces:
// the flag records how the mapper was built, not what the sizes say.  The
// numbering is still that of the single-basis creator.
TEST(per_component_equal_copies_are_declared_distinct)
{
    gsMultiBasis<real_t> mb = twoPatchBasis();
    std::vector<gsMultiBasis<real_t> > copies(2, mb);

    const gsDofMapper a = createMapper(copies, gsBoundaryConditions<real_t>(), 0, true, true);
    const gsDofMapper b = createMapper(mb, 2, true, true);

    CHECK(a.hasDistinctComponentSpaces());
    CHECK(!a.hasUniformComponents());
    CHECK(!b.hasDistinctComponentSpaces());
    for (index_t c = 0; c != 2; ++c)
        checkComponentAgainst(a, c, createMapper(mb, 1, true, true));
    CHECK_EQUAL(b.size(), a.size());
}

// A condition selecting a component that does not exist is rejected, for
// every kind, in every build type, and whether or not the input is
// delegated.  One for another unknown is not looked at.
TEST(per_component_rejects_invalid_component_selection)
{
    const std::vector<gsMultiBasis<real_t> > rt = rtPair(4, 2);
    gsMultiBasis<real_t> mb = twoPatchBasis();
    std::vector<const gsFunctionSet<real_t>*> same(2, &mb);

    const index_t bad[] = {2, 5, -2};
    for (size_t b = 0; b != sizeof(bad)/sizeof(bad[0]); ++b)
        for (index_t kind = 0; kind != 6; ++kind)
        {
            Selection s = {{-1, -1, -1, -1, -1, -1}};
            s.sel[kind] = bad[b];
            const gsBoundaryConditions<real_t> bc = allKinds(s);
            CHECK_THROW(createMapper(rt, bc, 0, true, true), std::runtime_error);
            CHECK_THROW(createMapper(same, mb.topology(), bc, 0, true, true), std::runtime_error);
            // filtered out by the unknown: accepted
            createMapper(rt, bc, 1, true, true);
        }
}

// The single-basis creator keeps ignoring a Clamped, Collapsed, coupled or
// corner condition whose component selects none of its components, as it
// always has, while a Dirichlet one is an error.
TEST(single_basis_keeps_ignoring_unmatched_component_selection)
{
    gsMultiBasis<real_t> mb = twoPatchBasis();
    const gsDofMapper plain = createMapper(mb, 2, true, true);

    for (index_t kind = 0; kind != 6; ++kind)
    {
        Selection s = {{omit, omit, omit, omit, omit, omit}};
        s.sel[kind] = 5;
        const gsBoundaryConditions<real_t> bc = allKinds(s);
        if (kind < 2)
            CHECK_THROW(createMapper(mb, bc, 2, 0, true, true), std::runtime_error);
        else
            checkSameMapper(createMapper(mb, bc, 2, 0, true, true), plain);
    }
}

// Mapped bases have no per-component meaning here, whether every component
// or only some are mapped -- also for one component, which would otherwise
// be delegated -- and in every dimension a gsMappedBasis is instantiated
// for.  The single-basis mapped path is unchanged (mapped_basis_with_bcs).
TEST(per_component_rejects_mapped_bases)
{
    auto geom = gsNurbsCreator<real_t>::BSplineSquare(2);
    gsMultiBasis<real_t> mb(geom->basis());
    mb.basis(0).uniformRefine(2);
    const index_t sz = mb.basis(0).size();
    gsSparseMatrix<real_t> ident(sz, sz);
    ident.setIdentity();
    gsMappedBasis<2, real_t> mapB(mb, ident);

    std::vector<const gsFunctionSet<real_t>*> allMapped(2, &mapB), mixed, one(1, &mapB);
    mixed.push_back(&mb);
    mixed.push_back(&mapB);

    const gsBoundaryConditions<real_t> none;
    CHECK_THROW(createMapper(allMapped, mb.topology(), none, 0, false, true), std::runtime_error);
    CHECK_THROW(createMapper(mixed    , mb.topology(), none, 0, false, true), std::runtime_error);
    CHECK_THROW(createMapper(one      , mb.topology(), none, 0, false, true), std::runtime_error);

    // 1D: a line, not a square
    gsKnotVector<real_t> kv(0.0, 1.0, 2, 3);
    const gsBSplineBasis<real_t> line(kv);
    gsMultiBasis<real_t> mb1(line);
    const index_t sz1 = mb1.basis(0).size();
    gsSparseMatrix<real_t> ident1(sz1, sz1);
    ident1.setIdentity();
    gsMappedBasis<1, real_t> mapB1(mb1, ident1);

    std::vector<const gsFunctionSet<real_t>*> one1(1, &mapB1), mixed1;
    mixed1.push_back(&mb1);
    mixed1.push_back(&mapB1);
    CHECK_THROW(createMapper(one1  , mb1.topology(), none, 0, false, true), std::runtime_error);
    CHECK_THROW(createMapper(mixed1, mb1.topology(), none, 0, false, true), std::runtime_error);
}

// Function sets that cannot describe the components of one variable.
TEST(per_component_rejects_incompatible_function_sets)
{
    const gsBoundaryConditions<real_t> none;
    gsMultiBasis<real_t> two = twoPatchBasis();
    gsMultiBasis<real_t> one(two.basis(0));

    std::vector<const gsFunctionSet<real_t>*> v;
    CHECK_THROW(createMapper(v, two.topology(), none, 0, true, true), std::runtime_error);
    CHECK_THROW(createMapper(std::vector<gsMultiBasis<real_t> >(), none, 0, true, true),
                std::runtime_error);

    v.push_back(&two);
    v.push_back(nullptr);
    CHECK_THROW(createMapper(v, two.topology(), none, 0, true, true), std::runtime_error);

    v[1] = &one;    // 2 patches against 1
    CHECK_THROW(createMapper(v, two.topology(), none, 0, true, true), std::runtime_error);

    gsMultiBasis<real_t> empty;
    v[0] = &empty;  // no patches at all
    v[1] = &empty;
    CHECK_THROW(createMapper(v, gsBoxTopology(), none, 0, true, true), std::runtime_error);

    // same patch count, different domain dimension
    auto cubeGeo = gsNurbsCreator<real_t>::BSplineCube(1);
    gsMultiBasis<real_t> cube(cubeGeo->basis());
    v[0] = &one;
    v[1] = &cube;
    CHECK_THROW(createMapper(v, gsBoxTopology(), none, 0, false, true), std::runtime_error);
}

// The gsMultiBasis overload takes the topology from component 0, so every
// component must have the same one.  Interface orientation and labels do
// not matter; a missing interface, a different interaction type or a
// different boundary does.
TEST(per_component_topologies_must_agree)
{
    const gsBoundaryConditions<real_t> none;
    const std::vector<gsMultiBasis<real_t> > rt = rtPair(4, 2);
    const gsBoxTopology & topo = rt[0].topology();
    const boundaryInterface bi = *topo.iBegin();

    std::vector<gsMultiBasis<real_t> > v = rt;
    gsBoxTopology t = topo;
    // the same interface stored the other way round
    t.interfaces()[0] = bi.getInverse();
    v[1].setTopology(t);
    checkSameMapper(createMapper(v, none, 0, true, true),
                    createMapper(pointers(rt), topo, none, 0, true, true));

    t = topo;
    t.interfaces().clear();
    v[1].setTopology(t);
    CHECK_THROW(createMapper(v, none, 0, true, true), std::runtime_error);

    t = topo;
    t.interfaces()[0].setAsContact();
    v[1].setTopology(t);
    CHECK_THROW(createMapper(v, none, 0, true, true), std::runtime_error);

    t = topo;
    t.boundaries().pop_back();
    v[1].setTopology(t);
    CHECK_THROW(createMapper(v, none, 0, true, true), std::runtime_error);
}

// References to patches, sides or corners that do not exist are rejected in
// every build type, by both creators.
TEST(creators_reject_missing_patches)
{
    gsFunctionExpr<real_t> g("0", 2);
    const std::vector<gsMultiBasis<real_t> > rt = rtPair(4, 2);
    gsMultiBasis<real_t> mb = twoPatchBasis();

    gsBoundaryConditions<real_t> bcPatch, bcCorner, bcCoupled, bcCornerIndex, bcCornerZero;
    bcPatch  .addCondition(2, boundary::west, condition_type::dirichlet, &g);
    bcCorner .addCornerValue(boundary::southwest, 0.0, 7, 0);
    bcCoupled.addCoupled(0, boundary::west, 3, boundary::east, 2, 0);
    // a 2D patch has the corners 1..4
    bcCornerIndex.addCornerValue(boxCorner(5), 0.0, 0, 0);
    bcCornerZero .addCornerValue(boxCorner(0), 0.0, 0, 0);
    const gsBoundaryConditions<real_t> * bcs[5] =
        {&bcPatch, &bcCorner, &bcCoupled, &bcCornerIndex, &bcCornerZero};
    for (index_t i = 0; i != 5; ++i)
    {
        CHECK_THROW(createMapper(rt, *bcs[i], 0, true, true), std::runtime_error);
        CHECK_THROW(createMapper(mb, *bcs[i], 1, 0, true, true), std::runtime_error);
    }

    // an interface on a patch the function sets do not have
    gsBoxTopology topo(2, 3);
    topo.addInterface(0, boundary::east, 2, boundary::west);
    const gsBoundaryConditions<real_t> none;
    CHECK_THROW(createMapper(pointers(rt), topo, none, 0, true, true), std::runtime_error);
    CHECK_THROW(createMapper(mb, topo, 1, true, true), std::runtime_error);
    // not visited without the conforming loop
    createMapper(pointers(rt), topo, none, 0, false, true);
    createMapper(mb, topo, 1, false, true);
}

// A coupled condition between sides with different numbers of dofs in some
// component cannot be applied, and says so.
TEST(per_component_coupled_sides_must_match)
{
    const std::vector<gsMultiBasis<real_t> > rt = rtPair(4, 2);
    // component 0 has e1+2 = 4 dofs on a west side, e0+3 = 7 on a south one
    CHECK_EQUAL(4, rt[0].basis(0).boundary(boundary::west ).rows());
    CHECK_EQUAL(7, rt[0].basis(0).boundary(boundary::south).rows());

    gsBoundaryConditions<real_t> bc;
    bc.addCoupled(0, boundary::west, 0, boundary::south, 2, 0, 0);
    CHECK_THROW(createMapper(rt, bc, 0, true, true), std::runtime_error);
}

// Distinct per-component bases are matched component by component, which is
// undefined across an interface that permutes the parametric directions:
// with patch 1 attached by its south side to patch 0's east side, the
// component normal to the interface on patch 0 is tangential on patch 1.
// Such an interface is rejected by name before any trace is matched.
TEST(per_component_rejects_interfaces_permuting_directions)
{
    const std::vector<gsMultiBasis<real_t> > rt = rtPair(3, 4);
    gsBoxTopology rotated(2, 2);
    rotated.addInterface(0, boundary::east, 1, boundary::south);
    const gsBoundaryConditions<real_t> none;
    const char * const reason = "requires every interface to keep each direction";

    // Two distinct copies of one basis: under the rotation every
    // component's traces are equally long (e1+2 == e0+3 == 6), so nothing
    // but the direction map tells the case apart, and without the check the
    // interface would be glued without any error.
    const std::vector<gsMultiBasis<real_t> > copies(2, rt[0]);
    CHECK_EQUAL(copies[0].basis(0).boundary(boundary::east ).rows(),
                copies[1].basis(1).boundary(boundary::south).rows());
    {
        CerrCapture err;
        CHECK_THROW(createMapper(pointers(copies), rotated, none, 0, true, true),
                    std::runtime_error);
        CHECK(err.contains(reason));
    }

    // The Raviart-Thomas pair is rejected for the same reason, not for its
    // component 1 traces differing in size (e1+3 against e0+2).
    {
        CerrCapture err;
        CHECK_THROW(createMapper(pointers(rt), rotated, none, 0, true, true),
                    std::runtime_error);
        CHECK(err.contains(reason));
    }

    // Not visited without the conforming loop.
    createMapper(pointers(rt), rotated, none, 0, false, true);

    // One shared function set is matched as by the single-basis overload.
    const std::vector<const gsFunctionSet<real_t>*> shared(2, &rt[0]);
    const gsDofMapper m = createMapper(shared, rotated, none, 0, true, true);
    CHECK(!m.hasDistinctComponentSpaces());
    CHECK(m.coupledSize() > 0);
    checkSameMapper(m, createMapper(rt[0], rotated, none, 2, 0, true, true));
}

// Contact interfaces are not glued, and the finalize flag is honoured, as in
// the single-basis creator.
TEST(per_component_contact_and_finalize)
{
    const std::vector<gsMultiBasis<real_t> > rt = rtPair(4, 4);
    const gsBoundaryConditions<real_t> none;

    gsBoxTopology contact = rt[0].topology();
    contact.interfaces()[0].setAsContact();
    const gsDofMapper m = createMapper(pointers(rt), contact, none, 0, true, true);
    CHECK_EQUAL(0, m.coupledSize());
    CHECK_EQUAL(4*42, m.size());

    gsDofMapper u = createMapper(pointers(rt), rt[0].topology(), none, 0, true, false);
    CHECK(!u.isFinalized());
    CHECK(u.hasDistinctComponentSpaces());
    u.finalize();
    checkSameMapper(u, createMapper(rt, none, 0, true, true));
}

// Conditions are filtered by unknown as in the single-basis creator:
// side conditions by their own unknown, coupled and corner conditions with
// unknown -1 for every unknown.
TEST(per_component_unknown_filter)
{
    gsFunctionExpr<real_t> g("0", 2);
    const std::vector<gsMultiBasis<real_t> > rt = rtPair(4, 2);
    const gsDofMapper plain = createMapper(rt, gsBoundaryConditions<real_t>(), 0, true, true);

    gsBoundaryConditions<real_t> bc;
    bc.addCondition(0, boundary::west, condition_type::dirichlet, &g, /*unknown=*/1, false, 1);
    bc.addCornerValue(boundary::northeast, 0.0, 1, /*unknown=*/-1, 0);

    const gsDofMapper m0 = createMapper(rt, bc, 0, true, true);
    CHECK_EQUAL(1, m0.boundarySize());                    // the wildcard corner only
    const gsDofMapper m1 = createMapper(rt, bc, 1, true, true);
    CHECK_EQUAL(1 + rt[1].basis(0).boundary(boundary::west).rows(), m1.boundarySize());
    const gsDofMapper mAll = createMapper(rt, bc, -1, true, true);
    CHECK_EQUAL(m1.boundarySize(), mAll.boundarySize());
    CHECK_EQUAL(plain.size(), m1.size());
}

}
