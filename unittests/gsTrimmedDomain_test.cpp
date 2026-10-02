/** @file gsTrimmedDomain_test.cpp

    @brief Behaviour lock for gsTrimmedDomain's cut-cell classification.

    PURPOSE. gsTrimmedDomain classifies every leaf of its kd-tree as interior
    (-1), exterior (+1) or cut (0) by sampling the level set on a Lobatto grid.
    Four init() overloads feed that machinery, and they do NOT share a code
    path: init(tbasis,.) and init(bbox,.) run _classifyTree(), while
    init(htbasis,.) and the size-based init() run _classifyTreeAdaptive() with
    different stopping predicates. Any change to the classification loops --
    reordering it, skipping work, or parallelising it -- must leave all four
    producing the identical tree and the identical signs.

    WHAT IS AND IS NOT EVIDENCE HERE. Element/leaf COUNTS alone are weak: a
    change that swapped two leaves' signs, or that classified the right number
    of leaves at the wrong subdivision level, reproduces every count exactly.
    So each test pins a FINGERPRINT -- an order-sensitive FNV-1a hash over
    (sign, level, lowerCorner, upperCorner) of every leaf in leaf-iterator
    order. That is the real invariant; it is built only from integers, so it
    carries no floating-point portability risk. The counts are asserted
    alongside it purely to LOCALIZE a failure: if a count moves you know which
    sign class broke, whereas the hash alone only says "something did".

    The expected values were captured from the pre-refactor implementation and
    must not be "refreshed" to match new behaviour: that would invert the
    purpose of the file. A deliberate behaviour change should edit the values
    in the same commit that justifies it in the message.

    The level set is a genuine signed distance (exactly 1-Lipschitz), not an
    implicit surface with arbitrary gradient scaling, so these fixtures stay
    valid for any future Lipschitz-based classification shortcut.

    This suite deliberately depends on no optional module.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst
 **/

#include "gismo_unittest.h"       // Brings in G+Smo and the UnitTest++ framework

#include <gsDomain/gsImplicitTrimmedDomain.h>

namespace {

/// Signed distance to a sphere, with a FINITE support() so that the
/// size-based init() overload (which reads support() to get its root box) is
/// exercisable -- gsFunctionExpr reports unbounded support and cannot drive
/// that path at all.
template<short_t d>
class gsSphereDist : public gsFunction<real_t>
{
public:
    GISMO_CLONE_FUNCTION(gsSphereDist)

    gsSphereDist(real_t cx, real_t radius) : m_c(cx), m_r(radius) {}

    short_t domainDim() const override { return d; }
    short_t targetDim() const override { return 1; }

    gsMatrix<real_t> support() const override
    {
        gsMatrix<real_t> s(d, 2);
        s.col(0).setZero();
        s.col(1).setOnes();
        return s;
    }

    void eval_into(const gsMatrix<real_t> & u, gsMatrix<real_t> & result) const override
    {
        result.resize(1, u.cols());
        for (index_t k = 0; k < u.cols(); ++k)
        {
            real_t s = 0;
            for (short_t j = 0; j < d; ++j)
                s += (u(j,k) - m_c) * (u(j,k) - m_c);
            result(0,k) = math::sqrt(s) - m_r;
        }
    }

private:
    real_t m_c, m_r;
};

/// Order-sensitive fingerprint of the classified tree: FNV-1a over
/// (sign, level, lowerCorner, upperCorner) of every leaf, in leaf-iterator
/// order. Integer-only, so it is bit-reproducible across platforms.
template<short_t d, class Domain>
unsigned long long treeFingerprint(const Domain & dom, size_t & nLeaves)
{
    unsigned long long h = 14695981039346656037ull;
    const auto mix = [&h](long long v)
    {
        h ^= static_cast<unsigned long long>(v);
        h *= 1099511628211ull;
    };

    auto it = dom.tree().beginLeafIterator();
    nLeaves = 0;
    while (it.good())
    {
        mix(it.data().sign());
        mix(it.data().level());
        for (short_t j = 0; j < d; ++j)
        {
            mix(static_cast<long long>(it.data().lowerCorner()[j]));
            mix(static_cast<long long>(it.data().upperCorner()[j]));
        }
        ++nLeaves;
        it.next();
    }
    return h;
}

/// Every leaf must carry one of the three legal signs, and the sign classes
/// must partition the element count exactly. Note AllSign is `s <= 0`, i.e.
/// the ACTIVE cells (interior + cut) -- not the whole tree; AnySign is the
/// whole tree (gsTrimmedDomain.h:25-33). Getting these two confused is the
/// easy mistake, so both identities are asserted.
///
/// This runs on every path and is independent of the pinned constants below,
/// so it keeps teeth even if someone does refresh them.
template<class Domain>
void checkPartition(const Domain & dom)
{
    auto it = dom.tree().beginLeafIterator();
    while (it.good())
    {
        const short_t s = it.data().sign();
        CHECK(s == -1 || s == 0 || s == 1);
        it.next();
    }
    const size_t interior = dom.template numElements<InteriorSign>();
    const size_t boundary = dom.template numElements<BoundarySign>();
    const size_t exterior = dom.template numElements<ExteriorSign>();

    CHECK_EQUAL(dom.template numElements<AllSign>(), interior + boundary);
    CHECK_EQUAL(dom.template numElements<AnySign>(), interior + boundary + exterior);
}

} // anonymous namespace

SUITE(gsTrimmedDomain_test)       // The suite should have the same name as the file
{

/// PATH 1 -- init(tbasis, samples) -> _classifyTree().
///
/// This is the path the immersed FCM examples take. The tree starts as one
/// leaf spanning every level-0 element and is bisected until each leaf holds a
/// single element; classification happens along the way.
///
/// A 3D 8^3 grid with a sphere of radius 0.3 centred in the unit cube gives a
/// genuinely three-way split (interior, cut shell, exterior bulk) rather than
/// the degenerate all-cut or all-exterior cases a coarse grid can produce.
TEST(TensorBasisPath)
{
    gsSphereDist<3> phi(0.5, 0.3);
    gsKnotVector<real_t> kv(0, 1, 7, 2);          // 8 elements per direction
    gsTensorBSplineBasis<3,real_t> tb(kv, kv, kv);

    gsImplicitTrimmedDomain<3,real_t> dom(phi, tb, 5);

    size_t nLeaves = 0;
    const unsigned long long h = treeFingerprint<3>(dom, nLeaves);

    checkPartition(dom);
    CHECK_EQUAL(512u, (unsigned)dom.numElements<AnySign>());
    CHECK_EQUAL(136u, (unsigned)dom.numElements<AllSign>());
    CHECK_EQUAL(8u,   (unsigned)dom.numElements<InteriorSign>());
    CHECK_EQUAL(128u, (unsigned)dom.numElements<BoundarySign>());
    CHECK_EQUAL(376u, (unsigned)dom.numElements<ExteriorSign>());
    CHECK_EQUAL(512u, (unsigned)nLeaves);
    CHECK_EQUAL(17270146533631227653ull, h);
}

/// PATH 2 -- init(bbox, numCells, samples, deg) -> _classifyTree().
///
/// Same classification loop as PATH 1 but reached without a basis, and on a
/// DELIBERATELY non-square box with a DELIBERATELY non-symmetric cell count:
/// a transposed direction index or a collapse to a single count cannot hide
/// behind symmetry here, whereas it could in the cubic PATH 1 fixture.
TEST(BoundingBoxPath)
{
    gsSphereDist<2> phi(0.5, 0.3);
    gsMatrix<real_t> bbox(2,2);
    bbox << 0.0, 0.0,        // row 0 = lower corners  (2 x d layout)
            1.0, 0.8;        // row 1 = upper corners
    gsVector<index_t,2> nc;
    nc << 6, 5;

    gsImplicitTrimmedDomain<2,real_t> dom(phi, bbox, nc, 5, 1);

    size_t nLeaves = 0;
    const unsigned long long h = treeFingerprint<2>(dom, nLeaves);

    checkPartition(dom);
    CHECK_EQUAL(30u, (unsigned)dom.numElements<AnySign>());
    CHECK_EQUAL(16u, (unsigned)dom.numElements<AllSign>());
    CHECK_EQUAL(30u, (unsigned)nLeaves);
    CHECK_EQUAL(4u,  (unsigned)dom.numElements<InteriorSign>());
    CHECK_EQUAL(12u, (unsigned)dom.numElements<BoundarySign>());
    CHECK_EQUAL(14u, (unsigned)dom.numElements<ExteriorSign>());
    CHECK_EQUAL(12097644341665813797ull, h);
}

/// PATH 3 -- init(htbasis, samples) -> _classifyTreeAdaptive(), always-stop.
///
/// The hierarchical structure is already fully specified by the basis, so the
/// stopping predicate is the constant `true` and NO refinement may occur: the
/// classified tree must mirror the THB tree leaf for leaf. The locally refined
/// box makes the tree genuinely multi-level, so a change that mishandled leaf
/// levels (e.g. reading break vectors at the wrong level) would move the
/// fingerprint even where the counts survived.
TEST(HTensorBasisPath)
{
    // Small, off-centre circle on an 8-element base with a refined lower-left
    // quadrant: this keeps some leaves entirely clear of the interface, so the
    // fixture exercises all three signs instead of collapsing to all-cut (a
    // centred circle on a coarse base classifies every large leaf as 0, which
    // would make the sign part of the fingerprint vacuous).
    gsSphereDist<2> phi(0.3, 0.12);
    gsKnotVector<real_t> kv(0, 1, 7, 2);          // 8 elements per direction
    gsTensorBSplineBasis<2,real_t> tb(kv, kv);
    gsTHBSplineBasis<2,real_t> thb(tb);

    std::vector<index_t> box;                      // {level, lo_0, lo_1, up_0, up_1}
    box.push_back(1);
    box.push_back(0); box.push_back(0);
    box.push_back(8); box.push_back(8);
    thb.refineElements(box);

    gsImplicitTrimmedDomain<2,real_t> dom(phi, thb, 5);

    size_t nLeaves = 0;
    const unsigned long long h = treeFingerprint<2>(dom, nLeaves);

    checkPartition(dom);
    // The mirrored THB tree is only 3 large kd-leaves, so no leaf fits strictly
    // inside a radius-0.12 disc and the interior class is legitimately empty
    // here; PATH 4 is the adaptive-path fixture that covers sign -1.
    CHECK_EQUAL(3u,   (unsigned)nLeaves);
    CHECK_EQUAL(112u, (unsigned)dom.numElements<AnySign>());
    CHECK_EQUAL(0u,   (unsigned)dom.numElements<InteriorSign>());
    CHECK_EQUAL(64u,  (unsigned)dom.numElements<BoundarySign>());
    CHECK_EQUAL(48u,  (unsigned)dom.numElements<ExteriorSign>());
    CHECK_EQUAL(5044972201710581758ull, h);
}

/// PATH 4 -- init(maxElementSize, minElementSize, samples, deg)
///           -> _classifyTreeAdaptive() with a REFINING predicate.
///
/// The only path where the classification loop actually mutates the tree while
/// walking it: leaves are split and promoted to the next dyadic level until
/// they meet the size criterion, with cut leaves driven to the finer
/// minElementSize. This is the path that constrains any restructuring of that
/// loop into rounds, so the multi-level fingerprint matters most here.
TEST(SizeBasedPath)
{
    // This ctor hardcodes init(maxElementSize=1, minElementSize=0.1,
    // samples=10, deg), so the fixture can only be steered through the level
    // set. A radius of 0.4 makes the disc large enough that whole refined
    // leaves land strictly inside it, so interior cells actually appear --
    // at 0.3 the interior class comes out empty and the sign part of the
    // fingerprint loses most of its discriminating power.
    gsSphereDist<2> phi(0.5, 0.4);
    gsImplicitTrimmedDomain<2,real_t> dom(phi, 1);

    size_t nLeaves = 0;
    const unsigned long long h = treeFingerprint<2>(dom, nLeaves);

    checkPartition(dom);
    CHECK(dom.numLevels() > 1);                    // refinement really happened
    // Non-cut leaves stop refining immediately (their volume already meets
    // maxElementSize=1), so at r=0.4 the disc leaves no leaf strictly outside
    // and the exterior class is legitimately empty here; PATH 3 covers sign +1
    // on the adaptive path.
    CHECK_EQUAL(16u,  (unsigned)nLeaves);
    CHECK_EQUAL(256u, (unsigned)dom.numElements<AnySign>());
    CHECK_EQUAL(64u,  (unsigned)dom.numElements<InteriorSign>());
    CHECK_EQUAL(192u, (unsigned)dom.numElements<BoundarySign>());
    CHECK_EQUAL(0u,   (unsigned)dom.numElements<ExteriorSign>());
    CHECK_EQUAL(10947057143385013373ull, h);
}

/// Classification must not depend on how many threads happen to run it.
/// Once the classification loops are parallelised this is the test that
/// catches a race or a schedule-dependent result; before that it is a cheap
/// tautology, which is the point -- it is in place BEFORE the change.
TEST(ThreadCountInvariance)
{
    gsSphereDist<3> phi(0.5, 0.3);
    gsKnotVector<real_t> kv(0, 1, 7, 2);
    gsTensorBSplineBasis<3,real_t> tb(kv, kv, kv);

    const int saved = omp_get_max_threads();

    omp_set_num_threads(1);
    size_t n1 = 0;
    gsImplicitTrimmedDomain<3,real_t> dom1(phi, tb, 5);
    const unsigned long long h1 = treeFingerprint<3>(dom1, n1);

    omp_set_num_threads(4);
    size_t n4 = 0;
    gsImplicitTrimmedDomain<3,real_t> dom4(phi, tb, 5);
    const unsigned long long h4 = treeFingerprint<3>(dom4, n4);

    omp_set_num_threads(saved);

    CHECK_EQUAL(n1, n4);
    CHECK_EQUAL(h1, h4);
}

}
