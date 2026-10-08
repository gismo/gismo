/** @file gsComputeTopology_test.cpp

    @brief Tests that the spatially binned gsMultiPatch::computeTopology
    reproduces the exhaustive pairwise search exactly.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "gismo_unittest.h"

#include <cmath>
#include <cstdint>
#include <vector>

using namespace gismo;

namespace {

// Re-exports the protected pairwise reference; never instantiated.
struct gsTopoProbe : public gsMultiPatch<real_t>
{
    using gsMultiPatch<real_t>::computeTopologyPairwise;
};

// &gsTopoProbe::f names the base member, so this is a member pointer of gsMultiPatch<real_t>.
const auto pairwise = &gsTopoProbe::computeTopologyPairwise;

// Fixed-seed 64-bit LCG (deterministic across platforms and stdlibs).
struct Lcg
{
    uint64_t s;
    explicit Lcg(uint64_t seed) : s(seed) {}
    // 31 high-quality bits
    uint64_t next()
    {
        s = s * 6364136223846793005ULL + 1442695040888963407ULL;
        return s >> 33;
    }
    // uniform in [-1,1)
    real_t sym() { return real_t(2) * real_t(next()) / real_t(1ULL << 31) - real_t(1); }
};

struct Snapshot
{
    std::vector<boundaryInterface> ifc;
    std::vector<patchSide>         bdr;
};

Snapshot snap(const gsMultiPatch<real_t> & mp)
{
    Snapshot s;
    s.ifc = mp.interfaces();
    s.bdr = mp.boundaries();
    return s;
}

bool sameSide(const patchSide & a, const patchSide & b)
{
    return a.patch == b.patch && int(a.side()) == int(b.side());
}

// Entry-wise comparison: same lists, same order, same sides, dirMap and
// dirOrientation. Prints the first mismatch.
bool sameTopology(const Snapshot & a, const Snapshot & b, bool verbose = true)
{
    if (a.ifc.size() != b.ifc.size())
    {
        if (verbose)
            gsInfo << "interface count " << a.ifc.size() << " vs " << b.ifc.size() << "\n";
        return false;
    }
    for (size_t i = 0; i != a.ifc.size(); ++i)
    {
        const boundaryInterface & x = a.ifc[i];
        const boundaryInterface & y = b.ifc[i];
        bool same = sameSide(x.first(), y.first()) && sameSide(x.second(), y.second())
            && x.dirMap().size() == y.dirMap().size()
            && x.dirOrientation().size() == y.dirOrientation().size();
        for (index_t k = 0; same && k != x.dirMap().size(); ++k)
            same = x.dirMap()(k) == y.dirMap()(k);
        for (index_t k = 0; same && k != x.dirOrientation().size(); ++k)
            same = x.dirOrientation()(k) == y.dirOrientation()(k);
        if (!same)
        {
            if (verbose)
                gsInfo << "interface " << i << " differs:\n  " << x << "\n  " << y << "\n";
            return false;
        }
    }
    if (a.bdr.size() != b.bdr.size())
    {
        if (verbose)
            gsInfo << "boundary count " << a.bdr.size() << " vs " << b.bdr.size() << "\n";
        return false;
    }
    for (size_t i = 0; i != a.bdr.size(); ++i)
        if (!sameSide(a.bdr[i], b.bdr[i]))
        {
            if (verbose)
                gsInfo << "boundary " << i << " differs:\n  " << a.bdr[i] << "\n  " << b.bdr[i] << "\n";
            return false;
        }
    return true;
}

bool counts(const Snapshot & s, size_t nIfc, size_t nBdr)
{
    const bool ok = s.ifc.size() == nIfc && s.bdr.size() == nBdr;
    if (!ok)
        gsInfo << "counts: interfaces " << s.ifc.size() << " (expected " << nIfc
               << "), boundaries " << s.bdr.size() << " (expected " << nBdr << ")\n";
    return ok;
}

// Runs the binned (public) and the exhaustive search on the same object.
// True iff both return true and produce identical topologies; ref receives
// the exhaustive result.
bool agree(gsMultiPatch<real_t> & mp, real_t tol, bool cornersOnly, Snapshot & ref)
{
    const bool r1 = mp.computeTopology(tol, cornersOnly);
    const Snapshot hashed = snap(mp);
    const bool r2 = (mp.*pairwise)(tol, cornersOnly);
    ref = snap(mp);
    return r1 && r2 && sameTopology(hashed, ref);
}

// agree() plus the derived counts of the exhaustive result.
bool agreeWithCounts(gsMultiPatch<real_t> & mp, real_t tol, bool cornersOnly,
                     size_t nIfc, size_t nBdr)
{
    Snapshot ref;
    const bool ok = agree(mp, tol, cornersOnly, ref);
    return counts(ref, nIfc, nBdr) && ok;
}

// agreeWithCounts() in both matching modes.
bool agreeBothModes(gsMultiPatch<real_t> & mp, real_t tol, size_t nIfc, size_t nBdr)
{
    const bool a = agreeWithCounts(mp, tol, false, nIfc, nBdr);
    const bool b = agreeWithCounts(mp, tol, true, nIfc, nBdr);
    return a && b;
}

// Whether some interface has a tangential dirMap entry that is not the
// identity (dm), and whether some has a reversed tangential direction (dor).
void orientationFlags(const Snapshot & s, bool & dm, bool & dor)
{
    dm = false;
    dor = false;
    for (size_t i = 0; i != s.ifc.size(); ++i)
    {
        const boundaryInterface & x = s.ifc[i];
        for (index_t k = 0; k != x.dirMap().size(); ++k)
        {
            if (x.dirMap()(k) != k)
                dm = true;
            if (k != x.first().direction() && !x.dirOrientation()(k))
                dor = true;
        }
    }
}

// Applies one LCG-chosen reparameterisation to a tensor B-spline patch: one
// of the D reverse(k) and D(D-1)/2 swapDirections(i,j) ops, or nothing. The
// draw is uniform on [0,of), so the patch is left alone with probability
// 1 - nOps/of, where nOps = D + D(D-1)/2.
template<short_t D>
void randomReparam(gsTensorBSpline<D, real_t> & g, Lcg & rng, uint64_t of)
{
    const uint64_t c = rng.next() % of;
    if (c < D)
    {
        g.reverse(unsigned(c));
        return;
    }
    uint64_t q = D;
    for (unsigned i = 0; i != D; ++i)
        for (unsigned j = i + 1; j != D; ++j, ++q)
            if (q == c)
            {
                g.swapDirections(i, j);
                return;
            }
}

// n x m grid of unit squares: optionally inserted in shuffled order, a third
// of the patches reparameterised, all patches rotated by 'angle' about the origin.
gsMultiPatch<real_t> squareGrid(int n, int m, bool shuffle, bool reparam,
                                real_t angle, uint64_t seed = 1)
{
    Lcg rng(seed);
    std::vector<int> order(size_t(n) * m);
    for (size_t k = 0; k != order.size(); ++k)
        order[k] = int(k);
    if (shuffle)
        for (size_t k = order.size() - 1; k > 0; --k)
            std::swap(order[k], order[rng.next() % (k + 1)]);
    gsMultiPatch<real_t> mp;
    for (size_t k = 0; k != order.size(); ++k)
    {
        gsNurbsCreator<>::TensorBSpline2Ptr p =
            gsNurbsCreator<>::BSplineSquare(real_t(1), real_t(order[k] % n), real_t(order[k] / n));
        if (reparam)
            randomReparam<2>(*p, rng, 9);
        mp.addPatch(*p);
    }
    if (angle != real_t(0))
        for (size_t k = 0; k != mp.nPatches(); ++k)
            mp.patch(k).rotate(angle);
    return mp;
}

// n x n x n grid of unit cubes, every patch randomly reparameterised, all
// patches rotated about an oblique axis.
gsMultiPatch<real_t> cubeGridReparam(int n, uint64_t seed)
{
    Lcg rng(seed);
    gsMultiPatch<real_t> mp;
    for (int k = 0; k != n; ++k)
        for (int j = 0; j != n; ++j)
            for (int i = 0; i != n; ++i)
            {
                gsNurbsCreator<>::TensorBSpline3Ptr p =
                    gsNurbsCreator<>::BSplineCube(real_t(1), real_t(i), real_t(j), real_t(k));
                randomReparam<3>(*p, rng, 7);
                mp.addPatch(*p);
            }
    gsVector<real_t, 3> axis;
    axis << 1, 2, 3;
    for (size_t k = 0; k != mp.nPatches(); ++k)
        mp.patch(k).rotate(real_t(0.7), axis);
    return mp;
}

// K separate pairs of unit squares (D=2) or cubes (D=3). Pair j: patch A at
// (3j, 0, ..) shifted by the phase j*2*tol/K on every axis, so the pairs
// sample all positions relative to any grid of cell width <= 2*tol; patch B is
// A shifted by (1+extraX) along x and, if 'offsets', by d_j. The offsets d_j
// cycle through the diagonals (+-0.4 tol per axis, Euclidean norm 0.57 tol in
// 2D, 0.69 tol in 3D) and the single-axis shifts +-0.9 tol; all are below tol.
gsMultiPatch<real_t> straddlePairs(int D, real_t tol, int K, real_t extraX, bool offsets)
{
    std::vector<gsVector<real_t> > off;
    for (int mask = 0; mask != (1 << D); ++mask)
    {
        gsVector<real_t> v(D);
        for (int a = 0; a != D; ++a)
            v(a) = ((mask >> a) & 1 ? real_t(0.4) : real_t(-0.4)) * tol;
        off.push_back(v);
    }
    for (int a = 0; a != D; ++a)
        for (int sgn = -1; sgn <= 1; sgn += 2)
        {
            gsVector<real_t> v = gsVector<real_t>::Zero(D);
            v(a) = real_t(sgn) * real_t(0.9) * tol;
            off.push_back(v);
        }

    gsMultiPatch<real_t> mp;
    for (int j = 0; j != K; ++j)
    {
        const real_t phase = real_t(j) * real_t(2) * tol / real_t(K);
        gsVector<real_t> a = gsVector<real_t>::Constant(D, phase);
        a(0) += real_t(3 * j);
        gsVector<real_t> b = a;
        b(0) += real_t(1) + extraX;
        if (offsets)
            b += off[size_t(j) % off.size()];
        for (int side = 0; side != 2; ++side)
        {
            const gsVector<real_t> & q = side == 0 ? a : b;
            if (D == 2)
                mp.addPatch(*gsNurbsCreator<>::BSplineSquare(real_t(1), q(0), q(1)));
            else
                mp.addPatch(*gsNurbsCreator<>::BSplineCube(real_t(1), q(0), q(1), q(2)));
        }
    }
    return mp;
}

// k exactly coincident unit squares.
gsMultiPatch<real_t> coincidentSquares(int k)
{
    gsMultiPatch<real_t> mp;
    for (int i = 0; i != k; ++i)
        mp.addPatch(*gsNurbsCreator<>::BSplineSquare());
    return mp;
}

} // namespace

SUITE(gsComputeTopology_test)
{
    // Each perturbation of a pairwise snapshot must be rejected by the
    // comparator used in all other tests.
    TEST(ComparatorCanFail)
    {
        gsMultiPatch<real_t> mp = squareGrid(3, 3, false, false, 0);
        Snapshot ref;
        CHECK(agreeWithCounts(mp, real_t(1e-4), false, 12, 12));
        CHECK(agree(mp, real_t(1e-4), false, ref));

        Snapshot same = ref;
        CHECK(sameTopology(ref, same, false));

        Snapshot swapped = ref;
        std::swap(swapped.ifc[0], swapped.ifc[1]);
        CHECK(!sameTopology(ref, swapped, false));

        Snapshot flipped = ref;
        {
            const boundaryInterface & b = flipped.ifc[0];
            gsVector<bool> o = b.dirOrientation();
            const index_t t = 1 - b.first().direction();
            o(t) = !o(t);
            flipped.ifc[0] = boundaryInterface(b.first(), b.second(), b.dirMap(), o);
        }
        CHECK(!sameTopology(ref, flipped, false));

        Snapshot dropped = ref;
        dropped.bdr.pop_back();
        CHECK(!sameTopology(ref, dropped, false));

        Snapshot fewer = ref;
        fewer.ifc.pop_back();
        CHECK(!sameTopology(ref, fewer, false));
    }

    // Counts (n-1)m + n(m-1) / 2(n+m) for squares and 3n^2(n-1) / 6n^2 for
    // cubes, read from the topology the creators compute themselves.
    TEST(KnownGridCounts)
    {
        gsMultiPatch<real_t> sq = gsNurbsCreator<>::BSplineSquareGrid(3, 4);
        CHECK_EQUAL(17u, sq.nInterfaces());
        CHECK_EQUAL(14u, sq.nBoundary());
        CHECK(agreeBothModes(sq, real_t(1e-4), 17, 14));

        gsMultiPatch<real_t> cu = gsNurbsCreator<>::BSplineCubeGrid(4, 4, 4);
        CHECK_EQUAL(144u, cu.nInterfaces());
        CHECK_EQUAL(96u, cu.nBoundary());
        CHECK(agreeBothModes(cu, real_t(1e-4), 144, 96));
    }

    // Shuffled insertion order and mixed orientations: a wrong dirMap or
    // dirOrientation, or a missed match anywhere in the 256 patches, differs
    // from the exhaustive result.
    TEST(LargeSquareGridShuffledRotated)
    {
        gsMultiPatch<real_t> mp = squareGrid(16, 16, true, true, real_t(0.3), 12345);
        CHECK_EQUAL(256u, mp.nPatches());
        for (int co = 0; co != 2; ++co)
        {
            Snapshot ref;
            CHECK(agree(mp, real_t(1e-4), co == 1, ref));
            CHECK(counts(ref, 480, 64));
            bool dm, dor;
            orientationFlags(ref, dm, dor);
            CHECK(dm);
            CHECK(dor);
        }
    }

    TEST(CubeGridReparameterised)
    {
        gsMultiPatch<real_t> mp = cubeGridReparam(4, 777);
        CHECK_EQUAL(64u, mp.nPatches());
        for (int co = 0; co != 2; ++co)
        {
            Snapshot ref;
            CHECK(agree(mp, real_t(1e-4), co == 1, ref));
            CHECK(counts(ref, 144, 96));
            bool dm, dor;
            orientationFlags(ref, dm, dor);
            CHECK(dm);
            CHECK(dor);
        }
    }

    // Matching pairs whose key points lie in different but neighbouring cells
    // (across a cell face, edge or corner): a same-cell-only search, a
    // one-sided {0,+1}^d neighbourhood, a face-only neighbourhood, or cells
    // much narrower than tol miss some pair.
    TEST(CellBoundaryStraddle2D)
    {
        const real_t tol = 1e-3;
        gsMultiPatch<real_t> mp = straddlePairs(2, tol, 16, 0, true);
        CHECK(agreeBothModes(mp, tol, 16, 96));
    }

    TEST(CellBoundaryStraddle3D)
    {
        const real_t tol = 1e-3;
        gsMultiPatch<real_t> mp = straddlePairs(3, tol, 16, 0, true);
        CHECK(agreeBothModes(mp, tol, 16, 160));
    }

    // The key points of a pair 0.99 tol apart are close, and of a pair
    // 1.01 tol apart are not; the candidate sets are filtered by the exact
    // test either way.
    TEST(NearMissAtTol)
    {
        const real_t tol = 1e-3;
        gsMultiPatch<real_t> in = straddlePairs(2, tol, 16, real_t(0.99) * tol, false);
        CHECK(agreeBothModes(in, tol, 16, 96));
        gsMultiPatch<real_t> out = straddlePairs(2, tol, 16, real_t(1.01) * tol, false);
        CHECK(agreeBothModes(out, tol, 0, 128));
        gsMultiPatch<real_t> in3 = straddlePairs(3, tol, 16, real_t(0.99) * tol, false);
        CHECK(agreeBothModes(in3, tol, 16, 160));
        gsMultiPatch<real_t> out3 = straddlePairs(3, tol, 16, real_t(1.01) * tol, false);
        CHECK(agreeBothModes(out3, tol, 0, 192));
    }

    TEST(CornersOnlyKeyPoint)
    {
        const real_t tol = 1e-4;
        {
            // Degree-2 squares A=[0,1]^2 and B=[1,2]x[0,1]. B's west middle
            // control point (u-fastest coefficient row 3) moves by 0.2, which
            // moves the west side centre by 0.2/2 = 0.1 (Bezier weight 1/2) and
            // leaves the corners in place.
            gsNurbsCreator<>::TensorBSpline2Ptr a = gsNurbsCreator<>::BSplineSquareDeg(2);
            gsNurbsCreator<>::TensorBSpline2Ptr b = gsNurbsCreator<>::BSplineSquareDeg(2);
            gsVector<real_t> sh(2);
            sh << 1, 0;
            b->translate(sh);
            CHECK_EQUAL(9, b->coefs().rows());
            CHECK_CLOSE(1.0, b->coef(3, 0), 1e-14);
            CHECK_CLOSE(0.5, b->coef(3, 1), 1e-14);
            b->coef(3, 0) += real_t(0.2);
            gsMultiPatch<real_t> mp;
            mp.addPatch(*a);
            mp.addPatch(*b);
            CHECK(agreeWithCounts(mp, tol, false, 0, 8));
            CHECK(agreeWithCounts(mp, tol, true, 1, 6));
        }
        {
            // Shared line x=1: the corner minima coincide at (1,0), the
            // corners do not, so the cornersOnly candidate must be rejected.
            gsMultiPatch<real_t> mp;
            mp.addPatch(*gsNurbsCreator<>::BSplineSquare());
            gsMatrix<real_t> box(2, 2);
            box << 1, 2, 0, 1.5;
            mp.addPatch(*gsNurbsCreator<>::BSplineSquare(box));
            CHECK(agreeBothModes(mp, tol, 0, 8));
        }
    }

    // Two trilinear patches share a face with bitwise-identical corners, but
    // the second patch has u and v swapped, so the two patches list the face
    // corners in a different order. The corners are ~1e6 with tol = 1e-11
    // (below their resolution) and the non-shared faces are mirrored, so a
    // key formed from the corner sum cancels to O(0.1), passes the overflow
    // guard, and rounds differently per patch: it lands two cells away from
    // its partner and the interface is missed. The cornersOnly key must be
    // independent of the corner order.
    TEST(CornersOnlyKeyCornerOrder3D)
    {
        const real_t X[4][3] = { { real_t(-1000000.5209384176), 0, 0 },
                                 { real_t( 1000000.393255095 ), 0, 0 },
                                 { real_t(-300000.4896935205 ), 1, 0 },
                                 { real_t( 300000.02957496396), 1, 0 } };
        gsKnotVector<real_t> kv(0, 1, 0, 2);
        gsTensorBSplineBasis<3, real_t> basis(kv, kv, kv);
        gsMultiPatch<real_t> mp;
        for (int swapUV = 0; swapUV < 2; ++swapUV)
        {
            gsMatrix<real_t> C(8, 3);
            for (int k = 0; k < 2; ++k)
                for (int j = 0; j < 2; ++j)
                    for (int i = 0; i < 2; ++i)
                    {
                        const real_t * p = swapUV ? X[j + 2*i] : X[i + 2*j];
                        const bool shared = swapUV ? (k == 0) : (k == 1);
                        const real_t sg = shared ? real_t(1) : real_t(-1);
                        const int row = i + 2*j + 4*k;
                        C(row, 0) = sg * p[0];
                        C(row, 1) = sg * p[1];
                        C(row, 2) = p[2] + real_t(swapUV ? k : k - 1);
                    }
            mp.addPatch(gsTensorBSpline<3, real_t>(basis, C));
        }
        const real_t tol = real_t(1e-11);
        CHECK(agreeWithCounts(mp, tol, true, 1, 10));
        CHECK(agreeWithCounts(mp, tol, false, 0, 12));
    }

    // Three physical coordinates for two parametric directions; the fold keeps
    // the 3-D shared edges at x=2.
    TEST(EmbeddedSurface)
    {
        gsMultiPatch<real_t> mp = squareGrid(4, 2, false, false, 0);
        mp.embed(3);
        CHECK_EQUAL(3, mp.geoDim());
        CHECK_EQUAL(2, mp.parDim());
        gsVector<real_t> back(3), fwd(3);
        back << -2, 0, 0;
        fwd << 2, 0, 0;
        gsVector<real_t, 3> yAxis;
        yAxis << 0, 1, 0;
        for (size_t k = 0; k != mp.nPatches(); ++k)
        {
            if (k % 4 < 2) // patches with x < 2
                continue;
            mp.patch(k).translate(back);
            mp.patch(k).rotate(real_t(EIGEN_PI / 2), yAxis);
            mp.patch(k).translate(fwd);
        }
        gsVector<real_t, 3> axis;
        axis << 1, -2, 0.5;
        for (size_t k = 0; k != mp.nPatches(); ++k)
            mp.patch(k).rotate(real_t(0.9), axis);
        CHECK(agreeBothModes(mp, real_t(1e-4), 10, 12));
    }

    // Sides with several matches: the interface order follows ascending
    // (side, other) only if every side collects all candidates and sorts them.
    TEST(CoincidentDuplicates)
    {
        {
            gsMultiPatch<real_t> mp = coincidentSquares(2);
            CHECK(agreeBothModes(mp, real_t(1e-4), 4, 0));
        }
        {
            gsMultiPatch<real_t> mp = coincidentSquares(3);
            CHECK(agreeBothModes(mp, real_t(1e-4), 12, 0));
        }
        {
            gsMultiPatch<real_t> mp = squareGrid(2, 2, false, false, 0);
            mp.addPatch(mp.patch(0));
            CHECK(agreeBothModes(mp, real_t(1e-4), 10, 6));
        }
        {
            const real_t tol = 1e-3;
            Lcg rng(99);
            gsMultiPatch<real_t> mp;
            for (int c = 0; c != 3; ++c)
            {
                gsVector<real_t> jitter(2);
                jitter << real_t(0.17) * tol * rng.sym(), real_t(0.17) * tol * rng.sym();
                gsMultiPatch<real_t> g = squareGrid(2, 2, false, false, 0);
                for (size_t k = 0; k != g.nPatches(); ++k)
                {
                    g.patch(k).translate(jitter);
                    mp.addPatch(g.patch(k));
                }
            }
            CHECK(agreeBothModes(mp, tol, 84, 0));
        }
    }

    // tol <= 0 matches nothing; there is no cell size, so only the exhaustive
    // path can be taken. Only crashes, asserts and a changed result are caught.
    TEST(NonPositiveTolFallback)
    {
        gsMultiPatch<real_t> mp = squareGrid(4, 4, false, false, 0);
        CHECK(agreeBothModes(mp, real_t(0), 0, 64));
        CHECK(agreeBothModes(mp, real_t(-1), 0, 64));
    }

    // Coordinates of 1e300: the cell indices are not representable and the
    // squared distances overflow. A missing guard usually only over-collects
    // candidates, so this catches crashes, asserts, traps and a changed
    // result, not a silently absent guard.
    TEST(OverflowFallback)
    {
        gsMultiPatch<real_t> mp = squareGrid(3, 3, false, false, 0);
        const real_t s = real_t(std::ldexp(1.0, 996));
        for (size_t k = 0; k != mp.nPatches(); ++k)
            mp.patch(k).scale(s);
        CHECK(agreeBothModes(mp, real_t(1e-4), 12, 12));
    }

    TEST(EmptyMultiPatch)
    {
        gsMultiPatch<real_t> mp;
        CHECK(mp.computeTopology());
        CHECK_EQUAL(0u, mp.nInterfaces());
        CHECK_EQUAL(0u, mp.nBoundary());
        CHECK((mp.*pairwise)(real_t(1e-4), false));
        CHECK_EQUAL(0u, mp.nInterfaces());
        CHECK_EQUAL(0u, mp.nBoundary());
        CHECK(mp.computeTopology(real_t(1e-4), true));
        CHECK((mp.*pairwise)(real_t(1e-4), true));
        CHECK_EQUAL(0u, mp.nInterfaces() + mp.nBoundary());
    }
}
