/** @file gsSpaceFillingCurve_test.cpp

    @brief Tests the Morton/Hilbert space-filling curve encoders in gsSpaceFillingCurve.h

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst
**/

#include "gismo_unittest.h"

// ---------------------------------------------------------------------------
// File-local helpers (plain free functions, not fixtures)
// ---------------------------------------------------------------------------

/// L1 distance of two lattice points, computed without signed casts.
static uint64_t l1dist(const uint64_t * a, const uint64_t * b, short_t dim)
{
    uint64_t d = 0;
    for (short_t i = 0; i != dim; ++i)
        d += (a[i] > b[i] ? a[i] - b[i] : b[i] - a[i]);
    return d;
}

/// Enumerate the full lattice: writes the \a idx-th of the 2^(dim*bits)
/// lattice points into \a c (axis 0 varies slowest).
static void latticePoint(uint64_t idx, short_t dim, unsigned bits, uint64_t * c)
{
    const uint64_t L = 1ULL << bits;
    for (short_t i = dim; i-- > 0; ) { c[i] = idx % L; idx /= L; }
}

static bool sameCoords(const uint64_t * a, const uint64_t * b, short_t dim)
{
    for (short_t i = 0; i != dim; ++i) if (a[i] != b[i]) return false;
    return true;
}

// Small box builders (column 0 = lower corner, column 1 = upper corner).
static gsMatrix<real_t> box1(real_t x0, real_t x1)
{
    gsMatrix<real_t> b(1,2);
    b(0,0) = x0; b(0,1) = x1;
    return b;
}
static gsMatrix<real_t> box2(real_t x0, real_t y0, real_t x1, real_t y1)
{
    gsMatrix<real_t> b(2,2);
    b(0,0) = x0; b(1,0) = y0; b(0,1) = x1; b(1,1) = y1;
    return b;
}
static gsMatrix<real_t> box3(real_t x0, real_t y0, real_t z0,
                             real_t x1, real_t y1, real_t z1)
{
    gsMatrix<real_t> b(3,2);
    b(0,0) = x0; b(1,0) = y0; b(2,0) = z0;
    b(0,1) = x1; b(1,1) = y1; b(2,1) = z1;
    return b;
}

// --- helpers used inside CHECK_THROW ---------------------------------------
// Never CHECK_THROW on a bare constructor expression: `gsSpaceFillingCurve(box);`
// is the most-vexing-parse (it declares a variable `box` and calls the default
// ctor, which never throws), and a two-argument ctor call would additionally be
// split by the macro on its comma.  Each helper returns a value so that the
// construction cannot be elided.

static unsigned bitsPerAxisOf(short_t dim)
{ return gsSpaceFillingCurve::bitsPerAxis(dim); }

static short_t buildAndGetGeoDim(const gsMatrix<real_t> & box)
{ gsSpaceFillingCurve sfc(box); return sfc.geoDim(); }

static unsigned buildAndGetBits(const gsMatrix<real_t> & box, real_t relTol)
{ gsSpaceFillingCurve sfc(box, gsSpaceFillingCurve::Hilbert, relTol); return sfc.bits(); }

/// A geometric point as a gsVector
static gsVector<real_t> pt3(real_t x, real_t y, real_t z)
{
    gsVector<real_t> p(3);
    p[0] = x; p[1] = y; p[2] = z;
    return p;
}
static gsVector<real_t> pt2(real_t x, real_t y)
{
    gsVector<real_t> p(2);
    p[0] = x; p[1] = y;
    return p;
}
static gsVector<real_t> pt1(real_t x)
{
    gsVector<real_t> p(1);
    p[0] = x;
    return p;
}

// The four (dim,bits) lattices exercised by the looped tests
static const short_t  s_dims [4] = { 2, 3, 2, 3 };
static const unsigned s_bits [4] = { 3, 3, 8, 5 };

SUITE(gsSpaceFillingCurve_test)
{

// ---------------------------------------------------------------------------
TEST(BitsPerAxis)
{
    // Catches: a wrong per-axis bit budget, i.e. a packed key that no longer
    // fits into a single uint64_t.
    CHECK_EQUAL(0u,  gsSpaceFillingCurve::bitsPerAxis(0));
    CHECK_EQUAL(63u, gsSpaceFillingCurve::bitsPerAxis(1));
    CHECK_EQUAL(31u, gsSpaceFillingCurve::bitsPerAxis(2));
    CHECK_EQUAL(21u, gsSpaceFillingCurve::bitsPerAxis(3));

    // The packed key must fit in 64 bits: 0, 63, 62, 63
    for (short_t d = 0; d <= 3; ++d)
        CHECK( (unsigned)d * gsSpaceFillingCurve::bitsPerAxis(d) <= 64u );

    // Out-of-range dimensions are rejected by GISMO_ENSURE (live in every build)
    CHECK_THROW(bitsPerAxisOf(4),  std::runtime_error);
    CHECK_THROW(bitsPerAxisOf(-1), std::runtime_error);
}

// ---------------------------------------------------------------------------
TEST(MortonRoundTripSmall)
{
    // Catches: a broken interleave/deinterleave inverse pair, and a Morton
    // encoder that is not exactly the bit interleaving.
    uint64_t c[3], r[3];

    // --- fixed vectors pinning the bit convention -------------------------
    // MSB of X[0] is the MSB of the key; bit j of X[i] sits at key bit dim*j+(dim-1-i)
    c[0] = 2; c[1] = 1;                                   // dim=2, bits=2
    CHECK_EQUAL((uint64_t)9, gsSpaceFillingCurve::mortonEncode(c, (short_t)2, 2u));
    gsSpaceFillingCurve::mortonDecode((uint64_t)9, (short_t)2, 2u, r);
    CHECK_EQUAL((uint64_t)2, r[0]);
    CHECK_EQUAL((uint64_t)1, r[1]);

    c[0] = 1; c[1] = 0; c[2] = 0;                         // dim=3, bits=21
    CHECK_EQUAL((uint64_t)4, gsSpaceFillingCurve::mortonEncode(c, (short_t)3, 21u));
    c[0] = 0; c[1] = 1; c[2] = 0;
    CHECK_EQUAL((uint64_t)2, gsSpaceFillingCurve::mortonEncode(c, (short_t)3, 21u));
    c[0] = 0; c[1] = 0; c[2] = 1;
    CHECK_EQUAL((uint64_t)1, gsSpaceFillingCurve::mortonEncode(c, (short_t)3, 21u));

    c[0] = 1; c[1] = 0;                                   // dim=2, bits=31
    CHECK_EQUAL((uint64_t)2, gsSpaceFillingCurve::mortonEncode(c, (short_t)2, 31u));
    c[0] = 0; c[1] = 1;
    CHECK_EQUAL((uint64_t)1, gsSpaceFillingCurve::mortonEncode(c, (short_t)2, 31u));
    c[0] = 2; c[1] = 0;
    CHECK_EQUAL((uint64_t)8, gsSpaceFillingCurve::mortonEncode(c, (short_t)2, 31u));

    // dim==1 interleaving is the identity
    const uint64_t ones[5] = { 0, 1, 7, 123456789, (1ULL<<62) };
    for (int i = 0; i != 5; ++i)
    {
        c[0] = ones[i];
        CHECK_EQUAL(ones[i], gsSpaceFillingCurve::mortonEncode(c, (short_t)1, 63u));
    }

    // --- full-lattice round trip on all four lattices ---------------------
    for (int t = 0; t != 4; ++t)
    {
        const short_t  dim  = s_dims[t];
        const unsigned bits = s_bits[t];
        const uint64_t N    = 1ULL << ((unsigned)dim * bits);
        uint64_t bad = 0, firstBad = 0, visited = 0;
        for (uint64_t idx = 0; idx != N; ++idx)
        {
            latticePoint(idx, dim, bits, c);
            ++visited;
            const uint64_t key = gsSpaceFillingCurve::mortonEncode(c, dim, bits);
            gsSpaceFillingCurve::mortonDecode(key, dim, bits, r);
            if (!sameCoords(c, r, dim)) { if (0 == bad) firstBad = idx; ++bad; }
        }
        CHECK_EQUAL((uint64_t)0, bad);
        CHECK_EQUAL(N, visited);
        if (0 != bad)
            gsInfo << "Morton round trip: first violation at lattice index "
                   << firstBad << " (dim=" << dim << ", bits=" << bits << ")\n";
    }

    // --- Morton must be *exactly* interleave ------------------------------
    {
        const short_t  dim  = 2;
        const unsigned bits = 3;
        const uint64_t N    = 1ULL << ((unsigned)dim * bits);
        uint64_t bad = 0;
        for (uint64_t idx = 0; idx != N; ++idx)
        {
            latticePoint(idx, dim, bits, c);
            if ( gsSpaceFillingCurve::mortonEncode(c, dim, bits) !=
                 gsSpaceFillingCurve::interleave  (c, dim, bits) ) ++bad;
        }
        CHECK_EQUAL((uint64_t)0, bad);
    }
}

// ---------------------------------------------------------------------------
TEST(HilbertRoundTripSmall)
{
    // Catches: a mis-transcribed Skilling AxestoTranspose/TransposetoAxes pair.
    uint64_t c[3], r[3];

    // --- orientation-pinning values at dim=2, bits=1 ----------------------
    // These pin the curve itself, not merely its bijectivity.
    const uint64_t cx[4] = { 0, 0, 1, 1 };
    const uint64_t cy[4] = { 0, 1, 1, 0 };
    const uint64_t ck[4] = { 0, 1, 2, 3 };
    for (int i = 0; i != 4; ++i)
    {
        c[0] = cx[i]; c[1] = cy[i];
        CHECK_EQUAL(ck[i], gsSpaceFillingCurve::hilbertEncode(c, (short_t)2, 1u));
        gsSpaceFillingCurve::hilbertDecode(ck[i], (short_t)2, 1u, r);
        CHECK_EQUAL(cx[i], r[0]);
        CHECK_EQUAL(cy[i], r[1]);
    }

    // --- the origin maps to key 0 at production widths --------------------
    const short_t  zdim [3] = { 1, 2, 3 };
    const unsigned zbits[3] = { 63u, 31u, 21u };
    for (int i = 0; i != 3; ++i)
    {
        c[0] = c[1] = c[2] = 0;
        CHECK_EQUAL((uint64_t)0, gsSpaceFillingCurve::hilbertEncode(c, zdim[i], zbits[i]));
        r[0] = r[1] = r[2] = 12345;
        gsSpaceFillingCurve::hilbertDecode((uint64_t)0, zdim[i], zbits[i], r);
        for (short_t k = 0; k != zdim[i]; ++k) CHECK_EQUAL((uint64_t)0, r[k]);
    }

    // --- full-lattice round trip on all four lattices ---------------------
    for (int t = 0; t != 4; ++t)
    {
        const short_t  dim  = s_dims[t];
        const unsigned bits = s_bits[t];
        const uint64_t N    = 1ULL << ((unsigned)dim * bits);
        uint64_t bad = 0, firstBad = 0, visited = 0;
        for (uint64_t idx = 0; idx != N; ++idx)
        {
            latticePoint(idx, dim, bits, c);
            ++visited;
            const uint64_t key = gsSpaceFillingCurve::hilbertEncode(c, dim, bits);
            gsSpaceFillingCurve::hilbertDecode(key, dim, bits, r);
            if (!sameCoords(c, r, dim)) { if (0 == bad) firstBad = idx; ++bad; }
        }
        CHECK_EQUAL((uint64_t)0, bad);
        CHECK_EQUAL(N, visited);
        if (0 != bad)
            gsInfo << "Hilbert round trip: first violation at lattice index "
                   << firstBad << " (dim=" << dim << ", bits=" << bits << ")\n";
    }
}

// ---------------------------------------------------------------------------
TEST(Bijectivity)
{
    // Catches: collisions or gaps in the key map, i.e. an encoder that loses
    // information.  The key set must be exactly [0, 2^(dim*bits)).
    const gsSpaceFillingCurve::Curve curves[2] =
        { gsSpaceFillingCurve::Morton, gsSpaceFillingCurve::Hilbert };

    uint64_t c[3];
    for (int ci = 0; ci != 2; ++ci)
        for (int t = 0; t != 4; ++t)
        {
            const short_t  dim  = s_dims[t];
            const unsigned bits = s_bits[t];
            const uint64_t N    = 1ULL << ((unsigned)dim * bits);

            std::vector<char> seen((std::size_t)N, 0);
            uint64_t outOfRange = 0, collisions = 0, firstCollision = 0;
            uint64_t pointsVisited = 0;

            for (uint64_t idx = 0; idx != N; ++idx)
            {
                latticePoint(idx, dim, bits, c);
                ++pointsVisited;
                const uint64_t key = gsSpaceFillingCurve::encode(curves[ci], c, dim, bits);
                if (key >= N) { ++outOfRange; continue; }
                if (0 != seen[(std::size_t)key])
                { if (0 == collisions) firstCollision = key; ++collisions; }
                seen[(std::size_t)key] = 1;
            }
            uint64_t gaps = 0;
            for (uint64_t k = 0; k != N; ++k) if (0 == seen[(std::size_t)k]) ++gaps;

            CHECK_EQUAL((uint64_t)0, outOfRange);
            CHECK_EQUAL((uint64_t)0, collisions);
            CHECK_EQUAL((uint64_t)0, gaps);
            // non-vacuity: a mis-written loop bound must not pass by visiting nothing
            CHECK_EQUAL(N, pointsVisited);
            if (0 != collisions)
                gsInfo << "Bijectivity: first collision at key " << firstCollision
                       << " (curve=" << (int)curves[ci] << ", dim=" << dim
                       << ", bits=" << bits << ")\n";
        }
}

// ---------------------------------------------------------------------------
TEST(HilbertAdjacencyMortonNot)
{
    // Catches: the two encoders wired to the same implementation, and any
    // Hilbert transcription that is a bijection but not a Hilbert curve.
    uint64_t a[3], b[3];

    for (int t = 0; t != 4; ++t)
    {
        const short_t  dim  = s_dims[t];
        const unsigned bits = s_bits[t];
        const uint64_t N    = 1ULL << ((unsigned)dim * bits);

        // Hilbert: consecutive keys are L1-neighbours, everywhere.
        uint64_t hilbertViolations = 0, firstBad = 0;
        for (uint64_t k = 0; k + 1 < N; ++k)
        {
            gsSpaceFillingCurve::hilbertDecode(k,   dim, bits, a);
            gsSpaceFillingCurve::hilbertDecode(k+1, dim, bits, b);
            if (1 != l1dist(a, b, dim))
            { if (0 == hilbertViolations) firstBad = k; ++hilbertViolations; }
        }
        CHECK_EQUAL((uint64_t)0, hilbertViolations);
        if (0 != hilbertViolations)
            gsInfo << "Hilbert adjacency: first violation at key " << firstBad
                   << " (dim=" << dim << ", bits=" << bits << ")\n";

        // Morton: must NOT be adjacent everywhere (only meaningful for dim>=2;
        // at dim==1 both curves are the identity, so Morton *is* adjacent).
        if (dim >= 2)
        {
            uint64_t mortonViolations = 0;
            for (uint64_t k = 0; k + 1 < N; ++k)
            {
                gsSpaceFillingCurve::mortonDecode(k,   dim, bits, a);
                gsSpaceFillingCurve::mortonDecode(k+1, dim, bits, b);
                if (1 != l1dist(a, b, dim)) ++mortonViolations;
            }
            CHECK(mortonViolations > 0);
            gsInfo << "Morton adjacency violations at dim=" << dim << ", bits="
                   << bits << ": " << mortonViolations << " of " << N-1 << " pairs\n";
        }
    }
}

// ---------------------------------------------------------------------------
TEST(EncodersAreDistinct)
{
    // Catches: both curves wired to a single implementation, and a swapped
    // Curve dispatch -- as a direct statement, independent of the adjacency
    // argument of the previous test.
    uint64_t c[3];

    c[0] = 1; c[1] = 1;
    CHECK_EQUAL((uint64_t)3, gsSpaceFillingCurve::mortonEncode (c, (short_t)2, 1u));
    CHECK_EQUAL((uint64_t)2, gsSpaceFillingCurve::hilbertEncode(c, (short_t)2, 1u));
    CHECK( gsSpaceFillingCurve::mortonEncode (c, (short_t)2, 1u) !=
           gsSpaceFillingCurve::hilbertEncode(c, (short_t)2, 1u) );

    c[0] = 1; c[1] = 0;
    CHECK_EQUAL((uint64_t)2, gsSpaceFillingCurve::mortonEncode (c, (short_t)2, 1u));
    CHECK_EQUAL((uint64_t)3, gsSpaceFillingCurve::hilbertEncode(c, (short_t)2, 1u));
    CHECK( gsSpaceFillingCurve::mortonEncode (c, (short_t)2, 1u) !=
           gsSpaceFillingCurve::hilbertEncode(c, (short_t)2, 1u) );

    // The Curve dispatch must not be swapped
    const uint64_t px[3] = { 0, 3, 5 };
    const uint64_t py[3] = { 1, 2, 7 };
    for (int i = 0; i != 3; ++i)
    {
        c[0] = px[i]; c[1] = py[i];
        CHECK_EQUAL( gsSpaceFillingCurve::mortonEncode(c, (short_t)2, 3u),
                     gsSpaceFillingCurve::encode(gsSpaceFillingCurve::Morton,  c, (short_t)2, 3u) );
        CHECK_EQUAL( gsSpaceFillingCurve::hilbertEncode(c, (short_t)2, 3u),
                     gsSpaceFillingCurve::encode(gsSpaceFillingCurve::Hilbert, c, (short_t)2, 3u) );
    }
}

// ---------------------------------------------------------------------------
TEST(MaxBitsSpotValues)
{
    // Catches: a `1 << b` written on int/unsigned instead of 1ULL -- an
    // overflow the small-`bits` loops above cannot see, while the production
    // code paths use exactly bits = 63 / 31 / 21.
    uint64_t c[3], r[3];

    // --- dim = 1, bits = 63 ------------------------------------------------
    {
        const uint64_t vals[4] = { 0, 1, (1ULL<<62), (1ULL<<63)-1 };
        for (int i = 0; i != 4; ++i)
        {
            c[0] = vals[i];
            CHECK_EQUAL(vals[i], gsSpaceFillingCurve::mortonEncode (c, (short_t)1, 63u));
            CHECK_EQUAL(vals[i], gsSpaceFillingCurve::hilbertEncode(c, (short_t)1, 63u));
            r[0] = 7;
            gsSpaceFillingCurve::mortonDecode(vals[i], (short_t)1, 63u, r);
            CHECK_EQUAL(vals[i], r[0]);
            r[0] = 7;
            gsSpaceFillingCurve::hilbertDecode(vals[i], (short_t)1, 63u, r);
            CHECK_EQUAL(vals[i], r[0]);
        }
    }

    // --- dim = 2, bits = 31 ------------------------------------------------
    {
        const uint64_t cx[5] = { 0, 1, 0, (1ULL<<31)-1, 12345 };
        const uint64_t cy[5] = { 0, 0, 1, (1ULL<<31)-1, 67890 };
        for (int i = 0; i != 5; ++i)
        {
            c[0] = cx[i]; c[1] = cy[i];
            const uint64_t km = gsSpaceFillingCurve::mortonEncode (c, (short_t)2, 31u);
            const uint64_t kh = gsSpaceFillingCurve::hilbertEncode(c, (short_t)2, 31u);
            CHECK(km < (1ULL<<62));
            CHECK(kh < (1ULL<<62));
            gsSpaceFillingCurve::mortonDecode(km, (short_t)2, 31u, r);
            CHECK(sameCoords(c, r, (short_t)2));
            gsSpaceFillingCurve::hilbertDecode(kh, (short_t)2, 31u, r);
            CHECK(sameCoords(c, r, (short_t)2));
        }
        c[0] = (1ULL<<31)-1; c[1] = (1ULL<<31)-1;
        CHECK_EQUAL((1ULL<<62)-1, gsSpaceFillingCurve::mortonEncode(c, (short_t)2, 31u));
    }

    // --- dim = 3, bits = 21 ------------------------------------------------
    {
        const uint64_t cx[6] = { 0, 1, 0, 0, (1ULL<<21)-1, 1000 };
        const uint64_t cy[6] = { 0, 0, 1, 0, (1ULL<<21)-1, 2000 };
        const uint64_t cz[6] = { 0, 0, 0, 1, (1ULL<<21)-1, 3000 };
        for (int i = 0; i != 6; ++i)
        {
            c[0] = cx[i]; c[1] = cy[i]; c[2] = cz[i];
            const uint64_t km = gsSpaceFillingCurve::mortonEncode (c, (short_t)3, 21u);
            const uint64_t kh = gsSpaceFillingCurve::hilbertEncode(c, (short_t)3, 21u);
            CHECK(km < (1ULL<<63));
            CHECK(kh < (1ULL<<63));
            gsSpaceFillingCurve::mortonDecode(km, (short_t)3, 21u, r);
            CHECK(sameCoords(c, r, (short_t)3));
            gsSpaceFillingCurve::hilbertDecode(kh, (short_t)3, 21u, r);
            CHECK(sameCoords(c, r, (short_t)3));
        }
        c[0] = (1ULL<<21)-1; c[1] = (1ULL<<21)-1; c[2] = (1ULL<<21)-1;
        CHECK_EQUAL((1ULL<<63)-1, gsSpaceFillingCurve::mortonEncode(c, (short_t)3, 21u));
    }
}

// ---------------------------------------------------------------------------
TEST(Dim1IsIdentity)
{
    // Catches: Skilling's loops being applied at n == 1.  This is the ONLY
    // test that can catch it: a Skilling-transformed 1-D map is still a
    // bijection and still round-trips, so MortonRoundTripSmall,
    // HilbertRoundTripSmall and Bijectivity all pass on the broken version.
    uint64_t c[1], r[1];

    // --- static level, bits = 10 over the full range ----------------------
    {
        uint64_t bad = 0, firstBad = 0;
        const uint64_t N = 1ULL << 10;
        for (uint64_t v = 0; v != N; ++v)
        {
            c[0] = v;
            const uint64_t kh = gsSpaceFillingCurve::hilbertEncode(c, (short_t)1, 10u);
            const uint64_t km = gsSpaceFillingCurve::mortonEncode (c, (short_t)1, 10u);
            r[0] = 0;
            gsSpaceFillingCurve::hilbertDecode(kh, (short_t)1, 10u, r);
            if (kh != v || km != kh || r[0] != v)
            { if (0 == bad) firstBad = v; ++bad; }
        }
        CHECK_EQUAL((uint64_t)0, bad);
        if (0 != bad) gsInfo << "Dim1IsIdentity: first violation at c = " << firstBad << "\n";
    }

    // --- static level, bits = 63 on spot values ---------------------------
    {
        const uint64_t vals[5] = { 0, 1, (1ULL<<31), (1ULL<<62), (1ULL<<63)-1 };
        for (int i = 0; i != 5; ++i)
        {
            c[0] = vals[i];
            CHECK_EQUAL(vals[i], gsSpaceFillingCurve::hilbertEncode(c, (short_t)1, 63u));
            CHECK_EQUAL(gsSpaceFillingCurve::mortonEncode(c, (short_t)1, 63u),
                        gsSpaceFillingCurve::hilbertEncode(c, (short_t)1, 63u));
            r[0] = 0;
            gsSpaceFillingCurve::hilbertDecode(vals[i], (short_t)1, 63u, r);
            CHECK_EQUAL(vals[i], r[0]);
        }
    }

    // --- geometric level: 1-D box [-1,3], both curves ---------------------
    {
        const gsMatrix<real_t> b = box1(-1.0, 3.0);
        gsSpaceFillingCurve sfcM(b, gsSpaceFillingCurve::Morton);
        gsSpaceFillingCurve sfcH(b, gsSpaceFillingCurve::Hilbert);

        CHECK_EQUAL((short_t)1, sfcH.geoDim());
        CHECK_EQUAL((short_t)1, sfcH.curveDim());
        CHECK_EQUAL(63u, sfcH.bits());
        CHECK_EQUAL((std::size_t)1, sfcH.activeAxes().size());
        CHECK_EQUAL((short_t)0, sfcH.activeAxes()[0]);
        CHECK_EQUAL((short_t)1, sfcM.curveDim());
        CHECK_EQUAL(63u, sfcM.bits());

        const real_t xs[6] = { -1.0, 0.0, 0.5, 1.0, 2.0, 3.0 };
        uint64_t keys[6];
        for (int i = 0; i != 6; ++i)
        {
            keys[i] = sfcH.encode(pt1(xs[i]));
            // dim==1 is the identity, so both curves must agree exactly
            CHECK_EQUAL(keys[i], sfcM.encode(pt1(xs[i])));
        }
        // the ordering is the identity: strictly increasing in x
        for (int i = 0; i != 5; ++i)
            CHECK(keys[i] < keys[i+1]);

        // exact endpoints and midpoint (all exactly representable)
        CHECK_EQUAL((uint64_t)0,        keys[0]);   // x = -1  -> t = 0
        CHECK_EQUAL((uint64_t)(1ULL<<62), keys[3]); // x =  1  -> t = 0.5
        CHECK_EQUAL((uint64_t)((1ULL<<63)-1), keys[5]); // x = 3 -> t = 1, clamped to L-1

        // out-of-box clamping
        CHECK_EQUAL((uint64_t)0,              sfcH.encode(pt1(-100.0)));
        CHECK_EQUAL((uint64_t)((1ULL<<63)-1), sfcH.encode(pt1( 100.0)));
    }
}

// ---------------------------------------------------------------------------
TEST(BoxDegenerateAxis)
{
    // Catches: division by a zero extent, and -- the case that matters for
    // geometric partitioning -- sizing the bit budget from geoDim() instead
    // of curveDim().
    const gsSpaceFillingCurve::Curve curves[2] =
        { gsSpaceFillingCurve::Morton, gsSpaceFillingCurve::Hilbert };

    for (int ci = 0; ci != 2; ++ci)
    {
        // --- Case A: flat plate in 3D, lo=(0,0,0), hi=(1,1,0) -------------
        {
            gsSpaceFillingCurve sfc(box3(0,0,0, 1,1,0), curves[ci]);
            CHECK_EQUAL((short_t)3, sfc.geoDim());
            CHECK_EQUAL((short_t)2, sfc.curveDim());
            // 31 bits, NOT 21: the budget follows curveDim(), not geoDim()
            CHECK_EQUAL(31u, sfc.bits());
            CHECK_EQUAL(gsSpaceFillingCurve::bitsPerAxis(sfc.curveDim()), sfc.bits());
            CHECK_EQUAL((std::size_t)2, sfc.activeAxes().size());
            CHECK_EQUAL((short_t)0, sfc.activeAxes()[0]);
            CHECK_EQUAL((short_t)1, sfc.activeAxes()[1]);

            // the degenerate axis quantises to 0 for any input, with no
            // exception and no NaN (no division by the zero extent)
            CHECK_EQUAL((uint64_t)0, sfc.quantize( 0.0,   (short_t)2));
            CHECK_EQUAL((uint64_t)0, sfc.quantize( 1e30,  (short_t)2));
            CHECK_EQUAL((uint64_t)0, sfc.quantize(-1e30,  (short_t)2));

            // encode ignores the degenerate coordinate
            gsMatrix<real_t> pts(3,3);
            pts(0,0)=0.3; pts(1,0)=0.7; pts(2,0)= 0.0;
            pts(0,1)=0.3; pts(1,1)=0.7; pts(2,1)= 5.0;
            pts(0,2)=0.3; pts(1,2)=0.7; pts(2,2)=-5.0;
            const uint64_t k0 = sfc.encode(pts, 0);
            CHECK_EQUAL(k0, sfc.encode(pts, 1));
            CHECK_EQUAL(k0, sfc.encode(pts, 2));
            CHECK_EQUAL(k0, sfc.encode(pt3(0.3, 0.7, 0.0)));
            CHECK(k0 < (1ULL<<62));

            gsVector<real_t> p;
            sfc.decode(k0, p);
            CHECK_EQUAL((index_t)3, p.size());
            CHECK_EQUAL(0.0, p[2]);              // inactive axis gets lo[2] verbatim
            CHECK_CLOSE(0.3, p[0], 1e-9);        // half cell = 0.5/2^31 = 2.3e-10
            CHECK_CLOSE(0.7, p[1], 1e-9);
            for (index_t i = 0; i != p.size(); ++i) CHECK(math::isfinite(p[i]));
        }

        // --- Case B: degenerate MIDDLE axis, lo=(0,0,0), hi=(2,0,3) -------
        {
            gsSpaceFillingCurve sfc(box3(0,0,0, 2,0,3), curves[ci]);
            CHECK_EQUAL((short_t)2, sfc.curveDim());
            CHECK_EQUAL(31u, sfc.bits());
            // assert the ENTRIES: this proves activeAxes() is consulted at
            // encode time, rather than the first curveDim() axes being used
            CHECK_EQUAL((std::size_t)2, sfc.activeAxes().size());
            CHECK_EQUAL((short_t)0, sfc.activeAxes()[0]);
            CHECK_EQUAL((short_t)2, sfc.activeAxes()[1]);

            const uint64_t k0 = sfc.encode(pt3(0.5, 0.0, 1.5));
            CHECK_EQUAL(k0, sfc.encode(pt3(0.5,  7.0, 1.5)));
            CHECK_EQUAL(k0, sfc.encode(pt3(0.5, -7.0, 1.5)));
        }

        // --- Case C: relTol is really used, hi=(1,1,1e-20) ----------------
        {
            const gsMatrix<real_t> b = box3(0,0,0, 1,1,1e-20);
            gsSpaceFillingCurve sfcDefault(b, curves[ci]);              // relTol = 1e-12
            CHECK_EQUAL((short_t)2, sfcDefault.curveDim());
            CHECK_EQUAL(31u, sfcDefault.bits());

            gsSpaceFillingCurve sfcZero(b, curves[ci], 0.0);            // relTol = 0
            CHECK_EQUAL((short_t)3, sfcZero.curveDim());
            CHECK_EQUAL(21u, sfcZero.bits());
        }

        // --- Case D: degenerate axis in 2D, lo=(0,0), hi=(1,0) ------------
        {
            gsSpaceFillingCurve sfc(box2(0,0, 1,0), curves[ci]);
            CHECK_EQUAL((short_t)2, sfc.geoDim());
            CHECK_EQUAL((short_t)1, sfc.curveDim());
            CHECK_EQUAL(63u, sfc.bits());
            CHECK_EQUAL((std::size_t)1, sfc.activeAxes().size());
            CHECK_EQUAL((short_t)0, sfc.activeAxes()[0]);
        }
    }
}

// ---------------------------------------------------------------------------
TEST(BoxFullyDegenerate)
{
    // Catches: a division by zero / NaN key when every extent is zero.
    {
        gsSpaceFillingCurve sfc(box3(1.5,-2.0,7.0, 1.5,-2.0,7.0));
        CHECK_EQUAL((short_t)3, sfc.geoDim());
        CHECK_EQUAL((short_t)0, sfc.curveDim());
        CHECK_EQUAL(0u, sfc.bits());
        CHECK(sfc.activeAxes().empty());

        CHECK_EQUAL((uint64_t)0, sfc.encode(pt3(1.5, -2.0, 7.0)));
        CHECK_EQUAL((uint64_t)0, sfc.encode(pt3(1e30, -1e30, 0.0)));
        CHECK_EQUAL((uint64_t)0, sfc.encode(pt3(0.0, 0.0, 0.0)));
        for (short_t i = 0; i != 3; ++i)
        {
            CHECK_EQUAL((uint64_t)0, sfc.quantize( 0.0,  i));
            CHECK_EQUAL((uint64_t)0, sfc.quantize( 1e30, i));
        }

        gsVector<real_t> p;
        sfc.decode(0, p);
        CHECK_EQUAL((index_t)3, p.size());
        CHECK_EQUAL( 1.5, p[0]);
        CHECK_EQUAL(-2.0, p[1]);
        CHECK_EQUAL( 7.0, p[2]);
        gsVector<real_t> q;
        sfc.decode(123456789, q);            // curveDim 0 returns early
        CHECK_EQUAL((index_t)3, q.size());
        CHECK_EQUAL( 1.5, q[0]);
        CHECK_EQUAL(-2.0, q[1]);
        CHECK_EQUAL( 7.0, q[2]);
    }

    // 1-D single-point box
    {
        gsSpaceFillingCurve sfc(box1(4.0, 4.0));
        CHECK_EQUAL((short_t)0, sfc.curveDim());
        CHECK_EQUAL(0u, sfc.bits());
        CHECK_EQUAL((uint64_t)0, sfc.encode(pt1(4.0)));
        CHECK_EQUAL((uint64_t)0, sfc.encode(pt1(1e30)));
    }

    // default-constructed curve
    {
        gsSpaceFillingCurve sfc;
        CHECK_EQUAL((short_t)0, sfc.geoDim());
        CHECK_EQUAL((short_t)0, sfc.curveDim());
        CHECK_EQUAL(0u, sfc.bits());
        CHECK(sfc.activeAxes().empty());
        CHECK_EQUAL((uint64_t)0, sfc.encode(pt3(1.0, 2.0, 3.0)));
        CHECK_EQUAL((uint64_t)0, sfc.encode(pt1(0.0)));
    }
}

// ---------------------------------------------------------------------------
TEST(BoxRoundTrip)
{
    // Catches: a quantise/decode pair that is not each other's inverse to
    // within half a lattice cell.
    const gsSpaceFillingCurve::Curve curves[2] =
        { gsSpaceFillingCurve::Morton, gsSpaceFillingCurve::Hilbert };

    for (int ci = 0; ci != 2; ++ci)
    {
        // --- 2-D box [0,1] x [0,2], bits = 31 -----------------------------
        // half cell = 0.5*extent/2^31 = 2.3e-10 (x) and 4.7e-10 (y): the
        // 1e-9 tolerance below is lattice-derived, not a magic constant.
        {
            gsSpaceFillingCurve sfc(box2(0,0, 1,2), curves[ci]);
            CHECK_EQUAL(31u, sfc.bits());
            const real_t xs[4] = { 0.0, 0.1, 0.5, 0.9 };
            const real_t ys[4] = { 0.0, 0.2, 1.0, 1.9 };
            for (int i = 0; i != 4; ++i)
            {
                gsVector<real_t> p;
                sfc.decode(sfc.encode(pt2(xs[i], ys[i])), p);
                CHECK_EQUAL((index_t)2, p.size());
                CHECK_CLOSE(xs[i], p[0], 1e-9);
                CHECK_CLOSE(ys[i], p[1], 1e-9);
            }
        }

        // --- 3-D box [-1,1]^3, bits = 21 ----------------------------------
        // half cell = 0.5*2/2^21 = 4.8e-7: hence the 1e-6 tolerance.
        {
            gsSpaceFillingCurve sfc(box3(-1,-1,-1, 1,1,1), curves[ci]);
            CHECK_EQUAL(21u, sfc.bits());
            const real_t xs[3] = { 0.0,  0.3, -0.9 };
            const real_t ys[3] = { 0.0, -0.4,  0.9 };
            const real_t zs[3] = { 0.0,  0.5,  0.1 };
            for (int i = 0; i != 3; ++i)
            {
                gsVector<real_t> p;
                sfc.decode(sfc.encode(pt3(xs[i], ys[i], zs[i])), p);
                CHECK_EQUAL((index_t)3, p.size());
                CHECK_CLOSE(xs[i], p[0], 1e-6);
                CHECK_CLOSE(ys[i], p[1], 1e-6);
                CHECK_CLOSE(zs[i], p[2], 1e-6);
            }

            // clamping: a point far outside encodes without throwing and
            // decodes to a point inside the box
            gsVector<real_t> p;
            sfc.decode(sfc.encode(pt3(1e6, -1e6, 42.0)), p);
            for (index_t i = 0; i != p.size(); ++i)
            {
                CHECK(math::isfinite(p[i]));
                CHECK(p[i] >= -1.0);
                CHECK(p[i] <=  1.0);
            }
        }
    }
}

// ---------------------------------------------------------------------------
TEST(BoxMalformedThrows)
{
    // Catches: a malformed box silently producing garbage keys.  The box
    // preconditions are GISMO_ENSURE (not GISMO_ASSERT), so they throw
    // std::runtime_error in EVERY build -- hence plain CHECK_THROW.

    // hi < lo on some axis
    CHECK_THROW(buildAndGetGeoDim(box3(0,0,0, 1,-1,1)), std::runtime_error);
    CHECK_THROW(buildAndGetGeoDim(box1(1.0, 0.0)),      std::runtime_error);

    // box.cols() != 2
    {
        gsMatrix<real_t> b31(3,1); b31.setZero();
        CHECK_THROW(buildAndGetGeoDim(b31), std::runtime_error);
        gsMatrix<real_t> b33(3,3); b33.setZero();
        CHECK_THROW(buildAndGetGeoDim(b33), std::runtime_error);
    }

    // box.rows() outside [1,3]
    {
        gsMatrix<real_t> b42(4,2); b42.setZero();
        CHECK_THROW(buildAndGetGeoDim(b42), std::runtime_error);
        gsMatrix<real_t> b02(0,2);
        CHECK_THROW(buildAndGetGeoDim(b02), std::runtime_error);
    }

    // relTol < 0
    CHECK_THROW(buildAndGetBits(box3(0,0,0, 1,1,1), -1.0), std::runtime_error);

    // sanity: a well-formed box does NOT throw
    CHECK_EQUAL((short_t)3, buildAndGetGeoDim(box3(0,0,0, 1,1,1)));
}

} // SUITE(gsSpaceFillingCurve_test)
