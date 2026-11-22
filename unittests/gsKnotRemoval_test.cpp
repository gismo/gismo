/** @file gsKnotRemoval_test.cpp

    @brief Unit tests for knot removal (Tiller's algorithm).

    Tests verify exact round-trip: insert one or more knots into a B-spline
    (or tensor-product B-spline), then remove the same knots and check that
    the original coefficients and knot vector are recovered.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H. Verhelst
**/

#include "gismo_unittest.h"

SUITE(gsKnotRemoval_test)
{

// ---------------------------------------------------------------------------
// Helper: insert knot `val` r times into a 1D B-spline (avoids overload ambiguity)
// ---------------------------------------------------------------------------
static void insertKnot1D(gsBSpline<real_t>& spl, real_t val, int r)
{
    for (int i = 0; i < r; ++i)
    {
        std::vector<real_t> k(1, val);
        spl.insertKnots(k.begin(), k.end());
    }
}

// ---------------------------------------------------------------------------
// Helper: build a degree-p clamped B-spline on [0,1] with n interior knots
// and random-ish control points.
// ---------------------------------------------------------------------------
static gsBSpline<real_t> makeCubicSpline()
{
    // degree 3, knot vector [0 0 0 0  0.25 0.5 0.75  1 1 1 1]
    gsKnotVector<real_t> kv(0, 1, 3, 4); // (start, end, interior_knots, degree+1)
    // 7 control points
    gsMatrix<real_t> coefs(7, 2);
    coefs << 0.0, 0.0,
             0.2, 0.5,
             0.4, 0.8,
             0.6, 0.6,
             0.8, 0.3,
             0.9, 0.1,
             1.0, 0.0;
    return gsBSpline<real_t>(kv, coefs);
}

// ---------------------------------------------------------------------------
TEST(InsertThenRemove_1D_SingleKnot)
{
    gsBSpline<real_t> orig = makeCubicSpline();
    gsBSpline<real_t> spl  = orig;

    const real_t knotVal = 0.5; // must be an existing interior knot for exact round-trip

    // Insert once
    insertKnot1D(spl, knotVal, 1);
    CHECK_EQUAL(orig.numCoefs() + 1, spl.numCoefs());

    // Remove once — must succeed
    const index_t removed = spl.removeKnot(knotVal, 1);
    CHECK_EQUAL(1, removed);

    // Coefficient count must be restored
    CHECK_EQUAL(orig.numCoefs(), spl.numCoefs());

    // Coefficients must match original
    CHECK_ARRAY_CLOSE(orig.coefs().data(), spl.coefs().data(),
                      orig.numCoefs() * orig.geoDim(), 1e-10);

    // Knot vectors must match
    CHECK_EQUAL(orig.knots().size(), spl.knots().size());
    for (index_t i = 0; i < (index_t)orig.knots().size(); ++i)
        CHECK_CLOSE(orig.knots()[i], spl.knots()[i], 1e-14);
}

// ---------------------------------------------------------------------------
TEST(InsertThenRemove_1D_MultipleDistinctKnots)
{
    gsBSpline<real_t> orig = makeCubicSpline();
    gsBSpline<real_t> spl  = orig;

    const real_t k1 = 0.25, k2 = 0.75; // must be existing interior knots

    insertKnot1D(spl, k1, 1);
    insertKnot1D(spl, k2, 1);
    CHECK_EQUAL(orig.numCoefs() + 2, spl.numCoefs());

    CHECK_EQUAL(1, spl.removeKnot(k1, 1));
    CHECK_EQUAL(1, spl.removeKnot(k2, 1));

    CHECK_EQUAL(orig.numCoefs(), spl.numCoefs());
    CHECK_ARRAY_CLOSE(orig.coefs().data(), spl.coefs().data(),
                      orig.numCoefs() * orig.geoDim(), 1e-10);
}

// ---------------------------------------------------------------------------
TEST(InsertThenRemove_1D_Multiplicity2)
{
    gsBSpline<real_t> orig = makeCubicSpline();
    gsBSpline<real_t> spl  = orig;

    const real_t knotVal = 0.5; // already present once; degree 3 allows up to mult 3

    // Insert 2 more times (total multiplicity becomes 3)
    insertKnot1D(spl, knotVal, 2);
    CHECK_EQUAL(orig.numCoefs() + 2, spl.numCoefs());

    // Remove 2 times
    const index_t removed = spl.removeKnot(knotVal, 2);
    CHECK_EQUAL(2, removed);

    CHECK_EQUAL(orig.numCoefs(), spl.numCoefs());
    CHECK_ARRAY_CLOSE(orig.coefs().data(), spl.coefs().data(),
                      orig.numCoefs() * orig.geoDim(), 1e-10);
}

// ---------------------------------------------------------------------------
TEST(RemoveFails_OnFreshSpline)
{
    // A fresh spline has no "artificially inserted" knots: exact removal
    // of an interior knot should fail (return 0) because the geometry
    // cannot be represented with fewer control points exactly.
    gsBSpline<real_t> spl = makeCubicSpline();

    // Insert 0.5 once to get multiplicity 2, then perturb a control point
    // in the affected span so removal becomes geometrically impossible.
    // For knot 0.5 with s=2: first=r-p=6-3=3, last=r-s=6-2=4, sz=4 (even).
    // The feasibility check compares two distinct temp entries — a perturbation
    // of coef 3 or 4 will break it.
    insertKnot1D(spl, 0.5, 1);
    spl.coef(3) += gsVector<real_t,2>::Ones() * 0.123;

    // Try to remove the (perturbed) 0.5 knot — should fail
    const index_t removed = spl.removeKnot(0.5, 1);
    CHECK_EQUAL(0, removed);
}

// ---------------------------------------------------------------------------
// Tensor-product (2D) tests
// ---------------------------------------------------------------------------
static gsTensorBSpline<2, real_t> makeBilinearPatch()
{
    // degree (2,2), knot vectors [0 0 0  0.5  1 1 1]
    gsKnotVector<real_t> kv(0, 1, 1, 3); // 1 interior knot, degree 2 => 4 basis funcs

    gsMatrix<real_t> coefs(16, 2); // 4x4 control points
    coefs << 0.0, 0.0,   0.3, 0.1,   0.7, 0.1,   1.0, 0.0,
             0.1, 0.3,   0.3, 0.4,   0.7, 0.4,   0.9, 0.3,
             0.1, 0.7,   0.3, 0.8,   0.7, 0.8,   0.9, 0.7,
             0.0, 1.0,   0.3, 0.9,   0.7, 0.9,   1.0, 1.0;

    return gsTensorBSpline<2, real_t>(kv, kv, coefs);
}

// ---------------------------------------------------------------------------
TEST(InsertThenRemove_Tensor2D_Dir0)
{
    gsTensorBSpline<2, real_t> orig = makeBilinearPatch();
    gsTensorBSpline<2, real_t> spl  = orig;

    const real_t knotVal = 0.25;

    spl.insertKnot(knotVal, 0, 1);
    CHECK_EQUAL(orig.coefsSize() + 4, spl.coefsSize()); // 4 extra (4 rows in dir 1)

    const int removed = spl.removeKnot(knotVal, 0, 1);
    CHECK_EQUAL(1, removed);

    CHECK_EQUAL(orig.coefsSize(), spl.coefsSize());
    CHECK_ARRAY_CLOSE(orig.coefs().data(), spl.coefs().data(),
                      orig.coefsSize() * orig.geoDim(), 1e-10);
}

// ---------------------------------------------------------------------------
TEST(InsertThenRemove_Tensor2D_Dir1)
{
    gsTensorBSpline<2, real_t> orig = makeBilinearPatch();
    gsTensorBSpline<2, real_t> spl  = orig;

    const real_t knotVal = 0.75;

    spl.insertKnot(knotVal, 1, 1);
    const int removed = spl.removeKnot(knotVal, 1, 1);
    CHECK_EQUAL(1, removed);

    CHECK_EQUAL(orig.coefsSize(), spl.coefsSize());
    CHECK_ARRAY_CLOSE(orig.coefs().data(), spl.coefs().data(),
                      orig.coefsSize() * orig.geoDim(), 1e-10);
}

// ---------------------------------------------------------------------------
// Tolerance behaviour
// ---------------------------------------------------------------------------

// The feasibility threshold of Tiller's algorithm is RELATIVE to the largest
// control-point entry.  A removal that succeeds on an O(1) control net must
// therefore still succeed after the net is scaled up: the round-off residual
// scales with the net, so an absolute threshold would reject it.
TEST(Tolerance_ScaleInvariance_Removable)
{
    const real_t scale = 1e9;

    gsBSpline<real_t> small = makeCubicSpline();
    insertKnot1D(small, 0.5, 1);

    gsBSpline<real_t> big = makeCubicSpline();
    big.coefs() *= scale;
    insertKnot1D(big, 0.5, 1);

    CHECK_EQUAL(1, small.removeKnot(0.5, 1));
    CHECK_EQUAL(1, big.removeKnot(0.5, 1));
}

// Converse of the above: a relative threshold must not become vacuous on a
// large control net.  A perturbation that is large RELATIVE to the net must
// still be rejected, however big the net is.
TEST(Tolerance_ScaleInvariance_NotRemovable)
{
    const real_t scale = 1e9;

    gsBSpline<real_t> big = makeCubicSpline();
    big.coefs() *= scale;
    insertKnot1D(big, 0.5, 1);
    big.coef(3) += gsVector<real_t,2>::Ones() * (scale * 1e-3);

    CHECK_EQUAL(0, big.removeKnot(0.5, 1));
}

// The tolerance is a parameter, not a hardcoded constant: the same perturbed
// spline is irremovable at the default tolerance and removable at a loose one.
TEST(Tolerance_ExplicitArgument)
{
    gsBSpline<real_t> perturbed = makeCubicSpline();
    insertKnot1D(perturbed, 0.5, 1);
    perturbed.coef(3) += gsVector<real_t,2>::Ones() * 1e-6;

    gsBSpline<real_t> tight = perturbed;
    gsBSpline<real_t> loose = perturbed;

    CHECK_EQUAL(0, tight.removeKnot(0.5, 1));                    // default 1e-10
    CHECK_EQUAL(1, loose.removeKnot(0.5, 1, (real_t)1e-3));
}

// ---------------------------------------------------------------------------
// Tensor-product (3D) tests
// ---------------------------------------------------------------------------
static gsTensorBSpline<3, real_t> makeTrivariatePatch()
{
    // degree (2,2,2), knot vectors [0 0 0  0.5  1 1 1] in every direction
    gsKnotVector<real_t> kv(0, 1, 1, 3);
    gsTensorBSplineBasis<3, real_t> basis(kv, kv, kv);

    // A mildly distorted identity map, so that no coefficient fiber is
    // degenerate and the removal has something to reproduce.
    gsMatrix<real_t> coefs = basis.anchors().transpose();
    coefs.col(0) += 0.10 * coefs.col(1);
    coefs.col(1) -= 0.05 * coefs.col(2);
    coefs.col(2) += 0.07 * coefs.col(0);

    return gsTensorBSpline<3, real_t>(basis, give(coefs));
}

TEST(InsertThenRemove_Tensor3D_AllDirections)
{
    const gsTensorBSpline<3, real_t> orig = makeTrivariatePatch();
    const real_t knotVal = 0.25;

    for (short_t dir = 0; dir != 3; ++dir)
    {
        gsTensorBSpline<3, real_t> spl = orig;

        spl.insertKnot(knotVal, dir, 1);
        CHECK_EQUAL(orig.coefsSize() + 16, spl.coefsSize()); // 4x4 slab

        const index_t removed = spl.removeKnot(knotVal, dir, 1);
        CHECK_EQUAL(1, removed);

        CHECK_EQUAL(orig.coefsSize(), spl.coefsSize());
        CHECK_ARRAY_CLOSE(orig.coefs().data(), spl.coefs().data(),
                          orig.coefsSize() * orig.geoDim(), 1e-10);
    }
}

TEST(InsertThenRemove_Tensor3D_Multiplicity2)
{
    const gsTensorBSpline<3, real_t> orig = makeTrivariatePatch();
    const real_t knotVal = 0.5; // present once; degree 2 allows up to mult 2

    gsTensorBSpline<3, real_t> spl = orig;
    spl.insertKnot(knotVal, 1, 1);

    CHECK_EQUAL(1, spl.removeKnot(knotVal, 1, 1));
    CHECK_EQUAL(orig.coefsSize(), spl.coefsSize());
    CHECK_ARRAY_CLOSE(orig.coefs().data(), spl.coefs().data(),
                      orig.coefsSize() * orig.geoDim(), 1e-10);
}

TEST(Tolerance_Tensor3D_ScaleInvariance)
{
    const real_t scale = 1e9;

    gsTensorBSpline<3, real_t> big = makeTrivariatePatch();
    big.coefs() *= scale;
    big.insertKnot(0.25, 2, 1);

    CHECK_EQUAL(1, big.removeKnot(0.25, 2, 1));

    // ... and a relatively large perturbation is still rejected.
    gsTensorBSpline<3, real_t> bad = makeTrivariatePatch();
    bad.coefs() *= scale;
    bad.insertKnot(0.25, 2, 1);
    bad.coefs().row(5) += gsVector<real_t,3>::Ones().transpose() * (scale * 1e-3);

    CHECK_EQUAL(0, bad.removeKnot(0.25, 2, 1));
}

} // SUITE
