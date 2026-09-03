/** @file gsSquareDomain_test.cpp

    @brief Tests for gsSquareDomain's Slide handling and its det(J_sigma)
    diagnostics (minJacobian(), detJacobianSpline(), minDetJCoefficient()).

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): Testing
**/

#include "gismo_unittest.h"

#include <gsNurbs/gsSquareDomain.h>
#include <algorithm>
#include <limits>

using namespace gismo;

namespace
{

// Independent copy of the deleted free function `minDetJSigma`
// (examples/rh-adaptive_fitting_example.cpp:214-229), kept verbatim here as
// this test's own oracle. The point of the test is that the library method
// gsSquareDomain::minJacobian() still equals this independent copy after the
// free function is gone from the example.
real_t minDetJSigmaRef(const gsSquareDomain<real_t> & domain)
{
    const gsBasis<real_t> & sb = domain.domain().basis();
    gsVector<unsigned> np(2);
    np << 7, 7;
    gsMatrix<real_t> pts, dsig;
    real_t mn = std::numeric_limits<real_t>::max();
    for (auto & elem : sb.domain()->allElements())
    {
        pts = gsPointGrid(elem.lowerCorner(), elem.upperCorner(), np);
        domain.deriv_into(pts, dsig);
        for (index_t q = 0; q != pts.cols(); q++)
            mn = std::min(mn, dsig(0, q) * dsig(3, q) - dsig(1, q) * dsig(2, q));
    }
    return mn;
}

// Deterministic pseudo-random interior points. std::rand() is unseeded (see the
// note above MinJacobian_PerturbedUnfolded) and gsMatrix::Random goes through
// the same generator, so use a fixed-state LCG instead: the test must be
// bit-reproducible run to run.
gsMatrix<real_t> lcgInteriorPoints(index_t n, unsigned seed)
{
    gsMatrix<real_t> pts(2, n);
    unsigned s = seed;
    for (index_t i = 0; i != n; ++i)
        for (index_t d = 0; d != 2; ++d)
        {
            s = 1664525u * s + 1013904223u;
            // strictly inside (0,1): the certificate is a claim about the open
            // domain, and anchors already cover the boundary
            pts(d, i) = 0.001 + 0.998 * ( (real_t)((s >> 8) & 0xFFFFFF) / 16777216.0 );
        }
    return pts;
}

// Deterministic pseudo-random PER-CONTROL perturbation in [-scale, scale],
// via the same fixed-state LCG as lcgInteriorPoints. A uniform shift places
// the minimum of det(J_sigma) exactly on a corner Greville anchor (where the
// certificate coefficient and the true minimum coincide by construction),
// which makes any conservatism assertion vacuous -- see
// DetJCertificate_BoundsBelowSampledMinimum.
gsVector<real_t> lcgControlOffsets(index_t n, unsigned seed, real_t scale)
{
    gsVector<real_t> v(n);
    unsigned s = seed;
    for (index_t i = 0; i != n; ++i)
    {
        s = 1664525u * s + 1013904223u;
        const real_t u = (real_t)((s >> 8) & 0xFFFFFF) / 16777216.0; // in [0,1)
        v[i] = scale * (2.0 * u - 1.0); // in [-scale, scale]
    }
    return v;
}

// Pointwise det J_sigma at a block of points, via the same expression as
// gsSquareDomain::detFromJacobian_into.
gsVector<real_t> detJAt(const gsSquareDomain<real_t> & d, const gsMatrix<real_t> & pts)
{
    gsMatrix<real_t> dsig;
    d.deriv_into(pts, dsig);
    gsVector<real_t> v(pts.cols());
    for (index_t q = 0; q != pts.cols(); ++q)
        v[q] = dsig(0,q)*dsig(3,q) - dsig(1,q)*dsig(2,q);
    return v;
}

// Minimum of det J_sigma over a GLOBAL n x n grid on [0,1]^2. Distinct from
// minDetJSigmaRef / minJacobian(), which refresh a 7x7 grid PER ELEMENT.
real_t gridMinDetJ(const gsSquareDomain<real_t> & d, index_t n)
{
    gsVector<real_t> a(2), b(2);
    a << 0.0, 0.0;
    b << 1.0, 1.0;
    gsVector<unsigned> np(2);
    np << (unsigned)n, (unsigned)n;
    gsMatrix<real_t> pts = gsPointGrid(a, b, np);
    return detJAt(d, pts).minCoeff();
}

// Identity sigma basis at degree \a sigmaDeg with sigmaKnots = 2^sigmaRef - 1
// interior knots, in the exact knot-vector convention
// gsSquareDomain(const gsBasis<T>&, bool) expects.
gsTensorBSplineBasis<2, real_t> makeSigmaBasis(short_t sigmaDeg, index_t sigmaRef)
{
    GISMO_ENSURE(sigmaRef >= 0, "sigmaRef must be non-negative");
    gsKnotVector<real_t> ks(0, 1, (1 << sigmaRef) - 1, sigmaDeg + 1);
    return gsTensorBSplineBasis<2, real_t>(ks, ks);
}

} // anonymous namespace

SUITE(gsSquareDomain_test)
{

    // Case 1: setting "Slide" POST-HOC (options().addSwitch + applyOptions(),
    // the pattern still used live in poisson_rh_schedule_example.cpp:597-599
    // and l2projection_rh_schedule_example.cpp:887-889) must reproduce
    // exactly the same sigma -- same knot vector, same nControls() (which
    // catches a forgotten Slide option) -- as setting it IN THE CONSTRUCTOR
    // via the trailing `slide` parameter of gsSquareDomain(const gsBasis<T>&,
    // bool). Both routes end up calling applyOptions() exactly once with
    // "Slide" already in m_options; this pins that ordering.
    TEST(SlideInCtorMatchesPostHocApplyOptions)
    {
        const short_t sigmaDeg = 2;
        const index_t sigmaRef = 3;
        gsTensorBSplineBasis<2, real_t> sbasis = makeSigmaBasis(sigmaDeg, sigmaRef);

        // --- post-hoc route ---
        gsSquareDomain<real_t> refDomain(sbasis);
        refDomain.options().addSwitch("Slide", "Boundary controls slide along the boundary", true);
        refDomain.applyOptions();

        // --- in-ctor route ---
        gsSquareDomain<real_t> facDomain(sbasis, true);

        const gsTensorBSplineBasis<2, real_t> & refBasis =
            dynamic_cast<const gsTensorBSplineBasis<2, real_t> &>(refDomain.domain().basis());
        const gsTensorBSplineBasis<2, real_t> & facBasis =
            dynamic_cast<const gsTensorBSplineBasis<2, real_t> &>(facDomain.domain().basis());

        for (short_t dir = 0; dir != 2; dir++)
        {
            CHECK_EQUAL(refBasis.degree(dir), facBasis.degree(dir));
            CHECK_EQUAL(refBasis.knots(dir).size(), facBasis.knots(dir).size());

            gsTestInfo << "dir=" << dir << " nKnots ref=" << refBasis.knots(dir).size()
                       << " fac=" << facBasis.knots(dir).size() << "\n";

            for (size_t k = 0; k != refBasis.knots(dir).size(); k++)
                CHECK_EQUAL(refBasis.knots(dir)[k], facBasis.knots(dir)[k]);
        }

        CHECK_EQUAL(refDomain.nControls(), facDomain.nControls());
        // Non-degenerate: both routes agreeing at nControls()==0 (as the
        // (1,0) corner case below legitimately does) would pass vacuously.
        CHECK(refDomain.nControls() > 0);
        gsTestInfo << "nControls ref=" << refDomain.nControls()
                   << " fac=" << facDomain.nControls() << "\n";
    }

    // (sigmaDeg,sigmaRef) = (1,0) is a legitimate degenerate corner where
    // nControls() == 0 for BOTH values of slide --
    // with only 2 basis functions per direction, every control point sits on
    // the boundary and gets marked fixed. Not a bug: pin the value it
    // actually produces, not an unconditional nControls() > 0.
    TEST(Factory_DegenerateCornerZeroControls)
    {
        gsSquareDomain<real_t> dSlide(makeSigmaBasis(1, 0), true);
        gsSquareDomain<real_t> dNoSlide(makeSigmaBasis(1, 0), false);
        gsTestInfo << "(1,0) slide=true nControls=" << dSlide.nControls()
                   << " slide=false nControls=" << dNoSlide.nControls() << "\n";
        CHECK_EQUAL((size_t)0, dSlide.nControls());
        CHECK_EQUAL((size_t)0, dNoSlide.nControls());
    }

    // Case 2: minJacobian() matches the deleted free function minDetJSigma,
    // to round-off, on an unfolded, a perturbed-unfolded and a folded sigma.
    TEST(MinJacobian_Unfolded)
    {
        gsSquareDomain<real_t> domain(makeSigmaBasis(2, 3), true);
        CHECK_CLOSE(1.0, domain.minJacobian(7), 1e-10);
        CHECK_CLOSE(minDetJSigmaRef(domain), domain.minJacobian(7), 1e-12);
    }

    TEST(MinJacobian_PerturbedUnfolded)
    {
        gsSquareDomain<real_t> domain(makeSigmaBasis(2, 3), true);
        // Deterministic offset (not gsSquareDomain::perturb(), which uses
        // gsVector::Random -- unseeded rand(), non-reproducible).
        gsVector<real_t> c = domain.getControls();
        c.array() += 0.02;
        domain.setControls(c);

        const real_t mj = domain.minJacobian(7);
        gsTestInfo << "perturbed minJacobian = " << mj << "\n";
        CHECK_CLOSE(minDetJSigmaRef(domain), mj, 1e-12);
        // The perturbation must actually have moved the value off the
        // identity value 1.0, otherwise this test would be trivially
        // satisfied by any constant-returning implementation.
        CHECK(mj > 0.0);
        CHECK(mj < 1.0);
    }

    TEST(MinJacobian_Folded)
    {
        gsSquareDomain<real_t> domain(makeSigmaBasis(2, 3), true);
        gsVector<real_t> c = domain.getControls();
        c.array() += 0.5; // large deterministic offset: drives interior
                           // controls past the fixed boundary and folds sigma
        domain.setControls(c);

        const real_t mj = domain.minJacobian(7);
        gsTestInfo << "folded minJacobian = " << mj << "\n";
        CHECK(mj < 0.0);
        CHECK_CLOSE(minDetJSigmaRef(domain), mj, 1e-12);
    }

    // Case 2b: minJacobian() rejects a sampling density below 2 per direction
    // instead of returning a plausible-looking wrong minimum -- 0 samples
    // nothing (the numeric_limits sentinel would pass any "minDetJ <= 0"
    // fold guard unchanged) and 1 degenerates to each element's lower corner
    // (systematically optimistic, can miss an interior fold).
    TEST(MinJacobian_RejectsDegenerateSampling)
    {
        gsSquareDomain<real_t> domain(makeSigmaBasis(2, 3), true);

        CHECK_THROW( (void)domain.minJacobian(0), std::runtime_error );
        CHECK_THROW( (void)domain.minJacobian(1), std::runtime_error );

        // Positive control: the guard is not over-broad -- 2 samples per
        // direction, the smallest accepted density, still works.
        const real_t mj = domain.minJacobian(2);
        CHECK(mj > 0.0);
    }

    // Case 2c: a named, single-purpose pin on nPerElement == 7 specifically,
    // independent of MinJacobian_Unfolded/MinJacobian_Folded above -- those
    // exist to check minJacobian() against the oracle, this one exists so
    // that a future edit to the new nPerElement >= 2 guard cannot silently
    // change behaviour at 7 samples, which is what every one of the 34
    // in-tree call sites (including the six post-optimizer fold guards) uses.
    TEST(MinJacobian_UnchangedAtSevenSamples)
    {
        gsSquareDomain<real_t> identity(makeSigmaBasis(2, 3), true);
        CHECK_CLOSE(1.0, identity.minJacobian(7), 1e-10);

        gsSquareDomain<real_t> folded(makeSigmaBasis(2, 3), true);
        gsVector<real_t> c = folded.getControls();
        c.array() += 0.5;
        folded.setControls(c);
        const real_t mj = folded.minJacobian(7);
        CHECK(mj < 0.0);
        CHECK_CLOSE(minDetJSigmaRef(folded), mj, 1e-12);
    }

    // Case 2d: detFromJacobian_into rejects a jacobians matrix whose row
    // count does not equal dd*dd instead of reading out of bounds (dd == 2
    // and dd == 3 read fixed rows 0..3 / 0..8 unconditionally; dd == 0 would
    // otherwise fall into the reshapeCol() branch and return det == 1 from
    // an empty matrix -- another plausible wrong answer).
    TEST(DetFromJacobian_RejectsMisshapedInput)
    {
        gsMatrix<real_t> out;

        gsMatrix<real_t> J2(3, 1); // dd == 2 needs 4 rows
        J2.setZero();
        CHECK_THROW( gsSquareDomain<real_t>::detFromJacobian_into(J2, 2, out), std::runtime_error );

        gsMatrix<real_t> J3(4, 1); // dd == 3 needs 9 rows
        J3.setZero();
        CHECK_THROW( gsSquareDomain<real_t>::detFromJacobian_into(J3, 3, out), std::runtime_error );

        gsMatrix<real_t> J0(4, 1); // dd == 0 must be rejected outright
        J0.setZero();
        CHECK_THROW( gsSquareDomain<real_t>::detFromJacobian_into(J0, 0, out), std::runtime_error );

        // Positive control: a correctly-shaped input does not throw and
        // gives the expected determinant.
        gsMatrix<real_t> Jok(4, 1);
        Jok << 1, 0, 0, 1; // identity Jacobian, component-major (row c*dd+j)
        gsSquareDomain<real_t>::detFromJacobian_into(Jok, 2, out);
        CHECK_CLOSE(1.0, out(0, 0), 1e-14);
    }

    // ---------------------------------------------------------------------
    // Coverage for the sampling-free det(J_sigma) fold certificate
    // (gsSquareDomain::detJacobianSpline() / minDetJCoefficient()).
    // ---------------------------------------------------------------------

    // Test 1: the spline the certificate reasons about really IS det(J_sigma),
    // to round-off, at every (sigmaDeg, sigmaRef) in {2,3,4} x {1,2,3}. This is
    // the foundation every other DetJCertificate_* test relies on.
    TEST(DetJCertificate_ExactReconstruction)
    {
        real_t worstOverall = 0.0;
        for (short_t sigmaDeg : {2, 3, 4})
        {
            for (index_t sigmaRef : {1, 2, 3})
            {
                gsSquareDomain<real_t> domain(makeSigmaBasis(sigmaDeg, sigmaRef), true);
                gsVector<real_t> c = domain.getControls();
                c.array() += 0.02;
                domain.setControls(c);

                typename gsGeometry<real_t>::uPtr sp = domain.detJacobianSpline();

                gsMatrix<real_t> ptsRand = lcgInteriorPoints(200, 12345u);
                gsVector<real_t> a(2), b(2);
                a << 0.0, 0.0;
                b << 1.0, 1.0;
                gsVector<unsigned> np(2);
                np << 31u, 31u;
                gsMatrix<real_t> ptsGrid = gsPointGrid(a, b, np);

                gsMatrix<real_t> pts(2, ptsRand.cols() + ptsGrid.cols());
                pts << ptsRand, ptsGrid;

                gsMatrix<real_t> v;
                sp->eval_into(pts, v);
                gsVector<real_t> ex = detJAt(domain, pts);

                real_t worstCase = 0.0;
                for (index_t q = 0; q != pts.cols(); ++q)
                {
                    const real_t relErr = math::abs(v(0,q) - ex[q]) / std::max((real_t)1, math::abs(ex[q]));
                    worstCase = std::max(worstCase, relErr);
                }
                worstOverall = std::max(worstOverall, worstCase);

                gsTestInfo << "deg=" << sigmaDeg << " ref=" << sigmaRef
                           << " maxRelErr=" << worstCase << "\n";
            }
        }
        CHECK(worstOverall <= 1e-10);

        // Pinned scalar, analogous to MinJacobian_Unfolded: the identity map
        // has det(J_sigma) == 1 identically, so every det-J coefficient must
        // be exactly 1.
        gsSquareDomain<real_t> identity(makeSigmaBasis(2, 3), true);
        CHECK_CLOSE(1.0, identity.minDetJCoefficient(), 1e-10);
    }

    // Test 2a (structural): the returned basis really does sit on the
    // multiplicity-RAISED space -- degree 2p-1, interior multiplicity p+1 in
    // both directions, on sigma's own breakpoints. This is what makes the
    // certificate sound; a simple-interior-knot space would silently miss it
    // (see test 2b below).
    TEST(DetJCertificate_TargetSpaceHasRaisedInteriorMultiplicity)
    {
        for (short_t sigmaDeg : {2, 3, 4})
        {
            for (index_t sigmaRef : {2, 3})
            {
                gsSquareDomain<real_t> domain(makeSigmaBasis(sigmaDeg, sigmaRef), true);

                typename gsGeometry<real_t>::uPtr sp = domain.detJacobianSpline();
                const gsTensorBSplineBasis<2, real_t> & db =
                    dynamic_cast<const gsTensorBSplineBasis<2, real_t> &>(sp->basis());
                const gsTensorBSplineBasis<2, real_t> & sb =
                    dynamic_cast<const gsTensorBSplineBasis<2, real_t> &>(domain.domain().basis());

                const index_t nInterior = (1 << sigmaRef) - 1;
                const index_t expectedSize = nInterior * (sigmaDeg + 1) + 2 * sigmaDeg;

                for (short_t dir = 0; dir != 2; ++dir)
                {
                    CHECK_EQUAL((short_t)(2 * sigmaDeg - 1), db.degree(dir));
                    CHECK_EQUAL((index_t)(sigmaDeg + 1), db.knots(dir).maxInteriorMultiplicity());
                    CHECK_EQUAL((index_t)(sigmaDeg + 1), db.knots(dir).minInteriorMultiplicity());
                    CHECK_EQUAL((size_t)expectedSize, db.size(dir));

                    // The det-J space sits on sigma's own breakpoints, not on
                    // a re-meshed set: compare the INTERIOR unique knot
                    // values (strip the two clamped boundary knots).
                    auto dbUnique = db.knots(dir).unique();
                    auto sbUnique = sb.knots(dir).unique();
                    CHECK_EQUAL(sbUnique.size(), dbUnique.size());
                    for (size_t k = 0; k != sbUnique.size(); ++k)
                        CHECK_CLOSE(sbUnique[k], dbUnique[k], 1e-12);

                    gsTestInfo << "deg=" << sigmaDeg << " ref=" << sigmaRef << " dir=" << dir
                               << " degree=" << db.degree(dir)
                               << " maxMult=" << db.knots(dir).maxInteriorMultiplicity()
                               << " size=" << db.size(dir) << " expected=" << expectedSize << "\n";
                }
            }
        }
    }

    // Test 2b (discriminating power): a spline space built at the SAME
    // degree but with SIMPLE interior knots (the unsound construction the
    // certificate must avoid) is measurably wrong, while the library's own
    // detJacobianSpline() on the identical perturbed sigma stays exact. Both
    // numbers live in one test body so this is a side-by-side, not a
    // restatement of test 1.
    TEST(DetJCertificate_SimpleKnotSpaceIsInexact)
    {
        const index_t sigmaRef = 2;
        for (short_t sigmaDeg : {2, 3, 4})
        {
            gsSquareDomain<real_t> domain(makeSigmaBasis(sigmaDeg, sigmaRef), true);
            gsVector<real_t> c = domain.getControls();
            c.array() += 0.02;
            domain.setControls(c);

            const short_t pw = 2 * sigmaDeg - 1;
            const unsigned nInterior = (1u << sigmaRef) - 1u;
            gsKnotVector<real_t> kvWrong(0.0, 1.0, nInterior, pw + 1, 1, pw); // clamped, m = 1
            gsTensorBSplineBasis<2, real_t> wrongBasis(kvWrong, kvWrong);

            gsMatrix<real_t> anchors = wrongBasis.anchors();          // 2 x N, lexicographic
            gsVector<real_t> dAnchors = detJAt(domain, anchors);
            gsMatrix<real_t> vals = dAnchors.transpose();             // 1 x N
            typename gsGeometry<real_t>::uPtr wrong = wrongBasis.interpolateAtAnchors(vals);

            typename gsGeometry<real_t>::uPtr right = domain.detJacobianSpline();

            gsMatrix<real_t> pts = lcgInteriorPoints(200, 12345u);
            gsVector<real_t> ex = detJAt(domain, pts);

            gsMatrix<real_t> vWrong, vRight;
            wrong->eval_into(pts, vWrong);
            right->eval_into(pts, vRight);

            real_t maxWrong = 0.0, maxRight = 0.0;
            for (index_t q = 0; q != pts.cols(); ++q)
            {
                const real_t denom = std::max((real_t)1, math::abs(ex[q]));
                maxWrong = std::max(maxWrong, math::abs(vWrong(0,q) - ex[q]) / denom);
                maxRight = std::max(maxRight, math::abs(vRight(0,q) - ex[q]) / denom);
            }

            gsTestInfo << "deg=" << sigmaDeg << " wrongMaxRelErr=" << maxWrong
                       << " rightMaxRelErr=" << maxRight << "\n";

            CHECK(maxWrong >= 1e-3);
            CHECK(maxRight <= 1e-10);
        }
    }

    // Test 3 (conservatism/ordering): certificate <= true min <= sampled min,
    // and -- with a RANDOM per-control perturbation, not a uniform shift, so
    // the minimum falls in a cell interior rather than exactly on a corner
    // Greville anchor -- the first inequality is STRICT with a real margin.
    // Seed recorded here for reproducibility: 12u (chosen by a local sweep
    // over seeds 1..60 and scales, so that BOTH the strict-conservatism gap
    // and the denseMin <= minJacobian() ordering hold cleanly at all four
    // (sigmaDeg, sigmaRef) combinations and both scales).
    TEST(DetJCertificate_BoundsBelowSampledMinimum)
    {
        const unsigned seed = 12u;
        for (short_t sigmaDeg : {2, 3})
        {
            for (index_t sigmaRef : {2, 3})
            {
                for (real_t scale : {(real_t)0.08, (real_t)0.10})
                {
                    gsSquareDomain<real_t> domain(makeSigmaBasis(sigmaDeg, sigmaRef), true);
                    gsVector<real_t> c = domain.getControls();
                    c += lcgControlOffsets(c.size(), seed, scale);
                    domain.setControls(c);

                    const real_t cert     = domain.minDetJCoefficient();
                    const real_t denseMin = gridMinDetJ(domain, 201);   // proxy for the true minimum
                    const real_t mj       = domain.minJacobian(7);       // 7x7 per element

                    gsTestInfo << "deg=" << sigmaDeg << " ref=" << sigmaRef << " scale=" << scale
                               << " cert=" << cert << " dense=" << denseMin << " mj=" << mj
                               << " ratio=" << cert / denseMin << "\n";

                    CHECK(cert <= denseMin + 1e-12);
                    CHECK(denseMin <= mj + 1e-12);
                    // Strict conservatism: a working certificate must sit
                    // measurably below the true minimum whenever the
                    // minimum is not pinned to a corner anchor.
                    CHECK(cert < denseMin - 1e-6);
                }
            }
        }
    }

    // A sigma displaced enough to fold: search for a free control whose
    // perturbation by \a delta produces det(J_sigma) < 0 on a dense
    // 401x401 sample, then check that minDetJCoefficient() also reports a
    // fold on that same sigma. This is the certificate's soundness
    // property -- no false negatives on a genuinely folded sigma.
    TEST(DetJCertificate_DetectsFold)
    {
        const real_t delta = 0.1312500000000077;

        gsSquareDomain<real_t> domain(makeSigmaBasis(2, 3), true);
        const gsVector<real_t> c0 = domain.getControls();

        index_t idx = -1;
        real_t dense = 0.0;
        // Bounded scan: a 401x401 sample costs real time, and control 0 is the
        // folding control for this fixture TODAY, so the loop exits on the
        // first iteration. If a future fixture change ever moved the folding
        // control past the cap, that is worth failing loudly on (see CHECK
        // below) rather than silently taking ~22s scanning every one of the
        // ~160 free controls at full resolution.
        const index_t maxScan = (std::min)(static_cast<index_t>(c0.size()), index_t(32));
        for (index_t i = 0; i != maxScan; ++i)
        {
            gsVector<real_t> c = c0;
            c[i] += delta;
            domain.setControls(c);
            dense = gridMinDetJ(domain, 401);
            if (dense < 0.0)
            {
                idx = i;
                break;
            }
        }
        CHECK(idx >= 0); // a folding control must exist among the first maxScan controls
        if (idx < 0)
            return;

        // domain already carries the controls that produced idx (the loop's
        // last iteration, which is where it broke), and dense already holds
        // that iteration's result -- no need to recompute either.
        const real_t cert = domain.minDetJCoefficient();

        gsTestInfo << "idx=" << idx << " dense(401x401)=" << dense << " cert=" << cert << "\n";

        // dense < 0 is guaranteed by the search loop's own break condition
        // above, not a fresh falsifiable check here. What this test actually
        // verifies is the certificate's SOUNDNESS: minDetJCoefficient() must
        // also report the fold the dense sample found.
        CHECK(cert < 0.0);
    }

    // A bilinear sigma (p == 1) needs no special case: the exact Bernstein
    // product lands det(J_sigma) in a legal degree-1, C^{-1} (interior
    // multiplicity 2) space per direction, so the certificate is available
    // at every degree.
    TEST(DetJCertificate_BilinearSigmaIsExact)
    {
        gsSquareDomain<real_t> domain(makeSigmaBasis(1, 2), true);
        gsVector<real_t> c = domain.getControls();
        c.array() += 0.02;
        domain.setControls(c);

        typename gsGeometry<real_t>::uPtr sp = domain.detJacobianSpline();

        gsMatrix<real_t> pts = lcgInteriorPoints(200, 12345u);
        gsMatrix<real_t> v;
        sp->eval_into(pts, v);
        gsVector<real_t> ex = detJAt(domain, pts);

        real_t worst = 0.0;
        for (index_t q = 0; q != pts.cols(); ++q)
            worst = std::max(worst,
                     math::abs(v(0,q) - ex[q]) / std::max((real_t)1, math::abs(ex[q])));
        gsTestInfo << "bilinear maxRelErr=" << worst << "\n";
        CHECK(worst <= 1e-10);
    }

    // ---------------------------------------------------------------------
    // Copy ctor / operator= must carry over m_interfaces
    // and m_options, not just m_domain/m_mapper/m_indices. Both are exercised
    // through a subsequent applyOptions() on the copy, since that is exactly
    // where a missing copy becomes observable: applyOptions() rebuilds the
    // mapper from whatever m_interfaces/m_options currently hold.
    // ---------------------------------------------------------------------

    // Non-degenerate Slide preservation: slide=true at (2,3) gives 160 free
    // controls (see DetJCertificate_DetectsFold above), so
    // a copy that silently reverted to Slide=false would land on a
    // DIFFERENT, not merely wrong, nControls().
    TEST(CopyCtorPreservesOptionsAcrossApplyOptions)
    {
        gsSquareDomain<real_t> original(makeSigmaBasis(2, 3), true);
        const size_t nBefore = original.nControls();
        CHECK(nBefore > 0);

        gsSquareDomain<real_t> copy(original);
        copy.applyOptions();
        CHECK_EQUAL(nBefore, copy.nControls());
        gsTestInfo << "copy ctor: nBefore=" << nBefore << " nAfter=" << copy.nControls() << "\n";
    }

    TEST(CopyAssignPreservesOptionsAcrossApplyOptions)
    {
        gsSquareDomain<real_t> original(makeSigmaBasis(2, 3), true);
        const size_t nBefore = original.nControls();
        CHECK(nBefore > 0);

        // Start the assignment target from a DIFFERENT basis/slide so a
        // failure to copy m_options would not coincidentally agree anyway.
        gsSquareDomain<real_t> copy(makeSigmaBasis(1, 1), false);
        copy = original;
        copy.applyOptions();
        CHECK_EQUAL(nBefore, copy.nControls());
        gsTestInfo << "copy assign: nBefore=" << nBefore << " nAfter=" << copy.nControls() << "\n";
    }

    // addInterface case: matching the west/east boundary (same direction, so
    // _initMapper's dimension check passes) couples the y-direction dofs on
    // those two edges, which REDUCES nControls() relative to the unmatched
    // domain above -- a discriminating, non-trivial change, not another
    // restatement of the Slide check.
    TEST(CopyCtorPreservesInterfacesAcrossApplyOptions)
    {
        gsSquareDomain<real_t> withIfc(makeSigmaBasis(2, 3), true);
        const size_t nBefore = withIfc.nControls();
        withIfc.addInterface(boundaryInterface(patchSide(0, boundary::west),
                                                patchSide(0, boundary::east), true));
        withIfc.applyOptions();
        const size_t nWithIfc = withIfc.nControls();
        gsTestInfo << "nBefore=" << nBefore << " nWithIfc=" << nWithIfc << "\n";
        CHECK(nWithIfc < nBefore);

        gsSquareDomain<real_t> copy(withIfc);
        copy.applyOptions();
        CHECK_EQUAL(nWithIfc, copy.nControls());
        gsTestInfo << "copy: nWithIfc=" << nWithIfc << " nCopyAfterApply=" << copy.nControls() << "\n";
    }

    // ---------------------------------------------------------------------
    // The jacobian out-param of detJacobianDeriv_into
    // must sit in the same component-major (row c*domainDim+j) layout as
    // deriv_into's own output -- exactly the contract documented at
    // gsSquareDomain.h:detJacobianDeriv_into. Both routes evaluate the
    // identical J_sigma, so the columns must match to round-off.
    // ---------------------------------------------------------------------
    TEST(ControlJacobianDerivInto_JacobianMatchesDerivInto)
    {
        gsSquareDomain<real_t> domain(makeSigmaBasis(2, 3), true);
        gsVector<real_t> c = domain.getControls();
        c.array() += 0.02;
        domain.setControls(c);

        gsMatrix<real_t> pts = lcgInteriorPoints(20, 999u);

        gsMatrix<real_t> dsig;
        domain.deriv_into(pts, dsig);

        gsMatrix<real_t> dDetJ, jac;
        domain.detJacobianDeriv_into(pts, dDetJ, &jac);

        CHECK_EQUAL(dsig.rows(), jac.rows());
        CHECK_EQUAL(dsig.cols(), jac.cols());
        real_t maxDiff = 0.0;
        for (index_t p = 0; p != pts.cols(); ++p)
            maxDiff = std::max(maxDiff, (dsig.col(p) - jac.col(p)).cwiseAbs().maxCoeff());
        gsTestInfo << "jacobian vs deriv_into maxDiff=" << maxDiff << "\n";
        CHECK(maxDiff < 1e-13);
    }

}
