/** @file gsMomentRule_test.cpp

    @brief Tests the core solve-free nodal moment fitting rule gsMomentRule.

    WARNING on what is and is not evidence here. The nodal Lagrange basis is
    evaluated in the *normalized* (second) barycentric form, i.e. it is
    self-normalising: sum_i L_i(x) == 1 for every x, in every coordinate frame
    and for every index order. Therefore
        sum_i w_i = sum_k omega_k sum_i L_i(x_k) = sum_k omega_k
    is an ALGEBRAIC IDENTITY: every "total area/volume/mass is preserved" check
    is NECESSARY BUT NOT SUFFICIENT -- it cannot distinguish a correct rule from
    one that attaches the correct weights to the WRONG POINTS (a broken affine
    x -> u map, a transposed lexicographic index, an off-by-one in the output
    order all reproduce the area exactly). Every such assertion below carries an
    inline necessary-not-sufficient comment.

    The correctness burden is carried by exactly two mechanisms:
      1. ReproducesGauss / AnisotropicOrders: node-for-node AND weight-for-weight
         identity with gsGaussRule on a box forced through the compression path
         (this locks the affine frame of the Lagrange evaluation);
      2. PolynomialExactness: per-element polynomial exactness on NON-CONSTANT,
         ASYMMETRIC integrands x^a y^b with a != b (this locks the index order).

    This suite is deliberately free of any gsAlgoim dependency: moment fitting
    is not algoim-specific and must compile and pass with that module disabled.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst
 **/

#include "gismo_unittest.h"       // Brings in G+Smo and the UnitTest++ framework

#include <gsDomain/gsImplicitTrimmedDomain.h>

namespace {

/// x^a, by repeated multiplication (an exact polynomial, unlike std::pow).
inline real_t mono(real_t x, index_t a)
{
    real_t v = 1;
    for (index_t i = 0; i < a; ++i) v *= x;
    return v;
}

/// Machine-precision level of real_t, so tolerances do not hardcode double.
inline real_t mrEps() { return std::numeric_limits<real_t>::epsilon(); }

/// Node/weight comparison tolerance, scaled with real_t.
inline real_t mrTol() { return math::max((real_t)1e-13, (real_t)100 * mrEps()); }

/// The unit circle level set (negative inside).
inline gsFunctionExpr<real_t> circle2D()
{ return gsFunctionExpr<real_t>("sqrt(x^2+y^2)-1", 2); }

/// A DELIBERATELY non-square bounding box around the unit circle, in the layout
/// gsTrimmedDomain::init(bbox,...) expects: a 2 x d matrix whose ROW 0 holds the
/// lower corners and ROW 1 the upper corners (NOT the d x 2 layout of
/// gsFunction::support()); see gsTrimmedDomain.h:446-448.
inline gsMatrix<real_t> circleBBox()
{
    gsMatrix<real_t> bbox(2,2);
    bbox << -1.2, -1.0,     // lower corners
             1.2,  1.0;     // upper corners
    return bbox;
}

/// A DELIBERATELY non-symmetric background cell grid (5 x 4): a transposed
/// index or a collapse to a single direction count cannot hide behind symmetry.
inline gsVector<index_t,2> circleCells()
{
    gsVector<index_t,2> nc;
    nc << 5, 4;
    return nc;
}

inline gsVector<index_t> constNodes(index_t d, index_t n)
{ return gsVector<index_t>::Constant(d, n); }

} // anonymous namespace

SUITE(gsMomentRule_test)          // The suite should have the same name as the file
{

/// TEST 1 -- THE FRAME LOCK.
///
/// A plain gsGaussRule with 6 nodes/direction (36 points) is wrapped and
/// compressed onto a 3x3 output grid. 36 > 9, so the pass-through guard does
/// NOT fire and the Lagrange accumulation really runs (asserted through Stats).
///
/// 6 is chosen so that NO underlying point coincides with an output node: the
/// 6-point Gauss abscissae contain neither 0 nor +-sqrt(3/5), which are exactly
/// the 3-point abscissae. This forces _lagrange1d's barycentric quotient -- the
/// part that depends on the affine map x -> u -- to execute for every point. Do
/// NOT lower this to 5: the 5-point rule shares the centre node with the
/// 3-point rule and the "hit" shortcut would fire there.
///
/// Since the underlying rule covers the whole box and is exact to degree 11 >=
/// deg(L_i) = 2, we have w_i = sum_k omega_k L_i(x_k) = int_box L_i = w_i^Gauss
/// exactly (Gauss is interpolatory), so the compressed rule must reproduce
/// gsGaussRule node-by-node AND weight-by-weight.
///
/// DETECTS: a wrong affine x -> u map, a wrong output grid, a wrong barycentric
/// weight, a wrong element frame.
/// DOES NOT DETECT: a transposed lexicographic index. Tensor-Gauss weights
/// factorise, w_ij = hx*hy*a_i*a_j, which is symmetric under (i,j) -> (j,i)
/// EVEN on a non-square box; if the transposition is applied consistently to
/// nodes and weights, the comparison is against a permuted-but-equal set. That
/// burden falls on TEST 2 (PolynomialExactness), via asymmetric monomials.
TEST(ReproducesGauss)
{
    const index_t n = 3;
    typename gsQuadRule<real_t>::uPtr inner(gsGaussRule<real_t>::make(constNodes(2,6)));
    gsMomentRule<real_t> mr(give(inner), constNodes(2,n));

    gsVector<real_t> lo(2), hi(2);
    lo << -0.30, -0.10;
    hi <<  0.10,  0.25;          // non-square, off-centre

    gsMatrix<real_t> nodes;
    gsVector<real_t> wts;
    mr.mapTo(lo, hi, nodes, wts);

    // The compression path really ran (no pass-through, no zero-column exit).
    CHECK_EQUAL(0, mr.stats().nPassThroughQPs);
    CHECK_EQUAL(0, mr.stats().nPassThroughElements);
    CHECK_EQUAL(n*n, mr.stats().nOutputQPs);
    CHECK_EQUAL(36,  mr.stats().nUnderlyingQPs);

    gsGaussRule<real_t> gauss(constNodes(2,n));
    gsMatrix<real_t> gnodes;
    gsVector<real_t> gwts;
    gauss.mapTo(lo, hi, gnodes, gwts);

    CHECK_EQUAL(n*n, nodes.cols());
    CHECK_EQUAL(gnodes.cols(), nodes.cols());
    CHECK_EQUAL(gwts.size(), wts.size());
    CHECK((nodes - gnodes).array().abs().maxCoeff() < mrTol());
    CHECK((wts   - gwts  ).array().abs().maxCoeff() < mrTol());
}

/// TEST 2 -- THE ACCURACY LOCK (and the index-order lock).
///
/// For any f that the output grid interpolates exactly (tensor degree n-1 per
/// direction),
///     sum_i f(g_i) w_i = sum_k omega_k sum_i f(g_i) L_i(x_k)
///                      = sum_k omega_k f(x_k),
/// i.e. the compressed rule integrates f EXACTLY as the underlying rule does,
/// independently of how accurate the underlying rule itself is. That identity
/// is what breaks under a wrong u-mapping or a permuted lexicographic index.
///
/// Two deliberate design points:
///  * the comparison is PER ELEMENT. The circle is symmetric in x <-> y, so a
///    transposed index order would cancel between mirror elements in a global
///    sum and go unnoticed.
///  * the exponent bound is written literally as n-1 per direction (n being the
///    output order) because moment fitting reproduces polynomials only to that
///    degree -- an off-by-one in the output order cannot hide inside a loose
///    loop bound. The pairs with a != b -- (2,0), (0,2), (2,1), (1,2) -- are
///    included by construction and are the transposition detector.
TEST(PolynomialExactness)
{
    const index_t n = 3;                    // output points per direction
    gsFunctionExpr<real_t> phi = circle2D();
    gsImplicitTrimmedDomain<2,real_t> domain(phi, circleBBox(), circleCells(), 5, 2);

    // Wrapped rule: octree cut-cell, 2 levels, over a 4-point Gauss leaf rule.
    typename gsQuadRule<real_t>::uPtr innerA(
        new gsOctreeCutCellRule<real_t>(gsGaussRule<real_t>::make(constNodes(2,4)), phi, 2));
    gsMomentRule<real_t> mr(give(innerA), constNodes(2,n));

    // An independent, identically configured copy to compare against.
    gsOctreeCutCellRule<real_t> und(gsGaussRule<real_t>::make(constNodes(2,4)), phi, 2);

    gsMatrix<real_t> mfPts, uPts;
    gsVector<real_t> mfW, uW;
    std::vector<real_t> globMf(n*n, (real_t)0), globUnd(n*n, (real_t)0);
    real_t worst = 0;                       // worst scaled per-element deviation
    index_t nCompressed = 0;

    for (auto & elem : domain.allElements())
    {
        mr .mapTo(elem.lowerCorner(), elem.upperCorner(), mfPts, mfW);
        und.mapTo(elem.lowerCorner(), elem.upperCorner(), uPts,  uW );

        // Comparing a rule with itself proves nothing: skip the elements the
        // pass-through guard forwarded verbatim (and the empty ones).
        const bool compressed = (uPts.cols() > n*n);
        if (compressed) ++nCompressed;

        // Bound written as n-1 per direction, NOT as a literal.
        for (index_t a = 0; a <= n-1; ++a)
        for (index_t b = 0; b <= n-1; ++b)
        {
            real_t Imf = 0, Iund = 0, scale = 0;
            for (index_t i = 0; i < mfW.size(); ++i)
            {
                const real_t f = mono(mfPts(0,i), a) * mono(mfPts(1,i), b);
                Imf   += mfW[i] * f;
                scale += math::abs(mfW[i] * f);
            }
            for (index_t k = 0; k < uW.size(); ++k)
            {
                const real_t f = mono(uPts(0,k), a) * mono(uPts(1,k), b);
                Iund  += uW[k] * f;
                scale += math::abs(uW[k] * f);
            }
            globMf [a*n+b] += Imf;
            globUnd[a*n+b] += Iund;
            if (compressed && scale > 0)   // cancellation-aware forward error
                worst = math::max(worst, math::abs(Imf - Iund) / scale);
        }
    }

    // The discriminating check: per element, asymmetric monomials included.
    CHECK(worst < (real_t)1e3 * mrEps());

    // Compression must actually have happened.
    CHECK(nCompressed > 0);
    CHECK(mr.stats().nOutputQPs > 0);
    CHECK_EQUAL(nCompressed * n * n, mr.stats().nOutputQPs);

    // Global identity over all monomials, including the constant a = b = 0.
    // NECESSARY, NOT SUFFICIENT for a != b entries too: mirror elements can
    // cancel a transposition error here. Kept as a cheap smoke test only; the
    // per-element check above is the load-bearing one.
    real_t globWorst = 0;
    for (index_t m = 0; m != n*n; ++m)
        globWorst = math::max(globWorst, math::abs(globMf[m] - globUnd[m])
                              / math::max((real_t)1, math::abs(globUnd[m])));
    CHECK(globWorst < math::max((real_t)1e-12, (real_t)1e4 * mrEps()));
}

/// TEST 3 -- ANISOTROPIC ORDERS, numNodes = (3,5).
///
/// Unlike TEST 1 this test DOES discriminate a transposed index, because a 3x5
/// grid and a 5x3 grid differ in shape: they place 15 nodes at DIFFERENT
/// positions, so the node-by-node comparison in part (a) reacts (a 5x3 grid
/// would also survive a collapse to numNodes.maxCoeff() = 25 or .minCoeff() = 9
/// only by changing the count). That is a bonus; TEST 2 remains the primary
/// defence against transposition.
///
/// Two parts, because they need DIFFERENT underlying rules:
///
/// (a) node/weight identity with gsGaussRule(3,5), wrapping a 6-point/direction
///     Gauss rule on non-square boxes (36 > 15, so the compression runs). Note
///     that with a full-box underlying rule the compressed weights collapse to
///     w_i = int L_i = w_i^Gauss, i.e. the compressed rule IS the (3,5) Gauss
///     rule -- which is exactly what makes it a reference, but also means such
///     a configuration can say nothing about the n-1 exactness bound (it would
///     be exact to 2n-1 there). Hence part (b).
///
/// (b) the per-direction exactness bound nx-1 / ny-1, measured where it is real:
///     around an octree cut-cell rule, whose point cloud is NOT a Gauss rule of
///     the element, so sum_i w_i f(g_i) = sum_k omega_k f(x_k) genuinely holds
///     only up to tensor degree (nx-1, ny-1). Compared PER ELEMENT.
TEST(AnisotropicOrders)
{
    const index_t nx = 3, ny = 5;
    gsVector<index_t> nn(2); nn << nx, ny;

    // ---- (a) node-for-node / weight-for-weight identity --------------------
    {
        typename gsQuadRule<real_t>::uPtr inner(gsGaussRule<real_t>::make(constNodes(2,6)));
        gsMomentRule<real_t> mr(give(inner), nn);

        // Several non-square, non-concentric boxes: the checks are per box.
        std::vector<std::pair<gsVector<real_t>,gsVector<real_t> > > boxes;
        {
            gsVector<real_t> lo(2), hi(2);
            lo << -0.30, -0.10;  hi << 0.10, 0.25;  boxes.push_back(std::make_pair(lo,hi));
            lo <<  0.50,  1.30;  hi << 1.70, 1.55;  boxes.push_back(std::make_pair(lo,hi));
            lo << -2.00, -0.05;  hi << -1.10, 0.90; boxes.push_back(std::make_pair(lo,hi));
        }

        gsMatrix<real_t> nodes, gnodes;
        gsVector<real_t> wts, gwts;
        gsGaussRule<real_t> gauss(nn);

        for (size_t e = 0; e != boxes.size(); ++e)
        {
            const gsVector<real_t> & lo = boxes[e].first;
            const gsVector<real_t> & hi = boxes[e].second;

            mr.mapTo(lo, hi, nodes, wts);

            // 3*5 -- not 5*3 (same count, different nodes), not 25, not 9.
            CHECK_EQUAL(nx*ny, nodes.cols());
            CHECK_EQUAL(nx*ny, wts.size());

            gauss.mapTo(lo, hi, gnodes, gwts);
            CHECK_EQUAL(gnodes.cols(), nodes.cols());
            CHECK((nodes - gnodes).array().abs().maxCoeff() < mrTol());
            CHECK((wts   - gwts  ).array().abs().maxCoeff() < mrTol());
        }

        CHECK_EQUAL(0, mr.stats().nPassThroughQPs);
        CHECK_EQUAL((index_t)boxes.size() * nx * ny, mr.stats().nOutputQPs);
    }

    // ---- (b) the per-direction exactness bound, on a genuine point cloud ---
    {
        gsFunctionExpr<real_t> phi = circle2D();
        gsImplicitTrimmedDomain<2,real_t> domain(phi, circleBBox(), circleCells(), 5, 2);

        typename gsQuadRule<real_t>::uPtr inner(
            new gsOctreeCutCellRule<real_t>(gsGaussRule<real_t>::make(constNodes(2,4)), phi, 2));
        gsMomentRule<real_t> mr(give(inner), nn);

        gsOctreeCutCellRule<real_t> und(gsGaussRule<real_t>::make(constNodes(2,4)), phi, 2);

        gsMatrix<real_t> mfPts, uPts;
        gsVector<real_t> mfW, uW;
        real_t worst = 0;
        index_t nCompressed = 0;

        for (auto & elem : domain.allElements())
        {
            mr .mapTo(elem.lowerCorner(), elem.upperCorner(), mfPts, mfW);
            und.mapTo(elem.lowerCorner(), elem.upperCorner(), uPts,  uW );

            if (uPts.cols() <= nx*ny)     // pass-through / empty: nothing to prove
                continue;
            ++nCompressed;
            CHECK_EQUAL(nx*ny, mfPts.cols());

            // Bounds written as nx-1 and ny-1, NOT as literals.
            for (index_t a = 0; a <= nx-1; ++a)
            for (index_t b = 0; b <= ny-1; ++b)
            {
                real_t Imf = 0, Iund = 0, scale = 0;
                for (index_t i = 0; i < mfW.size(); ++i)
                {
                    const real_t f = mono(mfPts(0,i), a) * mono(mfPts(1,i), b);
                    Imf   += mfW[i] * f;
                    scale += math::abs(mfW[i] * f);
                }
                for (index_t k = 0; k < uW.size(); ++k)
                {
                    const real_t f = mono(uPts(0,k), a) * mono(uPts(1,k), b);
                    Iund  += uW[k] * f;
                    scale += math::abs(uW[k] * f);
                }
                if (scale > 0)
                    worst = math::max(worst, math::abs(Imf - Iund) / scale);
            }
        }

        CHECK(nCompressed > 0);
        CHECK(worst < (real_t)1e3 * mrEps());
        CHECK_EQUAL(nCompressed * nx * ny, mr.stats().nOutputQPs);
    }
}

/// TEST 4 -- ALPHA, against a CLOSED-FORM oracle (no cut geometry involved).
///
/// alpha is the paper's multiplicative field INSIDE the integrand (not the FCM
/// fictitious-domain blend fraction of gsAlgoimMomentFittingRule -- these are
/// two distinct quantities and there is deliberately no range guard here).
///
/// On [0,1]^2 with alpha(x,y) = x + 2y (non-constant and ASYMMETRIC) and
/// f(x,y) = x^2 y:
///     int_0^1 int_0^1 (x+2y) x^2 y dx dy
///        = (int x^3)(int y) + 2 (int x^2)(int y^2)
///        = (1/4)(1/2) + 2(1/3)(1/3) = 1/8 + 2/9 = 25/72.
/// The compression reproduces f exactly (per-direction degree 2 and 1, both
/// <= n-1 = 3), so the result equals the underlying Gauss integral of alpha*f,
/// whose per-direction degree is <= 4 <= 2*8-1 and which is therefore itself
/// exact. Both sides are exact, so the oracle is the closed form.
TEST(AlphaWeight)
{
    const index_t n = 4;
    gsFunctionExpr<real_t> alpha("x+2*y", 2);

    gsVector<real_t> lo(2), hi(2);
    lo << 0.0, 0.0;
    hi << 1.0, 1.0;

    const real_t exact = (real_t)25 / (real_t)72;
    const real_t tol   = math::max((real_t)1e-13, (real_t)1e3 * mrEps());

    // (a) 8-point/direction Gauss wrapped: 64 > 16, the guard cannot fire.
    {
        typename gsQuadRule<real_t>::uPtr inner(gsGaussRule<real_t>::make(constNodes(2,8)));
        gsMomentRule<real_t> mr(give(inner), constNodes(2,n), &alpha);

        gsMatrix<real_t> nodes;
        gsVector<real_t> wts;
        mr.mapTo(lo, hi, nodes, wts);

        CHECK_EQUAL(n*n, nodes.cols());
        CHECK_EQUAL(0, mr.stats().nPassThroughQPs);

        real_t I = 0;
        for (index_t i = 0; i < wts.size(); ++i)
            I += wts[i] * mono(nodes(0,i),2) * mono(nodes(1,i),1);
        CHECK_CLOSE(exact, I, tol);
    }

    // (b) 2-point/direction Gauss wrapped: only 4 points, so 4 <= 16 WOULD
    // trigger the pass-through guard -- but the guard is conditioned on "no
    // alpha", because forwarding the points unchanged would silently drop the
    // alpha factor. The rule must still emit 16 compressed points, and the
    // closed form must still come out: alpha*f has per-direction degree <= 3,
    // which the 2-point Gauss (exact to 2*2-1 = 3) integrates exactly.
    {
        typename gsQuadRule<real_t>::uPtr inner(gsGaussRule<real_t>::make(constNodes(2,2)));
        gsMomentRule<real_t> mr(give(inner), constNodes(2,n), &alpha);

        gsMatrix<real_t> nodes;
        gsVector<real_t> wts;
        mr.mapTo(lo, hi, nodes, wts);

        CHECK_EQUAL(n*n, nodes.cols());
        CHECK_EQUAL(0, mr.stats().nPassThroughQPs);
        CHECK_EQUAL(0, mr.stats().nPassThroughElements);
        CHECK_EQUAL(n*n, mr.stats().nOutputQPs);

        real_t I = 0;
        for (index_t i = 0; i < wts.size(); ++i)
            I += wts[i] * mono(nodes(0,i),2) * mono(nodes(1,i),1);
        CHECK_CLOSE(exact, I, tol);
    }
}

/// TEST 5 -- THE CAPABILITY CLAIM: moment fitting around the CORE cut-cell
/// rules, with NO gsAlgoim in the build at all. This is only possible since the
/// rules take ownership of the wrapped rule through a uPtr: stored by value
/// they would be sliced down to the base affine map and the octree subdivision
/// would silently disappear.
///
/// The load-bearing assertion is NOT the area (see the file header) but the
/// per-element integral of the ASYMMETRIC monomial x^2 y against the underlying
/// rule's own value.
TEST(WrapsCoreCutCellRules)
{
    const index_t n = 3;
    gsFunctionExpr<real_t> phi = circle2D();
    gsImplicitTrimmedDomain<2,real_t> domain(phi, circleBBox(), circleCells(), 5, 2);

    for (index_t variant = 0; variant != 2; ++variant)
    {
        typename gsQuadRule<real_t>::uPtr inner;
        typename gsQuadRule<real_t>::uPtr undPtr;
        if (0 == variant)   // (a) plain cut-cell rule, 6 Gauss nodes/direction
        {
            inner.reset(new gsCutCellRule<real_t>(
                            gsGaussRule<real_t>::make(constNodes(2,6)), phi));
            undPtr.reset(new gsCutCellRule<real_t>(
                            gsGaussRule<real_t>::make(constNodes(2,6)), phi));
        }
        else                // (b) octree cut-cell rule, 2 levels, 4 nodes/dir
        {
            inner.reset(new gsOctreeCutCellRule<real_t>(
                            gsGaussRule<real_t>::make(constNodes(2,4)), phi, 2));
            undPtr.reset(new gsOctreeCutCellRule<real_t>(
                            gsGaussRule<real_t>::make(constNodes(2,4)), phi, 2));
        }

        gsMomentRule<real_t> mr(give(inner), constNodes(2,n));

        gsMatrix<real_t> mfPts, uPts;
        gsVector<real_t> mfW, uW;
        real_t worst = 0, areaMf = 0, areaUnd = 0;
        index_t nCompressed = 0;

        for (auto & elem : domain.allElements())
        {
            mr     .mapTo(elem.lowerCorner(), elem.upperCorner(), mfPts, mfW);
            undPtr->mapTo(elem.lowerCorner(), elem.upperCorner(), uPts,  uW );

            areaMf  += mfW.sum();
            areaUnd += uW.sum();

            if (uPts.cols() <= n*n)     // pass-through / empty: nothing to prove
                continue;
            ++nCompressed;

            CHECK_EQUAL(n*n, mfPts.cols());

            // The discriminating check: an ASYMMETRIC monomial, per element.
            real_t Imf = 0, Iund = 0, scale = 0;
            for (index_t i = 0; i < mfW.size(); ++i)
            {
                const real_t f = mono(mfPts(0,i),2) * mono(mfPts(1,i),1);
                Imf   += mfW[i] * f;
                scale += math::abs(mfW[i] * f);
            }
            for (index_t k = 0; k < uW.size(); ++k)
            {
                const real_t f = mono(uPts(0,k),2) * mono(uPts(1,k),1);
                Iund  += uW[k] * f;
                scale += math::abs(uW[k] * f);
            }
            if (scale > 0)
                worst = math::max(worst, math::abs(Imf - Iund) / scale);
        }

        CHECK(nCompressed > 0);
        CHECK(worst < (real_t)1e3 * mrEps());

        // SMOKE TEST ONLY -- NECESSARY, NOT SUFFICIENT. sum_i w_i == sum_k
        // omega_k is an algebraic identity of the Lagrange partition of unity,
        // so reproducing the underlying rule's area says nothing about where
        // the points went. The absolute bounds against pi are sanity bounds on
        // the UNDERLYING rule's own (staircase) accuracy, not accuracy claims:
        // (a) is a crude one-level staircase (measured area 3.13258, error
        // 9.0e-3), (b) has 2 octree levels (measured 3.14203, error 4.4e-4).
        // The bounds below sit ~5x/~10x above those measurements so that a
        // different real_t or platform cannot make them flap; do not read them
        // as convergence statements.
        CHECK_CLOSE(areaUnd, areaMf, math::max((real_t)1e-12, (real_t)1e4*mrEps()));
        CHECK_CLOSE(EIGEN_PI, areaMf, 0 == variant ? (real_t)5e-2 : (real_t)5e-3);
    }
}

/// TEST 6 -- THE ZERO-COLUMN CONTRACT.
///
/// gsOctreeCutCellRule returns immediately (0 columns) for a box lying entirely
/// outside the level set. gsMomentRule must forward that: emit NOTHING.
///
/// Why this matters: gsExprAssembler skips an element whose points() has zero
/// columns. Emitting n^d zero-weight points instead would silently multiply the
/// assembly cost of every exterior element in the bounding box.
TEST(ZeroColumnContract)
{
    gsFunctionExpr<real_t> phi = circle2D();

    gsVector<real_t> lo(2), hi(2);
    lo << 1.2, 1.2;             // entirely outside the unit circle
    hi << 1.5, 1.5;

    // Premise: the wrapped rule really does emit nothing here.
    {
        gsOctreeCutCellRule<real_t> oc(gsGaussRule<real_t>::make(constNodes(2,4)), phi, 2);
        gsMatrix<real_t> uPts;
        gsVector<real_t> uW;
        oc.mapTo(lo, hi, uPts, uW);
        CHECK_EQUAL(0, uPts.cols());
        CHECK_EQUAL(0, uW.size());
    }

    typename gsQuadRule<real_t>::uPtr inner(
        new gsOctreeCutCellRule<real_t>(gsGaussRule<real_t>::make(constNodes(2,4)), phi, 2));
    gsMomentRule<real_t> mr(give(inner), constNodes(2,3));

    gsMatrix<real_t> nodes;
    gsVector<real_t> wts;
    mr.mapTo(lo, hi, nodes, wts);

    CHECK_EQUAL(0, nodes.cols());
    CHECK_EQUAL(0, wts.size());
    CHECK_EQUAL(1, mr.stats().nElements);
    CHECK_EQUAL(0, mr.stats().nOutputQPs);
    CHECK_EQUAL(0, mr.stats().nPassThroughQPs);
}

/// TEST 7 -- THE PASS-THROUGH GUARD.
///
/// Compressing FEWER points onto MORE points only adds interpolation error, so
/// a wrapped rule delivering at most numNodes.prod() points is forwarded
/// verbatim (when no alpha is set). This guard is what makes uncut elements
/// free without any inside/outside classification at all.
///
/// "Unchanged" is asserted with ZERO tolerance: a tolerance would let an actual
/// (harmless-looking but wrong) compression sneak through.
TEST(PassThroughGuard)
{
    typename gsQuadRule<real_t>::uPtr inner(gsGaussRule<real_t>::make(constNodes(2,2)));
    gsMomentRule<real_t> mr(give(inner), constNodes(2,3));   // 9 >= 4

    gsGaussRule<real_t> und(constNodes(2,2));

    gsVector<real_t> lo(2), hi(2);
    lo << -0.30, -0.10;
    hi <<  0.10,  0.25;          // non-square

    gsMatrix<real_t> nodes, gnodes;
    gsVector<real_t> wts, gwts;
    mr .mapTo(lo, hi, nodes, wts);
    und.mapTo(lo, hi, gnodes, gwts);

    CHECK_EQUAL(4, nodes.cols());
    CHECK_EQUAL(gnodes.cols(), nodes.cols());
    CHECK_EQUAL(gwts.size(), wts.size());
    CHECK_EQUAL((real_t)0, (nodes - gnodes).array().abs().maxCoeff());
    CHECK_EQUAL((real_t)0, (wts   - gwts  ).array().abs().maxCoeff());

    CHECK_EQUAL(4, mr.stats().nPassThroughQPs);
    CHECK_EQUAL(1, mr.stats().nPassThroughElements);
    CHECK_EQUAL(0, mr.stats().nOutputQPs);
}

/// TEST 8 -- THE EXACTNESS ORDER.
///
/// Moment fitting reproduces polynomials only to tensor degree n-1 per
/// direction, so a mass-like integrand of tensor degree 2*p needs
/// n >= 2*p+1 per direction. exactnessOrder() must return exactly that, in
/// every direction (an off-by-one here under-integrates: measured on a p = 2
/// mass matrix, n = 2*p = 4 is not enough and n = 3 is indefinite at every
/// refinement).
TEST(ExactnessOrder)
{
    {
        const gsVector<index_t> v = gsMomentRule<real_t>::exactnessOrder(1, 2);
        CHECK_EQUAL(2, v.size());
        CHECK_EQUAL(3, v[0]);
        CHECK_EQUAL(3, v[1]);
    }
    {
        const gsVector<index_t> v = gsMomentRule<real_t>::exactnessOrder(2, 2);
        CHECK_EQUAL(2, v.size());
        CHECK_EQUAL(5, v[0]);
        CHECK_EQUAL(5, v[1]);
    }
    {
        const gsVector<index_t> v = gsMomentRule<real_t>::exactnessOrder(3, 3);
        CHECK_EQUAL(3, v.size());
        CHECK_EQUAL(7, v[0]);
        CHECK_EQUAL(7, v[1]);
        CHECK_EQUAL(7, v[2]);
    }
    // Generic form, so a special-cased implementation cannot pass by table.
    for (short_t p = 0; p != 5; ++p)
        for (short_t d = 1; d != 4; ++d)
        {
            const gsVector<index_t> v = gsMomentRule<real_t>::exactnessOrder(p, d);
            CHECK_EQUAL((index_t)d, v.size());
            for (index_t i = 0; i != v.size(); ++i)
                CHECK_EQUAL(2*(index_t)p + 1, v[i]);
        }
}

} // SUITE
