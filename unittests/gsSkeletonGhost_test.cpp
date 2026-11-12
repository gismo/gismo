/** @file gsSkeletonGhost_test.cpp

    @brief Tests for skeleton/ghost face iteration, the dnk(u,G,k) expression,
    face assembly/integration and the two-sided jump()/avg() symbols.

    Covers skeleton/ghost iteration and integration, dnk(u,G,k), face
    assembly (also into one block of a multi-block system), the two-sided
    jump()/avg() symbols and the guards that reject them outside a face loop.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst
*/

#include "gismo_unittest.h"

#include <cmath>
#include <functional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

using namespace gismo;

namespace {

/// (d-1)-measure of a face: the product of the extents in every direction
/// except the face normal, in which the box is degenerate.
real_t faceMeasure(const gsDomainIteratorWrapper<real_t> & it)
{
    const gsVector<real_t> lc = it.lowerCorner();
    const gsVector<real_t> uc = it.upperCorner();
    const index_t dir = it.side().direction();
    real_t m = 1.0;
    for (index_t j = 0; j != lc.rows(); ++j)
        if (j != dir) m *= uc[j] - lc[j];
    return m;
}

/// Coefficients of the 1D factor: sin-sequence on [pad, n-pad), zero on the
/// first and last \a pad entries. The zero padding is what makes the closed
/// form below exact, since it keeps every coefficient involved in the
/// interior, uniform-spacing region.
void separableCoefs(index_t n, index_t pad, std::vector<real_t> & a)
{
    a.assign((size_t)n, 0.0);
    for (index_t i = pad; i < n - pad; ++i)
        a[i] = math::sin(1.7*(real_t)i + 0.3);
}

/// (k+1)-th forward difference of \a a, built by k+1 passes of first
/// differences: on return d[i] = Delta^{k+1} a_i, size a.size()-(k+1).
void forwardDiff(const std::vector<real_t> & a, index_t k, std::vector<real_t> & d)
{
    d = a;
    for (index_t m = 0; m <= k; ++m)
    {
        for (size_t i = 0; i + 1 < d.size(); ++i) d[i] = d[i+1] - d[i];
        d.pop_back();
    }
}

/// Closed-form skeleton integral of [[d^k s/dn^k]]^2 for the separable spline
/// with 1D coefficients \a a on an N^d grid of spacing h: for a uniform
/// clamped degree-k B-spline s(x)=sum a_i B_i(x), s^{(k)} is piecewise
/// constant with element value Delta^k a_m/h^k, so the jump across the
/// interior knot x=j*h is Delta^{k+1} a_{j-1}/h^k. The sum below is that jump
/// squared, over the N-1 interior knot lines, times the (d-1)-measure
/// (N*h)^(d-1) of one knot line/plane.
real_t closedFormJumpSum(const std::vector<real_t> & a, index_t k,
                         real_t h, index_t N, index_t d)
{
    std::vector<real_t> delta;
    forwardDiff(a, k, delta);
    const real_t line = math::pow(N*h, d-1);
    real_t s = 0.0;
    for (index_t j = 1; j != N; ++j)
    {
        const real_t v = delta[j-1] / math::pow(h, k);
        s += v*v;
    }
    return s*line;
}

/// Redirects std::cerr into a buffer for the lifetime of the object.
/// GISMO_ERROR/GISMO_ENSURE throw std::runtime_error with a fixed message and
/// stream the actual message and condition to std::cerr, so the captured text
/// is the only way to tell which check fired.
struct CerrCapture
{
    std::ostringstream buf;
    std::streambuf * old;
    CerrCapture() : buf(), old(std::cerr.rdbuf(buf.rdbuf())) {}
    ~CerrCapture() { std::cerr.rdbuf(old); }
    std::string str() const { return buf.str(); }
};

/// Runs \a f under a fresh stderr capture. Returns true iff \a f threw a
/// std::runtime_error and the captured text contains \a needle; any other
/// exception propagates.
bool throwsWithCerr(const std::function<void()> & f, const std::string & needle)
{
    CerrCapture cap;
    bool threw = false;
    try { f(); }
    catch (const std::runtime_error &) { threw = true; }
    return threw && std::string::npos != cap.str().find(needle);
}

/// Unit square of degree \a deg, uniformly refined \a nRefine times. The
/// degree is elevated before refining, which keeps C^{deg-1} continuity
/// across the interior knots.
void unitSquare(gsMultiPatch<real_t> & mp, gsMultiBasis<real_t> & dbasis,
                index_t deg, index_t nRefine)
{
    mp.addPatch(gsNurbsCreator<real_t>::BSplineSquare(1,0,0));
    dbasis = gsMultiBasis<real_t>(mp);
    if (deg > 1) dbasis.degreeElevate(deg-1);
    dbasis.uniformRefine(nRefine);
}

/// Deterministic, non-polynomial coefficient vector. The caller keeps the
/// matrix alive: getSolution stores a pointer to it.
void sinCoefs(gsMatrix<real_t> & c, index_t n)
{
    c.resize(n, 1);
    for (index_t i = 0; i != n; ++i)
        c(i,0) = math::sin(1.7*(real_t)i + 0.3);
}

/// Runs dnk's space path (u.jump(), assembled) and solution path (sol.jump(),
/// evaluated) on the skeleton faces of \a mp, refined \a nRefine times. Both
/// loops are serial, so a thrown GISMO_ENSURE reaches the caller instead of
/// terminating inside an OpenMP region.
void runDnkMapGuardArms(const gsMultiPatch<real_t> & mp, index_t nRefine, bool expectThrow)
{
    gsMultiBasis<real_t> dbasis(mp);
    dbasis.uniformRefine(nRefine);
    gsBoundaryConditions<real_t> bc; bc.setGeoMap(mp);

    gsExprAssembler<real_t> A(1,1);
    A.setIntegrationElements(dbasis); auto G = A.getMap(mp);
    auto u = A.getSpace(dbasis); u.setup(bc, dirichlet::interpolation, 0);
    A.initSystem();
    auto dJ = dnk(u.jump(), G.left(), 1);

    // Fresh assembler: a throw in the space path may leave A's shared
    // gsExprHelper half-initialised.
    gsExprAssembler<real_t> B(1,1);
    B.setIntegrationElements(dbasis); auto Gb = B.getMap(mp);
    auto v  = B.getSpace(dbasis); v.setup(bc, dirichlet::interpolation, 0);
    B.initSystem();
    gsMatrix<real_t> coefs;
    sinCoefs(coefs, B.numDofs());
    auto sol = B.getSolution(v, coefs);
    gsExprEvaluator<real_t> ev(B);

    if (expectThrow)
    {
        CHECK(throwsWithCerr([&]{ A.assembleSkeleton( dJ*dJ.tr() ); }, "dnk_offDiagMax"));
        CHECK(throwsWithCerr([&]{ ev.integralSkeleton( dnk(sol.jump(), Gb.left(), 1).sqNorm() ); },
                             "dnk_offDiagMax"));
    }
    else
    {
        A.assembleSkeleton( dJ*dJ.tr() );
        CHECK(A.matrix().norm() > 1e-3);

        const real_t integral = ev.integralSkeleton( dnk(sol.jump(), Gb.left(), 1).sqNorm() );
        CHECK(integral > 1e-6);
    }
}

} // namespace

SUITE(gsSkeletonGhost_test)
{

// Tensor grids: closed-form counts and per-face invariants
TEST(Skeleton_UntrimmedCount_2D)
{
    gsKnotVector<real_t> kvx(0.0, 1.0, 3, 2);   // 4 elements, h = 0.25
    gsKnotVector<real_t> kvy(0.0, 1.0, 2, 2);   // 3 elements, h = 1/3
    gsTensorBSplineBasis<2,real_t> tbs(kvx, kvy);
    gsDomain<real_t>::Ptr dom = tbs.domain();   // non-owning: tbs must outlive dom

    CHECK_EQUAL((size_t)17, dom->numSkeletonFaces());   // (N-1)*M + N*(M-1) = 9 + 8

    gsDomain<real_t>::iterator endIt = dom->endSkeleton();
    size_t visited = 0;
    real_t total = 0.0;
    index_t prevDir = 0;
    size_t first = 0, second = 0;
    for (gsDomain<real_t>::iterator it = dom->beginSkeleton(); it < endIt; ++it)
    {
        // Element ids and signs live on the face iterator behind the wrapper.
        const gsDomainFaceIterator<real_t> & f =
            static_cast<const gsDomainFaceIterator<real_t>&>(*it);
        if (0 == visited) first = f.leftElementId();
        if (1 == visited) second = f.leftElementId();
        ++visited;
        total += faceMeasure(it);

        const index_t dir = it.side().direction();
        CHECK(dir >= prevDir);
        prevDir = dir;

        CHECK_CLOSE(0.0, (it.upperCorner()-it.lowerCorner())[dir], 1e-14);
        CHECK(it.side().parameter());

        const real_t hExpect = (0==dir) ? 0.25 : 1.0/3.0;
        CHECK_CLOSE(hExpect, it.getPerpendicularCellSize(), 1e-12);
        CHECK_CLOSE(hExpect, it.getPerpendicularCellSizeRight(), 1e-12);

        const size_t strideExpect = (0==dir) ? 1 : 4;
        CHECK_EQUAL(strideExpect, f.rightElementId() - f.leftElementId());

        CHECK_EQUAL(-1, f.leftSign());
        CHECK_EQUAL(-1, f.rightSign());
    }
    CHECK_EQUAL((size_t)17, visited);
    CHECK_CLOSE(5.0, total, 1e-12);          // 9*(1/3) + 8*(1/4)
    CHECK_EQUAL((size_t)0, first);
    CHECK_EQUAL((size_t)1, second);          // component 0 fastest, not the transposed 0,4

    CHECK_EQUAL((size_t)0, dom->numGhostFaces());
    gsDomain<real_t>::iterator endGhost = dom->endGhost();
    size_t nGhost = 0;
    for (gsDomain<real_t>::iterator it = dom->beginGhost(); it < endGhost; ++it) ++nGhost;
    CHECK_EQUAL((size_t)0, nGhost);
}

TEST(Skeleton_UntrimmedCount_3D)
{
    gsKnotVector<real_t> kv(0.0, 1.0, 2, 2);    // 3 elements per direction, h = 1/3
    gsTensorBSplineBasis<3,real_t> tbs(kv, kv, kv);
    gsDomain<real_t>::Ptr dom = tbs.domain();

    CHECK_EQUAL((size_t)54, dom->numSkeletonFaces());   // d*N^(d-1)*(N-1) = 3*9*2

    gsDomain<real_t>::iterator endIt = dom->endSkeleton();
    size_t visited = 0;
    real_t total = 0.0;
    index_t prevDir = 0;
    const real_t h = 1.0/3.0;
    for (gsDomain<real_t>::iterator it = dom->beginSkeleton(); it < endIt; ++it)
    {
        const gsDomainFaceIterator<real_t> & f =
            static_cast<const gsDomainFaceIterator<real_t>&>(*it);
        ++visited;
        total += faceMeasure(it);

        const index_t dir = it.side().direction();
        CHECK(dir >= prevDir);
        prevDir = dir;

        CHECK_CLOSE(0.0, (it.upperCorner()-it.lowerCorner())[dir], 1e-14);
        CHECK(it.side().parameter());
        CHECK_CLOSE(h, it.getPerpendicularCellSize(), 1e-12);
        CHECK_CLOSE(h, it.getPerpendicularCellSizeRight(), 1e-12);

        const size_t strideExpect = (0==dir) ? 1 : ((1==dir) ? 3 : 9);
        CHECK_EQUAL(strideExpect, f.rightElementId() - f.leftElementId());

        CHECK_EQUAL(-1, f.leftSign());
        CHECK_EQUAL(-1, f.rightSign());
    }
    CHECK_EQUAL((size_t)54, visited);
    CHECK_CLOSE(6.0, total, 1e-12);          // 54 * (1/3)^2

    CHECK_EQUAL((size_t)0, dom->numGhostFaces());
}

TEST(BaseDomain_NoSkeletonConcept_Throws)
{
    // Skeleton and ghost faces exist only on tensor domains of dimension >= 2;
    // the hierarchical and the 1D knot-vector domains inherit the throwing base.
    gsKnotVector<real_t> kv2(0,1,3,3);
    gsTensorBSplineBasis<2,real_t> tb(kv2,kv2);
    gsTHBSplineBasis<2,real_t> thb(tb);
    gsDomain<real_t>::Ptr hd = thb.domain();
    const gsHDomain<2,real_t> * hdp = dynamic_cast<const gsHDomain<2,real_t>*>(hd.get());
    CHECK(nullptr != hdp);

    gsKnotVector<real_t> kv(0,1,3,2);
    gsBSplineBasis<real_t> b1(kv);
    gsDomain<real_t>::Ptr kd = b1.domain();
    CHECK(nullptr != dynamic_cast<const gsKnotVector<real_t>*>(kd.get()));

    CHECK(throwsWithCerr([&]{ hd->numSkeletonFaces(); },
                         "gsDomain::numSkeletonFaces: this domain has no skeleton concept"));
    CHECK(throwsWithCerr([&]{ hd->numGhostFaces(); },
                         "gsDomain::numGhostFaces: this domain has no skeleton concept"));
    CHECK(throwsWithCerr([&]{ kd->numSkeletonFaces(); },
                         "gsDomain::numSkeletonFaces: this domain has no skeleton concept"));
    CHECK(throwsWithCerr([&]{ kd->numGhostFaces(); },
                         "gsDomain::numGhostFaces: this domain has no skeleton concept"));
}

// gsExprEvaluator: face measure and normal orientation
TEST(Evaluator_FaceMeasureAndNormal)
{
    gsMultiPatch<real_t> mp;
    gsMultiBasis<real_t> dbasis;
    unitSquare(mp, dbasis, 1, 3);                // 4 x 4 elements, J = I

    gsExprEvaluator<real_t> ev;
    ev.setIntegrationElements(dbasis);
    ev.getVariable(dbasis);                       // registers the basis degree
    auto G = ev.getMap(mp);

    gsFunctionExpr<real_t> ex("1", "0", 2);       // NAMED local (outlives ev)
    auto c = ev.getVariable(ex);

    CHECK_CLOSE(6.0, ev.integralSkeleton( nv(G.left()).norm() ), 1e-12);
    CHECK_EQUAL((size_t)24, ev.nValues());

    CHECK_CLOSE(3.0, ev.integralSkeleton( (c.tr()*nv(G.left())).val() ), 1e-12);

    CHECK_CLOSE(0.0, ev.integralGhost( nv(G.left()).norm() ), 1e-14);
}

// dnk(u,G,k) against igrad*n and n^T*ihess*n
TEST(Dnk1_EqualsIgradNormal_Assembled)
{
    gsMultiPatch<real_t> mp; gsMultiBasis<real_t> dbasis;
    unitSquare(mp, dbasis, 1, 3);
    gsBoundaryConditions<real_t> bc; bc.setGeoMap(mp);

    // matrix 1: dnk(.,.,1), all four terms with plus signs
    gsExprAssembler<real_t> A(1,1);
    A.setIntegrationElements(dbasis); auto G = A.getMap(mp);
    auto u = A.getSpace(dbasis); u.setup(bc, dirichlet::interpolation, 0);
    A.initSystem();
    CHECK_EQUAL(25, A.numDofs());                       // 5 x 5 degree-1 functions

    auto dL = dnk(u.left(),  G.left(),  1);
    auto dR = dnk(u.right(), G.right(), 1);
    A.assembleSkeleton( dL*dL.tr(), dL*dR.tr(), dR*dL.tr(), dR*dR.tr() );
    gsSparseMatrix<real_t> Mdnk = A.matrix();

    // matrix 2, fresh assembler: igrad(.,.)*nv(G.left()), same all-plus pattern
    gsExprAssembler<real_t> B(1,1);
    B.setIntegrationElements(dbasis); auto G2 = B.getMap(mp);
    auto v = B.getSpace(dbasis); v.setup(bc, dirichlet::interpolation, 0);
    B.initSystem();

    auto gL = igrad(v.left(),  G2.left())  * nv(G2.left());
    auto gR = igrad(v.right(), G2.right()) * nv(G2.left());     // nv(G2.left()) for both sides
    B.assembleSkeleton( gL*gL.tr(), gL*gR.tr(), gR*gL.tr(), gR*gR.tr() );
    gsSparseMatrix<real_t> Mig = B.matrix();

    CHECK(Mig.norm() > 1e-3);
    CHECK((Mdnk - Mig).norm() < 1e-10 * Mig.norm());
}

TEST(Dnk2_EqualsNormalHessian_OnFace)
{
    gsMultiPatch<real_t> mp; gsMultiBasis<real_t> dbasis;
    unitSquare(mp, dbasis, 2, 3);
    gsBoundaryConditions<real_t> bc; bc.setGeoMap(mp);

    gsExprAssembler<real_t> A(1,1);
    A.setIntegrationElements(dbasis); auto G = A.getMap(mp);
    auto u = A.getSpace(dbasis); u.setup(bc, dirichlet::interpolation, 0);
    A.initSystem();

    gsMatrix<real_t> coefs;
    sinCoefs(coefs, A.numDofs());
    auto sol = A.getSolution(u, coefs);

    gsExprEvaluator<real_t> ev(A);                 // shares A's gsExprHelper

    auto d1 = dnk(sol, G.left(), 1).val();
    auto d2 = dnk(sol, G.left(), 2).val();
    auto r1 = (igrad(sol,G.left())*nv(G.left())).val();
    auto r2 = (nv(G.left()).tr() * ihess(sol,G.left()) * nv(G.left())).val();

    const real_t n1 = ev.integralSkeleton( r1*r1 );
    const real_t e1 = ev.integralSkeleton( (d1-r1)*(d1-r1) );
    const real_t n2 = ev.integralSkeleton( r2*r2 );
    const real_t e2 = ev.integralSkeleton( (d2-r2)*(d2-r2) );

    CHECK_EQUAL((size_t)24, ev.nValues());
    CHECK(n1 > 1e-6);              // non-vacuity: the k=1 oracle itself is non-trivial
    CHECK(e1 < 1e-10 * n1);
    CHECK(n2 > 1e-6);              // non-vacuity: the k=2 oracle itself is non-trivial
    CHECK(e2 < 1e-10 * n2);
}

TEST(Dnk3_ThirdDerivative_3D_Dir1)
{
    // The test is 3D and singles out direction()==1 faces because a 2D k=3
    // test covers nothing: there the pure-derivative compositions (3,0) and
    // (0,3) are the first and last entries of the composition walk, i.e. both
    // endpoints of the order-3 block of evalAllDers_into, which almost any
    // indexing bug gets right by accident. Only in 3D does dir=1 reach the
    // strictly interior composition (0,3,0) of the 10 compositions of order 3.
    const index_t N = 8;
    const index_t k = 3;
    const index_t pad = k+1;
    const real_t h = 1.0/(real_t)N;

    gsMultiPatch<real_t> mp;
    mp.addPatch(gsNurbsCreator<real_t>::BSplineCube(1,0,0,0));
    gsMultiBasis<real_t> dbasis(mp);
    dbasis.degreeElevate(2);      // degree 3, elevate before refining
    dbasis.uniformRefine(N-1);    // N elements per direction
    gsBoundaryConditions<real_t> bc; bc.setGeoMap(mp);

    CHECK_EQUAL(3, dbasis.basis(0).degree(0));
    CHECK_EQUAL(3, dbasis.basis(0).degree(1));
    CHECK_EQUAL(3, dbasis.basis(0).degree(2));

    gsExprAssembler<real_t> A(1,1);
    A.setIntegrationElements(dbasis); auto G = A.getMap(mp);
    auto u = A.getSpace(dbasis); u.setup(bc, dirichlet::interpolation, 0);
    A.initSystem();

    const gsBasis<real_t> & b = dbasis.basis(0);
    CHECK_EQUAL(b.size(), A.numDofs());

    const gsTensorBSplineBasis<3,real_t> & tb =
        static_cast<const gsTensorBSplineBasis<3,real_t>&>(dbasis.basis(0));

    // y-only layout: s(x,y,z) = s_1D(y) by partition of unity in x and z, so
    // d^3_x s = d^3_z s = 0 identically and only direction()==1 faces carry
    // anything.
    std::vector<real_t> a;
    separableCoefs(tb.size(1), pad, a);
    gsMatrix<real_t> coefs(A.numDofs(), 1);
    for (index_t i = 0; i != tb.size(0); ++i)
        for (index_t j = 0; j != tb.size(1); ++j)
            for (index_t l = 0; l != tb.size(2); ++l)
                coefs(tb.index(i,j,l), 0) = a[j];
    auto sol = A.getSolution(u, coefs);

    // Random coefficients give a genuinely nonzero 3rd normal-derivative jump
    // (the space is only C^2), so the vanishing checks below cannot pass on a
    // silent zero.
    gsMatrix<real_t> coefsRandom;
    sinCoefs(coefsRandom, A.numDofs());
    auto solRandom = A.getSolution(u, coefsRandom);

    gsFunctionExpr<real_t> poly("y^3 + x*y^2 + z^3", 3);
    gsMatrix<real_t> vals;
    poly.eval_into(b.anchors(), vals);
    memory::unique_ptr< gsGeometry<real_t> > gpoly = b.interpolateAtAnchors(vals);
    gsMatrix<real_t> coefsPoly = gpoly->coefs();
    auto solPoly = A.getSolution(u, coefsPoly);

    gsExprEvaluator<real_t> ev(A);

    // Element values Delta^3 a_m/h^3 of d^3_y s on the two sides of the knot
    // plane y = j0*h.
    const index_t j0 = 5;
    const real_t expL = (a[j0+2] - 3.0*a[j0+1] + 3.0*a[j0] - a[j0-1]) / (h*h*h);  // element j0-1
    const real_t expR = (a[j0+3] - 3.0*a[j0+2] + 3.0*a[j0+1] - a[j0]) / (h*h*h);  // element j0
    CHECK(math::abs(expL) > 1e-3);
    CHECK(math::abs(expR) > 1e-3);

    // dnk itself, per face, zipped against the domain the evaluator's face
    // loop walks: the patch's subdomain is the same stored object each call.
    ev.integralSkeleton( dnk(sol.left(), G.left(), 3).val() );
    const gsDomain<real_t>::Ptr domPtr = A.domain().subdomain(0);
    const gsDomain<real_t> & dom = *domPtr;
    gsDomain<real_t>::iterator endIt = dom.endSkeleton();
    size_t i = 0, countDir1 = 0;
    for (gsDomain<real_t>::iterator it = dom.beginSkeleton(); it < endIt; ++it, ++i)
    {
        if (1 != it.side().direction())
            // Identically zero by construction (d^3_x s = d^3_z s = 0).
            CHECK(math::abs(ev.elementwise()[i]) < 1e-12*math::abs(expL)*h*h);
        else
        {
            const gsVector<real_t> lc = it.lowerCorner();
            if (math::abs(lc[1] - j0*h) < 1e-12)
            {
                CHECK_CLOSE(expL, ev.elementwise()[i]/faceMeasure(it), 1e-8*math::abs(expL));
                ++countDir1;
            }
        }
    }
    CHECK_EQUAL(ev.elementwise().size(), i);
    CHECK_EQUAL((size_t)(N*N), countDir1);

    // Closed form over the whole skeleton. The (N*h)^(d-1) factor is a knot
    // plane's area in 3D, (N*h)^2 = 1, not a knot line's length.
    const real_t ref = closedFormJumpSum(a, k, h, N, 3);
    CHECK(ref > 1e-3);
    const real_t J = ev.integralSkeleton(
          dnk(sol.jump(), G.left(), k) * dnk(sol.jump(), G.left(), k) );
    CHECK(math::abs(J - ref) < 1e-10*ref);

    // A global cubic has no 3rd normal-derivative jump. The integrand is a
    // square, so the whole-skeleton bound controls every face.
    const real_t Jrandom = ev.integralSkeleton(
          dnk(solRandom.jump(), G.left(), k) * dnk(solRandom.jump(), G.left(), k) );
    CHECK(Jrandom > 1e-3);

    const real_t Jpoly = ev.integralSkeleton(
          dnk(solPoly.jump(), G.left(), k) * dnk(solPoly.jump(), G.left(), k) );
    CHECK(math::abs(Jpoly) < 1e-10*Jrandom);
}

TEST(Dnk_RejectsNonAxisAlignedMap_2D)
{
    gsMultiPatch<real_t> mpId; mpId.addPatch(gsNurbsCreator<real_t>::BSplineSquare(1,0,0));
    runDnkMapGuardArms(mpId, 3, false);

    // Axis-aligned, non-identity: off-diagonal Jacobian entries are exactly 0.
    // Catches a guard over-tightened to J == I.
    gsMultiPatch<real_t> mpAA; mpAA.addPatch(gsNurbsCreator<real_t>::BSplineSquare(1,0,0));
    gsVector<real_t> s(2); s << 2.0, 0.5;
    gsVector<real_t> t(2); t << -1.0, 0.25;
    mpAA.patch(0).scale(s);
    mpAA.patch(0).translate(t);
    runDnkMapGuardArms(mpAA, 3, false);

    gsMultiPatch<real_t> mpRot; mpRot.addPatch(gsNurbsCreator<real_t>::BSplineSquare(1,0,0));
    mpRot.patch(0).rotate(0.3);
    runDnkMapGuardArms(mpRot, 3, true);
}

// gsExprAssembler::assembleSkeleton and the pattern pass
TEST(AssembleSkeleton_JumpForm_Deg1)
{
    gsMultiPatch<real_t> mp; gsMultiBasis<real_t> dbasis;
    unitSquare(mp, dbasis, 1, 3);
    gsBoundaryConditions<real_t> bc; bc.setGeoMap(mp);

    gsExprAssembler<real_t> A(1,1);
    A.setIntegrationElements(dbasis); auto G = A.getMap(mp);
    auto u = A.getSpace(dbasis); u.setup(bc, dirichlet::interpolation, 0);
    A.initSystem();
    CHECK_EQUAL(25, A.numDofs());

    auto dL = dnk(u.left(),  G.left(),  1);
    auto dR = dnk(u.right(), G.right(), 1);
    A.assembleSkeleton( dL*dL.tr(), -dL*dR.tr(), -dR*dL.tr(), dR*dR.tr() );
    const gsSparseMatrix<real_t> & M = A.matrix();

    // Degree-1 splines are C^0, so the normal derivative genuinely jumps.
    // This also proves that two distinct sides were evaluated: if left and
    // right resolved to the same data, all four terms would cancel.
    CHECK(M.norm() > 1e-3);
    CHECK((M - gsSparseMatrix<real_t>(M.transpose())).norm() < 1e-12);

    // Polynomial kernel: a global degree-<=1 function has a continuous normal
    // derivative. No dof is eliminated (empty BCs), so anchor order is dof order.
    gsMatrix<real_t> anch = dbasis.basis(0).anchors();      // 2 x 25
    gsMatrix<real_t> c1 = gsMatrix<real_t>::Ones(25,1);
    gsMatrix<real_t> cx = anch.row(0).transpose();
    gsMatrix<real_t> cy = anch.row(1).transpose();
    CHECK((M*c1).norm() < 1e-10);
    CHECK((M*cx).norm() < 1e-10);
    CHECK((M*cy).norm() < 1e-10);

    // Cross-element band: with 5 functions per direction and x fastest,
    // functions 0 and 2 have supports [0,0.25] and [0.25,0.75] and share no
    // element, yet the face at x=0.25 couples them. Functions 0 and 3 are
    // coupled by no face at all.
    CHECK(math::abs(M.coeff(0,2)) > 1e-6);      // cross block present
    CHECK_CLOSE(0.0, M.coeff(0,3), 1e-14);      // pattern is not simply dense
}

// gsFeSolution::right() must route to the mirrored (right) side of a face
TEST(SolutionRight_RoutesToMirror)
{
    gsMultiPatch<real_t> mp; gsMultiBasis<real_t> dbasis;
    unitSquare(mp, dbasis, 2, 3);                  // 4 x 4 elements, C^1
    gsBoundaryConditions<real_t> bc; bc.setGeoMap(mp);

    gsExprAssembler<real_t> A(1,1);
    A.setIntegrationElements(dbasis); auto G = A.getMap(mp);
    auto u = A.getSpace(dbasis); u.setup(bc, dirichlet::interpolation, 0);
    A.initSystem();

    const gsBasis<real_t> & b = dbasis.basis(0);
    CHECK_EQUAL(b.size(), A.numDofs());            // anchor order == dof order (empty BCs)

    gsMatrix<real_t> coefsRandom;
    sinCoefs(coefsRandom, A.numDofs());
    auto solRandom = A.getSolution(u, coefsRandom);

    // Interpolant of the global quadratic x^2+x*y: continuous 2nd normal
    // derivative across every face, hence a zero-jump oracle.
    gsFunctionExpr<real_t> poly("x^2 + x*y", 2);
    gsMatrix<real_t> vals;
    poly.eval_into(b.anchors(), vals);
    memory::unique_ptr< gsGeometry<real_t> > gpoly = b.interpolateAtAnchors(vals);
    gsMatrix<real_t> coefsPoly = gpoly->coefs();
    auto solPoly = A.getSolution(u, coefsPoly);

    gsExprEvaluator<real_t> ev(A);                  // shares A's gsExprHelper

    // Four-term jump [[.]]^2 = LL - LR - RL + RR, evaluated on the solution.
    // If .right() silently re-evaluated the left side, every term would
    // resolve to the same data and the sum would cancel to exactly 0.
    const real_t Jrandom = ev.integralSkeleton(
          dnk(solRandom.left(), G.left(), 2)*dnk(solRandom.left(), G.left(), 2)
        - dnk(solRandom.left(), G.left(), 2)*dnk(solRandom.right(),G.right(),2)
        - dnk(solRandom.right(),G.right(),2)*dnk(solRandom.left(), G.left(), 2)
        + dnk(solRandom.right(),G.right(),2)*dnk(solRandom.right(),G.right(),2) );

    const real_t Jpoly = ev.integralSkeleton(
          dnk(solPoly.left(), G.left(), 2)*dnk(solPoly.left(), G.left(), 2)
        - dnk(solPoly.left(), G.left(), 2)*dnk(solPoly.right(),G.right(),2)
        - dnk(solPoly.right(),G.right(),2)*dnk(solPoly.left(), G.left(), 2)
        + dnk(solPoly.right(),G.right(),2)*dnk(solPoly.right(),G.right(),2) );

    CHECK(Jrandom > 1e-6);                  // non-vacuity: the jump is genuinely nonzero
    CHECK(math::abs(Jpoly) < 1e-10 * Jrandom);   // a global quadratic has no 2nd normal-derivative jump

    // Tie the solution path to the space path: the same four-term form
    // assembled as a bilinear form on u reproduces Jrandom as c^T*M*c.
    auto dL = dnk(u.left(),  G.left(),  2);
    auto dR = dnk(u.right(), G.right(), 2);
    A.assembleSkeleton( dL*dL.tr(), -dL*dR.tr(), -dR*dL.tr(), dR*dR.tr() );
    const gsSparseMatrix<real_t> & M = A.matrix();

    const real_t quadForm = (coefsRandom.transpose()*(M*coefsRandom)).value();
    CHECK(math::abs(quadForm - Jrandom) < 1e-9 * math::abs(Jrandom));
}

// u.jump() / u.avg(): stacked two-sided symbols vs. the four-term forms
TEST(Jump_MatchesFourTerm_Deg2)
{
    gsMultiPatch<real_t> mp; gsMultiBasis<real_t> dbasis;
    unitSquare(mp, dbasis, 2, 3);                  // 4 x 4 elements, C^1
    gsBoundaryConditions<real_t> bc; bc.setGeoMap(mp);

    gsExprAssembler<real_t> A(1,1);
    A.setIntegrationElements(dbasis); auto G = A.getMap(mp);
    auto u = A.getSpace(dbasis); u.setup(bc, dirichlet::interpolation, 0);
    A.initSystem();
    CHECK_EQUAL(dbasis.basis(0).size(), A.numDofs());          // anchor order == dof order
    auto dL = dnk(u.left(),  G.left(),  2);
    auto dR = dnk(u.right(), G.right(), 2);
    A.assembleSkeleton( dL*dL.tr(), -dL*dR.tr(), -dR*dL.tr(), dR*dR.tr() );
    gsSparseMatrix<real_t> Mref = A.matrix();

    gsExprAssembler<real_t> B(1,1);
    B.setIntegrationElements(dbasis); auto G2 = B.getMap(mp);
    auto v  = B.getSpace(dbasis); v.setup(bc, dirichlet::interpolation, 0);
    B.initSystem();
    auto dJ = dnk(v.jump(), G2.left(), 2);
    B.computePatternSkeleton( dJ*dJ.tr() );
    B.assembleSkeleton( dJ*dJ.tr() );
    gsSparseMatrix<real_t> Mj = B.matrix();

    // If the derivative order did not reach the stacked entry and the second
    // derivative came back as zeros, Mj.norm() would collapse.
    CHECK(Mref.norm() > 1e-3);
    CHECK((Mj-Mref).norm() < 1e-12 * Mref.norm());
    CHECK((Mj - gsSparseMatrix<real_t>(Mj.transpose())).norm() < 1e-12 * Mj.norm());
}

TEST(Avg_MatchesFourTerm)
{
    gsMultiPatch<real_t> mp; gsMultiBasis<real_t> dbasis;
    unitSquare(mp, dbasis, 1, 3);                  // 4 x 4 elements
    gsBoundaryConditions<real_t> bc; bc.setGeoMap(mp);

    // reference: the all-plus four-term form LL+LR+RL+RR
    gsExprAssembler<real_t> A(1,1);
    A.setIntegrationElements(dbasis); auto G = A.getMap(mp);
    auto u = A.getSpace(dbasis); u.setup(bc, dirichlet::interpolation, 0);
    A.initSystem();
    auto dL = dnk(u.left(),  G.left(),  1);
    auto dR = dnk(u.right(), G.right(), 1);
    A.assembleSkeleton( dL*dL.tr(), dL*dR.tr(), dR*dL.tr(), dR*dR.tr() );
    gsSparseMatrix<real_t> Mall = A.matrix();

    gsExprAssembler<real_t> B(1,1);
    B.setIntegrationElements(dbasis); auto G2 = B.getMap(mp);
    auto v  = B.getSpace(dbasis); v.setup(bc, dirichlet::interpolation, 0);
    B.initSystem();
    auto aJ = dnk(v.avg(), G2.left(), 1);
    B.assembleSkeleton( aJ*aJ.tr() );
    gsSparseMatrix<real_t> Mavg = B.matrix();

    // avg scales both halves by 1/2, so the outer product carries 1/4 on
    // every one of the four blocks: a missing 1/2 shows up as a factor-4
    // discrepancy, not as noise.
    gsSparseMatrix<real_t> Mq = Mall;
    Mq *= 0.25;
    CHECK(Mall.norm() > 1e-3);
    CHECK((Mavg - Mq).norm() < 1e-12 * Mall.norm());
}

TEST(JumpAvg_SignConvention)
{
    // No squared-jump oracle can see a sign flip ((L-R)^2 == (R-L)^2), so the
    // convention jump = left-right, avg = (left+right)/2 needs its own test.
    // The evaluator path on a solution involves no scatter and is exact to
    // machine precision; sol.right() is verified by SolutionRight_RoutesToMirror.
    gsMultiPatch<real_t> mp; gsMultiBasis<real_t> dbasis;
    unitSquare(mp, dbasis, 2, 3);                  // 4 x 4 elements, C^1
    gsBoundaryConditions<real_t> bc; bc.setGeoMap(mp);

    gsExprAssembler<real_t> A(1,1);
    A.setIntegrationElements(dbasis); auto G = A.getMap(mp);
    auto u = A.getSpace(dbasis); u.setup(bc, dirichlet::interpolation, 0);
    A.initSystem();

    gsMatrix<real_t> coefsRandom;
    sinCoefs(coefsRandom, A.numDofs());
    auto sol = A.getSolution(u, coefsRandom);

    gsExprEvaluator<real_t> ev(A);                  // shares A's gsExprHelper

    const real_t Jsq = ev.integralSkeleton(
          dnk(sol.jump(), G.left(), 2) * dnk(sol.jump(), G.left(), 2) );
    CHECK(Jsq > 1e-6);                       // non-vacuity of the scale

    const real_t resJ = ev.integralSkeleton(
          ( dnk(sol.jump(),  G.left(),  2)
          - dnk(sol.left(),  G.left(),  2)
          + dnk(sol.right(), G.right(), 2) )
        * ( dnk(sol.jump(),  G.left(),  2)
          - dnk(sol.left(),  G.left(),  2)
          + dnk(sol.right(), G.right(), 2) ) );
    CHECK(resJ < 1e-12 * Jsq);               // jump == left - right

    const real_t Asq = ev.integralSkeleton(
          dnk(sol.avg(), G.left(), 2) * dnk(sol.avg(), G.left(), 2) );
    CHECK(Asq > 1e-6);

    // 2*avg - left - right == 0; a flipped/unscaled avg gives 4*Asq here.
    const real_t resA = ev.integralSkeleton(
          ( dnk(sol.avg(),   G.left(),  2) + dnk(sol.avg(),  G.left(),  2)
          - dnk(sol.left(),  G.left(),  2) - dnk(sol.right(), G.right(), 2) )
        * ( dnk(sol.avg(),   G.left(),  2) + dnk(sol.avg(),  G.left(),  2)
          - dnk(sol.left(),  G.left(),  2) - dnk(sol.right(), G.right(), 2) ) );
    CHECK(resA < 1e-12 * Asq);
}

TEST(ValueJump_IsOfOrderFaceShift)
{
    // With the faceShift, k-th derivative jumps are exact but jumps of VALUES
    // carry an O(faceShift) error, because the quadrature points sit inside
    // the two elements rather than on the face. The face loop shifts the
    // left/right quadrature points in the parameter direction and the map
    // here is the identity, h_L = h_R = h, so for a C^1 solution
    //   [[u]] = u_L - u_R = -(s*h_L + s*h_R)*du/dn + O(s^2)
    //   int_skel [[u]]^2 ~= 4*s^2*h^2*int_skel(du/dn)^2
    gsMultiPatch<real_t> mp; gsMultiBasis<real_t> dbasis;
    unitSquare(mp, dbasis, 2, 3);                  // 4 x 4 elements, h = 0.25
    gsBoundaryConditions<real_t> bc; bc.setGeoMap(mp);

    gsExprAssembler<real_t> A(1,1);
    A.setIntegrationElements(dbasis); auto G = A.getMap(mp);
    auto u = A.getSpace(dbasis); u.setup(bc, dirichlet::interpolation, 0);
    A.initSystem();

    gsMatrix<real_t> coefsRandom;
    sinCoefs(coefsRandom, A.numDofs());
    auto sol = A.getSolution(u, coefsRandom);

    gsExprEvaluator<real_t> ev(A);
    const real_t h = 0.25;

    // The oracle: the left-side normal derivative, squared, over the skeleton.
    const real_t n1 = ev.integralSkeleton(
          (igrad(sol,G.left())*nv(G.left())).val() * (igrad(sol,G.left())*nv(G.left())).val() );
    CHECK(n1 > 1e-6);

    // (a) value jump is tiny at the default shift.
    const real_t J0 = ev.integralSkeleton( sol.jump()*sol.jump() );
    CHECK(J0 < 1e-8);

    // (b) the average reproduces the one-sided value to O(s). This must run
    // before the faceShift is widened below: at s = 1e-3 the avg-vs-left
    // agreement degrades to O(2*s*h*u'/u) ~ 5e-4 and 1e-4 would fail.
    const real_t Aavg  = ev.integralSkeleton( sol.avg()*sol.avg()  );
    const real_t Aleft = ev.integralSkeleton( sol.left()*sol.left() );
    CHECK(Aleft > 1e-6);
    CHECK(math::abs(Aavg - Aleft) < 1e-4 * Aleft);

    // (c) the O(s) law, measured where round-off cannot reach it.
    const real_t s = 1e-3;
    ev.options().setReal("faceShift", s);
    const real_t J1 = ev.integralSkeleton( sol.jump()*sol.jump() );
    const real_t pred = 4.0*s*s*h*h*n1;
    CHECK_CLOSE(1.0, J1/pred, 0.05);          // absolute prediction, 5%
    CHECK_CLOSE(1.0e6, J1/J0, 0.05e6);        // quadratic scaling in s: (1e-3/1e-6)^2

    // (d) shifts outside (0,1), or too small/close to the far knot line to
    // move the points in floating point, are rejected.
    const std::string rangeMsg = "faceShift must lie in (0,1)";
    ev.options().setReal("faceShift", 0.0);
    CHECK(throwsWithCerr([&]{ ev.integralSkeleton( sol.jump()*sol.jump() ); }, rangeMsg));
    ev.options().setReal("faceShift", 1.0);
    CHECK(throwsWithCerr([&]{ ev.integralSkeleton( sol.jump()*sol.jump() ); }, rangeMsg));
    ev.options().setReal("faceShift", 1e-17);
    CHECK(throwsWithCerr([&]{ ev.integralSkeleton( sol.jump()*sol.jump() ); },
                         "below the floating-point resolution"));
    ev.options().setReal("faceShift", std::nextafter((real_t)1, (real_t)0));
    CHECK(throwsWithCerr([&]{ ev.integralSkeleton( sol.jump()*sol.jump() ); },
                         "opposite knot line"));
    ev.options().setReal("faceShift", 1e-6);
}

// jump()/avg() outside a two-sided loop
TEST(Jump_OneSidedLoop_Throws)
{
    // jump()/avg() are only defined on a skeleton/ghost face. Outside such a
    // loop they must throw catchably, and the guarded entry points reject
    // them before their OpenMP region (an exception cannot leave one). Each
    // check gets its own stderr capture so one entry point's message cannot
    // satisfy another's.
    gsMultiPatch<real_t> mp; gsMultiBasis<real_t> dbasis;
    unitSquare(mp, dbasis, 1, 3);
    gsBoundaryConditions<real_t> bc; bc.setGeoMap(mp);

    // assemble(): the pattern pass would otherwise run inside the guarded
    // call, so it is done up front and the guard alone decides the outcome.
    gsExprAssembler<real_t> A(1,1);
    A.setIntegrationElements(dbasis);
    auto u = A.getSpace(dbasis); u.setup(bc, dirichlet::interpolation, 0);
    A.initSystem();
    A.computePattern(u*u.tr());
    CHECK(throwsWithCerr([&]{ A.assemble(u.jump()*u.jump().tr()); },
                         "gsExprAssembler: jump()/avg()"));

    // assembleJacobian() builds its quadrature rule only when SameQuadrature
    // is off.
    gsExprAssembler<real_t> B(1,1);
    B.setIntegrationElements(dbasis);
    auto v = B.getSpace(dbasis); v.setup(bc, dirichlet::interpolation, 0);
    B.initSystem();
    gsMatrix<real_t> coefsB;
    sinCoefs(coefsB, B.numDofs());
    auto s = B.getSolution(v, coefsB);
    B.options().addSwitch("SameQuadrature", "", false);
    CHECK(throwsWithCerr([&]{ B.assembleJacobian(v*s.jump(), s); },
                         "gsExprAssembler: jump()/avg()"));

    // The six public evaluator reductions.
    gsExprAssembler<real_t> C(1,1);
    C.setIntegrationElements(dbasis);
    auto w = C.getSpace(dbasis); w.setup(bc, dirichlet::interpolation, 0);
    C.initSystem();
    gsMatrix<real_t> coefsC;
    sinCoefs(coefsC, C.numDofs());
    auto sc = C.getSolution(w, coefsC);
    gsExprEvaluator<real_t> ev(C);
    const std::string evMsg = "gsExprEvaluator: jump()/avg()";
    CHECK(throwsWithCerr([&]{ ev.integral(sc.jump()); },     evMsg));
    CHECK(throwsWithCerr([&]{ ev.integralElWise(sc.jump()); }, evMsg));
    CHECK(throwsWithCerr([&]{ ev.max(sc.jump()); },          evMsg));
    CHECK(throwsWithCerr([&]{ ev.maxElWise(sc.jump()); },    evMsg));
    CHECK(throwsWithCerr([&]{ ev.min(sc.jump()); },          evMsg));
    CHECK(throwsWithCerr([&]{ ev.minElWise(sc.jump()); },    evMsg));

    // The serial point evaluation is unguarded and stops at the helper's own
    // registration check.
    gsVector<real_t> pt(2); pt << 0.5, 0.5;
    CHECK(throwsWithCerr([&]{ ev.eval(sc.jump()*sc.jump(), pt, 0); },
                         "A jump()/avg() symbol is registered"));
}

// Multi-block scatter
TEST(MultiBlock_SkeletonScatter_HonoursBlockOffset)
{
    // Every space's gsDofMapper is shifted once, at initSystem(); the scatter
    // loop then reads the row/column indices straight off that mapper.
    // Assembling on space 1 only and comparing against a standalone
    // single-block assembly of the identical expression therefore isolates
    // the mapper shift: the face loop is serial and visits faces in a fixed
    // order, so the two local matrices are bit-identical and an exact
    // comparison is valid.
    gsMultiPatch<real_t> mp; gsMultiBasis<real_t> dbasis;
    unitSquare(mp, dbasis, 1, 3);                  // 4 x 4 elements, degree 1
    gsBoundaryConditions<real_t> bc; bc.setGeoMap(mp);

    // Standalone reference: one scalar space, block (0,0) of a 1x1 system.
    gsExprAssembler<real_t> B(1,1);
    B.setIntegrationElements(dbasis); auto Gb = B.getMap(mp);
    auto v  = B.getSpace(dbasis); v.setup(bc, dirichlet::interpolation, 0);
    B.initSystem();
    const index_t n = B.numDofs();
    CHECK_EQUAL(dbasis.basis(0).size(), n);
    auto dJb = dnk(v.jump(), Gb.left(), 1);
    B.computePatternSkeleton( dJb*dJb.tr() );
    B.assembleSkeleton( dJb*dJb.tr() );
    const gsMatrix<real_t> Dref = B.matrix().toDense();
    const real_t refMax = Dref.cwiseAbs().maxCoeff();
    CHECK(refMax > 1e-6);                                    // non-vacuity

    // Two-block system: identical expression, assembled on space 1 only.
    gsExprAssembler<real_t> A(2,2);
    A.setIntegrationElements(dbasis); auto G  = A.getMap(mp);
    auto u0 = A.getSpace(dbasis, 1, 0);
    auto u1 = A.getSpace(dbasis, 1, 1); u0.setup(bc, dirichlet::interpolation, 0);
    u1.setup(bc, dirichlet::interpolation, 0);
    A.initSystem();
    CHECK_EQUAL(2*n, A.numDofs());
    CHECK_EQUAL(n, u1.mapper().firstIndex());
    auto dJ = dnk(u1.jump(), G.left(), 1);
    A.computePatternSkeleton( dJ*dJ.tr() );
    A.assembleSkeleton( dJ*dJ.tr() );
    const gsMatrix<real_t> D2 = A.matrix().toDense();

    // cwiseAbs().maxCoeff() catches every entry, not just the aggregate
    // Frobenius scale, and ignores the structural zeros of the pattern pass.
    const real_t diag11Max = (D2.block(n,n,n,n) - Dref).cwiseAbs().maxCoeff();
    const real_t diag00Max = D2.block(0,0,n,n).cwiseAbs().maxCoeff();
    const real_t off01Max  = D2.block(0,n,n,n).cwiseAbs().maxCoeff();
    const real_t off10Max  = D2.block(n,0,n,n).cwiseAbs().maxCoeff();
    CHECK_EQUAL(0.0, diag11Max);
    CHECK_EQUAL(0.0, diag00Max);
    CHECK_EQUAL(0.0, off01Max);
    CHECK_EQUAL(0.0, off10Max);
}

// A face call leaves no two-sided state behind
TEST(FaceThenVolume_SameInstance)
{
    // A one-sided call after a face call on the same instance must not trip
    // the two-sided guards: a face call leaves no two-sided registration
    // behind.
    gsMultiPatch<real_t> mp; gsMultiBasis<real_t> dbasis;
    unitSquare(mp, dbasis, 2, 3);
    gsBoundaryConditions<real_t> bc; bc.setGeoMap(mp);

    gsExprAssembler<real_t> A(1,1);
    A.setIntegrationElements(dbasis);
    auto u = A.getSpace(dbasis); u.setup(bc, dirichlet::interpolation, 0);
    A.initSystem();
    gsMatrix<real_t> coefs;
    sinCoefs(coefs, A.numDofs());
    auto sol = A.getSolution(u, coefs);
    gsExprEvaluator<real_t> ev(A);

    const real_t Js = ev.integralSkeleton( sol.jump()*sol.jump() );
    CHECK(Js >= 0);
    real_t I = 0;
    CHECK(!throwsWithCerr([&]{ I = ev.integral(sol*sol); }, "jump()/avg()"));
    CHECK(I > 0);

    gsMultiBasis<real_t> dbasis1(mp);
    dbasis1.uniformRefine(3);
    gsExprAssembler<real_t> B(1,1);
    B.setIntegrationElements(dbasis1); auto G = B.getMap(mp);
    auto v = B.getSpace(dbasis1); v.setup(bc, dirichlet::interpolation, 0);
    B.initSystem();
    auto dJ = dnk(v.jump(), G.left(), 1);
    B.assembleSkeleton( dJ*dJ.tr() );
    B.initSystem();
    CHECK(!throwsWithCerr([&]{ B.assemble(v*v.tr()); }, "jump()/avg()"));
    CHECK(B.matrix().norm() > 0);
}

} // SUITE(gsSkeletonGhost_test)
