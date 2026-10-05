/** @file gsDirichletValues_test.cpp

    @brief Tests for strongly imposed Dirichlet values.

    On tensor bases the prescribed values produced by dirichlet::interpolation
    are compared against those produced by dirichlet::l2Projection on a
    manufactured solution of degree one. A polynomial of degree one lies in
    every spline space used here, so the trace lies in the boundary space and
    both strategies are exact: their coefficients must agree to round-off, and
    the Galerkin solution of -Laplace(u) = 0 must reproduce the exact solution.
    Any deviation is a defect in the construction of the prescribed values and
    not a discretization error.

    Tensor and rational tensor bases use dirichlet::interpolation (anchors of
    the boundary basis); that method rejects every other basis with an error.
    dirichlet::quasiInterpolation handles truncated hierarchical bases (THB,
    including rational THB; Speleers and Manni, Numer. Math. 132 (2016)
    155-184), which reproduces the THB boundary space, and non-truncated
    hierarchical bases (HB): if the Dirichlet data is the trace of a
    hierarchical spline, the prescribed values are its coefficients on the
    boundary functions. No warning is printed.
    Anchor collocation is not usable on hierarchical bases, because two
    boundary functions of different levels can share a Greville point and make
    the collocation matrix singular (see
    thb_corner_box_anchor_collocation_singular).

    A physical (non-parametric) condition needs a geometry map
    (bc.setGeoMap); a missing map is reported by an exception.

    Covers regressions of four defects:

      - anchor collocation on a hierarchical boundary basis can be rank deficient
        (a wrong Dirichlet lift without any complaint);
      - the prescribed values were read with linear (column-major) indexing,
        so every component of a vector-valued unknown received the first
        component of the boundary function;
      - a wrong Dirichlet lift on THB and rational THB patches;
      - a null geometry map dereferenced for physical conditions.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst
*/

#include "gismo_unittest.h"
#include <gsMSplines/gsMappedBasis.h>

#include <sstream>

namespace {

const index_t degree = 2;
const index_t numRef = 2;
const real_t  tol    = 1e-10;

/// Body of the warning of an l2Projection fallback, which Dirichlet setup must
/// never print. The "Warning: " prefix is not matched: colour builds put escape
/// codes between the prefix and the body.
const std::string fallbackWarning = "using dirichlet::l2Projection instead";

/// Message of the error thrown by gsDirichletValuesByTPInterpolation for a
/// basis that is not a tensor basis.
const std::string tensorOnlyMessage = "only implemented for tensor bases";

/// Remedies named by the errors of dirichlet::interpolation (hierarchical
/// bases) and dirichlet::quasiInterpolation (everything else unsupported).
const std::string qiHint = "dirichlet::quasiInterpolation";
const std::string l2Hint = "dirichlet::l2Projection";

/// Redirects std::cout (where gsWarn writes) into a buffer for the lifetime of
/// the object. The original stream buffer is restored by release() or, at the
/// latest, by the destructor, so a throw cannot leave std::cout redirected.
class CoutCapture
{
public:
    CoutCapture() : m_old(std::cout.rdbuf(m_buf.rdbuf())) { }
    ~CoutCapture() { release(); }

    /// Restores std::cout and returns everything written meanwhile
    std::string release()
    {
        if (m_old != nullptr)
        {
            std::cout.rdbuf(m_old);
            m_old = nullptr;
        }
        return m_buf.str();
    }

private:
    CoutCapture(const CoutCapture &);
    CoutCapture & operator=(const CoutCapture &);

    std::ostringstream m_buf;
    std::streambuf   * m_old;
};

/// Redirects std::cerr (where GISMO_ENSURE writes its message) into a buffer
/// for the lifetime of the object; see CoutCapture.
class CerrCapture
{
public:
    CerrCapture() : m_old(std::cerr.rdbuf(m_buf.rdbuf())) { }
    ~CerrCapture() { release(); }

    /// Restores std::cerr and returns everything written meanwhile
    std::string release()
    {
        if (m_old != nullptr)
        {
            std::cerr.rdbuf(m_old);
            m_old = nullptr;
        }
        return m_buf.str();
    }

private:
    CerrCapture(const CerrCapture &);
    CerrCapture & operator=(const CerrCapture &);

    std::ostringstream m_buf;
    std::streambuf   * m_old;
};

/// Number of non-overlapping occurrences of \a needle in \a text.
index_t countOccurrences(const std::string & text, const std::string & needle)
{
    index_t n = 0;
    for (std::string::size_type pos = text.find(needle);
         pos != std::string::npos;
         pos = text.find(needle, pos + needle.size()))
        ++n;
    return n;
}

/// Tensor B-spline basis of degree \a degree on the unit cube of dimension
/// \a d, uniformly refined \a numRef times.
template <short_t d>
typename gsTensorBSplineBasis<d,real_t>::uPtr unitCubeBasis()
{
    gsKnotVector<real_t> kv(0.0, 1.0, 0, degree+1);
    std::vector< gsBSplineBasis<real_t>* > cb(d);
    for (short_t i=0; i!=d; ++i) cb[i] = new gsBSplineBasis<real_t>(kv);
    typename gsTensorBSplineBasis<d,real_t>::uPtr tb(
        new gsTensorBSplineBasis<d,real_t>(cb) );
    for (index_t r=0; r!=numRef; ++r) tb->uniformRefine();
    return tb;
}

/// Discretization basis over \a tb. The hierarchical variant is refined with
/// two nested boxes anchored at the origin, so that the refinement cuts the
/// west/south (/front) sides and the boundary bases there are of mixed level.
template <short_t d>
typename gsBasis<real_t>::uPtr discretizationBasis(
    const gsTensorBSplineBasis<d,real_t> & tb, bool useTHB)
{
    if (!useTHB) return tb.clone();

    gsTHBSplineBasis<d,real_t> * thb = new gsTHBSplineBasis<d,real_t>(tb);
    gsMatrix<real_t> boxes(d, 4);
    boxes.col(0).setConstant(0.00); boxes.col(1).setConstant(0.50);
    boxes.col(2).setConstant(0.00); boxes.col(3).setConstant(0.25);
    thb->refine(boxes);
    return typename gsBasis<real_t>::uPtr(thb);
}

/// Number of points a tensor grid of component anchors would carry on \a side.
template <short_t d>
index_t tensorAnchorGridSize(const gsBasis<real_t> & b, boxSide side)
{
    index_t n = 1;
    for (short_t i=0; i!=d; ++i)
        if (i != side.direction()) n *= b.component(i).size();
    return n;
}

/// Exact solution of degree one, one expression per target component.
gsFunctionExpr<real_t> exactSolution(short_t d, index_t nComp)
{
    const std::string e0 = (2==d) ? "1 + 2*x + 3*y" : "1 + 2*x + 3*y + 4*z";
    const std::string e1 = (2==d) ? "4 + 5*x + 6*y" : "5 + 6*x + 7*y + 8*z";
    if (1==nComp) return gsFunctionExpr<real_t>(e0, d);
    if (2==nComp) return gsFunctionExpr<real_t>(e0, e1, d);
    return gsFunctionExpr<real_t>(e0, e1, "9 - x - y - z", d);
}

gsBoundaryConditions<real_t> allSidesDirichlet(const gsMultiPatch<real_t> & mp,
                                               gsFunctionExpr<real_t> & g,
                                               short_t d)
{
    gsBoundaryConditions<real_t> bc;
    for (index_t s=1; s<=2*d; ++s)
        bc.addCondition(0, boxSide(s), condition_type::dirichlet, &g);
    bc.setGeoMap(mp);
    return bc;
}

/// Largest relative difference between the prescribed values obtained with
/// \a strategy and with dirichlet::l2Projection, for an unknown
/// with \a nComp components.
template <short_t d>
real_t prescribedValueDifference(const gsMultiPatch<real_t> & mp,
                                 const gsMultiBasis<real_t> & dbasis,
                                 const dirichlet::values strategy,
                                 index_t nComp)
{
    gsFunctionExpr<real_t> g = exactSolution(d, nComp);
    gsBoundaryConditions<real_t> bc = allSidesDirichlet(mp, g, d);

    gsExprAssembler<real_t> A(1,1);
    A.setIntegrationDomain(dbasis.domain());
    gsExprAssembler<real_t>::geometryMap G = A.getMap(mp);
    GISMO_UNUSED(G);
    gsExprAssembler<real_t>::space u = A.getSpace(dbasis, nComp);

    u.setup(bc, strategy, 0);
    const gsMatrix<real_t> interp = u.fixedPart();
    u.setup(bc, dirichlet::l2Projection, 0);
    const gsMatrix<real_t> proj = u.fixedPart();

    // A vacuous comparison of two empty vectors must not be read as a pass
    CHECK(interp.size() > 0);
    CHECK_EQUAL(interp.size(), proj.size());

    return (interp - proj).cwiseAbs().maxCoeff()
         / math::max( (real_t)1, proj.cwiseAbs().maxCoeff() );
}

/// L2 error of the Galerkin solution of -Laplace(u) = f against the exact
/// solution \a g, with the Dirichlet values built by \a strategy. If
/// \a setupOutput is given, it receives what the Dirichlet setup (and nothing
/// else) wrote to std::cout.
template <short_t d>
real_t poissonError(const gsMultiPatch<real_t> & mp,
                    const gsMultiBasis<real_t> & dbasis,
                    dirichlet::values strategy,
                    gsFunctionExpr<real_t> & g,
                    const gsFunctionExpr<real_t> & fSource,
                    std::string * setupOutput = nullptr)
{
    gsBoundaryConditions<real_t> bc = allSidesDirichlet(mp, g, d);

    gsExprAssembler<real_t> A(1,1);
    A.setIntegrationDomain(dbasis.domain());
    gsExprEvaluator<real_t> ev(A);

    gsExprAssembler<real_t>::geometryMap G = A.getMap(mp);
    gsExprAssembler<real_t>::space u = A.getSpace(dbasis);
    auto f = A.getCoeff(fSource, G);
    auto g_ex = ev.getVariable(g, G);

    gsMatrix<real_t> solVector;
    gsExprAssembler<real_t>::solution u_sol = A.getSolution(u, solVector);

    {
        CoutCapture capture;
        u.setup(bc, strategy, 0);
        const std::string out = capture.release();
        if (setupOutput != nullptr) *setupOutput = out;
    }
    A.initSystem();
    A.assemble( igrad(u,G) * igrad(u,G).tr() * meas(G), u * f * meas(G) );

    gsSparseSolver<real_t>::SimplicialLDLT solver;
    solver.compute( A.matrix() );
    solVector = solver.solve( A.rhs() );

    return math::sqrt( ev.integral( (u_sol-g_ex).sqNorm() * meas(G) ) );
}

/// L2 error of the Galerkin solution of -Laplace(u) = 0 against the exact
/// solution of degree one.
template <short_t d>
real_t poissonError(const gsMultiPatch<real_t> & mp,
                    const gsMultiBasis<real_t> & dbasis,
                    dirichlet::values strategy,
                    std::string * setupOutput = nullptr)
{
    gsFunctionExpr<real_t> g = exactSolution(d, 1);
    const gsFunctionExpr<real_t> zero("0", d);
    return poissonError<d>(mp, dbasis, strategy, g, zero, setupOutput);
}

/// Unit square with the identity geometry, a p=3 tensor B-spline basis with
/// 4x4 elements, and its THB counterpart refined once in the corner box
/// [0,0.5]^2.
///
/// On the west side the boundary basis then has 9 functions; the level-0
/// function B_2 and the level-1 function B_3 have the same Greville anchor
/// 0.25, so collocation at the anchors has rank 8 of 9. The east side is
/// untouched by the refinement (7 functions, one level).
struct CornerBoxCase
{
    CornerBoxCase()
        : kv(0.0, 1.0, 3, 4), tb(kv, kv), thb(tb)
    {
        mp.addPatch( gsTensorBSpline<2,real_t>(tb, tb.anchors().transpose()) );

        gsMatrix<real_t> box(2,2);
        box << 0, 0.5,
               0, 0.5;
        thb.refine(box);

        mb.reset(new gsMultiBasis<real_t>(thb));
        mb->setTopology(mp);
    }

    gsKnotVector<real_t>               kv;
    gsTensorBSplineBasis<2,real_t>     tb;
    gsTHBSplineBasis<2,real_t>         thb;
    gsMultiPatch<real_t>               mp;
    std::unique_ptr< gsMultiBasis<real_t> > mb;

private:
    CornerBoxCase(const CornerBoxCase &);
    CornerBoxCase & operator=(const CornerBoxCase &);
};

/// Manufactured solution and source on the corner-box case. The cubic trace
/// is contained in the level-0 space and hence in the THB space, so a correct
/// Dirichlet lift reproduces it to solver precision. Its second derivative
/// depends on y, which keeps the lift on the west/east sides non-trivial.
gsFunctionExpr<real_t> cornerBoxSolution() { return gsFunctionExpr<real_t>("1 + y - 2*y^2 + y^3", 2); }
gsFunctionExpr<real_t> cornerBoxSource()   { return gsFunctionExpr<real_t>("4 - 6*y", 2); }

/// Deterministic values sin(a*i + b), i = 0..n-1, in [-1,1]. A fixed sequence
/// keeps filtered and full-suite runs identical (gsMatrix::Random shares the
/// std::rand state of the whole binary).
gsMatrix<real_t> spreadValues(index_t n, real_t a, real_t b)
{
    gsMatrix<real_t> v(n, 1);
    for (index_t i=0; i!=n; ++i)
        v(i,0) = math::sin(a * (real_t)i + b);
    return v;
}

/// Prescribes the trace of \a spline (a function over \a basis's parameter
/// domain) as Dirichlet data on every side of a single patch and returns the
/// largest |fixed dof - coefs(i)| over all i in basis.boundary(side), all
/// sides. The geometry map \a mp is always set: the l2 projection evaluates
/// its measure even for parametric data. \a out receives what the setup wrote
/// to std::cout, \a nCompared the number of compared values.
template <short_t d>
real_t traceCoefError(const gsBasis<real_t> & basis, const gsMultiPatch<real_t> & mp,
                      const gsFunction<real_t> & spline, const gsMatrix<real_t> & coefs,
                      dirichlet::values strategy, bool parametric,
                      std::string & out, index_t & nCompared)
{
    gsMultiBasis<real_t> mb(basis);
    gsExprAssembler<real_t> A(1,1);
    A.setIntegrationElements(mb);
    gsExprAssembler<real_t>::space u = A.getSpace(mb, 1, 0);
    gsBoundaryConditions<real_t> bc;
    bc.setGeoMap(mp);
    for (boxSide s = boxSide::getFirst(d); s < boxSide::getEnd(d); ++s)
        bc.addCondition(0, s, condition_type::dirichlet, spline, 0, parametric, -1);
    { CoutCapture capture; u.setup(bc, strategy, 0); out = capture.release(); }

    real_t err = 0;
    nCompared = 0;
    for (boxSide s = boxSide::getFirst(d); s < boxSide::getEnd(d); ++s)
    {
        const gsMatrix<index_t> idx = basis.boundary(s);
        for (index_t l = 0; l != idx.rows(); ++l)
        {
            const real_t fixed = u.fixedPart().at( u.mapper().bindex(idx(l,0), 0, 0) );
            err = math::max(err, math::abs(fixed - coefs(idx(l,0), 0)));
            ++nCompared;
        }
    }
    return err;
}

/// Runs every check for one dimension and one basis type.
template <short_t d>
void runCase(bool useTHB)
{
    typename gsTensorBSplineBasis<d,real_t>::uPtr tb = unitCubeBasis<d>();

    gsMultiPatch<real_t> mp;
    mp.addPatch( gsTensorBSpline<d,real_t>(*tb, tb->anchors().transpose()) );

    typename gsBasis<real_t>::uPtr dbas = discretizationBasis<d>(*tb, useTHB);
    gsMultiBasis<real_t> dbasis(*dbas);
    dbasis.setTopology(mp);

    // The hierarchical case is only a regression test if the refinement
    // actually reaches a Dirichlet side: otherwise every boundary basis is a
    // plain tensor basis and the defect cannot be reproduced.
    if (useTHB)
    {
        bool anyHierarchical = false;
        for (index_t s=1; s<=2*d; ++s)
            if ( dbasis.basis(0).boundaryBasis(boxSide(s))->size()
                 != tensorAnchorGridSize<d>(dbasis.basis(0), boxSide(s)) )
                anyHierarchical = true;
        CHECK(anyHierarchical);
    }

    // The boundary basis must number its functions like
    // basis.boundary(side), which is what the prescribed values are
    // scattered with.
    for (index_t s=1; s<=2*d; ++s)
        CHECK_EQUAL( dbasis.basis(0).boundary(boxSide(s)).size(),
                     dbasis.basis(0).boundaryBasis(boxSide(s))->size() );

    const dirichlet::values method = useTHB ? dirichlet::quasiInterpolation : dirichlet::interpolation;

    // Scalar unknown
    CHECK_CLOSE( 0.0, prescribedValueDifference<d>(mp, dbasis, method, 1), tol );

    std::string interpOut;
    CHECK_CLOSE( 0.0, poissonError<d>(mp, dbasis, method, &interpOut), tol );
    CHECK_CLOSE( 0.0, poissonError<d>(mp, dbasis, dirichlet::l2Projection ), tol );
    CHECK_EQUAL( 0, countOccurrences(interpOut, fallbackWarning) );

    // Vector-valued unknown: component r must be taken from component r
    // of the boundary function
    CHECK_CLOSE( 0.0, prescribedValueDifference<d>(mp, dbasis, method, d), tol );
}

} // anonymous namespace

SUITE(gsDirichletValues_test)
{

TEST(tensor_2d)  { runCase<2>(false); }
TEST(hierarchical_2d) { runCase<2>(true ); }
TEST(tensor_3d)  { runCase<3>(false); }
TEST(hierarchical_3d) { runCase<3>(true ); }

// Why anchor collocation is unusable on THB boundaries (and quasi-
// interpolation is used instead): on the west side of the corner-box THB basis
// two functions of different levels share the Greville anchor 0.25, so the
// square collocation matrix is singular. The east side has a single level and
// full rank.
TEST(thb_corner_box_anchor_collocation_singular)
{
    CornerBoxCase c;

    gsBasis<real_t>::uPtr west = c.thb.boundaryBasis(boundary::west);
    CHECK_EQUAL( 9, west->size() );
    {
        const gsMatrix<real_t> anchors = west->anchors();
        index_t atQuarter = 0;
        for (index_t i=0; i!=anchors.cols(); ++i)
            if ( math::abs(anchors(0,i) - 0.25) < 1e-14 ) ++atQuarter;
        CHECK_EQUAL( 2, atQuarter );

        const gsMatrix<real_t> C = west->collocationMatrix(anchors).toDense();
        gsEigen::FullPivLU< gsMatrix<real_t> > lu(C);
        CHECK_EQUAL( 8, lu.rank() );
    }

    gsBasis<real_t>::uPtr east = c.thb.boundaryBasis(boundary::east);
    CHECK_EQUAL( 7, east->size() );
    {
        const gsMatrix<real_t> anchors = east->anchors();
        const gsMatrix<real_t> C = east->collocationMatrix(anchors).toDense();
        gsEigen::FullPivLU< gsMatrix<real_t> > lu(C);
        CHECK_EQUAL( 7, lu.rank() );
    }
}

// The exact solution lies in the THB space, so a correct Dirichlet lift
// reproduces it to solver precision. THB quasi-interpolation prints no
// warning, and neither does l2 projection.
TEST(thb_corner_box_quasiInterpolation_exact)
{
    CornerBoxCase c;
    gsFunctionExpr<real_t> g = cornerBoxSolution();
    const gsFunctionExpr<real_t> f = cornerBoxSource();

    std::string interpOut, projOut;
    const real_t errInterp = poissonError<2>(c.mp, *c.mb, dirichlet::quasiInterpolation, g, f, &interpOut);
    const real_t errProj   = poissonError<2>(c.mp, *c.mb, dirichlet::l2Projection , g, f, &projOut  );

    CHECK_CLOSE( 0.0, errInterp, tol );
    CHECK_EQUAL( 0, countOccurrences(interpOut, fallbackWarning) );

    CHECK_CLOSE( 0.0, errProj, tol );
    CHECK_EQUAL( 0, countOccurrences(projOut, fallbackWarning) );
}

// Same problem through the legacy assembler, whose quasi-interpolation is the
// one of gsDirichletValues.
TEST(legacy_assembler_thb_quasiInterpolation_exact)
{
    CornerBoxCase c;
    gsFunctionExpr<real_t> g = cornerBoxSolution();
    const gsFunctionExpr<real_t> f = cornerBoxSource();

    gsBoundaryConditions<real_t> bc = allSidesDirichlet(c.mp, g, 2);

    // Identity geometry: parametric and physical coordinates coincide
    const index_t n = 11;
    gsMatrix<real_t> pts(2, n*n);
    for (index_t i=0; i!=n; ++i)
        for (index_t j=0; j!=n; ++j)
        {
            pts(0, i*n+j) = (real_t)i / (n-1);
            pts(1, i*n+j) = (real_t)j / (n-1);
        }

    real_t maxErr = -1;
    std::string failure, out;
    try
    {
        CoutCapture capture;
        gsPoissonAssembler<real_t> PA(c.mp, *c.mb, bc, f);
        PA.options().setInt("DirichletValues", dirichlet::quasiInterpolation);
        PA.assemble();

        gsSparseSolver<real_t>::SimplicialLDLT solver;
        solver.compute( PA.matrix() );
        const gsMatrix<real_t> solVector = solver.solve( PA.rhs() );

        gsMultiPatch<real_t> sol;
        PA.constructSolution(solVector, sol);
        out = capture.release();

        const gsMatrix<real_t> vals = sol.patch(0).eval(pts);
        const gsMatrix<real_t> ex   = g.eval(pts);
        maxErr = (vals - ex).cwiseAbs().maxCoeff();
    }
    catch (const std::exception & e) { failure = e.what(); }
    catch (...)                      { failure = "non-standard exception"; }

    CHECK_EQUAL( std::string(), failure );
    CHECK_CLOSE( 0.0, maxErr, tol );
    CHECK_EQUAL( 0, countOccurrences(out, fallbackWarning) );
}

// A Dirichlet function whose target dimension cannot supply the requested
// component is rejected rather than silently read as its first component.
TEST(component_mismatch_throws)
{
    typename gsTensorBSplineBasis<2,real_t>::uPtr tb = unitCubeBasis<2>();
    gsMultiPatch<real_t> mp;
    mp.addPatch( gsTensorBSpline<2,real_t>(*tb, tb->anchors().transpose()) );
    gsMultiBasis<real_t> dbasis(*tb);
    dbasis.setTopology(mp);

    // Scalar function prescribed for a two-component unknown, no component
    // selected on the condition
    gsFunctionExpr<real_t> g("1 + 2*x + 3*y", 2);
    gsBoundaryConditions<real_t> bc = allSidesDirichlet(mp, g, 2);

    gsExprAssembler<real_t> A(1,1);
    A.setIntegrationDomain(dbasis.domain());
    gsExprAssembler<real_t>::geometryMap G = A.getMap(mp);
    GISMO_UNUSED(G);
    gsExprAssembler<real_t>::space u = A.getSpace(dbasis, 2);

    CHECK_THROW( u.setup(bc, dirichlet::interpolation, 0), std::runtime_error );
}

// Tensor and rational tensor bases use anchor interpolation: the boundary
// indices, face anchors and interpolateAtAnchors all resolve through
// gsRationalBasis::source(). THB and rational THB bases use
// dirichlet::quasiInterpolation without a warning. A mapped basis is computed
// with dirichlet::l2Projection, whereas gsDirichletValuesByTPInterpolation
// itself rejects every non-tensor basis.
// Every case must reproduce the exact trace, either directly or through an
// equivalent tensor reference.
TEST(dirichletInterpolationBasisGuard)
{
    // Case 1: genuine NURBS (non-unit weights). A rational basis is a
    // partition of unity, so interpolating a constant must reproduce it
    // exactly in every boundary coefficient.
    {
        gsMultiPatch<> mp( *gsNurbsCreator<>::NurbsQuarterAnnulus(1,2) );
        gsMultiBasis<> mb(mp);
        mb.uniformRefine();

        CHECK(mb.basis(0).isRational());
        const gsTensorNurbsBasis<2,real_t> * nurbsBasis =
            dynamic_cast<const gsTensorNurbsBasis<2,real_t>*>(&mb.basis(0));
        CHECK(nullptr != nurbsBasis);

        if (nurbsBasis != nullptr)
        {
            CHECK( (nurbsBasis->weights().array() - 1).abs().maxCoeff() > 1e-3 );

            // The curved side is not known a priori from the construction
            // alone; find it by scanning the boundary weights instead of
            // assuming which side is curved.
            boundary::side chosenSide = boundary::none;
            real_t chosenDev = 0;
            for (boundary::side s : {boundary::west, boundary::south})
            {
                const gsMatrix<index_t> idx = mb.basis(0).boundary(s);
                real_t dev = 0;
                for (index_t i = 0; i != idx.rows(); ++i)
                    dev = std::max(dev, std::abs(nurbsBasis->weights()(idx(i,0),0) - 1));
                gsInfo << "[dirichletNurbs] side "<<boxSide(s)<<" weight deviation = "<<dev<<"\n";
                if (dev > 1e-3 && chosenSide == boundary::none)
                    chosenSide = s;
                if (dev > 1e-3 && chosenDev < dev) chosenDev = dev;
            }
            gsInfo << "[dirichletNurbs] chosen side = "<<boxSide(chosenSide)
                   << ", deviation = "<<chosenDev<<"\n";
            CHECK(chosenSide != boundary::none);

            if (chosenSide != boundary::none)
            {
                gsExprAssembler<> A(1,1);
                A.setIntegrationElements(mb);
                gsExprAssembler<>::space u = A.getSpace(mb, 1, 0);
                gsBoundaryConditions<> bc;
                bc.setGeoMap(mp);
                gsFunctionExpr<> f("3", 2);
                bc.addCondition(0, chosenSide, condition_type::dirichlet, &f, 0, false, -1);

                bool threw = false;
                std::string out;
                {
                    CoutCapture capture;
                    try { u.setup(bc, dirichlet::interpolation, 0); }
                    catch (const std::runtime_error &) { threw = true; }
                    out = capture.release();
                }
                gsInfo << "[dirichletNurbs] case1 threw = "<<threw<<"\n";
                CHECK(!threw);
                CHECK_EQUAL( 0, countOccurrences(out, fallbackWarning) );

                if (!threw)
                {
                    CHECK(u.mapper().boundarySize() > 0);
                    const gsMatrix<index_t> idx = mb.basis(0).boundary(chosenSide);
                    for (index_t i = 0; i != idx.rows(); ++i)
                        CHECK_CLOSE(3.0, u.fixedPart().at( u.mapper().bindex(idx(i,0), 0, 0) ), 1e-8);
                }
            }
        }
    }

    // Fixture shared by cases 2-5: a plain tensor B-spline basis on the
    // unit square. The parametric trace of g on the north side, 3+x^2,
    // is not constant (so a permuted boundary numbering changes the
    // prescribed values) and lies in every boundary space used here (so
    // interpolation and L2 projection must agree).
    gsKnotVector<> kv(0, 1, 3, 3);
    gsTensorBSplineBasis<2> tb(kv, kv);

    gsMultiPatch<> mp;
    mp.addPatch( *gsNurbsCreator<>::BSplineSquare(1.0, 0.0, 0.0) );
    mp.computeTopology();

    gsFunctionExpr<> g("1+x*x+2*y", 2);

    // Prescribed values of g on the north side of mbX with the given
    // strategy; threw reports whether setup raised, out receives what setup
    // wrote to std::cout and full the coefficients of the trace with respect
    // to all functions of mbX (zero away from the north side).
    auto northValues = [&](const gsMultiBasis<> & mbX, const index_t method, bool & threw,
                           std::string & out, gsMatrix<> & full)
    {
        gsExprAssembler<> AX(1,1);
        AX.setIntegrationElements(mbX);
        gsExprAssembler<>::space uX = AX.getSpace(mbX, 1, 0);
        gsBoundaryConditions<> bcX;
        bcX.setGeoMap(mp);
        bcX.addCondition(0, boundary::north, condition_type::dirichlet, &g, 0, true, -1);
        threw = false;
        {
            CoutCapture capture;
            try { uX.setup(bcX, method, 0); }
            catch (const std::runtime_error &) { threw = true; }
            out = capture.release();
        }
        full.setZero(mbX.basis(0).size(), 1);
        if (!threw)
        {
            const gsMatrix<index_t> idx = mbX.basis(0).boundary(boundary::north);
            for (index_t i = 0; i != idx.rows(); ++i)
                full(idx(i,0), 0) = uX.fixedPart().at( uX.mapper().bindex(idx(i,0), 0, 0) );
        }
        return gsMatrix<>(uX.fixedPart());
    };

    // Largest deviation of the trace with coefficients full from the exact
    // north trace 3 + x^2, sampled along the north side.
    auto northTraceError = [&](const gsBasis<> & b, const gsMatrix<> & full)
    {
        gsMatrix<> pts(2, 21);
        for (index_t i = 0; i != pts.cols(); ++i)
        {
            pts(0,i) = (real_t)i / (pts.cols()-1);
            pts(1,i) = 1.0;
        }
        const gsMatrix<> vals = b.makeGeometry(full)->eval(pts);
        real_t err = 0;
        for (index_t i = 0; i != pts.cols(); ++i)
            err = math::max(err, math::abs(vals(0,i) - (3.0 + pts(0,i)*pts(0,i))));
        return err;
    };

    auto sameValues = [](const gsMatrix<> & a, const gsMatrix<> & b)
    {
        return a.rows() == b.rows() && a.cols() == b.cols() &&
               (a - b).cwiseAbs().maxCoeff() < 1e-8;
    };

    bool threwB = false;
    std::string outB;
    gsMatrix<> fullB;
    const gsMatrix<> valsB = northValues(gsMultiBasis<>(tb), dirichlet::interpolation, threwB, outB, fullB);
    CHECK(!threwB);
    CHECK_EQUAL( 0, countOccurrences(outB, fallbackWarning) );
    CHECK( valsB.size() > 0 && valsB.maxCoeff() - valsB.minCoeff() > 0.1 );

    // Case 2: unit-weight NURBS spans the tensor B-spline space, and
    // interpolation is unique, so the rational path (generic collocation
    // + BiCGSTABILUT) must match the tensor path to solver tolerance.
    {
        gsTensorNurbsBasis<2> nb(tb);
        bool threwN = false;
        std::string outN;
        gsMatrix<> fullN;
        const gsMatrix<> valsN = northValues(gsMultiBasis<>(nb), dirichlet::interpolation, threwN, outN, fullN);
        gsInfo << "[dirichletGuard] case2 threwN = "<<threwN<<"\n";
        CHECK(!threwN);
        CHECK_EQUAL( 0, countOccurrences(outN, fallbackWarning) );
        if (!threwB && !threwN) CHECK( sameValues(valsB, valsN) );
    }

    // Case 3: a single-level THB basis spans the tensor space and is not a
    // tensor basis: the values come from quasi-interpolation (no warning) and
    // must equal the tensor interpolation.
    gsTHBSplineBasis<2> thb(tb);
    {
        bool threwH = false;
        std::string outH;
        gsMatrix<> fullH;
        const gsMatrix<> valsH = northValues(gsMultiBasis<>(thb), dirichlet::quasiInterpolation, threwH, outH, fullH);
        gsInfo << "[dirichletGuard] case3 threwH = "<<threwH<<"\n";
        CHECK(!threwH);
        CHECK_EQUAL( 0, countOccurrences(outH, fallbackWarning) );
        if (!threwB && !threwH) CHECK( sameValues(valsB, valsH) );
    }

    // Case 4: THB refined in [0,0.5]x[0.5,1], so the north boundary basis
    // is of mixed level. THB and unit-weight rational THB are both
    // quasi-interpolated without a warning; the trace 3 + x^2 lies in the
    // boundary space, so the prescribed values must reproduce it.
    thb.refine( (gsMatrix<>(2,2) << 0, 0.5, 0.5, 1).finished() );
    CHECK( thb.boundaryBasis(boundary::north)->size() > tb.boundaryBasis(boundary::north)->size() );
    {
        bool threwI = false, threwR = false;
        std::string outI, outR;
        gsMatrix<> fullI, fullR;
        const gsMatrix<> valsI = northValues(gsMultiBasis<>(thb), dirichlet::quasiInterpolation, threwI, outI, fullI);
        gsRationalTHBSplineBasis<2> rthb(thb);
        const gsMatrix<> valsR = northValues(gsMultiBasis<>(rthb), dirichlet::quasiInterpolation, threwR, outR, fullR);
        gsInfo << "[dirichletGuard] case4 threwI = "<<threwI<<", threwR = "<<threwR<<"\n";
        CHECK(!threwI);
        CHECK(!threwR);
        CHECK_EQUAL( 0, countOccurrences(outI, fallbackWarning) );
        CHECK_EQUAL( 0, countOccurrences(outR, fallbackWarning) );
        CHECK( valsI.size() > 0 && valsR.size() > 0 );
        if (!threwI) CHECK_CLOSE( 0.0, northTraceError(thb , fullI), 1e-8 );
        if (!threwR) CHECK_CLOSE( 0.0, northTraceError(rthb, fullR), 1e-8 );
    }

    // Case 5: a mapped basis is not a tensor basis. Through
    // gsDirichletValues with dirichlet::l2Projection and, since the identity
    // mapping spans the tensor space, it reproduces the tensor values;
    // dirichlet::interpolation itself rejects it.
    {
        gsMultiBasis<> mbM(tb);
        const index_t sz = mbM.basis(0).size();
        gsSparseMatrix<> ident(sz, sz);
        ident.setIdentity();
        gsMappedBasis<2, real_t> mapB(mbM, ident);

        gsExprAssembler<> AM(1,1);
        AM.setIntegrationElements(mbM);
        gsExprAssembler<>::space uM = AM.getSpace(mapB, 1, 0);
        gsBoundaryConditions<> bcM;
        bcM.setGeoMap(mp);
        bcM.addCondition(0, boundary::north, condition_type::dirichlet, &g, 0, true, -1);

        // (i) l2 projection through gsDirichletValues (dirichlet::interpolation
        // on a mapped basis throws)
        bool threwM = false;
        std::string outM;
        {
            CoutCapture capture;
            try { uM.setup(bcM, dirichlet::l2Projection, 0); }
            catch (const std::runtime_error &) { threwM = true; }
            outM = capture.release();
        }
        gsInfo << "[dirichletGuard] case5 threwM = "<<threwM<<"\n";
        CHECK(!threwM);
        CHECK_EQUAL( 0, countOccurrences(outM, fallbackWarning) );
        if (!threwM && !threwB) CHECK( sameValues(valsB, gsMatrix<>(uM.fixedPart())) );

        // (ii) the interpolation routine rejects it. GISMO_ENSURE throws the
        // literal "GISMO_ENSURE" and writes its message to std::cerr, so the
        // message is captured there.
        uM.setup(bcM, dirichlet::l2Projection, 0);
        std::ostringstream err;
        std::streambuf * cerrBuf = std::cerr.rdbuf(err.rdbuf());
        bool threwTP = false;
        try { gsDirichletValuesByTPInterpolation(uM, bcM); }
        catch (const std::runtime_error &) { threwTP = true; }
        std::cerr.rdbuf(cerrBuf);
        gsInfo << "[dirichletGuard] case5 threwTP = "<<threwTP<<", stderr = "<<err.str()<<"\n";
        CHECK(threwTP);
        CHECK( err.str().find(tensorOnlyMessage) != std::string::npos );
    }

    // Case 6: the direct tensor interpolation routine also rejects a refined
    // THB basis.
    {
        gsMultiBasis<> mbT(thb);
        gsExprAssembler<> AT(1,1);
        AT.setIntegrationElements(mbT);
        gsExprAssembler<>::space uT = AT.getSpace(mbT, 1, 0);
        gsBoundaryConditions<> bcT;
        bcT.setGeoMap(mp);
        bcT.addCondition(0, boundary::north, condition_type::dirichlet, &g, 0, true, -1);

        uT.setup(bcT, dirichlet::l2Projection, 0);
        std::ostringstream err;
        std::streambuf * cerrBuf = std::cerr.rdbuf(err.rdbuf());
        bool threw = false;
        try { gsDirichletValuesByTPInterpolation(uT, bcT); }
        catch (const std::runtime_error &) { threw = true; }
        std::cerr.rdbuf(cerrBuf);
        CHECK(threw);
        CHECK( err.str().find(tensorOnlyMessage) != std::string::npos );
    }
}


namespace {

/// Mixed-level discriminator: the refinement must enlarge the boundary basis
/// of \a side beyond that of the tensor basis \a tb.
bool boundaryGrown(const gsBasis<real_t> & basis, const gsBasis<real_t> & tb, boxSide side)
{
    return basis.boundaryBasis(side)->size() > tb.boundaryBasis(side)->size();
}

} // anonymous namespace

// The trace of a THB spline with non-polynomial (sinusoidal) coefficients lies
// in the boundary space but not in any polynomial space, so only an
// interpolant that reproduces the THB boundary space returns its coefficients.
TEST(dirichletTrace_thb_parametric)
{
    CornerBoxCase c;
    const gsMatrix<real_t> pc = spreadValues(c.thb.size(), 0.7, 0.3);
    gsGeometry<real_t>::uPtr spline = c.thb.makeGeometry(pc);

    CHECK( boundaryGrown(c.thb, c.tb, boundary::west) );

    std::string out; index_t n = 0;
    const real_t err = traceCoefError<2>(c.thb, c.mp, *spline, pc,
                                         dirichlet::quasiInterpolation, true, out, n);
    CHECK( n > 0 );
    CHECK_EQUAL( 0, countOccurrences(out, fallbackWarning) );
    CHECK( err <= 1e-12 );
}

// Same data prescribed as a physical condition; the identity geometry makes
// the parametric and physical traces coincide, so the lift must be identical.
TEST(dirichletTrace_thb_physical)
{
    CornerBoxCase c;
    const gsMatrix<real_t> pc = spreadValues(c.thb.size(), 0.7, 0.3);
    gsGeometry<real_t>::uPtr spline = c.thb.makeGeometry(pc);

    CHECK( boundaryGrown(c.thb, c.tb, boundary::west) );

    std::string out; index_t n = 0;
    const real_t err = traceCoefError<2>(c.thb, c.mp, *spline, pc,
                                         dirichlet::quasiInterpolation, false, out, n);
    CHECK( n > 0 );
    CHECK_EQUAL( 0, countOccurrences(out, fallbackWarning) );
    CHECK( err <= 1e-12 );
}

// Refinement one level-1 element thick along the west side: the west trace may
// consist of level-1 functions only, whereas the south side is of mixed level.
TEST(dirichletTrace_thb_thinStrip)
{
    CornerBoxCase c;
    gsTHBSplineBasis<2,real_t> thb(c.tb);
    gsMatrix<real_t> box(2,2);
    box << 0, 0.125,
           0, 1;
    thb.refine(box);
    const gsMatrix<real_t> pc = spreadValues(thb.size(), 0.7, 0.3);
    gsGeometry<real_t>::uPtr spline = thb.makeGeometry(pc);

    CHECK( boundaryGrown(thb, c.tb, boundary::west) );
    CHECK( boundaryGrown(thb, c.tb, boundary::south) );

    std::string out; index_t n = 0;
    const real_t err = traceCoefError<2>(thb, c.mp, *spline, pc,
                                         dirichlet::quasiInterpolation, true, out, n);
    CHECK( n > 0 );
    CHECK_EQUAL( 0, countOccurrences(out, fallbackWarning) );
    CHECK( err <= 1e-12 );
}

// Three-dimensional counterpart with a corner box refined once.
TEST(dirichletTrace_thb_3d)
{
    gsKnotVector<real_t> kv(0.0, 1.0, 3, 4);
    std::vector< gsKnotVector<real_t> > kvs(3, kv);
    gsTensorBSplineBasis<3,real_t> tb(kvs);
    gsMultiPatch<real_t> mp;
    mp.addPatch( gsTensorBSpline<3,real_t>(tb, tb.anchors().transpose()) );

    gsTHBSplineBasis<3,real_t> thb(tb);
    gsMatrix<real_t> box(3,2);
    box << 0, 0.5,
           0, 0.5,
           0, 0.5;
    thb.refine(box);
    const gsMatrix<real_t> pc = spreadValues(thb.size(), 0.7, 0.3);
    gsGeometry<real_t>::uPtr spline = thb.makeGeometry(pc);

    CHECK( boundaryGrown(thb, tb, boundary::west) );

    std::string out; index_t n = 0;
    const real_t err = traceCoefError<3>(thb, mp, *spline, pc,
                                         dirichlet::quasiInterpolation, true, out, n);
    CHECK( n > 0 );
    CHECK_EQUAL( 0, countOccurrences(out, fallbackWarning) );
    CHECK( err <= 1e-12 );
}

// Rational THB with weights in [0.5,2]: the value of the spline is
// sum_i c_i w_i T_i / W, so the prescribed values of its trace are the c_i.
// Unit weights would not distinguish a rational-aware lift from a polynomial
// one.
TEST(dirichletTrace_rationalTHB_genuineWeights)
{
    CornerBoxCase c;
    gsRationalTHBSplineBasis<2,real_t> rthb(c.thb);
    const gsMatrix<real_t> W = 1.25 + 0.75 * spreadValues(c.thb.size(), 1.3, 0.1).array();
    rthb.setWeights(W);
    CHECK( (W.array() - 1).abs().maxCoeff() > 0.1 );

    const gsMatrix<real_t> pc = spreadValues(c.thb.size(), 0.7, 0.3);
    gsGeometry<real_t>::uPtr spline = rthb.makeGeometry(pc);

    CHECK( boundaryGrown(rthb, c.tb, boundary::west) );

    std::string out; index_t n = 0;
    const real_t err = traceCoefError<2>(rthb, c.mp, *spline, pc,
                                         dirichlet::quasiInterpolation, true, out, n);
    CHECK( n > 0 );
    CHECK_EQUAL( 0, countOccurrences(out, fallbackWarning) );
    CHECK( err <= 1e-11 );
}

// The HB trace data are reproduced by the level-by-level quasi-interpolation,
// without a warning.
TEST(dirichletTrace_hb_quasiInterpolation)
{
    CornerBoxCase c;
    gsHBSplineBasis<2,real_t> hb(c.tb);
    gsMatrix<real_t> box(2,2);
    box << 0, 0.5,
           0, 0.5;
    hb.refine(box);
    const gsMatrix<real_t> pc = spreadValues(hb.size(), 0.7, 0.3);
    gsGeometry<real_t>::uPtr spline = hb.makeGeometry(pc);

    CHECK( boundaryGrown(hb, c.tb, boundary::west) );

    std::string out; index_t n = 0;
    const real_t err = traceCoefError<2>(hb, c.mp, *spline, pc,
                                         dirichlet::quasiInterpolation, true, out, n);
    CHECK( n > 0 );
    CHECK_EQUAL( 0, countOccurrences(out, fallbackWarning) );
    CHECK( err <= 1e-8 );
}

/// Prescribed values on \a side alone, in the order of basis.boundary(side).
/// \a out receives what the setup wrote to std::cout.
gsMatrix<real_t> prescribedOnSide(const gsBasis<real_t> & basis, const gsMultiPatch<real_t> & mp,
                                  const gsFunction<real_t> & data, bool parametric,
                                  boxSide side, dirichlet::values strategy, std::string & out)
{
    gsMultiBasis<real_t> mb(basis);
    gsExprAssembler<real_t> A(1,1);
    A.setIntegrationElements(mb);
    gsExprAssembler<real_t>::space u = A.getSpace(mb, 1, 0);
    gsBoundaryConditions<real_t> bc;
    bc.setGeoMap(mp);
    bc.addCondition(0, side, condition_type::dirichlet, data, 0, parametric, -1);
    { CoutCapture capture; u.setup(bc, strategy, 0); out = capture.release(); }

    const gsMatrix<index_t> idx = basis.boundary(side);
    gsMatrix<real_t> vals(idx.rows(), 1);
    for (index_t l = 0; l != idx.rows(); ++l)
        vals(l,0) = u.fixedPart().at( u.mapper().bindex(idx(l,0), 0, 0) );
    return vals;
}

/// Largest entry-wise difference of two equally sized arrays.
real_t maxDiff(const gsMatrix<real_t> & a, const gsMatrix<real_t> & b)
{
    CHECK_EQUAL( a.rows(), b.rows() );
    CHECK_EQUAL( a.cols(), b.cols() );
    return (a - b).cwiseAbs().maxCoeff();
}

// Physical data on a curved patch. The patch maps the west side onto a curve
// x = 0.1 sin(pi y), so a lift that evaluated g at the parametric points
// instead of the physical ones would prescribe different values.
TEST(dirichletTrace_thb_curvedPatch_physical)
{
    CornerBoxCase c;
    gsMatrix<real_t> cp = c.tb.anchors().transpose();
    for (index_t i = 0; i != cp.rows(); ++i)
    {
        const real_t x = cp(i,0), y = cp(i,1);
        cp(i,0) = x + 0.1 * math::sin(EIGEN_PI * y) * (1 + x);
        cp(i,1) = y + 0.05 * (1 + x) * math::sin(EIGEN_PI * y);
    }
    gsMultiPatch<real_t> curved;
    curved.addPatch( gsTensorBSpline<2,real_t>(c.tb, cp) );

    gsFunctionExpr<real_t> g("sin(2*x) + y^2", 2);
    CHECK( boundaryGrown(c.thb, c.tb, boundary::west) );

    const gsGeometrySlice<real_t> westCurve(&curved.patch(0), 0, 0.0);
    const gsComposedFunction<real_t> physicalTrace(&westCurve, &g);
    gsMatrix<real_t> oracle;
    gsQuasiInterpolate<real_t>::localIntpl(*c.thb.boundaryBasis(boundary::west), physicalTrace, oracle);

    // The premise: the physical oracle is not the one of g at parametric points
    const gsGeometrySlice<real_t> westParam(&g, 0, 0.0);
    gsMatrix<real_t> parametricOracle;
    gsQuasiInterpolate<real_t>::localIntpl(*c.thb.boundaryBasis(boundary::west), westParam, parametricOracle);
    CHECK( maxDiff(oracle, parametricOracle) > 1e-2 );

    std::string out;
    const gsMatrix<real_t> fixed = prescribedOnSide(c.thb, curved, g, false, boundary::west,
                                                    dirichlet::quasiInterpolation, out);
    CHECK( fixed.size() > 0 );
    CHECK_EQUAL( 0, countOccurrences(out, fallbackWarning) );
    CHECK( maxDiff(fixed, oracle) <= 1e-12 );
}

// Data outside the spline space: the prescribed values are those of hierarchical
// quasi-interpolation, which differ from the l2 projection.
TEST(dirichletTrace_thb_outOfSpace_matchesQI)
{
    CornerBoxCase c;
    gsFunctionExpr<real_t> g("exp(y)*cos(3*y)", 2);
    CHECK( boundaryGrown(c.thb, c.tb, boundary::west) );

    const gsGeometrySlice<real_t> westTrace(&g, 0, 0.0);
    gsMatrix<real_t> oracle;
    gsQuasiInterpolate<real_t>::localIntpl(*c.thb.boundaryBasis(boundary::west), westTrace, oracle);

    std::string out, outL2;
    const gsMatrix<real_t> fixed = prescribedOnSide(c.thb, c.mp, g, true, boundary::west,
                                                    dirichlet::quasiInterpolation, out);
    const gsMatrix<real_t> l2 = prescribedOnSide(c.thb, c.mp, g, true, boundary::west,
                                                 dirichlet::l2Projection, outL2);
    CHECK( fixed.size() > 0 );
    CHECK( maxDiff(fixed, oracle) <= 1e-12 );
    CHECK( maxDiff(l2, oracle) > 1e-6 );
    CHECK_EQUAL( 0, countOccurrences(out, fallbackWarning) );
}

/// THB basis of degree three on the parameter domain [0,4]^2 (unit knot
/// spans), refined in [2,4]^2 and once more in [3,4]^2, so that the east and
/// north boundary bases are of mixed level whereas the west and south ones
/// are not.
struct NonUnitDomainCase
{
    NonUnitDomainCase()
        : kv(0.0, 4.0, 3, 4), tb(kv, kv), thb(tb)
    {
        mp.addPatch( gsTensorBSpline<2,real_t>(tb, tb.anchors().transpose()) );
        gsMatrix<real_t> box(2,2);
        box << 2, 4,
               2, 4;
        thb.refine(box);
        box << 3, 4,
               3, 4;
        thb.refine(box);
        mb.reset(new gsMultiBasis<real_t>(thb));
        mb->setTopology(mp);
    }

    gsKnotVector<real_t>               kv;
    gsTensorBSplineBasis<2,real_t>     tb;
    gsTHBSplineBasis<2,real_t>         thb;
    gsMultiPatch<real_t>               mp;
    std::unique_ptr< gsMultiBasis<real_t> > mb;

private:
    NonUnitDomainCase(const NonUnitDomainCase &);
    NonUnitDomainCase & operator=(const NonUnitDomainCase &);
};

// Boundary coefficients of a THB spline on a non-unit parameter domain, through
// the parametric lift with both strategies.
TEST(dirichletTrace_thb_nonUnitDomain)
{
    NonUnitDomainCase c;
    const gsMatrix<real_t> pc = spreadValues(c.thb.size(), 0.7, 0.3);
    gsGeometry<real_t>::uPtr spline = c.thb.makeGeometry(pc);

    CHECK( boundaryGrown(c.thb, c.tb, boundary::east) );
    CHECK( boundaryGrown(c.thb, c.tb, boundary::north) );

    std::string out; index_t n = 0;
    const real_t err = traceCoefError<2>(c.thb, c.mp, *spline, pc,
                                         dirichlet::quasiInterpolation, true, out, n);
    CHECK( n > 0 );
    CHECK_EQUAL( 0, countOccurrences(out, fallbackWarning) );
    CHECK( err <= 1e-12 );

    std::string outL2; index_t nL2 = 0;
    const real_t errL2 = traceCoefError<2>(c.thb, c.mp, *spline, pc,
                                           dirichlet::l2Projection, true, outL2, nL2);
    CHECK( nL2 > 0 );
    CHECK( errL2 <= 1e-8 );
}

// The same problem through the legacy assembler: a solution constructed from
// the assembled system carries the Dirichlet values as coefficients of the
// boundary functions.
TEST(legacy_assembler_thb_nonUnitDomain_quasiInterpolation)
{
    NonUnitDomainCase c;
    const gsMatrix<real_t> pc = spreadValues(c.thb.size(), 0.7, 0.3);
    gsGeometry<real_t>::uPtr spline = c.thb.makeGeometry(pc);
    const gsFunctionExpr<real_t> f("0", 2);

    gsBoundaryConditions<real_t> bc;
    bc.setGeoMap(c.mp);
    for (boxSide s = boxSide::getFirst(2); s < boxSide::getEnd(2); ++s)
        bc.addCondition(0, s, condition_type::dirichlet, *spline, 0, true, -1);

    real_t err = -1;
    index_t nCompared = 0;
    std::string failure, out;
    try
    {
        CoutCapture capture;
        gsPoissonAssembler<real_t> PA(c.mp, *c.mb, bc, f);
        PA.options().setInt("DirichletValues", dirichlet::quasiInterpolation);
        PA.assemble();

        gsSparseSolver<real_t>::SimplicialLDLT solver;
        solver.compute( PA.matrix() );
        const gsMatrix<real_t> solVector = solver.solve( PA.rhs() );

        gsMultiPatch<real_t> sol;
        PA.constructSolution(solVector, sol);
        out = capture.release();

        err = 0;
        for (boxSide s = boxSide::getFirst(2); s < boxSide::getEnd(2); ++s)
        {
            const gsMatrix<index_t> idx = c.thb.boundary(s);
            for (index_t l = 0; l != idx.rows(); ++l, ++nCompared)
                err = math::max(err, math::abs(sol.patch(0).coefs()(idx(l,0), 0) - pc(idx(l,0), 0)));
        }
    }
    catch (const std::exception & e) { failure = e.what(); }
    catch (...)                      { failure = "non-standard exception"; }

    CHECK_EQUAL( std::string(), failure );
    CHECK( nCompared > 0 );
    CHECK( err >= 0 && err <= 1e-12 );
    CHECK_EQUAL( 0, countOccurrences(out, fallbackWarning) );
}

namespace {

/// Curved single patch over the tensor basis \a tb: the identity map bent by
/// smooth sinusoidal offsets, so that the physical image of a boundary point
/// differs from its parameter by O(0.1).
gsMultiPatch<real_t> curvedPatch(const gsTensorBSplineBasis<2,real_t> & tb)
{
    gsMatrix<real_t> cp = tb.anchors().transpose();
    for (index_t i = 0; i != cp.rows(); ++i)
    {
        const real_t x = cp(i,0), y = cp(i,1);
        cp(i,0) = x + 0.1 * math::sin(EIGEN_PI * y) * (1 + x);
        cp(i,1) = y + 0.05 * (1 + x) * math::sin(EIGEN_PI * y);
    }
    gsMultiPatch<real_t> curved;
    curved.addPatch( gsTensorBSpline<2,real_t>(tb, cp) );
    return curved;
}

/// Dirichlet dofs of the legacy gsPoissonAssembler (source zero) with
/// \a data prescribed on all four sides of the single patch \a mp, read back
/// as the coefficients of the constructed solution on the boundary functions
/// of \a basis. \a vals receives the full coefficient vector of the solution,
/// \a out what the assembly wrote to std::cout and \a failure the text of an
/// exception, if one was thrown.
void legacyDirichletCoefs(const gsBasis<real_t> & basis, const gsMultiPatch<real_t> & mp,
                          const gsFunction<real_t> & data, bool parametric,
                          dirichlet::values strategy,
                          gsMatrix<real_t> & vals, std::string & out, std::string & failure)
{
    gsMultiBasis<real_t> mb(basis);
    mb.setTopology(mp);
    gsBoundaryConditions<real_t> bc;
    bc.setGeoMap(mp);
    for (boxSide s = boxSide::getFirst(2); s < boxSide::getEnd(2); ++s)
        bc.addCondition(0, s, condition_type::dirichlet, data, 0, parametric, -1);
    const gsFunctionExpr<real_t> f("0", 2);

    try
    {
        CoutCapture capture;
        gsPoissonAssembler<real_t> PA(mp, mb, bc, f);
        PA.options().setInt("DirichletValues", strategy);
        PA.assemble();

        gsSparseSolver<real_t>::SimplicialLDLT solver;
        solver.compute( PA.matrix() );
        const gsMatrix<real_t> solVector = solver.solve( PA.rhs() );

        gsMultiPatch<real_t> sol;
        PA.constructSolution(solVector, sol);
        out = capture.release();
        vals = sol.patch(0).coefs();
    }
    catch (const std::exception & e) { failure = e.what(); }
    catch (...)                      { failure = "non-standard exception"; }
}

/// Largest |vals(i) - coefs(i)| over the boundary functions i of \a basis on
/// all sides; \a nCompared receives the number of compared entries.
real_t boundaryCoefError(const gsBasis<real_t> & basis, const gsMatrix<real_t> & vals,
                         const gsMatrix<real_t> & coefs, index_t & nCompared)
{
    real_t err = 0;
    nCompared = 0;
    for (boxSide s = boxSide::getFirst(2); s < boxSide::getEnd(2); ++s)
    {
        const gsMatrix<index_t> idx = basis.boundary(s);
        for (index_t l = 0; l != idx.rows(); ++l, ++nCompared)
            err = math::max(err, math::abs(vals(idx(l,0), 0) - coefs(idx(l,0), 0)));
    }
    return err;
}

/// Identity patch over the unit square on the p=3 tensor basis \a tb.
gsMultiPatch<real_t> identityPatch(const gsTensorBSplineBasis<2,real_t> & tb)
{
    gsMultiPatch<real_t> mp;
    mp.addPatch( gsTensorBSpline<2,real_t>(tb, tb.anchors().transpose()) );
    return mp;
}

} // anonymous namespace

// Parametric data on a curved patch through the legacy l2 projection: the
// data are a function of the parameter, so the prescribed values of the trace
// of a spline over the boundary space are its coefficients. Interpolation,
// which handles parametric data, serves as the control.
TEST(legacy_l2_parametric_tensor_curvedPatch)
{
    gsKnotVector<real_t> kv(0.0, 1.0, 3, 4);
    gsTensorBSplineBasis<2,real_t> tb(kv, kv);
    const gsMultiPatch<real_t> curved = curvedPatch(tb);
    const gsMatrix<real_t> pc = spreadValues(tb.size(), 0.7, 0.3);
    const gsTensorBSpline<2,real_t> spline(tb, pc);

    // Premise: the physical image of the boundary is far from the parameter
    gsMatrix<real_t> probe(2,1);
    probe << 0.0, 0.5;
    CHECK( (curved.patch(0).eval(probe) - probe).norm() > 1e-2 );

    for (int k = 0; k != 2; ++k)
    {
        const dirichlet::values strategy = (0 == k) ? dirichlet::interpolation
                                                    : dirichlet::l2Projection;
        gsMatrix<real_t> vals;
        std::string out, failure;
        legacyDirichletCoefs(tb, curved, spline, true, strategy, vals, out, failure);

        index_t n = 0;
        CHECK_EQUAL( std::string(), failure );
        CHECK( vals.size() == pc.size() );
        const real_t err = (vals.size() == pc.size()) ? boundaryCoefError(tb, vals, pc, n) : -1;
        CHECK( n > 0 );
        CHECK( err >= 0 && err <= tol );
    }
}

// Part of the east side is refined. A quadrature point on the side, mapped
// through the identity geometry, must still evaluate the parametric data at
// the parameter itself: rounding of the mapped point to u = 1 + eps would
// leave the THB support, where the spline vanishes.
TEST(legacy_l2_parametric_thb_partlyRefinedSide)
{
    gsKnotVector<real_t> kv(0.0, 1.0, 3, 4);
    gsTensorBSplineBasis<2,real_t> tb(kv, kv);
    const gsMultiPatch<real_t> mp = identityPatch(tb);

    gsTHBSplineBasis<2,real_t> thb(tb);
    gsMatrix<real_t> box(2,2);
    box << 0.75, 1.0,
           0.5,  1.0;
    thb.refine(box);
    CHECK( boundaryGrown(thb, tb, boundary::east) );

    const gsMatrix<real_t> pc = spreadValues(thb.size(), 0.7, 0.3);
    gsGeometry<real_t>::uPtr spline = thb.makeGeometry(pc);

    gsMatrix<real_t> vals;
    std::string out, failure;
    legacyDirichletCoefs(thb, mp, *spline, true, dirichlet::l2Projection, vals, out, failure);

    index_t n = 0;
    CHECK_EQUAL( std::string(), failure );
    CHECK( vals.size() == pc.size() );
    const real_t err = (vals.size() == pc.size()) ? boundaryCoefError(thb, vals, pc, n) : -1;
    CHECK( n > 0 );
    CHECK( err >= 0 && err <= tol );
}

// Non-truncated HB basis: the legacy quasi-interpolation reproduces the HB
// trace, without a warning.
TEST(legacy_qi_parametric_hb)
{
    gsKnotVector<real_t> kv(0.0, 1.0, 3, 4);
    gsTensorBSplineBasis<2,real_t> tb(kv, kv);
    const gsMultiPatch<real_t> mp = identityPatch(tb);

    gsHBSplineBasis<2,real_t> hb(tb);
    gsMatrix<real_t> box(2,2);
    box << 0.75, 1.0,
           0.5,  1.0;
    hb.refine(box);
    CHECK( boundaryGrown(hb, tb, boundary::east) );

    const gsMatrix<real_t> pc = spreadValues(hb.size(), 0.7, 0.3);
    gsGeometry<real_t>::uPtr spline = hb.makeGeometry(pc);

    gsMatrix<real_t> vals;
    std::string out, failure;
    legacyDirichletCoefs(hb, mp, *spline, true, dirichlet::quasiInterpolation, vals, out, failure);

    index_t n = 0;
    CHECK_EQUAL( std::string(), failure );
    CHECK_EQUAL( 0, countOccurrences(out, fallbackWarning) );
    CHECK( vals.size() == pc.size() );
    const real_t err = (vals.size() == pc.size()) ? boundaryCoefError(hb, vals, pc, n) : -1;
    CHECK( n > 0 );
    CHECK( err >= 0 && err <= tol );
}

// Physical data on a curved patch through the legacy l2 projection agree with
// the expression-path projection of the same data. g = 1 + 2x - 3y composed
// with the geometry lies in the patch's spline space, so the projection is
// also the coefficient vector 1 + 2 cx - 3 cy of the geometry.
TEST(legacy_l2_physical_curvedPatch_control)
{
    gsKnotVector<real_t> kv(0.0, 1.0, 3, 4);
    gsTensorBSplineBasis<2,real_t> tb(kv, kv);
    const gsMultiPatch<real_t> curved = curvedPatch(tb);
    const gsFunctionExpr<real_t> g("1 + 2*x - 3*y", 2);

    gsMatrix<real_t> vals;
    std::string out, failure;
    legacyDirichletCoefs(tb, curved, g, false, dirichlet::l2Projection, vals, out, failure);
    CHECK_EQUAL( std::string(), failure );
    CHECK( vals.rows() == tb.size() );

    gsMultiBasis<real_t> mb(tb);
    mb.setTopology(curved);
    gsBoundaryConditions<real_t> bc;
    bc.setGeoMap(curved);
    for (boxSide s = boxSide::getFirst(2); s < boxSide::getEnd(2); ++s)
        bc.addCondition(0, s, condition_type::dirichlet, g, 0, false, -1);
    gsExprAssembler<real_t> A(1,1);
    A.setIntegrationElements(mb);
    gsExprAssembler<real_t>::geometryMap G = A.getMap(curved);
    GISMO_UNUSED(G);
    gsExprAssembler<real_t>::space u = A.getSpace(mb, 1, 0);
    u.setup(bc, dirichlet::l2Projection, 0);

    const gsMatrix<real_t> & cg = curved.patch(0).coefs();
    real_t errExpr = 0, errGeo = 0;
    index_t n = 0;
    for (boxSide s = boxSide::getFirst(2); s < boxSide::getEnd(2); ++s)
    {
        const gsMatrix<index_t> idx = tb.boundary(s);
        for (index_t l = 0; l != idx.rows() && vals.rows() == tb.size(); ++l, ++n)
        {
            const index_t i = idx(l,0);
            const real_t fixed = u.fixedPart().at( u.mapper().bindex(i, 0, 0) );
            errExpr = math::max(errExpr, math::abs(vals(i,0) - fixed));
            errGeo  = math::max(errGeo,  math::abs(vals(i,0) - (1 + 2*cg(i,0) - 3*cg(i,1))));
        }
    }
    CHECK( n > 0 );
    CHECK( errExpr <= tol );
    CHECK( errGeo <= tol );
}

namespace {

/// Dirichlet dofs of a scalar unknown over a multibasis: fixedPart() and, per
/// patch, the value of the fixed dof of every local function (zero where the
/// function is free).
struct Lift
{
    gsMatrix<real_t>                 fixed;
    std::vector< gsMatrix<real_t> >  local;
    std::vector< std::vector<char> > isFixed;

    index_t numFixed(index_t k) const
    {
        index_t n = 0;
        for (size_t i = 0; i != isFixed[k].size(); ++i) n += isFixed[k][i];
        return n;
    }
};

/// Sets up the unknown over \a mb with \a bc and \a method. \a out receives
/// what the setup wrote to std::cout.
Lift liftDirichlet(const gsMultiBasis<real_t> & mb, const gsBoundaryConditions<real_t> & bc,
                   index_t method, std::string * out = nullptr)
{
    gsExprAssembler<real_t> A(1,1);
    A.setIntegrationElements(mb);
    gsExprAssembler<real_t>::space u = A.getSpace(mb, 1, 0);
    {
        CoutCapture capture;
        u.setup(bc, method, 0);
        const std::string o = capture.release();
        if (out != nullptr) *out = o;
    }

    Lift L;
    L.fixed = u.fixedPart();
    L.local.resize(mb.nBases());
    L.isFixed.resize(mb.nBases());
    for (size_t k = 0; k != mb.nBases(); ++k)
    {
        const index_t n = mb.basis(k).size();
        L.local[k].setZero(n, 1);
        L.isFixed[k].assign(n, 0);
        for (index_t i = 0; i != n; ++i)
            if (u.mapper().is_boundary(i, k, 0))
            {
                L.isFixed[k][i] = 1;
                L.local[k](i,0) = u.fixedPart().at( u.mapper().bindex(i, k, 0) );
            }
    }
    return L;
}

/// Largest entry-wise difference of two arrays; infinity if their shapes
/// differ, so that a `== 0` check on it cannot pass on mismatched sizes.
real_t exactDifference(const gsMatrix<real_t> & a, const gsMatrix<real_t> & b)
{
    if (a.rows() != b.rows() || a.cols() != b.cols())
        return std::numeric_limits<real_t>::infinity();
    return (a - b).cwiseAbs().maxCoeff();
}

/// Dirichlet condition \a data (whose piece k is used) on every side of patch \a k.
void addAllSides(gsBoundaryConditions<real_t> & bc, index_t k, const gsFunctionSet<real_t> & data,
                 bool parametric, short_t d = 2)
{
    for (boxSide s = boxSide::getFirst(d); s < boxSide::getEnd(d); ++s)
        bc.addCondition(k, s, condition_type::dirichlet, data, 0, parametric, -1);
}

/// Cubic-in-each-direction HB counterpart of CornerBoxCase::thb: tensor basis
/// \a tb refined in [0,0.5]^2.
gsHBSplineBasis<2,real_t> cornerBoxHB(const gsTensorBSplineBasis<2,real_t> & tb)
{
    gsHBSplineBasis<2,real_t> hb(tb);
    gsMatrix<real_t> box(2,2);
    box << 0, 0.5,
           0, 0.5;
    hb.refine(box);
    return hb;
}

/// Two unit squares side by side without a shared interface: patch 0 is the
/// identity map over the tensor basis \a tb (p=3, 4x4 elements), patch 1 the
/// same map translated by (2,0) and carrying the THB basis \a thb of \a tb
/// refined in [0,0.5]^2.
struct MixedTensorTHBCase
{
    MixedTensorTHBCase()
        : kv(0.0, 1.0, 3, 4), tb(kv, kv), thb(tb)
    {
        mp.addPatch( gsTensorBSpline<2,real_t>(tb, tb.anchors().transpose()) );
        gsMatrix<real_t> cp = tb.anchors().transpose();
        cp.col(0).array() += 2.0;
        mp.addPatch( gsTensorBSpline<2,real_t>(tb, cp) );
        mp.computeTopology();

        gsMatrix<real_t> box(2,2);
        box << 0, 0.5,
               0, 0.5;
        thb.refine(box);

        mb.addBasis(tb.clone().release());
        mb.addBasis(thb.clone().release());
        mb.setTopology(mp);
    }

    gsKnotVector<real_t>           kv;
    gsTensorBSplineBasis<2,real_t> tb;
    gsTHBSplineBasis<2,real_t>     thb;
    gsMultiPatch<real_t>           mp;
    gsMultiBasis<real_t>           mb;

private:
    MixedTensorTHBCase(const MixedTensorTHBCase &);
    MixedTensorTHBCase & operator=(const MixedTensorTHBCase &);
};

/// Quarter annulus (NURBS with non-unit weights), refined once.
struct AnnulusCase
{
    AnnulusCase()
        : mp( *gsNurbsCreator<real_t>::NurbsQuarterAnnulus(1,2) ), mb(mp)
    {
        mb.uniformRefine();
    }

    gsMultiPatch<real_t> mp;
    gsMultiBasis<real_t> mb;

private:
    AnnulusCase(const AnnulusCase &);
    AnnulusCase & operator=(const AnnulusCase &);
};

/// End values computed by the Dirichlet method \a method with the 1-D basis \a b on the unit interval:
/// data cos(x)+x^2 on both ends, either parametric or in physical coordinates of
/// the geometry with basis \a b and coefficients 2 + 3 anchors (about
/// xi -> 2 + 3 xi). A clamped basis is 1 at the end point of its end function,
/// so the geometry takes the coefficient of that function there. Returns the
/// largest deviation of the two end dofs from the data value at the end point;
/// \a separation receives |g(parametric west end) - g(physical west end)|.
/// \a method is the Dirichlet method under test.
real_t oneDimEndError(const gsBasis<real_t> & b, bool parametric, real_t & separation,
                      index_t method = dirichlet::quasiInterpolation)
{
    gsMatrix<real_t> coefs = (2 + 3 * b.anchors().array()).matrix().transpose();
    gsMultiPatch<real_t> mp;
    mp.addPatch( b.makeGeometry(coefs) );

    gsFunctionExpr<real_t> g("cos(x)+x^2", 1);
    gsBoundaryConditions<real_t> bc;
    bc.setGeoMap(mp);
    bc.addCondition(0, boundary::west, condition_type::dirichlet, &g, 0, parametric, -1);
    bc.addCondition(0, boundary::east, condition_type::dirichlet, &g, 0, parametric, -1);

    gsMultiBasis<real_t> mb(b);
    const Lift L = liftDirichlet(mb, bc, method);

    const index_t w = b.boundary(boundary::west)(0,0), e = b.boundary(boundary::east)(0,0);
    const real_t cw = coefs(w,0), ce = coefs(e,0);
    const real_t gWestParam = math::cos(0.0) + 0.0, gEastParam = math::cos(1.0) + 1.0;
    const real_t gWestPhys  = math::cos(cw) + cw*cw, gEastPhys  = math::cos(ce) + ce*ce;
    separation = math::abs(gWestParam - gWestPhys);

    return math::max( math::abs(L.local[0](w,0) - (parametric ? gWestParam : gWestPhys)),
                      math::abs(L.local[0](e,0) - (parametric ? gEastParam : gEastPhys)) );
}

/// Coefficients of the legacy gsPoissonAssembler (see legacyDirichletCoefs) with
/// the DirichletValues option set to \a method, or left at its default if
/// \a method is negative.
void legacyCoefsWithMethod(const gsBasis<real_t> & basis, const gsMultiPatch<real_t> & mp,
                           const gsFunction<real_t> & data, bool parametric, index_t method,
                           gsMatrix<real_t> & vals, std::string & out, std::string & failure)
{
    gsMultiBasis<real_t> mb(basis);
    mb.setTopology(mp);
    gsBoundaryConditions<real_t> bc;
    bc.setGeoMap(mp);
    for (boxSide s = boxSide::getFirst(2); s < boxSide::getEnd(2); ++s)
        bc.addCondition(0, s, condition_type::dirichlet, data, 0, parametric, -1);
    const gsFunctionExpr<real_t> f("0", 2);

    try
    {
        CoutCapture capture;
        gsPoissonAssembler<real_t> PA(mp, mb, bc, f);
        if (method >= 0) PA.options().setInt("DirichletValues", method);
        PA.assemble();

        gsSparseSolver<real_t>::SimplicialLDLT solver;
        solver.compute( PA.matrix() );
        const gsMatrix<real_t> solVector = solver.solve( PA.rhs() );

        gsMultiPatch<real_t> sol;
        PA.constructSolution(solVector, sol);
        out = capture.release();
        vals = sol.patch(0).coefs();
    }
    catch (const std::exception & e) { failure = e.what(); }
    catch (...)                      { failure = "non-standard exception"; }
}

/// Free coefficients of the L2 projection of \a f onto \a mb with the Dirichlet
/// condition \a bc, assembled exactly as gsProjection<ProjectionNorm::L2>::project
/// does but with the Dirichlet values computed by \a method.
gsMatrix<real_t> projectionWithMethod(const gsMultiBasis<real_t> & mb, const gsMultiPatch<real_t> & mp,
                                      const gsFunction<real_t> & f, const gsBoundaryConditions<real_t> & bc,
                                      index_t method)
{
    gsExprAssembler<real_t> A(1,1);
    A.setIntegrationDomain(mb.domain());
    gsExprAssembler<real_t>::space u = A.getSpace(mb, f.targetDim());
    auto ff = A.getCoeff(f);
    gsExprAssembler<real_t>::geometryMap G = A.getMap(mp);
    u.setup(bc, method, -1);
    A.initSystem();
    A.assemble( 1.0 * (u * u.tr()) * meas(G) );
    A.assemble( 1.0 * (u * ff) * meas(G) );
    gsSparseSolver<real_t>::uPtr solver = gsSparseSolver<real_t>::get("SimplicialLDLT");
    solver->compute(A.matrix());
    return solver->solve(A.rhs());
}

/// The clamped 1-D bases of the end-point tests, with the label of each: B-spline,
/// NURBS with non-unit weights, and THB and HB bases of the B-spline refined at
/// the west or at the east end.
std::vector< std::pair<std::string, gsBasis<real_t>::uPtr> > oneDimBases()
{
    gsKnotVector<real_t> kv(0.0, 1.0, 3, 3);
    gsBSplineBasis<real_t> bs(kv);
    gsMatrix<real_t> w(bs.size(), 1);
    for (index_t i = 0; i != w.rows(); ++i) w(i,0) = 1 + 0.5*i;

    gsMatrix<real_t> boxWest(1,2), boxEast(1,2);
    boxWest << 0.0,  0.25;
    boxEast << 0.75, 1.0;

    std::vector< std::pair<std::string, gsBasis<real_t>::uPtr> > out;
    out.push_back( std::make_pair(std::string("bspl"),
                   gsBasis<real_t>::uPtr(new gsBSplineBasis<real_t>(bs))) );
    out.push_back( std::make_pair(std::string("nurbs"),
                   gsBasis<real_t>::uPtr(new gsNurbsBasis<real_t>(kv, w))) );
    for (int east = 0; east != 2; ++east)
    {
        gsTHBSplineBasis<1,real_t> thb(bs);
        thb.refine(east ? boxEast : boxWest);
        out.push_back( std::make_pair(std::string(east ? "thb_east" : "thb_west"),
                       gsBasis<real_t>::uPtr(new gsTHBSplineBasis<1,real_t>(thb))) );
    }
    for (int east = 0; east != 2; ++east)
    {
        gsHBSplineBasis<1,real_t> hb(bs);
        hb.refine(east ? boxEast : boxWest);
        out.push_back( std::make_pair(std::string(east ? "hb_east" : "hb_west"),
                       gsBasis<real_t>::uPtr(new gsHBSplineBasis<1,real_t>(hb))) );
    }
    return out;
}

/// Dirichlet dofs at the two ends of the 1-D basis \a b, set up through the
/// expression assembler with \a method for the data cos(x)+x^2 on a patch
/// with basis \a b and coefficients 2 + 3 anchors (see oneDimEndError).
Lift oneDimLift(const gsBasis<real_t> & b, bool parametric, index_t method)
{
    gsMatrix<real_t> coefs = (2 + 3 * b.anchors().array()).matrix().transpose();
    gsMultiPatch<real_t> mp;
    mp.addPatch( b.makeGeometry(coefs) );
    gsFunctionExpr<real_t> g("cos(x)+x^2", 1);
    gsBoundaryConditions<real_t> bc;
    bc.setGeoMap(mp);
    bc.addCondition(0, boundary::west, condition_type::dirichlet, &g, 0, parametric, -1);
    bc.addCondition(0, boundary::east, condition_type::dirichlet, &g, 0, parametric, -1);
    gsMultiBasis<real_t> mb(b);
    return liftDirichlet(mb, bc, method);
}

/// Dirichlet end values of the 1-D basis \a b on the patch \a mp through the
/// legacy gsPoissonAssembler with the DirichletValues option \a method, for
/// \a data on both ends. Returns true if the assembly threw; \a err then holds
/// what was written to std::cerr. Otherwise \a west and \a east are the end dofs.
bool legacyOneDimEnds(const gsBasis<real_t> & b, const gsMultiPatch<real_t> & mp,
                      const gsFunction<real_t> & data, bool parametric, index_t method,
                      real_t & west, real_t & east, std::string & err)
{
    gsMultiBasis<real_t> mb(b);
    gsBoundaryConditions<real_t> bc;
    bc.setGeoMap(mp);
    for (boxSide s = boxSide::getFirst(1); s < boxSide::getEnd(1); ++s)
        bc.addCondition(0, s, condition_type::dirichlet, data, 0, parametric, -1);
    const gsFunctionExpr<real_t> f("0", 1);

    bool threw = false;
    CoutCapture cout_capture;
    CerrCapture cerr_capture;
    try
    {
        gsPoissonAssembler<real_t> pa(mp, mb, bc, f);
        pa.options().setInt("DirichletValues", method);
        pa.refresh();
        pa.assemble();
        const gsDofMapper & map = pa.system().colMapper(0);
        const index_t iw = b.boundary(boundary::west)(0,0), ie = b.boundary(boundary::east)(0,0);
        west = pa.fixedDofs(0)(map.bindex(iw, 0, 0), 0);
        east = pa.fixedDofs(0)(map.bindex(ie, 0, 0), 0);
    }
    catch (const std::exception &) { threw = true; }
    err = cerr_capture.release();
    return threw;
}

} // anonymous namespace

// Quasi-interpolation of the trace of a tensor B-spline reproduces its
// coefficients, parametric and physical, in two and three dimensions.
TEST(quasiInterpolation_tensor_2d)
{
    CornerBoxCase c;
    const gsMatrix<real_t> pc = spreadValues(c.tb.size(), 0.7, 0.3);
    gsGeometry<real_t>::uPtr spline = c.tb.makeGeometry(pc);

    for (int param = 1; param >= 0; --param)
    {
        std::string out; index_t n = 0;
        const real_t err = traceCoefError<2>(c.tb, c.mp, *spline, pc,
                                             dirichlet::quasiInterpolation, param != 0, out, n);
        CHECK( n > 0 );
        CHECK( err <= 1e-12 );
    }
}

TEST(quasiInterpolation_tensor_3d)
{
    gsKnotVector<real_t> kv(0.0, 1.0, 3, 4);
    std::vector< gsKnotVector<real_t> > kvs(3, kv);
    gsTensorBSplineBasis<3,real_t> tb(kvs);
    gsMultiPatch<real_t> mp;
    mp.addPatch( gsTensorBSpline<3,real_t>(tb, tb.anchors().transpose()) );
    const gsMatrix<real_t> pc = spreadValues(tb.size(), 0.7, 0.3);
    gsGeometry<real_t>::uPtr spline = tb.makeGeometry(pc);

    std::string out; index_t n = 0;
    const real_t err = traceCoefError<3>(tb, mp, *spline, pc,
                                         dirichlet::quasiInterpolation, true, out, n);
    CHECK( n > 0 );
    CHECK( err <= 1e-12 );
}

// On a tensor patch the method is quasi-interpolation and not anchor
// interpolation: for data outside the boundary space the two prescribed
// vectors differ, and the former is the bulk quasi-interpolant.
TEST(quasiInterpolation_tensor_isNotAnchorInterpolation)
{
    CornerBoxCase c;
    gsFunctionExpr<real_t> g("exp(y)*cos(3*y)", 2);

    const gsGeometrySlice<real_t> westTrace(&g, 0, 0.0);
    gsMatrix<real_t> oracle;
    gsQuasiInterpolate<real_t>::localIntpl(*c.tb.boundaryBasis(boundary::west), westTrace, oracle);

    std::string out;
    const gsMatrix<real_t> qi = prescribedOnSide(c.tb, c.mp, g, true, boundary::west,
                                                 dirichlet::quasiInterpolation, out);
    const gsMatrix<real_t> ai = prescribedOnSide(c.tb, c.mp, g, true, boundary::west,
                                                 dirichlet::interpolation, out);
    CHECK( qi.size() > 0 );
    CHECK( maxDiff(qi, oracle) <= 1e-12 );
    CHECK( maxDiff(qi, ai) > 1e-6 );
}

// A NURBS patch with non-unit weights: the trace of a NURBS spline with
// coefficients c has prescribed values c.
TEST(quasiInterpolation_nurbs_genuineWeights)
{
    AnnulusCase c;
    CHECK( c.mb.basis(0).isRational() );
    const gsTensorNurbsBasis<2,real_t> * nurbsBasis =
        dynamic_cast<const gsTensorNurbsBasis<2,real_t>*>(&c.mb.basis(0));
    CHECK( nurbsBasis != nullptr );
    if (nurbsBasis == nullptr) return;
    CHECK( (nurbsBasis->weights().array() - 1).abs().maxCoeff() > 1e-3 );

    const gsMatrix<real_t> pc = spreadValues(nurbsBasis->size(), 0.7, 0.3);
    gsGeometry<real_t>::uPtr spline = nurbsBasis->makeGeometry(pc);

    std::string out; index_t n = 0;
    const real_t err = traceCoefError<2>(*nurbsBasis, c.mp, *spline, pc,
                                         dirichlet::quasiInterpolation, true, out, n);
    CHECK( n > 0 );
    CHECK( err <= 1e-11 );
}

// HB in three dimensions and in physical coordinates, and the public routine
// gsDirichletValuesByQuasiInterpolation called directly.
TEST(quasiInterpolation_hb_3d_and_physical)
{
    {
        gsKnotVector<real_t> kv(0.0, 1.0, 3, 4);
        std::vector< gsKnotVector<real_t> > kvs(3, kv);
        gsTensorBSplineBasis<3,real_t> tb(kvs);
        gsMultiPatch<real_t> mp;
        mp.addPatch( gsTensorBSpline<3,real_t>(tb, tb.anchors().transpose()) );

        gsHBSplineBasis<3,real_t> hb(tb);
        gsMatrix<real_t> box(3,2);
        box << 0, 0.5,
               0, 0.5,
               0, 0.5;
        hb.refine(box);
        CHECK( boundaryGrown(hb, tb, boundary::west) );

        const gsMatrix<real_t> pc = spreadValues(hb.size(), 0.7, 0.3);
        gsGeometry<real_t>::uPtr spline = hb.makeGeometry(pc);
        std::string out; index_t n = 0;
        const real_t err = traceCoefError<3>(hb, mp, *spline, pc,
                                             dirichlet::quasiInterpolation, true, out, n);
        CHECK( n > 0 );
        CHECK( err <= 1e-11 );
    }

    CornerBoxCase c;
    const gsHBSplineBasis<2,real_t> hb = cornerBoxHB(c.tb);
    CHECK( boundaryGrown(hb, c.tb, boundary::west) );
    const gsMatrix<real_t> pc = spreadValues(hb.size(), 0.7, 0.3);
    gsGeometry<real_t>::uPtr spline = hb.makeGeometry(pc);
    {
        std::string out; index_t n = 0;
        const real_t err = traceCoefError<2>(hb, c.mp, *spline, pc,
                                             dirichlet::quasiInterpolation, false, out, n);
        CHECK( n > 0 );
        CHECK( err <= 1e-11 );
    }

    // The enumerated method and the routine called on a zeroed unknown agree
    // bitwise, since every coefficient is computed independently.
    gsMultiBasis<real_t> mbHB(hb);
    mbHB.setTopology(c.mp);
    gsBoundaryConditions<real_t> bc;
    bc.setGeoMap(c.mp);
    addAllSides(bc, 0, *spline, true);

    gsExprAssembler<real_t> A(1,1);
    A.setIntegrationElements(mbHB);
    gsExprAssembler<real_t>::space u = A.getSpace(mbHB, 1, 0);
    u.setup(bc, dirichlet::quasiInterpolation, 0);
    const gsMatrix<real_t> viaSetup = u.fixedPart();
    u.setup(bc, dirichlet::homogeneous, 0);
    CHECK_EQUAL( viaSetup.size(), u.fixedPart().size() );
    CHECK( u.fixedPart().cwiseAbs().maxCoeff() == 0 );
    gsDirichletValuesByQuasiInterpolation(u, bc);
    CHECK( viaSetup.size() > 0 );
    CHECK( viaSetup.cwiseAbs().maxCoeff() > 0.1 );
    CHECK( exactDifference(viaSetup, u.fixedPart()) == 0 );
}

// The HB trace of data outside the boundary space is the bulk quasi-interpolant
// on the boundary basis, which is not the l2 projection.
TEST(dirichletTrace_hb_outOfSpace_matchesQI)
{
    CornerBoxCase c;
    const gsHBSplineBasis<2,real_t> hb = cornerBoxHB(c.tb);
    gsFunctionExpr<real_t> g("exp(y)*cos(3*y)", 2);
    CHECK( boundaryGrown(hb, c.tb, boundary::west) );

    const gsGeometrySlice<real_t> westTrace(&g, 0, 0.0);
    gsMatrix<real_t> oracle;
    gsQuasiInterpolate<real_t>::localIntpl(*hb.boundaryBasis(boundary::west), westTrace, oracle);

    std::string out, outL2;
    const gsMatrix<real_t> fixed = prescribedOnSide(hb, c.mp, g, true, boundary::west,
                                                    dirichlet::quasiInterpolation, out);
    const gsMatrix<real_t> l2 = prescribedOnSide(hb, c.mp, g, true, boundary::west,
                                                 dirichlet::l2Projection, outL2);
    CHECK( fixed.size() > 0 );
    CHECK( maxDiff(fixed, oracle) <= 1e-12 );
    CHECK( maxDiff(l2, oracle) > 1e-6 );
}

// On a 1-D patch a side is an end point and the dof is the data value there.
TEST(quasiInterpolation_1d_parametric)
{
    gsKnotVector<real_t> kv(0.0, 1.0, 3, 4);
    gsBSplineBasis<real_t> b(kv);
    real_t separation = 0;
    CHECK( oneDimEndError(b, true, separation) <= 1e-13 );
    CHECK( separation > 0.1 );
}

TEST(quasiInterpolation_1d_physical)
{
    gsKnotVector<real_t> kv(0.0, 1.0, 3, 4);
    gsBSplineBasis<real_t> b(kv);
    real_t separation = 0;
    CHECK( oneDimEndError(b, false, separation) <= 1e-13 );
    CHECK( separation > 0.1 );
}

// 1-D NURBS (non-unit weights) and 1-D hierarchical bases refined at either
// end, taken as the boundary bases of 2-D bases.
TEST(quasiInterpolation_1d_nurbs_and_hierarchical)
{
    CornerBoxCase c;

    AnnulusCase an;
    const gsTensorNurbsBasis<2,real_t> * nurbs =
        dynamic_cast<const gsTensorNurbsBasis<2,real_t>*>(&an.mb.basis(0));
    CHECK( nurbs != nullptr );
    std::vector< gsBasis<real_t>::uPtr > bases;
    if (nurbs != nullptr)
    {
        bases.push_back( nurbs->boundaryBasis(boundary::west) );
        const gsNurbsBasis<real_t> * nb1 = dynamic_cast<const gsNurbsBasis<real_t>*>(bases.back().get());
        CHECK( nb1 != nullptr );
        if (nb1 != nullptr)
            CHECK( (nb1->weights().array() - 1).abs().maxCoeff() > 1e-3 );
    }

    gsMatrix<real_t> boxEast(2,2), boxWest(2,2);
    boxEast << 0.75, 1.0,
               0.0,  1.0;
    boxWest << 0.0,  0.25,
               0.0,  1.0;
    for (int hb = 0; hb != 2; ++hb)
        for (int east = 0; east != 2; ++east)
        {
            gsBasis<real_t>::uPtr h;
            if (hb)
            {
                gsHBSplineBasis<2,real_t> basis(c.tb);
                basis.refine(east ? boxEast : boxWest);
                h = basis.boundaryBasis(boundary::south);
            }
            else
            {
                gsTHBSplineBasis<2,real_t> basis(c.tb);
                basis.refine(east ? boxEast : boxWest);
                h = basis.boundaryBasis(boundary::south);
            }
            CHECK( h->size() > c.tb.boundaryBasis(boundary::south)->size() );
            bases.push_back( give(h) );
        }

    CHECK_EQUAL( 5u, bases.size() );
    for (size_t i = 0; i != bases.size(); ++i)
        for (int param = 0; param != 2; ++param)
        {
            real_t separation = 0;
            CHECK( oneDimEndError(*bases[i], param != 0, separation) <= 1e-13 );
        }
}

// On every clamped 1-D basis (B-spline, NURBS, THB and HB refined at either
// end) dirichlet::interpolation and dirichlet::automatic give the data value at
// the end point, parametric and physical, and the same end dofs as
// quasi-interpolation.
TEST(interpolation_automatic_1d_allBases)
{
    const std::vector< std::pair<std::string, gsBasis<real_t>::uPtr> > bases = oneDimBases();
    CHECK_EQUAL( 6u, bases.size() );

    index_t visited = 0;
    for (size_t i = 0; i != bases.size(); ++i)
    {
        const gsBasis<real_t> & b = *bases[i].second;
        const std::string & name = bases[i].first;
        const bool west = name.find("west") != std::string::npos;
        const bool east = name.find("east") != std::string::npos;
        if (west || east)
        {
            const gsHTensorBasis<1,real_t> * h = dynamic_cast<const gsHTensorBasis<1,real_t>*>(&b);
            CHECK( h != nullptr );
            if (h != nullptr)
            {
                CHECK( b.size() > 6 );
                const index_t end = b.boundary(west ? boundary::west : boundary::east)(0,0);
                CHECK_EQUAL( 1, h->levelOf(end) );
            }
        }

        for (int param = 0; param != 2; ++param)
        {
            const index_t errMethods[2] = { dirichlet::interpolation, dirichlet::automatic };
            for (int m = 0; m != 2; ++m)
            {
                real_t separation = 0;
                CHECK( oneDimEndError(b, param != 0, separation, errMethods[m]) <= 1e-13 );
                CHECK( separation > 0.1 );
            }

            const index_t w = b.boundary(boundary::west)(0,0), e = b.boundary(boundary::east)(0,0);
            gsMatrix<real_t> ends[3];
            const index_t methods[3] = { dirichlet::interpolation, dirichlet::quasiInterpolation,
                                         dirichlet::automatic };
            for (int m = 0; m != 3; ++m)
            {
                const Lift L = oneDimLift(b, param != 0, methods[m]);
                ends[m].resize(2,1);
                ends[m] << L.local[0](w,0), L.local[0](e,0);
            }
            CHECK( exactDifference(ends[0], ends[1]) == 0 );
            CHECK( exactDifference(ends[2], ends[1]) == 0 );
        }
        ++visited;
    }
    CHECK_EQUAL( 6, visited );
}

// The legacy gsPoissonAssembler gives the same end dofs for the interpolation
// and automatic Dirichlet methods on every clamped 1-D basis.
TEST(legacy_interpolation_automatic_1d_allBases)
{
    const std::vector< std::pair<std::string, gsBasis<real_t>::uPtr> > bases = oneDimBases();
    CHECK_EQUAL( 6u, bases.size() );

    gsFunctionExpr<real_t> g("cos(x)+x^2", 1);
    index_t runs = 0;
    for (size_t i = 0; i != bases.size(); ++i)
    {
        const gsBasis<real_t> & b = *bases[i].second;
        const gsMatrix<real_t> coefs = (2 + 3 * b.anchors().array()).matrix().transpose();
        gsMultiPatch<real_t> mp;
        mp.addPatch( b.makeGeometry(coefs) );
        const index_t w = b.boundary(boundary::west)(0,0), e = b.boundary(boundary::east)(0,0);
        const real_t cw = coefs(w,0), ce = coefs(e,0);

        for (int param = 0; param != 2; ++param)
            for (int m = 0; m != 2; ++m)
            {
                const index_t method = m ? dirichlet::automatic : dirichlet::interpolation;
                real_t west = 0, east = 0;
                std::string err;
                const bool threw = legacyOneDimEnds(b, mp, g, param != 0, method, west, east, err);
                CHECK( !threw );
                const real_t expW = param ? math::cos(0.0) : math::cos(cw) + cw*cw;
                const real_t expE = param ? math::cos(1.0) + 1.0 : math::cos(ce) + ce*ce;
                CHECK( math::abs(west - expW) <= 1e-13 );
                CHECK( math::abs(east - expE) <= 1e-13 );
                ++runs;
            }
    }
    CHECK_EQUAL( 24, runs );
}

// The L2 projection reproduces the data at both ends of a 1-D basis: for the
// linear data 2+3x on x -> 2x+1 the end values are 2 and 5 (parametric), 5 and
// 11 (physical). B-spline, NURBS with non-unit weights and THB.
TEST(l2Projection_1d_endValues_exact)
{
    gsKnotVector<real_t> kv(0.0, 1.0, 3, 3);
    gsBSplineBasis<real_t> bs(kv);
    gsMatrix<real_t> wts(bs.size(), 1);
    for (index_t i = 0; i != wts.rows(); ++i) wts(i,0) = 1 + 0.5*i;
    gsNurbsBasis<real_t> nb(kv, wts);
    gsTHBSplineBasis<1,real_t> thb(bs);
    gsMatrix<real_t> box(1,2);
    box << 0, 0.25;
    thb.refine(box);

    gsKnotVector<real_t> lkv(0.0, 1.0, 0, 2);
    gsBSplineBasis<real_t> lb(lkv);
    gsMatrix<real_t> c(2,1);
    c << 1, 3;
    gsMultiPatch<real_t> mp;
    mp.addPatch( gsBSpline<real_t>(lb, c) );
    gsFunctionExpr<real_t> g("2+3*x", 1);

    const gsBasis<real_t> * bases[3] = { &bs, &nb, &thb };
    index_t runs = 0;
    for (int k = 0; k != 3; ++k)
        for (int param = 0; param != 2; ++param)
        {
            gsBoundaryConditions<real_t> bc;
            bc.setGeoMap(mp);
            bc.addCondition(0, boundary::west, condition_type::dirichlet, &g, 0, param != 0, -1);
            bc.addCondition(0, boundary::east, condition_type::dirichlet, &g, 0, param != 0, -1);
            gsMultiBasis<real_t> mb(*bases[k]);
            const Lift L = liftDirichlet(mb, bc, dirichlet::l2Projection);

            const real_t expW = param ? 2.0 : 5.0, expE = param ? 5.0 : 11.0;
            const index_t w = bases[k]->boundary(boundary::west)(0,0);
            const index_t e = bases[k]->boundary(boundary::east)(0,0);
            CHECK( math::abs(L.local[0](w,0) - expW) <= 1e-12 * expW );
            CHECK( math::abs(L.local[0](e,0) - expE) <= 1e-12 * expE );
            ++runs;
        }
    CHECK_EQUAL( 6, runs );
}

// Same end values through the legacy gsPoissonAssembler.
TEST(legacy_l2Projection_1d_endValues_exact)
{
    gsKnotVector<real_t> kv(0.0, 1.0, 3, 3);
    gsBSplineBasis<real_t> bs(kv);
    gsMatrix<real_t> wts(bs.size(), 1);
    for (index_t i = 0; i != wts.rows(); ++i) wts(i,0) = 1 + 0.5*i;
    gsNurbsBasis<real_t> nb(kv, wts);
    gsTHBSplineBasis<1,real_t> thb(bs);
    gsMatrix<real_t> box(1,2);
    box << 0, 0.25;
    thb.refine(box);

    gsKnotVector<real_t> lkv(0.0, 1.0, 0, 2);
    gsBSplineBasis<real_t> lb(lkv);
    gsMatrix<real_t> c(2,1);
    c << 1, 3;
    gsMultiPatch<real_t> mp;
    mp.addPatch( gsBSpline<real_t>(lb, c) );
    gsFunctionExpr<real_t> g("2+3*x", 1);

    const gsBasis<real_t> * bases[3] = { &bs, &nb, &thb };
    index_t runs = 0;
    for (int k = 0; k != 3; ++k)
        for (int param = 0; param != 2; ++param)
        {
            real_t west = 0, east = 0;
            std::string err;
            const bool threw = legacyOneDimEnds(*bases[k], mp, g, param != 0,
                                                dirichlet::l2Projection, west, east, err);
            CHECK( !threw );
            const real_t expW = param ? 2.0 : 5.0, expE = param ? 5.0 : 11.0;
            CHECK( math::abs(west - expW) <= 1e-12 * expW );
            CHECK( math::abs(east - expE) <= 1e-12 * expE );
            ++runs;
        }
    CHECK_EQUAL( 6, runs );
}

// Two disconnected patches, a tensor and a THB one: every Dirichlet dof of
// either patch is the coefficient of the data spline over its own basis.
TEST(quasiInterpolation_mixed_tensor_thb_noInterface)
{
    MixedTensorTHBCase c;
    CHECK( c.mp.nInterfaces() == 0 );
    CHECK( boundaryGrown(c.thb, c.tb, boundary::west) );

    const gsMatrix<real_t> pc0 = spreadValues(c.tb.size(), 0.7, 0.3);
    const gsMatrix<real_t> pc1 = spreadValues(c.thb.size(), 0.7, 0.3);
    gsGeometry<real_t>::uPtr s0 = c.tb.makeGeometry(pc0);
    gsGeometry<real_t>::uPtr s1 = c.thb.makeGeometry(pc1);

    gsMultiPatch<real_t> data;
    data.addPatch( give(s0) );
    data.addPatch( give(s1) );

    gsBoundaryConditions<real_t> bc;
    bc.setGeoMap(c.mp);
    addAllSides(bc, 0, data, true);
    addAllSides(bc, 1, data, true);
    const Lift L = liftDirichlet(c.mb, bc, dirichlet::quasiInterpolation);

    const gsMatrix<real_t> * pcs[2] = { &pc0, &pc1 };
    const gsBasis<real_t> * bs[2]   = { &c.tb, &c.thb };
    for (index_t k = 0; k != 2; ++k)
    {
        index_t n = 0;
        real_t err = 0;
        for (boxSide s = boxSide::getFirst(2); s < boxSide::getEnd(2); ++s)
        {
            const gsMatrix<index_t> idx = bs[k]->boundary(s);
            for (index_t l = 0; l != idx.rows(); ++l, ++n)
                err = math::max(err, math::abs(L.local[k](idx(l,0),0) - (*pcs[k])(idx(l,0),0)));
        }
        CHECK( n > 0 );
        CHECK( L.numFixed(k) > 0 );
        CHECK( err <= 1e-12 );
    }
}

// A dof shared by two Dirichlet sides keeps the value of the side processed
// last, i.e. the one added to the boundary conditions last.
TEST(quasiInterpolation_tensor_corner_lastWriteWins)
{
    CornerBoxCase c;
    gsFunctionExpr<real_t> g("exp(x+y)*cos(3*y+2*x)", 2);

    const gsMatrix<index_t> iw = c.tb.boundary(boundary::west);
    const gsMatrix<index_t> is = c.tb.boundary(boundary::south);
    index_t corner = -1, lw = -1, ls = -1, nCorner = 0;
    for (index_t l = 0; l != iw.rows(); ++l)
        for (index_t m = 0; m != is.rows(); ++m)
            if (iw(l,0) == is(m,0)) { corner = iw(l,0); lw = l; ls = m; ++nCorner; }
    CHECK_EQUAL( 1, nCorner );
    if (nCorner != 1) return;

    const gsGeometrySlice<real_t> westTrace(&g, 0, 0.0), southTrace(&g, 1, 0.0);
    gsMatrix<real_t> oW, oS;
    gsQuasiInterpolate<real_t>::localIntpl(*c.tb.boundaryBasis(boundary::west), westTrace, oW);
    gsQuasiInterpolate<real_t>::localIntpl(*c.tb.boundaryBasis(boundary::south), southTrace, oS);
    CHECK_EQUAL( iw.rows(), oW.rows() );
    CHECK_EQUAL( is.rows(), oS.rows() );
    CHECK( math::abs(oW(lw,0) - oS(ls,0)) > 1e-6 );

    gsMultiBasis<real_t> mbT(c.tb);
    for (int westFirst = 1; westFirst >= 0; --westFirst)
    {
        gsBoundaryConditions<real_t> bc;
        bc.setGeoMap(c.mp);
        if (westFirst)
        {
            bc.addCondition(0, boundary::west , condition_type::dirichlet, &g, 0, true, -1);
            bc.addCondition(0, boundary::south, condition_type::dirichlet, &g, 0, true, -1);
        }
        else
        {
            bc.addCondition(0, boundary::south, condition_type::dirichlet, &g, 0, true, -1);
            bc.addCondition(0, boundary::west , condition_type::dirichlet, &g, 0, true, -1);
        }
        const Lift L = liftDirichlet(mbT, bc, dirichlet::quasiInterpolation);

        const real_t winner = westFirst ? oS(ls,0) : oW(lw,0);
        const real_t loser  = westFirst ? oW(lw,0) : oS(ls,0);
        CHECK( math::abs(L.local[0](corner,0) - winner) <= 1e-12 );
        CHECK( math::abs(L.local[0](corner,0) - loser)  > 1e-6 );

        real_t err = 0;
        for (index_t l = 0; l != iw.rows(); ++l)
            if (l != lw) err = math::max(err, math::abs(L.local[0](iw(l,0),0) - oW(l,0)));
        for (index_t m = 0; m != is.rows(); ++m)
            if (m != ls) err = math::max(err, math::abs(L.local[0](is(m,0),0) - oS(m,0)));
        CHECK( err <= 1e-12 );
    }
}

// dirichlet::automatic is anchor interpolation on tensor and NURBS patches:
// bitwise the values of dirichlet::interpolation, and not those of
// quasi-interpolation.
TEST(automatic_tensor_equals_interpolation)
{
    CornerBoxCase c;
    gsMultiBasis<real_t> mbT(c.tb);
    mbT.setTopology(c.mp);
    gsFunctionExpr<real_t> g("exp(x)*cos(3*y)", 2);
    gsBoundaryConditions<real_t> bc;
    bc.setGeoMap(c.mp);
    addAllSides(bc, 0, g, false);

    const Lift a = liftDirichlet(mbT, bc, dirichlet::automatic);
    const Lift i = liftDirichlet(mbT, bc, dirichlet::interpolation);
    const Lift q = liftDirichlet(mbT, bc, dirichlet::quasiInterpolation);
    CHECK( a.fixed.size() > 0 );
    CHECK( exactDifference(a.fixed, i.fixed) == 0 );
    CHECK( exactDifference(q.fixed, i.fixed) > 1e-6 );
}

TEST(automatic_nurbs_equals_interpolation)
{
    AnnulusCase c;
    CHECK( c.mb.basis(0).isRational() );
    gsFunctionExpr<real_t> g("exp(x)*cos(3*y)", 2);
    gsBoundaryConditions<real_t> bc;
    bc.setGeoMap(c.mp);
    addAllSides(bc, 0, g, false);

    const Lift a = liftDirichlet(c.mb, bc, dirichlet::automatic);
    const Lift i = liftDirichlet(c.mb, bc, dirichlet::interpolation);
    const Lift q = liftDirichlet(c.mb, bc, dirichlet::quasiInterpolation);
    CHECK( a.fixed.size() > 0 );
    CHECK( exactDifference(a.fixed, i.fixed) == 0 );
    CHECK( exactDifference(q.fixed, i.fixed) > 1e-6 );
}

// dirichlet::automatic is quasi-interpolation on THB and HB patches.
TEST(automatic_thb_hb_equals_quasiInterpolation)
{
    CornerBoxCase c;
    gsFunctionExpr<real_t> g("exp(x)*cos(3*y)", 2);
    gsMultiBasis<real_t> mbHB(cornerBoxHB(c.tb));
    mbHB.setTopology(c.mp);

    const gsMultiBasis<real_t> * mbs[2] = { c.mb.get(), &mbHB };
    for (int k = 0; k != 2; ++k)
    {
        const gsBoundaryConditions<real_t> bc = allSidesDirichlet(c.mp, g, 2);
        const Lift a = liftDirichlet(*mbs[k], bc, dirichlet::automatic);
        const Lift q = liftDirichlet(*mbs[k], bc, dirichlet::quasiInterpolation);
        CHECK( a.fixed.size() > 0 );
        CHECK( a.fixed.cwiseAbs().maxCoeff() > 0.1 );
        CHECK( exactDifference(a.fixed, q.fixed) == 0 );
    }
}

// On a mixed multipatch dirichlet::automatic resolves per patch.
TEST(automatic_mixed_tensor_thb_perPatch)
{
    MixedTensorTHBCase c;
    gsFunctionExpr<real_t> g("exp(x)*cos(3*y)", 2);

    gsBoundaryConditions<real_t> bc;
    bc.setGeoMap(c.mp);
    addAllSides(bc, 0, g, true);
    addAllSides(bc, 1, g, true);
    const Lift mixed = liftDirichlet(c.mb, bc, dirichlet::automatic);

    gsMultiPatch<real_t> mpOne = identityPatch(c.tb);
    gsBoundaryConditions<real_t> bcOne;
    bcOne.setGeoMap(mpOne);
    addAllSides(bcOne, 0, g, true);
    gsMultiBasis<real_t> mbTb(c.tb), mbThb(c.thb);
    const Lift tbInterp = liftDirichlet(mbTb,  bcOne, dirichlet::interpolation);
    const Lift tbQI     = liftDirichlet(mbTb,  bcOne, dirichlet::quasiInterpolation);
    const Lift thbQI    = liftDirichlet(mbThb, bcOne, dirichlet::quasiInterpolation);

    CHECK( mixed.numFixed(0) > 0 );
    CHECK( mixed.numFixed(1) > 0 );
    CHECK( exactDifference(tbQI.local[0], tbInterp.local[0]) > 1e-6 );
    CHECK( exactDifference(mixed.local[0], tbInterp.local[0]) == 0 );
    CHECK( exactDifference(mixed.local[1], thbQI.local[0]) == 0 );
}

// A mapped basis is neither interpolated nor quasi-interpolated: automatic
// computes the l2 projection and says nothing about it.
TEST(automatic_mapped_equals_l2Projection_silently)
{
    gsKnotVector<> kv(0, 1, 3, 3);
    gsTensorBSplineBasis<2> tb(kv, kv);
    gsMultiPatch<> mp;
    mp.addPatch( *gsNurbsCreator<>::BSplineSquare(1.0, 0.0, 0.0) );
    mp.computeTopology();
    gsFunctionExpr<> g("exp(x)*cos(3*y)", 2);

    // Premise: on the tensor space of the mapped basis the three methods give
    // different north values, so that equality with l2 identifies it.
    std::string out;
    const gsMatrix<real_t> intpl = prescribedOnSide(tb, mp, g, true, boundary::north, dirichlet::interpolation, out);
    const gsMatrix<real_t> l2Ref = prescribedOnSide(tb, mp, g, true, boundary::north, dirichlet::l2Projection, out);
    const gsMatrix<real_t> qiRef = prescribedOnSide(tb, mp, g, true, boundary::north, dirichlet::quasiInterpolation, out);
    CHECK( maxDiff(intpl, l2Ref) > 1e-6 );
    CHECK( maxDiff(qiRef, l2Ref) > 1e-6 );

    gsMultiBasis<> mbM(tb);
    const index_t sz = mbM.basis(0).size();
    gsSparseMatrix<> ident(sz, sz);
    ident.setIdentity();
    gsMappedBasis<2, real_t> mapB(mbM, ident);

    gsExprAssembler<> AM(1,1);
    AM.setIntegrationElements(mbM);
    gsExprAssembler<>::space uM = AM.getSpace(mapB, 1, 0);
    gsBoundaryConditions<> bcM;
    bcM.setGeoMap(mp);
    bcM.addCondition(0, boundary::north, condition_type::dirichlet, &g, 0, true, -1);

    bool threw = false;
    std::string coutText, cerrText;
    {
        CoutCapture co;
        CerrCapture ce;
        try { uM.setup(bcM, dirichlet::automatic, 0); }
        catch (const std::exception &) { threw = true; }
        coutText = co.release();
        cerrText = ce.release();
    }
    CHECK( !threw );
    CHECK_EQUAL( 0, countOccurrences(coutText, fallbackWarning) );
    CHECK_EQUAL( 0, countOccurrences(cerrText, fallbackWarning) );
    CHECK_EQUAL( 0, countOccurrences(coutText, "arning") );
    CHECK_EQUAL( 0, countOccurrences(cerrText, "arning") );
    const gsMatrix<real_t> automaticFixed = uM.fixedPart();

    uM.setup(bcM, dirichlet::l2Projection, 0);
    const gsMatrix<real_t> l2Fixed = uM.fixedPart();
    CHECK( automaticFixed.size() > 0 );
    CHECK_EQUAL( automaticFixed.size(), l2Fixed.size() );
    if (automaticFixed.size() == l2Fixed.size())
        CHECK( (automaticFixed - l2Fixed).cwiseAbs().maxCoeff()
               <= 1e-12 * math::max((real_t)1, l2Fixed.cwiseAbs().maxCoeff()) );
}

TEST(automatic_is_default)
{
    CHECK_EQUAL( (index_t)dirichlet::automatic, gsAssembler<real_t>::defaultOptions().getInt("DirichletValues") );
    CHECK_EQUAL( (index_t)dirichlet::automatic, gsPoissonAssembler<real_t>::defaultOptions().getInt("DirichletValues") );
    CHECK_EQUAL( (index_t)dirichlet::automatic, gsExprAssembler<real_t>::defaultOptions().getInt("DirichletValues") );
}

// The default of the legacy assembler accepts a THB patch (interpolation would
// throw) and reproduces data in the boundary space. The l2 projection is exact
// here as well (about 1e-15), so this test does not tell the two apart.
TEST(legacy_assembler_defaultOptions_thb)
{
    gsKnotVector<real_t> kv(0.0, 1.0, 3, 4);
    gsTensorBSplineBasis<2,real_t> tb(kv, kv);
    const gsMultiPatch<real_t> mp = identityPatch(tb);

    gsTHBSplineBasis<2,real_t> thb(tb);
    gsMatrix<real_t> box(2,2);
    box << 0.75, 1.0,
           0.5,  1.0;
    thb.refine(box);
    CHECK( boundaryGrown(thb, tb, boundary::east) );

    const gsMatrix<real_t> pc = spreadValues(thb.size(), 0.7, 0.3);
    gsGeometry<real_t>::uPtr spline = thb.makeGeometry(pc);

    gsMatrix<real_t> vals;
    std::string out, failure;
    legacyCoefsWithMethod(thb, mp, *spline, true, -1, vals, out, failure);

    index_t n = 0;
    CHECK_EQUAL( std::string(), failure );
    CHECK( vals.size() == pc.size() );
    const real_t err = (vals.size() == pc.size()) ? boundaryCoefError(thb, vals, pc, n) : -1;
    CHECK( n > 0 );
    CHECK( err >= 0 && err <= 1e-12 );
}

// Legacy assembler: dirichlet::automatic on a tensor patch gives the anchor
// interpolation values bitwise, which differ from quasi-interpolation.
TEST(legacy_assembler_automatic_tensor_equals_interpolation)
{
    CornerBoxCase c;
    gsFunctionExpr<real_t> g("exp(x)*cos(3*y)", 2);

    gsMatrix<real_t> vA, vI, vQ;
    std::string out, failure;
    legacyCoefsWithMethod(c.tb, c.mp, g, false, dirichlet::automatic, vA, out, failure);
    CHECK_EQUAL( std::string(), failure );
    legacyCoefsWithMethod(c.tb, c.mp, g, false, dirichlet::interpolation, vI, out, failure);
    CHECK_EQUAL( std::string(), failure );
    legacyCoefsWithMethod(c.tb, c.mp, g, false, dirichlet::quasiInterpolation, vQ, out, failure);
    CHECK_EQUAL( std::string(), failure );

    CHECK( vA.size() == c.tb.size() && vI.size() == c.tb.size() && vQ.size() == c.tb.size() );
    if (vA.size() == c.tb.size() && vI.size() == c.tb.size() && vQ.size() == c.tb.size())
    {
        index_t nAI = 0, nQI = 0;
        CHECK( boundaryCoefError(c.tb, vA, vI, nAI) == 0 );
        CHECK( nAI > 0 );
        CHECK( boundaryCoefError(c.tb, vQ, vI, nQI) > 1e-6 );
    }
}

// Projection onto a THB space with Dirichlet conditions needs no tensor
// structure, and a spline of the space is reproduced.
TEST(projection_dirichlet_thb)
{
    CornerBoxCase c;
    const gsMatrix<real_t> pc = spreadValues(c.thb.size(), 0.7, 0.3);
    gsGeometry<real_t>::uPtr spline = c.thb.makeGeometry(pc);

    gsBoundaryConditions<real_t> bc;
    bc.setGeoMap(c.mp);
    addAllSides(bc, 0, *spline, true);

    gsMatrix<real_t> coefs;
    real_t err = -1;
    std::string failure;
    try { err = gsProjection<ProjectionNorm::L2, real_t>::project(*c.mb, c.mp, *spline, coefs, bc); }
    catch (const std::exception & e) { failure = e.what(); }
    catch (...)                      { failure = "non-standard exception"; }

    CHECK_EQUAL( std::string(), failure );
    CHECK( coefs.size() > 0 );
    CHECK( err >= 0 && err <= 1e-10 );
}

// On a tensor basis the projection takes interpolation-equivalent Dirichlet
// values. Out-of-space data separates them from quasi-interpolation.
TEST(projection_dirichlet_tensor_usesInterpolation)
{
    CornerBoxCase c;
    gsMultiBasis<real_t> mbT(c.tb);
    mbT.setTopology(c.mp);
    gsFunctionExpr<real_t> g("exp(x)*sin(2*y)+y^3", 2);

    gsBoundaryConditions<real_t> bc;
    bc.setGeoMap(c.mp);
    bc.addCondition(0, boundary::west, condition_type::dirichlet, &g, 0, true, -1);

    gsMatrix<real_t> coefs;
    gsProjection<ProjectionNorm::L2, real_t>::project(mbT, c.mp, g, coefs, bc);

    const gsMatrix<real_t> refI = projectionWithMethod(mbT, c.mp, g, bc, dirichlet::interpolation);
    const gsMatrix<real_t> refQ = projectionWithMethod(mbT, c.mp, g, bc, dirichlet::quasiInterpolation);
    CHECK( refI.size() > 0 );
    CHECK_EQUAL( refI.size(), coefs.size() );
    CHECK( maxDiff(refI, refQ) > 1e-8 );

    // gsExprAssembler accumulates with OpenMP atomic updates, so two assemblies
    // of the same system agree only to round-off and not bitwise.
    const real_t diff = maxDiff(coefs, refI);
    gsInfo << "[projection_dirichlet_tensor_usesInterpolation] |coefs - refI|_max = " << diff << "\n";
    CHECK( diff <= 1e-12 * math::max((real_t)1, refI.cwiseAbs().maxCoeff()) );
}

// The guard tests below run last: on a violated guard they read out of
// bounds or dereference a null pointer. UnitTest++ reports a SIGSEGV as a
// test failure, but the process state after it is undefined, and SIGABRT is
// not caught at all.

// A scalar function prescribed for a two-component unknown on a THB patch is
// rejected instead of reading a second component that does not exist.
TEST(dirichletGuard_thb_componentMismatch)
{
    CornerBoxCase c;
    gsFunctionExpr<real_t> g("1 + 2*x + 3*y", 2);
    gsBoundaryConditions<real_t> bc = allSidesDirichlet(c.mp, g, 2);

    gsExprAssembler<real_t> A(1,1);
    A.setIntegrationElements(*c.mb);
    gsExprAssembler<real_t>::space u = A.getSpace(*c.mb, 2);

    std::string err;
    {
        CerrCapture capture;
        CHECK_THROW( u.setup(bc, dirichlet::quasiInterpolation, 0), std::runtime_error );
        err = capture.release();
    }
    CHECK( err.find("cannot supply component") != std::string::npos );
}

/// West condition on \a mb without a geometry map, set up with \a strategy.
/// The condition is given in physical coordinates unless \a parametric. Returns
/// whether the setup threw; \a err receives what it wrote to std::cerr.
bool setupWithoutGeoMap(const gsMultiBasis<real_t> & mb, dirichlet::values strategy,
                        bool parametric, std::string & err)
{
    gsFunctionExpr<real_t> g = cornerBoxSolution();
    gsBoundaryConditions<real_t> bc;
    bc.addCondition(0, boundary::west, condition_type::dirichlet, &g, 0, parametric, -1);

    gsExprAssembler<real_t> A(1,1);
    A.setIntegrationElements(mb);
    gsExprAssembler<real_t>::space u = A.getSpace(mb, 1, 0);

    bool threw = false;
    {
        CoutCapture out;
        CerrCapture capture;
        try { u.setup(bc, strategy, 0); }
        catch (const std::runtime_error &) { threw = true; }
        err = capture.release();
    }
    return threw;
}

const std::string missingGeoMapMessage = "call bc.setGeoMap(";

// A physical condition without bc.setGeoMap is reported by an exception that
// names the remedy.
TEST(dirichletGuard_missingGeoMap_thbQuasiInterpolation)
{
    CornerBoxCase c;
    std::string err;
    CHECK( setupWithoutGeoMap(*c.mb, dirichlet::quasiInterpolation, false, err) );
    CHECK( err.find(missingGeoMapMessage) != std::string::npos );
}

TEST(dirichletGuard_missingGeoMap_thbL2Projection)
{
    CornerBoxCase c;
    std::string err;
    CHECK( setupWithoutGeoMap(*c.mb, dirichlet::l2Projection, false, err) );
    CHECK( err.find(missingGeoMapMessage) != std::string::npos );
}

TEST(dirichletGuard_missingGeoMap_tensorInterpolation)
{
    CornerBoxCase c;
    std::string err;
    CHECK( setupWithoutGeoMap(gsMultiBasis<real_t>(c.tb), dirichlet::interpolation, false, err) );
    CHECK( err.find(missingGeoMapMessage) != std::string::npos );
}

// L2 projection on a non-truncated HB patch needs the geometry map for its
// measure even for parametric data. Parametric quasi-interpolation does not.
TEST(dirichletGuard_hb_missingGeoMap_parametric)
{
    CornerBoxCase c;
    gsHBSplineBasis<2,real_t> hb(c.tb);
    gsMatrix<real_t> box(2,2);
    box << 0, 0.5,
           0, 0.5;
    hb.refine(box);

    std::string err;
    CHECK( setupWithoutGeoMap(gsMultiBasis<real_t>(hb), dirichlet::l2Projection, true, err) );
    CHECK( err.find(missingGeoMapMessage) != std::string::npos );
}

namespace {

/// Non-truncated hierarchical basis with rational weights; only what the
/// Dirichlet routines inspect (the type of the basis) is implemented.
class RationalHB2 : public gsRationalBasis<gsHBSplineBasis<2,real_t> >
{
public:
    typedef gsRationalBasis<gsHBSplineBasis<2,real_t> > Base;
    RationalHB2(gsHBSplineBasis<2,real_t> * src, gsMatrix<real_t> w) : Base(src, give(w)) { }
    GISMO_CLONE_FUNCTION(RationalHB2)
    gsGeometry<real_t>::uPtr makeGeometry(gsMatrix<real_t>) const override { GISMO_NO_IMPLEMENTATION }
    std::ostream & print(std::ostream & os) const override { return os << "RationalHB2\n"; }
};

/// Condition on every side of every patch of \a mb with the data
/// exp(x)cos(3y) in physical coordinates, set up with \a method: the setup must
/// throw, and what it wrote to std::cerr must contain \a needle.
void expectSetupThrows(const gsMultiBasis<real_t> & mb, const gsMultiPatch<real_t> & mp,
                       index_t method, const std::string & needle)
{
    gsFunctionExpr<real_t> g("exp(x)*cos(3*y)", 2);
    gsBoundaryConditions<real_t> bc;
    bc.setGeoMap(mp);
    for (size_t k = 0; k != mb.nBases(); ++k)
        addAllSides(bc, k, g, false);

    gsExprAssembler<real_t> A(1,1);
    A.setIntegrationElements(mb);
    gsExprAssembler<real_t>::space u = A.getSpace(mb, 1, 0);

    std::string err;
    {
        CoutCapture out;
        CerrCapture capture;
        CHECK_THROW( u.setup(bc, method, 0), std::runtime_error );
        err = capture.release();
    }
    CHECK( err.find(needle) != std::string::npos );
}

} // anonymous namespace

// Anchor interpolation rejects THB and HB bases and names the remedy for them.
TEST(dirichletGuard_interpolation_hierarchical_namesQI)
{
    CornerBoxCase c;
    expectSetupThrows(*c.mb, c.mp, dirichlet::interpolation, qiHint);

    gsMultiBasis<real_t> mbHB(cornerBoxHB(c.tb));
    mbHB.setTopology(c.mp);
    expectSetupThrows(mbHB, c.mp, dirichlet::interpolation, qiHint);
}

// A multipatch with one unsupported patch is rejected as a whole.
TEST(dirichletGuard_interpolation_mixed_namesQI)
{
    MixedTensorTHBCase c;
    expectSetupThrows(c.mb, c.mp, dirichlet::interpolation, qiHint);
}

// A mapped basis is rejected by both interpolations, each naming the remedy
// available for it.
TEST(dirichletGuard_mapped_namesRemedy)
{
    gsKnotVector<> kv(0, 1, 3, 3);
    gsTensorBSplineBasis<2> tb(kv, kv);
    gsMultiPatch<> mp;
    mp.addPatch( *gsNurbsCreator<>::BSplineSquare(1.0, 0.0, 0.0) );
    mp.computeTopology();
    gsFunctionExpr<> g("exp(x)*cos(3*y)", 2);

    gsMultiBasis<> mbM(tb);
    const index_t sz = mbM.basis(0).size();
    gsSparseMatrix<> ident(sz, sz);
    ident.setIdentity();
    gsMappedBasis<2, real_t> mapB(mbM, ident);

    gsExprAssembler<> AM(1,1);
    AM.setIntegrationElements(mbM);
    gsExprAssembler<>::space uM = AM.getSpace(mapB, 1, 0);
    gsBoundaryConditions<> bcM;
    bcM.setGeoMap(mp);
    bcM.addCondition(0, boundary::north, condition_type::dirichlet, &g, 0, true, -1);

    for (int k = 0; k != 2; ++k)
    {
        const index_t method = (0 == k) ? dirichlet::interpolation : dirichlet::quasiInterpolation;
        const std::string & needle = (0 == k) ? qiHint : l2Hint;
        std::string err;
        {
            CoutCapture out;
            CerrCapture capture;
            CHECK_THROW( uM.setup(bcM, method, 0), std::runtime_error );
            err = capture.release();
        }
        CHECK( err.find(needle) != std::string::npos );
    }
}

// Quasi-interpolation does not cover a rational basis over HB; the remedy is
// the l2 projection.
TEST(dirichletGuard_quasiInterpolation_rationalHB_namesL2)
{
    CornerBoxCase c;
    gsHBSplineBasis<2,real_t> * hb = new gsHBSplineBasis<2,real_t>( cornerBoxHB(c.tb) );
    const gsMatrix<real_t> W = 1.25 + 0.75 * spreadValues(hb->size(), 1.3, 0.1).array();
    CHECK( (W.array() - 1).abs().maxCoeff() > 0.1 );
    RationalHB2 rb(hb, W);

    gsMultiBasis<real_t> mbR(rb);
    expectSetupThrows(mbR, c.mp, dirichlet::quasiInterpolation, l2Hint);
}

// The legacy assembler rejects interpolation on THB with the same hint, and
// quasi-interpolation on a rational HB basis likewise.
TEST(dirichletGuard_legacy_namesRemedy)
{
    CornerBoxCase c;
    const gsMatrix<real_t> pc = spreadValues(c.thb.size(), 0.7, 0.3);
    gsGeometry<real_t>::uPtr spline = c.thb.makeGeometry(pc);

    {
        gsMatrix<real_t> vals;
        std::string out, failure, err;
        {
            CerrCapture capture;
            legacyDirichletCoefs(c.thb, c.mp, *spline, true, dirichlet::interpolation, vals, out, failure);
            err = capture.release();
        }
        CHECK( !failure.empty() );
        CHECK( err.find(qiHint) != std::string::npos );
    }

    {
        gsHBSplineBasis<2,real_t> * hb = new gsHBSplineBasis<2,real_t>( cornerBoxHB(c.tb) );
        const gsMatrix<real_t> W = 1.25 + 0.75 * spreadValues(hb->size(), 1.3, 0.1).array();
        RationalHB2 rb(hb, W);

        const gsFunctionExpr<real_t> g("exp(x)*cos(3*y)", 2);
        gsMatrix<real_t> vals;
        std::string out, failure, err;
        {
            CerrCapture capture;
            legacyCoefsWithMethod(rb, c.mp, g, false, dirichlet::quasiInterpolation, vals, out, failure);
            err = capture.release();
        }
        CHECK( !failure.empty() );
        CHECK( err.find(l2Hint) != std::string::npos );
    }
}

// A 1-D basis whose end function is not 1 at the end point (unclamped knot
// vector) is rejected by quasi-interpolation.
TEST(dirichletGuard_quasiInterpolation_1d_unclamped)
{
    gsKnotVector<real_t> kv(0.0, 1.0, 5, 1, 1, 2);
    gsBSplineBasis<real_t> b(kv);
    gsMultiBasis<real_t> mb(b);
    gsFunctionExpr<real_t> g("cos(x)+x^2", 1);
    gsBoundaryConditions<real_t> bc;
    bc.addCondition(0, boundary::west, condition_type::dirichlet, &g, 0, true, -1);

    gsExprAssembler<real_t> A(1,1);
    A.setIntegrationElements(mb);
    gsExprAssembler<real_t>::space u = A.getSpace(mb, 1, 0);

    std::string err;
    {
        CoutCapture out;
        CerrCapture capture;
        CHECK_THROW( u.setup(bc, dirichlet::quasiInterpolation, 0), std::runtime_error );
        err = capture.release();
    }
    CHECK( err.find("open (clamped) knot vector") != std::string::npos );
}

// A 1-D basis whose end function is not 1 at the end point (unclamped knot
// vector) is rejected by dirichlet::interpolation and dirichlet::automatic, in
// the expression assembler and in the legacy assembler.
TEST(dirichletGuard_interpolation_automatic_1d_unclamped)
{
    gsKnotVector<real_t> kv(0.0, 1.0, 5, 1, 1, 2);
    gsBSplineBasis<real_t> b(kv);
    gsMultiBasis<real_t> mb(b);
    gsFunctionExpr<real_t> g("cos(x)+x^2", 1);
    const std::string hint = "open (clamped) knot vector";

    const gsMatrix<real_t> coefs = (2 + 3 * b.anchors().array()).matrix().transpose();
    gsMultiPatch<real_t> mp;
    mp.addPatch( b.makeGeometry(coefs) );

    index_t checked = 0;
    for (int m = 0; m != 2; ++m)
    {
        const index_t method = m ? dirichlet::automatic : dirichlet::interpolation;

        gsBoundaryConditions<real_t> bc;
        bc.addCondition(0, boundary::west, condition_type::dirichlet, &g, 0, true, -1);
        gsExprAssembler<real_t> A(1,1);
        A.setIntegrationElements(mb);
        gsExprAssembler<real_t>::space u = A.getSpace(mb, 1, 0);
        std::string err;
        {
            CoutCapture out;
            CerrCapture capture;
            CHECK_THROW( u.setup(bc, method, 0), std::runtime_error );
            err = capture.release();
        }
        CHECK( err.find(hint) != std::string::npos );

        real_t west = 0, east = 0;
        std::string legacyErr;
        CHECK( legacyOneDimEnds(b, mp, g, true, method, west, east, legacyErr) );
        CHECK( legacyErr.find(hint) != std::string::npos );
        ++checked;
    }
    CHECK_EQUAL( 2, checked );
}

}
