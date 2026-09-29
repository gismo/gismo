/** @file gsDirichletValues_test.cpp

    @brief Tests for strongly imposed Dirichlet values.

    The prescribed values produced by dirichlet::interpolation are compared
    against those produced by dirichlet::l2Projection on a manufactured
    solution of degree one. A polynomial of degree one lies in every spline
    space used here, so the trace lies in the boundary space and both
    strategies are exact: their coefficients must agree to round-off, and the
    Galerkin solution of -Laplace(u) = 0 must reproduce the exact solution.
    Any deviation is a defect in the construction of the prescribed values and
    not a discretization error.

    Covers regressions of two defects in gsDirichletValuesByTPInterpolation:

      - the interpolation nodes were taken from a tensor grid of
        basis.component(i).anchors(), which does not match the size of the
        boundary basis once the basis is hierarchical;
      - the prescribed values were read with linear (column-major) indexing,
        so every component of a vector-valued unknown received the first
        component of the boundary function.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst
*/

#include "gismo_unittest.h"
#include <gsMSplines/gsMappedBasis.h>

namespace {

const index_t degree = 2;
const index_t numRef = 2;
const real_t  tol    = 1e-10;

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
/// dirichlet::interpolation and with dirichlet::l2Projection, for an unknown
/// with \a nComp components.
template <short_t d>
real_t prescribedValueDifference(const gsMultiPatch<real_t> & mp,
                                 const gsMultiBasis<real_t> & dbasis,
                                 index_t nComp)
{
    gsFunctionExpr<real_t> g = exactSolution(d, nComp);
    gsBoundaryConditions<real_t> bc = allSidesDirichlet(mp, g, d);

    gsExprAssembler<real_t> A(1,1);
    A.setIntegrationDomain(dbasis.domain());
    gsExprAssembler<real_t>::geometryMap G = A.getMap(mp);
    GISMO_UNUSED(G);
    gsExprAssembler<real_t>::space u = A.getSpace(dbasis, nComp);

    u.setup(bc, dirichlet::interpolation, 0);
    const gsMatrix<real_t> interp = u.fixedPart();
    u.setup(bc, dirichlet::l2Projection, 0);
    const gsMatrix<real_t> proj = u.fixedPart();

    // A vacuous comparison of two empty vectors must not be read as a pass
    CHECK(interp.size() > 0);
    CHECK_EQUAL(interp.size(), proj.size());

    return (interp - proj).cwiseAbs().maxCoeff()
         / math::max( (real_t)1, proj.cwiseAbs().maxCoeff() );
}

/// L2 error of the Galerkin solution of -Laplace(u) = 0 against the exact
/// solution of degree one, with the Dirichlet values built by \a strategy.
template <short_t d>
real_t poissonError(const gsMultiPatch<real_t> & mp,
                    const gsMultiBasis<real_t> & dbasis,
                    dirichlet::values strategy)
{
    gsFunctionExpr<real_t> g = exactSolution(d, 1);
    gsFunctionExpr<real_t> zero("0", d);
    gsBoundaryConditions<real_t> bc = allSidesDirichlet(mp, g, d);

    gsExprAssembler<real_t> A(1,1);
    A.setIntegrationDomain(dbasis.domain());
    gsExprEvaluator<real_t> ev(A);

    gsExprAssembler<real_t>::geometryMap G = A.getMap(mp);
    gsExprAssembler<real_t>::space u = A.getSpace(dbasis);
    auto f = A.getCoeff(zero, G);
    auto g_ex = ev.getVariable(g, G);

    gsMatrix<real_t> solVector;
    gsExprAssembler<real_t>::solution u_sol = A.getSolution(u, solVector);

    u.setup(bc, strategy, 0);
    A.initSystem();
    A.assemble( igrad(u,G) * igrad(u,G).tr() * meas(G), u * f * meas(G) );

    gsSparseSolver<real_t>::SimplicialLDLT solver;
    solver.compute( A.matrix() );
    solVector = solver.solve( A.rhs() );

    return math::sqrt( ev.integral( (u_sol-g_ex).sqNorm() * meas(G) ) );
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

    // The boundary basis must number its functions like basis.boundary(side),
    // which is what the prescribed values are scattered with.
    for (index_t s=1; s<=2*d; ++s)
        CHECK_EQUAL( dbasis.basis(0).boundary(boxSide(s)).size(),
                     dbasis.basis(0).boundaryBasis(boxSide(s))->size() );

    // Scalar unknown
    CHECK_CLOSE( 0.0, prescribedValueDifference<d>(mp, dbasis, 1), tol );
    CHECK_CLOSE( 0.0, poissonError<d>(mp, dbasis, dirichlet::interpolation), tol );
    CHECK_CLOSE( 0.0, poissonError<d>(mp, dbasis, dirichlet::l2Projection ), tol );

    // Vector-valued unknown: component r must be taken from component r of the
    // boundary function
    CHECK_CLOSE( 0.0, prescribedValueDifference<d>(mp, dbasis, d), tol );
}

} // anonymous namespace

SUITE(gsDirichletValues_test)
{

TEST(tensor_2d)  { runCase<2>(false); }
TEST(hierarchical_2d) { runCase<2>(true ); }
TEST(tensor_3d)  { runCase<3>(false); }
TEST(hierarchical_3d) { runCase<3>(true ); }

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

// gsDirichletValuesByTPInterpolation's basis guard admits tensor and
// hierarchical bases and their rational counterparts: the boundary
// indices, face anchors and interpolateAtAnchors all resolve through
// gsRationalBasis::source(). Every admitted basis must produce the same
// prescribed values as an equivalent reference basis; any other basis
// must be rejected by the guard itself.
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
                try { u.setup(bc, dirichlet::interpolation, 0); }
                catch (const std::runtime_error &) { threw = true; }
                gsInfo << "[dirichletNurbs] case1 threw = "<<threw<<"\n";
                CHECK(!threw);

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
    // strategy; threw reports whether setup raised.
    auto northValues = [&](const gsMultiBasis<> & mbX, const index_t method, bool & threw)
    {
        gsExprAssembler<> AX(1,1);
        AX.setIntegrationElements(mbX);
        gsExprAssembler<>::space uX = AX.getSpace(mbX, 1, 0);
        gsBoundaryConditions<> bcX;
        bcX.setGeoMap(mp);
        bcX.addCondition(0, boundary::north, condition_type::dirichlet, &g, 0, true, -1);
        threw = false;
        try { uX.setup(bcX, method, 0); }
        catch (const std::runtime_error &) { threw = true; }
        return gsMatrix<>(uX.fixedPart());
    };

    auto sameValues = [](const gsMatrix<> & a, const gsMatrix<> & b)
    {
        return a.rows() == b.rows() && a.cols() == b.cols() &&
               (a - b).cwiseAbs().maxCoeff() < 1e-8;
    };

    bool threwB = false;
    const gsMatrix<> valsB = northValues(gsMultiBasis<>(tb), dirichlet::interpolation, threwB);
    CHECK(!threwB);
    CHECK( valsB.size() > 0 && valsB.maxCoeff() - valsB.minCoeff() > 0.1 );

    // Case 2: unit-weight NURBS spans the tensor B-spline space, and
    // interpolation is unique, so the rational path (generic collocation
    // + BiCGSTABILUT) must match the tensor path to solver tolerance.
    {
        gsTensorNurbsBasis<2> nb(tb);
        bool threwN = false;
        const gsMatrix<> valsN = northValues(gsMultiBasis<>(nb), dirichlet::interpolation, threwN);
        gsInfo << "[dirichletGuard] case2 threwN = "<<threwN<<"\n";
        CHECK(!threwN);
        if (!threwB && !threwN) CHECK( sameValues(valsB, valsN) );
    }

    // Case 3: a single-level THB basis spans the tensor space.
    gsTHBSplineBasis<2> thb(tb);
    {
        bool threwH = false;
        const gsMatrix<> valsH = northValues(gsMultiBasis<>(thb), dirichlet::interpolation, threwH);
        gsInfo << "[dirichletGuard] case3 threwH = "<<threwH<<"\n";
        CHECK(!threwH);
        if (!threwB && !threwH) CHECK( sameValues(valsB, valsH) );
    }

    // Case 4: THB refined in [0,0.5]x[0.5,1], so the north boundary basis
    // is of mixed level. Interpolation must match L2 projection, and the
    // unit-weight rational THB basis must match the THB basis.
    thb.refine( (gsMatrix<>(2,2) << 0, 0.5, 0.5, 1).finished() );
    CHECK( thb.boundaryBasis(boundary::north)->size() > tb.boundaryBasis(boundary::north)->size() );
    {
        bool threwI = false, threwL = false, threwR = false;
        const gsMatrix<> valsI = northValues(gsMultiBasis<>(thb), dirichlet::interpolation, threwI);
        const gsMatrix<> valsL = northValues(gsMultiBasis<>(thb), dirichlet::l2Projection , threwL);
        gsRationalTHBSplineBasis<2> rthb(thb);
        const gsMatrix<> valsR = northValues(gsMultiBasis<>(rthb), dirichlet::interpolation, threwR);
        gsInfo << "[dirichletGuard] case4 threwI = "<<threwI<<", threwL = "<<threwL
               << ", threwR = "<<threwR<<"\n";
        CHECK(!threwI);
        CHECK(!threwL);
        CHECK(!threwR);
        if (!threwI && !threwL) CHECK( sameValues(valsI, valsL) );
        if (!threwI && !threwR) CHECK( sameValues(valsI, valsR) );
    }

    // Case 5: a mapped basis is neither tensor nor hierarchical and must
    // be rejected by the guard, not by some later failure.
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

        // GISMO_ENSURE throws the literal "GISMO_ENSURE" and writes its
        // message to std::cerr, so the message is captured there.
        std::ostringstream err;
        std::streambuf * cerrBuf = std::cerr.rdbuf(err.rdbuf());
        bool threwM = false;
        try { uM.setup(bcM, dirichlet::interpolation, 0); }
        catch (const std::runtime_error &) { threwM = true; }
        std::cerr.rdbuf(cerrBuf);
        gsInfo << "[dirichletGuard] case5 threwM = "<<threwM<<", stderr = "<<err.str()<<"\n";
        CHECK(threwM);
        CHECK( err.str().find("only implemented for tensor and hierarchical bases") != std::string::npos );
    }
}

}
