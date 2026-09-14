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

}
