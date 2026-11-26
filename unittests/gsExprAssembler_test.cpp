/** @file gsExprAssembler_test.cpp

    @brief Tests for gsExprAssembler

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): R. Schneckenleitner
*/

#include "gismo_unittest.h"

#include <atomic>

namespace
{

/// Deliberately simple non-G+Smo rule used to verify that the assembler accepts
/// arbitrary gsQuadRule implementations, rather than merely another quRule
/// option.
template <class T>
class MidpointRule : public gismo::gsQuadRule<T>
{
public:
    explicit MidpointRule(index_t dim)
    {
        this->m_nodes.setZero(dim, 1);
        this->m_weights.resize(1);
        this->m_weights[0] = T(1);
        for (index_t i = 0; i < dim; ++i)
            this->m_weights[0] *= T(2);
    }
};

/// Poisson assembler whose single unknown has two components, so that the
/// legacy Dirichlet code sees an unknown size U = 2. A plain
/// gsPoissonAssembler always has U = 1; a vector-valued unknown exists only
/// through a custom gsSparseSystem, which this refresh() sets up.
class TwoComponentPoissonAssembler : public gismo::gsPoissonAssembler<real_t>
{
public:
    TwoComponentPoissonAssembler(const gismo::gsPoissonPde<real_t> & pde,
                                 const gismo::gsMultiBasis<real_t> & mb)
    : gismo::gsPoissonAssembler<real_t>(pde, mb) { }

    void refresh() override
    {
        gismo::gsDofMapper map = gismo::createMapper(m_bases[0], m_pde_ptr->bc(), 1, 0, true, true);
        std::vector<gismo::gsDofMapper> maps(2, map);
        gismo::gsVector<index_t> dims(1);
        dims << 2;
        m_system = gismo::gsSparseSystem<real_t>(maps, dims);
    }
};

} // anonymous namespace


SUITE(gsExprAssembler_test)
{
    // Verifies that a custom quadrature factory is called with the patch and
    // side context, changes the assembled mass matrix, and that clearing it
    // restores the default rule. The reference is the matrix assembled
    // without a factory in the same test; the one-point midpoint rule cannot
    // integrate the cubic mass matrix exactly, so the matrices must differ.
    TEST(CustomQuadratureFactory)
    {
        gsBSplineBasis<real_t> bb(0.0, 1.0, 3, 3);
        gsMultiBasis<real_t> mb(bb);
        gsBoundaryConditions<real_t> bcs;

        gsExprAssembler<real_t> standard(1, 1);
        standard.setIntegrationElements(mb);
        auto uStandard = standard.getSpace(mb);
        uStandard.setup(bcs, dirichlet::homogeneous, 0);
        standard.initSystem();
        standard.assemble(uStandard * uStandard.tr());
        const gsSparseMatrix<real_t> Mstandard = standard.matrix();

        gsExprAssembler<real_t> custom(1, 1);
        custom.setIntegrationElements(mb);
        auto uCustom = custom.getSpace(mb);
        uCustom.setup(bcs, dirichlet::homogeneous, 0);
        custom.initSystem();

        std::atomic<index_t> calls(0);
        std::atomic<bool> contextIsCorrect(true);
        custom.setQuadratureFactory(
            [&calls, &contextIsCorrect](const gsBasis<real_t> & basis,
                                        const gsOptionList &,
                                        index_t patch,
                                        short_t fixedDirection)
                -> gsExprAssembler<real_t>::QuadratureRulePtr
            {
                ++calls;
                if (patch != 0 || fixedDirection != -1)
                    contextIsCorrect = false;
                return gsExprAssembler<real_t>::QuadratureRulePtr(
                    new MidpointRule<real_t>(basis.dim()));
            });

        CHECK(custom.hasCustomQuadrature());
        custom.assemble(uCustom * uCustom.tr());
        const gsSparseMatrix<real_t> Mcustom = custom.matrix();
        
        const index_t current_calls = calls.load();
        CHECK(current_calls > 0);
        CHECK(contextIsCorrect.load());
        CHECK((Mstandard - Mcustom).norm() > 1e-8);

        custom.clearQuadratureFactory();
        CHECK(!custom.hasCustomQuadrature());
        custom.clearMatrix();
        custom.assemble(uCustom * uCustom.tr());
        const gsSparseMatrix<real_t> Mrestored = custom.matrix();

        CHECK_EQUAL(current_calls, calls.load());
        CHECK_EQUAL(Mstandard.rows(), Mrestored.rows());
        CHECK_EQUAL(Mstandard.cols(), Mrestored.cols());
        CHECK((Mstandard - Mrestored).norm() < 1e-14);
    }

    TEST(MultiSpaceBlockDims)
    {
        // _blockDims sizes the row/column blocks and resetDimensions sets the
        // shift applied to later blocks. mapper.freeSize() already reports the
        // component-inclusive dof count of a space, so multiplying it by that
        // space's dim counts the dimension twice. A single space of dim 1
        // cannot expose that -- the erroneous factor is 1 -- so this needs two
        // spaces of different, non-trivial dimension in one assembler.
        //
        // A 2D geometry is required: the assemble() arm below needs meas(G)
        // to build a genuine bilinear form per space, rather than only
        // inspecting mappers.
        gsMultiPatch<real_t> mp(*gsNurbsCreator<real_t>::BSplineSquare());
        gsMultiBasis<real_t> dbasis(mp);
        dbasis.degreeElevate(2);
        dbasis.uniformRefine(5);
        gsBoundaryConditions<real_t> bcs;
        const index_t n = dbasis.basis(0).size();

        // Control arm: two SCALAR spaces (dim 1 each). Since dim==1 makes
        // the erroneous factor in the fixed formula equal to 1, this arm
        // cannot itself fail on the bug -- it exists to pin that the shift
        // between blocks is otherwise sane, so a failure in the main arm
        // below can be attributed to the dim>1 handling rather than to the
        // fixture or to block bookkeeping in general.
        {
            gsExprAssembler<real_t> A(2, 2);
            A.setIntegrationElements(dbasis);
            auto u0 = A.getSpace(dbasis, 1, 0);
            auto u1 = A.getSpace(dbasis, 1, 1);
            u0.setup(bcs, dirichlet::homogeneous, 0);
            u1.setup(bcs, dirichlet::homogeneous, 0);
            A.initSystem();
            CHECK_EQUAL(2*n, A.numDofs());
            CHECK_EQUAL(n,   u1.mapper().firstIndex());
        }

        // Main arm: a vector-valued space of dimension d sharing an
        // assembler with a scalar space. Looping d = 2 and d = 3 is
        // essential -- the erroneous shift was dim()*(dim()*freeSize())
        // versus the correct dim()*freeSize(), i.e. an extra factor of
        // dim(). At a single d, a formula quadratic in dim() is
        // indistinguishable from a linear one with a different constant;
        // the pair of values pins the exponent.
        for (index_t d = 2; d <= 3; ++d)
        {
            gsExprAssembler<real_t> A(2, 2);
            A.setIntegrationElements(dbasis);
            auto G = A.getMap(mp);
            auto v = A.getSpace(dbasis, d, 0); // vector-valued space, dim d
            auto p = A.getSpace(dbasis, 1, 1); // scalar space, dim 1
            v.setup(bcs, dirichlet::homogeneous, 0);
            p.setup(bcs, dirichlet::homogeneous, 0);
            A.initSystem();

            CHECK_EQUAL((d+1)*n, A.numDofs());
            CHECK_EQUAL(0,       v.mapper().firstIndex());
            CHECK_EQUAL(d*n,     v.mapper().freeSize());
            CHECK_EQUAL(d*n,     p.mapper().firstIndex());
            CHECK_EQUAL(n,       p.mapper().freeSize());

            // Two distinct-block terms in one assemble() call: this is what
            // actually writes into both diagonal blocks of the system
            // matrix, so a wrong offset for p's block (d*(d*n) instead of
            // d*n) would either write out of range or leave part of v's
            // block untouched.
            A.assemble(v*v.tr()*meas(G), p*p.tr()*meas(G));

            // Counting entirely-zero rows separates "wrong total" from
            // "structurally broken": pre-fix, each block was internally
            // self-consistent and merely sat at the wrong offset, so
            // numDofs() alone does not reveal that (d^2-d)*n rows exist
            // that nothing ever writes to.
            const gsSparseMatrix<real_t> & M = A.matrix();
            gsVector<bool> rowTouched(M.rows());
            rowTouched.setZero();
            for (index_t c = 0; c != M.cols(); ++c)
                for (gsSparseMatrix<real_t>::InnerIterator it(M, c); it; ++it)
                    rowTouched(it.row()) = true;
            const index_t zeroRows = M.rows() - rowTouched.array().count();
            CHECK_EQUAL(0, zeroRows);

            // matrixBlockView() is the only assertion here reaching
            // _blockDims directly (numDofs()/firstIndex() above reach the
            // same defect only through resetDimensions): the block sizes
            // must match the per-space free dof counts exactly.
            auto view = A.matrixBlockView();
            CHECK_EQUAL(d*n, view(0,0).rows());
            CHECK_EQUAL(d*n, view(0,0).cols());
            CHECK_EQUAL(n,   view(1,1).rows());
            CHECK_EQUAL(n,   view(1,1).cols());
        }
    }

    // getSolution() takes the solution vector of the whole system. With two
    // spaces, each space's dofs lie at its mapper's shift in that vector.
    TEST(SolutionOfShiftedSpace)
    {
        // Set-up: the unit square, with uniformRefine(2) (two knots inserted
        // per span). Two scalar spaces are unknowns 0 and 1 of one system, and
        // there are no boundary conditions.
        gsMultiPatch<> mp(*gsNurbsCreator<>::BSplineSquare());
        gsMultiBasis<> mb(mp);
        mb.uniformRefine(2);
        gsBoundaryConditions<> bcs;

        gsExprAssembler<> A(2,2);
        A.setIntegrationElements(mb);
        auto u0 = A.getSpace(mb, 1, 0);
        auto u1 = A.getSpace(mb, 1, 1);
        u0.setup(bcs, dirichlet::homogeneous, 0);
        u1.setup(bcs, dirichlet::homogeneous, 0);
        A.initSystem();
        const index_t n = mb.basis(0).size();
        CHECK_EQUAL(2*n, A.numDofs());

        // Partition of unity: the solutions equal 1 and 2
        gsMatrix<> solVector = gsMatrix<>::Zero(A.numDofs(), 1);
        auto s0 = A.getSolution(u0, solVector);
        auto s1 = A.getSolution(u1, solVector);
        s0.setComponent(0, 1.0);
        s1.setComponent(0, 2.0);
        CHECK(s0.check());
        CHECK(s1.check());

        // A vector holding only the first space is too short for the second
        gsMatrix<> firstPart = solVector.topRows(n);
        CHECK( A.getSolution(u0, firstPart).check());
        CHECK(!A.getSolution(u1, firstPart).check());

        // eval() asserts check() inside the evaluator's parallel region
        if (s0.check() && s1.check())
        {
            // Integrals over the parameter domain (no meas) of the constants 1
            // and 2, and of the gradient of a constant. Gauss quadrature is
            // exact for them and the derivatives of the partition of unity sum
            // to 0, so only round-off remains.
            gsExprEvaluator<> ev(A);
            CHECK_CLOSE(1.0, ev.integral(s0), 1e-12);
            CHECK_CLOSE(2.0, ev.integral(s1), 1e-12);
            CHECK_CLOSE(0.0, ev.integral(grad(s1).sqNorm()), 1e-12);
        }
    }

    // matrix() is `m_modified ? makeMatrix() : m_matrix`, and m_matrix is only
    // populated from the fiber matrix by makeMatrix(). clearMatrix() must
    // therefore invalidate the cache on every path, including the one that
    // resizes rather than zeroes: initSystem() reaches exactly that path, so
    // without the flag matrix() reports the default-constructed 0x0 m_matrix
    // for a system whose dimensions are already known.
    TEST(MatrixSizedAfterInitSystem)
    {
        gsBSplineBasis<real_t> bb(0.0, 1.0, 3, 3);
        gsMultiBasis<real_t> mb(bb);
        gsBoundaryConditions<real_t> bcs;

        gsExprAssembler<real_t> A(1, 1);
        A.setIntegrationElements(mb);
        auto u = A.getSpace(mb, 1, 0);
        u.setup(bcs, dirichlet::homogeneous, 0);
        A.initSystem();

        // No assemble() here: sizing must hold from initSystem() alone.
        CHECK_EQUAL(A.numTestDofs(), A.matrix().rows());
        CHECK_EQUAL(A.numDofs(),     A.matrix().cols());
        CHECK(A.matrix().rows() > 0);
    }

    // lapl(u) of a space must be usable in an assembled form; it is the
    // parametric Laplacian, i.e. hess(u).trace()
    TEST(LaplacianOfSpace)
    {
        // Set-up: two unit squares side by side, degree elevated by 2 and
        // refined once.
        gsMultiPatch<> mp = gsNurbsCreator<>::BSplineSquareGrid(2,1,1.0);
        gsMultiBasis<> mb(mp);
        mb.degreeElevate(2);
        mb.uniformRefine(1);
        gsFunctionExpr<> ff("sin(x)*y", 2);

        gsExprAssembler<> A(1,1);
        A.setIntegrationElements(mb);
        auto G = A.getMap(mp);
        auto u = A.getSpace(mb);
        auto f = A.getCoeff(ff, G);
        u.setup(0);

        // Assembly with lapl(u): matrix M and right-hand side b.
        A.initSystem();
        A.assemble(lapl(u) * lapl(u).tr() * meas(G) + (lapl(u) + u) * u.tr() * meas(G),
                   lapl(u) * f * meas(G));
        const gsSparseMatrix<> M = A.matrix();
        const gsMatrix<>       b = A.rhs();

        // Second assembly of the same forms, with hess(u).trace() in place of
        // lapl(u): the reference.
        A.initSystem();
        A.assemble(hess(u).trace() * hess(u).trace().tr() * meas(G)
                   + (hess(u).trace() + u) * u.tr() * meas(G),
                   hess(u).trace() * f * meas(G));

        // The two lower bounds are non-vacuity guards: the relative checks
        // below would pass trivially for a zero matrix or vector. The relative
        // tolerance 1e-12 holds because the two forms are the same quantity,
        // evaluated with the same quadrature, so they agree to round-off.
        CHECK(M.norm() > 1);
        CHECK(b.norm() > 1e-3);
        CHECK((M - A.matrix()).norm() <= 1e-12 * M.norm());
        CHECK((b - A.rhs()).norm()    <= 1e-12 * b.norm());
    }

    // Verifies that re-registering a test space returns a handle to that same
    // space. The expected id 0 and the mapper identity follow from the
    // registration order: the test space of unknown 0 was registered first.
    TEST(TestSpaceReRegistrationKeepsItsData)
    {
        gsBSplineBasis<real_t> bb(0.0, 1.0, 3, 3);
        gsMultiBasis<real_t> mb(bb);

        gsExprAssembler<real_t> A(2, 2);
        auto u0 = A.getSpace(mb, 1, 0);
        A.getTestSpace(u0, mb);   // distinct test space for unknown 0
        A.getSpace(mb, 2, 1);     // registers more space data after it

        // Re-registering test space 0 updates its existing data; the
        // returned handle must refer to that data, not to whichever
        // space data happened to be registered last (unknown 1).
        auto v0 = A.getTestSpace(mb, 1, 0);
        CHECK_EQUAL(0, v0.id());
        CHECK(&v0.mapper() == &A.testSpace(0).mapper());
    }

    // For a space of dimension d, a Dirichlet condition with d-component
    // data must give component r of the space the r-th data component (also
    // when the condition constrains a single component), and scalar data
    // must be applied to every component the condition constrains. The
    // oracle is the same set of conditions split into one scalar condition
    // per component, in the same order (later conditions overwrite shared
    // corner dofs the same way).
    TEST(DirichletValuesVectorData)
    {
        // Set-up: two unit squares side by side, degree elevated by 1 and
        // refined once. g is vector data, h scalar data. bcs holds the
        // conditions as given; split holds the same conditions as scalar
        // conditions per component, the oracle.
        gsMultiPatch<> mp = gsNurbsCreator<>::BSplineSquareGrid(2,1,1.0);
        gsMultiBasis<> mb(mp);
        mb.degreeElevate(1);
        mb.uniformRefine(1);

        gsFunctionExpr<> g ("sin(1.5*x) + y", "x*y - cos(y)", 2);
        gsFunctionExpr<> g0("sin(1.5*x) + y", 2);
        gsFunctionExpr<> g1("x*y - cos(y)", 2);
        gsFunctionExpr<> h ("x + 2*y*y", 2);

        gsBoundaryConditions<> bcs, split;
        bcs.setGeoMap(mp);
        split.setGeoMap(mp);
        // vector data, all components
        bcs  .addCondition(0, boundary::west , condition_type::dirichlet, &g , 0, false, -1);
        split.addCondition(0, boundary::west , condition_type::dirichlet, &g0, 0, false,  0);
        split.addCondition(0, boundary::west , condition_type::dirichlet, &g1, 0, false,  1);
        // vector data, component 1 only
        bcs  .addCondition(1, boundary::south, condition_type::dirichlet, &g , 0, false,  1);
        split.addCondition(1, boundary::south, condition_type::dirichlet, &g1, 0, false,  1);
        // scalar data, all components
        bcs  .addCondition(1, boundary::east , condition_type::dirichlet, &h , 0, false, -1);
        split.addCondition(1, boundary::east , condition_type::dirichlet, &h , 0, false,  0);
        split.addCondition(1, boundary::east , condition_type::dirichlet, &h , 0, false,  1);
        // scalar data, component 0 only
        bcs  .addCondition(0, boundary::north, condition_type::dirichlet, &h , 0, false,  0);
        split.addCondition(0, boundary::north, condition_type::dirichlet, &h , 0, false,  0);

        // Both Dirichlet methods are checked: interpolation and L2 projection.
        const index_t methods[2] = {dirichlet::interpolation, dirichlet::l2Projection};
        for (index_t m = 0; m != 2; ++m)
        {
            // u carries the vector-data conditions, w the scalar split of the
            // same conditions, both with the same Dirichlet method.
            gsExprAssembler<> A(1,1), B(1,1);
            A.setIntegrationElements(mb);
            B.setIntegrationElements(mb);
            auto u = A.getSpace(mb, 2);
            auto w = B.getSpace(mb, 2);
            u.setup(bcs,   methods[m], 0);
            w.setup(split, methods[m], 0);
            CHECK_EQUAL(w.mapper().boundarySize(), u.mapper().boundarySize());

            // Compare every constrained dof of every component, and track the
            // largest difference and the largest value.
            real_t maxDiff = 0, maxVal = 0;
            index_t numFixed = 0;
            for (size_t k = 0; k != mp.nPatches(); ++k)
                for (index_t i = 0; i != mb.basis(k).size(); ++i)
                    for (index_t r = 0; r != 2; ++r)
                    {
                        CHECK_EQUAL(w.mapper().is_boundary(i,k,r), u.mapper().is_boundary(i,k,r));
                        if (!u.mapper().is_boundary(i,k,r)) continue;
                        const real_t a = u.fixedPart().at(u.mapper().bindex(i,k,r));
                        const real_t b = w.fixedPart().at(w.mapper().bindex(i,k,r));
                        maxDiff = math::max(maxDiff, math::abs(a-b));
                        maxVal  = math::max(maxVal , math::abs(b));
                        ++numFixed;
                    }
            // The data are O(1): on the east side of patch 1 (x = 2), h = x +
            // 2y^2 >= 2, so the relative check below is not against near-zero
            // values. Both set-ups compute the fixed dofs from the same data with
            // the same method (interpolation at the same anchors, or the same
            // boundary L2 system solved to machine-epsilon tolerance), so they
            // differ by round-off. 1e-10 (about 4.5e5 eps) leaves room for the
            // conditioning of the small boundary systems.
            CHECK(numFixed > 0);
            CHECK(maxVal > 0.5);
            CHECK(maxDiff <= 1e-10 * maxVal);
        }
    }

    // The side of a curve is its end point, where the L2 projection of
    // Dirichlet data reduces to the value there: the dof of the end basis
    // function takes the data value, whatever the tangent of the curve.
    // G(t) = t^2, on the line and as a curve in the plane, has a vanishing
    // tangent at t = 0; x = 2t is a regular curve. The expression assembler
    // projects for all three curves; gsAssembler (through gsPoissonAssembler)
    // only for the two curves on the line, because the loop skips the planar
    // curve with `if (1 != dim) continue;`.
    TEST(DirichletL2ProjectionAtCurveEnd)
    {
        // Constant Dirichlet data at the west and east end of each curve. The
        // three degree-2 curves, on one span, are x = t^2, the planar
        // (t^2, 2t^2), and x = 2t.
        const real_t gW = 0.75, gE = -1.25;
        gsKnotVector<> kv(0, 1, 0, 3); // degree 2, one span
        gsMatrix<> cps[3];
        cps[0].resize(3, 1); cps[0] << 0, 0, 1;           // x = t^2
        cps[1].resize(3, 2); cps[1] << 0, 0,  0, 0,  1, 2; // (t^2, 2t^2)
        cps[2].resize(3, 1); cps[2] << 0, 1, 2;           // x = 2t
        for (index_t c = 0; c != 3; ++c)
        {
            gsMultiPatch<> mp;
            mp.addPatch(gsBSpline<>(kv, cps[c]));
            gsMultiBasis<> mb(mp);
            // uniformRefine(2) inserts two knots per span
            mb.uniformRefine(2);
            const index_t n = mb.basis(0).size();
            const short_t dim = mp.geoDim();

            gsConstantFunction<> fW(gW, dim), fE(gE, dim);
            gsBoundaryConditions<> bc;
            bc.addCondition(0, boundary::west, condition_type::dirichlet, &fW);
            bc.addCondition(0, boundary::east, condition_type::dirichlet, &fE);
            bc.setGeoMap(mp);

            gsExprAssembler<> A(1,1);
            A.setIntegrationElements(mb);
            auto u = A.getSpace(mb);
            u.setup(bc, dirichlet::l2Projection, 0);
            CHECK_EQUAL(2, u.mapper().boundarySize());
            // On an end point each boundary dof has a 1x1 equation
            // N(end)^2 x = N(end) g, with N(end) = 1 and measure 1. The CG solver
            // solves this diagonal system to round-off, and 1e-10 is a round-off
            // bound.
            CHECK_CLOSE(gW, u.fixedPart().at(u.mapper().bindex(0  , 0)), 1e-10);
            CHECK_CLOSE(gE, u.fixedPart().at(u.mapper().bindex(n-1, 0)), 1e-10);

            // The same projection through gsAssembler, here gsPoissonAssembler,
            // for the curves on the line only. The tolerance is that of the
            // expression assembler above.
            if (1 != dim) continue;
            gsConstantFunction<> f(0, 1);
            gsPoissonAssembler<> P(mp, mb, bc, f);
            P.options().setInt("DirichletValues", dirichlet::l2Projection);
            P.computeDirichletDofs();
            const gsDofMapper & m = P.system().colMapper(0);
            CHECK_EQUAL(2, m.boundarySize());
            CHECK_CLOSE(gW, P.fixedDofs().at(m.bindex(0  , 0)), 1e-10);
            CHECK_CLOSE(gE, P.fixedDofs().at(m.bindex(n-1, 0)), 1e-10);
        }
    }

    // Dirichlet data with fewer components than the condition addresses is
    // refused when the space is set up (a runtime error, in every build type),
    // for both ways of computing the values. A 3-component space with
    // 2-component data (x, y) is refused for a condition on all components
    // (com = -1 needs 3) and on component 2 (needs 3), and accepted for
    // component 1 (needs 2), which reads component 1 of the data, y. The
    // accepted case also shows the data is used: the boundary values hold y,
    // which reaches 1 at the top corner of the unit square, so a maximum above
    // 0.5 (half of the exact maximum 1) rules out values left at zero.
    TEST(DirichletValuesTooFewComponents)
    {
        gsMultiPatch<> mp = gsNurbsCreator<>::BSplineSquareGrid(2,1,1.0);
        gsMultiBasis<> mb(mp);
        mb.degreeElevate(1);
        mb.uniformRefine(1);
        gsFunctionExpr<> g("x", "y", 2);

        const index_t methods[2] = {dirichlet::interpolation, dirichlet::l2Projection};
        const index_t comps[3]   = {-1, 2, 1};
        for (index_t m = 0; m != 2; ++m)
            for (index_t c = 0; c != 3; ++c)
            {
                gsBoundaryConditions<> bcs;
                bcs.setGeoMap(mp);
                bcs.addCondition(0, boundary::west, condition_type::dirichlet, &g, 0, false, comps[c]);

                gsExprAssembler<> A(1,1);
                A.setIntegrationElements(mb);
                auto u = A.getSpace(mb, 3);
                if (c < 2)
                    CHECK_THROW(u.setup(bcs, methods[m], 0), std::runtime_error);
                else
                {
                    u.setup(bcs, methods[m], 0);
                    CHECK(u.mapper().boundarySize() > 0);
                    CHECK(u.fixedPart().cwiseAbs().maxCoeff() > 0.5);
                }
            }
    }

    // The legacy gsAssembler accepts Dirichlet data of one component for an
    // unknown of U components with a single right-hand side (the value is used
    // for every component), and refuses any other mismatch of the number of
    // data components with U times the number of right-hand sides, in both the
    // interpolation and the L2 projection. Oracle for the accepted case:
    // gsDirichletValues (through gsFeSpace::setup) on the same condition and a
    // 2-component space, which also broadcasts scalar data. Both compute the
    // same thing: the interpolation uses the anchors of the side's boundary
    // basis, the L2 projection the Gauss rule of that basis with 3 points
    // per direction (degree 2 after one elevation, exact up to degree 5), so
    // the cubic datum 1 + x + y^3 is integrated exactly by both, and the two
    // CG solves of a few unknowns agree to far below the tolerance
    // 1e-10 * maxVal, a relative tolerance on values of order 1. On the west
    // side x = 0, so the datum is 1 + y^3, between 1 and 2; maxVal > 0.5
    // (half the datum's minimum) rules out values left at zero.
    TEST(LegacyDirichletScalarDataOnVectorUnknown)
    {
        gsMultiPatch<> mp = gsNurbsCreator<>::BSplineSquareGrid(1,1,1.0);
        gsMultiBasis<> mb(mp);
        mb.degreeElevate(1);

        gsFunctionExpr<> h ("1 + x + y*y*y", 2);
        gsFunctionExpr<> h3("x", "y", "x*y", 2);
        gsFunctionExpr<> rhs1("1", 2);
        gsFunctionExpr<> rhs2("1", "2", 2);

        const index_t methods[2] = {dirichlet::interpolation, dirichlet::l2Projection};
        for (index_t method : methods)
        {
            gsBoundaryConditions<> bc;
            bc.setGeoMap(mp);
            bc.addCondition(0, boundary::west, condition_type::dirichlet, &h);

            // accepted: scalar data, U = 2, one right-hand side
            gsPoissonPde<> pde(mp, bc, rhs1);
            TwoComponentPoissonAssembler A(pde, mb);
            A.options().setInt("DirichletValues", method);
            A.refresh();
            A.computeDirichletDofs(0);

            gsExprAssembler<> E(1,1);
            E.setIntegrationElements(mb);
            auto u = E.getSpace(mb, 2);
            u.setup(bc, method, 0);

            const gsDofMapper & m = A.system().colMapper(0);
            real_t maxDiff = 0, maxVal = 0;
            index_t numFixed = 0;
            for (index_t i = 0; i != mb.basis(0).size(); ++i)
            {
                if (!m.is_boundary(i, 0)) continue;
                for (index_t r = 0; r != 2; ++r)
                {
                    CHECK(u.mapper().is_boundary(i, 0, r));
                    const real_t a = A.fixedDofs()(m.bindex(i, 0), r);
                    const real_t b = u.fixedPart().at(u.mapper().bindex(i, 0, r));
                    maxDiff = math::max(maxDiff, math::abs(a-b));
                    maxVal  = math::max(maxVal , math::abs(b));
                    ++numFixed;
                }
            }
            CHECK(numFixed > 0);
            CHECK(maxVal > 0.5);
            CHECK(maxDiff <= 1e-10 * maxVal);

            // refused: 3 components of data for U = 2
            gsBoundaryConditions<> bc3;
            bc3.setGeoMap(mp);
            bc3.addCondition(0, boundary::west, condition_type::dirichlet, &h3);
            gsPoissonPde<> pde3(mp, bc3, rhs1);
            TwoComponentPoissonAssembler A3(pde3, mb);
            A3.options().setInt("DirichletValues", method);
            A3.refresh();
            CHECK_THROW(A3.computeDirichletDofs(0), std::runtime_error);

            // refused: scalar data with two right-hand sides, 1 != U * R = 4
            gsPoissonPde<> pde2(mp, bc, rhs2);
            TwoComponentPoissonAssembler A2(pde2, mb);
            A2.options().setInt("DirichletValues", method);
            A2.refresh();
            CHECK_THROW(A2.computeDirichletDofs(0), std::runtime_error);
        }
    }

    // Verifies the domain and interface integrals of a function with a jump:
    // f = sign(x) is odd in x and both the domain and the interface are
    // symmetric about x = 0 (the grid is translated by -0.5), so the exact
    // integrals of f and of f.left() + f.right() are 0. The tests compare the
    // squared value with 1e-10.
    TEST(InterfaceExpression)
    {
        const index_t numRef = 2;
        gsVector<> translation(2);
        translation << -0.5, -1;
        gsFunctionExpr<> ff("if(x>0,1,-1)", 2); // function with jump
        gsMultiPatch<> patches = gsNurbsCreator<>::BSplineSquareGrid(1,2,1);
        patches.patch(0).translate(translation);
        patches.patch(1).translate(translation);

        gsMultiBasis<> mb(patches);

        for(index_t i = 0; i < numRef; i++)
            mb.uniformRefine();

        gsExprEvaluator<> ev;
        ev.setIntegrationDomain(mb.domain());
        gsExprEvaluator<>::geometryMap G = ev.getMap(patches);
        auto f = ev.getVariable(ff, G);

        ev.integral(f);
        const real_t v = ev.value();
        CHECK( v*v < 1e-10 );

        CHECK(1==patches.interfaces().size());
        ev.integralInterface(f.left() + f.right() , patches.interfaces());
        const real_t w = ev.value();
        CHECK( w*w < 1e-10 );
    }

    // A solution taken with right() must be evaluated on the second patch of
    // the interface, and left() on the first one
    TEST(InterfaceSolution)
    {
        gsMultiPatch<> patches = gsNurbsCreator<>::BSplineSquareGrid(1,2,1);
        gsMultiBasis<> mb(patches);
        mb.uniformRefine();

        gsExprAssembler<> A(1,1);
        A.setIntegrationElements(mb);
        auto u = A.getSpace(mb);
        u.setup(-1); // no coupling: the solution may jump across the interface

        // Partition of unity: the solution equals 1+p on patch p
        gsMatrix<> solVector = gsMatrix<>::Zero(u.mapper().freeSize(), 1);
        auto s = A.getSolution(u, solVector);
        s.setComponent(0, 1.0, 0);
        s.setComponent(0, 2.0, 1);

        CHECK(1==patches.interfaces().size());
        const boundaryInterface & iFace = patches.interfaces().front();
        const real_t first  = 1.0 + iFace.first ().patch;
        const real_t second = 1.0 + iFace.second().patch;

        CHECK(!s.isAcross());
        CHECK( s.right().isAcross());
        CHECK(!s.right().left().isAcross());

        // The interface has parametric length one
        gsExprEvaluator<> ev(A);
        CHECK_CLOSE(first , ev.integralInterface(s,                 patches.interfaces()), 1e-12);
        CHECK_CLOSE(first , ev.integralInterface(s.left(),          patches.interfaces()), 1e-12);
        CHECK_CLOSE(second, ev.integralInterface(s.right(),         patches.interfaces()), 1e-12);
        CHECK_CLOSE(first , ev.integralInterface(s.right().left(),  patches.interfaces()), 1e-12);
        CHECK_CLOSE(second - first,
                    ev.integralInterface(s.right() - s.left(), patches.interfaces()), 1e-12);
    }

    // assembleJacobianIfc() differentiates an interface residual by finite
    // differences. For a residual linear in the solution, the Jacobian is the
    // matrix assembleIfc() builds from the corresponding bilinear form.
    TEST(InterfaceJacobian)
    {
        // Set-up: two stacked unit squares, degree 2, uniformRefine(2) (two
        // knots per span), unglued (setup(-1)), and a random solution vector.
        gsMultiPatch<> patches = gsNurbsCreator<>::BSplineSquareGrid(1,2,1);
        gsMultiBasis<> mb(patches);
        mb.degreeElevate(1);
        mb.uniformRefine(2);

        gsExprAssembler<> A(1,1);
        A.setIntegrationElements(mb);
        auto G = A.getMap(patches);
        auto u = A.getSpace(mb);
        u.setup(-1);
        A.initSystem();

        gsMatrix<> solVector = gsMatrix<>::Random(A.numDofs(), 1);
        auto s = A.getSolution(u, solVector);

        // J: the finite-difference Jacobian of the residual u * s.
        A.assembleJacobianIfc(patches.interfaces(), u * s * nv(G).norm(), s);
        const gsSparseMatrix<> J = A.matrix();

        // M: assembleIfc() of the bilinear form u * u.tr() the residual
        // linearizes to.
        A.initSystem();
        A.assembleIfc(patches.interfaces(), u * u.tr() * nv(G).norm());
        const gsSparseMatrix<> M = A.matrix();

        // The tolerance 1e-8 is that of the finite-difference Jacobian with the
        // default fixed interface (movingInterface == false): the five-point
        // central-difference stencil with step 1e-5 is exact for a residual
        // that is a polynomial of degree at most 4 in the perturbed dof, and
        // this one is linear, so only cancellation round-off remains, about
        // eps/step, 2e-11 relative. 1e-8 leaves about three orders of magnitude
        // of margin.
        CHECK(M.norm() > 0);
        CHECK((J - M).norm() <= 1e-8 * M.norm());
    }

    // The interface Jacobian also differentiates with respect to the
    // solution on the second side (s.right()), and has the rows of the test
    // functions of either side: compared with assembleIfc() of the
    // linearized forms, for unglued patches and for a C0 interface whose
    // dofs both sides share.
    TEST(InterfaceJacobianBothSides)
    {
        gsMultiPatch<> patches = gsNurbsCreator<>::BSplineSquareGrid(1,2,1);
        gsMultiBasis<> mb(patches);
        mb.degreeElevate(1);
        mb.uniformRefine(2);

        for (index_t glued = 0; glued != 2; ++glued)
        {
            gsExprAssembler<> A(1,1);
            A.setIntegrationElements(mb);
            auto G = A.getMap(patches);
            auto u = A.getSpace(mb);
            if (glued)
                u.setup(gsBoundaryConditions<>(), dirichlet::homogeneous, 0);
            else
                u.setup(-1);
            A.initSystem();

            gsMatrix<> solVector = gsMatrix<>::Random(A.numDofs(), 1);
            auto s = A.getSolution(u, solVector);

            // J of the residual minus M of the linearized form, relative to M.
            // The tolerance 1e-8 is that of the finite-difference Jacobian with
            // the default fixed interface (movingInterface == false): the
            // five-point central-difference stencil with step 1e-5 is exact for
            // a residual that is a polynomial of degree at most 4 in the
            // perturbed dof, so only cancellation round-off remains, about
            // eps/step, 2e-11 relative. The quadratic case u * s_R^2 is still
            // within the exact range of the stencil.
            gsSparseMatrix<> J, M;
            auto relErr = [&]() { return (J - M).norm() / M.norm(); };

            // u * s.right() against u * u.right().tr(): the dependence on the
            // dofs of the second side.
            A.initSystem();
            A.assembleJacobianIfc(patches.interfaces(), u * s.right() * nv(G).norm(), s);
            J = A.matrix();
            A.initSystem();
            A.assembleIfc(patches.interfaces(), u * u.right().tr() * nv(G).norm());
            M = A.matrix();
            CHECK(relErr() <= 1e-8);

            // u.right() * s against u.right() * u.tr(): the test functions of
            // the second side.
            A.initSystem();
            A.assembleJacobianIfc(patches.interfaces(), u.right() * s * nv(G).norm(), s);
            J = A.matrix();
            A.initSystem();
            A.assembleIfc(patches.interfaces(), u.right() * u.tr() * nv(G).norm());
            M = A.matrix();
            CHECK(relErr() <= 1e-8);

            // nonlinear: d/ds (u s_R^2) = 2 s_R u u_R^T
            A.initSystem();
            A.assembleJacobianIfc(patches.interfaces(),
                                  u * (s.right() * s.right()) * nv(G).norm(), s);
            J = A.matrix();
            A.initSystem();
            A.assembleIfc(patches.interfaces(),
                          2 * u * s.right() * u.right().tr() * nv(G).norm());
            M = A.matrix();
            CHECK(relErr() <= 1e-8);

            // the normal derivative from the second side also depends on its
            // dofs off the interface, glued or not
            A.initSystem();
            A.assembleJacobianIfc(patches.interfaces(),
                                  u * (igrad(s.right(), G.right()) * nv(G)), s);
            J = A.matrix();
            A.initSystem();
            A.assembleIfc(patches.interfaces(),
                          u * (igrad(u.right(), G.right()) * nv(G)).tr());
            M = A.matrix();
            CHECK(relErr() <= 1e-8);
        }
    }

    // Whether assembleJacobianIfc() perturbs the dofs of the second side
    // depends on the residual only, not on the interface operations the
    // assembler did before. For a residual on the first side, perturbing the
    // second side's dofs leaves finite-difference round-off at structurally
    // zero positions, so the sparsity would change between calls.
    TEST(InterfaceJacobianRepeatable)
    {
        gsMultiPatch<> patches = gsNurbsCreator<>::BSplineSquareGrid(1,2,1);
        gsMultiBasis<> mb(patches);
        mb.degreeElevate(1);
        mb.uniformRefine(2);

        gsMatrix<> solVector;
        // The Jacobian of u * s from the last of \a calls calls on a new
        // assembler, after an assembleIfc() if \a ifcFirst. Each call starts
        // from the same solution vector, as the perturbations leave
        // round-off in it.
        const auto jacobian = [&](index_t calls, bool ifcFirst)
        {
            gsExprAssembler<> A(1,1);
            A.setIntegrationElements(mb);
            auto G = A.getMap(patches);
            auto u = A.getSpace(mb);
            u.setup(-1);
            A.initSystem();
            if (0 == solVector.size())
                solVector = gsMatrix<>::Random(A.numDofs(), 1);
            gsMatrix<> sol = solVector;
            auto s = A.getSolution(u, sol);

            if (ifcFirst)
                A.assembleIfc(patches.interfaces(), u * u.tr() * nv(G).norm());
            for (index_t k = 0; k != calls; ++k)
            {
                sol = solVector;
                A.initSystem();
                A.assembleJacobianIfc(patches.interfaces(), u * s * nv(G).norm(), s);
            }
            return gsSparseMatrix<>(A.matrix());
        };

        // J1 is the first call. A second call on the same assembler, and a call
        // after an assembleIfc(), must give the sparsity and the values of the
        // first call. Every call restarts from the same vector (sol =
        // solVector) and runs the same finite differences on the same elements,
        // so the Jacobians differ at most by round-off: 1e-12 relative is a
        // round-off bound, far tighter than the 1e-8 with which a
        // finite-difference Jacobian is compared with an exact matrix.
        const gsSparseMatrix<> J1 = jacobian(1, false);
        CHECK(J1.norm() > 0);
        const gsSparseMatrix<> later[2] = {jacobian(2, false), jacobian(1, true)};
        for (const gsSparseMatrix<> & J : later)
        {
            CHECK_EQUAL(J1.nonZeros(), J.nonZeros());
            CHECK((J - J1).norm() <= 1e-12 * J1.norm());
        }
    }

    // On a non-matching interface, the quadrature points of one element of
    // the first side lie in several elements of the second side, so the
    // functions active there differ between the points. The interface
    // Jacobian has the union of them as columns (for the dofs of the second
    // side) and as rows (for test functions on the second side): compared
    // with assembleIfc() of the linearized forms, for scalar and vector
    // valued spaces.
    TEST(InterfaceJacobianNonMatching)
    {
        gsMultiPatch<> patches = gsNurbsCreator<>::BSplineSquareGrid(1,2,1);
        CHECK(1==patches.interfaces().size());
        const boundaryInterface & iFace = patches.interfaces().front();
        gsMultiBasis<> mb(patches);
        mb.degreeElevate(1);
        mb.uniformRefine();
        // refine the second side along the interface only
        const short_t along = 1 - iFace.second().side().direction();
        mb.basis(iFace.second().patch).uniformRefine(1, 1, along);
        CHECK(mb.basis(iFace.second().patch).component(along).numElements() >
              mb.basis(iFace.first ().patch).component(along).numElements());

        // The number of columns of \a m above \a tol; finite differences
        // leave round-off in the columns of dofs the residual does not
        // depend on
        const auto usedCols = [](const gsSparseMatrix<> & m, real_t tol)
        {
            index_t n = 0;
            for (index_t j = 0; j != m.cols(); ++j)
                n += (m.col(j).norm() > tol);
            return n;
        };

        for (index_t dim = 1; dim != 3; ++dim)
        {
            gsExprAssembler<> A(1,1);
            A.setIntegrationElements(mb);
            auto G = A.getMap(patches);
            auto u = A.getSpace(mb, dim);
            u.setup(-1);
            A.initSystem();

            gsMatrix<> solVector = gsMatrix<>::Random(A.numDofs(), 1);
            auto s = A.getSolution(u, solVector);

            // J of the residual minus M of the linearized form, relative to M
            gsSparseMatrix<> J, M;
            // The tolerance 1e-8 is that of the finite-difference Jacobian with
            // the default fixed interface (movingInterface == false), which
            // uses the five-point central-difference stencil with step 1e-5. The
            // stencil is exact for a residual that is a polynomial of degree at
            // most 4 in the perturbed dof, and every residual here is linear or
            // quadratic, so there is no truncation error. What remains is
            // cancellation round-off, about eps/step, 2e-11 relative, so 1e-8
            // leaves about three orders of magnitude of margin. As a column
            // threshold, 1e-8 * |M| separates the columns the residual depends
            // on (entries of the size of M's) from columns holding only
            // finite-difference round-off.
            const auto check = [&]()
            {
                CHECK(M.norm() > 0);
                CHECK((J - M).norm() <= 1e-8 * M.norm());
                CHECK_EQUAL(usedCols(M, 1e-8 * M.norm()), usedCols(J, 1e-8 * M.norm()));
            };

            // u * s.right(): the columns come from the second side.
            A.initSystem();
            A.assembleJacobianIfc(patches.interfaces(), u * s.right() * nv(G).norm(), s);
            J = A.matrix();
            A.initSystem();
            A.assembleIfc(patches.interfaces(), u * u.right().tr() * nv(G).norm());
            M = A.matrix();
            check();

            // u.right() * s: the rows are those of second-side test functions.
            A.initSystem();
            A.assembleJacobianIfc(patches.interfaces(), u.right() * s * nv(G).norm(), s);
            J = A.matrix();
            A.initSystem();
            A.assembleIfc(patches.interfaces(), u.right() * u.tr() * nv(G).norm());
            M = A.matrix();
            check();

            // u.right() * (s - s.right()): both, with its two-form assembleIfc
            // as the reference.
            A.initSystem();
            A.assembleJacobianIfc(patches.interfaces(),
                                  u.right() * (s - s.right()) * nv(G).norm(), s);
            J = A.matrix();
            A.initSystem();
            A.assembleIfc(patches.interfaces(),
                          u.right() * u.tr() * nv(G).norm(),
                          (-1) * u.right() * u.right().tr() * nv(G).norm());
            M = A.matrix();
            check();

            if (1 != dim) continue;
            // the normal derivative from the second side, which depends on
            // its dofs off the interface too
            A.initSystem();
            A.assembleJacobianIfc(patches.interfaces(),
                                  u.right() * (igrad(s.right(), G.right()) * nv(G)), s);
            J = A.matrix();
            A.initSystem();
            A.assembleIfc(patches.interfaces(),
                          u.right() * (igrad(u.right(), G.right()) * nv(G)).tr());
            M = A.matrix();
            check();
        }
    }

    // assembleJacobianIfc() must perturb the dofs of the second side only
    // when the residual depends on them. A residual that reaches across the
    // interface only through the geometry (G.right()) and not through the
    // solution space is linear in the dofs of u on the first side. Then:
    //  - J equals the matrix assembleIfc() builds from the linearized form,
    //    u * u.tr() with the same geometry weight; the finite-difference
    //    stencil is exact for a residual linear in the dofs, so only
    //    round-off, of order epsilon * |R| / delta (delta = 1e-5, about 1e-11
    //    relative in double precision), remains, below the tolerance 1e-8
    //    that the sibling tests use;
    //  - for unglued patches (setup(-1)) the dofs of the second patch belong
    //    to the second side only, so their columns are structurally zero and
    //    nothing may be stored in them. Perturbing them would store the
    //    round-off of the stencil (8a - a - 8a + a)/(12 delta), which does not
    //    cancel exactly in floating point, as stored entries. The count of
    //    stored entries is exact (0), because a perturbation that was
    //    skipped leaves no entry at all.
    TEST(InterfaceJacobianGeometryAcross)
    {
        gsMultiPatch<> patches = gsNurbsCreator<>::BSplineSquareGrid(1,2,1);
        gsMultiBasis<> mb(patches);
        mb.degreeElevate(1);
        mb.uniformRefine(2);

        gsExprAssembler<> A(1,1);
        A.setIntegrationElements(mb);
        auto G = A.getMap(patches);
        auto u = A.getSpace(mb);
        u.setup(-1);
        A.initSystem();

        gsMatrix<> solVector = gsMatrix<>::Random(A.numDofs(), 1);
        auto s = A.getSolution(u, solVector);

        A.assembleJacobianIfc(patches.interfaces(),
                              u * s * nv(G).norm() * (G.right().tr()*G.right()).val(), s);
        const gsSparseMatrix<> J = A.matrix();

        A.initSystem();
        A.assembleIfc(patches.interfaces(),
                      u * u.tr() * nv(G).norm() * (G.right().tr()*G.right()).val());
        const gsSparseMatrix<> M = A.matrix();

        CHECK(M.norm() > 0);
        CHECK((J - M).norm() <= 1e-8 * M.norm());

        // stored entries in the columns of the second patch
        const index_t p2 = patches.interfaces().front().second().patch;
        index_t numCols = 0, numStored = 0;
        for (index_t i = 0; i != mb.basis(p2).size(); ++i)
        {
            ++numCols;
            for (gsSparseMatrix<>::InnerIterator it(J, u.mapper().index(i, p2)); it; ++it)
                ++numStored;
        }
        CHECK(numCols > 0);
        CHECK_EQUAL(0, numStored);
    }

    // A mapper installed by setupMapper() must be the one initSystem() uses.
    // initSystem() rebuilds any mapper whose component count differs from the
    // space dimension, discarding its eliminated dofs, so setupMapper() has to
    // reject such a mapper instead of accepting it.
    TEST(SetupMapperComponents)
    {
        gsMultiPatch<> patches = gsNurbsCreator<>::BSplineSquareGrid(1,1,1);
        gsMultiBasis<> mb(patches);
        mb.uniformRefine();

        const auto boundaryEliminated = [&mb](index_t nComp)
        {
            gsDofMapper m = createMapper(mb, nComp);
            for (index_t c = 0; c != nComp; ++c)
                m.markBoundary(0, mb.basis(0).allBoundary(), c);
            m.finalize();
            return m;
        };

        gsExprAssembler<> A(1, 1);
        A.setIntegrationElements(mb);
        auto u = A.getSpace(mb, 2);

        CHECK_THROW(u.setupMapper(boundaryEliminated(3)), std::runtime_error);
        CHECK_THROW(u.setupMapper(boundaryEliminated(1)), std::runtime_error);

        const gsDofMapper matching = boundaryEliminated(2);
        u.setupMapper(matching);
        A.initSystem();
        CHECK_EQUAL(matching.freeSize(), u.mapper().freeSize());
        CHECK_EQUAL(matching.freeSize(), A.numDofs());
    }

    // Verifies boundary, area and element-area integrals on the unit disk
    // (curve and surface in 3D, then the planar disk). The references are the
    // analytic circumference 2*pi and area pi, and for el.area(G) the integral
    // over the parameter domain of each element's area, i.e. the disk area pi
    // times the parametric element size 1/16 (16 elements after two
    // refinements) = 2*pi/32.
    TEST(BoundaryIntegral)
    {
        // Create a circle
        gsMultiPatch<> mp;
        mp.addPatch(gsNurbsCreator<>::NurbsDisk(1.0));
        mp.computeTopology();
        mp.embed(3);
        mp.uniformRefine(1);
        mp.uniformRefine(1);

        // Rotate it 30 degrees
        gsVector<real_t,3> rv;
        rv.setZero(); rv[0]=1;
        real_t angle = EIGEN_PI/3;
        mp.patch(0).rotate(angle, rv);

        // Create evaluator
        gsMultiBasis<> mb(mp);
        gsExprEvaluator<> ev;
        ev.setIntegrationDomain(mb.domain());
        ev.options().addReal("quA","",2); // added precision to approx pi
        ev.options().addInt("quB","", 2); // added precision to approx pi
        auto G = ev.getMap(mp);
        typedef gsExprAssembler<>::element element;
        element el = ev.getElement();

        // Test measure of a curve in 3D
        CHECK_CLOSE(ev.integralBdr(meas(G)), 2*EIGEN_PI, EPSILON);
        // Test tangent of a curve in 3D
        CHECK(math::abs(ev.integralBdr(tv(G).norm())-2*EIGEN_PI) < 1e-10);

        // Test measure of a surface in 3D
        CHECK(math::abs(ev.integral(meas(G))-EIGEN_PI) < 1e-10);
        // Test surface normal of a surface in 3D
        CHECK(math::abs(ev.integral(sn(G).norm())-EIGEN_PI) < 1e-10);
        // XXXX
        CHECK(math::abs(ev.integral(el.area(G))-2*EIGEN_PI/32) < 1e-10);

        mp.patch(0).rotate(-angle, rv);
        mp.embed(2);

        // Test measure of a curve in 3D
        CHECK(math::abs(ev.integralBdr(meas(G))-2*EIGEN_PI) < 1e-10);
        // Test tangent of a curve in 3D
        CHECK(math::abs(ev.integralBdr(tv(G).norm())-2*EIGEN_PI) < 1e-10);
        // Test outer normal of a curve in 3D
        CHECK(math::abs(ev.integralBdr(nv(G).norm())-2*EIGEN_PI) < 1e-10);
        //
        CHECK(math::abs(ev.integral(el.area(G))-2*EIGEN_PI/32) < 1e-10);
    }

    // extractFull() is component-major: rows are ordered by component, then
    // patch offset, then patch-local index. The expected entries are those of
    // the per-patch extract() result, which does not depend on the full ordering.
    TEST(ExtractFullComponentMajor)
    {
        // Two patches, one above the other, with a Dirichlet condition on the west side of patch 0, so the
        // mapper has coupled and eliminated dofs; a vector-valued space (dimension 2) makes the
        // component-major order observable
        gsMultiPatch<> patches = gsNurbsCreator<>::BSplineSquareGrid(1,2,1);
        gsMultiBasis<> mb(patches);
        mb.uniformRefine();
        gsBoundaryConditions<> bcs;
        bcs.addCondition(0, boundary::west, condition_type::dirichlet,
                         static_cast<gsFunctionSet<real_t>*>(nullptr));

        gsExprAssembler<> A(1,1);
        A.setIntegrationElements(mb);
        auto u = A.getSpace(mb, 2);
        u.setup(bcs, dirichlet::homogeneous, 0);
        const gsDofMapper & dm = u.mapper();

        // Free coefficient j is 1+j and eliminated coefficient j is -1-j: all values are distinct and
        // free and fixed ones have opposite signs, so a misplaced entry cannot match by accident
        gsMatrix<> sol(dm.freeSize(), 1);
        for (index_t j = 0; j < sol.rows(); ++j) sol(j) = 1 + j;
        for (index_t j = 0; j < (index_t)dm.boundarySize(); ++j) u.fixedPart()(j) = -1 - j;

        auto s = A.getSolution(u, sol);
        gsMatrix<> full;
        s.extractFull(full);

        // The fixture has two patches with both eliminated and coupled dofs, and the full vector has
        // one row per mapper entry (patch-local dofs of all patches and components)
        CHECK_EQUAL(2, (index_t)dm.numPatches());
        CHECK(dm.boundarySize() > 0);
        CHECK(dm.coupledSize() > 0);
        CHECK_EQUAL((index_t)dm.mapSize(), full.rows());

        // For patch p and component c, row (component block + patch offset + local index) of the full
        // vector equals column c of the patch-local extraction. The block of component c starts after
        // the totalSize() rows of components 0..c-1. Both entries are copies of the same stored
        // coefficient (free solution entry or Dirichlet value), with no arithmetic applied, so they
        // are compared exactly and not within a tolerance.
        for (index_t p = 0; p < 2; ++p)
        {
            gsMatrix<> cf;
            s.extract(cf, p);
            for (index_t c = 0; c < 2; ++c)
            {
                index_t base = 0;
                for (index_t q = 0; q < c; ++q) base += (index_t)dm.totalSize(q);
                for (index_t i = 0; i < (index_t)dm.patchSize(p,c); ++i)
                {
                    const index_t row = base + (index_t)dm.offset(p) + i;
                    CHECK(row < full.rows());
                    if (row < full.rows())
                        CHECK_EQUAL(cf(i,c), full(row,0));
                }
            }
        }
    }

    // getSolution() takes a vector with the structure of the whole system,
    // so check() must accept it for every space of a multi-space assembler.
    // A vector with too few rows must be rejected, also when it has several
    // columns so that its total entry count still exceeds the system size.
    TEST(FeSolutionCheckWholeSystem)
    {
        // Two scalar unknowns on one square patch without boundary conditions: each space has all
        // mb.basis(0).size() basis functions as free dofs
        gsMultiPatch<> mp(*gsNurbsCreator<>::BSplineSquare());
        gsMultiBasis<> mb(mp);
        mb.uniformRefine();
        gsBoundaryConditions<> bcs;

        gsExprAssembler<> A(2,2);
        A.setIntegrationElements(mb);
        auto u0 = A.getSpace(mb, 1, 0);
        auto u1 = A.getSpace(mb, 1, 1);
        u0.setup(bcs, dirichlet::homogeneous, 0);
        u1.setup(bcs, dirichlet::homogeneous, 0);
        A.initSystem();

        // Whole-system vector: the unknowns are stacked, so it has numDofs() rows, and numDofs() is
        // twice the basis size because both spaces use the same basis without eliminated dofs
        gsMatrix<> x(A.numDofs(), 1);
        x.setZero();
        auto s0 = A.getSolution(u0, x);
        auto s1 = A.getSolution(u1, x);

        CHECK_EQUAL(2*mb.basis(0).size(), A.numDofs());
        CHECK(s0.check());
        CHECK(s1.check());

        // The second space ends at index numDofs(), so numDofs()-1 rows is the longest vector that is
        // still too short for it. check() reports this through a warning.
        gsMatrix<> y(A.numDofs() - 1, 1);
        auto t1 = A.getSolution(u1, y);
        UnitTest::deactivate_output();
        CHECK(!t1.check());
        UnitTest::reactivate_output();

        // Two columns with numDofs()-1 rows hold 2*(numDofs()-1) > numDofs() entries, so only the
        // row count reveals that the vector is too short
        gsMatrix<> y2(A.numDofs() - 1, 2);
        y2.setZero();
        CHECK(y2.size() > A.numDofs());
        auto t2 = A.getSolution(u1, y2);
        UnitTest::deactivate_output();
        CHECK(!t2.check());
        UnitTest::reactivate_output();
    }
}
