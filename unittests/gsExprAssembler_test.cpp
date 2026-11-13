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

} // anonymous namespace


SUITE(gsExprAssembler_test)
{
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
            // matrix, so a wrong offset for p's block (computed pre-fix as
            // d*(d*n) instead of d*n) would either write out of range or
            // leave part of v's block untouched.
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

        A.initSystem();
        A.assemble(lapl(u) * lapl(u).tr() * meas(G) + (lapl(u) + u) * u.tr() * meas(G),
                   lapl(u) * f * meas(G));
        const gsSparseMatrix<> M = A.matrix();
        const gsMatrix<>       b = A.rhs();

        A.initSystem();
        A.assemble(hess(u).trace() * hess(u).trace().tr() * meas(G)
                   + (hess(u).trace() + u) * u.tr() * meas(G),
                   hess(u).trace() * f * meas(G));

        CHECK(M.norm() > 1);
        CHECK(b.norm() > 1e-3);
        CHECK((M - A.matrix()).norm() <= 1e-12 * M.norm());
        CHECK((b - A.rhs()).norm()    <= 1e-12 * b.norm());
    }

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

        const index_t methods[2] = {dirichlet::interpolation, dirichlet::l2Projection};
        for (index_t m = 0; m != 2; ++m)
        {
            gsExprAssembler<> A(1,1), B(1,1);
            A.setIntegrationElements(mb);
            B.setIntegrationElements(mb);
            auto u = A.getSpace(mb, 2);
            auto w = B.getSpace(mb, 2);
            u.setup(bcs,   methods[m], 0);
            w.setup(split, methods[m], 0);
            CHECK_EQUAL(w.mapper().boundarySize(), u.mapper().boundarySize());

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
            CHECK(numFixed > 0);
            CHECK(maxVal > 0.5);
            CHECK(maxDiff <= 1e-10 * maxVal);
        }
    }

    // The side of a curve is its end point, where the L2 projection of
    // Dirichlet data reduces to the value there: the dof of the end basis
    // function takes the data value, whatever the tangent of the curve.
    // G(t) = t^2, on the line and as a curve in the plane, has a vanishing
    // tangent at t = 0; x = 2t is a regular curve. Both the expression
    // assembler and gsAssembler project.
    TEST(DirichletL2ProjectionAtCurveEnd)
    {
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
            CHECK_CLOSE(gW, u.fixedPart().at(u.mapper().bindex(0  , 0)), 1e-10);
            CHECK_CLOSE(gE, u.fixedPart().at(u.mapper().bindex(n-1, 0)), 1e-10);

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

        A.assembleJacobianIfc(patches.interfaces(), u * s * nv(G).norm(), s);
        const gsSparseMatrix<> J = A.matrix();

        A.initSystem();
        A.assembleIfc(patches.interfaces(), u * u.tr() * nv(G).norm());
        const gsSparseMatrix<> M = A.matrix();

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

            // J of the residual minus M of the linearized form, relative to M
            gsSparseMatrix<> J, M;
            auto relErr = [&]() { return (J - M).norm() / M.norm(); };

            A.initSystem();
            A.assembleJacobianIfc(patches.interfaces(), u * s.right() * nv(G).norm(), s);
            J = A.matrix();
            A.initSystem();
            A.assembleIfc(patches.interfaces(), u * u.right().tr() * nv(G).norm());
            M = A.matrix();
            CHECK(relErr() <= 1e-8);

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
            const auto check = [&]()
            {
                CHECK(M.norm() > 0);
                CHECK((J - M).norm() <= 1e-8 * M.norm());
                CHECK_EQUAL(usedCols(M, 1e-8 * M.norm()), usedCols(J, 1e-8 * M.norm()));
            };

            A.initSystem();
            A.assembleJacobianIfc(patches.interfaces(), u * s.right() * nv(G).norm(), s);
            J = A.matrix();
            A.initSystem();
            A.assembleIfc(patches.interfaces(), u * u.right().tr() * nv(G).norm());
            M = A.matrix();
            check();

            A.initSystem();
            A.assembleJacobianIfc(patches.interfaces(), u.right() * s * nv(G).norm(), s);
            J = A.matrix();
            A.initSystem();
            A.assembleIfc(patches.interfaces(), u.right() * u.tr() * nv(G).norm());
            M = A.matrix();
            check();

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
}
