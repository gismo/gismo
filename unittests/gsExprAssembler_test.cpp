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

/// Rule that maps an element to zero nodes when the midpoint of its
/// parametric box lies below \a cut in direction \a dir, and to a tensor
/// Gauss rule otherwise -- models an immersed rule on a fully exterior cell.
class EmptyBelowRule : public gismo::gsQuadRule<real_t>
{
public:
    EmptyBelowRule(const gsVector<index_t> & numNodes, short_t dir, real_t cut)
    : m_gauss(numNodes), m_dir(dir), m_cut(cut) { }

    using gismo::gsQuadRule<real_t>::mapTo;

    void mapTo(const gsVector<real_t> & lower, const gsVector<real_t> & upper,
               gsMatrix<real_t> & nodes, gsVector<real_t> & weights) const override
    {
        if (0.5*(lower[m_dir] + upper[m_dir]) < m_cut)
        {
            nodes.resize(lower.size(), 0);
            weights.resize(0);
            return;
        }
        m_gauss.mapTo(lower, upper, nodes, weights);
    }
private:
    gsGaussRule<real_t> m_gauss;
    short_t m_dir;
    real_t  m_cut;
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

    // A vector-valued expression equal to the geometry map must reproduce the
    // stored physical points node by node
    TEST(EvaluatorEvalAtInterface)
    {
        gsMultiPatch<> mp = gsNurbsCreator<>::BSplineSquareGrid(1,2,1);
        gsMultiBasis<> mb(mp);
        mb.uniformRefine();
        mb.uniformRefine();

        gsExprEvaluator<> ev;
        ev.setIntegrationDomain(mb.domain());
        gsExprEvaluator<>::geometryMap G = ev.getMap(mp);

        CHECK(1==mp.interfaces().size());
        ev.evalAtInterface(G, G, mp.interfaces());

        // 2 refinements: 4 elements along the interface, 2 Gauss nodes per
        // direction for degree 1
        const index_t nElem = 4;
        const index_t nNodes = 2;
        const gsMatrix<> pts = ev.allPoints();
        CHECK_EQUAL(2, pts.rows());
        CHECK_EQUAL(nElem * nNodes, pts.cols());
        CHECK_EQUAL((size_t)(2 * nElem * nNodes), ev.nValues());
        if (2 != pts.rows() || (size_t)pts.size() != ev.nValues())
            return;

        const gsMatrix<> vals = ev.allValues(2, pts.cols());
        CHECK_ARRAY_CLOSE(pts.data(), vals.data(), pts.size(), 1e-12);
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

    // An element whose custom quadrature rule has no nodes contributes
    // nothing and must not be precomputed: compute_impl, reached through
    // integral()/integralElWise()/min(), is the volume loop that must skip it.
    TEST(EvaluatorEmptyRuleVolume)
    {
        gsMultiPatch<> mp(*gsNurbsCreator<>::BSplineSquare()); // identity map on [0,1]^2
        gsMultiBasis<> mb(mp);
        mb.uniformRefine();
        mb.uniformRefine();
        CHECK_EQUAL((size_t)16, mb.basis(0).numElements());

        gsExprEvaluator<> ev;
        ev.setIntegrationElements(mb);
        gsExprEvaluator<>::geometryMap G = ev.getMap(mp);

        gsFunctionExpr<> fexpr("x^2*y", 2);
        gsFunctionExpr<> fxexpr("x", 2);
        auto f  = ev.getVariable(fexpr,  G);
        auto fx = ev.getVariable(fxexpr, G);

        // Control arm: standard quadrature, no empty elements.
        CHECK_CLOSE(1.0, ev.integral(meas(G)), 1e-13);
        CHECK(ev.min(fx) < 0.5);

        std::atomic<index_t> calls(0);
        std::atomic<bool> contextIsCorrect(true);
        ev.setQuadratureFactory(
            [&calls, &contextIsCorrect](const gsDomain<real_t> & domain,
                                        const gsOptionList &,
                                        index_t patch,
                                        short_t fixedDirection)
                -> gsExprEvaluator<real_t>::QuadratureRulePtr
            {
                ++calls;
                if (domain.numElements() != (size_t)16 ||
                    fixedDirection != -1 || patch != 0)
                    contextIsCorrect = false;
                gsVector<index_t> nn;
                nn.setConstant(domain.dim(), 3);
                if (fixedDirection >= 0) nn[fixedDirection] = 1;
                return gsExprEvaluator<real_t>::QuadratureRulePtr(
                    new EmptyBelowRule(nn, 0, 0.5));
            });
        CHECK(ev.hasCustomQuadrature());

        // Custom arm: elements with midpoint x in {0.125, 0.375} are empty
        // (8 of the 16 elements); the rest is a plain 3x3 Gauss rule.
        CHECK_CLOSE(0.5, ev.integral(meas(G)), 1e-13);
        CHECK_CLOSE(7.0/48.0, ev.integral(f*meas(G)), 1e-13);

        ev.integralElWise(f*meas(G));
        CHECK_EQUAL((size_t)16, ev.elementwise().size());
        index_t zeroCount = 0, positiveCount = 0;
        real_t sum = 0;
        for (const real_t & v : ev.elementwise())
        {
            if (v == 0.0) ++zeroCount;
            else if (v > 0.0) ++positiveCount;
            sum += v;
        }
        CHECK_EQUAL(8, zeroCount);
        CHECK_EQUAL(8, positiveCount);
        CHECK_CLOSE(7.0/48.0, sum, 1e-13);
        CHECK_CLOSE(7.0/48.0, ev.value(), 1e-13);

        // Every surviving node has x > 0.5, unlike the control arm's minimum.
        CHECK(ev.min(fx) > 0.5);

        // One factory call per thread per computation: integral, integral,
        // integralElWise and min.
#ifdef _OPENMP
        CHECK(calls.load() >= 4);
#else
        CHECK_EQUAL(4, calls.load());
#endif
        CHECK(contextIsCorrect.load());
    }

    // computeBdr_impl, reached through integralBdr(expr, bdrlist), is the
    // boundary loop that must skip an element whose rule has no nodes.
    TEST(EvaluatorEmptyRuleBoundary)
    {
        gsMultiPatch<> mp(*gsNurbsCreator<>::BSplineSquare());
        gsMultiBasis<> mb(mp);
        mb.uniformRefine();
        mb.uniformRefine();
        CHECK_EQUAL((size_t)16, mb.basis(0).numElements());

        gsExprEvaluator<> ev;
        ev.setIntegrationElements(mb);
        gsExprEvaluator<>::geometryMap G = ev.getMap(mp);

        gsFunctionExpr<> fexpr("x^2*y", 2);
        auto f = ev.getVariable(fexpr, G);

        gsExprEvaluator<>::bContainer sides;
        sides.push_back(patchSide(0, boundary::west));
        sides.push_back(patchSide(0, boundary::east));
        sides.push_back(patchSide(0, boundary::south));
        sides.push_back(patchSide(0, boundary::north));

        // Control arm: standard quadrature over all four sides.
        CHECK_CLOSE(4.0, ev.integralBdr(meas(G), sides), 1e-13);

        std::atomic<index_t> calls(0);
        std::atomic<bool> contextIsCorrect(true);
        ev.setQuadratureFactory(
            [&calls, &contextIsCorrect](const gsDomain<real_t> & domain,
                                        const gsOptionList &,
                                        index_t patch,
                                        short_t fixedDirection)
                -> gsExprEvaluator<real_t>::QuadratureRulePtr
            {
                ++calls;
                if (domain.numElements() != (size_t)16 || patch != 0 ||
                    (fixedDirection != 0 && fixedDirection != 1))
                    contextIsCorrect = false;
                gsVector<index_t> nn;
                nn.setConstant(domain.dim(), 3);
                if (fixedDirection >= 0) nn[fixedDirection] = 1;
                return gsExprEvaluator<real_t>::QuadratureRulePtr(
                    new EmptyBelowRule(nn, 0, 0.5));
            });
        CHECK(ev.hasCustomQuadrature());

        // West (u=0) is empty on every element, east (u=1) is kept in full,
        // south/north keep only x in [0.5,1].
        CHECK_CLOSE(2.0, ev.integralBdr(meas(G), sides), 1e-13);
        CHECK_EQUAL(4, calls.load());

        CHECK_CLOSE(19.0/24.0, ev.integralBdr(f*meas(G), sides), 1e-13);
        CHECK_EQUAL(8, calls.load());

        CHECK(contextIsCorrect.load());
    }

    // computeInterface_impl, reached through integralInterface(expr, iFaces),
    // is the interface loop that must skip an element whose rule has no
    // nodes.
    TEST(EvaluatorEmptyRuleInterface)
    {
        gsMultiPatch<> mp = gsNurbsCreator<>::BSplineSquareGrid(2,1,0.5);
        CHECK(1==mp.interfaces().size());

        gsMultiBasis<> mb(mp);
        mb.uniformRefine();
        mb.uniformRefine();
        CHECK_EQUAL((size_t)16, mb.basis(0).numElements());

        gsExprEvaluator<> ev;
        ev.setIntegrationElements(mb);
        gsExprEvaluator<>::geometryMap G = ev.getMap(mp);

        gsFunctionExpr<> oneExpr("1", 2);
        gsFunctionExpr<> fyExpr ("y", 2);
        auto one = ev.getVariable(oneExpr, G);
        auto fy  = ev.getVariable(fyExpr,  G);

        // Control arm: standard quadrature on the single interface.
        CHECK_CLOSE(1.0,  ev.integralInterface(one.left(), mp.interfaces()), 1e-13);
        CHECK_CLOSE(0.25, ev.integralInterface(fy.left(),  mp.interfaces()), 1e-13);

        std::atomic<index_t> calls(0);
        std::atomic<bool> contextIsCorrect(true);
        ev.setQuadratureFactory(
            [&calls, &contextIsCorrect](const gsDomain<real_t> & domain,
                                        const gsOptionList &,
                                        index_t patch,
                                        short_t fixedDirection)
                -> gsExprEvaluator<real_t>::QuadratureRulePtr
            {
                ++calls;
                if (domain.numElements() != (size_t)16 || fixedDirection != 0 ||
                    (patch != 0 && patch != 1))
                    contextIsCorrect = false;
                gsVector<index_t> nn;
                nn.setConstant(domain.dim(), 3);
                if (fixedDirection >= 0) nn[fixedDirection] = 1;
                return gsExprEvaluator<real_t>::QuadratureRulePtr(
                    new EmptyBelowRule(nn, 1, 0.5));
            });
        CHECK(ev.hasCustomQuadrature());

        // Interface elements with v-midpoint in {0.125, 0.375} are empty.
        CHECK_CLOSE(0.5,    ev.integralInterface(one.left(), mp.interfaces()), 1e-13);
        CHECK_CLOSE(0.1875, ev.integralInterface(fy.left(),  mp.interfaces()), 1e-13);

        // One factory call per interface per computation (first side only):
        // two integralInterface calls on the single interface.
        CHECK_EQUAL(2, calls.load());
        CHECK(contextIsCorrect.load());
    }

    // The elementwise vector of an evaluator with an empty-rule factory is
    // indexed by the element id of the integration domain: skipped elements
    // hold exactly zero for every reduction (sum, min, max) and the global
    // value excludes them.
    TEST(EvaluatorEmptyRuleElementwiseAlignment)
    {
        gsMultiPatch<> mp(*gsNurbsCreator<>::BSplineSquare()); // identity map on [0,1]^2
        gsMultiBasis<> mb(mp);
        mb.uniformRefine();
        mb.uniformRefine();

        gsDomain<real_t>::Ptr dom = mb.domain();
        CHECK_EQUAL((size_t)16, dom->numElements());

        gsExprEvaluator<> ev;
        ev.setIntegrationDomain(dom);
        gsExprEvaluator<>::geometryMap G = ev.getMap(mp);

        gsFunctionExpr<> fexpr("x^2*y", 2);
        gsFunctionExpr<> fxexpr("x", 2);
        auto f  = ev.getVariable(fexpr,  G);
        auto fx = ev.getVariable(fxexpr, G);

        std::atomic<bool> contextIsCorrect(true);
        ev.setQuadratureFactory(
            [&contextIsCorrect](const gsDomain<real_t> & domain,
               const gsOptionList &,
               index_t patch,
               short_t fixedDirection)
                -> gsExprEvaluator<real_t>::QuadratureRulePtr
            {
                if (domain.numElements() != (size_t)16 ||
                    fixedDirection != -1 || patch != 0)
                    contextIsCorrect = false;
                gsVector<index_t> nn;
                nn.setConstant(domain.dim(), 3);
                if (fixedDirection >= 0) nn[fixedDirection] = 1;
                return gsExprEvaluator<real_t>::QuadratureRulePtr(
                    new EmptyBelowRule(nn, 0, 0.5));
            });

        const real_t I = ev.integralElWise(f*meas(G));
        const std::vector<real_t> eInt = ev.elementwise();
        CHECK_EQUAL((size_t)16, ev.nValues());

        ev.minElWise(fx);
        const real_t mn = ev.value();
        const std::vector<real_t> eMin = ev.elementwise();
        CHECK_EQUAL((size_t)16, ev.nValues());

        ev.maxElWise(fx);
        const real_t mx = ev.value();
        const std::vector<real_t> eMax = ev.elementwise();
        CHECK_EQUAL((size_t)16, ev.nValues());

        CHECK_EQUAL((size_t)16, eInt.size());
        CHECK_EQUAL((size_t)16, eMin.size());
        CHECK_EQUAL((size_t)16, eMax.size());

        index_t skipped = 0;
        real_t sum = 0;
        for (auto & el : dom->allElements())
        {
            const index_t id = el.id();
            const real_t a = el.lowerCorner()[0], b = el.upperCorner()[0];
            const real_t c = el.lowerCorner()[1], d = el.upperCorner()[1];
            if (0.5*(a + b) < 0.5)
            {
                ++skipped;
                CHECK_EQUAL((real_t)0, eInt[id]);
                CHECK_EQUAL((real_t)0, eMin[id]);
                CHECK_EQUAL((real_t)0, eMax[id]);
            }
            else
            {
                // Exact element integral of x^2*y over [a,b]x[c,d].
                CHECK_CLOSE((b*b*b - a*a*a)/3 * (d*d - c*c)/2, eInt[id], 1e-14);
                // Gauss nodes are interior to the element.
                CHECK(a < eMin[id] && eMin[id] < b);
                CHECK(a < eMax[id] && eMax[id] < b);
            }
            sum += eInt[id];
        }
        CHECK_EQUAL(8, skipped);
        CHECK_CLOSE(7.0/48.0, I, 1e-13);
        CHECK_CLOSE(I, sum, 1e-14);
        // The global extrema exclude the skipped zeros.
        CHECK(mn > 0.5);
        CHECK(mx < 1.0);
        CHECK(contextIsCorrect.load());
    }

    // Skeleton and ghost face loops always use plain Gauss rules: neither
    // the assembler's nor the evaluator's quadrature factory is consulted.
    TEST(FaceLoopsIgnoreQuadratureFactories)
    {
        gsMultiPatch<real_t> mp;
        mp.addPatch(gsNurbsCreator<real_t>::BSplineSquare(1,0,0));
        gsMultiBasis<real_t> mb(mp);
        mb.uniformRefine(3); // 4 x 4 elements, degree 1
        gsBoundaryConditions<real_t> bc;
        bc.setGeoMap(mp);

        // Reference arm: no factories.
        gsExprAssembler<real_t> A(1,1);
        A.setIntegrationElements(mb);
        auto G = A.getMap(mp);
        auto u = A.getSpace(mb);
        u.setup(bc, dirichlet::interpolation, 0);
        auto dJ = dnk(u.jump(), G.left(), 1); // C^0 basis: the k=1 jump is non-zero
        A.initSystem();
        A.assembleSkeleton(dJ*dJ.tr());
        const gsSparseMatrix<real_t> K0 = A.matrix();
        CHECK(K0.norm() > 1e-3);

        gsExprEvaluator<real_t> ev;
        ev.setIntegrationElements(mb);
        auto Ge = ev.getMap(mp);
        const real_t s0 = ev.integralSkeleton(nv(Ge.left()).norm());
        CHECK_CLOSE(6.0, s0, 1e-12); // 24 interior faces of length 1/4
        const std::vector<real_t> e0 = ev.elementwise();
        CHECK_EQUAL((size_t)24, e0.size());

        // Factory arm: fresh instances with both factories installed. The
        // positive controls run before any two-sided expression is registered.
        gsExprAssembler<real_t> A1(1,1);
        A1.setIntegrationElements(mb);
        auto G1 = A1.getMap(mp);
        auto u1 = A1.getSpace(mb);
        u1.setup(bc, dirichlet::interpolation, 0);
        auto dJ1 = dnk(u1.jump(), G1.left(), 1);

        gsExprEvaluator<real_t> ev1;
        ev1.setIntegrationElements(mb);
        auto Ge1 = ev1.getMap(mp);

        std::atomic<index_t> aCalls(0), eCalls(0);
        std::atomic<bool> eContextIsCorrect(true);
        A1.setQuadratureFactory(
            [&aCalls](const gsBasis<real_t> & basis,
                      const gsOptionList &,
                      index_t,
                      short_t)
                -> gsExprAssembler<real_t>::QuadratureRulePtr
            {
                ++aCalls;
                return gsExprAssembler<real_t>::QuadratureRulePtr(
                    new MidpointRule<real_t>(basis.dim()));
            });
        ev1.setQuadratureFactory(
            [&eCalls, &eContextIsCorrect](const gsDomain<real_t> & domain,
                      const gsOptionList &,
                      index_t patch,
                      short_t fixedDirection)
                -> gsExprEvaluator<real_t>::QuadratureRulePtr
            {
                ++eCalls;
                if (domain.numElements() != (size_t)16 ||
                    fixedDirection != -1 || patch != 0)
                    eContextIsCorrect = false;
                gsVector<index_t> nn;
                nn.setConstant(domain.dim(), 3);
                if (fixedDirection >= 0) nn[fixedDirection] = 1;
                return gsExprEvaluator<real_t>::QuadratureRulePtr(
                    new EmptyBelowRule(nn, 0, 0.5));
            });

        // Both factories are live for volume loops.
        A1.initSystem();
        A1.assemble(u1*u1.tr());
        CHECK(aCalls.load() > 0);
        CHECK_CLOSE(0.5, ev1.integral(meas(Ge1)), 1e-13);
        CHECK(eCalls.load() > 0);

        const index_t a0 = aCalls.load(), c0 = eCalls.load();

        A1.initSystem();
        A1.assembleSkeleton(dJ1*dJ1.tr());
        CHECK_EQUAL(a0, aCalls.load());
        CHECK_EQUAL((real_t)0, (A1.matrix() - K0).norm());

        const real_t s1 = ev1.integralSkeleton(nv(Ge1.left()).norm());
        CHECK_EQUAL(s0, s1);
        CHECK(ev1.elementwise() == e0);
        CHECK_EQUAL(c0, eCalls.load());
        CHECK(eContextIsCorrect.load());
    }

    // assembleJacobian under default options must evaluate each element on
    // its own patch; the residual is linear in the solution, so its Jacobian
    // is the stiffness matrix.
    TEST(JacobianDefaultOptionsMultipatch)
    {
        gsMultiPatch<> mp;
        mp.addPatch(gsNurbsCreator<>::BSplineSquare(1.0, 0.0, 0.0));
        mp.addPatch(gsNurbsCreator<>::BSplineQuarterAnnulus());
        gsMultiBasis<> mb(mp);
        mb.uniformRefine();

        gsExprAssembler<> A(1,1);
        A.setIntegrationElements(mb);
        auto G = A.getMap(mp);
        auto v = A.getSpace(mb);
        gsBoundaryConditions<> bc;
        v.setup(bc, dirichlet::interpolation, 0);
        A.initSystem();

        A.assemble(igrad(v,G)*igrad(v,G).tr()*meas(G));
        const gsMatrix<> K = A.matrix().toDense();
        CHECK(K.norm() > 0);

        gsMatrix<> coefs(A.numDofs(), 1);
        for (index_t i = 0; i != A.numDofs(); ++i)
            coefs(i,0) = math::sin(1.7*(real_t)i + 0.3);
        auto s = A.getSolution(v, coefs);

        A.initSystem();
        A.assembleJacobian(igrad(v,G)*igrad(s,G).tr()*meas(G), s);
        gsMatrix<> J = A.matrix().toDense();
        CHECK(K.rows() == J.rows() && K.cols() == J.cols());
        CHECK((J - K).norm() <= 1e-8 * K.norm());

        A.options().addSwitch("SameQuadrature", "", false);
        A.initSystem();
        A.assembleJacobian(igrad(v,G)*igrad(s,G).tr()*meas(G), s);
        J = A.matrix().toDense();
        CHECK((J - K).norm() <= 1e-8 * K.norm());
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
