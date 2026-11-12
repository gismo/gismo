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
