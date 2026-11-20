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

// Largest |a - b| over all entries, or +infinity when the sizes differ or either operand holds a
// non-finite entry. A size mismatch or a NaN then fails a CHECK on the result instead of aborting
// in Eigen's size assertion or vanishing in cwiseAbs().maxCoeff().
real_t maxAbsDiff(const gismo::gsMatrix<real_t> & a, const gismo::gsMatrix<real_t> & b)
{
    if (a.rows() != b.rows() || a.cols() != b.cols() || !a.allFinite() || !b.allFinite())
        return std::numeric_limits<real_t>::infinity();
    if (a.size() == 0)
        return real_t(0);
    return (a - b).cwiseAbs().maxCoeff();
}

// Runs the enclosing scope on one OpenMP thread and restores the previous count when the scope
// ends. Declare it before any assembler of the scope: util::gsThreaded sizes its per-thread
// storage when an assembler is constructed, and locals are destroyed in reverse order, so every
// assembler is gone before the count is raised again.
struct OneThreadScope
{
    OneThreadScope() : m_threads(omp_get_max_threads()) { omp_set_num_threads(1); }
    ~OneThreadScope() { omp_set_num_threads(m_threads); }
    const int m_threads;
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

    // Property: the boundary loops build their Gauss rule from the degree of the
    // integration domain, not from the degree of the trial basis.
    //
    // Fixture: two unit squares (affine geometry), a degree-2 trial basis mbT and a degree-3
    // integration basis mbI on identical 4x4 elements per patch. A Gauss rule has
    // cast<index_t>(quA*degree + quB + 0.5) nodes per direction; with quA = 1, quB = 0 that is
    // 3 nodes on the integration degree and 2 on the trial degree. The boundary mass of a
    // degree-2 basis has degree 4 along a side: 3 Gauss nodes are exact (up to degree 5), 2 are
    // not (up to degree 3).
    //
    // Oracle: the reference assembler R integrates over mbT with the default options (quA = 1,
    // quB = 1: 3 nodes at degree 2). A and R then evaluate the same 3 nodes per element in the
    // same element order, so their matrices agree bitwise. The witness W (integration over mbT,
    // the quB = 0 options, hence 2 nodes) must differ from R: it shows that the fixture tells the
    // two node counts apart. The entries of R sum to the length of the boundary (partition of
    // unity, 6 unit sides).
    TEST(RuleSourceBoundaryFollowsIntegrationDegree)
    {
        gsMultiPatch<real_t> mp = gsNurbsCreator<real_t>::BSplineSquareGrid(2, 1, 1.0);
        gsMultiBasis<real_t> mbT(mp);
        mbT.degreeElevate(1); mbT.uniformRefine(); mbT.uniformRefine();
        gsMultiBasis<real_t> mbI(mp);
        mbI.degreeElevate(2); mbI.uniformRefine(); mbI.uniformRefine();
        CHECK_EQUAL(2, mbT.basis(0).degree(0));
        CHECK_EQUAL(3, mbI.basis(0).degree(0));
        CHECK_EQUAL(6, static_cast<index_t>(mp.nBoundary()));

        // Neumann data on every boundary side; assembleBdr takes the boundary source from it.
        gsFunctionExpr<real_t> g("1+x+2*y", 2);
        gsBoundaryConditions<real_t> bc;
        for (gsMultiPatch<real_t>::const_biterator bit = mp.bBegin(); bit != mp.bEnd(); ++bit)
            bc.addCondition(*bit, condition_type::neumann, &g);
        bc.setGeoMap(mp);

        gsOptionList optI = gsExprAssembler<real_t>::defaultOptions();
        optI.setReal("quA", 1.0);
        optI.setInt("quB", 0);

        // Tested assembler: trial degree 2, integration degree 3, 3 nodes.
        gsExprAssembler<real_t> A(1, 1);
        A.setOptions(optI);
        A.setIntegrationElements(mbI);
        auto GA = A.getMap(mp);
        auto uA = A.getSpace(mbT);
        uA.setup(-1);

        // Reference: integration over the trial basis, default options, 3 nodes.
        gsExprAssembler<real_t> R(1, 1);
        R.setIntegrationElements(mbT);
        auto GR = R.getMap(mp);
        auto uR = R.getSpace(mbT);
        uR.setup(-1);

        // Witness: integration over the trial basis with quB = 0, 2 nodes.
        gsExprAssembler<real_t> W(1, 1);
        W.setOptions(optI);
        W.setIntegrationElements(mbT);
        auto GW = W.getMap(mp);
        auto uW = W.getSpace(mbT);
        uW.setup(-1);

        const real_t eps = math::limits::epsilon();
        const real_t inf = std::numeric_limits<real_t>::infinity();

        // Boundary mass through the list of Neumann sides.
        A.initSystem();
        A.assembleBdr(bc.get("Neumann"), uA * uA.tr() * meas(GA));
        const gsMatrix<real_t> KA = A.matrix().toDense();
        R.initSystem();
        R.assembleBdr(bc.get("Neumann"), uR * uR.tr() * meas(GR));
        const gsMatrix<real_t> KR = R.matrix().toDense();
        W.initSystem();
        W.assembleBdr(bc.get("Neumann"), uW * uW.tr() * meas(GW));
        const gsMatrix<real_t> KW = W.matrix().toDense();

        const real_t maxAbs = KR.cwiseAbs().maxCoeff();
        CHECK(maxAbs > 0);
        CHECK_CLOSE(real_t(6), KR.sum(), 1000 * eps);
        CHECK_EQUAL(KR.rows(), KA.rows());
        CHECK_EQUAL(KR.cols(), KA.cols());
        CHECK_EQUAL(real_t(0), maxAbsDiff(KA, KR));
        // 2 nodes are inexact for this mass matrix, 3 are exact.
        CHECK_EQUAL(KR.rows(), KW.rows());
        const real_t dW = maxAbsDiff(KW, KR);
        CHECK(dW > 1e-6 * maxAbs && dW < inf);

        // Boundary mass through the container of all boundary sides.
        auto eA = uA * uA.tr() * meas(GA);
        A.initSystem();
        A.assembleBdr(mp.boundaries(), eA);
        const gsMatrix<real_t> LA = A.matrix().toDense();
        auto eR = uR * uR.tr() * meas(GR);
        R.initSystem();
        R.assembleBdr(mp.boundaries(), eR);
        const gsMatrix<real_t> LR = R.matrix().toDense();

        CHECK(LR.cwiseAbs().maxCoeff() > 0);
        CHECK_CLOSE(real_t(6), LR.sum(), 1000 * eps);
        CHECK_EQUAL(LR.rows(), LA.rows());
        CHECK_EQUAL(LR.cols(), LA.cols());
        CHECK_EQUAL(real_t(0), maxAbsDiff(LA, LR));
    }

    // Property: the interface loop builds its Gauss rule from the degree of the integration
    // domain, not from the degree of the trial basis.
    //
    // Fixture and oracle as in RuleSourceBoundaryFollowsIntegrationDegree: trial degree 2,
    // integration degree 3 on identical elements (3 nodes against 2 for quA = 1, quB = 0), and a
    // reference assembler that integrates over the trial basis with the default options (also
    // 3 nodes). The interface mass of a degree-2 basis has degree 4 along the interface, so only
    // the 3-node rule is exact. The single interface has unit length, which the entries of the
    // reference sum to (partition of unity of the left patch).
    TEST(RuleSourceInterfaceFollowsIntegrationDegree)
    {
        gsMultiPatch<real_t> mp = gsNurbsCreator<real_t>::BSplineSquareGrid(2, 1, 1.0);
        gsMultiBasis<real_t> mbT(mp);
        mbT.degreeElevate(1); mbT.uniformRefine(); mbT.uniformRefine();
        gsMultiBasis<real_t> mbI(mp);
        mbI.degreeElevate(2); mbI.uniformRefine(); mbI.uniformRefine();
        CHECK_EQUAL(2, mbT.basis(0).degree(0));
        CHECK_EQUAL(3, mbI.basis(0).degree(0));
        CHECK_EQUAL(1, static_cast<index_t>(mp.nInterfaces()));

        gsOptionList optI = gsExprAssembler<real_t>::defaultOptions();
        optI.setReal("quA", 1.0);
        optI.setInt("quB", 0);

        gsExprAssembler<real_t> A(1, 1);
        A.setOptions(optI);
        A.setIntegrationElements(mbI);
        auto GA = A.getMap(mp);
        auto uA = A.getSpace(mbT);
        uA.setup(-1);

        gsExprAssembler<real_t> R(1, 1);
        R.setIntegrationElements(mbT);
        auto GR = R.getMap(mp);
        auto uR = R.getSpace(mbT);
        uR.setup(-1);

        A.initSystem();
        A.assembleIfc(mp.interfaces(), uA.left() * uA.left().tr() * nv(GA.left()).norm());
        const gsMatrix<real_t> KA = A.matrix().toDense();
        R.initSystem();
        R.assembleIfc(mp.interfaces(), uR.left() * uR.left().tr() * nv(GR.left()).norm());
        const gsMatrix<real_t> KR = R.matrix().toDense();

        CHECK(KR.cwiseAbs().maxCoeff() > 0);
        CHECK_CLOSE(real_t(1), KR.sum(), 1000 * math::limits::epsilon());
        CHECK_EQUAL(KR.rows(), KA.rows());
        CHECK_EQUAL(KR.cols(), KA.cols());
        CHECK_EQUAL(real_t(0), maxAbsDiff(KA, KR));
    }

    // Property: the finite-difference Jacobian assembly with default options (no
    // "SameQuadrature" option) builds its Gauss rule from the degree of the integration
    // domain, and does so again when the loop moves to the second patch.
    //
    // Fixture and nodes as in RuleSourceBoundaryFollowsIntegrationDegree: trial degree 2,
    // integration degree 3 on identical elements, quA = 1, quB = 0 (3 nodes against 2).
    // Oracle: the residual igrad(u) * igrad(us)^T * meas(G) is linear in the solution us, so
    // its Jacobian is the stiffness matrix igrad(u) * igrad(u)^T * meas(G). The reference
    // assembles that stiffness matrix over the trial basis with the default options (3 nodes).
    // The stiffness integrand has degree 4 in the second direction: 3 nodes are exact, 2 are
    // not. The Jacobian is a finite-difference quotient, so it matches only up to the
    // difference error; 1e-6 of the largest entry is far below the 2-node error and far above
    // the difference error.
    TEST(RuleSourceJacobianDefaultOptions)
    {
        // The finite-difference Jacobian assembly aborts on an assertion (a thread-active
        // check) when it runs on two or more threads, and an abort ends the whole test binary.
        OneThreadScope serial;

        gsMultiPatch<real_t> mp = gsNurbsCreator<real_t>::BSplineSquareGrid(2, 1, 1.0);
        gsMultiBasis<real_t> mbT(mp);
        mbT.degreeElevate(1); mbT.uniformRefine(); mbT.uniformRefine();
        gsMultiBasis<real_t> mbI(mp);
        mbI.degreeElevate(2); mbI.uniformRefine(); mbI.uniformRefine();
        CHECK_EQUAL(2, mbT.basis(0).degree(0));
        CHECK_EQUAL(3, mbI.basis(0).degree(0));

        gsOptionList optI = gsExprAssembler<real_t>::defaultOptions();
        optI.setReal("quA", 1.0);
        optI.setInt("quB", 0);

        gsExprAssembler<real_t> A(1, 1);
        A.setOptions(optI);
        A.setIntegrationElements(mbI);
        auto GA = A.getMap(mp);
        auto uA = A.getSpace(mbT);
        uA.setup(-1);

        gsExprAssembler<real_t> R(1, 1);
        R.setIntegrationElements(mbT);
        auto GR = R.getMap(mp);
        auto uR = R.getSpace(mbT);
        uR.setup(-1);

        A.initSystem();
        gsMatrix<real_t> sv = gsVector<real_t>::LinSpaced(A.numDofs(), 0.1, 1.0);
        auto us = A.getSolution(uA, sv);
        A.assembleJacobian(igrad(uA, GA) * igrad(us, GA).tr() * meas(GA), us);
        const gsMatrix<real_t> J = A.matrix().toDense();

        R.initSystem();
        R.assemble(igrad(uR, GR) * igrad(uR, GR).tr() * meas(GR));
        const gsMatrix<real_t> K = R.matrix().toDense();

        const real_t maxAbs = K.cwiseAbs().maxCoeff();
        CHECK(maxAbs > 0);
        CHECK_EQUAL(K.rows(), J.rows());
        CHECK_EQUAL(K.cols(), J.cols());
        CHECK(maxAbsDiff(J, K) <= 1e-6 * maxAbs);
    }

    // Property: quPointsWeights iterates the elements of the integration domain and builds its
    // Gauss rule from the integration degree, one point set per patch.
    //
    // Fixture: as in RuleSourceBoundaryFollowsIntegrationDegree (trial degree 2, integration
    // degree 3, 4x4 = 16 elements per patch, quA = 1, quB = 0).
    // Oracle: the integration degree gives 3 nodes per direction, i.e. 16 * 3 * 3 = 144 points
    // per patch; the trial degree would give 16 * 2 * 2 = 64. The weight sum alone cannot tell
    // the two rules apart (every rule integrates the parametric area 1 exactly), so the point
    // count carries the test.
    TEST(RuleSourceQuPointsWeightsFollowsIntegrationDomain)
    {
        gsMultiPatch<real_t> mp = gsNurbsCreator<real_t>::BSplineSquareGrid(2, 1, 1.0);
        gsMultiBasis<real_t> mbT(mp);
        mbT.degreeElevate(1); mbT.uniformRefine(); mbT.uniformRefine();
        gsMultiBasis<real_t> mbI(mp);
        mbI.degreeElevate(2); mbI.uniformRefine(); mbI.uniformRefine();
        CHECK_EQUAL(2, mbT.basis(0).degree(0));
        CHECK_EQUAL(3, mbI.basis(0).degree(0));

        gsOptionList optI = gsExprAssembler<real_t>::defaultOptions();
        optI.setReal("quA", 1.0);
        optI.setInt("quB", 0);

        gsExprAssembler<real_t> A(1, 1);
        A.setOptions(optI);
        A.setIntegrationElements(mbI);
        auto uA = A.getSpace(mbT);
        uA.setup(-1);
        A.initSystem();

        std::vector<gsMatrix<real_t> > P;
        std::vector<gsVector<real_t> > W;
        A.quPointsWeights(P, W);

        CHECK_EQUAL(2, static_cast<index_t>(W.size()));
        CHECK_EQUAL(2, static_cast<index_t>(P.size()));
        for (size_t p = 0; p != W.size() && p != P.size(); ++p)
        {
            CHECK(W[p].allFinite());
            CHECK_EQUAL(144, W[p].size());
            CHECK_EQUAL(W[p].size(), P[p].cols());
            CHECK_EQUAL(2, P[p].rows());
            CHECK_CLOSE(real_t(1), W[p].sum(), 1000 * math::limits::epsilon());
        }
    }

    // Property: on a hierarchical basis with a separate integration basis, boundary assembly
    // with quRule = 4 (gsHPatchRule) equals Gauss, the factory of gsHPatchRule::makeFactory
    // gives bitwise the matrix of the option path, and quPointsWeights iterates the elements
    // of the integration basis.
    //
    // The hierarchical fixture has one patch: hierarchical domains do not report the patch
    // index of their elements, so a multi-patch hierarchical volume loop would integrate every
    // element on patch 0.
    //
    // Fixture: one truncated hierarchical patch (quadratic, refined to level 2 in a corner) on the
    // identity geometry, with the integration basis of gsHPatchRule::makeIntegrationBasis and
    // Neumann data on all 4 sides.
    // Assemblers: Ag integrates over the trial basis with default Gauss (3 nodes, exact for the
    // products of degree at most 4 here), Ah integrates over the integration basis with quRule = 4
    // (target degree 4, regularity 1), and Af integrates over the integration basis with
    // inexact Gauss (quA = 1, quB = 0: 2 nodes) plus a counting factory that forwards to
    // gsHPatchRule::makeFactory(integMB, 4, 1, false) and ignores its own options argument.
    // (a) is a regression guard: both rules are exact, so it does not depend on which basis the
    //     rule is built on. (b) compares the factory path with the option path. With Gauss options
    //     on Af the check also fails when the factory is ignored.
    TEST(RuleSourceHPatchOnSeparateIntegrationBasis)
    {
        gsKnotVector<real_t> kv(0, 1, 0, 3, 1);
        kv.uniformRefine(); kv.uniformRefine();
        gsTensorBSplineBasis<2, real_t> tb(kv, kv);
        gsTHBSplineBasis<2, real_t> thb(tb);
        thb.uniformRefine();
        std::vector<index_t> box = {2, 0,0, 4,4};
        thb.refineElements(box);
        const gsTHBSplineBasis<2, real_t> integ = gsHPatchRule<real_t>::makeIntegrationBasis<2>(thb);

        // Identity geometry on the level-0 tensor basis.
        gsTensorBSplineBasis<2, real_t> & l0 = thb.tensorLevel(0);
        gsMatrix<real_t> coefs = l0.anchors().transpose();
        gsTensorBSpline<2, real_t> geom(l0, coefs);
        gsMultiPatch<real_t> mp;
        mp.addPatch(geom);
        mp.computeTopology();
        CHECK_EQUAL(4, static_cast<index_t>(mp.nBoundary()));
        gsMultiBasis<real_t> mb;
        mb.addBasis(thb.clone().release());
        mb.setTopology(mp);
        gsMultiBasis<real_t> integMB;
        integMB.addBasis(integ.clone().release());
        integMB.setTopology(mp);

        gsFunctionExpr<real_t> g("1+x+2*y", 2);
        gsBoundaryConditions<real_t> bc;
        for (gsMultiPatch<real_t>::const_biterator bit = mp.bBegin(); bit != mp.bEnd(); ++bit)
            bc.addCondition(*bit, condition_type::neumann, &g);
        bc.setGeoMap(mp);

        gsExprAssembler<real_t> Ag(1, 1);
        Ag.setIntegrationElements(mb);
        auto Gg = Ag.getMap(mp);
        auto ug = Ag.getSpace(mb);
        ug.setup(-1);

        gsOptionList optH = gsExprAssembler<real_t>::defaultOptions();
        optH.setInt("quRule", 4);
        optH.setReal("quA", 4.0);
        optH.setInt("quB", 1);
        gsExprAssembler<real_t> Ah(1, 1);
        Ah.setOptions(optH);
        Ah.setIntegrationElements(integMB);
        auto Gh = Ah.getMap(mp);
        auto uh = Ah.getSpace(mb);
        uh.setup(-1);

        gsOptionList optF = gsExprAssembler<real_t>::defaultOptions();
        optF.setReal("quA", 1.0);
        optF.setInt("quB", 0);
        gsExprAssembler<real_t> Af(1, 1);
        Af.setOptions(optF);
        Af.setIntegrationElements(integMB);
        auto Gf = Af.getMap(mp);
        auto uf = Af.getSpace(mb);
        uf.setup(-1);
        const gsExprAssembler<real_t>::QuadratureFactory hpatch =
            gsHPatchRule<real_t>::makeFactory(integMB, 4, 1, false);
        std::atomic<index_t> volCalls(0);
        std::atomic<index_t> bdrCalls(0);
        Af.setQuadratureFactory(
            [hpatch, &volCalls, &bdrCalls](const gsBasis<real_t> & basis,
                                           const gsOptionList & options,
                                           index_t patch,
                                           short_t fixedDirection)
                -> gsExprAssembler<real_t>::QuadratureRulePtr
            {
                if (fixedDirection == -1) ++volCalls; else ++bdrCalls;
                return hpatch(basis, options, patch, fixedDirection);
            });
        CHECK(Af.hasCustomQuadrature());

        const real_t eps = math::limits::epsilon();

        // (a) Boundary mass and boundary right-hand side: quRule = 4 equals Gauss.
        Ag.initSystem();
        Ag.assembleBdr(bc.get("Neumann"), ug * ug.tr() * meas(Gg));
        const gsMatrix<real_t> Mg = Ag.matrix().toDense();
        Ah.initSystem();
        Ah.assembleBdr(bc.get("Neumann"), uh * uh.tr() * meas(Gh));
        const gsMatrix<real_t> Mh = Ah.matrix().toDense();
        const real_t maxAbsM = Mg.cwiseAbs().maxCoeff();
        CHECK(maxAbsM > 0);
        CHECK_EQUAL(Mg.rows(), Mh.rows());
        CHECK_EQUAL(Mg.cols(), Mh.cols());
        CHECK(maxAbsDiff(Mh, Mg) <= 500 * eps * maxAbsM);

        auto gNg = Ag.getBdrFunction(Gg);
        Ag.initSystem();
        Ag.assembleBdr(bc.get("Neumann"), ug * gNg * meas(Gg));
        const gsMatrix<real_t> bg = Ag.rhs();
        auto gNh = Ah.getBdrFunction(Gh);
        Ah.initSystem();
        Ah.assembleBdr(bc.get("Neumann"), uh * gNh * meas(Gh));
        const gsMatrix<real_t> bh = Ah.rhs();
        const real_t maxAbsB = bg.cwiseAbs().maxCoeff();
        CHECK(maxAbsB > 0);
        CHECK_EQUAL(bg.rows(), bh.rows());
        CHECK_EQUAL(bg.cols(), bh.cols());
        CHECK(maxAbsDiff(bh, bg) <= 500 * eps * maxAbsB);

        // (b) The factory path is the option path, bitwise, and the factory is asked once per
        // boundary side (4 sides, fixed direction >= 0), never for a volume rule.
        index_t vol0 = volCalls.load();
        index_t bdr0 = bdrCalls.load();
        Af.initSystem();
        Af.assembleBdr(bc.get("Neumann"), uf * uf.tr() * meas(Gf));
        const gsMatrix<real_t> Mf = Af.matrix().toDense();
        CHECK_EQUAL(4, bdrCalls.load() - bdr0);
        CHECK_EQUAL(0, volCalls.load() - vol0);
        CHECK_EQUAL(Mh.rows(), Mf.rows());
        CHECK_EQUAL(Mh.cols(), Mf.cols());
        CHECK_EQUAL(real_t(0), maxAbsDiff(Mf, Mh));

        auto eh = uh * uh.tr() * meas(Gh);
        Ah.initSystem();
        Ah.assembleBdr(mp.boundaries(), eh);
        const gsMatrix<real_t> Nh = Ah.matrix().toDense();
        auto ef = uf * uf.tr() * meas(Gf);
        vol0 = volCalls.load();
        bdr0 = bdrCalls.load();
        Af.initSystem();
        Af.assembleBdr(mp.boundaries(), ef);
        const gsMatrix<real_t> Nf = Af.matrix().toDense();
        CHECK_EQUAL(4, bdrCalls.load() - bdr0);
        CHECK_EQUAL(0, volCalls.load() - vol0);
        CHECK(Nh.cwiseAbs().maxCoeff() > 0);
        CHECK_EQUAL(Nh.rows(), Nf.rows());
        CHECK_EQUAL(Nh.cols(), Nf.cols());
        CHECK_EQUAL(real_t(0), maxAbsDiff(Nf, Nh));

        // (c) quPointsWeights runs over the leaves of the integration basis; the factory is
        // asked once for the single patch (fixed direction -1) and gives the same points.
        Ah.initSystem();
        std::vector<gsMatrix<real_t> > P;
        std::vector<gsVector<real_t> > W;
        Ah.quPointsWeights(P, W);
        CHECK_EQUAL(1, static_cast<index_t>(W.size()));
        CHECK_EQUAL(1, static_cast<index_t>(P.size()));
        if (W.size() == 1 && P.size() == 1)
        {
            CHECK(W[0].allFinite());
            CHECK_EQUAL(532, W[0].size()); // all nodes of patch 0 under the merged partition
            CHECK_EQUAL(W[0].size(), P[0].cols());
            CHECK_CLOSE(real_t(1), W[0].sum(), 1000 * eps);

            std::vector<gsMatrix<real_t> > Pf;
            std::vector<gsVector<real_t> > Wf;
            vol0 = volCalls.load();
            bdr0 = bdrCalls.load();
            Af.initSystem();
            Af.quPointsWeights(Pf, Wf);
            CHECK_EQUAL(1, volCalls.load() - vol0);
            CHECK_EQUAL(0, bdrCalls.load() - bdr0);
            CHECK_EQUAL(1, static_cast<index_t>(Wf.size()));
            CHECK_EQUAL(1, static_cast<index_t>(Pf.size()));
            if (Wf.size() == 1 && Pf.size() == 1)
            {
                CHECK_EQUAL(W[0].size(), Wf[0].size());
                CHECK_EQUAL(P[0].cols(), Pf[0].cols());
                CHECK_EQUAL(real_t(0), maxAbsDiff(gsMatrix<real_t>(Wf[0]), gsMatrix<real_t>(W[0])));
                CHECK_EQUAL(real_t(0), maxAbsDiff(Pf[0], P[0]));
            }
        }
    }

    // Property: an installed quadrature factory supplies the rule in every loop (volume,
    // boundary list, boundary container, interface, Jacobian, quPointsWeights), and the
    // default-option Jacobian assembly works with it.
    //
    // Tensor-product bases are used because hierarchical domains do not report the patch index
    // of their elements, so a two-patch hierarchical volume loop would integrate every element
    // on patch 0. One thread is forced for two reasons: the finite-difference Jacobian assembly
    // aborts on an assertion at two or more threads, and the call counts of the volume loops
    // are literal only on one thread (they ask for a rule once per thread and patch).
    //
    // Fixture: two unit squares, trial = integration = degree 2 on 4x4 elements per patch.
    // Ao uses Lobatto options (quRule = 2) and no factory. Af uses default Gauss options plus a
    // factory that returns the Lobatto rule of Ao's options and ignores the options it is
    // given. Gauss and 3-node Lobatto differ on the degree-4 mass, so a loop that ignored the
    // factory would reproduce the Gauss matrix and miss Ao's. Per loop the matrices of Af and
    // Ao agree bitwise (same nodes, same element order, same finite differences), and the
    // counter of the loop's kind grows by a literal number of rule requests: volume and
    // Jacobian loops ask per patch (fixed direction -1), boundary and interface loops ask per
    // side or interface (fixed direction >= 0), quPointsWeights asks per patch.
    // The volume, boundary, interface and quPointsWeights rows are regression guards for the
    // precedence of the factory, which the option path cannot disturb; the Jacobian row also
    // needs a rule at all under default options.
    TEST(RuleSourceFactoryInEveryLoop)
    {
        OneThreadScope serial;

        gsMultiPatch<real_t> mp = gsNurbsCreator<real_t>::BSplineSquareGrid(2, 1, 1.0);
        gsMultiBasis<real_t> mbT(mp);
        mbT.degreeElevate(1); mbT.uniformRefine(); mbT.uniformRefine();
        CHECK_EQUAL(2, mbT.basis(0).degree(0));
        CHECK_EQUAL(6, static_cast<index_t>(mp.nBoundary()));
        CHECK_EQUAL(1, static_cast<index_t>(mp.nInterfaces()));

        gsFunctionExpr<real_t> g("1+x+2*y", 2);
        gsBoundaryConditions<real_t> bc;
        for (gsMultiPatch<real_t>::const_biterator bit = mp.bBegin(); bit != mp.bEnd(); ++bit)
            bc.addCondition(*bit, condition_type::neumann, &g);
        bc.setGeoMap(mp);

        gsOptionList lobattoOpts = gsExprAssembler<real_t>::defaultOptions();
        lobattoOpts.setInt("quRule", 2);

        gsExprAssembler<real_t> Ao(1, 1);
        Ao.setOptions(lobattoOpts);
        Ao.setIntegrationElements(mbT);
        auto Go = Ao.getMap(mp);
        auto uo = Ao.getSpace(mbT);
        uo.setup(-1);

        gsExprAssembler<real_t> Af(1, 1);
        Af.setIntegrationElements(mbT);
        auto Gf = Af.getMap(mp);
        auto uf = Af.getSpace(mbT);
        uf.setup(-1);

        const real_t inf = std::numeric_limits<real_t>::infinity();

        // Volume mass. Without the factory Af integrates with Gauss; the witness shows that
        // Gauss and Lobatto differ on this mass matrix.
        Ao.initSystem();
        Ao.assemble(uo * uo.tr() * meas(Go));
        const gsMatrix<real_t> Ko = Ao.matrix().toDense();
        Af.initSystem();
        Af.assemble(uf * uf.tr() * meas(Gf));
        const gsMatrix<real_t> Kgauss = Af.matrix().toDense();
        const real_t maxAbs = Ko.cwiseAbs().maxCoeff();
        CHECK(maxAbs > 0);
        CHECK_EQUAL(Ko.rows(), Kgauss.rows());
        const real_t dGauss = maxAbsDiff(Kgauss, Ko);
        CHECK(dGauss > 1e-6 * maxAbs && dGauss < inf);

        std::atomic<index_t> volCalls(0);
        std::atomic<index_t> bdrCalls(0);
        Af.setQuadratureFactory(
            [lobattoOpts, &volCalls, &bdrCalls](const gsBasis<real_t> & basis,
                                                const gsOptionList &,
                                                index_t,
                                                short_t fixedDirection)
                -> gsExprAssembler<real_t>::QuadratureRulePtr
            {
                if (fixedDirection == -1) ++volCalls; else ++bdrCalls;
                return gsQuadrature::getPtr(basis, lobattoOpts, fixedDirection);
            });

        index_t vol0 = volCalls.load();
        index_t bdr0 = bdrCalls.load();
        Af.initSystem();
        Af.assemble(uf * uf.tr() * meas(Gf));
        const gsMatrix<real_t> Kf = Af.matrix().toDense();
        CHECK_EQUAL(2, volCalls.load() - vol0);
        CHECK_EQUAL(0, bdrCalls.load() - bdr0);
        CHECK_EQUAL(Ko.rows(), Kf.rows());
        CHECK_EQUAL(Ko.cols(), Kf.cols());
        CHECK_EQUAL(real_t(0), maxAbsDiff(Kf, Ko));

        // Boundary mass over the list of Neumann sides.
        Ao.initSystem();
        Ao.assembleBdr(bc.get("Neumann"), uo * uo.tr() * meas(Go));
        const gsMatrix<real_t> Bo = Ao.matrix().toDense();
        vol0 = volCalls.load();
        bdr0 = bdrCalls.load();
        Af.initSystem();
        Af.assembleBdr(bc.get("Neumann"), uf * uf.tr() * meas(Gf));
        const gsMatrix<real_t> Bf = Af.matrix().toDense();
        CHECK_EQUAL(6, bdrCalls.load() - bdr0);
        CHECK_EQUAL(0, volCalls.load() - vol0);
        CHECK(Bo.cwiseAbs().maxCoeff() > 0);
        CHECK_EQUAL(Bo.rows(), Bf.rows());
        CHECK_EQUAL(Bo.cols(), Bf.cols());
        CHECK_EQUAL(real_t(0), maxAbsDiff(Bf, Bo));

        // Boundary mass over the container of all boundary sides.
        auto eo = uo * uo.tr() * meas(Go);
        Ao.initSystem();
        Ao.assembleBdr(mp.boundaries(), eo);
        const gsMatrix<real_t> Co = Ao.matrix().toDense();
        auto ef = uf * uf.tr() * meas(Gf);
        vol0 = volCalls.load();
        bdr0 = bdrCalls.load();
        Af.initSystem();
        Af.assembleBdr(mp.boundaries(), ef);
        const gsMatrix<real_t> Cf = Af.matrix().toDense();
        CHECK_EQUAL(6, bdrCalls.load() - bdr0);
        CHECK_EQUAL(0, volCalls.load() - vol0);
        CHECK(Co.cwiseAbs().maxCoeff() > 0);
        CHECK_EQUAL(Co.rows(), Cf.rows());
        CHECK_EQUAL(Co.cols(), Cf.cols());
        CHECK_EQUAL(real_t(0), maxAbsDiff(Cf, Co));

        // Interface mass.
        Ao.initSystem();
        Ao.assembleIfc(mp.interfaces(), uo.left() * uo.left().tr() * nv(Go.left()).norm());
        const gsMatrix<real_t> Io = Ao.matrix().toDense();
        vol0 = volCalls.load();
        bdr0 = bdrCalls.load();
        Af.initSystem();
        Af.assembleIfc(mp.interfaces(), uf.left() * uf.left().tr() * nv(Gf.left()).norm());
        const gsMatrix<real_t> If = Af.matrix().toDense();
        CHECK_EQUAL(1, bdrCalls.load() - bdr0);
        CHECK_EQUAL(0, volCalls.load() - vol0);
        CHECK(Io.cwiseAbs().maxCoeff() > 0);
        CHECK_EQUAL(Io.rows(), If.rows());
        CHECK_EQUAL(Io.cols(), If.cols());
        CHECK_EQUAL(real_t(0), maxAbsDiff(If, Io));

        // Finite-difference Jacobian under default options.
        Ao.initSystem();
        gsMatrix<real_t> svo = gsVector<real_t>::LinSpaced(Ao.numDofs(), 0.1, 1.0);
        auto uso = Ao.getSolution(uo, svo);
        Ao.assembleJacobian(igrad(uo, Go) * igrad(uso, Go).tr() * meas(Go), uso);
        const gsMatrix<real_t> Jo = Ao.matrix().toDense();
        vol0 = volCalls.load();
        bdr0 = bdrCalls.load();
        Af.initSystem();
        gsMatrix<real_t> svf = gsVector<real_t>::LinSpaced(Af.numDofs(), 0.1, 1.0);
        auto usf = Af.getSolution(uf, svf);
        Af.assembleJacobian(igrad(uf, Gf) * igrad(usf, Gf).tr() * meas(Gf), usf);
        const gsMatrix<real_t> Jf = Af.matrix().toDense();
        CHECK_EQUAL(2, volCalls.load() - vol0);
        CHECK_EQUAL(0, bdrCalls.load() - bdr0);
        CHECK(Jo.cwiseAbs().maxCoeff() > 0);
        CHECK_EQUAL(Jo.rows(), Jf.rows());
        CHECK_EQUAL(Jo.cols(), Jf.cols());
        CHECK_EQUAL(real_t(0), maxAbsDiff(Jf, Jo));

        // Quadrature points and weights, one set per patch.
        Ao.initSystem();
        std::vector<gsMatrix<real_t> > Po;
        std::vector<gsVector<real_t> > Wo;
        Ao.quPointsWeights(Po, Wo);
        vol0 = volCalls.load();
        bdr0 = bdrCalls.load();
        Af.initSystem();
        std::vector<gsMatrix<real_t> > Pf;
        std::vector<gsVector<real_t> > Wf;
        Af.quPointsWeights(Pf, Wf);
        CHECK_EQUAL(2, volCalls.load() - vol0);
        CHECK_EQUAL(0, bdrCalls.load() - bdr0);
        CHECK_EQUAL(2, static_cast<index_t>(Wo.size()));
        CHECK_EQUAL(2, static_cast<index_t>(Po.size()));
        CHECK_EQUAL(Wo.size(), Wf.size());
        CHECK_EQUAL(Po.size(), Pf.size());
        for (size_t p = 0; p != Wo.size() && p != Wf.size() && p != Po.size() && p != Pf.size(); ++p)
        {
            CHECK(Wo[p].size() > 0);
            CHECK_EQUAL(Wo[p].size(), Wf[p].size());
            CHECK_EQUAL(Po[p].cols(), Pf[p].cols());
            CHECK_EQUAL(real_t(0), maxAbsDiff(gsMatrix<real_t>(Wf[p]), gsMatrix<real_t>(Wo[p])));
            CHECK_EQUAL(real_t(0), maxAbsDiff(Pf[p], Po[p]));
        }
    }
}
