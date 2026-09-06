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
        // Regression test for a469c2d04: _blockDims/resetDimensions used
        // dim()*mapper.freeSize() for block sizes, but freeSize() already
        // includes dim(), so a space with dim>1 doubled up on its own
        // dimension in the row/col block sizes and in the shift applied to
        // later blocks. A single space of dim 1 cannot expose this (the
        // erroneous factor is 1), so this test needs two spaces of
        // different, non-trivial dimension sharing one assembler, matching
        // the "vector space v and scalar space q" example from the fix.
        gsBSplineBasis<real_t> bb(0.0, 1.0, 3, 3);
        gsMultiBasis<real_t> mb(bb);
        gsBoundaryConditions<real_t> bcs;

        gsExprAssembler<real_t> A(2, 2);
        A.setIntegrationElements(mb);
        auto v = A.getSpace(mb, 3, 0); // vector-valued space, dim 3
        auto q = A.getSpace(mb, 1, 1); // scalar space, dim 1
        v.setup(bcs, dirichlet::homogeneous, 0);
        q.setup(bcs, dirichlet::homogeneous, 0);
        A.initSystem();

        // Expectation computed directly from the per-space dof mappers,
        // independent of both _blockDims and numDofs()/matrix() sizing:
        // freeSize() already reports the total (component-inclusive) dof
        // count for that space, so the system size is simply their sum.
        const index_t expected = v.mapper().freeSize() + q.mapper().freeSize();

        // Measured: 28 here (3*7 + 7). Before the fix numDofs() reported 70,
        // the shift for q's block having been computed as 3*(3*7) instead of
        // 3*7.
        CHECK_EQUAL(expected, A.numDofs());

        // The system matrix is sized by MatrixSizedAfterInitSystem below; the
        // check here is confined to the block arithmetic. _blockDims itself,
        // which feeds blockView(), stays uncovered -- the check above reaches
        // the same defect through resetDimensions.
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
