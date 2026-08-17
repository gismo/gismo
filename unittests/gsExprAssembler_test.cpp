/** @file gsExprAssembler_test.cpp

    @brief Tests for gsExprAssembler

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): R. Schneckenleitner
*/

#include "gismo_unittest.h"

#include <gsDomain/gsImplicitTrimmedDomain.h>


SUITE(gsExprAssembler_test)
{
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

    /// Test 4 -- gsExprAssembler integrating over a gsImplicitTrimmedDomain
    /// whose stored m_deg is 1, with a degree-3 space, sizes quadrature from
    /// the space (via the degrees forwarded through setIntegrationDomain's
    /// path), not from the trimmed domain's own (wrong) guessed degree.
    TEST(TrimmedDomainUsesSpaceDegree)
    {
        typedef gsExprAssembler<>::geometryMap geometryMap;
        typedef gsExprAssembler<>::space       space;

        gsMultiPatch<> patches = gsNurbsCreator<>::BSplineSquareGrid(1,1,1); // unit square
        gsMultiBasis<> mb(patches);
        mb.degreeElevate(2);        // degree 1 -> degree 3
        mb.uniformRefine();         // 2 elements/direction
        mb.uniformRefine();         // 4 elements/direction

        // Assert the element grid assumed below, so a change in
        // gsNurbsCreator/uniformRefine semantics fails loudly rather than
        // silently mismatching the trimmed grid.
        CHECK_EQUAL(3, mb.basis(0).degree(0));
        CHECK_EQUAL(3, mb.basis(0).degree(1));
        CHECK_EQUAL((size_t)16, mb.basis(0).numElements());

        // Trimmed domain over the SAME 4x4 grid, m_deg = 1, fully-interior
        // level set (constant -1 everywhere -> no cell is trimmed).
        gsFunctionExpr<real_t> phiConst("-1", 2);
        gsMatrix<real_t> bbox(2,2);
        bbox << 0.0, 0.0,   // lower corners
                1.0, 1.0;   // upper corners
        gsVector<index_t,2> nc; nc << 4, 4;
        memory::shared_ptr<gsImplicitTrimmedDomain<2,real_t> > trDom =
            memory::make_shared(new gsImplicitTrimmedDomain<2,real_t>(phiConst, bbox, nc, 5, /*deg=*/1));
        CHECK_EQUAL((short_t)1, trDom->degree(0));
        CHECK_EQUAL((short_t)1, trDom->degree(1));

        gsMatrix<> A_trim, A_fallback, A_untrimmed;

        // A_trim: trimmed domain, default quA=1, quB=1 -> 4 nodes/direction
        // once the space degree (3) is wired through.
        {
            gsExprAssembler<> A(1,1);
            A.setIntegrationDomain(trDom);
            geometryMap G = A.getMap(patches);
            space u = A.getSpace(mb);
            A.initSystem();
            A.assemble(u * u.tr() * meas(G));
            A_trim = A.matrix().toDense();
        }

        // A_fallback: same trimmed domain, but quA=0, quB=2 -> exactly the
        // 2 nodes/direction that the old domain-degree fallback (deg=1,
        // 1*1+1=2) produced under the DEFAULT options. Reproducing the old
        // count through the new API keeps the test self-contained.
        {
            gsExprAssembler<> A(1,1);
            A.setIntegrationDomain(trDom);
            A.options().setReal("quA", 0.0);
            A.options().setInt ("quB", 2);
            geometryMap G = A.getMap(patches);
            space u = A.getSpace(mb);
            A.initSystem();
            A.assemble(u * u.tr() * meas(G));
            A_fallback = A.matrix().toDense();
        }

        // A_untrimmed: the plain basis domain, default quA=1, quB=1.
        {
            gsExprAssembler<> A(1,1);
            A.setIntegrationElements(mb);
            geometryMap G = A.getMap(patches);
            space u = A.getSpace(mb);
            A.initSystem();
            A.assemble(u * u.tr() * meas(G));
            A_untrimmed = A.matrix().toDense();
        }

        // A_trim must differ from A_fallback: 2-point Gauss cannot integrate
        // the degree-6 mass integrand exactly, 4-point can. If they are
        // equal, the degree wiring is not live and this check must fail.
        const real_t diffFallback = (A_trim - A_fallback).norm() / A_fallback.norm();
        CHECK(diffFallback > 1e-8);

        // A_trim must equal A_untrimmed: both integrate the same 4x4 element
        // grid with the same 4-node rule, and the constant -1 level set
        // makes every cell interior, so there is no trimming.
        const real_t diffUntrimmed = (A_trim - A_untrimmed).norm() / A_untrimmed.norm();
        CHECK(diffUntrimmed < 1e-12);

        // Cheap and strong: partition of unity on the unit square with an
        // identity geometry map -> total mass is the domain area (1.0).
        CHECK_CLOSE(1.0, A_trim.sum(), 1e-12);
    }
}
