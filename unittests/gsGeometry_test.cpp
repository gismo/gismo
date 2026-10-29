/** @file gsGeometry_test.cpp

    @brief Tsting gsGeometry

    == BASIC REFERENCE ==
         - TEST(NAME_OF_TEST) { body_of_test }
         - TEST_FIXTURE(NAME_OF_FIXTURE,NAME_OF_TEST){ body_of_test }

    == CHECK MACRO REFERENCE ==
         - CHECK(EXPR);
         - CHECK_EQUAL(EXPECTED,ACTUAL);
         - CHECK_CLOSE(EXPECTED,ACTUAL,EPSILON);
         - CHECK_ARRAY_EQUAL(EXPECTED,ACTUAL,LENGTH);
         - CHECK_ARRAY_CLOSE(EXPECTED,ACTUAL,LENGTH,EPSILON);
         - CHECK_ARRAY2D_EQUAL(EXPECTED,ACTUAL,ROWCOUNT,COLCOUNT);
         - CHECK_ARRAY2D_CLOSE(EXPECTED,ACTUAL,ROWCOUNT,COLCOUNT,EPSILON);
         - CHECK_THROW(EXPR,EXCEPTION_TYPE_EXPECTED);

    == TIME CONSTRAINTS ==
         - UNITTEST_TIME_CONSTRAINT(TIME_IN_MILLISECONDS);
         - UNITTEST_TIME_CONSTRAINT_EXEMPT();

    == MORE INFO ==
         See: https://unittest-cpp.github.io/

    Author(s): A. Mantzaflaris,  H. Weiner
 **/

#include "gismo_unittest.h"       // Brings in G+Smo and the UnitTest++ framework

SUITE(gsGeometry_test)                 // The suite should have the same name as the file
{
    gsGeometry<>::uPtr g = gsReadFile<>("surfaces/simple.xml");   
    gsMatrix<> uv = gsPointGrid<>(g->support(), 5);
    gsMatrix<> xyz = g->eval(uv);
    
    TEST(recoverPoints)
    {
        gsMatrix<> puv, pxy = xyz;
        g->recoverPoints(pxy, puv, 2, 1e-8);
        CHECK(  (puv-uv ).norm() < 1e-6 );
        CHECK(  (pxy-xyz).norm() < 1e-6 );
    }

    TEST(rotate2D)
    {
        gsGeometry<>::uPtr g = gsNurbsCreator<>::BSplineSquare();
        gsMatrix<> orig = g->coefs();
        g->rotate(EIGEN_PI / 2);
        // 90-degree CCW rotation around origin: (x, y) -> (-y, x)
        for (index_t i = 0; i < orig.rows(); ++i)
        {
            CHECK_CLOSE(g->coefs()(i, 0), -orig(i, 1), EPSILON);
            CHECK_CLOSE(g->coefs()(i, 1),  orig(i, 0), EPSILON);
        }
    }

    TEST(rotate3D)
    {
        gsGeometry<>::uPtr g = gsNurbsCreator<>::BSplineCube();
        gsMatrix<> orig = g->coefs();
        gsVector<real_t, 3> axis;
        axis << 1, 0, 0;
        g->rotate(EIGEN_PI / 2, axis);
        // 90-degree rotation around X-axis: (x, y, z) -> (x, -z, y)
        for (index_t i = 0; i < orig.rows(); ++i)
        {
            CHECK_CLOSE(g->coefs()(i, 0),  orig(i, 0), EPSILON);
            CHECK_CLOSE(g->coefs()(i, 1), -orig(i, 2), EPSILON);
            CHECK_CLOSE(g->coefs()(i, 2),  orig(i, 1), EPSILON);
        }
    }


    // Quarter circle of radius 2 from (2,0) to (0,2), a curve in 2D whose
    // speed |G'| varies along it and is not 1 at the end points.
    static gsNurbs<real_t> quarterCircle()
    {
        const real_t w = 1 / math::sqrt(real_t(2));
        gsMatrix<real_t> weights(3, 1), cp(3, 2);
        weights << 1, w, 1;
        cp << 2, 0,  2, 2,  0, 2;
        return gsNurbs<real_t>(gsKnotVector<real_t>(0, 1, 0, 3), weights, cp);
    }

    static gsMapData<real_t> endPointData(const gsGeometry<real_t> & geo, boxSide side)
    {
        gsMapData<real_t> md(NEED_VALUE | NEED_MEASURE | NEED_OUTER_NORMAL | NEED_GRAD_TRANSFORM);
        md.side   = side;
        md.points = geo.support().col(boundary::west == side ? 0 : 1);
        geo.computeMap(md);
        return md;
    }

    // The boundary of a curve is a point: measure 1, and the outer normal is
    // the unit tangent, pointing out of the parameter interval.
    TEST(curveEndPointMapData)
    {
        const real_t tol = 100 * math::limits::epsilon();
        const boxSide sides[2] = {boundary::west, boundary::east};

        // curve in 2D: outward is -y at (2,0) and -x at (0,2)
        const gsNurbs<real_t> arc = quarterCircle();
        const real_t arcNormal[2][2] = {{0, -1}, {-1, 0}};
        for (index_t i = 0; i != 2; ++i)
        {
            const gsMapData<real_t> md = endPointData(arc, sides[i]);
            CHECK_EQUAL(1, md.measures.size());
            CHECK_CLOSE(1, md.measures(0, 0), tol);
            CHECK_EQUAL(2, md.outNormals.rows());
            CHECK_CLOSE(arcNormal[i][0], md.outNormals(0, 0), tol);
            CHECK_CLOSE(arcNormal[i][1], md.outNormals(1, 0), tol);
        }

        // curve in 3D: the line (t, 2t, 2t), speed 3
        gsMatrix<real_t> cp3(2, 3);
        cp3 << 0, 0, 0,  1, 2, 2;
        const gsBSpline<real_t> line(gsKnotVector<real_t>(0, 1, 0, 2), cp3);
        for (index_t i = 0; i != 2; ++i)
        {
            const gsMapData<real_t> md = endPointData(line, sides[i]);
            const real_t sgn = (0 == i ? -1 : 1);
            CHECK_CLOSE(1, md.measures(0, 0), tol);
            CHECK_EQUAL(3, md.outNormals.rows());
            CHECK_CLOSE(sgn / 3,     md.outNormals(0, 0), tol);
            CHECK_CLOSE(sgn * 2 / 3, md.outNormals(1, 0), tol);
            CHECK_CLOSE(sgn * 2 / 3, md.outNormals(2, 0), tol);
        }

        // 1D -> 1D, orientation reversing: x = 2 - 2t, so the outer normal
        // is +1 at t = 0 (x = 2) and -1 at t = 1 (x = 0); the other
        // requested data are filled as well
        gsMatrix<real_t> cp1(2, 1);
        cp1 << 2, 0;
        const gsBSpline<real_t> seg(gsKnotVector<real_t>(0, 1, 0, 2), cp1);
        for (index_t i = 0; i != 2; ++i)
        {
            const gsMapData<real_t> md = endPointData(seg, sides[i]);
            CHECK_EQUAL(1, md.measures.size());
            CHECK_CLOSE(1, md.measures(0, 0), tol);
            CHECK_CLOSE(0 == i ? 1 : -1, md.outNormals(0, 0), tol);
            CHECK_EQUAL(1, md.jacInvTr.size());
            CHECK_CLOSE(-0.5, md.jacInvTr(0, 0), tol);
        }

        // a vanishing tangent at the end point has no direction: refused
        gsMatrix<real_t> cpd(3, 2);
        cpd << 0, 0,  0, 0,  1, 0;
        const gsBSpline<real_t> cusp(gsKnotVector<real_t>(0, 1, 0, 3), cpd);
        CHECK_THROW(endPointData(cusp, boundary::west), std::runtime_error);
        CHECK_EQUAL(1, endPointData(cusp, boundary::east).measures.size());
    }

    // Where the derivative of a curve vanishes at its end point, the measure
    // there is still 1, and a 1D -> 1D map still has an outer normal: the
    // orientation of the map.
    TEST(curveEndPointVanishingTangent)
    {
        const real_t tol = 100 * math::limits::epsilon();
        const auto westData = [](const gsGeometry<real_t> & geo, unsigned flags)
        {
            gsMapData<real_t> md(NEED_VALUE | flags);
            md.side   = boundary::west;
            md.points = geo.support().col(0);
            geo.computeMap(md);
            return md;
        };

        // x = t^2 (orientation preserving) and x = 1 - t^2 (reversing)
        const real_t outward[2] = {-1, 1};
        for (index_t i = 0; i != 2; ++i)
        {
            gsMatrix<real_t> cp(3, 1);
            if (0 == i) cp << 0, 0, 1;
            else        cp << 1, 1, 0;
            const gsBSpline<real_t> line(gsKnotVector<real_t>(0, 1, 0, 3), cp);
            const gsMapData<real_t> md =
                westData(line, NEED_MEASURE | NEED_OUTER_NORMAL);
            CHECK_CLOSE(1, md.measures(0, 0), tol);
            CHECK_EQUAL(1, md.outNormals.size());
            CHECK_CLOSE(outward[i], md.outNormals(0, 0), tol);
        }

        // a curve in 2D with the same cusp: its measure alone is available
        gsMatrix<real_t> cpd(3, 2);
        cpd << 0, 0,  0, 0,  1, 0;
        const gsBSpline<real_t> cusp(gsKnotVector<real_t>(0, 1, 0, 3), cpd);
        CHECK_CLOSE(1, westData(cusp, NEED_MEASURE).measures(0, 0), tol);
    }

    // A point Neumann term at the end of a curve is g times the end basis
    // function, whichever spelling of the boundary measure is used.
    TEST(curveEndPointNeumann)
    {
        const real_t tol = 100 * math::limits::epsilon();
        gsMultiPatch<real_t> mp(quarterCircle());
        mp.uniformRefine();
        mp.uniformRefine();
        gsMultiBasis<real_t> mb(mp);
        const index_t n = mb.basis(0).size();
        gsFunctionExpr<real_t> g("sin(x) + y", 2);
        const real_t gEnd[2] = {math::sin(real_t(2)), 2}; // at (2,0) and (0,2)

        const boxSide sides[2] = {boundary::west, boundary::east};
        for (index_t i = 0; i != 2; ++i)
        {
            gsBoundaryConditions<real_t> bc;
            bc.addCondition(0, sides[i], condition_type::neumann, &g);

            gsExprAssembler<real_t> A(1,1);
            A.setIntegrationElements(mb);
            auto G = A.getMap(mp);
            auto u = A.getSpace(mb);
            auto h = A.getBdrFunction(G);
            u.setup(gsBoundaryConditions<real_t>(), dirichlet::homogeneous, 0);

            gsMatrix<real_t> expected = gsMatrix<real_t>::Zero(n, 1);
            expected(u.mapper().index(0 == i ? 0 : n-1, 0), 0) = gEnd[i];

            A.initSystem();
            A.assembleBdr(bc.get("Neumann"), u * h * meas(G));
            CHECK((A.rhs() - expected).norm() < tol);

            A.initSystem();
            A.assembleBdr(bc.get("Neumann"), u * h * nv(G).norm());
            CHECK((A.rhs() - expected).norm() < tol);
        }

        // the boundary of a curve is its two end points
        gsExprEvaluator<real_t> ev;
        ev.setIntegrationElements(mb);
        auto G = ev.getMap(mp);
        CHECK_CLOSE(2, ev.integralBdr(meas(G)), tol);
        CHECK_CLOSE(2, ev.integralBdr(nv(G).norm()), tol);
    }

}
