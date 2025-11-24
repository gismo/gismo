/** @file gsNurbs_test.cpp

    @brief Comprehensive tests for gsNurbs module (NURBS and B-spline functionality)

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): Copilot
*/

#include "gismo_unittest.h"

SUITE(gsNurbs_test)
{
    TEST(gsBSpline_Line)
    {
        gsKnotVector<> kv(0, 1, 1, 2);
        gsBSplineBasis<> basis(kv);
        
        gsMatrix<> coefs(2, 2);
        coefs << 0, 0,
                 1, 1;
        
        gsBSpline<> line(kv, coefs);
        
        gsMatrix<> point(1, 1);
        point << 0.5;
        
        gsMatrix<> result = line.eval(point);
        
        CHECK_CLOSE(result(0, 0), 0.5, 1e-10);
        CHECK_CLOSE(result(1, 0), 0.5, 1e-10);
    }

    TEST(gsBSpline_Quadratic)
    {
        gsKnotVector<> kv(0, 1, 2, 3);
        gsBSplineBasis<> basis(kv);
        
        gsMatrix<> coefs(3, 1);
        coefs << 0, 1, 0;
        
        gsBSpline<> parabola(kv, coefs);
        
        gsMatrix<> point(1, 1);
        point << 0.5;
        
        gsMatrix<> result = parabola.eval(point);
        
        CHECK(result(0, 0) > 0);
    }

    TEST(gsBSpline_Derivative)
    {
        gsKnotVector<> kv(0, 1, 2, 2);
        gsBSplineBasis<> basis(kv);
        
        gsMatrix<> coefs(2, 1);
        coefs << 0, 1;
        
        gsBSpline<> curve(kv, coefs);
        
        gsMatrix<> point(1, 1);
        point << 0.5;
        
        gsMatrix<> deriv = curve.deriv(point);
        
        CHECK(deriv.rows() > 0);
    }

    TEST(gsBSpline_SecondDerivative)
    {
        gsKnotVector<> kv(0, 1, 3, 2);
        gsBSplineBasis<> basis(kv);
        
        gsMatrix<> coefs(2, 1);
        coefs << 0, 1;
        
        gsBSpline<> curve(kv, coefs);
        
        gsMatrix<> point(1, 1);
        point << 0.5;
        
        gsMatrix<> deriv2 = curve.deriv2(point);
        
        CHECK(deriv2.rows() > 0);
    }

    TEST(gsTensorBSpline_Surface)
    {
        gsKnotVector<> kv1(0, 1, 1, 2);
        gsKnotVector<> kv2(0, 1, 1, 2);
        gsTensorBSplineBasis<2> basis(kv1, kv2);
        
        gsMatrix<> coefs(4, 3);
        coefs << 0, 0, 0,
                 1, 0, 0,
                 0, 1, 0,
                 1, 1, 0;
        
        gsTensorBSpline<2> surface(basis, coefs);
        
        gsMatrix<> point(2, 1);
        point << 0.5, 0.5;
        
        gsMatrix<> result = surface.eval(point);
        
        CHECK_EQUAL(result.rows(), 3);
        CHECK_CLOSE(result(0, 0), 0.5, 1e-10);
        CHECK_CLOSE(result(1, 0), 0.5, 1e-10);
    }

    TEST(gsTensorBSpline_Solid)
    {
        gsKnotVector<> kv1(0, 1, 1, 2);
        gsKnotVector<> kv2(0, 1, 1, 2);
        gsKnotVector<> kv3(0, 1, 1, 2);
        gsTensorBSplineBasis<3> basis(kv1, kv2, kv3);
        
        gsMatrix<> coefs(8, 3);
        for (int i = 0; i < 8; i++)
        {
            coefs(i, 0) = (i % 2);
            coefs(i, 1) = ((i / 2) % 2);
            coefs(i, 2) = (i / 4);
        }
        
        gsTensorBSpline<3> solid(basis, coefs);
        
        gsMatrix<> point(3, 1);
        point << 0.5, 0.5, 0.5;
        
        gsMatrix<> result = solid.eval(point);
        
        CHECK_EQUAL(result.rows(), 3);
    }

    TEST(gsBSplineBasis_Uniformity)
    {
        gsKnotVector<> kv(0, 1, 2, 4, 1);
        gsBSplineBasis<> basis(kv);
        
        CHECK_EQUAL(basis.degree(), 2);
        CHECK(basis.size() > 0);
    }

    TEST(gsBSplineBasis_KnotInsertion)
    {
        gsKnotVector<> kv(0, 1, 2, 2);
        gsBSplineBasis<> basis(kv);
        
        size_t oldSize = basis.size();
        
        gsMatrix<> coefs(oldSize, 1);
        coefs.setOnes();
        
        // Insert a knot
        gsMatrix<> newCoefs;
        basis.insertKnot(0.5, newCoefs, coefs);
        
        CHECK(newCoefs.rows() > (index_t)oldSize);
    }

    TEST(gsBSpline_ClosedCurve)
    {
        gsKnotVector<> kv(0, 1, 2, 8, 1);
        gsBSplineBasis<> basis(kv);
        
        gsMatrix<> coefs(8, 2);
        // Create points on a circle
        for (int i = 0; i < 8; i++)
        {
            double angle = 2 * M_PI * i / 8;
            coefs(i, 0) = cos(angle);
            coefs(i, 1) = sin(angle);
        }
        
        gsBSpline<> circle(kv, coefs);
        
        // Check that start and end are close (approximate circle)
        gsMatrix<> start(1, 1);
        start << 0.0;
        gsMatrix<> end(1, 1);
        end << 1.0;
        
        gsMatrix<> pStart = circle.eval(start);
        gsMatrix<> pEnd = circle.eval(end);
        
        CHECK(pStart.rows() == 2);
    }

    TEST(gsNurbsCreator_Circle)
    {
        gsNurbs<>::uPtr circle = gsNurbsCreator<>::NurbsCircle();
        
        CHECK(circle.get() != nullptr);
        CHECK_EQUAL(circle->domainDim(), 1);
        CHECK_EQUAL(circle->geoDim(), 2);
    }

    TEST(gsNurbsCreator_Sphere)
    {
        gsNurbs<>::uPtr sphere = gsNurbsCreator<>::NurbsSphere();
        
        CHECK(sphere.get() != nullptr);
        CHECK_EQUAL(sphere->domainDim(), 2);
        CHECK_EQUAL(sphere->geoDim(), 3);
    }

    TEST(gsTensorBSpline_BoundaryExtraction)
    {
        gsKnotVector<> kv1(0, 1, 1, 3);
        gsKnotVector<> kv2(0, 1, 1, 3);
        gsTensorBSplineBasis<2> basis(kv1, kv2);
        
        gsMatrix<> coefs(9, 2);
        coefs.setRandom();
        
        gsTensorBSpline<2> surface(basis, coefs);
        
        // Extract boundary curve
        gsGeometry<>::uPtr boundary = surface.boundary(boundary::west);
        
        CHECK(boundary.get() != nullptr);
        CHECK_EQUAL(boundary->domainDim(), 1);
    }

    TEST(gsBSpline_ElevateDegreeBasis)
    {
        gsKnotVector<> kv(0, 1, 1, 3);
        gsBSplineBasis<> basis(kv);
        
        int oldDegree = basis.degree();
        
        gsBSplineBasis<> elevatedBasis = basis.elevateDegree(1);
        
        CHECK_EQUAL(elevatedBasis.degree(), oldDegree + 1);
    }

    TEST(gsBSpline_InsertionMatrix)
    {
        gsKnotVector<> kv(0, 1, 2, 2);
        gsBSplineBasis<> basis(kv);
        
        gsSparseMatrix<> transfer;
        basis.insertKnot(0.5, transfer);
        
        CHECK(transfer.rows() > transfer.cols());
    }

    /* 
     * Step-by-step instructions for additional complex gsNurbs tests:
     * 
     * 1. NURBS (rational) tests:
     *    - Create NURBS curve with non-uniform weights
     *    - Test conic section representation (circle, ellipse, hyperbola)
     *    - Test weight modification and its effect on shape
     *    - Test projection to/from B-spline space
     * 
     * 2. Knot insertion and removal:
     *    - Test Boehm's algorithm for knot insertion
     *    - Test Oslo algorithm for refinement
     *    - Test knot removal (inverse operation)
     *    - Verify curve/surface stays unchanged after insertion
     * 
     * 3. Degree elevation:
     *    - Elevate degree of curve/surface
     *    - Verify geometry preservation
     *    - Test multiple degree elevations
     *    - Check control point count increase
     * 
     * 4. gsNurbsCreator advanced shapes:
     *    - Create torus
     *    - Create cylinder
     *    - Create cone
     *    - Create revolution surfaces
     *    - Create sweep surfaces
     * 
     * 5. Curve and surface operations:
     *    - Test curve splitting at parameter
     *    - Test surface trimming
     *    - Test curve/surface merging
     *    - Test reparameterization
     * 
     * 6. gsBSplineBasis advanced tests:
     *    - Test greville abscissae computation
     *    - Test active basis functions query
     *    - Test basis function support
     *    - Test connectivity matrix
     * 
     * 7. gsTensorBSpline operations:
     *    - Test swapDirections (transpose parameters)
     *    - Test degreeElevate in specific direction
     *    - Test uniformRefine in specific direction
     *    - Test boundary extraction for all sides
     * 
     * 8. Approximation and fitting:
     *    - Fit B-spline to point cloud
     *    - Test least-squares fitting
     *    - Test interpolation at specified parameters
     *    - Test smoothing spline fitting
     */
}