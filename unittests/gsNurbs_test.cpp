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
    TEST(gsBSpline)
    {
        // Evaluation Test
        gsKnotVector<> kvLine(0, 1, 0, 2);
        gsInfo<<"kvLine: " << kvLine << "\n";
        gsBSplineBasis<> basisLine(kvLine);
        gsMatrix<> coefsLine(2, 2);
        coefsLine << 0, 0,
                     1, 1;
        gsBSpline<> line(kvLine, coefsLine);
        gsMatrix<> point(1, 1);
        point << 0.5;
        gsMatrix<> resultLine; line.eval_into(point, resultLine);
        CHECK_CLOSE(resultLine(0, 0), 0.5, 1e-10);
        CHECK_CLOSE(resultLine(1, 0), 0.5, 1e-10);
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
        auto sphere = gsNurbsCreator<>::NurbsSphere();
        
        CHECK(sphere.get() != nullptr);
        CHECK_EQUAL(sphere->domainDim(), 2);
        CHECK_EQUAL(sphere->geoDim(), 3);
    }

    TEST(gsTensorBSpline_BoundaryExtraction)
    {
        auto surface = gsNurbsCreator<>::BSplineSquare();

        // Extract boundary curve
        gsGeometry<>::uPtr boundary = surface->boundary(boundary::west);
        
        CHECK(boundary.get() != nullptr);
        CHECK_EQUAL(boundary->domainDim(), 1);
    }

    TEST(gsBSpline_ElevateDegreeBasis)
    {
        gsKnotVector<> kv(0, 1, 1, 3);
        gsBSplineBasis<> basis(kv);
        
        int oldDegree = basis.degree();
        basis.setDegree(oldDegree + 1);
        
        CHECK_EQUAL(basis.degree(), oldDegree + 1);
    }

    TEST(gsBSpline_InsertionMatrix)
    {
        gsKnotVector<> kv(0, 1, 2, 2);
        gsBSplineBasis<> basis(kv);
        
        gsSparseMatrix<double,RowMajor> transfer;
        std::vector<double> knots; knots.push_back(0.5);
        basis.refine_withTransfer(transfer, knots);
        
        CHECK(transfer.rows() > transfer.cols());
    }
}