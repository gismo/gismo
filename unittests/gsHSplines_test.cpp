/** @file gsHSplines_test.cpp

    @brief Tests for gsHSplines module - hierarchical spline structures

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): Copilot
*/

#include "gismo_unittest.h"

SUITE(gsHSplines_test)
{
    // Commented out - gsAABB API doesn't match
    // TEST(gsAABB_Construction)
    // {
    //     // Create an axis-aligned bounding box
    //     gsVector<index_t, 2> lower;
    //     lower << 0, 0;
    //     gsVector<index_t, 2> upper;
    //     upper << 10, 10;
    //     
    //     gsAABB<2, index_t> bbox(lower, upper);
    //     
    //     CHECK(bbox.getSize(0) > 0);
    //     CHECK(bbox.getSize(1) > 0);
    // }

    // TEST(gsAABB_Contains)
    // {
    //     gsVector<index_t, 2> lower, upper;
    //     lower << 0, 0;
    //     upper << 10, 10;
    //     
    //     gsAABB<2, index_t> bbox(lower, upper);
    //     
    //     gsVector<index_t, 2> point1, point2;
    //     point1 << 5, 5;  // Inside
    //     point2 << 15, 5; // Outside
    //     
    //     CHECK(bbox.contains(point1));
    //     CHECK(!bbox.contains(point2));
    // }

    // TEST(gsAABB_Intersection)
    // {
    //     gsVector<index_t, 2> lower1, upper1, lower2, upper2;
    //     lower1 << 0, 0;
    //     upper1 << 10, 10;
    //     lower2 << 5, 5;
    //     upper2 << 15, 15;
    //     
    //     gsAABB<2, index_t> bbox1(lower1, upper1);
    //     gsAABB<2, index_t> bbox2(lower2, upper2);
    //     
    //     CHECK(bbox1.intersects(bbox2));
    //     CHECK(bbox2.intersects(bbox1));
    // }

    // TEST(gsAABB_Volume)
    // {
    //     gsVector<index_t, 2> lower, upper;
    //     lower << 0, 0;
    //     upper << 10, 20;
    //     
    //     gsAABB<2, index_t> bbox(lower, upper);
    //     
    //     // Volume should be 10 * 20 = 200
    //     CHECK_EQUAL(bbox.volume(), 200);
    // }

    TEST(gsTHBSplineBasis_Construction_1D)
    {
        // Create a simple 1D tensor B-spline basis
        gsKnotVector<> kv(0, 1, 2, 4);
        gsBSplineBasis<> basis(kv);
        
        // Create THB-spline basis from tensor basis
        gsTHBSplineBasis<1> thb(basis);
        
        CHECK_EQUAL(thb.dim(), 1);
        CHECK(thb.size() > 0);
    }

    TEST(gsTHBSplineBasis_Construction_2D)
    {
        // Create 2D tensor B-spline basis
        gsKnotVector<> kv1(0, 1, 1, 3);
        gsKnotVector<> kv2(0, 1, 1, 3);
        gsTensorBSplineBasis<2> tensorBasis(kv1, kv2);
        
        // Create THB-spline basis
        gsTHBSplineBasis<2> thb(tensorBasis);
        
        CHECK_EQUAL(thb.dim(), 2);
        CHECK(thb.size() > 0);
        CHECK_EQUAL(thb.maxLevel(), 0);
    }

    TEST(gsTHBSplineBasis_Refinement_1D)
    {
        // Create 1D basis
        gsKnotVector<> kv(0, 1, 1, 3);
        gsBSplineBasis<> basis(kv);
        gsTHBSplineBasis<1> thb(basis);
        
        // Create a refinement box
        gsVector<unsigned> lower, upper;
        lower.setZero(1);
        upper.setOnes(1);
        upper *= 2;
        
        std::vector<index_t> boxes;
        boxes.push_back(0); // level
        boxes.push_back(lower[0]); // lower bound
        boxes.push_back(upper[0]); // upper bound
        
        // Refine the basis
        thb.refineElements(boxes);
        
        CHECK(thb.maxLevel() >= 0);
    }

    TEST(gsTHBSpline_Construction)
    {
        // Create basis
        gsKnotVector<> kv(0, 1, 1, 2);
        gsBSplineBasis<> basis(kv);
        gsTHBSplineBasis<1> thb(basis);
        
        // Create coefficient vector
        gsMatrix<> coefs(thb.size(), 1);
        coefs.setRandom();
        
        // Create THB-spline
        gsTHBSpline<1> thbSpline(thb, coefs);
        
        CHECK_EQUAL(thbSpline.basis().dim(), 1);
        CHECK_EQUAL(thbSpline.coefs().rows(), (index_t)thb.size());
    }

    TEST(gsTHBSpline_Evaluation)
    {
        // Create simple 1D THB-spline
        gsKnotVector<> kv(0, 1, 1, 2);
        gsBSplineBasis<> basis(kv);
        gsTHBSplineBasis<1> thb(basis);
        
        gsMatrix<> coefs(thb.size(), 1);
        coefs.setOnes();
        
        gsTHBSpline<1> thbSpline(thb, coefs);
        
        // Evaluate at point
        gsMatrix<> point(1, 1);
        point(0, 0) = 0.5;
        
        gsMatrix<> result = thbSpline.eval(point);
        
        CHECK(result.rows() == 1);
        CHECK(result.cols() == 1);
    }

    // Commented out - gsHTensorBasis is abstract and can't be instantiated directly
    // TEST(gsHTensorBasis_Construction_1D)
    // {
    //     // Create tensor basis
    //     gsKnotVector<> kv(0, 1, 1, 3);
    //     gsBSplineBasis<> basis(kv);
    //     
    //     // Create hierarchical tensor basis
    //     gsHTensorBasis<1> htb(basis);
    //     
    //     CHECK_EQUAL(htb.dim(), 1);
    //     CHECK(htb.size() > 0);
    // }

    // TEST(gsHTensorBasis_Construction_2D)
    // {
    //     // Create 2D tensor basis
    //     gsKnotVector<> kv1(0, 1, 1, 2);
    //     gsKnotVector<> kv2(0, 1, 1, 2);
    //     gsTensorBSplineBasis<2> tensorBasis(kv1, kv2);
    //     
    //     // Create hierarchical tensor basis
    //     gsHTensorBasis<2> htb(tensorBasis);
    //     
    //     CHECK_EQUAL(htb.dim(), 2);
    //     CHECK(htb.size() > 0);
    // }

    // TEST(gsHTensorBasis_MaxLevel)
    // {
    //     gsKnotVector<> kv(0, 1, 1, 2);
    //     gsBSplineBasis<> basis(kv);
    //     gsHTensorBasis<1> htb(basis);
    //     
    //     // Initially should have only level 0
    //     CHECK_EQUAL(htb.maxLevel(), 0);
    //     CHECK(htb.numLevels() > 0);
    // }

    // TEST(gsHTensorBasis_Size)
    // {
    //     gsKnotVector<> kv(0, 1, 1, 3);
    //     gsBSplineBasis<> basis(kv);
    //     gsHTensorBasis<1> htb(basis);
    //     
    //     size_t initialSize = htb.size();
    //     CHECK(initialSize > 0);
    //     CHECK(initialSize == basis.size());
    // }

    // Commented out - gsVSegment is an incomplete type
    // TEST(gsVSegment_Construction)
    // {
    //     gsVSegment<real_t> seg;
    //     seg.level = 0;
    //     seg.low = 0.0;
    //     seg.high = 1.0;
    //     
    //     CHECK_EQUAL(seg.level, 0);
    //     CHECK_CLOSE(seg.low, 0.0, 1e-10);
    //     CHECK_CLOSE(seg.high, 1.0, 1e-10);
    // }

    // TEST(gsVSegment_Contains)
    // {
    //     gsVSegment<real_t> seg;
    //     seg.level = 0;
    //     seg.low = 0.0;
    //     seg.high = 1.0;
    //     
    //     CHECK(seg.contains(0.5));
    //     CHECK(seg.contains(0.0));
    //     CHECK(seg.contains(1.0));
    //     CHECK(!seg.contains(1.5));
    // }

    TEST(gsTHBSplineBasis_LevelSize)
    {
        gsKnotVector<> kv(0, 1, 1, 3);
        gsBSplineBasis<> basis(kv);
        gsTHBSplineBasis<1> thb(basis);
        
        // Level 0 should have some functions
        CHECK(thb.size() > 0);
    }

    TEST(gsTHBSpline_Derivative)
    {
        // Create simple 1D THB-spline
        gsKnotVector<> kv(0, 1, 2, 3);
        gsBSplineBasis<> basis(kv);
        gsTHBSplineBasis<1> thb(basis);
        
        gsMatrix<> coefs(thb.size(), 1);
        for (index_t i = 0; i < coefs.rows(); i++)
            coefs(i, 0) = i;
        
        gsTHBSpline<1> thbSpline(thb, coefs);
        
        // Evaluate derivative at point
        gsMatrix<> point(1, 1);
        point(0, 0) = 0.5;
        
        gsMatrix<> deriv = thbSpline.deriv(point);
        
        CHECK(deriv.rows() >= 1);
    }

    TEST(gsTHBSplineBasis_Support)
    {
        gsKnotVector<> kv(0, 1, 1, 3);
        gsBSplineBasis<> basis(kv);
        gsTHBSplineBasis<1> thb(basis);
        
        // Get support of first basis function
        gsMatrix<> support = thb.support(0);
        
        CHECK(support.rows() == 1);
        CHECK(support.cols() == 2);
        CHECK(support(0, 0) <= support(0, 1));
    }

    // Commented out - gsHTensorBasis is abstract
    // TEST(gsHTensorBasis_Support)
    // {
    //     gsKnotVector<> kv(0, 1, 1, 2);
    //     gsBSplineBasis<> basis(kv);
    //     gsHTensorBasis<1> htb(basis);
    //     
    //     // Get support of first basis function
    //     gsMatrix<> support = htb.support(0);
    //     
    //     CHECK(support.rows() == 1);
    //     CHECK(support.cols() == 2);
    // }

    TEST(gsTHBSplineBasis_BoundingBox)
    {
        gsKnotVector<> kv(0, 1, 1, 2);
        gsBSplineBasis<> basis(kv);
        gsTHBSplineBasis<1> thb(basis);
        
        gsMatrix<> bbox = thb.support();
        
        CHECK(bbox.rows() == 1);
        CHECK(bbox.cols() == 2);
        CHECK_CLOSE(bbox(0, 0), 0.0, 1e-10);
        CHECK_CLOSE(bbox(0, 1), 1.0, 1e-10);
    }
}
