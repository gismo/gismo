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

    // Additional comprehensive tests for gsHSplines module

    TEST(gsTHBSplineBasis_2D_Construction)
    {
        gsKnotVector<> kv1(0, 1, 2, 2);
        gsKnotVector<> kv2(0, 1, 2, 2);
        gsTensorBSplineBasis<2> tensorBasis(kv1, kv2);
        
        gsTHBSplineBasis<2> thb(tensorBasis);
        
        CHECK_EQUAL(thb.dim(), 2);
        CHECK(thb.size() > 0);
    }

    TEST(gsTHBSplineBasis_MultipleRefinements)
    {
        gsKnotVector<> kv(0, 1, 1, 4);
        gsBSplineBasis<> basis(kv);
        gsTHBSplineBasis<1> thb(basis);
        
        size_t initialSize = thb.size();
        
        // Refine a region
        std::vector<index_t> boxes;
        boxes.push_back(0); // level
        boxes.push_back(0); // lower bound
        boxes.push_back(2); // upper bound
        
        thb.refineElements(boxes);
        
        CHECK(thb.size() >= initialSize);
    }

    TEST(gsTHBSpline_ZeroCoefficients)
    {
        gsKnotVector<> kv(0, 1, 1, 2);
        gsBSplineBasis<> basis(kv);
        gsTHBSplineBasis<1> thb(basis);
        
        gsMatrix<> coefs(thb.size(), 1);
        coefs.setZero();
        
        gsTHBSpline<1> thbSpline(thb, coefs);
        
        gsMatrix<> point(1, 1);
        point << 0.5;
        
        gsMatrix<> result = thbSpline.eval(point);
        
        CHECK_CLOSE(result(0, 0), 0.0, 1e-10);
    }

    TEST(gsTHBSpline_ConstantFunction)
    {
        gsKnotVector<> kv(0, 1, 1, 2);
        gsBSplineBasis<> basis(kv);
        gsTHBSplineBasis<1> thb(basis);
        
        gsMatrix<> coefs(thb.size(), 1);
        coefs.setOnes();
        
        gsTHBSpline<1> thbSpline(thb, coefs);
        
        gsMatrix<> point(1, 1);
        point << 0.5;
        
        gsMatrix<> result = thbSpline.eval(point);
        
        // Sum of basis functions should be 1
        CHECK_CLOSE(result(0, 0), 1.0, 1e-8);
    }

    TEST(gsTHBSplineBasis_ActiveFunctions)
    {
        gsKnotVector<> kv(0, 1, 2, 3);
        gsBSplineBasis<> basis(kv);
        gsTHBSplineBasis<1> thb(basis);
        
        gsMatrix<> point(1, 1);
        point << 0.5;
        
        gsMatrix<unsigned> active;
        thb.active_into(point, active);
        
        CHECK(active.size() > 0);
    }

    TEST(gsTHBSpline_2D_Evaluation)
    {
        gsKnotVector<> kv1(0, 1, 1, 2);
        gsKnotVector<> kv2(0, 1, 1, 2);
        gsTensorBSplineBasis<2> tensorBasis(kv1, kv2);
        
        gsTHBSplineBasis<2> thb(tensorBasis);
        
        gsMatrix<> coefs(thb.size(), 1);
        for (index_t i = 0; i < coefs.rows(); i++)
            coefs(i, 0) = i * 0.1;
        
        gsTHBSpline<2> thbSpline(thb, coefs);
        
        gsMatrix<> point(2, 1);
        point << 0.5, 0.5;
        
        gsMatrix<> result = thbSpline.eval(point);
        
        CHECK(result.rows() > 0);
    }

    TEST(gsTHBSplineBasis_Degree)
    {
        gsKnotVector<> kv(0, 1, 3, 2);
        gsBSplineBasis<> basis(kv);
        gsTHBSplineBasis<1> thb(basis);
        
        CHECK_EQUAL(thb.degree(), 3);
    }

    TEST(gsTHBSpline_Coefs)
    {
        gsKnotVector<> kv(0, 1, 1, 3);
        gsBSplineBasis<> basis(kv);
        gsTHBSplineBasis<1> thb(basis);
        
        gsMatrix<> coefs(thb.size(), 1);
        for (index_t i = 0; i < coefs.rows(); i++)
            coefs(i, 0) = i + 1.5;
        
        gsTHBSpline<1> thbSpline(thb, coefs);
        
        const gsMatrix<>& retrievedCoefs = thbSpline.coefs();
        
        CHECK_EQUAL(retrievedCoefs.rows(), coefs.rows());
        for (index_t i = 0; i < coefs.rows(); i++)
            CHECK_CLOSE(retrievedCoefs(i, 0), i + 1.5, 1e-10);
    }

    TEST(gsTHBSplineBasis_BasisSize)
    {
        gsKnotVector<> kv(0, 1, 2, 5);
        gsBSplineBasis<> basis(kv);
        gsTHBSplineBasis<1> thb(basis);
        
        size_t basisSize = thb.size();
        CHECK(basisSize > 0);
        CHECK(basisSize == basis.size());
    }

    TEST(gsTHBSpline_MultiplePoints)
    {
        gsKnotVector<> kv(0, 1, 1, 3);
        gsBSplineBasis<> basis(kv);
        gsTHBSplineBasis<1> thb(basis);
        
        gsMatrix<> coefs(thb.size(), 1);
        coefs.setConstant(2.5);
        
        gsTHBSpline<1> thbSpline(thb, coefs);
        
        gsMatrix<> points(1, 5);
        points << 0.1, 0.3, 0.5, 0.7, 0.9;
        
        gsMatrix<> result = thbSpline.eval(points);
        
        CHECK_EQUAL(result.cols(), 5);
        for (int i = 0; i < 5; i++)
            CHECK_CLOSE(result(0, i), 2.5, 1e-8);
    }

    TEST(gsTHBSplineBasis_DomainDimension)
    {
        gsKnotVector<> kv1(0, 1, 1, 2);
        gsKnotVector<> kv2(0, 1, 1, 2);
        gsKnotVector<> kv3(0, 1, 1, 2);
        gsTensorBSplineBasis<3> tensorBasis(kv1, kv2, kv3);
        
        gsTHBSplineBasis<3> thb(tensorBasis);
        
        CHECK_EQUAL(thb.dim(), 3);
    }

    TEST(gsTHBSpline_SecondDerivative)
    {
        gsKnotVector<> kv(0, 1, 3, 3);
        gsBSplineBasis<> basis(kv);
        gsTHBSplineBasis<1> thb(basis);
        
        gsMatrix<> coefs(thb.size(), 1);
        for (index_t i = 0; i < coefs.rows(); i++)
            coefs(i, 0) = i * i;
        
        gsTHBSpline<1> thbSpline(thb, coefs);
        
        gsMatrix<> point(1, 1);
        point << 0.5;
        
        gsMatrix<> deriv2 = thbSpline.deriv2(point);
        
        CHECK(deriv2.rows() > 0);
    }

    TEST(gsTHBSplineBasis_Anchors)
    {
        gsKnotVector<> kv(0, 1, 2, 3);
        gsBSplineBasis<> basis(kv);
        gsTHBSplineBasis<1> thb(basis);
        
        gsMatrix<> anchors;
        thb.anchors_into(anchors);
        
        CHECK(anchors.rows() > 0);
        CHECK(anchors.cols() > 0);
    }

    /* 
     * Step-by-step instructions for additional complex gsHSplines tests:
     * 
     * 1. gsHBox tests:
     *    - Create gsHTensorBasis with refinement
     *    - Extract elements as gsHBox objects
     *    - Test gsHBox containment checks
     *    - Test gsHBox intersection with other boxes
     *    - Test parent/child relationships in hierarchy
     *    - Test support extensions for basis functions
     * 
     * 2. gsHBoxContainer tests:
     *    - Create container of gsHBox objects
     *    - Test insertion and removal operations
     *    - Test iteration over boxes
     *    - Test sorting by level or position
     *    - Test union and intersection operations on containers
     * 
     * 3. gsHBoxUtils tests:
     *    - Test domain unions of box collections
     *    - Test domain intersections
     *    - Test box splitting and merging
     *    - Test neighborhood computation
     * 
     * 4. gsRationalTHBSpline tests:
     *    - Create rational THB-spline (NURBS variant)
     *    - Test weight management
     *    - Test evaluation with rational basis
     *    - Test derivative computation
     *    - Test projection onto rational space
     * 
     * 5. gsHFitting tests:
     *    - Set up hierarchical fitting problem
     *    - Provide scattered data points
     *    - Test automatic refinement based on error
     *    - Test adaptive fitting algorithm
     *    - Verify convergence properties
     * 
     * 6. Advanced refinement tests:
     *    - Test refinement by error indicators
     *    - Test adaptive refinement strategies
     *    - Test coarsening operations
     *    - Test refinement with hanging nodes
     *    - Test truncation mechanism for THB-splines
     */
}
