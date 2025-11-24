/** @file gsDomain_test.cpp

    @brief Tests for gsDomain classes and domain iterators

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): Copilot
*/

#include "gismo_unittest.h"

SUITE(gsDomain_test)
{
    TEST(gsPointDomain_Construction)
    {
        // Create a matrix of test points
        gsMatrix<> points(2, 5);
        points << 0.0, 0.25, 0.5, 0.75, 1.0,
                  0.0, 0.25, 0.5, 0.75, 1.0;
        
        gsPointDomain<> pd(points);
        
        // Test basic properties
        CHECK_EQUAL(pd.dim(), 2);
        CHECK_EQUAL(pd.numElements(), 5);
        CHECK_EQUAL(pd.numElementsBdr(boundary::none), 0);
    }

    TEST(gsPointDomain_Iterator)
    {
        // Create a simple point domain
        gsMatrix<> points(2, 3);
        points << 0.0, 0.5, 1.0,
                  0.0, 0.5, 1.0;
        
        gsPointDomain<> pd(points);
        
        // Test iterator
        int count = 0;
        for (auto it = pd.beginAll(); it != pd.endAll(); ++it)
        {
            count++;
        }
        CHECK_EQUAL(count, 3);
    }

    TEST(gsPointDomain_Points)
    {
        gsMatrix<> points(3, 4);
        points << 1.0, 2.0, 3.0, 4.0,
                  5.0, 6.0, 7.0, 8.0,
                  9.0, 10.0, 11.0, 12.0;
        
        gsPointDomain<> pd(points);
        
        const gsMatrix<>& pts = pd.points();
        CHECK_EQUAL(pts.rows(), 3);
        CHECK_EQUAL(pts.cols(), 4);
        CHECK_CLOSE(pts(0, 1), 2.0, 1e-10);
    }

    TEST(gsTensorDomain_1D)
    {
        // Create a 1D knot vector for testing
        gsKnotVector<> kv(0, 1, 2, 4); // interval [0,1], degree 2, 4 elements
        
        std::vector<gsDomain<>::Ptr> kvs;
        kvs.push_back(memory::make_shared_not_owned(&kv));
        
        gsTensorDomain<real_t, 1> td(kvs);
        
        CHECK_EQUAL(td.dim(), 1);
        CHECK_EQUAL(td.degree(0), 2);
        CHECK(td.numElements() > 0);
    }

    TEST(gsTensorDomain_2D)
    {
        // Create 2D tensor domain
        gsKnotVector<> kv1(0, 1, 1, 3);
        gsKnotVector<> kv2(0, 1, 1, 3);
        
        std::vector<gsDomain<>::Ptr> kvs;
        kvs.push_back(memory::make_shared_not_owned(&kv1));
        kvs.push_back(memory::make_shared_not_owned(&kv2));
        
        gsTensorDomain<real_t, 2> td(kvs);
        
        CHECK_EQUAL(td.dim(), 2);
        CHECK_EQUAL(td.degree(0), 1);
        CHECK_EQUAL(td.degree(1), 1);
        
        size_t nElem = kv1.numElements() * kv2.numElements();
        CHECK_EQUAL(td.numElements(), nElem);
    }

    TEST(gsTensorDomain_BoundingBox)
    {
        gsKnotVector<> kv1(0, 1, 1, 2);
        gsKnotVector<> kv2(-1, 2, 1, 2);
        
        std::vector<gsDomain<>::Ptr> kvs;
        kvs.push_back(memory::make_shared_not_owned(&kv1));
        kvs.push_back(memory::make_shared_not_owned(&kv2));
        
        gsTensorDomain<real_t, 2> td(kvs);
        
        gsMatrix<> bbox = td.boundingBox();
        CHECK_EQUAL(bbox.rows(), 2);
        CHECK_EQUAL(bbox.cols(), 2);
        CHECK_CLOSE(bbox(0, 0), 0.0, 1e-10);
        CHECK_CLOSE(bbox(0, 1), 1.0, 1e-10);
        CHECK_CLOSE(bbox(1, 0), -1.0, 1e-10);
        CHECK_CLOSE(bbox(1, 1), 2.0, 1e-10);
    }

    TEST(gsTensorDomain_Iterator)
    {
        gsKnotVector<> kv(0, 1, 0, 5);
        
        std::vector<gsDomain<>::Ptr> kvs;
        kvs.push_back(memory::make_shared_not_owned(&kv));
        
        gsTensorDomain<real_t, 1> td(kvs);
        
        int count = 0;
        for (auto it = td.beginAll(); it != td.endAll(); ++it)
        {
            count++;
        }
        CHECK_EQUAL(count, (int)td.numElements());
    }

    TEST(gsTensorDomain_BoundaryElements)
    {
        gsKnotVector<> kv1(0, 1, 1, 3);
        gsKnotVector<> kv2(0, 1, 1, 4);
        
        std::vector<gsDomain<>::Ptr> kvs;
        kvs.push_back(memory::make_shared_not_owned(&kv1));
        kvs.push_back(memory::make_shared_not_owned(&kv2));
        
        gsTensorDomain<real_t, 2> td(kvs);
        
        // Boundary on side with direction 0 should have kv2.numElements() elements
        size_t nBdr0 = td.numElementsBdr(boxSide(1));
        CHECK_EQUAL(nBdr0, kv2.numElements());
        
        // Boundary on side with direction 1 should have kv1.numElements() elements
        size_t nBdr1 = td.numElementsBdr(boxSide(3));
        CHECK_EQUAL(nBdr1, kv1.numElements());
    }

    TEST(gsTensorDomain_ComponentAccess)
    {
        gsKnotVector<> kv1(0, 1, 2, 3);
        gsKnotVector<> kv2(0, 2, 1, 4);
        
        std::vector<gsDomain<>::Ptr> kvs;
        kvs.push_back(memory::make_shared_not_owned(&kv1));
        kvs.push_back(memory::make_shared_not_owned(&kv2));
        
        gsTensorDomain<real_t, 2> td(kvs);
        
        auto comp0 = td.component(0);
        auto comp1 = td.component(1);
        
        CHECK(comp0.get() != nullptr);
        CHECK(comp1.get() != nullptr);
    }

    TEST(gsDomain_Subdomain)
    {
        gsMatrix<> points(2, 3);
        points << 0.0, 0.5, 1.0,
                  0.0, 0.5, 1.0;
        
        gsPointDomain<> pd(points);
        
        // Single piece domain should return itself for subdomain 0
        auto sub = pd.subdomain(0);
        CHECK_EQUAL(sub->numElements(), pd.numElements());
        CHECK_EQUAL(pd.nPieces(), 1);
    }

    TEST(gsDomain_Print)
    {
        gsMatrix<> points(2, 5);
        points << 0.0, 0.25, 0.5, 0.75, 1.0,
                  0.0, 0.25, 0.5, 0.75, 1.0;
        
        gsPointDomain<> pd(points);
        
        std::ostringstream oss;
        oss << pd;
        std::string output = oss.str();
        
        // Check that output contains expected information
        CHECK(output.find("Domain") != std::string::npos);
        CHECK(output.find("2") != std::string::npos); // dimension
    }

    TEST(gsPointDomain_BoundaryIterator)
    {
        gsMatrix<> points(2, 3);
        points << 0.0, 0.5, 1.0,
                  0.0, 0.5, 1.0;
        
        gsPointDomain<> pd(points);
        
        // Point domain should have no boundary elements
        auto it = pd.beginBdr(boundary::none);
        CHECK(it == pd.endBdr(boundary::none));
    }

    TEST(gsTensorDomain_3D)
    {
        gsKnotVector<> kv1(0, 1, 0, 2);
        gsKnotVector<> kv2(0, 1, 0, 3);
        gsKnotVector<> kv3(0, 1, 0, 2);
        
        std::vector<gsDomain<>::Ptr> kvs;
        kvs.push_back(memory::make_shared_not_owned(&kv1));
        kvs.push_back(memory::make_shared_not_owned(&kv2));
        kvs.push_back(memory::make_shared_not_owned(&kv3));
        
        gsTensorDomain<real_t, 3> td(kvs);
        
        CHECK_EQUAL(td.dim(), 3);
        size_t expectedElements = kv1.numElements() * kv2.numElements() * kv3.numElements();
        CHECK_EQUAL(td.numElements(), expectedElements);
    }

    // Additional comprehensive tests for gsDomain module

    TEST(gsPointDomain_1D)
    {
        gsMatrix<> points(1, 10);
        for (int i = 0; i < 10; i++)
            points(0, i) = i * 0.1;
        
        gsPointDomain<> pd(points);
        CHECK_EQUAL(pd.dim(), 1);
        CHECK_EQUAL(pd.numElements(), 10);
    }

    TEST(gsPointDomain_3D)
    {
        gsMatrix<> points(3, 8);
        for (int i = 0; i < 8; i++)
        {
            points(0, i) = (i % 2) * 0.5;
            points(1, i) = ((i / 2) % 2) * 0.5;
            points(2, i) = (i / 4) * 0.5;
        }
        
        gsPointDomain<> pd(points);
        CHECK_EQUAL(pd.dim(), 3);
        CHECK_EQUAL(pd.numElements(), 8);
    }

    TEST(gsPointDomain_IteratorWithLargeSet)
    {
        gsMatrix<> points(2, 100);
        for (int i = 0; i < 100; i++)
        {
            points(0, i) = i * 0.01;
            points(1, i) = i * 0.01;
        }
        
        gsPointDomain<> pd(points);
        
        int count = 0;
        for (auto it = pd.beginAll(); it != pd.endAll(); ++it)
            count++;
        
        CHECK_EQUAL(count, 100);
    }

    TEST(gsTensorDomain_BoundaryIterator_West)
    {
        gsKnotVector<> kv1(0, 1, 1, 3);
        gsKnotVector<> kv2(0, 1, 1, 4);
        
        std::vector<gsDomain<>::Ptr> kvs;
        kvs.push_back(memory::make_shared_not_owned(&kv1));
        kvs.push_back(memory::make_shared_not_owned(&kv2));
        
        gsTensorDomain<real_t, 2> td(kvs);
        
        int count = 0;
        for (auto it = td.beginBdr(boundary::west); it != td.endBdr(boundary::west); ++it)
            count++;
        
        CHECK(count > 0);
    }

    TEST(gsTensorDomain_BoundaryIterator_AllSides)
    {
        gsKnotVector<> kv1(0, 1, 1, 2);
        gsKnotVector<> kv2(0, 1, 1, 2);
        
        std::vector<gsDomain<>::Ptr> kvs;
        kvs.push_back(memory::make_shared_not_owned(&kv1));
        kvs.push_back(memory::make_shared_not_owned(&kv2));
        
        gsTensorDomain<real_t, 2> td(kvs);
        
        // Test all boundary sides
        for (int side = 1; side <= 4; side++)
        {
            size_t nBdr = td.numElementsBdr(boxSide(side));
            CHECK(nBdr > 0);
        }
    }

    TEST(gsTensorDomain_MultipleComponents)
    {
        gsKnotVector<> kv1(0, 1, 2, 3);
        gsKnotVector<> kv2(0, 1, 1, 4);
        gsKnotVector<> kv3(0, 1, 3, 2);
        
        std::vector<gsDomain<>::Ptr> kvs;
        kvs.push_back(memory::make_shared_not_owned(&kv1));
        kvs.push_back(memory::make_shared_not_owned(&kv2));
        kvs.push_back(memory::make_shared_not_owned(&kv3));
        
        gsTensorDomain<real_t, 3> td(kvs);
        
        CHECK_EQUAL(td.dim(), 3);
        CHECK_EQUAL(td.degree(0), 2);
        CHECK_EQUAL(td.degree(1), 1);
        CHECK_EQUAL(td.degree(2), 3);
    }

    TEST(gsPointDomain_EmptyPoints)
    {
        gsMatrix<> points(2, 0);
        gsPointDomain<> pd(points);
        
        CHECK_EQUAL(pd.dim(), 2);
        CHECK_EQUAL(pd.numElements(), 0);
    }

    TEST(gsTensorDomain_UniformKnotVector)
    {
        gsKnotVector<> kv(0, 1, 2, 5, 1); // Uniform knots
        
        std::vector<gsDomain<>::Ptr> kvs;
        kvs.push_back(memory::make_shared_not_owned(&kv));
        
        gsTensorDomain<real_t, 1> td(kvs);
        
        CHECK_EQUAL(td.numElements(), 5);
        CHECK_EQUAL(td.degree(0), 2);
    }

    TEST(gsTensorDomain_NonUniformKnotVector)
    {
        gsKnotVector<> kv;
        kv.initClamped(0, 1, 3, 2); // degree 3, 2 interior knots
        
        std::vector<gsDomain<>::Ptr> kvs;
        kvs.push_back(memory::make_shared_not_owned(&kv));
        
        gsTensorDomain<real_t, 1> td(kvs);
        
        CHECK_EQUAL(td.degree(0), 3);
        CHECK(td.numElements() > 0);
    }

    TEST(gsTensorDomain_LargeDimensionBoundingBox)
    {
        gsKnotVector<> kv1(-5, 5, 1, 3);
        gsKnotVector<> kv2(-10, 10, 1, 4);
        gsKnotVector<> kv3(0, 100, 1, 2);
        
        std::vector<gsDomain<>::Ptr> kvs;
        kvs.push_back(memory::make_shared_not_owned(&kv1));
        kvs.push_back(memory::make_shared_not_owned(&kv2));
        kvs.push_back(memory::make_shared_not_owned(&kv3));
        
        gsTensorDomain<real_t, 3> td(kvs);
        
        gsMatrix<> bbox = td.boundingBox();
        CHECK_EQUAL(bbox.rows(), 3);
        CHECK_EQUAL(bbox.cols(), 2);
        CHECK_CLOSE(bbox(0, 0), -5.0, 1e-10);
        CHECK_CLOSE(bbox(0, 1), 5.0, 1e-10);
        CHECK_CLOSE(bbox(1, 0), -10.0, 1e-10);
        CHECK_CLOSE(bbox(1, 1), 10.0, 1e-10);
        CHECK_CLOSE(bbox(2, 0), 0.0, 1e-10);
        CHECK_CLOSE(bbox(2, 1), 100.0, 1e-10);
    }

    TEST(gsPointDomain_SinglePoint)
    {
        gsMatrix<> points(3, 1);
        points << 1.0, 2.0, 3.0;
        
        gsPointDomain<> pd(points);
        
        CHECK_EQUAL(pd.dim(), 3);
        CHECK_EQUAL(pd.numElements(), 1);
        
        const gsMatrix<>& pts = pd.points();
        CHECK_CLOSE(pts(0, 0), 1.0, 1e-10);
        CHECK_CLOSE(pts(1, 0), 2.0, 1e-10);
        CHECK_CLOSE(pts(2, 0), 3.0, 1e-10);
    }

    TEST(gsTensorDomain_HighDegree)
    {
        gsKnotVector<> kv(0, 1, 5, 3); // degree 5
        
        std::vector<gsDomain<>::Ptr> kvs;
        kvs.push_back(memory::make_shared_not_owned(&kv));
        
        gsTensorDomain<real_t, 1> td(kvs);
        
        CHECK_EQUAL(td.degree(0), 5);
        CHECK(td.numElements() > 0);
    }

    // Tests for domain iterator behavior
    TEST(gsTensorDomain_IteratorIncrement)
    {
        gsKnotVector<> kv(0, 1, 1, 4);
        
        std::vector<gsDomain<>::Ptr> kvs;
        kvs.push_back(memory::make_shared_not_owned(&kv));
        
        gsTensorDomain<real_t, 1> td(kvs);
        
        auto it = td.beginAll();
        auto it2 = it;
        ++it2;
        
        // Ensure iterator advanced
        CHECK(it != it2);
    }

    TEST(gsTensorDomain_IteratorEquality)
    {
        gsKnotVector<> kv(0, 1, 1, 3);
        
        std::vector<gsDomain<>::Ptr> kvs;
        kvs.push_back(memory::make_shared_not_owned(&kv));
        
        gsTensorDomain<real_t, 1> td(kvs);
        
        auto it1 = td.beginAll();
        auto it2 = td.beginAll();
        
        CHECK(it1 == it2);
    }

    TEST(gsPointDomain_SubdomainCheck)
    {
        gsMatrix<> points(2, 5);
        for (int i = 0; i < 5; i++)
        {
            points(0, i) = i * 0.2;
            points(1, i) = i * 0.3;
        }
        
        gsPointDomain<> pd(points);
        
        // Test that subdomain returns valid pointer
        auto sub = pd.subdomain(0);
        CHECK(sub.get() != nullptr);
        CHECK_EQUAL(sub->nPieces(), 1);
    }

    TEST(gsTensorDomain_AllElementsRange)
    {
        gsKnotVector<> kv1(0, 1, 1, 5);
        gsKnotVector<> kv2(0, 1, 1, 5);
        
        std::vector<gsDomain<>::Ptr> kvs;
        kvs.push_back(memory::make_shared_not_owned(&kv1));
        kvs.push_back(memory::make_shared_not_owned(&kv2));
        
        gsTensorDomain<real_t, 2> td(kvs);
        
        // Test allElements() range
        int count = 0;
        for (auto elem : td.allElements())
        {
            count++;
        }
        
        CHECK(count > 0);
    }

    /* 
     * Step-by-step instructions for additional complex domain tests:
     * 
     * 1. gsHDomain tests (hierarchical domains):
     *    - Create gsHTensorBasis with multiple levels of refinement
     *    - Test gsHDomain construction from the basis
     *    - Test gsHDomainIterator for iterating over hierarchical elements
     *    - Test gsHDomainBoundaryIterator for boundary iteration
     *    - Test gsHDomainLeafIter for iterating over leaf nodes only
     *    - Verify correct level assignment for each element
     *    - Test element support queries on different levels
     * 
     * 2. gsCompositeDomain tests:
     *    - Create multiple sub-domains (tensor domains)
     *    - Combine them into a gsCompositeDomain
     *    - Test iteration over all pieces
     *    - Test subdomain access by index
     *    - Verify total element count across all pieces
     * 
     * 3. gsDomainIterator advanced tests:
     *    - Test element center calculation
     *    - Test element volume/area calculation  
     *    - Test Jacobian computation at element centers
     *    - Test element support queries
     *    - Test iterator arithmetic (+=, -=, distance)
     * 
     * 4. gsBreaksIterator tests:
     *    - Create a knot vector with repeated knots
     *    - Test iteration over break points (unique knots)
     *    - Verify multiplicities are correctly identified
     * 
     * 5. gsKdNode and gsHTree tests:
     *    - Build a k-d tree for spatial indexing
     *    - Test insertion and search operations
     *    - Test neighbor finding
     *    - Test tree balance properties
     */
}
