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
    TEST(gsPointDomain)
    {
        // Construction test
        gsMatrix<> points(2, 5);
        points << 0.0, 0.25, 0.5, 0.75, 1.0,
                  0.0, 0.25, 0.5, 0.75, 1.0;
        gsPointDomain<> pd(points);
        CHECK_EQUAL(pd.dim(), 2);
        CHECK_EQUAL(pd.numElements(), 5);

        // Iterator test
        int count = 0;
        for (auto it = pd.beginAll(); it != pd.endAll(); ++it)
            count++;
        CHECK_EQUAL(count, 5);

        // Points access test
        const gsMatrix<>& pts = pd.points();
        CHECK_EQUAL(pts.rows(), 2);
        CHECK_EQUAL(pts.cols(), 5);
        CHECK_CLOSE(pts(0, 1), 0.25, 1e-10);

        // Additional tests for 1D, 3D, and empty points
        gsMatrix<> points1D(1, 10);
        for (int i = 0; i < 10; i++)
            points1D(0, i) = i * 0.1;
        gsPointDomain<> pd1D(points1D);
        CHECK_EQUAL(pd1D.dim(), 1);
        CHECK_EQUAL(pd1D.numElements(), 10);

        gsMatrix<> points3D(3, 8);
        for (int i = 0; i < 8; i++)
        {
            points3D(0, i) = (i % 2) * 0.5;
            points3D(1, i) = ((i / 2) % 2) * 0.5;
            points3D(2, i) = (i / 4) * 0.5;
        }
        gsPointDomain<> pd3D(points3D);
        CHECK_EQUAL(pd3D.dim(), 3);
        CHECK_EQUAL(pd3D.numElements(), 8);

        gsMatrix<> emptyPoints(2, 0);
        gsPointDomain<> pdEmpty(emptyPoints);
        CHECK_EQUAL(pdEmpty.dim(), 2);
        CHECK_EQUAL(pdEmpty.numElements(), 0);
    }

    TEST(gsTensorDomain)
    {
        // 1D Tensor Domain
        gsKnotVector<> kv(0, 1, 2, 4);
        std::vector<gsDomain<>::Ptr> kvs1D;
        kvs1D.push_back(memory::make_shared_not_owned(&kv));
        gsTensorDomain<real_t, 1> td1D(kvs1D);
        CHECK_EQUAL(td1D.dim(), 1);
        CHECK_EQUAL(td1D.degree(0), kv.degree());
        CHECK(td1D.numElements() > 0);

        // 2D Tensor Domain
        gsKnotVector<> kv1(0, 1, 1, 3);
        gsKnotVector<> kv2(0, 1, 1, 3);
        std::vector<gsDomain<>::Ptr> kvs2D;
        kvs2D.push_back(memory::make_shared_not_owned(&kv1));
        kvs2D.push_back(memory::make_shared_not_owned(&kv2));
        gsTensorDomain<real_t, 2> td2D(kvs2D);
        CHECK_EQUAL(td2D.dim(), 2);
        CHECK_EQUAL(td2D.degree(0), kv1.degree());
        CHECK_EQUAL(td2D.degree(1), kv2.degree());

        // Additional tests for bounding box, iterators, and components
        gsMatrix<> bbox = td2D.boundingBox();
        CHECK_EQUAL(bbox.rows(), 2);
        CHECK_EQUAL(bbox.cols(), 2);
        CHECK_CLOSE(bbox(0, 0), 0.0, 1e-10);
        CHECK_CLOSE(bbox(0, 1), 1.0, 1e-10);
        CHECK_CLOSE(bbox(1, 0), 0.0, 1e-10);
        CHECK_CLOSE(bbox(1, 1), 1.0, 1e-10);

        int count = 0;
        for (auto it = td2D.beginAll(); it != td2D.endAll(); ++it)
            count++;
        CHECK_EQUAL(count, (int)td2D.numElements());

        auto comp0 = td2D.component(0);
        auto comp1 = td2D.component(1);
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
        for (auto & elem : td.allElements())
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
