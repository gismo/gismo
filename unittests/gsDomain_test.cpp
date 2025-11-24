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
}
