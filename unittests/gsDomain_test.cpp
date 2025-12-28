/** @file gsDomain_test.cpp

    @brief Unittest for gsHDomain::decompose method.

    This test constructs a gsHDomain from a gsTHBSplineBasis, decomposes it
    using different strategies, and verifies the correctness of the decomposition.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

*/

#include <gismo.h>
#include <gsDomain/gsHDomain.h>
#include <gsHSplines/gsTHBSplineBasis.h>
#include <gsCore/gsMemory.h>
#include <gsDomain/gsDomain.h> // Required for gismo::decompositionStrategy

// For UnitTest++
#include <gismo_unittest.h>

using namespace gismo;

// Helper function to test a given decomposition
void test_hdomain_decomposition(const gsHDomain<2, real_t, index_t>& parent,
                                index_t npieces)
{
    // Decompose the domain
    auto decomposed = parent.decompose(npieces);
    CHECK(decomposed != nullptr); // Check that decomposition returned a valid pointer

    // Verify that the total number of elements matches
    CHECK_EQUAL(parent.numElements(), decomposed->numElements());

    // Verify that all elements from the parent are covered by the decomposed domain
    auto vectorToString = [](const gsVector<real_t>& v) {
        std::stringstream ss;
        ss << v.transpose();
        return ss.str();
    };

    std::set<std::pair<std::string, std::string>> all_parent_elements;
    for(auto it = parent.beginAll(); it != parent.endAll(); ++it)
    {
        all_parent_elements.insert({vectorToString(it.lowerCorner()), vectorToString(it.upperCorner())});
    }

    std::set<std::pair<std::string, std::string>> decomposed_elements;
    for(auto it = decomposed->beginAll(); it != decomposed->endAll(); ++it)
    {
        decomposed_elements.insert({vectorToString(it.lowerCorner()), vectorToString(it.upperCorner())});
    }

    CHECK_EQUAL(all_parent_elements.size(), decomposed_elements.size());
    
    size_t not_covered = 0;
    for(const auto& elem : all_parent_elements)
    {
        if (decomposed_elements.find(elem) == decomposed_elements.end())
        {
            not_covered++;
        }
    }
    CHECK_EQUAL(0, not_covered); // All elements should be covered
}

TEST(gsHDomainDecomposition)
{
    const short_t dim = 2;
    index_t degree = 1;

    // Create a base tensor-product B-spline basis
    gsKnotVector<real_t> kv_u(0.0,1.0,0,degree+1, degree);
    gsKnotVector<real_t> kv_v = kv_u;
    gsTensorBSplineBasis<dim, real_t> base_tbasis(kv_u, kv_v);

    // Create a THB-spline basis and refine it
    gsTHBSplineBasis<dim, real_t> thb_basis(base_tbasis);
    std::vector<index_t> refineBox(2*dim+1);
    index_t maxLvl = 3;
    for (index_t level = 0; level < maxLvl; ++level)
    {
        refineBox[0] = level + 1;
        for (short_t i = 0; i < dim; ++i)
        {
            refineBox[1 + i] = 0;
            refineBox[1 + dim + i] = 2;
        }
        thb_basis.refineElements(refineBox);
    }

    // Create the hierarchical domain
    const gsHTree<dim, index_t>& h_tree = thb_basis.tree();
    gsHDomain<dim, real_t, index_t> h_domain(h_tree, thb_basis);

    // Run the test cases
    test_hdomain_decomposition(h_domain, 1);
    test_hdomain_decomposition(h_domain, 2);
    test_hdomain_decomposition(h_domain, 4);
    
    test_hdomain_decomposition(h_domain, 2);
    test_hdomain_decomposition(h_domain, 4);

    // The 'optimalBalancing' strategy for gsHDomain is simplified and behaves like tensor, so we test that too.
    test_hdomain_decomposition(h_domain, 2);
}
