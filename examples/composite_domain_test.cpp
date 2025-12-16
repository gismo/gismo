/** @file composite_domain_test.cpp

    @brief Example to test gsCompositeDomain and its iterators.

    This example creates a 4-patch gsMultiPatch geometry and then tests the
    gsCompositeDomain::decompose method for 1, 2, 4, and 8 pieces.
    For each case, it iterates through the resulting subdomains and prints the
    patch() ID from the iterators to verify their correctness.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

*/

#include <gismo.h>

using namespace gismo;

// Helper function to test a given decomposition
void test_decomposition(const gsCompositeDomain<real_t>& parent, 
                        const gsMultiPatch<real_t>& mp,
                        index_t npieces,
                        gismo::decompositionStrategy strategy)
{
    std::string strategy_str = (strategy == gismo::decompositionStrategy::tensor ? "tensor" :
                                (strategy == gismo::decompositionStrategy::localOptimalBalancing ? "localOptimal" : "optimal"));
    
    gsInfo << "\n=======================================================\n";
    gsInfo << "Testing " << strategy_str << " decomposition into " << npieces << " piece(s)...\n";
    gsInfo << "=======================================================\n";

    auto decomposed = parent.decompose(npieces, strategy);
    gsInfo << "Decomposition resulted in " << decomposed->nPieces() << " subdomains.\n\n";

    // 1. Iterate through each subdomain individually
    gsInfo << "--- Iterating through each subdomain separately ---\n";
    for (index_t i = 0; i < decomposed->nPieces(); ++i)
    {
        auto sub = decomposed->subdomain(i);
        gsInfo << "Subdomain " << i << " has " << sub->numElements() << " elements.\n";
    }

    // 2. Plot the result
    std::string filename = "decomposed_" + strategy_str + "_" + std::to_string(npieces);
    gsWriteParaview(*decomposed, mp, filename, 2);
    gsInfo << "Plot of decomposed domain saved to " << filename << ".pvd\n";

    // 3. Verify that all elements are covered
    gsInfo << "\n--- Verifying element coverage ---\n";
    if (parent.numElements() != decomposed->numElements())
    {
        gsWarn << "Error: Number of elements in parent and decomposed domains do not match!\n";
        gsWarn << "Parent has " << parent.numElements() << " elements.\n";
        gsWarn << "Decomposed has " << decomposed->numElements() << " elements.\n";
    }
    else
    {
        gsInfo << "✓ Total number of elements is correct.\n";
    }

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

    if (all_parent_elements.size() != decomposed_elements.size())
    {
        gsWarn << "Error: Number of unique elements in parent and decomposed domains do not match!\n";
        gsWarn << "Parent has " << all_parent_elements.size() << " unique elements.\n";
        gsWarn << "Decomposed has " << decomposed_elements.size() << " unique elements.\n";
    }
    
    size_t not_covered = 0;
    for(const auto& elem : all_parent_elements)
    {
        if (decomposed_elements.find(elem) == decomposed_elements.end())
        {
            gsWarn << "Error: Element with corners [" << elem.first << "] and [" << elem.second << "] is not covered by the decomposition!\n";
            not_covered++;
        }
    }

    if (not_covered == 0)
    {
        gsInfo << "✓ All elements are covered by the decomposition.\n";
    }
    else
    {
        gsWarn << not_covered << " elements not covered.\n";
    }
}

int main(int argc, char* argv[])
{
    gsInfo << "Creating a 3x3 single-patch geometry...";
    gsMultiPatch<real_t> singlepatch = gsNurbsCreator<real_t>::BSplineSquareGrid(1,1);
    singlepatch.uniformRefine(3); 
    gsCompositeDomain<real_t> parent_domain_single(singlepatch);
    gsInfo << "Parent composite domain created with " << parent_domain_single.numElements() << " total elements.\n";

    test_decomposition(parent_domain_single, singlepatch, 2, decompositionStrategy::tensor);
    test_decomposition(parent_domain_single, singlepatch, 2, decompositionStrategy::localOptimalBalancing);

    gsInfo << "\nCreating a 4-patch geometry (2x2 grid of b-spline squares)...";

    // Use gsNurbsCreator to make a simple 2x2 multi-patch geometry
    gsMultiPatch<real_t> multipatch = gsNurbsCreator<real_t>::BSplineSquareGrid(2,2);
    multipatch.uniformRefine(2); // Refine for more elements

    gsInfo << "Geometry created with " << multipatch.nPatches() << " patches.\n";
    for (index_t i = 0; i < multipatch.nPatches(); ++i)
    {
        gsInfo << "  Patch " << i << " has " << multipatch.basis(i).domain()->numElements() << " elements.\n";
    }

    // Create a composite domain from the multipatch
    gsCompositeDomain<real_t> parent_domain_multi(multipatch);
    gsInfo << "Parent composite domain created with " << parent_domain_multi.numElements() << " total elements.\n";

    gsInfo << "\n--- Parent domain elements (physical corners) ---\n";
    for(auto it = parent_domain_multi.beginAll(); it != parent_domain_multi.endAll(); ++it)
    {
        const auto& patch = multipatch.patch(it.patch());
        gsMatrix<real_t> corners(2,2);
        corners.col(0) = it.lowerCorner();
        corners.col(1) = it.upperCorner();
        gsMatrix<real_t> phys_corners = patch.eval(corners);

        gsInfo << "Patch: " << it.patch() << ", LocalID: " << it.localId() 
               << ", Lower: " << phys_corners.col(0).transpose() 
               << ", Upper: " << phys_corners.col(1).transpose() << "\n";
    }

    // Run the test cases
    test_decomposition(parent_domain_multi, multipatch, 8, decompositionStrategy::tensor);
    test_decomposition(parent_domain_multi, multipatch, 8, decompositionStrategy::localOptimalBalancing);
    test_decomposition(parent_domain_multi, multipatch, 8, decompositionStrategy::optimalBalancing);

    gsInfo << "\n✓ Composite domain test completed.\n";

    return 0;
}

