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

    Author(s): H.M. Verhelst

*/

#include <gismo.h>

using namespace gismo;

// Helper function to test a given decomposition
void test_decomposition(const gsCompositeDomain<real_t>& parent, 
                        const gsMultiPatch<real_t>& mp,
                        index_t npieces)
{
    gsInfo << "\n=======================================================\n";
    gsInfo << "Testing decomposition into " << npieces << " piece(s)...\n";
    gsInfo << "=======================================================\n";

    auto decomposed = parent.decompose(npieces);
    gsInfo << "Decomposition resulted in " << decomposed->nPieces() << " subdomains.\n\n";

    // 1. Iterate through each subdomain individually
    gsInfo << "--- Iterating through each subdomain separately ---\n";
    for (index_t i = 0; i < decomposed->nPieces(); ++i)
    {
        auto sub = decomposed->subdomain(i);
        gsInfo << "Subdomain " << i << " has " << sub->numElements() << " elements.\n";
    }

    // 2. Plot the result
    std::string filename = "decomposed";
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
    std::string input;
    std::string output = ".";

    index_t numRef = 2;
    index_t numElev= 1;

    index_t numDecompositions = 4;

    gsCmdLine cmd("Domain decomposition test");
    cmd.addInt("r", "refinements", "Number of uniform refinements to apply to the geometry", numRef);
    cmd.addInt("e", "elevations", "Number of degree elevations to apply to the geometry", numElev);
    cmd.addInt("d", "decompositions", "Number of decomposition tests to run (tensor, localOptimal, optimal)", numDecompositions);
    cmd.addString("i", "input", "Input geometry file name", input);
    cmd.addString("o", "output", "Output directory for results", output);
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    GISMO_ENSURE(!input.empty(), "No input file specified.");
    GISMO_ENSURE(gsFileManager::fileExists(input), "Input file "+ input + " does not exist.");

    gsFileData<> fd(input);
    gsMultiPatch<> mp;
    fd.getFirst<gsMultiPatch<>>(mp);
    gsInfo << "Loaded geometry with " << mp.nPatches() << " patches.\n";

    // Apply refinements and elevations
    mp.degreeElevate(numElev);
    for (index_t i = 0; i < numRef; ++i)
        mp.uniformRefine();

    gsCompositeDomain<real_t> domain(mp);
    gsInfo << "Parent composite domain created with " << domain.numElements() << " total elements.\n";

    // Run the test cases
    test_decomposition(domain, mp, numDecompositions);

    gsInfo << "\n✓ Composite domain test completed.\n";

    return 0;
}

