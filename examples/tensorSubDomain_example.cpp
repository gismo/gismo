/** @file tensorSubDomain_example.cpp

    @brief Example of using gsTensorSubDomain for efficient tensor-structured subdomains.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst
*/

#include <gismo.h>

using namespace gismo;

int main(int argc, char* argv[])
{
    // Create a 2D tensor domain with 5x6 elements
    gsInfo << "Creating a 2D tensor domain (5x6 elements)...\n";
    
    auto kvU = memory::make_shared(new gsKnotVector<real_t>(0, 1, 5, 1));
    auto kvV = memory::make_shared(new gsKnotVector<real_t>(0, 1, 6, 1));
    
    std::vector<typename gsKnotVector<real_t>::Ptr> knotVectors = {kvU, kvV};
    gsTensorDomain<real_t, 2> parentDomain(knotVectors);
    
    gsInfo << "Parent domain has " << parentDomain.numElements() << " elements (5x6).\n";
    
    // Create a tensor subdomain: u in [1,3], v in [2,5]
    // This creates a 2x3 = 6 element subdomain
    std::vector<gsTensorSubDomain<real_t, 2>::Range> ranges = {
        {1, 3},  // u direction: elements 1-2
        {2, 5}   // v direction: elements 2-4
    };
    
    gsTensorSubDomain<real_t, 2> tensorSubdomain(parentDomain, ranges, 0);
    
    gsInfo << "\nTensor subdomain [u ∈ [1,3), v ∈ [2,5)) contains " 
           << tensorSubdomain.numElements() << " elements (2x3).\n";
    
    // Get all element indices
    auto allIndices = tensorSubdomain.elementIndices();
    gsInfo << "\nElement indices:\n";
    for (size_t i = 0; i < allIndices.size(); ++i)
    {
        gsInfo << "  " << allIndices[i];
        if ((i + 1) % 3 == 0) gsInfo << "\n";
    }
    
    // Iterate over subdomain elements
    gsInfo << "\nIterating over tensor subdomain:\n";
    auto it = tensorSubdomain.beginAll();
    auto endIt = tensorSubdomain.endAll();
    index_t count = 0;
    for (; it != endIt; ++it)
    {
        gsInfo << "  Iteration " << count << ": local element index = " << it.localId() << "\n";
        ++count;
    }
    
    // Test containment
    gsInfo << "\nTesting containment:\n";
    for (index_t i = 0; i < 10; ++i)
    {
        gsInfo << "  Element " << i << ": "
               << (tensorSubdomain.contains(i) ? "in subdomain" : "not in subdomain") << "\n";
    }
    
    gsInfo << "\n✓ gsTensorSubDomain example completed successfully.\n";

    gsInfo << "\n----------------------------------------\n";
    gsInfo << "Partitioning Demonstration\n";
    gsInfo << "----------------------------------------\n";

    // Create first tensor subdomain: left half [u ∈ [0,2), v ∈ [0,6)]
    std::vector<gsTensorSubDomain<real_t, 2>::Range> ranges1 = {
        {0, 2},  // u: elements 0-1
        {0, 6}   // v: all elements
    };
    gsTensorSubDomain<real_t, 2> subdomain1(parentDomain, ranges1, 0);
    
    gsInfo << "Subdomain 1 [u ∈ [0,2), v ∈ [0,6)]: "
           << subdomain1.numElements() << " elements.\n";
    
    // Create second tensor subdomain: right half [u ∈ [2,5), v ∈ [0,6)]
    std::vector<gsTensorSubDomain<real_t, 2>::Range> ranges2 = {
        {2, 5},  // u: elements 2-4
        {0, 6}   // v: all elements
    };
    gsTensorSubDomain<real_t, 2> subdomain2(parentDomain, ranges2, 1);
    
    gsInfo << "Subdomain 2 [u ∈ [2,5), v ∈ [0,6)]: "
           << subdomain2.numElements() << " elements.\n";

    // Verify partition completeness
    gsInfo << "\nPartition analysis:\n";
    int inSub1 = 0, inSub2 = 0, inBoth = 0, inNeither = 0;
    
    for (size_t i = 0; i < parentDomain.numElements(); ++i)
    {
        bool in1 = subdomain1.contains(i);
        bool in2 = subdomain2.contains(i);
        
        if (in1 && in2) inBoth++;
        else if (in1) inSub1++;
        else if (in2) inSub2++;
        else inNeither++;
    }
    
    gsInfo << "Elements in subdomain 1 only: " << inSub1 << "\n";
    gsInfo << "Elements in subdomain 2 only: " << inSub2 << "\n";
    gsInfo << "Elements in both (overlap): " << inBoth << "\n";
    gsInfo << "Elements in neither (gap): " << inNeither << "\n";
    gsInfo << "Total coverage: " << (inSub1 + inSub2 + inBoth) << " / " 
           << parentDomain.numElements() << " elements\n";
    
    return 0;
}
