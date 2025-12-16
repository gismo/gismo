/** @file decompose_example.cpp

    @brief Example of domain decomposition using the decompose() method.

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
    gsInfo << "Domain Decomposition Example\n\n";
    gsInfo << "This example demonstrates precise domain decomposition\n";
    gsInfo << "of different domain types (Tensor, Multi-patch, Hierarchical).\n\n";

    // Example 1: Tensor domain decomposition with exact piece count
    gsInfo << "=== Example 1: Tensor Domain Decomposition ===\n";
    {
        // Create a 2D tensor domain with 5x4 = 20 elements
        gsKnotVector<real_t> kvU(0, 1, 5, 2);  // 5 elements
        gsKnotVector<real_t> kvV(0, 1, 4, 2);  // 4 elements
        
        gsTensorBSplineBasis<2, real_t> basis(kvU, kvV);
        auto domain = basis.domain();
        
        gsInfo << "Original tensor domain:\n";
        gsInfo << "  - Dimension: " << domain->dim() << "\n";
        gsInfo << "  - Number of elements: " << domain->numElements() << " (5x4)\n\n";
        
        // Test decomposition with exact piece counts
        std::vector<index_t> test_npieces = {1, 2, 4, 5, 8, 10};
        for (index_t npieces : test_npieces) {
            auto decomposed = domain->decompose(npieces);
            if (auto composite = std::dynamic_pointer_cast<gsCompositeDomain<real_t>>(decomposed)) {
                gsInfo << "Requested " << npieces << " pieces: Got " 
                       << composite->subdomains().size() << " subdomains\n";
                index_t totalElems = 0;
                for (size_t i = 0; i < composite->subdomains().size(); ++i) {
                    totalElems += composite->subdomain(i)->numElements();
                }
                gsInfo << "  Total elements: " << totalElems << "\n";
            }
        }
    }

    // Example 2: Multi-patch domain decomposition
    gsInfo << "\n=== Example 2: Multi-Patch Domain Decomposition ===\n";
    {
        // Create two patches using basis domains
        gsKnotVector<real_t> kvU1(0, 1, 3, 2);  // 3 elements
        gsKnotVector<real_t> kvV1(0, 1, 3, 2);  // 3 elements
        gsTensorBSplineBasis<2, real_t> basis1(kvU1, kvV1);
        auto patchA = basis1.domain();
        
        gsKnotVector<real_t> kvU2(0, 1, 4, 2);  // 4 elements
        gsKnotVector<real_t> kvV2(0, 1, 4, 2);  // 4 elements
        gsTensorBSplineBasis<2, real_t> basis2(kvU2, kvV2);
        auto patchB = basis2.domain();
        
        auto multiPatch = memory::make_shared(new gsCompositeDomain<real_t>());
        multiPatch->addDomain(patchA);
        multiPatch->addDomain(patchB);
        
        gsInfo << "Multi-patch domain (2 patches):\n";
        gsInfo << "  Patch A: 3x3 = 9 elements\n";
        gsInfo << "  Patch B: 4x4 = 16 elements\n";
        gsInfo << "  Total: " << multiPatch->numElements() << " elements\n\n";
        
        // Test cases
        gsInfo << "Case 1 - npieces == npatches (2):\n";
        {
            auto decomp = multiPatch->decompose(2);
            if (auto comp = std::dynamic_pointer_cast<gsCompositeDomain<real_t>>(decomp)) {
                gsInfo << "  Result: " << comp->subdomains().size() << " subdomains\n";
                for (size_t i = 0; i < comp->subdomains().size(); ++i) {
                    gsInfo << "    Subdomain " << i << ": " 
                           << comp->subdomain(i)->numElements() << " elements\n";
                }
            }
        }
        
        gsInfo << "\nCase 2 - npieces < npatches (1):\n";
        {
            auto decomp = multiPatch->decompose(1);
            if (auto comp = std::dynamic_pointer_cast<gsCompositeDomain<real_t>>(decomp)) {
                gsInfo << "  Result: " << comp->subdomains().size() << " subdomains\n";
                index_t totalElems = 0;
                for (size_t i = 0; i < comp->subdomains().size(); ++i) {
                    totalElems += comp->subdomain(i)->numElements();
                }
                gsInfo << "  Total elements: " << totalElems << "\n";
            }
        }
        
        gsInfo << "\nCase 3 - npieces > npatches (3):\n";
        {
            auto decomp = multiPatch->decompose(3);
            if (auto comp = std::dynamic_pointer_cast<gsCompositeDomain<real_t>>(decomp)) {
                gsInfo << "  Result: " << comp->subdomains().size() << " subdomains\n";
                for (size_t i = 0; i < comp->subdomains().size(); ++i) {
                    gsInfo << "    Subdomain " << i << ": " 
                           << comp->subdomain(i)->numElements() << " elements\n";
                }
            }
        }
    }

    // Example 3: Index subdomains
    gsInfo << "\n=== Example 3: Index Subdomains ===\n";
    {
        // Create a tensor domain
        gsKnotVector<real_t> kvU(0, 1, 3, 1);
        gsKnotVector<real_t> kvV(0, 1, 3, 1);
        
        gsTensorBSplineBasis<2, real_t> basis(kvU, kvV);
        auto domain = basis.domain();
        
        gsInfo << "Original domain has " << domain->numElements() << " elements\n";
        
        // Manually create index subdomains
        gsIndexSubDomain<real_t> sub1(*domain, {0, 1, 2, 3});
        gsIndexSubDomain<real_t> sub2(*domain, {4, 5, 6, 7});
        
        gsInfo << "Created 2 index subdomains:\n";
        gsInfo << "  - Subdomain 0 has " << sub1.numElements() << " elements\n";
        gsInfo << "  - Subdomain 1 has " << sub2.numElements() << " elements\n";
        
        // Iterate over first subdomain
        gsInfo << "\nIterating over first subdomain:\n";
        auto it = sub1.beginAll();
        auto endIt = sub1.endAll();
        index_t count = 0;
        for (; it != endIt; ++it, ++count)
        {
            if (count < 4)  // Print first 4
                gsInfo << "  Element " << count << ": iteration id=" << it.id() 
                       << ", local id=" << it.localId() << "\n";
        }
    }

    // Example 4: Tensor subdomain ranges
    gsInfo << "\n=== Example 4: Tensor Subdomain Ranges ===\n";
    {
        // Create a 2D tensor domain
        gsKnotVector<real_t> kvU(0, 1, 6, 2);
        gsKnotVector<real_t> kvV(0, 1, 4, 2);
        
        gsTensorBSplineBasis<2, real_t> basis(kvU, kvV);
        auto domain = basis.domain();
        
        auto tensorDomainPtr = std::dynamic_pointer_cast<gsTensorDomain<real_t, 2>>(domain);
        if (!tensorDomainPtr)
        {
            gsWarn << "Failed to cast to gsTensorDomain\n";
            return 1;
        }
        
        auto& tensorDomain = *tensorDomainPtr;
        
        gsInfo << "Tensor domain size: " << tensorDomain.numElements() << " elements (6x4)\n";
        
        // Create a subdomain covering half the domain in u direction
        std::vector<gsTensorSubDomain<real_t, 2>::Range> ranges = {
            {0, 3},  // First half in u direction (3 elements)
            {0, 4}   // All elements in v direction (4 elements)
        };
        
        gsTensorSubDomain<real_t, 2> subdom(tensorDomain, ranges, 0);
        
        gsInfo << "Created tensor subdomain:\n";
        gsInfo << "  - Range in u: [0, 3) = 3 elements\n";
        gsInfo << "  - Range in v: [0, 4) = 4 elements\n";
        gsInfo << "  - Number of elements: " << subdom.numElements() << "\n";
        gsInfo << "  - Expected: 3 * 4 = 12 elements\n";
    }

    gsInfo << "\n✓ Domain decomposition example completed.\n";

    return 0;
}
