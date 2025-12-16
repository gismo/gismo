/** @file indexSubDomain_example.cpp

    @brief Example of using gsIndexSubDomain to work with a subset of elements.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst
*/

#include <gismo.h>
#include <gsDomain/gsIndexSubDomain.h>

using namespace gismo;

int main(int argc, char* argv[])
{
    // Create a simple 2D tensor domain
    gsInfo << "Creating a 2D tensor domain...\n";
    
    // Create knot vectors for [0,1] with 4 elements in each direction
    auto kvU = memory::make_shared(new gsKnotVector<real_t>(0, 1, 4, 1));
    auto kvV = memory::make_shared(new gsKnotVector<real_t>(0, 1, 4, 1));
    
    std::vector<typename gsKnotVector<real_t>::Ptr> knotVectors = {kvU, kvV};
    gsTensorDomain<real_t, 2> parentDomain(knotVectors);
    
    gsInfo << "Parent domain has " << parentDomain.numElements() << " elements.\n";
    
    // Create a subdomain with specific elements: {0, 2, 5, 7}
    std::vector<index_t> indices = {0, 2, 5, 7};
    gsIndexSubDomain<real_t> subdomain(parentDomain, indices);
    
    gsInfo << "\nIndex subdomain contains " << subdomain.numElements() 
           << " elements: ";
    for (auto idx : subdomain.elementIndices())
        gsInfo << idx << " ";
    gsInfo << "\n";
    
    // Test iteration
    gsInfo << "\nIterating over subdomain elements:\n";
    auto it = subdomain.beginAll();
    auto endIt = subdomain.endAll();
    for (; it != endIt; ++it)
    {
        gsInfo << "  Iteration ID " << it.id() << ": global element index = " 
               << it.localId() << "\n";
    }
    
    // Test containment
    gsInfo << "\nTesting containment:\n";
    for (index_t i = 0; i < 10; ++i)
    {
        gsInfo << "  Element " << i << ": " 
               << (subdomain.contains(i) ? "in subdomain" : "not in subdomain") << "\n";
    }
    
    gsInfo << "\n✓ gsIndexSubDomain example completed successfully.\n";
    
    return 0;
}
