/** @file hDomain_example.cpp

    @brief Example of iterating through the elements of a hierarchical domain (gsHDomain).

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst
*/

#include <gismo.h>
#include <gsDomain/gsHDomainIterator.h>

using namespace gismo;

int main(int argc, char* argv[])
{
    gsInfo << "Hierarchical Domain Example\n";
    gsInfo << "Demonstrates iterating over the elements of a gsHDomain from a THB-spline basis.\n\n";
    
    // Create a simple 2D THB-spline basis manually by refining a coarse B-spline basis
    gsInfo << "Creating a 2D THB-spline basis with 3 levels of refinement...\n";
    
    int p = 2;  // degree
    int n = 4;  // initial basis functions per direction
    
    gsKnotVector<real_t> kvU(0, 1, n, p);
    gsKnotVector<real_t> kvV(0, 1, n, p);
    
    gsTensorBSplineBasis<2, real_t> coarseBasis(kvU, kvV);
    gsTHBSplineBasis<2, real_t> hBasis(coarseBasis);
    hBasis.refineElements({0, 1, 2, 3}); // Refine some elements on level 0
    hBasis.refineElements({4, 5});      // Refine some elements on level 1
    hBasis.refineElements({10});        // Refine an element on level 2

    // Get the hierarchical domain from the basis
    auto domain = hBasis.domain();
    
    // Check if it's a gsHDomain and iterate through its elements
    if (auto hdom = std::dynamic_pointer_cast<gsHDomain<2, real_t>>(domain)) 
    {
        gsInfo << "\nSuccessfully obtained a gsHDomain with " << hdom->numElements() << " active elements.\n";
        
        // To access hierarchical-specific information like an element's level,
        // we must use a traditional iterator loop and dynamic_cast the iterator
        // to the derived gsHDomainIterator type.
        gsInfo << "Iterating over all active elements in the hierarchical domain:\n";
        index_t count = 0;
        for (auto it = hdom->beginAll(); it != hdom->endAll(); ++it)
        {
            if (auto h_it = dynamic_cast<const gsHDomainIterator<real_t, 2>*>(it.get()))
            {
                if (count < 10)
                {
                    gsInfo << "  Element " << count << ": id=" << h_it->id() 
                           << " (level=" << h_it->getLevel() << ")"
                           << " (patch=" << h_it->patch() << ")" << "\n";
                }
                count++;
            }
        }

        if (count > 10)
            gsInfo << "  ... (total " << count << " elements)\n";
    }
    else
    {
        gsInfo << "Could not cast basis domain to gsHDomain.\n";
    }
    
    gsInfo << "\n✓ Hierarchical Domain example completed successfully.\n";
    
    return 0;
}
