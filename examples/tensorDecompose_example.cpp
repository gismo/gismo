#include <gismo.h>
#include <gsDomain/gsTensorDomain.h>
#include <gsDomain/gsCompositeDomain.h>

using namespace gismo;

int main(int argc, char *argv[])
{
    // Create a simple 2D tensor domain using a tensor basis
    gsKnotVector<real_t> kv_u(0, 1, 6, 1);  // 6 elements
    gsKnotVector<real_t> kv_v(0, 1, 8, 1);  // 8 elements
    
    gsTensorBSplineBasis<2, real_t> basis(kv_u, kv_v);
    auto domain = basis.domain();
    
    gsInfo << "Creating a 2D tensor domain (" << kv_u.numElements() << "x" 
           << kv_v.numElements() << " elements)...\n";
    gsInfo << "Total elements: " << domain->numElements() << "\n\n";
    
    // Test decomposition into 2 pieces
    gsInfo << "Decomposing into 2 pieces:\n";
    auto decomposed_2 = domain->decompose(2);
    if (auto composite = std::dynamic_pointer_cast<gsCompositeDomain<real_t>>(decomposed_2))
    {
        gsInfo << "Number of subdomains: " << composite->nPieces() << "\n";
        for (size_t i = 0; i < composite->nPieces(); ++i) {
            gsInfo << "  Subdomain " << (i+1) << ": " 
                   << composite->subdomain(i)->numElements() << " elements\n";
        }
    }
    
    // Test decomposition into 4 pieces
    gsInfo << "\nDecomposing into 4 pieces:\n";
    auto decomposed_4 = domain->decompose(4);
    if (auto composite = std::dynamic_pointer_cast<gsCompositeDomain<real_t>>(decomposed_4))
    {
        gsInfo << "Number of subdomains: " << composite->nPieces() << "\n";
        for (size_t i = 0; i < composite->nPieces(); ++i) {
            gsInfo << "  Subdomain " << (i+1) << ": " 
                   << composite->subdomain(i)->numElements() << " elements\n";
        }
    }
    
    // Test decomposition into 6 pieces (should adjust since we can't get exactly 6 with 2D factors)
    gsInfo << "\nDecomposing into 6 pieces:\n";
    auto decomposed_6 = domain->decompose(6);
    if (auto composite = std::dynamic_pointer_cast<gsCompositeDomain<real_t>>(decomposed_6))
    {
        gsInfo << "Number of subdomains: " << composite->nPieces() << "\n";
        for (size_t i = 0; i < composite->nPieces(); ++i) {
            gsInfo << "  Subdomain " << (i+1) << ": " 
                   << composite->subdomain(i)->numElements() << " elements\n";
        }
    }
    
    gsInfo << "\n✓ Tensor domain decompose example completed.\n";
    return 0;
}
