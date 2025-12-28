#include <gismo.h>
#include <gsDomain/gsIndexSubDomain.h>
#include <gsCore/gsMemory.h>
#include <gsNurbs/gsKnotVector.h>      // Explicitly include for KnotVector
#include <gsDomain/gsTensorDomain.h>   // Explicitly include for TensorDomain

using namespace gismo;

int main(int argc, char* argv[])
{
    gsInfo << "Starting gsIndexSubDomain example..." << std::endl;

    // 1. Create a simple parent domain (e.g., a tensor domain from a single patch)
    // Degree for the basis
    short_t p = 1; 

    gsInfo << "Before kv_u construction." << std::endl;
    // Create knot vectors with explicit values for more elements
    std::vector<real_t> knots_u = {0.0, 0.0, 0.5, 1.0, 1.0}; // Degree 1, 2 elements
    gsKnotVector<real_t> kv_u(knots_u, p);
    gsInfo << "After kv_u construction." << std::endl;

    gsInfo << "Before kv_v construction." << std::endl;
    std::vector<real_t> knots_v = {0.0, 0.0, 0.5, 1.0, 1.0}; // Degree 1, 2 elements
    gsKnotVector<real_t> kv_v(knots_v, p);
    gsInfo << "After kv_v construction." << std::endl;

    gsInfo << "Before parentDomain construction." << std::endl;
    gsTensorDomain<real_t, 2> parentDomain(kv_u, kv_v); 
    gsInfo << "After parentDomain construction." << std::endl;

    gsInfo << "Parent domain has " << parentDomain.numElements() << " elements." << std::endl;

    // 2. Define a subset of elements for the index subdomain
    // Let's take the first half of the elements for the subdomain.
    std::vector<index_t> sub_indices;
    size_t num_elements = parentDomain.numElements();
    for (size_t i = 0; i < num_elements / 2; ++i) {
        sub_indices.push_back(i);
    }

    gsInfo << "Before indexSubDomain construction." << std::endl;
    gsIndexSubDomain<real_t> indexSubDomain(parentDomain, std::move(sub_indices));
    gsInfo << "After indexSubDomain construction." << std::endl;

    gsInfo << "Index subdomain has " << indexSubDomain.numElements() << " elements." << std::endl;

    gsInfo << "Iterating over index subdomain elements:" << std::endl;
    size_t count = 0;
    for (auto it = indexSubDomain.beginAll(); it != indexSubDomain.endAll(); ++it) {
        gsInfo << "  Subdomain Element " << count << ": Local ID = " << it.localId() 
               << ", Parent Global ID = " << it.id() << std::endl;
        count++;
    }

    gsInfo << "Bounding box of the first element in the index subdomain (if available):" << std::endl;
    if (indexSubDomain.numElements() > 0) {
        auto it = indexSubDomain.beginAll();
        gsInfo << "  Lower Corner: " << it.lowerCorner().transpose() << std::endl;
        gsInfo << "  Upper Corner: " << it.upperCorner().transpose() << std::endl;
    }

    gsInfo << "gsIndexSubDomain example finished successfully." << std::endl;

    return 0;
}
