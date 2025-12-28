#include <gismo.h>
#include <gsDomain/gsHDomain.h>
#include <gsHSplines/gsHTensorBasis.h>    // Base class declaration
#include <gsHSplines/gsTHBSplineBasis.h> // Concrete implementation
#include <gsNurbs/gsKnotVector.h>        // Knot vectors
#include <gsNurbs/gsTensorBSplineBasis.h> // Base B-spline for THB
#include <gsCore/gsMemory.h>
#include <gsDomain/gsHDomainIterator.h> // Explicitly include for direct access to gsHDomainIterator

using namespace gismo;

int main(int argc, char* argv[])
{
    const short_t dim = 2;
    index_t degree = 1;
    index_t npieces = 4;
    gsCmdLine cmd("Example demonstrating gsHDomain creation and decomposition.");
    cmd.addInt("d", "degree", "Degree of the base B-spline basis", degree);
    cmd.addInt("n", "numpieces", "Number of pieces to decompose into", npieces);
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    gsInfo << "Starting gsHDomain example..." << std::endl;

    gsKnotVector<real_t> kv_u(0.0,1.0,0,degree+1, degree);
    gsKnotVector<real_t> kv_v = kv_u;
    gsTensorBSplineBasis<dim, real_t> base_tbasis(kv_u, kv_v);

    gsTHBSplineBasis<dim, real_t> thb_basis(base_tbasis);
    std::vector<index_t> refineBox1(2*dim+1);
    std::vector<index_t> refineBox2(2*dim+1);
    index_t maxLvl = 5;
    for (index_t level = 0; level < maxLvl; ++level)
    {
        refineBox1[0] = level + 1;
        refineBox2[0] = level + 1;
        for (short_t i = 0; i < dim; ++i)
        {
            refineBox1[1 + i] = 0; // lower corner
            refineBox1[1 + dim + i] = 2; // upper corner
            refineBox2[1 + i] = 2*math::pow(2, level)-1; // lower corner
            refineBox2[1 + dim + i] = 2*math::pow(2, level); // upper corner
        }
        thb_basis.refineElements(refineBox1);
        if (level == maxLvl - 1)
            break;
        thb_basis.refineElements(refineBox2);
    }

    const gsHTree<dim, index_t>& h_tree = thb_basis.tree();
    gsHDomain<dim, real_t, index_t> h_domain(h_tree, thb_basis);
    gsInfo << "Hierarchical domain has " << h_domain.numElements() << " elements at the finest level.\n";

    gsInfo << "\nAttempting to decompose hierarchical domain into " << npieces << " pieces...\n";
    typename gsDomain<real_t>::Ptr decomposed_h_domain = h_domain.decompose(static_cast<size_t>(npieces));

    if (decomposed_h_domain) {
        gsInfo << "Decomposition resulted in " << decomposed_h_domain->nPieces() << " subdomains.\n";
        gsInfo << "Total elements covered by decomposed domain: " << decomposed_h_domain->numElements() << "\n";

        for (size_t i = 0; i < decomposed_h_domain->nPieces(); ++i)
        {
            gsInfo << "  Subdomain " << i << " has " << decomposed_h_domain->subdomain(i)->numElements() << " elements.\n";
        }
    } else {
        gsInfo << "Decomposition failed or returned null.\n";
    }

    gsWriteParaview(*decomposed_h_domain, "hDomain_decomposed");

    gsInfo << "gsHDomain example finished successfully." << std::endl;

    return 0;
}