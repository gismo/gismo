/** @file elasticityAutodiff_rev_simple_test.cpp

    @brief Test reverse-mode AD through a simple elasticity problem
    
    Incrementally test var type support with minimal code
*/

#include <gismo.h>
#include <gsAutoDiff/gsAutoDiff2.h>
#include <gsAutoDiff/gsAutoDiffVar.h>

#include <gsElasticity/src/gsElasticityAssembler.h>
#include <gsElasticity/src/gsIterative.h>
#include <gsElasticity/src/gsLinearMaterial.h>

using namespace gismo;
using namespace autodiff;

int main()
{
    gsInfo << "=== Testing var type with elasticity ===\n";

    std::string filename = "cooks.xml";
    
    // Test 1: Can we create var-typed material and assembler objects?
    gsInfo << "\nTest 1: Creating var-typed material...\n";
    try {
        var E(240.565e6);
        var nu(0.4);
        gsLinearMaterial<var> mat(E, nu, 2);
        gsInfo << "✓ Material created successfully\n";
    } catch (const std::exception& e) {
        gsInfo << "✗ Failed: " << e.what() << "\n";
        return 1;
    }

    // Test 2: Can we create var-typed geometry and basis?
    gsInfo << "\nTest 2: Loading var-typed geometry and creating basis...\n";
    try {
        var E(240.565e6);
        var nu(0.4);
        
        gsMultiPatch<var> geometry;
        gsReadFile<var>(filename, geometry);
        gsInfo << "✓ Geometry loaded\n";
        
        gsMultiBasis<var> basis(geometry);
        basis.uniformRefine();
        gsInfo << "✓ Basis created and refined\n";
    } catch (const std::exception& e) {
        gsInfo << "✗ Failed: " << e.what() << "\n";
        return 1;
    }

    // Test 3: Can we assemble with var types?
    gsInfo << "\nTest 3: Assembling system with var types...\n";
    try {
        var E(240.565e6);
        var nu(0.4);
        
        gsMultiPatch<var> geometry;
        gsReadFile<var>(filename, geometry);
        gsMultiBasis<var> basis(geometry);
        basis.uniformRefine();
        
        // BCs
        gsConstantFunction<var> f(var(0.), var(625e4), 2);
        gsBoundaryConditions<var> bcInfo;
        for (index_t d = 0; d < 2; ++d)
            bcInfo.addCondition(0, boundary::west, condition_type::dirichlet, nullptr, d);
        bcInfo.addCondition(0, boundary::east, condition_type::neumann, &f);
        
        gsConstantFunction<var> g(var(0.), var(0.), 2);
        
        gsLinearMaterial<var> mat(E, nu, 2);
        gsElasticityAssembler<var> assembler(geometry, basis, bcInfo, g, &mat);
        
        gsInfo << "✓ Assembler created with " << assembler.numDofs() << " dofs\n";
        
        // Try to assemble
        gsInfo << "Attempting to assemble...\n";
        assembler.assemble();
        gsInfo << "✓ Assembly succeeded\n";
        
    } catch (const std::exception& e) {
        gsInfo << "✗ Failed: " << e.what() << "\n";
        return 1;
    }

    gsInfo << "\n=== All basic tests passed ===\n";
    return 0;
}
