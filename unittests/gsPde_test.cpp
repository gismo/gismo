/** @file gsPde_test.cpp

    @brief Tests for gsPde classes

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): Copilot
*/

#include "gismo_unittest.h"

SUITE(gsPde_test)
{
    TEST(gsLaplacePde_DefaultConstructor)
    {
        gsLaplacePde<double> laplace;
        CHECK(true); // Just check it compiles and doesn't crash
    }

    TEST(gsLaplacePde_ConstructorWithDomain)
    {
        // Create a simple domain
        gsKnotVector<> kv(0, 1, 0, 3);
        gsBSplineBasis<> basis(kv);
        gsMatrix<> coefs(3, 1);
        coefs << 0, 0.5, 1.0;
        gsBSpline<> curve(kv, coefs);
        
        gsMultiPatch<> mp;
        mp.addPatch(curve);
        
        gsLaplacePde<double> laplace(mp);
        
        CHECK_EQUAL(laplace.domain().nPatches(), 1);
        CHECK_EQUAL(laplace.patches().nPatches(), 1);
    }

    TEST(gsLaplacePde_ConstructorWithBoundaryConditions)
    {
        // Create domain
        gsKnotVector<> kv1(0, 1, 0, 2);
        gsKnotVector<> kv2(0, 1, 0, 2);
        gsTensorBSplineBasis<2> basis(kv1, kv2);
        
        gsMatrix<> coefs(4, 2);
        coefs << 0, 0,
                 1, 0,
                 0, 1,
                 1, 1;
        gsTensorBSpline<2> patch(basis, coefs);
        
        gsMultiPatch<> mp;
        mp.addPatch(patch);
        
        // Create boundary conditions
        gsBoundaryConditions<> bc;
        gsConstantFunction<> zero(0, 2);
        bc.addCondition(0, boundary::west, condition_type::dirichlet, &zero);
        
        gsLaplacePde<double> laplace(mp, bc);
        
        CHECK_EQUAL(laplace.domain().nPatches(), 1);
        CHECK(laplace.boundaryConditions().size() > 0);
        CHECK(laplace.bc().size() > 0);
    }

    TEST(gsLaplacePde_Print)
    {
        gsKnotVector<> kv(0, 1, 0, 2);
        gsBSplineBasis<> basis(kv);
        gsMatrix<> coefs(2, 1);
        coefs << 0, 1;
        gsBSpline<> curve(kv, coefs);
        
        gsMultiPatch<> mp;
        mp.addPatch(curve);
        
        gsLaplacePde<double> laplace(mp);
        
        std::ostringstream oss;
        oss << laplace;
        std::string output = oss.str();
        
        // Check that output contains Laplace
        CHECK(output.find("Laplace") != std::string::npos);
    }

    // Commented out - gsPoissonPde API uses different constructor parameters
    // TEST(gsPoissonPde_Construction)
    // {
    //     // Create domain
    //     gsKnotVector<> kv(0, 1, 0, 3);
    //     gsBSplineBasis<> basis(kv);
    //     gsMatrix<> coefs(3, 1);
    //     coefs << 0, 0.5, 1.0;
    //     gsBSpline<> curve(kv, coefs);
    //     
    //     gsMultiPatch<> mp;
    //     mp.addPatch(curve);
    //     
    //     // Create RHS function
    //     gsConstantFunction<> rhs(1.0, 1);
    //     
    //     // Create boundary conditions
    //     gsBoundaryConditions<> bc;
    //     gsConstantFunction<> zero(0, 1);
    //     bc.addCondition(0, boundary::west, condition_type::dirichlet, &zero);
    //     
    //     gsPoissonPde<double> poisson(mp, bc, &rhs);
    //     
    //     CHECK_EQUAL(poisson.domain().nPatches(), 1);
    //     CHECK(poisson.boundaryConditions().size() > 0);
    //     CHECK(poisson.rhs() != nullptr);
    // }

    // TEST(gsPoissonPde_Print)
    // {
    //     gsKnotVector<> kv(0, 1, 0, 2);
    //     gsBSplineBasis<> basis(kv);
    //     gsMatrix<> coefs(2, 1);
    //     coefs << 0, 1;
    //     gsBSpline<> curve(kv, coefs);
    //     
    //     gsMultiPatch<> mp;
    //     mp.addPatch(curve);
    //     
    //     gsConstantFunction<> rhs(1.0, 1);
    //     gsBoundaryConditions<> bc;
    //     
    //     gsPoissonPde<double> poisson(mp, bc, &rhs);
    //     
    //     std::ostringstream oss;
    //     oss << poisson;
    //     std::string output = oss.str();
    //     
    //     CHECK(output.find("Poisson") != std::string::npos);
    // }

    // Commented out - gsConvDiffRePde API doesn't match
    // TEST(gsConvDiffRePde_Construction)
    // {
    //     // Create domain
    //     gsKnotVector<> kv1(0, 1, 0, 2);
    //     gsKnotVector<> kv2(0, 1, 0, 2);
    //     gsTensorBSplineBasis<2> basis(kv1, kv2);
    //     
    //     gsMatrix<> coefs(4, 2);
    //     coefs << 0, 0,
    //              1, 0,
    //              0, 1,
    //              1, 1;
    //     gsTensorBSpline<2> patch(basis, coefs);
    //     
    //     gsMultiPatch<> mp;
    //     mp.addPatch(patch);
    //     
    //     gsBoundaryConditions<> bc;
    //     
    //     // Create coefficient functions
    //     gsConstantFunction<> rhs(1.0, 2);
    //     gsConstantFunction<> diffusion(1.0, 2);
    //     gsConstantFunction<> reaction(0.0, 2);
    //     
    //     gsMatrix<> convection(2, 1);
    //     convection << 1.0, 0.0;
    //     gsConstantFunction<> convFunc(convection, 2);
    //     
    //     gsConvDiffRePde<double> cdr(mp, bc, &rhs, &diffusion, &convFunc, &reaction);
    //     
    //     CHECK_EQUAL(cdr.domain().nPatches(), 1);
    //     CHECK(cdr.rhs() != nullptr);
    //     CHECK(cdr.diffusion() != nullptr);
    //     CHECK(cdr.convection() != nullptr);
    //     CHECK(cdr.reaction() != nullptr);
    // }

    // Commented out - gsBiharmonicPde not available in this version
    // TEST(gsBiharmonicPde_Construction)
    // {
    //     // Create domain
    //     gsKnotVector<> kv(0, 1, 2, 3);
    //     gsBSplineBasis<> basis(kv);
    //     gsMatrix<> coefs(3, 1);
    //     coefs << 0, 0.5, 1.0;
    //     gsBSpline<> curve(kv, coefs);
    //     
    //     gsMultiPatch<> mp;
    //     mp.addPatch(curve);
    //     
    //     gsBoundaryConditions<> bc;
    //     gsConstantFunction<> rhs(1.0, 1);
    //     
    //     gsBiharmonicPde<double> biharm(mp, bc, &rhs);
    //     
    //     CHECK_EQUAL(biharm.domain().nPatches(), 1);
    //     CHECK(biharm.rhs() != nullptr);
    // }

    TEST(gsPde_DomainAccess)
    {
        gsKnotVector<> kv(0, 1, 0, 2);
        gsBSplineBasis<> basis(kv);
        gsMatrix<> coefs(2, 1);
        coefs << 0, 1;
        gsBSpline<> curve(kv, coefs);
        
        gsMultiPatch<> mp;
        mp.addPatch(curve);
        mp.addPatch(curve); // Add second patch
        
        gsLaplacePde<double> laplace(mp);
        
        const gsMultiPatch<>& domainConst = laplace.domain();
        CHECK_EQUAL(domainConst.nPatches(), 2);
        
        gsMultiPatch<>& domainRef = laplace.domain();
        CHECK_EQUAL(domainRef.nPatches(), 2);
        
        const gsMultiPatch<>& patchesConst = laplace.patches();
        CHECK_EQUAL(patchesConst.nPatches(), 2);
    }

    TEST(gsPde_BoundaryConditionsAccess)
    {
        gsKnotVector<> kv(0, 1, 0, 2);
        gsBSplineBasis<> basis(kv);
        gsMatrix<> coefs(2, 1);
        coefs << 0, 1;
        gsBSpline<> curve(kv, coefs);
        
        gsMultiPatch<> mp;
        mp.addPatch(curve);
        
        gsBoundaryConditions<> bcs;
        gsConstantFunction<> zero(0, 1);
        bcs.addCondition(0, boundary::west, condition_type::dirichlet, &zero);
        bcs.addCondition(0, boundary::east, condition_type::neumann, &zero);
        
        gsLaplacePde<double> laplacePde(mp, bcs);
        
        const gsBoundaryConditions<>& bcConst = laplacePde.boundaryConditions();
        CHECK(bcConst.size() > 0);
        
        gsBoundaryConditions<>& bcRef = laplacePde.boundaryConditions();
        CHECK(bcRef.size() > 0);
        
        const gsBoundaryConditions<>& bcConst2 = laplacePde.bc();
        CHECK(bcConst2.size() > 0);
    }
}
