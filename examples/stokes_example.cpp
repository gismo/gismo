/** @file stokes_example.cpp

    @brief Assembling the Stokes equations with the expression assembler

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): S. Takacs
*/

#include <ctime>
#include <gismo.h>

using namespace gismo;

int main(int argc, char *argv[])
{
    /************** Define command line options *************/

    std::string geometry("domain2d/yeti_mp2.xml");
    index_t splitPatches = 1;
    real_t stretchGeometry = 1;
    index_t refinements = 1;
    index_t degree = 2;
    std::string boundaryConditions("d");

    gsCmdLine cmd("Assembles the Stokes system.");
    cmd.addString("g", "Geometry",              "Geometry file", geometry);
    cmd.addInt   ("",  "SplitPatches",          "Split every patch that many times in 2^d patches", splitPatches);
    cmd.addReal  ("",  "StretchGeometry",       "Stretch geometry in x-direction by the given factor", stretchGeometry);
    cmd.addInt   ("r", "Refinements",           "Number of uniform h-refinement steps to perform before solving", refinements);
    cmd.addInt   ("p", "Degree",                "Degree of the B-spline discretization space", degree);
    cmd.addString("b", "BoundaryConditions",    "Boundary conditions", boundaryConditions);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    if ( ! gsFileManager::fileExists(geometry) )
    {
        gsInfo << "Geometry file could not be found.\n";
        gsInfo << "I was searching in the current directory and in: " << gsFileManager::getSearchPaths() << "\n";
        return EXIT_FAILURE;
    }

    gsInfo << "Run ieti_example with options:\n" << cmd << std::endl;

    /******************* Define geometry ********************/

    gsInfo << "Define geometry... " << std::flush;

    gsMultiPatch<>::uPtr mpPtr = gsReadFile<>(geometry);
    if (!mpPtr)
    {
        gsInfo << "No geometry found in file " << geometry << ".\n";
        return EXIT_FAILURE;
    }
    gsMultiPatch<>& mp = *mpPtr;

    for (index_t i=0; i<splitPatches; ++i)
    {
        gsInfo << "split patches uniformly... " << std::flush;
        mp = mp.uniformSplit();
    }

    if (stretchGeometry!=1)
    {
       gsInfo << "and stretch it... " << std::flush;
       for (size_t i=0; i!=mp.nPatches(); ++i)
           const_cast<gsGeometry<>&>(mp[i]).scale(stretchGeometry,0);
       // Const cast is allowed since the object itself is not const. Stretching the
       // overall domain keeps its topology.
    }

    gsInfo << "done.\n";

    /************** Define boundary conditions **************/

    gsInfo << "Define right-hand-side and boundary conditions... " << std::flush;

    // Right-hand-side
    gsFunctionExpr<> f1( "2*sin(x)*cos(y)", mp.geoDim() );
    gsFunctionExpr<> f2( "2*sin(x)*cos(y)", mp.geoDim() );

    // Dirichlet function
    gsFunctionExpr<> gD1( "sin(x)*cos(y)", mp.geoDim() );
    gsFunctionExpr<> gD2( "sin(x)*cos(y)", mp.geoDim() );

    // Neumann
    gsConstantFunction<> gN( 1.0, mp.geoDim() );

    gsBoundaryConditions<> bc0, bc1, bc2;

    for (gsMultiPatch<>::const_biterator it = mp.bBegin(); it < mp.bEnd(); ++it)
    {
        bc1.addCondition( *it, condition_type::dirichlet, &gD1 );
        bc2.addCondition( *it, condition_type::dirichlet, &gD2 );
    }

    /************ Setup bases and adjust degree *************/

    gsMultiBasis<> mbU(mp);
    gsMultiBasis<> mbV(mp);
    gsMultiBasis<> mbP(mp);

    gsInfo << "Setup bases and adjust degree... " << std::flush;

    for ( size_t i = 0; i < mbU.nBases(); ++ i )
    {
        mbU[i].setDegreePreservingMultiplicity(degree+1);
        mbV[i].setDegreePreservingMultiplicity(degree+1);
        mbP[i].setDegreePreservingMultiplicity(degree);
        mbU[i].reduceContinuity(1);
        mbV[i].reduceContinuity(1);
    }

    for ( index_t i = 0; i < refinements; ++i )
    {
        mbU.uniformRefine();
        mbV.uniformRefine();
        mbP.uniformRefine();
    }

    gsInfo << "done.\n";

    /********* Setup assembler and assemble matrix **********/

    gsInfo << "Setup assembler and assemble matrix... " << std::flush;

    // The usual stuff for the expression assembler
    typedef gsExprAssembler<>::geometryMap geometryMap;
    typedef gsExprAssembler<>::space       space;
    typedef gsExprAssembler<>::solution    solution;

    // We set up the assembler
    gsExprAssembler<> assembler(3,3);

    // Elements used for numerical integration
    assembler.setIntegrationElements(mbU);  // All grids are the same, so we can choose any
    gsExprEvaluator<> ev(assembler);

    // Set the geometry map
    geometryMap G = assembler.getMap(mp);

    // Set the discretization spaces
    // We define different spaces for the two velocity components since we want
    // to be able to have different spaces (like for Nedelec elements) in future.
    space u = assembler.getSpace(mbU,1,0);
    space v = assembler.getSpace(mbV,1,1);
    space p = assembler.getSpace(mbP,1,2);

    // Incorporate Dirichlet BC
    bc1.setGeoMap(mp);
    bc2.setGeoMap(mp);
    bc0.setGeoMap(mp);
    u.setup(bc1, dirichlet::interpolation, 0);
    v.setup(bc2, dirichlet::interpolation, 0);
    p.setup(bc0, dirichlet::interpolation, 0);

    // Set the source term
    auto ff1 = assembler.getCoeff(f1, G);
    auto ff2 = assembler.getCoeff(f2, G);

    // Initialize the system
    assembler.initSystem();

    // Compute the system matrix and right-hand side
    assembler.assemble(
        igrad(u, G) * igrad(u, G).tr() * meas(G),
        igrad(v, G) * igrad(v, G).tr() * meas(G),
        igrad(u, G)[0] * p.tr() * meas(G),
        igrad(v, G)[1] * p.tr() * meas(G),
        p * igrad(u, G)[0].tr() * meas(G),
        p * igrad(v, G)[1].tr() * meas(G),
        u * ff1 * meas(G),
        v * ff2 * meas(G)
    );

    gsInfo << "done.\n";

    /************ Solve resulting linear system *************/

    index_t size = assembler.matrix().rows();
    gsInfo << size << " dofs.\n";
    gsSparseMatrix<> mat = assembler.matrix();
    mat(size-1,size-1)=1;
    
    gsInfo << mat << "\n\n";

    gsInfo << "Solve resulting linear system... " << std::flush;

    // Solve
    gsSparseSolver<>::LU solver( mat );
    gsMatrix<> sol = solver.solve( assembler.rhs() );

    gsInfo << "done.\n\n";

    // Print the solution
    gsInfo << "Solution:\n" << sol.transpose() << "\n\n";

    return EXIT_SUCCESS;
}
