/** @file poisson2_example.cpp

    @brief Tutorial on how to use expression assembler to solve the Poisson equation

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris
*/

//! [Include namespace]
#include <gismo.h>

using namespace gismo;
//! [Include namespace]

int main(int argc, char *argv[])
{
    //! [Parse command line]
    bool plot = false;
    bool save = false;
    index_t numRefine  = 0;
    index_t numElevate = 0;
    bool last = false;
    std::string fn("pde/poisson2d_bvp.xml");

    gsCmdLine cmd("Tutorial on solving a Poisson problem.");
    cmd.addInt( "e", "degreeElevation",
                "Number of degree elevation steps to perform before solving (0: equalize degree in all directions)", numElevate );
    cmd.addInt( "r", "uniformRefine", "Number of Uniform h-refinement steps to perform before solving",  numRefine );
    cmd.addString( "f", "file", "Input XML file", fn );
    cmd.addSwitch("plot", "Create a ParaView visualization file with the solution", plot);
    cmd.addSwitch("save", "Save gismo XML file", save);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }
    //! [Parse command line]

    //! [Read input file]

    gsFileData<> fd(fn);
    gsInfo << "Loaded file "<< fd.lastPath() <<"\n";

    gsMultiPatch<> mp;
    fd.getId(0, mp); // id=0: Multipatch domain

    gsFunctionExpr<> f;
    fd.getId(1, f); // id=1: source function
    gsInfo<<"Source function "<< f << "\n";

    gsBoundaryConditions<> bc;
    fd.getId(2, bc); // id=2: boundary conditions
    bc.setGeoMap(mp);
    gsInfo<<"Boundary conditions:\n"<< bc <<"\n";

    gsOptionList Aopt;
    fd.getId(4, Aopt); // id=4: assembler options

    //! [Read input file]

    //! [Refinement]
    gsMultiBasis<> dbasis(mp);

    // Elevate and p-refine the basis to order p + numElevate
    // where p is the highest degree in the bases
    dbasis.setDegree( dbasis.maxCwiseDegree() + numElevate);

    // h-refine each basis
    for (int r =0; r < numRefine-1; ++r)
        dbasis.uniformRefine();

    gsInfo << "Patches: "<< mp.nPatches() <<", degree: "<< dbasis.minCwiseDegree() <<"\n";
    //! [Refinement]

    //! [Problem setup]
    gsExprAssembler<> A(1,1);
    A.setOptions(Aopt);

    //gsInfo<<"Active options:\n"<< A.options() <<"\n";
    typedef gsExprAssembler<>::geometryMap geometryMap;
    typedef gsExprAssembler<>::variable    variable;
    typedef gsExprAssembler<>::space       space;
    typedef gsExprAssembler<>::solution    solution;

    // Elements used for numerical integration
    A.setIntegrationElements(dbasis);
    gsExprEvaluator<> ev(A);

    // Set the geometry map
    geometryMap G = A.getMap(mp);

    // Set the discretization space
    space u = A.getSpace(dbasis);

    //next two steps are now moved to setup(.)
    // u.setInterfaceCont(0);
    // u.addBc( bc.get("Dirichlet") );

    // Set the source term
    variable ff = A.getCoeff(f, G);

    // Solution vector and solution variable
    gsMatrix<> solVector;
    solution u_sol = A.getSolution(u, solVector);

    gsSparseSolver<>::CGDiagonal solver;

    //! [Problem setup]

    //! [Solver loop]
    gsInfo<< "(dot1=assembled, dot2=solved)\n"
        "\nDoFs: ";

    //labels: Dirichlet, CornerValues, Collapsed, Clamped
    u.setup(bc, dirichlet::interpolation, 0);

    // Initialize the system
    A.initSystem(false);

    gsInfo<< A.numDofs() <<std::flush;

    // Compute the system matrix and right-hand side
    A.assemble( igrad(u, G) * igrad(u, G).tr() * meas(G), u * ff * meas(G) );

    // Enforce Neumann conditions to right-hand side
    variable g_N = A.getBdrFunction();
    A.assembleRhsBc(u * g_N.val() * nv(G).norm(), bc.neumannSides() );
    //gsInfo<<"Sparse Matrix:\n"<< A.matrix().toDense() <<"\n";
    //gsInfo<<"Rhs vector:\n"<< A.rhs().transpose() <<"\n";

    gsInfo<< "." <<std::flush;// Assemblying done

    solver.compute( A.matrix() );
    solVector = solver.solve(A.rhs());

    gsInfo<< ".Done. \n" <<std::flush; // Linear solving done

    //! [Solver loop]
    if (save)
    {
        gsInfo<<"Writing to XML...\n";
        gsMultiPatch<> mp_export;
        u_sol.extract(mp_export);
        gsWrite(mp,"geometry");
        gsWrite(mp_export,"solution");
    }


    //! [Export visualization in ParaView]
    if (plot)
    {
        gsInfo<<"Plotting in Paraview...\n";
        ev.options().setSwitch("plot.elements", true);
        ev.writeParaview( u_sol   , G, "solution");
        //ev.writeParaview( u, G, "aa");

        gsFileManager::open("solution.pvd");
    }
    else
        gsInfo << "No output created, re-run with --plot to get a ParaView "
                  "file containing the solution.\n";
    //! [Export visualization in ParaView]

    return EXIT_SUCCESS;

}// end main
