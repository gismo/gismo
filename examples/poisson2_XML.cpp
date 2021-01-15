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
    bool verbose = false;
    index_t numRefine  = 0;
    index_t numElevate = 0;
    std::string fn("pde/poisson3d_bvp.xml");

    std::string save;
    std::string geometryRef;
    std::string geometry;

    gsCmdLine cmd("Tutorial on solving a Poisson problem.");
    cmd.addInt( "e", "degreeElevation",
                "Number of degree elevation steps to perform before solving (0: equalize degree in all directions)", numElevate );
    cmd.addInt( "r", "uniformRefine", "Number of Uniform h-refinement steps to perform before solving",  numRefine );
    cmd.addString( "f", "file", "Input XML file", fn );
    cmd.addSwitch("plot", "Create a ParaView visualization file with the solution", plot);
    cmd.addString("S","save", "Save solution to gismo XML file", save);
    cmd.addString("R","rgeom", "Save refined geometry to gismo XML file", geometryRef);
    cmd.addString("G","geom", "Save initial geometry to XML file", geometry);
    cmd.addSwitch("verbose","Verbose output", verbose);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }
    //! [Parse command line]

    //! [Read input file]

    gsFileData<> fd(fn);
    if (verbose)
        gsInfo << "Loaded file "<< fd.lastPath() <<"\n";

    index_t num;
    gsMultiPatch<> mp;
    num = fd.template count<gsMultiPatch<real_t>>(); // id=0: Multipatch domain
    GISMO_ENSURE(num==1,"Number of multipatch objects in XML should be 1, but is "<<num);
    fd.template getFirst<gsMultiPatch<real_t>>(mp); // Multipatch domain

    if (!geometry.empty())
    {
        if (verbose)
        {
            gsInfo<<"Writing geometry to XML...\n";
            gsInfo<<"Path = "<<geometry<<"\n";
        }
        gsWrite(mp,geometry);
    }

    gsFunctionExpr<> f;
    fd.getId(1, f); // id=1: source function
    if (verbose)
        gsInfo<<"Source function "<< f << "\n";

    gsBoundaryConditions<> bc;
    num = fd.template count<gsBoundaryConditions<>>(); // id=0: Multipatch domain
    GISMO_ENSURE(num==1,"Number of boundary condition objects in XML should be 1, but is "<<num);
    fd.template getFirst<gsBoundaryConditions<>>(bc); // Multipatch domain
    bc.setGeoMap(mp);
    if (verbose)
        gsInfo<<"Boundary conditions:\n"<< bc <<"\n";

    gsOptionList Aopt;
    num = fd.template count<gsOptionList>(); // id=0: Multipatch domain
    GISMO_ENSURE(num==1,"Number of options objects in XML should be 1, but is "<<num);
    fd.template getFirst<gsOptionList>(Aopt); // Multipatch domain

    //! [Read input file]

    //! [Refinement]
    gsMultiBasis<> dbasis(mp);

    // Elevate and p-refine the basis to order p + numElevate
    // where p is the highest degree in the bases
    mp.degreeElevate(numElevate);

    // h-refine each basis
    for (int r =0; r < numRefine-1; ++r)
        mp.uniformRefine();

    if (verbose)
    {
        gsInfo << "Patches: "<< mp.nPatches() <<", degree: "<< dbasis.minCwiseDegree() <<"\n";
        for (size_t k=0; k!=mp.nPatches(); k++)
        {
            gsInfo<<"------------------------------------------------------\n";
            gsInfo<<"Basis"<<k<<":\n";
            gsInfo<<mp.basis(k)<<"\n";
        }
        gsInfo<<"------------------------------------------------------\n";
    }
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
    //labels: Dirichlet, CornerValues, Collapsed, Clamped
    u.setup(bc, dirichlet::interpolation, 0);

    // Initialize the system
    A.initSystem(false);

    if (verbose)
        gsInfo<< "(dot1=assembled, dot2=solved)\n"
        "\nDoFs: "<<A.numDofs() <<std::flush;

    // Compute the system matrix and right-hand side
    A.assemble( igrad(u, G) * igrad(u, G).tr() * meas(G), u * ff * meas(G) );

    // Enforce Neumann conditions to right-hand side
    variable g_N = A.getBdrFunction();
    A.assembleRhsBc(u * g_N.val() * nv(G).norm(), bc.neumannSides() );
    //gsInfo<<"Sparse Matrix:\n"<< A.matrix().toDense() <<"\n";
    //gsInfo<<"Rhs vector:\n"<< A.rhs().transpose() <<"\n";

    if (verbose)
        gsInfo<< "." <<std::flush;// Assemblying done

    solver.compute( A.matrix() );
    solVector = solver.solve(A.rhs());

    if (verbose)
        gsInfo<< ".Done. \n" <<std::flush; // Linear solving done

    //! [Solver loop]
    if (!save.empty())
    {
        if (verbose)
        {
            gsInfo<<"Writing to XML...\n";
            gsInfo<<"Path = "<<save<<"\n";
        }
        gsMultiPatch<> mp_export;
        u_sol.extract(mp_export);
        gsWrite(mp_export,save);
    }
    if (!geometryRef.empty())
    {
        if (verbose)
        {
            gsInfo<<"Writing refined geometry to XML...\n";
            gsInfo<<"Path = "<<geometryRef<<"\n";
        }
        gsWrite(mp,geometryRef);
    }
    //! [Export visualization in ParaView]
    if (plot)
    {
        if (verbose)
            gsInfo<<"Plotting in Paraview...\n";

        ev.options().setSwitch("plot.elements", true);
        ev.writeParaview( u_sol   , G, "solution");
        //ev.writeParaview( u, G, "aa");

        gsFileManager::open("solution.pvd");
    }
    //! [Export visualization in ParaView]

    return EXIT_SUCCESS;

}// end main
