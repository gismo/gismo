/** @file poisson2_example.cpp

    @brief Poisson example with command line arguments.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris, J. Sogn
*/

# include <gismo.h>

using namespace gismo;

int main(int argc, char *argv[])
{

    /////////////////// Input ///////////////////
    std::string fn_pde("");
    bool Nitsche  = false;
    bool dG = false;
    bool plot = false;
    index_t plot_pts = 1000;
    index_t numElevate = -1;
    index_t numRefine = 2;
    std::string fn_basis("");
    std::string fn("");

    gsCmdLine cmd( "Solves Poisson's equation with an isogeometric discretization." );
    cmd.addString( "p", "pde", "File containing a poisson PDE (.xml)", fn_pde );
    cmd.addSwitch( "nitsche", "Use the Nitsche's method for Dirichlet sides", Nitsche );
    cmd.addSwitch( "discGalerkin", "Use Discontiouous Galerkin method for patch interfaces", dG );
    cmd.addSwitch( "plot", "Plot result in ParaView format", plot );
    cmd.addInt( "s", "plotSamples", "Number of sample points to use for plotting", plot_pts );
    cmd.addInt( "e", "degreeElevation",
                "Number of degree elevation steps to perform before solving (0: equalize degree in all directions)",
                numElevate );
    cmd.addInt( "r", "uniformRefine", "Number of Uniform h-refinement steps to perform before solving",
                numRefine );
    cmd.addString( "b", "basis", "File containing basis for discretization (.xml)", fn_basis );
    cmd.addString( "g", "geometry", "File containing Geometry (.xml, .axl, .txt)", fn );
    cmd.getValues(argc,argv);

    gsMultiBasis<> bases;
    if ( ! fn_basis.empty() )
    {
        gsBasis<>::uPtr bb = gsReadFile<>( fn_basis );
        gsInfo << "Got basis: "<< * bb<<"\n";
        bases.addBasis(bb.release());
    }

    if (numRefine<0)
    {
        gsInfo << "Number of refinements must be non-negative.\n";
        return -1;
    }

    if (numElevate<-1)
    {
        gsInfo << "Number of elevations must be non-negative.\n";
        return -1;
    }

    gsMultiPatch<>::uPtr geo;
    if ( fn_pde.empty() )
    {
        fn_pde = GISMO_DATA_DIR;
        if ( !fn.empty() )
        {
            geo = gsReadFile<>( fn );
            if ( !geo )
            {
                gsWarn<< "Did not find any geometry in "<<fn<<", quitting.\n";
                return false;
            }
            switch ( geo->geoDim() )
            {
            case 1:
                fn_pde+="/pde/poisson1d_sin.xml";
                break;
            case 2:
                fn_pde+="/pde/poisson2d_sin.xml";
                break;
            case 3:
                fn_pde+="/pde/poisson3d_sin.xml";
                break;
            default:
                return -1;
            }
        }
        else
        fn_pde+="/pde/poisson2d_sin.xml";
    }

    memory::unique_ptr< gsPoissonPde<> > ppde = gsReadFile<>(fn_pde);

    gsFunctionExpr<>::uPtr exactSol = gsReadFile<>(fn_pde, 100);
    if ( !ppde )
    {
        gsWarn<< "Did not find any PDE in "<< fn<<", quitting.\n";
        return -1;
    }

    if ( fn.empty() )
    {
        fn = GISMO_DATA_DIR;
        switch ( ppde->m_compat_dim )
        {
        case 1:
            fn+= "domain1d/bspline1d_01.xml";
            break;
        case 2:
            fn+= "domain2d/square.xml";
            break;
        case 3:
            fn+= "domain3d/cube.xml";
            break;
        default:
            return false;
        }
    }

    geo = gsReadFile<>( fn );
    if ( !geo )
    {
        gsInfo << "Did not find any geometries in "<< fn<<", quitting.\n";
        return -1;
    }

    /////////////////// Print info ///////////////////
    gsInfo<<"Domain: "<< *geo <<"\n";
    gsInfo<<"Number of patches are " << geo->nPatches() << "\n";
    gsInfo<<"Source function "<< *ppde->rhs() << "\n";
    gsInfo<<"Exact solution "<< *exactSol <<".\n" << "\n";
    gsInfo<<"p-refinent steps before solving: "<< numElevate <<"\n";
    gsInfo<<"h-refinent steps before solving: "<< numRefine <<"\n";

    gsInfo<< * ppde <<"\n";

    /////////////////// Setup boundary conditions ///////////////////
    // Define Boundary conditions
    gsBoundaryConditions<> bcInfo;

    ppde->patches() = *geo;

    // Create Dirichlet boundary conditions for all boundaries
    for (gsMultiPatch<>::const_biterator
         bit = geo->bBegin(); bit != geo->bEnd(); ++bit)
    {
        bcInfo.addCondition( *bit, condition_type::dirichlet, exactSol.get() );
    }

    ppde->boundaryConditions() = bcInfo;

    /////////////////// Refinement h and p ///////////////////
    if ( bases.nBases() == 0 )
        bases = gsMultiBasis<>(*geo);

    // Elevate and refine the solution space
    if ( numElevate > -1 )
    {
        // get maximum degree
        int tmp = bases.maxDegree(0);

        // Elevate all degrees uniformly
        tmp += numElevate;
        for (size_t j = 0; j < bases.nBases(); ++j )
                bases[j].setDegree(tmp);
    }

    // Refining the basis
    for (size_t j = 0; j < bases.nBases(); ++j )
        for (int i = 0; i < numRefine; ++i)
            bases[j].uniformRefine();

    gsInfo << "Discrete. Space 0: "<< bases[0] << "\n";


    /////////////////// Setup solver ///////////////////
    //Initialize Solver

    gsPoissonAssembler<real_t> PoissonAssembler;

    gsOptionList options = PoissonAssembler.defaultOptions();
    //Use Nitsche's method for Dirichlet boundaries
    if ( Nitsche )
    {
        gsInfo<<"Using Nitsche's method for Dirichlet boundaries.\n";
        options.setInt("DirichletStrategy", dirichlet::nitsche);
    }

    if ( dG )
    {
        gsInfo<<"Using dG method for patch interfaces.\n";
        options.setInt("InterfaceStrategy", iFace::dg);
    }

    PoissonAssembler.initialize(*ppde, bases, options);

    // Generate system matrix and load vector
    gsInfo<<"Assembling...\n";
    PoissonAssembler.assemble();

    // gsDebugVar(PoissonAssembler.matrix().rows());
    // gsDebugVar(PoissonAssembler.matrix().cols());
    // gsDebugVar(PoissonAssembler.rhs().rows());

    // Initialize the conjugate gradient solver
    gsInfo<<"Solving...\n";
    gsSparseSolver<>::CGDiagonal solver( PoissonAssembler.matrix() );
    gsMatrix<> solVector = solver.solve( PoissonAssembler.rhs()    );

    // Construct the solution as a scalar field
    gsMultiPatch<> mpsol;
    PoissonAssembler.constructSolution(solVector, mpsol);
    gsField<> sol( PoissonAssembler.patches(), mpsol);

    gsInfo <<"Sol:"<< mpsol <<"\n";

    // Plot solution in paraview
    int result = 0;
    if (plot)
    {
        // Write approximate and exact solution to paraview files
        gsInfo<<"Plotting in Paraview...\n";
        gsWriteParaview<>(sol, "poisson2d", plot_pts);
        gsField<> exact(*geo, *exactSol, false );
        gsWriteParaview<>(exact, "poisson2d_exact", plot_pts);

        // Run paraview
        result = system("paraview poisson2d.pvd &");
    }

    gsInfo << "Test is done: Exiting" << "\n";

    return result;
}
