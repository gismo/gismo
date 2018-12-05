/** @file multiGrid_example.cpp

    @brief Provides an examples for the multigrid solver.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): C. Hofreither, S. Takacs
*/

#include <gismo.h>

using namespace gismo;

int main(int argc, char *argv[])
{
    /************** Define command line options *************/

    std::string geometry("domain2d/yeti_mp2.xml");
    index_t refinements = 3;
    index_t degree = 2;
    bool dg = false;
    index_t levels = -1;
    index_t cycles = 1;
    index_t presmooth = 1;
    index_t postsmooth = 1;
    std::string smoother("GaussSeidel");
    real_t damping = -1;
    real_t scaling = 0.12;
    real_t tolerance = 1.e-8;
    index_t maxIterations = 100;
    bool plot = false;
    std::string boundary_conditions("d");

    gsCmdLine cmd("Solves a PDE with an isogeometric discretization using a multigrid solver.");
    cmd.addString("g", "Geometry",              "Geometry file", geometry);
    cmd.addInt   ("r", "Refinements",           "Number of uniform h-refinement steps to perform before solving", refinements);
    cmd.addInt   ("p", "Degree",                "Degree of the B-spline discretization space", degree);
    cmd.addSwitch("",  "DG",                    "Use a discontinuous Galerkin discretization", dg);
    cmd.addInt   ("l", "MG.Levels",             "Number of levels to use for multigrid iteration", levels);
    cmd.addInt   ("c", "MG.Cycles",             "Number of multi-grid cycles", cycles);
    cmd.addInt   ("",  "MG.Presmooth",          "Number of pre-smoothing steps", presmooth);
    cmd.addInt   ("",  "MG.Postsmooth",         "Number of post-smoothing steps", postsmooth);
    cmd.addString("s", "MG.Smoother",           "Smoothing method", smoother);
    cmd.addReal  ("",  "MG.Damping",            "Damping factor for the smoother", damping);
    cmd.addReal  ("",  "MG.Scaling",            "Scaling factor for the subspace corrected mass smoother", scaling);
    cmd.addReal  ("t", "CG.Tolerance",          "Stopping criterion for cg", tolerance);
    cmd.addInt   ("",  "CG.MaxIterations",      "Stopping criterion for cg", maxIterations);
    cmd.addString("b", "BoundaryConditions",    "Boundary conditions", boundary_conditions);
    cmd.addSwitch("",  "Plot",                  "Plot the result with Paraview", plot);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    gsOptionList opt = cmd.getOptionList();

    // Default case is levels:=refinements, so replace invalid default accordingly
    if (levels <0) { levels = refinements; opt.setInt( "MG.Levels", levels ); }
    // The smoothers know their defaults, so remove the invalid default
    if (damping<0) { opt.remove( "MG.Damping" ); }

    // Define assembler options
    opt.remove( "DG" );
    opt.addInt( "Ass.InterfaceStrategy", "", (index_t)( dg ? iFace::dg : iFace::conforming )  );
    opt.addInt( "Ass.DirichletStrategy", "", (index_t) dirichlet::elimination                 );

    if ( ! gsFileManager::fileExists(geometry) )
    {
        gsInfo << "Geometry file could not be found.\n";
        gsInfo << "I was searching in the current directory and in: " << gsFileManager::getSearchPaths() << "\n";
        return EXIT_FAILURE;
    }

    gsInfo << "Run multiGrid_example with options:\n" << opt << std::endl;

    /******************* Define geometry ********************/

    gsInfo << "Define geometry... " << std::flush;

    gsMultiPatch<>::uPtr mpPtr = gsReadFile<>(geometry);
    if (!mpPtr)
    {
        gsInfo << "No geometry found in file " << geometry << ".\n";
        return EXIT_FAILURE;
    }
    gsMultiPatch<>& mp = *mpPtr;

    gsInfo << "done.\n";

    /************** Define boundary conditions **************/

    gsInfo << "Define boundary conditions... " << std::flush;

    gsConstantFunction<> one(1.0, mp.geoDim());

    gsBoundaryConditions<> bc;
    {
        const index_t len = boundary_conditions.length();
        index_t i = 0;
        for (gsMultiPatch<>::const_biterator it = mp.bBegin(); it < mp.bEnd(); ++it)
        {
            char b_local;
            if ( len == 1 )
                b_local = boundary_conditions[0];
            else if ( i < len )
                b_local = boundary_conditions[i];
            else
            {
                gsInfo << "\nNot enough boundary conditions given.\n";
                return EXIT_FAILURE;
            }

            if ( b_local == 'd' )
                bc.addCondition( *it, condition_type::dirichlet, &one );
            else if ( b_local == 'n' )
                bc.addCondition( *it, condition_type::neumann, &one );
            else
            {
                gsInfo << "\nInvalid boundary condition given; only 'd' (Dirichlet) and 'n' (Neumann) are supported.\n";
                return EXIT_FAILURE;
            }

            ++i;
        }
        if ( len > i )
            gsInfo << "\nToo much boundary conditions have been specified. Ingnoring the remaining ones.\n";
        gsInfo << "done. "<<i<<" boundary conditions set.\n";
    }


    /************ Setup bases and adjust degree *************/

    gsInfo << "Setup bases and adjust degree... " << std::flush;

    gsMultiBasis<> mb(mp);

    for ( size_t i = 0; i < mb.nBases(); ++ i )
        mb[i].setDegreePreservingMultiplicity(degree);

    for ( index_t i = 0; i < refinements; ++i )
        mb.uniformRefine();

    gsInfo << "done.\n";

    /********* Setup assembler and assemble matrix **********/

    gsInfo << "Setup assembler and assemble matrix... " << std::flush;

    gsPoissonAssembler<> assembler(
        mp,
        mb,
        bc,
        gsConstantFunction<>(1,mp.geoDim()),
        (dirichlet::strategy) opt.getInt("Ass.DirichletStrategy"),
        (iFace::strategy)     opt.getInt("Ass.InterfaceStrategy")
    );
    assembler.assemble();

    gsInfo << "done.\n";

    /**************** Setup solver and solve ****************/

    gsInfo << "Setup solver and solve... " << std::flush;

    std::vector< gsSparseMatrix<real_t,RowMajor> > transferMatrices;
    std::vector< gsMultiBasis<real_t> > multiBases; // only needed for subspace corrected mass smoother

    gsGridHierarchy<>::buildByCoarsening(give(mb), bc, opt.getGroup("MG"))
        .moveMultiBasesTo(multiBases)
        .moveTransferMatricesTo(transferMatrices)
        .clear();

    gsMultiGridOp<>::Ptr mg = gsMultiGridOp<>::make( assembler.matrix(), transferMatrices );
    mg->setOptions( opt.getGroup("MG") );

    for (index_t i = 1; i < mg->numLevels(); ++i)
    {
        gsPreconditionerOp<>::Ptr smootherOp;
        if ( opt.getString("MG.Smoother") == "Richardson" )
            smootherOp = makeRichardsonOp(mg->matrix(i));
        else if ( opt.getString("MG.Smoother") == "Jacobi" )
            smootherOp = makeJacobiOp(mg->matrix(i));
        else if ( opt.getString("MG.Smoother") == "GaussSeidel" )
            smootherOp = makeGaussSeidelOp(mg->matrix(i));
        else if ( opt.getString("MG.Smoother") == "SubspaceCorrectedMassSmoother" )
        {
            if (multiBases[i].nBases() != 1)
            {
                gsInfo << "The chosen smoother only works for single-patch geometries.\n\n";
                return EXIT_FAILURE;
            }
            smootherOp = gsPreconditionerFromOp<>::make(
                mg->underlyingOp(i),
                gsPatchPreconditionersCreator<>::subspaceCorrectedMassSmootherOp(multiBases[i][0],bc,opt.getGroup("Ass"),scaling)
            );
        }
        else
        {
            gsInfo << "The chosen smoother is unknown.\n\nKnown are:\n  Richardson\n  Jacobi\n  GaussSeidel"
                      "\n  SubspaceCorrectedMassSmoother\n\n";
            return EXIT_FAILURE;
        }
        smootherOp->setOptions( opt.getGroup("MG") );
        mg->setSmoother(i, smootherOp);
    }

    gsMatrix<> x, errorHistory;
    x.setRandom( assembler.matrix().rows(), 1 );

    gsConjugateGradient<>( assembler.matrix(), mg )
        .setOptions( opt.getGroup("CG") )
        .solveDetailed( assembler.rhs(), x, errorHistory );

    gsInfo << "done.\n\n";

    /******************** Print end Exit ********************/

    const index_t iter = errorHistory.rows()-1;
    const bool success = errorHistory(iter,0) < tolerance;
    if (success)
        gsInfo << "Reached desired tolerance after " << iter << " iterates:\n";
    else
        gsInfo << "Did not reach desired tolerance after " << iter << " iterates:\n";

    if (errorHistory.rows() < 20)
        gsInfo << errorHistory.transpose() << "\n\n";
    else
        gsInfo << errorHistory.topRows(5).transpose() << " ... " << errorHistory.bottomRows(5).transpose()  << "\n\n";

    if (plot)
    {
        // Construct the solution as a scalar field
        gsMultiPatch<> mpsol;
        assembler.constructSolution(x, mpsol);
        gsField<> sol( assembler.patches(), mpsol );

        // Write approximate and exact solution to paraview files
        gsInfo << "Plotting in Paraview.\n";
        gsWriteParaview<>(sol, "multiGrid_result", 1000);
        gsFileManager::open("multiGrid_result.pvd");
    }
    else
    {
        gsInfo << "Done. No output created, re-run with --Plot to get a ParaView "
                  "file containing the solution.\n";
    }
    return success ? EXIT_SUCCESS : EXIT_FAILURE;
};
