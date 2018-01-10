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
    real_t tolerance = 1.e-8;
    index_t maxIterations = 100;

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
    cmd.addReal  ("",  "MG.Damping",            "Damping factor for the smoother (handed over to smoother)", damping);
    cmd.addReal  ("t", "CG.Tolerance",          "Stopping criterion for cg", tolerance);
    cmd.addInt   ("",  "CG.MaxIterations",      "Stopping criterion for cg", maxIterations);

    cmd.getValues(argc,argv);

    gsOptionList opt = cmd.getOptionList();

    // Handle some non-trivial standards
    if (levels <0)  { levels = refinements; opt.setInt( "MG.Levels", refinements );                             }
    if (damping<0)  { opt.remove( "MG.Damping" );                                                               }
    if (dg)         { opt.addInt( "MG.InterfaceStrategy", "", (index_t)iFace::dg         ); opt.remove( "DG" ); }
    else            { opt.addInt( "MG.InterfaceStrategy", "", (index_t)iFace::conforming ); opt.remove( "DG" ); }

    if ( ! gsFileManager::fileExists(geometry) )
    {
        gsInfo << "Geometry file could not be found.\n";
        gsInfo << "I was searching in the current directory and in: " << gsFileManager::getSearchPaths() << "\n";
        return EXIT_FAILURE;
    }

    gsInfo << "Run gsMultiGridTutorial with options:\n" << opt << std::endl;

    /******************* Define geometry ********************/

    gsInfo << "Define geometry... " << std::flush;

    gsMultiPatch<> mp;

    gsFileData<> fileData(geometry);
    if (!fileData.has< gsMultiPatch<> >())
    {
        gsInfo << "No multipatch object found in file " << geometry << ".\n";
        return EXIT_FAILURE;
    }
    fileData.getFirst< gsMultiPatch<> >(mp);

    gsInfo << "done.\n";

    /************** Define boundary conditions **************/

    gsInfo << "Define boundary conditions... " << std::flush;

    gsConstantFunction<> one(1.0, mp.geoDim());

    gsBoundaryConditions<> bc;
    bc.addCondition( boundary::west,  condition_type::dirichlet, &one );
    bc.addCondition( boundary::east,  condition_type::dirichlet, &one );
    if (mp.geoDim() >= 2)
    {
        bc.addCondition( boundary::south, condition_type::dirichlet, &one );
        bc.addCondition( boundary::north, condition_type::dirichlet, &one );
    }
    if (mp.geoDim() >= 3)
    {
        bc.addCondition( boundary::front, condition_type::dirichlet, &one );
        bc.addCondition( boundary::back,  condition_type::dirichlet, &one );
    }

    gsInfo << "done.\n";

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
        gsConstantFunction<>(0,mp.geoDim()),
        (dirichlet::strategy) opt.askInt("MG.DirichletStrategy",dirichlet::elimination),
        (iFace::strategy) opt.askInt("MG.InterfaceStrategy",iFace::conforming)
    );
    assembler.assemble();

    gsInfo << "done.\n";

    /**************** Setup solver and solve ****************/

    gsInfo << "Setup solver and solve... " << std::flush;

    std::vector< gsSparseMatrix<real_t,RowMajor> > transferMatrices;

    gsGridHierarchy<>::buildByCoarsening(give(mb), bc, opt.getGroup("MG"))
        .moveTransferMatricesTo(transferMatrices)
        .clear();

    gsMultiGridOp<>::Ptr mg = gsMultiGridOp<>::make( assembler.matrix(), transferMatrices );
    mg->setOptions( opt.getGroup("MG") );

    for (index_t i = 1; i < mg->numLevels(); ++i)
    {
        gsPreconditionerOp<>::Ptr smootherOp;
        if ( opt.getString("MG.Smoother") == "Richardson" )
            smootherOp = makeRichardsonOp(mg->matrix(i),(real_t)1/2);
        else if ( opt.getString("MG.Smoother") == "Jacobi" )
            smootherOp = makeJacobiOp(mg->matrix(i),(real_t)1/2);
        else if ( opt.getString("MG.Smoother") == "GaussSeidel" )
            smootherOp = makeGaussSeidelOp(mg->matrix(i));
        else
        {
            gsInfo << "The chosen smoother is unknown.\n\n";
            return EXIT_FAILURE;
        }
        smootherOp->setOptions( opt.getGroup("MG") );
        mg->setSmoother(i, smootherOp);
    }

    gsMatrix<> x, history;
    x.setRandom( assembler.matrix().rows(), 1 );

    gsConjugateGradient<>( assembler.matrix(), mg )
        .setOptions( opt.getGroup("CG") )
        .solveDetailed( assembler.rhs(), x, history );

    gsInfo << "done.\n";

    /******************** Print end Exit ********************/

    gsInfo << "\n\nResidual:\n" << history << "\n\n";

    if (history(history.rows()-1,0) < tolerance)
        gsInfo << "Reached desired tolerance after " << history.rows()-1 << " iterates.\n";
    else
        gsInfo << "Did not reach desired tolerance after " << history.rows()-1 << " iterates.\n";


    return EXIT_SUCCESS;
};
