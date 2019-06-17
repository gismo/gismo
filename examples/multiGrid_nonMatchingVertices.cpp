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

gsPreconditionerOp<>::Ptr setupSubspaceCorrectedMassSmoother(const gsSparseMatrix<>&, const gsMultiBasis<>&, const gsMultiPatch<>&,
                                                             const gsBoundaryConditions<>&, const gsOptionList&, bool useGeo, bool autoDamping);

int main(int argc, char *argv[])
{
    /************** Define command line options *************/

    std::string geometry("domain2d/yeti_mp2.xml");
    index_t refinements = 3;
    index_t degree = 2;
    bool nonMatching = false;
    bool dg = false;
    index_t levels = -1;
    index_t cycles = 1;
    index_t presmooth = 1;
    index_t postsmooth = 1;
    bool extrasmooth = false;
    std::string smoother("GaussSeidel");
    real_t damping = -1;
    bool autoDamping = false;
    real_t scaling = 0.12;
    std::string iterativeSolver("cg");
    real_t tolerance = 1.e-8;
    index_t maxIterations = 100;
    bool plot = false;
    bool doNURBS = false;
    std::string boundary_conditions("d");
    index_t c = 1;

    gsCmdLine cmd("Solves a PDE with an isogeometric discretization using a multigrid solver.");
    cmd.addString("g", "Geometry",              "Geometry file", geometry);
    cmd.addInt   ("r", "Refinements",           "Number of uniform h-refinement steps to perform before solving", refinements);
    cmd.addInt   ("p", "Degree",                "Degree of the B-spline discretization space", degree);
    cmd.addSwitch("",  "NonMatching",           "Set up a non-matching multi-patch discretization", nonMatching);
    cmd.addSwitch("",  "DG",                    "Use a discontinuous Galerkin discretization", dg);
    cmd.addInt   ("l", "MG.Levels",             "Number of levels to use for multigrid iteration", levels);
    cmd.addInt   ("c", "MG.Cycles",             "Number of multi-grid cycles", cycles);
    cmd.addInt   ("",  "MG.Presmooth",          "Number of pre-smoothing steps", presmooth);
    cmd.addInt   ("",  "MG.Postsmooth",         "Number of post-smoothing steps", postsmooth);
    cmd.addSwitch("",  "MG.Extrasmooth",        "Doubles the number of smoothing steps for each coarser level", extrasmooth);
    cmd.addString("s", "MG.Smoother",           "Smoothing method", smoother);
    cmd.addReal  ("",  "MG.Damping",            "Damping factor for the smoother", damping);
    cmd.addSwitch(     "MG.AutoDamping",        "Adjeusts the damping of the subspace corrected mass smoother such that "
                                                "rho(tau P*A) == damping.", autoDamping);
    cmd.addReal  ("",  "MG.Scaling",            "Scaling factor for the subspace corrected mass smoother", scaling);
    cmd.addString("i", "IterativeSolver",       "Iterative solver: apply multigrid directly (d) or as a preconditioner for conjugate gradient (cg)", iterativeSolver);
    cmd.addReal  ("t", "Solver.Tolerance",      "Stopping criterion for linear solver", tolerance);
    cmd.addInt   ("",  "Solver.MaxIterations",  "Stopping criterion for linear solver", maxIterations);
    cmd.addString("b", "BoundaryConditions",    "Boundary conditions", boundary_conditions);
    cmd.addSwitch("",  "Plot",                  "Plot the result with Paraview", plot);
    cmd.addSwitch("",  "NURBS",                 "Use NURBS as Ansatz functions or BSplines", doNURBS);
    cmd.addInt   ("",  "Case",                  "Choose a test case", c);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    gsOptionList opt = cmd.getOptionList();

    // Default case is levels:=refinements, so replace invalid default accordingly
    if (levels <0) { levels = refinements; opt.setInt( "MG.Levels", levels ); }
    // The smoothers know their defaults, so remove the invalid default
    if (damping<0) { if (autoDamping) opt.setReal("MG.Damping", 1); else opt.remove( "MG.Damping" ); }

    // Define assembler options
    opt.remove( "DG" );
    opt.addInt( "MG.InterfaceStrategy", "", (index_t)( dg ? iFace::dg : iFace::conforming )  );
    opt.addInt( "MG.DirichletStrategy", "", (index_t) dirichlet::elimination                 );

    if(c == 1)
        if ( ! gsFileManager::fileExists(geometry) )
        {
            gsInfo << "Geometry file could not be found.\n";
            gsInfo << "I was searching in the current directory and in: " << gsFileManager::getSearchPaths() << "\n";
            return EXIT_FAILURE;
        }

    gsInfo << "Run multiGrid_example with options:\n" << opt << std::endl;

    if (autoDamping && ( smoother == "Richardson" || smoother == "r"  || smoother == "Jacobi" || smoother == "j"
        || smoother == "GaussSeidel" || smoother == "gs" ) )
    {
        gsInfo << "--MG.AutoDamping is only available for the mass smoothers.\n";
        return EXIT_FAILURE;
    }

    /******************* Define geometry ********************/

    gsInfo << "Define geometry... " << std::flush;

    gsMultiPatch<>::uPtr mpPtr;

    switch (c)
    {
        case 1:
        {
            mpPtr = gsReadFile<>(geometry);
            if (!mpPtr)
            {
                gsInfo << "No geometry found in file " << geometry << ".\n";
                return EXIT_FAILURE;
            }
            break;
        }
        case 2:
        {
            mpPtr = memory::unique_ptr<gsMultiPatch<> >(  new gsMultiPatch<>(gsNurbsCreator<>::BSplineSquareGrid(2, 2, 0.5)) );
            //mpPtr->computeTopology();
            break;
        }
        case 3:
        {
            mpPtr = memory::unique_ptr<gsMultiPatch<> >(  new gsMultiPatch<>(gsNurbsCreator<>::BSplineSquareGrid(2, 1, 0.5)) );
            gsNurbsCreator<>::TensorBSpline2Ptr temp2 = (gsNurbsCreator<>::BSplineRectangle(0, 0.5, 1, 1));

            //mpPtr->addPatch(give(temp1[0]));
            //mpPtr->addPatch(give(temp1[1]));
            mpPtr->addPatch(give(temp2));

            mpPtr->clearTopology();
            mpPtr->addInterface(0, 2, 1, 1);
            mpPtr->addInterface(2, 3, 0, 4);
            mpPtr->addInterface(2, 3, 1, 4);

            mpPtr->addBoundary(0, 3);
            mpPtr->addBoundary(0, 1);
            mpPtr->addBoundary(1, 3);
            mpPtr->addBoundary(1, 2);
            mpPtr->addBoundary(2, 4);
            mpPtr->addBoundary(2, 1);
            mpPtr->addBoundary(2, 2);
            break;
        }
        case 4:
        {
            mpPtr = memory::unique_ptr<gsMultiPatch<> >(  new gsMultiPatch<>(gsNurbsCreator<>::BSplineSquareGrid(2, 1, 3)) );
            gsMultiPatch<> temp2 = (gsNurbsCreator<>::BSplineSquareGrid(3, 1, 2, 0, 3));

            mpPtr->addPatch(temp2[0]);
            mpPtr->addPatch(temp2[1]);
            mpPtr->addPatch(temp2[2]);

            mpPtr->clearTopology();
            mpPtr->addInterface(0, 2, 1, 1);
            mpPtr->addInterface(2, 2, 3, 1);
            mpPtr->addInterface(3, 2, 4, 1);
            mpPtr->addInterface(2, 3, 0, 4);
            mpPtr->addInterface(3, 3, 0, 4);
            mpPtr->addInterface(3, 3, 1, 4);
            mpPtr->addInterface(4, 3, 1, 4);

            mpPtr->addBoundary(0, 3);
            mpPtr->addBoundary(0, 1);
            mpPtr->addBoundary(1, 3);
            mpPtr->addBoundary(1, 2);
            mpPtr->addBoundary(2, 4);
            mpPtr->addBoundary(2, 1);
            mpPtr->addBoundary(3, 4);
            mpPtr->addBoundary(4, 4);
            mpPtr->addBoundary(4, 2);
            break;
        }
        case 5: // circular ring
        {
            real_t phi = 10.;
            mpPtr = memory::unique_ptr<gsMultiPatch<> >( new gsMultiPatch<>() );
            gsGeometry<real_t>::Ptr geo1 = ( gsNurbsCreator<>::NurbsQuarterAnnulus(1.0, 1.5) );
            gsGeometry<real_t>::Ptr geo2 = ( gsNurbsCreator<>::NurbsQuarterAnnulus(1.5, 2.) );

            mpPtr->addPatch( *gsNurbsCreator<>::NurbsQuarterAnnulus(1.0, 1.5) );
            mpPtr->addPatch( *gsNurbsCreator<>::NurbsQuarterAnnulus(1.5, 2. ) );

            size_t size = mpPtr->nPatches();

            for(size_t rot = 1; rot < 4; rot++)
            {
                for (size_t i = 0; i < size; i++)
                {
                    gsGeometry<real_t>::Ptr dub1 = mpPtr->patch(i).clone();
                    dub1->rotate(rot * M_PI / 2.);
                    mpPtr->addPatch(*dub1);
                }

            }
            //gsInfo << "number of patches: " << mpPtr->detail() << "\n";

            mpPtr->addInterface(0, 4, 2, 3);
            mpPtr->addInterface(0, 2, 1, 1);
            mpPtr->addInterface(0, 3, 6, 4);
            mpPtr->addInterface(1, 3, 7, 4);
            mpPtr->addInterface(1, 4, 3, 3);
            mpPtr->addInterface(2, 2, 3, 1);
            mpPtr->addInterface(2, 4, 4, 3);
            mpPtr->addInterface(3, 4, 5, 3);
            mpPtr->addInterface(4, 2, 5, 1);
            mpPtr->addInterface(4, 4, 6, 3);
            mpPtr->addInterface(5, 4, 7, 3);
            mpPtr->addInterface(6, 2, 7, 1);

            mpPtr->addBoundary(0, 1);
            mpPtr->addBoundary(1, 2);
            mpPtr->addBoundary(2, 1);
            mpPtr->addBoundary(3, 2);
            mpPtr->addBoundary(4, 1);
            mpPtr->addBoundary(5, 2);
            mpPtr->addBoundary(6, 1);
            mpPtr->addBoundary(7, 2);

            // additionally rotate the inner ring
            mpPtr->patch(0).rotate(phi*M_PI/180.);
            mpPtr->patch(2).rotate(phi*M_PI/180.);
            mpPtr->patch(4).rotate(phi*M_PI/180.);
            mpPtr->patch(6).rotate(phi*M_PI/180.);

            // need a few additional interfaces when we rotate the inner ring
            if(phi != 0)
            {
                mpPtr->addInterface(0, 2, 3, 1);
                mpPtr->addInterface(2, 2, 5, 1);
                mpPtr->addInterface(4, 2, 7, 1);
                mpPtr->addInterface(6, 2, 1, 1);
            }

            break;
        }
        default:
            GISMO_ERROR("Unknown case makes gismo unhappy.");
    }
    gsMultiPatch<> & mp = *mpPtr;

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

    gsMultiBasis<> mb;
    if(doNURBS)
        mb.patchBases() = mp.basesCopy(false);//continue with NURBSBasis
    else
        mb.patchBases() = mp.basesCopy(true);//continue with B-Spline Basis
    mb.setTopology(mp.topology());

    if (nonMatching)
    {
        gsInfo << "Option NonMatching: Make uniform refinement for every third patch... " << std::flush;
        for (size_t i = 0; i < mb.nBases(); ++i)
            if ( i%3 == 0 )
                mb[i].uniformRefine();
        gsInfo << "done.\n";
        --refinements;

        gsInfo << "Option NonMatching: Increase spline degree for every other third patch... " << std::flush;
        for (size_t i = 0; i < mb.nBases(); ++i)
            if ( i%3 == 1 )
                mb[i].setDegreePreservingMultiplicity(degree+1);
        gsInfo << "done.\n";

        if (!dg)
        {
            gsInfo << "\nThe option --NonMatching does not allow a conforming discretization. Thus, option --DG is required.\n";
            return EXIT_FAILURE;
        }

    }

    gsInfo << "Setup bases and adjust degree... " << std::flush;

    if(!doNURBS)
        for ( size_t i = 0; i < mb.nBases(); ++ i )
            mb[i].setDegreePreservingMultiplicity(degree);
    else
        gsWarn << "\nsetDegreePreservingMultiplicity(...) is not implemented yet for NURBS!\n";

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
        (dirichlet::strategy) opt.getInt("MG.DirichletStrategy"),
        (iFace::strategy)     opt.getInt("MG.InterfaceStrategy")
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
        if ( smoother == "Richardson" || smoother == "r" )
            smootherOp = makeRichardsonOp(mg->matrix(i));
        else if ( smoother == "Jacobi" || smoother == "j" )
            smootherOp = makeJacobiOp(mg->matrix(i));
        else if ( smoother == "GaussSeidel" || smoother == "gs" )
            smootherOp = makeGaussSeidelOp(mg->matrix(i));
        else if ( smoother == "SubspaceCorrectedMassSmoother" || smoother == "scms" || smoother == "Hybrid" || smoother == "hyb"
            || smoother == "SubspaceCorrectedMassSmootherGeo" || smoother == "scmsg" || smoother == "HybridGeo" || smoother == "hybg"
            )
        {
            bool useGeo = smoother == "SubspaceCorrectedMassSmootherGeo" || smoother == "scmsg" || smoother == "HybridGeo" || smoother == "hybg";

            if (multiBases[i].nBases() == 1)
            {
                smootherOp = gsPreconditionerFromOp<>::make(
                    mg->underlyingOp(i),
                    gsPatchPreconditionersCreator<>::subspaceCorrectedMassSmootherOp(multiBases[i][0],bc,opt.getGroup("MG"),scaling,0,useGeo?&mp[0]:NULL)
                );

                if (autoDamping)
                {
                    opt.setReal( "MG.Damping", opt.getReal( "MG.Damping" ) / smootherOp->estimateLargestEigenvalue() );
                }

            }
            else
            {
                smootherOp = setupSubspaceCorrectedMassSmoother( mg->matrix(i), multiBases[i], mp, bc, opt.getGroup("MG"), useGeo, autoDamping);
            }

            if ( smoother == "Hybrid" || smoother == "hyb" || smoother == "HybridGeo" || smoother == "hybg" )
            {

                smootherOp->setOptions( opt.getGroup("MG") );
                smootherOp = gsCompositePrecOp<>::make( makeGaussSeidelOp(mg->matrix(i)), smootherOp );
            }
        }
        else
        {
            gsInfo << "\n\nThe chosen smoother is unknown.\n\nKnown are:\n  Richardson (r)\n  Jacobi (j)\n  GaussSeidel (gs)"
                      "\n  SubspaceCorrectedMassSmoother (scms)\n  Hybrid (hyb)\n  SubspaceCorrectedMassSmootherGeo (scmsg)\n  HybridGeo (hybg)\n\n";
            return EXIT_FAILURE;
        }

        smootherOp->setOptions( opt.getGroup("MG") );

        // Handle the extra-smooth option. On the finest grid level, there is nothing to handle.
        if (extrasmooth && i < mg->numLevels()-1 )
        {
            smootherOp->setNumOfSweeps( 1 << (mg->numLevels()-1-i) );
            smootherOp = gsPreconditionerFromOp<>::make(mg->underlyingOp(i),smootherOp);
        }
        mg->setSmoother(i, smootherOp);

    }

    gsMatrix<> x, errorHistory;
    x.setRandom( assembler.matrix().rows(), 1 );

    if (iterativeSolver=="cg")
        gsConjugateGradient<>( assembler.matrix(), mg )
            .setOptions( opt.getGroup("Solver") )
            .solveDetailed( assembler.rhs(), x, errorHistory );
    else if (iterativeSolver=="d")
        gsGradientMethod<>( assembler.matrix(), mg )
            .setOptions( opt.getGroup("Solver") )
            .solveDetailed( assembler.rhs(), x, errorHistory );
    else if (iterativeSolver=="gm")
        gsGMRes<>( assembler.matrix(), mg )
            .setOptions( opt.getGroup("Solver") )
            .solveDetailed( assembler.rhs(), x, errorHistory );
    else
    {
        gsInfo << "\n\nThe chosen iterative solver is unknown.\n\nKnown are:\n  conjugate gradient (cg)\n  direct (d)\n\n";
        return EXIT_FAILURE;
    }

    gsInfo << "done.\n\n";

    /******************** Print end Exit ********************/

    const index_t iter = errorHistory.rows()-1;

    gsInfo << "Number of dofs: " << assembler.matrix().rows() << "\n";
    const bool success = errorHistory(iter,0) < tolerance;
    if (success)
        gsInfo << "Reached desired tolerance after " << iter << " iterations:\n";
    else
        gsInfo << "Did not reach desired tolerance after " << iter << " iterations:\n";

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
        // Run paraview
        system("paraview multiGrid_result.pvd  &");    }
    else
    {
        gsInfo << "Done. No output created, re-run with --Plot to get a ParaView "
                  "file containing the solution.\n";
    }
    return success ? EXIT_SUCCESS : EXIT_FAILURE;
}


gsPreconditionerOp<>::Ptr setupSubspaceCorrectedMassSmoother(
    const gsSparseMatrix<>& matrix,
    const gsMultiBasis<>& mb,
    const gsMultiPatch<>& mp,
    const gsBoundaryConditions<>& bc,
    const gsOptionList& opt,
    bool useGeo,
    bool autoDamping
)
{
    const short_t dim = mb.topology().dim();

    // Setup dof mapper
    gsDofMapper dm;
    mb.getMapper(
       (dirichlet::strategy)opt.askInt("DirichletStrategy",11),
       (iFace    ::strategy)opt.askInt("InterfaceStrategy", 1),
       bc,
       dm,
       0
    );
    const index_t nTotalDofs = dm.freeSize();

    // Decompose the whole domain into components
    std::vector< std::vector<patchComponent> > components = mb.topology().allNonMatchingComponents(false);
    //std::vector< std::vector<patchComponent> > components = mb.topology().allComponents(false);
    const index_t nr_components = components.size();
    

    gsInfo << "Got "<<nr_components<<" components: \n";
    for (index_t i=0; i<nr_components; ++i)
    {
        for (size_t j=0; j<components[i].size(); ++j)
        {
            gsInfo << "[" << components[i][j].patch() << ":" << components[i][j].index() << "] ";
            //gsInfo << "corner  " << components[i][j].containedCorners().size() << "\n";
        }
        gsInfo << "\n";
    }


    // Setup Dirichlet boundary conditions
    gsBoundaryConditions<> dir_bc;
    for( index_t ps=0; ps < 2*dim; ++ps )
        dir_bc.addCondition( 0, 1+ps, condition_type::dirichlet, NULL );

    // Setup transfer matrices and local preconditioners
    std::vector< gsSparseMatrix<real_t,RowMajor> > transfers;
    transfers.reserve(nr_components);

    std::vector< gsLinearOperator<>::Ptr > ops;
    ops.reserve(nr_components);

    for (index_t i=0; i<nr_components; ++i)
    {
        gsMatrix<unsigned> indices;
        std::vector<gsBasis<>::uPtr> bases = mb.componentBasis_withIndices(components[i],dm,indices,true);

        index_t sz = indices.rows();
        gsSparseEntries<> se;
        se.reserve(sz);
        for (index_t i=0; i<sz; ++i)
            se.add(indices(i,0),i,real_t(1));
        gsSparseMatrix<real_t,RowMajor> transfer(nTotalDofs,sz);
        transfer.setFrom(se);

        if (sz>0)
        {
            if ((bases[0]->dim() == dim))
            {
                GISMO_ASSERT ( components[i].size() == 1, "Only one basis is expected for each patch." );
                ops.push_back(
                    gsPatchPreconditionersCreator<>::subspaceCorrectedMassSmootherOp(
                        *(bases[0]),
                        dir_bc,
                        gsOptionList(),
                        opt.getReal("Scaling"),
                        0,
                        useGeo ? &mp[components[i][0].patch()] : NULL
                    )
                );
            }
            else
            {
                gsSparseMatrix<> mat = transfer.transpose() * matrix * transfer;
                ops.push_back( makeSparseCholeskySolver(mat) );
            }

            transfers.push_back(give(transfer));
        }
    }

    gsAdditiveOp<>::uPtr prec = gsAdditiveOp<>::make(transfers, ops);
    if (autoDamping) prec->setRelativeScaling(matrix);
    return gsPreconditionerFromOp<>::make( makeMatrixOp(matrix), give(prec) );
}
