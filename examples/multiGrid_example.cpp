/** @file multiGrid_example.cpp

    @brief Provides an examples for the multigrid solver.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): C. Hofreither, S. Takacs
*/

#include <ctime>

#include <gismo.h>

using namespace gismo;

gsPreconditionerOp<>::Ptr setupSubspaceCorrectedMassSmoother(index_t, index_t, const gsSparseMatrix<>&,
    const gsMultiBasis<>&, const gsBoundaryConditions<>&, const gsOptionList&, std::vector<real_t>& );

int main(int argc, char *argv[])
{
    /************** Define command line options *************/

    std::string geometry("domain2d/yeti_mp2.xml");
    index_t splitPatches = 0;
    real_t stretchGeometry = 1;
    index_t xRefine = 0;
    index_t refinements = 3;
    index_t degree = 2;
    bool nonMatching = false;
    bool dg = false;
    bool nitsche = false;
    index_t levels = -1;
    index_t cycles = 1;
    index_t presmooth = 1;
    index_t postsmooth = 1;
    bool extrasmooth = false;
    std::string smoother("GaussSeidel");
    real_t damping = -1;
    real_t scaling = 0.12;
    std::string iterativeSolver("cg");
    real_t tolerance = 1.e-8;
    index_t maxIterations = 100;
    std::string boundaryConditions("d");
    std::string out;
    bool plot = false;

    gsCmdLine cmd("Solves a PDE with an isogeometric discretization using a multigrid solver.");
    cmd.addString("g", "Geometry",              "Geometry file", geometry);
    cmd.addInt   ("",  "SplitPatches",          "Split every patch that many times in 2^d patches", splitPatches);
    cmd.addReal  ("",  "StretchGeometry",       "Stretch geometry in x-direction by the given factor", stretchGeometry);
    cmd.addInt   ("",  "XRefine",               "Refine in x-direction by adding that many knots to every knot span", xRefine);
    cmd.addInt   ("r", "Refinements",           "Number of uniform h-refinement steps to perform before solving", refinements);
    cmd.addInt   ("p", "Degree",                "Degree of the B-spline discretization space", degree);
    cmd.addSwitch("",  "NonMatching",           "Set up a non-matching multi-patch discretization", nonMatching);
    cmd.addSwitch("",  "DG",                    "Use discontinuous Galerkin (SIPG) for coupling the patches", dg);
    cmd.addSwitch("",  "Nitsche",               "Use Nitsche method for Dirichlet boundary conditions", nitsche);
    cmd.addInt   ("l", "MG.Levels",             "Number of levels to use for multigrid iteration", levels);
    cmd.addInt   ("c", "MG.NumCycles",          "Number of multi-grid cycles", cycles);
    cmd.addInt   ("",  "MG.NumPreSmooth",       "Number of pre-smoothing steps", presmooth);
    cmd.addInt   ("",  "MG.NumPostSmooth",      "Number of post-smoothing steps", postsmooth);
    cmd.addSwitch("",  "MG.Extrasmooth",        "Doubles the number of smoothing steps for each coarser level", extrasmooth);
    cmd.addString("s", "MG.Smoother",           "Smoothing method", smoother);
    cmd.addReal  ("",  "MG.Damping",            "Damping factor for the smoother", damping);
    cmd.addReal  ("",  "MG.Scaling",            "Scaling factor for the subspace corrected mass smoother", scaling);
    cmd.addString("i", "IterativeSolver",       "Iterative solver: apply multigrid directly (d) or as a preconditioner for "
                                                "conjugate gradient (cg)", iterativeSolver);
    cmd.addReal  ("t", "Solver.Tolerance",      "Stopping criterion for linear solver", tolerance);
    cmd.addInt   ("",  "Solver.MaxIterations",  "Stopping criterion for linear solver", maxIterations);
    cmd.addString("b", "BoundaryConditions",    "Boundary conditions", boundaryConditions);
    cmd.addString("",  "out",                   "Write solution and used options to file", out);
    cmd.addSwitch(     "plot",                  "Plot the result with Paraview", plot);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    // Default case is levels:=refinements, so replace invalid default accordingly
    if (levels <0) { levels = refinements; cmd.setInt( "MG.Levels", levels ); }
    // The smoothers know their defaults, so remove the invalid default
    if (damping<0) { cmd.remove( "MG.Damping" ); }

    // Define assembler options
    cmd.remove( "DG" );
    cmd.addInt( "MG.InterfaceStrategy", "", (index_t)( dg      ? iFace::dg          : iFace::conforming      ) );
    cmd.addInt( "MG.DirichletStrategy", "", (index_t)( nitsche ? dirichlet::nitsche : dirichlet::elimination ) );

    if ( ! gsFileManager::fileExists(geometry) )
    {
        gsInfo << "Geometry file could not be found.\n";
        gsInfo << "I was searching in the current directory and in: " << gsFileManager::getSearchPaths() << "\n";
        return EXIT_FAILURE;
    }

    gsInfo << "Run multiGrid_example with options:\n" << cmd << "\n";

    /******************* Define geometry ********************/

    gsInfo << "Define geometry... " << std::flush;

    //! [Define Geometry]
    gsMultiPatch<>::uPtr mpPtr = gsReadFile<>(geometry);
    //! [Define Geometry]
    if (!mpPtr)
    {
        gsInfo << "No geometry found in file " << geometry << ".\n";
        return EXIT_FAILURE;
    }
    gsMultiPatch<>& mp = *mpPtr;

    //! [Define Geometry2]
    for (index_t i=0; i<splitPatches; ++i)
    {
        gsInfo << "split patches uniformly... " << std::flush;
        mp = mp.uniformSplit();
    }
    //! [Define Geometry2]

    // This allows to strech a single-patch geometry in x-direction.
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

    gsInfo << "Define boundary conditions... " << std::flush;

    //! [Define Source]
    // Right-hand-side
    gsFunctionExpr<> f( "2*sin(x)*cos(y)", mp.geoDim() );

    // Dirichlet function
    gsFunctionExpr<> gD( "sin(x)*cos(y)", mp.geoDim() );

    // Neumann
    gsConstantFunction<> gN( 1.0, mp.geoDim() );

    gsBoundaryConditions<> bc;
    //! [Define Source]
    {
        const index_t len = boundaryConditions.length();
        index_t i = 0;
        for (gsMultiPatch<>::const_biterator it = mp.bBegin(); it < mp.bEnd(); ++it)
        {
            char b_local;
            if ( len == 1 )
                b_local = boundaryConditions[0];
            else if ( i < len )
                b_local = boundaryConditions[i];
            else
            {
                gsInfo << "\nNot enough boundary conditions given.\n";
                return EXIT_FAILURE;
            }

            if ( b_local == 'd' )
                bc.addCondition( *it, condition_type::dirichlet, &gD );
            else if ( b_local == 'n' )
                bc.addCondition( *it, condition_type::neumann, &gN );
            else
            {
                gsInfo << "\nInvalid boundary condition given; only 'd' (Dirichlet) and 'n' (Neumann) are supported.\n";
                return EXIT_FAILURE;
            }

            ++i;
        }
        if ( len > i )
            gsInfo << "\nToo many boundary conditions have been specified. Ignoring the remaining ones.\n";
        gsInfo << "done. "<<i<<" boundary conditions set.\n";
    }


    /************ Setup bases and adjust degree *************/

    //! [Define Basis]
    gsMultiBasis<> mb(mp);
    //! [Define Basis]

    if (xRefine)
    {
        if (mp.nPatches() > 1)
        {
            gsInfo << "\nThe option --XRefine is only available for single-patch geometries..\n";
            return EXIT_FAILURE;
        }
        gsInfo << "Option --XRefine: Refine the grid in x-direction by interting " << xRefine
               << " knots in every knot span... " << std::flush;
        mb[0].component(0).uniformRefine(xRefine,1);
        gsInfo << "done.\n";
    }

    gsInfo << "Setup bases and adjust degree... " << std::flush;

    //! [Set degree and refine]
    for ( size_t i = 0; i < mb.nBases(); ++ i )
        mb[i].setDegreePreservingMultiplicity(degree);

    for ( index_t i = 0; i < refinements; ++i )
        mb.uniformRefine();
    //! [Set degree and refine]

    gsInfo << "done.\n";

    //! [Define Non Matching]
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
    //! [Define Non Matching]

    /********* Setup assembler and assemble matrix **********/

    gsInfo << "Setup assembler and assemble matrix... " << std::flush;

    //! [Assemble]
    gsPoissonAssembler<> assembler(
        mp,
        mb,
        bc,
        f,
        (dirichlet::strategy) cmd.getInt("MG.DirichletStrategy"),
        (iFace::strategy)     cmd.getInt("MG.InterfaceStrategy")
    );
    assembler.assemble();
    //! [Assemble]

    gsInfo << "done.\n";

    /**************** Setup solver and solve ****************/

    gsInfo << "Setup solver and solve... " << std::flush;

    //! [Define vectors]
    std::vector< gsSparseMatrix<real_t,RowMajor> > transferMatrices;
    std::vector< gsMultiBasis<real_t> > multiBases;  // Needed for setupSubspaceCorrectedMassSmoother
    std::vector<real_t> patchLocalDampingParameters; // Needed for setupSubspaceCorrectedMassSmoother
    //! [Define vectors]

    // Setup grid hiearachy by coarsening of the given matrix
    // We move the constructed hiearchy of multi bases into a variable (only required for the subspace smoother)
    // Then we move the transfer matrices into a variable
    //! [Setup grid hierarchy]
    gsGridHierarchy<>::buildByCoarsening(give(mb), bc, cmd.getGroup("MG"))
        .moveMultiBasesTo(multiBases)
        .moveTransferMatricesTo(transferMatrices);
    //! [Setup grid hierarchy]

    // Setup the multigrid solver
    //! [Setup multigrid]
    gsMultiGridOp<>::Ptr mg = gsMultiGridOp<>::make( assembler.matrix(), transferMatrices );
    mg->setOptions( cmd.getGroup("MG") );
    //! [Setup multigrid]

    // Since we are solving a symmetric positive definite problem,we can use a Cholesky solver
    // (instead of the LU solver that would be created by default).
    //
    // mg->matrix(0) gives the matrix for the coarsest grid level (=level 0).
    //! [Define coarse solver]
    mg->setCoarseSolver( makeSparseCholeskySolver( mg->matrix(0) ) );
    //! [Define coarse solver]

    // Set up of the smoothers
    // This has to be done for each grid level separately
    //! [Define smoothers]
    for (index_t i = 1; i < mg->numLevels(); ++i)
    {
        gsPreconditionerOp<>::Ptr smootherOp;
        if ( smoother == "Richardson" || smoother == "r" )
            smootherOp = makeRichardsonOp(mg->matrix(i));
        else if ( smoother == "Jacobi" || smoother == "j" )
            smootherOp = makeJacobiOp(mg->matrix(i));
        else if ( smoother == "GaussSeidel" || smoother == "gs" )
            smootherOp = makeGaussSeidelOp(mg->matrix(i));
        else if ( smoother == "IncompleteLU" || smoother == "ilu" )
            smootherOp = makeIncompleteLUOp(mg->matrix(i));
        else if ( smoother == "SubspaceCorrectedMassSmoother" || smoother == "scms" )
            smootherOp = setupSubspaceCorrectedMassSmoother( i, mg->numLevels(), mg->matrix(i),
                multiBases[i], bc, cmd.getGroup("MG"), patchLocalDampingParameters );
        else if ( smoother == "Hybrid" || smoother == "hyb" )
            smootherOp = gsCompositePrecOp<>::make(
                makeGaussSeidelOp(mg->matrix(i)),
                setupSubspaceCorrectedMassSmoother( i, mg->numLevels(), mg->matrix(i),
                    multiBases[i], bc, cmd.getGroup("MG"), patchLocalDampingParameters )
                );
        //! [Define smoothers]
        else
        {
            gsInfo << "\n\nThe chosen smoother is unknown.\n\nKnown are:\n  Richardson (r)\n  Jacobi (j)\n  GaussSeidel (gs)"
                      "\n  IncompleteLU (ilu)\n  SubspaceCorrectedMassSmoother (scms)\n  Hybrid (hyb)\n\n";
            return EXIT_FAILURE;
        }

        //! [Define smoothers2]
        smootherOp->setOptions( cmd.getGroup("MG") );
        //! [Define smoothers2]

        // Handle the extra-smooth option. On the finest grid level, there is nothing to handle.
        if (extrasmooth && i < mg->numLevels()-1)
        {
            smootherOp->setNumOfSweeps( 1 << (mg->numLevels()-1-i) );
            smootherOp = gsPreconditionerFromOp<>::make(mg->underlyingOp(i),smootherOp);
        }

    //! [Define smoothers3]
        mg->setSmoother(i, smootherOp);
    } // end for
    //! [Define smoothers3]

    gsMatrix<> errorHistory;

    //! [Initial guess]
    gsMatrix<> x;
    x.setRandom( assembler.matrix().rows(), 1 );
    //! [Initial guess]

    //! [Solve]
    if (iterativeSolver=="cg")
        gsConjugateGradient<>( assembler.matrix(), mg )
            .setOptions( cmd.getGroup("Solver") )
            .solveDetailed( assembler.rhs(), x, errorHistory );
    else if (iterativeSolver=="d")
        gsGradientMethod<>( assembler.matrix(), mg )
            .setOptions( cmd.getGroup("Solver") )
            .solveDetailed( assembler.rhs(), x, errorHistory );
    //! [Solve]
    else
    {
        gsInfo << "\n\nThe chosen iterative solver is unknown.\n\nKnown are:\n  conjugate gradient (cg)\n  direct (d)\n\n";
        return EXIT_FAILURE;
    }

    gsInfo << "done.\n\n";

    /******************** Print end Exit ********************/

    const index_t iter = errorHistory.rows()-1;
    const bool success = errorHistory(iter,0) < tolerance;
    if (success)
        gsInfo << "Reached desired tolerance after " << iter << " iterations:\n";
    else
        gsInfo << "Did not reach desired tolerance after " << iter << " iterations:\n";

    if (errorHistory.rows() < 20)
        gsInfo << errorHistory.transpose() << "\n\n";
    else
        gsInfo << errorHistory.topRows(5).transpose() << " ... " << errorHistory.bottomRows(5).transpose()  << "\n\n";

    if (!out.empty())
    {
        gsFileData<> fd;
        std::time_t time = std::time(NULL);
        fd.add(cmd);
        fd.add(x);
        fd.addComment(std::string("multiGrid_example   Timestamp:")+std::ctime(&time));
        fd.save(out);
        gsInfo << "Write solution to file " << out << "\n";
    }

    if (plot)
    {
        // Construct the solution as a scalar field
        gsMultiPatch<> mpsol;
        assembler.constructSolution(x, mpsol);
        gsField<> sol( assembler.patches(), mpsol );

        // Write solution to paraview files
        gsInfo << "Write Paraview data to file multiGrid_result.pvd\n";
        gsWriteParaview<>(sol, "multiGrid_result", 1000);
        gsFileManager::open("multiGrid_result.pvd");
    }
    if (!plot&&out.empty())
    {
        gsInfo << "Done. No output created, re-run with --plot to get a ParaView "
                  "file containing the solution or --out to write solution to xml file.\n";
    }
    return success ? EXIT_SUCCESS : EXIT_FAILURE;
}


gsPreconditionerOp<>::Ptr setupSubspaceCorrectedMassSmoother(
    index_t level,
    index_t nrLevels,
    const gsSparseMatrix<>& matrix,
    const gsMultiBasis<>& mb,
    const gsBoundaryConditions<>& bc,
    const gsOptionList& opt,
    std::vector<real_t>& patchLocalDampingParameters
)
{
    // This function sets up a patchwise version of the subspace corrected mass smoother.

    const short_t dim = mb.topology().dim();

    const iFace::strategy iFaceStrategy = (iFace::strategy)opt.askInt("InterfaceStrategy", 1);
    GISMO_ASSERT( iFaceStrategy == iFace::dg || iFaceStrategy == iFace::conforming,
      "Unknown interface strategy." );

    // Setup dof mapper
    gsDofMapper dm;
    mb.getMapper(
       (dirichlet::strategy)opt.askInt("DirichletStrategy",11),
       iFaceStrategy,
       bc,
       dm,
       0
    );
    const index_t nTotalDofs = dm.freeSize();

    // Decompose the whole domain into components
    std::vector< std::vector<patchComponent> > components = mb.topology().allComponents(true);
    const index_t nr_components = components.size();

    if (patchLocalDampingParameters.size() == 0)
        patchLocalDampingParameters.resize(nr_components);

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
        const index_t cdim = components[i][0].dim();
        gsMatrix<index_t> indices;
        std::vector<gsBasis<>::uPtr> bases = mb.componentBasis_withIndices(components[i],dm,indices,true);
        index_t sz = indices.rows();
        gsSparseEntries<> se;
        se.reserve(sz);
        for (index_t j=0; j<sz; ++j)
            se.add(indices(j,0),j,(real_t)1);
        gsSparseMatrix<real_t,RowMajor> transfer(nTotalDofs,sz);
        transfer.setFrom(se);

        if (sz>0)  // It might be the case that there is no active dof in some components.
        {
            gsSparseMatrix<> localMatrix = transfer.transpose() * matrix * transfer;

            if (cdim >= 2 && cdim >= dim-1)
            {
                // If the component has dimension of at least 2, we have to do something fancy.
                // For dim > 3, we ignore interaces of dimensionality dim-2 and smaller (for simplicity).
                index_t nrDifferentBases;
                real_t stiffFactor, massFactor;
                if ( cdim == dim )
                {
                    nrDifferentBases = 1;
                    stiffFactor = 1;
                    massFactor = 0;
                }
                else // cdim == dim-1
                {
                    if ( iFaceStrategy == iFace::dg && components[i].size() == 2 )
                    {
                        // components[i].size() == 2 means that we consider an interface, not a boundary
                        const index_t patch0 = components[i][0].patch(),
                                      patch1 = components[i][1].patch();
                        const real_t h = math::min( mb[patch0].getMinCellLength(), mb[patch1].getMinCellLength() );
                        const index_t p = math::max( mb[patch0].maxDegree(), mb[patch1].maxDegree() );
                        nrDifferentBases = 2;
                        stiffFactor = h/p;
                        massFactor = p/h;
                    }
                    else
                    {
                        const index_t patch = components[i][0].patch();
                        const real_t h  = mb[patch].getMinCellLength();
                        const index_t p = mb[patch].maxDegree();
                        nrDifferentBases = 1;
                        stiffFactor = h/p;
                        massFactor = p/h;
                    }
                }

                gsLinearOperator<>::uPtr localOperators[2];
                for (index_t j=0; j!=nrDifferentBases; ++j)
                {
                    localOperators[j] = gsPatchPreconditionersCreator<>::subspaceCorrectedMassSmootherOp(
                        *(bases[j]),
                        dir_bc,
                        gsOptionList(),
                        opt.getReal("Scaling"),
                        massFactor,
                        stiffFactor
                    );
                }

                if ( nrDifferentBases > 1 )
                {
                    // In the case of dG, we use a conjugate gradient solver with a block-diagonal
                    // preconditioner
                    gsLinearOperator<>::Ptr underlying = makeMatrixOp( localMatrix.moveToPtr() );

                    gsBlockOp<>::uPtr blockPreconder = gsBlockOp<>::make(nrDifferentBases, nrDifferentBases);
                    for (index_t j=0; j!=nrDifferentBases; ++j)
                        blockPreconder->addOperator(j,j,give(localOperators[j]));

                    ops.push_back( gsIterativeSolverOp< gsConjugateGradient<> >::make(underlying, give(blockPreconder)) );
                }
                else
                {
                    // Otherwiese, we just scale the local solvers properly
                    gsLinearOperator<>::Ptr underlying = makeMatrixOp( localMatrix.moveToPtr() );
                    gsPreconditionerFromOp<>::uPtr pc = gsPreconditionerFromOp<>::make(underlying,give(localOperators[0]));
                    if (patchLocalDampingParameters[i] == 0)
                    {
                        const real_t damping = 1/pc->estimateLargestEigenvalueOfPreconditionedSystem(10);
                        pc->setDamping(damping);

                        // We store the local damping parameter if we can expect that it does not change
                        // too much any more. When the number of inner knots is p or smaller, the subspace
                        // corrected mass smoother is an exact solver. Thus, these cases are not comparable.
                        bool saveLambda = (level > nrLevels - 5);
                        for (index_t j=0; j!=bases[0]->dim(); ++j)
                            saveLambda &= bases[0]->component(j).numElements()
                                              > static_cast<size_t>( bases[0]->component(j).maxDegree() );
                        if ( saveLambda )
                            patchLocalDampingParameters[i] = damping;
                    }
                    else
                        pc->setDamping(patchLocalDampingParameters[i]);
                    ops.push_back(give(pc));
                }

            }
            else
            {
                // If the component has dimension 0 or 1, we can just use direct solves
                ops.push_back( makeSparseCholeskySolver(localMatrix) );
            }

            transfers.push_back(give(transfer));
        }
    }

    gsPreconditionerFromOp<>::Ptr result
        = gsPreconditionerFromOp<>::make(
            makeMatrixOp(matrix),
            gsAdditiveOp<>::make(transfers, ops)
        );
    result->setOptions(opt);
    return result;
}
