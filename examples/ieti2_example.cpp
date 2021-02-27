/** @file ieti2_example.cpp

    @brief Provides an example for the ieti solver

    Here, MINRES solves the saddle point formulation. For solving
    the Schur complement formulation with CG, see ieti_example.cpp.

    This class uses the gsPoissonAssembler, for the expression
    assembler, see ieti_example.cpp.

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
    real_t tolerance = 1.e-8;
    index_t maxIterations = 100;
    std::string out;
    bool plot = false;

    gsCmdLine cmd("Solves a PDE with an isogeometric discretization using an isogeometric tearing and interconnecting (IETI) solver.");
    cmd.addString("g", "Geometry",              "Geometry file", geometry);
    cmd.addInt   ("",  "SplitPatches",          "Split every patch that many times in 2^d patches", splitPatches);
    cmd.addReal  ("",  "StretchGeometry",       "Stretch geometry in x-direction by the given factor", stretchGeometry);
    cmd.addInt   ("r", "Refinements",           "Number of uniform h-refinement steps to perform before solving", refinements);
    cmd.addInt   ("p", "Degree",                "Degree of the B-spline discretization space", degree);
    cmd.addString("b", "BoundaryConditions",    "Boundary conditions", boundaryConditions);
    cmd.addReal  ("t", "Solver.Tolerance",      "Stopping criterion for linear solver", tolerance);
    cmd.addInt   ("",  "Solver.MaxIterations",  "Stopping criterion for linear solver", maxIterations);
    cmd.addString("",  "out",                   "Write solution and used options to file", out);
    cmd.addSwitch(     "plot",                  "Plot the result with Paraview", plot);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    gsOptionList opt = cmd.getOptionList();

    if ( ! gsFileManager::fileExists(geometry) )
    {
        gsInfo << "Geometry file could not be found.\n";
        gsInfo << "I was searching in the current directory and in: " << gsFileManager::getSearchPaths() << "\n";
        return EXIT_FAILURE;
    }

    gsInfo << "Run ieti_example with options:\n" << opt << std::endl;

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
    gsFunctionExpr<> f( "2*sin(x)*cos(y)", mp.geoDim() );

    // Dirichlet function
    gsFunctionExpr<> gD( "sin(x)*cos(y)", mp.geoDim() );

    // Neumann
    gsConstantFunction<> gN( 1.0, mp.geoDim() );

    gsBoundaryConditions<> bc;
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

    gsMultiBasis<> mb(mp);

    gsInfo << "Setup bases and adjust degree... " << std::flush;

    for ( size_t i = 0; i < mb.nBases(); ++ i )
        mb[i].setDegreePreservingMultiplicity(degree);

    for ( index_t i = 0; i < refinements; ++i )
        mb.uniformRefine();

    gsInfo << "done.\n";

    /********* Setup assembler and assemble matrix **********/

    gsInfo << "Setup assembler and assemble matrix... " << std::flush;

    const index_t nPatches = mp.nPatches();

    gsIetiMapper<> ietiMapper;

    // We start by setting up a global FeSpace that allows us to
    // obtain a dof mapper and the Dirichlet data
    {
      gsPoissonAssembler<> assembler(
            mp,
            mb,
            bc,
            f,
            dirichlet::elimination,
            iFace::glue
        );
        assembler.computeDirichletDofs();
        ietiMapper.init( mb, assembler.system().rowMapper(0), assembler.fixedDofs() );
    }

    // Compute the jump matrices
    bool fullyRedundant = true,
         noLagrangeMultipliersForCorners = true;
    ietiMapper.computeJumpMatrices(fullyRedundant, noLagrangeMultipliersForCorners);

    // We tell the ieti mapper which primal constraints we want; calling
    // more than one such function is possible.
    ietiMapper.cornersAsPrimals();

    // The ieti system does not have a special treatment for the
    // primal dofs. They are just one more subdomain
    gsIetiSystem<> ieti;
    ieti.reserve(nPatches+1);

    // The scaled Dirichlet preconditioner is independent of the
    // primal dofs.
    gsScaledDirichletPrec<> prec;
    prec.reserve(nPatches);

    // Setup the primal system, which needs to know the number of primal dofs.
    gsPrimalSystem<> primal(ietiMapper.nPrimalDofs());

    // Setup of the block-diagonal preconditioner for the saddle point problem
    // First, we need to know its size
    const index_t bdPrecSz = nPatches + 1 + (ietiMapper.nPrimalDofs()>0?1:0);
    gsBlockOp<>::Ptr bdPrec = gsBlockOp<>::make(bdPrecSz,bdPrecSz);

    for (index_t k=0; k<nPatches; ++k)
    {
        // We use the local variants of everything
        gsBoundaryConditions<> bc_local;
        bc.getConditionsForPatch(k,bc_local);
        gsMultiPatch<> mp_local = mp[k];
        gsMultiBasis<> mb_local = mb[k];

        // Setup assembler
        gsPoissonAssembler<> assembler(
            mp_local,
            mb_local,
            bc_local,
            f,
            dirichlet::elimination,
            iFace::glue
        );

        // This provides a new dof mapper and the Dirichlet data
        // This is necessary since it might happen that a 2d-patch touches the
        // Dirichlet boundary just with a corner or that a 3d-patch touches the
        // Dirichlet boundary with a corner or an edge. These cases are not
        // covered by bc.getConditionsForPatch
        gsDofMapper tmp = ietiMapper.dofMapperLocal(k);
        assembler.system() = gsSparseSystem<>(tmp); // gsSparseSystem takes per reference and steals contets
        assembler.setFixedDofVector(ietiMapper.fixedPart(k));

        // Assemble
        assembler.assemble();

        // Fetch data
        gsSparseMatrix<real_t, RowMajor> jumpMatrix  = ietiMapper.jumpMatrix(k);
        gsSparseMatrix<>                 localMatrix = assembler.matrix();
        gsMatrix<>                       localRhs    = assembler.rhs();

        // Add the patch to the scaled Dirichlet preconditioner
        //
        // This can be done using gsScaledDirichletPrec<>::restrictToSkeleton
        // as in ieti_example. Here, we call the underlying commands directly
        // to show how one can choose an alternative solver.
        std::vector<index_t> skeletonDofs = ietiMapper.skeletonDofs(k);

        gsScaledDirichletPrec<>::Blocks blocks
            = gsScaledDirichletPrec<>::matrixBlocks(localMatrix, skeletonDofs);

        prec.addSubdomain(
            prec.restrictJumpMatrix(jumpMatrix, skeletonDofs).moveToPtr(),
            gsScaledDirichletPrec<>::schurComplement( blocks, makeSparseCholeskySolver(blocks.A11) )
        );

        // Now, we handle the primal constraints.
        //
        // This can be done using primal.handleConstraints as in ieti_example.
        // Here, we call the underlying commands directly to show how one can
        // choose an alternative solver.
        gsSparseMatrix<>  modifiedLocalMatrix, localEmbedding, embeddingForBasis;
        gsMatrix<>        rhsForBasis;

        const bool eliminatePointwiseDofs = true;

        gsPrimalSystem<>::incorporateConstraints(
            ietiMapper.primalConstraints(k),
            eliminatePointwiseDofs,
            localMatrix,
            modifiedLocalMatrix,
            localEmbedding,
            embeddingForBasis,
            rhsForBasis
        );

        gsLinearOperator<>::Ptr localSolver = makeSparseLUSolver(modifiedLocalMatrix);

        primal.addContribution(
            jumpMatrix, localMatrix, localRhs,
            gsPrimalSystem<>::primalBasis(
                localSolver, embeddingForBasis, rhsForBasis, ietiMapper.primalDofIndices(k), primal.nPrimalDofs()
            )
        );
        gsMatrix<>                       modifiedLocalRhs     = localEmbedding.transpose() * localRhs;
        gsSparseMatrix<real_t, RowMajor> modifiedJumpMatrix   = jumpMatrix * localEmbedding;


        // Register the local solver to the block preconditioner. We use
        // a sparse LU solver since the local saddle point problem is not
        // positive definite.
        bdPrec->addOperator(k,k,localSolver);

        // Add the patch to the Ieti system
        ieti.addSubdomain(
            modifiedJumpMatrix.moveToPtr(),
            makeMatrixOp(modifiedLocalMatrix.moveToPtr()),
            give(modifiedLocalRhs)
        );
    }

    // Add the primal problem if there are primal constraints
    if (ietiMapper.nPrimalDofs()>0)
    {
        // Register the local solver to the block preconditioner
        bdPrec->addOperator(nPatches, nPatches, makeSparseCholeskySolver(primal.localMatrix()));

        // Add to IETI system
        ieti.addSubdomain(
            primal.jumpMatrix().moveToPtr(),
            makeMatrixOp(primal.localMatrix().moveToPtr()),
            give(primal.localRhs())
        );
    }

    gsInfo << "done.\n";

    /**************** Setup solver and solve ****************/

    gsInfo << "Setup solver and solve... \n"
        "    Setup multiplicity scaling... " << std::flush;

    // Tell the preconditioner to set up the scaling
    prec.setupMultiplicityScaling();

    // The scaled Dirichlet preconditioner is in the last block
    gsLinearOperator<>::Ptr sdPrec = prec.preconditioner();
    bdPrec->addOperator(bdPrecSz-1,bdPrecSz-1,sdPrec);

    gsInfo << "done.\n    Setup minres solver and solve... " << std::flush;
    // Initial guess
    gsMatrix<> x;
    x.setRandom( bdPrec->rows(), 1 );

    // This is the main cg iteration
    gsMatrix<> errorHistory;
    gsMinimalResidual<>( ieti.saddlePointProblem(), bdPrec )
        .setOptions( opt.getGroup("Solver") )
        .solveDetailed( ieti.rhsForSaddlePoint(), x, errorHistory );

    gsInfo << "done.\n    Reconstruct solution from Lagrange multipliers... " << std::flush;

    // Now, we want to have the global solution for u
    gsMatrix<> uVec = ietiMapper.constructGlobalSolutionFromLocalSolutions(
        primal.distributePrimalSolution(
            ieti.constructSolutionFromSaddlePoint(x)
        )
    );
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
        fd.add(opt);
        fd.add(uVec);
        fd.addComment(std::string("ieti2_example   Timestamp:")+std::ctime(&time));
        fd.save(out);
        gsInfo << "Write solution to file " << out << "\n";
    }

    if (plot)
    {
        gsInfo << "Write Paraview data to file multiGrid_result.pvd\n";
        gsPoissonAssembler<> assembler(
            mp,
            mb,
            bc,
            f,
            dirichlet::elimination,
            iFace::glue
        );
        assembler.computeDirichletDofs();
        gsMultiPatch<> mpsol;
        assembler.constructSolution(uVec, mpsol);
        gsField<> sol( assembler.patches(), mpsol );
        gsWriteParaview<>(sol, "ieti_result", 1000);
        //gsFileManager::open("ieti_result.pvd");
    }
    if (!plot&&out.empty())
    {
        gsInfo << "Done. No output created, re-run with --plot to get a ParaView "
                  "file containing the solution or --out to write solution to xml file.\n";
    }
    return success ? EXIT_SUCCESS : EXIT_FAILURE;
}
