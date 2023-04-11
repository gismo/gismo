/** @file ietidG2_example.cpp

    @brief Provides an example for the ieti solver for a dG setting

    This class uses the gsPoissonAssembler. We use CG to solve for the
    Schur complement formulation.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): S. Takacs, R. Schneckenleitner
*/
#include <ctime>

#ifdef NDEBUG
#define DEBUGVAR(v)
#define DEBUGMAT(m)
#else
#define DEBUGVAR(v) gsInfo << #v << ": " << (v) << " (line nr. " << __LINE__ << ")\n"
#define DEBUGMAT(m) gsInfo << #m << ": " << (m).rows() << "x" << (m).cols() << " (line nr. " << __LINE__ << ")\n"
#endif

#include <gismo.h>
#include <gsAssembler/gsVisitorDg.h>
#include <gsIeti/gsIetidGMapper.h>
#include <gsSolver/gsPreconditioner.h>
#include <unsupported/src/gsSolver/gsTimedOp.h>

using namespace gismo;

gsMultiPatch<real_t> approximateQuarterAnnulus(index_t deg)
{
    gsGeometry<>::uPtr quann = gsNurbsCreator<>::NurbsQuarterAnnulus();

    gsKnotVector<> KV1(0,1, 0, deg+1);        // no interior knots in x direction
    gsKnotVector<> KV2(0,1, 0, deg+1);        // no interior knot in y direction

    gsTensorBSplineBasis<2> tbsp (give(KV1), give(KV2));
    gsMatrix<real_t> eval = quann->eval(tbsp.anchors());
    gsGeometry<>::uPtr approxGeom = tbsp.interpolateAtAnchors( eval );
    gsMultiPatch<real_t> mp(*approxGeom);

    return mp;

}

/// struct to measure timings related to the setup of and solving with the IETI method
struct timings
{
    real_t assembling_time = 0;
    real_t setup_time = 0;
    real_t solving_time = 0;

    void print(std::ostream& out)
    {
        std::ios  state(NULL);
        state.copyfmt(std::cout);

        //out<<std::setprecision(4)<<"\n\n";
        out<<std::setw(12)<<"Assembling time: "<<assembling_time<<"\n";
        out<<std::setw(12)<<"Setup: "<<setup_time<<"\n";
        out<<std::setw(12)<<"Solving time: "<<solving_time<<"\n";

        out<<"\n";
        std::cout.copyfmt(state);
    }
};

int main(int argc, char *argv[])
{
    /************** Define command line options *************/

    std::string geometry("domain2d/square.xml");
    index_t splitPatches = 1;
    real_t stretchGeometry = 1;
    index_t refinements = 1;
    index_t degree = 2;
    bool nonMatching = false;
    real_t alpha = 1;
    real_t beta = 1;
    real_t penalty = 5;
    real_t sigma = 1;
    std::string boundaryConditions("d");
    std::string primals("c");
    bool eliminatePointwiseDofs = !true; // true for MFD
    real_t tolerance = 1.e-6;
    index_t maxIterations = 20000;
    std::string out;
    bool plot = false;

    gsCmdLine cmd("Solves a PDE with an isogeometric discretization using an isogeometric tearing and interconnecting (IETI) solver.");
    cmd.addString("g", "Geometry",              "Geometry file", geometry);
    cmd.addInt   ("",  "SplitPatches",          "Split every patch that many times in 2^d patches", splitPatches);
    cmd.addReal  ("",  "StretchGeometry",       "Stretch geometry in x-direction by the given factor", stretchGeometry);
    cmd.addInt   ("r", "Refinements",           "Number of uniform h-refinement steps to perform before solving", refinements);
    cmd.addInt   ("p", "Degree",                "Degree of the B-spline discretization space", degree);
    cmd.addSwitch("",  "NonMatching",           "Set up a non-matching multi-patch discretization", nonMatching);
    cmd.addReal  ("",  "DG.Alpha",              "Parameter alpha for dG scheme; use 1 for SIPG and NIPG.", alpha );
    cmd.addReal  ("",  "DG.Beta",               "Parameter beta for dG scheme; use 1 for SIPG and -1 for NIPG", beta );
    cmd.addReal  ("",  "DG.Penalty",            "Penalty parameter delta for dG scheme; if negative, default 4(p+d)(p+1) is used.", penalty );
    cmd.addReal  ("",  "Sigma",                 "Scaling parameter sigma for scaled Dirichlet preconder", sigma );
    cmd.addString("b", "BoundaryConditions",    "Boundary conditions", boundaryConditions);
    cmd.addString("c", "Primals",               "Primal constraints (c=corners, e=edges, f=faces)", primals);
    //cmd.addSwitch("e", "EliminateCorners",      "Eliminate corners (if they are primals)", eliminatePointwiseDofs);
    cmd.addReal  ("t", "Solver.Tolerance",      "Stopping criterion for linear solver", tolerance);
    cmd.addInt   ("",  "Solver.MaxIterations",  "Maximum number of iterations for linear solver", maxIterations);
    cmd.addString("",  "out",                   "Write solution and used options to file", out);
    cmd.addSwitch(     "plot",                  "Plot the result with Paraview", plot);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    //if ( ! gsFileManager::fileExists(geometry) )
    //{
    //    gsInfo << "Geometry file could not be found.\n";
    //    gsInfo << "I was searching in the current directory and in: " << gsFileManager::getSearchPaths() << "\n";
    //    return EXIT_FAILURE;
    //}

    gsInfo << "Run ietidG4_example with options:\n" << cmd << std::endl;

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
    //! [Define Geometry2]
    gsMultiPatch<>& mp = *mpPtr;
    {

        /*std::vector<gsGeometry<>*> ptch;
        for(size_t np = 0; np< mp.nPatches(); ++np)
        {
            std::vector<gsGeometry<>* >  ptch_ =  mp.patch(np).uniformSplit(1);
            ptch.insert(ptch.end(),ptch_.begin(),ptch_.end());
        }
        mp = gsMultiPatch<real_t>(ptch);*/

        for (index_t i=0; i<splitPatches; ++i)
        {
            gsInfo << "split patches uniformly... " << std::flush;
            mp = mp.uniformSplit();
        }
        mp.computeTopology();
    }
    //! [Define Geometry2]

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

    //! [Define Source]
    // Right-hand-side
    gsFunctionExpr<> f( "2*pi^2*sin(pi*x)*sin(pi*y)", mp.geoDim() );

    // Dirichlet function
    gsFunctionExpr<> gD( "2*pi^2*sin(pi*x)*sin(pi*y)", mp.geoDim() );

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

    gsInfo << "Setup bases and adjust degree... " << std::flush;

    //! [Define Basis]
    gsMultiBasis<> mb(mp);
    //! [Define Basis]

    //! [Set degree and refine]
    for ( size_t i = 0; i < mb.nBases(); ++ i )
        mb[i].setDegreePreservingMultiplicity(degree);

    for ( index_t i = 0; i < refinements; ++i )
        mb.uniformRefine();

    // We might want to create a non-matching setup such that we are not in the
    // special case where the function spaces on the interfaces actually agree
    // (as this would be the case for the Yeti footprint, the default domain.)
    //gsPiecewiseFunction<> piece(mb.size());
    //gsFunctionExpr<> one("1.0", 2);
    //gsFunctionExpr<> two("2.0", 2);
    //gsFunctionExpr<> zero("0.0", 2);
    if (nonMatching)
    {
        gsInfo << "\n  Option NonMatching: Make uniform refinement for every third patch... " << std::flush;
        for (size_t i = 0; i < mb.nBases(); ++i)
            if ( i%3 == 0 ) {
                mb[i].uniformRefine();
                //piece.addPiece(zero);
            }
        gsInfo << "done.\n";
        --refinements;

        gsInfo << "  Option NonMatching: Increase spline degree for every other third patch... " << std::flush;
        for (size_t i = 0; i < mb.nBases(); ++i)
            if ( i%3 == 1 ) {
                mb[i].setDegreePreservingMultiplicity(degree+1);
                //piece.addPiece(one);
            }
        gsInfo << "done.\n";

    }
    //! [Set degree and refine]

    gsInfo << "done.\n";

    /************** Compute penalty parameter **************/
    if(penalty < 0)
        penalty *= mb.maxCwiseDegree();

    /********* Setup assembler and assemble matrix **********/

    gsStopwatch timer;
    timings times;

    gsInfo << "Setup assembler and assemble matrix... " << std::flush;

    ///timings setup_time begin
    timer.restart();

    const index_t nPatches = mp.nPatches();
    //! [Define global mapper]
    gsIetidGMapper<> ietiMapper;
    {
        // We start by setting up a global assembler that allows us to
        // obtain a dof mapper and the Dirichlet data
        gsOptionList assemblerOptions = gsGenericAssembler<>::defaultOptions();
        assemblerOptions.setInt("DirichletStrategy", dirichlet::elimination);
        assemblerOptions.setInt("InterfaceStrategy", iFace::dg);
        gsGenericAssembler<> assembler(
            mp,
            mb,
            assemblerOptions,
            &bc
        );
        assembler.computeDirichletDofs();
        //! [Define global mapper]

        //! [Define Ieti Mapper]
        ietiMapper.init(
            mb,
            assembler.system().rowMapper(0),
            assembler.fixedDofs(),
            gsIetidGMapper<>::allArtificialIfaces(mb)
        );
    }
    //! [Define Ieti Mapper]

    // Which primal dofs should we choose?
    bool cornersAsPrimals = false, edgesAsPrimals = false, facesAsPrimals = false;
    for (size_t i=0; i<primals.length(); ++i)
        switch (primals[i])
        {
            case 'c': cornersAsPrimals = true;   break;
            case 'e': edgesAsPrimals = true;     break;
            case 'f': facesAsPrimals = true;     break;
            default:
                gsInfo << "\nUnkown type of primal constraint: \"" << primals[i] << "\"\n";
                return EXIT_FAILURE;
        }

    // Compute the jump matrices
    bool fullyRedundant = true,
         noLagrangeMultipliersForCorners = cornersAsPrimals;
    //! [Define jumps]
    ietiMapper.computeJumpMatrices(fullyRedundant, noLagrangeMultipliersForCorners);
    //! [Define jumps]

    // We tell the ieti mapper which primal constraints we want; calling
    // more than one such function is possible.
    //! [Define primals]
    if (cornersAsPrimals)
        ietiMapper.cornersAsPrimals();

    if (edgesAsPrimals)
        ietiMapper.interfaceAveragesAsPrimals(mp,1);

    if (facesAsPrimals)
        ietiMapper.interfaceAveragesAsPrimals(mp,2);
    //! [Define primals]

    //! [Setup]
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


    if (eliminatePointwiseDofs)
        primal.setEliminatePointwiseConstraints(true);
    //! [Setup]


    // We set up the assembler
    gsOptionList assemblerOptions = gsGenericAssembler<>::defaultOptions();
    assemblerOptions.setInt("DirichletStrategy", dirichlet::elimination);
    assemblerOptions.setInt("InterfaceStrategy", iFace::dg);
    assemblerOptions.setSwitch("DG.OneSided", true);
    assemblerOptions.setReal("DG.Alpha", alpha);
    assemblerOptions.setReal("DG.Beta", beta);
    assemblerOptions.setReal("DG.Penalty", penalty);

    ///timings setup_time end
    times.setup_time += timer.stop();

    //! [Assemble]
    for (index_t k=0; k<nPatches; ++k)
    {
        //time assembling_time start
        timer.restart();

        DEBUGVAR(k);
        DEBUGVAR(nPatches);

        // We use the local variants of everything
        gsBoundaryConditions<> bc_local;
        bc.getConditionsForPatch(k,bc_local);
        gsMultiPatch<> mp_local;
        gsMultiBasis<> mb_local;
        ietiMapper.localSpaces(mp,k,mp_local,mb_local);
/*
        // only required for the yeti footprint
        if(k == 4 || k == 12 || k == 36 || k == 60)
        {
            bc_local.addCornerValue(1,0,0);
            bc_local.addCornerValue(1,0,0);
        } else if(k == 5 || k == 13)
        {
            bc_local.addCornerValue(3, 0, 0);
        } else if(k == 38 || k == 46 || k == 54 || k == 62)
        {
            bc_local.addCornerValue(2, 0, 0);
        } else if(k == 47 || k == 55)
        {
            bc_local.addCornerValue(4, 0, 0);
        }
*/

        gsGenericAssembler<> assembler(
            mp_local,
            mb_local,
            assemblerOptions,
            &bc_local
        );

        // This function provides a new dof mapper and the Dirichlet data
        // This is necessary since it might happen that a 2d-patch touches the
        // Dirichlet boundary just with a corner or that a 3d-patch touches the
        // Dirichlet boundary with a corner or an edge. These cases are not
        // covered by bc.getConditionsForPatch
        assembler.refresh(ietiMapper.augmentedDofMapperLocal(k));
        assembler.setFixedDofVector(ietiMapper.augmentedFixedPart(k));
        assembler.system().reserve(mb_local, assemblerOptions, 1);

        // Assemble
        assembler.assembleStiffness(0,false);

        assembler.assembleMoments(f,0,false);
        gsBoundaryConditions<>::bcContainer neumannSides = bc_local.neumannSides();
        for (gsBoundaryConditions<>::bcContainer::const_iterator it = neumannSides.begin();
                it!= neumannSides.end(); ++it)
            assembler.assembleNeumann(*it,false);

        for (size_t i=0; i<ietiMapper.artificialIfaces(k).size(); ++i)
        {
            patchSide side1(0,ietiMapper.artificialIfaces(k)[i].assignedTo.side());
            patchSide side2(i+1,ietiMapper.artificialIfaces(k)[i].takenFrom.side());
            boundaryInterface bi(side1, side2, mp.geoDim());
            assembler.assembleDG(bi,false);
        }

        // Fetch data
        gsSparseMatrix<real_t, RowMajor> jumpMatrix  = ietiMapper.jumpMatrix(k);
        gsMatrix<>                       localRhs    = assembler.rhs();
        gsSparseMatrix<>                 localMatrix = assembler.matrix();
        //! [Assemble]
        DEBUGMAT(jumpMatrix);
        DEBUGMAT(localRhs);
        DEBUGMAT(localMatrix);

        ///timings assembling_time end
        times.assembling_time += timer.stop();

        ///timings setup_time start
        timer.restart();

        //! [Patch to preconditioner]
        // Add the patch to the scaled Dirichlet preconditioner
        std::vector<index_t> skeletonDofs = ietiMapper.skeletonDofs(k);

        gsScaledDirichletPrec<>::Blocks blocks
                = gsScaledDirichletPrec<>::matrixBlocks(localMatrix, skeletonDofs);

        prec.addSubdomain(
                prec.restrictJumpMatrix(jumpMatrix, skeletonDofs).moveToPtr(),
                gsScaledDirichletPrec<>::schurComplement(blocks, makeSparseCholeskySolver(blocks.A11))
        );
        //! [Patch to preconditioner]

        // This function writes back to jumpMatrix, localMatrix, and localRhs,
        // so it must be called after prec.addSubdomain().
        //! [Patch to primals]
        gsSparseMatrix<> modifiedLocalMatrix, localEmbedding, embeddingForBasis;
        gsMatrix<> rhsForBasis;
        gsPrimalSystem<>::incorporateConstraints(
                ietiMapper.primalConstraints(k),
                eliminatePointwiseDofs,
                localMatrix,
                modifiedLocalMatrix,
                localEmbedding,
                embeddingForBasis,
                rhsForBasis
        );
        DEBUGMAT(localMatrix);
        DEBUGMAT(modifiedLocalMatrix);
        DEBUGMAT(localEmbedding);
        DEBUGMAT(embeddingForBasis);
        DEBUGMAT(rhsForBasis);

        gsMatrix<>                       modifiedLocalRhs     = localEmbedding * localRhs;
        gsSparseMatrix<real_t, RowMajor> modifiedJumpMatrix   = jumpMatrix * localEmbedding.transpose();

        DEBUGMAT(modifiedLocalRhs);
        DEBUGMAT(modifiedJumpMatrix);

        //! [Patch to primals]
        gsLinearOperator<>::Ptr localSolver = makeSparseLUSolver(modifiedLocalMatrix);
        DEBUGVAR(localEmbedding.rows()-localEmbedding.cols());

        gsSparseMatrix<> basisDel = gsPrimalSystem<>::primalBasis(
                localSolver,
                embeddingForBasis,
                rhsForBasis,
                ietiMapper.primalDofIndices(k),
                primal.nPrimalDofs()
        );
        DEBUGMAT(basisDel);
        //DEBUGVAR(basisDel.toDense());

        primal.addContribution(
                jumpMatrix,
                localMatrix,
                localRhs,
                basisDel
        );
        // Register the local solver to the block preconditioner. We use
        // a sparse LU solver since the local saddle point problem is not
        // positive definite.
        bdPrec->addOperator(k,k, localSolver);

        // Add the patch to the Ieti system
        //! [Patch to system]
        ieti.addSubdomain(
            modifiedJumpMatrix.moveToPtr(),
            makeMatrixOp(modifiedLocalMatrix.moveToPtr()),
            give(modifiedLocalRhs)
        );
        #ifndef NDEBUG
        //gsInfo << "Assembling fin.\n";
        #endif
        ///timings setup_time end
        times.setup_time += timer.stop();

        //! [Patch to system]
    //! [End of assembling loop]
    } // end for
    //! [End of assembling loop]

    ///timings setup_time start
    timer.restart();


    // Add the primal problem if there are primal constraints
    //! [Primal to system]
    if (ietiMapper.nPrimalDofs()>0)
    {
        // It is not required to provide a local solver to .addSubdomain,
        // since a sparse LU solver would be set up on the fly if required.
        // Here, we make use of the fact that we can use a Cholesky solver
        // because the primal problem is symmetric and positive definite:
        bdPrec->addOperator(nPatches, nPatches, makeSparseLUSolver(primal.localMatrix())); //TODO: Cholesky?

        // Add to IETI system
        ieti.addSubdomain(
            primal.jumpMatrix().moveToPtr(),
            makeMatrixOp(primal.localMatrix().moveToPtr()),
            give(primal.localRhs())
        );
    }
    //! [Primal to system]

    gsInfo << "done. " << ietiMapper.nPrimalDofs() << " primal dofs.\n";

    /**************** Setup solver and solve ****************/

    gsInfo << "Setup solver and solve... \n"
        "    Setup multiplicity scaling... " << std::flush;

    // Tell the preconditioner to set up the scaling
    //! [Setup scaling]
    prec.setupMultiplicityScaling();

    // The scaled Dirichlet preconditioner is in the last block
    gsLinearOperator<>::Ptr sdPrec = prec.preconditioner();
    bdPrec->addOperator(bdPrecSz-1,bdPrecSz-1, gsScaledOp<>::make(sdPrec,sigma));
    //! [Setup scaling]

    ///timings setup_time end
    times.setup_time += timer.stop();

    gsInfo << "done.\n    Setup gmres solver and solve... \n" << std::flush;
    // Initial guess

    ///timings solving_time start
    timer.restart();

    //! [Define initial guess]
    gsMatrix<> x;
    //TODO: Not possible: x.setRandom( bdPrec->rows(), 1 );
    x.setZero( bdPrec->rows(), 1 );

    gsInfo << "[" << bdPrec->rows() << " dofs of spp]";

    //! [Define initial guess]

    //gsInfo << ieti.rhsForSaddlePoint().transpose() << "\n";
    //gsInfo << x.transpose() << "\n";

    gsMatrix<> errorHistory;

    // This is the main cg iteration
    //! [Solve]
    //gsMinimalResidual //TODO
    gsGMRes<>( ieti.saddlePointProblem(), bdPrec )
            .setOptions( cmd.getGroup("Solver") )
            .solveDetailed( ieti.rhsForSaddlePoint(), x, errorHistory );
    //! [Solve]

    gsInfo << "done.\n    Reconstruct solution from Lagrange multipliers... " << std::flush;
    // Now, we want to have the global solution for u
    //! [Recover]
    gsMatrix<> uVec = ietiMapper.constructGlobalSolutionFromLocalSolutions(
        primal.distributePrimalSolution(
            ieti.constructSolutionFromSaddlePoint(x)
        )
    );
    //! [Recover]
    gsInfo << "done.\n\n";
    ///timings setup_time end
    times.solving_time += timer.stop();

    /******************** Print end Exit ********************/

    gsInfo << "Print timings ... \n";
    times.print(gsInfo);

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
        //fd.add(cmd);
        //fd.add(uVec);
        //gsMatrix<> mat; ieti.saddlePointProblem()->toMatrix(mat); fd.add(mat);
        fd.addComment(std::string("ietidG3_example   Timestamp:")+std::ctime(&time));
        fd.save(out);
        gsInfo << "Write solution to file " << out << "\n";
    }

    if (plot)
    {
        gsInfo << "Write Paraview data to file ieti_result.pvd\n";
        gsOptionList assemblerOptions = gsGenericAssembler<>::defaultOptions();
        assemblerOptions.setInt("DirichletStrategy", dirichlet::elimination);
        assemblerOptions.setInt("InterfaceStrategy", iFace::dg);
        gsGenericAssembler<> assembler(
            mp,
            mb,
            assemblerOptions,
            &bc
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
