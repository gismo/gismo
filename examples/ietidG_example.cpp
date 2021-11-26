/** @file ietidG_example.cpp

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

#define DEBUGVAR(v) gsInfo << #v << ": " << (v) << "\n"
#define DEBUGMAT(m) gsInfo << #m << ": " << (m).rows() << "x" << (m).cols() << "\n"

#include <gismo.h>
#include <gsAssembler/gsVisitorDg.h>
#include <gsIeti/gsArtificialIfaces.h>

using namespace gismo;


struct ctr {
   gsSparseVector<real_t> c;
   index_t k;
   index_t i;

   bool operator <( const ctr& other ) const
   {
       if (i!=other.i) return i<other.i;
       return k<other.k;
   }
};

void printPrimalConstraints( const gsIetiMapper<>& ietiMapper )
{

//c=ietiMapper.m_primalConstraints
//ii=ietiMapper.m_primalDofIndices

    std::vector<ctr> data;
    for (index_t i=0; i<ietiMapper.nPrimalDofs(); ++i)
         for (size_t j=0; j<ietiMapper.primalConstraints(i).size(); ++j)
         {
	     ctr d;
	     d.c=ietiMapper.primalConstraints(i)[j];
	     d.k=i;
	     d.i=ietiMapper.primalDofIndices(i)[j];
	     data.push_back(d);
	 }
    std::sort(data.begin(), data.end());
    for (size_t i=0; i<data.size(); ++i)
    {
        std::cout << data[i].i << " [" << data[i].k << "]" << data[i].c.transpose() << "\n";
    }
}

int main(int argc, char *argv[])
{
    /************** Define command line options *************/

    std::string geometry("domain2d/yeti_mp2.xml");
    index_t splitPatches = 1;
    real_t stretchGeometry = 1;
    index_t refinements = 1;
    index_t degree = 2;
    real_t alpha = 1;
    real_t beta = 1;
    real_t penalty = -1;
    std::string boundaryConditions("d");
    std::string primals("c");
    bool eliminateCorners = false;
    real_t tolerance = 1.e-8;
    index_t maxIterations = 100;
    bool calcEigenvalues = false;
    std::string out;
    bool plot = false;

    gsCmdLine cmd("Solves a PDE with an isogeometric discretization using an isogeometric tearing and interconnecting (IETI) solver.");
    cmd.addString("g", "Geometry",              "Geometry file", geometry);
    cmd.addInt   ("",  "SplitPatches",          "Split every patch that many times in 2^d patches", splitPatches);
    cmd.addReal  ("",  "StretchGeometry",       "Stretch geometry in x-direction by the given factor", stretchGeometry);
    cmd.addInt   ("r", "Refinements",           "Number of uniform h-refinement steps to perform before solving", refinements);
    cmd.addInt   ("p", "Degree",                "Degree of the B-spline discretization space", degree);
    cmd.addReal  ("",  "DG.Alpha",              "Parameter alpha for dG scheme; use 1 for SIPG and NIPG.", alpha );
    cmd.addReal  ("",  "DG.Beta",               "Parameter beta for dG scheme; use 1 for SIPG and -1 for NIPG", beta );
    cmd.addReal  ("",  "DG.Penalty",            "Penalty parameter delta for dG scheme; if negative, default 4(p+d)(p+1) is used.", penalty );
    cmd.addString("b", "BoundaryConditions",    "Boundary conditions", boundaryConditions);
    cmd.addString("c", "Primals",               "Primal constraints (c=corners, e=edges, f=faces)", primals);
    cmd.addSwitch("e", "EliminateCorners",      "Eliminate corners (if they are primals)", eliminateCorners);
    cmd.addReal  ("t", "Solver.Tolerance",      "Stopping criterion for linear solver", tolerance);
    cmd.addInt   ("",  "Solver.MaxIterations",  "Maximum number of iterations for linear solver", maxIterations);
    cmd.addSwitch("",  "Solver.CalcEigenvalues","Estimate eigenvalues based on Lanczos", calcEigenvalues);
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

    for (index_t i=0; i<splitPatches; ++i)
    {
        gsInfo << "split patches uniformly... " << std::flush;
        mp = mp.uniformSplit();
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

    gsInfo << "Setup bases and adjust degree... " << std::flush;

    //! [Set degree and refine]
    for ( size_t i = 0; i < mb.nBases(); ++ i )
        mb[i].setDegreePreservingMultiplicity(degree);

    for ( index_t i = 0; i < refinements; ++i )
        mb.uniformRefine();
    //! [Set degree and refine]

    gsInfo << "done.\n";

    /********* Setup assembler and assemble matrix **********/

    gsInfo << "Setup assembler and assemble matrix... " << std::flush;

    const index_t nPatches = mp.nPatches();

    // We start by setting up a global FeSpace that allows us to
    // obtain a dof mapper and the Dirichlet data
    gsPoissonAssembler<> assembler(
        mp,
        mb,
        bc,
        f,
        dirichlet::elimination,
        iFace::dg
    );
    assembler.computeDirichletDofs();

    gsInfo << "Register artificial interfaces ... "<< std::flush;
    gsArtificialIfaces<> ai( mb, assembler.system().rowMapper(0), assembler.fixedDofs() );
    ai.registerAllArtificialIfaces();
    ai.finalize();
    gsInfo << "done\n";

    gsIetiMapper<> ietiMapper( mb, ai.dofMapperMod(), assembler.fixedDofs() );

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


printPrimalConstraints(ietiMapper);

    // TODO: Jump matrices before or after primals?
    // Compute the jump matrices
    bool fullyRedundant = true,
         noLagrangeMultipliersForCorners = eliminateCorners;
    //! [Define jumps]
    ietiMapper.computeJumpMatrices(fullyRedundant, noLagrangeMultipliersForCorners);
    //! [Define jumps]

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
    //! [Setup]

    //! [Assemble]
    for (index_t k=0; k<nPatches; ++k)
    {
        // We use the local variants of everything
        gsBoundaryConditions<> bc_local;
        bc.getConditionsForPatch(k,bc_local);

        gsOptionList assemblerOptions = gsGenericAssembler<>::defaultOptions();
        assemblerOptions.setInt("DirichletStrategy", dirichlet::elimination);
        assemblerOptions.setInt("InterfaceStrategy", iFace::dg);
        assemblerOptions.setSwitch("DG.OneSided", true);
        assemblerOptions.setReal("DG.Alpha", alpha);
        assemblerOptions.setReal("DG.Beta", beta);
        assemblerOptions.setReal("DG.Penalty", penalty);

        const std::vector<gsArtificialIfaces<>::ArtificialIface>& artIfaces = ai.artificialIfaces(k);

        std::vector<gsGeometry<>* > mp_local_data;
        mp_local_data.reserve(1+2*mp.geoDim());
        std::vector<gsBasis<>* > mb_local_data;
        mp_local_data.reserve(1+2*mp.geoDim());

        mp_local_data.push_back(mp[k].clone().release());
        mb_local_data.push_back(mb[k].clone().release());
        for (size_t i=0; i<artIfaces.size(); ++i)
        {
             mp_local_data.push_back(mp[artIfaces[i].artificialIface.patch].clone().release());
             mb_local_data.push_back(mb[artIfaces[i].artificialIface.patch].clone().release());
        }
        gsMultiPatch<> mp_local(mp_local_data);
        mp_local.computeTopology();
        gsMultiBasis<> mb_local(mb_local_data,mp_local);

        gsGenericAssembler<> gAssembler(
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
        gAssembler.refresh(ai.dofMapperLocal(k));
        gsMatrix<> fixedPart(ai.dofMapperLocal(k).boundarySize(), 1 );
        fixedPart.setZero();
        fixedPart.topRows(ietiMapper.fixedPart(k).rows()) = ietiMapper.fixedPart(k); // TODO
        gAssembler.setFixedDofVector(fixedPart);
        gAssembler.system().reserve(mb_local, assemblerOptions, 1);

        // Assemble
        gAssembler.assembleStiffness(0,false);
        gAssembler.assembleMoments(f,0,false);
        gsBoundaryConditions<>::bcContainer neumannSides = bc_local.neumannSides();
        for (gsBoundaryConditions<>::bcContainer::const_iterator it = neumannSides.begin();
                it!= neumannSides.end(); ++it)
            gAssembler.assembleNeumann(*it,false);

        for (size_t i=0; i<artIfaces.size(); ++i)
        {
            patchSide side1(0,artIfaces[i].realIface.side());
            patchSide side2(i+1,artIfaces[i].artificialIface.side());
            boundaryInterface bi(side1, side2, mp.geoDim());
            gAssembler.assembleDG(bi,false);
        }

        // Fetch data
        gsSparseMatrix<real_t, RowMajor> jumpMatrix  = ietiMapper.jumpMatrix(k);
        gsMatrix<>                       localRhs    = gAssembler.rhs();
        gsSparseMatrix<>                 localMatrix = gAssembler.matrix();
        //! [Assemble]

        // Add the patch to the scaled Dirichlet preconditioner
        //! [Patch to preconditioner]
        prec.addSubdomain(
            gsScaledDirichletPrec<>::restrictToSkeleton(
                jumpMatrix,
                localMatrix,
                ietiMapper.skeletonDofs(k)
            )
        );
        //! [Patch to preconditioner]
        // This function writes back to jumpMatrix, localMatrix, and localRhs,
        // so it must be called after prec.addSubdomain().
        //! [Patch to primals]
        primal.handleConstraints(
            ietiMapper.primalConstraints(k),
            ietiMapper.primalDofIndices(k),
            jumpMatrix,
            localMatrix,
            localRhs
        );
        //! [Patch to primals]
DEBUGMAT(jumpMatrix);
        // Add the patch to the Ieti system
        //! [Patch to system]
        ieti.addSubdomain(
            jumpMatrix.moveToPtr(),
            makeMatrixOp(localMatrix.moveToPtr()),
            give(localRhs)
        );
        //! [Patch to system]
    //! [End of assembling loop]
    } // end for
    //! [End of assembling loop]

    gsInfo << "All patches are assembled\nNow handle primal system..." << std::flush;
    // Add the primal problem if there are primal constraints
    //! [Primal to system]
    if (ietiMapper.nPrimalDofs()>0)
    {
        // It is not required to provide a local solver to .addSubdomain,
        // since a sparse LU solver would be set up on the fly if required.
        // Here, we make use of the fact that we can use a Cholesky solver
        // because the primal problem is symmetric and positive definite:
        gsLinearOperator<>::Ptr localSolver
            = makeSparseCholeskySolver(primal.localMatrix());

        // Add to IETI system
        ieti.addSubdomain(
            primal.jumpMatrix().moveToPtr(),
            makeMatrixOp(primal.localMatrix().moveToPtr()),
            give(primal.localRhs()),
            localSolver
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
    //! [Setup scaling]

    gsInfo << "done.\n    Setup rhs... " << std::flush;
    // Compute the Schur-complement contribution for the right-hand-side
    //! [Setup rhs]
    gsMatrix<> rhsForSchur = ieti.rhsForSchurComplement();
    //! [Setup rhs]

    gsInfo << "done.\n    Setup cg solver for Lagrange multipliers and solve... " << std::flush;
    // Initial guess
    //! [Define initial guess]
    gsMatrix<> lambda;
    lambda.setRandom( ieti.nLagrangeMultipliers(), 1 );
    //! [Define initial guess]

    gsMatrix<> errorHistory;

    // This is the main cg iteration
    //! [Solve]
    gsConjugateGradient<> PCG( ieti.schurComplement(), prec.preconditioner() );
    PCG.setOptions( opt.getGroup("Solver") ).solveDetailed( rhsForSchur, lambda, errorHistory );
    //! [Solve]

    gsInfo << "done.\n    Reconstruct solution from Lagrange multipliers... " << std::flush;
    // Now, we want to have the global solution for u
    //! [Recover]
    gsMatrix<> uVec = ietiMapper.constructGlobalSolutionFromLocalSolutions(
        primal.distributePrimalSolution(
            ieti.constructSolutionFromLagrangeMultipliers(lambda)
        )
    );
    //! [Recover]
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

    if (calcEigenvalues)
        gsInfo << "Estimated condition number: " << PCG.getConditionNumber() << "\n";

    if (!out.empty())
    {
        gsFileData<> fd;
        std::time_t time = std::time(NULL);
        fd.add(opt);
        fd.add(uVec);
        gsMatrix<> mat; ieti.saddlePointProblem()->toMatrix(mat); fd.add(mat);
        fd.addComment(std::string("ietidG_example   Timestamp:")+std::ctime(&time));
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
