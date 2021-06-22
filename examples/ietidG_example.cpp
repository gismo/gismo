/** @file ietidG_example.cpp

    @brief Provides an example for the ieti solver for a dG setting

    Here, MINRES solves the saddle point formulation. For solving
    the Schur complement formulation with CG, see ieti_example.cpp.

    This class uses the gsPoissonAssembler, for the expression
    assembler, see ieti_example.cpp.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): S. Takacs, R. Schneckenleitner
*/

#include <ctime>

#define DEBUGVAR(v) gsInfo << #v << ": " << v << "\n";
#define DEBUGMAT(m) gsInfo << #m << ": " << m.rows() << "x" << m.cols() << "\n";

#include <gismo.h>
#include <gsAssembler/gsVisitorDg.h>
#include <gsIeti/gsArtificialIfaces.h>

using namespace gismo;


void adddGInterfaceContributions(
    const gsArtificialIfaces<>& ai,
    const gsMultiPatch<>& mp,
    const gsMultiBasis<>& mb,
          gsOptionList options,
    const index_t patch,
    gsSparseMatrix<>& localMatrix,
    gsMatrix<>& localRhs
);

int main(int argc, char *argv[])
{
    /************** Define command line options *************/

    std::string geometry("domain2d/yeti_mp2.xml");
    index_t splitPatches = 1;
    real_t stretchGeometry = 1;
    index_t refinements = 1;
    index_t degree = 2;
    std::string boundaryConditions("d");
    std::string primals("c");
    real_t tolerance = 1.e-8;
    index_t maxIterations = 100;
    std::string fn;
    bool plot = false;

    gsCmdLine cmd("Solves a PDE with an isogeometric discretization using an isogeometric tearing and interconnecting (IETI) solver.");
    cmd.addString("g", "Geometry",              "Geometry file", geometry);
    cmd.addInt   ("",  "SplitPatches",          "Split every patch that many times in 2^d patches", splitPatches);
    cmd.addReal  ("",  "StretchGeometry",       "Stretch geometry in x-direction by the given factor", stretchGeometry);
    cmd.addInt   ("r", "Refinements",           "Number of uniform h-refinement steps to perform before solving", refinements);
    cmd.addInt   ("p", "Degree",                "Degree of the B-spline discretization space", degree);
    cmd.addString("b", "BoundaryConditions",    "Boundary conditions", boundaryConditions);
    cmd.addString("c", "Primals",               "Primal constraints (c=corners, e=edges, f=faces)", primals);
    cmd.addReal  ("t", "Solver.Tolerance",      "Stopping criterion for linear solver", tolerance);
    cmd.addInt   ("",  "Solver.MaxIterations",  "Maximum number of iterations for linear solver", maxIterations);
    cmd.addString("" , "fn",                    "Write solution and used options to file", fn);
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
       // overall mp keeps its topology.
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

    // Compute the jump matrices
    bool fullyRedundant = true,
         noLagrangeMultipliersForCorners = true;
    ietiMapper.computeJumpMatrices(fullyRedundant, noLagrangeMultipliersForCorners);

    // We tell the ieti mapper which primal constraints we want; calling
    // more than one such function is possible.
    ietiMapper.cornersAsPrimals();

    // The ieti system does not have a special treatment for the
    // primal dofs. They are just one more submp
    gsIetiSystem<> ieti;
    ieti.reserve(nPatches+1);

    // The scaled Dirichlet preconditioner is independent of the
    // primal dofs.
    gsScaledDirichletPrec<> prec;
    prec.reserve(nPatches);

    // Setup the primal system, which needs to know the number of primal dofs.
    gsPrimalSystem<> primal(ietiMapper.nPrimalDofs());

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
            iFace::dg
        );

        // This provides a new dof mapper and the Dirichlet data
        // This is necessary since it might happen that a 2d-patch touches the
        // Dirichlet boundary just with a corner or that a 3d-patch touches the
        // Dirichlet boundary with a corner or an edge. These cases are not
        // covered by bc.getConditionsForPatch
        gsDofMapper tmp = ietiMapper.dofMapperLocal(k);
        assembler.system() = gsSparseSystem<>(tmp); // gsSparseSystem takes per reference and steals contents
        assembler.setFixedDofVector(ietiMapper.fixedPart(k));

        // Assemble
        assembler.assemble();

        // Fetch data
        gsSparseMatrix<real_t, RowMajor> jumpMatrix  = ietiMapper.jumpMatrix(k);
        gsSparseMatrix<>                 localMatrix = assembler.matrix();
        gsMatrix<>                       localRhs    = assembler.rhs();

        gsInfo << "\n patch " << k<< "\n";
        localMatrix.uncompress();
        const index_t diff = ai.dofMapperLocal(k).freeSize() - localMatrix.rows();
        if (diff!=0)
        {
            DEBUGVAR(ai.dofMapperLocal(k).freeSize());
            DEBUGVAR(ai.dofMapperLocal(k).size());
            //DEBUGVAR(ai.dofMapperLocal(k).fixedSize());
            DEBUGVAR(ietiMapper.dofMapperLocal(k).freeSize());
            DEBUGVAR(ietiMapper.dofMapperLocal(k).size());
            //DEBUGVAR(ietiMapper.dofMapperLocal(k).fixedSize());
            DEBUGMAT(localMatrix);
            DEBUGMAT(jumpMatrix);
        }

        /*if (diff>0)
        {
            gsInfo << "Enlarge local system by " << diff << " dofs. Why??\n";
            localMatrix.conservativeResize(ai.dofMapperLocal(k).freeSize(), ai.dofMapperLocal(k).freeSize());
            localRhs.conservativeResize(ai.dofMapperLocal(k).freeSize(),1);
            localRhs.bottomRows(diff).setZero();
        }*/
        GISMO_ENSURE (diff==0, "diff: " <<diff);

        adddGInterfaceContributions(
            ai,
            mp,
            mb,
            assembler.options(),
            k,
            localMatrix,
            localRhs
        );
        localMatrix.makeCompressed();
        gsInfo << "\n matrix \n" << localMatrix.toDense() << "\n";

        // Add the patch to the scaled Dirichlet preconditioner
        prec.addSubdomain(
            gsScaledDirichletPrec<>::restrictToSkeleton(
                jumpMatrix,
                localMatrix,
                ietiMapper.skeletonDofs(k)
            )
        );

        // This function writes back to jumpMatrix, localMatrix, and localRhs,
        // so it must be called after prec.addSubmp().
        primal.handleConstraints(
            ietiMapper.primalConstraints(k),
            ietiMapper.primalDofIndices(k),
            jumpMatrix,
            localMatrix,
            localRhs
        );

        // Add the patch to the Ieti system
        ieti.addSubdomain(
            jumpMatrix.moveToPtr(),
            makeMatrixOp(localMatrix.moveToPtr()),
            give(localRhs)
        );
    }
    gsInfo << "All patches are assembled\nNow handle primal system..." << std::flush;
    // Add the primal problem if there are primal constraints
    if (ietiMapper.nPrimalDofs()>0)
    {
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

    gsInfo << "done.\n    Setup rhs... " << std::flush;
    // Compute the Schur-complement contribution for the right-hand-side
    gsMatrix<> rhsForSchur = ieti.rhsForSchurComplement();

    gsInfo << "done.\n    Setup cg solver for Lagrange multipliers and solve... " << std::flush;
    // Initial guess
    gsMatrix<> lambda;
    lambda.setRandom( ieti.nLagrangeMultipliers(), 1 );

    gsMatrix<> errorHistory;

    // This is the main cg iteration
    gsConjugateGradient<> PCG( ieti.schurComplement(), prec.preconditioner() );
    PCG.solveDetailed( rhsForSchur, lambda, errorHistory );

    gsInfo << "done.\n    Reconstruct solution from Lagrange multipliers... " << std::flush;
    // Now, we want to have the global solution for u
    gsMatrix<> uVec = ietiMapper.constructGlobalSolutionFromLocalSolutions(
        primal.distributePrimalSolution(
            ieti.constructSolutionFromLagrangeMultipliers(lambda)
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

    if (!fn.empty())
    {
        gsFileData<> fd;
        std::time_t time = std::time(NULL);
        fd.add(opt);
        fd.add(uVec);
        fd.addComment(std::string("ietidG_example   Timestamp:")+std::ctime(&time));
        fd.save(fn);
        gsInfo << "Write solution to file " << fn << "\n";
    }
    /*
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
    */
    if (!plot&&fn.empty())
    {
        gsInfo << "Done. No output created, re-run with --plot to get a ParaView "
                  "file containing the solution or --fn to write solution to xml file.\n";
    }
    return success ? EXIT_SUCCESS : EXIT_FAILURE;
}


void adddGInterfaceContributions(
    const gsArtificialIfaces<>& ai,
    const gsMultiPatch<>& mp,
    const gsMultiBasis<>& mb,
          gsOptionList options,
    const index_t patch,
    gsSparseMatrix<>& localMatrix,
    gsMatrix<>& localRhs
)
{
    options.addSwitch("DG.NonSymmetric", "", true);
    const std::vector<gsArtificialIfaces<>::ArtificialIface>& artIfaces = ai.artificialIfaces(patch);

    gsDofMapper localMapperTempCopy = ai.dofMapperLocal(patch);
    gsSparseSystem<real_t> sparseSystem (localMapperTempCopy); // Yummi. Mapper tastes good.
    sparseSystem.reserve(7,1); // TODO: first parameter should be number of non-zeros per row

    std::vector< gsMatrix<> > fixedPart(1, ai.fixedDofs(patch) );

    gsMatrix<> quNodes1, quNodes2;// Mapped nodes
    gsVector<> quWeights;         // Mapped weights

    gsQuadRule<real_t> QuRule;

    for (size_t l = 0; l < artIfaces.size(); ++l)
    {

        const patchSide& side1     = artIfaces[l].realIface;
        const patchSide& side2     = artIfaces[l].artificialIface;
        const gsBasis<>& basis1    = mb[side1.patch];
        const gsBasis<>& basis2    = mb[side2.patch];
        const gsGeometry<>& geo1   = mp[side1.patch];
        const gsGeometry<>& geo2   = mp[side2.patch];

        gsVisitorDg<real_t> dg;

        boundaryInterface bi(side1, side2, mp.domainDim());
        gsRemapInterface<real_t> interfaceMap(mp, mb, bi);

        // Initialize
        dg.initialize(basis1, basis2, bi, options, QuRule);

        // Initialize mp element iterators
        // TODO: What should we do to tell the system which basis is to be used for iterating?
        typename gsBasis<>::domainIter domIt = interfaceMap.makeDomainIterator();

        // The following iterators are required such that the visitor knows the correct penality term
        typename gsBasis<>::domainIter domIt_helper1 = basis1.makeDomainIterator(bi.first().side());
        typename gsBasis<>::domainIter domIt_helper2 = basis2.makeDomainIterator(bi.second().side());

        // iterate over all boundary grid cells on the "left"
        for (; domIt->good(); domIt->next() )
        {
            QuRule.mapTo( domIt->lowerCorner(), domIt->upperCorner(), quNodes1, quWeights);
            interfaceMap.eval_into(quNodes1,quNodes2);

            // Perform required evaluations on the quadrature nodes
            dg.evaluate(basis1, geo1, basis2, geo2, quNodes1, quNodes2);

            // Assemble on element
            dg.assemble(*domIt_helper1, *domIt_helper2, quWeights);
            // Here, the iterators are only used for the calculation of the cell size

            // Map it back
            dg.localToGlobal(0, l+1, fixedPart, sparseSystem);

        }

    }

    localRhs +=       sparseSystem.rhs();
    localMatrix +=    sparseSystem.matrix();

}
