/** @file biharmonic_multigrid_example.cpp

    @brief A Biharmonic example for a single patch with a multigrid solver.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): P. Weinmueller & S. Takacs
*/

# include <gismo.h>

using namespace gismo;


void setMapperForBiharmonic(const gsBoundaryConditions<> & bc, gsMappedBasis<2, real_t> & mappedBasis, gsDofMapper & mapper)
{
    mapper.setIdentity(mappedBasis.nPatches(), mappedBasis.size(), 1);
    gsMatrix<index_t> bnd;
    for (typename gsBoundaryConditions<real_t>::const_iterator
                 it = bc.begin("Dirichlet"); it != bc.end("Dirichlet"); ++it)
    {
        bnd = mappedBasis.basis(it->ps.patch).boundary(it->ps.side());
        mapper.markBoundary(it->ps.patch, bnd, 0);
    }

    for (typename gsBoundaryConditions<real_t>::const_iterator
                 it = bc.begin("Neumann"); it != bc.end("Neumann"); ++it)
    {
        bnd = mappedBasis.basis(it->ps.patch).boundaryOffset(it->ps.side(),1);
        mapper.markBoundary(it->ps.patch, bnd, 0);
    }
    mapper.finalize();
}

void setDirichletNeumannValuesL2Projection(gsMultiPatch<> & mp, gsMultiBasis<> & mb,
                                           gsMappedBasis<2, real_t> & mappedBasis,
                                           gsBoundaryConditions<> & bc, const expr::gsFeSpace<real_t> & u) {
    const gsDofMapper &mapper = u.mapper();

    gsMatrix<index_t> bnd = mapper.findFree(mapper.numPatches() - 1);
    gsDofMapper mapperBdy;
    mapperBdy.setIdentity(mappedBasis.nPatches(), mappedBasis.size(), 1);  // bb2.nPatches() == 1
    mapperBdy.markBoundary(0, bnd, 0);
    mapperBdy.finalize();

    gsExprAssembler<real_t> A(1, 1);
    A.setIntegrationElements(mb);

    auto G = A.getMap(mp);
    auto uu = A.getSpace(mappedBasis);
    auto g_bdy = A.getBdrFunction(G);

    uu.setupMapper(mapperBdy);
    gsMatrix<real_t> &fixedDofs_A = const_cast<expr::gsFeSpace<real_t> &>(uu).fixedPart();
    fixedDofs_A.setZero(uu.mapper().boundarySize(), 1);

    real_t lambda = 1e-5;

    A.initSystem();
    A.assembleBdr(bc.get("Dirichlet"), uu * uu.tr() * meas(G));
    A.assembleBdr(bc.get("Dirichlet"), uu * g_bdy * meas(G));
    A.assembleBdr(bc.get("Neumann"),
                  lambda * (igrad(uu, G) * nv(G).normalized()) * (igrad(uu, G) * nv(G).normalized()).tr() * meas(G));
    A.assembleBdr(bc.get("Neumann"),
                  lambda * (igrad(uu, G) * nv(G).normalized()) * (g_bdy.tr() * nv(G).normalized()) * meas(G));

    gsSparseSolver<real_t>::SimplicialLDLT solver;
    solver.compute(A.matrix());
    gsMatrix<real_t> &fixedDofs = const_cast<expr::gsFeSpace<real_t> & >(u).fixedPart();
    fixedDofs = solver.solve(A.rhs());
}


namespace {
    template <typename T>
    struct take_first {
        T operator() (const T& a, const T&) { return a; }
    };
}

// Copied from gsMultiBasis (not changed yet)
void combineTransferMatrices(
        const std::vector< gsSparseMatrix<real_t, RowMajor> >& localTransferMatrices,
        const gsDofMapper& coarseMapper,
        const gsDofMapper& fineMapper,
        gsSparseMatrix<real_t, RowMajor>& transferMatrix
)
{
    const index_t nBases = localTransferMatrices.size();

    index_t nonzeros = 0;
    index_t trRows = 0;
    index_t trCols = 0;

    for (index_t j=0; j<nBases; ++j)
    {
        nonzeros += localTransferMatrices[j].nonZeros();
        trRows += localTransferMatrices[j].rows();
        trCols += localTransferMatrices[j].cols();
    }

    gsSparseEntries<real_t> entries;
    entries.reserve( nonzeros );

    for (index_t j=0; j<nBases; ++j)
    {
        for (index_t k=0; k < localTransferMatrices[j].outerSize(); ++k)
        {
            for (typename gsSparseMatrix<real_t, RowMajor>::iterator it(localTransferMatrices[j],k); it; ++it)
            {
                const index_t coarse_dof_idx = coarseMapper.index(it.col(),j);
                const index_t   fine_dof_idx = fineMapper.index(it.row(),j);

                if (coarseMapper.is_free_index(coarse_dof_idx) && fineMapper.is_free_index(fine_dof_idx))
                    entries.add(fine_dof_idx, coarse_dof_idx, it.value());
            }
        }
    }

    transferMatrix.resize(fineMapper.freeSize(), coarseMapper.freeSize());
    transferMatrix.setFromTriplets(entries.begin(), entries.end(), take_first<real_t>());
    transferMatrix.makeCompressed();
}


template <typename BasisType>
void uniformCoarsen_withTransfer(
        BasisType & mBasis,
        gsSparseMatrix<real_t, RowMajor>& transferMatrix,
        const gsBoundaryConditions<real_t>& boundaryConditions,
        const gsOptionList & assemblerOptions,
        int numKnots = 1,
        index_t unk = 0
)
{
    // Get fine mapper
    gsDofMapper fineMapper;
    setMapperForBiharmonic(boundaryConditions, mBasis, fineMapper);

    // Refine
    std::vector< gsSparseMatrix<real_t, RowMajor> > localTransferMatrices(1);
    for (size_t k = 0; k < mBasis.nPatches(); ++k)
    {
        // TODO: Make projection
        mBasis.getBase(k).uniformCoarsen_withTransfer(localTransferMatrices[k], numKnots);
    }
    // Get coarse mapper
    gsDofMapper coarseMapper;
    setMapperForBiharmonic(boundaryConditions, mBasis, coarseMapper);

    // restrict to free dofs
    combineTransferMatrices( localTransferMatrices, coarseMapper, fineMapper, transferMatrix );
}

template <typename BasisType>
void buildByCoarsening(
    BasisType & mBasis,
    const gsBoundaryConditions<real_t>& boundaryConditions,
    const gsOptionList& options,
    index_t levels,
    index_t degreesOfFreedom,
    index_t unk,
    std::vector<gsSparseMatrix<real_t, RowMajor>> & transferMatrices
    )
{
    index_t lastSize = mBasis.size();

    for (index_t i = 0; i < levels-1 && lastSize > degreesOfFreedom; ++i)
    {
        gsSparseMatrix<real_t, RowMajor> transferMatrix;
        uniformCoarsen_withTransfer(mBasis, transferMatrix, boundaryConditions, options, 1, unk);

        index_t newSize = mBasis.size();
        // If the number of dofs could not be decreased, then cancel. However, if only the number
        // of levels was specified, then this should be ignored (the caller might need to have a
        // fixed number of levels).
        if (lastSize <= newSize && degreesOfFreedom > 0)
             break;
        lastSize = newSize;

        transferMatrices.push_back(give(transferMatrix));
    }

    std::reverse( transferMatrices.begin(), transferMatrices.end() );

}


template <typename BasisType>
void buildByCoarsening(
    BasisType & mBasis,
    const gsBoundaryConditions<real_t>& boundaryConditions,
    const gsOptionList& options,
        std::vector<gsSparseMatrix<real_t, RowMajor>> & transferMatrices
)
{
    return buildByCoarsening(
        mBasis,
        boundaryConditions,
        options,
        options.askInt( "Levels", 3 ),
        options.askInt( "DegreesOfFreedom", 0 ),
        0,
        transferMatrices
    );
}


int main(int argc, char *argv[])
{
    //! [Parse command line]
    bool plot = false;

    index_t numRefine  = 5;
    index_t degree = 3;
    index_t smoothness = 2;
    bool last = true;
    bool second = false;
    std::string fn;

    index_t levels = -1;
    index_t cycles = 2;
    index_t presmooth = 1;
    index_t postsmooth = 1;
    std::string smoother("GaussSeidel");
    real_t damping = -1;

    std::string iterativeSolver("cg");
    real_t tolerance = 1.e-8;
    index_t maxIterations = 100;


    gsCmdLine cmd("Example for solving the biharmonic problem (single patch only).");
    cmd.addInt("p", "degree","Set discrete polynomial degree", degree);
    cmd.addInt("", "smoothness", "Set discrete regularity",  smoothness);
    cmd.addInt("r", "refinementLoop", "Number of refinement steps",  numRefine);
    cmd.addString("f", "file", "Input geometry file (with .xml)", fn);
    //cmd.addSwitch("last", "Solve problem only on the last level of h-refinement", last);
    cmd.addSwitch("plot", "Create a ParaView visualization file with the solution", plot);

    cmd.addSwitch("second", "Compute second biharmonic problem with u = g1 and Delta u = g2 "
                            "(default first biharmonic problem: u = g1 and partial_n u = g2)", second);

    cmd.addInt   ("l", "MG.Levels",             "Number of levels to use for multigrid iteration", levels);
    cmd.addInt   ("c", "MG.NumCycles",          "Number of multi-grid cycles", cycles);
    cmd.addInt   ("",  "MG.NumPreSmooth",       "Number of pre-smoothing steps", presmooth);
    cmd.addInt   ("",  "MG.NumPostSmooth",      "Number of post-smoothing steps", postsmooth);
    cmd.addString("s", "MG.Smoother",           "Smoothing method", smoother);
    cmd.addReal  ("",  "MG.Damping",            "Damping factor for the smoother", damping);

    cmd.addString("i", "IterativeSolver",       "Iterative solver: apply multigrid directly (d) or as a preconditioner for "
                                                "conjugate gradient (cg)", iterativeSolver);
    cmd.addReal  ("t", "Solver.Tolerance",      "Stopping criterion for linear solver", tolerance);
    cmd.addInt   ("",  "Solver.MaxIterations",  "Stopping criterion for linear solver", maxIterations);


    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }
    //! [Parse command line]

    // Default case is levels:=refinements, so replace invalid default accordingly
    if (levels <0) { levels = numRefine; cmd.setInt( "MG.Levels", levels ); }
    // The smoothers know their defaults, so remove the invalid default
    if (damping<0) { cmd.remove( "MG.Damping" ); }


    //! [Read geometry]
    gsMultiPatch<> mp;
    if (fn.empty())
        mp = gsMultiPatch<>( *gsNurbsCreator<>::BSplineSquare(1,1,1) );
    else
    {
        gsInfo << "Filedata: " << fn << "\n";
        gsReadFile<>(fn, mp);
    }
    mp.clearTopology();
    mp.computeTopology();

    if (mp.nPatches() != 1)
    {
        gsInfo << "The geometry has more than one patch. Run the code with a single patch!\n";
        return EXIT_FAILURE;
    }
    //! [Read geometry]

    gsFunctionExpr<>f("256*pi*pi*pi*pi*(4*cos(4*pi*x)*cos(4*pi*y) - cos(4*pi*x) - cos(4*pi*y))",2);
    gsInfo << "Source function: " << f << "\n";

    gsFunctionExpr<> ms("(cos(4*pi*x) - 1) * (cos(4*pi*y) - 1)",2);
    gsInfo << "Exact function: " << ms << "\n";

    //! [Refinement]
    gsMultiBasis<> basis(mp, false);//true: poly-splines (not NURBS)

    // Elevate and p-refine the basis to order p + numElevate
    // where p is the highest degree in the bases
    basis.setDegree(degree); // preserve smoothness

    // h-refine each basis
    if (last)
    {
        for (index_t r =0; r < numRefine; ++r)
            basis.uniformRefine(1, degree-smoothness);
        numRefine = 0;
    }
    //! [Refinement]

    //! [Boundary condition]
    gsBoundaryConditions<> bc;
    for (gsMultiPatch<>::const_biterator bit = mp.bBegin(); bit != mp.bEnd(); ++bit)
    {
        // Laplace
        gsFunctionExpr<> laplace ("-16*pi*pi*(2*cos(4*pi*x)*cos(4*pi*y) - cos(4*pi*x) - cos(4*pi*y))",2);

        // Neumann
        gsFunctionExpr<> sol1der("-4*pi*(cos(4*pi*y) - 1)*sin(4*pi*x)",
                                 "-4*pi*(cos(4*pi*x) - 1)*sin(4*pi*y)", 2);

        bc.addCondition(*bit, condition_type::dirichlet, ms);
        if (second)
            bc.addCondition(*bit, condition_type::laplace, laplace);
        else
            bc.addCondition(*bit, condition_type::neumann, sol1der);
    }
    bc.setGeoMap(mp);
    gsInfo << "Boundary conditions:\n" << bc << "\n";
    //! [Boundary condition]

    //! [Problem setup]
    gsExprAssembler<real_t> A(1,1);
    //gsInfo<<"Active options:\n"<< A.options() <<"\n";

    // Elements used for numerical integration
    A.setIntegrationElements(basis);
    gsExprEvaluator<real_t> ev(A);

    // Set the geometry map
    auto G = A.getMap(mp);

    // Set the source term
    auto ff = A.getCoeff(f, G); // Laplace example

    // Set the discretization space
    gsSparseMatrix<> global2local;
    gsMappedBasis<2, real_t> mappedBasis;
    auto u = A.getSpace(mappedBasis);

    // Solution vector and solution variable
    gsMatrix<real_t> solVector;
    auto u_sol = A.getSolution(u, solVector);

    // Recover manufactured solution
    auto u_ex = ev.getVariable(ms, G);
    //! [Problem setup]

#ifdef _OPENMP
    gsInfo << "Available threads: "<< omp_get_max_threads() <<"\n";
#endif

    //! [Solver loop]
    gsSparseSolver<real_t>::SimplicialLDLT solver;

    gsVector<real_t> l2err(numRefine+1), h1err(numRefine+1), h2err(numRefine+1),
            dofs(numRefine+1), meshsize(numRefine+1);
    gsInfo<< "(dot1=assembled, dot2=solved, dot3=got_error)\n"
             "\nDoFs: ";
    double setup_time(0), ma_time(0), slv_time(0), err_time(0);
    gsStopwatch timer;

    // Refine uniform once
    basis.uniformRefine(1,degree -smoothness);

    const index_t r=0;

    meshsize[r] = basis.basis(0).getMaxCellLength();

    global2local.resize(basis.size(), basis.size());
    global2local.setIdentity();
    mappedBasis.init(basis,global2local);

    // Setup the Mapper
    gsDofMapper map;
    setMapperForBiharmonic(bc, mappedBasis, map);

    // Setup the system
    u.setupMapper(map);
    setDirichletNeumannValuesL2Projection(mp, basis, mappedBasis, bc, u);

    // Initialize the system
    A.initSystem();
    setup_time += timer.stop();

    gsInfo<< A.numDofs() <<std::flush;

    timer.restart();
    // Compute the system matrix and right-hand side
    A.assemble(ilapl(u, G) * ilapl(u, G).tr() * meas(G),u * ff * meas(G));

    // Enforce Laplace conditions to right-hand side
    auto g_L = A.getBdrFunction(G); // Set the laplace bdy value
    //auto g_L = A.getCoeff(laplace, G);
    A.assembleBdr(bc.get("Laplace"), (igrad(u, G) * nv(G)) * g_L.tr() );

    dofs[r] = A.numDofs();
    ma_time += timer.stop();
    gsInfo << "." << std::flush;// Assemblying done

    timer.restart();

    /**************** Setup solver and solve ****************/

    gsInfo << "Setup solver and solve... " << std::flush;

    std::vector< gsSparseMatrix<real_t,RowMajor> > transferMatrices;

    // Setup grid hiearachy by coarsening of the given matrix
    // We move the constructed hiearchy of multi bases into a variable (only required for the subspace smoother)
    // Then we move the transfer matrices into a variable
    //! [Setup grid hierarchy]
    buildByCoarsening(mappedBasis, bc, cmd.getGroup("MG"), transferMatrices);
    //! [Setup grid hierarchy]

    // Setup the multigrid solver
    //! [Setup multigrid]
    gsMultiGridOp<>::Ptr mg = gsMultiGridOp<>::make( A.matrix(), transferMatrices );
    mg->setOptions( cmd.getGroup("MG") );
    //! [Setup multigrid]

    // Since we are solving a symmetric positive definite problem, we can use a Cholesky solver
    // (instead of the LU solver that would be created by default).
    //
    // mg->matrix(0) gives the matrix for the coarsest grid level (=level 0).
    //! [Define coarse solver]
    mg->setCoarseSolver( makeSparseCholeskySolver( mg->matrix(0) ) );
    //! [Define coarse solver]

    // Set up of the smoothers
    // This has to be done for each grid level separately
    //! [Define smoothers
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
        else
        {
            gsInfo << "\n\nThe chosen smoother is unknown.\n\nKnown are:\n  Richardson (r)\n  Jacobi (j)\n  GaussSeidel (gs)"
                      "\n  IncompleteLU (ilu)\n\n";
            return EXIT_FAILURE;
        }

        smootherOp->setOptions( cmd.getGroup("MG") );

        mg->setSmoother(i, smootherOp);
    } // end for
    //! [Define smoothers]

    gsMatrix<> errorHistory;

    //! [Initial guess]
    gsMatrix<> x;
    x.setRandom( A.matrix().rows(), 1 );
    //! [Initial guess]

    //! [Solve]
    if (iterativeSolver=="cg")
        gsConjugateGradient<>( A.matrix(), mg )
            .setOptions( cmd.getGroup("Solver") )
            .solveDetailed( A.rhs(), x, errorHistory );
    else if (iterativeSolver=="d")
        gsGradientMethod<>( A.matrix(), mg )
            .setOptions( cmd.getGroup("Solver") )
            .solveDetailed( A.rhs(), x, errorHistory );
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


    //linferr[r] = ev.max( f-s ) / ev.max(f);

    //l2err[r]= math::sqrt( ev.integral( (u_ex - u_sol).sqNorm() * meas(G) ) ); // / ev.integral(ff.sqNorm()*meas(G)) );
    //h1err[r]= l2err[r] +
    //          math::sqrt(ev.integral( ( igrad(u_ex) - igrad(u_sol,G) ).sqNorm() * meas(G) )); // /ev.integral( igrad(ff).sqNorm()*meas(G) ) );

    //h2err[r]= h1err[r] +
    //          math::sqrt(ev.integral( ( ihess(u_ex) - ihess(u_sol,G) ).sqNorm() * meas(G) )); // /ev.integral( ihess(ff).sqNorm()*meas(G) )

    /*err_time += timer.stop();
    gsInfo << ". " << std::flush; // Error computations done
    //! [Solver loop]

    timer.stop();
    gsInfo << "\n\nTotal time: "<< setup_time+ma_time+slv_time+err_time <<"\n";
    gsInfo << "     Setup: "<< setup_time <<"\n";
    gsInfo << "  Assembly: "<< ma_time    <<"\n";
    gsInfo << "   Solving: "<< slv_time   <<"\n";
    gsInfo << "     Norms: "<< err_time   <<"\n";
    gsInfo << " Mesh-size: "<< meshsize.transpose() << "\n";
    gsInfo << "      Dofs: "<<dofs.transpose() << "\n";

    //! [Error and convergence rates]
    gsInfo << "\nL2 error: "<<std::scientific<<std::setprecision(3)<<l2err.transpose()<<"\n";
    gsInfo << "H1 error: "<<std::scientific<<h1err.transpose()<<"\n";
    gsInfo << "H2 error: "<<std::scientific<<h2err.transpose()<<"\n";

    if (numRefine>0)
    {
        gsInfo<< "\nEoC (L2): " << std::fixed<<std::setprecision(2)
              <<  ( l2err.head(numRefine).array()  /
                    l2err.tail(numRefine).array() ).log().transpose() / std::log(2.0)
              <<"\n";

        gsInfo<<   "EoC (H1): "<< std::fixed<<std::setprecision(2)
              <<( h1err.head(numRefine).array() /
                  h1err.tail(numRefine).array() ).log().transpose() / std::log(2.0) <<"\n";

        gsInfo<<   "EoC (H2): "<< std::fixed<<std::setprecision(2)
              <<( h2err.head(numRefine).array() /
                  h2err.tail(numRefine).array() ).log().transpose() / std::log(2.0) <<"\n";
    }
    //! [Error and convergence rates]

    //! [Export visualization in ParaView]
    if (plot)
    {
        // Write approximate and exact solution to paraview files
        gsInfo << "Plotting in Paraview...\n";
        ev.options().setSwitch("plot.elements", true);
        ev.writeParaview( u_sol   , G, "solution");
        gsInfo << "Saved with solution.pvd \n";
    }
    else
        gsInfo << "Done. No output created, re-run with --plot to get a ParaView "
                  "file containing the solution.\n";
    //! [Export visualization in ParaView]
    */
    return  EXIT_SUCCESS;
}
