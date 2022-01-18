/** @file pMultiGrid_example.cpp

    @brief The p-multigrid shows use of p-multigrid methods

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): R. Tielen, S. Takacs
*/

#include <gismo.h>
#include <string>

using namespace gismo;

gsPreconditionerOp<>::Ptr setupSubspaceCorrectedMassSmoother(
    const gsSparseMatrix<>&,
    const gsMultiBasis<>&,
    const gsBoundaryConditions<>&,
    const gsOptionList&
);

gsPreconditionerOp<>::Ptr setupBlockILUT(
    const gsSparseMatrix<>&,
    const gsMultiBasis<>&,
    const gsBoundaryConditions<>&,
    const gsOptionList&
);

gsMatrix<> assembleLumpedMass(
    const gsMultiPatch<>& mp,
    const gsMultiBasis<>& basis,
    const gsBoundaryConditions<>& bcInfo,
    const gsOptionList&
);

gsSparseMatrix<> assembleMass(
    const gsMultiPatch<>& mp,
    const gsMultiBasis<>& basis,
    const gsBoundaryConditions<>& bcInfo,
    const gsOptionList&
);

gsSparseMatrix<> assembleMixedMass(
    const gsMultiPatch<>& mp,
    const gsMultiBasis<>& basisU,
    const gsMultiBasis<>& basisV,
    const gsBoundaryConditions<>& bcInfo,
    const gsOptionList&
);

int main(int argc, char* argv[])
{
    index_t numDegree = 2;
    index_t numRefine = 6;
    index_t numSmoothing = 1;
    index_t numLevels = 2;
    index_t numBenchmark = 3;
    index_t numSplits = 0;
    index_t typeSolver = 2;
    index_t typeCycle_p = 1;
    index_t typeCycle_h = 2;
    index_t typeBCHandling = 2;
    index_t typeLumping = 1;
    index_t typeProjection = 2;
    index_t typeSmoother = 1;
    index_t typeCoarseOperator = 2;
    std::string typeCoarsening;
    index_t maxIter = 20;
    real_t tol = 1e-8;
    real_t dampingSCMS = 1;

    // Command line argument parser
    gsCmdLine cmd("This programm solves the CDR-equation with a p-multigrid or h-multigrid method");

    // Add command line arguments
    cmd.addInt("p", "Degree", "Spline degree for finest grid", numDegree);
    cmd.addInt("r", "Refinement", "Number of global refinements to obtain finest grid", numRefine);
    cmd.addInt("v", "Smoothing", "Number of pre/post smoothing steps", numSmoothing);
    cmd.addInt("l", "Levels", "Number of levels in multigrid method", numLevels);
    cmd.addInt("b", "Benchmark", "Number of the benchmark", numBenchmark);
    cmd.addInt("P", "PatchSplits", "Number of patch splittings: split each pach that many times into 2^d patches", numSplits);
    cmd.addInt("s", "Solver", "Type of solver: (1) mg as stand-alone solver (2) BiCGStab prec. with mg (3) CG prec. with mg", typeSolver);
    cmd.addInt("m", "Cycle_p", "Type of cycle where p or hp-refinement is applied: (1) V-cycle or (2) W-cycle", typeCycle_p);
    cmd.addInt("M", "Cycle_h", "Type of cycle where h-refinement is applied: (1) V-cycle or (2) W-cycle", typeCycle_h);
    cmd.addInt("d", "BCHandling", "Handles Dirichlet BC's by (1) elimination or (2) Nitsche's method", typeBCHandling);
    cmd.addInt("L", "Lumping", "Restriction and Prolongation performed with the (1) lumped or (2) consistent mass matrix", typeLumping);
    cmd.addInt("D", "Projection", "A p or z coarsening refers to (1) changing the degree directly to 1 or (2) reducing the degree by 1", typeProjection);
    cmd.addInt("S", "Smoother", "Type of smoother: (1) ILUT (2) Gauss-Seidel (3) SCMS or (4) Block ILUT", typeSmoother);
    cmd.addInt("G", "CoarseOperator", "Type of coarse operator in h-multigrid: (1) rediscretization (2) Galerkin projection", typeCoarseOperator);
    cmd.addString("z", "Coarsening", "Coarsening strategy for each of the levels: (h) h-refinement, (p) p-refinement or (z) both", typeCoarsening);
    cmd.addInt("", "MaxIter", "Maximum number of iterations for iterative solver", maxIter);
    cmd.addReal("", "Tolerance", "Threshold for iterative solver", tol);
    cmd.addReal("", "DampingSCMS", "Damping for subspace corrected mass smoother (otherwise ignored)", dampingSCMS);

    // Read parameters from command line
    try { cmd.getValues(argc,argv);  } catch (int rv) { return rv; }

    gsOptionList opt;

    // Check validity of input
    if (numDegree<=0)
    {
        gsInfo << "--Degree has to be a positive value.\n";
        return EXIT_FAILURE;
    }

    if (numRefine<0)
    {
        gsInfo << "--Refinement has to be a non-negative value.\n";
        return EXIT_FAILURE;
    }

    if (numSmoothing<0)
    {
        gsInfo << "--Smoothing has to be a non-negative value.\n";
        return EXIT_FAILURE;
    }

    if (numLevels<=0)
    {
        gsInfo << "--Levels has to be a positive value.\n";
        return EXIT_FAILURE;
    }

    if (numSplits<0)
    {
        gsInfo << "--PatchSplits has to be a non-negative value.\n";
        return EXIT_FAILURE;
    }

    if (typeCycle_p<=0 || typeCycle_h<=0)
    {
        gsInfo << "--Cycle_p and Cycle_h have to be a positvie value.\n";
        return EXIT_FAILURE;
    }

    switch (typeBCHandling)
    {
        case 1:
            opt.addInt("DirichletStrategy","",dirichlet::elimination); break;
        case 2:
            opt.addInt("DirichletStrategy","",dirichlet::nitsche    ); break;
        default:
            gsInfo << "Unknown boundary condition handling type chosen.\n";
            return EXIT_FAILURE;
    }

    if (typeLumping < 1 || typeLumping > 2)
    {
        gsInfo << "Unknown lumping type chosen.\n";
        return -1;
    }

    if (typeProjection < 1 || typeProjection > 2)
    {
        gsInfo << "Unknown projection type chosen.\n";
        return -1;
    }

    if (typeCoarseOperator < 1 || typeCoarseOperator > 2)
    {
        gsInfo << "Unknown coarse operator chosen.\n";
        return -1;
    }

    // Initialize solution, rhs and geometry
    //! [Define benchmark problem]
    gsFunctionExpr<> sol_exact, rhs_exact, coeff_diff, coeff_conv, coeff_reac;
    gsMultiPatch<> mp;
    gsInfo << "|| Benchmark information ||\n";
    switch (numBenchmark)
    {
        case 1:
            gsInfo << "1) CDR-equation the unit square\n";
            mp = gsMultiPatch<>(*gsNurbsCreator<>::BSplineSquare(1.0, 0.0, 0.0));
            sol_exact = gsFunctionExpr<>("sin(pi*x)*sin(pi*y)",2);
            rhs_exact = gsFunctionExpr<>("1.2*pi*pi*sin(pi*x)*sin(pi*y)+0.9*pi*pi*sin(pi*x)*sin(pi*y)+0.7*pi*pi*cos(pi*x)*cos(pi*y)"
                            " + 0.4*pi*pi*cos(pi*x)*cos(pi*y) +0.4*pi*cos(pi*x)*sin(pi*y)-0.2*pi*sin(pi*x)*cos(pi*y)+0.3*sin(pi*x)*sin(pi*y)", 2);
            coeff_diff = gsFunctionExpr<>("1.2","-0.7","-0.4","0.9",2);
            coeff_conv = gsFunctionExpr<>("0.4","-0.2",2);
            coeff_reac = gsFunctionExpr<>("0.3",2);
            break;

        case 2:
            gsInfo << "2) Poisson equation on the quarter annulus\n";
            mp = gsMultiPatch<>(*gsNurbsCreator<>::BSplineFatQuarterAnnulus(1.0, 2.0));
            sol_exact = gsFunctionExpr<>( "-(x*x+y*y-1)*(x*x+y*y-4)*x*y*y", 2);
            rhs_exact = gsFunctionExpr<>( "2*x*(22*x*x*y*y+21*y*y*y*y-45*y*y+x*x*x*x-5*x*x+4)", 2);
            coeff_diff = gsFunctionExpr<>("1","0","0","1",2);
            coeff_conv = gsFunctionExpr<>("0","0",2);
            coeff_reac = gsFunctionExpr<>("0",2);
            break;

        case 3:
            gsInfo << "3) Poisson equation on the quarter annulus\n";
            mp = gsMultiPatch<>(*gsNurbsCreator<>::BSplineFatQuarterAnnulus(1.0, 2.0));
            sol_exact = gsFunctionExpr<>( "(x^2+y^2-3*sqrt(x^2+y^2)+2)*sin(2*atan(y/x))", 2);
            rhs_exact = gsFunctionExpr<>( "(8-9*sqrt(x^2 + y^2))*sin(2*atan(y/x))/(x^2+y^2)", 2);
            coeff_diff = gsFunctionExpr<>("1","0","0","1",2);
            coeff_conv = gsFunctionExpr<>("0","0",2);
            coeff_reac = gsFunctionExpr<>("0",2);
            break;

        case 4:
            gsInfo << "4) Poisson equation on an L-shaped domain\n";
            mp = gsMultiPatch<>(*gsNurbsCreator<>::BSplineLShape_p1());
            sol_exact = gsFunctionExpr<>( "if ( y>0, ( (x^2+y^2)^(1.0/3.0) )*sin( (2*atan2(y,x) - pi)/3.0 ), ( (x^2+y^2)^(1.0/3.0) )*sin( (2*atan2(y,x) +3*pi)/3.0 ) )", 2);
            rhs_exact = gsFunctionExpr<>( "0", 2);
            coeff_diff = gsFunctionExpr<>("1","0","0","1",2);
            coeff_conv = gsFunctionExpr<>("0","0",2);
            coeff_reac = gsFunctionExpr<>("0",2);
            break;

        case 5:
            gsInfo << "5) Poisson equation on the unit cube\n";
            mp = gsMultiPatch<>(*gsNurbsCreator<>::BSplineCube(1));
            sol_exact = gsFunctionExpr<>( "sin(pi*x)*sin(pi*y)*sin(pi*z)", 3);
            rhs_exact = gsFunctionExpr<>( "(3*pi^2)*sin(pi*x)*sin(pi*y)*sin(pi*z)", 3);
            coeff_diff = gsFunctionExpr<>("1","0","0","0","1","0","0","0","1",3);
            coeff_conv = gsFunctionExpr<>("0","0","0",3);
            coeff_reac = gsFunctionExpr<>("0",3);
            break;

        case 6:
            gsInfo << "6) Poisson's equation on Yeti footprint\n";
            mp = *static_cast<gsMultiPatch<>::uPtr>(gsReadFile<>("domain2d/yeti_mp2.xml"));
            sol_exact = gsFunctionExpr<>( "sin(5*pi*x)*sin(5*pi*y)", 2);
            rhs_exact = gsFunctionExpr<>( "(50*pi^2)*sin(5*pi*x)*sin(5*pi*y)", 2);
            coeff_diff = gsFunctionExpr<>("1","0","0","1",2);
            coeff_conv = gsFunctionExpr<>("0","0",2);
            coeff_reac = gsFunctionExpr<>("0",2);
            break;

        default:
            gsInfo << "Unknown benchmark case chosen.\n";
            return EXIT_FAILURE;
    }
    //! [Define benchmark problem]

    // Print information about benchmark
    gsInfo << "Exact solution: " << sol_exact << "\n";
    gsInfo << "Right hand side: " << rhs_exact << "\n";
    gsInfo << "Handling of boundary conditions: " << (typeBCHandling==1 ? "elimination" : "Nitsche") << "\n";

    // Handle the uniform splitting
    //! [Splitting of patches]
    for (index_t i=0; i<numSplits; ++i)
        mp = mp.uniformSplit();
    //! [Splitting of patches]
    gsInfo << "Number of patches: " << mp.nPatches() << "\n\n";

    // Define boundary conditions
    //! [Boundary conditions]
    gsBoundaryConditions<> bcInfo;
    gsFunctionExpr<> zero("0", mp.geoDim());
    for (gsMultiPatch<>::const_biterator bit = mp.bBegin(); bit != mp.bEnd(); ++bit)
    {
        bcInfo.addCondition(*bit, condition_type::dirichlet, &sol_exact);
    }
    bcInfo.setGeoMap(mp);
    //! [Boundary conditions]

    // Setup of vector of bases
    //! [Default coarsening]
    if (typeCoarsening.empty())
    {
        // Just setup some defaults
        if (typeProjection==1)
        {
            // If we have a transition from degree p to degree 1, we only need
            // a single refinement where the degree is changed. We apply this
            // between the finest and the second finest levels.
            typeCoarsening = std::string(numLevels-2, 'h') + "z";
        }
        else //if (typeProjection==2)
        {
            // If we have a transition where the degree is changed by only one,
            // we need up to p-1 transitions. The remainder needs to by of type h.
            if (numLevels>numDegree)
                typeCoarsening = std::string(numLevels-numDegree, 'h') + std::string(numDegree-1, 'z');
            else
                // It can happen that the coarsest problem does not have degree 1 (if there
                // are not enough levels).
                typeCoarsening = std::string(numLevels-1, 'z');
        }
    }
    //! [Default coarsening]
    if (typeCoarsening.size() != (size_t)numLevels-1)
    {
        gsInfo << "The string provided to --Coarsening should have length " << numLevels-1 << "\n";
        return EXIT_FAILURE;
    }
    // Now, we count the coarsening type informations
    index_t numRefH = 0, numRefP = 0, numRefZ = 0;
    for ( index_t i = 0; i < numLevels-1 ; ++i)
    {
        if ( typeCoarsening[i] == 'h')
            numRefH++;
        else if ( typeCoarsening[i] == 'p')
            numRefP++;
        else if ( typeCoarsening[i] == 'z')
            numRefZ++;
        else
        {
            gsInfo << "Unknown type given in --Coarsening.\n";
            return EXIT_FAILURE;
        }
    }

    if (typeProjection==1 && numRefP+numRefZ > 1)
    {
        gsInfo << "If --Projection 1, there should be only one coarsening of type p or z in --Coarsening.\n";
        return EXIT_FAILURE;
    }

    if (numRefine - numRefH - numRefZ<0)
    {
        gsInfo << "Cannot have more h or z refinements in --Coarsening than overall number of refinement steps.\n";
        return EXIT_FAILURE;
    }

    // Vector of multibasis objects
    //! [Coarsest grid]
    // Index 0 refers to coarsest level; we first setup that basis
    std::vector<gsMultiBasis<>> bases;
    bases.emplace_back(mp);

    // Adjust spline degree for coarsest level
    const index_t degreeOffset = numDegree - (numRefP+numRefZ) * (typeProjection==1?(numDegree-1):1) - bases[0].degree();
    if (degreeOffset>0)
        bases[0].degreeIncrease(degreeOffset);
    else if (degreeOffset<0)
        bases[0].degreeReduce(-degreeOffset);

    // Apply refinement in h for the coarses level
    for (index_t i = 0; i < numRefine - numRefH - numRefZ; ++i)
        bases[0].uniformRefine();
    //! [Coarsest grid]

    // Now, we setup the remaining levels
    //! [Finer grids]
    for (index_t i = 1; i < numLevels; i++)
    {
        // New basis object, which is just a copy
        bases.push_back(bases.back());
        if ( typeCoarsening[i-1] == 'p' )
        {
            if (typeProjection == 1)
                bases.back().degreeIncrease(numDegree-1);
            else
                bases.back().degreeIncrease();
        }
        else if ( typeCoarsening[i-1] == 'h' )
        {
            bases.back().uniformRefine();
        }
        else //if ( typeCoarsening[i-1] == 'z' )
        {
            bases.back().uniformRefine();
            if (typeProjection == 1)
                bases.back().degreeIncrease(numDegree-1);
            else
                bases.back().degreeIncrease();
        }
    }
    //! [Finer grids]

    gsStopwatch clock;

    // Determine prolongation/restriction operators
    std::vector<gsSparseMatrix<real_t,RowMajor>> prolongation_H(numLevels-1);
    std::vector<gsLinearOperator<>::Ptr> restriction(numLevels-1);
    std::vector<gsLinearOperator<>::Ptr> prolongation(numLevels-1);

    clock.restart();
    //! [Intergrid]
    gsInfo << "|| Multigrid hierarchy ||\n";
    for (index_t i = 1; i < numLevels; i++)
    {
        if (typeCoarsening[i-1] != 'h')
        {
            if (typeLumping == 1)
            {
                gsInfo << "Restriction and prolongation between levels " << i << " and " << i-1 << ": L2-projection with lumped mass.\n";
                gsSparseMatrix<> mixedMass = assembleMixedMass(mp, bases[i], bases[i-1], bcInfo, opt);
                gsSparseMatrix<real_t> prolongationMatrix
                      = assembleLumpedMass(mp, bases[i], bcInfo, opt).asDiagonal().inverse()
                      * mixedMass;
                gsSparseMatrix<real_t> restrictionMatrix
                      = assembleLumpedMass(mp, bases[i-1], bcInfo, opt).asDiagonal().inverse()
                      * mixedMass.transpose();
                prolongation[i-1] = makeMatrixOp(prolongationMatrix.moveToPtr());
                restriction[i-1] = makeMatrixOp(restrictionMatrix.moveToPtr());
            }
            else
            {
                gsInfo << "Restriction and prolongation between levels " << i << " and " << i-1 << ": exact L2-projection.\n";
                gsSparseMatrix<> prolongationP =  assembleMixedMass(mp, bases[i], bases[i-1], bcInfo, opt);
                gsSparseMatrix<> restrictionP =  prolongationP.transpose();
                gsSparseMatrix<> prolongationM = assembleMass(mp, bases[i], bcInfo, opt);
                gsSparseMatrix<> restrictionM = assembleMass(mp, bases[i-1], bcInfo, opt);

                prolongation[i-1] = makeLinearOp(
                    [prolongationP=give(prolongationP),prolongationM=give(prolongationM)]
                    (const gsMatrix<>& Xcoarse, gsMatrix<>& Xfine)
                    {
                        gsConjugateGradient<> CGSolver(prolongationM);
                        CGSolver.setTolerance(1e-12);
                        CGSolver.solve(prolongationP*Xcoarse,Xfine);
                    },
                    prolongationP.rows(), prolongationP.cols()
                );
                restriction[i-1] = makeLinearOp(
                    [restrictionP=give(restrictionP),restrictionM=give(restrictionM)]
                    (const gsMatrix<>& Xfine, gsMatrix<>& Xcoarse)
                    {
                        gsConjugateGradient<> CGSolver(restrictionM);
                        CGSolver.setTolerance(1e-12);
                        CGSolver.solve(restrictionP*Xfine,Xcoarse);
                    },
                    restrictionP.rows(), restrictionP.cols()
                );

            }
        }
        else //if (typeCoarsening[i-1] == 'h')
        {
            gsInfo << "Restriction and prolongation between levels " << i << " and " << i-1 << ": canonical embedding and its transpose.\n";
            bases[i].clone()->uniformCoarsen_withTransfer(prolongation_H[i-1],bcInfo,opt);
            prolongation[i-1] = makeMatrixOp(prolongation_H[i-1]);
            restriction[i-1] = makeMatrixOp(prolongation_H[i-1].transpose());

        }
    }
    //! [Intergrid]
    double Time_Transfer = clock.stop();

    // Determine stiffness matrices (assembling or Galerkin projection)
    gsInfo << "\n|| Assembling ||\n";
    std::vector<gsSparseMatrix<>> matrices(numLevels);
    std::vector<gsLinearOperator<>::Ptr> operators(numLevels);
    gsMatrix<> rhs;
    double Time_Assembly = 0, Time_Galerkin = 0;
    //! [Assembling]
    for (index_t i = numLevels-1; i >= 0; --i)
    {
        if (typeCoarseOperator == 1 || i == numLevels-1 || typeCoarsening[i] != 'h')
        {
            gsInfo << "Assemble on level " << i << ": ";
            clock.restart();
            gsCDRAssembler<real_t> assembler(
                mp,
                bases[i],
                bcInfo,
                rhs_exact,
                coeff_diff, coeff_conv, coeff_reac,
                (dirichlet::strategy)opt.getInt("DirichletStrategy"),
                iFace::glue
            );
            assembler.assemble();
            matrices[i] = assembler.matrix();
            if (i==numLevels-1)
                rhs = assembler.rhs();
            Time_Assembly += clock.stop();
            gsInfo << "Degree: " << bases[i].degree()  << ", Ndof: " << matrices[i].rows() << "\n";
        }
        else
        {
            gsInfo << "Apply Galerkin projection on level " << i << ": ";
            clock.restart();
            matrices[i] = prolongation_H[i].transpose() * matrices[i+1] * prolongation_H[i];
            Time_Galerkin += clock.stop();
            gsInfo << "Ndof: " << matrices[i].rows() << "\n";
        }
        operators[i] = makeMatrixOp(matrices[i]);
    }
    gsSparseMatrix<>& matrix = matrices.back();
    //! [Assembling]

    // Setup of multigrid object
    //! [Setup]
    gsMultiGridOp<>::Ptr mg = gsMultiGridOp<>::make(operators, prolongation, restriction);
    //! [Setup]

    // Setup of solver for coarsest grid level
    clock.restart();
    //! [Coarse solver]
    mg->setCoarseSolver(makeSparseLUSolver(matrices[0]));
    //! [Coarse solver]
    double Time_Coarse_Solver_Setup = clock.stop();

    // Specify the number of cycles
    // For coarsest level, we stick to 1 in any case (exact solver is never cycled)
    //! [Set cycle]
    if (typeCycle_p==typeCycle_h)
        gsInfo << "Multigrid cycle type: " << typeCycle_h << "\n";
    else
        gsInfo << "Multigrid cycle type: " << typeCycle_h << " for h-refinement and " << typeCycle_p << " for p- and z-refinement\n";
    for (index_t i=1; i<numLevels-1; i++)
    {
        if (typeCoarsening[i] != 'h')
            mg->setNumCycles(i,typeCycle_p);
        else
            mg->setNumCycles(i,typeCycle_h);
    }
    //! [Set cycle]

    // Setup of smoothers
    gsInfo << "\n|| Smoother ||\n";
    opt.addReal("Scaling","",0.12);           // only used for SCMS
    opt.addReal("Damping","",dampingSCMS);    // only used for SCMS
    clock.restart();
    //! [Smoother]
    for (index_t i = 1; i < numLevels; i++)
    {
        if (bases[i].degree()==1)
        {
            mg->setSmoother(i,makeGaussSeidelOp(matrices[i]));
            gsInfo << "Smoother for level " << i << ": Gauss-Seidel\n";
        }
        else switch (typeSmoother)
        {
            case 1:
                mg->setSmoother(i,makeIncompleteLUOp(matrices[i]));
                gsInfo << "Smoother for level " << i << ": Incomplete LU\n";
                break;
            case 2:
                mg->setSmoother(i,makeGaussSeidelOp(matrices[i]));
                gsInfo << "Smoother for level " << i << ": Gauss-Seidel\n";
                break;
            case 3:
                mg->setSmoother(i,setupSubspaceCorrectedMassSmoother(matrices[i], bases[i], bcInfo, opt));
                gsInfo << "Smoother for level " << i << ": Subspace corrected mass smoother (damping = "<<dampingSCMS<<")\n";
                break;
            case 4:
                mg->setSmoother(i,setupBlockILUT(matrices[i], bases[i], bcInfo, opt));
                gsInfo << "Smoother for level " << i << ": Blockwise Incomplete LU\n";
                break;
            default:
                gsInfo << "Unknown smoother chosen.\n";
                return EXIT_FAILURE;
        }
    }
    gsInfo << "Number of smoothing steps: " << numSmoothing << "\n";
    mg->setNumPreSmooth(numSmoothing);
    mg->setNumPostSmooth(numSmoothing);
    // For the conjugate Gradient solver, the post smoother needs to be the transpose of the
    // pre smoother. For all other iterative solvers, it does not matter. So, we choose
    // pre smoothing and post smoothing to be the same.
    gsInfo << "Post smoothing is transpose of pre smoothing: " << (typeSolver == 3?"yes":"no") << "\n";
    mg->setSymmSmooth(typeSolver == 3);
    //! [Smoother]
    double Time_Smoother_Setup = clock.stop();

    gsInfo << "\n|| Setup Timings || \n";
    gsInfo << "Total transfer setup time: " << Time_Transfer << "\n";
    gsInfo << "Total Assembly time: " << Time_Assembly << "\n";
    gsInfo << "Total Galerkin projections: " << Time_Galerkin << "\n";
    gsInfo << "Coarse solver setup time: " << Time_Coarse_Solver_Setup << "\n";
    gsInfo << "Smoother setup time: " << Time_Smoother_Setup << "\n";
    gsInfo << "Total setup time: " << Time_Assembly + Time_Galerkin + Time_Transfer
        + Time_Coarse_Solver_Setup + Time_Smoother_Setup << "\n";

    // Setup of iterative solver
    //! [Solver]
    gsIterativeSolver<>::Ptr solver;
    switch (typeSolver)
    {
        case 1:
            gsInfo << "\n|| Solver information ||\np-multigrid is applied as stand-alone solver\n";
            // The preconditioned gradient method is noting but applying p-multigrid as a stand-alone solver.
            // We have to set step size = 1 to deactivate automatic stepsize control.
            solver = gsGradientMethod<>::make(matrix, mg, 1);
            break;
        case 2:
            gsInfo << "\n|| Solver information ||\nBiCGStab is applied as solver, p-multigrid as a preconditioner\n";
            solver = gsBiCgStab<>::make(matrix, mg);
            break;
        case 3:
            gsInfo << "\n|| Solver information ||\nCG is applied as solver, p-multigrid as a preconditioner\n";
            solver = gsConjugateGradient<>::make(matrix, mg);
            break;
        default:
            gsInfo << "Unknown iterative solver chosen.\n";
            return EXIT_FAILURE;
    }
    //! [Solver]

    //! [Initial guess]
    gsMatrix<> x = gsMatrix<>::Random(matrix.rows(),1);
    //! [Initial guess]

    // Unfortunately, the stopping criterion is relative to the rhs not to the initial residual (not yet configurable)
    const real_t rhs_norm = rhs.norm();
    solver->setTolerance( tol * (rhs-matrix*x).norm() / rhs_norm );
    solver->setMaxIterations(maxIter);
    gsMatrix<> error_history;
    clock.restart();
    //! [Solve]
    solver->solveDetailed( rhs, x, error_history );
    //! [Solve]
    double Time_Solve = clock.stop();

    gsInfo << "   0  |  Residual norm: " << std::setprecision(6) << (error_history(0,0) * rhs_norm) << "\n";
    for (index_t i=1; i<error_history.rows(); ++i)
        gsInfo << std::right << std::setw(4) << i
               << "  |  Residual norm: "
               << std::setprecision(6) << std::left << std::setw(15) << (error_history(i,0) * rhs_norm)
               << "  reduction:  1 / "
               << std::setprecision(3) << (error_history(i-1,0)/error_history(i,0))
               << "\n";

    const bool success = solver->error() <= solver->tolerance();
    if (success)
        gsInfo << "Solver reached accuracy goal after " << solver->iterations() << " iterations.\n";
    else
        gsInfo << "Solver did not reach accuracy goal within " << solver->iterations() << " iterations.\n";
    gsInfo << "The iteration took " << Time_Solve << " seconds.\n";
    // Determine residual and l2 error
    gsInfo << "Residual after solving: "  << (rhs-matrix*x).norm() << "\n";

    return success ? EXIT_SUCCESS : EXIT_FAILURE;
}

// Create the subspace corrected mass smoother
gsPreconditionerOp<>::Ptr setupSubspaceCorrectedMassSmoother(
    const gsSparseMatrix<>& matrix,
    const gsMultiBasis<>& mb,
    const gsBoundaryConditions<>& bc,
    const gsOptionList& opt
)
{
    const short_t dim = mb.topology().dim();

    // Setup dof mapper
    gsDofMapper dm;
    mb.getMapper(
       (dirichlet::strategy)opt.getInt("DirichletStrategy"),
       iFace::glue,
       bc,
       dm,
       0
    );
    const index_t nTotalDofs = dm.freeSize();

    // Decompose the whole domain into components
    std::vector< std::vector<patchComponent> > components = mb.topology().allComponents(true);
    const index_t nr_components = components.size();

    // Setup Dirichlet boundary conditions
    gsBoundaryConditions<> dir_bc;
    for ( index_t ps=0; ps < 2*dim; ++ps )
        dir_bc.addCondition( 0, 1+ps, condition_type::dirichlet, NULL );

    // Setup transfer matrices and local preconditioners
    std::vector< gsSparseMatrix<real_t,RowMajor> > transfers;
    transfers.reserve(nr_components);
    std::vector< gsLinearOperator<>::Ptr > ops;
    ops.reserve(nr_components);

    for (index_t i=0; i<nr_components; ++i)
    {
        gsMatrix<index_t> indices;
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
            if (bases[0]->dim() == dim)
            {
                GISMO_ASSERT ( bases.size() == 1, "Only one basis is expected for each patch." );
                ops.push_back(
                    gsPatchPreconditionersCreator<>::subspaceCorrectedMassSmootherOp(
                        *(bases[0]),
                        dir_bc,
                        gsOptionList(),
                        opt.getReal("Scaling")
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
    gsPreconditionerOp<>::uPtr result = gsPreconditionerFromOp<>::make(makeMatrixOp(matrix), gsAdditiveOp<>::make(transfers, ops));
    result->setOptions(opt);
    return result;
}

gsPreconditionerOp<>::Ptr setupBlockILUT(
    const gsSparseMatrix<>& A,
    const gsMultiBasis<>& mb,
    const gsBoundaryConditions<>& bc,
    const gsOptionList& opt
)
{
    const index_t nPatches = mb.nPieces();

    // Setup dof mapper
    gsDofMapper dm;
    mb.getMapper(
       (dirichlet::strategy)opt.getInt("DirichletStrategy"),
       iFace::glue,
       bc,
       dm,
       0
    );
    // Subdivide the overal dofs into
    //   a) the patch-local dofs (l=0,...,nPatches-1)
    //   b) the coupled dofs (l=nPatches)
    std::vector<index_t> sizes(nPatches+1);
    std::vector<index_t> shifts(nPatches+1);
    shifts[0] = 0;
    for (index_t k=0; k<nPatches; k++)
    {
        sizes[k] = dm.findFreeUncoupled(k).rows();
        shifts[k+1] = shifts[k] + sizes[k];
    }
    sizes[nPatches] = A.rows() - shifts[nPatches];

    // Vector of factorized operators
    std::vector< Eigen::PermutationMatrix<Dynamic,Dynamic,index_t> > P(nPatches+1);
    // Vector of factorized operators
    std::vector< Eigen::PermutationMatrix<Dynamic,Dynamic,index_t> > Pinv(nPatches+1);
    // Define the A_aprox matrix
    // Inefficient since it does not use a sparse matrix
    gsMatrix<> A_aprox(A.rows(),A.cols());
    A_aprox.setZero();

    gsSparseMatrix<> S = A.block(shifts[nPatches],shifts[nPatches],sizes[nPatches],sizes[nPatches]);

    for (index_t k=0 ; k<nPatches; k++)
    {
        // Diagonal entries
        gsSparseMatrix<> block = A.block(shifts[k],shifts[k],sizes[k],sizes[k]);
        Eigen::IncompleteLUT<real_t> ilu;
        ilu.setFillfactor(1);
        ilu.compute(block);
        P[k] = ilu.fillReducingPermutation();
        Pinv[k] = ilu.inversePermutation();
        A_aprox.block(shifts[k],shifts[k],sizes[k],sizes[k]) = ilu.factors();

        // Schur complement and off-diagonal contributions
        if (sizes[nPatches]>0)
        {
            gsMatrix<> ddB = A.block(shifts[nPatches],shifts[k],sizes[nPatches],sizes[k]);
            gsMatrix<> ddC = A.block(shifts[k],shifts[nPatches],sizes[k],sizes[nPatches]);
            gsMatrix<> ddBtilde = ilu.factors().triangularView<Eigen::Upper>().transpose().solve((ddB*P[k]).transpose()).transpose();
            gsMatrix<> ddCtilde = ilu.factors().triangularView<Eigen::UnitLower>().solve(Pinv[k]*ddC);

            S -= (ddBtilde*ddCtilde).sparseView();
            A_aprox.block(shifts[k],shifts[nPatches],sizes[k],sizes[nPatches]) = ddCtilde;
            A_aprox.block(shifts[nPatches],shifts[k],sizes[nPatches],sizes[k]) = ddBtilde;
        }
    }

    // If there is no coupling (for example if there is only one patch), the work
    // is done. We should not try to compute the ilu factorization then.
    if (sizes[nPatches]>0)
    {
        // Preform ILUT on the S-matrix.
        Eigen::IncompleteLUT<real_t> ilu;
        ilu.setFillfactor(1);
        ilu.compute(S);
        P[nPatches] = ilu.fillReducingPermutation();
        Pinv[nPatches] = ilu.inversePermutation();
        A_aprox.block(shifts[nPatches],shifts[nPatches],sizes[nPatches],sizes[nPatches]) = ilu.factors();
    }

    return gsPreconditionerFromOp<>::make(
        makeMatrixOp(A),
        makeLinearOp(
            [m_A_aprox=give(A_aprox), m_P=give(P), m_Pinv=give(Pinv),m_shifts=give(shifts),m_sizes=give(sizes),
                m_nPatches=nPatches](const gsMatrix<>& rhs, gsMatrix<>& x)
            {
                x.setZero(rhs.rows(),rhs.cols());
                for (index_t k=0; k<=m_nPatches; ++k)
                    x.block(m_shifts[k],0,m_sizes[k],1) = m_Pinv[k]*rhs.block(m_shifts[k],0,m_sizes[k],1);
                x = m_A_aprox.template triangularView<Eigen::UnitLower>().solve(x);
                x = m_A_aprox.template triangularView<Eigen::Upper>().solve(x);
                for (index_t k=0; k<=m_nPatches; ++k)
                    x.block(m_shifts[k],0,m_sizes[k],1) = m_P[k]*x.block(m_shifts[k],0,m_sizes[k],1);
            },
            A.rows(),
            A.cols()
        )
    );
}

/// @brief Determine \f$ \int_\Omega p_i dx \f$ with basis functions \f$ p_i \f$ as vector
/// The entries of this vector are the row-sums of the mass matrix
gsMatrix<> assembleLumpedMass(
    const gsMultiPatch<>& mp,
    const gsMultiBasis<>& basis,
    const gsBoundaryConditions<>& bcInfo,
    const gsOptionList& opt
)
{
    typedef gsExprAssembler<real_t>::geometryMap geometryMap;
    typedef gsExprAssembler<real_t>::variable variable;
    typedef gsExprAssembler<real_t>::space space;

    gsExprAssembler<real_t> ex(1,1);
    geometryMap G = ex.getMap(mp);
    space w_n = ex.getSpace(basis, 1, 0);
    w_n.setInterfaceCont(0);
    if ((dirichlet::strategy)opt.getInt("DirichletStrategy") == dirichlet::elimination)
    {
        w_n.setup(bcInfo, dirichlet::interpolation, 0);
    }
    ex.setIntegrationElements(basis);
    ex.initSystem();
    ex.assemble(w_n * meas(G));
    return ex.rhs();
}

/// @brief Determines the mass matrix
gsSparseMatrix<> assembleMass(
    const gsMultiPatch<>& mp,
    const gsMultiBasis<>& basis,
    const gsBoundaryConditions<>& bcInfo,
    const gsOptionList& opt
)
{
    typedef gsExprAssembler<real_t>::geometryMap geometryMap;
    typedef gsExprAssembler<real_t>::variable variable;
    typedef gsExprAssembler<real_t>::space space;

    gsExprAssembler<real_t> ex(1,1);
    geometryMap G = ex.getMap(mp);
    space w_n = ex.getSpace(basis, 1, 0);
    if ((dirichlet::strategy)opt.getInt("DirichletStrategy") == dirichlet::elimination)
    {
        w_n.setup(bcInfo, dirichlet::interpolation, 0);
    }
    ex.setIntegrationElements(basis);
    ex.initSystem();
    ex.assemble(w_n * meas(G) * w_n.tr());
    return ex.matrix();
}

/// @brief Determines the mass matrix with different bases for trial and test functions
/// The evaluation of the integrals is baed on the trial space (basisU)
gsSparseMatrix<> assembleMixedMass(
    const gsMultiPatch<>& mp,
    const gsMultiBasis<>& basisU,
    const gsMultiBasis<>& basisV,
    const gsBoundaryConditions<>& bcInfo,
    const gsOptionList& opt
)
{
    typedef gsExprAssembler<real_t>::geometryMap geometryMap;
    typedef gsExprAssembler<real_t>::variable variable;
    typedef gsExprAssembler<real_t>::space space;

    gsExprAssembler<real_t> ex(1,1);
    geometryMap G = ex.getMap(mp);
    space v_n = ex.getSpace(basisU, 1, 0);
    space u_n = ex.getTestSpace(v_n, basisV);
    if ((dirichlet::strategy)opt.getInt("DirichletStrategy") == dirichlet::elimination)
    {
        v_n.setup(bcInfo, dirichlet::interpolation, 0);
        u_n.setup(bcInfo, dirichlet::interpolation, 0);
    }
    ex.setIntegrationElements(basisU);
    ex.initSystem();
    ex.assemble(u_n * meas(G) * v_n.tr());
    return ex.matrix().transpose();
}
