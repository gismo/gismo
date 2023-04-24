/** @file poisson_performance_example.cpp

    @brief Performance tuning based on Poisson problem

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): M. Moller

    This example is meant to benckmark the performance of assembly and
    solver for the Poisson equation.

    It supports the following command line arguments:

    -a  type of assembler (0 : expression-based, 1: visitor-based)
    -d  spatial dimension (1-4)
    -e  number of uniform order-elevation steps
    -r  number of uniform h-refinement steps
    -s  type of solver (0 : )
    -t  number of patches in t-direction
    -x  number of patches in x-direction
    -y  number of patches in y-direction
    -z  number of patches in z-direction

    --plot  Write solution to file
    --solve Solve the linear system and compute the error

    Example usage:
    ./bin/poisson_performance_example -d 2 -x 2 -y 4 -e 1 -r 1 -a 0 --solve

    This will solve a 2D Poisson equation ('-d 2') on a 2x4
    multi-patch topology ('-x 2 -y 4') with quadratic B-Splines ('-e
    1', order elevation by one) on a once refined basis ('-r 1'). The
    assembler used is gsExprAssembler ('-a 0'). This can be changed to
    the classical visitor-based assembler by using ('-a 1'). If the
    argument '--solve' is left out, only the system matrix and the
    right-hand side are assembled but the linear system is not solved.

    Comment on performance:

    So far, the assembly is parallelized over the number of degrees of
    freedoms per patch but not over the patches. 
*/

//! [Include namespace]
#include <gismo.h>

using namespace gismo;
//! [Include namespace]

int main(int argc, char *argv[])
{
    //! [Parse command line]
    bool    solve      = false;
    bool    plot       = false;
    short_t assembler  = 0;
    short_t solver     = 0;
    short_t numDim     = 3;
    index_t numRefine  = 0;
    index_t numElevate = 0;
    index_t numX       = 1;
    index_t numY       = 1;
    index_t numZ       = 1;
    index_t numT       = 1;
    
    gsCmdLine cmd("Performance tuning based on Poisson problem.");
    cmd.addInt( "a", "assembler", "Type of assembler(0: expression, 1: visitor)",  assembler );
    cmd.addInt( "d", "spatialDimension", "Number of spatial dimension", numDim );
    cmd.addInt( "e", "degreeElevation",
                "Number of degree elevation steps to perform before solving (0: equalize degree in all directions)", numElevate );
    cmd.addInt( "r", "uniformRefine", "Number of uniform h-refinement loops",  numRefine );
    cmd.addInt( "s", "solver", "Type of solver (0: Eigen-SimplicialLLT, 1: Eigen-SimplicialLDLT, 2: Eigen-QR, 3: Eigen-LU, 4: PARDISO-LLT, 5: PARDISO-LDLT, 6: PARDISO-LU, 7: Eigen-CG-Identity, 8: Eigen-CG-Diagonal, 9: Eigen-BiCGStab-Identity, 10: Eigen-BiCGStab-Diagonal, 11: Eigen-BiCGStabILUT, 12: CG, 13: BiCGStab, 14: GMRES)", solver );
    cmd.addInt( "t", "patchesTdir", "Number of patches in t-direction", numT );
    cmd.addInt( "x", "patchesXdir", "Number of patches in x-direction", numX );
    cmd.addInt( "y", "patchesYdir", "Number of patches in y-direction", numY );
    cmd.addInt( "z", "patchesZdir", "Number of patches in z-direction", numZ );

    cmd.addSwitch("solve", "Solve the linear system (default: no)", solve);
    cmd.addSwitch("plot", "Create a ParaView visualization file with the solution", plot);
    
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }
    //! [Parse command line]

    double setup_time(0), assembler_time(0), solver_time(0), error_time(0);
    gsStopwatch timer;
    timer.restart();
    
    //! [Create multi-patch topology]
    gsMultiPatch<real_t> mp;
    switch(numDim) {
    case 2:
      for (int i = 0; i < numX; ++i)
        for (int j = 0; j < numY; ++j)
          mp.addPatch(gsNurbsCreator<real_t>::BSplineSquare(1, (real_t)(i), (real_t)(j)));
      break;
    case 3:
      for (int i = 0; i < numX; ++i)
        for (int j = 0; j < numY; ++j)
          for (int k = 0; k < numZ; ++k)
            mp.addPatch(gsNurbsCreator<real_t>::BSplineCube(1, (real_t)(i), (real_t)(j), (real_t)(k)));
      break;
    case 4:
      for (int i = 0; i < numX; ++i)
        for (int j = 0; j < numY; ++j)
          for (int k = 0; k < numZ; ++k)
            for (int l = 0; l < numT; ++l)
              mp.addPatch(gsNurbsCreator<real_t>::BSplineHyperCube(1, (real_t)(i), (real_t)(j), (real_t)(k), (real_t)(l)));
      break;
    default:
      GISMO_ERROR("Unsupported spatial dimension");
    }
    mp.computeTopology();
    //! [Create multi-patch topology]

    //! [Source function]
    gsFunctionExpr<real_t> f;

    switch(numDim){
    case 2:
      f = gsFunctionExpr<real_t>("2*pi^2*sin(pi*x)*sin(pi*y)", 2);
      break;
    case 3:
      f = gsFunctionExpr<real_t>("3*pi^2*sin(pi*x)*sin(pi*y)*sin(pi*z)", 3);
      break;
    case 4:
        f = gsFunctionExpr<real_t>("4*pi^2*sin(pi*x)*sin(pi*y)*sin(pi*z)*sin(pi*w)", 4);//note: time var = w
      break;
    default:
      GISMO_ERROR("Unsupported spatial dimension");
    }
    //! [Source function]

    //! [Boundary conditions]
    gsFunctionExpr<real_t> g;

    switch(numDim){
    case 2:
      g = gsFunctionExpr<real_t>("sin(pi*x)*sin(pi*y)", 2);
      break;
    case 3:
      g = gsFunctionExpr<real_t>("sin(pi*x)*sin(pi*y)*sin(pi*z)", 3);
      break;
    case 4:
      g = gsFunctionExpr<real_t>("sin(pi*x)*sin(pi*y)*sin(pi*z)*sin(pi*w)", 4);//note: time var = w
      break;
    default:
      GISMO_ERROR("Unsupported spatial dimension");
    }
    
    gsBoundaryConditions<real_t> bc;
    for (gsMultiPatch<>::const_biterator bit = mp.bBegin(); bit != mp.bEnd(); ++bit)
      bc.addCondition( *bit, condition_type::dirichlet, &g );

    bc.setGeoMap(mp);
    //! [Boundary conditions]

    //! [Create multi-basis]
    gsMultiBasis<> mb(mp, true);
    mb.setDegree(mb.maxCwiseDegree() + numElevate);
    for (int r=0; r < numRefine; ++r)
      mb.uniformRefine();
    //! [Create multi-basis]
    
    gsInfo << "Patches: " << mp.nPatches()
           << ", degree: " << mb.minCwiseDegree()
           << ", unknowns: " << mb.totalSize()
           << "\n";
#ifdef _OPENMP
    gsInfo<< "Available threads: "<< omp_get_max_threads() <<"\n";
#endif

    // Solution vector
    gsMatrix<real_t> solVector;
    
    switch(assembler) {      
    case 0: {
      // Expression assembler

      //! [Problem setup]
      gsExprAssembler<real_t> A(1,1);
      
      typedef gsExprAssembler<real_t>::geometryMap geometryMap;
      typedef gsExprAssembler<real_t>::variable    variable;
      typedef gsExprAssembler<real_t>::space       space;
      typedef gsExprAssembler<real_t>::solution    solution;
      
      // Elements used for numerical integration
      A.setIntegrationElements(mb);
      gsExprEvaluator<> ev(A);
      
      // Set the geometry map
      geometryMap G = A.getMap(mp);
      
      // Set the discretization space
      space u = A.getSpace(mb);
      
      // Set the source term
      auto ff = A.getCoeff(f, G);
      
      // Solution vector and solution variable
      solution u_sol = A.getSolution(u, solVector);

      // Setup Dirichlet boundary conditions
      u.setup(bc, dirichlet::l2Projection, 0);

      // Initialize the system
      A.initSystem();
      setup_time += timer.stop();
      //! [Problem setup]

      //! [Assembly]
      timer.restart();
      // Compute the system matrix and right-hand side
      A.assemble(igrad(u, G) * igrad(u, G).tr() * meas(G), //matrix
                 u * ff * meas(G) ); //rhs vector
      assembler_time += timer.stop();
      //! [Assembly]

      //! [Solver]
      if (solve) {
        timer.restart();

        switch(solver) {
        case 0: { // Eigen-SimplicialLLT
          gsSparseSolver<>::SimplicialLLT solver;
          solver.compute(A.matrix());
          solVector = solver.solve(A.rhs());
          break;
        }
        case 1: { // Eigen-SimplicialLDLT
          gsSparseSolver<>::SimplicialLDLT solver;
          solver.compute(A.matrix());
          solVector = solver.solve(A.rhs());
          break;
        }
        case 2: { // Eigen-QR
          gsSparseSolver<>::QR solver;
          solver.compute(A.matrix());
          solVector = solver.solve(A.rhs());
          break;
        }
        case 3: { // Eigen-LU
          gsSparseSolver<>::LU solver;
          solver.compute(A.matrix());
          solVector = solver.solve(A.rhs());
          break;
        }
#ifdef GISMO_WITH_PARDISO
        case 4: { // PARDISO-LLT
          gsSparseSolver<>::PardisoLLT solver;
          solver.compute(A.matrix());
          solVector = solver.solve(A.rhs());
          break;
        }
        case 5: { // PARDISO-LDLT
          gsSparseSolver<>::PardisoLDLT solver;
          solver.compute(A.matrix());
          solVector = solver.solve(A.rhs());
          break;
        }
        case 6: { // PARDISO-LU
          gsSparseSolver<>::PardisoLU solver;
          solver.compute(A.matrix());
          solVector = solver.solve(A.rhs());
          break;
        }
#endif
        case 7: { // Eigen-CG-Identity
          gsSparseSolver<>::CGIdentity solver;
          solver.compute(A.matrix());
          solVector = solver.solve(A.rhs());
          break;
        }
        case 8: { // Eigen-CG-Diagonal
          gsSparseSolver<>::CGDiagonal solver;
          solver.compute(A.matrix());
          solVector = solver.solve(A.rhs());
          break;
        }
        case 9: { // Eigen-BiCGStab-Identity
          gsSparseSolver<>::BiCGSTABIdentity solver;
          solver.compute(A.matrix());
          solVector = solver.solve(A.rhs());
          break;
        }
        case 10: { // Eigen-BiCGStab-Diagonal
          gsSparseSolver<>::BiCGSTABDiagonal solver;
          solver.compute(A.matrix());
          solVector = solver.solve(A.rhs());
          break;
        }
        case 11: { // Eigen-BiCGStab-ILUT
          gsSparseSolver<>::BiCGSTABILUT solver;
          solver.compute(A.matrix());
          solVector = solver.solve(A.rhs());
          break;
        }
        case 12: { // CG
          gsConjugateGradient<> solver(A.matrix());
          solver.solve(A.rhs(), solVector);
          break;
        }
        case 13: { // BiCGStab
          gsBiCgStab<> solver(A.matrix());
          solver.solve(A.rhs(), solVector);
          break;
        }
        case 14: { // GMRES
          gsGMRes<> solver(A.matrix());
          solver.solve(A.rhs(), solVector);
          break;
        }
        default:
          GISMO_ERROR("Unsupported solver");
        }
        
        solver_time += timer.stop();
      }
      //! [Solver]

      //! [Error computation]
      if (solve) {
        timer.restart();
        
        // Set the exact solution      
        auto u_ex = A.getCoeff(g, G);
        
        // Compute L2-error
        real_t l2err = math::sqrt( ev.integral( (u_ex - u_sol).sqNorm() * meas(G) ) );
        
        real_t h1err = l2err +
          math::sqrt(ev.integral( ( igrad(u_ex) - igrad(u_sol,G) ).sqNorm() * meas(G) ));
        error_time += timer.stop();

        gsInfo << "\nL2 error: "
               << std::scientific << std::setprecision(3)
               << l2err<<"\n";
        gsInfo << "H1 error: "
               <<std::scientific << h1err << "\n";

        //! [ParaView visualization]
        if (plot) {
          ev.options().setSwitch("plot.elements", true);
          ev.writeParaview(u_sol, G, "poisson");
        }
        //! [ParaView visualization]
      }
      //! [Error computation]
     
      break;
    }
      
    case 1: {
      // Visitor-based assembler
      
      //! [Problem setup]
      gsMultiPatch<real_t> mpsol;
      gsPoissonAssembler<real_t> poisson_assembler(mp, mb, bc, f);
      poisson_assembler.options().setInt("DirichletValues", dirichlet::l2Projection);
      setup_time += timer.stop();
      //! [Problem setup]

      //! [Assembly]
      timer.restart();
      // Compute the system matrix and right-hand side
      poisson_assembler.assemble();
      assembler_time += timer.stop();
      //! [Assembly]

      //! [Solver]
      if (solve) {
        timer.restart();

        switch(solver) {
        case 0: { // Eigen-SimplicialLLT
          gsSparseSolver<>::SimplicialLLT solver;
          solver.compute(poisson_assembler.matrix());
          solVector = solver.solve(poisson_assembler.rhs());
          break;
        }
        case 1: { // Eigen-SimplicialLDLT
          gsSparseSolver<>::SimplicialLDLT solver;
          solver.compute(poisson_assembler.matrix());
          solVector = solver.solve(poisson_assembler.rhs());
          break;
        }
        case 2: { // Eigen-QR
          gsSparseSolver<>::QR solver;
          solver.compute(poisson_assembler.matrix());
          solVector = solver.solve(poisson_assembler.rhs());
          break;
        }
        case 3: { // Eigen-LU
          gsSparseSolver<>::LU solver;
          solver.compute(poisson_assembler.matrix());
          solVector = solver.solve(poisson_assembler.rhs());
          break;
        }
#ifdef GISMO_WITH_PARDISO
        case 4: { // PARDISO-LLT
          gsSparseSolver<>::PardisoLLT solver;
          solver.compute(poisson_assembler.matrix());
          solVector = solver.solve(poisson_assembler.rhs());
          break;
        }
        case 5: { // PARDISO-LDLT
          gsSparseSolver<>::PardisoLDLT solver;
          solver.compute(poisson_assembler.matrix());
          solVector = solver.solve(poisson_assembler.rhs());
          break;
        }
        case 6: { // PARDISO-LU
          gsSparseSolver<>::PardisoLU solver;
          solver.compute(poisson_assembler.matrix());
          solVector = solver.solve(poisson_assembler.rhs());
          break;
        }
#endif
        case 7: { // Eigen-CG-Identity
          gsSparseSolver<>::CGIdentity solver;
          solver.compute(poisson_assembler.matrix());
          solVector = solver.solve(poisson_assembler.rhs());
          break;
        }
        case 8: { // Eigen-CG-Diagonal
          gsSparseSolver<>::CGDiagonal solver;
          solver.compute(poisson_assembler.matrix());
          solVector = solver.solve(poisson_assembler.rhs());
          break;
        }
        case 9: { // Eigen-BiCGStab-Identity
          gsSparseSolver<>::BiCGSTABIdentity solver;
          solver.compute(poisson_assembler.matrix());
          solVector = solver.solve(poisson_assembler.rhs());
          break;
        }
        case 10: { // Eigen-BiCGStab-Diagonal
          gsSparseSolver<>::BiCGSTABDiagonal solver;
          solver.compute(poisson_assembler.matrix());
          solVector = solver.solve(poisson_assembler.rhs());
          break;
        }
        case 11: { // Eigen-BiCGStab-ILUT
          gsSparseSolver<>::BiCGSTABILUT solver;
          solver.compute(poisson_assembler.matrix());
          solVector = solver.solve(poisson_assembler.rhs());
          break;
        }
        case 12: { // CG
          gsConjugateGradient<> solver(poisson_assembler.matrix());
          solver.solve(poisson_assembler.rhs(), solVector);
          break;
        }
        case 13: { // BiCGStab
          gsBiCgStab<> solver(poisson_assembler.matrix());
          solver.solve(poisson_assembler.rhs(), solVector);
          break;
        }
        case 14: { // GMRES
          gsGMRes<> solver(poisson_assembler.matrix());
          solver.solve(poisson_assembler.rhs(), solVector);
          break;
        }
        default:
          GISMO_ERROR("Unsupported solver");
        }
        
        solver_time += timer.stop();
      }
      //! [Solver]

      //! [Error computation]
      if (solve) {
        timer.restart();

        // Construct the solution
        poisson_assembler.constructSolution(solVector, mpsol);

        gsExprAssembler<real_t> A(1,1);
        
        typedef gsExprAssembler<real_t>::geometryMap geometryMap;
        typedef gsExprAssembler<real_t>::variable    variable;
        typedef gsExprAssembler<real_t>::space       space;
        typedef gsExprAssembler<real_t>::solution    solution;

        // Elements used for numerical integration
        A.setIntegrationElements(mb);
        gsExprEvaluator<> ev(A);
        
        // Set the geometry map
        geometryMap G = A.getMap(mp);

        // Set the exact solution      
        auto u_ex = A.getCoeff(g, G);
                
        gsExprEvaluator<real_t>::variable u_sol = ev.getVariable(mpsol);

        // Compute L2-error
        real_t l2err = math::sqrt( ev.integral( (u_ex - u_sol).sqNorm() * meas(G) ) );
        
        real_t h1err = l2err +
          math::sqrt(ev.integral( ( igrad(u_ex) - igrad(u_sol,G) ).sqNorm() * meas(G) ));
        error_time += timer.stop();

        gsInfo << "\nL2 error: "
               << std::scientific << std::setprecision(3)
               << l2err<<"\n";
        gsInfo << "H1 error: "
               <<std::scientific << h1err << "\n";

        //! [ParaView visualization]
        if (plot) {
          ev.options().setSwitch("plot.elements", true);
          ev.writeParaview(u_sol, G, "poisson");
        }
        //! [ParaView visualization]
      }
      //! [Error computation]
      
      break;
    }
    default:
      GISMO_ERROR("Unsupported assembler");
    }

    //! [Output statistics]
    gsInfo << "\n\nTotal time: "<< setup_time + assembler_time + solver_time + error_time <<"\n";
    gsInfo << "     Setup: " << setup_time << "\n";
    gsInfo << " Assembler: " << assembler_time <<"\n";
    gsInfo << "    Solver: " << solver_time <<"\n";
    gsInfo << "     Error: " << error_time <<"\n";
    //! [Output statistics]
    
    return EXIT_SUCCESS;

}// end main
