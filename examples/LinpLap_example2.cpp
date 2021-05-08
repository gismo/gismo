/** @file poisson_example.cpp
brief Tutorial on how to use G+Smo to solve the Poisson equation,
see the \ref PoissonTutorial
This file is part of the G+Smo library.
This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this
file, You can obtain one at http://mozilla.org/MPL/2.0/.
Author(s):
*/

//! [Include namespace]
#include<string>
#include <gismo.h>
#include <ctime>

using namespace gismo;
//! [Include namespace]

gsMatrix<real_t> projectL2(const gsMultiPatch<real_t> &mp, const gsMultiBasis<real_t> &mb, const gsFunction<real_t> &g);
gsMatrix<real_t> addDirVal(const gsAssembler<real_t> &a, const gsMatrix<real_t> &solVector);
gsMatrix<real_t> reduceDirichlet(const gsAssembler<real_t> &a, const gsMatrix<real_t> &w);

int main(int argc, char* argv[])
{
    real_t eps = 1;
    real_t p = 1.8;            //p-Laplace parameter
    index_t problemId = 2;
    index_t k = 1;            //spline degree
    index_t reduceCont = 0;
    index_t maxiter = 100;
    real_t tol = 1e-6;        //residual error tolerance
    real_t tol0 = 0;
    index_t num = 8;        //number of refinements
    index_t str = 2;
    bool require_fin = false;
    index_t subdiv = 1;
    index_t bc = 1;
    index_t startrefine = 1;
    real_t lambda = 0;
    real_t alpha = p-2; //default Eigenvalue problem (if f=0)
    real_t gamma = 1;
    index_t initial = 1;
    std::string solver_type("lu");

    gsCmdLine cmd("Linearized p-Laplace example");

    cmd.addReal("e", "eps", "variable for eps", eps);
    cmd.addReal("p", "pow", "p-Laplace Parameter", p);
    cmd.addInt("", "problem", "which problem to solve", problemId);
    cmd.addInt("k", "degree", "degree of basis", k);
    cmd.addInt("c", "reduceCont", "reduce continuity of basis by...", reduceCont);
    cmd.addInt("i", "maxiter", "maximal iterations", maxiter);
    cmd.addReal("t", "tol", "Residual error tolerance", tol);
    cmd.addReal("", "tol0", "Residual error tolerance for coarsest grid (if 0: like --tol)", tol0);
    cmd.addInt("r", "numRefine", "number of refinements of the mesh", num);
    cmd.addInt("s", "strat", "Method for Dirichlet Imposition", str);
    cmd.addInt("", "sub", "Number of subdivision of an element for quadrature nodes", subdiv);
    cmd.addSwitch("", "fin", "After computation, wait until button is pressed", require_fin);
    cmd.addInt("", "bc", "Type of Boundary conditions on all boundaries", bc);
    cmd.addInt("", "startrefine", "Number of refinement steps before entering nested iteration", startrefine);
    cmd.addReal("","lambda","Parameter for lambda",lambda);
    cmd.addReal("","alpha","Parameter for alpha",alpha);
    cmd.addReal("","gamma","Parameter for gamma",gamma);
    cmd.addInt("","initial","Choice for initial guess u_0",initial);
    cmd.addString("", "solver", "Solver to be used (lu, cg, cg-mg, gmres, gmres-mg)", solver_type);
    try { cmd.getValues(argc, argv); }
    catch (int rv) { return rv; }

    if (tol0==0) tol0 = tol;

    gsInfo << "Printing command line arguments:\n"
          << "eps               = " << eps << "\n"
          << "pow               = " << p << "\n"
          << "problem           = " << problemId << "\n"
          << "degree            = " << k << "\n"
          << "reduceCont        = " << reduceCont << "\n"
          << "maxiter           = " << maxiter << "\n"
          << "tol               = " << tol << "\n"
          << "tol0              = " << tol0 << "\n"
          << "numRefine         = " << num << "\n"
          << "subdiv            = " << subdiv << "\n"
          << "solver_type       = " << solver_type << "\n"
          << "BC                = " << bc << "\n";
    gsOptionList opt = gsAssembler<>::defaultOptions();
    //opt.setInt("DirichletValues", dirichlet::l2Projection);

    if (str == 1)
    {
        gsInfo << "DirichletStrategy = dirichlet::elimination\n\n";
        opt.setInt("DirichletStrategy", dirichlet::elimination);
    }
    else if (str == 2)
    {
        gsInfo << "DirichletStrategy = dirichlet::nitsche\n\n";
        opt.setInt("DirichletStrategy", dirichlet::nitsche);
    }
    else
    {
        gsInfo << "DirichletStrategy unknown.";
        return 1;
    }
    
   	gsOptionList opt2 = opt;
	  opt2.setInt("DirichletStrategy", dirichlet::nitsche);

    real_t e_0 = 0;
    real_t e_F = 0;
    real_t e_0old;
    real_t e_Fold;

    real_t Lp_rate;
    real_t F_rate;

    //! [Function data]

    real_t omega = 2;
    index_t l = 3;

    gsFunctionExpr<> u; // Exact solution
    gsFunctionExpr<> u_derEast;
    gsFunctionExpr<> u_derWest;
    gsFunctionExpr<> u_derNorth;
    gsFunctionExpr<> u_derSouth;
    gsFunctionExpr<> f; // Source function

    if (problemId == 1)
    {
        f = 	gsFunctionExpr<>("-((" + std::to_string(gamma) + "^2*(x^2 + y^2)^(" + std::to_string(gamma) + "/2)*(" + std::to_string(eps) + "^2 + " + std::to_string(gamma) + "^2*(x^2 + y^2)^(-1 + " + std::to_string(gamma) + "))^(" + std::to_string(p) + "/2)*(" + std::to_string(eps) + "^2*(x^2 + y^2) + " + std::to_string(gamma) + "*(2 + " + std::to_string(gamma) + "*(-1 + " + std::to_string(p) + ") - " + std::to_string(p) + ")*(x^2 + y^2)^" + std::to_string(gamma) + "))/(" + std::to_string(eps) + "^2*(x^2 + y^2) + " + std::to_string(gamma) + "^2*(x^2 + y^2)^" + std::to_string(gamma) + ")^2)", 2);
        u = gsFunctionExpr<>("(x^2+y^2)^(" + std::to_string(gamma) + "/2)", 2);

    }
    else if (problemId == 2)
    {
        f	= gsFunctionExpr<>("2*" + std::to_string(omega) + "^2*pi^2*(" + std::to_string(eps*eps) + "+2*" + std::to_string(omega) + "^2*pi^2*cos(" + std::to_string(omega) + "*pi*(x+y))^2)^((" + std::to_string(p) + "-4)/2)*(" + std::to_string(eps*eps) + "+2*" + std::to_string(omega) + "^2*(" + std::to_string(p) + "-1)*pi^2*cos(" + std::to_string(omega) + "*pi*(x+y))^2)*sin(" + std::to_string(omega) + "*pi*(x+y))+ " + std::to_string(lambda) + "*abs(sin(" + std::to_string(omega) + "*pi*(x+y)))^" + std::to_string(alpha) +"*sin(" + std::to_string(omega) + "*pi*(x+y))", 2);
        
        u =	gsFunctionExpr<>("sin(" + std::to_string(omega) + "*pi*(x+y))", 2);
        
        u_derEast = gsFunctionExpr<>("(" + std::to_string(eps*eps) + "+2*cos(" + std::to_string(omega) + "*pi*(x+y))^2*" + std::to_string(omega) + "^2*pi^2)^((" + std::to_string(p) + "-2)/2)*(" + std::to_string(omega) + "*pi*cos(pi*" + std::to_string(omega) + "*(x+y)))", 2);
        u_derWest = gsFunctionExpr<>("(" + std::to_string(eps*eps) + "+2*cos(" + std::to_string(omega) + "*pi*(x+y))^2*" + std::to_string(omega) + "^2*pi^2)^((" + std::to_string(p) + "-2)/2)*(-" + std::to_string(omega) + "*pi*cos(pi*" + std::to_string(omega) + "*(x+y)))", 2);
        u_derNorth = gsFunctionExpr<>("(" + std::to_string(eps*eps) + "+2*cos(" + std::to_string(omega) + "*pi*(x+y))^2*" + std::to_string(omega) + "^2*pi^2)^((" + std::to_string(p) + "-2)/2)*(" + std::to_string(omega) + "*pi*cos(pi*" + std::to_string(omega) + "*(x+y)))", 2);
        u_derSouth = gsFunctionExpr<>("(" + std::to_string(eps*eps) + "+2*cos(" + std::to_string(omega) + "*pi*(x+y))^2*" + std::to_string(omega) + "^2*pi^2)^((" + std::to_string(p) + "-2)/2)*(-" + std::to_string(omega) + "*pi*cos(pi*" + std::to_string(omega) + "*(x+y)))", 2);
        
    }
    else if (problemId == 3)
    {
        f = gsFunctionExpr<>("8*pi^2*(" + std::to_string(eps*eps) + "+2*pi^2+pi^2*(-(" + std::to_string(p) + "-2)*cos(4*pi*y)-cos(4*pi*x)*(" + std::to_string(p) + "-2+2*(" + std::to_string(p) + "-1)*cos(4*pi*y))))*(" + std::to_string(eps*eps) + "+2*pi^2-pi^2*(cos(4*pi*(x-y))+cos(4*pi*(x+y))))^((" + std::to_string(p) + "-4)/2)*(sin(2*pi*x)*sin(2*pi*y))+" + std::to_string(lambda) + "*abs(sin(2*pi*x)*sin(2*pi*y))^" + std::to_string(alpha) +"*sin(2*pi*x)*sin(2*pi*y)", 2);
        
        u = gsFunctionExpr<>("sin(2*pi*x)*sin(2*pi*y)", 2);
    }
    else if (problemId == 4)
    {
        f = gsFunctionExpr<>("(" + std::to_string(eps*eps) + "+cos(x)^2)^(" + std::to_string(p) + "/2-2)*(" + std::to_string(eps*eps) + "+(" + std::to_string(p) + "-1)*cos(x)^2)*sin(x) + " + std::to_string(lambda) + "* abs(sin(x))^" + std::to_string(alpha) +"*sin(x)", 2);
        
        u = gsFunctionExpr<>("sin(x)", 2);
    }
    else if (problemId == 5)
    {
        f = gsFunctionExpr<>("1", 2);
    }
    else
    {
        gsInfo << "problemId unknown.";
        return 1;
    }

    gsFunctionExpr<> u0;
  
    gsFunctionExpr<> Z("0",2);
  
    if(initial==1)
    {
	      u0 = gsFunctionExpr<>("0", 2);
    }
    else if(initial==2)
    {
        u0 = gsFunctionExpr<>("x+y", 2);
    }
    else if(initial==3)
    {
        u0 = gsFunctionExpr<>("(" + u.expression() + ")*exp(x*(1-x)*y*(1-y))", 2);
    }

    // Print out source function and solution
    gsInfo << "Source function " << f << "\n";
    gsInfo << "Exact solution " << u << "\n";
    gsInfo << "Initial guess" << u0 << "\n\n";
    //! [Function data]


    //! [Geometry data]
   	gsMultiPatch<> patch = gsMultiPatch<>(*gsNurbsCreator<real_t>::BSplineSquare(1,0,0));
    //! [Geometry data]

    //! [Boundary conditions]
   	gsBoundaryConditions<> bcInfo;
	  gsBoundaryConditions<> hbcInfo;
    if (bc == 1)
	  {
		    bcInfo.addCondition(0, boundary::west, condition_type::dirichlet, &u);
		    bcInfo.addCondition(0, boundary::east, condition_type::dirichlet, &u);
		    bcInfo.addCondition(0, boundary::north, condition_type::dirichlet, &u);
		    bcInfo.addCondition(0, boundary::south, condition_type::dirichlet, &u);

		    hbcInfo.addCondition(0, boundary::west, condition_type::dirichlet, &Z);
		    hbcInfo.addCondition(0, boundary::east, condition_type::dirichlet, &Z);
		    hbcInfo.addCondition(0, boundary::north, condition_type::dirichlet, &Z);
		    hbcInfo.addCondition(0, boundary::south, condition_type::dirichlet, &Z);
	  }
	  else if (bc==2)
	  {
		    bcInfo.addCondition(0, boundary::west, condition_type::neumann, &u_derWest);
		    bcInfo.addCondition(0, boundary::east, condition_type::neumann, &u_derEast);
		    bcInfo.addCondition(0, boundary::north, condition_type::neumann, &u_derNorth);
		    bcInfo.addCondition(0, boundary::south, condition_type::neumann, &u_derSouth);
		    //bcInfo.addCornerValue(boundary::southwest, 0);

		    hbcInfo.addCondition(0, boundary::west, condition_type::neumann, &Z);
		    hbcInfo.addCondition(0, boundary::east, condition_type::neumann, &Z);
		    hbcInfo.addCondition(0, boundary::north, condition_type::neumann, &Z);
		    hbcInfo.addCondition(0, boundary::south, condition_type::neumann, &Z);
		    //hbcInfo.addCornerValue(boundary::southwest,0);
    }
    else
    {
        bcInfo.addCondition(0, boundary::west, condition_type::dirichlet, &u);
		    bcInfo.addCondition(0, boundary::east, condition_type::neumann, &u_derEast);
		    bcInfo.addCondition(0, boundary::north, condition_type::neumann, &u_derNorth);
		    bcInfo.addCondition(0, boundary::south, condition_type::dirichlet, &u);
		    //bcInfo.addCornerValue(boundary::southwest, 0);

		    hbcInfo.addCondition(0, boundary::west, condition_type::dirichlet, &Z);
		    hbcInfo.addCondition(0, boundary::east, condition_type::neumann, &Z);
		    hbcInfo.addCondition(0, boundary::north, condition_type::neumann, &Z);
		    hbcInfo.addCondition(0, boundary::south, condition_type::dirichlet, &Z);
		    //hbcInfo.addCornerValue(boundary::southwest,0);
    }
    //! [Boundary conditions]

    //! [Refinement]
    // Copy basis from the geometry

    gsMultiBasis<> basis(patch);
    for (index_t i = 0; i < startrefine; i++)
          basis.uniformRefine();
    basis[0].setDegreePreservingMultiplicity(k);
    basis.reduceContinuity(reduceCont);
    
    // Setup for solving
    if (str == 2) { hbcInfo = bcInfo; }
    gsMatrix<real_t> w = projectL2(patch, basis, u0);  // initial guess on coarsest grid
    gsLinpLapPde<real_t> pde(patch, bcInfo, f, eps, p, w, lambda, alpha);

    gsInfo << "eps = " << eps << " , p = " << p << " , k = " << k << " , lambda = " << lambda << "\n";
    gsInfo << "Dofs       &CPU time  &L_p error& L_p rate  &F error  & F rate    &iter      &p         &eps       &||rh||_{\ell^2} \n";

    gsMatrix<> solVector;
    
    gsLinpLapAssembler<real_t> A;
    A.initialize(pde, eps, basis, opt, subdiv);

    std::vector< gsSparseMatrix<real_t, RowMajor> > transfers;
    for (index_t i = startrefine; i < num; i++)
    {
        gsSparseMatrix<real_t, RowMajor> transfer;
        basis.uniformRefine_withTransfer(transfer, bcInfo, opt2);
        pde.w = transfer * pde.w; // get initial guess from coaeser grid // Here, the Dirichlet values are not the exact L2-projection -- ST
        transfers.push_back(give(transfer));

        std::clock_t c_start = std::clock();

        gsLinpLapAssembler<real_t> A;
        A.initialize(pde, eps, basis, opt, subdiv);
        A.assemble();

        solVector = reduceDirichlet(A, pde.w);

        gsMatrix<> rh = A.matrix() * solVector - A.rhs();
        index_t iter = 0;

        real_t final_res_sq = (rh.transpose()*rh).value();
        real_t initial_res_sq = final_res_sq;

        real_t tol_k = i==startrefine ? tol0 : tol;

/*
        gsInfo << "   " << 0 << " (initial res "
               << std::setprecision(4) << std::scientific << final_res_sq << ")\n";
*/

        do
        {
            index_t ls_iter = 0;
            real_t cg_cond = 0;
            if (solver_type=="lu")
            {
                gsSparseSolver<>::LU solver(A.matrix());
                solVector += solver.solve(-rh);
            } else if (solver_type=="cg") {
                gsConjugateGradient<> solver(A.matrix());
                solver.setCalcEigenvalues(true);
                //solVector.setZero();
                solver.solve(A.rhs(), solVector);
                ls_iter = solver.iterations();
                cg_cond = solver.getConditionNumber();
            } else if (solver_type=="cg-mg") {
                gsMultiGridOp<>::Ptr mg=gsMultiGridOp<>::make(A.matrix(), transfers);
                for (index_t i=1; i<mg->numLevels(); ++i)
                    mg->setSmoother(i,makeGaussSeidelOp(mg->matrix(i)));
                gsConjugateGradient<> solver(A.matrix());
                solver.setCalcEigenvalues(true);
                //solVector.setZero();
                solver.solve(A.rhs(), solVector);
                ls_iter = solver.iterations();
                cg_cond = solver.getConditionNumber();
            } else if (solver_type=="gmres") {
                gsGMRes<> solver(A.matrix());
                //solVector.setZero();
                solver.solve(A.rhs(), solVector);
                ls_iter = solver.iterations();
            } else if (solver_type=="gmres-mg") {
                gsMultiGridOp<>::Ptr mg=gsMultiGridOp<>::make(A.matrix(), transfers);
                for (index_t i=1; i<mg->numLevels(); ++i)
                    mg->setSmoother(i,makeGaussSeidelOp(mg->matrix(i)));
                gsGMRes<> solver(A.matrix());
                //solVector.setZero();
                solver.solve(A.rhs(), solVector);
                ls_iter = solver.iterations();
            } else {
                GISMO_ENSURE(0, "Unknown solver");
            }

            pde.w = addDirVal(A, solVector); // Add Dirichlet values to current solution and set as new w. // Here, the exact L2-projection of the Dirichlet values is used -- ST
            A.initialize(pde, eps, basis, opt);
            A.assemble();
            rh = A.matrix() * solVector - A.rhs();
            final_res_sq = (rh.transpose()*rh).value();

            iter++;

/*
            if (solver_type=="lu")
                gsInfo << "   " << iter << " (lin solver: obtained res "
                       << std::setprecision(4) << std::scientific << final_res_sq << ")\n";
            else if (solver_type=="cg"||solver_type=="cg-mg")
                gsInfo << "   " << iter << " (lin solver: " << ls_iter << " iterations to reach res "
                       << std::setprecision(4) << std::scientific << final_res_sq << ". "
                       << " cond: " << cg_cond << ")\n";
            else if (solver_type=="gmres"||solver_type=="gmres-mg")
                gsInfo << "   " << iter << " (lin solver: " << ls_iter << " iterations to reach res "
                       << std::setprecision(4) << std::scientific << final_res_sq << ")\n";
*/
        } while (iter < maxiter && final_res_sq>tol_k*tol_k*initial_res_sq);

        std::clock_t c_end = std::clock();
        double time = (c_end - c_start) / CLOCKS_PER_SEC;

        gsMultiPatch<> mpsol;
        A.constructSolution(solVector, mpsol); //construct solution from the free DoFs via the assembler that is set to elimination.
        gsField<> sol(A.patches(), mpsol);

        e_0old = e_0;
        e_Fold = e_F;

        e_0 = sol.distanceLp(u, basis, p, false);
        e_F = sol.distanceF(u, basis, eps, p, false);

        if (i == startrefine)
        {
            gsInfo << std::setw(8) << A.numDofs() << "   &"
              << std::setprecision(2) << std::fixed << std::setw(7) << time << "s  &"
              << std::setprecision(2) << std::scientific << std::setw(8) << e_0 << " &  "
              << "   -     &"
              << std::setprecision(2) << std::scientific << std::setw(8) << e_F << " & "
              << "   -      &"
              << std::setw(8) << iter << "  &"
              << std::setprecision(2) << std::fixed << std::setw(8) << p << "  &"
              << std::setprecision(2) << std::fixed << std::setw(8) << eps << "  &"
              << std::setprecision(2) << std::scientific << std::setw(8) << math::sqrt(final_res_sq) << "\n";
        }
        else
        {
            Lp_rate = math::log(e_0 / e_0old) / math::log(0.5);
            F_rate = math::log(e_F / e_Fold) / math::log(0.5);
            gsInfo << std::setw(8) << A.numDofs() << "   &"
              << std::setprecision(2) << std::fixed << std::setw(7) << time << "s  &"
              << std::setprecision(2) << std::scientific << std::setw(8) << e_0 << " &"
              << std::setprecision(4) << std::fixed << std::setw(9) << Lp_rate << "  &"
              << std::setprecision(2) << std::scientific << std::setw(8) << e_F << " &"
              << std::setprecision(4) << std::fixed << std::setw(9) << F_rate << "  &"
              << std::setw(8) << iter << "  &"
              << std::setprecision(2) << std::fixed << std::setw(8) << p << "  &"
              << std::setprecision(2) << std::fixed << std::setw(8) << eps << "  &"
              << std::setprecision(2) << std::scientific << std::setw(8) << math::sqrt(final_res_sq) << "\n";
        }
    }

    if (require_fin)
    {
        gsInfo << "fin";
        std::cin.get();
    }
    return EXIT_SUCCESS;

}// end main


//-------------------------------------------------------------------------------------------------------------------------------

gsMatrix<real_t> projectL2(const gsMultiPatch<real_t> &mp, const gsMultiBasis<real_t> &mb, const gsFunction<real_t> &g)
{
    gsGenericAssembler<> assembler(mp,mb);
    gsSparseMatrix<> mass = assembler.assembleMass();
    gsMatrix<> moments = assembler.assembleMoments(g);
    gsSparseSolver<real_t>::LU solver(mass);
    return solver.solve(moments);
}

/*
add Dirichlet values to the solution in the same manner as gsAssembler<T>::constructSolution does
*/
gsMatrix<real_t> addDirVal(const gsAssembler<real_t> &a, const gsMatrix<real_t> &solVector)
{
      gsDofMapper mapper = a.system().colMapper(0); //DoF mapper
      size_t n = a.multiBasis(0).size();
      gsMatrix<real_t> solVector_new(n, 1);
      for (size_t i = 0; i < n; ++i)
      {
            if (mapper.is_free(i, 0)) // DoF value is in the solVector, unknown is 0 since we only have 1.
            {
                  solVector_new.row(i) = solVector.row(mapper.index(i, 0));
            }
            else // eliminated DoF: fill with Dirichlet data
            {
                  solVector_new.row(i) = a.fixedDofs(0).row(mapper.bindex(i, 0)).head(1);
            }
      }
      return solVector_new;
}

/*
delete Dirichlet values from the vector.
*/
gsMatrix<real_t> reduceDirichlet(const gsAssembler<real_t> &a, const gsMatrix<real_t> &w)
{
    gsDofMapper mapper = a.system().colMapper(0);
    const size_t n = a.multiBasis(0).size();
    gsMatrix<real_t> w_new(a.numDofs(), 1);
    index_t k=0;
      for (size_t i = 0; i < n; ++i)
      {
            if (mapper.is_free(i, 0)) // not part of the Dirichlet Boundary
                  w_new.row(i - k) = w.row(i);
        else
            ++k;
      }

      return w_new;
}
