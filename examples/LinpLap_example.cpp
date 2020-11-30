/** @file poisson_example.cpp

@brief Tutorial on how to use G+Smo to solve the Poisson equation,
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

gsMatrix<real_t> projectL2(gsMultiPatch<real_t> mp, gsMultiBasis<real_t> mb, gsFunction<real_t> &g);
gsMatrix<real_t> addDirVal(gsAssembler<real_t> a, gsMatrix<real_t> solVector);
gsMatrix<real_t> reduceDirichlet(gsAssembler<real_t> a, gsMatrix<real_t> w_);

int main(int argc, char* argv[])
{
	real_t eps = 1;		
	real_t p = 1.8;			//p-Laplace parameter
	index_t k = 1;			//spline degree
	index_t maxiter = 100;
	real_t TOL = 1e-12;		//residual error tolerance
	index_t num = 3;		//number of refinements
	index_t str = 1;
	bool require_fin = false;
	gsCmdLine cmd("Linearized p-Laplace example");

	cmd.addReal("e", "eps", "variable for eps", eps);
	cmd.addReal("p", "pow", "p-Laplace Parameter", p);
	cmd.addInt("k", "degree", "degree of basis", k);
	cmd.addInt("i", "maxiter", "maximal iterations", maxiter);
	cmd.addInt("r", "numRefine", "number of refinements of the mesh", num);
	cmd.addInt("s", "strat", "Method for Dirichlet Imposition", str);
	cmd.addSwitch("fin", "After computation, wait until button is pressed", require_fin);
	try { cmd.getValues(argc, argv); } catch (int rv) { return rv; }
	
	gsInfo << "Printing command line arguments:\n\n\n"
		<< "eps = " << eps << "\n\n"
		<< "p = " << p << "\n\n"
		<< "Degree = " << k << "\n\n"
		<< "Maxiter = " << maxiter << "\n\n"
		<< "Number of refinements = " << num << "\n\n";

	gsOptionList opt = gsAssembler<>::defaultOptions();
	//opt.setInt("DirichletValues", dirichlet::l2Projection);

	if (str == 1)
	{
		gsInfo << "DirichletStrategy = dirichlet::elimination\n\n";
		opt.setInt("DirichletStrategy", dirichlet::elimination);
	}
	else if(str == 2)
	{
		gsInfo << "DirichletStrategy = dirichlet::nitsche\n\n";
		opt.setInt("DirichletStrategy", dirichlet::nitsche);
	}

	opt.setReal("quA", 1.0); //number of quadrature points = quA * deg + quB;
	opt.setInt("quB", 1);
	gsInfo << opt << "\n";

	real_t e_0 = 0;
	real_t e_F = 0;
	real_t e_0old;
	real_t e_Fold;

	real_t Lp_rate;
	real_t F_rate;

	//! [Function data]

	double gamma = 10;

	// Define source function
	gsFunctionExpr<> f1("-4*(" + std::to_string(eps*eps) + "+4*x^2+4*y^2)^(" + std::to_string(p) + "/2-1)-8*(" + std::to_string(p) + "-2)*(" + std::to_string(eps*eps) + "+4*x^2+4*y^2)^(" + std::to_string(p) + "/2-2)*(x^2+y^2)", 2);
	//gsFunctionExpr<> f1("-4*(" + std::to_string(eps*eps) + "+4*(x^2+y^2))^((" + std::to_string(p) + "-4)/2)*(" + std::to_string(eps*eps) + "+2*" + std::to_string(p) + "*(x^2+y^2))", 2);
	gsFunctionExpr<> f2("2*" + std::to_string(gamma) + "^2*pi^2*(" + std::to_string(eps*eps) + "+2*" + std::to_string(gamma) + "^2*pi^2*cos(" + std::to_string(gamma) + "*pi*(x+y))^2)^((" + std::to_string(p) + "-4)/2)*(" + std::to_string(eps*eps) + "+2*" + std::to_string(gamma) + "^2*(" + std::to_string(p) + "-1)*pi^2*cos(" + std::to_string(gamma) + "*pi*(x+y))^2)*sin(" + std::to_string(gamma) + "*pi*(x+y))", 2);
	gsFunctionExpr<> f3("8*pi^2*(" + std::to_string(eps*eps) + "+2*pi^2+pi^2*(-(" + std::to_string(p) + "-2)*cos(4*pi*y)-cos(4*pi*x)*(" + std::to_string(p) + "-2+2*(" + std::to_string(p) + "-1)*cos(4*pi*y))))*(" + std::to_string(eps*eps) + "+2*pi^2-pi^2*(cos(4*pi*(x-y))+cos(4*pi*(x+y))))^((" + std::to_string(p) + "-4)/2)*(sin(2*pi*x)*sin(2*pi*y))", 2);
	gsFunctionExpr<> f4("(" + std::to_string(eps*eps) + "+cos(x)^2)^(" + std::to_string(p) + "/2-2)*(" + std::to_string(eps*eps) + "+(" + std::to_string(p) + "-1)*cos(x)^2)*sin(x)", 2);

	// Define exact solution (optional)
	gsFunctionExpr<> u1("x^2+y^2", 2);
	gsFunctionExpr<> u2("sin(" + std::to_string(gamma) + "*pi*(x+y))", 2);
	gsFunctionExpr<> u3("sin(2*pi*x)*sin(2*pi*y)", 2);
	gsFunctionExpr<> u4("sin(x)", 2);

	gsFunctionExpr<> f = f2;
	gsFunctionExpr<> u = u2;
	gsFunctionExpr<> u0("("+u.expression()+")*exp(x*(1-x)*y*(1-y))", 2);

	// Print out source function and solution
	gsInfo  << "Source function " << f << "\n";
	gsInfo  << "Exact solution " << u << "\n";
	gsInfo << "Initial guess" << u0 << "\n\n";
	//! [Function data]

	
	//! [Geometry data]
	// Define Geometry, must be a gsMultiPatch object
	//gsMultiPatch<> patches;

	gsMultiPatch<> patch = gsMultiPatch<>(*gsNurbsCreator<real_t>::BSplineSquareDeg(k));

	//! [Geometry data]

	//! [Boundary conditions]
	gsBoundaryConditions<> bcInfo;
	bcInfo.addCondition(0, boundary::west, condition_type::dirichlet, &u);
	bcInfo.addCondition(0, boundary::east, condition_type::dirichlet, &u);
	bcInfo.addCondition(0, boundary::north, condition_type::dirichlet, &u);
	bcInfo.addCondition(0, boundary::south, condition_type::dirichlet, &u);

	//! [Boundary conditions]

	//! [Refinement]
	// Copy basis from the geometry

	int startrefine = 0;

	gsMultiBasis<> refine_basis(patch);
	for (int i = 0; i < startrefine; i++)
	{
		refine_basis.uniformRefine();
	}

	////////////// Setup solver and solve //////////////
	// Initialize Solver
	// Setup method for handling Dirichlet boundaries, options:
	//
	// * elimination: Eliminate the Dirichlet DoFs from the linear system.
	//
	// * nitsche: Keep the Dirichlet DoFs and enforce the boundary
	//
	// condition weakly by a penalty term.
	// Setup method for handling patch interfaces, options:
	//
	// * glue:Glue patches together by merging DoFs across an interface into one.
	//   This only works for conforming interfaces
	//
	// * dg: Use discontinuous Galerkin-like coupling between adjacent patches.
	//       (This option might not be available yet)

	//not used

	int n = refine_basis.size();

	//starting solution as 0
	//gsMatrix<real_t> w_ = projectL2(patch, refine_basis, g0);
	gsMatrix<real_t> w_ = gsMatrix<real_t>::Zero(n,1);

	gsLinpLapPde<real_t> pde(patch, bcInfo, f, eps, p, w_);

	gsLinpLapAssembler<real_t> assembler; 
	assembler.initialize(pde, refine_basis, opt);

	gsMatrix<real_t> solVector = w_;

	gsSparseMatrix<real_t, 1> transfer;
	gsSparseMatrix<real_t, 1> Kh;
	//gsSparseMatrix<real_t, 1> Mh;

	gsMatrix<real_t> fh;
	gsMatrix<real_t> rh;
	//gsMatrix<real_t> vh

	gsField<> sol;
	//gsField<> solnew;

	//variables for errors

	gsInfo <<"eps = " << eps << " , p = " << p << " , k = " << k  << "\n";
	gsInfo << "mesh size & CPU time & L_p error & L_p rate & F error & F rate & N_max \n";

	for (int i = startrefine; i < num; i++) 
	{
		//transfer recent solution to finer mesh. with elimination it only transfers free DoFs and not Dirichlet values.
		refine_basis.uniformRefine_withTransfer(transfer, bcInfo, opt);

		//refine_basis.uniformRefine();

		/*
		gsInfo << refine_basis.size() << "\n";
		gsInfo << transfer.rows() << " x " << transfer.cols() << "\n";
		gsInfo << pde.w.size() << "\n";
		*/
		
		int n = refine_basis.size();

		//start with 0 solution for every mesh for now, later with transfer.
		pde.w = gsMatrix<real_t>::Zero(n, 1);
		//pde.w = projectL2(patch, refine_basis, u0);
		//pde.w = transfer * pde.w; //update w

		//gsLinpLapAssembler<real_t> assembler;

		assembler.initialize(pde, refine_basis, opt);

		assembler.assemble();
		
		Kh = assembler.matrix();

		fh = assembler.rhs(); //set rhs, which is constant for current mesh --> independant of w

		int iter = 0;

		//gsInfo << "------------------------------------------------" << "\n";

		std::clock_t c_start = std::clock();

		do
		{
			gsSparseSolver<>::LU solver(Kh);
			solVector = solver.solve(fh);

			pde.w = addDirVal(assembler,solVector); //add Dirichlet values to current solution and set as new w.

			//gsInfo << pde.w << "\n";
			//gsInfo << pde.w.size() << "\n";

		//	gsLinpLapAssembler<real_t> assembler;
			assembler.initialize(pde, refine_basis, opt);

			assembler.assemble();

			Kh = assembler.matrix(); //compute new lhs matrix to compute residuum of the nonlinear problem --> Kh(uh)*uh-fh
			fh = assembler.rhs();

			rh = Kh * solVector - fh;

			iter++;
		} while (iter < maxiter && rh.norm()>TOL);

		std::clock_t c_end = std::clock();
		double time = (c_end - c_start) / CLOCKS_PER_SEC;

		
		//transfer coarse solutions to finer meshes for better computation of the norms. neglected right now for simplicity.
		/*
		gsMultiBasis<> norm_bases = refine_basis;
		gsMatrix<> norm_solVector = solVector;

		for (int k = i; k <= 5; k++)
		{
			gsSparseMatrix<real_t, 1> transfer2;
			norm_bases.uniformRefine_withTransfer(transfer2, bcInfo, assembler.options());

			norm_solVector = transfer2 * norm_solVector;
		}

		gsLinpLapAssembler<real_t> norm_assembler(pde, norm_bases,
			dirichlet::nitsche, iFace::glue);

		gsMultiPatch<> mpsol;
		norm_assembler.constructSolution(norm_solVector, mpsol);
		gsField<> sol(norm_assembler.patches(), mpsol);
		*/
		
		gsMultiPatch<> mpsol;
		assembler.constructSolution(solVector, mpsol); //construct solution from the free DoFs via the assembler that is set to elimination.
		gsField<> sol(assembler.patches(), mpsol);

		e_0old = e_0;
		e_Fold = e_F;

		e_0 = sol.distanceLp(u, refine_basis, p, false);
		e_F = sol.distanceF(u, refine_basis, eps, p, false);

		if (i == startrefine)
		{
			gsInfo << 1. / sqrt(refine_basis.totalElements()) << " & " << time << "s & " << e_0 << " & - & " << e_F << " & - & " << iter << "\n";
		}
		else
		{
			Lp_rate = math::log(e_0 / e_0old) / math::log(0.5);
			F_rate = math::log(e_F / e_Fold) / math::log(0.5);
			gsInfo << 1. / sqrt(refine_basis.totalElements()) << " & " << time << "s & " << e_0 << " & " << Lp_rate << " & " << e_F << " & " << F_rate << " & " << iter << "\n";
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
/*
	Project Function g to discrete space with basis mb by solving Mh*uh=fh for uh 
	where fh_i=<g,phi_i> the moment vector and Mh_ij=<phi_j,phi_i> the mass matrix
*/
gsMatrix<real_t> projectL2(gsMultiPatch<real_t> mp, gsMultiBasis<real_t> mb, gsField<real_t> g)
{
	gsSparseSolver<real_t>::LU solver;
	gsAssembler<> MA;

	gsOptionList opt = gsAssembler<>::defaultOptions();
	opt.setInt("DirichletValues", dirichlet::l2Projection);
	opt.setInt("DirichletStrategy", dirichlet::elimination);
	opt.setInt("InterfaceStrategy", iFace::conforming);

	gsBoundaryConditions<real_t> bcInfo;
	
	bcInfo.addCondition(0, boundary::west, condition_type::neumann, 0);
	bcInfo.addCondition(0, boundary::east, condition_type::neumann, 0);
	bcInfo.addCondition(0, boundary::north, condition_type::neumann, 0);
	bcInfo.addCondition(0, boundary::south, condition_type::neumann, 0);
	
	gsPoissonPde<real_t> pde(mp, bcInfo, g.function());

	MA.initialize(pde, mb, opt);

	gsDofMapper mapper; // Gets the indices mapped from Basis --> Matrix

	mb.getMapper((dirichlet::strategy)opt.getInt("DirichletStrategy"),
		(iFace::strategy)opt.getInt("InterfaceStrategy"),
		bcInfo, mapper, 0);

	gsSparseSystem<> sys(mapper);
	sys.reserve(MA.multiBasis(0), MA.options(), MA.pde().numRhs()); // reserving enough space is crutial for performance!
	MA.setSparseSystem(sys);


	MA.push<gsVisitorMass<real_t>>();

	gsSparseMatrix<real_t, 1> Mh = MA.matrix();

	MA.push<gsVisitorMoments<real_t>>(gsVisitorMoments<real_t>(g.function()));

	MA.finalize();

	gsMatrix<real_t> fh = MA.rhs();

	solver.compute(Mh);
	return solver.solve(fh);
}

gsMatrix<real_t> projectL2(gsMultiPatch<real_t> mp, gsMultiBasis<real_t> mb, gsFunction<real_t> &g)
{
	gsSparseSolver<real_t>::LU solver;
	gsAssembler<> MA;

	gsOptionList opt = gsAssembler<>::defaultOptions();
	opt.setInt("DirichletValues", dirichlet::l2Projection);
	opt.setInt("DirichletStrategy", dirichlet::elimination);
	opt.setInt("InterfaceStrategy", iFace::conforming);

	gsBoundaryConditions<real_t> bcInfo;

	bcInfo.addCondition(0, boundary::west, condition_type::neumann, 0);
	bcInfo.addCondition(0, boundary::east, condition_type::neumann, 0);
	bcInfo.addCondition(0, boundary::north, condition_type::neumann, 0);
	bcInfo.addCondition(0, boundary::south, condition_type::neumann, 0);

	gsPoissonPde<real_t> pde(mp, bcInfo, g);

	MA.initialize(pde, mb, opt);

	gsDofMapper mapper; // Gets the indices mapped from Basis --> Matrix

	mb.getMapper((dirichlet::strategy)opt.getInt("DirichletStrategy"),
		(iFace::strategy)opt.getInt("InterfaceStrategy"),
		bcInfo, mapper, 0);

	gsSparseSystem<> sys(mapper);
	sys.reserve(MA.multiBasis(0), MA.options(), MA.pde().numRhs()); // reserving enough space is crutial for performance!
	MA.setSparseSystem(sys);


	MA.push<gsVisitorMass<real_t>>();

	gsSparseMatrix<real_t, 1> Mh = MA.matrix();

	MA.push<gsVisitorMoments<real_t>>(gsVisitorMoments<real_t>(g));

	MA.finalize();

	gsMatrix<real_t> fh = MA.rhs();

	solver.compute(Mh);
	return solver.solve(fh);
}
/*
	add Dirichlet values to the solution in the same manner as gsAssembler<T>::constructSolution does
*/
gsMatrix<real_t> addDirVal(gsAssembler<real_t> a, gsMatrix<real_t> solVector)
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
	delete Dirichlet values from the vector. not needed right now.
*/
gsMatrix<real_t> reduceDirichlet(gsAssembler<real_t> a, gsMatrix<real_t> w_)
{
	gsDofMapper mapper = a.system().colMapper(0);

	size_t n = a.multiBasis(0).size();

	gsMatrix<real_t> w_new(a.numDofs(), 1);

	int k = 0;

	for (size_t i = 0; i < n; ++i)
	{
		if (mapper.is_free(i, 0)) // DoF value is in the solVector
		{
			w_new.row(i-k) = w_.row(mapper.index(i, 0));
			k = 0;
		}
		else // eliminated DoF: fill with Dirichlet data
		{
			k++;
		}
	}

	return w_new;
}






