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

//#define wolfe_powell_debug 1
#define wolf_powell_alt 0

using namespace gismo;
//! [Include namespace]

gsMatrix<real_t> projectL2(gsMultiPatch<real_t> mp, gsMultiBasis<real_t> mb, gsFunction<real_t> &g);
gsMatrix<real_t> addDirVal(gsAssembler<real_t> a, gsMatrix<real_t> solVector);
gsMatrix<real_t> reduceDirichlet(gsAssembler<real_t> a, gsMatrix<real_t> w_);
real_t stepsize(gsSparseMatrix<real_t, 1> &Kh, gsMatrix<real_t> &fh, gsLinpLapPde<real_t> &pde,real_t epsR, gsMultiBasis<real_t> basis, gsOptionList opt, gsMatrix<real_t> u_, gsMatrix<real_t> s_, gsMatrix<real_t> &rh, real_t &Jh, real_t mu, real_t sigma, real_t tau_min, real_t tau_max, index_t subdiv);

int main(int argc, char* argv[])
{
	real_t eps = 1;
	real_t epsR = eps;
	real_t p = 1.8;			//p-Laplace parameter
	index_t k = 1;			//spline degree
	index_t maxiter = 100;
	real_t TOL = 1e-6;		//residual error tolerance
	index_t num = 8;		//number of refinements
	index_t str = 2;
	bool require_fin = true;
	real_t mu = 0.7;
	real_t sigma = 0.3;
	real_t tau_min = 1;
	real_t tau_max = 1;
	index_t subdiv = 1;
	bool prec = false;
	gsCmdLine cmd("Linearized p-Laplace example");

	cmd.addReal("e", "eps", "variable for eps", eps);
	//cmd.addReal("f", "eps_", "variable for eps_", eps_);
	cmd.addReal("p", "pow", "p-Laplace Parameter", p);
	cmd.addInt("k", "degree", "degree of basis", k);
	cmd.addInt("i", "maxiter", "maximal iterations", maxiter);
	cmd.addInt("r", "numRefine", "number of refinements of the mesh", num);
	cmd.addInt("s", "strat", "Method for Dirichlet Imposition", str);
	cmd.addInt("", "sub", "Number of subdivision of an element for quadrature nodes", subdiv);
	cmd.addReal("", "mu", "Wolfe-Powell-Parameter mu", mu);
	cmd.addReal("", "sigma", "Wolfe-Powell-Parameter sigma", sigma);
	cmd.addReal("", "tau_min", "Minimum stepsize", tau_min);
	cmd.addReal("", "tau_max", "Maximum stepsize", tau_max);
	cmd.addReal("t", "tol", "Residual error tolerance", TOL);
	cmd.addSwitch("", "fin", "After computation, wait until button is pressed", require_fin);
	cmd.addSwitch("", "prec", "Preconditioning switch", prec);
	cmd.addReal("", "epsR", "regularizing parameter", epsR);
	try { cmd.getValues(argc, argv); }
	catch (int rv) { return rv; }

	real_t eps_ = eps;
	real_t p_=p;
	real_t ip = 1.8;
	real_t ieps = 1;

	gsInfo << "Printing command line arguments:\n"
		<< "eps               = " << eps << "\n"
		<< "pow               = " << p << "\n"
		<< "degree            = " << k << "\n"
		<< "maxiter           = " << maxiter << "\n"
		<< "numRefine         = " << num << "\n"
		<< "mu                = " << mu << "\n"
		<< "sigma             = " << sigma << "\n"
		<< "tau_min           = " << tau_min << "\n"
		<< "tau_max           = " << tau_min << "\n"
		<< "tol               = " << TOL << "\n"
		<< "subdiv            = " << subdiv << "\n"
		<< "prec              = " << prec << "\n"
		<< "epsR              = " << epsR << "\n";
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

	double gamma = 2;
	int l = 3;
	double lambda = 2.2;//l - 2. / p + 0.01;

						// Define source function
	gsFunctionExpr<> f1("-(2*" + std::to_string(lambda) + "^2*(0.5+(x-1)*x+(y-1)*y)^(" + std::to_string(lambda) + "/2)*(" + std::to_string(eps*eps) + "+" + std::to_string(lambda) + "^2*(0.5+(x-1)*x+(y-1)*y)^(" + std::to_string(lambda) + "-1))^(" + std::to_string(p) + "/2)*(2*" + std::to_string(lambda) + "*(2+" + std::to_string(lambda) + "*(" + std::to_string(p) + "-1)-" + std::to_string(p) + ")*(0.5+(x-1)*x+(y-1)*y)^(" + std::to_string(lambda) + ")+" + std::to_string(eps*eps) + "*(1+2*(x-1)*x+2*(y-1)*y)))/(2*" + std::to_string(lambda) + "^2*(0.5+(x-1)*x+(y-1)*y)^(" + std::to_string(lambda) + ")+" + std::to_string(eps*eps) + "*(1+2*(x-1)*x+2*(y-1)*y))^2", 2);
	//gsFunctionExpr<> f1("-4*(" + std::to_string(eps*eps) + "+4*(x^2+y^2))^((" + std::to_string(p) + "-4)/2)*(" + std::to_string(eps*eps) + "+2*" + std::to_string(p) + "*(x^2+y^2))", 2);
	gsFunctionExpr<> f2("2*" + std::to_string(gamma) + "^2*pi^2*(" + std::to_string(eps*eps) + "+2*" + std::to_string(gamma) + "^2*pi^2*cos(" + std::to_string(gamma) + "*pi*(x+y))^2)^((" + std::to_string(p) + "-4)/2)*(" + std::to_string(eps*eps) + "+2*" + std::to_string(gamma) + "^2*(" + std::to_string(p) + "-1)*pi^2*cos(" + std::to_string(gamma) + "*pi*(x+y))^2)*sin(" + std::to_string(gamma) + "*pi*(x+y))", 2);
	gsFunctionExpr<> f3("8*pi^2*(" + std::to_string(eps*eps) + "+2*pi^2+pi^2*(-(" + std::to_string(p) + "-2)*cos(4*pi*y)-cos(4*pi*x)*(" + std::to_string(p) + "-2+2*(" + std::to_string(p) + "-1)*cos(4*pi*y))))*(" + std::to_string(eps*eps) + "+2*pi^2-pi^2*(cos(4*pi*(x-y))+cos(4*pi*(x+y))))^((" + std::to_string(p) + "-4)/2)*(sin(2*pi*x)*sin(2*pi*y))", 2);
	gsFunctionExpr<> f4("(" + std::to_string(eps*eps) + "+cos(x)^2)^(" + std::to_string(p) + "/2-2)*(" + std::to_string(eps*eps) + "+(" + std::to_string(p) + "-1)*cos(x)^2)*sin(x)", 2);

	gsFunctionExpr<> f5("1", 2);

	// Define exact solution (optional)
	gsFunctionExpr<> u1("((x-0.5)^2+(y-0.5)^2)^(" + std::to_string(lambda) + "/2)", 2);
	gsFunctionExpr<> u2("sin(" + std::to_string(gamma) + "*pi*(x+y))", 2);
	gsFunctionExpr<> u3("sin(2*pi*x)*sin(2*pi*y)", 2);
	gsFunctionExpr<> u4("sin(x)", 2);

	gsFunctionExpr<> u2_derEast("(" + std::to_string(eps*eps) + "+2*" + std::to_string(gamma) + "^2*pi^2*cos(" + std::to_string(gamma) + "*pi*(1+y))^2)^((" + std::to_string(p) + "-2)/2)*(" + std::to_string(gamma) + "*pi*cos(pi*" + std::to_string(gamma) + "*(1+y)))", 2);
	gsFunctionExpr<> u2_derWest("-(" + std::to_string(eps*eps) + "+2*" + std::to_string(gamma) + "^2*pi^2*cos(" + std::to_string(gamma) + "*pi*y)^2)^((" + std::to_string(p) + "-2)/2)*(" + std::to_string(gamma) + "*pi*cos(pi*" + std::to_string(gamma) + "*y))", 2);
	gsFunctionExpr<> u2_derNorth("(" + std::to_string(eps*eps) + "+2*" + std::to_string(gamma) + "^2*pi^2*cos(" + std::to_string(gamma) + "*pi*(x+1))^2)^((" + std::to_string(p) + "-2)/2)*(" + std::to_string(gamma) + "*pi*cos(pi*" + std::to_string(gamma) + "*(x+1)))", 2);
	gsFunctionExpr<> u2_derSouth("-(" + std::to_string(eps*eps) + "+2*" + std::to_string(gamma) + "^2*pi^2*cos(" + std::to_string(gamma) + "*pi*x)^2)^((" + std::to_string(p) + "-2)/2)*(" + std::to_string(gamma) + "*pi*cos(pi*" + std::to_string(gamma) + "*x))", 2);

	gsFunctionExpr<> f = f2;
	gsFunctionExpr<> u = u2;
	gsFunctionExpr<> B("(" + u.expression() + ")*exp(x*(1-x)*y*(1-y))", 2);
	gsFunctionExpr<> Z("0", 2);
	gsFunctionExpr<> L("x+y", 2);

	gsFunctionExpr<> u0 = B;


	// Print out source function and solution
	gsInfo << "Source function " << f << "\n";
	gsInfo << "Exact solution " << u << "\n";
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

	gsBoundaryConditions<> hbcInfo;
	hbcInfo.addCondition(0, boundary::west, condition_type::dirichlet, &Z);
	hbcInfo.addCondition(0, boundary::east, condition_type::dirichlet, &Z);
	hbcInfo.addCondition(0, boundary::north, condition_type::dirichlet, &Z);
	hbcInfo.addCondition(0, boundary::south, condition_type::dirichlet, &Z);

	//! [Boundary conditions]

	//! [Refinement]
	// Copy basis from the geometry

	int startrefine = 1;

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

	//int n = refine_basis.size();

	gsMatrix<real_t> w_ = projectL2(patch, refine_basis, L);
	//gsMatrix<real_t> w_ = gsMatrix<real_t>::Zero(n, 1);

	if (str == 2) { hbcInfo = bcInfo; }

	gsLinpLapPde<real_t> pde(patch, bcInfo, f, eps, p, w_);
	gsLinpLapPde<real_t> pde_(patch, hbcInfo, f, eps_, p, w_);

	gsLinpLapAssembler<real_t> A;
	A.initialize(pde, epsR, refine_basis, opt, subdiv);

	gsLinpLapAssembler<real_t> rA;
	rA.initialize(pde_, epsR, refine_basis, opt, subdiv, prec);

	A.initialize(pde, epsR, refine_basis, opt, subdiv);
	rA.initialize(pde_, epsR, refine_basis, opt, subdiv, prec);

	gsMatrix<real_t> solVector = w_;
	gsMatrix<real_t> step;

	gsSparseMatrix<real_t, 1> transfer;
	gsSparseMatrix<real_t, 1> Kh;
	gsSparseMatrix<real_t, 1> Kh_;
	//gsSparseMatrix<real_t, 1> Mh;

	gsMatrix<real_t> fh;
	gsMatrix<real_t> fh_;
	gsMatrix<real_t> rh;
	gsMatrix<real_t> r0;
	//gsMatrix<real_t> vh

	real_t Jh;

	gsField<> sol;
	//gsField<> solnew;

	gsInfo << "eps = " << eps << " , p = " << p << " , k = " << k << " , lambda = " << lambda << "\n";
	gsInfo << "Dofs & CPU time & L_p error & L_p rate & F error & F rate & N_max \n";

	for (int i = startrefine; i < num; i++)
	{
		if (prec)
		{
			p_ = math::max(p,ip - (i - startrefine)*(ip-p) / (5-startrefine));
			eps_ = math::max(eps,ieps - (i - startrefine)*(ieps - eps) / (5 - startrefine)); //reach the true parameter after 5 iterations
		}

		//transfer recent solution to finer mesh. with elimination it only transfers free DoFs and not Dirichlet values.
		refine_basis.uniformRefine_withTransfer(transfer, bcInfo, opt2);

		//refine_basis.uniformRefine();
		//A.initialize(pde, refine_basis, opt, subdiv);
		//rA.initialize(pde_, refine_basis, opt, subdiv, prec);

		/*
		gsInfo << A.numDofs() << "\n";
		gsInfo << transfer.rows() << " x " << transfer.cols() << "\n";
		gsInfo << pde.w.size() << "\n";
		*/

		//int n = refine_basis.size();

		//start with 0 solution for every mesh for now, later with transfer.
		//pde.w = gsMatrix<real_t>::Zero(n, 1);
		//pde.w = projectL2(patch, refine_basis, u0);

		pde.w = transfer * pde.w; //update w
		pde_.w = pde.w;
		pde_.p = p_;
		pde_.eps = eps_;

		A.initialize(pde_, epsR, refine_basis, opt, subdiv);
		rA.initialize(pde_, epsR, refine_basis, opt, subdiv, prec);

		A.assemble();
		rA.assemble();

		Kh = A.matrix();
		Kh_ = rA.matrix();

		fh = A.rhs();
		//fh_ = rA.rhs();

		//Jh = A.energy();

		/*
		gsInfo << Kh.rows() << " x " << Kh.cols() << "\n";
		gsInfo << solVector.size() << "\n";
		*/

		solVector = reduceDirichlet(A, pde.w);

		rh = Kh * solVector - fh;
		r0 = rh;
		int iter = 0;

		//gsInfo << "------------------------------------------------" << "\n";

		std::clock_t c_start = std::clock();

		do
		{
			gsSparseSolver<>::LU solver(Kh_);
			step = solver.solve(-rh);

			//gsInfo << (rh.transpose()*step).value()/(rh.norm()*step.norm()) << "\n";

			//std::cin.get();

			real_t tau =2*eps_/epsR;// stepsize(Kh, fh, pde, epsR, refine_basis, opt, solVector, step, rh, Jh, mu, sigma, tau_min, tau_max, subdiv);

			solVector = solVector + tau * step;

			pde.w = addDirVal(A, solVector); //add Dirichlet values to current solution and set as new w.
			pde_.w = pde.w;

			A.initialize(pde_, epsR, refine_basis, opt);
			rA.initialize(pde_, epsR, refine_basis, opt, subdiv, prec);

			A.assemble();
			rA.assemble();

			Kh = A.matrix(); //compute new lhs matrix to compute residuum of the nonlinear problem --> Kh(uh)*uh-fh
			Kh_ = rA.matrix();

			fh = A.rhs();
			//fh_ = rA.rhs();

			//Jh = A.energy();

			rh = Kh * solVector - fh;

			//gsInfo << " " << iter << ": rh = " << rh.norm() << ", tau = " << tau << "\n";

			iter++;
		} while (iter < maxiter && (rh.transpose()*rh).value()>TOL*TOL*(r0.transpose()*r0).value());

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

		pde.w = addDirVal(A, solVector);
		pde_.w = pde.w;

		gsMultiPatch<> mpsol;
		A.constructSolution(solVector, mpsol); //construct solution from the free DoFs via the assembler that is set to elimination.
		gsField<> sol(A.patches(), mpsol);

		e_0old = e_0;
		e_Fold = e_F;

		e_0 = sol.distanceLp(u, refine_basis, p, false);
		e_F = sol.distanceF(u, refine_basis, eps, p, false);

		if (i == startrefine)
		{
			gsInfo << A.numDofs() << " & " << time << "s & " << e_0 << " & - & " << e_F << " & - & " << iter << " & " << p_ << " & " << eps_ << "\n";
		}
		else
		{
			Lp_rate = math::log(e_0 / e_0old) / math::log(0.5);
			F_rate = math::log(e_F / e_Fold) / math::log(0.5);
			gsInfo << A.numDofs() << " & " << time << "s & " << e_0 << " & " << Lp_rate << " & " << e_F << " & " << F_rate << " & " << iter << " & " << p_ << " & " << eps_ << "\n";
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
		if (mapper.is_free(i, 0)) // not part of the Dirichlet Boundary
		{
			w_new.row(i - k) = w_.row(i);
		}
		else // eliminated DoF: neglect
		{
			k++;
		}
	}

	return w_new;
}
/*
calculate stepsize of the iteration
*/
real_t stepsize(gsSparseMatrix<real_t, 1> &Kh, gsMatrix<real_t> &fh, gsLinpLapPde<real_t> &pde,real_t epsR, gsMultiBasis<> basis, gsOptionList opt, gsMatrix<real_t> u_, gsMatrix<real_t> s_, gsMatrix<real_t> &rh, real_t &Jh, real_t mu, real_t sigma, real_t tau_min, real_t tau_max, index_t subdiv)
{
	real_t tau = 1;
	//int iter = 0;
	//int maxiter = 10;

	//gsMatrix<real_t> fh;
	gsMatrix<real_t> rh_new;
	//gsSparseMatrix<real_t, 1> Kh;
	real_t Jh_new;

	const real_t downscale = 0.8;
	const real_t upscale = 1.2;

	gsLinpLapAssembler<real_t> A;
	A.initialize(pde, epsR, basis, opt, subdiv);
	A.assemble();

	pde.w = addDirVal(A, u_ + tau * s_);

	A.initialize(pde,epsR, basis, opt, subdiv);
	A.assemble();

	Kh = A.matrix();
	fh = A.rhs();
	rh_new = Kh * (u_ + tau * s_) - fh;
	Jh_new = A.energy();

#ifndef wolf_powell_alt
	bool c1 = Jh_new <= Jh + tau * mu*(rh.transpose()*s_).value();
#else
	bool c1 = (rh_new.transpose() * rh_new).value() <= (rh.transpose() * rh).value() + mu * tau * (rh.transpose()*Kh*s_).value();
#endif
	bool c2 = (rh_new.transpose()*s_).value() >= sigma * (rh.transpose()*s_).value();
#ifdef wolfe_powell_debug
#ifndef wolf_powell_alt
	gsInfo << "   WP1: " << Jh_new << " <= " << Jh + tau * mu*(rh.transpose()*s_).value()
		<< " = " << Jh << "+" << tau << "*" << mu << "*" << rh.norm() << "*" << s_.norm() << "*" <<
		((rh.transpose()*s_).value() / (rh.norm()*s_.norm())) << "     " << (c1 ? "true" : "false") << "\n";
#else
	gsInfo << "   WP1: " << (rh_new.transpose() * rh_new).value() << " <= "
		<< (rh.transpose() * rh).value() + mu * tau * (rh.transpose()*Kh*s_).value()
		<< " = " << (rh.transpose() * rh).value() << "+" << mu << "*" << tau
		<< "*" << (rh.transpose()*Kh*s_).value() << "     " << (c1 ? "true" : "false") << "\n";
#endif
	gsInfo << "   WP2: " << (rh_new.transpose()*s_).value() << " >= " << sigma * (rh.transpose()*s_).value()
		<< "     " << (c2 ? "true" : "false") << "\n";
#endif
	if (!c1 && tau*downscale >= tau_min)
	{
		do
		{
			tau *= downscale;

			pde.w = addDirVal(A, u_ + tau * s_);
			A.initialize(pde, epsR, basis, opt, subdiv);
			A.assemble();

			Kh = A.matrix();
			fh = A.rhs();
			rh_new = Kh * (u_ + tau * s_) - fh;
			Jh_new = A.energy();

#ifndef wolf_powell_alt
			c1 = Jh_new <= Jh + tau * mu*(rh.transpose()*s_).value();
#else
			c1 = (rh_new.transpose() * rh_new).value() <= (rh.transpose() * rh).value() + mu * tau * (rh.transpose()*Kh*s_).value();
#endif
#ifdef wolfe_powell_debug
			c2 = (rh_new.transpose()*s_).value() >= sigma * (rh.transpose()*s_).value();
#ifndef wolf_powell_alt
			gsInfo << "   WP1: " << Jh_new << " <= " << Jh + tau * mu*(rh.transpose()*s_).value()
				<< " = " << Jh << "+" << tau << "*" << mu << "*" << rh.norm() << "*" << s_.norm() << "*" <<
				((rh.transpose()*s_).value() / (rh.norm()*s_.norm())) << "     " << (c1 ? "true" : "false") << "\n";
#else
			gsInfo << "   WP1: " << (rh_new.transpose() * rh_new).value() << " <= "
				<< (rh.transpose() * rh).value() + mu * tau * (rh.transpose()*Kh*s_).value()
				<< " = " << (rh.transpose() * rh).value() << "+" << mu << "*" << tau
				<< "*" << (rh.transpose()*Kh*s_).value() << "     " << (c1 ? "true" : "false") << "\n";
#endif
			gsInfo << "   WP2: " << (rh_new.transpose()*s_).value() << " >= " << sigma * (rh.transpose()*s_).value()
				<< "     " << (c2 ? "true" : "false") << "\n";
#endif
			//iter++;
		} while (!c1 && tau*downscale >= tau_min);
	}
	else if (c1 && !c2 && tau*upscale <= tau_max)
	{
		do
		{
			tau *= upscale;

			pde.w = addDirVal(A, u_ + tau * s_);
			A.initialize(pde, epsR, basis, opt, subdiv);
			A.assemble();

			Kh = A.matrix();
			fh = A.rhs();
			rh_new = Kh * (u_ + tau * s_) - fh;
			Jh_new = A.energy();

#ifndef wolf_powell_alt
			c1 = Jh_new <= Jh + tau * mu*(rh.transpose()*s_).value();
#else
			c1 = (rh_new.transpose() * rh_new).value() <= (rh.transpose() * rh).value() + mu * tau * (rh.transpose()*Kh*s_).value();
#endif
			c2 = (rh_new.transpose()*s_).value() >= sigma * (rh.transpose()*s_).value();
#ifdef wolfe_powell_debug
#ifndef wolf_powell_alt
			gsInfo << "   WP1: " << Jh_new << " <= " << Jh + tau * mu*(rh.transpose()*s_).value()
				<< " = " << Jh << "+" << tau << "*" << mu << "*" << rh.norm() << "*" << s_.norm() << "*" <<
				((rh.transpose()*s_).value() / (rh.norm()*s_.norm())) << "     " << (c1 ? "true" : "false") << "\n";
#else
			gsInfo << "   WP1: " << (rh_new.transpose() * rh_new).value() << " <= "
				<< (rh.transpose() * rh).value() + mu * tau * (rh.transpose()*Kh*s_).value()
				<< " = " << (rh.transpose() * rh).value() << "+" << mu << "*" << tau
				<< "*" << (rh.transpose()*Kh*s_).value() << "     " << (c1 ? "true" : "false") << "\n";
#endif
			gsInfo << "   WP2: " << (rh_new.transpose()*s_).value() << " >= " << sigma * (rh.transpose()*s_).value()
				<< "     " << (c2 ? "true" : "false") << "\n";
#endif
			//iter++;
		} while (c1 && !c2 && tau*upscale <= tau_max);
	}
	rh = rh_new;
	Jh = Jh_new;

	return tau;
}
