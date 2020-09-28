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
#include <gismo.h>
#include <ctime>
#include <Eigen/SparseCholesky>
//#include <unsupported/Eigen/MatrixFunctions>
#include <Eigen/Sparse>



using namespace gismo;
using namespace expr;
//! [Include namespace]

real_t distanceF(gsField<real_t> &f1, gsFunctionSet<real_t> &f2, gsMultiBasis<real_t> mb, real_t eps, real_t p, bool isFunc_param);
gsMatrix<real_t> projectL2(gsMultiPatch<real_t> mp, gsMultiBasis<real_t> mb, gsFunctionExpr<real_t> g);

int main(int argc, char *argv[])
{
	real_t eps = 1;
	real_t p = 1.8;
	int k = 2;
	int maxiter = 30;
	real_t TOL = 1e-9;

	gsMultiPatch<> patch_ = gsMultiPatch<>(*gsNurbsCreator<real_t>::BSplineSquareDeg(k));
	gsMultiBasis<> basis_(patch_);

	for (int i = 0; i < 5; i++)
	{
		basis_.uniformRefine();
	}

	gsFunctionExpr<> f3("(2-" + std::to_string(p) + ")*4*x*(" + std::to_string(eps*eps) + "+4*x^2)^((" + std::to_string(p) + "-4)/2)", 2);
	gsFunctionExpr<> w("x^2", 2);
	gsFunctionExpr<> u("x+y", 2);

	gsMatrix<real_t> _w = projectL2(patch_, basis_, w);

	gsBoundaryConditions<> bcInfo_;
	bcInfo_.addCondition(0, boundary::west, condition_type::dirichlet, &u);
	bcInfo_.addCondition(0, boundary::east, condition_type::dirichlet, &u);
	bcInfo_.addCondition(0, boundary::north, condition_type::dirichlet, &u);
	bcInfo_.addCondition(0, boundary::south, condition_type::dirichlet, &u);
	
	gsLinpLapPde<real_t> pde_(patch_, bcInfo_, f3, eps, p, _w);

	gsLinpLapAssembler<real_t> assembler_(pde_, basis_,
		dirichlet::elimination, iFace::glue);

	gsLinpLapAssembler<real_t> assembler2_(pde_, basis_,
		dirichlet::nitsche, iFace::glue);

	assembler_.assemble();

	gsSparseSolver<real_t>::LU solver_;
	solver_.compute(assembler_.matrix());
	gsMatrix<real_t> solVector_ = solver_.solve(assembler_.rhs());

	gsMultiPatch<> mpsol_;
	assembler_.constructSolution(solVector_, mpsol_);

	gsMultiPatch<> mpsol2_;
	assembler2_.constructSolution(_w, mpsol2_);

	gsField<> sol_(assembler_.patches(), mpsol_);
	gsField<> approx(assembler2_.patches(), mpsol2_);

	gsInfo << sol_.distanceL2(u) << "\n";

	gsMatrix<real_t> points(2,3);
	points << 0.5, 0.7, 1,
			  0.5, 0.33,0.67;
			
	gsInfo << approx.function(0).eval(points) << "\n";
	gsInfo << w.eval(points) << "\n";

	//! [Function data]
	// Define source function
	gsFunctionExpr<> f1("-4*(" + std::to_string(eps*eps) + "+4*x^2+4*y^2)^(" + std::to_string(p) + "/2-1)-8*(" + std::to_string(p) + "-2)*(" + std::to_string(eps*eps) + "+4*x^2+4*y^2)^(" + std::to_string(p) + "/2-2)*(x^2+y^2)", 2);
	gsFunctionExpr<> f2("(" + std::to_string(eps*eps) + "+200*pi^2*cos(10*pi*(x+y))^2)^((" + std::to_string(p) + "-2)/2)*sin(10*pi*(x+y))*200*pi^2+(" + std::to_string(p) + "-2)*(" + std::to_string(eps*eps) + "+200*pi^2*cos(10*pi*(x+y))^2)^((" + std::to_string(p) + "-4)/2)*40000*pi^4*cos(10*pi*(x+y))^2*sin(10*pi*(x+y))", 2);
	// Define exact solution (optional)
	gsFunctionExpr<> g1("x^2+y^2", 2);
	gsFunctionExpr<> g2("sin(10*pi*(x+y))", 2);

	gsFunctionExpr<> f = f2;
	gsFunctionExpr<> g = g2;

	// Print out source function and solution
	gsInfo  << "Source function " << f << "\n";
	gsInfo  << "Exact solution " << g << "\n\n";
	//! [Function data]

	
	//! [Geometry data]
	// Define Geometry, must be a gsMultiPatch object
	//gsMultiPatch<> patches;

	gsMultiPatch<> patch = gsMultiPatch<>(*gsNurbsCreator<real_t>::BSplineSquareDeg(k));

	//! [Geometry data]

	//! [Boundary conditions]
	gsBoundaryConditions<> bcInfo;
	bcInfo.addCondition(0, boundary::west, condition_type::dirichlet, &g);
	bcInfo.addCondition(0, boundary::east, condition_type::dirichlet, &g);
	bcInfo.addCondition(0, boundary::north, condition_type::dirichlet, &g);
	bcInfo.addCondition(0, boundary::south, condition_type::dirichlet, &g);

	//! [Boundary conditions]

	//! [Refinement]
	// Copy basis from the geometry
	gsMultiBasis<> refine_basis(patch);

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

	gsAssembler<real_t> MA;

	gsOptionList opt = gsAssembler<>::defaultOptions();
	opt.setInt("DirichletValues", dirichlet::l2Projection);
	opt.setInt("DirichletStrategy", dirichlet::elimination);
	opt.setInt("InterfaceStrategy", iFace::conforming);

	int n = refine_basis.size();
	gsMatrix<real_t> w_ = gsMatrix<real_t>::Zero(n, 1);

	gsLinpLapPde<real_t> pde(patch, bcInfo, f, eps, p, w_);

	gsLinpLapAssembler<real_t> assembler(pde, refine_basis,
		dirichlet::nitsche, iFace::glue);

	gsMatrix<real_t> solVector;
	gsSparseMatrix<real_t, 1> Kh;
	//gsSparseMatrix<real_t, 1> Mh;

	gsMatrix<real_t> fh;
	gsMatrix<real_t> rh;
	//gsMatrix<real_t> vh;

	gsSparseSolver<real_t>::LU solver;

	gsField<> sol;
	//gsField<> solnew;

	real_t e_0 = 0;
	real_t e_F = 0;
	real_t e_0old;
	real_t e_Fold;

	real_t L2_rate;
	real_t F_rate;

	gsSparseMatrix<real_t, 1> transfer;

	//real_t step;

	gsInfo <<"eps = " << eps << " , p = " << p << " , k = " << k  << "\n";

	for (int i = 1; i <= 6; i++) 
	{
		gsSparseMatrix<real_t, 1> transfer;
		refine_basis.uniformRefine_withTransfer(transfer, bcInfo, assembler.options());
		//patches.uniformRefine();

		pde = gsLinpLapPde<real_t>(patch, bcInfo, f, eps, p, transfer*pde.w);

		assembler = gsLinpLapAssembler<real_t>(pde, refine_basis,
			dirichlet::nitsche, iFace::glue);

		assembler.assemble();

		Kh = assembler.matrix();

		fh = assembler.rhs();
		/*
		MA.initialize(pde, refine_basis, opt);

		gsDofMapper mapper; // Gets the indices mapped from Basis --> Matrix

		refine_basis.getMapper((dirichlet::strategy)opt.getInt("DirichletStrategy"),
			(iFace::strategy)opt.getInt("InterfaceStrategy"),
			bcInfo, mapper, 0);

		gsSparseSystem<> sys(mapper);
		sys.reserve(MA.multiBasis(0), MA.options(), MA.pde().numRhs()); // reserving enough space is crutial for performance!
		MA.setSparseSystem(sys);

		//MA.computeDirichletDofs();

		MA.push<gsVisitorMass<real_t> >();

		//MA.push<gsVisitorNitsche<real_t> >(MA.pde().bc().dirichletSides());

		MA.finalize();
		gsSparseMatrix<real_t> Mh = MA.matrix();
		
		gsInfo << Mh << "\n";
		*/

		/*
		Eigen::SimplicialLLT<gsSparseMatrix<real_t,1>> chol(Mh); // compute the Cholesky decomposition of Mh
		gsSparseMatrix<real_t,1> L = chol.matrixL(); // retrieve factor L  in the decomposition
		*/

		//gsInfo << L.toDense() << "\n";
		//gsInfo << (Mh - L * L).norm << "\n";

		int iter = 0;

		//gsInfo << "------------------------------------------------" << "\n";

		std::clock_t c_start = std::clock();

		do
		{
			/*
			assembler = gsLinpLapAssembler<real_t>(pde, refine_basis,
				dirichlet::nitsche, iFace::glue);

			assembler.assemble();
			*/

			solver.compute(Kh);

			solVector = solver.solve(fh);

			pde.w = solVector;

			assembler = gsLinpLapAssembler<real_t>(pde, refine_basis,
				dirichlet::nitsche, iFace::glue);

			assembler.assemble();

			Kh = assembler.matrix();

			rh = Kh * solVector - fh;

			/*
			solver.compute(Mh);
			vh = solver.solve(rh);
			*/

			iter++;
		} while (iter < maxiter && rh.norm()>TOL);

		std::clock_t c_end = std::clock();
		double time = 1000.0*(c_end - c_start) / CLOCKS_PER_SEC;

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

		e_0old = e_0;
		e_Fold = e_F;

		e_0 = sol.distanceLp(g, norm_bases, p, false);
		e_F = distanceF(sol, g, norm_bases, eps, p, false);

		if (i == 1)
		{
			gsInfo << 1. / sqrt(refine_basis.totalElements()) << " & " << time << "ms & " << e_0 << " & - & " << e_F << " & - & " << iter << "\n";
		}
		else
		{
			L2_rate = math::log(e_0 / e_0old) / math::log(0.5);
			F_rate = math::log(e_F / e_Fold) / math::log(0.5);
			gsInfo << 1. / sqrt(refine_basis.totalElements()) << " & " << time << "ms & " << e_0 << " & " << L2_rate << " & " << e_F << " & " << F_rate << " & " << iter << "\n";
		}
	}

	gsInfo << "fin";
	std::cin.get();
	return EXIT_SUCCESS;
	
}// end main

real_t distanceF(gsField<real_t> &f1, gsFunctionSet<real_t> &f2,gsMultiBasis<real_t> mb, real_t eps, real_t p, bool isFunc_param)
{
	gsExprEvaluator<> ev;
	ev.setIntegrationElements(mb);

	typename gsExprEvaluator<real_t>::geometryMap G = ev.getMap(f1.patches());
	typename gsExprEvaluator<real_t>::variable f_1 =
		(f1.isParametric() ? ev.getVariable(f1.fields()) : ev.getVariable(f1.fields(), G));
	typename gsExprEvaluator<real_t>::variable f_2 = 
		(isFunc_param ? ev.getVariable(f2) : ev.getVariable(f2, G));

	//value_expr<_expr<real_t>> eps_expr(_expr<real_t>(eps*eps));

	if (f1.isParametric() && isFunc_param)
		return math::sqrt(ev.integral((((eps*eps+igrad(f_1,G).sqNorm())^((p-2)/4)) *igrad(f_1,G) - ((eps*eps+igrad(f_2,G).sqNorm())^((p-2)/4)) *igrad(f_2,G)).sqNorm()*meas(G)));
	if (f1.isParametric())
		return math::sqrt(ev.integral((((eps*eps+igrad(f_1,G).sqNorm())^((p-2)/4)) *igrad(f_1,G) - ((eps*eps+igrad(f_2).sqNorm())^((p-2)/4)) *igrad(f_2)).sqNorm()*meas(G)));
	if (isFunc_param)
		return math::sqrt(ev.integral((((eps*eps+igrad(f_1).sqNorm())^((p-2)/4)) *igrad(f_1) - ((eps*eps+igrad(f_2,G).sqNorm())^((p-2)/4)) *igrad(f_2,G)).sqNorm()*meas(G)));

	return math::sqrt(ev.integral((((eps*eps+igrad(f_1).sqNorm())^((p-2)/4)) *igrad(f_1) - ((eps*eps+(igrad(f_2).sqNorm()))^((p-2)/4)) *igrad(f_2)).sqNorm()*meas(G)));
}

gsMatrix<real_t> projectL2(gsMultiPatch<real_t> mp, gsMultiBasis<real_t> mb, gsFunctionExpr<real_t> g)
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

	//MA.computeDirichletDofs();

	MA.push<gsVisitorMass<real_t>>();

	gsSparseMatrix<real_t, 1> Mh = MA.matrix();

	MA.push<gsVisitorMoments<real_t>>(gsVisitorMoments<real_t>(g));

	//MA.push<gsVisitorNitsche<real_t> >(MA.pde().bc().dirichletSides());

	MA.finalize();

	gsMatrix<real_t> fh = MA.rhs();

	solver.compute(Mh);
	return solver.solve(fh);
}




