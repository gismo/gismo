/** @file poisson2_example.cpp

    @brief Tutorial on how to use expression assembler to solve a non-linear Poisson equation

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris
*/

//! [Include namespace]
#include <gismo.h>

using namespace gismo;
//! [Include namespace]

int main(int argc, char *argv[])
{
    //! [Parse command line]
    bool plot = false;
    index_t numRefine  = 5;
    index_t numElevate = 0;
    index_t maxIter = 100;
    index_t method = 1;
    bool last = false;

    gsCmdLine cmd("Tutorial on solving a non-linear Poisson problem.");
    cmd.addInt("m","method","Method to use: 0: Newton with automated Jacobian, 1: Newton with Precomputed Jacobian, 2: Picard iteration", method);
    cmd.addInt("i", "iter", "Maximum number of iterations for the iterative method", maxIter);
    cmd.addInt( "e", "degreeElevation",
                "Number of degree elevation steps to perform before solving (0: equalize degree in all directions)", numElevate );
    cmd.addInt( "r", "uniformRefine", "Number of Uniform h-refinement loops",  numRefine );
    //cmd.addString( "f", "file", "Input XML file", fn );
    cmd.addSwitch("last", "Solve solely for the last level of h-refinement", last);
    cmd.addSwitch("plot", "Create a ParaView visualization file with the solution", plot);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }
    //! [Parse command line]

    //! [Read input file]

    //gsFileData<> fd(fn);
    //gsInfo << "Loaded file "<< fd.lastPath() <<"\n";

    gsMultiPatch<> mp = gsNurbsCreator<>::BSplineSquareGrid(1,1,1, .3, .3);

    // Manufactured solition
    gsFunctionExpr<> s("sin(2*pi*y)*sin(2*pi*x)",2);

    // Right-hand side function
    gsFunctionExpr<> f("8*sin(2*pi*y)*pi^2*sin(2*pi*x)*(sin(2*pi*y)^2*sin(2*pi*x)^2-sin(2*pi*y)^2*cos(2*pi*x)^2-sin(2*pi*x)^2*cos(2*pi*y)^2+1)",2);
    gsInfo<<"Source function "<< f << "\n";

    gsFunctionExpr<> bfunc("0",2);
    gsBoundaryConditions<> bc;
    bc.setGeoMap(mp);
    bc.addCondition(0,1, condition_type::dirichlet, &s,0,false);
    bc.addCondition(0,2, condition_type::dirichlet, &s,0,false);
    bc.addCondition(0,3, condition_type::dirichlet, &s,0,false);
    bc.addCondition(0,4, condition_type::dirichlet, &s,0,false);
    //gsDebugVar( bc.allConditions()[0].parametric() );
    gsInfo<<"Boundary conditions:\n"<< bc <<"\n";

    //gsOptionList Aopt;

    //! [Read input file]

    //! [Refinement]
    gsMultiBasis<> dbasis(mp, true);//true: poly-splines (not NURBS)

    // Elevate and p-refine the basis to order p + numElevate
    // where p is the highest degree in the bases
    dbasis.setDegree( dbasis.maxCwiseDegree() + numElevate);

    // h-refine each basis
    if (last)
    {
        for (int r =0; r < numRefine; ++r)
            dbasis.uniformRefine();
        numRefine = 0;
    }

    gsInfo << "Patches: "<< mp.nPatches() <<", degree: "<< dbasis.minCwiseDegree() <<"\n";
#ifdef _OPENMP
    gsInfo<< "Available threads: "<< omp_get_max_threads() <<"\n";
#endif
    //! [Refinement]

    //! [Problem setup]
    gsExprAssembler<> A(1,1);
    //A.setOptions(Aopt);

    gsInfo<<"Active options:\n"<< A.options() <<"\n";

    typedef gsExprAssembler<>::geometryMap geometryMap;
    typedef gsExprAssembler<>::variable    variable;
    typedef gsExprAssembler<>::space       space;
    typedef gsExprAssembler<>::solution    solution;

    // Elements used for numerical integration
    A.setIntegrationElements(dbasis);
    gsExprEvaluator<> ev(A);

    // Set the geometry map
    geometryMap G = A.getMap(mp);

    // Set the discretization space
    space u = A.getSpace(dbasis);

    // Set the source term
    auto ff = A.getCoeff(f, G);

    // Recover manufactured solution
    auto u_ex = ev.getVariable(s, G);

    // Solution vector and solution variable
    gsMatrix<> solVector;
    solution u_sol = A.getSolution(u, solVector);

    //! [Problem setup]

    //! [Solver loop]
    gsSparseSolver<>::CGDiagonal solver;

    gsVector<> l2err(numRefine+1), h1err(numRefine+1);
    gsInfo<< "(dot1=assembled, dot2=solved, dot3=got_error)\n"
        "\nDoFs: ";
    double setup_time(0), ma_time(0), slv_time(0), err_time(0);
    gsStopwatch timer;
    for (int r=0; r<=numRefine; ++r)
    {
        dbasis.uniformRefine();

        u.setup(bc, dirichlet::interpolation, 0); //todo

        // Initialize the system
        A.initSystem();
        setup_time += timer.stop();

        gsInfo<< A.numDofs() <<std::flush;

        A.options().setInt("DirichletStrategy", 1);
        timer.restart();
        // Initial solution (to do: use previous level)
        A.assemble(
            igrad(u, G) * igrad(u, G).tr() * meas(G) //matrix
            ,
            u * ff * meas(G) //rhs vector
            );

        // Compute the Neumann terms defined on physical space
        //auto g_N = A.getBdrFunction(G);
        //A.assembleBdr(bc.get("Neumann"), u * g_N.tr() * nv(G) );

        ma_time += timer.stop();

        gsInfo<< "." <<std::flush;// Assemblying done

        timer.restart();
        solver.compute( A.matrix() );
        solVector = solver.solve( A.rhs() );
        slv_time += timer.stop();

        gsInfo<< "." <<std::flush; // Linear solving done
        int l;

        if (2==method)
        {
            //! Picard method
            gsMatrix<> sv0;
            for (l=0; l<maxIter; ++l)
            {
                sv0 = solVector;
                timer.restart();
                A.clearMatrix();
                A.clearRhs();
                A.assemble( (1+u_sol.val()*
                             u_sol.val())*igrad(u, G) * igrad(u, G).tr() * meas(G),
                            u * ff * meas(G) );
                ma_time += timer.stop();

                timer.restart();
                solver.compute( A.matrix() );
                solVector = solver.solve(A.rhs());
                slv_time += timer.stop();
                if ( (solVector-sv0).norm()<1e-8 ) break;
            }
            //! Picard method
        }
        else
        {
            //! Newton method
            A.options().setInt("DirichletStrategy", 0);// swich off elimination
            auto residual = (1+(u_sol*u_sol).val()) * igrad(u,G) *
                igrad(u_sol, G).tr() * meas(G) - u * ff * meas(G);

            for (l=0; l<maxIter; ++l)
            {
                timer.restart();
                if (1==method) // Precomputed Jacobian
                {
                    A.clearMatrix();
                    A.clearRhs();
                    A.assemble(
                        (1+(u_sol*u_sol).val())*igrad(u,G)*igrad(u,G).tr() * meas(G)
                        + (2*u_sol.val())*igrad(u, G) * igrad(u_sol, G).tr() * u.tr() * meas(G),
                        residual );
                }
                else // Jacobian using finite differences
                {
                    GISMO_ASSERT(0==method, "Invalid method");
                    A.assembleJacobian(residual, u_sol);
                }
                //else // TODO:  Jacobian using automatic differentiation
                ma_time += timer.stop();

                timer.restart();                
                solver.compute( A.matrix() );
                auto du = solver.solve(A.rhs());
                solVector -= du;
                slv_time += timer.stop();
                if ( du.norm() < 1e-8 ) break;
            }
        }
        gsDebug<< "Iterations:"<< l <<" \n";
        //! Newton method

        gsInfo<< "." <<std::flush; // Non-linear iteration done

        timer.restart();
        l2err[r]= math::sqrt( ev.integral( (u_ex - u_sol).sqNorm() * meas(G) ) );
        
        h1err[r]= l2err[r] +
            math::sqrt(ev.integral( ( igrad(u_ex) - igrad(u_sol,G) ).sqNorm() * meas(G) ));
        err_time += timer.stop();
        gsInfo<< ". " <<std::flush; // Error computations done

    } //for loop

    //! [Solver loop]

    timer.stop();
    gsInfo<<"\n\nTotal time: "<< setup_time+ma_time+slv_time+err_time <<"\n";
    gsInfo<<"     Setup: "<< setup_time <<"\n";
    gsInfo<<"  Assembly: "<< ma_time    <<"\n";
    gsInfo<<"   Solving: "<< slv_time   <<"\n";
    gsInfo<<"     Norms: "<< err_time   <<"\n";

    //! [Error and convergence rates]
    gsInfo<< "\nL2 error: "<<std::scientific<<std::setprecision(3)<<l2err.transpose()<<"\n";
    gsInfo<< "H1 error: "<<std::scientific<<h1err.transpose()<<"\n";

    if (!last && numRefine>0)
    {
        gsInfo<< "\nEoC (L2): " << std::fixed<<std::setprecision(2)
              <<  ( l2err.head(numRefine).array()  /
                   l2err.tail(numRefine).array() ).log().transpose() / std::log(2.0)
                   <<"\n";

        gsInfo<<   "EoC (H1): "<< std::fixed<<std::setprecision(2)
              <<( h1err.head(numRefine).array() /
                  h1err.tail(numRefine).array() ).log().transpose() / std::log(2.0) <<"\n";
    }
    //! [Error and convergence rates]

    //! [Export visualization in ParaView]
    if (plot)
    {
        gsInfo<<"Plotting in Paraview...\n";

        gsParaviewCollection collection("ParaviewOutput/solution", &ev);
        collection.options().setSwitch("plotElements", true);
        collection.options().setInt("plotElements.resolution", 16);
        collection.newTimeStep(&mp);
        collection.addField(u_sol,"numerical solution");
        collection.addField(u_ex, "exact solution");
        collection.addField((u_ex-u_sol).norm(), "error");
        collection.saveTimeStep();
        collection.save();


        gsFileManager::open("ParaviewOutput/solution.pvd");
    }
    else
        gsInfo << "Done. No output created, re-run with --plot to get a ParaView "
                  "file containing the solution.\n";
    //! [Export visualization in ParaView]

    return EXIT_SUCCESS;

}// end main
