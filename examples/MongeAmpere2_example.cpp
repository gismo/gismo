/** @file MongeAmpere2_example.cpp

    @brief Tutorial on how to use expression assembler to solve the Monge-Ampere equation

                -det(Hess(u))       = f, <==> -div(grad(u)) + epsilon u= -( div(grad(u))^2 + d!(f- det(Hess(u)) )  )

                            u=g & eps=0 or \partial_n u = g & eps = 1e-8 : eps is the pinalization coefficient
                            u is convex (u in H^2(\Omega))

    N. A simple test of the solver using the nonlinear Poisson equation equation is provided
                        -div((1+u**2)grad(u)) = f

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris M. Bahari
*/

//! [Include namespace]
#include <gismo.h>

using namespace gismo;
//! [Include namespace]

int main(int argc, char *argv[])
{
    //! [Parse command line]
    bool plot = false;
    index_t numRefine  = 4;
    index_t numElevate = 0;
    index_t maxIter = 30;
    double eps{1e-7}; /// pinalization coefficient
    double l2errRes{0.}, tolerancePicard{1e-6};
    bool last{false}, export_b64{false}, adaptiveMesh{false};
    std::string fn("pde/MongeAmpere2d_bvp.xml");

    gsCmdLine cmd("Tutorial on solving a nonlinear Poisson(Monge-Ampere) problem.");
    cmd.addInt( "e", "degreeElevation",
                "Number of degree elevation steps to perform before solving (0: equalize degree in all directions)", numElevate );
    cmd.addInt( "r", "uniformRefine", "Number of Uniform h-refinement loops",  numRefine );
    cmd.addString( "f", "file", "Input XML file", fn );
    cmd.addSwitch("last", "Solve solely for the last level of h-refinement",
                  last);
    cmd.addSwitch(
        "plot", "Create a ParaView visualization file with the solution", plot);
    cmd.addSwitch("binary", "Use B64 encoding for Paraview", export_b64);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }
    //! [Parse command line]

    //! [Read input file]

    gsFileData<> fd(fn);
    gsInfo << "Loaded file "<< fd.lastPath() <<"\n";

    gsMultiPatch<> mp;
    fd.getId(0, mp); // id=0: Multipatch domain

    gsFunctionExpr<> f;
    // id=1: source function for Monge-Ampere equation
    fd.getId(1, f); // id=5: source function for Nonlinear Poisson equation
    gsInfo<<"Source function "<< f << "\n";

    gsBoundaryConditions<> bc;
    fd.getId(2, bc); // id=2: boundary conditions
    bc.setGeoMap(mp);
    gsInfo<<"Boundary conditions:\n"<< bc <<"\n";

    gsFunctionExpr<> ms;
    fd.getId(3, ms); // id=3: reference solution

    gsOptionList Aopt;
    fd.getId(4, Aopt); // id=4: assembler options

    //! [Read input file]

    //gsOptionList Aopt;

    //! [Read input file]

    //! [Refinement]
    gsMultiBasis<> dbasis(mp, true);//true: poly-splines (not NURBS)

    // Elevate and p-refine the basis to order p + numElevate
    // where p is the highest degree in the bases
    dbasis.setDegree( dbasis.maxCwiseDegree() + numElevate);
    dbasis.degreeElevate(1);

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
    auto u_ex = ev.getVariable(ms, G);

    // Solution vector and solution variable
    gsMatrix<> solVector;
    solution u_sol = A.getSolution(u, solVector);

    //! [Problem setup]

    //! [Solver loop]
    gsSparseSolver<>::CGDiagonal solver;

    gsVector<>  h1err(numRefine+1); //l2err(numRefine+1),
    gsInfo<< "(dot1=assembled, dot2=solved, dot3=nonlinear_loop,dot4=got_error)\n"
        "\nDoFs: ";
    double setup_time(0), ma_time(0), slv_time(0), err_time(0);    
    gsStopwatch timer;
    for (int r=0; r<=numRefine; ++r)
    {
        dbasis.uniformRefine();
        mp.uniformRefine();

//        u.setup(bc, dirichlet::interpolation, 0);
//        u.setup(bc, dirichlet::l2Projection, 0);

        // Compute the system matrix and right-hand side

        // Initialize the system : start Computing the conductivity coeffeicient ...
        // Compute the Neumann terms defined on physical space
        auto g_N = A.getBdrFunction(G);
        auto Neumann_Int{ev.integralBdrBc(bc.get("Neumann"), g_N.tr() * nv(G) )};
        // ...
        auto CoeffConductivity{Neumann_Int/ev.integral(pow(2.+2. * ff.val(), 0.5) * meas(G))};
        //... end 

        // Initialize the system
        A.initSystem();
        setup_time += timer.stop();

        gsInfo<< A.numDofs() <<std::flush;

        timer.restart();

        A.assemble(
           igrad(u, G) * igrad(u, G).tr() * meas(G) + eps * u *u.tr()* meas(G) //matrix
           ,
           u* (-1.)*pow(2.+2. * ff.val(), 0.5) *CoeffConductivity* meas(G) //rhs vector
           );
        
        // Compute the Neumann terms defined on physical space
        //auto g_N = A.getBdrFunction(G);
        A.assembleBdr(bc.get("Neumann"), u * g_N.tr() * nv(G) );

        ma_time += timer.stop();

        // gsDebugVar(A.matrix().toDense());
        // gsDebugVar(A.rhs().transpose()   );

        gsInfo<< "." <<std::flush;// Assemblying done

        timer.restart();
        solver.compute( A.matrix() );
        solVector = solver.solve(A.rhs());

        slv_time += timer.stop();

        gsInfo<< "." <<std::flush; // Linear solving done

        // Picard loop
        index_t NiterPicard{0};
        gsMatrix<> sv0; //
        for(int ip{0}; ip<maxIter; ++ip)
        {
            if(adaptiveMesh)
            {
                gsMultiPatch<> UU;
                u_sol.extract(UU);
                gsWrite(UU, "U_solution");
                auto u_s = A.getCoeff(UU);

                //gsMultiBasis<> gbasis(dbasis);
                //gbasis.reduceContinuity(1);
                space v = A.getSpace(dbasis);
                gsMatrix<> vsolVector;
                solution v_sol = A.getSolution(v, vsolVector);
                A.initSystem(2);

                //gsVector<> pt(2); pt.setConstant(0.5);
                //ev.testEval( v, pt );
                //ev.testEval( igrad(u_sol,G), pt );

                // Obtain control points for the gradient of Psi
                A.assemble( v * v.tr() , v * igrad(u_s,G) );
                vsolVector = solver.compute(A.matrix()).solve(A.rhs());
                
                gsMultiPatch<> Psi;
                v_sol.extract(Psi);
                geometryMap PP = A.getMap(Psi);
                auto fp = A.getCoeff(f,PP);

                //::::::::::::::::::::      mesh adaptation solver         :::::::::::::::::::::::::
                sv0 = solVector;
        //        u.setup(bc, dirichlet::interpolation, 0);
        //        u.setup(bc, dirichlet::l2Projection, 0);
                
                solution u_sol = A.getSolution(u, solVector);

                // Initialize the system
                A.initSystem();
                setup_time += timer.stop();

                //gsInfo<< A.numDofs() <<std::flush;

                timer.restart();
                // Compute the system matrix and right-hand side
                // just for test
                // A.assemble(
                //     igrad(u, G) * (1.+pow(u_sol,2))*igrad(u, G).tr() * meas(G) + eps*u*u.tr() * meas(G) //matrix
                //     ,
                //     u * ff* meas(G) //rhs vector
                //     );
                //... Monge-Ampere eqaution

                // .. update Coeffeicient of conductivity
                CoeffConductivity = Neumann_Int/ev.integral(pow( (ilapl(u_sol,G)*ilapl(u_sol,G).tr()).val() + 2.*(fp.val() - ihess(u_sol,G).det()), 0.5) * meas(G));

                A.assemble(
                igrad(u, G) * igrad(u, G).tr() * meas(G) +  eps * u * u.tr() * meas(G) //matrix
                ,
                u * (-1.) * pow( (ilapl(u_sol,G)*ilapl(u_sol,G).tr()).val() + 2.*(fp.val() - ihess(u_sol,G).det()), 0.5) *CoeffConductivity * meas(G) //rhs vector
                );

                // Compute the Neumann terms defined on physical space
                auto g_N = A.getBdrFunction(G);
                A.assembleBdr(bc.get("Neumann"), u * g_N.tr() * nv(G) );

                ma_time += timer.stop();

                // gsDebugVar(A.matrix().toDense());
                // gsDebugVar(A.rhs().transpose() );
                

                gsInfo<< " ." <<std::flush;// Assemblying done

                timer.restart();
                solver.compute( A.matrix() );
                solVector = solver.solve(A.rhs());
                // gsMatrix<> ToUse();
                // solVector = solVector -  ToUse * solVector.sum()/A.numDofs();
                //gsInfo << "Error in Picard algorithm " << math::sqrt( ev.integral( (u_ex - u_sol).sqNorm() * meas(G) ) ) << "\n";

                slv_time += timer.stop();

                gsInfo<< "." <<std::flush; // Linear solving done

                // omp_set_dynamic(0);     // Explicitly disable dynamic teams
                // omp_set_num_threads(1); // Use these threads for later parallel regions

                ++NiterPicard;
                l2errRes = (solVector-sv0).norm();// TODO
                if ( l2errRes < tolerancePicard ) break; // TODO
            }
            else
            {
            sv0 = solVector;
    //        u.setup(bc, dirichlet::interpolation, 0);
    //        u.setup(bc, dirichlet::l2Projection, 0);
            
            solution u_sol = A.getSolution(u, solVector);

            // Initialize the system
            A.initSystem();
            setup_time += timer.stop();

            //gsInfo<< A.numDofs() <<std::flush;

            timer.restart();
            // Compute the system matrix and right-hand side
            // just for test
            // A.assemble(
            //     igrad(u, G) * (1.+pow(u_sol,2))*igrad(u, G).tr() * meas(G) + eps*u*u.tr() * meas(G) //matrix
            //     ,
            //     u * ff* meas(G) //rhs vector
            //     );
            //... Monge-Ampere eqaution

            // .. update Coeffeicient of conductivity
            CoeffConductivity = Neumann_Int/ev.integral(pow( (ilapl(u_sol,G)*ilapl(u_sol,G).tr()).val() + 2.*(ff.val() - ihess(u_sol,G).det()), 0.5) * meas(G));

            A.assemble(
               igrad(u, G) * igrad(u, G).tr() * meas(G) +  eps * u * u.tr() * meas(G) //matrix
               ,
               u * (-1.) * pow( (ilapl(u_sol,G)*ilapl(u_sol,G).tr()).val() + 2.*(ff.val() - ihess(u_sol,G).det()), 0.5) *CoeffConductivity * meas(G) //rhs vector
               );

            // Compute the Neumann terms defined on physical space
            auto g_N = A.getBdrFunction(G);
            A.assembleBdr(bc.get("Neumann"), u * g_N.tr() * nv(G) );

            ma_time += timer.stop();

            // gsDebugVar(A.matrix().toDense());
            // gsDebugVar(A.rhs().transpose()   );
            

            gsInfo<< " ." <<std::flush;// Assemblying done

            timer.restart();
            solver.compute( A.matrix() );
            solVector = solver.solve(A.rhs());
            // gsMatrix<> ToUse();
            // solVector = solVector -  ToUse * solVector.sum()/A.numDofs();
            //gsInfo << "Error in Picard algorithm " << math::sqrt( ev.integral( (u_ex - u_sol).sqNorm() * meas(G) ) ) << "\n";

            slv_time += timer.stop();

            gsInfo<< "." <<std::flush; // Linear solving done

            // omp_set_dynamic(0);     // Explicitly disable dynamic teams
            // omp_set_num_threads(1); // Use these threads for later parallel regions

            ++NiterPicard;
            l2errRes = (solVector-sv0).norm();// TODO
            if ( l2errRes < tolerancePicard ) break; // TODO
            }
        } //for loop
        // ! end Picard loop
        gsInfo<< "\n Niter in Picard : " << NiterPicard << " L2 residual : "<<std::scientific<<l2errRes<<"\n";
        // omp_set_dynamic(0);     // Explicitly disable dynamic teams
        // omp_set_num_threads(1); // Use these threads for later parallel regions

        timer.restart();
        //l2err[r]= math::sqrt( ev.integral( (u_ex - u_sol).sqNorm() * meas(G) ) );

        h1err[r]= math::sqrt(ev.integral( ( igrad(u_ex) - igrad(u_sol,G) ).sqNorm() * meas(G) ));
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
    //gsInfo<< "\nL2 error: "<<std::scientific<<std::setprecision(3)<<l2err.transpose()<<"\n";
    gsInfo<< "H1 error: "<<std::scientific<<h1err.transpose()<<"\n";

    if (!last && numRefine>0)
    {
        // gsInfo<< "\nEoC (L2): " << std::fixed<<std::setprecision(2)
        //       <<  ( l2err.head(numRefine).array()  /
        //            l2err.tail(numRefine).array() ).log().transpose() / std::log(2.0)
        //            <<"\n";

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
        collection.options().setSwitch("base64", export_b64);
        collection.options().setInt("plotElements.resolution", 16);
        collection.newTimeStep(&mp);
        collection.addField(u_sol,"numerical solution");
        collection.addField(igrad(u_sol,G),"gradient_numerical solution");
        if(adaptiveMesh)
        collection.addField(ff, "density function");
        else
        collection.addField(u_ex, "exact solution");
        collection.saveTimeStep();
        collection.save();


        gsFileManager::open("ParaviewOutput/solution.pvd");


        gsMultiPatch<> UU;
        u_sol.extract(UU);
        gsWrite(UU, "U_solution");
        auto u_s = A.getCoeff(UU);

        gsMultiBasis<> gbasis(dbasis);
        gbasis.reduceContinuity(1);
        space v = A.getSpace(gbasis);
        gsMatrix<> vsolVector;
        solution v_sol = A.getSolution(v, vsolVector);
        A.initSystem(2);

        //gsVector<> pt(2); pt.setConstant(0.5);
        //ev.testEval( v, pt );
        //ev.testEval( igrad(u_sol,G), pt );

        // Obtain control points for the gradient of Psi
        A.assemble( v * v.tr() , v * igrad(u_s,G) );
        vsolVector = solver.compute(A.matrix()).solve(A.rhs());
        gsMultiPatch<> Psi;
        v_sol.extract(Psi);
        //geometryMap PP = A.getMap(Psi);
        //auto fp = A.getCoeff(f,PP);

        gsWrite(Psi, "Psi_mapping");
        gsInfo << "Result written in Psi_mapping.xml \n";
    }
    else
        gsInfo << "Done. No output created, re-run with --plot to get a ParaView "
                  "file containing the solution.\n";
    //! [Export visualization in ParaView]

    return EXIT_SUCCESS;


}// end main
