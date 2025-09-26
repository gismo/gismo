/** @file Monge_Ampere_example.cpp

    @brief Tutorial on how to use expression assembler to builde multipatch-solver from anayltic density function

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): M. BAHARI
*/

//! [Include namespace]
#include <gismo.h>

using namespace gismo;
//! [Include namespace]

int main(int argc, char *argv[])
{
    //! [Parse command line]
    bool plot          = false;
    index_t numRefine  = 3;
    index_t numElevate = 0;
    index_t maxIter    = 50;
    bool last          = false;
    // Specify the file path
    //std::string fn("pde/circle.xml");
    //std::string fn("pde/mhd.xml");
    std::string fn("surfaces/cylinder.xml"); 
    // load the file
    gsFileData<> fd(fn);


    gsCmdLine cmd("Tutorial on solving a non-linear Monge-Ampere problem.");
    cmd.addInt("i", "iter", "Maximum number of iterations for the iterative Picard", maxIter);
    cmd.addInt( "e", "degreeElevation",
                "Number of degree elevation steps to perform before solving (0: equalize degree in all directions)", numElevate );
    cmd.addInt( "r", "uniformRefine", "Number of Uniform h-refinement loops",  numRefine );
    //cmd.addString( "f", "file", "Input XML file", fn );
    cmd.addInt("quRule",
                 "Quadrature rule [1:GaussLegendre,2:GaussLobatto,3:PatchRule]",
                 2);
    cmd.addSwitch("last", "Solve solely for the last level of h-refinement", last);
    cmd.addSwitch("plot", "Create a ParaView visualization file with the solution", plot);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    //gsFileData<> fd(fn);
    //gsInfo << "Loaded file "<< fd.lastPath() <<"\n";
    // .... one single patch
   gsMultiPatch<> mp = gsNurbsCreator<>::BSplineSquareGrid(1,1,1, 0.0, 0.0);
   //... patch 1
   //mp.addPatch(gsNurbsCreator<>::BSplineSquare(1, 0.,1.0));
//    //... patch 2 (L-shape)
    // mp.addPatch(gsNurbsCreator<>::BSplineSquare(1, 1.,0.0));
   //mp.addInterface(0,2,2,1);
   
//    // ... patch 0-1
//    gsMultiPatch<> mp = gsNurbsCreator<>::BSplineSquareGrid(1,2,1, -1.0, -1.0);
//    // ... patch 2
//    mp.addPatch(gsNurbsCreator<>::BSplineSquare(1,-1.,-2.0));
//    mp.addInterface(0,3,2,4);
//    // ... patch 3
//    mp.addPatch(gsNurbsCreator<>::BSplineSquare(1,0.,-2.0));
//    mp.addInterface(2,2,3,1);
//    // ... patch 4
//    mp.addPatch(gsNurbsCreator<>::BSplineSquare(1, 1.,-2.0));
//    mp.addInterface(3,2,4,1);
//    // ... patch 5
//    mp.addPatch(gsNurbsCreator<>::BSplineSquare(1, 1.,-1.0));
//    mp.addInterface(4,4,5,3);
//    // ... patch 6
//    mp.addPatch(gsNurbsCreator<>::BSplineSquare(1, 1., 0.0));
//    mp.addInterface(5,4,6,3);
//    // ... patch 7
//    mp.addPatch(gsNurbsCreator<>::BSplineSquare(1, 0.0, 0.0));
//    mp.addInterface(6,1,7,2);
//    mp.addInterface(1,2,7,1);
   // Get all interfaces and boundaries:
    mp.computeTopology();

    //identity mapping in each direction
    gsFunctionExpr<> sx("x",2);
    gsFunctionExpr<> sy("y",2);
    // Manufactured identity mapping
    gsFunctionExpr<> sN("x","y",2);
    // Right-hand side function : Analytical density function (det(H(u))=f= sigma/rho)
    //gsFunctionExpr<> f("1./(2.+cos(8.*pi*sqrt((x-0.5-0.25*0.)**2+(y-0.5)**2)))",2);
    //gsFunctionExpr<> f("1.+6.*( 1/(1.+exp((y -x  - 0.3)/0.01)) - 1/(1.+exp((y - x  - 0.1)/0.01)) )",2);
    //Manufactured density function
    //gsFunctionExpr<> f("1.+6.*( 1/(1.+exp(( (y-0.98)**2+(x-0.899)**2  - 0.002)/0.001)) + 1/(1.+exp((y -x  - 0.25)/0.001)) - 1./(1.+exp((y - x  - 0.15)/0.001)) +  0/(1.+exp((y - 1.0)/0.001)) - 0./(1.+exp((y - 0.975)/0.001))  +  1/(1.+exp((x - 1.0)/0.001)) - 1./(1.+exp((x - 0.95)/0.001)) )",2);    
    gsFunctionExpr<> f("(1.+ 9./(1.+(10.*sqrt((x-0.7-0.25*0.)**2+(y-0.5)**2)*cos(atan2(y-0.5,x-0.7-0.25*0.) -20.*((x-0.7-0.25*0.)**2+(y-0.5)**2)))**2) )",2);
    //gsFunctionExpr<> f("( 1.+ 5.*exp(-50.*abs((x-0.5-0.25*cos(2.*pi*0.25))**2-(y-0.5-0.5 *sin(2.*pi*0.25))**2- 0.01)))",2);
    //gsFunctionExpr<> f("(1. + 5./cosh( 5.*((x-sqrt(3)/2)**2+(y-0.5)**2 - (pi/2)**2) )**2 + 5./cosh( 5.*((x+sqrt(3)/2)**2+(y-0.5)**2 - (pi/2)**2) )**2)",2);
    //gsFunctionExpr<> f("(1.+10.(1.+0.*x*(abs(y-0.5)<0.3))/cosh( 10.*(x -0.5 ) )**2)",2);
    gsInfo<<"Source function "<< f << "\n";

    gsInfo<<"The domain is "<< mp.detail() << "\n";
    
    gsMultiPatch<> mpLeft;// Initial geometry
    fd.getId(1,mpLeft);

    gsBoundaryConditions<> bc_x;
    bc_x.setGeoMap(mp);
   for ( gsMultiPatch<>::const_biterator
            bit = mp.bBegin(); bit != mp.bEnd(); ++bit)
   {
       bc_x.addCondition( *bit, condition_type::dirichlet, &sx );
   }
    gsBoundaryConditions<> bc_y;
    bc_y.setGeoMap(mp);
   for ( gsMultiPatch<>::const_biterator
            bit = mp.bBegin(); bit != mp.bEnd(); ++bit)
   {
       bc_y.addCondition( *bit, condition_type::dirichlet, &sy );
   }

    //! [Refinement]
    gsMultiBasis<> dbasis(mp, true);//true: poly-splines (not NURBS)

    // Elevate and p-refine the basis to order p + numElevate
    // where p is the highest degree in the bases
    dbasis.setDegree( dbasis.maxCwiseDegree() + numElevate);
    dbasis.degreeElevate(1);

    gsInfo << "Patches: "<< mp.nPatches() <<", degree: "<< dbasis.minCwiseDegree() <<"\n";
#ifdef _OPENMP
    gsInfo<< "Available threads: "<< omp_get_max_threads() <<"\n";
#endif
    //! [Refinement]

    //! [Problem setup]
    gsExprAssembler<> A(1,1);
    //A.setOptions(Aopt);
    // A.options().addInt("quRule",
    //              "Quadrature rule [1:GaussLegendre,2:GaussLobatto,3:PatchRule]",
    //              1);
    // A.options().addInt("InterfaceStrategy", "Interface strategy conforming", iFace::none);
    // A.options().setSwitch("SameElement",false);
    gsInfo<<"Active options:\n"<< A.options() <<"\n";

    typedef gsExprAssembler<>::geometryMap geometryMap;
    typedef gsExprAssembler<>::variable    variable;
    typedef gsExprAssembler<>::space       space;
    typedef gsExprAssembler<>::solution    solution;

    // Elements used for numerical integration
    A.setIntegrationElements(dbasis);
    gsExprEvaluator<> ev(A);
    // ev.options().addInt("quRule",
    //              "Quadrature rule [1:GaussLegendre,2:GaussLobatto,3:PatchRule]",
    //             1);
    // ev.options().addInt("InterfaceStrategy", "Interface strategy conforming", iFace::none);
    // A.options().setSwitch("SameElement",false); 
    // Set the geometry map
    geometryMap G = A.getMap(mp);

    // Set the discretization space
    space u = A.getSpace(dbasis);
    
    // Set the source term
    auto ff = A.getCoeff(f, G);

    // idensity mapping
    auto u_I = ev.getVariable(sN, G);
    auto u_x = ev.getVariable(sx, G);
    auto u_y = ev.getVariable(sy, G);

    // Solution vector and solution variable in x direction
    gsMatrix<> solVectorx;
    solution u_solx = A.getSolution(u, solVectorx);
    // Solution vector and solution variable in y direction
    gsMatrix<> solVectory;
    solution u_soly = A.getSolution(u, solVectory);

    // ... Adaptive mapping
    gsMultiPatch<> Psix;
    gsMultiPatch<> Psiy;
    gsMultiPatch<> Psi;

    // geometry adaptive map
    geometryMap PP = A.getMap(Psi);

    auto ff_Psi = A.getCoeff(f, G);

    //! [Solver loop]
    gsSparseSolver<>::CGDiagonal solver;
    // h-refine each basis
    if (last)
    {
        for (int r =0; r < numRefine; ++r){
            dbasis.uniformRefine();
            //mpLeft.uniformRefine();
        }
        numRefine = 0;
    }
    //... nromalisation of density function
    auto CoeffDensity = 10.;
    gsVector<>  h1err(numRefine+1); //l2err(numRefine+1) : The solution exists up to an additive constant.
    gsInfo<< "(dot1=assembled, dot2=solved, dot3=nonlinear_loop,dot4=got_error)\n"
        "\nDoFs: ";
    double setup_time(0), ma_time(0), slv_time(0), err_time(0);    
    gsStopwatch timer;

    //::::::::::::::::::::      mesh adaptation solver         :::::::::::::::::::::::::
    for (int r=0; r<=numRefine; ++r)
    {
        dbasis.uniformRefine();
        mp.uniformRefine();

        u.setup(bc_x, dirichlet::interpolation, 0);
        // Initialize the system :  identity mapping as initial guess
        A.initSystem();
        setup_time += timer.stop();

        gsInfo<< A.numDofs() <<std::flush;
        timer.restart();
        A.assemble(
        igrad(u, G) * igrad(u, G).tr() *meas(G) //matrix
        ,
         (1.- ff.val()/CoeffDensity) * igrad(u, G) * grad(u_x).tr() *meas(G) //rhs vector
        );

        ma_time += timer.stop();

        // gsDebugVar(A.matrix().toDense());
        // gsDebugVar(A.rhs().transpose()   );

        gsInfo<< "." <<std::flush;// Assemblying done

        timer.restart();
        solVectorx  = solver.compute(A.matrix()).solve(A.rhs());

        u_solx.extract(Psix);
        Psi = Psix;

        // in y direction 
        u.setup(bc_y, dirichlet::interpolation, 0);
        // Initialize the system :  identity mapping as initial guess
        A.initSystem();
        setup_time += timer.stop();

        gsInfo<< A.numDofs() <<std::flush;
        timer.restart();
        A.assemble(
        igrad(u, G) * igrad(u, G).tr() *meas(G) //matrix
        ,
         (1.- ff.val()/CoeffDensity) * igrad(u, G) * grad(u_y).tr() *meas(G) //rhs vector
        );
        ma_time += timer.stop();
        gsInfo<< "." <<std::flush;// Assemblying done

        timer.restart();
        solVectory  = solver.compute(A.matrix()).solve(A.rhs());

        slv_time += timer.stop();
        gsInfo<< "." <<std::flush; // Linear solving done

        Psi.patch(0).embed(2);
        u_soly.extract(Psiy);
        Psi.patch(0).coefs().col(1) = Psiy.patch(0).coefs();
        // ..===================================================================
        // Picard loop
        index_t NiterPicard{0};
        gsMatrix<> sv0; //
        solution u_lsol = A.getSolution(u, sv0);
        gsInfo<< A.numDofs() <<std::flush;
        for(int ip{0}; ip<=maxIter; ++ip)
        {
            gsInfo<<std::flush;
            sv0        = solVectorx;
            // .. computes the composition
            auto    ff = A.getCoeff(f,PP);
            //====================================== Update 
            u.setup(bc_x, dirichlet::interpolation, 0);
            // Initialize the system :  identity mapping as initial guess
            A.initSystem();
            setup_time += timer.stop();

            timer.restart();
            A.assemble(
            igrad(u, G) * igrad(u, G).tr() *meas(G) //matrix
            ,
            (1.- ff.val()/CoeffDensity) * igrad(u, G) * grad(u_solx).tr() *meas(G) //rhs vector
            );

            ma_time += timer.stop();

            // gsDebugVar(A.matrix().toDense());
            // gsDebugVar(A.rhs().transpose()   );

            gsInfo<< "." <<std::flush;// Assemblying done

            timer.restart();
            solVectorx  = solver.compute(A.matrix()).solve(A.rhs());
            u_solx.extract(Psix);

            // in y direction 
            u.setup(bc_y, dirichlet::interpolation, 0);
            // Initialize the system :  identity mapping as initial guess
            A.initSystem();
            setup_time += timer.stop();

            gsInfo<< A.numDofs() <<std::flush;
            timer.restart();
            A.assemble(
            igrad(u, G) * igrad(u, G).tr() *meas(G) //matrix
            ,
            (1.- ff.val()/CoeffDensity) * igrad(u, G) * grad(u_soly).tr() *meas(G) //rhs vector
            );
            ma_time += timer.stop();
            gsInfo<< "." <<std::flush;// Assemblying done

            timer.restart();
            solVectory  = solver.compute(A.matrix()).solve(A.rhs());

            slv_time += timer.stop();
            gsInfo<< "." <<std::flush; // Linear solving done
            u_soly.extract(Psiy);

            // Update the mapping
            Psi.patch(0).coefs().col(0) = Psix.patch(0).coefs();
            Psi.patch(0).coefs().col(1) = Psiy.patch(0).coefs();
            // ..===================================================================
            // Check convergence
            
            ++NiterPicard;
            auto l2errRes = math::sqrt(ev.integral( (u_lsol - u_solx).sqNorm() ) );
            auto Ddet     = ev.min(jac(PP).det());
            if ( l2errRes < 1e-8 || ip == maxIter ){
                // ! end Picard loop
                gsInfo<< "\n Niter in Picard : " << ip
                        << ".. L2 MAE residual : "<<std::scientific<<l2errRes
                        << ".. min JAcobian : "<<Ddet<<"..";
                break;
                } //

        }//for loop
        // omp_set_dynamic(0);     // Explicitly disable dynamic teams
        // omp_set_num_threads(1); // Use these threads for later parallel regions

        timer.restart();
        h1err[r]= math::sqrt(ev.integral( ( u_x- u_solx ).sqNorm() + ( u_y- u_soly ).sqNorm() ));
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
        collection.options().setSwitch("base64", false);
        collection.options().setInt("plotElements.resolution", 16);
        collection.options().setInt("numPoints", 100000);
        collection.newTimeStep(&Psi);
        collection.addField(ff_Psi, "density function");
        // collection.addField( ff.val()/CoeffDensity*meas(PP), "MAE_rhs");
        collection.saveTimeStep();
        collection.save();
        gsFileManager::open("ParaviewOutput/solution.pvd");

        // gsWrite(Psi, "Psi_mapping");
        // gsInfo << "Result written in Psi_mapping.xml \n";
    }
    else
        gsInfo << "Done. No output created, re-run with --plot to get a ParaView "
                  "file containing the solution.\n";
    //! [Export visualization in ParaView]

    return EXIT_SUCCESS;


}// end main