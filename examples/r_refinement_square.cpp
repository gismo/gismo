/** @file r_refinement.cpp

    @brief Tutorial on how to use expression assembler to solve a Poisson equation in adaptive mesh.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris & M. BAHARI
*/

//! [Include namespace]
#include <gismo.h>
# include <gsAssembler/gsAdaptiveRefUtils.h>

using namespace std;
using namespace gismo;
//! [Include namespace]


void ProjectionNormalCPoints(gsMultiPatch<>& Psi, gsMultiPatch<> mp){
    // Projection normal of control points (exact geometry)
    int boxMaxNumber = mp.nBoxes();
    for (int boxNumber = 0; boxNumber < boxMaxNumber; ++boxNumber)
    {
        // test if the boundary interface is not an inner interface between patches
        if(!mp.isInterface( patchSide(boxNumber,1) ) ){
        for (int i_x =0; i_x < Psi.patch(boxNumber).basis().boundary(1).size(); ++i_x) // x=0 control points be like (0,:) in this case
        {
            Psi.patch(boxNumber).coef( Psi.patch(boxNumber).basis().boundary(1).at(i_x) ).array()[0] = mp.patch(boxNumber).coef( mp.patch(boxNumber).basis().boundary(1).at(0) ).array()[0];
        }
        }

        if(!mp.isInterface( patchSide(boxNumber,2) ) ){
        for (int i_x =0; i_x < Psi.patch(boxNumber).basis().boundary(2).size(); ++i_x)// x=1 control points be like (1,:) in this case
        {
        Psi.patch(boxNumber).coef( Psi.patch(boxNumber).basis().boundary(2).at(i_x) ).array()[0] = mp.patch(boxNumber).coef( mp.patch(boxNumber).basis().boundary(1).at(0) ).array()[0] + 1.;
        }
        }

        if(!mp.isInterface( patchSide(boxNumber,3) ) ){
        for (int i_x =0; i_x < Psi.patch(boxNumber).basis().boundary(3).size(); ++i_x) // y=0 control points be like (:,0) in this case
        {
        Psi.patch(boxNumber).coef( Psi.patch(boxNumber).basis().boundary(3).at(i_x) ).array()[1] = mp.patch(boxNumber).coef( mp.patch(boxNumber).basis().boundary(1).at(0) ).array()[1];
        }
        }
        if(!mp.isInterface( patchSide(boxNumber,4) ) ){
        for (int i_x =0; i_x < Psi.patch(boxNumber).basis().boundary(4).size(); ++i_x)// y=1 control points be like (:,1) in this case
        {
        Psi.patch(boxNumber).coef( Psi.patch(boxNumber).basis().boundary(4).at(i_x) ).array()[1] = mp.patch(boxNumber).coef( mp.patch(boxNumber).basis().boundary(1).at(0) ).array()[1]+1.;
        }
        }
        }
};


int main(int argc, char *argv[])
{
    //! [Parse command line]
    bool plot          = false;
    index_t numRefine  = 4;// for local refinement:  0 means no local h-refinement
    index_t UnifRefine = 3;// initial refinement: for MAE resolution take at least >=3 for Bejictive mapping 
    index_t DegElevate = 2; // degree Elevation
    index_t NumArMarEl = 1; // Number of ring of cells around marked elements
    index_t maxIter    = 30;
    double eps         = 1e-5; // pinalization coefficient
    double tolPicard   = 1e-8;
    bool last = false, export_b64 =false;
    // ...PNormalCP: Correct the normal part of the mapping and CornersLshape: adjust the corners of the three patches that form L.
    bool PNormalCP{true};
    // --------------- adaptive refinement ---------------
    // Specify cell-marking strategy...
    MarkingStrategy adaptRefCrit = PUCA;
    //MarkingStrategy adaptRefCrit = GARU;
    //MarkingStrategy adaptRefCrit = errorFraction;

    // ... and parameter.
    const real_t adaptRefParam = 0.5;


    gsCmdLine cmd("Tutorial on solving a non-linear Monge-Ampere problem.");
    cmd.addInt("i", "iter", "Maximum number of iterations for the iterative Picard", maxIter);
    cmd.addInt( "e", "degreeElevation",
                "Number of degree elevation steps to perform before solving (0: equalize degree in all directions)", DegElevate);
    cmd.addInt( "r", "uniformRefine", "Number of Uniform h-refinement loops",  numRefine );
    //cmd.addString( "f", "file", "Input XML file", fn );
    cmd.addSwitch("last", "Solve solely for the last level of h-refinement", last);
    cmd.addSwitch("plot", "Create a ParaView visualization file with the solution", plot);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    //gsFileData<> fd(fn);
    //gsInfo << "Loaded file "<< fd.lastPath() <<"\n";
    // .... one single patch
    gsMultiPatch<> mp,  mptp = gsNurbsCreator<>::BSplineSquareGrid(1,1,1, 0.0, 0.0);
    // //... patch 2 (L-shape)
    // gsMultiPatch<> mp,  mptp = gsNurbsCreator<>::BSplineSquareGrid(2,1,1, -1.0, 0.0);
    // mptp.addPatch(gsNurbsCreator<>::BSplineSquare(1, 0.,1.0));

    // ... need regularity to be at least C^1
    mptp.degreeElevate(DegElevate);
    for(size_t i =0; i<mptp.nPatches(); ++i)
        mp.addPatch(gsTHBSpline<2>( dynamic_cast<const gsTensorBSpline<2>&>(mptp.patch(i)) ));
    mp.addAutoBoundaries();

    //..... Test 2
    // Manufactured solition
    gsFunctionExpr<> s("1./(1.+exp((y - x  - 0.2)/0.01))",2);
    // // Manufactured Grad solition
    //gsFunctionExpr<> sP("2.06115362243856e-7*exp(-100.0*x + 100.0*y)/(2.06115362243856e-9*exp(-100.0*x + 100.0*y) + 1.0)**2","-2.06115362243856e-7*exp(-100.0*x + 100.0*y)/(2.06115362243856e-9*exp(-100.0*x + 100.0*y) + 1.0)**2",2);
    // // Right-hand side function
    gsFunctionExpr<> SourceFunc("4.12230724487712e-5*exp(-100.0*x + 100.0*y)/(2.06115362243856e-9*exp(-100.0*x + 100.0*y) + 1.0)**2 - 1.69934170211664e-13*exp(-200.0*x + 200.0*y)/(2.06115362243856e-9*exp(-100.0*x + 100.0*y) + 1.0)**3",2);

    //..... Test 2
    // convex function
    //gsFunctionExpr<> s("0.5*(x**2 + y**2)",2);
    // // Manufactured identity mapping
    gsFunctionExpr<> sN("x","y",2);
    // Right-hand side function : Analytical density function (det(H(u))=f= sigma/rho)
    //gsFunctionExpr<> f("(1./(2.+cos(8.*pi*sqrt((x-0.5-0.25*0.)**2+(y-0.5)**2))))",2);
    //
    //gsFunctionExpr<> f("(1.+ 9./(1.+(10.*sqrt((x-0.7-0.25*0.)**2+(y-0.5)**2)*cos(atan2(y-0.5,x-0.7-0.25*0.) -20.*((x-0.7-0.25*0.)**2+(y-0.5)**2)))**2) )",2);
    //gsFunctionExpr<> f("( 1.+ 5.*exp(-50.*abs((x-0.5-0.25*cos(2.*pi*0.25))**2-(y-0.5-0.5 *sin(2.*pi*0.25))**2- 0.01)))",2);
    gsFunctionExpr<> f("1.+6.*( 1/(1.+exp((y -x  - 0.3)/0.01)) - 1/(1.+exp((y - x  - 0.1)/0.01)) )",2);
    //gsFunctionExpr<> f("(1. + 5./cosh( 5.*((x-sqrt(3)/2)**2+(y-0.5)**2 - (pi/2)**2) )**2 + 5./cosh( 5.*((x+sqrt(3)/2)**2+(y-0.5)**2 - (pi/2)**2) )**2)",2);
    gsInfo<<"Source function "<< f << "\n";// + 5./cosh( 10.*((x-0.2)**2 - 0.9) )

    gsInfo<<"The Initial domain is "<< mp.detail() << "\n";

    gsBoundaryConditions<> bc_mae;
    bc_mae.setGeoMap(mp);
    // For simplicity, set Neumann boundary conditions
    for ( gsMultiPatch<>::const_biterator
            bit = mp.bBegin(); bit != mp.bEnd(); ++bit)
    {
       bc_mae.addCondition( *bit, condition_type::neumann, &sN,0, false);
    }
    //gsDebugVar( bc.allConditions()[0].parametric() );
    gsInfo<<"Boundary conditions:\n"<< bc_mae <<"\n";

    gsBoundaryConditions<> bc;

    //gsOptionList Aopt;

    //! [Refinement]
    gsMultiBasis<> dbasis(mp, true);//true: poly-splines (not NURBS)

    // Elevate and p-refine the basis to order p + numElevate
    // where p is the highest degree in the bases
    //dbasis.setDegree( dbasis.maxCwiseDegree() + numElevate);
    //gsInfo << dbasis.degree(0) << " degree  \n";

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


    //! [Solver loop]
    gsSparseSolver<>::CGDiagonal solver;

    // Elements used for numerical integration
    A.setIntegrationElements(dbasis);
    gsExprEvaluator<> ev(A);

    // Set the geometry map
    geometryMap G = A.getMap(mp);

    // Set pow for BFO method
    auto IGdim     = G.domainDim();
    // Set factor for BFO method
    auto gammaMAE = factorial(G.domainDim());

    // Set the discretization space
    space u = A.getSpace(dbasis);

    // Set the source term
    auto ff = A.getCoeff(f, G);

    // Solution vector and solution variable
    gsMatrix<> solVector;
    solution u_sol = A.getSolution(u, solVector);

    //! [Problem setup] ***-------------------Initailisation for adaptive mapping ---------------------***
    gsMultiPatch<> Psi;
    geometryMap PP = A.getMap(Psi);    
    // Set the source term for Poisson equation
    auto SFunc = A.getCoeff(SourceFunc, PP);

    // Recover manufactured solution for Poisson equation
    auto u_ex = ev.getVariable(s, PP);

    // Recover manufactured density function for Poisson equation
    auto fp = ev.getVariable(f, PP);

    // Set the discretization space // different boundary condition !
    space ru = A.getSpace(dbasis);
    
    // Solution vector and solution variable
    gsMatrix<> rsolVector;
    solution ru_sol = A.getSolution(ru, rsolVector);

    Psi.addAutoBoundaries();
    bc.setGeoMap(Psi);
    // For simplicity, set Dirichlet boundary conditions
    for ( gsMultiPatch<>::const_biterator
            bit = mp.bBegin(); bit != mp.bEnd(); ++bit)
    {
       bc.addCondition( *bit, condition_type::dirichlet, &s,0, false);
    }
    //... END Initailisation

    gsVector<>  h1err(numRefine+1), l2err(numRefine+1);
    gsVector<int>  DoFPDE(numRefine+1);
    gsInfo<< "(dot1=assembled, dot2=solved, dot3=nonlinear_loop,dot4=got_error)\n";
    double setup_time(0), ma_time(0), slv_time(0), err_time(0);
    gsStopwatch timer;
    //... 
    mp.uniformRefine();
    for (int r=0; r<=UnifRefine; ++r)
    {
        dbasis.uniformRefine();
        //mp.uniformRefine();
        Psi.uniformRefine();
    }
    //::::::::::::::::::::      mesh adaptation solver         :::::::::::::::::::::::::
    for (int r=0; r<=numRefine; ++r)
    {
        // Elements used for numerical integration
        A.setIntegrationElements(dbasis);
        gsExprEvaluator<> ev(A);

        // Set the geometry map
        geometryMap G = A.getMap(mp);

        // Set the discretization space
        space u = A.getSpace(dbasis);

        // Set the geometry optimal map
        geometryMap PP = A.getMap(Psi);    
        // Set the discretization space // different boundary condition !
        space ru = A.getSpace(dbasis);
        if (r==0){
            //*********************************************************//

            //dbasis.uniformRefine();
            // mp.uniformRefine();
            //Psi.uniformRefine();
            // Compute the system matrix and right-hand side

            //... nromalisation of density function
            auto CoeffDensity{ev.integral(ff.val() * meas(G))};
            // Initialize the system : start Computing the conductivity coeffeicient ...
            // Compute the Neumann terms defined on physical space
            auto g_N = A.getBdrFunction(G);
            auto Neumann_Int{ev.integralBdrBc(bc_mae.get("Neumann"), g_N.tr() * nv(G) )};
            // ...
            auto CoeffConductivity{Neumann_Int/ev.integral(pow(gammaMAE* CoeffDensity/ff.val(), 1./IGdim) * meas(G))};
            //... end  G.domainDim()+

            // Initialize the system :  identity mapping as initial guess
            A.initSystem();
            setup_time += timer.stop();

            gsInfo<< "\nDoFs: " << A.numDofs() <<std::flush << "\n";

            timer.restart();

            A.assemble(
            igrad(u, G) * igrad(u, G).tr() * meas(G) + eps * u *u.tr()* meas(G) //matrix
            ,
            u*  CoeffConductivity * (-1.)*pow(gammaMAE * CoeffDensity/ff.val(), 1./IGdim) * meas(G) //rhs vector
            );
            
            // Compute the Neumann terms defined on physical space
            //auto g_N = A.getBdrFunction(G);
            A.assembleBdr(bc_mae.get("Neumann"), u * g_N.tr() * nv(G) );

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
            solution u_lsol = A.getSolution(u, sv0);
            for(int ip{0}; ip<=maxIter; ++ip)
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

                // Obtain control points for the gradient of Psi
                A.assemble( v * v.tr() , v * igrad(u_s,G) );
                vsolVector = solver.compute(A.matrix()).solve(A.rhs());
                
                v_sol.extract(Psi);                
                
                auto fp = A.getCoeff(f,PP);

                // ...  0  dirichlet for boundaries
                
                sv0 = solVector;

                // Initialize the system
                A.initSystem();
                setup_time += timer.stop();

                //gsInfo<< A.numDofs() <<std::flush;

                timer.restart();
                // Compute the system matrix and right-hand side ... Monge-Ampere eqaution .....
                
                // .. update Coeffeicient of conductivity
                CoeffConductivity = Neumann_Int/ev.integral(pow( (ilapl(u_sol,G)*ilapl(u_sol,G).tr()).val() + gammaMAE*(CoeffDensity/fp.val() - ihess(u_sol,G).det()), 1./IGdim) * meas(G));

                // MAE system
                A.assemble(
                igrad(u, G) * igrad(u, G).tr() * meas(G) +  eps * u * u.tr()* meas(G) //matrix
                ,
                u * CoeffConductivity * (-1.) * pow( (ilapl(u_sol,G)*ilapl(u_sol,G).tr()).val() + gammaMAE*(CoeffDensity/fp.val() - ihess(u_sol,G).det()), 1./IGdim) * meas(G) //rhs vector
                );
                // Compute the Neumann terms defined on physical space
                auto g_N = A.getBdrFunction(G);
                A.assembleBdr(bc_mae.get("Neumann"), u * g_N.tr() * nv(G) );

                ma_time += timer.stop();

                // gsDebugVar(A.matrix().toDense());
                // gsDebugVar(A.rhs().transpose()   );
                

                gsInfo<< " ." <<std::flush;// Assemblying done

                timer.restart();
                solver.compute( A.matrix() );
                solVector = solver.solve(A.rhs());
                
                slv_time += timer.stop();

                gsInfo<< "." <<std::flush; // Linear solving done

                // omp_set_dynamic(0);     // Explicitly disable dynamic teams
                // omp_set_num_threads(1); // Use these threads for later parallel regions

                ++NiterPicard;
                auto l2errRes = math::sqrt(ev.integral( ( igrad(u_lsol,G) - igrad(u_sol,G) ).sqNorm() * meas(G) ));
                if ( l2errRes < tolPicard || ip == maxIter ){
                    // ! end Picard loop
                    gsInfo<< "\n Niter in Picard : " << NiterPicard << ".. L2 residual : "<<std::scientific<<l2errRes<<"\n";
                    break; 
                    } // 
            }//for loop
            // omp_set_dynamic(0);     // Explicitly disable dynamic teams
            // omp_set_num_threads(1); // Use these threads for later parallel regions

            timer.restart();

            // ... correct boundary
            if (PNormalCP)
                ProjectionNormalCPoints(Psi, mp);
            if(mp.nPatches()>1){
            Psi.addInterface(0,2,1,1);
            Psi.addInterface(1,4,2,3);
            }
            Psi.addAutoBoundaries();
            gsInfo<<"The PDE domain is "<< Psi.detail() << "\n";
            gsInfo<<"Boundary conditions:\n"<< bc <<"\n";
        }
        gsInfo << "Patches: "<< Psi.nPatches() <<", degree: "<< dbasis.minCwiseDegree() <<"\n";
        //::::::::::::::::::::   Poisson equation - (manufactured exact solution)         :::::::::::::::::::::::::

        ru.setup(bc, dirichlet::l2Projection, 0);

        // Compute the system matrix and right-hand side
        // Initialize the system
        A.initSystem();
        setup_time += timer.stop();

        gsInfo<< "Solving PDEs " <<std::flush;
        gsInfo<< A.numDofs() <<std::flush;

        timer.restart();

        A.assemble(
        igrad(ru, PP) * igrad(ru, PP).tr() * meas(PP) //matrix
        ,
        ru * SFunc * meas(PP) //rhs vector
        );

        ma_time += timer.stop();

        // gsDebugVar(A.matrix().toDense());
        // gsDebugVar(A.rhs().transpose()   );

        gsInfo<< "." <<std::flush;// Assemblying done

        timer.restart();
        solver.compute( A.matrix() );
        rsolVector = solver.solve(A.rhs());

        slv_time += timer.stop();

        gsInfo<< "." <<std::flush; // Linear solving done

        // omp_set_dynamic(0);     // Explicitly disable dynamic teams
        // omp_set_num_threads(1); // Use these threads for later parallel regions
        DoFPDE[r] = A.numDofs();

        timer.restart();
        l2err[r]= math::sqrt( ev.integral( (u_ex - ru_sol).sqNorm() * meas(PP) ) );

        h1err[r]= math::sqrt(ev.integral( ( igrad(u_ex) - igrad(ru_sol,PP) ).sqNorm() * meas(PP) ));
        err_time += timer.stop();
        gsInfo<< ". " <<std::flush; // Error computations done

        if(r < numRefine){
        //! [beginRefLoop]
            gsInfo << "====== Loop " << r << " of "
                    <<numRefine<< " ======" << "\n";
        // --------------- error estimation/computation ---------------
        // Get the element-wise norms.
        ev.integralElWise( ( ilapl(ru_sol, PP)+ SFunc ).sqNorm()*meas(PP) );
        const std::vector<real_t> & eltErrs  = ev.elementwise();
        //! [errorComputation]

        //! [adaptRefinementPart]
        // Mark elements for refinement, based on the computed local errors and
        // the refinement-criterion and -parameter.
        std::vector<bool> elMarked( eltErrs.size() );
        gsMarkElementsForRef( eltErrs, adaptRefCrit, adaptRefParam, elMarked);
        gsInfo <<"Marked "<< std::count(elMarked.begin(), elMarked.end(), true) <<" elements.\n";
        // Refine the marked elements with a 1-ring of cells around marked elements
        gsRefineMarkedElements( dbasis, elMarked, NumArMarEl);
        gsRefineMarkedElements( Psi, elMarked, NumArMarEl);
        }
    }
    //! [Solver loop]    


    timer.stop();
    gsInfo<<"\n\nTotal time: "<< setup_time+ma_time+slv_time+err_time <<"\n";
    gsInfo<<"     Setup: "<< setup_time <<"\n";
    gsInfo<<"  Assembly: "<< ma_time    <<"\n";
    gsInfo<<"   Solving: "<< slv_time   <<"\n";
    gsInfo<<"     Norms: "<< err_time   <<"\n";
    

    //! [Error and convergence rates]
    gsInfo<< "\nDoF_PDE = "<<std::scientific<<DoFPDE.transpose()<<"\n";
    gsInfo<< "L2_error = "<<std::scientific<<std::setprecision(3)<<l2err.transpose()<<"\n";
    gsInfo<< "H1_error= "<<std::scientific<<std::setprecision(3)<<h1err.transpose()<<"\n";

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
        collection.options().setSwitch("base64", export_b64);
        collection.options().setInt("plotElements.resolution", 16);
        collection.options().setInt("numPoints", 10000);
        collection.newTimeStep(&Psi);
        collection.addField(ru_sol,"numerical solution");
        collection.addField(igrad(ru_sol,PP),"gradient_numerical solution");
        collection.addField(jac(PP).det(), "Jacobian function");
        collection.addField(u_ex, "exact solution");
        collection.addField(fp,"Density function");
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