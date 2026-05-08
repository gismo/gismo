/** @file multipatch_r_refinement.cpp TODO

    @brief Tutorial on how to use expression assembler to builde multipatch adaptive-mesh solver from anayltic density function

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
// Project control points following  normal direction at the boundaries for square domain (Exact square recovery after refinement)
// void aforcing_c0(gsMultiPatch<>& mp, gsMultiPatch<>& Psi)
// {
//     //TODO correct control points of Psi to be c0 with respect to multipatch interfaces for Psi based on the interfaces of mp (mp only give us connection not points), while Psi is multipatch but only unit-square (0,1)^2 means if interface 1/2 Psi control poitns cpts at interface 1(x=0) should have same as for Psi cpts right (x=1) which leads to correct only in y direction
//     // normal Projection of control points (exact geometry)
//     auto interfaces = mp.interfaces();
//     #pragma omp parallel for
//     for (auto interface : interfaces)
//     {
//         index_t bndIter  = 1;
//         auto leftbox = interface.first();
//         auto secbox = interface.second();
//         auto i_dir = leftbox.direction();
//         // x=0 control points be like (0,:) in this case
//         for (int i_x =0; i_x < Psi.patch(leftbox).basis().boundary(bndIter).size(); ++i_x) 
//         {
//             // moyen
//             auto lval = Psi.patch(leftbox).coef( Psi.patch(leftbox).basis().boundary(bndIter).at(i_x) ).array()[i_dir];
//             lval = lval+ Psi.patch(secbox).coef( Psi.patch(secbox).basis().boundary(bndIter).at(i_x) ).array()[i_dir];
//             // correction 
//             Psi.patch(leftbox).coef( Psi.patch(leftbox).basis().boundary(bndIter).at(i_x) ).array()[i_dir] = lval;
//             Psi.patch(secbox).coef( Psi.patch(secbox).basis().boundary(bndIter).at(i_x) ).array()[i_dir] = lval;
//         }
//     }
// }

int main(int argc, char *argv[])
{
    //! [Parse command line]
    bool plot          = false;
    double Intensity   = 0.;
    index_t numRefine  = 3;
    index_t numElevate = 0;
    index_t maxIter    = 50;
    bool last          = false;
    // Specify the file path
    //std::string fn("pde/circle.xml");
    // std::string fn("surfaces/cylinder.xml"); 
    std::string fn( "pde/solovev_relaxed.xml" );
    // load the file
    gsFileData<> fd(fn);
    gsInfo << "Loaded file "<< fd.lastPath() <<"\n";
    gsMultiPatch<> mpLeft;// = gsNurbsCreator<>::BSplineSquareGrid(1,1,1, 0.0, 0.0);// Initial geometry
    fd.getId(5,mpLeft);
    //.. density 
    // std::string frho("density_function.xml");
    // gsFileData<> strho(frho);
    // gsMultiPatch<> monitor_fct;
    // strho.getId(5,monitor_fct);

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

    // .... one single patch
   gsMultiPatch<> mp = gsNurbsCreator<>::BSplineSquareGrid(1,1,1, 0.0, 0.0);

    //identity mapping in each direction
    gsFunctionExpr<> sx("x",2);
    gsFunctionExpr<> sy("y",2);
    gsFunctionExpr<> s0("0.",2);
    // Manufactured identity mapping
    // gsFunctionExpr<> sN("x","y",2);
    // Right-hand side function : Analytical density function (det(H(u))=f= sigma/rho)
    gsFunctionExpr<> f("1./(2.+cos(8.*pi*sqrt((x-0.5-0.25*0.)**2+(y-0.5)**2)))",2);
    //gsFunctionExpr<> f("1.+6.*( 1/(1.+exp((y -x  - 0.3)/0.01)) - 1/(1.+exp((y - x  - 0.1)/0.01)) )",2);
    //Manufactured density function
    //gsFunctionExpr<> f("1.+6.*( 1/(1.+exp(( (y-0.98)**2+(x-0.899)**2  - 0.002)/0.001)) + 1/(1.+exp((y -x  - 0.25)/0.001)) - 1./(1.+exp((y - x  - 0.15)/0.001)) +  0/(1.+exp((y - 1.0)/0.001)) - 0./(1.+exp((y - 0.975)/0.001))  +  1/(1.+exp((x - 1.0)/0.001)) - 1./(1.+exp((x - 0.95)/0.001)) )",2);    
    //gsFunctionExpr<> f("(1.+ 9./(1.+(10.*sqrt((x-0.7-0.25*0.)**2+(y-0.5)**2)*cos(atan2(y-0.5,x-0.7-0.25*0.) -20.*((x-0.7-0.25*0.)**2+(y-0.5)**2)))**2) )",2);
    //gsFunctionExpr<> f("( 1.+ 5.*exp(-50.*abs((x-0.5-0.25*cos(2.*pi*0.25))**2-(y-0.5-0.5 *sin(2.*pi*0.25))**2- 0.01)))",2);
    //gsFunctionExpr<> f("(1. + 5./cosh( 5.*((x-sqrt(3)/2)**2+(y-0.5)**2 - (pi/2)**2) )**2 + 5./cosh( 5.*((x+sqrt(3)/2)**2+(y-0.5)**2 - (pi/2)**2) )**2)",2);
    // gsFunctionExpr<> f("(1.+10.(1.+0.*x*(abs(y-0.5)<0.3))/cosh( 10.*(x -0.5 ) )**2)",2);
    gsInfo<<"Source function "<< f << "\n";

    gsInfo<<"The domain is "<< mp.detail() << "\n";

    gsBoundaryConditions<> bc_x;
    bc_x.setGeoMap(mp);
    bc_x.addCondition(0,1, condition_type::dirichlet, &sx,0,false);
    bc_x.addCondition(0,2, condition_type::dirichlet, &sx,0,false);
    bc_x.addCondition(0,3, condition_type::neumann, &s0,0,false);
    bc_x.addCondition(0,4, condition_type::neumann, &s0,0,false);
//     for ( gsMultiPatch<>::const_biterator
//             bit = mp.bBegin(); bit != mp.bEnd(); ++bit)
//    {
//        bc_x.addCondition( *bit, condition_type::dirichlet, &sx );
//    }
    gsBoundaryConditions<> bc_y;
    bc_y.setGeoMap(mp);
    bc_y.addCondition(0,1, condition_type::neumann, &s0,0,false);
    bc_y.addCondition(0,2, condition_type::neumann, &s0,0,false);
    bc_y.addCondition(0,3, condition_type::dirichlet, &sy,0,false);
    bc_y.addCondition(0,4, condition_type::dirichlet, &sy,0,false);
//     for ( gsMultiPatch<>::const_biterator
//             bit = mp.bBegin(); bit != mp.bEnd(); ++bit)
//    {
//        bc_y.addCondition( *bit, condition_type::dirichlet, &sy );
//    }

    //! [Refinement]
    gsMultiBasis<> dbasis(mp, true);//true: poly-splines (not NURBS)

    // Elevate and p-refine the basis to order p + numElevate
    // where p is the highest degree in the bases
    dbasis.setDegree( dbasis.maxCwiseDegree() + numElevate);

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
    A.options().setSwitch("SameElement",false);
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

    // idensity mapping
    // auto u_I = ev.getVariable(sN, G);
    auto u_x = ev.getVariable(sx, G);
    auto u_y = ev.getVariable(sy, G);

    // density vector and density variable
    gsMatrix<> solVectorho;
    solution u_solho = A.getSolution(u, solVectorho);
    // Solution vector and solution variable in x direction
    gsVector<gsMatrix<>> solVectorxMP(mpLeft.nPatches());
    gsVector<gsMatrix<>> solVectoryMP(mpLeft.nPatches());

    gsMatrix<> solVectorx;
    solution u_solx = A.getSolution(u, solVectorx);
    // Solution vector and solution variable in y direction
    gsMatrix<> solVectory;
    solution u_soly = A.getSolution(u, solVectory);

    // ... Adaptive mapping
    gsMultiPatch<> Psix;
    gsMultiPatch<> Psiy;

    gsMultiPatch<> rho, Psi, newLeft, TargetPsi;

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
    double CoeffDensity = 20.; // maximum of the density function for normalization

    gsVector<>  h1err(numRefine+1); //l2err(numRefine+1).
    gsInfo<< "(dot1=assembled, dot2=solved, dot3=nonlinear_loop,dot4=got_error)\n"
        "\n ";
    double setup_time(0), ma_time(0), slv_time(0), err_time(0);
    gsStopwatch timer;

    //::::::::::::::::::::      mesh adaptation solver         :::::::::::::::::::::::::
    for (int r=0; r<=numRefine; ++r)
    {
        dbasis.uniformRefine();
        mp.uniformRefine();
        newLeft.clear();
        TargetPsi.clear();
        rho.clear();
        // we first build a identity mapping in the parametric domain, then we compose it with the geometry map to get the new geometry map
        for (size_t n_patch = 0; n_patch < mpLeft.nPatches(); ++n_patch){
            gsMatrix<> intGrid             = dbasis.basis(0).anchors();
            // Evaluate f at the Greville points
            gsGeometry<>::uPtr interpolant = dbasis.basis(0).interpolateData(intGrid, intGrid);
            // extract the mapping
            TargetPsi.addPatch(give(interpolant));   
            //---------------------------------------------
            // Set density function in computational domain
            //---------------------------------------------
            gsMultiPatch<> rhotmp;
            geometryMap GLeft = A.getMap(mpLeft.patch(n_patch));
            auto ff = A.getCoeff(f, GLeft);
            // ---- manipulation of density function ----
            auto empldensity = (ev.max(abs(ff.val()))-ev.min(abs(ff.val())));
            double  int_uh_0 = 1.;
            double  int_uh_1 = 1.;

            if (empldensity < 1e-5|| Intensity <= 1. )
            {
                gsInfo << " rho = 1.~~";
            }
            else{
                int_uh_0     = (Intensity-1.)/empldensity;
                int_uh_1     = (ev.max(abs(ff.val()))-Intensity*ev.min(abs(ff.val())))/empldensity;
            }
            // end of manipulation
            u.setup(bc_x, dirichlet::interpolation, -1);

            A.initSystem();
            A.assemble(
            u * u.tr() * meas(G) //matrix
            ,
            u * (int_uh_0*ff[0].val()+int_uh_1) * meas(G) //rhs vector
            );        
            solVectorho  = solver.compute(A.matrix()).solve(A.rhs());
            u_solho.extract(rhotmp);
            rho.addPatch(rhotmp.patch(0));
        }
        for (size_t n_patch = 0; n_patch < mpLeft.nPatches(); ++n_patch){
            gsInfo << "patch number = "<< n_patch<< "----------------\n";
            //----------------------------------------------
            // New density function defined in square
            //----------------------------------------------
            auto frho = A.getCoeff(rho.patch(n_patch), G);

            //... nromalisation of density function
            CoeffDensity = std::max(CoeffDensity, ev.max(frho)+1.);
            //----------------------------------------------
            //Start MMPDE refinement
            //----------------------------------------------
            u.setup(bc_x, dirichlet::interpolation, 0);
            // Initialize the system :  identity mapping as initial guess
            A.initSystem();
            setup_time += timer.stop();

            gsInfo<< "DoFs:" << A.numDofs() <<std::flush<< "\n" ;
            timer.restart();
            A.assemble(
            grad(u) * grad(u).tr() //matrix
            ,
            (1.- frho.val()/CoeffDensity) * grad(u) * grad(u_x).tr() //rhs vector
            );

            ma_time += timer.stop();
            // gsInfo<< A.matrix()  << "." <<std::flush<< "\n" ;

            // gsDebugVar(A.matrix().toDense());
            // gsDebugVar(A.rhs().transpose()   );

            gsInfo<< "." <<std::flush;// Assemblying done

            timer.restart();
            solVectorx  = solver.compute(A.matrix()).solve(A.rhs());
            solVectorxMP[n_patch] = solVectorx;

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
            grad(u) * grad(u).tr() //matrix
            ,
            (1.- frho.val()/CoeffDensity) * grad(u) * grad(u_y).tr() //rhs vector
            );
            ma_time += timer.stop();
            gsInfo<< "." <<std::flush;// Assemblying done

            timer.restart();
            solVectory  = solver.compute(A.matrix()).solve(A.rhs());
            solVectoryMP[n_patch] = solVectory;

            slv_time += timer.stop();
            gsInfo<< "." <<std::flush; // Linear solving done
            u_soly.extract(Psiy);
            for(size_t Mp=0; Mp<Psi.nPatches(); ++Mp){
                Psi.patch(Mp).embed(2);
                Psi.patch(Mp).coefs().col(1) = Psiy.patch(Mp).coefs();
            }
            // extract the multipatch mampping in computational domain
            TargetPsi.patch(n_patch).coefs() = Psi.patch(0).coefs(); 
        }
            // ..===================================================================
            // .. correct interfaces to be conforming : C-1-> C0 -> C1
            // ..===================================================================
        
        // Picard loop
        gsInfo << "Step two in r-refinement" << "\n";
        index_t NiterPicard{0};
        gsMatrix<> sv0; //
        solution u_lsol = A.getSolution(u, sv0);
        gsInfo<< A.numDofs() <<std::flush;
        std::vector<int> patch_iter(mpLeft.nPatches());//list of patches

        for (size_t n_patch = 0; n_patch < mpLeft.nPatches(); ++n_patch)
                patch_iter[n_patch] = n_patch;
        for(int ip{0}; ip<=maxIter; ++ip)
        {
        for (size_t n_patch : patch_iter){
                gsInfo<<std::flush;
                sv0        = solVectorx;
                // Set the geometry map
                geometryMap PP = A.getMap(TargetPsi.patch(n_patch));
                // .. computes the composition
                auto    frho = A.getCoeff(rho.patch(n_patch), PP);
                // update solutions
                solVectorx = solVectorxMP[n_patch];
                solVectory = solVectoryMP[n_patch];

                //====================================== Update 
                u.setup(bc_x, dirichlet::interpolation, 0);
                // Initialize the system :  identity mapping as initial guess
                A.initSystem();
                setup_time += timer.stop();

                timer.restart();
                A.assemble(
                grad(u) * grad(u).tr() //matrix
                ,
                (1.- frho.val()/CoeffDensity) * grad(u) * grad(u_solx).tr() //rhs vector
                );

                ma_time += timer.stop();

                // gsDebugVar(A.matrix().toDense());
                // gsDebugVar(A.rhs().transpose()   );

                gsInfo<< "." <<std::flush;// Assemblying done

                timer.restart();
                solVectorx  = solver.compute(A.matrix()).solve(A.rhs());
                u_solx.extract(Psix);
                solVectorxMP[n_patch] = solVectorx;
                // in y direction 
                u.setup(bc_y, dirichlet::interpolation, 0);
                // Initialize the system :  identity mapping as initial guess
                A.initSystem();
                setup_time += timer.stop();

                // gsInfo<< A.numDofs() <<std::flush;
                timer.restart();
                A.assemble(
                grad(u) * grad(u).tr() //matrix
                ,
                (1.- frho.val()/CoeffDensity) * grad(u) * grad(u_soly).tr() //rhs vector
                );
                ma_time += timer.stop();
                gsInfo<< "." <<std::flush;// Assemblying done

                timer.restart();
                solVectory  = solver.compute(A.matrix()).solve(A.rhs());
                solVectoryMP[n_patch] = solVectory;

                slv_time += timer.stop();
                gsInfo<< "." <<std::flush; // Linear solving done
                u_soly.extract(Psiy);

                // Update the mapping
                Psi.patch(0).coefs().col(0) = Psix.patch(0).coefs();
                Psi.patch(0).coefs().col(1) = Psiy.patch(0).coefs();
                // extract the multipatch mampping in computational domain
                TargetPsi.patch(n_patch).coefs() = Psi.patch(0).coefs(); 
                // ..===================================================================
                // Check convergence                
                ++NiterPicard;
                auto l2errRes = math::sqrt(ev.integral( (u_lsol - u_solx).sqNorm() ) );
                auto Ddet     = ev.min(jac(PP).det());
                if ( l2errRes < 1e-8 || ip == maxIter ){
                    // ! end Picard loop
                    // patch_iter.erase(patch_iter.begin() + n_patch);
                    gsInfo<< "\n Niter in Picard : " << ip
                            << ".. L2 MAE residual : "<<std::scientific<<l2errRes
                            << ".. min JAcobian : "<<Ddet<<"..";
                    break;
                    } //
            }//for loop
            // ..===================================================================
            // .. correct interfaces to be conforming : C-1-> C0 -> C1
            // ..===================================================================

            // omp_set_dynamic(0);     // Explicitly disable dynamic teams
            // omp_set_num_threads(1); // Use these threads for later parallel regions

            timer.restart();
            h1err[r]= math::sqrt(ev.integral( ( u_x- u_solx ).sqNorm() + ( u_y- u_soly ).sqNorm() ));
            err_time += timer.stop();
            gsInfo<< ". " <<std::flush; // Error computations done

            gsStopwatch timer;
            timer.restart();
        }   
    } //for loop
    //! [Solver loop]    
    // ----------------------
    //.. composition 
    gsInfo<<"<Col> computes composition";
    for(size_t n_patch = 0; n_patch < mpLeft.nPatches(); ++n_patch){
        gsMatrix<> intGrid             = dbasis.basis(0).anchors();
        // Evaluate f at the Greville points
        gsMatrix<> intfavlues          = TargetPsi.patch(n_patch).eval(intGrid);
        intfavlues                     = intfavlues.cwiseMax(0).cwiseMin(1);
        gsMatrix<> fValues             = mpLeft.patch(n_patch).eval(intfavlues);
        gsGeometry<>::uPtr interpolant = dbasis.basis(0).interpolateData(fValues, intGrid);
        // extract the mapping
        newLeft.addPatch(give(interpolant));
    }

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
        //! [Refinement]
        gsMultiBasis<> dbasis(newLeft, true);//true: poly-splines (not NURBS)
        gsExprAssembler<> A(1,1);
        A.setIntegrationElements(dbasis);
        gsExprEvaluator<> ev(A);
        // geometry adaptive map
        geometryMap PPFinal = A.getMap(newLeft);

        auto ff_Psi = A.getCoeff(f, PPFinal);
        gsInfo<<"Plotting in Paraview...\n";
        gsParaviewCollection collection("ParaviewOutput/solution", &ev);
        collection.options().setSwitch("plotElements", true);
        collection.options().setSwitch("base64", false);
        collection.options().setInt("plotElements.resolution", 16);
        collection.options().setInt("numPoints", 100000);
        collection.newTimeStep(&newLeft);
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