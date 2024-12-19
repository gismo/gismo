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
    bool plot          = false;
    index_t numRefine  = 3;// for local refinement:  0 means no local h-refinement
    index_t UnifRefine = 3;// initial refinement: for MAE resolution take at least >=3 for Bejictive mapping 
    index_t DegElevate = 2; // degree Elevation
    index_t NumArMarEl = 1; // Number of ring of cells around marked elements
    index_t maxIter    = 30;
    double eps         = 1e-5; // pinalization coefficient
    double tolPicard   = 1e-8;
    bool last = false, export_b64 =false;
    // ...PNormalCP: Correct the normal part of the mapping and CornersLshape: adjust the corners of the three patches that form L.
    bool PNormalCP     = true;
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
    gsMultiPatch<> Psi, mp, mptp = gsNurbsCreator<>::BSplineSquareGrid(1,1,1, 0.0, 0.0);
    // ... need regularity to be at least C^1
    mptp.degreeElevate(DegElevate);
    for(size_t i =0; i<mptp.nPatches(); ++i)
        mp.addPatch(gsTHBSpline<2>( dynamic_cast<const gsTensorBSpline<2>&>(mptp.patch(i)) ));
    mp.addAutoBoundaries();
    
    // // Identity mapping stay fix
    gsFunctionExpr<> sN("x","y",2);

    /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    //..... Test 1 : POISSON EQUATION
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
    // Define Stabilization method
    auto Stabilizationtype = stabilizerCDR::none;
    // convection coefficient
    gsFunctionExpr<> coeff_conv("0.","0.",2);
    // diffusion coefficient:
    gsFunctionExpr<> coeff_diff("1.","0","0","1.",2);
    // For a posterior error estimate
    gsFunctionExpr<> coeff_diffMax("1.",2);
    // reaction coefficient:
    gsFunctionExpr<> coeff_reac("0",2);
    /* *********************** Test 1 *********************** */
    // // Define  Dirichlet boundary conditions
    // gsFunctionExpr<> Dg("1./(1.+exp((y - x  - 0.2)/0.01))", 2);
    // // Manufactured solition
    // gsFunctionExpr<> s("1./(1.+exp((y - x  - 0.2)/0.01))",2);
    // // // Right-hand side function
    // gsFunctionExpr<> SourceFunc("4.12230724487712e-5*exp(-100.0*x + 100.0*y)/(2.06115362243856e-9*exp(-100.0*x + 100.0*y) + 1.0)**2 - 1.69934170211664e-13*exp(-200.0*x + 200.0*y)/(2.06115362243856e-9*exp(-100.0*x + 100.0*y) + 1.0)**3",2);
    // // Manufactured density function
    // gsFunctionExpr<> f("1.+6.*( 1/(1.+exp((y -x  - 0.3)/0.01)) - 1/(1.+exp((y - x  - 0.1)/0.01))  )",2);
    /* *********************** Test 2 *********************** */
    // Define  Dirichlet boundary conditions
    gsFunctionExpr<> Dg("1./(1.+exp(100 * ( x**2 + (y-0.5)**2-0.75*sin(pi*y)) ))", 2);
    // Manufactured solition
    gsFunctionExpr<> s("1./(1.+exp(100 * ( x**2 + (y-0.5)**2-0.75*sin(pi*y)) ))",2);
    // // Right-hand side function
    gsFunctionExpr<> SourceFunc("40000.0*x**2*exp(100*x**2 + 100*(y - 0.5)**2 - 75.0*sin(pi*y))/(exp(100*x**2 + 100*(y - 0.5)**2 - 75.0*sin(pi*y)) + 1.0)**2 - 80000.0*x**2*exp(200*x**2 + 200*(y - 0.5)**2 - 150.0*sin(pi*y))/(exp(100*x**2 + 100*(y - 0.5)**2 - 75.0*sin(pi*y)) + 1.0)**3 + 1.0*(75.0*pi**2*sin(pi*y) + 200)*exp(100*x**2 + 100*(y - 0.5)**2 - 75.0*sin(pi*y))/(exp(100*x**2 + 100*(y - 0.5)**2 - 75.0*sin(pi*y)) + 1.0)**2 + 1.0*(40000*(y - 0.375*pi*cos(pi*y) - 0.5)**2)*exp(100*x**2 + 100*(y - 0.5)**2 - 75.0*sin(pi*y))/(exp(100*x**2 + 100*(y - 0.5)**2 - 75.0*sin(pi*y)) + 1.0)**2 + 200.0*exp(100*x**2 + 100*(y - 0.5)**2 - 75.0*sin(pi*y))/(exp(100*x**2 + 100*(y - 0.5)**2 - 75.0*sin(pi*y)) + 1.0)**2 - 80000.0*(y - 0.375*pi*cos(pi*y) - 0.5)**2*exp(200*x**2 + 200*(y - 0.5)**2 - 150.0*sin(pi*y))/(exp(100*x**2 + 100*(y - 0.5)**2 - 75.0*sin(pi*y)) + 1.0)**3",2);
    // Manufactured density function
    gsFunctionExpr<> f("1.+ 0.*exp(-100 * ( x**2 + (y-0.5)**2-0.75*sin(pi*y))**2 )",2);

    /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // //..... Test 2 ADVECTION DUFFFUSION
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
    // // Define Stabilization method
    //auto Stabilizationtype = stabilizerCDR::SUPG;
    // // Define  Dirichlet boundary conditions
    // gsFunctionExpr<> Dg("if( y <= 0.2*(1.-x), 1,0)", 2);
    // // Manufactured solition
    // gsFunctionExpr<> s("1./(1.+exp((y - x  - 0.2)/0.001))",2);
    // // convection coefficient:
    // gsFunctionExpr<> coeff_conv("cos(pi/4)","sin(pi/4)",2);
    // // diffusion coefficient:
    // gsFunctionExpr<> coeff_diff("0.000001","0","0","0.000001",2);
    // // For a posterior error estimate
    //gsFunctionExpr<> coeff_diffMax("0.000001",2);
    // // reaction coefficient:
    // gsFunctionExpr<> coeff_reac("0",2);
    // // // Right-hand side function
    // gsFunctionExpr<> SourceFunc("0.",2);
    // //Manufactured density function
    // gsFunctionExpr<> f("1.+6.*( 1/(1.+exp((y -x  - 0.3)/0.01)) - 1/(1.+exp((y - x  - 0.1)/0.01)) + 1/(1.+exp((0.9-x)/0.01)) )",2);


    gsInfo<<"Source function "<< f << "\n";

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

    typedef gsExprAssembler<>::geometryMap geometryMap;
    typedef gsExprAssembler<>::variable    variable;
    typedef gsExprAssembler<>::space       space;
    typedef gsExprAssembler<>::solution    solution;

    //! [Solver loop]
    gsSparseSolver<>::CGDiagonal solver;

    gsVector<>  h1err(numRefine+1), l2err(numRefine+1);
    gsVector<int>  DoFPDE(numRefine+1);
    gsInfo<< "(dot1=assembled, dot2=solved, dot3=nonlinear_loop,dot4=got_error)\n";
    double setup_time(0), ma_time(0), slv_time(0), err_time(0);    
    gsStopwatch timer;

     //! [Refinement]
    mp.uniformRefine();
    for (int r=0; r<=UnifRefine; ++r)
    {
        dbasis.uniformRefine();
        //mp.uniformRefine();
        Psi.uniformRefine();
    }
    //::::::::::::::::::::      mesh adaptation solver         :::::::::::::::::::::::::

    //! [Problem setup]
    gsExprAssembler<> A(1,1);
    //A.setOptions(Aopt);
    gsInfo<<"Active options:\n"<< A.options() <<"\n";

    // Elements used for numerical integration
    A.setIntegrationElements(dbasis);
    gsExprEvaluator<> ev(A);

    // Set the discretization space
    space u = A.getSpace(dbasis); 

    // Set the geometry map
    geometryMap G = A.getMap(mp);

    // Set pow for BFO method
    auto IGdim     = G.domainDim();
    // Set factor for BFO method
    auto gammaMAE = factorial(G.domainDim());

    // Set the source term
    auto ff = A.getCoeff(f, G);

    // Solution vector and solution variable
    gsMatrix<> solVector;
    solution u_sol = A.getSolution(u, solVector);
    
    /* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ###   Step 2: r_refinement
     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
    //... nromalisation of density function
    auto CoeffDensity{ev.integral(ff.val() * meas(G))};
    // Initialize the system : start Computing the conductivity coeffeicient ...
    // Compute the Neumann terms defined on physical space
    auto g_N = A.getBdrFunction(G);
    auto Neumann_Int{ev.integralBdrBc(bc_mae.get("Neumann"), g_N.tr() * nv(G) )};
    // ...
    auto CoeffConductivity{Neumann_Int/ev.integral(pow(IGdim*IGdim+gammaMAE * CoeffDensity/ff.val(), 1./IGdim) * meas(G))};
    //... end 

    //u.setup(bc_mae, dirichlet::l2Projection, 0);
    // Initialize the system :  identity mapping as initial guess
    A.initSystem();
    setup_time += timer.stop();

    gsInfo<< "\nDoFs: " << A.numDofs() <<std::flush << "\n";

    timer.restart();

    A.assemble(
    igrad(u, G) * igrad(u, G).tr() * meas(G) + eps * u *u.tr()* meas(G) //matrix
    ,
    u*  CoeffConductivity * (-1.)*pow(IGdim*IGdim+gammaMAE * CoeffDensity/ff.val(), 1./IGdim) * meas(G) //rhs vector
    );
    
    // Compute the Neumann terms defined on physical space
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
        space v = A.getSpace(dbasis);
        gsMatrix<> vsolVector;
        solution v_sol = A.getSolution(v, vsolVector);
        A.initSystem(2);

        // Obtain control points for the gradient of Psi
        A.assemble( v * v.tr() , v * igrad(u_s,G) );
        vsolVector = solver.compute(A.matrix()).solve(A.rhs());
        
        v_sol.extract(Psi);                
        
        // Set the geometry optimal map
        geometryMap PP = A.getMap(Psi);
        auto ff = A.getCoeff(f,PP);

        // ...  0  dirichlet for boundaries
        sv0 = solVector;
        
        //u.setup(bc_mae, dirichlet::l2Projection, 0);

        // Initialize the system
        A.initSystem();
        setup_time += timer.stop();
        timer.restart();
        // Compute the system matrix and right-hand side ... Monge-Ampere eqaution .....
        
        // .. update Coeffeicient of conductivity
        CoeffConductivity = Neumann_Int/ev.integral(pow( (ilapl(u_sol,G)*ilapl(u_sol,G).tr()).val() + gammaMAE*(CoeffDensity/ff.val() - ihess(u_sol,G).det()), 1./IGdim) * meas(G));

        // MAE system
        A.assemble(
        igrad(u, G) * igrad(u, G).tr() * meas(G) +  eps * u * u.tr()* meas(G) //matrix
        ,
        u * CoeffConductivity * (-1.) * pow( (ilapl(u_sol,G)*ilapl(u_sol,G).tr()).val() + gammaMAE*(CoeffDensity/ff.val() - ihess(u_sol,G).det()), 1./IGdim) * meas(G) //rhs vector
        );

        // Compute the Neumann terms defined on physical space
        auto g_N = A.getBdrFunction(G);
        A.assembleBdr(bc_mae.get("Neumann"), u * g_N.tr() * nv(G) );

        ma_time += timer.stop();            

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
    /* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ###   Step 3: Correct boundary
     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
    if (PNormalCP)
        ProjectionNormalCPoints(Psi, mp);
    //  construct boundaries and interfaces
    if(mp.nPatches()>1){
    Psi.addInterface(0,2,1,1);
    Psi.addInterface(0,4,2,3);
    Psi.addAutoBoundaries();
    //Psi.computeTopology();
    }
    Psi.addAutoBoundaries();

    /* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ###   Step 4: Define hierarchical adaptive mapping
     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
    //gsMultiPatch<> Psi = Psi;
    // for(size_t i =0; i<Psi.nPatches(); ++i)
    //     Psi.addPatch(gsTHBSpline<2>( dynamic_cast<const gsTensorBSpline<2>&>(Psi.patch(i)) ));
    // //Psi.addAutoBoundaries();
    // Psi.computeTopology();
    gsInfo<<"The PDE domain is "<< Psi.detail() << "\n";

    // Set Dirichlet boundary conditions
    gsBoundaryConditions<> bc;
    bc.setGeoMap(Psi);
    // given by exact solution Dg on all boundaries:
    for ( gsMultiPatch<>::const_biterator
                bit = Psi.bBegin(); bit != Psi.bEnd(); ++bit)
    {
        bc.addCondition( *bit, condition_type::dirichlet, &Dg );
    }
    gsInfo<<"Boundary conditions:\n"<< bc <<"\n";

    //! [Refinement]
    //gsMultiBasis<> hdbasis(Psi, true);//true: poly-splines (not NURBS)

    /* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ###   Step 5: local h-refinement
     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
    gsInfo << "Patches: "<< Psi.nPatches() <<", degree: "<< dbasis.minCwiseDegree() <<"\n";

    // // Set the geometry optimal map
    geometryMap PP    = ev.getMap(Psi);  

    // Recover manufactured solution for Poisson equation
    auto u_ex         = ev.getVariable(s, PP);
    // Recover rhs for Poisson equation
    auto SFunc        = ev.getVariable(SourceFunc, PP);
    // Coeffs for advection-reaction diffusion equation
    auto coeff_convPP = ev.getVariable(coeff_conv, PP);
    auto coeff_diffPP = ev.getVariable(coeff_diffMax, PP);
    auto coeff_reacPP = ev.getVariable(coeff_reac, PP);
    // numerical solutionas Vector
    gsMatrix<> rsolVector;
    for (int r=0; r<=numRefine; ++r)
    {
        // --------------- define Pde ---------------
        timer.restart();
        //! [definePde]
        gsConvDiffRePde<real_t> cdrPde(Psi, bc, & coeff_diff,& coeff_conv, & coeff_reac, & SourceFunc);
        //! [definePde]
        //! [constructAssembler]
        // Construct assembler
        gsCDRAssembler<real_t> cdrAss( cdrPde, dbasis);
        // Set stabilization flag to 1 = SUPG
        cdrAss.options().setInt("Stabilization", Stabilizationtype);
        // Compute Dirichlet values by L2-projection
        // Caution: Interpolation does not work for locally refined (T)HB-splines!
        cdrAss.options().setInt("DirichletValues",dirichlet::l2Projection);
        //! [constructAssembler]
        gsInfo<< "." <<std::flush;// Assemblying done
        ma_time += timer.stop();
        // Generate system matrix and load vector
        cdrAss.assemble();
        // Solution vector and solution variable
        timer.restart();
        // Solve the system
        rsolVector = gsSparseSolver<>::BiCGSTABILUT( cdrAss.matrix() ).solve( cdrAss.rhs() );

        slv_time += timer.stop();
        gsInfo<< "DoFs in PDEs " << cdrAss.numDofs() <<std::flush;
        DoFPDE[r] = cdrAss.numDofs();
        gsField<> solField;        // Construct the solution as a scalar field
        solField = cdrAss.constructSolution(rsolVector);

        ev.setIntegrationElements(dbasis);

        gsExprEvaluator<>::variable ru_sol = ev.getVariable(solField.fields());

        // omp_set_dynamic(0);     // Explicitly disable dynamic teams
        // omp_set_num_threads(1); // Use these threads for later parallel regions

        //! [errorComputation]
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
        // for test :ev.integralElWise( ( ilapl(ru_sol, PP) + SFunc ).sqNorm()*meas(PP) );
        // Get the element-wise norms.
        ev.integralElWise( ( coeff_diffPP * ilapl(ru_sol,PP) - igrad(ru_sol, PP)*coeff_convPP - coeff_reacPP * ru_sol + SFunc).sqNorm()*meas(PP) );
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

       // Call repair interfaces to make sure that the new meshes
       // match along patch interfaces.
       dbasis.repairInterfaces( Psi.interfaces() );
       //! [refreshAssembler]
       cdrAss.refresh();
        }
    if(plot && r == numRefine){
        // gsInfo<<"Storing paraview...\n";
        // // Write the computed solution to paraview files
        // gsWriteParaview<>( solField, "adaptRef", 1000, true);
        gsInfo<<"Making in Paraview...\n";
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
        collection.addField(ff,"Density function");
        collection.saveTimeStep();
        collection.save();
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
       //! [Export to Paraview]
       // Export the final solution
    if(plot){
        //------------------------------------
        gsInfo<<"Plotting in Paraview...\n";
        // Run paraview
        gsFileManager::open("ParaviewOutput/solution.pvd");
    }
    else
        gsInfo << "Done. No output created, re-run with --plot to get a ParaView "
                  "file containing the solution.\n";
    //! [Export visualization in ParaView]

    return EXIT_SUCCESS;


}// end main