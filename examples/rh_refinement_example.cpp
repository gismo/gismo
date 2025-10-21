/** @file rh_refinement_example.cpp

    @brief Tutorial on how to use expression assembler to solve a Ellepitc PDE with adaptive refinement

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris & M. BAHARI
*/

//! [Include namespace]
#include <gismo.h>
#include <fstream> // For file operations
#include <gsAssembler/gsAdaptiveMultiPatchBuilder.h>  // Include the new class of r_refinement

int main(int argc, char *argv[])
{
    //! [Parse command line]
    bool plot             = false;
    index_t numRefine     = 4;
    index_t numLRefine    = 3;
    index_t numElevate    = 0;
    index_t maxIter       = 50;
    index_t NumArMarEl    = 0; // Number of ring of cells around marked elements
    double IntensityMAE   = 12.;
    bool export_b64       = false;
    bool errorsave        = false;
    bool basisRation      = true; //true: poly-splines (not NURBS)
    // --------------- adaptive refinement ---------------
    // Specify cell-marking strategy... 
    index_t adaptRefCrit  = 2;  // 1: GARU, 2: PUCA, 3: BULK, 4: PBULK
    real_t  adaptRefParam = 0.; // ... adapt parameter.
    // Specify the file path
    // std::string fn("pde/quart_annulus.xml");
    std::string fn("pde/circle.xml");
    // std::string fn("pde/lshape.xml");
    // std::string fn("domain2d/lake.xml");
    
    gsCmdLine cmd("Tutorial on solving a non-linear Monge-Ampere problem.");
    cmd.addReal( "a",    "adaptRefParam",    "parameter for local h-refinement loops",                                  adaptRefParam );
    cmd.addInt( "c",     "NumArMarEl",       "augement NumArMarEl with such quantity in local h-refinement loops",      NumArMarEl );
    cmd.addString( "d",  "file",             "Input XML file data",                                                     fn );
    cmd.addInt( "e",     "degreeElevation",  "Number of degree elevation steps to perform before solving",              numElevate );
    cmd.addReal( "f",    "IntensityMAE",     "Intensity of density function",                                           IntensityMAE);
    cmd.addInt("i",      "iter",             "Maximum number of iterations for the iterative Picard",                   maxIter);
    cmd.addInt( "l",     "numLRefine",       "Number of local h-refinement loops",                                      numLRefine );
    cmd.addInt( "r",     "adaptRefCrit",     "Adaptive refinement criterion [1:GARU,2:PUCA,3:BULK,4:PBULK]",            adaptRefCrit );
    cmd.addInt( "u",     "uniformRefine",    "Number of Uniform h-refinement loops",                                    numRefine );

    cmd.addInt("quRule",                     "Quadrature rule [1:GaussLegendre,2:GaussLobatto,3:PatchRule]",            1);
    cmd.addSwitch("plot",                    "Create a ParaView visualization file with the solution",                  plot);
    cmd.addSwitch("errorsave",               "Create a file in ... and save errors",                                    errorsave);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    // Load the file
    gsFileData<> fd(fn);
    gsInfo << "Loaded file " << fd.lastPath() << "\n";
    // Create a gsMultipatch and add the loaded geometry
    gsMultiPatch<> mpLeft;// = gsNurbsCreator<>::BSplineSquareGrid(1,1,1, 0.0, 0.0);
    fd.getId(1,mpLeft);
    // Elevate and p-refine the basis to order p + numElevate
    // where p is the highest degree in the bases
    mpLeft.degreeElevate(numElevate);
    mpLeft.computeTopology();

    //..... Test 1
    // convection coefficient:
    gsMatrix<> coeff_conv{1,2};
    // diffusion coefficient:
    //double coeff_diff = 1.;
    // reaction coefficient:
    //double coeff_reac = 0.;

    // source term: and manufactured solution
    gsFunctionExpr<> s;
    fd.getId(2000, s);
    gsFunctionExpr<> rhs;
    fd.getId(2001, rhs);

    //! [Refinement]
    gsMultiBasis<double> dbasis(mpLeft, basisRation);
    gsInfo << "Patches: "<< mpLeft.nPatches() <<", degree: "<< dbasis.minCwiseDegree() <<"\n";
#ifdef _OPENMP
    gsInfo<< "Available threads: "<< omp_get_max_threads() <<"\n";
#endif
    //! [Problem setup]
    gsExprAssembler<> A(1,1);
    //A.setOptions(Aopt);
    gsInfo<<"Active options:\n"<< A.options() <<"\n";

    A.setIntegrationElements(dbasis);
    gsExprEvaluator<> ev(A);

    typedef gsExprAssembler<>::geometryMap geometryMap;
    typedef gsExprAssembler<>::variable    variable;
    typedef gsExprAssembler<>::space       space;
    typedef gsExprAssembler<>::solution    solution;


    //! [Solver loop]
    gsSparseSolver<>::CGDiagonal solver;

    for (int r=0; r<=numRefine; ++r)
    {
        dbasis.uniformRefine();
        mpLeft.uniformRefine();
    }
    gsBoundaryConditions<> bc;
    bc.setGeoMap(mpLeft);
    // For simplicity, set Dirichlet boundary conditions
    for ( gsMultiPatch<>::const_biterator
            bit = mpLeft.bBegin(); bit != mpLeft.bEnd(); ++bit)
    {
       bc.addCondition( *bit, condition_type::dirichlet, &s,0, false);
    }
    geometryMap GLeft = A.getMap(mpLeft);
    gsStopwatch timer;

    // Set the discretization space // different boundary condition !
    space ru   = A.getSpace(dbasis);

    // Set the source term for Poisson equation
    auto rhs_f = A.getCoeff(rhs, GLeft);

    // Solution vector and solution variable
    gsMatrix<> rsolVector;
    solution u_sol = A.getSolution(ru, rsolVector);

    /* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ###   Step 0: Computes the initial solution of the PDEs 
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
    ru.setup(bc, dirichlet::l2Projection, 0);
    // Initialize the system
    A.initSystem();
    gsInfo<< "Solving PDEs " <<std::flush;
    gsInfo<< A.numDofs() <<std::flush;    
    //auto h_Tau =  m_h/(2.*coeff_conv.squaredNorm()+m_h);
    timer.restart();
    A.assemble(
    igrad(ru, GLeft) * igrad(ru, GLeft).tr() * meas(GLeft) //matrix
    ,
    ru * rhs_f * meas(GLeft) //rhs vector
    );
    gsInfo<< "." <<std::flush;// Assemblying done
    timer.restart();
    solver.compute( A.matrix() );
    rsolVector     = solver.solve(A.rhs());
    ev.integralElWise( (ilapl(u_sol, GLeft) +rhs_f).sqNorm()*meas(GLeft) );
    //ev.integralElWise( igrad(u_sol, GLeft).sqNorm() );
    auto elwise    = ev.elementwise();


    /* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ###   Step 1-2 : Computes the density function
    ###         and the multipatch adaptive mapping
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
    gsAdaptiveMultiPatchBuilder MAE = gsAdaptiveMultiPatchBuilder(dbasis, mpLeft, maxIter, IntensityMAE);
    gsMultiPatch<> Psitp;
    if (IntensityMAE <= 1.)
    {
        gsInfo << "IntensityMAE < 1, so mpLeft is used.\n";
        // Use the original multipatch
        Psitp = mpLeft;
    }  
    else{
    //.. mark elements location
    std::vector<bool> eldensityMarked( elwise.size() );
    gsMarkElementsForRef( elwise, adaptRefCrit, 0.8, eldensityMarked);
    auto density   = MAE.buildDensity(dbasis, eldensityMarked);
    MAE.buildMultiPatch(density);// compute adaptive mapping
    Psitp     = MAE.buildCompMultiPatch(dbasis);// computes the composition mapping mpLeft o Psitp    
    }

    /* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ###   Step 3: Define hierarchical adaptive mapping
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
    gsMultiPatch<> Psi;
    if (Psitp.basis(0).weights().any()){
    for(size_t i =0; i<Psitp.nPatches(); ++i)
        Psi.addPatch(gsRationalTHBSpline<2>( dynamic_cast<const gsTensorNurbs<2>&>(Psitp.patch(i)) ));
    }
    else{
    for(size_t i =0; i<Psitp.nPatches(); ++i)
        Psi.addPatch(gsTHBSpline<2>( dynamic_cast<const gsTensorBSpline<2>&>(Psitp.patch(i)) ));
    }
    Psi.addAutoBoundaries();
    Psi.computeTopology();

    /* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ###   Step 3: Start hierarchical refinement
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
    //::::::::::::::::::::   Poisson equation - (manufactured exact solution)         :::::::::::::::::::::::::
    if (true){
    gsMultiPatch<> sol_restr; // restricted solution
    //boubdary conditions
    gsBoundaryConditions<> bc;
    bc.setGeoMap(Psi);
    // For simplicity, set Dirichlet boundary conditions
    for ( gsMultiPatch<>::const_biterator
            bit = Psi.bBegin(); bit != Psi.bEnd(); ++bit)
    {
       bc.addCondition( *bit, condition_type::dirichlet, &s,0, false);
    }

    // gsFunctionExpr<> g("0.*x","0.*y",2);
    // bc.addCondition(0,2, condition_type::neumann, &g,0,false);
    // bc.addCondition(0,4, condition_type::dirichlet, &s,0,false);
    // bc.addCondition(0,3, condition_type::dirichlet, &s,0,false);
    // bc.addCondition(0,1, condition_type::neumann, &g,0,false);
    gsMultiBasis<> dbasis(Psi, true);//true: poly-splines (not NURBS)

    gsInfo<<"The PDE domain is "<< Psi.detail() << "\n";
    gsInfo<<"Source function is "<< rhs << "\n";
    gsInfo<<"Boundary conditions:\n"<< bc <<"\n";
    //! [Problem setup]
    gsExprAssembler<> A(1,1);
    A.setIntegrationElements(dbasis);
    gsExprEvaluator<> ev(A);

    geometryMap PP  = A.getMap(Psi);
    
    gsStopwatch timer;
    // Elements used for numerical integration
    A.setIntegrationElements(dbasis);
    
    // Set the discretization space // different boundary condition !
    space ru        = A.getSpace(dbasis);

    // Set the source term for Poisson equation
    auto rhs_f      = A.getCoeff(rhs, PP);

    // Recover manufactured solution for Poisson equation
    auto u_ex       = ev.getVariable(s, PP);

    // Solution vector and solution variable
    gsMatrix<> rsolVector;
    solution ru_sol = A.getSolution(ru, rsolVector);

    gsVector<>  h1err(numLRefine+1), l2err(numLRefine+1);
    gsVector<int>  DoFPDE(numLRefine+1);

    gsInfo<< "(dot1=assembled, dot2=solved)\n";
    double setup_time(0), ma_time(0), slv_time(0), err_time(0);
    for (int r=0; r<=numLRefine; ++r)
    {
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
        ru * rhs_f * meas(PP) //rhs vector
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

        // Compute the global error indicators.
        DoFPDE[r] = A.numDofs();
        timer.restart();
        l2err[r]= math::sqrt( ev.integral( (u_ex - ru_sol).sqNorm() * meas(PP) ) );
        h1err[r]= math::sqrt(ev.integral( ( grad(u_ex) - igrad(ru_sol,PP) ).sqNorm() * meas(PP) ));
        err_time += timer.stop();
        gsInfo<< ". " <<std::flush; // Error computations done
        if(r < numLRefine){
        //! [beginRefLoop]
            gsInfo << "====== Loop " << r << " of "
                    <<numLRefine<< " ====adapt Parameter ="<< adaptRefParam << " ======" << "\n";

            // --------------- error estimation/computation ---------------
            // Get the element-wise norms.
            // ev.integralElWise( (  igrad(ru_sol, PP) ).sqNorm() );
            ev.integralElWise( (  ilapl(ru_sol, PP) + rhs_f ).sqNorm() );
            //! [errorComputation]
            std::vector<real_t> eltErrs  = ev.elementwise();            
            if (IntensityMAE >1.){
                std::vector<bool> eldensityMarked( eltErrs.size() );
                gsMarkElementsForRef( eltErrs, adaptRefCrit, 0.8, eldensityMarked);                 
                auto density   = MAE.buildDensity( dbasis, eldensityMarked);
                MAE.buildMultiPatch(density);// compute adaptive mapping
                Psi            = MAE.buildCompMultiPatch(dbasis);// computes the composition mapping mpLeft o Psitp
            }

            // /* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            std::vector<bool> elMarked( eltErrs.size() );
            // //! [adaptRefinementPart]
            // Mark elements for refinement, based on the computed local errors and
            // the refinement-criterion and -parameter.
            gsMarkElementsForRef( eltErrs, adaptRefCrit, adaptRefParam, elMarked);
            if (IntensityMAE >1.){
                // -----------------
                double Minvalue     = *std::max_element(eltErrs.begin(), eltErrs.end());                
                for(size_t i=0; i<eltErrs.size(); ++i)
                {
                    if (elMarked[i] == true and Minvalue > eltErrs[i]) // Avoid numerical issues
                    {
                        Minvalue = eltErrs[i];
                    }
                }
                for(size_t i=0; i<eltErrs.size(); ++i)
                {
                    auto POWMinvalue = Minvalue/pow(DoFPDE[r],adaptRefParam);
                    if (POWMinvalue < eltErrs[i]) // Avoid numerical issues
                    {
                        elMarked[i] = true;
                    }
                }
            }
            gsInfo <<"Marked "<< std::count(elMarked.begin(), elMarked.end(), true) <<" elements.\n";
            // Refine the marked elements with a 1-ring of cells around marked elements
            gsRefineMarkedElements( dbasis, elMarked, NumArMarEl);
            gsRefineMarkedElements( Psi, elMarked, NumArMarEl);
            //
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

    if (errorsave)
    {
    // Assuming DoFPDE, l2err, and h1err are gsMatrix or similar types
    std::ofstream outFile("error_analysis.txt", std::ios::app); // Open file in append mode
    if (outFile.is_open())
    {
        outFile << "#DoF_PDE: " << adaptRefParam <<" "<< NumArMarEl <<" " << IntensityMAE << " \n"<< std::scientific << DoFPDE.transpose() << "\n";
        outFile << "#L2_error: \n" << std::scientific << std::setprecision(3) << l2err.transpose() << "\n";
        outFile << "#H1_error: \n" << std::scientific << std::setprecision(3) << h1err.transpose() << "\n";
        outFile << "#-------------------------------------------------------------------------------\n"; // Optional separator for readability
        outFile.close(); // Close the file after writing
    }
    else
    {
        gsInfo << "Error: Unable to open file for writing : error_analysis.txt.\n";
    }
    }
    else
    {
        gsInfo << "Errors are not saved. To save them, try with --errorsave.\n";
    }

    //! [Error and convergence rates]
    if (numLRefine>0 && false)
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
        collection.addField((  ilapl(ru_sol, PP)+ rhs_f ).sqNorm()*meas(PP),"indecator");
        collection.addField(jac(PP).det(), "Jacobian function");
        collection.addField(u_ex, "exact solution");
        collection.saveTimeStep();
        collection.save();
        gsFileManager::open("ParaviewOutput/solution.pvd");
    }
    else
        gsInfo << "Done. No output created, re-run with --plot to get a ParaView "
                  "file containing the solution.\n";
    //! [Export visualization in ParaView]
    }

    return EXIT_SUCCESS;


}// end main