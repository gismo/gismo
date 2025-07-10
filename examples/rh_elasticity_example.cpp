/** @file rh_elasticity_example.cpp

    @brief Tutorial on how to use expression assembler to solve a linear elasticity equation

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris & M. BAHARI
*/

//! [Include namespace]
#include <gismo.h>
#include <gsElasticity/src/gsElasticityAssembler.h>
#include <gsElasticity/src/gsWriteParaviewMultiPhysics.h>
#include <gsAssembler/gsAdaptiveMultiPatchBuilder.h>  // Include the new class of r_refinement

using namespace gismo;
//! [Include namespace]
// Method to Project normal control points
void CorrecNormalCPoints(const gsMultiPatch<>& mpLeft, gsMultiPatch<>& Psi, index_t boxMaxNumber = 1)
{
    // Projection normal of control points (exact geometry)
    for (int boxNumber = 0; boxNumber < boxMaxNumber; ++boxNumber)
    {
        // test if the boundary interface is not an inner interface between patches
        for (int i_x =0; i_x < Psi.patch(boxNumber).basis().boundary(1).size(); ++i_x) // x=0 control points be like (0,:) in this case
        {
            Psi.patch(boxNumber).coef( Psi.patch(boxNumber).basis().boundary(1).at(i_x) ).array()[1] = mpLeft.patch(boxNumber).coef( mpLeft.patch(boxNumber).basis().boundary(1).at(i_x) ).array()[1];
        }

        for (int i_x =0; i_x < Psi.patch(boxNumber).basis().boundary(2).size(); ++i_x)// x=1 control points be like (1,:) in this case
        {
        Psi.patch(boxNumber).coef( Psi.patch(boxNumber).basis().boundary(2).at(i_x) ).array()[0] = mpLeft.patch(boxNumber).coef( mpLeft.patch(boxNumber).basis().boundary(2).at(i_x) ).array()[0];
        }

        for (int i_x =0; i_x < Psi.patch(boxNumber).basis().boundary(3).size(); ++i_x) // y=0 control points be like (:,0) in this case
        {
        Psi.patch(boxNumber).coef( Psi.patch(boxNumber).basis().boundary(3).at(i_x) ).array()[1] = mpLeft.patch(boxNumber).coef( mpLeft.patch(boxNumber).basis().boundary(3).at(i_x) ).array()[1];
        Psi.patch(boxNumber).coef( Psi.patch(boxNumber).basis().boundary(3).at(i_x) ).array()[0] = mpLeft.patch(boxNumber).coef( mpLeft.patch(boxNumber).basis().boundary(3).at(i_x) ).array()[0];
        }
        // x= 4 y = -4
        for (int i_x =0; i_x < Psi.patch(boxNumber).basis().boundary(4).size(); ++i_x)// y=1 control points be like (:,1) in this case
        {
        Psi.patch(boxNumber).coef( Psi.patch(boxNumber).basis().boundary(4).at(i_x) ).array()[1] = mpLeft.patch(boxNumber).coef( mpLeft.patch(boxNumber).basis().boundary(4).at(i_x) ).array()[1];
        Psi.patch(boxNumber).coef( Psi.patch(boxNumber).basis().boundary(4).at(i_x) ).array()[0] = mpLeft.patch(boxNumber).coef( mpLeft.patch(boxNumber).basis().boundary(4).at(i_x) ).array()[0];
        }
    }
}

int main(int argc, char *argv[])
{
    //! [Parse command line]
    bool plot             = false;
    index_t numRefine     = 2;
    index_t numLRefine    = 0;
    index_t numElevate    = 0;
    index_t maxIter       = 30;
    index_t NumArMarEl    = 0; // Number of ring of cells around marked elements
    double IntensityMAE   = 12.;
    //bool export_b64     = false;
    bool errorsave        = false;
    real_t adaptRefParam  = 0.;     // ... adapt parameter.
    index_t FactRefPar    = 0;    // ... adapt parameter : adaptRefParam += FactRefPar in each iter
    // Specify the file path
    std::string fn("pde/infinit_plate.xml");
    // --------------- adaptive refinement ---------------
    index_t adaptRefCrit  = 2;  // 1: GARU, 2: PUCA, 3: BULK, 4: PBULK

    gsCmdLine cmd("Tutorial on solving a non-linear Monge-Ampere problem.");
    cmd.addInt("i", "iter", "Maximum number of iterations for the iterative Picard", maxIter);
    cmd.addInt( "e", "degreeElevation",
                "Number of degree elevation steps to perform before solving (0: equalize degree in all directions)", numElevate );
    cmd.addInt( "u", "uniformRefine", "Number of Uniform h-refinement loops",  numRefine );
    cmd.addInt( "l", "numLRefine", "Number of local h-refinement loops",  numLRefine );
    cmd.addString( "d", "file", "Input XML file data", fn );
    cmd.addInt("quRule",
                 "Quadrature rule [1:GaussLegendre,2:GaussLobatto,3:PatchRule]",
                 1);
    cmd.addReal( "a", "adaptRefParam", "parameter for local h-refinement loops",  adaptRefParam );
    cmd.addInt( "p", "FactRefPar", "augement adaptRefParam with such quantity in local h-refinement loops",  FactRefPar );
    cmd.addSwitch("plot", "Create a ParaView visualization file with the solution", plot);
    cmd.addSwitch("errorsave", "Create a file in ... and save errors", errorsave);
    cmd.addReal( "f", "IntensityMAE", "Intensity of density function",  IntensityMAE);
    cmd.addInt( "c", "NumArMarEl", "augement NumArMarEl with such quantity in local h-refinement loops",  NumArMarEl );
    cmd.addInt( "r", "adaptRefCrit", "Adaptive refinement criterion [1:GARU,2:PUCA,3:BULK,4:PBULK]",  adaptRefCrit );

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    // Load the file
    gsFileData<> fd(fn);
    gsInfo << "Loaded file " << fd.lastPath() << "\n";
    // Create a gsMultipatch and add the loaded geometry
    gsMultiPatch<> mpLeft;
    fd.getId(1,mpLeft);
    // Elevate and p-refine the basis to order p + numElevate
    // where p is the highest degree in the bases
    mpLeft.degreeElevate(numElevate);
    mpLeft.computeTopology();

    // Right-hand side function : Analytical density function (det(H(u))=f= sigma/rho)
    // Load the file
    gsFunctionExpr<> f;
    fd.getId(2003, f);
    gsFunctionExpr<> rhs;
    fd.getId(2001, rhs);
    gsInfo<<"Density function "<< f << "\n";

    //=============================================//
        // Setting loads and boundary conditions //
    //=============================================//

    gsFunctionExpr<> analyticalStresses("1-1/(x^2+y^2)*(3/2*cos(2*atan2(y,x)) + cos(4*atan2(y,x))) + 3/2/(x^2+y^2)^2*cos(4*atan2(y,x))",
                                        "-1/(x^2+y^2)*(1/2*cos(2*atan2(y,x)) - cos(4*atan2(y,x))) - 3/2/(x^2+y^2)^2*cos(4*atan2(y,x))",
                                        "-1/(x^2+y^2)*(1/2*sin(2*atan2(y,x)) + sin(4*atan2(y,x))) + 3/2/(x^2+y^2)^2*sin(4*atan2(y,x))",2);
    // boundary load neumann BC
    gsFunctionExpr<> traction("(1-1/(x^2+y^2)*(3/2*cos(2*atan2(y,x)) + cos(4*atan2(y,x))) + 3/2/(x^2+y^2)^2*cos(4*atan2(y,x))) * (x/sqrt(x^2+y^2)) +"
                              "(-1/(x^2+y^2)*(1/2*sin(2*atan2(y,x)) + sin(4*atan2(y,x))) + 3/2/(x^2+y^2)^2*sin(4*atan2(y,x))) * (y/sqrt(x^2+y^2))",
                              "(-1/(x^2+y^2)*(1/2*sin(2*atan2(y,x)) + sin(4*atan2(y,x))) + 3/2/(x^2+y^2)^2*sin(4*atan2(y,x))) * (x/sqrt(x^2+y^2)) +"
                              "(-1/(x^2+y^2)*(1/2*cos(2*atan2(y,x)) - cos(4*atan2(y,x))) - 3/2/(x^2+y^2)^2*cos(4*atan2(y,x))) * (y/sqrt(x^2+y^2))",2);
    // material parameters
    real_t youngsModulus = 1.0e3;
    real_t poissonsRatio = 0.3;

    // boundary conditions
    gsBoundaryConditions<> bcInfo;
    bcInfo.addCondition(0,boundary::north,condition_type::neumann,&traction);
    bcInfo.addCondition(0,boundary::west,condition_type::dirichlet,nullptr,1); // last number is a component (coordinate) number
    bcInfo.addCondition(0,boundary::east,condition_type::dirichlet,nullptr,0);

    // source function, rhs
    gsConstantFunction<> g(0.,0.,2);

    //! [Refinement]
    gsMultiBasis<> dbasis(mpLeft, true);//true: poly-splines (not NURBS)
    
    gsInfo << "Patches: "<< mpLeft.nPatches() <<", degree: "<< dbasis.minCwiseDegree() <<"\n";
#ifdef _OPENMP
    gsInfo<< "Available threads: "<< omp_get_max_threads() <<"\n";
#endif
    //::::::::::::::::::::      mesh adaptation solver         :::::::::::::::::::::::::
    for (int r=0; r<=numRefine; ++r)
    {
        dbasis.uniformRefine();
        mpLeft.uniformRefine();
    }

    //=================================================//
         // Assembling & solving for the initial guess//
    //=================================================//

    // creating assembler
    gsElasticityAssembler<real_t> assembler(mpLeft,dbasis,bcInfo,g);
    assembler.options().setReal("YoungsModulus",youngsModulus);
    assembler.options().setReal("PoissonsRatio",poissonsRatio);
    gsInfo<<"Assembling...\n";
    gsStopwatch clock;
    clock.restart();
    assembler.assemble();
    gsInfo << "Assembled a system (matrix and load vector) with "
           << assembler.numDofs() << " dofs in " << clock.stop() << "s.\n";

    gsInfo << "Solving...\n";
    clock.restart();

#ifdef GISMO_WITH_PARDISO
    gsSparseSolver<>::PardisoLLT solver(assembler.matrix());
    gsVector<> solVector = solver.solve(assembler.rhs());
    gsInfo << "Solved the system with PardisoLDLT solver in " << clock.stop() <<"s.\n";
#else
    gsSparseSolver<>::SimplicialLDLT solver(assembler.matrix());
    gsVector<> solVector = solver.solve(assembler.rhs());
    gsInfo << "Solved the system with EigenLDLT solver in " << clock.stop() <<"s.\n";
#endif

    //=============================================//
                      // Output //
    //=============================================//

    // constructing displacement as an IGA function
    gsMultiPatch<> solution;
    assembler.constructSolution(solVector,assembler.allFixedDofs(),solution);
    // constructing stress tensor
    gsPiecewiseFunction<> stresses;
    assembler.constructCauchyStresses(solution,stresses,stress_components::all_2D_vector);

    // constructing an IGA field (geometry + solution) for displacement
    gsField<> solutionField(assembler.patches(),solution);
    // constructing an IGA field (geometry + solution) for stresses
    gsField<> stressField(assembler.patches(),stresses,true);

    gsExprEvaluator<> evi;
    evi.setIntegrationElements(assembler.multiBasis());
    gsExprEvaluator<>::geometryMap PPi = evi.getMap(mpLeft);
    //... error computation
    auto istress = evi.getVariable(stressField.fields());

    // ev.integralElWise( idiv(istress, PP) * meas(PP) )
    gsInfo << "Stress: min "<< istress.ppart() << "\n";
    evi.integralElWise( idiv(istress, PPi).sqNorm());
    // ev.integralBdrBc(bcInfo.get("nuemann"),(istress* nv(PP)).sqNorm());
    auto elwise = evi.elementwise();

    /* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ###   Step 1-2 : Computes the density function
    ###        and the multipatch adaptove mapping
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
    gsAdaptiveMultiPatchBuilder MAE = gsAdaptiveMultiPatchBuilder(dbasis, mpLeft, numElevate, maxIter, IntensityMAE);
    auto density                    = MAE.buildStrategyDensity( elwise, 0.8);

    if (IntensityMAE>0){
    // auto density                 = MAE.buildAnalyticDensity( f);
    auto geometrytp                 = MAE.buildMultiPatch(density);
    CorrecNormalCPoints(mpLeft, geometrytp);
    index_t numPaches               = geometrytp.nPatches();
    for( index_t i=0; i<numPaches; ++i)
    {
    index_t coefsNum                = geometrytp.patch(i).coefsSize();
    for ( index_t j=0; j<coefsNum; ++j)
    {
        mpLeft.patch(i).coef(j)     = geometrytp.patch(i).coef(j);
    }
    }
    }
    // mpLeft.coefs().swap(geometrytp.coefs());
    // gsInfo<<"Making geometry"<< (mpLeft.coefs()-geometrytp.coefs())<< "\n";
    /* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ###   Step 4: Define hierarchical adaptive mapping
     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
    gsMultiPatch<> geometry;	
    for(size_t i =0; i<mpLeft.nPatches(); ++i)
        geometry.addPatch(gsRationalTHBSpline<2>( dynamic_cast<const gsTensorNurbs<2>&>(mpLeft.patch(i)) ));
        //geometry.addPatch(gsTHBSpline<2>( dynamic_cast<const gsTensorBSpline<2>&>(geometrytp.patch(i)) ));
    geometry.addAutoBoundaries();
    geometry.computeTopology();
    //gsWrite(geometry, "geometry_mapping");
    //#-++++++++++++++++++++++++ End of sharing part of any geometry------------------------------

    //::::::::::::::::::::  Elasticity equation - (manufactured exact solution)         :::::::::::::::::::::::::
    if (true){
    // creating basis
    gsMultiBasis<> basis(geometry, true);//true: poly-splines (not NURBS)
    gsVector<>     l2err(numLRefine+1), h1err(numLRefine+1);
    gsVector<int>  DoFPDE(numLRefine+1);

    //=============================================//
              // Assembling & solving //
    //=============================================//

    // creating assembler
    gsElasticityAssembler<real_t> assembler(geometry,basis,bcInfo,g);
    assembler.options().setReal("YoungsModulus",youngsModulus);
    assembler.options().setReal("PoissonsRatio",poissonsRatio);
    
    for (int r=0; r<=numLRefine; ++r)
    {

        //=============================================//
                // Assembling & solving //
        //=============================================//
        // creating assembler
        gsElasticityAssembler<real_t> assembler(geometry,basis,bcInfo,g);
        assembler.options().setReal("YoungsModulus",youngsModulus);
        assembler.options().setReal("PoissonsRatio",poissonsRatio);
        gsInfo<<"Assembling...\n";
        gsStopwatch clock;
        clock.restart();
        assembler.assemble();
        gsInfo << "Assembled a system (matrix and load vector) with "
            << assembler.numDofs() << " dofs in " << clock.stop() << "s.\n";

        gsInfo << "Solving...\n";
        clock.restart();

    #ifdef GISMO_WITH_PARDISO
        gsSparseSolver<>::PardisoLLT solver(assembler.matrix());
        gsVector<> solVector = solver.solve(assembler.rhs());
        gsInfo << "Solved the system with PardisoLDLT solver in " << clock.stop() <<"s.\n";
    #else
        gsSparseSolver<>::SimplicialLDLT solver(assembler.matrix());
        gsVector<> solVector = solver.solve(assembler.rhs());
        gsInfo << "Solved the system with EigenLDLT solver in " << clock.stop() <<"s.\n";
    #endif

        //=============================================//
                        // Output //
        //=============================================//
        // // constructing displacement as an IGA function
        gsMultiPatch<> solution;
        assembler.constructSolution(solVector,assembler.allFixedDofs(),solution);
        // constructing stress tensor
        gsPiecewiseFunction<> stresses;
        assembler.constructCauchyStresses(solution,stresses,stress_components::all_2D_vector);

        // constructing an IGA field (geometry + solution) for displacement
        gsField<> solutionField(assembler.patches(),solution);
        // constructing an IGA field (geometry + solution) for stresses
        gsField<> stressField(assembler.patches(),stresses,true);
        // analytical stresses
        gsField<> analyticalStressField(assembler.patches(),analyticalStresses,false);

        gsExprEvaluator<> ev;
        ev.setIntegrationElements(assembler.multiBasis());
        gsExprEvaluator<>::geometryMap PP = ev.getMap(geometry);
        auto sigm_ex      = ev.getVariable(analyticalStresses.piece(0), PP);
        //... error computation
        auto istress      = ev.getVariable(stressField.fields().piece(0));
        // eval stress at the top of the circular cut
        gsMatrix<> A(2,1);
        A << 1.,0.; // parametric coordinates for the isogeometric solution
        gsMatrix<> res;
        stresses.piece(0).eval_into(A,res);
        A << 0., 1.; // spatial coordinates for the analytical solution
        gsMatrix<> analytical;
        analyticalStresses.eval_into(A,analytical);
        //...
        DoFPDE[r]         = assembler.numDofs();
        l2err[r]          = math::sqrt( ev.integral( ( sigm_ex - istress).sqNorm() * meas(PP) ));
        h1err[r]          = math::sqrt(pow(res.at(0)-analytical.at(0),2));
        //
        gsInfo << "err: " << math::sqrt( ev.integral( ( sigm_ex - istress).sqNorm() * meas(PP) )) << " (errc), " << math::sqrt(pow(res.at(0)-analytical.at(0),2)) << " \n";
        gsInfo << "XX-stress at the top of the circle: " << res.at(0) << " (computed), " << analytical.at(0) << " (analytical)\n";
        gsInfo << "YY-stress at the top of the circle: " << res.at(1) << " (computed), " << analytical.at(1) << " (analytical)\n";
        gsInfo << "XY-stress at the top of the circle: " << res.at(2) << " (computed), " << analytical.at(2) << " (analytical)\n";
        if(r < numLRefine){
        //! [beginRefLoop]
            gsInfo << "====== Loop " << r << " of "
                    <<numLRefine<< " ====adapt Parameter ="<< adaptRefParam << " ======" << "\n";
            // --------------- error estimation/computation ---------------
            // Compute the error indicators
            //ev.integralElWise( ff );
            ev.integralElWise( idiv(istress,PP).sqNorm());
            //ev.integralBdrBc(bcInfo.get("neumann"),(istress* nv(PP)).sqNorm());

            std::vector<real_t> eltErrs  = ev.elementwise();
            //! [adaptRefinementPart]
            // Mark elements for refinement, based on the computed local errors and
            // the refinement-criterion and -parameter.
            std::vector<bool> elMarked( eltErrs.size() );
            gsMarkElementsForRef( eltErrs, adaptRefCrit, adaptRefParam, elMarked);
            if (IntensityMAE >1.){
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
                    if (Minvalue/pow(DoFPDE[r],adaptRefParam) < eltErrs[i]) // Avoid numerical issues
                    {
                        elMarked[i] = true;
                    }
                }
            }
            gsInfo <<"Marked "<< std::count(elMarked.begin(), elMarked.end(), true) <<" elements.\n";

            // Refine the marked elements with a 1-ring of cells around marked elements
            gsRefineMarkedElements( basis, elMarked, NumArMarEl);
            gsRefineMarkedElements( assembler.multiBasis(), elMarked, NumArMarEl);
            gsRefineMarkedElements( geometry, elMarked, NumArMarEl);
            gsInfo << "assemble refined\n";
            
            NumArMarEl = NumArMarEl + FactRefPar;
            assembler.refresh();
            }
    //! [Export visualization in ParaView]
    if (plot and r ==numLRefine)
    {
        // Write the computed solution to paraview files
        gsInfo<<"Making in Paraview...\n";
        gsParaviewCollection collection("ParaviewOutput/solution", &ev);
        collection.options().setSwitch("plotElements", true);
        collection.options().setSwitch("base64", false);
        collection.options().setInt("plotElements.resolution", 16);
        collection.options().setInt("numPoints", 10000);
        collection.newTimeStep(&geometry);
        collection.addField(istress,"numerical stress");
        // collection.addField(jac(PP).det(), "Jacobian function");
        collection.addField(sigm_ex, "exact stress");
        // collection.addField(ff_Ggeometry,"Density function");
        collection.saveTimeStep();
        collection.save();
    }
    }
    //! [Solver loop]    

    //! [Error and convergence rates]
    gsInfo<< "\nDoF_PDE = "<<std::scientific<<DoFPDE.transpose()<<"\n";
    gsInfo<< "L2_error_x = "<<std::scientific<<std::setprecision(3)<<l2err.transpose()<<"\n";
    gsInfo<< "L2_error_y= "<<std::scientific<<std::setprecision(3)<<h1err.transpose()<<"\n";

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
        //------------------------------------
        gsInfo<<"Plotting in Paraview...\n";
        // Run paraview
        gsFileManager::open("ParaviewOutput/solution.pvd");
    }
    else
        gsInfo << "Done. No output created, re-run with --plot to get a ParaView "
                  "file containing the solution.\n";
    //! [Export visualization in ParaView]
    }

    return EXIT_SUCCESS;


}// end main