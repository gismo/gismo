/// This is the 2D linear elasticity benchmark "Infinite plate with circular hole"
/// as described in V.P.Nguyen, C.Anitescu, S.P.A.Bordas, T.Rabczuk, 2015
/// "Isogeometric analysis: An overview and computer implementation aspects".
///
/// Author: A.Shamanskiy (2016 - ...., TU Kaiserslautern)
#include <gismo.h>
#include <gsElasticity/gsElasticityAssembler.h>
#include <gsElasticity/gsWriteParaviewMultiPhysics.h>
#include <gsAssembler/gsAdaptiveMultiPatchBuilder.h>  // Include the new class of r_refinement

using namespace gismo;

int main(int argc, char* argv[]){

    gsInfo << "This is the 2D linear elasticity benchmark: infinite plate with circular hole.\n";

    //=====================================//
                // Input //
    //=====================================//
    std::string filename = ("pde/infinit_plate.xml");
    //! [Parse command line]
    bool plot             = true;
    index_t numRefine     = 3;
    index_t numLRefine    = 0;
    index_t numElevate    = 0;
    index_t maxIter       = 30;
    index_t NumArMarEl    = 0; // Number of ring of cells around marked elements
    double IntensityMAE   = 9.;
    //bool export_b64     = false;
    bool errorsave        = false;
    real_t adaptRefParam  = 0.;     // ... adapt parameter.
    index_t FactRefPar    = 0;    // ... adapt parameter : adaptRefParam += FactRefPar in each iter
    // Specify the file path
    std::string fn("pde/infinit_plate.xml");
    // --------------- adaptive refinement ---------------
    // Specify cell-marking strategy... 
    //MarkingStrategy adaptRefCrit = PUCA;
    //MarkingStrategy adaptRefCrit = GARU;
    //MarkingStrategy adaptRefCrit = errorFraction;
    // Elements used for numerical integration

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

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    //=============================================//
        // Scanning geometry and creating bases //
    //=============================================//
    // scanning geometry
    gsFileData<> fd(filename);
    gsInfo << "Loaded file " << fd.lastPath() << "\n";
    // Create a gsMultipatch and add the loaded geometry
    gsMultiPatch<> mpLeft;
    fd.getId(1,mpLeft);
    // Right-hand side function : Analytical density function (det(H(u))=f= sigma/rho)
    // Load the file
    gsFunctionExpr<> f;
    fd.getId(2003, f);

    // creating basis
    gsMultiBasis<> dbasis(mpLeft);
    for (index_t i = 0; i < numElevate; ++i)
        dbasis.degreeElevate();
    for (index_t i = 0; i < numRefine; ++i)
        dbasis.uniformRefine();    

    /* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ###   Step 1-2 : Computes the density function
    ###         and the multipatch adaptove mapping
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
    gsAdaptiveMultiPatchBuilder MAE = gsAdaptiveMultiPatchBuilder(dbasis, mpLeft, numElevate, maxIter, IntensityMAE);
    auto density    = MAE.buildAnalyticDensity( f);
    auto geometry   = MAE.buildMultiPatch(density);

    gsWrite(geometry, "geometry_mapping");
    /* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ###   Step 4: Define hierarchical adaptive mapping
     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
    //  gsMultiPatch<> geometry;
    //  for(size_t i =0; i<geometrytp.nPatches(); ++i)
    //      geometry.addPatch(gsTHBSpline<2>( dynamic_cast<const gsTensorBSpline<2>&>(geometrytp.patch(i)) ));
    //  geometry.addAutoBoundaries();
    //  geometry.computeTopology();

    // creating basis
    gsMultiBasis<> basis(geometry);
    // for (index_t i = 0; i < numElevate; ++i)
    //     basis.degreeElevate();
    // for (index_t i = 0; i < numRefine; ++i)
    //     basis.uniformRefine();

    gsVector<>   l2err(numLRefine+1);//h1err(numLRefine+1),
    gsVector<int>  DoFPDE(numLRefine+1);
    gsInfo<< "(dot1=assembled, dot2=solved)\n";

    //=============================================//
        // Setting loads and boundary conditions //
    //=============================================//

    gsFunctionExpr<> analyticalStresses("1-1/(x^2+y^2)*(3/2*cos(2*atan2(y,x)) + cos(4*atan2(y,x))) + 3/2/(x^2+y^2)^2*cos(4*atan2(y,x))",
                                        "-1/(x^2+y^2)*(1/2*cos(2*atan2(y,x)) - cos(4*atan2(y,x))) - 3/2/(x^2+y^2)^2*cos(4*atan2(y,x))",
                                        "-1/(x^2+y^2)*(1/2*sin(2*atan2(y,x)) + sin(4*atan2(y,x))) + 3/2/(x^2+y^2)^2*sin(4*atan2(y,x))",2);
    // boundary load neumann BC
    gsFunctionExpr<> traction("(-1+1/(x^2+y^2)*(3/2*cos(2*atan2(y,x)) + cos(4*atan2(y,x))) - 3/2/(x^2+y^2)^2*cos(4*atan2(y,x))) * (x==-4) +"
                              "(-1/(x^2+y^2)*(1/2*sin(2*atan2(y,x)) + sin(4*atan2(y,x))) + 3/2/(x^2+y^2)^2*sin(4*atan2(y,x))) * (y==4)",
                              "(1/(x^2+y^2)*(1/2*sin(2*atan2(y,x)) + sin(4*atan2(y,x))) - 3/2/(x^2+y^2)^2*sin(4*atan2(y,x))) * (x==-4) +"
                              "(-1/(x^2+y^2)*(1/2*cos(2*atan2(y,x)) - cos(4*atan2(y,x))) - 3/2/(x^2+y^2)^2*cos(4*atan2(y,x))) * (y==4)",2);
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

    // constructing displacement as an IGA function
    gsMultiPatch<> solution;
    assembler.constructSolution(solVector,assembler.allFixedDofs(),solution);
    // constructing stress tensor
    gsPiecewiseFunction<> stresses;
    assembler.constructCauchyStresses(solution,stresses,stress_components::all_2D_vector);
    // eval stress at the top of the circular cut
    gsMatrix<> A(2,1);
    A << 1.,0.; // parametric coordinates for the isogeometric solution
    gsMatrix<> res;
    stresses.piece(0).eval_into(A,res);
    A << 0., 1.; // spatial coordinates for the analytical solution
    gsMatrix<> analytical;
    analyticalStresses.eval_into(A,analytical);
    gsInfo << "XX-stress at the top of the circle: " << res.at(0) << " (computed), " << analytical.at(0) << " (analytical)\n";
    gsInfo << "YY-stress at the top of the circle: " << res.at(1) << " (computed), " << analytical.at(1) << " (analytical)\n";
    gsInfo << "XY-stress at the top of the circle: " << res.at(2) << " (computed), " << analytical.at(2) << " (analytical)\n";

    gsExprEvaluator<> ev;
    ev.setIntegrationElements(assembler.multiBasis());
    gsExprEvaluator<>::geometryMap PP = ev.getMap(geometry);
    auto sigm_ex   = ev.getVariable(analyticalStresses, PP);
    auto ff   = ev.getVariable(f, PP);

    // constructing an IGA field (geometry + solution) for displacement
    gsField<> solutionField(assembler.patches(),solution);
    // constructing an IGA field (geometry + solution) for stresses
    gsField<> stressField(assembler.patches(),stresses,true);
    // analytical stresses
    gsField<> analyticalStressField(assembler.patches(),analyticalStresses,false);
    //... error computation
    auto istress = ev.getVariable(stressField.fields());

    // omp_set_num_threads(1); // Use these threads for later parallel regions
    DoFPDE[0] = assembler.numDofs();
    l2err[0]  = math::sqrt( ev.integral( ( sigm_ex - istress).sqNorm() * meas(PP) ));
    gsInfo << " min Jacobian function and maximum " << ev.min(jac(PP).det())<< " " << ev.max(jac(PP).det())<<"\n";


    //! [Error and convergence rates]
    gsInfo<< "\nDoF_PDE = "<<std::scientific<<DoFPDE.transpose()<<"\n";
    gsInfo<< "L2_error_x = "<<std::scientific<<std::setprecision(3)<<l2err.transpose()<<"\n";
    //gsInfo<< "L2_error_y= "<<std::scientific<<std::setprecision(3)<<h1err.transpose()<<"\n";
    //! [Export visualization in ParaView]
    if (plot)
    {
        // constructing an IGA field (geometry + solution) for displacement
        gsField<> solutionField(assembler.patches(),solution);
        // constructing an IGA field (geometry + solution) for stresses
        gsField<> stressField(assembler.patches(),stresses,true);
        // analytical stresses
        gsField<> analyticalStressField(assembler.patches(),analyticalStresses,false);
        //... density function
        gsField<> densityField(assembler.patches(),f,false);

        // creating a container to plot all fields to one Paraview file
        std::map<std::string,const gsField<> *> fields;
        fields["Deformation"] = &solutionField;
        fields["Stress"] = &stressField;
        fields["StressAnalytical"] = &analyticalStressField;
        fields["DensityAnalytical"] = &densityField;
        gsWriteParaviewMultiPhysics(fields,"infinit_plate",10000,plot);
        gsInfo << "Open \"infinit_plate.pvd\" in Paraview for visualization. Stress wiggles on the left side are caused by "
                    "a singularity in the parametrization.\n";
        gsFileManager::open("infinit_plate.pvd");

        // gsInfo<<"Storing paraview...\n";
        // // Write the computed solution to paraview files
        // gsInfo<<"Making in Paraview...\n";
        // gsParaviewCollection collection("ParaviewOutput/solution", &ev);
        // collection.options().setSwitch("plotElements", true);
        // collection.options().setSwitch("base64", false);
        // collection.options().setInt("plotElements.resolution", 16);
        // collection.options().setInt("numPoints", 10000);
        // collection.newTimeStep(&geometry);
        // // collection.addField(istress,"numerical stress");
        // // collection.addField(jac(PP).det(), "Jacobian function");
        // // collection.addField(sigm_ex, "exact stress");
        // // collection.addField(ff_Ggeometry,"Density function");
        // collection.saveTimeStep();
        // collection.save();
        // //------------------------------------
        // gsInfo<<"Plotting in Paraview...\n";
        // // Run paraview
        // gsFileManager::open("ParaviewOutput/solution.pvd");
    }
    else
        gsInfo << "Done. No output created, re-run with --plot to get a ParaView "
                    "file containing the solution.\n";
    //! [Export visualization in ParaView]

    return EXIT_SUCCESS;
}// end main
