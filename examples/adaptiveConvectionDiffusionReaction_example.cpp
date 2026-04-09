/** @file adaptiveConvectionDiffusionReaction_example.cpp

    @brief Tutorial on how to use G+Smo to solve a convection-diffusion-reaction problem.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): S. Kleiss
*/

//! [Include namespace]
# include <gismo.h>
# include <gsAssembler/gsAdaptiveRefUtils.h>

using namespace std;
using namespace gismo;
//! [Include namespace]

int main(int argc, char *argv[])
{
   //! [Parse command line]
   bool plot = false;

   gsCmdLine cmd("Example for solving a convection-diffusion problem.");
   cmd.addSwitch("plot", "Create a ParaView visualization file with the solution", plot);
   try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }
   //! [Parse command line]

   //! [Function data] TEST   2------------------------------------------------------------------ SOLOVEV
//    gsFunctionExpr<> g("sqrt( 4.*y**2 + ( x**2-1.)**2 )", 2);
//    // Define source function
//    gsFunctionExpr<> rhs("-1e0*(( ( (6.0*x**2 + 2.0) / sqrt(4.0*y**2 + (x**2 - 1.0)**2)) - ((4.0*x**2*(x**2 - 1.0)**2 + 16.0*y**2) / (sqrt(4.0*y**2 + (x**2 - 1.0)**2)**3)) ))",2);

//    // diffusion coefficient:
//    gsFunctionExpr<> coeff_diff("1e0+(1e+4-1e0)*(y**2)/( y**2 + 0.25 * x**2 * (x**2 - 1.0)**2)",
//                                "    (1e+4-1e0)*(-0.5 * x * (x**2 - 1.0) * y)/( y**2 + 0.25 * x**2 * (x**2 - 1.0)**2)",
//                                "    (1e+4-1e0)*(-0.5 * x * (x**2 - 1.0) * y)/( y**2 + 0.25 * x**2 * (x**2 - 1.0)**2)",
//                                "1e0+(1e+4-1e0)*(0.25 * x**2 * (x**2 - 1.0)**2)/( y**2 + 0.25 * x**2 * (x**2 - 1.0)**2)",2);
   //! [Function data] TEST   2------------------------------------------------------------------SOVINEC
   gsFunctionExpr<> g("cos(pi*x)*cos(pi*y)", 2);
   // Define source function
   gsFunctionExpr<> rhs("2*pi*pi*cos(pi*x)*cos(pi*y)",2);

   // diffusion coefficient:
gsFunctionExpr<> coeff_diff(
    "if(x**2+y**2 <=1e-16, 1e0, 1e0 + (1e8-1e0)*(cos(pi*x)^2*sin(pi*y)^2)/(sin(pi*x)^2*cos(pi*y)^2 + cos(pi*x)^2*sin(pi*y)^2))",
    "if(x**2+y**2 <=1e-16, 0,         (1e8-1e0)*(-sin(pi*x)*cos(pi*y)*cos(pi*x)*sin(pi*y))/(sin(pi*x)^2*cos(pi*y)^2 + cos(pi*x)^2*sin(pi*y)^2))",
    "if(x**2+y**2 <=1e-16, 0,         (1e8-1e0)*(-sin(pi*x)*cos(pi*y)*cos(pi*x)*sin(pi*y))/(sin(pi*x)^2*cos(pi*y)^2 + cos(pi*x)^2*sin(pi*y)^2))",
    "if(x**2+y**2 <=1e-16, 1e0, 1e0 + (1e8-1e0)*(sin(pi*x)^2*cos(pi*y)^2)/(sin(pi*x)^2*cos(pi*y)^2 + cos(pi*x)^2*sin(pi*y)^2+1e-16))", 2);
    // convection coefficient:
   gsFunctionExpr<> coeff_conv("0.","0.",2);
   // reaction coefficient:
   gsFunctionExpr<> coeff_reac("0",2);
   //! [Function data]
   // Print out source function and solution
   gsInfo<<"Source function " << rhs << "\n";
   gsInfo<<"Dirichlet boundary conditions "  << g << "\n\n";


   // --------------- read geometry from file ---------------

   // Read geometry from file
   //! [GetGeometryData]
   // Read xml and create gsMultiPatch
//    string fileSrc( "planar/lshape2d_3patches_thb.xml" );
//    gsMultiPatch<real_t> patches;
//    gsReadFile<real_t>( fileSrc, patches);
// Test 2 
   string fileSrc( "pde/unitsquare.xml" );
//    string fileSrc( "pde/solovev_relaxed.xml" );
//    string fileSrc( "pde/solovev_relaxed.xml" );
//    string fileSrc( "pde/annulus2d_bvp.xml" );
    gsFileData<> fd(fileSrc);
    gsMultiPatch<> Psi, patches;
    fd.getId(1,Psi);
    if ( Psi.basis(0).weights().any()){
        gsInfo<<Psi.basis(0).weights().any()<<" Rational mapping \n";
    for(size_t i =0; i<Psi.nPatches(); ++i)
        patches.addPatch(gsRationalTHBSpline<2>( dynamic_cast<const gsTensorNurbs<2>&>(Psi.patch(i)) ));
    }
    else{
        gsInfo<<"nonRational mapping \n";
    for(size_t i =0; i<Psi.nPatches(); ++i)
        patches.addPatch(gsTHBSpline<2>( dynamic_cast<const gsTensorBSpline<2>&>(Psi.patch(i)) ));            
    }
   //! [computeTopology]
   // Get all interfaces and boundaries:
   patches.computeTopology();
   //! [computeTopology]

   //! [GetGeometryData]
   gsInfo << "The domain is a "<< patches <<"\n";

   // --------------- add bonudary conditions ---------------
   //! [Boundary conditions]
   gsBoundaryConditions<> bcInfo;

    // fd.getId(2, bcInfo); // id=2: boundary conditions
    // bcInfo.setGeoMap(patches);
    // gsInfo<<"Boundary conditions:\n"<< bcInfo <<"\n";
   // For simplicity, set Dirichlet boundary conditions
   // given by exact solution g on all boundaries:
   for ( gsMultiPatch<>::const_biterator
            bit = patches.bBegin(); bit != patches.bEnd(); ++bit)
   {
       bcInfo.addCondition( *bit, condition_type::dirichlet, &g );
   }
   //! [Boundary conditions]

   // --------------- define Pde ---------------
   //! [definePde]
   gsConvDiffRePde<real_t> cdrPde(patches, bcInfo, & coeff_diff,& coeff_conv, & coeff_reac, & rhs);
   //! [definePde]


   // --------------- set up basis ---------------
   // Specify cell-marking strategy...
   MarkingStrategy adaptRefCrit = PUCA;
   //MarkingStrategy adaptRefCrit = GARU;
   //MarkingStrategy adaptRefCrit = errorFraction;

   //! [GetBasisFromTHB]
   // Copy basis from the geometry
   gsMultiBasis<> bases( patches , true);
   bases.uniformRefine(4);
//    bases.degreeElevate(1);
   //! [GetBasisFromTHB]


   //! [initialRefinements]
   // Number of initial uniform refinement steps:
//    int numInitUniformRefine  = 1;

//    for (int i = 0; i < numInitUniformRefine; ++i){
//      bases.uniformRefine();
//     }
   //! [initialRefinements]


   // --------------- set up adaptive refinement loop ---------------
   // ... and parameter.
   const real_t adaptRefParam = 0.7;

   //! [adaptRefSettings]
   // Number of refinement loops to be done
   int numRefinementLoops = 7;


   //! [constructAssembler]
   // Construct assembler
   gsCDRAssembler<real_t> cdrAss( cdrPde, bases);
   // Set stabilization flag to 1 = SUPG
   cdrAss.options().setInt("Stabilization", stabilizerCDR::none);
   // Compute Dirichlet values by L2-projection
   // Caution: Interpolation does not work for locally refined (T)HB-splines!
   cdrAss.options().setInt("DirichletValues",dirichlet::l2Projection);
   //! [constructAssembler]

   // --------------- adaptive refinement loop ---------------
    gsVector<>     l2err(numRefinementLoops+1), h1err(numRefinementLoops+1),l1err(numRefinementLoops+1);
    gsVector<int>  DoFPDE(numRefinementLoops+1);
   //! [beginRefLoop]
   for( int refLoop = 0; refLoop <= numRefinementLoops; refLoop++)
   {
   //! [beginRefLoop]
       gsInfo << "====== Loop " << refLoop << " of "
              <<numRefinementLoops<< " ======" << "\n";

       // --------------- solving ---------------

       //! [solverPart]
       // Generate system matrix and load vector
       cdrAss.assemble();
        gsInfo << "System matrix has " << cdrAss.matrix().nonZeros() << " non-zero entries.\n";
       // Solve the system
       gsMatrix<real_t> solVector =
           gsSparseSolver<>::CGDiagonal( cdrAss.matrix() ).solve( cdrAss.rhs());

       // Construct the solution as a scalar field
       gsField<> solField;
       solField = cdrAss.constructSolution(solVector);
       //! [solverPart]

       // --------------- error estimation/computation ---------------

       //! [errorComputation]
       // Compute the H1-seminorm of the computed solution
       // ( which is, at least, equivalent to the energy norm in this example )
       // using the known exact solution.
       gsExprEvaluator<> ev;
       ev.setIntegrationDomain(cdrAss.multiBasis().domain());
       gsExprEvaluator<>::geometryMap Gm = ev.getMap(patches);
       gsExprEvaluator<>::variable is = ev.getVariable(solField.fields());
       auto ms = ev.getVariable(g, Gm);

       // Get the element-wise norms.
       gsInfo << " error L2 norm: "<< sqrt(ev.integral( (is-ms).sqNorm() )) << "\n";
       gsInfo << " error H1 seminorm: "<< sqrt(ev.integral( ( igrad(is,Gm) - igrad(ms,Gm)).sqNorm() )) << "\n";
        //...
        gsMatrix<> A(2,1);
        A << 0.5,0.5; // parametric coordinates for the isogeometric solution
        gsInfo<< ev.eval(ms, A) <<"\n";
        gsInfo<< ev.eval(is, A, 0) <<"\n";
        DoFPDE[refLoop]         = cdrAss.multiBasis().totalSize();
        l1err[refLoop]          = abs(1.-1./ev.eval(is, A).value());
        l2err[refLoop]          = sqrt(ev.integral( (is-ms).sqNorm() ));
        h1err[refLoop]          = sqrt(ev.integral( ( igrad(is,Gm) - grad(ms)).sqNorm() ));

       // --------------- adaptive refinement ---------------
       ev.integralElWise( ( is- ms).sqNorm()*meas(Gm) );
        //------------------
       const std::vector<real_t> & eltErrs  = ev.elementwise();
       //! [errorComputation]
       //! [adaptRefinementPart]
       // Mark elements for refinement, based on the computed local errors and
       // the refinement-criterion and -parameter.
       std::vector<bool> elMarked( eltErrs.size() );
       gsMarkElementsForRef( eltErrs, adaptRefCrit, adaptRefParam, elMarked);
       gsInfo <<"Marked "<< std::count(elMarked.begin(), elMarked.end(), true) <<" elements.\n";

       // Refine the marked elements with a 1-ring of cells around marked elements
        gsRefineMarkedElements( cdrAss.multiBasis(), elMarked, 1 );
        gsInfo <<"Refined mesh has degrees of freedom.\n";
        //! [adaptRefinementPart]


        //! [repairInterfaces]
        // Call repair interfaces to make sure that the new meshes
        // match along patch interfaces.
        // cdrAss.multiBasis().repairInterfaces( patches.interfaces() );
        //! [repairInterfaces]
        // gsInfo <<"Refined mesh has degrees of freedom.\n";

        //! [refreshAssembler]
        cdrAss.refresh();
        gsInfo <<"Refined mesh has degrees of freedom.\n";

       //! [refreshAssembler]

       //! [Export to Paraview]
       // Export the final solution
       if( plot && refLoop == numRefinementLoops )
       {
           // Write the computed solution to paraview files
           gsWriteParaview<>( solField, "adaptRef", 1000, true);
       }
       //! [Export to Paraview]

   }
    //! [Error and convergence rates]
    gsInfo<< "\nDoF_PDE = "<<std::scientific<<DoFPDE.transpose()<<"\n";
    gsInfo<< "L1_error = "<<std::scientific<<std::setprecision(3)<<l1err.transpose()<<"\n";
    gsInfo<< "L2_error = "<<std::scientific<<std::setprecision(3)<<l2err.transpose()<<"\n";
    gsInfo<< "H1_error= "<<std::scientific<<std::setprecision(3)<<h1err.transpose()<<"\n";
    //! [Error and convergence rates]
    if (numRefinementLoops>0)
    {
        gsInfo<< "\nEoC (L1): " << std::fixed<<std::setprecision(2)
              <<  ( l1err.head(numRefinementLoops).array()  /
                   l1err.tail(numRefinementLoops).array() ).log().transpose() / std::log(2.0)
                   <<"\n";

        gsInfo<< "\nEoC (L2): " << std::fixed<<std::setprecision(2)
              <<  ( l2err.head(numRefinementLoops).array()  /
                   l2err.tail(numRefinementLoops).array() ).log().transpose() / std::log(2.0)
                   <<"\n";

        gsInfo<<   "EoC (H1): "<< std::fixed<<std::setprecision(2)
             <<( h1err.head(numRefinementLoops).array() /
                  h1err.tail(numRefinementLoops).array() ).log().transpose() / std::log(2.0) <<"\n";
    }
    // Assuming DoFPDE, l2err, and h1err are gsMatrix or similar types
    std::ofstream outFile("error_analysis.txt", std::ios::app); // Open file in append mode
    if (outFile.is_open())
    {
        outFile << "#DoF_PDE: " << bases.degree(0) << " \n"<< std::scientific << DoFPDE.transpose() << "\n";
        outFile << "#L1_error: \n" << std::scientific << std::setprecision(3) << l1err.transpose() << "\n";
        outFile << "#L2_error: \n" << std::scientific << std::setprecision(3) << l2err.transpose() << "\n";
        outFile << "#H1_error: \n" << std::scientific << std::setprecision(3) << h1err.transpose() << "\n";
        outFile << "#-------------------------------------------------------------------------------\n"; // Optional separator for readability
        outFile.close(); // Close the file after writing
    }
    else
    {
        gsInfo << "Error: Unable to open file for writing : error_analysis.txt.\n";
    }

   //! [Plot in Paraview]
   if( plot )
   {
       // Run paraview
       gsFileManager::open("adaptRef.pvd");
   }
   //! [Plot in Paraview]
   else
   {
       gsInfo<<"Done. No output created, re-run with --plot to get a ParaView "
               "file containing Plotting image data.\n";
   }
   return EXIT_SUCCESS;

}// end main