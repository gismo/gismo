/** @file r_refinement_example.cpp

    @brief Tutorial on how to use compute the Monge-Ampere mapping for adaptive mesh in two and three dimensions.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris & M. BAHARI
*/

//! [Include namespace]
#include <gismo.h>
#include <fstream>  // For file operations
#include <gsAssembler/gsAdaptiveMultiPatchBuilder.h>  // Include the new class of r_refinement

using namespace gismo;
//! [Include namespace]

int main(int argc, char *argv[])
{
    //! [Parse command line]
    bool plot           = false;
    index_t numRefine   = 2;
    index_t numLRefine  = 1;
    index_t numRefineMAE= 3;
    index_t numReduce   = 0;
    index_t numElevate  = 0;
    index_t maxIter     = 50;
    index_t elevDegree  = 0; // degree elevation for the composition of geometry maps
    double IntensityMAE = 9.;
    double quadValue    = 2.0;
    bool bs_nrbs        = false;
    bool last           = true;
    bool colloc         = false;
    bool fit            = false;
    bool L2proj         = true;

    // gsStopwatch timer;
    // timer.restart();

    // Specify the file path
    // std::string fn("pde/example3D.xml");
    // std::string fn("volumes/GshapedVolume.xml");
    // Specify the file path
    std::string fn("pde/bsimple.xml");
    // std::string fn("pde/infinit_plate.xml");
    // std::string fn("pde/circle.xml");
    // std::string fn("surfaces/egg.xml"); 
    // std::string fn("domain2d/lake.xml");

    gsCmdLine cmd("Tutorial on solving a non-linear Monge-Ampere problem.");
    cmd.addInt("i", "iter", "Maximum number of iterations for the iterative Picard", 
                maxIter);
    cmd.addInt( "e", "degreeElevation", "Number of degree Elevation steps to perform before solving (0: equalize degree in all directions)", 
                numElevate );
    cmd.addInt( "v", "elevDegree", "Number of degree Elevation steps to perform fro the composition (0: equalize degree in all directions)", 
                elevDegree );
    cmd.addInt( "r", "degreeRedution", "Number of degree Reduction steps to perform before solving (0: equalize degree in all directions)", 
                numReduce );
    cmd.addInt( "u", "uniformRefine", "Number of Uniform h-refinement loops",  
                numRefine );
    cmd.addInt( "m", "uniformRefineMAE", "Number of Uniform h-refinement loops for MAE mapping",  
                numRefineMAE );
    cmd.addInt( "l", "numLRefine", "Number of local h-refinement loops",  
                numLRefine );
    cmd.addReal( "f", "IntensityMAE", "Intensity of density function",  
                IntensityMAE);
    cmd.addReal("q","quadValue", "Quadrature rule number of  points in each direction", 
                quadValue);
    cmd.addString( "d", "file", "Input XML file data", 
                fn );
    cmd.addInt("quRule", "Quadrature rule [1:GaussLegendre,2:GaussLobatto,3:PatchRule]",
                1);
    cmd.addSwitch("plot", "Create a ParaView visualization file with the solution", 
                plot);
    cmd.addSwitch("nrb", "basis use nrubs or bspline format",
                bs_nrbs);
    cmd.addSwitch("last", "Solve solely for the last level of h-refinement",
                last);
    cmd.addSwitch("colloc", "Compute the the compodition using collocation method",
                colloc);                
    cmd.addSwitch("fit", "Use fitting to compute the composition",
                fit);
    cmd.addSwitch("L2proj", "Use L2-projection to compute the composition",
                L2proj);
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    if(fit)
      bs_nrbs = true;// fitting doesn't support nurbs format
    // Load the file
    gsFileData<> fd(fn);
    // ...
    gsMultiPatch<> mpLeft, mpPsi;// Initial geometry and the resluted adaptive mapping
    fd.getId(1,mpLeft);
    gsInfo << "Loaded file " << fd.lastPath() << "\n";
    // Create a gsMultipatch and add the loaded geometry
    // gsMultiPatch<> mpLeft; mpLeft.addPatch( gsNurbsCreator<>::BSplineCube(1,0,0,0) );
    //gsMultiPatch<> mpLeft; mpLeft.addPatch( gsNurbsCreator<>::NurbsSphere(1.,0.,0.,0.));
    // mpLeft = gsNurbsCreator<>::BSplineSquareGrid(1,1,1, 0.0, 0.0);
    // mpLeft = gsNurbsCreator<>::BSplineCubeGrid(1,1,1,1.,-0.5,-0.5,-0.5);

    // Elevate and p-refine the basis to order p + numElevate
    // where p is the highest degree in the bases
    mpLeft.computeTopology();
    
    // Right-hand side function : Analytical density function rho_1
    // Load the file
    gsFunctionExpr<> f;
    fd.getId(2003, f);
    gsInfo<<"Density function "<< f << "\n";

    //! [Refinement]
    gsMultiBasis<double> dbasis(mpLeft, bs_nrbs);//true: poly-splines (not NURBS)
    // Elevate and p-refine the basis to order p + numElevate
    // where p is the highest degree in the bases
    dbasis.setDegree( dbasis.maxCwiseDegree() + numElevate);
    //dbasis.degreeIncrease(numElevate);
    
    //! [Problem setup]
    gsExprAssembler<> A(1,1);
    A.options().setReal("quA", quadValue);
    // A.options().setInt("quB", 1);
    A.options().setSwitch("SameElement",false);
    gsInfo<<"Active options:\n"<< A.options() <<"\n";

    // Elements used for numerical integration
    A.setIntegrationElements(dbasis);
    gsExprEvaluator<> ev(A);

    typedef gsExprAssembler<>::geometryMap geometryMap;

    /* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ###                                  Step r* : Computes the density function
    ###                                     and the multipatch adaptive mapping from a given mesh
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
    gsAdaptiveMultiPatchBuilder MAE = gsAdaptiveMultiPatchBuilder(mpLeft, numRefineMAE, maxIter, IntensityMAE, numReduce);
    auto density        = MAE.buildAnalyticDensity(f); // build the density function (we avoid composing rho o F o Psi here)
    MAE.buildMultiPatch(density, 1e-8);// build the adaptive mapping
    // //------------------------------------
    geometryMap G       = A.getMap(mpLeft);
    geometryMap PP      = A.getMap(MAE.MAmapping);
    auto Cmp            = A.getCoeff(mpLeft, PP);

    // h-refine each basis
    if (last)
    {
        for (int r =0; r < numRefine; ++r){
            dbasis.uniformRefine();
        }
        numRefine = 0;
    }
    // ... 
    gsVector<>  Bdrerr(numRefine+1), Volerr(numRefine+1)// L2 projection errors
                , CVolerr(numRefine+1)
                , IBdrerr(numRefine+1), IVolerr(numRefine+1)// Interpolation errors
                , FBdrerr(numRefine+1), FVolerr(numRefine+1);// Fitting errors
    gsVector<int>  DoFPDE(numRefine+1);
    for (int r=0; r<= numRefine; ++r)
    {
    dbasis.uniformRefine();

    //... some infos on the computational domain
    gsInfo << r <<"th iter:{ numElement " << dbasis.basis(0).numElements() << " degree " << dbasis.degree() 
            <<" dim " <<dbasis.dim()<<" Geodim " << mpLeft.geoDim() 
            <<"}------------------------------------------------------\n";

    DoFPDE[r]               = dbasis.basis(0).size();
    CVolerr[r]              = std::abs(abs(ev.integral( meas(G)  )) - abs( ev.integral( jac(Cmp).det()*jac(PP).det() ) ));

    //----------------------------------------------------------------------
    // ... computes the composition of geometry maps L^2-projection method !
    //----------------------------------------------------------------------
    if(L2proj)
    {
    mpPsi           = MAE.buildCompMultiPatch(dbasis, quadValue);
    geometryMap GPi = A.getMap(mpPsi);
    // ... Error analysis
    Volerr[r]               = std::abs(abs(ev.integral( meas(G)  )) - abs( ev.integral(meas(GPi)) ));
    Bdrerr[r]               = math::sqrt(ev.integral((GPi-Cmp).sqNorm()));
    }

    if (colloc){
    //---------------------------------------------------------- 
    //...Interpolation of the mapping by collocation method !
    //----------------------------------------------------------
    gsMultiPatch<> mpPsi    = MAE.buildColCompMultiPatch(dbasis);
    geometryMap PGI         = A.getMap(mpPsi);
    // ... Error using Interpolation method
    IVolerr[r]              = std::abs(abs(ev.integral( meas(G)  )) - abs( ev.integral(meas(PGI)) ));
    IBdrerr[r]              = math::sqrt(ev.integral((PGI-Cmp).sqNorm())); 
    }

    if(fit){
    //---------------------------------------------------------- 
    //...Interpolation of the mapping by fit method !
    //----------------------------------------------------------
    gsMultiPatch<> mpPsi    = MAE.buildFitCompMultiPatch(dbasis,50);
    geometryMap PGF         = A.getMap(mpPsi);
    // ... Error analysis
    FVolerr[r]              = std::abs(abs(ev.integral( meas(G)  )) - abs( ev.integral(meas(PGF)) ));
    FBdrerr[r]              = math::sqrt(ev.integral((PGF-Cmp).sqNorm())); 
    }
    }

    // Assuming DoFPDE, l2err, and h1err are gsMatrix or similar types
    std::ofstream outFile("errorGeometry_analysis.txt", std::ios::app); // Open file in append mode
    if (outFile.is_open())
    {
        outFile << "#DoF_PDE: q"<< quadValue<<"pPr"<< dbasis.basis(0).maxDegree()<<"pPsi"<< mpLeft.basis(0).maxDegree()-numReduce <<"\n"
                << std::scientific << DoFPDE.transpose() << "\n";
        if (L2proj){
        outFile << "#V_error: \n" << std::scientific << std::setprecision(3) << Volerr.transpose() << "\n";
        outFile << "#L_error: \n" << std::scientific << std::setprecision(3) << Bdrerr.transpose() << "\n";
        }
        if (colloc){
        outFile << "#INTERPOL V_error: \n" << std::scientific << std::setprecision(3) << IVolerr.transpose() << "\n";
        outFile << "#INTERPOL L_error: \n" << std::scientific << std::setprecision(3) << IBdrerr.transpose() << "\n";
        }
        if (fit){
        outFile << "#FITTING V_error: \n" << std::scientific << std::setprecision(3) << FVolerr.transpose() << "\n";
        outFile << "#FITTING L_error: \n" << std::scientific << std::setprecision(3) << FBdrerr.transpose() << "\n";
        }
        outFile << "#C_error:  "<< quadValue << ": "<< std::scientific << std::setprecision(3) << CVolerr.transpose() << "\n";
        outFile << "#-------------------------------------------------------------------------------\n"; // Optional separator for readability
        outFile.close(); // Close the file after writing
        gsInfo << "Error analysis results appended to errorGeometry_analysis.txt.\n";
    }
    else
    {
        gsInfo << "Error: Unable to open file for writing : error_analysis.txt.\n";
    }
    //------------------------------------    
    //! [Error and convergence rates]
    // --- Print header ---
    gsInfo << std::setw(12) << "DoFs" << " & "
         << std::setw(13) << "V(F∘Ψ)"     << " & "
         << std::setw(13) << "V(Πp(F∘Ψ))"  << " & "
         << std::setw(13) << "L2(Πp(F∘Ψ))" << " & "
         << std::setw(6)  << "EOcL2"      << "\n";
    gsInfo << std::string(50, '-') << "\n";
    // --- Print table row by row ---
    auto orderofConv = ( Bdrerr.head(numRefine).array() /
                  Bdrerr.tail(numRefine).array() ).log().transpose() / std::log(2.0);
    gsInfo << std::setw(12) << DoFPDE[0] << " & "
            << std::setw(12) <<std::setprecision(3)<<std::scientific<< CVolerr[0] << " & "
            << std::setw(12) <<std::setprecision(3)<<std::scientific<< Volerr[0] << " & "
            << std::setw(12) <<std::setprecision(3)<<std::scientific<< Bdrerr[0] << "&"
            << std::setw(12) << "--" << "\n";
    for (int i = 1; i <= numRefine; i++) {
        gsInfo << std::setw(12) << DoFPDE[i] << " & "
             << std::setw(12) <<std::setprecision(3)<<std::scientific<< CVolerr[i] << " & "
             << std::setw(12) <<std::setprecision(3)<<std::scientific<< Volerr[i] << " & "
             << std::setw(12) <<std::setprecision(3)<<std::scientific<< Bdrerr[i] << "&"
             << std::setw(12) <<std::fixed<<std::setprecision(2)<< orderofConv[i-1] << "\n";
            }

    //! [Export visualization in ParaView] 
    if (plot)
    {
        gsMultiPatch<> Psi;
        if ( mpPsi.basis(0).weights().any()){
            gsInfo<<"Rational mapping \n";
        if (mpLeft.dim()== 3){
        for(size_t i =0; i<mpPsi.nPatches(); ++i)
            Psi.addPatch(gsRationalTHBSpline<3>( dynamic_cast<const gsTensorNurbs<3>&>(mpPsi.patch(i)) ));
        }
        else{
        for(size_t i =0; i<mpPsi.nPatches(); ++i)
            Psi.addPatch(gsRationalTHBSpline<2>( dynamic_cast<const gsTensorNurbs<2>&>(mpPsi.patch(i)) ));
        }
        }
        else{
            gsInfo<<"nonRational mapping \n";
        if (mpLeft.dim()== 3){
        for(size_t i =0; i<mpPsi.nPatches(); ++i)
            Psi.addPatch(gsTHBSpline<3>( dynamic_cast<const gsTensorBSpline<3>&>(mpPsi.patch(i)) ));
        }
        else{
        for(size_t i =0; i<mpPsi.nPatches(); ++i)
            Psi.addPatch(gsTHBSpline<2>( dynamic_cast<const gsTensorBSpline<2>&>(mpPsi.patch(i)) ));            
        }
        }
        Psi.addAutoBoundaries();
        Psi.computeTopology();

        //Psi.uniformRefine();
        gsMultiBasis<> basis(Psi, true);//true: poly-splines (not NURBS)
        // Elements used for numerical integration
        A.setIntegrationElements(basis);
        gsExprEvaluator<> ev(A);
        // ...
        geometryMap PPsi = A.getMap(Psi);
        auto ff_TG       = A.getCoeff(f, PPsi);
        // --------------- adaptive refinement ---------------
        // Specify cell-marking strategy...
        MarkingStrategy adaptRefCrit = PUCA;
        //MarkingStrategy adaptRefCrit = GARU;
        //MarkingStrategy adaptRefCrit = errorFraction;
        real_t adaptRefParam = 0.7;

        for (int r=0; r<=numLRefine; ++r)
        {
        //! [beginRefLoop]
            gsInfo << "====== Loop " << r << " of "
                    <<numLRefine<< " ====adapt Parameter ="<< adaptRefParam << " ======" << "\n";
            // --------------- error estimation/computation ---------------
            // Get the element-wise norms.
            ev.integralElWise( ff_TG.val() );
            //ev.integralElWise( 1/jac(GPi).absDet() );

            const std::vector<real_t> eltErrs  = ev.elementwise();
            //! [errorComputation]

            //! [adaptRefinementPart]
            // Mark elements for refinement, based on the computed local errors and
            // the refinement-criterion and -parameter.
            std::vector<bool> elMarked( eltErrs.size() );
            gsMarkElementsForRef( eltErrs, adaptRefCrit, adaptRefParam, elMarked);
            gsInfo <<"Marked "<< std::count(elMarked.begin(), elMarked.end(), true) <<" elements.\n";
            // Refine the marked elements with a 1-ring of cells around marked elements
            gsRefineMarkedElements( basis, elMarked, 1);
            gsRefineMarkedElements( Psi, elMarked, 1);
        }

        //::::::::::::::::::::      end       :::::::::::::::::::::::::
        gsInfo<<"Plotting in Paraview...\n";
        gsParaviewCollection collection("ParaviewOutput/solution", &ev);
        collection.options().setSwitch("plotElements", true);
        collection.options().setSwitch("base64", false);
        collection.options().setInt("plotElements.resolution", 16);
        collection.options().setInt("numPoints", 10000);
        collection.newTimeStep(&Psi);
        collection.addField(ff_TG, "density function");
        collection.addField(jac(PPsi).det(), "Jacobian function");
        collection.saveTimeStep();
        collection.save();

        gsFileManager::open("ParaviewOutput/solution.pvd");
        // gsWrite(Psi, "Psi_mapping");
        // gsInfo << "Result written in Psi_mapping.xml \n";
    }
    else
        gsInfo << "Done. No output created, re-run with --plot to get a ParaView "
                  "file containing the solution.\n";

    return EXIT_SUCCESS;


}// end main