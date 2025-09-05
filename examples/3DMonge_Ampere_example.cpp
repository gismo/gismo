/** @file 3DMonge_Ampere_example.cpp

    @brief Tutorial on how to use expression assembler to solve a non-linear Monge-Ampere equation in 3D

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
    index_t numRefine   = 3;
    index_t numLRefine  = 1;
    index_t numReduce   = 0;
    index_t numElevate  = 0;
    index_t maxIter     = 30;
    double IntensityMAE = 9.;
    double quadValue    = 4.0;
    bool export_b64     = false;
    bool last           = false;
    index_t coef_V      = 1.;
    // Specify the file path
    // std::string fn("pde/example3D.xml");
    //std::string fn("volumes/GshapedVolume.xml");
    // Specify the file path
    //std::string fn("pde/quart_annulus.xml");
    //std::string fn("pde/infinit_plate.xml");
    std::string fn("pde/circle.xml");
    //std::string fn("surfaces/cylinder.xml"); 
    //std::string fn("domain2d/lake.xml");

    gsCmdLine cmd("Tutorial on solving a non-linear Monge-Ampere problem.");
    cmd.addInt("i", "iter", "Maximum number of iterations for the iterative Picard", maxIter);
    cmd.addInt( "e", "degreeElevation",
                "Number of degree Elevation steps to perform before solving (0: equalize degree in all directions)", numElevate );
    cmd.addInt( "r", "degreeRedution",
                "Number of degree Reduction steps to perform before solving (0: equalize degree in all directions)", numReduce );
    cmd.addInt( "u", "uniformRefine", "Number of Uniform h-refinement loops",  numRefine );
    cmd.addInt( "l", "numLRefine", "Number of local h-refinement loops",  numLRefine );
    cmd.addInt("c", "coef_V", "volume coef mult", coef_V);

    cmd.addString( "d", "file", "Input XML file data", fn );
    cmd.addInt("quRule",
                 "Quadrature rule [1:GaussLegendre,2:GaussLobatto,3:PatchRule]",
                 1);
    cmd.addSwitch("plot", "Create a ParaView visualization file with the solution", plot);
    cmd.addSwitch("last", "Solve solely for the last level of h-refinement",
                  last);
    cmd.addReal( "f", "IntensityMAE", "Intensity of density function",  IntensityMAE);
    cmd.addReal("q","quadValue", "Quadrature rule number of  points in each direction", quadValue);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    // Load the file
    gsFileData<> fd(fn);
    gsInfo << "Loaded file " << fd.lastPath() << "\n";
    // Create a gsMultipatch and add the loaded geometry
    // gsMultiPatch<> mpLeft; mpLeft.addPatch( gsNurbsCreator<>::BSplineCube(1,0,0,0) );
    // gsMultiPatch<> mpLeft; mpLeft.addPatch( gsNurbsCreator<>::NurbsSphere(1.,0.,0.,0.));
    // gsMultiPatch<> mpLeft = gsNurbsCreator<>::BSplineSquareGrid(1,1,1, 0.0, 0.0);
    // gsMultiPatch<> mpLeft = gsNurbsCreator<>::BSplineCubeGrid(1,1,1,1.,-0.5,-0.5,-0.5);
    // ...
    gsMultiPatch<> mpLeft;
    fd.getId(1,mpLeft);
    // Elevate and p-refine the basis to order p + numElevate
    // where p is the highest degree in the bases
    mpLeft.patch(0).coefs() *= coef_V;
    mpLeft.degreeElevate(numElevate);
    mpLeft.computeTopology();
    
    // Right-hand side function : Analytical density function rho_1
    // Load the file
    gsFunctionExpr<> f;
    fd.getId(2003, f);
    gsInfo<<"Density function "<< f << "\n";

    //! [Refinement]
    gsMultiBasis<double> dbasis(mpLeft, false);//true: poly-splines (not NURBS)

    //! [Problem setup]
    gsExprAssembler<> A(1,1);
    A.options().setReal("quA", quadValue);
    // A.options().setInt("quB", 1);
    A.options().setSwitch("SameElement",false);
    gsInfo<<"Active options:\n"<< A.options() <<"\n";

    typedef gsExprAssembler<>::geometryMap geometryMap;
    typedef gsExprAssembler<>::variable    variable;
    typedef gsExprAssembler<>::space       space;
    typedef gsExprAssembler<>::solution    solution;

    //::::::::::::::::::::      mesh adaptation solver         :::::::::::::::::::::::::
    while (dbasis.basis(0).numElements()<1e2)
    {
        dbasis.uniformRefine();
    }
    // h-refine each basis
    if (last)
    {
        for (int r =0; r < numRefine; ++r)
            dbasis.uniformRefine();
        numRefine = 0;
    }
    // ... 
    gsVector<>  h1err(numRefine+1), l2err(numRefine+1);
    gsVector<int>  DoFPDE(numRefine+1);
    for (int r=0; r<= numRefine; ++r)
    {
    dbasis.uniformRefine();
    // //... condition for the convergence 
    // while (dbasis.basis(0).numElements()<1e3)
    // {
    //     dbasis.uniformRefine();
    //     numRefine++;
    // }    
    // Elements used for numerical integration
    A.setIntegrationElements(dbasis);
    gsExprEvaluator<> ev(A);

    gsStopwatch timer;
    timer.restart();

    //... some infos on the computational domain
    auto corners         = dbasis.basis(0).support();
    gsInfo << "numElement :" << dbasis.basis(0).numElements() << " degree " << dbasis.degree() <<" dim " <<dbasis.dim()<<" Geodim " << mpLeft.geoDim() <<"\n";
    gsInfo << "corners : (";
    for (index_t i = 0; i < dbasis.dim()-1; i++)
    {
      gsInfo << corners.at(i)<< ",";
    }
    gsInfo << corners.at(dbasis.dim()-1)<<"),(";
    for (index_t i = 0; i < dbasis.dim()-1; i++)
    {
      gsInfo << corners.at(dbasis.dim()+i)<< ",";            
    }
    gsInfo <<  corners.at(2*dbasis.dim()-1)<< ")\n";

    /* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ###                                  Step 1-2 : Computes the density function
    ###                                     and the multipatch adaptive mapping
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
    gsAdaptiveMultiPatchBuilder MAE = gsAdaptiveMultiPatchBuilder(dbasis, mpLeft, maxIter, IntensityMAE, numReduce = numReduce);
    auto density   = MAE.buildAnalyticDensity(f);
    auto Psitp     = MAE.buildMultiPatch(density);// build the adaptive mapping
    // //------------------------------------
    geometryMap G  = A.getMap(mpLeft);
    geometryMap PP = A.getMap(Psitp);
    auto comp      = A.getCoeff(mpLeft, PP);
    gsMultiPatch<> PsiF  = MAE.buildCompMultiPatch(Psitp, quadValue); //composition of geometry maps
    // gsInfo << "Composition of geometry maps computed\n";
    geometryMap PPF = A.getMap(PsiF);
    // ...
    DoFPDE[r] = dbasis.basis(0).size();
    gsInfo << "DOF of the PDE space: "<< DoFPDE[r] <<"\n";
    timer.restart();
    l2err[r]  = EIGEN_PI*coef_V*coef_V + ev.integral(jac(PPF).det() );
    h1err[r]  = math::sqrt(ev.integralBdr((PPF-comp).sqNorm()));
    // ...
    std::cout << std::setprecision(15);
    std::cout << ", G :   volume Initial: "<< ev.integral( jac(G).det()  ) <<"\n";
    std::cout << ", G :   geometry volume: "<< EIGEN_PI*coef_V*coef_V + ev.integral( jac(G).det()  ) <<"\n";
    std::cout << ", Comp: geometry volume: "<< EIGEN_PI*coef_V*coef_V + ev.integral( (jac(PP).det() *jac(comp).det()) ) <<"\n";
    // ...
    std::cout << ", GPP : geometry volume: "<< l2err[r] <<"\n";
    std::cout << ", GPP : geometry boundary: "<<  h1err[r] <<"\n";
    //------------------------------------    
    //! [Export visualization in ParaView] 
    if (plot)
    {
        gsMultiPatch<> Psi;
        if ( dbasis.basis(0).IsRational){
        if (mpLeft.dim()== 3){
        for(size_t i =0; i<PsiF.nPatches(); ++i)
            Psi.addPatch(gsRationalTHBSpline<3>( dynamic_cast<const gsTensorNurbs<3>&>(PsiF.patch(i)) ));
            // Psi.addPatch(gsTHBSpline<3>( dynamic_cast<const gsTensorBSpline<3>&>(PsiF.patch(i)) ));
        }
        else{
        for(size_t i =0; i<PsiF.nPatches(); ++i)
            Psi.addPatch(gsRationalTHBSpline<2>( dynamic_cast<const gsTensorNurbs<2>&>(PsiF.patch(i)) ));
            //Psi.addPatch(gsTHBSpline<2>( dynamic_cast<const gsTensorBSpline<2>&>(PsiF.patch(i)) ));            
        }
        }
        else{
        if (mpLeft.dim()== 3){
        for(size_t i =0; i<PsiF.nPatches(); ++i)
            Psi.addPatch(gsRationalTHBSpline<3>( dynamic_cast<const gsTensorNurbs<3>&>(PsiF.patch(i)) ));
            // Psi.addPatch(gsTHBSpline<3>( dynamic_cast<const gsTensorBSpline<3>&>(PsiF.patch(i)) ));
        }
        else{
        for(size_t i =0; i<PsiF.nPatches(); ++i)
            Psi.addPatch(gsRationalTHBSpline<2>( dynamic_cast<const gsTensorNurbs<2>&>(PsiF.patch(i)) ));
            //Psi.addPatch(gsTHBSpline<2>( dynamic_cast<const gsTensorBSpline<2>&>(PsiF.patch(i)) ));            
        }
        }
        Psi.addAutoBoundaries();
        Psi.computeTopology();

        //Psi.uniformRefine();
        gsMultiBasis<> dbasis(Psi, true);//true: poly-splines (not NURBS)

        geometryMap PPF = A.getMap(Psi);
        auto ff_TG      = A.getCoeff(f, PPF);
        // --------------- adaptive refinement ---------------
        // Specify cell-marking strategy...
        MarkingStrategy adaptRefCrit = PUCA;
        //MarkingStrategy adaptRefCrit = GARU;
        //MarkingStrategy adaptRefCrit = errorFraction;
        real_t adaptRefParam = 0.7;
        // Elements used for numerical integration
        A.setIntegrationElements(dbasis);
        gsExprEvaluator<> ev(A);

        for (int r=0; r<=numLRefine; ++r)
        {
        //! [beginRefLoop]
            gsInfo << "====== Loop " << r << " of "
                    <<numLRefine<< " ====adapt Parameter ="<< adaptRefParam << " ======" << "\n";
            // --------------- error estimation/computation ---------------
            // Get the element-wise norms.
            ev.integralElWise( ( ff_TG ).sqNorm() );
            //ev.integralElWise( 1/jac(PPF).absDet() );

            const std::vector<real_t> eltErrs  = ev.elementwise();
            //! [errorComputation]

            //! [adaptRefinementPart]
            // Mark elements for refinement, based on the computed local errors and
            // the refinement-criterion and -parameter.
            std::vector<bool> elMarked( eltErrs.size() );
            gsMarkElementsForRef( eltErrs, adaptRefCrit, adaptRefParam, elMarked);
            gsInfo <<"Marked "<< std::count(elMarked.begin(), elMarked.end(), true) <<" elements.\n";
            // Refine the marked elements with a 1-ring of cells around marked elements
            gsRefineMarkedElements( dbasis, elMarked, 1);
            gsRefineMarkedElements( Psi, elMarked, 1);
            }

        //::::::::::::::::::::      end       :::::::::::::::::::::::::
        gsInfo<<"Plotting in Paraview...\n";
        gsParaviewCollection collection("ParaviewOutput/solution", &ev);
        collection.options().setSwitch("plotElements", true);
        collection.options().setSwitch("base64", export_b64);
        collection.options().setInt("plotElements.resolution", 16);
        collection.options().setInt("numPoints", 10000);
        collection.newTimeStep(&Psi);
        collection.addField(ff_TG, "density function");
        collection.addField(jac(PPF).det(), "Jacobian function");
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
    }
    // Assuming DoFPDE, l2err, and h1err are gsMatrix or similar types
    std::ofstream outFile("errorGeometry_analysis.txt", std::ios::app); // Open file in append mode
    if (outFile.is_open())
    {
        outFile << "#DoF_PDE:  \n"<< std::scientific << DoFPDE.transpose() << "\n";
        outFile << "#V_error: \n" << std::scientific << std::setprecision(3) << l2err.transpose() << "\n";
        outFile << "#B_error: \n" << std::scientific << std::setprecision(3) << h1err.transpose() << "\n";
        outFile << "#-------------------------------------------------------------------------------\n"; // Optional separator for readability
        outFile.close(); // Close the file after writing
    }
    else
    {
        gsInfo << "Error: Unable to open file for writing : error_analysis.txt.\n";
    }
    return EXIT_SUCCESS;


}// end main