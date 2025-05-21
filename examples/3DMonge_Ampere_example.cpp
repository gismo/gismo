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
    index_t numElevate  = 0;
    index_t maxIter     = 30;
    double IntensityMAE = 12.;
    bool export_b64     = false;

    // Specify the file path
    //std::string fn("pde/example3D.xml");
    //std::string fn("volumes/GshapedVolume.xml");
    // Specify the file path
    //std::string fn("pde/quart_annulus.xml");
    //std::string fn("pde/mhd.xml");
    //std::string fn("pde/infinit_plate.xml");
    //std::string fn("pde/circle.xml");
    std::string fn("surfaces/simple.xml"); 
    //std::string fn("surfaces/egg.xml");

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
    cmd.addSwitch("plot", "Create a ParaView visualization file with the solution", plot);
    cmd.addReal( "f", "IntensityMAE", "Intensity of density function",  IntensityMAE);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    // Load the file
    gsFileData<> fd(fn);
    gsInfo << "Loaded file " << fd.lastPath() << "\n";
    // Create a gsMultipatch and add the loaded geometry
    // gsMultiPatch<> mpLeft; mpLeft.addPatch( gsNurbsCreator<>::BSplineCube(1,0,0,0) );
    //gsMultiPatch<> mpLeft; mpLeft.addPatch( gsNurbsCreator<>::NurbsSphere(1.,0.,0.,0.));
    // gsMultiPatch<> mpLeft = gsNurbsCreator<>::BSplineSquareGrid(1,1,1, 0.0, 0.0);
    // gsMultiPatch<> mpLeft = gsNurbsCreator<>::BSplineCubeGrid(1,1,1,1.,-0.5,-0.5,-0.5);
    // ...
    gsMultiPatch<> mpLeft;
    fd.getId(1,mpLeft);
    // Elevate and p-refine the basis to order p + numElevate
    // where p is the highest degree in the bases
    // auto kv1 =  static_cast<gsTensorNurbs<2> &>( mpLeft.patch(0)).knots(0);
    // auto kv2 =  static_cast<gsTensorNurbs<2> &>( mpLeft.patch(0)).knots(1);
    mpLeft.degreeElevate(numElevate);
    mpLeft.computeTopology();

    // Right-hand side function : Analytical density function rho_1
    // Load the file
    gsFunctionExpr<> f;
    fd.getId(2003, f);
    gsInfo<<"Density function "<< f << "\n";

    //! [Refinement]
    gsMultiBasis<double> dbasis(mpLeft, true);//true: poly-splines (not NURBS)
    // TODO : build Hdiv solver
    //gsMultiBasis<double> Hdivbasis(mpLeft, true);//true: poly-splines (not NURBS)
    //Hdivbasis.degreeElevate(1);
    //dbasis.setDegree(2);


    //! [Problem setup]
    gsExprAssembler<> A(1,1);
    // A.options().setReal("quA", 2.0);
    // A.options().setInt("quB", 2);
    A.options().setSwitch("SameElement",false);
    gsInfo<<"Active options:\n"<< A.options() <<"\n";

    typedef gsExprAssembler<>::geometryMap geometryMap;
    typedef gsExprAssembler<>::variable    variable;
    typedef gsExprAssembler<>::space       space;
    typedef gsExprAssembler<>::solution    solution;

    // Elements used for numerical integration
    A.setIntegrationElements(dbasis);
    gsExprEvaluator<> ev(A);

    gsStopwatch timer;
    timer.restart();
    //::::::::::::::::::::      mesh adaptation solver         :::::::::::::::::::::::::
    for (int r=0; r<=numRefine; ++r)
    {
        dbasis.uniformRefine();
        // mpLeft.uniformRefine();
    }
    //... condition for the convergence 
    while (dbasis.basis(0).numElements()<1e3)
    {
        dbasis.uniformRefine();
        numRefine++;
    }    
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
    gsAdaptiveMultiPatchBuilder MAE = gsAdaptiveMultiPatchBuilder(dbasis, mpLeft, numElevate, maxIter, IntensityMAE);
    auto density = MAE.buildAnalyticDensity(f);
    // //------------------------------------
    auto Psitp   = MAE.buildMultiPatch(density, true);//if true : composition of geometry maps
    gsInfo << "MultiPatch geometry computed\n";

    //! [Export visualization in ParaView]
    if (plot)
    {
        gsMultiPatch<> Psi;
        if (mpLeft.dim()== 3){
        for(size_t i =0; i<Psitp.nPatches(); ++i)
            Psi.addPatch(gsTHBSpline<3>( dynamic_cast<const gsTensorBSpline<3>&>(Psitp.patch(i)) ));
        }
        else{
        for(size_t i =0; i<Psitp.nPatches(); ++i)
            Psi.addPatch(gsTHBSpline<2>( dynamic_cast<const gsTensorBSpline<2>&>(Psitp.patch(i)) ));            
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

    return EXIT_SUCCESS;


}// end main