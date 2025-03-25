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
    index_t numRefine   = 4;
    index_t numLRefine  = 3;
    index_t numElevate  = 0;
    index_t numrRefine  = -1; // number of composition bewteen adaptive mappings ()
    index_t maxIter     = 30;
    double IntensityMAE = 9.;
    bool export_b64     = false;
    bool errorsave      = false;

    // Specify the file path
    std::string fn("pde/example3D.xml");

    gsCmdLine cmd("Tutorial on solving a non-linear Monge-Ampere problem.");
    cmd.addInt("i", "iter", "Maximum number of iterations for the iterative Picard", maxIter);
    cmd.addInt( "e", "degreeElevation",
                "Number of degree elevation steps to perform before solving (0: equalize degree in all directions)", numElevate );
    cmd.addInt( "u", "uniformRefine", "Number of Uniform h-refinement loops",  numRefine );
    cmd.addInt( "l", "numLRefine", "Number of local h-refinement loops",  numLRefine );
    cmd.addInt( "r", "numrRefine", "Number of local r-refinement compostion loops",  numrRefine);

    cmd.addString( "d", "file", "Input XML file data", fn );
    cmd.addInt("quRule",
                 "Quadrature rule [1:GaussLegendre,2:GaussLobatto,3:PatchRule]",
                 1);
    cmd.addSwitch("plot", "Create a ParaView visualization file with the solution", plot);
    cmd.addSwitch("errorsave", "Create a file in ... and save errors", errorsave);
    cmd.addReal( "f", "IntensityMAE", "Intensity of density function",  IntensityMAE);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    // Load the file
    gsFileData<> fd(fn);
    gsInfo << "Loaded file " << fd.lastPath() << "\n";
    // Create a gsMultipatch and add the loaded geometry
    gsMultiPatch<> mpLeft; mpLeft.addPatch( gsNurbsCreator<>::BSplineCube(1,0,0,0) );
    //gsMultiPatch<> mpLeft = gsNurbsCreator<>::BSplineCubeGrid(1,1,1,1.,0.,0.,0.);
    // fd.getId(1,mpLeft);
    // Elevate and p-refine the basis to order p + numElevate
    // where p is the highest degree in the bases
    mpLeft.degreeElevate(numElevate);
    mpLeft.computeTopology();
    gsInfo << "INFO IN PARAMETRIC DOMAIN "<< mpLeft.dim() << mpLeft.parDim() <<"\n";

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

    // Set pow for BFO method dim in parameteric domain
    auto IGdim     = dbasis.dim();

    gsStopwatch timer;
    timer.restart();
    //::::::::::::::::::::      mesh adaptation solver         :::::::::::::::::::::::::
    for (int r=0; r<=numRefine; ++r)
    {
        dbasis.uniformRefine();
        // mp.uniformRefine();
        // mpLeft.uniformRefine();
    }
    timer.restart();
    /* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ###   Step 1-2 : Computes the density function
    ###         and the multipatch adaptive mapping
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
    gsAdaptiveMultiPatchBuilder MAE = gsAdaptiveMultiPatchBuilder(dbasis, mpLeft, numElevate, maxIter, IntensityMAE);
    auto density = MAE.buildAnalyticDensity(f);
    auto Psitp   = MAE.buildMultiPatch(density, false);

    //! [Export visualization in ParaView]
    if (plot)
    {
        space v        = A.getSpace(dbasis);
        gsMatrix<> vsolVector;
        solution v_sol = A.getSolution(v, vsolVector);

        for (index_t i = 0; i < numrRefine; i++){
            /// compose adaptive mappings : not working after 2nd composition.(can be used for the blew-ip prblem)
            geometryMap PP = A.getMap(Psitp);
            auto  comp = PP(Psitp);
            A.initSystem(IGdim);
            //Obtain control points for the gradient of mpLeft.comp(Psi)
            A.assemble( v * v.tr() , v * comp.tr() );// blocked by this one
            vsolVector = MAE.Poisson.L2ProjectVec(A.rhs());
            v_sol.extract(Psitp);
            Psitp.addAutoBoundaries();
            Psitp.computeTopology();
        }
        //::::::::::::::::::::    Compute the composition of geometry maps      :::::::::::::::::::::::::
        geometryMap PP = A.getMap(Psitp);
        // auto comp  = A.getCoeff(mpLeft, PP);
        // A.initSystem(IGdim);
        // //Obtain control points for the gradient of mpLeft.comp(Psi)
        // A.assemble( v * v.tr() , v * comp.tr() );// blocked by this one
        // vsolVector = MAE.Poisson.L2ProjectVec(A.rhs());
        // v_sol.extract(Psitp);
        // Psitp.addAutoBoundaries();
        Psitp.computeTopology();
        gsInfo << "end of adaptive mapping computation\n" << Psitp<< "\n";
        gsMultiPatch<> Psi;
        for(size_t i =0; i<Psitp.nPatches(); ++i)
            Psi.addPatch(gsTHBSpline<3>( dynamic_cast<const gsTensorBSpline<3>&>(Psitp.patch(i)) ));
        Psi.addAutoBoundaries();
        Psi.computeTopology();

        //Psi.uniformRefine();
        gsMultiBasis<> dbasis(Psi, true);//true: poly-splines (not NURBS)

        geometryMap PPF = A.getMap(Psi);
        auto ff_TG      = A.getCoeff(density, PPF);
        // --------------- adaptive refinement ---------------
        // Specify cell-marking strategy...
        MarkingStrategy adaptRefCrit = PUCA;
        //MarkingStrategy adaptRefCrit = GARU;
        //MarkingStrategy adaptRefCrit = errorFraction;
        real_t adaptRefParam = 0.75;
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