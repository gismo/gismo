/** @file poissonHB_refinement.cpp

    @brief Example for refinement the HB spline in a different way.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): P. Weinmüller
*/

#include <gismo.h>
#include <gsAssembler/gsAdaptiveRefUtils.h>

using namespace gismo;

int main(int argc, char *argv[])
{
    // Number of initial uniform mesh refinements
    index_t initUnifRef;
    // Number of adaptive refinement loops
    index_t RefineLoopMax;

    // Degree to use for discretization
    index_t degree;

    // Flag whether final mesh should be plotted in ParaView
    bool plot = false;

    RefineLoopMax = 2;
    initUnifRef = 2;
    degree = 2;

    gsCmdLine cmd("Solving a PDE with adaptive refinement using THB-splines.");
    cmd.addSwitch("plot", "Plot resulting mesh in ParaView", plot);
    cmd.addInt("k", "refine", "Maximum number of adaptive refinement steps to perform",
               RefineLoopMax);
    cmd.addInt("i", "initial-ref", "Initial number of uniform refinement steps to perform",
               initUnifRef);
    cmd.addInt("p", "degree", "Spline degree of the THB basis", degree);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    // ****** Prepared test examples ******
    //
    // f       ... source term
    // g       ... exact solution
    // patches ... the computational domain given as object of gsMultiPatch
    //


    //! [Function data]
    // Define source function
    gsFunctionExpr<> f ("16*pi*pi*(2*cos(4*pi*x)*cos(4*pi*y) - cos(4*pi*x) - cos(4*pi*y))",2);
    gsFunctionExpr<> g("(1 - cos(4*pi*x)) * (1 - cos(4*pi*y))",2);

    // Print out source function and solution
    gsInfo<<"Source function "<< f << "\n";
    gsInfo<<"Exact solution "<< g <<"\n\n";
    //! [Function data]

    gsFileData<> fd("planar/quarterCircle.xml");
    gsInfo << "Loaded file "<< fd.lastPath() <<"\n";

    gsMultiPatch<> multiPatch;
    fd.getId(0, multiPatch); // id=0: Multipatch domain
    multiPatch.computeTopology();

    multiPatch.degreeElevate(degree);

    for (int i = 0; i < initUnifRef; ++i)
        multiPatch.uniformRefine();

    gsWriteParaview(multiPatch,"geometry",5000,true);


    gsInfo<<"Source function "<< f << "\n";
    gsInfo<<"Exact solution "<< g <<".\n" << "\n";

    // Define Boundary conditions
    gsBoundaryConditions<> bcInfo;
    // Dirichlet BCs
    bcInfo.addCondition( boundary::west,  condition_type::dirichlet, &g );
    bcInfo.addCondition( boundary::east,  condition_type::dirichlet, &g );
    bcInfo.addCondition( boundary::north, condition_type::dirichlet, &g );
    bcInfo.addCondition( boundary::south, condition_type::dirichlet, &g );

    gsTensorBSpline<2,real_t> * geo = dynamic_cast< gsTensorBSpline<2,real_t> * >( & multiPatch.patch(0) );
    gsInfo << " --- Geometry:\n" << *geo << "\n";
    gsInfo << "Number of patches: " << multiPatch.nPatches() << "\n";

    gsTensorBSplineBasis<2,real_t> tbb = geo->basis();

    gsInfo << "\nCoarse discretization basis:\n" << tbb << "\n";

    // With this gsTensorBSplineBasis, it's possible to call the THB-Spline constructor
    //gsTHBSplineBasis<2,real_t> THB( tbb );
    gsHBSplineBasis<2,real_t> HB( tbb );

    // Finally, create a vector (of length one) of this gsTHBSplineBasis
    gsMultiBasis<real_t> bases(HB);

    gsMultiPatch<> mpsol; // holds computed solution
    gsPoissonAssembler<real_t> pa(multiPatch,bases,bcInfo,f);// constructs matrix and rhs
    pa.options().setInt("DirichletValues", dirichlet::l2Projection);


    gsMatrix<> mat(4,4);
    mat.setOnes();
    if (mat(1,1) * mat(1,0) > 10e-15)
        mat(1,1) = 0;

    real_t h_mesh = 1;

    gsVector<> l2err(RefineLoopMax+1), h1err(RefineLoopMax+1);

    // So, ready to start the adaptive refinement loop:
    for( int RefineLoop = 1; RefineLoop <= RefineLoopMax ; RefineLoop++ )
    {
        h_mesh = 0.5 * h_mesh;
        gsInfo << "\n ====== Loop " << RefineLoop << " of " << RefineLoopMax << " ======" << "\n" << "\n";
        gsInfo << "\n ====== Mesh size = " << h_mesh << "\n" << "\n";

        gsHBSplineBasis<2, real_t> hbSplineBasis = dynamic_cast<gsHBSplineBasis<2, real_t> & > (pa.multiBasis().basis(0));
        std::vector< gsSortedVector< unsigned > > Xmatrix =  hbSplineBasis.getXmatrix();
        for (index_t i = 0; i< Xmatrix.size(); i++)
            for (index_t j = 0; j < Xmatrix[i].size(); j++)
                gsInfo <<"Basis: "<< Xmatrix[i][j] <<"\n";

        // Assemble matrix and rhs
        gsInfo << "Assembling... " << std::flush;
        pa.assemble();
        gsInfo << "done." << "\n";

        gsInfo << "DIM: " << pa.matrix().dim().first << " : " << pa.matrix().dim().second << "\n";

        // Solve system
        gsInfo << "Solving... " << std::flush;
        gsMatrix<> solVector = gsSparseSolver<>::CGDiagonal(pa.matrix() ).solve( pa.rhs() );
        gsInfo << "done." << "\n";

        // Construct the solution for plotting the mesh later
        pa.constructSolution(solVector, mpsol);
        gsField<> sol(pa.patches(), mpsol);

        /*
        // Set up and compute the L2-error to the known exact solution...
        gsExprEvaluator<> ev;
        ev.setIntegrationElements(pa.multiBasis());
        gsExprEvaluator<>::geometryMap Gm = ev.getMap(patches);
        gsExprEvaluator<>::variable f1 = ev.getVariable(mpsol);
        gsExprEvaluator<>::variable ff = ev.getVariable(f, Gm);

        // The vector with element-wise local error estimates.
        ev.integralElWise( (ilapl(f1,Gm) + ff).sqNorm() * meas(Gm) );
        const std::vector<real_t> & elErrEst = ev.elementwise();

        // Get the vector with element-wise local (known in this case) errors...
        //gsExprEvaluator<>::variable gg = ev.getVariable(g, Gm);
        //ev.integralElWise( (f1 - gg).sqNorm() * meas(Gm) );
        //const std::vector<real_t> & elErrEst = ev.elementwise();
*/
        // Mark elements for refinement, based on the computed local errors and
        // refCriterion and refParameter.
        std::vector<bool> elMarked( pa.multiBasis().basis(0).numElements() );
        //gsMarkElementsForRef( elErrEst, refCriterion, refParameter, elMarked);

        gsInfo << "Elements : " <<  pa.multiBasis().basis(0).numElements() << "\n";
        //elMarked.at(pa.multiBasis().basis(0).numElements()-1) = true;
        elMarked.at(0) = true;

        // for all elements in patch pn
/*        index_t globalCounter = 0;
        typename gsBasis<real_t>::domainIter domIt = pa.multiBasis().basis(0).makeDomainIterator();

        for (; domIt->good(); domIt->next())
        {
            gsMapData<real_t> mapData;
            mapData.flags = NEED_MEASURE | NEED_VALUE;
            gsMatrix<> points(2,2);
            points.col(0) = domIt->lowerCorner();
            points.col(1) = domIt->upperCorner();

            mapData.points = points;

            multiPatch.patch(0).computeMap(mapData);

            gsMatrix<> cellSize = mapData.eval(1) - mapData.eval(0);
            gsInfo << "Point : " << cellSize  << " : " << (domIt->getMinCellLength() * mapData.measures) [0] <<  "\n";

            //if (cellSize.array().abs().maxCoeff() > (h_mesh - 10e-6) )
            if (cellSize[1] > (h_mesh - 10e-6) )
            //if ((domIt->getMinCellLength() * mapData.measures) [0] > (h_mesh - 10e-6) )
            {

                elMarked.at(globalCounter) = true;
                //gsInfo << globalCounter << "\n";
            }

            globalCounter++;
        }
*/
        gsInfo <<"Marked "<< std::count(elMarked.begin(), elMarked.end(), true) << "\n";

        // Refine the elements of the mesh, based on elMarked.
        gsRefineMarkedElements( pa.multiBasis(), elMarked);

        //Contruct the H2 norm, part by part.
        real_t errorH1Semi = sol.distanceH1(g, false);
        real_t errorL2 = sol.distanceL2(g, false);
        real_t errorH1 = math::sqrt(errorH1Semi*errorH1Semi + errorL2*errorL2);

        gsInfo << "The L2 error of the solution is : " << errorL2 << "\n";
        gsInfo << "The H1 error of the solution is : " << errorH1 << "\n";

        l2err[RefineLoop-1] = errorL2;
        h1err[RefineLoop-1] = errorH1;

        // Refresh the assembler, since basis is now changed
        pa.refresh();

        gsWriteParaview(pa.multiBasis().basis(0),"ParameterSpace",5000,true);

        if ( (RefineLoop == RefineLoopMax) && plot)
        {
            // Write approximate solution to paraview files
            gsInfo<<"Plotting in Paraview...\n";
            gsWriteParaview<>(sol, "p2d_adaRef_sol", 5001, true);
            // Run paraview and plot the last mesh
            //gsFileManager::open("p2d_adaRef_sol.pvd");
        }

    }


    //! [Error and convergence rates]
    gsInfo<< "\n\nL2 error: "<<std::scientific<<std::setprecision(3)<<l2err.transpose()<<"\n";
    gsInfo<< "H1 error: "<<std::scientific<<h1err.transpose()<<"\n";

    if (RefineLoopMax>0)
    {
        gsInfo<< "\nEoC (L2): " << std::fixed<<std::setprecision(2)
              << ( l2err.head(RefineLoopMax).array() /
                  l2err.tail(RefineLoopMax).array() ).log().transpose() / std::log(2.0) <<"\n";

        gsInfo<<   "EoC (H1): "<< std::fixed<<std::setprecision(2)
              <<( h1err.head(RefineLoopMax).array() /
                  h1err.tail(RefineLoopMax).array() ).log().transpose() / std::log(2.0) <<"\n";
    }
    //! [Error and convergence rates]

    gsInfo << "\nFinal basis: " << pa.multiBasis().basis(0) <<"\n";

    return EXIT_SUCCESS;
}
