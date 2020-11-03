/** @file poisson_exampleThb.cpp

    @brief Example for using the gsPoissonSolver with adaptive refinement with THB-splines.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): S. Kleiss
*/

#include <gismo.h>
#include <gsAssembler/gsAdaptiveRefUtils.h>

using namespace gismo;

//S.Kleiss
//
//This is a test example for a illustrating the adaptive
//refinement procedure implemented in the gsPoissonAssembler
//
//Flags, parameters, geometry and prescribed exact solution
//are specified within the main() function

int main(int argc, char *argv[])
{
    // Number of initial uniform mesh refinements
    index_t initUnifRef;
    // Number of adaptive refinement loops
    index_t RefineLoopMax;

    // Flag for refinemet criterion
    // (see doxygen documentation of the free function
    // gsMarkElementsForRef explanation)
    index_t refCriterion;   // MarkingStrategy
    // Parameter for computing adaptive refinement threshold
    // (see doxygen documentation of the free function
    // gsMarkElementsForRef explanation)
    real_t refParameter;  // ...specified below with the examples

    // Degree to use for discretization
    index_t degree;

    // Flag whether final mesh should be plotted in ParaView
    bool plot = false;
    bool dump;

    RefineLoopMax = 2;
    initUnifRef = 2;
    degree = 2;
    refCriterion = PUCA;
    refParameter = 0.85;
    dump = false;

    gsCmdLine cmd("Solving a PDE with adaptive refinement using THB-splines.");
    cmd.addSwitch("plot", "Plot resulting mesh in ParaView", plot);
    cmd.addInt("r", "refine", "Maximum number of adaptive refinement steps to perform",
            RefineLoopMax);
    cmd.addInt("i", "initial-ref", "Initial number of uniform refinement steps to perform",
            initUnifRef);
    cmd.addInt("", "degree", "Spline degree of the THB basis", degree);
    cmd.addInt("c", "criterion",  "Criterion to be used for adaptive refinement (1-3, see documentation)",
            refCriterion);
    cmd.addReal("p", "parameter", "Parameter for adaptive refinement", refParameter);
    cmd.addSwitch("dump", "Write geometry and sequence of bases into XML files",
                dump);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    // ****** Prepared test examples ******
    //
    // f       ... source term
    // g       ... exact solution
    // patches ... the computational domain given as object of gsMultiPatch
    //


    // ------ Example 1 ------

    // --- Unit square, with a spike of the source function at (0.25, 0.6)
    gsFunctionExpr<>  f("if( (x-0.25)^2 + (y-0.6)^2 < 0.2^2, 1, 0 )",2);
    //gsFunctionExpr<>  f("if( (x-0.25)^2 + (y-1.6)^2 < 0.2^2, 1, 0 )",2);
    gsFunctionExpr<>  g("0",2);
    gsMultiPatch<> patches( *gsNurbsCreator<>::BSplineRectangle(0.0,0.0,2.0,1.0) );
    //gsMultiPatch<> patches( *gsNurbsCreator<>::BSplineFatQuarterAnnulus(1.0, 2.0) );

    //RefineLoopMax = 6;
    //refParameter = 0.6;

    // ^^^^^^ Example 1 ^^^^^^

    /*
    // ------ Example 2 ------

    // The classical example associated with the L-Shaped domain.
    // The exact solution has a singularity at the re-entrant corner at the origin.

    gsFunctionExpr<>  f("0",2);
    gsFunctionExpr<>  g("if( y>0, ( (x^2+y^2)^(1.0/3.0) )*sin( (2*atan2(y,x) - pi)/3.0 ), ( (x^2+y^2)^(1.0/3.0) )*sin( (2*atan2(y,x) +3*pi)/3.0 ) )",2);

    //gsMultiPatch<> patches (*gsNurbsCreator<>::BSplineLShape_p2C0());
    gsMultiPatch<> patches (*gsNurbsCreator<>::BSplineLShape_p1());

    //RefineLoopMax = 8;
    //refParameter = 0.85;

    // ^^^^^^ Example 2 ^^^^^^
    */

    gsInfo<<"Source function "<< f << "\n";
    gsInfo<<"Exact solution "<< g <<".\n" << "\n";

    // Define Boundary conditions
    gsBoundaryConditions<> bcInfo;
    // Dirichlet BCs
    bcInfo.addCondition( boundary::west,  condition_type::dirichlet, &g );
    bcInfo.addCondition( boundary::east,  condition_type::dirichlet, &g );
    bcInfo.addCondition( boundary::north, condition_type::dirichlet, &g );
    bcInfo.addCondition( boundary::south, condition_type::dirichlet, &g );

    gsTensorBSpline<2,real_t> * geo = dynamic_cast< gsTensorBSpline<2,real_t> * >( & patches.patch(0) );
    gsInfo << " --- Geometry:\n" << *geo << "\n";
    gsInfo << "Number of patches: " << patches.nPatches() << "\n";

    if (dump)
        gsWrite(*geo, "adapt_geo.xml");

    gsTensorBSplineBasis<2,real_t> tbb = geo->basis();
    tbb.setDegree(degree);

    gsInfo << "\nCoarse discretization basis:\n" << tbb << "\n";

    // With this gsTensorBSplineBasis, it's possible to call the THB-Spline constructor
    gsTHBSplineBasis<2,real_t> THB( tbb );

    // Finally, create a vector (of length one) of this gsTHBSplineBasis
    gsMultiBasis<real_t> bases(THB);

    for (int i = 0; i < initUnifRef; ++i)
        bases.uniformRefine();

    gsMultiPatch<> mpsol; // holds computed solution
    gsPoissonAssembler<real_t> pa(patches,bases,bcInfo,f);// constructs matrix and rhs
    pa.options().setInt("DirichletValues", dirichlet::l2Projection);

    if (dump)
        gsWrite(bases[0], "adapt_basis_0.xml");

    // So, ready to start the adaptive refinement loop:
    for( int RefineLoop = 1; RefineLoop <= RefineLoopMax ; RefineLoop++ )
    {
        gsInfo << "\n ====== Loop " << RefineLoop << " of " << RefineLoopMax << " ======" << "\n" << "\n";

        gsInfo <<"Basis: "<< pa.multiBasis() <<"\n";

        // Assemble matrix and rhs
        gsInfo << "Assembling... " << std::flush;
        pa.assemble();
        gsInfo << "done." << "\n";

        // Solve system
        gsInfo << "Solving... " << std::flush;
        gsMatrix<> solVector = gsSparseSolver<>::CGDiagonal(pa.matrix() ).solve( pa.rhs() );
        gsInfo << "done." << "\n";

        // Construct the solution for plotting the mesh later
        pa.constructSolution(solVector, mpsol);
        gsField<> sol(pa.patches(), mpsol);

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

        // Mark elements for refinement, based on the computed local errors and
        // refCriterion and refParameter.
        std::vector<bool> elMarked( elErrEst.size() );
        gsMarkElementsForRef( elErrEst, refCriterion, refParameter, elMarked);

        gsInfo <<"Marked "<< std::count(elMarked.begin(), elMarked.end(), true);

        // Refine the elements of the mesh, based on elMarked.
        gsRefineMarkedElements( pa.multiBasis(), elMarked);

        // Refresh the assembler, since basis is now changed
        pa.refresh();

        if (dump)
        {
            std::stringstream ss;
            ss << "adapt_basis_" << RefineLoop << ".xml";
            gsWrite(bases[0], ss.str());
        }

        if ( (RefineLoop == RefineLoopMax) && plot)
        {
            // Write approximate solution to paraview files
            gsInfo<<"Plotting in Paraview...\n";
            gsWriteParaview<>(sol, "p2d_adaRef_sol", 5001, true);
            // Run paraview and plot the last mesh
            gsFileManager::open("p2d_adaRef_sol.pvd");
        }

    }

    gsInfo << "\nFinal basis: " << bases[0] << "\n";

    return EXIT_SUCCESS;
}
