/** @file ViscoEl.cpp

    @Replicating viscoelasticity results

    Author(s): A. Barion
*/

//! [Include namespace]
#include <gismo.h>

using namespace gismo;
//! [Include namespace]

int main(int argc, char* argv[])
{
    //! [Parse command line]
    bool plot = false;
    index_t numRefine = 2;
    index_t numElevate = 0;
    bool last = false;

    real_t endTime = 1.0;
    int numSteps = 100;
    real_t dt = endTime / numSteps;

    gsCmdLine cmd("Viscoelasticity on unit square problem.");
    cmd.addInt("e", "degreeElevation", "Number of degree elevation steps to perform before solving (0: equalize degree in all directions)", numElevate);
    cmd.addInt("r", "uniformRefine", "Number of Uniform h-refinement steps to perform before solving", numRefine);
    cmd.addSwitch("last", "Solve solely for the last level of h-refinement", last);
    cmd.addSwitch("plot", "Create a ParaView visualization file with the solution", plot);

    try { cmd.getValues(argc, argv); }
    catch (int rv) { return rv; }
    //! [Parse command line]

    //! [Read input file]

    // Define Geometry, must be a gsMultiPatch object
    gsMultiPatch<> patch;

    // Create one patch as the unit square [0,2]^2
    patch = gsNurbsCreator<>::BSplineSquareGrid(1, 1, 2);
    gsInfo << "The domain is a " << patch << "\n";

    //Source function for velocities
    //gsFunctionExpr<> f("if ( sqrt((x-6)*(x-6)+(y-6)*(y-6))<=1,-(x-6),0)", "if ( sqrt((x-6)*(x-6)+(y-6)*(y-6)) <= 1, -(y-6), 0)", 2);
    gsFunctionExpr<> f("0.05", 2);

    gsInfo << "Source function for heat reads: " << f << "\n";

    // Assembler options
    //gsOptionList opt = gsAssembler<>::defaultOptions();
    //opt.setInt("DirichletValues", 101);
    //opt.setInt("DirichletStrategy", dirichlet::nitsche);
    //opt.setReal("quA", 1.0);
    //opt.setInt("quB", 1);
    //opt.setReal("bdA", 2.0);
    //opt.setInt("bdB", 1);
    //opt.setInt("quRule", 1);
    //gsInfo << "Assembler " << opt << "\n";


    //! [Refinement]

    // p-refine
    // Elevate and p-refine the basis to order p + numElevate
    // where p is the highest degree in the bases
    if (numElevate != 0)
        patch.degreeElevate(numElevate);

    // h-refine
    for (int r = 0; r < numRefine; ++r)
        patch.uniformRefine();

    gsMultiBasis<> dbasis(patch);

    gsWriteParaview(patch, "mp", 1000, true);

    gsInfo << "Patches: " << patch.nPatches() << ", degree: " << dbasis.minCwiseDegree() << "\n";
    gsInfo << dbasis.basis(0) << "\n";

    gsBoundaryConditions<> bc;
    bc.addCondition(boundary::north, condition_type::dirichlet, 0, 0);
    bc.addCondition(boundary::east, condition_type::dirichlet, 0, 0);
    bc.addCondition(boundary::south, condition_type::dirichlet, 0, 0);
    bc.addCondition(boundary::west, condition_type::dirichlet, 0, 0);
    gsInfo << "Boundary conditions:\n" << bc << "\n";

    //! [Refinement]

    //! [Problem setup]
    gsExprAssembler<> A(1, 1);
    gsOptionList opt = gsAssembler<>::defaultOptions();
    //opt.setInt("DirichletValues", dirichlet::l2Projection);
    //opt.setInt("DirichletStrategy", dirichlet::nitsche);
    opt.setReal("quA", 1.0);
    opt.setInt("quB", 1);
    opt.setReal("bdA", 2.0);
    opt.setInt("bdB", 1);
    //gsInfo << "Assembler " << opt;
    A.setOptions(opt);

    gsInfo << "Active options:\n" << A.options() << "\n";
    typedef gsExprAssembler<>::geometryMap geometryMap;
    typedef gsExprAssembler<>::variable    variable;
    typedef gsExprAssembler<>::space       space;
    typedef gsExprAssembler<>::solution    solution;

    // Elements used for numerical integration
    A.setIntegrationElements(dbasis);
    gsExprEvaluator<> ev(A);


    // Set the geometry map
    geometryMap G = A.getMap(patch);

    // Set the discretization space
    space u = A.getSpace(dbasis);

    u.setInterfaceCont(0); // todo: 1 (smooth basis)
    u.addBc(bc.get("Dirichlet")); // (!) must be called only once


    variable ff = A.getCoeff(f, G);

    // Solution vector and solution variable
    
    A.initSystem();

    gsSparseSolver<>::CGDiagonal solver;
    
    gsInfo << "Number of degrees of freedom: " << A.numDofs() << "\n";

    A.assemble(u * u.tr() * meas(G));
    gsSparseMatrix<> Mass;
    Mass = A.matrix();
    gsInfo << "M:\n" << Mass << "\n";

    int ndof = A.numDofs();
    gsMatrix<> solVector;
    solution u_sol = A.getSolution(u, solVector);

    A.initSystem();
    gsInfo << "test\n";
    A.assemble(igrad(u, G) * igrad(u, G).tr() * meas(G), u * ff * meas(G));
    gsSparseMatrix<> S;
    S = A.matrix();

    gsMatrix<> RHS;
    RHS = A.rhs();

    gsMatrix<> F;

    gsMatrix<> auxVec;
    gsMatrix<> temp;
    
    auxVec.setZero(ndof, 1);
    solVector.setZero(ndof, 1);


    gsInfo << "Sum of Mass*auxVec:" << (Mass * auxVec).sum() << "\n";
    gsInfo << "auxVec:\n" << auxVec << "\n";
    gsSparseMatrix<> solveMat;
    solveMat = Mass + dt * S;
    solver.compute(solveMat);
    gsInfo << "computed S:\n" << S << "\n";
    gsInfo << "computed Mass + dt * S:\n" << solveMat << "\n";

    for( index_t it=0; it<numSteps; ++it){
        F = Mass * solVector + dt * RHS;
        //gsInfo << "F:\n" << F << "\n";
        solVector = solver.solve(F);
        gsInfo << "solved timestep: " << it << "\n";
        //gsInfo << "auxVec:\n" << auxVec << "\n";
        //gsInfo << "temp to assign:\n" << temp << "\n";
        //auxVec = temp;
        //gsInfo << "auxVec:\n" << auxVec << "\n";
        //gsInfo << it << "\n";
    }

    //solVector = auxVec;

    gsInfo << "auxVec:\n" << auxVec << "\n";
    gsInfo << "F:\n" << F << "\n";
    gsInfo << "solVector:\n" << solVector << "\n";

    if (plot){
    gsInfo << "Plotting in Paraview...\n";

    
    ev.options().setSwitch("plot.elements", true);
    ev.writeParaview(u_sol, G, "solution");
    //ev.writeParaview( u_ex    , G, "solution_ex");
    //ev.writeParaview( u, G, "aa");
    //gsFileManager::open("solution.pvd");
    }
    return EXIT_SUCCESS;

}// end main
