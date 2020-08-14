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
    gsFunctionExpr<> f("0", 2);

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


    //gsFunction<> *boundCond;
    //gsFunction<> test(1);
    //boundCond = &test;
    gsBoundaryConditions<> bc;
    bc.addCondition(boundary::north, condition_type::dirichlet, 0, 0);
    bc.addCondition(boundary::east, condition_type::dirichlet, 0, 0);
    bc.addCondition(boundary::south, condition_type::dirichlet, 0, 0);
    bc.addCondition(boundary::west, condition_type::dirichlet, 0, 0);
    bc.addCondition(boundary::north, condition_type::dirichlet, 0, 1);
    bc.addCondition(boundary::east, condition_type::dirichlet, 0, 1);
    bc.addCondition(boundary::south, condition_type::dirichlet, 0, 1);
    bc.addCondition(boundary::west, condition_type::dirichlet, 0, 1);
    gsInfo << "Boundary conditions:\n" << bc << "\n";

    //! [Refinement]

    //! [Problem setup]
    gsExprAssembler<> A(2, 2);
    gsExprAssembler<> EpsA(1,1);
    gsOptionList opt = gsAssembler<>::defaultOptions();
    //opt.setInt("DirichletValues", dirichlet::l2projection);
    opt.setInt("DirichletStrategy", dirichlet::nitsche);
    opt.setReal("quA", 1.0);
    opt.setInt("quB", 1);
    opt.setReal("bdA", 2.0);
    opt.setInt("bdB", 1);
    //gsInfo << "Assembler " << opt;
    A.setOptions(opt);
    EpsA.setOptions(opt);

    gsInfo << "Active options:\n" << A.options() << "\n";
    typedef gsExprAssembler<>::geometryMap geometryMap;
    typedef gsExprAssembler<>::variable    variable;
    typedef gsExprAssembler<>::space       space;
    typedef gsExprAssembler<>::solution    solution;

    // Elements used for numerical integration
    A.setIntegrationElements(dbasis);
    gsExprEvaluator<> ev(A);
    EpsA.setIntegrationElements(dbasis);
    gsExprEvaluator<> evEps(EpsA);


    // Set the geometry map
    geometryMap G = A.getMap(patch);
    geometryMap EpsG = EpsA.getMap(patch);

    // Set the discretization space
    space u = A.getSpace(dbasis,1,0);
    space v = A.getSpace(dbasis,1,1);
    space eps = EpsA.getSpace(dbasis);

    gsFunctionExpr<> firstComp("1","0",2);
    variable firstCoeff = A.getCoeff(firstComp, G);
    gsInfo << "firstcomp: " << firstComp << "\n";
    gsFunctionExpr<> secondComp("0","1",2);
    variable secondCoeff = A.getCoeff(secondComp, G);
    gsInfo << "secondcomp: " << secondComp << "\n";

    EpsA.initSystem();
    //EpsA.assemble(eps * eps.tr() * meas(EpsG));
    //gsMatrix<> MassEps = EpsA.matrix();
    //gsInfo << "MassEps: \n" << MassEps << "\n";
    EpsA.initSystem();
    EpsA.assemble(eps * eps.tr() * (nv(EpsG)*firstCoeff) * nv(EpsG).norm());

    u.setInterfaceCont(0); // todo: 1 (smooth basis)
    u.addBc(bc.get("Dirichlet",0)); // (!) must be called only once
    v.setInterfaceCont(0); // todo: 1 (smooth basis)
    v.addBc(bc.get("Dirichlet",1)); // (!) must be called only once

    variable ff = A.getCoeff(f, G);

    // Solution vector and solution variable


    //Sxx
    A.initSystem();
    real_t a = 2.0;
    A.assemble((igrad(u,G)*firstCoeff)*(igrad(u,G)*firstCoeff).tr()*meas(G));
    //A.assemble((igrad(u,G)*firstCoeff)*(igrad(u,G)*firstCoeff).tr()*meas(G));//will sum this to precedent values
    gsMatrix<> Temp;
    gsMatrix<> Temp2;
    Temp = A.matrix();
    //gsInfo << "Temp: " << Temp.toDense() << "\n";
    //Syy
    A.initSystem();
    A.assemble((igrad(u,G)*secondCoeff)*(igrad(u,G)*secondCoeff).tr()*meas(G));
    Temp2 = A.matrix();

    A.initSystem();
    A.assemble(igrad(u,G)*igrad(u,G).tr()*meas(G));
    gsInfo << "Diff: " << (Temp+Temp2-A.matrix().toDense()).norm() << "\n";


    //Ix
    A.initSystem();
    A.assemble(u*(igrad(u,G)*firstCoeff).tr()*meas(G));////this should be transpose of Ix
    Temp = A.matrix();
    //gsInfo << "Temp: " << Temp << "\n";
    A.initSystem();
    A.assemble((igrad(u,G)*firstCoeff)*u.tr()*meas(G));//should be this
    Temp2 = A.matrix();
    //gsInfo << "Temp2: " << Temp2 << "\n";
    gsInfo << "Norm of difference: "<< (Temp-Temp2.transpose()).norm() << "\n";


    //Iy
    A.initSystem();
    A.assemble(u*(igrad(u,G)*secondCoeff).tr()*meas(G));//this should be transpose of Iy
    Temp = A.matrix();
    //gsInfo << "Temp: " << Temp << "\n";
    A.initSystem();
    A.assemble((igrad(u,G)*secondCoeff)*u.tr()*meas(G));//should be this
    Temp2 = A.matrix();
    //gsInfo << "Temp2: " << Temp2 << "\n";
    gsInfo << "Norm of difference: "<< (Temp-Temp2.transpose()).norm() << "\n";


    gsSparseSolver<>::CGDiagonal solver;
    
    gsInfo << "Number of degrees of freedom: " << A.numDofs() << "\n";

    A.assemble(u * u.tr() * meas(G));
    A.assemble(v * v.tr() * meas(G));
    A.assembleRhsBc((igrad(u,G)*secondCoeff) * u.tr() * nv(G).norm(), bc.dirichletSides() );

    gsInfo << "RHS:\n" << A.rhs() << "\n";
    gsSparseMatrix<> Mass;
    Mass = A.matrix();
    //gsInfo << "M:\n" << Mass.toDense() << "\n";

    int ndof = A.numDofs();
    gsMatrix<> solVector;
    solution u_sol = A.getSolution(u, solVector);
    solution v_sol = A.getSolution(v, solVector);

    solVector.setRandom(ndof,1);
    gsInfo << "solVector: \n" << solVector << "\n";

    if (plot){
    gsInfo << "Plotting in Paraview...\n";

    
    ev.options().setSwitch("plot.elements", true);
    ev.writeParaview(u_sol, G, "solutionU");
    ev.writeParaview(v_sol, G, "solutionV");
    //ev.writeParaview( u_ex    , G, "solution_ex");
    //ev.writeParaview( u, G, "aa");
    //gsFileManager::open("solution.pvd");
    }
    return EXIT_SUCCESS;

}// end main
