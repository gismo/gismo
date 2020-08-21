/** @file MyViscoEl.cpp

    @Replicating viscoelasticity results

    Author(s): A. Barion
*/

//! [Include namespace]
#include <gismo.h>
#include <chrono>

using namespace gismo;
//! [Include namespace]

int main(int argc, char* argv[])
{
    //! [Parse command line]
    bool plot = false;
    index_t numRefine = 2;
    index_t numElevate = 0;
    bool last = false;

    real_t endTime = 5.0;
    index_t numSteps = 50;
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
    real_t nu = 0.48;
    real_t E = 31.0;
    real_t mu1 = 100.0;
    real_t mu2 = 100.0;
    real_t rho = 1.02;
    real_t G_tilde = E*sqrt(rho)/(1+nu);

    // Define Geometry, must be a gsMultiPatch object
    gsMultiPatch<> patch;

    // Create one patch as the unit square [0,12]^2
    patch = gsNurbsCreator<>::BSplineSquareGrid(1, 1, 12);
    gsInfo << "The domain is a " << patch << "\n";

    //Source function for velocities
    gsFunctionExpr<> f1("if ( (x<=7.5) and (x>=4.5) and (y<=7.5) and (y>=4.5), -sgn(x-6), 0)", 2);
    gsFunctionExpr<> f2("0", 2);
    //gsFunctionExpr<> f("0", 2);
    

    gsInfo << "Source function for first velocity component: " << f1 << "\n";
    gsInfo << "Source function for second velocity component: " << f2 << "\n";

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


    gsConstantFunction<> fun(1,2);
    // gsFunctionExpr<> fun("sin(x)",2);

    gsBoundaryConditions<> bc;
    // bc.addCondition(boundary::north, condition_type::dirichlet, &fun, 0);
    bc.addCondition(boundary::north, condition_type::dirichlet, 0, 3);
    bc.addCondition(boundary::east, condition_type::dirichlet, 0, 3);
    bc.addCondition(boundary::south, condition_type::dirichlet, 0, 3);
    bc.addCondition(boundary::west, condition_type::dirichlet, 0, 3);
    bc.addCondition(boundary::north, condition_type::dirichlet, 0, 4);
    bc.addCondition(boundary::east, condition_type::dirichlet, 0, 4);
    bc.addCondition(boundary::south, condition_type::dirichlet, 0, 4);
    bc.addCondition(boundary::west, condition_type::dirichlet, 0, 4);
    gsInfo << "Boundary conditions:\n" << bc << "\n";


    //! [Refinement]

    //! [Problem setup]
    gsExprAssembler<> A(5, 5);                              // velocity assembler

    gsOptionList opt = gsAssembler<>::defaultOptions();
    //opt.setInt("DirichletValues", dirichlet::l2projection);
    opt.setInt("DirichletStrategy", dirichlet::nitsche);
    opt.setReal("quA", 1.0);
    opt.setInt("quB", 1);
    opt.setReal("bdA", 2.0);
    opt.setInt("bdB", 1);
    //gsInfo << "Assembler " << opt;
    A.setOptions(opt);
    

    //gsInfo << "Active options:\n" << A.options() << "\n";
    typedef gsExprAssembler<>::geometryMap geometryMap;
    typedef gsExprAssembler<>::variable    variable;
    typedef gsExprAssembler<>::space       space;
    typedef gsExprAssembler<>::solution    solution;

    //typedef typename gsExprHelper<real_t>::nullExpr    nullExpr;


    // Elements used for numerical integration
    A.setIntegrationElements(dbasis);
    gsExprEvaluator<> ev(A);

    // Set the geometry map
    geometryMap G = A.getMap(patch);

    // Set the discretization space
    space eps11 = A.getSpace(dbasis,1,0);
    space eps22 = A.getSpace(dbasis,1,1);
    space eps12 = A.getSpace(dbasis,1,2);
    space u = A.getSpace(dbasis,1,3);
    space v = A.getSpace(dbasis,1,4);

    u.setInterfaceCont(0); // todo: 1 (smooth basis)
    u.addBc(bc.get("Dirichlet",3)); // (!) must be called only once
    v.setInterfaceCont(0); // todo: 1 (smooth basis)
    v.addBc(bc.get("Dirichlet",4)); // (!) must be called only once

    gsFunctionExpr<> firstComp("1","0",2);
    variable firstCoeff = A.getCoeff(firstComp, G);
    gsInfo << "firstcomp: " << firstComp << "\n";
    gsFunctionExpr<> secondComp("0","1",2);
    variable secondCoeff = A.getCoeff(secondComp, G);
    gsInfo << "secondcomp: " << secondComp << "\n";

    //EpsA.assemble(eps * eps.tr() * meas(EpsG));
    //gsMatrix<> MassEps = EpsA.matrix();
    //gsInfo << "MassEps: \n" << MassEps << "\n";


    gsVector<> pt(2);
    pt<<0.25,0.25;
    // pt.setConstant(0.25);
    
    variable ff1 = A.getCoeff(f1, G);
    variable ff2 = A.getCoeff(f2, G);

    // Solution vector and solution variable

    //Sxx_u
    
    // A.assemble(Sxx);
    // //A.assemble((igrad(u,G)*firstCoeff)*(igrad(u,G)*firstCoeff).tr()*meas(G));//will sum this to precedent values
    // gsMatrix<> Temp;
    // gsMatrix<> Temp2;
    // Temp = A.matrix();
    // //gsInfo << "Temp: " << Temp.toDense() << "\n";
    // //Syy
    // A.initSystem();
    // A.assemble(Syy);
    // Temp2 = A.matrix();

    // A.initSystem();
    // A.assemble(igrad(u,G)*igrad(u,G).tr()*meas(G));
    // gsInfo << "Diff: " << (Temp+Temp2-A.matrix().toDense()).norm() << "\n";

    // gsExprEvaluator<> evA(A);
    // // gsDebug<<"\n"<<evA.eval(Sxx,pt)<<"\n";

    // //Ix
    // A.initSystem();
    // A.assemble(u*(igrad(u,G)*firstCoeff).tr()*meas(G));////this should be transpose of Ix
    // Temp = A.matrix();
    // //gsInfo << "Temp: " << Temp << "\n";
    // A.initSystem();
    // A.assemble((igrad(u,G)*firstCoeff)*u.tr()*meas(G));//should be this
    // Temp2 = A.matrix();
    // //gsInfo << "Temp2: " << Temp2 << "\n";
    // gsInfo << "Norm of difference: "<< (Temp-Temp2.transpose()).norm() << "\n";


    // //Iy
    // A.initSystem();
    // A.assemble(u*(igrad(u,G)*secondCoeff).tr()*meas(G));//this should be transpose of Iy
    // Temp = A.matrix();
    // //gsInfo << "Temp: " << Temp << "\n";
    // A.initSystem();
    // A.assemble((igrad(u,G)*secondCoeff)*u.tr()*meas(G));//should be this
    // Temp2 = A.matrix();
    // //gsInfo << "Temp2: " << Temp2 << "\n";
    // gsInfo << "Norm of difference: "<< (Temp-Temp2.transpose()).norm() << "\n";
    A.initSystem();

    gsSparseSolver<>::BiCGSTABDiagonal solver;

    gsInfo << "Number of degrees of freedom: " << A.numDofs() << "\n";

    index_t ndof = A.numDofs();
    gsMatrix<> solVector;
    solution u_sol = A.getSolution(u, solVector);
    solution v_sol = A.getSolution(v, solVector);
    solution eps11_sol = A.getSolution(eps11, solVector);
    solution eps22_sol = A.getSolution(eps22, solVector);
    solution eps12_sol = A.getSolution(eps12, solVector);

    solVector.setZero(ndof,1);

    gsMatrix<> aux;
    gsVector<> displ1;
    gsVector<> displ2;
    gsMatrix<> prevTimestep;
    gsMatrix<> displContrib;
    gsVector<> displVec;
    displVec.setZero(ndof);
    eps11.getCoeffs(solVector,aux);
    displ1.setZero(aux.size());
    displ2 = displ1;
    index_t displ_dof = displ1.size();
    gsMatrix<> oldMass;
    gsMatrix<> u_aux;
    gsMatrix<> v_aux; 
    u.getCoeffs(solVector,u_aux);
    v.getCoeffs(solVector,v_aux);
    real_t half = 0.5;
    real_t damp;

    gsMatrix<> F;
    gsSparseMatrix<> solveMat;

    auto Sxx_u = (igrad(u,G)*firstCoeff)*(igrad(u,G)*firstCoeff).tr()*meas(G);
    auto Syy_u = (igrad(u,G)*secondCoeff)*(igrad(u,G)*secondCoeff).tr()*meas(G);
    auto Sxx_v = (igrad(v,G)*firstCoeff)*(igrad(v,G)*firstCoeff).tr()*meas(G);
    auto Syy_v = (igrad(v,G)*secondCoeff)*(igrad(v,G)*secondCoeff).tr()*meas(G);

    for(index_t it = 1; it <= 1; ++it){

        auto start = std::chrono::high_resolution_clock::now();
        A.initSystem();
        A.assemble(rho * u * u.tr() * meas(G));
        A.assemble(rho * v * v.tr() * meas(G));
        oldMass = A.matrix();
        prevTimestep = oldMass * solVector;

        u.getCoeffs(solVector,u_aux);
        v.getCoeffs(solVector,v_aux);

        //gsInfo << "u_aux: \n" << u_aux << "\n";

        patch.patch(0).coefs().col(0) += dt*u_aux;
        patch.patch(0).coefs().col(1) += dt*v_aux;

        displ1 = displ1 + dt * u_aux;
        displ2 = displ2 + dt * v_aux;
        displVec.segment(0,displ_dof) = displ1;
        displVec.segment(displ_dof, displ_dof) = displ2;
        
        //gsInfo << "displVec:" << displVec << "\n";
        A.initSystem();
        A.assemble(eps11*(igrad(eps11,G)*firstCoeff).tr()*meas(G));
        
        A.assemble(eps22*(igrad(eps22,G)*secondCoeff).tr()*meas(G));

        A.assemble((half*eps12)*(igrad(eps11,G)*secondCoeff).tr()*meas(G) +
                    (half*eps12)*(igrad(eps22,G)*firstCoeff).tr()*meas(G));
        
        
        displContrib = A.matrix() * displVec;

        //gsInfo << "displContrib:\n" << displContrib << "\n";
        
        A.initSystem();
        damp = 4.2*(1-exp(-4*(dt*it-0.1)/(20-0.1)));
        //gsInfo << "damp: " << damp << "\n";

        //Mass matrices for strains and velocity components
        A.assemble(eps11 * eps11.tr() * meas(G));
        A.assemble(eps22 * eps22.tr() * meas(G));
        A.assemble(eps12 * eps12.tr() * meas(G));
        A.assemble(rho * u * u.tr() * meas(G), damp * u * ff1 * meas(G));
        A.assemble(rho * v * v.tr() * meas(G), damp * v * ff2 * meas(G));

        //Strain contributions to vel comps
        //eps11-->u
        //A.assemble((-dt)*G_tilde*(1-nu)/(1-2*nu)*(u*(nv(G).tr()*firstCoeff)) * eps11.tr() * nv(G).norm());
        A.assemble((dt)*G_tilde*(1-nu)/(1-2*nu)*((igrad(u,G)*firstCoeff))*eps11.tr()*meas(G));
        
        //eps22-->u
        //A.assemble((-dt)*G_tilde*(nu)/(1-2*nu)*(u*(nv(G).tr()*firstCoeff)) * eps22.tr() * nv(G).norm());
        A.assemble((dt)*G_tilde*(nu)/(1-2*nu)*((igrad(u,G)*firstCoeff))*eps22.tr()*meas(G));
        
        //eps12-->u
        //A.assemble((-dt)*G_tilde*(u*(nv(G).tr()*secondCoeff)) * eps12.tr() * nv(G).norm());
        A.assemble((dt)*G_tilde*((igrad(u,G)*secondCoeff))*eps12.tr()*meas(G));

        //eps11-->v
        //A.assemble((-dt)*G_tilde*(nu)/(1-2*nu)*(v*(nv(G).tr()*secondCoeff)) * eps11.tr() * nv(G).norm());
        A.assemble((dt)*G_tilde*(nu)/(1-2*nu)*((igrad(v,G)*secondCoeff))*eps11.tr()*meas(G));

        //eps22-->v
        //A.assemble((-dt)*G_tilde*(1-nu)/(1-2*nu)*(v*(nv(G).tr()*secondCoeff)) * eps22.tr() * nv(G).norm());
        A.assemble((dt)*G_tilde*(1-nu)/(1-2*nu)*((igrad(v,G)*secondCoeff))*eps22.tr()*meas(G));

        //eps12-->v
        //A.assemble((-dt)*G_tilde*(u*(nv(G).tr()*firstCoeff)) * eps12.tr() * nv(G).norm());
        A.assemble((dt)*G_tilde*((igrad(v,G)*firstCoeff))*eps12.tr()*meas(G));

        // gsDebug<<"(firstCoeff) = \n"<<evEps.eval(firstCoeff_Eps,pt)<<"\n";
        //gsDebug<<"igrad(v,G) = \n"<<ev.eval(igrad(v,G),pt)<<"\n";
        //gsDebug<<"(nv(G).tr()*firstCoeff/nv(G).norm()) = \n"<<ev.eval(nv(G).tr()*firstCoeff/nv(G).norm(),pt)<<"\n";


        //Velocity contributions
       
        //u-->u
        A.assemble(dt*(mu1+mu2)*Sxx_u);
        A.assemble(dt*mu1/2*Syy_u);

        //u-->v
        A.assemble(dt*mu2*(igrad(u,G)*firstCoeff)*(igrad(v,G)*secondCoeff).tr()*meas(G));
        A.assemble(dt*mu1/2*(igrad(u,G)*secondCoeff)*(igrad(v,G)*firstCoeff).tr()*meas(G));

        //v-->v
        A.assemble(dt*(mu1+mu2)*Syy_v);
        A.assemble(dt*mu1/2*Sxx_v);

        //v-->u
        A.assemble(dt*mu1/2*(igrad(v,G)*firstCoeff)*(igrad(u,G)*secondCoeff).tr()*meas(G));
        A.assemble(dt*mu2*(igrad(v,G)*secondCoeff)*(igrad(u,G)*firstCoeff).tr()*meas(G));

        auto finish = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = finish - start;
        gsInfo << "Elapsed time: " << elapsed.count() << " s\n";
        //gsInfo << "A.matrix():\n" << A.matrix() << "\n";
        solver.compute(A.matrix());
        F = prevTimestep + displContrib + dt*A.rhs();

        //gsInfo << "F: " << F << "\n";
        solVector = solver.solve(F);
        real_t zero = 0;

        A.initSystem();
        A.assembleLhsRhsBc((eps11*(nv(G).tr()*firstCoeff)) * eps11.tr(),
                              zero*eps11,
                              bc.dirichletSides());
        gsInfo << "BndMatrix: \n" << A.matrix().toDense() << "\n";
        //gsInfo << "solVector: \n" << solVector << "\n";
        gsInfo << "Finished t= " << dt*it << "\n";
        //gsInfo << "dof: " << ndof << "\n";
    }

    if (plot){
    gsInfo << "Plotting in Paraview...\n";


    ev.options().setSwitch("plot.elements", true);
    ev.writeParaview(u_sol, G, "solutionU");
    ev.writeParaview(v_sol, G, "solutionV");
    ev.writeParaview(eps11_sol, G, "solutionEps11");
    ev.writeParaview(eps22_sol, G, "solutionEps22");
    ev.writeParaview(eps12_sol, G, "solutionEps12");
    //ev.writeParaview( u_ex    , G, "solution_ex");
    //ev.writeParaview( u, G, "aa");
    //gsFileManager::open("solution.pvd");
    }
    return EXIT_SUCCESS;

}// end main
