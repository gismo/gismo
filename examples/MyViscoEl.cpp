/** @file MyViscoEl.cpp

    @Replicating viscoelasticity results

    Author(s): A. Barion
*/

//! [Include namespace]
#include <gismo.h>
#include <chrono>

using namespace gismo;
//! [Include namespace]
void writeToCSVfile(std::string name, gsMatrix<> matrix)
{
    std::ofstream file(name.c_str());
    for(int  i = 0; i < matrix.rows(); i++){
        for(int j = 0; j < matrix.cols(); j++){
           std::string str = std::to_string(matrix(i,j));
           if(j+1 == matrix.cols()){
               file<<str;
           }else{
               file<<str<<',';
           }
        }
        file<<'\n';
    }
  }

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
    real_t zero = 0.0;

    // Define Geometry, must be a gsMultiPatch object
    gsMultiPatch<> patch;

    // Create one patch as the unit square [0,12]^2
    patch = gsNurbsCreator<>::BSplineSquareGrid(1, 1, 12);
    gsInfo << "The domain is a " << patch << "\n";

    //Source function for velocities
    gsFunctionExpr<> f1("if ((abs(x-6.0)<=1.5) and (abs(y-6.0)<=1.5), -sgn(x-6.0), 0.0)", 2);
    gsFunctionExpr<> f2("0.0", 2);
    
    

    gsInfo << "Source function for first velocity component: " << f1 << "\n";
    gsInfo << "Source function for second velocity component: " << f2 << "\n";


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
    pt<<0,0.25;
    // pt.setConstant(0.25);
    
    variable ff1 = A.getCoeff(f1, G);
    variable ff2 = A.getCoeff(f2, G);

    // Solution vector and solution variable

    A.initSystem();

    gsSparseSolver<>::BiCGSTABILUT solver;

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
    gsMatrix<> RHS;
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

    for(index_t it = 1; it <= numSteps; ++it){

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

        //gsInfo << "u_aux:\n" << u_aux << "\n";
        //gsInfo << "v_aux:\n" << v_aux << "\n";

        displ1 = displ1 + dt * u_aux;
        displ2 = displ2 + dt * v_aux;
        displVec.segment(0,displ_dof) = displ1;
        displVec.segment(displ_dof, displ_dof) = displ2;
        
        //gsInfo << "displDof:" << displ_dof << "\n";
        A.initSystem();
        A.assemble(eps11*(igrad(eps11,G)*firstCoeff).tr()*meas(G));
        
        A.assemble(eps22*(igrad(eps22,G)*secondCoeff).tr()*meas(G));

        A.assemble((half*eps12)*(igrad(eps11,G)*secondCoeff).tr()*meas(G));
        A.assemble((half*eps12)*(igrad(eps22,G)*firstCoeff).tr()*meas(G));
        
        
        displContrib = A.matrix() * displVec;

        //gsInfo << "mat:\n" << A.matrix() << "\n";
        
        A.initSystem();
        damp = 4.2*(1.0-exp(-4.0*(dt*it-0.1)/(20.0-0.1)));
        //gsInfo << "damp: " << damp << "\n";

        //Mass matrices for strains and velocity components
        A.assemble(eps11 * eps11.tr() * meas(G));
        A.assemble(eps22 * eps22.tr() * meas(G));
        A.assemble(eps12 * eps12.tr() * meas(G));
        A.assemble(rho * u * u.tr() * meas(G), damp * u * ff1 * meas(G));
        A.assemble(rho * v * v.tr() * meas(G), damp * v * ff2 * meas(G));

        //Strain contributions to vel comps
        //eps11-->u
        A.assembleLhsRhsBc((-dt)*G_tilde*(1-nu)/(1-2*nu)*(u*(nv(G).tr()*firstCoeff)) * eps11.tr(),
                              zero*u,
                              bc.reducedContainer(bc.dirichletSides(), 3));
        A.assemble((dt)*G_tilde*(1-nu)/(1-2*nu)*((igrad(u,G)*firstCoeff))*eps11.tr()*meas(G));
        
        //eps22-->u
        A.assembleLhsRhsBc((-dt)*G_tilde*(nu)/(1-2*nu)*(u*(nv(G).tr()*firstCoeff)) * eps22.tr(),
                               zero*u,
                               bc.reducedContainer(bc.dirichletSides(), 3));
        A.assemble((dt)*G_tilde*(nu)/(1-2*nu)*((igrad(u,G)*firstCoeff))*eps22.tr()*meas(G));
        
        //eps12-->u
        A.assembleLhsRhsBc((-dt)*G_tilde*(u*(nv(G).tr()*secondCoeff)) * eps12.tr(),
                               zero*u,
                               bc.reducedContainer(bc.dirichletSides(), 3));
        A.assemble((dt)*G_tilde*((igrad(u,G)*secondCoeff))*eps12.tr()*meas(G));

        //eps11-->v
        A.assembleLhsRhsBc((-dt)*G_tilde*(nu)/(1-2*nu)*(v*(nv(G).tr()*secondCoeff)) * eps11.tr(),
                               zero*v,
                               bc.reducedContainer(bc.dirichletSides(), 4));
        A.assemble((dt)*G_tilde*(nu)/(1-2*nu)*((igrad(v,G)*secondCoeff))*eps11.tr()*meas(G));

        //eps22-->v
        A.assembleLhsRhsBc((-dt)*G_tilde*(1-nu)/(1-2*nu)*(v*(nv(G).tr()*secondCoeff)) * eps22.tr(),
                               zero*v,
                               bc.reducedContainer(bc.dirichletSides(), 4));
        A.assemble((dt)*G_tilde*(1-nu)/(1-2*nu)*((igrad(v,G)*secondCoeff))*eps22.tr()*meas(G));

        //eps12-->v
        A.assembleLhsRhsBc((-dt)*G_tilde*(v*(nv(G).tr()*firstCoeff)) * eps12.tr(),
                               zero*v,
                               bc.reducedContainer(bc.dirichletSides(), 4));
        A.assemble((dt)*G_tilde*((igrad(v,G)*firstCoeff))*eps12.tr()*meas(G));

        // gsDebug<<"(firstCoeff) = \n"<<evEps.eval(firstCoeff_Eps,pt)<<"\n";
        //gsDebug<<"igrad(v,G) = \n"<<ev.eval(igrad(v,G),pt)<<"\n";
        //gsDebug<<"(nv(G).tr()*firstCoeff/nv(G).norm()) = \n"<<ev.eval(nv(G).tr()*firstCoeff/nv(G).norm(),pt)<<"\n";


        //Velocity contributions
       
        //u,u
        A.assemble(dt*(mu1+mu2)*Sxx_u);
        A.assemble(dt*mu1/2*Syy_u);
        A.assembleLhsRhsBc((-dt)*(mu1+mu2)*(u*(nv(G).tr()*firstCoeff))*(igrad(u,G)*firstCoeff).tr(),
                               zero*u,
                               bc.reducedContainer(bc.dirichletSides(), 3));
        A.assembleLhsRhsBc((-dt)*mu1/2*(u*(nv(G).tr()*secondCoeff))*(igrad(u,G)*secondCoeff).tr(),
                               zero*u,
                               bc.reducedContainer(bc.dirichletSides(), 3));

        //u,v
        A.assemble(dt*mu2*(igrad(u,G)*firstCoeff)*(igrad(v,G)*secondCoeff).tr()*meas(G));
        A.assemble(dt*mu1/2*(igrad(u,G)*secondCoeff)*(igrad(v,G)*firstCoeff).tr()*meas(G));
        A.assembleLhsRhsBc((-dt)*mu2*(u*(nv(G).tr()*firstCoeff))*(igrad(v,G)*secondCoeff).tr(),
                                zero*u,
                                bc.reducedContainer(bc.dirichletSides(), 3));
        A.assembleLhsRhsBc((-dt)*mu1/2*(u*(nv(G).tr()*secondCoeff))*(igrad(v,G)*firstCoeff).tr(),
                                zero*u,
                                bc.reducedContainer(bc.dirichletSides(), 3));


        //v,v
        A.assemble(dt*(mu1+mu2)*Syy_v);
        A.assemble(dt*mu1/2*Sxx_v);
        A.assembleLhsRhsBc((-dt)*(mu1+mu2)*(v*(nv(G).tr()*secondCoeff))*(igrad(v,G)*secondCoeff).tr(),
                               zero*v,
                               bc.reducedContainer(bc.dirichletSides(), 4));
        A.assembleLhsRhsBc((-dt)*mu1/2*(v*(nv(G).tr()*firstCoeff))*(igrad(v,G)*firstCoeff).tr(),
                               zero*v,
                               bc.reducedContainer(bc.dirichletSides(), 4));


        //v,u
        A.assemble(dt*mu1/2*(igrad(v,G)*firstCoeff)*(igrad(u,G)*secondCoeff).tr()*meas(G));
        A.assemble(dt*mu2*(igrad(v,G)*secondCoeff)*(igrad(u,G)*firstCoeff).tr()*meas(G));
        A.assembleLhsRhsBc((-dt)*mu1/2*(v*(nv(G).tr()*firstCoeff))*(igrad(u,G)*secondCoeff).tr(),
                               zero*v,
                               bc.reducedContainer(bc.dirichletSides(), 4));
        A.assembleLhsRhsBc((-dt)*mu2*(v*(nv(G).tr()*secondCoeff))*(igrad(u,G)*firstCoeff).tr(),
                               zero*v,
                               bc.reducedContainer(bc.dirichletSides(), 4));


        auto finish = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = finish - start;
        gsInfo << "Elapsed time: " << elapsed.count() << " s\n";
        //gsInfo << "A.matrix():\n" << A.matrix() << "\n";

        solveMat = A.matrix();
        solveMat.makeCompressed();
        solver.compute(solveMat);
        
        RHS = A.rhs();
        //gsInfo << "rhs:\n " << RHS << "\n";

        
        F = prevTimestep + displContrib + dt*RHS;
        //gsInfo << "rhs:\n " << RHS << "\n";
        //writeToCSVfile("matrix.csv",A.matrix().toDense());
        
        solVector = solver.solve(F);
        //solVector = (solveMat.toDense()).colPivHouseholderQr().solve(F); //dense linear algebra
        
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
