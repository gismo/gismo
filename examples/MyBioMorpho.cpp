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
    real_t E = 10.0;
    real_t mu1 = 10.0;
    real_t mu2 = 10.0;
    real_t rho_t = 1.02;
    real_t G_tilde = E*sqrt(rho_t)/(1+nu); // TO BE CHANGED TO rho
    real_t zeta = 0.3;
    real_t t_ap = 5.0;
    real_t tau_c = 20.0;
    real_t D_F = 1.0;
    real_t chi_F = 1.0;
    real_t q = 1-0.42;
    real_t r_F = 1.0;
    real_t r_F_max = 1.0;
    real_t a_c_I = 1.0;
    real_t kappa_F = 1.0;


    // Define Geometry, must be a gsMultiPatch object
    gsMultiPatch<> patch;

    // Create one patch as the unit square [0,12]^2
    patch = gsNurbsCreator<>::BSplineSquareGrid(1, 1, 12);
    gsInfo << "The domain is a " << patch << "\n";

    //Source function for velocities
    gsFunctionExpr<> f1("if ((abs(x-6)<=1.5) and (abs(y-6)<=1.5), -sgn(x-6), 0)", 2);
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

    gsBoundaryConditions<> bc;
    gsBoundaryConditions<> bc_bio;
    // bc.addCondition(boundary::north, condition_type::dirichlet, &fun, 0);
    bc.addCondition(boundary::north, condition_type::dirichlet, 0, 3);
    bc.addCondition(boundary::east, condition_type::dirichlet, 0, 3);
    bc.addCondition(boundary::south, condition_type::dirichlet, 0, 3);
    bc.addCondition(boundary::west, condition_type::dirichlet, 0, 3);
    bc.addCondition(boundary::north, condition_type::dirichlet, 0, 4);
    bc.addCondition(boundary::east, condition_type::dirichlet, 0, 4);
    bc.addCondition(boundary::south, condition_type::dirichlet, 0, 4);
    bc.addCondition(boundary::west, condition_type::dirichlet, 0, 4);

    gsFunctionExpr<> N_bar("10000",2);
    gsFunctionExpr<> M_bar("0",2);
    gsFunctionExpr<> c_bar("0",2);
    gsFunctionExpr<> rho_bar("0.1",2);


    bc_bio.addCondition(boundary::north, condition_type::dirichlet, &N_bar, 0);
    bc_bio.addCondition(boundary::east, condition_type::dirichlet, &N_bar, 0);
    bc_bio.addCondition(boundary::south, condition_type::dirichlet, &N_bar, 0);
    bc_bio.addCondition(boundary::west, condition_type::dirichlet, &N_bar, 0);
    bc_bio.addCondition(boundary::north, condition_type::dirichlet, &M_bar, 1);
    bc_bio.addCondition(boundary::east, condition_type::dirichlet, &M_bar, 1);
    bc_bio.addCondition(boundary::south, condition_type::dirichlet, &M_bar, 1);
    bc_bio.addCondition(boundary::west, condition_type::dirichlet, &M_bar, 1);
    bc_bio.addCondition(boundary::north, condition_type::dirichlet, &c_bar, 2);
    bc_bio.addCondition(boundary::east, condition_type::dirichlet, &c_bar, 2);
    bc_bio.addCondition(boundary::south, condition_type::dirichlet, &c_bar, 2);
    bc_bio.addCondition(boundary::west, condition_type::dirichlet, &c_bar, 2);
    bc_bio.addCondition(boundary::north, condition_type::dirichlet, &rho_bar, 3);
    bc_bio.addCondition(boundary::east, condition_type::dirichlet, &rho_bar, 3);
    bc_bio.addCondition(boundary::south, condition_type::dirichlet, &rho_bar, 3);
    bc_bio.addCondition(boundary::west, condition_type::dirichlet, &rho_bar, 3);
    
    gsInfo << "Boundary conditions:\n" << bc << "\n";
    gsInfo << "Boundary conditions:\n" << bc_bio << "\n";

    //! [Problem setup]
    gsExprAssembler<> A(5, 5);                              // velocity+strain assembler
    gsExprAssembler<> B(5, 5);                              // cell constituents assembler

    gsOptionList opt = gsAssembler<>::defaultOptions();
    //opt.setInt("DirichletValues", dirichlet::l2projection);
    opt.setInt("DirichletStrategy", dirichlet::nitsche);
    opt.setReal("quA", 1.0);
    opt.setInt("quB", 1);
    opt.setReal("bdA", 2.0);
    opt.setInt("bdB", 1);
    //gsInfo << "Assembler " << opt;
    A.setOptions(opt);
    B.setOptions(opt);
    

    //gsInfo << "Active options:\n" << A.options() << "\n";
    typedef gsExprAssembler<>::geometryMap geometryMap;
    typedef gsExprAssembler<>::variable    variable;
    typedef gsExprAssembler<>::space       space;
    typedef gsExprAssembler<>::solution    solution;

    //typedef typename gsExprHelper<real_t>::nullExpr    nullExpr;


    // Elements used for numerical integration
    A.setIntegrationElements(dbasis);
    B.setIntegrationElements(dbasis);
    gsExprEvaluator<> ev(A);
    gsExprEvaluator<> evB(B);


    // Set the geometry map
    geometryMap G = A.getMap(patch);
    geometryMap G_bio = B.getMap(patch);

    // Set the discretization space
    space eps11 = A.getSpace(dbasis,1,0);
    space eps22 = A.getSpace(dbasis,1,1);
    space eps12 = A.getSpace(dbasis,1,2);
    space u = A.getSpace(dbasis,1,3);
    space v = A.getSpace(dbasis,1,4);

    space N = B.getSpace(dbasis,1,0);
    space M = B.getSpace(dbasis,1,1);
    space c = B.getSpace(dbasis,1,2);
    space rho = B.getSpace(dbasis,1,3);


    u.setInterfaceCont(0); // todo: 1 (smooth basis)
    u.addBc(bc.get("Dirichlet",3)); // (!) must be called only once
    v.setInterfaceCont(0); // todo: 1 (smooth basis)
    v.addBc(bc.get("Dirichlet",4)); // (!) must be called only once

    N.setInterfaceCont(0); // todo: 1 (smooth basis)
    N.addBc(bc.get("Dirichlet",0)); // (!) must be called only once
    M.setInterfaceCont(0); // todo: 1 (smooth basis)
    M.addBc(bc.get("Dirichlet",1)); // (!) must be called only once
    c.setInterfaceCont(0); // todo: 1 (smooth basis)
    c.addBc(bc.get("Dirichlet",2)); // (!) must be called only once
    rho.setInterfaceCont(0); // todo: 1 (smooth basis)
    rho.addBc(bc.get("Dirichlet",3)); // (!) must be called only once

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
    B.initSystem();

    gsSparseSolver<>::BiCGSTABILUT solver;

    gsInfo << "Number of degrees of freedom: " << A.numDofs() << "\n";

    index_t ndof = A.numDofs();
    index_t ndof_bio = B.numDofs();

    gsMatrix<> solVector;
    gsMatrix<> solVector_bio;
    
    solution u_sol = A.getSolution(u, solVector);
    solution v_sol = A.getSolution(v, solVector);
    solution eps11_sol = A.getSolution(eps11, solVector);
    solution eps22_sol = A.getSolution(eps22, solVector);
    solution eps12_sol = A.getSolution(eps12, solVector);

    solution N_sol = A.getSolution(N, solVector_bio);
    solution M_sol = A.getSolution(M, solVector_bio);
    solution c_sol = A.getSolution(c, solVector_bio);
    solution rho_sol = A.getSolution(rho, solVector_bio);
    

    //Initial conditions
    solVector.setZero(ndof,1);
    solVector_bio.setZero(ndof_bio,1); //TO DO!!!!!!!!!!!!!!!!!!!!

    gsMatrix<> aux;
    gsVector<> displ1;
    gsVector<> displ2;
    gsMatrix<> prevTimestep;

    eps11.getCoeffs(solVector,aux);
    displ1.setZero(aux.size());
    displ2 = displ1;
    gsSparseMatrix<> oldMass;
    gsMatrix<> u_aux;
    gsMatrix<> v_aux; 
    u.getCoeffs(solVector,u_aux);
    v.getCoeffs(solVector,v_aux);
    real_t half = 0.5;
    real_t zero = 0;
    real_t one = 1;
    real_t damp;
    real_t tol = pow(10,-4);

    gsMatrix<> F;
    gsSparseMatrix<> solveMat;

    auto Sxx_u = (igrad(u,G)*firstCoeff)*(igrad(u,G)*firstCoeff).tr()*meas(G);
    auto Syy_u = (igrad(u,G)*secondCoeff)*(igrad(u,G)*secondCoeff).tr()*meas(G);
    auto Sxx_v = (igrad(v,G)*firstCoeff)*(igrad(v,G)*firstCoeff).tr()*meas(G);
    auto Syy_v = (igrad(v,G)*secondCoeff)*(igrad(v,G)*secondCoeff).tr()*meas(G);

    gsSparseMatrix<> M_Linear;
    gsMatrix<> picardVector;
    gsMatrix<> aux_bio;
    aux_bio = solVector_bio;
    picardVector = solVector;
    solution N_old = A.getSolution(N, aux_bio);
    solution M_old = A.getSolution(M, aux_bio);
    solution c_old = A.getSolution(c, aux_bio);
    solution rho_old = A.getSolution(rho, aux_bio);
    solution u_old = A.getSolution(u, picardVector);
    solution v_old = A.getSolution(v, picardVector);
    real_t diff;

    //plotting
    gsParaviewCollection collectionU("solutionU");
    gsParaviewCollection collectionV("solutionV");
    gsParaviewCollection collectionE11("solutionE11");
    gsParaviewCollection collectionE22("solutionE22");
    gsParaviewCollection collectionE12("solutionE12");
    ev.options().setSwitch("plot.elements", true);
    //! plotting

    for(index_t it = 1; it <= 1; ++it){

        B.initSystem();
        B.assemble(N * N.tr() * meas(G_bio));
        B.assemble(M * M.tr() * meas(G_bio));
        B.assemble(c * c.tr() * meas(G_bio));
        B.assemble(rho * rho.tr() * meas(G_bio));
        oldMass = B.matrix();
        prevTimestep = oldMass * solVector_bio;

        B.initSystem();
        //N fluxes: diffusive and convective
        B.assemble(-D_F*(N_old+M_old)*igrad(N,G_bio)*igrad(N,G_bio)*meas(G_bio)); //diffusive
        B.assemble(chi_F*(igrad(N,G_bio)*grad(c_old).tr())*N.tr()*meas(G_bio)); //convective

        //N bndFluxes
        B.assembleLhsRhsBc(-D_F*(N_old+M_old)*N*(igrad(N,G_bio)*nv(G_bio)).tr(),
                              zero*N,
                              bc_bio.reducedContainer(bc_bio.dirichletSides(), 0)); //diffusive
        B.assembleLhsRhsBc(chi_F*(N*(grad(c_old)*nv(G_bio)))*N.tr(),
                              zero*N,
                              bc_bio.reducedContainer(bc_bio.dirichletSides(), 0)); //diffusive

        //N Forcing
        B.assemble(r_F*(one+(r_F_max*c_old)/(a_c_I+c_old.val()))*(-kappa_F*pow((N_old.val()),1.0+q))*N*N.tr()*meas(G_bio),
                    r_F*(one+(r_F_max*c_old)/(a_c_I+c_old.val()))*pow((N_old.val()),1.0+q)*N.tr()*meas(G_bio)+
                    r_F*(one+(r_F_max*c_old)/(a_c_I+c_old.val()))*(-kappa_F*M_old*pow((N_old.val()),1.0+q))*N.tr()*meas(G_bio));



        //auto start = std::chrono::high_resolution_clock::now();
        A.initSystem();
        A.assemble(eps11 * eps11.tr() * meas(G));
        A.assemble(eps22 * eps22.tr() * meas(G));
        A.assemble(eps12 * eps12.tr() * meas(G));
        A.assemble(rho_t * u * u.tr() * meas(G));
        A.assemble(rho_t * v * v.tr() * meas(G));
        oldMass = A.matrix();
        prevTimestep = oldMass * solVector;

        u.getCoeffs(solVector,u_aux);
        v.getCoeffs(solVector,v_aux);

        patch.patch(0).coefs().col(0) += dt*u_aux;
        patch.patch(0).coefs().col(1) += dt*v_aux;

        displ1 = displ1 + dt * u_aux;
        displ2 = displ2 + dt * v_aux;
        
        A.initSystem();
        A.assemble(eps11*(igrad(eps11,G)*firstCoeff).tr()*meas(G));
        
        A.assemble(eps22*(igrad(eps22,G)*secondCoeff).tr()*meas(G));

        A.assemble((half*eps12)*(igrad(eps11,G)*secondCoeff).tr()*meas(G) +
                    (half*eps12)*(igrad(eps22,G)*firstCoeff).tr()*meas(G));
        
        A.initSystem();
        
        damp=0;
        if(dt*it < t_ap){
            damp=1;
        }
        

        //Mass matrices for strains and velocity components
        A.assemble((1 + zeta * dt) * eps11 * eps11.tr() * meas(G));
        A.assemble((1 + zeta * dt) * eps22 * eps22.tr() * meas(G));
        A.assemble((1 + zeta * dt) * eps12 * eps12.tr() * meas(G));
        A.assemble(rho_t * u * u.tr() * meas(G), damp * tau_c * u * ff1 * meas(G));
        A.assemble(rho_t * v * v.tr() * meas(G), damp * tau_c * v * ff2 * meas(G));

        //Vel comps contributions to strain
        //u-->eps11
        A.assemble((-dt)*eps11*(igrad(u,G)*firstCoeff).tr()*meas(G));

        //v-->eps22
        A.assemble((-dt)*eps22*(igrad(v,G)*secondCoeff).tr()*meas(G));

        //u,v-->eps12
        A.assemble((-dt/2)*eps12*(igrad(u,G)*secondCoeff).tr()*meas(G));
        A.assemble((-dt/2)*eps12*(igrad(v,G)*firstCoeff).tr()*meas(G));


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

        //auto finish = std::chrono::high_resolution_clock::now();
        //std::chrono::duration<double> elapsed = finish - start;
        //gsInfo << "Elapsed time: " << elapsed.count() << " s\n";
        
        M_Linear = A.matrix();
        F = prevTimestep + dt*A.rhs();

        // ! Linear assembly
             

        // Picard loop
        for(index_t pic=0; pic < 5; ++pic){

            gsInfo << "Picard iteration: " << pic+1 << "\n";

            A.initSystem();
            //eps11
            A.assemble((-dt) * eps11 * eps11.tr() * ( grad(v_old)*secondCoeff ).val() * meas(G) );
            A.assemble((dt) * eps11 * eps22.tr() * ( grad(u_old)*firstCoeff ).val() * meas(G) );
            A.assemble((-dt) * eps11 * eps12.tr() * ( grad(u_old)*secondCoeff - grad(v_old) * firstCoeff).val() * meas(G) );
            
            //eps22
            A.assemble((-dt) * eps22 * eps22.tr() * ( grad(u_old)*firstCoeff ).val() * meas(G) );
            A.assemble((dt) * eps22 * eps11.tr() * ( grad(v_old)*secondCoeff ).val() * meas(G) );
            A.assemble((-dt) * eps22 * eps12.tr() * ( grad(v_old)*firstCoeff - grad(u_old) * secondCoeff).val() * meas(G) );

            //eps12
            A.assemble((-dt) * eps12 * eps12.tr() * ( grad(u_old)*firstCoeff + grad(v_old) * secondCoeff ).val() * meas(G) );
            A.assemble((dt) * eps12 * eps11.tr() * ( grad(u_old)*secondCoeff ).val() * meas(G) );
            A.assemble((dt) * eps12 * eps22.tr() * ( grad(v_old)*firstCoeff ).val() * meas(G) );

            solveMat = M_Linear + A.matrix();
            //gsInfo << "Mat to solve: " << solveMat << "\n"; //check if sum remains sparse
            solver.compute(solveMat);

            picardVector = solver.solve(F);

            diff = ((solVector-picardVector).cwiseAbs()).maxCoeff();

            solVector = picardVector; //update solution to last picard iteration

            if(diff<tol){
                gsInfo << "Stopped at picard iteration: " << pic+1 << "\n";
                gsInfo << "Diff = " << diff << "\n";
                break;
            }
        }

        // A.initSystem();
        // gsDebug << "nv(G)/nv(G).norm()" << ev.eval(nv(G).norm(), pt) << "\n";
        // A.assembleLhsRhsBc((eps11*((nv(G)/nv(G).norm()).tr()*firstCoeff)) * eps11.tr(),
        //                       zero*eps11,
        //                       bc.dirichletSides());
        // gsInfo << "BndMatrix: \n" << A.matrix().toDense() << "\n";
        //gsInfo << "solVector: \n" << solVector << "\n";
        gsInfo << "Finished t= " << dt*it << "\n";
        //gsInfo << "dof: " << ndof << "\n";


        if(plot){
            ev.writeParaview(u_sol, G, "solutionU_"+ util::to_string(it));
            ev.writeParaview(v_sol, G, "solutionV_" + util::to_string(it));
            ev.writeParaview(eps11_sol, G, "solutionEps11_" + util::to_string(it));
            ev.writeParaview(eps22_sol, G, "solutionEps22_" + util::to_string(it));
            ev.writeParaview(eps12_sol, G, "solutionEps12_" + util::to_string(it));

            collectionU.addTimestep("solutionU_" + util::to_string(it),it,"0.vts");
            collectionU.addTimestep("solutionU_" + util::to_string(it),it,"0_mesh.vtp");

            collectionV.addTimestep("solutionV_" + util::to_string(it),it,"0.vts");
            collectionV.addTimestep("solutionV_" + util::to_string(it),it,"0_mesh.vtp");

            collectionE11.addTimestep("solutionEps11_" + util::to_string(it),it,"0.vts");
            collectionE11.addTimestep("solutionEps11_" + util::to_string(it),it,"0_mesh.vtp");

            collectionE22.addTimestep("solutionEps22_" + util::to_string(it),it,"0.vts");
            collectionE22.addTimestep("solutionEps22_" + util::to_string(it),it,"0_mesh.vtp");

            collectionE12.addTimestep("solutionEps12_" + util::to_string(it),it,"0.vts");
            collectionE12.addTimestep("solutionEps12_" + util::to_string(it),it,"0_mesh.vtp");
        }
    }

    if(plot){
        collectionU.save();
        collectionV.save();
        collectionE11.save();
        collectionE22.save();
        collectionE12.save();
    }

    // gsInfo << "Plotting in Paraview...\n";


    // ev.options().setSwitch("plot.elements", true);
    // ev.writeParaview(u_sol, G, "solutionU");
    // ev.writeParaview(v_sol, G, "solutionV");
    // ev.writeParaview(eps11_sol, G, "solutionEps11");
    // ev.writeParaview(eps22_sol, G, "solutionEps22");
    // ev.writeParaview(eps12_sol, G, "solutionEps12");

    return EXIT_SUCCESS;

}// end main
