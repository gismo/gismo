/** @file rh_elasticity_example.cpp

    @brief Tutorial on how to use expression assembler to solve a linear elasticity equation

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris & M. BAHARI
*/

//! [Include namespace]
#include <gismo.h>
#include <fstream> // For file operations
#include <gsElasticity/gsElasticityAssembler.h>
#include <gsElasticity/gsWriteParaviewMultiPhysics.h>

using namespace gismo;
//! [Include namespace]

void ProjectionNormalCPoints(gsMultiPatch<>& Psi, gsMultiPatch<> mpLeft, bool Left){
    int boxMaxNumber = 1;
    // Projection normal of control points (exact geometry)
    if (Left){
    //gsInfo << "ProjectionNormalCPoints for left geometry\n";
    for (int boxNumber = 0; boxNumber < boxMaxNumber; ++boxNumber)
    {
        // test if the boundary interface is not an inner interface between patches
        auto lVal = 0.; //int(1.1*Psi.patch(boxNumber).coef( Psi.patch(boxNumber).basis().boundary(1).at(0) ).array()[0]);
        auto hVal = 1.; //int(1.1*Psi.patch(boxNumber).coef( Psi.patch(boxNumber).basis().boundary(2).at(0) ).array()[0]);
        for (int i_x =0; i_x < Psi.patch(boxNumber).basis().boundary(1).size(); ++i_x) // x=0 control points be like (0,:) in this case
        {
            Psi.patch(boxNumber).coef( Psi.patch(boxNumber).basis().boundary(1).at(i_x) ).array()[0] = lVal;
        }

        for (int i_x =0; i_x < Psi.patch(boxNumber).basis().boundary(2).size(); ++i_x)// x=1 control points be like (1,:) in this case
        {
        Psi.patch(boxNumber).coef( Psi.patch(boxNumber).basis().boundary(2).at(i_x) ).array()[0] = hVal;
        }

        lVal = 0.; //int(1.1*Psi.patch(boxNumber).coef( Psi.patch(boxNumber).basis().boundary(3).at(0) ).array()[1]);
        hVal = 1.; //int(1.1*Psi.patch(boxNumber).coef( Psi.patch(boxNumber).basis().boundary(4).at(0) ).array()[1]);
        for (int i_x =0; i_x < Psi.patch(boxNumber).basis().boundary(3).size(); ++i_x) // y=0 control points be like (:,0) in this case
        {
        Psi.patch(boxNumber).coef( Psi.patch(boxNumber).basis().boundary(3).at(i_x) ).array()[1] = lVal;
        }
        for (int i_x =0; i_x < Psi.patch(boxNumber).basis().boundary(4).size(); ++i_x)// y=1 control points be like (:,1) in this case
        {
        Psi.patch(boxNumber).coef( Psi.patch(boxNumber).basis().boundary(4).at(i_x) ).array()[1] = hVal;
        }
    }
}
    else{
    gsInfo << "ProjectionNormalCPoints for right geometry\n";
    for (int boxNumber = 0; boxNumber < boxMaxNumber; ++boxNumber)
    {
        // test if the boundary interface is not an inner interface between patches
        for (int i_x =0; i_x < Psi.patch(boxNumber).basis().boundary(1).size(); ++i_x) // x=0 control points be like (0,:) in this case
        {
            //auto lVal = mpLeft.patch(boxNumber).coef( mpLeft.patch(boxNumber).basis().boundary(1).at(i_x) ).array()[0];
            auto hVal = mpLeft.patch(boxNumber).coef( mpLeft.patch(boxNumber).basis().boundary(1).at(i_x) ).array()[1];
            //Psi.patch(boxNumber).coef( Psi.patch(boxNumber).basis().boundary(1).at(i_x) ).array()[0] = lVal;
            Psi.patch(boxNumber).coef( Psi.patch(boxNumber).basis().boundary(1).at(i_x) ).array()[1] = hVal;
        }

        for (int i_x =0; i_x < Psi.patch(boxNumber).basis().boundary(2).size(); ++i_x)// x=1 control points be like (1,:) in this case
        {
            auto lVal = mpLeft.patch(boxNumber).coef( mpLeft.patch(boxNumber).basis().boundary(2).at(i_x) ).array()[0];
            //auto hVal = mpLeft.patch(boxNumber).coef( mpLeft.patch(boxNumber).basis().boundary(2).at(i_x) ).array()[1];
            Psi.patch(boxNumber).coef( Psi.patch(boxNumber).basis().boundary(2).at(i_x) ).array()[0] = lVal;
            //Psi.patch(boxNumber).coef( Psi.patch(boxNumber).basis().boundary(2).at(i_x) ).array()[1] = hVal;
        }
        for (int i_x =0; i_x < Psi.patch(boxNumber).basis().boundary(3).size(); ++i_x) // y=0 control points be like (:,0) in this case
        {
            auto lVal = mpLeft.patch(boxNumber).coef( mpLeft.patch(boxNumber).basis().boundary(3).at(i_x) ).array()[0];
            auto hVal = mpLeft.patch(boxNumber).coef( mpLeft.patch(boxNumber).basis().boundary(3).at(i_x) ).array()[1];
            Psi.patch(boxNumber).coef( Psi.patch(boxNumber).basis().boundary(3).at(i_x) ).array()[0] = lVal;
            Psi.patch(boxNumber).coef( Psi.patch(boxNumber).basis().boundary(3).at(i_x) ).array()[1] = hVal;
        }
        for (int i_x =0; i_x < Psi.patch(boxNumber).basis().boundary(4).size(); ++i_x)// y=1 control points be like (:,1) in this case
        {
            auto lVal = mpLeft.patch(boxNumber).coef( mpLeft.patch(boxNumber).basis().boundary(4).at(i_x) ).array()[0];
            auto hVal = mpLeft.patch(boxNumber).coef( mpLeft.patch(boxNumber).basis().boundary(4).at(i_x) ).array()[1];
            Psi.patch(boxNumber).coef( Psi.patch(boxNumber).basis().boundary(4).at(i_x) ).array()[0] = lVal;
            Psi.patch(boxNumber).coef( Psi.patch(boxNumber).basis().boundary(4).at(i_x) ).array()[1] = hVal;
        }
    }
    }
};

int main(int argc, char *argv[])
{
    //! [Parse command line]
    bool plot             = true;
    index_t numRefine     = 3;
    index_t numLRefine    = 0;
    index_t numElevate    = 0;
    index_t maxIter       = 50;
    index_t NumArMarEl    = 0; // Number of ring of cells around marked elements
    double eps            = 1e-7; // pinalization coefficient
    double tolPicard      = 1e-8;
    double IntensityMAE   = 6.;
    //bool export_b64     = false;
    bool errorsave        = false;
    real_t adaptRefParam  = 0.;     // ... adapt parameter.
    index_t FactRefPar    = 0;    // ... adapt parameter : adaptRefParam += FactRefPar in each iter
    // Specify the file path
    std::string fn("pde/infinit_plate.xml");

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
    cmd.addReal( "a", "adaptRefParam", "parameter for local h-refinement loops",  adaptRefParam );
    cmd.addInt( "p", "FactRefPar", "augement adaptRefParam with such quantity in local h-refinement loops",  FactRefPar );
    cmd.addSwitch("plot", "Create a ParaView visualization file with the solution", plot);
    cmd.addSwitch("errorsave", "Create a file in ... and save errors", errorsave);
    cmd.addReal( "f", "IntensityMAE", "Intensity of density function",  IntensityMAE);
    cmd.addInt( "c", "NumArMarEl", "augement NumArMarEl with such quantity in local h-refinement loops",  NumArMarEl );

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    // Load the file
    gsFileData<> fd(fn);
    gsInfo << "Loaded file " << fd.lastPath() << "\n";
    // Create a gsMultipatch and add the loaded geometry
    gsMultiPatch<> mpLeft;
    fd.getId(1,mpLeft);
    // Elevate and p-refine the basis to order p + numElevate
    // where p is the highest degree in the bases
    mpLeft.degreeElevate(numElevate);
    mpLeft.computeTopology();

    // .... one single patch
    gsMultiPatch<> mp = gsNurbsCreator<>::BSplineSquareGrid(1,1,1, 0.0, 0.0);
    //Get all interfaces and boundaries:
    mp.degreeElevate(numElevate);
    mp.computeTopology();
    //mp.addAutoBoundaries();

    // Manufactured identity mapping
    gsFunctionExpr<> sN("x","y",2);
    // Right-hand side function : Analytical density function (det(H(u))=f= sigma/rho)
    // Load the file
    gsFunctionExpr<> f;
    fd.getId(2003, f);
    gsFunctionExpr<> rhs;
    fd.getId(2001, rhs);
    gsInfo<<"Density function "<< f << "\n";

    gsBoundaryConditions<> bc_mae;
    bc_mae.setGeoMap(mp);
    // For simplicity, set Neumann boundary conditions
   for ( gsMultiPatch<>::const_biterator
            bit = mp.bBegin(); bit != mp.bEnd(); ++bit)
   {
       bc_mae.addCondition( *bit, condition_type::neumann, &sN );
   }
    gsInfo<<"Boundary conditions:\n"<< bc_mae <<"\n";

    //! [Refinement]
    gsMultiBasis<> dbasis(mpLeft, true);//true: poly-splines (not NURBS)
    
    gsInfo << "Patches: "<< mp.nPatches() <<", degree: "<< dbasis.minCwiseDegree() <<"\n";
#ifdef _OPENMP
    gsInfo<< "Available threads: "<< omp_get_max_threads() <<"\n";
#endif
    //! [Refinement]
    
    //! [Problem setup]
    gsExprAssembler<> A(1,1);
    //A.setOptions(Aopt);
    gsInfo<<"Active options:\n"<< A.options() <<"\n";

    typedef gsExprAssembler<>::geometryMap geometryMap;
    typedef gsExprAssembler<>::variable    variable;
    typedef gsExprAssembler<>::space       space;
    typedef gsExprAssembler<>::solution    solution;

    // Elements used for numerical integration
    A.setIntegrationElements(dbasis);
    gsExprEvaluator<> ev(A);
    // Set the geometry map
    geometryMap G = A.getMap(mp);

    // Set the Target geometry map
    geometryMap GLeft = A.getMap(mpLeft);

    // Set pow for BFO method
    auto IGdim     = G.domainDim();
    // Set factor for BFO method
    auto gammaMAE = factorial(G.domainDim());

    // Set the discretization space
    space u = A.getSpace(dbasis);
    
    // Set the source term with respect to target geometry
    auto ff = A.getCoeff(f, G);

    //gsFunctionExpr<> sI("0.5*(x**2+y**2)+x*y",2);
    auto u_I = ev.getVariable(sN, G);

    // Solution vector and solution variable
    gsMatrix<> solVector;
    solution u_sol = A.getSolution(u, solVector);
    //! [Problem setup]

    //! [Solver loop]
    gsSparseSolver<>::CGDiagonal solver;

    gsInfo<< "(dot1=assembled, dot2=solved)\n"
        "\nDoFs: ";
    double setup_time(0), ma_time(0), slv_time(0);    
    gsStopwatch timer;
    timer.restart();
    //::::::::::::::::::::      mesh adaptation solver         :::::::::::::::::::::::::
    for (int r=0; r<=numRefine; ++r)
    {
        dbasis.uniformRefine();
        mp.uniformRefine();
        mpLeft.uniformRefine();
    }
    //Initialization of Fast diagonalization solver
    gsPatchPreconditionersCreator<double>::Poisson_FastDiag Poisson(dbasis.basis(0), bc_mae, A.options(), eps);

    u.setup(bc_mae, dirichlet::l2Projection, 0);
    // Compute the system matrix and right-hand side

    // Initialize the system :  identity mapping as initial guess
    A.initSystem();
    // Initialize the system : start Computing the conductivity coeffeicient ...
    // Compute the Neumann terms defined on physical space
    auto g_N = A.getBdrFunction(G);
    auto Neumann_Int{ev.integralBdrBc(bc_mae.get("Neumann"), g_N.tr() * nv(G) )};
    //... nromalisation of density function
    auto CoeffDensity{ev.integral((1.+IntensityMAE*ff.val())* meas(G))};
    auto CoeffConductivity{Neumann_Int/ev.integral(pow(pow(IGdim,IGdim)+gammaMAE * CoeffDensity/(1.+IntensityMAE*ff.val()), 1./IGdim) )};

    setup_time += timer.stop();

    gsInfo<< A.numDofs() <<std::flush;

    timer.restart();
    A.assemble(
    igrad(u, G) * igrad(u, G).tr() * meas(G) + eps * u *u.tr()* meas(G) //matrix
    ,
    u*  CoeffConductivity * (-1.)*pow(pow(IGdim,IGdim)+gammaMAE* CoeffDensity/(1.+IntensityMAE*ff.val()), 1./IGdim) * meas(G) //rhs vector
    );
    
    // Compute the Neumann terms defined on physical space
    //auto g_N = A.getBdrFunction(G);
    A.assembleBdr(bc_mae.get("Neumann"), u * g_N.tr() * nv(G) );
    A.assembleIfc(mp.interfaces(), u.left() * (u_I.tr() * nv(G.left())));
    A.assembleIfc(mp.interfaces(), u.right() * (u_I.tr() * nv(G.right())));

    ma_time += timer.stop();

    gsInfo<< "." <<std::flush;// Assemblying done

    timer.restart();
    // solver.compute( A.matrix() );
    // solVector = solver.solve(A.rhs());
    solVector = Poisson.solve(A.rhs());
    slv_time += timer.stop();

    gsInfo<< "." <<std::flush; // Linear solving done

    // Picard loop
    gsMatrix<> sv0; //
    solution u_lsol = A.getSolution(u, sv0);
    for(int ip{0}; ip<=maxIter; ++ip)
    {
        timer.restart();
        gsMultiPatch<> UU;
        u_sol.extract(UU);
        gsWrite(UU, "U_solution");
        auto u_s = A.getCoeff(UU);

        space v = A.getSpace(dbasis);
        gsMatrix<> vsolVector;
        solution v_sol = A.getSolution(v, vsolVector);
        A.initSystem(2);

        // Obtain control points for the gradient of Psi
        A.assemble( v * v.tr() , v * grad(u_s) );
        //gsInfo <<"rhs vec = " << A.rhs().size() << "\n";
        vsolVector = Poisson.L2ProjectVec(A.rhs());
        
        gsMultiPatch<> Psi;
        v_sol.extract(Psi);

        // ... correct boundary
        ProjectionNormalCPoints(Psi, mpLeft, true);
        Psi.addAutoBoundaries();
        Psi.computeTopology();
        geometryMap PP    = A.getMap(Psi);
        auto ff      = A.getCoeff(f, GLeft, PP);/// Exact composition

        // ...  0  dirichlet for boundaries
        sv0 = solVector;
        u.setup(bc_mae, dirichlet::l2Projection, 0);
    
        solution u_sol = A.getSolution(u, solVector);

        // Initialize the system
        A.initSystem();
        setup_time += timer.stop();

        //gsInfo<< A.numDofs() <<std::flush;

        timer.restart();
        // Compute the system matrix and right-hand side ... Monge-Ampere eqaution .....
        
        // .. update Coeffeicient of conductivity
        auto  ExprMAE     = pow( pow(div(PP).val(),IGdim) + gammaMAE*(CoeffDensity/(1.+IntensityMAE*ff.val()) - jac(PP).det()), 1./IGdim);
        auto IntegDensity = ev.integral(ExprMAE);
        CoeffConductivity = Neumann_Int/IntegDensity;
        //.. Assemble the Monge-Ampere equation rhs and matrix        
        A.assemble(
        grad(u) * grad(u).tr()  +  eps * u * u.tr()//matrix
        ,
        u * CoeffConductivity * (-1.) * ExprMAE  //rhs vector
        );
        //gsInfo << "End Assemnles \n";

        // Compute the Neumann terms defined on physical space
        auto g_N = A.getBdrFunction(G);
        A.assembleBdr(bc_mae.get("Neumann"), u * g_N.tr() * nv(G) );
        A.assembleIfc(mp.interfaces(), u.left() * (u_I.tr() * nv(G.left())));
        A.assembleIfc(mp.interfaces(), u.right() * (u_I.tr() * nv(G.right())));
        ma_time += timer.stop();

        gsInfo<< " ." <<std::flush;// Assemblying done

        timer.restart();
        solVector = Poisson.solve(A.rhs());
        slv_time += timer.stop();

        gsInfo<< "." <<std::flush; // Linear solving done

        // omp_set_dynamic(0);     // Explicitly disable dynamic teams
        // omp_set_num_threads(1); // Use these threads for later parallel regions

        auto l2errRes = math::sqrt(ev.integral( ( grad(u_lsol) - grad(u_sol) ).sqNorm()  ));
        auto L2MAERes = math::sqrt(ev.integral( pow( CoeffDensity - (1.+IntensityMAE*ff.val())*jac(PP).det(),2)  ));
        if ( l2errRes < tolPicard || ip == maxIter ){
            // ! end Picard loop
            gsInfo<< "\n Niter in Picard : " << ip << ".. H1 residual : "<<std::scientific<<l2errRes<< ".. L2 MAE residual : "<<std::scientific<<L2MAERes<<"\n";
            break; 
            } // 
    }//for loop
    // omp_set_dynamic(0);     // Explicitly disable dynamic teams
    // omp_set_num_threads(1); // Use these threads for later parallel regions
    gsMultiPatch<> UU;
    u_sol.extract(UU);
    gsWrite(UU, "U_solution");
    auto u_s = A.getCoeff(UU);

    //gsMultiBasis<> gbasis(dbasis);
    //gbasis.reduceContinuity(1);
    space v = A.getSpace(dbasis);
    gsMatrix<> vsolVector;
    solution v_sol = A.getSolution(v, vsolVector);
    A.initSystem(2);

    //gsVector<> pt(2); pt.setConstant(0.5);
    //ev.testEval( v, pt );
    //ev.testEval( igrad(u_sol,G), pt );

    // Obtain control points for the gradient of Psi
    A.assemble( v * v.tr() , v * igrad(u_s,G) );
    vsolVector = solver.compute(A.matrix()).solve(A.rhs());
    gsMultiPatch<> Psi, Psitp;
    v_sol.extract(Psitp);
    //... correct the boundary
    ProjectionNormalCPoints(Psitp, mpLeft, true);

    //::::::::::::::::::::    Compute the composition of geometry maps      :::::::::::::::::::::::::
    // Psi.addAutoBoundaries();
    geometryMap PPi = A.getMap(Psitp);
    auto  comp = PPi(mpLeft);
    A.initSystem(2);
    //Obtain control points for the gradient of mpLeft.comp(Psi)
    A.assemble( v * v.tr() , v * comp.tr() );// blocked by this one
    vsolVector = solver.compute(A.matrix()).solve(A.rhs());
    v_sol.extract(Psitp);
    ProjectionNormalCPoints(Psitp, mpLeft, false);
    Psitp.addAutoBoundaries();
    Psitp.computeTopology();
    gsInfo << "end of adaptive mapping computation\n" << Psitp<< "\n";

    /* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ###   Step 4: Define hierarchical adaptive mapping
     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
    for(size_t i =0; i<Psitp.nPatches(); ++i)
        Psi.addPatch(gsTHBSpline<2>( dynamic_cast<const gsTensorBSpline<2>&>(Psitp.patch(i)) ));
    Psi.addAutoBoundaries();
    Psi.computeTopology();
    gsWrite(Psi, "Psi_mapping");
    //#-++++++++++++++++++++++++ End of sharing part of any geometry------------------------------
    timer.stop();
    gsInfo<<"\n\nTotal time: "<< setup_time+ma_time+slv_time<<"\n";
    gsInfo<<"     Setup: "<< setup_time <<"\n";
    gsInfo<<"  Assembly: "<< ma_time    <<"\n";
    gsInfo<<"   Solving: "<< slv_time   <<"\n";
    //::::::::::::::::::::  Elasticity equation - (manufactured exact solution)         :::::::::::::::::::::::::
    if (true){
    // creating basis
    gsMultiBasis<> basis(Psi);
    gsVector<>   l2err(numLRefine+1);//h1err(numLRefine+1),
    gsVector<int>  DoFPDE(numLRefine+1);
    gsInfo<< "(dot1=assembled, dot2=solved)\n";
    double setup_time(0), ma_time(0), slv_time(0), err_time(0);
    //=============================================//
        // Setting loads and boundary conditions //
    //=============================================//

    gsFunctionExpr<> analyticalStresses("1-1/(x^2+y^2)*(3/2*cos(2*atan2(y,x)) + cos(4*atan2(y,x))) + 3/2/(x^2+y^2)^2*cos(4*atan2(y,x))",
                                        "-1/(x^2+y^2)*(1/2*cos(2*atan2(y,x)) - cos(4*atan2(y,x))) - 3/2/(x^2+y^2)^2*cos(4*atan2(y,x))",
                                        "-1/(x^2+y^2)*(1/2*sin(2*atan2(y,x)) + sin(4*atan2(y,x))) + 3/2/(x^2+y^2)^2*sin(4*atan2(y,x))",2);
    // boundary load neumann BC
    gsFunctionExpr<> traction("(-1+1/(x^2+y^2)*(3/2*cos(2*atan2(y,x)) + cos(4*atan2(y,x))) - 3/2/(x^2+y^2)^2*cos(4*atan2(y,x))) * (x==-4) +"
                              "(-1/(x^2+y^2)*(1/2*sin(2*atan2(y,x)) + sin(4*atan2(y,x))) + 3/2/(x^2+y^2)^2*sin(4*atan2(y,x))) * (y==4)",
                              "(1/(x^2+y^2)*(1/2*sin(2*atan2(y,x)) + sin(4*atan2(y,x))) - 3/2/(x^2+y^2)^2*sin(4*atan2(y,x))) * (x==-4) +"
                              "(-1/(x^2+y^2)*(1/2*cos(2*atan2(y,x)) - cos(4*atan2(y,x))) - 3/2/(x^2+y^2)^2*cos(4*atan2(y,x))) * (y==4)",2);
    // material parameters
    real_t youngsModulus = 1.0e3;
    real_t poissonsRatio = 0.3;

    // boundary conditions
    gsBoundaryConditions<> bcInfo;
    bcInfo.addCondition(0,boundary::north,condition_type::neumann,&traction);
    bcInfo.addCondition(0,boundary::west,condition_type::dirichlet,nullptr,1); // last number is a component (coordinate) number
    bcInfo.addCondition(0,boundary::east,condition_type::dirichlet,nullptr,0);

    // source function, rhs
    gsConstantFunction<> g(0.,0.,2);

    //=============================================//
              // Assembling & solving //
    //=============================================//

    // creating assembler
    gsElasticityAssembler<real_t> assembler(Psi,basis,bcInfo,g);
    assembler.options().setReal("YoungsModulus",youngsModulus);
    assembler.options().setReal("PoissonsRatio",poissonsRatio);
    gsInfo<<"Assembling...\n";
    gsStopwatch clock;
    clock.restart();
    assembler.assemble();
    gsInfo << "Assembled a system (matrix and load vector) with "
           << assembler.numDofs() << " dofs in " << clock.stop() << "s.\n";

    gsInfo << "Solving...\n";
    clock.restart();

#ifdef GISMO_WITH_PARDISO
    gsSparseSolver<>::PardisoLLT solver(assembler.matrix());
    gsVector<> solVector = solver.solve(assembler.rhs());
    gsInfo << "Solved the system with PardisoLDLT solver in " << clock.stop() <<"s.\n";
#else
    gsSparseSolver<>::SimplicialLDLT solver(assembler.matrix());
    gsVector<> solVector = solver.solve(assembler.rhs());
    gsInfo << "Solved the system with EigenLDLT solver in " << clock.stop() <<"s.\n";
#endif

    //=============================================//
                      // Output //
    //=============================================//

    // constructing displacement as an IGA function
    gsMultiPatch<> u_disp;
    assembler.constructSolution(solVector,assembler.allFixedDofs(),u_disp);
    // constructing stress tensor
    gsPiecewiseFunction<> stresses;
    assembler.constructCauchyStresses(u_disp,stresses,stress_components::all_2D_vector);

    // constructing an IGA field (geometry + solution) for displacement
    gsField<> solutionField(assembler.patches(),u_disp);
    // constructing an IGA field (geometry + solution) for stresses
    gsField<> stressField(assembler.patches(),stresses,true);
    // analytical stresses
    gsField<> analyticalStressField(assembler.patches(),analyticalStresses,false);
    // --------------- adaptive refinement ---------------
    // Specify cell-marking strategy... 
    MarkingStrategy adaptRefCrit = PUCA;
    //MarkingStrategy adaptRefCrit = GARU;
    //MarkingStrategy adaptRefCrit = errorFraction;
    // Elements used for numerical integration

    gsExprEvaluator<> ev;
    ev.setIntegrationElements(assembler.multiBasis());
    gsExprEvaluator<>::geometryMap PP = ev.getMap(Psi);
    auto sigm_ex   = ev.getVariable(analyticalStresses, PP);
    auto ff_GPsi   = ev.getVariable(f, PP);
    //... error computation
    auto istress = ev.getVariable(stressField.fields());
        // omp_set_num_threads(1); // Use these threads for later parallel regions
    DoFPDE[0] = assembler.numDofs();
    timer.restart();
    l2err[0]= math::sqrt( ev.integral( ( sigm_ex - istress).sqNorm() * meas(PP) ));
    //h1err[0]= math::sqrt(ev.max( (sigm_ex - istress).sqNorm()));
    // eval stress at the top of the circular cut
    gsMatrix<> A(2,1);
    A << 1.,0.; // parametric coordinates for the isogeometric solution
    gsMatrix<> res;
    stresses.piece(0).eval_into(A,res);
    A << 0., 1.; // spatial coordinates for the analytical solution
    gsMatrix<> analytical;
    analyticalStresses.eval_into(A,analytical);
    gsInfo << "XX-stress at the top of the circle: " << res.at(0) << " (computed), " << analytical.at(0) << " (analytical)\n";
    gsInfo << "YY-stress at the top of the circle: " << res.at(1) << " (computed), " << analytical.at(1) << " (analytical)\n";
    gsInfo << "XY-stress at the top of the circle: " << res.at(2) << " (computed), " << analytical.at(2) << " (analytical)\n";
    for (int r=1; r<=numLRefine; ++r)
    {
        // if(r < numLRefine){
        //! [beginRefLoop]
            gsInfo << "====== Loop " << r << " of "
                    <<numLRefine<< " ====adapt Parameter ="<< adaptRefParam << " ======" << "\n";
            // --------------- error estimation/computation ---------------
            // Compute the error indicators
            ev.integralElWise( ff_GPsi );

            const std::vector<real_t> eltErrs  = ev.elementwise();
            //! [errorComputation]

            //! [adaptRefinementPart]
            // Mark elements for refinement, based on the computed local errors and
            // the refinement-criterion and -parameter.
            std::vector<bool> elMarked( eltErrs.size() );
            gsMarkElementsForRef( eltErrs, adaptRefCrit, adaptRefParam, elMarked);
            gsInfo <<"Marked "<< std::count(elMarked.begin(), elMarked.end(), true) <<" elements.\n";

            // Refine the marked elements with a 1-ring of cells around marked elements
            gsRefineMarkedElements( basis, elMarked, NumArMarEl);
            gsRefineMarkedElements( assembler.multiBasis(), elMarked, NumArMarEl);
            // assembler.multiBasis().repairInterfaces( Psi.interfaces() );
            gsRefineMarkedElements( Psi, elMarked, NumArMarEl);
            
            if (r%2==0)
            NumArMarEl = NumArMarEl + FactRefPar;
            // }
        //=============================================//
                // Assembling & solving //
        //=============================================//
        // basis.uniformRefine();
        // Psi.uniformRefine();
        // creating assembler
        gsElasticityAssembler<real_t> assembler(Psi,basis,bcInfo,g);
        assembler.options().setReal("YoungsModulus",youngsModulus);
        assembler.options().setReal("PoissonsRatio",poissonsRatio);
        gsInfo<<"Assembling...\n";
        gsStopwatch clock;
        clock.restart();
        assembler.assemble();
        gsInfo << "Assembled a system (matrix and load vector) with "
            << assembler.numDofs() << " dofs in " << clock.stop() << "s.\n";

        gsInfo << "Solving...\n";
        clock.restart();

    #ifdef GISMO_WITH_PARDISO
        gsSparseSolver<>::PardisoLLT solver(assembler.matrix());
        gsVector<> solVector = solver.solve(assembler.rhs());
        gsInfo << "Solved the system with PardisoLDLT solver in " << clock.stop() <<"s.\n";
    #else
        gsSparseSolver<>::SimplicialLDLT solver(assembler.matrix());
        gsVector<> solVector = solver.solve(assembler.rhs());
        gsInfo << "Solved the system with EigenLDLT solver in " << clock.stop() <<"s.\n";
    #endif

        //=============================================//
                        // Output //
        //=============================================//
        // // constructing displacement as an IGA function
        gsMultiPatch<> u_disp;
        assembler.constructSolution(solVector,assembler.allFixedDofs(),u_disp);
        // constructing stress tensor
        gsPiecewiseFunction<> stresses;
        assembler.constructCauchyStresses(u_disp,stresses,stress_components::all_2D_vector);

        // constructing an IGA field (geometry + solution) for displacement
        gsField<> solutionField(assembler.patches(),u_disp);
        // constructing an IGA field (geometry + solution) for stresses
        gsField<> stressField(assembler.patches(),stresses,true);
        // analytical stresses
        gsField<> analyticalStressField(assembler.patches(),analyticalStresses,false);

        gsExprEvaluator<> ev;
        ev.setIntegrationElements(assembler.multiBasis());
        gsExprEvaluator<>::geometryMap PP = ev.getMap(Psi);
        auto sigm_ex   = ev.getVariable(analyticalStresses, PP);
        auto ff_GPsi   = ev.getVariable(f, PP);
        //... error computation
        auto istress = ev.getVariable(stressField.fields());
            // omp_set_num_threads(1); // Use these threads for later parallel regions
        DoFPDE[r] = assembler.numDofs();
        timer.restart();
        l2err[r]  = math::sqrt( ev.integral( ( sigm_ex - istress).sqNorm() * meas(PP) ));
        //h1err[0]= math::sqrt(ev.max( (sigm_ex - istress).sqNorm()));
        // eval stress at the top of the circular cut
        gsMatrix<> A(2,1);
        A << 1.,0.; // parametric coordinates for the isogeometric solution
        gsMatrix<> res;
        stresses.piece(0).eval_into(A,res);
        A << 0., 1.; // spatial coordinates for the analytical solution
        gsMatrix<> analytical;
        analyticalStresses.eval_into(A,analytical);
        gsInfo << "XX-stress at the top of the circle: " << res.at(0) << " (computed), " << analytical.at(0) << " (analytical)\n";
        gsInfo << "YY-stress at the top of the circle: " << res.at(1) << " (computed), " << analytical.at(1) << " (analytical)\n";
        gsInfo << "XY-stress at the top of the circle: " << res.at(2) << " (computed), " << analytical.at(2) << " (analytical)\n";
        assembler.refresh();
    }
    //! [Solver loop]    


    timer.stop();
    gsInfo<<"\n\nTotal time: "<< setup_time+ma_time+slv_time+err_time <<"\n";
    gsInfo<<"     Setup: "<< setup_time <<"\n";
    gsInfo<<"  Assembly: "<< ma_time    <<"\n";
    gsInfo<<"   Solving: "<< slv_time   <<"\n";
    gsInfo<<"     Norms: "<< err_time   <<"\n";
    

    //! [Error and convergence rates]
    gsInfo<< "\nDoF_PDE = "<<std::scientific<<DoFPDE.transpose()<<"\n";
    gsInfo<< "L2_error_x = "<<std::scientific<<std::setprecision(3)<<l2err.transpose()<<"\n";
    //gsInfo<< "L2_error_y= "<<std::scientific<<std::setprecision(3)<<h1err.transpose()<<"\n";

    if (errorsave)
    {
    // Assuming DoFPDE, l2err, and h1err are gsMatrix or similar types
    std::ofstream outFile("error_analysis.txt", std::ios::app); // Open file in append mode
    if (outFile.is_open())
    {
        outFile << "#DoF_PDE: " << adaptRefParam <<" "<< NumArMarEl <<" " << IntensityMAE << " \n"<< std::scientific << DoFPDE.transpose() << "\n";
        outFile << "#L2_error: \n" << std::scientific << std::setprecision(3) << l2err.transpose() << "\n";
        //outFile << "#H1_error: \n" << std::scientific << std::setprecision(3) << h1err.transpose() << "\n";
        outFile << "#-------------------------------------------------------------------------------\n"; // Optional separator for readability
        outFile.close(); // Close the file after writing
    }
    else
    {
        gsInfo << "Error: Unable to open file for writing : error_analysis.txt.\n";
    }
    }
    else
    {
        gsInfo << "Errors are not saved. To save them, try with --errorsave.\n";
    }

    //! [Error and convergence rates]
    if (numLRefine>0 && false)
    {
        gsInfo<< "\nEoC (L2): " << std::fixed<<std::setprecision(2)
              <<  ( l2err.head(numRefine).array()  /
                   l2err.tail(numRefine).array() ).log().transpose() / std::log(2.0)
                   <<"\n";

        //gsInfo<<   "EoC (H1): "<< std::fixed<<std::setprecision(2)
        //     <<( h1err.head(numRefine).array() /
        //          h1err.tail(numRefine).array() ).log().transpose() / std::log(2.0) <<"\n";
    }
    //! [Error and convergence rates]
    //! [Export visualization in ParaView]
    if (plot)
    {
        // constructing an IGA field (geometry + solution) for displacement
        gsField<> solutionField(assembler.patches(),u_disp);
        // constructing an IGA field (geometry + solution) for stresses
        gsField<> stressField(assembler.patches(),stresses,true);
        // analytical stresses
        gsField<> analyticalStressField(assembler.patches(),analyticalStresses,false);
        //... density function
        gsField<> densityField(assembler.patches(),f,false);

        // creating a container to plot all fields to one Paraview file
        std::map<std::string,const gsField<> *> fields;
        fields["Deformation"] = &solutionField;
        fields["Stress"] = &stressField;
        fields["StressAnalytical"] = &analyticalStressField;
        fields["DensityAnalytical"] = &densityField;
        gsWriteParaviewMultiPhysics(fields,"infinit_plate",10000,plot);
        gsInfo << "Open \"infinit_plate.pvd\" in Paraview for visualization. Stress wiggles on the left side are caused by "
                  "a singularity in the parametrization.\n";
        gsFileManager::open("infinit_plate.pvd");

        gsInfo<<"Storing paraview...\n";
        // Write the computed solution to paraview files
        gsInfo<<"Making in Paraview...\n";
        gsParaviewCollection collection("ParaviewOutput/solution", &ev);
        collection.options().setSwitch("plotElements", true);
        collection.options().setSwitch("base64", false);
        collection.options().setInt("plotElements.resolution", 16);
        collection.options().setInt("numPoints", 10000);
        collection.newTimeStep(&Psi);
        // collection.addField(istress,"numerical stress");
        // collection.addField(jac(PP).det(), "Jacobian function");
        // collection.addField(sigm_ex, "exact stress");
        // collection.addField(ff_GPsi,"Density function");
        collection.saveTimeStep();
        collection.save();
        //------------------------------------
        gsInfo<<"Plotting in Paraview...\n";
        // Run paraview
        gsFileManager::open("ParaviewOutput/solution.pvd");
    }
    else
        gsInfo << "Done. No output created, re-run with --plot to get a ParaView "
                  "file containing the solution.\n";
    //! [Export visualization in ParaView]
    }

    return EXIT_SUCCESS;


}// end main