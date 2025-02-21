/** @file Monge_Ampere_example.cpp

    @brief Tutorial on how to use expression assembler to solve a non-linear Monge-Ampere equation

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris & M. BAHARI
*/

//! [Include namespace]
#include <gismo.h>
#include <fstream> // For file operations

using namespace gismo;
//! [Include namespace]

void ProjectionNormalCPoints(gsMultiPatch<>& Psi, int boxMaxNumber = 1){
    // Projection normal of control points (exact geometry)
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
};

int main(int argc, char *argv[])
{
    //! [Parse command line]
    bool plot             = false;
    index_t numRefine     = 3;
    index_t numLRefine    = 3;
    index_t numElevate    = 0;
    index_t maxIter       = 50;
    index_t NumArMarEl    = 0; // Number of ring of cells around marked elements
    double eps            = 1e-7; // pinalization coefficient
    double tolPicard      = 1e-8;
    double IntensityMAE   = 9.;
    bool export_b64       = false;
    bool errorsave        = false;
    // --------------- adaptive refinement ---------------
    // Specify cell-marking strategy... 
    index_t adaptRefCrit  = 2; // 1: GARU, 2: PUCA, 3: BULK, 4: PBULK
    real_t adaptRefParam  = 0.;     // ... adapt parameter.
    index_t FactRefPar    = 0;    // ... adapt parameter : adaptRefParam += FactRefPar in each iter
    // Specify the file path
    //std::string fn("pde/quart_annulus.xml");
    std::string fn("pde/circle.xml");


    gsCmdLine cmd("Tutorial on solving a non-linear Monge-Ampere problem.");
    cmd.addReal( "a", "adaptRefParam", "parameter for local h-refinement loops",  adaptRefParam );
    cmd.addInt( "c", "NumArMarEl", "augement NumArMarEl with such quantity in local h-refinement loops",  NumArMarEl );
    cmd.addString( "d", "file", "Input XML file data", fn );
    cmd.addInt( "e", "degreeElevation",
                "Number of degree elevation steps to perform before solving (0: equalize degree in all directions)", numElevate );
    cmd.addReal( "f", "IntensityMAE", "Intensity of density function",  IntensityMAE);
    cmd.addInt("i", "iter", "Maximum number of iterations for the iterative Picard", maxIter);
    cmd.addInt( "l", "numLRefine", "Number of local h-refinement loops",  numLRefine );
    cmd.addInt( "p", "FactRefPar", "augement adaptRefParam with such quantity in local h-refinement loops",  FactRefPar );
    cmd.addInt( "r", "adaptRefCrit", "Adaptive refinement criterion [1:GARU,2:PUCA,3:BULK,4:PBULK]",  adaptRefCrit );
    cmd.addInt( "u", "uniformRefine", "Number of Uniform h-refinement loops",  numRefine );
    cmd.addInt("quRule",
                 "Quadrature rule [1:GaussLegendre,2:GaussLobatto,3:PatchRule]",
                 1);

    cmd.addSwitch("plot", "Create a ParaView visualization file with the solution", plot);
    cmd.addSwitch("errorsave", "Create a file in ... and save errors", errorsave);


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

    //..... Test 1
    // convection coefficient:
    gsMatrix<> coeff_conv{1,2};
    // diffusion coefficient:
    //double coeff_diff = 1.;
    // reaction coefficient:
    //double coeff_reac = 0.;

    // source term: and manufactured solution
    gsFunctionExpr<> s;
    fd.getId(2000, s);
    gsFunctionExpr<> rhs;
    fd.getId(2001, rhs);



    //================================== -------------------------- =====================================
    //==================================   Independent of the geometry ==================================
    //================================== This part is fixed for all    ==================================
    //================================== -------------------------- =====================================
    // .... one single patch
    gsMultiPatch<> mp = gsNurbsCreator<>::BSplineSquareGrid(1,1,1, 0.0, 0.0);
    //Get all interfaces and boundaries:
    mp.degreeElevate(numElevate);
    mp.computeTopology();
    //mp.addAutoBoundaries();

    // Manufactured identity mapping
    gsFunctionExpr<> sN("x","y",2);
    // Right-hand side function : Analytical density function rho_1
    // Load the file
    gsFunctionExpr<> f;
    fd.getId(2003, f);
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
    gsMultiBasis<double> dbasis(mpLeft, true);//true: poly-splines (not NURBS)
    // TODO : build Hdiv solver
    //gsMultiBasis<double> Hdivbasis(mpLeft, true);//true: poly-splines (not NURBS)
    //Hdivbasis.degreeElevate(1);

    gsInfo << "Patches: "<< mpLeft.nPatches() <<", degree: "<< dbasis.minCwiseDegree() <<"\n";
#ifdef _OPENMP
    gsInfo<< "Available threads: "<< omp_get_max_threads() <<"\n";
#endif
    //! [Refinement]

    //! [Problem setup]
    gsExprAssembler<> A(1,1);
    // A.options().setReal("quA", 2.0);
    // A.options().setInt("quB", 2);
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

    // Set pow for BFO method dim in parameteric domain
    auto IGdim     = G.domainDim();

    // Set dimension of target geometry
    // auto ITdim     = mpLeft.geoDim();

    // Set factor for BFO method
    auto gammaMAE = factorial(G.domainDim());

    // Set the discretization space
    space u = A.getSpace(dbasis);

    // Set the source term with respect to target geometry
    auto ff = A.getCoeff(f, GLeft);

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
    timer.restart();
    gsPatchPreconditionersCreator<double>::Poisson_FastDiag Poisson(dbasis.basis(0), bc_mae, A.options(), eps);
    slv_time += timer.stop();

    // ......... INITIALIZE THE SYSTEM BY COMPUTIONG A Appr-DENSITY IN UNIT-SQUARE .........
    // Solution vector and solution variable
    gsMatrix<> densityVector;
    solution density_sol = A.getSolution(u, densityVector);
    u.setup(bc_mae, dirichlet::l2Projection, 0);
    A.initSystem();
    A.assemble(
    u *u.tr() //matrix
    ,
    u* ff.val()  //rhs vector
    );
    densityVector = Poisson.L2ProjectScalar(A.rhs());
    gsMultiPatch<> density;
    density_sol.extract(density);
    auto rho = A.getCoeff(density, G);
    // ... manipulation of density function
    auto empldensity = (ev.max(abs(rho.val()))-ev.min(abs(rho.val())));
    double  int_uh_0 = 0.;
    double  int_uh_1  = 1.;
    if (empldensity < 1e-5|| IntensityMAE <= 1. )
    {
        gsInfo << "Density function is constant in the domain rho = 1.\n";
    }
    else{
        int_uh_0  = (IntensityMAE-1.)/empldensity;
        int_uh_1  = (1.*ev.max(abs(rho.val()))-IntensityMAE*ev.min(abs(rho.val())))/empldensity;
        gsInfo << "Density function is not constant in the domain\n";
    }
    gsInfo << "Density functio: min "<< ev.min(int_uh_0*abs(rho.val()) + int_uh_1)<<"/ max " << ev.max(int_uh_0*abs(rho.val()) + int_uh_1) << "\n";
    // ......... End initialization for density.........

    // ......... Start solving the Monge-Ampere equation .........
    u.setup(bc_mae, dirichlet::l2Projection, 0);
    // Compute the system matrix and right-hand side

    // Initialize the system :  identity mapping as initial guess
    A.initSystem();
    // Initialize the system : start Computing the conductivity coeffeicient ...
    // Compute the Neumann terms defined on physical space
    auto g_N = A.getBdrFunction(G);
    auto Neumann_Int{ev.integralBdrBc(bc_mae.get("Neumann"), g_N.tr() * nv(G) )};
    //... nromalisation of density function
    auto CoeffDensity{ev.integral((int_uh_0*abs(rho.val()) + int_uh_1))};
    auto CoeffConductivity{Neumann_Int/ev.integral(pow(pow(IGdim,IGdim)-gammaMAE+gammaMAE * CoeffDensity/(int_uh_0*abs(rho.val()) + int_uh_1), 1./IGdim) )};

    setup_time += timer.stop();

    timer.restart();
    A.assemble(
    grad(u) * grad(u).tr()  + eps * u *u.tr() //matrix
    ,
    u*  CoeffConductivity * (-1.)*pow(pow(IGdim,IGdim)-gammaMAE+gammaMAE* CoeffDensity/(int_uh_0*abs(rho.val()) + int_uh_1), 1./IGdim)  //rhs vector
    );

    // Compute the Neumann terms defined on physical space
    //auto g_N = A.getBdrFunction(G);
    A.assembleBdr(bc_mae.get("Neumann"), u * g_N.tr() * nv(G) );
    A.assembleIfc(mp.interfaces(), u.left() * (u_I.tr() * nv(G.left())));
    A.assembleIfc(mp.interfaces(), u.right() * (u_I.tr() * nv(G.right())));

    ma_time += timer.stop();

    gsInfo<< "." <<std::flush;// Assemblying done
    solVector = A.rhs();
    timer.restart();
    solVector = Poisson.solve(A.rhs());

    // solver.compute( A.matrix() );
    // solVector = solver.solve(A.rhs());
    slv_time += timer.stop();

    gsInfo<< "." <<std::flush; // Linear solving done
    //! [Solver loop]
    gsInfo<< A.numDofs() <<std::flush;
    // Picard loop
    gsVector<>  h1Res(maxIter+1), l2err(maxIter+1), Iter_mae(maxIter+1);
    gsMatrix<> sv0; //
    solution u_lsol = A.getSolution(u, sv0);
    for(int ip{0}; ip<=maxIter; ++ip)
    {
        timer.restart();
        gsMultiPatch<> UU;
        u_sol.extract(UU);
        auto u_s = A.getCoeff(UU);

        space v = A.getSpace(dbasis);
        gsMatrix<> vsolVector;
        solution v_sol = A.getSolution(v, vsolVector);
        A.initSystem(2);

        // Obtain control points for the gradient of Psi
        A.assemble( v * v.tr() , v * grad(u_s) );
        vsolVector = Poisson.L2ProjectVec(A.rhs());

        gsMultiPatch<> Psi, PsiLoc;
        v_sol.extract(Psi);
        v_sol.extract(PsiLoc);

        // ... correct boundary
        ProjectionNormalCPoints(Psi);
        Psi.addAutoBoundaries();
        Psi.computeTopology();
        geometryMap PP    = A.getMap(Psi);
        geometryMap PPrho = A.getMap(Psi);
        auto rho = PPrho(density);
        // auto rho = A.getCoeff(density, PP);

        // ...  0  dirichlet for boundaries
        sv0 = solVector;
        u.setup(bc_mae, dirichlet::l2Projection, 0);

        solution u_sol = A.getSolution(u, solVector);

        // Initialize the system
        A.initSystem();
        setup_time += timer.stop();

        timer.restart();
        // Compute the system matrix and right-hand side ... Monge-Ampere eqaution .....
        // .. update Coeffeicient of conductivity
        auto  ExprMAE     = pow( abs(pow(div(PP).val(),IGdim) - gammaMAE*jac(PP).det())+ gammaMAE*CoeffDensity/(int_uh_0*abs(rho.val()) + int_uh_1), 1./IGdim);
        // if (dbasis.minCwiseDegree() > 2)
        // auto  ExprMAE     = pow( abs(pow(lapl(u_sol).val(),IGdim) - gammaMAE*hess(u_sol).det())+ gammaMAE*CoeffDensity/(int_uh_0*abs(rho.val()) + int_uh_1), 1./IGdim);
        auto IntegDensity = ev.integral(ExprMAE);
        CoeffConductivity = Neumann_Int/IntegDensity;
        // MAE system
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
        auto L2MAERes = math::sqrt(ev.integral( pow( CoeffDensity - (int_uh_0*abs(rho.val()) + int_uh_1)*jac(PP).det(),2)  ));
        Iter_mae[ip]  = ip;
        h1Res[ip]     = l2errRes;// Compute the H1 residual
        l2err[ip]     = L2MAERes;// Compute the L2 error in MA equation
        if ( l2errRes < tolPicard || ip == maxIter ){
            // ! end Picard loop
            gsInfo<< "\n Niter in Picard : " << ip
                    << ".. H1 residual : "<<std::scientific<<l2errRes
                    << ".. L2 MAE residual : "<<std::scientific<<L2MAERes<<"\n";
            break;
            } //
    }//for loop
    // omp_set_dynamic(0);     // Explicitly disable dynamic teams
    // omp_set_num_threads(1); // Use these threads for later parallel regions
    gsMultiPatch<> UU;
    u_sol.extract(UU);
    //gsWrite(UU, "U_solution");
    auto u_s = A.getCoeff(UU);

    //gsMultiBasis<> gbasis(dbasis);
    //gbasis.reduceContinuity(1);
    space v = A.getSpace(dbasis);
    gsMatrix<> vsolVector;
    solution v_sol = A.getSolution(v, vsolVector);
    A.initSystem(2);
    // Obtain control points for the gradient of Psi
    A.assemble( v * v.tr() , v * igrad(u_s,G) );
    vsolVector = solver.compute(A.matrix()).solve(A.rhs());
    gsMultiPatch<> Psi, Psitp;
    v_sol.extract(Psitp);
    //... correct the boundary
    ProjectionNormalCPoints(Psitp);

    //::::::::::::::::::::    Compute the composition of geometry maps      :::::::::::::::::::::::::
    // Psi.addAutoBoundaries();
    geometryMap PP = A.getMap(Psitp);
    auto  comp = PP(mpLeft);
    A.initSystem(2);
    //Obtain control points for the gradient of mpLeft.comp(Psi)
    A.assemble( v * v.tr() , v * comp.tr() );// blocked by this one
    vsolVector = solver.compute(A.matrix()).solve(A.rhs());
    v_sol.extract(Psitp);
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
    //#-++++++++++++++++++++++++ End of sharing part of any geometry------------------------------
    timer.stop();
    gsInfo<<"\n\nTotal time: "<< setup_time+ma_time+slv_time<<"\n";
    gsInfo<<"     Setup: "<< setup_time <<"\n";
    gsInfo<<"  Assembly: "<< ma_time    <<"\n";
    gsInfo<<"   Solving: "<< slv_time   <<"\n";

    //::::::::::::::::::::   Poisson equation - (manufactured exact solution)         :::::::::::::::::::::::::
    if (true){
    gsBoundaryConditions<> bc;
    bc.setGeoMap(Psi);
    // For simplicity, set Dirichlet boundary conditions
    for ( gsMultiPatch<>::const_biterator
            bit = mp.bBegin(); bit != mp.bEnd(); ++bit)
    {
       bc.addCondition( *bit, condition_type::dirichlet, &s,0, false);
    }
    gsInfo<<"The PDE domain is "<< Psi.detail() << "\n";
    gsInfo<<"Source function is "<< rhs << "\n";
    gsInfo<<"Boundary conditions:\n"<< bc <<"\n";


    dbasis.clear();
    gsMultiBasis<> dbasis(Psi, true);//true: poly-splines (not NURBS)

    geometryMap PP = A.getMap(Psi);
    auto ff_GPsi   = A.getCoeff(f, PP);
    
    // Elements used for numerical integration
    A.setIntegrationElements(dbasis);
    gsExprEvaluator<> ev(A);
    // Set the discretization space // different boundary condition !
    space ru = A.getSpace(dbasis);

    // Set the source term for Poisson equation
    auto SFunc      = A.getCoeff(rhs, PP);

    // Recover manufactured solution for Poisson equation
    auto u_ex       = ev.getVariable(s, PP);

    // Solution vector and solution variable
    gsMatrix<> rsolVector;
    solution ru_sol = A.getSolution(ru, rsolVector);

    gsVector<>  h1err(numLRefine+1), l2err(numLRefine+1);
    gsVector<int>  DoFPDE(numLRefine+1);
    gsInfo<< "(dot1=assembled, dot2=solved)\n";
    double setup_time(0), ma_time(0), slv_time(0), err_time(0);
    for (int r=0; r<=numLRefine; ++r)
    {
        //::::::::::::::::::::   Poisson equation - (manufactured exact solution)         :::::::::::::::::::::::::
        ru.setup(bc, dirichlet::l2Projection, 0);

        // Compute the system matrix and right-hand side
        // Initialize the system
        A.initSystem();
        setup_time += timer.stop();

        gsInfo<< "Solving PDEs " <<std::flush;
        gsInfo<< A.numDofs() <<std::flush;
        
        //auto h_Tau =  m_h/(2.*coeff_conv.squaredNorm()+m_h);

        timer.restart();
        A.assemble(
        igrad(ru, PP) * igrad(ru, PP).tr() * meas(PP) //matrix
        ,
        ru * SFunc * meas(PP) //rhs vector
        );

        ma_time += timer.stop();

        // gsDebugVar(A.matrix().toDense());
        // gsDebugVar(A.rhs().transpose()   );

        gsInfo<< "." <<std::flush;// Assemblying done

        timer.restart();
        solver.compute( A.matrix() );
        rsolVector = solver.solve(A.rhs());

        slv_time += timer.stop();

        gsInfo<< "." <<std::flush; // Linear solving done

        // omp_set_dynamic(0);     // Explicitly disable dynamic teams
        // omp_set_num_threads(1); // Use these threads for later parallel regions
        DoFPDE[r] = A.numDofs();

        timer.restart();
        l2err[r]= math::sqrt( ev.integral( (u_ex - ru_sol).sqNorm() * meas(PP) ) );

        h1err[r]= math::sqrt(ev.integral( ( igrad(u_ex) - igrad(ru_sol,PP) ).sqNorm() * meas(PP) ));
        err_time += timer.stop();
        gsInfo<< ". " <<std::flush; // Error computations done
        if(r < numLRefine){
        //! [beginRefLoop]
            gsInfo << "====== Loop " << r << " of "
                    <<numLRefine<< " ====adapt Parameter ="<< adaptRefParam << " ======" << "\n";
            // --------------- error estimation/computation ---------------
            // Get the element-wise norms.
            //ev.integralElWise( (  ilapl(ru_sol, PP)+ SFunc ).sqNorm() );
            //if (IntensityMAE > 1.)
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
            gsRefineMarkedElements( dbasis, elMarked, NumArMarEl);
            gsRefineMarkedElements( Psi, elMarked, NumArMarEl);
            if (r%2==0)
            NumArMarEl = NumArMarEl + FactRefPar;
            }
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
    gsInfo<< "L2_error = "<<std::scientific<<std::setprecision(3)<<l2err.transpose()<<"\n";
    gsInfo<< "H1_error= "<<std::scientific<<std::setprecision(3)<<h1err.transpose()<<"\n";

    if (errorsave)
    {
    // Assuming DoFPDE, l2err, and h1err are gsMatrix or similar types
    std::ofstream outFile("error_analysis.txt", std::ios::app); // Open file in append mode
    if (outFile.is_open())
    {
        outFile << "#DoF_PDE: " << adaptRefParam <<" "<< NumArMarEl <<" " << IntensityMAE << " \n"<< std::scientific << DoFPDE.transpose() << "\n";
        outFile << "#L2_error: \n" << std::scientific << std::setprecision(3) << l2err.transpose() << "\n";
        outFile << "#H1_error: \n" << std::scientific << std::setprecision(3) << h1err.transpose() << "\n";
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

        gsInfo<<   "EoC (H1): "<< std::fixed<<std::setprecision(2)
              <<( h1err.head(numRefine).array() /
                  h1err.tail(numRefine).array() ).log().transpose() / std::log(2.0) <<"\n";
    }
    //! [Error and convergence rates]
    //! [Export visualization in ParaView]
    if (plot)
    {
        gsInfo<<"Plotting in Paraview...\n";
        gsParaviewCollection collection("ParaviewOutput/solution", &ev);
        collection.options().setSwitch("plotElements", true);
        collection.options().setSwitch("base64", export_b64);
        collection.options().setInt("plotElements.resolution", 16);
        collection.options().setInt("numPoints", 10000);
        collection.newTimeStep(&Psi);
        collection.addField(ru_sol,"numerical solution");
        collection.addField(igrad(ru_sol,PP),"gradient_numerical solution");
        collection.addField(jac(PP).det(), "Jacobian function");
        collection.addField(u_ex, "exact solution");
        collection.addField(ff_GPsi,"Density function");
        collection.saveTimeStep();
        collection.save();
        gsFileManager::open("ParaviewOutput/solution.pvd");
    }
    else
        gsInfo << "Done. No output created, re-run with --plot to get a ParaView "
                  "file containing the solution.\n";
    //! [Export visualization in ParaView]
    }

    return EXIT_SUCCESS;


}// end main