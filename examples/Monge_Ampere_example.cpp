/** @file Monge_Ampere_example.cpp

    @brief Tutorial on how to use expression assembler to solve a non-linear Monge-Ampere equation

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): M. BAHARI
*/

//! [Include namespace]
#include <gismo.h>
#include <fstream>  // For file operations

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
    bool plot           = false; // plot solution or not
    index_t numRefine   = 4;     // number of uniform refinement
    index_t numIRefine  = 3;     // number for initial refinement
    index_t numElevate  = 1;     // number of degree elevation
    index_t maxIter     = 50;    // maximum number of Picard iterations
    double eps          = 1e-6;  // pinalization coefficient
    double tolPicard    = 1e-12; // tolerance for Picard iterations
    double IntensityMAE = 9.;    // Intensity of density function
    bool export_b64     = false; // export in base64 format
    bool last           = false; // solve solely for the last level of h-refinement
    bool onlyH          = false; // only local h-refinement after initial refinement
    bool locrh          = true;  // local h-refinement or not

    // --------------- adaptive refinement ---------------
    // Specify cell-marking strategy...
    MarkingStrategy adaptRefCrit = PUCA;
    // MarkingStrategy adaptRefCrit = GARU;
    //MarkingStrategy adaptRefCrit = errorFraction;
    real_t adaptRefParam = 0.75;

    // Specify the file path
    //std::string fn("pde/quart_annulus.xml");
    //std::string fn("pde/mhd.xml");
    //std::string fn("pde/infinit_plate.xml");
    std::string fn("pde/circle.xml");
    //std::string fn("surfaces/egg.xml");
    //std::string fn("domain2d/lake.xml");

    gsCmdLine cmd("Tutorial on solving a non-linear Monge-Ampere problem.");
    cmd.addInt("i",       "iter",            "Maximum number of iterations for Picard",    maxIter);
    cmd.addInt( "e",      "degreeElevation", "Number of degree elevation ",                numElevate );
    cmd.addInt( "u",      "uniformRefine",    "Number of Uniform h-refinement loops",      numRefine );
    cmd.addInt( "g",      "G-uniformRefine",    "Number of Uniform global-refinement",     numIRefine);
    cmd.addString("d",    "file",             "Input XML file data",                       fn );
    cmd.addReal( "l",     "adaptRefParam",     "percentage of local h-refinement loops",   adaptRefParam );
    cmd.addReal( "f",     "IntensityMAE",      "Intensity of density function",            IntensityMAE);
    cmd.addReal( "t",     "tolPicard",         "the tolerance for Picard iterations",      tolPicard);
    cmd.addInt("quRule",  "Quadrature rule [1:GaussLegendre,2:GaussLobatto,3:PatchRule]",  1);
    cmd.addSwitch("plot", "Create a ParaView visualization file with the solution",        plot);
    cmd.addSwitch("last", "Solve solely for the last level of h-refinement",               last);
    cmd.addSwitch("b64",  "Export the solution in base64 format",                         export_b64);
    cmd.addSwitch("onlyH","Perform only local h-refinement after initial refinement",      onlyH);
    cmd.addSwitch("locrh","Perform local h-refinement",                                    locrh);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    // Load the file
    gsFileData<> fd(fn);
    gsInfo << "Loaded file " << fd.lastPath() << "\n";
    // Create a gsMultipatch and add the loaded geometry
    gsMultiPatch<> mpLeft; //= gsNurbsCreator<>::BSplineSquareGrid(1,1,1, 0.0, 0.0);
    fd.getId(1,mpLeft);
    // Elevate and p-refine the basis to order p + numElevate
    // where p is the highest degree in the bases
    mpLeft.degreeElevate(numElevate);
    mpLeft.computeTopology();
    gsInfo << "INFO IN PARAMETRIC DOMAIN "<< mpLeft.dim() << mpLeft.parDim() <<"\n";
    //================================== -------------------------- =====================================
    //==================================   Independent of the geometry ==================================
    //==================================    This part is fixed for all ==================================
    //================================== -------------------------- =====================================
    // .... one single patch
    gsMultiPatch<> mp = gsNurbsCreator<>::BSplineSquareGrid(1,1,1, 0.0, 0.0);
    //Get all interfaces and boundaries:
    mp.degreeElevate(numElevate);
    mp.computeTopology();
    //mp.addAutoBoundaries();
    gsMultiPatch<> Psi;
    if (locrh){
    for(size_t i =0; i<mp.nPatches(); ++i)
        Psi.addPatch(gsTHBSpline<2>( dynamic_cast<const gsTensorBSpline<2>&>(mp.patch(i)) ));
    }
    else
        Psi = mp;
    Psi.addAutoBoundaries();
    Psi.computeTopology();
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
    gsMultiBasis<double> dbasis(Psi, true);//true: poly-splines (not NURBS)

    //.. initial refinement to have the same number of elements as the adaptive geometry
    for (int r=0; r<numIRefine; ++r)
    {
        dbasis.uniformRefine();
        Psi.uniformRefine();
    }
    if( last){
    for (int r=0; r<numRefine; ++r)
    {
        dbasis.uniformRefine();
        Psi.uniformRefine();
    }
    numRefine = 0;
    }

    gsInfo << "Patches: "<< mp.nPatches() <<", degree: "<< dbasis.minCwiseDegree() <<" nbasis: " << dbasis.basis(0).size() <<"\n";
#ifdef _OPENMP
    gsInfo<< "Available threads: "<< omp_get_max_threads() <<"\n";
#endif
    //! [Refinement]

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
    // Set the square geometry map
    geometryMap G     = A.getMap(mp);

    // Set the Physical geometry map
    geometryMap GLeft = A.getMap(mpLeft);

    // Set the Target geometry adaptive map
    geometryMap PP    = A.getMap(Psi);

    // Set pow for BFO method dim in parameteric domain
    auto IGdim        = G.domainDim();

    // Set factor for BFO method
    auto gammaMAE     = factorial(G.domainDim());

    // Set the discretization space
    space u           = A.getSpace(dbasis);

    //... solution vector and solution variable for gradient of potential function
    gsMatrix<> vsolVector;
    solution v_sol    = A.getSolution(u, vsolVector);

    // Set the source term with respect to target geometry
    auto ff           = A.getCoeff(f, GLeft);
    // Set the source term with respect to target geometry
    auto fP           = A.getCoeff(f,  GLeft, PP);

    //gsFunctionExpr<> sI("0.5*(x**2+y**2)+x*y",2);
    auto u_I          = ev.getVariable(sN, G);

    // Solution vector and solution variable
    gsMatrix<> solVector;
    solution u_sol    = A.getSolution(u, solVector);
    //! [Problem setup]

    gsInfo<< "(dot1=assembled, dot2=solved)\n"
        "\nDoFs: ";
    double setup_time(0), ma_time(0), slv_time(0);
    gsStopwatch timer;
    timer.restart();
    //::::::::::::::::::::      mesh adaptation solver         :::::::::::::::::::::::::
    gsVector<>  DoFs_sys(numRefine+1), Solv_time(numRefine+1);;
    gsVector<>  l2Gdisp(numRefine+1), error_mae(numRefine+1); //    
    // gsPatchPreconditionersCreator<double>::Poisson_FastDiag Poisson(dbasis.basis(0), bc_mae, A.options(), eps);
    gsSparseSolver<>::CGDiagonal solver; // exact solver
    gsSparseSolver<>::CGDiagonal Ssolver;// relaxation solver
    // Ssolver.setMaxIterations(20);
    Ssolver.setTolerance(1e-4);
    timer.restart();
    // -------------------- for projection --------------------
    u.setup(bc_mae, dirichlet::l2Projection, 0);
    A.initSystem();
    A.assemble(u *u.tr());//matrix
    auto MProj = A.matrix();
    solver.compute( MProj);
    // -------------------- for MAE system ----------------
    u.setup(bc_mae, dirichlet::l2Projection, 0);
    A.initSystem();
    A.assemble(grad(u) * grad(u).tr() + eps * u *u.tr());//matrix
    auto MMAe = A.matrix();
    Ssolver.compute( MMAe );
    setup_time += timer.stop();

    // ......... INITIALIZE THE SYSTEM BY COMPUTIONG A Appr-DENSITY IN UNIT-SQUARE .........
    // Solution vector and solution variable
    gsMatrix<> densityVector;
    solution density_sol = A.getSolution(u, densityVector);
    u.setup(bc_mae, dirichlet::l2Projection, 0);
    timer.restart();
    A.initSystem();
    A.assemble(u* ff.val());//rhs vector
    ma_time += timer.stop();
    timer.restart();
    // densityVector = Poisson.L2ProjectScalar(A.rhs());
    densityVector = solver.solve(A.rhs());
    slv_time += timer.stop();
    for(index_t ds= 0; ds<densityVector.size(); ++ds){
        if (densityVector.coeff(ds)<=0.)
            densityVector.coeffRef(ds) = 0.;
    }
    gsMultiPatch<> density, u_phi;
    density_sol.extract(density);
    // gsWrite(density, "density");
    // std::string fr("pde/density_hand.xml");
    // gsMultiPatch<> density;
    // gsFileData<> fdr(fr);
    // fdr.getId(1, density);
    auto rho = A.getCoeff(density, PP);
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
        gsInfo << "Density functio: min "<< ev.min(int_uh_0*abs(rho.val()) + int_uh_1)<<"/ max " << ev.max(int_uh_0*abs(rho.val()) + int_uh_1) << "\n";
    }

    // ......... End initialization for density.........
        
    // ......... Start solving the Monge-Ampere equation .........
    // u.setup(bc_mae, dirichlet::l2Projection, 0);
    // Compute the system matrix and right-hand side

    // Initialize the system :  identity mapping as initial guess
    A.initSystem();
    // Initialize the system : start Computing the conductivity coeffeicient ...
    // Compute the Neumann terms defined on physical space
    auto g_N      = A.getBdrFunction(G);
    auto Neumann_Int{ev.integralBdrBc(bc_mae.get("Neumann"), g_N.tr() * nv(G) )};
    //... normalisation of a density function
    auto CoeffDensity{ev.integral((int_uh_0*abs(rho.val()) + int_uh_1))};
    auto ExprMAE     = pow(pow(IGdim,IGdim)-gammaMAE+gammaMAE * CoeffDensity/(int_uh_0*abs(rho.val()) + int_uh_1), 1./IGdim);
    auto CoeffConductivity{Neumann_Int/ev.integral(ExprMAE)};

    A.assemble(u*  CoeffConductivity * (-1.)*ExprMAE//rhs vector
    );
    // Compute the Neumann terms defined on physical space
    A.assembleBdr(bc_mae.get("Neumann"), u * g_N.tr() * nv(G) );
    gsInfo<< "." <<std::flush;// Assemblying done
    timer.restart();
    // solVector = Poisson.solve(A.rhs());
    solVector = Ssolver.solve(A.rhs());
    slv_time += timer.stop();
    u_sol.extract(u_phi);
    gsInfo<< "." <<std::flush; // Linear solving done

    //! [Solver loop]
    gsInfo<< A.numDofs() <<std::flush;

    // ... compute the projection of gradient of potential function
    timer.restart();
    A.initSystem(IGdim);
    // Obtain control points for the gradient of Psi
    A.assemble(u * igrad(u_sol,G));
    ma_time += timer.stop();
    for(index_t Mp=0; Mp<mp.dim(); ++Mp){
    gsMultiPatch<> PsiPsitp_temp;
    timer.restart();
    // vsolVector = Poisson.L2ProjectScalar(A.rhs().col(Mp));
    vsolVector = solver.solve(A.rhs().col(Mp));
    slv_time += timer.stop();
    v_sol.extract(PsiPsitp_temp);
    Psi.patch(0).coefs().col(Mp) = PsiPsitp_temp.patch(0).coefs().col(0);
    }
    // // ... correct boundary
    ProjectionNormalCPoints(Psi);
    Psi.addAutoBoundaries();
    Psi.computeTopology();
    //! [beginRefLoop]
    gsInfo << "====== Loop without h-refinement ======\n";
    for (int r=0; r<=numRefine; ++r)
    {
        // initialize time counters
        setup_time=0.; ma_time=0.; slv_time=0.;
        // ......... End initialization for density.........

        if( r > 0 && locrh){
        //! [beginRefLoop]
        gsInfo << "====== Loop " << r << " of "
                <<numRefine<< " ====adapt Parameter ="<< adaptRefParam << " ======\n";
        // --------------- error estimation/computation ---------------
        // Get the element-wise norms.
        // ev.integralElWise(int_uh_0*rho.val()+int_uh_1 );
        // ev.integralElWise( pow( 1. - (int_uh_0*abs(rho.val()) + int_uh_1)*jac(PP).det()/CoeffDensity,2)   );
        ev.integralElWise( fP.val() );

        //! [errorComputation]
        const std::vector<real_t> eltErrs  = ev.elementwise();
        //! [adaptRefinementPart]
        // Mark elements for refinement, based on the computed local errors and
        // the refinement-criterion and -parameter.
        std::vector<bool> elMarked( eltErrs.size() );
        gsMarkElementsForRef( eltErrs, adaptRefCrit, adaptRefParam, elMarked);
        gsInfo <<"Marked "<< std::count(elMarked.begin(), elMarked.end(), true) <<" elements.\n";
        // Refine the marked elements with a 1-ring of cells around marked elements
        gsRefineMarkedElements( dbasis, elMarked, 0);
        gsRefineMarkedElements( Psi, elMarked, 0);
        gsRefineMarkedElements( u_phi, elMarked, 0);
        //gsRefineMarkedElements( density, elMarked, 0);
        }
        else if (r>0)
        {
            dbasis.uniformRefine();
            Psi.uniformRefine();
            u_phi.uniformRefine();
        }
        gsInfo<< "." <<std::flush; // h-refinement done
        // ... refinement
        if (r > 0 && onlyH){
            continue;
        }

        // -------------------- It solver for projection --------------------
        // gsPatchPreconditionersCreator<double>::Poisson_FastDiag Poisson(dbasis.basis(0), bc_mae, A.options(), eps);
        u.setup(bc_mae, dirichlet::l2Projection, 0);
        A.initSystem();
        gsInfo<< A.numDofs() <<std::flush;
        DoFs_sys[r] = A.numDofs();
        A.assemble(u *u.tr());//matrix
        auto MProj = A.matrix();
        solver.compute( MProj);
        gsInfo<< "." <<std::flush; // Initialize Projection Iterative solver done        
        // -------------------- It solver for MAE system ----------------
        u.setup(bc_mae, dirichlet::l2Projection, 0);
        A.initSystem();
        A.assemble(grad(u) * grad(u).tr() + eps * u *u.tr());//matrix
        auto MMAe = A.matrix();
        Ssolver.compute( MMAe );
        gsInfo<< "." <<std::flush; // Initialize Iterative solver done


        // ......... INITIALIZE THE SYSTEM BY COMPUTIONG A Appr-DENSITY IN UNIT-SQUARE .........
        u.setup(bc_mae, dirichlet::l2Projection, 0);
        timer.restart();
        A.initSystem();
        A.assemble(u* ff.val());//rhs vector
        ma_time += timer.stop();
        timer.restart();
        // densityVector = Poisson.L2ProjectScalar(A.rhs());
        densityVector = solver.solve(A.rhs());
        slv_time += timer.stop();
        for(index_t ds= 0; ds<densityVector.size(); ++ds){
            if (densityVector.coeff(ds)<=0.)
                densityVector.coeffRef(ds) = 0.;
        }
        density_sol.extract(density);
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
            std::cout << std::setprecision(1);
            gsInfo << ".rho~"<< ev.min(int_uh_0*abs(rho.val()) + int_uh_1)<<"/" << ev.max(int_uh_0*abs(rho.val()) + int_uh_1) << "."<<std::flush;
        }
        // .... Picard loop  .... ...........................
        auto  sv0 = u_phi.patch(0).coefs(); // initial guess set as last solution
        solution u_lsol = A.getSolution(u, sv0);

        for(int ip{0}; ip<=maxIter; ++ip)
        {
            // ...  0  dirichlet for boundaries
            timer.restart();
            //u.setup(bc_mae, dirichlet::l2Projection, 0);
            // Initialize the system
            A.initSystem();
            // .. update Coeffeicient of conductivity
            auto  ExprMAE     = pow( abs(pow(div(PP).val(),IGdim) - gammaMAE*jac(PP).det())+ gammaMAE*CoeffDensity/(int_uh_0*abs(rho.val()) + int_uh_1), 1./IGdim);
            auto IntegDensity = ev.integral(ExprMAE);
            CoeffConductivity = Neumann_Int/IntegDensity;
            // MAE system
            A.assemble(u * CoeffConductivity * (-1.) * ExprMAE);//rhs vector

            // Compute the Neumann terms defined on physical space
            auto g_N = A.getBdrFunction(G);
            A.assembleBdr(bc_mae.get("Neumann"), u * g_N.tr() * nv(G) );
            ma_time += timer.stop();

            gsInfo<< " ." <<std::flush;// Assemblying done
            timer.restart();
            // solVector = Poisson.solve(A.rhs());
            // solVector = Ssolver.solve(A.rhs());
            solVector = Ssolver.solveWithGuess(A.rhs(), sv0);
            slv_time += timer.stop();
            gsInfo<< "." <<std::flush; // Linear solving done

            // omp_set_dynamic(0);     // Explicitly disable dynamic teams
            // omp_set_num_threads(1); // Use these threads for later parallel regions
            auto h1Res  = math::sqrt(ev.integral( ( grad(u_sol)-grad(u_lsol)).sqNorm()  ));
            sv0         = solVector;

            if ( h1Res < tolPicard || ip == maxIter ){
                // ! end Picard loop
                gsInfo<< "\n Niter in Picard : " << ip
                        << ".. H1 residual : "<<std::scientific<<h1Res
                        << ".. min JAcobian : "<<ev.min(jac(PP).det())<<"..";
                break;
                } //
            timer.restart();
            A.initSystem(IGdim);
            // Obtain control points for the gradient of Psi
            A.assemble(u * grad(u_sol));
            ma_time += timer.stop();
            for(index_t Mp=0; Mp<mp.dim(); ++Mp){
            gsMultiPatch<> PsiPsitp_temp;
            timer.restart();
            // vsolVector = Poisson.L2ProjectScalar(A.rhs().col(Mp));
            // vsolVector = solver.solve(A.rhs().col(Mp));
            vsolVector = solver.solveWithGuess(A.rhs().col(Mp), Psi.patch(0).coefs().col(Mp));
            slv_time += timer.stop();
            v_sol.extract(PsiPsitp_temp);
            Psi.patch(0).coefs().col(Mp) = PsiPsitp_temp.patch(0).coefs().col(0);
            }
            // // ... correct boundary
            ProjectionNormalCPoints(Psi);
            Psi.addAutoBoundaries();
            Psi.computeTopology();
        }//for loop
        // omp_set_dynamic(0);     // Explicitly disable dynamic teams
        // omp_set_num_threads(1); // Use these threads for later parallel regions
        error_mae[r] = math::sqrt(ev.integral( pow( 1. - (int_uh_0*abs(rho.val()) + int_uh_1)*jac(PP).det()/CoeffDensity,2)  ));
        l2Gdisp[r]   = math::sqrt(ev.integral( ( PP-G).sqNorm()  ));
        Solv_time[r] = slv_time;// cpu time for a solver
    gsInfo<<"\n\nTotal time: "<< setup_time+ma_time+slv_time<<"\n";
    gsInfo<<"     Setup: "<< setup_time <<"\n";
    gsInfo<<"  Assembly: "<< ma_time    <<"\n";
    gsInfo<<"   Solving: "<< slv_time <<"\n";
    }
    timer.stop();

    //! [Error and convergence rates]
    gsInfo<< "itermax   :"<< IntensityMAE << "\n";
    gsInfo<< "\n#DoFs  : \n"<<std::scientific<<std::setprecision(3)<<DoFs_sys.transpose()<<"\n";
    gsInfo<< "#CPU-time: \n"<<std::scientific<<Solv_time.transpose()<<"\n";
    gsInfo<< "#L2 error: \n"<<std::scientific<<error_mae.transpose()<<"\n";
    gsInfo<< "#L2 disp : \n"<<std::scientific<<l2Gdisp.transpose()<<"\n";
    //! [Export visualization in ParaView]
    if (plot)
    {
        //::::::::::::::::::::    Compute the composition of geometry maps      :::::::::::::::::::::::::
        // -------------------- for projection --------------------
        u.setup(bc_mae, dirichlet::l2Projection, 0);
        A.initSystem();
        A.assemble(u *u.tr());//matrix
        auto MProj = A.matrix();
        solver.compute( MProj);
        // // //::::::::::::::::::::...
        auto comp  = A.getCoeff(mpLeft, PP);
        A.initSystem(mpLeft.geoDim());
        A.assemble(u * comp.tr() );// blocked by this one
        for(index_t Mp=0; Mp<mpLeft.geoDim(); ++Mp){
        gsMultiPatch<> PsiPsitp_temp;
        timer.restart();
        vsolVector = solver.solve(A.rhs().col(Mp));
        slv_time += timer.stop();
        v_sol.extract(PsiPsitp_temp);

        Psi.patch(0).coefs().col(Mp) = PsiPsitp_temp.patch(0).coefs().col(0);
        }

        if (mpLeft.dim() < mpLeft.geoDim() ){
            /// special case where a mpping is surface in three dimensions
            gsInfo << "surface in three dimensions\n";
            gsMultiPatch<> PsiPsitp_temp;
            vsolVector = solver.solve(A.rhs().col(2));
            v_sol.extract(PsiPsitp_temp);
            Psi.embed(3);
            Psi.patch(0).coefs().col(2) = PsiPsitp_temp.patch(0).coefs().col(0);
        }
        Psi.computeTopology();

        //::::::::::::::::::::      end       :::::::::::::::::::::::::
        gsInfo<<"Plotting in Paraview...\n";
        gsParaviewCollection collection("ParaviewOutput/solution", &ev);
        collection.options().setSwitch("plotElements", true);
        collection.options().setSwitch("base64", export_b64);
        collection.options().setInt("plotElements.resolution", 16);
        collection.options().setInt("numPoints", 10000);
        collection.newTimeStep(&Psi);
        collection.addField(fP, "density function");
        collection.addField(jac(PP).det(), "Jacobian function");
        collection.saveTimeStep();
        collection.save();

        gsFileManager::open("ParaviewOutput/solution.pvd");
    }
    else
        gsInfo << "Done. No output created, re-run with --plot to get a ParaView "
                  "file containing the solution.\n";
    //! [Export visualization in ParaView]

    return EXIT_SUCCESS;

}// end main