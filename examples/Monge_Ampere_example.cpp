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
    bool plot           = false;
    index_t numRefine   = 3;
    index_t numLRefine  = 3;
    index_t numElevate  = 1;
    index_t numrRefine  = -1; // number of composition bewteen adaptive mappings ()
    index_t maxIter     = 50;
    double eps          = 1e-6; // pinalization coefficient
    double tolPicard    = 1e-8;
    double IntensityMAE = 6.;
    bool export_b64     = false;
    // Specify the file path
    std::string fn("pde/quart_annulus.xml");
    //std::string fn("pde/infinit_plate.xml");
    //std::string fn("pde/circle.xml");
    //std::string fn("surfaces/egg.xml");
    //std::string fn("domain2d/lake.xml");

    gsCmdLine cmd("Tutorial on solving a non-linear Monge-Ampere problem.");
    cmd.addInt("i", "iter", "Maximum number of iterations for the iterative Picard", maxIter);
    cmd.addInt( "e", "degreeElevation",
                "Number of degree elevation steps to perform before solving (0: equalize degree in all directions)", numElevate );
    cmd.addInt( "u", "uniformRefine", "Number of Uniform h-refinement loops",  numRefine );
    cmd.addInt( "l", "numLRefine", "Number of local h-refinement loops",  numLRefine );
    cmd.addInt( "r", "numrRefine", "Number of local r-refinement compostion loops",  numrRefine);
    cmd.addString( "d", "file", "Input XML file data", fn );
    cmd.addInt("quRule",
                 "Quadrature rule [1:GaussLegendre,2:GaussLobatto,3:PatchRule]",
                 1);
    cmd.addSwitch("plot", "Create a ParaView visualization file with the solution", plot);
    cmd.addReal( "f", "IntensityMAE", "Intensity of density function",  IntensityMAE);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    // Load the file
    gsFileData<> fd(fn);
    gsInfo << "Loaded file " << fd.lastPath() << "\n";
    // Create a gsMultipatch and add the loaded geometry
    gsMultiPatch<> mpLeft; // gsNurbsCreator<>::BSplineSquareGrid(1,1,1, 0.0, 0.0);
    fd.getId(1,mpLeft);
    // Elevate and p-refine the basis to order p + numElevate
    // where p is the highest degree in the bases
    mpLeft.degreeElevate(numElevate);
    mpLeft.computeTopology();

    // .... one single patch
    //gsInfo << "INFO IN PARAMETRIC DOMAIN "<< mpLeft.dim() << mpLeft.parDim() <<"\n";
    gsMultiPatch<> mp = gsNurbsCreator<>::BSplineSquareGrid(1,1,1, 0.0, 0.0);
    //Get all interfaces and boundaries:
    mp.degreeElevate(numElevate);
    mp.computeTopology();
    //mp.addAutoBoundaries();

    //..... Test 2
    // Manufactured identity mapping
    gsFunctionExpr<> sN("x","y",2);
    // Right-hand side function : Analytical density function (det(H(u))=f= sigma/rho)
    // Load the file
    gsFunctionExpr<> f;
    fd.getId(2003, f);
    gsInfo<<"Density function "<< f << "\n";

    gsBoundaryConditions<> bc;
    bc.setGeoMap(mp);
    // For simplicity, set Neumann boundary conditions
   for ( gsMultiPatch<>::const_biterator
            bit = mp.bBegin(); bit != mp.bEnd(); ++bit)
   {
       bc.addCondition( *bit, condition_type::neumann, &sN );
   }
    gsInfo<<"Boundary conditions:\n"<< bc <<"\n";

    //! [Refinement]
    gsMultiBasis<double> dbasis(mpLeft, true);//true: poly-splines (not NURBS)
    
    gsInfo << "Patches: "<< mp.nPatches() <<", degree: "<< dbasis.minCwiseDegree() <<"\n";
#ifdef _OPENMP
    gsInfo<< "Available threads: "<< omp_get_max_threads() <<"\n";
#endif
    //! [Refinement]
    
    //! [Problem setup]
    gsExprAssembler<> A(1,1);
    //A.setOptions(Aopt);
    // A.options().addInt("quRule",
    //              "Quadrature rule [1:GaussLegendre,2:GaussLobatto,3:PatchRule]",
    //              1);
    A.options().setReal("quA", 2.0);
    A.options().setInt("quB", 2);
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
    auto ITdim     = mpLeft.geoDim();

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
    //auto Poisson  = gsPatchPreconditionersCreator<>::fastDiagonalizationOp(dbasis.basis(0),bc,A.options(), 1.,eps,0.);
    gsPatchPreconditionersCreator<double>::Poisson_FastDiag Poisson(dbasis.basis(0), bc, A.options(), eps);
    slv_time += timer.stop();

    u.setup(bc, dirichlet::l2Projection, 0);
    // Compute the system matrix and right-hand side

    // Initialize the system :  identity mapping as initial guess
    A.initSystem();
    // Initialize the system : start Computing the conductivity coeffeicient ...
    // Compute the Neumann terms defined on physical space
    gsInfo << "evaluate integral " << ev.integral(ff.val()) << "\n";

    auto g_N = A.getBdrFunction(G);
    auto Neumann_Int{ev.integralBdrBc(bc.get("Neumann"), g_N.tr() * nv(G) )};
    //... nromalisation of density function
    auto CoeffDensity{ev.integral((1.+IntensityMAE*ff.val()))};
    auto CoeffConductivity{Neumann_Int/ev.integral(pow(pow(IGdim,IGdim)-gammaMAE+gammaMAE * CoeffDensity/(1.+IntensityMAE*ff.val()), 1./IGdim) )};

    setup_time += timer.stop();

    gsInfo<< A.numDofs() <<std::flush;

    timer.restart();
    A.assemble(
    igrad(u, G) * igrad(u, G).tr()  + eps * u *u.tr() //matrix
    ,
    u*  CoeffConductivity * (-1.)*pow(pow(IGdim,IGdim)-gammaMAE+gammaMAE* CoeffDensity/(1.+IntensityMAE*ff.val()), 1./IGdim)  //rhs vector
    );
    
    // Compute the Neumann terms defined on physical space
    //auto g_N = A.getBdrFunction(G);
    A.assembleBdr(bc.get("Neumann"), u * g_N.tr() * nv(G) );
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

    gsInfo<< "." << solVector.size() <<std::flush; // Linear solving done
 
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
        
        gsMultiPatch<> Psi, PsiLoc;
        v_sol.extract(Psi);
        v_sol.extract(PsiLoc);

        // ... correct boundary
        ProjectionNormalCPoints(Psi);
        ProjectionNormalCPoints(PsiLoc);
        PsiLoc.addAutoBoundaries();
        PsiLoc.computeTopology();
        Psi.addAutoBoundaries();
        Psi.computeTopology();
        geometryMap PP    = A.getMap(Psi);
        geometryMap PPLoc = A.getMap(PsiLoc);
        //::::::::::::::::::::    Compute the composition of geometry maps      :::::::::::::::::::::::::
        PPLoc(mpLeft);
        auto comp0 = A.getCoeff(mpLeft, PP);

        gsInfo << std::setprecision(14) <<"\n\n";
        auto int_ex = ev.integral( (GLeft).sqNorm()/meas(PP)); /// Exact integral
        gsInfo << "integral getCoef ERR "<< int_ex - ev.integral( (PP).sqNorm()) <<"\n";
        gsInfo << "integral getCoef ERR "<< int_ex - ev.integral( (comp0).sqNorm()) <<"\n";//Strange behavior
        gsInfo << "integral Comp ERR "<< int_ex - ev.integral( (PPLoc).sqNorm()) <<"\n";//Strange behavior
        gsInfo <<"max error Quadrature " << ev.max( (comp0-PPLoc).norm() ) <<"\n";// Strange behavior they are the same

        A.initSystem(ITdim);
        //Obtain control points for the gradient of mpLeft.comp(Psi)
        A.assemble( v * v.tr() , v * PPLoc.tr() );// blocked by this one
        vsolVector = Poisson.L2ProjectVec(A.rhs());
        //vsolVector = solver.compute(A.matrix()).solve(A.rhs());
        v_sol.extract(PsiLoc);// Psiloc is reparameterization of Gleft o Psi by L2projection
        //::::::::::::::::::::      end       ::::::::::::::::::::::::: 
        geometryMap PPfLoc = A.getMap(PsiLoc);
        auto ff     = A.getCoeff(f, PPfLoc);/// Gleft o Psi
        auto ffG     = A.getCoeff(f, GLeft, PP);/// Exact composition


        // gsMatrix<> ptst(2, 1); // Create a 2x1 matrix (2D point)
        // ptst(0, 0) = 0.385;     // Set the x-coordinate
        // ptst(1, 0) = 1.;     // Set the y-coordinate
        // auto CCff  =  f.eval(mpLeft.piece(0).eval( Psi.piece(0).eval(ptst) ));
        // gsInfo <<"foFoPsi exact : "<< CCff <<"\n";
        // gsInfo <<" pt3-Stcompo : "<< ev.eval(ffG, ptst) <<"\n";
        // gsInfo <<" pt3-getCoeff : "<< ev.eval(ff, ptst) <<"\n";
        // gsInfo << "Error_" << ev.integral((PP-grad(u_sol).tr()).sqNorm() ) << "...\n";

        // auto ff = A.getCoeff(f, GLeft, PP);
        // auto ffG = A.getCoeff(f, GLeft);
        gsInfo << "ff" << ev.integral(ff.val()) << "...\n";
        gsInfo << "ffG " << ev.integral(ffG.val()) << "...\n";
        gsInfo << "ffDiff " << ev.integral(ff.val()-ffG.val()) << "...\n";

        // ...  0  dirichlet for boundaries
        sv0 = solVector;
        u.setup(bc, dirichlet::l2Projection, 0);
    
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
        // MAE system
        gsInfo << " Err Mae RHS :  "<< IntegDensity -ev.integral(pow( pow(div(PP).val(),IGdim) + gammaMAE*(CoeffDensity/(1.+IntensityMAE*ffG.val()) - jac(PP).det()), 1./IGdim))<< "\n";
        
        A.assemble(
        grad(u) * grad(u).tr()  +  eps * u * u.tr()//matrix
        ,
        u * CoeffConductivity * (-1.) * ExprMAE  //rhs vector
        );
        //gsInfo << "End Assemnles \n";

        // Compute the Neumann terms defined on physical space
        auto g_N = A.getBdrFunction(G);
        A.assembleBdr(bc.get("Neumann"), u * g_N.tr() * nv(G) );
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
        gsInfo<< "\n Niter in Picard : " << ip << ".. H1 residual : "<<std::scientific<<l2errRes<< ".. L2 MAE residual : "<<std::scientific<<L2MAERes<<"\n";
        if ( l2errRes < tolPicard || ip == maxIter ){
            // ! end Picard loop
            break; 
            } // 
    }//for loop
        // omp_set_dynamic(0);     // Explicitly disable dynamic teams
        // omp_set_num_threads(1); // Use these threads for later parallel regions
 
    timer.stop();
    gsInfo<<"\n\nTotal time: "<< setup_time+ma_time+slv_time<<"\n";
    gsInfo<<"     Setup: "<< setup_time <<"\n";
    gsInfo<<"  Assembly: "<< ma_time    <<"\n";
    gsInfo<<"   Solving: "<< slv_time   <<"\n";

    //! [Export visualization in ParaView]
    if (plot)
    {
        gsMultiPatch<> UU;
        u_sol.extract(UU);
        gsWrite(UU, "U_solution");
        auto u_s = A.getCoeff(UU);

        //gsMultiBasis<> gbasis(dbasis);
        //gbasis.reduceContinuity(1);
        space v = A.getSpace(dbasis);
        gsMatrix<> vsolVector;
        solution v_sol = A.getSolution(v, vsolVector);
        A.initSystem(IGdim);

        //gsVector<> pt(2); pt.setConstant(0.5);
        //ev.testEval( v, pt );
        //ev.testEval( igrad(u_sol,G), pt );

        // Obtain control points for the gradient of Psi
        A.assemble( v * v.tr() , v * igrad(u_s,G));
        vsolVector = Poisson.L2ProjectVec(A.rhs());
        gsMultiPatch<> Psi, Psitp;
        v_sol.extract(Psitp);
        //... correct the boundary
        ProjectionNormalCPoints(Psitp);

        for (index_t i = 0; i < numrRefine; i++){
            // Psi.addAutoBoundaries();
            geometryMap PP = A.getMap(Psitp);
            //gsInfo<< " Int  "<< ev.integral(PP.sqNorm()) << "\n";

            auto  comp = PP(Psitp);
            //auto comp = A.getCoeff(Psitp, PP);
            A.initSystem(IGdim);
            //Obtain control points for the gradient of mpLeft.comp(Psi)
            A.assemble( v * v.tr() , v * comp.tr() );// blocked by this one
            vsolVector = Poisson.L2ProjectVec(A.rhs());
            v_sol.extract(Psitp);
            Psitp.addAutoBoundaries();
            Psitp.computeTopology();
        }
        //::::::::::::::::::::    Compute the composition of geometry maps      :::::::::::::::::::::::::
        // Psi.addAutoBoundaries();
        geometryMap PP = A.getMap(Psitp);
        gsInfo<< " Int  "<< ev.integral(PP.sqNorm()) << "\n";

        PP(mpLeft);
        //auto comp = A.getCoeff(mpLeft, PP);
        A.initSystem(ITdim);
        //Obtain control points for the gradient of mpLeft.comp(Psi)
        A.assemble( v * v.tr() , v * PP.tr() );// blocked by this one
        vsolVector = Poisson.L2ProjectVec(A.rhs());
        v_sol.extract(Psitp);
        Psitp.addAutoBoundaries();
        Psitp.computeTopology();
        gsInfo << "end of adaptive mapping computation\n" << Psitp<< "\n";

        for(size_t i =0; i<Psitp.nPatches(); ++i)
            Psi.addPatch(gsTHBSpline<2>( dynamic_cast<const gsTensorBSpline<2>&>(Psitp.patch(i)) ));
        Psi.addAutoBoundaries();
        Psi.computeTopology();

        //Psi.uniformRefine();
        gsMultiBasis<> dbasis(Psi, true);//true: poly-splines (not NURBS)

        geometryMap PPF = A.getMap(Psi);
        auto ff_TG      = A.getCoeff(f, PPF);
        // --------------- adaptive refinement ---------------
        // Specify cell-marking strategy...
        MarkingStrategy adaptRefCrit = PUCA;
        //MarkingStrategy adaptRefCrit = GARU;
        //MarkingStrategy adaptRefCrit = errorFraction;
        real_t adaptRefParam = 0.7;
        // Elements used for numerical integration
        A.setIntegrationElements(dbasis);
        gsExprEvaluator<> ev(A);

        for (int r=0; r<=numLRefine; ++r)
        {
        //! [beginRefLoop]
            gsInfo << "====== Loop " << r << " of "
                    <<numLRefine<< " ====adapt Parameter ="<< adaptRefParam << " ======" << "\n";
            // --------------- error estimation/computation ---------------
            // Get the element-wise norms.
            ev.integralElWise( ( ff_TG ).sqNorm() );
            const std::vector<real_t> eltErrs  = ev.elementwise();
            //! [errorComputation]

            //! [adaptRefinementPart]
            // Mark elements for refinement, based on the computed local errors and
            // the refinement-criterion and -parameter.
            std::vector<bool> elMarked( eltErrs.size() );
            gsMarkElementsForRef( eltErrs, adaptRefCrit, adaptRefParam, elMarked);
            gsInfo <<"Marked "<< std::count(elMarked.begin(), elMarked.end(), true) <<" elements.\n";
            // Refine the marked elements with a 1-ring of cells around marked elements
            gsRefineMarkedElements( dbasis, elMarked, 1);
            gsRefineMarkedElements( Psi, elMarked, 1);
            }

        //::::::::::::::::::::      end       :::::::::::::::::::::::::   
        gsInfo<<"Plotting in Paraview...\n";
        gsParaviewCollection collection("ParaviewOutput/solution", &ev);
        collection.options().setSwitch("plotElements", true);
        collection.options().setSwitch("base64", export_b64);
        collection.options().setInt("plotElements.resolution", 16);
        collection.options().setInt("numPoints", 10000);
        collection.newTimeStep(&Psi);
        collection.addField(ff_TG, "density function");
        collection.addField(jac(PPF).det(), "Jacobian function");
        collection.saveTimeStep();
        collection.save();

        gsFileManager::open("ParaviewOutput/solution.pvd");
        // gsWrite(Psi, "Psi_mapping");
        // gsInfo << "Result written in Psi_mapping.xml \n";
    }
    else
        gsInfo << "Done. No output created, re-run with --plot to get a ParaView "
                  "file containing the solution.\n";
    //! [Export visualization in ParaView]

    return EXIT_SUCCESS;


}// end main