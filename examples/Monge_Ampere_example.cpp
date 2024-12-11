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

void ProjectionNormalCPointsAll(gsMultiPatch<>& Psi, gsMultiPatch<> mp){
    // Projection normal of control points (exact geometry)
    int boxMaxNumber = mp.nBoxes();
    for (int boxNumber = 0; boxNumber < boxMaxNumber; ++boxNumber)
    {
        // test if the boundary interface is not an inner interface between patches
        for (int i_x =0; i_x < Psi.patch(boxNumber).basis().boundary(1).size(); ++i_x) // x=0 control points be like (0,:) in this case
        {
            Psi.patch(boxNumber).coef( Psi.patch(boxNumber).basis().boundary(1).at(i_x) ).array()[0] = mp.patch(boxNumber).coef( mp.patch(boxNumber).basis().boundary(1).at(0) ).array()[0];
        }
        for (int i_x =0; i_x < Psi.patch(boxNumber).basis().boundary(2).size(); ++i_x)// x=1 control points be like (1,:) in this case
        {
        Psi.patch(boxNumber).coef( Psi.patch(boxNumber).basis().boundary(2).at(i_x) ).array()[0] = mp.patch(boxNumber).coef( mp.patch(boxNumber).basis().boundary(1).at(0) ).array()[0] + 1.;
        }
        for (int i_x =0; i_x < Psi.patch(boxNumber).basis().boundary(3).size(); ++i_x) // y=0 control points be like (:,0) in this case
        {
        Psi.patch(boxNumber).coef( Psi.patch(boxNumber).basis().boundary(3).at(i_x) ).array()[1] = mp.patch(boxNumber).coef( mp.patch(boxNumber).basis().boundary(1).at(0) ).array()[1];
        }
        for (int i_x =0; i_x < Psi.patch(boxNumber).basis().boundary(4).size(); ++i_x)// y=1 control points be like (:,1) in this case
        {
        Psi.patch(boxNumber).coef( Psi.patch(boxNumber).basis().boundary(4).at(i_x) ).array()[1] = mp.patch(boxNumber).coef( mp.patch(boxNumber).basis().boundary(1).at(0) ).array()[1]+1.;
        }
        }
};

void ProjectionNormalCPoints(gsMultiPatch<>& Psi, gsMultiPatch<> mp){
    // Projection normal of control points (exact geometry)
    int boxMaxNumber = mp.nBoxes();
    for (int boxNumber = 0; boxNumber < boxMaxNumber; ++boxNumber)
    {
        // test if the boundary interface is not an inner interface between patches
        if(!mp.isInterface( patchSide(boxNumber,1) ) ){
        for (int i_x =0; i_x < Psi.patch(boxNumber).basis().boundary(1).size(); ++i_x) // x=0 control points be like (0,:) in this case
        {
            Psi.patch(boxNumber).coef( Psi.patch(boxNumber).basis().boundary(1).at(i_x) ).array()[0] = mp.patch(boxNumber).coef( mp.patch(boxNumber).basis().boundary(1).at(0) ).array()[0];
        }
        }

        if(!mp.isInterface( patchSide(boxNumber,2) ) ){
        for (int i_x =0; i_x < Psi.patch(boxNumber).basis().boundary(2).size(); ++i_x)// x=1 control points be like (1,:) in this case
        {
        Psi.patch(boxNumber).coef( Psi.patch(boxNumber).basis().boundary(2).at(i_x) ).array()[0] = mp.patch(boxNumber).coef( mp.patch(boxNumber).basis().boundary(1).at(0) ).array()[0] + 1.;
        }
        }

        if(!mp.isInterface( patchSide(boxNumber,3) ) ){
        for (int i_x =0; i_x < Psi.patch(boxNumber).basis().boundary(3).size(); ++i_x) // y=0 control points be like (:,0) in this case
        {
        Psi.patch(boxNumber).coef( Psi.patch(boxNumber).basis().boundary(3).at(i_x) ).array()[1] = mp.patch(boxNumber).coef( mp.patch(boxNumber).basis().boundary(1).at(0) ).array()[1];
        }
        }
        if(!mp.isInterface( patchSide(boxNumber,4) ) ){
        for (int i_x =0; i_x < Psi.patch(boxNumber).basis().boundary(4).size(); ++i_x)// y=1 control points be like (:,1) in this case
        {
        Psi.patch(boxNumber).coef( Psi.patch(boxNumber).basis().boundary(4).at(i_x) ).array()[1] = mp.patch(boxNumber).coef( mp.patch(boxNumber).basis().boundary(1).at(0) ).array()[1]+1.;
        }
        }
        }
};

void MatchtangentialCPoints(gsMultiPatch<>& Psi, int boxNumber, int boxNumberR, int s){
    // Projecting tangantial of control points (matching geometry)
    // looking for boundary interface
    if( s==1){
    for (int i_x =0; i_x < Psi.patch(boxNumber).basis().boundary(1).size(); ++i_x) // x=0 control points be like (0,:) in this case
    {
        Psi.patch(boxNumber).coef( Psi.patch(boxNumber).basis().boundary(1).at(i_x) ).array()[1]  += Psi.patch(boxNumberR).coef( Psi.patch(boxNumberR).basis().boundary(2).at(i_x) ).array()[1];
        Psi.patch(boxNumber).coef( Psi.patch(boxNumber).basis().boundary(1).at(i_x) ).array()[1]  *= 0.5;
        Psi.patch(boxNumberR).coef( Psi.patch(boxNumberR).basis().boundary(2).at(i_x) ).array()[1] = Psi.patch(boxNumber).coef( Psi.patch(boxNumber).basis().boundary(1).at(i_x) ).array()[1];
    }
    }
    if(s==2 ){
    for (int i_x =0; i_x < Psi.patch(boxNumber).basis().boundary(2).size(); ++i_x)// x=1 control points be like (1,:) in this case
    {
    Psi.patch(boxNumber).coef( Psi.patch(boxNumber).basis().boundary(2).at(i_x) ).array()[1]  += Psi.patch(boxNumberR).coef( Psi.patch(boxNumberR).basis().boundary(1).at(i_x) ).array()[1];
    Psi.patch(boxNumber).coef( Psi.patch(boxNumber).basis().boundary(2).at(i_x) ).array()[1]  *= 0.5;
    Psi.patch(boxNumberR).coef( Psi.patch(boxNumberR).basis().boundary(1).at(i_x) ).array()[1] = Psi.patch(boxNumber).coef( Psi.patch(boxNumber).basis().boundary(2).at(i_x) ).array()[1];
    }
    }

    if(s==3){
    for (int i_x =0; i_x < Psi.patch(boxNumber).basis().boundary(3).size(); ++i_x) // y=0 control points be like (:,0) in this case
    {
    Psi.patch(boxNumber).coef( Psi.patch(boxNumber).basis().boundary(3).at(i_x) ).array()[0]  += Psi.patch(boxNumberR).coef( Psi.patch(boxNumberR).basis().boundary(4).at(i_x) ).array()[0];
    Psi.patch(boxNumber).coef( Psi.patch(boxNumber).basis().boundary(3).at(i_x) ).array()[0]  *= 0.5;
    Psi.patch(boxNumberR).coef( Psi.patch(boxNumberR).basis().boundary(4).at(i_x) ).array()[0] = Psi.patch(boxNumber).coef( Psi.patch(boxNumber).basis().boundary(3).at(i_x) ).array()[0];
    }
    }
    if(s==4){
    for (int i_x =0; i_x < Psi.patch(boxNumber).basis().boundary(4).size(); ++i_x)// y=1 control points be like (:,1) in this case
    {
        Psi.patch(boxNumber).coef( Psi.patch(boxNumber).basis().boundary(4).at(i_x) ).array()[0]  += Psi.patch(boxNumberR).coef( Psi.patch(boxNumberR).basis().boundary(3).at(i_x) ).array()[0];
        Psi.patch(boxNumber).coef( Psi.patch(boxNumber).basis().boundary(4).at(i_x) ).array()[0]  *= 0.5;
        Psi.patch(boxNumberR).coef( Psi.patch(boxNumberR).basis().boundary(3).at(i_x) ).array()[0] = Psi.patch(boxNumber).coef( Psi.patch(boxNumber).basis().boundary(4).at(i_x) ).array()[0];
    }
    }
};


int main(int argc, char *argv[])
{
    //! [Parse command line]
    bool plot          = false;
    index_t numRefine  = 3;
    index_t numElevate = 0;
    index_t maxIter    = 30;
    double eps         = 1e-5; // pinalization coefficient
    double tolPicard   = 1e-8;
    bool last = false, export_b64{false}, adaptiveMesh{true};
    // ...PNormalCP: Correct the normal part of the mapping and CornersLshape: adjust the corners of the three patches that form L.
    bool PNormalCP{false};
    if(adaptiveMesh){
        PNormalCP = true;
    }

    gsCmdLine cmd("Tutorial on solving a non-linear Monge-Ampere problem.");
    cmd.addInt("i", "iter", "Maximum number of iterations for the iterative Picard", maxIter);
    cmd.addInt( "e", "degreeElevation",
                "Number of degree elevation steps to perform before solving (0: equalize degree in all directions)", numElevate );
    cmd.addInt( "r", "uniformRefine", "Number of Uniform h-refinement loops",  numRefine );
    //cmd.addString( "f", "file", "Input XML file", fn );
    cmd.addInt("quRule",
                 "Quadrature rule [1:GaussLegendre,2:GaussLobatto,3:PatchRule]",
                 2);
    cmd.addSwitch("last", "Solve solely for the last level of h-refinement", last);
    cmd.addSwitch("plot", "Create a ParaView visualization file with the solution", plot);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    //gsFileData<> fd(fn);
    //gsInfo << "Loaded file "<< fd.lastPath() <<"\n";
    // .... one single patch
   gsMultiPatch<> mp = gsNurbsCreator<>::BSplineSquareGrid(1,1,1, 0.0, 0.0);
   //... patch 1
   //mp.addPatch(gsNurbsCreator<>::BSplineSquare(1, 0.,1.0));
//    //... patch 2 (L-shape)
    // mp.addPatch(gsNurbsCreator<>::BSplineSquare(1, 1.,0.0));
   //mp.addInterface(0,2,2,1);
   
//    // ... patch 0-1
//    gsMultiPatch<> mp = gsNurbsCreator<>::BSplineSquareGrid(1,2,1, -1.0, -1.0);
//    // ... patch 2
//    mp.addPatch(gsNurbsCreator<>::BSplineSquare(1,-1.,-2.0));
//    mp.addInterface(0,3,2,4);
//    // ... patch 3
//    mp.addPatch(gsNurbsCreator<>::BSplineSquare(1,0.,-2.0));
//    mp.addInterface(2,2,3,1);
//    // ... patch 4
//    mp.addPatch(gsNurbsCreator<>::BSplineSquare(1, 1.,-2.0));
//    mp.addInterface(3,2,4,1);
//    // ... patch 5
//    mp.addPatch(gsNurbsCreator<>::BSplineSquare(1, 1.,-1.0));
//    mp.addInterface(4,4,5,3);
//    // ... patch 6
//    mp.addPatch(gsNurbsCreator<>::BSplineSquare(1, 1., 0.0));
//    mp.addInterface(5,4,6,3);
//    // ... patch 7
//    mp.addPatch(gsNurbsCreator<>::BSplineSquare(1, 0.0, 0.0));
//    mp.addInterface(6,1,7,2);
//    mp.addInterface(1,2,7,1);
   // Get all interfaces and boundaries:
    mp.computeTopology();
    //mp.addAutoBoundaries();

    //..... Test 2
    // // Manufactured solition
    // gsFunctionExpr<> s("exp(0.5*(x**2 + y**2))",2);
    // // Manufactured Grad solition
    // gsFunctionExpr<> sN("x*exp(0.5*(x**2 + y**2))","y*exp(0.5*(x**2 + y**2))",2);
    // // Right-hand side function
    // gsFunctionExpr<> f("(1.+x**2+y**2)*exp(x**2 + y**2)",2);

    //..... Test 2
    //convex function
    gsFunctionExpr<> s("0.5*(x**2 + y**2)",2);
    // Manufactured identity mapping
    gsFunctionExpr<> sN("x","y",2);
    // Right-hand side function : Analytical density function (det(H(u))=f= sigma/rho)
    //gsFunctionExpr<> f("1./(2.+cos(8.*pi*sqrt((x-0.5-0.25*0.)**2+(y-0.5)**2)))",2);
    //gsFunctionExpr<> f("1.+6.*( 1/(1.+exp((y -x  - 0.3)/0.01)) - 1/(1.+exp((y - x  - 0.1)/0.01)) )",2);
    // Manufactured density function
    gsFunctionExpr<> f("1.+6.*( 1/(1.+exp(( (y-0.98)**2+(x-0.899)**2  - 0.002)/0.001)) + 1/(1.+exp((y -x  - 0.25)/0.001)) - 1./(1.+exp((y - x  - 0.15)/0.001)) +  0/(1.+exp((y - 1.0)/0.001)) - 0./(1.+exp((y - 0.975)/0.001))  +  1/(1.+exp((x - 1.0)/0.001)) - 1./(1.+exp((x - 0.95)/0.001)) )",2);    
    //gsFunctionExpr<> f("(1.+ 9./(1.+(10.*sqrt((x-0.7-0.25*0.)**2+(y-0.5)**2)*cos(atan2(y-0.5,x-0.7-0.25*0.) -20.*((x-0.7-0.25*0.)**2+(y-0.5)**2)))**2) )",2);
    //gsFunctionExpr<> f("( 1.+ 5.*exp(-50.*abs((x-0.5-0.25*cos(2.*pi*0.25))**2-(y-0.5-0.5 *sin(2.*pi*0.25))**2- 0.01)))",2);
    //gsFunctionExpr<> f("(1. + 5./cosh( 5.*((x-sqrt(3)/2)**2+(y-0.5)**2 - (pi/2)**2) )**2 + 5./cosh( 5.*((x+sqrt(3)/2)**2+(y-0.5)**2 - (pi/2)**2) )**2)",2);
    gsInfo<<"Source function "<< f << "\n";

    gsInfo<<"The domain is "<< mp.detail() << "\n";

    gsBoundaryConditions<> bc;
    bc.setGeoMap(mp);
    // For simplicity, set Neumann boundary conditions
    //bc.addCornerValue(4,1.,0,0);
    //bc.addCornerValue(2,1.,1,0);
    //bc.addCoupled(0,4,1,3,2,0);
   for ( gsMultiPatch<>::const_biterator
            bit = mp.bBegin(); bit != mp.bEnd(); ++bit)
   {
       bc.addCondition( *bit, condition_type::neumann, &sN );
   }
    gsInfo<<"Boundary conditions:\n"<< bc <<"\n";

    //! [Refinement]
    gsMultiBasis<> dbasis(mp, true);//true: poly-splines (not NURBS)

    // Elevate and p-refine the basis to order p + numElevate
    // where p is the highest degree in the bases
    dbasis.setDegree( dbasis.maxCwiseDegree() + numElevate);
    dbasis.degreeElevate(2);

    gsInfo << "Patches: "<< mp.nPatches() <<", degree: "<< dbasis.minCwiseDegree() <<"\n";
#ifdef _OPENMP
    gsInfo<< "Available threads: "<< omp_get_max_threads() <<"\n";
#endif
    //! [Refinement]

    //! [Problem setup]
    gsExprAssembler<> A(1,1);
    //A.setOptions(Aopt);
    A.options().addInt("quRule",
                 "Quadrature rule [1:GaussLegendre,2:GaussLobatto,3:PatchRule]",
                 1);
    A.options().addInt("InterfaceStrategy", "Interface strategy conforming", iFace::none);
    gsInfo<<"Active options:\n"<< A.options() <<"\n";

    typedef gsExprAssembler<>::geometryMap geometryMap;
    typedef gsExprAssembler<>::variable    variable;
    typedef gsExprAssembler<>::space       space;
    typedef gsExprAssembler<>::solution    solution;

    // Elements used for numerical integration
    A.setIntegrationElements(dbasis);
    gsExprEvaluator<> ev(A);
    ev.options().addInt("quRule",
                 "Quadrature rule [1:GaussLegendre,2:GaussLobatto,3:PatchRule]",
                1);
    ev.options().addInt("InterfaceStrategy", "Interface strategy conforming", iFace::none);
    // Set the geometry map
    geometryMap G = A.getMap(mp);

    // Set pow for BFO method
    auto IGdim     = 1./G.domainDim();
    // Set factor for BFO method
    auto gammaMAE = factorial(G.domainDim());

    // Set the discretization space
    space u = A.getSpace(dbasis);
    
    // Set the source term
    auto ff = A.getCoeff(f, G);

    // Recover manufactured solution
    auto u_ex = ev.getVariable(s, G);

    //gsFunctionExpr<> sI("0.5*(x**2+y**2)+x*y",2);
    auto u_I = ev.getVariable(sN, G);

    // Solution vector and solution variable
    gsMatrix<> solVector;
    solution u_sol = A.getSolution(u, solVector);
    //! [Problem setup]

    //! [Solver loop]
    gsSparseSolver<>::CGDiagonal solver;

    gsVector<>  h1err(numRefine+1); //l2err(numRefine+1) : The solution exists up to an additive constant.
    gsInfo<< "(dot1=assembled, dot2=solved, dot3=nonlinear_loop,dot4=got_error)\n"
        "\nDoFs: ";
    double setup_time(0), ma_time(0), slv_time(0), err_time(0);    
    gsStopwatch timer;
    // Initialize the system : start Computing the conductivity coeffeicient ...
    // Compute the Neumann terms defined on physical space
    auto g_N = A.getBdrFunction(G);
    auto Neumann_Int{ev.integralBdrBc(bc.get("Neumann"), g_N.tr() * nv(G) )};
    //... nromalisation of density function
    auto CoeffDensity{ev.integral(ff.val() * meas(G))};
    auto CoeffConductivity{Neumann_Int/ev.integral(pow(gammaMAE * CoeffDensity/ff.val(), IGdim) * meas(G))};
    if(adaptiveMesh)
    {
        //::::::::::::::::::::      mesh adaptation solver         :::::::::::::::::::::::::
        for (int r=0; r<=numRefine; ++r)
        {
            dbasis.uniformRefine();
            mp.uniformRefine();
            if( last && r != numRefine)
                continue;
            u.setup(bc, dirichlet::l2Projection, 0);
            // Compute the system matrix and right-hand side
            // ...
            //... end 

            // Initialize the system :  identity mapping as initial guess
            A.initSystem();
            setup_time += timer.stop();

            gsInfo<< A.numDofs() <<std::flush;
//Optionally, assemble Nitsche
        gsMatrix<> mu_interfaces;
        //auto m_penalty = 1;
        index_t i = 0;
        for ( typename gsMultiPatch<>::const_iiterator it = mp.iBegin(); it != mp.iEnd(); ++it, ++i)
        {
            auto stab     = 4 * ( dbasis.maxCwiseDegree() + dbasis.dim() ) * ( dbasis.maxCwiseDegree() + 1 );
           auto m_h      = dbasis.basis(0).getMinCellLength(); // m_basis.basis(0).getMinCellLength();
           //auto mu       = 2 * stab / m_h;
            auto alpha   = 2 * stab / m_h;

            // //mu = penalty_init == -1.0 ? mu : penalty_init / m_h;
            // if (m_penalty == -1)
            //     mu = mu_interfaces(i,0) / m_h;
            // else
            //     mu = m_penalty / m_h;

            std::vector<boundaryInterface> iFace;
            iFace.push_back(*it);
            A.assembleIfc(iFace,
                    // //B11
                        //   +mu * u.left() *
                        //   u.left().tr() * nv(G.left()).norm(),
                        //  -0.5*alpha *
                        //   (igrad(u.left(), G) * tv(G.left()).normalized() * u.left().tr()).tr() *
                        //   nv(G.left()).norm(),
                        //  -0.5*alpha *
                        //     (u.left()*(igrad(u.left(), G) * tv(G.left()).normalized()).tr()).tr() *
                        //   nv(G.left()).norm(),
                    // //B12
                        //  -mu * u.left() *
                        //   u.right().tr() * nv(G.left()).norm(),
                        //  +0.5*alpha *
                        //   (igrad(u.left(), G) * tv(G.left()).normalized() * u.right().tr()).tr() *
                        //   nv(G.left()).norm(),
                        //  -0.5*alpha *
                        //     (u.left()*(igrad(u.right(), G) * tv(G.right()).normalized()).tr()).tr() *
                        //   nv(G.left()).norm(),
                    // //B21
                        //  - mu * u.right() *
                        //   u.left().tr() * nv(G.left()).norm(),
                        //  -0.5*alpha *
                        //   (igrad(u.right(), G) * tv(G.right()).normalized() * u.left().tr()).tr() *
                        //   nv(G.left()).norm(),
                        //  +0.5*alpha *
                        //     (u.right()*(igrad(u.left(), G) * tv(G.left()).normalized()).tr()).tr() *
                        //   nv(G.left()).norm(),
                    // //B22
                        //   + mu * u.right() *
                        //   u.right().tr() * nv(G.left()).norm(),
                        //  +0.5*alpha *
                        //   (igrad(u.right(), G) * tv(G.right()).normalized() * u.right().tr()).tr() *
                        //   nv(G.left()).norm(),
                        //  +0.5*alpha *
                        //     (u.right()*(igrad(u.right(), G) * tv(G.right()).normalized()).tr()).tr() *
                        //   nv(G.left()).norm()
                          
                    // // E11
                    //      0.5*alpha * igrad(u.left(), G.left()) * nv(G.left()).normalized() *
                    //       (igrad(u.left(), G.left()) * nv(G.left()).normalized()).tr() * nv(G.left()).norm(),
                    // //-E12
                    //       -0.5*alpha *(igrad(u.left(), G.left()) * nv(G.left()).normalized()) *
                    //       (igrad(u.right(), G.right()) * nv(G.left()).normalized()).tr() * nv(G.left()).norm(),
                    // //-E21
                    //       -0.5*alpha *(igrad(u.right(), G.right()) * nv(G.left()).normalized()) *
                    //       (igrad(u.left(), G.left()) * nv(G.left()).normalized()).tr() * nv(G.left()).norm(),
                    // // E22
                    //      0.5*alpha * igrad(u.right(), G.right()) * nv(G.left()).normalized() *
                    //       (igrad(u.right(), G.right()) * nv(G.left()).normalized()).tr() * nv(G.left()).norm()

                    // //E11
                    //     0.5* alpha * igrad(u.left(), G.left()) * tv(G.left()).normalized() *
                    //       (igrad(u.left(), G.left()) * tv(G.left()).normalized()).tr() * nv(G.left()).norm(),
                    // //-E12
                    //       -0.5* alpha *(igrad(u.left(), G.left()) * tv(G.left()).normalized()) *
                    //       (igrad(u.right(), G.right()) * tv(G.right()).normalized()).tr() * nv(G.left()).norm(),
                    // //-E21
                    //       - 0.5*alpha *(igrad(u.right(), G.right()) * tv(G.right()).normalized()) *
                    //       (igrad(u.left(), G.left()) * tv(G.left()).normalized()).tr() * nv(G.left()).norm(),
                    // //E22
                    //       0.5*alpha * igrad(u.right(), G.right()) * tv(G.right()).normalized() *
                    //       (igrad(u.right(), G.right()) * tv(G.right()).normalized()).tr() * nv(G.left()).norm()
                    
                     //E11  Nitsche or Penality method for forcing tangential part
                        - u.left() *
                          (igrad(u.left(), G.left()) * tv(G.left()).normalized()).tr() * meas(G),
                    //-E12
                          u.left() *
                          (igrad(u.right(), G.right()) * tv(G.left()).normalized()).tr() * meas(G),
                    //-E21
                        -  u.right() *
                          (igrad(u.left(), G.left()) * tv(G.left()).normalized()).tr() * meas(G),
                    //E22
                          u.right() *
                          (igrad(u.right(), G.right()) * tv(G.left()).normalized()).tr() * meas(G),
                    //E11
                          alpha * u.left() *
                          (igrad(u.left(), G.left()) * tv(G.left()).normalized()).tr() * nv(G.left()).norm(),
                    //-E12
                          -alpha *u.left() *
                          (igrad(u.right(), G.right()) * tv(G.left()).normalized()).tr() * nv(G.left()).norm(),
                    //-E21
                          alpha * u.right() *
                          (igrad(u.left(), G.left()) * tv(G.left()).normalized()).tr() * nv(G.left()).norm(),
                    //E22
                         - alpha *u.right() *
                          (igrad(u.right(), G.right()) * tv(G.left()).normalized()).tr() * nv(G.left()).norm()
            );
        }
            timer.restart();
            A.assemble(
            igrad(u, G) * igrad(u, G).tr() * meas(G) + eps * u *u.tr()* meas(G) //matrix
            ,
            u*  CoeffConductivity * (-1.)*pow(gammaMAE* CoeffDensity/ff.val(), IGdim) * meas(G) //rhs vector
            );
            
            // Compute the Neumann terms defined on physical space
            //auto g_N = A.getBdrFunction(G);
            A.assembleBdr(bc.get("Neumann"), u * g_N.tr() * nv(G) );
            A.assembleIfc(mp.interfaces(), u.left() * (u_I.tr() * nv(G.left())));
            A.assembleIfc(mp.interfaces(), u.right() * (u_I.tr() * nv(G.right())));

            ma_time += timer.stop();

            // gsDebugVar(A.matrix().toDense());
            // gsDebugVar(A.rhs().transpose()   );

            gsInfo<< "." <<std::flush;// Assemblying done

            timer.restart();
            solver.compute( A.matrix() );
            solVector = solver.solve(A.rhs());

            // gsMultiPatch<> UUs;
            // u_sol.extract(UUs);
            // gsInfo << "dim of solution---" << UUs.patch(0).coefsSize() << " ---\n";
            // gsInfo << "nBases---size\n" << UUs.patch(0).basis().boundary(4).size() << " ---\n";
            // gsInfo << "nBases---1\n" << UUs.patch(0).basis().boundary(1) << " ---\n";
            // gsInfo << "nBases---2\n" << UUs.patch(0).basis().boundary(2) << " ---\n";
            // gsInfo << "nBases---3\n" << UUs.patch(0).basis().boundary(3) << " ---\n";
            // gsInfo << "nBases---4\n" << UUs.patch(0).basis().boundary(4) << " ---\n";
            // gsInfo << "Im of solution---" << UUs.patch(1).coef( UUs.patch(1).basis().boundary(4).at(UUs.patch(1).basis().boundary(4).size()-1) ).array() << " ---\n";

            slv_time += timer.stop();

            gsInfo<< "." <<std::flush; // Linear solving done

            // Picard loop
            index_t NiterPicard{0};
            gsMatrix<> sv0; //
            solution u_lsol = A.getSolution(u, sv0);
            for(int ip{0}; ip<=maxIter; ++ip)
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
                A.initSystem(2);

                //gsVector<> pt(2); pt.setConstant(0.5);
                //ev.testEval( v, pt );
                //ev.testEval( igrad(u_sol,G), pt );

                // Obtain control points for the gradient of Psi
                A.assemble( v * v.tr() , v * igrad(u_s,G) );
                vsolVector = solver.compute(A.matrix()).solve(A.rhs());
                
                gsMultiPatch<> Psi;
                v_sol.extract(Psi);

                // ... correct boundary
                //if (PNormalCP)
                //    ProjectionNormalCPoints(Psi, mp);
                //if (CornersLshape)
                //    CorrectCornersLshape(Psi, mp);
                
                geometryMap PP = A.getMap(Psi);
                auto ff = A.getCoeff(f,PP);

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
                CoeffConductivity = Neumann_Int/ev.integral(pow( (ilapl(u_sol,G)*ilapl(u_sol,G).tr()).val() + gammaMAE*(CoeffDensity/ff.val() - ihess(u_sol,G).det()), IGdim) * meas(G));
//Optionally, assemble Nitsche
        gsMatrix<> mu_interfaces;
        //auto m_penalty = 1;
        index_t i = 0;
        for ( typename gsMultiPatch<>::const_iiterator it = mp.iBegin(); it != mp.iEnd(); ++it, ++i)
        {
            //auto stab     = 4 * ( dbasis.maxCwiseDegree() + dbasis.dim() ) * ( dbasis.maxCwiseDegree() + 1 );
           //auto m_h      = dbasis.basis(0).getMinCellLength(); // m_basis.basis(0).getMinCellLength();
           // auto mu       = 10.;//2 * stab / m_h;
            auto alpha = 1e+4;

            // //mu = penalty_init == -1.0 ? mu : penalty_init / m_h;
            // if (m_penalty == -1)
            //     mu = mu_interfaces(i,0) / m_h;
            // else
            //     mu = m_penalty / m_h;

            std::vector<boundaryInterface> iFace;
            iFace.push_back(*it);
            solution u_lsol = A.getSolution(u, sv0);
            A.assembleIfc(iFace,
                    // //B11
                    //       +mu * u.left() *
                    //       u.left().tr() * nv(G.left()).norm(),
                    //      -0.5*alpha *
                    //       (igrad(u.left(), G) * nv(G.left()).normalized() * u.left().tr()).tr() *
                    //       nv(G.left()).norm(),
                    //      -0.5*alpha *
                    //         (u.left()*(igrad(u.left(), G) * nv(G.left()).normalized()).tr()).tr() *
                    //       nv(G.left()).norm(),
                    // //B12
                    //      -mu * u.left() *
                    //       u.right().tr() * nv(G.left()).norm(),
                    //      +0.5*alpha *
                    //       (igrad(u.left(), G) * nv(G.left()).normalized() * u.right().tr()).tr() *
                    //       nv(G.left()).norm(),
                    //      -0.5*alpha *
                    //         (u.left()*(igrad(u.right(), G) * nv(G.right()).normalized()).tr()).tr() *
                    //       nv(G.left()).norm(),
                    // //B21
                    //      - mu * u.right() *
                    //       u.left().tr() * nv(G.left()).norm(),
                    //      -0.5*alpha *
                    //       (igrad(u.right(), G) * nv(G.right()).normalized() * u.left().tr()).tr() *
                    //       nv(G.left()).norm(),
                    //      +0.5*alpha *
                    //         (u.right()*(igrad(u.left(), G) * nv(G.left()).normalized()).tr()).tr() *
                    //       nv(G.left()).norm(),
                    // //B22
                    //       + mu * u.right() *
                    //       u.right().tr() * nv(G.left()).norm(),
                    //      +0.5*alpha *
                    //       (igrad(u.right(), G) * nv(G.right()).normalized() * u.right().tr()).tr() *
                    //       nv(G.left()).norm(),
                    //      +0.5*alpha *
                    //         (u.right()*(igrad(u.right(), G) * nv(G.right()).normalized()).tr()).tr() *
                    //       nv(G.left()).norm()
                          
                    // // E11
                    //      0.5*alpha * igrad(u.left(), G.left()) * nv(G.left()).normalized() *
                    //       (igrad(u.left(), G.left()) * nv(G.left()).normalized()).tr() * nv(G.left()).norm(),
                    // //-E12
                    //       -0.5*alpha *(igrad(u.left(), G.left()) * nv(G.left()).normalized()) *
                    //       (igrad(u.right(), G.right()) * nv(G.left()).normalized()).tr() * nv(G.left()).norm(),
                    // //-E21
                    //       -0.5*alpha *(igrad(u.right(), G.right()) * nv(G.left()).normalized()) *
                    //       (igrad(u.left(), G.left()) * nv(G.left()).normalized()).tr() * nv(G.left()).norm(),
                    // // E22
                    //      0.5*alpha * igrad(u.right(), G.right()) * nv(G.left()).normalized() *
                    //       (igrad(u.right(), G.right()) * nv(G.left()).normalized()).tr() * nv(G.left()).norm()

                    // //E11
                    //     0.5* alpha * igrad(u.left(), G.left()) * tv(G.left()).normalized() *
                    //       (igrad(u.left(), G.left()) * tv(G.left()).normalized()).tr() * nv(G.left()).norm(),
                    // //-E12
                    //       -0.5* alpha *(igrad(u.left(), G.left()) * tv(G.left()).normalized()) *
                    //       (igrad(u.right(), G.right()) * tv(G.right()).normalized()).tr() * nv(G.left()).norm(),
                    // //-E21
                    //       - 0.5*alpha *(igrad(u.right(), G.right()) * tv(G.right()).normalized()) *
                    //       (igrad(u.left(), G.left()) * tv(G.left()).normalized()).tr() * nv(G.left()).norm(),
                    // //E22
                    //       0.5*alpha * igrad(u.right(), G.right()) * tv(G.right()).normalized() *
                    //       (igrad(u.right(), G.right()) * tv(G.right()).normalized()).tr() * nv(G.left()).norm()
                    
                     //E11  Nitsche or Penality method for forcing tangential part
                        - u.left() *
                          (igrad(u.left(), G.left()) * tv(G.left()).normalized()).tr() * nv(G.left()).norm(),
                    //-E12
                          u.left() *
                          (igrad(u.right(), G.right()) * tv(G.left()).normalized()).tr() * nv(G.left()).norm(),
                    //-E21
                        -  u.right() *
                          (igrad(u.left(), G.left()) * tv(G.left()).normalized()).tr() * nv(G.left()).norm(),
                    //E22
                          u.right() *
                          (igrad(u.right(), G.right()) * tv(G.left()).normalized()).tr() * nv(G.left()).norm(),
                    //E11
                          alpha * u.left() *
                          (igrad(u.left(), G.left()) * tv(G.left()).normalized()).tr() * nv(G.left()).norm(),
                    //-E12
                          -alpha *u.left() *
                          (igrad(u.right(), G.right()) * tv(G.left()).normalized()).tr() * nv(G.left()).norm(),
                    //-E21
                          -alpha * u.right() *
                          (igrad(u.left(), G.left()) * tv(G.left()).normalized()).tr() * nv(G.left()).norm(),
                    //E22
                          alpha *u.right() *
                          (igrad(u.right(), G.right()) * tv(G.left()).normalized()).tr() * nv(G.left()).norm()

                    // // //
                    //      +0.5*alpha *
                    //       (igrad(u.left(), G) * tv(G.left()).normalized()) * u.left().tr() *
                    //       nv(G.left()).norm(),
                    //     //
                    //      -0.5*alpha *
                    //       (igrad(u.left(), G) * tv(G.left()).normalized()) * u.right().tr() *
                    //        nv(G.left()).norm(),
                    //     //
                    //      -0.5*alpha *
                    //       (igrad(u.right(), G) * tv(G.left()).normalized()) * u.left().tr() *
                    //       nv(G.left()).norm(),
                    //     //
                    //      +0.5*alpha *
                    //       (igrad(u.right(), G) * tv(G.left()).normalized() )* u.right().tr() *
                    //       nv(G.left()).norm()
            );
        }
                // MAE system
                A.assemble(
                igrad(u, G) * igrad(u, G).tr() * meas(G) +  eps * u * u.tr()* meas(G)//matrix
                ,
                u * CoeffConductivity * (-1.) * pow( (ilapl(u_sol,G)*ilapl(u_sol,G).tr()).val() + gammaMAE*(CoeffDensity/ff.val() - ihess(u_sol,G).det()), IGdim) * meas(G) //rhs vector
                );

                // Compute the Neumann terms defined on physical space
                auto g_N = A.getBdrFunction(G);
                A.assembleBdr(bc.get("Neumann"), u * g_N.tr() * nv(G) );
                A.assembleIfc(mp.interfaces(), u.left() * (u_I.tr() * nv(G.left())));
                A.assembleIfc(mp.interfaces(), u.right() * (u_I.tr() * nv(G.right())));



                ma_time += timer.stop();

                // gsDebugVar(A.matrix().toDense());
                // gsDebugVar(A.rhs().transpose()   );
                

                gsInfo<< " ." <<std::flush;// Assemblying done

                timer.restart();
                solver.compute( A.matrix() );
                solVector = solver.solve(A.rhs());
                
                slv_time += timer.stop();

                gsInfo<< "." <<std::flush; // Linear solving done

                // omp_set_dynamic(0);     // Explicitly disable dynamic teams
                // omp_set_num_threads(1); // Use these threads for later parallel regions

                ++NiterPicard;
                auto l2errRes = math::sqrt(ev.integral( ( igrad(u_lsol,G) - igrad(u_sol,G) ).sqNorm() * meas(G) ));
                if ( l2errRes < tolPicard || ip == maxIter ){
                    // ! end Picard loop
                    gsInfo<< "\n Niter in Picard : " << NiterPicard << ".. L2 residual : "<<std::scientific<<l2errRes<<"\n";
                    break; 
                    } // 
            }//for loop
            // omp_set_dynamic(0);     // Explicitly disable dynamic teams
            // omp_set_num_threads(1); // Use these threads for later parallel regions

            timer.restart();
            //l2err[r]= math::sqrt( ev.integral( (u_ex - u_sol).sqNorm() * meas(G) ) );

            h1err[r]= math::sqrt(ev.integral( ( igrad(u_ex) - igrad(u_sol,G) ).sqNorm() * meas(G) ));
            err_time += timer.stop();
            gsInfo<< ". " <<std::flush; // Error computations done

        } //for loop
        //! [Solver loop]

    }
    else
    {
        //::::::::::::::::::::     manufactured exact solution         :::::::::::::::::::::::::
        for (int r=0; r<=numRefine; ++r)
        {
            dbasis.uniformRefine();
            mp.uniformRefine();
            if( last && r != numRefine)
                continue;
            u.setup(bc, dirichlet::l2Projection, 0);
            // Compute the system matrix and right-hand side

            // Initialize the system : start Computing the conductivity coeffeicient ...
            // Compute the Neumann terms defined on physical space
            auto g_N = A.getBdrFunction(G);
            auto Neumann_Int{ev.integralBdrBc(bc.get("Neumann"), g_N.tr() * nv(G) )};
            // ...
            auto CoeffConductivity{Neumann_Int/ev.integral(pow(2.+gammaMAE * ff.val(), IGdim) * meas(G))};
            //... end 

            // Initialize the system
            A.initSystem();
            setup_time += timer.stop();

            gsInfo<< A.numDofs() <<std::flush;

            timer.restart();

            A.assemble(
            igrad(u, G) * igrad(u, G).tr() * meas(G) + eps * u *u.tr()* meas(G) //matrix
            ,
            u*  CoeffConductivity * (-1.)*pow(2.+gammaMAE* ff.val(), IGdim) * meas(G) //rhs vector
            );
            
            // Compute the Neumann terms defined on physical space
            //auto g_N = A.getBdrFunction(G);
            A.assembleBdr(bc.get("Neumann"), u * g_N.tr() * nv(G) );

            ma_time += timer.stop();

            // gsDebugVar(A.matrix().toDense());
            // gsDebugVar(A.rhs().transpose()   );

            gsInfo<< "." <<std::flush;// Assemblying done

            timer.restart();
            solver.compute( A.matrix() );
            solVector = solver.solve(A.rhs());

            slv_time += timer.stop();

            gsInfo<< "." <<std::flush; // Linear solving done

            // Picard loop
            index_t NiterPicard{0};
            gsMatrix<> sv0; //
            solution u_lsol = A.getSolution(u, sv0);
            for(int ip{0}; ip<=maxIter; ++ip)
            {
                sv0 = solVector;
        //        u.setup(bc, dirichlet::interpolation, 0);
                u.setup(bc, dirichlet::l2Projection, 0);
                
                solution u_sol = A.getSolution(u, solVector);

                // Initialize the system
                A.initSystem();
                setup_time += timer.stop();

                //gsInfo<< A.numDofs() <<std::flush;

                timer.restart();
                // Compute the system matrix and right-hand side

                // .. update Coeffeicient of conductivity
                CoeffConductivity = Neumann_Int/ev.integral(pow( (ilapl(u_sol,G)*ilapl(u_sol,G).tr()).val() + gammaMAE*(ff.val() - ihess(u_sol,G).det()), IGdim) * meas(G));

                // MAE system
                A.assemble(
                igrad(u, G) * igrad(u, G).tr() * meas(G) +  eps * u * u.tr()  * meas(G) //matrix
                ,
                u * CoeffConductivity * (-1.) * pow( (ilapl(u_sol,G)*ilapl(u_sol,G).tr()).val() + gammaMAE*(ff.val() - ihess(u_sol,G).det()), IGdim) *meas(G) //rhs vector
                );

                // Compute the Neumann terms defined on physical space
                auto g_N = A.getBdrFunction(G);
                A.assembleBdr(bc.get("Neumann"), u * g_N.tr() * nv(G) );

                ma_time += timer.stop();

                // gsDebugVar(A.matrix().toDense());
                // gsDebugVar(A.rhs().transpose()   );
                

                gsInfo<< " ." <<std::flush;// Assemblying done

                timer.restart();
                solver.compute( A.matrix() );
                solVector = solver.solve(A.rhs());
                slv_time += timer.stop();

                gsInfo<< "." <<std::flush; // Linear solving done

                // omp_set_dynamic(0);     // Explicitly disable dynamic teams
                // omp_set_num_threads(1); // Use these threads for later parallel regions

                ++NiterPicard;
                auto l2errRes = math::sqrt(ev.integral( ( igrad(u_lsol,G) - igrad(u_sol,G) ).sqNorm() * meas(G) ));
                if ( l2errRes < tolPicard  || ip == maxIter ){
                    // ! end Picard loop
                    gsInfo<< "\n Niter in Picard : " << NiterPicard << ".. L2 residual : "<<std::scientific<<l2errRes<<"\n";
                    break; 
                    } // 
            }//for loop
        // omp_set_dynamic(0);     // Explicitly disable dynamic teams
        // omp_set_num_threads(1); // Use these threads for later parallel regions

        timer.restart();
        //l2err[r]= math::sqrt( ev.integral( (u_ex - u_sol).sqNorm() * meas(G) ) );

        h1err[r]= math::sqrt(ev.integral( ( igrad(u_ex) - igrad(u_sol,G) ).sqNorm() * meas(G) ));
        err_time += timer.stop();
        gsInfo<< ". " <<std::flush; // Error computations done
        } //for loop
        //! [Solver loop]

    } // end of solver
    


    timer.stop();
    gsInfo<<"\n\nTotal time: "<< setup_time+ma_time+slv_time+err_time <<"\n";
    gsInfo<<"     Setup: "<< setup_time <<"\n";
    gsInfo<<"  Assembly: "<< ma_time    <<"\n";
    gsInfo<<"   Solving: "<< slv_time   <<"\n";
    gsInfo<<"     Norms: "<< err_time   <<"\n";

    //! [Error and convergence rates]
    //gsInfo<< "\nL2 error: "<<std::scientific<<std::setprecision(3)<<l2err.transpose()<<"\n";
    gsInfo<< "H1 error: "<<std::scientific<<h1err.transpose()<<"\n";

    if (!last && numRefine>0)
    {
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
        collection.options().setInt("numPoints", 100000);
        collection.newTimeStep(&mp);
        collection.addField(u_sol,"numerical solution");
        collection.addField(igrad(u_sol,G),"gradient_numerical solution");
        if(adaptiveMesh){
        collection.addField(ff, "density function");
        collection.addField(ihess(u_sol,G).det(), "Jacobian function");}
        else{
        collection.addField(u_ex, "exact solution");
        }
        if(maxIter == 0)
        collection.addField(CoeffConductivity * (-1.)*pow(2.+gammaMAE * CoeffDensity/ff.val(), IGdim) * meas(G), "MAE_rhs");
        else
        collection.addField(CoeffConductivity * (-1.) * pow( (ilapl(u_sol,G)*ilapl(u_sol,G).tr()).val() + gammaMAE*(CoeffDensity/ff.val() - ihess(u_sol,G).det()), IGdim) * meas(G), "MAE_rhs");
        collection.saveTimeStep();
        collection.save();


        gsFileManager::open("ParaviewOutput/solution.pvd");


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
        gsMultiPatch<> Psi;
        v_sol.extract(Psi);
        
        //... correct the boundary
        if (PNormalCP)
            ProjectionNormalCPoints(Psi, mp);
        //MatchtangentialCPoints(Psi, 0, 1, 4);
        //MatchtangentialCPoints(Psi, 0, 2, 2);
        //geometryMap PP = A.getMap(Psi);
        //auto fp = A.getCoeff(f,PP);

        gsWrite(Psi, "Psi_mapping");
        gsInfo << "Result written in Psi_mapping.xml \n";
    }
    else
        gsInfo << "Done. No output created, re-run with --plot to get a ParaView "
                  "file containing the solution.\n";
    //! [Export visualization in ParaView]

    return EXIT_SUCCESS;


}// end main