/** @file flow-over-heated-plate.cpp

    @brief Heat equation participant for the PreCICE example "flow over heated plate"

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst
*/

#include <gismo.h>

using namespace gismo;

int main(int argc, char *argv[])
{
    //! [Parse command line]
    bool plot = false;
    index_t numRefine  = 5;
    index_t numElevate = 0;
    index_t numRefineT = 0;
    index_t numElevateT= 0;
    index_t steps = 10;
    index_t testCase = 0;
    bool last{false}, export_b64{false};
    real_t alpha = 5e-3;
    real_t beta = 1;

    gsCmdLine cmd("Tutorial on solving a Poisson problem.");
    cmd.addReal( "a", "alpha", "Heat conduction coefficient", alpha );
    cmd.addInt( "e", "degreeElevation",
                "Number of degree elevation steps to perform before solving (0: equalize degree in all directions)", numElevate );
    cmd.addInt( "r", "uniformRefine", "Number of Uniform h-refinement loops",  numRefine );
    cmd.addInt( "R", "uniformRefineT", "Number of Uniform h-refinement loops in time",  numRefineT );
    cmd.addInt( "E", "degreeElevationT",
                "Number of degree elevation steps to perform in time before solving (0: equalize degree in all directions)", numElevateT );
    cmd.addInt( "N", "steps", "Number of time steps", steps );
    cmd.addInt( "t", "testCase", "Test case (0: , 1: )", testCase );
    cmd.addSwitch("last", "Solve solely for the last level of h-refinement",
                  last);
    cmd.addSwitch(
        "plot", "Create a ParaView visualization file with the solution", plot);
    cmd.addSwitch("binary", "Use B64 encoding for Paraview", export_b64);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }
    //! [Parse command line]

    //! [Read input file]

    real_t t_max = 10;

    // directions 0 and 1 are spatial
    // direction 2 is time
    gsMultiPatch<> mpBSpline;
    mpBSpline.addPatch(*gsNurbsCreator<>::BSplineCube());
    // mpBSpline.patch(0).coefs().col(2)*=t_max;

    gsBoundaryConditions<> bc;

    std::string source;
    std::string initial;
    if (testCase==0)
    {
        source = "0";
        std::string r = "10*sqrt((x-0.5)^2+(y-0.5)^2)";
        initial = "if("+r+" <= 1,(1-("+r+")^2)^2,0)";
        bc.addCondition(boundary::west, condition_type::dirichlet, nullptr); //&ms
        bc.addCondition(boundary::east, condition_type::dirichlet, nullptr); //&ms
        bc.addCondition(boundary::south, condition_type::dirichlet, nullptr); //&ms
        bc.addCondition(boundary::north, condition_type::dirichlet, nullptr); //&ms
    }
    else if (testCase==1)
    {
        std::string x0 = "(0.25*cos(2*pi*z)+0.5)";
        std::string y0 = "(0.25*sin(2*pi*z)+0.5)";
        std::string r = "10*sqrt((x-"+x0+")^2+(y-"+y0+")^2)";
        source = "if("+r+" <= 1,(1-("+r+")^2)^2,0)";
        initial = "0";
    }
    else if (testCase==2)
    {
        alpha = 1.0;         //[W/mm/K] ;k in the paper
        beta = 1.0;          //[J/kg/K]*[kg/mm3]; rho*C_h in the paper
        real_t P = 9e5;      //[W]
        real_t eta = 0.33;   //[-]
        real_t r_h = 100e-3/10; //[mm] HMV: DIVIDED BY 10 BECAUSE THE DOMAIN IS SMALLER!
        std::string x0 = "0.5+0.5*cos(0.5*pi*z+pi/4)";
        std::string y0 = "0.3+0.5*sin(0.5*pi*z+pi/4)";
        source = util::to_string(P*eta)+"*exp(-((x-("+x0+"))^2+(y-("+y0+"))^2)/"+util::to_string(pow(r_h,2))+")";
        initial = "20";      //[K]
    }
    else
        GISMO_ERROR("Unknown test case");

    gsFunctionExpr<> f(source,3);
    gsInfo<<"Source function "<< f << "\n";
    gsFunctionExpr<> u0(initial,3);
    gsInfo<<"Initial condition function "<< u0 << "\n";

    // front: initial condition
    bc.addCondition(boundary::front, condition_type::dirichlet, &u0);
    gsInfo<<"Boundary conditions:\n"<< bc <<"\n";
    bc.setGeoMap(mpBSpline);


    mpBSpline.degreeElevate(numElevate,0);
    mpBSpline.degreeElevate(numElevate,1);
    mpBSpline.degreeElevate(numElevateT,2);
    // h-refine each basis
    for (int r =0; r < numRefine; ++r)
    {
        mpBSpline.uniformRefine(1,1,0);
        mpBSpline.uniformRefine(1,1,1);
    }
    for (int r =0; r < numRefineT; ++r)
        mpBSpline.uniformRefine(1,1,2);

    // Convert to THB splines
    gsMultiPatch<> mp;
    for (size_t k=0; k!=mpBSpline.nPatches(); ++k)
    {
        gsTensorBSpline<3,real_t> *geo = dynamic_cast< gsTensorBSpline<3,real_t> * > (&mpBSpline.patch(k));
        gsTHBSpline<3,real_t> thb(*geo);
        mp.addPatch(thb);
    }


    //! [Refinement]
    gsMultiBasis<> dbasis(mp, true);//true: poly-splines (not NURBS)


    // Set heat conduction coefficient
    typedef gsExprAssembler<>::geometryMap geometryMap;
    typedef gsExprAssembler<>::space       space;
    typedef gsExprAssembler<>::solution    solution;

    gsExprAssembler<> A(1,1);

    gsInfo<<"Active options:\n"<< A.options() <<"\n";

    A.setIntegrationElements(dbasis);

    gsExprEvaluator<> ev(A);

    // Set the geometry map
    geometryMap G = A.getMap(mp);

    // Set the discretization space
    space u = A.getSpace(dbasis);

    // Set the source term
    auto ff = A.getCoeff(f, G);

    // Set the solution
    gsMatrix<> solVector;
    solution theta = A.getSolution(u, solVector);

    // Assemble mass matrix
    auto dudt = subgrad(u,2,0);
    auto dudX = subgrad(u,0,1);
    auto dthetadt = subgrad(theta,2,0);
    auto dthetadX = subgrad(theta,0,1);
    auto d2thetadX = sublapl(theta,0,1);

#ifdef GISMO_WITH_PARDISO
    gsSparseSolver<>::PardisoLDLT solver;
#else
    // gsSparseSolver<>::CGDiagonal solver;
    gsSparseSolver<>::QR solver;
#endif

    gsAdaptiveMeshing<3,real_t> mesher(dbasis);
    mesher.options().setInt("RefineRule",3);
    mesher.options().setInt("CoarsenRule",3);
    mesher.options().setSwitch("Admissible",true);
    mesher.options().setReal("RefineParam",0.2);
    // mesher.options().setReal("CoarsenParam",0.001);
    mesher.options().setInt("MaxLevel",8);
    mesher.getOptions();
    gsHBoxContainer<3,real_t> refine,coarsen;

    gsParaviewCollection collection("ParaviewOutput/solution", &ev);
    if (plot)
    {
        collection.options().setSwitch("plotElements", true);
        collection.options().setSwitch("base64", export_b64);
        collection.options().setInt("plotElements.resolution", 16);
        collection.options().setInt("numPoints",1e5);
    }

    index_t it=0;
    real_t E_I, E_I_prev, E_T, E_T_prev;
    E_I_prev = E_T_prev = 1;
    real_t error, initialError;
    error = initialError = 1;
    std::vector<real_t> errors;

    gsInfo<<std::setw(12)<<"Iteration"<<std::setw(12)<<"Error"<<std::setw(12)<<"E_I"<<std::setw(12)<<"E_T"<<std::setw(12)<<"Rel. error"<<std::setw(12)<<"Rel. E_I"<<std::setw(12)<<"Rel. E_T"<<std::setw(12)<<"DoFs"<<"\n";
    while (error/initialError > 1e-3 && it < 20)
    {
        u.setup(bc, dirichlet::l2Projection, 0);
        A.initSystem();
        A.assemble( (u*dudt.tr() + alpha*dudX*dudX.tr()) * meas(G), u*ff*meas(G));
        solVector = solver.compute(A.matrix()).solve(A.rhs());

        error = ev.integralElWise((ff-beta*dthetadt+alpha*d2thetadX).sqNorm()*meas(G));
        errors = ev.elementwise();
        E_I = 0.5*ev.integral(alpha*dthetadX*dthetadX.tr()*meas(G));
        E_T = E_I + 0.5*ev.integral(beta*theta*dthetadt*meas(G));
        gsInfo  <<std::setw(12)<<it
                <<std::setw(12)<<error
                <<std::setw(12)<<E_I
                <<std::setw(12)<<E_T
                <<std::setw(12)<<error/initialError
                <<std::setw(12)<<(E_I-E_I_prev)/E_I
                <<std::setw(12)<<(E_T-E_T_prev)/E_T
                <<std::setw(12)<<A.numDofs()
                <<"\n";

        if (plot)
        {
            collection.newTimeStep(&mp);
            collection.addField(theta,"numerical solution");
            collection.addField(ff,"source");
            collection.addField((ff-beta*dthetadt+alpha*d2thetadX).sqNorm(),"residual");
            collection.saveTimeStep();
        }

        if (error/initialError < 1e-3)
        {
            gsInfo<<"Converged\n";
            break;
        }

        if (it + 1 == 20)
        {
            gsInfo<<"Did not converge\n";
            break;
        }

        mesher.markRef_into(errors,refine);
        // mesher.markCrs_into(errors,refine,coarsen);
        mesher.refine(refine);
        // mesher.unrefine(coarsen);
        mesher.rebuild();

        E_I_prev = E_I;
        E_T_prev = E_T;
        if (it==0)
            initialError = error;
        it++;

    }

    if (plot)
        collection.save();





    return  EXIT_SUCCESS;
}
