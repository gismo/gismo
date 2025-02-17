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
    bool last{false}, export_b64{false};

    gsCmdLine cmd("Tutorial on solving a Poisson problem.");
    cmd.addInt( "e", "degreeElevation",
                "Number of degree elevation steps to perform before solving (0: equalize degree in all directions)", numElevate );
    cmd.addInt( "r", "uniformRefine", "Number of Uniform h-refinement loops",  numRefine );
    cmd.addInt( "R", "uniformRefineT", "Number of Uniform h-refinement loops in time",  numRefineT );
    cmd.addInt( "E", "degreeElevationT",
                "Number of degree elevation steps to perform in time before solving (0: equalize degree in all directions)", numElevateT );
    cmd.addInt( "N", "steps", "Number of time steps", steps );
    cmd.addSwitch("last", "Solve solely for the last level of h-refinement",
                  last);
    cmd.addSwitch(
        "plot", "Create a ParaView visualization file with the solution", plot);
    cmd.addSwitch("binary", "Use B64 encoding for Paraview", export_b64);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }
    //! [Parse command line]

    //! [Read input file]

    real_t t_max = 10;
    real_t lambda = 1/(32*pow(EIGEN_PI,2));

    // directions 0 and 1 are spatial
    // direction 2 is time
    gsMultiPatch<> mp;
    mp.addPatch(*gsNurbsCreator<>::BSplineCube());
    // mp.patch(0).coefs().col(2)*=t_max;

    // std::string r = "10*sqrt((x-0.5)^2+(y-0.5)^2)";
    // gsFunctionExpr<> initial("if("+r+" <= 1,(1-("+r+")^2)^2,0)",3);
    gsFunctionExpr<> initial("0.1 * cos(2*pi*x) * cos(2*pi*y)",3);
    gsInfo<<"Initial condition function "<< initial << "\n";

    // front: initial condition
    gsBoundaryConditions<> bc;
    bc.addCondition(boundary::front, condition_type::dirichlet, &initial);
    bc.addCondition(boundary::north, condition_type::clamped, nullptr,0,false,0);
    bc.addCondition(boundary::east , condition_type::clamped, nullptr,0,false,1);
    bc.addCondition(boundary::south, condition_type::clamped, nullptr,0,false,0);
    bc.addCondition(boundary::west , condition_type::clamped, nullptr,0,false,1);

    gsInfo<<"Boundary conditions:\n"<< bc <<"\n";
    bc.setGeoMap(mp);


    mp.degreeElevate(numElevate,0);
    mp.degreeElevate(numElevate,1);
    mp.degreeElevate(numElevateT,2);
    // h-refine each basis
    for (int r =0; r < numRefine; ++r)
    {
        mp.uniformRefine(1,1,0);
        mp.uniformRefine(1,1,1);
    }
    for (int r =0; r < numRefineT; ++r)
        mp.uniformRefine(1,1,2);

    //! [Refinement]
    gsMultiBasis<> dbasis(mp, true);//true: poly-splines (not NURBS)
    gsInfo<<"Basis:\n"<< dbasis <<"\n";
    gsInfo<<"Size: "<< dbasis.size() <<"\n";

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
    space w = A.getSpace(dbasis);
    w.setup(bc, dirichlet::l2Projection, 0);

    auto dwdt = subgrad(w,2,0);
    auto dwdX = subgrad(w,0,1);
    auto d2wdX= sublapl(w,0,1);
    // auto d2wdX= lapl(w);

    gsMatrix<> C;
    auto c  = A.getSolution(w, C);
    // auto dc = A.getSolution(u, dC);

    gsWarn<<"How to handle ijac?\n";
    auto dcdt = subgrad(c,2,0);
    auto dcdX = subgrad(c,0,1);
    auto d2cdX= sublapl(c,0,1);

    // Derivatives of the double well potential (Gomez et al., 2008)
    auto dmu_c = - 1.0 + 3.0 * (c*c).val(); // f_2 (second derivative of double well)
    auto ddmu_c = 6*c.val(); // f_3 (third derivative of double well)

    // Mobility
    auto M_c  = 1.0 + 0.0*c.val(); // replace with const_expr(1.0) instead of using 0*c
    // auto dM_c = 0.0 * igrad(c,G); // replace with const_expr(1.0) instead of using 0*c!!

    auto lam  = lambda + 0.0*c.val(); // replace with const_expr(lambda) instead of using 0*c

/*
    auto residual =
                    w*dcdt.tr()
                    +
                    alpha*dwdX  * dcdX.tr()
                    +
                    d2wdX * d2cdX.val()
                    ;
    auto jacobian = //
                    w*dwdt.tr()
                    +
                    alpha*dwdX * dwdX.tr()
                    +
                    d2wdX * d2wdX.tr()
                    ;
 */



    auto residual = //
                    w*dcdt.tr()
                    +
                    dwdX  * dmu_c * dcdX.tr()
                    + // F_bar
                    d2wdX * lam   * d2cdX.val()
                    ; // K_laplacian

    auto jacobian = //
                    w*dwdt.tr()
                    +
                    dwdX * dmu_c  * dwdX.tr()
                    + // K_f1
                    dwdX * ddmu_c * dcdX.tr() * w.tr()
                    + // K_f2
                    d2wdX* lam    * d2wdX.tr()
                    ; // K_laplacian

    gsMatrix<> dC;
    gsMatrix<> Q;

    // Define linear solver (install SuperLUMT-devel)
#ifdef GISMO_WITH_PARDISO
    gsSparseSolver<>::PardisoLDLT solver;
#else
    // gsSparseSolver<>::QR solver;
    // gsSparseSolver<>::CGDiagonal solver;
    gsSparseSolver<>::LU solver;
#endif

    A.initSystem();
    C.setZero(A.numDofs(),1);
    // C.setRandom(A.numDofs(),1);
    // C *= 0.01;

    // // Linear solve
    // A.assemble(jacobian * meas(G), w*0);
    // solver.compute(A.matrix());
    // C = solver.solve(A.rhs());
    // ev.writeParaview(c,G,"initial");

    // w.setup(bc, dirichlet::homogeneous, 0);

    index_t maxIter = 50;
    real_t R0 = 0;
    gsInfo<<"It.\t|R|/|R0|\t\t|U|\t\t|dU|/|U|\n";
    for (index_t iter = 0; iter < maxIter; ++iter)
    {
        A.clearMatrix();
        A.clearRhs();

        // Assemble residual
        A.assemble(residual * meas(G));
        if (iter==0) R0 = A.rhs().norm();
        if (A.rhs().norm()/R0 < 1e-6 || iter == maxIter-1)
        {
            gsInfo<<"Converged with residual = "<<Q.norm()<<"\n";
            break;
        }


        // Assemble Jacobian
        // A.assembleJacobian(residual * meas(G), c);
        A.assemble(jacobian * meas(G));
        // Solve
        solver.compute(A.matrix());
        dC = solver.solve(-A.rhs());
        C += dC;
        gsInfo<<iter<<"\t"<<A.rhs().norm()<<"\t"<<C.norm()<<"\t"<<dC.norm()/C.norm()<<"\n";
        if ((dC.norm()/C.norm() < 1e-3 && A.rhs().norm()/R0 < 1e-6) || iter == maxIter-1)
        {
            gsInfo<<"Converged with update = "<<dC.norm()<<"\n";
            break;
        }
    }

    gsDebugVar(C.transpose());

    // w.setup(bc, dirichlet::l2Projection, 0);

    if (plot)
    {
        gsInfo<<"Plotting in Paraview...\n";

        gsParaviewCollection collection("ParaviewOutput/solution", &ev);
        collection.options().setSwitch("plotElements", true);
        collection.options().setSwitch("base64", export_b64);
        collection.options().setInt("plotElements.resolution", 16);
        collection.options().setInt("numPoints",1e5);
        collection.newTimeStep(&mp);
        collection.addField(c,"numerical solution");
        collection.saveTimeStep();
        collection.save();
    }
    return  EXIT_SUCCESS;
}
