/** @file composed_domain_poisson.cpp

    @brief Tutorial on how to use expression assembler to solve the Poisson equation

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris
*/

//! [Include namespace]
#include <gismo.h>
#include <gsNurbs/gsSquareDomain.h>

using namespace gismo;
//! [Include namespace]

int main(int argc, char *argv[])
{
    //! [Parse command line]
    bool plot = false;
    bool plotbasis = false;
    index_t numRefine  = 1;
    index_t numElevate = 0;

    gsCmdLine cmd("Tutorial on solving a Poisson problem.");
    cmd.addInt( "e", "degreeElevation",
                "Number of degree elevation steps to perform before solving (0: equalize degree in all directions)", numElevate );
    cmd.addInt( "r", "uniformRefine", "Number of Uniform h-refinement loops",  numRefine );
    cmd.addSwitch("plot", "Create a ParaView visualization file with the solution", plot);
    cmd.addSwitch("plotB", "Create a ParaView visualization file with the solution", plotbasis);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }
    //! [Parse command line]

    gsMultiPatch<> mp0;
    mp0.addPatch(gsNurbsCreator<>::BSplineSquare());
    mp0.patch(0).coefs().array() -= 0.5;
    mp0.embed(3);

    if (numElevate!=0)
        mp0.degreeElevate(numElevate);

    // h-refine
    for (int r =0; r < numRefine; ++r)
        mp0.uniformRefine();


    // Make composed geometry and basis
    const gsBasis<> & tbasis = mp0.basis(0); // basis(u,v) -> deriv will give dphi/du ,dphi/dv
    const gsGeometry<> & tgeom = mp0.patch(0); //G(u,v) -> deriv will give dG/du, dG/dv

    // The domain sigma
    // gsSquareDomain<2,real_t> domain;

    // gsMatrix<> pars = domain.controls();
    // pars *= 0.99;
    // // pars(0,0) -= 0.1;
    // domain.controls() = pars.col(0);
    // domain.updateGeom();

    gsFunctionExpr<> domain("(x)^(2)","(y)^(2)",2);
    // gsFunctionExpr<> domain("(x","y",2);


    // Define a composite basis and composite geometry
    // The basis is composed by the square domain
    gsComposedBasis<real_t> cbasis(domain,tbasis); // basis(u,v) = basis(sigma(xi,eta)) -> deriv will give dphi/dxi, dphi/deta
    // The geometry is defined using the composite basis and some coefficients
    gsComposedGeometry<real_t> cgeom(domain, tgeom); // G(u,v) = G(sigma(xi,eta))  -> deriv will give dG/dxi, dG/deta

    if (plotbasis)
    {
        gsWriteParaview(cbasis,"cbasis");
        gsWriteParaview(tbasis,"tbasis");
    }

    // const gsBasis<> & cbasis = tbasis; // basis(u,v) -> deriv will give dphi/du ,dphi/dv
    // const gsGeometry<> & cgeom = tgeom;


    gsMultiPatch<> mp;
    mp.addPatch(cgeom);

    gsMultiBasis<> dbasis(cbasis);

    // gsMultiBasis<> dbasis(mp, true);

    //! [Refinement]

     // Source function:
     gsFunctionExpr<> f("((tanh(20*(x^2 + y^2)^(1/2) - 5)^2 - 1)*(20*x^2 + 20*y^2)*(40*tanh(20*(x^2 + y^2)^(1/2) - 5)*(x^2 + y^2)^(1/2) - 1))/(x^2 + y^2)^(3/2)",3);

     // Exact solution
     gsFunctionExpr<> ms("tanh((0.25-sqrt(x^2+y^2))/0.05)+1",3);

    // Source function:
//    gsFunctionExpr<> f("2*pi^2*cos(pi*x)*cos(pi*y)",2);

    // Exact solution
//    gsFunctionExpr<> ms("cos(pi*x)*cos(pi*y)",2);

    gsBoundaryConditions<> bc;
    bc.addCondition(boundary::side::west ,condition_type::dirichlet,&ms);
    bc.addCondition(boundary::side::east ,condition_type::dirichlet,&ms);
    bc.addCondition(boundary::side::south,condition_type::dirichlet,&ms);
    bc.addCondition(boundary::side::north,condition_type::dirichlet,&ms);
    bc.setGeoMap(mp);

    //! [Problem setup]
    gsExprAssembler<> A(1,1);

    typedef gsExprAssembler<>::geometryMap geometryMap;
    typedef gsExprAssembler<>::variable    variable;
    typedef gsExprAssembler<>::space       space;
    typedef gsExprAssembler<>::solution    solution;

    // Elements used for numerical integration
    A.setIntegrationElements(dbasis);
    gsExprEvaluator<> ev(A);

    // Set the geometry map
    geometryMap G = A.getMap(mp);

    // Set the discretization space
    space u = A.getSpace(dbasis);

    // Set the source term
    auto ff = A.getCoeff(f, G);

    // Recover manufactured solution
    auto u_ex = ev.getVariable(ms, G);

    // Solution vector and solution variable
    gsMatrix<> solVector;
    solution u_sol = A.getSolution(u, solVector);

    //! [Problem setup]

    gsSparseSolver<>::CGDiagonal solver;

    u.setup(bc, dirichlet::l2Projection, 0);

    // Initialize the system
    A.initSystem();

    gsInfo<< "A.numDofs() = " << A.numDofs() <<std::flush;

    // Compute the system matrix and right-hand side
   A.assemble(
       igrad(u, G) * igrad(u, G).tr() * meas(G) //matrix
       ,
       u * ff * meas(G) //rhs vector
       );

// //  grad(u)*jac(G).ginv()
//   A.assemble(
//       (grad(u)*(jac(G).ginv().tr())) * (grad(u)*(jac(G).ginv().tr())).tr() * meas(G) //matrix
//       ,
//       u * ff * meas(G) //rhs vector
//   );

    solver.compute( A.matrix() );
    solVector = solver.solve(A.rhs());

    // Compute the error
    real_t L2err = math::sqrt( ev.integral( (u_ex - u_sol).sqNorm() * meas(G) ) );
    gsInfo<<"\nL2 error = "<<L2err<<"\n";

    //! [Export visualization in ParaView]
    if (plot)
    {
        gsInfo<<"Plotting in Paraview...\n";

        gsParaviewCollection collection("ParaviewOutput/solution", &ev);
        collection.options().setSwitch("plotElements", true);
        collection.options().setInt("plotElements.resolution", 16);
        collection.options().setInt("numPoints", 1000);
        collection.options().setInt("precision", 12);
//        collection.options().setInt("plotElements.resolution", 16);
        collection.newTimeStep(&mp);
        collection.addField(u_sol,"numerical solution");
        collection.addField(u_ex, "exact solution");
        collection.addField((u_ex-u_sol).sqNorm(), "error");
        collection.saveTimeStep();
        collection.save();

        gsWriteParaview(cgeom,"cgeom");

        // gsFileManager::open("ParaviewOutput/solution.pvd");
    }
    else
        gsInfo << "Done. No output created, re-run with --plot to get a ParaView "
                  "file containing the solution.\n";
    //! [Export visualization in ParaView]




//    gsDebug<<"GEOMETRY==================================================================\n";
//
//    gsMatrix<> point(2,1);
//    point.col(0) << 0.5, 0.25;
//    gsDebugVar(ev.eval(jac(G), point));
//
////    gsDebugVar( ev.integral(jac(G).det()) );
////    gsDebugVar( ev.integral( meas(G) ) );
//
//    gsVector<> pp = point.col(0);
//    gsMatrix<> ev1, ev2;
//    gsMatrix<> der;
//
//    real_t delta = 1e-6;
//    gsVector<> pt1 = pp, pt2 = pp;
//    pt1.at(1) += delta;
//    pt2.at(1) -= delta;
//    cgeom.eval_into(pt1,ev1);
//    cgeom.eval_into(pt2,ev2);
//    gsVector<> der_eta = (ev1-ev2)/(2*delta);
//    gsDebugVar(der_eta);
//
//    pt1 = pp, pt2 = pp;
//    pt1.at(0) += delta;
//    pt2.at(0) -= delta;
//    cgeom.eval_into(pt1, ev1);
//    cgeom.eval_into(pt2, ev2);
//    gsVector<> der_xi = (ev1-ev2)/(2*delta);
//    gsDebugVar(der_xi);
//
//    cgeom.deriv_into(pp,der);
////  gsDebugVar(der);
//    gsDebugVar(der.reshape(2,2));

//    gsDebug<<"BASIS==================================================================\n";
//    gsMatrix<index_t> act;
//    cbasis.active_into(pp,act);
//    cbasis.evalSingle_into(act(0,0),pt1,ev1);
//    cbasis.evalSingle_into(act(0,0),pt2,ev2);
//
//    gsVector<> der2 = (ev1-ev2)/(2*delta);
//    gsDebugVar(der2);
//
//    cbasis.derivSingle_into(act(0,0),pp,der);
//    gsDebugVar(der);

    return EXIT_SUCCESS;

}// end main
