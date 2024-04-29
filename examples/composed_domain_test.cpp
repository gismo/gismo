/** @file composed_domain_test.cpp

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
    // a
    bool b_morph = false;   // b, composition on the basis
    // c
    index_t spaceDim = 2;   // d, dimension of the physical solution
    index_t numElevate = 0; // e, deg elevation
    // f
    bool g_morph = false;   // g, composition on the geometry
    // i, j, k, l, m, n, o
    index_t problem = 0;    // p, problem: 0 --> poisson, 1 -->  L2, 2 --> another pde;
    // q
    index_t numRefine  = 1; // r, uniform refiment
    // s, t, u, v, w, x, y, z
    bool ref_last = false;  // last, perform last refinement only
    bool plot = false;
    bool plotbasis = false;


    gsCmdLine cmd("Tutorial on solving a Poisson problem.");
    cmd.addInt( "e", "degreeElevation",
                "Number of degree elevation steps to perform before solving (0: equalize degree in all directions)", numElevate );
    cmd.addInt("d", "dimension", "specify the problem dimension", spaceDim);
    cmd.addInt("p", "equation", "specify the problem", problem);
    cmd.addInt( "r", "uniformRefine", "Number of Uniform h-refinement loops",  numRefine );
    cmd.addSwitch("plot", "Create a ParaView visualization file with the solution", plot);
    cmd.addSwitch("plotB", "Create a ParaView visualization file with the solution", plotbasis);
    cmd.addSwitch("b", "basis", "Apply composition to the basis", b_morph);
    cmd.addSwitch("g", "geom", "Apply composition to the geometry", g_morph);
    cmd.addSwitch("last", "one-shot uniform refinement", ref_last);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }
    //! [Parse command line]

    gsStopwatch gsTime;
    time_t now = time(0);

    std::ofstream file_out;
    file_out.open(std::to_string(now)+"composedDomain_results.csv");
    file_out << "problem, cbasis, cgeom, deg, ref, dofs, L2err\n";

    std::string problem_name = "";
    if(problem == 0)
    {
      problem_name = "poisson";
      gsInfo << "Solve Poisson equation.\n";
    }
    else if(problem == 1)
    {
      problem_name = "L2_projection";
      gsInfo << "Solve L2 projection.\n";
    }
    else
    {
      gsWarn << "Unknown problem, exiting." << std::endl;
      return -1;
    }


    gsMultiPatch<> mp0;
    mp0.addPatch(gsNurbsCreator<>::BSplineSquare());
    mp0.patch(0).coefs().array() -= 0.5;
    if (spaceDim > 2)
      mp0.embed(spaceDim);

    if (numElevate!=0)
        mp0.degreeElevate(numElevate);

    gsInfo << "degree: (" << mp0.patch(0).degree(0) << "," << mp0.patch(0).degree(1) << ")\n";

    for(index_t refCount = 1; refCount <= numRefine; refCount++)
    {

        if(ref_last)
        {
          refCount = numRefine+1;
          index_t numKnts = pow(2, numRefine) -1;
          mp0.uniformRefine(numKnts);
        }
        else
          mp0.uniformRefine();

        // Make geometry and basis
        const gsBasis<> & tbasis = mp0.basis(0); // basis(u,v) -> deriv will give dphi/du ,dphi/dv
        const gsGeometry<> & tgeom = mp0.patch(0); //G(u,v) -> deriv will give dG/du, dG/dv

        // The domain sigma
        gsSquareDomain<2,real_t> domain(1, 1);
        gsMatrix<> pars = domain.controls();

        gsInfo << domain << "\n";
        gsInfo << "Domain coefficients:\n" << pars << "\n";

        if(b_morph || g_morph)
        {
          gsInfo << "Apply sigma.\n";
          pars *= 0.8;
          // pars(0,0) -= 0.1;
          domain.controls() = pars.col(0);
        }
        domain.updateGeom();

        //gsFunctionExpr<> domain("(x)^(2)","(y)^(2)",2);
        //gsFunctionExpr<> domain("(x","y",2);

        gsBasis<>::Ptr cbasis;
        // Define a composite basis and composite geometry
        if(b_morph)
        {
          gsInfo << "Composed basis.\n";
          // gsComposedBasis<real_t> cbasis(domain,tbasis); // basis(u,v) = basis(sigma(xi,eta)) -> deriv will give dphi/dxi, dphi/deta
          cbasis = memory::make_shared(new gsComposedBasis<real_t>(domain,tbasis));
        }
        else
        {
          cbasis = memory::make_shared_not_owned(&tbasis);
          // const gsBasis<> & cbasis = tbasis; // basis(u,v) -> deriv will give dphi/du ,dphi/dv
        }


        // const gsGeometry<> & cgeom = tgeom;
        gsGeometry<>::Ptr cgeom;
        if(g_morph)
        {
          gsInfo << "Composed geometry.\n";
          cgeom = memory::make_shared(new gsComposedGeometry<real_t>(domain,tgeom));
          //gsComposedGeometry<real_t> cgeom(domain, tgeom); // G(u,v) = G(sigma(xi,eta))  -> deriv will give dG/dxi, dG/deta
        }
        else
          cgeom = memory::make_shared_not_owned(&tgeom);


        if (plotbasis)
        {
            gsWriteParaview(*cbasis,"cbasis");
            gsWriteParaview(tbasis,"tbasis");
        }


        gsMultiPatch<> mp;
        mp.addPatch(*cgeom);

        gsMultiBasis<> dbasis(*cbasis);

        // gsMultiBasis<> dbasis(mp, true);

        //! [Refinement]

         // Source function:
         gsFunctionExpr<> f("((tanh(20*(x^2 + y^2)^(1/2) - 5)^2 - 1)*(20*x^2 + 20*y^2)*(40*tanh(20*(x^2 + y^2)^(1/2) - 5)*(x^2 + y^2)^(1/2) - 1))/(x^2 + y^2)^(3/2)",spaceDim);

         // Exact solution
         gsFunctionExpr<> ms("tanh((0.25-sqrt(x^2+y^2))/0.05)+1",spaceDim);

        // Source function:
    //    gsFunctionExpr<> f("2*pi^2*cos(pi*x)*cos(pi*y)",spaceDim);

        // Exact solution
    //    gsFunctionExpr<> ms("cos(pi*x)*cos(pi*y)",spaceDim);

        gsBoundaryConditions<> bc;
        if (problem == 0)
        {
          gsInfo << "Add bc for Poisson.\n";
          bc.addCondition(boundary::side::west ,condition_type::dirichlet,&ms);
          bc.addCondition(boundary::side::east ,condition_type::dirichlet,&ms);
          bc.addCondition(boundary::side::south,condition_type::dirichlet,&ms);
          bc.addCondition(boundary::side::north,condition_type::dirichlet,&ms);
        }
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

        gsInfo<< "A.numDofs() = " << A.numDofs() << "\n";

        // Compute the system matrix and right-hand side
        if (problem == 0)
        {
          gsInfo << "Assemble poisson.\n";
          A.assemble( igrad(u, G) * igrad(u, G).tr() * meas(G) //matrix
           ,
           u * ff * meas(G) //rhs vector
           );
        }
        else if(problem == 1)
        {
          gsInfo << "Assemble L2 projection.\n";
          A.assemble(u * u.tr() * meas(G), u * u_ex * meas(G));
        }
        else
        {
          gsWarn << "Unknown problem, exiting." << std::endl;
          return -1;
        }


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

        // file_out << "problem, morph, deg, ref, dofs, L2err\n";
        file_out << problem_name << "," << b_morph << "," << g_morph << ","
                 << std::max(mp0.patch(0).degree(0),mp0.patch(0).degree(1))   << "," << refCount << ","
                 << A.numDofs()  << "," << L2err     << "\n";

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

            gsWriteParaview(*cgeom,"cgeom");

            // gsFileManager::open("ParaviewOutput/solution.pvd");
        }
        else
            gsInfo << "Done. No output created, re-run with --plot to get a ParaView "
                      "file containing the solution.\n";
        //! [Export visualization in ParaView]
    } // refCount
    file_out.close();


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
