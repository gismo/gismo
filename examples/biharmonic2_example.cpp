/** @file biharmonic2_example.cpp

    @brief Tutorial on how to use expression assembler to solve the Poisson equation

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris
*/

//! [Include namespace]
#include <gismo.h>

#include <gsUnstructuredSplines/gsApproxC1Spline.h>

using namespace gismo;
//! [Include namespace]

int main(int argc, char *argv[])
{
    //! [Parse command line]
    bool plot = false;
    index_t numRefine  = 3;
    index_t numElevate = 0;
    bool last = false;
    std::string fn;

    index_t geometry = 1000;

    gsCmdLine cmd("Tutorial on solving a Poisson problem.");
    cmd.addInt( "e", "degreeElevation",
                "Number of degree elevation steps to perform before solving (0: equalize degree in all directions)", numElevate );
    cmd.addInt( "r", "uniformRefine", "Number of Uniform h-refinement loops",  numRefine );
    cmd.addString( "f", "file", "Input geometry file", fn );
    cmd.addInt( "g", "geometry", "Which geometry",  geometry );
    cmd.addSwitch("last", "Solve solely for the last level of h-refinement", last);
    cmd.addSwitch("plot", "Create a ParaView visualization file with the solution", plot);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }
    //! [Parse command line]

    //! [Read geometry]
    std::string string_geo;
    if (fn.empty())
        string_geo = "planar/geometries/g" + util::to_string(geometry) + ".xml";
    else
        string_geo = fn;

    gsMultiPatch<> mp;
    gsInfo << "Filedata: " << string_geo << "\n";
    gsReadFile<>(string_geo, mp);
    mp.clearTopology();
    mp.computeTopology();

    std::string string_file = "planar/biharmonic_pde/bvp1.xml";
    gsFileData<> fd(string_file);
    gsFunctionExpr<> f, laplace, ms;
    fd.getId(0,f);
    gsInfo<<"Source function "<< f << "\n";

    fd.getId(1,ms); // Exact solution
    gsInfo<<"Exact function "<< ms << "\n";

    fd.getId(2,laplace); // Laplace for the bcs
    gsInfo<<"Finished\n";

    //! [Boundary condition]
    gsBoundaryConditions<> bc, bc2;
    for (gsMultiPatch<>::const_biterator bit = mp.bBegin(); bit != mp.bEnd(); ++bit)
    {
        bc.addCondition(*bit, condition_type::dirichlet, &ms);
        bc2.addCondition(*bit, condition_type::neumann, &laplace); // Is not the usually neumann condition
    }
    bc.setGeoMap(mp);
    bc2.setGeoMap(mp);
    //! [Boundary condition]
    gsInfo<<"First boundary conditions:\n"<< bc <<"\n";
    gsInfo<<"Second boundary conditions:\n"<< bc2 <<"\n";

    gsOptionList optionList;
    fd.getId(100, optionList); // id=100: assembler options
    gsInfo << "OptionList: " << optionList << "\n";
    //! [Read input file]

    //! [Refinement]
    gsMultiBasis<> dbasis(mp, false);//true: poly-splines (not NURBS)

    // Elevate and p-refine the basis to order p + numElevate
    // where p is the highest degree in the bases
    if (dbasis.maxCwiseDegree() + numElevate < 3)
        numElevate = 3 - dbasis.maxCwiseDegree();
    dbasis.setDegree( dbasis.maxCwiseDegree() + numElevate);

    // h-refine each basis
    if (last)
    {
        for (int r =0; r < numRefine; ++r)
            dbasis.uniformRefine();
        numRefine = 0;
    }

    gsInfo << "Patches: "<< mp.nPatches() <<", degree: "<< dbasis.minCwiseDegree() <<"\n";
#ifdef _OPENMP
    gsInfo<< "Available threads: "<< omp_get_max_threads() <<"\n";
#endif
    //! [Refinement]

    //! [Problem setup]
    gsExprAssembler<> A(1,1);
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

    // Set the source term
    auto ff = A.getCoeff(f, G); // Laplace example

    // Set the discretization space
    gsMappedBasis<2,real_t> bb2;
    space u = A.getSpace(bb2);

    // Solution vector and solution variable
    gsMatrix<> solVector;
    solution u_sol = A.getSolution(u, solVector);

    // Recover manufactured solution
    auto u_ex = ev.getVariable(ms, G);
    //! [Problem setup]

    //! [Solver loop]
    gsSparseSolver<>::CGDiagonal solver;

    gsVector<> l2err(numRefine+1), h1err(numRefine+1), h2err(numRefine+1);
    gsInfo<< "(dot1=assembled, dot2=solved, dot3=got_error)\n"
        "\nDoFs: ";
    double setup_time(0), ma_time(0), slv_time(0), err_time(0);
    gsStopwatch timer;
    for (int r=0; r<=numRefine; ++r)
    {
        dbasis.uniformRefine();

        // Compute the approx C1 space
        gsApproxC1Spline<2,real_t> approxC1(mp,dbasis);
        //approxC1.options().setSwitch("info",info);
        approxC1.options().setSwitch("plot",plot);

        approxC1.init();
        approxC1.compute();

        gsSparseMatrix<> global2local = approxC1.getSystem();
        global2local = global2local.transpose();
        global2local.pruned(1,1e-10);
        gsMultiBasis<> dbasis_temp;
        approxC1.getMultiBasis(dbasis_temp);
        bb2.init(dbasis_temp,global2local);

        // Setup the system
        u.setup(bc, dirichlet::homogeneous, -1);

        // Initialize the system
        A.initSystem();
        setup_time += timer.stop();

        gsInfo<< A.numDofs() <<std::flush;

        timer.restart();
        // Compute the system matrix and right-hand side
        A.assemble(
                ilapl(u, G) * ilapl(u, G).tr()
                * meas(G)
                ,
                u * ff
                * meas(G)
        );

        // Enforce Laplace conditions to right-hand side
        auto g_L = A.getCoeff(laplace, G); // Set the laplace bdy value
        A.assembleRhsBc( (igrad(u, G) * nv(G)) * g_L.tr(), bc2.neumannSides() );

        ma_time += timer.stop();
        gsInfo<< "." <<std::flush;// Assemblying done

        timer.restart();
        solver.compute( A.matrix() );
        solVector = solver.solve(A.rhs());

        slv_time += timer.stop();
        gsInfo<< "." <<std::flush; // Linear solving done

        timer.restart();
        l2err[r]= math::sqrt( ev.integral( (u_ex - u_sol).sqNorm() * meas(G) ) );
        
        h1err[r]= l2err[r] +
            math::sqrt(ev.integral( ( igrad(u_ex) - igrad(u_sol,G) ).sqNorm() * meas(G) ));

        //h2err[r]= l2err[r] + h1err[r] +
        //          math::sqrt(ev.integral( ( ihess(u_ex) - ihess(u_sol,G) ).sqNorm() * meas(G) ));

        err_time += timer.stop();
        gsInfo<< ". " <<std::flush; // Error computations done

    } //for loop

    //! [Solver loop]

    timer.stop();
    gsInfo<<"\n\nTotal time: "<< setup_time+ma_time+slv_time+err_time <<"\n";
    gsInfo<<"     Setup: "<< setup_time <<"\n";
    gsInfo<<"  Assembly: "<< ma_time    <<"\n";
    gsInfo<<"   Solving: "<< slv_time   <<"\n";
    gsInfo<<"     Norms: "<< err_time   <<"\n";

    //! [Error and convergence rates]
    gsInfo<< "\nL2 error: "<<std::scientific<<std::setprecision(3)<<l2err.transpose()<<"\n";
    gsInfo<< "H1 error: "<<std::scientific<<h1err.transpose()<<"\n";
    gsInfo<< "H2 error: "<<std::scientific<<h2err.transpose()<<"\n";

    if (!last && numRefine>0)
    {
        gsInfo<< "\nEoC (L2): " << std::fixed<<std::setprecision(2)
              <<  ( l2err.head(numRefine).array()  /
                   l2err.tail(numRefine).array() ).log().transpose() / std::log(2.0)
                   <<"\n";

        gsInfo<<   "EoC (H1): "<< std::fixed<<std::setprecision(2)
              <<( h1err.head(numRefine).array() /
                  h1err.tail(numRefine).array() ).log().transpose() / std::log(2.0) <<"\n";

        gsInfo<<   "EoC (H2): "<< std::fixed<<std::setprecision(2)
              <<( h2err.head(numRefine).array() /
                  h2err.tail(numRefine).array() ).log().transpose() / std::log(2.0) <<"\n";
    }
    //! [Error and convergence rates]

    //! [Export visualization in ParaView]
    if (plot)
    {
        gsInfo<<"Plotting in Paraview...\n";
        ev.options().setSwitch("plot.elements", false);
        ev.writeParaview( u_sol   , G, "solution");
        //ev.writeParaview( u_ex    , G, "solution_ex");
        //ev.writeParaview( u, G, "aa");

        //gsFileManager::open("solution.pvd");
    }
    else
        gsInfo << "Done. No output created, re-run with --plot to get a ParaView "
                  "file containing the solution.\n";
    //! [Export visualization in ParaView]

    return EXIT_SUCCESS;

}// end main
