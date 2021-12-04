/** @file biharmonic2_example.cpp

    @brief Tutorial on how to use expression assembler and the (approx.) C1 basis function
                to solve the Biharmonic equation

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): P. Weinmueller
*/

//! [Include namespace]
#include <gismo.h>

#include <gsUnstructuredSplines/gsApproxC1Spline.h>
//#include <gsUnstructuredSplines/gsDPatch.h>

#include <gsIO/gsCSVWriter.h>

using namespace gismo;
//! [Include namespace]

int main(int argc, char *argv[])
{
    //! [Parse command line]
    bool plot = false;
    index_t smoothing = 2;

    index_t numRefine  = 3;
    index_t discreteDegree = 3;
    index_t discreteRegularity = 2;
    bool last = false;
    bool info = false;
    bool neumann = false;
    std::string fn;

    index_t geometry = 1000;

    gsCmdLine cmd("Tutorial on solving a Poisson problem.");
    cmd.addInt( "s", "smoothing","Smoothing", smoothing );
    cmd.addInt( "p", "discreteDegree","Which discrete degree?", discreteDegree );
    cmd.addInt( "r", "discreteRegularity", "Number of discreteRegularity",  discreteRegularity );
    cmd.addInt( "l", "refinementLoop", "Number of refinementLoop",  numRefine );
    cmd.addString( "f", "file", "Input geometry file", fn );
    cmd.addInt( "g", "geometry", "Which geometry",  geometry );
    cmd.addSwitch("last", "Solve solely for the last level of h-refinement", last);
    cmd.addSwitch("plot", "Create a ParaView visualization file with the solution", plot);
    cmd.addSwitch("info", "Getting the information inside of Approximate C1 basis functions", info);

    cmd.addSwitch("neumann", "Neumann", neumann);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }
    //! [Parse command line]

    if (discreteDegree - discreteRegularity > 1)
        gsInfo << "ERROR, does not work for C^1: Choose C^2 \n";

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

    // Neumann
    gsFunctionExpr<> sol1der("-4*pi*(cos(4*pi*y) - 1)*sin(4*pi*x)",
                                         "-4*pi*(cos(4*pi*x) - 1)*sin(4*pi*y)",2);

    //! [Boundary condition]
    gsBoundaryConditions<> bc;
    if (geometry == 1000 || geometry == 1100)
        fd.getId(3, bc); // id=2: boundary conditions
    else
        for (gsMultiPatch<>::const_biterator bit = mp.bBegin(); bit != mp.bEnd(); ++bit)
        {
            bc.addCondition(*bit, condition_type::dirichlet, &ms);
            if (neumann)
                bc.addCondition(*bit, condition_type::neumann, &sol1der);
            else
                bc.addCondition(*bit, condition_type::laplace, &laplace);
        }
    bc.setGeoMap(mp);
    //! [Boundary condition]
    gsInfo<<"Boundary conditions:\n"<< bc <<"\n";
    gsInfo<<"Finished\n";

    gsOptionList optionList;
    fd.getId(100, optionList); // id=100: assembler options
    gsInfo << "OptionList: " << optionList << "\n";
    //! [Read input file]

    //! [Refinement]
    gsMultiBasis<> dbasis(mp, false);//true: poly-splines (not NURBS)

    // Elevate and p-refine the basis to order p + numElevate
    // where p is the highest degree in the bases
    dbasis.setDegree( discreteDegree);

    // h-refine each basis
    if (last)
    {
        for (int r =0; r < numRefine; ++r)
            dbasis.uniformRefine(1, discreteDegree-discreteRegularity);
        numRefine = 0;
    }

    // Assume that the condition holds for each patch TODO
    // Refine once
    if (dbasis.basis(0).numElements() < 4)
        dbasis.uniformRefine(1, discreteDegree-discreteRegularity);


    gsInfo << "Patches: "<< mp.nPatches() <<", degree: "<< dbasis.minCwiseDegree() <<"\n";
#ifdef _OPENMP
    gsInfo<< "Available threads: "<< omp_get_max_threads() <<"\n";
#endif
    //! [Refinement]

    //! [Problem setup]
    gsExprAssembler<> A(1,1);
    //gsInfo<<"Active options:\n"<< A.options() <<"\n";

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
    gsSparseSolver<>::SimplicialLDLT solver;

    gsVector<> l2err(numRefine+1), h1err(numRefine+1), h2err(numRefine+1);
    gsInfo<< "(dot1=approxC1construction, dot2=assembled, dot3=solved, dot4=got_error)\n"
        "\nDoFs: ";
    double setup_time(0), ma_time(0), slv_time(0), err_time(0);
    gsStopwatch timer;
    for (int r=0; r<=numRefine; ++r)
    {
        dbasis.uniformRefine(1,discreteDegree -discreteRegularity);

        //gsDebugVar(dbasis.basis(0));

        gsSparseMatrix<real_t> global2local;

        gsApproxC1Spline<2,real_t> approxC1(mp,dbasis);
        approxC1.options().setSwitch("info",info);
        approxC1.options().setSwitch("plot",plot);

        approxC1.init();
        approxC1.compute();

        global2local = approxC1.getSystem();
        global2local = global2local.transpose();
        //global2local.pruned(1,1e-10);
        gsMultiBasis<> dbasis_temp;
        approxC1.getMultiBasis(dbasis_temp);
        bb2.init(dbasis_temp,global2local);
        // Compute the approx C1 space

        gsInfo<< "." <<std::flush; // Approx C1 construction done

        // Setup the system
        u.setup(bc, dirichlet::l2Projection, -1);

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
        A.assemble(bc.get("Laplace"), (igrad(u, G) * nv(G)) * g_L.tr() );

        // Enforce Neumann conditions to right-hand side
        //A.assembleRhsBc(ilapl(u, G) * (igrad(u_ex) * nv(G)).tr(), bc.neumannSides() );

        ma_time += timer.stop();
        gsInfo<< "." <<std::flush;// Assemblying done

        timer.restart();
        solver.compute( A.matrix() );
        solVector = solver.solve(A.rhs());

        slv_time += timer.stop();
        gsInfo<< "." <<std::flush; // Linear solving done

        timer.restart();
        //linferr[r] = ev.max( f-s ) / ev.max(f);

        l2err[r]= math::sqrt( ev.integral( (u_ex - u_sol).sqNorm() * meas(G) ) ); // / ev.integral(f.sqNorm()*meas(G)) );
        
        h1err[r]= l2err[r] +
            math::sqrt(ev.integral( ( igrad(u_ex) - igrad(u_sol,G) ).sqNorm() * meas(G) )); // /ev.integral( igrad(f).sqNorm()*meas(G) ) );

        h2err[r]= h1err[r] +
                 math::sqrt(ev.integral( ( ihess(u_ex) - ihess(u_sol,G) ).sqNorm() * meas(G) )); // /ev.integral( ihess(f).sqNorm()*meas(G) )



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
        ev.options().setInt   ("plot.npts"    , 1000);
        ev.writeParaview( u_sol   , G, "solution");
        //ev.writeParaview( u_ex    , G, "solution_ex");
        //ev.writeParaview( grad(s), G, "solution_grad");
        //ev.writeParaview( grad(f), G, "solution_ex_grad");
        //ev.writeParaview( (f-s), G, "error_pointwise");
    }
    else
        gsInfo << "Done. No output created, re-run with --plot to get a ParaView "
                  "file containing the solution.\n";
    //! [Export visualization in ParaView]

    //! [Export data to CSV]
    std::string cmdName = "-g1000-p3-r2-l5";
    gsCSVOutput<2, real_t> csvOutput(mp, dbasis, cmdName);
    csvOutput.flags = GEOMETRY | MESH | ERROR;
    csvOutput.saveCSVFile(cmdName);
    //! [Export data to CSV]

    return EXIT_SUCCESS;

}// end main
