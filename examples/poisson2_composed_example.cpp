/** @file composed_domain_test.cpp

    @brief Tutorial on how to use expression assembler to solve the Poisson equation
           using composite maps

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst, S. Imperatore, A. Mantzaflaris
*/

//! [Include namespace]
#include <gismo.h>
#include <gsNurbs/gsSquareDomain.h>
#include <gsAssembler/gsAdaptiveParametrization.h>
#include <gsNurbs/gsMobiusDomain.h>

using namespace gismo;
//! [Include namespace]

int main(int argc, char *argv[])
{
    //! [Parse command line]
    bool plot = false;
    index_t numRefine  = 5;
    index_t numElevate = 0;
    index_t numRefineC = 1;

    bool last{false}, export_b64{false}, write{false};
    std::string fn("pde/poisson2d_bvp.xml");

    gsCmdLine cmd("Tutorial on solving a Poisson problem.");
    cmd.addInt( "e", "degreeElevation",
                "Number of degree elevation steps to perform before solving (0: equalize degree in all directions)", numElevate );
    cmd.addInt( "r", "uniformRefine", "Number of Uniform h-refinement loops",  numRefine );
    cmd.addString( "f", "file", "Input XML file", fn );
    cmd.addSwitch("last", "Solve solely for the last level of h-refinement",
                  last);
    cmd.addSwitch("plot", "Create a ParaView visualization file with the solution", plot);
    cmd.addSwitch("write", "Write the convergence data to a file", write);
    cmd.addSwitch("binary", "Use B64 encoding for Paraview", export_b64);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }
    //! [Parse command line]

    // std::ofstream file_out;
    // file_out.open("composedDomain_results_r"+util::to_string(numRefine)+"e"+util::to_string(numElevate)+"R"+util::to_string(numRefineC)+"E"+util::to_string(numElevateC)+".csv");
    // file_out << "problem, deg, ref, dofs, L2err\n";

    //! [Read input file]

    gsFileData<> fd(fn);
    gsInfo << "Loaded file "<< fd.lastPath() <<"\n";

    gsMultiPatch<> mp, cmp;
    fd.getId(0, cmp); // id=0: Multipatch domain
    GISMO_ENSURE(cmp.nPatches()==1, "This example works only with one patch");
    GISMO_ENSURE(nullptr!=dynamic_cast<gsComposedGeometry<real_t>*>(&cmp.patch(0)),"This example works only with composed geometries");

    gsComposedGeometry<real_t> & cgeom = static_cast<gsComposedGeometry<real_t> &>(cmp.patch(0));
    gsComposedBasis<real_t>    & cbasis = static_cast<gsComposedBasis<real_t> &>(cgeom.basis());
    gsFunction<>::uPtr domain = cbasis.composition().clone();
    gsBasis<>::uPtr    basis  = cbasis.basis().clone();
    gsGeometry<>::uPtr geom   = basis->makeGeometry(cgeom.coefs());

    // gsComposedGeometry<real_t> & cgeom = static_cast<gsComposedGeometry<real_t> &>(cmp.patch(0));
    // gsFunction<> & domain = cgeom.composition();
    // gsComposedBasis<real_t> & cbasis = static_cast<gsComposedBasis<real_t> &>(cgeom.basis());
    // gsBasis<>               & basis  = cbasis.basis();
    // mp.addPatch(basis.makeGeometry(cgeom.coefs()));

    gsFunctionExpr<> f;
    fd.getId(1, f); // id=1: source function
    gsInfo<<"Source function "<< f << "\n";

    gsBoundaryConditions<> bc;
    gsBoundaryConditions<> cbc;
    fd.getId(2, bc); // id=2: boundary conditions
    bc.setGeoMap(*geom);
    gsInfo<<"Boundary conditions:\n"<< bc <<"\n";

    std::vector<gsFunction<real_t> *> functions(bc.size());
    std::vector<gsComposedFunction<real_t>> cfunctions(bc.size());
    index_t fIndex=0;
    for (typename gsBoundaryConditions<real_t>::const_iterator it = bc.begin("Dirichlet"); it!=bc.end("Dirichlet"); ++it, fIndex++)
    {
        functions[fIndex] = dynamic_cast<gsFunction<real_t> *>(it->function().get());
        cfunctions[fIndex] = gsComposedFunction<real_t>(&cgeom,functions[fIndex]);
        cbc.addCondition(it->side(),it->type(),&cfunctions[fIndex],it->unknown(),true,it->unkComponent());
    }
    for (typename gsBoundaryConditions<real_t>::const_iterator it = bc.begin("Neumann"); it!=bc.end("Neumann"); ++it, fIndex++)
    {
        functions[fIndex] = dynamic_cast<gsFunction<real_t> *>(it->function().get());
        cfunctions[fIndex] = gsComposedFunction<real_t>(&cgeom,functions[fIndex]);
        cbc.addCondition(it->side(),it->type(),&cfunctions[fIndex],it->unknown(),it->parametric(),it->unkComponent());
    }
    for (typename gsBoundaryConditions<real_t>::citerator it = bc.cornerBegin(); it!=bc.cornerEnd(); ++it)
        cbc.addCornerValue(it->corner,it->value,it->patch,it->unknown,it->component);

    cbc.setGeoMap(cgeom);

    gsFunctionExpr<> ms;
    fd.getId(3, ms); // id=3: reference solution

    gsOptionList Aopt;
    fd.getId(4, Aopt); // id=4: assembler options

    //! [Read input file]



    //! [Refinement]
    if (numElevate!=0)
        basis->degreeElevate(numElevate);
        // mp.degreeElevate(numElevate);

    // h-refine each basis
    if (last)
    {
        for (int r =0; r < numRefine; ++r)
            basis->uniformRefine();
        numRefine = 0;
    }
#ifdef _OPENMP
    gsInfo<< "Available threads: "<< omp_get_max_threads() <<"\n";
#endif
    //! [Refinement]

    //! [Problem setup]
    gsExprAssembler<> A(1,1);
    A.setOptions(Aopt);
    A.options().setSwitch("SameElement",false); // HMV: needed
    gsExprEvaluator<> ev(A);
    // gsExprEvaluator takes the assembler's DATA but its own default options
    // (gsExprEvaluator.h:74-76), and SameElement is not among them, so
    // askSwitch defaults it to TRUE. For a composed map the quadrature points
    // of one integration cell do NOT share an analysis element, so the error
    // evaluation must switch it off explicitly -- otherwise the reported error
    // is corrupted for any non-identity sigma (identity sigma is unaffected,
    // which is why a null gate alone does not catch this).
    ev.options().addSwitch("SameElement","",false);

    gsInfo<<"Active options:\n"<< A.options() <<"\n";

    typedef gsExprAssembler<>::geometryMap geometryMap;
    typedef gsExprAssembler<>::variable    variable;
    typedef gsExprAssembler<>::space       space;
    typedef gsExprAssembler<>::solution    solution;

    //! [Problem setup]

    //! [Solver loop]
#ifdef GISMO_WITH_PARDISO
    gsSparseSolver<>::PardisoLDLT solver;
#else
    gsSparseSolver<>::CGDiagonal solver;
#endif

    gsParaviewCollection collection("ParaviewOutput/solution", &ev);
    collection.options().setSwitch("plotElements", true);
    collection.options().setSwitch("base64", export_b64);
    collection.options().setInt("plotElements.resolution", 16);

    gsVector<> l2err(numRefine+1), h1err(numRefine+1);
    gsInfo<< "(dot1=assembled, dot2=solved, dot3=got_error)\n"
        "\nDoFs: ";
    for (int r=0; r<=numRefine; ++r)
    {
        basis->uniformRefine();

        // Integration basis: the UNION of the analysis knots and sigma's knots
        // at the product degree. Integrating on the analysis mesh alone leaves
        // sigma's knot lines inside elements and loses 1-2 orders of accuracy.
        const gsGeometry<real_t> * sigmaGeo =
            dynamic_cast<const gsGeometry<real_t>*>(domain.get());
        GISMO_ENSURE(nullptr!=sigmaGeo,"The composition must be a geometry.");
        const gsTensorBSplineBasis<2,real_t> & sigmaBasis =
            static_cast<const gsTensorBSplineBasis<2,real_t>&>(sigmaGeo->basis());

        const gsTensorBSplineBasis<2,real_t> * anaTB =
            dynamic_cast<const gsTensorBSplineBasis<2,real_t>*>(basis.get());
        const gsTensorBSplineBasis<2,real_t> anaBasis = (nullptr!=anaTB) ? *anaTB
            : static_cast<const gsTensorNurbsBasis<2,real_t>&>(*basis).source();

        gsMultiBasis<> ibasis(
            gsAdaptiveParametrization<real_t,MonitorMode::ValueBased>::
                makeIntegrationBasis<2>(anaBasis,sigmaBasis));
        // Analysis basis
        gsComposedBasis<real_t> cbasis(*domain,*basis);
        gsMultiBasis<> cmb(cbasis);

        // Composed function
        gsComposedFunction<real_t> cf(cgeom,f);

        // Composed ms
        gsComposedFunction<real_t> cms(cgeom,ms);

        // gsBoundaryConditions<> cbc;
        // cbc.addCondition(boundary::side::west ,condition_type::dirichlet,&cms,0,true);
        // cbc.addCondition(boundary::side::east ,condition_type::dirichlet,&cms,0,true);
        // cbc.addCondition(boundary::side::south,condition_type::dirichlet,&cms,0,true);
        // cbc.addCondition(boundary::side::north,condition_type::dirichlet,&cms,0,true);
        // cbc.setGeoMap(cgeom);

        // Elements used for numerical integration
        A.setIntegrationElements(ibasis);

        // Set the geometry map
        geometryMap G = A.getMap(cgeom);

        // Set the discretization space
        space u = A.getSpace(cmb);

        // Set the source term
        auto ff = A.getCoeff(cf);

        // Recover manufactured solution
        // Bind the PHYSICAL exact solution to the composed map, exactly as
        // poisson2_example.cpp:118 does. Using the pre-composed cms without a
        // map makes igrad(u_ex) the parametric gradient and corrupts the H1
        // error (the L2 value is unaffected).
        auto u_ex = ev.getVariable(ms, G);

        // Solution vector and solution variable
        gsMatrix<> solVector;
        solution u_sol = A.getSolution(u, solVector);

        // Setup the space \a u with strongly imposed Dirichlet part
        //u.setup(bc, dirichlet::interpolation, 0);
        u.setup(cbc, dirichlet::l2Projection, 0);

        // Initialize the system
        A.initSystem();
        // Compute sparsity patter: this is done automatically - but
        // is needed if assemble(.) is called twice
        A.computePattern( igrad(u) * igrad(u).tr() );
        gsInfo<< A.numDofs() <<std::flush;

        // Compute the system matrix and right-hand side
        A.assemble(
            igrad(u, G) * igrad(u, G).tr() * meas(G) //matrix
            ,
            u * ff * meas(G) //rhs vector
            );

        // Compute the Neumann terms defined on physical space
        auto g_N = A.getBdrFunction(G);
        A.assembleBdr(bc.get("Neumann"), u * g_N.tr() * nv(G) );
        gsInfo<< "." <<std::flush;// Assemblying done

        solver.compute( A.matrix() );
        solVector = solver.solve(A.rhs());
        gsInfo<< "." <<std::flush; // Linear solving done

        // Compute the L2 and H1 error, based on manufactured solution
        l2err[r]= math::sqrt( ev.integral( (u_ex - u_sol).sqNorm() * meas(G) ) );
        h1err[r]= l2err[r] +
            math::sqrt(ev.integral( ( igrad(u_ex) - igrad(u_sol,G) ).sqNorm() * meas(G) ));
        gsInfo<< ". " <<std::flush; // Error computations done

        if (plot)
        {
            collection.newTimeStep(&cmp);
            collection.addField(u_sol,"numerical solution");
            collection.addField(u_ex, "exact solution");
            collection.addField((u_ex-u_sol).sqNorm(), "L2-error");
            collection.addField((u_ex-u_sol).sqNorm() + ( igrad(u_ex) - igrad(u_sol) ).sqNorm(), "H1-error");
            collection.saveTimeStep();
        }

    } //for loop
    //! [Solver loop]

    //! [Error and convergence rates]
    gsInfo<< "\nL2 error: "<<std::scientific<<std::setprecision(3)<<l2err.transpose()<<"\n";
    gsInfo<< "H1 error: "<<std::scientific<<h1err.transpose()<<"\n";

    if (!last && numRefine>0)
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
        collection.save();
    }
    else
        gsInfo << "Done. No output created, re-run with --plot to get a ParaView "
                  "file containing the solution.\n";
    //! [Export visualization in ParaView]

    return EXIT_SUCCESS;

}// end main
