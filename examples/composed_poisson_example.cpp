/** @file monitor_poisson_composed_r-adaptivity.cpp

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
#include <gsOptim/gsOptim.h>
#include <gsAssembler/gsAdaptiveParametrization.h>

template <short_t _DIM, class T>
gsTensorBSplineBasis<_DIM,T> integrationBasis(const gsTensorBSplineBasis<_DIM,T> & basis1,
                                           const gsTensorBSplineBasis<_DIM,T> & basis2)
{
    gsTensorBSplineBasis<_DIM,T> ibasis(basis1);
    // Integration basis: parent basis with knots of composition basis inserted, and the degree is the sum of the two degrees (?)
    index_t targetDegree;
    for (size_t d = 0; d!=_DIM; d++)
    {
        // 1. Insert interior knots of composition basis
        for (typename gsKnotVector<T>::uiterator it = std::next(basis2.knots(d).ubegin());
                                                    it!= std::prev(basis2.knots(d).uend());
                                                    ++it)
            {
                if (ibasis.knots(d).has(*it))
                    continue;
                ibasis.insertKnot(*it,d);
            }
        // 2. Increase the degree
        targetDegree = ibasis.degree(d) * basis2.degree(d);
        ibasis.degreeIncrease(targetDegree-ibasis.degree(d),d);
    }
    return ibasis;
}

//! [Include namespace]
int main(int arg, char *argv[])
{
    //! [Parse command line]
    bool plot = false;
    index_t numRefine  = 2;
    index_t numElevate = 0;
    std::string geometry_file;

    gsCmdLine cmd("Tutorial on solving a Poisson problem.");
    cmd.addInt( "e", "elevAnalysis","Number of degree elevation steps to perform for the analysis", numElevate );
    cmd.addInt( "r", "refAnalysis", "Number of Uniform h-refinement loops for the analysis",  numRefine );
    cmd.addSwitch("plot", "Create a ParaView visualization file with the solution", plot);
    cmd.addPlainString( "geometry", "Geometry file", geometry_file );

    try { cmd.getValues(arg,argv); } catch (int rv) { return rv; }
    //! [Parse command line]

    GISMO_ASSERT(!geometry_file.empty(),"Please provide a geometry file.");

    std::string dirname = "r-adaptivity_composed_poisson_e"+util::to_string(numElevate)+"_r"+util::to_string(numRefine);
    gsFileManager::mkdir(dirname);
    dirname += gsFileManager::getNativePathSeparator();

    gsFileData<> fd(geometry_file);
    GISMO_ASSERT(fd.hasAny<gsGeometry<real_t>>(),"The input file must contain a geometry.");

    gsGeometry<>::uPtr geom = fd.getFirst<gsGeometry<>>();

    GISMO_ASSERT((dynamic_cast<gsComposedGeometry<real_t>*>(geom.get())),"The geometry must be a composed geometry.");
    gsComposedGeometry<>::uPtr cgeom = memory::make_unique(static_cast<gsComposedGeometry<real_t>*>(geom.release()));
    gsComposedBasis<>::uPtr    cbasis= memory::make_unique(cgeom->basis().clone().release());

    gsFunction<>::uPtr composition = memory::make_unique(cgeom->composition().clone().release());
    // gsGeometry<>::uPtr geometry    = memory::make_unique(cgeom->geometry().clone().release());
    gsBasis<>::uPtr    basis       = memory::make_unique(cbasis->basis().clone().release());
    gsGeometry<>::uPtr geometry    = basis->makeGeometry(cgeom->coefs());

    // Basis for the analysis
    basis->degreeElevate(numElevate);
    for (index_t i = 0; i < numRefine; i++)
        basis->uniformRefine();

    gsInfo<<"Analysis basis:\n"<<*basis<<"\n";
    GISMO_ASSERT((dynamic_cast<gsGeometry<real_t>*>(composition.get())),"The composition must be a geometry.");
    gsInfo<<"Mapper basis:\n"<<static_cast<gsGeometry<real_t>*>(composition.get())->basis()<<"\n";

    ///////////////////////////////////////////////////////
    // Solve poisson
    ///////////////////////////////////////////////////////

    gsMultiPatch<> mp;
    mp.addPatch(*geometry);
    gsMultiPatch<> cmp;
    cmp.addPatch(*cgeom);

    // Define the basis
    gsMultiBasis<> mb(mp);
    gsMultiBasis<> cmb(cmp);

    gsDebugVar(cmp.patch(0));
    gsDebugVar(cmb.basis(0));

    // Source function:
    // const index_t dimexpr = mp.geoDim();
    const index_t dimexpr = mp.domainDim();
    std::string fstring = "-9*pi^2*(2*x - 1.0)^2*cos(3*pi*(sqrt((x - 0.5)^2 + (y - 0.5)^2) - 1/6))^2*exp(-162*(sqrt((x - 0.5)^2 + (y - 0.5)^2) - 1/3)^2)/(2*((x - 0.5)^2 + (y - 0.5)^2)) + 972*sin(3*pi*(sqrt((x - 0.5)^2 + (y - 0.5)^2) - 1/6))*(sqrt((x - 0.5)^2 + (y - 0.5)^2) - 1/3)*(2*x - 1.0)^2*exp(-162*(sqrt((x - 0.5)^2 + (y - 0.5)^2) - 1/3)^2)*pi*cos(3*pi*(sqrt((x - 0.5)^2 + (y - 0.5)^2) - 1/6))/((x - 0.5)^2 + (y - 0.5)^2) + (3*sin(3*pi*(sqrt((x - 0.5)^2 + (y - 0.5)^2) - 1/6))*exp(-162*(sqrt((x - 0.5)^2 + (y - 0.5)^2) - 1/3)^2)*pi*(2*x - 1.0)^2*cos(3*pi*(sqrt((x - 0.5)^2 + (y - 0.5)^2) - 1/6)))/(2*((x - 0.5)^2 + (y - 0.5)^2)^(3/2)) - 12*sin(3*pi*(sqrt((x - 0.5)^2 + (y - 0.5)^2) - 1/6))*exp(-162*(sqrt((x - 0.5)^2 + (y - 0.5)^2) - 1/3)^2)*pi*cos(3*pi*(sqrt((x - 0.5)^2 + (y - 0.5)^2) - 1/6))/sqrt((x - 0.5)^2 + (y - 0.5)^2) + 9*sin(3*pi*(sqrt((x - 0.5)^2 + (y - 0.5)^2) - 1/6))^2*exp(-162*(sqrt((x - 0.5)^2 + (y - 0.5)^2) - 1/3)^2)*pi^2*(2*x - 1.0)^2/(2*((x - 0.5)^2 + (y - 0.5)^2)) + 81*sin(3*pi*(sqrt((x - 0.5)^2 + (y - 0.5)^2) - 1/6))^2*(2*x - 1.0)^2*exp(-162*(sqrt((x - 0.5)^2 + (y - 0.5)^2) - 1/3)^2)/((x - 0.5)^2 + (y - 0.5)^2) - 81*sin(3*pi*(sqrt((x - 0.5)^2 + (y - 0.5)^2) - 1/6))^2*(sqrt((x - 0.5)^2 + (y - 0.5)^2) - 1/3)*(2*x - 1.0)^2*exp(-162*(sqrt((x - 0.5)^2 + (y - 0.5)^2) - 1/3)^2)/((x - 0.5)^2 + (y - 0.5)^2)^(3/2) + 648*sin(3*pi*(sqrt((x - 0.5)^2 + (y - 0.5)^2) - 1/6))^2*(sqrt((x - 0.5)^2 + (y - 0.5)^2) - 1/3)*exp(-162*(sqrt((x - 0.5)^2 + (y - 0.5)^2) - 1/3)^2)/sqrt((x - 0.5)^2 + (y - 0.5)^2) - 26244*sin(3*pi*(sqrt((x - 0.5)^2 + (y - 0.5)^2) - 1/6))^2*(sqrt((x - 0.5)^2 + (y - 0.5)^2) - 1/3)^2*(2*x - 1.0)^2*exp(-162*(sqrt((x - 0.5)^2 + (y - 0.5)^2) - 1/3)^2)/((x - 0.5)^2 + (y - 0.5)^2) - 9*pi^2*(2*y - 1.0)^2*cos(3*pi*(sqrt((x - 0.5)^2 + (y - 0.5)^2) - 1/6))^2*exp(-162*(sqrt((x - 0.5)^2 + (y - 0.5)^2) - 1/3)^2)/(2*((x - 0.5)^2 + (y - 0.5)^2)) + 972*sin(3*pi*(sqrt((x - 0.5)^2 + (y - 0.5)^2) - 1/6))*(sqrt((x - 0.5)^2 + (y - 0.5)^2) - 1/3)*(2*y - 1.0)^2*exp(-162*(sqrt((x - 0.5)^2 + (y - 0.5)^2) - 1/3)^2)*pi*cos(3*pi*(sqrt((x - 0.5)^2 + (y - 0.5)^2) - 1/6))/((x - 0.5)^2 + (y - 0.5)^2) + (3*sin(3*pi*(sqrt((x - 0.5)^2 + (y - 0.5)^2) - 1/6))*exp(-162*(sqrt((x - 0.5)^2 + (y - 0.5)^2) - 1/3)^2)*pi*(2*y - 1.0)^2*cos(3*pi*(sqrt((x - 0.5)^2 + (y - 0.5)^2) - 1/6)))/(2*((x - 0.5)^2 + (y - 0.5)^2)^(3/2)) + 9*sin(3*pi*(sqrt((x - 0.5)^2 + (y - 0.5)^2) - 1/6))^2*exp(-162*(sqrt((x - 0.5)^2 + (y - 0.5)^2) - 1/3)^2)*pi^2*(2*y - 1.0)^2/(2*((x - 0.5)^2 + (y - 0.5)^2)) + 81*sin(3*pi*(sqrt((x - 0.5)^2 + (y - 0.5)^2) - 1/6))^2*(2*y - 1.0)^2*exp(-162*(sqrt((x - 0.5)^2 + (y - 0.5)^2) - 1/3)^2)/((x - 0.5)^2 + (y - 0.5)^2) - 81*sin(3*pi*(sqrt((x - 0.5)^2 + (y - 0.5)^2) - 1/6))^2*(sqrt((x - 0.5)^2 + (y - 0.5)^2) - 1/3)*(2*y - 1.0)^2*exp(-162*(sqrt((x - 0.5)^2 + (y - 0.5)^2) - 1/3)^2)/((x - 0.5)^2 + (y - 0.5)^2)^(3/2) - 26244*sin(3*pi*(sqrt((x - 0.5)^2 + (y - 0.5)^2) - 1/6))^2*(sqrt((x - 0.5)^2 + (y - 0.5)^2) - 1/3)^2*(2*y - 1.0)^2*exp(-162*(sqrt((x - 0.5)^2 + (y - 0.5)^2) - 1/3)^2)/((x - 0.5)^2 + (y - 0.5)^2)";
    std::string msstring = "sin(3*pi*(sqrt((x - 0.5)^2 + (y - 0.5)^2) - 1/6))^2*exp(-162*(sqrt((x - 0.5)^2 + (y - 0.5)^2) - 1/3)^2)";

    gsFunctionExpr<> f(fstring, dimexpr);
    gsFunctionExpr<> ms(msstring,dimexpr);

    gsDebugVar(f.domainDim());
    gsDebugVar(ms.domainDim());
    gsDebugVar(composition->domainDim());
    gsComposedFunction<real_t> cf(*composition,f);
    gsComposedFunction<real_t> cms(*composition,ms);

    gsBoundaryConditions<> bc;
    bc.addCondition(boundary::side::west ,condition_type::dirichlet,&cms,0,true);
    bc.addCondition(boundary::side::east ,condition_type::dirichlet,&cms,0,true);
    bc.addCondition(boundary::side::south,condition_type::dirichlet,&cms,0,true);
    bc.addCondition(boundary::side::north,condition_type::dirichlet,&cms,0,true);
    bc.setGeoMap(cmp);

    gsMultiBasis<> ib = mb; // integration basis

    // Expression assembler and evaluator
    gsExprAssembler<> A(1,1);

    typedef typename gsExprAssembler<>::geometryMap geometryMap;
    typedef typename gsExprAssembler<>::variable    variable;
    typedef typename gsExprAssembler<>::space       space;
    typedef typename gsExprAssembler<>::solution    solution;

    // Elements used for numerical integration
    A.setIntegrationElements(ib);
    gsExprEvaluator<> ev(A);
    gsMatrix<> solVector;

    A.options().setReal("quA",2.0);
    A.options().setInt("quB",2.0);
    A.options().setSwitch("SameElement",false);

    // Set the geometry map
    geometryMap G = A.getMap(cmp);

    // Set the discretization space
    space u = A.getSpace(cmb);

    // Set the source term
    auto ff = A.getCoeff(cf);

    // Recover manufactured solution
    auto u_ex = ev.getVariable(cms);

    // Solution vector and solution variable
    solution u_sol = A.getSolution(u, solVector);

    //! [Problem setup]

    gsParaviewCollection collection(dirname+"solution",&ev);
    collection.options().setSwitch("plotElements", true);
    collection.options().setInt("plotElements.resolution", 4);
    collection.options().setInt("numPoints", 1000);
    typename gsSparseSolver<>::CGDiagonal solver;
    for (index_t i = 0; i < 5; i++)
    {
        u.setup(bc, dirichlet::l2Projection, 0);

        // Initialize the system
        A.initSystem();
        gsInfo<<"Number of analysis degrees of freedom: "<<A.numDofs()<<"\n";

        // Compute the system matrix and right-hand side
        A.assemble(
            igrad(u, G) * igrad(u, G).tr() * meas(G) //matrix
            ,
            u * ff * meas(G) //rhs vector
            );


        solver.compute( A.matrix() );
        solVector = solver.solve(A.rhs());

        // Convert solution to gsMultiPatch
        gsMultiPatch<> fun;
        u_sol.extract(fun);

        gsDebug<<"error = "<<math::sqrt( ev.integral( (u_ex - u_sol).sqNorm() * meas(G) ) )<<"\n";

        collection.newTimeStep(&mp);
        collection.addField(u_sol,"numerical solution");
        collection.addField(u_ex, "exact solution");
        collection.addField((u_ex-u_sol).norm(), "error");
        collection.saveTimeStep();
    }

    collection.save();

    return 0;
}// end main
