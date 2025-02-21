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
    bool write = false;
    bool last = false;
    bool rational = false;
    index_t numRefine  = 2;
    index_t numElevate = 0;
    std::string geometry_file;
    std::string output;

    gsCmdLine cmd("Tutorial on solving a Poisson problem.");
    cmd.addInt( "e", "elevAnalysis","Number of degree elevation steps to perform for the analysis", numElevate );
    cmd.addInt( "r", "refAnalysis", "Number of Uniform h-refinement loops for the analysis",  numRefine );
    cmd.addSwitch("plot", "Create a ParaView visualization file with the solution", plot);
    cmd.addSwitch("write", "Write the convergence data to a file", write);
    cmd.addPlainString( "geometry", "Geometry file", geometry_file );
    cmd.addString("o", "output", "Output directory", output);
    cmd.addSwitch("last", "Solve solely for the last level of h-refinement", last);
    cmd.addSwitch("rational", "Use rational NURBS in gsMultiBasis", rational);
    try { cmd.getValues(arg,argv); } catch (int rv) { return rv; }
    //! [Parse command line]

    GISMO_ASSERT(!geometry_file.empty(),"Please provide a geometry file.");

    std::string dirname = (output.empty()) ? "composed_poisson_e"+util::to_string(numElevate)+"_r"+util::to_string(numRefine) : output;
    if (plot || write)
        gsFileManager::mkdir(dirname);

    dirname += gsFileManager::getNativePathSeparator();

    gsFileData<> fd(geometry_file);
    GISMO_ASSERT(fd.hasAny<gsGeometry<real_t>>(),"The input file must contain a geometry.");

    gsGeometry<>::uPtr geom = fd.getFirst<gsGeometry<>>();

    // For NURBS, we want to activate the NoRational flag of gsMultiBasis
    // Therefore, we need to extract the geometry and the tensor basis
    // The following naming is used:
    // cgeometry: composed geometry (NURBS)
    // composition: the composition of the NURBS
    // composition_basis: the basis of the composition (possibly NURBS)
    // geometry: non-composed geometry (NURBS)
    // cbasis: composed basis (NOT NURBS)
    // basis: non-composed basis (NOT NURBS)
    GISMO_ASSERT((dynamic_cast<gsComposedGeometry<real_t>*>(geom.get())),"The geometry must be a composed geometry.");
    gsComposedGeometry<>::uPtr cgeometry = memory::make_unique(static_cast<gsComposedGeometry<real_t>*>(geom.release()));
    GISMO_ASSERT((dynamic_cast<gsGeometry<real_t>*>(&cgeometry->composition())),"The composition must be a geometry.");
    gsGeometry<>::uPtr composition = memory::make_unique(static_cast<gsGeometry<real_t>*>(&cgeometry->composition())->clone().release());
    gsGeometry<>::uPtr geometry    = cgeometry->basis().basis().makeGeometry(cgeometry->coefs()); // always NURBS if cgeometry is NURBS
    GISMO_ASSERT((dynamic_cast<gsTensorBSplineBasis<2>*>(&composition->basis())),"The composition must be a tensor basis.");
    gsTensorBSplineBasis<2> * composition_basis = static_cast<gsTensorBSplineBasis<2>*>(&composition->basis());

    // Obtain the non-composed basis
    gsBasis<>::uPtr basis;
    if (rational)
        basis = memory::make_unique(cgeometry->basis().basis().clone().release());
    else
        basis = memory::make_unique(cgeometry->basis().basis().source().clone().release());
    // Obtain the composed basis
    gsComposedBasis<> cbasis(*composition,*basis);

    gsInfo<<"Mapper basis:\n"<<*composition_basis<<"\n";
    gsInfo<<"Analysis basis:\n"<<cbasis<<"\n";

    ///////////////////////////////////////////////////////
    // Solve poisson
    ///////////////////////////////////////////////////////

    gsMultiPatch<> mp;
    mp.addPatch(*geometry);
    gsMultiPatch<> cmp;
    cmp.addPatch(*cgeometry);

    // Define the basis
    gsMultiBasis<> mb(*basis,true); // Takes NURBS if requested
    gsMultiBasis<> cmb(cbasis,true); // Takes NURBS if requested

    // Basis for the analysis
    cmb.degreeElevate(numElevate);
    if (last)
    {
        for (index_t i = 0; i < numRefine; i++)
            cmb.uniformRefine();
        numRefine = 0;
    }

    gsDebugVar(cmp.patch(0));
    gsDebugVar(cmb.basis(0));

    // Source function:
    const index_t dimexpr = mp.geoDim();
    std::string fstring = "u:=atan2(y,x);v:=atan2(sqrt(x^2+y^2),z);w:=sqrt(x^2+y^2+z^2);t:=100;-(-2*cos(v)*t*(v - 1/2)*exp(-t*(v - 1/2)^2) - 2*sin(v)*t*exp(-t*(v - 1/2)^2) + 4*sin(v)*t^2*(v - 1/2)^2*exp(-t*(v - 1/2)^2))/(w^2*sin(v))";
    std::string msstring = "u:=atan2(y,x);v:=atan2(sqrt(x^2+y^2),z);w:=sqrt(x^2+y^2+z^2);t:=100;exp(-t*(v - 1/2)^2)";

    gsFunctionExpr<> f(fstring, dimexpr);
    gsFunctionExpr<> ms(msstring,dimexpr);

    gsComposedFunction<real_t> cf(*cgeometry,f);
    gsComposedFunction<real_t> cms(*cgeometry,ms);

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

    gsVector<> l2err(numRefine+1), h1err(numRefine+1), dofs(numRefine+1), bsize(numRefine+1);
    for (index_t i = 0; i <= numRefine; i++)
    {
        cmb.uniformRefine();
        gsComposedBasis<> & cmb_basis = dynamic_cast<gsComposedBasis<> &>(cmb.basis(0));
        if (const gsTensorBSplineBasis<2,real_t> * tbasis = dynamic_cast<const gsTensorBSplineBasis<2,real_t> *>(&cmb_basis.basis()))
            ib = integrationBasis(*tbasis,*composition_basis);
        else if (const gsTensorNurbsBasis<2,real_t> * nbasis = dynamic_cast<const gsTensorNurbsBasis<2,real_t> *>(&cmb_basis.basis()))
            ib = integrationBasis(nbasis->source(),*composition_basis);
        else
            GISMO_ERROR("The basis must be either a tensor or a NURBS basis.");

        A.setIntegrationElements(ib);
        ev.setIntegrationElements(ib);
        u.setup(bc, dirichlet::l2Projection, 0);

        // Initialize the system
        A.initSystem();
        gsInfo<<"Number of analysis degrees of freedom: "<<A.numDofs()<<"\n";

        // Compute the system matrix and right-hand side
        A.computePattern( igrad(u) * igrad(u).tr() );
        A.assemble(
            igrad(u, G) * igrad(u, G).tr() * meas(G) //matrix
            ,
            u * ff * meas(G) //rhs vector
            );

        solver.compute( A.matrix() );
        solVector = solver.solve(A.rhs());

        // Compute the L2 and H1 error, based on manufactured solution
        l2err[i]= math::sqrt( ev.integral( (u_ex - u_sol).sqNorm() * meas(G) ) );
        h1err[i]= l2err[i] +
            math::sqrt(ev.integral( ( igrad(u_ex, G) - igrad(u_sol,G) ).sqNorm() * meas(G) ));
        dofs [i]= A.numDofs();
        bsize[i]= cmb_basis.size();

        if (plot)
        {
            collection.newTimeStep(&cmp);
            collection.addField(u_sol,"numerical solution");
            collection.addField(u_ex, "exact solution");
            collection.addField((u_ex-u_sol).norm(), "error");
            collection.saveTimeStep();
        }
    }

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

    if (write)
    {
        std::ofstream file;
        file.open(dirname+"convergence.txt");
        file<<"numDofs, basisSize, L2 error, H1 error\n";
        for (index_t i = 0; i <= numRefine; i++)
            file<<dofs[i]<<", "<<bsize[i]<<", "<<l2err[i]<<", "<<h1err[i]<<"\n";
        file.close();
    }

    if (plot)
        collection.save();

    return 0;
}// end main
