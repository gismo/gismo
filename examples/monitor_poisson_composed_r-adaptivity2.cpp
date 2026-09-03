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

//! [Include namespace]

template <class T>
class gsL2Difference : public gsFunction<T>
{
    typedef memory::shared_ptr<gsL2Difference<T>> Ptr;
    typedef memory::unique_ptr<gsL2Difference<T>> uPtr;

public:
    gsL2Difference(const gsFunction<T> & sol, const gsFunction<T> & ex, T scaling = 1)
    :
    m_sol(sol),
    m_ex(ex),
    m_scaling(scaling)
    {
        GISMO_ASSERT(m_sol.domainDim() == m_ex.domainDim(),"The solution and the exact solution must have the same domain dimension.");
    }

    GISMO_CLONE_FUNCTION(gsL2Difference)

    short_t domainDim() const override {return m_sol.domainDim();}

    void eval_into(const gsMatrix<T> & u, gsMatrix<T> & result) const override
    {
        result.resize(1,u.cols());
        gsMatrix<T> tmp;
        m_sol.eval_into(u,tmp);
        tmp -= m_ex.eval(u);
        result.row(0) = tmp.colwise().norm();
        result *= m_scaling;
    }

    void deriv_into(const gsMatrix<T> & u, gsMatrix<T> & result) const override
    {
        // Result is a 2xN matrix
        // result.col(i) = (d(sol-ex)) / ||sol-ex||
        gsMatrix<T> num, denom;
        m_sol.eval_into(u,num);
        num -= m_ex.eval(u);
        denom = num.colwise().norm();
        m_sol.deriv_into(u,result);
        result -= m_ex.deriv(u);
        result.array().rowwise() *= num.array();
        result.array().rowwise() /= denom.array();
    }

protected:
    const gsFunctionSet<T> & m_sol;
    const gsFunctionSet<T> & m_ex;
    T m_scaling;
};

template <class T>
gsMultiPatch<T> solvePoisson(const gsMultiPatch<T> & mp,
                             const gsMultiBasis<T> & mb,
                             const gsMultiBasis<T> & ib,
                             const gsFunction<T> & f,
                             const gsFunction<T> & ms,
                             const gsBoundaryConditions<T> & bc)
{
        //! [Problem setup]
    gsExprAssembler<T> A(1,1);

    typedef typename gsExprAssembler<T>::geometryMap geometryMap;
    typedef typename gsExprAssembler<T>::variable    variable;
    typedef typename gsExprAssembler<T>::space       space;
    typedef typename gsExprAssembler<T>::solution    solution;

    // Elements used for numerical integration
    A.setIntegrationElements(ib);
    gsExprEvaluator<T> ev(A);
    gsMatrix<T> solVector;

    A.options().setReal("quA",2.0);
    A.options().setInt("quB",2.0);
    A.options().setSwitch("SameElement",false);

    // Set the geometry map
    geometryMap G = A.getMap(mp);
    // geometryMap cG = A.getMap(cmp);

    // Set the discretization space
    space u = A.getSpace(mb);

    // Set the source term
    auto ff = A.getCoeff(f);

    // Recover manufactured solution
    auto u_ex = ev.getVariable(ms);

    // Solution vector and solution variable
    solution u_sol = A.getSolution(u, solVector);

    //! [Problem setup]

    typename gsSparseSolver<T>::CGDiagonal solver;

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
    gsMultiPatch<T> fun;
    u_sol.extract(fun);

    gsDebug<<"error = "<<math::sqrt( ev.integral( (u_ex - u_sol).sqNorm() * meas(G) ) )<<"\n";

    return fun;
}

int main(int arg, char *argv[])
{
    //! [Parse command line]
    bool plot = false;
    index_t numRefine  = 2;
    index_t numRefineD = 2;
    index_t numElevate = 0;
    index_t numElevateD= 0;
    index_t maxIt = 100;
    real_t tol_g = 5e-5;
    real_t eps = 1e-4;
    bool slide = true;
    index_t expression = 0;

    std::string fn("pde/poisson2d_bvp.xml");
    std::string geo;

    gsCmdLine cmd("Tutorial on solving a Poisson problem.");
    cmd.addInt( "e", "elevAnalysis","Number of degree elevation steps to perform for the analysis", numElevate );
    cmd.addInt( "E", "elevDomain","Number of degree elevation steps to perform for the domain", numElevateD );
    cmd.addInt( "r", "refAnalysis", "Number of Uniform h-refinement loops for the analysis",  numRefine );
    cmd.addInt( "R", "refDomain", "Number of Uniform h-refinement loops for the domain",  numRefineD );
    cmd.addSwitch("plot", "Create a ParaView visualization file with the solution", plot);
    cmd.addReal("g", "tolG", "relative tol", tol_g);
    cmd.addInt( "i", "maxIt", "max num iterations",  maxIt );
    cmd.addReal("", "eps", "eps",  eps );
    cmd.addSwitch("noslide", "Do not slide the boundaries",  slide );
    cmd.addInt( "f", "expr", "Problem to be solved: 0 default peak in the center of the domain; 1 peak in the bottom left corner of the domain.",  expression );
    cmd.addString( "f", "file", "Input XML file", fn );
    cmd.addString( "g", "geometry", "Geometry file", geo );

    try { cmd.getValues(arg,argv); } catch (int rv) { return rv; }
    //! [Parse command line]

    std::string dirname = "r-adaptivity_composed_poisson_e"+util::to_string(numElevate)+"_r"+util::to_string(numRefine)+"_E"+util::to_string(numElevateD)+"_R"+util::to_string(numRefineD);
    gsFileManager::mkdir(dirname);
    dirname += gsFileManager::getNativePathSeparator();

    gsFileData<> fd(fn);
    gsInfo << "Loaded file "<< fd.lastPath() <<"\n";

    gsFunctionExpr<> f;
    fd.getId(1, f); // id=1: source function
    gsInfo<<"Source function "<< f << "\n";

    gsBoundaryConditions<> bc;
    fd.getId(2, bc); // id=2: boundary conditions
    bc.setGeoMap(mp);
    gsInfo<<"Boundary conditions:\n"<< bc <<"\n";

    gsFunctionExpr<> ms;
    fd.getId(3, ms); // id=3: reference solution

    gsOptionList Aopt;
    fd.getId(4, Aopt); // id=4: assembler options

    if (!geo.empty())
        fd.read(geo);

    gsMultiPatch<> mp, cmp;
    gsComposedGeometry<real_t>::uPtr cspline;
    gsComposedBasis<real_t>::uPtr cbasis;
    gsFunction<real_t>::uPtr composition;
    gsGeometry<real_t>::uPtr spline;
    gsBasis<real_t>::uPtr basis;
    if (fd.hasAny<gsMultiPatch<real_t>>())
    {
        fd.getId(0, cmp); // id=0: Multipatch domain
        GISMO_ASSERT((dynamic_cast<gsComposedGeometry<T>*>(&cmp.patch(0))),"The geometry must be a composed geometry.");
        cspline = memory::make_unique(static_cast<gsComposedGeometry<real_t>*>(&cmp.patch(0)));
    }
    else if (fd.hasAny<gsGeometry<real_t>>())
    {
        gsGeometry<real_t>::uPtr patch = fd.getId<gsGeometry<real_t>>(0);
        GISMO_ASSERT((dynamic_cast<gsComposedGeometry<T>*>(patch.get())),"The geometry must be a composed geometry.");
        cspline = memory::make_unique(static_cast<gsComposedGeometry<real_t>*>(patch.release()));
        cmp.addPatch(cspline);
    }
    else
        GISMO_ERROR("The input file must contain a geometry.");

    cbasis= memory::make_unique(cgeom->basis().clone().release());
    basis = memory::make_unique(cbasis->basis().clone().release());
    spline= basis->makeGeometry(cgeom->coefs());
    mp.addPatch(spline);

    // Basis for the analysis
    basis->degreeElevate(numElevate);
    for (index_t i = 0; i < numRefine; i++)
        basis->uniformRefine();

    gsInfo<<"Analysis basis:\n"<<*basis<<"\n";


    ///////////////////////////////////////////////////////
    // Solve poisson
    ///////////////////////////////////////////////////////

    // // Define the domain (which is a composed geometry)
    // gsComposedGeometry<real_t> cspline(domain,tbspline);

    // gsMultiPatch<> mp;
    // mp.addPatch(tbspline);
    // gsMultiPatch<> cmp;
    // cmp.addPatch(cspline);

    // Define the basis
    gsMultiBasis<> mb(mp);
    gsMultiBasis<> cmb(cmp);

    // Source function:
    gsComposedFunction<real_t> cf(cmp.patch(0),f);
    gsComposedFunction<real_t> cms(cmp.patch(0),ms);

    gsBoundaryConditions<> bc;
    bc.addCondition(boundary::side::west ,condition_type::dirichlet,&cms,0,true);
    bc.addCondition(boundary::side::east ,condition_type::dirichlet,&cms,0,true);
    bc.addCondition(boundary::side::south,condition_type::dirichlet,&cms,0,true);
    bc.addCondition(boundary::side::north,condition_type::dirichlet,&cms,0,true);
    bc.setGeoMap(cmp);

    gsParaviewCollection errcol(dirname+"errors");
    gsParaviewCollection solcol(dirname+"solution");

    gsMultiPatch<> sol;
    real_t totalError, maxError;

    index_t i = 0;
    /*
        Functions/FunctionSets:
        mp     := the original geometry
        cmp    := the composed geometry
        mb     := the original basis
        cmb    := the composed basis
        f      := the source function
        cf     := the composed source function
        ms     := the manufactured solution
        cms    := the composed manufactured solution
        sol    := the numerical solution defined using a composed basis
        parsol := the numerical solution defined using the original basis and the coefficients of sol
        err    := the L2 error between the numerical solution and the manufactured solution

        Expressions:
        G      := the geometry map of the original geometry
        cG     := the geometry map of the composed geometry
        u_sol  := the numerical solution defined using the composed basis
        u_ex   := the exact solution defined on the composed parametric domain
        */

    // Write the mesh and control net to the paraview collections
    writeSingleCompMesh(cmp.basis(0),cmp.patch(0),dirname+"domain_mesh_"+util::to_string(i));
    errcol.addTimestep("domain_mesh_"+util::to_string(i),i,".vtp");
    solcol.addTimestep("domain_mesh_"+util::to_string(i),i,".vtp");
    writeSingleControlNet(domain.domain(),dirname+"domain_cnet_"+util::to_string(i));
    errcol.addTimestep("domain_cnet_"+util::to_string(i),i,".vtp");
    solcol.addTimestep("domain_cnet_"+util::to_string(i),i,".vtp");

    // Solve the poisson equation
    sol = solvePoisson(cmp,cmb,mb,cf,cms,bc);

    // Construct parsol and the error in the original parametric domain
    gsMultiPatch<> parsol;
    parsol.addPatch(mb.basis(0).makeGeometry(sol.patch(0).coefs()));
    gsL2Difference<real_t> err(ms,parsol.patch(0));

    // Prepare the error computation
    gsExprEvaluator<real_t> ev;
    ev.setIntegrationElements(mb);
    // gsWarn<<"Construct an integration basis instead of hard-coding the quadrature rule\n";
    ev.options().setReal("quA",4.0);
    ev.options().setInt("quB",2.0);
    auto G     = ev.getMap(mp);
    auto cG    = ev.getMap(cmp);
    auto u_sol = ev.getVariable(sol);
    auto u_psol = ev.getVariable(parsol);
    auto u_ex  = ev.getVariable(cms);
    auto E     = ev.getVariable(err);

    // Evaluate the error(s) and maximum error
    gsInfo<<"ev.integral((u_sol-u_ex).norm()*meas(G)) = "<<ev.integral((u_sol-u_ex).norm()*meas(cG))<<"\n";
    gsInfo<<"ev.integral(E*meas(G))                   = "<<ev.integral(E*meas(G))<<"\n";
    totalError = ev.integral((u_sol-u_ex).norm()*meas(cG));
    gsInfo<<"total error after iteration "<<i<<": "<<totalError<<"\n";
    maxError = ev.max((u_sol-u_ex).norm());

    // Write the error(s) and solution to paraview
    writeSinglePatchField(cmp.patch(0),err,true,dirname+"error_"+util::to_string(i),1000);
    errcol.addTimestep("error_"+util::to_string(i),i,".vts");
    writeSinglePatchField(cmp.patch(0),parsol.patch(0),true,dirname+"sol_"+util::to_string(i),1000);
    solcol.addTimestep("sol_"+util::to_string(i),i,".vts");

    ev.writeParaview((u_sol),cG,dirname+"u_sol");
    ev.writeParaview((u_psol),cG,dirname+"u_psol");
    ev.writeParaview((u_ex),cG,dirname+"u_ex");
    // ev.writeParaview((u_ex-u_sol).norm(),cG,dirname+"u_err");

    // Solve r-adaptivity
    gsInfo<<"R-Adaptivity iteration "<<i<<"\n";
    gsL2Difference<real_t> err_scaled(ms,parsol.patch(0),10./maxError);
    // gsAdaptiveParametrization<real_t,MonitorMode::ValueBased> relocator(domain,mp.patch(0),err_scaled,mb.basis(0),optimizer,true);
    // gsAdaptiveParametrization<real_t,MonitorMode::ValueBased> relocator(domain,mp.patch(0),parsol.patch(0),mb.basis(0),optimizer,true);
    // gsAdaptiveParametrization<real_t,MonitorMode::GradientBased> relocator(domain,mp.patch(0),parsol.patch(0),mb.basis(0),optimizer,true);
    // relocator.solve();

    gsVector<> controls = domain.getControls();
    controls *= 0.5;
    domain.setControls(controls);

    gsInfo<<"Domain Basis:\n"<<domain.domain().basis()<<"\n";
    gsInfo<<"Domain Coefficients:\n"<<domain.domain().coefs()<<"\n";

    // Solve the poisson equation
    sol = solvePoisson(cmp,cmb,mb,cf,cms,bc);
    parsol.patch(0).coefs() = sol.patch(0).coefs();

    ev.writeParaview((u_ex),cG,dirname+"u_ex____2");


    ++i;
    // Write the mesh and control net to the paraview collections
    writeSingleCompMesh(cmp.basis(0),cmp.patch(0),dirname+"domain_mesh_"+util::to_string(i));
    errcol.addTimestep("domain_mesh_"+util::to_string(i),i,".vtp");
    solcol.addTimestep("domain_mesh_"+util::to_string(i),i,".vtp");
    writeSingleControlNet(domain.domain(),dirname+"domain_cnet_"+util::to_string(i));
    errcol.addTimestep("domain_cnet_"+util::to_string(i),i,".vtp");
    solcol.addTimestep("domain_cnet_"+util::to_string(i),i,".vtp");
    // Write the error(s) and solution to paraview
    writeSinglePatchField(cmp.patch(0),err,true,dirname+"error_"+util::to_string(i),1000);
    errcol.addTimestep("error_"+util::to_string(i),i,".vts");
    writeSinglePatchField(cmp.patch(0),parsol.patch(0),true,dirname+"sol_"+util::to_string(i),1000);
    solcol.addTimestep("sol_"+util::to_string(i),i,".vts");

    errcol.save();
    solcol.save();
    gsWriteParaview(domain.domain(),dirname+"domain",1000,true,true);




    return 0;
}// end main
