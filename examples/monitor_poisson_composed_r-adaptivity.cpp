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

    // Set the geometry map
    geometryMap G = A.getMap(mp);
    // geometryMap cG = A.getMap(cmp);

    // Set the discretization space
    space u = A.getSpace(mb);

    // Set the source term
    auto ff = A.getCoeff(f);

    // // Recover manufactured solution
    // auto u_ex = ev.getVariable(ms, G);

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

    try { cmd.getValues(arg,argv); } catch (int rv) { return rv; }
    //! [Parse command line]

    std::string dirname = "r-adaptivity_composed_poisson_e"+util::to_string(numElevate)+"_r"+util::to_string(numRefine)+"_E"+util::to_string(numElevateD)+"_R"+util::to_string(numRefineD);
    gsFileManager::mkdir(dirname);
    dirname += gsFileManager::getNativePathSeparator();

    // Basis for the analysis
    gsKnotVector<> kv({0,0,1,1},1);
    gsTensorBSplineBasis<2> tbasis(kv,kv);

    tbasis.degreeElevate(numElevate);
    for (index_t i = 0; i < numRefine; i++)
        tbasis.uniformRefine();

    gsInfo<<"Analysis basis:\n"<<tbasis<<"\n";

    gsMatrix<> coefs = tbasis.anchors().transpose();
    gsTensorBSpline<2> tbspline(tbasis,coefs);

    // Basis for the square domain
    gsKnotVector<> kv2({0,0,1,1},1);
    gsTensorBSplineBasis<2> dbasis(kv2,kv2);
    dbasis.degreeElevate(numElevateD);
    for (index_t i = 0; i < numRefineD; i++)
        dbasis.uniformRefine();

    gsInfo<<"Mapper basis:\n"<<dbasis<<"\n";

    gsSquareDomain<2,real_t> domain(dbasis);
    domain.options().addSwitch("Slide","",slide);
    domain.applyOptions();

    ///////////////////////////////////////////////////////
    // Solve poisson
    ///////////////////////////////////////////////////////

    // Define the domain (which is a composed geometry)
    gsComposedGeometry<real_t> cspline(domain,tbspline);

    gsMultiPatch<> mp;
    mp.addPatch(tbspline);
    gsMultiPatch<> cmp;
    cmp.addPatch(cspline);

    // Define the basis
    gsMultiBasis<> mb(mp);
    gsMultiBasis<> cmb(cmp);

    // Source function:
    index_t dimexpr = mp.geoDim();
    std::string fstring = "";
    std::string msstring = "";

    switch (expression) {
    case 1:
        // bottom left corner
        fstring = "((tanh(20*(x^2 + y^2)^(1/2) - 5)^2 - 1)*(20*x^2 + 20*y^2)*(40*tanh(20*(x^2 + y^2)^(1/2) - 5)*(x^2 + y^2)^(1/2) - 1))/(x^2 + y^2)^(3/2)";
        msstring = "tanh((0.25-sqrt(x^2+y^2))/0.05)+1";
        break;
    default:
        // center of the domain
        fstring = "((tanh(20*((x-0.5)^2 + (y-0.5)^2)^(1/2) - 5)^2 - 1)*(20*(x-0.5)^2 + 20*(y-0.5)^2)*(40*tanh(20*((x-0.5)^2 + (y-0.5)^2)^(1/2) - 5)*((x-0.5)^2 + (y-0.5)^2)^(1/2) - 1))/((x-0.5)^2 + (y-0.5)^2)^(3/2)";
        msstring = "tanh((0.25-sqrt((x-0.5)^2+(y-0.5)^2))/0.05)+1";
        break;
    }

    gsFunctionExpr<> f(fstring, dimexpr);
    gsFunctionExpr<> ms(msstring,dimexpr);
    gsComposedFunction<real_t> cf(domain,f);
    gsComposedFunction<real_t> cms(domain,ms);

    gsBoundaryConditions<> bc;
    bc.addCondition(boundary::side::west ,condition_type::dirichlet,&cms,0,true);
    bc.addCondition(boundary::side::east ,condition_type::dirichlet,&cms,0,true);
    bc.addCondition(boundary::side::south,condition_type::dirichlet,&cms,0,true);
    bc.addCondition(boundary::side::north,condition_type::dirichlet,&cms,0,true);
    bc.setGeoMap(cmp);

    ///////////////////////////////////////////////////////
    // Perform R-Adaptivity
    ///////////////////////////////////////////////////////

    gsOptim<real_t>::LBFGS optimizer;
    optimizer.options().setInt("MaxIterations",maxIt);
    optimizer.options().setInt("Verbose",1);
    optimizer.options().setReal("GradErrTol",tol_g);

    gsParaviewCollection errcol(dirname+"errors");
    gsParaviewCollection solcol(dirname+"solution");

    gsMultiPatch<> sol;
    real_t totalError, maxError;
    for (index_t i = 0; i < 5; i++)
    {
        writeSingleCompMesh(cmp.basis(0),cmp.patch(0),dirname+"domain_mesh_"+util::to_string(i));
        errcol.addTimestep("domain_mesh_"+util::to_string(i),i,".vtp");
        solcol.addTimestep("domain_mesh_"+util::to_string(i),i,".vtp");
        writeSingleControlNet(domain.domain(),dirname+"domain_cnet_"+util::to_string(i));
        errcol.addTimestep("domain_cnet_"+util::to_string(i),i,".vtp");
        solcol.addTimestep("domain_cnet_"+util::to_string(i),i,".vtp");


        sol = solvePoisson(cmp,cmb,mb,cf,bc);
        gsMultiPatch<> parsol;
        parsol.addPatch(sol.basis(0).makeGeometry(sol.patch(0).coefs()));


        gsExprEvaluator<real_t> ev;
        ev.setIntegrationElements(mb);
        ev.options().setReal("quA",4.0);
        ev.options().setInt("quB",2.0);
        auto G     = ev.getMap(mp);
        auto cG    = ev.getMap(cmp);
        auto u_sol = ev.getVariable(sol);
        auto u_ex  = ev.getVariable(cms);


        gsL2Difference<real_t> err(ms,parsol.patch(0));
        gsL2Difference<real_t> cerr(cms,sol.patch(0));
        // gsL2Difference<real_t> err(cms,sol.patch(0));
        auto E     = ev.getVariable(err);
        gsInfo<<"ev.integral((u_sol-u_ex).norm()*meas(G)) = "<<ev.integral((u_sol-u_ex).norm()*meas(cG))<<"\n";
        gsInfo<<"ev.integral(E*meas(G))                   = "<<ev.integral(E*meas(cG))<<"\n";
        totalError = ev.integral((u_sol-u_ex).norm()*meas(cG));
        gsInfo<<"total error after iteration "<<i<<": "<<totalError<<"\n";


        maxError = ev.max((u_sol-u_ex).norm());

        writeSinglePatchField(cmp.patch(0),cerr,true,dirname+"error_"+util::to_string(i),1000);
        errcol.addTimestep("error_"+util::to_string(i),i,".vts");
        writeSinglePatchField(cmp.patch(0),parsol.patch(0),true,dirname+"sol_"+util::to_string(i),1000);
        solcol.addTimestep("sol_"+util::to_string(i),i,".vts");

        gsWriteParaview(mp,parsol,"parsol_mp");
        gsWriteParaview(cmp,parsol,"parsol_cmp");
        gsWriteParaview(mp,sol,"sol_mp");
        gsWriteParaview(cmp,sol,"sol_cmp");

        // writeSinglePatchField(domain.domain(),sol.patch(0),true,dirname+"sol_"+util::to_string(i),1000);
        // errcol.addTimestep("sol_"+util::to_string(i),i,".vts");

        ev.writeParaview((u_sol),cG,dirname+"u_sol");
        ev.writeParaview((u_ex),cG,dirname+"u_ex");
        ev.writeParaview((u_ex-u_sol).norm(),cG,dirname+"u_err");


        // POSSIBLE PROBLEM: err IS DOUBLE COMPOSED (Once here and once in the gsAdaptiveParametrization class)


        // gsMultiPatch<> tmp = mp, tmp_sol;
        // tmp.patch(0).coefs() = domain.domain().coefs();
        // tmp_sol = solvePoisson(tmp,mb,f,bc);

        // gsExprEvaluator<real_t> tmp_ev;
        // tmp_ev.setIntegrationElements(mb);
        // auto tmp_G     = tmp_ev.getMap(tmp);
        // auto tmp_u_sol = tmp_ev.getVariable(tmp_sol);
        // auto tmp_u_ex  = tmp_ev.getVariable(ms);

        // totalError = tmp_ev.integral((tmp_u_sol-tmp_u_ex).norm()*meas(tmp_G));
        // gsInfo<<"total error after iteration "<<i<<": "<<totalError<<"\n";

        gsInfo<<"R-Adaptivity iteration "<<i<<"\n";
        gsL2Difference<real_t> err_scaled(ms,parsol.patch(0),10./maxError);
        gsAdaptiveParametrization<real_t,MonitorMode::ValueBased> relocator(domain,mp.patch(0),err_scaled,mb.basis(0),optimizer,true);
        // gsAdaptiveParametrization<real_t,MonitorMode::ValueBased> relocator(domain,mp.patch(0),parsol.patch(0),mb.basis(0),optimizer,true);
        // gsAdaptiveParametrization<real_t,MonitorMode::GradientBased> relocator(domain,mp.patch(0),parsol.patch(0),mb.basis(0),optimizer,true);
        relocator.solve();
    }

    errcol.save();
    solcol.save();

    // gsMultiPatch<> sol = solvePoisson(cmp,cmb,mb,cf,bc);

    // ///////////////////////////////////////////////////////
    // // ERROR ESTIMATION
    // ///////////////////////////////////////////////////////

    // gsExprEvaluator<real_t> ev;
    // ev.setIntegrationElements(mb);
    // auto G     = ev.getMap(cmp);
    // auto u_sol = ev.getVariable(sol);
    // auto u_ex  = ev.getVariable(ms);


    // ev.writeParaview(u_sol,G,dirname+"solution");

    // real_t totalError = ev.integral((u_sol-u_ex).norm()*meas(G));
    // gsDebugVar(totalError);

    // ev.writeParaview((u_sol-u_ex).norm()/totalError,G,dirname+"error");
    // // ev.writeParaview(ijac(u_sol,G),G,dirname+"solution_gradient");
    // // ev.writeParaview(ijac(u_sol,G).sqNorm(),G,dirname+"solution_gradient_sqNorm");

    // gsL2Difference<real_t> err(ms.piece(0),sol.patch(0),1./totalError);
    // gsWriteParaview(mp,err,"error",3000);

    // gsComposedFunction<real_t> cerr = gsComposedFunction<real_t>(domain,err);
    // gsComposedFunction<real_t> cserr = gsComposedFunction<real_t>(cspline,err);
    // gsWriteParaview(mp,cerr,"cerror",3000);
    // gsWriteParaview(mp,cserr,"cserror",3000);


    // gsAdaptiveParametrization<real_t,MonitorMode::ValueBased> relocator(domain,mp.patch(0),err,mb.basis(0),optimizer,true);
    // relocator.solve();


    // gsL2Difference<real_t> err2(ms.piece(0),sol.patch(0),1./totalError);
    // gsComposedFunction<real_t> cerr2 = gsComposedFunction<real_t>(domain,err2);
    // gsComposedFunction<real_t> cserr2 = gsComposedFunction<real_t>(cspline,err2);
    // gsWriteParaview(mp,cerr2,"cerror2",3000);
    // gsWriteParaview(mp,cserr2,"cserror2",3000);

    // gsDebugVar(domain.domain().coefs());


    // // gsWriteParaview(cspline,"cspline",1000,true);
    // // gsWriteParaview(cspline.basis(),"cbasis",1000);
    gsWriteParaview(domain.domain(),dirname+"domain",1000,true,true);

    // gsInfo<<"Area = "<<ev.integral(meas(G))<<"\n";
    // ev.writeParaview(jac(G).det(),G,dirname+"jacobian_determinant");



    return 0;
}// end main
