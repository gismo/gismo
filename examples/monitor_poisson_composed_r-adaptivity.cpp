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

using namespace gismo;

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
    std::string geometry_file;

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
    cmd.addString( "G", "geometry", "Geometry file", geometry_file );

    try { cmd.getValues(arg,argv); } catch (int rv) { return rv; }
    //! [Parse command line]

    GISMO_ASSERT(!geometry_file.empty(),"Please provide a geometry file.");

    std::string dirname = "r-adaptivity_composed_poisson_e"+util::to_string(numElevate)+"_r"+util::to_string(numRefine)+"_E"+util::to_string(numElevateD)+"_R"+util::to_string(numRefineD);
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

    // Source function:
    const index_t dimexpr = mp.geoDim();
    std::string fstring = "-9*pi^2*(2*x - 1.0)^2*cos(3*pi*(sqrt((x - 0.5)^2 + (y - 0.5)^2) - 1/6))^2*exp(-162*(sqrt((x - 0.5)^2 + (y - 0.5)^2) - 1/3)^2)/(2*((x - 0.5)^2 + (y - 0.5)^2)) + 972*sin(3*pi*(sqrt((x - 0.5)^2 + (y - 0.5)^2) - 1/6))*(sqrt((x - 0.5)^2 + (y - 0.5)^2) - 1/3)*(2*x - 1.0)^2*exp(-162*(sqrt((x - 0.5)^2 + (y - 0.5)^2) - 1/3)^2)*pi*cos(3*pi*(sqrt((x - 0.5)^2 + (y - 0.5)^2) - 1/6))/((x - 0.5)^2 + (y - 0.5)^2) + (3*sin(3*pi*(sqrt((x - 0.5)^2 + (y - 0.5)^2) - 1/6))*exp(-162*(sqrt((x - 0.5)^2 + (y - 0.5)^2) - 1/3)^2)*pi*(2*x - 1.0)^2*cos(3*pi*(sqrt((x - 0.5)^2 + (y - 0.5)^2) - 1/6)))/(2*((x - 0.5)^2 + (y - 0.5)^2)^(3/2)) - 12*sin(3*pi*(sqrt((x - 0.5)^2 + (y - 0.5)^2) - 1/6))*exp(-162*(sqrt((x - 0.5)^2 + (y - 0.5)^2) - 1/3)^2)*pi*cos(3*pi*(sqrt((x - 0.5)^2 + (y - 0.5)^2) - 1/6))/sqrt((x - 0.5)^2 + (y - 0.5)^2) + 9*sin(3*pi*(sqrt((x - 0.5)^2 + (y - 0.5)^2) - 1/6))^2*exp(-162*(sqrt((x - 0.5)^2 + (y - 0.5)^2) - 1/3)^2)*pi^2*(2*x - 1.0)^2/(2*((x - 0.5)^2 + (y - 0.5)^2)) + 81*sin(3*pi*(sqrt((x - 0.5)^2 + (y - 0.5)^2) - 1/6))^2*(2*x - 1.0)^2*exp(-162*(sqrt((x - 0.5)^2 + (y - 0.5)^2) - 1/3)^2)/((x - 0.5)^2 + (y - 0.5)^2) - 81*sin(3*pi*(sqrt((x - 0.5)^2 + (y - 0.5)^2) - 1/6))^2*(sqrt((x - 0.5)^2 + (y - 0.5)^2) - 1/3)*(2*x - 1.0)^2*exp(-162*(sqrt((x - 0.5)^2 + (y - 0.5)^2) - 1/3)^2)/((x - 0.5)^2 + (y - 0.5)^2)^(3/2) + 648*sin(3*pi*(sqrt((x - 0.5)^2 + (y - 0.5)^2) - 1/6))^2*(sqrt((x - 0.5)^2 + (y - 0.5)^2) - 1/3)*exp(-162*(sqrt((x - 0.5)^2 + (y - 0.5)^2) - 1/3)^2)/sqrt((x - 0.5)^2 + (y - 0.5)^2) - 26244*sin(3*pi*(sqrt((x - 0.5)^2 + (y - 0.5)^2) - 1/6))^2*(sqrt((x - 0.5)^2 + (y - 0.5)^2) - 1/3)^2*(2*x - 1.0)^2*exp(-162*(sqrt((x - 0.5)^2 + (y - 0.5)^2) - 1/3)^2)/((x - 0.5)^2 + (y - 0.5)^2) - 9*pi^2*(2*y - 1.0)^2*cos(3*pi*(sqrt((x - 0.5)^2 + (y - 0.5)^2) - 1/6))^2*exp(-162*(sqrt((x - 0.5)^2 + (y - 0.5)^2) - 1/3)^2)/(2*((x - 0.5)^2 + (y - 0.5)^2)) + 972*sin(3*pi*(sqrt((x - 0.5)^2 + (y - 0.5)^2) - 1/6))*(sqrt((x - 0.5)^2 + (y - 0.5)^2) - 1/3)*(2*y - 1.0)^2*exp(-162*(sqrt((x - 0.5)^2 + (y - 0.5)^2) - 1/3)^2)*pi*cos(3*pi*(sqrt((x - 0.5)^2 + (y - 0.5)^2) - 1/6))/((x - 0.5)^2 + (y - 0.5)^2) + (3*sin(3*pi*(sqrt((x - 0.5)^2 + (y - 0.5)^2) - 1/6))*exp(-162*(sqrt((x - 0.5)^2 + (y - 0.5)^2) - 1/3)^2)*pi*(2*y - 1.0)^2*cos(3*pi*(sqrt((x - 0.5)^2 + (y - 0.5)^2) - 1/6)))/(2*((x - 0.5)^2 + (y - 0.5)^2)^(3/2)) + 9*sin(3*pi*(sqrt((x - 0.5)^2 + (y - 0.5)^2) - 1/6))^2*exp(-162*(sqrt((x - 0.5)^2 + (y - 0.5)^2) - 1/3)^2)*pi^2*(2*y - 1.0)^2/(2*((x - 0.5)^2 + (y - 0.5)^2)) + 81*sin(3*pi*(sqrt((x - 0.5)^2 + (y - 0.5)^2) - 1/6))^2*(2*y - 1.0)^2*exp(-162*(sqrt((x - 0.5)^2 + (y - 0.5)^2) - 1/3)^2)/((x - 0.5)^2 + (y - 0.5)^2) - 81*sin(3*pi*(sqrt((x - 0.5)^2 + (y - 0.5)^2) - 1/6))^2*(sqrt((x - 0.5)^2 + (y - 0.5)^2) - 1/3)*(2*y - 1.0)^2*exp(-162*(sqrt((x - 0.5)^2 + (y - 0.5)^2) - 1/3)^2)/((x - 0.5)^2 + (y - 0.5)^2)^(3/2) - 26244*sin(3*pi*(sqrt((x - 0.5)^2 + (y - 0.5)^2) - 1/6))^2*(sqrt((x - 0.5)^2 + (y - 0.5)^2) - 1/3)^2*(2*y - 1.0)^2*exp(-162*(sqrt((x - 0.5)^2 + (y - 0.5)^2) - 1/3)^2)/((x - 0.5)^2 + (y - 0.5)^2)";
    std::string msstring = "sin(3*pi*(sqrt((x - 0.5)^2 + (y - 0.5)^2) - 1/6))^2*exp(-162*(sqrt((x - 0.5)^2 + (y - 0.5)^2) - 1/3)^2)";

    gsFunctionExpr<> f(fstring, dimexpr);
    gsFunctionExpr<> ms(msstring,dimexpr);
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
    auto u_ex = ev.getVariable(cms, G);

    // Solution vector and solution variable
    solution u_sol = A.getSolution(u, solVector);

    //! [Problem setup]

    gsParaviewCollection collection(dirname+"solution",&ev);
    collection.options().setSwitch("plotElements", true);
    const bool export_b64 = false; // was a CLI switch in poisson2_example, lost on copy
    collection.options().setSwitch("base64", export_b64);
    collection.options().setInt("plotElements.resolution", 4);
    collection.options().setInt("plot.npts", 1000);
    gsSparseSolver<real_t>::CGDiagonal solver;
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
        collection.save();
    }

    collection.save();

    return 0;
}// end main
