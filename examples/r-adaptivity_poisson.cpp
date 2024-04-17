/** @file r-adaptivity_poisson.cpp

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

template <typename T>
class solvePoisson : public gsOptProblem<T>
//! [OptProblemExample Class]
{
public:

    solvePoisson(   const index_t numElevate,
                    const index_t numRefine)
    :
    m_numElevate(numElevate),
    m_numRefine(numRefine),
    m_domain(gsSquareDomain<2,T>())
    {
        // Number of design variables
        m_numDesignVars = m_domain.nControls();

        // Number of constraints
        m_numConstraints = 0;

        // Number of non-zeros in the Jacobian of the constraints
        m_numConJacNonZero = 0;

        // we initialize x in bounds, in the upper right quadrant
        m_curDesign.resize(m_numDesignVars,1);
        m_curDesign.col(0) = m_domain.controls();


        m_desLowerBounds.resize(m_numDesignVars);
        m_desUpperBounds.resize(m_numDesignVars);

        m_desLowerBounds.array() = std::numeric_limits<T>::min();
        m_desUpperBounds.array() = std::numeric_limits<T>::max();

        // m_desLowerBounds.array() = m_curDesign.col(0).array()-0.1;
        // m_desUpperBounds.array() = m_curDesign.col(0).array()+0.1;

        m_conJacRows.resize(m_numConJacNonZero);
        m_conJacCols.resize(m_numConJacNonZero);
    }

    T evalObj( const gsAsConstVector<T> & C ) const
    {
        // gsDebugVar(C.transpose());
        // gsDebugVar(m_curDesign.transpose());

        gsMultiPatch<T> mp0;
        mp0.addPatch(gsNurbsCreator<T>::BSplineSquare());
        mp0.patch(0).coefs().array() -= 0.5;

        if (m_numElevate!=0)
            mp0.degreeElevate(m_numElevate);

        // h-refine
        for (int r =0; r < m_numRefine; ++r)
            mp0.uniformRefine();


        // Make composed geometry and basis
        const gsBasis<T> & tbasis = mp0.basis(0); // basis(u,v) -> deriv will give dphi/du ,dphi/dv
        const gsGeometry<T> & tgeom = mp0.patch(0); //G(u,v) -> deriv will give dG/du, dG/dv

        // The domain sigma
        m_domain.controls() = C;
        m_domain.updateGeom();

        // Define a composite basis and composite geometry
        // The basis is composed by the square domain
        gsComposedBasis<real_t> cbasis(m_domain,tbasis); // basis(u,v) = basis(sigma(xi,eta)) -> deriv will give dphi/dxi, dphi/deta
        // The geometry is defined using the composite basis and some coefficients
        // gsComposedGeometry<real_t> cgeom(cbasis,tgeom.coefs()); // G(u,v) = G(sigma(xi,eta))  -> deriv will give dG/dxi, dG/deta
        gsComposedGeometry<real_t> cgeom(m_domain,tgeom);

        gsMultiPatch<T> mp;
        mp.addPatch(cgeom);

        gsMultiBasis<T> dbasis(mp,true);


        // Source function:
        gsFunctionExpr<T> f("((tanh(20*(x^2 + y^2)^(1/2) - 5)^2 - 1)*(20*x^2 + 20*y^2)*(40*tanh(20*(x^2 + y^2)^(1/2) - 5)*(x^2 + y^2)^(1/2) - 1))/(x^2 + y^2)^(3/2)",2);

        // Exact solution
        gsFunctionExpr<T> ms("tanh((0.25-sqrt(x^2+y^2))/0.05)+1",2);

        gsBoundaryConditions<T> bc;
        bc.addCondition(boundary::side::west ,condition_type::dirichlet,&ms);
        bc.addCondition(boundary::side::east ,condition_type::dirichlet,&ms);
        bc.addCondition(boundary::side::south,condition_type::dirichlet,&ms);
        bc.addCondition(boundary::side::north,condition_type::dirichlet,&ms);
        bc.setGeoMap(mp);

        //! [Problem setup]
        gsExprAssembler<T> A(1,1);

        typedef typename gsExprAssembler<T>::geometryMap geometryMap;
        typedef typename gsExprAssembler<T>::variable    variable;
        typedef typename gsExprAssembler<T>::space       space;
        typedef typename gsExprAssembler<T>::solution    solution;

        // Elements used for numerical integration
        A.setIntegrationElements(dbasis);
        gsExprEvaluator<T> ev(A);

        // Set the geometry map
        geometryMap G = A.getMap(mp);

        // Set the discretization space
        space u = A.getSpace(dbasis);

        // Set the source term
        auto ff = A.getCoeff(f, G);

        // Recover manufactured solution
        auto u_ex = ev.getVariable(ms, G);

        // Solution vector and solution variable
        solution u_sol = A.getSolution(u, m_solution);

        //! [Problem setup]

        typename gsSparseSolver<T>::CGDiagonal solver;

        u.setup(bc, dirichlet::l2Projection, 0);

        // Initialize the system
        A.initSystem();

        // Compute the system matrix and right-hand side
        A.assemble(
            igrad(u, G) * igrad(u, G).tr() * meas(G) //matrix
            ,
            u * ff * meas(G) //rhs vector
            );

        solver.compute( A.matrix() );
        m_solution = solver.solve(A.rhs());

        return math::sqrt( ev.integral( (u_ex - u_sol).sqNorm() * meas(G) ) );
    }

    // void gradObj_into( const gsAsConstVector<T> & u, gsAsVector<T> & result) const
    // {}

    void evalCon_into( const gsAsConstVector<T> & u, gsAsVector<T> & result) const
    {
        result.resize(m_numConstraints,1);
    }

    // void jacobCon_into( const gsAsConstVector<T> & u, gsAsVector<T> & result) const
    // {

    // }

    gsMatrix<T> solution() const {return m_solution;}

    gsMatrix<T> currentDesign() const {return m_curDesign;}

private:

    index_t m_numElevate, m_numRefine;
    mutable gsMatrix<T> m_solution;
    mutable gsSquareDomain<2,T> m_domain;

    // Lastly, we forward the memebers of the base clase gsOptProblem
    using gsOptProblem<T>::m_numDesignVars;
    using gsOptProblem<T>::m_numConstraints;
    using gsOptProblem<T>::m_numConJacNonZero;

    using gsOptProblem<T>::m_desLowerBounds;
    using gsOptProblem<T>::m_desUpperBounds;

    using gsOptProblem<T>::m_conLowerBounds;
    using gsOptProblem<T>::m_conUpperBounds;

    using gsOptProblem<T>::m_conJacRows;
    using gsOptProblem<T>::m_conJacCols;

    using gsOptProblem<T>::m_curDesign;
};


int main(int argc, char *argv[])
{
    //! [Parse command line]
    bool plot = false;
    index_t numRefine  = 1;
    index_t numElevate = 0;

    gsCmdLine cmd("Tutorial on solving a Poisson problem.");
    cmd.addInt( "e", "degreeElevation",
                "Number of degree elevation steps to perform before solving (0: equalize degree in all directions)", numElevate );
    cmd.addInt( "r", "uniformRefine", "Number of Uniform h-refinement loops",  numRefine );
    cmd.addSwitch("plot", "Create a ParaView visualization file with the solution", plot);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }
    //! [Parse command line]

    solvePoisson<real_t> problem(numElevate,numRefine);
    gsMatrix<> currentDesign = problem.currentDesign();
    gsHLBFGS<real_t> solver(&problem);
    solver.options().setInt("Verbose",2);
    solver.options().setReal("MinGradientLength", 1e-9); // 1e-6 : more or less as refinement tolerance; there should be a balance between the two;
    solver.options().setInt("LBFGSUpdates", 20);
    solver.options().setReal("MinStepLength", 1e-12);
    solver.options().setInt("MaxIterations",100);

    solver.solve(currentDesign);

    /**
     * ONLY FOR RECONSTRUCTION OF THE SOLUTION
     */
    gsMultiPatch<> mp0;
    mp0.addPatch(gsNurbsCreator<>::BSplineSquare());
    mp0.patch(0).coefs().array() -= 0.5;

    if (numElevate!=0)
        mp0.degreeElevate(numElevate);

    // h-refine
    for (int r =0; r < numRefine; ++r)
        mp0.uniformRefine();


    // Make composed geometry and basis
    const gsBasis<> & tbasis = mp0.basis(0); // basis(u,v) -> deriv will give dphi/du ,dphi/dv
    const gsGeometry<> & tgeom = mp0.patch(0); //G(u,v) -> deriv will give dG/du, dG/dv

    // The domain sigma
    gsSquareDomain<2,real_t> domain;
    domain.controls() = solver.currentDesign();
    domain.updateGeom();

    gsDebugVar(solver.currentDesign());


    // Define a composite basis and composite geometry
    // The basis is composed by the square domain
    gsComposedBasis<real_t> cbasis(domain,tbasis); // basis(u,v) = basis(sigma(xi,eta)) -> deriv will give dphi/dxi, dphi/deta
    // The geometry is defined using the composite basis and some coefficients
    gsComposedGeometry<real_t> cgeom(cbasis,tgeom.coefs()); // G(u,v) = G(sigma(xi,eta))  -> deriv will give dG/dxi, dG/deta

    gsMultiPatch<> mp;
    mp.addPatch(cgeom);

    gsMultiBasis<> dbasis(mp,true);

    gsFunctionExpr<> ms("tanh((0.25-sqrt(x^2+y^2))/0.05)+1",2);

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
    u.setup(bc, dirichlet::l2Projection, 0);

    // Recover manufactured solution
    auto u_ex = ev.getVariable(ms, G);

    // Solution vector and solution variable
    gsMatrix<> sol = problem.solution();
    solution u_sol = A.getSolution(u, sol);

    //! [Export visualization in ParaView]
    if (plot)
    {
        gsInfo<<"Plotting in Paraview...\n";

        gsParaviewCollection collection("ParaviewOutput/solution", &ev);
        collection.options().setSwitch("plotElements", true);
        collection.options().setInt("plotElements.resolution", 100);
        collection.newTimeStep(&mp);
        collection.addField(u_sol,"numerical solution");
        collection.addField(u_ex, "exact solution");
        collection.addField((u_ex-u_sol).sqNorm(), "error");
        collection.saveTimeStep();
        collection.save();


        // gsFileManager::open("ParaviewOutput/solution.pvd");
    }
    else
        gsInfo << "Done. No output created, re-run with --plot to get a ParaView "
                  "file containing the solution.\n";
    //! [Export visualization in ParaView]

    return EXIT_SUCCESS;

}// end main
