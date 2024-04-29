/** @file r-adaptivity_test.cpp

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
// class solvePoisson : public gsOptProblem<T>
class solvePde : public gsOptProblem<T>
//! [OptProblemExample Class]
{
public:

    // solvePoisson(   const index_t numElevate,
    //                 const index_t numRefine)
    solvePde(   const index_t pdeDef,
                const index_t spaceDim,
                const index_t numElevate,
                const index_t numRefine,
                const bool composeBasis,
                const bool composeGeom)
    :
    m_problemType(pdeDef),
    m_Dim(spaceDim),
    m_numElevate(numElevate),
    m_numRefine(numRefine),
    m_Bmorph(composeBasis),
    m_Gmorph(composeGeom),
    m_domain(gsSquareDomain<2,T>(1,1))
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

        // gsMultiPatch<T> mp0;
        // mp0.addPatch(gsNurbsCreator<T>::BSplineSquare());
        // mp0.patch(0).coefs().array() -= 0.5; // This I do not remember why we put here.
        //
        // if (m_numElevate!=0)
        //     mp0.degreeElevate(m_numElevate);

        gsMultiPatch<T> mp0;
        mp0.addPatch(gsNurbsCreator<>::BSplineSquare());
        mp0.patch(0).coefs().array() -= 0.5;
        if (m_Dim > 2)
          mp0.embed(m_Dim);

          if (m_numElevate!=0)
              mp0.degreeElevate(m_numElevate);

        // h-refine
        for (int r =0; r < m_numRefine; ++r) // to be removed.
            mp0.uniformRefine();


        // Make composed geometry and basis
        const gsBasis<T> & tbasis = mp0.basis(0); // basis(u,v) -> deriv will give dphi/du ,dphi/dv
        const gsGeometry<T> & tgeom = mp0.patch(0); //G(u,v) -> deriv will give dG/du, dG/dv



        // apply sigma
        if(m_Bmorph || m_Gmorph)
          m_domain.controls() = C;

        m_domain.updateGeom();

        // Define a composite basis and composite geometry

        /* OLD
        gsComposedBasis<real_t> cbasis(m_domain,tbasis); // basis(u,v) = basis(sigma(xi,eta)) -> deriv will give dphi/dxi, dphi/deta
        */

        typename gsBasis<T>::Ptr cbasis;
        if(m_Bmorph)
          cbasis = memory::make_shared(new gsComposedBasis<T>(m_domain,tbasis)); // The basis is composed by the square domain
        else
          cbasis = memory::make_shared_not_owned(&tbasis); // The basis is NOT composed


        // The geometry is defined using the composite basis and some coefficients

        // gsComposedGeometry<real_t> cgeom(m_domain,tgeom);


        typename gsGeometry<T>::Ptr cgeom;
        if(m_Gmorph)
          cgeom = memory::make_shared(new gsComposedGeometry<T>(m_domain,tgeom));
        else
          cgeom = memory::make_shared_not_owned(&tgeom);


        gsMultiPatch<T> mp;
        mp.addPatch(*cgeom);

        // gsMultiBasis<T> dbasis(mp,true); // ??? composed ???
        gsMultiBasis<T> dbasis(*cbasis);


        // Source function:
        gsFunctionExpr<T> f("((tanh(20*(x^2 + y^2)^(1/2) - 5)^2 - 1)*(20*x^2 + 20*y^2)*(40*tanh(20*(x^2 + y^2)^(1/2) - 5)*(x^2 + y^2)^(1/2) - 1))/(x^2 + y^2)^(3/2)",m_Dim);

        // Exact solution
        gsFunctionExpr<T> ms("tanh((0.25-sqrt(x^2+y^2))/0.05)+1",m_Dim);

        gsBoundaryConditions<T> bc;
        if (m_problemType == 0)
        {
          bc.addCondition(boundary::side::west ,condition_type::dirichlet,&ms);
          bc.addCondition(boundary::side::east ,condition_type::dirichlet,&ms);
          bc.addCondition(boundary::side::south,condition_type::dirichlet,&ms);
          bc.addCondition(boundary::side::north,condition_type::dirichlet,&ms);
        }
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


        if (m_problemType == 0)
        {
          //gsInfo << "Assemble poisson.\n";
          // A.assemble( igrad(u, G) * igrad(u, G).tr() * meas(G) //matrix
          //  ,
          //  u * ff * meas(G) //rhs vector
          //  );
          A.assemble(
            (grad(u)*(jac(G).ginv().tr())) * (grad(u)*(jac(G).ginv().tr())).tr() * meas(G) //matrix
            ,
            u * ff * meas(G) //rhs vector
          );
        }
        else if(m_problemType == 1)
        {
          //gsInfo << "Assemble L2 projection.\n";
          A.assemble(
            u * u.tr() * meas(G),
            u * u_ex * meas(G)); // to check
        }
        else
        {
          gsWarn << "Unknown problem, exiting." << std::endl;
          return -1;
        }

        // A.assemble(
        //   (grad(u)*(jac(G).ginv().tr())) * (grad(u)*(jac(G).ginv().tr())).tr() * meas(G) //matrix
        //   ,
        //   u * ff * meas(G) //rhs vector
        // );

        solver.compute( A.matrix() );
        m_solution = solver.solve(A.rhs());

        return math::sqrt( ev.integral( (u_ex - u_sol).sqNorm() * meas(G) ) );
    }

     void gradObj_into( const gsAsConstVector<T> & u, gsAsVector<T> & result) const override
     {
       const index_t n = u.rows();
       //GISMO_ASSERT((index_t)m_numDesignVars == n*m, "Wrong design.");

       gsMatrix<T> uu = u;//copy
       gsAsVector<T> tmp(uu.data(), n);
       gsAsConstVector<T> ctmp(uu.data(), n);
       index_t c = 0;

       const T e0 = this->evalObj(ctmp);
       // for all partial derivatives (column-wise)
       for ( index_t i = 0; i!=n; i++ )
       {
         // to do: add m_desLowerBounds m_desUpperBounds check
         tmp[i]  += T(0.00001);
         const T e1 = this->evalObj(ctmp);
         tmp[i]   = u[i];
         result[c++]= ( e1 - e0 ) / T(0.00001);
//         tmp[i]   = u[i] + T(0.00002);
//         tmp[i]   = u[i] + T(0.00002);
//         const T e3 = this->evalObj(ctmp);
//         tmp[i]   = u[i] - T(0.00001);
//         const T e2 = this->evalObj(ctmp);
//         tmp[i]   = u[i] - T(0.00002);
//         const T e4 = this->evalObj(ctmp);
//         tmp[i]   = u[i];
//         result[c++]= ( 8 * (e1 - e2) + e4 - e3 ) / T(0.00012);
       }
     }

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

    const index_t m_problemType;
    const index_t m_Dim;
    index_t m_numElevate, m_numRefine;
    const bool m_Bmorph, m_Gmorph;
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
    // a
    bool b_morph = false;   // b, composition on the basis
    // c
    index_t spaceDim = 2;   // d, dimension of the physical solution
    index_t numElevate = 0; // e, degree elevation
    // f
    bool g_morph = false;   // g, composition on the geometry
    // i, j, k, l, m, n, o
    index_t pdeDef = 0; // p, definition of the pde to be solved: 0 for poisso; 1 for L2 projection;
    // q
    index_t numRefine  = 1; // r, uniform refiment
    // s, t u, v, w, x, y, z
    bool plot = false;
    bool ref_last = false;  // last, perform last refinement only

    gsCmdLine cmd("Tutorial on solving a Poisson problem.");
    cmd.addInt( "e", "degreeElevation",
                "Number of degree elevation steps to perform before solving (0: equalize degree in all directions)", numElevate );
    cmd.addInt( "r", "uniformRefine", "Number of Uniform h-refinement loops",  numRefine );
    cmd.addInt("p", "pde", "specify the pde problem", pdeDef);
    cmd.addSwitch("b", "basis", "Apply composition to the basis", b_morph);
    cmd.addSwitch("g", "geom", "Apply composition to the geometry", g_morph);
    cmd.addSwitch("last", "one-shot uniform refinement", ref_last);
    cmd.addSwitch("plot", "Create a ParaView visualization file with the solution", plot);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }
    //! [Parse command line]

    gsStopwatch gsTime;
    time_t now = time(0);

    std::ofstream file_out;
    file_out.open(std::to_string(now)+"radaptivity_results.csv");
    file_out << "problem, cbasis, cgeom, deg, ref, dofs, L2err\n";


    std::string problem_name = "";
    if(pdeDef == 0)
    {
      problem_name = "poisson";
      gsInfo << "Solve Poisson equation.\n";
    }
    else if(pdeDef == 1)
    {
      problem_name = "L2_projection";
      gsInfo << "Solve L2 projection.\n";
    }
    else
    {
      gsWarn << "Unknown problem, exiting." << std::endl;
      return -1;
    }


    // solvePoisson<real_t> problem(numElevate,numRefine);
    solvePde<real_t> problem(pdeDef, spaceDim, numElevate, numRefine, b_morph, g_morph);
    gsMatrix<> currentDesign = problem.currentDesign();
    gsHLBFGS<real_t> solver(&problem);
    solver.options().setInt("Verbose",2);
    solver.options().setReal("MinGradientLength", 1e-6); // 1e-6 : more or less as refinement tolerance; there should be a balance between the two;
    solver.options().setInt("LBFGSUpdates", 20);
    solver.options().setReal("MinStepLength", 1e-12);
    solver.options().setInt("MaxIterations",1000); // this can help

    solver.solve(currentDesign);

    /**
     * ONLY FOR RECONSTRUCTION OF THE SOLUTION
     // why do we need to repeat everything here?
     // Is there a more efficient way to do it?
     */
    gsMultiPatch<> mp0;
    mp0.addPatch(gsNurbsCreator<>::BSplineSquare());
    mp0.patch(0).coefs().array() -= 0.5;
    if (spaceDim > 2)
      mp0.embed(3);

    if (numElevate!=0)
        mp0.degreeElevate(numElevate);

    // h-refine
    for (int r =0; r < numRefine; ++r)
        mp0.uniformRefine();


    // Make composed geometry and basis
    const gsBasis<> & tbasis = mp0.basis(0); // basis(u,v) -> deriv will give dphi/du ,dphi/dv
    const gsGeometry<> & tgeom = mp0.patch(0); //G(u,v) -> deriv will give dG/du, dG/dv

    // The domain sigma
    // gsSquareDomain<2,real_t> domain;
    gsSquareDomain<2,real_t> domain(1, 1);
    domain.controls() = solver.currentDesign();
    domain.updateGeom();

    gsDebugVar(solver.currentDesign().minCoeff());
    gsDebugVar(solver.currentDesign().maxCoeff());


    // Define a composite basis and composite geometry

    gsBasis<>::Ptr cbasis;
    if(b_morph)
    {
      gsInfo << "Composed basis.\n";
      cbasis = memory::make_shared(new gsComposedBasis<real_t>(domain,tbasis));
    }
    else
      cbasis = memory::make_shared_not_owned(&tbasis);

    gsGeometry<>::Ptr cgeom;
    if(g_morph)
    {
      gsInfo << "Composed geometry.\n";
      cgeom = memory::make_shared(new gsComposedGeometry<real_t>(domain,tgeom));
    }
    else
      cgeom = memory::make_shared_not_owned(&tgeom);


    gsMultiPatch<> mp;
    mp.addPatch(*cgeom);

    gsMultiBasis<> dbasis(*cbasis);

    gsFunctionExpr<> ms("tanh((0.25-sqrt(x^2+y^2))/0.05)+1",spaceDim);

    gsBoundaryConditions<> bc;
    if(pdeDef == 0)
    {
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
//        collection.options().setSwitch("plotElements", true);
//        collection.options().setInt("plotElements.resolution", 100);
        collection.newTimeStep(&mp);
        collection.addField(u_sol,"numerical solution");
        collection.addField(u_ex, "exact solution");
        collection.addField((u_ex-u_sol).sqNorm(), "error");
        collection.saveTimeStep();
        collection.save();

        gsComposedGeometry<real_t> cgeom2(domain,tgeom);
        gsWriteParaview(*cbasis,"cbasis",2000);
        gsWriteParaview(cgeom2,"cgeom",1000);

        // gsFileManager::open("ParaviewOutput/solution.pvd");
    }
    else
        gsInfo << "Done. No output created, re-run with --plot to get a ParaView "
                  "file containing the solution.\n";
    //! [Export visualization in ParaView]
    file_out.close();
    return EXIT_SUCCESS;

}// end main
