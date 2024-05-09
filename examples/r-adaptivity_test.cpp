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
#include <gsNurbs/gsMobiusDomain.h>
#include <gsOptimizer/gsGradientDescent.h>

using namespace gismo;
//! [Include namespace]

template <typename T>
// class solvePoisson : public gsOptProblem<T>
class solvePde : public gsOptProblem<T>
//! [OptProblemExample Class]
{

    typedef typename gsExprAssembler<T>::geometryMap geometryMapT;
    typedef typename gsExprAssembler<T>::variable    variableT;
    typedef typename gsExprAssembler<T>::space       spaceT;
    typedef typename gsExprAssembler<T>::solution    solutionT;

public:
    solvePde(   const gsBasis<T> &              basis,
                const gsGeometry<T> &           geom,
                gsFunction<T> &                 composition,
                const gsBoundaryConditions<T> & BCs,
                const gsFunction<T> *           exact,
                const bool                      bmorph = false,
                const bool                      gmorph = false)
    :
    m_basis(basis),
    m_geom(geom),
    m_composition(composition),
    m_BCs(BCs),
    m_source(nullptr),
    m_exact(exact),
    m_Bmorph(bmorph),
    m_Gmorph(gmorph)
    {
        // Compute the number of design variables
        m_numDesignVars = m_composition.nControls();

        // Number of constraints
        m_numConstraints = 0;

        // Number of non-zeros in the Jacobian of the constraints
        m_numConJacNonZero = 0;

        // we initialize x in bounds, in the upper right quadrant
        m_curDesign.resize(m_numDesignVars,1);
        m_curDesign = m_composition.controls();

        m_desLowerBounds.resize(m_numDesignVars);
        m_desUpperBounds.resize(m_numDesignVars);

        // m_desLowerBounds.array() = std::numeric_limits<T>::min();
        // m_desUpperBounds.array() = std::numeric_limits<T>::max();

        m_desLowerBounds.array() = m_curDesign.col(0).array()-0.1;
        m_desUpperBounds.array() = m_curDesign.col(0).array()+0.1;

        m_conJacRows.resize(m_numConJacNonZero);
        m_conJacCols.resize(m_numConJacNonZero);
    }

    void setSource(const gsFunction<T> & source)
    {
        m_source = &source;
    }

    void setModes(bool bmorph, bool gmorph)
    {
        m_Bmorph = bmorph;
        m_Gmorph = gmorph;
    }

    T evalObj( const gsAsConstVector<T> & C) const
    {
        return evalObj(C,false,false);
    }

    T evalObj( const gsAsConstVector<T> & C, bool computeErrors, bool plot) const
    {
        // apply sigma
        if(m_Bmorph || m_Gmorph)
            m_composition.controls() = C;
        // m_composition.updateGeom();

        typename gsBasis<T>::Ptr cbasis;
        if(m_Bmorph)
          cbasis = memory::make_shared(new gsComposedBasis<T>(m_composition,m_basis)); // The basis is composed by the square domain
        else
          cbasis = memory::make_shared_not_owned(&m_basis); // The basis is NOT composed

        typename gsGeometry<T>::Ptr cgeom;
        if(m_Gmorph)
          cgeom = memory::make_shared(new gsComposedGeometry<T>(m_composition,m_geom));
        else
          cgeom = memory::make_shared_not_owned(&m_geom);

        gsMultiPatch<T> mp;
        mp.addPatch(*cgeom);

        gsMultiBasis<T> dbasis(*cbasis);

        // gsBoundaryConditions<T> BCs = m_BCs;
        m_BCs.setGeoMap(mp);

        //! [Problem setup]
        gsExprAssembler<T> A(1,1);

        // Elements used for numerical integration
        A.setIntegrationElements(dbasis);
        gsExprEvaluator<T> ev(A);

        // Set the geometry map
        geometryMapT G = A.getMap(mp);

        // Set the discretization space
        spaceT u = A.getSpace(dbasis);

        // Recover manufactured solution
        auto u_ex = ev.getVariable(*m_exact, G);

        // Solution vector and solution variable
        solutionT u_sol = A.getSolution(u, m_solution);

        //! [Problem setup]

        typename gsSparseSolver<T>::CGDiagonal solver;

        u.setup(m_BCs, dirichlet::l2Projection, 0);

        // Initialize the system
        A.initSystem();
        m_numDofs = A.numDofs();

        if (m_source != nullptr)
        {
            // Set the source term
            auto ff = A.getCoeff(*m_source, G);

            A.assemble(
                        igrad(u,G) * igrad(u,G).tr() * meas(G) //matrix
                        ,
                        u * ff * meas(G) //rhs vector
                        );
        }
        else
        {
            //gsInfo << "Assemble L2 projection.\n";
            A.assemble(
                        u * u.tr() * meas(G),
                        u * u_ex * meas(G)
                        ); // to check
        }

        solver.compute( A.matrix() );
        m_solution = solver.solve(A.rhs());

        if (computeErrors)
        {
            m_errors.resize(0);
            // if (m_source != nullptr)
            // {
            //     auto ff = A.getCoeff(*m_source, G);
            //     ev.integralElWise( (ilapl(u_sol,G) + ff).sqNorm() * meas(G) );
            // }
            // else
                ev.integralElWise(  (u_ex - u_sol).sqNorm() * meas(G) );

            m_errors = ev.elementwise();
        }
        if (plot)
        {
            m_spline = u_sol.extractPiece(0);
            m_finalGeom = cgeom->clone();
        }

        return math::sqrt( ev.integral( (u_ex - u_sol).sqNorm() * meas(G) ) );

        // Improvement:
        // if (m_source != nullptr)
        //     return ev.integral( (ilapl(u_sol,G) + ff).sqNorm() * meas(G) );
        // else
        //    return math::sqrt( ev.integral( (u_ex - u_sol).sqNorm() * meas(G) ) );
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

    // gsMatrix<T> solution() const {return m_solution;}

    gsMatrix<T> currentDesign() const {return m_curDesign;}

    const std::vector<T> & errors() {return m_errors;}

    // Return the number of DoFs of the analysis basis (!)
    index_t numDofs() const {return m_numDofs;}

    typename gsGeometry<T>::Ptr geometry() { return m_finalGeom; }
    typename gsGeometry<T>::Ptr solution() { return m_spline; }

private:

    const gsBasis<T> & m_basis;
    const gsGeometry<T> & m_geom;
    gsFunction<T> & m_composition;
    mutable gsBoundaryConditions<T> m_BCs;
    const gsFunction<T> * m_source;
    const gsFunction<T> * m_exact;
    const bool m_Bmorph;
    const bool m_Gmorph;

    mutable gsMatrix<T> m_solution;

    mutable std::vector<T> m_errors;

    mutable index_t m_numDofs;

    mutable typename gsGeometry<T>::Ptr m_spline;
    mutable typename gsGeometry<T>::Ptr m_finalGeom;


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
    index_t numSquareElev  = 1; // q, uniform refiment
    index_t numRefine  = 1; // r, uniform refiment
    index_t numGeoRefine  = 1; // R, uniform refiment
    index_t numSquareRefine  = 1; // s, uniform refiment
    //t, u, v, w, x, y, z
    bool plot = false;
    bool ref_last = false;  // last, perform last refinement only
    bool writeData = false;

    gsCmdLine cmd("Tutorial on solving a Poisson problem.");
    cmd.addInt( "e", "degreeElevation",
                "Number of degree elevation steps to perform before solving (0: equalize degree in all directions)", numElevate );
    cmd.addInt( "r", "uniformRefine", "Number of Uniform h-refinement loops",  numRefine );
    cmd.addInt( "R", "uniformRefineGeom", "Number of uniform refinements of the geometry",  numGeoRefine );
    cmd.addInt( "s", "uniformRefineSquare", "Number of uniform refinements of the sigma square",  numSquareRefine );
    cmd.addInt( "q", "degElevSquare", "Number of degree elevation of the sigma square",  numSquareElev );
    cmd.addInt("p", "pde", "specify the pde problem", pdeDef);
    cmd.addSwitch("b", "basis", "Apply composition to the basis", b_morph);
    cmd.addSwitch("g", "geom", "Apply composition to the geometry", g_morph);
    cmd.addSwitch("last", "one-shot uniform refinement", ref_last);
    cmd.addSwitch("plot", "Create a ParaView visualization file with the solution", plot);
    cmd.addSwitch("csv", "Save the data in .csv file", writeData);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }
    //! [Parse command line]

    gsStopwatch gsTime;
    time_t now = time(0);

    std::ofstream file_out;
    if (writeData)
    {
      file_out.open(std::to_string(now)+"radaptivity_results.csv");
      file_out << "problem, cbasis, cgeom, a-deg, a-ref, a-dofs, opt-deg, opt-dofs, geo-deg, geo-dofs, L2err\n";
    }


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


    gsFunctionExpr<> f("((tanh(20*(x^2 + y^2)^(1/2) - 5)^2 - 1)*(20*x^2 + 20*y^2)*(40*tanh(20*(x^2 + y^2)^(1/2) - 5)*(x^2 + y^2)^(1/2) - 1))/(x^2 + y^2)^(3/2)",spaceDim);
    gsFunctionExpr<> ms("tanh((0.25-sqrt(x^2+y^2))/0.05)+1",spaceDim);

    gsBoundaryConditions<> bc;
    if(pdeDef == 0)
    {
        bc.addCondition(boundary::side::west ,condition_type::dirichlet,&ms);
        bc.addCondition(boundary::side::east ,condition_type::dirichlet,&ms);
        bc.addCondition(boundary::side::south,condition_type::dirichlet,&ms);
        bc.addCondition(boundary::side::north,condition_type::dirichlet,&ms);
    }

    // gsSquareDomain<2,real_t> composition(numSquareElev,numSquareRefine);
    gsMobiusDomain<2,real_t> composition;

    gsMultiPatch<> mp0;
    mp0.addPatch(gsNurbsCreator<>::BSplineSquare());
    mp0.patch(0).coefs().array() -= 0.5;
    if (spaceDim > 2)
        mp0.embed(3);

    if (numElevate!=0)
        mp0.degreeElevate(numElevate);

    gsMultiBasis<> basis(mp0);

    for (index_t k=0; k!=numGeoRefine; k++)
        mp0.uniformRefine();


    solvePde<real_t> problem(basis.basis(0),mp0.patch(0),composition,bc,&ms,b_morph,g_morph);
    if (pdeDef == 0)
        problem.setSource(f);

    gsMatrix<> currentDesign;

    gsHLBFGS<real_t> solver(&problem);
    solver.options().setInt("Verbose",2);
    solver.options().setReal("MinGradientLength", 1e-3); // 1e-6 : more or less as refinement tolerance; there should be a balance between the two;
    solver.options().setInt("LBFGSUpdates", 20);
    solver.options().setReal("MinStepLength", 1e-12);
    solver.options().setInt("MaxIterations",1000); // this can help

    // gsGradientDescent<real_t> solver(&problem);
    // solver.options().setInt("Verbose",1);
    // solver.options().setReal("MinGradientLength", 1e-6); // 1e-6 : more or less as refinement tolerance; there should be a balance between the two;
    // solver.options().setReal("MinStepLength", 1e-6);
    // solver.options().setInt("MaxIterations",1000); // this can help

    for(index_t refCount = 1; refCount <= numRefine; refCount++)
    {
        if(ref_last)
        {
            refCount = numRefine+1;
            index_t numKnts = pow(2, numRefine) -1;
            basis.uniformRefine(numKnts);
        }
        else
        {
            basis.uniformRefine();
        }

        gsDebugVar(mp0.patch(0).basis().size());
        currentDesign = problem.currentDesign();

        real_t L2err;
        if (b_morph || g_morph)
        {
            gsInfo<<"Solving r-adaptivity..."<<std::flush;
            solver.solve(currentDesign);
            currentDesign = solver.currentDesign();
            gsAsConstVector<real_t> curr(currentDesign.data(),currentDesign.rows());
            L2err = problem.evalObj(curr,true,true);
            gsInfo<<"Finished\n";
        }
        else
        {
            gsInfo<<"Solving without r-adaptivity..."<<std::flush;
            L2err = problem.evalObj(gsAsConstVector<real_t>(currentDesign.data(),currentDesign.rows()),true,true);
            gsInfo<<"Finished\n";
        }

        //! [Export visualization in ParaView]
        std::cout << "\nL2 error = " << std::setprecision(10) << std::scientific << L2err<<"\n";

        // file_out << "problem, cbasis, cgeom, a-deg, a-ref, a-dofs, opt-deg, opt-dofs, geo-deg, geo-dofs, L2err\n";
        if (writeData)
        {
          file_out << problem_name << ","
                 << b_morph << ","
                 << g_morph << ","
                 << std::max(mp0.patch(0).degree(0),mp0.patch(0).degree(1))   << "," // analysis degree
                 << refCount << "," // analysis unif refinement
                 << problem.numDofs()  << "," // analysis dofs
                 // << composition.maxDegree() << "," // opt degree
                 << "" << ","
                 << composition.nControls() << "," // optDim
                 << mp0.patch(0).basis().maxDegree() << ","
                 << mp0.patch(0).basis().size() << ","
                 << std::setprecision(10) << std::scientific << L2err<<"\n";
        }


        if (plot)
        {
            gsField<> field(*problem.geometry(),*problem.solution());
            gsWriteParaview(field,"radaptivity_solution",1000,true);
            gsWriteParaview(*problem.geometry(),"radaptivity_geometry",1000,true);
            composition.controls() = currentDesign;
            // composition.updateGeom();
            // gsWriteParaview(composition.domain(),"radaptivity_composition",1000,true,true);

        }

    }
    
    if (writeData)
      file_out.close();
    return EXIT_SUCCESS;

}// end main
