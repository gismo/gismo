/** @file poisson2_example.cpp

    @brief Tutorial on how to use expression assembler to solve the Poisson equation

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris
*/

//! [Include namespace]
#include <gismo.h>

#  define MatExprType  auto

namespace gismo{ 
namespace expr{

template<class E>
class mygrad_expr : public _expr<mygrad_expr<E> >
{
    typename E::Nested_t _u;
public:
    enum{ Space = E::Space };

    typedef typename E::Scalar Scalar;

    mygrad_expr(const E & u) : _u(u)
    { GISMO_ASSERT(1==u.dim(),"grad(.) requires 1D variable, use jac(.) instead.");}

    MatExprType eval(const index_t k) const
    {
        // numActive x dim
        return _u.data().values[1].reshapeCol(k, cols(), rows()).transpose();
    }

    index_t rows() const
    {
        //return _u.data().values[0].rows();
        return _u.data().values[1].rows() / cols();
    }
    //index_t rows() const { return _u.data().actives.size(); }
    //index_t rows() const { return _u.rows(); }

    //index_t rows() const { return _u.source().targetDim() is wrong }
    index_t cols() const { return _u.source().domainDim(); }

    void setFlag() const
    {
        _u.data().flags |= NEED_GRAD;
        if (_u.composed() )
            _u.mapData().flags |= NEED_VALUE;
    }

    void parse(gsSortedVector<const gsFunctionSet<Scalar>*> & evList) const
    {
        //GISMO_ASSERT(NULL!=m_fd, "FeVariable: FuncData member not registered");
        evList.push_sorted_unique(&_u.source());
        _u.data().flags |= NEED_GRAD;
        if (_u.composed() )
            _u.mapData().flags |= NEED_VALUE;
    }

    const gsFeSpace<Scalar> & rowVar() const { return _u.rowVar(); }
    const gsFeSpace<Scalar> & colVar() const
    {return gsNullExpr<Scalar>::get();}

    static constexpr bool rowSpan() {return E::rowSpan(); }
    static constexpr bool colSpan() {return false;}

    void print(std::ostream &os) const { os << "grad("; _u.print(os); os <<")"; }
};

template<class E> EIGEN_STRONG_INLINE
mygrad_expr<E> mygrad(const E & u) { return mygrad_expr<E>(u); }

}}

using namespace gismo;

//! [Include namespace]



int main(int argc, char *argv[])
{
    //! [Parse command line]
    bool plot = false;
    index_t numRefine  = 5;
    index_t numElevate = 0;
    bool last = false;
    std::string fn("pde/poisson2d_bvp.xml");

    gsCmdLine cmd("Tutorial on solving a Poisson problem.");
    cmd.addInt( "e", "degreeElevation",
                "Number of degree elevation steps to perform before solving (0: equalize degree in all directions)", numElevate );
    cmd.addInt( "r", "uniformRefine", "Number of Uniform h-refinement steps to perform before solving",  numRefine );
    cmd.addString( "f", "file", "Input XML file", fn );
    cmd.addSwitch("last", "Solve solely for the last level of h-refinement", last);
    cmd.addSwitch("plot", "Create a ParaView visualization file with the solution", plot);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }
    //! [Parse command line]

    //! [Read input file]

    gsFileData<> fd(fn);
    gsInfo << "Loaded file "<< fd.lastPath() <<"\n";

    gsMultiPatch<> mp;
    fd.getId(0, mp); // id=0: Multipatch domain

    gsFunctionExpr<> f;
    fd.getId(1, f); // id=1: source function
    gsInfo<<"Source function "<< f << "\n";

    gsBoundaryConditions<> bc;
    fd.getId(2, bc); // id=2: boundary conditions
    gsInfo<<"Boundary conditions:\n"<< bc <<"\n";

    gsOptionList Aopt;
    fd.getId(4, Aopt); // id=4: assembler options

    //! [Read input file]

    //! [Refinement]
    gsMultiBasis<> dbasis(mp);

    // Elevate and p-refine the basis to order p + numElevate
    // where p is the highest degree in the bases
    dbasis.setDegree( dbasis.maxCwiseDegree() + numElevate);

    // h-refine each basis
    if (last)
    {
        for (int r =0; r < numRefine-1; ++r)
            dbasis.uniformRefine();
        numRefine = 0;
    }

    gsInfo << "Patches: "<< mp.nPatches() <<", degree: "<< dbasis.minCwiseDegree() <<"\n";
    //! [Refinement]

    //! [Problem setup]
    gsExprAssembler<> A(1,1);
    A.setOptions(Aopt);

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

    // Set the discretization space
    space u = A.getSpace(dbasis);
    u.setInterfaceCont(0);
    u.addBc( bc.get("Dirichlet") );

    // Set the source term
    variable ff = A.getCoeff(f, G);

    // Recover manufactured solution
    gsFunctionExpr<> ms;
    fd.getId(3, ms); // id=3: reference solution
    //gsInfo<<"Exact solution: "<< ms << "\n";
    variable u_ex = ev.getVariable(ms, G);

    // Solution vector and solution variable
    gsMatrix<> solVector;
    solution u_sol = A.getSolution(u, solVector);

    gsSparseSolver<>::CGDiagonal solver;

    //! [Problem setup]

    //! [Solver loop]
    gsVector<> l2err(numRefine+1), h1err(numRefine+1);
    gsInfo<< "(dot1=assembled, dot2=solved, dot3=got_error)\n"
        "\nDoFs: ";
    for (int r=0; r<=numRefine; ++r)
    {
        dbasis.uniformRefine();

        // Initialize the system
        A.initSystem();

        gsInfo<< A.numDofs() <<std::flush;

        // Compute the system matrix and right-hand side
        A.assemble( mygrad(u)*jac(G).ginv() * (mygrad(u)*jac(G).ginv()).tr() * meas(G), u * ff * meas(G) );

        // Enforce Neumann conditions to right-hand side
        variable g_N = A.getBdrFunction();
        A.assembleRhsBc(u * g_N.val() * nv(G).norm(), bc.neumannSides() );
        //gsInfo<<"Sparse Matrix:\n"<< A.matrix().toDense() <<"\n";
        //gsInfo<<"Rhs vector:\n"<< A.rhs().transpose() <<"\n";

        gsInfo<< "." <<std::flush;// Assemblying done

        solver.compute( A.matrix() );
        solVector = solver.solve(A.rhs());

        gsInfo<< "." <<std::flush; // Linear solving done

        l2err[r]= math::sqrt( ev.integral( (u_ex - u_sol).sqNorm() * meas(G) ) );
        h1err[r]= l2err[r] +
        math::sqrt(ev.integral( ( igrad(u_ex) - grad(u_sol)*jac(G).inv() ).sqNorm() * meas(G) ));

        gsInfo<< ". " <<std::flush; // Error computations done

    } //for loop

    //! [Solver loop]

    //! [Error and convergence rates]
    gsInfo<< "\n\nL2 error: "<<std::scientific<<std::setprecision(3)<<l2err.transpose()<<"\n";
    gsInfo<< "H1 error: "<<std::scientific<<h1err.transpose()<<"\n";

    if (!last && numRefine>0)
    {
        gsInfo<< "\nEoC (L2): " << std::fixed<<std::setprecision(2)
              << ( l2err.head(numRefine).array() /
                   l2err.tail(numRefine).array() ).log().transpose() / std::log(2.0) <<"\n";

        gsInfo<<   "EoC (H1): "<< std::fixed<<std::setprecision(2)
              <<( h1err.head(numRefine).array() /
                  h1err.tail(numRefine).array() ).log().transpose() / std::log(2.0) <<"\n";
    }
    //! [Error and convergence rates]

    // if (save)

    //! [Export visualization in ParaView]
    if (plot)
    {
        gsInfo<<"Plotting in Paraview...\n";
        ev.options().setSwitch("plot.elements", true);
        ev.writeParaview( u_sol   , G, "solution");
        //ev.writeParaview( u_ex    , G, "solution_ex");
        //ev.writeParaview( u, G, "aa");

        gsFileManager::open("solution.pvd");
    }
    //! [Export visualization in ParaView]

    return EXIT_SUCCESS;

}// end main
