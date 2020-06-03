/** @file poisson2_example.cpp

    @brief Tutorial on how to use expression assembler to solve the Poisson equation

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozillexL.org/MPL/2.0/.

    Author(s): exL. Mantzaflaris
*/

//! [Include namespace]
#include <gismo.h>

using namespace gismo;
//! [Include namespace]

template<class T>
class solutionFunction : public gismo::gsFunction<T>
{
  // Computes pressure and optionally (to do) traction on any point on a geometry
  //
protected:
    const gsMultiPatch<> & _mp;
    const gsExprAssembler<> & _A; // Expression Assembler
    const gsExprAssembler<>::solution & _soln;
    mutable index_t m_pIndex;

public:
    /// Shared pointer for solutionFunction
    typedef memory::shared_ptr< solutionFunction > Ptr;

    /// Unique pointer for solutionFunction
    typedef memory::unique_ptr< solutionFunction > uPtr;

    solutionFunction(   const gsMultiPatch<> & mp,
                        const gsExprAssembler<> & A,
                     gsExprAssembler<>::solution & soln
                     ) : _mp(mp), _A(A), _soln(soln), _mm_piece(nullptr){ }

    GISMO_CLONE_FUNCTION(solutionFunction)

    short_t domainDim() const {return 2;}

    short_t targetDim() const {return 1;}

    mutable solutionFunction<T> * _mm_piece; // todo: improve the way pieces are accessed

    const gsFunction<T> & piece(const index_t p) const
    {
        // delete m_piece;
        // m_piece = new gsMaterialMatrix(m_patches->piece(k), *m_thickness, *m_YoungsModulus, *m_PoissonRatio);
        _mm_piece = new solutionFunction(*this);
        _mm_piece->setPatch(p);
        return *_mm_piece;
    }

    ~solutionFunction() { delete _mm_piece; }

    void setPatch(index_t p) {m_pIndex = p; }

    void eval_into(const gsMatrix<>& u, gsMatrix<>& res) const
    {
        // u is an d x n matrix of points with dimension d and n points
        // res2 is temporary vector for pressure components
        // res3 is vector with pressure and traction.
        // Data is stored as: [ tractionX1,tractionX2,tractionX3,...
        //                      tractionY1,tractionY2,tractionY3,...
        //                      pressure1, pressure2, pressure3 ,... ]

        gsExprEvaluator<T> _ev(_A);

        // Pre-allocate memory for pressure.
        res.resize(this->targetDim(), u.cols());
        res.setZero();
        // Create objects for parametric points (v) and for pressure-ONLY result (res2)
        gsVector <T> v;
        gsVector <index_t> patches;
        gsMatrix<T> points;
        _mp.locatePoints(u, patches, points);

        for( index_t k = 0; k != u.cols(); k++)
        {
            v = points.col(k);
            res.col(k) = _ev.eval(_soln, v, patches.at(k));
        }
    }
};


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

    // gsMultiPatch<> mp;
    // mp = *gsNurbsCreator<>::BSplineRectangle(-0.5,-0.5,0.5,0.5);

    // gsFunctionExpr<> f("0",2);

    // gsBoundaryConditions<> bc;
    // bc.addCondition(boundary::north, condition_type::dirichlet, 0 );
    // bc.addCondition(boundary::east, condition_type::dirichlet, 0 );
    // bc.addCondition(boundary::south, condition_type::dirichlet, 0 );
    // bc.addCondition(boundary::west, condition_type::dirichlet, 0 );

    gsMultiPatch<> mp;
    fd.getId(0, mp); // id=0: Multipatch domain

    gsFunctionExpr<> f;
    fd.getId(1, f); // id=1: source function
    gsInfo<<"Source function "<< f << "\n";

    gsBoundaryConditions<> bc;
    fd.getId(2, bc); // id=2: boundary conditions
    // bc.setMap(mp);
    gsInfo<<"Boundary conditions:\n"<< bc <<"\n";

    gsOptionList Aopt;
    fd.getId(4, Aopt); // id=4: assembler options

    //! [Read input file]

    //! [Refinement]
    gsMultiBasis<> basisL(mp);

    // Elevate and p-refine the basis to order p + numElevate
    // where p is the highest degree in the bases
    basisL.setDegree( basisL.maxCwiseDegree() + numElevate);

    // h-refine each basis
    for (int r =0; r < numRefine-1; ++r)
        basisL.uniformRefine();

    numRefine = 0;

    gsInfo<<"Basis Primal: "<<basisL.basis(0)<<"\n";


    gsInfo << "Patches: "<< mp.nPatches() <<", degree: "<< basisL.minCwiseDegree() <<"\n";
    //! [Refinement]

    //! [Problem setup]
    gsExprAssembler<> exL(1,1);
    exL.setOptions(Aopt);

    //gsInfo<<"Active options:\n"<< exL.options() <<"\n";
    typedef gsExprAssembler<>::geometryMap geometryMap;
    typedef gsExprAssembler<>::variable    variable;
    typedef gsExprAssembler<>::space       space;
    typedef gsExprAssembler<>::solution    solution;

    // Elements used for numerical integration
    exL.setIntegrationElements(basisL);
    gsExprEvaluator<> evL(exL);

    // Set the geometry map
    geometryMap G = exL.getMap(mp);

    // Set the discretization space
    space u = exL.getSpace(basisL);

    u.setInterfaceCont(0);
    u.addBc( bc.get("Dirichlet") );

    // Set the source term
    variable ff = exL.getCoeff(f, G);

    // Recover manufactured solution
    // gsFunctionExpr<> ms("((x-0.5)*(x+0.5)*(y-0.5)*(y+0.5)) / ((x+2/3)*y+3/4)",2);
    //gsInfo<<"Exact solution: "<< ms << "\n";
    gsFunctionExpr<> ms;
    fd.getId(3, ms); // id=3: reference solution
    variable u_exA = evL.getVariable(ms, G);

    // Solution vector and solution variable
    gsMatrix<> uVectorL, uVectorL2, uVectorH, zVectorH, phiVectorL, phiVectorH ;

    // Solutions and their projections
        // PDE solution on low-order mesh
    solution uL = exL.getSolution(u, uVectorL);

    gsSparseSolver<>::CGDiagonal solver;

    exL.initSystem();



    //! [Problem setup]

    //Treat labels: Dirichlet, CornerValues, Collapsed, Clamped
    // u.setup(bc.get("Dirichlet"), dirichlet::interpolation, 0); // def=-1
    //u.setupAsInteriorOnly(0); // def=-1

    // Initialize the system
    gsInfo<< "NumDofs Primal: "<<exL.numDofs() <<std::flush;

    // [ORIGINAL PROBLEM]
    // Compute the system matrix and right-hand side
    exL.assemble( igrad(u, G) * igrad(u, G).tr() * meas(G), u * ff * meas(G) );

    // Enforce Neumann conditions to right-hand side
    variable g_N = exL.getBdrFunction();
    exL.assembleRhsBc(u * g_N.val() * nv(G).norm(), bc.neumannSides() );
    //gsInfo<<"Sparse Matrix:\n"<< exL.matrix().toDense() <<"\n";
    //gsInfo<<"Rhs vector:\n"<< exL.rhs().transpose() <<"\n";

    // gsDebugVar(exL.matrix().toDense());
    // gsDebugVar(exL.rhs().transpose());

    gsInfo<< ".\n" <<std::flush;// Assemblying done

    solver.compute( exL.matrix() );
    uVectorL = solver.solve(exL.rhs());

    if (plot)
    {
        gsInfo<<"Plotting in Paraview...\n";
        // ev.options().setSwitch("plot.elements", true);
        evL.writeParaview( uL   , G, "solution");
    }

    gsInfo<<"Objective function errors J(u)-J(u_h)\n";
    // gsInfo<<"\t exact: "<<ev.integral(u_ex*meas(G))-ev.integral(u_sol*meas(G))<<"\n";
    gsInfo<<"\t exact: "<<evL.integral((u_exA-uL)*meas(G))<<"\n";
    // [!ORIGINAL PROBLEM]

    if (plot)
    {
        gsInfo<<"Plotting in Paraview...\n";
        // ev.options().setSwitch("plot.elements", true);
        evL.writeParaview( uL   , G, "solutionL");
    }

    // Let us construct a vector of coefficients including the boundaries.
    uL.extractFull(uVectorH);

    // Set the discretization space
    space v = exL.getSpace(basisL);
    exL.initSystem();

    // PDE solution projection on high-order mesh
    solution uH = exL.getSolution(v, uVectorH);


    gsDebugVar(uVectorH);
    gsDebugVar(uH.coefs());

    if (plot)
    {
        gsInfo<<"Plotting in Paraview...\n";
        // ev.options().setSwitch("plot.elements", true);
        evL.writeParaview( uH   , G, "solutionH");
    }

    return EXIT_SUCCESS;

}// end main
