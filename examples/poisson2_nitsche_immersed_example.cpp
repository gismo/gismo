/** @file poisson2_quarter_circle_immersed_example.cpp

    @brief Poisson on an immersed quarter-circle domain inside [0,1]^2.

    Embedding box: [0,1]x[0,1]
    Physical domain: Omega = {(x,y): x>0, y>0, x^2+y^2<1}
    ./build/bin/poisson2_nitsche_immersed_example -r 7 -e 1 --plot 

    This file is part of the G+Smo library.
*/

#include <gismo.h>

using namespace gismo;

int main(int argc, char *argv[])
{
    bool plot = false;
    index_t numRefine  = 5;
    index_t numElevate = 0;
    real_t penalty = 1e3;

    gsCmdLine cmd("Poisson immersed quarter-circle in [0,1]^2.");
    cmd.addInt("e", "degreeElevation", "Number of degree elevation steps", numElevate);
    cmd.addInt("r", "uniformRefine", "Number of uniform refinement steps (before solve)", numRefine);
    cmd.addReal("p", "penalty", "Weak Dirichlet penalty", penalty);
    cmd.addSwitch("plot", "Create a ParaView file", plot);
    try { cmd.getValues(argc, argv); } catch (int rv) { return rv; }

    // Embedding box [0,1]^2
    gsMultiPatch<> mp(*gsNurbsCreator<>::BSplineSquare());

    gsFunctionExpr<> u_exact("sin(pi*x)*sin(pi*y)", 2);
    gsFunctionExpr<> f_rhs("2*pi^2*sin(pi*x)*sin(pi*y)", 2);

    // Quarter-circle level-set: inside when max(-x,-y,x^2+y^2-1) < 0
    // gsFunctionExpr<> impl_fun("max(max(-x,-y),x^2+y^2-1)", 2);
    // gsFunctionExpr<> impl_fun("1 - x^2 - y^2", 2);
    gsFunctionExpr<> impl_fun("-0.16 + (x-0.5)^2 + (y-0.5)^2", 2);
    // Square
    gsFunctionExpr<> impl_square("x*(x-1)*y*(y-1)",2);
    // gsFunctionExpr<> normal_immersed("2*(x-0.5)","2*(y-0.5)",2);
    gsFunctionExpr<> normal_immersed_unit(
        "2*(x-0.5)/sqrt(4*(x-0.5)^2 + 4*(y-0.5)^2)",
        "2*(y-0.5)/sqrt(4*(x-0.5)^2 + 4*(y-0.5)^2)",
        2
    );


    gsBoundaryConditions<> bc;
    // bc.addCondition(0, boundary::west, condition_type::weak_dirichlet, &u_exact);
    // bc.addCondition(0, boundary::south, condition_type::weak_dirichlet, &u_exact);
    // bc.addCondition(0, boundary::none,  condition_type::weak_dirichlet, &u_exact); // immersed interface
    bc.setGeoMap(mp);

    gsMultiBasis<> dbasis(mp, true);
    dbasis.setDegree(dbasis.maxCwiseDegree() + numElevate);
    for (index_t i = 0; i < numRefine; ++i)
        dbasis.uniformRefine();

    gsTensorBSplineBasis<2,real_t> * tbsPtr =
        dynamic_cast<gsTensorBSplineBasis<2,real_t>*>(&dbasis.basis(0));
    GISMO_ENSURE(tbsPtr, "Basis is not a tensor B-spline basis");

    gsImplicitTrimmedDomain<2,real_t> tr_domain(impl_fun, *tbsPtr);

    // Debug export of trimmed-domain element classification:
    // sign < 0 : interior cell, sign == 0 : cut/boundary cell.
    gsMatrix<> boxesInterior(2,0), boxesBoundary(2,0), boxesAll(2,0);
    std::vector<real_t> classValues;
    classValues.reserve(tr_domain.numElements());

    auto appendBox = [](gsMatrix<> & boxes,
                        const gsVector<> & lo,
                        const gsVector<> & hi)
    {
        const index_t c = boxes.cols();
        boxes.conservativeResize(2, c + 2);
        boxes.col(c)     = lo;
        boxes.col(c + 1) = hi;
    };

    // internal elmnts
    for (auto it = tr_domain.beginInterior(); it != tr_domain.end<InteriorSign>(); ++it)
    {
        appendBox(boxesInterior, it.lowerCorner(), it.upperCorner());
        appendBox(boxesAll,      it.lowerCorner(), it.upperCorner());
        classValues.push_back(real_t(1));
    }
    // cutelements
    for (auto it = tr_domain.beginBdr(boundary::none); it != tr_domain.endBdr(boundary::none); ++it)
    {
        appendBox(boxesBoundary, it.lowerCorner(), it.upperCorner());
        appendBox(boxesAll,      it.lowerCorner(), it.upperCorner());
        classValues.push_back(real_t(2));
    }

    typedef gsExprAssembler<>::geometryMap geometryMap;
    typedef gsExprAssembler<>::space       space;
    typedef gsExprAssembler<>::solution    solution;

    gsExprAssembler<> A(1, 1);
    A.options().setInt("quRule", gsQuadrature::AlgoimRule);
    A.setIntegrationDomain(memory::make_shared_not_owned(&tr_domain));
    gsExprEvaluator<> ev(A);

    // auto begin = tr_domain.begin<BoundarySign>();
    // auto end = tr_domain.end<BoundarySign>();
    // auto begin = tr_domain.beginBdr(boundary::none);
    // auto end = tr_domain.endBdr(boundary::none);


    std::vector<patchSide> bdr_immersed(1); // vector of size 1
    bdr_immersed[0]= patchSide(0,boundary::none); // assign the first and only element
    // std::vector<patchSide> bdr_immersed = { patchSide(0, boundary::none) };

    // for (;begin!=end;begin++)
    // {
    //     gsInfo<< "upper corner of boundary element iterator: "<< begin.upperCorner() <<"\n";
    // }

    // for (gsBoxTopology::const_biterator it = bdr_immersed.begin();
    //      it != bdr_immersed.end(); ++it )
    // {
    //     auto begin = tr_domain.beginBdr(it->side());
    //     auto end = tr_domain.endBdr(it->side());
        
    //     for (;begin!=end;begin++)
    //     {
    //         gsInfo<< "upper corner of boundary element iterator: "<< begin.upperCorner() <<"\n";
    //     }
    // }

    // A.assembleBdr(bdr_immersed,expr);
    
    // Plot elements in the boundary (with container)
    // Verify that elements are correct! (!!)
    // Dont know the normal vector (normals known in the straight 

    // A.assembleBdr(bdr_immersed, expr);
    // Immerse a full circle 

    // Partially immersed partially non immersed

    // Problems 
    // gsFunctionExpr (normnal expression of hte circle)
    // Normal computation and boundary iteration 

    // Meshing library could give you the normal (CGAL)

    geometryMap G = A.getMap(mp);
    space u = A.getSpace(dbasis);

    auto ff    = A.getCoeff(f_rhs, G);
    auto n_imm = A.getCoeff(normal_immersed_unit, G); // grad(phi) — normalized?
    auto u_ex = ev.getVariable(u_exact, G);
    auto impl = A.getCoeff(impl_fun,G);
    gsMatrix<> solVector;
    solution u_sol = A.getSolution(u, solVector);

#ifdef GISMO_WITH_PARDISO
    gsSparseSolver<>::PardisoLDLT solver;
#else
    gsSparseSolver<>::CGDiagonal solver;
#endif

    gsInfo << "(dot1=assembled, dot2=solved, dot3=error)\n\nDoFs: ";

    // Determine maximum mesh size
    real_t hmax = 0;
    for (size_t p=0; p!=dbasis.nBases(); p++)
    {
        hmax = math::max(hmax, dbasis.basis(p).getMaxCellLength());
    }

    const index_t deg = dbasis.maxCwiseDegree();
    const real_t gamma = (penalty < 0)? real_t(2.5) * (deg + 2) * (deg + 1): penalty;

    u.setup(bc, dirichlet::none, 0);
    A.initSystem();
    A.computePattern(igrad(u) * igrad(u).tr());
    // A.getCoeff()

    A.assemble(
        igrad(u, G) * igrad(u, G).tr() * meas(G),
        u * ff * meas(G)
    );

    // auto g_D = A.getBdrFunction(G); //???? is it the immersed boundary? value? 0?
    auto g_D = A.getCoeff(u_exact, G); // prescribe solution at the immersed boundary!

    // A.assembleBdr(
    //     bc.get("Weak Dirichlet"),
    //     penalty * u * u.tr(),
    //     penalty * u * g_D
    // );

    // Symmetric Nitsche terms (matrix part)
    A.assembleBdr(
        bdr_immersed,
        - (igrad(u, G) * igrad(impl,G)) * u.tr()
        - u * (igrad(u, G) * igrad(impl,G)).tr()
        + gamma / hmax * u * u.tr() * meas(G)
    );

    // Symmetric Nitsche terms (rhs part)
    // unit outward normal: n_imm.normalized()
    A.assembleBdr(
        bdr_immersed,
        - (igrad(u, G) * igrad(impl,G)) * g_D * meas(G)
        + gamma / hmax * u * g_D * meas(G)
    );


    gsInfo << A.numDofs() << "." << std::flush;
    gsInfo<< "Matrix norm after assembly: " << A.matrix().norm() <<"\n";

    solver.compute(A.matrix());
    solVector = solver.solve(A.rhs());
    gsInfo << "." << std::flush;

    const real_t l2err = math::sqrt(ev.integral((u_ex - u_sol).sqNorm() * meas(G)));
    const real_t h1err = l2err +
        math::sqrt(ev.integral((igrad(u_ex) - igrad(u_sol, G)).sqNorm() * meas(G)));

    gsInfo << ".\n\nL2 error: " << std::scientific << std::setprecision(3) << l2err << "\n";
    gsInfo << "H1 error: " << std::scientific << h1err << "\n";

    if (plot)
    {
        const std::string debugDir = "ParaviewOutput/trimmed_domain_debug";
        gsFileManager::mkdir(debugDir);
        const std::string base = debugDir + "/poisson2_quarter_circle_immersed_circle0_mesh";

        if (boxesInterior.cols() > 0)
            gsWriteParaview(boxesInterior, base + "_inside", real_t(1));

        if (boxesBoundary.cols() > 0)
            gsWriteParaview(boxesBoundary, base + "_boundarysign", real_t(2));

        A.setIntegrationDomain(dbasis.domain());
        ev.options().setSwitch("plot.elements", true);
        ev.options().setInt("plot.npts", 100000);
        ev.writeParaview(u_sol, G, "poisson2_quarter_circle_immersed_circle");
    }

    return EXIT_SUCCESS;
}
