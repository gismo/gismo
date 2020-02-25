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
    bool fullL2 = false;
    bool plot = false;
    index_t numRefine  = 5;
    index_t numElevate = 0;
    bool last = false;
    std::string fn("pde/poisson2d_bvp2.xml");

    gsCmdLine cmd("Tutorial on solving a Poisson problem.");
    cmd.addInt( "e", "degreeElevation",
                "Number of degree elevation steps to perform before solving (0: equalize degree in all directions)", numElevate );
    cmd.addInt( "r", "uniformRefine", "Number of Uniform h-refinement steps to perform before solving",  numRefine );
    cmd.addString( "f", "file", "Input XML file", fn );
    cmd.addSwitch("last", "Solve solely for the last level of h-refinement", last);
    cmd.addSwitch("full", "full L2 projection", fullL2);
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
    gsMultiBasis<> basisH(mp);

    // Elevate and p-refine the basis to order p + numElevate
    // where p is the highest degree in the bases
    basisL.setDegree( basisL.maxCwiseDegree() + numElevate);
    basisH.setDegree( basisH.maxCwiseDegree() + numElevate);

    // h-refine each basis
    for (int r =0; r < numRefine-1; ++r)
    {
        basisL.uniformRefine();
        basisH.uniformRefine();
    }
    basisH.degreeElevate(1);

    numRefine = 0;

    gsInfo<<"Basis Primal: "<<basisL.basis(0)<<"\n";
    gsInfo<<"Basis Dual:   "<<basisH.basis(0)<<"\n";


    gsInfo << "Patches: "<< mp.nPatches() <<", degree: "<< basisL.minCwiseDegree() <<"\n";
    //! [Refinement]

    //! [Problem setup]
    gsExprAssembler<> exL(1,1);
    gsExprAssembler<> exH(1,1);
    exL.setOptions(Aopt);
    exH.setOptions(Aopt);

    //gsInfo<<"Active options:\n"<< exL.options() <<"\n";
    typedef gsExprAssembler<>::geometryMap geometryMap;
    typedef gsExprAssembler<>::variable    variable;
    typedef gsExprAssembler<>::space       space;
    typedef gsExprAssembler<>::solution    solution;

    // Elements used for numerical integration
    exL.setIntegrationElements(basisL);
    gsExprEvaluator<> evL(exL);

    exH.setIntegrationElements(basisH);
    gsExprEvaluator<> evH(exH);

    // Set the geometry map
    geometryMap G = exL.getMap(mp);
    geometryMap H = exH.getMap(mp);

    // Set the discretization space
    space u = exL.getSpace(basisL);
    space v = exH.getSpace(basisH);

    u.setInterfaceCont(0);
    u.addBc( bc.get("Dirichlet") );

    v.setInterfaceCont(0);
    v.addBc( bc.get("Dirichlet") );

    // Set the source term
    variable ff = exL.getCoeff(f, G);
    variable gg = exH.getCoeff(f, H);

    // Recover manufactured solution
    // gsFunctionExpr<> ms("((x-0.5)*(x+0.5)*(y-0.5)*(y+0.5)) / ((x+2/3)*y+3/4)",2);
    //gsInfo<<"Exact solution: "<< ms << "\n";
    gsFunctionExpr<> msP;
    fd.getId(3, msP); // id=3: reference solution primal
    gsFunctionExpr<> msD;
    fd.getId(9, msD); // id=9: reference solution dual

    variable primal_exL = evL.getVariable(msP, G);
    variable dual_exL   = evL.getVariable(msD, G);
    variable primal_exH = evH.getVariable(msP, H);
    variable dual_exH   = evH.getVariable(msD, H);

    // Solution vector and solution variable
    gsMatrix<> primalL, primalH, primalLp, dualL, dualLp, dualH, phiVectorL, phiVectorH ;

    // Solutions and their projections
        // PDE solution on low-order mesh
    solution uL = exL.getSolution(u, primalL);
        // PDE solution projection on high-order mesh
    solution uH = exH.getSolution(v, primalH);

        // Dual solution on low-order mesh
    solution zL = exL.getSolution(u, dualL);
        // Dual solution projection on high-order mesh
    solution zH = exH.getSolution(v, dualH);

    solution phiL = exL.getSolution(u, phiVectorL);
    solution phiH = exH.getSolution(v, phiVectorH);

    gsSparseSolver<>::CGDiagonal solver;

    exL.initSystem();
    exH.initSystem();



    //! [Problem setup]

    //Treat labels: Dirichlet, CornerValues, Collapsed, Clamped
    // u.setup(bc.get("Dirichlet"), dirichlet::interpolation, 0); // def=-1
    //u.setupAsInteriorOnly(0); // def=-1

    // Initialize the system
    exL.initSystem();

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
    primalL = solver.solve(exL.rhs());

    if (plot)
    {
        gsInfo<<"Plotting in Paraview...\n";
        // ev.options().setSwitch("plot.elements", true);
        evL.writeParaview( uL   , G, "solution");
        evL.writeParaview( primal_exL   , G, "solution_exact");
    }

    gsInfo<<"Objective function errors J(u)-J(u_h)\n";
    // gsInfo<<"\t exact: "<<ev.integral(u_ex*meas(G))-ev.integral(u_sol*meas(G))<<"\n";
    gsInfo<<"\t exact: "<<evL.integral((primal_exL-uL)*meas(G))<<"\n";
    // [!ORIGINAL PROBLEM]




    // [Low-order Dual PROBLEM]
    // Compute the system matrix and right-hand side
    exL.initSystem();
    exL.assemble( igrad(u, G) * igrad(u, G).tr() * meas(G), u * meas(G) );

    // Enforce Neumann conditions to right-hand side
    // variable g_N = exL.getBdrFunction();
    // exL.assembleRhsBc(u * g_N.val() * nv(G).norm(), bc.neumannSides() );

    solver.compute( exL.matrix() );
    dualL = solver.solve(exL.rhs());

    if (plot)
    {
        gsInfo<<"Plotting in Paraview...\n";
        // ev.options().setSwitch("plot.elements", true);
        evL.writeParaview( zL   , G, "dualL");
        evL.writeParaview( dual_exL   , G, "dual_exact");

    }

    // [!Low-order Dual PROBLEM]



    // gsInfo<<"pos\tpatch\tdof\tfree\tbnd\tcpl\n";
    // for (index_t j=0; j!=2; j++)
    // for (index_t i=0; i!=u.mapper().patchSize(j); i++)
    // {
    //     gsInfo<<i<<"\t"<<j<<"\t"<<u.mapper().index(i,j)<<"\t"<<u.mapper().is_free(i,j)<<"\t"<<u.mapper().is_boundary(i,j)<<"\t"<<u.mapper().is_coupled(i,j)<<"\n";
    // }

    // gsInfo<<"pos\tpatch\tdof\tfree\tbnd\tcpl\n";
    // for (index_t j=0; j!=2; j++)
    // for (index_t i=0; i!=v.mapper().patchSize(j); i++)
    // {
    //     gsInfo<<i<<"\t"<<j<<"\t"<<v.mapper().index(i,j)<<"\t"<<v.mapper().is_free(i,j)<<"\t"<<v.mapper().is_boundary(i,j)<<"\t"<<v.mapper().is_coupled(i,j)<<"\n";
    // }

    // Initialize the system
    exH.initSystem();

    gsInfo<< "NumDofs Dual: "<<exH.numDofs() <<std::flush;

    // [DUAL PROBLEM]
    // Compute the system matrix and right-hand side
    exH.assemble( igrad(v, H) * igrad(v, H).tr() * meas(H), v * meas(H) );

    gsInfo<< ".\n" <<std::flush;// Assemblying done

    solver.compute( exH.matrix() );
    dualH = solver.solve(exH.rhs());

    gsMatrix<> dualH0;
    space v0 = exH.getSpace(basisH);
    solution zH0 = exH.getSolution(v0,dualH0);
    if (fullL2)
    {
        zH.extractFull(dualH0);
        exH.initSystem();

        if (plot)
        {
            gsInfo<<"Plotting in Paraview...\n";
            // ev.options().setSwitch("plot.elements", true);
            evH.writeParaview( zH0   , H, "dualH");
        }
    }
    else
    {
        if (plot)
        {
            gsInfo<<"Plotting in Paraview...\n";
            // ev.options().setSwitch("plot.elements", true);
            evH.writeParaview( zH   , H, "dualH");
        }
    }
    // [!DUAL PROBLEM]

    // [Project the dual]
    // FROM strong dual basis (H) TO strong primal basis (L)
    gsExprAssembler<> exM(1,1);
    exM.setOptions(Aopt);
    geometryMap GM = exM.getMap(mp);

    space v_n = exM.getSpace(basisH);
    space u_n = exM.getTestSpace(v_n , basisL);


    if (!fullL2)
    {
        u_n.setInterfaceCont(0);
        u_n.addBc(bc.get("Dirichlet"));
        v_n.setInterfaceCont(0);
        v_n.addBc(bc.get("Dirichlet"));
    }
    exM.setIntegrationElements(basisH);
    exM.initSystem();
    exM.assemble(u_n * meas(GM)* v_n.tr()); // * meas(G));
    gsSparseMatrix<> M_HL = exM.matrix(); // from H to L
    gsSparseMatrix<> M_LH = M_HL.transpose(); // from L to H

    // exM.cleanUp();
    // geometryMap G2 = exM.getMap(mp);

    space w1_n = exM.getSpace(basisH);
    if (!fullL2)
    {
        w1_n.setInterfaceCont(0);
        w1_n.addBc(bc.get("Dirichlet"));
    }
    exM.setIntegrationElements(basisH);
    exM.initSystem();
    exM.assemble(w1_n * meas(GM)* w1_n.tr()); // * meas(G));
    gsSparseMatrix<> M_HH = exM.matrix();

    // exM.cleanUp();
    // geometryMap G3 = exM.getMap(mp);

    space w2_n = exM.getSpace(basisL);
    if (!fullL2)
    {
        w2_n.setInterfaceCont(0);
        w2_n.addBc(bc.get("Dirichlet"));
    }
    exM.setIntegrationElements(basisL);
    exM.initSystem();

    exM.assemble(w2_n * meas(GM)* w2_n.tr()); // * meas(G));
    gsSparseMatrix<> M_LL = exM.matrix();

    // gsDebugVar(M_HH.toDense());
    // gsDebugVar(M_LL.toDense());
    // gsDebugVar(M_LH.toDense());
    // gsDebugVar(M_HL.toDense());


    gsDebug<<"M_HH "<<M_HH.rows()<<"x"<<M_HH.cols()<<"\n";
    gsDebug<<"M_LL "<<M_LL.rows()<<"x"<<M_LL.cols()<<"\n";
    gsDebug<<"M_LH "<<M_LH.rows()<<"x"<<M_LH.cols()<<"\n";
    gsDebug<<"M_HL "<<M_HL.rows()<<"x"<<M_HL.cols()<<"\n";

    gsMatrix<> dualL0, primalL0;
    if (fullL2)
    {
        zL.extractFull(dualL0);
        uL.extractFull(primalL0);
    }

    solver.compute(M_HH);
    dualLp = solver.solve(M_LH*dualL0);

    solver.compute(M_HH);
    primalLp = solver.solve(M_LH*primalL0);

    if (fullL2)
    {
        solution zLp = exH.getSolution(v0, dualLp);
        solution uLp = exH.getSolution(v0, primalLp);
        if (plot)
        {
            gsInfo<<"Plotting in Paraview...\n";
            // ev.options().setSwitch("plot.elements", true);
            evH.writeParaview( zLp   , H, "dualLp");
            evH.writeParaview( uLp   , H, "primalLp");
        }
        gsVector<> pt(2);
        pt.setConstant(0.5);
        // gsInfo<<evH.eval( (igrad(zH, H)),pt)<<"\n\n";
        // gsInfo<<evH.eval( (grad(zH)*jac(H).ginv()),pt)<<"\n\n";
        // gsInfo<<evH.eval( igrad(uLp, H).tr(),pt )<<"\n\n";
        // gsInfo<<evH.eval( (igrad(zH, H) - igrad(zLp, H)) * igrad(uLp, H).tr() ,pt)<<"\n\n";
        // gsInfo<<evH.eval( (igrad(zH, H) - igrad(zLp, H)) * igrad(uLp, H).tr() * meas(H),pt )<<"\n\n";
        // gsInfo<<evH.eval( (zH-zLp) * ff * meas(H) ,pt)<<"\n\n";
        // // gsInfo<<evH.integral( (igrad(zH, G) - igrad(zLp, G)) * igrad(uLp, G).tr()  )<<"\n";
        // gsInfo<<evH.integral( grad(uLp)*jac(H).ginv() * (grad(uLp)*jac(H).ginv()).tr()  )<<"\n";
        // gsInfo<<evH.integral( (grad(zH)*jac(H).ginv() - grad(zLp)*jac(H).ginv()) * (grad(zH)*jac(H).ginv() - grad(zLp)*jac(H).ginv()).tr()  )<<"\n";
        // // gsInfo<<evH.integral( (zH-zLp)  )<<"\n";

        // variable g_H = exH.getBdrFunction();

        gsInfo<<"Objective function errors J(u)-J(u_h)\n";

        real_t error = evL.integral((primal_exL)*meas(G)) - evL.integral((uL)*meas(G));
        real_t errest = evH.integral( (zH-zLp) * gg * meas(H)-(((igrad(zH) - igrad(zLp))*igrad(uLp).tr()) ) * meas(H));
        // gsInfo<<"\t exact: "<<ev.integral(u_ex*meas(G))-ev.integral(u_sol*meas(G))<<"\n";
        gsInfo<<"\texact:\t"   <<error<<"\n";

        // gsInfo<<evH.integral( (zH-zLp) * gg * meas(H)-((grad(zH) - grad(zLp)) * (grad(uLp)).tr()) * meas(H))<<"\n";
        gsInfo<<"\testimate:\t"<<errest<<"\n";

        gsInfo<<"\teff:\t"<<errest/error<<"\n";
                    // - evH.integralBdr(uLp * g_H.val() * nv(H).norm())<<"\n"; //, bc.neumannSides()

        // gsInfo<<evH.eval(jac(H).ginv().tr()*( hess(uLp) - summ(grad(uLp)*jac(H).ginv(),hess(H)) ) * jac(H).ginv(),pt);
        // gsInfo<<evH.eval(ilapl,pt);
        // gsInfo<<evH.integralElWise( gg - lapl*meas(H))<<"\n";

        gsInfo<<"Computation errors (L2-norm)\n";
        // gsInfo<<"\t exact: "<<ev.integral(u_ex*meas(G))-ev.integral(u_sol*meas(G))<<"\n";
        real_t err1, err2;
        err1 = math::sqrt(evL.integral((primal_exL-uL).sqNorm()*meas(G)));
        err2 = math::sqrt(evH.integral((primal_exH-uLp).sqNorm()*meas(H)));
        gsInfo<<"\t (1) primal exact: "<<err1<<"\n";
        gsInfo<<"\t (2) primal proje: "<<err2<<"\n";
        gsInfo<<"\t ((1)-(2)/(1)    : "<<(err1-err2)/err1<<"\n";

        err1 = math::sqrt(evL.integral((dual_exL-zL).sqNorm()*meas(G)));
        err2 = math::sqrt(evH.integral((dual_exH-zLp).sqNorm()*meas(H)));
        gsInfo<<"\t (1) dual (L) exact: "<<err1<<"\n";
        gsInfo<<"\t (2) dual (L) proje: "<<err2<<"\n";
        gsInfo<<"\t ((1)-(2)/(1)      : "<<(err1-err2)/err1<<"\n";

        gsInfo<<"\t dual (H) exact: "<<math::sqrt(evH.integral((dual_exH-zH).sqNorm()*meas(H)))<<"\n";



        gsInfo<<evH.integralElWise((ff - lapl(uLp)) * (zH - zLp)*meas(H))<<"\n";
        // gsInfo<<evH.eval( hess(uLp),pt)<<"\n";
        // gsInfo<<evH.eval( grad(uLp),pt)<<"\n";
        // gsInfo<<evH.eval( lapl(uLp),pt)<<"\n";
        // gsInfo<<evH.eval( zH,pt)<<"\n";
        gsInfo<<evH.eval( (ff - lapl(uLp)) * (zH - zLp)*meas(H),pt)<<"\n";


    }
    else
    {
        solution zLp = exH.getSolution(v, dualLp);
        solution uLp = exH.getSolution(v, primalLp);
        if (plot)
        {
            gsInfo<<"Plotting in Paraview...\n";
            // ev.options().setSwitch("plot.elements", true);
            evH.writeParaview( zLp   , H, "dualLp");
            evH.writeParaview( uLp   , H, "primalLp");
        }
        gsVector<> pt(2);
        pt.setConstant(0.5);
        // gsInfo<<evH.eval(zH-zLp,pt)<<"\n";
        gsInfo<<evH.integral( -((grad(zH)*jac(H).ginv() - grad(zLp)*jac(H).ginv()) * (grad(uLp)*jac(H).ginv()).tr()) * meas(H) + (zH-zLp) * gg * meas(H))<<"\n";
    }


    // // [L2 PROJECTION]
    // solutionFunction<real_t> solDual(mp,D,zH);
    // variable ud = exL.getCoeff(solDual, G);

    // gsVector<> pt(2);
    // pt.setConstant(0.5);

    // exL.initSystem();
    // exL.assemble(u*u.tr(),u*uexH.val());

    // solver.compute( exL.matrix() );
    // // gsDebugVar(exL.matrix().toDense());
    // // gsDebugVar(exL.rhs().transpose());
    // projVector = solver.solve(exL.rhs());

    // if (plot)
    // {
    //     gsInfo<<"Plotting in Paraview...\n";
    //     // ev.options().setSwitch("plot.elements", true);
    //     evL.writeParaview( u_dual_proj  , G, "dual_proj");
    // }
    // solutionFunction<real_t> solDualP(mp,A,u_dual_proj);
    // variable udp = exL.getCoeff(solDualP, G);

    // gsInfo<<exL.integral(igrad(ud-udp, G) * igrad(u_sol, G).tr() * meas(G) - (ud-udp) * ff * meas(G))<<"\n";
    // gsInfo<<exL.integral((igrad(ud,G)-igrad(udp, G)) * igrad(u_sol, G).tr() * meas(G) - (ud-udp) * ff * meas(G))<<"\n";

    // gsInfo<<exL.eval(u_dual_proj,pt)<<"\n";
    // gsInfo<<exH.eval(u_dual,pt)<<"\n";

    // gsInfo<<exL.eval(igrad(u_dual_proj,G),pt)<<"\n";
    // gsInfo<<exH.eval(igrad(u_dual,G),pt)<<"\n";

    // //! [Export visualization in ParaView]
    // if (plot)
    // {
    //     gsInfo<<"Plotting in Paraview...\n";
    //     // ev.options().setSwitch("plot.elements", true);
    //     evH.writeParaview( u_exD   , H, "solution_ex");
    //     //ev.writeParaview( u, G, "aa");
    // }
    // else
    //     gsInfo << "Done. No output created, re-run with --plot to get a ParaView "
    //               "file containing the solution.\n";
    // //! [Export visualization in ParaView]

    return EXIT_SUCCESS;

}// end main
