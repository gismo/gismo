/** @file pMultiGrid_example.cpp

    @brief The p-multigrid shows use of p-multigrid methods

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): R. Tielen
*/

#include <gismo.h>
#include <string>

using namespace gismo;

gsPreconditionerOp<>::Ptr setupSubspaceCorrectedMassSmoother(
    const gsSparseMatrix<>&,
    const gsMultiBasis<>&,
    const gsBoundaryConditions<>&,
    const gsOptionList&
);

gsPreconditionerOp<>::Ptr setupBlockILUT(
    const gsSparseMatrix<>&,
    const gsMultiBasis<>&,
    const gsBoundaryConditions<>&,
    const gsOptionList&
);


/** @brief The p-multigrid class implements a generic p-multigrid solver
 *  that can be customized by passing assembler and coarse
 *  solver as template arguments.
 *
 *  @note: This implementation assumes that all required prolongation/
 *  restriction operators are generated internally. Therefore, a
 *  problem-specific assembler has to be passed as template argument.
 */
struct pMultigrid
{
private:

    /// Vector of multi-basis objects
    std::vector<memory::shared_ptr<gsMultiBasis<>>> m_basis;

    /// Vector of prolongation matrices for h-refinement
    std::vector<gsSparseMatrix<real_t,RowMajor>> m_prolongation_H;

    /// Vector of restriction matrices for h-refinement
    std::vector<gsSparseMatrix<>> m_restriction_H;

    /// Vector of restriction operators
    std::vector<gsLinearOperator<>::Ptr> m_restriction;

    /// Vector of prolongation operators
    std::vector<gsLinearOperator<>::Ptr> m_prolongation;

    /// Vector of stiffness matrices
    std::vector<memory::shared_ptr<gsSparseMatrix<>>> m_matrices;

    /// Vector of linear operators
    std::vector<gsLinearOperator<>::Ptr> m_operators;

    /// Right-hand-side for finest grid
    gsMatrix<> m_rhs;

    double m_Time_Assembly, m_Time_Transfer, m_Time_Assembly_Galerkin;
public:

    /// @brief Set-up
    pMultigrid(
         const gsMultiPatch<> & mp,
         const gsMultiBasis<> & basis,
         const gsBoundaryConditions<> & bcInfo,
         const gsFunctionExpr<> & rhs,
         index_t numLevels,
         index_t numDegree,
         index_t typeBCHandling,
         index_t typeLumping,
         const gsMatrix<index_t>& hp,
         index_t typeProjection,
         index_t typeCoarseOperator,
         const gsFunctionExpr<> coeff_diff,
         const gsFunctionExpr<> coeff_conv,
         const gsFunctionExpr<> coeff_reac
        )
    {
        m_basis.push_back(memory::make_shared_not_owned(&basis));

        for (index_t i = 1; i < numLevels; i++)
        {
            m_basis.push_back(give(m_basis.back()->clone()));
            switch ( hp(i-1,0) )
            {
                case 0 : (typeProjection == 1 ? m_basis.back()->degreeIncrease(numDegree-1) : m_basis.back()->degreeIncrease()); break;

                case 1 : m_basis.back()->uniformRefine(); break;

                case 2 : m_basis.back()->uniformRefine(); m_basis.back()->degreeIncrease(); break;
            }
        }

        // Resize vectors of operators
        m_prolongation_H.resize(numLevels-1);
        m_restriction_H.resize(numLevels-1);
        m_restriction.resize(numLevels-1);
        m_prolongation.resize(numLevels-1);

        // Assemble operators at finest level
        gsStopwatch clock;

        // Determine prolongation/restriction operators
        clock.restart();
        gsOptionList options;
        options.addInt("DirichletStrategy","",typeBCHandling == 1 ? dirichlet::elimination : dirichlet::nitsche);
        for (index_t i = 1; i < numLevels; i++)
        {
            if (hp(i-1,0) != 1)
            {
                if (typeLumping == 1)
                {
                    gsSparseMatrix<> mixedMass = assembleMixedMass(mp, *m_basis[i], *m_basis[i-1], bcInfo, typeBCHandling);
                    gsSparseMatrix<real_t> prolongation
                          = assembleLumpedMass(mp, *m_basis[i], bcInfo, typeBCHandling).asDiagonal().inverse()
                          * mixedMass;
                    gsSparseMatrix<real_t> restriction
                          = assembleLumpedMass(mp, *m_basis[i-1], bcInfo, typeBCHandling).asDiagonal().inverse()
                          * mixedMass.transpose();
                    m_prolongation[i-1] = makeMatrixOp(prolongation.moveToPtr());
                    m_restriction[i-1] = makeMatrixOp(restriction.moveToPtr());
                }
                else
                {
                    gsSparseMatrix<> prolongationP =  assembleMixedMass(mp, *m_basis[i], *m_basis[i-1], bcInfo, typeBCHandling);
                    gsSparseMatrix<> restrictionP =  prolongationP.transpose();
                    gsSparseMatrix<> prolongationM = assembleMass(mp, *m_basis[i], bcInfo, typeBCHandling);
                    gsSparseMatrix<> restrictionM = assembleMass(mp, *m_basis[i-1], bcInfo, typeBCHandling);

                    m_prolongation[i-1] = makeLinearOp(
                        [m_prolongationP=give(prolongationP),m_prolongationM=give(prolongationM)]
                        (const gsMatrix<>& Xcoarse, gsMatrix<>& Xfine)
                        {
                            gsConjugateGradient<> CGSolver(m_prolongationM);
                            CGSolver.setTolerance(1e-12);
                            CGSolver.solve(m_prolongationP*Xcoarse,Xfine);
                        },
                        prolongationP.rows(), prolongationP.cols()
                    );
                    m_restriction[i-1] = makeLinearOp(
                        [m_restrictionP=give(restrictionP),m_restrictionM=give(restrictionM)]
                        (const gsMatrix<>& Xfine, gsMatrix<>& Xcoarse)
                        {
                            gsConjugateGradient<> CGSolver(m_restrictionM);
                            CGSolver.setTolerance(1e-12);
                            CGSolver.solve(m_restrictionP*Xfine,Xcoarse);
                        },
                        restrictionP.rows(), restrictionP.cols()
                    );

                }
            }
            else //if (hp(i-1,0) == 1)
            {
                gsMultiBasis<> basis_copy = *m_basis[i];
                basis_copy.uniformCoarsen_withTransfer(m_prolongation_H[i-1],bcInfo,options);
                m_restriction_H[i-1] = m_prolongation_H[i-1].transpose();

                m_prolongation[i-1] = makeMatrixOp(m_prolongation_H[i-1]);
                m_restriction[i-1] = makeMatrixOp(m_restriction_H[i-1]);

            }
        }
        m_Time_Transfer = clock.stop();

        gsInfo << "|| Multigrid hierarchy ||\n";
        m_matrices.resize(numLevels);
        m_operators.resize(m_matrices.size());
        m_Time_Assembly = 0;
        m_Time_Assembly_Galerkin = 0;
        for (index_t i = numLevels-1; i >= 0; --i)
        {
            if (typeCoarseOperator == 1 || i == numLevels-1 || hp(i,0) != 1)
            {
                gsInfo << "Assemble on level " << i << ": ";
                clock.restart();
                gsCDRAssembler<real_t> assembler(
                    mp,
                    *(m_basis[i]),
                    bcInfo,
                    rhs,
                    coeff_diff, coeff_conv, coeff_reac,
                    (typeBCHandling == 1 ? dirichlet::elimination : dirichlet::nitsche),
                    iFace::glue
                );
                assembler.assemble();
                m_matrices[i] = gsSparseMatrix<>(assembler.matrix()).moveToPtr();
                if (i==numLevels-1)
                    m_rhs = give(assembler.rhs());
                m_Time_Assembly += clock.stop();
                gsInfo << "Degree: " << m_basis[i]->degree()  << ", Ndof: " << m_matrices[i]->rows() << "\n";
            }
            else
            {
                gsInfo << "Apply Galerkin projection on level " << i << ": ";
                clock.restart();
                m_matrices[i] = gsSparseMatrix<>(m_restriction_H[i] * *(m_matrices[i+1]) * m_prolongation_H[i]).moveToPtr();
                m_Time_Assembly_Galerkin += clock.stop();
                gsInfo << "Ndof: " << m_matrices[i]->rows() << "\n";
            }
            m_operators[i] = makeMatrixOp(m_matrices[i]);
        }

    }

private:

    /// @brief Determine \f$ \int_\Omega p_i dx \f$ with basis functions \f$ p_i \f$ as vector
    /// The entries of this vector are the row-sums of the mass matrix
    static gsMatrix<> assembleLumpedMass(
        const gsMultiPatch<>& mp,
        const gsMultiBasis<>& basis,
        const gsBoundaryConditions<>& bcInfo,
        index_t typeBCHandling
    )
    {
        typedef gsExprAssembler<real_t>::geometryMap geometryMap;
        typedef gsExprAssembler<real_t>::variable variable;
        typedef gsExprAssembler<real_t>::space space;

        gsExprAssembler<real_t> ex(1,1);
        geometryMap G = ex.getMap(mp);
        space w_n = ex.getSpace(basis, 1, 0);
        w_n.setInterfaceCont(0);
        if (typeBCHandling == 1)
        {
            w_n.setup(bcInfo, dirichlet::interpolation, 0);
        }
        ex.setIntegrationElements(basis);
        ex.initSystem();
        ex.assemble(w_n * meas(G));
        return ex.rhs();
    }

    /// @brief Determines the mass matrix
    static gsSparseMatrix<> assembleMass(
        const gsMultiPatch<>& mp,
        const gsMultiBasis<>& basis,
        const gsBoundaryConditions<>& bcInfo,
        index_t typeBCHandling
    )
    {
        typedef gsExprAssembler<real_t>::geometryMap geometryMap;
        typedef gsExprAssembler<real_t>::variable variable;
        typedef gsExprAssembler<real_t>::space space;

        gsExprAssembler<real_t> ex(1,1);
        geometryMap G = ex.getMap(mp);
        space w_n = ex.getSpace(basis, 1, 0);
        if (typeBCHandling == 1)
        {
            w_n.setup(bcInfo, dirichlet::interpolation, 0);
        }
        ex.setIntegrationElements(basis);
        ex.initSystem();
        ex.assemble(w_n * meas(G) * w_n.tr());
        return ex.matrix();
    }

    /// @brief Determines the mass matrix with different bases for trial and test functions
    /// The evaluation of the integrals is baed on the trial space (basisU)
    static gsSparseMatrix<> assembleMixedMass(
        const gsMultiPatch<>& mp,
        const gsMultiBasis<>& basisU,
        const gsMultiBasis<>& basisV,
        const gsBoundaryConditions<>& bcInfo,
        index_t typeBCHandling
    )
    {
        typedef gsExprAssembler<real_t>::geometryMap geometryMap;
        typedef gsExprAssembler<real_t>::variable variable;
        typedef gsExprAssembler<real_t>::space space;

        gsExprAssembler<real_t> ex(1,1);
        geometryMap G = ex.getMap(mp);
        space v_n = ex.getSpace(basisU, 1, 0);
        space u_n = ex.getTestSpace(v_n, basisV);
        if (typeBCHandling == 1)
        {
            v_n.setup(bcInfo, dirichlet::interpolation, 0);
            u_n.setup(bcInfo, dirichlet::interpolation, 0);
        }
        ex.setIntegrationElements(basisU);
        ex.initSystem();
        ex.assemble(u_n * meas(G) * v_n.tr());
        return ex.matrix().transpose();
    }

public:
    const gsMatrix<>& rhs() const { return m_rhs; }
    const gsSparseMatrix<>& matrix(index_t level) const { return *(m_matrices[level]); }
    const gsSparseMatrix<>& matrix() const { return *(m_matrices.back()); }
    const gsMultiBasis<>& basis(index_t level) const { return *(m_basis[level]); }

    const std::vector<gsLinearOperator<>::Ptr>& restriction() const { return m_restriction; }
    const std::vector<gsLinearOperator<>::Ptr>& prolongation() const { return m_prolongation; }
    const std::vector<gsLinearOperator<>::Ptr>& operators() const { return m_operators; }

    double TimeAssembly() const { return m_Time_Assembly; }
    double TimeAssemblyGalerkin() const { return m_Time_Assembly_Galerkin; }
    double TimeTransfer() const { return m_Time_Transfer; }

};

int main(int argc, char* argv[])
{
    index_t numDegree = 2;
    index_t numRefine = 6;
    index_t numSmoothing = 1;
    index_t numLevels = 2;
    index_t numBenchmark = 3;
    index_t numSplits = 0;
    index_t typeSolver = 1;
    index_t typeCycle_p = 1;
    index_t typeCycle_h = 2;
    index_t typeBCHandling = 2;
    index_t typeLumping = 1;
    index_t typeProjection = 2;
    index_t typeSmoother = 1;
    index_t typeCoarseOperator = 1;
    std::string typeCoarsening = "h";
    index_t maxIter = 20;
    real_t tol = 1e-8;
    real_t dampingSCMS = 1;

    // Command line argument parser
    gsCmdLine cmd("This programm solves the CDR-equation with a p-multigrid or h-multigrid method");

    // Add command line arguments
    cmd.addInt("p", "Degree", "Number of order elevation steps", numDegree);
    cmd.addInt("r", "Refinement", "Number of global refinements", numRefine);
    cmd.addInt("v", "Smoothing", "Number of pre/post smoothing steps", numSmoothing);
    cmd.addInt("l", "Levels", "Number of levels in multigrid method", numLevels);
    cmd.addInt("b", "Benchmark", "Number of the benchmark", numBenchmark);
    cmd.addInt("P", "PatchSplits", "Number of patch splittings: split each pach that many times into 2^d patches", numSplits);
    cmd.addInt("s", "Solver", "Type of solver: (1) mg as stand-alone solver (2) BiCGStab prec. with mg (3) CG prec. with mg", typeSolver);
    cmd.addInt("m", "Cycle_p", "Type of cycle where p or hp-refinement is applied: (1) V-cycle or (2) W-cycle", typeCycle_p);
    cmd.addInt("M", "Cycle_h", "Type of cycle where h-refinement is applied: (1) V-cycle or (2) W-cycle", typeCycle_h);
    cmd.addInt("d", "BCHandling", "Handles Dirichlet BC's by (1) elimination or (2) Nitsche's method", typeBCHandling);
    cmd.addInt("L", "Lumping", "Restriction and Prolongation performed with the (1) lumped or (2) consistent mass matrix", typeLumping);
    cmd.addInt("D", "Projection", "Direct projection on (1) coarsest level or (2) via all other levels", typeProjection);
    cmd.addInt("S", "Smoother", "Type of smoother: (1) ILUT (2) Gauss-Seidel (3) SCMS or (4) Block ILUT", typeSmoother);
    cmd.addInt("G", "CoarseOperator", "Type of coarse operator in h-multigrid: (1) Rediscretization (2) Galerkin Projection", typeCoarseOperator);
    cmd.addString("z", "Coarsening", "Expression that defines coarsening strategy", typeCoarsening);
    cmd.addInt("", "MaxIter", "Maximum number of iterations for iterative solver", maxIter);
    cmd.addReal("", "Tolerance", "Threshold for iterative solver", tol);
    cmd.addReal("", "DampingSCMS", "Damping for subspace corrected mass smoother (otherwise ignored)", dampingSCMS);

    // Read parameters from command line
    try { cmd.getValues(argc,argv);  } catch (int rv) { return rv; }

    gsOptionList opt;

    if (typeBCHandling < 1 || typeBCHandling > 2)
    {
        gsInfo << "Unknown boundary condition handling type chosen.\n";
        return -1;
    }
    opt.addInt("DirichletStrategy","",typeBCHandling == 1 ? dirichlet::elimination : dirichlet::nitsche);

    if (typeLumping < 1 || typeLumping > 2)
    {
        gsInfo << "Unknown lumping type chosen.\n";
        return -1;
    }

    if (typeProjection < 1 || typeProjection > 2)
    {
        gsInfo << "Unknown projection type chosen.\n";
        return -1;
    }

    if (typeCoarseOperator < 1 || typeCoarseOperator > 4)
    {
        gsInfo << "Unknown smoother chosen.\n";
        return -1;
    }

    // Initialize solution, rhs and geometry
    gsFunctionExpr<> sol_exact, rhs_exact, coeff_diff, coeff_conv, coeff_reac;
    gsMultiPatch<> mp;
    gsInfo << "|| Benchmark information ||\n";
    switch (numBenchmark)
    {
        case 1:
            gsInfo << "CDR-equation the unit square\n";
            mp = gsMultiPatch<>(*gsNurbsCreator<>::BSplineSquare(1.0, 0.0, 0.0));
            sol_exact = gsFunctionExpr<>("sin(pi*x)*sin(pi*y)",2);
            rhs_exact = gsFunctionExpr<>("1.2*pi*pi*sin(pi*x)*sin(pi*y)+0.9*pi*pi*sin(pi*x)*sin(pi*y)+0.7*pi*pi*cos(pi*x)*cos(pi*y) + 0.4*pi*pi*cos(pi*x)*cos(pi*y) +0.4*pi*cos(pi*x)*sin(pi*y)-0.2*pi*sin(pi*x)*cos(pi*y)+0.3*sin(pi*x)*sin(pi*y)", 2);
            coeff_diff = gsFunctionExpr<>("1.2","-0.7","-0.4","0.9",2);
            coeff_conv = gsFunctionExpr<>("0.4","-0.2",2);
            coeff_reac = gsFunctionExpr<>("0.3",2);
            break;

        case 2:
            gsInfo << "Poisson equation on the quarter annulus (1)\n";
            mp = gsMultiPatch<>(*gsNurbsCreator<>::BSplineFatQuarterAnnulus(1.0, 2.0));
            sol_exact = gsFunctionExpr<>( "-(x*x+y*y-1)*(x*x+y*y-4)*x*y*y", 2);
            rhs_exact = gsFunctionExpr<>( "2*x*(22*x*x*y*y+21*y*y*y*y-45*y*y+x*x*x*x-5*x*x+4)", 2);
            coeff_diff = gsFunctionExpr<>("1","0","0","1",2);
            coeff_conv = gsFunctionExpr<>("0","0",2);
            coeff_reac = gsFunctionExpr<>("0",2);
            break;

        case 3:
            gsInfo << "Poisson equation on the quarter annulus (2)\n";
            mp = gsMultiPatch<>(*gsNurbsCreator<>::BSplineFatQuarterAnnulus(1.0, 2.0));
            sol_exact = gsFunctionExpr<>( "(x^2+y^2-3*sqrt(x^2+y^2)+2)*sin(2*atan(y/x))", 2);
            rhs_exact = gsFunctionExpr<>( "(8-9*sqrt(x^2 + y^2))*sin(2*atan(y/x))/(x^2+y^2)", 2);
            coeff_diff = gsFunctionExpr<>("1","0","0","1",2);
            coeff_conv = gsFunctionExpr<>("0","0",2);
            coeff_reac = gsFunctionExpr<>("0",2);
            break;

        case 4:
            gsInfo << "Poisson equation on an L-shaped domain\n";
            mp = gsMultiPatch<>(*gsNurbsCreator<>::BSplineLShape_p1());
            sol_exact = gsFunctionExpr<>( "if ( y>0, ( (x^2+y^2)^(1.0/3.0) )*sin( (2*atan2(y,x) - pi)/3.0 ), ( (x^2+y^2)^(1.0/3.0) )*sin( (2*atan2(y,x) +3*pi)/3.0 ) )", 2);
            rhs_exact = gsFunctionExpr<>( "0", 2);
            coeff_diff = gsFunctionExpr<>("1","0","0","1",2);
            coeff_conv = gsFunctionExpr<>("0","0",2);
            coeff_reac = gsFunctionExpr<>("0",2);
            break;

        case 5:
            gsInfo << "Poisson equation on the unit cube\n";
            mp = gsMultiPatch<>(*gsNurbsCreator<>::BSplineCube(1));
            sol_exact = gsFunctionExpr<>( "sin(pi*x)*sin(pi*y)*sin(pi*z)", 3);
            rhs_exact = gsFunctionExpr<>( "(3*pi^2 )*sin(pi*x)*sin(pi*y)*sin(pi*z)", 3);
            coeff_diff = gsFunctionExpr<>("1","0","0","0","1","0","0","0","1",3);
            coeff_conv = gsFunctionExpr<>("0","0","0",3);
            coeff_reac = gsFunctionExpr<>("0",3);
            break;

        case 6:
            gsInfo << "Poisson's equation on Yeti footprint\n";
            mp = *static_cast<gsMultiPatch<>::uPtr>(gsReadFile<>("domain2d/yeti_mp2.xml"));
            sol_exact = gsFunctionExpr<>( "sin(5*pi*x)*sin(5*pi*y)", 2);
            rhs_exact = gsFunctionExpr<>( "(50*pi^2 )*sin(5*pi*x)*sin(5*pi*y)", 2);
            coeff_diff = gsFunctionExpr<>("1","0","0","1",2);
            coeff_conv = gsFunctionExpr<>("0","0",2);
            coeff_reac = gsFunctionExpr<>("0",2);
            break;

        default:
            gsInfo << "Unknown benchmark case chosen.\n";
            return EXIT_FAILURE;
    }

    // Print information about benchmark
    gsInfo << "Exact solution: " << sol_exact << "\n";
    gsInfo << "Right hand side: " << rhs_exact << "\n";

    // Handle the uniform splitting
    for (index_t i=0; i<numSplits; ++i)
        mp = mp.uniformSplit();

    gsInfo << "Number of patches: " << mp.nPatches() << "\n\n";

    // Multi basis
    gsMultiBasis<> basis(mp);

    gsMatrix<index_t> hp = gsMatrix<index_t>::Zero(numLevels-1,1);

    // Read string from command line
    index_t numRefH = 0;
    index_t numRefP = 0;
    index_t numRefZ = 0;

    // Convert input string to array
    if (typeCoarsening.size() != static_cast<size_t>(numLevels-1))
    {
        gsInfo << "The string provided to --Coarsening should have length " << numLevels-1 << "\n";
        return EXIT_FAILURE;
    }
    for ( index_t i = 0; i < numLevels-1 ; ++i)
    {
        if ( typeCoarsening[i] == 'h')
        {
            hp(i,0) = 1;
            numRefH = numRefH + 1;
        }
        else if ( typeCoarsening[i] == 'p')
        {
            hp(i,0) = 0;
            numRefP = numRefP + 1;
        }
        else if ( typeCoarsening[i] == 'z')
        {
            hp(i,0) = 2;
            numRefZ = numRefZ + 1;
        }
        else
        {
            gsInfo << "Unknown coarsening type given.\n";
            return EXIT_FAILURE;
        }
    }

    // Apply refinement in p for coarse level
    if (numRefP + numRefZ == numDegree)
    {
        basis.degreeReduce(1);
    }
    else
    {
        basis.degreeIncrease(numDegree-numRefP-numRefZ-1);
    }

    // Apply refinement in h for coarse and fine level
    for (index_t i = 0; i < numRefine - numRefH - numRefZ; ++i)
    {
        basis.uniformRefine();
    }

    // Define boundary conditions
    gsBoundaryConditions<> bcInfo;
    gsFunctionExpr<> zero("0", mp.geoDim());
    for (gsMultiPatch<>::const_biterator bit = mp.bBegin(); bit != mp.bEnd(); ++bit)
    {
        bcInfo.addCondition(*bit, condition_type::dirichlet, &sol_exact);
    }
    bcInfo.setGeoMap(mp);

    // Generate sequence of bases on all levels
    if (typeProjection == 1)
    {
        numLevels = numLevels - numDegree + 2;
    }

    const bool symmSmoothing = typeSolver == 3;

    // Setup of p-mg object
    pMultigrid My_MG(
        mp,
        basis,
        bcInfo,
        rhs_exact,
        numLevels,
        numDegree,
        typeBCHandling,
        typeLumping,
        hp,
        typeProjection,
        typeCoarseOperator,
        coeff_diff,
        coeff_conv,
        coeff_reac
    );

    gsMultiGridOp<>::Ptr mg = gsMultiGridOp<>::make(My_MG.operators(), My_MG.prolongation(), My_MG.restriction());

    gsStopwatch clock;
    clock.restart();
    mg->setCoarseSolver(makeSparseLUSolver(mg->matrix(0)));
    double Time_Coarse_Solver_Setup = clock.stop();

    // Setup of smoothers
    opt.addReal("Scaling","",0.12);           // only used for SCMS
    opt.addReal("Damping","",dampingSCMS);    // only used for SCMS
    clock.restart();
    for (index_t i = 0; i < numLevels; i++)
    {
        switch (typeSmoother)
        {
            case 1:
                if (typeProjection == 2 || i == numLevels-1)
                    mg->setSmoother(i,makeIncompleteLUOp(mg->matrix(i)));
                else
                    mg->setSmoother(i,makeGaussSeidelOp(mg->matrix(i)));
                break;
            case 2:
                mg->setSmoother(i,makeGaussSeidelOp(mg->matrix(i)));
                break;
            case 3:
                mg->setSmoother(i,setupSubspaceCorrectedMassSmoother(mg->matrix(i), My_MG.basis(i), bcInfo, opt));
                break;
            case 4:
                if (typeProjection == 2 || i == numLevels-1)
                    mg->setSmoother(i,setupBlockILUT(mg->matrix(i), My_MG.basis(i), bcInfo, opt));
                else
                    mg->setSmoother(i,makeGaussSeidelOp(mg->matrix(i)));
                break;
            default:
                gsInfo << "Unknown smoother chosen.\n";
                return EXIT_FAILURE;
        }
    }
    double Time_Smoother_Setup = clock.stop();

    mg->setNumPreSmooth(numSmoothing);
    mg->setNumPostSmooth(numSmoothing);
    mg->setSymmSmooth(symmSmoothing);
    // For coarsest level, we stick to 1 in any case (exact solver is never cycled)
    for (index_t i=1; i<numLevels-1; i++)
    {
        if (hp(i,0) == 0) // offset?
            mg->setNumCycles(i,typeCycle_p);
        else
            mg->setNumCycles(i,typeCycle_h);
    }

    gsInfo << "\n|| Setup Timings || \n";
    gsInfo << "Total Assembly time: " << My_MG.TimeAssembly() << "\n";
    gsInfo << "Total Assembly time (Galerkin): " << My_MG.TimeAssemblyGalerkin() << "\n";
    gsInfo << "Total transfer setup time: " << My_MG.TimeTransfer() << "\n";
    gsInfo << "Coarse solver setup time: " << Time_Coarse_Solver_Setup << "\n";
    gsInfo << "Smoother setup time: " << Time_Smoother_Setup << "\n";
    gsInfo << "Total setup time: " << My_MG.TimeAssembly() + My_MG.TimeAssemblyGalerkin() + My_MG.TimeTransfer() 
        + Time_Coarse_Solver_Setup + Time_Smoother_Setup << "\n";

    // Setup of iterative solver
    gsIterativeSolver<>::Ptr solver;
    switch (typeSolver)
    {
        case 1:
            gsInfo << "\n|| Solver information ||\np-multigrid is applied as stand-alone solver\n";
            // The preconditioned gradient method is noting but applying p-multigrid as a stand-alone solver.
            // We have to set step size = 1 to deactivate automatic stepsize control.
            solver = gsGradientMethod<>::make(mg->underlyingOp(), mg, 1);
            break;
        case 2:
            gsInfo << "\n|| Solver information ||\nBiCGStab is applied as solver, p-multigrid as a preconditioner\n";
            solver = gsBiCgStab<>::make(mg->underlyingOp(), mg);
            break;
        case 3:
            gsInfo << "\n|| Solver information ||\nCG is applied as solver, p-multigrid as a preconditioner\n";
            solver = gsConjugateGradient<>::make(mg->underlyingOp(), mg);
            break;
        default:
            gsInfo << "Unknown iterative solver chosen.\n";
            return EXIT_FAILURE;
    }

    gsMatrix<> x = gsMatrix<>::Random(mg->underlyingOp()->rows(),1);

    // Unfortunately, the stopping criterion is relative to the rhs not to the initial residual (not yet configurable)
    const real_t rhs_norm = My_MG.rhs().norm();
    solver->setTolerance(tol * (My_MG.rhs()-My_MG.matrix()*x).norm() / rhs_norm );
    solver->setMaxIterations(maxIter);
    gsMatrix<> error_history;
    clock.restart();
    solver->solveDetailed( My_MG.rhs(), x, error_history );
    real_t Time_Solve = clock.stop();

    for (index_t i=1; i<error_history.rows(); ++i)
    {
        gsInfo << std::right << std::setw(4) << i
               << "  |  Residual norm: "
               << std::setprecision(6) << std::left << std::setw(15) << (error_history(i,0) * rhs_norm)
               << "  reduction:  1 / "
               << std::setprecision(3) << (error_history(i-1,0)/error_history(i,0))
               << "\n";
    }
    if (solver->error() <= solver->tolerance())
        gsInfo << "Solver reached accuracy goal after " << solver->iterations() << " iterations.\n";
    else
        gsInfo << "Solver did not reach accuracy goal within " << solver->iterations() << " iterations.\n";
    gsInfo << "The iteration took " << Time_Solve << " seconds.\n";
    // Determine residual and l2 error
    gsInfo << "Residual after solving: "  << (My_MG.rhs()-My_MG.matrix()*x).norm() << "\n";

    return 0;
}

// Create the subspace corrected mass smoother
gsPreconditionerOp<>::Ptr setupSubspaceCorrectedMassSmoother(
    const gsSparseMatrix<>& matrix,
    const gsMultiBasis<>& mb,
    const gsBoundaryConditions<>& bc,
    const gsOptionList& opt
)
{
    const short_t dim = mb.topology().dim();

    // Setup dof mapper
    gsDofMapper dm;
    mb.getMapper(
       (dirichlet::strategy)opt.askInt("DirichletStrategy",11),
       (iFace    ::strategy)opt.askInt("InterfaceStrategy", 1),
       bc,
       dm,
       0
    );
    const index_t nTotalDofs = dm.freeSize();

    // Decompose the whole domain into components
    std::vector< std::vector<patchComponent> > components = mb.topology().allComponents(true);
    const index_t nr_components = components.size();

    // Setup Dirichlet boundary conditions
    gsBoundaryConditions<> dir_bc;
    for ( index_t ps=0; ps < 2*dim; ++ps )
        dir_bc.addCondition( 0, 1+ps, condition_type::dirichlet, NULL );

    // Setup transfer matrices and local preconditioners
    std::vector< gsSparseMatrix<real_t,RowMajor> > transfers;
    transfers.reserve(nr_components);
    std::vector< gsLinearOperator<>::Ptr > ops;
    ops.reserve(nr_components);

    for (index_t i=0; i<nr_components; ++i)
    {
        gsMatrix<index_t> indices;
        std::vector<gsBasis<>::uPtr> bases = mb.componentBasis_withIndices(components[i],dm,indices,true);
        index_t sz = indices.rows();
        gsSparseEntries<> se;
        se.reserve(sz);
        for (index_t i=0; i<sz; ++i)
            se.add(indices(i,0),i,real_t(1));
        gsSparseMatrix<real_t,RowMajor> transfer(nTotalDofs,sz);
        transfer.setFrom(se);
        if (sz>0)
        {
            if (bases[0]->dim() == dim)
            {
                GISMO_ASSERT ( bases.size() == 1, "Only one basis is expected for each patch." );
                ops.push_back(
                    gsPatchPreconditionersCreator<>::subspaceCorrectedMassSmootherOp(
                        *(bases[0]),
                        dir_bc,
                        gsOptionList(),
                        opt.getReal("Scaling")
                    )
                );
            }
            else
            {
                gsSparseMatrix<> mat = transfer.transpose() * matrix * transfer;
                ops.push_back( makeSparseCholeskySolver(mat) );
            }
            transfers.push_back(give(transfer));
        }
    }
    gsPreconditionerOp<>::uPtr result = gsPreconditionerFromOp<>::make(makeMatrixOp(matrix), gsAdditiveOp<>::make(transfers, ops));
    result->setOptions(opt);
    return result;
}

gsPreconditionerOp<>::Ptr setupBlockILUT(
    const gsSparseMatrix<>& A,
    const gsMultiBasis<>& mb,
    const gsBoundaryConditions<>& bc,
    const gsOptionList& opt
)
{
    const index_t nPatches = mb.nPieces();

    // Setup dof mapper
    gsDofMapper dm;
    mb.getMapper(
       (dirichlet::strategy)opt.askInt("DirichletStrategy",11),
       (iFace    ::strategy)opt.askInt("InterfaceStrategy", 1),
       bc,
       dm,
       0
    );
    // Subdivide the overal dofs into
    //   a) the patch-local dofs (l=0,...,nPatches-1)
    //   b) the coupled dofs (l=nPatches)
    std::vector<index_t> sizes(nPatches+1);
    std::vector<index_t> shifts(nPatches+1);
    shifts[0] = 0;
    for (index_t k=0; k<nPatches; k++)
    {
        sizes[k] = dm.findFreeUncoupled(k).rows();
        shifts[k+1] = shifts[k] + sizes[k];
    }
    sizes[nPatches] = A.rows() - shifts[nPatches];

    // Vector of factorized operators
    std::vector< Eigen::PermutationMatrix<Dynamic,Dynamic,index_t> > P(nPatches+1);
    // Vector of factorized operators
    std::vector< Eigen::PermutationMatrix<Dynamic,Dynamic,index_t> > Pinv(nPatches+1);
    // Define the A_aprox matrix
    // Inefficient since it does not use a sparse matrix
    gsMatrix<> A_aprox(A.rows(),A.cols());
    A_aprox.setZero();

    gsSparseMatrix<> S = A.block(shifts[nPatches],shifts[nPatches],sizes[nPatches],sizes[nPatches]);

    for (index_t k=0 ; k<nPatches; k++)
    {
        // Diagonal entries
        gsSparseMatrix<> block = A.block(shifts[k],shifts[k],sizes[k],sizes[k]);
        Eigen::IncompleteLUT<real_t> ilu;
        ilu.setFillfactor(1);
        ilu.compute(block);
        P[k] = ilu.fillReducingPermutation();
        Pinv[k] = ilu.inversePermutation();
        A_aprox.block(shifts[k],shifts[k],sizes[k],sizes[k]) = ilu.factors();

        // Schur complement and off-diagonal contributions
        if (sizes[nPatches]>0)
        {
            gsMatrix<> ddB = A.block(shifts[nPatches],shifts[k],sizes[nPatches],sizes[k]);
            gsMatrix<> ddC = A.block(shifts[k],shifts[nPatches],sizes[k],sizes[nPatches]);
            gsMatrix<> ddBtilde = ilu.factors().triangularView<Eigen::Upper>().transpose().solve((ddB*P[k]).transpose()).transpose();
            gsMatrix<> ddCtilde = ilu.factors().triangularView<Eigen::UnitLower>().solve(Pinv[k]*ddC);

            S -= (ddBtilde*ddCtilde).sparseView();
            A_aprox.block(shifts[k],shifts[nPatches],sizes[k],sizes[nPatches]) = ddCtilde;
            A_aprox.block(shifts[nPatches],shifts[k],sizes[nPatches],sizes[k]) = ddBtilde;
        }
    }

    // If there is no coupling (for example if there is only one patch), the work
    // is done. We should not try to compute the ilu factorization then.
    if (sizes[nPatches]>0)
    {
        // Preform ILUT on the S-matrix.
        Eigen::IncompleteLUT<real_t> ilu;
        ilu.setFillfactor(1);
        ilu.compute(S);
        P[nPatches] = ilu.fillReducingPermutation();
        Pinv[nPatches] = ilu.inversePermutation();
        A_aprox.block(shifts[nPatches],shifts[nPatches],sizes[nPatches],sizes[nPatches]) = ilu.factors();
    }

    return gsPreconditionerFromOp<>::make(
        makeMatrixOp(A),
        makeLinearOp(
            [m_A_aprox=give(A_aprox), m_P=give(P), m_Pinv=give(Pinv),m_shifts=give(shifts),m_sizes=give(sizes),
                m_nPatches=nPatches](const gsMatrix<>& rhs, gsMatrix<>& x)
            {
                x.setZero(rhs.rows(),rhs.cols());
                for (index_t k=0; k<=m_nPatches; ++k)
                    x.block(m_shifts[k],0,m_sizes[k],1) = m_Pinv[k]*rhs.block(m_shifts[k],0,m_sizes[k],1);
                x = m_A_aprox.template triangularView<Eigen::UnitLower>().solve(x);
                x = m_A_aprox.template triangularView<Eigen::Upper>().solve(x);
                for (index_t k=0; k<=m_nPatches; ++k)
                    x.block(m_shifts[k],0,m_sizes[k],1) = m_P[k]*x.block(m_shifts[k],0,m_sizes[k],1);
            },
            A.rows(),
            A.cols()
        )
    );
}
