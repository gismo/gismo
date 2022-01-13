/** @file pMultigrid_example.cpp

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
using std::max;
using std::min;

gsPreconditionerOp<>::Ptr setupSubspaceCorrectedMassSmoother(
    const gsSparseMatrix<>&,
    const gsMultiBasis<>&,
    const gsBoundaryConditions<>&,
    const gsOptionList&,
    const index_t&
);

gsPreconditionerOp<>::Ptr setupBlockILUT(
    const gsSparseMatrix<>&,
    const gsMultiPatch<>&,
    const gsMultiBasis<>&
);


/** @brief The p-multigrid class implements a generic p-multigrid solver
 *  that can be customized by passing assembler and coarse
 *  solver as template arguments.
 *
 *  @note: This implementation assumes that all required prolongation/
 *  restriction operators are generated internally. Therefore, a
 *  problem-specific assembler has to be passed as template argument.
 */
template<class Assembler>
struct pMultigrid
{
private:
    /// Shared pointer to multi-patch geometry
    memory::shared_ptr<gsMultiPatch<>> m_mp_ptr;

    /// Shared pointer to boundary conditions
    memory::shared_ptr<gsBoundaryConditions<>> m_bcInfo_ptr;

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
    std::vector<gsSparseMatrix<>> m_matrices;

    /// Vector of assembler objects
    std::vector<Assembler> m_assembler;

    double m_Time_Assembly, m_Time_Transfer, m_Time_Assembly_Galerkin;
public:

    ///  @brief Set-up p-multigrid solver
    gsMultiGridOp<>::uPtr setup(
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
        m_mp_ptr = memory::make_shared_not_owned(&mp);
        m_bcInfo_ptr = memory::make_shared_not_owned(&bcInfo);
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

        // Generate sequence of assembler objects and assemble
        for (std::vector<memory::shared_ptr<gsMultiBasis<>>>::iterator it = m_basis.begin();
        it != m_basis.end(); ++it)
        {
            m_assembler.push_back(Assembler(*m_mp_ptr,*(*it).get(),*m_bcInfo_ptr,rhs, coeff_diff , coeff_conv , coeff_reac ,(typeBCHandling == 1 ? dirichlet::elimination : dirichlet::nitsche), iFace::glue));
        }

        // Resize vectors of operators
        m_matrices.resize(numLevels);
        m_prolongation_H.resize(numLevels-1);
        m_restriction_H.resize(numLevels-1);
        m_restriction.resize(numLevels-1);
        m_prolongation.resize(numLevels-1);

        // Assemble operators at finest level
        gsStopwatch clock;
        gsInfo << "|| Multigrid hierarchy ||\n";
        for (index_t i = 0; i < numLevels; i++)
        {
            gsInfo << "Level " << i+1 << " " ;
            if (typeCoarseOperator == 1)
            {
                m_assembler[i].assemble();
                m_matrices[i] = m_assembler[i].matrix();
                gsInfo << "Degree: " << m_basis[i]->degree()  << ", Ndof: " << m_basis[i]->totalSize() << "\n";
            }
            else
            {
                if (hp(min(i,hp.rows()-1),0) == 0 || i == numLevels-1)
                {
                    m_assembler[i].assemble();
                    m_matrices[i] = m_assembler[i].matrix();
                    gsInfo << "\nDegree of the basis: " << m_basis[i]->degree() << "\n";
                    gsInfo << "Size of the basis functions: " << m_basis[i]->totalSize() << "\n";
                }
            }
        }
        m_Time_Assembly = clock.stop();

        // Determine prolongation/restriction operators
        clock.restart();
        gsOptionList options;
        options.addInt("DirichletStrategy","",typeBCHandling == 1 ? dirichlet::elimination : dirichlet::nitsche);
        for (index_t i = 1; i < numLevels; i++)
        {
            if (hp(i-1,0) == 0)
            {
                if (typeLumping == 1)
                {
                    gsSparseMatrix<> prolongationP = prolongation_P(mp, *m_basis[i], *m_basis[i-1], bcInfo, typeBCHandling);
                    gsSparseMatrix<> restrictionP = prolongationP.transpose();
                    gsMatrix<> prolongationM = prolongation_M(mp, *m_basis[i], bcInfo, typeBCHandling);
                    gsMatrix<> restrictionM = restriction_M(mp, *m_basis[i-1], bcInfo, typeBCHandling);

                    m_prolongation[i-1] = makeLinearOp(
                        [=](const gsMatrix<>& Xcoarse, gsMatrix<>& Xfine)
                        {
                            gsVector<> temp = prolongationP*Xcoarse;
                            gsMatrix<> Minv = prolongationM.array().inverse();
                            Xfine = Minv.cwiseProduct(temp);
                        },
                        prolongationP.rows(), prolongationP.cols()
                    );
                    m_restriction[i-1] = makeLinearOp(
                        [=](const gsMatrix<>& Xfine, gsMatrix<>& Xcoarse)
                        {
                            gsVector<> temp = restrictionP*Xfine;
                            gsMatrix<> Minv = restrictionM.array().inverse();
                            Xcoarse = Minv.cwiseProduct(temp);
                        },
                        restrictionP.rows(), restrictionP.cols()
                    );
                }
                else
                {
                    gsSparseMatrix<> prolongationP =  prolongation_P(mp, *m_basis[i], *m_basis[i-1], bcInfo, typeBCHandling);
                    gsSparseMatrix<> restrictionP =  prolongationP.transpose();
                    gsSparseMatrix<> prolongationM2 = prolongation_M2(mp, *m_basis[i], bcInfo, typeBCHandling);
                    gsSparseMatrix<> restrictionM2 = restriction_M2(mp, *m_basis[i-1], bcInfo, typeBCHandling);

                    m_prolongation[i-1] = makeLinearOp(
                        [=](const gsMatrix<>& Xcoarse, gsMatrix<>& Xfine)
                        {
                            gsMatrix<> temp = prolongationP*Xcoarse;
                            gsConjugateGradient<> CGSolver(prolongationM2);
                            CGSolver.setTolerance(1e-12);
                            CGSolver.solve(temp,Xfine);
                        },
                        prolongationP.rows(), prolongationP.cols()
                    );
                    m_restriction[i-1] = makeLinearOp(
                        [=](const gsMatrix<>& Xfine, gsMatrix<>& Xcoarse)
                        {
                            gsMatrix<> temp = restrictionP*Xfine;
                            gsConjugateGradient<> CGSolver(restrictionM2);
                            CGSolver.setTolerance(1e-12);
                            CGSolver.solve(temp,Xcoarse);
                        },
                        restrictionP.rows(), restrictionP.cols()
                    );

                }
            }
            else if (hp(i-1,0) == 1)
            {
                gsMultiBasis<> basis_copy = *m_basis[i];
                basis_copy.uniformCoarsen_withTransfer(m_prolongation_H[i-1],*m_bcInfo_ptr,options);
                m_restriction_H[i-1] = m_prolongation_H[i-1].transpose();

                m_prolongation[i-1] = makeMatrixOp(m_prolongation_H[i-1]);
                m_restriction[i-1] = makeMatrixOp(m_restriction_H[i-1]);

            }
        }
        m_Time_Transfer = clock.stop();

        // Obtain operators with Galerkin projection
        clock.restart();
        if (typeCoarseOperator == 2)
        {
            for (index_t i = numLevels-1; i > -1; i--)
            {
                if (hp(hp.rows()-1,0) == 0)
                {
                    if (hp(min(i,hp.rows()-1),0) == 1)
                    {
                        m_matrices[i] = m_restriction_H[i]*m_matrices[i+1]*m_prolongation_H[i];
                    }
                }
                else
                {
                    if (hp(min(i,hp.rows()-1),0) == 1 && i > 0)
                    {
                        m_matrices[i-1] = m_restriction_H[i-1]*m_matrices[i]*m_prolongation_H[i-1];
                    }
                }
            }
        }
        m_Time_Assembly_Galerkin = clock.stop();

        std::vector<gsLinearOperator<>::Ptr> operators(m_matrices.size());
        for (size_t i=0; i<m_matrices.size(); i++)
        {
            operators[i] = makeMatrixOp(m_matrices[i].moveToPtr());
        }

        return gsMultiGridOp<>::make(operators, m_prolongation, m_restriction);
    }

private:

    /// @brief Construct prolongation operator at level numLevels
    static gsMatrix<> prolongation_M(const gsMultiPatch<>& mp, const gsMultiBasis<>& basisH, const gsBoundaryConditions<>& bcInfo, index_t typeBCHandling)
    {
        typedef gsExprAssembler<real_t>::geometryMap geometryMap;
        typedef gsExprAssembler<real_t>::variable variable;
        typedef gsExprAssembler<real_t>::space space;

        // Determine matrix M (high_order * high_order)
        gsExprAssembler<real_t> ex(1,1);
        geometryMap G = ex.getMap(mp);
        space w_n = ex.getSpace(basisH, 1, 0);
        w_n.setInterfaceCont(0);
        if (typeBCHandling == 1)
        {
            w_n.setup(bcInfo, dirichlet::interpolation, 0);
        }
        ex.setIntegrationElements(basisH);
        ex.initSystem();
        ex.assemble(w_n * meas(G));
        return ex.rhs();
    }

    static gsSparseMatrix<> prolongation_M2(const gsMultiPatch<>& mp, const gsMultiBasis<>& basisH, const gsBoundaryConditions<>& bcInfo, index_t typeBCHandling)
    {
        typedef gsExprAssembler<real_t>::geometryMap geometryMap;
        typedef gsExprAssembler<real_t>::variable variable;
        typedef gsExprAssembler<real_t>::space space;

        // Determine matrix M (high_order * high_order)
        gsExprAssembler<real_t> ex(1,1);
        geometryMap G = ex.getMap(mp);
        space w_n = ex.getSpace(basisH, 1, 0);
        if (typeBCHandling == 1)
        {
            w_n.setup(bcInfo, dirichlet::interpolation, 0);
        }
        ex.setIntegrationElements(basisH);
        ex.initSystem();
        ex.assemble(w_n * meas(G) * w_n.tr());
        return ex.matrix();
    }

    /// @brief Construct prolongation operator at level numLevels
    static gsSparseMatrix<> prolongation_P(const gsMultiPatch<>& mp, const gsMultiBasis<>& basisH, const gsMultiBasis<>& basisL, const gsBoundaryConditions<>& bcInfo, index_t typeBCHandling)
    {
        typedef gsExprAssembler<real_t>::geometryMap geometryMap;
        typedef gsExprAssembler<real_t>::variable variable;
        typedef gsExprAssembler<real_t>::space space;

        // Determine matrix P (high_order * low_order)
        gsExprAssembler<real_t> ex(1,1);
        geometryMap G = ex.getMap(mp);
        space v_n = ex.getSpace(basisH, 1, 0);
        space u_n = ex.getTestSpace(v_n, basisL);
        if (typeBCHandling == 1)
        {
            v_n.setup(bcInfo, dirichlet::interpolation, 0);
            u_n.setup(bcInfo, dirichlet::interpolation, 0);
        }
        ex.setIntegrationElements(basisH);
        ex.initSystem();
        ex.assemble(u_n * meas(G) * v_n.tr());
        return ex.matrix().transpose();
    }

    /// @brief Construct restriction operator at level numLevels
    static gsMatrix<> restriction_M(const gsMultiPatch<>& mp, const gsMultiBasis<>& basisL, const gsBoundaryConditions<>& bcInfo, index_t typeBCHandling)
    {
        typedef gsExprAssembler<real_t>::geometryMap geometryMap;
        typedef gsExprAssembler<real_t>::variable variable;
        typedef gsExprAssembler<real_t>::space    space;

        // Determine matrix M (low_order * low_order)
        gsExprAssembler<real_t> ex(1,1);
        geometryMap G = ex.getMap(mp);
        space w_n = ex.getSpace(basisL ,1, 0);
        if (typeBCHandling == 1)
        {
            w_n.setup(bcInfo, dirichlet::interpolation, 0);
        }
        ex.setIntegrationElements(basisL);
        ex.initSystem();
        ex.assemble(w_n * meas(G));
        return ex.rhs();
    }

    /// @brief Construct restriction operator at level numLevels
    static gsSparseMatrix<> restriction_M2(const gsMultiPatch<>& mp, const gsMultiBasis<>& basisL, const gsBoundaryConditions<>& bcInfo, index_t typeBCHandling)
    {
        typedef gsExprAssembler<real_t>::geometryMap geometryMap;
        typedef gsExprAssembler<real_t>::variable variable;
        typedef gsExprAssembler<real_t>::space space;

        // Determine matrix M (low_order * low_order)
        gsExprAssembler<real_t> ex(1,1);
        geometryMap G = ex.getMap(mp);
        space w_n = ex.getSpace(basisL, 1, 0);
        if (typeBCHandling == 1)
        {
            w_n.setup(bcInfo, dirichlet::interpolation, 0);
        }
        ex.setIntegrationElements(basisL);
        ex.initSystem();
        ex.assemble(w_n * meas(G) * w_n.tr());
        return ex.matrix();
    }

public:
    const gsMatrix<>& rhs(index_t level) const { return m_assembler[level].rhs(); }
    const gsSparseMatrix<>& matrix(index_t level) const { return m_assembler[level].matrix(); }
    const gsMultiBasis<>& basis(index_t level) const { return *(m_basis[level]); }

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
    index_t numPatches = 1;
    index_t typeSolver = 1;
    index_t typeCycle_p = 1;
    index_t typeCycle_h = 2;
    index_t typeBCHandling = 2;
    index_t typeLumping = 1;
    index_t typeProjection = 2;
    index_t typeSmoother = 1;
    index_t typeCoarseOperator = 1;
    std::string typeCoarsening = "h";

    // Command line argument parser
    gsCmdLine cmd("This programm solves the CDR-equation with a p-multigrid or h-multigrid method");

    // Add command line arguments
    cmd.addInt("p", "Degree", "Number of order elevation steps", numDegree);
    cmd.addInt("r", "Refinement", "Number of global refinements", numRefine);
    cmd.addInt("v", "Smoothing", "Number of pre/post smoothing steps", numSmoothing);
    cmd.addInt("l", "Levels", "Number of levels in multigrid method", numLevels);
    cmd.addInt("b", "Benchmark", "Number of the benchmark", numBenchmark);
    cmd.addInt("P", "Patches", "Number of patch splittings (1) no splitting, (2) split each patch into 2^d patches, (3) split each patch into 4^d patches, etc.", numPatches);
    cmd.addInt("s", "Solver", "Type of solver: (1) mg as stand-alone solver (2) BiCGStab prec. with mg (3) CG prec. with mg", typeSolver);
    cmd.addInt("m", "Cycle_p", "Type of cycle, eather V-cycle (1) or W-cycle (2)", typeCycle_p);
    cmd.addInt("M", "Cycle_h", "Type of cycle, eather V-cycle (1) or W-cycle (2)", typeCycle_h);
    cmd.addInt("d", "BCHandling", "Handles Dirichlet BC's by elimination (1) or Nitsche's method (2)", typeBCHandling);
    cmd.addInt("L", "Lumping", "Restriction and Prolongation performed with the lumped (1) or consistent (2) mass matrix", typeLumping);
    cmd.addInt("D", "Projection", "Direct projection on coarsest level (1) or via all other levels (2)", typeProjection);
    cmd.addInt("S", "Smoother", "Type of smoother: (1) ILUT (2) Gauss-Seidel (3) SCMS or (4) Block ILUT", typeSmoother);
    cmd.addInt("G", "CoarseOperator", "Type of coarse operator in h-multigrid: (1) Rediscretization (2) Galerkin Projection", typeCoarseOperator);
    cmd.addString("z", "Coarsening", "Expression that defines coarsening strategy", typeCoarsening);

    // Read parameters from command line
    try { cmd.getValues(argc,argv);  } catch (int rv) { return rv; }

    if (typeSolver < 1 || typeSolver > 3)
    {
        gsInfo << "Unknown solver chosen.\n";
        return -1;
    }

    if (typeBCHandling < 1 || typeBCHandling > 2)
    {
        gsInfo << "Unknown boundary condition handling type chosen.\n";
        return -1;
    }

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

    if (typeSmoother < 1 || typeSmoother > 4)
    {
        gsInfo << "Unknown smoother chosen.\n";
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
            return -1;
    }

    // Print information about benchmark
    gsInfo << "Exact solution: " << sol_exact << "\n";
    gsInfo << "Right hand side: " << rhs_exact << "\n";

    // Handle the uniform splitting
    for (index_t i=0; i<numPatches-1; ++i)
    {
        mp = mp.uniformSplit();
    }

    gsInfo << "Number of patches: " << mp.nPatches() << "\n\n";

    // Construct two bases (coarse and fine level)
    gsMultiBasis<> basisL(mp);
    gsMultiBasis<> basisH(mp);

    gsMatrix<index_t> hp = gsMatrix<index_t>::Zero(numLevels-1,1);

    // Read string from command line
    real_t numRefH = 0;
    real_t numRefP = 0;
    real_t numRefZ = 0;

    // Convert input string to array
    if (typeCoarsening.size() != static_cast<size_t>(numLevels-1))
    {
        gsInfo << "The string provided to --Coarsening should have length " << numLevels-1 << "\n";
        return -1;
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
        else
        {
            hp(i,0) = 2;
            numRefZ = numRefZ + 1;
        }
    }

    // Apply refinement in p for coarse level
    if (numRefP + numRefZ == numDegree)
    {
        basisL.degreeReduce(1);
    }
    else
    {
        basisL.degreeIncrease(numDegree-numRefP-numRefZ-1);
    }

    // Apply refinement in h for coarse and fine level
    for (index_t i = 0; i < numRefine - numRefH - numRefZ; ++i)
    {
        basisL.uniformRefine();
    }
    for (index_t i = 0; i < numRefine ; ++i)
    {
        basisH.uniformRefine();
    }

    // Apply refinement in p for fine level
    basisH.degreeIncrease(numDegree-1);

    // Define boundary conditions
    gsBoundaryConditions<> bcInfo;
    gsFunctionExpr<> zero("0", mp.geoDim());
    for (gsMultiPatch<>::const_biterator bit = mp.bBegin(); bit != mp.bEnd(); ++bit)
    {
        bcInfo.addCondition( *bit, condition_type::dirichlet, &sol_exact);
    }
    bcInfo.setGeoMap(mp);

    // Generate sequence of bases on all levels
    if (typeProjection == 1)
    {
        numLevels = numLevels - numDegree + 2;
    }

    const bool symmSmoothing = typeSolver == 3;

    // Setup of p-mg object
    pMultigrid< gsCDRAssembler<real_t> > My_MG;
    gsMultiGridOp<>::Ptr mg = My_MG.setup(
        mp,
        basisL,
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

    gsStopwatch clock;
    clock.restart();
    mg->setCoarseSolver(makeSparseLUSolver(mg->matrix(0)));
    double Time_Coarse_Solver_Setup = clock.stop();

    // Setup of smoothers
    clock.restart();
    for (index_t i = 0; i < numLevels; i++)
    {
        switch (typeSmoother)
        {
            case 1:
            {
                if (typeProjection == 2 || i == numLevels-1)
                    mg->setSmoother(i,makeIncompleteLUOp(mg->matrix(i)));
                else
                    mg->setSmoother(i,makeGaussSeidelOp(mg->matrix(i)));
                break;
            }
            case 2:
            {
                mg->setSmoother(i,makeGaussSeidelOp(mg->matrix(i)));
                break;
            }
            case 3:
            {
                gsOptionList opt;
                opt.addReal("Scaling","",0.12);
                mg->setSmoother(i,setupSubspaceCorrectedMassSmoother(mg->matrix(i), My_MG.basis(i), bcInfo, opt, typeBCHandling));
                break;
            }
            case 4:
            {
                if (typeProjection == 2 || i == numLevels-1)
                    mg->setSmoother(i,setupBlockILUT(mg->matrix(i), mp, My_MG.basis(i)));
                else
                    mg->setSmoother(i,makeGaussSeidelOp(mg->matrix(i)));
                break;
            }
            default:;
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
    gsInfo << "Total setup time: " << My_MG.TimeAssembly() + My_MG.TimeAssemblyGalerkin() + My_MG.TimeTransfer() + Time_Coarse_Solver_Setup + Time_Smoother_Setup << "\n";

    // Setup of solver
    gsIterativeSolver<>::Ptr solver;
    if (typeSolver == 1)
    {
        gsInfo << "\n|| Solver information ||\np-multigrid is applied as stand-alone solver\n";
        // The preconditioned gradient method is noting but applying p-multigrid as a stand-alone solver.
        // We have to set step size = 1 to deactivate automatic stepsize control.
        solver = gsGradientMethod<>::make(mg->underlyingOp(), mg, 1);
    }
    else if (typeSolver == 2)
    {
        gsInfo << "\n|| Solver information ||\nBiCGStab is applied as solver, p-multigrid as a preconditioner\n";
        solver = gsBiCgStab<>::make(mg->underlyingOp(), mg);
    }
    else if (typeSolver == 3)
    {
        gsInfo << "\n|| Solver information ||\nCG is applied as solver, p-multigrid as a preconditioner\n";
        solver = gsConjugateGradient<>::make(mg->underlyingOp(), mg);
    }

    gsMatrix<> x = gsMatrix<>::Random(mg->underlyingOp()->rows(),1);

    real_t maxIter = 20; // mg->underlyingOp()->rows();
    real_t tol = 1e-8;

    // Unfortunately, the stopping criterion is relative to the rhs not to the initial residual (not yet configurable)
    const real_t rhs_norm = My_MG.rhs(numLevels-1).norm();
    solver->setTolerance(tol * (My_MG.rhs(numLevels-1)-My_MG.matrix(numLevels-1)*x).norm() / rhs_norm );
    solver->setMaxIterations(maxIter);
    gsMatrix<> error_history;
    clock.restart();
    solver->solveDetailed( My_MG.rhs(numLevels-1), x, error_history );
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
    gsInfo << "Residual after solving: "  << (My_MG.rhs(numLevels-1)-My_MG.matrix(numLevels-1)*x).norm() << "\n";

    return 0;
}

// Create the subspace corrected mass smoother
gsPreconditionerOp<>::Ptr setupSubspaceCorrectedMassSmoother(
    const gsSparseMatrix<>& matrix,
    const gsMultiBasis<>& mb,
    const gsBoundaryConditions<>& bc,
    const gsOptionList& opt,
    const index_t &typeBCHandling
)
{
    const short_t dim = mb.topology().dim();

    // Setup dof mapper
    gsDofMapper dm;
    mb.getMapper(
       typeBCHandling == 1 ? (dirichlet::strategy)opt.askInt("DirichletStrategy",11) : (dirichlet::strategy)opt.askInt("DirichletStrategy",14),
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
    return gsPreconditionerFromOp<>::make(makeMatrixOp(matrix), gsAdditiveOp<>::make(transfers, ops));
}

class gsBlockILUT : public gsPreconditionerOp<real_t>
{
public:
    typedef memory::unique_ptr<gsBlockILUT> uPtr;

    gsBlockILUT( const gsSparseMatrix<>& A, const gsMultiPatch<>& mp, const gsMultiBasis<>& mb ) : m_A(A)
    {
        index_t numPatches = mp.nPatches();

        gsDofMapper dm;
        mb.getMapper(true,dm,false); // This is what partition would do; but what about the Nitsche case?
        dm.finalize();

        // Subdivide the overal dofs into
        //   a) the patch-local dofs (l=0,...,numPatches-1)
        //   b) the coupled dofs (l=numPatches)
        m_sizes.resize(numPatches+1);
        m_shifts.resize(numPatches+1);
        m_shifts[0] = 0;
        for (index_t k=0; k<numPatches; k++)
        {
            m_sizes[k] = dm.findFreeUncoupled(k).rows();
            m_shifts[k+1] = m_shifts[k] + m_sizes[k];
        }
        m_sizes[numPatches] = m_A.rows() - m_shifts[numPatches];

        // Vector of factorized operators
        m_P.resize(numPatches+1);
        // Vector of factorized operators
        m_Pinv.resize(numPatches+1);
        // Define the A_aprox matrix
        m_A_aprox = gsSparseMatrix<>(m_A.rows(),m_A.cols());

        gsSparseMatrix<> S = m_A.block(m_shifts[numPatches],m_shifts[numPatches],m_sizes[numPatches],m_sizes[numPatches]);

        for (index_t k=0 ; k<numPatches; k++)
        {
            // Diagonal entries
            const gsSparseMatrix<> block = m_A.block(m_shifts[k],m_shifts[k],m_sizes[k],m_sizes[k]);
            Eigen::IncompleteLUT<real_t> ilu;
            ilu.setFillfactor(1);
            ilu.compute(block);
            m_P[k] = ilu.fillReducingPermutation();
            m_Pinv[k] = ilu.inversePermutation();
            m_A_aprox.block(m_shifts[k],m_shifts[k],m_sizes[k],m_sizes[k]) = ilu.factors();

            // Schur complement and off-diagonal contributions
            if (m_sizes[numPatches]>0)
            {
                gsSparseMatrix<real_t,RowMajor> ddB = m_A.block(m_shifts[numPatches],m_shifts[k],m_sizes[numPatches],m_sizes[k]);
                gsSparseMatrix<real_t> ddC = m_A.block(m_shifts[k],m_shifts[numPatches],m_sizes[k],m_sizes[numPatches]);
                gsMatrix<> ddBtilde(m_sizes[k],m_sizes[numPatches]);
                gsMatrix<> ddCtilde(m_sizes[k],m_sizes[numPatches]);
                for (index_t j=0; j<m_sizes[numPatches]; j++)
                {
                    gsMatrix<> Brhs = ddB.row(j)*m_P[k];
                    gsMatrix<> Crhs = m_Pinv[k]*ddC.col(j);
                    ddBtilde.col(j) = ilu.factors().triangularView<Eigen::Upper>().transpose().solve(Brhs.transpose());
                    ddCtilde.col(j) = ilu.factors().triangularView<Eigen::UnitLower>().solve(Crhs);
                }
                S -= (ddBtilde.transpose()*ddCtilde).sparseView();

                m_A_aprox.block(m_shifts[k],m_shifts[numPatches],m_sizes[k],m_sizes[numPatches]) = ddCtilde;
                m_A_aprox.block(m_shifts[numPatches],m_shifts[k],m_sizes[numPatches],m_sizes[k]) = ddBtilde.transpose();
            }
        }

        // If there is no coupling (for example if there is only one patch), the work
        // is done. We should not try to compute the ilu factorization then.
        if (m_sizes[numPatches]==0)
          return;

        // Preform ILUT on the S-matrix!
        Eigen::IncompleteLUT<real_t> ilu;
        ilu.setFillfactor(1);
        ilu.compute(S);
        m_P[numPatches] = ilu.fillReducingPermutation();
        m_Pinv[numPatches] = ilu.inversePermutation();
        m_A_aprox.block(m_shifts[numPatches],m_shifts[numPatches],m_sizes[numPatches],m_sizes[numPatches]) = ilu.factors();
    }

    void step(const gsMatrix<>& rhs, gsMatrix<>& x) const
    {
        const index_t sz = m_shifts.size();
        gsMatrix<> d = rhs-m_A*x;
        for (index_t k=0; k<sz; ++k)
            d.block(m_shifts[k],0,m_sizes[k],1) = m_Pinv[k]*d.block(m_shifts[k],0,m_sizes[k],1);
        d = m_A_aprox.template triangularView<Eigen::UnitLower>().solve(d);
        d = m_A_aprox.template triangularView<Eigen::Upper>().solve(d);
        for (index_t k=0; k<sz-1; ++k)
            d.block(m_shifts[k],0,m_sizes[k],1) = m_P[k]*d.block(m_shifts[k],0,m_sizes[k],1);
        x += d;
    }

    void stepT(const gsMatrix<>& rhs, gsMatrix<>& x) const
    {
        const index_t sz = m_shifts.size();
        gsMatrix<> d = rhs-m_A*x;
        for (index_t k=0; k<sz; ++k)
            d.block(m_shifts[k],0,m_sizes[k],1) = m_Pinv[k]*d.block(m_shifts[k],0,m_sizes[k],1);
        d = m_A_aprox.template triangularView<Eigen::UnitLower>().solve(d);
        d = m_A_aprox.template triangularView<Eigen::Upper>().solve(d);
        for (index_t k=0; k<sz-1; ++k)
            d.block(m_shifts[k],0,m_sizes[k],1) = m_P[k]*d.block(m_shifts[k],0,m_sizes[k],1);
        x += d;
    }

    index_t rows() const { return m_A.rows(); }
    index_t cols() const { return m_A.cols(); }
    gsLinearOperator<real_t>::Ptr underlyingOp() const { return makeMatrixOp(m_A); }

private:
    /// Underlying matrix
    gsSparseMatrix<> m_A;

    /// Block operator object
    gsMatrix<> m_A_aprox;

    /// Vector of factorized operators
    std::vector< Eigen::PermutationMatrix<Dynamic,Dynamic,index_t> > m_P;

    /// Vector of factorized operators
    std::vector< Eigen::PermutationMatrix<Dynamic,Dynamic,index_t> > m_Pinv;

    std::vector<index_t> m_sizes;
    std::vector<index_t> m_shifts;


};

gsPreconditionerOp<>::Ptr setupBlockILUT(
    const gsSparseMatrix<>& op,
    const gsMultiPatch<>& mp,
    const gsMultiBasis<>& mb
)
{
    return gsBlockILUT::uPtr(new gsBlockILUT(op, mp, mb));
}
