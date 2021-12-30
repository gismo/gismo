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


/** @brief The p-multigrid base class provides the basic
 *  methods (smoothing, prolongation, restriction) for
 *  implementing p-multigrid methods
 */

template<class T>
struct pMultigridBase
{

public:

    /// @brief Apply p-multigrid solver to given right-hand side on level l
    virtual void solve(
        const gsMatrix<T> & rhs,
        std::vector<memory::shared_ptr<gsMultiBasis<T> > > m_basis,
        gsMatrix<T>& x,
        const int& numLevels,
        const int& numCoarsening,
        const int& numRefine,
        const int& numSmoothing,
        int& numCoarseCycles,
        const int& typeCycle_p,
        int& typeCycle_h,
        const int& typeSolver,
        const int& typeBCHandling,
        gsBoundaryConditions<T> bcInfo,
        gsMultiPatch<> mp,
        gsGeometry<>::Ptr geo,
        const int& typeLumping,
        const int& typeProjection,
        const int& typeSmoother,
        std::vector<gsSparseMatrix<T>>& m_prolongation_P,
        std::vector<gsSparseMatrix<T>>& m_restriction_P,
        std::vector<gsMatrix<T>>& m_prolongation_M,
        std::vector<gsMatrix<T>>& m_restriction_M,
        std::vector<gsSparseMatrix<T>>& m_prolongation_H,
        std::vector<gsSparseMatrix<T>>& m_restriction_H,
        const gsMatrix<>& hp
        )
    {
        if ( numLevels == 1)
        {
            solvecoarse(rhs, x, numLevels);
            return;
        }

        const index_t typeCycle = (hp(max(numLevels-2,0),0) == 0) ? typeCycle_p : typeCycle_h;

        gsMatrix<T> fineRes, coarseRes, fineCorr, coarseCorr;
        presmoothing(rhs, x, numLevels, numSmoothing);
        residual(fineRes, rhs, x, numLevels);
        restriction(fineRes, coarseRes, numLevels, numCoarsening, m_basis, typeLumping,
            typeBCHandling, bcInfo, mp, geo, typeProjection, m_prolongation_P, m_restriction_P,
            m_prolongation_M, m_restriction_M, m_prolongation_H, m_restriction_H, hp);
        coarseCorr.setZero(coarseRes.rows(),1);
        for ( index_t j = 0 ; j < typeCycle ; j++)
        {
            solve(coarseRes, m_basis, coarseCorr, numLevels-1, numCoarsening, numRefine, numSmoothing,
                numCoarseCycles, typeCycle_p, typeCycle_h, typeSolver, typeBCHandling, bcInfo, mp, geo,
                typeLumping, typeProjection, typeSmoother, m_prolongation_P, m_restriction_P,
                m_prolongation_M, m_restriction_M, m_prolongation_H, m_restriction_H, hp);
        }
        prolongation(coarseCorr, fineCorr, numLevels, numCoarsening, m_basis, typeLumping,
            typeBCHandling, bcInfo, mp, geo, typeProjection, m_prolongation_P, m_restriction_P,
            m_prolongation_M, m_restriction_M, m_prolongation_H, m_restriction_H, hp);

        const real_t alpha = 1;
        x -= alpha * fineCorr;
        postsmoothing(rhs, x, numLevels, numSmoothing, typeSolver);
    }

    /// @brief Setup p-multigrid to given linear system
    virtual void setup(const gsMatrix<T> & rhs,
        std::vector<memory::shared_ptr<gsMultiBasis<T> > > m_basis ,
        gsMatrix<T>& x,
        const int& numLevels,
        const int& numCoarsening,
        const int& numRefine,
        const int& numSmoothing,
        int& numCoarseCycles,
        const int& typeCycle_p,
        const int& typeCycle_h,
        const int& typeSolver,
        const int& typeBCHandling,
        gsBoundaryConditions<T> bcInfo,
        gsMultiPatch<> mp,
        gsGeometry<>::Ptr geo,
        const int& typeLumping,
        const int& typeProjection,
        const int& typeSmoother,
        std::vector<gsSparseMatrix<T>>& m_prolongation_P,
        std::vector<gsSparseMatrix<T>>& m_restriction_P,
        std::vector<gsMatrix<T>>& m_prolongation_M,
        std::vector<gsMatrix<T>>& m_restriction_M,
        std::vector<gsSparseMatrix<T>>& m_prolongation_H,
        std::vector<gsSparseMatrix<T>>& m_restriction_H,
        const gsMatrix<>& hp)
    {}

    /// @brief Get residual (pure virtual method)
    virtual void residual(gsMatrix<T>& res, const gsMatrix<T>& rhs, const gsMatrix<T>& x, const int& numLevels) = 0;

    /// @brief Apply fixed number of smoothing steps (pure virtual method)
    virtual void presmoothing(const gsMatrix<T>& rhs, gsMatrix<T>& x, const int& numLevels,
        const int& numSmoothing) = 0;

    /// @brief Apply fixed number of smoothing steps (pure virtual method)
    virtual void postsmoothing(const gsMatrix<T>& rhs, gsMatrix<T>& x, const int& numLevels,
        const int& numSmoothing, const int& typeSolver) = 0;

    /// @brief Apply coarse solver (pure virtual method)
    virtual void solvecoarse(const gsMatrix<T>& rhs, gsMatrix<T>& x, const int& numLevels) = 0;

    /// @brief Prolongate coarse space function to fine space
    virtual gsSparseMatrix<T> prolongation_P(const int& numLevels, std::vector<memory::shared_ptr<gsMultiBasis<T> > > m_basis,
        const int& typeLumping, const int& typeBCHandling, gsGeometry<>::Ptr geo, const int& typeProjection) = 0;

    /// @brief Prolongate coarse space function to fine space
    virtual gsSparseMatrix<T> restriction_P(const int& numLevels, std::vector<memory::shared_ptr<gsMultiBasis<T> > > m_basis,
        const int& typeLumping, const int& typeBCHandling, gsGeometry<>::Ptr geo, const int& typeProjection) = 0;

    /// @brief Prolongate coarse space function to fine space
    virtual gsMatrix<T> prolongation_M(const int& numLevels, std::vector<memory::shared_ptr<gsMultiBasis<T> > > m_basis,
        const int& typeLumping, const int& typeBCHandling, gsGeometry<>::Ptr geo, const int& typeProjection) = 0;

    /// @brief Prolongate coarse space function to fine space
    virtual gsMatrix<T> restriction_M(const int& numLevels, std::vector<memory::shared_ptr<gsMultiBasis<T> > > m_basis,
        const int& typeLumping, const int& typeBCHandling, gsGeometry<>::Ptr geo, const int& typeProjection) = 0;

    /// @brief Prolongate coarse space function to fine space
    virtual void prolongation(const gsMatrix<T>& Xcoarse, gsMatrix<T>& Xfine, const int& numLevels, const int& numCoarsening,
        std::vector<memory::shared_ptr<gsMultiBasis<T> > > m_basis, const int& typeLumping, const int& typeBCHandling,
        gsBoundaryConditions<T> bcInfo, gsMultiPatch<> mp, gsGeometry<>::Ptr geo, const int& typeProjection,
        std::vector<gsSparseMatrix<T>>& m_prolongation_P, std::vector<gsSparseMatrix<T>>& m_restriction_P,
        std::vector<gsMatrix<T>>& m_prolongation_M, std::vector<gsMatrix<T>>& m_restriction_M,
        std::vector<gsSparseMatrix<T>>& m_prolongation_H, std::vector<gsSparseMatrix<T>>& m_restriction_H,
        const gsMatrix<>& hp)
    {
        if (hp(numLevels-2,0) == 1)
        {
            Xfine = m_prolongation_H[numLevels-2]*Xcoarse;
        }
        else
        {
            if (typeLumping == 1)
            {
                gsVector<> temp = m_prolongation_P[numLevels-2]*Xcoarse;
                gsMatrix<> M_L_inv = (m_prolongation_M[numLevels-2]).array().inverse();
                Xfine = (M_L_inv).cwiseProduct(temp);
            }
            else
            {
                // Define the low and high order basis
                gsMultiBasis<> basisL = *m_basis[numLevels-2];
                gsMultiBasis<> basisH = *m_basis[numLevels-1];
                typedef gsExprAssembler<real_t>::geometryMap geometryMap;
                typedef gsExprAssembler<real_t>::variable variable;
                typedef gsExprAssembler<real_t>::space space;
                // Determine matrix M (high_order * high_order)
                gsExprAssembler<real_t> ex2(1,1);
                geometryMap G2 = ex2.getMap(mp);
                space w_n = ex2.getSpace(basisH ,1, 0);
                //w_n.setInterfaceCont(0);
                if (typeBCHandling == 1)
                {
                    w_n.setup(bcInfo, dirichlet::interpolation, 0);
                }
                ex2.setIntegrationElements(basisH);
                ex2.initSystem();
                ex2.assemble(w_n * meas(G2) * w_n.tr());

                // Prolongate Xcoarse to Xfine
                gsVector<> temp = m_prolongation_P[numLevels-2]*Xcoarse;
                gsSparseMatrix<> M = ex2.matrix();
                gsConjugateGradient<> CGSolver(M);
                CGSolver.setTolerance(1e-12);
                CGSolver.solve(temp,Xfine);
            }
        }
    }

    /// @brief Restrict fine space function to coarse space
    virtual void restriction(const gsMatrix<T>& Xfine, gsMatrix<T>& Xcoarse, const int& numLevels, const int& numCoarsening,
        std::vector<memory::shared_ptr<gsMultiBasis<T> > > m_basis, const int& typeLumping, const int& typeBCHandling,
        const gsBoundaryConditions<T> & bcInfo, gsMultiPatch<> mp, gsGeometry<>::Ptr geo, const int& typeProjection,
        std::vector<gsSparseMatrix<T>>& m_prolongation_P, std::vector<gsSparseMatrix<T>>& m_restriction_P,
        std::vector<gsMatrix<T>>& m_prolongation_M, std::vector<gsMatrix<T>>& m_restriction_M,
        std::vector<gsSparseMatrix<T>>& m_prolongation_H, std::vector<gsSparseMatrix<T>>& m_restriction_H, const gsMatrix<>& hp)
    {
        if (hp(numLevels-2,0) == 1)
        {
            Xcoarse = m_restriction_H[numLevels-2]*Xfine;
        }
        else
        {
            if (typeLumping == 1)
            {
                // Standard way
                gsVector<> temp = m_restriction_P[numLevels-2]*Xfine;
                gsMatrix<> M_L_inv = (m_restriction_M[numLevels-2]).array().inverse();
                Xcoarse = (M_L_inv).cwiseProduct(temp);
            }
            else
            {
                // Define the low and high order basis
                gsMultiBasis<> basisL = *m_basis[numLevels-2];
                gsMultiBasis<> basisH = *m_basis[numLevels-1];
                typedef gsExprAssembler<real_t>::geometryMap geometryMap;
                typedef gsExprAssembler<real_t>::variable variable;
                typedef gsExprAssembler<real_t>::space space;

                // Determine matrix M (low_order * low_order)
                gsExprAssembler<real_t> ex2(1,1);
                geometryMap G2 = ex2.getMap(mp);
                space w_n = ex2.getSpace(basisL ,1, 0);
                //w_n.setInterfaceCont(0);
                if (typeBCHandling == 1)
                {
                    w_n.setup(bcInfo, dirichlet::interpolation, 0);
                }
                ex2.setIntegrationElements(basisL);
                ex2.initSystem();
                ex2.assemble(w_n * meas(G2) * w_n.tr());

                // Restrict Xfine to Xcoarse
                gsMatrix<> temp = m_restriction_P[numLevels-2]*Xfine;
                gsSparseMatrix<> M = ex2.matrix();
                gsConjugateGradient<> CGSolver(M);
                CGSolver.setTolerance(1e-12);
                CGSolver.solve(temp,Xcoarse);
            }
        }
    }
};


/** @brief The p-multigrid class implements a generic p-multigrid solver
 *  that can be customized by passing assembler and coarse
 *  solver as template arguments.
 *
 *  @note: This implementation assumes that all required prolongation/
 *  restriction operators are generated internally. Therefore, a
 *  problem-specific assembler has to be passed as template argument.
 */
template<class T, class CoarseSolver, class Assembler>
struct pMultigrid : public pMultigridBase<T>
{
private:

    /// Base class type
    typedef pMultigridBase<T> Base;

    /// Shared pointer to multi-patch geometry
    memory::shared_ptr<gsMultiPatch<T> > m_mp_ptr;

    /// Shared pointer to boundary conditions
    memory::shared_ptr<gsBoundaryConditions<T> > m_bcInfo_ptr;

    /// Vector of multi-basis objects
    std::vector<memory::shared_ptr<gsMultiBasis<T> > > m_basis;

    /// Vector of prolongation operators
    std::vector< gsSparseMatrix<T> > m_prolongation_P;

    /// Vector of restriction operators
    std::vector< gsSparseMatrix<T> > m_restriction_P;

    /// Vector of prolongation operators
    std::vector< gsMatrix<T> > m_prolongation_M;

    /// Vector of restriction operators
    std::vector< gsMatrix<T> > m_restriction_M;

    /// Vector of prolongation operators
    std::vector< gsSparseMatrix<T> > m_prolongation_H;

    /// Vector of restriction operators
    std::vector< gsSparseMatrix<T> > m_restriction_H;

    /// Vector of smoother objects
    std::vector< gsPreconditionerOp<>::Ptr > m_smoother;

    /// Vector of operator objects
    std::vector< gsSparseMatrix<T> > m_operator;

    /// Vector of assembler objects
    std::vector<Assembler> m_assembler;

    /// The coarse solver
    CoarseSolver m_csolver;
public:

    // Constructor
    pMultigrid(const gsMultiPatch<T> & mp, const gsMultiBasis<T> & basis, const gsBoundaryConditions<T> & bcInfo)
    {
        m_mp_ptr = memory::make_shared_not_owned(&mp);
        m_bcInfo_ptr = memory::make_shared_not_owned(&bcInfo);
        m_basis.push_back(memory::make_shared_not_owned(&basis));
    }

public:

    ///  @brief Set-up p-multigrid solver
    void setup(const gsFunctionExpr<T> & rhs, const gsFunctionExpr<T> & sol_exact, gsMatrix<T>& x, const int& numSmoothing, gsMatrix<T> f,const int& typeSolver, int& iterTot, int& typeCycle_p, int& typeCycle_h, int numLevels, const int& numCoarsening, const int& numDegree, const int& numRefine, const int& numBenchmark, const int& typeMultigrid, const int& typeBCHandling, gsGeometry<>::Ptr geo, const int& typeLumping, const gsMatrix<>& hp, const int& typeProjection,const int& typeSmoother, const int& typeCoarseOperator, const gsFunctionExpr<> coeff_diff, const gsFunctionExpr<> coeff_conv,const gsFunctionExpr<> coeff_reac)
    {
        for (int i = 1; i < numLevels; i++)
        {
            m_basis.push_back(give(m_basis.back()->clone()));
            switch((int) hp(i-1,0) )
            {
                case 0 : (typeProjection == 1 ? m_basis.back()->degreeIncrease(numDegree-1) : m_basis.back()->degreeIncrease()); break;

                case 1 : m_basis.back()->uniformRefine();  break;

                case 2:  m_basis.back()->uniformRefine();
                        m_basis.back()->degreeIncrease(); break;
            }
        }

        // Generate sequence of assembler objects and assemble
        for (typename std::vector<memory::shared_ptr<gsMultiBasis<T> > >::iterator it = m_basis.begin();
        it != m_basis.end(); ++it)
        {
            m_assembler.push_back(Assembler(*m_mp_ptr,*(*it).get(),*m_bcInfo_ptr,rhs, coeff_diff , coeff_conv , coeff_reac ,(typeBCHandling == 1 ? dirichlet::elimination : dirichlet::nitsche), iFace::glue));
        }

        // Resize vectors of operators
        m_operator.resize(numLevels); m_prolongation_P.resize(numLevels-1); m_prolongation_M.resize(numLevels-1); m_prolongation_H.resize(numLevels-1); m_restriction_P.resize(numLevels-1); m_restriction_M.resize(numLevels-1); m_restriction_H.resize(numLevels-1);

        // Assemble operators at finest level
        gsStopwatch clock;
        gsInfo << "|| Multigrid hierarchy ||\n";
        for (int i = 0; i < numLevels; i++)
        {
            gsInfo << "Level " << i+1 << " " ;
            if (typeCoarseOperator == 1)
            {
                m_assembler[i].assemble();
                m_operator[i] = m_assembler[i].matrix();
                gsInfo << "Degree: " << m_basis[i]->degree()  << ", Ndof: " << m_basis[i]->totalSize() << "\n";
            }
            else
            {
                if (hp(min(i,hp.rows()-1),0) == 0 || i == numLevels-1)
                {
                m_assembler[i].assemble();
                m_operator[i] = m_assembler[i].matrix();
                gsInfo << "\nDegree of the basis: " << m_basis[i]->degree() << "\n";
                gsInfo << "Size of the basis functions: " << m_basis[i]->totalSize() << "\n";
                }
            }
        }
        real_t Time_Assembly = clock.stop();

        // Determine prolongation/restriction operators in p
        clock.restart();
        for (int i = 1; i < numLevels; i++)
        {
            if (hp(i-1,0) == 0)
            {
                m_prolongation_P[i-1] =  prolongation_P(i+1, m_basis, typeLumping, typeBCHandling, geo, typeProjection);
                m_restriction_P[i-1] =  m_prolongation_P[i-1].transpose(); //restriction_P(i+1, m_basis, typeLumping, typeBCHandling, geo, typeProjection);
                m_prolongation_M[i-1] =  prolongation_M(i+1, m_basis, typeLumping, typeBCHandling, geo, typeProjection);
                m_restriction_M[i-1] = restriction_M(i+1, m_basis, typeLumping, typeBCHandling, geo, typeProjection);
            }
        }

        // Determine prolongation/restriction operators in h
        gsSparseMatrix<real_t, RowMajor> transferMatrix;
        gsOptionList options;
        typeBCHandling == 1 ? options.addInt("DirichletStrategy","",dirichlet::elimination) : options.addInt("DirichletStrategy","",dirichlet::nitsche);
        for (int i = 1; i < numLevels; i++)
        {
            if (hp(i-1,0) == 1)
            {
                gsMultiBasis<T> m_basis_copy = *m_basis[i];
                m_basis_copy.uniformCoarsen_withTransfer(transferMatrix,*m_bcInfo_ptr,options);
                m_prolongation_H[i-1] = transferMatrix;
                m_restriction_H[i-1] = m_prolongation_H[i-1].transpose();
            }
        }
        real_t Time_Transfer = clock.stop();

        // Obtain operators with Galerkin projection
        clock.restart();
        if (typeCoarseOperator == 2)
        {
            for (int i = numLevels-1; i > -1; i--)
            {
                if (hp(hp.rows()-1,0) == 0)
                {
                    if (hp(min(i,hp.rows()-1),0) == 1)
                    {
                        m_operator[i] = m_restriction_H[i]*m_operator[i+1]*m_prolongation_H[i];
                    }
                }
                else
                {
                    if (hp(min(i,hp.rows()-1),0) == 1 && i > 0)
                    {
                        m_operator[i-1] = m_restriction_H[i-1]*m_operator[i]*m_prolongation_H[i-1];
                    }
                }
            }
        }
        real_t Time_Assembly_Galerkin = clock.stop();


        // Setting up the all the smoothers
        m_smoother.resize(numLevels);

        // Determine ILUT factorizations at each level
        clock.restart();

        if (typeSmoother == 1)
        {
            // Generate factorizations (ILUT)
            for (int i = 0; i < numLevels; i++)
            {
                if (typeProjection == 2 || i == numLevels-1)
                {
                    m_smoother[i] = makeIncompleteLUOp(m_operator[i]);
                }
                else
                {
                    // Use Gauss-Seidel on all the other levels
                    m_smoother[i] = makeGaussSeidelOp(m_operator[i]);
                }
            }
        }
        real_t Time_ILUT_Factorization = clock.stop();
        if (typeSmoother == 2)
        {
            for (int i = 0; i < numLevels; i++)
            {
                m_smoother[i] = makeGaussSeidelOp(m_operator[i]);
            }
        }
        clock.restart();
        if (typeSmoother == 3)
        {
            // Generate sequence of SCM smoothers
            gsOptionList opt;
            opt.addReal("Scaling","",0.12);
            for (int i = 0 ; i < numLevels ; i++)
            {
                m_smoother[i] = setupSubspaceCorrectedMassSmoother(m_operator[i], *m_basis[i], *m_bcInfo_ptr, opt, typeBCHandling);
            }
        }
        real_t Time_SCMS = clock.stop();

        clock.restart();
        if (typeSmoother == 5)
        {
            for (int i = 0; i < numLevels; i++)
            {
                if (typeProjection == 2 || i == numLevels-1)
                {
                    m_smoother[i] = setupBlockILUT(m_operator[i], *m_mp_ptr, *(m_basis[i]));
                }
                else
                {
                    // Use Gauss-Seidel on all the other levels
                    m_smoother[i] = makeGaussSeidelOp(m_operator[i]);
                }
            }
        }
        real_t Time_Block_ILUT_Factorization = clock.stop();

        clock.restart();
        m_csolver.analyzePattern(m_operator[0]);
        m_csolver.factorize(m_operator[0]);
        real_t Time_Coarse_Solver_Setup = clock.stop();

        gsInfo << "\n|| Setup Timings || \n";
        gsInfo << "Total Assembly time: " << Time_Assembly << "\n";
        gsInfo << "Total Assembly time (Galerkin): " << Time_Assembly_Galerkin << "\n";
        gsInfo << "Total ILUT factorization time: " << Time_ILUT_Factorization << "\n";
        gsInfo << "Total block ILUT factorization time: " << Time_Block_ILUT_Factorization << "\n";
        gsInfo << "Total SCMS time: " << Time_SCMS << "\n";
        gsInfo << "Total Coarse solver setup time: " << Time_Coarse_Solver_Setup << "\n";
        gsInfo << "Total setup time: " << Time_Assembly_Galerkin + Time_Assembly + Time_Transfer + Time_ILUT_Factorization + Time_SCMS + Time_Coarse_Solver_Setup << "\n";
    }

    ///  @brief Apply p-multigrid solver to given right-hand side on level l
    void solve(const gsFunctionExpr<T> & rhs, const gsFunctionExpr<T> & sol_exact, gsMatrix<T>& x, const int& numSmoothing, gsMatrix<T> f,const int& typeSolver, int& iterTot, int& typeCycle_p, int& typeCycle_h, int numLevels, const int& numCoarsening, const int& numDegree, const int& numRefine, const int& numBenchmark, const int& typeMultigrid, const int& typeBCHandling, gsGeometry<>::Ptr geo, const int& typeLumping, const gsMatrix<>& hp, const int& typeProjection,const int& typeSmoother, const int& typeCoarseOperator)
    {
        gsStopwatch clock;

        if (typeSolver == 1)
        {
            x = gsMatrix<>::Random(m_operator[numLevels-1].rows(),1);
        }

        gsMatrix<> b;
        typeSolver == 1 ? b = m_assembler.back().rhs() : b = f;


        // Determine residual and L2 error
        real_t r0 = (m_operator[numLevels-1]*x - b).norm();
        real_t r = r0;
        real_t tol = 1e-8;
        int iter = 1;
        int numCoarseCycles = 0;

        // Solve with p-multigrid method
        real_t r_old = r0;
        clock.restart();
        while( (typeSolver == 1 || typeSolver == 5) ? r/r0 > tol && iter < 100000 : iter < 2)
        {
            // Call solver from base class
            Base::solve(b, m_basis,  x, numLevels, numCoarsening, numRefine, numSmoothing, numCoarseCycles, typeCycle_p, typeCycle_h,
                typeSolver, typeBCHandling, *m_bcInfo_ptr, *m_mp_ptr, geo, typeLumping, typeProjection, typeSmoother, m_prolongation_P, m_restriction_P,
                m_prolongation_M, m_restriction_M, m_prolongation_H, m_restriction_H, hp);
            numCoarseCycles = 0;
            r = (m_operator[numLevels-1]*x - b).norm();
            if ( r_old < r)
            {
                gsInfo << "Residual increased during solving!!! \n";
            }
            r_old = r;
            //gsInfo << "Residual after cycle " << iter << " equals: " << r << "\n";
            iter++;
            iterTot++;
        }
        real_t Time_Solve = clock.stop();
        gsInfo << "\n|| Solver information || \n";
        gsInfo << "Solver converged in " << Time_Solve << " seconds!\n";
        gsInfo << "Solver converged in " << iter-1 << " iterations!\n";

        if (typeSolver == 1)
        {
            // Determine residual and L2 errpr
            gsInfo << "Residual after solving: "  << (b-m_operator[numLevels-1]*x).norm() << "\n";
        }
    }

private:

    /// @brief Apply coarse solver
    virtual void solvecoarse(const gsMatrix<T>& rhs, gsMatrix<T>& x, const int& numLevels)
    {
        x = m_csolver.solve(rhs);
    }

    /// @brief Construct prolongation operator at level numLevels
    virtual gsMatrix<T> prolongation_M(const int& numLevels, std::vector<memory::shared_ptr<gsMultiBasis<T> > > m_basis, const int& typeLumping, const int& typeBCHandling, gsGeometry<>::Ptr geo, const int& typeProjection)
    {
        // Define the low and high order basis
        gsMultiBasis<> basisL = *m_basis[numLevels-2];
        gsMultiBasis<> basisH = *m_basis[numLevels-1];

        // Determine matrix M (high_order * high_order)
        typedef gsExprAssembler<real_t>::geometryMap geometryMap;
        typedef gsExprAssembler<real_t>::variable variable;
        typedef gsExprAssembler<real_t>::space    space;
        gsExprAssembler<real_t> ex2(1,1);
        geometryMap G2 = ex2.getMap(*m_mp_ptr);
        space w_n = ex2.getSpace(basisH ,1, 0);
        w_n.setInterfaceCont(0);
        if (typeBCHandling == 1)
        {
            w_n.setup(*m_bcInfo_ptr, dirichlet::interpolation, 0);
        }
        ex2.setIntegrationElements(basisH);
        ex2.initSystem();
        ex2.assemble(w_n * meas(G2) );
        return ex2.rhs();
    }

    /// @brief Construct prolongation operator at level numLevels
    virtual gsSparseMatrix<T> prolongation_P(const int& numLevels, std::vector<memory::shared_ptr<gsMultiBasis<T> > > m_basis, const int& typeLumping, const int& typeBCHandling, gsGeometry<>::Ptr geo, const int& typeProjection)
    {
        // Define the low and high order basis
        gsMultiBasis<> basisL = *m_basis[numLevels-2];
        gsMultiBasis<> basisH = *m_basis[numLevels-1];

        // Determine matrix P (high_order * low_order)
        typedef gsExprAssembler<real_t>::geometryMap geometryMap;
        gsExprAssembler<real_t> ex(1,1);
        geometryMap G = ex.getMap(*m_mp_ptr);
        typedef gsExprAssembler<real_t>::variable variable;
        typedef gsExprAssembler<real_t>::space    space;
        space v_n = ex.getSpace(basisH ,1, 0);
        //v_n.setInterfaceCont(0);
        space u_n = ex.getTestSpace(v_n , basisL);
        //u_n.setInterfaceCont(0);
        if (typeBCHandling == 1)
        {
            v_n.setup(*m_bcInfo_ptr, dirichlet::interpolation, 0);
            u_n.setup(*m_bcInfo_ptr, dirichlet::interpolation, 0);
        }
        ex.setIntegrationElements(basisH);
        ex.initSystem();
        ex.assemble(u_n*meas(G) * v_n.tr());
        gsSparseMatrix<> P = ex.matrix().transpose();
        return P;
    }

    /// @brief Construct restriction operator at level numLevels
    virtual gsMatrix<T> restriction_M(const int& numLevels, std::vector<memory::shared_ptr<gsMultiBasis<T> > > m_basis, const int& typeLumping, const int& typeBCHandling, gsGeometry<>::Ptr geo, const int& typeProjection)
    {
        // Define the low and high order basis
        gsMultiBasis<> basisL = *m_basis[numLevels-2];
        gsMultiBasis<> basisH = *m_basis[numLevels-1];

        // Determine matrix M (low_order * low_order)
        typedef gsExprAssembler<real_t>::geometryMap geometryMap;
        typedef gsExprAssembler<real_t>::variable variable;
        typedef gsExprAssembler<real_t>::space    space;
        gsExprAssembler<real_t> ex2(1,1);
        geometryMap G2 = ex2.getMap(*m_mp_ptr);
        space w_n = ex2.getSpace(basisL ,1, 0);
        //w_n.setInterfaceCont(0);
        if (typeBCHandling == 1)
        {
            w_n.setup(*m_bcInfo_ptr, dirichlet::interpolation, 0);
        }
        ex2.setIntegrationElements(basisL);
        ex2.initSystem();
        ex2.assemble(w_n * meas(G2) );
        return ex2.rhs();
    }

    /// @brief Construct restriction operator at level numLevels
    virtual gsSparseMatrix<T> restriction_P(const int& numLevels, std::vector<memory::shared_ptr<gsMultiBasis<T> > > m_basis, const int& typeLumping, const int& typeBCHandling, gsGeometry<>::Ptr geo, const int& typeProjection)
    {
        // Define the low and high order basis
        gsMultiBasis<> basisL = *m_basis[numLevels-2];
        gsMultiBasis<> basisH = *m_basis[numLevels-1];

        // Determine matrix P (high_order * low_order)
        gsExprAssembler<real_t> ex(1,1);
        typedef gsExprAssembler<real_t>::geometryMap geometryMap;
        geometryMap G = ex.getMap(*m_mp_ptr);

        typedef gsExprAssembler<real_t>::variable variable;
        typedef gsExprAssembler<real_t>::space    space;
        space v_n = ex.getSpace(basisH ,1, 0);
        //v_n.setInterfaceCont(0);
        space u_n = ex.getTestSpace(v_n , basisL);
        //u_n.setInterfaceCont(0);
        if ( typeBCHandling == 1)
        {
            u_n.setup(*m_bcInfo_ptr, dirichlet::interpolation, 0);
            v_n.setup(*m_bcInfo_ptr, dirichlet::interpolation, 0);
        }
        ex.setIntegrationElements(basisH);
        ex.initSystem();
        ex.assemble(u_n * meas(G)* v_n.tr());
        gsSparseMatrix<> P = ex.matrix();
        return P;
    }

    /// @brief Get residual (pure virtual method)
    virtual void residual(gsMatrix<T>& res, const gsMatrix<T>& rhs, const gsMatrix<T>& x, const int& numLevels)
    {
        res = m_operator[numLevels-1]*x - rhs;
    }

    /// @brief Apply fixed number of presmoothing steps
    virtual void presmoothing(const gsMatrix<T>& rhs, gsMatrix<T>& x, const int& numLevels, const int& numSmoothing)
    {
        for (index_t i = 0 ; i < numSmoothing ; i++)
        {
            m_smoother[numLevels-1]->step(rhs,x);
        }
    }

    /// @brief Apply fixed number of postsmoothing steps
    virtual void postsmoothing(const gsMatrix<T>& rhs, gsMatrix<T>& x, const int& numLevels, const int& numSmoothing, const int& typeSolver)
    {
        for (index_t i = 0 ; i < numSmoothing ; i++)
        {
            if (typeSolver==3)
                m_smoother[numLevels-1]->stepT(rhs,x);
            else
                m_smoother[numLevels-1]->step(rhs,x);
        }
    }
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
    index_t typeMultigrid = 1;
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
    cmd.addInt("b", "Benchmark", "Number of the benchmark",numBenchmark);
    cmd.addInt("P", "Patches", "Number of patches (1, 4, 16 or 64)", numPatches);
    cmd.addInt("s", "Solver", "Type of solver: (1) mg as stand-alone solver (2) BiCGStab prec. with mg (3) CG prec. with mg", typeSolver);
    cmd.addInt("t", "Multigrid", "p-multigrid (1) or h-multigrid (2)", typeMultigrid);
    cmd.addInt("m", "Cycle_p", "Type of cycle, eather V-cycle (1) or W-cycle (2)", typeCycle_p);
    cmd.addInt("M", "Cycle_h", "Type of cycle, eather V-cycle (1) or W-cycle (2)", typeCycle_h);
    cmd.addInt("d", "BCHandling", "Handles Dirichlet BC's by elimination (1) or Nitsche's method (2)", typeBCHandling);
    cmd.addInt("L", "Lumping", "Restriction and Prolongation performed with the lumped (1) or consistent (2) mass matrix", typeLumping);
    cmd.addInt("D", "Projection", "Direct projection on coarsest level (1) or via all other levels (2)", typeProjection);
    cmd.addInt("S", "Smoother", "Type of smoother: (1) ILUT (2) Gauss-Seidel (3) SCMS or (5) Block Gauss-Seidel", typeSmoother);
    cmd.addInt("G", "CoarseOperator", "Type of coarse operator in h-multigrid: (1) Rediscretization (2) Galerkin Projection", typeCoarseOperator);
    cmd.addString("z", "Coarsening", "Expression that defines coarsening strategy", typeCoarsening);

    // Read parameters from command line
    try { cmd.getValues(argc,argv);  } catch (int rv) { return rv; }

    GISMO_ENSURE(typeSmoother == 1||typeSmoother == 2||typeSmoother == 3||typeSmoother == 5,
        "Unknown smoother chosen.");

    // Initialize solution, rhs and geometry
    std::string solution_exact,rhs_exact;
    gsGeometry<>::Ptr geo;
    gsInfo << "|| Benchmark information ||\n";
    switch(numBenchmark)
    {
      case 1 : gsInfo << "CDR-equation the unit square\n";
               solution_exact = "sin(pi*x)*sin(pi*y)";
               rhs_exact = "1.2*pi*pi*sin(pi*x)*sin(pi*y)+0.9*pi*pi*sin(pi*x)*sin(pi*y)+0.7*pi*pi*cos(pi*x)*cos(pi*y) + 0.4*pi*pi*cos(pi*x)*cos(pi*y) +0.4*pi*cos(pi*x)*sin(pi*y)-0.2*pi*sin(pi*x)*cos(pi*y)+0.3*sin(pi*x)*sin(pi*y)";
               geo = gsNurbsCreator<>::BSplineSquare(1.0, 0.0, 0.0); break;

      case 2:  gsInfo << "Poisson equation on the quarter annulus (1)\n";
               solution_exact = "-(x*x+y*y-1)*(x*x+y*y-4)*x*y*y";
               rhs_exact = "2*x*(22*x*x*y*y+21*y*y*y*y-45*y*y+x*x*x*x-5*x*x+4)";
               geo = gsNurbsCreator<>::BSplineFatQuarterAnnulus(1.0, 2.0); break;

      case 3:  gsInfo << "Poisson equation on the quarter annulus (2)\n";
               solution_exact = "(x^2+y^2-3*sqrt(x^2+y^2)+2)*sin(2*atan(y/x))";
               rhs_exact = "(8-9*sqrt(x^2 + y^2))*sin(2*atan(y/x))/(x^2+y^2)";
               geo = gsNurbsCreator<>::BSplineFatQuarterAnnulus(1.0, 2.0); break;

      case 4:  gsInfo << "Poisson equation on an L-shaped domain\n";
               solution_exact = "if ( y>0, ( (x^2+y^2)^(1.0/3.0) )*sin( (2*atan2(y,x) - pi)/3.0 ), ( (x^2+y^2)^(1.0/3.0) )*sin( (2*atan2(y,x) +3*pi)/3.0 ) )";
               rhs_exact = "0";
               geo = gsNurbsCreator<>::BSplineLShape_p1(); break;

      case 5:  gsInfo << "Poisson equation on the unit cube\n";
               solution_exact = "sin(pi*x)*sin(pi*y)*sin(pi*z)";
               rhs_exact = "(3*pi^2 )*sin(pi*x)*sin(pi*y)*sin(pi*z)";
               geo = gsNurbsCreator<>::BSplineCube(1); break;

      case 6: gsInfo << "Poisson's equation on Yeti footprint\n";
               solution_exact = "sin(5*pi*x)*sin(5*pi*y)";
               rhs_exact = "(50*pi^2 )*sin(5*pi*x)*sin(5*pi*y)";
               geo = gsReadFile<>("domain2d/yeti_mp2.xml"); break;
    }

    // Print information about benchmark
    gsInfo << "Exact solution: " << solution_exact << "\n";
    gsInfo << "Right hand side: " << rhs_exact << "\n";

    gsMultiPatch<> mp;
    switch(numPatches)
    {
        case 1:
        {
            gsMultiPatch<> mp0(*geo);
            mp = mp0;
            break;
        }
        case 2:
        {
            gsMultiPatch<> mp1(*geo);
            mp = mp1.uniformSplit();
            break;
        }
        case 3:
        {
            gsMultiPatch<> mp2(*geo);
            gsMultiPatch<> mp3 = mp2.uniformSplit();
            mp = mp3.uniformSplit();
            break;
        }
        case 4:
        {
            gsMultiPatch<> mp4(*geo);
            gsMultiPatch<> mp5 = mp4.uniformSplit();
            gsMultiPatch<> mp6 = mp5.uniformSplit();
            mp = mp6.uniformSplit();
            break;
        }
    }

    // To read the entire Yeti footprint
    if (numBenchmark == 10)
    {
        std::string geometry("domain2d/yeti_mp2.xml");
        gsMultiPatch<>::uPtr mpPtr = gsReadFile<>(geometry);
        mp = *mpPtr;
    }
    gsInfo << "Number of patches: " << mp.nPatches() << "\n\n";
    gsFunctionExpr<> sol_exact(solution_exact,mp.geoDim());
    gsFunctionExpr<> f(rhs_exact,mp.geoDim());

    // Construct two bases (coarse and fine level)
    gsMultiBasis<> basisL(mp);
    gsMultiBasis<> basisH(mp);

    // Geometry of benchmark 10 is describe by quadratic B-splines
    if (numBenchmark == 10)
    {
        basisL.degreeReduce();
        basisH.degreeReduce();
    }
    // Geometry of benchmark 7 is describe by quadratic B-splines
    if (numBenchmark == 7)
    {
        basisL.degreeReduce(numDegree-1);
        basisH.degreeReduce(numDegree-1);
    }
    gsMatrix<> hp = gsMatrix<>::Zero(numLevels-1,1);

    // Read string from command line
    real_t numRefH = 0;
    real_t numRefP = 0;
    real_t numRefZ = 0;

    // Convert input string to array
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

    index_t iterTot = 1;
    index_t numCoarsening = numRefH+1;

    pMultigrid<real_t, gsSparseSolver<real_t>::LU , gsCDRAssembler<real_t> > My_MG(mp, basisL, bcInfo);
    gsMatrix<real_t> x;

    // Generate sequence of bases on all levels
    if (typeProjection == 1)
    {
        numLevels = numLevels - numDegree + 2;
    }

    // Define coefficients CDR equation (= Poisson)
    gsFunctionExpr<> coeff_diff, coeff_conv, coeff_reac;
    switch(sol_exact.domainDim())
    {
      case 1: coeff_diff = gsFunctionExpr<>("1",1);
              coeff_conv = gsFunctionExpr<>("0",1);
              coeff_reac = gsFunctionExpr<>("0",1); break;
      case 2: coeff_diff = gsFunctionExpr<>("1","0","0","1",2);
              coeff_conv = gsFunctionExpr<>("0","0",2);
              coeff_reac = gsFunctionExpr<>("0",2); break;
      case 3: coeff_diff = gsFunctionExpr<>("1","0","0","0","1","0","0","0","1",3);
              coeff_conv = gsFunctionExpr<>("0","0","0",3);
              coeff_reac = gsFunctionExpr<>("0",3);
    }
    if (numBenchmark == 1)
    {
        coeff_diff = gsFunctionExpr<>("1.2","-0.7","-0.4","0.9",2);
        coeff_conv = gsFunctionExpr<>("0.4","-0.2",2);
        coeff_reac = gsFunctionExpr<>("0.3",2);
    }

    // Apply p-Multigrid as stand-alone solver
    if (typeSolver == 1)
    {
        gsInfo << "p-multigrid is applied as stand-alone solver\n\n";
        My_MG.setup(f, sol_exact, x, numSmoothing, x, typeSolver, iterTot, typeCycle_p,typeCycle_h, numLevels, numCoarsening, numDegree, numRefine, numBenchmark, typeMultigrid, typeBCHandling, geo, typeLumping, hp, typeProjection, typeSmoother, typeCoarseOperator, coeff_diff, coeff_conv, coeff_reac);
        My_MG.solve(f, sol_exact, x, numSmoothing, x, typeSolver, iterTot, typeCycle_p,typeCycle_h, numLevels, numCoarsening, numDegree, numRefine, numBenchmark, typeMultigrid, typeBCHandling, geo, typeLumping, hp, typeProjection, typeSmoother, typeCoarseOperator);
        return 0;
    }

    // Assemble matrix fine level (p-mg as preconditioner)
    gsCDRAssembler<real_t> pa(mp, basisH, bcInfo, sol_exact, coeff_diff, coeff_conv, coeff_reac, (typeBCHandling == 1 ? dirichlet::elimination : dirichlet::nitsche), iFace::glue);
    pa.assemble();

    // Apply BiCGStab or CG
    x = gsMatrix<>::Random(pa.matrix().rows(),1);
    gsVector <> r = pa.rhs() - pa.matrix() * x;
    gsVector<> r0 = r;
    real_t maxIter = pa.matrix().rows();
    real_t tol = 1e-8;
    index_t i = 1;

    // Set up, determine initial L2 error
    gsField<> solKrylov;
    gsField<> sol = pa.constructSolution(x);
    real_t oldResNorm = r0.norm();

    // Perform BiCGStab
    if (typeSolver == 2)
    {
      gsInfo << "BiCGStab is applied as solver, p-multigrid as a preconditioner\n\n";
      // Define vectors needed in BiCGStab
      gsVector<> t = gsVector<>::Zero(pa.matrix().rows());
      gsVector<> s = gsVector<>::Zero(pa.matrix().rows());
      gsVector<> p = gsVector<>::Zero(pa.matrix().rows());
      gsVector<> v = gsVector<>::Zero(pa.matrix().rows());
      gsMatrix<> y = gsMatrix<>::Zero(pa.matrix().rows(),1);
      gsMatrix<> z = gsMatrix<>::Zero(pa.matrix().rows(),1);
      real_t alp = 1; real_t rho = 1; real_t w = 1;

      // Set residual norm and #restarts
      real_t r0_sqnorm = r0.dot(r0);
      index_t restarts = 0;

      // Construct P-Multigrid objects
      pMultigrid<real_t, gsSparseSolver<real_t>::LU,gsCDRAssembler<real_t> > My_MG1(mp, basisL, bcInfo);
      My_MG1.setup(f, sol_exact, y, numSmoothing, p, typeSolver, iterTot, typeCycle_p, typeCycle_h, numLevels, numCoarsening, numDegree, numRefine, numBenchmark, typeMultigrid, typeBCHandling, geo, typeLumping, hp, typeProjection, typeSmoother, typeCoarseOperator, coeff_diff, coeff_conv, coeff_reac);
      pMultigrid<real_t, gsSparseSolver<real_t>::LU,gsCDRAssembler<real_t> > My_MG2(mp, basisL, bcInfo);
      My_MG2.setup(f, sol_exact, z, numSmoothing, s, typeSolver, iterTot, typeCycle_p, typeCycle_h, numLevels, numCoarsening, numDegree, numRefine, numBenchmark, typeMultigrid,typeBCHandling, geo, typeLumping, hp, typeProjection, typeSmoother, typeCoarseOperator, coeff_diff, coeff_conv, coeff_reac);

      // Perform BiCGStab
      while(r.norm()/r0.norm() > tol && i < maxIter)
      {
        real_t rho_old = rho;
        rho = r0.dot(r);
        if (abs(rho) < 1e-32*r0_sqnorm)
        {
              gsInfo << "Residual too orthogonal, restart with new r0 \n";
              r = pa.rhs() - pa.matrix()*x;
              r0 = r;
              rho = r0_sqnorm = r.dot(r);
              if (restarts++ == 0)
              {
                Eigen::IncompleteLUT<real_t> ilu;
                ilu.setFillfactor(1);
                ilu.compute(pa.matrix());
                i=0;
              }
        }
        real_t beta = (rho/rho_old)*(alp/w);
        p = r + beta*(p - w*v);

        // Apply preconditioning by solving Ay = p
        y.setZero();
        My_MG1.solve(f, sol_exact, y, numSmoothing, p, typeSolver, iterTot, typeCycle_p, typeCycle_h, numLevels, numCoarsening, numDegree, numRefine, numBenchmark, typeMultigrid, typeBCHandling, geo, typeLumping, hp, typeProjection, typeSmoother, typeCoarseOperator);
        v = pa.matrix()*y;
        alp = rho/(r0.dot(v));
        s = r - alp*v;

        // Apply preconditioning by solving Az = s
        z.setZero();
        My_MG2.solve(f, sol_exact, z, numSmoothing, s, typeSolver, iterTot, typeCycle_p, typeCycle_h, numLevels, numCoarsening, numDegree, numRefine, numBenchmark, typeMultigrid,typeBCHandling, geo, typeLumping, hp, typeProjection, typeSmoother, typeCoarseOperator);
        t = pa.matrix()*z;
        if (t.dot(t) > 0)
          w = t.dot(s)/t.dot(t);
        else
          w = 0;
        x= x + alp*y + w*z;
        r = s - w*t;

        // Print information about residual and L2 error
        gsInfo << "BiCGStab iteration: " << i << "      |   Residual norm: "   << std::left << std::setw(15) << r.norm()
               << "           reduction:  1 / " << std::setprecision(3) << (oldResNorm/r.norm()) <<        std::setprecision  (6) << "\n";
        oldResNorm = r.norm();
        solKrylov = pa.constructSolution(x);

        ++i;
      }
  }
  else if (typeSolver == 3)
  {
      gsInfo << "CG is applied as solver, p-multigrid as a preconditioner\n\n";

      // Apply preconditioner
      pMultigrid<real_t, gsSparseSolver<real_t>::LU,gsCDRAssembler<real_t> > My_MG2(mp, basisL, bcInfo);
      gsMatrix<> z1 = gsMatrix<>::Zero(pa.matrix().rows(),1);
      My_MG2.setup(f, sol_exact, z1, numSmoothing, r0, typeSolver, iterTot, typeCycle_p, typeCycle_h, numLevels, numCoarsening, numDegree, numRefine, numBenchmark, typeMultigrid, typeBCHandling, geo, typeLumping, hp, typeProjection, typeSmoother, typeCoarseOperator, coeff_diff, coeff_conv, coeff_reac);
      My_MG2.solve(f, sol_exact, z1, numSmoothing, r0, typeSolver, iterTot, typeCycle_p, typeCycle_h, numLevels, numCoarsening, numDegree, numRefine, numBenchmark, typeMultigrid, typeBCHandling, geo, typeLumping, hp, typeProjection, typeSmoother, typeCoarseOperator);
      gsVector<> z = z1;
      gsVector<> p = z;
      real_t alpha, beta;
      pMultigrid<real_t, gsSparseSolver<real_t>::LU,gsCDRAssembler<real_t> > My_MG3(mp, basisL, bcInfo);
      gsMatrix<> z2 = gsMatrix<>::Zero(pa.matrix().rows(),1);
      My_MG3.setup(f, sol_exact, z2, numSmoothing, z2, typeSolver, iterTot, typeCycle_p, typeCycle_h, numLevels, numCoarsening, numDegree, numRefine, numBenchmark, typeMultigrid, typeBCHandling, geo,typeLumping, hp, typeProjection, typeSmoother, typeCoarseOperator, coeff_diff, coeff_conv, coeff_reac);

      while(r.norm()/r0.norm() >  tol && i < maxIter)
      {
        // Determine alpha
        alpha = r.transpose()*z;
        alpha = alpha/(p.transpose()*pa.matrix()*p);

        // Update solution and residual
        x = x + alpha*p;
        gsVector<> r_new;
        r_new = r - alpha*pa.matrix()*p;

        // Obtain new values
        gsMatrix<> z2 = gsMatrix<>::Zero(pa.matrix().rows(),1);
        My_MG3.solve(f, sol_exact, z2, numSmoothing,r_new, typeSolver, iterTot, typeCycle_p, typeCycle_h, numLevels, numCoarsening, numDegree, numRefine, numBenchmark, typeMultigrid, typeBCHandling, geo,typeLumping, hp, typeProjection, typeSmoother, typeCoarseOperator);
        gsVector<> z3 = z2;

        // Determine beta
        beta = z3.transpose()*r_new;
        beta = beta/(z.transpose()*r) ;
        p = z3 + beta*p;
        z = z3;
        r = r_new;

        // Print information about residual and L2 error
        gsInfo << "CG iteration: " << i << "       |  Residual norm: "   << std::left << std::setw(15) << r.norm()
               << "            reduction:  1 / " << std::setprecision(3) << (oldResNorm/r.norm()) <<  std::setprecision (6) << "\n";
        oldResNorm = r.norm();
        solKrylov = pa.constructSolution(x);
        ++i;
      }
  }
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

    gsBlockILUT( const gsSparseMatrix<>& op, const gsMultiPatch<>& mp, const gsMultiBasis<>& mb ) : m_op(op)
    {
        index_t numPatch = mp.nPatches();

        index_t shift0 = 0;

        // Use of partition functions
        std::vector<gsVector<index_t> > interior, boundary;
        std::vector<std::vector<gsVector<index_t> > > interface;
        std::vector<gsMatrix<index_t> >  global_interior, global_boundary;
        std::vector<std::vector<gsMatrix<index_t> > > global_interface;
        mb.partition(interior,boundary,interface,global_interior,global_boundary,global_interface);
        // Vector of vector of shift objects
        std::vector<index_t> shift(numPatch+1);
        for (index_t l=0; l< numPatch; l++)
        {
            shift[l] = global_interior[l].rows();
        }
        shift[numPatch] = 0;
        shift[numPatch] = m_op.rows() - accumulate(shift.begin(),shift.end(),0);

        // Vector of factorized operators
        std::vector< gsSparseMatrix<> > ILUT(numPatch+1);
        // Vector of factorized operators
        std::vector< Eigen::PermutationMatrix<Dynamic,Dynamic,index_t> > P(numPatch+1);
        // Vector of factorized operators
        std::vector< Eigen::PermutationMatrix<Dynamic,Dynamic,index_t> > Pinv(numPatch+1);
        for (index_t j = 0 ; j < numPatch ; j++)
        {
            const gsSparseMatrix<> block = m_op.block(shift0,shift0,shift[j],shift[j]);
            Eigen::IncompleteLUT<real_t> ilu;
            ilu.setFillfactor(1);
            ilu.compute(block);
            ILUT[j] = ilu.factors();
            P[j] = ilu.fillReducingPermutation();
            Pinv[j] = ilu.inversePermutation();
            shift0 = shift0 + shift[j];
        }

        shift0 = 0;
        // Obtain the blocks of the matrix
        std::vector< gsSparseMatrix<> > ddB(numPatch+1);
        std::vector< gsSparseMatrix<> > ddC(numPatch+1);

        for (index_t j = 0 ; j < numPatch+1 ; j++)
        {
            ddB[j] = m_op.block(m_op.rows()-shift[numPatch],shift0,shift[numPatch],shift[j]);
            ddC[j] = m_op.block(shift0,m_op.cols()-shift[numPatch],shift[j],shift[numPatch]);
            shift0 = shift0 + shift[j];
        }
        shift0 = 0;

        // Define the A_aprox matrix
        m_A_aprox = gsSparseMatrix<>(m_op.rows(),m_op.cols());

        // Retrieve a block of each patch
        for (index_t k = 0; k < numPatch; k++)
        {
            m_A_aprox.block(shift0,shift0,shift[k],shift[k]) = ILUT[k];
            shift0 = shift0 + shift[k];
        }
        shift0 = 0;
        std::vector< gsMatrix<> > ddBtilde(numPatch);
        std::vector< gsMatrix<> > ddCtilde(numPatch);

        for (index_t j=0 ; j < numPatch ; j++)
        {
            ddBtilde[j] = gsSparseMatrix<>(shift[j],shift[numPatch]);
            ddCtilde[j] = gsSparseMatrix<>(shift[j],shift[numPatch]);
            for (index_t k=0 ; k < shift[numPatch]; k++)
            {
                gsMatrix<> Brhs = ddC[j].col(k);
                gsMatrix<> Crhs = ddC[j].col(k);
                ddBtilde[j].col(k) = ILUT[j].triangularView<Eigen::Upper>().transpose().solve(Brhs);
                ddCtilde[j].col(k) = ILUT[j].triangularView<Eigen::UnitLower>().solve(Crhs);
            }
        }
        // Define matrix S
        gsSparseMatrix<> S = ddC[numPatch];
        for (index_t l = 0 ; l < numPatch ; l++)
        {
            S -= ddBtilde[l].transpose()*ddCtilde[l];
        }

        // Fill matrix A_aprox
        for (index_t m = 0 ; m < numPatch ; m++)
        {
          m_A_aprox.block(shift0,m_A_aprox.rows() - shift[numPatch],shift[m],shift[numPatch]) = ddCtilde[m];
          m_A_aprox.block(m_A_aprox.rows() - shift[numPatch],shift0,shift[numPatch],shift[m]) = ddBtilde[m].transpose();
          shift0 += shift[m];
        }
        shift0 = 0;

        // Preform ILUT on the S-matrix!
        Eigen::IncompleteLUT<real_t> ilu;
        ilu.setFillfactor(1);
        gsSparseMatrix<> II = S;
        ilu.compute(II);
        m_A_aprox.block(
            m_A_aprox.rows() - shift[numPatch],
            m_A_aprox.rows() - shift[numPatch],
            shift[numPatch],
            shift[numPatch]
        ) = ilu.factors();
    }

    void step(const gsMatrix<>& rhs, gsMatrix<>& x) const
    {
        gsMatrix<> d = rhs-m_op*x;
        gsMatrix<> e = m_A_aprox.template triangularView<Eigen::UnitLower>().solve(d);
        x += m_A_aprox.template triangularView<Eigen::Upper>().solve(e);
    }

    void stepT(const gsMatrix<>& rhs, gsMatrix<>& x) const
    {
        gsMatrix<> d = rhs-m_op*x;
        gsMatrix<> e = m_A_aprox.template triangularView<Eigen::UnitLower>().solve(d);
        x += m_A_aprox.template triangularView<Eigen::Upper>().solve(e);
    }

    index_t rows() const { return m_op.rows(); }
    index_t cols() const { return m_op.cols(); }
    gsLinearOperator<real_t>::Ptr underlyingOp() const { return makeMatrixOp(m_op); }

private:
    /// Underlying matrix
    gsSparseMatrix<> m_op;

    /// Block operator object
    gsMatrix<> m_A_aprox;
};

gsPreconditionerOp<>::Ptr setupBlockILUT(
    const gsSparseMatrix<>& op,
    const gsMultiPatch<>& mp,
    const gsMultiBasis<>& mb
)
{
    return gsBlockILUT::uPtr(new gsBlockILUT(op, mp, mb));
}
