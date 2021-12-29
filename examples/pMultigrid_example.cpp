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

using namespace std;
using namespace gismo;

gsPreconditionerOp<>::Ptr setupSubspaceCorrectedMassSmoother(
    const gsSparseMatrix<>&,
    const gsMultiBasis<>&,
    const gsBoundaryConditions<>&,
    const gsOptionList&,
    const int&
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
        vector<memory::shared_ptr<gsMultiBasis<T> > > m_basis,
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
        vector<gsSparseMatrix<T>>& m_prolongation_P,
        vector<gsSparseMatrix<T>>& m_restriction_P,
        vector<gsMatrix<T>>& m_prolongation_M,
        vector<gsMatrix<T>>& m_restriction_M,
        vector<gsSparseMatrix<T>>& m_prolongation_H,
        vector<gsSparseMatrix<T>>& m_restriction_H,
        const gsMatrix<>& hp
        )
    {
        if( numLevels == 1)
        {
            solvecoarse(rhs, x, numLevels);
            return;
        }

        if(hp(max(numLevels-2,0),0) == 0 )
        {
            gsMatrix<T> fineRes, coarseRes, fineCorr, coarseCorr, postRes;
            presmoothing(rhs, x, numLevels, numSmoothing, fineRes, numRefine, typeSmoother,hp);
            restriction(fineRes, coarseRes, numLevels, numCoarsening, m_basis, typeLumping,
                typeBCHandling, bcInfo, mp, geo, typeProjection, m_prolongation_P, m_restriction_P,
                m_prolongation_M, m_restriction_M, m_prolongation_H, m_restriction_H, hp);
            coarseCorr.setZero(coarseRes.rows(),1);
            for( int j = 0 ; j < (typeCycle_p == 2 ? 2 : 1) ; j++)
            {
                  solve(coarseRes, m_basis, coarseCorr, numLevels-1, numCoarsening, numRefine, numSmoothing,
                      numCoarseCycles, typeCycle_p, typeCycle_h, typeSolver, typeBCHandling, bcInfo, mp, geo,
                      typeLumping, typeProjection, typeSmoother, m_prolongation_P, m_restriction_P, m_prolongation_M,
                      m_restriction_M,m_prolongation_H, m_restriction_H, hp);
            }
            prolongation(coarseCorr, fineCorr, numLevels, numCoarsening, m_basis, typeLumping,
                typeBCHandling, bcInfo, mp, geo, typeProjection, m_prolongation_P, m_restriction_P,
                m_prolongation_M, m_restriction_M, m_prolongation_H, m_restriction_H, hp);
            postsmoothing(rhs,x, numLevels, numSmoothing, fineCorr, postRes, typeSolver, numRefine, typeSmoother,hp);
        }

        if(hp(max(numLevels-2,0),0) == 1 )
        {
            gsMatrix<T> fineRes, coarseRes, fineCorr, coarseCorr, postRes;
            presmoothing(rhs, x, numLevels, numSmoothing, fineRes, numRefine, typeSmoother,hp);
            restriction(fineRes, coarseRes, numLevels, numCoarsening, m_basis, typeLumping,
                typeBCHandling, bcInfo, mp, geo, typeProjection, m_prolongation_P, m_restriction_P,
                m_prolongation_M, m_restriction_M, m_prolongation_H, m_restriction_H, hp);
            coarseCorr.setZero(coarseRes.rows(),1);
            for( int i = 0 ; i < (typeCycle_h == 2 ? 2 : 1) ; i++)
            {
                solve(coarseRes, m_basis, coarseCorr, numLevels-1, numCoarsening, numRefine, numSmoothing,
                    numCoarseCycles, typeCycle_p, typeCycle_h, typeSolver, typeBCHandling, bcInfo, mp, geo,
                    typeLumping, typeProjection, typeSmoother, m_prolongation_P, m_restriction_P,
                    m_prolongation_M, m_restriction_M, m_prolongation_H, m_restriction_H, hp);
            }
            prolongation(coarseCorr, fineCorr, numLevels, numCoarsening, m_basis, typeLumping,
                typeBCHandling, bcInfo, mp, geo, typeProjection, m_prolongation_P, m_restriction_P,
                m_prolongation_M, m_restriction_M, m_prolongation_H, m_restriction_H, hp);
            postsmoothing(rhs,x, numLevels, numSmoothing, fineCorr, postRes, typeSolver, numRefine, typeSmoother,hp);
        }
    }

    /// @brief Setup p-multigrid to given linear system
    virtual void setup(const gsMatrix<T> & rhs,
        vector<memory::shared_ptr<gsMultiBasis<T> > > m_basis ,
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
        vector<gsSparseMatrix<T>>& m_prolongation_P,
        vector<gsSparseMatrix<T>>& m_restriction_P,
        vector<gsMatrix<T>>& m_prolongation_M,
        vector<gsMatrix<T>>& m_restriction_M,
        vector<gsSparseMatrix<T>>& m_prolongation_H,
        vector<gsSparseMatrix<T>>& m_restriction_H,
        const gsMatrix<>& hp)
    {}

    /// @brief Apply fixed number of smoothing steps (pure virtual method)
    virtual void presmoothing(const gsMatrix<T>& rhs, gsMatrix<T>& x, const int& numLevels,
        const int& numSmoothing, gsMatrix<T> & fineRes , const int& numRefine, const int& typeSmoother,
        const gsMatrix<>& hp) = 0;

    /// @brief Apply fixed number of smoothing steps (pure virtual method)
    virtual void postsmoothing(const gsMatrix<T>& rhs, gsMatrix<T>& x, const int& numLevels,
        const int& numSmoothing, gsMatrix<T> & fineCorr, gsMatrix<T> & postRes, const int& typeSolver,
        const int& numRefine, const int& typeSmoother, const gsMatrix<>& hp) = 0;

    /// @brief Apply coarse solver (pure virtual method)
    virtual void solvecoarse(const gsMatrix<T>& rhs, gsMatrix<T>& x, const int& numLevels) = 0;

    /// @brief Prolongate coarse space function to fine space
    virtual gsSparseMatrix<T> prolongation_P(const int& numLevels, vector<memory::shared_ptr<gsMultiBasis<T> > > m_basis,
        const int& typeLumping, const int& typeBCHandling, gsGeometry<>::Ptr geo, const int& typeProjection) = 0;

    /// @brief Prolongate coarse space function to fine space
    virtual gsSparseMatrix<T> restriction_P(const int& numLevels, vector<memory::shared_ptr<gsMultiBasis<T> > > m_basis,
        const int& typeLumping, const int& typeBCHandling, gsGeometry<>::Ptr geo, const int& typeProjection) = 0;

    /// @brief Prolongate coarse space function to fine space
    virtual gsMatrix<T> prolongation_M(const int& numLevels, vector<memory::shared_ptr<gsMultiBasis<T> > > m_basis,
        const int& typeLumping, const int& typeBCHandling, gsGeometry<>::Ptr geo, const int& typeProjection) = 0;

    /// @brief Prolongate coarse space function to fine space
    virtual gsMatrix<T> restriction_M(const int& numLevels, vector<memory::shared_ptr<gsMultiBasis<T> > > m_basis,
        const int& typeLumping, const int& typeBCHandling, gsGeometry<>::Ptr geo, const int& typeProjection) = 0;

    /// @brief Prolongate coarse space function to fine space
    virtual void prolongation(const gsMatrix<T>& Xcoarse, gsMatrix<T>& Xfine, const int& numLevels, const int& numCoarsening,
        vector<memory::shared_ptr<gsMultiBasis<T> > > m_basis, const int& typeLumping, const int& typeBCHandling,
        gsBoundaryConditions<T> bcInfo, gsMultiPatch<> mp, gsGeometry<>::Ptr geo, const int& typeProjection,
        vector<gsSparseMatrix<T>>& m_prolongation_P, vector<gsSparseMatrix<T>>& m_restriction_P,
        vector<gsMatrix<T>>& m_prolongation_M, vector<gsMatrix<T>>& m_restriction_M,
        vector<gsSparseMatrix<T>>& m_prolongation_H, vector<gsSparseMatrix<T>>& m_restriction_H,
        const gsMatrix<>& hp)
    {
        if(hp(numLevels-2,0) == 1)
        {
            Xfine = m_prolongation_H[numLevels-2]*Xcoarse;
        }
        else
        {
            if(typeLumping == 1)
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
                if(typeBCHandling == 1)
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
        vector<memory::shared_ptr<gsMultiBasis<T> > > m_basis, const int& typeLumping, const int& typeBCHandling,
        const gsBoundaryConditions<T> & bcInfo, gsMultiPatch<> mp, gsGeometry<>::Ptr geo, const int& typeProjection,
        vector<gsSparseMatrix<T>>& m_prolongation_P,vector<gsSparseMatrix<T>>& m_restriction_P,
        vector<gsMatrix<T>>& m_prolongation_M, vector<gsMatrix<T>>& m_restriction_M,
        vector<gsSparseMatrix<T>>& m_prolongation_H, vector<gsSparseMatrix<T>>& m_restriction_H, const gsMatrix<>& hp)
    {
        if(hp(numLevels-2,0) == 1)
        {
            Xcoarse = m_restriction_H[numLevels-2]*Xfine;
        }
        else
        {
            if(typeLumping == 1)
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
                if(typeBCHandling == 1)
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
  vector<memory::shared_ptr<gsMultiBasis<T> > > m_basis;

  /// Vector of prolongation operators
  vector< gsSparseMatrix<T> > m_prolongation_P;

  /// Vector of restriction operators
  vector< gsSparseMatrix<T> > m_restriction_P;

  /// Vector of prolongation operators
  vector< gsMatrix<T> > m_prolongation_M;

  /// Vector of restriction operators
  vector< gsMatrix<T> > m_restriction_M;

  /// Vector of prolongation operators
  vector< gsSparseMatrix<T> > m_prolongation_H;

  /// Vector of restriction operators
  vector< gsSparseMatrix<T> > m_restriction_H;

  /// Vector of factorized operators
  vector< vector< gsSparseMatrix<T> > > m_ILUT;

  /// Vector of factorized operators
  vector< vector < Eigen::PermutationMatrix<Dynamic,Dynamic,index_t> > > m_P;

  /// Vector of factorized operators
  vector < vector < Eigen::PermutationMatrix<Dynamic,Dynamic,index_t> > > m_Pinv;

  /// Vector of SCM smoother object
  vector< gsPreconditionerOp<>::Ptr > m_SCMS;

  /// Vector of operator objects
  vector< gsSparseMatrix<T> > m_operator;

  /// Vector of vector of block operator objects
  vector < vector< gsSparseMatrix<T> > > m_block_operator;

  /// Vector of vector of block operator objects
  vector < vector  < gsSparseMatrix<T> > > m_ddB;

  /// Vector of vector of block operator objects
  vector < vector  < gsSparseMatrix<T> > > m_ddC;

  /// Vector of vector of block operator objects
  vector < vector <  gsMatrix<T>  > > m_ddBtilde;

  /// Vector of vector of block operator objects
  vector < vector <  gsMatrix<T>  > > m_ddCtilde;

  /// Vector of vector of block operator objects
  vector <  gsMatrix<T> > m_A_aprox;

  /// Vector of vector of block operator objects
  vector <  gsSparseMatrix<T> > m_S;

  /// Vector of vector of shift objects
  vector < std::vector< int > > m_shift;

  /// Vector of assembler objects
  vector<Assembler> m_assembler;

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
    for (typename vector<memory::shared_ptr<gsMultiBasis<T> > >::iterator it = m_basis.begin();
       it != m_basis.end(); ++it)
    {
      m_assembler.push_back(Assembler(*m_mp_ptr,*(*it).get(),*m_bcInfo_ptr,rhs, coeff_diff , coeff_conv , coeff_reac ,(typeBCHandling == 1 ? dirichlet::elimination : dirichlet::nitsche), iFace::glue));
    }

    // Resize vectors of operators
    m_operator.resize(numLevels); m_prolongation_P.resize(numLevels-1); m_prolongation_M.resize(numLevels-1); m_prolongation_H.resize(numLevels-1); m_restriction_P.resize(numLevels-1); m_restriction_M.resize(numLevels-1); m_restriction_H.resize(numLevels-1);

    // Assemble operators at finest level
    gsStopwatch clock;
    gsInfo << "|| Multigrid hierarchy ||" << endl;
    for (int i = 0; i < numLevels; i++)
    {
      gsInfo << "Level " << i+1 << " " ;
      if(typeCoarseOperator == 1)
      {
        m_assembler[i].assemble();
        m_operator[i] = m_assembler[i].matrix();
        gsInfo << "Degree: " << m_basis[i]->degree()  << ", Ndof: " << m_basis[i]->totalSize() << endl;
      }
      else
      {
        if(hp(min(i,hp.rows()-1),0) == 0 || i == numLevels-1)
        {
          m_assembler[i].assemble();
          m_operator[i] = m_assembler[i].matrix();
          gsInfo << "\nDegree of the basis: " << m_basis[i]->degree() << endl;
          gsInfo << "Size of the basis functions: " << m_basis[i]->totalSize() << endl;
        }
      }
    }
    real_t Time_Assembly = clock.stop();

    // Determine prolongation/restriction operators in p
    clock.restart();
    for (int i = 1; i < numLevels; i++)
    {
      if(hp(i-1,0) == 0)
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
    for(int i = 1; i < numLevels; i++)
    {
      if(hp(i-1,0) == 1)
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
    if(typeCoarseOperator == 2)
    {
      for (int i = numLevels-1; i > -1; i--)
      {
        if(hp(hp.rows()-1,0) == 0)
        {
          if(hp(min(i,hp.rows()-1),0) == 1)
          {
            m_operator[i] = m_restriction_H[i]*m_operator[i+1]*m_prolongation_H[i];
          }
        }
        else
        {
          if(hp(min(i,hp.rows()-1),0) == 1 && i > 0)
          {
            m_operator[i-1] = m_restriction_H[i-1]*m_operator[i]*m_prolongation_H[i-1];
          }
        }
      }
    }
    real_t Time_Assembly_Galerkin = clock.stop();


    // Setting up the subspace corrected mass smoother
    clock.restart();
    if(typeSmoother == 3)
    {
      // Generate sequence of SCM smoothers
      m_SCMS.resize(numLevels);
      gsOptionList opt;
      opt.addReal("Scaling","",0.12);
      for(int i = 0 ; i < numLevels ; i++)
      {
        m_SCMS[i] = setupSubspaceCorrectedMassSmoother(m_operator[i], *m_basis[i], *m_bcInfo_ptr, opt, typeBCHandling);
      }
    }
    real_t Time_SCMS = clock.stop();

    // Determine ILUT factorizations at each level
    clock.restart();
    int numPatch = m_mp_ptr->nPatches();

    if(typeSmoother == 1)
    {
      // Generate factorizations (ILUT)
      m_ILUT.resize(numLevels);
      m_P.resize(numLevels);
      m_Pinv.resize(numLevels);
      for(int i = 0; i < numLevels; i++)
      {
        m_ILUT[i].resize(1);
        m_P[i].resize(1);
        m_Pinv[i].resize(1);
        if(typeProjection == 2)
        {
            Eigen::IncompleteLUT<real_t> ilu;
            ilu.setFillfactor(1);
            ilu.compute(m_operator[i]);
            m_ILUT[i][0] = ilu.factors();
            m_P[i][0] = ilu.fillReducingPermutation();
            m_Pinv[i][0] = ilu.inversePermutation();
        }
        else
        {
            if(i == numLevels-1) // Only at finest level
            {
              Eigen::IncompleteLUT<real_t> ilu;
              ilu.setFillfactor(1);
              ilu.compute(m_operator[i]);
              m_ILUT[i][0] = ilu.factors();
              m_P[i][0] = ilu.fillReducingPermutation();
              m_Pinv[i][0] = ilu.inversePermutation();
            }
        }
      }
    }
    real_t Time_ILUT_Factorization = clock.stop();
    clock.restart();
    if(typeSmoother == 5)
    {
      int shift0 = 0;
      m_ddB.resize(numLevels);
      m_ddC.resize(numLevels);
      m_ddBtilde.resize(numLevels);
      m_ddCtilde.resize(numLevels);

      m_ILUT.resize(numLevels);
      m_P.resize(numLevels);
      m_Pinv.resize(numLevels);
      m_shift.resize(numLevels);
      m_S.resize(numLevels);

      for(int i = 0 ; i < numLevels ; i++)
      {
        m_shift[i].resize(numPatch+1);
        m_ILUT[i].resize(numPatch+1);
        m_P[i].resize(numPatch+1);
        m_Pinv[i].resize(numPatch+1);

        // Use of partition functions
        std::vector<gsVector<index_t> > interior, boundary;
        std::vector<std::vector<gsVector<index_t> > > interface;
        std::vector<gsMatrix<index_t> >  global_interior, global_boundary;
        std::vector<std::vector<gsMatrix<index_t> > > global_interface;
        m_basis[i]->partition(interior,boundary,interface,global_interior,global_boundary,global_interface);
        for(int l=0; l< numPatch; l++)
        {
          m_shift[i][l] = global_interior[l].rows();
        }
        m_shift[i][numPatch] = 0;
        m_shift[i][numPatch] = m_operator[i].rows() - accumulate(m_shift[i].begin(),m_shift[i].end(),0);

        // Put shift on zero
        shift0 = 0;
        for(int j = 0 ; j < numPatch ; j++)
        {
          const gsSparseMatrix<> block = m_operator[i].block(shift0,shift0,m_shift[i][j],m_shift[i][j]);
          Eigen::IncompleteLUT<real_t> ilu;
          ilu.setFillfactor(1);
          ilu.compute(block);
          m_ILUT[i][j] = ilu.factors();

          m_P[i][j] = ilu.fillReducingPermutation();
          m_Pinv[i][j] = ilu.inversePermutation();
         shift0 = shift0 + m_shift[i][j];

        }

        shift0 = 0;
        // Obtain the blocks of the matrix
        m_ddB[i].resize(numPatch+1);
        m_ddC[i].resize(numPatch+1);

        for(int j = 0 ; j < numPatch+1 ; j++)
        {
            m_ddB[i][j] = m_operator[i].block(m_operator[i].rows()-m_shift[i][numPatch],shift0,m_shift[i][numPatch],m_shift[i][j]);
            m_ddC[i][j] = m_operator[i].block(shift0,m_operator[i].cols()-m_shift[i][numPatch],m_shift[i][j],m_shift[i][numPatch]);
            shift0 = shift0 + m_shift[i][j];
        }
        shift0 = 0;
      }

      m_A_aprox.resize(numLevels);
      for(int i = 0 ; i < numLevels ; i++)
      {
         // Define the A_aprox matrix
        m_A_aprox[i] = gsSparseMatrix<>(m_operator[i].rows(),m_operator[i].cols());

        // Retrieve a block of each patch
        for(int k=0; k< numPatch; k++)
        {
          m_A_aprox[i].block(shift0,shift0,m_shift[i][k],m_shift[i][k]) = m_ILUT[i][k];
          shift0 = shift0 + m_shift[i][k];
        }
        shift0 = 0;
        m_ddBtilde[i].resize(numPatch);
        m_ddCtilde[i].resize(numPatch);

        for(int j=0 ; j < numPatch ; j ++)
        {
          m_ddBtilde[i][j] = gsSparseMatrix<>(m_shift[i][j],m_shift[i][numPatch]);
          m_ddCtilde[i][j] = gsSparseMatrix<>(m_shift[i][j],m_shift[i][numPatch]);
          for(int k=0 ; k < m_shift[i][numPatch]; k++)
          {
            gsMatrix<> Brhs = m_ddC[i][j].col(k);
            gsMatrix<> Crhs = m_ddC[i][j].col(k);
            m_ddBtilde[i][j].col(k) = m_ILUT[i][j].template triangularView<Eigen::Upper>().transpose().solve(Brhs);
            m_ddCtilde[i][j].col(k) = m_ILUT[i][j].template triangularView<Eigen::UnitLower>().solve(Crhs);
          }
        }

        // Define matrix S
        m_S[i] = m_ddC[i][numPatch];
        for(int l = 0 ; l < numPatch ; l++)
        {
          m_S[i] = m_S[i] - m_ddBtilde[i][l].transpose()*m_ddCtilde[i][l];
        }

        // Fill matrix A_aprox
        for(int m = 0 ; m < numPatch ; m++)
        {
          m_A_aprox[i].block(shift0,m_A_aprox[i].rows() - m_shift[i][numPatch],m_shift[i][m],m_shift[i][numPatch]) = m_ddCtilde[i][m];
          m_A_aprox[i].block(m_A_aprox[i].rows() - m_shift[i][numPatch],shift0,m_shift[i][numPatch],m_shift[i][m]) = m_ddBtilde[i][m].transpose();
          shift0 = shift0 + m_shift[i][m];
        }
        shift0 = 0;

        // Preform ILUT on the S-matrix!
        Eigen::IncompleteLUT<real_t> ilu;
        ilu.setFillfactor(1);
        gsSparseMatrix<> II = m_S[i];
        ilu.compute(II);
        m_A_aprox[i].block(m_A_aprox[i].rows() - m_shift[i][numPatch],m_A_aprox[i].rows() - m_shift[i][numPatch],m_shift[i][numPatch],m_shift[i][numPatch]) = ilu.factors();
      }
    }

    real_t Time_Block_ILUT_Factorization = clock.stop();
    gsInfo << "\n|| Setup Timings || " << endl;
    gsInfo << "Total Assembly time: " << Time_Assembly << endl;
    gsInfo << "Total Assembly time (Galerkin): " << Time_Assembly_Galerkin << endl;
    gsInfo << "Total ILUT factorization time: " << Time_ILUT_Factorization << endl;
    gsInfo << "Total block ILUT factorization time: " << Time_Block_ILUT_Factorization << endl;
    gsInfo << "Total SCMS time: " << Time_SCMS << endl;
    gsInfo << "Total setup time: " << Time_Assembly_Galerkin + Time_Assembly + Time_Transfer + Time_ILUT_Factorization + Time_SCMS << endl;
  }

  ///  @brief Apply p-multigrid solver to given right-hand side on level l
  void solve(const gsFunctionExpr<T> & rhs, const gsFunctionExpr<T> & sol_exact, gsMatrix<T>& x, const int& numSmoothing, gsMatrix<T> f,const int& typeSolver, int& iterTot, int& typeCycle_p, int& typeCycle_h, int numLevels, const int& numCoarsening, const int& numDegree, const int& numRefine, const int& numBenchmark, const int& typeMultigrid, const int& typeBCHandling, gsGeometry<>::Ptr geo, const int& typeLumping, const gsMatrix<>& hp, const int& typeProjection,const int& typeSmoother, const int& typeCoarseOperator)
  {
      gsStopwatch clock;

      if(typeSolver == 1)
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
          Base::solve(b, m_basis,  x, numLevels, numCoarsening, numRefine, numSmoothing, numCoarseCycles, typeCycle_p, typeCycle_h
            , typeSolver, typeBCHandling, *m_bcInfo_ptr, *m_mp_ptr, geo, typeLumping, typeProjection, typeSmoother, m_prolongation_P, m_restriction_P,m_prolongation_M, m_restriction_M,m_prolongation_H, m_restriction_H, hp);
          numCoarseCycles = 0;
          r = (m_operator[numLevels-1]*x - b).norm();
          if( r_old < r)
          {
            gsInfo << "Residual increased during solving!!! " << endl;
          }
          r_old = r;
          //gsInfo << "Residual after cycle " << iter << " equals: " << r << endl;
          iter++;
          iterTot++;
      }
      real_t Time_Solve = clock.stop();
      gsInfo << "\n|| Solver information || " << endl;
      gsInfo << "Solver converged in " << Time_Solve << " seconds!" << endl;
      gsInfo << "Solver converged in " << iter-1 << " iterations!" << endl;

      if(typeSolver == 1)
      {
          // Determine residual and L2 errpr
          gsInfo << "Residual after solving: "  << (b-m_operator[numLevels-1]*x).norm() << endl;
      }
  }

private:

  /// @brief Apply coarse solver
  virtual void solvecoarse(const gsMatrix<T>& rhs, gsMatrix<T>& x, const int& numLevels)
  {
      //gsInfo << "Coarse solver is applied! " << endl;

      // Direct solver (LU factorization)
      CoarseSolver solver;
      solver.analyzePattern(m_operator[0]);
      solver.factorize(m_operator[0]);
      x = solver.solve(rhs);
  }

  /// @brief Construct prolongation operator at level numLevels
  virtual gsMatrix<T> prolongation_M(const int& numLevels, vector<memory::shared_ptr<gsMultiBasis<T> > > m_basis, const int& typeLumping, const int& typeBCHandling, gsGeometry<>::Ptr geo, const int& typeProjection)
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
      if(typeBCHandling == 1)
      {
          w_n.setup(*m_bcInfo_ptr, dirichlet::interpolation, 0);
      }
      ex2.setIntegrationElements(basisH);
      ex2.initSystem();
      ex2.assemble(w_n * meas(G2) );
      return ex2.rhs();
  }

  /// @brief Construct prolongation operator at level numLevels
  virtual gsSparseMatrix<T> prolongation_P(const int& numLevels, vector<memory::shared_ptr<gsMultiBasis<T> > > m_basis, const int& typeLumping, const int& typeBCHandling, gsGeometry<>::Ptr geo, const int& typeProjection)
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
      if(typeBCHandling == 1)
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
  virtual gsMatrix<T> restriction_M(const int& numLevels, vector<memory::shared_ptr<gsMultiBasis<T> > > m_basis, const int& typeLumping, const int& typeBCHandling, gsGeometry<>::Ptr geo, const int& typeProjection)
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
      if(typeBCHandling == 1)
      {
          w_n.setup(*m_bcInfo_ptr, dirichlet::interpolation, 0);
      }
      ex2.setIntegrationElements(basisL);
      ex2.initSystem();
      ex2.assemble(w_n * meas(G2) );
      return ex2.rhs();
  }

  /// @brief Construct restriction operator at level numLevels
  virtual gsSparseMatrix<T> restriction_P(const int& numLevels, vector<memory::shared_ptr<gsMultiBasis<T> > > m_basis, const int& typeLumping, const int& typeBCHandling, gsGeometry<>::Ptr geo, const int& typeProjection)
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
      if( typeBCHandling == 1)
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

  /// @brief Apply fixed number of presmoothing steps
  virtual void presmoothing(const gsMatrix<T>& rhs,  gsMatrix<T>& x, const int& numLevels, const int& numSmoothing, gsMatrix<T> & fineRes, const int& numRefine, const int& typeSmoother, const gsMatrix<>& hp)
  {
    for(int i = 0 ; i < numSmoothing ; i++)
    {
      if(typeSmoother == 1)
      {
        if(hp(numLevels-2,0) == 1 && hp(hp.rows()-1,0) == 0)
        {
          internal::gaussSeidelSweep(m_operator[numLevels-1],x,rhs);
        }
        else
        {
          gsMatrix<> e;
          gsMatrix<> d = rhs-m_operator[numLevels-1]*x;
          e = m_Pinv[numLevels-1][0]*d;
          e = m_ILUT[numLevels-1][0].template triangularView<Eigen::UnitLower>().solve(e);
          e = m_ILUT[numLevels-1][0].template triangularView<Eigen::Upper>().solve(e);
          e = m_P[numLevels-1][0]*e;
          x = x + e;
        }
      }
      if(typeSmoother == 2)
      {
        internal::gaussSeidelSweep(m_operator[numLevels-1],x,rhs);
      }
      if(typeSmoother == 3)
      {
          m_SCMS[numLevels-1]->step(rhs,x);
      }
      if(typeSmoother == 5)
      {
        if(hp(numLevels-2,0) == 1 && hp(hp.rows()-1,0) == 0)
        {
          internal::gaussSeidelSweep(m_operator[numLevels-1],x,rhs);
        }
        else
        {
          gsMatrix<> e;
          gsMatrix<> d = rhs-m_operator[numLevels-1]*x;
          e = m_A_aprox[numLevels-1].template triangularView<Eigen::UnitLower>().solve(d);
          e = m_A_aprox[numLevels-1].template triangularView<Eigen::Upper>().solve(e);
          x = x + e;
        }
      }
    }
    fineRes = m_operator[numLevels-1]*x - rhs;
  }

  /// @brief Apply fixed number of postsmoothing steps
  virtual void postsmoothing(const gsMatrix<T>& rhs,  gsMatrix<T>& x, const int& numLevels, const int& numSmoothing, gsMatrix<T> & fineCorr, gsMatrix<T> & postRes, const int& typeSolver, const int& numRefine, const int& typeSmoother, const gsMatrix<>& hp)
  {
    real_t alpha = 1;
    x = x - alpha*fineCorr;
    for(int i = 0 ; i < numSmoothing ; i++)
    {
      if(typeSmoother == 1)
      {
        if(hp(numLevels-2,0) == 1 && hp(hp.rows()-1,0) == 0)
        {
         ( typeSolver == 3 ? internal::reverseGaussSeidelSweep(m_operator[numLevels-1],x,rhs) : internal::gaussSeidelSweep(m_operator[numLevels-1],x,rhs));
        }
        else
        {
          gsMatrix<> e;
          gsMatrix<> d = rhs-m_operator[numLevels-1]*x;
          e = m_Pinv[numLevels-1][0]*d;
          e = m_ILUT[numLevels-1][0].template triangularView<Eigen::UnitLower>().solve(e);
          e = m_ILUT[numLevels-1][0].template triangularView<Eigen::Upper>().solve(e);
          e = m_P[numLevels-1][0]*e;
          x = x + e;
        }
      }
      if(typeSmoother == 2)
      {
        ( typeSolver == 3 ? internal::reverseGaussSeidelSweep(m_operator[numLevels-1],x,rhs) : internal::gaussSeidelSweep(m_operator[numLevels-1],x,rhs));
      }
      if(typeSmoother == 3)
      {
         m_SCMS[numLevels-1]->step(rhs,x);
      }
      if(typeSmoother == 5)
      {
        if(hp(numLevels-2,0) == 1 && hp(hp.rows()-1,0) == 0)
        {
         ( typeSolver == 3 ? internal::reverseGaussSeidelSweep(m_operator[numLevels-1],x,rhs) : internal::gaussSeidelSweep(m_operator[numLevels-1],x,rhs));
        }
        else
        {
          gsMatrix<> e;
          gsMatrix<> d = rhs-m_operator[numLevels-1]*x;
          e = m_A_aprox[numLevels-1].template triangularView<Eigen::UnitLower>().solve(d);
          e = m_A_aprox[numLevels-1].template triangularView<Eigen::Upper>().solve(e);
          x = x + e;
        }
      }
      postRes = rhs - m_operator[numLevels-1]*x;
    }
  }
};

/** @brief The p-multigrid class implements a generic p-multigrid solver
 *  that can be customized by passing assembler and coarse
 *  solver as template arguments.
 *
 *  @note: This implementation assumes that all required prolongation/
 *  restriction operators are generated externally and provided as
 *  constant references through the constructor. Therefore, no assembler
 *  is passed as template parameter.
 */
template<class T, class CoarseSolver>
struct pMultigrid<T, CoarseSolver, void> : public pMultigridBase<T>
{
  // Default constructor
  pMultigrid()
  {
      gsInfo << "The specific case";
  }
};

int main(int argc, char* argv[])
{
    int numDegree = 2;
    int numRefine = 6;
    int numSmoothing = 1;
    int numLevels = 2;
    int numBenchmark = 3;
    int numPatches = 1;
    int typeSolver = 1;
    int typeCycle_p = 1;
    int typeCycle_h = 2;
    int typeMultigrid = 1;
    int typeBCHandling = 2;
    int typeLumping = 1;
    int typeProjection = 2;
    int typeSmoother = 1;
    int typeCoarseOperator = 1;
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
    cmd.addInt("S", "Smoother", "Type of smoother: (1) ILUT (2) Gauss-Seidel (3) SCMS (4) Block ILUT (5) Block Gauss-Seidel", typeSmoother);
    cmd.addInt("G", "CoarseOperator", "Type of coarse operator in h-multigrid: (1) Rediscretization (2) Galerkin Projection", typeCoarseOperator);
    cmd.addString("z", "Coarsening", "Expression that defines coarsening strategy", typeCoarsening);

    // Read parameters from command line
    try { cmd.getValues(argc,argv);  } catch (int rv) { return rv; }

    // Initialize solution, rhs and geometry
    std::string solution_exact,rhs_exact;
    gsGeometry<>::Ptr geo;
    gsInfo << "|| Benchmark information ||" << endl;
    switch(numBenchmark)
    {
      case 1 : gsInfo << "CDR-equation the unit square" << endl;
               solution_exact = "sin(pi*x)*sin(pi*y)";
               rhs_exact = "1.2*pi*pi*sin(pi*x)*sin(pi*y)+0.9*pi*pi*sin(pi*x)*sin(pi*y)+0.7*pi*pi*cos(pi*x)*cos(pi*y) + 0.4*pi*pi*cos(pi*x)*cos(pi*y) +0.4*pi*cos(pi*x)*sin(pi*y)-0.2*pi*sin(pi*x)*cos(pi*y)+0.3*sin(pi*x)*sin(pi*y)";
               geo = gsNurbsCreator<>::BSplineSquare(1.0, 0.0, 0.0); break;

      case 2:  gsInfo << "Poisson equation on the quarter annulus (1)" << endl;
               solution_exact = "-(x*x+y*y-1)*(x*x+y*y-4)*x*y*y";
               rhs_exact = "2*x*(22*x*x*y*y+21*y*y*y*y-45*y*y+x*x*x*x-5*x*x+4)";
               geo = gsNurbsCreator<>::BSplineFatQuarterAnnulus(1.0, 2.0); break;

      case 3:  gsInfo << "Poisson equation on the quarter annulus (2)" << endl;
               solution_exact = "(x^2+y^2-3*sqrt(x^2+y^2)+2)*sin(2*atan(y/x))";
               rhs_exact = "(8-9*sqrt(x^2 + y^2))*sin(2*atan(y/x))/(x^2+y^2)";
               geo = gsNurbsCreator<>::BSplineFatQuarterAnnulus(1.0, 2.0); break;

      case 4:  gsInfo << "Poisson equation on an L-shaped domain" << endl;
               solution_exact = "if( y>0, ( (x^2+y^2)^(1.0/3.0) )*sin( (2*atan2(y,x) - pi)/3.0 ), ( (x^2+y^2)^(1.0/3.0) )*sin( (2*atan2(y,x) +3*pi)/3.0 ) )";
               rhs_exact = "0";
               geo = gsNurbsCreator<>::BSplineLShape_p1(); break;

      case 5:  gsInfo << "Poisson equation on the unit cube" << endl;
               solution_exact = "sin(pi*x)*sin(pi*y)*sin(pi*z)";
               rhs_exact = "(3*pi^2 )*sin(pi*x)*sin(pi*y)*sin(pi*z)";
               geo = gsNurbsCreator<>::BSplineCube(1); break;

      case 6: gsInfo << "Poisson's equation on Yeti footprint" << endl;
               solution_exact = "sin(5*pi*x)*sin(5*pi*y)";
               rhs_exact = "(50*pi^2 )*sin(5*pi*x)*sin(5*pi*y)";
               geo = gsReadFile<>("domain2d/yeti_mp2.xml"); break;
    }

    // Print information about benchmark
    gsInfo << "Exact solution: " << solution_exact << endl;
    gsInfo << "Right hand side: " << rhs_exact << endl;

    gsMultiPatch<> mp;
    switch(numPatches)
    {
      case 1:
      {
        gsMultiPatch<> mp0(*geo);
        mp = mp0; break;
      }
      case 2:
      {
        gsMultiPatch<> mp1(*geo);
        mp = mp1.uniformSplit(); break;
      }
      case 3:
      {
        gsMultiPatch<> mp2(*geo);
        gsMultiPatch<> mp3 = mp2.uniformSplit();
        mp = mp3.uniformSplit(); break;
      }
      case 4:
      {
        gsMultiPatch<> mp4(*geo);
        gsMultiPatch<> mp5 = mp4.uniformSplit();
        gsMultiPatch<> mp6 = mp5.uniformSplit();
        mp = mp6.uniformSplit(); break;
      }
    }

    // To read the entire Yeti footprint
    if(numBenchmark == 10)
    {
      std::string geometry("domain2d/yeti_mp2.xml");
      gsMultiPatch<>::uPtr mpPtr = gsReadFile<>(geometry);
      mp = *mpPtr;
    }
    gsInfo << "Number of patches: " << mp.nPatches() << "\n" << endl;
    gsFunctionExpr<> sol_exact(solution_exact,mp.geoDim());
    gsFunctionExpr<> f(rhs_exact,mp.geoDim());

    // Construct two bases (coarse and fine level)
    gsMultiBasis<> basisL(mp);
    gsMultiBasis<> basisH(mp);

    // Geometry of benchmark 10 is describe by quadratic B-splines
    if(numBenchmark == 10)
    {
      basisL.degreeReduce();
      basisH.degreeReduce();
    }
    // Geometry of benchmark 7 is describe by quadratic B-splines
    if(numBenchmark == 7)
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
    for( int i = 0; i < numLevels-1 ; ++i)
    {
      if( typeCoarsening[i] == 'h')
      {
        hp(i,0) = 1;
        numRefH = numRefH + 1;
      }
      else if( typeCoarsening[i] == 'p')
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
    if((numRefP + numRefZ) == numDegree)
    {
      basisL.degreeReduce(1);
    }
    else
    {
      basisL.degreeIncrease(numDegree-numRefP-numRefZ-1);
    }

    // Apply refinement in h for coarse and fine level
    for (int i = 0; i < numRefine - numRefH - numRefZ; ++i)
    {
      basisL.uniformRefine();
    }
    for (int i = 0; i < numRefine ; ++i)
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

    int iterTot = 1;
    int numCoarsening = numRefH+1;

    pMultigrid<real_t, gsSparseSolver<real_t>::LU , gsCDRAssembler<real_t> > My_MG(mp, basisL, bcInfo);
    gsMatrix<real_t> x;

    // Generate sequence of bases on all levels
    if(typeProjection == 1)
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
    if(numBenchmark == 1)
    {
      coeff_diff = gsFunctionExpr<>("1.2","-0.7","-0.4","0.9",2);
      coeff_conv = gsFunctionExpr<>("0.4","-0.2",2);
      coeff_reac = gsFunctionExpr<>("0.3",2);
    }


    // Apply p-Multigrid as stand-alone solver
    if(typeSolver == 1)
    {
          gsInfo << "p-multigrid is applied as stand-alone solver" << "\n" << endl;
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
    int i = 1;

    // Set up, determine initial L2 error
    gsField<> solKrylov;
    gsField<> sol = pa.constructSolution(x);
    real_t oldResNorm = r0.norm();

    // Perform BiCGStab
    if(typeSolver == 2)
    {
      gsInfo << "BiCGStab is applied as solver, p-multigrid as a preconditioner" << "\n" << endl;
      // Define vectors needed in BiCGStab
      gsVector<> t = gsVector<>::Zero(pa.matrix().rows()); gsVector<> s = gsVector<>::Zero(pa.matrix().rows()); gsVector<> p = gsVector<>::Zero(pa.matrix().rows()); gsVector<> v = gsVector<>::Zero(pa.matrix().rows()); gsMatrix<> y = gsMatrix<>::Zero(pa.matrix().rows(),1); gsMatrix<> z = gsMatrix<>::Zero(pa.matrix().rows(),1);
      real_t alp = 1; real_t rho = 1; real_t w = 1;

      // Set residual norm and #restarts
      real_t r0_sqnorm = r0.dot(r0);
      int restarts = 0;

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
        gsInfo << "BiCGStab iteration: " << i << "      |   Residual norm: "   << left << setw(15) << r.norm() << "           reduction:  1 / " << setprecision(3) << (oldResNorm/r.norm()) <<        setprecision  (6) << "\n";
        oldResNorm = r.norm();
        solKrylov = pa.constructSolution(x);

        ++i;
      }
  }
  else if(typeSolver == 3)
  {
      gsInfo << "CG is applied as solver, p-multigrid as a preconditioner" << "\n" << endl;

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
        gsInfo << "CG iteration: " << i << "       |  Residual norm: "   << left << setw(15) << r.norm() << "            reduction:  1 / " << setprecision(3) << (oldResNorm/r.norm()) <<  setprecision (6) << "\n";
        oldResNorm = r.norm();
        solKrylov = pa.constructSolution(x);
        ++i;
      }
  }
 return 0;
}

// Create the subspace corrected mass smoother
gsPreconditionerOp<>::Ptr setupSubspaceCorrectedMassSmoother(const gsSparseMatrix<>& matrix, const gsMultiBasis<>& mb, const gsBoundaryConditions<>& bc, const gsOptionList& opt, const int &typeBCHandling)
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
    for( index_t ps=0; ps < 2*dim; ++ps )
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
