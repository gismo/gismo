#include <gismo.h>
#include <string>

namespace gismo {

  /** @brief The p-multigrid base class provides the basic
   *  methods (smoothing, prolongation, restriction) for
   *  implementing p-multigrid methods
   */

  template<class T>
  struct gsXBraidMultigridBase
  {
  protected:
    int maxIter;
    int numLevels;
    int numSmoothing;
    int typeBCHandling;
    int typeCycle_h;
    int typeCycle_p;
    int typeLumping;
    int typeProjection;
    int typeSmoother;
    gsMatrix<> hp;
    T   tol;
    
  public:
    /// @brief Constructor
    gsXBraidMultigridBase()
      : maxIter(100000),
        numLevels(1),
        numSmoothing(1),
        typeBCHandling(1),
        typeCycle_h(2),
        typeCycle_p(1),
        typeLumping(1),
        typeProjection(1),
        typeSmoother(1),
        tol(1e-8)
    {}

    void setMaxIter(int maxIter)
    { this->maxIter = maxIter; }

    void setTolerance(T tol)
    { this->tol = tol; }
    
    void setNumLevels(int numLevels, int typeProjection, int numDegree)
    { 
     if(typeProjection == 1)
     {
      this->numLevels = numLevels - numDegree + 2;
     }
     else
     { 
      this->numLevels = numLevels; 
     }
    }

    void setNumSmoothing(int numSmoothing)
    { this->numSmoothing = numSmoothing; }

    void setTypeBCHandling(int typeBCHandling)
    { this->typeBCHandling = typeBCHandling; }

    void setTypeCycle_h(int typeCycle_h)
    { this->typeCycle_h = typeCycle_h; }

    void setTypeCycle_p(int typeCycle_p)
    { this->typeCycle_p = typeCycle_p; }

    void setTypeLumping(int typeLumping)
    { this->typeLumping = typeLumping; }

    void setTypeProjection(int typeProjection)
    { this->typeProjection = typeProjection; }

    void setTypeSmoother(int typeSmoother)
    { this->typeSmoother = typeSmoother; }

    void setCoarsening(gsMatrix<> hp)
    { this->hp = hp; }

    virtual gsXBraidMultigridBase& compute(const gsSparseMatrix<T>& mat, const T tstep, const int& numDegree, index_t typeMethod)
    {
      // Get arguments explicitly
      gsMatrix<T> x = gsMatrix<>::Zero(mat.rows(),1);
      gsMatrix<T> b = gsMatrix<>::Zero(mat.rows(),1);
      gsFunctionExpr<> rhs("1",2);
      int iterTot = 1;
      int typeMultigrid = 2;
      int typeCoarseOperator = 1;
      
      ///  @brief Set-up p-multigrid solver 
           setup(rhs,
                 x,
                 b,
                 iterTot,
                 numLevels,
                 numDegree,
                 typeMultigrid,
                 hp,
                 typeCoarseOperator,
                 tstep,
                 typeMethod);
 

      return *this; }

    virtual gsMatrix<T> solveWithGuess(const gsMatrix<T>& b,
                                       const gsMatrix<T>& x0)
    {
      // Get arguments explicitly
      gsMatrix<T> x(x0);
      x = x0;

      gsFunctionExpr<> rhs("1",2);
      int iterTot = 1;
      int typeMultigrid = 2;
      int typeCoarseOperator = 1;
      
      ///  @brief Apply p-multigrid solver to given right-hand side on level l
            solve(rhs,
                  x,
                  b,
                  iterTot,
                  numLevels,
                  typeMultigrid,
                  hp,
                  typeCoarseOperator);
      return x;
    }
    
    /// @brief Apply p-multigrid solver to given right-hand side on level l
    virtual void solveMG(const gsMatrix<T> & rhs,
                       std::vector<memory::shared_ptr<gsMultiBasis<T> > > m_bases,
                       gsMatrix<T>& x,
                       const int& numLevels,
                       gsBoundaryConditions<T> bcInfo,
                       gsMultiPatch<T> mp,
                       std::vector<gsSparseMatrix<T> >& m_prolongation_P,
                       std::vector<gsSparseMatrix<T> >& m_restriction_P,
                       std::vector<gsMatrix<T> >& m_prolongation_M,
                       std::vector<gsMatrix<T> >& m_restriction_M,
                       std::vector<gsSparseMatrix<T> >& m_prolongation_H,
                       std::vector<gsSparseMatrix<T> >& m_restriction_H,
                       const gsMatrix<T>& hp)
    {
      if ( numLevels == 1)
      {
          solvecoarse(rhs, x, numLevels);
          return;
      }

      if (hp(std::max(numLevels-2,0),0) == 0 )
        {
          gsMatrix<T> fineRes, coarseRes, fineCorr, coarseCorr, postRes;
          presmoothing(rhs, x, numLevels, fineRes, hp); 
          restriction(fineRes, coarseRes, numLevels, m_bases,
                      bcInfo, mp, 
                      m_prolongation_P, m_restriction_P,
                      m_prolongation_M, m_restriction_M,
                      m_prolongation_H, m_restriction_H, hp);
          //coarseRes.setZero(coarseRes.rows(),1);
          coarseCorr.setZero(coarseRes.rows(),1);
          for( int j = 0 ; j < (typeCycle_p == 2 ? 2 : 1) ; j++)
            {   
              solveMG(coarseRes, m_bases, coarseCorr, numLevels-1,
                    bcInfo, mp,
                    m_prolongation_P, m_restriction_P,
                    m_prolongation_M, m_restriction_M,
                    m_prolongation_H, m_restriction_H, hp);  
            }
          prolongation(coarseCorr, fineCorr, numLevels, m_bases,
                       bcInfo, mp,
                       m_prolongation_P, m_restriction_P,
                       m_prolongation_M, m_restriction_M,
                       m_prolongation_H, m_restriction_H, hp);
          postsmoothing(rhs, x, numLevels, fineCorr, postRes,
                        hp);
        }
   
      if (hp(std::max(numLevels-2,0),0) == 1 )
        {
          gsMatrix<T> fineRes, coarseRes, fineCorr, coarseCorr, postRes;
          presmoothing(rhs, x, numLevels, fineRes, hp); 
          restriction(fineRes, coarseRes, numLevels, m_bases,
                      bcInfo, mp,
                      m_prolongation_P, m_restriction_P,
                      m_prolongation_M, m_restriction_M,
                      m_prolongation_H, m_restriction_H, hp);
          //coarseRes.setZero(coarseRes.rows(),1);
          coarseCorr.setZero(coarseRes.rows(),1);
          for( int i = 0 ; i < (typeCycle_h == 2 ? 2 : 1) ; i++)
            {  
              solveMG(coarseRes, m_bases, coarseCorr, numLevels-1,
                    bcInfo, mp,
                    m_prolongation_P, m_restriction_P,
                    m_prolongation_M, m_restriction_M,
                    m_prolongation_H, m_restriction_H, hp);  
            }   
          prolongation(coarseCorr, fineCorr, numLevels, m_bases,
                       bcInfo, mp,
                       m_prolongation_P, m_restriction_P,
                       m_prolongation_M, m_restriction_M,
                       m_prolongation_H, m_restriction_H,  hp);
          postsmoothing(rhs,x, numLevels, fineCorr, postRes,
                        hp);
        }
    }

    virtual void setup(const gsFunctionExpr<T> & rhs,
               gsMatrix<T>& x,
               gsMatrix<T> f,
               const int& iterTot,
               const int& numLevels,
               const int& numDegree,
               const int& typeMultigrid,
               const gsMatrix<T>& hp,
               const int& typeCoarseOperator,
               T tstep,
               index_t typeMethod){}

   virtual void solve(const gsFunctionExpr<T> & rhs,
               gsMatrix<T>& x,
               gsMatrix<T> f,
               const int& iterTot,
               const int& numLevels,
               const int& typeMultigrid,
               const gsMatrix<T>& hp,
               const int& typeCoarseOperator){}
  
    /// @brief Apply fixed number of smoothing steps (pure virtual method)
    virtual void presmoothing(const gsMatrix<T>& rhs,
                              gsMatrix<T>& x,
                              const int& numLevels,
                              gsMatrix<T> & fineRes ,
                              const gsMatrix<T>& hp) = 0;

    /// @brief Apply fixed number of smoothing steps (pure virtual method)
    virtual void postsmoothing(const gsMatrix<T>& rhs,
                               gsMatrix<T>& x,
                               const int& numLevels,
                               gsMatrix<T> & fineCorr,
                               gsMatrix<T> & postRes,
                               const gsMatrix<T>& hp) = 0;

    /// @brief Apply coarse solver (pure virtual method)
    virtual void solvecoarse(const gsMatrix<T>& rhs,
                             gsMatrix<T>& x,
                             const int& numLevels) = 0;
   
    /// @brief Prolongate coarse space function to fine space
    virtual gsSparseMatrix<T> prolongation_P(const int& numLevels,
                                             std::vector<memory::shared_ptr<gsMultiBasis<T> > > m_bases) = 0;
   
    /// @brief Prolongate coarse space function to fine space
    virtual gsSparseMatrix<T> restriction_P(const int& numLevels,
                                            std::vector<memory::shared_ptr<gsMultiBasis<T> > > m_bases) = 0;
  
    /// @brief Prolongate coarse space function to fine space
    virtual gsMatrix<T> prolongation_M(const int& numLevels,
                                       std::vector<memory::shared_ptr<gsMultiBasis<T> > > m_bases) = 0;
   
    /// @brief Prolongate coarse space function to fine space
    virtual gsMatrix<T> restriction_M(const int& numLevels,
                                      std::vector<memory::shared_ptr<gsMultiBasis<T> > > m_bases) = 0;
 
    /// @brief Prolongate coarse space function to fine space
    virtual void prolongation(const gsMatrix<T>& Xcoarse,
                              gsMatrix<T>& Xfine,
                              const int& numLevels,
                              std::vector<memory::shared_ptr<gsMultiBasis<T> > > m_bases,
                              gsBoundaryConditions<T> bcInfo,
                              gsMultiPatch<T> mp,
                              std::vector<gsSparseMatrix<T> >& m_prolongation_P,
                              std::vector<gsSparseMatrix<T> >& m_restriction_P,
                              std::vector<gsMatrix<T> >& m_prolongation_M,
                              std::vector<gsMatrix<T> >& m_restriction_M,
                              std::vector<gsSparseMatrix<T> >& m_prolongation_H,
                              std::vector<gsSparseMatrix<T> >& m_restriction_H,
                              const gsMatrix<T>& hp)
    {
      if (hp(numLevels-2,0) == 1)
        {
          Xfine = m_prolongation_H[numLevels-2]*Xcoarse;
        }
      else
        {
          if (typeLumping == 1)
            {
              gsMatrix<T> temp = m_prolongation_P[numLevels-2]*Xcoarse;
              gsMatrix<T> M_L_inv = (m_prolongation_M[numLevels-2]).array().inverse();
              Xfine = (M_L_inv).cwiseProduct(temp);
            }
          else
            {
              // Define the low and high order bases
              gsMultiBasis<T> basesL = *m_bases[numLevels-2];
              gsMultiBasis<T> basesH = *m_bases[numLevels-1];
              typedef gsExprAssembler<real_t>::geometryMap geometryMap;
              typedef gsExprAssembler<real_t>::variable variable;
              typedef gsExprAssembler<real_t>::space space;
      
              // Determine matrix M (high_order * high_order)
              gsExprAssembler<real_t> ex2(1,1);
              geometryMap G2 = ex2.getMap(mp);
              space w_n = ex2.getSpace(basesH ,1, 0);
              w_n.setInterfaceCont(0);
              if (typeBCHandling == 1)
                {
		  w_n.setup(bcInfo, dirichlet::l2Projection, 0);
                  //#w_n.addBc(bcInfo.get("Dirichlet"));
                }
              ex2.setIntegrationElements(basesH);
              ex2.initSystem();
              ex2.assemble(w_n * meas(G2) * w_n.tr()); 
        
              // Prolongate Xcoarse to Xfine
              gsMatrix<T> temp = m_prolongation_P[numLevels-2]*Xcoarse;
              gsSparseMatrix<T> M = ex2.matrix();  
              gsConjugateGradient<T> CGSolver(M);
              CGSolver.setTolerance(1e-12);
              CGSolver.solve(temp,Xfine);        
            }
        }
    }

    /// @brief Restrict fine space function to coarse space
    virtual void restriction(const gsMatrix<T>& Xfine,
                             gsMatrix<T>& Xcoarse,
                             const int& numLevels,
                             std::vector<memory::shared_ptr<gsMultiBasis<T> > > m_bases,
                             gsBoundaryConditions<T> bcInfo,
                             gsMultiPatch<T> mp,
                             std::vector<gsSparseMatrix<T> >& m_prolongation_P,
                             std::vector<gsSparseMatrix<T> >& m_restriction_P,
                             std::vector<gsMatrix<T> >& m_prolongation_M,
                             std::vector<gsMatrix<T> >& m_restriction_M,
                             std::vector<gsSparseMatrix<T> >& m_prolongation_H,
                             std::vector<gsSparseMatrix<T> >& m_restriction_H,
                             const gsMatrix<T>& hp)
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
              gsMatrix<T> temp = m_restriction_P[numLevels-2]*Xfine;
              gsMatrix<T> M_L_inv = (m_restriction_M[numLevels-2]).array().inverse();
              Xcoarse = (M_L_inv).cwiseProduct(temp);
            }
          else
            {
              // Define the low and high order bases
              gsMultiBasis<T> basesL = *m_bases[numLevels-2];
              gsMultiBasis<T> basesH = *m_bases[numLevels-1];
              typedef gsExprAssembler<real_t>::geometryMap geometryMap;
              typedef gsExprAssembler<real_t>::variable variable;
              typedef gsExprAssembler<real_t>::space space;
      
              // Determine matrix M (low_order * low_order)
              gsExprAssembler<real_t> ex2(1,1);
              geometryMap G2 = ex2.getMap(mp);
              space w_n = ex2.getSpace(basesL, 1, 0);
              w_n.setInterfaceCont(0);
              if (typeBCHandling == 1)
                {
		  w_n.setup(bcInfo, dirichlet::l2Projection, 0);
                  //#w_n.addBc(bcInfo.get("Dirichlet"));
                }
              ex2.setIntegrationElements(basesL);
              ex2.initSystem();
              ex2.assemble(w_n * meas(G2) * w_n.tr());
        
              // Restrict Xfine to Xcoarse
              gsMatrix<T> temp = m_restriction_P[numLevels-2]*Xfine;
              gsSparseMatrix<T> M = ex2.matrix();  
              gsConjugateGradient<T> CGSolver(M);
              CGSolver.setTolerance(1e-12);
              CGSolver.solve(temp, Xcoarse);      
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
  template<class T, class CoarseSolver>
  struct gsXBraidMultigrid : public gsXBraidMultigridBase<T>
  {
  private:

    /// Base class type
    typedef gsXBraidMultigridBase<T> Base;

    /// Shared pointer to multi-patch geometry
    memory::shared_ptr<gsMultiPatch<T> > m_mp_ptr;

    /// Shared pointer to boundary conditions
    memory::shared_ptr<gsBoundaryConditions<T> > m_bcInfo_ptr;
 
    /// std::vector of multi-basis objects
    std::vector<memory::shared_ptr<gsMultiBasis<T> > > m_bases;

    /// std::vector of prolongation operators
    std::vector< gsSparseMatrix<T> > m_prolongation_P;

    /// std::vector of restriction operators
    std::vector< gsSparseMatrix<T> > m_restriction_P;
  
    /// std::vector of prolongation operators
    std::vector< gsMatrix<T> > m_prolongation_M;

    /// std::vector of restriction operators
    std::vector< gsMatrix<T> > m_restriction_M;

    /// std::vector of prolongation operators
    std::vector< gsSparseMatrix<T> > m_prolongation_H;

    /// std::vector of restriction operators
    std::vector< gsSparseMatrix<T> > m_restriction_H;

    /// std::vector of factorized operators
    std::vector< std::vector< gsSparseMatrix<T> > > m_ILUT;

    /// std::vector of factorized operators
    std::vector< std::vector < Eigen::PermutationMatrix<Dynamic,Dynamic,index_t> > > m_P;
  
    /// std::vector of factorized operators
    std::vector < std::vector < Eigen::PermutationMatrix<Dynamic,Dynamic,index_t> > > m_Pinv;

    /// std::vector of SCM smoother object
    std::vector< typename gsPreconditionerOp<T>::Ptr > m_SCMS;
  
    /// std::vector of operator objects
    std::vector< gsSparseMatrix<T> > m_operator;

    /// std::vector of std::vector of block operator objects
    std::vector < std::vector< gsSparseMatrix<T> > > m_block_operator;

    /// std::vector of std::vector of block operator objects
    std::vector < std::vector  < gsSparseMatrix<T> > > m_ddB;

    /// std::vector of std::vector of block operator objects
    std::vector < std::vector  < gsSparseMatrix<T> > > m_ddC;

    /// std::vector of std::vector of block operator objects
    std::vector < std::vector <  gsMatrix<T>  > > m_ddBtilde;

    /// std::vector of std::vector of block operator objects
    std::vector < std::vector <  gsMatrix<T>  > > m_ddCtilde;

    /// std::vector of std::vector of block operator objects
    std::vector <  gsMatrix<T> > m_A_aprox;

    /// std::vector of std::vector of block operator objects
    std::vector <  gsSparseMatrix<T> > m_S;
  
    /// std::vector of std::vector of shift objects
    std::vector < std::vector< int > > m_shift;

  public:

    // Constructor
    gsXBraidMultigrid(const gsMultiPatch<T> & mp,
                      const gsMultiBasis<T> & bases,
                      const gsBoundaryConditions<T> & bcInfo)
    {
      m_mp_ptr = memory::make_shared_not_owned(&mp);
      m_bcInfo_ptr = memory::make_shared_not_owned(&bcInfo);
      m_bases.push_back(memory::make_shared_not_owned(&bases));
    }

    virtual ~gsXBraidMultigrid() {}
    
  public:

    ///  @brief Set-up p-multigrid solver 
    void setup(const gsFunctionExpr<T> & rhs,
               gsMatrix<T>& x,
               gsMatrix<T> f,
               const int& iterTot,
               const int& numLevels,
               const int& numDegree,
               const int& typeMultigrid,
               const gsMatrix<T>& hp,
               const int& typeCoarseOperator,
               T tstep,
               index_t typeMethod)
    {
      for (int i = 1; i < numLevels; i++)
        {
          m_bases.push_back(give(m_bases.back()->clone()));
          switch((int) hp(i-1,0) )
            {
            case 0 : (Base::typeProjection == 1 ?
                      m_bases.back()->degreeIncrease(numDegree-1) :
                      m_bases.back()->degreeIncrease()); break;

            case 1 : m_bases.back()->uniformRefine();  break;

            case 2:  m_bases.back()->uniformRefine();
              m_bases.back()->degreeIncrease(); break;
            }
        }

     // Generate sequence of matrix K and M 
      m_operator.resize(numLevels);
      gsStopwatch clock;
      //gsInfo << "|| Multigrid hierarchy ||" <<std::endl;
     
      for (int i = 0; i < numLevels; i++)
      {  
        //gsInfo << "Level " << i+1 << " " ;
        //gsInfo << "Degree: " << m_bases[i]->degree() << ", Ndof: " << m_bases[i]->totalSize() <<std::endl; 

        typedef typename gsExprAssembler<T>::geometryMap geometryMap;
        typedef typename gsExprAssembler<T>::variable    variable;
        typedef typename gsExprAssembler<T>::space       space;
        typedef typename gsExprAssembler<T>::solution    solution;

        gsExprAssembler<T> K, M; 

        // Set the bases
        K.setIntegrationElements(*m_bases[i]);
        M.setIntegrationElements(*m_bases[i]);

        // Set the geometry map
        geometryMap G_K = K.getMap(*m_mp_ptr);
        geometryMap G_M = M.getMap(*m_mp_ptr);

        // Set the discretization space
        space u_K = K.getSpace(*m_bases[i]);
        space u_M = M.getSpace(*m_bases[i]);
        u_K.setInterfaceCont(0);
        u_M.setInterfaceCont(0);
	u_K.setup(*m_bcInfo_ptr, dirichlet::l2Projection, 0);
	u_M.setup(*m_bcInfo_ptr, dirichlet::l2Projection, 0);
        //#u_K.addBc( m_bcInfo_ptr->get("Dirichlet") );
        //#u_M.addBc( m_bcInfo_ptr->get("Dirichlet") );

        // Set the source term
        auto ff_K = K.getCoeff(rhs, G_K);
        auto ff_M = M.getCoeff(rhs, G_M);

        // Initialize and assemble the system matrix
        K.initSystem();
        K.assemble( igrad(u_K, G_K) * igrad(u_K, G_K).tr() * meas(G_K), u_K * ff_K * meas(G_K) );

        // Initialize and assemble the mass matrix
        M.initSystem();
        M.assemble( u_M * u_M.tr() * meas(G_M), u_M * ff_M * meas(G_M) );


        m_operator[i] = M.matrix() + tstep*K.matrix();
        switch(typeMethod)
          {
            case 0: m_operator[i] = M.matrix(); break;
            case 1: m_operator[i] = M.matrix() + tstep*K.matrix(); break;
            case 2: m_operator[i] = M.matrix() + 0.5*tstep*K.matrix();
          }

      }
      real_t Time_Assembly = clock.stop();
      GISMO_UNUSED(Time_Assembly);
 

      // Resize vector of operators
      m_prolongation_P.resize(numLevels-1);
      m_prolongation_M.resize(numLevels-1);
      m_prolongation_H.resize(numLevels-1);
      m_restriction_P.resize(numLevels-1);
      m_restriction_M.resize(numLevels-1);
      m_restriction_H.resize(numLevels-1);
    
      // Determine prolongation/restriction operators in p
      clock.restart();
      for (int i = 1; i < numLevels; i++)
        {
          if (hp(i-1,0) == 0)
            {
              m_prolongation_P[i-1] =  prolongation_P(i+1, m_bases);
              m_restriction_P[i-1] =  m_prolongation_P[i-1].transpose(); //restriction_P(i+1, m_bases);
              m_prolongation_M[i-1] =  prolongation_M(i+1, m_bases);
              m_restriction_M[i-1] = restriction_M(i+1, m_bases);
            }
        }

      // Determine prolongation/restriction operators in h
      gsSparseMatrix<real_t, RowMajor> transferMatrix;    
      gsOptionList options;
      Base::typeBCHandling == 1 ? options.addInt("DirichletStrategy","",dirichlet::elimination) : options.addInt("DirichletStrategy","",dirichlet::nitsche);
      for(int i = 1; i < numLevels; i++)
        {
          if (hp(i-1,0) == 1)
            {
              gsMultiBasis<T> m_bases_copy = *m_bases[i]; 
              m_bases_copy.uniformCoarsen_withTransfer(transferMatrix,*m_bcInfo_ptr,options); 
              m_prolongation_H[i-1] = transferMatrix;
              m_restriction_H[i-1] = m_prolongation_H[i-1].transpose();
            }  
        }
      real_t Time_Transfer = clock.stop();
      GISMO_UNUSED(Time_Transfer);
    
      // Obtain operators with Galerkin projection (TO DO)
      clock.restart();
      if (typeCoarseOperator == 2)
        {
          for (int i = numLevels-1; i > -1; i--)
            {
              if (hp(hp.rows()-1,0) == 0)
                {
                  if (hp(std::min(i,hp.rows()-1),0) == 1)
                    {
                      m_operator[i] = m_restriction_H[i]*m_operator[i+1]*m_prolongation_H[i];  
                    }
                }
              else
                {
                  if (hp(std::min(i,hp.rows()-1),0) == 1 && i > 0)
                    {
                      m_operator[i-1] = m_restriction_H[i-1]*m_operator[i]*m_prolongation_H[i-1];    
                    }
                }
            }
        }
      real_t Time_Assembly_Galerkin = clock.stop();
      GISMO_UNUSED(Time_Assembly_Galerkin);

      // Setting up the subspace corrected mass smoother
      clock.restart();
      if (Base::typeSmoother == 3)
        {
          // Generate sequence of SCM smoothers
          m_SCMS.resize(numLevels);
          gsOptionList opt;
          opt.addReal("Scaling","",0.12);
          for(int i = 0 ; i < numLevels ; i++)
            {
              m_SCMS[i] = setupSubspaceCorrectedMassSmoother(m_operator[i], *m_bases[i], *m_bcInfo_ptr, opt, Base::typeBCHandling);
            }
        }
      real_t Time_SCMS = clock.stop();
      GISMO_UNUSED(Time_SCMS);

      // Determine ILUT factorizations at each level
      clock.restart();  
      int numPatch = m_mp_ptr->nPatches();
      
      if (Base::typeSmoother == 1)
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
              if (Base::typeProjection == 2)
                {
                  Eigen::IncompleteLUT<real_t> ilu;
                  ilu.setFillfactor(1);
                  ilu.compute(m_operator[i]);
                  m_ILUT[i][0] = ilu.m_lu;
                  m_P[i][0] = ilu.m_P;
                  m_Pinv[i][0] = ilu.m_Pinv;
                }
              else
                {
                  if (i == numLevels-1) // Only at finest level
                    {
                      Eigen::IncompleteLUT<real_t> ilu;
                      ilu.setFillfactor(1);
                      ilu.compute(m_operator[i]);
                      m_ILUT[i][0] = ilu.m_lu;
                      m_P[i][0] = ilu.m_P;
                      m_Pinv[i][0] = ilu.m_Pinv;
                    }
                } 
            }
        } 
      real_t Time_ILUT_Factorization = clock.stop();
      GISMO_UNUSED(Time_ILUT_Factorization);
      
      clock.restart();   
      if (Base::typeSmoother == 5)
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
              //m_bases[i]->partition(interior,boundary,interface,global_interior,global_boundary,global_interface);
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
                  const gsSparseMatrix<T> block = m_operator[i].block(shift0,shift0,m_shift[i][j],m_shift[i][j]);
                  Eigen::IncompleteLUT<real_t> ilu;
                  ilu.setFillfactor(1);
                  ilu.compute(block);
                  m_ILUT[i][j] = ilu.m_lu;

                  m_P[i][j] = ilu.m_P;
                  m_Pinv[i][j] = ilu.m_Pinv;
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
              m_A_aprox[i] = gsSparseMatrix<T>(m_operator[i].rows(),m_operator[i].cols());

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
                  m_ddBtilde[i][j] = gsSparseMatrix<T>(m_shift[i][j],m_shift[i][numPatch]);
                  m_ddCtilde[i][j] = gsSparseMatrix<T>(m_shift[i][j],m_shift[i][numPatch]);
                  for(int k=0 ; k < m_shift[i][numPatch]; k++)
                    {
                      gsMatrix<T> Brhs = m_ddC[i][j].col(k);
                      gsMatrix<T> Crhs = m_ddC[i][j].col(k);
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
        
              // Perform ILUT on the S-matrix!
              Eigen::IncompleteLUT<real_t> ilu;
              ilu.setFillfactor(1);
              gsSparseMatrix<T> II = m_S[i];
              ilu.compute(II);
              m_A_aprox[i].block(m_A_aprox[i].rows() - m_shift[i][numPatch],m_A_aprox[i].rows() - m_shift[i][numPatch],m_shift[i][numPatch],m_shift[i][numPatch]) = ilu.m_lu;
            }
        }
      
      real_t Time_Block_ILUT_Factorization = clock.stop();
      GISMO_UNUSED(Time_Block_ILUT_Factorization);
      
      // gsInfo << "\n|| Setup Timings || " <<std::endl;
      // gsInfo << "Total Assembly time: " << Time_Assembly <<std::endl;
      // gsInfo << "Total ILUT factorization time: " << Time_ILUT_Factorization <<std::endl;
      // gsInfo << "Total block ILUT factorization time: " << Time_Block_ILUT_Factorization <<std::endl;
      // gsInfo << "Total SCMS time: " << Time_SCMS <<std::endl;
      // gsInfo << "Total setup time: " << Time_Assembly_Galerkin + Time_Assembly + Time_Transfer + Time_ILUT_Factorization + Time_SCMS <<std::endl;
    }

    ///  @brief Apply p-multigrid solver to given right-hand side on level l
    void solve(const gsFunctionExpr<T> & rhs,
               gsMatrix<T>& x,
               gsMatrix<T> f,
               const int& iterTot,
               const int& numLevels,
               const int& typeMultigrid,
               const gsMatrix<T>& hp,
               const int& typeCoarseOperator)
    {
      gsStopwatch clock;
      gsMatrix<T> b = f;

      // Determine residual and L2 error
      real_t r0 = (m_operator[numLevels-1]*x - b).norm();
      real_t r = r0;
      int iter = 1;

      // Solve with p-multigrid method 
      real_t r_old = r0;
      clock.restart();
      // Adjusted stopping criterion!!
       while( r/b.norm() > Base::tol && iter < Base::maxIter )   
        {
          // Call solver from base class
          Base::solveMG(b, m_bases, x, numLevels,
                      *m_bcInfo_ptr, *m_mp_ptr,
                      m_prolongation_P, m_restriction_P,
                      m_prolongation_M, m_restriction_M,
                      m_prolongation_H, m_restriction_H, hp);
          
          r = (m_operator[numLevels-1]*x - b).norm();
          if ( r_old < r)
            {
              gsInfo << "Residual increased during solving!!! " <<std::endl;
            }
          r_old = r;
          //gsInfo << "Residual after cycle " << iter << " equals: " << r <<std::endl;
          iter++;
        }
      real_t Time_Solve = clock.stop();
      GISMO_UNUSED(Time_Solve);

    // gsInfo << "\n|| Solver information || " <<std::endl;
    // gsInfo << "Solver converged in " << Time_Solve << " seconds!" <<std::endl;
    // gsInfo << "Solver converged in: " << iter-1 << " iterations!" <<std::endl;           
  }

  private:

    /// @brief Apply coarse solver
    virtual void solvecoarse(const gsMatrix<T>& rhs,
                             gsMatrix<T>& x,
                             const int& numLevels)
    {
      //Direct solver (LU factorization)
      CoarseSolver solver;
      solver.analyzePattern(m_operator[numLevels-1]);
      solver.factorize(m_operator[numLevels-1]);
      x = solver.solve(rhs); 
    }
  
    /// @brief Construct prolongation operator at level numLevels
    virtual gsMatrix<T> prolongation_M(const int& numLevels,
                                       std::vector<memory::shared_ptr<gsMultiBasis<T> > > m_bases)
    {
      // Define the low and high order bases
      gsMultiBasis<T> basesL = *m_bases[numLevels-2];
      gsMultiBasis<T> basesH = *m_bases[numLevels-1];

      // Determine matrix M (high_order * high_order)
      typedef gsExprAssembler<real_t>::geometryMap geometryMap;
      typedef gsExprAssembler<real_t>::variable variable;
      typedef gsExprAssembler<real_t>::space    space;
      gsExprAssembler<real_t> ex2(1,1);
      geometryMap G2 = ex2.getMap(*m_mp_ptr);
      space w_n = ex2.getSpace(basesH ,1, 0);
      w_n.setInterfaceCont(0);
      if (Base::typeBCHandling == 1)
        {
	  w_n.setup(*m_bcInfo_ptr, dirichlet::l2Projection, 0);
          //#w_n.addBc(m_bcInfo_ptr->get("Dirichlet"));
        }
      ex2.setIntegrationElements(basesH);
      ex2.initSystem();
      ex2.assemble(w_n * meas(G2) );   
      return ex2.rhs();
    }

    /// @brief Construct prolongation operator at level numLevels
    virtual gsSparseMatrix<T> prolongation_P(const int& numLevels,
                                             std::vector<memory::shared_ptr<gsMultiBasis<T> > > m_bases)
    {
      // Define the low and high order bases
      gsMultiBasis<T> basesL = *m_bases[numLevels-2];
      gsMultiBasis<T> basesH = *m_bases[numLevels-1];

      // Determine matrix P (high_order * low_order)
      typedef gsExprAssembler<real_t>::geometryMap geometryMap;
      gsExprAssembler<real_t> ex(1,1);
      geometryMap G = ex.getMap(*m_mp_ptr);
      typedef gsExprAssembler<real_t>::variable variable;
      typedef gsExprAssembler<real_t>::space    space;
      space v_n = ex.getSpace(basesH ,1, 0);
      v_n.setInterfaceCont(0);
      space u_n = ex.getTestSpace(v_n , basesL);
      u_n.setInterfaceCont(0);
      if (Base::typeBCHandling == 1)
        {
	  v_n.setup(*m_bcInfo_ptr, dirichlet::l2Projection, 0);
	  u_n.setup(*m_bcInfo_ptr, dirichlet::l2Projection, 0);
          //#v_n.addBc(m_bcInfo_ptr->get("Dirichlet"));
          //#u_n.addBc(m_bcInfo_ptr->get("Dirichlet"));
        }
      ex.setIntegrationElements(basesH);
      ex.initSystem();
      ex.assemble(u_n*meas(G) * v_n.tr()); 
      gsSparseMatrix<T> P = ex.matrix().transpose();
      return P;    
    }

    /// @brief Construct restriction operator at level numLevels
    virtual gsMatrix<T> restriction_M(const int& numLevels,
                                      std::vector<memory::shared_ptr<gsMultiBasis<T> > > m_bases)
    {
      // Define the low and high order bases
      gsMultiBasis<T> basesL = *m_bases[numLevels-2];
      gsMultiBasis<T> basesH = *m_bases[numLevels-1];
      
      // Determine matrix M (low_order * low_order)
      typedef gsExprAssembler<real_t>::geometryMap geometryMap;
      typedef gsExprAssembler<real_t>::variable variable;
      typedef gsExprAssembler<real_t>::space    space;
      gsExprAssembler<real_t> ex2(1,1);
      geometryMap G2 = ex2.getMap(*m_mp_ptr);
      space w_n = ex2.getSpace(basesL ,1, 0);
      w_n.setInterfaceCont(0);
      if (Base::typeBCHandling == 1)
        {
	  w_n.setup(*m_bcInfo_ptr, dirichlet::l2Projection, 0);
          //#w_n.addBc(m_bcInfo_ptr->get("Dirichlet"));
        }
      ex2.setIntegrationElements(basesL);
      ex2.initSystem();
      ex2.assemble(w_n * meas(G2) ); 
      return ex2.rhs();     
    }

    /// @brief Construct restriction operator at level numLevels
    virtual gsSparseMatrix<T> restriction_P(const int& numLevels,
                                            std::vector<memory::shared_ptr<gsMultiBasis<T> > > m_bases)
    {
      // Define the low and high order bases
      gsMultiBasis<T> basesL = *m_bases[numLevels-2];
      gsMultiBasis<T> basesH = *m_bases[numLevels-1];
      
      // Determine matrix P (high_order * low_order)
      gsExprAssembler<real_t> ex(1,1);
      typedef gsExprAssembler<real_t>::geometryMap geometryMap;
      geometryMap G = ex.getMap(*m_mp_ptr);
      
      typedef gsExprAssembler<real_t>::variable variable;
      typedef gsExprAssembler<real_t>::space    space;
      space v_n = ex.getSpace(basesH ,1, 0);
      v_n.setInterfaceCont(0);
      space u_n = ex.getTestSpace(v_n , basesL);
      u_n.setInterfaceCont(0);
      if (Base::typeBCHandling == 1)
        {
	  u_n.setup(*m_bcInfo_ptr, dirichlet::l2Projection, 0);
	  v_n.setup(*m_bcInfo_ptr, dirichlet::l2Projection, 0);
          //#u_n.addBc(m_bcInfo_ptr->get("Dirichlet"));
          //#v_n.addBc(m_bcInfo_ptr->get("Dirichlet"));
        }
      ex.setIntegrationElements(basesH);
      ex.initSystem();
      ex.assemble(u_n * meas(G)* v_n.tr());
      gsSparseMatrix<T> P = ex.matrix();
      return P;   
    }
  
    /// @brief Apply fixed number of presmoothing steps
    virtual void presmoothing(const gsMatrix<T>& rhs,
                              gsMatrix<T>& x,
                              const int& numLevels,
                              gsMatrix<T> & fineRes,
                              const gsMatrix<T>& hp)
    { 
      //gsInfo << "Residual before presmoothing: " << (rhs-m_operator[numLevels-1]*x).norm() << " at level " << numLevels <<std::endl;
      for(int i = 0 ; i < Base::numSmoothing ; i++)
        {
          if (Base::typeSmoother == 1)
            {
              if (hp(numLevels-2,0) == 1 && hp(hp.rows()-1,0) == 0)
                {
                  internal::gaussSeidelSweep(m_operator[numLevels-1],x,rhs);
                }
              else
                {   
                  gsMatrix<T> e; 
                  gsMatrix<T> d = rhs-m_operator[numLevels-1]*x;
                  e = m_Pinv[numLevels-1][0]*d;
                  e = m_ILUT[numLevels-1][0].template triangularView<Eigen::UnitLower>().solve(e);
                  e = m_ILUT[numLevels-1][0].template triangularView<Eigen::Upper>().solve(e);
                  e = m_P[numLevels-1][0]*e;
                  x = x + e; 
                }     
            }
          if (Base::typeSmoother == 2)
            {
              internal::gaussSeidelSweep(m_operator[numLevels-1],x,rhs); 
            }
          if (Base::typeSmoother == 3)
            {  
              m_SCMS[numLevels-1]->step(rhs,x);
            }
          if (Base::typeSmoother == 5)
            {
              if (hp(numLevels-2,0) == 1 && hp(hp.rows()-1,0) == 0)
                {
                  internal::gaussSeidelSweep(m_operator[numLevels-1],x,rhs);
                }
              else
                {
                  gsMatrix<T> e; 
                  gsMatrix<T> d = rhs-m_operator[numLevels-1]*x;
                  e = m_A_aprox[numLevels-1].template triangularView<Eigen::UnitLower>().solve(d);
                  e = m_A_aprox[numLevels-1].template triangularView<Eigen::Upper>().solve(e);
                  x = x + e;   
                }     
            }
        }          
      //   gsInfo << "Residual after presmoothing: " << (rhs-m_operator[numLevels-1]*x).norm() << " at level " << numLevels <<std::endl;
      fineRes = m_operator[numLevels-1]*x - rhs;
    }

    /// @brief Apply fixed number of postsmoothing steps
    virtual void postsmoothing(const gsMatrix<T>& rhs,
                               gsMatrix<T>& x,
                               const int& numLevels,
                               gsMatrix<T> & fineCorr,
                               gsMatrix<T> & postRes,
                               const gsMatrix<T>& hp)
    {
      real_t alpha = 1;
      x = x - alpha*fineCorr;
      //gsInfo << "Residual before postsmoothing: " << (rhs-m_operator[numLevels-1]*x).norm() << " at level " << numLevels <<std::endl;

      for(int i = 0 ; i < Base::numSmoothing ; i++)
        {
          if (Base::typeSmoother == 1)
            {
              if (hp(numLevels-2,0) == 1 && hp(hp.rows()-1,0) == 0)
                { 
                   internal::gaussSeidelSweep(m_operator[numLevels-1],x,rhs);      
                }
              else
                { 
                  gsMatrix<T> e; 
                  gsMatrix<T> d = rhs-m_operator[numLevels-1]*x;
                  e = m_Pinv[numLevels-1][0]*d;
                  e = m_ILUT[numLevels-1][0].template triangularView<Eigen::UnitLower>().solve(e);
                  e = m_ILUT[numLevels-1][0].template triangularView<Eigen::Upper>().solve(e);
                  e = m_P[numLevels-1][0]*e;
                  x = x + e;
                }
            }
          if (Base::typeSmoother == 2)
            {
                internal::gaussSeidelSweep(m_operator[numLevels-1],x,rhs);      
            }
          if (Base::typeSmoother == 3)
            {
              m_SCMS[numLevels-1]->step(rhs,x);
            }
          if (Base::typeSmoother == 5)
            {
              if (hp(numLevels-2,0) == 1 && hp(hp.rows()-1,0) == 0)
                {
                  internal::gaussSeidelSweep(m_operator[numLevels-1],x,rhs); 
                }
              else
                { 
                  gsMatrix<T> e; 
                  gsMatrix<T> d = rhs-m_operator[numLevels-1]*x;
                  e = m_A_aprox[numLevels-1].template triangularView<Eigen::UnitLower>().solve(d);
                  e = m_A_aprox[numLevels-1].template triangularView<Eigen::Upper>().solve(e);
                  x = x + e;    
                }
            }  
          postRes = rhs - m_operator[numLevels-1]*x;        
          //      gsInfo << "Residual after postsmoothing: " << (rhs-m_operator[numLevels-1]*x).norm() << " at level " << numLevels <<std::endl;
        }
    }
  };

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


} // namespace gismo
