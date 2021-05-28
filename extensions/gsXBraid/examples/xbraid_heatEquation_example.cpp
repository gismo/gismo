/** @file xbraid_example.cpp

    @brief XBraid integration

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris, M. Moeller
*/

#include <gismo.h>
#include <gsXBraid/gsXBraid.h>
#include "gsXBraidMultigrid.h"

using namespace gismo;

#ifdef GISMO_WITH_XBRAID

namespace gismo {

  enum class gsXBraid_typeMethod
  {
    FE_FE = 0, // forward Euler   (all grids)
    BE_BE = 1, // backward Euler  (all grids)
    CN_CN = 2, // Crank-Nicholson (all grids)
    FE_BE = 3, // forward Euler   (fine grid), backward Euler (coarser grids)
    CN_BE = 4  // Crank-Nicholson (fine grid), backward Euler (coarser grids)      
  };

/**
   \brief Derived class implementing the XBraid wrapper for the heat equation
*/
template<typename T>
class gsXBraid_app : public gsXBraid< gsMatrix<T> >
{
private:
  // Spatial discretisation parameters
  index_t numRefine, numElevate, numIncrease;

  // Temporal discretisation parameters
  index_t numSteps, typeMethod;
  T tstart, tstop, tstep;

  // Spatial discretizations
  gsMultiPatch<T> mp;
  gsMultiBasis<T> basesH, basesL;
  
  // Boundary conditions
  gsBoundaryConditions<T> bc;

  // Assembler options
  gsOptionList Aopt, Sopt, Topt;
  
  // Expression assembler
  gsExprAssembler<T> K, M;
  gsFunctionExpr<T> f, u0, ms;
  
  // Solution
  gsMatrix<T> sol;
  
  // Multigrid solver
  typedef gsXBraidMultigrid<T, typename gsSparseSolver<T>::LU > solver_mg;
  std::vector< solver_mg* > m_solver;    
    
  typedef typename gsExprAssembler<T>::geometryMap geometryMap;
  typedef typename gsExprAssembler<T>::variable    variable;
  typedef typename gsExprAssembler<T>::space       space;
  typedef typename gsExprAssembler<T>::solution    solution;
  
 public:
  /// Contructor
  gsXBraid_app(const gsMpiComm& comm,
               const T&         tstart,
               const T&         tstop,
               index_t          typeMethod,
               index_t          numSteps,
               index_t          numRefine,
               index_t          numElevate,
               index_t          numIncrease,
               std::string&     fn)
    : gsXBraid< gsMatrix<T> >::gsXBraid(comm, tstart, tstop, (int)numSteps),
      numRefine(numRefine),
      numElevate(numElevate),
      numIncrease(numIncrease),
      numSteps(numSteps),
      typeMethod(typeMethod),
      tstart(tstart),
      tstop(tstop),
      tstep( (tstop-tstart)/numSteps ),
      K(1,1), M(1,1)
  {
    /////////////////////////////////////////////////////////////////////////////////////////////
    //                           Code for heat equation starts here                            //
    /////////////////////////////////////////////////////////////////////////////////////////////
    
    gsFileData<T> fd(fn);
    if (this->id() == 0) gsInfo << "Loaded file " << fd.lastPath() << "\n";

    fd.getId(0, mp); // id=0: Multipatch domain
    basesH = gsMultiBasis<T>(mp);
    basesL = gsMultiBasis<T>(mp);
    
    fd.getId(1, f); // id=1: right-hand side function
    if (this->id() == 0) gsInfo << "Source function " << f << "\n";
    
    fd.getId(2, bc); // id=2: boundary conditions
    if (this->id() == 0) gsInfo << "Boundary conditions:\n" << bc << "\n";

    fd.getId(3, u0); // id=3: initial conditions
    if (this->id() == 0) gsInfo << "Initial conditions:\n" << u0 << "\n";
    
    fd.getId(4, ms); // id=4: manufactured solution
    if (this->id() == 0) gsInfo << "Manufactured solution:\n" << ms << "\n";
    
    fd.getId(5, Aopt); // id=5: assembler options
    if (this->id() == 0) gsInfo << "Assembler options:\n" << Aopt << "\n";
    K.setOptions(Aopt);
    M.setOptions(Aopt);

    fd.getId(6, Topt); // id=6: multigrid-in-time options
    if (this->id() == 0) gsInfo << "Multigrid-in-time options:\n" << Topt << "\n";

    this->SetCFactor(Topt.getInt("CFactor"));
    this->SetMaxIter(Topt.getInt("maxIter"));
    this->SetMaxLevels(Topt.getInt("maxLevel"));
    this->SetMaxRefinements(Topt.getInt("numMaxRef"));
    this->SetMinCoarse(Topt.getInt("minCLevel"));
    this->SetNFMG(Topt.getInt("numFMG"));
    this->SetNFMGVcyc(Topt.getInt("numFMGVcyc"));
    this->SetNRelax(Topt.getInt("numRelax"));
    this->SetAccessLevel(Topt.getInt("access"));
    this->SetPrintLevel(Topt.getInt("print"));
    this->SetStorage(Topt.getInt("numStorage"));
    this->SetTemporalNorm(Topt.getInt("norm"));

    if (Topt.getSwitch("fmg"))           this->SetFMG();
    if (Topt.getSwitch("incrMaxLevels")) this->SetIncrMaxLevels();
    if (Topt.getSwitch("periodic"))      this->SetPeriodic(1); else this->SetPeriodic(0);
    if (Topt.getSwitch("refine"))        this->SetRefine(1);   else this->SetRefine(0);
    if (Topt.getSwitch("sequential"))    this->SetSeqSoln(1);  else this->SetSeqSoln(0);
    if (Topt.getSwitch("skip"))          this->SetSkip(1);     else this->SetSkip(0);
    if (Topt.getSwitch("spatial"))       this->SetSpatialCoarsenAndRefine();
    if (Topt.getSwitch("tol"))           this->SetAbsTol(Topt.getReal("absTol"));
    else                                 this->SetRelTol(Topt.getReal("relTol"));    
    
    fd.getId(7, Sopt); // id=6: spatial solver options
    if (this->id() == 0) gsInfo << "Spatial solver options:\n" << Sopt << "\n";
        
    std::string typeCoarsening = Sopt.getString("coarseStrategy");
    gsMatrix<> hp = gsMatrix<>::Zero(Sopt.getInt("numLevels")-1,1);

    // Read string from command line
    real_t numRefH = 0;
    real_t numRefP = 0;
    real_t numRefZ = 0;
    
    // Convert input string to array
    for( int i = 0; i < Sopt.getInt("numLevels")-1 ; ++i)
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
    if((numRefP + numRefZ) == numIncrease )
    {
      basesL.degreeReduce(1);
    }    
    else
    {
      basesL.degreeIncrease(numIncrease-numRefP-numRefZ-1); 
    }
 
    // Apply refinement in h for coarse and fine level
    for (int i = 0; i < numRefine - numRefH - numRefZ; ++i)
    {
      basesL.uniformRefine();
    }
    for (int i = 0; i < numRefine ; ++i)
    {
      basesH.uniformRefine();
    }

    // Apply refinement in p for fine level
    basesH.degreeIncrease(numIncrease-1);

    // Set the bases
    K.setIntegrationElements(basesH);
    M.setIntegrationElements(basesH);   

    // Set the geometry map
    geometryMap G_K = K.getMap(mp);
    geometryMap G_M = M.getMap(mp);

    // Set the discretization space
    space u_K = K.getSpace(basesH);
    space u_M = M.getSpace(basesH);
    u_K.setInterfaceCont(0);
    u_M.setInterfaceCont(0);
    u_K.addBc( bc.get("Dirichlet") );
    u_M.addBc( bc.get("Dirichlet") );

    // Set the source term
    variable ff_K = K.getCoeff(f, G_K);
    variable ff_M = M.getCoeff(f, G_M);

    // Initialize and assemble the system matrix
    K.initSystem();
    K.assemble( igrad(u_K, G_K) * igrad(u_K, G_K).tr() * meas(G_K), u_K * ff_K * meas(G_K) );

    // Initialize and assemble the mass matrix
    M.initSystem();
    M.assemble( u_M * u_M.tr() * meas(G_M), u_M * ff_M * meas(G_M) );

    // Enforce Neumann conditions to right-hand side
    variable g_Neumann = K.getBdrFunction();
    K.assembleRhsBc(u_K * g_Neumann.val() * nv(G_K).norm(), bc.neumannSides() );

    // Determine MGRIT levels a priori
    int numMGRITLevels = 1;
    int StepsLevel = numSteps;
    for(int i = 1 ; i < 10000; i++){
        StepsLevel = StepsLevel/Topt.getInt("CFactor"); 
        if(StepsLevel < Topt.getInt("minCLevel")) 
          break;
        numMGRITLevels = numMGRITLevels + 1;
    } 
    
    m_solver.resize(numMGRITLevels);
    real_t tstep_level = tstep;
    for(int i = 0 ; i < numMGRITLevels ; i++)
    {
      m_solver[i] = new solver_mg(mp, basesL, bc);
      m_solver[i]->setMaxIter(Sopt.getInt("maxIter"));
      m_solver[i]->setTolerance(Sopt.getReal("tol"));
      m_solver[i]->setNumLevels(Sopt.getInt("numLevels"),Sopt.getInt("projection"),numIncrease);
      m_solver[i]->setNumSmoothing(Sopt.getInt("numSmoothing"));
      m_solver[i]->setTypeBCHandling(Sopt.getInt("bcHandling"));
      m_solver[i]->setTypeCycle_h(Sopt.getInt("cycle_h"));
      m_solver[i]->setTypeCycle_p(Sopt.getInt("cycle_p"));
      m_solver[i]->setTypeLumping(Sopt.getInt("lumping"));
      m_solver[i]->setTypeProjection(Sopt.getInt("projection"));
      m_solver[i]->setTypeSmoother(Sopt.getInt("smoother"));
      m_solver[i]->setCoarsening(hp);
      if(typeMethod > 2 && i == 0)
      {
        m_solver[i]->compute(M.matrix(),tstep_level,numIncrease,typeMethod);
      }
      else
      {
        // Apple Backward Euler at coarser levels (FE_BE and CN_BE)
        m_solver[i]->compute(M.matrix(),tstep_level,numIncrease,1);  
      }
      tstep_level = tstep_level*Topt.getInt("CFactor");
    }
    
    if (this->id() == 0) {
 
      gsStopwatch clock;
      clock.restart();
      
      sol.setZero(M.numDofs());
 
      switch((gsXBraid_typeMethod)typeMethod) {
      case gsXBraid_typeMethod::FE_FE:
      case gsXBraid_typeMethod::FE_BE:
        // Forward Euler method
        
        for ( int i = 1; i<=numSteps; ++i) // for all timesteps
          // Compute the system for the timestep i (rhs is assumed constant wrt time)
          sol = m_solver[0]->solveWithGuess(tstep*K.rhs() +
                                                   (M.matrix()-tstep*K.matrix())*sol,
                                                   sol);
        break;
        
      case gsXBraid_typeMethod::BE_BE:
        // Backward Euler method
        
       for ( int i = 1; i<=numSteps; ++i) // for all timesteps
          // Compute the system for the timestep i (rhs is assumed constant wrt time)
          sol = m_solver[0]->solveWithGuess(tstep*K.rhs() + (M.matrix())*sol, sol);
       break;
        
      case gsXBraid_typeMethod::CN_CN:
      case gsXBraid_typeMethod::CN_BE:
        // Crank-Nicholson method
        for ( int i = 1; i<=numSteps; ++i) // for all timesteps
          // Compute the system for the timestep i (rhs is assumed constant wrt time)
          sol = m_solver[0]->solveWithGuess(tstep*K.rhs() +
                                                   (M.matrix()-tstep*0.5*K.matrix())*sol,
                                                   sol);
        break;
        
      default:
        throw std::runtime_error("Unsupported time-stepping method");
      }
      
      gsInfo << "wall time = " << clock.stop() << "\n"
             << "L2 norm of the solution  = " << sol.norm() << "\n";             
      
      // gsExprEvaluator<T> ev(M);
      // solution u_sol = M.getSolution(u_M, sol);
      // variable u_ex  = ev.getVariable(ms, G_M);
      // T l2err = math::sqrt( ev.integral( (u_ex - u_sol).sqNorm() * meas(G_M) ) );
      // T h1err = l2err +
      //   math::sqrt(ev.integral( ( igrad(u_ex) - grad(u_sol)*jac(G_M).inv() ).sqNorm() * meas(G_M) ));

      // gsInfo << "L2 error of the solution = " << l2err << "\n"
      //        << "H1 error of the solution = " << h1err << std::flush;
    }
  }

  /// Destructor
  virtual ~gsXBraid_app()
  {

  }
  
  /// Creates instance from command line argument
  static inline gsXBraid_app create(const gsMpiComm& comm,
                                    int              argc,
                                    char**           argv)
  {
    // Problem parameters
    std::string fn(XBRAID_DATA_DIR"pde/heat2d_square_ibvp1.xml");
    
    // Spatial discretisation parameters
    index_t numRefine     = 2;
    index_t numElevate    = 0;
    index_t numIncrease   = 0;
    
    // Temporal discretisation parameters
    index_t numSteps      = 40;
    index_t typeMethod    = (index_t)gsXBraid_typeMethod::BE_BE;
    T       tfinal        = 0.1;
        
    gsCmdLine cmd("Tutorial on solving a Heat equation problem using parallel-in-time multigrid.");

    // Problem parameters
    cmd.addString( "f", "file", "Input XML file", fn );
    
    // Spatial discretisation parameters
    cmd.addInt( "e", "degreeElevation",
                "Number of degree elevation steps to perform before solving (0: equalize degree in all directions)", numElevate );
    cmd.addInt( "i", "degreeIncrease",
                "Number of degree increase steps to perform before solving (0: equalize degree in all directions)", numIncrease );
    cmd.addInt( "r", "uniformRefine", "Number of uniform h-refinement steps to perform before solving",  numRefine );

    // Temporal diescretisation parameters
    cmd.addInt( "n",  "numSteps", "Number of parallel-in-time steps", numSteps );
    cmd.addInt( "T", "typeMethod", "Time-stepping scheme", typeMethod);
    cmd.addReal( "t", "tfinal", "Final time", tfinal );
    
    cmd.getValues(argc,argv);

    // Create instance
    gsXBraid_app<T> app(comm, 0.0, tfinal, typeMethod, numSteps, numRefine, numElevate, numIncrease, fn);
    
    return app;
  }

  /// Initializes a vector
  braid_Int Init(braid_Real    t,
                 braid_Vector *u_ptr)
#if __cplusplus >= 201103L || _MSC_VER >= 1600
    override
#endif
  {
    gsMatrix<T>* u = new gsMatrix<T>(M.numDofs(), 1);
    
    if (t != tstart) {
      // Intermediate solution
      u->setZero(M.numDofs());
    } else {
      // Initial solution
      u->setZero(M.numDofs());
    }

    *u_ptr = (braid_Vector) u;
    return braid_Int(0);
  }
  
  /// Performs a single step of the parallel-in-time multigrid
  braid_Int Step(braid_Vector    u,
                 braid_Vector    ustop,
                 braid_Vector    fstop,
                 BraidStepStatus &status)
#if __cplusplus >= 201103L || _MSC_VER >= 1600
    override
#endif
  {
    gsMatrix<T>* u_ptr = (gsMatrix<T>*) u;
    gsMatrix<T>* ustop_ptr = (gsMatrix<T>*) ustop;

    // XBraid forcing
    if (fstop != NULL) {
      gsMatrix<T>* fstop_ptr = (gsMatrix<T>*) fstop;
      *u_ptr += *fstop_ptr;
    }
    
    // Get time step information
    std::pair<braid_Real, braid_Real> time =
      static_cast<gsXBraidStepStatus&>(status).timeInterval();
    T tstep(time.second - time.first);

    switch((gsXBraid_typeMethod)typeMethod) {
    case gsXBraid_typeMethod::FE_FE:
      // Forward Euler method (all grids)
      *u_ptr = m_solver[static_cast<gsXBraidStepStatus&>(status).level()]->solveWithGuess(tstep*K.rhs() +
                                                  (M.matrix()-tstep*K.matrix())*(*u_ptr),
                                                  *ustop_ptr);
      break;
      
    case gsXBraid_typeMethod::FE_BE:
      if (static_cast<gsXBraidStepStatus&>(status).level() == 0) {
        // Forward Euler method (fine grid)
        *u_ptr = m_solver[static_cast<gsXBraidStepStatus&>(status).level()]->solveWithGuess(tstep*K.rhs() +
                                                    (M.matrix()-tstep*K.matrix())*(*u_ptr),
                                                    *ustop_ptr);
      } else {
        // Backward Euler method (coarse grids)
        *u_ptr = m_solver[static_cast<gsXBraidStepStatus&>(status).level()]->solveWithGuess(tstep*K.rhs() +
                                                    (M.matrix())*(*u_ptr),
                                                    *ustop_ptr);
      }
      break;

    case gsXBraid_typeMethod::BE_BE: {
      // Backward Euler method (all grids)
       *u_ptr = m_solver[static_cast<gsXBraidStepStatus&>(status).level()]->solveWithGuess(tstep*K.rhs() + (M.matrix())*(*u_ptr), *ustop_ptr);
      } break;

    case gsXBraid_typeMethod::CN_CN:
      // Crank-Nicholson method (all grids)
      *u_ptr = m_solver[static_cast<gsXBraidStepStatus&>(status).level()]->solveWithGuess(tstep*K.rhs() +
                                                    (M.matrix()-tstep*0.5*K.matrix())*(*u_ptr),
                                                    *ustop_ptr);
      break;
      
    case gsXBraid_typeMethod::CN_BE:
      if (static_cast<gsXBraidStepStatus&>(status).level() == 0) {
        *u_ptr = m_solver[static_cast<gsXBraidStepStatus&>(status).level()]->solveWithGuess(tstep*K.rhs() +
                                                    (M.matrix()-tstep*0.5*K.matrix())*(*u_ptr),
                                                    *ustop_ptr);
      } else {
        // Backward Euler method (coarse grids)
        *u_ptr = m_solver[static_cast<gsXBraidStepStatus&>(status).level()]->solveWithGuess(tstep*K.rhs() +
                                                    (M.matrix())*(*u_ptr),
                                                    *ustop_ptr);
      }
      break;

    default:
      throw std::runtime_error("Unsupported time-stepping method");
    }
      
    // Carry out adaptive refinement in time
    if (static_cast<gsXBraidStepStatus&>(status).level() == 0) {
      braid_Real error = static_cast<gsXBraidStepStatus&>(status).error();
      if (error != braid_Real(-1.0)) {
        braid_Int rfactor = (braid_Int) std::ceil( std::sqrt( error / 1e-3) );
        status.SetRFactor(rfactor);
      } else
        status.SetRFactor(1);
    }
    
    return braid_Int(0);
  }

  /// Sets the size of the MPI communication buffer
  braid_Int BufSize(braid_Int         *size_ptr,
                    BraidBufferStatus &status)
#if __cplusplus >= 201103L || _MSC_VER >= 1600
    override
#endif
  {
    *size_ptr = sizeof(T)*(M.numDofs()+2);
    return braid_Int(0);
  }

  /// Handles access for input/output
  braid_Int Access(braid_Vector       u,
                   BraidAccessStatus &status)
#if __cplusplus >= 201103L || _MSC_VER >= 1600
    override
#endif
  {
    if (static_cast<gsXBraidAccessStatus&>(status).done() &&
        static_cast<gsXBraidAccessStatus&>(status).timeIndex() ==
        static_cast<gsXBraidAccessStatus&>(status).times()) {
      gsMatrix<T>* u_ptr = (gsMatrix<T>*) u;
      gsInfo << "norm of the solution = " << u_ptr->norm() << std::endl;    
    }
    return braid_Int(0);
  }

  /// Performs spatial coarsening
  braid_Int Coarsen(braid_Vector           fu,
                    braid_Vector          *cu_ptr,
                    BraidCoarsenRefStatus &status)
#if __cplusplus >= 201103L || _MSC_VER >= 1600
    override
#endif
  {
    // gsInfo << "Coarsen on level = "
    //        << static_cast<gsXBraidCoarsenRefStatus&>(status).level()
    //        << " of "
    //        << static_cast<gsXBraidCoarsenRefStatus&>(status).levels()
    //        << "\n";
    gsMatrix<T> *fu_ptr = (gsMatrix<T>*) fu;    
    gsMatrix<T>* cu     = new gsMatrix<T>();
    *cu = *fu_ptr;
    *cu_ptr = (braid_Vector) cu;
    return braid_Int(0);
  }
  
  // Performs spatial refinement
  braid_Int Refine(braid_Vector           cu,
                   braid_Vector          *fu_ptr,
                   BraidCoarsenRefStatus &status)
#if __cplusplus >= 201103L || _MSC_VER >= 1600
    override
#endif
  {
    // gsInfo << "Refine on level = "
    //        << static_cast<gsXBraidCoarsenRefStatus&>(status).level()
    //        << " of "
    //        << static_cast<gsXBraidCoarsenRefStatus&>(status).levels()
    //        << "\n";
    gsMatrix<T> *cu_ptr = (gsMatrix<T>*) cu;    
    gsMatrix<T>* fu     = new gsMatrix<T>();
    *fu = *cu_ptr;
    *fu_ptr = (braid_Vector) fu;
    return braid_Int(0);
  }
};
  
} // ending namespace gismo

#endif

int main(int argc, char**argv)
{
#ifdef GISMO_WITH_XBRAID
  
  // Initialize the MPI environment and obtain the world communicator
  gsMpiComm comm = gsMpi::init(argc, argv).worldComm();

  // Set up app structure
  gsXBraid_app<real_t> app = gsXBraid_app<real_t>::create(comm, argc, argv);

  // Perform parallel-in-time multigrid
  app.solve();

#else

  gsInfo << "\n";
 
#endif

  return 0;
  
}
