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

using namespace gismo;

#ifdef GISMO_WITH_XBRAID

namespace gismo {

  enum class gsXBraid_typeMethod
  {
    FE_FE = 1, // forward Euler   (all grids)
    BE_BE = 2, // backward Euler  (all grids)
    CN_CN = 3, // Crank-Nicholson (all grids)
    FE_BE = 4, // forward Euler   (fine grid), backward Euler (coarser grids)
    CN_BE = 5  // Crank-Nicholson (fine grid), backward Euler (coarser grids)      
  };

/**
   \brief Derived class implementing the XBraid wrapper for the heat equation
*/
template<typename T>
class gsXBraid_app : public gsXBraid< gsVector<T> >
{
private:
  // Spatial discretisation parameters
  index_t numRefine, numElevate, numIncrease;

  // Temporal discretisation parameters
  index_t numSteps, typeMethod;
  T tstart, tstop, tstep;

  // Spatial discretization
  gsMultiPatch<T> mp;
  gsMultiBasis<T> bases;
  
  // Boundary conditions
  gsBoundaryConditions<T> bc;

  // Assembler options
  gsOptionList Aopt;
  
  // Expression assembler
  gsExprAssembler<T> K, M;
  gsFunctionExpr<T> f, u0, ms;
  
  // Solution
  gsVector<T> sol;

  // Solver
  typedef typename gsSparseSolver<>::CGDiagonal solver;
  solver* m_solver;
  
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
    : gsXBraid< gsVector<T> >::gsXBraid(comm, tstart, tstop, (int)numSteps),
      numRefine(numRefine),
      numElevate(numElevate),
      numIncrease(numIncrease),
      numSteps(numSteps),
      typeMethod(typeMethod),
      tstart(tstart),
      tstop(tstop),
      tstep( (tstop-tstart)/numSteps ),
      K(1,1), M(1,1),
      m_solver(new solver)
  {
    /////////////////////////////////////////////////////////////////////////////////////////////
    //                           Code for heat equation starts here                            //
    /////////////////////////////////////////////////////////////////////////////////////////////

    gsFileData<T> fd(fn);
    if (this->id() == 0) gsInfo << "Loaded file " << fd.lastPath() << "\n";

    fd.getId(0, mp); // id=0: Multipatch domain
    bases = gsMultiBasis<T>(mp);
    
    fd.getId(1, f); // id=1: right-hand side function
    if (this->id() == 0) gsInfo << "Source function " << f << "\n";
    
    fd.getId(2, bc); // id=2: boundary conditions
    if (this->id() == 0) gsInfo << "Boundary conditions:\n" << bc << "\n";

    fd.getId(3, u0); // id=3: initial conditions
    if (this->id() == 0) gsInfo << "Initial conditions:\n" << u0 << "\n";
    
    fd.getId(4, ms); // id=4: manufactured solution
    if (this->id() == 0) gsInfo << "Manufactured solution:\n" << ms << "\n";
    
    fd.getId(5, Aopt); // id=5: assembler options
    K.setOptions(Aopt);
    M.setOptions(Aopt);

    // Elevate and p-refine the basis to order k + numElevate
    // where k is the highest degree in the bases
    if ( numElevate > -1 )
    {
        // Find maximum degree with respect to all the variables
        int tmp = bases.maxDegree(0);
        for (short_t j = 1; j < mp.parDim(); ++j )
            if ( tmp < bases.maxDegree(j) )
                tmp = bases.maxDegree(j);

        // Elevate all degrees uniformly
        tmp += numElevate;
        bases.setDegree(tmp);
    }

    // Increase and p-refine the basis
    if (numIncrease > 0)
      bases.degreeIncrease(numIncrease);
      
    // h-refine the basis
    for (int i = 0; i < numRefine; ++i)
        bases.uniformRefine();
    
    // Set the basis
    K.setIntegrationElements(bases);
    M.setIntegrationElements(bases);   

    // Set the geometry map
    geometryMap G_K = K.getMap(mp);
    geometryMap G_M = M.getMap(mp);

    // Set the discretization space
    space u_K = K.getSpace(bases);
    space u_M = M.getSpace(bases);
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
          sol = m_solver->compute(M.matrix()
                                  ).solve(tstep*K.rhs() +
                                          (M.matrix()-tstep*K.matrix())*sol);
        break;
        
      case gsXBraid_typeMethod::BE_BE:
        // Backward Euler method
        
        for ( int i = 1; i<=numSteps; ++i) // for all timesteps
          // Compute the system for the timestep i (rhs is assumed constant wrt time)
          sol = m_solver->compute(M.matrix() +
                                  tstep*K.matrix()
                                  ).solve(tstep*K.rhs() +
                                          (M.matrix())*sol);        
        break;
        
      case gsXBraid_typeMethod::CN_CN:
      case gsXBraid_typeMethod::CN_BE:
        // Crank-Nicholson method
        
        for ( int i = 1; i<=numSteps; ++i) // for all timesteps
          // Compute the system for the timestep i (rhs is assumed constant wrt time)
          sol = m_solver->compute(M.matrix() +
                                  tstep*0.5*K.matrix()
                                  ).solve(tstep*K.rhs() +
                                          (M.matrix()-tstep*0.5*K.matrix())*sol);        
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
    delete m_solver;
  }
  
  /// Creates instance from command line argument
  static inline gsXBraid_app create(const gsMpiComm& comm,
                                    int              argc,
                                    char**           argv)
  {
    // Problem parameters
    std::string fn(XBRAID_DATA_DIR"pde/heat2d_square_ibvp.xml");
    
    // Spatial discretisation parameters
    index_t numRefine     = 2;
    index_t numElevate    = 0;
    index_t numIncrease   = 0;
    
    // Temporal discretisation parameters
    index_t numSteps      = 40;
    index_t typeMethod    = (index_t)gsXBraid_typeMethod::CN_CN;
    T       tfinal        = 0.1;
    
    // Parallel-in-time multigrid parameters
    index_t CFactor       = 2;
    index_t access        = 1;
    index_t maxIter       = 100;
    index_t maxLevel      = 30;
    index_t minCLevel     = 2;
    index_t numFMG        = 1;
    index_t numFMGVcyc    = 1;
    index_t numMaxRef     = 1;
    index_t numRelax      = 1;
    index_t numStorage    =-1;
    index_t print         = 2;
    index_t tnorm         = 2; // 1-norm, 2-norm, inf-norm
    
    T       absTol        = 1e-10;
    T       relTol        = 1e-3;

    bool    fmg           = false;
    bool    incrMaxLevels = false;
    bool    periodic      = false;
    bool    refine        = false;
    bool    sequential    = false;
    bool    skip          = true;
    
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

    // Parallel-in-time multigrid parameters
    cmd.addInt( "",   "numStorage", "Number of storage of the parallel-in-time multigrid solver", numStorage );
    cmd.addInt( "A",  "access", "Access level (neve [=0], =after finished [=1(default)], each iteration [=2]", access );
    cmd.addInt( "C",  "CFactor", "Coarsening factor of the parallel-in-time multigrid solver", CFactor );
    cmd.addInt( "F",  "numFMG", "Number of full multigrid steps of the parallel-in-time multigrid solver", numFMG );
    cmd.addInt( "L",  "maxLevel", "Maximum numbers of parallel-in-time multigrid levels", maxLevel );
    cmd.addInt( "M",  "maxIter", "Maximum iteration numbers  of the parallel-in-time multigrid solver", maxIter );
    cmd.addInt( "N",  "norm", "Temporal norm of the parallel-in-time multigrid solver (1-norm [=1], 2-norm [=2(default)], inf-norm [=3])", tnorm );
    cmd.addInt( "P",  "print", "Print level (no output [=0], =runtime inforation [=1], run statistics [=2(default)], debug [=3])", print );
    cmd.addInt( "R",  "numMaxRef", "Maximum number of refinements of the parallel-in-time multigrid solver", numMaxRef );
    cmd.addInt( "V",  "numFMGVcyc", "Number of full multigrid V-cycles of the parallel-in-time multigrid solver", numFMGVcyc );
    cmd.addInt( "X",  "numRelax", "Number of relaxation steps of the parallel-in-time multigrid solver", numRelax );
    cmd.addInt( "l",  "minCLevel", "Minimum level of the parallel-in-time multigrid solver", minCLevel );
    
    cmd.addReal( "",  "absTol", "Absolute tolerance of the parallel-in-time multigrid solver", absTol );
    cmd.addReal( "",  "relTol", "Relative tolerance of the parallel-in-time multigrid solver", relTol );

    cmd.addSwitch( "fmg" , "Perform full multigrid (default is off)", fmg);
    cmd.addSwitch( "incrMaxLevels" , "Increase the maximum number of parallel-in-time multigrid levels after performing a refinement (default is off)", incrMaxLevels);
    cmd.addSwitch( "periodic" , "Periodic time grid (default is off)", periodic);
    cmd.addSwitch( "refine" , "Perform refinement in time (default off)", refine);
    cmd.addSwitch( "sequential", "Set the initial guess of the parallel-in-time multigrid solver as the sequential time stepping solution (default is off)", sequential);
    cmd.addSwitch( "skip" , "Skip all work on the first down cycle of the parallel-in-time multigrid solver (default on)", skip);
    
    cmd.getValues(argc,argv);

    // Create instance
    gsXBraid_app<T> app(comm, 0.0, tfinal, typeMethod, numSteps, numRefine, numElevate, numIncrease, fn);

    if (absTol != 1e-10)
      app.SetAbsTol(absTol);
    else if (relTol != 1e-3)
      app.SetRelTol(relTol);
    else
      app.SetAbsTol(absTol);

    app.SetCFactor(CFactor);
    app.SetMaxIter(maxIter);
    app.SetMaxLevels(maxLevel);
    app.SetMaxRefinements(numMaxRef);
    app.SetMinCoarse(minCLevel);
    app.SetNFMG(numFMG);
    app.SetNFMGVcyc(numFMGVcyc);
    app.SetNRelax(numRelax);
    app.SetAccessLevel(access);
    app.SetPrintLevel(print);
    app.SetStorage(numStorage);
    app.SetTemporalNorm(tnorm);

    if (fmg)           app.SetFMG();
    if (incrMaxLevels) app.SetIncrMaxLevels();
    if (periodic)      app.SetPeriodic(1); else app.SetPeriodic(0);
    if (refine)        app.SetRefine(1);   else app.SetRefine(0);
    if (sequential)    app.SetSeqSoln(1);  else app.SetSeqSoln(0);
    if (skip)          app.SetSkip(1);     else app.SetSkip(0);

    //app.SetSpatialCoarsenAndRefine();
    
    return app;
  }

  /// Initializes a vector
  braid_Int Init(braid_Real    t,
                 braid_Vector *u_ptr)
#if __cplusplus >= 201103L || _MSC_VER >= 1600
    override
#endif
  {
    gsVector<T>* u = new gsVector<T>(M.numDofs());
    
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
                 BraidStepStatus &pstatus)
#if __cplusplus >= 201103L || _MSC_VER >= 1600
    override
#endif
  {
    gsVector<T>* u_ptr = (gsVector<T>*) u;
    gsVector<T>* ustop_ptr = (gsVector<T>*) ustop;

    // XBraid forcing
    if (fstop != NULL) {
      gsVector<T>* fstop_ptr = (gsVector<T>*) fstop;
      *u_ptr += *fstop_ptr;
    }
    
    // Get time step information
    std::pair<braid_Real, braid_Real> time =
      static_cast<gsXBraidStepStatus&>(pstatus).timeInterval();
    T tstep(time.second - time.first);
    
    switch((gsXBraid_typeMethod)typeMethod) {
    case gsXBraid_typeMethod::FE_FE:
      // Forward Euler method (all grids)
      *u_ptr = m_solver->compute(M.matrix()
                                 ).solveWithGuess(tstep*K.rhs() +
                                                  (M.matrix()-tstep*K.matrix())*(*u_ptr),
                                                  *ustop_ptr);
      break;
      
    case gsXBraid_typeMethod::FE_BE:
      if (static_cast<gsXBraidStepStatus&>(pstatus).level() == 0) {
        // Forward Euler method (fine grid)
        *u_ptr = m_solver->compute(M.matrix()
                                   ).solveWithGuess(tstep*K.rhs() +
                                                    (M.matrix()-tstep*K.matrix())*(*u_ptr),
                                                    *ustop_ptr);
      } else {
        // Backward Euler method (coarse grids)
        *u_ptr = m_solver->compute(M.matrix() +
                                   tstep*K.matrix()
                                   ).solveWithGuess(tstep*K.rhs() +
                                                    (M.matrix())*(*u_ptr),
                                                    *ustop_ptr);
      }
      break;

    case gsXBraid_typeMethod::BE_BE:
      // Backward Euler method (all grids)
      *u_ptr = m_solver->compute(M.matrix() +
                                 tstep*K.matrix()
                                 ).solveWithGuess(tstep*K.rhs() +
                                                  (M.matrix())*(*u_ptr),
                                                  *ustop_ptr);
      break;

    case gsXBraid_typeMethod::CN_CN:
      // Crank-Nicholson method (all grids)
      *u_ptr = m_solver->compute(M.matrix() +
                                 tstep*0.5*K.matrix()
                                 ).solveWithGuess(tstep*K.rhs() +
                                                  (M.matrix()-tstep*0.5*K.matrix())*(*u_ptr),
                                                  *ustop_ptr);
      break;
      
    case gsXBraid_typeMethod::CN_BE:
      if (static_cast<gsXBraidStepStatus&>(pstatus).level() == 0) {
        *u_ptr = m_solver->compute(M.matrix() +
                                   tstep*0.5*K.matrix()
                                   ).solveWithGuess(tstep*K.rhs() +
                                                    (M.matrix()-tstep*0.5*K.matrix())*(*u_ptr),
                                                    *ustop_ptr);
      } else {
        // Backward Euler method (coarse grids)
        *u_ptr = m_solver->compute(M.matrix() +
                                   tstep*K.matrix()
                                   ).solveWithGuess(tstep*K.rhs() +
                                                    (M.matrix())*(*u_ptr),
                                                    *ustop_ptr);
      }
      break;

    default:
      throw std::runtime_error("Unsupported time-stepping method");
    }
      
    // Carry out adaptive refinement in time
    if (static_cast<gsXBraidStepStatus&>(pstatus).level() == 0) {
      braid_Real error = static_cast<gsXBraidStepStatus&>(pstatus).error();
      if (error != braid_Real(-1.0)) {
        braid_Int rfactor = (braid_Int) std::ceil( std::sqrt( error / 1e-3) );
        pstatus.SetRFactor(rfactor);
      } else
        pstatus.SetRFactor(1);
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
                   BraidAccessStatus &astatus)
#if __cplusplus >= 201103L || _MSC_VER >= 1600
    override
#endif
  {
    if(static_cast<gsXBraidAccessStatus&>(astatus).done() &&
       static_cast<gsXBraidAccessStatus&>(astatus).timeIndex() ==
       static_cast<gsXBraidAccessStatus&>(astatus).times()) {
      gsVector<T>* u_ptr = (gsVector<T>*) u;
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
    gsInfo << "Coarsen\n";
    gsVector<T> *fu_ptr = (gsVector<T>*) fu;    
    gsVector<T>* cu     = new gsVector<T>();
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
    gsInfo << "Refine\n";
    gsVector<T> *cu_ptr = (gsVector<T>*) cu;    
    gsVector<T>* fu     = new gsVector<T>();
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
