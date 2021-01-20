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

/**
   \brief Derived class implementing the XBraid wrapper for the heat equation
*/
template<typename T>
class gsXBraid_app : public gsXBraid< gsVector<T> >
{
private:
  // Spatial discretisation parameters
  index_t numRefine, numElevate;

  // Temporal discretisation parameters
  index_t numTime;
  T tstart, tstop, theta, tstep;

  // Spatial discretization
  gsMultiPatch<T> patches;
  gsMultiBasis<T> bases;
  
  // Boundary conditions
  gsBoundaryConditions<T> bcInfo;
  gsConstantFunction<T> g_D, g_N;

  // Expression assembler
  gsExprAssembler<T> K, M;
  gsConstantFunction<T> f;
  
  // Solution
  gsVector<T> sol;
  
  typedef typename gsExprAssembler<T>::geometryMap geometryMap;
  typedef typename gsExprAssembler<T>::variable    variable;
  typedef typename gsExprAssembler<T>::space       space;
  typedef typename gsExprAssembler<T>::solution    solution;
  
 public:
  /// Contructor
  gsXBraid_app(const gsMpiComm& comm,
               const T&         tstart,
               const T&         tstop,
               index_t          numTime,
               index_t          numRefine,
               index_t          numElevate)
    : gsXBraid< gsVector<T> >::gsXBraid(comm, tstart, tstop, (int)numTime),
      numRefine(numRefine),
      numElevate(numElevate),
      numTime(numTime),
      tstart(tstart),
      tstop(tstop),
      theta(0.0),
      tstep( (tstop-tstart)/numTime ),
      patches(*gsNurbsCreator<>::BSplineSquareDeg(2)),
      bases(patches),
      g_D(0,2), g_N(1,2),
      K(1,1), M(1,1), f(1,2)
  {
    /////////////////////////////////////////////////////////////////////////////////////////////
    //                           Code for heat equation starts here                            //
    /////////////////////////////////////////////////////////////////////////////////////////////
    
    // Define geometry, must be a gsMultiPatch object
    patches.computeTopology();

    // Boundary conditions    
    bcInfo.addCondition(0, boundary::west,  condition_type::neumann  , &g_N);
    bcInfo.addCondition(0, boundary::east,  condition_type::dirichlet, &g_D);
    bcInfo.addCondition(0, boundary::north, condition_type::dirichlet, &g_D);
    bcInfo.addCondition(0, boundary::south, condition_type::dirichlet, &g_D);

    // Elevate and p-refine the basis to order k + numElevate
    // where k is the highest degree in the bases
    if ( numElevate > -1 )
    {
        // Find maximum degree with respect to all the variables
        int tmp = bases.maxDegree(0);
        for (short_t j = 1; j < patches.parDim(); ++j )
            if ( tmp < bases.maxDegree(j) )
                tmp = bases.maxDegree(j);

        // Elevate all degrees uniformly
        tmp += numElevate;
        bases.setDegree(tmp);
    }

    // h-refine the basis
    for (int i = 0; i < numRefine; ++i)
        bases.uniformRefine();
    
    // Set the basis
    K.setIntegrationElements(bases);
    M.setIntegrationElements(bases);   

    // Set the geometry map
    geometryMap G_K = K.getMap(patches);
    geometryMap G_M = M.getMap(patches);

    // Set the discretization space
    space u_K = K.getSpace(bases);
    space u_M = M.getSpace(bases);
    u_K.setInterfaceCont(0);
    u_M.setInterfaceCont(0);
    u_K.addBc( bcInfo.get("Dirichlet") );
    u_M.addBc( bcInfo.get("Dirichlet") );

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
    K.assembleRhsBc(u_K * g_Neumann.val() * nv(G_K).norm(), bcInfo.neumannSides() );

    if (this->id() == 0) {
      gsStopwatch clock;
      clock.restart();
      
      gsSparseSolver<>::CGDiagonal solver;
      sol.setZero(M.numDofs());
      
      for ( int i = 1; i<=numTime; ++i) // for all timesteps
        // Compute the system for the timestep i (rhs is assumed constant wrt time)
        sol = solver.compute(M.matrix() +
                             tstep*theta*K.matrix()
                             ).solve(tstep*K.rhs() +
                                     (M.matrix()-tstep*(1.0-theta)*K.matrix())*sol);
      
      gsInfo << "norm of the solution = " << sol.norm() << "\n"
             << "wall time = " << clock.stop() << std::endl;
    }
  }

  /// Destructor
  virtual ~gsXBraid_app() {}
  
  /// Creates instance from command line argument
  static inline gsXBraid_app create(const gsMpiComm& comm,
                                    int              argc,
                                    char**           argv)
  {
    // Problem parameters
    std::string fn("pde/poisson2d_bvp.xml");
    
    // Spatial discretisation parameters
    index_t numRefine  = 2;
    index_t numElevate = 0;
    
    // Temporal discretisation parameters
    index_t numTime    = 40;
    T       tfinal     = 0.1;
    
    // Parallel-in-time multigrid parameters
    index_t CFactor    = 2;
    index_t info       = 2;
    index_t maxIter    = 100;
    index_t maxLevel   = 30;
    index_t minCLevel  = 2;
    index_t numFMG     = 1;
    index_t numFMGVcyc = 1;
    index_t numMaxRef  = 1;
    index_t numRelax   = 1;
    index_t numStorage =-1;
    index_t tnorm      = 2; // 1-norm, 2-norm, inf-norm
    
    T       absTol     = 1e-10;
    T       relTol     = 1e-3;

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
    cmd.addInt( "r", "uniformRefine", "Number of uniform h-refinement steps to perform before solving",  numRefine );

    // Temporal diescretisation parameters
    cmd.addInt( "n",  "timeSteps", "Number of parallel-in-time steps", numTime );
    cmd.addReal( "t", "time", "Final time", tfinal );

    // Parallel-in-time multigrid parameters
    cmd.addInt( "C",  "CFactor", "Coarsening factor of the parallel-in-time multigrid solver", CFactor );
    cmd.addInt( "I",  "info", "Print level (no output [=0], =runtime inforation [=1], run statistics [=2(default)], debug [=3])", info );
    cmd.addInt( "M",  "maxIter", "Maximum iteration numbers  of the parallel-in-time multigrid solver", maxIter );
    cmd.addInt( "L",  "maxLevel", "Maximum numbers of parallel-in-time multigrid levels", maxLevel );
    cmd.addInt( "l",  "minCLevel", "Minimum level of the parallel-in-time multigrid solver", minCLevel );
    cmd.addInt( "F",  "numFMG", "Number of full multigrid steps of the parallel-in-time multigrid solver", numFMG );
    cmd.addInt( "V",  "numFMGVcyc", "Number of full multigrid V-cycles of the parallel-in-time multigrid solver", numFMGVcyc );
    cmd.addInt( "R",  "numMaxRef", "Maximum number of refinements of the parallel-in-time multigrid solver", numMaxRef );
    cmd.addInt( "X",  "numRelax", "Number of relaxation steps of the parallel-in-time multigrid solver", numRelax );
    cmd.addInt( "",  "numStorage", "Number of storage of the parallel-in-time multigrid solver", numStorage );
    cmd.addInt( "T",  "tnorm", "Temporal norm of the parallel-in-time multigrid solver (1-norm [=1], 2-norm [=2(default)], inf-norm [=3])", tnorm );
    
    cmd.addReal( "", "absTol", "Absolute tolerance of the parallel-in-time multigrid solver", absTol );
    cmd.addReal( "", "relTol", "Relative tolerance of the parallel-in-time multigrid solver", relTol );

    cmd.addSwitch( "fmg" , "Perform full multigrid (default is off)", fmg);
    cmd.addSwitch( "incrMaxLevels" , "Increase the maximum number of parallel-in-time multigrid levels after performing a refinement (default is off)", incrMaxLevels);
    cmd.addSwitch( "periodic" , "Periodic time grid (default is off)", periodic);
    cmd.addSwitch( "refine" , "Perform refinement in time (default off)", refine);
    cmd.addSwitch( "sequential", "Set the initial guess of the parallel-in-time multigrid solver as the sequential time stepping solution (default is off)", sequential);
    cmd.addSwitch( "skip" , "Skip all work on the first down cycle of the parallel-in-time multigrid solver (default on)", skip);
    
    cmd.getValues(argc,argv);

    // Create instance
    gsXBraid_app<T> app(comm, 0.0, tfinal, numTime, numRefine, numElevate);

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
    app.SetPrintLevel(info);
    app.SetStorage(numStorage);
    app.SetTemporalNorm(tnorm);

    if (fmg) app.SetFMG();
    if (incrMaxLevels) app.SetIncrMaxLevels();
    if (periodic) app.SetPeriodic(1); else app.SetPeriodic(0);
    if (refine) app.SetRefine(1); else app.SetRefine(0);
    if (sequential) app.SetSeqSoln(1); else app.SetSeqSoln(0);
    if (skip) app.SetSkip(1); else app.SetSkip(0);
   
    return app;
  }

  /// Initializes a vector
  braid_Int Init(braid_Real    t,
                 braid_Vector *u_ptr)
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
  {
    gsVector<T>* u_ptr = (gsVector<T>*) u;
    
    // Get time step information
    std::pair<braid_Real, braid_Real> time =
      static_cast<gsXBraidStepStatus&>(pstatus).timeInterval();
    T tstep(time.second - time.first);
    
    // Solve spatial problem
    gsSparseSolver<>::CGDiagonal solver;
    *u_ptr = solver.compute(M.matrix() +
                            tstep*theta*K.matrix()
                            ).solve(tstep*K.rhs() +
                                    (M.matrix()-tstep*(1.0-theta)*K.matrix())*(*u_ptr));    
    // no refinement
    pstatus.SetRFactor(1);
    return braid_Int(0);
  }

  /// Sets the size of the MPI communication buffer
  braid_Int BufSize(braid_Int         *size_ptr,
                    BraidBufferStatus &status)
  {
    *size_ptr = sizeof(T)*(M.numDofs()+2);
    return braid_Int(0);
  }

  /// Handles access for input/output
  braid_Int Access(braid_Vector       u,
                   BraidAccessStatus &astatus)
  {
    if(static_cast<gsXBraidAccessStatus&>(astatus).done() &&
       static_cast<gsXBraidAccessStatus&>(astatus).timeIndex() ==
       static_cast<gsXBraidAccessStatus&>(astatus).times()) {
      gsVector<T>* u_ptr = (gsVector<T>*) u;
      gsInfo << "norm of the solution = " << u_ptr->norm() << std::endl;    
    }
    return braid_Int(0);
  }
  
  // Not needed in this example
  braid_Int Residual(braid_Vector     u,
                     braid_Vector     r,
                     BraidStepStatus &pstatus)
  {
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
