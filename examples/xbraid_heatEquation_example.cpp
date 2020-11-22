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
class gsXBraid_app : public gsXBraid<T>
{
 public:
  /// Inherit all constructors from base class
  using gsXBraid<T>::gsXBraid;

  /// Creates instance from command line argument
  static inline gsXBraid_app create(const gsMpiComm& comm,
                                    int              argc,
                                    char**           argv)
  {
    // Problem parameters
    std::string fn("pde/poisson2d_bvp.xml");
    
    // Spatial discretisation parameters
    index_t numRefine  = 5;
    index_t numElevate = 0;

    // Temporal discretisation parameters
    index_t numTime    = 1;
    T       tfinal     = 1.0;
    
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
    gsXBraid_app<T> app(comm, 0.0, tfinal, numTime);

    app.SetAbsTol(absTol);
    app.SetRelTol(relTol);

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
  
  /// Destructor
  ~gsXBraid_app() override
    {}

  int Step(braid_Vector    u,
           braid_Vector    ustop,
           braid_Vector    fstop,
           BraidStepStatus &pstatus) override
  {}
  
  int Clone(braid_Vector  u,
            braid_Vector *v_ptr) override
  {}
  
  int Init(T             t,
           braid_Vector *u_ptr) override
  {}
  
  int Free(braid_Vector u) override
  {}
  
  int Sum(T            alpha,
          braid_Vector x,
          T            beta,
          braid_Vector y) override
  {}
  
  int SpatialNorm(braid_Vector  u,
                  T            *norm_ptr) override
  {}
  
  int BufSize(index_t           *size_ptr,
              BraidBufferStatus &status) override
  {}
  
  int BufPack(braid_Vector       u,
              void              *buffer,
              BraidBufferStatus &status) override
  {}
  
  int BufUnpack(void              *buffer,
                braid_Vector      *u_ptr,
                BraidBufferStatus &status) override
  {}
  
  int Access(braid_Vector       u,
             BraidAccessStatus &astatus) override
  {}
  
  // Not needed in this example
  int Residual(braid_Vector     u,
               braid_Vector     r,
               BraidStepStatus &pstatus) override
  {}
  
  // Not needed in this example
  int Coarsen(braid_Vector           fu,
              braid_Vector          *cu_ptr,
              BraidCoarsenRefStatus &status) override
  {}
  
  // Not needed in this example
  int Refine(braid_Vector           cu,
             braid_Vector          *fu_ptr,
             BraidCoarsenRefStatus &status) override
  {}
};

} // ending namespace gismo

#endif

int main(int argc, char**argv)
{
  // Initialize the MPI environment and obtain the world communicator
  gsMpiComm comm = gsMpi::init(argc, argv).worldComm();
  
#ifdef GISMO_WITH_XBRAID

  // Set up app structure
  gsXBraid_app<real_t> app = gsXBraid_app<real_t>::create(comm, argc, argv);

  // Perform parallel-in-time multigrid
  app.solve();
  
#endif

  return 0;
  
}
