/** @file gsXBraid.h

    @brief Provides declarations of the XBraid wrapper

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): M. Moller
*/

#pragma once

#include <gsCore/gsConfig.h>
#include <gsMpi/gsMpi.h>

#if !defined(GISMO_WITH_MPI)
#define braid_SEQUENTIAL 1
#endif

#include <braid.hpp>

namespace gismo {
  
  /**
     \brief Class defining the XBraid wrapper
  */

  template<typename T>
  class gsXBraid : public BraidApp
  {
  public:
    /// Constructor
    gsXBraid(const gsMpiComm& comm,
             const T&         tstart,
             const T&         tstop,
             int              ntime);
       
    /// Destructor
    virtual ~gsXBraid();

    // Performs one time step
    virtual int Step(braid_Vector    u,
                     braid_Vector    ustop,
                     braid_Vector    fstop,
                     BraidStepStatus &pstatus) = 0;

    // Clones the given vectors
    virtual int Clone(braid_Vector  u,
                      braid_Vector *v_ptr) = 0 ;

    // Initializes the given vector
    virtual int Init(T             t,
                     braid_Vector *u_ptr) = 0;

    // Fianlizes the given vector
    virtual int Free(braid_Vector u) = 0;

    // Computes the weighted sum of two given vectors
    virtual int Sum(T            alpha,
                    braid_Vector x,
                    T            beta,
                    braid_Vector y) = 0;

    // Computes the spatial norm of the given vector
    virtual int SpatialNorm(braid_Vector  u,
                            T            *norm_ptr) = 0;

    // Computes the buffer size
    virtual int BufSize(index_t           *size_ptr,
                        BraidBufferStatus &status) = 0;

    // Packes the given vector into the given buffer
    virtual int BufPack(braid_Vector       u,
                        void              *buffer,
                        BraidBufferStatus &status) = 0;

    // Unpacks the given buffer into the given vector
    virtual int BufUnpack(void              *buffer,
                          braid_Vector      *u_ptr,
                          BraidBufferStatus &status) = 0;

    // Accesses the given vector
    virtual int Access(braid_Vector       u,
                       BraidAccessStatus &astatus) = 0;

    // Calculates the residual
    virtual int Residual(braid_Vector     u,
                         braid_Vector     r,
                         BraidStepStatus &pstatus) = 0;

    /// Performs coarsening in time
    virtual int Coarsen(braid_Vector           fu,
                        braid_Vector          *cu_ptr,
                        BraidCoarsenRefStatus &status) = 0;

    /// Performs refinement in time
    virtual int Refine(braid_Vector           cu,
                       braid_Vector          *fu_ptr,
                       BraidCoarsenRefStatus &status) = 0;

    /// Runs the parallel-in-time multigrid solver
    void solve() { core.Drive(); }
    
  public:
    // Sets the maximum number of multigrid levels.
    void SetMaxLevels(int max_levels) { core.SetMaxLevels(max_levels); }
    
    // Increases the max number of multigrid levels after performing a refinement.
    void SetIncrMaxLevels() { core.SetIncrMaxLevels(); }

    // Sets whether to skip all work on the first down cycle (skip = 1).  On by default.
    void SetSkip(int skip) { core.SetSkip(skip); }

    // Sets the minimum allowed coarse grid size. gsXBraid stops
    // coarsening whenever creating the next coarser grid will result
    // in a grid smaller than min_coarse. The maximum possible coarse
    // grid size will be min_coarse*coarsening_factor.
    void SetMinCoarse(int min_coarse) { core.SetMinCoarse(min_coarse); }

    // Sets the number of relaxation sweeps *nrelax* on grid
    // *level*. Level 0 is the finest grid. One sweep is a CF
    // relaxation sweep.
    void SetNRelax(int level, int nrelax) { core.SetNRelax(level, nrelax); }

    // Sets the number of relaxation sweeps *nrelax* on all grid
    // levels. One sweep is a CF relaxation sweep.
    void SetNRelax(int nrelax) { core.SetNRelax(-1, nrelax); }
        
    // Sets absolute stopping tolerance.
    void SetAbsTol(T tol) { core.SetAbsTol(tol); }

    // Sets relative stopping tolerance.
    void SetRelTol(T tol) { core.SetRelTol(tol); }

    // Sets the temporal norm: 1-norm (1), 2-norm (2:default), inf-norm (3)
    void SetTemporalNorm(int tnorm) { core.SetTemporalNorm(tnorm); }

    // Sets the coarsening factor *cfactor* on grid *level* (default is 2)
    void SetCFactor(int level, int cfactor) { core.SetCFactor(level, cfactor); }

    // Sets the coarsening factor *cfactor* on all grid levels
    void SetCFactor( int cfactor) { core.SetCFactor(-1, cfactor); }

    // Sets periodic time grid (default is 0)
    void SetPeriodic(int periodic) { core.SetPeriodic(periodic); }

    // Sets max number of multigrid iterations.
    void SetMaxIter(int max_iter) { core.SetMaxIter(max_iter); }

    // Sets the print level for runtime print message.
    // - Level 0: no output
    // - Level 1: print runtime information like the residual history 
    // - Level 2: level 1 output, plus post-Braid run statistics (default)
    // - Level 3: level 2 output, plus debug level output.
    void SetPrintLevel(int print_level) { core.SetPrintLevel(print_level); }

    // Sets the output file for runtime print message.
    void SetPrintFile(const char *printfile_name) { core.SetPrintFile(printfile_name); }
    
    // Sets the initial guess to gsXBraid as the sequential time stepping solution.
    // - 0: The user's Init() function initializes the state vector (default)
    // - 1: Sequential time stepping, with the user's initial condition from
    //      Init(t=0) initializes the state vector
    void SetSeqSoln(int use_seq_soln) { core.SetSeqSoln(use_seq_soln); }

    // Sets the acces level for gsXBraid. This controls how often the
    // user's access routine is called.
    // - Level 0:  Never call the user's access routine
    // - Level 1:  Only call the user's access routine after gsXBraid is finished (default)
    // - Level 2:  Call the user's access routine every iteration and on every level.
    //             This is during _braid_FRestrict, during the down-cycle part of a
    //             gsXBraid iteration. 
    void SetAccessLevel(int access_level) { core.SetAccessLevel(access_level); }

    // Sets FMG (F-cycle)
    void SetFMG() { core.SetFMG(); }

    // Sets the number of initial F-cycles to do before switching to V-cycles
    void SetNFMG(int k) { core.SetNFMG(k); }

    // Sets the number of V-cycles to do at each FMG level (default is 1)
    void SetNFMGVcyc(int nfmg_Vcyc) { core.SetNFMGVcyc(nfmg_Vcyc); }

    // Sets the storage properties of the code.
    //  -1     : Default, store only C-points
    //   0     : Full storage of C- and F-Points on all levels
    //   x > 0 : Full storage on all levels >= x 
    void SetStorage(int storage) { core.SetStorage(storage); }

    // Turns time refinement on (refine = 1) or off (refine = 0).
    void SetRefine(int refine) {core.SetRefine(refine);}

    // Sets the max number of time grid refinement levels allowed.
    void SetMaxRefinements(int max_refinements) {core.SetMaxRefinements(max_refinements);}

    // Turns on built-in Richardson-based error estimation and/or
    // extrapolation with gsXBraid.  When enabled, the Richardson
    // extrapolation (RE) option (richardson == 1) is used to improve
    // the accuracy of the solution at the C-points on the finest
    // level.  When the built-in error estimate option is turned on
    // (est_error == 1), RE is used to estimate the local truncation
    // error at each point. These estimates can be accessed through
    // StepStatus and AccessStatus functions.  The last parameter is
    // local_order, which represents the LOCAL order of the* time
    // integration scheme. e.g. local_order = 2 for Backward Euler.
    // Also, the Richardson error estimate is only available after
    // roughly 1 Braid iteration.  The estimate is given a dummy value
    // of -1.0, until an actual estimate is available.  Thus after an
    // adaptive refinement, and a new hierarchy is formed, another
    // iteration must pass before the error estimates are available
    // again.
    void SetRichardsonEstimation(int est_error, int richardson, int local_order) { core.SetRichardsonEstimation(est_error, richardson, local_order); }

  public:
    // Sets user-defined residual routine.
    void SetResidual() { core.SetResidual(); }
    
    // Sets user-defined coarsening and refinement routine.
    void SetSpatialCoarsenAndRefine() { core.SetSpatialCoarsenAndRefine(); }
    
    // Sets user-defined sync routine.
    void SetSync() { core.SetSync(); }
    
  public:
    void GetNumIter(int *niter_ptr) { core.GetNumIter(niter_ptr); }
    
    void GetRNorms(int *nrequest_ptr, double *rnorms) { core.GetRNorms(nrequest_ptr, rnorms); }
    
    void GetNLevels(int *nlevels_ptr) { core.GetNLevels(nlevels_ptr); }

    int iterations() {
      int niter;
      GetNumIter(&niter);
      return niter;
    }

    int levels() {
      int nlevels;
      GetNLevels(&nlevels);
      return nlevels;
    }
    
  protected:
    /// Braid Core object
    BraidCore core;    
  };
  
}// namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsXBraid.hpp)
#endif
