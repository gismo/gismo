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
    void SetMaxLevels(int max_levels) { core.SetMaxLevels(max_levels); }
    
    void SetIncrMaxLevels() { core.SetIncrMaxLevels(); }
    
    void SetSkip(int skip) { core.SetSkip(skip); }
    
    void SetMinCoarse(int min_coarse) { core.SetMinCoarse(min_coarse); }

    void SetNRelax(int level, int nrelax) { core.SetNRelax(level, nrelax); }

    void SetAbsTol(T tol) { core.SetAbsTol(tol); }
    
    void SetRelTol(T tol) { core.SetRelTol(tol); }
    
    void SetTemporalNorm(int tnorm) { core.SetTemporalNorm(tnorm); }
    
    void SetCFactor(int level, int cfactor) { core.SetCFactor(level, cfactor); }

    void SetAggCFactor(int cfactor0) { core.SetAggCFactor(cfactor0); }

    void SetSpatialCoarsenAndRefine() { core.SetSpatialCoarsenAndRefine(); }

    void SetPeriodic(int periodic) { core.SetPeriodic(periodic); }
    
    void SetSync() { core.SetSync(); }
    
    void SetResidual() { core.SetResidual(); }
    
    void SetMaxIter(int max_iter) { core.SetMaxIter(max_iter); }
    
    void SetPrintLevel(int print_level) { core.SetPrintLevel(print_level); }
    
    void SetSeqSoln(int use_seq_soln) { core.SetSeqSoln(use_seq_soln); }
    
    void SetPrintFile(const char *printfile_name) { core.SetPrintFile(printfile_name); }
    
    void SetAccessLevel(int access_level) { core.SetAccessLevel(access_level); }
    
    void SetFMG() { core.SetFMG(); }
    
    void SetNFMG(int k) { core.SetNFMG(k); }
    
    void SetNFMGVcyc(int nfmg_Vcyc) { core.SetNFMGVcyc(nfmg_Vcyc); }
    
    void SetStorage(int storage) { core.SetStorage(storage); }
    
    void SetRefine(int refine) {core.SetRefine(refine);}
    
    void SetMaxRefinements(int max_refinements) {core.SetMaxRefinements(max_refinements);}
    
    void SetRichardsonEstimation(int est_error, int richardson, int local_order) { core.SetRichardsonEstimation(est_error, richardson, local_order); }
    
    void GetNumIter(int *niter_ptr) { core.GetNumIter(niter_ptr); }
    
    void GetRNorms(int *nrequest_ptr, double *rnorms) { core.GetRNorms(nrequest_ptr, rnorms); }
    
    void GetNLevels(int *nlevels_ptr) { core.GetNLevels(nlevels_ptr); }
    
  protected:
    /// Braid Core object
    BraidCore core;    
  };
  
}// namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsXBraid.hpp)
#endif
