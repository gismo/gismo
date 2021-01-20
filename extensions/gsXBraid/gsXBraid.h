/** @file gsXBraid.h

    @brief Provides declarations of the XBraid wrapper

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): M. Moller
*/

#pragma once

#include <gismo.h>

#if !defined(GISMO_WITH_MPI)
#define braid_SEQUENTIAL 1
#endif

#include <braid.hpp>

namespace gismo {

  class gsXBraidAccessStatus;
  class gsXBraidSyncStatus;
  class gsXBraidStepStatus;
  class gsXBraidCoarsenRefStatus;
  class gsXBraidBufferStatus;
  class gsXBraidObjectiveStatus;
  
  /**
     \brief Class defining the XBraid wrapper

     The gsXBraid class wraps the BraidApp class provided by the
     XBraid project and adds a set of commodity functions.

     In order to implement an XBraid application the user has to
     implement a derived class

     \code{.cpp}
     template<typename T>
     class gsXBraid_app : public gsXBraid<T>
     { ... };
     \endcode
     
     and implement the following application-specific functions:

     \code{.cpp}
     braid_Int Access(...)
     braid_Int BufPack(...)
     braid_Int BufSize(...)
     braid_Int BufUnpack(...)
     braid_Int Clone(...)
     braid_Int Free(...)
     braid_Int Init(...)
     braid_Int Residual(...)
     braid_Int SpatialNorm(...)
     braid_Int Step(...)
     \endcode
     
     which are declared as (pure) virtual functions in BraidApp.

     The generic implementation of the gsXBraid class leaves all of
     these methods unimplemented. We also provide specializations for
     gsXBraid< gsMatrix<T> > and gsXBraid< std::vector< gsMatrix<T> >
     > which assume that the data type for storing the solution
     (passed as braid_Vector) is of type gsMatrix<T> and std::vector<
     gsMatrix<T> >, respectively. The latter can be used to pass a
     hierarchy of matrices/vectors in a multi-level setup.
  */

  template<typename T>
  class gsXBraid : public BraidApp
  {
  public:
    /// Constructor
    gsXBraid(const gsMpiComm& comm,
             const braid_Real tstart,
             const braid_Real tstop,
             braid_Int        ntime);
       
    /// Destructor
    virtual ~gsXBraid();

    /// Free
    virtual braid_Int Free(braid_Vector u) { return braid_Int(0); }
    
    /// Runs the parallel-in-time multigrid solver
    void solve() { core.Drive(); }
    
  public:
    /// Sets the maximum number of multigrid levels.
    void SetMaxLevels(braid_Int max_levels) { core.SetMaxLevels(max_levels); }
    
    /// Increases the max number of multigrid levels after performing a refinement.
    void SetIncrMaxLevels() { core.SetIncrMaxLevels(); }

    /// Sets whether to skip all work on the first down cycle (skip = 1).  On by default.
    void SetSkip(braid_Int skip) { core.SetSkip(skip); }

    /// Sets the minimum allowed coarse grid size. gsXBraid stops
    /// coarsening whenever creating the next coarser grid will result
    /// in a grid smaller than min_coarse. The maximum possible coarse
    /// grid size will be min_coarse*coarsening_factor.
    void SetMinCoarse(braid_Int min_coarse) { core.SetMinCoarse(min_coarse); }

    /// Sets the number of relaxation sweeps *nrelax* on grid
    /// *level*. Level 0 is the finest grid. One sweep is a CF
    /// relaxation sweep.
    void SetNRelax(braid_Int level, braid_Int nrelax) { core.SetNRelax(level, nrelax); }

    /// Sets the number of relaxation sweeps *nrelax* on all grid
    /// levels. One sweep is a CF relaxation sweep.
    void SetNRelax(braid_Int nrelax) { core.SetNRelax(-1, nrelax); }
        
    /// Sets absolute stopping tolerance.
    void SetAbsTol(braid_Real tol) { core.SetAbsTol(tol); }

    /// Sets relative stopping tolerance.
    void SetRelTol(braid_Real tol) { core.SetRelTol(tol); }

    /// Sets the temporal norm: 1-norm (1), 2-norm (2:default), inf-norm (3)
    void SetTemporalNorm(braid_Int tnorm) { core.SetTemporalNorm(tnorm); }

    /// Sets the coarsening factor *cfactor* on grid *level* (default is 2)
    void SetCFactor(braid_Int level, braid_Int cfactor) { core.SetCFactor(level, cfactor); }

    /// Sets the coarsening factor *cfactor* on all grid levels
    void SetCFactor(braid_Int cfactor) { core.SetCFactor(-1, cfactor); }

    /// Sets periodic time grid (default is 0)
    void SetPeriodic(braid_Int periodic) { core.SetPeriodic(periodic); }

    /// Sets max number of multigrid iterations.
    void SetMaxIter(braid_Int max_iter) { core.SetMaxIter(max_iter); }

    /// Sets the print level for runtime print message.
    /// - Level 0: no output
    /// - Level 1: print runtime information like the residual history 
    /// - Level 2: level 1 output, plus post-Braid run statistics (default)
    /// - Level 3: level 2 output, plus debug level output.
    void SetPrintLevel(braid_Int print_level) { core.SetPrintLevel(print_level); }

    /// Sets the output file for runtime print message.
    void SetPrintFile(const char *printfile_name) { core.SetPrintFile(printfile_name); }
    
    /// Sets the initial guess to gsXBraid as the sequential time stepping solution.
    /// - 0: The user's Init() function initializes the state vector (default)
    /// - 1: Sequential time stepping, with the user's initial condition from
    ///      Init(t=0) initializes the state vector
    void SetSeqSoln(braid_Int use_seq_soln) { core.SetSeqSoln(use_seq_soln); }

    /// Sets the acces level for gsXBraid. This controls how often the
    /// user's access routine is called.
    /// - Level 0:  Never call the user's access routine
    /// - Level 1:  Only call the user's access routine after gsXBraid is finished (default)
    /// - Level 2:  Call the user's access routine every iteration and on every level.
    ///             This is during _braid_FRestrict, during the down-cycle part of a
    ///             gsXBraid iteration. 
    void SetAccessLevel(braid_Int access_level) { core.SetAccessLevel(access_level); }

    /// Sets FMG (F-cycle)
    void SetFMG() { core.SetFMG(); }

    /// Sets the number of initial F-cycles to do before switching to V-cycles
    void SetNFMG(braid_Int k) { core.SetNFMG(k); }

    /// Sets the number of V-cycles to do at each FMG level (default is 1)
    void SetNFMGVcyc(braid_Int nfmg_Vcyc) { core.SetNFMGVcyc(nfmg_Vcyc); }

    /// Sets the storage properties of the code.
    ///  -1     : Default, store only C-points
    ///   0     : Full storage of C- and F-Points on all levels
    ///   x > 0 : Full storage on all levels >= x 
    void SetStorage(braid_Int storage) { core.SetStorage(storage); }

    /// Turns time refinement on (refine = 1) or off (refine = 0).
    void SetRefine(braid_Int refine) {core.SetRefine(refine);}

    /// Sets the max number of time grid refinement levels allowed.
    void SetMaxRefinements(braid_Int max_refinements) {core.SetMaxRefinements(max_refinements);}

    /// Turns on built-in Richardson-based error estimation and/or
    /// extrapolation with gsXBraid.  When enabled, the Richardson
    /// extrapolation (RE) option (richardson == 1) is used to improve
    /// the accuracy of the solution at the C-points on the finest
    /// level.  When the built-in error estimate option is turned on
    /// (est_error == 1), RE is used to estimate the local truncation
    /// error at each point. These estimates can be accessed through
    /// StepStatus and AccessStatus functions.  The last parameter is
    /// local_order, which represents the LOCAL order of the* time
    /// integration scheme. e.g. local_order = 2 for Backward Euler.
    /// Also, the Richardson error estimate is only available after
    /// roughly 1 Braid iteration.  The estimate is given a dummy value
    /// of -1.0, until an actual estimate is available.  Thus after an
    /// adaptive refinement, and a new hierarchy is formed, another
    /// iteration must pass before the error estimates are available
    /// again.
    void SetRichardsonEstimation(braid_Int est_error, braid_Int richardson, braid_Int local_order)
    { core.SetRichardsonEstimation(est_error, richardson, local_order); }

  public:
    /// Sets user-defined residual routine.
    void SetResidual() { core.SetResidual(); }
    
    /// Sets user-defined coarsening and refinement routine.
    void SetSpatialCoarsenAndRefine() { core.SetSpatialCoarsenAndRefine(); }
    
    /// Sets user-defined sync routine.
    void SetSync() { core.SetSync(); }
    
  public:
    /// Gets the number of iterations (XBraid style)
    void GetNumIter(braid_Int *niter_ptr) { core.GetNumIter(niter_ptr); }

    /// Gets the residual norm (XBraid style)
    void GetRNorms(braid_Int *nrequest_ptr, braid_Real *rnorms) { core.GetRNorms(nrequest_ptr, rnorms); }

    /// Gets the total number of levels (XBraid style)
    void GetNLevels(braid_Int *nlevels_ptr) { core.GetNLevels(nlevels_ptr); }

    /// Returns the number of iterations
    braid_Int iterations() {
      braid_Int niter;
      GetNumIter(&niter);
      return niter;
    }

    /// Returns the residual norm
    braid_Real norm(braid_Int nrequest) {
      braid_Real rnorm;
      GetRNorms(&nrequest, &rnorm);
      return rnorm;
    }
    
    /// Returns the total number of levels
    braid_Int levels() {
      braid_Int nlevels;
      GetNLevels(&nlevels);
      return nlevels;
    }
    
  protected:
    /// Braid Core object
    BraidCore core;    
  };

  
  /**
     \brief Specializations for gsXBraid< gsMatrix<T> >
  */
  template<typename T>
  class gsXBraid< gsMatrix<T> > : public gsXBraid<T>
  {
  public:
    /// Constructor
    gsXBraid(const gsMpiComm& comm,
             const braid_Real tstart,
             const braid_Real tstop,
             braid_Int        ntime);
    
    /// Destructor
    virtual ~gsXBraid();
    
    /// Clones a given vector
    virtual braid_Int Clone(braid_Vector  u,
                            braid_Vector *v_ptr)
    {
      gsMatrix<T>* _u = (gsMatrix<T>*) u;
      gsMatrix<T>*  v = new gsMatrix<T>();
      (*v) = (*_u);
      *v_ptr = (braid_Vector) v;
      return braid_Int(0);
    }
    
    /// Frees a given vector
    virtual braid_Int Free(braid_Vector u)
    {
      gsMatrix<T>* _u = (gsMatrix<T>*) u;
      delete _u;
      return braid_Int(0);
    }
    
    /// Computes the sum of two given vectors
    virtual braid_Int Sum(braid_Real   alpha,
                          braid_Vector x,
                          braid_Real   beta,
                          braid_Vector y)
    {
      gsMatrix<T>* _x = (gsMatrix<T>*) x;
      gsMatrix<T>* _y = (gsMatrix<T>*) y;
      *_y = (T)alpha * (*_x) + (T)beta * (*_y);
      return braid_Int(0);
    }
    
    /// Computes the spatial norm of a given vector
    virtual braid_Int SpatialNorm(braid_Vector  u,
                                  braid_Real   *norm_ptr)
    {
      gsMatrix<T> *_u = (gsMatrix<T>*) u;    
      *norm_ptr = _u->norm();
      return braid_Int(0);
    }
    
    /// Packs the given vector into the MPI communication buffer
    virtual braid_Int BufPack(braid_Vector       u,
                              void              *buffer,
                              BraidBufferStatus &status)
    {
      gsMatrix<T> *_u = (gsMatrix<T>*) u;
      T* _buffer      = (T*) buffer;
      T* _data        = _u->data();
      index_t size    = _u->rows()*_u->cols();
      
      _buffer[0] = _u->rows();
      _buffer[1] = _u->cols();
      for (index_t idx = 0; idx < size; ++idx)
        _buffer[idx+2] = _data[idx];
      
      status.SetSize(sizeof(T)*(size+2));
      return braid_Int(0);
    }

    /// Unpacks a vector from the MPI communication buffer
    virtual braid_Int BufUnpack(void              *buffer,
                                braid_Vector      *u_ptr,
                                BraidBufferStatus &status)
    {
      T* _buffer     = (T*) buffer;
      index_t rows   = _buffer[0];
      index_t cols   = _buffer[1];
      gsMatrix<T>* u = new gsMatrix<T>(rows,cols);
      T* _data       = u->data();

    for (index_t idx = 0; idx < rows*cols; ++idx)
      _data[idx] = _buffer[idx+2];
    
    *u_ptr = (braid_Vector) u;
    return braid_Int(0);
    }
  };


  /**
     \brief Specializations for gsXBraid< std::vector< gsMatrix<T> > >
  */
  template<typename T>
  class gsXBraid< std::vector< gsMatrix<T> > > : public gsXBraid<T>
  {
  public:
    /// Constructor
    gsXBraid(const gsMpiComm& comm,
             const braid_Real tstart,
             const braid_Real tstop,
             braid_Int        ntime);
    
    /// Destructor
    virtual ~gsXBraid();
    
    /// Clones a given vector
    virtual braid_Int Clone(braid_Vector  u,
                            braid_Vector *v_ptr)
    {
      std::vector< gsMatrix<T> >* _u = (std::vector< gsMatrix<T> >*) u;
      std::vector< gsMatrix<T> >*  v = new std::vector< gsMatrix<T> >();

      for (typename std::vector< gsMatrix<T> >::const_iterator it = _u->cbegin();
           it != _u->cend(); ++it)
        v->push_back( *it );
      *v_ptr = (braid_Vector) v;
      return braid_Int(0);
    }
    
    /// Frees a given vector
    virtual braid_Int Free(braid_Vector u)
    {
      gsMatrix<T>* _u = (gsMatrix<T>*) u;
      delete _u;
      return braid_Int(0);
    }
    
    /// Computes the sum of two given vectors
    virtual braid_Int Sum(braid_Real   alpha,
                          braid_Vector x,
                          braid_Real   beta,
                          braid_Vector y)
    {
      gsMatrix<T>* _x = (gsMatrix<T>*) x;
      gsMatrix<T>* _y = (gsMatrix<T>*) y;
      *_y = (T)alpha * (*_x) + (T)beta * (*_y);
      return braid_Int(0);
    }
    
    /// Computes the spatial norm of a given vector
    virtual braid_Int SpatialNorm(braid_Vector  u,
                                  braid_Real   *norm_ptr)
    {
      gsMatrix<T> *_u = (gsMatrix<T>*) u;    
      *norm_ptr = _u->norm();
      return braid_Int(0);
    }
    
    /// Packs the given vector into the MPI communication buffer
    virtual braid_Int BufPack(braid_Vector       u,
                              void              *buffer,
                              BraidBufferStatus &status)
    {
      gsMatrix<T> *_u = (gsMatrix<T>*) u;
      T* _buffer      = (T*) buffer;
      T* _data        = _u->data();
      index_t size    = _u->rows()*_u->cols();
      
      _buffer[0] = _u->rows();
      _buffer[1] = _u->cols();
      for (index_t idx = 0; idx < size; ++idx)
        _buffer[idx+2] = _data[idx];
      
      status.SetSize(sizeof(T)*(size+2));
      return braid_Int(0);
    }

    /// Unpacks a vector from the MPI communication buffer
    virtual braid_Int BufUnpack(void              *buffer,
                                braid_Vector      *u_ptr,
                                BraidBufferStatus &status)
    {
      T* _buffer     = (T*) buffer;
      index_t rows   = _buffer[0];
      index_t cols   = _buffer[1];
      gsMatrix<T>* u = new gsMatrix<T>(rows,cols);
      T* _data       = u->data();

    for (index_t idx = 0; idx < rows*cols; ++idx)
      _data[idx] = _buffer[idx+2];
    
    *u_ptr = (braid_Vector) u;
    return braid_Int(0);
    }
  };
  
  
  /**
     \brief Class defining the XBraid access status wrapper

     The wrapper provides all functionality of the BraidAccessStatus
     class plus some functions that return the information by value
  */
  class gsXBraidAccessStatus : public BraidAccessStatus
  {
  public:    
    /// Returns the number of iterations
    braid_Int iterations() {
      braid_Int iter;
      GetIter(&iter);
      return iter;
    }

    /// Returns the current multigrid level
    braid_Int level() {
      braid_Int level;
      GetLevel(&level);
      return level;
    }

    /// Returns the total number of multigrid levels
    braid_Int levels() {
      braid_Int nlevels;
      GetNLevels(&nlevels);
      return nlevels;
    }

    /// Returns the total number of refinements
    braid_Int refines() {
      braid_Int nref;
      GetNRefine(&nref);
      return nref;
    }

    /// Returns the total number of time instances
    braid_Int times() {
      braid_Int ntpoints;
      GetNTPoints(&ntpoints);
      return ntpoints;
    }
    
    /// Returns true if XBraid has completed
    bool done() {
      braid_Int status;
      GetDone(&status);
      return bool(status);
    }

    /// ???
    braid_Int callingFunction() {
      braid_Int callingfcn;
      GetCallingFunction(&callingfcn);
      return callingfcn;
    }

    /// Returns the current time instance
    braid_Real time() {
      braid_Real t;
      GetT(&t);
      return t;
    }
    
    /// Returns the index of the time instance
    braid_Int timeIndex() {
      braid_Int tindex;
      GetTIndex(&tindex);
      return tindex;
    }

    /// ???
    braid_Int test() {
      braid_Int wtest;
      GetWrapperTest(&wtest);
      return wtest;
    }

    /// Returns the residual norm
    braid_Real norm() {
      braid_Real rnorm;
      GetResidual(&rnorm);
      return rnorm;
    }

    /// Returns the estimated error
    braid_Real error() {
      braid_Real errorest;
      GetSingleErrorEstAccess(&errorest);
      return errorest;
    }
  };

  /**
     \brief Class defining the XBraid sync status wrapper

     The wrapper provides all functionality of the BraidSyncStatus
     class plus some functions that return the information by value
  */
  class gsXBraidSyncStatus : public BraidSyncStatus
  {
  public:
    /// Returns the number of iterations
    braid_Int iterations() {
      braid_Int iter;
      GetIter(&iter);
      return iter;
    }

    /// Returns the current multigrid level
    braid_Int level() {
      braid_Int level;
      GetLevel(&level);
      return level;
    }

    /// Returns the total number of multigrid levels
    braid_Int levels() {
      braid_Int nlevels;
      GetNLevels(&nlevels);
      return nlevels;
    }

    /// Returns the total number of refinements
    braid_Int refines() {
      braid_Int nref;
      GetNRefine(&nref);
      return nref;
    }
    
    /// Returns the total number of time instances
    braid_Int times() {
      braid_Int ntpoints;
      GetNTPoints(&ntpoints);
      return ntpoints;
    }
    
    /// Returns true if XBraid is completed
    bool done() {
      braid_Int status;
      GetDone(&status);
      return bool(status);
    }

    /// ???
    braid_Int callingFunction() {
      braid_Int callingfcn;
      GetCallingFunction(&callingfcn);
      return callingfcn;
    }

    /// Returns the estimated errors
    braid_Real errors() {
      braid_Real errorest;
      GetAllErrorEst(&errorest);
      return errorest;
    }

    /// Returns the number of estimated errors
    braid_Int nerrors() {
      braid_Int numerrorest;
      GetNumErrorEst(&numerrorest);
      return numerrorest;
    }
  };
  
}// namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsXBraid.hpp)
#endif
