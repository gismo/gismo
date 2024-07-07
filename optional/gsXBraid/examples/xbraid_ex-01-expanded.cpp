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

#ifdef gsXBraid_ENABLED

namespace gismo {

/**
   \brief Derived class implementing the XBraid wrapper for the heat equation
*/
template<typename T>
class gsXBraid_app : public gsXBraid< gsMatrix<T> >
{
private:

 public:

  gsXBraid_app(const gsMpiComm& comm,
               const T&         tstart,
               const T&         tstop,
               index_t          numSteps)
  :
  gsXBraid< gsMatrix<T> >::gsXBraid(comm, tstart, tstop, (int)numSteps)
  {}

  /// Destructor
  virtual ~gsXBraid_app()
  {

  }

  /// Initializes a vector
  braid_Int Init(braid_Real    t,
                 braid_Vector *u_ptr)
#if __cplusplus >= 201103L || _MSC_VER >= 1600
    override
#endif
  {
    gsMatrix<T> * u = new gsMatrix<T>(1, 1);

    if (t == this->tstart) /* Initial condition */
      *u<<1.0;
    else /* All other time points set to arbitrary value */
      *u<<0.456;

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
    if (fstop != NULL)
    {
      gsMatrix<T>* fstop_ptr = (gsMatrix<T>*) fstop;
      *u_ptr += *fstop_ptr;
    }

    // Get time step information
    std::pair<braid_Real, braid_Real> time =
      static_cast<gsXBraidStepStatus&>(status).timeInterval();
    T tstep(time.second - time.first);

    *u_ptr = 1./(1. + tstep)*(*u_ptr);

    status.SetRFactor(1);
    return braid_Int(0);
  }

  /// Sets the size of the MPI communication buffer
  braid_Int BufSize(braid_Int         *size_ptr,
                    BraidBufferStatus &status)
#if __cplusplus >= 201103L || _MSC_VER >= 1600
    override
#endif
  {
    *size_ptr = sizeof(T)*(1);
    return braid_Int(0);
  }

  /// Handles access for input/output
  braid_Int Access(braid_Vector       u,
                   BraidAccessStatus &status)
#if __cplusplus >= 201103L || _MSC_VER >= 1600
    override
#endif
  {
    if (static_cast<gsXBraidAccessStatus&>(status).done()// &&
        //static_cast<gsXBraidAccessStatus&>(status).timeIndex() ==
        //static_cast<gsXBraidAccessStatus&>(status).times()
        )
    {
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

  /// Performs spatial coarsening
  braid_Int Residual(braid_Vector           u,
                    braid_Vector           r,
                    BraidStepStatus       &status)
#if __cplusplus >= 201103L || _MSC_VER >= 1600
    override
#endif
  {
    // Get time step information
    std::pair<braid_Real, braid_Real> time =
      static_cast<gsXBraidStepStatus&>(status).timeInterval();
    T tstep(time.second - time.first);

    gsMatrix<T>* u_ptr = (gsMatrix<T>*) u;
    gsMatrix<T>* r_ptr = (gsMatrix<T>*) r;

    /* compute A(u_i, u_{i-1}) */
    *r_ptr = (1. + tstep)*(*u_ptr) - (*r_ptr);

    gsDebugVar(*r_ptr);

    return 0;
  }
};

} // ending namespace gismo

#endif

int main(int argc, char**argv)
{
#ifdef gsXBraid_ENABLED

  index_t           ntime         = 10;
  index_t           max_levels    = 2;
  index_t           nrelax        = 1;
  index_t           nrelax0       = -1;
  index_t           nrelaxc       = 5;
  real_t            tol           = 1.0e-06;
  index_t           cfactor       = 2;
  index_t           max_iter      = 100;
  index_t           timings       = 0;
  index_t           fmg           = 0;
  index_t           skip          = 1;
  index_t           res           = 0;
  index_t           mydt          = 0;
  index_t           sync          = 0;
  index_t           periodic      = 0;
  index_t           relax_only_cg = 0;
  index_t           bufalloc      = 0;

  gsCmdLine cmd("\nExample 1: Solve a scalar ODE \n\n");
  cmd.addInt("n","ntime","set num time points",ntime);
  cmd.addInt("l","ml","set max levels",max_levels);
  cmd.addInt("f","nu","set num F-C relaxations",nrelax);
  cmd.addInt("F","nu0","set num F-C relaxations on level 0",nrelax0);
  cmd.addInt("C","nuc","set num F-C relaxations on coarsest grid",nrelaxc);
  cmd.addReal("t","tol","set stopping tolerance",tol);
  cmd.addInt("a","cf","set coarsening factor",cfactor);
  cmd.addInt("M","mi","set max iterations",max_iter);
  cmd.addInt("s","skip","set skip relaxations on first down-cycle; 0: no skip;  1: skip",skip);
  cmd.addInt("T","timings","turn XBraid internal timings on/off",timings);
  cmd.addInt("g","fmg","use FMG cycling",fmg);
  cmd.addInt("r","res","use my residual",res);
  cmd.addInt("S","sync","enable calls to the sync function",sync);
  cmd.addInt("P","periodic","solve a periodic problem",periodic);
  cmd.addInt("B","bufalloc","user-defined MPI buffer allocation",bufalloc);
  // cmd.addInt("","tg","use user-specified time grid as global fine time grid, options are\n");
  // cmd.addInt("","                   1 - uniform time grid\n");
  // cmd.addInt("","                   2 - nonuniform time grid, where dt*0.5 for n = 1, ..., nt/2; dt*1.5 for n = nt/2+1, ..., nt\n\n");

  try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }


  // Initialize the MPI environment and obtain the world communicator
  gsMpiComm comm = gsMpi::init(argc, argv).worldComm();

  // Print MPI/OpenMP configuration
  if (comm.rank() == 0)
  {
    gsInfo << "Number of MPI processes    : " << comm.size() << std::endl;
#ifdef _OPENMP
    gsInfo << "Number of OpenMP processes : " << omp_get_num_procs() << std::endl;
#endif
  }

  // Set up app structure
  gsXBraid_app<real_t> app(comm, 0,5,10);
  app.SetMaxLevels(max_levels);
  app.SetNRelax(-1, nrelax);
  if (nrelax0 > -1)
  {
     app.SetNRelax( 0, nrelax0);
  }

  if (res)
  {
     app.SetResidual( );
  }
  app.SetNRelax(max_levels-1, nrelaxc);
  app.SetAbsTol(tol);
  app.SetCFactor(-1, cfactor);
  app.SetMaxIter(max_iter);
  // app.SetTimings(timings);
  app.SetSkip(skip);
  if (fmg)
  {
     app.SetFMG();
  }
  if (periodic)
  {
     app.SetPeriodic(periodic);
  }
  // if (relax_only_cg)
  // {
  //    app.SetRelaxOnlyCG(relax_only_cg);
  // }
  // if (bufalloc)
  // {
  //    app.SetBufAllocFree(my_BufAlloc, my_BufFree);
  // }

  // Perform parallel-in-time multigrid
  app.solve();

#else

  gsInfo << "\n";

#endif

  return 0;

}
