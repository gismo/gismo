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
    index_t numRefine  = 5;
    index_t numElevate = 0;
    index_t numTime    = 1;
    T       tfinal     = 1.0;    
    std::string fn("pde/poisson2d_bvp.xml");
    
    gsCmdLine cmd("Tutorial on solving a Heat equation problem using parallel-in-time multigrid.");

    cmd.addInt( "e", "degreeElevation",
                "Number of degree elevation steps to perform before solving (0: equalize degree in all directions)", numElevate );
    cmd.addInt( "r", "uniformRefine", "Number of uniform h-refinement steps to perform before solving",  numRefine );
    cmd.addInt( "n", "timeSteps", "Number of parallel-in-time steps", numTime );
    cmd.addString( "f", "file", "Input XML file", fn );
    cmd.addReal( "t", "time", "Final time", tfinal );

    cmd.getValues(argc,argv);

    return gsXBraid_app<T>(comm, 0.0, tfinal, numTime);
  }
  
  /// Destructor
  ~gsXBraid_app()
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
