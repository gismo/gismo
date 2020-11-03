/** @file gsXBraid.h

    @brief Provides declarations of the XBraid wrapper

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): M. Moller
*/

#include <braid.hpp>

namespace gismo {
  
  /**
     \brief Class defining the XBraid wrapper
  */

  class gsXBraid : public BraidApp
  {
  public:
    
    // Default constructor
    gsXBraid() = delete;

    // Constructor
    gsXBraid(MPI_Comm comm,
             int      rank,
             double   start,
             double   stop,
             int      timesteps);
       
    // Destructor
    virtual ~gsXBraid();

    // Define all the Braid Wrapper routines
    virtual int Step(braid_Vector    u_,
                     braid_Vector    ustop_,
                     braid_Vector    fstop_,
                     BraidStepStatus &pstatus) = 0;
    
    virtual int Clone(braid_Vector  u_,
                      braid_Vector *v_ptr) = 0 ;
    
    virtual int Init(double        t,
                     braid_Vector *u_ptr) = 0;
    
    virtual int Free(braid_Vector u_) = 0;
    
    virtual int Sum(double       alpha,
                    braid_Vector x_,
                    double       beta,
                    braid_Vector y_) = 0;

    virtual int SpatialNorm(braid_Vector  u_,
                            double       *norm_ptr) = 0;

    virtual int BufSize(int *size_ptr,
                        BraidBufferStatus  &status) = 0;

    virtual int BufPack(braid_Vector  u_,
                        void         *buffer,
                        BraidBufferStatus  &status) = 0;

    virtual int BufUnpack(void         *buffer,
                          braid_Vector *u_ptr,
                          BraidBufferStatus  &status) = 0;

    virtual int Access(braid_Vector       u_,
                       BraidAccessStatus &astatus) = 0;

    // Not needed in this example
    virtual int Residual(braid_Vector     u_,
                         braid_Vector     r_,
                         BraidStepStatus &pstatus) = 0;

    // Not needed in this example
    virtual int Coarsen(braid_Vector   fu_,
                        braid_Vector  *cu_ptr,
                        BraidCoarsenRefStatus &status) = 0;

    // Not needed in this example
    virtual int Refine(braid_Vector   cu_,
                       braid_Vector  *fu_ptr,
                       BraidCoarsenRefStatus &status) = 0;
    
  protected:
    int rank;
    
  };
  
}// namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsXBraid.hpp)
#endif
