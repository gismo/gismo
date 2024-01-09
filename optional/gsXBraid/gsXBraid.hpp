/** @file gsXBraid.hpp

    @brief Provides implementations of the XBraid wrapper.

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): M. Moller
*/

#pragma once

#include <iostream>

namespace gismo {

  // Constructor
  template <typename T>
  gsXBraid<T>::gsXBraid(const gsMpiComm& comm,
                        const braid_Real tstart,
                        const braid_Real tstop,
                        braid_Int        ntime)
    : BraidApp(static_cast<MPI_Comm>(comm), tstart, tstop, ntime),
      core(static_cast<MPI_Comm>(comm), this)
  {}

  // Destructor
  template <typename T>
  gsXBraid<T>::~gsXBraid()
  {}

  // Constructor
  template <typename T>
  gsXBraid< gsMatrix<T> >::gsXBraid(const gsMpiComm& comm,
                                    const braid_Real tstart,
                                    const braid_Real tstop,
                                    braid_Int        ntime)
    : gsXBraid<T>(comm, tstart, tstop, ntime)
  {}
  
  // Destructor
  template <typename T>
  gsXBraid< gsMatrix<T> >::~gsXBraid()
  {}
  
  // Constructor
  template <typename T>
  gsXBraid< gsVector<T> >::gsXBraid(const gsMpiComm& comm,
                                    const braid_Real tstart,
                                    const braid_Real tstop,
                                    braid_Int        ntime)
    : gsXBraid<T>(comm, tstart, tstop, ntime)
  {}
  
  // Destructor
  template <typename T>
  gsXBraid< gsVector<T> >::~gsXBraid()
  {}

}// namespace gismo
