/** @file gsTensorOptions.hpp

    @brief Provides implementation of the gsTensorOptions class which is a
    wrapper of the LibTorch torch::TensorOptions class (https://pytorch.org)

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): M. Moller
*/

#include <gsCore/gsDebug.h>

namespace gismo
{
  // Setter functions
  template<typename T> gsTensorOptions<T> gsTensorOptions<T>::setActive()
  { return this->requires_grad(true); }

  template<typename T> gsTensorOptions<T> gsTensorOptions<T>::setPassive()
  { return this->requires_grad(false); }
  
  template<typename T> gsTensorOptions<T> gsTensorOptions<T>::setStrided()
  { return this->layout(torch::kStrided); }

  template<typename T> gsTensorOptions<T> gsTensorOptions<T>::setSparse()
  { return this->layout(torch::kSparse); }

  template<typename T> gsTensorOptions<T> gsTensorOptions<T>::setCPU()
  { return this->device(torch::kCPU); }

  template<typename T> gsTensorOptions<T> gsTensorOptions<T>::setCUDA(int device)
  { return this->device(torch::kCUDA, device); }

  template<typename T> gsTensorOptions<T> gsTensorOptions<T>::setPinnedMemory()
  { return this->pinned_memory(true); }

  template<typename T> gsTensorOptions<T> gsTensorOptions<T>::setNonPinnedMemory()
  { return this->pinned_memory(false); }
  
} // end namespace gismo
