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
  template<typename T> gsTensorOptions<T>& gsTensorOptions<T>::setRequiresGrad()
  {
    *this = requires_grad(true);
    return *this;
  }

  template<typename T> gsTensorOptions<T>& gsTensorOptions<T>::unsetRequiresGrad()
  {
    *this = requires_grad(false);
    return *this;
  }
  
  template<typename T> gsTensorOptions<T>& gsTensorOptions<T>::setStrided()
  {
    *this = layout(torch::kStrided);
    return *this;
  }

  template<typename T> gsTensorOptions<T>& gsTensorOptions<T>::setSparse()
  {
    *this = layout(torch::kSparse);
    return *this;
  }

  template<typename T> gsTensorOptions<T>& gsTensorOptions<T>::setSparseCsr()
  {
    *this = layout(torch::kSparseCsr);
    return *this;
  }

  template<typename T> gsTensorOptions<T>& gsTensorOptions<T>::setCPU()
  {
    *this = device(torch::kCPU);
    return *this;
  }

  template<typename T> gsTensorOptions<T>& gsTensorOptions<T>::setCUDA(int index)
  {
    *this = device(torch::kCUDA, index);
    return *this;
  }

  template<typename T> gsTensorOptions<T>& gsTensorOptions<T>::setPinnedMemory()
  {
    *this = pinned_memory(true);
    return *this;
  }

  template<typename T> gsTensorOptions<T>& gsTensorOptions<T>::unsetPinnedMemory()
  {
    *this = pinned_memory(false);
    return *this;
  }

  template<typename T> gsTensorOptions<T>& gsTensorOptions<T>::setMemoryFormatPreserve()
  {
    *this = memory_format(at::MemoryFormat::Preserve);
    return *this;
  }
  
  template<typename T> gsTensorOptions<T>& gsTensorOptions<T>::setMemoryFormatContiguous()
  {
    *this = memory_format(at::MemoryFormat::Contiguous);
    return *this;
  }
  
  template<typename T> gsTensorOptions<T>& gsTensorOptions<T>::setMemoryFormatChannelsLast()
  {
    *this = memory_format(at::MemoryFormat::ChannelsLast);
    return *this;
  }
  
  template<typename T> gsTensorOptions<T>& gsTensorOptions<T>::setMemoryFormatChannelsLast3d()
  {
    *this = memory_format(at::MemoryFormat::ChannelsLast3d);
    return *this;
  }

  template<typename T> gsTensorOptions<T>& gsTensorOptions<T>::unsetMemoryFormat()
  {
    *this = memory_format(c10::nullopt);
    return *this;
  }

  // Comparison operators
  template<typename T> bool gsTensorOptions<T>::operator==(const gsTensorOptions& other) const
  {
    return (
            (dtype() == other.dtype()) &&
            (requires_grad() == other.requires_grad()) &&
            (device() == other.device()) &&
            (layout() == other.layout()) &&
            (pinned_memory() == other.pinned_memory())
            );
  }

  template<typename T> bool gsTensorOptions<T>::operator!=(const gsTensorOptions& other) const
  {
    return !(*this==other);
  }  
  
} // end namespace gismo
