/** @file gsTensorOptions.h

    @brief Provides declarations of the gsTensorOptions class which is a
    wrapper of the LibTorch torch::TensorOptions class (https://pytorch.org)

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): M. Moller
*/

#pragma once

#pragma push_macro("index_t")
#undef index_t
#include <torch/torch.h>
#pragma pop_macro("index_t")

namespace gismo {

  /// @brief
  template <typename T>
  class gsTensorOptions : public torch::TensorOptions
  {
  public:
    // Default onstructor
    gsTensorOptions()
      : torch::TensorOptions(caffe2::TypeMeta::Make<T>()) {}

    // Copy constructor (from torch::TensorOptions)
    gsTensorOptions(const torch::TensorOptions& other)
      : torch::TensorOptions(other.dtype(caffe2::TypeMeta::Make<T>())) {}

    // Move constructor (from torch::TensorOptions)
    gsTensorOptions(torch::TensorOptions&& other)
      : torch::TensorOptions(other.dtype(caffe2::TypeMeta::Make<T>())) {}

    // Copy constructor (from gsTensorOptions)
    gsTensorOptions(const gsTensorOptions&) = default;

    // Move constructor (from gsTensorOptions)
    gsTensorOptions(gsTensorOptions&&) = default;

  public:
    // Copy assignment operator (from torch::TensorOptions)
    gsTensorOptions& operator=(const torch::TensorOptions& other)
    {
      torch::TensorOptions::operator=(other.dtype(caffe2::TypeMeta::Make<T>()));
      return *this;
    }

    // Move assignment operator (from torch::TensorOptions)
    gsTensorOptions& operator=(torch::TensorOptions&& other)
    {
      torch::TensorOptions::operator=(other.dtype(caffe2::TypeMeta::Make<T>()));
      return *this;
    }

    // Copy assignment operator (from gsTensorOptions)
    gsTensorOptions& operator=(const gsTensorOptions& other)
    {
      torch::TensorOptions::operator=(other);
      return *this;
    }

    // Move assignment operator (from gsTensorOptions)
    gsTensorOptions& operator=(gsTensorOptions&& other)
    {
      torch::TensorOptions::operator=(other);
      return *this;
    }
    
    
  public:
    // Setter functions
    gsTensorOptions setActive();
    gsTensorOptions setPassive();
    gsTensorOptions setStrided();
    gsTensorOptions setSparse();
    gsTensorOptions setCPU();
    gsTensorOptions setCUDA(int device=0);
    gsTensorOptions setPinnedMemory();
    gsTensorOptions setNonPinnedMemory();
  };
  
} // end namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsTensorOptions.hpp)
#endif

