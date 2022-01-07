/** @file gsTensor.h

    @brief Provides declarations of the gsTensor class which is a
    wrapper of the LibTorch torch::Tensor class (https://pytorch.org)

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): M. Moller
*/

#pragma once

#include <gsTorch/gsTensorOptions.h>

namespace gismo {

  template <typename T>
  class gsTensor : public torch::Tensor
  {
  public:
    // Default constructor
    gsTensor()
      : torch::Tensor(torch::empty({}, caffe2::TypeMeta::Make<T>())) {}

    // Copy constructor (from torch::Tensor)
    gsTensor(const torch::Tensor& other)
      : torch::Tensor(other.to(caffe2::TypeMeta::Make<T>())) {}

    // Move constructor (from torch::Tensor)
    gsTensor(torch::Tensor&& other)
      : torch::Tensor(other.to(caffe2::TypeMeta::Make<T>())) {}

    // Copy constructor (from gsTensor)
    gsTensor(const gsTensor&) = default;

    // Move constructor (from gsTensor)
    gsTensor(gsTensor&&) = default;

    // Constructor from dimension array
    gsTensor(at::IntArrayRef size)
      : torch::Tensor(torch::empty(size, caffe2::TypeMeta::Make<T>())) {}

  public:
    // Copy assignment operator (from torch::Tensor)
    gsTensor& operator=(const torch::Tensor& other)
    {
      torch::Tensor::operator=(other.to(caffe2::TypeMeta::Make<T>()));
      return *this;
    }

    // Move assignment operator (from torch::Tensor)
    gsTensor& operator=(torch::Tensor&& other)
    {
      torch::Tensor::operator=(other.to(caffe2::TypeMeta::Make<T>()));
      return *this;
    }

    // Copy assignment operator (from torch::Tensor)
    gsTensor& operator=(const gsTensor& other)
    {
      torch::Tensor::operator=(other);
      return *this;
    }

    // Move assignment operator (from torch::Tensor)
    gsTensor& operator=(gsTensor&& other)
    {
      torch::Tensor::operator=(other);
      return *this;
    }
    
  public:
    // Query functions
    bool isActive() const;
    bool isPassive() const;
    bool isStrided() const;
    bool isSparse() const;
    bool isCPU() const;
    bool isCUDA(int index=0) const;
    bool isPinnedMemory() const;
    bool isNonPinnedMemory() const;

    // Setter functions
    gsTensor& setActive();
    gsTensor& setPassive();
  };

} // end namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsTensor.hpp)
#endif
