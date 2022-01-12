/** @file gsTensor.hpp

    @brief Provides implementation of the gsTensor class which is a
    wrapper of the LibTorch torch::Tensor class (https://pytorch.org)

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): M. Moller
*/

namespace gismo
{
  // Query functions
  template<typename T> bool gsTensor<T>::isStrided() const
  { return this->layout() == torch::kStrided; }

  template<typename T> bool gsTensor<T>::isSparse() const
  { return this->layout() == torch::kSparse; }

  template<typename T> bool gsTensor<T>::isSparseCsr() const
  { return this->layout() == torch::kSparseCsr; }

  template<typename T> bool gsTensor<T>::isCPU() const
  { return this->device().type() == torch::kCPU; }

  template<typename T> bool gsTensor<T>::isCUDA(int index) const
  { return this->device().type() == torch::kCUDA
      &&   this->device().index() == index; }

  template<typename T> bool gsTensor<T>::isPinnedMemory() const
  { return this->is_pinned(); }
  
  // Setter functions
  template<typename T> gsTensor<T>& gsTensor<T>::setRequiresGrad()
  { this->requires_grad_(true);
    return *this; }

  template<typename T> gsTensor<T>& gsTensor<T>::unsetRequiresGrad()
  { this->requires_grad_(false);
    return *this; }

} // end namespace gismo
