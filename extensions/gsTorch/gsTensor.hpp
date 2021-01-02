/** @file gsTensor.hpp

    @brief Provides implementation of the gsTensor class which is a
    wrapper of the PyTorch torch::tensor class (https://pytorch.org)

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): M. Moller
*/

namespace gismo
{  
  // Specializations for floating-point numbers
  template<> gsTensorOptions<double>::gsTensorOptions() : torch::TensorOptions(torch::kFloat64) {}
  template<> gsTensorOptions<float>::gsTensorOptions()  : torch::TensorOptions(torch::kFloat32) {}

  // Specializations for integers
  template<> gsTensorOptions<long long>::gsTensorOptions()  : torch::TensorOptions(torch::kInt64) {}
  template<> gsTensorOptions<long>::gsTensorOptions()       : torch::TensorOptions(sizeof(long)  == 32 ? torch::kInt32 : torch::kInt64) {}
  template<> gsTensorOptions<int>::gsTensorOptions()        : torch::TensorOptions(sizeof(int)   == 16 ? torch::kInt16 : torch::kInt32) {}
  template<> gsTensorOptions<short>::gsTensorOptions()      : torch::TensorOptions(sizeof(short) == 16 ? torch::kInt16 : torch::kInt32) {}
  
  template<typename T> gsTensorOptions<T> gsTensorOptions<T>::setActive()        { return this->requires_grad(true); }
  template<typename T> gsTensorOptions<T> gsTensorOptions<T>::setPassive()       { return this->requires_grad(false); }
  template<typename T> gsTensorOptions<T> gsTensorOptions<T>::setStrided()       { return this->layout(torch::kStrided); }
  template<typename T> gsTensorOptions<T> gsTensorOptions<T>::setSparse()        { return this->layout(torch::kSparse); }
  template<typename T> gsTensorOptions<T> gsTensorOptions<T>::setCPU()           { return this->device(torch::kCPU); }
  template<typename T> gsTensorOptions<T> gsTensorOptions<T>::setGPU(int device) { return this->device(torch::kCUDA, device); }
  
} // end namespace gismo
