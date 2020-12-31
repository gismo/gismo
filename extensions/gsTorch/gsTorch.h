/** @file gsTorch.h

    @brief Provides declarations of Torch integration routines

    https://pytorch.org

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): M. Moller
*/

#pragma once

#ifdef index_t
#define _index_t index_t
#undef index_t
#endif

#include <torch/torch.h>

#ifdef _index_t
#define index_t _index_t
#undef _index_t
#endif

#include <gsCore/gsConfig.h>

namespace gismo {

  namespace math {

    // LibTorch types cannot be be used as native types in G+Smo and
    // are therefore not included in gsCore/gsMath.h
    // However, to make math functions accessible under the namespace
    // gismo::math::func we import them here.
    using torch::abs;
    using torch::acos;
    using torch::asin;
    using torch::atan2;
    using torch::atan;
    using torch::ceil;
    using torch::cos;
    using torch::cosh;
    using torch::exp;
    using torch::floor;
    //using torch::frexp;
    //using torch::ldexp;
    using torch::log10;
    using torch::log;
    using torch::max;
    using torch::min;
    using torch::pow;
    using torch::sin;
    using torch::sinh;
    using torch::sqrt;
    using torch::tan;
    using torch::tanh;
    using torch::real;
    using torch::imag;
    using torch::conj;

    // dummy
    inline torch::Tensor frexp(const torch::Tensor & a, int* b) { return a; }
    inline torch::Tensor ldexp(const torch::Tensor & a, int* b) { return a; }
    
  } // end namespace math
  
  // Import factory functions from namespace torch
  using torch::arange;
  using torch::bartlett_window;
  using torch::blackman_window;
  using torch::empty;
  using torch::empty_like;
  using torch::empty_meta;
  using torch::empty_strided;
  using torch::eye;
  using torch::from_blob;
  using torch::from_file;
  using torch::full;
  using torch::full_like;
  using torch::hamming_window;
  using torch::hann_window;
  using torch::kaiser_window;
  using torch::linspace;
  using torch::logspace;
  using torch::normal;
  using torch::ones;
  using torch::ones_like;
  using torch::rand;
  using torch::rand_like;
  using torch::randint;
  using torch::randint_like;
  using torch::randn;
  using torch::randn_like;
  using torch::randperm;
  using torch::range;
  using torch::requires_grad;
  using torch::scalar_tensor;
  using torch::sparse_coo_tensor;
  using torch::tril_indices;
  using torch::triu_indices;
  using torch::zeros;
  using torch::zeros_like;

  using torch::adaptive_avg_pool1d;
  using torch::adaptive_avg_pool2d;
  using torch::adaptive_avg_pool3d;
  using torch::adaptive_avg_pool3d_backward;
  using torch::adaptive_max_pool1d;
  using torch::adaptive_max_pool2d;
  using torch::adaptive_max_pool2d_backward;
  using torch::adaptive_max_pool3d;
  using torch::adaptive_max_pool3d_backward;

  using torch::add;
  using torch::addbmm;
  using torch::addcdiv;
  using torch::addcmul;
  using torch::addmm;
  using torch::addmv;
  using torch::addr;

  using torch::affine_grid_generator;
  using torch::affine_grid_generator_backward;

  using torch::alias;
  using torch::align_tensors;
  //using torch::align_to;
  using torch::all;
  using torch::allclose;
  using torch::alpha_dropout;

  using torch::amax;
  using torch::amin;
  using torch::angle;
  using torch::any;

  using torch::arccos;
  using torch::arcsin;
  using torch::arccosh;
  using torch::arcsinh;
  using torch::arctan;
  using torch::arctanh;

  using torch::argmax;
  using torch::argmin;
  using torch::argsort;

  using torch::as_strided;

  using torch::atleast_1d;
  using torch::atleast_2d;
  using torch::atleast_3d;

  using torch::avg_pool1d;
  using torch::avg_pool2d;
  using torch::avg_pool2d_backward;
  using torch::avg_pool3d;
  using torch::avg_pool3d_backward;

  //using torch::backward;

  using torch::baddbmm;
  using torch::batch_norm;
  using torch::batch_norm_backward_elemt;
  using torch::batch_norm_backward_reduce;
  using torch::batch_norm_elemt;
  using torch::batch_norm_gather_stats;
  using torch::batch_norm_gather_stats;
  using torch::batch_norm_gather_stats_with_counts;
  using torch::batch_norm_stats;
  using torch::batch_norm_update_stats;
  using torch::bernoulli;
  using torch::bilinear;
  using torch::binary_cross_entropy;
  using torch::binary_cross_entropy_backward;
  using torch::binary_cross_entropy_with_logits;
  using torch::binary_cross_entropy_with_logits_backward;
  using torch::bincount;
  using torch::binomial;
  
  

  
  
  
  
  // Create alias to torch::Tensor (do not use at:Tensor which is not differentiable)
  template <typename T> using gsTensor = torch::Tensor;
  
  // Create alias to torch::TensorOptions
  template <typename T>
  class gsTensorOptions : public torch::TensorOptions
  {
  public:
    gsTensorOptions() : torch::TensorOptions() {}
    gsTensorOptions(torch::TensorOptions other) : torch::TensorOptions(other) {}
    gsTensorOptions(const gsTensorOptions&) = default;
    gsTensorOptions(gsTensorOptions&&)      = default;

    // Core functionality
    gsTensorOptions setActive()          { return this->requires_grad(true); }
    gsTensorOptions setPassive()         { return this->requires_grad(false); }
    gsTensorOptions setStrided()         { return this->layout(torch::kStrided); }
    gsTensorOptions setSparse()          { return this->layout(torch::kSparse); }
    gsTensorOptions setCPU()             { return this->device(torch::kCPU); }
    gsTensorOptions setGPU(int device=0) { return this->device(torch::kCUDA, device); }
  };

  // Specializations for floating-point numbers
  template<> gsTensorOptions<double>::gsTensorOptions() : torch::TensorOptions(torch::kFloat64) {}
  template<> gsTensorOptions<float>::gsTensorOptions()  : torch::TensorOptions(torch::kFloat32) {}

  // Specializations for integers
  template<> gsTensorOptions<long long>::gsTensorOptions()  : torch::TensorOptions(torch::kInt64) {}
  template<> gsTensorOptions<long>::gsTensorOptions()       : torch::TensorOptions(sizeof(long)  == 32 ? torch::kInt32 : torch::kInt64) {}
  template<> gsTensorOptions<int>::gsTensorOptions()        : torch::TensorOptions(sizeof(int)   == 16 ? torch::kInt16 : torch::kInt32) {}
  template<> gsTensorOptions<short>::gsTensorOptions()      : torch::TensorOptions(sizeof(short) == 16 ? torch::kInt16 : torch::kInt32) {}

} // end namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsTorch.hpp)
#endif
