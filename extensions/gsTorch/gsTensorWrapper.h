/** @file gsTensorWrapper.h

    @brief Provides declarations of the gsTensorWrapper class
    
    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): M. Moller
*/

#pragma once

#include <gsTorch/gsTensor.h>

namespace gismo {

  template <typename T>
  class gsTensorWrapper : public torch::autograd::Function<gsTensorWrapper<T> >
  {
  private:
  public:
    static torch::Tensor forward(torch::autograd::AutogradContext *ctx,
                                 torch::Tensor input,
                                 torch::Tensor weight,
                                 torch::Tensor bias = torch::Tensor())
    {
      ctx->save_for_backward({input, weight, bias});
      auto output = input.mm(weight.t());
      if (bias.defined()) {
        output += bias.unsqueeze(0).expand_as(output);
      }
      return output;
    }
    
    static torch::autograd::tensor_list backward(torch::autograd::AutogradContext *ctx,
                                                 torch::autograd::tensor_list grad_outputs)
    {
      auto saved  = ctx->get_saved_variables();
      auto input  = saved[0];
      auto weight = saved[1];
      auto bias   = saved[2];
      
      auto grad_output = grad_outputs[0];
      auto grad_input  = grad_output.mm(weight);
      auto grad_weight = grad_output.t().mm(input);
      auto grad_bias   = torch::Tensor();
      if (bias.defined()) {
        grad_bias = grad_output.sum(0);
      }
      
      return {grad_input, grad_weight, grad_bias};
    }
  };
  
} // end namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsTensorWrapper.hpp)
#endif

