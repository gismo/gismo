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

/**
   \brief Class defining a tensor
*/
template <typename T>
class gsTensor : public torch::Tensor
{
public:
  using torch::Tensor::Tensor;
};

using torch::eye;

} //namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsTorch.hpp)
#endif
