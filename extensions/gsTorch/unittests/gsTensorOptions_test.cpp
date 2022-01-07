/** @file gsTensorOptions_test.cpp

    @brief test gsTensorOptions

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): M. Moller
**/

#include "gismo_unittest.h"
#include <gsTorch/gsTorch.h>

using namespace gismo;

SUITE(gsTensorOptions_test)
{
  TEST(construction_test)
  {
    {
      // Check default constructor
      gsTensorOptions<double> to_double;
      CHECK(to_double.dtype() == caffe2::TypeMeta::Make<double>());
      
      gsTensorOptions<float> to_float;
      CHECK(to_float.dtype() == caffe2::TypeMeta::Make<float>());

      gsTensorOptions<long> to_long;
      CHECK(to_long.dtype() == caffe2::TypeMeta::Make<long>());
      
      gsTensorOptions<int> to_int;
      CHECK(to_int.dtype() == caffe2::TypeMeta::Make<int>());
      
      gsTensorOptions<short> to_short;
      CHECK(to_short.dtype() == caffe2::TypeMeta::Make<short>());
    }
    
    {
      // Check copy constructor (from torch::TensorOptions)
      torch::TensorOptions to_src(torch::kDouble);
      gsTensorOptions<short> to_dst(to_src);
      CHECK(to_dst.dtype() == caffe2::TypeMeta::Make<short>());
    }

    {
      // Check move constructor (from torch::TensorOptions)
      gsTensorOptions<short> to_dst(torch::TensorOptions{}.dtype(torch::kDouble));
      CHECK(to_dst.dtype() == caffe2::TypeMeta::Make<short>());
    }

    {
      // Check copy constructor (from gsTensorOptions)
      gsTensorOptions<short> to_src;
      gsTensorOptions<short> to_dst(to_src);
      CHECK(to_dst.dtype() == caffe2::TypeMeta::Make<short>());
    }

    {
      // Check move constructor (from torch::TensorOptions)
      gsTensorOptions<short> to_dst( gsTensorOptions<short>{} );
      CHECK(to_dst.dtype() == caffe2::TypeMeta::Make<short>());
    }
  }

  TEST(operator_test)
  {
    {
      // Check copy assignment (from torch::TensorOptions)
      torch::TensorOptions to_src(torch::kDouble);
      gsTensorOptions<short> to_dst; to_dst = to_src;
      CHECK(to_dst.dtype() == caffe2::TypeMeta::Make<short>());
    }
  }
}
