/** @file gsTensor_test.cpp

    @brief test gsTensor

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): M. Moller
**/

#include "gismo_unittest.h"
#include <gsTorch/gsTorch.h>

using namespace gismo;

SUITE(gsTensor_test)
{
  TEST(construction_test)
  {
    {
      // Check default constructor
      gsTensor<double> tensor_double;
      CHECK(tensor_double.dtype() == caffe2::TypeMeta::Make<double>());

      gsTensor<float> tensor_float;
      CHECK(tensor_float.dtype() == caffe2::TypeMeta::Make<float>());

      // Not supported yet
      // gsTensor<long> tensor_long;
      // CHECK(tensor_long.dtype() == caffe2::TypeMeta::Make<long>());

      gsTensor<int> tensor_int;
      CHECK(tensor_int.dtype() == caffe2::TypeMeta::Make<int>());

      gsTensor<short> tensor_short;
      CHECK(tensor_short.dtype() == caffe2::TypeMeta::Make<short>());
    }

    {
      // Check copy constructor (from torch::Tensor)
      torch::Tensor tensorSrc = torch::ones(1,torch::TensorOptions(torch::kDouble));
      gsTensor<short> tensorDst(tensorSrc);
      CHECK(tensorDst.dtype() == caffe2::TypeMeta::Make<short>());
    }

    {
      // Check move constructor (from torch::Tensor)
      gsTensor<short> tensorDst(torch::ones(1,torch::TensorOptions(torch::kDouble)));
      CHECK(tensorDst.dtype() == caffe2::TypeMeta::Make<short>());
    }

    {
      // Check copy constructor (from gsTensor)
      gsTensor<short> tensorSrc;
      gsTensor<short> tensorDst(tensorSrc);
      CHECK(tensorDst.dtype() == caffe2::TypeMeta::Make<short>());
    }

    {
      // Check move constructor (from gsTensor)
      gsTensor<short> tensorDst(gsTensor<short>{});
      CHECK(tensorDst.dtype() == caffe2::TypeMeta::Make<short>());
    }
  }

  TEST(operator_test)
  {
    {
      // Check copy assignment operator (from torch::Tensor)
      torch::Tensor tensorSrc = torch::ones(1,torch::TensorOptions(torch::kDouble));
      gsTensor<short> tensorDst;
      tensorDst = tensorSrc;
      CHECK(tensorDst.dtype() == caffe2::TypeMeta::Make<short>());
    }

    {
      // Check move assignment operator (from torch::TensorOptions)
      gsTensor<short> tensorDst;
      tensorDst = torch::ones(1,torch::TensorOptions(torch::kDouble));
      CHECK(tensorDst.dtype() == caffe2::TypeMeta::Make<short>());
    }

    {
      // Check copy assignment operator (from gsTensorOptions)
      gsTensor<short> tensorSrc;
      gsTensor<short> tensorDst; tensorDst = tensorSrc;
      CHECK(tensorDst.dtype() == caffe2::TypeMeta::Make<short>());
    }

    {
      // Check move assignment operator (from gsTensorOptions)
      gsTensor<short> tensorDst;
      tensorDst = gsTensor<short>{};
      CHECK(tensorDst.dtype() == caffe2::TypeMeta::Make<short>());
    }
  }

  TEST(setter_test)
  {
    gsTensor<double> tensor;

    tensor.setRequiresGrad();
    CHECK(tensor.requires_grad() == true);
    tensor.unsetRequiresGrad();
    CHECK(tensor.requires_grad() == false);
  }
  
  TEST(xml_test)
  {
    gsTensor<double> tensor_out = torch::ones({2,2}), tensor_in;
    CHECK( (tensor_out-tensor_in).norm().item<double>() > 0 );
    
    gsFileData<> fd_out;
    fd_out << tensor_out;
    fd_out.save(gsFileManager::getTempPath()+"output.xml");

    gsFileData<> fd_in(gsFileManager::getTempPath()+"output.xml");
    fd_in.getId(0, tensor_in);
    CHECK( (tensor_out-tensor_in).norm().item<double>() < 1e-12 );
  }
}
