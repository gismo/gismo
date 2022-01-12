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
      gsTensorOptions<double> tensorOpt_double;
      CHECK(tensorOpt_double.dtype() == caffe2::TypeMeta::Make<double>());

      gsTensorOptions<float> tensorOpt_float;
      CHECK(tensorOpt_float.dtype() == caffe2::TypeMeta::Make<float>());

      gsTensorOptions<long> tensorOpt_long;
      CHECK(tensorOpt_long.dtype() == caffe2::TypeMeta::Make<long>());

      gsTensorOptions<int> tensorOpt_int;
      CHECK(tensorOpt_int.dtype() == caffe2::TypeMeta::Make<int>());

      gsTensorOptions<short> tensorOpt_short;
      CHECK(tensorOpt_short.dtype() == caffe2::TypeMeta::Make<short>());
    }

    {
      // Check copy constructor (from torch::TensorOptions)
      torch::TensorOptions tensorOptSrc(torch::kDouble);
      gsTensorOptions<short> tensorOptDst(tensorOptSrc);
      CHECK(tensorOptDst.dtype() == caffe2::TypeMeta::Make<short>());
    }

    {
      // Check move constructor (from torch::TensorOptions)
      gsTensorOptions<short> tensorOptDst(torch::TensorOptions{}.dtype(torch::kDouble));
      CHECK(tensorOptDst.dtype() == caffe2::TypeMeta::Make<short>());
    }

    {
      // Check copy constructor (from gsTensorOptions)
      gsTensorOptions<short> tensorOptSrc;
      gsTensorOptions<short> tensorOptDst(tensorOptSrc);
      CHECK(tensorOptDst.dtype() == caffe2::TypeMeta::Make<short>());
    }

    {
      // Check move constructor (from gsTensorOptions)
      gsTensorOptions<short> tensorOptDst(gsTensorOptions<short>{});
      CHECK(tensorOptDst.dtype() == caffe2::TypeMeta::Make<short>());
    }
  }

  TEST(operator_test)
  {
    {
      // Check copy assignment operator (from torch::TensorOptions)
      torch::TensorOptions tensorOptSrc(torch::kDouble);
      gsTensorOptions<short> tensorOptDst;
      tensorOptDst = tensorOptSrc;
      CHECK(tensorOptDst.dtype() == caffe2::TypeMeta::Make<short>());
    }

    {
      // Check move assignment operator (from torch::TensorOptions)
      gsTensorOptions<short> tensorOptDst;
      tensorOptDst = torch::TensorOptions{}.dtype(torch::kDouble);
      CHECK(tensorOptDst.dtype() == caffe2::TypeMeta::Make<short>());
    }

    {
      // Check copy assignment operator (from gsTensorOptions)
      gsTensorOptions<short> tensorOptSrc;
      gsTensorOptions<short> tensorOptDst;
      tensorOptDst = tensorOptSrc;
      CHECK(tensorOptDst.dtype() == caffe2::TypeMeta::Make<short>());
    }

    {
      // Check move assignment operator (from gsTensorOptions)
      gsTensorOptions<short> tensorOptDst;
      tensorOptDst = gsTensorOptions<short>{};
      CHECK(tensorOptDst.dtype() == caffe2::TypeMeta::Make<short>());
    }
  }

  TEST(setter_test)
  {
    gsTensorOptions<double> to;

    to.setRequiresGrad();
    CHECK(to.requires_grad() == true);
    to.unsetRequiresGrad();
    CHECK(to.requires_grad() == false);

    to.setStrided();
    CHECK(to.layout() == torch::kStrided);
    to.setSparse();
    CHECK(to.layout() == torch::kSparse);
    to.setSparseCsr();
    CHECK(to.layout() == torch::kSparseCsr);

    if (torch::cuda::is_available()) {
      to.setCUDA();
      CHECK(to.device() == torch::kCUDA);
    }
    to.setCPU();
    CHECK(to.device() == torch::kCPU);

    to.setPinnedMemory();
    CHECK(to.pinned_memory() == true);
    to.unsetPinnedMemory();
    CHECK(to.pinned_memory() == false);
  }

  TEST(comparison_test)
  {
    gsTensorOptions<double> tensorOptA, tensorOptB;
    CHECK(tensorOptA == tensorOptB);

    tensorOptA.unsetRequiresGrad();
    tensorOptB.setRequiresGrad();
    CHECK(tensorOptA != tensorOptB);
  }

  TEST(xml_test)
  {
    gsTensorOptions<double> tensorOptOut, tensorOptIn;
    tensorOptOut.setStrided().setRequiresGrad();
    CHECK(tensorOptOut != tensorOptIn);

    gsFileData<> fdOut;
    fdOut << tensorOptOut;
    fdOut.save(gsFileManager::getTempPath()+"output.xml");

    gsFileData<> fdIn(gsFileManager::getTempPath()+"output.xml");
    fdIn.getId(0, tensorOptIn);
    CHECK(tensorOptOut == tensorOptIn);
  }
}
