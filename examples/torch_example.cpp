/** @file gsTorch.cpp

    @brief Demonstrate use of Torch extension

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): M. Moller
*/

#include <iostream>

#include <gismo.h>
#include <gsTorch/gsTorch.h>

using namespace gismo;

int main(int argc, char *argv[])
{
  torch::Tensor s = torch::eye(3);
  gsTensor<real_t> t;

  gsInfo << t << std::endl;
  
  return 0;
}
