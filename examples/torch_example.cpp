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
#include <gsCoDiPack/gsCoDiPack.h>

using namespace gismo;

int main(int argc, char *argv[])
{
  {
    gsTensor<real_t> x = ones({2, 2}, gsTensorOptions<real_t>().setActive());
    auto y = x + x + x;
    auto z = torch::mm(y,y) * 3;
    auto f = z.mean();
    
    f.backward();
    gsInfo << "f(x)  : " << f << std::endl;
    gsInfo << "df/dx : " << x.grad() << std::endl;
  }

  {
    codi::RealReverse::TapeType& tape = codi::RealReverse::getGlobalTape();
    tape.setActive();

    gsMatrix<codi::RealReverse> x(2,2); x << 1,1,1,1;
    tape.registerInput(x(0,0));
    tape.registerInput(x(0,1));
    tape.registerInput(x(1,0));
    tape.registerInput(x(1,1));
    
    auto y = x + x + x;
    auto z = y * y * 3;
    codi::RealReverse f = z.mean();
    tape.registerOutput(f);
    
    tape.setPassive();
    f.setGradient(1.0);
    tape.evaluate();

    gsInfo << "f(x)  : " << f << std::endl;
    gsInfo << "df/dx : "
           << x(0,0).getGradient() << "," << x(0,1).getGradient() << ";"
           << x(1,0).getGradient() << "," << x(1,1).getGradient()
           << std::endl;
  }
  
  return 0;
}
