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
//#include <gsCoDiPack/gsCoDiPack.h>

using namespace gismo;

int main(int argc, char *argv[])
{
  {
    gsTensorOptions<short> ao;
    gsInfo << ao << std::endl;

    torch::TensorOptions bo;

    bo = bo.requires_grad(true).pinned_memory(true);
    gsInfo << bo << std::endl;

    gsTensorOptions<short> co(bo);
    gsInfo << co << std::endl;

    
    
    return 0;
    
    gsTensor<short> a({10,20});
    gsInfo << a.dtype() << std::endl;
    
    auto b(a);      
    gsInfo << b.dtype() << std::endl;

    torch::Tensor c = torch::ones({2, 2});
    gsInfo << c.dtype() << std::endl;

    gsTensor<short> d(c);
    gsInfo << d.dtype() << std::endl;
    
    return 0;
    
    gsTensor<real_t> x = torch::ones({2, 2}, gsTensorOptions<real_t>().setActive());

    gsInfo << x.dtype() << std::endl;
    
    auto y = x + x + x;
    auto z = torch::mm(y,y) * 3;
    auto f = z.mean();
    
    f.backward();
    gsInfo << "f(x)  : " << f << std::endl;
    gsInfo << "df/dx : " << x.grad() << std::endl;

    gsInfo << "Active " << static_cast<gsTensor<real_t>>(f).isActive() << "\n";
    gsInfo << "Passive " << static_cast<gsTensor<real_t>>(f).isPassive() << "\n";
    gsInfo << "Strided " << static_cast<gsTensor<real_t>>(f).isStrided() << "\n";
    gsInfo << "Sparse " << static_cast<gsTensor<real_t>>(f).isSparse() << "\n";
    gsInfo << "CPU " << static_cast<gsTensor<real_t>>(f).isCPU() << "\n";
    gsInfo << "CUDA " << static_cast<gsTensor<real_t>>(f).isCUDA() << "\n";
  }

#if 0
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
#endif
  
  return 0;
}
