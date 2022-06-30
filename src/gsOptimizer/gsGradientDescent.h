/** @file gsGradientDescent.h

    @brief Provides declaration of the gradient descent method.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris
*/

#pragma once

#include <gsCore/gsLinearAlgebra.h>
#include "gdcpp.h"
//#include "lsqcpp.h"


namespace gismo
{

template<typename Objective,
         typename T = real_t,
         typename StepSize=gdc::BarzilaiBorwein<T>,
         typename Callback=gdc::NoCallback<T>,
         typename FiniteDifferences=gdc::CentralDifferences<T> >
using gsGradientDescent = gdc::GradientDescent<T, Objective, StepSize, Callback, FiniteDifferences>;

} //namespace gismo
