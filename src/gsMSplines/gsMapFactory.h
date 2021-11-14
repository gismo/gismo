/** @file gsMapFactory.h

    @brief Provides declaration of abstract interface of gsMapFactory.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): F. Buchegger
*/

#pragma once

#include <gsMSplines/gsWeightMapper.h>

namespace gismo
{

/** @brief
      class to construct objects of type gsWeightMapper
  */

class gsMapFactory
{
public:

    gsMapFactory(){ }

    virtual ~gsMapFactory() { }

    virtual gsWeightMapper<real_t> * makeMapper() const = 0;

}; // gsMapFactory

}
