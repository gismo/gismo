/** @file exprtk_universal_forward.hpp

    @brief Provides an exprtk adaptor for the Unum Posit arithmetic
    type

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): M. Moller
*/

#pragma once

#include <gsUniversal/gsUniversal.h>

typedef sw::universal::posit<256,5> posit_256_5;
typedef sw::universal::posit<128,4> posit_128_4;
typedef sw::universal::posit< 64,3> posit_64_3;
typedef sw::universal::posit< 32,2> posit_32_2;
typedef sw::universal::posit< 16,1> posit_16_1;
typedef sw::universal::posit<  8,1> posit_8_1;
typedef sw::universal::posit<  8,0> posit_8_0;
typedef sw::universal::posit<  4,0> posit_4_0;
typedef sw::universal::posit<  3,1> posit_3_1;
typedef sw::universal::posit<  3,0> posit_3_0;
typedef sw::universal::posit<  2,0> posit_2_0;

#define UNIVERSAL_TYPE posit_256_5
#include "exprtk_universal_forward.h"
#undef UNIVERSAL_TYPE

#define UNIVERSAL_TYPE posit_128_4
#include "exprtk_universal_forward.h"
#undef UNIVERSAL_TYPE

#define UNIVERSAL_TYPE posit_64_3
#include "exprtk_universal_forward.h"
#undef UNIVERSAL_TYPE

#define UNIVERSAL_TYPE posit_32_2
#include "exprtk_universal_forward.h"
#undef UNIVERSAL_TYPE

#define UNIVERSAL_TYPE posit_16_1
#include "exprtk_universal_forward.h"
#undef UNIVERSAL_TYPE

#define UNIVERSAL_TYPE posit_8_1
#include "exprtk_universal_forward.h"
#undef UNIVERSAL_TYPE

#define UNIVERSAL_TYPE posit_8_0
#include "exprtk_universal_forward.h"
#undef UNIVERSAL_TYPE

#define UNIVERSAL_TYPE posit_4_0
#include "exprtk_universal_forward.h"
#undef UNIVERSAL_TYPE

#define UNIVERSAL_TYPE posit_3_1
#include "exprtk_universal_forward.h"
#undef UNIVERSAL_TYPE

#define UNIVERSAL_TYPE posit_3_0
#include "exprtk_universal_forward.h"
#undef UNIVERSAL_TYPE

#define UNIVERSAL_TYPE posit_2_0
#include "exprtk_universal_forward.h"
#undef UNIVERSAL_TYPE
