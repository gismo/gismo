/** @file exprtk_unum_posit_adaptor.hpp

    @brief Provides an exprtk adaptor for the Unum Posit arithmetic
    type

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): M. Moller
*/

#pragma once

#include <string>
#include <gsUnum/gsUnum.h>

#include "exprtk_unum_posit_forward.hpp"
#include "exprtk.hpp"


namespace exprtk
{
namespace details
{

namespace constant_unum_posit
{
static const double e       =  2.71828182845904523536028747135266249775724709369996;
static const double pi      =  3.14159265358979323846264338327950288419716939937510;
static const double pi_2    =  1.57079632679489661923132169163975144209858469968755;
static const double pi_4    =  0.78539816339744830961566084581987572104929234984378;
static const double pi_180  =  0.01745329251994329576923690768488612713442871888542;
static const double _1_pi   =  0.31830988618379067153776752674502872406891929148091;
static const double _2_pi   =  0.63661977236758134307553505349005744813783858296183;
static const double _180_pi = 57.29577951308232087679815481410517033240547246656443;
static const double log2    =  0.69314718055994530941723212145817656807550013436026;
static const double sqrt2   =  1.41421356237309504880168872420969807856967187537695;
} // namespace constant_unum_posit  
}
}


#define UNUM_TYPE posit_256_5
#include "exprtk_unum_posit_adaptor.h"
#undef UNUM_TYPE

#define UNUM_TYPE posit_128_4
#include "exprtk_unum_posit_adaptor.h"
#undef UNUM_TYPE

#define UNUM_TYPE posit_64_3
#include "exprtk_unum_posit_adaptor.h"
#undef UNUM_TYPE

#define UNUM_TYPE posit_32_2
#include "exprtk_unum_posit_adaptor.h"
#undef UNUM_TYPE

#define UNUM_TYPE posit_16_1
#include "exprtk_unum_posit_adaptor.h"
#undef UNUM_TYPE

#define UNUM_TYPE posit_8_1
#include "exprtk_unum_posit_adaptor.h"
#undef UNUM_TYPE

#define UNUM_TYPE posit_8_0
#include "exprtk_unum_posit_adaptor.h"
#undef UNUM_TYPE

#define UNUM_TYPE posit_4_0
#include "exprtk_unum_posit_adaptor.h"
#undef UNUM_TYPE

#define UNUM_TYPE posit_3_1
#include "exprtk_unum_posit_adaptor.h"
#undef UNUM_TYPE

#define UNUM_TYPE posit_3_0
#include "exprtk_unum_posit_adaptor.h"
#undef UNUM_TYPE

#define UNUM_TYPE posit_2_0
#include "exprtk_unum_posit_adaptor.h"
#undef UNUM_TYPE
