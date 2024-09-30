/** @file exprtk_codi_adaptor.hpp

    @brief Provides an exprtk adaptor for CoDiPack
    arithmetic types of autodiff

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): M. Moller
*/

#pragma once

#include <string>
#include <gsCoDiPack/gsCoDiPack.h>

#include "exprtk_codi_forward.hpp"
#include "exprtk.hpp"

namespace exprtk
{
namespace details
{

namespace constant_codi
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
} // namespace constant_codi
}
}

#define CODI_TYPE codi_real_forward_t
#include "exprtk_codi_adaptor.h"
#undef CODI_TYPE

#define CODI_TYPE codi_real_reverse_t
#include "exprtk_codi_adaptor.h"
#undef CODI_TYPE

#define CODI_TYPE codi_real_reverse_index_t
#include "exprtk_codi_adaptor.h"
#undef CODI_TYPE

#define CODI_TYPE codi_real_reverse_primal_t
#include "exprtk_codi_adaptor.h"
#undef CODI_TYPE

#define CODI_TYPE codi_real_reverse_primal_index_t
#include "exprtk_codi_adaptor.h"
#undef CODI_TYPE

// The unchecked versions lead to code redefinition

// #define CODI_TYPE codi_real_reverse_unchecked_t
// #include "exprtk_codi_adaptor.h"
// #undef CODI_TYPE

// #define CODI_TYPE codi_real_reverse_index_unchecked_t
// #include "exprtk_codi_adaptor.h"
// #undef CODI_TYPE

// #define CODI_TYPE codi_real_reverse_primal_unchecked_t
// #include "exprtk_codi_adaptor.h"
// #undef CODI_TYPE

// #define CODI_TYPE codi_real_reverse_primal_index_unchecked_t
// #include "exprtk_codi_adaptor.h"
// #undef CODI_TYPE
