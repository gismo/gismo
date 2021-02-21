/** @file exprtk_codi_forward.hpp

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

typedef codi::RealForwardGen<codi::codi_real<real_t>::type>                codi_real_forward_t;
typedef codi::RealReverseGen<codi::codi_real<real_t>::type>                codi_real_reverse_t;
typedef codi::RealReverseIndexGen<codi::codi_real<real_t>::type>           codi_real_reverse_index_t;
typedef codi::RealReverseIndexUncheckedGen<codi::codi_real<real_t>::type>  codi_real_reverse_index_unchecked_t;
typedef codi::RealReverseUncheckedGen<codi::codi_real<real_t>::type>       codi_real_reverse_unchecked_t;

#define CODI_TYPE codi_real_forward_t
#include "exprtk_codi_forward.h"
#undef CODI_TYPE

#define CODI_TYPE codi_real_reverse_t
#include "exprtk_codi_forward.h"
#undef CODI_TYPE

#define CODI_TYPE codi_real_reverse_index_t
#include "exprtk_codi_forward.h"
#undef CODI_TYPE

#define CODI_TYPE codi_real_reverse_index_unchecked_t
#include "exprtk_codi_forward.h"
#undef CODI_TYPE

#define CODI_TYPE codi_real_reverse_unchecked_t
#include "exprtk_codi_forward.h"
#undef CODI_TYPE
