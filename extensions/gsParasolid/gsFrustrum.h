/** @file gsFrustrum.h

    @brief Defines the Parasolud frustrim

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): A. Mantzaflaris
    Based on Parasolid templates
*/

#pragma once

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <time.h>

#ifdef PS_AIX
#include <signal.h>
#include <sys/param.h>
#endif

#ifdef _WIN32
#include <windows.h>
#endif


#include "frustrum_ifails.h"
#include "frustrum_tokens.h"

#include <parasolid_kernel.h>
// #include <kernel_interface.h>

namespace gismo {

namespace extensions {

int register_frustrum ();

}//extensions

}//gismo

#define PARASOLID_ERROR(name, err) \
    if (err) gsWarn<< "Parasolid " #name ": "<<err<<".\n" ;
