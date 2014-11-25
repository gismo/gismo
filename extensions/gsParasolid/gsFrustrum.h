
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

int register_frustrum ();

#define PARASOLID_ERROR(name, err) \
    if (err) gsWarn<< "Parasolid " #name ": "<<err<<".\n" ;
