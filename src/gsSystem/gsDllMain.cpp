/** @file gsDllMain.cpp

    @brief Initialization of the dynamic library: loads the modules

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): A. Mantzaflaris
*/

#include <gsSystem/gsModuleManager.h>

namespace
{
// Loads every module found in the search paths of the module manager
void gsLoadModules()
{ gismo::gsModuleManager::get().loadAll(); }
}

#ifdef _WIN32

#include <windows.h>

// Any library initialization code goes here
BOOL APIENTRY
DllMain( HANDLE /* hModule */,
         DWORD  ul_reason_for_call,
         LPVOID /* lpReserved */  )
{
    if ( DLL_PROCESS_ATTACH == ul_reason_for_call )
        gsLoadModules();
    return TRUE;
}

#else

// Runs when the library is loaded, before main() of a linked executable
__attribute__((constructor))
static void gsLibraryInit()
{ gsLoadModules(); }

#endif

// Note: Microsoft diadvises calling LoadLibrary from DllMain
// https://learn.microsoft.com/en-us/windows/win32/dlls/dllmain
