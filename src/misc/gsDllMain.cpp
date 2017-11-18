/** @file gsDllMain.cpp

    @brief Required function for MS windows DLL

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): A. Mantzaflaris
*/

#ifdef _WIN32

#include <windows.h>

// Any library initialization code goes here
BOOL APIENTRY 
DllMain( HANDLE /* hModule */,
         DWORD  /* ul_reason_for_call */, 
         LPVOID /* lpReserved */  )
{ return TRUE; }
#endif
