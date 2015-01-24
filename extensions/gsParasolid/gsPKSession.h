/** @file gsPKSession.h

    @brief Manages starting and stopping Parasolid session

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): A. Mantzaflaris
*/

#pragma once

namespace gismo {

namespace extensions {

class gsPKSession  
{
public:
    /// Starts Parasolid session
	static bool start();

    /// Stops Parasolid session
	static bool stop();

	// static void ReturnMemory(int * nBytes, char * * memory, int * ifail);
	// static void GetMemory(int * nBytes, char * * memory, int * ifail);
	// static void StopFrustrum( int * ifail );
	// static void StartFrustrum( int * ifail );
	// static PK_ERROR_code_t PKerrorHandler( PK_ERROR_sf_t* error );
};

}//extensions

}//gismo

