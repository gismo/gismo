
#pragma once

#include <gsParasolid/gsFrustrum.h>

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
