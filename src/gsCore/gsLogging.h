
#pragma once

#include <iostream>

/** Logging messages:
 *  If the LOGGING variable is defined then gismo::gsLogging will log to the
 *  screen. Otherwise it just ignores everything that gets streamed to it.
 */

#ifdef GISMO_LOGGING_WARN
#define gsWarn std::cout<<"Warning: "
#else
#define gsWarn if (0) std::cerr
#endif

#ifdef GISMO_LOGGING_INFO
#define gsInfo std::cout  //<< "Info: "
#else
#define gsInfo std::cerr
#endif

#ifndef  NDEBUG 
//#ifdef GISMO_LOGGING_DEBUG
  #define gsDebug std::cout<<"GISMO_DEBUG: "
  #define gsDebugVar(variable) gsDebug <<"L"<<__LINE__<< ", "#variable": "<<(variable)<<"\n"
#else
  #define gsDebug if (0) std::cerr
  #define gsDebugVar(variable)
#endif

