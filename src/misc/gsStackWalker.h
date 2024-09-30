/** @file gsStackWalker.h

    @brief Provides declaration of gsStackwalker.

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): A. Mantzaflaris
*/

#pragma once

#include <gsCore/gsConfig.h>
#include <gsCore/gsExport.h>

#include <signal.h>

#if defined(_WIN32) && !defined(NOMINMAX)
  #define NOMINMAX 1
  #include <windows.h>
  #undef NOMINMAX
#endif


namespace gismo {

namespace internal {

#if defined(__GNUC__) && !defined(_WIN32)

    /// Prints the call stack
    void gsStackWalker(void * context);

    /// Exception hook 
    void gsExceptionHook(int sig, siginfo_t * siginfo, void * context);

    /// Prints demangled GCC function names from mangled symbols
    void printGccDemangled(char * symbol);

#endif


#if defined(_WIN32)

    /// Prints the call stack
    GISMO_EXPORT void gsStackWalker(CONTEXT * context);

    /// Exception hook 
    GISMO_EXPORT LONG WINAPI gsExceptionHook(EXCEPTION_POINTERS * ExceptionInfo);

    /// Used to take control of all exceptions
    GISMO_EXPORT BOOL PreventSetUnhandledExceptionFilter();


    // see https://msdn.microsoft.com/en-us/library/5at7yxcs.aspx
    //static const int gismoCrtDbgFlag = 
    //_CrtSetDbgFlag ( _CRTDBG_ALLOC_MEM_DF | _CRTDBG_LEAK_CHECK_DF );

#endif


/// Exception handler
GISMO_EXPORT bool gsExceptionHandler();


/// Initialize exception handler for stack backtrace
static const bool gismoExceptionHandler = gsExceptionHandler();

} // end namespace internal

} // end namespace gismo


#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsStackWalker.cpp)
#endif
