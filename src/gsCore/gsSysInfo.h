/** @file gsSysInfo.h

    @brief Provides system information.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): M. Moller
*/

#pragma once

#include <gsCore/gsExport.h>
#include <gsCore/gsForwardDeclarations.h>
#include <gsCore/gsLinearAlgebra.h>

#include <string>

#if defined(_WIN32) || defined(_WIN64)
#   include <windows.h>
#elif __APPLE__
#   include <sys/utsname.h>
#   include <sys/sysctl.h>
#elif __linux__
#   include <unistd.h>
#   if defined(__x86_64__) && ( defined(__GNUC__) || defined(__clang__) || defined(__INTEL_COMPILER) || define(__SUNCC_PRO))
#      include <cpuid.h>
#   else
#      include <limits.h>
#   endif
#elif __unix__
#endif

namespace gismo
{

  class GISMO_EXPORT gsSysInfo
    {
    public:
    
    /// Returns the version of G+Smo
    static std::string getGismoVersion();

    /// Returns the version of Eigen
    static std::string getEigenVersion();

    /// Returns the version of the compiler
    static std::string getCompilerVersion();

    /// Returns the version of the C++ standard
    static std::string getCppVersion();

    /// Returns the version of the standard library
    static std::string getStdLibVersion();

    /// Returns the version of extra libraries
    static std::string getExtraLibsVersion();

    /// Returns CPU information
    static std::string getCpuInfo();

    /// Returns memory information
    static std::string getMemoryInfo();
    }; // class gsSysInfo
  
} // namespace gismo
