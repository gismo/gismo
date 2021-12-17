/** @file gsSysInfo.cpp

    @brief Provides implemementation of system information.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): M. Moller
*/

#include <gsCore/gsSysInfo.h>
#include <gsCore/gsLinearAlgebra.h>

#include <string>

#if defined(_WIN32) || defined(_WIN64)
#   include <windows.h>
#elif __APPLE__
#   include <sys/utsname.h>
#   include <sys/sysctl.h>
#elif __linux__
#   include <unistd.h>
#   if defined(__x86_64__) && ( defined(__GNUC__) || defined(__clang__) || defined(__INTEL_COMPILER) || defined(__SUNCC_PRO))
#      include <cpuid.h>
#   else
#      include <limits.h>
#   endif
#elif __unix__
#endif

namespace gismo
{

  std::string gsSysInfo::getGismoVersion()
  {
    return util::to_string(GISMO_VERSION);
  }

  std::string gsSysInfo::getEigenVersion()
  {
    return util::to_string(EIGEN_WORLD_VERSION)+"."
      +    util::to_string(EIGEN_MAJOR_VERSION)+"."
      +    util::to_string(EIGEN_MINOR_VERSION);
  }

  std::string gsSysInfo::getCompilerVersion()
  {
    // This code is copied from the CMakeCXXCompilerId.cpp file that was
    // automatically generated with CMake 3.21.4

    // The following two macros have been modified as we do not want to
    // return the compiler version in the specific CMake format
#define DEC(n) n
#define HEX(n) n

    /* Version number components: V=Version, R=Revision, P=Patch
       Version date components:   YYYY=Year, MM=Month,   DD=Day  */

#if defined(__COMO__)
# define COMPILER_ID "Comeau"
    /* __COMO_VERSION__ = VRR */
# define COMPILER_VERSION_MAJOR DEC(__COMO_VERSION__ / 100)
# define COMPILER_VERSION_MINOR DEC(__COMO_VERSION__ % 100)

#elif defined(__INTEL_COMPILER) || defined(__ICC)
# define COMPILER_ID "Intel"
# if defined(_MSC_VER)
#  define SIMULATE_ID "MSVC"
# endif
# if defined(__GNUC__)
#  define SIMULATE_ID "GNU"
# endif
    /* __INTEL_COMPILER = VRP prior to 2021, and then VVVV for 2021 and later,
       except that a few beta releases use the old format with V=2021.  */
# if __INTEL_COMPILER < 2021 || __INTEL_COMPILER == 202110 || __INTEL_COMPILER == 202111
#  define COMPILER_VERSION_MAJOR DEC(__INTEL_COMPILER/100)
#  define COMPILER_VERSION_MINOR DEC(__INTEL_COMPILER/10 % 10)
#  if defined(__INTEL_COMPILER_UPDATE)
#   define COMPILER_VERSION_PATCH DEC(__INTEL_COMPILER_UPDATE)
#  else
#   define COMPILER_VERSION_PATCH DEC(__INTEL_COMPILER   % 10)
#  endif
# else
#  define COMPILER_VERSION_MAJOR DEC(__INTEL_COMPILER)
#  define COMPILER_VERSION_MINOR DEC(__INTEL_COMPILER_UPDATE)
    /* The third version component from --version is an update index,
       but no macro is provided for it.  */
#  define COMPILER_VERSION_PATCH DEC(0)
# endif
# if defined(__INTEL_COMPILER_BUILD_DATE)
    /* __INTEL_COMPILER_BUILD_DATE = YYYYMMDD */
#  define COMPILER_VERSION_TWEAK DEC(__INTEL_COMPILER_BUILD_DATE)
# endif
# if defined(_MSC_VER)
    /* _MSC_VER = VVRR */
#  define SIMULATE_VERSION_MAJOR DEC(_MSC_VER / 100)
#  define SIMULATE_VERSION_MINOR DEC(_MSC_VER % 100)
# endif
# if defined(__GNUC__)
#  define SIMULATE_VERSION_MAJOR DEC(__GNUC__)
# elif defined(__GNUG__)
#  define SIMULATE_VERSION_MAJOR DEC(__GNUG__)
# endif
# if defined(__GNUC_MINOR__)
#  define SIMULATE_VERSION_MINOR DEC(__GNUC_MINOR__)
# endif
# if defined(__GNUC_PATCHLEVEL__)
#  define SIMULATE_VERSION_PATCH DEC(__GNUC_PATCHLEVEL__)
# endif

#elif (defined(__clang__) && defined(__INTEL_CLANG_COMPILER)) || defined(__INTEL_LLVM_COMPILER)
# define COMPILER_ID "IntelLLVM"
#if defined(_MSC_VER)
# define SIMULATE_ID "MSVC"
#endif
#if defined(__GNUC__)
# define SIMULATE_ID "GNU"
#endif
    /* __INTEL_LLVM_COMPILER = VVVVRP prior to 2021.2.0, VVVVRRPP for 2021.2.0 and
     * later.  Look for 6 digit vs. 8 digit version number to decide encoding.
     * VVVV is no smaller than the current year when a version is released.
     */
#if __INTEL_LLVM_COMPILER < 1000000L
# define COMPILER_VERSION_MAJOR DEC(__INTEL_LLVM_COMPILER/100)
# define COMPILER_VERSION_MINOR DEC(__INTEL_LLVM_COMPILER/10 % 10)
# define COMPILER_VERSION_PATCH DEC(__INTEL_LLVM_COMPILER    % 10)
#else
# define COMPILER_VERSION_MAJOR DEC(__INTEL_LLVM_COMPILER/10000)
# define COMPILER_VERSION_MINOR DEC(__INTEL_LLVM_COMPILER/100 % 100)
# define COMPILER_VERSION_PATCH DEC(__INTEL_LLVM_COMPILER     % 100)
#endif
#if defined(_MSC_VER)
    /* _MSC_VER = VVRR */
# define SIMULATE_VERSION_MAJOR DEC(_MSC_VER / 100)
# define SIMULATE_VERSION_MINOR DEC(_MSC_VER % 100)
#endif
#if defined(__GNUC__)
# define SIMULATE_VERSION_MAJOR DEC(__GNUC__)
#elif defined(__GNUG__)
# define SIMULATE_VERSION_MAJOR DEC(__GNUG__)
#endif
#if defined(__GNUC_MINOR__)
# define SIMULATE_VERSION_MINOR DEC(__GNUC_MINOR__)
#endif
#if defined(__GNUC_PATCHLEVEL__)
# define SIMULATE_VERSION_PATCH DEC(__GNUC_PATCHLEVEL__)
#endif

#elif defined(__PATHCC__)
# define COMPILER_ID "PathScale"
# define COMPILER_VERSION_MAJOR DEC(__PATHCC__)
# define COMPILER_VERSION_MINOR DEC(__PATHCC_MINOR__)
# if defined(__PATHCC_PATCHLEVEL__)
#  define COMPILER_VERSION_PATCH DEC(__PATHCC_PATCHLEVEL__)
# endif

#elif defined(__BORLANDC__) && defined(__CODEGEARC_VERSION__)
# define COMPILER_ID "Embarcadero"
# define COMPILER_VERSION_MAJOR HEX(__CODEGEARC_VERSION__>>24 & 0x00FF)
# define COMPILER_VERSION_MINOR HEX(__CODEGEARC_VERSION__>>16 & 0x00FF)
# define COMPILER_VERSION_PATCH DEC(__CODEGEARC_VERSION__     & 0xFFFF)

#elif defined(__BORLANDC__)
# define COMPILER_ID "Borland"
    /* __BORLANDC__ = 0xVRR */
# define COMPILER_VERSION_MAJOR HEX(__BORLANDC__>>8)
# define COMPILER_VERSION_MINOR HEX(__BORLANDC__ & 0xFF)

#elif defined(__WATCOMC__) && __WATCOMC__ < 1200
# define COMPILER_ID "Watcom"
    /* __WATCOMC__ = VVRR */
# define COMPILER_VERSION_MAJOR DEC(__WATCOMC__ / 100)
# define COMPILER_VERSION_MINOR DEC((__WATCOMC__ / 10) % 10)
# if (__WATCOMC__ % 10) > 0
#  define COMPILER_VERSION_PATCH DEC(__WATCOMC__ % 10)
# endif

#elif defined(__WATCOMC__)
# define COMPILER_ID "OpenWatcom"
    /* __WATCOMC__ = VVRP + 1100 */
# define COMPILER_VERSION_MAJOR DEC((__WATCOMC__ - 1100) / 100)
# define COMPILER_VERSION_MINOR DEC((__WATCOMC__ / 10) % 10)
# if (__WATCOMC__ % 10) > 0
#  define COMPILER_VERSION_PATCH DEC(__WATCOMC__ % 10)
# endif

#elif defined(__SUNPRO_CC)
# define COMPILER_ID "SunPro"
# if __SUNPRO_CC >= 0x5100
    /* __SUNPRO_CC = 0xVRRP */
#  define COMPILER_VERSION_MAJOR HEX(__SUNPRO_CC>>12)
#  define COMPILER_VERSION_MINOR HEX(__SUNPRO_CC>>4 & 0xFF)
#  define COMPILER_VERSION_PATCH HEX(__SUNPRO_CC    & 0xF)
# else
    /* __SUNPRO_CC = 0xVRP */
#  define COMPILER_VERSION_MAJOR HEX(__SUNPRO_CC>>8)
#  define COMPILER_VERSION_MINOR HEX(__SUNPRO_CC>>4 & 0xF)
#  define COMPILER_VERSION_PATCH HEX(__SUNPRO_CC    & 0xF)
# endif

#elif defined(__HP_aCC)
# define COMPILER_ID "HP"
    /* __HP_aCC = VVRRPP */
# define COMPILER_VERSION_MAJOR DEC(__HP_aCC/10000)
# define COMPILER_VERSION_MINOR DEC(__HP_aCC/100 % 100)
# define COMPILER_VERSION_PATCH DEC(__HP_aCC     % 100)

#elif defined(__DECCXX)
# define COMPILER_ID "Compaq"
    /* __DECCXX_VER = VVRRTPPPP */
# define COMPILER_VERSION_MAJOR DEC(__DECCXX_VER/10000000)
# define COMPILER_VERSION_MINOR DEC(__DECCXX_VER/100000  % 100)
# define COMPILER_VERSION_PATCH DEC(__DECCXX_VER         % 10000)

#elif defined(__IBMCPP__) && defined(__COMPILER_VER__)
# define COMPILER_ID "zOS"
    /* __IBMCPP__ = VRP */
# define COMPILER_VERSION_MAJOR DEC(__IBMCPP__/100)
# define COMPILER_VERSION_MINOR DEC(__IBMCPP__/10 % 10)
# define COMPILER_VERSION_PATCH DEC(__IBMCPP__    % 10)

#elif defined(__ibmxl__) && defined(__clang__)
# define COMPILER_ID "XLClang"
# define COMPILER_VERSION_MAJOR DEC(__ibmxl_version__)
# define COMPILER_VERSION_MINOR DEC(__ibmxl_release__)
# define COMPILER_VERSION_PATCH DEC(__ibmxl_modification__)
# define COMPILER_VERSION_TWEAK DEC(__ibmxl_ptf_fix_level__)


#elif defined(__IBMCPP__) && !defined(__COMPILER_VER__) && __IBMCPP__ >= 800
# define COMPILER_ID "XL"
    /* __IBMCPP__ = VRP */
# define COMPILER_VERSION_MAJOR DEC(__IBMCPP__/100)
# define COMPILER_VERSION_MINOR DEC(__IBMCPP__/10 % 10)
# define COMPILER_VERSION_PATCH DEC(__IBMCPP__    % 10)

#elif defined(__IBMCPP__) && !defined(__COMPILER_VER__) && __IBMCPP__ < 800
# define COMPILER_ID "VisualAge"
    /* __IBMCPP__ = VRP */
# define COMPILER_VERSION_MAJOR DEC(__IBMCPP__/100)
# define COMPILER_VERSION_MINOR DEC(__IBMCPP__/10 % 10)
# define COMPILER_VERSION_PATCH DEC(__IBMCPP__    % 10)

#elif defined(__NVCOMPILER)
# define COMPILER_ID "NVHPC"
# define COMPILER_VERSION_MAJOR DEC(__NVCOMPILER_MAJOR__)
# define COMPILER_VERSION_MINOR DEC(__NVCOMPILER_MINOR__)
# if defined(__NVCOMPILER_PATCHLEVEL__)
#  define COMPILER_VERSION_PATCH DEC(__NVCOMPILER_PATCHLEVEL__)
# endif

#elif defined(__PGI)
# define COMPILER_ID "PGI"
# define COMPILER_VERSION_MAJOR DEC(__PGIC__)
# define COMPILER_VERSION_MINOR DEC(__PGIC_MINOR__)
# if defined(__PGIC_PATCHLEVEL__)
#  define COMPILER_VERSION_PATCH DEC(__PGIC_PATCHLEVEL__)
# endif

#elif defined(_CRAYC)
# define COMPILER_ID "Cray"
# define COMPILER_VERSION_MAJOR DEC(_RELEASE_MAJOR)
# define COMPILER_VERSION_MINOR DEC(_RELEASE_MINOR)

#elif defined(__TI_COMPILER_VERSION__)
# define COMPILER_ID "TI"
    /* __TI_COMPILER_VERSION__ = VVVRRRPPP */
# define COMPILER_VERSION_MAJOR DEC(__TI_COMPILER_VERSION__/1000000)
# define COMPILER_VERSION_MINOR DEC(__TI_COMPILER_VERSION__/1000   % 1000)
# define COMPILER_VERSION_PATCH DEC(__TI_COMPILER_VERSION__        % 1000)

#elif defined(__CLANG_FUJITSU)
# define COMPILER_ID "FujitsuClang"
# define COMPILER_VERSION_MAJOR DEC(__FCC_major__)
# define COMPILER_VERSION_MINOR DEC(__FCC_minor__)
# define COMPILER_VERSION_PATCH DEC(__FCC_patchlevel__)
# define COMPILER_VERSION_INTERNAL_STR __clang_version__


#elif defined(__FUJITSU)
# define COMPILER_ID "Fujitsu"
# if defined(__FCC_version__)
#   define COMPILER_VERSION __FCC_version__
# elif defined(__FCC_major__)
#   define COMPILER_VERSION_MAJOR DEC(__FCC_major__)
#   define COMPILER_VERSION_MINOR DEC(__FCC_minor__)
#   define COMPILER_VERSION_PATCH DEC(__FCC_patchlevel__)
# endif
# if defined(__fcc_version)
#   define COMPILER_VERSION_INTERNAL DEC(__fcc_version)
# elif defined(__FCC_VERSION)
#   define COMPILER_VERSION_INTERNAL DEC(__FCC_VERSION)
# endif


#elif defined(__ghs__)
# define COMPILER_ID "GHS"
    /* __GHS_VERSION_NUMBER = VVVVRP */
# ifdef __GHS_VERSION_NUMBER
# define COMPILER_VERSION_MAJOR DEC(__GHS_VERSION_NUMBER / 100)
# define COMPILER_VERSION_MINOR DEC(__GHS_VERSION_NUMBER / 10 % 10)
# define COMPILER_VERSION_PATCH DEC(__GHS_VERSION_NUMBER      % 10)
# endif

#elif defined(__SCO_VERSION__)
# define COMPILER_ID "SCO"

#elif defined(__ARMCC_VERSION) && !defined(__clang__)
# define COMPILER_ID "ARMCC"
#if __ARMCC_VERSION >= 1000000
    /* __ARMCC_VERSION = VRRPPPP */
# define COMPILER_VERSION_MAJOR DEC(__ARMCC_VERSION/1000000)
# define COMPILER_VERSION_MINOR DEC(__ARMCC_VERSION/10000 % 100)
# define COMPILER_VERSION_PATCH DEC(__ARMCC_VERSION     % 10000)
#else
    /* __ARMCC_VERSION = VRPPPP */
# define COMPILER_VERSION_MAJOR DEC(__ARMCC_VERSION/100000)
# define COMPILER_VERSION_MINOR DEC(__ARMCC_VERSION/10000 % 10)
# define COMPILER_VERSION_PATCH DEC(__ARMCC_VERSION    % 10000)
#endif


#elif defined(__clang__) && defined(__apple_build_version__)
# define COMPILER_ID "AppleClang"
# if defined(_MSC_VER)
#  define SIMULATE_ID "MSVC"
# endif
# define COMPILER_VERSION_MAJOR DEC(__clang_major__)
# define COMPILER_VERSION_MINOR DEC(__clang_minor__)
# define COMPILER_VERSION_PATCH DEC(__clang_patchlevel__)
# if defined(_MSC_VER)
    /* _MSC_VER = VVRR */
#  define SIMULATE_VERSION_MAJOR DEC(_MSC_VER / 100)
#  define SIMULATE_VERSION_MINOR DEC(_MSC_VER % 100)
# endif
# define COMPILER_VERSION_TWEAK DEC(__apple_build_version__)

#elif defined(__clang__) && defined(__ARMCOMPILER_VERSION)
# define COMPILER_ID "ARMClang"
# define COMPILER_VERSION_MAJOR DEC(__ARMCOMPILER_VERSION/1000000)
# define COMPILER_VERSION_MINOR DEC(__ARMCOMPILER_VERSION/10000 % 100)
# define COMPILER_VERSION_PATCH DEC(__ARMCOMPILER_VERSION     % 10000)
# define COMPILER_VERSION_INTERNAL DEC(__ARMCOMPILER_VERSION)

#elif defined(__clang__)
# define COMPILER_ID "Clang"
# if defined(_MSC_VER)
#  define SIMULATE_ID "MSVC"
# endif
# define COMPILER_VERSION_MAJOR DEC(__clang_major__)
# define COMPILER_VERSION_MINOR DEC(__clang_minor__)
# define COMPILER_VERSION_PATCH DEC(__clang_patchlevel__)
# if defined(_MSC_VER)
    /* _MSC_VER = VVRR */
#  define SIMULATE_VERSION_MAJOR DEC(_MSC_VER / 100)
#  define SIMULATE_VERSION_MINOR DEC(_MSC_VER % 100)
# endif

#elif defined(__GNUC__) || defined(__GNUG__)
# define COMPILER_ID "GNU"
# if defined(__GNUC__)
#  define COMPILER_VERSION_MAJOR DEC(__GNUC__)
# else
#  define COMPILER_VERSION_MAJOR DEC(__GNUG__)
# endif
# if defined(__GNUC_MINOR__)
#  define COMPILER_VERSION_MINOR DEC(__GNUC_MINOR__)
# endif
# if defined(__GNUC_PATCHLEVEL__)
#  define COMPILER_VERSION_PATCH DEC(__GNUC_PATCHLEVEL__)
# endif

#elif defined(_MSC_VER)
# define COMPILER_ID "MSVC"
    /* _MSC_VER = VVRR */
# define COMPILER_VERSION_MAJOR DEC(_MSC_VER / 100)
# define COMPILER_VERSION_MINOR DEC(_MSC_VER % 100)
# if defined(_MSC_FULL_VER)
#  if _MSC_VER >= 1400
    /* _MSC_FULL_VER = VVRRPPPPP */
#   define COMPILER_VERSION_PATCH DEC(_MSC_FULL_VER % 100000)
#  else
    /* _MSC_FULL_VER = VVRRPPPP */
#   define COMPILER_VERSION_PATCH DEC(_MSC_FULL_VER % 10000)
#  endif
# endif
# if defined(_MSC_BUILD)
#  define COMPILER_VERSION_TWEAK DEC(_MSC_BUILD)
# endif

#elif defined(__VISUALDSPVERSION__) || defined(__ADSPBLACKFIN__) || defined(__ADSPTS__) || defined(__ADSP21000__)
# define COMPILER_ID "ADSP"
#if defined(__VISUALDSPVERSION__)
    /* __VISUALDSPVERSION__ = 0xVVRRPP00 */
# define COMPILER_VERSION_MAJOR HEX(__VISUALDSPVERSION__>>24)
# define COMPILER_VERSION_MINOR HEX(__VISUALDSPVERSION__>>16 & 0xFF)
# define COMPILER_VERSION_PATCH HEX(__VISUALDSPVERSION__>>8  & 0xFF)
#endif

#elif defined(__IAR_SYSTEMS_ICC__) || defined(__IAR_SYSTEMS_ICC)
# define COMPILER_ID "IAR"
# if defined(__VER__) && defined(__ICCARM__)
#  define COMPILER_VERSION_MAJOR DEC((__VER__) / 1000000)
#  define COMPILER_VERSION_MINOR DEC(((__VER__) / 1000) % 1000)
#  define COMPILER_VERSION_PATCH DEC((__VER__) % 1000)
#  define COMPILER_VERSION_INTERNAL DEC(__IAR_SYSTEMS_ICC__)
# elif defined(__VER__) && (defined(__ICCAVR__) || defined(__ICCRX__) || defined(__ICCRH850__) || defined(__ICCRL78__) || defined(__ICC430__) || defined(__ICCRISCV__) || defined(__ICCV850__) || defined(__ICC8051__) || defined(__ICCSTM8__))
#  define COMPILER_VERSION_MAJOR DEC((__VER__) / 100)
#  define COMPILER_VERSION_MINOR DEC((__VER__) - (((__VER__) / 100)*100))
#  define COMPILER_VERSION_PATCH DEC(__SUBVERSION__)
#  define COMPILER_VERSION_INTERNAL DEC(__IAR_SYSTEMS_ICC__)
# endif


    /* These compilers are either not known or too old to define an
       identification macro.  Try to identify the platform and guess that
       it is the native compiler.  */
#elif defined(__hpux) || defined(__hpua)
# define COMPILER_ID "HP"

#else /* unknown compiler */
# define COMPILER_ID "Unknown-Compiler"
#endif

    return util::to_string(COMPILER_ID)
#ifdef COMPILER_VERSION
      +" "+util::to_string(COMPILER_VERSION);
#elif defined(COMPILER_VERSION_MAJOR)
    +" "+util::to_string(COMPILER_VERSION_MAJOR)
# ifdef COMPILER_VERSION_MINOR
      +"."+util::to_string(COMPILER_VERSION_MINOR)
#  ifdef COMPILER_VERSION_PATCH
      +"."+util::to_string(COMPILER_VERSION_PATCH)
#   ifdef COMPILER_VERSION_TWEAK
      +"."+util::to_string(COMPILER_VERSION_TWEAK)
#   endif
#  endif
# endif
      ;
#endif

#undef DEC
#undef HEX
#undef COMPILER_ID
#undef COMPILER_VERSION
#undef COMPILER_VERSION_MAJOR
#undef COMPILER_VERSION_MINOR
#undef COMPILER_VERSION_PATCH
#undef COMPILER_VERSION_TWEAK
#undef SIMULATE_VERSION_MAJOR
#undef SIMULATE_VERSION_MINOR
#undef SIMULATE_VERSION_PATCH
#undef SIMULATE_VERSION_TWEAK
  }

  std::string gsSysInfo::getCppVersion()
  {
#if defined(_MSC_VER) && _MSC_VER < 1600
    return "C++ 199711L";
#elsif _MSC_VER >= 1900
    return "C++ "+util::to_string(_MSVC_LANG);
#elsif _MSC_VER >= 1600
    return "C++ 201103L";
#else
    return "C++ "+util::to_string(__cplusplus);
#endif
  }

  std::string gsSysInfo::getStdLibVersion()
  {
#ifdef _LIBCPP_VERSION
    return "libc++ "+util::to_string(_LIBCPP_VERSION);
#  elif defined(__GLIBCXX__)
    return "glibc++ "+util::to_string(__GLIBCXX__);
#  elif defined(__GLIBCPP__)
    return "glibc++ "+util::to_string(__GLIBCPP__);
#elif defined(__LIBCOMO__)
    return "Comeau STL "+util::to_string(__LIBCOMO__);
#  elif defined(__STL_CONFIG_H)
    return "SGI STL";
#  elif defined(__MSL_CPP__)
    return "MSL standard lib";
#  elif defined(__IBMCPP__)
    return "VACPP STL";
#  elif defined(MSIPL_COMPILE_H)
    return "Modena C++ STL";
#  elif (defined(_YVALS) && !defined(__IBMCPP__)) || defined(_CPPLIB_VER)
    return "Dinkumware STL "+util::to_string(_CPPLIB_VER);
#  elif defined(__STD_RWCOMPILER_H__) || defined(_RWSTD_VER)
    return "Rogue Wave lib "+util::to_string(_RWSTD_VER);
#else
    return "Unknown-STD";
#endif
  }

  std::string gsSysInfo::getExtraLibsVersion()
  {
    std::string s("");

    // CoDiPack extension
#if defined(CODI_VERSION)
    if (!s.empty()) s+= ", ";
    s += "CoDiPack "+util::to_string(CODI_VERSION);
#elif defined(CODI_MAJOR_VERSION) && \
      defined(CODI_MINOR_VERSION) && \
      defined(CODI_BUILD_VERSION)
        if (!s.empty()) s+= ", ";
    s += "CoDiPack "+util::to_string(CODI_MAJOR_VERSION)
      +          "."+util::to_string(CODI_MINOR_VERSION)
      +          "."+util::to_string(CODI_BUILD_VERSION);
#endif

    // GMP library
#if defined(__GNU_MP_VERSION)       &&   \
    defined(__GNU_MP_VERSION_MINOR) &&   \
    defined(__GNU_MP_VERSION_PATCHLEVEL)
    if (!s.empty()) s+= ", ";
    s += "gmp "+util::to_string(__GNU_MP_VERSION)
      +     "."+util::to_string(__GNU_MP_VERSION_MINOR)
      +     "."+util::to_string(__GNU_MP_VERSION_PATCHLEVEL);
#endif

    // IpOpt library
#if defined(IPOPT_VERSION)
    if (!s.empty()) s+= ", ";
    s += "IpOpt "+util::to_string(IPOPT_VERSION);
#elif defined(IPOPT_VERSION_MAJOR) && \
      defined(IPOPT_VERSION_MINOR) && \
      defined(IPOPT_VERSION_RELEASE)
    if (!s.empty()) s+= ", ";
    s += "IpOpt "+util::to_string(IPOPT_VERSION_MAJOR)
      +       "."+util::to_string(IPOPT_VERSION_MINOR)
      +       "."+util::to_string(IPOPT_VERSION_RELEASE);
#endif

    // Intel MKL library
#if defined(INTEL_MKL_VERSION)
    if (!s.empty()) s+= ", ";
    s += "MKL "+util::to_string(INTEL_MKL_VERSION);
#endif

    // MPFR library
#if defined(MPFR_VERSION_STRING)
    if (!s.empty()) s+= ", ";
    s += "mpfr "+util::to_string(MPFR_VERSION_STRING);
#elif defined(MPFR_VERSION_MAJOR) && \
      defined(MPFR_VERSION_MINOR) && \
      defined(MPFR_VERSION_PATCHLEVEL)
    if (!s.empty()) s+= ", ";
    s += "mpfr "+util::to_string(MPFR_VERSION_MAJOR)
      +      "."+util::to_string(MPFR_VERSION_MINOR)
      +      "."+util::to_string(MPFR_VERSION_PATCHLEVEL);
#endif

    // OpenCascade
#if defined(OCC_VERSION_COMPLETE)
    if (!s.empty()) s+= ", ";
    s += "occ "+util::to_string(OCC_VERSION_COMPLETE);
#elif defined(OCC_VERSION_MAJOR) && \
      defined(OCC_VERSION_MINOR) && \
      defined(OCC_VERSION_MAINTENANCE)
    if (!s.empty()) s+= ", ";
    s += "occ "+util::to_string(OCC_VERSION_MAJOR)
      +     "."+util::to_string(OCC_VERSION_MINOR)
      +     "."+util::to_string(OCC_VERSION_MAINTENANCE);
#endif

    // OpenNurbs
#if defined(OPENNURBS_VERSION)
    if (!s.empty()) s+= ", ";
    s += "onurbs "+util::to_string(OPENNURBS_VERSION);
#endif

    // Spectra library
#if defined(SPECTRA_MAJOR_VERSION) && \
    defined(SPECTRA_MINOR_VERSION) && \
    defined(SPECTRA_PATCH_VERSION)
    if (!s.empty()) s+= ", ";
    s += "spectra "+util::to_string(SPECTRA_MAJOR_VERSION)
      +         "."+util::to_string(SPECTRA_MINOR_VERSION)
      +         "."+util::to_string(SPECTRA_PATCH_VERSION);
#endif

    return s;
  }

  std::string gsSysInfo::getCpuInfo()
  {
#if defined(_WIN32) || defined(_WIN64)

    int CPUInfo[4] = {-1};
    unsigned   nExIds, i =  0;
    char CPUBrandString[0x40];

    __cpuid(CPUInfo, 0x80000000);
    nExIds = CPUInfo[0];

    for (i=0x80000000; i<=nExIds; ++i) {
      __cpuid(CPUInfo, i);
      if  (i == 0x80000002)
        memcpy(CPUBrandString, CPUInfo, sizeof(CPUInfo));
      else if  (i == 0x80000003)
        memcpy(CPUBrandString + 16, CPUInfo, sizeof(CPUInfo));
      else if  (i == 0x80000004)
        memcpy(CPUBrandString + 32, CPUInfo, sizeof(CPUInfo));
    }

    return CPUBrandString;

#elif __APPLE__

    std::string CPUBrandString;
    std::size_t size = 32;

    // Supply an oversized buffer, and avoid an extra call to sysctlbyname.
    CPUBrandString.resize(size);
    if (sysctlbyname("machdep.cpu.brand_string", &CPUBrandString[0], &size, NULL, 0) == 0 && size > 0) {
      if (CPUBrandString[size-1] == '\0')
        size--;
      CPUBrandString.resize(size);
      return CPUBrandString;
    }

#elif __linux__
#   if defined(__x86_64__) && ( defined(__GNUC__) || defined(__clang__) || defined(__INTEL_COMPILER) || defined(__SUNCC_PRO))

    char CPUBrandString[0x40];
    unsigned int CPUInfo[4] = {0,0,0,0};

    __cpuid(0x80000000, CPUInfo[0], CPUInfo[1], CPUInfo[2], CPUInfo[3]);
    unsigned int nExIds = CPUInfo[0];

    memset(CPUBrandString, 0, sizeof(CPUBrandString));

    for (unsigned int i = 0x80000000; i <= nExIds; ++i)
      {
        __cpuid(i, CPUInfo[0], CPUInfo[1], CPUInfo[2], CPUInfo[3]);

        if (i == 0x80000002)
          memcpy(CPUBrandString, CPUInfo, sizeof(CPUInfo));
        else if (i == 0x80000003)
          memcpy(CPUBrandString + 16, CPUInfo, sizeof(CPUInfo));
        else if (i == 0x80000004)
          memcpy(CPUBrandString + 32, CPUInfo, sizeof(CPUInfo));
      }

    return CPUBrandString;

#   else

    char hostname[HOST_NAME_MAX + 1];
    gethostname(hostname, HOST_NAME_MAX + 1);

    return "Unknown-CPU ["+hostname+"]";

#   endif

#elif __unix__

    // No generic implementation yet

#endif

    return "Unknown-CPU";
  }

  std::string gsSysInfo::getMemoryInfo()
  {
    uint64_t memsize = gsSysInfo::getMemoryInBytes();
    if (memsize>0) {
      if (memsize<1024)
        return util::to_string(memsize)+" B";
      else if (memsize<1024*1024)
        return util::to_string(memsize/1024)+" KB";
      else if (memsize<1024*1024*1024)
        return util::to_string(memsize/(1024*1024))+" MB";
      else
        return util::to_string(memsize/(1024*1024*1024))+" GB";
    }
    else
      return "Unknown-Memory";
  }

  uint64_t gsSysInfo::getMemoryInBytes()
  {
#if defined(_WIN32) || defined(_WIN64)

    MEMORYSTATUSEX status;
    status.dwLength = sizeof(status);
    GlobalMemoryStatusEx(&status);
    return (uint64_t)status.ullTotalPhys;

#elif __APPLE__

    int64_t memsize;
    std::size_t size = sizeof(memsize);

    if (sysctlbyname("hw.memsize", &memsize, &size, NULL, 0) == 0) {
      return (uint64_t)memsize;
    }

#elif __linux__

    long pages = sysconf(_SC_PHYS_PAGES);
    long page_size = sysconf(_SC_PAGE_SIZE);
    return (uint64_t)(pages * page_size);

#elif __unix__

    // No generic implementation yet

#endif

    return 0;
  }

} // namespace gismo
