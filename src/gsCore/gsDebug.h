/** @file gsDebug.h

    @brief This file contains the debugging and messaging system of G+Smo. 

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): A. Mantzaflaris
*/


// Start DEBUG_GROUP of Doxygen
/** @{ */

#pragma once

#include <iostream>
#include <sstream>
#include <iomanip>
#include <stdexcept>
#include <typeinfo>

// See also about memory leak detection:
// http://msdn.microsoft.com/en-us/library/e5ewb1h3%28v=vs.90%29.aspx
#if _MSC_VER //>= 1400
#include <crtdbg.h>
//#include <errno.h>
#endif

//#ifdef GISMO_EXTRA_DEBUG
//  #include <misc/gsStackWalker.h>
//#endif

namespace gismo {

/** Logging messages:
 *  gsInfo is ment to be the standard output stream, like for the output of the
 *  executables. In general, the library should not write to gsInfo.
 */
#define gsInfo std::cout

/** Logging messages:
 *  gsWarn is for warnings, eg, for missing functionality or problem in the input.
 *
 *  Note that gsWarn cannot be given as a parameter to another function.
 */
#define gsWarn std::cout<<"Warning: "
//#define gsWarn std::cerr

/** Logging messages:
 *  gsDebug and gsDebugVar(.) are for debugging messages and are enabled in debug
 *  mode only.
 *
 *  Note that gsDebug cannot be given as a parameter to another function.
 */
#ifndef  NDEBUG 

    #define gsDebug std::cout<<"GISMO_DEBUG: "

    #define gsDebugVar(variable) gsDebug << (strrchr(__FILE__, '/') ?          \
                             strrchr(__FILE__, '/') + 1 : __FILE__) <<":"<<    \
                              __LINE__<< ", "#variable": \n"<<(variable)<<"\n"
#define gsDebugIf(cond,variable) if (cond) gsDebug <<"[ "#cond" ] -- "<<       \
              (strrchr(__FILE__, '/') ? strrchr(__FILE__, '/') + 1 : __FILE__) \
               <<":"<<__LINE__<< ", "#variable": \n"<<(variable)<<"\n"
#else
    #define gsDebug if (0) std::cout
    #define gsDebugVar(variable)
    #define gsDebugIf(cond,variable)
#endif

/** 
 *  Used for optional inclusion of .hpp header files in the .h files.
 *  Allows to drop dependencies on the .hpp files when using
 *  GISMO_BUILD_LIB for compiling a library instance.  
 *  When compiling as a pure template library CMake will detect
 *  dependency on .hpp files.
 */
#define GISMO_HPP_HEADER(x) #x

/**  
 *  Runtime assertions which display a message  
 *
 */
#ifndef NDEBUG
#   define GISMO_ASSERT(cond, message) do if(!(cond)) {std::stringstream _m_;\
       _m_<<"GISMO_ASSERT `"<<#cond<<"` "<<message<<"\n"<<__FILE__<<", line "\
        <<__LINE__<<" ("<<__FUNCTION__<<")";                                 \
       throw std::logic_error(_m_.str()); } while(false)
#else
#   define GISMO_ASSERT(condition, message)
#endif

/**  
 *  Runtime check and display error message. This command is the same as
 *  GISMO_ASSERT but it is executed in release builds as well.
 *
 */
#define GISMO_ENSURE(cond, message) do if(!(cond)) {std::stringstream _m_;   \
    _m_<<"GISMO_ENSURE `"<<#cond<<"` "<<message<<"\n"<<__FILE__<<", line "   \
     <<__LINE__<<" ("<< __FUNCTION__<< ")";                                  \
    throw std::runtime_error(_m_.str());} while(false)

/**  
 *  Denote a variable as unused, used to silence warnings in release
 *  mode builds.
 *
 */
#define GISMO_UNUSED(x)  static_cast<void>(x)

/**
 *  Runtime error message
 *
 */
#define GISMO_ERROR(message) do {std::stringstream _m_; _m_<<"GISMO_ERROR "  \
    <<message<<"\n"<<__FILE__<<", line " <<__LINE__<<" ("<<__FUNCTION__<<")";\
    throw std::runtime_error(_m_.str());} while(false)

/**  
 *  Runtime "no implementation" error happens when the user calls a

 *  virtual member function without a default implementation.
 */
 
// TO DO: for GCC __PRETTY_FUNC__ is better
#define GISMO_NO_IMPLEMENTATION {std::stringstream _m_;                            \
    _m_<<"Virtual member function `"<<__FUNCTION__<<"` has not been implemented\n" \
     <<__FILE__<<", line "<<__LINE__<<"\n"<<typeid(*this).name();                  \
    throw std::runtime_error(_m_.str());}

/*
#ifdef _MSC_VER
#include <float.h>
template <typename T> int isnan   (T a) {return _isnan(a); }
template <typename T> int isfinite(T a){return _finite(a);}
template <typename T> bool isinf(T a) {return (_FPCLASS_PINF|_FPCLASS_NINF) & _fpclass(a);}
 #else
#ifdef _INTEL_COMPILER
#include <mathimf.h>
#else
using std::isnan;
using std::isfinite;
using std::isinf;
#endif
*/
/**
   Check if a floating point number is different than NAN (not a number)

   See https://en.wikipedia.org/wiki/Floating_point#Special_values
   and https://en.wikipedia.org/wiki/NaN
 */
template <typename T> bool gsIsnumber(T a) {return a == a;}
template <typename T> bool gsIsnan   (T a) {return a != a;}
/**
   Check if a flaoting point number is different than INF

   See https://en.wikipedia.org/wiki/Floating_point#Special_values
 */
template <typename T> bool gsIsfinite(T a) {return  (a - a) == (a - a);}

}//namespace gismo

/*
  Disable debug/abort popup windows on MS Windows
  
  See http://msdn.microsoft.com/en-us/library/1y71x448.aspx

  You might also need to disable "error reporting" on your windows
  system for popup-free runs.
*/
#if _MSC_VER //>= 1400 
static const int    gismo_CrtSetReportMode = _CrtSetReportMode( 
    _CRT_ASSERT, _CRTDBG_MODE_FILE   );
static const _HFILE gismo_CrtSetReportFile = _CrtSetReportFile( 
    _CRT_ASSERT, _CRTDBG_FILE_STDERR );
static const int  gismo_set_abort_behavior = _set_abort_behavior( 
    0x0, _WRITE_ABORT_MSG | _CALL_REPORTFAULT);
#endif

/*
  Disable some Warnings
*/

#ifdef _MSC_VER
// 4100 - unreferenced formal parameter
// 4101 - unreferenced local variable
// 4127 - conditional expression is constant (triggered by assertion macros)
// 4146 - unary minus operator applied to unsigned type, result still unsigned
// 4181 - qualifier applied to reference type ignored
// 4211 - nonstandard extension used : redefined extern to static
// 4231 - nonstandard extension used: extern before template explicit instantiation
// 4244 - 'argument' : conversion from 'type1' to 'type2', possible loss of data
// 4251 - needs to have dll-interface to be used by clients of class
// 4273 - QtAlignedMalloc, inconsistent DLL linkage
// 4275 - non dll-interface base 
// 4324 - structure was padded due to declspec(align())
// 4428 - universal-character-name encountered in source
// 4503 - decorated name length exceeded
// 4505 - unreferenced local function has been removed
// 4512 - assignment operator could not be generated
// 4522 - 'class' : multiple assignment operators specified
// 4566 - character represented by universal-character-name cannot be represented in the current code page
// 4661 - no definition available
// 4700 - uninitialized local variable 'xyz' used
// 4702 - unreachable code
// 4714 - function marked as __forceinline not inlined
// 4717 - 'function' : recursive on all control paths, function will cause runtime stack overflow
// 4789 - destination of memory copy is too small (for Eigen)
// 4996 - 'sprintf': This function or variable may be unsafe. Consider using sprintf_s instead.
// 4510 - default constructor could not be generated
// 4610 - user defined constructor required
  #pragma warning( push )
  #pragma warning( disable : 4100 4127 4146 4231 4251 4428 4275 4503 4505 4512 4566 4661 4714 4789 4996 4510 4610)

#elif defined __INTEL_COMPILER
// 2196 - routine is both "inline" and "noinline" ("noinline" assumed)
//        ICC 12 generates this warning even without any inline keyword, when defining class methods 'inline' i.e. inside of class body
//        typedef that may be a reference type.
// 279  - controlling expression is constant
//        ICC 12 generates this warning on assert(constant_expression_depending_on_template_params) and frankly this is a legitimate use case.
// 161  - unrecognized pragma
// 175  - subscript out of range
//        to avoid warnings on #pragma GCC diagnostic          
  #pragma warning push
  #pragma warning disable 2196 279 161 175

#elif defined __clang__
// -Wconstant-logical-operand - warning: use of logical && with constant operand; switch to bitwise & or remove constant
// -Wbind-to-temporary-copy - warning: Warn about an unusable copy constructor when binding a reference to a temporary
  #pragma clang diagnostic push
  #pragma clang diagnostic ignored "-Wconstant-logical-operand"
  #pragma clang diagnostic ignored "-Wbind-to-temporary-copy"

#elif defined __GNUC__ // major version >=4
// typedef locally defined but not used [-Wunused-local-typedefs]
#if ( __GNUC__>4 || (__GNUC__==4 && __GNUC_MINOR__>7) )
#pragma GCC diagnostic ignored "-Wunused-local-typedefs"
#endif

#if (__cplusplus < 201703L && __GNUC__>6)
// mangled name will change in C++17 because the exception
// specification is part of a function type [-Wnoexcept-type]
#pragma GCC diagnostic ignored "-Wnoexcept-type"
#endif

#endif


/*
   Compile-time assertions: 
  
  - in GISMO_STATIC_ASSERT(CONDITION,MSG) the parameter CONDITION
     must be a compile time boolean expression, and MSG an error
     message (string)
 
   - define GISMO_NO_STATIC_ASSERT to disable them (and save
     compilation time) in that case, the static assertion is
     converted to a runtime assertion.
 
  - GISMO_STATIC_ASSERT can only be used in function scope
 */
#ifndef GISMO_NO_STATIC_ASSERT

  #if defined(__GXX_EXPERIMENTAL_CXX0X__) || _MSC_VER >= 1600

    // Native static_assert is available
    #define GISMO_STATIC_ASSERT(X,MSG) static_assert(X,#MSG);

  #else // not C++11

    namespace gismo {
    namespace internal {

    template<bool condition> struct static_assertion {};
    template<> struct static_assertion<true> { enum { STATIC_ASSERTION_FAILED }; };

    } // end namespace internal
    } // end namespace gismo

    // Specialized implementation for MSVC to avoid "conditional
    // expression is constant" warnings.  This implementation doesn't
    // appear to work under GCC, hence the multiple implementations.
    #ifdef _MSC_VER

      #define GISMO_STATIC_ASSERT(CONDITION,MSG) \
        {gismo::internal::static_assertion<bool(CONDITION)>::STATIC_ASSERTION_FAILED;}

    #else

      #define GISMO_STATIC_ASSERT(CONDITION,MSG) \
        if (gismo::internal::static_assertion<bool(CONDITION)>::STATIC_ASSERTION_FAILED) {}

    #endif

  #endif // not C++11

#else // GISMO_NO_STATIC_ASSERT

#define GISMO_STATIC_ASSERT(CONDITION,MSG) GISMO_ASSERT(CONDITION, #MSG);

#endif // GISMO_NO_STATIC_ASSERT



#ifdef GISMO_WARNINGS
    //#pragma message("G+Smo Warnings ON")
  #ifdef __GNUC__
    #define GISMO_DEPRECATED __attribute__((deprecated))
  #elif defined(_MSC_VER)
    #define GISMO_DEPRECATED __declspec(deprecated)
  #else
    #define GISMO_DEPRECATED
  #endif
#else
  #define GISMO_DEPRECATED
#endif


// Next line closes the DEBUG_GROUP of Doxygen
/** @} */
